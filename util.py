"""
util.py — unified utility module for MotifDiff and MotifScore
=============================================================
Merged from:
  - util.py      (SegmentDataBed, SegmentDataSeq — used by MotifScore)
  - util(1).py   (vcfData, gzip support, strand_specific — used by MotifDiff)

N / gap-character bug fix
--------------------------
Problem: returnonehot() encodes N (and '-') as the all-zero vector [0,0,0,0].
Because PWM kernels store log(p_motif) ≤ 0, any window that overlaps an N
scores exactly 0.0 — which is *above* the maximum achievable score for a real
sequence.  The probNorm spline therefore maps those positions to CDF(0) ≈ 1.0,
producing spurious binding probabilities of ~1.

Fix applied here:
  1. returnonehot() now treats BOTH 'N' and '-' as missing (zero-encoded, same
     as before), but also returns an optional boolean mask of invalid columns
     when called with return_invalid=True.
  2. get_invalid_mask(seq, kernel_len, dinucleotide) — new helper that returns
     a 1-D boolean array of length (len(seq) - kernel_len + 1) marking every
     output position whose receptive field overlaps at least one N or '-'.
  3. apply_invalid_mask(scores, invalid_mask) — sets contaminated positions in
     the conv1d output tensor to -inf so that max-pool and mc_spline never see
     them.  mc_spline itself is unchanged; callers are responsible for masking
     before calling it.

Usage in MotifDiff / MotifScore (see those files for the full integration):
    scores = F.conv1d(one_hot_batch, kernels)          # shape (B, 2*M, L)
    for b in range(B):
        inv = get_invalid_mask(raw_seq[b], kernel_len, dinucleotide)
        apply_invalid_mask(scores[b], inv)             # in-place, -inf
    # now max-pool / mc_spline as usual
"""

import numpy as np
import torch
import pandas as pd
from pysam import FastaFile
from Bio.Seq import Seq
import time
import itertools
import xml.etree.ElementTree as ET
import os
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator
import regex as re
import itertools as itt
from collections import namedtuple
import gzip

torch.set_printoptions(precision=8)
np.set_printoptions(precision=8)


# ── I/O helpers ───────────────────────────────────────────────────────────────

def open_maybe_gzip(path):
    """Open a plain or gzip-compressed text file transparently."""
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path, "r")


def number_of_headers(filename):
    header = 0
    with open_maybe_gzip(filename) as file:
        while True:
            line = file.readline()
            if line.startswith("#"):
                header += 1
            else:
                break
    return header


def kmers_count(seq, k=2):
    lookup = {"".join(i): 0 for i in itertools.product(["A", "C", "G", "T"], repeat=k)}
    mers = [seq[i:i+2] for i in range(len(seq)-k+1)]
    for i in mers:
        if i in lookup:
            lookup[i] += 1
    for i in lookup:
        lookup[i] /= (len(seq)-k+1)
    return list(lookup.values())


def kmers(k=2):
    return ["".join(i) for i in itertools.product(["A", "C", "G", "T"], repeat=k)]


def logit(x, a, b):
    return 1 / (1 + np.exp(-a * x - b))


def logit_torch(x, a, b):
    return 1 / (1 + torch.exp(-a * x - b))


def init_dist(dmin, dmax, dp, weights, probs):
    out = np.zeros(int(np.round((dmax - dmin) / dp) + 1))
    ii = np.array(np.round((weights - dmin) / dp), dtype=int)
    for i in range(len(probs)):
        out[ii[i]] = out[ii[i]] + probs[i]
    return out


def scoreDist(pwm, nucleotide_prob=None, gran=None, size=1000):
    if nucleotide_prob is None:
        nucleotide_prob = [np.ones(4) / 4] * pwm.shape[0]
    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pwm) - np.min(pwm)) / (size - 1)
    pwm = np.round(pwm / gran) * gran
    pwm_max, pwm_min = pwm.max(axis=1), pwm.min(axis=1)
    distribution = init_dist(pwm_min[0], pwm_max[0], gran, pwm[0], nucleotide_prob[0])
    for i in range(1, pwm.shape[0]):
        kernel = init_dist(pwm_min[i], pwm_max[i], gran, pwm[i], nucleotide_prob[i])
        distribution = np.convolve(distribution, kernel)
    support_min = pwm_min.sum()
    ii = np.where(distribution > 0)[0]
    support = support_min + (ii) * gran
    return support, distribution[ii]


# ── Dinucleotide helpers ───────────────────────────────────────────────────────

class diNucMat:
    def __init__(self, values, colnames) -> None:
        assert colnames == ["".join(i) for i in itt.product(["A", "C", "G", "T"], repeat=2)]
        self.values = values
        self.colnames = colnames

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        assert len(values.shape) == 2
        assert values.shape[1] == 16
        assert values.shape[0] > 0
        self._values = values
        self._colMats = np.reshape(self.values, newshape=(values.shape[0], 4, 4))

    @property
    def colMats(self):
        return self._colMats

    @property
    def colnames(self):
        return self._colnames

    @colnames.setter
    def colnames(self, colnames):
        self._colnames = colnames


class diNucProbMat(diNucMat):
    @diNucMat.values.setter
    def values(self, values):
        assert len(values.shape) == 2
        assert values.shape[1] == 16
        assert values.shape[0] > 0
        assert np.all(values >= 0)
        assert np.all(values <= 1)
        assert np.all(np.isclose(np.sum(values, axis=1), 1))
        self._values = values
        self._colMats = np.reshape(self.values, newshape=(self.values.shape[0], 4, 4))
        self._trnMats = self._colMats / np.sum(self._colMats, axis=2, keepdims=True)

    @property
    def trnMats(self) -> np.array:
        return self._trnMats


def diNucMotDist(pssm, prob, gran=None, size=1000):
    assert isinstance(pssm, diNucMat)
    assert isinstance(prob, diNucProbMat)
    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pssm) - np.min(pssm)) / (size - 1)

    def vals2inds(vals, mnscore, gran):
        return np.rint(((vals - mnscore) / gran)).astype(int)

    mnscore = np.floor(np.sum(np.min(pssm.values, axis=1)))
    mxscore = np.ceil(np.sum(np.max(pssm.values, axis=1)))
    if mxscore < 0:
        mxscore = 0
    nscores = int(np.rint((mxscore - mnscore) / gran) + 1)

    SD = np.zeros((4, nscores))
    for i in range(4):
        for j in range(4):
            SD[i, vals2inds(pssm.colMats[0, j, i], mnscore, gran)] += prob.colMats[0, j, i]

    SD_tmp = SD.copy()
    for pos in range(1, pssm.values.shape[0]):
        for nuc in range(4):
            tvec = np.zeros(nscores)
            scores = pssm.colMats[pos, :, nuc]
            shifts = np.rint(scores / gran).astype(int)
            for i in range(4):
                tvec += np.roll(SD_tmp[i, :], shifts[i]) * prob.trnMats[pos, i, nuc]
            SD[nuc, :] = tvec
        SD_tmp = SD.copy()

    x = np.arange(mnscore, mxscore + gran, gran)
    y = np.sum(SD, axis=0)
    ii = np.where(y > 0)[0]
    support = mnscore + (ii) * gran
    x = support
    y = y[ii]
    Dist = namedtuple('Dist', ['x', 'y'])
    return Dist(x, y)


def mono2di(ppm):
    num_rows = ppm.shape[0]
    num_cols = ppm.shape[1] ** 2
    ppm_di = np.zeros((num_rows - 1, num_cols))
    for i in range(num_rows - 1):
        for j in range(ppm.shape[1]):
            ppm_di[i, 4*j:4*j+4] = ppm[i, j] * ppm[i+1, :]
    return ppm_di


def scoreDistDinuc(pwm, gran=None, size=1000):
    tmp = pwm
    cn = ["".join(i) for i in itt.product(["A", "C", "G", "T"], repeat=2)]
    pssm = diNucMat(tmp, cn)
    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pwm) - np.min(pwm)) / (size - 1)
    tmp = np.exp(tmp)
    tmp = tmp / np.sum(tmp, axis=1, keepdims=True)
    prob = diNucProbMat(tmp, cn)
    sd_mot = diNucMotDist(pssm, prob, gran=0.01)
    avg_dinuc_freqs = prob.values.mean(axis=0)
    iid_prob_values = np.repeat(avg_dinuc_freqs, pwm.shape[0]).reshape((pwm.shape[0], 16), order='F')
    prob_bg1 = diNucProbMat(iid_prob_values, cn)
    sd_bg1 = diNucMotDist(pssm, prob_bg1, gran=0.01)
    return (sd_mot.x, sd_bg1.y, sd_mot.y)


# ── Spline fitting & application ───────────────────────────────────────────────

def MCspline_fitting(pwms, nucleotide_prob=None, gran=None, size=1000, nuc="mono", method="motif_based"):
    spline_list = []
    for i in range(0, pwms.shape[0], 2):
        pwm = pwms[i].numpy().T
        pwm = pwm[pwm.sum(axis=1) != 0, :]
        nucleotide_prob = np.exp(pwm) / np.sum(np.exp(pwm), axis=1, keepdims=True)
        if nuc == "mono":
            s, d = scoreDist(pwm, nucleotide_prob, gran, size)
            spl = PchipInterpolator(s, np.cumsum(d))
        if nuc == "di":
            s, d_iid, d_m = scoreDistDinuc(pwm, gran=gran, size=size)
            if method == "iid":
                spl = PchipInterpolator(s, np.cumsum(d_iid))
            if method == "motif_based":
                spl = PchipInterpolator(s, np.cumsum(d_m))
            if method == "mixture":
                spl = PchipInterpolator(s, np.cumsum(0.5*d_m + 0.25*d_iid + 0.25*1/len(d_iid)))
        spline_list.append(spl)
    return spline_list


def mc_spline(mat, spline_list):
    out = torch.empty_like(mat)
    assert mat.shape[1] == len(spline_list)
    for i in range(len(spline_list)):
        spl = spline_list[i]
        out_i = spl(mat[:, i])
        out_i[out_i > 1] = 1
        out_i[out_i < 0] = 0
        out[:, i] = torch.tensor(out_i)
    return out


# ── N / gap masking (THE BUG FIX) ─────────────────────────────────────────────

# Characters that should be treated as "missing" / ambiguous.
_INVALID_CHARS = frozenset({'N', '-'})


def returnonehot(string, dinucleotide=False, return_invalid=False):
    """
    Convert a DNA string to a one-hot (or dinucleotide) tensor.

    N and '-' (gap) characters are treated as missing: the corresponding
    column(s) in the output are left as all-zeros.  This is correct for the
    one-hot matrix itself, but callers that compute PWM scores must mask those
    positions afterwards (see get_invalid_mask / apply_invalid_mask).

    Parameters
    ----------
    string : str
        DNA sequence (may contain N, '-', upper or lower case).
    dinucleotide : bool
        If True, return a (16, len-1) dinucleotide one-hot matrix.
    return_invalid : bool
        If True, also return a boolean 1-D array marking columns that contain
        at least one N or '-' character (True = invalid).

    Returns
    -------
    out : np.ndarray, shape (4 or 16, L)
    invalid_cols : np.ndarray of bool, shape (L,)  — only if return_invalid=True
    """
    string = string.upper()
    tmp = np.array(list(string))

    if dinucleotide:
        lookup = {"".join(i): n for n, i in enumerate(itertools.product(["A", "C", "G", "T"], repeat=2))}
        # positions that are N or '-'
        bad = np.where(np.isin(tmp, list(_INVALID_CHARS)))[0]
        # a dinucleotide position is invalid if either of its two bases is bad
        bad_di = np.unique(np.clip(np.concatenate([bad, bad - 1]), 0, len(tmp) - 2))
        valid_di = np.where(~np.isin(np.arange(len(tmp) - 1), bad_di))[0]
        tmp_di = np.array([tmp[i] + tmp[i+1] for i in range(len(tmp) - 1)])
        irow = np.array([lookup[t] for t in tmp_di[valid_di]])
        out = np.zeros((16, len(tmp) - 1), dtype=np.float32)
        if len(valid_di) > 0:
            out[irow, valid_di] = 1
        if return_invalid:
            invalid_cols = np.isin(np.arange(len(tmp) - 1), bad_di)
            return out, invalid_cols
        return out

    else:
        lookup = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        # positions that are N or '-' (or any other non-ACGT character)
        bad = np.where(np.isin(tmp, list(_INVALID_CHARS)))[0]
        valid = np.where(~np.isin(tmp, list(_INVALID_CHARS)))[0]
        # further filter: only keep positions whose character is in lookup
        valid = np.array([i for i in valid if tmp[i] in lookup])
        irow = np.array([lookup[tmp[i]] for i in valid])
        out = np.zeros((4, len(tmp)), dtype=np.float32)
        if len(valid) > 0:
            out[irow, valid] = 1
        if return_invalid:
            invalid_cols = np.isin(np.arange(len(tmp)), list(_INVALID_CHARS) + [c for c in np.unique(tmp) if c not in lookup])
            # simpler: a column is invalid if it is all-zero AND the original char was not ACGT
            invalid_cols = np.zeros(len(tmp), dtype=bool)
            invalid_cols[bad] = True
            # also mark any character not in lookup (e.g. ambiguity codes other than N)
            other_bad = np.where(~np.isin(tmp, ['A', 'C', 'G', 'T', 'N', '-']))[0]
            invalid_cols[other_bad] = True
            return out, invalid_cols
        return out


def get_invalid_mask(seq, kernel_len, dinucleotide=False, motif_len=None):
    """
    Return a boolean array of length n_output_positions marking every conv1d
    output position whose receptive field overlaps at least one N or '-'.

    For a sequence of length L and a padded kernel of length K:
      n_output_positions = L - K + 1   (mono)
      n_output_positions = (L-1) - K + 1  (dinucleotide)

    The contamination radius is determined by the *actual* motif length, not
    the padded kernel length.  Pass motif_len explicitly when the kernel is
    zero-padded (as in MEME_probNorm / gimme PFM files).  If motif_len is None
    it defaults to kernel_len.

    Parameters
    ----------
    seq : str
        Raw DNA sequence (may contain N, '-').
    kernel_len : int
        Padded kernel length (= kernels.shape[2]).  Used only to compute
        n_output_positions.
    dinucleotide : bool
    motif_len : int or None
        Actual motif length (number of non-zero columns in the kernel).
        Defaults to kernel_len when None.

    Returns
    -------
    invalid : np.ndarray of bool, shape (n_output_positions,)
        True where the output position is contaminated by N or '-'.
    """
    if motif_len is None:
        motif_len = kernel_len

    seq = seq.upper()
    arr = np.array(list(seq))
    bad_pos = np.where(np.isin(arr, list(_INVALID_CHARS)))[0]

    if dinucleotide:
        seq_len = len(seq) - 1          # di-encoded length
    else:
        seq_len = len(seq)

    n_out = seq_len - kernel_len + 1
    if n_out <= 0:
        return np.ones(max(0, n_out), dtype=bool)

    invalid = np.zeros(n_out, dtype=bool)
    if len(bad_pos) == 0:
        return invalid

    for bp in bad_pos:
        if dinucleotide:
            # dinucleotide position i covers bases i and i+1
            # output position j covers di-positions j .. j+motif_len-1
            # di-position d is bad if base d or d+1 is bad
            bad_di = set()
            if bp > 0:
                bad_di.add(bp - 1)
            bad_di.add(bp)
            for d in bad_di:
                # output positions j where d is in [j, j+motif_len-1]
                j_min = max(0, d - motif_len + 1)
                j_max = min(n_out - 1, d)
                if j_min <= j_max:
                    invalid[j_min:j_max + 1] = True
        else:
            # output position j covers bases j .. j+motif_len-1
            j_min = max(0, bp - motif_len + 1)
            j_max = min(n_out - 1, bp)
            if j_min <= j_max:
                invalid[j_min:j_max + 1] = True

    return invalid


def apply_invalid_mask(scores, invalid_mask):
    """
    Set contaminated positions in a conv1d output to -inf IN PLACE.

    Parameters
    ----------
    scores : torch.Tensor, shape (n_kernels, n_positions)
        Output of F.conv1d for a single sequence (one element of the batch).
        n_positions must equal len(invalid_mask).
    invalid_mask : np.ndarray of bool, shape (n_positions,)
        True where the position is contaminated by N or '-'.
    """
    if not invalid_mask.any():
        return
    idx = torch.from_numpy(np.where(invalid_mask)[0]).long()
    # Guard against index out of bounds (e.g. mask computed with wrong kernel_len)
    idx = idx[idx < scores.shape[1]]
    if len(idx) > 0:
        scores[:, idx] = -torch.inf


# ── VCF / BED / sequence data loaders ─────────────────────────────────────────

def readvcf(filename):
    nh = number_of_headers(filename)
    compression = "gzip" if filename.endswith((".gz", ".bgz")) else None
    if nh > 1:
        data = pd.read_csv(filename, skiprows=nh, header=None, sep="\t", compression=compression)
    elif nh == 1:
        data = pd.read_csv(filename, skiprows=1, header=None, sep="\t", compression=compression)
    else:
        data = pd.read_csv(filename, header=None, sep="\t", compression=compression)
    return data


def readbed(filename, up):
    data = pd.read_csv(filename, sep="\t", header=None)
    chrs = data[0].to_numpy()
    start = data[1].to_numpy(dtype=int)
    end = data[2].to_numpy(dtype=int)
    if data.shape[1] > 3:
        gene = data[3].to_numpy(dtype=str)
        if data.shape[1] > 5:
            print("Strand detected")
            up = int(np.floor(up))
            strand = data[5].to_numpy()
            start = start - (strand == "+") * up
            end = end + (strand == "-") * up
    else:
        gene = np.array([None] * len(chrs))
    return chrs, start, end, gene


def returnmask(i, mask, windowsize, start, end, dinucleotide):
    if dinucleotide:
        tmp = np.zeros(mask.shape[2] + 1)
        tmp[int(windowsize-1):int(end-start-windowsize+1)] = 1
        mask[i, :, :] = torch.from_numpy(np.convolve(tmp, [1, 1], mode="valid"))
    else:
        mask[i, :, int(windowsize-1):int(end-start-windowsize+1)] = 1


def countlowercase(arr):
    return sum([1 for c in arr if c.islower()])


def stringstats(string):
    lowercaseratio = countlowercase(string) / len(string)
    string = string.upper()
    tmp = np.array(list(string))
    gccount = np.sum(np.logical_or(tmp == 'C', tmp == 'G')) / len(tmp)
    patterns = kmers_count(string)
    return np.array([gccount, lowercaseratio, *patterns], dtype=np.float32)


# ── PWM / kernel utilities ─────────────────────────────────────────────────────

def read_pwm(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    values = []
    for line in lines:
        if not line.startswith(">"):
            values.append(line.strip().split("\t"))
    values = np.array(values, dtype=float)
    if np.min(values) >= 0:
        values = values / values.sum(axis=1, keepdims=True)
    return np.array(values, dtype=float)


def transform_kernel(kernel, smoothing, background):
    if np.min(kernel) < 0:
        out = kernel
    else:
        out = np.log(kernel / background + smoothing)
    c = out.max(axis=1)
    out = out - c[:, np.newaxis]
    norm = out.min(axis=1).sum()
    return out, norm


# ── MEME parsers ───────────────────────────────────────────────────────────────

class MEME_probNorm():
    def __init__(self, precision=1e-7, smoothing=0.02, background=None):
        self.version = 0
        self.alphabet = ""
        self.strands = ""
        self.background = []
        self.names = []
        self.nmotifs = 0
        self.precision = precision
        self.smoothing = smoothing
        self.background_prob = background

    def parse(self, text, nuc="mono", transform=False, strand_specific=False):
        if nuc == "mono":
            if self.background_prob is None:
                background_prob = np.ones(4) / 4
            else:
                background_prob = self.background

            if text.endswith(".pfm") or text.endswith(".ppm"):
                print("motif is pfm or ppm format")
                with open(text, 'r') as file:
                    data = file.read()
                self.names = re.findall(r"(>GM\.5\.0\.\S+)", data)
                self.synonames = re.findall(r"(#GM\.5\.0\.\S+)", data)
                letter_probs = re.findall(
                    r"(>GM.*\n((?:[ \t]*\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*"
                    r"[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]*\n)+))", data)
                assert len(letter_probs) == len(self.names)
                self.nmotifs = len(self.names)
                out_channels = self.nmotifs * 2
                in_channels = 4
                matrices = []
                length = 0
                for i in range(len(letter_probs)):
                    matrix = letter_probs[i][0].split("\n")
                    matrix = matrix[1:-1] if len(matrix[-1]) == 0 else matrix[1:]
                    matrices.append(np.array([r.split() for r in matrix], dtype=float))
                    if matrices[-1].shape[0] > length:
                        length = matrices[-1].shape[0]

            if text.endswith(".meme"):
                with open(text, 'r') as file:
                    data = file.read()
                self.version = re.compile(r'MEME version ([\d+\.*]+)').match(data).group(1)
                self.names = re.findall(r"MOTIF (.*)\n", data)
                self.background = re.findall(
                    r"Background letter frequencies.*\n(A .* C .* G .* T .*)\n", data)[0]
                self.strands = re.findall(r"strands: (.*)\n", data)[0].strip()
                self.alphabet = re.findall(r"ALPHABET=(.*)\n", data)[0].strip()
                letter_probs = re.findall(
                    r"letter-probability.*\n((?:[ \t]*\d*\.?\d+[eE]?-?\d*[ \t]+"
                    r"\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]+"
                    r"\d*\.?\d+[eE]?-?\d*[ \t]*\n)+)", data)
                assert len(letter_probs) == len(self.names)
                self.nmotifs = len(letter_probs)
                out_channels = self.nmotifs * 2
                in_channels = 4
                matrices = []
                length = 0
                for i in range(len(letter_probs)):
                    matrix = letter_probs[i].split("\n")
                    matrix = matrix[:-1] if len(matrix[-1]) == 0 else matrix
                    matrices.append(np.array([r.split() for r in matrix], dtype=float))
                    if matrices[-1].shape[0] > length:
                        length = matrices[-1].shape[0]

            if os.path.isdir(text):
                self.names = os.listdir(text)
                self.nmotifs = len(self.names)
                in_channels = 4
                out_channels = self.nmotifs * 2
                matrices = []
                length = 0
                for k, i in enumerate(self.names):
                    if i.endswith(".pcm") or i.endswith(".pwm"):
                        matrix = read_pwm(os.path.join(text, i))
                        matrices.append(matrix)
                        if matrix.shape[0] > length:
                            length = matrix.shape[0]

        if nuc == "di":
            if self.background_prob is None:
                background_prob = np.ones(16) / 16
            else:
                background_prob = self.background_prob

            if text.endswith(".meme"):
                with open(text, 'r') as file:
                    data = file.read()
                self.version = re.compile(r'MEME version ([\d+\.*]+)').match(data).group(1)
                self.names = re.findall(r"MOTIF (.*)\n", data)
                self.background = re.findall(
                    r"Background letter frequencies.*\n(A .* C .* G .* T .*)\n", data)[0]
                self.strands = re.findall(r"strands: (.*)\n", data)[0].strip()
                self.alphabet = re.findall(r"ALPHABET=(.*)\n", data)[0].strip()
                letter_probs = re.findall(
                    r"(letter-probability.*\n([ \t]*\d+\.?\d*[ \t]+\d+\.?\d*"
                    r"[ \t]+\d+\.?\d*[ \t]+\d+\.?\d*[ \t]*\n)+)", data)
                assert len(letter_probs) == len(self.names)
                self.nmotifs = len(letter_probs)
                out_channels = self.nmotifs * 2
                in_channels = 16
                matrices = []
                length = 0
                for i in range(len(letter_probs)):
                    matrix = letter_probs[i][0].split("\n")
                    matrix = matrix[1:-1] if len(matrix[-1]) == 0 else matrix[1:]
                    m = np.array([r.split() for r in matrix], dtype=float)
                    if m.shape[1] == 4:
                        m = mono2di(m)
                    matrices.append(m)
                    if matrices[-1].shape[0] > length:
                        length = matrices[-1].shape[0]
            else:
                self.names = os.listdir(text)
                self.nmotifs = len(self.names)
                in_channels = 16
                out_channels = self.nmotifs * 2
                matrices = []
                length = 0
                for k, i in enumerate(self.names):
                    if i.endswith(".dpcm") or i.endswith(".dpwm"):
                        matrix = read_pwm(os.path.join(text, i))
                        matrices.append(matrix)
                        if matrix.shape[0] > length:
                            length = matrix.shape[0]

        out = np.zeros((out_channels, in_channels, length), dtype=np.float32)
        mask = torch.zeros((out_channels, 1, length), dtype=torch.uint8)
        for k, kernel in enumerate(matrices):
            if transform:
                kernel, _ = transform_kernel(kernel, self.smoothing, background_prob)
            else:
                if np.min(kernel) < 0:
                    kernel = kernel
                else:
                    kernel[kernel == 0] = self.precision
                    kernel = np.log(kernel)

            if strand_specific:
                out[2*k,   :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel.T
            else:
                out[2*k,   :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel[::-1, ::-1].T
            mask[2*k,   :, :kernel.shape[0]] = 1
            mask[2*k+1, :, :kernel.shape[0]] = 1

        return torch.from_numpy(out), mask


class MEME_FABIAN():
    def __init__(self, precision=1e-7, smoothing=0.02, background=None):
        self.version = 0
        self.alphabet = ""
        self.strands = ""
        self.background = []
        self.names = []
        self.nmotifs = 0
        self.precision = precision
        self.smoothing = smoothing
        self.background_prob = background

    def parse(self, text, nuc="mono", strand_specific=False):
        if nuc == "mono":
            if self.background_prob is None:
                background_prob = np.ones(4) / 4
            else:
                background_prob = self.background_prob
            with open(text, 'r') as file:
                data = file.read()
            self.version = re.compile(r'MEME version ([\d+\.*]+)').match(data).group(1)
            self.names = re.findall(r"MOTIF (.*)\n", data)
            self.background = re.findall(
                r"Background letter frequencies.*\n(A .* C .* G .* T .*)\n", data)[0]
            self.strands = re.findall(r"strands: (.*)\n", data)[0].strip()
            self.alphabet = re.findall(r"ALPHABET=(.*)\n", data)[0].strip()
            letter_probs = re.findall(
                r"(letter-probability.*\n([ \t]*\d+\.?\d*[ \t]+\d+\.?\d*"
                r"[ \t]+\d+\.?\d*[ \t]+\d+\.?\d*[ \t]*\n)+)", data)
            assert len(letter_probs) == len(self.names)
            self.nmotifs = len(letter_probs)
            out_channels = self.nmotifs * 2
            in_channels = 4
            matrices = []
            length = 0
            for i in range(len(letter_probs)):
                matrix = letter_probs[i][0].split("\n")
                matrix = matrix[1:-1] if len(matrix[-1]) == 0 else matrix[1:]
                matrices.append(np.array([r.split() for r in matrix], dtype=float))
                if matrices[-1].shape[0] > length:
                    length = matrices[-1].shape[0]

        if nuc == "di":
            if self.background_prob is None:
                background_prob = np.ones(16) / 16
            else:
                background_prob = self.background_prob
            self.names = os.listdir(text)
            self.nmotifs = len(self.names)
            in_channels = 16
            out_channels = self.nmotifs * 2
            matrices = []
            length = 0
            for i in self.names:
                if i.endswith(".dpcm") or i.endswith(".dpwm"):
                    matrix = read_pwm(os.path.join(text, i))
                    matrices.append(matrix)
                    if matrix.shape[0] > length:
                        length = matrix.shape[0]

        out = np.zeros((out_channels, in_channels, length), dtype=np.float32)
        mask = torch.zeros((out_channels, 1, length), dtype=torch.uint8)
        motif_norms = np.zeros(self.nmotifs, dtype=np.float32)
        for k, kernel in enumerate(matrices):
            kernel, motif_norms[k] = transform_kernel(kernel, self.smoothing, background_prob)
            if strand_specific:
                out[2*k,   :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel.T
            else:
                out[2*k,   :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel[::-1, ::-1].T
            mask[2*k,   :, :kernel.shape[0]] = 1
            mask[2*k+1, :, :kernel.shape[0]] = 1
        return torch.from_numpy(out), mask, motif_norms


# ── Data loaders ───────────────────────────────────────────────────────────────

class vcfData:
    """
    Loads variant data from a VCF file and returns reference / alternate
    one-hot tensors for each variant.  Used by MotifDiff.
    """
    def __init__(self, vcf, batchsize, genome, windowsize, dinucleotide=False, strand='+'):
        data = readvcf(vcf)
        self.headers = data.columns.to_list()
        self.strand = strand

        self.ref = data.iloc[:, 3].to_numpy()
        self.alt = data.iloc[:, 4].to_numpy()

        f = np.vectorize(len)
        self.reflength = f(self.ref)
        self.altlength = f(self.alt)

        self.chrs = data.iloc[:, 0].to_numpy()
        self.refstarts = data.iloc[:, 1].to_numpy() - int(windowsize)
        self.refends   = data.iloc[:, 1].to_numpy() + self.reflength - 1 + int(windowsize) - 1
        self.altstarts = data.iloc[:, 1].to_numpy() - int(windowsize)
        self.altends   = data.iloc[:, 1].to_numpy() + self.altlength - 1 + int(windowsize) - 1
        self.pos = data.iloc[:, 1].to_numpy()
        self.variant_names = data.iloc[:, 2].to_numpy()

        self.batchsize = batchsize
        self.n = data.shape[0]
        self.seqs = FastaFile(genome)
        self.windowsize = windowsize
        refs = self.seqs.references
        lengths = self.seqs.lengths
        self.limits = {refs[i]: lengths[i] for i in range(len(refs))}
        self.out = open("coordinatesUsed.bed", "w")
        self.lookup = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.dinucleotide = dinucleotide
        # Store raw sequences for N-masking in MotifDiff
        self._raw_ref_seqs = {}
        self._raw_alt_seqs = {}

    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))

    def names(self):
        return self.variant_names

    def __getitem__(self, i):
        i1, i2 = i * self.batchsize, (i + 1) * self.batchsize
        if i2 >= self.n:
            i2 = self.n
        batchsize = int(i2 - i1)
        targetlength = max(np.max(self.reflength[i1:i2]), np.max(self.altlength[i1:i2]))
        if self.dinucleotide:
            offset = 1
            height = (self.windowsize - 1) * 2 + targetlength - 1
            width = 16
        else:
            offset = 0
            height = (self.windowsize - 1) * 2 + targetlength
            width = 4
        batch    = np.zeros((batchsize, width, height), dtype=np.float32)
        mask     = torch.zeros((batchsize, 1, height), dtype=torch.uint8)
        altbatch = np.zeros((batchsize, width, height), dtype=np.float32)
        altmask  = torch.zeros((batchsize, 1, height), dtype=torch.uint8)
        stats    = np.empty((batchsize, 4))
        raw_ref_seqs = []
        raw_alt_seqs = []

        for j, (c, refs, refe, alts, alte, r, a, lenr, lena) in enumerate(zip(
                self.chrs[i1:i2], self.refstarts[i1:i2], self.refends[i1:i2],
                self.altstarts[i1:i2], self.altends[i1:i2],
                self.ref[i1:i2], self.alt[i1:i2],
                self.reflength[i1:i2], self.altlength[i1:i2])):
            if refs > 0 and refe < self.limits[c]:
                seg = self.seqs.fetch(c, refs, refe)
                if self.strand == '+':
                    seg = seg.upper()
                else:
                    seg = str(Seq(seg).reverse_complement()).upper()
                    r = str(Seq(r).reverse_complement())
                    a = str(Seq(a).reverse_complement())
                assert seg[self.windowsize-1:-(self.windowsize-1)] == r
                ref_seg = seg
                alt_seg = seg[:self.windowsize-1] + a + seg[-(self.windowsize-1):]
            else:
                ref_seg = "N" * height
                alt_seg = "N" * height

            raw_ref_seqs.append(ref_seg)
            raw_alt_seqs.append(alt_seg)
            batch[j,    :, :int(refe-refs-offset)] = returnonehot(ref_seg, self.dinucleotide)
            returnmask(j, mask, self.windowsize, refs, refe, self.dinucleotide)
            altbatch[j, :, :int(alte-alts-offset)] = returnonehot(alt_seg, self.dinucleotide)
            returnmask(j, altmask, self.windowsize, alts, alte, self.dinucleotide)

        return (torch.from_numpy(batch), mask,
                torch.from_numpy(altbatch), altmask,
                raw_ref_seqs, raw_alt_seqs)


class SegmentDataSeq:
    """
    Loads sequences from individual .fa files in a directory.  Used by MotifScore.
    """
    def __init__(self, seq, batchsize, windowsize, up, dinucleotide=False):
        self.files = os.listdir(seq)
        self.id = [file for file in self.files if file.endswith(".fa")]
        self.n = len(self.id)
        self.batchsize = batchsize
        self.additional = 4 * 4 + 2
        self.dinucleotide = dinucleotide
        self.dir = seq

    def names(self):
        return self.id

    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))

    def __getitem__(self, i):
        i1, i2 = i * self.batchsize, (i + 1) * self.batchsize
        if i2 >= self.n:
            i2 = self.n
        batchsize = int(i2 - i1)
        if self.dinucleotide:
            height = 500
            width = 16
        else:
            height = 501
            width = 4
        batch = np.zeros((batchsize, width, height), dtype=np.float32)
        stats = np.empty((batchsize, self.additional), dtype=np.float32)
        raw_seqs = []

        for j in range(i2 - i1):
            if self.id[j + i1].endswith(".fa"):
                sequence = FastaFile(f"{self.dir}/{self.id[j+i1]}")
                for ref in sequence.references:
                    seg = sequence.fetch(ref)
                stats[j] = stringstats(seg)
                raw_seqs.append(seg.upper())
                if self.dinucleotide:
                    batch[j, :, :len(seg)-1] = returnonehot(seg, dinucleotide=True)
                else:
                    batch[j, :, :len(seg)] = returnonehot(seg)
            else:
                print("not a readable fasta file!")
                raw_seqs.append("")
        return torch.from_numpy(batch), stats, raw_seqs


class SegmentDataBed:
    """
    Extracts sequences from a genome FASTA given a BED file of coordinates.
    Used by MotifScore.
    """
    def __init__(self, bed, batchsize, genome, windowsize, up, dinucleotide=False):
        self.chrs, self.starts, self.ends, self.gene = readbed(bed, up)
        self.id = ["_".join([c, str(s), str(e), str(g)])
                   for c, s, e, g in zip(self.chrs, self.starts, self.ends, self.gene)]
        self.midpoints = np.asarray(np.ceil((self.starts + self.ends) / 2), dtype=int)
        self.seqs = FastaFile(genome)
        refs = self.seqs.references
        lengths = self.seqs.lengths
        self.new_starts = self.midpoints - windowsize
        self.new_ends   = self.midpoints + windowsize
        self.batchsize = batchsize
        self.n = len(self.chrs)
        self.padding = windowsize
        self.additional = 4 * 4 + 2
        self.limits = {refs[i]: lengths[i] for i in range(len(refs))}
        self.out = open("coordinatesUsed.bed", "w")
        self.dinucleotide = dinucleotide

    def names(self):
        return self.id

    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))

    def __getitem__(self, i):
        i1, i2 = i * self.batchsize, (i + 1) * self.batchsize
        if i2 >= self.n:
            i2 = self.n
        batchsize = int(i2 - i1)
        if self.dinucleotide:
            height = np.max(self.new_ends[i1:i2] - self.new_starts[i1:i2]) - 1
            width = 16
        else:
            height = np.max(self.new_ends[i1:i2] - self.new_starts[i1:i2])
            width = 4
        batch = np.zeros((batchsize, width, height), dtype=np.float32)
        stats = np.empty((batchsize, self.additional), dtype=np.float32)
        raw_seqs = []

        for j, (c, new_s, new_e) in enumerate(zip(
                self.chrs[i1:i2], self.new_starts[i1:i2], self.new_ends[i1:i2])):
            self.out.write(c + "\t" + str(new_s) + "\t" + str(new_e) + "\n")
            if new_s > 0 and new_e < self.limits[c]:
                seg = self.seqs.fetch(c, new_s, new_e)
            else:
                seg = "N" * (self.padding * 2)
            stats[j] = stringstats(seg)
            raw_seqs.append(seg.upper())
            if self.dinucleotide:
                batch[j, :, :(new_e - new_s) - 1] = returnonehot(seg, dinucleotide=True)
            else:
                batch[j, :, :(new_e - new_s)] = returnonehot(seg)
        return torch.from_numpy(batch), stats, raw_seqs

    def __del__(self):
        pass


if __name__ == "__main__":
    # Quick smoke-test of the N-masking fix
    import torch.nn.functional as F

    # Synthetic 6-bp motif kernel (all-log-prob, max real score < 0)
    kernel = torch.zeros(1, 4, 6)
    kernel[0, 0, :] = -0.1   # slight preference for A at every position

    seq_clean = "AAAAAA"
    seq_n     = "AANNAAA"   # contains N at positions 2,3

    oh_clean = torch.from_numpy(returnonehot(seq_clean)).unsqueeze(0)
    oh_n     = torch.from_numpy(returnonehot(seq_n)).unsqueeze(0)

    # scores shape: (1, 1, n_positions) → keep as (1, n_positions) for apply_invalid_mask
    score_clean = F.conv1d(oh_clean, kernel).squeeze(0)   # shape (1, L)
    score_n     = F.conv1d(oh_n,     kernel).squeeze(0)   # shape (1, L)

    print("Clean scores:", score_clean[0].tolist())
    print("N-seq scores (before mask):", score_n[0].tolist())

    inv = get_invalid_mask(seq_n, kernel_len=6, dinucleotide=False)
    print("Invalid mask:", inv)
    apply_invalid_mask(score_n, inv)   # score_n is (n_kernels, n_positions)
    print("N-seq scores (after mask):", score_n[0].tolist())

    # Verify: no non-inf score should be > 0
    finite_scores = [s for s in score_n[0].tolist() if not np.isinf(s)]
    assert all(s <= 0 for s in finite_scores), \
        f"Spurious score > 0 found: {finite_scores}"
    # Verify: masked positions are -inf
    assert all(np.isinf(score_n[0, i].item()) for i in range(len(inv)) if inv[i]), \
        "Masked positions should be -inf"
    print("Smoke test passed.")

    # Additional test: gap character '-'
    seq_gap = "AAA-AAA"
    oh_gap  = torch.from_numpy(returnonehot(seq_gap)).unsqueeze(0)
    score_gap = F.conv1d(oh_gap, kernel).squeeze(0)
    print("\nGap-seq scores (before mask):", score_gap[0].tolist())
    inv_gap = get_invalid_mask(seq_gap, kernel_len=6, dinucleotide=False)
    print("Gap invalid mask:", inv_gap)
    apply_invalid_mask(score_gap, inv_gap)
    print("Gap-seq scores (after mask):", score_gap[0].tolist())
    finite_gap = [s for s in score_gap[0].tolist() if not np.isinf(s)]
    assert all(s <= 0 for s in finite_gap), f"Spurious gap score > 0: {finite_gap}"
    print("Gap smoke test passed.")
