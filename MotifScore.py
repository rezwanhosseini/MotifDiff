#!/usr/bin/env python
"""
MotifScore.py — per-sequence PWM binding score computation
===========================================================
N / gap bug fix
---------------
After F.conv1d, any output position whose receptive field overlaps an N or '-'
in the input sequence is set to -inf using apply_invalid_mask().  This prevents
the zero-vector dot-product artefact (score = 0.0 > all real scores) from
propagating through mc_spline / max-pool and producing spurious binding
probabilities of ~1.
"""
import os
from torch.nn import functional as F
import pandas as pd
from util import (SegmentDataBed, SegmentDataSeq,
                  MEME_probNorm, MEME_FABIAN,
                  MCspline_fitting, mc_spline, kmers,
                  get_invalid_mask, apply_invalid_mask)
import torch
import numpy as np
from enum import Enum
import typer
import pickle
import time
import resource
import matplotlib.pyplot as plt
from pysam import FastaFile


app = typer.Typer()

torch.backends.cudnn.deterministic = True


def write_output_motif_features(filename, mat, names, index=None):
    if index is None:
        df = pd.DataFrame(
            mat,
            columns=names + ["GC_ratio", "Masked_ratio"] + [i + "_pattern" for i in kmers()])
        print("writing to csv with no index...")
        df.to_csv(filename, sep="\t", float_format='%.3f', index=False)
    else:
        df = pd.DataFrame(
            mat,
            columns=names + ["GC_ratio", "Masked_ratio"] + [i + "_pattern" for i in kmers()],
            index=index)
        print("writing to csv with index...")
        df.to_csv(filename, sep="\t", float_format='%.3f', index=True)


class mode_type(str, Enum):
    max = "max"
    average = "average"


class kernel_type(str, Enum):
    TFFM = "TFFM"
    PWM = "PWM"


class norm_type(str, Enum):
    iid = "iid"
    motif_based = "motif_based"
    mixture = "mixture"
    from_mono = "from_mono"


class score_type(str, Enum):
    FABIAN = "FABIAN"
    probNorm = "probNorm"
    NONE = "NONE"


class nucleotide_type(str, Enum):
    mono = "mono"
    di = "di"


def _mask_batch(scores, raw_seqs, kernel_len, dinucleotide, motif_len=None):
    """
    Apply N/gap invalid mask to a batch of conv1d scores IN PLACE.

    Parameters
    ----------
    scores : torch.Tensor, shape (B, 2*M, L)
    raw_seqs : list of str, length B
    kernel_len : int   — padded kernel length (kernels.shape[2])
    dinucleotide : bool
    motif_len : int or None — actual motif length for contamination radius
    """
    for b, seq in enumerate(raw_seqs):
        if not seq:
            continue
        inv = get_invalid_mask(seq, kernel_len, dinucleotide, motif_len=motif_len)
        if inv.any():
            apply_invalid_mask(scores[b], inv)


@app.command()
def variantdiff(
        genome: str = typer.Option(None, help="fasta file for the genome"),
        motif_file: str = typer.Option(..., "--motif", help="meme file for the motifs"),
        seqs: str = typer.Option(
            ...,
            help="path to sequences (one folder per gene containing the orthologous sequences "
                 "in fasta format), or a bed file of coordinates to extract the sequence from "
                 "the given genome."),
        diff_score: score_type = typer.Option(
            score_type.NONE, "--method",
            help="how to calculate the diff score (FABIAN/probNorm/NONE)"),
        max_scale: bool = typer.Option(False, "--MaxScale",
                                       help="Apply max transformation or not"),
        nucleotide: nucleotide_type = typer.Option(
            nucleotide_type.mono, "--nuc",
            help="length of the nucleotides in the motifs (mono/di)"),
        normalization_file: str = typer.Option(
            None, "--norm",
            help="file including normalization params. should be consistent with the transform option"),
        up: int = typer.Option(0, help="add upstream"),
        mode: mode_type = typer.Option(mode_type.max,
                                       help="Operation mode for the pooling layer (max/average)"),
        batch: int = typer.Option(128, help="batch size"),
        out_file: str = typer.Option(..., "--out", help="output directory"),
        window: int = typer.Option(500, "--window", help="window size"),
        kernel: kernel_type = typer.Option(
            kernel_type.PWM,
            help="Choose between PWM (4 dimensional) or TFFM (16 dimensional) (default = PWM)."),
        bin: int = typer.Option(1, help="number of bins")
):
    s = time.time()

    kernel = kernel.value
    if kernel == "PWM":
        if diff_score == "FABIAN":
            motif = MEME_FABIAN()
            print(f"Reading the motifs from {motif_file}")
            kernels, kernel_mask, kernel_norms = motif.parse(motif_file, nuc=nucleotide)
        else:
            if normalization_file is None:
                motif = MEME_probNorm()
                kernels, kernel_mask = motif.parse(
                    motif_file, nuc=nucleotide, transform=max_scale)
                if diff_score == "probNorm":
                    print("normalization params not given. normalization process running...")
                    spline_list = MCspline_fitting(kernels, nuc=nucleotide)
            else:
                if "max_scaled" in normalization_file:
                    max_scale = True
                else:
                    max_scale = False
                motif = MEME_probNorm()
                kernels, kernel_mask = motif.parse(
                    motif_file, nuc=nucleotide, transform=max_scale)
                with open(normalization_file, "rb") as fp:
                    spline_list = pickle.load(fp)

    kernel_len   = kernels.shape[2]   # padded kernel length
    dinucleotide = (nucleotide == "di")
    # Maximum actual motif length across all motifs (conservative masking radius)
    motif_lens = []
    for ki in range(0, kernels.shape[0], 2):
        col_norms = kernels[ki].abs().sum(dim=0)
        motif_lens.append(int((col_norms != 0).sum().item()))
    max_motif_len = max(motif_lens) if motif_lens else kernel_len

    if seqs.endswith(".bed"):
        print("positions given.")
        if genome is None:
            raise ValueError(
                "genome is required to extract the sequence. use --genome argument.")
        segments = SegmentDataBed(seqs, batch, genome, int((window + up) / 2), up,
                                  dinucleotide=dinucleotide)
    if os.path.isdir(seqs):
        print("sequences given.")
        segments = SegmentDataSeq(seqs, batch, int((window + up) / 2), up,
                                  dinucleotide=dinucleotide)

    out = np.empty((segments.n, bin * motif.nmotifs + segments.additional))

    print(f"Batch size: {batch}")
    print("Calculating convolutions")

    for i in range(len(segments)):
        print(f"Batch {i+1}:")
        i1, i2 = i * batch, (i + 1) * batch
        if i2 >= segments.n:
            i2 = segments.n

        # SegmentData loaders now return raw_seqs as the 3rd element
        mat, out[i1:i2, bin * motif.nmotifs:], raw_seqs = segments[i]

        tmp = F.conv1d(mat, kernels)

        # ── N / gap masking (bug fix) ──────────────────────────────────────────
        _mask_batch(tmp, raw_seqs, kernel_len, dinucleotide, motif_len=max_motif_len)
        # ──────────────────────────────────────────────────────────────────────

        if diff_score == "NONE":
            tmp = tmp.view(tmp.shape[0], tmp.shape[1] // 2, tmp.shape[2] * 2)
            tmp = F.max_pool1d(tmp, tmp.shape[2])
            tmp = np.squeeze(tmp).numpy()

        if diff_score == "probNorm":
            if mode == "average":
                tmp = tmp.view(tmp.shape[0], tmp.shape[1] // 2, tmp.shape[2] * 2)
                tmp = mc_spline(tmp, spline_list)
                tmp = torch.nan_to_num(tmp, nan=0.0, posinf=0.0, neginf=0.0)
                tmp = F.avg_pool1d(torch.tensor(tmp), tmp.shape[2])
            if mode == "max":
                tmp = tmp.view(tmp.shape[0], tmp.shape[1] // 2, tmp.shape[2] * 2)
                # Replace -inf with a very negative finite value so max-pool works
                tmp = torch.nan_to_num(tmp, nan=-1e9, posinf=-1e9, neginf=-1e9)
                tmp = F.max_pool1d(torch.tensor(tmp), tmp.shape[2])
                tmp = mc_spline(tmp, spline_list)

        if diff_score == "FABIAN":
            tmp = F.max_pool1d(tmp, tmp.shape[2]).numpy()
            tmp = np.max(tmp.reshape(tmp.shape[0], -1, 2), axis=2)

        out[i1:i2, :bin * motif.nmotifs] = tmp.reshape(
            tmp.shape[0], tmp.shape[1] * tmp.shape[2])

    names = motif.names
    new_names = []
    for i in range(bin * len(names)):
        if i % bin > int(bin / 2):
            new_names.append(f'{names[int(i/bin)]} - right - {i%bin}')
        if i % bin == int(bin / 2):
            new_names.append(f'{names[int(i/bin)]} - middle')
        if i % bin < int(bin / 2):
            new_names.append(f'{names[int(i/bin)]} - left - {i%bin}')

    write_output_motif_features(
        out_file + f"_{nucleotide}_{diff_score}_{mode}_{window}",
        out, new_names, segments.names())

    e = time.time()
    print(f"real runtime for {out_file}_{nucleotide}_{diff_score}_{mode}_{window} = ", e - s)


if __name__ == "__main__":
    app()
