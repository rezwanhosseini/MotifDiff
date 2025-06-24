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

torch.set_printoptions(precision=8)
np.set_printoptions(precision=8)


def number_of_headers(filename):
    header=0
    with open(filename,"r") as file:      
        while True:
            line = file.readline()      
            if line.startswith("#"):
                header=header+1
            else:
                break
    return header

def kmers_count(seq, k=2):
    lookup = {"".join(i):0 for i in itertools.product(["A","C","G","T"], repeat=k)}
    mers = [seq[i:i+2] for i in range(len(seq)-k+1)]
    for i in mers:
        if i in lookup:
            lookup[i] += 1
    for i in lookup:
        lookup[i] /= (len(seq)-k+1)
    return list(lookup.values())

def kmers(k=2):
    return ["".join(i) for i in itertools.product(["A","C","G","T"], repeat=k)]

def logit(x, a, b):
    return 1/(1 + np.exp(-a * x - b))

def logit_torch(x, a, b):
    return 1/(1 + torch.exp(-a * x - b))

def init_dist(dmin, dmax, dp, weights, probs):
    out = np.zeros(int(np.round((dmax-dmin)/dp)+1))
    ii = np.array(np.round((weights-dmin)/dp), dtype=int)
    for i in range(len(probs)):
        out[ii[i]] = out[ii[i]] + probs[i]
    return out

def scoreDist(pwm, nucleotide_prob=None, gran=None, size=1000):
    if nucleotide_prob is None:
        nucleotide_prob = [np.ones(4)/4]*pwm.shape[0]
    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pwm) - np.min(pwm))/(size - 1)
    pwm = np.round(pwm/gran)*gran
    pwm_max, pwm_min = pwm.max(axis=1), pwm.min(axis=1)
    distribution = init_dist(pwm_min[0], pwm_max[0], gran, pwm[0], nucleotide_prob[0])   
    for i in range(1, pwm.shape[0]):
        kernel = init_dist(pwm_min[i], pwm_max[i], gran, pwm[i], nucleotide_prob[i])
        distribution = np.convolve(distribution, kernel)
    support_min = pwm_min.sum()
    ii = np.where(distribution > 0)[0]
    support = support_min + (ii) * gran
    return support, distribution[ii]

class diNucMat:
    def __init__(self, values, colnames) -> None:
        #- check if colnames is a 1-D array with 16 elements, each a dinucleotide in the correct order
        assert colnames == ["".join(i) for i in itt.product(["A","C","G","T"], repeat=2)]
        self.values   = values
        self.colnames = colnames      
    @property
    def values(self):
        return self._values    
    @values.setter
    def values(self, values):
        #- check if probs is a 2-D array with 16 columns
        assert len(values.shape) == 2
        assert values.shape[1] == 16
        assert values.shape[0] > 0
        self._values = values
        #- make the column matrices, rely on correct column order
        #  therefore need the assert in init...
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
#==============================================================================
    #- assert pssm  is a ov class diNucMat and 
    #  prob is of class diNucProbMat
    assert isinstance(pssm, diNucMat)    
    assert isinstance(prob, diNucProbMat)

    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pssm) - np.min(pssm))/(size - 1)  
        
    # utility function for converting a score to an index (of the distribution)
    def vals2inds(vals, mnscore, gran):
        return np.rint(((vals - mnscore)/gran)).astype(int)
    
    # discretization and score range
    mnscore = np.floor(np.sum(np.min(pssm.values, axis=1))) #- lower bound
    mxscore = np.ceil(np.sum(np.max(pssm.values, axis=1))) #- upper bound
    if mxscore<0: mxscore=0
    nscores = int(np.rint((mxscore - mnscore) / gran) + 1) #- shouldn't really need to round
    
    # make the score distribution
    #- INITIALIZATION
    SD = np.zeros((4, nscores))
    #- FIRST POSITION
    for i in range(4):
        #- update SD at the right indices
        for j in range(4):
            SD[i, vals2inds(pssm.colMats[0, j, i], mnscore, gran)] += prob.colMats[0, j, i]

    #- need a copy of SD for updating
    SD_tmp = SD.copy()
    #- ITERATE through rest of the motif
    for pos in range(1, pssm.values.shape[0]):
        #- nuc is ending nucleotide (second of the two)
        for nuc in range(4):
            tvec   = np.zeros(nscores)
            scores = pssm.colMats[pos, :, nuc]
            shifts = np.rint(scores / gran).astype(int)
            #- i is the starting nucleotide
            for i in range(4):
                tvec += np.roll(SD_tmp[i, :], shifts[i]) * prob.trnMats[pos, i, nuc] 
            #- overwriteing SD here, so thats why we need to use SD_tmp above
            SD[nuc, :] = tvec
        #- need to update SD_tmp for the next iteration
        SD_tmp = SD.copy()

    # RETURN final motif distribution

    x = np.arange(mnscore, mxscore + gran, gran)
    y = np.sum(SD, axis=0)
    ii = np.where(y > 0)[0]
    support = mnscore + (ii) * gran
    
    x = support
    y = y[ii]
    
    #- results as a named tuple
    Dist = namedtuple('Dist', ['x', 'y'])
    return Dist(x, y)



def mono2di(ppm): 
    num_rows = ppm.shape[0]
    num_cols = ppm.shape[1]** 2
    ppm_di = np.zeros((num_rows-1, num_cols))
    for i in range(num_rows-1):
        for j in range(ppm.shape[1]):
            ppm_di[i,4*j:4*j+4] = ppm[i,j]*ppm[i+1,:]
    return ppm_di    

def scoreDistDinuc(pwm, gran=None, size=1000):
    tmp  = pwm 
    cn   = ["".join(i) for i in itt.product(["A","C","G","T"], repeat=2)]
    pssm = diNucMat(tmp, cn)

    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pwm) - np.min(pwm))/(size - 1)    
            
    tmp  = np.exp(tmp)
    tmp  = tmp / np.sum(tmp, axis=1, keepdims=True)   

    prob = diNucProbMat(tmp, cn)
    # calculating the score dist
    sd_mot = diNucMotDist(pssm, prob, gran=0.01)

    #- now let's make a dinucleotide model with iid columns
    avg_dinuc_freqs = prob.values.mean(axis=0)
    iid_prob_values = np.repeat(avg_dinuc_freqs,pwm.shape[0]).reshape((pwm.shape[0],16), order = 'F')
    prob_bg1        = diNucProbMat(iid_prob_values, cn)
    sd_bg1          = diNucMotDist(pssm, prob_bg1, gran=0.01)
    
    return(sd_mot.x, sd_bg1.y, sd_mot.y)



#def return_coef_for_normalization(pwms, nucleotide_prob=None, gran=None, size=1000, nuc="mono"):
#    params = []
#    for i in range(0,pwms.shape[0],2):
#        pwm = pwms[i].numpy().T      
#        pwm = pwm[pwm.sum(axis=1) != 0, :]
#        nucleotide_prob = np.exp(pwm) / np.sum(np.exp(pwm), axis=1, keepdims=True)
#       if nuc=="mono":
#           s, d = scoreDist(pwm, nucleotide_prob, gran, size)
#        if nuc=="di":
#            s, d = scoreDistDinuc(pwm, nucleotide_prob, gran=gran, size=size)
#        param, _ = curve_fit(logit, s, np.cumsum(d), maxfev=5000)
#        #f = interp1d(np.exp(s), np.cumsum(d))
#        #print(curve_fit(logit, np.exp(s), np.cumsum(d), maxfev=5000))
#        #params.append(param)
#        params.append(param)
#    return params

def MCspline_fitting(pwms, nucleotide_prob=None, gran=None, size=1000, nuc="mono", method="motif_based"):
    spline_list = []
    for i in range(0,pwms.shape[0],2):
        pwm = pwms[i].numpy().T      
        pwm = pwm[pwm.sum(axis=1) != 0, :]
        nucleotide_prob = np.exp(pwm) / np.sum(np.exp(pwm), axis=1, keepdims=True)
        if nuc=="mono":
            s, d = scoreDist(pwm, nucleotide_prob, gran, size)
            spl = PchipInterpolator(s, np.cumsum(d))
        if nuc=="di":
            s, d_iid, d_m = scoreDistDinuc(pwm, gran=gran, size=size)  
            if method=="iid":
                spl = PchipInterpolator(s, np.cumsum(d_iid))
            if method=="motif_based":
                spl = PchipInterpolator(s, np.cumsum(d_m))
            if method=="mixture":
                spl = PchipInterpolator(s, np.cumsum(0.5*d_m + 0.25*d_iid + 0.25*1/len(d_iid)))
        spline_list.append(spl)
    return spline_list

#def return_coef_for_normalization_diff(pwms, nucleotide_prob=None, gran=None, size=1000, length_correction=1):
#   params = []
#    for i in range(0,pwms.shape[0],2):
#        pwm = pwms[i].numpy().T
#        pwm = pwm[pwm.sum(axis=1) != 0, :]
#        #prob = pwm.sum(axis=0)/pwm.sum()
#       prob = np.sum(np.exp(pwm) / np.exp(pwm).sum(axis=1).reshape(-1,1), axis=0)/np.sum(np.exp(pwm) / np.exp(pwm).sum(axis=1).reshape(-1,1))
#       s, d = scoreDist(pwm, prob, gran, size)#, diff=True)
#       param, _ = curve_fit(logit, s, np.power(np.cumsum(d), length_correction))
#       params.append(param)
#   return params

#def normalize_mat(mat, params):
#    out = torch.empty_like(mat)
#    assert mat.shape[1] == len(params)
#   for i in range(len(params)):
#       #out[:,i] = logit(mat[:,i], *params[i])
#       #tmp = np.clip(mat[:,i],params[i].x.min(), params[i].x.max())
#       #tmp = params[i](tmp)
#       out[:,i] = logit_torch(mat[:,i], *params[i])
#   return out

def mc_spline (mat, spline_list):
    out = torch.empty_like(mat)
    #print(mat.shape)
    #print(len(spline_list))
    assert mat.shape[1] == len(spline_list)
    for i in range(len(spline_list)):
        spl = spline_list[i]
        out_i = spl(mat[:,i])
        out_i[out_i>1]=1
        out_i[out_i<0]=0
        out[:,i] = torch.tensor(out_i)
    return out

#def readvcf(filename):
#    nh = number_of_headers(filename)
#    if nh > 1:
#        data = pd.read_csv(filename, header=list(range(nh)), sep="\t")
#        data.columns = pd.MultiIndex.from_tuples([tuple(i[1:] for i in data.columns[0])] +list(data.columns)[1:])
#    elif nh == 1:
#        data = pd.read_csv(filename, header=0, sep="\t")
#        data.columns = [data.columns[0][1:]] + data.columns.to_list()[1:]
#    else:
#        data = pd.read_csv(filename, header=None, sep="\t")
#    return data  

def readvcf(filename):
    nh = number_of_headers(filename)
    if nh > 1:
        #print(nh, " headers in the vcf file.")
        data = pd.read_csv(filename, skiprows=nh, header=None, sep="\t")
        #data.columns = pd.MultiIndex.from_tuples([tuple(i[1:] for i in data.columns[0])] +list(data.columns)[1:])
    elif nh == 1:
        data = pd.read_csv(filename, skiprows=1, header=None, sep="\t")
        #data.columns = [data.columns[0][1:]] + data.columns.to_list()[1:]
    else:
        #print("no header")
        data = pd.read_csv(filename, header=None, sep="\t")
    return data  

def readbed(filename, up):
    data = pd.read_csv(filename, sep = "\t", header = None)
    chrs = data[0].to_numpy()
    start = data[1].to_numpy(dtype=int)
    end = data[2].to_numpy(dtype=int)
    if (data.shape[1]>3):
        peaks = data[3].to_numpy(dtype=str)
        if(data.shape[1] > 5): #get the strand   
            print("Strand detected")
            up = int(np.floor(up))
            strand = data[5].to_numpy()
            #adjust the regions to acccount for strand and up
            start = start - (strand == "+") * up #[start[i]-up if strand[i]=="+" else start[i] for i in range(len(start))]
            end = end + (strand == "-") * up #[end[i]+up if strand[i]=="-" else end[i] for i in range(len(start))]
    else:
        peaks=np.array([None]*len(chrs))
    return chrs, start, end, peaks

def returnmask(i, mask, windowsize, start, end, dinucleotide):
    if dinucleotide:
        tmp = np.zeros(mask.shape[2]+1)
        tmp[int(windowsize-1):int(end-start-windowsize+1)] = 1
        mask[i,:,:] = torch.from_numpy(np.convolve(tmp, [1,1], mode="valid"))
    else:
        mask[i, :, int(windowsize-1):int(end-start-windowsize+1)] = 1


def returnonehot(string, dinucleotide=False):
    string = string.upper()
    tmp = np.array(list(string))

    if dinucleotide:
        lookup = {"".join(i):n for n,i in enumerate(itertools.product(["A","C","G","T"], repeat=2))}
        icol = np.where(tmp == 'N')[0]
        #icol = np.unique(icol // 2)
        #icol = np.where(np.logical_not(np.isin(np.arange(len(tmp)//2), icol)))[0]
        icol = np.unique(np.clip(np.concatenate([icol, icol-1]), 0, len(tmp)-2))
        icol = np.where(np.logical_not(np.isin(np.arange(len(tmp)-1), icol)))[0]
        tmp = np.array([tmp[i] + tmp[i+1] for i in range(len(tmp)-1)])
        irow = np.array([lookup[i] for i in tmp[icol]])
    else:
        lookup = {'A':0, 'C':1, 'G':2, 'T':3}
        icol = np.where(tmp != 'N')[0]
        irow = np.array([lookup[i] for i in tmp[icol]])

    out = np.zeros((len(lookup),len(tmp)), dtype = np.float32)

    if len(icol)>0:
        out[irow,icol] = 1

    return out

#def read_TFFM(file):
#    tree = ET.parse(file)
#    root = tree.getroot()
#    data = []
#    for state in root[0].iterfind("state"):
#        discrete = state[0]
#       if "order" in discrete.attrib:
#           data.append(discrete.text.split(","))
#   return np.array(data, dtype=float)

def read_pwm(filename):
    with open(filename,'r') as file:
        lines = file.readlines()
    values = []
    for line in lines:
        if not line.startswith(">"):
            values.append(line.strip().split("\t"))
    values = np.array(values, dtype=float)
    if np.min(values)>=0:
        values = values/values.sum(axis=1, keepdims=True)
    return np.array(values, dtype=float)

def transform_kernel(kernel, smoothing, background):
    if np.min(kernel)<0: #(if kernels are already log transformed)
        out=kernel
    else: 
        out = np.log(kernel / background + smoothing)
    c = out.max(axis=1)
    out = out - c[:, np.newaxis]
    norm = out.min(axis=1).sum()
    return out, norm

class MEME_probNorm():
    def __init__(self, precision=1e-7, smoothing=0.02, background=None):
        self.version = 0
        self.alphabet = ""
        self.strands = ""
        #self.headers = []
        self.background = []
        self.names = []
        self.nmotifs = 0
        self.precision=1e-7
        self.smoothing = smoothing
        self.background_prob = background

    def parse(self, text, nuc="mono", transform=False, strand_specific=False):
        if nuc == "mono":  
            if self.background_prob is None:
                background_prob = np.ones(4)/4
            else:
                background_prob = self.background
            
            if text.endswith(".pfm") or text.endswith(".ppm"):
                print("motif is pfm or ppm format")
                with open(text,'r') as file:
                    data = file.read()
                self.names = re.findall(r"(>GM\.5\.0\.\S+)", data)
                self.synonames = re.findall(r"(#GM\.5\.0\.\S+)", data)
                letter_probs = re.findall(r"(>GM.*\n((?:[ \t]*\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]*\n)+))", data)
                assert len(letter_probs) == len(self.names)
                self.nmotifs = len(self.names)
                out_channels = self.nmotifs * 2
                in_channels = 4
                matrices = []
                length = 0
                for i in range(len(letter_probs)):
                    matrix = letter_probs[i][0].split("\n")
                    if len(matrix[-1]) == 0:
                        matrix = matrix[1:-1]
                    else:
                        matrix = matrix[1:]
                    matrices.append(np.array([i.split() for i in matrix], dtype=float))
                    if matrices[-1].shape[0] > length:
                        length = matrices[-1].shape[0]
                        
            if text.endswith(".meme"):
                with open(text,'r') as file:
                    data = file.read()
                self.version = re.compile(r'MEME version ([\d+\.*]+)').match(data).group(1)
                self.names = re.findall(r"MOTIF (.*)\n", data)
                self.background = re.findall(r"Background letter frequencies.*\n(A .* C .* G .* T .*)\n", data)[0]
                self.strands = re.findall(r"strands: (.*)\n", data)[0].strip()
                self.alphabet = re.findall(r"ALPHABET=(.*)\n", data)[0].strip()
                letter_probs = re.findall(r"letter-probability.*\n((?:[ \t]*\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]+\d*\.?\d+[eE]?-?\d*[ \t]*\n)+)", data)
                assert len(letter_probs) == len(self.names)
                self.nmotifs = len(letter_probs)
                out_channels = self.nmotifs * 2
                in_channels = 4
                matrices = []
                length = 0
                for i in range(len(letter_probs)):
                    matrix = letter_probs[i].split("\n")
                    if len(matrix[-1]) == 0:
                        matrix = matrix[1:-1]
                    else:
                        matrix = matrix[1:]
                    matrices.append(np.array([i.split() for i in matrix], dtype=float))
                    if matrices[-1].shape[0] > length:
                        length = matrices[-1].shape[0]
            
            if os.path.isdir(text):
                self.names = os.listdir(text)
                self.nmotifs = len(self.names)
                in_channels = 4
                out_channels = self.nmotifs * 2
                matrices = []
                length = 0
                for k,i in enumerate(self.names):
                    if i.endswith(".pcm") or i.endswith(".pwm"):
                        matrix = read_pwm(os.path.join(text, i))
                        matrices.append(matrix)
                        if matrix.shape[0]>length:
                            length = matrix.shape[0] 
            
        if nuc == "di":
            if self.background_prob is None:
                background_prob = np.ones(16)/16
            else:
                background_prob = self.background_prob
            
            if text.endswith(".meme"):
                with open(text,'r') as file:
                    data = file.read()
                self.version = re.compile(r'MEME version ([\d+\.*]+)').match(data).group(1)
                self.names = re.findall(r"MOTIF (.*)\n", data)
                self.background = re.findall(r"Background letter frequencies.*\n(A .* C .* G .* T .*)\n", data)[0]
                self.strands = re.findall(r"strands: (.*)\n", data)[0].strip()
                self.alphabet = re.findall(r"ALPHABET=(.*)\n", data)[0].strip()
                letter_probs = re.findall(r"(letter-probability.*\n([ \t]*\d+\.?\d*[ \t]+\d+\.?\d*[ \t]+\d+\.?\d*[ \t]+\d+\.?\d*[ \t]*\n)+)", data)
                assert len(letter_probs) == len(self.names)
                self.nmotifs = len(letter_probs)
                out_channels = self.nmotifs * 2
                in_channels = 16
                matrices = []
                length = 0
                for i in range(len(letter_probs)):
                    matrix = letter_probs[i][0].split("\n")
                    if len(matrix[-1]) == 0:
                        matrix = matrix[1:-1]
                    else:
                        matrix = matrix[1:]
                    m = np.array([i.split() for i in matrix], dtype=float)
                    if m.shape[1]==4:
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
                for k,i in enumerate(self.names):
                    if i.endswith(".dpcm") or i.endswith(".dpwm"):
                        matrix = read_pwm(os.path.join(text, i))
                        matrices.append(matrix)
                        if matrix.shape[0]>length:
                            length = matrix.shape[0]              
        
        out = np.zeros((out_channels, in_channels, length), dtype=np.float32)
        mask = torch.zeros((out_channels, 1, length), dtype=torch.uint8)
        for k, kernel in enumerate(matrices):
            #if transform == "constant":
            #    bg=np.repeat(0.25, in_channels).reshape(1,4)
            #if transform == "local":
            #    bg=np.average(kernel,0).reshape(1,4)
            #if transform != "none":
            #   offset=np.min(kernel[kernel>0])
            #    bgMat=np.tile(bg,(kernel.shape[0],1))
            #    kernel=np.log((kernel+offset)/bgMat)
            
            #if k==144: print("real pwm", kernel)
            if transform:
                kernel, _ = transform_kernel(kernel, self.smoothing, background_prob)
            else:    
                if np.min(kernel)<0:
                    #print( "it's already the log likelihood, no need to do the log transform")
                    kernel = kernel 
                else:
                    kernel[kernel == 0] = self.precision
                    kernel = np.log(kernel)
            #if k==144: print("after transformation transpose", kernel.T)
            
            if strand_specific:
                out[2*k  , :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel.T
                mask[2*k  , :, :kernel.shape[0]] = 1
                mask[2*k+1, :, :kernel.shape[0]] = 1
            else:
                out[2*k  , :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel[::-1, ::-1].T
                mask[2*k  , :, :kernel.shape[0]] = 1
                mask[2*k+1, :, :kernel.shape[0]] = 1
            #if k==144: print(torch.from_numpy(out[2*k  , :, :kernel.shape[0]]))
        return torch.from_numpy(out), mask
    
    def names(self):
        return self.names
    
    #def Names(self, text):
    #    if text.endswith(".meme"):
    #        with open(text,'r') as file:
    #            data = file.read()
    #        names = re.findall(r"MOTIF (.*)\n", data)
    #    else:
    #        names = os.listdir(text)
    #    return names

class MEME_FABIAN():
    def __init__(self, precision=1e-7, smoothing=0.02, background=None):
        self.version = 0
        self.alphabet = ""
        self.strands = ""
        #self.headers = []
        self.background = []
        self.names = []
        self.nmotifs = 0
        self.precision=1e-7
        self.smoothing = smoothing
        self.background_prob = background
            
    def parse(self, text, nuc="mono", strand_specific=False):
        if nuc == "mono":
            if self.background_prob is None:
                background_prob = np.ones(4)/4
            else:
                background_prob = self.background_prob   
            with open(text,'r') as file:
                data = file.read()
            self.version = re.compile(r'MEME version ([\d+\.*]+)').match(data).group(1)
            self.names = re.findall(r"MOTIF (.*)\n", data)
            self.background = re.findall(r"Background letter frequencies.*\n(A .* C .* G .* T .*)\n", data)[0]
            self.strands = re.findall(r"strands: (.*)\n", data)[0].strip()
            self.alphabet = re.findall(r"ALPHABET=(.*)\n", data)[0].strip()
            letter_probs = re.findall(r"(letter-probability.*\n([ \t]*\d+\.?\d*[ \t]+\d+\.?\d*[ \t]+\d+\.?\d*[ \t]+\d+\.?\d*[ \t]*\n)+)", data)
            assert len(letter_probs) == len(self.names)
            self.nmotifs = len(letter_probs)
            out_channels = self.nmotifs * 2
            in_channels = 4
            matrices = []
            length = 0
            for i in range(len(letter_probs)):
                matrix = letter_probs[i][0].split("\n")
                if len(matrix[-1]) == 0:
                    matrix = matrix[1:-1]
                else:
                    matrix = matrix[1:]
                matrices.append(np.array([i.split() for i in matrix], dtype=float))
                if matrices[-1].shape[0] > length:
                    length = matrices[-1].shape[0]
        
        if nuc == "di":
            if self.background_prob is None:
                background_prob = np.ones(16)/16
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
                    if matrix.shape[0]>length:
                        length = matrix.shape[0]             
        out = np.zeros((out_channels, in_channels, length), dtype=np.float32)
        mask = torch.zeros((out_channels, 1, length), dtype=torch.uint8)
        motif_norms = np.zeros(self.nmotifs, dtype=np.float32)
        for k, kernel in enumerate(matrices):
            kernel, motif_norms[k] = transform_kernel(kernel, self.smoothing, background_prob)
            if strand_specific:
                out[2*k  , :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel.T
                mask[2*k  , :, :kernel.shape[0]] = 1
                mask[2*k+1, :, :kernel.shape[0]] = 1
            else:
                out[2*k  , :, :kernel.shape[0]] = kernel.T
                out[2*k+1, :, :kernel.shape[0]] = kernel[::-1, ::-1].T
                mask[2*k  , :, :kernel.shape[0]] = 1
                mask[2*k+1, :, :kernel.shape[0]] = 1
        return torch.from_numpy(out), mask, motif_norms

#class TFFM():
#    def __init__(self):
#       self.names = []
#       self.nmotifs = 0
#
#    def parse(self, directory):
#        self.names = os.listdir(directory)
#        self.nmotifs = len(self.names)
#       in_channels = 16
#        out_channels = self.nmotifs
#       data = []
#        height = 0
#        for i in self.names:
#            tffm = read_TFFM(os.path.join(directory, i))
#            data.append(tffm)
#            if tffm.shape[0] > height:
#                height = tffm.shape[0]
#        out = np.zeros((out_channels, in_channels, height), dtype=np.float32)
#        mask = torch.zeros((out_channels, 1 , height), dtype=torch.uint8)
#        for n, tffm in enumerate(data):
#           out[n, :, :tffm.shape[0]] = tffm.T
#            mask[n, :, :tffm.shape[0]] = 1
#        return torch.from_numpy(out), mask

#class TFFM_with_Transformation():
#    def __init__(self, precision=1e-7, smoothing=0.02, background=None):
#        self.names = []
#        self.nmotifs = 0
#        self.precision=1e-7
#        self.smoothing = smoothing
#        self.background = []
#        if background is None:
#            self.background_prob = np.ones(16)*0.0625
#        else:
#            self.background_prob = background
#    def parse(self, directory):
#        self.names = os.listdir(directory)
#        self.nmotifs = len(self.names)
#        in_channels = 16
#        out_channels = self.nmotifs * 2
#        data = []
#        height = 0
#        for i in self.names:
#            if i.endswith(".dpcm") or i.endswith(".dpwm"):
#                tffm = read_pwm(os.path.join(directory, i))
#                data.append(tffm)
#                if tffm.shape[0]>height:
#                    height = tffm.shape[0]               
#            else:
#                tffm = read_TFFM(os.path.join(directory, i))
#                data.append(tffm)
#                if tffm.shape[0] > height:
#                    height = tffm.shape[0]
#        #print(data)
#        out = np.zeros((out_channels, in_channels, height), dtype=np.float32)
#        mask = torch.zeros((out_channels, 1 , height), dtype=torch.uint8)
#        motif_norms = np.zeros(self.nmotifs, dtype=np.float32)
#        for n, tffm in enumerate(data):
#            tffm, motif_norms[n] = transform_kernel(tffm, self.smoothing, self.background_prob)
#            out[2*n  , :, :tffm.shape[0]] = tffm.T
#            out[2*n+1, :, :tffm.shape[0]] = tffm[::-1, ::-1].T
#            mask[2*n , :, :tffm.shape[0]] = 1
#            mask[2*n+1,:, :tffm.shape[0]] = 1
#        return torch.from_numpy(out), mask, motif_norms



class vcfData:
    def __init__(self, vcf, batchsize, genome, windowsize, dinucleotide = False, strand='+'):
        data = readvcf(vcf)
        #print(data)
        #print(data.shape)
        self.headers = data.columns.to_list()

        self.strand = strand
        
        self.ref = data.iloc[:,3].to_numpy()
        self.alt = data.iloc[:,4].to_numpy()

        f = np.vectorize(len)

        self.reflength = f(self.ref)
        self.altlength = f(self.alt)

        self.chrs = data.iloc[:,0].to_numpy()

        self.refstarts = data.iloc[:,1].to_numpy() - int(windowsize)
        self.refends = data.iloc[:,1].to_numpy() + self.reflength - 1 + int(windowsize) - 1

        self.altstarts = data.iloc[:,1].to_numpy() - int(windowsize)
        self.altends = data.iloc[:,1].to_numpy() + self.altlength - 1 + int(windowsize) - 1

        self.pos = data.iloc[:,1].to_numpy()

        self.variant_names = data.iloc[:, 2].to_numpy()

        self.batchsize = batchsize
        self.n = data.shape[0] 
        self.seqs = FastaFile(genome)
        self.windowsize = windowsize
        refs = self.seqs.references
        lengths = self.seqs.lengths
        self.limits = {refs[i]: lengths[i] for i in range(len(refs))}
        self.out = open("coordinatesUsed.bed", "w")
        self.lookup = {'A':0, 'C':1, 'G':2, 'T':3}
        self.dinucleotide = dinucleotide
        
    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))

    def names(self):
        return self.variant_names

    def __getitem__(self, i):
        i1, i2 = i*self.batchsize, (i+1)*self.batchsize
        if i2 >= self.n: i2 = self.n
        batchsize = int(i2 - i1)
        targetlength = max(np.max(self.reflength[i1:i2]), np.max(self.altlength[i1:i2]))
        if self.dinucleotide:
            offset = 1
            height = (self.windowsize-1)*2 + targetlength - 1 #np.max(self.ends[i1:i2] - self.starts[i1:i2])# + self.padding
            width = 16 
        else:
            offset=0
            height = (self.windowsize-1)*2 + targetlength #np.max(self.ends[i1:i2] - self.starts[i1:i2])# + self.padding
            width = 4
        batch = np.zeros((batchsize, width, height), dtype=np.float32) 
        mask = torch.zeros((batchsize, 1, height), dtype=torch.uint8)
        altbatch = np.zeros((batchsize, width, height), dtype=np.float32) 
        altmask = torch.zeros((batchsize, 1, height), dtype=torch.uint8)
        stats = np.empty((batchsize, 4))
        for i, c, refs, refe, alts, alte, r, a, lenr, lena in zip(range(i2-i1), self.chrs[i1:i2], self.refstarts[i1:i2], self.refends[i1:i2], self.altstarts[i1:i2], self.altends[i1:i2], self.ref[i1:i2], self.alt[i1:i2], self.reflength[i1:i2], self.altlength[i1:i2]):
            if refs>0 and refe<self.limits[c]:
                seg = self.seqs.fetch(c, refs, refe)
                #print(seg)
                #print(self.strand)
                if self.strand=='+':
                    seg=seg
                    seg=seg.upper()
                else:
                    seg = str(Seq(seg).reverse_complement())
                    seg=seg.upper()
                    revcomp_r = str(Seq(r).reverse_complement())
                    revcomp_a = str(Seq(a).reverse_complement())
                    r = revcomp_r
                    a = revcomp_a
                #print(seg)
                #print(seg[self.windowsize-1:-(self.windowsize-1)], r)
                assert(seg[self.windowsize-1:-(self.windowsize-1)]==r)

                
                #if i==0: 
                #    if self.dinucleotide:
                #        #seg="CTGCATAAACCGTCGGCAACGTTGGCCACCAGGGGGCGCCATGCACGTGGG"
                #        seg="GCATAAACCGTCGGCAACGTTGGCCACCAGGGGGCGCCATGCACGTG"
                #    else:
                #        seg="GGCAGCAGAGGGAATGCAGATGGCCACCAGGGGGCGCCAGGAGTCAGCA"
                #   print(seg, len(seg))
                #    print(self.windowsize-1, -(self.windowsize-1), seg[self.windowsize-1:-(self.windowsize-1)], r)
                #    print(r, a)
                #print(f"Sequence: {seg[:self.windowsize-1]} {seg[self.windowsize-1]} {seg[self.windowsize:]}")
                #print(f"a: ({a}, {self.lookup[a]}), r: ({r}, {self.lookup[r]}), Target: {seg[self.windowsize-1]}")
                #assert(seg[self.windowsize-1]==r or len(a)!=1 or len(r)!=1)
                #print(i, c, refs+int(self.windowsize), seg[self.windowsize-1:-(self.windowsize-1)], r)
                
                
                batch[i, :, :int(refe-refs-offset)] = returnonehot(seg, self.dinucleotide)
                returnmask(i, mask, self.windowsize, refs, refe, self.dinucleotide)
                #print(f"{seg[:self.windowsize-1]} + {a} + {seg[-(self.windowsize-1):]}, {self.dinucleotide}")
                altbatch[i, :, :int(alte-alts-offset)] = returnonehot(seg[:self.windowsize-1] + a + seg[-(self.windowsize-1):], self.dinucleotide)
                returnmask(i, altmask, self.windowsize, alts, alte, self.dinucleotide)
                #if i in range(2):
                    #print(a, revcomp_a, '--', r, revcomp_r)
                    #print('ref onehot:', seg)
                    #print('alt onehot:', seg[:self.windowsize-1] + a + seg[-(self.windowsize-1):])
                    #print(i, batch[i, :, :int(refe-refs-offset)]-altbatch[i, :, :int(alte-alts-offset)])
        return torch.from_numpy(batch), mask, torch.from_numpy(altbatch), altmask #torch.from_numpy(batch)

def countlowercase(arr):
    return sum([1 for c in arr if c.islower()])

def stringstats(string):
    lowercaseratio = countlowercase(string)/len(string)
    string = string.upper()
    tmp = np.array(list(string))
    gccount = np.sum(np.logical_or(tmp == 'C', tmp == 'G'))/len(tmp)
    #gcpattern = string.count("GC")/(len(tmp)-1)
    #cgpattern = string.count("CG")/(len(tmp)-1)
    patterns = kmers_count(string)
    return np.array([gccount, lowercaseratio, *patterns], dtype=np.float32)

#class SegmentData:
#    def __init__(self, bed, batchsize, genome, windowsize, up, dinucleotide=False):
#        self.chrs, self.starts, self.ends, self.peaks = readbed(bed, up)
#        self.id = ["_".join([c, str(s), str(e)]) for c, s, e in zip(self.chrs, self.starts, self.ends)]
#        self.midpoints = np.asarray(np.ceil((self.starts + self.ends)/2),dtype=int)
#        self.seqs = FastaFile(genome)
#        refs = self.seqs.references
#        lengths = self.seqs.lengths
#        if windowsize>(min(lengths)/2):
#            self.new_starts = self.starts
#            self.new_ends=self.ends
#        else:
#            self.new_starts = self.midpoints - windowsize
#            self.new_ends = self.midpoints + windowsize
#        self.batchsize = batchsize
#        self.n = len(self.chrs)
#        self.padding = windowsize
#        self.additional = 4 * 4 + 2
#        self.limits = {refs[i]: lengths[i] for i in range(len(refs))}
#        self.out = open("coordinatesUsed.bed", "w")
#        self.dinucleotide = dinucleotide

#    def names(self):
#        return self.id

#    def __len__(self):
#        return int(np.ceil(self.n / self.batchsize))

#    def __getitem__(self, i):
#        i1, i2 = i*self.batchsize, (i+1)*self.batchsize
#        if i2 >= self.n: i2 = self.n
#        batchsize = int(i2 - i1)
#        if self.dinucleotide:
#            height = np.max(self.new_ends[i1:i2] - self.new_starts[i1:i2])-1# + self.padding
#            width = 16
#        else:
#            height = np.max(self.new_ends[i1:i2] - self.new_starts[i1:i2])# + self.padding
#            width = 4
#        batch = np.zeros((batchsize, width, height), dtype=np.float32) 
#        stats = np.empty((batchsize, self.additional), dtype=np.float32)
#        for i, c, p, s, e, new_s, new_e in zip(range(i2-i1), self.chrs[i1:i2], self.peaks[i1:i2], self.starts[i1:i2], self.ends[i1:i2], self.new_starts[i1:i2], self.new_ends[i1:i2]):
#            self.out.write(c+"\t"+str(new_s)+"\t"+str(new_e)+"\n")
#            if all(self.peaks!=None) and all('peak' in string for string in self.peaks):
#                if i==0: print('peaks available')
#                seg = self.seqs.fetch(p, new_s-s, new_e-s)
#            else:
#                if i==0: print('peaks not available')
#                if new_s>0 and new_e<self.limits[c]:
#                    seg = self.seqs.fetch(c, new_s, new_e)
#                else:
#                    seg = "N"*(self.padding*2)

#            stats[i] = stringstats(seg)
#            if self.dinucleotide:
#                batch[i, :, :(new_e-new_s)-1] = returnonehot(seg, dinucleotide=True)
#            else:
#                batch[i, :, :(new_e-new_s)] = returnonehot(seg)
#        return torch.from_numpy(batch), stats

#    def __del__(self):
#        pass
#        #self.out.close()

if __name__ == "__main__":
    motif = MEME_FABIAN()
    kernels = motif.parse("TestInput/Test.meme", "none")
    segments = vcfData("TestInput/TestHg38.vcf", 128, "/data/genomes/human/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa", kernels.shape[2])
    print(segments.headers)
    start = time.time()
    for i in range(len(segments)):
        orig, alt = segments[i]
    end = time.time()
    print(f"test took {end-start}")
