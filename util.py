import numpy as np
import torch
import pandas as pd
from pysam import FastaFile
import time
import itertools
import xml.etree.ElementTree as ET
import os
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator
import regex as re


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
        nucleotide_prob = np.ones(4)/4
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


def scoreDistDinuc(pssm, prob, gran=None, size=1000):
    
    nucleotides = ['A', 'C', 'G', 'T']
    nms = [a + b for a in nucleotides for b in nucleotides]
    pssm = pd.DataFrame(pssm)
    pssm.columns = nms
      
    if prob is None:
        #prob = dict(zip(['A', 'C', 'G', 'T'], np.ones(4)/4))
        bg_prob = pd.DataFrame(np.repeat(np.ones(4)/4, pssm.shape[0]).reshape(pssm.shape[0],4), columns=nucleotides)
    else:
        prob = pd.DataFrame(prob, columns=nms)
        bg_prob = pd.DataFrame(np.zeros((prob.shape[0],4)), columns=nucleotides)
        for pos in range(0,prob.shape[0]):
            bg_prob.iloc[pos,:] = [prob.iloc[pos,:][[f'{nuc1}{nuc2}' for nuc2 in nucleotides]].sum() for nuc1 in nucleotides]
    
    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pssm) - np.min(pssm))/(size - 1)            
    
    pssm = np.round(pssm / gran)
    pssm = pssm * gran
    mnscore = np.min(pssm, axis=1).sum()
    mxscore = np.max(pssm, axis=1).sum()
    nscores = int(np.round((mxscore - mnscore) / gran)) + 1
    
    def s2ind(s):
        return int((s / gran - mnscore / gran).round())
    
    SD = np.zeros((4, nscores))
    nucleotides = ['A', 'C', 'G', 'T']
    SD = pd.DataFrame(SD, index=nucleotides)

    pssm_ainds = [i for i, name in enumerate(pssm.columns) if name.endswith('A')]
    pssm_cinds = [i for i, name in enumerate(pssm.columns) if name.endswith('C')]
    pssm_ginds = [i for i, name in enumerate(pssm.columns) if name.endswith('G')]
    pssm_tinds = [i for i, name in enumerate(pssm.columns) if name.endswith('T')]

    ascores = pssm.iloc[0, pssm_ainds].values
    cscores = pssm.iloc[0, pssm_cinds].values
    gscores = pssm.iloc[0, pssm_ginds].values
    tscores = pssm.iloc[0, pssm_tinds].values

    
    ascores_i = {name: s2ind(score) for name, score in zip(pssm.columns[pssm_ainds], ascores)}
    cscores_i = {name: s2ind(score) for name, score in zip(pssm.columns[pssm_cinds], cscores)}
    gscores_i = {name: s2ind(score) for name, score in zip(pssm.columns[pssm_ginds], gscores)}
    tscores_i = {name: s2ind(score) for name, score in zip(pssm.columns[pssm_tinds], tscores)}

    # Probabilities of initial dinucleotides
    def ffu(nms):
        return [np.prod([bg_prob.iloc[0,:][x] for x in nm]) for nm in nms]
        
    aprobs = ffu(ascores_i.keys())
    cprobs = ffu(cscores_i.keys())
    gprobs = ffu(gscores_i.keys())
    tprobs = ffu(tscores_i.keys())
    
    for i in range(4):
        SD.loc['A', ascores_i[list(ascores_i.keys())[i]]] += aprobs[i]
        SD.loc['C', cscores_i[list(cscores_i.keys())[i]]] += cprobs[i]
        SD.loc['G', gscores_i[list(gscores_i.keys())[i]]] += gprobs[i]
        SD.loc['T', tscores_i[list(tscores_i.keys())[i]]] += tprobs[i]

    
    def update_dist(nuc, pssm_inds, pos):
        vals = pssm.iloc[pos,pssm_inds]
        shifts = (vals/gran).round().astype(int) # shouldn't we use s2ind function here too?
        tvec = np.zeros(SD.shape[1])
        for i in range(4):
            tvec += np.roll(SD.iloc[i, :], shifts[i]) * bg_prob.iloc[pos,:][nuc]
        return(tvec)

    for pos in range(1, pssm.shape[0]):
        t1 = update_dist('A', pssm_ainds, pos)
        t2 = update_dist('C', pssm_cinds, pos)
        t3 = update_dist('G', pssm_ginds, pos)
        t4 = update_dist('T', pssm_tinds, pos)
        SD.iloc[0,:] = t1
        SD.iloc[1,:] = t2
        SD.iloc[2,:] = t3
        SD.iloc[3,:] = t4

    #print(mxscore)
    #print(mxscore+gran)
    x = np.arange(mnscore,mxscore,gran)
    y = np.sum(SD, axis=0)
    if len(x)==len(y):
        #Tprint("lengths are equal")
        #print(x[len(x)-3:len(x)])
        return(tuple((x,y)))
    if len(x)<len(y):
        #print("lenagth are NOT equal")
        x = np.arange(mnscore,mxscore+gran-1e-5,gran)
        y = np.sum(SD, axis=0)
        #print(x[len(x)-3:len(x)])
        return(tuple((x,y)))


def return_coef_for_normalization(pwms, nucleotide_prob=None, gran=None, size=1000, nuc="mono"):
    params = []
    for i in range(0,pwms.shape[0],2):
        pwm = pwms[i].numpy().T      
        pwm = pwm[pwm.sum(axis=1) != 0, :]
        nucleotide_prob = np.exp(pwm) / np.sum(np.exp(pwm), axis=1, keepdims=True)
        if nuc=="mono":
            s, d = scoreDist(pwm, nucleotide_prob, gran, size)
        if nuc=="di":
            s, d = scoreDistDinuc(pwm, nucleotide_prob, gran=gran, size=size)
        param, _ = curve_fit(logit, s, np.cumsum(d), maxfev=5000)
        #f = interp1d(np.exp(s), np.cumsum(d))
        #print(curve_fit(logit, np.exp(s), np.cumsum(d), maxfev=5000))
        #params.append(param)
        params.append(param)
    return params

def MCspline_fitting(pwms, nucleotide_prob=None, gran=None, size=1000, nuc="mono"):
    spline_list = []
    for i in range(0,pwms.shape[0],2):
        pwm = pwms[i].numpy().T      
        pwm = pwm[pwm.sum(axis=1) != 0, :]
        nucleotide_prob = np.exp(pwm) / np.sum(np.exp(pwm), axis=1, keepdims=True)
        if nuc=="mono":
            s, d = scoreDist(pwm, nucleotide_prob, gran, size)
        if nuc=="di":
            s, d = scoreDistDinuc(pwm, nucleotide_prob, gran=gran, size=size)  
        spl = PchipInterpolator(s, np.cumsum(d))
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

def normalize_mat(mat, params):
    out = torch.empty_like(mat)
    assert mat.shape[1] == len(params)
    for i in range(len(params)):
        #out[:,i] = logit(mat[:,i], *params[i])
        #tmp = np.clip(mat[:,i],params[i].x.min(), params[i].x.max())
        #tmp = params[i](tmp)
        out[:,i] = logit_torch(mat[:,i], *params[i])
    return out

def mc_spline (mat, spline_list):
    out = torch.empty_like(mat)
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
        data = pd.read_csv(filename, skiprows=nh, header=None, sep="\t")
        #data.columns = pd.MultiIndex.from_tuples([tuple(i[1:] for i in data.columns[0])] +list(data.columns)[1:])
    elif nh == 1:
        data = pd.read_csv(filename, skiprows=1, header=None, sep="\t")
        #data.columns = [data.columns[0][1:]] + data.columns.to_list()[1:]
    else:
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

def read_TFFM(file):
    tree = ET.parse(file)
    root = tree.getroot()
    data = []
    for state in root[0].iterfind("state"):
        discrete = state[0]
        if "order" in discrete.attrib:
            data.append(discrete.text.split(","))
    return np.array(data, dtype=float)

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

    def parse(self, text, nuc="mono", transform=False):
        if nuc == "mono":  
            if self.background_prob is None:
                background_prob = np.ones(4)/4
            else:
                background_prob = self.background
            
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
            else:
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
                    matrices.append(np.array([i.split() for i in matrix], dtype=float))
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
            
            if transform:
                kernel, _ = transform_kernel(kernel, self.smoothing, background_prob)
            else:    
                if np.min(kernel)<0:
                    #print( "it's already the log likelihood, no need to do the log transform")
                    kernel = kernel 
                else:
                    kernel[kernel == 0] = self.precision
                    kernel = np.log(kernel)
            
            
            out[2*k  , :, :kernel.shape[0]] = kernel.T
            out[2*k+1, :, :kernel.shape[0]] = kernel[::-1, ::-1].T
            mask[2*k  , :, :kernel.shape[0]] = 1
            mask[2*k+1, :, :kernel.shape[0]] = 1
        
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
            
    def parse(self, text, nuc="mono"):
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
            out[2*k  , :, :kernel.shape[0]] = kernel.T
            out[2*k+1, :, :kernel.shape[0]] = kernel[::-1, ::-1].T
            mask[2*k  , :, :kernel.shape[0]] = 1
            mask[2*k+1, :, :kernel.shape[0]] = 1
        return torch.from_numpy(out), mask, motif_norms

class TFFM():
    def __init__(self):
        self.names = []
        self.nmotifs = 0

    def parse(self, directory):
        self.names = os.listdir(directory)
        self.nmotifs = len(self.names)
        in_channels = 16
        out_channels = self.nmotifs
        data = []
        height = 0
        for i in self.names:
            tffm = read_TFFM(os.path.join(directory, i))
            data.append(tffm)
            if tffm.shape[0] > height:
                height = tffm.shape[0]
        out = np.zeros((out_channels, in_channels, height), dtype=np.float32)
        mask = torch.zeros((out_channels, 1 , height), dtype=torch.uint8)
        for n, tffm in enumerate(data):
            out[n, :, :tffm.shape[0]] = tffm.T
            mask[n, :, :tffm.shape[0]] = 1
        return torch.from_numpy(out), mask

class TFFM_with_Transformation():
    def __init__(self, precision=1e-7, smoothing=0.02, background=None):
        self.names = []
        self.nmotifs = 0
        self.precision=1e-7
        self.smoothing = smoothing
        self.background = []
        if background is None:
            self.background_prob = np.ones(16)*0.0625
        else:
            self.background_prob = background
    def parse(self, directory):
        self.names = os.listdir(directory)
        self.nmotifs = len(self.names)
        in_channels = 16
        out_channels = self.nmotifs * 2
        data = []
        height = 0
        for i in self.names:
            if i.endswith(".dpcm") or i.endswith(".dpwm"):
                tffm = read_pwm(os.path.join(directory, i))
                data.append(tffm)
                if tffm.shape[0]>height:
                    height = tffm.shape[0]               
            else:
                tffm = read_TFFM(os.path.join(directory, i))
                data.append(tffm)
                if tffm.shape[0] > height:
                    height = tffm.shape[0]
        #print(data)
        out = np.zeros((out_channels, in_channels, height), dtype=np.float32)
        mask = torch.zeros((out_channels, 1 , height), dtype=torch.uint8)
        motif_norms = np.zeros(self.nmotifs, dtype=np.float32)
        for n, tffm in enumerate(data):
            tffm, motif_norms[n] = transform_kernel(tffm, self.smoothing, self.background_prob)
            out[2*n  , :, :tffm.shape[0]] = tffm.T
            out[2*n+1, :, :tffm.shape[0]] = tffm[::-1, ::-1].T
            mask[2*n , :, :tffm.shape[0]] = 1
            mask[2*n+1,:, :tffm.shape[0]] = 1
        return torch.from_numpy(out), mask, motif_norms



class vcfData:
    def __init__(self, vcf, batchsize, genome, windowsize, dinucleotide = False):
        data = readvcf(vcf)
        self.headers = data.columns.to_list()
        
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
                seg=seg.upper()
                #print(f"Sequence: {seg[:self.windowsize-1]} {seg[self.windowsize-1]} {seg[self.windowsize:]}")
                #print(f"a: ({a}, {self.lookup[a]}), r: ({r}, {self.lookup[r]}), Target: {seg[self.windowsize-1]}")
                #assert(seg[self.windowsize-1]==r or len(a)!=1 or len(r)!=1)
                #print(i, c, refs+int(self.windowsize), seg[self.windowsize-1:-(self.windowsize-1)], r)
                assert(seg[self.windowsize-1:-(self.windowsize-1)]==r)
                batch[i, :, :int(refe-refs-offset)] = returnonehot(seg, self.dinucleotide)
                returnmask(i, mask, self.windowsize, refs, refe, self.dinucleotide)
                #print(f"{seg[:self.windowsize-1]} + {a} + {seg[-(self.windowsize-1):]}, {self.dinucleotide}")
                altbatch[i, :, :int(alte-alts-offset)] = returnonehot(seg[:self.windowsize-1] + a + seg[-(self.windowsize-1):], self.dinucleotide)
                returnmask(i, altmask, self.windowsize, alts, alte, self.dinucleotide)
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

class SegmentData:
    def __init__(self, bed, batchsize, genome, windowsize, up, dinucleotide=False):
        self.chrs, self.starts, self.ends, self.peaks = readbed(bed, up)
        self.id = ["_".join([c, str(s), str(e)]) for c, s, e in zip(self.chrs, self.starts, self.ends)]
        self.midpoints = np.asarray(np.ceil((self.starts + self.ends)/2),dtype=int)
        self.seqs = FastaFile(genome)
        refs = self.seqs.references
        lengths = self.seqs.lengths
        if windowsize>(min(lengths)/2):
            self.new_starts = self.starts
            self.new_ends=self.ends
        else:
            self.new_starts = self.midpoints - windowsize
            self.new_ends = self.midpoints + windowsize
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
        i1, i2 = i*self.batchsize, (i+1)*self.batchsize
        if i2 >= self.n: i2 = self.n
        batchsize = int(i2 - i1)
        if self.dinucleotide:
            height = np.max(self.new_ends[i1:i2] - self.new_starts[i1:i2])-1# + self.padding
            width = 16
        else:
            height = np.max(self.new_ends[i1:i2] - self.new_starts[i1:i2])# + self.padding
            width = 4
        batch = np.zeros((batchsize, width, height), dtype=np.float32) 
        stats = np.empty((batchsize, self.additional), dtype=np.float32)
        for i, c, p, s, e, new_s, new_e in zip(range(i2-i1), self.chrs[i1:i2], self.peaks[i1:i2], self.starts[i1:i2], self.ends[i1:i2], self.new_starts[i1:i2], self.new_ends[i1:i2]):
            self.out.write(c+"\t"+str(new_s)+"\t"+str(new_e)+"\n")
            if all(self.peaks!=None) and all('peak' in string for string in self.peaks):
                if i==0: print('peaks available')
                seg = self.seqs.fetch(p, new_s-s, new_e-s)
            else:
                if i==0: print('peaks not available')
                if new_s>0 and new_e<self.limits[c]:
                    seg = self.seqs.fetch(c, new_s, new_e)
                else:
                    seg = "N"*(self.padding*2)

            stats[i] = stringstats(seg)
            if self.dinucleotide:
                batch[i, :, :(new_e-new_s)-1] = returnonehot(seg, dinucleotide=True)
            else:
                batch[i, :, :(new_e-new_s)] = returnonehot(seg)
        return torch.from_numpy(batch), stats

    def __del__(self):
        pass
        #self.out.close()

if __name__ == "__main__":
    motif = MEME_FABIAN()
    kernels = motif.parse("TestInput/Test.meme", "none")
    segments = vcfData("TestInput/TestHg38.vcf", 128, "/data/genomes/human/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa", kernels.shape[2])
    print(segments.headers)
    start = time.time()
    for i in range(len(segments)):
        orig, alt = segments[i]
    end = time.time()
    print(f"other took {end-start}")
