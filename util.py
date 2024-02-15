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
import itertools as itt
import diNucMat as dnma
import diNucMotDist as dnsd

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
    pssm = dnma.diNucMat(tmp, cn)

    if gran is None:
        if size is None:
            raise ValueError("provide either gran or size. Both missing.")
        gran = (np.max(pwm) - np.min(pwm))/(size - 1)    
            
    tmp  = np.exp(tmp)
    tmp  = tmp / np.sum(tmp, axis=1, keepdims=True)   

    prob = dnma.diNucProbMat(tmp, cn)
    # calculating the score dist
    sd_mot = dnsd.score_dist(pssm, prob, gran=0.01)

    #- now let's make a dinucleotide model with iid columns
    avg_dinuc_freqs = prob.values.mean(axis=0)
    iid_prob_values = np.repeat(avg_dinuc_freqs,pwm.shape[0]).reshape((pwm.shape[0],16), order = 'F')
    prob_bg1        = dnma.diNucProbMat(iid_prob_values, cn)
    sd_bg1          = dnsd.score_dist(pssm, prob_bg1, gran=0.01)
    
    return(sd_mot.x, sd_bg1.y, sd_mot.y)


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
            #if k==144: print(torch.from_numpy(out[2*k  , :, :kernel.shape[0]]))
        return torch.from_numpy(out), mask
    
    def names(self):
        return self.names
    

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
