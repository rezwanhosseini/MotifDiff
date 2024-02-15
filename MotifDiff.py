#!/usr/bin/env python
from torch.nn import functional as F
import pandas as pd
from utilNEW import vcfData, MEME_probNorm, MEME_FABIAN, MCspline_fitting, mc_spline #
import torch
import numpy as np
from enum import Enum
import typer
import pickle
import time
import resource
import matplotlib.pyplot as plt
#import pyarrow.feather as feather


app = typer.Typer()

torch.backends.cudnn.deterministic = True

def write_output_diff(filename, mat, names, index=None):
    if index is None:
        pd.DataFrame(mat, columns = names).to_csv(filename, header=True, index=False, sep="\t")
        #feather.write_feather(pd.DataFrame(mat, columns = names), f"{filename}.feather")
    else:
        pd.DataFrame(mat, columns = names, index=index).to_csv(filename, header=True, index=True, sep="\t")
        #feather.write_feather(pd.DataFrame(mat, columns = names, index=index), f"{filename}.feather")


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
    
@app.command()
def variantdiff(genome: str = typer.Option(..., help="fasta file for the genome"),
                motif_file: str = typer.Option(..., "--motif", help="meme file for the motifs"),
                vcf: str = typer.Option(..., help="vcf file"), 
                diff_score: score_type = typer.Option(score_type.NONE, "--method", help="how to calculate the diff score (FABIAN/probNorm/NONE)"),
                max_scale: bool = typer.Option(False, "--MaxScale", help="Apply max transformation or not"),
                nucleotide: nucleotide_type = typer.Option(nucleotide_type.mono, "--nuc", help="length of the nucleotides in the motifs (mono/di)"),
                norm_method: norm_type = typer.Option(None, "--methodNorm", help="what method to use for the normalization (iid/motifbased/mixture)"),
                normalization_file: str = typer.Option(None, "--norm", help="file including normalization params. should be consistent with the transform option"),
                mode: mode_type = typer.Option(mode_type.max, help="Operation mode for the pooling layer (max/average)"),                 
                batch: int = typer.Option(128, help="batch size"),
                out_file: str = typer.Option(..., "--out", help="output directory"),
                window:int = typer.Option(0, help="window size"), # change this in a way that if it gets a value use that value and if not the default will be kernel size (instead of setting 0 as the default, make it none or not receiving an input as the default)
                kernel: kernel_type = typer.Option(kernel_type.PWM, help="Choose between PWM (4 dimensional) or TFFM (16 dimensional) (default = PWM).")
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
                kernels, kernel_mask = motif.parse(motif_file, nuc=nucleotide, transform=max_scale)
                if diff_score=="probNorm":
                    print("normalization params not given. normalization process running...")
                    spline_list = MCspline_fitting(kernels, nuc=nucleotide)

                    
            else:
                if "max_scaled" in normalization_file:
                    max_scale=True
                else:
                    max_scale=False
                motif = MEME_probNorm()
                kernels, kernel_mask = motif.parse(motif_file, nuc=nucleotide, transform=max_scale)
                with open(normalization_file, "rb") as fp:
                    spline_list = pickle.load(fp)
        
        if window==0:
            windowsize = kernels.shape[2]
        else:
            windowsize=window


    
    segments = vcfData(vcf, batch, genome, windowsize, dinucleotide=(nucleotide=="di"))
    
    outRef = np.empty((segments.n, motif.nmotifs), dtype=np.float32)
    outAlt = np.empty((segments.n, motif.nmotifs), dtype=np.float32)
    alpha = 0.1
    print(f"Batch size: {batch}")
    print("Calculating convolutions")
    for i in range(1):#len(segments)):
        print(f"Batch {i+1}:")
        i1, i2 = i*batch, (i+1)*batch
        if i2 >= segments.n: i2 = segments.n
        matref, maskref, matalt, maskalt = segments[i]

        ref = F.conv1d(matref, kernels)
        alt = F.conv1d(matalt, kernels)
        
        if window==0:
            motif_mask = F.conv1d(maskref, kernel_mask)
            ref[motif_mask == 0] = -torch.inf
            motif_mask = F.conv1d(maskalt, kernel_mask)        
            alt[motif_mask == 0] = -torch.inf
        
        if diff_score == "NONE":
            ref = ref.view(ref.shape[0], ref.shape[1]//2, ref.shape[2]*2)
            alt = alt.view(alt.shape[0], alt.shape[1]//2, alt.shape[2]*2)
            ref = F.max_pool1d(ref, ref.shape[2])
            alt = F.max_pool1d(alt, alt.shape[2])
            ref = np.squeeze(ref).numpy()
            alt = np.squeeze(alt).numpy()
            
        if diff_score == "probNorm":
            if mode == "average":
                ref = ref.view(ref.shape[0], ref.shape[1]//2, ref.shape[2]*2)
                alt = alt.view(alt.shape[0], alt.shape[1]//2, alt.shape[2]*2)
                ref = mc_spline(ref, spline_list)
                alt = mc_spline(alt, spline_list)           
                ref = np.nan_to_num(ref, nan=0)
                alt = np.nan_to_num(alt, nan=0)            
                ref = F.avg_pool1d(torch.tensor(ref), ref.shape[2])
                alt = F.avg_pool1d(torch.tensor(alt), alt.shape[2])   
            if mode == "max":
                ref = ref.view(ref.shape[0], ref.shape[1]//2, ref.shape[2]*2)
                alt = alt.view(alt.shape[0], alt.shape[1]//2, alt.shape[2]*2)
                ref = np.nan_to_num(ref, nan=0)
                alt = np.nan_to_num(alt, nan=0)
                ref = F.max_pool1d(torch.tensor(ref), ref.shape[2])
                alt = F.max_pool1d(torch.tensor(alt), alt.shape[2])
                ref = mc_spline(ref, spline_list)
                alt = mc_spline(alt, spline_list)
            ref = np.squeeze(ref).numpy()
            alt = np.squeeze(alt).numpy()
            
        if diff_score == "FABIAN":                 
            ref = F.max_pool1d(ref, ref.shape[2]).numpy()
            alt = F.max_pool1d(alt, alt.shape[2]).numpy()
            ref = np.max(ref.reshape(ref.shape[0],-1,2), axis=2) # separates the convolutions from the original kernel and the reverse complement kernel into two differetn columns AND then keeps the maximum between those two
            alt = np.max(alt.reshape(alt.shape[0],-1,2), axis=2)
        
        outRef[i1:i2, :motif.nmotifs] = ref 
        outAlt[i1:i2, :motif.nmotifs] = alt

    print("calculating diff...")
    if diff_score == "FABIAN":
        outRef = 1 - outRef/kernel_norms[np.newaxis,:]
        outAlt = 1 - outAlt/kernel_norms[np.newaxis,:]
        mask = outRef > outAlt
        f = np.empty_like(outRef)
        f[np.where(mask)] = -(1-outAlt[mask]+alpha)/(1-outRef[mask]+alpha)+1
        mask = outAlt >= outRef
        f[np.where(mask)] = (1-outRef[mask]+alpha)/(1-outAlt[mask]+alpha)-1
        outDiff = 2/(1+np.power(2, -2*f)) - 1
        
    
    if diff_score == "probNorm" or diff_score == "NONE":
        outDiff = outAlt-outRef
    
    motif_names = motif.names
    
    
    print(f"Writing the results to {out_file}_{nucleotide}_{diff_score}_{mode}")
    write_output_diff(out_file+f"_{nucleotide}_{diff_score}_{mode}.alt", outAlt, motif_names, segments.names())
    write_output_diff(out_file+f"_{nucleotide}_{diff_score}_{mode}.ref", outRef, motif_names, segments.names())
    write_output_diff(out_file+f"_{nucleotide}_{diff_score}_{mode}.diff", outDiff, motif_names, segments.names())
    
    e = time.time()
    print(f"real runtime for {out_file}_{nucleotide}_{diff_score}_{mode} = ", e-s)


if __name__ == "__main__":
    app()
