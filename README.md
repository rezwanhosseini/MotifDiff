# MotifDiff
MotifDiff is a novel computational tool, designed to quantify variant effects using mono and dinucleotide PWMs. It offers several key advantages, including scalability to score millions of variants within minutes, implementation of various normalization strategies for optimal performance, and support for both dinucleotide and mononucleotide models for transcription factor binding and also considering the effect of variants both inside and outside of the binding site. More details including MotifDiff's efficacy across diverse ground truth datasets and also the mathematical details behind the tool can be found in the paper [**Ultra-fast variant effect prediction using biophysical transcription factor binding models**](https://www.biorxiv.org/content/10.1101/2024.06.26.600873v2)

## Installation
motifDiff can be both installed by pip or cloning the repository.

### pip installation
```
pip install MotifDiff-pkg
```
### cloning the repository
```
git clone ...
```
then navigate to the directory and create the virtual environment to install the required dependencies.
```
cd MotifDiff/
conda env create -f env.yml
conda activate motifdiff_env
```

* if you're using Spyder, you might see a dependency version conflict between spyder, pyqt, and pyqtwebengine. we recommend using a compatible version of spyder, otherwise, the error does not interfere with the installation, and the installation will be done regardless.  
## How to use 
### pip
If installed by pip, motifDiff can be run by the following command:
```
getDiff --genome <path-to-reference-genome> --motif <path-to-motif-file> --vcf <path-to-vcf> --out <output-name>
```
### cloned
If installed by cloning the repository, motifDiff can be run by the following command:
```
python3 MotifDiff.py --genome <path-to-reference-genome> --motif <path-to-motif-file> --vcf <path-to-vcf> --out <output-name>
```
where,
- ```--genome``` takes the path to the FASTA file for the reference genome, which could be either hg38.fa or hg19.fa depending on your variants.
- ```--motif```  takes the path to the file containing the motif models. motif models can be of the following formats: .meme, .pwm, .pcm, .pfm, .ppm. If the format is .meme, .pfm, or .ppm, it is expected that all the motifs are in a single file. If the formate is .pwm, or .pcm, it is expected that each motif is in a separate file and the input to this argument must be the directory containing those files.
our evaluations in the [paper](https://www.biorxiv.org/content/10.1101/2024.06.26.600873v2) were based on HOCOMOCO file in this repository. any other versions of HOCOMOCO motifs can lead to different results. 
- ```--vcf```  takes the path to a file containing the variants. This input doesn't have to be in the proper vcf format. it only needs the following columns in the mentioned order: chr, pos, ID, ref, and alt. Look at TestHg38.vcf in the repository for an example.
- ```--out``` takes the name you'd want your output to be saved as. if not passed, "output" will be used as default. There will be multiple outputs saved in the directory including:
  a file ending with .ref, containing the binding scores of given motifs on the reference (WT) sequence,
  a file ending with .alt, containing the binding scores of given motifs on the alternative (MT) sequence,
  a file endig with .diff, containg the difference in the binding score on the sequence before and after the variant. if <0, the variant is considered destruptive, decreasing the binding probability. if >0, the variant is considered constructive, increasing the binding probability.
  two additional files will also be saved ending with .refMaxPos and .altMaxPos which contain the location on the sequence (ref, alt respectively) where it matches the given motif the most. 

### more options to use:
- ```--method```  by default is "probNorm" and will use probNorm normalization to calculate the effect. Other options are "FABIAN", which re-implements the heuristic normalization of [FABIAN-variant](https://pmc.ncbi.nlm.nih.gov/articles/PMC9252790/) to calculate the variant effects, or "No normalization", which provides the raw log-odd difference of the ref and alt sequence as scores of matching the binding, without any normalization.
- ```--mode```  by default is "max" and is the mode for pooling the mapped scores after the convolution. "max" will capture the best site (larger effects) on the sequence that can be bound by a TF and is recommended when looking at more specific or distinguished binding sites. "average" captures the average binding affinity of the whole sequence (accouting for smaller effects as well) and is recommended when looking at more general binding site.
- ```--nuc```  by default is "mono" and is to specify whether what nucleotide models are used, could be "mono", or "di". This option has to be compatible with the given motif model to --motif.
- ```--MaxScale``` by default is True and scales the PWM values up to the max value at each position to increase the highest matching value when scanning the sequence by each motif. "False" can used if no scaling is desired.
- ```--batch``` is the batch size and by default is 128, can be set to higher or lower values depending on the resources.
- ```--window``` by default is twice the length of longest given motif containing the variant in the center of the sequence. Can be set to larger values centered around variants to capture the effect of variants further away from the binding site.
- ```--strand``` specifies the strand on which the TF is expected to bind. if "+"/"-", it only looks at the +/- strand and will calculate the effect of the variant on the binding of the TF on that strand. if not passed (default), it will caluclate the effect on both strands and takes the max value.

