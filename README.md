# MotifDiff
MotifDiff is a novel computational tool designed to quantify variant effects using mono and dinucleotide PWMs. It offers several key advantages, including scalability to score millions of variants within minutes, implementation of various normalization strategies for optimal performance, and support for both dinucleotide and mononucleotide models for transcription factor binding and also considering the effect of variants both inside and outside of the binding site. More details including MotifDiff's efficacy across diverse ground truth datasets and also the mathematical details behind the tool can be found in the paper **Ultra-fast variant effect prediction using biophysical transcription factor binding models**

## Installation
```
pip install MotifDiff-pkg
```
## How to use
The very basic run with the required inputs will be:
```
getDiff --genome hg38/hg19 --motif motif_file --vcf vcf_file --out output_file
```
where,
- ```--genome``` is the FASTA file for hg38 or hg19 depending on the variants.
- ```motif_file```  is either a meme format of PPMs or .pwm or .pcm format files.
- ```vcf_file```  is a vcf file of variants. the order of the columns in the file should be: chr, pos, ID, ref, alt. .csv files are accepted too.
- ```out_file``` is a name/path for the score files to be saved. there will be three outputs ending with .alt which is the scores for the alt variants, .ref which is the scores for ref variants, and .diff which is the diff scores for the variants.
  This will calculate the simple log-odd difference as the scores for all motifs provided in the ```motif_file``` for each variant.

  Following options could be added to get normalized scores from different methods:
- ```--method```  is the method for normalizing the scores from scanning the sequence. could be "FABIAN", or "probNorm".
- ```--nuc```  is to specify whether mono nucleotide models are used or di nucleotide, could be "mono", or "di".
- ```--MaxScale``` is to specify whether to do the Max Scaling or not. Just passing the ```--MaxScale``` is considered True and not passing it is considered as False.
- ```--batch``` is the batch size. default is 128.
- ```--window``` is the window size for the sequence around variants
