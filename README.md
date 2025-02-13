# MotifDiff
MotifDiff is a novel computational tool, designed to quantify variant effects using mono and dinucleotide PWMs. It offers several key advantages, including scalability to score millions of variants within minutes, implementation of various normalization strategies for optimal performance, and support for both dinucleotide and mononucleotide models for transcription factor binding and also considering the effect of variants both inside and outside of the binding site. More details including MotifDiff's efficacy across diverse ground truth datasets and also the mathematical details behind the tool can be found in the paper [**Ultra-fast variant effect prediction using biophysical transcription factor binding models**](https://www.biorxiv.org/content/10.1101/2024.06.26.600873v1)

## Installation
```
pip install MotifDiff-pkg
```
* if you're using Spyder, you might see a dependency version conflict between spyder, pyqt, and pyqtwebengine. we recommend using a compatible version of spyder, otherwise, the error does not interfere with the installation, and the installation will be done regardless.  
## How to use
The basic run with the required inputs will be:
```
getDiff --genome hg38/hg19 --motif motif_file.meme --vcf vcf_file.vcf --out output_file
```
where,
- ```--genome``` is the FASTA file for hg38 or hg19 depending on the variants.
- ```motif_file```  is a file containing the PPMs of motifs and could be any of these formats: .meme, .pwm, .pcm, .pfm, .ppm. our evaluations were based on HOCOMOCO motifs which can be found here https://hocomoco11.autosome.org/downloads_v11 for both mono nucleotides and di nucleotides
- ```vcf_file```  is a vcf file of variants. the order of the columns in the file should be chr, pos, ID, ref, and alt. .csv or .txt files are accepted too.
- ```out_file``` is a name/path for the score files to be saved. there will be three outputs ending with .alt, the scores for the alt variants, .ref, the scores for ref variants, and .diff, the diff scores for the variants.
  
  The following options could be added to get normalized scores from different methods:
- ```--method```  is the method for normalizing the scores from scanning the sequence. could be "FABIAN", or "probNorm". The default is "No normalization" which provides the log-odd difference as the scores.
- ```--mode```  is the mode for pooling the mapped scores after normalization. could be "max", or "average".
- ```--nuc```  is to specify whether mono nucleotide models or di nucleotide models are used, could be "mono", or "di".
- ```--MaxScale``` is to specify whether to do the Max Scaling or not. Just passing the ```--MaxScale``` is considered True and not passing it is considered as False.
- ```--norm``` is a file containing the normalizing parameters, already been prepared for HOCOMOCOv11 Human TFBS PWMs, and can be used here for faster implementation. The parameters for both mono nucleotide and di nucleotide models can be found here: https://www.dropbox.com/scl/fo/sqkltsnhdn2il7olab7vr/h?rlkey=0sj46r3xsq8dbudstezzc1khv&dl=0. it is recommended to use these parameters only if you're using the ``` HOCOMOCOv11_full_HUMAN_mono_meme_format.meme``` file from this repository. if downloading other versions of motifs, using these parameters might lead to errors caused by difference in length of motifs or might not be acurate since the given parameters are speciically calculated based on this version. if this argument is not passed the model will automatically calculate normalizing paramters for the given motifs.
- ```--batch``` is the batch size. default is 128.
- ```--window``` is the window size for the sequence around variants. the default is twice the length of each motif which will only consider the variants inside the binding site. but it can also be set to larger sizes to include variants further away from the binding site.

## test run
```
getDiff --genome hg38.fa --motif HOCOMOCOv11_full_HUMAN_mono_meme_format.meme --norm HOCOMOCOv11_HUMAN_mono_params --nuc mono --vcf vcf_test.csv --method probNorm --mode average --out test
```
