## BioMM
BioMM: Biological-informed Multi-stage Machine learning framework for phenotype prediction using omics data

## Features of BioMM in a nutshell

1. Applicability for various omics data modalities (e.g. methylome, transcriptomics, genomics).   
2. Various biological stratification strategies.    
3. Prioritizing outcome-associated functional patterns.   
4. End-to-end prediction at the individual level based on biological stratified patterns.   
4. Possibility for an extension to machine learning models of interest.   
6. Parallel computing. 

## Installation 

BioMM installation from Github
```{r eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("transbioZI/BioMM", build_vignettes=TRUE)
``` 

BioMM has been incorporated into the [Bioconductor](http://www.bioconductor.org/packages/devel/bioc//html/BioMM.html).
To install this package from BioConductor, start R (version "4.0") and enter:

```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BioMM", version = "4.0")
``` 

## Tutorial 

The detailed instructions on how to use this package are explained in the most updated [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/BioMM/inst/doc/BioMMtutorial.html). 

## Citation

NIPS ML4H submission: Chen, J. and Schwarz, E., 2017. BioMM: Biologically-informed Multi-stage Machine learning for identification of epigenetic fingerprints. arXiv preprint arXiv:1712.00336.

Chen, Junfang, et al. "Association of a Reproducible Epigenetic Risk Profile for Schizophrenia With Brain Methylation and Function." JAMA psychiatry (2020).