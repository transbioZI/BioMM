## BioMM
BioMM: Biological-informed Multi-stage Machine learning framework for phenotype prediction using omics data

## Features of BioMM in a nutshell

1. Applicability for various omics data modalities.   
2. Different biological stratification strategies.    
3. Prioritizing outcome-associated functional patterns.   
4. End-to-end prediction at the individual level based on biological stratified patterns.   
5. Possibility for extension to machine learning models of interest.   
6. Parallel computing. 


## Installation 

BioMM has been incorporated into the [Bioconductor](http://www.bioconductor.org/packages/devel/bioc//html/BioMM.html).
To install this package, start R (version "3.6") and enter:

```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BioMM", version = "3.9")
``` 

## Tutorial 

The detailed instructions on how to use this package are explained in this [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/BioMM/inst/doc/BioMMtutorial.html). 

## Citation

NIPS ML4H submission: Chen, J. and Schwarz, E., 2017. BioMM: Biologically-informed Multi-stage Machine learning for identification of epigenetic fingerprints. arXiv preprint arXiv:1712.00336.

