## BioMM
BioMM: Biological-informed Multi-stage Machine learning framework for phenotype prediction using omics data

## Features of BioMM in a nutshell

1. Applicability for various omics data modalities.     
2. Prioritizing outcome-associated functional patterns.   
3. End-to-end prediction at the individual level based on biological stratified patterns.   
4. Possibility for extension to machine learning models of interest.   
5. Parallel computing. 


## Installation 

BioMM has been incorporated into the [Bioconductor version: Development (3.10)](http://www.bioconductor.org/packages/devel/bioc//html/BioMM.html).
To install this package, start R (version "3.6") and enter:

```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# The following initializes usage of Bioc devel
BiocManager::install(version='devel')
BiocManager::install("BioMM")
``` 

Alternative installation from Github, start R (version "3.5") and enter:

```{r eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("transbioZI/BioMM")
```

## Tutorial 

The detailed instructions on how to use this package are explained in this [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/BioMM/inst/doc/BioMMtutorial.html). 

## Citation

NIPS ML4H submission: Chen, J. and Schwarz, E., 2017. BioMM: Biologically-informed Multi-stage Machine learning for identification of epigenetic fingerprints. arXiv preprint arXiv:1712.00336.
