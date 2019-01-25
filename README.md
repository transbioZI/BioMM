# BioMM
BioMM: Biological-informed Multi-stage Machine learning framework for phenotype prediction using omics data


# Getting started  

## Installation 
Development version from Github:

2.) Install BioMM in R  
```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
BiocManager::install("BioMM")
``` 

## Tutorial 

The best view of the tutorial is in HTML format with a table of contents by executing the following R codes after you have downloaded the tutorial (BioMMtutorial.Rmd in vignettes directory).

```{r eval=FALSE}
install.packages("rmarkdown")
library("rmarkdown")
render("BioMMtutorial.Rmd")
```

