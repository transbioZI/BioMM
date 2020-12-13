## BioMM
BioMM: Biological-informed Multi-stage Machine learning framework for phenotype prediction using omics data


# Overview
## Motivation
Identifying reproducible and interpretable biological patterns from high-dimensional omics data is a critical factor in understanding the risk mechanism of complex disease. As such, explainable machine learning can offer biological insight in addition to personalized risk scoring.

## Deliverables
We have implemented a biologically informed multi-stage machine learning framework termed __BioMM__ [1] specifically for phenotype prediction using omics-scale data based on prior biological information including gene ontological (GO) and/or KEGG pathways.   

**Features of BioMM in a nutshell**:   

1. Applicability for multiple omics data modalities (e.g. methylome, genomics). 
2. Various biological stratification strategies.    
3. Prioritizing outcome-associated functional patterns.   
4. Personalized scoring based on biological stratified patterns.   
5. Possibility for an extension to learning models of interest.   

## Installation 

BioMM installation from Github
```{r eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("transbioZI/BioMM", build_vignettes=FALSE)
``` 

* Load required libraries
```{r loadPkg, eval=TRUE, results="hide"}
library(BioMM)
library(BiocParallel) 
library(ranger)
library(rms)
library(glmnet)
library(e1071)
library(precrec)
library(vioplot)
library(CMplot)
library(imager)
library(topGO)
library(xlsx)
```

## DNA methylation data
For a better understanding of the BioMM framework, we used one small examplar datasets: one genome-wide DNA methylation data consisting of 40 subjects and 26486 CpGs for demonstration.  

```{r studyData, eval=TRUE}
## DNA methylation data 
methylData <- readRDS(system.file("extdata", "/methylData.rds", package="BioMM"))
# The first column is the label, and the rest are the features (e.g., CpGs) 
head(methylData[,1:4])
# 0: control and 1: patient
table(methylData[,1]) 
dim(methylData)

``` 
## Feature mapping
Features like CpGs, genes or SNPs can be mapped into pathways based on genomic location and pathway annotation, as implemented in the function `omics2pathlist()`. The examples of pathway databases are gene ontology (GO), KEGG and Reactome. Gene ontological and KEGG pathways are used in this tutorial.

```{r annotationFile, eval=TRUE}
## Load feature annotation data
featureAnno <- readRDS(system.file("extdata", "cpgAnno.rds", package="BioMM"))
# The mapping between CpGs and genes (i.e. entrezID or gene symbol)
head(featureAnno)
# total number of CpGs under investigation
str(unique(featureAnno[,1]))

## Reprocessed Gene ontological pathway annotation with 10 and 200 genes for each pathway
golist <- readRDS(system.file("extdata", "goDB.rds", package="BioMM")) 
## Number of investigated biological processes
length(golist)
str(golist[1:3])

## Reprocessed KEGG pathway annotation with 10 and 200 genes for each pathway
kegglist <- readRDS(system.file("extdata", "keggDB.rds", package="BioMM"))  
## Number of investigated KEGG pathways 
length(kegglist)
str(kegglist[1:3]) 

``` 

To annotate pathways, we demonstrate the usage of `omics2pathlist()` function based on two different pathway databases and two data modalities as follows.

```{r pathlist, eval=TRUE} 
## Feature annotation to pathways 
## Use 100 pathways to reduce the runtime for the downstream analysis. 
numPath <- 100

# GO pathway mapping using DNA methylation data
golistSub <- golist[seq_len(numPath)]
methylGOlist <- omics2pathlist(data=methylData, pathlistDB=golistSub, 
                               featureAnno=featureAnno, 
                               restrictUp=200, restrictDown=10, minPathSize=10) 
# KEGG pathway mapping using DNA methylation data
kegglistSub <- kegglist[seq_len(numPath)]
methylKEGGlist <- omics2pathlist(data=methylData, pathlistDB=kegglistSub, 
                                 featureAnno=featureAnno, 
                                 restrictUp=200, restrictDown=10, minPathSize=10)

``` 

# BioMM framework
## Recapitulation
The BioMM model framework mainly consists of two learning stages [1]. During the first stage, the stage-1 model (either supervised or unsupervised learning model) is aimed to learn the 'latent variables' (i.e., stage-2 features) based on the original dataset incorporating biological meta-information with a resampling strategy. In the second stage, a supervised model (stage-2 model) is built using the stage-2 data with non-negative outcome-associated features for final prediction.

## Choice of model parameters
The end-to-end prediction is performed using `BioMM()` fit on the training dataset to provide predictions on the test data set. `pathlistDB` indicates the type of the stratification of predictors using biological information. Both supervised and unsupervised learning are implemented. `predMode` indicates the prediction type. Generic resampling methods defined using the argument `resample1="CV"` or `resample1="BS"` are used for reconstructing stage-2 features. For more details regarding parameter options, please check `BioMM()` in the manual. The 'biological' parameters such as choice of pathway databases, gene length within a pathway, and targeted CpGs should be reasonably prespecified, with human involvement. In terms of other parameters such as the number of features at both stages, the classifier, and the hyperparameters within each classifier can be informed using cross-validation or bootstrapping procedure.  

#### An example of disease outcome prediction
To apply BioMM with the Random Forest model, we use the argument `supervisedStage1=TRUE` and `classifier=randForest` in `BioMM()`. DNA methylation data mapping to GO pathways are used.

```{r BioMMrandForest4methylGO, eval=TRUE}
## random split the data into training and test
set.seed(1)
trainIndex <- sample(nrow(methylData), nrow(methylData)*0.6)
trainData <- methylData[trainIndex,]
testData <- methylData[-trainIndex,] 
testDataY <- testData[,1] 
table(testDataY)

## Prespecificed parameters
supervisedStage1=TRUE
FSmethod1 <- NULL
cutP1 <- 0.05
FSmethod2 <- "positive"
cutP2 <- NULL
classifier <- "randForest"
predMode <- "probability"
paramlist <- list(ntree=100, nthreads=20)   
core <- MulticoreParam(workers = 10) 

set.seed(123)
result <- BioMM(trainData=trainData, testData=testData, 
                pathlistDB=golistSub, featureAnno, 
                restrictUp=200, restrictDown=10, minPathSize=10, 
                supervisedStage1, typePCA="regular", 
                resample1="BS", dataMode="allTrain",
                repeatA1=50, repeatA2=1, repeatB1=10, repeatB2=1, 
                nfolds=10, FSmethod1, FSmethod2, 
                cutP1, cutP2, fdr2=NULL, 
                FScore=MulticoreParam(workers = 1), 
                classifier, predMode, paramlist, innerCore=core)

metricTest <- getMetrics(dataY = testDataY, predY = result)
message("Test set prediction performance:")
print(metricTest)


``` 

## Stage-2 feature investigation 
### Generation of stage-2 feature
Stage-2 feature is reconstructed using resampling prediction on the training set or on the test set if the argument `testData` is provided. Here we use BioMM with a set of predetermined/optimal parameters to create stage-2 pathway-level data.  

```{r stage2dataAprep, eval=TRUE} 
pathType <- "GO"
# pathType <- "KEGG"
if (pathType == "GO"){
    studylist <- methylGOlist
} else if (pathType == "KEGG"){
    studylist <- methylKEGGlist
} else {
    stop("Wrong specified pathType!")
} 

## Model parameters 
supervisedStage1=TRUE
FSmethod <- NULL
cutP <- 0.05
classifier <- "randForest"
predMode <- "probability"
paramlist <- list(ntree=100, nthreads=20)   
core <- MulticoreParam(workers = 10) 

set.seed(123)
stage2dataA <- reconBySupervised(trainDataList=studylist, 
			                     testDataList=NULL,
			                     resample="BS", dataMode="allTrain",
			                     repeatA=50, repeatB=1, nfolds=10,
			                     FSmethod, cutP, fdr=NULL, 
			                     FScore=MulticoreParam(workers = 1),
			                     classifier, predMode, paramlist,
			                     innerCore=core, outFileA=NULL, outFileB=NULL)
## Check stage-2 data
dim(stage2dataA)
print(table(stage2dataA[,1]))
head(stage2dataA[,1:4])

``` 

### Feature Visualization
#### Explained variation of stage-2 data
The distribution of the proportion of variance explained for the individual generated feature of stage-2 data for the classification task is illustrated `plotVarExplained()` below. Nagelkerke pseudo R-squared measure is used to compute the explained variance. The argument `posF=TRUE` indicates that only positively outcome-associated features are plotted since negative associations likely reflect random effects in the underlying data.

``` {r stage2dataViz, eval=TRUE}
core <- MulticoreParam(workers = 10)   
fileName <- paste0(pathType, "_featuresVarExplained.png")
plotVarExplained(data=stage2dataA, posF=TRUE, binarize=FALSE, core=core, 
                 pathTitle=paste0(pathType, " pathways"), fileName)

plot(load.image(fileName)) 

``` 

#### Prioritization of outcome-associated functional patterns 
`plotRankedFeature()` is employed to rank and visualize the outcome-associated features from stage-2 data. The argument `topF=10` and `posF=TRUE` are used to define the top 10 positively outcome-associated features. The negative log P value using logistic regression is utilized to evaluate the importance of the ranked features as indicated by the argument `rankMetric="negPlogit"`. Other metrics including Nagelkerke pseudo R-squared "R2", and Z score "Zscore" measure are also available (see `plotRankedFeature` in the manual for details). The size of the respective pathway is pictured as the argument `colorMetric="size"`. 

``` {r topPathFeatures, eval=TRUE, fig.show="hold"} 
core <- MulticoreParam(workers = 1)   
rankMetric <- "negPlogit" 
filePrefix <- paste0(pathType, "_topPath_", rankMetric)
topPath <- plotRankedFeature(data=stage2dataA, 
                             posF=TRUE, topF=10, binarize=FALSE, 
                             blocklist=studylist,  
                             rankMetric=rankMetric, 
                             colorMetric="size",  core, 
                             pathTitle=paste0(pathType, " pathways"), 
                             fileName=paste0(filePrefix, ".png"))
plot(load.image(paste0(filePrefix, ".png")))   

``` 

The statistical metrics and descriptions of these above top pathways are shown below:

``` {r reportTopPath, eval=TRUE} 
## Report the top pathways
if (pathType == "GO"){
  goterms = unlist(Term(GOTERM))  
  topGOID = gsub("\\.", ":", rownames(topPath))
  subTerm = goterms[is.element(names(goterms), topGOID)] 
  topIDmatch = subTerm[match(topGOID, names(subTerm))]  
  topPath <- data.frame(topPath, Description=topIDmatch)
  
} else if (pathType == "KEGG"){
  ## A matching list between KEGG ID and names. Data freezes on Aug 2020
  keggID2name <- readRDS(system.file("extdata", "/keggID2name202008.rds", 
                                     package="BioMM"))  
  keggSub <- keggID2name[is.element(keggID2name[,"ID"], rownames(topPath)),]
  topIDmatch <- keggSub[match(rownames(topPath), keggSub[,"ID"]),] 
  topPath <- data.frame(topPath, Description=topIDmatch[,"name"])
}

print(topPath) 
write.xlsx(topPath,file=paste0(filePrefix, ".xlsx"))
# write.table(topPath,file=paste0(filePrefix, ".txt"), sep="\t")

``` 

#### The significance of CpGs in pathways of interest
`cirPlot4pathway()` illustrates the significance of the individual CpGs (for DNA methylation data) or genes (for gene expression data) falling into a set of pathways. Here the top 10 outcome-associated pathways are investigated. Negative log P value is used to define the significance of each CpG or genes within these pathways.

``` {r cirPlot, eval=TRUE, fig.show="hold"}  
core <- MulticoreParam(workers = 10)   
pathID <- gsub("\\.", ":", rownames(topPath))
## The number of top pathways must be bigger than overall investigated pathways
pathSet <- studylist[is.element(names(studylist), pathID)]
pathMatch <- pathSet[match(pathID, names(pathSet))]
fileName <- paste0(pathType, "_SigRankBy_", rankMetric)

cirPlot4pathway(datalist=pathMatch, topPathID=names(pathMatch), core, fileName)
``` 

## Computational consideration
BioMM with supervised models at both stages will take longer to run than unsupervised approaches. But the prediction can be more powerful. Given the potential of BioMM beyond its prediction performance, the computational load may be negligible. Furthermore, more efficient programming languages or environments can help reduce the runtime vastly, and the adoption of next-generation technology (e.g., 5G) is pushing advances in computational storage and speed. 

The stability of BioMM prediction is often facilitated with the increasing number of resampling repetitions and some other related model parameters such as the number of trees used in the Random Forest model. Due to the runtime, we only showcased the smaller examples and models with less computation in this vignette.


## Citation

NIPS ML4H submission: Chen, J. and Schwarz, E., 2017. BioMM: Biologically-informed Multi-stage Machine learning for identification of epigenetic fingerprints. arXiv preprint arXiv:1712.00336.
