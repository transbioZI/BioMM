test_BaseModels <- function() {

    library(RUnit)
    library(BioMM) 
    ## Load data  
    methylfile <- system.file("extdata", "methylData.rds", package="BioMM")  
    methylData <- readRDS(methylfile)  
    dataY <- methylData[,1]
    ## test a subset of genome-wide methylation data at random
    methylSub <- data.frame(label=dataY, methylData[,c(2:1001)])  
    trainIndex <- sample(nrow(methylSub), 12)
    trainData = methylSub[trainIndex,]
    testData = methylSub[-trainIndex,]
    library(ranger)
    predMode <- "classification"
    paramlist <- list(ntree=300, nthreads=10)
    predY <- baseRandForest(trainData, testData, predMode, paramlist)  
    ## less than two levels for the binary classification predicted output
    checkTrue( nlevels(factor(predY)) <= 2 ) 

}
