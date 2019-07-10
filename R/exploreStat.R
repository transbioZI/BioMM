
# File : exploreStat.R Author : Junfang Chen


############################################################################### 

#' Return the data after feature selection 
#'
#' @description 
#' Identify and select a subset of outcome-associated or predictive features in 
#' the training data based on filtering methods. Return the same set of selected 
#' features for the test data if it is available. 

#' @param trainData The input training dataset. The first column is the label.  
#' @param testData The input test dataset. The first column is the label.   
#' @param FSmethod Feature selection methods. Available options are 
#' c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 'posWilcox', 
#' or 'top10pCor'). 'positive' is the positively outcome-associated features
#' using the Pearson correlation method. 'posWilcox' is the positively 
#' outcome-associated features using Pearson correlation method together 
#' with 'wilcox.text' method. 'top10pCor' is the top 10% of   
#' outcome-associcated features. This is helpful when no features can be 
#' picked during stringent feature selection procedure.
#' @param cutP The cutoff used for p value thresholding. It can be any value 
#' between 0 and 1. Commonly used cutoffs are c(0.5, 0.1, 0.05, 0.01, etc.). 
#' The default is 0.1. 
#' @param fdr Multiple testing correction method. Available options are 
#' c(NULL, 'fdr', 'BH', 'holm' etc). 
#' See also \code{\link[stats]{p.adjust}}. The default is NULL.
#' @param FScore The number of cores used for some feature selection methods.
#' The default is 10.

#' @return Both training and test data (if provided) with pre-selected 
#' features are returned if feature selection method is applied. 
#' If no feature can be selected during feature selection procedure, then the 
#' output is NULL. 
#' @details Parallel computing is helpful if your input data is high  
#' dimensional. For 'cutP', a soft thresholding of 0.1 may be favorable than 
#' more stringent p value cutoff because the features with small effect size 
#' can be taken into consideration for downstream analysis. However, for high 
#' dimensional (e.g. p > 10,000) data, many false positive features may 
#' exist, thus, rigorous p value thresholding should be applied.  
#' The choice of feature selection method depends on the characteristics of 
#' the input data.
#' @export 
#' @import BiocParallel
#' @author Junfang Chen
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' trainIndex <- sample(nrow(methylData), 20)
#' trainData = methylData[trainIndex,]
#' testData = methylData[-trainIndex,]
#' ## Feature selection
#' library(BiocParallel)
#' param <- MulticoreParam(workers = 10)
#' ## Select outcome-associated features based on the Wilcoxon test (P<0.1)
#' datalist <- getDataAfterFS(trainData, testData, FSmethod="wilcox.test", 
#'                            cutP=0.1, fdr=NULL, FScore=param)
#' trainDataSub <- datalist[[1]] 
#' testDataSub <- datalist[[2]] 
#' print(dim(trainData))
#' print(dim(trainDataSub))


getDataAfterFS <- function(trainData, testData, FSmethod, cutP = 0.1, 
    fdr = NULL, FScore = MulticoreParam()) {
    
    if (colnames(trainData)[1] != "label") {
        stop("The first column of the 'trainData' must be the 'label'!")
    }
    trainX <- trainData[, -1]
    trainY <- trainData[, 1]
    featureNames <- colnames(trainX)
    if (is.factor(trainY)) {
        trainY <- as.numeric(trainY) - 1
    }
    if (!is.null(testData)) {
        if (colnames(testData)[1] != "label") {
            stop("The first column of the 'testData' must be the 'label'!")
        }
        testX <- testData[, -1]
        testY <- testData[, 1]
    }
    
    if (is.null(FSmethod)) {
        ## no FS;
        selFeature <- seq_len(ncol(trainX))
    } else if (FSmethod == "positive") {
        ## only positively correlated use 'which' to keep the remaining index
        selFeature <- which(cor(trainX, trainY) > 0)
    } else if (FSmethod == "wilcox.test") {
        featurelist <- as.list(seq_len(ncol(trainX)))
        pvTrain <- unlist(bplapply(featurelist, function(i) {
            wilcox.test(trainX[, i] ~ as.factor(trainY))$p.value
        }, BPPARAM = FScore))
        selFeature <- which(pvTrain < cutP)
    } else if (FSmethod == "cor.test") {
        featurelist <- as.list(seq_len(ncol(trainX)))
        pvTrain <- unlist(bplapply(featurelist, function(i) {
            cor.test(trainX[, i], trainY)$p.value
        }, BPPARAM = FScore))
        selFeature <- which(pvTrain < cutP)
    } else if (FSmethod == "chisq.test") {
        featurelist <- as.list(seq_len(ncol(trainX)))
        pvTrain <- unlist(bplapply(featurelist, function(i) {
            if (length(table(trainX[, i])) != 1) {
                pv <- chisq.test(trainX[, i], as.factor(trainY))$p.value
                names(pv) <- featureNames[i]
                pv
            }
        }, BPPARAM = FScore))
        indexNew <- match(featureNames, names(pvTrain))
        pvTrain2 <- pvTrain[indexNew]
        selFeature <- which(pvTrain2 < cutP)
    } else if (FSmethod == "posWilcox") {
        whPos <- which(cor(trainX, trainY) > 0)
        featurelist <- as.list(seq_len(ncol(trainX)))
        pvTrain <- unlist(bplapply(featurelist, function(i) {
            wilcox.test(trainX[, i] ~ as.factor(trainY))$p.value
        }, BPPARAM = FScore))
        varTmp <- which(pvTrain < cutP)
        selFeature <- intersect(whPos, varTmp)
    } else if (FSmethod == "top10pCor") { 
        topN <- round(ncol(trainX)*0.1) 
        featurelist <- as.list(seq_len(ncol(trainX)))
        pvTrain <- unlist(mclapply(featurelist, function(i){
                          wilcox.test(trainX[,i]~as.factor(trainY))$p.value}, 
                          mc.cores=FScore))  
        selFeature <- order(pvTrain, decreasing=FALSE)[seq_len(topN)]   
        if (length(selFeature) <= 2){
            selFeature <- order(pvTrain, decreasing=FALSE)[seq_len(2)]
        }     
        print("top 10% sigWil Features")
    }
    
    if (!is.null(fdr)) {
        ## For FS based on fdr
        if (FSmethod != "positive" && !is.null(FSmethod)) {
            selFpAdj <- which(p.adjust(pvTrain, method = fdr) < cutP)
            selFeature <- intersect(selFeature, selFpAdj)
        }
    }
    
    if (length(selFeature) == 0) {
        ## return NULL
        message("Warning: No feature selected!")
        subTrain <- NULL
        subTest <- NULL
    } else if (length(selFeature) == 1) {
        ## in case only 1 feature is selected
        featName <- featureNames[selFeature]
        subTrain <- data.frame(label = trainY, one = trainX[, selFeature])
        colnames(subTrain) <- c("label", featName)
        if (!is.null(testData)) {
            subTest <- data.frame(label = testY, one = testX[, selFeature])
            colnames(subTest) <- c("label", featName)
        }
    } else if (length(selFeature) >= 2) {
        subTrain <- data.frame(label = trainY, trainX[, selFeature])
        if (!is.null(testData)) {
            subTest <- data.frame(label = testY, testX[, selFeature])
        }
    }
    
    if (!is.null(testData)) {
        result <- list(subTrain, subTest)
    } else {
        result <- subTrain
    }
    
}
