
# File : BioMM.R Author: Junfang Chen


############################################################################### 

#' Prediction by random forest 
#'
#' @description
#' Prediction by random forest with different settings: 'probability', 
#' 'classification' and 'regression'.

#' @param trainData The input training dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param predMode The prediction mode. Available options are 
#' c('probability', 'classification', 'regression').
#' @param paramlist A set of model parameters defined in an R list object. 
#' The valid option: list(ntree, nthreads). 'ntree' is the number of trees 
#' used. The defaul is 2000. 'nthreads' is the number of threads used for 
#' computation. The default is 20.
#' @return The predicted output for the test data.
#' @import ranger
#' @export 
#' @author Junfang Chen
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' dataY <- methylData[,1]
#' ## test a subset of genome-wide methylation data at random
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' trainIndex <- sample(nrow(methylSub), 30)
#' trainData = methylSub[trainIndex,]
#' testData = methylSub[-trainIndex,]
#' library(ranger)
#' predY <- baseRandForest(trainData, testData, 
#'                         predMode='classification', 
#'                         paramlist=list(ntree=300, nthreads=20)) 
#' testY <- testData[,1]
#' accuracy <- classifiACC(dataY=testY, predY=predY)
#' print(accuracy) 


baseRandForest <- function(trainData, testData, predMode = c("classification", 
    "probability", "regression"), paramlist = list(ntree = 2000, nthreads = 20)) {
    
    predMode <- match.arg(predMode)
    ## input parameters
    if (!is.element("ntree", names(paramlist))) {
        stop(" 'ntree' is missing in the 'paramlist'!")
    }
    if (!is.element("nthreads", names(paramlist))) {
        stop(" 'nthreads' is missing in the 'paramlist'!")
    }
    ntree <- paramlist$ntree
    nthreads <- paramlist$nthreads
    if (predMode == "probability") {
        trainData$label <- as.factor(trainData$label)
        model <- ranger(label ~ ., data = trainData, num.trees = ntree,
            num.threads = nthreads, write.forest = TRUE, probability = TRUE)
        yhatClassProb <- predictions(predict(model, testData))
        predTest <- round(yhatClassProb[, 2], 3)
    } else if (predMode == "classification") {
        trainData$label <- as.factor(trainData$label)
        model <- ranger(label ~ ., data = trainData, num.trees = ntree,
            num.threads = nthreads, write.forest = TRUE, classification = TRUE)
        predTest <- predictions(predict(model, testData))
        predTest <- as.numeric(predTest) - 1
    } else if (predMode == "regression") {
        model <- ranger(label ~ ., data = trainData, num.trees = ntree, 
            num.threads = nthreads)
        predTest <- predictions(predict(model, testData))
        predTest <- round(predTest, 3)
    }

}


############################################################################### 

#' Prediction by SVM
#'
#' @description
#' Prediction by support vector machine (SVM) with two different settings: 
#' 'classification' and 'regression'.

#' @param trainData The input training dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param predMode The prediction mode. Available options are 
#' c('classification', 'probability', 'regression').
#' @param paramlist A set of model parameters defined in an R list object. 
#' The valid option: list(kernel, gamma, cost, tuneP).  
#' \enumerate{
#'   \item 'tuneP':  a logical value indicating if hyperparameter tuning 
#'                   should be conducted or not. The default is FALSE.
#'   \item 'kernel': options are c('linear', 'polynomial', 'radial', 
#'                   'sigmoid'). The defaut is 'radial'.
#'   \item 'gamma': the parameter needed for all kernels except 'linear'. 
#'         If tuneP is TRUE, more than one value is suggested.  
#'   \item 'cost': is the cost of constraints violation.
#'         If tuneP is TRUE, more than one value is suggested. 
#' } 
#' @return The predicted output for the test data.
#' @details Hyperparameter tuning is recommended in many biological data 
#' mining applications. The best parameters can be determined via an internal 
#' cross validation.
#' @import e1071
#' @export 
#' @author Junfang Chen  
##' @seealso \code{\link[e1071]{svm}}
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' dataY <- methylData[,1]
#' ## select a subset of genome-wide methylation data at random
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' trainIndex <- sample(nrow(methylSub), 30)
#' trainData = methylSub[trainIndex,]
#' testData = methylSub[-trainIndex,]
#' library(e1071)
#' predY <- baseSVM(trainData, testData, 
#'                  predMode='classification', 
#'                  paramlist=list(tuneP=FALSE, kernel='radial', 
#'                                 gamma=10^(-3:-1), cost=10^(-3:1))) 
#' testY <- testData[,1]
#' accuracy <- classifiACC(dataY=testY, predY=predY)
#' print(accuracy) 


baseSVM <- function(trainData, testData, 
    predMode = c("classification", "probability", "regression"), 
    paramlist = list(tuneP = TRUE, kernel = "radial", gamma = 10^(-3:-1), 
    cost = 10^(-2:2))) {
    
    predMode <- match.arg(predMode)
    ## input parameters
    if (!is.element("tuneP", names(paramlist))) {
        stop(" 'tuneP' is missing in the 'paramlist'!")
    }
    if (!is.element("kernel", names(paramlist))) {
        stop(" 'kernel' is missing in the 'paramlist'!")
    }
    if (!is.element("gamma", names(paramlist))) {
        stop(" 'gamma' is missing in the 'paramlist'!")
    }
    if (!is.element("cost", names(paramlist))) {
        stop(" 'cost' is missing in the 'paramlist'!")
    }
    tuneP <- paramlist$tuneP
    kernel <- paramlist$kernel
    gamma <- paramlist$gamma
    cost <- paramlist$cost
    if (tuneP == TRUE) {
        kernel <- paramlist$kernel
    }
    if (tuneP == TRUE) {
        gamma <- paramlist$gamma
    }
    if (tuneP == TRUE) {
        cost <- paramlist$cost
    }
    
    if (predMode == "classification" || predMode == "probability") {
        type <- "C-classification"  ## 
        trainData$label <- as.factor(trainData$label)
    } else if (predMode == "regression") {
        type <- "eps-regression"
        if (is.factor(trainData$label)) {
            trainData$label <- as.numeric(trainData$label) - 1
        }
    }
    
    if (predMode == "probability") {
        if (tuneP == TRUE) {
            ## use internal CV
            tuneOut <- tune.svm(label ~ ., data = trainData, type, kernel, 
                gamma, cost, probability = TRUE)
            model <- tuneOut$best.model
        } else {
            ## no internal CV required
            model <- svm(label ~ ., data = trainData, probability = TRUE)
        }
        predTest2 <- predict(model, testData, probability = TRUE)
        predTest <- round(attr(predTest2, "probabilities")[, 2], 3)
        
    } else {
        if (tuneP == TRUE) {
            ## use internal CV
            tuneOut <- tune.svm(label ~ ., data = trainData, type, kernel, 
                gamma, cost)
            model <- tuneOut$best.model
        } else {
            ## no CV required
            model <- svm(label ~ ., data = trainData)
        }
        predTest <- predict(model, testData)
        if (predMode == "classification") {
            predTest <- as.numeric(predTest) - 1
        }
    }

}


############################################################################### 

#' Prediction by generalized linear regression models
#'
#' @description
#' Prediction by generalized regression models with lasso or elastic net 
#' regularization.

#' @param trainData The input training dataset. The first column is named 
#' the 'label'.
#' @param testData The input test dataset. The first column is named 
#' the 'label'.
#' @param predMode The prediction mode. Available options are 
#' c('classification', 'probability', 'regression').
#' @param paramlist A set of model parameters defined in an R list object. 
#' The valid option: list(family, alpha, typeMeasure, typePred).  
#' \enumerate{
#'   \item 'family': Response type: 'gaussian','binomial','poisson',
#'                'multinomial','cox','mgaussian'. (Default: 'binomial')
#'   \item 'alpha': The elastic net mixing parameter, with 0<=alpha<= 1. 
#'   \item 'typeMeasure': error metrics for internal cross-validation.
#'                        'mse' uses squared loss; 
#'                        'deviance' uses actual deviance; 
#'                        'mae' uses mean absolute error; 
#'                        'class' gives misclassification error;
#'                        'auc' (for two-class logistic regression ONLY) 
#'                        gives area under the ROC curve.
#'   \item 'typePred': The type of prediction: 'response' and 'class'. 
#'                     (Default: 'class' for binary classification)
#' } 
#' @return The predicted output for the test data.
#' @export 
#' @import glmnet
#' @author Junfang Chen 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' dataY <- methylData[,1]
#' ## select a subset of genome-wide methylation data at random
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' trainIndex <- sample(nrow(methylSub), 30)
#' trainData = methylSub[trainIndex,]
#' testData = methylSub[-trainIndex,]
#' library(glmnet)
#' ## classification
#' predY <- baseGLMnet(trainData, testData, 
#'                     predMode='classification', 
#'                     paramlist=list(family='binomial', alpha=0.5, 
#'                                    typeMeasure='mse', typePred='class')) 
#' testY <- testData[,1]
#' accuracy <- classifiACC(dataY=testY, predY=predY)
#' print(accuracy)  


baseGLMnet <- function(trainData, testData, 
    predMode = c("classification", "probability", "regression"), 
    paramlist = list(family = "binomial", alpha = 0.5, typeMeasure = "mse", 
    typePred = "class")) {
    
    predMode <- match.arg(predMode)
    ## input parameters
    if (!is.element("family", names(paramlist))) {
        stop(" 'family' is missing in the 'paramlist'!")
    }
    if (!is.element("alpha", names(paramlist))) {
        stop(" 'alpha' is missing in the 'paramlist'!")
    }
    if (!is.element("typeMeasure", names(paramlist))) {
        stop(" 'typeMeasure' is missing in the 'paramlist'!")
    }
    if (!is.element("typePred", names(paramlist))) {
        stop(" 'typePred' is missing in the 'paramlist'!")
    }
    family <- paramlist$family
    alpha <- paramlist$alpha
    typeMeasure <- paramlist$typeMeasure
    typePred <- paramlist$typePred
    trainDataY <- trainData$label
    ## To avoid 'X should be a matrix with 2 or more cols'
    if (ncol(trainData) == 2) {
        trainData$one <- rep(1, nrow(trainData))
        testData$one <- rep(1, nrow(testData))
        trainDataX <- trainData[, -1]
        testDataX <- testData[, -1]
    } else {
        trainDataX <- trainData[, -1]
        testDataX <- testData[, -1]
    }
    cvfit <- cv.glmnet(as.matrix(trainDataX), trainDataY, family = family, 
        type.measure = typeMeasure, nfolds = 10)
    model <- glmnet(as.matrix(trainDataX), trainDataY, family = family, 
        alpha = alpha, lambda = cvfit$lambda.min)
    yhat <- predict(model, as.matrix(testDataX), type = typePred)
    
    if (predMode == "classification") {
        ## set type='class'
        predTest = as.numeric(yhat)
    } else if (predMode == "probability" || predMode == "regression") {
        predTest = round(yhat[, 1], 3)
    }

}


############################################################################### 

#' Base supervised machine learning models for prediction
#'
#' @description
#' Prediction using different supervised machine learning models.

#' @param trainData The input training dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param classifier Machine learning classifiers. 
#' Available options are c('randForest', 'SVM', 'glmnet').
#' @param predMode The prediction mode. Available options are
#' c('classification', 'probability', 'regression'). 
#' 'probability' is currently only for 'randForest'. 
#' @param paramlist A set of model parameters defined in an R list object. 
#' See more details for each individual model. 

#' @return The predicted output for the test data.
#' @export 

#' @author Junfang Chen 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' dataY <- methylData[,1]
#' ## select a subset of genome-wide methylation data at random
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' trainIndex <- sample(nrow(methylSub), 30)
#' trainData = methylSub[trainIndex,]
#' testData = methylSub[-trainIndex,]
#' library(ranger) 
#' set.seed(123)
#' predY <- baseModel(trainData, testData, 
#'                    classifier='randForest',  
#'                    predMode='classification', 
#'                    paramlist=list(ntree=300, nthreads=20)) 
#' print(table(predY)) 
#' testY <- testData[,1]
#' accuracy <- classifiACC(dataY=testY, predY=predY)
#' print(accuracy)  

baseModel <- function(trainData, testData, 
    classifier = c("randForest", "SVM", "glmnet"), 
    predMode = c("classification", "probability", "regression"), paramlist) {
    
    if (colnames(trainData)[1] != "label") {
        stop("The first column of the 'trainData' must be the 'label'!")
    }
    if (colnames(testData)[1] != "label") {
        stop("The first column of the 'testData' must be the 'label'!")
    }
    classifier <- match.arg(classifier)
    predMode <- match.arg(predMode)
    if (classifier == "randForest") {
        predTest <- baseRandForest(trainData, testData, predMode, paramlist)
    } else if (classifier == "SVM") {
        predTest <- baseSVM(trainData, testData, predMode, paramlist)
    } else if (classifier == "glmnet") {
        predTest <- baseGLMnet(trainData, testData, predMode, paramlist)
    }

}


############################################################################### 

#' Prediction by supervised machine learning along with feature selection
#'
#' @description
#' Prediction by supervised machine learning along with feature selection.

#' @param trainData The input training dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param FSmethod Feature selection methods. Available options are 
#' c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 'posWilcox', 
#' or 'top10pCor').
#' @param cutP The cutoff used for p value thresholding.  
#' Commonly used cutoffs are c(0.5, 0.1, 0.05, etc.). The default is 0.05.
#' @param fdr Multiple testing correction method. Available options are 
#' c(NULL, 'fdr', 'BH', 'holm', etc.). 
#' See also \code{\link[stats]{p.adjust}}. The default is NULL.
#' @param FScore The number of cores used for feature selection.  
#' @param classifier Machine learning classifiers. 
#' Available options are c('randForest', 'SVM', 'glmnet').
#' @param predMode The prediction mode. Available options are 
#' c('probability', 'classification', 'regression'). 
#' @param paramlist A set of model parameters defined in an R list object. 

#' @return The predicted output for the test data.
#' @details If no feature selected or just one selected feature, 
#' then top 10% outcome-associated features will be selected by default.
#' @export 

#' @author Junfang Chen 
#' @seealso \code{\link{getDataAfterFS}} 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' dataY <- methylData[,1]
#' ## select a subset of genome-wide methylation data at random
#' methylSub <- data.frame(label=dataY, methylData[,c(2:501)])  
#' trainIndex <- sample(nrow(methylSub), 30)
#' trainData = methylSub[trainIndex,]
#' testData = methylSub[-trainIndex,]
#' library(ranger) 
#' library(BiocParallel)
#' param <- MulticoreParam(workers = 10)
#' predY <- predByFS(trainData, testData, 
#'                   FSmethod='cor.test', cutP=0.1, 
#'                   fdr=NULL, FScore=param, 
#'                   classifier='randForest',
#'                   predMode='classification', 
#'                   paramlist=list(ntree=300, nthreads=20))  
#' testY <- testData[,1]
#' accuracy <- classifiACC(dataY=testY, predY=predY)
#' print(accuracy)  


predByFS <- function(trainData, testData, FSmethod, cutP, fdr, 
    FScore = MulticoreParam(), classifier, predMode, paramlist) {
    
    datalist <- getDataAfterFS(trainData, testData, FSmethod, cutP, fdr, FScore)
    ## include the label
    trainDataSub <- datalist[[1]]
    testDataSub <- datalist[[2]]
    
    ## If no selected features or just one selected feature.
    if (is.null(trainDataSub) || ncol(trainDataSub) == 2) {
        datalist <- getDataAfterFS(trainData, testData, FSmethod = "top10pCor", 
            cutP, fdr, FScore)
        trainDataSub <- datalist[[1]]
        testDataSub <- datalist[[2]]
        predTest <- baseModel(trainData = trainDataSub, testData = testDataSub, 
            classifier, predMode, paramlist)
        
    } else if (!is.null(trainDataSub)) {
        predTest <- baseModel(trainData = trainDataSub, testData = testDataSub, 
            classifier, predMode, paramlist)
        
    } else {
        message("Error: Input data is wrong!")
    }
    
    if (is.factor(predTest)) {
        predTest <- as.numeric(predTest) - 1
    }

    return(predTest)
}



############################################################################### 

#' Bootstrap resampling prediction via supervised machine learning with 
#' feature selection
#'
#' @description
#' Prediction via supervised machine learning using bootstrap resampling 
#' along with feature selection methods.

#' @param trainData The input training dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param dataMode The input training data mode for model training.
#' It is used only if 'testData' is present. It can be a subset of 
#' the whole training data or the entire training data. 'subTrain' 
#' is the given for subsetting and 'allTrain' for the entire training
#' dataset. 
#' @param repeats The number of repeats used for boostrapping.
#' @param FSmethod Feature selection methods. Available options are 
#' c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 'posWilcox', 
#' or 'top10pCor').
#' @param cutP The cutoff used for p value thresholding.  
#' Commonly used cutoffs are c(0.5, 0.1, 0.05, 0.01, etc). The default is 0.05.
#' @param fdr Multiple testing correction method. Available options are 
#' c(NULL, 'fdr', 'BH', 'holm', etc). 
#' See also \code{\link[stats]{p.adjust}}. The default is NULL.
#' @param FScore The number of cores used for feature selection 
#' if parallel computing needed.
#' @param classifier Machine learning classifiers. 
#' @param predMode The prediction mode. Available options are 
#' c('probability', 'classification', 'regression').
#' @param paramlist A set of model parameters defined in an R list object. 
#' @param innerCore The number of cores used for computation.

#' @return The predicted output for the test data. 
#' @export  
#' @import BiocParallel 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' dataY <- methylData[,1]
#' ## select a subset of genome-wide methylation data at random
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' trainIndex <- sample(nrow(methylSub), 30)
#' trainData = methylSub[trainIndex,]
#' testData = methylSub[-trainIndex,]
#' library(ranger) 
#' library(BiocParallel)
#' param1 <- MulticoreParam(workers = 1)
#' param2 <- MulticoreParam(workers = 20)
#' predY <- predByBS(trainData, testData, 
#'                   dataMode='allTrain', repeats=50,
#'                   FSmethod=NULL, cutP=0.1, 
#'                   fdr=NULL, FScore=param1, 
#'                   classifier='randForest',
#'                   predMode='classification', 
#'                   paramlist=list(ntree=300, nthreads=10),
#'                   innerCore=param2)  
#' testY <- testData[,1]
#' accuracy <- classifiACC(dataY=testY, predY=predY)
#' print(accuracy)  


predByBS <- function(trainData, testData, dataMode, repeats, FSmethod, cutP, 
    fdr, FScore = MulticoreParam(), classifier, predMode, paramlist, 
    innerCore = MulticoreParam()) {
    
    if (is.null(testData)) {
        ## BS only applied for training data
        data <- trainData
        predTmp <- rep(NA, nrow(data))
        replist <- seq_len(repeats)
        predTmpList <- bplapply(replist, function(reps) {
            trainIndex <- unique(sample(nrow(data), replace = TRUE))
            testIndex <- setdiff(seq_len(nrow(data)), trainIndex)
            trainData <- data[trainIndex, ]
            testData <- data[testIndex, ]
            ## OOB predicted scores
            predTmp[testIndex] <- predByFS(trainData, testData, FSmethod, cutP, 
                fdr, FScore, classifier, predMode, paramlist)
            predTmp
        }, BPPARAM = innerCore)
        predTest <- round(rowMeans(do.call(cbind, predTmpList), na.rm = TRUE), 
            3)
        naCount <- sum(is.na(predTest))
        ## If not enough bootstrapping repeats
        if (naCount > 0) {
            message(paste0("Warning: ", naCount, " NA(s) in predicted testY!"))
        }
        
    } else if (!is.null(testData)) {
        ## repeated prediction on test data
        if (dataMode == "subTrain") {
            data <- trainData
            predTmp <- rep(NA, nrow(data))
            replist <- seq_len(repeats)
            predTmpList <- bplapply(replist, function(reps) {
                trainIndex <- unique(sample(nrow(data), replace = TRUE))
                trainData <- data[trainIndex, ]
                predTmp <- predByFS(trainData, testData, FSmethod, cutP, fdr, 
                  FScore, classifier, predMode, paramlist)
                predTmp
            }, BPPARAM = innerCore)
            predTest <- round(rowMeans(do.call(cbind, predTmpList)), 3)
        } else if (dataMode == "allTrain") {
            predTmp <- rep(NA, nrow(testData))
            replist <- seq_len(repeats)
            predTmpList <- bplapply(replist, function(reps) {
                predTmp <- predByFS(trainData, testData, FSmethod, cutP, fdr, 
                  FScore, classifier, predMode, paramlist)
                predTmp
            }, BPPARAM = innerCore)
            predTest <- round(rowMeans(do.call(cbind, predTmpList)), 3)
        }
    } else {
        message("Input data is wrong!")
    }
    
    if (predMode == "classification") {
        predTest <- ifelse(predTest >= 0.5, 1, 0)
    }
    return(predTest)
}


############################################################################### 

#' Cross validation prediction by supervised machine learning and feature 
#' selection

#' @description
#' Prediction by supervised machine learning models using cross validation 
#' along with feature selection methods.

#' @param data The input dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member. 
#' @param repeats The number of repeats used for cross validation.
#' Repeated cross validation is performed if N >=2.
#' @param nfolds The number of folds is defined for cross validation.
#' @param FSmethod Feature selection methods. Available options are 
#' c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 'posWilcox', 
#' or 'top10pCor').
#' @param cutP The cutoff used for p value thresholding.  
#' Commonly used cutoffs are c(0.5, 0.1, 0.05, 0.01, etc). The default is 0.05.
#' @param fdr Multiple testing correction method. Available options are 
#' c(NULL, 'fdr', 'BH', 'holm', etc). 
#' See also \code{\link[stats]{p.adjust}}. The default is NULL.
#' @param FScore The number of cores used for feature selection if parallel 
#' computing needed. 
#' @param classifier Machine learning classifiers. 
#' @param predMode The prediction mode. Available options are 
#' c('probability', 'classification', 'regression').
#' @param paramlist A set of model parameters defined in an R list object. 
#' @param innerCore The number of cores used for computation.

#' @return The predicted cross validation output. 
#' @export 
#' @import BiocParallel 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)   
#' dataY <- methylData[,1]
#' ## select a subset of genome-wide methylation data at random
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' library(ranger) 
#' library(BiocParallel)
#' param1 <- MulticoreParam(workers = 1)
#' param2 <- MulticoreParam(workers = 20)
#' predY <- predByCV(methylSub, repeats=1, nfolds=10,   
#'                   FSmethod=NULL, cutP=0.1, 
#'                   fdr=NULL, FScore=param1, 
#'                   classifier='randForest',
#'                   predMode='classification', 
#'                   paramlist=list(ntree=300, nthreads=1),
#'                   innerCore=param2)
#' dataY <- methylData[,1]
#' accuracy <- classifiACC(dataY=dataY, predY=predY)
#' print(accuracy)  


predByCV <- function(data, repeats, nfolds, FSmethod, cutP, fdr, FScore = MulticoreParam(), 
    classifier, predMode, paramlist, innerCore = MulticoreParam()) {
    
    cvMat <- c()
    replist <- seq_len(repeats)
    for (reps in replist) {
        foldlists <- split(sample(nrow(data)), rep(seq_len(nfolds), 
            length = nrow(data)))
        cvY <- rep(NA, nrow(data))
        cvEstimate <- bplapply(foldlists, function(fold) {
            trainData <- data[-fold, ]
            testData <- data[fold, ]
            predTest <- predByFS(trainData, testData, FSmethod, cutP, fdr, 
                FScore, classifier, predMode, paramlist)
        }, BPPARAM = innerCore)
        cvY[unlist(foldlists)] <- unlist(cvEstimate)
        cvMat <- cbind(cvMat, cvY)
    }
    
    if (repeats == 1) {
        predTest <- cvMat[, 1]
    } else {
        ## average over the repeats
        predTest <- round(rowMeans(cvMat, na.rm = TRUE), 3)
    }
    
    if (predMode == "classification") {
        predTest <- ifelse(predTest >= 0.5, 1, 0)
    }

}


############################################################################### 

#' Reconstruct stage-2 data by supervised machine learning prediction

#' @description
#' Reconstruct stage-2 data by supervised machine learning prediction.

#' @param trainDataList The input training data list containing ordered 
#' collections of matrices.
#' @param testDataList The input test data list containing ordered collections 
#' of matrices. 
#' @param resample The resampling methods. Valid options are 'CV' and 'BS'. 
#' 'CV' for cross validation and 'BS' for bootstrapping resampling.
#' The default is 'BS'.
#' @param dataMode The mode of data used. 'subTrain' or 'allTrain'. 
#' @param repeatA The number of repeats N is used during resampling procedure.
#' Repeated cross validation or multiple boostrapping is performed if N >=2. 
#' One can choose 10 repeats for 'CV' and 100 repeats for 'BS'.
#' @param repeatB The number of repeats N is used for generating test data 
#' prediction scores. 
#' @param nfolds The number of folds is defined for cross validation.
#' @param FSmethod Feature selection methods. Available options are 
#' c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 'posWilcox', 
#' or 'top10pCor').
#' @param cutP The cutoff used for p value thresholding.  
#' Commonly used cutoffs are c(0.5, 0.1, 0.05, 0.01, etc). The default is 0.05.
#' @param fdr Multiple testing correction method. Available options are 
#' c(NULL, 'fdr', 'BH', 'holm', etc). 
#' See also \code{\link[stats]{p.adjust}}. The default is NULL.
#' @param FScore The number of cores used for feature selection, if parallel 
#' computing needed. 
#' @param classifier Machine learning classifiers. 
#' @param predMode The prediction mode. Available options are 
#' c('probability', 'classification', 'regression').
#' @param paramlist A set of model parameters defined in an R list object. 
#' @param innerCore The number of cores used for computation.
#' @param outFileA The file name of stage-2 training data with the '.rds' 
#' file extension. If it's provided, then the result will be saved in this 
#' file. The default is NULL.
#' @param outFileB The file name of stage-2 training data with the '.rds' 
#' file extension. If it's provided, then the result will be saved in this 
#' file. The default is NULL.

#' @return The predicted stage-2 training data and also stage-2 test data, if 
#' 'testDataList' provided. If outFileA and outFileB are provided, 
#' then the results will be stored in the files.
#' @details Stage-2 training data can be learned either using bootstrapping 
#' or cross validation resampling methods. Stage-2 test data is learned via 
#' independent test set prediction.
#' @export  
#' @import stats
#' @import BiocParallel 
#' @author Junfang Chen 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)  
#' ## Annotation files for Mapping CpGs into chromosome  
#' probeAnnoFile <- system.file('extdata', 'cpgAnno.rds', package='BioMM')  
#' probeAnno <- readRDS(file=probeAnnoFile)  
#' ## Mapping CpGs into Chromosome
#' dataList <- omics2chrlist(data=methylData, probeAnno)
#' length(dataList)
#' library(ranger) 
#' library(BiocParallel)
#' param1 <- MulticoreParam(workers = 1)
#' param2 <- MulticoreParam(workers = 20)
#' ## Not Run
#' ## stage2data <- BioMMreconData(trainDataList=dataList, testDataList=NULL, 
#' ##                             resample='CV', dataMode='allTrain', 
#' ##                             repeatA=1, repeatB=1, nfolds=10, 
#' ##                             FSmethod=NULL, cutP=0.1, 
#' ##                             fdr=NULL, FScore=param1, 
#' ##                             classifier='randForest',
#' ##                             predMode='classification', 
#' ##                             paramlist=list(ntree=300, nthreads=20),
#' ##                             innerCore=param2, outFileA=NULL, outFileB=NULL) 
#' ## print(dim(stage2data))
#' ## print(head(stage2data[,1:5]))


BioMMreconData <- function(trainDataList, testDataList, resample = "BS", 
    dataMode, repeatA, repeatB, nfolds, FSmethod, cutP, fdr, FScore = MulticoreParam(), classifier, 
    predMode, paramlist, innerCore = MulticoreParam(), outFileA = NULL, outFileB = NULL) {
    
    reconDataAx <- c()
    reconDataBx <- c()
    for (i in seq_along(trainDataList)) {
        ## for each block
        trainData = trainDataList[[i]]
        dataAy <- trainData[, 1]
        if (resample == "CV") {
            predA <- predByCV(data = trainData, repeats = repeatA, nfolds, 
                FSmethod, cutP, fdr, FScore, classifier, predMode, paramlist, 
                innerCore)
        } else if (resample == "BS") {
            predA <- predByBS(trainData, testData = NULL, dataMode, 
                repeats = repeatA, FSmethod, cutP, fdr, FScore, classifier, 
                predMode, paramlist, innerCore)
        }
        reconDataAx <- cbind(reconDataAx, predA)
        
        if (!is.null(testDataList)) {
            ## ind. prediction
            testData = testDataList[[i]]
            dataBy <- testData[, 1]
            predB <- predByBS(trainData, testData, dataMode, repeats = repeatB, 
                FSmethod, cutP, fdr, FScore, classifier, predMode, paramlist, 
                innerCore)
            reconDataBx <- cbind(reconDataBx, predB)
        }
    }
    
    row.names(reconDataAx) <- row.names(trainData)  ## 
    reconDataA <- data.frame(label = dataAy, reconDataAx)
    colnames(reconDataA) <- c("label", names(trainDataList))
    ## Provide fileName
    if (!is.null(outFileA)) {
        saveRDS(reconDataA, file = outFileA)
    }
    if (!is.null(testDataList)) {
        ## ind. prediction
        row.names(reconDataBx) <- row.names(testData)  ## 
        reconDataB <- data.frame(label = dataBy, reconDataBx)
        colnames(reconDataB) <- c("label", names(testDataList))
        if (!is.null(outFileB)) {
            saveRDS(reconDataB, file = outFileB)
        }
        result <- list(reconDataA, reconDataB)
    } else {
        result <- reconDataA
    }

}



############################################################################### 


#' Prediction performance for stage-2 data using supervised machine learning

#' @description
#' Prediction performance for reconstructed stage-2 data using supervised 
#' machine learning with feature selection methods.

#' @param trainData The input training dataset (stage-2 data). The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset (stage-2 data). The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param resample The resampling methods. Valid options are 'CV' and 'BS'. 
#' 'CV' for cross validation and 'BS' for bootstrapping resampling.
#' The default is 'CV'.
#' @param dataMode The mode of data used. 'subTrain' or 'allTrain'.
#' @param repeatA The number of repeats N is used during resampling prediction. 
#' The default is 1.
#' @param repeatB The number of repeats N is used for test data prediction. 
#' The default is 1. 
#' @param nfolds The number of folds is defined for cross validation.
#' @param FSmethod Feature selection methods. Available options are 
#' c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 'posWilcox', 
#' or 'top10pCor').
#' @param cutP The cutoff used for p value thresholding.  
#' Commonly used cutoffs are c(0.5, 0.1, 0.05, 0.01, etc). The default is 0.05.
#' @param fdr Multiple testing correction method. Available options are 
#' c(NULL, 'fdr', 'BH', 'holm', etc). 
#' See also \code{\link[stats]{p.adjust}}. The default is NULL.
#' @param FScore The number of cores used for feature selection if parallel 
#' computing needed. 
#' @param classifier Machine learning classifiers. 
#' @param predMode The prediction mode. Available options are 
#' c('probability', 'classification', 'regression').
#' @param paramlist A set of model parameters defined in an R list object. 
#' @param innerCore The number of cores used for computation. 

#' @return The CV or BS prediction performance for stage-2 training data and 
#' test set prediction performance for stage-2 test data if the test set is 
#' given.

#' @details Stage-2 prediction is performed typically using positively  
#' correlated features. Since negative associations likely reflect random 
#' effects in the underlying data
#' @export  
#' @import stats
#' @import utils  
#' @references Perlich, C., & Swirszcz, G. (2011). On cross-validation and  
#' stacking: Building seemingly predictive models on random data. ACM SIGKDD    
#' Explorations Newsletter, 12(2), 11-15. 
#' @author Junfang Chen 


BioMMstage2pred <- function(trainData, testData, resample = "CV", dataMode, 
    repeatA = 1, repeatB = 1, nfolds, FSmethod, cutP, fdr, FScore = MulticoreParam(), classifier, 
    predMode, paramlist, innerCore = MulticoreParam()) {
    
    resample <- match.arg(resample)
    if (!is.null(resample)) {
        if (is.factor(trainData$label)) {
            trainData$label <- as.numeric(trainData$label) - 1
        }
        trainDataY <- trainData$label
        if (resample == "CV") {
            predY <- predByCV(data = trainData, repeats = repeatA, nfolds, 
                FSmethod, cutP, fdr, FScore, classifier, predMode, paramlist, 
                innerCore)
            message("CrossValidation >>> ")
        } else if (resample == "BS") {
            predY <- predByBS(trainData, testData = NULL, dataMode, 
                repeats = repeatA, FSmethod, cutP, fdr, FScore, classifier, 
                predMode, paramlist, innerCore)
            message("Bootstrapping >>> ")
        }
        
        if (predMode == "probability") {
            predY <- ifelse(predY >= 0.5, 1, 0)
            metricCV <- getMetrics(dataY = trainDataY, predY)
        } else if (predMode == "classification") {
            metricCV <- getMetrics(dataY = trainDataY, predY)
        } else if (predMode == "regression") {
            metricCV <- cor(trainDataY, predY) 
        }
    }
    
    if (!is.null(testData)) {
        ## ind. prediction
        if (is.factor(testData$label)) {
            testData$label <- as.numeric(testData$label) - 1
        }
        testY <- testData$label
        predTest <- predByBS(trainData, testData, dataMode, repeats = repeatB, 
            FSmethod, cutP, fdr, FScore, classifier, predMode, paramlist, 
            innerCore)
        ## Prediction performance for the ind. test performance
        message(paste0("Test set performance: "))
        if (predMode == "probability") {
            predTest <- ifelse(predTest >= 0.5, 1, 0)
            metricTest <- getMetrics(dataY = testY, predTest)
        } else if (predMode == "classification") {
            metricTest <- getMetrics(dataY = testY, predTest)
        } else if (predMode == "regression") {
            metricTest <- cor(testY, predTest)
        }
        result <- list(metricCV, metricTest)
    } else {
        result <- metricCV
    }

}




############################################################################### 

#' Reconstruct stage-2 data by PCA 

#' @description
#' Stage-2 data reconstruction by regular or sparse constrained principal 
#' component analysis (PCA).  

#' @param trainDataList The input training data list containing ordered 
#' collections of matrices.
#' @param testDataList The input test data list containing ordered collections 
#' of matrices.
#' @param typeMode The type of PCA prediction mode. Available options are 
#' c('regular', 'sparse'). 
#' (Default: regular)
#' @param topPC The number of top PCs selected. The default is 1, 
#' i.e. the first PC. 
#' @param innerCore The number of cores used for computation. 
#' @param outFileA The file name of stage-2 training data with the '.rds' file 
#' extension.
#' If it's provided, then the result will be saved in this file. 
#' The default is NULL.
#' @param outFileB The file name of stage-2 training data with the '.rds' file 
#' extension. If it's provided, then the result will be saved in this file. 
#' The default is NULL.

#' @return The predicted stage-2 training data and also stage-2 test data if 
#' 'testDataList' provided. If outFileA and outFileB are provided then the 
#' results will be stored in the files.

#' @export 
#' @import nsprcomp
#' @import stats
#' @import BiocParallel

#' @author Junfang Chen   
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)    
#' ## Annotation files for Mapping CpGs into chromosome  
#' probeAnnoFile <- system.file('extdata', 'cpgAnno.rds', package='BioMM')  
#' probeAnno <- readRDS(file=probeAnnoFile)  
#' ## Mapping CpGs into Chromosome
#' dataList <- omics2chrlist(data=methylData, probeAnno)
#' length(dataList) 
#' library(BiocParallel)
#' param <- MulticoreParam(workers = 10) 
#' stage2data <- BioMMstage1pca(trainDataList=dataList, testDataList=NULL, 
#'                              typeMode='regular', topPC=1,  
#'                              innerCore=param, outFileA=NULL, outFileB=NULL) 
#' print(dim(stage2data))
#' print(head(stage2data[,1:5]))


BioMMstage1pca <- function(trainDataList, testDataList, typeMode = "regular", 
    topPC = 1, innerCore = MulticoreParam(), outFileA = NULL, outFileB = NULL) {
    
    reconDataAx <- c()  ## for stage-2 training data   
    trainData <- trainDataList[[1]]
    trainDataY <- trainData[, 1]
    if (!is.null(testDataList)) {
        reconDataBx <- c()
        testData <- testDataList[[1]]
        testDataY <- testData[, 1]
    }
    
    numlist <- seq_along(trainDataList)
    predMat <- bplapply(numlist, function(i) {
        trainData <- trainDataList[[i]]
        trainDataX = trainData[, -1]
        if (!is.null(testDataList)) {
            testData <- testDataList[[i]]
            testDataX = testData[, -1]
        }
        ## remove zero variance columns from trainData
        whNon0var <- apply(trainDataX, 2, var) != 0
        if (sum(whNon0var) == 0) {
            ## if all features constant..
            predA <- trainDataX[, 1]
            if (!is.null(testDataList)) {
                predB <- testDataX[, 1]
            }
        } else {
            trainDataX2 <- trainDataX[, whNon0var]
            if (!is.null(testDataList)) {
                testX2 <- as.matrix(testDataX[, whNon0var])
            }
            if (typeMode == "regular") {
                pca = prcomp(trainDataX2, center = TRUE, scale = TRUE)
                predA <- pca$x[, seq_len(min(ncol(pca$x), topPC))]
                if (!is.null(testDataList)) {
                  predB <- predict(pca, testX2)[, seq_len(min(ncol(pca$x), topPC))]
                }
            } else if (typeMode == "sparse") {
                nspc <- nsprcomp(trainDataX2, ncomp = topPC)
                predA <- nspc$x
                if (!is.null(testDataList)) {
                  predB <- predict(nspc, testX2)[, topPC]
                }
            }
        }
        predA <- round(predA, 3)
        pred <- predA
        if (!is.null(testDataList)) {
            predB <- round(predB, 3)
            pred <- list(predA, predB)
        }
        pred
    }, BPPARAM = innerCore)
    
    if (!is.null(testDataList)) {
        ## PC1
        predMatA <- lapply(predMat, function(data) {
            data[[1]]
        })
        reconDataAx <- matrix(unlist(predMatA), ncol = length(predMatA))
        row.names(reconDataAx) <- row.names(trainData)  ## 
        colnames(reconDataAx) <- names(trainDataList)
        reconDataA = data.frame(label = trainDataY, reconDataAx)
        if (!is.null(outFileA)) {
            saveRDS(reconDataA, file = outFileA)
        }
        ## dataB
        predMatB <- lapply(predMat, function(data) {
            data[[2]]
        })
        reconDataBx <- matrix(unlist(predMatB), ncol = length(predMatB))
        row.names(reconDataBx) <- row.names(testData)  ## 
        colnames(reconDataBx) <- names(trainDataList)
        reconDataB = data.frame(label = testDataY, reconDataBx)
        if (!is.null(outFileB)) {
            saveRDS(reconDataB, file = outFileB)
        }
        result <- list(reconDataA, reconDataB)
    } else {
        ## Note 'predMat' is different from above
        reconDataAx <- matrix(unlist(predMat), ncol = length(predMat))
        row.names(reconDataAx) <- row.names(trainData)  ## 
        colnames(reconDataAx) <- names(trainDataList)
        reconDataA = data.frame(label = trainDataY, reconDataAx)
        if (!is.null(outFileA)) {
            saveRDS(reconDataA, file = outFileA)
        }
        result <- reconDataA
    }

}



############################################################################### 

#' BioMM end-to-end prediction 

#' @description
#' End-to-end prediction by BioMM framework using either supervised or 
#' unsupervised learning at stage-1, then supervised learning at stage-2. 

#' @param trainData The input training dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset. The first column
#' is the label or the output. For binary classes, 
#' 0 and 1 are used to indicate the class member. 
#' @param pathlistDB A list of pathways with pathway IDs and their 
#' corresponding genes ('entrezID' is used). This is only used for 
#' pathway-based stratification (only \code{stratify} is 'pathway'). 
#' @param featureAnno The annotation data stored in a data.frame for probe 
#' mapping. It must have at least two columns named 'ID' and 'entrezID'.  
#' If it's NULL, then the input probe is from the transcriptomic data. 
#' (Default: NULL)
#' @param restrictUp The upper-bound of the number of probes or genes in each  
#' biological stratified block. 
#' @param restrictDown The lower-bound of the number of probes or genes in 
#' each biological stratified block.
#' @param minPathSize The minimal defined pathway size after mapping your 
#' own data to GO database. This is only used for 
#' pathway-based stratification (only \code{stratify} is 'pathway'). 
#' @param supervisedStage1 A logical value. If TRUE, then supervised learning 
#' models are applied; if FALSE, unsupervised learning.
#' @param typePCA the type of PCA. Available options are c('regular', 
#' 'sparse').
#' @param resample1 The resampling methods at stage-1. Valid options are 
#' 'CV' and 'BS'. 'CV' for cross validation and 'BS' for bootstrapping 
#' resampling. The default is 'BS'.
#' @param resample2 The resampling methods at stage-2. Valid options are 'CV' 
#' and 'BS'. 'CV' for cross validation and 'BS' for bootstrapping resampling.
#' The default is 'CV'. 
#' @param repeatA1 The number of repeats N is used during resampling procedure.
#' Repeated cross validation or multiple boostrapping is performed if N >=2. 
#' One can choose 10 repeats for 'CV' and 100 repeats for 'BS'.
#' @param repeatA2 The number of repeats N is used during resampling 
#' prediction. The default is 1 for 'CV'.  
#' @param repeatB1 The number of repeats N is used for generating stage-2 test 
#' data prediction scores. The default is 20.
#' @param repeatB2 The number of repeats N is used for test data prediction. 
#' The default is 1. 
#' @param nfolds The number of folds is defined for cross validation. 
#' The default is 10.
#' @param FSmethod1 Feature selection methods at stage-1. Available options 
#' are c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 
#' 'posWilcox'). 
#' @param FSmethod2 Feature selection methods at stage-2. Available options  
#' are c(NULL, 'positive', 'wilcox.test', 'cor.test', 'chisq.test', 
#' 'posWilcox').
#' @param cutP1 The cutoff used for p value thresholding at stage-1.  
#' Commonly used cutoffs are c(0.5, 0.1, 0.05, 0.01, etc). 
#' @param cutP2 The cutoff used for p value thresholding at stage-2.   
#' @param fdr2 Multiple testing correction method at stage-2. 
#' Available options are c(NULL, 'fdr', 'BH', 'holm', etc). 
#' See also \code{\link[stats]{p.adjust}}. The default is NULL.
#' @param FScore The number of cores used for feature selection.
#' @param classifier Machine learning classifiers at both stages.  
#' @param predMode The prediction mode at both stages. Available options are 
#' c('probability', 'classification', 'regression'). 
#' @param paramlist A list of model parameters at both stages.  
#' @param innerCore The number of cores used for computation.

#' @return The CV or BS prediction performance for the training data and 
#' test set prediction performance if \code{testData} is given.
#' @details Stage-2 training data can be learned either using bootstrapping 
#' or cross validation resampling methods in the supervised learning settting.
#' Stage-2 test data is learned via independent test set prediction.
#' @export  

#' @references Chen, J., & Schwarz, E. (2017). BioMM: Biologically-informed 
#' Multi-stage Machine learning for identification of epigenetic fingerprints. 
#' arXiv preprint arXiv:1712.00336.
#' @references Perlich, C., & Swirszcz, G. (2011). On cross-validation and  
#' stacking: Building seemingly predictive models on random data. ACM SIGKDD    
#' Explorations Newsletter, 12(2), 11-15. 

#' @seealso \code{\link{BioMMreconData}}; \code{\link{BioMMstage1pca}}; 
#' \code{\link{BioMMstage2pred}}

#' @examples  
#' ## Load data    
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)    
#' ## Annotation files for Mapping CpGs into chromosome  
#' probeAnnoFile <- system.file('extdata', 'cpgAnno.rds', package='BioMM')  
#' probeAnno <- readRDS(file=probeAnnoFile)   
#' supervisedStage1=TRUE
#' classifier <- 'randForest'
#' predMode <- 'classification'
#' paramlist <- list(ntree=300, nthreads=30)   
#' library(BiocParallel)
#' library(ranger)
#' param1 <- MulticoreParam(workers = 2)
#' param2 <- MulticoreParam(workers = 20)
#' ## Not Run 
#' ## result <- BioMM(trainData=methylData, testData=NULL,
#' ##                 pathlistDB, featureAnno=probeAnno, 
#' ##                 restrictUp=10, restrictDown=200, minPathSize=10, 
#' ##                 supervisedStage1, typePCA='regular', 
#' ##                 resample1='BS', resample2='CV', 
#' ##                 repeatA1=20, repeatA2=1, repeatB1=20, repeatB2=1, 
#' ##                 nfolds=10, FSmethod1=NULL, FSmethod2=NULL, 
#' ##                 cutP1=0.1, cutP2=0.1, fdr2=NULL, FScore=param1, 
#' ##                 classifier, predMode, paramlist, innerCore=param2)


BioMM <- function(trainData, testData, pathlistDB, featureAnno, 
    restrictUp, restrictDown, minPathSize, supervisedStage1 = TRUE, 
    typePCA, resample1 = "BS", resample2 = "CV",  
    repeatA1 = 100, repeatA2 = 1, repeatB1 = 20, repeatB2 = 1, 
    nfolds = 10, FSmethod1, FSmethod2, 
    cutP1, cutP2, fdr2, FScore = MulticoreParam(), classifier, 
    predMode, paramlist, innerCore = MulticoreParam()) {

    trainDataList <- omics2pathlist(data=trainData, pathlistDB, 
                                    featureAnno, restrictUp, 
                                    restrictDown, minPathSize)  
    if (!is.null(testData)){ 
        testDataList <- omics2pathlist(data=testData, pathlistDB, 
                                       featureAnno, restrictUp, 
                                       restrictDown, minPathSize) 
    } else{
        testDataList <- NULL
    }
    
    ## generation of stage-2 data
    if (supervisedStage1 == TRUE) {
        stage2data <- BioMMreconData(trainDataList = trainDataList, 
            testDataList = testDataList,  resample = resample1, dataMode, 
            repeatA = repeatA1, repeatB = repeatB1, nfolds, 
            FSmethod = FSmethod1, cutP = cutP1, fdr = NULL, FScore, 
            classifier = classifier, predMode = predMode, 
            paramlist = paramlist, innerCore, 
            outFileA = NULL, outFileB = NULL)
    } else {
        stage2data <- BioMMstage1pca(trainDataList = trainDataList, 
            testDataList = testDataList, typeMode = typePCA, topPC = 1, 
            innerCore, outFileA = NULL, outFileB = NULL)
    }
    
    if (is.null(testDataList)) {
        trainData2 <- stage2data
        message("Stage-2: >>> ")
        message(paste0("Number of blocks: ", ncol(trainData2) - 1))
        trainPos2 <- getDataAfterFS(trainData = trainData2, testData = NULL, 
            FSmethod = "positive", cutP = 0.1, fdr = NULL, FScore)
        message(paste0("Number of positive blocks: ", ncol(trainPos2) - 1))
        testPos2 <- NULL
    } else {
        ## if testData provided
        trainData2 <- stage2data[[1]]
        testData2 <- stage2data[[2]]
        datalist <- getDataAfterFS(trainData2, testData2, FSmethod = "positive", 
            cutP = 0.1, fdr = NULL, FScore)
        ## include the label
        trainPos2 <- datalist[[1]]
        testPos2 <- datalist[[2]]
        message(paste0("Number of positive blocks: ", ncol(trainPos2) - 1))
    }
    
    ## If no positive features
    if (is.null(trainPos2)) {
        message("Warning: no positive features!!")
        result <- data.frame(pv = 1, cor = 0, AUC = 0.5, ACC = 0.5, R2 = 0)
    } else {
        ## make prediction
        result <- BioMMstage2pred(trainData = trainPos2, testData = testPos2, 
            resample = resample2, dataMode, repeatA = repeatA2, 
            repeatB = repeatB2, nfolds, FSmethod = FSmethod2, cutP = cutP2, 
            fdr = fdr2, FScore, classifier = classifier, 
            predMode = predMode, paramlist = paramlist, innerCore)
    }
    
}
