# File : resultReport.R Author : Junfang Chen




############################################################################### 
#' Compute the classification accuracy

#' @description
#' Compute the classification accuracy for the binary classification problem.

#' @param dataY The observed outcome.
#' @param predY The predicted outcome. 

#' @return The classification accuracy in terms of percentage.
#' @export 
#' @author Junfang Chen 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile) 
#' dataY <- methylData[,1]
#' methylSub <- data.frame(label=dataY, methylData[,c(2:1001)])  
#' library(ranger)  
#' library(BiocParallel)
#' param1 <- MulticoreParam(workers = 1) 
#' param2 <- MulticoreParam(workers = 10) 
#' predY <- predByCV(methylSub, repeats=1, nfolds=10,   
#'                   FSmethod=NULL, cutP=0.1, 
#'                   fdr=NULL, FScore=param1, 
#'                   classifier='randForest',
#'                   predMode='classification', 
#'                   paramlist=list(ntree=300, nthreads=1),
#'                   innerCore=param2)  
#' accuracy <- classifiACC(dataY=dataY, predY=predY)
#' print(accuracy)  


classifiACC <- function(dataY, predY) {
    
    tab <- table(dataY, predY)
    num1 <- sum(diag(tab))
    denom1 <- sum(tab)
    signif(num1/denom1, 3)
}


#' Compute the evaluation metrics

#' @description
#' Compute the evaluation metrics in the classification setting: 
#' area under curve (AUC), classification accuracy (ACC) and 
#' the pseudo R square (R2).

#' @param dataY The observed outcome.
#' @param predY The predicted outcome.
#' @details If all samples are predicted into one class, then R2 is 0.

#' @return A set of metrics for model evaluation: AUC, ACC and R2.
#' @export 
#' @import rms
#' @import glmnet   
#' @author Junfang Chen 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)   
#' dataY <- methylData[,1]
#' methylSub <- data.frame(label=dataY, methylData[,c(2:1001)])  
#' library(ranger) 
#' library(rms)
#' library(BiocParallel) 
#' param1 <- MulticoreParam(workers = 1) 
#' param2 <- MulticoreParam(workers = 10)  
#' predY <- predByCV(methylSub, repeats=1, nfolds=10,   
#'                   FSmethod=NULL, cutP=0.1, 
#'                   fdr=NULL, FScore=param1, 
#'                   classifier='randForest',
#'                   predMode='classification', 
#'                   paramlist=list(ntree=300, nthreads=20),
#'                   innerCore=param2)   
#' accuracy <- getMetrics(dataY=dataY, predY=predY)
#' print(accuracy)  

getMetrics <- function(dataY, predY){
    
    auc <- roc(dataY, predY)$auc      
    predY <- ifelse(predY>=.5, 1, 0) 
    cat("\n Levels of predicted Y =", nlevels(factor(predY)),"\n\n") 
    ACC <- classifiACC(dataY, predY)  
    if (nlevels(factor(predY)) > 1){
        R2 <- lrm(dataY ~ predY)$stats["R2"]
    } else {
        R2 <- 0
    }
    eMat <- data.frame(AUC=round(auc,3), ACC=round(ACC,3), R2 = round(R2, 3))  
    print(eMat) 
    return(eMat)
} 


############################################################################### 

#' Plot data summary statistics 

#' @description Plot data summary statistics in terms of the proportion of 
#' variance explained.

#' @param data The input dataset (either data.frame or matrix). 
#' Rows are the samples, columns are the probes/genes, except that 
#' the first column is the label (the outcome).
#' @param posF A logical value indicating if only positively outcome-associated
#' features should be used. (Default: TRUE)
#' @param stratify A string. Pathway based stratification method to generate 
#' \code{blocklist}. 
#' @param core The number of cores used for computation. (Default: 1)
#' @param fileName The file name specified for the plot. If it is not NULL,
#' then the plot will be generated. The plot will project the data on the 
#' first two components. (Default: 'R2explained.png') 
#' @return An output image file with '.png' format.  
#' @references Yu, Guangchuang, et al. 'clusterProfiler: an R package for 
#' comparing biological themes among gene clusters.' Omics: a journal of 
#' integrative biology 16.5 (2012): 284-287. 
#' @references Perlich, C., & Swirszcz, G. (2011). On cross-validation and  
#' stacking: Building seemingly predictive models on random data. ACM SIGKDD    
#' Explorations Newsletter, 12(2), 11-15. 

#' @import BiocParallel
#' @import variancePartition
#' @import grDevices
#' @export  


plotVarExplained <- function(data, posF = TRUE, stratify = c("pathway"), 
    core = MulticoreParam(), fileName = NULL) {
    
    if (colnames(data)[1] != "label") {
        stop("The first column of the 'data' must be the 'label'!")
    }
    dataX <- data[, -1]
    ## convert prob. to integer
    dataX <- apply(dataX, 2, round)
    dataY <- data[, 1]
    if (is.factor(dataY)) {
        dataY <- as.numeric(dataY) - 1
    }
    if (posF) {
        corr <- cor(dataX, dataY)
        nPos <- sum(corr > 0)
        message(paste0("posFeature: ", nPos))
        if (nPos == 0) {
            stop("No positively outcome-associated features!")
        }
        ## 'NA' may appear, use is.element to avoid NA.
        dataXsub <- dataX[, is.element(corr > 0, TRUE)]
    } else {
        dataX <- dataXsub
    }
    featurelist <- seq_len(ncol(dataXsub))

    r2mat <- unlist(bplapply(featurelist, function(i) {
        r2 <- lrm(dataXsub[, i] ~ dataY)$stats["R2"]
    }, BPPARAM = core))
    
    r2plot <- data.frame(Stage2data = r2mat)
    stratify <- match.arg(stratify)
    colnames(r2plot) <- paste0("Reconstructed ", stratify, "s")
    rownames(r2plot) <- colnames(dataXsub)
    head(r2plot)
    if (is.null(fileName)) {
        fileName <- "R2explained.png"
    }
    png(fileName, width = 5, height = 6, units = "in", res = 330)
    print(plotVarPart(r2plot, label.angle = 10, ylab = "Variance explained (%)", 
        convertToPercent = FALSE))
    dev.off()
}


############################################################################### 

#' Plot top outcome-associated features

#' @description Plot top ranked outcome-associated features from stage-2 data. 
#' The ranking criteria are based on metrics such as Nagelkerke pseudo 
#' R-square. 
#' @param data The input stage-2 data (either data.frame or matrix). 
#' Rows are the samples, columns are pathway names,  
#' except that the first column is the label (the outcome).
#' @param posF A logical value indicating if only positively outcome-associated
#' features should be used. (Default: TRUE)
#' @param topF The top ranked number of features at stage-2 (\code{topF} >= 2).
#' (Default: 10)
#' @param blocklist A list of matrices with block IDs as the associated list 
#' member names. The block IDs identical to the stage-2 feature names.  
#' For each matrix, rows are the samples and columns are the probe names,  
#' except that the first column is named 'label'. See also 
#' \code{\link{omics2pathlist}}.
#' @param stratify A string. Pathway based stratification method to generate 
#' \code{blocklist}. 
#' @param rankMetric A string representing the metrics used for ranking. 
#' Valid options are c('AUC', 'ACC', 'R2', 'size').
#' 'size' is the block size.
#' @param colorMetric A string representing the metric used to color the plot. 
#' Valid options are c('AUC', 'ACC', 'R2', 'size').
#' 'size' is the block size.
#' @param core The number of cores used for computation. (Default: 10)
#' @param fileName The plot file name. (Default: 'plottopF.png') 

#' @return An output image file. 

#' @details If the argument \code{posF} is TRUE, 
#' and no positively outcome-associated features are present in stage-2 data 
#' , then an error is reported. In addition, if \code{topF} is bigger than
#' the number of positively outcome-associated features, an error is returned. 

#' @references Perlich, C., & Swirszcz, G. (2011). On cross-validation and  
#' stacking: Building seemingly predictive models on random data. ACM SIGKDD    
#' Explorations Newsletter, 12(2), 11-15. 

#' @import BiocParallel 
#' @import lattice  
#' @import ggplot2 
#' @export  
#' @seealso  \code{\link{omics2pathlist}}.


plotRankedFeature <- function(data, posF = TRUE, topF = 10, blocklist, 
    stratify = c("pathway"), 
    rankMetric = c("AUC", "ACC", "R2", "size"), 
    colorMetric = c("AUC", "ACC", "R2", "size"), 
    core = MulticoreParam(), fileName = NULL) {
    
    .getBlockSize <- function(blocklist, stratify = c("pathway")) {
    
        stratify <- match.arg(stratify)
        ID <- gsub("\\:", ".", names(blocklist))
        size <- unlist(lapply(blocklist, function(d) {
            ncol(d) - 1
        }))
        blockSize <- data.frame(ID, size, stringsAsFactors = FALSE)        
        return(blockSize)
    }
    
    
    if (colnames(data)[1] != "label") {
        stop("The first column of the 'data' must be the 'label'!")
    }
    dataX <- data[, -1]
    ## convert prob. to integer
    dataX <- apply(dataX, 2, round)
    dataY <- data[, 1]
    if (is.factor(dataY)) {
        dataY <- as.numeric(dataY) - 1
    }
    
    if (posF) {
        corr <- cor(dataX, dataY)
        nPos <- sum(corr > 0)
        message(paste0("posFeature: ", nPos))
        if (nPos == 0) {
            stop("No positively outcome-associated features!")
        }
        if (topF > nPos) {
            stop("'topF' bigger than # of positively associated features!")
        }
        ## 'NA' may appear, use is.element to avoid NA.
        dataXsub <- dataX[, is.element(corr>0, TRUE)] 
        featurelist <- as.list(seq_len(ncol(dataXsub)))
        metrics <- unlist(bplapply(featurelist, function(i) {
            invisible(capture.output(eMat <- getMetrics(dataXsub[, i], dataY)))
            eMat
        },  BPPARAM = core))
    } else {
        featurelist <- as.list(seq_len(ncol(dataX)))
        metrics <- unlist(bplapply(featurelist, function(i) {
            invisible(capture.output(eMat <- getMetrics(dataX[, i], dataY)))
            eMat
        }, BPPARAM = core))
        dataXsub <- dataX
    }
    
    eMat <- matrix(unlist(metrics), nrow = ncol(dataXsub), byrow = TRUE)
    colnames(eMat) <- c("AUC", "ACC", "R2")
    if (stratify == "pathway") {
        goID <- gsub("\\:", ".", colnames(dataXsub))
        rownames(eMat) <- goID
    } else {
        rownames(eMat) <- colnames(dataXsub)
    }
    ## checking the 'blocklist'
    blockSize <- .getBlockSize(blocklist, stratify)
    ## double check the overlapping IDs
    sharedID <- intersect(rownames(eMat), blockSize[, "ID"])
    eMatSub <- eMat[is.element(rownames(eMat), sharedID), ]
    eMat2 <- eMatSub[match(rownames(eMatSub), sharedID), ]
    blockSub <- blockSize[is.element(blockSize[, "ID"], sharedID), ]
    blockMatch <- blockSub[match(rownames(eMatSub), blockSub[, "ID"]), ]
    ## attached block Size
    blockInfo <- data.frame(eMat2, blockMatch, stringsAsFactors = FALSE)
    ## ranking
    rankMetric <- match.arg(rankMetric)
    blockInfo2 <- blockInfo[order(blockInfo[, rankMetric], decreasing = TRUE),]
    topPat <- head(blockInfo2, topF)
    topPat$ID <- factor(topPat$ID, levels = rev(unique(topPat$ID))) 
    x <- "ID"
    y <- rankMetric 
    colorby <- match.arg(colorMetric)
    stratify <- match.arg(stratify) 
    if (stratify == "pathway") {
        subtitle <- "Pathways"
    } 
    title <- paste0("Top ", topF, " ", subtitle)
    if (is.null(fileName)) {
        fileName <- paste0("plotTopF", topF, "_", rankMetric, "_", stratify, 
            ".png")
    }
    png(fileName, width = 5, height = 6, units = "in", res = 330)
    print(ggplot(topPat, aes_string(x = x, y = y, fill = colorby)) + 
        scale_fill_continuous(low = "red", high = "blue", name = colorby, 
            guide = guide_colorbar(reverse = TRUE)) + 
        geom_bar(stat = "identity") + coord_flip() + ggtitle(title) + 
        xlab(NULL) + ylab(y))
    dev.off()
}
