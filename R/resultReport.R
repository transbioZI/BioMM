# File : resultReport.R   Author : Junfang Chen




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


#' Compute the machine learning evaluation metrics

#' @description
#' Compute the evaluation metrics in the classification setting: 
#' area under curve (AUC), the area under the Precision-Recall curve, 
#' classification accuracy (ACC) and the pseudo R square (R2).

#' @param dataY The observed outcome.
#' @param predY The predicted outcome.
#' @details If all samples are predicted into one class, then R2 is 0.

#' @return A set of metrics for model evaluation: AUC, AUCPR, ACC and R2.
#' @export 
#' @import rms
#' @import precrec
#' @author Junfang Chen 
#' @examples  
#' ## Load data  
#' methylfile <- system.file('extdata', 'methylData.rds', package='BioMM')  
#' methylData <- readRDS(methylfile)   
#' dataY <- methylData[,1]
#' methylSub <- data.frame(label=dataY, methylData[,c(2:1001)])  
#' library(ranger) 
#' library(precrec)
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
#' metrics <- getMetrics(dataY=dataY, predY=predY)
#' print(metrics)  

 

getMetrics <- function(dataY, predY){
    
    sscurves <- evalmod(scores=predY, labels=dataY)  
    AUC <- attr(sscurves[[1]][[1]], "auc")
    AUCPR <- attr(sscurves[[2]][[1]], "auc") 
    predYbinary <- ifelse(predY>=.5, 1, 0) 
    ACC <- classifiACC(dataY, predYbinary)  
    if (nlevels(factor(predY)) > 1){
        R2 <- lrm(dataY ~ predY)$stats["R2"]
    } else {
        R2 <- 0
    }
    eMat <- data.frame(AUC=round(AUC,3), AUCPR=round(AUCPR,3), 
                       ACC=round(ACC,3), R2 = round(R2, 3))  

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
#' @param binarize A logical value indicating if the individual features under
#' investigation should be binarized. The default is FALSE, which provides the 
#' estimated class probabilities for each pathway-level feature. If TRUE, then
#' the binary output is given for each feature.
#' @param core The number of cores used for computation. (Default: 1)
#' @param pathTitle A string indicating the name of pathway under investigation.
#' This will be displayed as the name of y-axis.
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
#' @import grDevices
#' @import vioplot
#' @export  


plotVarExplained <- function(data, posF = TRUE, binarize = FALSE, 
    core = MulticoreParam(), pathTitle = "GO pathways", fileName = NULL) {
    
    if (colnames(data)[1] != "label") {
        stop("The first column of the 'data' must be the 'label'!")
    }
    dataX <- data[, -1]
    if (binarize == TRUE){
        ## convert prob. to integer
        dataX <- apply(dataX, 2, round)  
    }
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
        r2 <- lrm(dataY ~ dataXsub[,i])$stats["R2"]
    }, BPPARAM = core))
    
    r2plot <- data.frame(Stage2data = r2mat) 
    colnames(r2plot) <- pathTitle
    rownames(r2plot) <- colnames(dataXsub)
    head(r2plot)
    if (is.null(fileName)) {
        fileName <- "R2explained.png"
        vioplot(r2plot, names=pathTitle, main = "", col="#F8766D", 
            horizontal=TRUE, xlab="Variance explained (%)",  areaEqual=TRUE)
    } else {
        png(fileName, width = 7, height = 5, units = "in", res = 330)
        vioplot(r2plot, names=pathTitle, main = "", col="#F8766D", 
            horizontal=TRUE, xlab="Variance explained (%)",  areaEqual=TRUE)
        dev.off()
    }
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
#' @param binarize A logical value indicating if the individual features under
#' investigation should be binarized. The default is FALSE, which provides the 
#' estimated class probabilities for each pathway-level feature. If TRUE, then
#' the binary output is given for each feature.
#' @param rankMetric A string representing the metrics used for ranking. 
#' Valid options are c("AUC", "R2", "Zscore", "negPlogit", "negPwilcox").
#' "negPlogit" denotes the negative log P value from the logistic regression 
#' and "negPwilcox" means the negative log P value based on the Wilcoxon test. 
#' "size" is the block size.
#' @param colorMetric A string representing the metric used to color the plot. 
#' Valid options are c("AUC", "R2", "Zscore", "negPlogit", "negPwilcox").
#' "negPlogit" denotes the negative log P value from the logistic regression 
#' and "negPwilcox" means the negative log P value based on wilcoxon test. 
#' "size" is the block size.
#' @param core The number of cores used for computation. (Default: 10)
#' @param pathTitle A string indicating the name of pathway under investigation.
#' @param fileName The plot file name. (Default: 'plottopF.png') 

#' @return An output image file and the summary statistics of the top pathways. 

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


plotRankedFeature <- function(data, posF = TRUE, topF = 10, blocklist, binarize=FALSE, 
    rankMetric = c("AUC", "R2", "Zscore", "negPlogit", "negPwilcox", "size"), 
    colorMetric = c("AUC", "R2", "Zscore", "negPlogit", "negPwilcox", "size"),
    core = MulticoreParam(), pathTitle = "GO pathways", fileName = NULL) {
    
    .getBlockSize <- function(blocklist) {
     
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
    if (binarize == TRUE){
        dataX <- apply(dataX, 2, round)    
    }
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
    } else {  
        dataXsub <- dataX
    }

    featurelist <- as.list(seq_len(ncol(dataXsub)))
    metrics <- unlist(bplapply(featurelist, function(i) { 
        predY <- dataXsub[,i]
        eMatTmp <- getMetrics(dataY, predY)
        eMatTmp <- eMatTmp[c(1,4)]
        fit <- glm(dataY~predY, family="binomial")
        zPvMat <- summary(fit)$coefficients[,3:4]
        zPv <- zPvMat[is.element(rownames(zPvMat), "predY"),]   
        names(zPv) <- c("Zscore", "negPlogit")
        Zscore <- round(zPv[1], 3)
        negPlogit <- round(-log10(zPv[2]), 3)
        pWilcox <- wilcox.test(predY ~ dataY)$p.value  
        negPwilcox <- round(-log10(pWilcox), 3) 
        metrics <- unlist(c(eMatTmp, Zscore, negPlogit, negPwilcox=negPwilcox))
        metrics
    },  BPPARAM = core))
 
    eMat <- matrix(unlist(metrics), nrow = ncol(dataXsub), byrow = TRUE) 
    colnames(eMat) <- names(metrics)[seq_len(ncol(eMat))]
    goID <- gsub("\\:", ".", colnames(dataXsub))
    rownames(eMat) <- goID 
    ## checking the 'blocklist'
    blockSize <- .getBlockSize(blocklist)
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
    title <- paste0("Top ", topF, " ", pathTitle)
    if (is.null(fileName)) {
        fileName <- paste0("plotTopF", topF, "_", rankMetric, "_pathway.png")
    }
    png(fileName, width = 5, height = 6, units = "in", res = 330)
    p <- ggplot(topPat, aes_string(x = x, y = y, fill = colorby)) + 
        scale_fill_continuous(low = "gray", high = "#F8766D", name = colorby, 
            guide = guide_colorbar(reverse = TRUE)) + 
        geom_bar(stat = "identity") + coord_flip() + ggtitle(title) + 
        xlab(NULL) + ylab(y)
    if (rankMetric == "negPlogit" | rankMetric == "negPwilcox"){
        p <- p + geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red")
        p <- p + labs(y = "-logP")   
    }    
    print(p)
    dev.off()
 
    return(topPat)
}


############################################################################### 

#' Circular plot for a set of pathways

#' @description The individual CpGs or genes within a given set of pathways are  
#' displayed as the dots in the resulting plot. The significance of the CpGs or genes   
#' are illustrated by the negative log P value.

#' @param datalist The input data list containing ordered 
#' collections of matrices.
#' @param topPathID A predefined pathway IDs.  
#' @param core The number of cores used for computation. (Default: 10)
#' @param fileName Add a character to the output file name. 
#' (Default: 'Circular-Manhattan.pval.jpeg') 

#' @return An output image file. 

#' @details Top 10 or 20 pathways are usually suggested to be visualized.
#' The significant features (if any) are highlighted using filled diamond. 
#' The significance line is set as 0.05 marked as dashed red line.
#' @import BiocParallel
#' @import CMplot
#' @export
#' @seealso  \code{\link{omics2pathlist}}.


cirPlot4pathway <- function(datalist, topPathID, core = MulticoreParam(), fileName = NULL){
  
    sublistV0 <- datalist[is.element(names(datalist), topPathID)]  
    sublist <- sublistV0[match(topPathID, names(sublistV0))]   
    message(paste0("Checking top ", length(sublist), " pathways >> "))

    cpgPvByGO <- list()
    cpg <- c()
    pathway <- c()
    pos <- c()
    pval <- c()

    for (i in seq_along(sublist)){
        print(i)
        data <- sublist[[i]] 
        dataX <- data[,-1]
        dataY <- data[,1]
        cpgTmp <- colnames(dataX)
        pathwayTmp <- rep(i, ncol(dataX))
        posTmp <- seq_len(ncol(dataX)) 

        featurelist <- as.list(seq_len(ncol(dataX)))
        pvalTmp <- unlist(bplapply(featurelist, function(i){ 
            # wilcox.test(dataX[,i] ~ dataY)$p.value  
            predY <- dataX[,i] 
            fit <- glm(dataY~predY, family="binomial")
            zPvMat <- summary(fit)$coefficients[,3:4]
            zPv <- zPvMat[is.element(rownames(zPvMat), "predY"),]   
            pval <- zPv[2] 
        }, BPPARAM = core))    
        names(pvalTmp)  = colnames(dataX) 

        cpgPvByGO[[i]] <- pvalTmp
        cpg <- c(cpg, cpgTmp)
        pathway <- c(pathway, pathwayTmp)
        pos <- c(pos, posTmp)
        pval <- c(pval, pvalTmp)   
    } 
 
    cpgResults <- data.frame(cpg, pathway, pos, pval, stringsAsFactors=FALSE)
   
    ## -log10 was used
    CMplot(cpgResults, plot.type="c", r=2, 
            cir.legend=TRUE, cex.axis=0.9,
            outward=TRUE, cir.legend.col="black", cir.chr.h=0.8,
            threshold=0.05, signal.pch=18,
            chr.den.col="white", file="jpg",
            memo=fileName, dpi=400, chr.labels=names(sublist))
 
 }