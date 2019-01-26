# File   : 4.resultReport.R
# Author : Junfang Chen




############################################################################### 
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
#' methylfile <- system.file("extdata", "methylData.rds", package="BioMM")  
#' methylData <- readRDS(methylfile) 
#' dataY <- methylData[,1]
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' library(ranger) 
#' predY <- predByCV(methylSub, repeats=1, nfolds=10,   
#'                   FSmethod=NULL, cutP=0.1, 
#'                   fdr=NULL, FScore=1, 
#'                   classifier="randForest",
#'                   predMode="classification", 
#'                   paramlist=list(ntree=300, nthreads=1),
#'                   innerCore=1)  
#' accuracy <- classifiACC(dataY=dataY, predY=predY)
#' print(accuracy)  


classifiACC <- function(dataY, predY){ 

    tab <- table(dataY, predY)
    num1 <- sum(diag(tab))
    denom1 <- sum(tab)
    signif(num1/denom1,3)
}


#' Compute the evaluation metrics

#' @description
#' Compute the evaluation metrics in the classification setting: 
#' P value based on chi-square test (pv), pearson correlation coefficient 
#' (cor), area under curve (AUC), classification accuracy (ACC) and 
#' the pseudo R square (R2).

#' @param dataY The observed outcome.
#' @param predY The predicted outcome.
#' @details If all samples are predicted into one class, then we assign 
#' R2=0, cor=0, and AUC=0.5.

#' @return A set of metrics for model evaluation: pv, cor, AUC, ACC and R2.
#' @export 
#' @import rms
#' @import glmnet   
#' @author Junfang Chen 
#' @examples  
#' ## Load data  
#' methylfile <- system.file("extdata", "methylData.rds", package="BioMM")  
#' methylData <- readRDS(methylfile)   
#' dataY <- methylData[,1]
#' methylSub <- data.frame(label=dataY, methylData[,c(2:2001)])  
#' library(ranger) 
#' library(rms)
#' predY <- predByCV(methylSub, repeats=1, nfolds=10,   
#'                   FSmethod=NULL, cutP=0.1, 
#'                   fdr=NULL, FScore=1, 
#'                   classifier="randForest",
#'                   predMode="classification", 
#'                   paramlist=list(ntree=300, nthreads=1),
#'                   innerCore=1)   
#' accuracy <- getMetrics(dataY=dataY, predY=predY)
#' print(accuracy)  


getMetrics <- function(dataY, predY){
    
    cat("\n Levels of predicted Y =", nlevels(factor(predY)),"\n\n") 
    ACC <- classifiACC(dataY, predY) 
    pv <- chisq.test(table(dataY, predY))$p.value 
    if (nlevels(factor(predY)) > 1){     
        Cor <- cor(predY, dataY)  
        AUC <- auc(predY, dataY)   
        R2 <- lrm(predY~dataY)$stats["R2"] 
        eMat <- data.frame(pv=pv, cor=round(Cor,2),
                    AUC=round(AUC,2), ACC=round(ACC,2), R2=round(R2,3))  
    } else {
        message("Warning: all predicted samples in one class!") 
        eMat <- data.frame(pv=pv, cor=0, AUC=0.5, ACC=round(ACC,2), R2=0) 
    }   
    print(eMat) 
    return(eMat)
} 


###############################################################################
###############################################################################

#' Plot data summary statistics 

#' @description Plot data summary statistics in terms of the proportion of 
#' variance explained.

#' @param data The input dataset (either data.frame or matrix). 
#' Rows are the samples, columns are the probes/genes, except that 
#' the first column is the label (the outcome).
#' @param posF A logical value indicating if only positively outcome-associated
#' features should be used. (Default: TRUE)
#' @param stratify A string. The applied stratification method to generate 
#' \code{blocklist}. Valid options are c("gene", "pathway", "chromosome").
#' @param core The number of cores used for computation. (Default: 1)
#' @param fileName The file name specified for the plot. If it is not NULL,
#' then the plot will be generated. The plot will project the data on the 
#' first two components. (Default: "R2explained.png") 
#' @return An output image file with ".png" format.  
#' @references Yu, G., Wang, L.G., Dall'Olio, G., Yu, M.G. and GOSemSim, A., 
#' 2018. Package ‘clusterProfiler’.
#' @references Claudia Perlich and Grzegorz Swirszcz. On cross-validation and 
#' stacking: Building seemingly predictive models on random data. ACM SIGKDD 
#' Explorations Newsletter, 12(2):11–15, 2011.

#' @import BiocParallel
#' @import variancePartition
#' @import grDevices
#' @export  


plotVarExplained <- function(data, posF=TRUE,
                            stratify=c("gene", "pathway", "chromosome"),
                            core=1, fileName=NULL){ 

    if ( colnames(data)[1] != "label" ){
        stop("The first column of the 'data' must be the 'label'!") 
    }
    dataX <- data[,-1]
    ## convert prob. to integer
    dataX <- apply(dataX, 2, round) 
    dataY <- data[,1]  
    if (is.factor(dataY)){dataY <- as.numeric(dataY)-1}   
    if (posF){ 
        corr <- cor(dataX, dataY)
        nPos <- length(which(corr > 0))
        message(paste0("posFeature: ", nPos))
        if (nPos == 0){
            stop("No positively outcome-associated features!")    
        }
        ## use 'which' to avoid possible 'NAs'
        dataXsub <- dataX[,which(corr > 0)] 
    } else { 
        dataX <- dataXsub
    }       
    featurelist <- seq_len(ncol(dataXsub))    
    r2mat <- unlist(bplapply(featurelist, function(i){  
            r2 <-  lrm(dataXsub[,i]~dataY)$stats["R2"]  
    }, BPPARAM=SnowParam(workers=core)))  

    r2plot <- data.frame(Stage2data=r2mat)
    stratify <- match.arg(stratify) 
    colnames(r2plot) <- paste0("Reconstructed ", stratify, "s")
    rownames(r2plot) <- colnames(dataXsub)
    head(r2plot) 
    if (is.null(fileName)){fileName <- "R2explained.png"}
    png(fileName, width=5, height=6, units='in', res=330)  
    print(plotVarPart(r2plot, label.angle=10,
                        ylab="Variance explained (%)",  
                        convertToPercent=FALSE))
    dev.off() 
}


###############################################################################
###############################################################################

#' Plot top outcome-associated features

#' @description Plot top ranked outcome-associated features from stage-2 data. 
#' The ranking criteria are based on metrics such as Nagelkerke pseudo 
#' R-square. 
#' @param data The input stage-2 data (either data.frame or matrix). 
#' Rows are the samples, columns are gene IDs, or pathway names or chromosome 
#' IDs, except that the first column is the label (the outcome).
#' @param posF A logical value indicating if only positively outcome-associated
#' features should be used. (Default: TRUE)
#' @param topF The top ranked number of features at stage-2 (\code{topF} >= 2).
#' (Default: 10)
#' @param blocklist A list of matrices with block IDs as the associated list 
#' member names. The block IDs identical to the stage-2 feature names. 
#' The block can be gene, pathway or chromosome. 
#' For each matrix, rows are the samples and columns are the probe names,  
#' except that the first column is named "label". See also 
#' \code{\link{omics2genelist}}; \code{\link{omics2pathlist}}; 
#' \code{\link{omics2chrlist}}
#' @param stratify A string. The applied stratification method to generate 
#' \code{blocklist}. Valid options are c("gene", "pathway", "chromosome"). 
#' @param rankMetric A string representing the metrics used for ranking. 
#' Valid options are c("cor", "AUC", "ACC", "R2", "size").
#' "size" is the block size.
#' @param colorMetric A string representing the metric used to color the plot. 
#' Valid options are c("cor", "AUC", "ACC", "R2", "size").
#' "size" is the block size.
#' @param core The number of cores used for computation. (Default: 10)
#' @param fileName The plot file name. (Default: "plottopF.png") 

#' @return An output image file. 

#' @details If the argument \code{posF} is TRUE, 
#' and no positively outcome-associated features are present in stage-2 data 
#' , then an error is reported. In addition, if \code{topF} is bigger than
#' the number of positively outcome-associated features, an error is returned. 

#' @references Claudia Perlich and Grzegorz Swirszcz. On cross-validation and 
#' stacking: Building seemingly predictive models on random data. ACM SIGKDD 
#' Explorations Newsletter, 12(2):11–15, 2011.

#' @import BiocParallel 
#' @import lattice  
#' @import ggplot2 
#' @export  
#' @seealso \code{\link{omics2genelist}}; \code{\link{omics2pathlist}}; 
#' \code{\link{omics2chrlist}}


plotRankedFeature <- function(data, posF=TRUE, topF=10,
                            blocklist,
                            stratify=c("gene", "pathway", "chromosome"),
                            rankMetric=c("cor", "AUC", "ACC", "R2", "size"),
                            colorMetric=c("cor", "AUC", "ACC", "R2", "size"),
                            core=10, fileName=NULL){  

    .getBlockSize <- function(blocklist,
                            stratify=c("gene", "pathway", "chromosome")){

        stratify <- match.arg(stratify)
        ## blocks of sub-datasets preparation
        if (stratify == "gene" || stratify == "chromosome"){  
            ID <- names(blocklist)
            ## exclude the label (minus 1)
            size <- unlist(lapply(blocklist, function(d){ncol(d)-1}))
            blockSize <- data.frame(ID, size, stringsAsFactors=FALSE) 
        } else if (stratify == "pathway") {
            ID <- gsub("\\:", ".", names(blocklist)) 
            size <- unlist(lapply(blocklist, function(d){ncol(d)-1}))
            blockSize <- data.frame(ID, size, stringsAsFactors=FALSE)  
        } else {
            stop("Wrong stratification method.")
        }  
        
        return(blockSize)
    }


    if ( colnames(data)[1] != "label" ){
        stop("The first column of the 'data' must be the 'label'!") 
    }
    dataX <- data[,-1]
    ## convert prob. to integer
    dataX <- apply(dataX, 2, round) 
    dataY <- data[,1]  
    if (is.factor(dataY)){dataY <- as.numeric(dataY)-1}  
    if (posF){ 
        corr <- cor(dataX, dataY)
        nPos <- length(which(corr > 0))
        message(paste0("posFeature: ", nPos))
        if (nPos == 0){
            stop("No positively outcome-associated features!")
        }
        if (topF > nPos ){
            stop("'topF' bigger than # of positively associated features!")
        }
        ## use 'which' to avoid possible 'NAs'
        dataXsub <- dataX[,which(corr > 0)]
        featurelist <- as.list(seq_len(ncol(dataXsub))) 
        metrics <- unlist(bplapply(featurelist, function(i){   
                eMat <- getMetrics(dataXsub[,i], dataY)
        }, BPPARAM=SnowParam(workers=core)))  
    } else { 
        featurelist <- as.list(seq_len(ncol(dataX))) 
        metrics <- unlist(bplapply(featurelist, function(i){   
                eMat <- getMetrics(dataX[,i], dataY)
        }, BPPARAM=SnowParam(workers=core))) 
        dataXsub <- dataX
    }    
    
    eMat <- matrix(unlist(metrics), nrow=ncol(dataXsub), byrow=TRUE)
    met1 <- getMetrics(dataXsub[,1], dataY)
    colnames(eMat) <- colnames(met1)
    if (stratify == "pathway"){ 
        goID <- gsub("\\:", ".", colnames(dataXsub))
        rownames(eMat) <- goID
    } else {
        rownames(eMat) <- colnames(dataXsub)
    }  
    ## checking the 'blocklist'
    blockSize <- .getBlockSize(blocklist, stratify)
    ## double check the overlapping IDs
    sharedID <- intersect(rownames(eMat), blockSize[,"ID"])
    eMat2 <- eMat[match(rownames(eMat), sharedID), ]
    blockMatch <- blockSize[match(rownames(eMat), sharedID), ] 
    ## attached block Size  
    blockInfo <- data.frame(eMat2, blockMatch, stringsAsFactors=FALSE)
    ## ranking 
    rankMetric <- match.arg(rankMetric)
    print(rankMetric)
    blockInfo2 <- blockInfo[order(blockInfo[, rankMetric], decreasing=TRUE),]  
    topPat <- head(blockInfo2, topF)  
    topPat$ID <- factor(topPat$ID, levels=rev(unique(topPat$ID))) 
    print(head(topPat)) 
    x <- "ID"
    y <- rankMetric
    print(y)
    colorby <- match.arg(colorMetric) 
    stratify <- match.arg(stratify) 
    if (stratify == "gene"){subtitle <- "Genes"}
    if (stratify == "pathway"){subtitle <- "Pathways"} 
    if (stratify == "chromosome"){subtitle <- "Chromosomes"} 
    title <- paste0("Top ", topF, " ", subtitle)
    if (is.null(fileName)){
        fileName <- paste0("plotTopF", topF, "_",
                        rankMetric, "_", stratify, ".png")
    }
    png(fileName, width=5, height=6, units='in', res=330) 
    print(ggplot(topPat, aes_string(x=x, y=y, fill=colorby)) +
        scale_fill_continuous(low="red", high="blue", name=colorby,
            guide=guide_colorbar(reverse=TRUE)) +
        geom_bar(stat="identity") + coord_flip() +
        ggtitle(title) + xlab(NULL) + ylab(y))
    dev.off() 
} 
