

############################################################################### 
###############################################################################

#' Map individual probes into gene

#' @description
#' Map a set of individual probes from different omics 
#' (i.e. SNPs and CpGs etc.) into genes.

#' @param data The input dataset (either data.frame or matrix). 
#' Rows are the samples, columns are the probes, except that 
#' the first column is the label.  
#' @param featureAnno The annotation data stored in a data.frame for probe 
#' mapping. It must have at least two columns named "ID" and "entrezID".  
#' @param restrictUp The upper-bound of the number of probes in each gene.
#' The default is 500.
#' @param restrictDown The lower-bound of the number of probes in each gene. 
#' The default is 5.
#' @return A list of matrices with gene IDs as the associated list member  
#' names. For each matrix, rows are the samples and columns are the probe 
#' names, except that the first column is named "label".
#' 
#' @details 
#' This function is not applicable for transcriptomic data. The data types  
#' including DNA methylation and GWAS often have several probes within one 
#' gene, therefore gene-based stratification is feasible 

#' @export 
#' @examples  
#' ## Load data from DNA methylation
#' methylfile <- system.file("extdata", "methylData.rds", package="BioMM")  
#' methylData <- readRDS(methylfile)  
#' ## Annotation files for Mapping CpGs into genes  
#' featureAnnoFile <- system.file("extdata", "cpgAnno.rds", package="BioMM")  
#' featureAnno <- readRDS(file=featureAnnoFile)  
#' ## Mapping CpGs into gene list
#' ## Not run 
#' ## dataList <- omics2genelist(data=methylData, featureAnno, 
#' ##                            restrictUp=100, restrictDown=5)
#' ## length(dataList)


omics2genelist <- function(data, featureAnno, restrictUp=500, restrictDown=5){

    if ( colnames(data)[1] != "label" ){
        stop("The first column of the 'data' must be the 'label'!")
    }
    dataX <- data[,-1] 
    dataY <- data[,1] 
    probeName <- colnames(dataX) 
    featureAnno <- as.data.frame(featureAnno, stringsAsFactors=FALSE) 
    colnameAnno <- colnames(featureAnno) 
    if ( length(grep("^ID", colnameAnno)) == 0 ){ 
        stop("No 'ID' column in 'featureAnno'! ")
    }
    if ( length(grep("^entrezID", colnameAnno)) == 0) {
        stop("No 'entrezID' column in 'featureAnno'! ")
    }
    probeAnno <- featureAnno[is.element(featureAnno[,'ID'], probeName),
                            c('ID', 'entrezID')]
    if ( nrow(probeAnno) == 0 ) {
        stop("Wrong 'featureAnno' provided!")
    }    
    genes <- unique(probeAnno[,"entrezID"])  
    matlist <- list()
    for (i in seq_along(genes)) {   
        annoSub <- probeAnno[probeAnno[,"entrezID"] == genes[i],]  
        IDperGene <- annoSub[,"ID"] 
        dataXsub <- dataX[,intersect(probeName, IDperGene)] 
        if (length(IDperGene) == 1){
            mat <- cbind(label=dataY, dataXsub)
            colnames(mat) <- c("label", IDperGene)
        } else { 
            mat <- cbind(label=dataY, dataXsub)
        }     
        matlist[[i]] <- mat 
    } 
    names(matlist) <- genes  
    ## exclude the label
    geneSize <- unlist(lapply(matlist, function(i) {ncol(i)-1}))  
    geneFilteredIndex <- geneSize >= restrictDown & geneSize <= restrictUp
    subGeneSize <- geneSize[geneFilteredIndex]
    matlistSub <- matlist[geneFilteredIndex]
    message("# Genes summary: ")  
    print(summary(geneSize))   
    message("# Filtered Genes summary: ") 
    print(summary(subGeneSize))   
    message(paste0('Total number of original genes: ', length(matlist))) 
    message(paste0('Total number of filtered genes: ', length(matlistSub)))  
    return(matlistSub) 
}


############################################################################### 
###############################################################################

#' Map individual probes into pathway

#' @description
#' Map a set of individual probes from different omics (i.e. SNPs, gene 
#' expression probes, CpGs etc.) into pathway such as Gene Ontology (GO) 
#' categories and KEGG.
#' @param data The input dataset (either data.frame or matrix). 
#' Rows are the samples, columns are the probes/genes, except that 
#' the first column is the label. If it's transcriptomic data, gene ID is 
#' the "entrezID".
#' @param pathlistDB A list of pathways with pathway IDs and their 
#' corresponding genes ("entrezID" is used). 
#' @param featureAnno The annotation data stored in a data.frame for probe 
#' mapping. It must have at least two columns named "ID" and "entrezID".  
#' If it's NULL, then the input probe is from transcriptomic data.  
#' @param restrictUp The upper-bound of the number of genes in each pathway. 
#' The default is 200.
#' @param restrictDown The lower-bound of the number of genes in each pathway. 
#' The default is 10.
#' @param minPathSize The minimal required number of probes in each pathway 
#' after mapping the input data to \code{pathlistDB}.
#' @return A list of matrices with pathway IDs as the associated list member
#' names. For each matrix, rows are the samples and columns are the probe   
#' names, except that the first column is named "label".
#' 
#' @details 
#' If gene expression \code{data} is the input, then \code{featureAnno} is
#' NULL,since the gene IDs are already defined as column names of the 
#' \code{data}. Since online database is updated from time to time, 
#' it is adivsed to make sure that the study database (e.g. \code{pathlistDB})  
#' is frozen at particular time for reproducing the results. 
#' The number of genes in each pathway can be restricted for downstream 
#' analysis because too small pathways are sparsely distributed, and too large 
#' pathways are often computationally intensive, and likely nonspecific. 

#' @export 
#' @examples  
#' ## Load data from DNA methylation
#' methylfile <- system.file("extdata", "methylData.rds", package="BioMM")  
#' methylData <- readRDS(methylfile)  
#' ## Annotation files for Mapping CpGs into pathways 
#' pathlistDBfile <- system.file("extdata", "goDB.rds", package="BioMM")
#' featureAnnoFile <- system.file("extdata", "cpgAnno.rds", package="BioMM") 
#' pathlistDB <- readRDS(file=pathlistDBfile)
#' featureAnno <- readRDS(file=featureAnnoFile)  
#' ## To reduce runtime
#' pathlistDB <- pathlistDB[1:20]
#' ## Mapping CpGs into pathway list 
#' dataList <- omics2pathlist(data=methylData, 
#'                                 pathlistDB, featureAnno, 
#'                                 restrictUp=100, restrictDown=20, 
#'                                 minPathSize=10)
#' length(dataList)


omics2pathlist <- function(data, pathlistDB, featureAnno=NULL,
                            restrictUp=200, restrictDown=10, minPathSize=2){

    if ( colnames(data)[1] != "label" ){ 
        stop("The first column of the 'data' must be the 'label'!")
    }
    dataX <- data[,-1] 
    dataY <- data[,1]
    # Restrict the pathways of size x-xx gene for downstream analysis   
    genePerPath <- lapply(pathlistDB, function(i){length(i)})  
    if (!is.null(restrictUp)) { 
        pathSizeIndex <- genePerPath >= restrictDown &
                            genePerPath <= restrictUp
    } else { 
        pathSizeIndex <- genePerPath >= restrictDown
    } 
    pathlistSub <- pathlistDB[pathSizeIndex]
    
    probeName <- colnames(dataX)   
    if (!is.null(featureAnno)){ 
        colnameAnno <- colnames(featureAnno) 
        if ( length(grep("^ID", colnameAnno)) == 0 ){ 
            stop("No 'ID' column in 'featureAnno'! ")
        }
        if ( length(grep("^entrezID", colnameAnno)) == 0){
            stop("No 'entrezID' column in 'featureAnno'! ")
        }
        probeAnno <- featureAnno[is.element(featureAnno[,'ID'], probeName),
                                c('ID', 'entrezID')]
        if ( nrow(probeAnno) == 0 ){ 
            stop("Wrong matching between 'data' and 'featureAnno'!")
        }    
    }
    pathlist <- list()
    for (i in seq_along(pathlistSub)) {  
        entrezIDpath <- pathlistSub[[i]]
        if (!is.null(featureAnno)){  
            probePerPath <- probeAnno[is.element(probeAnno[,'entrezID'],
                                        entrezIDpath), 'ID'] 
            dataXsub <- dataX[,intersect(probeName, probePerPath)]
        } else { ## for gene expression data 
            dataXsub <- dataX[,intersect(probeName, entrezIDpath)]
        }    
        pathMat <- cbind(label=dataY, dataXsub) 
        pathlist[[i]] <- pathMat 
    } 

    names(pathlist) <- names(pathlistSub)  
    ## exclude the label
    probeNperPath <- unlist(lapply(pathlist, function(i) {ncol(i)-1})) 
    # use 'which' to keep the remaining index
    minPathIndex <- which(probeNperPath < minPathSize) 
    ## remove pathways with too small size
    if (length(minPathIndex)!=0) {pathlist <- pathlist[-minPathIndex]}
    message("Summary stat of # probes in each mapped pathway: ")
    print(summary(unlist(lapply(pathlist, function(i) {ncol(i)-1} )))) 
    
    message(paste0('Total number of original pathways: ', length(pathlistDB))) 
    message(paste0('Total number of filtered original pathways: ',
                length(pathlistSub))) 
    message('Remove pathways of too small size if any')
    message(paste0('Total number of retained pathways: ', length(pathlist))) 
    return(pathlist) 
}


############################################################################### 
###############################################################################

#' Map individual probes into chromosome

#' @description
#' Map a set of individual probes from different omics (i.e. SNPs, gene 
#' expression probes, CpGs etc.) into the corresponding chromosomes.

#' @param data The input dataset (either data.frame or matrix). 
#' Rows are the samples, columns are the probes/genes, except that 
#' the first column is the label. 
#' @param probeAnno The annotation data stored in a data.frame for probe 
#' mapping. It must have at least two columns named "ID" and "chr".  
#' @return A list of matrices with chromosome IDs as the associated list 
#' member names. For each matrix, rows are the samples and columns are the  
#' probe names, except that the first column is named "label".
#' 
#' @export 
#' @examples  
#' ## Load data from DNA methylation
#' methylfile <- system.file("extdata", "methylData.rds", package="BioMM")  
#' methylData <- readRDS(methylfile)  
#' ## Annotation files for Mapping CpGs into chromosome  
#' probeAnnoFile <- system.file("extdata", "cpgAnno.rds", package="BioMM")  
#' probeAnno <- readRDS(file=probeAnnoFile)  
#' ## Mapping CpGs into Chromosome
#' dataList <- omics2chrlist(data=methylData, probeAnno)
#' length(dataList)


omics2chrlist <- function(data, probeAnno){

    if ( colnames(data)[1] != "label" ){
        stop("The first column of the 'data' must be the 'label'!")
    }
    dataX <- data[,-1] 
    dataY <- data[,1] 

    probeName <- colnames(dataX)    
    colnameAnno <- colnames(probeAnno) 
    if ( length(grep("^ID", colnameAnno)) == 0 ){
        stop("No 'ID' column in 'probeAnno'! ")
    }
    if ( length(grep("^chr", colnameAnno)) == 0){
        stop("No 'chr' column in 'probeAnno'! ")
    }
    probeAnno <- probeAnno[is.element(probeAnno[,'ID'], probeName),
                            c('ID', 'chr')]
    if ( nrow(probeAnno) == 0 ){ 
        stop("Wrong 'probeAnno' provided!")
    }    
    chrDist <- table(probeAnno[,"chr"])
    currentChr <- names(chrDist)
    message("Chromosomes: ")  
    print(currentChr) 

    matlist <- list()
    for (i in currentChr){  
        annoSub <- probeAnno[probeAnno[,"chr"] == i,]
        entrezIDperChr <- annoSub[,"ID"] 
        dataXsub <- dataX[,intersect(probeName, entrezIDperChr)] 
        mat <- cbind(label=dataY, dataXsub) 
        matlist[[i]] <- mat 
    } 

    names(matlist) <- paste0("chr", currentChr)
    ## exclude the label
    probeNperChr <- unlist(lapply(matlist, function(i) {ncol(i)-1}))  
    message("Summary stat of # probes in each Chr: ")
    print(summary(probeNperChr))  
    return(matlist) 
}