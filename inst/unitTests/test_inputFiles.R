test_inputFiles <- function() {

    library(RUnit)
    library(BioMM)
    ## Input data with the column "label"
    methylfile <- system.file("extdata", "methylData.rds", package="BioMM")  
    methylData <- readRDS(methylfile)   
    checkTrue(colnames(methylData)[1]=="label") 

    ## Annotation files for Mapping CpGs  
    featureAnnoFile <- system.file("extdata", "cpgAnno.rds", package="BioMM") 
    featureAnno <- readRDS(file=featureAnnoFile)   
    colnameAnno <- colnames(featureAnno) 
    checkTrue(any(grep("^ID", colnameAnno))) 
    checkTrue(any(grep("^entrezID", colnameAnno)))  

}
