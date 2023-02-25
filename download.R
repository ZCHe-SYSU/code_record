#BIocMAnager::install("remote")
#BiocManager::install("ExperimentHub")
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
#download from TCGA
a <- "LUSC"
a <- commandArgs(T)
library(TCGAbiolinks)
library(stringr)
library(SummarizedExperiment)
library(data.table)
factor2 <- paste('TCGA-',a,sep='')
query.counts <- GDCquery(project = factor2 ,
	data.category = 'Transcriptome Profiling' ,
	data.type='Gene Expression Quantification' ,
	workflow.type ='STAR - Counts')
#GDCdownload(query.counts)
prepare.counts <-GDCprepare(query = query.counts)
count <- data.table(assay(prepare.counts),keep.rownames="ID")
fwrite(count,"counts.csv")
query <- GDCquery(project = factor2 ,data.category = 'Clinical' ,file.type = 'xml')
GDCdownload(query)

