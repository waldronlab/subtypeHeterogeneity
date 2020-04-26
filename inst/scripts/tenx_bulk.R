############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-10-30 17:27:47
# 
# descr: analysis of bulk RNA-seq for four 10X tumors 
# 
############################################################

library(EnrichmentBrowser)
library(consensusOV)

setwd("/Users/ludwig/Documents/PostDoc/NewYork/ovc/scRNA_seq")

dat <- read.delim("Bulk_RNAseq.txt", as.is=TRUE)
dat2 <- read.delim("Bulk_RNAseq2.txt", as.is=TRUE)

cpm.ind <- seq_len(which(dat[,1] == "")[1] - 1)
cpm.dat <- dat[cpm.ind, 1:8]
raw.dat <- dat[,9:16]

cpm.sym <- cpm.dat[,1]
raw.sym <- raw.dat[,1]

rownames(cpm.dat) <- cpm.dat[,2] 
cpm.dat <- cpm.dat[,3:8]
colnames(cpm.dat) <- paste0("T", c(59, 77, 76, 89, 90, 91)) 
cpm.dat <- as.matrix(cpm.dat)
mode(cpm.dat) <- "numeric"
cpm.se <- SummarizedExperiment(assays = list(cpm = cpm.dat))
cpm.se <- EnrichmentBrowser::idMap(cpm.se, org="hsa", from="ENSEMBL", to="ENTREZID")
cst <- consensusOV::get.consensus.subtypes(assay(cpm.se), names(cpm.se))

$consensusOV.subtypes
[1] IMR_consensus PRO_consensus MES_consensus IMR_consensus PRO_consensus
[6] DIF_consensus
Levels: IMR_consensus DIF_consensus PRO_consensus MES_consensus

$rf.probs
    IMR_consensus DIF_consensus PRO_consensus MES_consensus
T59         0.780         0.050         0.012         0.158
T77         0.288         0.254         0.332         0.126
T76         0.022         0.018         0.046         0.914
T89         0.838         0.062         0.012         0.088
T90         0.112         0.120         0.452         0.316
T91         0.030         0.946         0.018         0.006
attr(,"class")
[1] "matrix" "array"  "votes" 



rownames(raw.dat) <- raw.dat[,2] 
raw.dat <- raw.dat[,3:8]
colnames(raw.dat) <- paste0("T", c(59, 77, 76, 89, 90, 91)) 
raw.dat <- as.matrix(raw.dat)
mode(raw.dat) <- "numeric"
raw.dat <- edgeR::cpm(raw.dat, log=TRUE)
raw.se <- SummarizedExperiment(assays = list(raw = raw.dat))
raw.se <- EnrichmentBrowser::idMap(raw.se, org="hsa", from="ENSEMBL", to="ENTREZID")
cst2 <- consensusOV::get.consensus.subtypes(assay(raw.se), names(raw.se))

##
$consensusOV.subtypes
[1] IMR_consensus PRO_consensus MES_consensus IMR_consensus PRO_consensus
[6] DIF_consensus
Levels: IMR_consensus DIF_consensus PRO_consensus MES_consensus

$rf.probs
    IMR_consensus DIF_consensus PRO_consensus MES_consensus
T59         0.768         0.052         0.016         0.164
T77         0.298         0.266         0.338         0.098
T76         0.004         0.000         0.020         0.976
T89         0.836         0.078         0.008         0.078
T90         0.082         0.124         0.504         0.290
T91         0.016         0.962         0.016         0.006
attr(,"class")
[1] "matrix" "array"  "votes" 

