############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-11-08 14:54:00
# 
# descr: CN signatures of Macintyre et al. (Nat Genet, 2018) 
# 
############################################################

setwd("~/Documents/github/cnsignatures")

source("helper_functions.R")
source("main_functions.R")
setwd("manuscript_Rmarkdown")

library(Biobase)
library(flexmix)
library(NMF)
library(doMC)

CN_components <- readRDS("data/CN_components.rds")
sigs <- readRDS("data/sigs.rds")
tcga_sigs <- readRDS("data/tcga_sigs.rds")

reord_components<-c(11:13,24:31,17:23,32:36,14:16,1:10)

#britroc feat_sig matrix
nsig<-7
feat_sig_mat<-basis(sigs)
reord_britroc<-as.integer(c(2,6,5,4,7,3,1))
names(reord_britroc)<-paste0("s",1:7)
feat_sig_mat<-feat_sig_mat[,reord_britroc]
colnames(feat_sig_mat)<-paste0("s",1:nsig)
sig_feat_mat<-t(feat_sig_mat)

#tcga feat_sig matrix
feat_sig_mat_tcga<-basis(tcga_sigs)[,]
sig_feat_mat_tcga<-t(feat_sig_mat_tcga)
colnames(feat_sig_mat_tcga)<-paste0("s",1:nsig)

reord_tcga<-apply(feat_sig_mat,2,function(x){which.max(apply(feat_sig_mat_tcga,2,cor,x,method="s"))})
sigcor_tcga<-apply(feat_sig_mat,2,function(x){max(apply(feat_sig_mat_tcga,2,cor,x,method="s"))})

sig_feat_mat_tcga<-data.frame(sig_feat_mat_tcga[reord_tcga,],stringsAsFactors = F)
rownames(sig_feat_mat_tcga)<-paste0("s",1:nsig)
feat_sig_mat_tcga<-t(sig_feat_mat_tcga)


sig_thresh<-0.01

sig_pat_mat_tcga<-scoef(tcga_sigs)
rownames(sig_pat_mat_tcga)<-paste0("s",1:nsig)
sig_pat_mat_tcga<-sig_pat_mat_tcga[reord_tcga,]
rownames(sig_pat_mat_tcga)<-paste0("s",1:nsig)
sig_pat_mat_tcga<-normaliseMatrix(sig_pat_mat_tcga,sig_thresh)

saveRDS(sig_pat_mat_tcga, file="~/Documents/github/cnsignatures/manuscript_Rmarkdown/data/sig_pat_mat_tcga.rds")

# Get subtype annotations
# Get the TCGA OV microarray data:
library(curatedTCGAData)
ov.ctd <- curatedTCGAData(diseaseCode="OV", assays="mRNAArray_huex*", dry.run=FALSE)
se <- ov.ctd[[1]]
assays(se) <- list(as.matrix(assay(se)))
res.file <- file.path(data.dir, "consensus_subtypes_TCGA_marray_huex.rds")

# Assign consensus subtypes:
data.dir <- system.file("extdata", package="subtypeHeterogeneity") 
res.file <- file.path(data.dir, "consensus_subtypes_TCGA_marray_huex.rds")
cst <- readRDS(res.file)
sts <- as.vector(cst$consensusOV.subtypes)
sts <- sub("_consensus$", "", sts)
names(sts) <- substring(colnames(se), 1, 15)

colnames(sig_pat_mat_tcga) <- substring(colnames(sig_pat_mat_tcga), 1, 15)
isect <- intersect(colnames(sig_pat_mat_tcga), names(sts))
sts <- sts[isect]
sig_pat_mat_tcga <- sig_pat_mat_tcga[,isect]


st <- as.factor(sts)
df <- data.frame(Subtype = st)
scols <- list(Subtype=stcols)
ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col=scols)

ComplexHeatmap::Heatmap(sig_pat_mat_tcga, name="exposure", top_annotation = ha,
                            show_row_names=TRUE, show_column_names=FALSE,
                            column_title="Tumors", row_title="Signatures")


# background probabilities
bprobs <- table(sts) / sum(table(sts))
sts
      DIF       IMR       MES       PRO 
0.3063725 0.2647059 0.2450980 0.1838235 

# summary stats per sig
m <- apply(sig_pat_mat_tcga, 1, summary)
chisq.sig <- function(sig, stat = "Mean")
{
    ind <- sig_pat_mat_tcga[sig,] > m[stat,sig]
    chisq.test(table(sts[ind]), p = bprobs)$p.value
}
vapply(1:7, chisq.test, numeric(1))
[1] 0.3868094879 0.0751184679 0.0390187603 0.0201684821 0.5478753726
[6] 0.0001193329 0.0957830840

library(ggpubr)
SUBTYPES <- c("DIF", "IMR", "MES", "PRO")
bpType <- function(stat="Mean")
{   
    getPerc <- function(sig) 
    {
        ind <- sig_pat_mat_tcga[sig,] > m[stat,sig]
        tab <- table(sts[ind]) 
        tab / sum(tab) * 100
    }
    ttab <- vapply(1:7, getPerc, numeric(4))
    colnames(ttab) <- paste0("s", 1:7)
    ttab <- cbind(ttab, bprobs * 100)
    colnames(ttab)[8] <- "bg"
    df <- reshape2::melt(t(ttab))
    colnames(df) <- c("SIGNATURES", "SUBTYPE", "value")
    
    type <- "SUBTYPE"
    df[[type]] <- factor(df[[type]], levels = rev(SUBTYPES))
    pal <- rev(get_palette("npg", 4)[c(4,2,3,1)])

    ggbarplot(df, "SIGNATURES", "value", fill = type, legend = "right",
                color = type, palette = pal , ylab = "TUMORS [%]")
}
bpType()
df.stars <- data.frame(x=1:7, y=101, label=c("", "", "*", "*", "", "***", ""))
ggtext(df.stars, x, y, label)

