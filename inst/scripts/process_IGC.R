############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-08-15 17:11:30
# 
# descr: processing data on chemoresistant OV 
#        (from Patch et al, Nature, 2015)  
#
# data downloaded from: https://dcc.icgc.org 
# 
############################################################

# RNA-seq data
dat <- read.delim("exp_seq.tsv", as.is=TRUE)

extractMatrix <- function(dat, what=c("raw", "normalized"))
{
    what <- match.arg(what)
    col <- paste(what, "read_count", sep="_")
    mat <- split(dat[,col], dat[,"submitted_sample_id"])
    mat <- do.call(cbind, mat)
    colnames(mat) <- unique(dat[,"submitted_sample_id"])
    rownames(mat) <- unique(dat[,"gene_id"])
    return(mat)
}

raw.mat <- extractMatrix(dat, "raw")
norm.mat <- extractMatrix(dat, "normalized")

# clinical data

# code: 
# 1 Primary tumour - solid tissue
# 2 Recurrent tumour - solid tissue
# 4 Recurrent tumour - solid tissue
# 9 Normal - blood derived
# 12 Primary tumour - other
# 13 Recurrent tumour - other
# 14 Recurrent tumour - other
# 16 Metastatic tumour - metastasis to distant location
# 17 Metastatic tumour - metastasis to distant location

# sample data
sampl.dat <- read.delim("sample.tsv", as.is=TRUE)
rownames(sampl.dat) <- sampl.dat[,"submitted_specimen_id"] 

# donor data
donor.dat <- read.delim("donor.tsv", as.is=TRUE)
rownames(donor.dat) <- donor.dat[,"submitted_donor_id"] 

sample2donor <- sampl.dat[,"submitted_donor_id"]
names(sample2donor) <- sampl.dat[,"submitted_sample_id"]

don <- sample2donor[colnames(raw.mat)]  
don.cols <- c("disease_status_last_followup", "donor_age_at_diagnosis",
                "donor_relapse_interval",  "donor_tumour_stage_at_diagnosis", 
                    "donor_survival_time")
don.cols <- donor.dat[don,don.cols]

# specimen data
spec.dat <- read.delim("specimen.tsv", as.is=TRUE)
rownames(spec.dat) <- spec.dat[,"submitted_specimen_id"] 

sample2spec <- sampl.dat[,"submitted_specimen_id"]
names(sample2spec) <- sampl.dat[,"submitted_sample_id"]

spec <- sample2spec[colnames(raw.mat)]  
spec.cols <- c("specimen_type", "specimen_type_other", 
                "specimen_donor_treatment_type", "specimen_donor_treatment_type_other")
spec.cols <- spec.dat[spec,spec.cols]

# combined clinical data
cdat <- cbind(spec, don, spec.cols, don.cols)
colnames(cdat) <- sub("^specimen_", "", colnames(cdat))
colnames(cdat) <- sub("^donor_", "", colnames(cdat))

# create coordinated data object
library(SummarizedExperiment)
se <- SummarizedExperiment(assays = list(raw = raw.mat, norm=norm.mat),
                            colData = cdat)
# map IDs
library(EnrichmentBrowser)
se <- idMap(se, org="hsa", from="ENSEMBL", to="ENTREZID")

# consensus subtypes
library(consensusOV)
cst <- get.consensus.subtypes(log(assays(se)$norm + 0.1, base=2), names(se))

# alternative voom processing
dge <- DGEList(assays(se)$raw)
keep <- filterByExpr(dge)
dge <- dge[keep,]
dge <- calcNormFactors(dge)
voom.mat <- limma::voom(dge)$E
cst.voom <- get.consensus.subtypes(voom.mat, rownames(voom.mat))

# find pairings: primary <-> recurrent
se$spec <- as.vector(se$spec)
se$don <- as.vector(se$don)
scode <- vapply(se$spec, function(s) unlist(strsplit(s, "-"))[[3]], character(1))

table(se$don)[table(se$don) > 1]

# AOCS-093-2-2: Primary tumour
# AOCS-093-10-1: primary ascitic fluid
# AOCS-093-4-X:  recurrent ascitic fluid

# AOCS-135 has two recurrent ascites, compare first against second?
# AOCS-141 has two recurrent ascites, compare first against second?
# AOCS-150 has two recurrent ascites, compare first against second?
# AOCS-170 has primary tumor and primary ascite - compare?
# AOCS-171 has primary tumor and primary ascite - compare?


## compare primary tumor and primary ascite
####################
> cst$rf.probs[c("AOCS-171-2-4", "AOCS-171-4-1"),]
             IMR_consensus DIF_consensus PRO_consensus MES_consensus
AOCS-171-2-4         0.050         0.046         0.404         0.500
AOCS-171-4-1         0.082         0.746         0.150         0.022
#
> cst$rf.probs[c("AOCS-170-2-1", "AOCS-170-4-9"),]
             IMR_consensus DIF_consensus PRO_consensus MES_consensus
AOCS-170-2-1          0.01         0.002         0.274         0.714
AOCS-170-4-9          0.33         0.542         0.034         0.094
#
> cst$rf.probs[c("AOCS-093-2-2", "AOCS-093-10-1"),]
              IMR_consensus DIF_consensus PRO_consensus MES_consensus
AOCS-093-2-2           0.06         0.802         0.094         0.044
AOCS-093-10-1          0.32         0.334         0.226         0.120
######################

> cst.voom$rf.probs[c("AOCS-171-2-4", "AOCS-171-4-1"),]
             IMR_consensus DIF_consensus PRO_consensus MES_consensus
AOCS-171-2-4         0.042         0.018         0.472         0.468
AOCS-171-4-1         0.108         0.636         0.226         0.030
> cst.voom$rf.probs[c("AOCS-170-2-1", "AOCS-170-4-9"),]
             IMR_consensus DIF_consensus PRO_consensus MES_consensus
AOCS-170-2-1         0.010         0.014         0.312         0.664
AOCS-170-4-9         0.346         0.534         0.054         0.066
> cst.voom$rf.probs[c("AOCS-093-2-2", "AOCS-093-10-1"),]
              IMR_consensus DIF_consensus PRO_consensus MES_consensus
AOCS-093-2-2          0.072         0.778         0.086         0.064
AOCS-093-10-1         0.296         0.400         0.156         0.148
######################

primary.ids <- c("AOCS-034-2-4",
                    "AOCS-064-2-X",
                    "AOCS-065-2-2",
                    "AOCS-086-2-9",
                    "AOCS-088-2-4",
                    "AOCS-090-2-4",
                    "AOCS-091-2-7",
                    "AOCS-092-2-X",
                    "AOCS-093-2-2",
                    "AOCS-093-10-1",
                    "AOCS-094-2-5",
                    "AOCS-095-2-8",
                    "AOCS-137-4-0") 
                
recurrent.ids <- c("AOCS-034-4-1",
                    "AOCS-064-4-7",
                    "AOCS-065-4-X",
                    "AOCS-086-4-6",
                    "AOCS-088-4-1",
                    "AOCS-090-4-1",
                    "AOCS-091-4-4",
                    "AOCS-092-4-7",
                    "AOCS-093-4-X",    
                    "AOCS-093-4-X",
                    "AOCS-094-4-2",
                    "AOCS-095-4-5",
                    "AOCS-137-10-2")

ptumor.ids <- c("AOCS-171-2-4", 
                "AOCS-170-2-1", 
                "AOCS-093-2-2")

pascite.ids <- c("AOCS-171-4-1",
                "AOCS-170-4-9",
                "AOCS-093-10-1")

# viz

cb.pink <- "#CC79A7"
cb.darkred <- "#B42F32"
cb.red <- "#D55E00"
cb.lightred <- "#DF6747"
cb.blue <- "#0072B2"
cb.yellow <- "#F0E442"
cb.green <- "#009E73"
cb.lightblue <- "#56B4E9"
cb.lightorange <- "#FAAC77"
cb.orange <- "#E69F00"
cb.darkorange <- "#F6893D"
cb.lightgrey <- "#C9C9BD"
cb.darkgrey <- "#878D92"

stcols <- c(cb.lightblue, cb.green, cb.orange, cb.pink) 
names(stcols) <- c("PRO", "MES", "DIF", "IMR")

# compare primary tumor and recurrent ascite
pr.ids <- list(primary.tumor = primary.ids, recurrent.ascite = recurrent.ids)

# compare primary tumor and primary ascite
ta.ids <- list(primary.tumor = ptumor.ids, primary.ascite = pascite.ids)


stHeatmap <- function(se, cst, ids, title = "")
{
    st <- sub("_consensus$", "", cst$consensusOV.subtypes)

    pind <- match(ids[[1]], colnames(se))
    prim <- cst$rf.probs[pind,]
    colnames(prim) <- sub("_consensus$", "", colnames(prim))
    rownames(prim) <- substring(rownames(prim), 1, 8)
    pst <- as.factor(st[pind])

    rind <- match(ids[[2]], colnames(se))
    recur <- cst$rf.probs[rind,]
    colnames(recur) <- sub("_consensus$", "", colnames(prim))
    rownames(recur) <- substring(rownames(recur), 1, 8) 
    rst <- as.factor(st[rind]) 

    scol <- list()
    scol[[names(ids)[1]]] <- scol[[names(ids)[2]]] <- stcols
   
    df <- data.frame(pst, rst)
    colnames(df) <- names(ids)

    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, 
        col=scol, 
        show_legend = c(TRUE, FALSE))

    names(ids)[2] <- "ascite" 
    Heatmap(rbind(t(prim), t(recur)), name = "ClassProb", cluster_rows=FALSE,
        cluster_columns=FALSE, top_annotation = ha, column_title=title,    
        row_split = factor(rep(names(ids), each=4), levels=names(ids)))
}   
 
h1 <- stHeatmap(se, cst.voom, pr.ids, title = "A) primary tumor / recurrent ascite")
h2 <- stHeatmap(se, cst.voom, ta.ids, title = "B) primary tumor / primary ascite")
ht_list <- h1 + h2
draw(ht_list, ht_gap = unit(4, "cm"))











