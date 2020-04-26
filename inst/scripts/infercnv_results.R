# Run with R-release, there are issues with ComplexHeatmap in devel 

#tumor 59: 3,500 epithelial cells vs. 9,500 reference (non-epithelial) cells
#tumor 61: 1,200 epithelial cells vs. 7,400 reference (non-epithelial) cells
#tumor 76: 210 epithelial cells vs. 13,500 reference (non-epithelial) cells
#tumor 77: 2,900 epithelial cells vs. 4,000 reference (non-epithelial) cells
#tumor 89: 300 epithelial cells vs. 4,600 reference (non-epithelial) cells
#tumor 90: 900 epithelial cells vs. 2,700 reference (non-epithelial) cells
#
library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(subtypeHeterogeneity)

SUBTYPES <- c("DIF", "IMR", "MES", "PRO")
CELL.TYPES <- c("EPI", "LYMPH", "MYE", "STROM", "ENDO")

cb.pink <- "#CC79A7"
cb.red <- "#D55E00"
cb.blue <- "#0072B2"
cb.yellow <- "#F0E442"
cb.green <- "#009E73"
cb.lightblue <- "#56B4E9"
cb.orange <- "#E69F00"

stcols <- c(cb.lightblue, cb.green, cb.orange, cb.pink)
names(stcols) <- SUBTYPES

sample.dir <- "/Users/ludwig/Documents/PostDoc/NewYork/ovc/scRNA_seq/"
#sample.dir <- "/nobackup/16tb_b/scRNA/10x_genomics/"
setwd(sample.dir)

# reading the sces
tumors <- c(59, 61, 76, 77, 89, "90_PM39") 
readSCE <- function(i) readRDS(gsub("XX", i, "tumorXX/sampleXX_sce.rds")) 
#readSCE <- function(i) readRDS(gsub("XX", i, "XX/sampleXX_sce.rds")) 
sces <- lapply(tumors, readSCE)

for(i in seq_along(sces)) colnames(sces[[i]]) <- sces[[i]]$Barcode

readObs <- function(i)
{
    obs.file <- file.path(paste0("tumor", i), "infercnv.observations.txt")
    #obs.file <- file.path(i, "infercnv.observations.txt")
    obs <- readr::read_delim(obs.file, delim = " ")
    obs.mat <- as.matrix(obs[,2:ncol(obs)])
    rownames(obs.mat) <- obs$GENE
    obs.mat
}
obs.mats <- lapply(tumors, readObs)

tumors <- c(59, 77, 90) 

for(i in seq_along(tumors))
{
    sce <- sces[[i]]
    obs.mat <- obs.mats[[i]]
    
    sce <- sce[,colnames(obs.mat)]
    ind <- sce$subtype %in% c("DIF", "PRO")
    obs.mat <- obs.mat[,ind]
    sce <- sce[,ind]

    sces[[i]] <- sce
    obs.mats[[i]] <- obs.mat
}

# intersection heatmap
isect <- Reduce(intersect, lapply(obs.mats, rownames))
obs.mat <- do.call(cbind, lapply(obs.mats, function(o) o[isect,]))

# row annotation: Subtype + number cells 
getCellAnnotation <- function(sces)
{
    sts <- do.call(c, lapply(sces, function(sce) sce$subtype))
    margins <- do.call(c, lapply(sces, function(sce) sce$margin))
    df <- data.frame(Subtype = sts, Margin = margins)

    margin.ramp <- circlize::colorRamp2(
                    seq(quantile(margins, 0.01), quantile(margins,
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))

    cols <- list(Subtype = stcols, Margin = margin.ramp)
    ComplexHeatmap::HeatmapAnnotation(df = df, col = cols, which = "row")
}
ra <- getCellAnnotation(sces)

# 2. column annotation: genomic ranges
getEGenes <- function(obs.mat)
{
    eids <- AnnotationDbi::mapIds(org.Hs.eg.db, column = "ENTREZID", 
                                keys = rownames(obs.mat), keytype = "SYMBOL")
    egenes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)[unname(eids)]
    egenes <- unstrand(egenes)
    names(egenes) <- names(eids)
    sort(egenes)
}
egenes <- getEGenes(obs.mat)
obs.mat <- obs.mat[names(egenes),]
chr <- as.character(seqnames(egenes))

#df <- data.frame(Chromosome = sts)
#ta <- ComplexHeatmap::HeatmapAnnotation(df = df, col = cols)

# 3. highlight specific recurrent regions with genes 
getCGenes <- function(obs.mat)
{
    cnv.genes <- getCnvGenesFromTCGA()
    eids <- AnnotationDbi::mapIds(org.Hs.eg.db, column = "ENTREZID", 
                                keys = names(cnv.genes), keytype = "SYMBOL")
    cgenes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)
    cgenes <- cgenes[unname(eids)] 
    cgenes <- unlist(GRangesList(lapply(cgenes, function(g) g[1])))
    names(cgenes) <- names(eids)
    clabs <- vector("character", length = nrow(obs.mat)) 
    ind <- rownames(obs.mat) %in% names(cgenes)
    clabs[ind] <- rownames(obs.mat)[ind]
    clabs
}
clabs <- getCGenes(obs.mat)

# 4. GISTIC
getGISTIC <- function(egenes, obs.mat)
{
    gisticOV <- gistic2RSE(ctype="OV", peak="wide")
    #gdf <- as.data.frame(rowRanges(gisticOV))
    #gstr <- paste(gdf[,1], paste(gdf[,2], gdf[,3], sep = "-"), sep = ":")
    #cat(gstr, file = "gistic_hg19.txt", sep = "\n")
    gtype <- rowData(gisticOV)$type  

    not.mapped <- c(22, 31, 43, 47, 67)
    gstr <- scan("gistic_hg38_liftover.txt", what = "character", sep = "\n")
    gstr <- gstr[-not.mapped]
    gistic <- GRanges(gstr)
    gtype <- gtype[-not.mapped]

    stopifnot(all(names(egenes) == rownames(obs.mat)))
    olaps <- findOverlaps(egenes, gistic)
    sh <- subjectHits(olaps)
    qh <- queryHits(olaps)
    ampdel <- rep("Neutral", nrow(obs.mat))
    ampdel[qh] <- gtype[sh]
    ampdel <- substring(ampdel, 1, 3) 
    df <- data.frame(GISTIC2 = ampdel)
    col <- list(GISTIC2 = c(Amp = "red", Del = "blue", Neu = "white"))
    ComplexHeatmap::HeatmapAnnotation(df = df, col = col)
}
ta <- getGISTIC(egenes, obs.mat)

ComplexHeatmap::Heatmap(t(obs.mat), 
    #top_annotation = ta,
    right_annotation = ra,
    name = "Expression", 
    col = colorRamp2(c(0.895, 1, 1.105), c("blue", "white", "red")),
    row_title_rot = 0, column_title_rot = 90,
    show_row_names = FALSE, show_column_names = TRUE, 
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = rep(paste0("T", tumors), sapply(obs.mats, ncol)),
    column_spl = factor(chr, levels = unique(chr)), 
    column_gap = unit(0, "mm"),
    column_title_gp = gpar(fontsize = 6),
    column_labels = clabs,
    column_names_gp = gpar(fontsize = 6),
    border = TRUE)

# union heatmap
un <- Reduce(union, lapply(obs.mats, rownames))
for(i in seq_along(obs.mats))
{
    om <- obs.mats[[i]]
    mgenes <- setdiff(un, rownames(om))
    m <- matrix(1, nrow = length(mgenes), ncol = ncol(om))
    rownames(m) <- mgenes
    colnames(m) <- colnames(om)
    obs.mats[[i]] <- rbind(om, m)
}
obs.mat <- do.call(cbind, lapply(obs.mats, function(o) o[un,]))

ra <- getCellAnnotation(sces)
egenes <- getEGenes(obs.mat)
clabs <- getCGenes(obs.mat)
ta <- getGISTIC(egenes, obs.mat)



refs <- readr::read_delim("infercnv.references.txt", delim = " ")
ref.mat <- as.matrix(refs[,2:ncol(refs)])
rownames(ref.mat) <- refs$GENE


