suppressPackageStartupMessages({
    library(scater)
    library(SingleCellExperiment)
    library(consensusOV)
    library(EnrichmentBrowser)
    library(SingleR)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "annotate.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

bp <- BiocParallel::registered()[[1]]

## 1. marker genes
cat("Plotting marker genes ...:\n", file=log.file)

pdf(file.path(plot.dir, "marker_genes.pdf"))

# Adopted from Shih et al., PLoS One (2018)
epithelial = c('EPCAM', 'KRT8', 'KRT18', 'KRT19')
lymphocyte = c('PTPRC', 'CD3E', 'CD19', 'MS4A1')
endothelial = c('PECAM1', 'CD34')
fibroblast = c('ACTA2', 'DCN', 'ACTB')
stromal = c('THY1', 'ENG', 'VIM', 'CD44')

### Ovarian cancer marker expression
ov.epithelial = c('WFDC2', 'CLDN4', 'FXYD3', 'CD24', 'ELF3', 'CLDN3', 
                        'MUC1', 'SPINT2', 'KRT8', 'SLPI', 'KRT18', 'KRT19')

### Cell origin marker expression
origin = c("PAX8", "CALB2", "COL11A1")

fontsize <- theme(axis.text=element_text(size=8), axis.title=element_text(size=10))

multiplot.markers <- function(markers, type)
{
    mlist <- vector("list", length=length(markers))
    for(i in seq_along(markers)) 
        mlist[[i]] <- plotTSNE(sce, colour_by = markers[i]) + fontsize + ggtitle(type)
    multiplot(plot.list = mlist, cols = 2)
}

### Epithelial marker expression
multiplot.markers(epithelial, "epithelial")
multiplot.markers(lymphocyte, "lymphocyte")
multiplot.markers(endothelial, "endothelial")
multiplot.markers(fibroblast, "fibroblast")
multiplot.markers(stromal, "stromal")
multiplot.markers(ov.epithelial, "ov.epithelial")
multiplot.markers(origin, "origin")
dev.off()

## 2. consensus subtypes
cat("Computing subtypes ...\n", file=log.file, append=TRUE)
margin <- function (rf.probs) 
{
    subtr <- apply(rf.probs, 1, function(row) sort(row)[3])
    pred.margins <- Biobase::rowMax(rf.probs) - subtr
    return(pred.margins)
}


sce.entrez <- idMap(sce, org="hsa", from="SYMBOL", to="ENTREZID")
cst <- get.consensus.subtypes(as.matrix(assays(sce.entrez)$logcounts), names(sce.entrez))

sts <- sub("_consensus$", "", as.vector(cst$consensusOV.subtypes))
sce$subtype <- sts
sce$margin <- margin(cst$rf.probs)  

s <- summary(sce$margin)

cat("Margins:\n", file=log.file, append=TRUE)
cat(names(s), file=log.file, append=TRUE)
cat("\n", file=log.file, append=TRUE)
cat(s, file=log.file, append=TRUE)


pdf(file.path(plot.dir, "subtypes.pdf"))
plotTSNE(sce, colour_by = "subtype")
plotTSNE(sce, colour_by = "margin")
dev.off()


## 3. annotate cell types

cat("\n\n Annotate cell types ... \n", file=log.file, append=TRUE)

annotateCellType <- function(sce, ref=c("hpca", "encode"))
{ 
    ref <- match.arg(ref)
    if(ref == "hpca") se <- HumanPrimaryCellAtlasData()
    else se <- BlueprintEncodeData()
    
    common <- intersect(rownames(sce), rownames(se))
    se <- se[common,]
    sce <- sce[common,]

    pred <- SingleR(test = sce, ref = se,
                        labels = se$label.main,
                        assay.type.ref = "logcounts",
                        BPPARAM = BiocParallel::registered()[[1]])

    fname <- paste0("sample", snr, "_celltypes_", ref, ".rds")
    saveRDS(pred, file = file.path(sample.dir, fname))
    
    return(pred)
}

# annotate and plot: can be done with laptop
hpc <- annotateCellType(sce, "hpca")
encode <- annotateCellType(sce, "encode")

colData(sce)$hpca.celltype <- hpc$labels 
colData(sce)$encode.celltype <- encode$labels

maxScore <- function(res)
    vapply(1:nrow(res), function(i) res$scores[i, res[i,"labels"]], numeric(1))

colData(sce)$hpca.celltype.score <- maxScore(hpc)
colData(sce)$encode.celltype.score <- maxScore(hpc)

plotMainCellTypes <- function(sce, col)
{
    main.cats <- sort(table(sce[[col]]), decreasing=TRUE)
    main.cats <- names(main.cats)

    sce.sub <- sce[,sce[[col]] %in% main.cats[1:6]]
    sce.sub[[col]] <- factor(as.vector(sce.sub[[col]]), levels=main.cats[1:6])
    plotTSNE(sce.sub, colour_by=col, shape_by=col)
}

pdf(file.path(plot.dir, "cell_types.pdf"))
plotMainCellTypes(sce, "hpca.celltype")
plotMainCellTypes(sce, "encode.celltype")
# investigate ambiguity of cell type assignments
plotScoreHeatmap(hpc)
plotScoreHeatmap(encode)
dev.off()

saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
