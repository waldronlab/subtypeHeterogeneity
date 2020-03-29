suppressPackageStartupMessages({
    library(scater)
    library(SingleCellExperiment)
    library(SingleR)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

bp <- BiocParallel::registered()[[1]]

## annotate cell types

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
