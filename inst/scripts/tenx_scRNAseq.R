############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-08-22 16:02:45
# 
# descr: annotate cell types with SingleR
# 
############################################################

library(scater)

# install SingleR from https://github.com/LTLA/SingleR 
library(SingleR)

# exute the actual annotation on wallabe4, to leverage parallel computation
sce <- readRDS("sample59_sce_13k_cells.rda")

hpca.se <- HumanPrimaryCellAtlasData()
encode.se <- BlueprintEncodeData()

common <- intersect(rownames(sce), rownames(hpca.se))
hpca.se <- hpca.se[common,]
sce <- sce[common,]

pred.hpca <- SingleR(test = sce, ref = encode.se,
                        labels = hpca.se$label.main,
                        assay.type.ref = "logcounts",
                        BPPARAM = BiocParallel::registered()[[1]])

saveRDS(pred.hpca, file="sample59_celltypes_annotated_encode.rds")

# annotate and plot: can be done with laptop
hpc <- readRDS("sample59_celltypes_annotated_hpca.rds")
encode <- readRDS("sample59_celltypes_annotated_encode.rds")

colData(sce)$hpca.celltype <- hpc$labels 
colData(sce)$encode.celltype <- encode$labels

maxScore <- function(res) 
    vapply(1:nrow(res), function(i) res$scores[i, res[i,"labels"]], numeric(1))

colData(sce)$hpca.celltype.score <- maxScore(hpc) 
colData(sce)$encode.celltype.score <- maxScore(hpc) 

plotMainCellTypes <- function(sce, col)
{
    main.cats <- sort(table(sce[[col]]), decreasing=TRUE)
    main.cats <- names(main.cats)[main.cats > 99]

    sce.sub <- sce[,sce[[col]] %in% main.cats[1:6]]
    sce.sub[[col]] <- factor(as.vector(sce.sub[[col]]), levels=main.cats[1:6])
    plotTSNE(sce.sub, colour_by=col, shape_by=col)
}

# investigate ambiguity of cell type assignments
plotScoreHeatmap(hpc)
plotScoreHeatmap(encode)

