suppressPackageStartupMessages({
    library(DropletUtils)
    library(scater)
    library(scran)
    library(pheatmap)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "cluster.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

bp <- BiocParallel::registered()[[1]]
# 7. Modelling the mean-variance trend
dec.pbmc <- modelGeneVarByPoisson(x=sce, BPPARAM=bp)
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1)

pdf(file.path(plot.dir, "sc_cluster.pdf"))
plot(dec.pbmc$mean, dec.pbmc$total, pch=16, cex=0.5,
    xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec.pbmc)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

## 8. Dimensionality reduction
#Denoise log-expression data by removing principal components corresponding to technical noise.

sce <- denoisePCA(sce, subset.row=top.pbmc, technical=dec.pbmc, BPPARAM = bp)

cat("Reduced dims:\n", file=log.file)
cat(ncol(reducedDim(sce, "PCA")), file=log.file, append=TRUE)

plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
    ylab="Proportion of variance explained")
abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")

sce <- runTSNE(sce, use_dimred="PCA", perplexity=30, BPPARAM = bp)
plotTSNE(sce, colour_by="sizeFactor")

# TODO: run UMAP

## 9. Clustering with graph-based methods
snn.gr <- buildSNNGraph(sce, use.dimred="PCA", BPPARAM = bp, k = 25)
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)

cat("\n\n Clusters:\n", file=log.file, append=TRUE)
cat(table(sce$Cluster), file=log.file, append=TRUE)

cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.weights=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)

pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
    color=colorRampPalette(c("white", "blue"))(100))

plotTSNE(sce, colour_by="Cluster")
dev.off()

saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
