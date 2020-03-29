suppressPackageStartupMessages({
    library(DropletUtils)
    library(scater)
    library(scran)
    library(EnsDb.Hsapiens.v86)
})

snr <- commandArgs()[8]

## 2. Setting up the data
##### Reading in a sparse matrix
sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "preproc.log")
dir.create(plot.dir)

bp <- BiocParallel::registered()[[1]]

sce = read10xCounts(sample.dir)

cat("Initial dims:\n", file=log.file)
cat(dim(sce), file=log.file, append=TRUE)

##### Annotating the rows
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
    column="SEQNAME", keytype="GENEID")

##### Testing for deviations from ambient expression
bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)

pdf(file.path(plot.dir, "barcode_ranks.pdf"))
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
    xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=bcrank$inflection, col="darkgreen", lty=2)
abline(h=bcrank$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
    col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
dev.off()

## 4. Quality control on the cells
df <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")),
    BPPARAM = bp)

pdf(file.path(plot.dir, "cellQC.pdf"))
par(mfrow=c(1,3))

# hist(sce$log10_total_counts, breaks=20, col="grey80",
#     xlab="Log-total UMI count")
# hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
#     xlab="Log-total number of expressed features")
# hist(sce$pct_counts_Mito, breaks=20, col="grey80",
#     xlab="Proportion of reads in mitochondrial genes")

hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Mito, xlab="Mitochondrial proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")
dev.off()

high.mito <- isOutlier(df$subsets_Mito_percent, nmads = 3, type = "higher")
libsize.drop <- isOutlier(df$sum, nmads = 1, type = "lower", log = TRUE)
feature.drop <- isOutlier(df$detected, nmads = 1, type = "lower", log = TRUE)
sce <- sce[,!(high.mito | libsize.drop | feature.drop)]

df <- data.frame(ByHighMito=sum(high.mito),
           ByLibSize=sum(libsize.drop), 
           ByFeature=sum(feature.drop), 
           Remaining=ncol(sce))

cat("\n\n CellQC:\n", file=log.file, append=TRUE)
cat(paste("\n ByHighMito:", sum(high.mito)), file=log.file, append=TRUE)
cat(paste("\n ByLibSize:", sum(libsize.drop)), file=log.file, append=TRUE)
cat(paste("\n ByFeature:", sum(feature.drop)), file=log.file, append=TRUE)
cat(paste("\n Remaining:", ncol(sce)), file=log.file, append=TRUE)

cat("\n\n Dims after CellQC:\n", file=log.file, append=TRUE)
cat(dim(sce), file=log.file, append=TRUE)

## 5. Examining gene expression
ave <- calculateAverage(sce, BPPARAM = bp)
pdf(file.path(plot.dir, "average_expression_hist.pdf"))
hist(log10(ave), breaks=100, main="", col="grey",
    xlab=expression(Log[10]~"average count"))
dev.off()

#remove genes that have average counts of zero, as this means that they are not expressed in any cell
rowData(sce)$AveCount <- ave
to.keep = ave > 0.001
sce = sce[to.keep,]

cat(paste("\n\n Exluding", sum(to.keep), 
    "genes due to insufficient expression\n\n"), file=log.file, append=TRUE)

cat("Dims after low expression filter:\n", file=log.file, append=TRUE)
cat(dim(sce), file=log.file, append=TRUE)


## 6. Normalizing for cell-specific biases
clusters <- quickCluster(sce, method="igraph", min.mean=0.1, BPPARAM = bp)
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters, BPPARAM = bp)

pdf(file.path(plot.dir, "sizeFactors.pdf"))
plot(sce$total_counts, sizeFactors(sce), log="xy")
dev.off()

sce <- scater::normalize(sce)
saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
