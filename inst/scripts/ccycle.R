suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(SingleCellExperiment)
    library(SingleR)
    library(scater)
    library(scran)
    library(scRNAseq)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "ccycle.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))


## 1. Cyclins
cat("Cyclins ...:\n", file = log.file)
cyclin.genes <- grep("^CCN[ABDE][0-9]$", names(sce), value = TRUE)

am <- as.matrix(assay(sce, "logcounts"))
am <- am[sort(cyclin.genes),]
ind <- colSums(am) > 0 & !is.na(sce$celltype) & sce$celltype == "EPI"

df <- data.frame(Subtype = colData(sce)[ind, "subtype"])
ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col = stcols)
                                
Heatmap(am[,ind], top_annotation = ha, cluster_rows = FALSE)


## 2. SingleR
cat("SingleR ...:\n", file = log.file, append = TRUE)
leng <- LengESCData()

## 3. Cyclone
cat("Cyclone ...:\n", file = log.file, append = TRUE)

## 4. 
cat("Run infercnv (sample mode) ...:\n", file = log.file)
scl.dir <- file.path(sample.dir, "subcluster")
res2 <- infercnv::run(infercnv_obj,
                                cutoff = 0.1,
                                out_dir = scl.dir, 
                                num_threads = 20,
                                #cluster_by_groups = TRUE, 
                                analysis_mode = "subclusters",
                                denoise = TRUE,
                                HMM = TRUE)

res.file <- file.path(sample.dir, paste0("sample", snr, "_infercnv_subcluster_res.rds"))
saveRDS(res, res.file)

#saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
