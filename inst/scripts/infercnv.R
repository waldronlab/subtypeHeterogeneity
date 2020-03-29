suppressPackageStartupMessages({
    library(infercnv)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(SingleCellExperiment)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "cnv.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))


## 1. Writing input files
cat("Writing input files ...:\n", file=log.file)

gorder.file <- file.path(sample.dir, paste0("sample", snr, "_gene_order.tab"))
eids <- AnnotationDbi::mapIds(org.Hs.eg.db, column = "ENTREZID", keys = rownames(sce), keytype = "SYMBOL")

ind <- !is.na(eids)
sce <- sce[ind,]
eids <- eids[ind]
ind <- unname(eids) %in% names(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
sce <- sce[ind,]
eids <- eids[ind]
egenes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)[unname(eids)]
names(egenes) <- names(eids)
egenes <- sort(egenes)
egenes <- keepStandardChromosomes(egenes, pruning.mode = "coarse")
sce <- sce[names(egenes), ]

egenes <- as.data.frame(egenes)
write.table(egenes[,1:3], file = gorder.file, sep = "\t", 
            quote = FALSE, col.names = FALSE)

cmat.file <- file.path(sample.dir, paste0("sample", snr, "_count_matrix.tab"))
cmat <- as.matrix(assay(sce))
colnames(cmat) <- sce$Barcode
write.table(cmat, file = cmat.file, sep = "\t", quote = FALSE)

sanno.file <- file.path(sample.dir, paste0("sample", snr, "_sample_anno.tab"))
sce <- sce[,!is.na(sce$celltype)]
sanno <- colData(sce)[,c("Barcode", "celltype")]
#sanno[,2] <- ifelse(sanno[,2] == "EPI", "malignant", "normal")
write.table(sanno, file = sanno.file, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)


## 2. Creating infercnv object
cat("Creating infercnv object ...:\n", file=log.file)
ref.names <- unique(sce$celltype)
ref.names <- ref.names[ref.names != "EPI"]
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = cmat.file,
                                    annotations_file = sanno.file,
                                    delim = "\t",
                                    gene_order_file = gorder.file,
                                    ref_group_names = ref.names) 
#                                    ref_group_names = "normal") 

## 3. Run infercnv (sample mode)
cat("Run infercnv (sample mode) ...:\n", file=log.file)
res <- infercnv::run(infercnv_obj,
                             cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir = sample.dir, 
                             num_threads = 20,
                             #cluster_by_groups = TRUE, 
                             denoise = TRUE,
                             HMM = TRUE)

res.file <- file.path(sample.dir, paste0("sample", snr, "_infercnv_sample_res.rds"))
saveRDS(res, res.file)

## 4. Run infercnv (subcluster mode)
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
