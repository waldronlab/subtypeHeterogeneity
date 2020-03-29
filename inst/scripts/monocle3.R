suppressPackageStartupMessages({
    library(scater)
    library(SingleCellExperiment)
    library(monocle3)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "monocle3.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

bp <- BiocParallel::registered()[[1]]

## 1. marker genes
cat("Running monocle ...:\n", file=log.file)

# try: full and filter mat
# try: only first two PCs
# try: clusterLabels and celltype labels

cds <- new_cell_data_set(assays(sce)$counts, 
		      	cell_metadata = as.data.frame(colData(sce)),
			gene_metadata = as.data.frame(rowData(sce)))

cat("Reducing dimensions + clustering ...:\n", file=log.file)

cds <- preprocess_cds(cds, num_dim=50)

pdf(file.path(plot.dir, "monocle3.pdf"))
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "subtype")
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "hpca.celltype")

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cat("Learn trajectory graph ...:\n", file=log.file, append=TRUE)

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "subtype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE)

plot_cells(cds,
           color_cells_by = "hpca.celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE)
dev.off()
saveRDS(cds, file = file.path(sample.dir, paste0("sample", snr, "_monocle3.rds")))

# cds <- readRDS(file.path(sample.dir, paste0("sample", snr, "_monocle3.rds")))
# cds <- order_cells(cds)

# annotate to sce
#sce$monocle3.pseudotime <- cds$Pseudotime 
#sce$monocle3.state <- cds$State
#
#pdf(file.path(plot.dir, "monocle2.pdf"))
#plotTSNE(sce, colour_by="monocle.pseudotime")
#plotTSNE(sce, colour_by="monocle.state")
#plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1) + scale_color_viridis_c()
#plot_cell_trajectory(cds, color_by = "State", cell_size = 1)
#plot_cell_trajectory(cds, color_by = "subtype", cell_size = 1)
#dev.off() 
#
#saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
