suppressPackageStartupMessages({
    library(scater)
    library(SingleCellExperiment)
    library(monocle)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "monocle2.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

bp <- BiocParallel::registered()[[1]]

## 1. marker genes
cat("Running monocle ...:\n", file=log.file)

#!!! Monocle 2 only infers one trajectory for the entire dataset, so non-neuronal cells like endothelial cells and erythrocytes may be mistaken as highly differentiated cells from the neuronal lineage. So we will remove cell types not of the neural or glial lineages. Cell types are also helpful to orient the trajectory; neuronal progenitor cells must come before neurons. 
# restrict to epithelial cells

# try: full and filter mat
# try: only first two PCs
# try: clusterLabels and celltype labels

cds <- newCellDataSet(assays(sce)$counts, 
		      	phenoData = new("AnnotatedDataFrame", data=as.data.frame(colData(sce))),
			featureData = new("AnnotatedDataFrame", data=as.data.frame(rowData(sce))))

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

disp_table <- dispersionTable(cds)
clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, clustering_genes$gene_id)

cat("Reducing dimensions + clustering ...:\n", file=log.file)

cds <- reduceDimension(cds, num_dim = 40, reduction_method = 'tSNE')
cds <- clusterCells(cds, method = "louvain", append=TRUE)

#TODO: + subtype
diff_genes <- differentialGeneTest(cds, fullModelFormulaStr = "~ Cluster",
                                   cores = 20)

ordering_genes <- row.names(subset(diff_genes, qval < 1e-3))[order(diff_genes$qval)][1:3000]
cds <- setOrderingFilter(cds, ordering_genes)

#TODO: increase nr components
cat("Ordering cells ...:\n", file=log.file, append=TRUE)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

saveRDS(cds, file = file.path(sample.dir, paste0("sample", snr, "_monocle2.rds")))
# cds <- readRDS(file.path(sample.dir, paste0("sample", snr, "_monocle2.rds")))

# annotate to sce
sce$monocle2.pseudotime <- cds$Pseudotime 
sce$monocle2.state <- cds$State

pdf(file.path(plot.dir, "monocle2.pdf"))
plotTSNE(sce, colour_by="monocle2.pseudotime")
plotTSNE(sce, colour_by="monocle2.state")
plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1) + scale_color_viridis_c()
plot_cell_trajectory(cds, color_by = "State", cell_size = 1)
plot_cell_trajectory(cds, color_by = "subtype", cell_size = 1)
dev.off() 

saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
