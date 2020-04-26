# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3440647/
# https://www.ncbi.nlm.nih.gov/pubmed/?term=Search+Results+Web+results++A+study+of+CCND1+with+epithelial+ovarian+cancer+cell
# https://en.wikipedia.org/wiki/G1_phase#In_cancer
# https://www.cyclacel.com/research_science_cell-cycle.shtml

suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(SingleCellExperiment)
    library(SingleR)
    library(scater)
    library(scran)
    library(scRNAseq)
    library(org.Hs.eg.db)
    library(EnrichmentBrowser)
    library(ggpubr)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/nobackup/16tb_b/scRNA/10x_genomics/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "ccycle.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

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

subsetEPI <- function(sce)
{
    cell.type <- "EPI" 
    subtype <- c("DIF", "PRO")
 
    sce <- sce[,!is.na(sce$celltype)]
    sce <- sce[,sce$celltype %in% cell.type]
    sce <- sce[,sce$subtype %in% subtype]
    ind <- do.call(order, colData(sce)[,c("subtype", "margin")])
    sce[,ind]
}
sce <- subsetEPI(sce) 

cat(paste0("Total EPI-DIF:", sum(sce$subtype == "DIF"), "\n"), 
    file = log.file)
cat(paste0("Total EPI-PRO:", sum(sce$subtype == "PRO"), "\n"), 
    file = log.file, append = TRUE)

## 1. Cyclins
cat("Cyclins ...:\n", file = log.file, append = TRUE)
cyclin.genes <- grep("^CCN[ABDE][0-9]$", names(sce), value = TRUE)

am <- as.matrix(assay(sce, "logcounts"))
am <- am[sort(cyclin.genes),]
ind <- colSums(am) > 0

rowMeans(am[,sce$subtype == "DIF"])
rowMeans(am[,sce$subtype == "PRO"])
boxplot(am["CCND1", tail(which(sce$subtype == "DIF"), 100)], am["CCND1", tail(which(sce$subtype == "PRO"), 100)])
t.test(am["CCND1", tail(which(sce$subtype == "DIF"), 100)], am["CCND1", tail(which(sce$subtype == "PRO"), 100)])

## ggboxplot of cyclin expression between DIF and PRO
df <- reshape2::melt(am)
df$subtype <- sce$subtype[df[,2]]
df$margin <- sce$margin[df[,2]]
ind.dif <- tail(which(sce$subtype == "DIF"), 100)
ind.pro <- tail(which(sce$subtype == "PRO"), 100)
df$top <- ifelse(df[,2] %in% c(ind.dif, ind.pro), "Top100", "All")
df2 <- df[df$top == "Top100",]
df2$top <- "All"
df <- rbind(df, df2) 
colnames(df)[c(1,3)] <- c("Cyclins", "log2CPM")
all.lab <- paste0("All epithelial cells (", 
                    sum(sce$subtype == "DIF"), " DIF, ", 
                    sum(sce$subtype == "PRO"), " PRO", ")")
top.lab <- paste0("Top 100 highest margin cells (", 
                    length(ind.dif), " DIF, ", 
                    length(ind.pro), " PRO", ")")

ggboxplot(df, "Cyclins", "log2CPM", color = "subtype", 
            palette = stcols[subtype], facet.by = "top",
            nrow = 2, ncol = 1, 
            panel.labs = list(top = c(all.lab, top.lab)))

cat(paste0("Cyclin EPI-DIF:", sum(sce$subtype[ind] == "DIF"), "\n"), 
    file = log.file, append = TRUE)
cat(paste0("Cyclin EPI-PRO:", sum(sce$subtype[ind] == "PRO"), "\n"), 
    file = log.file, append = TRUE)


## cyclin heatmap
sts <- colData(sce)[ind, "subtype"]
margins <- colData(sce)[ind, "margin"]
df <- data.frame(Subtype = sts, Margin = margins)

margin.ramp <- circlize::colorRamp2(
                    seq(quantile(margins, 0.01), quantile(margins, 
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))

ha <- ComplexHeatmap::HeatmapAnnotation(df = df, 
                                        col = list(Subtype = stcols, Margin = margin.ramp))

pdf(file.path(sample.dir, "cyclins.pdf"))                                
draw(Heatmap(am[,ind], name = "log2TPM", top_annotation = ha, 
        cluster_rows = FALSE, cluster_columns = FALSE))
dev.off()

## Marker heatmap
# markers (most DE) between PRO and DIF cells
markers <- findMarkers(sce, groups = sce$subtype)

# TODO: split heatmap into top 10 up in DIF and top 10 up in PRO
# TODO: check whether these markers correspond to Verhaak genes (annotation bar)

# markers for T59
> head(markers$DIF)
DataFrame with 6 rows and 4 columns
               Top               p.value                   FDR
         <integer>             <numeric>             <numeric>
PDZK1IP1         1 7.67750820550576e-280 1.38801670847338e-275
LCN2             2 6.04990167538199e-183 5.46880861946149e-179
RARRES1          3 5.17450949826119e-128 3.11833190730212e-124
NNMT             4  2.32691294228646e-92  1.05170647708991e-88
SLPI             5  6.69699645399425e-72  2.42149997783525e-68
CLU              6  2.45787462560593e-63  7.40598589272154e-60
                 logFC.PRO
                 <numeric>
PDZK1IP1  1.04129135790315
LCN2       1.4421170032941
RARRES1   1.30003069747262
NNMT     0.966470113525509
SLPI      1.04623530516874


Lipocalin2 Expressions Correlate Significantly With Tumor Differentiation in Epithelial Ovarian Cancer

> markers$DIF["CCND1",]
DataFrame with 1 row and 4 columns
            Top              p.value                  FDR          logFC.PRO
      <integer>            <numeric>            <numeric>          <numeric>
CCND1        69 1.30891435385047e-17 3.42954530482067e-15 -0.411589958912849


> head(same.markers$DIF[same.markers$DIF$logFC.PRO < 0,])
DataFrame with 6 rows and 4 columns
               Top              p.value                  FDR          logFC.PRO
         <integer>            <numeric>            <numeric>          <numeric>
STMN1           21 5.90138550762085e-36 5.08053088534654e-33 -0.483397822380303
MARCKSL1        31 8.51181647903995e-29 4.96403645563109e-26 -0.444470919223906
MFAP2           38 1.28208849992324e-25 6.09970473424007e-23 -0.375557828173884
BMP7            41 2.67123703429496e-24  1.1778852278785e-21 -0.305717911076186
ODC1            53 3.48572777936481e-21  1.1890277834554e-18 -0.388359574494035
SOX4            56 1.17823515216135e-19 3.80380594927233e-17 -0.457772730094846

Stathmin gene silencing suppresses proliferation, migration and invasion of gastric cancer cells via AKT/sCLU and STAT3 signaling. Shu F, et al. Int J Oncol, 2019 Mar. PMID 30628664
MiR-101 inhibits cell proliferation and invasion of pancreatic cancer through targeting 
STMN1. Zhu L, et al. Cancer Biomark, 2018. PMID 30198871 
MARCKSL1 promotes the proliferation, migration and invasion of lung adenocarcinoma cells


Heatmap(am[rownames(same.markers$DIF)[1:25],ind], name = "log2TPM", top_annotation = ha, cluster_rows = FALSE, cluster_columns = FALSE)

# DE
# wilcox-test
markers <- findMarkers(sce, subset.row=cyclin.genes, 
                        group=sce$subtype, test.type="wilcox", direction="up")
for(n in names(markers))
{
    cat(paste(n, "\n"), file = log.file, append = TRUE) 
    write.table(markers[[n]], file = log.file,
                sep = "\t", quote = FALSE, append = TRUE)
}

# t-test
markers <- findMarkers(sce, subset.row=cyclin.genes, group=sce$subtype)
for(n in names(markers))
{ 
    cat(paste(n, "\n"), file = log.file, append = TRUE) 
    write.table(markers[[n]], file = log.file,
                sep = "\t", quote = FALSE, append = TRUE)
}

## 2. SingleR
cat("SingleR ...:\n", file = log.file, append = TRUE)
sce.ref <- LengESCData()
names(assays(sce.ref)) <- "counts"
sce.ref <- logNormCounts(sce.ref)

# get genes annotated to cell cycle
cycle.anno <- select(org.Hs.eg.db, keytype="GOALL", 
                        keys="GO:0007049", columns="SYMBOL")[,"SYMBOL"]

# Find genes that are present in both datasets and are cell cycle-related
candidates <- Reduce(intersect, 
    list(rownames(sce.ref), rownames(sce), cycle.anno))

# Identifying markers between cell cycle phases
phase.stats <- pairwiseWilcox(logcounts(sce.ref), sce.ref$Phase, 
    direction="up", subset.row=candidates)

cycle.markers <- getTopMarkers(phase.stats[[1]], phase.stats[[2]])

sce.ref <- sce.ref[,!is.na(sce.ref$Phase)]
sr.assignments <- SingleR(sce, ref=sce.ref, label=sce.ref$Phase, genes=cycle.markers)
saveRDS(sr.assignments, file = file.path(sample.dir, "SingleR.rds"))
tab <- table(sr.assignments$labels, sce$subtype)
cat("SingleR ccycle:\n", file = log.file, append = TRUE) 
write.table(tab, file = log.file, sep = "\t", quote = FALSE, append = TRUE)
cat(paste("\n", chisq.test(tab)$p.value, "\n"), file = log.file, append = TRUE)

## 3. Cyclone
cat("Cyclone ...:\n", file = log.file, append = TRUE)

sce.ensembl <- idMap(sce, org = "hsa", from = "SYMBOL", to = "ENSEMBL")
hm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
cy.assignments <- cyclone(sce.ensembl, hm.pairs)
saveRDS(cy.assignments, file = file.path(sample.dir, "cyclone.rds"))

pdf(file.path(sample.dir, "cyclone.pdf"))
plot(cy.assignments$score$G1, cy.assignments$score$G2M,
    xlab="G1 score", ylab="G2/M score", pch=16)
tab <- table(cy.assignments$phases, sce$subtype)
dev.off()

write.table(tab, file = log.file, sep = "\t", quote = FALSE, append = TRUE)
cat(paste("\n", chisq.test(tab)$p.value, "\n"), file = log.file, append = TRUE)

#saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
