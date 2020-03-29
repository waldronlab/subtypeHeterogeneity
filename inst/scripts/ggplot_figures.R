############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-10-17 17:01:07
# 
# descr: 
# 
############################################################

library(ggpubr)
library(scater)
library(viridis)

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

# facet subtype
tumors <- c(59, 61, 76, 77, 89) 
readSCE <- function(i) readRDS(gsub("XX", i, "tumorXX/sampleXX_sce.rds")) 
sces <- lapply(tumors, readSCE)

facetTumors <- function(sces, col = "subtype", pal = stcols, cont = FALSE)
{
    psces <- lapply(sces, plotTSNE, colour_by = col)
    psces <- lapply(psces, function(p) p$data)
    for(i in tumors) psces[[match(i, tumors)]]$tumor <-  paste0("Tumor", i)
    df <- do.call(rbind, psces)

    colnames(df)[3] <- col
    df <- df[!is.na(df[[col]]),]
    p <- ggplot(df, aes_string(x = "X", y = "Y", fill = col)) + 
            geom_point(alpha = 0.75, shape=21, color="grey", size=1) + 
        facet_wrap(~tumor, nrow=1, ncol=5) + 
        xlab("Dimension 1") +
        ylab("Dimension 2") +
        labs(fill = col) +
        theme_bw(base_size = 12) #+
        #guides(fill = guide_legend(override.aes = list(size=7))) 
    
    is.cont <- col %in% names(sces[[1]]) || cont
    if(is.cont) p <- p + scale_fill_viridis() 
    else p <- p + scale_fill_manual(values=pal)
    p
}

stp <- facetTumors(sces, col = "subtype")
pal <- get_palette("npg", 5)
names(pal) <- CELL.TYPES
ctp <- facetTumors(sces, col = "celltype", pal = pal)
ggarrange(stp, ctp, nrow=2, align = "hv")

# summarize
annoCellType <- function(sce)
{
    hpca <- sce$hpca.celltype
    encode <- sce$encode.celltype

    .decide <- function(hp, en)
    {
        epithelial <- c("Epithelial_cells", "Epithelial cells", "Embryonic_stem_cells",
                        "MSC", "iPS_cells", "Keratinocytes", "Neuroepithelial_cell",
                        "Neurons")
        myeloid <- c("Macrophage", "Monocyte", "Macrophages", "Monocytes", "DC",
                        "Erythrocytes", "Neutrophils")
        lymphocyte <- c("T_cells", "B_cell", "Pre-B_cell_CD34-", "NK_cell", "CD8+ T-cells",
                        "CD4+ T-cells", "HSC", "B−cells", "NK cells", "Pro-B_cell_CD34+")
        endothelial <- c("Endothelial_cells", "Endothelial cells")
        stromal <- c("Fibroblasts", "Tissue_stem_cells", "Smooth_muscle_cells",
                        "Osteoblasts", "Chondrocytes", "Myocytes", "Fibroblasts", 
                        "Pericytes", "Chondrocytes", "Smooth muscle", "Adipocytes")
        
        if(hp %in% epithelial  || en %in% epithelial) ct <- "EPI"
        else if(hp %in% myeloid || en %in% myeloid) ct <- "MYE"
        else if(hp %in% lymphocyte || en %in% lymphocyte) ct <- "LYMPH"
        else if(hp %in% endothelial || en %in% endothelial) ct <- "ENDO"
        else if(hp %in% stromal || en %in% stromal) ct <- "STROM"
        else ct <- NA_character_
        return(ct)
    }
    
    sce$celltype <- vapply(seq_len(ncol(sce)), 
                            function(i) .decide(hpca[i], encode[i]), 
                            character(1)) 
    return(sce)
}

sces <- lapply(sces, annoCellType)

##
### arrange subtype and celltype side-by-side
##
plotType <- function(sce, type=c("celltype", "subtype"))
{
    pal <- "lancet"
    
    type <- match.arg(type)
    p <- plotTSNE(sce, colour_by = type)
    df <- p$data
    type <- toupper(type)
    colnames(df)[3] <- type 
    df <- df[!is.na(df[[type]]),]
    
    if(type == "CELLTYPE") df[[type]] <- factor(df[[type]], levels = CELL.TYPES)
    else pal <- get_palette("npg", 4)[c(4,2,3,1)]
    
    p <- ggscatter(df, x = "X", y = "Y", shape = 21, size = 1.5, color = "grey", 
                    alpha = 0.8, fill = type, palette = pal, ggtheme = theme_bw())
    ggpar(p, xlab = "Dimension 1", ylab = "Dimension 2")
}

p2a <- plotType(sces[[1]])
p2b <- plotType(sces[[1]], type="subtype")

p2 <- ggarrange(p2b + theme_bw() + theme(plot.margin = margin(r = 1)), 
          p2a + theme_bw() + theme(axis.text.y = element_blank(), 
                         axis.ticks.y = element_blank(), 
                         axis.title.y = element_blank(),
                         plot.margin = margin(l = 1, r=1)), 
            legend="top", align="h", labels = c("A", "B"))

# or
p3 <- ggarrange(p2b + theme_bw()+ theme(axis.text.x = element_blank(),  
                         axis.ticks.x = element_blank(), 
                         axis.title.x = element_blank()),  
            p2a + theme_bw(),   
            nrow=2, legend="none", align="hv", labels = c("A", "B"))

##
### barplot percentage of subtype / cell type for all 5 tenx tumors
##
getPerc <- function(sce, col="celltype", perc=TRUE)
{
    col <- sce[[col]]
    col <- col[!is.na(col)]
    ctab <- table(col)
    if(perc) ctab <- round(ctab / sum(ctab) * 100) 
    return(ctab)
} 

cttab <- vapply(sces, getPerc, numeric(5), perc=FALSE)
total <- colSums(cttab)
prop.test(cttab["EPI",], total)
prop.test(cttab["STROM",], total)

bpType <- function(sces, type=c("celltype", "subtype"))
{
    type <- match.arg(type)
    n <- ifelse(type == "celltype", 5, 4)
    cttab <- vapply(sces, getPerc, numeric(n), col = type)
    colnames(cttab) <- paste0("T", tumors)
    df <- reshape2::melt(t(cttab))
    type <- toupper(type)
    colnames(df) <- c("TUMORS", type, "value")
    
    if(type == "CELLTYPE") 
    {
        df[[type]] <- factor(df[[type]], levels = rev(CELL.TYPES))
        pal <- rev(get_palette("lancet", 5))
    }
    else 
    {
        df[[type]] <- factor(df[[type]], levels = rev(SUBTYPES))
        pal <- rev(get_palette("npg", 4)[c(4,2,3,1)])
    }

    ggbarplot(df, "TUMORS", "value", fill = type, 
                color = type, palette = pal , ylab = "CELLS [%]")
}

stbp <- bpType(sces, "subtype")
ctbp <- bpType(sces) 

p1 <- ggarrange(stbp + theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank()),
                        # axis.ticks.x = element_blank(),
                        # axis.line.x = element_blank()), 
          ctbp,
          ncol = 1,
          legend="left", align="hv", labels = c("C", "D"))

##
### CT vs ST matrix
##
cts <- unlist(lapply(sces, function(sce) sce$celltype))
sts <- unlist(lapply(sces, function(sce) sce$subtype))
nna.ind <- !is.na(cts)
cts <- cts[nna.ind]
sts <- sts[nna.ind]
spl <- split(sts, cts)
mat <- vapply(spl, table, numeric(4))

# barplot instead
perc <- round(mat / rowSums(mat) * 100, digits=1)
df <- reshape2::melt(perc)
type <- "CELL.TYPE"; colnames(df) <- c("SUBTYPE", type, "value")
df[[type]] <- factor(df[[type]], levels = rev(CELL.TYPES))
pal <- rev(get_palette("lancet", 5))
ctstbp <- ggbarplot(df, "SUBTYPE", "value", fill = type, color = type, 
                        palette = pal , ylab = "CELLS [%]", legend = "none", xlab = "")


# or matrix plot
mat <- round(mat / length(cts) * 100)
mat <- mat[,CELL.TYPES]

cnames <- c("SUBTYPE", "CELL.TYPE")
high.color <- "darkgreen"

matrixPlot <- function(mat, cnames, high.color = "red", bw.thresh = 50, gthresh = 0)
{
    lev1 <- colnames(mat)
    lev2 <- rownames(mat)
    df <- reshape2::melt(mat)
    col <- ifelse(df$value > bw.thresh, "white", "black")
    df$color <- ifelse(df$value > gthresh, col, "darkgrey")

    colnames(df)[1:2] <- cnames

    ggplot(df, aes(get(cnames[1]), get(cnames[2]))) +
        geom_tile(data=df, aes(fill=value), color="white") +
        scale_fill_gradient2(low="white", high=high.color, guide=FALSE) +
        geom_text(aes(label=value), size=5, color=df$color, fontface="bold") +
        theme(text = element_text(size = 14), 
            axis.text = element_text(color = "black", size=12),
            axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    #coord_equal() + 
    xlab(cnames[1]) + ylab(cnames[2])  
}

##
### Margin scores
##
spl.margin <- function(sce, col = "celltype")
{ 
    spl <- split(sce$margin, sce[[col]])
    df <- reshape2::melt(spl)
    
}
psces <- lapply(sces, spl.margin) 
for(i in tumors) psces[[match(i, tumors)]]$tumor <-  paste0("Tumor", i)
df <- do.call(rbind, psces)
colnames(df) <- c("MARGIN", "CELL.TYPE", "TUMOR")
df[[2]] <- factor(df[[2]], levels = CELL.TYPES)

# TODO: supplementary figure facet margin tsne + boxplot 
bp <- ggboxplot(df, x = "CELL.TYPE", y = "MARGIN", width = 0.8, notch = TRUE, 
            fill = "CELL.TYPE", palette = get_palette("lancet", 5),
            facet.by = "TUMOR", nrow = 1, 
            x.text.angle = 45, legend="none", ggtheme = theme_bw(base_size = 12)) 
fp <- facetTumors(sces, "margin", cont = TRUE)
ggarrange(fp, bp, align="hv", nrow=2)

# only epi
df <- subset(df, CELL.TYPE == "EPI")
df <- subset(df, TUMOR != "T61")
bp <- ggboxplot(df, x = "TUMOR", y = "MARGIN", width = 0.8, notch = TRUE, fill = cb.red, 
            legend="none", ggtheme = theme_bw(), ylab = "MARGIN (EPI CELLS)",
            font.x = 14, font.y = 14, font.tickslab = c(13, "plain", "black")) 
            # ylim = c(0,1),x.text.angle = 45, 
df.points <- data.frame(x=1:4, y=c(0.376, 0.766, 0.182, 0.646))
bp <- bp + geom_point(data=df.points, aes(x=x, y=y), shape = 8, size = 3) +
           geom_text(data=df.points[2,], aes(x=x, y=y), 
                         label="Bulk", hjust=0, nudge_x = 0.1, size=5)


# with bulk
mat <- t(rf.probs)
rownames(mat) <- sub("_consensus", "", rownames(mat))
mat <- round(mat, digits=2)
mat <- mat[c(3,4,1,2),c(1,3,2,4)]
mp <- matrixPlot(t(mat), cnames=c("TUMORS", "CLASS.PROB"), bw.thresh=0.5, gthresh=0.2)
mp <- mp + theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank(),
                          axis.ticks.x = element_blank(),
                           axis.line.x = element_blank(),
                            plot.margin = margin(b = 0))
mp <- mp + ylab("BULK")
p0 <- ggarrange(ctstbp, mp, bp, nrow=3, align="v", heights = c(2,1,3), labels = c("E", "F"))

##
### Heatmaps
##
library(consensusOV)
library(EnrichmentBrowser)
getClassProbs <- function(i)
{
    message(i)
    sce <- sces[[match(i, tumors)]]
    sce.entrez <- idMap(sce, org="hsa", from="SYMBOL", to="ENTREZID")
    cst <- get.consensus.subtypes(as.matrix(assays(sce.entrez)$logcounts), names(sce.entrez))
    probs <- cst$rf.probs
    colData(sce)[,colnames(probs)] <- probs
    saveRDS(sce, file=gsub("XX", i, "tumorXX/sampleXX_sce.rds")) 
    sce
}
sces <- lapply(tumors, getClassProbs)

library(ComplexHeatmap)
scHeatmap <- function(sce, cell.type = "EPI", subtype = c("DIF", "PRO"), 
        name = "log2TPM", title = "Cells", add = FALSE, fsize = 12, legend=TRUE)
{
    #sce <- idMap(sce, org="hsa", from="SYMBOL", to="ENTREZID")
    sce <- sce[,!is.na(sce$celltype)]
    sce <- sce[,sce$celltype %in% cell.type]
    sce <- sce[,sce$subtype %in% subtype]
    
    clgenes.sym <- AnnotationDbi::mapIds(org.Hs.eg.db,
                    column="SYMBOL", keys=clgenes, keytype="ENTREZID")
    sce <- sce[intersect(clgenes.sym, names(sce)),]
    
    # cell type
    ccol <- rev(get_palette("lancet", length(cell.type)))
    ct <- as.factor(sce$celltype)
    names(ccol) <- levels(ct)

    # subtype calls
    st <- as.factor(sce$subtype)

    # class probs per subtype        
    dat <- colData(sce)[, paste0(SUBTYPES, "_consensus")]
    dat <- as.matrix(dat)
    colnames(dat) <- sub("_consensus$", "", colnames(dat))
    if(add) colnames(dat) <- paste0(colnames(dat), "2") 

    cp.ramp <- circlize::colorRamp2(
                    seq(quantile(dat, 0.01), quantile(dat, 
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))
 
    # put together
    scol <- list(Subtype=stcols[subtype], CellType=ccol)
    for(i in seq_len(ncol(dat)))
    {   
        s <- colnames(dat)[i]
        scol[[s]] <- cp.ramp 
    }   

    df <- data.frame(CellType = ct, Subtype = st, dat)
    if(add) names(scol)[3] <- colnames(df)[3] <- "ClassProb"   
 
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, 
            col=scol, 
            show_legend = c(rep(TRUE, ifelse(add,3,2)), rep(FALSE,ifelse(add,3,4))),
            show_annotation_name = c(rep(TRUE,ifelse(add,3,0)), rep(FALSE,ifelse(add,3,6))),
            annotation_name_offset = unit(2, "mm"),
            gap = unit(c(0, 2, 0, 0, 0), "mm"))

    expr <- as.matrix(assays(sce)$logcounts)
    #rsums <- rowSums(expr)
    #top50 <- names(rsums)[tail(order(rsums), n=50)] 
    #expr <- expr[top50, ]

    #ComplexHeatmap::draw(
        ComplexHeatmap::Heatmap(expr, name=name, top_annotation = ha,
                            show_row_names=TRUE, show_column_names=FALSE,
                            column_title=title, row_title="Genes", 
                            row_names_gp = gpar(fontsize = fsize),
                            show_heatmap_legend = legend)
    #)

    #if(!add)
    #    for(i in seq_len(ncol(dat)))
    #        ComplexHeatmap::decorate_annotation(colnames(df)[i+2],
    #            {grid::grid.text(colnames(dat)[i], grid::unit(-2, "mm"), just = "left")})
}

h1 <- scHeatmap(sces[[1]], cell.type = "EPI", 
                subtype = c("DIF", "PRO"), title = "Cells (T59)", fsize=8)
h2 <- scHeatmap(sces[[4]], cell.type = "EPI", subtype = c("DIF", "PRO"), 
                title = "Cells (T77)", fsize=8, name="h2", add=TRUE, legend=FALSE)
for(st in names(stcols)) ComplexHeatmap::decorate_annotation(st, 
                            {grid::grid.text(st, grid::unit(-2, "mm"), just = "right")})
