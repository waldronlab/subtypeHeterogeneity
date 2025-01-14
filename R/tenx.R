############################################################
# 
# author: Ludwig Geistlinger
# date: 2020-04-18 12:22:34
# 
# descr: analysis and visualization of 10X tumors
# 
############################################################

SUBTYPES <- c("DIF", "IMR", "MES", "PRO")
CELL.TYPES <- c("EPI", "LYMPH", "MYE", "STROM", "ENDO")

#
# analysis
#
preprocessSCE <- function(sce)
{
    bp <- BiocParallel::registered()[[1]]
    ave <- scater::calculateAverage(sce, BPPARAM = bp)
    rowData(sce)$AveCount <- ave
    to.keep <- ave > 0.001
    sce <- sce[to.keep,]

    clusters <- scran::quickCluster(sce, method = "igraph",
                                    min.mean = 0.1, BPPARAM = bp)
    sce <- scran::computeSumFactors(sce, min.mean = 0.1,
                                    cluster = clusters, BPPARAM = bp)
    sce <- scater::logNormCounts(sce)
    dec.pbmc <- scran::modelGeneVarByPoisson(sce, BPPARAM = bp)
    top.pbmc <- scran::getTopHVGs(dec.pbmc, prop = 0.1)
    sce <- scran::denoisePCA(sce, subset.row = top.pbmc,
                             technical = dec.pbmc, BPPARAM = bp)

    sce <- scater::runTSNE(sce, dimred = "PCA", perplexity = 30, BPPARAM = bp)
    sce$sizeFactor <- sizeFactors(sce)
    snn.gr <- scran::buildSNNGraph(sce, use.dimred = "PCA", BPPARAM = bp, k = 25)
    clusters <- igraph::cluster_walktrap(snn.gr)
    sce$Cluster <- factor(clusters$membership)
    sce
}   

compileMarkers <- function(sc.markers, bulk.markers, nr.markers = 10)
{
    grid <- seq_len(nr.markers)
    
    sc.dif.markers <- rownames(sc.markers$DIF[sc.markers$DIF$logFC.PRO > 0,])
    sc.pro.markers <- rownames(sc.markers$DIF[sc.markers$DIF$logFC.PRO < 0,])
    
    sc.dif.markers <- sc.dif.markers[sc.dif.markers %in% bulk.markers$DIF]
    sc.pro.markers <- sc.pro.markers[sc.pro.markers %in% bulk.markers$PRO]

    sc.markers$DIF[c(sc.dif.markers[grid], sc.pro.markers[grid]),]         
}

getCommonMarkers <- function(marker.list, subtype = c("DIF", "PRO"), n = 5)
{
    subtype <- match.arg(subtype)
    .subsetMarkers <- function(markers)
    {
        mst <- markers[["DIF"]]
        if(subtype == "DIF") ind <- mst$logFC.PRO > 0
        else ind <- mst$logFC.PRO < 0
        rownames(mst[ind,])
    }
    st.markers <- lapply(marker.list, .subsetMarkers)
    rn <- Reduce(intersect, st.markers)
    ranks <- lapply(st.markers, function(markers) match(rn, markers))
    cn <- do.call(cbind, ranks)
    rmeans <- rowMeans(cn)
    head(rn[order(rmeans)], n) 
}


generateInferCNVInput <- function(sce, sample.dir)
{
    # (1) gene order file
    gorder.file <- paste0("sample", snr, "_gene_order.tab")
    gorder.file <- file.path(sample.dir, gorder.file)

    orgdb <- org.Hs.eg.db::org.Hs.eg.db
    eids <- AnnotationDbi::mapIds(orgdb, 
                                  column = "ENTREZID",
                                  keys = rownames(sce),
                                  keytype = "SYMBOL")

    ind <- !is.na(eids)
    sce <- sce[ind,]
    eids <- eids[ind]
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    ind <- unname(eids) %in% names(GenomicFeatures::genes(txdb))
    sce <- sce[ind,]
    eids <- eids[ind]

    egenes <- GenomicFeatures::genes(txdb)[unname(eids)]
    names(egenes) <- names(eids)
    egenes <- sort(egenes)
    egenes <- GenomeInfoDb::keepStandardChromosomes(egenes, pruning.mode = "coarse")
    sce <- sce[names(egenes), ]
    
    egenes <- as.data.frame(egenes)
    write.table(egenes[,1:3], file = gorder.file, sep = "\t", 
                quote = FALSE, col.names = FALSE)


    # (2) count matrix file
    cmat.file <- paste0("sample", snr, "_count_matrix.tab")
    cmat.file <- file.path(sample.dir, cmat.file)
    cmat <- as.matrix(assay(sce))
    colnames(cmat) <- sce$Barcode
    write.table(cmat, file = cmat.file, sep = "\t", quote = FALSE)
    
    # (3) sample annotation file   
    sanno.file <- paste0("sample", snr, "_sample_anno.tab")
    sanno.file <- file.path(sample.dir, sanno.file)
    sce <- sce[,!is.na(sce$celltype)]
    sanno <- colData(sce)[,c("Barcode", "celltype")]
    #sanno[,2] <- ifelse(sanno[,2] == "EPI", "malignant", "normal")
    write.table(sanno, file = sanno.file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)

    # (4) generate input object
    ref.names <- unique(sce$celltype)
    ref.names <- ref.names[ref.names != "EPI"]
    infercnv::CreateInfercnvObject(raw_counts_matrix = cmat.file,
                                   annotations_file = sanno.file,
                                   delim = "\t",
                                   gene_order_file = gorder.file,
                                   ref_group_names = ref.names)
#                                  ref_group_names = "normal")
}

subsetObservations <- function(obs.mat, sce)
{
    sce <- sce[,colnames(obs.mat)]
    ind <- sce$subtype %in% c("DIF", "PRO")
    obs.mat[,ind]
}

getRecurrentDriverGenes <- function(symbols)
{
    symbols[match("ANKRD22", symbols)] <- "PTEN" 
    symbols[match("TRIT1", symbols)] <- "MYCL"
    cgenes <- getCnvGenesFromTCGA(excl = FALSE, build = "hg38")
    clabs <- vector("character", length = length(symbols))
    ind <- symbols %in% names(cgenes)
    clabs[ind] <- symbols[ind]
    clabs
}

getGISTIC <- function(egenes, obs.mat)
{
    gisticOV <- gistic2RSE(ctype="OV", peak="wide")
    #gdf <- as.data.frame(rowRanges(gisticOV))
    #gstr <- paste(gdf[,1], paste(gdf[,2], gdf[,3], sep = "-"), sep = ":")
    #cat(gstr, file = "gistic_hg19.txt", sep = "\n")
    gtype <- rowData(gisticOV)$type  

    not.mapped <- c(22, 31, 43, 47, 67) 
    gfile <- system.file("extdata", package = "subtypeHeterogeneity")
    gfile <- file.path(gfile, "gisticOV_ranges_hg38.txt")
    gstr <- scan(gfile, what = "character", sep = "\n")
    #gstr <- gstr[-not.mapped]
    gistic <- GenomicRanges::GRanges(gstr)
    gtype <- gtype[-not.mapped]

    stopifnot(all(names(egenes) == rownames(obs.mat)))
    olaps <- GenomicRanges::findOverlaps(egenes, gistic)
    sh <- S4Vectors::subjectHits(olaps)
    qh <- S4Vectors::queryHits(olaps)
    ampdel <- rep("Neutral", nrow(obs.mat))
    ampdel[qh] <- gtype[sh]
    ampdel <- substring(ampdel, 1, 3)  
    df <- data.frame(GISTIC2 = ampdel)
    col <- list(GISTIC2 = c(Amp = "red", Del = "blue", Neu = "white"))
    ComplexHeatmap::HeatmapAnnotation(df = df, col = col)
}

symbols2ranges <- function(symbols)
{  
    orgdb <- org.Hs.eg.db::org.Hs.eg.db
    eids <- AnnotationDbi::mapIds(orgdb, 
                                    column = "ENTREZID", 
                                    keys = symbols, 
                                    keytype = "SYMBOL")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    egenes <- GenomicFeatures::genes(txdb)[unname(eids)]
    egenes <- BiocGenerics::unstrand(egenes)
    names(egenes) <- names(eids)
    sort(egenes)
}

extendObservations <- function(om, un)
{
    mgenes <- setdiff(un, rownames(om))
    m <- matrix(1, nrow = length(mgenes), ncol = ncol(om))
    rownames(m) <- mgenes
    colnames(m) <- colnames(om)
    rbind(om, m)
}

getCellAnnotation <- function(sces, obs.mats, bw = FALSE)
{   
    for(i in seq_along(sces)) 
        sces[[i]] <- sces[[i]][,colnames(obs.mats[[i]])]
    sts <- do.call(c, lapply(sces, function(sce) sce$subtype))
    margins <- do.call(c, lapply(sces, function(sce) sce$margin))
    df <- data.frame(Subtype = sts, Margin = margins)

    if(bw)
    { 
        stcols <- bw.stcols
        margin.ramp <- circlize::colorRamp2(
                    seq(quantile(margins, 0.01), quantile(margins,
                    0.99), length = 2), c("white", "black"))
    }
    else margin.ramp <- circlize::colorRamp2(
                    seq(quantile(margins, 0.01), quantile(margins,
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))

    cols <- list(Subtype = stcols, Margin = margin.ramp)
    ComplexHeatmap::HeatmapAnnotation(df = df, col = cols, which = "row")
}

readObservations <- function(i)
{
    obs.file <- file.path(paste0("tumor", i), "infercnv.observations.txt")
    obs <- readr::read_delim(obs.file, delim = " ")
    obs.mat <- as.matrix(obs[,2:ncol(obs)])
    rownames(obs.mat) <- obs$GENE
    obs.mat
}

subsetByCellType <- function(sce, cell.type = "EPI", subtype = c("DIF", "PRO"))
{
    sce <- sce[,!is.na(sce$celltype)]
    sce <- sce[,sce$celltype %in% cell.type]
    sce <- sce[,sce$subtype %in% subtype]
    ind <- do.call(order, as.data.frame(colData(sce)[,c("subtype", "margin")]))
    sce[,ind]
}

getClassProbs <- function(i, tumors, save2file = FALSE)
{   
    message(i)
    sce <- sces[[match(i, tumors)]]
    sce.entrez <- EnrichmentBrowser::idMap(sce, org="hsa", from="SYMBOL", to="ENTREZID")
    am <- as.matrix(assays(sce.entrez)$logcounts)
    cst <- consensusOV::get.consensus.subtypes(am, names(sce.entrez))
    probs <- cst$rf.probs
    colData(sce)[,colnames(probs)] <- probs
    if(save2file) saveRDS(sce, file=gsub("XX", i, "tumorXX/sampleXX_sce.rds"))
    sce
}

getCelltypeSubtypeMatrix <- function(sces)
{
    cts <- unlist(lapply(sces, function(sce) sce$celltype))
    sts <- unlist(lapply(sces, function(sce) sce$subtype))
    nna.ind <- !is.na(cts)
    cts <- cts[nna.ind]
    sts <- sts[nna.ind]
    spl <- split(sts, cts)
    vapply(spl, table, numeric(4))
}

getPerc <- function(sce, col="celltype", perc=TRUE)
{
    col <- sce[[col]]
    col <- col[!is.na(col)]
    ctab <- table(col)
    if(perc) ctab <- round(ctab / sum(ctab) * 100)
    return(ctab)
}

transferCellType <- function(sce, ref = c("hpca", "encode"))
{ 
    ref <- match.arg(ref)
    if(ref == "hpca") se <- SingleR::HumanPrimaryCellAtlasData()
    else se <- SingleR::BlueprintEncodeData()
        
    common <- intersect(rownames(sce), rownames(se))
    se <- se[common,]
    sce <- sce[common,]

    SingleR::SingleR(test = sce, ref = se, 
                     labels = se$label.main,
                     assay.type.ref = "logcounts",
                     BPPARAM = BiocParallel::registered()[[1]])
}


annoCellType <- function(sce)
{
    hpca <- sce$hpca.celltype
    encode <- sce$encode.celltype
            
    sce$celltype <- vapply(seq_len(ncol(sce)),
                            function(i) .decideCellType(hpca[i], encode[i]),
                            character(1)) 
    return(sce)                 
}

.decideCellType <- function(hp, en) 
{   
    epithelial <- c("Epithelial_cells", "Epithelial cells", "Embryonic_stem_cells",
                    "MSC", "iPS_cells", "Keratinocytes", "Neuroepithelial_cell",
                    "Neurons")
    myeloid <- c("Macrophage", "Monocyte", "Macrophages", "Monocytes", "DC",
                    "Erythrocytes", "Neutrophils")
    lymphocyte <- c("T_cells", "B_cell", "Pre-B_cell_CD34-", "NK_cell", "CD8+ T-cells",
                    "CD4+ T-cells", "HSC", "B-cells", "NK cells", "Pro-B_cell_CD34+")
    endothelial <- c("Endothelial_cells", "Endothelial cells")
    stromal <- c("Fibroblasts", "Tissue_stem_cells", "Smooth_muscle_cells",
                    "Osteoblasts", "Chondrocytes", "Myocytes", "Fibroblasts", 
                    "Pericytes", "Chondrocytes", "Smooth muscle", "Adipocytes")
                        
    if(hp %in% epithelial || en %in% epithelial) ct <- "EPI"
    else if(hp %in% myeloid || en %in% myeloid) ct <- "MYE"
    else if(hp %in% lymphocyte || en %in% lymphocyte) ct <- "LYMPH"
    else if(hp %in% endothelial || en %in% endothelial) ct <- "ENDO"
    else if(hp %in% stromal || en %in% stromal) ct <- "STROM"
    else ct <- NA_character_
    return(ct)
}

#
# visualization
#
facetTumors <- function(sces, col = "subtype", pal = stcols, cont = FALSE, shape = 21)
{   
    if(!cont) for(i in seq_along(sces)) sces[[i]][[col]] <- factor(sces[[i]][[col]], 
                                                         levels = names(pal)) 
    psces <- lapply(sces, scater::plotTSNE, colour_by = col)
    psces <- lapply(psces, function(p) p$data)
    for(i in tumors) psces[[match(i, tumors)]]$tumor <-  paste0("Tumor", i)
    df <- do.call(rbind, psces)
    
    colnames(df)[3] <- col
    df <- df[!is.na(df[[col]]),] 
    if(!cont) shape <- col
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "X", y = "Y", fill = col, shape = 21, color = col)) + 
            ggplot2::geom_point(alpha = 0.7, size = 0.75, shape = 21) +
        ggplot2::facet_wrap(~tumor, nrow=1, ncol=5) +
        ggplot2::xlab("Dimension 1") +
        ggplot2::ylab("Dimension 2") +
        ggplot2::labs(fill = col) +
        ggplot2::theme_bw(base_size = 12)

    is.cont <- col %in% names(sces[[1]]) || cont
    if(is.cont) p <- p + viridis::scale_fill_viridis()
    else p <- p + ggplot2::scale_fill_manual(values=pal) +
                  ggplot2::scale_color_manual(values=pal) +
                  guides(fill = guide_legend(override.aes = list(size=7, alpha =1))) 
    p
}

## arrange subtype and celltype side-by-side
plotType <- function(sce, type = c("celltype", "subtype"))
{   
    pal <- "lancet"
    
    type <- match.arg(type)
    p <- scater::plotTSNE(sce, colour_by = type)
    df <- p$data
    type <- toupper(type)
    colnames(df)[3] <- type 
    df <- df[!is.na(df[[type]]),]
    
    if(type == "CELLTYPE") df[[type]] <- factor(df[[type]], levels = CELL.TYPES)
    else pal <- ggpubr::get_palette("npg", 4)[c(4,2,3,1)]
    
    p <- ggpubr::ggscatter(df, x = "X", y = "Y", shape = 21, size = 1.5, color = "grey", 
                    alpha = 0.8, fill = type, palette = pal, ggtheme = ggplot2::theme_bw())
    ggpubr::ggpar(p, xlab = "Dimension 1", ylab = "Dimension 2")
}

bpType <- function(sces, type=c("celltype", "subtype"), pal)
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
        #pal <- rev(ggpubr::get_palette("lancet", 5))
    }
    else
    {
        df[[type]] <- factor(df[[type]], levels = rev(SUBTYPES))
        #pal <- rev(ggpubr::get_palette("npg", 4)[c(4,2,3,1)])
    }

    ggpubr::ggbarplot(df, "TUMORS", "value", fill = type,
                color = type, palette = pal , ylab = "CELLS [%]")
}

bpCrossType <- function(mat, pal)
{
    perc <- round(mat / rowSums(mat) * 100, digits=1)
    df <- reshape2::melt(perc)
    type <- "CELL.TYPE"
    colnames(df) <- c("SUBTYPE", type, "value")
    df[[type]] <- factor(df[[type]], levels = rev(CELL.TYPES))
    #pal <- rev(ggpubr::get_palette("lancet", 5))
    ggpubr::ggbarplot(df, "SUBTYPE", "value", fill = type, color = type, 
                        palette = pal, legend = "none",
                        ylab = "CELLS [%]", xlab = "SUBTYPES")
}

matrixPlot <- function(mat, cnames, high.color = "red", bw.thresh = 50, gthresh = 0)
{  
    lev1 <- colnames(mat)
    lev2 <- rownames(mat)
    df <- reshape2::melt(mat)
    col <- ifelse(df$value > bw.thresh, "white", "black")
    df$color <- ifelse(df$value > gthresh, col, "darkgrey")
    
    colnames(df)[1:2] <- cnames

    ggplot2::ggplot(df, ggplot2::aes(get(cnames[1]), get(cnames[2]))) +
        ggplot2::geom_tile(data=df, ggplot2::aes(fill=value), color="white") +
        ggplot2::scale_fill_gradient2(low="white", high=high.color, guide=FALSE) +
        ggplot2::geom_text(aes(label=value), size=4, color=df$color, fontface="bold") +
        ggplot2::theme(text = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(color = "black", size=10),
            axis.text.x = ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    #coord_equal() + 
    ggplot2::xlab(cnames[1]) + ggplot2::ylab(cnames[2])
}

scHeatmap2 <- function(sce, cell.type = "EPI", subtype = c("DIF", "PRO"),
        name = "log2TPM", title = "Cells", add = "", fsize = 12, legend = TRUE,
        last = FALSE)
{
    #sce <- idMap(sce, org="hsa", from="SYMBOL", to="ENTREZID")
    sce <- sce[,!is.na(sce$celltype)]
    sce <- sce[,sce$celltype %in% cell.type]
    sce <- sce[,sce$subtype %in% subtype]

    clgenes <- getClassifierGenes()    
    clgenes.sym <- AnnotationDbi::mapIds(org.Hs.eg.db,
                    column="SYMBOL", keys=clgenes, keytype="ENTREZID")
    sce <- sce[intersect(clgenes.sym, names(sce)),]

    # cell type
    ccol <- rev(ggpubr::get_palette("lancet", length(cell.type)))
    ct <- as.factor(sce$celltype)
    names(ccol) <- levels(ct)

    # subtype calls
    st <- as.factor(sce$subtype)

    # class probs per subtype        
    dat <- colData(sce)[, paste0(SUBTYPES, "_consensus")]
    dat <- as.matrix(dat)
    colnames(dat) <- sub("_consensus$", "", colnames(dat))
    colnames(dat) <- paste0(colnames(dat), add)

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
    if(last) names(scol)[3] <- colnames(df)[3] <- "ClassProb"

    ha <- ComplexHeatmap::HeatmapAnnotation(df = df,
            col=scol,
            show_legend = c(rep(TRUE, ifelse(last,3,2)), rep(FALSE,ifelse(last,3,4))),
            show_annotation_name = c(rep(TRUE,ifelse(last,3,0)), rep(FALSE,ifelse(last,3,6))),
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

bpCyclins <- function(sce, top.only = FALSE, nr.cells = 100, bw = FALSE)
{
    sce <- subsetByCellType(sce)
    message(paste0("Total EPI-DIF:", sum(sce$subtype == "DIF")))
    message(paste0("Total EPI-PRO:", sum(sce$subtype == "PRO")))

    cyclin.genes <- grep("^CCN[ABDE][0-9]$", names(sce), value = TRUE)
    am <- as.matrix(assay(sce, "logcounts"))
    am <- am[sort(cyclin.genes),]  
    rownames(am) <- sub("^CCN", "", rownames(am))

    df <- reshape2::melt(am)
    df$subtype <- sce$subtype[df[,2]]
    df$margin <- sce$margin[df[,2]]
    ind.dif <- tail(which(sce$subtype == "DIF"), nr.cells)
    ind.pro <- tail(which(sce$subtype == "PRO"), nr.cells)
    top <-  paste0("Top", nr.cells)
    df$top <- ifelse(df[,2] %in% c(ind.dif, ind.pro), top, "All")
    colnames(df)[c(1,3)] <- c("Cyclins", "log2CPM")

    if(bw) stcols <- bw.stcols    

    if(top.only)
    { 
        df <- df[df$top == top,]
        ggpubr::ggboxplot(df, "Cyclins", "log2CPM", legend = "none", ggtheme = theme_bw(),
                          color = "subtype", palette = stcols[c("DIF", "PRO")])
    }
    else
    {
        df2 <- df[df$top == top,]
        df2$top <- "All"
        df <- rbind(df, df2)
        all.lab <- paste0("All epithelial cells (",
                            sum(sce$subtype == "DIF"), " DIF, ",
                            sum(sce$subtype == "PRO"), " PRO", ")")
        top.lab <- paste0("Top ", nr.cells, " highest margin cells (",
                            length(ind.dif), " DIF, ",
                            length(ind.pro), " PRO", ")")
 
        ggpubr::ggboxplot(df, "Cyclins", "log2CPM", color = "subtype",
                palette = stcols[c("DIF", "PRO")], facet.by = "top",
                nrow = 2, ncol = 1,
                panel.labs = list(top = c(all.lab, top.lab)))  
    }
}

markerHeatmap <- function(sce, markers, row.split = FALSE)
{
    sce <- subsetByCellType(sce)
    
    st <- sce$subtype
    ind <- c(rev(which(st == "DIF")), which(st == "PRO"))
    sce <- sce[,ind]
    
    am <- as.matrix(assay(sce, "logcounts"))
    am <- am[markers,] 
    ind <- colSums(am) > 0

    sts <- colData(sce)[ind, "subtype"]
    margins <- colData(sce)[ind, "margin"]
    df <- data.frame(Subtype = sts, Margin = margins)

    margin.ramp <- circlize::colorRamp2(
                    seq(quantile(margins, 0.01), quantile(margins,
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))

    col <- list(Subtype = stcols, Margin = margin.ramp)
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col = col)

    if(row.split)
    { 
        ststr <- rep(c("DIF", "PRO"), each = length(markers) / 2)
        row.split <- paste(ststr, "bulk markers")
    }
    else row.split <- NULL   

    am <- t(scale(t(am[markers,ind])))
    ComplexHeatmap::Heatmap(am, name = "Expression", top_annotation = ha,
                            cluster_rows = FALSE, cluster_columns = FALSE, 
                            row_split = row.split)
}

markerHeatmapList <- function(sces, markers, row.split = FALSE, bw = FALSE)
{
    sces <- lapply(sces, subsetByCellType)
    
    .orderCellsBySubtype <- function(sce)
    {
        st <- sce$subtype
        ind <- c(rev(which(st == "DIF")), which(st == "PRO"))
        sce[,ind]
    }
    sces <- lapply(sces, .orderCellsBySubtype)
    
    ams <- lapply(sces, function(sce) as.matrix(assay(sce, "logcounts"))[markers,])
    am <- do.call(cbind, ams)
    for(i in seq_along(sces)) sces[[i]] <- sces[[i]][,colSums(ams[[i]]) > 0]

    sts <- lapply(sces, function(sce) sce$subtype)
    sts <- unlist(sts)
    margins <- lapply(sces, function(sce) sce$margin)
    margins <- unlist(margins)
    df <- data.frame(Subtype = sts, Margin = margins)

    if(bw)
    { 
        stcols <- bw.stcols
        margin.ramp <- circlize::colorRamp2(
                    seq(quantile(margins, 0.01), quantile(margins,
                    0.99), length = 2), c("white", "black"))
        

    }
    else margin.ramp <- circlize::colorRamp2(
                    seq(quantile(margins, 0.01), quantile(margins,
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))

    col <- list(Subtype = stcols, Margin = margin.ramp)
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col = col, show_legend = FALSE)

    if(row.split)
        row.split <- rep(c("DIF", "PRO"), each = length(markers) / 2)
    else row.split <- NULL   

    col.split <- rep(names(sces), vapply(sces, ncol, integer(1)))

    am <- t(scale(t(am)))
    if(bw) main.ramp <- circlize::colorRamp2(
                    seq(quantile(as.vector(am), 0.01), quantile(as.vector(am),
                    0.99), length = 2), c("white", "black"))
    else main.ramp <- margin.ramp <- circlize::colorRamp2(
                    seq(quantile(as.vector(am), 0.01), quantile(as.vector(am),
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))

    ComplexHeatmap::Heatmap(am, name = "Expression", col = main.ramp, top_annotation = ha,
                            cluster_rows = FALSE, cluster_columns = FALSE, 
                            row_split = row.split, column_split = col.split,
                            show_heatmap_legend = FALSE)
}

plotMainCellTypes <- function(sce, col, nr.types = 6)
{
    main.cats <- sort(table(sce[[col]]), decreasing=TRUE)
    main.cats <- names(main.cats)

    grid <- seq_len(nr.types)
    sce.sub <- sce[,sce[[col]] %in% main.cats[grid]]
    sce.sub[[col]] <- factor(as.vector(sce.sub[[col]]), levels = main.cats[grid])
    scater::plotTSNE(sce.sub, colour_by = col, shape_by = col)
}

pseudotimeHeatmap <- function(sce, lineage = 1, nr.genes = 100, 
                              cluster.rows = FALSE, scale = TRUE)
{
    t <- sce[[paste("slingPseudotime", lineage, sep = "_")]]
    grid <- seq_len(nr.genes)
    Y <- as.matrix(assay(sce, "logcounts"))
    vartop <- names(sort(apply(Y, 1, var), decreasing = TRUE))[grid]

    Y <- Y[vartop,]
    
    gam.pval <- apply(Y, 1, 
        function(z)
        {
            d <- data.frame(z=z, t=t)
            suppressWarnings({
                lo <- gam::lo
                tmp <- suppressWarnings(gam::gam(z ~ lo(t), data=d))
            })
            summary(tmp)[3][[1]][2,3]
        })

    topgenes <- names(sort(gam.pval, decreasing = FALSE))[grid]
    heatdata <- Y[topgenes, order(t, na.last = NA)]
    heatclus <- sce$subtype[order(t, na.last = NA)]
    
    # plot the heatmap
    df <- data.frame(Subtype = heatclus)
    col <- list(Subtype = stcols)
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col = col, show_legend = FALSE)
    am <- as.matrix(heatdata)
    if(scale) am <- t(scale(t(am)))
    ComplexHeatmap::Heatmap(am, 
        name = "Expression", top_annotation = ha,
        show_heatmap_legend = FALSE,
        show_row_names = FALSE,
        cluster_rows = cluster.rows, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6), row_title = "Genes",
        column_title = "Cells")
}
