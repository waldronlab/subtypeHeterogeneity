############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-10-22 17:31:07
# 
# descr: functionality for scRNA-seq analysis 
# 
############################################################

# get Supplementary Table 8A from Verhaak et al, JCI, 2013
getExtendedVerhaakSignature <- function(verhaak.file, 
                                        nr.genes.per.subtype = 200,
                                        by.subtype = FALSE)
{
    dat <- read.delim(verhaak.file)
    dat <- dat[-c(1:2),1:13]
    rownames(dat) <- as.vector(dat[,1])
    dat <- dat[,-1]
    sts <- c("DIF", "IMR", "MES", "PRO")
    de <- c("F", "FC", "Q")
    nam <- vapply(sts, function(s) paste(s, de, sep="."), character(3))
    colnames(dat) <- as.vector(nam)
    dat <- as.matrix(dat)
    dat <- gsub(",", ".", dat)
    mode(dat) <- "numeric"

    .getGenes <- function(st)
    {
        rcol <- paste0(st, ".F")
        ind <- order(dat[,rcol], decreasing=TRUE)
        genes <- rownames(dat)[ind][seq_len(nr.genes.per.subtype)]
        return(genes)
    }

    st.genes <- lapply(sts, .getGenes)
    if(by.subtype) names(st.genes) <- sts
    else st.genes <- sort(unique(unlist(st.genes))) 
    return(st.genes)
}

scHeatmap <- function(sc.expr, subtypes, c2t.file="scRNAseq_celltypes.txt")
{
    # cell type: stromal / epithelial
    c2t.file <- file.path(data.dir, c2t.file)
    c2t <- read.delim(c2t.file, as.is=TRUE)
    n <- paste0("X", c2t[,1])
    c2t <- c2t[,2]
    names(c2t) <- n

    # cell type
    ccol <- c("#2AB68C", "#B62A84")
    names(ccol) <- unique(unname(c2t))
    ct <- as.factor(c2t[colnames(sc.expr)])

    # subtype calls
    st <- sub("_consensus$", "", subtypes$consensusOV.subtypes)
    st <- as.factor(st)

    # class probs per subtype        
    dat <- subtypes$rf.probs
    colnames(dat) <- sub("_consensus$", "", colnames(dat))
   
    cp.ramp <- circlize::colorRamp2(
                    seq(quantile(dat, 0.01), quantile(dat, 
                    0.99), length = 3), c("blue", "#EEEEEE", "red"))
 
    # put together
    scol <- list(Subtype=stcols, CellType=ccol)
    for(i in seq_len(ncol(dat)))
    {
        s <- colnames(dat)[i]
        scol[[s]] <- cp.ramp 
    }
    df <- data.frame(CellType=ct, Subtype = st, dat)
    names(scol)[3] <- colnames(df)[3] <- "ClassProb"
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, 
            col=scol, 
            show_legend = c(rep(TRUE,3), rep(FALSE,3)),
            show_annotation_name = c(rep(TRUE,3), rep(FALSE,3)),
            annotation_name_offset = unit(2, "mm"),
            gap = unit(c(0, 2, 0, 0, 0), "mm"))

    expr <- log(sc.expr + 1, base=2)
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(expr, name="log2TPM", top_annotation = ha,
                            show_row_names=FALSE, show_column_names=FALSE,
                            column_title="Cells", row_title="Genes"))
    
    for(i in seq_len(ncol(dat)))
        ComplexHeatmap::decorate_annotation(colnames(df)[i+2], 
            {grid::grid.text(colnames(dat)[i], grid::unit(-2, "mm"), just = "right")})
}

margin <- function(rf.probs)
{
    subtr <- apply(rf.probs, 1, function(row) sort(row)[3])
    pred.margins <- Biobase::rowMax(rf.probs) - subtr
    return(pred.margins)
}

getClassifierGenes <- function()
{
    fnames <- lapply(
        consensusOV:::esets.rescaled.classified.filteredgenes,
        Biobase::featureNames)
    clgenes <- sort(unique(unlist(fnames)))
    sub("^geneid.", "", clgenes)
}



