############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-10-22 17:31:07
# 
# descr: functionality for scRNA-seq analysis 
# 
############################################################

# get Supplementary Table 8A from Verhaak et al, JCI, 2013
getExtendedVerhaakSignature <- function(verhaak.file, nr.genes.per.subtype=200)
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
    st.genes <- sort(unique(unlist(st.genes))) 
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

    # subtype
    st <- sub("_consensus$", "", subtypes)
    st <- as.factor(st)

    # cell type
    ccol <- c("#2AB68C", "#B62A84")
    names(ccol) <- unique(unname(c2t))
    ct <- as.factor(c2t[colnames(sc.expr)])

    scol <- list(Subtype=stcols, CellType=ccol)
    df <- data.frame(Subtype = st, CellType=ct)
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col=scol)

    expr <- log(sc.expr + 1, base=2)
    ComplexHeatmap::Heatmap(expr, name="log2TPM", top_annotation = ha,
                            show_row_names=FALSE, show_column_names=FALSE,
                            column_title="Cells", row_title="Genes")
}

margin <- function(rf.probs)
{
    subtr <- apply(rf.probs, 1, function(row) sort(row)[3])
    pred.margins <- Biobase::rowMax(rf.probs) - subtr
    return(pred.margins)
}





