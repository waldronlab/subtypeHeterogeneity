############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-10-22 17:31:07
# 
# descr: 
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
