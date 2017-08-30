############################################################
# 
# author: Ludwig Geistlinger and Levi Waldron
# date: 2017-08-12 22:31:59
# 
# descr: read and process ABSOLUTE calls
# 
############################################################

## SETUP

if(!require(Biobase) || packageVersion("Biobase") < "2.37")
  stop("Please install the development version of Bioconductor - this is cutting-edge stuff!")
if(!require(RTCGAToolbox) || packageVersion("RTCGAToolbox") < "2.7.5")
  BiocInstaller::biocLite("LiNk-NY/RTCGAToolbox")
if(!require(TCGAutils) || packageVersion("TCGAutils") < "0.3.28")
  BiocInstaller::biocLite("waldronlab/TCGAutils")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(RaggedExperiment)
  library(MultiAssayExperiment)
  library(dplyr)
  library(TCGAutils)
  library(readr)
})

## cancer types
ctypes <- getFirehoseDatasets() 

## Read ABSOLUTE calls and transform to GRangesList

### if already processed
absGRL <- readRDS("data/ABSOLUTE/ABSOLUTE_grangeslist.rds")
grlOV <- readRDS("data/ABSOLUTE/ABSOLUTE_OV_grangeslist.rds")

### processing raw data
isSubClonal <- function(x) as.integer(x$Subclonal_HSCN_a1 | x$Subclonal_HSCN_a2)

totalCN <- function(x) x$Modal_HSCN_1 + Modal_HSCN_2

# @returns: a GRangesList with per-sample ABSOLUTE calls
readAbsolute <- function(file="data/ABSOLUTE/TCGA_mastercalls.abs_segtabs.fixed.txt.bz2")
{
    con <- bzfile(file)
    df <- readr::read_tsv(con, col_names = TRUE)
    df <-  filter(df, !is.na(Chromosome))
    df <-  filter(df, Chromosome != 23)
    df <- as.data.frame(df, stringsAsFactors=FALSE)
    df[,"Chromosome"] <- paste0("chr", df[,"Chromosome"])

    ind2n <- which(df[,"Modal_HSCN_1"] == 1 & df[,"Modal_HSCN_2"] == 1)
    df <- df[-ind2n,]
    df$score <- isSubClonal(df)

    rel.cols <- c("Sample", "Chromosome", "Start", 
        "End", "Modal_HSCN_1", "Modal_HSCN_2", "score")
    df <- df[,rel.cols]
    grl <- makeGRangesListFromDataFrame(df, split.field = "Sample", keep.extra.columns = TRUE) 
    saveRDS(grl, file="data/ABSOLUTE/ABSOLUTE_grangeslist.rds")
    return(grl)
}

## GISTIC
# @returns: a RangedSummarizedExperiment
gistic2RSE <- function(ctype="OV", peak=c("wide", "narrow", "full"))
{

    BROAD.URL <- paste0("http://gdac.broadinstitute.org/",
                "runs/analyses__2016_01_28/data/ctype/20160128")
 
    GISTIC.FILE <- paste0("gdac.broadinstitute.org_ctype-TP.",
                    "CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz")
    
    # download the tar
    url <- file.path(BROAD.URL, GISTIC.FILE)
    url <- gsub("ctype", ctype, url)
    if(ctype == "LAML") url <- sub("TP", "TB", url)
    else if(ctype == "SKCM") url <- sub("TP", "TM", url)
    dest.file <- paste0("data/GISTIC/", ctype, ".tar.gz")
    download.file(url, dest.file)
    files <- untar(dest.file, list=TRUE)
    
    # extract the lesions file
    basef <- "all_lesions.conf_99.txt"
    basef <- files[grepl(basef, files)]
    untar(dest.file, files=basef)

    # read the lesion file
    gistic <- read.delim(basef, as.is=TRUE)
    file.remove(dest.file)    
    unlink(dirname(basef), recursive=TRUE, force=TRUE)

    # transform to RangedSummarizedExperiment
    rel.rows <- grepl("Peak +[0-9]+$", gistic[,1])
    gistic <- gistic[rel.rows,]

    # (a) get the ranges from chosen peaks
    peak <- match.arg(peak)
    peak.col <- grep(peak, c("wide", "narrow", "full")) + 2
    ranges <- gistic[rel.rows,peak.col]
    ranges <- sub("\\(probes [0-9]+:[0-9]+\\) *$", "", ranges)   
    ranges <- as(ranges, "GRanges")
    ind <- seqnames(ranges) != "chrX"
    ranges <- ranges[ind]
    gistic <- gistic[as.logical(ind),]
    ind <- orderSeqlevels(seqlevels(ranges))
    seqlevels(ranges) <- seqlevels(ranges)[ind]
    ind <- order(ranges)
    ranges <- ranges[ind]
    gistic <- gistic[ind, ]
  
    # (b) get the peak type (amplification / deletion) 
    peak.type <- sapply(gistic[,1], 
        function(x) unlist(strsplit(x, " "))[1], USE.NAMES=FALSE)
    
    # (c) create the SE
    rel.cols <- grepl("^TCGA", colnames(gistic))
    gistic <- gistic[rel.cols]
    gistic <- as.matrix(gistic)
    gisticSE <- SummarizedExperiment(assays=list(counts=gistic))
    rowRanges(gisticSE) <- ranges
    rowData(gisticSE)$type <- peak.type 

    # ensure consistent naming with subtypes and absolute
    colnames(gisticSE)<- gsub("\\.", "-", colnames(gisticSE))
    colnames(gisticSE) <- TCGAutils::TCGAbarcode(colnames(gisticSE), sample=TRUE)
    colnames(gisticSE)<- sub("A$", "", colnames(gisticSE))

    return(gisticSE)
}


## Broad subtypes
# @returns: a matrix with sample IDs as rownames and at least a column "cluster"
getBroadSubtypes <- function(ctype="OV", clust.alg=c("CNMF", "Consensus_Plus"))
{
    BROAD.URL <- paste0("http://gdac.broadinstitute.org/runs/analyses__latest",
        "/reports/cancer/ctype-TP/mRNA_Clustering_calg/ctype-TP.bestclus.txt")

    # insert selected cancer type
    url <- gsub("ctype", ctype, BROAD.URL)
    if(ctype == "LAML") url <- sub("TP", "TB", url)
    else if(ctype == "SKCM") url <- sub("TP", "TM", url)

    # insert selected cluster algorithm
    clust.alg <- match.arg(clust.alg)
    url <- sub("calg", clust.alg, url)
    
    # is mRNA cLustering available?
    if(!RCurl::url.exists(url))
    {
        warning(paste("mRNA clustering not available for", 
            ctype, "- mRNAseq clustering is taken instead"))
        url <- sub("mRNA", "mRNAseq", url)
    }

    subtys <- read.delim(url, skip=1, as.is=TRUE)
    subtys[,1] <- TCGAutils::TCGAbarcode(subtys[, 1], sample=TRUE)
    subtys[,1] <- sub("A$", "", subtys[,1])
    rownames(subtys) <- subtys[,1]
    subtys <- subtys[,-1] 
    return(subtys)
}


# subtys: a matrix with sample IDs as rownames and at least a column 'cluster'
# absGRL: a GRangesList storing the per-sample ABSOLUTE calls
# @returns: a RaggedExperiment storing the ABSOLUTE calls that 
getMatchedAbsoluteCalls <- function(absGRL, subtys, save=FALSE, ctype="OV")
{
    subtys <- subtys[rownames(subtys) %in% names(absGRL), ]
    absGRLmatched <- absGRL[match(rownames(subtys), names(absGRL))]
    ra <- RaggedExperiment(absGRLmatched, colData=subtys)

    if(save) 
    {
        outf <- paste("ABSOLUTE", ctype, "RaggedExp.rds",  sep="_") 
        outf <- file.path("data/ABSOLUTE/", outf)
        saveRDS(ra, file=outf)
    }
    return(ra)
}

# gistic: a RangedSummarizedExperiment
# subtys: a matrix with sample IDs as rownames and at least a column 'cluster'
testSubtypes <- function(gistic, subtys, 
    test.type=c("chisq", "perm"), padj.method="BH")
{
    test.type <- match.arg(test.type)
    subtys <- subtys[rownames(subtys) %in% colnames(gistic), ]
    gistic <- gistic[,match(rownames(subtys), colnames(gistic))]
    subtys <- subtys$cluster    
    slot <- ifelse(test.type == "perm", "statistic", "p.value")   
 
    res <- apply(assay(gistic), 1, 
        function(x) chisq.test(x, subtys)[[slot]])
    if(test.type == "chisq") res <- p.adjust(res, method=padj.method)       
    return(res)
}


# two ways of summarizing calls:
# 1:
maxScore <- function(scores, ranges, qranges) max(scores, na.rm=TRUE)
# 2:
wmean <- function(scores, ranges, qranges)
{
    isects <- pintersect(ranges, qranges)
    s <- sum(scores * width(isects)) / sum(width(isects))
    return(round(s))
}

# @ra: RaggedExperiment
# @query: GRanges
querySubclonality <- function(ra, query, sum.method=c("any", "wmean"))
{
    sum.method <- match.arg(sum.method)
    sum.method <- ifelse(sum.method == "wmean", wmean, maxScore)
    qa <- RaggedExperiment::qreduceAssay(ra, query, 
                simplifyReduce=sum.method, i="score", background=0)
    return(qa)
}

perm.test <- function()
replicate(1:1000,
    { 
        ind <- sample(nrow(ovsubs))
        ovsubs.perm <- ovsubs[ind,]
        rownames(ovsubs.perm) <- rownames(ovsubs)
        stat
    })

