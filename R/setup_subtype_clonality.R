############################################################
# 
# author: Ludwig Geistlinger and Levi Waldron
# date: 2017-08-12 22:31:59
# 
# descr: read and process ABSOLUTE calls
# 
############################################################

## SETUP
ctypes <- RTCGAToolbox::getFirehoseDatasets()

## cancer types
.isSubClonal <- function(x) as.integer(x$Subclonal_HSCN_a1 | x$Subclonal_HSCN_a2)
.totalCN <- function(x) x$Modal_HSCN_1 + x$Modal_HSCN_2

### processing raw data
# @returns: a GRangesList with per-sample ABSOLUTE calls
# infile = TCGA_mastercalls.abs_segtabs.fixed.txt.bz2
# outfile = ABSOLUTE_grangeslist.rds
.readAbsolute <- function(infile, outfile)
{
    con <- bzfile(infile)
    df <-  read.delim(con, as.is=TRUE)
    df <-  df[!is.na(df[,"Chromosome"]),]
    df <-  df[df[,"Chromosome"] != 23,]
    df[,"Chromosome"] <- paste0("chr", df[,"Chromosome"])

    ind2n <- which(df[,"Modal_HSCN_1"] == 1 & df[,"Modal_HSCN_2"] == 1)
    df <- df[-ind2n,]
    df$score <- .isSubClonal(df)

    rel.cols <- c("Sample", "Chromosome", "Start", 
        "End", "Modal_HSCN_1", "Modal_HSCN_2", "score")
    df <- df[,rel.cols]
    grl <- makeGRangesListFromDataFrame(df, 
        split.field="Sample", keep.extra.columns=TRUE) 
    saveRDS(grl, file=outfile)
    return(grl)
}


# @returns: a matrix with sample IDs as rownames and columns 'purity', 'ploidy',
#           'Genome.doublings',  and 'Subclonal.genome.fraction'  
.readPurityPloidy <- function(infile, outfile)
{
    con <- bzfile(infile)
    dat <- read.delim(infile, as.is=TRUE)
    dat <- dat[,-c(2,3,7,8,10)]
    rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    saveRDS(dat, file=outfile)
    return(dat)
}


## GISTIC
# @returns: a RangedSummarizedExperiment
gistic2RSE <- function(ctype=RTCGAToolbox::getFirehoseDatasets(), 
    peak=c("wide", "narrow", "full"), cache=TRUE)
{
    ctype <- match.arg(ctype)
    peak <- match.arg(peak)
    rname <- paste("gistic", ctype, peak, sep="_")
    
    if(cache)
    {
        gisticSE <- GSEABenchmarkeR:::.getResourceFromCache(
            rname, update.value=NA, ucdir="subtypeHeterogeneity")
        if(!is.null(gisticSE)) return(gisticSE)
    }
    
    # download the tar
    BROAD.URL <- paste0("http://gdac.broadinstitute.org/",
                "runs/analyses__2016_01_28/data/ctype/20160128")
 
    GISTIC.FILE <- paste0("gdac.broadinstitute.org_ctype-TP.",
                    "CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz")

    url <- file.path(BROAD.URL, GISTIC.FILE)
    url <- gsub("ctype", ctype, url)
    if(ctype == "LAML") url <- sub("TP", "TB", url)
    else if(ctype == "SKCM") url <- sub("TP", "TM", url)
    
    dest.file <- paste0(ctype, "_gistic2.tar.gz")
    out.dir <- system.file("extdata", package="subtypeHeterogeneity")
    dest.file <- file.path(out.dir, dest.file)
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

    GSEABenchmarkeR::cacheResource(gisticSE, rname, ucdir="subtypeHeterogeneity")

    return(gisticSE)
}


## Broad subtypes
# @returns: a matrix with sample IDs as rownames and at least a column "cluster"
getBroadSubtypes <- function(ctype=RTCGAToolbox::getFirehoseDatasets(), 
    clust.alg=c("CNMF", "Consensus_Plus"), cache=TRUE)
{
    ctype <- match.arg(ctype)
    clust.alg <- match.arg(clust.alg)
    rname <- paste("subtypes", ctype, clust.alg, sep="_")
    if(cache)
    {
        subtys <- GSEABenchmarkeR:::.getResourceFromCache(
            rname, update.value=NA, ucdir="subtypeHeterogeneity")
        if(!is.null(subtys)) return(subtys)
    }

    BROAD.URL <- paste0("http://gdac.broadinstitute.org/runs/analyses__latest",
        "/reports/cancer/ctype-TP/mRNA_Clustering_calg/ctype-TP.bestclus.txt")

    # insert selected cancer type
    url <- gsub("ctype", ctype, BROAD.URL)
    if(ctype == "LAML") url <- gsub("TP", "TB", url)
    else if(ctype == "SKCM") url <- gsub("TP", "TM", url)

    # insert selected cluster algorithm
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
    
    GSEABenchmarkeR::cacheResource(subtys, rname, ucdir="subtypeHeterogeneity")
    return(subtys)
}

# @args
#  subtys: a matrix with sample IDs as rownames and at least a column 'cluster'
#  absGRL: a GRangesList storing the per-sample ABSOLUTE calls
#
# @returns: a RaggedExperiment storing the ABSOLUTE calls that 
getMatchedAbsoluteCalls <- function(absGRL, subtys, max.cn=Inf)
{
    subtys <- subtys[rownames(subtys) %in% names(absGRL), ]
    absGRLmatched <- absGRL[match(rownames(subtys), names(absGRL))]
    
    if(max.cn < Inf)
        absGRLmatched <- endoapply(absGRLmatched, 
            function(x) x[.totalCN(x) <= max.cn])

    subtys <- S4Vectors::DataFrame(subtys)
    RaggedExperiment::RaggedExperiment(absGRLmatched, colData = subtys)
}

# @args
#  broad.subtys: a matrix with sample IDs as rownames and a column 'cluster'
#  pooled.subtys: data.frame from consensusOV manuscript 
#                   (OV subtypes from different studies)  
#
# @returns: a 4x4 contingeny table 
#            (Broad: {1,2,3,4} vs Verhaak: {PRO, MES, DIF, IMR})
mapOVSubtypes <- function(broad.subtys, pooled.subtys)
{
    ind <- pooled.subtys[,"data.source"] == "TCGA"
    pooledTCGA <- pooled.subtys[ind, "Verhaak"]
    ids <- rownames(pooled.subtys)[ind]
    ids <- sub("^TCGA.TCGA_", "", ids)
    ids <- gsub("\\.", "-", ids)
    
    broad.ids <- suppressWarnings( TCGAutils::TCGAbarcode(rownames(broad.subtys)) )
    ind <- match(ids, broad.ids)   
    broad.subtys <- broad.subtys[ind,"cluster"]
       
    spl <- split(pooledTCGA, broad.subtys)
    tab <- sapply(spl, table)
    return(tab)   
}

# requires a column named 'band' in cbands
annotateCytoBands <- function(gistic, cbands)
{
    olaps <- GenomicRanges::findOverlaps(rowRanges(gistic), cbands)
    sh <- S4Vectors::subjectHits(olaps)
    qh <- S4Vectors::queryHits(olaps)
    olist <- GRangesList(split(cbands[sh], qh))
    isects <- GenomicRanges::pintersect(olist, rowRanges(gistic))
    ind <- IntegerList(lapply(width(isects), which.max)) 
    bands <- unlist(isects[ind])$band
    mcols(gistic)$band <- bands
    return(gistic)
}

tissueScore <- function(emat)
{
    # get the signature 
    sig.file <- system.file("extdata/FT_vs_OSE_signature.txt", 
                                package="subtypeHeterogeneity")

    sig <- read.delim(sig.file, as.is=TRUE)
    sig <- split(sig[,"GENE"], sig[,"TISSUE"])
    
    # warn if certain signature genes are not present in the data
    not.ft <- sum(!(sig$FT %in% rownames(emat)))
    if(not.ft) warning(paste(not.ft, "of", 
                                length(sig$FT), "FT genes not present"))

    not.ose <- sum(!(sig$OSE %in% rownames(emat)))    
    if(not.ose) warning(paste(not.ose, "of", 
                                length(sig$OSE), "OSE genes not present"))

    # scoring for one sample at a time
    .scoreSample <- function(evals)
    {
        evals <- sort(evals)
        .getRank <- function(genes) match(genes, names(evals)) / length(evals)
        r <- lapply(sig, .getRank)
        nna <- sum(!is.na(unlist(r)))
        rs <-  sum(r[[1]], na.rm=TRUE) + sum(1 - r[[2]], na.rm=TRUE) 
        s <- 2 * rs / nna - 1
        return(s)
    }

    # apply to each sample & assign tissue of origin by thresholding
    scores <- apply(emat, 2, .scoreSample)
    tissues <- ifelse(scores > 0, "FT", "OSE")
    tscores <- S4Vectors::DataFrame(tissues = tissues, scores = scores)

    return(tscores)
} 
