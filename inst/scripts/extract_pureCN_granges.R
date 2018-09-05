############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-09-04 15:17:28
# 
# descr: extract pureCN calls from PureCN WES run 
# 
############################################################

map.file <- "/data/16tb/CNVworkflow/S0293689/sampleMap_S0293689.csv"
sample.map <- read.csv(map.file)
n <- as.vector(sample.map[,"filename"])
n <- sub(".bam$", "", n)
sample.map <- as.vector(sample.map[,"barcode"])
sample.map <- TCGAutils::TCGAbarcode(sample.map, sample=TRUE)
sample.map <- substring(sample.map, 1, 15) 
names(sample.map) <- n


data.dir <- "/data/16tb/CNVworkflow/S0293689/purecn_output/S0293689_PureCN/matching_normal"
ddirs <- list.dirs(data.dir, full.names=FALSE, recursive=FALSE)

rds.files <- file.path(data.dir, ddirs, paste0(ddirs, ".rds"))
ddirs <- ddirs[file.exists(rds.files)]

extractCalls <- function(d)
{
    message(d)
    rds.file <- paste0(d, ".rds") 
    rds.file <- file.path(data.dir, d, rds.file)  
    dat <- readRDS(rds.file)          

    pos <- dat$results[[1]]$seg[,c("chrom", "loc.start", "loc.end")]
    colnames(pos)[2:3] <- sub("^loc.", "", colnames(pos)[2:3])
    colnames(pos)[1] <- "seqnames"
    pos[,1] <- paste0("chr", pos[,1])
    mlc <- dat$results[[1]]$C.posterior$ML.C
    mlsubcl <- dat$results[[1]]$C.posterior$ML.Subclonal
    mlsubcl <- as.integer(mlsubcl)
    pos <- data.frame(pos, CN=mlc, score=mlsubcl)
    gr <- GenomicRanges::makeGRangesFromDataFrame(pos, keep.extra.columns=TRUE) 
    return(gr)
}

grl <- GenomicRanges::GRangesList(lapply(ddirs, extractCalls))
names(grl) <- unname(sample.map[ddirs])

###


