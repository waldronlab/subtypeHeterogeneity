############################################################
# 
# author: Ludwig Geistlinger
# date: 2017-08-22 17:34:32
# 
# descr: 
# 
############################################################

## Create RaggedExperiment for OVC samples
source("R/setup_subtype_clonality.R")

# access RaggedExperiment: assay(ra[1:5,1:5], "score")

## summarize absolute calls in gistic regions
maxScore <- function(scores, ranges, qranges) max(scores, na.rm=TRUE)

weightedmean2 <- function(scores, ranges, qranges)
{
    isects <- pintersect(ranges, qranges)
    sum(scores * width(isects)) / sum(width(isects))
}


querySubclonality <- function(ra, query)
{
    qa <- RaggedExperiment::qreduceAssay(ra, query=gistic, 
                simplifyReduce=maxScore, i="score", background=0)

}

# find gistic regions that are overlapped by more than one absolute call
g <- sapply(names(grlOV), 
    function(n)
    { 
        message(match(n, names(grlOV)))
        x <- grlOV[[n]]
        olaps <- findOverlaps(rowRanges(gisticOV), x)
        qh <- queryHits(olaps)
        sh <- subjectHits(olaps)
        dups <- duplicated(qh)
        dups <- unique(qh[dups])

        res <- sapply(dups, 
            function(d)
            {
                shits <- sh[qh==d]
                scores <- x[shits]$score
                is.consistent <- all(scores == scores[1])
                return(is.consistent)
            })
        return(res)
    }, simplify=FALSE)

## summarize individual absolute calls to population regions 
source("R/population_ranges.R")
