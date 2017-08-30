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

# find gistic regions that are overlapped by more than one absolute call
getAmbigiousCalls <- function(grl, query)
{
    g <- sapply(names(grl), 
        function(n)
        { 
            message(match(n, names(grl)))
            x <- grl[[n]]
            olaps <- findOverlaps(query, x)
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
}


## summarize individual absolute calls to population regions 
# source("R/population_ranges.R")
