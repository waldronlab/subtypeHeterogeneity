############################################################
# 
# author: Ludwig Geistlinger
# date: 2017-08-22 17:34:32
# 
# descr: analyzing subclonality of expression subtypes
# 
############################################################

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

# @gistic: a RangedSummarizedExperiment
# @subtys: a matrix with sample IDs as rownames and at least a column 'cluster'
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


perm.test <- function()
replicate(1:1000,
    { 
        ind <- sample(nrow(ovsubs))
        ovsubs.perm <- ovsubs[ind,]
        rownames(ovsubs.perm) <- rownames(ovsubs)
        stat
    })

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
