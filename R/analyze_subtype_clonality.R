############################################################
# 
# author: Ludwig Geistlinger
# date: 2017-08-22 17:34:32
# 
# descr: analyzing subclonality of expression subtypes
# 
############################################################

# @gistic: a RangedSummarizedExperiment
# @subtys: a matrix with sample IDs as rownames and at least a column 'cluster'
testSubtypes <- function(gistic, subtys, stat.only=FALSE, padj.method="BH")
{
    subtys <- subtys[rownames(subtys) %in% colnames(gistic), ]
    gistic <- gistic[,match(rownames(subtys), colnames(gistic))]
    subtys <- subtys$cluster    
    slot <- ifelse(stat.only, "statistic", "p.value")   
 
    res <- apply(assay(gistic), 1, 
        function(x) chisq.test(x, subtys)[[slot]])
    if(!stat.only) res <- p.adjust(res, method=padj.method)       
    return(res)
}

#
# testing the significance of Spearman's correlation
# (see https://en.wikipedia.org/wiki/Spearman
#       %27s_rank_correlation_coefficient#Determining_significance) 
#

# permutation test
corPermTest <- function(gistic, subtys, subcl.score, nperm=1000)
{
    obs.stat <- testSubtypes(gistic, subtys, stat.only=TRUE) 
    obs.cor <- abs(cor(obs.stat, subcl.score, method="spearman"))    
    times.greater <- 0

    x <- replicate(nperm,
    { 
        ind <- sample(nrow(subtys))
        subs.perm <- subtys[ind,]
        rownames(subs.perm) <- rownames(subtys)
        perm.stat <- testSubtypes(gistic, subs.perm, stat.only=TRUE)
        perm.cor <- cor(perm.stat, subcl.score, method="spearman")
        is.greater <- abs(perm.cor) >= obs.cor
        times.greater <<- times.greater + is.greater
    })
    p <- (times.greater + 1) / (nperm + 1)
    return(p)
}
# access RaggedExperiment: assay(ra[1:5,1:5], "score")

## summarize absolute calls in gistic regions
.maxScore <- function(scores, ranges, qranges) max(scores, na.rm=TRUE)

.wmean <- function(scores, ranges, qranges)
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
    sum.method <- ifelse(sum.method == "wmean", .wmean, .maxScore)
    qa <- RaggedExperiment::qreduceAssay(ra, query, 
                simplifyReduce=sum.method, i="score", background=NA)
    return(qa)
}

stratifyByPurity <- function(subtys, puri.ploi, method=c("quintile", "equal.bin"))
{
    cids <- intersect(rownames(subtys), rownames(puri.ploi))
    
    method <- match.arg(method)
    if(method == "quintile")
        cps <- quantile(puri.ploi[cids,1], seq(0, 1, 0.2), na.rm=TRUE) 
    else cps <- 5

    n <- levels(cut(puri.ploi[cids,1], cps))
    bins <- cut(puri.ploi[cids,1], cps, labels=FALSE) 
    spl <- split(cids, bins)
    names(spl) <- n
    return(spl)
}

analyzeStrata <- function(ids, absGRL, gistic, subtys)
{
    sabs <- absGRL[intersect(ids, names(absGRL))]
    sgis <- gistic[,intersect(ids, colnames(gistic))]
    ssub <- subtys[intersect(ids, rownames(subtys)),]
    
    ra <- getMatchedAbsoluteCalls(sabs, ssub)
    subcl <- querySubclonality(ra, query=rowRanges(sgis))
    subcl.score <- rowMeans(subcl, na.rm=TRUE)
    assoc.score <- suppressWarnings(
                        testSubtypes(sgis, ssub, stat.only=TRUE))
    rho <- cor(assoc.score, subcl.score, method="spearman")
    
    return(rho) 
}

analyzeCancerType <- function(ctype, absGRL)
{
    suppressWarnings({
    gistic <- gistic2RSE(ctype)
    subs <- getBroadSubtypes(ctype)
    ra <- getMatchedAbsoluteCalls(absGRL, subs)
    subcl <- querySubclonality(ra, query=rowRanges(gistic))
    subcl.score <- rowMeans(subcl, na.rm=TRUE)
    assoc.score <- testSubtypes(gistic, subs, stat.only=TRUE)
    rho <- cor(assoc.score, subcl.score, method="spearman")
    p <- cor.test(assoc.score, subcl.score, method="spearman", exact=FALSE)$p.value
    return(list(
        rho=rho, p=p, 
        assoc.score=assoc.score, 
        subcl.score=subcl.score,
        gistic=gistic, subtypes=subs
        ))
    
    })
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
