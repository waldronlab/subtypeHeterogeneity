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
testSubtypes <- function(gistic, subtys, 
    what=c("p.value", "statistic", "subtype", "full"), padj.method="BH")
{
    subtys <- subtys[rownames(subtys) %in% colnames(gistic), ]
    gistic <- gistic[,match(rownames(subtys), colnames(gistic))]
    subtys <- subtys$cluster  
    
    what <- match.arg(what)  
    res <- apply(assay(gistic), 1, function(x) chisq.test(x, subtys))
    if(what == "full") return(res)

    if(what == "subtype")
    {
        res <- sapply(res,
            function(x)
            {
                diff <- x$observed - x$expected
                csums <- colSums(diff^2 / x$expected)
                return(which.max(csums))
            }, USE.NAMES=FALSE)
    }
    else 
    {
        res <- sapply(res, function(x) x[[what]])
        if(what == "p.value") res <- p.adjust(res, method=padj.method)
    } 
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
    obs.stat <- testSubtypes(gistic, subtys, what="statistic") 
    obs.cor <- abs(cor(obs.stat, subcl.score, method="spearman"))    
    times.greater <- 0

    x <- replicate(nperm,
    { 
        ind <- sample(nrow(subtys))
        subs.perm <- subtys[ind,]
        rownames(subs.perm) <- rownames(subtys)
        perm.stat <- testSubtypes(gistic, subs.perm, what="statistic")
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
querySubclonality <- function(ra, query, 
    sum.method=c("any", "wmean"), ext.range=0)
{
    sum.method <- match.arg(sum.method)
    sum.method <- ifelse(sum.method == "wmean", .wmean, .maxScore)

    if(ext.range)
    {
        start(query) <- sapply(start(query),
            function(s) ifelse(s < ext.range, 1, s-ext.range))
        end(query) <- sapply(end(query), function(s) s + ext.range)
    }

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
    subcl <- querySubclonality(ra, query=rowRanges(sgis), ext=500000)
    subcl.score <- rowMeans(subcl, na.rm=TRUE)
    assoc.score <- suppressWarnings(
                        testSubtypes(sgis, ssub, what="statistic"))
    rho <- cor(assoc.score, subcl.score, method="spearman")
    
    return(rho) 
}

analyzeCooccurence <- function(gistic)
{
    coocc <- matrix(0, nrow=nrow(gistic), ncol=nrow(gistic))
    for(i in seq_len(nrow(gistic)))
    {
        curri <- assay(gistic)[i,]
        for(j in seq_len(nrow(gistic)))
        {
            currj <- assay(gistic)[j,]
            coocc[i,j] <- sum((curri > 0) & (currj > 0))
        }
    }
    return(coocc)
}

stratifySubclonality <- function(gistic, subcl, subtys, tests=NULL, 
    st.names=names(stcols))
{
    isamples <- intersect(colnames(subcl), colnames(gistic))
    gistic <- gistic[,isamples]
    subcl <- subcl[,isamples]
    subtys <- subtys[colnames(subcl), "cluster"]
    .subclFract <- function(i)
    {
        # split by subtype
        subclspl <- split(subcl[i,], subtys)
        gisticspl <- split(assay(gistic)[i,], subtys)
        
        # split by copy number
        .subclFractCN <- function(j)
        { 
            spl <- split(subclspl[[j]], gisticspl[[j]])
            f2c <- vapply(spl, mean, numeric(1), na.rm=TRUE)
            return(f2c)
        }

        # subcl fraction per subtype
        f2s <- lapply(seq_along(st.names), .subclFractCN)
        
        # check length (sometimes one CN state is not present for one subtype)
        ls <- lengths(f2s)
        uls <- unique(ls)
        if(length(uls) > 1)
        {
            m <- which.max(ls)
            n <- names(f2s[[m]])
            for(k in seq_along(st.names))
            { 
                fk <- f2s[[k]]
                lfk <- length(fk)
                if(lfk < ls[m]) 
                {
                    ind <- setdiff(n, names(fk))  
                    f2s[[k]] <- c(fk, rep(0,length(ind)))
                    fk <- f2s[[k]]
                    lfk <- length(fk)
                    names(f2s[[k]])[(lfk-length(ind)+1):lfk] <- ind
                    f2s[[k]] <- f2s[[k]][n]
                }
            }
        } 
        
        f2s <- do.call(cbind, f2s)
        return(f2s)
    }
    stratL <- lapply(seq_len(nrow(subcl)), .subclFract)
    if(is.null(tests)) return(stratL)

    # extrapolate?
    .extrapolateSubcl <- function(i)
    {
        ot <- tests[[i]]$observed
        rt <- round(ot * stratL[[i]])
        fl <- rbind(ot[2,] - rt[2,], rt[2,])
        if(nrow(ot) > 2) fl <- rbind(fl, ot[3,] - rt[3,], rt[3,])
        ot <- rbind(ot[1,], fl)
        colnames(ot) <- st.names
        return(ot)
    }
    extraL <- lapply(seq_len(nrow(subcl)), .extrapolateSubcl)
    return(extraL)
}

analyzeCancerType <- function(ctype, absGRL)
{
    suppressWarnings({
    gistic <- gistic2RSE(ctype)
    subs <- getBroadSubtypes(ctype)
    ra <- getMatchedAbsoluteCalls(absGRL, subs)
    subcl <- querySubclonality(ra, query=rowRanges(gistic))
    subcl.score <- rowMeans(subcl, na.rm=TRUE)
    assoc.score <- testSubtypes(gistic, subs, what="statistic")
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
