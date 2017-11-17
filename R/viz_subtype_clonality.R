############################################################
# 
# author: Ludwig Geistlinger and Levi Waldron
# date: 2017-08-12 22:31:59
# 
# descr: functionality for visualization of subtype clonality
# 
############################################################

# @args: gistic ... RangedSummarizedExperiment
circosSubtypeAssociation <- function(gistic)
{
    gistic <- as.data.frame(rowRanges(gistic))
    gistic$alterationTypeColor <- ifelse(gistic$type=="Deletion", "blue", "red")
    circlize::circos.initializeWithIdeogram(species="hg19", chr=paste0("chr",1:22), plotType=NULL)
    df <- circlize::read.chromInfo(species="hg19")
    circlize::circos.genomicInitialize(df$df[df$df[,1] %in% paste0("chr",1:22),], plotType=NULL)
    circlize::circos.trackPlotRegion(ylim=c(0, 1), 
        panel.fun =
            function(x, y)
            {   
                chr <- circlize::get.cell.meta.data("sector.index")
                xlim <- circlize::get.cell.meta.data("xlim")
                ylim <- circlize::get.cell.meta.data("ylim")
                circlize::circos.rect(xlim[1], 0, xlim[2], 0.5, col="white")
                circlize::circos.text(mean(xlim), 0.9, chr, 
                    cex=0.5, facing="clockwise", niceFacing=TRUE)       
            
                # add cnvrs
                ccnvrs <- gistic[gistic[,1]==chr,]
                if(nrow(ccnvrs))
                    for(i in seq_len(nrow(ccnvrs)))
                        circlize::circos.rect(ccnvrs[i,2], 0, 
                            ccnvrs[i,3], 0.5, 
                            col=ccnvrs[i,"alterationTypeColor"], 
                            border=ccnvrs[i,"alterationTypeColor"])
            }, bg.border = NA)
    circlize::circos.clear()
}

volcanoCorrelation <- function(rho, p)
{
    plot(x=rho, y=-log(p, base=10), 
        ylab="-log10(p)", xlab="Spearman correlation", col="white")
    text(x=rho, y=-log(p, base=10), names(allCT), cex=0.8)
}

plotSubclonalityDistributions <- function(subcl.scores)
{
    meds <- sapply(subcl.scores, median)
    ind <- order(meds)
    subcl.scores <- subcl.scores[ind]
    par(las=2)
    boxplot(subcl.scores, ylab="subclonality score")
}

plotNrSubtypesVsNrGisticRegions <- function(nr.sts, nr.gregs)
{
    plot(nr.sts, nr.gregs, col="white", xlab="#subtypes", ylab="#GISTIC.regions")
    text(nr.sts, nr.gregs, names(nr.sts), cex=0.8)
}

# plot nr of samples available for each cancer type
plotNrSamples <- function(gistic, subtypes, absolute)
{
    dat <- rbind(gistic, subtypes, absolute) 
    dat <- dat[,order(subtypes)]

    par(las=2)
    barplot(dat, beside=TRUE, ylab="#samples")
    legend("topleft", lwd=3, 
        legend=c("gistic", "subtype", "absolute"),  
        col=c("black", "darkgrey", "lightgrey"))
}

# scatterplot correlation between 
#   (1) subtype association score, and 
#   (2) subclonality score
plotCorrelation <- function(assoc.score, subcl.score)
{
    par(pch=20)
    plot(assoc.score, subcl.score, 
        xlab="subtype association score", ylab="subclonality score")
}
.extractUpper <- function(x)
{
    spl <- unlist(strsplit(x, ","))[2]
    spl <- substring(spl, 1, nchar(spl)-1)
    spl <- as.numeric(spl)
    return(spl)
}

plotPurityStrata <- function(subtys, puri.ploi)
{
    cids <- intersect(rownames(subtys), rownames(puri.ploi))    

    par(pch=20)
    plot(puri.ploi[cids,1], puri.ploi[cids,2], xlab="purity", ylab="ploidy")

    sebin <- stratifyByPurity(ovsubs, puri.ploi, method="equal.bin")
    sex <- sapply(names(sebin)[1:4], .extractUpper) 
    abline(v=sex, col="green", lty=3)

    squint <- stratifyByPurity(ovsubs, puri.ploi, method="quintile")
    sqx <- sapply(names(squint)[1:4], .extractUpper)
    abline(v=sqx, col="blue", lty=3)

    #legend("topright", lty=3, 
    #    legend=c("equal.bin", "quintile"), col=c("green", "blue"))
}



#@subtys: a matrix with sample IDs as rownames and at least a column 'cluster'
#@puri.ploi: a matrix with sample IDs as rownames and at least two columns
#               named 'purity' and 'ploidy'
#@returns: plots to a graphics device
plotSubtypePurityPloidy <- function(subtys, puri.ploi)
{
    cids <- intersect(rownames(subtys), rownames(puri.ploi)) 
    subtys <- subtys[cids,"cluster"] 
    puri.ploi <- puri.ploi[cids,]
    par(mfrow=c(2,2))
    for(i in c(1:2,4,3))
    {   
        title <- colnames(puri.ploi)[i]
        vlist <- split(puri.ploi[,i], subtys)
        if(i != 3) suppressWarnings( boxplot(vlist, 
                    ylab=title, varwidth=TRUE, notch=TRUE, main=title) )
        else
        {
            tlist <- sapply(vlist,table)
            barplot(tlist, ylab=title, main=title)
        }
    }
}

# GVIZ
#library(Gviz)
#gen <- "hg19"
#chr <- "1"
#pstart <- 1
##full
#pend <- 249250621
#
#col0n <- "steelblue"
#col1n <- "darkseagreen3"
#col3n <- "coral"
#col4n <- "indianred"
#colM <- "gray24"
#colB <- "gray94"
#
#ideoTrack <- IdeogramTrack(genome = gen, chromosome = chr)
#axisTrack <- GenomeAxisTrack(fontsize=15, littleTicks=F, from=start, to=end)
#
#covTrack <-   DataTrack(range=absCov, genome=gen, name="#Samples",  
#                        chromosome=chr, type="histogram", start=pstart, end=pend, 
#                        col = col0n,
#                        fill.histogram=col0n, 
#                        ylim=c(0,516), 
#                        cex.title=1, cex.axis=1, font.axis=2, 
#                        col.title=colM, col.axis=colM, background.title=colB)
#
#atracks <- sapply(1:100, function(i) 
#    AnnotationTrack(grlmatched[[i]], stacking="dense", overlapping=TRUE))
#
#plotTracks(c(ideoTrack, axisTrack, covTrack), from=pstart, to=pend, chromosome=chr)
