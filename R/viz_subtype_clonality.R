############################################################
# 
# author: Ludwig Geistlinger and Levi Waldron
# date: 2017-08-12 22:31:59
# 
# descr: functionality for visualization of subtype clonality
# 
############################################################

bang_wong_colors <-
    c(
        "#CC79A7",
        "#D55E00",
        "#0072B2",
        "#F0E442",
        "#009E73",
        "#56B4E9",
        "#E69F00",
        "#000000"
    )

stcols <- c("steelblue", "darkseagreen3", "coral", "firebrick")
names(stcols) <- c("PRO", "MES", "DIF", "IMR")

# as from Figure 1b in TCGA OVC paper, Nature 2011
getCnvGenesFromTCGA <- function()
{
    data.dir <- system.file("extdata", package="subtypeHeterogeneity")
    cnv.genes.file <- file.path(data.dir, "ovc_cnv_genes.txt") 
    cnv.genes <- scan(cnv.genes.file, what="character")

    excl.genes <- c("XPR1", "PAX8", "SKP2", "PRIM2", "SOX17", "ERBB2",
                        "ERBB3", "DEPTOR", "CDKN2A", "ZMYND8", "ANKRD11")
    cnv.genes <- setdiff(cnv.genes, excl.genes)

    edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    filtr <- AnnotationFilter::GenenameFilter(cnv.genes)
    gr <- GenomicFeatures::genes(edb, filter=list(filtr))
    names(gr) <- gr$symbol
    gr <- as.data.frame(gr)
    gr <- gr[,1:3]
    gr[,1] <- as.integer(as.vector(gr[,1]))
    gr <- gr[do.call(order, gr),]
    gr[,1] <- paste0("chr", gr[,1])
    gr <- makeGRangesFromDataFrame(gr)
    return(gr) 
}

# @args: gistic ... RangedSummarizedExperiment
circosSubtypeAssociation <- function(gistic, cnv.genes)
{
    stcols <- c("steelblue", "darkseagreen3", "coral", "firebrick")
    gistic <- as.data.frame(rowRanges(gistic))
    gistic$alterationTypeColor <- ifelse(gistic$type=="Deletion", "blue", "red")
    gistic$subtypeColor <- stcols[gistic$subtype]
    circlize::circos.initializeWithIdeogram(species="hg19", chr=paste0("chr",1:22), plotType=NULL)
    df <- circlize::read.chromInfo(species="hg19")
    circlize::circos.genomicInitialize(df$df[df$df[,1] %in% paste0("chr",1:22),], plotType=NULL)

    # outer circle: GISTIC alteration type (amplification / deletion)
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

    # inner circle: subtype association
    circlize::circos.trackPlotRegion(ylim=c(0,1),
        panel.fun=
            function(x, y)
            {
                # chromosomal layout
                chr <- circlize::get.cell.meta.data("sector.index")
                xlim <- circlize::get.cell.meta.data("xlim")
                ylim <- circlize::get.cell.meta.data("ylim")
                circlize::circos.rect(xlim[1], 0, xlim[2], 0.5, col="white")
                
                # add cnvrs
                ccnvrs <- gistic[gistic[,1]==chr,]
                if(nrow(ccnvrs))
                    for(i in seq_len(nrow(ccnvrs)))
                    {
                        circlize::circos.rect(ccnvrs[i,2], 0,  
                            ccnvrs[i,3], 0.5, 
                            col=ccnvrs[i,"subtypeColor"], 
                            border=ccnvrs[i,"subtypeColor"])
                        
                       circlize::circos.text(mean(unlist(ccnvrs[i,2:3])), 0.2, 
                            ccnvrs[i,"significance"], cex=0.7)
                         #   #, facing="clockwise", niceFacing=TRUE)
                    }

                # add genes
                cnv.genes <- as.data.frame(cnv.genes)
                cgenes <- cnv.genes[cnv.genes[,1]==chr,]
                if(nrow(cgenes))
                    for(i in seq_len(nrow(cgenes)))
                    {
                        circlize::circos.text(mean(unlist(cgenes[i,2:3])), 0.8, 
                            substring(rownames(cgenes)[i],1,5), cex=0.4 
                                , facing="clockwise", niceFacing=TRUE)
            
                        #circlize    
                    }

            }, bg.border=NA)
    circlize::circos.clear()

    legend("bottomleft", legend=c("deletion", "amplification"), 
        col=c("blue", "red"), lwd=2, cex=0.6, title="Outer circle")

    legend("bottomright", legend=c("PRO", "MES", "DIF", "IMR"), 
        col=stcols, lwd=2, cex=0.6, title="Inner circle")
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
    plot(assoc.score, subcl.score, col="firebrick",
        xlab=expression(paste("Subtype association score ", italic(S[A]))), 
        ylab=expression(paste("Subclonality score ", italic(S[C]))))
    
    rho <- cor(assoc.score, subcl.score, method="spearman")
    rho <- paste("cor", round(rho, digits=3), sep=" = ")

    p <- cor.test(assoc.score, subcl.score, method="spearman", exact=FALSE)
    p <- p$p.value
    p <- paste("   p", round(p, digits=3), sep=" = ")
    
    legend("topright", legend=c(rho, p))
}

# barplot summarizing associated CNAs per subtype
plotNrCNAsPerSubtype <- function(type, subtype)
{
    spl <- split(type, subtype)
    m <- vapply(spl, table, integer(2))
    m <- m[2:1,]     
    colnames(m) <- names(stcols)
    m <- m[,order(colSums(m))]

    bp <- barplot(m, col=c("blue", "red"), border=NA, ylab="#CNAs")
    for(i in 1:4) 
        rect( xleft=bp[i]-0.5, ybottom=0, 
                xright=bp[i]+0.5, ytop=sum(m[,i]), 
                border=stcols[c(2:4,1)][i], lwd=3)
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
    subtys <- subtys[cids, "subtype"] 
    puri.ploi <- puri.ploi[cids,]
    par(mfrow=c(2,2))
    par(pch=20)
    par(cex.axis=1.1)
    par(cex=1.1)
    par(cex.lab=1.1)
    for(i in c(1:2,4,3))
    {   
        title <- colnames(puri.ploi)[i]
        title <- gsub("\\.", " ", title)
        vlist <- split(puri.ploi[,i], subtys)
        if(i != 3) suppressWarnings( boxplot(vlist, col=stcols[names(vlist)], 
                    ylab=title, varwidth=TRUE, notch=TRUE, main="") )
        else
        {
            cols <- c("plum3", "mistyrose2", "moccasin")
            tlist <- sapply(vlist, table)
            barplot(tlist, ylab=title, ylim=c(0,170), main="", col=cols)
            legend("topleft", legend=0:2, lwd=4, col=cols, horiz=TRUE)
            
        }
    }
}
