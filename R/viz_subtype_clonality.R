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
cb.pink <- "#CC79A7"
cb.darkred <- "#B42F32"
cb.red <- "#D55E00"
cb.lightred <- "#DF6747"
cb.blue <- "#0072B2"
cb.yellow <- "#F0E442"
cb.green <- "#009E73"
cb.lightblue <- "#56B4E9"
cb.lightorange <- "#FAAC77"
cb.orange <- "#E69F00"
cb.darkorange <- "#F6893D"
cb.lightgrey <- "#C9C9BD"
cb.darkgrey <- "#878D92"


stcols <- c(cb.lightblue, cb.green, cb.orange, cb.pink) 
bw.stcols <- grey(c(0.8, 0.6, 0.4, 0.2)) 
SUBTYPES <- c("DIF", "IMR", "MES", "PRO")
names(stcols) <- names(bw.stcols) <- SUBTYPES 

# as from Figure 1b in TCGA OVC paper, Nature 2011
getCnvGenesFromTCGA <- function(excl = TRUE, build = c("hg19", "hg38"))
{
    build <- match.arg(build)
    data.dir <- system.file("extdata", package="subtypeHeterogeneity")
    cnv.genes.file <- file.path(data.dir, "ovc_cnv_genes.txt") 
    cnv.genes <- scan(cnv.genes.file, what="character")

    if(excl) excl.genes <- c("XPR1", "PAX8", "SKP2", "PRIM2", "SOX17", "ERBB2",
                             "ERBB3", "DEPTOR", "CDKN2A", "ZMYND8", "ANKRD11")
    else excl.genes <- c("XPR1", "CD47", "TACC3", "DEPTOR", "ERBB2", "ERBB3",
                         "METTL17", "ANKRD11", "MAP2K4", "ZMYND8", "SC5D")
    cnv.genes <- setdiff(cnv.genes, excl.genes)

    if(build == "hg19") edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    else
    {
        ah <- AnnotationHub::AnnotationHub() 
        # query(ah, "EnsDb for Homo sapiens")
        edb <- ah[["AH78783"]] 
    }
    filtr <- AnnotationFilter::GeneNameFilter(cnv.genes)
    gr <- GenomicFeatures::genes(edb, filter = list(filtr))
    gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
    gr <- gr[gr$gene_biotype == "protein_coding"]
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
circosSubtypeAssociation <- function(gistic, cnv.genes, bw = FALSE)
{
    gistic <- as.data.frame(rowRanges(gistic))
    gistic$alterationTypeColor <- ifelse(gistic$type=="Deletion", cb.blue, cb.red)
    if(bw)
    {
        gistic$alterationTypeColor <- ifelse(gistic$type=="Deletion", gray(0.8), gray(0.1))
        stcols <- rev(gray(c(0.1, 0.3, 0.6, 0.9)))
        names(stcols) <- c("DIF", "IMR", "MES", "PRO")
    }
    gistic$subtype <- c("PRO", "MES", "DIF", "IMR")[gistic$subtype]
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
        col=c(gray(0.8), gray(0.1)), lwd=2, cex=0.6, title="Outer circle")

    legend("bottomright", legend = names(stcols), 
        col=stcols, lwd=2, cex=0.6, title="Inner circle")
}

gvizRegion <- function(r, g, ov.abscalls, sarc.abscalls=NULL, 
    ov.rel=NULL, sarc.rel=NULL, type=c("+", "-"), pcn.calls=NULL, weight=FALSE)
{
    colM <- "gray24"
    colB <- "gray94"

    genome <- GenomeInfoDb::genome(r)
    chr <- as.character(GenomicRanges::seqnames(r))
    pstart <- GenomicRanges::start(r) - 30000
    pend <- GenomicRanges::end(r) + 30000

    # Ideogram & Genome axis
    ideoTrack <- Gviz::IdeogramTrack(genome=genome, chromosome=chr, fontsize=15)
    axisTrack <- Gviz::GenomeAxisTrack(fontsize=15, littleTicks=FALSE, 
                                        from=pstart, 
                                        to=pend)
    # CNV region
    cnvrTrack <- Gviz::AnnotationTrack(r, name="GISTIC2", fill=cb.red, 
        col.title=colM, col.axis=colM, background.title=colB, cex.title=0.8)

    # Genes contained in the CNV region
    geneTrack <- Gviz::AnnotationTrack(g, name="", group=names(g), just.group="right",
            col.title=colM, col.axis=colM, background.title=colB, cex.title=0.7)
        

    # OV: Coverage of calls in the CNV region
    rel <- table(ovsubs$cluster)
    type <- match.arg(type)
    totCN <- ov.abscalls$Modal_HSCN_1 + ov.abscalls$Modal_HSCN_2
    ind <- if(type == "+") totCN > 2 else totCN < 2
    ov.abscalls <- ov.abscalls[ind]
    cl.olaps <- ov.abscalls[ov.abscalls$score == 0]
    subcl.olaps <- ov.abscalls[ov.abscalls$score == 1]
    cl.covTrack <- .constrCovTrack(r, cl.olaps, 
        trans=TRUE, weight=weight, relative=ov.rel, title="OV clonal")
    subcl.covTrack <- .constrCovTrack(r, subcl.olaps, 
        trans=TRUE, weight=weight, relative=ov.rel, title="OV subclonal")
    tracks <- c(ideoTrack, axisTrack, cnvrTrack, geneTrack, cl.covTrack, subcl.covTrack) 

    # SARC: Coverage of calls in the CNV region
    if(!is.null(sarc.abscalls))
    {
        cl.olaps <- sarc.abscalls[sarc.abscalls$score == 0]
        subcl.olaps <- sarc.abscalls[sarc.abscalls$score == 1]
        scl.covTrack <- .constrCovTrack(r, cl.olaps, 
            trans=TRUE, weight=weight, relative=sarc.rel, title="SARC clonal")
        ssubcl.covTrack <- .constrCovTrack(r, subcl.olaps, 
            trans=TRUE, weight=weight, relative=sarc.rel, title="SARC subclonal")
        tracks <- c(tracks, scl.covTrack, ssubcl.covTrack)
    } 
    # Put it all together and plot it
    Gviz::plotTracks(tracks, from=pstart, to=pend,
            groupAnnotation="group", chromosome=chr, littleTicks=TRUE)
    
}

.constrCovTrack <- function(r, calls, type="l", weight=FALSE, relative=NULL,
    delimit=TRUE, trans=FALSE, minCov=0, stateOrder=NULL, title="#samples")
{
    colM <- "gray24"
    colB <- "gray94"

    gen <- GenomeInfoDb::genome(r)
    chr <- as.character(GenomicRanges::seqnames(r))
    gstart <- GenomicRanges::start(r)
    gend <- GenomicRanges::end(r)
    pstart <- gstart - 30000
    pend <- gend + 30000

    .getFreq <- function(gcalls, state)
    {   
        scalls <- subset(gcalls, State == state)
        w <- 1L
        if(weight) w <- scalls$Modal_HSCN_1 + scalls$Modal_HSCN_2
        x <- GenomicRanges::coverage(scalls, weight=w)[[chr]]
        cs <- cumsum(S4Vectors::runLength(x))
        len <- length(cs)
        grid <- seq_len(len-1)
        starts <-  cs[grid] + 1
        ends <- cs[grid + 1]
        cov <- S4Vectors::runValue(x)[grid + 1]
        df <- data.frame(seqnames=chr, start=starts, end=ends)
        gr <- GenomicRanges::makeGRangesFromDataFrame(df)
        gr$freq <- as.integer(cov)
        if(!is.null(relative)) 
            gr$freq <- round(gr$freq / relative[state] * 100, digits=2)
        gr <- gr[gr$freq >= minCov]
        if(delimit) gr <- GenomicRanges::restrict(gr,
                            start=GenomicRanges::start(r), 
                            end=GenomicRanges::end(r))
        return(gr)
    }

    ustates <- sort(unique(calls$State))
    if(!is.null(stateOrder)) ustates <- ustates[stateOrder]
    freqs <- lapply(ustates, function(s) .getFreq(calls, s))
    ind <- lengths(freqs) > 0
    ustates <- ustates[ind]
    freqs <- freqs[ind]
    maxSamples <- vapply(freqs, function(f) max(f$freq), numeric(1))
    maxSamples <- max(maxSamples)

    covTracks <- list()
    covTracks[[1]] <-   Gviz::DataTrack(range=freqs[[1]], genome=gen, name=title,
                        chromosome=chr, type=type, start=pstart, end=pend, 
                        col=stcols[1],
                        fill.histogram=stcols[1],
                        ylim=c(0, maxSamples),
                        cex.title=1, cex.axis=1, font.axis=2,
                        col.title=colM, col.axis=colM, background.title=colB)

    if(length(ustates) > 1)
    {
        grid <- seq_len(length(ustates) - 1) + 1
        for(i in grid)
        {
            covTracks[[i]] <- Gviz::DataTrack(range=freqs[[i]], genome=gen,
                        chromosome=chr, start=pstart, end=pend, 
                        type=type, 
                        col=stcols[i],
                        fill.histogram=stcols[i],
                        ylim=c(0, maxSamples),
                        col.title=colM, col.axis=colM, background.title=colB)
        }
    }

    if(trans) 
        for(i in seq_along(ustates)) 
            Gviz::displayPars(covTracks[[i]]) <- list(alpha = 0.75)
    
    covTrack <- Gviz::OverlayTrack(covTracks, background.title=colB)
    return(covTrack)
}

bpSubtypeStrata <- function(x, type=c("gain", "loss"))
{
    type <- match.arg(type)
    
    col <- .getStateColors(type) 
    ltext <- .getStateLegend(type)
    dens <- c(-1,-1, 40, -1, 40)

    if(nrow(x) < 5)
    {
        col <- col[1:3]
        ltext <- ltext[1:3]
        dens <- dens[1:3]
    }

    barplot(x, ylab="#tumors",
        col=col, 
        legend.text=ltext, 
        args.legend=c(x=ifelse(ncol(x) == 3, "topright", "top")), 
        density=dens, 
        border=NA)
}

ggplotSubtypeStrata <- function(x, type=c("gain", "loss"))
{
    .single <- function(y, ty=type)
    {
        ty <- match.arg(ty)
        col <- .getStateColors(ty) 
        ltext <- .getStateLegend(ty)
        alpha <- c(1,1,0.5,1,0.5)   
 
        if(nrow(y) < 5)
        {
            col <- col[1:3]
            ltext <- ltext[1:3]
            alpha <- alpha[1:3]
        }

        df <- reshape2::melt(y)
        colnames(df) <- c("state", "subtype", "tumors")
        df$state <- factor(ltext[df$state], levels=rev(ltext))
        df$subtype <- factor(as.vector(df$subtype), levels=sort(levels(df$subtype)))
        df$color <- rev(col)[df$state]
        df$alpha <- ifelse(grepl("subclonal", df$state), 0.5, 1)

        return(df)
    }
    
    if(is.list(x))
    {
        genes <- names(x)
        df <- lapply(seq_along(x), function(i) .single(x[[i]], type[i]))
        nr <- vapply(df, nrow, integer(1))
        df <- do.call(rbind, df)
        df$state <- factor(as.vector(df$state), levels=sort(levels(df$state))[c(2,1,4,3,6,5,7)])
        df$gene <- factor(rep(genes, times=nr), levels=genes)
    }
    else df <- .single(x, type)

    p <- ggplot2::ggplot(data=df, ggplot2::aes(x=subtype, y=tumors, fill=state, alpha=state)) + 
                ggplot2::geom_col() + 
                ggplot2::scale_fill_manual(values=rev(df$color[!duplicated(df$state)])) + 
                ggplot2::scale_alpha_manual(values=rev(df$alpha[!duplicated(df$state)])) 

    if(is.list(x)) p + ggplot2::facet_wrap( ~ gene, 
                       ncol = ifelse(length(x) %% 2 == 0, length(x) / 2, length(x)))
    else p
}

.getStateColors <- function(type=c("gain", "loss"))
{
    type <- match.arg(type)
    if(type == "gain")
    {
        n1.col <- cb.orange
        n2.col <- cb.red
    }
    else
    {
        n1.col <- cb.green
        n2.col <- cb.blue
    }
    col <- c("grey", rep(n1.col, 2), rep(n2.col, 2))
    return(col)
}

.getStateLegend <- function(type=c("gain", "loss"))
{
    type <- match.arg(type)
    ltext <- c("normal", 
                paste0("1-copy xxxx (", c("", "sub"), "clonal)"),
                paste0("2-copy xxxx (", c("", "sub"), "clonal)"))
    ltext <- sub("xxxx", type, ltext)
    if(type == "gain")
    {
        #ltext[4] <- substitute(pre>=suff, list(pre="", suff=ltext[4]))
        #ltext[5] <- substitute(pre>=suff, list(pre="", suff=ltext[5])) 
        ltext[4:5] <- paste0(">=", ltext[4:5])
    }
    return(ltext)
}


plotlyPie <- function(labels, values, colors, out.file=NULL)
{
    require(dplyr)
    p <- plotly::plot_ly(
        labels = labels, 
        values = values, 
        type = 'pie',
        textposition = 'inside',
        textinfo = 'label+percent',
        insidetextfont = list(color = '#FFFFFF', size=30),
        pull = 0.01,
        hole = 0.01, 
        marker = list(colors = colors,
                     line = list(color = '#FFFFFF', width = 1)),
       showlegend = FALSE) %>% 
 
    plotly::layout(
        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

    if(!is.null(out.file))
    {
        Sys.setenv('MAPBOX_TOKEN' = "pk.eyJ1IjoibHVkd2lnZyIsImEiOiJjamx3cXFwMm8xOGJ3M2tvZGV5amozNG5rIn0.-EpukkiU2QxcbUvFtq1QOw")
        plotly::orca(p, "pie-plot.pdf")
    }
    return(p)
}


plotCommonSCNAs <- function(ovranges, sarcranges)
{
    # restrict to significant subtye assocation
    ovsig <- ovranges[ovranges$adj.pval < 0.1]
    sarcsig <- sarcranges[sarcranges$adj.pval < 0.1]
    
    hits <- GenomicRanges::findOverlaps(ovsig, sarcsig)
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    ovsig <- ovsig[qh]
    sarcsig <- sarcsig[sh]
    
    is.type <- ovsig$type == sarcsig$type
    len <- length(hits)
    col <- rep("grey", len)
    for(i in seq_len(len)) 
        if(is.type[i]) col[i] <- ifelse(ovsig$type[i] == "Deletion", cb.blue, cb.red)

    par(pch=20)
    par(cex=1.1)
    plot(ovsig$subcl, sarcsig$subcl, 
            xlab="OV subclonality", ylab="SARC subclonality", col=col)
    text(x=0.3, y=0.249, "chr13:q14.2 (RB1)", pos=4)
    text(x=0.56, y=0.625, "chr20:q13.33 (EEF1A2)", pos=2)
    text(x=0.6, y=0.678, "chr8:q24.21 (MYC)", pos=2)
    text(x=0.463, y=0.351, "chr3:q13.31 (TUSC7)", pos=4)
    text(x=0.483, y=0.455, "chr2:q37.3 (ING5)", pos=4) 
    text(x=0.483, y=0.512, "chr2:q37.3 (TWIST2)", pos=4) 
    text(x=0.305, y=0.625, "chr8:p23.2 (CSMD1)", pos=4)
    text(x=0.265, y=0.58, "chr15:q15.1 (MGA)", pos=4)
    text(x=0.29, y=0.564, "chr15:q11.2 (Prader-Willi)", pos=4)
    text(x=0.47, y=0.5, "chr1:q42.3 (ARID4B)", pos=2)
    text(x=0.477, y=0.53, "chr6:p22.3 (RNF144B)", pos=4)

    legend("bottomright", lwd=2, col=c(cb.blue, cb.red, "grey"), legend=c("deletion", "amplification", "amp-OV/del-SARC"))
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
plotCorrelation <- function(assoc.score, subcl.score, 
    subtypes=NULL, stcols=NULL, xlim=NULL, lpos="bottomright")
{
    if(is.null(subtypes))
    { 
        col <- cb.red
        pch <- 20
    }
    else
    { 
        col <- stcols[subtypes]
        pch <- c(20, 15, 17, 8)
        pch <- pch[subtypes]
    }
    if(is.null(xlim)) xlim <- c(0, max(assoc.score))
    plot(assoc.score, subcl.score, col=col, xlim=xlim, pch = pch,
        xlab=expression(paste("Subtype association score ", italic(S[A]))), 
        ylab=expression(paste("Subclonality score ", italic(S[C]))))
    
    rho <- cor(assoc.score, subcl.score, method="spearman")
    rho <- paste("cor", round(rho, digits=3), sep=" = ")

    p <- cor.test(assoc.score, subcl.score, method="spearman", exact=FALSE)
    p <- p$p.value
    p <- paste("   p", round(p, digits=3), sep=" = ")
    
    legend("topright", legend=c(rho, p))
    if(!is.null(subtypes))
    {
        ind <- order(names(stcols))
        stcols <- stcols[ind]
        pch <- c(20, 15, 17, 8)
        legend(lpos, legend=names(stcols), col=stcols, pch = pch[ind])
    }

    abline(lm(subcl.score ~ assoc.score), lty=2, col="grey")
}

# barplot summarizing associated CNAs per subtype
plotNrCNAsPerSubtype <- function(type, subtype, bw = FALSE)
{
    spl <- split(type, subtype)
    .countType <- 
        function(s)
        {
            nr.amp <- sum(s == "Amplification")
            nr.del <- length(s) - nr.amp
            res <- c(nr.amp, nr.del)
            names(res) <- c("Amplification", "Deletion")
            return(res)
        }
    m <- vapply(spl, .countType, integer(2))
    m <- m[2:1,]     
    m <- m[,order(colSums(m))]
    colnames(m) <- c("MES", "DIF", "IMR", "PRO")

    col <- if(bw) gray(c(0.8, 0.1)) else c(cb.blue, cb.red)  
    if(bw)
    {
        stcols <- rev(gray(c(0.1, 0.3, 0.6, 0.9)))
        names(stcols) <- c("DIF", "IMR", "MES", "PRO")
    }
    bcol <- stcols[colnames(m)]
 
    bp <- barplot(m, col=col, border=NA, ylab="#CNAs")
    for(i in 1:4) 
        rect( xleft=bp[i]-0.5, ybottom=0, 
                xright=bp[i]+0.5, ytop=sum(m[,i]), 
                border=bcol[i], lwd=3)
}

.extractUpper <- function(x)
{
    spl <- unlist(strsplit(x, ","))[2]
    spl <- substring(spl, 1, nchar(spl)-1)
    spl <- as.numeric(spl)
    return(spl)
}

.extractLower <- function(x)
{
    spl <- unlist(strsplit(x, ","))[1]
    spl <- substring(spl, 2, nchar(spl))
    spl <- as.numeric(spl)
    return(spl)
}

plotPurityStrata <- function(puri.ploi, subtys, gistic, absGRL,
    method=c("equal.bin", "quintile"))
{
    method <- match.arg(method)
    puri.ploi <- puri.ploi[!is.na(puri.ploi[,1]),]
    cids <- intersect(rownames(subtys), rownames(puri.ploi))    

    sebin <- stratifyByPurity(ovsubs, puri.ploi, method=method)
    sex <- vapply(names(sebin), .extractUpper, numeric(1))
    sel <- vapply(names(sebin), .extractLower, numeric(1))
    if(method=="equal.bin") sebin <- sebin[-1] 
    strata.cor <- vapply(sebin, function(ids) 
                    analyzeStrata(ids, absGRL, gistic, subtys), numeric(1))
    
    par(pch=20)
    hcols <- rev(grey(2:6/8))
    col.ind <- vapply(puri.ploi[cids,1], 
                    function(p) min(which(p <= sex)), integer(1)) 
    plot(puri.ploi[cids,1], puri.ploi[cids,2], 
            col=hcols[col.ind], xlab="purity", ylab="ploidy")
    
    abline(v=sex[1:4], col=hcols[2:5], lty=3)
    if(method=="equal.bin")
    {
        sex <- sex[-1]
        sel <- sel[-1]
    } 
    xm <- sel + (sex - sel) / 2
    text(x=xm, y=8, labels=round(strata.cor, digits=2), col=hcols, font=2)
    axis(3, labels=lengths(sebin), at=xm)
    mtext("#tumors", line=3)
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
            cols <- c("mistyrose2", "moccasin", "plum3")
            tlist <- sapply(vlist, table)
            barplot(tlist, ylab="#tumors", ylim=c(0,170), main="", col=cols)
            legend("topleft", legend=0:2, lwd=4, col=cols, horiz=TRUE)
            
        }
    }
}
