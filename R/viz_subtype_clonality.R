############################################################
# 
# author: Ludwig Geistlinger and Levi Waldron
# date: 2017-08-12 22:31:59
# 
# descr: functionality for visualization of subtype clonality
# 
############################################################

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
