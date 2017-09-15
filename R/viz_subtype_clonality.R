############################################################
# 
# author: Ludwig Geistlinger and Levi Waldron
# date: 2017-08-12 22:31:59
# 
# descr: functionality for visualization of subtype clonality
# 
############################################################

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
