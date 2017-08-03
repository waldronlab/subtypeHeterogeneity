library(dplyr)
library(GenomicRanges)
library(MultiAssayExperiment)


# Load the subtypes from Ovcsubtypes repo - part of another validation project, but used here to get
# the sample IDs of patient and their subtypes (as classified by Verhaak, Helland or Konecny or Consensus classifier) 

load("../OvcSubtypes_Files/pooled.subtypes.updated.RData")
TCGA_subtypes <- pooled.subtypes[pooled.subtypes$data.source == "TCGA",]
rownames(TCGA_subtypes) <-
  gsub("\\.","-",
       substring(row.names(TCGA_subtypes),first= 11, last = 22))
TCGA_subtypes$individual_id <- row.names(TCGA_subtypes)

subtype.margin_purity <- merge(Pathology_ABSOLUTE,TCGA_subtypes, by = "individual_id")
subtype.margin_purity.filtered <- subtype.margin_purity[subtype.margin_purity$individual_id %in% row.names(TCGA_subtypes),]


subtype.margin_purity.filtered$Consensus.genepairs.rf <- factor(sub("_consensus", "",
                                                                    as.character(subtype.margin_purity.filtered$Consensus.genepairs.rf)))

subtype.margin_purity.filtered$Consensus.genepairs.rf <- factor(subtype.margin_purity.filtered$Consensus.genepairs.rf, levels = c("IMR", "DIF", "PRO", "MES"))

subtype.margin_purity.filtered$Consensus.genepairs.rf 

# Note other methods like Verhaak, Konecny or Helland may also be loaded instead of the Consensus classifier - part of another project  

### Processing ABSOLUTE pancan CNV file
#ABSOLUTE_hetro_scores <- read.table("/Users/lavanyakannan/Documents/ovc_purecn/ABSOLUTE/pancan12_absolute.segtab.txt",header = TRUE, sep="\t", as.is = TRUE)
ABSOLUTE_hetro_scores <- readr::read_tsv("/Users/lavanyakannan/Documents/ovc_purecn/ABSOLUTE/pancan12_absolute.segtab.txt")

### Processing the gistic Ranges (146 of them taken from Firebrowser)
cnv_lesions <- read.table("/Users/lavanyakannan/Documents/ovc_purecn/ABSOLUTE/firebrowseOV_all_lesions.conf_99.txt", 
                              header = TRUE, sep = "\t",
                              as.is=TRUE)


## Extract TCGA samples


names(cnv_lesions)
sample_names <- substr(names(cnv_lesions),start = 1, stop = 12)
names(cnv_lesions) <- gsub("\\.", "-", sample_names)

# table(gsub("\\.", "-", sample_names) %in% TCGA_subtypes$individual_id)

TCGA_individual_ids <- subtype.margin_purity.filtered$individual_id[!is.na(subtype.margin_purity.filtered$Consensus.genepairs.rf)]
tcga_samples_with_subtypes <- names(cnv_lesions)[names(cnv_lesions) %in% TCGA_individual_ids]
cnv_lesions.filtered <- cnv_lesions[,c("Peak-Limits",tcga_samples_with_subtypes)]
## Why are the CNV values in decimal points for rows 74 - 146?

### convert cnv_lesions.filtered into granges
n_gistic_ranges <- dim(cnv_lesions.filtered)[1]
split_lessions <- sapply(1:n_gistic_ranges, function (i) strsplit(cnv_lesions.filtered$"Peak-Limits"[i], ":"))
cnv_lesions.filtered$chr <- sapply(1:n_gistic_ranges, function (l) split_lessions[[l]][1])
split_lessions_start <- sapply(1:n_gistic_ranges, function (i) strsplit(split_lessions[[i]][[2]], "-"))
cnv_lesions.filtered$start <- sapply(1:n_gistic_ranges, function (l) split_lessions_start[[l]][1])
split_lessions_end <- sapply(1:n_gistic_ranges, function (i) strsplit(split_lessions_start[[i]][2], "[(]")) 
cnv_lesions.filtered$end <- sapply(1:n_gistic_ranges, function (l) split_lessions_end[[l]][1])

### get associations between granges (rows 1:73 for which the cnv are integer) and subtypes
subtypes <- subtype.margin_purity.filtered[match(names(cnv_lesions.filtered), TCGA_individual_ids),"Consensus.genepairs.rf"]
subtypes <- subtypes[!is.na(subtypes)]
p_values.with.gistic.cnvs <- lapply(1:73, function (i) chisq.test(cnv_lesions.filtered[i,-c(1,438, 437, 436)] , subtypes)$p.value) %>% unlist
hist(p_values.with.gistic.cnvs)

### create a GRanges object with gistic ranges 
gistic_granges <- makeGRangesFromDataFrame(df = cnv_lesions.filtered, seqinfo = NULL,keep.extra.columns = FALSE,seqnames.field = "chr", start.field = "start", end.field = "end")

##### fyhtyjbukttgfix the problematic ones 

makeDisjointWithScore <- function(x, FUN=mean){
  dj <- disjoin(x, with.revmap=TRUE)
  dj$score <- sapply(dj$revmap, function(i) FUN(x$score[i]))
  return( dj )
}


p_values_chi_sq_cnv_vs_subtype <- function (score_type = "Modal", restrict_to_gistic = TRUE) 
{

 
ABSOLUTE_hetro_scores$score <- 
  if (score_type == "Modal") (ABSOLUTE_hetro_scores$Modal_HSCN_1 + ABSOLUTE_hetro_scores$Modal_HSCN_2) else
    if (score_type == "subclonal") (ABSOLUTE_hetro_scores$Subclonal_HSCN_a1 + ABSOLUTE_hetro_scores$Subclonal_HSCN_a2 ) else
      if (score_type == "Modal_minus_subclonal") (ABSOLUTE_hetro_scores$Modal_HSCN_1 + ABSOLUTE_hetro_scores$Modal_HSCN_2) -
  (ABSOLUTE_hetro_scores$Subclonal_HSCN_a1 + ABSOLUTE_hetro_scores$Subclonal_HSCN_a2 )

  

ABSOLUTE_hetro_scores$score <- ABSOLUTE_hetro_scores$score %>% as.numeric
############# Extract only the Ovarian tumors whose subtypes are known 
load("/Users/lavanyakannan/Documents/ovc_purecn/OvcSubtypes_Files/pooled.subtypes.updated.RData")
TCGA_subtypes <- pooled.subtypes[pooled.subtypes$data.source == "TCGA",]
rownames(TCGA_subtypes) <-
  gsub("\\.","-",
       substring(row.names(TCGA_subtypes),first= 11, last = 22))
TCGA_subtypes$individual_id <- row.names(TCGA_subtypes)

ovarian_subtypes <- substr(x = ABSOLUTE_hetro_scores$Sample,start = 1, stop = 12) %in% rownames(TCGA_subtypes)

############# Convert ABSOLUTE CNV scores into grangeslist
## all unnecessary columns are removed in df1 and subsequent containers
df1 <- ABSOLUTE_hetro_scores[ovarian_subtypes,c("Sample", "Chromosome","Start", "End", "score")]
ABSOLUTE_GRanges_list <- makeGRangesListFromDataFrame(df = df1, split.field = "Sample", keep.extra.columns=TRUE)
#ABSOLUTE_GRanges_list <- endoapply(ABSOLUTE_GRanges_list, function (sample) {sample$Sample <- NULL
                                                   # sample})

##### find the GRanges that are disjoint

idj <- sapply(ABSOLUTE_GRanges_list, isDisjoint)
summary(idj)

idj[!idj]



#### reassemble the absolute_GRanges replacing the problematic ones by the fixed ones 
# # below is an attempt to only replace those Grangs that have intersecting segments - not working
# problems <- ABSOLUTE_GRanges_list[!idj]
# fixed <- endoapply(problems, makeDisjointWithScore)
# gr1 <- ABSOLUTE_GRanges_list[idj]
# gr2 <- endoapply(fixed, function (sample) {sample$revmap <- NULL
#                                            sample})
# ABSOLUTE_GRanges_list_fixed <- GRangesList(splitAsList(c(gr1,gr2), seq_along(gr1)))
# GRangesList(gr1,gr2)
# , seq_along(ABSOLUTE_GRanges_list[idj])))



ABSOLUTE_GRanges_list_fixed <- endoapply(ABSOLUTE_GRanges_list, makeDisjointWithScore)
  
############# form the Ranged Ragged Assay
ABSOLUTE_RRA <- RangedRaggedAssay(ABSOLUTE_GRanges_list_fixed)

############ return a matrix of copy number scores
djr <- disjoin(unlist(ABSOLUTE_RRA))
cnv_profile_matrix <- assay(ABSOLUTE_RRA, ranges=djr, background = 2)
# cnv_profile_matrix2 <- assay(ABSOLUTE_RRA, background = 2)
# below line is not working because assay method does not work well for independently supplied ranges
# should this be fixed now/ here? - needed if gistic ranges are to be used 
# cnv_profile_gistic_ranges <- assay(ABSOLUTE_RRA, ranges=gistic_granges, background = 2)



# View(cnv_profile_matrix)
# View(cnv_profile_gistic_ranges)

# change genomic ranges formatt of granges
gr <- unlist(ABSOLUTE_RRA)
newStyle <- mapSeqlevels(seqlevels(gr), "UCSC")
gr <- renameSeqlevels(gr, newStyle)



subtypes <- subtype.margin_purity.filtered[match(substring(colnames(cnv_profile_matrix), first = 1, last=12), TCGA_individual_ids),"Consensus.genepairs.rf"]
# subtypes <- subtypes[!is.na(subtypes)]
if (restrict_to_gistic){
mtch <- findOverlaps(gr, gistic_granges, type = "within")
p_values <- lapply(queryHits(mtch), function (i) chisq.test(cnv_profile_matrix[i,], subtypes)$p.value) %>% unlist
} else {
p_values <- lapply(1:(cnv_profile_matrix %>% dim), function (i) chisq.test(cnv_profile_matrix[i,], subtypes)$p.value) %>% unlist
}
p_values
}

piper
magrittr


modal.gistic_ranges <- p_values_chi_sq_cnv_vs_subtype(score_type = "Modal", restrict_to_gistic = TRUE)
subclonal.gistic_ranges <- p_values_chi_sq_cnv_vs_subtype(score_type = "subclonal", restrict_to_gistic = TRUE) 
modal_minus_subclonal.gistic_ranges <- p_values_chi_sq_cnv_vs_subtype(score_type = "Modal_minus_subclonal", restrict_to_gistic = TRUE)

modal <- p_values_chi_sq_cnv_vs_subtype(score_type = "Modal", restrict_to_gistic = FALSE)
subclonal <- p_values_chi_sq_cnv_vs_subtype(score_type = "subclonal", restrict_to_gistic = FALSE) 
modal_minus_subclonal <- p_values_chi_sq_cnv_vs_subtype(score_type = "Modal_minus_subclonal", restrict_to_gistic = FALSE)
par(mfrow = c(2,3))
hist(modal.gistic_ranges)
hist(subclonal.gistic_ranges)
hist(modal_minus_subclonal.gistic_ranges)
hist(modal)
hist(subclonal)
hist(modal_minus_subclonal)
