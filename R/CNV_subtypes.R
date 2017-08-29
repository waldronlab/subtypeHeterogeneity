# Load the subtypes from Ovcsubtypes repo

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


Number of losses  Number of gains
Sample 1
Sample 2
...

##########


cnv_by_genes <- read.table("ABSOLUTE/OV/all_thresholded.by_genes.OV.txt.bz2", 
           header = TRUE, sep = "\t",
           as.is=TRUE)

cnv_lesions <- read.table("ABSOLUTE/OV/all_lesions.conf_95.OV.txt.bz2", 
                           header = TRUE, sep = "\t",
                           as.is=TRUE)

## Extract TCGA samples


names(cnv_lesions)
sample_names <- substr(names(cnv_lesions),start = 1, stop = 12)
names(cnv_lesions) <- gsub("\\.", "-", sample_names)


table(gsub("\\.", "-", sample_names) %in% TCGA_subtypes$individual_id)

TCGA_individual_ids <- subtype.margin_purity.filtered$individual_id[!is.na(subtype.margin_purity.filtered$Consensus.genepairs.rf)]
tcga_samples_with_subtypes <- names(cnv_lesions)[names(cnv_lesions) %in% TCGA_individual_ids]
cnv_lesions.filtered <- cnv_lesions[1:66,c("Peak-Limits",tcga_samples_with_subtypes)]
#cnv_lesions.filtered <- cnv_lesions.filtered[rowSums(abs(cnv_lesions.filtered)) >= 200,]
subtypes <- subtype.margin_purity.filtered[match(names(cnv_lesions.filtered), TCGA_individual_ids),"Consensus.genepairs.rf"]
subtypes <- subtypes[!is.na(subtypes)]
p_values <- lapply(1:66, function (i) chisq.test(cnv_lesions.filtered[i,-1], subtypes)$p.value) %>% unlist
hist(p_values)
as.matrix(cnv_by_genes[,c(-1,-2,-3)]) + 2

### convert cnv_lesions.filtered into granges
split_lessions <- sapply(1:66, function (i) strsplit(cnv_lesions.filtered$"Peak-Limits"[i], ":"))
cnv_lesions.filtered$chr <- sapply(1:66, function (l) split_lessions[[l]][1])
split_lessions_start <- sapply(1:66, function (i) strsplit(split_lessions[[i]][[2]], "-"))
cnv_lesions.filtered$start <- sapply(1:66, function (l) split_lessions_start[[l]][1])
split_lessions_end <- sapply(1:66, function (i) strsplit(split_lessions_start[[i]][2], "[(]")) 
cnv_lesions.filtered$end <- sapply(1:66, function (l) split_lessions_end[[l]][1])
  
library(GenomicRanges)
gistic_granges <- makeGRangesFromDataFrame(df = cnv_lesions.filtered, seqinfo = NULL,keep.extra.columns = FALSE,seqnames.field = "chr", start.field = "start", end.field = "end")
