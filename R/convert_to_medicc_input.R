### this file will convert a copy number profile from a txt file to the 

library(dplyr)

## sample input file
cnv_lesions <- read.table("/Users/lavanyakannan/Documents/ovc_purecn/ABSOLUTE/OV/all_lesions.conf_95.OV.txt", 
                          header = TRUE, sep = "\t",
                          as.is=TRUE)






## create desc file
unlink("cnv_lesions_desc.txt")
sink(file="cnv_lesions_desc.txt")
cat("Peaks cnv_lesions_major.fasta cnv_lesions_minor.fasta")
sink()
file.show("cnv_lesions_desc.txt")



## create major allele
unlink("cnv_lesions_major.fasta")
sink(file="cnv_lesions_major.fasta")
for (i in names(cnv_lesions)[-c(1:9,573)]){
cat(paste(">", i,"\n",sep=''))
cat(paste0(cnv_lesions[1:66,i] %>% as.character,collapse =''))
cat("\n")
}
sink(type = "message")
sink()
file.show("cnv_lesions_major.fasta")

## create minor allele
minor_allele_dummy <- matrix(1, 66, length(names(cnv_lesions)[-c(1:9,573)]))
unlink("cnv_lesions_minor.fasta")
sink(file="cnv_lesions_minor.fasta")
for (i in 1:length(names(cnv_lesions)[-c(1:9,573)])){
  cat(paste(">", names(cnv_lesions)[-c(1:9,573)][i],"\n",sep=''))
  cat(paste0(minor_allele_dummy[,i] %>% as.character,collapse =''))
  cat("\n")
}
sink(type = "message")
sink()
file.show("cnv_lesions_minor.fasta")


