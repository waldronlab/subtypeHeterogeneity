scriptdir = "/home/lwaldron/git/ovc_purecn"
data.dir = "/media/sphadmin/16tb/ovcnormals"

java7.exec = "/usr/local/share/jre1.7.0_79/bin/java"
mutect.jar = "/usr/local/share/mutect/1.1.7/bin/mutect-1.1.7.jar"
gatk.jar = "/usr/local/share/gatk/3.6/GenomeAnalysisTK.jar"

REFERENCE="~/mutect/hg38bundle/Homo_sapiens_assembly38.fasta"
DBSNP_VCF="~/mutect/hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz"
COSMIC_VCF="~/mutect/CosmicCodingMuts.vcf.gz"
COVERAGE_FILE="~/mutect/SeqCapEZ_Exome_v3_gcgene_file.txt"

out.dir = file.path(scriptdir, "output")
manifest.file = file.path(scriptdir, "output/ovc_manifest.tsv")
bed.dir = file.path(scriptdir, "data")

manifest = read.delim(manifest.file, as.is=TRUE)

manifest <- manifest[!is.na(manifest$bedfiles), ]

allcalls <- lapply(seq_along(manifest$id), function(i){
  uuid <- manifest[i, "id"]
  fname <- manifest[i, "filename"]
  fullname <- file.path(data.dir, uuid, fname)
  bed <- manifest[i, "bedfiles"]
  patient_id <- manifest[i, "patient_id"]
  mycall <- paste(java7.exec,
                  "-Xmx6g -jar", 
                  mutect.jar, 
                  "--analysis_type MuTect",
                  "-R", REFERENCE,
                  "--artifact_detection_mode", 
                  "--dbsnp", DBSNP_VCF,
                  "--cosmic", COSMIC_VCF,
                  "-dt None",
                  "-I:tumor", fullname,
                  "-o", file.path(out.dir, paste0(patient_id, "_artifact_stats.txt")),
                  "-vcf", file.path(out.dir, paste0(patient_id, "_artifact_mutect.vcf")))
})


allcalls <- allcalls[!file.exists(file.path(out.dir, paste0(manifest[,"patient_id"], "_artifact_mutect.vcf")))]

if(length(allcalls) > 0){
  library(BiocParallel)
  res <- bplapply(allcalls, function(mycall) system(mycall), BPPARAM = MulticoreParam(12))
}

## Merge normal VCFs into one:
merged.file = file.path(out.dir, "normals.merged.min2.vcf")
merge.cmd = paste("java -Xmx6g -jar", gatk.jar, "-T CombineVariants --minimumN 2 --genotypemergeoption UNSORTED -R", REFERENCE)
merge.cmd = paste(merge.cmd, 
                  paste("--variant", dir(out.dir, pattern="^.*_artifact_mutect\\.vcf$", full.names = TRUE), collapse=" "))
merge.cmd = paste(merge.cmd, ">", merged.file)
print(merge.cmd)
system(merge.cmd)

# for VCF in $OUT/*bwa_mutect_artifact_detection_mode.vcf;
# do
# CMD="$CMD --variant $VCF"
# done
# CMD="$CMD -o $OUT/normals.merged.min2.vcf"
# echo $CMD > $OUT/merge_normals_min2.sh

