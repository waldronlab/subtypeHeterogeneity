# this script is to analyze the demographic data stratified by the methods and the subtypes; any interesting things about this data for subtypes may be reported in the upcoming JNCI paper
# To Do: ANOVA for subtype ages 

load("../OvcSubtypes_Files/esets.not.rescaled.classified.RData") 
load("../OvcSubtypes_Files/classification.vals.genepairs.default.cutoff.rf.RData")

esets <- esets.not.rescaled.classified
for (i in 1:length(esets))
  pData(esets[[i]]) <- pData(esets[[i]])[, 
                                         c("age_at_initial_pathologic_diagnosis","days_to_death", "vital_status",
                                           "Konecny.subtypes", "Helland.subtypes", "Verhaak.subtypes"
                                           )
                                         ]


### below code prints all ages  


for (i in 1:length(esets)){
  eset <- esets[[i]]
  accessionID <- sub("_eset", "", names(esets)[i])
  eset.name <- names(esets)[i]
  ##pdata.nonblank will contain pdata columns with any non-missing values:
  pdata.nonblank <- pData(eset)
  pdata.nonblank <- pdata.nonblank[,apply(pdata.nonblank,2,function(x) sum(!is.na(x)) > 0)]
  cat(paste("GEO Accession ID:", accessionID))
  cat("\n")
  ages <- eset$age_at_initial_pathologic_diagnosis
  print("Average age (years):")
  cat("\n")
  print(mean(ages,na.rm = TRUE))
  cat("\n")
  Consensus.subtypes <- classification.vals.genepairs.default.cutoff.rf[[names(esets[i])]]
  Konecny.mean.ages <- sapply(split(ages, eset$Konecny.subtypes), 
                              function (subgroup_ages) mean(subgroup_ages,na.rm=TRUE)) 
  Helland.mean.ages <-  sapply(split(ages, factor(eset$Helland.subtypes,levels = c("C2","C4","C5","C1"))), 
                               function (subgroup_ages) mean(subgroup_ages,na.rm=TRUE))  
  Verhaak.mean.ages <- sapply(split(ages, eset$Verhaak.subtypes), 
                              function (subgroup_ages) mean(subgroup_ages,na.rm=TRUE))
  Consensus.rf.mean.ages <- sapply(split(ages, Consensus.subtypes), 
                              function (subgroup_ages) mean(subgroup_ages,na.rm=TRUE))
  
  print("Konecny subtype ages \n")
  print(Konecny.mean.ages)
  print("Helland subtype ages \n")
  print(Helland.mean.ages)
  print("Verhaak subtype ages \n")
  print(Verhaak.mean.ages)
  print("Consensus subtype ages \n")
  print(Consensus.rf.mean.ages)

}

subtype_ages <- lapply(1:length(esets), 
                       function (i) 
                       {
                         eset <- esets[[i]]
                        Consensus.subtypes <- classification.vals.genepairs.default.cutoff.rf[[names(esets[i])]]
                        ages <- eset$age_at_initial_pathologic_diagnosis
                         Konecny <- split(ages, eset$Konecny.subtypes)  
                         Helland <- split(ages, factor(eset$Helland.subtypes,levels = c("C2","C4","C5","C1")))
                         Verhaak <- split(ages, eset$Verhaak.subtypes)
                         Consensus <- split(ages, Consensus.subtypes)
                        list(Konecny,Helland,Verhaak,Consensus)
                       }
                        )


Konecny_subtypes <- lapply(1:4,function (j) lapply(1:15, function (i) subtype_ages[[i]][[1]][[j]]) %>% unlist )
Helland_subtypes <- lapply(1:4,function (j) lapply(1:15, function (i) subtype_ages[[i]][[1]][[j]]) %>% unlist )
Verhaak_subtypes <- lapply(1:4,function (j) lapply(1:15, function (i) subtype_ages[[i]][[1]][[j]]) %>% unlist )
Consensus_subtypes <- lapply(1:4,function (j) lapply(1:15, function (i) subtype_ages[[i]][[1]][[j]]) %>% unlist )
