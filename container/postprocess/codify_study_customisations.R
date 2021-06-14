library(tidyverse)
library(stringr)
library(data.table)
source("../container/downstream/DIA_downstream_datacarpentry.R")

# The functions tie together file load of annotation and rda, 
# and preprocess the dataframe such that it can be used in the 
# anlaysis R code of `postprocess` 
# ! They are not generalised since the annotation is not uniform
# ! (or was not during MSstats), or MSstats object name is not  
# ! uniform and some extra massaging is needed to get the
# ! dataframe in an (postprocessing) acceptable shape.

# It basically follows this scheme in sequence:
#   1. group_consistency
#   2. clean_protein_names
#   3. normalise 
#   4. merge with annotations 
#   5. kick out superfluous cols
#   6. dont return a grouped dataframe
#   extras: you can build in a group filter argument or 
#           a normalisation method selector argument
  

tally_protpepirt <- function(rda, annot) {
  irt_peptides <- nrow(rda$ProcessedData %>% 
                         dplyr::filter(PROTEIN=="1/iRT_protein") %>% 
                         dplyr::mutate(irtp =  gsub("_.*","",PEPTIDE) ) %>% dplyr::distinct(irtp)
                       )
  prot_all <- length(levels(rda$RunlevelData$Protein))
  
  prot_50 <- nrow(rda$RunlevelData %>% 
                    dplyr::filter(!more50missing == TRUE) %>% 
                    dplyr::distinct(Protein)
                  )

  table_procdat <- data.frame(
    Analysis = c("Total num of proteins", 
                 "Total num of proteins with <50% missing values",
                 "Total num of iRT peptides",
                 "Total num of peptides"),
    measure = c(prot_all,
                prot_50,
                irt_peptides,
                length(levels(rda$ProcessedData$PEPTIDE)))
  )
  return(table_procdat)
  
}

# PXD014194
ours_PXD014194 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD014194_annotation_corrected.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD014194.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
           group_consistency(annot) %>%
           clean_protein_names() %>% 
           dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
           {if(norm == "median") run_median_norm(.) else .} %>%
           {if(norm == "quantile") run_quantile_norm(.) else .} %>%
           {if(norm == "zscore") run_zscore_norm(.) else .} %>%
           {if(norm == "minmax") run_minmax_norm(.) else .} %>%
           dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
           dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                         feature_missingrate_per_group,feature_missingrate_per_run,technical_replicate) %>%
           dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                              dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                              dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD010912
ours_PXD010912 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD010912_annotation.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD010912.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run) %>%
    dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD004873
ours_PXD004873 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD004873_annotation_corrected.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD004873.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run,technical_replicate) %>%
    dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD004691
ours_PXD004691 <- function(groups, norm="median", tally=FALSE){
  # annot <- read.delim("../inputs/annotations/PXD004691_annotation.txt") %>% dplyr::select(Run,Condition,BioReplicate,technical_replicate,sample_name)
  annot <- read.delim("../inputs/annotations/PXD004691_annotation_corrected.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD004691.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run,technical_replicate) %>%
    dplyr::ungroup()

  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD004589
ours_PXD004589 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD004589_annotation.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD004589.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run) %>%
    dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD003497
ours_PXD003497 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD003497_annotation_corrected.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD003497.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::rename(sample_id = sample_name) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run,technical_replicate,sample_id) %>%
    dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD000672
ours_PXD000672 <- function(groups, norm="median", tally=FALSE){
  # annot <- read.delim("../inputs/annotations/PXD000672_annotation.txt")
  # rda <- "../inputs/rdas/fdr1/PXD000672.rda"
  annot <- read.delim("../inputs/annotations/PXD000672_annotation_corrected.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD000672.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    dplyr::mutate(GROUP_ORIGINAL = str_trim(GROUP_ORIGINAL)) %>%  # needs to be after groups filter BUT BEWARE! groups argument must be original names (incl. trailing whitespaces)
    dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run) %>%
    dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD001064
ours_PXD001064 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD001064_serum_annotation.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD001064_serum.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run) %>%
    dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD003539
ours_PXD003539 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD003539_annotation.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD003539.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run) %>%
    dplyr::ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD014943
ours_PXD014943 <- function(groups, norm="median", tally=FALSE){
  annot <- read.delim("../inputs/annotations/PXD014943_annotation_corrected_norecalc.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD014943.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  our<-goldstandard.proposed$RunlevelData %>% 
    clean_protein_names() %>% 
    group_consistency(annot) %>%
    filter(GROUP_ORIGINAL %in% groups) %>%
    {if(norm == "median") run_median_norm(.) else .} %>%
    {if(norm == "quantile") run_quantile_norm(.) else .} %>%
    {if(norm == "zscore") run_zscore_norm(.) else .} %>%
    {if(norm == "minmax") run_minmax_norm(.) else .} %>%
    full_join(annot, by = c("originalRUN" = "Run")) %>%
    dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                  feature_missingrate_per_group,feature_missingrate_per_run) %>%
    ungroup()
  
  if(tally) {
    procdat <- tally_protpepirt(goldstandard.proposed, annot)   
    rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
          c("Total num of Proteins <50% features missing in group", nrow(our %>% 
                                                                           dplyr::filter(feature_missingrate_per_group < 50) %>% 
                                                                           dplyr::distinct(Protein))))
  }
  else {
    return(our)
  }
}

# PXD000
ours_test <- function(diff, groups, norm="median") {
  annot <- read.delim("../inputs/annotations/PXD014943_annotation_corrected_norecalc.txt")
  rda <- "../inputs/intermediate-results/msstats_rdas/fdr1_top3_inference/PXD014943.rda"
  load(rda)
  
  if(missing(groups)) {
    groups <- as.character((goldstandard.proposed$RunlevelData %>% distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  
  if(missing(diff)) {
    return(our<-goldstandard.proposed$RunlevelData %>% 
             clean_protein_names() %>% 
             group_consistency() %>%
             dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
             {if(norm == "median") run_median_norm(.) else .} %>%
             {if(norm == "quantile") run_quantile_norm(.) else .} %>%
             {if(norm == "zscore") run_zscore_norm(.) else .} %>%
             {if(norm == "minmax") run_minmax_norm(.) else .} %>%
             dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
             dplyr::select(Protein,LogIntensities,NormLogIntensities,originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                           feature_missingrate_per_group,feature_missingrate_per_run,pct_batch,tissue_region) %>%
             dplyr::ungroup())
  } else {
    if (diff == "msstatsquant") {
      return(our<-quantification(goldstandard.proposed$RunlevelData,  type="Sample", format="long") %>% 
               clean_protein_names() %>% 
               group_consistency(annot) %>%
               run_median_norm() %>%
               full_join(annot, by = c("originalRUN" = "Run")) %>%
               dplyr::select(Protein,LogIntensities,matches("NormLogIntensities"),originalRUN,GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                             feature_missingrate_per_group,feature_missingrate_per_run,pct_batch,tissue_region)
      )
    } else {
      return(our<-goldstandard.proposed$RunlevelData %>% 
               clean_protein_names() %>% 
               run_median_norm() %>%
               rename(OriLogIntensities = LogIntensities, LogIntensities = NormLogIntensities)
      )
    }
  }
}
