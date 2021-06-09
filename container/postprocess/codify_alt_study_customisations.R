source("../container/downstream/DIA_downstream_datacarpentry.R")

# must start with ours_...
#FDR 1%, 'all' inference
# PXD004691
ours_all_PXD004691 <- function(groups, norm="median"){
  annot <- read.delim("../inputs/annotations/PXD004691_annotation_corrected.txt")
  rda <- "../inputs/rdas/fdr1_all_inference/PXD004691_corrected.rda"
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
}

# PXD000672
ours_all_PXD000672 <- function(groups, norm="median"){
  annot <- read.delim("../inputs/annotations/PXD000672_annotation_corrected.txt")
  rda <- "../inputs/rdas/fdr1_all_inference/PXD000672_corrected.rda"
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
}

# PXD014943
ours_all_PXD014943 <- function(groups, norm="median"){
  annot <- read.delim("../inputs/annotations/PXD014943_annotation_corrected_norecalc.txt")
  rda <- "../inputs/rdas/fdr1_all_inference/PXD014943.rda"
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
}
