source("../container/postprocess/DIA_postprocess_variation.R")
source("../container/postprocess/codify_study_customisations.R")
source("../container/postprocess/codify_originalresult_integration.R")

pxds <- c("PXD003497","PXD014194","PXD004873")
for (p in pxds ){
  print(p)
  df1 <- get(paste("ours_",p,sep=""))() %>% filter(feature_missingrate_per_group < 50) 
  assign(paste("Variances_reanalysis",p,sep="_"), analyse_cv(df1))
  df2 <- get(paste("theirs_",p,sep=""))() %>% rename(NormLogIntensities = LogIntensities) # we have to assume their intensities are normalised and filtered
  assign(paste("Variances_original",p,sep="_"), analyse_cv(df2))
}

(Variances_reanalysis_PXD003497$p_cv | Variances_original_PXD003497$p_cv ) / 
  (Variances_reanalysis_PXD004873$p_cv | Variances_original_PXD004873$p_cv ) / 
  (Variances_reanalysis_PXD014194$p_cv | Variances_original_PXD014194$p_cv ) / 
  plot_annotation(tag_levels = 'A')

Variances_reanalysis_PXD004691 <- analyse_cv(ours_PXD004691() %>% filter(feature_missingrate_per_group < 50) )
# No reproduction of variance from original for PXD004691 possible, data starts at sample level (tech_rep already merged)
