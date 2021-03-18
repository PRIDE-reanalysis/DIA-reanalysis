source("../container/postprocess/DIA_postprocess_correlation.R")
source("../container/postprocess/codify_study_customisations.R")
source("../container/postprocess/codify_originalresult_integration.R")

pxds <- c("PXD003497","PXD014194","PXD004873")
for (p in pxds ){
  print(p)
  df1 <- get(paste("ours_",p,sep=""))() %>% filter(feature_missingrate_per_group < 50) 
  assign(paste("Correlation_reanalysis",p,sep="_"), correlate_techreps(df1))
  df2 <- get(paste("theirs_",p,sep=""))() %>% mutate(NormLogIntensities = LogIntensities) # we have to assume their intensities are normalised and filtered
  assign(paste("Correlation_original",p,sep="_"), correlate_techreps(df2))
}


(Correlation_reanalysis_PXD003497 | Correlation_original_PXD003497 ) / 
  (Correlation_reanalysis_PXD004873 | Correlation_original_PXD004873 ) / 
  (Correlation_reanalysis_PXD014194 | Correlation_original_PXD014194 ) / 
  plot_annotation(tag_levels = 'A')
