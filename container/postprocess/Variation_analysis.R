### Fig. 3

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



panel_cv <- (Variances_reanalysis_PXD003497$p_cv + 
    labs(title = "PXD003497", subtitle = "(Reanalysis)") + 
    theme(legend.position = "none") + 
    scale_x_discrete(labels=c("Total","AAC", "BPH", "DAC" )) | 
    Variances_original_PXD003497$p_cv + 
    labs(title = "PXD003497", subtitle = "(Original)") + 
    theme(legend.position = "none", legend.title=element_blank()) + 
    scale_x_discrete(labels=c("Total","AAC", "BPH", "DAC" )) ) / 
  (Variances_reanalysis_PXD004873$p_cv + 
     labs(title = "PXD004873", subtitle = "(Reanalysis)") + 
     theme(legend.position = "none") + 
     scale_x_discrete(labels=c("Normal", "Tumour", "Total")) | 
     Variances_original_PXD004873$p_cv + 
     labs(title = "PXD004873", subtitle = "(Original)") + 
     theme(legend.position = "none", legend.title=element_blank())  + 
     scale_x_discrete(labels=c("Normal", "Tumour", "Total" )) )/ 
  (Variances_reanalysis_PXD014194$p_cv + labs(title = "PXD014194", subtitle = "(Reanalysis)") + theme(legend.position = "none") | 
     Variances_original_PXD014194$p_cv + labs(title = "PXD014194", subtitle = "(Original)") + theme(legend.position = "none", legend.title=element_blank()) ) / 
  plot_annotation(tag_levels = 'a')
  #plot_annotation(tag_levels = 'a', title = 'Coefficient of variation in technical replicates')

ggsave(
  "dia_paper_cv_panel_top3update.pdf",
  plot = panel_cv,
  device = "pdf",
  scale = 1,
  width = 38,
  height = 38,
  units = "cm" ,
  dpi = 300,
  limitsize = TRUE
)


Variances_reanalysis_PXD004691 <- analyse_cv(ours_PXD004691() %>% filter(feature_missingrate_per_group < 50) )
# No reproduction of variance from original for PXD004691 possible, data starts at sample level (tech_rep already merged)
