### Fig. 4

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
   plot_annotation(tag_levels = 'a')

ggsave(
  "dia_paper_correlation_panel_top3update.pdf",
  plot = panel_cr,
  device = "pdf",
  scale = 1,
  width = 38,
  height = 38,
  units = "cm" ,
  dpi = 300,
  limitsize = TRUE
)

### Fig. 4 with fixed axis (all plots in panel same reference axis), outlier (< log2(0.001)) treated, plot point size adjusted 

pxds <- c("PXD003497","PXD014194","PXD004873")
for (p in pxds ){
  print(p)
  df1 <- get(paste("ours_",p,sep=""))() %>% filter(feature_missingrate_per_group < 50) %>% filter(LogIntensities > log2(0.001))
  assign(paste("Correlation_reanalysis",p,sep="_"), correlate_techreps(df1, pntsz=0.75))
  df2 <- get(paste("theirs_",p,sep=""))() %>% mutate(NormLogIntensities = LogIntensities) # we have to assume their intensities are normalised and filtered
  assign(paste("Correlation_original",p,sep="_"), correlate_techreps(df2, pntsz=0.75))
}

panel_cr <- (Correlation_reanalysis_PXD003497 + xlim(0,25) + ylim(0,25) + labs(title = "PXD003497", subtitle = "(Reanalysis)") | 
               Correlation_original_PXD003497 + xlim(0,25) + ylim(0,25) + labs(title = "PXD003497", subtitle = "(Original)") ) / 
  (Correlation_reanalysis_PXD004873 + xlim(0,25) + ylim(0,25) + labs(title = "PXD004873", subtitle = "(Reanalysis)") | 
     Correlation_original_PXD004873 + xlim(0,25) + ylim(0,25) + labs(title = "PXD004873", subtitle = "(Original)") ) / 
  (Correlation_reanalysis_PXD014194 + xlim(0,25) + ylim(0,25) + labs(title = "PXD014194", subtitle = "(Reanalysis)") | 
     Correlation_original_PXD014194 + xlim(0,25) + ylim(0,25) + labs(title = "PXD014194", subtitle = "(Original)") )/ 
  plot_annotation(tag_levels = 'a')
  #plot_annotation(tag_levels = 'a', title = 'Protein correlation in technical replicates')

ggsave(
  "dia_paper_correlation_panel_top3update_axisfix_outliertreated.pdf",
  plot = panel_cr,
  device = "pdf",
  scale = 1,
  width = 38,
  height = 38,
  units = "cm" ,
  dpi = 300,
  limitsize = TRUE
)

### Suppl. Fig. S3

Correlation_runs_PXD004873 <- correlate_study(ours_PXD004873(),theirs_PXD004873())
panel_ic <- Correlation_runs_PXD004873 + xlab("Replicate 1 (log2 intensity)") + ylab("Replicate 2 (log2 intensity)")

ggsave(
  "dia_paper_PXD004873_ProtIntCorrel_top3update_outliertreated.pdf",
  plot = panel_ic,
  device = "pdf",
  scale = 1,
  width = 38,
  height = 76,
  units = "cm" ,
  dpi = 300,
  limitsize = TRUE
)

### Suppl. Fig. S4

source("FC_variance_results_datacarpentry.R")
source("../container/postprocess/DIA_postprocess_correlation.R")
source("../container/postprocess/codify_study_customisations.R")
source("../container/postprocess/codify_originalresult_integration.R")

pxds <- c("PXD003497","PXD014194","PXD004873")
for (p in pxds ){
  print(p)
  df1 <- get(paste("ours_1_1_top3_",p,sep=""))() %>% filter(feature_missingrate_per_group < 50) %>% filter(LogIntensities > log2(0.001))
  assign(paste("Correlation_reanalysis_top3",p,sep="_"), correlate_techreps(df1))
  df2 <- get(paste("theirs_",p,sep=""))() %>% mutate(NormLogIntensities = LogIntensities) # we have to assume their intensities are normalised and filtered
  assign(paste("Correlation_original",p,sep="_"), correlate_techreps(df2))
}

(Correlation_reanalysis_top3_PXD003497 + xlim(0,25) + ylim(0,25) + labs(title = "PXD003497", subtitle = "(Reanalysis)") | 
    Correlation_original_PXD003497 + xlim(0,25) + ylim(0,25) + labs(title = "PXD003497", subtitle = "(Original)") ) / 
  (Correlation_reanalysis_top3_PXD004873 + xlim(0,25) + ylim(0,25) + labs(title = "PXD004873", subtitle = "(Reanalysis)") | 
     Correlation_original_PXD004873  + xlim(0,25) + ylim(0,25) + labs(title = "PXD004873", subtitle = "(Original)") ) / 
  (Correlation_reanalysis_top3_PXD014194 + xlim(0,25) + ylim(0,25) + labs(title = "PXD014194", subtitle = "(Reanalysis)") | 
     Correlation_original_PXD014194 + xlim(0,25) + ylim(0,25) + labs(title = "PXD014194", subtitle = "(Original)") )/ 
  plot_annotation(tag_levels = 'a', title = 'Protein correlation in technical replicates')


