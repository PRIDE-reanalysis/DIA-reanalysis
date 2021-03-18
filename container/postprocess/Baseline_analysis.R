library(tidyverse)

source("../container/downstream/DIA_downstream_datacarpentry.R")
source("../container/postprocess/codify_study_customisations.R")

pxds <- c("PXD004873", "PXD000672", "PXD004691", "PXD014943", "PXD003497", "PXD004589", "PXD014194", "PXD003539", "PXD001064", "PXD010912")

for ( i in pxds ) {
  data <- get(paste("ours_",i,sep=""))(norm="median") %>% 
  filter(feature_missingrate_per_group < 50) %>%
  mutate(NormIntensities = `^`(2,NormLogIntensities)) %>%
  select(-LogIntensities,-NormLogIntensities,-GROUP_ORIGINAL, -SUBJECT_ORIGINAL, -feature_missingrate_per_run, -feature_missingrate_per_group) %>%
  select(-matches("technical_replicate"), -matches("sample_id"), -matches("sample_name")) %>%
  spread(originalRUN, NormIntensities) 

  write.table(data,
              paste(i,"_baseline_mediannorm_untransformed.tsv",sep=""),
              quote = FALSE,
              sep = "\t",
              row.names = F)
}
