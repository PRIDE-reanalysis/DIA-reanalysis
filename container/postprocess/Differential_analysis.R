source("../container/downstream/DIA_downstream_datacarpentry.R")
source("../container/postprocess/DIA_postprocess_differential.R")
source("../container/postprocess/codify_study_customisations.R")

pxds <- list("PXD014943" =  list(c("eDLBCL", "PCNSL") , c("IVL", "eDLBCL")),
             "PXD004691" = list(c("F-N", "F-T"), c("P-N","P-T")),
             "PXD000672" = list(c("normal_ccRCC", "ccRCC"), c("ccRCC", "pRCC")) 
)

for (name in names(pxds)) {
  assign(paste("Contrast",name,sep="_"), calc_contrasts(name, pxds[[name]]) )
}

# fix headers
Contrast_PXD004691$`F-N-F-T` <- Contrast_PXD004691$`F-N-F-T` + labs(title="PXD004691 normal(ff)-PrC(ff)")
Contrast_PXD004691$`P-N-P-T` <- Contrast_PXD004691$`P-N-P-T` + labs(title="PXD004691 normal(pe)-PrC(pe)")
Contrast_PXD000672$`normal_ccRCC-ccRCC` <- Contrast_PXD000672$`normal_ccRCC-ccRCC` + labs(title="PXD000672 normal-ccRCC")

# make figure panel
panel_de <- (Contrast_PXD014943$`eDLBCL-PCNSL`| Contrast_PXD014943$`IVL-eDLBCL`) / 
  (Contrast_PXD004691$`F-N-F-T` | Contrast_PXD004691$`P-N-P-T` ) / 
  (Contrast_PXD000672$`normal_ccRCC-ccRCC` | Contrast_PXD000672$`ccRCC-pRCC`) / 
  plot_annotation(tag_levels = 'A')
# panel_de

# write protein tables
# not working :/
# for (dataset in c(Contrast_PXD014943,Contrast_PXD004691,Contrast_PXD000672)) {
#   for (name in names(dataset)) {
#     print(name)    
#     if (startsWith(name, "contrast")) {
#       data <- dataset[[as.name(name)]]
#       print(paste(paste(deparse(substitute(dataset)), name, sep=" "), ".tsv", sep=""))
#       # write.table(data %>% dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
#       #                                                        if_else(log.fc <= (-1*fc_threshold), "right", 
#       #                                                                if_else(log.fc >= fc_threshold, "left" , "center")
#       #                                                        ))) %>% filter(!threshold == "center") ,
#       #           paste(paste(deparse(substitute(dataset)), name, sep=" "), ".tsv", sep=""),
#       #           quote = FALSE, row.names = FALSE,
#       #           sep = "\t")
#     }
#   }
# }
# manually then
fc_threshold <- 1
write.table(Contrast_PXD014943$`contrast-eDLBCL-PCNSL` %>% 
              dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
                  if_else(log.fc <= (-1*fc_threshold), "right", 
                          if_else(log.fc >= fc_threshold, "left" , "center")
                  ))) %>% filter(!threshold == "center") ,
            "PXD014943_DLBCL-PCNSL_log2fc>1.tsv",
            quote = FALSE, row.names = FALSE,
            sep = "\t")
	    
write.table(Contrast_PXD014943$`contrast-IVL-eDLBCL` %>% 
              dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
                  if_else(log.fc <= (-1*fc_threshold), "right", 
                          if_else(log.fc >= fc_threshold, "left" , "center")
                  ))) %>% filter(!threshold == "center") ,
            "PXD014943_IVL-eDLBCL_log2fc>1.tsv",
            quote = FALSE, row.names = FALSE,
            sep = "\t")


write.table(Contrast_PXD004691$`contrast-F-N-F-T` %>% 
              dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
                  if_else(log.fc <= (-1*fc_threshold), "right", 
                          if_else(log.fc >= fc_threshold, "left" , "center")
                  ))) %>% filter(!threshold == "center") ,
            "PXD004691_F-N-F-T_log2fc>1.tsv",
            quote = FALSE, row.names = FALSE,
            sep = "\t")
	    
write.table(Contrast_PXD004691$`contrast-P-N-P-T` %>% 
              dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
                  if_else(log.fc <= (-1*fc_threshold), "right", 
                          if_else(log.fc >= fc_threshold, "left" , "center")
                  ))) %>% filter(!threshold == "center") ,
            "PXD004691_P-N-P-T_log2fc>1.tsv",
            quote = FALSE, row.names = FALSE,
            sep = "\t")


write.table(Contrast_PXD000672$`contrast-normal_ccRCC-ccRCC` %>% 
              dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
                  if_else(log.fc <= (-1*fc_threshold), "right", 
                          if_else(log.fc >= fc_threshold, "left" , "center")
                  ))) %>% filter(!threshold == "center") ,
            "PXD000672_normal_ccRCC-ccRCC_log2fc>1.tsv",
            quote = FALSE, row.names = FALSE,
            sep = "\t")
	    
write.table(Contrast_PXD000672$`contrast-ccRCC-pRCC` %>% 
              dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
                  if_else(log.fc <= (-1*fc_threshold), "right", 
                          if_else(log.fc >= fc_threshold, "left" , "center")
                  ))) %>% filter(!threshold == "center") ,
            "PXD000672ccRCC-pRCC_log2fc>1.tsv",
            quote = FALSE, row.names = FALSE,
            sep = "\t")

