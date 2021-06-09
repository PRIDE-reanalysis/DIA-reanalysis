source("../container/downstream/DIA_downstream_datacarpentry.R")
source("../container/postprocess/DIA_postprocess_differential.R")
source("../container/postprocess/codify_study_customisations.R")

pxds <- list("PXD014943" = list(c("eDLBCL", "PCNSL") , c("eDLBCL", "IVL")),
             "PXD004691" = list(c("F-N", "F-T"), c("P-N","P-T")),
             "PXD000672" = list(c("ccRCC", "normal_ccRCC"), c("ccRCC", "pRCC")) 
)

for (name in names(pxds)) {
  assign(paste("Contrast",name,sep="_"), calc_contrasts(name, pxds[[name]]) )
}

# fix headers
Contrast_PXD004691$`F-N-F-T` <- Contrast_PXD004691$`F-N-F-T` + labs(title="PXD004691 normal(ff)-PrC(ff)")
Contrast_PXD004691$`P-N-P-T` <- Contrast_PXD004691$`P-N-P-T` + labs(title="PXD004691 normal(pe)-PrC(pe)")
Contrast_PXD000672$`ccRCC-normal_ccRCC` <- Contrast_PXD000672$`ccRCC-normal_ccRCC` + labs(title="PXD000672 ccRCC-normal")

# make figure panel
panel_de <- (Contrast_PXD014943$`eDLBCL-PCNSL` + ylim(0,30) | Contrast_PXD014943$`eDLBCL-IVL`+ ylim(0,30) ) / 
  (Contrast_PXD004691$`F-N-F-T` + ylim(0,25) | Contrast_PXD004691$`P-N-P-T` + ylim(0,25) ) / 
  (Contrast_PXD000672$`ccRCC-normal_ccRCC` + ylim(0,15) | Contrast_PXD000672$`ccRCC-pRCC`+ ylim(0,15) ) / 
  plot_annotation(tag_levels = 'A', title = 'Volcano plots for differential expression')
# panel_de
# point size should be 2

ggsave(
  "dia_paper_volcano_panel_top3update.pdf",
  plot = panel_de,
  device = "pdf",
  scale = 1,
  width = 38,
  height = 38,
  units = "cm" ,
  dpi = 300,
  limitsize = TRUE
)


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


# recalc differential expressions with top3 imputation
source("FC_variance_results_datacarpentry.R")

PXD000672_1_1_top3 <- calc_contrasts("1_1_top3_PXD000672", list(c("ccRCC", "normal_ccRCC"), c("ccRCC", "pRCC")) )
PXD004691_1_1_top3 <- calc_contrasts("1_1_top3_PXD004691", list(c("F-N", "F-T"), c("P-N","P-T")) )
PXD014943_1_1_top3 <- calc_contrasts("1_1_top3_PXD014943", list(c("eDLBCL", "PCNSL") , c("eDLBCL", "IVL")) )

# fix headers
PXD004691_1_1_top3$`F-N-F-T` <- PXD004691_1_1_top3$`F-N-F-T` + labs(title="PXD004691 normal(ff)-PrC(ff)")
PXD004691_1_1_top3$`P-N-P-T` <- PXD004691_1_1_top3$`P-N-P-T` + labs(title="PXD004691 normal(pe)-PrC(pe)")
PXD000672_1_1_top3$`ccRCC-normal_ccRCC` <- PXD000672_1_1_top3$`ccRCC-normal_ccRCC` + labs(title="PXD000672 ccRCC-normal")
PXD000672_1_1_top3$`ccRCC-pRCC` <- PXD000672_1_1_top3$`ccRCC-pRCC` + labs(title="PXD000672 ccRCC-pRCC")
PXD014943_1_1_top3$`eDLBCL-PCNSL` <- PXD014943_1_1_top3$`eDLBCL-PCNSL` + labs(title="PXD014943 eDLBCL-PCNSL")
PXD014943_1_1_top3$`eDLBCL-IVL` <- PXD014943_1_1_top3$`eDLBCL-IVL` + labs(title="PXD014943 eDLBCL-IVL")

# make figure panel
panel_de <- (PXD014943_1_1_top3$`eDLBCL-PCNSL` + ylim(0,30) | PXD014943_1_1_top3$`eDLBCL-IVL`+ ylim(0,30) ) / 
  (PXD004691_1_1_top3$`F-N-F-T` + ylim(0,25) | PXD004691_1_1_top3$`P-N-P-T` + ylim(0,25) ) / 
  (PXD000672_1_1_top3$`ccRCC-normal_ccRCC` + ylim(0,15) | PXD000672_1_1_top3$`ccRCC-pRCC`+ ylim(0,15) ) / 
  plot_annotation(tag_levels = 'A', title = 'Volcano plots for differential expression')

