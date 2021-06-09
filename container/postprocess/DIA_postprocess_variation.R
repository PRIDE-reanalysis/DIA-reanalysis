library(tidyverse)
library(broom)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

calc_cv_by_tr <- function(df){
  cvdf <- df %>%
    dplyr::select(GROUP_ORIGINAL,SUBJECT_ORIGINAL,Protein,NormLogIntensities,technical_replicate) %>%
    dplyr::group_by(GROUP_ORIGINAL,SUBJECT_ORIGINAL,Protein) %>%
    dplyr::mutate(Intensities = `^`(2,NormLogIntensities)) %>%
    dplyr::summarise(tr_sd = sd(Intensities), tr_mean = mean(Intensities), tr_count = n(), .groups = "keep") %>%
    dplyr::filter(tr_count > 1,) %>%
    dplyr::filter(!is.na(tr_sd)) %>%
    dplyr::mutate(cv = tr_sd/tr_mean*100)
}

analyse_cv <- function(df) {
  cv_tab <- calc_cv_by_tr(df)
  cv_all <- cv_tab %>% dplyr::mutate(GROUP_ORIGINAL = "Total")
  
  # Groups
  groups <- cv_tab %>% dplyr::group_by(GROUP_ORIGINAL) %>%
    dplyr::summarise(cv_gr = median(cv), .groups = "keep")
  
  # Subjects
  subjects <- cv_tab %>% dplyr::group_by(SUBJECT_ORIGINAL) %>%
    dplyr::summarise(cv_gr = median(cv), .groups = "keep")
  
  # Overall 
  oamedian <- median(cv_tab$cv)
  
  # Violin plots
  plot <- ggplot(rbind(cv_tab,cv_all) %>% dplyr::group_by(GROUP_ORIGINAL, .groups = "keep" ),
                 aes(x = GROUP_ORIGINAL, y = cv, fill=GROUP_ORIGINAL) ) + 
    geom_violin(trim=TRUE) +
    geom_boxplot(width=0.1) + 
    labs(title = "Coefficient of variation", x = "Types", y = "CV [%]") + 
    ylim(0,100) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return_list <- list("t_groups_cv" = groups, "t_subjects_cv" = subjects, "median_cv" = oamedian, "p_cv" = plot)
  return(return_list) 
}
