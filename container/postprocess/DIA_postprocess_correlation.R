library(tidyverse)
library(stringr)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(patchwork)

library(MASS, mask.ok = FALSE) 

# Calculates the density for a given number of 2D points
# x: vector of 1st dimension coordinate items 
# y: vector of 2nd dimension coordinate items
get_density <- function(x, y, ...) {
  # creds to slowkow(dot)com/notes/ggplot2-color-by-density/
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Creates correlation plot for all technical replicates pairs in a study
# df: df long form, needs cols `c(LogIntensities,Protein,originalRUN,technical_replicate or sample_id)`
# ! technical_replicate tuples must be named a/b/...
# ! correlates only a with b (c.. is ignored)
correlate_techreps <- function(df) {
  # in case different technical replicates are made, an optional sample_id column can be used to (uniquely) identifiy TR pairs
  if("sample_id" %in% names(df)) {
    corr <- df %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(Protein,sample_id,sep="_")) %>%
      dplyr::select(technical_replicate,LogIntensities,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>% 
      tidyr::drop_na(c(a,b)) 
  } else {
    corr <- df %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(GROUP_ORIGINAL,SUBJECT_ORIGINAL,sep="")) %>%
      dplyr::select(technical_replicate,LogIntensities,Protein,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>% 
      tidyr::drop_na(c(a,b)) 
  }
  corr$density <- get_density(corr$a, corr$b, n = 100)
  
  plot <- ggplot(corr, aes(x = a, y = b, color = density)) + 
    geom_point() + 
    #scale_fill_gradient(low="blue", high="red") +
    scale_color_viridis_c(option = "magma", guide=FALSE) + 
    geom_smooth(method=lm) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + 
    labs(title = "Protein correlation in technical replicates", x = "Replicate 1 (log2 intensity)", y = "Replicate 2 (log2 intensity)") 
  
  return(plot)
}

# Creates correlation plots for all runs in a study analysis method x vs analysis method y
# ours: df argument, long form, needs cols `c(LogIntensities,Protein,originalRUN)`
# theirs: df argument, long form, needs cols `c(LogIntensities,Protein,originalRUN)`
correlate_study <- function(ours,theirs) {
  # merge ours with theirs
  corr <- ours %>% 
    dplyr::mutate(merge_id = paste(Protein,originalRUN,sep="")) %>%
    dplyr::select(LogIntensities,Protein,merge_id) %>% 
    dplyr::inner_join(
      theirs %>% 
        dplyr::mutate(merge_id = paste(Protein,originalRUN,sep="")) %>%
        dplyr::select(LogIntensities,Protein,originalRUN,merge_id), 
      by="merge_id"
    ) %>% tidyr::drop_na() %>%
    dplyr::select(-merge_id) %>%
    dplyr::rename(run = originalRUN) %>%
    dplyr::filter(is.finite(LogIntensities.x)) %>%
    dplyr::filter(is.finite(LogIntensities.y))
  
  outlier_threshold <- log2(0.05)
  corr_w_out <- corr %>% mutate(outliers =  LogIntensities.x < outlier_threshold | LogIntensities.y < outlier_threshold) %>% 
                        group_by(run) %>% summarise(outliers = sum(outliers)) %>%
                        filter(outliers > 0)
  print("Removed this many 'close to zero' intensity pairs")
  print(corr_w_out)
    
  # crude outlier removal 
  corr <- corr %>%
    dplyr::filter(LogIntensities.x >= outlier_threshold) %>%
    dplyr::filter(LogIntensities.y >= outlier_threshold)
  
  corr$density <- get_density(corr$LogIntensities.x, corr$LogIntensities.y, n = 100)
  
  plot <- ggplot(corr, aes(x = LogIntensities.x, y = LogIntensities.y, color = density)) + 
    geom_point() + 
    #scale_fill_gradient(low="blue", high="red") +
    scale_color_viridis_c(option = "magma", guide=FALSE) + 
    geom_smooth(method=lm) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + 
    labs(title = "Protein correlation between study analyses", x = "Reanalysis", y = "Original") + 
    facet_wrap(vars(run), ncol=4)
  return(plot)
} 


correlate_study_trs <- function(ours,theirs) {
  
  if("sample_id" %in% names(ours)) {
    corr_o <- df %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(Protein,sample_id,sep="_")) %>%
      dplyr::select(technical_replicate,LogIntensities,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>%
      dplyr::rename(a_ours = a, b_ours = b)
  } else {
    corr_o <- ours %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(GROUP_ORIGINAL,SUBJECT_ORIGINAL,sep="")) %>%
      dplyr::select(technical_replicate,LogIntensities,Protein,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>%
      dplyr::rename(a_ours = a, b_ours = b)
  }
  if("sample_id" %in% names(theirs)) {
    corr_t <- df %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(Protein,sample_id,sep="_")) %>%
      dplyr::select(technical_replicate,LogIntensities,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>%
      dplyr::rename(a_theirs = a, b_theirs = b)
  } else {
    corr_t <- theirs %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(GROUP_ORIGINAL,SUBJECT_ORIGINAL,sep="")) %>%
      dplyr::select(technical_replicate,LogIntensities,Protein,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>%
      dplyr::rename(a_theirs = a, b_theirs = b)
  }
  
  # merge ours with theirs
  corr_ot <- corr_o %>% inner_join(corr_t, by=c("sample_id","Protein")) %>% 
    dplyr::select(a_ours, b_theirs) %>%
    dplyr::rename(a = a_ours, b = b_theirs) %>%
    tidyr::drop_na(c(a,b)) 
  corr_to <- corr_t %>% inner_join(corr_o, by=c("sample_id","Protein")) %>% 
    dplyr::select(a_theirs, b_ours) %>%
    dplyr::rename(a = a_theirs, b = b_ours) %>%
    tidyr::drop_na(c(a,b)) 
  
  # paste ot and to to one and calc density for all three
  corr <- rbind(corr_ot, corr_to %>% dplyr::rename(a = b, b = a))
  corr$density <- get_density(corr$a, corr$b, n = 100)
  corr_ot$density <- get_density(corr_ot$a, corr_ot$b, n = 100)
  corr_to$density <- get_density(corr_to$a, corr_to$b, n = 100)
  
  plot_ot <- ggplot(corr_ot, aes(x = a, y = b, color = density)) + 
    geom_point() + 
    #scale_fill_gradient(low="blue", high="red") +
    scale_color_viridis_c(option = "magma", guide=FALSE) + 
    geom_smooth(method=lm) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + 
    labs(title = "Protein correlation in technical replicates \n reanalysis vs. original", 
         x = "Replicate 1 (reanalysis)", y = "Replicate 2 (original)") 
  
  plot_to <- ggplot(corr_to, aes(x = a, y = b, color = density)) + 
    geom_point() + 
    #scale_fill_gradient(low="blue", high="red") +
    scale_color_viridis_c(option = "magma", guide=FALSE) + 
    geom_smooth(method=lm) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + 
    labs(title = "Protein correlation in technical replicates \n reanalysis vs. original", 
         x = "Replicate 1 (original)", y = "Replicate 2 (reanalysis)") 
  
  plot <- ggplot(corr, aes(x = a, y = b, color = density)) + 
    geom_point() + 
    #scale_fill_gradient(low="blue", high="red") +
    scale_color_viridis_c(option = "magma", guide=FALSE) + 
    geom_smooth(method=lm) + 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + 
    labs(title = "Protein correlation in technical replicates \n reanalysis vs. original", 
         x = "Replicate 1", y = "Replicate 2") 
  
  return_list <- list("corr_R1_ours_R2_theirs" = plot_ot,
                      "corr_R1_theirs_R2_ours" = plot_to,
                      "corr_R1_ours_R2_theirs_and_R1_theirs_R2_ours" = plot)
  return(return_list) 
}

# Example
# ours <- ours_PXD004873()
# theirs <- theirs_PXD004873()
# p_corr_ours <- correlate_techreps(ours)
# p_corr_vs <- correlate_study(ours,theirs)

correlate_techreps_logratios <- function(df) {
  # in case different technical replicates are made, an optional sample_id column can be used to (uniquely) identifiy TR pairs
  if("sample_id" %in% names(df)) {
    corr <- df %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(Protein,sample_id,sep="_")) %>%
      dplyr::select(technical_replicate,LogIntensities,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>% 
      tidyr::drop_na(c(a,b)) 
  } else {
    corr <- df %>% 
      dplyr::filter(technical_replicate %in% c('a','b')) %>% 
      dplyr::mutate(sample_id = paste(GROUP_ORIGINAL,SUBJECT_ORIGINAL,sep="")) %>%
      dplyr::select(technical_replicate,LogIntensities,Protein,sample_id) %>%
      dplyr::filter(is.finite(LogIntensities)) %>%
      tidyr::spread(technical_replicate, LogIntensities) %>% 
      tidyr::drop_na(c(a,b)) 
  }
  
  corr$ab <- corr$a - corr$b
  corr$density <- get_density(corr$ab, corr$b, n = 100)
  
  plot <- ggplot(corr, aes(x = b, y = ab, color = density)) + 
    geom_point() + 
    #scale_fill_gradient(low="blue", high="red") +
    scale_color_viridis_c(option = "magma", guide=FALSE) + 
    geom_smooth(method=lm) + 
    #stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + 
    labs(title = "Protein correlation in technical replicates", y = "Log2(Replicate 1/Replicate 2)", x = "Replicate 2") +
    geom_hline(yintercept=0, linetype="dashed", color = "red")
  box <- ggplot(corr, aes(x=0, y = ab)) + 
    geom_boxplot(alpha = 0.80) +
    scale_x_discrete() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines"),
          plot.background = element_blank()) +
    geom_hline(yintercept=0, linetype="dashed", color = "red")
  
  return((plot | box ) + plot_layout(widths = c(10, 1)))
}
