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
    labs(title = "Protein correlation in technical replicates", x = "Replicate 1", y = "Replicate 2") 
  
  return(plot)
}

# Creates correlation plots for all runs in a study analysis method x vs analysis method y
# ours: df argument, long form, needs cols `c(LogIntensities,Protein,originalRUN)`
# theirs: df argument, long form, needs cols `c(LogIntensities,Protein,originalRUN)`
correlate_study <- function(ours,theirs) {
  # merge ours with theirs
  corr <- ours %>% mutate(merge_id = paste(Protein,originalRUN,sep="")) %>%
    select(LogIntensities,Protein,merge_id) %>% 
    inner_join(
      theirs %>% 
        mutate(merge_id = paste(Protein,originalRUN,sep="")) %>%
        select(LogIntensities,Protein,originalRUN,merge_id), 
      by="merge_id"
    ) %>% drop_na() %>%
    select(-merge_id) %>%
    rename(run = originalRUN) %>%
    filter(is.finite(LogIntensities.x)) %>%
    filter(is.finite(LogIntensities.y))
  
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

# Example
# ours <- ours_PXD004873()
# theirs <- theirs_PXD004873()
# p_corr_ours <- correlate_techreps(ours)
# p_corr_vs <- correlate_study(ours,theirs)
