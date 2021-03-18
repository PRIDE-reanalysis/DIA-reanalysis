library(tidyverse)
library(broom)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

# `data` df argument: long form w/ cols: LogIntensities, NormLogIntensities (median), 
#   feature_missingrate_per_group, originalRUN, `GROUP_ORIGINAL (see colourcol)`
#   i.e. as it comes from `codify_study_customisations.R`
# name: argument should contain dataset name
# colourcol: argument should contain the column name that is used for colour coding (defaults to GROUP_ORIGINAL)
norm_plots <- function(data, name, colourcol, ...) {
  if(missing(colourcol)) {
    colcol <- "GROUP_ORIGINAL"
  }
  
  brf <- ggplot(df %>% 
                  filter(feature_missingrate_per_group < 50) ,   
                aes(x = LogIntensities , y = originalRUN, fill = get(colcol) ) ) +
    geom_boxplot(outlier.size = 0.5) + 
    theme(legend.position = "none") + 
    labs(title=paste("Raw intensities filtered", name, sep = " "),x="LogIntensity", y = "Runs")
  bnf <- ggplot(data %>% 
                  filter(feature_missingrate_per_group < 50) ,  
                aes(x = NormLogIntensities , y = originalRUN , fill = get(colcol)) ) +
    geom_boxplot(outlier.size = 0.5) + 
    theme(legend.position = "none") + 
    labs(title=paste("Median normalised filtered", name, sep = " "),x="LogIntensity", y = "Runs")
  
  bru <- ggplot(data, 
                aes(x = LogIntensities , y = originalRUN , fill = get(colcol)) ) +
    geom_boxplot(outlier.size = 0.5) + 
    theme(legend.position = "none") + 
    labs(title=paste("Raw intensities unfiltered", name, sep = " "),x="LogIntensity", y = "Runs")
  bnu <- ggplot(data, 
                aes(x = NormLogIntensities , y = originalRUN , fill = get(colcol)) ) +
    geom_boxplot(outlier.size = 0.5) + 
    theme(legend.position = "none") + 
    labs(title=paste("Median normalised unfiltered", name, sep = " "),x="LogIntensity", y = "Runs")
  
  
  drf <- ggplot(data %>% 
                  filter(feature_missingrate_per_group < 50) ,  
                aes(x = LogIntensities , color = get(colcol)) ) +
    geom_density() + 
    theme(legend.position = "none") + 
    labs(title=paste("Raw intensities filtered", name, sep = " "), x="LogIntensity", y="Density")
  dnf <- ggplot(data %>% 
                  filter(feature_missingrate_per_group < 50) , 
                aes(x = NormLogIntensities , color = get(colcol)) ) +
    geom_density() +
    theme(legend.position = "none") + 
    labs(title=paste("Median normalised filtered", name, sep = " "), x="LogIntensity", y="Density")
  
  dru <- ggplot(data, 
                aes(x = LogIntensities , color = get(colcol)) ) +
    geom_density() + 
    theme(legend.position = "none") + 
    labs(title=paste("Raw intensities unfiltered", name, sep = " "), x="LogIntensity", y="Density")
  dnu <- ggplot(data, 
                aes(x = NormLogIntensities , color = get(colcol)) ) +
    geom_density() +
    theme(legend.position = "none") + 
    labs(title=paste("Median normalised unfiltered", name, sep = " "), x="LogIntensity", y="Density")
  
  return_list <- list("box_normalised_filtered" = bnf,
                      "box_normalised_unfiltered" = bnu,
                      "box_raw_filtered" = brf,
                      "box_raw_unfiltered" = bru,
                      "density_normalised_filtered" = dnf,
                      "density_normalised_unfiltered" = dnu,
                      "density_raw_filtered" = drf,
                      "density_raw_unfiltered" = dru)
  return(return_list) 
}