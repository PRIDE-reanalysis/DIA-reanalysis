library(tidyverse)
library(broom)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

# Welch t-test, but not for each protein
# t.test(NormLogIntensities ~ GROUP_ORIGINAL, data %>% filter(GROUP_ORIGINAL %in% c("eDLBCL","PCNSL") ))
# from dataframe: 
# example <- data_select %>% filter(Protein == "A0FGR8") 
# t_test_results <- t.test(NormLogIntensities ~ GROUP_ORIGINAL,
#                                                          data = example,
#                                                          alternative = "two.sided",
#                                                          mu = 0,
#                                                          paired = FALSE,
#                                                          var.equal = FALSE,
#                                                          conf.level = 0.95)

# `data` df argument: long form w/ cols: log.fc,log.pval
# name: argument should contain dataset name and 'left condition'-'right condition' names
# fc_threshold: argument abs. number, default 1 - assumes log2 fc
# *_colour: arguments col. hex or name, left will show up on the right and vice versa, because 
#   left/right refers to contrast order (and thus name order)
fc_plot <- function(data, name, left_colour, right_colour, center_colour, fc_threshold, ...) {
  
  if(missing(left_colour) | missing(right_colour)) {
    left_colour <- brewer.pal(n = 8, name = "Dark2")[1]
    right_colour <- brewer.pal(n = 8, name = "Dark2")[3]
  }
  if(missing(center_colour)) {
    center_colour <- "black"
  }
  if(missing(fc_threshold)) {
    fc_threshold <- 1
  }

  # Plot with points coloured according to the threshold
  #   threshold (-) right, neg. FC, overexpr. in 'right' side of contrast
  #   threshold (+) left, pos. FC, overexpr. in 'left' side of contrast
  #   threshold (<>) center, FC of pval not over threshold 
  plot <- ggplot( data %>%
                    dplyr::mutate(threshold = if_else(log.pval < 1.3 , "center" , 
                                                      if_else(log.fc <= (-1*fc_threshold), "right", 
                                                              if_else(log.fc >= fc_threshold, "left" , "center")
                                                      )
                                              )
                    ), 
                  aes(log.fc,log.pval, colour = threshold, size = 2)) +
    geom_point(alpha = 0.75, show.legend = FALSE) +
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + # threshold indicators, semi-transparent
    geom_vline(xintercept = fc_threshold, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = (-1*fc_threshold), linetype = 2, alpha = 0.5) +
    scale_x_continuous(breaks=c(0,seq(-5,5,2)), limits=c(-5.5, 5.5)) + 
    scale_y_continuous(expand = c(0,0), breaks=seq(1,100,2)) + 
    scale_colour_manual(values = c("left"= left_colour, "center"= center_colour, "right"= right_colour)) + # colour by threshold
    labs(title=name,x="log2 fold change", y = "-log10 p-value") +
    #geom_text_repel(data=data %>% filter(threshold == "I"), aes(label=Protein), show.legend = FALSE) +
    #theme(legend.position="none") + # Hide the legend
    #guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE) + # double-down on hiding the legend
    theme_classic() # Set the theme
  return(plot)
}

# `data` df annotations:
do_fc <-function(data, grp1, grp1_names, grp2, grp2_names){
  # first filter for proteins eligible for fold-change tests
  data_select <- data %>% 
    # consistency filter 
    dplyr::filter(feature_missingrate_per_group < 50) %>%
    # group contrast filter
    dplyr::filter(GROUP_ORIGINAL %in% c(grp1,grp2)) %>% 
    # also filter per protein to be at least detected in two runs per group
    dplyr::ungroup() %>%
    dplyr::group_by(Protein, GROUP_ORIGINAL) %>% 
    dplyr::filter(n() >1) %>%
    # and filter for proteins in both groups
    dplyr::ungroup() %>%
    dplyr::group_by(Protein) %>%
    dplyr::filter(n_distinct(GROUP_ORIGINAL) >1)
  
  # needs to be group_by(Protein)
  data_pvalues <- do(data_select, tidy(t.test(.$NormLogIntensities ~ .$GROUP_ORIGINAL,
                                              alternative = "two.sided",
                                              mu = 0,
                                              paired = FALSE,
                                              var.equal = FALSE,
                                              conf.level = 0.95
  ))) 
  # immediately adjust for multiple tests
  data_pvalues$p.adjust <- p.adjust(data_pvalues$p.value, method="BH")
  # convert to log.pval and get rid of the broom generated columns we don't need
  data_pvalues <- data_pvalues %>% 
    dplyr::mutate(log.pval = -1*log10(p.adjust)) %>% 
    dplyr::select(-starts_with("estim"),-statistic,-parameter) 
  
  # calculate the fold change
  data_fc <- data_select %>% 
    dplyr::group_by(Protein,GROUP_ORIGINAL) %>% 
    dplyr::mutate(group_mean = mean(NormLogIntensities)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Protein,GROUP_ORIGINAL,group_mean) %>%
    dplyr::distinct(Protein,GROUP_ORIGINAL,group_mean,.keep_all = FALSE) %>% 
    dplyr::group_by(Protein) %>% spread(key=GROUP_ORIGINAL, value=group_mean ) %>%
    tidyr::drop_na() %>% mutate(log.fc = (!!as.name(grp1))-(!!as.name(grp2))) %>%
    dplyr::inner_join(data_pvalues, by = "Protein")
  
  # combine data_fc and data_pvalues for plotting purposes
  data_plot <- data_fc %>%
    dplyr::select(Protein,log.fc,log.pval,p.adjust,p.value) %>% 
    dplyr::distinct(Protein, .keep_all = TRUE)
  
  return(data_plot)
}

# name: descriptive name for the study - This would be the PXD for us, usually
# contrast_list: list of vectors with string pairs - each string must correspond 
#   to a condition as used in the annotation given to MSstats (mind stray `\s`!)
# ! Also needs a function available to acquire a processed (MSstats) dataframe  
# ! with the study quant data; function name must `ours_<name>` and have arguments
# ! for the `groups`(for filtering) and `norm` (for median normalisation)
# ! see codify_study_customisations.R
# With datastructure constructed like this, you can loop over `names(pxds)` and 
# get the name argument (iterated) and the contrast list argument `pxds[[name]]`:
# pxds <- list("PXD014943" = list(c("eDLBCL", "PCNSL") , c("IVL", "eDLBCL")),
#              "PXD004691" = list(c("F-N", "F-T"), c("P-N","P-T")),
#              "PXD000672" = list(c("normal_ccRCC", "ccRCC"), c("ccRCC", "pRCC")) )
calc_contrasts <- function(name, contrast_list){
  return_list <- list()
  print(name)
  for (c in contrast_list ) {
    grp1 <- c[[1]]
    grp2 <- c[[2]]
    data <- get(paste("ours_",name,sep=""))(groups=c(grp1,grp2),norm="median")
    
    print(paste(grp1,grp2,sep="-"))
    
    grp1_names <- data %>% filter(GROUP_ORIGINAL == grp1) %>% distinct(originalRUN) %>% pull(originalRUN) %>% as.character()
    grp2_names <- data %>% filter(GROUP_ORIGINAL == grp2) %>% distinct(originalRUN) %>% pull(originalRUN) %>% as.character()
    
    print( paste(length(grp1_names),"vs",length(grp2_names), sep = " ") )
    
    plot_data <- do_fc(data, grp1 = grp1, grp1_names = grp1_names, grp2 = grp2, grp2_names = grp2_names)
    return_list[[as.name(paste("contrast",grp1,grp2,sep="-"))]] <- plot_data
    print(paste(nrow(plot_data), "Proteins considered",sep = " "))
    
    plt <- fc_plot(plot_data, paste(name, paste(grp1,grp2,sep="-"), sep=" "))
    return_list[[as.name(paste(grp1,grp2,sep="-"))]] <- plt
  }    
  return(return_list)
}

calc_ext_contrasts <- function(name, contrast_list){
  return_list <- list()
  print(name)
  for (c in contrast_list ) {
    grp1 <- c[[1]]
    grp2 <- c[[2]]
    data <- get(paste("theirs_",name,sep=""))()
    
    print(paste(grp1,grp2,sep="-"))
    
    data <- data %>% dplyr::filter(GROUP_ORIGINAL %in% c(grp1, grp2)) %>% 
      mutate(feature_missingrate_per_group = 1) %>%  # assume already consistency filtered
      mutate(NormLogIntensities = LogIntensities) %>% # assume already normalised
      dplyr::filter(GROUP_ORIGINAL %in% c(grp1, grp2))
    
    grp1_names <- data %>% filter(GROUP_ORIGINAL == grp1) %>% distinct(originalRUN) %>% pull(originalRUN) %>% as.character()
    grp2_names <- data %>% filter(GROUP_ORIGINAL == grp2) %>% distinct(originalRUN) %>% pull(originalRUN) %>% as.character()
    
    print( paste(length(grp1_names),"vs",length(grp2_names), sep = " ") )
    
    plot_data <- do_fc(data, grp1 = grp1, grp1_names = grp1_names, grp2 = grp2, grp2_names = grp2_names)
    return_list[[as.name(paste("contrast",grp1,grp2,sep="-"))]] <- plot_data
    print(paste(nrow(plot_data), "Proteins considered",sep = " "))
    
    plt <- fc_plot(plot_data, paste(name, paste(grp1,grp2,sep="-"), sep=" "))
    return_list[[as.name(paste(grp1,grp2,sep="-"))]] <- plt
  }    
  return(return_list)
}
