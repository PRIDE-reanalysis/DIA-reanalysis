source("../container/postprocess/DIA_postprocess_normfilter.R")
source("../container/postprocess/codify_study_customisations.R")


###
pxds <- c("PXD004873", "PXD000672", "PXD004691", "PXD014943", "PXD003497", "PXD004589", "PXD014194", "PXD003539", "PXD001064", "PXD010912")

for (n in pxds ){
  df <- get(paste("ours",n,sep="_"))()
  assign(paste("Norm_plot",n,sep="_"), norm_plots(df, "")) # fix headers for panel
}

for (n in pxds ){

  panel_normbox <- ( get(paste("Norm_plot",n,sep="_"))$box_raw_unfiltered | get(paste("Norm_plot",n,sep="_"))$box_raw_filtered ) / 
    ( get(paste("Norm_plot",n,sep="_"))$box_normalised_unfiltered | get(paste("Norm_plot",n,sep="_"))$box_normalised_filtered ) /  
    plot_annotation(tag_levels = 'A',
                    title = paste('Intensities overview for',n,sep=" "),
                    #subtitle = 'These 3 plots will reveal yet-untold secrets about our beloved data-set',
                    caption = 'Filter: feature_missingrate_per_group < 50%, Colours: by condition')
  panel_normdens <- ( get(paste("Norm_plot",n,sep="_"))$density_raw_unfiltered | get(paste("Norm_plot",n,sep="_"))$density_raw_filtered ) / 
    ( get(paste("Norm_plot",n,sep="_"))$density_normalised_unfiltered | get(paste("Norm_plot",n,sep="_"))$density_normalised_filtered ) /  
      plot_annotation(tag_levels = 'A',
                      title = paste('Intensities overview for',n,sep=" "),
                      #subtitle = 'These 3 plots will reveal yet-untold secrets about our beloved data-set',
                      caption = 'Filter: feature_missingrate_per_group < 50%, Colours: by condition')

  # panel_normbox
  # panel_normdens
  
}


    