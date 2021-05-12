library(tidyverse)
library(stringr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

source("../container/downstream/DIA_downstream_datacarpentry.R")
source("../container/postprocess/codify_study_customisations.R")

# This needs a nifty table that has all the numbers collected
protein_counts_all_projects <- read.delim("../R/protein_counts_all_projects.tsv")

# define order manually
PXD_order <- c("PXD004873",
               "PXD000672",
               "PXD004691",
               "PXD014943",
               "PXD003497",
               "PXD004589",
               "PXD014194",
               "PXD003539",
               "PXD001064",
               "PXD010912")

protein_counts_all_projects <- protein_counts_all_projects %>% 
  mutate(PXD = factor(PXD, levels=PXD_order %>% 
                        stringr::str_replace('\\*', ''))
  )

fdr <- ggplot(protein_counts_all_projects %>% 
                filter(!FDR == "0.5%") %>% 
                mutate_at(vars(FDR), list(factor)) , aes(x=FDR, y=Total.proteins, fill=FDR)) +
  geom_col(position = "dodge") +
  xlab("FDR level") + ylab("Total distinct proteins") + 
  theme_bw() + theme(legend.position = "none",
                     panel.grid.major.x = element_blank()) + scale_x_discrete(limits = rev) +
  facet_wrap(~PXD,scales = "free_x") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 10000, 8*250), minor_breaks = seq(0, 10000, 250)) + 
  scale_fill_brewer(palette="Dark2")


protein_comp <- read.delim("../R/protein_comparisons_where_available.tsv") %>% 
  #dplyr::rename(Reanalysis_runfiltered = Proteins..50..missing) %>%
  #dplyr::rename(Original_Data = Original...data) %>%
  dplyr::rename(`Original data\n (Publication)` = Original...after.filter) %>%
  dplyr::rename(`Reanalysis\n (unfiltered)` = Reanalysis.proteins ) %>%
  dplyr::rename(`Reanalysis\n (consistencyfilter)` = Reanalysis.proteins...50..per.group. ) %>%
  dplyr::select(PXD,`Reanalysis\n (unfiltered)`,`Reanalysis\n (consistencyfilter)`,`Original data\n (Publication)`) %>% 
  tidyr::gather(Source, count, c(`Reanalysis\n (unfiltered)`,`Reanalysis\n (consistencyfilter)`,`Original data\n (Publication)`)) %>%
  #mutate(count =  na_if(count,"N/A")) %>%
  #drop_na() %>%
  dplyr::mutate(count = ifelse(count == "N/A", 0, count)) %>% mutate(count = as.numeric(count))

# manually order the factors
protein_comp <- protein_comp %>% mutate(order = case_when(
  Source == "Reanalysis\n (unfiltered)" ~ 3,
  Source == "Reanalysis\n (consistencyfilter)" ~ 2,
  Source == "Original data\n (Publication)" ~ 1)) %>%
  mutate(Source = fct_reorder(as.factor(Source), desc(order))) %>%
  mutate(PXD = factor(PXD, levels=PXD_order))

discovery <- ggplot(protein_comp, aes(x=PXD, y=count, fill=Source)) +
  geom_bar(stat='identity', position='dodge') +
  ylab("Protein count") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                     panel.grid.major.x = element_blank(),
                     legend.key.height = unit(30, "pt") ) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 10000, 4*250), minor_breaks = seq(0, 10000, 250)) + 
  scale_fill_brewer(palette="Dark2")


fdr + discovery + plot_annotation(tag_levels = 'A')
# ggsave('dia_paper_numbers_panel.png', fdr + discovery + plot_annotation(tag_levels = 'A') ,
#        width = 160,
#        height = 106, 
#        units = "mm", dpi = 300
# )
