#!/usr/bin/env Rscript
library(optparse)
library(MSstats)
library(tidyverse)
library(stringr)
library(data.table)

tally_protpepirt <- function(rda, annot) {
  irt_peptides <- nrow(rda$ProcessedData %>% 
                         dplyr::filter(PROTEIN=="1/iRT_protein") %>% 
                         dplyr::mutate(irtp =  gsub("_.*","",PEPTIDE) ) %>% dplyr::distinct(irtp)
  )
  prot_all <- length(levels(rda$RunlevelData$Protein))
  
  prot_50 <- nrow(rda$RunlevelData %>% 
                    dplyr::filter(!more50missing == TRUE) %>% 
                    dplyr::distinct(Protein)
  )
  
  table_procdat <- data.frame(
    Analysis = c("Total num of proteins", 
                 "Total num of proteins with <50% missing values",
                 "Total num of iRT peptides",
                 "Total num of peptides"),
    measure = c(prot_all,
                prot_50,
                irt_peptides,
                length(levels(rda$ProcessedData$PEPTIDE)))
  )
  return(table_procdat)
  
}

option_list = list(
  make_option(c("-r", "--rda"), type="character", default=NULL, 
              help="full path to the MSstats file", metavar="character"),
  make_option(c("-a", "--ann"), type="character", default=NULL, 
              help="full path to the annotation file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="full path to the output tsv file to be written", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$rda) || is.null(opt$ann) ){
  print_help(opt_parser)
  stop("Required input arguments must be supplied.", call.=FALSE)
}

# This script assumes the MSstats dataframe object is called goldstandard.proposed
# ./CMDL_MSstats_rda_count.R --rda <...> --ann <...> --out <...>  and dont forget to chmod u+x

annot <- read.delim(opt$ann)
load(opt$rda)

procdat <- tally_protpepirt(goldstandard.proposed, annot)   

prot_group <- nrow(goldstandard.proposed$RunlevelData %>% 
                     dplyr::group_by(Protein,GROUP_ORIGINAL) %>% 
                     dplyr::mutate( sumNumMeasuredFeature = sum(NumMeasuredFeature), sumNumImputedFeature = sum(NumImputedFeature) ) %>% 
                     dplyr::mutate( feature_missingrate_per_group = (100/(sumNumMeasuredFeature+sumNumImputedFeature)) * sumNumImputedFeature ) %>%
                     dplyr::ungroup() %>%
                     dplyr::mutate(feature_missingrate_per_run = MissingPercentage * 100) %>%
                     dplyr::select(-MissingPercentage) %>%
                     dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
                     dplyr::ungroup() %>% 
                     dplyr::filter(feature_missingrate_per_group < 50) %>% 
                     dplyr::distinct(Protein))

tab <- rbind(procdat %>% dplyr::mutate(Analysis = as.character(Analysis)),
      c("Total num of Proteins <50% features missing in group", prot_group))

write.table(tab, opt$out, quote = FALSE, sep = "\t")

# Analysis measure
# 1                                Total num of proteins    3530
# 2       Total num of proteins with <50% missing values    2134
# 3                            Total num of iRT peptides       1
# 4                                Total num of peptides   22231
# 5 Total num of Proteins <50% features missing in group    3392
# Analysis measure
# 1                                Total num of proteins    4053
# 2       Total num of proteins with <50% missing values    2128
# 3                            Total num of iRT peptides       1
# 4                                Total num of peptides   31888
# 5 Total num of Proteins <50% features missing in group    3962
# Analysis measure
# 1                                Total num of proteins    2872
# 2       Total num of proteins with <50% missing values    1513
# 3                            Total num of iRT peptides       1
# 4                                Total num of peptides   17943
# 5 Total num of Proteins <50% features missing in group    2686
# Analysis measure
# 1                                Total num of proteins    5946
# 2       Total num of proteins with <50% missing values    3488
# 3                            Total num of iRT peptides       1
# 4                                Total num of peptides   45986
# 5 Total num of Proteins <50% features missing in group    5568
# Analysis measure
# 1                                Total num of proteins    2754
# 2       Total num of proteins with <50% missing values    1799
# 3                            Total num of iRT peptides       1
# 4                                Total num of peptides   20296
# 5 Total num of Proteins <50% features missing in group    2704
# Analysis measure
# 1                                Total num of proteins    3703
# 2       Total num of proteins with <50% missing values    1891
# 3                            Total num of iRT peptides       1
# 4                                Total num of peptides   31727
# 5 Total num of Proteins <50% features missing in group    3298
# Analysis measure
# 1                                Total num of proteins    2239
# 2       Total num of proteins with <50% missing values     881
# 3                            Total num of iRT peptides       1
# 4                                Total num of peptides   11593
# 5 Total num of Proteins <50% features missing in group    2145
      