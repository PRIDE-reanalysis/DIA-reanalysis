require(tidyverse)
require(stringr)

run_median_norm <- function(data, ...) {
  data %>% group_by(originalRUN) %>% 
    mutate(run_sd = sd(LogIntensities), run_mdn = median(LogIntensities), run_count = n()) %>%
    mutate(NormLogIntensities = LogIntensities-run_mdn) %>%
    ungroup() %>%
    dplyr::select(-run_sd, -run_mdn, -run_count)
}

run_quantile_norm <- function(data,... ){
  matrix <- data %>% dplyr::select(Protein, originalRUN, LogIntensities) %>%
    tidyr::spread(originalRUN, LogIntensities) %>% 
    column_to_rownames("Protein") %>%
  drop_na()  # alt. 0 imputation

  # substitution function (according to rank)
  index_lookup <- function(i, m){
    return(m[i])
  }
  
  data_rank <- map_df(matrix,rank,ties.method="average")  # 1. rank
  data_sorted <- map_df(matrix,sort)  # 2. sort
  data_mean <- rowMeans(data_sorted)  # 3. mean
  data_norm <- map_df(data_rank, index_lookup, m=data_mean)  # 4. substitute
  
  final <- data_norm %>% #bind_cols(row.names(matrix),.) %>% rename(Protein = ...1) %>% 
    bind_cols(matrix %>% rownames_to_column("Protein") %>% dplyr::select(Protein), . ) %>%
    tidyr::gather(key= "originalRUN", value= "NormLogIntensities", -Protein)
  
  return(data %>% # dplyr::select(-NormLogIntensities) %>%  # in case 
           dplyr::mutate(mergecol = paste(Protein,originalRUN,sep="_")) %>% 
           dplyr::inner_join(final %>% 
                      dplyr::mutate(mergecol= paste(Protein,originalRUN, sep="_")) %>%
                      dplyr::select(-originalRUN,-Protein), 
                    by= "mergecol") %>%
           dplyr::select(-mergecol)
  )
}

run_minmax_norm <- function(data, ...) {
  data %>% group_by(originalRUN) %>% 
    dplyr::mutate(run_min = min(LogIntensities), run_max = max(LogIntensities), run_count = n()) %>%
    dplyr::mutate(NormLogIntensities = (LogIntensities-run_min)/(run_max-run_min)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-run_min, -run_max, -run_count)
}

run_zscore_norm <- function(data, ...) {
  # normalisation on log values
  # https://www.researchgate.net/post/log_transformation_and_standardization_which_should_come_first
  data %>% group_by(originalRUN) %>% 
    dplyr::mutate(run_sd = sd(LogIntensities), run_mean = mean(LogIntensities), run_count = n()) %>%
    dplyr::mutate(NormLogIntensities = (LogIntensities-run_mean)/run_sd) %>%
    dplyr::ungroup() %>%
    dplyr::select(-run_sd, -run_mean, -run_count)
}

group_consistency <- function(data, ...){
  data %>%  
    dplyr::group_by(Protein,GROUP_ORIGINAL) %>% 
    dplyr::mutate( sumNumMeasuredFeature = sum(NumMeasuredFeature), sumNumImputedFeature = sum(NumImputedFeature) ) %>% 
    dplyr::mutate( feature_missingrate_per_group = (100/(sumNumMeasuredFeature+sumNumImputedFeature)) * sumNumImputedFeature ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(feature_missingrate_per_run = MissingPercentage * 100) %>%
    dplyr::select(-MissingPercentage)
}

clean_protein_names <- function(data, ...) {
  data %>%
    dplyr::mutate(Protein = str_replace(Protein, "^\\d+/", "")) %>% 
    dplyr::mutate(Protein = str_replace(Protein, "/$", "")) %>%
    dplyr::filter(Protein != "")
}

process_msstats_rda <- function(rda_obj, annot_obj, groups, norm="median") {
  if(missing(groups)) {
    groups <- as.character((rda_obj$RunlevelData %>% dplyr::distinct(GROUP_ORIGINAL))$GROUP_ORIGINAL)
  }
  return(our<-goldstandard.proposed$RunlevelData %>% 
           group_consistency(annot) %>%
           clean_protein_names() %>% 
           dplyr::filter(GROUP_ORIGINAL %in% groups) %>%
           {if(norm == "median") run_median_norm(.) else .} %>%
           {if(norm == "quantile") run_quantile_norm(.) else .} %>%
           {if(norm == "zscore") run_zscore_norm(.) else .} %>%
           {if(norm == "minmax") run_minmax_norm(.) else .} %>%
           dplyr::full_join(annot, by = c("originalRUN" = "Run")) %>%
           dplyr::select(Protein,dplyr::matches("Intensities"),originalRUN,
                         GROUP_ORIGINAL,SUBJECT_ORIGINAL,
                         dplyr::matches("missingrate") ) %>%
           dplyr::ungroup())
}



