# 09/24

library(dplyr)

prepare_fa_dta = function(subplot_dta, trait_col, qc_pass_envts = NULL, extra_cols = NULL) {

  dta = subplot_dta

  if (!is.null(qc_pass_envts)) {
    data <- data %>%
      inner_join(
        qc_pass_envts %>% distinct(.data$season, .data$location),
        by = c('season', 'location')
      )
  }

  dta = dta %>% filter(!is.na(.data[[trait_col]]))

  out = dta %>%
    transmute(
      Genotype = as.character(Genotype),
      Env = paste(location, season),
      Block = as.character(BlockNo),
      y = .data[[trait_col]]
    )

  # 10/26 to pass Age
  if (!is.null(extra_cols)) {
    missing_cols = setdiff(extra_cols, names(dta))
    if (length(missing_cols) > 0) {
      stop(
        "prepare_fa_dta(): extra_cols not found in subplot_dta: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    out = dplyr::bind_cols(out, dta %>% dplyr::select(dplyr::all_of(extra_cols)))
  }

  out

}

prepare_nf_fa_dta = function(dta, trait_col, genotype_col = 'genotype',
                                        block_col, location_col = 'location', season_col = 'season', extra_cols = NULL) {
  dta = dta %>% filter( !is.na(.data[[trait_col]]))

  out = dta %>%
    # filter(!is.na(.data[[trait_col]])) %>%
    transmute(
      Genotype = as.character(.data[[genotype_col]]),
      Env = paste(.data[[location_col]], .data[[season_col]]),
      Block = as.character(.data[[block_col]]),
      y = .data[[trait_col]]
    )
  
  if (!is.null(extra_cols)) {
    # print(names(dta))
    # print(dim(out))
    # print(dim(dta %>% filter( !is.na(.data[[trait_col]])) %>% dplyr::select(dplyr::all_of(extra_cols))))
    missing_cols = setdiff(extra_cols, names(dta))
    if (length(missing_cols) > 0) {
      stop(
        "prepare_fa_dta(): extra_cols not found in subplot_dta: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    out = dplyr::bind_cols(out, dta %>% dplyr::select(dplyr::all_of(extra_cols)))
  }

  out
}


# 10/16 quality traits 
aggregate_nc_quality_for_fa = function(trait_dta, value_col, block_col = 'subplot', extra_cols = NULL) {

  trait_dta = trait_dta %>% filter(!is.na(.data[[block_col]]))

  group_cols = c('season', 'location', 'plot', block_col, 'genotype')

  dta = trait_dta %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(!!value_col := mean(.data[[value_col]], na.rm = TRUE), .groups = 'drop') %>%
    mutate(!!value_col := ifelse(is.nan(.data[[value_col]]), NA_real_, .data[[value_col]]))

  dta
}


# used within fa model calls
prep_and_scale_dta = function(fa_dta){
  fa_dta = fa_dta %>%
    mutate(Genotype = factor(Genotype), Env = factor(Env), Block = factor(Block)) %>%
    arrange(Env)
 
  env_stats = fa_dta %>% group_by(Env) %>%
    summarize(mean_env = mean(y, na.rm = TRUE), sd_env = sd(y, na.rm = TRUE))
 
  fa_dta_scaled = fa_dta %>% left_join(env_stats, by = 'Env') %>%
    mutate(y_scaled = (y - mean_env) / sd_env)
}