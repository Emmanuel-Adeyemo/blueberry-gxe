# 09/24

library(dplyr)

prepare_fa_dta = function(subplot_dta, trait_col, qc_pass_envts = NULL) {

  dta = subplot_dta

  if (!is.null(qc_pass_envts)) {
    data <- data %>%
      inner_join(
        qc_pass_envts %>% distinct(.data$season, .data$location),
        by = c('season', 'location')
      )
  }

  dta %>% filter(!is.na(.data[[trait_col]])) %>%
    transmute(
      Genotype = as.character(Genotype),
      Env = paste(location, season),
      Block = as.character(BlockNo),
      y = .data[[trait_col]]
    )

}

prepare_nf_fa_dta = function(data, trait_col, genotype_col = 'genotype',
                                        block_col, location_col = 'location', season_col = 'season') {
  data %>%
    filter(!is.na(.data[[trait_col]])) %>%
    transmute(
      Genotype = as.character(.data[[genotype_col]]),
      Env = paste(.data[[location_col]], .data[[season_col]]),
      Block = as.character(.data[[block_col]]),
      y = .data[[trait_col]]
    )
}


# used within fa model calls
prep_and_scale_dta = function(fa_dta){
  fa_dta = fa_dta %>% mutate(Genotype = factor(Genotype), Env = factor(Env), Block = factor(Block))

  env_stats = fa_dta %>% group_by(Env) %>%
    summarize(mean_env = mean(y, na.rm = TRUE), sd_env = sd(y, na.rm = TRUE))

  fa_dta_scaled = fa_dta %>% left_join(env_stats, by = 'Env') %>%
    mutate(y_scaled = (y - mean_env) / sd_env)
}

