# 09/22

# season = 2023, it spans from sept 2022 to ~June 2023, based on my definition
# will mostly likely cause error and misalign phenology. actual bloom calendar date is needed


library(dplyr)

build_pls_dta = function(bloom_doy_dta, target_col = 'bloom5') {
  bloom_doy_dta %>%
    group_by(season) %>%
    summarize(bloom_date = mean(.data[[target_col]], na.rm = T), .groups = 'drop') %>%
    filter(!is.na(bloom_date)) %>%
    transmute(
      Year = as.integer(format(bloom_date, '%Y')),
      pheno = as.integer(format(bloom_date, '%j'))
    )
}


build_bloom_pheno_by_season = function(bloom_doy_data, target_col = "bloom5", location) {
  bloom_doy_data %>%
    group_by(.data$season) %>%
    summarise(bloom_date = mean(.data[[target_col]], na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(.data$bloom_date)) %>%
    rowwise() %>%
    mutate(pheno = as.numeric(.data$bloom_date - season_anchor_date(.data$season))) %>%
    ungroup() %>%
    transmute(location = location, season = .data$season, pheno = .data$pheno)
}