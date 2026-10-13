# 09/18

library(dplyr)


aggregate_bloom_to_subplot = function(bloom_dta, trait_cols) {
  bloom_dta %>%
    group_by(season, location, BlockNo, SubplotNo, Genotype) %>%
    summarise(across(all_of(trait_cols), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
}


aggregate_al_bloom_for_qc = function(al_bloom_long_clean) {
  doy_plant = interpolate_bloom_doy(al_bloom_long_clean)
  doy_subplot = aggregate_bloom_to_subplot(
    doy_plant, c('bloom5', 'bloom50', 'bloom100', 'bloom_duration')
  )
  doy_subplot %>%
    mutate(
      bloom5   = days_since_season_anchor(bloom5, season),
      bloom50  = days_since_season_anchor(bloom50, season),
      bloom100 = days_since_season_anchor(bloom100, season),
      bloom_duration = as.numeric(bloom_duration)
    )
}