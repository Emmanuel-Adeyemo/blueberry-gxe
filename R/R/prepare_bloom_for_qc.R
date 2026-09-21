# 09/18

library(dplyr)


aggregate_bloom_to_subplot = function(bloom_dta, trait_cols) {
  bloom_dta %>%
    group_by(season, location, BlockNo, SubplotNo, Genotype) %>%
    summarise(across(all_of(trait_cols), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
}