# 09/12

# future bloom date should not have a regression in bloom proression
# using Sushan's carry forward approach ie change error value to last known correct value

library(dplyr)

carry_forward_bloom_correction = function(bloom_dta){

  bloom_dta %>%
    arrange(season, location, PlantID, obs_date) %>%
    group_by(season, location, PlantID) %>%
    mutate(
      bloom_percent_raw = bloom_percent,
      running_max = cummax(bloom_percent),
      carried_forward = running_max !=  bloom_percent_raw,
      bloom_percent = running_max
    ) %>%
      ungroup()
}


# never bloomed - assumed dead
exclude_never_bloomed = function(bloom_dta) {
  plant_max = bloom_dta %>%
    group_by(PlantID, season) %>%
    summarise(max_bloom = max(bloom_percent), .groups = 'drop')

  dead_plants = plant_max %>% filter(max_bloom == 0)

  filtered = bloom_dta %>% anti_join(dead_plants, by = c('PlantID', 'season'))

  list(
    bloom_dta = filtered,
    excluded_plants = dead_plants %>% select(PlantID, season)
  )
}