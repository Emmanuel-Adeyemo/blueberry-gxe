# 09/12

# future bloom date should not have a regression in bloom progression
# using Sushan's carry forward approach ie change error value to last known correct value

library(dplyr)


#' @description corrects non progressing bloom-percent readings ie later observation lower than previous. Sushan's carry-forward approach: any lower obs is replaced with
#'   the last known (highest) correct value, and flagged via `carried_forward`
#' @param bloom_dta long-format bloom obs data, e.g. `al_bloom_long` output of `ingest_all_al_dta`
#' @output `bloom_dta` df
carry_forward_bloom_correction = function(bloom_dta){

  bloom_dta %>%
    arrange(season, location, PlantID, obs_date, Age) %>%
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
#' @description identifies and excludes plants that never reached any bloom percent ie `max(bloom_percent) = 0`
#' @param bloom_dta long format bloom obs e.g output of `carry_forward_bloom_correction()`
#' @output list with clean bloom df and excluded df
exclude_never_bloomed = function(bloom_dta) {
  plant_max = bloom_dta %>%
    group_by(PlantID, season) %>%
    summarise(max_bloom = max(bloom_percent), .groups = 'drop')

  dead_plants = plant_max %>% filter(max_bloom == 0)

  filtered = bloom_dta %>% anti_join(dead_plants, by = c('PlantID', 'season'))

  list(
    bloom_dta = filtered,
    excluded_plants = dead_plants %>% dplyr::select(PlantID, season)
  )
}