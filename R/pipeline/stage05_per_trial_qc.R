# 09/18

stage05_targets = list(
  tar_target(
    al_qc_environments,
    dplyr::bind_rows(
      al_yield_long_clean %>% dplyr::distinct(season, location),
      al_bloom_long_clean %>% dplyr::distinct(season, location)
    ) %>% dplyr::distinct(season, location)
  ),

  tar_target(
    al_per_trial_qc,
    al_qc_environments %>%
    dplyr::mutate(
      qc = purrr::pmap(list(season, location), function(season, location){
        yield_env = al_yield_long_clean[al_yield_long_clean$season == season & al_yield_long_clean$location == location, ]
        bloom_env = al_bloom_long_clean [al_bloom_long_clean$season == season & al_bloom_long_clean$location == location, ]

        build_per_trial_qc_dta(yield_env, bloom_env)$results
      })
    ) %>%
      tidyr::unnest(qc)
  )
)