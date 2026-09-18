

stage04_targets <- list(
  tar_target(
    al_bloom_environments,
    al_bloom_long_clean %>% dplyr::distinct(season, location)
  ),

  tar_target(
    al_bloom_doy,
    al_bloom_environments %>%
      dplyr::mutate(
        doy = purrr::pmap(list(season, location), function(season, location) {
          env_data = al_bloom_long_clean[
            al_bloom_long_clean$season == season & al_bloom_long_clean$location == location,
          ]
          interpolate_bloom_doy(env_data)
        })
      ) %>%
      tidyr::unnest(doy)
  ),


  tar_target(
    al_bloom_auc,
    al_bloom_environments %>%
      dplyr::mutate(
        auc = purrr::pmap(list(season, location), function(season, location) {
          env_data = al_bloom_long_clean[
            al_bloom_long_clean$season == season & al_bloom_long_clean$location == location,
          ]
          compute_bloom_auc(env_data)
        })
      ) %>%
      tidyr::unnest(auc)
  ),


  # next two not production, just keeping this here since I spent so much time making it work

  tar_target(
    al_bloom_doy_nlme,
    al_bloom_environments %>%
      dplyr::mutate(
        doy = purrr::pmap(list(season, location), function(season, location) {
          env_data = al_bloom_long_clean[
            al_bloom_long_clean$season == season & al_bloom_long_clean$location == location,
          ]
          harmonize_bloom_doy_nlme(env_data)
        })
      ) %>%
      tidyr::unnest(doy),
    error = "continue"
  ),


  tar_target(
    al_bloom_cv_validation,
    al_bloom_environments %>%
      dplyr::mutate(
        cv = purrr::pmap(list(season, location), function(season, location) {
          env_data = al_bloom_long_clean[
            al_bloom_long_clean$season == season & al_bloom_long_clean$location == location,
          ]
          leave_visits_out_cv(env_data)
        })
      ) %>%
      tidyr::unnest(cv),
    error = "continue"
  )
)
