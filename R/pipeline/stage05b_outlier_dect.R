# 10/20
# Outlier detection on per-trial 


stage05b_targets = list(

  tar_target(
    al_yield_trial_dta,
    qc_envts %>% filter(location %in% c('BW', 'EV')) %>%
      dplyr::mutate(
        trial_dta = purrr::pmap(list(season, location), function(season, location) {
          yield_env = al_yield_long_clean[al_yield_long_clean$season == season & al_yield_long_clean$location == location, ]
          bloom_env = al_bloom_long_clean[al_bloom_long_clean$season == season & al_bloom_long_clean$location == location, ]
          build_per_trial_qc_dta(yield_env, bloom_env)$dta$TotalYield
        })
      )
  ),

  tar_target(
    nc_yield_trial_dta,
    qc_envts %>% filter(location == 'NC') %>%
      dplyr::mutate(
        trial_dta = purrr::pmap(list(season, location), function(season, location) {
          yield_env = nc_yield_clean[nc_yield_clean$season == season & nc_yield_clean$location == location, ]
          bloom_env = nc_bloom_doy[nc_bloom_doy$season == season & nc_bloom_doy$location == location, ]
          build_per_trial_qc_dta_nc_fl(yield_env, bloom_env, block_col = 'subplot')$dta$TotalYield
        })
      )
  ),

  tar_target(
    fl_yield_trial_dta,
    qc_envts %>% filter(location == 'FL') %>%
      dplyr::mutate(
        trial_dta = purrr::pmap(list(season, location), function(season, location) {
          yield_env = fl_yield_clean[fl_yield_clean$season == season & fl_yield_clean$location == location, ]
          build_per_trial_qc_dta_nc_fl(yield_env, NULL, block_col = 'subplot')$dta$TotalYield
        })
      )
  ),

  tar_target(
    yield_trial_dta_all,
    dplyr::bind_rows(al_yield_trial_dta, nc_yield_trial_dta, fl_yield_trial_dta) %>%
      dplyr::mutate(Env = paste(location, season)) %>%
      dplyr::filter(!purrr::map_lgl(trial_dta, ~ is.null(.x) || nrow(.x) == 0))
  ),

  # named list keyed by Env -- what build_trait_blues_dta() expects
  tar_target(
    yield_trial_dta_named,
    setNames(yield_trial_dta_all$trial_dta, yield_trial_dta_all$Env)
  ),

  tar_target(
    yield_blues,
    build_trait_blues_dta(yield_trial_dta_named)
  ),

  tar_target(
    yield_outlier_fit,
    detect_blue_outliers_asreml(yield_blues, trt = 'yield_sqrt')
  ),

  tar_target(
    yield_outliers,
    yield_outlier_fit$outliers
  ),

  tar_target(
    yield_blues_clean,
    yield_outlier_fit$blues_clean
  )

)