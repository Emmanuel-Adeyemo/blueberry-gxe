# 09/18
# nc and fl added 10/03

stage05_targets = list(
  tar_target(
    qc_envts,
    dplyr::bind_rows(
      al_yield_long_clean %>% dplyr::distinct(season, location),
      al_bloom_long_clean %>% dplyr::distinct(season, location),
      nc_yield_clean %>% distinct(season, location),
      fl_yield_clean %>% distinct(season, location)
    ) %>% dplyr::distinct(season, location)
  ),

  tar_target(
    al_per_trial_qc,
    qc_envts %>% filter(location %in% c('BW', 'EV')) %>%
    dplyr::mutate(
      qc = purrr::pmap(list(season, location), function(season, location){
        yield_env = al_yield_long_clean[al_yield_long_clean$season == season & al_yield_long_clean$location == location, ]
        bloom_env = al_bloom_long_clean [al_bloom_long_clean$season == season & al_bloom_long_clean$location == location, ]

        build_per_trial_qc_dta(yield_env, bloom_env)$results
      })
    ) %>%
      tidyr::unnest(qc)
  ),

  # nc
  tar_target(
    nc_per_trial_qc,
    qc_envts %>% 
      filter(location == 'NC') %>%
      mutate(
        qc = purrr::pmap(list(season, location), function(season, location){
          yield_env = nc_yield_clean[nc_yield_clean$season == season & nc_yield_clean$location == location, ]
          bloom_env = nc_bloom_doy [nc_bloom_doy$season == season & nc_bloom_doy$location == location, ]

          build_per_trial_qc_dta_nc_fl(yield_env, bloom_env, block_col = 'subplot')$results
      })
    )%>%
      tidyr::unnest(qc)
  ),

  # fl
  tar_target(
    fl_per_trial_qc,
    qc_envts %>% 
      filter(location == 'FL') %>%
      mutate(
        qc = purrr::pmap(list(season, location), function(season, location){
          yield_env = fl_yield_clean[fl_yield_clean$season == season & fl_yield_clean$location == location, ]
          # fl bloom dta has no reps, model will fail.
          # bloom_env = fl_bloom_doy [fl_bloom_doy$season == season & fl_bloom_doy$location == location, ]

          build_per_trial_qc_dta_nc_fl(yield_env, NULL, block_col = 'site_detail')$results
      })
    )%>%
      tidyr::unnest(qc)
  ),

  tar_target(
    per_trial_qc,
    bind_rows(al_per_trial_qc, nc_per_trial_qc, fl_per_trial_qc)
  )
)