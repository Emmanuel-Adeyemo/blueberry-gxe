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

        # 10/16 quality traits 
        plant_height_env = al_plant_height_long_clean[al_plant_height_long_clean$season == season & al_plant_height_long_clean$location == location, ]
        avg_berry_weight_env = al_avg_berry_weight_long_clean[al_avg_berry_weight_long_clean$season == season & al_avg_berry_weight_long_clean$location == location, ]
        firmness_env = al_firmness_long_clean[al_firmness_long_clean$season == season & al_firmness_long_clean$location == location, ]
        brix_env = al_brix_long_clean[al_brix_long_clean$season == season & al_brix_long_clean$location == location, ]
        ta_env = al_ta_long_clean[al_ta_long_clean$season == season & al_ta_long_clean$location == location, ]

        build_per_trial_qc_dta(
          yield_env, bloom_env,
          plant_height_dta = plant_height_env,
          avg_berry_weight_dta = avg_berry_weight_env,
          firmness_dta = firmness_env,
          brix_dta = brix_env,
          ta_dta = ta_env
        )$results
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

          # 10/16 quality traits 
          avg_berry_weight_env = nc_avg_berry_weight_clean[nc_avg_berry_weight_clean$season == season & nc_avg_berry_weight_clean$location == location, ]
          firmness_env = nc_firmness_clean[nc_firmness_clean$season == season & nc_firmness_clean$location == location, ]
          brix_env = nc_brix_clean[nc_brix_clean$season == season & nc_brix_clean$location == location, ]
          ta_env = nc_ta_clean[nc_ta_clean$season == season & nc_ta_clean$location == location, ]

          build_per_trial_qc_dta_nc_fl(
            yield_env, bloom_env, block_col = 'subplot',
            avg_berry_weight_dta = avg_berry_weight_env,
            firmness_dta = firmness_env,
            brix_dta = brix_env,
            ta_dta = ta_env
          )$results
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

          build_per_trial_qc_dta_nc_fl(yield_env, NULL, block_col = 'subplot')$results
      })
    )%>%
      tidyr::unnest(qc)
  ),

  tar_target(
    per_trial_qc,
    bind_rows(al_per_trial_qc, nc_per_trial_qc, fl_per_trial_qc)
  )
)