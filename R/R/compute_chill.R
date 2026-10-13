# 09/22

library(dplyr)
library(chillR)

compute_daily_chill_heat = function(hour_temps_one_loc) {

  hourtemps_clean = hour_temps_one_loc %>% dplyr::select(-any_of('Location'))

  results = chillR::daily_chill(
    hourtemps = as.data.frame(hour_temps_one_loc),
    running_mean = 1,
    models = list(Chill_Portions = chillR::Dynamic_Model, GDH = chillR::GDH )
  )
  as_tibble(results$daily_chill)
}


add_cumm_chill_heat = function(daily_chill_dta, season) {
  anchor = season_anchor_date(season)

  date_col = as.Date(paste(daily_chill_dta$Year, daily_chill_dta$Month, daily_chill_dta$Day, sep = '-'))

  daily_chill_dta %>%
    mutate(obs_date = date_col) %>%
    filter(obs_date >= anchor) %>%
    arrange(obs_date) %>%
    mutate(chill_portions_cumm = cumsum(Chill_Portions), gdh_cumm = cumsum(GDH))
}


run_window_pls = function(daily_chill_dta, bio_dta, split_month = 8, ...){

  daily_chill_obj = list(object_type = 'daily_chill', daily_chill = as.data.frame(daily_chill_dta))

  fit = chillR::PLS_chill_force(
    daily_chill_obj = daily_chill_obj,
    bio_data_frame = as.data.frame(bio_dta),
    split_month = split_month,
    chill_models = 'Chill_Portions',
    heat_models = 'GDH',
    ...
  )

  pls_summary = fit$Chill_Portions$GDH$PLS_summary

  if (is.null(pls_summary)) {
    stop(
      "run_critical_window_pls(): chillR::PLS_chill_force()'s output did ", 'not have the expected result$Chill_Portions$GDH$PLS_summary element. '
    )
  }

  as_tibble(pls_summary)
}