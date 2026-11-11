# 09/23, reworked 10/29 to use BW/EV station data instead of NASA power
#


stage06_targets <- list(
  tar_target(ev_weather_csv_path, "data/raw/weather/WeatherEV_2022_2026.csv", format = "file"),
  tar_target(bw_weather_csv_path, "data/raw/weather/WeatherBW_2022_2026.csv", format = "file"),
  tar_target(nc_weather_xlsx_path, "data/raw/weather/NCHCRS temp precip hourly 23-25.xlsx", format = "file"),

  tar_target(
    ev_weather_prepared,
    fill_hourly_gaps(prepare_station_weather_for_chillr(ev_weather_csv_path, "EV"), "EV")
  ),
  tar_target(
    bw_weather_prepared,
    fill_hourly_gaps(prepare_station_weather_for_chillr(bw_weather_csv_path, "BW"), "BW")
  ),
  tar_target(
    nc_weather_prepared,
    prepare_nc_weather_for_chillr(nc_weather_xlsx_path)
  ),

 
  tar_target(
    weather_gap_report,
    dplyr::bind_rows(
      summarize_gap_overlap(ev_weather_prepared, "EV", 2023:2027),
      summarize_gap_overlap(bw_weather_prepared, "BW", 2023:2027)
    )
  ),


  tar_target(ev_daily_chill_heat, compute_daily_chill_heat(ev_weather_prepared)),
  tar_target(bw_daily_chill_heat, compute_daily_chill_heat(bw_weather_prepared)),
  tar_target(nc_daily_chill_heat, compute_daily_chill_heat(nc_weather_prepared)),

  
  tar_target(
    chill_heat_environments,
    dplyr::bind_rows(
      get_available_seasons(ev_weather_prepared, "EV", 2023:2027),
      get_available_seasons(bw_weather_prepared, "BW", 2023:2027),
      get_available_seasons(nc_weather_prepared, "NC", 2023:2027)
    )
  ),

  tar_target(
    fixed_window_covariates,
    chill_heat_environments %>%
      dplyr::rowwise() %>%
      dplyr::mutate(cov = list({
        daily = switch(location, EV = ev_daily_chill_heat, BW = bw_daily_chill_heat, NC = nc_daily_chill_heat)
        compute_fixed_window_covariate(daily, season)
      })) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(cov, names_repair = "unique") %>%
      dplyr::rename(season = season...2) %>%
      dplyr::select(-season...3)
  ),

  tar_target(
    frost_scores,
    chill_heat_environments %>%
      dplyr::rowwise() %>%
      dplyr::mutate(cov = list({
        weather = switch(location, EV = ev_weather_prepared, BW = bw_weather_prepared, NC = nc_weather_prepared)
        compute_frost_score(weather, season)
      })) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(cov, names_repair = "unique") %>%
      dplyr::rename(season = season...2) %>%
      dplyr::select(-season...3)
  ),

  # actual environmental clustering input
  tar_target(
    env_clustering_covariates,
    fixed_window_covariates %>%
      dplyr::inner_join(frost_scores, by = c("location", "season")) %>%
      dplyr::left_join(weather_gap_report, by = c("location", "season"))
  ),

  tar_target(envt_cluster, plot_envt_cluster(env_clustering_covariates))
)