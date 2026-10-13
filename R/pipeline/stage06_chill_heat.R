# 09/23

stage06_targets <- list(
  tar_target(nasa_ev_csv_path, "data/raw/weather/NASA_EV_POWER_Point_Hourly_20220101_20260628_UTC.csv", format = "file"),
  tar_target(nasa_bw_csv_path, "data/raw/weather/NASA_BW_POWER_Point_Hourly_20220101_20260628_UTC.csv", format = "file"),

  tar_target(
    ev_weather_prepared,
    prepare_weather_for_chillr(readr::read_csv(nasa_ev_csv_path, show_col_types = FALSE), "EV")
  ),
  tar_target(
    bw_weather_prepared,
    prepare_weather_for_chillr(readr::read_csv(nasa_bw_csv_path, show_col_types = FALSE), "BW")
  ),

  tar_target(ev_daily_chill_heat, compute_daily_chill_heat(ev_weather_prepared)),
  tar_target(bw_daily_chill_heat, compute_daily_chill_heat(bw_weather_prepared)),

  tar_target(
    al_chill_heat_environments,
    dplyr::bind_rows(
      tibble::tibble(location = "EV", season = al_bloom_long_clean %>% dplyr::filter(location == "EV") %>% dplyr::distinct(season) %>% dplyr::pull(season)),
      tibble::tibble(location = "BW", season = al_bloom_long_clean %>% dplyr::filter(location == "BW") %>% dplyr::distinct(season) %>% dplyr::pull(season))
    )
  ),
  tar_target(
    al_cumulative_chill_heat,
    al_chill_heat_environments %>%
      dplyr::mutate(cumulative = purrr::map2(location, season, function(location, season) {
        daily <- if (location == "EV") ev_daily_chill_heat else bw_daily_chill_heat
        add_cumm_chill_heat(daily, season)
      })) %>%
      tidyr::unnest(cumulative)
  ),

  tar_target(
    ev_bio_data_frame,
    build_pls_dta(
      interpolate_bloom_doy(al_bloom_long_clean %>% dplyr::filter(location == "EV")),
      "bloom5"
    )
  ),
  tar_target(
    bw_bio_data_frame,
    build_pls_dta(
      interpolate_bloom_doy(al_bloom_long_clean %>% dplyr::filter(location == "BW")),
      "bloom5"
    )
  ),

  tar_target(ev_critical_window_pls, run_window_pls(ev_daily_chill_heat, ev_bio_data_frame, split_month = 8)),
  tar_target(bw_critical_window_pls, run_window_pls(bw_daily_chill_heat, bw_bio_data_frame, split_month = 8)),

  tar_target(ev_critical_window, extract_critical_window(smooth_critical_window(ev_critical_window_pls))),
  tar_target(bw_critical_window, extract_critical_window(smooth_critical_window(bw_critical_window_pls))),

  tar_target(
    al_fixed_window_covariates,
    dplyr::bind_rows(
      al_chill_heat_environments %>% dplyr::filter(location == "EV") %>%
        dplyr::mutate(season = as.integer(season)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(cov = list(compute_fixed_window_covariate(ev_daily_chill_heat, season))) %>%
        dplyr::ungroup() %>% tidyr::unnest(cov, names_repair = "unique"),
      al_chill_heat_environments %>% dplyr::filter(location == "BW") %>%
        dplyr::mutate(season = as.integer(season)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(cov = list(compute_fixed_window_covariate(bw_daily_chill_heat, season))) %>%
        dplyr::ungroup() %>% tidyr::unnest(cov, names_repair = "unique")
    ) %>% rename(season = season...2) %>% dplyr::select(-season...3)
  ),
  tar_target(
    al_discovered_window_covariates,
    dplyr::bind_rows(
      al_chill_heat_environments %>% dplyr::filter(location == "EV") %>%
        dplyr::mutate(season = as.integer(season)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(cov = list(compute_discovered_window_covariate(ev_daily_chill_heat, season, ev_critical_window))) %>%
        dplyr::ungroup() %>% tidyr::unnest(cov, names_repair = "unique"),
      al_chill_heat_environments %>% dplyr::filter(location == "BW") %>%
        dplyr::mutate(season = as.integer(season)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(cov = list(compute_discovered_window_covariate(bw_daily_chill_heat, season, bw_critical_window))) %>%
        dplyr::ungroup() %>% tidyr::unnest(cov, names_repair = "unique")
    ) %>% rename(season = season...2) %>% dplyr::select(-season...3)
  ),

  tar_target(
    al_bloom_pheno_by_season,
    dplyr::bind_rows(
      build_bloom_pheno_by_season(
        interpolate_bloom_doy(al_bloom_long_clean %>% dplyr::filter(location == "EV")), "bloom5", "EV"),
      build_bloom_pheno_by_season(
        interpolate_bloom_doy(al_bloom_long_clean %>% dplyr::filter(location == "BW")), "bloom5", "BW")
    )
  ),

  tar_target(al_loyo_fixed_chill, run_loyo_comparison(al_fixed_window_covariates, al_bloom_pheno_by_season, "Chill_Portions_sum")),
  tar_target(al_loyo_fixed_heat, run_loyo_comparison(al_fixed_window_covariates, al_bloom_pheno_by_season, "GDH_sum")),
  tar_target(al_loyo_discovered_chill, run_loyo_comparison(al_discovered_window_covariates, al_bloom_pheno_by_season, "Chill_Portions_sum")),
  tar_target(al_loyo_discovered_heat, run_loyo_comparison(al_discovered_window_covariates, al_bloom_pheno_by_season, "GDH_sum")),

  tar_target(
    al_loyo_summary,
    dplyr::bind_rows(
      summarise_loyo_comparison(al_loyo_fixed_chill, "Fixed: Chill"),
      summarise_loyo_comparison(al_loyo_fixed_heat, "Fixed: Heat"),
      summarise_loyo_comparison(al_loyo_discovered_chill, "PLS: Chill"),
      summarise_loyo_comparison(al_loyo_discovered_heat, "PLS: Heat")
    )
  ),

  # plots
  tar_target(ev_critical_window_plot, plot_critical_window_diagnostic(smooth_critical_window(ev_critical_window_pls), ev_critical_window)),
  tar_target(bw_critical_window_plot, plot_critical_window_diagnostic(smooth_critical_window(bw_critical_window_pls), bw_critical_window)),
  tar_target(ev_fixed_vs_discovered_plot, plot_fixed_vs_discovered(ev_critical_window)),
  tar_target(bw_fixed_vs_discovered_plot, plot_fixed_vs_discovered(bw_critical_window)),
  tar_target(al_loyo_comparison_plot, plot_loyo_comparison(al_loyo_summary))
)
