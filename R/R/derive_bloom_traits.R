# 09/14

library(dplyr)
library(purrr)
library(tibble)


# reproduce Sushan's interpolation (BM50.calc())
.interpolate_threshold_date = function(dates, percents, threshold){
  if(all(is.na(percents))) {
    return(as.Date(NA))

  }
  r2 = which(percents >= threshold)[1]

  if(is.na(r2)) {
    return (max(dates, na.rm = TRUE))
  }
  if (r2 == 1){
    return (dates[r2])
  }

  r1 = r2 - 1
  d1 = dates[r1]
  d2 = dates[r2]
  b1 = percents[r1]
  b2 = percents[r2]
  d1 + (threshold - b1)* as.numeric(d2 - d1) / (b2 - b1)
}



interpolate_bloom_doy = function(bloom_dta, thresholds = c(bloom1 = 1, bloom5 = 5, bloom10 = 10, bloom50 = 50, bloom100 = 100)) {

  plant_meta = bloom_dta %>%
    distinct(season, location, PlantID, BlockNo, SubplotNo, Genotype)

  bloom_dta %>%
    arrange(season, location, PlantID, obs_date) %>%
    group_by(season, location, PlantID) %>%
    reframe(
      map_dfc(thresholds, ~ .interpolate_threshold_date(obs_date, bloom_percent, .x))
    ) %>%
    left_join(plant_meta, by = c('season', 'location', 'PlantID')) %>%
    relocate('BlockNo', 'SubplotNo', 'Genotype', .after = 'PlantID')
  
}


# frost impact metric: AUC
.trapezoidal_bloom_auc = function(dates, values, window_start, window_end) {
  ord = order(dates)
  dates = dates[ord]
  values = values[ord]

  if (dates[1] > window_start) {
    dates = c(window_start, dates)
    values = c(0, values)
  }
  last_i <- length(dates)
  if (dates[last_i] < window_end) {
    dates = c(dates, window_end)
    values = c(values, values[last_i])
  }

  x = as.numeric(dates)
  n = length(x)
  sum(diff(x) * (values[-n] + values[-1]) / 2)
}


# compute AUC per plant for one envt on both raw and carry forward bloom percents
compute_bloom_auc = function(bloom_data, few_distinct_threshold = 3) {
  season = unique(bloom_data$season)
  if (length(season) != 1) {
    stop("compute_bloom_auc(): bloom_data must be a single environment (one season), got: ", paste(season, collapse = ', '))
  }

  window_start = season_anchor_date(season)
  window_end = max(bloom_data$obs_date)
  window_length = as.numeric(window_end - window_start)

  bloom_data %>%
    arrange(PlantID, obs_date) %>%
    group_by(PlantID) %>%
    summarise(
      auc_raw_absolute = .trapezoidal_bloom_auc(obs_date, bloom_percent_raw, window_start, window_end),
      auc_corrected_absolute = .trapezoidal_bloom_auc(obs_date, bloom_percent, window_start, window_end),
      n_distinct_bloom_values = n_distinct(bloom_percent),
      .groups = "drop"
    ) %>%
    mutate(
      auc_raw_relative = auc_raw_absolute / (100 * window_length),
      auc_corrected_relative = auc_corrected_absolute / (100 * window_length),
      frost_impact_absolute = auc_corrected_absolute - auc_raw_absolute,
      frost_impact_relative = auc_corrected_relative - auc_raw_relative,
      few_distinct_bloom_values = n_distinct_bloom_values <= few_distinct_threshold
    )
}