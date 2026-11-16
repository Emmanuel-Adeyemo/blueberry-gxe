# 09/14

library(dplyr)
library(purrr)
library(tibble)


# reproduce Sushan's interpolation (BM50.calc())

#' @description linearly interpolates the date at which a single plant's bloom percent first crosses a given threshold, if the threshold is never reached, returns the last observed date
#' if already met at the first obs, returns that date; 
#' @param dates vec of obs dates for one plant in one season
#' @param percents vec of bloom-percent values aligned with `dates`, for the same plant
#' @param threshold single bloom-percent threshold to interpolate the crossing date for e.g. 1, 5, 50, 100 ...
#' @output a single `Date` or `NA` if `percents` is all `NA`
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



#' @description gets per-plant bloom-timing traits for one environment's AL bloom data, with `.interpolate_threshold_date()` per plant per threshold
#' @param bloom_dta long-format bloom observation data for one envt (season x location), e.g. `al_bloom_long_clean` 
#' @param thresholds numeric vec of bloom-percent thresholds to interpolate a date for, defaults to `bloom1`/`bloom5`/`bloom10`/`bloom50`/`bloom100` at 1/5/10/50/100%
#' @output df
interpolate_bloom_doy = function(bloom_dta, thresholds = c(bloom1 = 1, bloom5 = 5, bloom10 = 10, bloom50 = 50, bloom100 = 100)) {

  plant_meta = bloom_dta %>%
    distinct(season, location, PlantID, BlockNo, SubplotNo, Genotype, Age)

  bloom_dta %>%
    arrange(season, location, PlantID, obs_date, Age) %>%
    group_by(season, location, PlantID) %>%
    reframe(
      map_dfc(thresholds, ~ .interpolate_threshold_date(obs_date, bloom_percent, .x))
    ) %>%
    left_join(plant_meta, by = c('season', 'location', 'PlantID')) %>%
    relocate('BlockNo', 'SubplotNo', 'Genotype', 'Age', .after = 'PlantID') %>%
    # added bloom duration since I wont get that from AUC - 09/18
    mutate(bloom_duration = bloom100 - bloom5)
  
}


# frost impact metric: AUC
#' @description computes the trapezoidal area-under-curve of a bloom-percent-over-time series for one plant, every plant's AUC is computed over the same fixed window. 
#' called internally by `compute_bloom_auc()`
#' @param dates vec of obs dates for one plant
#' @param values bloom-percent values aligned with `dates` 
#' @param window_start date to anchor the start of the AUC window to if earlier than the first obs
#' @param window_end Date to pad the end of the AUC window to if later than the last obs
#' @output a single numeric: the trapezoidal area under the bloom percent vs date curve between `window_start` and `window_end`
.trapezoidal_bloom_auc = function(dates, values, window_start, window_end) {
  ord = order(dates)
  dates = dates[ord]
  values = values[ord]

  if (dates[1] > window_start) {
    dates = c(window_start, dates)
    values = c(0, values)
  }
  last_i = length(dates)
  if (dates[last_i] < window_end) {
    dates = c(dates, window_end)
    values = c(values, values[last_i])
  }

  x = as.numeric(dates)
  n = length(x)
  sum(diff(x) * (values[-n] + values[-1]) / 2)
}


# compute AUC per plant for one envt on both raw and carry forward bloom percents
#' @description computes bloom-percent AUC for one envt, on both the raw and carry-forward-corrected bloom series, over a window from that season's anchor date to the 
#' last observed date and derives a "frost impact" metric as the corrected-minus-raw difference 
#' @param bloom_dta long format bloom observation data for one envt 
#' @param few_distinct_threshold flags plants with few obs, def = 3
#' @output df
compute_bloom_auc = function(bloom_dta, few_distinct_threshold = 3) {

  plant_meta = bloom_dta %>%
    distinct(season, location, PlantID, BlockNo, SubplotNo, Genotype, Age)


  season = unique(bloom_dta$season)
  if (length(season) != 1) {
    stop('compute_bloom_auc(): bloom_dta must be a single environment , got: ', paste(season, collapse = ', '))
  }

  window_start = season_anchor_date(season)
  window_end = max(bloom_dta$obs_date)
  window_length = as.numeric(window_end - window_start)

  bloom_dta %>%
    arrange(season, location, PlantID, obs_date, Age) %>%
    group_by(season, location, PlantID) %>%
    summarise(
      auc_raw_absolute = .trapezoidal_bloom_auc(obs_date, bloom_percent_raw, window_start, window_end),
      auc_corrected_absolute = .trapezoidal_bloom_auc(obs_date, bloom_percent, window_start, window_end),
      n_distinct_bloom_values = n_distinct(bloom_percent),
      .groups = 'drop'
    ) %>%
    mutate(
      auc_raw_percent = auc_raw_absolute / (window_length),
      auc_corrected_percent = auc_corrected_absolute / (window_length),
      frost_impact_absolute = auc_corrected_absolute - auc_raw_absolute,
      frost_impact_percent = auc_corrected_percent - auc_raw_percent,
      few_distinct_bloom_values = n_distinct_bloom_values <= few_distinct_threshold
    )%>%
    left_join(plant_meta, by = c('season', 'location', 'PlantID')) %>%
    relocate('BlockNo', 'SubplotNo', 'Genotype', 'Age', .after = 'PlantID')
}