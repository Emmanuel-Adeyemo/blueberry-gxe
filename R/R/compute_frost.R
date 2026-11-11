# 10/29

library(dplyr)

compute_frost_score = function(weather_prepared, season, frost_threshold = -2,
                                scoring_window = c(start = 195, end = 211), frost_weight = 1) {
  
  anchor = season_anchor_date(season)
  obs_date = as.Date(paste(weather_prepared$Year, weather_prepared$Month, weather_prepared$Day, sep = "-"))
  days_since_anchor = as.numeric(obs_date - anchor)

  
  dta = weather_prepared %>%
    mutate(days_since_anchor = days_since_anchor) %>%
    filter(days_since_anchor >= 0, days_since_anchor <= 211)

  late_window = dta %>%
    filter(days_since_anchor >= scoring_window["start"], days_since_anchor <= scoring_window["end"])

  frost_occurred = any(late_window$Temp <= frost_threshold, na.rm = TRUE)

  tibble(
    season = season,
    n_frost_hours_available = nrow(dta),
    frost_late_window = frost_occurred,
    frost_late_window_min_temp = suppressWarnings(min(late_window$Temp, na.rm = TRUE)),
    frost_score = if (frost_occurred) frost_weight else 0
  )
}