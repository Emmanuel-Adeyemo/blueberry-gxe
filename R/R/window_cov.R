# 09/23

library(dplyr)


# from Sushan/Sakshi's script
# Chill: Oct1 - Feb 28
# Heat: Jan 1 - Feb 28
# Both windows are expressed as day-offsets from season_anchor_date() 
FIXED_CHILL_WINDOW = c(start = 30, end = 180)
FIXED_HEAT_WINDOW = c(start = 122, end = 180)


compute_window_covariate = function(daily_chill_table, season, start_offset, end_offset) {
  anchor = season_anchor_date(season)
  obs_date = as.Date(paste(daily_chill_table$Year, daily_chill_table$Month, daily_chill_table$Day, sep = "-"))
  days_since_anchor = as.numeric(obs_date - anchor)

  daily_chill_table %>%
    mutate(days_since_anchor = days_since_anchor) %>%
    filter(.data$days_since_anchor >= start_offset, .data$days_since_anchor <= end_offset) %>%
    summarise(
      season = season,
      Chill_Portions_sum = sum(.data$Chill_Portions, na.rm = TRUE),
      GDH_sum = sum(.data$GDH, na.rm = TRUE)
    )
}


compute_fixed_window_covariate = function(daily_chill_table, season) {
  chill_part = compute_window_covariate(daily_chill_table, season, FIXED_CHILL_WINDOW["start"], FIXED_CHILL_WINDOW["end"])
  heat_part = compute_window_covariate(daily_chill_table, season, FIXED_HEAT_WINDOW["start"], FIXED_HEAT_WINDOW["end"])
  tibble(season = season, Chill_Portions_sum = chill_part$Chill_Portions_sum, GDH_sum = heat_part$GDH_sum)
}


# 10/29
get_available_seasons = function(weather_prepared, location, candidate_seasons) {
  obs_date = as.Date(paste(weather_prepared$Year, weather_prepared$Month, weather_prepared$Day, sep = "-"))
  data_min = min(obs_date, na.rm = TRUE)
  data_max = max(obs_date, na.rm = TRUE)

  eligible = Filter(function(season) {
    anchor = season_anchor_date(season)
    data_min <= anchor && data_max >= (anchor + 211)
  }, candidate_seasons)

  tibble(location = location, season = eligible)
}