# 09/23

library(dplyr)


# from Sushan/Sakshi's script
# Chill: Oct1 - Feb 28
# Heat: Jan 1 - Feb 28
FIXED_CHILL_WINDOW = c(start = 30, end = 180)
FIXED_HEAT_WINDOW = c(start = 122, end = 180)


mmdd_to_days_since_anchor = function(mmdd, season) {
  month = mmdd %/% 100
  day = mmdd %% 100
  year = ifelse(month >= 9, season - 1L, season)
  real_date = as.Date(paste(year, month, day, sep = "-"))
  as.numeric(real_date - season_anchor_date(season))
}


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



compute_discovered_window_covariate = function(daily_chill_table, season, critical_window) {
  .sum_metric_over_stretches = function(stretches, metric_col) {
    if (nrow(stretches) == 0) {
      return(NA_real_)
    }
    values = vapply(seq_len(nrow(stretches)), function(i) {
      start_offset = mmdd_to_days_since_anchor(stretches$start_Date[i], season)
      end_offset = mmdd_to_days_since_anchor(stretches$end_Date[i], season)
      cov = compute_window_covariate(daily_chill_table, season, start_offset, end_offset)
      cov[[metric_col]]
    }, numeric(1))
    sum(values)
  }

  chill_stretches = critical_window %>% filter(.data$Type == "Chill")
  heat_stretches = critical_window %>% filter(.data$Type == "Heat")

  tibble(
    season = season,
    Chill_Portions_sum = .sum_metric_over_stretches(chill_stretches, "Chill_Portions_sum"),
    GDH_sum = .sum_metric_over_stretches(heat_stretches, "GDH_sum")
  )
}