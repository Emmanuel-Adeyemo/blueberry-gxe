# 10/02

library(dplyr)

process_nc_fl_bloom = function(bloom_dta, max_days_since_anchor = 365) {
  .to_doy_flagged = function(date_col, season) {
    doy = days_since_season_anchor(date_col, season)
    out_of_range = !is.na(doy) & (doy < 0 | doy > max_days_since_anchor)
    if (any(out_of_range)) {
      warning(
        'process_nc_fl_bloom(): ', sum(out_of_range), ' bloom dates are more ',
        'than ', max_days_since_anchor, ' set to NA'
      )
    }
    if_else(out_of_range, NA_real_, doy)
  }

  bloom_dta %>%
    mutate(
      bloom5 = .to_doy_flagged(bloom5, season),
      bloom50 = .to_doy_flagged(bloom50, season),
      bloom100 = .to_doy_flagged(bloom100, season),
      bloom_duration = bloom100 - bloom5
    )
}