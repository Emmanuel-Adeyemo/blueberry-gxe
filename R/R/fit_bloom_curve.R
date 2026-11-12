# 09/12
# two fitting approaches
# nested and independent per plot before aggregating
# validate to choose best model


library(dplyr)
library(nlme)
library(purrr)
library(tibble)

# season anchor is September 1 as per Sushan's suggestion
season_anchor_date = function(season) {
  as.Date(paste0(season - 1L, '-09-01'))
}

# days since season anchor
days_since_season_anchor = function(obs_date, season) {
  as.numeric(obs_date - season_anchor_date(season))
}

