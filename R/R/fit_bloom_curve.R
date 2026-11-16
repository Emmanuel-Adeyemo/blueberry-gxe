# 09/12
# two fitting approaches
# nested and independent per plot before aggregating
# validate to choose best model


library(dplyr)
library(nlme)
library(purrr)
library(tibble)

# season anchor is September 1 as per Sushan's suggestion
#' @description computes the fixed calendar "anchor" date for a given bloom/fruiting season, Sept 1 of preceeding year
#' @param season the season (year) to compute the anchor for, a single value or vec
#' @output `Date` (vectorized over `season`): `<season - 1>-09-01`.
season_anchor_date = function(season) {
  as.Date(paste0(season - 1L, '-09-01'))
}

# days since season anchor
#' @description converts a calendar date to "days since season anchor" 
#' @param obs_date vec of obs dates
#' @param season season/year corresponding to each `obs_date` 
#' @output numeric vec: days elapsed from that season's anchor date to `obs_date` 
days_since_season_anchor = function(obs_date, season) {
  as.numeric(obs_date - season_anchor_date(season))
}

