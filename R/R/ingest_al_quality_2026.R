# 10/16

library(dplyr)
library(readxl)
library(tidyr)
library(stringr)


.al_2026_quality_id_cols = function(raw) {
  # same fix as ingest_al_dta.R's .al_resolve_id_cols() for the BW
  # SubloptNo/SubplotNo misspelling
  if ('SubloptNo' %in% names(raw)) {
    raw = raw %>% rename(SubplotNo = SubloptNo)
  }
  intersect(AL_ID_COLS, names(raw))
}


#' @param raw - single EV26/BW26 "Final" xlsx, Sheet1
#' @param season - 2026L
#' @param location - 'EV' or 'BW'
#' @output list(brix_long, ta_long) for this single environment
.process_al_2026_quality_envt = function(raw, season, location) {

  if ('SubloptNo' %in% names(raw)) {
    raw = raw %>% rename(SubplotNo = SubloptNo)
  }

  id_cols = .al_2026_quality_id_cols(raw)

  triples = list(
    c(date = 'Date1', brix = 'Brix1', ta = 'TA1'),
    c(date = 'Date2', brix = 'Brix2', ta = 'TA2'),
    c(date = 'Date3', brix = 'Brix3', ta = 'TA3')
  )

  long_one_triple = function(cols) {
    if (!all(cols %in% names(raw))) return(NULL)

    raw %>%
      dplyr::select(all_of(c(id_cols, unname(cols)))) %>%
      dplyr::transmute(
        across(all_of(id_cols)),
        season = season,
        location = location,
        obs_date = as.Date(.data[[cols[['date']]]]),
        brix = as.numeric(.data[[cols[['brix']]]]),
        ta = as.numeric(.data[[cols[['ta']]]])
      )
  }

  stacked = dplyr::bind_rows(lapply(triples, long_one_triple))

  if (is.null(stacked) || nrow(stacked) == 0) {
    return(list(brix_long = NULL, ta_long = NULL))
  }

  brix_long = stacked %>%
    dplyr::filter(!is.na(brix)) %>%
    dplyr::select(-ta) %>%
    dplyr::relocate(season, location, obs_date, .before = 1)

  ta_long = stacked %>%
    dplyr::filter(!is.na(ta)) %>%
    dplyr::select(-brix) %>%
    dplyr::relocate(season, location, obs_date, .before = 1)

  list(brix_long = brix_long, ta_long = ta_long)
}


#' @description ingest 2026 AL Brix/TA from the two "Final" xlsx files
#' @param ev26_path - path to EV26_071726_Final.xlsx
#' @param bw26_path - path to BW26_071726_Final.xlsx
#' @output list(brix_long, ta_long) combined across EV/BW 2026
ingest_al_quality_2026 = function(ev26_path, bw26_path) {

  ev26_raw = readxl::read_excel(ev26_path, sheet = 'Sheet1')
  bw26_raw = readxl::read_excel(bw26_path, sheet = 'Sheet1')

  ev26 = .process_al_2026_quality_envt(ev26_raw, season = 2026L, location = 'EV')
  bw26 = .process_al_2026_quality_envt(bw26_raw, season = 2026L, location = 'BW')

  list(
    brix_long = dplyr::bind_rows(ev26$brix_long, bw26$brix_long),
    ta_long = dplyr::bind_rows(ev26$ta_long, bw26$ta_long)
  )
}