# two fitting approaches
# nested and independent per plot before aggregating
# validate to choose best model


library(dplyr)
library(nlme)
library(purrr)
library(tibble)
library(nlme)

# season anchor is September 1 as per Sushan's suggestion
season_anchor_date = function(season) {
  as.Date(paste0(season - 1L, '-09-01'))
}

# days since season anchor
days_since_season_anchor = function(obs_date, season) {
  as.numeric(obs_date - season_anchor_date(season))
}

# SSlogis's is failing to auto generate starting values, this should manually take care of it
.estimate_logistic_start_values = function(x,y){
  asym_guess = max(y, na.rm=TRUE)
  half_target = asym_guess / 2
  xmid_guess = x[which.min(abs(y - half_target))]
  scal_guess = max(diff(range(x, na.rm=TRUE)) / 6, 1)
  c(Asym = asym_guess, xmid = xmid_guess, scal = scal_guess)
}


.add_subplot_id = function(bloom_dta){
  bloom_dta %>%
    mutate(subplot_id = paste(BlockNo, SubplotNo, sep = '_'))
}

# 4 parameter logistic curve
invert_logistic_day = function(Asym, xmid, scal, target){

  unreachable = target >= Asym
  if (any(unreachable, na.rm = T)) {
    warning(
      'Likely a poorly-fit curve.'
    )
  }
  out = xmid - scal * log(Asym / target - 1)
  out[unreachable] = NA_real_
  out
}


# fits independent - no borrowing - logistic curve
.fit_bloom_independent_one = function(plant_df){
  x = days_since_season_anchor(plant_df$obs_date, plant_df$season[1])
  y = plant_df$bloom_percent
  start_vals = .estimate_logistic_start_values(x, y)

  fit = tryCatch(
    nls(y ~ Asym / (1 + exp((xmid - x) / scal)), data = data.frame(x, y), start = start_vals),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return (list(Asym = NA_real_, xmid = NA_real_, scal = NA_real_, converged = FALSE))
  }
  cf = coef(fit)
  list(Asym = cf[['Asym']], xmid = cf[['xmid']], scal = cf[['scal']], converged = TRUE)
}


# per plot logistic curve
fit_bloom_independent = function(bloom_dta) {
  bloom_dta %>%
    group_by(PlantID) %>%
    group_modify(function(rows ) as_tibble(.fit_bloom_independent_one(rows))) %>%
    ungroup()
}


# complex model - nested
fit_bloom_nlme = function(bloom_dta, control = nlmeControl(maxIter = 100, msMaxIter = 100)){

  df = .add_subplot_id(bloom_dta) %>% 
    mutate(
      x = days_since_season_anchor(obs_date, season)
    )
}


print(days_since_season_anchor(as.Date('2022-12-25'),2023))