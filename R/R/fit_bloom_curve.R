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

# SSlogis's is failing to auto generate starting values, this should manually take care of it
.estimate_logistic_start_values = function(x,y){
  asym_guess = max(y, na.rm=TRUE) * 1.02
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
.fit_bloom_independent_one = function(plant_df, min_visits = 4, timeout_secs = 30){

  failed <- list(Asym = NA_real_, xmid = NA_real_, scal = NA_real_, converged = FALSE, bound_pinned = NA)

  if (nrow(plant_df) < min_visits) {
    return(failed)
  }

  x = days_since_season_anchor(plant_df$obs_date, plant_df$season[1])
  y = plant_df$bloom_percent
  start_vals = as.list(.estimate_logistic_start_values(x, y))

  asym_upper = 105

  fit = tryCatch(
    {
      setTimeLimit(cpu = timeout_secs, elapsed = timeout_secs, transient = TRUE)
  
      nls(y ~ Asym / (1 + exp((xmid - x) / scal)), data = data.frame(x, y), start = start_vals, algorithm='port', 
        lower = c(Asym = 0, xmid = -Inf, scal = 0.5),
        upper = c(Asym = asym_upper, xmid = Inf, scal = 200))
  
    },
    error = function(e) NULL,
    finally = setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)
  )
    
  if (is.null(fit)) {
    return(failed)
  }
  cf = coef(fit)
  list(Asym = cf[['Asym']], xmid = cf[['xmid']], scal = cf[['scal']], converged = TRUE, 
  bound_pinned = isTRUE(all.equal(cf[['Asym']], asym_upper, tolerance = 0.01)))
}


# per plot logistic curve
fit_bloom_independent = function(bloom_dta) {
  bloom_dta %>%
    group_by(PlantID) %>%
    group_modify(function(rows, ... ) as_tibble(.fit_bloom_independent_one(rows))) %>%
    ungroup()
}


# complex model - nested
fit_bloom_nlme = function(bloom_dta, control = nlme::nlmeControl(maxIter = 300, msMaxIter = 300)){

  df = .add_subplot_id(bloom_dta) %>% 
    mutate(
      x = days_since_season_anchor(obs_date, season),
      # BlockNo  =factor(BlockNo), changing to 2 level nesting
      subplot_id = factor(subplot_id),
      PlantID = factor(PlantID)
    )
  
  indy = fit_bloom_independent(bloom_dta)
  clean_indy = indy %>% filter(converged, !bound_pinned)

  if (nrow(clean_indy) >= 5) {
    start_vals = c(
      Asym = median(clean_indy$Asym, na.rm = TRUE),
      xmid = median(clean_indy$xmid, na.rm = TRUE),
      log_scal = median(log(clean_indy$scal), na.rm = TRUE)
    )
  } else {
    pooled = .estimate_logistic_start_values(df$x, df$bloom_percent)
    start_vals = c(Asym = pooled[["Asym"]], xmid = pooled[["xmid"]], log_scal = log(pooled[["scal"]]))
  }
  # start_vals = .estimate_logistic_start_values(df$x, df$bloom_percent)
  fit_a = tryCatch(
    nlme(
      bloom_percent ~ Asym / (1 + exp((xmid - x) / exp(log_scal))),
      data = df,
      fixed = Asym + xmid + log_scal ~ 1,
      random = list(
        # BlockNo = nlme::pdDiag(Asym + xmid + scal ~1),
        subplot_id = nlme::pdDiag(Asym + xmid + log_scal ~1),
        PlantID = nlme::pdDiag(xmid ~1)
      ), 
      start = start_vals,
      control = control
    ),
    error = function(e) NULL
  )
  if(!is.null(fit_a)) {
    attr(fit_a, 'model_used') = 'A'
    return (fit_a)
  }

  # fall back model for ev25 and ev26
  fit_b = nlme(
    bloom_percent ~ Asym / (1 + exp((xmid - x) / exp(log_scal))),
    data = df,
    fixed = Asym + xmid + log_scal ~ 1,
    random = list(
      subplot_id = nlme::pdDiag(xmid+ log_scal ~ 1),
      PlantID = nlme::pdDiag((xmid + log_scal ~ 1))
    ),
    start = start_vals,
    control = control
  )
  attr(fit_b, 'model_used') = 'B'
  fit_b
  
}


# extract params from fitted nlme
extract_plant_coefs_nlme = function(fitted_model){
  plant_coefs = coef(fitted_model, level = 2) # changed levels to 2, after block effect carried no weight
  model_used = attr(fitted_model, 'model_used')
  if(is.null(model_used)) model_used = NA_character_

  tibble(
    PlantID = rownames(plant_coefs),
    Asym = plant_coefs[['Asym']],
    xmid = plant_coefs[['xmid']],
    scal = plant_coefs[['log_scal']],
    model_spec = model_used
  ) %>% 
    mutate(PlantID = str_extract_plant_id(PlantID))

}


# rownames at innermost grouping levels not just PlantID
str_extract_plant_id = function(x){
  vapply(strsplit(x, '/', fixed=TRUE), function(parts) parts[length(parts)], character(1))
}


# converts into needed output
doy_from_coefs = function(plant_coefs, low =5, mid = 50, full_frac = 0.95) {
  plant_coefs %>% 
    mutate(
      doy_5 = invert_logistic_day(Asym, xmid, scal, low),
      doy_50 = invert_logistic_day(Asym, xmid, scal, mid),
      doy_full = invert_logistic_day(Asym, xmid, scal, full_frac * Asym)
    ) %>% select(
      PlantID, doy_5, doy_50, doy_full, any_of('model_spec')
    )
}


# convert back to date
season_doy_to_date = function(season, day_since_anchor) {
  season_anchor_date(season) + round(day_since_anchor)
}


# flag date estimates that fall outside plant's actual visit window especially nls()
flag_extrapolated_day_estimates = function(doy_table, bloom_dta, buffer_days = 10) {
  plant_ranges = bloom_dta %>%
    mutate(x = days_since_season_anchor(obs_date, season)) %>% 
    group_by(PlantID) %>%
    summarise(min_x = min(x), max_x = max(x), .groups = 'drop')

  doy_table %>%
    left_join(plant_ranges, by = 'PlantID') %>%
    mutate(
      extrapolated = if_else(
        is.na(day_50), NA, doy_5 < min_x - buffer_Days | doy_full > max_x + buffer_days
      )
    ) %>%
    select(-'min_x', -'max_x')
}


# entry for nlme
harmonize_bloom_doy_nlme = function(bloom_dta, low = 5, mid = 50, full_frac = 0.95) {
  season = unique(bloom_dta$season)
  if (length(season) != 1){
    stop('bloom_dta must be a single envt ie one season, got: ', paste(season, collaspe = ", "))
  }

  fitted_model = fit_bloom_nlme(bloom_dta)
  extract_plant_coefs_nlme(fitted_model) %>%
    doy_from_coefs(low = low, mid = mid, full_frac = full_frac) %>%
    mutate(
      date_5 = season_doy_to_date(season, doy_5),
      date_50 = season_doy_to_date(season, doy_50),
      date_full = season_doy_to_date(season, doy_full)
    )
}


# entry for independent
harmonize_bloom_doy_independent = function(bloom_dta, low = 5, mid = 50, full_frac = 0.95) {
  season = unique(bloom_dta$season)
  if (length(season) != 1) {
    stop('bloom_dta must be a single envt ie one season, got: ', paste(season, collapse = ', '))
  }

  fit_bloom_independent(bloom_dta) %>%
    doy_from_coefs(low = low, mid = mid, full_frac = full_frac) %>%
    mutate(
      date_5 = season_doy_to_date(season, doy_5),
      date_50 = season_doy_to_date(season, doy_50),
      date_full = season_doy_to_date(season, doy_full)
    )
}


  # print(days_since_season_anchor(as.Date('2022-12-25'),2023))



