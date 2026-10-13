# 09/12

library(dplyr)
library(purrr)
library(tibble)

# reproduce Sushan's interpolation (BM50.calc())
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



calc_b50_interpolation = function(bloom_data, thresholds = c(bloom1 = 1, bloom5 = 5, bloom10 = 10, bloom50 = 50, bloom100 = 100)) {
  bloom_data %>%
    arrange(PlantID, obs_date) %>%
    group_by(PlantID) %>%
    reframe(
      map_dfc(thresholds, ~ .interpolate_threshold_date(obs_date, bloom_percent, .x))
    )
}

# calc_b50_interpolation <- function(bloom_dta, thresholds = c(bloom1 = 1, bloom5 = 5, bloom10 = 10, bloom50 = 50, bloom100 = 100)) {
#   bloom_dta %>%
#     arrange(PlantID, obs_date) %>%
#     group_by(PlantID) %>%
#     group_modify(function(rows, ...) {
#       as_tibble(map(thresholds, function(t) {
#         .interpolate_threshold_date(obs_date, bloom_percent, t)
#       }))
#     }) %>%
#     ungroup()
# }





# cross validation of the models - leave visit out.
leave_visits_out_cv = function(bloom_dta, min_visits = 6, drop_frac = 0.4, n_plants = 20, n_reps = 5, seed = 785) {
  set.seed(seed)

  visit_counts = bloom_dta %>% count(PlantID, name = 'n_visits')
  eligible = visit_counts %>% filter(n_visits >= min_visits)

  if(nrow(eligible) == 0) {
    warning('leave_visits_out_cv(): no plants meet min_visits = ', min_visits, 'output is empty.')
    return(tibble(
      PlantID = character(0), rep = integer(0), approach = character(0),
      doy_50_full = double(0), doy_50_thinned = double(0), error = double(0),
      converged = logical(0)
    ))
  }

  test_plants = eligible$PlantID
  if(!is.null(n_plants) && length(test_plants) > n_plants) {
    test_plants = sample(test_plants, n_plants)

  }

  full_independent = fit_bloom_independent(bloom_dta) %>%
    doy_from_coefs() %>%
    dplyr::select(PlantID, doy_50_full = 'doy_50')

  full_nlme_fit = fit_bloom_nlme(bloom_dta) 
  full_nlme = extract_plant_coefs_nlme(full_nlme_fit) %>%
    doy_from_coefs() %>%
    dplyr::select(PlantID, doy_50_full = 'doy_50')

  results = map_dfr(test_plants, function(pid) {
    plant_visits = bloom_dta %>% filter(PlantID == pid)
    n_keep = max(2, round(nrow(plant_visits) * (1 - drop_frac)))

    map_dfr(seq_len(n_reps), function(r) {
      kept_rows = plant_visits %>% slice_sample(n = n_keep) %>% arrange(obs_date)
      thinned_dta = bloom_dta %>%
        filter(PlantID != pid) %>%
        bind_rows(kept_rows)

      # indy fit
      indy_thin = .fit_bloom_independent_one(kept_rows)
      indy_doy50 = if(indy_thin$converged) {
        invert_logistic_day(indy_thin$Asym, indy_thin$xmid, indy_thin$scal, 50)
      } else{
        NA_real_
      }

      # nlme fit
      nlme_thin_fit = tryCatch(fit_bloom_nlme(thinned_dta), error = function(e) NULL)
      nlme_converged = !is.null(nlme_thin_fit)
      nlme_doy50 = if (nlme_converged) {
        coefs = extract_plant_coefs_nlme(nlme_thin_fit) %>% filter(PlantID == pid)
        if(nrow(coefs) == 1) invert_logistic_day(coefs$Asym, coefs$xmid, coefs$scal, 50) else NA_real_
      }else{
        NA_real_
      }

      bind_rows(
        tibble(
          PlantID = pid, rep = r, approach = 'independent',
          doy_50_full = full_independent$doy_50_full[full_independent$PlantID == pid],
          doy_50_thinned = indy_doy50, converged = indy_thin$converged
        ),
        tibble(
          PlantID = pid, rep = r, approach = 'nlme',
          doy_50_full = full_nlme$doy_50_full[full_nlme$PlantID == pid],
          day_50_thinned = nlme_doy50, converged = nlme_converged
        )
      )
    })
  })
  results %>% mutate(error = doy_50_thinned - doy_50_full)

}




# dta agreement
agreement_by_dta_density = function(bloom_dta) {

  visit_counts = bloom_dta %>% count(PlantID, name='n_visits')

  indy = fit_bloom_independent(bloom_dta) %>%
    doy_from_coefs() %>%
    dplyr::select(PlantID, day_50_independent = 'doy_50')

  nlme_fit = fit_bloom_nlme(bloom_dta)
  mm = extract_plant_coefs_nlme(nlme_fit) %>%
    doy_from_coefs() %>%
    dplyr::select(PlantID, doy_50_nlme = 'doy_50')

  visit_counts %>%
    left_join(indy, by = 'PlantID') %>%
    left_join(mm, by = 'PlantID') %>%
    mutate(abs_diff = abs(doy_50_nlme - doy_50_independent))
}




extract_variance_components = function(fitted_model) {
  vc = VarCorr(fitted_model)
  rn = rownames(vc)
  is_header = grepl('=', rn, fixed = TRUE)

  level = character(length(rn))
  current = NA_character_
  for (i in seq_along(rn)) {
    if(is_header[i]) {
      current = sub(' =$', '', rn[i])
    }
    level[i] = current
  }

  tibble (

    level = ifelse(rn == 'Residual', 'Residual', level),
    parameter = rn,
    variance = suppressWarnings(as.numeric(vc[, 'Variance']))
  ) %>%
    filter(!is_header)
}