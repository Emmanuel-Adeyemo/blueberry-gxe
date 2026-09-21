# 09/18


build_per_trial_qc_dta = function(yield_dta, bloom_dta){

  dta = list()

  results = list()

  # yield
  yield_subplot = aggregate_yield_for_qc(yield_dta)
  yield_trial_dta = yield_subplot %>%
    transmute(genotype = Genotype, block = BlockNo, y = TotalYield)

  dta$TotalYield = yield_trial_dta
  results$TotalYield = fit_single_trial_qc(yield_trial_dta, 'TotalYield')

  # bloom doy
  season = unique(bloom_dta$season)
  # print(length(season))
  doy_plant = interpolate_bloom_doy(bloom_dta)
  # print(colnames(doy_plant))
  doy_subplot = aggregate_bloom_to_subplot(doy_plant, c('bloom5', 'bloom50', 'bloom100', 'bloom_duration'))

  # print(unique(doy_subplot$season))
  # print(sapply(doy_subplot, class))
  # print(nrow(yield_dta))
  if (nrow(bloom_dta) > 0) {  for (doy_col in c('bloom5', 'bloom50', 'bloom100')) {
      trial_dta = doy_subplot %>%
        transmute(
          genotype = Genotype, block = BlockNo, y = days_since_season_anchor(.data[[doy_col]], .env$season)
        )
      
      dta[[doy_col]] = trial_dta
      results[[doy_col]] = fit_single_trial_qc(trial_dta, doy_col)
    }

    duration_trial_dta = doy_subplot %>%
      transmute(genotype = Genotype, block = BlockNo, y = as.numeric(bloom_duration))

    dta$bloom_duration = duration_trial_dta
    results$bloom_duration = fit_single_trial_qc(duration_trial_dta, 'bloom_duration')

    # bloom auc
    auc_plant = compute_bloom_auc(bloom_dta)
    # print(colnames(auc_plant))
    auc_subplot = aggregate_bloom_to_subplot(auc_plant, 'auc_corrected_percent')
    auc_trial_dta = auc_subplot %>%
      transmute(genotype = Genotype, block = BlockNo, y = auc_corrected_percent)

    dta$auc_corrected_percent = auc_trial_dta
    results$auc_corrected_percent = fit_single_trial_qc(auc_trial_dta, 'auc_corrected_percent')
}
  
  list(
    dta = dta,
    results = bind_rows(results)
  )
}