# 10/03

library(dplyr)


build_per_trial_qc_dta_nc_fl = function(yield_dta, bloom_dta, block_col) {
  dta = list()
  results = list()

  # yield
  yield_trial_dta = yield_dta %>%
    transmute(genotype = genotype, block = .data[[block_col]], y = y)

  dta$TotalYield = yield_trial_dta
  results$TotalYield = fit_single_trial_qc(yield_trial_dta, 'TotalYield')

  # only if bloom dta exist - found some missingness in NC
  if (!is.null(bloom_dta) && nrow(bloom_dta) > 0) {
    for (bloom_col in c('bloom5', 'bloom50', 'bloom100')) {
      trial_dta = bloom_dta %>%
        transmute(genotype = genotype, block = .data[[block_col]], y = .data[[bloom_col]])

      dta[[bloom_col]] = trial_dta
      results[[bloom_col]] = fit_single_trial_qc(trial_dta, bloom_col)
    }

    duration_trial_dta = bloom_dta %>%
      transmute(genotype = genotype, block = .data[[block_col]], y = bloom_duration)

    dta$bloom_duration = duration_trial_dta
    results$bloom_duration = fit_single_trial_qc(duration_trial_dta, 'bloom_duration')
  }
  
  # like with al, dta is just for my own debugging
  list(dta = dta, results = bind_rows(results))
}