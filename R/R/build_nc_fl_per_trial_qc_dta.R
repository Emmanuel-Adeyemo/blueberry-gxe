# 10/03

library(dplyr)


build_per_trial_qc_dta_nc_fl = function(yield_dta, bloom_dta, block_col,
                                         avg_berry_weight_dta = NULL, firmness_dta = NULL,
                                         brix_dta = NULL, ta_dta = NULL) {
  dta = list()
  results = list()

  # yield
  # 10/22 - yield now fit on sqrt(y)
  yield_trial_dta = yield_dta %>%
    transmute(genotype = genotype, block = .data[[block_col]], y = sqrt(y))

  dta$TotalYield = yield_trial_dta
  results$TotalYield = fit_single_trial_qc(
    yield_trial_dta, 'TotalYield_sqrt', scale = 'sqrt',
    raw_mean_y = mean(yield_dta$y, na.rm = TRUE)
  )

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

  # 10/16 quality traits (avg berry weight, firmness, brix, and ta). no ph.
  quality_traits = list(
    avg_berry_weight = avg_berry_weight_dta,
    firmness = firmness_dta,
    brix = brix_dta,
    ta = ta_dta
  )

  for (trait_name in names(quality_traits)) {
    trait_dta = quality_traits[[trait_name]]
    if (!is.null(trait_dta) && nrow(trait_dta) > 0) {
      has_block = block_col %in% names(trait_dta) && any(!is.na(trait_dta[[block_col]]))

      trial_dta = if (has_block) {
        trait_dta %>% transmute(genotype = genotype, block = .data[[block_col]], y = .data[[trait_name]])
      } else {
        trait_dta %>% transmute(genotype = genotype, y = .data[[trait_name]])
      }

      dta[[trait_name]] = trial_dta
      results[[trait_name]] = fit_single_trial_qc(trial_dta, trait_name)
    }
  }

  # like with al, dta is just for my own debugging
  list(dta = dta, results = bind_rows(results))
}