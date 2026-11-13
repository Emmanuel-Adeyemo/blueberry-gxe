# 10/03

library(dplyr)

#' @description NC/FL version of `build_per_trial_qc_dta()` builds single-trial qc results for NC and FL trials. Already in the `genotype`/`y` style. Plnat height dta not avail
#' @param yield_dta raw long-format yield data for one Env e.g. `nc_yield_clean` (output of `clean_nc_genotypes()` from stage03_cleaning.R) filtered to a single season/loc
#' @param bloom_dta raw long-format bloom observation data for the same env as `yield_dta` e.g. `nc_bloom_doy` utput of `process_nc_fl_bloom()`
#' @param block_col name of the locking column in the site e.g. `subplot` for NC/FL dta as opposed to `BlockNo`in AL data
#' @param avg_berry_weight_dta optional - raw long-format of average berry weight data for the same Env e.g. `nc_avg_berry_weight_clean`, also output from `clean_nc_genotypes()`
#' @param firmness_dta optional - raw long format firmness data for same env e.g. `nc_firmness_clean` output of `clean_nc_genotypes()`
#' @param brix_dta optional - raw long format brix/tss data for same env e.g. `nc_brix_clean` output of `clean_nc_genotypes()`
#' @param ta_dta optional - raw long format ta data e.g. `nc_ta_clean`, output of `clean_nc_genotypes()`
#' @output a list with `dta` ie data actually used to fit model (not targets output) and `results` - one df with every trait's qc results
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