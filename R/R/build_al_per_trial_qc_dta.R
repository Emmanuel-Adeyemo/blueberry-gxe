# 09/18


#' @description builds single-trial (one Env = one season x location) QC results for AL trials. Fits traits via `fit_single_trial_qc()`
#' @param yield_dta raw long-format yield data for one Env e.g. `al_yield_long_clean` (output of `clean_al_trait_long()` from stage03_cleaning.R) filtered to a single season/loc
#' @param bloom_dta raw long-format bloom observation data for the same env as `yield_dta` e.g. `al_bloom_long_clean` 
#' @param plant_height_dta optional - raw long-format plant height data for the same Env e.g. `al_plant_height_long_clean` which is the output from `clean_al_trait_long()` filtered to seas/loc
#' @param avg_berry_weight_dta optional - raw long-format of average berry weight data for the same Env e.g. `al_avg_berry_weight_long_clean`, also output from `clean_al_trait_long()`
#' @param firmness_dta optional - raw long format firmness data for same env e.g. `al_firmness_long_clean` output of `clean_al_trait_long()`
#' @param brix_dta optional - raw long format brix/tss data for same env e.g. `al_brix_long_clean` output of `clean_al_trait_long()`
#' @param ta_dta optional - raw long format ta data e.g. `al_ta_long_clean`, output of `clean_al_trait_long()`
#' @output a list with `dta` ie data actually used to fit model (not targets output) and `results` - one df with every trait's qc results
#' 
build_per_trial_qc_dta = function(yield_dta, bloom_dta,
                                   plant_height_dta = NULL, avg_berry_weight_dta = NULL,
                                   firmness_dta = NULL, brix_dta = NULL, ta_dta = NULL){

  dta = list()

  results = list()

  # yield
  # 10/22 - yield now fit on sqrt(TotalYield)
  yield_subplot = aggregate_yield_for_qc(yield_dta)
  yield_trial_dta = yield_subplot %>%
    transmute(genotype = Genotype, block = BlockNo, y = sqrt(TotalYield))

  dta$TotalYield = yield_trial_dta
  results$TotalYield = fit_single_trial_qc(
    yield_trial_dta, 'TotalYield_sqrt', scale = 'sqrt',
    raw_mean_y = mean(yield_subplot$TotalYield, na.rm = TRUE)
  )

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

  # 10/16 quality traits (plant height, avg berry weight, firmness, brix, and ta)
  quality_traits = list(
    plant_height = plant_height_dta,
    avg_berry_weight = avg_berry_weight_dta,
    firmness = firmness_dta,
    brix = brix_dta,
    ta = ta_dta
  )

  for (trait_name in names(quality_traits)) {
    trait_dta = quality_traits[[trait_name]]
    if (!is.null(trait_dta) && nrow(trait_dta) > 0) {
      trait_subplot = aggregate_quality_trait_for_qc(trait_dta, trait_name)
      trial_dta = trait_subplot %>%
        transmute(genotype = Genotype, block = BlockNo, y = .data[[trait_name]])

      dta[[trait_name]] = trial_dta
      results[[trait_name]] = fit_single_trial_qc(trial_dta, trait_name)
    }
  }

  # dta is for my own debugging, to explore the dta that goes into the model, not a tar_target output

  list(
    dta = dta,
    results = bind_rows(results)
  )
}