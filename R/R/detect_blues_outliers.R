# 09/13


library(dplyr)
library(purrr)
library(asreml)


build_trait_blues_dta = function(trial_dta_list, fit_fun = fit_single_trial_blue_asreml) {
  purrr::imap_dfr(trial_dta_list, function(trial_dta, env_name) {
    if (is.null(trial_dta) || nrow(trial_dta) == 0) return(NULL)
    blue = fit_fun(trial_dta)
    if (is.null(blue) || nrow(blue) == 0) return(NULL)
    blue %>% mutate(Env = env_name, .before = 1)
  })
}


detect_blue_outliers_asreml = function(blues_dta, trt) {

  dta = blues_dta %>%
    transmute(
      Genotype = factor(Genotype),
      Env = factor(Env),
      predicted_value = predicted_value,
      std_error = std_error,
      wts = 1 / std_error^2
    ) %>%
    as.data.frame()


  fit = asreml::asreml(
    fixed = predicted_value ~ Env + Genotype,
    random = ~Genotype:Env,
    data = dta,
    residual = ~dsum(~id(units) | Env),
    predict = asreml::predict.asreml(classify = 'Genotype', ignore = 'Genotype:Env'),
    weights = 'wts',
    aom = TRUE,
    sparse = ~Genotype,
    family = asreml::asr_gaussian(dispersion = 1)
  )

  Rinv = fit$aom$R[, 2]

  # alpha for a ~3-SD (99.73%) two-sided cutoff
  alpha = 1 - .9973
  critical_value = qt(alpha / 2, df = fit$nedf, lower.tail = FALSE)

  is_outlier = abs(Rinv) > critical_value
  is_outlier[is.na(is_outlier)] = FALSE

  outliers = droplevels(dta[is_outlier, ])
  colnames(outliers) = gsub('predicted_value', trt, colnames(outliers))
  colnames(outliers) = gsub('std_error', paste0(trt, '_SE'), colnames(outliers))
  colnames(outliers) = gsub('wts', paste0(trt, '_wt'), colnames(outliers))
  outliers = cbind(join = paste0(outliers$Genotype, outliers$Env), outliers)

  blues_clean = droplevels(dta[!is_outlier, ])

  list(
    fit = fit,
    critical_value = critical_value,
    outliers = outliers,
    blues_clean = blues_clean
  )
}