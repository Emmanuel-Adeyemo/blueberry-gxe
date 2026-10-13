09/19

library(lme4)
library(dplyr)

fit_single_trial_lme4 = function(trial_dta) {
  dta = trial_dta %>%
    mutate(genotype = factor(genotype), block = factor(block))

  fit = tryCatch(
    lmer(y ~ block + (1 | genotype), data = dta),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(vc_g = NA_real_, vc_e = NA_real_, mean_y = mean(dta$y, na.rm = TRUE), pev = NA_real_, converged = FALSE, singular = NA))
  }

  vc = as.data.frame(VarCorr(fit))
  vc_g = vc$vcov[vc$grp == 'genotype']
  vc_e = vc$vcov[vc$grp == 'Residual']

  re = ranef(fit, condVar = TRUE)
  pev = attr(re$genotype, 'postVar')[1, 1, ]

  list(vc_g = vc_g,
  vc_e = vc_e,
  mean_y = mean(dta$y, na.rm = T),
  pev = pev,
  converged = TRUE,
  singular = isSingular(fit)
  )
}