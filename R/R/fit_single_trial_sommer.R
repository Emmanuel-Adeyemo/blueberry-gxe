# 09/18

library(dplyr)
library(sommer)


fit_single_trial_sommer = function(trial_dta) {


  dta = trial_dta %>%
    mutate(genotype = factor(genotype), block = factor(block))

  fit = tryCatch(
    mmer(
      fixed = y ~ block,
      random = ~genotype,
      rcov = ~units,
      data = dta,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(sigma2_g = NA_real_, sigma2_e = NA_real_, mean_y = mean(df$y, na.rm = TRUE), pev = NA_real_, converged = FALSE))
  }


  sigma2_g = as.numeric(fit$sigma$genotype)
  sigma2_e = as.numeric(fit$sigma$units)

  pev = diag(fit$PevU$genotype$y)

  list(
    sigma2_g = sigma2_g,
    sigma2_e = sigma2_e,
    mean_y = mean(dta$y, na.rm = TRUE),
    pev = pev,
    converged = TRUE
  )
}
