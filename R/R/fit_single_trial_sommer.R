# 09/18

library(dplyr)
library(sommer)


fit_single_trial_sommer = function(trial_dta) {

  if ('block' %in% names(trial_dta)){

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
}else {
    fit = NULL
  }
  if (is.null(fit)) {
    return(list(vc_g = NA_real_, vc_e = NA_real_, mean_y = mean(trial_dta$y, na.rm = TRUE), pev = NA_real_, converged = FALSE))
  }


  vc_g = as.numeric(fit$sigma$genotype)
  vc_e = as.numeric(fit$sigma$units)

  pev = diag(fit$PevU$genotype$y)

  list(
    vc_g = vc_g,
    vc_e = vc_e,
    mean_y = mean(trial_dta$y, na.rm = TRUE),
    pev = pev,
    converged = TRUE
  )
}
