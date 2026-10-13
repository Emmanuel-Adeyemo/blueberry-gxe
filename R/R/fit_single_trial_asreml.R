library(dplyr)
library(asreml)


fit_single_trial_asreml = function(trial_dta) {

  if ('block' %in% names(trial_dta)){
    dta = trial_dta %>%
        mutate(genotype = factor(genotype), block = factor(block))


      dta = as.data.frame(dta)

      fit = tryCatch(
        asreml(
          fixed = y ~ block,
          random = ~genotype,
          residual = ~units,
          na.action = na.method(x = 'include'),
          data = dta
        ),
        error = function(e) NULL
      )
  }else{
    fit = NULL
  }
  

  if (is.null(fit)) {
    return(list(
      vc_g = NA_real_, vc_e = NA_real_,
      mean_y = mean(trial_dta$y, na.rm = TRUE), pev = NA_real_,
      converged = FALSE, singular = NA
    ))
  }
  fit

  vc = summary(fit)$varcomp
  vc_g = vc['genotype', 'component']
  vc_e = vc['units!R', 'component']  

  pred = predict(fit, classify = 'genotype', sed = TRUE)
  preds = pred$pvals
  # vdBLUP_mat = pred$sed^2
  # vdBLUP_avg = mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])
  pev = preds$std.error^2

  # h2_cullis = 1 - (vdBLUP_avg/(vc_g*2))

  list(
    vc_g = vc_g,
    vc_e = vc_e,
    mean_y = mean(trial_dta$y, na.rm = TRUE),
    pev = pev,
    # h2_cullis = h2_cullis,
    converged = TRUE,
    singular = isTRUE(vc_g < 1e-8)  
  )
}
