# 09/18



#' @description compues single-trial QC metrics from a fitted single trial model's variance
#' @param vc_g genetic variance component from a single-trial fit, e.g. `asreml_res$vc_g` from `fit_single_trial_asreml()`
#' @param vc_e residual variance component from the same fit, e.g. `asreml_res$vc_e`
#' @param mean_y trial mean of the response on the scale the model was fit on, e.g. `asreml_res$mean_y`, only used when `scale = 'raw'`
#' @param pev vec of genotype blup prediction error variances from the fit, e.g. `asreml_res$pev`
#' @param scale scale the response was fit on, controlling how CV is calc
#' @param raw_mean_y trial mean of the response on the original scale, used when `scale = 'sqrt'`
#' @output list with `H2`, `Re`, and `CV`.
compute_qc_metrics = function(vc_g, vc_e, mean_y, pev, scale = 'raw', raw_mean_y = NULL) {
  H2 = vc_g / (vc_g + vc_e)
  Re = max(0, 1 - mean(pev, na.rm = TRUE) / vc_g)

  CV = switch(
    scale,
    'raw' = (sqrt(vc_e) / mean_y) * 100,
    'log' = sqrt(vc_e) * 100,
    'sqrt' = {
      if (is.null(raw_mean_y)) stop("compute_qc_metrics(scale = 'sqrt') needs raw_mean_y")
      2 * sqrt(vc_e / raw_mean_y) * 100
    },
    stop("compute_qc_metrics: unknown scale '", scale, "'")
  )

  list(H2 = H2, Re = Re, CV = CV)
}


#' @description applies fixed pass/fail QC thresholds to a single-trial's H2/Re/CV metrics and reports whether all three pass
#' @param H2 broad-sense H2 for trial/trait, e.g. `decision_metrics$H2` 
#' @param Re reliability metric for trial/trait, e.g. `decision_metrics$Re` 
#' @param CV coeff of variation for trial/trait, e.g. `decision_metrics$CV`
#' @param cv_threshold max acceptable CV (%), defaults to 70.5.
#' @output list 
apply_thresholds = function(H2, Re, CV, cv_threshold = 70.5) {
  passes_H2 = H2 >= 0.1
  passes_Re = Re >= 0.2
  passes_CV = CV <= cv_threshold

  list(
    passes_H2 = passes_H2,
    passes_Re = passes_Re,
    passes_CV = passes_CV,
    passes_all = passes_H2 && passes_Re && passes_CV
  )
}