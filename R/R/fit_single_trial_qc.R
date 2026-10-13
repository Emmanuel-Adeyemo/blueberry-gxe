# 09/18


library(dplyr)
library(tibble)


fit_single_trial_qc = function(trial_data, trait_name) {
  asreml_res = fit_single_trial_sommer(trial_data)
  
  sommer_res = fit_single_trial_sommer(trial_data)
  

  asreml_metrics = if (isTRUE(asreml_res$converged)) {
    compute_qc_metrics(asreml_res$vc_g, asreml_res$vc_e, asreml_res$mean_y, asreml_res$pev)
  } else {
    list(H2 = NA_real_, Re = NA_real_, CV = NA_real_)
  }

  sommer_metrics = if (isTRUE(sommer_res$converged)) {
    compute_qc_metrics(sommer_res$vc_g, sommer_res$vc_e, sommer_res$mean_y, sommer_res$pev)
  } else {
    list(H2 = NA_real_, Re = NA_real_, CV = NA_real_)
  }


  if (isTRUE(asreml_res$converged)) {
    decision_metrics = asreml_metrics
    decision_source = 'asreml'

  } else if (isTRUE(sommer_res$converged)) {
    decision_metrics = sommer_metrics
    decision_source = 'sommer_fallback'
  } else {
    decision_metrics = list(H2 = NA_real_, Re = NA_real_, CV = NA_real_)
    decision_source = 'neither_converged'
  }

  thresholds = if (!anyNA(unlist(decision_metrics))) {
    apply_thresholds(decision_metrics$H2, decision_metrics$Re, decision_metrics$CV)
  } else {
    list(passes_H2 = NA, passes_Re = NA, passes_CV = NA, passes_all = NA)
  }

  tibble(
    trait = trait_name,
    asreml_H2 = asreml_metrics$H2, asreml_Re = asreml_metrics$Re, asreml_CV = asreml_metrics$CV,
    sommer_H2 = sommer_metrics$H2, sommer_Re = sommer_metrics$Re, sommer_CV = sommer_metrics$CV,
    H2_disagreement = abs(sommer_metrics$H2 - asreml_metrics$H2),
    asreml_converged = isTRUE(asreml_res$converged),
    sommer_converged = isTRUE(sommer_res$converged),
    asreml_singular = asreml_res$singular,
    decision_source = decision_source,
    passes_H2 = thresholds$passes_H2,
    passes_Re = thresholds$passes_Re,
    passes_CV = thresholds$passes_CV,
    passes_all = thresholds$passes_all
  )
}
