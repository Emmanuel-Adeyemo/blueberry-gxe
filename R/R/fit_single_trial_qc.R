# 09/18


library(dplyr)
library(tibble)


fit_single_trial_qc = function(trial_data, trait_name) {
  lme4_res = fit_single_trial_lme4(trial_data)
  
  sommer_res = fit_single_trial_sommer(trial_data)
  

  lme4_metrics = if (isTRUE(lme4_res$converged)) {
    compute_qc_metrics(lme4_res$sigma2_g, lme4_res$sigma2_e, lme4_res$mean_y, lme4_res$pev)
  } else {
    list(H2 = NA_real_, Ac = NA_real_, CV = NA_real_)
  }

  sommer_metrics = if (isTRUE(sommer_res$converged)) {
    compute_qc_metrics(sommer_res$sigma2_g, sommer_res$sigma2_e, sommer_res$mean_y, sommer_res$pev)
  } else {
    list(H2 = NA_real_, Ac = NA_real_, CV = NA_real_)
  }


  if (isTRUE(sommer_res$converged)) {
    decision_metrics = sommer_metrics
    decision_source = "sommer"
  } else if (isTRUE(lme4_res$converged)) {
    decision_metrics = lme4_metrics
    decision_source = "lme4_fallback"
  } else {
    decision_metrics = list(H2 = NA_real_, Ac = NA_real_, CV = NA_real_)
    decision_source = "neither_converged"
  }

  thresholds = if (!anyNA(unlist(decision_metrics))) {
    apply_thresholds(decision_metrics$H2, decision_metrics$Ac, decision_metrics$CV)
  } else {
    list(passes_H2 = NA, passes_Ac = NA, passes_CV = NA, passes_all = NA)
  }

  tibble(
    trait = trait_name,
    sommer_H2 = sommer_metrics$H2, sommer_Ac = sommer_metrics$Ac, sommer_CV = sommer_metrics$CV,
    lme4_H2 = lme4_metrics$H2, lme4_Ac = lme4_metrics$Ac, lme4_CV = lme4_metrics$CV,
    H2_disagreement = abs(sommer_metrics$H2 - lme4_metrics$H2),
    sommer_converged = isTRUE(sommer_res$converged),
    lme4_converged = isTRUE(lme4_res$converged),
    lme4_singular = lme4_res$singular,
    decision_source = decision_source,
    passes_H2 = thresholds$passes_H2,
    passes_Ac = thresholds$passes_Ac,
    passes_CV = thresholds$passes_CV,
    passes_all = thresholds$passes_all
  )
}
