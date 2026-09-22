# 09/18

compute_qc_metrics = function(sigma2_g, sigma2_e, mean_y, pev) {
  H2 = sigma2_g / (sigma2_g + sigma2_e)
  Ac = max(0, 1 - mean(pev, na.rm = TRUE) / sigma2_g)
  CV = (sqrt(sigma2_e) / mean_y) * 100

  list(H2 = H2, Ac = Ac, CV = CV)
}

# applies Bakare et al. 2022 threshold cutoffs: H² < 0.1, Ac < 0.4, CV > 40.5% -- 

apply_thresholds = function(H2, Ac, CV) {
  passes_H2 = H2 >= 0.1
  passes_Ac = Ac >= 0.4
  passes_CV = CV <= 40.5

  list(
    passes_H2 = passes_H2,
    passes_Ac = passes_Ac,
    passes_CV = passes_CV,
    passes_all = passes_H2 && passes_Ac && passes_CV
  )
}