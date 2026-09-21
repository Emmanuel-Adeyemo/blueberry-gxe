# 09/18

compute_qc_metrics = function(sigma2_g, sigma2_e, mean_y, pev) {
  H2 = sigma2_g / (sigma2_g + sigma2_e)
  Re = max(0, 1 - mean(pev, na.rm = TRUE) / sigma2_g)
  CV = (sqrt(sigma2_e) / mean_y) * 100

  list(H2 = H2, Re = Re, CV = CV)
}

# applies Bakare et al. 2022 threshold cutoffs: H² < 0.1, Re < 0.4^2, CV > 40.5% -- 

apply_thresholds = function(H2, Re, CV) {
  passes_H2 = H2 >= 0.1
  passes_Re = Re >= 0.4
  passes_CV = CV <= 40.5

  list(
    passes_H2 = passes_H2,
    passes_Re = passes_Re,
    passes_CV = passes_CV,
    passes_all = passes_H2 && passes_Re && passes_CV
  )
}