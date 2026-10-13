# 09/18

compute_qc_metrics = function(vc_g, vc_e, mean_y, pev) {
  H2 = vc_g / (vc_g + vc_e)
  Re = max(0, 1 - mean(pev, na.rm = TRUE) / vc_g)
  CV = (sqrt(vc_e) / mean_y) * 100

  list(H2 = H2, Re = Re, CV = CV)
}


apply_thresholds = function(H2, Re, CV) {
  passes_H2 = H2 >= 0.1
  passes_Re = Re >= 0.2
  passes_CV = CV <= 70.5

  list(
    passes_H2 = passes_H2,
    passes_Re = passes_Re,
    passes_CV = passes_CV,
    passes_all = passes_H2 && passes_Re && passes_CV
  )
}