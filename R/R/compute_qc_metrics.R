# 09/18

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