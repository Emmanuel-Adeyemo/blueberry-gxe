# 09/24

library(dplyr)

smooth_critical_window = function(pls_summary, window = 7) {
  if (window %% 2 == 0) { # has to do with centering of av, odd just trips it off
    stop("smooth_critical_window(): `window` must be odd, got ", window)
  }
  half = (window - 1) / 2

  .rolling_mean_centered = function(x, half) {
    n = length(x)
    out = rep(NA_real_, n)
    for (i in seq_len(n)) {
      lo = i - half
      hi = i + half
      if (lo >= 1 && hi <= n) out[i] = mean(x[lo:hi], na.rm = TRUE)
    }
    out
  }

  pls_summary %>%
    arrange(.data$Type, .data$JDay) %>%
    group_by(.data$Type) %>%
    mutate(
      VIP_smoothed = .rolling_mean_centered(.data$VIP, half),
      Coef_smoothed = .rolling_mean_centered(.data$Coef, half)
    ) %>%
    ungroup()
}



extract_critical_window = function(smoothed_pls_summary, vip_threshold = 0.8, min_run_length = 5) {
  flagged = smoothed_pls_summary %>%
    filter(!is.na(.data$VIP_smoothed)) %>%
    arrange(.data$Type, .data$JDay) %>%
    group_by(.data$Type) %>%
    mutate(
      above = .data$VIP_smoothed >= vip_threshold,
      run_id = cumsum(.data$above != lag(.data$above, default = FALSE))
    ) %>%
    ungroup() %>%
    filter(.data$above)

  
  # issue bcos dplyr just throws a warning with min/max when there are 0 vip days instead of just straight error
  if (nrow(flagged) == 0) {
    return(tibble(
      Type = character(0), start_JDay = numeric(0), end_JDay = numeric(0),
      start_Date = numeric(0), end_Date = numeric(0),
      length_days = integer(0), mean_VIP_smoothed = numeric(0)
    ))
  }

  flagged %>%
    group_by(.data$Type, .data$run_id) %>%
    summarise(
      start_JDay = min(.data$JDay), end_JDay = max(.data$JDay),
      start_Date = .data$Date[which.min(.data$JDay)], end_Date = .data$Date[which.max(.data$JDay)],
      length_days = n(), mean_VIP_smoothed = mean(.data$VIP_smoothed),
      .groups = "drop"
    ) %>%
    filter(.data$length_days >= min_run_length) %>%
    dplyr::select(-any_of("run_id")) %>%
    arrange(.data$Type, .data$start_JDay)
}