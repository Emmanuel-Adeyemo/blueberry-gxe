# 09/23

library(ggplot2)
library(dplyr)


plot_critical_window_diagnostic = function(smoothed_pls_summary, extracted_windows = NULL, vip_threshold = 0.8) {
  df = smoothed_pls_summary %>%
    mutate(approx_date = as.Date(.data$JDay, origin = '2022-12-31'))

  p = ggplot(df, aes(x = .data$approx_date, color = .data$Type))

  if (!is.null(extracted_windows) && nrow(extracted_windows) > 0) {
    shade_df <- extracted_windows %>%
      mutate(
        xmin = as.Date(.data$start_JDay, origin = '2022-12-31'),
        xmax = as.Date(.data$end_JDay, origin = '2022-12-31')
      )
    p = p + geom_rect(
      data = shade_df, inherit.aes = FALSE,
      aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf, fill = .data$Type),
      alpha = 0.12
    )
  }

  p = p +
    geom_hline(yintercept = vip_threshold, linetype = 'dashed', color = 'grey40') +
    geom_line(aes(y = .data$VIP), linewidth = 0.3, alpha = 0.35) +
    geom_line(aes(y = .data$VIP_smoothed), linewidth = 1) +
    scale_color_manual(values = c(Chill = '#3B7FB6', Heat = '#D9534F')) +
    scale_fill_manual(values = c(Chill = '#3B7FB6', Heat = '#D9534F'), guide = "none") +
    scale_x_date(date_labels = '%b %d') +
    labs(
      title = 'Critical-window diagnostic: raw vs. smoothed VIP',
      subtitle = 'Faint = raw (noisy) daily VIP; bold = smoothed; shaded = extracted critical window; dashed = VIP 0.8 threshold',
      x = NULL, y = 'VIP score', color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = 'top')

  p
}


plot_fixed_vs_discovered = function(discovered_windows,
                                      fixed_chill_window = FIXED_CHILL_WINDOW,
                                      fixed_heat_window = FIXED_HEAT_WINDOW) {
  offset_to_date = function(offset) as.Date('2022-09-01') + offset

  fixed_df = tibble(
    Type = c('Chill', 'Heat'), Approach = 'Fixed convention',
    xmin = offset_to_date(c(fixed_chill_window['start'], fixed_heat_window['start'])),
    xmax = offset_to_date(c(fixed_chill_window['end'], fixed_heat_window['end']))
  )

  discovered_df = discovered_windows %>%
    mutate(
      Approach = 'PLS-discovered',
      xmin = as.Date(.data$start_JDay, origin = '2022-12-31'),
      xmax = as.Date(.data$end_JDay, origin = '2022-12-31')
    ) %>%
    dplyr::select('Type', 'Approach', 'xmin', 'xmax')

  combined = bind_rows(fixed_df, discovered_df) %>%
    mutate(row_label = paste(.data$Approach, .data$Type))

  ggplot(combined, aes(y = .data$row_label, color = .data$Type)) +
    geom_segment(aes(x = .data$xmin, xend = .data$xmax, yend = .data$row_label), linewidth = 6, alpha = 0.7) +
    scale_color_manual(values = c(Chill = '#3B7FB6', Heat = '#D9534F')) +
    scale_x_date(date_labels = '%b %d') +
    labs(
      title = 'Fixed convention vs. PLS-discovered windows',
      x = NULL, y = NULL, color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = 'top')
}




plot_loyo_comparison = function(loyo_summaries) {
  plot_df = loyo_summaries %>%
    dplyr::select('label', 'MAE', 'RMSE') %>%
    tidyr::pivot_longer(cols = c('MAE', 'RMSE'), names_to = 'metric', values_to = 'value')

  ggplot(plot_df, aes(x = .data$label, y = .data$value, fill = .data$metric)) +
    geom_col(position = 'dodge') +
    geom_text(aes(label = round(.data$value, 1)), position = position_dodge(width = 0.9), vjust = -0.4, size = 3.5) +
    labs(
      title = 'Leave-one-year-out prediction error by approach',
      subtitle = 'Lower is better -- error in days between predicted and actual bloom date',
      x = NULL, y = 'Days', fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = 'top', axis.text.x = element_text(angle = 20, hjust = 1))
}
