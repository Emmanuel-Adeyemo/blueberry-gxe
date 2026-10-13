# 09/18

library(dplyr)
library(ggplot2)
library(tidyr)


.bloom_auc_point_sequence = function(dates, values, window_start, window_end) {
  ord = order(dates)
  dates = dates[ord]
  values = values[ord]

  if (dates[1] > window_start) {
    dates = c(window_start, dates)
    values = c(0, values)
  }
  last_i = length(dates)
  if (dates[last_i] < window_end) {
    dates = c(dates, window_end)
    values = c(values, values[last_i])
  }

  tibble(obs_date = dates, bloom_percent = values)
}



plot_bloom_auc = function(bloom_dta, plant_ids) {
  season = unique(bloom_dta$season)
  if (length(season) != 1) {
    stop('plot_bloom_auc(): bloom_dta must be a single environment (one season), got: ', paste(season, collapse = ', '))
  }

  window_start = season_anchor_date(season)
  window_end = max(bloom_dta$obs_date)

  selected = bloom_dta %>% filter(PlantID %in% plant_ids)
  if (nrow(selected) == 0) {
    stop('plot_bloom_auc(): none of the requested plant_ids were found in bloom_dta.')
  }

  point_seq = selected %>%
    group_by(PlantID) %>%
    group_modify(function(rows, ...) {
      raw_seq = .bloom_auc_point_sequence(rows$obs_date, rows$bloom_percent_raw, window_start, window_end) %>%
        mutate(series = 'Raw (frost-affected)')
      corrected_seq = .bloom_auc_point_sequence(rows$obs_date, rows$bloom_percent, window_start, window_end) %>%
        mutate(series = 'Corrected (carry-forward)')
      bind_rows(raw_seq, corrected_seq)
    }) %>%
    ungroup()

  auc_labels = compute_bloom_auc(bloom_dta) %>%
    filter(PlantID %in% plant_ids) %>%
    mutate(
      label = sprintf(
        'Corrected AUC: %.0f%%-days\nRaw AUC: %.0f%%-days\nFrost impact: %.0f%%-days',
        auc_corrected_percent, auc_raw_percent, frost_impact_percent
      )
    )

  ggplot(point_seq, aes(x = obs_date, y = bloom_percent, fill = series)) +
    geom_area(position = 'identity', alpha = 0.35, outline.type = 'upper') +
    geom_line(aes(color = series), linewidth = 0.9) +
    facet_wrap(~PlantID, ncol = 2) +
    scale_fill_manual(values = c('Raw (frost-affected)' = '#D55E00', 'Corrected (carry-forward)' = '#0072B2')) +
    scale_color_manual(values = c('Raw (frost-affected)' = '#D55E00', 'Corrected (carry-forward)' = '#0072B2')) +
    labs(
      title = 'Bloom progression and area under the curve: raw vs. frost-corrected',
      x = 'Date', y = 'Bloom (%)', fill = NULL, color = NULL
    ) +
    geom_text(
      data = auc_labels, aes(x = window_start, y = 70, label = .data$label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, color = 'grey30', lineheight = 1.1
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = 'top', plot.title = element_text(face = 'bold'))
}
