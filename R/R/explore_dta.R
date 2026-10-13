# 09/11

library(dplyr)
library(ggplot2)
library(tidyr)



describe_trait_by_group = function(data, trait_col, group_cols) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n = sum(!is.na(.data[[trait_col]])),
      n_missing = sum(is.na(.data[[trait_col]])),
      mean = mean(.data[[trait_col]], na.rm = TRUE),
      sd = sd(.data[[trait_col]], na.rm = TRUE),
      se = .data$sd / sqrt(.data$n),
      cv_pct = ifelse(.data$mean == 0, NA_real_, 100 * .data$sd / .data$mean),
      min = suppressWarnings(min(.data[[trait_col]], na.rm = TRUE)),
      median = median(.data[[trait_col]], na.rm = TRUE),
      max = suppressWarnings(max(.data[[trait_col]], na.rm = TRUE)),
      .groups = 'drop'
    )
}


describe_genotype_coverage = function(data, group_cols, genotype_col = 'genotype') {
  per_genotype_counts = data %>%
    group_by(across(all_of(c(group_cols, genotype_col)))) %>%
    summarise(n_reps = n(), .groups = 'drop')

  per_genotype_counts %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_genotypes = n(),
      n_rows = sum(.data$n_reps),
      mean_reps_per_genotype = mean(.data$n_reps),
      min_reps_per_genotype = min(.data$n_reps),
      max_reps_per_genotype = max(.data$n_reps),
      .groups = 'drop'
    )
}


summarise_missingness = function(data, trait_cols, group_cols) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(
        all_of(trait_cols),
        list(n_total = ~ length(.x), n_missing = ~ sum(is.na(.x))),
        .names = '{.col}_{.fn}"'
      ),
      .groups = 'drop'
    ) %>%
    pivot_longer(
      cols = -all_of(group_cols),
      names_to = c('trait', '.value'),
      names_sep = '_'
    ) %>%
    mutate(pct_missing = 100 * .data$n_missing / .data$n_total)
}


plot_trait_by_group = function(data, trait_col, group_col, title = NULL) {
  if (is.null(title)) title = paste(trait_col, 'by', group_col)

  ggplot(data, aes(x = .data[[group_col]], y = .data[[trait_col]])) +
    geom_boxplot(outlier.shape = NA, fill = 'grey90') +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    labs(title = title, x = group_col, y = trait_col) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}