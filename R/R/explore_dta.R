# 09/11

library(dplyr)
library(ggplot2)
library(tidyr)


#' @description general-purpose exploratory data analysis helper,  summarizes one trait's descriptive statistics. not wired into targets
#' @param dta any df containing `trait_col` e.g. a raw or cleaned trait long df such as `al_yield_long_clean`/`nc_yield_clean` or a per-trial qc df such as `per_trial_qc`
#' @param group_cols character vec of col names to group by e.g. `c('season', 'location')`
#' @param trait_col name of the numeric trait col to summarize
#' @output df 
describe_trait_by_group = function(dta, trait_col, group_cols) {
  dta %>%
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


#' @description general purpose EDA helper: for each combination of grouping columns, counts distinct genotypes and summarizes replication. not wired into targets
#' @param data Any df containing `genotype_col` and all of `group_cols`, e.g. a trait-long or per-trial data frame such as `al_yield_long_clean`/`nc_yield_clean`, or an
#'   FA-model-ready df such as `prep_yield`, output of `prepare_fa_dta()`
#' @param group_cols character vec of col names to group by (e.g. `c('season', 'location')` or `'Env'`
#' @param genotype_col name of the genotype cols, def to `'genotype'` 
#' @output df 
describe_genotype_coverage = function(dta, group_cols, genotype_col = 'genotype') {
  per_genotype_counts = dta %>%
    group_by(across(all_of(c(group_cols, genotype_col)))) %>%
    summarise(n_reps = n(), .groups = 'drop')

  per_genotype_counts %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_genotypes = n(),
      n_rows = sum(n_reps),
      mean_reps_per_genotype = mean(n_reps),
      min_reps_per_genotype = min(n_reps),
      max_reps_per_genotype = max(n_reps),
      .groups = 'drop'
    )
}


#' @description general-purpose EDA helper, tallies total and missing obs per trait within each grouping columns, and reshapes to one row per (group, trait). 
#' @param trait_cols character vec of trait col names to check for missingness.
#' @param group_cols character vec of col names to group by (e.g. `c('season', 'location')`).
#' @output df
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
    mutate(pct_missing = 100 * (n_missing / n_total))
}


#' @description general-purpose EDA helper: boxplot + jittered points of one trait across levels of one grouping col. not wired into targets
#' @param data Any df containing `trait_col` and `group_col`.
#' @param trait_col name of the numeric trait col to plot on the y-axis.
#' @param group_col name of the col to plot on the x-axis, e.g. `'Env'` or `'season'`.
#' @param title optional plot title; defaults to `"<trait_col> by <group_col>"` if `NULL`.
#' @output ggplot obj
plot_trait_by_group = function(dta, trait_col, group_col, title = NULL) {
  if (is.null(title)) title = paste(trait_col, 'by', group_col)

  ggplot(dta, aes(x = .data[[group_col]], y = .data[[trait_col]])) +
    geom_boxplot(outlier.shape = NA, fill = 'grey90') +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    labs(title = title, x = group_col, y = trait_col) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}