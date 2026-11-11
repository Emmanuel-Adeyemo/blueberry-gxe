# 09/13

library(dplyr)
library(asreml)

fit_single_trial_blue_asreml = function(trial_dta) {

  # has_block = 'block' %in% names(trial_dta) &&
  #   length(unique(trial_dta$block[!is.na(trial_dta$block)])) > 1

  dta = trial_dta %>% mutate(genotype = factor(genotype), factor(block))
  # if (has_block) dta = dta %>% mutate(block = factor(block))
  dta = as.data.frame(dta)

  # fixed_form = if (has_block) y ~ block + genotype else y ~ genotype

  fit = tryCatch(
    asreml(
      fixed = y ~ block + genotype,
      residual = ~units,
      na.action = na.method(x = 'include'),
      data = dta
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  pred = predict(fit, classify = 'genotype')$pvals

  pred %>%
    transmute(
      Genotype = as.character(genotype),
      predicted_value = predicted.value,
      std_error = std.error
    )
  
}