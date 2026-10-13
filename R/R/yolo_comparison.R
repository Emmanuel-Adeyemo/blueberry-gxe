# 09/23

library(dplyr)


run_loyo_comparison = function(covariates_df, bloom_pheno_df, predictor_col = "Chill_Portions_sum") {
  combined = covariates_df %>%
    inner_join(bloom_pheno_df, by = c("location", "season")) %>%
    filter(!is.na(.data[[predictor_col]]), !is.na(.data$pheno))

  n = nrow(combined)
  if (n < 3) {
    stop(
      "run_loyo_comparison(): we need at least 3 envs, got ", n
    )
  }

  purrr::map_dfr(seq_len(n), function(i) {
    train = combined[-i, ]
    test = combined[i, ]
    fit = lm(stats::reformulate(predictor_col, response = "pheno"), data = train)
    pred = as.numeric(predict(fit, newdata = test))
    tibble(
      location = test$location, season = test$season,
      actual = test$pheno, predicted = pred, error = abs(pred - test$pheno)
    )
  })
}


summarise_loyo_comparison = function(loyo_result, label) {
  tibble(
    label = label, n = nrow(loyo_result),
    MAE = mean(loyo_result$error), RMSE = sqrt(mean(loyo_result$error^2)),
    max_error = max(loyo_result$error)
  )
}