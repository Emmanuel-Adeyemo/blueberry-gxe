# 09/25

library(asreml)
library(dplyr)

fit_fa_asreml = function(fa_dta, Ainv, maxiter = 300) {
  # dta = fa_dta %>% mutate(Genotype = factor(Genotype), Env = factor(Env), Block = factor(Block))
  dta = prep_and_scale_dta(fa_dta)

  fit = asreml(
    fixed = y_scaled ~ Env,
    random = ~fa(Env, k = 1):vm(Genotype, Ainv) + Env:Block,
    residual = ~dsum(~units | Env),
    data = dta,
    maxit = maxiter,
    ai.loadings = 5,   
    ai.penalty = 15     
  )
  fit$call$maxit = maxiter
  update.asreml(fit)
}



fit_fa2_asreml = function(fa_dta, Ainv, maxiter = 300) {
  dta = prep_and_scale_dta(fa_dta)

  fit = asreml(
    fixed = y_scaled ~ Env,
    random = ~fa(Env, k = 2):vm(Genotype, Ainv) + Env:Block,
    residual = ~dsum(~units | Env),
    data = dta,
    maxit = maxiter,
    ai.loadings = 5,   
    ai.penalty = 15     
  )
  fit$call$maxit = maxiter
  update.asreml(fit)
  
}


fit_fa3_asreml = function(fa_dta, Ainv, maxiter = 300) {
  dta = prep_and_scale_dta(fa_dta)

  fit = asreml(
    fixed = y_scaled ~ Env,
    random = ~fa(Env, k = 3):vm(Genotype, Ainv) + Env:Block,
    residual = ~dsum(~units | Env),
    data = dta,
    maxit = maxiter,
    ai.loadings = 5,   
    ai.penalty  = 15     
  )
  fit$call$maxit = maxiter
  update.asreml(fit)
}

fit_cs_asreml = function(cs_dta, Ainv){

  dta = cs_dta %>% mutate(Genotype = factor(Genotype), Env = factor(Env), Block = factor(Block))

  fit = asreml(
    fixed = y ~ Env,
    # random = ~corh(Env):vm(Genotype, Ainv) + Env:Block,
    random = ~ vm(Genotype, Ainv) + idv(Env):vm(Genotype, Ainv) + Env:Block,
    residual = ~idv(units),
    data = dta
  )
}




extract_rg_matrix_asreml = function(varcomp, env_levels, n_factors) {
  var_rows = varcomp[grep("!var$", rownames(varcomp)), , drop = FALSE]
  if (nrow(var_rows) == 0) {
    stop(
      "extract_rg_matrix_asreml(): no rows matching '!var' "
    )
  }
  R = var_rows[, 1]

  loadings = lapply(seq_len(n_factors), function(k) {
    rows = varcomp[grep(paste0("!fa", k, "$"), rownames(varcomp)), , drop = FALSE]
    if (nrow(rows) == 0) {
      stop(
        "extract_rg_matrix_asreml(): no rows matching '!fa", k, "'",   
      )
    }
    rows[, 1]
  })

  L = do.call(cbind, loadings)
  # rownames(L) = env_levels 
  # print('here3')
  cov = L %*% t(L) + diag(R)
  dimnames(cov) = list(env_levels, env_levels)

  list(cov = cov, corr = stats::cov2cor(cov), loadings = L)
}