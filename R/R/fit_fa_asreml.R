# 09/25

library(asreml)
library(dplyr)


fit_fa_asreml = function(fa_dta, Gmat, ka, maxiter = 300, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef ))

  # print(fixed_formula)

  fit = asreml(
    fixed = fixed_formula,
    random = as.formula(paste0("~fa(Env, k = ", ka, "):vm(Genotype, Gmat) + Env:Block")),
    residual = ~dsum(~units | Env),
    data = dta,
    maxit = maxiter,
    ai.loadings = 5,   
    ai.penalty = 15     
  )
  fit$call$maxit = maxiter
  update.asreml(fit)
}


fit_cs_asreml = function(cs_dta, Gmat, extra_fixed = NULL){

  dta = cs_dta %>% mutate(Genotype = factor(Genotype), Env = factor(Env), Block = factor(Block))

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef ))

  fit = asreml(
    fixed = fixed_formula,
    # random = ~corh(Env):vm(Genotype, Gmat) + Env:Block,
    random = ~ vm(Genotype, Gmat) + idv(Env):vm(Genotype, Gmat) + Env:Block,
    residual = ~idv(units),
    data = dta
  )
}



# 10/25 - pooled non additive
fit_fa_ide_asreml = function(fa_dta, Gmat, maxiter = 300, Dmat, ka, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef ))

  fit = asreml(
    fixed = fixed_formula,
    random = as.formula(paste0("~fa(Env, k = ", ka, "):vm(Genotype, Gmat) + vm(Genotype, Dmat) + Env:Block")),
    residual = ~dsum(~units | Env),
    data = dta,
    maxit = maxiter,
    ai.loadings = 5,
    ai.penalty = 15
  )
  fit$call$maxit = maxiter
  update.asreml(fit)
}

# non add diag only
fit_fa_diag_asreml = function(fa_dta, Gmat, maxiter = 300, Dmat, ka, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef ))

  fit = asreml(
    fixed = fixed_formula,
    random = as.formula(paste0("~fa(Env, k = ", ka, "):vm(Genotype, Gmat) + diag(Env):vm(Genotype, Dmat) + Env:Block")),
    residual = ~dsum(~units | Env),
    data = dta,
    maxit = maxiter,
    ai.loadings = 5,
    ai.penalty = 15
  )
  fit$call$maxit = maxiter
  update.asreml(fit)
}

# 10/14 - additive + non-additive split with fa, per Smith et al 2025 eq. 4-9.
# did not work out, I dont need the two way split descriped in the article
#

fit_fa_ke1_asreml = function(fa_dta, Gmat, maxiter = 300, Dmat, ka, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef ))

  fit = asreml(
    fixed = fixed_formula,
    random = as.formula(paste0("~fa(Env, k = ", ka, "):vm(Genotype, Gmat) + fa(Env, k = 1):vm(Genotype, Dmat) + Env:Block")),
    residual = ~dsum(~units | Env),
    data = dta,
    maxit = maxiter,
    ai.loadings = 5,
    ai.penalty = 15
  )
  fit$call$maxit = maxiter
  update.asreml(fit)
}


# the non-additive term is unique by envt

extract_nonadd_scalar_asreml = function(fit, dta, Gmat, Dmat) {
  predict(fit, classify = "Env:Genotype", data = dta,
          onlyuse = ~diag(Env):vm(Genotype, Dmat))$pvals %>%
    mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
    transmute(Genotype, Env, g_e = predicted.value, se_e = std.error)
}

# extract non add term - one value per genotype across envts
extract_nonadd_scalar_pooled_asreml = function(fit, dta, Gmat, Dmat) {
  predict(fit, classify = "Genotype", data = dta,
          onlyuse = ~vm(Genotype, Dmat))$pvals %>%
    mutate(Genotype = as.character(Genotype)) %>%
    transmute(Genotype, g_e = predicted.value, se_e = std.error)
}


# gets loadings, corr, and cov from asreml fit
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