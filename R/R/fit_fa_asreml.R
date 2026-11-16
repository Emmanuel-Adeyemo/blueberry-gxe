# 09/25

library(asreml)
library(dplyr)


#' @description fits the core factor-analytic (FA) GxE model in ASReml-R. this is the main additive-only model used across all trait stages (yield, bloom, AUC, plant height, quality)
#'   passes relationship mat not the inverse.
#' @param fa_dta FA-model-ready long-format trait data (one row per plant/subplot observation) e.g. `prep_yield`, combined output of `prepare_fa_dta()`/`prepare_nf_fa_dta()` in prepare_fa_dta.R
#' @param Gmat additive genomic (or pedigree) relationship matrix, not inverted, `g_mat` loaded via `readRDS('data/processed/markers/G_mat.rds')`
#' @param ka num of additive FA factors to fit e.g. 1, 2, or 3 
#' @param maxiter max ASReml-R iterations, def to 300
#' @param extra_fixed optional additional fixed-effect term to append to `Env` in the fixed e.g. `factor(Age)`
#' @output fitted `asreml` model obj
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


# fit_cs_asreml = function(cs_dta, Gmat, extra_fixed = NULL){

#   dta = cs_dta %>% mutate(Genotype = factor(Genotype), Env = factor(Env), Block = factor(Block))

#   fixed_ef = 'Env'
#   if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
#   fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef ))

#   fit = asreml(
#     fixed = fixed_formula,
#     # random = ~corh(Env):vm(Genotype, Gmat) + Env:Block,
#     random = ~ vm(Genotype, Gmat) + idv(Env):vm(Genotype, Gmat) + Env:Block,
#     residual = ~idv(units),
#     data = dta
#   )
# }



# 10/25 - pooled non additive
#' @description fits fa model with a pooled non additive term (identity structured)
#' @param fa_dta FA-model-ready long-format trait data (one row per plant/subplot observation) e.g. `prep_yield`, combined output of `prepare_fa_dta()`/`prepare_nf_fa_dta()` in prepare_fa_dta.R
#' @param Gmat additive genomic (or pedigree) relationship matrix, not inverted, `g_mat` loaded via `readRDS('data/processed/markers/G_mat.rds')`. dim is geno x geno
#' @param Dmat dominance genomic relationship matrix, not inverted, `d_mat` loaded via `readRDS('data/processed/markers/D_mat.rds')` dim is geno x geno
#' @param ka num of additive FA factors to fit e.g. 1, 2, or 3 
#' @param maxiter max ASReml-R iterations, def to 300
#' @param extra_fixed optional additional fixed-effect term to append to `Env` in the fixed e.g. `factor(Age)`
#' @output fitted `asreml` model obj
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
#' @description fits fa model with a diagonal non additive term (envt specific but uncorrelated across envts)
#' @param fa_dta FA-model-ready long-format trait data (one row per plant/subplot observation) e.g. `prep_yield`, combined output of `prepare_fa_dta()`/`prepare_nf_fa_dta()` in prepare_fa_dta.R
#' @param Gmat additive genomic (or pedigree) relationship matrix, not inverted, `g_mat` loaded via `readRDS('data/processed/markers/G_mat.rds')`. dim is geno x geno
#' @param Dmat dominance genomic relationship matrix, not inverted, `d_mat` loaded via `readRDS('data/processed/markers/D_mat.rds')` dim is geno x geno
#' @param ka num of additive FA factors to fit e.g. 1, 2, or 3 
#' @param maxiter max ASReml-R iterations, def to 300
#' @param extra_fixed optional additional fixed-effect term to append to `Env` in the fixed e.g. `factor(Age)`
#' @output fitted `asreml` model obj
fit_fa_diag_asreml = function(fa_dta, Gmat, maxiter = 300, Dmat, ka, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef ))

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
#' @description fits fa model with a second 1-factor fa structure for the non additive term (envt specific + correlated across envts)
#' @param fa_dta FA-model-ready long-format trait data (one row per plant/subplot observation) e.g. `prep_yield`, combined output of `prepare_fa_dta()`/`prepare_nf_fa_dta()` in prepare_fa_dta.R
#' @param Gmat additive genomic (or pedigree) relationship matrix, not inverted, `g_mat` loaded via `readRDS('data/processed/markers/G_mat.rds')`. dim is geno x geno
#' @param Dmat dominance genomic relationship matrix, not inverted, `d_mat` loaded via `readRDS('data/processed/markers/D_mat.rds')` dim is geno x geno
#' @param ka num of additive FA factors to fit e.g. 1, 2, or 3 
#' @param maxiter max ASReml-R iterations, def to 300
#' @param extra_fixed optional additional fixed-effect term to append to `Env` in the fixed e.g. `factor(Age)`
#' @output fitted `asreml` model obj
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
#' @description extracts envt specific non additive effect prediction from the fitted diag non additive model
#' @param fit a model fitted with a diagonal non-additive term, e.g. `fit_fa_diag_asreml()`
#' @param dta same dta used for fit e.g. `prep_yield`, combined output of `prepare_fa_dta()`/`prepare_nf_fa_dta()` in prepare_fa_dta.R
#' @param Gmat additive genomic (or pedigree) relationship matrix, not inverted, `g_mat` loaded via `readRDS('data/processed/markers/G_mat.rds')`. dim is geno x geno
#' @param Dmat dominance genomic relationship matrix, not inverted, `d_mat` loaded via `readRDS('data/processed/markers/D_mat.rds')` dim is geno x geno
#' @output fitted `asreml` model obj
extract_nonadd_scalar_asreml = function(fit, dta, Gmat, Dmat) {
  predict(fit, classify = "Env:Genotype", data = dta,
          onlyuse = ~diag(Env):vm(Genotype, Dmat))$pvals %>%
    mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
    transmute(Genotype, Env, g_e = predicted.value, se_e = std.error)
}

# extract non add term - one value per genotype across envts
#' @description extracts pooled non additive effect prediction from the fitted identity non additive model
#' @param fit a model fitted with a pooled non-additive term, e.g. `fit_fa_ide_asreml()`
#' @param dta same dta used for fit e.g. `prep_yield`, combined output of `prepare_fa_dta()`/`prepare_nf_fa_dta()` in prepare_fa_dta.R
#' @param Gmat additive genomic (or pedigree) relationship matrix, not inverted, `g_mat` loaded via `readRDS('data/processed/markers/G_mat.rds')`. dim is geno x geno
#' @param Dmat dominance genomic relationship matrix, not inverted, `d_mat` loaded via `readRDS('data/processed/markers/D_mat.rds')` dim is geno x geno
#' @output fitted `asreml` model obj
extract_nonadd_scalar_pooled_asreml = function(fit, dta, Gmat, Dmat) {
  predict(fit, classify = "Genotype", data = dta,
          onlyuse = ~vm(Genotype, Dmat))$pvals %>%
    mutate(Genotype = as.character(Genotype)) %>%
    transmute(Genotype, g_e = predicted.value, se_e = std.error)
}




# gets loadings, corr, and cov from asreml fit
#' @description constructs the genetic covariance matrix and genetic correlation matrix from an ASReml-R FA fit's var component table
#' @param varcomp The `varcomp` component of `summary(fit)` for an FA model fitted with `fit_fa_asreml()` or similar, e.g. `summary(fit_fa3_asr)$varcomp` 
#' @param env_levels character vec of envt level names
#' @param n_factors num of FA factors the model was fitted with (i.e. the `ka` passed to `fit_fa_asreml()`).
#' @output list with `cov` (genetic covariance matrix, Env x Env), `corr` (genetic correlation  matrix, via `stats::cov2cor()`), and `loadings` (the Env x `n_factors` loading matrix `L`)

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