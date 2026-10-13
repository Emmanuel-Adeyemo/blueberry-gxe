# 09/24

library(dplyr)
library(sommer)

# mmes uses henderson's mixed model eq  faster with sparse rel mat. more than 50% of amat is 0
# fit diag model to ptoduce H mat

# Henderson == T need A inverse

# prep a inverse mat for sommer
prep_ainv_sommer = function(Ainv) {
  # stopifnot(!is.null(rownames(Ainv)), all(genotypes %in% rownames(Ainv)))
  Ainv = as(as(as(Ainv, 'dMatrix'), 'generalMatrix'), 'CsparseMatrix')
  attr(Ainv, 'inverse') = TRUE
  Ainv
}

# not needed anymore - may be I can use for base model insteead of CS
fit_diag_sommer = function(fa_dta, Ainv, nIters = 150) {
  dta = prep_and_scale_dta(fa_dta)
  mmes(
    y_scaled ~ Env,
    henderson = TRUE,   
    random = ~ vsm(dsm(Env), ism(Genotype), Gu = Ainv) + vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE
  )
}



fit_fa_sommer = function(fa_dta, Ainv, n_factors, nIters = 150) {
  dta = prep_and_scale_dta(fa_dta)
  mmes(
    y_scaled ~ Env,
    henderson = TRUE,
    random = ~ vsm(fam(Env, k = n_factors), ism(Genotype), Gu = Ainv) + vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE, computeCi=2

  )
}

# keep for now
fit_fa_chill_sommer = function(data, Ainv, n_factors, nIters = 150) {
  dta = prep_and_scale_dta(data)
  mmes(
    y_scaled ~ Env,
    henderson = TRUE,
    random = ~ vsm(fam(Env, k = n_factors), ism(Genotype), Gu = Ainv) +  # uns GxE
               vsm(ism(Chill_c), ism(Genotype), Gu = Ainv) +             # chill sense slope
               vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE
  )
}

extract_rg_matrix_sommer = function(fa_fit) {
  info = loadings_mmes(fa_fit)
  Sigma = info$sigma2 * (info$loadings %*% t(info$loadings) + diag(info$specific))
  envs = rownames(info$loadings)
  dimnames(Sigma) = list(envs, envs)
  list(cov = Sigma, corr = cov2cor(Sigma),
       loadings = info$loadings, specific = info$specific, sigma2 = info$sigma2)
}

