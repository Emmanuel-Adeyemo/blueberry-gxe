# 09/24

library(dplyr)
library(sommer)

# mmes uses henderson's mixed model eq  faster with sparse rel mat. more than 50% of amat is 0

# Henderson == T need A inverse

# prep a inverse mat for sommer
prep_ainv_sommer = function(Ainv) {
  # stopifnot(!is.null(rownames(Ainv)), all(genotypes %in% rownames(Ainv)))
  Ainv = as(as(as(Ainv, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  attr(Ainv, "inverse") = TRUE
  Ainv
}

# not needed anymore - could use for base model insteead of CS
fit_diag_sommer = function(fa_dta, Ainv, nIters = 150, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef))
  mmes(
    fixed_formula,
    henderson = TRUE,   
    random = ~ vsm(dsm(Env), ism(Genotype), Gu = Ainv) + vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE
  )
}



fit_fa_sommer = function(fa_dta, Ainv, n_factors, nIters = 150, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef))
  mmes(
    fixed_formula,
    henderson = TRUE,
    random = ~ vsm(fam(Env, k = n_factors), ism(Genotype), Gu = Ainv) + vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE, computeCi=2

  )
}



# keep for now
# fit_fa_chill_sommer <- function(data, Ainv, n_factors, nIters = 150) {
#   dta = prep_and_scale_dta(data)
#   mmes(
#     y_scaled ~ Env + Age,
#     henderson = TRUE,
#     random = ~ vsm(fam(Env, k = n_factors), ism(Genotype), Gu = Ainv) +  # uns GxE
#                vsm(ism(Chill_c), ism(Genotype), Gu = Ainv) +             # chill sense slope
#                vsm(ism(Env:Block)),
#     rcov = ~ vsm(dsm(Env), ism(units)),
#     nIters = nIters, data = dta, verbose = FALSE,
#     getPEV = TRUE
#   )
# }

extract_rg_matrix_sommer = function(fa_fit) {
  info = loadings_mmes(fa_fit)
  Sigma = info$sigma2 * (info$loadings %*% t(info$loadings) + diag(info$specific))
  envs = rownames(info$loadings)
  dimnames(Sigma) = list(envs, envs)
  list(cov = Sigma, corr = cov2cor(Sigma),
       loadings = info$loadings, specific = info$specific, sigma2 = info$sigma2)
}


# 10/24 
fit_fa3_min_ide_sommer = function(fa_dta, Ainv, Dinv, ka = 3, nIters = 150, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef))
  mmes(
    fixed_formula,
    henderson = TRUE,
    random = ~ vsm(fam(Env, k = ka), ism(Genotype), Gu = Ainv) + vsm(ism(Genotype), Gu = Dinv) + vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE, computeCi = 2
  )
}


fit_fa3_min_sommer = function(fa_dta, Ainv, ka = 3, nIters = 150, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef))
  mmes(
    fixed_formula,
    henderson = TRUE,
    random = ~ vsm(fam(Env, k = ka), ism(Genotype), Gu = Ainv) +
               vsm(dsm(Env), ism(Genotype)) + vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE, computeCi = 2
  )
}


extract_nonadd_scalar_sommer = function(fit) {
  nms = names(fit$uList)
  term = nms[grepl("ism\\(Genotype\\)", nms) & !grepl("fam\\(|Gu", nms)]
  stopifnot(length(term) == 1)

  g_e = fit$uList[[term]]
  pev = fit$uPevList[[term]]
  stopifnot(identical(dim(g_e), dim(pev)), identical(rownames(g_e), rownames(pev)))

  g_e_df = as.data.frame(g_e) %>%
    tibble::rownames_to_column("Genotype") %>%
    tidyr::pivot_longer(-Genotype, names_to = "Env", values_to = "g_e")

  pev_df = as.data.frame(pev) %>%
    tibble::rownames_to_column("Genotype") %>%
    tidyr::pivot_longer(-Genotype, names_to = "Env", values_to = "pev")

  dplyr::left_join(g_e_df, pev_df, by = c("Genotype", "Env")) %>%
    dplyr::mutate(se_e = sqrt(pmax(pev, 0))) %>%
    dplyr::select(Genotype, Env, g_e, se_e)
}


fit_fa3_add_nonadd_sommer = function(fa_dta, Ainv, ka = 3, ke = 1, nIters = 150, extra_fixed = NULL) {
  dta = prep_and_scale_dta(fa_dta)

  fixed_ef = 'Env'
  if(!is.null(extra_fixed)) fixed_ef = paste(fixed_ef, extra_fixed, sep = '+')
  fixed_formula = as.formula(paste0('y_scaled ~ ', fixed_ef))
  mmes(
    fixed_formula,
    henderson = TRUE,
    random = ~ vsm(fam(Env, k = ka), ism(Genotype), Gu = Ainv) + vsm(fam(Env, k = ke), ism(Genotype)) + vsm(ism(Env:Block)),
    rcov = ~ vsm(dsm(Env), ism(units)),
    nIters = nIters, data = dta, verbose = FALSE,
    getPEV = TRUE, computeCi = 2
  )
}


find_fa_terms_sommer = function(fit) {
  nms = names(fit$covStruct)
  fa_terms = nms[grepl("fam\\(", nms)]
  stopifnot(length(fa_terms) == 2)
  has_gu = grepl("Gu", fa_terms)
  stopifnot(sum(has_gu) == 1)
  list(additive = fa_terms[has_gu], nonadditive = fa_terms[!has_gu])
}


extract_rg_matrix_sommer_split = function(fa_fit) {
  terms = find_fa_terms_sommer(fa_fit)
  build = function(term) {
    info = loadings_mmes(fa_fit, term = term)
    Sigma = info$sigma2 * (info$loadings %*% t(info$loadings) + diag(info$specific))
    envs = rownames(info$loadings)
    dimnames(Sigma) = list(envs, envs)
    list(cov = Sigma, corr = cov2cor(Sigma),
         loadings = info$loadings, specific = info$specific, sigma2 = info$sigma2)
  }
  list(additive = build(terms$additive), nonadditive = build(terms$nonadditive))
}



extract_blups_sommer_split = function(fit) {
  terms = find_fa_terms_sommer(fit)

  build = function(term) {
    U = fit$uList[[term]]
    PEV = fit$uPevList[[term]]
    stopifnot(identical(dim(U), dim(PEV)), identical(rownames(U), rownames(PEV)))

    U_df = as.data.frame(U) %>%
      tibble::rownames_to_column("Genotype") %>%
      tidyr::pivot_longer(-Genotype, names_to = "Env", values_to = "blup")

    PEV_df = as.data.frame(PEV) %>%
      tibble::rownames_to_column("Genotype") %>%
      tidyr::pivot_longer(-Genotype, names_to = "Env", values_to = "pev")

    dplyr::left_join(U_df, PEV_df, by = c("Genotype", "Env")) %>%
      dplyr::mutate(se = sqrt(pmax(pev, 0))) %>%
      dplyr::select(Genotype, Env, blup, se)
  }

  list(additive = build(terms$additive), nonadditive = build(terms$nonadditive))
}





find_fa_pooled_terms_sommer = function(fit) {
  nms = names(fit$covStruct)
  fa_term = nms[grepl("fam\\(", nms)]
  stopifnot(length(fa_term) == 1)
  pooled_term = nms[grepl("ism\\(Genotype\\)", nms) & !grepl("fam\\(", nms) & grepl("Gu", nms)]
  stopifnot(length(pooled_term) == 1)
  list(additive = fa_term, nonadditive = pooled_term)
}



extract_nonadd_scalar_pooled_sommer = function(fit) {
  terms = find_fa_pooled_terms_sommer(fit)
  g_e = fit$uList[[terms$nonadditive]]
  pev = fit$uPevList[[terms$nonadditive]]
  stopifnot(identical(dim(g_e), dim(pev)))

  tibble::tibble(
    Genotype = rownames(g_e),
    g_e = as.numeric(g_e[, 1]),
    se_e = sqrt(pmax(as.numeric(pev[, 1]), 0))
  )
}