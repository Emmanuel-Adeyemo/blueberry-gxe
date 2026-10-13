# 10/07


library(tidyr)


run_lrt = function(df) {
  df %>%
    mutate(
      LRT_stat = 2 * (LogL - lag(LogL)),
      df_diff = N_params - lag(N_params),
      p_value = pchisq(LRT_stat, df = df_diff, lower.tail = FALSE),
      Signif = case_when(
        p_value < 0.001 ~ '***',
        p_value < 0.01  ~ '**',
        p_value < 0.05  ~ '*',
        TRUE            ~ 'ns'
      )
    )
}

 # this is easier than extracting and transforming from coef$random
extract_blups_asreml = function(fit, dta, Ainv) {
  pv = predict(fit, classify = 'Env:Genotype', data = dta)$pvals
  env_fixed = predict(fit, classify = 'Env', data = dta)$pvals
  pv %>%
    mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
    left_join(
      env_fixed %>% transmute(Env = as.character(Env), env_effect = predicted.value),
      by = 'Env'
    ) %>%
    transmute(Genotype, Env, blup = predicted.value - env_effect, se = std.error)
}


extract_blups_sommer = function(fit) {

  fa_name = grep('fam\\(', names(fit$uList), value = TRUE)
  # sanity check, fa_name should be 1
  stopifnot(length(fa_name) == 1)
  U = fit$uList[[fa_name]]                   

  pev_list = fit$uPevList
  pev_name = grep('fam\\(', names(pev_list), value = TRUE)
  stopifnot(length(pev_name) == 1)
  PEV = pev_list[[pev_name]]
 

  stopifnot(identical(dim(U), dim(PEV)), identical(rownames(U), rownames(PEV)))

  # blups already centered
  U_df = as.data.frame(U) %>%
    tibble::rownames_to_column('Genotype') %>%
    tidyr::pivot_longer(-Genotype, names_to = 'Env', values_to = 'blup')

  PEV_df = as.data.frame(PEV) %>%
    tibble::rownames_to_column('Genotype') %>%
    tidyr::pivot_longer(-Genotype, names_to = 'Env', values_to = 'pev')

  dplyr::left_join(U_df, PEV_df, by = c('Genotype', 'Env')) %>%
    dplyr::mutate(se = sqrt(pmax(pev, 0))) %>%
    dplyr::select(Genotype, Env, blup, se)
}




add_reliability = function(blup_long, cov) {
  gv = tibble(Env = rownames(cov), gen_var = diag(cov))
  blup_long %>%
    left_join(gv, by = 'Env') %>%
    mutate(pev = se^2, reliability = pmax(0, 1 - pev / gen_var))
}


back_transform_blups = function(blup_long, env_stats) {
  blup_long %>%
    left_join(env_stats, by = 'Env"') %>%
    mutate(blup_original = blup * sd_env, se_original = se * sd_env)
}


# implementation of FAST analysis proposed by Smith & Cullis 2018 - OP and RMSD

fa_stability = function(blup_long, loadings, tested = NULL, engine = 'sommer') {
  U = blup_long %>%
    dplyr::select(Genotype, Env, blup) %>%
    pivot_wider(names_from = Env, values_from = blup)
  geno = U$Genotype
  U = as.matrix(U[, -1]); rownames(U) = geno

  env_levs = levels(factor(blup_long$Env))
  # asreml has no rownames by default
  if(engine == 'asreml') rownames(loadings) = env_levs

  stopifnot(setequal(colnames(U), rownames(loadings)))
  loadings = loadings[colnames(U), , drop = FALSE]

  s = svd(loadings)
  L1 = s$u[, 1] * s$d[1]
  if (mean(L1) < 0) L1 = -L1
  names(L1) = colnames(U)
  var_expl = s$d^2 / sum(s$d^2)

  tol = 0.05 * max(abs(L1))
  if (any(L1 < -tol)) {
    warning('Not all PC1 loadings are positive (', paste(names(L1)[L1 < -tol], collapse = ', '),'). 
    Smith & Cullis" OP assumes they are positive or a few very small negatives.')
  }

  mask = matrix(TRUE, nrow(U), ncol(U), dimnames = dimnames(U))
  if (!is.null(tested)) {
    mask[] = FALSE
    ij = cbind(match(tested$Genotype, geno), match(tested$Env, colnames(U)))
    ij = ij[complete.cases(ij), , drop = FALSE]
    mask[ij] = TRUE
  }

  se_w = blup_long %>% dplyr::select(Genotype, Env, se) %>%
    pivot_wider(names_from = Env, values_from = se)
  SE = as.matrix(se_w[match(geno, se_w$Genotype), colnames(U)])

  # factor-1 score - f1 is same as latent-regression slope
  f1 = as.vector(U %*% L1) / sum(L1^2)     
  tab = tibble(
    Genotype = geno,
    f1 = f1,
    # OP - value of the line at the mean loading
    OP = mean(L1) * f1,                          
    n_tested = rowSums(mask),

    # not just tested but also predicted
    RMSD = sapply(seq_along(geno), function(i) {
      sqrt(mean((U[i, ] - L1 * f1[i])^2))
    }),
    # RMSD = sapply(seq_along(geno), function(i) {
    #   j <- mask[i, ]
    #   if (!any(j)) return(NA_real_)
    #   sqrt(mean((U[i, j] - L1[j] * f1[i])^2))
    # }),
    # mean_se_tested = sapply(seq_along(geno), function(i) mean(SE[i, mask[i, ]]))
    mean_se_tested = sapply(seq_along(geno), function(i) mean(SE))
  ) 

  op_cut = mean(tab$OP)
  rmsd_cut = median(tab$RMSD, na.rm = TRUE)
  tab = tab %>%
    mutate(group = case_when(
      OP >= op_cut & RMSD <= rmsd_cut ~ 'Above average & stable',
      OP >= op_cut & RMSD > rmsd_cut ~ 'Above average, unstable',
      OP < op_cut & RMSD <= rmsd_cut ~ 'Below average, stable',
      TRUE ~ 'Below average, unstable'
    )) %>%
    arrange(desc(OP))

  list(table = tab, pc1_loadings = L1, var_explained = var_expl,
       op_cut = op_cut, rmsd_cut = rmsd_cut)
}




# PC loadings for biplot plotting
fa_pc_scores = function(blup_long, loadings, tested = NULL, n_pc = 2, engine = 'sommer') {

  U = blup_long %>% dplyr::select(Genotype, Env, blup) %>%
    pivot_wider(names_from = Env, values_from = blup)
  geno = U$Genotype
  U = as.matrix(U[, -1]); rownames(U) = geno

  env_levs = levels(factor(blup_long$Env))
  # asreml has no rownames by default
  if(engine == 'asreml') rownames(loadings) = env_levs

  stopifnot(setequal(colnames(U), rownames(loadings)), ncol(loadings) >= n_pc)
  loadings = loadings[colnames(U), , drop = FALSE]

  s = svd(loadings)
  Lpc = (s$u %*% diag(s$d, nrow = length(s$d)))[, seq_len(n_pc), drop = FALSE]  # rotated loadings
  
  sgn = sapply(seq_len(n_pc), function(k)
    if (k == 1) sign(mean(Lpc[, 1])) else sign(Lpc[which.max(abs(Lpc[, k])), k]))
  sgn[sgn == 0] = 1
  Lpc = sweep(Lpc, 2, sgn, `*`)
  dimnames(Lpc) = list(colnames(U), paste0('PC', seq_len(n_pc)))

  # PC loading cols orthogonal, geno score is a simple projection
  # score_k = sum_j(loading_jk * BLUP_ij) / sum_j(loading_jk^2)
  scores = U %*% Lpc %*% diag(1 / colSums(Lpc^2), nrow = n_pc)
  colnames(scores) = colnames(Lpc)

  n_t = if (is.null(tested)) rep(ncol(U), nrow(U)) else {
    cnt = tested %>% filter(Env %in% colnames(U)) %>% count(Genotype)
    ifelse(is.na(cnt$n[match(geno, cnt$Genotype)]), 0L, cnt$n[match(geno, cnt$Genotype)])
  }
  list(scores = as_tibble(scores) %>% mutate(Genotype = geno, n_tested = n_t, .before = 1),
       loadings = Lpc,
       var_explained = (s$d^2 / sum(s$d^2))[seq_len(n_pc)])   
}