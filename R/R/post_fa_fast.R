# 10/07


library(tidyr)
library(dplyr)
library(purrr)
library(tibble)


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


# 10/24 
#
split_varcomp_asreml = function(varcomp,
                                  add_pattern = "Gmat",
                                  nonadd_pattern = "ide\\(Genotype|Dmat") {
  rn = rownames(varcomp)
  is_add = grepl(add_pattern, rn)
  is_nonadd = grepl(nonadd_pattern, rn) & !is_add
  stopifnot(any(is_add), any(is_nonadd))
  list(additive = varcomp[is_add, , drop = FALSE],
       nonadditive = varcomp[is_nonadd, , drop = FALSE])
}


# 10/24 
extract_rg_matrix_asreml_split = function(varcomp, env_levels, ka, ke,
                                            add_pattern = "Gmat",
                                            nonadd_pattern = "ide\\(Genotype|Dmat") {
  parts = split_varcomp_asreml(varcomp, add_pattern = add_pattern, nonadd_pattern = nonadd_pattern)
  list(
    additive = extract_rg_matrix_asreml(parts$additive, env_levels, ka),
    nonadditive = extract_rg_matrix_asreml(parts$nonadditive, env_levels, ke)
  )
}


 # this is easier than extracting and transforming from coef$random
extract_blups_asreml_yld = function(fit, dta, Gmat) {
  pv = predict(fit, classify = "Env:Genotype", data = dta, average = 'factor(Age)')$pvals
  env_fixed = predict(fit, classify = "Env", data = dta, average = 'factor(Age)')$pvals
  pv %>%
    mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
    left_join(
      env_fixed %>% transmute(Env = as.character(Env), env_effect = predicted.value),
      by = "Env"
    ) %>%
    transmute(Genotype, Env, blup = predicted.value - env_effect, se = std.error)
}


extract_blups_asreml = function(fit, dta, Gmat) {
  pv = predict(fit, classify = "Env:Genotype", data = dta)$pvals
  env_fixed = predict(fit, classify = "Env", data = dta)$pvals
  pv %>%
    mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
    left_join(
      env_fixed %>% transmute(Env = as.character(Env), env_effect = predicted.value),
      by = "Env"
    ) %>%
    transmute(Genotype, Env, blup = predicted.value - env_effect, se = std.error)
}

extract_blups_asreml_additive = function(fit, dta, Gmat, Dmat = NULL) {
  pv = predict(fit, classify = "Env:Genotype", data = dta,
                onlyuse = ~fa(Env, k = 2):vm(Genotype, Gmat))$pvals
  env_fixed = predict(fit, classify = "Env", data = dta)$pvals
  pv %>%
    mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
    left_join(
      env_fixed %>% transmute(Env = as.character(Env), env_effect = predicted.value),
      by = "Env"
    ) %>%
    transmute(Genotype, Env, blup = predicted.value - env_effect, se = std.error)
}


extract_blups_asreml_additive_yld = function(fit, dta, Gmat, Dmat = NULL) {
  pv = predict(fit, classify = "Env:Genotype", data = dta,
                onlyuse = ~fa(Env, k = 3):vm(Genotype, Gmat), average = 'factor(Age)')$pvals
  env_fixed = predict(fit, classify = "Env", data = dta, average = 'factor(Age)')$pvals
  pv %>%
    mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
    left_join(
      env_fixed %>% transmute(Env = as.character(Env), env_effect = predicted.value),
      by = "Env"
    ) %>%
    transmute(Genotype, Env, blup = predicted.value - env_effect, se = std.error)
}


# extract_blups_asreml_split = function(fit, dta, ka, ke) {
#   pv_a = predict(fit, classify = "Env:Genotype", data = dta,
#                    onlyuse = ~fa(Env, k = ka):vm(Genotype, Gmat))$pvals
#   pv_e = predict(fit, classify = "Env:Genotype", data = dta,
#                    onlyuse = ~fa(Env, k = ke):ide(Genotype))$pvals
 
#   tidy = function(pv) {
#     pv %>%
#       mutate(Env = as.character(Env), Genotype = as.character(Genotype)) %>%
#       transmute(Genotype, Env, blup = predicted.value, se = std.error)
#   }
#   list(additive = tidy(pv_a), nonadditive = tidy(pv_e))
# }


extract_blups_sommer = function(fit) {

  fa_name = grep("fam\\(", names(fit$uList), value = TRUE)
  # sanity check, fa_name should be 1
  stopifnot(length(fa_name) == 1)
  U = fit$uList[[fa_name]]

  pev_list = fit$uPevList
  pev_name = grep("fam\\(", names(pev_list), value = TRUE)
  stopifnot(length(pev_name) == 1)
  PEV = pev_list[[pev_name]]


  stopifnot(identical(dim(U), dim(PEV)), identical(rownames(U), rownames(PEV)))

  # blups already centered
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




add_reliability = function(blup_long, cov) {
  gv = tibble(Env = rownames(cov), gen_var = diag(cov))
  blup_long %>%
    left_join(gv, by = "Env") %>%
    mutate(pev = se^2, reliability = pmax(0, 1 - pev / gen_var))
}


# 10/21 
#

back_transform_blups = function(blup_long, env_stats, is_sqrt_scale = FALSE, env_var = NULL) {
  out = blup_long %>%
    left_join(env_stats, by = 'Env') %>%
    mutate(
      blup_sqrt = blup * sd_env,  # un-scale: deviation (sqrt-yield units) from the Env's own sqrt-yield mean
      se_sqrt = se * sd_env,
      # kept for backward compatibility
      blup_original = blup_sqrt,
      se_original = se_sqrt
    )

  if (isTRUE(is_sqrt_scale)) {
    out = out %>%
      mutate(
        blup_sqrt_abs = blup_sqrt + mean_env,  
        blup_median = pmax(0, blup_sqrt_abs^2)
      )

    if (!is.null(env_var)) {
      ev = tibble::tibble(Env = names(env_var), env_var = unname(env_var))
      out = out %>%
        left_join(ev, by = 'Env') %>%
        mutate(blup_mean = pmax(0, blup_sqrt_abs^2 + env_var))
    } else {
      out = out %>% mutate(blup_mean = NA_real_)
    }
  }

  out
}


# 10/22 
extract_env_var_sqrt_scale = function(cov, env_stats) {
  v = tibble::tibble(Env = rownames(cov), var_z = diag(cov)) %>%
    left_join(env_stats, by = 'Env') %>%
    mutate(env_var = var_z * sd_env^2)
  setNames(v$env_var, v$Env)
}


# 10/22 
compare_engine_rg = function(corr_asr, corr_sm) {
  stopifnot(setequal(rownames(corr_asr), rownames(corr_sm)))
  envs = rownames(corr_asr)
  corr_sm = corr_sm[envs, envs]
  diff = corr_asr - corr_sm
  off_diag = diff[upper.tri(diff)]
  list(
    max_abs_diff = max(abs(off_diag)),
    mean_abs_diff = mean(abs(off_diag)),
    diff = diff
  )
}

compare_engine_blups = function(blup_asr, blup_sm) {
  joined = blup_asr %>%
    dplyr::select(Genotype, Env, blup_asr = blup) %>%
    dplyr::inner_join(blup_sm %>% dplyr::select(Genotype, Env, blup_sm = blup), by = c('Genotype', 'Env'))
  list(
    n = nrow(joined),
    cor = stats::cor(joined$blup_asr, joined$blup_sm, use = 'complete.obs'),
    rmse = sqrt(mean((joined$blup_asr - joined$blup_sm)^2, na.rm = TRUE)),
    max_abs_diff = max(abs(joined$blup_asr - joined$blup_sm), na.rm = TRUE),
    dta = joined
  )
}


# 10/24 
compare_engine_scalar = function(scalar_asr, scalar_sm) {
  joined = scalar_asr %>%
    dplyr::select(Genotype, g_e_asr = g_e) %>%
    dplyr::inner_join(scalar_sm %>% dplyr::select(Genotype, g_e_sm = g_e), by = 'Genotype')
  list(
    n = nrow(joined),
    cor = stats::cor(joined$g_e_asr, joined$g_e_sm, use = 'complete.obs'),
    rmse = sqrt(mean((joined$g_e_asr - joined$g_e_sm)^2, na.rm = TRUE)),
    max_abs_diff = max(abs(joined$g_e_asr - joined$g_e_sm), na.rm = TRUE),
    dta = joined
  )
}


# 10/23 
flag_boundary_correlations = function(corr, threshold = 0.97) {
  corr_off_diag = corr
  diag(corr_off_diag) = NA
  flagged = which(abs(corr_off_diag) >= threshold, arr.ind = TRUE)
  if (nrow(flagged) == 0) return(tibble::tibble(Env1 = character(0), Env2 = character(0), corr = numeric(0)))
  tibble::tibble(
    Env1 = rownames(corr)[flagged[, 1]],
    Env2 = colnames(corr)[flagged[, 2]],
    corr = corr[flagged]
  ) %>% dplyr::distinct(corr, .keep_all = TRUE)
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
  k = ncol(loadings)
 
  # rotated loadings for all k factors 
  s = svd(loadings)
  Lpc = s$u %*% diag(s$d, nrow = length(s$d))
  if (mean(Lpc[, 1]) < 0) Lpc[, 1] = -Lpc[, 1]
  dimnames(Lpc) = list(colnames(U), paste0('fa', seq_len(k)))
  L1 = Lpc[, 1]
  var_expl = s$d^2 / sum(s$d^2)
 
  tol = 0.05 * max(abs(L1))
  if (any(L1 < -tol)) {
    warning('Not all PC1 loadings are positive (', paste(names(L1)[L1 < -tol], collapse = ', '),').
    Smith & Cullis" OP assumes they are positive or a few very small negatives.')
  }
 

  scores = U %*% Lpc %*% diag(1 / colSums(Lpc^2), nrow = k)
  colnames(scores) = colnames(Lpc)
  f1 = scores[, 1]
 
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
 
  tab = tibble(
    Genotype = geno,
    f1 = f1,
    # OP - value of the line at the mean loading
    OP = mean(L1) * f1,
    n_tested = rowSums(mask),
 
    
    RMSD = if (k > 1) {
      sqrt(rowSums(scores[, 2:k, drop = FALSE]^2) / ncol(U))
    } else {
      rep(0, nrow(scores))
    },
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


# 10/14 add non-additive scalar to blups
#
add_nonadd_scalar_to_stability = function(stab, nonadd_scalar) {
  nonadd_mean = nonadd_scalar %>%
    dplyr::group_by(Genotype) %>%
    dplyr::summarise(g_e = mean(g_e), se_e = sqrt(mean(se_e^2)), .groups = "drop")

  tab = stab$table %>%
    dplyr::rename(OP_additive = OP) %>%
    dplyr::inner_join(nonadd_mean, by = "Genotype") %>%
    dplyr::mutate(OP = OP_additive + g_e)

  op_cut = mean(tab$OP)
  rmsd_cut = stab$rmsd_cut
  tab = tab %>%
    dplyr::mutate(group = dplyr::case_when(
      OP >= op_cut & RMSD <= rmsd_cut ~ 'Above average & stable',
      OP >= op_cut & RMSD > rmsd_cut ~ 'Above average, unstable',
      OP < op_cut & RMSD <= rmsd_cut ~ 'Below average, stable',
      TRUE ~ 'Below average, unstable'
    )) %>%
    dplyr::arrange(desc(OP))

  list(table = tab, pc1_loadings = stab$pc1_loadings, var_explained = stab$var_explained,
       op_cut = op_cut, rmsd_cut = rmsd_cut)
}



# 10/21 
add_nonadd_scalar_to_blups = function(blup_long, nonadd_scalar) {
  blup_long %>%
    dplyr::rename(blup_additive = blup, se_additive = se) %>%
    dplyr::inner_join(nonadd_scalar, by = c("Genotype", "Env")) %>%
    dplyr::mutate(
      blup = blup_additive + g_e,
      se = sqrt(se_additive^2 + se_e^2)
    )
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


# 10/24 
add_pooled_nonadd_to_stability = function(stab, nonadd_scalar) {
  tab = stab$table %>%
    dplyr::rename(OP_additive = OP) %>%
    dplyr::inner_join(nonadd_scalar, by = "Genotype") %>%
    dplyr::mutate(OP = OP_additive + g_e)

  op_cut = mean(tab$OP)
  rmsd_cut = stab$rmsd_cut
  tab = tab %>%
    dplyr::mutate(group = dplyr::case_when(
      OP >= op_cut & RMSD <= rmsd_cut ~ 'Above average & stable',
      OP >= op_cut & RMSD > rmsd_cut ~ 'Above average, unstable',
      OP < op_cut & RMSD <= rmsd_cut ~ 'Below average, stable',
      TRUE ~ 'Below average, unstable'
    )) %>%
    dplyr::arrange(desc(OP))

  list(table = tab, pc1_loadings = stab$pc1_loadings, var_explained = stab$var_explained,
       op_cut = op_cut, rmsd_cut = rmsd_cut)
}


add_pooled_nonadd_to_blups = function(blup_long, nonadd_scalar) {
  blup_long %>%
    dplyr::rename(blup_additive = blup, se_additive = se) %>%
    dplyr::inner_join(nonadd_scalar, by = "Genotype") %>%
    dplyr::mutate(
      blup = blup_additive + g_e,
      se = sqrt(se_additive^2 + se_e^2)
    )
}


# 10/24 
#
flag_rank_shifts = function(stab_additive, stab_total, rank_shift_threshold = 5) {
  add_rank = stab_additive$table %>%
    dplyr::transmute(Genotype, OP_additive = OP, group_additive = group) %>%
    dplyr::mutate(rank_additive = dplyr::min_rank(dplyr::desc(OP_additive)))

  total_rank = stab_total$table %>%
    dplyr::transmute(Genotype, OP_total = OP, group_total = group) %>%
    dplyr::mutate(rank_total = dplyr::min_rank(dplyr::desc(OP_total)))

  dplyr::inner_join(add_rank, total_rank, by = "Genotype") %>%
    dplyr::mutate(
      rank_shift    = rank_additive - rank_total,
      group_changed = group_additive != group_total
    ) %>%
    dplyr::filter(abs(rank_shift) >= rank_shift_threshold | group_changed) %>%
    dplyr::arrange(dplyr::desc(abs(rank_shift)))
}



#' Combine per-trait OP/RMSD into one selection index.
#'
#' @param stab_list named list of `fa_stability()` outputs or the table in it
#' @param weights breeder weights for traits
#'
#' @return tibble: Genotype, n_traits, Index, RMSD_index, and group
build_selection_index = function(stab_list, weights) {

  stopifnot(setequal(names(stab_list), names(weights)))

  # takes raw fa_stability() output or the table in it
  get_table = function(x) if (is.data.frame(x)) x else x$table

  per_trait = purrr::imap_dfr(stab_list, function(x, trait) {
    tab = get_table(x)
    stopifnot(all(c('Genotype', 'OP', 'RMSD') %in% names(tab)))
    tab %>%
      transmute(
        Genotype,
        trait = trait,
        OP_z   = as.numeric(scale(OP)),
        RMSD_z = as.numeric(scale(RMSD))
      )
  })

  per_trait = per_trait %>%
    mutate(
      w = weights[trait],
      OP_contrib   = w * OP_z,
      RMSD_contrib = abs(w) * RMSD_z
    )

  index_tab = per_trait %>%
    group_by(Genotype) %>%
    summarise(
      n_traits   = dplyr::n(),
      Index      = sum(OP_contrib, na.rm = TRUE),
      RMSD_index = sum(RMSD_contrib, na.rm = TRUE),
      .groups = 'drop'
    )

  idx_cut  = mean(index_tab$Index)
  rmsd_cut = median(index_tab$RMSD_index, na.rm = TRUE)

  index_tab %>%
    mutate(group = case_when(
      Index >= idx_cut & RMSD_index <= rmsd_cut ~ 'Above average & stable',
      Index >= idx_cut & RMSD_index >  rmsd_cut ~ 'Above average, unstable',
      Index <  idx_cut & RMSD_index <= rmsd_cut ~ 'Below average, stable',
      TRUE                                      ~ 'Below average, unstable'
    )) %>%
    arrange(desc(Index))
}



