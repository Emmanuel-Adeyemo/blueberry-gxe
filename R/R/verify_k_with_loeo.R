# 10/07

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(purrr)


ALL_ENVS = c("BW 2023", "BW 2024", "BW 2025", "BW 2026",
              "EV 2023", "EV 2024", "EV 2025", "EV 2026")
TRIO = c("BW 2023", "EV 2023", "EV 2026")
WATCH_GENOTYPES = c("NC5289", "MS2162", "SanJoaquin")

loeo_specs = function(trio = TRIO, all_envs = ALL_ENVS) {
  specs = c(
    list(baseline = all_envs),
    setNames(lapply(trio, function(e) setdiff(all_envs, e)),
              paste0("drop_", gsub(" ", "", trio))),
    list(drop_trio = setdiff(all_envs, trio))
  )
  specs
}


leave_one_out_specs = function(all_envs) {
  c(
    list(baseline = all_envs),
    setNames(
      lapply(all_envs, function(e) setdiff(all_envs, e)),
      paste0("drop_", gsub(" ", "", all_envs))
    )
  )
}


.drop_stale_scaling_cols = function(dta) {
  dta %>% dplyr::select(-dplyr::any_of(c("mean_env", "sd_env", "y_scaled")))
}

fit_loeo_sommer = function(dta, ai, envs, orders = 2:3) {
  dta_envs = dta %>% .drop_stale_scaling_cols() %>% dplyr::filter(Env %in% envs)
  AINV = prep_ainv_sommer(ai)
  # diag_fit = fit_diag_sommer(dta_envs, AINV)

  fa_fits = set_names(orders) %>% map(function(k) {
    f = fit_fa_sommer(dta_envs, AINV, k)
    list(
      fit = f,
      cor = extract_rg_matrix_sommer(f)$corr,
      loadings = extract_rg_matrix_sommer(f)$loadings,   
      blups = extract_blups_sommer(f)
    )
  })
  list(engine = "sommer", envs = envs, fa = fa_fits)
}


fit_loeo_asreml = function(dta, ai, envs, orders = 2:3){
  dta_envs = dta %>%
    .drop_stale_scaling_cols() %>%
    dplyr::filter(Env %in% envs) %>%
    dplyr::mutate(Env = droplevels(factor(Env)))

  env_levels = levels(dta_envs$Env)

  fa_fits = set_names(orders) %>% map(function(k){
    if (k == 2) f = fit_fa2_asreml(dta_envs, ai)
    if (k == 3) f = fit_fa3_asreml(dta_envs, ai)

    vc = summary(f)$varcomp
    rg_info = extract_rg_matrix_asreml(varcomp = vc, env_levels = env_levels, n_factors = k)

    blup_df = extract_blups_asreml(f, prep_and_scale_dta(dta_envs), ai)
    loadings = rg_info$loadings
    rownames(loadings) = env_levels

    # TODO: loadings bug: does not drop years in blup -> DONE
    list(
      fit = f,
      cor = rg_info$corr,
      loadings = loadings,  
      blups = blup_df
    )
  })
  list(engine = "asreml", envs = envs, fa = fa_fits)
}

run_loeo = function(dta, ai, engine = c("asreml", "sommer"), specs = loeo_specs(), orders = 2:3) {
  engine = match.arg(engine)

  fit_fn = if (engine == "asreml") fit_loeo_asreml else fit_loeo_sommer
  imap(specs, function(envs, config_name) {
    message(engine, " :: ", config_name, " (", length(envs), " envs)")
    fit_fn(dta, ai, envs = envs, orders = orders)  
  })
}


trio_pair_for_config = function(config_name, trio = TRIO) {
dropped = gsub("drop_", "", config_name)
dropped = gsub("([A-Z]{2})(\\d{4})", "\\1 \\2", dropped)   # "BW2023" = "BW 2023"
if (!dropped %in% trio) return(NULL)
setdiff(trio, dropped)
}

summarize_trio_pairs = function(loeo_results, trio = TRIO) {
  baseline = loeo_results[['baseline']]

  map_dfr(names(loeo_results), function(cfg) {
    pair = trio_pair_for_config(cfg, trio)
    if (is.null(pair)) return(NULL)
    res = loeo_results[[cfg]]

    map_dfr(names(res$fa), function(ord) {
      cor_loeo = res$fa[[ord]]$cor[pair[1], pair[2]]
      cor_base = baseline$fa[[ord]]$cor[pair[1], pair[2]]
      tibble(config = cfg, dropped_env = setdiff(trio, pair), order = paste0("FA", ord),
             pair = paste(pair, collapse = " vs "),
             cor_baseline_with_all3 = cor_base, cor_with_third_dropped = cor_loeo,
             delta = cor_loeo - cor_base)
    })
  })
}

track_watch_genotypes = function(loeo_results, tested_full, order = "2", genotypes = WATCH_GENOTYPES) {

  map_dfr(names(loeo_results), function(cfg) {
    
    res = loeo_results[[cfg]]
    if (!order %in% names(res$fa)) return(NULL)
    envs_here = res$envs
    tested_here = tested_full %>% filter(Env %in% envs_here)
    
    loadings = res$fa[[order]]$loadings  

    U_envs = unique(res$fa[[order]]$blups$Env)
    #debug for asreml
    # cat("in loadings, not in blups:", setdiff(rownames(loadings), U_envs), "\n")
    # cat("in blups, not in loadings:", setdiff(U_envs, rownames(loadings)), "\n")
    # cat("n loadings rows:", nrow(loadings), " n blup envs:", length(U_envs), "\n")    
    
    stab = fa_stability(res$fa[[order]]$blups, loadings, tested = tested_here)
    stab$table %>% filter(Genotype %in% genotypes) %>%
      mutate(config = cfg, order = paste0("FA", order),
             rank_OP = rank(-stab$table$OP)[match(Genotype, stab$table$Genotype)],
             rank_RMSD = rank(stab$table$RMSD)[match(Genotype, stab$table$Genotype)]) %>%
      dplyr::select(config, order, Genotype, OP, RMSD, n_tested, rank_OP, rank_RMSD)
  })
}


summarize_genotype_shifts = function(loeo_results, tested_full, order = "2") {
  baseline = loeo_results[["baseline"]]
  if (!order %in% names(baseline$fa)) return(NULL)
 
  base_tab = fa_stability(
    baseline$fa[[order]]$blups, baseline$fa[[order]]$loadings,
    tested = tested_full %>% dplyr::filter(Env %in% baseline$envs)
  )$table %>% dplyr::select(Genotype, base_OP = OP, base_RMSD = RMSD)
 
  map_dfr(names(loeo_results), function(cfg) {
    res = loeo_results[[cfg]]
    if (!order %in% names(res$fa)) return(NULL)
    tested_here = tested_full %>% dplyr::filter(Env %in% res$envs)
    stab = fa_stability(res$fa[[order]]$blups, res$fa[[order]]$loadings, tested = tested_here)
    stab$table %>%
      dplyr::select(Genotype, OP, RMSD, n_tested) %>%
      dplyr::left_join(base_tab, by = "Genotype") %>%
      dplyr::mutate(config = cfg, order = paste0("FA", order),
                    delta_OP = OP - base_OP, delta_RMSD = RMSD - base_RMSD) %>%
      dplyr::select(config, order, Genotype, OP, delta_OP, RMSD, delta_RMSD, n_tested)
  }) %>%
    dplyr::arrange(dplyr::desc(abs(delta_OP)))
}
 
track_watch_genotypes = function(loeo_results, tested_full, order = "2", genotypes = WATCH_GENOTYPES) {
 
  map_dfr(names(loeo_results), function(cfg) {
    
    res = loeo_results[[cfg]]
    if (!order %in% names(res$fa)) return(NULL)
    envs_here = res$envs
    tested_here = tested_full %>% filter(Env %in% envs_here)
    
    loadings = res$fa[[order]]$loadings  
 
    U_envs = unique(res$fa[[order]]$blups$Env)
    #debug for asreml
    # cat("in loadings, not in blups:", setdiff(rownames(loadings), U_envs), "\n")
    # cat("in blups, not in loadings:", setdiff(U_envs, rownames(loadings)), "\n")
    # cat("n loadings rows:", nrow(loadings), " n blup envs:", length(U_envs), "\n")    
    
    stab = fa_stability(res$fa[[order]]$blups, loadings, tested = tested_here)
    stab$table %>% filter(Genotype %in% genotypes) %>%
      mutate(config = cfg, order = paste0("FA", order),
             rank_OP = rank(-stab$table$OP)[match(Genotype, stab$table$Genotype)],
             rank_RMSD = rank(stab$table$RMSD)[match(Genotype, stab$table$Genotype)]) %>%
      dplyr::select(config, order, Genotype, OP, RMSD, n_tested, rank_OP, rank_RMSD)
  })
}



# added 10/09
#

# 
#' @param loeo_results output of run_loeo()
#' @param raw_dta the UNSCALED fitting dta (prep_yield / prep_bloom)
#' @param orders which FA orders 
loeo_predictive_ability = function(loeo_results, raw_dta, orders = c('2', '3'), group_envs = NULL) {
  baseline_envs = loeo_results[['baseline']]$envs
  configs = setdiff(names(loeo_results), 'baseline')
 
  purrr::map_dfr(orders, function(order) {
    purrr::map_dfr(configs, function(cfg) {
      res = loeo_results[[cfg]]
      if (!order %in% names(res$fa)) return(NULL)
 
      dropped_env = setdiff(baseline_envs, res$envs)
      if (length(dropped_env) != 1) return(NULL)  # e.g. "drop_trio" - no hold out
      if (!is.null(group_envs) && !(dropped_env %in% group_envs)) return(NULL)
 
      blup_long = res$fa[[order]]$blups
 
      if (is.null(group_envs)) {
        tested_train = raw_dta %>%
          dplyr::distinct(Genotype, Env) %>%
          dplyr::filter(Env %in% res$envs)
        op_tab = fa_stability(blup_long, res$fa[[order]]$loadings,
                                tested = tested_train)$table %>%
          dplyr::select(Genotype, OP)
        n_predictor_envs = length(res$envs)
      } else {
        # grouped: mean blup over only the grp's other members
        predictor_envs = intersect(group_envs, res$envs)
        op_tab = blup_long %>%
          dplyr::filter(Env %in% predictor_envs) %>%
          dplyr::group_by(Genotype) %>%
          dplyr::summarise(OP = mean(blup, na.rm = TRUE), .groups = "drop")
        n_predictor_envs = length(predictor_envs)
      }
 
      obs = raw_dta %>% dplyr::filter(Env == dropped_env)
      y_star = obs %>%
        dplyr::group_by(Genotype) %>%
        dplyr::summarise(y_star = mean(y, na.rm = TRUE), .groups = "drop")
 
      joined = dplyr::inner_join(op_tab, y_star, by = "Genotype")
      n_pairs = nrow(joined)
 
      if (n_pairs < 3 || stats::sd(joined$OP) == 0) {
        return(tibble::tibble(config = cfg, order = paste0("FA", order),
                               dropped_env = dropped_env, n_genotypes = n_pairs,
                               n_predictor_envs = n_predictor_envs,
                               r_p = NA_real_, rmse = NA_real_, nrmse = NA_real_))
      }
 
      r_p = stats::cor(joined$OP, joined$y_star)
 
      # 
      env_mean = mean(obs$y, na.rm = TRUE)
      env_sd = stats::sd(obs$y, na.rm = TRUE)
      op_pred = as.numeric(scale(joined$OP)) * env_sd + env_mean
      rmse = sqrt(mean((joined$y_star - op_pred)^2))
      nrmse = if (env_sd > 0) rmse / env_sd else NA_real_
 
      tibble::tibble(config = cfg, order = paste0("FA", order), dropped_env = dropped_env,
                      n_genotypes = n_pairs, n_predictor_envs = n_predictor_envs,
                      r_p = r_p, rmse = rmse, nrmse = nrmse)
    })
  }) %>%
    dplyr::arrange(order, dplyr::desc(abs(r_p)))
}
 