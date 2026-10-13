# 10/07

stage07_targets = list(

  tar_target(
    ped_path,
    'data/processed/pedigree/SHB_Pedigree_ordered.csv',
    format = 'file'
  ),

  tar_target(
    al_prep_yield,
    prepare_fa_dta(aggregate_yield_for_qc(al_yield_long_clean), 'TotalYield')
  ),

  tar_target(
    nc_prep_yield, 
    prepare_nf_fa_dta(nc_yield_clean, 'y', block_col = 'subplot')
  ),

  tar_target(
    prep_yield, 
    rbind(al_prep_yield, nc_prep_yield)
  ),

  tar_target(
    a_mat,
    build_pedigree_a_matrix(readr::read_csv(ped_path))
  ),

  tar_target(
    a_inv,
    build_A_inverse(a_mat, unique(prep_yield$Genotype))
  ),

  tar_target(
    a_inv_sm,
    prep_ainv_sommer(a_inv)
  ),

  # FA sommer
  tar_target(
    fit_diag_sm,
    fit_diag_sommer(prep_yield, a_inv_sm)
  ),

  tar_target(
    fit_fa1_sm,
    fit_fa_sommer(prep_yield, a_inv_sm, 1)
  ),

  tar_target(
    fit_fa2_sm,
    fit_fa_sommer(prep_yield, a_inv_sm, 2)
  ),

  tar_target(
    fit_fa3_sm,
    fit_fa_sommer(prep_yield, a_inv_sm, 3)
  ),

  # FA AsRemlR
  tar_target(
    fit_fa1_asr,
    fit_fa_asreml(prep_yield, a_inv)
  ),

  tar_target(
    fit_fa2_asr,
    fit_fa2_asreml(prep_yield, a_inv)
  ),

  tar_target(
    fit_fa3_asr,
    fit_fa3_asreml(prep_yield, a_inv)
  ),


  # Model Comp
  # Sommer
  tar_target(
    mods_comp_sm,
    dplyr::bind_rows(
      tibble::as_tibble(fa1_fa2[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA1"),
      tibble::as_tibble(fa1_fa2[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
      tibble::as_tibble(fa2_fa3[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA3")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),
  tar_target(
    diag_fa1,
    anova.mmes(fit_diag_sm, fit_fa1_sm)
  ),

  tar_target(
    fa1_fa2,
    anova.mmes(fit_fa1_sm, fit_fa2_sm)
  ),

  tar_target(
    fa2_fa3,
    anova.mmes(fit_fa2_sm, fit_fa3_sm)
  ),

  # AsReml
  tar_target(
    mods_asr,
    list(FA1 = fit_fa1_asr, FA2 = fit_fa2_asr, FA3 = fit_fa3_asr)
  ),

  tar_target(
    mods_comp_asr,
    data.frame(
    Model    = names(mods_asr),
    LogL     = sapply(mods_asr, function(m) m$loglik),
    N_params = sapply(mods_asr, function(m) (m$nwv)),
    AIC      = sapply(mods_asr, function(m) summary(m)$aic),
    BIC      = sapply(mods_asr, function(m) summary(m)$bic)
  )
  ),

  # Correlations
  # sommer
  tar_target(
    rg_fa2_sm,
    extract_rg_matrix_sommer(fit_fa2_sm)
  ),

  tar_target(
    rg_fa3_sm,
    extract_rg_matrix_sommer(fit_fa3_sm)
  ),

  tar_target(
    env_levels,
    levels(factor(prep_yield$Env))
  ),

  # asreml
  tar_target(
    rg_fa2_asr,
    extract_rg_matrix_asreml(summary(fit_fa2_asr)$varcomp, env_levels, 2)
  ),

  tar_target(
    rg_fa3_asr,
    extract_rg_matrix_asreml(summary(fit_fa3_asr)$varcomp, env_levels, 3)
  ),

  # Corr plots
  tar_target(
    corr_fa2_sm, 
    plot_corr(round(rg_fa2_sm$corr, 2), specs = 'Sommer FA2')
  ),
  tar_target(
    corr_fa3_sm, 
    plot_corr(round(rg_fa3_sm$corr, 2), specs = 'Sommer FA3')
  ),
  tar_target(
    corr_fa2_asr, 
    plot_corr(round(rg_fa2_asr$corr, 2), specs = 'ASReml FA2')
  ),
  tar_target(
    corr_fa3_asr, 
    plot_corr(round(rg_fa3_asr$corr, 2), specs = 'ASReml FA3')
  ),

  # HCs
  tar_target(
    hc_fa2_sm,
    plot_hc(rg_fa2_sm$corr, 'Sommer FA2')
  ),

  tar_target(
    hc_fa3_sm,
    plot_hc(rg_fa3_sm$corr, 'Sommer FA3')
  ),

  tar_target(
    hc_fa2_asr,
    plot_hc(rg_fa2_asr$corr, 'ASReml FA2')
  ),

  tar_target(
    hc_fa3_asr,
    plot_hc(rg_fa3_asr$corr, 'ASReml FA3')
  ),

  # FAST ### From Smith and Cullis 2018
  tar_target(
    tested, prep_yield %>% distinct(Genotype, Env)
  ),

  tar_target(
    blup_fa2_sm, extract_blups_sommer(fit_fa2_sm)
  ),

  tar_target(
    blup_fa3_sm, extract_blups_sommer(fit_fa3_sm)
  ),

  tar_target(
    blup_fa2_asr, extract_blups_asreml(fit_fa2_asr, prep_and_scale_dta(prep_yield), a_inv)
  ),

  tar_target(
    blup_fa3_asr, extract_blups_asreml(fit_fa3_asr, prep_and_scale_dta(prep_yield), a_inv)
  ),

  tar_target(
    stab_fa2_sm, fa_stability(blup_fa2_sm, loadings_mmes(fit_fa2_sm)$loadings, tested)
  ),

  tar_target(
    stab_fa3_sm, fa_stability(blup_fa3_sm, loadings_mmes(fit_fa3_sm)$loadings, tested)
  ),

  tar_target(
    stab_fa2_asr, fa_stability(blup_fa2_asr, rg_fa2_asr$loadings, tested, engine = 'asreml')
  ),

  tar_target(
    stab_fa3_asr, fa_stability(blup_fa3_asr, rg_fa3_asr$loadings, tested, engine = 'asreml')
  ),

  tar_target(
    pc_scores_fa2_sm, fa_pc_scores(blup_fa2_sm, loadings_mmes(fit_fa2_sm)$loadings, tested, n_pc = 2)
  ),

  tar_target(
    pc_scores_fa3_sm, fa_pc_scores(blup_fa3_sm, loadings_mmes(fit_fa3_sm)$loadings, tested, n_pc = 3)
  ),

  tar_target(
    pc_scores_fa2_asr, fa_pc_scores(blup_fa2_asr, rg_fa2_asr$loadings, tested, n_pc = 2, engine = 'asreml')
  ),

  tar_target(
    pc_scores_fa3_asr, fa_pc_scores(blup_fa3_asr, rg_fa3_asr$loadings, tested, n_pc = 3, engine = 'asreml')
  ),

  ## PLOTS
  # OP vs RMSD
  tar_target(
    op_rmsd_fa2_sm, plot_op_rmsd(stab_fa2_sm, specs = 'Sommer FA2')
  ),
  tar_target(
    op_rmsd_fa3_sm, plot_op_rmsd(stab_fa3_sm, specs = 'Sommer FA3')
  ),

  tar_target(
    op_rmsd_fa2_asr, plot_op_rmsd(stab_fa2_asr, specs = 'ASReml FA2')
  ),
  tar_target(
    op_rmsd_fa3_asr, plot_op_rmsd(stab_fa3_asr, specs = 'ASReml FA3')
  ),

  # HEATMAPS
  tar_target(
    blup_heatmap_fa2_sm, plot_blup_heatmap(blup_fa2_sm, stab_fa2_sm, tested, specs = 'Sommer FA2')
  ),
  tar_target(
    blup_heatmap_fa3_sm, plot_blup_heatmap(blup_fa3_sm, stab_fa3_sm, tested, specs = 'Sommer FA3')
  ),

  tar_target(
    blup_heatmap_fa2_asr, plot_blup_heatmap(blup_fa2_asr, stab_fa2_asr, tested, specs = 'ASReml FA2')
  ),
  tar_target(
    blup_heatmap_fa3_asr, plot_blup_heatmap(blup_fa3_asr, stab_fa3_asr, tested, specs = 'ASReml FA3')
  ),

  # SENSITIVITY/RXN NORM
  tar_target(
    rxn_norm_fa2_sm, plot_reaction_norms(blup_fa2_sm, stab_fa2_sm, tested, specs = 'Sommer FA2')
  ),
  tar_target(
    rxn_norm_fa3_sm, plot_reaction_norms(blup_fa3_sm, stab_fa3_sm, tested, specs = 'Sommer FA3')
  ),

  tar_target(
    rxn_norm_fa2_asr, plot_reaction_norms(blup_fa2_asr, stab_fa2_asr, tested, specs = 'ASReml FA2')
  ),
  tar_target(
    rxn_norm_fa3_asr, plot_reaction_norms(blup_fa3_asr, stab_fa3_asr, tested, specs = 'ASReml FA3')
  ),

  # BIPLOTS
  tar_target(
    biplot_fa2_sm, plot_fa_biplot(pc_scores_fa2_sm, specs = 'Sommer FA2')
  ),
  tar_target(
    biplot_fa3_sm, plot_fa_biplot(pc_scores_fa3_sm, specs = 'Sommer FA3')
  ),

  tar_target(
    biplot_fa2_asr, plot_fa_biplot(pc_scores_fa2_asr, specs = 'ASReml FA2')
  ),
  tar_target(
    biplot_fa3_asr, plot_fa_biplot(pc_scores_fa3_asr, specs = 'ASReml FA3')
  ),

  ## uniplots added 10/10
  tar_target(
    uniplot_fa2_sm, plot_uniplot(pc_scores_fa2_sm$loadings, rg_fa2_sm$cov, specs = 'Sommer FA2')
  ),
  tar_target(
    uniplot_fa3_sm, plot_uniplot(pc_scores_fa3_sm$loadings, rg_fa3_sm$cov, specs = 'Sommer FA3')
  ),
  tar_target(
    uniplot_fa2_asr, plot_uniplot(pc_scores_fa2_asr$loadings, rg_fa2_asr$cov, specs = 'AsReml FA2')
  ),
  tar_target(
    uniplot_fa3_asr, plot_uniplot(pc_scores_fa3_asr$loadings, rg_fa3_asr$cov, specs = 'AsReml FA3')
  ),

  # Group Means
  tar_target(
    group_means_fa2_sm, plot_group_means(
      blup_fa2_sm,
      list(normal = c("BW 2024", "BW 2025", "EV 2024", "EV 2025", "NC 2024", "NC 2025"),
                      frost  = c("BW 2023", "EV 2023", "EV 2026")),
                       tested, specs = 'Sommer FA2')
  ),
  tar_target(
    group_means_fa3_sm, plot_group_means(
      blup_fa3_sm,
      list(normal = c("BW 2024", "BW 2025", "EV 2024", "EV 2025", "NC 2024", "NC 2025"),
                      frost  = c("BW 2023", "EV 2023", "EV 2026")),
                       tested, specs = 'Sommer FA3')
  ),

  tar_target(
    group_means_fa2_asr, plot_group_means(
      blup_fa2_asr,
      list(normal = c("BW 2024", "BW 2025", "EV 2024", "EV 2025", "NC 2024", "NC 2025"),
                      frost  = c("BW 2023", "EV 2023", "EV 2026")),
                       tested, specs = 'ASReml FA2')
  ),
  tar_target(
    group_means_fa3_asr, plot_group_means(
      blup_fa3_asr,
      list(normal = c("BW 2024", "BW 2025", "EV 2024", "EV 2025", "NC 2024", "NC 2025"),
                      frost  = c("BW 2023", "EV 2023", "EV 2026")),
                       tested, specs = 'ASReml FA3')
  ),


  # Diagnostics for fa2 & 3
  tar_target(
    loeo_sm, run_loeo(prep_yield, a_inv_sm, engine = 'sommer')
  ),
  tar_target(
    loeo_asr, run_loeo(prep_yield, a_inv, engine = 'asreml')
  ),

  tar_target(
    res_sm, summarize_trio_pairs(loeo_sm)
  ),
  tar_target(
    res_asr, summarize_trio_pairs(loeo_asr)
  ),

  tar_target(
    trio_sm, plot_trio_pair_survival(res_sm, specs = 'Sommer')
  ),
  tar_target(
    trio_asr, plot_trio_pair_survival(res_asr, specs = 'ASReml')
  ),


  # LOEO CV for predictive ability
  tar_target(
    loeo_env_specs, leave_one_out_specs(env_levels)
  ),
  tar_target(
    loeo_full_sm, run_loeo(prep_yield, a_inv_sm, engine = 'sommer', specs = loeo_env_specs)
  ),
  tar_target(
    loeo_full_asr, run_loeo(prep_yield, a_inv, engine = 'asreml', specs = loeo_env_specs)
  ),
 
  tar_target(
    pred_ability_sm, loeo_predictive_ability(loeo_full_sm, prep_yield)
  ),
  tar_target(
    pred_ability_asr, loeo_predictive_ability(loeo_full_asr, prep_yield)
  ),
 
  tar_target(
    pred_ability_plot_sm, plot_predictive_ability(pred_ability_sm, specs = 'Sommer')
  ),
  tar_target(
    pred_ability_plot_asr, plot_predictive_ability(pred_ability_asr, specs = 'AsReml')
  ),

  tar_target(
    frost_envs, TRIO
  ),
  tar_target(
    normal_envs, c("BW 2024", "BW 2025", "EV 2024", "EV 2025", "NC 2024", "NC 2025")
  ),
 
  tar_target(
    pred_ability_frost_sm, loeo_predictive_ability(loeo_full_sm, prep_yield, group_envs = frost_envs)
  ),
  tar_target(
    pred_ability_frost_asr, loeo_predictive_ability(loeo_full_asr, prep_yield, group_envs = frost_envs)
  ),
  tar_target(
    pred_ability_normal_sm, loeo_predictive_ability(loeo_full_sm, prep_yield, group_envs = normal_envs)
  ),
  tar_target(
    pred_ability_normal_asr, loeo_predictive_ability(loeo_full_asr, prep_yield, group_envs = normal_envs)
  ),
 
  tar_target(
    pred_ability_plot_frost_sm,
    plot_predictive_ability(pred_ability_frost_sm, specs = 'Sommer, frost-trio predictor only')
  ),
  tar_target(
    pred_ability_plot_frost_asr,
    plot_predictive_ability(pred_ability_frost_asr, specs = 'ASReml, frost-trio predictor only')
  ),
  tar_target(
    pred_ability_plot_normal_sm,
    plot_predictive_ability(pred_ability_normal_sm, specs = 'Sommer, normal-group predictor only')
  ),
  tar_target(
    pred_ability_plot_normal_asr,
    plot_predictive_ability(pred_ability_normal_asr, specs = 'ASReml, normal-group predictor only')
  )
 
 
)



