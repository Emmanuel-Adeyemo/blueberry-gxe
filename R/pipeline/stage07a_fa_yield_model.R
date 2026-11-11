# 10/07

stage07a_targets = list(

  # tar_target(
  #   ped_path,
  #   'data/processed/pedigree/SHB_Pedigree_ordered.csv',
  #   format = 'file'
  # ),


  tar_target(
    al_yield_long_clean_filter,
    al_yield_long_clean %>% filter(!PlantID %in% c('BW_MS2162_311_2', 'BW_Farthing_510_2', 'BW_NC5287_407_2'))
  ),

  tar_target(
    al_prep_yield,
    prepare_fa_dta(
      aggregate_yield_for_qc(al_yield_long_clean) %>%
        dplyr::mutate(TotalYield = sqrt(TotalYield)),
      'TotalYield',
      extra_cols = 'Age'
    )
  ),

  tar_target(
    nc_prep_yield,
    prepare_nf_fa_dta(
      nc_yield_clean %>% dplyr::mutate(y = sqrt(y)),
      'y', block_col = 'subplot', extra_cols = 'Age'
    )
  ),

  # tar_target(
  #   fl_prep_yield,
  #   prepare_nf_fa_dta(
  #     fl_yield_clean %>% dplyr::mutate(y = sqrt(y)),
  #     'y', block_col = 'subplot'
  #   )
  # ),

  tar_target(
    all_prep_yield, 
    rbind(al_prep_yield, nc_prep_yield)
  ),

  tar_target(prep_yield,
    all_prep_yield %>% filter(Genotype != 'TH-2937')
  ),

  tar_target(
    g_mat,
    readRDS('data/processed/markers/G_mat.rds')
  ),

  tar_target(
    d_mat,
    readRDS('data/processed/markers/D_mat.rds')
  ),

  # tar_target(
  #   g_inv,
  #   build_A_inverse(g_mat, unique(prep_yield$Genotype))
  # ),
  # tar_target(
  #   g_inv,
  #   g_mat
  # ),

  # tar_target(
  #   a_mat,
  #   build_pedigree_a_matrix(readr::read_csv(ped_path))
  # ),

  # tar_target(
  #   a_inv,
  #   build_A_inverse(a_mat, unique(prep_yield$Genotype))
  # ),

  # tar_target(
  #   a_inv,
  #   g_inv
  # ),

  tar_target(
    a_inv_sm,
    prep_ainv_sommer(build_A_inverse(g_mat, unique(prep_yield$Genotype)))
  ),

  # tar_target(
  #   d_inv, d_mat
  # ),

  tar_target(
    d_inv_sm,
    prep_ainv_sommer(build_A_inverse(d_mat, unique(prep_yield$Genotype)))
  ),

  # FA sommer
  tar_target(
    fit_diag_sm,
    fit_diag_sommer(prep_yield, a_inv_sm, extra_fixed ='factor(Age)')
  ),

  tar_target(
    fit_fa1_sm,
    fit_fa_sommer(prep_yield, a_inv_sm, 1, extra_fixed ='factor(Age)')
  ),

  tar_target(
    fit_fa2_sm,
    fit_fa_sommer(prep_yield, a_inv_sm, 2, extra_fixed ='factor(Age)')
  ),

  tar_target(
    fit_fa3_sm,
    fit_fa_sommer(prep_yield, a_inv_sm, 3, extra_fixed ='factor(Age)')
  ),

  # FA AsRemlR
  tar_target(
    fit_fa1_asr,
    fit_fa_asreml(prep_yield, g_mat, ka = 1, extra_fixed = 'factor(Age)')
  ),

  tar_target(
    fit_fa2_asr,
    fit_fa_asreml(prep_yield, g_mat, ka = 2, extra_fixed = 'factor(Age)')
  ),

  tar_target(
    fit_fa3_asr,
    fit_fa_asreml(prep_yield, g_mat, ka = 3, extra_fixed = 'factor(Age)')
  ),


  # Model Comp
  # Sommer
  # tar_target(
  #   diag_fa1,
  #   anova.mmes(fit_diag_sm, fit_fa1_sm)
  # ),

  # tar_target(
  #   fa1_fa2,
  #   anova.mmes(fit_fa1_sm, fit_fa2_sm)
  # ),

  # tar_target(
  #   fa2_fa3,
  #   anova.mmes(fit_fa2_sm, fit_fa3_sm)
  # ),

  # 
  # tar_target(
  #   mods_comp_sm,
  #   dplyr::bind_rows(
  #     tibble::as_tibble(fa1_fa2[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA1"),
  #     tibble::as_tibble(fa1_fa2[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
  #     tibble::as_tibble(fa2_fa3[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA3")
  #   ) %>%
  #     dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  # ),

  # AsReml
  tar_target(
    mods_asr,
    list(FA2 = fit_fa2_asr, FA3 = fit_fa3_asr)
  ),

  tar_target(
    mods_comp_asr,
    run_lrt(data.frame(
    Model    = names(mods_asr),
    LogL     = sapply(mods_asr, function(m) m$loglik),
    N_params = sapply(mods_asr, function(m) (m$nwv)),
    AIC      = sapply(mods_asr, function(m) summary(m)$aic),
    BIC      = sapply(mods_asr, function(m) summary(m)$bic)
  ))
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

  # AsReml
  tar_target(
    rg_fa2_asr,
    extract_rg_matrix_asreml(summary(fit_fa2_asr)$varcomp, env_levels, 2)
  ),

  tar_target(
    rg_fa3_asr,
    extract_rg_matrix_asreml(summary(fit_fa3_asr)$varcomp, env_levels, 3)
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

  tar_target(
    fa3_asr_corr_plot, plot_corr(rg_fa3_asr$corr, 'ASReml FA3')
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
    blup_fa2_asr, extract_blups_asreml_yld(fit_fa2_asr, prep_and_scale_dta(prep_yield), g_mat)
  ),

  tar_target(
    blup_fa3_asr, extract_blups_asreml_yld(fit_fa3_asr, prep_and_scale_dta(prep_yield), g_mat)
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

  
  tar_target(
    fit_fa3_diag_asr,
    fit_fa_diag_asreml(prep_yield, Gmat = g_mat, Dmat = d_mat, ka = 3, extra_fixed = 'factor(Age)')
  ),

  tar_target(
    fit_fa3_ke1_asr,
    fit_fa_ke1_asreml(prep_yield, Gmat = g_mat, Dmat = d_mat, ka = 3, extra_fixed = 'factor(Age)')
  ),

  tar_target(
    fit_fa3_min_sm,
    fit_fa3_min_sommer(prep_yield, a_inv_sm, ka = 3, extra_fixed = 'factor(Age)')
  ),

  tar_target(
    fit_fa3_ke1_sm,
    fit_fa3_add_nonadd_sommer(prep_yield, a_inv_sm, ka = 3, ke = 1, extra_fixed = 'factor(Age)')
  ),


  tar_target(
    mods_noadd_asr,
    list(FA3 = fit_fa3_asr, FA3_diag = fit_fa3_diag_asr, FA3_ke1 = fit_fa3_ke1_asr)
  ),

  tar_target(
    mods_noadd_comp_asr,
    run_lrt(data.frame(
      Model    = names(mods_noadd_asr),
      LogL     = sapply(mods_noadd_asr, function(m) m$loglik),
      N_params = sapply(mods_noadd_asr, function(m) (m$nwv)),
      AIC      = sapply(mods_noadd_asr, function(m) summary(m)$aic),
      BIC      = sapply(mods_noadd_asr, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa3_fa3min_sm,
    anova.mmes(fit_fa3_sm, fit_fa3_min_sm)
  ),

  tar_target(
    fa3min_fa3ke1_sm,
    anova.mmes(fit_fa3_min_sm, fit_fa3_ke1_sm)
  ),

  tar_target(
    mods_noadd_comp_sm,
    dplyr::bind_rows(
      tibble::as_tibble(fa3_fa3min_sm[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA3"),
      tibble::as_tibble(fa3_fa3min_sm[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA3_min"),
      tibble::as_tibble(fa3min_fa3ke1_sm[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA3_ke1")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
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

 
  tar_target(
    uniplot_fa2_sm, plot_uniplot(pc_scores_fa2_sm$loadings, rg_fa2_sm$cov, specs = 'Sommer FA2')
  ),
  tar_target(
    uniplot_fa3_sm, plot_uniplot(pc_scores_fa3_sm$loadings, rg_fa3_sm$cov, specs = 'Sommer FA3')
  ),
  tar_target(
    uniplot_fa2_asr, plot_uniplot(pc_scores_fa2_asr$loadings, rg_fa2_asr$cov, specs = 'ASReml FA2')
  ),
  tar_target(
    uniplot_fa3_asr, plot_uniplot(pc_scores_fa3_asr$loadings, rg_fa3_asr$cov, specs = 'ASReml FA3')
  ),

  # Group means
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

  # 10/22 gram backtransform
  tar_target(
    yield_env_stats,
    prep_and_scale_dta(prep_yield) %>% dplyr::distinct(Env, mean_env, sd_env)
  ),

  tar_target(
    blup_fa2_sm_grams,
    back_transform_blups(blup_fa2_sm, yield_env_stats, is_sqrt_scale = TRUE,
                          env_var = extract_env_var_sqrt_scale(rg_fa2_sm$cov, yield_env_stats))
  ),
  tar_target(
    blup_fa3_sm_grams,
    back_transform_blups(blup_fa3_sm, yield_env_stats, is_sqrt_scale = TRUE,
                          env_var = extract_env_var_sqrt_scale(rg_fa3_sm$cov, yield_env_stats))
  ),
  tar_target(
    blup_fa2_asr_grams,
    back_transform_blups(blup_fa2_asr, yield_env_stats, is_sqrt_scale = TRUE,
                          env_var = extract_env_var_sqrt_scale(rg_fa2_asr$cov, yield_env_stats))
  ),
  tar_target(
    blup_fa3_asr_grams,
    back_transform_blups(blup_fa3_asr, yield_env_stats, is_sqrt_scale = TRUE,
                          env_var = extract_env_var_sqrt_scale(rg_fa3_asr$cov, yield_env_stats))
  ),

  

  tar_target(
    engine_agreement_fa3_rg,
    compare_engine_rg(rg_fa3_asr$corr, rg_fa3_sm$corr)
  ),
  tar_target(
    engine_agreement_fa3_blup,
    compare_engine_blups(blup_fa3_asr, blup_fa3_sm)
  ),

 




  tar_target(
    fit_fa3_ide_asr,
    fit_fa_ide_asreml(prep_yield, Gmat = g_mat, Dmat = d_mat, ka = 3, extra_fixed = 'factor(Age)')
  ),

  tar_target(
    fit_fa3_min_ide_sm,
    fit_fa3_min_ide_sommer(prep_yield, a_inv_sm, d_inv_sm, ka = 3)
  ),


  ###

  tar_target(
    mods_all_asr,
    list(FA3_base = fit_fa3_asr, FA3_ide = fit_fa3_ide_asr, fA3_diag = fit_fa3_diag_asr, FA3_ke1 = fit_fa3_ke1_asr)
  ),

  tar_target(
    mods_comp_all_asr,
    run_lrt(data.frame(
    Model    = names(mods_all_asr),
    LogL     = sapply(mods_all_asr, function(m) m$loglik),
    N_params = sapply(mods_all_asr, function(m) (m$nwv)),
    AIC      = sapply(mods_all_asr, function(m) summary(m)$aic),
    BIC      = sapply(mods_all_asr, function(m) summary(m)$bic)
  ))
  ),
  ###
  tar_target(
    mods_ide_asr,
    list(FA3 = fit_fa3_asr, FA3_ide = fit_fa3_ide_asr)
  ),

  
  tar_target(
    mods_min_ide_comp_asr,
    run_lrt(data.frame(
      Model    = names(mods_ide_asr),
      LogL     = sapply(mods_ide_asr, function(m) m$loglik),
      N_params = sapply(mods_ide_asr, function(m) (m$nwv)),
      AIC      = sapply(mods_ide_asr, function(m) summary(m)$aic),
      BIC      = sapply(mods_ide_asr, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa3_fa3minide_sm,
    anova.mmes(fit_fa3_sm, fit_fa3_min_ide_sm)
  ),

  # same Self & Liang (1987) boundary-correction here too
  tar_target(
    mods_min_ide_comp_sm,
    dplyr::bind_rows(
      tibble::as_tibble(fa3_fa3minide_sm[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA3"),
      tibble::as_tibble(fa3_fa3minide_sm[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA3_min_ide")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),

  # additive-only blups/rg from FA3_ide itself 
  tar_target(
    blup_fa3_ide_asr,
    extract_blups_asreml_additive_yld(fit_fa3_ide_asr, prep_and_scale_dta(prep_yield), Gmat = g_mat, Dmat = d_mat)
  ),

  
  tar_target(
    rg_fa3_ide_asr,
    extract_rg_matrix_asreml(summary(fit_fa3_ide_asr)$varcomp, env_levels, 3)
  ),

  
  tar_target(
    blup_fa3_min_ide_sm,
    extract_blups_sommer(fit_fa3_min_ide_sm)
  ),

  tar_target(
    rg_fa3_min_ide_sm,
    extract_rg_matrix_sommer(fit_fa3_min_ide_sm)
  ),

  tar_target(
    stab_fa3_ide_asr, fa_stability(blup_fa3_ide_asr, rg_fa3_ide_asr$loadings, tested, engine = 'asreml')
  ),

  tar_target(
    stab_fa3_min_ide_sm, fa_stability(blup_fa3_min_ide_sm, rg_fa3_min_ide_sm$loadings, tested)
  ),

  tar_target(
    pc_scores_fa3_ide_asr, fa_pc_scores(blup_fa3_ide_asr, rg_fa3_ide_asr$loadings, tested, n_pc = 3, engine = 'asreml')
  ),

  tar_target(
    pc_scores_fa3_min_ide_sm, fa_pc_scores(blup_fa3_min_ide_sm, rg_fa3_min_ide_sm$loadings, tested, n_pc = 3)
  ),

  # pooled non-add scalar (one g_e per Genotype, no Env index)
  tar_target(
    nonadd_pooled_fa3_asr,
    extract_nonadd_scalar_pooled_asreml(fit_fa3_ide_asr, prep_and_scale_dta(prep_yield), g_mat, d_mat)
  ),

  tar_target(
    nonadd_pooled_fa3_sm,
    extract_nonadd_scalar_pooled_sommer(fit_fa3_min_ide_sm)
  ),

  # total - add + pooled non-add
  tar_target(
    stab_fa3_total_ide_asr,
    add_pooled_nonadd_to_stability(stab_fa3_ide_asr, nonadd_pooled_fa3_asr)
  ),

  tar_target(
    stab_fa3_total_ide_sm,
    add_pooled_nonadd_to_stability(stab_fa3_min_ide_sm, nonadd_pooled_fa3_sm)
  ),

  tar_target(
    blup_fa3_total_ide_asr,
    add_pooled_nonadd_to_blups(blup_fa3_ide_asr, nonadd_pooled_fa3_asr)
  ),

  tar_target(
    blup_fa3_total_ide_sm,
    add_pooled_nonadd_to_blups(blup_fa3_min_ide_sm, nonadd_pooled_fa3_sm)
  ),

  # rank-shift flags (flag_rank_shifts() 
  tar_target(
    rank_shifts_fa3_ide_asr,
    flag_rank_shifts(stab_fa3_ide_asr, stab_fa3_total_ide_asr)
  ),

  tar_target(
    rank_shifts_fa3_ide_sm,
    flag_rank_shifts(stab_fa3_min_ide_sm, stab_fa3_total_ide_sm)
  ),

  # grams back transform
  tar_target(
    blup_fa3_total_ide_asr_grams,
    back_transform_blups(blup_fa3_total_ide_asr, yield_env_stats, is_sqrt_scale = TRUE,
                          env_var = extract_env_var_sqrt_scale(rg_fa3_ide_asr$cov, yield_env_stats))
  ),

  tar_target(
    blup_fa3_total_ide_sm_grams,
    back_transform_blups(blup_fa3_total_ide_sm, yield_env_stats, is_sqrt_scale = TRUE,
                          env_var = extract_env_var_sqrt_scale(rg_fa3_min_ide_sm$cov, yield_env_stats))
  ),

  # asreml-vs-sommer agreement 
  tar_target(
    engine_agreement_fa3_ide_rg,
    compare_engine_rg(rg_fa3_ide_asr$corr, rg_fa3_min_ide_sm$corr)
  ),

  tar_target(
    engine_agreement_fa3_ide_blup,
    compare_engine_blups(blup_fa3_ide_asr, blup_fa3_min_ide_sm)
  ),

  tar_target(
    engine_agreement_fa3_ide_nonadd,
    compare_engine_scalar(nonadd_pooled_fa3_asr, nonadd_pooled_fa3_sm)
  ),

  # Plots
  tar_target(
    hc_fa3_ide_asr, plot_hc(rg_fa3_ide_asr$corr, 'ASReml FA3_ide (additive)')
  ),
  tar_target(
    hc_fa3_min_ide_sm, plot_hc(rg_fa3_min_ide_sm$corr, 'Sommer FA3_min_ide (additive)')
  ),

  tar_target(
    op_rmsd_fa3_ide_asr, plot_op_rmsd(stab_fa3_total_ide_asr, specs = 'ASReml FA3_ide (additive + pooled non-additive)')
  ),
  tar_target(
    op_rmsd_fa3_min_ide_sm, plot_op_rmsd(stab_fa3_total_ide_sm, specs = 'Sommer FA3_min_ide (additive + pooled non-additive)')
  ),

  tar_target(
    blup_heatmap_fa3_ide_asr, plot_blup_heatmap(blup_fa3_total_ide_asr, stab_fa3_total_ide_asr, tested, specs = 'ASReml FA3_ide (additive + pooled non-additive)')
  ),
  tar_target(
    blup_heatmap_fa3_min_ide_sm, plot_blup_heatmap(blup_fa3_total_ide_sm, stab_fa3_total_ide_sm, tested, specs = 'Sommer FA3_min_ide (additive + pooled non-additive)')
  ),

  tar_target(
    rxn_norm_fa3_ide_asr, plot_reaction_norms(blup_fa3_ide_asr, stab_fa3_ide_asr, tested, specs = 'ASReml FA3_ide (additive)')
  ),
  tar_target(
    rxn_norm_fa3_min_ide_sm, plot_reaction_norms(blup_fa3_min_ide_sm, stab_fa3_min_ide_sm, tested, specs = 'Sommer FA3_min_ide (additive)')
  ),

  tar_target(
    biplot_fa3_ide_asr, plot_fa_biplot(pc_scores_fa3_ide_asr, specs = 'ASReml FA3_ide (additive)')
  ),
  tar_target(
    biplot_fa3_min_ide_sm, plot_fa_biplot(pc_scores_fa3_min_ide_sm, specs = 'Sommer FA3_min_ide (additive)')
  ),

  tar_target(
    group_means_fa3_ide_asr, plot_group_means(
      blup_fa3_total_ide_asr,
      list(normal = c("BW 2024", "BW 2025", "EV 2024", "EV 2025", "NC 2024", "NC 2025"),
                      frost  = c("BW 2023", "EV 2023", "EV 2026")),
                       tested, specs = 'ASReml FA3_ide (additive + pooled non-additive)')
  ),
  tar_target(
    group_means_fa3_min_ide_sm, plot_group_means(
      blup_fa3_total_ide_sm,
      list(normal = c("BW 2024", "BW 2025", "EV 2024", "EV 2025", "NC 2024", "NC 2025"),
                      frost  = c("BW 2023", "EV 2023", "EV 2026")),
                       tested, specs = 'Sommer FA3_min_ide (additive + pooled non-additive)')
  )

)
