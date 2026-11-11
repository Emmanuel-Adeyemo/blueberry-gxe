# 10/10


stage07c_targets = list(

  tar_target(
    al_prep_auc,
    prepare_fa_dta(
      aggregate_bloom_to_subplot(al_bloom_auc, 'auc_corrected_percent'),
      'auc_corrected_percent'
    )
  ),

  tar_target(
    prep_auc,
    al_prep_auc %>% filter(Genotype != 'TH-2937')
  ),

  # tar_target(
  #   a_inv_auc,
  #   build_A_inverse(g_mat, unique(prep_auc$Genotype)) # changed to g_mat 10/14
  # ),

  tar_target(
    g_mat_auc, g_mat
  ),

  tar_target(
    a_inv_sm_auc,
    prep_ainv_sommer(build_A_inverse(g_mat, unique(prep_auc$Genotype)))
  ),

  tar_target(
    d_mat_auc, d_mat
  ),

  tar_target(
    d_inv_sm_auc,
    prep_ainv_sommer(build_A_inverse(d_mat, unique(prep_auc$Genotype)))
  ),

  # FA sommer
  tar_target(fit_diag_sm_auc, fit_diag_sommer(prep_auc, a_inv_sm_auc)),
  tar_target(fit_fa1_sm_auc, fit_fa_sommer(prep_auc, a_inv_sm_auc, 1)),
  tar_target(fit_fa2_sm_auc, fit_fa_sommer(prep_auc, a_inv_sm_auc, 2)),
  tar_target(fit_fa3_sm_auc, fit_fa_sommer(prep_auc, a_inv_sm_auc, 3)),

  # FA asreml
  tar_target(fit_fa1_asr_auc, fit_fa_asreml(prep_auc, g_mat_auc, ka = 1)),
  tar_target(fit_fa2_asr_auc, fit_fa_asreml(prep_auc, g_mat_auc, ka = 2)),
  tar_target(fit_fa3_asr_auc, fit_fa_asreml(prep_auc, g_mat_auc, ka = 3)),

  # model comparison - sommer
  tar_target(diag_fa1_auc, anova.mmes(fit_diag_sm_auc, fit_fa1_sm_auc)),
  tar_target(fa1_fa2_auc, anova.mmes(fit_fa1_sm_auc, fit_fa2_sm_auc)),
  tar_target(fa2_fa3_auc, anova.mmes(fit_fa2_sm_auc, fit_fa3_sm_auc)),

  tar_target(
    mods_comp_sm_auc,
    dplyr::bind_rows(
      tibble::as_tibble(fa1_fa2_auc[2, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA1'),
      tibble::as_tibble(fa1_fa2_auc[1, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA2'),
      tibble::as_tibble(fa2_fa3_auc[1, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA3')
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),
  tar_target(mods_comp_sm_lrt_auc, run_lrt(mods_comp_sm_auc)),

  # model comp - asreml
  tar_target(mods_asr_auc, list(FA1 = fit_fa1_asr_auc, FA2 = fit_fa2_asr_auc, FA3 = fit_fa3_asr_auc)),
  tar_target(
    mods_comp_asr_auc,
    data.frame(
      Model    = names(mods_asr_auc),
      LogL     = sapply(mods_asr_auc, function(m) m$loglik),
      N_params = sapply(mods_asr_auc, function(m) (m$nwv)),
      AIC      = sapply(mods_asr_auc, function(m) summary(m)$aic),
      BIC      = sapply(mods_asr_auc, function(m) summary(m)$bic)
    )
  ),
  tar_target(mods_comp_asr_lrt_auc, run_lrt(mods_comp_asr_auc)),

  # corr
  tar_target(rg_fa2_sm_auc, extract_rg_matrix_sommer(fit_fa2_sm_auc)),
  tar_target(rg_fa3_sm_auc, extract_rg_matrix_sommer(fit_fa3_sm_auc)),
  tar_target(env_levels_auc, levels(factor(prep_auc$Env))),
  tar_target(rg_fa2_asr_auc, extract_rg_matrix_asreml(summary(fit_fa2_asr_auc)$varcomp, env_levels_auc, 2)),
  tar_target(rg_fa3_asr_auc, extract_rg_matrix_asreml(summary(fit_fa3_asr_auc)$varcomp, env_levels_auc, 3)),

  # hc
  tar_target(hc_fa2_sm_auc, plot_hc(rg_fa2_sm_auc$corr, 'Sommer FA2 AUC')),
  tar_target(hc_fa3_sm_auc, plot_hc(rg_fa3_sm_auc$corr, 'Sommer FA3 AUC')),
  tar_target(hc_fa2_asr_auc, plot_hc(rg_fa2_asr_auc$corr, 'ASReml FA2 AUC')),
  tar_target(hc_fa3_asr_auc, plot_hc(rg_fa3_asr_auc$corr, 'ASReml FA3 AUC')),

  tar_target(corr_fa2_sm_auc, plot_corr(rg_fa2_sm_auc$corr, 'Sommer FA2 AUC')),
  tar_target(corr_fa3_sm_auc, plot_corr(rg_fa3_sm_auc$corr, 'Sommer FA3 AUC')),
  tar_target(corr_fa2_asr_auc, plot_corr(rg_fa2_asr_auc$corr, 'ASReml FA2 AUC')),
  tar_target(corr_fa3_asr_auc, plot_corr(rg_fa3_asr_auc$corr, 'ASReml FA3 AUC')),

  # FAST - Smith & Cullis 2018
  tar_target(tested_auc, prep_auc %>% dplyr::distinct(Genotype, Env)),

  tar_target(blup_fa2_sm_auc, extract_blups_sommer(fit_fa2_sm_auc)),
  tar_target(blup_fa3_sm_auc, extract_blups_sommer(fit_fa3_sm_auc)),
  tar_target(blup_fa2_asr_auc, extract_blups_asreml(fit_fa2_asr_auc, prep_and_scale_dta(prep_auc), g_mat_auc)),
  tar_target(blup_fa3_asr_auc, extract_blups_asreml(fit_fa3_asr_auc, prep_and_scale_dta(prep_auc), g_mat_auc)),

  tar_target(stab_fa2_sm_auc, fa_stability(blup_fa2_sm_auc, loadings_mmes(fit_fa2_sm_auc)$loadings, tested_auc)),
  tar_target(stab_fa3_sm_auc, fa_stability(blup_fa3_sm_auc, loadings_mmes(fit_fa3_sm_auc)$loadings, tested_auc)),
  tar_target(stab_fa2_asr_auc, fa_stability(blup_fa2_asr_auc, rg_fa2_asr_auc$loadings, tested_auc, engine = 'asreml')),
  tar_target(stab_fa3_asr_auc, fa_stability(blup_fa3_asr_auc, rg_fa3_asr_auc$loadings, tested_auc, engine = 'asreml')),

  tar_target(pc_scores_fa2_sm_auc, fa_pc_scores(blup_fa2_sm_auc, loadings_mmes(fit_fa2_sm_auc)$loadings, tested_auc, n_pc = 2)),
  tar_target(pc_scores_fa3_sm_auc, fa_pc_scores(blup_fa3_sm_auc, loadings_mmes(fit_fa3_sm_auc)$loadings, tested_auc, n_pc = 3)),
  tar_target(pc_scores_fa2_asr_auc, fa_pc_scores(blup_fa2_asr_auc, rg_fa2_asr_auc$loadings, tested_auc, n_pc = 2, engine = 'asreml')),
  tar_target(pc_scores_fa3_asr_auc, fa_pc_scores(blup_fa3_asr_auc, rg_fa3_asr_auc$loadings, tested_auc, n_pc = 3, engine = 'asreml')),

  ## plots
  tar_target(op_rmsd_fa2_sm_auc, plot_op_rmsd(stab_fa2_sm_auc, specs = 'Sommer FA2 AUC')),
  tar_target(op_rmsd_fa3_sm_auc, plot_op_rmsd(stab_fa3_sm_auc, specs = 'Sommer FA3 AUC')),
  tar_target(op_rmsd_fa2_asr_auc, plot_op_rmsd(stab_fa2_asr_auc, specs = 'ASReml FA2 AUC')),
  tar_target(op_rmsd_fa3_asr_auc, plot_op_rmsd(stab_fa3_asr_auc, specs = 'ASReml FA3 AUC')),

  tar_target(blup_heatmap_fa2_sm_auc, plot_blup_heatmap(blup_fa2_sm_auc, stab_fa2_sm_auc, tested_auc, specs = 'Sommer FA2 AUC')),
  tar_target(blup_heatmap_fa3_sm_auc, plot_blup_heatmap(blup_fa3_sm_auc, stab_fa3_sm_auc, tested_auc, specs = 'Sommer FA3 AUC')),
  tar_target(blup_heatmap_fa2_asr_auc, plot_blup_heatmap(blup_fa2_asr_auc, stab_fa2_asr_auc, tested_auc, specs = 'ASReml FA2 AUC')),
  tar_target(blup_heatmap_fa3_asr_auc, plot_blup_heatmap(blup_fa3_asr_auc, stab_fa3_asr_auc, tested_auc, specs = 'ASReml FA3 AUC')),

  tar_target(rxn_norm_fa2_sm_auc, plot_reaction_norms(blup_fa2_sm_auc, stab_fa2_sm_auc, tested_auc, specs = 'Sommer FA2 AUC')),
  tar_target(rxn_norm_fa3_sm_auc, plot_reaction_norms(blup_fa3_sm_auc, stab_fa3_sm_auc, tested_auc, specs = 'Sommer FA3 AUC')),
  tar_target(rxn_norm_fa2_asr_auc, plot_reaction_norms(blup_fa2_asr_auc, stab_fa2_asr_auc, tested_auc, specs = 'ASReml FA2 AUC')),
  tar_target(rxn_norm_fa3_asr_auc, plot_reaction_norms(blup_fa3_asr_auc, stab_fa3_asr_auc, tested_auc, specs = 'AsReml FA3 AUC')),

  tar_target(biplot_fa2_sm_auc, plot_fa_biplot(pc_scores_fa2_sm_auc, specs = 'Sommer FA2 AUC')),
  tar_target(biplot_fa3_sm_auc, plot_fa_biplot(pc_scores_fa3_sm_auc, specs = 'Sommer FA3 AUC')),
  tar_target(biplot_fa2_asr_auc, plot_fa_biplot(pc_scores_fa2_asr_auc, specs = 'ASReml FA2 AUC')),
  tar_target(biplot_fa3_asr_auc, plot_fa_biplot(pc_scores_fa3_asr_auc, specs = 'ASReml FA3 AUC')),

  tar_target(uniplot_fa2_sm_auc, plot_uniplot(pc_scores_fa2_sm_auc$loadings, rg_fa2_sm_auc$cov, specs = 'Sommer FA2 AUC')),
  tar_target(uniplot_fa3_sm_auc, plot_uniplot(pc_scores_fa3_sm_auc$loadings, rg_fa3_sm_auc$cov, specs = 'Sommer FA3 AUC')),
  tar_target(uniplot_fa2_asr_auc, plot_uniplot(pc_scores_fa2_asr_auc$loadings, rg_fa2_asr_auc$cov, specs = 'ASReml FA2 AUC')),
  tar_target(uniplot_fa3_asr_auc, plot_uniplot(pc_scores_fa3_asr_auc$loadings, rg_fa3_asr_auc$cov, specs = 'ASReml FA3 AUC')),

  # 10/25 
  tar_target(
    fit_fa2_min_asr_auc,
    fit_fa_diag_asreml(prep_auc, Gmat = g_mat_auc, Dmat = d_mat_auc, ka = 2)
  ),

  tar_target(
    fit_fa2_ke1_asr_auc,
    fit_fa_ke1_asreml(prep_auc, Gmat = g_mat_auc, Dmat = d_mat_auc, ka = 2)
  ),

  tar_target(
    fit_fa2_min_sm_auc,
    fit_fa3_min_sommer(prep_auc, a_inv_sm_auc, ka = 2)
  ),

  tar_target(
    fit_fa2_ke1_sm_auc,
    fit_fa3_add_nonadd_sommer(prep_auc, a_inv_sm_auc, ka = 2, ke = 1)
  ),

  tar_target(
    mods_noadd_asr_auc,
    list(FA2 = fit_fa2_asr_auc, FA2_min = fit_fa2_min_asr_auc, FA2_ke1 = fit_fa2_ke1_asr_auc)
  ),

  tar_target(
    mods_noadd_comp_asr_auc,
    run_lrt(data.frame(
      Model    = names(mods_noadd_asr_auc),
      LogL     = sapply(mods_noadd_asr_auc, function(m) m$loglik),
      N_params = sapply(mods_noadd_asr_auc, function(m) (m$nwv)),
      AIC      = sapply(mods_noadd_asr_auc, function(m) summary(m)$aic),
      BIC      = sapply(mods_noadd_asr_auc, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa2_fa2min_sm_auc,
    anova.mmes(fit_fa2_sm_auc, fit_fa2_min_sm_auc)
  ),

  tar_target(
    fa2min_fa2ke1_sm_auc,
    anova.mmes(fit_fa2_min_sm_auc, fit_fa2_ke1_sm_auc)
  ),

  tar_target(
    mods_noadd_comp_sm_auc,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa2min_sm_auc[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
      tibble::as_tibble(fa2_fa2min_sm_auc[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_min"),
      tibble::as_tibble(fa2min_fa2ke1_sm_auc[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_ke1")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),




  tar_target(
    fit_fa2_min_ide_asr_auc,
    fit_fa_ide_asreml(prep_auc, Gmat = g_mat_auc, Dmat = d_mat_auc, ka = 2)
  ),

  tar_target(
    fit_fa2_min_ide_sm_auc,
    fit_fa3_min_ide_sommer(prep_auc, a_inv_sm_auc, d_inv_sm_auc, ka = 2)
  ),

  tar_target(
    mods_min_ide_asr_auc,
    list(FA2 = fit_fa2_asr_auc, FA2_min_ide = fit_fa2_min_ide_asr_auc)
  ),


  # per Self & Liang (1987), halve chi-sq_1 LRT p-value 
  tar_target(
    mods_min_ide_comp_asr_auc,
    run_lrt(data.frame(
      Model    = names(mods_min_ide_asr_auc),
      LogL     = sapply(mods_min_ide_asr_auc, function(m) m$loglik),
      N_params = sapply(mods_min_ide_asr_auc, function(m) (m$nwv)),
      AIC      = sapply(mods_min_ide_asr_auc, function(m) summary(m)$aic),
      BIC      = sapply(mods_min_ide_asr_auc, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa2_fa2minide_sm_auc,
    anova.mmes(fit_fa2_sm_auc, fit_fa2_min_ide_sm_auc)
  ),

  tar_target(
    mods_min_ide_comp_sm_auc,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa2minide_sm_auc[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
      tibble::as_tibble(fa2_fa2minide_sm_auc[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_min_ide")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),

  tar_target(
    blup_fa2_min_ide_asr_auc,
    extract_blups_asreml_additive(fit_fa2_min_ide_asr_auc, prep_and_scale_dta(prep_auc), g_mat_auc, d_mat_auc)
  ),


  tar_target(
    rg_fa2_min_ide_asr_auc,
    extract_rg_matrix_asreml(summary(fit_fa2_min_ide_asr_auc)$varcomp, env_levels_auc, 2)
  ),

  tar_target(
    blup_fa2_min_ide_sm_auc,
    extract_blups_sommer(fit_fa2_min_ide_sm_auc)
  ),

  tar_target(
    rg_fa2_min_ide_sm_auc,
    extract_rg_matrix_sommer(fit_fa2_min_ide_sm_auc)
  ),

  tar_target(
    stab_fa2_min_ide_asr_auc, fa_stability(blup_fa2_min_ide_asr_auc, rg_fa2_min_ide_asr_auc$loadings, tested_auc, engine = 'asreml')
  ),

  tar_target(
    stab_fa2_min_ide_sm_auc, fa_stability(blup_fa2_min_ide_sm_auc, rg_fa2_min_ide_sm_auc$loadings, tested_auc)
  ),

  tar_target(
    pc_scores_fa2_min_ide_asr_auc, fa_pc_scores(blup_fa2_min_ide_asr_auc, rg_fa2_min_ide_asr_auc$loadings, tested_auc, n_pc = 2, engine = 'asreml')
  ),

  tar_target(
    pc_scores_fa2_min_ide_sm_auc, fa_pc_scores(blup_fa2_min_ide_sm_auc, rg_fa2_min_ide_sm_auc$loadings, tested_auc, n_pc = 2)
  ),

  tar_target(
    nonadd_pooled_fa2_asr_auc,
    extract_nonadd_scalar_pooled_asreml(fit_fa2_min_ide_asr_auc, prep_and_scale_dta(prep_auc), g_mat_auc, d_mat_auc)
  ),

  tar_target(
    nonadd_pooled_fa2_sm_auc,
    extract_nonadd_scalar_pooled_sommer(fit_fa2_min_ide_sm_auc)
  ),

  tar_target(
    stab_fa2_total_ide_asr_auc,
    add_pooled_nonadd_to_stability(stab_fa2_min_ide_asr_auc, nonadd_pooled_fa2_asr_auc)
  ),

  tar_target(
    stab_fa2_total_ide_sm_auc,
    add_pooled_nonadd_to_stability(stab_fa2_min_ide_sm_auc, nonadd_pooled_fa2_sm_auc)
  ),

  tar_target(
    blup_fa2_total_ide_asr_auc,
    add_pooled_nonadd_to_blups(blup_fa2_min_ide_asr_auc, nonadd_pooled_fa2_asr_auc)
  ),

  tar_target(
    blup_fa2_total_ide_sm_auc,
    add_pooled_nonadd_to_blups(blup_fa2_min_ide_sm_auc, nonadd_pooled_fa2_sm_auc)
  ),

  tar_target(
    rank_shifts_fa2_ide_asr_auc,
    flag_rank_shifts(stab_fa2_min_ide_asr_auc, stab_fa2_total_ide_asr_auc)
  ),

  tar_target(
    rank_shifts_fa2_ide_sm_auc,
    flag_rank_shifts(stab_fa2_min_ide_sm_auc, stab_fa2_total_ide_sm_auc)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_rg_auc,
    compare_engine_rg(rg_fa2_min_ide_asr_auc$corr, rg_fa2_min_ide_sm_auc$corr)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_blup_auc,
    compare_engine_blups(blup_fa2_min_ide_asr_auc, blup_fa2_min_ide_sm_auc)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_nonadd_auc,
    compare_engine_scalar(nonadd_pooled_fa2_asr_auc, nonadd_pooled_fa2_sm_auc)
  ),

  tar_target(
    hc_fa2_min_ide_asr_auc, plot_hc(rg_fa2_min_ide_asr_auc$corr, 'ASReml FA2_min_ide (additive) AUC')
  ),
  tar_target(
    hc_fa2_min_ide_sm_auc, plot_hc(rg_fa2_min_ide_sm_auc$corr, 'Sommer FA2_min_ide (additive) AUC')
  ),

  tar_target(
    op_rmsd_fa2_min_ide_asr_auc, plot_op_rmsd(stab_fa2_total_ide_asr_auc, specs = 'ASReml FA2_min_ide (additive + pooled non-additive) AUC')
  ),
  tar_target(
    op_rmsd_fa2_min_ide_sm_auc, plot_op_rmsd(stab_fa2_total_ide_sm_auc, specs = 'Sommer FA2_min_ide (additive + pooled non-additive) AUC')
  ),

  tar_target(
    blup_heatmap_fa2_min_ide_asr_auc, plot_blup_heatmap(blup_fa2_total_ide_asr_auc, stab_fa2_total_ide_asr_auc, tested_auc, specs = 'ASReml FA2_min_ide (additive + pooled non-additive) AUC')
  ),
  tar_target(
    blup_heatmap_fa2_min_ide_sm_auc, plot_blup_heatmap(blup_fa2_total_ide_sm_auc, stab_fa2_total_ide_sm_auc, tested_auc, specs = 'Sommer FA2_min_ide (additive + pooled non-additive) AUC')
  ),

  tar_target(
    rxn_norm_fa2_min_ide_asr_auc, plot_reaction_norms(blup_fa2_min_ide_asr_auc, stab_fa2_min_ide_asr_auc, tested_auc, specs = 'ASReml FA2_min_ide (additive) AUC')
  ),
  tar_target(
    rxn_norm_fa2_min_ide_sm_auc, plot_reaction_norms(blup_fa2_min_ide_sm_auc, stab_fa2_min_ide_sm_auc, tested_auc, specs = 'Sommer FA2_min_ide (additive) AUC')
  ),

  tar_target(
    biplot_fa2_min_ide_asr_auc, plot_fa_biplot(pc_scores_fa2_min_ide_asr_auc, specs = 'ASReml FA2_min_ide (additive) AUC')
  ),
  tar_target(
    biplot_fa2_min_ide_sm_auc, plot_fa_biplot(pc_scores_fa2_min_ide_sm_auc, specs = 'Sommer FA2_min_ide (additive) AUC')
  )

)


