# 10/16

stage07d_targets = list(

  tar_target(
    al_prep_ph,
    prepare_fa_dta(
      aggregate_quality_trait_for_qc(al_plant_height_long_clean, 'plant_height'),
      'plant_height'
    )
  ),

  tar_target(
    prep_ph,
    al_prep_ph %>% filter(Genotype != 'TH-2937')
  ),

  # tar_target(
  #   a_inv_ph,
  #   build_A_inverse(g_mat, unique(prep_ph$Genotype))
  # ),

  tar_target(
    g_mat_ph, g_mat
  ),

  tar_target(
    a_inv_sm_ph,
    prep_ainv_sommer(build_A_inverse(g_mat, unique(prep_ph$Genotype)))
  ),

  tar_target(
    d_mat_ph, d_mat
  ),

  tar_target(
    d_inv_sm_ph,
    prep_ainv_sommer(build_A_inverse(d_mat, unique(prep_ph$Genotype)))
  ),

  # FA sommer
  tar_target(fit_diag_sm_ph, fit_diag_sommer(prep_ph, a_inv_sm_ph)),
  tar_target(fit_fa1_sm_ph, fit_fa_sommer(prep_ph, a_inv_sm_ph, 1)),
  tar_target(fit_fa2_sm_ph, fit_fa_sommer(prep_ph, a_inv_sm_ph, 2)),
  tar_target(fit_fa3_sm_ph, fit_fa_sommer(prep_ph, a_inv_sm_ph, 3)),

  # FA asreml
  tar_target(fit_fa1_asr_ph, fit_fa_asreml(prep_ph, g_mat_ph, ka = 1)),
  tar_target(fit_fa2_asr_ph, fit_fa_asreml(prep_ph, g_mat_ph, ka = 2)),
  tar_target(fit_fa3_asr_ph, fit_fa_asreml(prep_ph, g_mat_ph, ka = 3)),

  # model comparison - sommer
  tar_target(diag_fa1_ph, anova.mmes(fit_diag_sm_ph, fit_fa1_sm_ph)),
  tar_target(fa1_fa2_ph, anova.mmes(fit_fa1_sm_ph, fit_fa2_sm_ph)),
  tar_target(fa2_fa3_ph, anova.mmes(fit_fa2_sm_ph, fit_fa3_sm_ph)),

  tar_target(
    mods_comp_sm_ph,
    dplyr::bind_rows(
      tibble::as_tibble(fa1_fa2_ph[2, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA1'),
      tibble::as_tibble(fa1_fa2_ph[1, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA2'),
      tibble::as_tibble(fa2_fa3_ph[1, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA3')
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),
  tar_target(mods_comp_sm_lrt_ph, run_lrt(mods_comp_sm_ph)),

  # model comp - asreml
  tar_target(mods_asr_ph, list(FA1 = fit_fa1_asr_ph, FA2 = fit_fa2_asr_ph, FA3 = fit_fa3_asr_ph)),
  tar_target(
    mods_comp_asr_ph,
    data.frame(
      Model    = names(mods_asr_ph),
      LogL     = sapply(mods_asr_ph, function(m) m$loglik),
      N_params = sapply(mods_asr_ph, function(m) (m$nwv)),
      AIC      = sapply(mods_asr_ph, function(m) summary(m)$aic),
      BIC      = sapply(mods_asr_ph, function(m) summary(m)$bic)
    )
  ),
  tar_target(mods_comp_asr_lrt_ph, run_lrt(mods_comp_asr_ph)),

  # corr
  tar_target(rg_fa2_sm_ph, extract_rg_matrix_sommer(fit_fa2_sm_ph)),
  tar_target(rg_fa3_sm_ph, extract_rg_matrix_sommer(fit_fa3_sm_ph)),
  tar_target(env_levels_ph, levels(factor(prep_ph$Env))),
  tar_target(rg_fa2_asr_ph, extract_rg_matrix_asreml(summary(fit_fa2_asr_ph)$varcomp, env_levels_ph, 2)),
  tar_target(rg_fa3_asr_ph, extract_rg_matrix_asreml(summary(fit_fa3_asr_ph)$varcomp, env_levels_ph, 3)),

  # hc
  tar_target(hc_fa2_sm_ph, plot_hc(rg_fa2_sm_ph$corr, 'Sommer FA2 Plant Height')),
  tar_target(hc_fa3_sm_ph, plot_hc(rg_fa3_sm_ph$corr, 'Sommer FA3 Plant Height')),
  tar_target(hc_fa2_asr_ph, plot_hc(rg_fa2_asr_ph$corr, 'ASReml FA2 Plant Height')),
  tar_target(hc_fa3_asr_ph, plot_hc(rg_fa3_asr_ph$corr, 'ASReml FA3 Plant Height')),

  tar_target(corr_fa2_sm_ph, plot_corr(rg_fa2_sm_ph$corr, 'Sommer FA2 Plant Height')),
  tar_target(corr_fa3_sm_ph, plot_corr(rg_fa3_sm_ph$corr, 'Sommer FA3 Plant Height')),
  tar_target(corr_fa2_asr_ph, plot_corr(rg_fa2_asr_ph$corr, 'ASReml FA2 Plant Height')),
  tar_target(corr_fa3_asr_ph, plot_corr(rg_fa3_asr_ph$corr, 'ASReml FA3 Plant Height')),

  # FAST - Smith & Cullis 2018
  tar_target(tested_ph, prep_ph %>% dplyr::distinct(Genotype, Env)),

  tar_target(blup_fa2_sm_ph, extract_blups_sommer(fit_fa2_sm_ph)),
  tar_target(blup_fa3_sm_ph, extract_blups_sommer(fit_fa3_sm_ph)),
  tar_target(blup_fa2_asr_ph, extract_blups_asreml(fit_fa2_asr_ph, prep_and_scale_dta(prep_ph), g_mat_ph)),
  tar_target(blup_fa3_asr_ph, extract_blups_asreml(fit_fa3_asr_ph, prep_and_scale_dta(prep_ph), g_mat_ph)),

  tar_target(stab_fa2_sm_ph, fa_stability(blup_fa2_sm_ph, loadings_mmes(fit_fa2_sm_ph)$loadings, tested_ph)),
  tar_target(stab_fa3_sm_ph, fa_stability(blup_fa3_sm_ph, loadings_mmes(fit_fa3_sm_ph)$loadings, tested_ph)),
  tar_target(stab_fa2_asr_ph, fa_stability(blup_fa2_asr_ph, rg_fa2_asr_ph$loadings, tested_ph, engine = 'asreml')),
  tar_target(stab_fa3_asr_ph, fa_stability(blup_fa3_asr_ph, rg_fa3_asr_ph$loadings, tested_ph, engine = 'asreml')),

  tar_target(pc_scores_fa2_sm_ph, fa_pc_scores(blup_fa2_sm_ph, loadings_mmes(fit_fa2_sm_ph)$loadings, tested_ph, n_pc = 2)),
  tar_target(pc_scores_fa3_sm_ph, fa_pc_scores(blup_fa3_sm_ph, loadings_mmes(fit_fa3_sm_ph)$loadings, tested_ph, n_pc = 3)),
  tar_target(pc_scores_fa2_asr_ph, fa_pc_scores(blup_fa2_asr_ph, rg_fa2_asr_ph$loadings, tested_ph, n_pc = 2, engine = 'asreml')),
  tar_target(pc_scores_fa3_asr_ph, fa_pc_scores(blup_fa3_asr_ph, rg_fa3_asr_ph$loadings, tested_ph, n_pc = 3, engine = 'asreml')),

  ## plots
  tar_target(op_rmsd_fa2_sm_ph, plot_op_rmsd(stab_fa2_sm_ph, specs = 'Sommer FA2 Plant Height')),
  tar_target(op_rmsd_fa3_sm_ph, plot_op_rmsd(stab_fa3_sm_ph, specs = 'Sommer FA3 Plant Height')),
  tar_target(op_rmsd_fa2_asr_ph, plot_op_rmsd(stab_fa2_asr_ph, specs = 'ASReml FA2 Plant Height')),
  tar_target(op_rmsd_fa3_asr_ph, plot_op_rmsd(stab_fa3_asr_ph, specs = 'ASReml FA3 Plant Height')),

  tar_target(blup_heatmap_fa2_sm_ph, plot_blup_heatmap(blup_fa2_sm_ph, stab_fa2_sm_ph, tested_ph, specs = 'Sommer FA2 Plant Height')),
  tar_target(blup_heatmap_fa3_sm_ph, plot_blup_heatmap(blup_fa3_sm_ph, stab_fa3_sm_ph, tested_ph, specs = 'Sommer FA3 Plant Height')),
  tar_target(blup_heatmap_fa2_asr_ph, plot_blup_heatmap(blup_fa2_asr_ph, stab_fa2_asr_ph, tested_ph, specs = 'ASReml FA2 Plant Height')),
  tar_target(blup_heatmap_fa3_asr_ph, plot_blup_heatmap(blup_fa3_asr_ph, stab_fa3_asr_ph, tested_ph, specs = 'ASReml FA3 Plant Height')),

  tar_target(rxn_norm_fa2_sm_ph, plot_reaction_norms(blup_fa2_sm_ph, stab_fa2_sm_ph, tested_ph, specs = 'Sommer FA2 Plant Height')),
  tar_target(rxn_norm_fa3_sm_ph, plot_reaction_norms(blup_fa3_sm_ph, stab_fa3_sm_ph, tested_ph, specs = 'Sommer FA3 Plant Height')),
  tar_target(rxn_norm_fa2_asr_ph, plot_reaction_norms(blup_fa2_asr_ph, stab_fa2_asr_ph, tested_ph, specs = 'ASReml FA2 Plant Height')),
  tar_target(rxn_norm_fa3_asr_ph, plot_reaction_norms(blup_fa3_asr_ph, stab_fa3_asr_ph, tested_ph, specs = 'AsReml FA3 Plant Height')),

  tar_target(biplot_fa2_sm_ph, plot_fa_biplot(pc_scores_fa2_sm_ph, specs = 'Sommer FA2 Plant Height')),
  tar_target(biplot_fa3_sm_ph, plot_fa_biplot(pc_scores_fa3_sm_ph, specs = 'Sommer FA3 Plant Height')),
  tar_target(biplot_fa2_asr_ph, plot_fa_biplot(pc_scores_fa2_asr_ph, specs = 'ASReml FA2 Plant Height')),
  tar_target(biplot_fa3_asr_ph, plot_fa_biplot(pc_scores_fa3_asr_ph, specs = 'ASReml FA3 Plant Height')),

  tar_target(uniplot_fa2_sm_ph, plot_uniplot(pc_scores_fa2_sm_ph$loadings, rg_fa2_sm_ph$cov, specs = 'Sommer FA2 Plant Height')),
  tar_target(uniplot_fa3_sm_ph, plot_uniplot(pc_scores_fa3_sm_ph$loadings, rg_fa3_sm_ph$cov, specs = 'Sommer FA3 Plant Height')),
  tar_target(uniplot_fa2_asr_ph, plot_uniplot(pc_scores_fa2_asr_ph$loadings, rg_fa2_asr_ph$cov, specs = 'ASReml FA2 Plant Height')),
  tar_target(uniplot_fa3_asr_ph, plot_uniplot(pc_scores_fa3_asr_ph$loadings, rg_fa3_asr_ph$cov, specs = 'ASReml FA3 Plant Height')),

  # 10/25
  tar_target(
    fit_fa2_min_asr_ph,
    fit_fa_diag_asreml(prep_ph, Gmat = g_mat_ph, Dmat = d_mat_ph, ka = 2)
  ),

  tar_target(
    fit_fa2_ke1_asr_ph,
    fit_fa_ke1_asreml(prep_ph, Gmat = g_mat_ph, Dmat = d_mat_ph, ka = 2)
  ),

  tar_target(
    fit_fa2_min_sm_ph,
    fit_fa3_min_sommer(prep_ph, a_inv_sm_ph, ka = 2)
  ),

  tar_target(
    fit_fa2_ke1_sm_ph,
    fit_fa3_add_nonadd_sommer(prep_ph, a_inv_sm_ph, ka = 2, ke = 1)
  ),

  tar_target(
    mods_noadd_asr_ph,
    list(FA2 = fit_fa2_asr_ph, FA2_min = fit_fa2_min_asr_ph, FA2_ke1 = fit_fa2_ke1_asr_ph)
  ),

  tar_target(
    mods_noadd_comp_asr_ph,
    run_lrt(data.frame(
      Model    = names(mods_noadd_asr_ph),
      LogL     = sapply(mods_noadd_asr_ph, function(m) m$loglik),
      N_params = sapply(mods_noadd_asr_ph, function(m) (m$nwv)),
      AIC      = sapply(mods_noadd_asr_ph, function(m) summary(m)$aic),
      BIC      = sapply(mods_noadd_asr_ph, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa2_fa2min_sm_ph,
    anova.mmes(fit_fa2_sm_ph, fit_fa2_min_sm_ph)
  ),

  tar_target(
    fa2min_fa2ke1_sm_ph,
    anova.mmes(fit_fa2_min_sm_ph, fit_fa2_ke1_sm_ph)
  ),

  tar_target(
    mods_noadd_comp_sm_ph,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa2min_sm_ph[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
      tibble::as_tibble(fa2_fa2min_sm_ph[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_min"),
      tibble::as_tibble(fa2min_fa2ke1_sm_ph[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_ke1")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),

  

  tar_target(
    fit_fa2_min_ide_asr_ph,
    fit_fa_ide_asreml(prep_ph, Gmat = g_mat_ph, Dmat = d_mat_ph, ka = 2)
  ),

  tar_target(
    fit_fa2_min_ide_sm_ph,
    fit_fa3_min_ide_sommer(prep_ph, a_inv_sm_ph, d_inv_sm_ph, ka = 2)
  ),

  tar_target(
    mods_min_ide_asr_ph,
    list(FA2 = fit_fa2_asr_ph, FA2_min_ide = fit_fa2_min_ide_asr_ph)
  ),

  
  tar_target(
    mods_min_ide_comp_asr_ph,
    run_lrt(data.frame(
      Model    = names(mods_min_ide_asr_ph),
      LogL     = sapply(mods_min_ide_asr_ph, function(m) m$loglik),
      N_params = sapply(mods_min_ide_asr_ph, function(m) (m$nwv)),
      AIC      = sapply(mods_min_ide_asr_ph, function(m) summary(m)$aic),
      BIC      = sapply(mods_min_ide_asr_ph, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa2_fa2minide_sm_ph,
    anova.mmes(fit_fa2_sm_ph, fit_fa2_min_ide_sm_ph)
  ),

  tar_target(
    mods_min_ide_comp_sm_ph,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa2minide_sm_ph[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
      tibble::as_tibble(fa2_fa2minide_sm_ph[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_min_ide")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),

  tar_target(
    blup_fa2_min_ide_asr_ph,
    extract_blups_asreml_additive(fit_fa2_min_ide_asr_ph, prep_and_scale_dta(prep_ph), g_mat_ph, d_mat_ph)
  ),

  # 
  tar_target(
    rg_fa2_min_ide_asr_ph,
    extract_rg_matrix_asreml(summary(fit_fa2_min_ide_asr_ph)$varcomp, env_levels_ph, 2)
  ),

  tar_target(
    blup_fa2_min_ide_sm_ph,
    extract_blups_sommer(fit_fa2_min_ide_sm_ph)
  ),

  tar_target(
    rg_fa2_min_ide_sm_ph,
    extract_rg_matrix_sommer(fit_fa2_min_ide_sm_ph)
  ),

  tar_target(
    stab_fa2_min_ide_asr_ph, fa_stability(blup_fa2_min_ide_asr_ph, rg_fa2_min_ide_asr_ph$loadings, tested_ph, engine = 'asreml')
  ),

  tar_target(
    stab_fa2_min_ide_sm_ph, fa_stability(blup_fa2_min_ide_sm_ph, rg_fa2_min_ide_sm_ph$loadings, tested_ph)
  ),

  tar_target(
    pc_scores_fa2_min_ide_asr_ph, fa_pc_scores(blup_fa2_min_ide_asr_ph, rg_fa2_min_ide_asr_ph$loadings, tested_ph, n_pc = 2, engine = 'asreml')
  ),

  tar_target(
    pc_scores_fa2_min_ide_sm_ph, fa_pc_scores(blup_fa2_min_ide_sm_ph, rg_fa2_min_ide_sm_ph$loadings, tested_ph, n_pc = 2)
  ),

  tar_target(
    nonadd_pooled_fa2_asr_ph,
    extract_nonadd_scalar_pooled_asreml(fit_fa2_min_ide_asr_ph, prep_and_scale_dta(prep_ph), g_mat_ph, d_mat_ph)
  ),

  tar_target(
    nonadd_pooled_fa2_sm_ph,
    extract_nonadd_scalar_pooled_sommer(fit_fa2_min_ide_sm_ph)
  ),

  tar_target(
    stab_fa2_total_ide_asr_ph,
    add_pooled_nonadd_to_stability(stab_fa2_min_ide_asr_ph, nonadd_pooled_fa2_asr_ph)
  ),

  tar_target(
    stab_fa2_total_ide_sm_ph,
    add_pooled_nonadd_to_stability(stab_fa2_min_ide_sm_ph, nonadd_pooled_fa2_sm_ph)
  ),

  tar_target(
    blup_fa2_total_ide_asr_ph,
    add_pooled_nonadd_to_blups(blup_fa2_min_ide_asr_ph, nonadd_pooled_fa2_asr_ph)
  ),

  tar_target(
    blup_fa2_total_ide_sm_ph,
    add_pooled_nonadd_to_blups(blup_fa2_min_ide_sm_ph, nonadd_pooled_fa2_sm_ph)
  ),

  tar_target(
    rank_shifts_fa2_ide_asr_ph,
    flag_rank_shifts(stab_fa2_min_ide_asr_ph, stab_fa2_total_ide_asr_ph)
  ),

  tar_target(
    rank_shifts_fa2_ide_sm_ph,
    flag_rank_shifts(stab_fa2_min_ide_sm_ph, stab_fa2_total_ide_sm_ph)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_rg_ph,
    compare_engine_rg(rg_fa2_min_ide_asr_ph$corr, rg_fa2_min_ide_sm_ph$corr)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_blup_ph,
    compare_engine_blups(blup_fa2_min_ide_asr_ph, blup_fa2_min_ide_sm_ph)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_nonadd_ph,
    compare_engine_scalar(nonadd_pooled_fa2_asr_ph, nonadd_pooled_fa2_sm_ph)
  ),

  tar_target(
    hc_fa2_min_ide_asr_ph, plot_hc(rg_fa2_min_ide_asr_ph$corr, 'ASReml FA2_min_ide (additive) Plant Height')
  ),
  tar_target(
    hc_fa2_min_ide_sm_ph, plot_hc(rg_fa2_min_ide_sm_ph$corr, 'Sommer FA2_min_ide (additive) Plant Height')
  ),

  tar_target(
    op_rmsd_fa2_min_ide_asr_ph, plot_op_rmsd(stab_fa2_total_ide_asr_ph, specs = 'ASReml FA2_min_ide (additive + pooled non-additive) Plant Height')
  ),
  tar_target(
    op_rmsd_fa2_min_ide_sm_ph, plot_op_rmsd(stab_fa2_total_ide_sm_ph, specs = 'Sommer FA2_min_ide (additive + pooled non-additive) Plant Height')
  ),

  tar_target(
    blup_heatmap_fa2_min_ide_asr_ph, plot_blup_heatmap(blup_fa2_total_ide_asr_ph, stab_fa2_total_ide_asr_ph, tested_ph, specs = 'ASReml FA2_min_ide (additive + pooled non-additive) Plant Height')
  ),
  tar_target(
    blup_heatmap_fa2_min_ide_sm_ph, plot_blup_heatmap(blup_fa2_total_ide_sm_ph, stab_fa2_total_ide_sm_ph, tested_ph, specs = 'Sommer FA2_min_ide (additive + pooled non-additive) Plant Height')
  ),

  tar_target(
    rxn_norm_fa2_min_ide_asr_ph, plot_reaction_norms(blup_fa2_min_ide_asr_ph, stab_fa2_min_ide_asr_ph, tested_ph, specs = 'ASReml FA2_min_ide (additive) Plant Height')
  ),
  tar_target(
    rxn_norm_fa2_min_ide_sm_ph, plot_reaction_norms(blup_fa2_min_ide_sm_ph, stab_fa2_min_ide_sm_ph, tested_ph, specs = 'Sommer FA2_min_ide (additive) Plant Height')
  ),

  tar_target(
    biplot_fa2_min_ide_asr_ph, plot_fa_biplot(pc_scores_fa2_min_ide_asr_ph, specs = 'ASReml FA2_min_ide (additive) Plant Height')
  ),
  tar_target(
    biplot_fa2_min_ide_sm_ph, plot_fa_biplot(pc_scores_fa2_min_ide_sm_ph, specs = 'Sommer FA2_min_ide (additive) Plant Height')
  )

)