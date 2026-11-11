# 10/16


library(targets)
library(tarchetypes)

quality_traits = tibble::tibble(trait = c('avg_berry_weight', 'firmness', 'brix', 'ta'))

stage07e_targets = tar_map(
  values = quality_traits,
  names = 'trait',

  tar_target(
    al_quality_source_qual,
    switch(
      trait,
      avg_berry_weight = al_avg_berry_weight_long_clean,
      firmness = al_firmness_long_clean,
      brix = al_brix_long_clean,
      ta = al_ta_long_clean
    )
  ),
  tar_target(
    nc_quality_source_qual,
    switch(
      trait,
      avg_berry_weight = nc_avg_berry_weight_clean,
      firmness = nc_firmness_clean,
      brix = nc_brix_clean,
      ta = nc_ta_clean
    )
  ),

  tar_target(
    al_prep_qual,
    prepare_fa_dta(aggregate_quality_trait_for_qc(al_quality_source_qual, trait), trait)
  ),
  tar_target(
    nc_prep_qual,
    prepare_nf_fa_dta(
      aggregate_nc_quality_for_fa(nc_quality_source_qual, trait, block_col = 'subplot'),
      trait,
      block_col = 'subplot'
    )
  ),
  tar_target(
    prep_qual,
    rbind(al_prep_qual, nc_prep_qual) %>% filter(Genotype != 'TH-2937')
  ),

  # tar_target(
  #   a_inv_qual,
  #   build_A_inverse(g_mat, unique(prep_qual$Genotype))
  # ),
  tar_target(
    g_mat_qual, g_mat
  ),
  tar_target(
    a_inv_sm_qual,
    prep_ainv_sommer(build_A_inverse(g_mat, unique(prep_qual$Genotype)))
  ),

  tar_target(
    d_mat_qual, d_mat
  ),

  tar_target(
    d_inv_sm_qual,
    prep_ainv_sommer(build_A_inverse(d_mat, unique(prep_qual$Genotype)))
  ),

  # FA sommer
  tar_target(fit_fa2_sm_qual, fit_fa_sommer(prep_qual, a_inv_sm_qual, 2)),
  tar_target(fit_fa3_sm_qual, fit_fa_sommer(prep_qual, a_inv_sm_qual, 3)),

  # FA asreml
  tar_target(fit_fa1_asr_qual, fit_fa_asreml(prep_qual, g_mat_qual, ka = 1)),
  tar_target(fit_fa2_asr_qual, fit_fa_asreml(prep_qual, g_mat_qual, ka = 2)),
  tar_target(fit_fa3_asr_qual, fit_fa_asreml(prep_qual, g_mat_qual, ka = 3)),

  # mod comp
  tar_target(fa2_fa3_qual, anova.mmes(fit_fa2_sm_qual, fit_fa3_sm_qual)),
  tar_target(
    mods_comp_sm_qual,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa3_qual[2, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA2'),
      tibble::as_tibble(fa2_fa3_qual[1, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA3')
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),

  # lrt
  tar_target(mods_comp_sm_lrt_qual, run_lrt(mods_comp_sm_qual)),

  # asreml
  tar_target(mods_asr_qual, list(FA1 = fit_fa1_asr_qual, FA2 = fit_fa2_asr_qual, FA3 = fit_fa3_asr_qual)),
  tar_target(
    mods_comp_asr_qual,
    data.frame(
      Model    = names(mods_asr_qual),
      LogL     = sapply(mods_asr_qual, function(m) m$loglik),
      N_params = sapply(mods_asr_qual, function(m) (m$nwv)),
      AIC      = sapply(mods_asr_qual, function(m) summary(m)$aic),
      BIC      = sapply(mods_asr_qual, function(m) summary(m)$bic)
    )
  ),
  tar_target(mods_comp_asr_lrt_qual, run_lrt(mods_comp_asr_qual)),

  tar_target(tested_qual, prep_qual %>% dplyr::distinct(Genotype, Env)),
  tar_target(env_levels_qual, levels(factor(prep_qual$Env))),

  # corrs
  tar_target(rg_fa2_sm_qual, extract_rg_matrix_sommer(fit_fa2_sm_qual)),
  tar_target(rg_fa3_sm_qual, extract_rg_matrix_sommer(fit_fa3_sm_qual)),
  tar_target(rg_fa2_asr_qual, extract_rg_matrix_asreml(summary(fit_fa2_asr_qual)$varcomp, env_levels_qual, 2)),
  tar_target(rg_fa3_asr_qual, extract_rg_matrix_asreml(summary(fit_fa3_asr_qual)$varcomp, env_levels_qual, 3)),

  # hc
  tar_target(hc_fa2_sm_qual, plot_hc(rg_fa2_sm_qual$corr, paste('Sommer FA2', trait))),
  tar_target(hc_fa3_sm_qual, plot_hc(rg_fa3_sm_qual$corr, paste('Sommer FA3', trait))),
  tar_target(hc_fa2_asr_qual, plot_hc(rg_fa2_asr_qual$corr, paste('ASReml FA2', trait))),
  tar_target(hc_fa3_asr_qual, plot_hc(rg_fa3_asr_qual$corr, paste('ASReml FA3', trait))),

  # corr plots
  tar_target(corr_fa2_sm_qual, plot_corr(rg_fa2_sm_qual$corr, paste('Sommer FA2', trait))),
  tar_target(corr_fa3_sm_qual, plot_corr(rg_fa3_sm_qual$corr, paste('Sommer FA3', trait))),
  tar_target(corr_fa2_asr_qual, plot_corr(rg_fa2_asr_qual$corr, paste('ASReml FA2', trait))),
  tar_target(corr_fa3_asr_qual, plot_corr(rg_fa3_asr_qual$corr, paste('ASReml FA3', trait))),

  # blups
  tar_target(blup_fa2_sm_qual, extract_blups_sommer(fit_fa2_sm_qual)),
  tar_target(blup_fa3_sm_qual, extract_blups_sommer(fit_fa3_sm_qual)),
  tar_target(blup_fa2_asr_qual, extract_blups_asreml(fit_fa2_asr_qual, prep_and_scale_dta(prep_qual), g_mat_qual)),
  tar_target(blup_fa3_asr_qual, extract_blups_asreml(fit_fa3_asr_qual, prep_and_scale_dta(prep_qual), g_mat_qual)),

  # Smith & Cullis FAST stability
  tar_target(stab_fa2_sm_qual, fa_stability(blup_fa2_sm_qual, loadings_mmes(fit_fa2_sm_qual)$loadings, tested_qual)),
  tar_target(stab_fa3_sm_qual, fa_stability(blup_fa3_sm_qual, loadings_mmes(fit_fa3_sm_qual)$loadings, tested_qual)),
  tar_target(stab_fa2_asr_qual, fa_stability(blup_fa2_asr_qual, rg_fa2_asr_qual$loadings, tested_qual, engine = 'asreml')),
  tar_target(stab_fa3_asr_qual, fa_stability(blup_fa3_asr_qual, rg_fa3_asr_qual$loadings, tested_qual, engine = 'asreml')),

  # pc scores
  tar_target(pc_scores_fa2_sm_qual, fa_pc_scores(blup_fa2_sm_qual, loadings_mmes(fit_fa2_sm_qual)$loadings, tested_qual, n_pc = 2)),
  tar_target(pc_scores_fa3_sm_qual, fa_pc_scores(blup_fa3_sm_qual, loadings_mmes(fit_fa3_sm_qual)$loadings, tested_qual, n_pc = 3)),
  tar_target(pc_scores_fa2_asr_qual, fa_pc_scores(blup_fa2_asr_qual, rg_fa2_asr_qual$loadings, tested_qual, n_pc = 2, engine = 'asreml')),
  tar_target(pc_scores_fa3_asr_qual, fa_pc_scores(blup_fa3_asr_qual, rg_fa3_asr_qual$loadings, tested_qual, n_pc = 3, engine = 'asreml')),

  ## plots
  # op vs rmsd
  tar_target(op_rmsd_fa2_sm_qual, plot_op_rmsd(stab_fa2_sm_qual, specs = paste('Sommer FA2', trait))),
  tar_target(op_rmsd_fa3_sm_qual, plot_op_rmsd(stab_fa3_sm_qual, specs = paste('Sommer FA3', trait))),
  tar_target(op_rmsd_fa2_asr_qual, plot_op_rmsd(stab_fa2_asr_qual, specs = paste('ASReml FA2', trait))),
  tar_target(op_rmsd_fa3_asr_qual, plot_op_rmsd(stab_fa3_asr_qual, specs = paste('ASReml FA3', trait))),

  # heatmaps
  tar_target(blup_heatmap_fa2_sm_qual, plot_blup_heatmap(blup_fa2_sm_qual, stab_fa2_sm_qual, tested_qual, specs = paste('Sommer FA2', trait))),
  tar_target(blup_heatmap_fa3_sm_qual, plot_blup_heatmap(blup_fa3_sm_qual, stab_fa3_sm_qual, tested_qual, specs = paste('Sommer FA3', trait))),
  tar_target(blup_heatmap_fa2_asr_qual, plot_blup_heatmap(blup_fa2_asr_qual, stab_fa2_asr_qual, tested_qual, specs = paste('ASReml FA2', trait))),
  tar_target(blup_heatmap_fa3_asr_qual, plot_blup_heatmap(blup_fa3_asr_qual, stab_fa3_asr_qual, tested_qual, specs = paste('ASReml FA3', trait))),

  # rxn norms
  tar_target(rxn_norm_fa2_sm_qual, plot_reaction_norms(blup_fa2_sm_qual, stab_fa2_sm_qual, tested_qual, specs = paste('Sommer FA2', trait))),
  tar_target(rxn_norm_fa3_sm_qual, plot_reaction_norms(blup_fa3_sm_qual, stab_fa3_sm_qual, tested_qual, specs = paste('Sommer FA3', trait))),
  tar_target(rxn_norm_fa2_asr_qual, plot_reaction_norms(blup_fa2_asr_qual, stab_fa2_asr_qual, tested_qual, specs = paste('ASReml FA2', trait))),
  tar_target(rxn_norm_fa3_asr_qual, plot_reaction_norms(blup_fa3_asr_qual, stab_fa3_asr_qual, tested_qual, specs = paste('ASReml FA3', trait))),

  # biplots
  tar_target(biplot_fa2_sm_qual, plot_fa_biplot(pc_scores_fa2_sm_qual, specs = paste('Sommer FA2', trait))),
  tar_target(biplot_fa3_sm_qual, plot_fa_biplot(pc_scores_fa3_sm_qual, specs = paste('Sommer FA3', trait))),
  tar_target(biplot_fa2_asr_qual, plot_fa_biplot(pc_scores_fa2_asr_qual, specs = paste('ASReml FA2', trait))),
  tar_target(biplot_fa3_asr_qual, plot_fa_biplot(pc_scores_fa3_asr_qual, specs = paste('ASReml FA3', trait))),

  # uniplots
  tar_target(uniplot_fa2_sm_qual, plot_uniplot(pc_scores_fa2_sm_qual$loadings, rg_fa2_sm_qual$cov, specs = paste('Sommer FA2', trait))),
  tar_target(uniplot_fa3_sm_qual, plot_uniplot(pc_scores_fa3_sm_qual$loadings, rg_fa3_sm_qual$cov, specs = paste('Sommer FA3', trait))),
  tar_target(uniplot_fa2_asr_qual, plot_uniplot(pc_scores_fa2_asr_qual$loadings, rg_fa2_asr_qual$cov, specs = paste('ASReml FA2', trait))),
  tar_target(uniplot_fa3_asr_qual, plot_uniplot(pc_scores_fa3_asr_qual$loadings, rg_fa3_asr_qual$cov, specs = paste('ASReml FA3', trait))),

  
  tar_target(
    fit_fa2_min_asr_qual,
    fit_fa_diag_asreml(prep_qual, Gmat = g_mat_qual, Dmat = d_mat_qual, ka = 2)
  ),

  tar_target(
    fit_fa2_ke1_asr_qual,
    fit_fa_ke1_asreml(prep_qual, Gmat = g_mat_qual, Dmat = d_mat_qual, ka = 2)
  ),

  tar_target(
    fit_fa2_min_sm_qual,
    fit_fa3_min_sommer(prep_qual, a_inv_sm_qual, ka = 2)
  ),

  tar_target(
    fit_fa2_ke1_sm_qual,
    fit_fa3_add_nonadd_sommer(prep_qual, a_inv_sm_qual, ka = 2, ke = 1)
  ),

  tar_target(
    mods_noadd_asr_qual,
    list(FA2 = fit_fa2_asr_qual, FA2_min = fit_fa2_min_asr_qual, FA2_ke1 = fit_fa2_ke1_asr_qual)
  ),

  tar_target(
    mods_noadd_comp_asr_qual,
    run_lrt(data.frame(
      Model    = names(mods_noadd_asr_qual),
      LogL     = sapply(mods_noadd_asr_qual, function(m) m$loglik),
      N_params = sapply(mods_noadd_asr_qual, function(m) (m$nwv)),
      AIC      = sapply(mods_noadd_asr_qual, function(m) summary(m)$aic),
      BIC      = sapply(mods_noadd_asr_qual, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa2_fa2min_sm_qual,
    anova.mmes(fit_fa2_sm_qual, fit_fa2_min_sm_qual)
  ),

  tar_target(
    fa2min_fa2ke1_sm_qual,
    anova.mmes(fit_fa2_min_sm_qual, fit_fa2_ke1_sm_qual)
  ),

  tar_target(
    mods_noadd_comp_sm_qual,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa2min_sm_qual[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
      tibble::as_tibble(fa2_fa2min_sm_qual[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_min"),
      tibble::as_tibble(fa2min_fa2ke1_sm_qual[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_ke1")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),


  

  tar_target(
    fit_fa2_min_ide_asr_qual,
    fit_fa_ide_asreml(prep_qual, Gmat = g_mat_qual, Dmat = d_mat_qual, ka = 2)
  ),

  tar_target(
    fit_fa2_min_ide_sm_qual,
    fit_fa3_min_ide_sommer(prep_qual, a_inv_sm_qual, d_inv_sm_qual, ka = 2)
  ),

  tar_target(
    mods_min_ide_asr_qual,
    list(FA2 = fit_fa2_asr_qual, FA2_min_ide = fit_fa2_min_ide_asr_qual)
  ),


  tar_target(
    mods_min_ide_comp_asr_qual,
    run_lrt(data.frame(
      Model    = names(mods_min_ide_asr_qual),
      LogL     = sapply(mods_min_ide_asr_qual, function(m) m$loglik),
      N_params = sapply(mods_min_ide_asr_qual, function(m) (m$nwv)),
      AIC      = sapply(mods_min_ide_asr_qual, function(m) summary(m)$aic),
      BIC      = sapply(mods_min_ide_asr_qual, function(m) summary(m)$bic)
    ))
  ),

  tar_target(
    fa2_fa2minide_sm_qual,
    anova.mmes(fit_fa2_sm_qual, fit_fa2_min_ide_sm_qual)
  ),

  tar_target(
    mods_min_ide_comp_sm_qual,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa2minide_sm_qual[2, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2"),
      tibble::as_tibble(fa2_fa2minide_sm_qual[1, c("Df", "loLik", "AIC", "BIC")]) %>% dplyr::mutate(Model = "FA2_min_ide")
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),

  tar_target(
    blup_fa2_min_ide_asr_qual,
    extract_blups_asreml_additive(fit_fa2_min_ide_asr_qual, prep_and_scale_dta(prep_qual), g_mat_qual, d_mat_qual)
  ),

  tar_target(
    rg_fa2_min_ide_asr_qual,
    extract_rg_matrix_asreml(summary(fit_fa2_min_ide_asr_qual)$varcomp, env_levels_qual, 2)
  ),

  tar_target(
    blup_fa2_min_ide_sm_qual,
    extract_blups_sommer(fit_fa2_min_ide_sm_qual)
  ),

  tar_target(
    rg_fa2_min_ide_sm_qual,
    extract_rg_matrix_sommer(fit_fa2_min_ide_sm_qual)
  ),

  tar_target(
    stab_fa2_min_ide_asr_qual, fa_stability(blup_fa2_min_ide_asr_qual, rg_fa2_min_ide_asr_qual$loadings, tested_qual, engine = 'asreml')
  ),

  tar_target(
    stab_fa2_min_ide_sm_qual, fa_stability(blup_fa2_min_ide_sm_qual, rg_fa2_min_ide_sm_qual$loadings, tested_qual)
  ),

  tar_target(
    pc_scores_fa2_min_ide_asr_qual, fa_pc_scores(blup_fa2_min_ide_asr_qual, rg_fa2_min_ide_asr_qual$loadings, tested_qual, n_pc = 2, engine = 'asreml')
  ),

  tar_target(
    pc_scores_fa2_min_ide_sm_qual, fa_pc_scores(blup_fa2_min_ide_sm_qual, rg_fa2_min_ide_sm_qual$loadings, tested_qual, n_pc = 2)
  ),

  tar_target(
    nonadd_pooled_fa2_asr_qual,
    extract_nonadd_scalar_pooled_asreml(fit_fa2_min_ide_asr_qual, prep_and_scale_dta(prep_qual), g_mat_qual, d_mat_qual)
  ),

  tar_target(
    nonadd_pooled_fa2_sm_qual,
    extract_nonadd_scalar_pooled_sommer(fit_fa2_min_ide_sm_qual)
  ),

  tar_target(
    stab_fa2_total_ide_asr_qual,
    add_pooled_nonadd_to_stability(stab_fa2_min_ide_asr_qual, nonadd_pooled_fa2_asr_qual)
  ),

  tar_target(
    stab_fa2_total_ide_sm_qual,
    add_pooled_nonadd_to_stability(stab_fa2_min_ide_sm_qual, nonadd_pooled_fa2_sm_qual)
  ),

  tar_target(
    blup_fa2_total_ide_asr_qual,
    add_pooled_nonadd_to_blups(blup_fa2_min_ide_asr_qual, nonadd_pooled_fa2_asr_qual)
  ),

  tar_target(
    blup_fa2_total_ide_sm_qual,
    add_pooled_nonadd_to_blups(blup_fa2_min_ide_sm_qual, nonadd_pooled_fa2_sm_qual)
  ),

  tar_target(
    rank_shifts_fa2_ide_asr_qual,
    flag_rank_shifts(stab_fa2_min_ide_asr_qual, stab_fa2_total_ide_asr_qual)
  ),

  tar_target(
    rank_shifts_fa2_ide_sm_qual,
    flag_rank_shifts(stab_fa2_min_ide_sm_qual, stab_fa2_total_ide_sm_qual)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_rg_qual,
    compare_engine_rg(rg_fa2_min_ide_asr_qual$corr, rg_fa2_min_ide_sm_qual$corr)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_blup_qual,
    compare_engine_blups(blup_fa2_min_ide_asr_qual, blup_fa2_min_ide_sm_qual)
  ),

  tar_target(
    engine_agreement_fa2_min_ide_nonadd_qual,
    compare_engine_scalar(nonadd_pooled_fa2_asr_qual, nonadd_pooled_fa2_sm_qual)
  ),

  tar_target(
    hc_fa2_min_ide_asr_qual, plot_hc(rg_fa2_min_ide_asr_qual$corr, paste('ASReml FA2_min_ide (additive)', trait))
  ),
  tar_target(
    hc_fa2_min_ide_sm_qual, plot_hc(rg_fa2_min_ide_sm_qual$corr, paste('Sommer FA2_min_ide (additive)', trait))
  ),

  tar_target(
    op_rmsd_fa2_min_ide_asr_qual, plot_op_rmsd(stab_fa2_total_ide_asr_qual, specs = paste('ASReml FA2_min_ide (additive + pooled non-additive)', trait))
  ),
  tar_target(
    op_rmsd_fa2_min_ide_sm_qual, plot_op_rmsd(stab_fa2_total_ide_sm_qual, specs = paste('Sommer FA2_min_ide (additive + pooled non-additive)', trait))
  ),

  tar_target(
    blup_heatmap_fa2_min_ide_asr_qual, plot_blup_heatmap(blup_fa2_total_ide_asr_qual, stab_fa2_total_ide_asr_qual, tested_qual, specs = paste('ASReml FA2_min_ide (additive + pooled non-additive)', trait))
  ),
  tar_target(
    blup_heatmap_fa2_min_ide_sm_qual, plot_blup_heatmap(blup_fa2_total_ide_sm_qual, stab_fa2_total_ide_sm_qual, tested_qual, specs = paste('Sommer FA2_min_ide (additive + pooled non-additive)', trait))
  ),

  tar_target(
    rxn_norm_fa2_min_ide_asr_qual, plot_reaction_norms(blup_fa2_min_ide_asr_qual, stab_fa2_min_ide_asr_qual, tested_qual, specs = paste('ASReml FA2_min_ide (additive)', trait))
  ),
  tar_target(
    rxn_norm_fa2_min_ide_sm_qual, plot_reaction_norms(blup_fa2_min_ide_sm_qual, stab_fa2_min_ide_sm_qual, tested_qual, specs = paste('Sommer FA2_min_ide (additive)', trait))
  ),

  tar_target(
    biplot_fa2_min_ide_asr_qual, plot_fa_biplot(pc_scores_fa2_min_ide_asr_qual, specs = paste('ASReml FA2_min_ide (additive)', trait))
  ),
  tar_target(
    biplot_fa2_min_ide_sm_qual, plot_fa_biplot(pc_scores_fa2_min_ide_sm_qual, specs = paste('Sommer FA2_min_ide (additive)', trait))
  )

)

