# 10/08


library(targets)
library(tarchetypes)  

bloom_traits <- tibble::tibble(trait = c('bloom5', 'bloom50', 'bloom100', 'bloom_duration'))

stage07b_targets <- tar_map(
  values = bloom_traits,
  names = 'trait',

  tar_target(
    al_prep_bloom,
    prepare_fa_dta(aggregate_al_bloom_for_qc(al_bloom_long_clean), trait)
  ),
  tar_target(
    nc_prep_bloom,
    prepare_nf_fa_dta(nc_bloom_doy, trait, block_col = 'subplot')
  ),
  tar_target(
    prep_bloom,
    rbind(al_prep_bloom, nc_prep_bloom)
  ),

  tar_target(
    a_inv_bloom,
    build_A_inverse(a_mat, unique(prep_bloom$Genotype))
  ),
  tar_target(
    a_inv_sm_bloom,
    prep_ainv_sommer(a_inv_bloom)
  ),

  # FA sommer
  tar_target(fit_fa2_sm_bloom, fit_fa_sommer(prep_bloom, a_inv_sm_bloom, 2)),
  tar_target(fit_fa3_sm_bloom, fit_fa_sommer(prep_bloom, a_inv_sm_bloom, 3)),

  # FA asreml
  tar_target(fit_fa2_asr_bloom, fit_fa2_asreml(prep_bloom, a_inv_bloom)),
  tar_target(fit_fa3_asr_bloom, fit_fa3_asreml(prep_bloom, a_inv_bloom)),

  # mod comp
  tar_target(fa2_fa3_bloom, anova.mmes(fit_fa2_sm_bloom, fit_fa3_sm_bloom)),
  tar_target(
    mods_comp_sm_bloom,
    dplyr::bind_rows(
      tibble::as_tibble(fa2_fa3_bloom[2, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA2'),
      tibble::as_tibble(fa2_fa3_bloom[1, c('Df', 'loLik', 'AIC', 'BIC')]) %>% dplyr::mutate(Model = 'FA3')
    ) %>%
      dplyr::transmute(Model, LogL = loLik, N_params = Df, AIC, BIC)
  ),
  
  # lrt
  tar_target(mods_comp_sm_lrt_bloom, run_lrt(mods_comp_sm_bloom)),

  # asreml
  tar_target(mods_asr_bloom, list(FA2 = fit_fa2_asr_bloom, FA3 = fit_fa3_asr_bloom)),
  tar_target(
    mods_comp_asr_bloom,
    data.frame(
      Model    = names(mods_asr_bloom),
      LogL     = sapply(mods_asr_bloom, function(m) m$loglik),
      N_params = sapply(mods_asr_bloom, function(m) (m$nwv)),
      AIC      = sapply(mods_asr_bloom, function(m) summary(m)$aic),
      BIC      = sapply(mods_asr_bloom, function(m) summary(m)$bic)
    )
  ),
  tar_target(mods_comp_asr_lrt_bloom, run_lrt(mods_comp_asr_bloom)),

  tar_target(tested_bloom, prep_bloom %>% dplyr::distinct(Genotype, Env)),
  tar_target(env_levels_bloom, levels(factor(prep_bloom$Env))),

  # corrs
  tar_target(rg_fa2_sm_bloom, extract_rg_matrix_sommer(fit_fa2_sm_bloom)),
  tar_target(rg_fa3_sm_bloom, extract_rg_matrix_sommer(fit_fa3_sm_bloom)),
  tar_target(rg_fa2_asr_bloom, extract_rg_matrix_asreml(summary(fit_fa2_asr_bloom)$varcomp, env_levels_bloom, 2)),
  tar_target(rg_fa3_asr_bloom, extract_rg_matrix_asreml(summary(fit_fa3_asr_bloom)$varcomp, env_levels_bloom, 3)),

  # hc
  tar_target(hc_fa2_sm_bloom, plot_hc(rg_fa2_sm_bloom$corr, paste('Sommer FA2', trait))),
  tar_target(hc_fa3_sm_bloom, plot_hc(rg_fa3_sm_bloom$corr, paste('Sommer FA3', trait))),
  tar_target(hc_fa2_asr_bloom, plot_hc(rg_fa2_asr_bloom$corr, paste('ASReml FA2', trait))),
  tar_target(hc_fa3_asr_bloom, plot_hc(rg_fa3_asr_bloom$corr, paste('ASReml FA3', trait))),

  # corr plots
  tar_target(corr_fa2_sm_bloom, plot_corr(rg_fa2_sm_bloom$corr, paste('Sommer FA2', trait))),
  tar_target(corr_fa3_sm_bloom, plot_corr(rg_fa3_sm_bloom$corr, paste('Sommer FA3', trait))),
  tar_target(corr_fa2_asr_bloom, plot_corr(rg_fa2_asr_bloom$corr, paste('ASReml FA2', trait))),
  tar_target(corr_fa3_asr_bloom, plot_corr(rg_fa3_asr_bloom$corr, paste('ASReml FA3', trait))),

  # blups
  tar_target(blup_fa2_sm_bloom, extract_blups_sommer(fit_fa2_sm_bloom)),
  tar_target(blup_fa3_sm_bloom, extract_blups_sommer(fit_fa3_sm_bloom)),
  tar_target(blup_fa2_asr_bloom, extract_blups_asreml(fit_fa2_asr_bloom, prep_and_scale_dta(prep_bloom), a_inv_bloom)),
  tar_target(blup_fa3_asr_bloom, extract_blups_asreml(fit_fa3_asr_bloom, prep_and_scale_dta(prep_bloom), a_inv_bloom)),

  # Smith & Cullis FAST stability
  tar_target(stab_fa2_sm_bloom, fa_stability(blup_fa2_sm_bloom, loadings_mmes(fit_fa2_sm_bloom)$loadings, tested_bloom)),
  tar_target(stab_fa3_sm_bloom, fa_stability(blup_fa3_sm_bloom, loadings_mmes(fit_fa3_sm_bloom)$loadings, tested_bloom)),
  tar_target(stab_fa2_asr_bloom, fa_stability(blup_fa2_asr_bloom, rg_fa2_asr_bloom$loadings, tested_bloom, engine = 'asreml')),
  tar_target(stab_fa3_asr_bloom, fa_stability(blup_fa3_asr_bloom, rg_fa3_asr_bloom$loadings, tested_bloom, engine = 'asreml')),

  # pc scores 
  tar_target(pc_scores_fa2_sm_bloom, fa_pc_scores(blup_fa2_sm_bloom, loadings_mmes(fit_fa2_sm_bloom)$loadings, tested_bloom, n_pc = 2)),
  tar_target(pc_scores_fa3_sm_bloom, fa_pc_scores(blup_fa3_sm_bloom, loadings_mmes(fit_fa3_sm_bloom)$loadings, tested_bloom, n_pc = 3)),
  tar_target(pc_scores_fa2_asr_bloom, fa_pc_scores(blup_fa2_asr_bloom, rg_fa2_asr_bloom$loadings, tested_bloom, n_pc = 2, engine = 'asreml')),
  tar_target(pc_scores_fa3_asr_bloom, fa_pc_scores(blup_fa3_asr_bloom, rg_fa3_asr_bloom$loadings, tested_bloom, n_pc = 3, engine = 'asreml')),

  ## plots
  # op vs rmsd
  tar_target(op_rmsd_fa2_sm_bloom, plot_op_rmsd(stab_fa2_sm_bloom, specs = paste('Sommer FA2', trait))),
  tar_target(op_rmsd_fa3_sm_bloom, plot_op_rmsd(stab_fa3_sm_bloom, specs = paste('Sommer FA3', trait))),
  tar_target(op_rmsd_fa2_asr_bloom, plot_op_rmsd(stab_fa2_asr_bloom, specs = paste('ASReml FA2', trait))),
  tar_target(op_rmsd_fa3_asr_bloom, plot_op_rmsd(stab_fa3_asr_bloom, specs = paste('ASReml FA3', trait))),

  # heatmaps
  tar_target(blup_heatmap_fa2_sm_bloom, plot_blup_heatmap(blup_fa2_sm_bloom, stab_fa2_sm_bloom, tested_bloom, specs = paste('Sommer FA2', trait))),
  tar_target(blup_heatmap_fa3_sm_bloom, plot_blup_heatmap(blup_fa3_sm_bloom, stab_fa3_sm_bloom, tested_bloom, specs = paste('Sommer FA3', trait))),
  tar_target(blup_heatmap_fa2_asr_bloom, plot_blup_heatmap(blup_fa2_asr_bloom, stab_fa2_asr_bloom, tested_bloom, specs = paste('ASReml FA', trait))),
  tar_target(blup_heatmap_fa3_asr_bloom, plot_blup_heatmap(blup_fa3_asr_bloom, stab_fa3_asr_bloom, tested_bloom, specs = paste('ASReml FA3', trait))),

  # rxn norms
  tar_target(rxn_norm_fa2_sm_bloom, plot_reaction_norms(blup_fa2_sm_bloom, stab_fa2_sm_bloom, tested_bloom, specs = paste('Sommer FA2', trait))),
  tar_target(rxn_norm_fa3_sm_bloom, plot_reaction_norms(blup_fa3_sm_bloom, stab_fa3_sm_bloom, tested_bloom, specs = paste('Sommer FA3', trait))),
  tar_target(rxn_norm_fa2_asr_bloom, plot_reaction_norms(blup_fa2_asr_bloom, stab_fa2_asr_bloom, tested_bloom, specs = paste('ASReml FA2', trait))),
  tar_target(rxn_norm_fa3_asr_bloom, plot_reaction_norms(blup_fa3_asr_bloom, stab_fa3_asr_bloom, tested_bloom, specs = paste('ASReml FA3', trait))),

  # biplots
  tar_target(biplot_fa2_sm_bloom, plot_fa_biplot(pc_scores_fa2_sm_bloom, specs = paste('Sommer FA2', trait))),
  tar_target(biplot_fa3_sm_bloom, plot_fa_biplot(pc_scores_fa3_sm_bloom, specs = paste('Sommer FA3', trait))),
  tar_target(biplot_fa2_asr_bloom, plot_fa_biplot(pc_scores_fa2_asr_bloom, specs = paste('ASReml FA2', trait))),
  tar_target(biplot_fa3_asr_bloom, plot_fa_biplot(pc_scores_fa3_asr_bloom, specs = paste('ASReml FA3', trait))),

  # uniplots
  tar_target(uniplot_fa2_sm_bloom, plot_uniplot(pc_scores_fa2_sm_bloom$loadings, rg_fa2_sm_bloom$cov, specs = paste('Sommer FA2', trait))),
  tar_target(uniplot_fa3_sm_bloom, plot_uniplot(pc_scores_fa3_sm_bloom$loadings, rg_fa3_sm_bloom$cov, specs = paste('Sommer FA3', trait))),
  tar_target(uniplot_fa2_asr_bloom, plot_uniplot(pc_scores_fa2_asr_bloom$loadings, rg_fa2_asr_bloom$cov, specs = paste('ASReml FA2', trait))),
  tar_target(uniplot_fa3_asr_bloom, plot_uniplot(pc_scores_fa3_asr_bloom$loadings, rg_fa3_asr_bloom$cov, specs = paste('ASReml FA3', trait))),

  # loeo
  tar_target(loeo_env_specs_bloom, leave_one_out_specs(sort(unique(prep_bloom$Env)))),
  tar_target(loeo_sm_bloom, run_loeo(prep_bloom, a_inv_sm_bloom, engine = 'sommer', specs = loeo_env_specs_bloom)),
  tar_target(loeo_asr_bloom, run_loeo(prep_bloom, a_inv_bloom, engine = 'asreml', specs = loeo_env_specs_bloom)),

  tar_target(geno_shift_sm_bloom, summarize_genotype_shifts(loeo_sm_bloom, tested_bloom)),
  tar_target(geno_shift_asr_bloom, summarize_genotype_shifts(loeo_asr_bloom, tested_bloom)),

  # pred ability
  tar_target(pred_ability_sm_bloom, loeo_predictive_ability(loeo_sm_bloom, prep_bloom)),
  tar_target(pred_ability_asr_bloom, loeo_predictive_ability(loeo_asr_bloom, prep_bloom)),

  tar_target(pred_ability_plot_sm_bloom, plot_predictive_ability(pred_ability_sm_bloom, specs = paste('Sommer', trait))),
  tar_target(pred_ability_plot_asr_bloom, plot_predictive_ability(pred_ability_asr_bloom, specs = paste('ASReml', trait)))

)