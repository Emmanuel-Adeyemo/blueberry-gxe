# 10/07

library(statgenGxE)


prepare_td = function(blup_long, use_col = 'yield', tested = NULL) {
  d = blup_long %>%
    tidyr::separate(Env, into = c('loc', 'year'), sep = ' ', remove = FALSE) %>%
    mutate(year = as.integer(year))

  if (!is.null(tested)) {
    d = d %>%
      left_join(tested %>% mutate(.tested = TRUE), by = c('Genotype', 'Env')) %>%
      mutate(trait = ifelse(is.na(.tested), NA_real_, .data[[use_col]])) %>%
      dplyr::select(-.tested)
  } else {
    d = d %>% mutate(trait = .data[[use_col]])
  }

  statgenSTA::createTD(
    data = as.data.frame(d),
    genotype = 'Genotype',
    trial = 'Env',
    loc = 'loc',
    year = 'year'
  )
}

# FW
run_fw = function(td, trait = 'yield') {
  fit = gxeFw(td, trait = trait)
  list(
    fit = fit,
    table = fit$estimates,   # genotype, intercept (b0), slope (b1), MSdeviation
    plot_scatter = function() plot(fit, plotType = 'scatter'),
    plot_line = function() plot(fit, plotType = 'line')
  )
}

# AMMI
run_ammi = function(td, trait = 'yield', n_pc = 3) {
  fit = gxeAmmi(td, trait = trait, nPC = n_pc)
  list(
    fit = fit,
    plot_ammi1 = function() plot(fit, plotType = 'AMMI1'),   # mean vs PC1
    plot_ammi2 = function() plot(fit, plotType = 'AMMI2')    # PC1 vs PC2 biplot
  )
}

# GGE
run_gge = function(td, trait = 'yield') {
  fit = gxeGGE(td, trait = trait)
  list(
    fit = fit,
    plot_biplot = function() plot(fit, plotType = 'AMMI2'),
    # plot_which_won = function() plot(fit, plotType = 'AMMI2', primePlot = TRUE)
    
  )
}


run_mega = function(td, trait = 'yield') {
  fit = gxeMegaEnv(td, trait = trait)
  list(
    fit = fit,
    plot_mega = function() plot(fit)
  )
}


run_stability = function(td, trait = 'yield') {
  fit = gxeStability(td, trait = trait)
  list(
    fit = fit,
    summary = summary(fit, pctGeno = 2),
    plot_stab = function() plot(fit)
  )
}


# high cor_f1_pc1 = ammi pc1 agrees with f1/OP ranking - both captures main effects so it should be high
# cor_rmsd_pc2: this is the stability marker. rmsd uses PC1 regression while ammi uses fixed effects residual SVD

cross_check_vs_fa = function(ammi_fit, stab_table) {
  ammi_scores = as.data.frame(ammi_fit$fit$genoScores) %>%
    tibble::rownames_to_column('Genotype')

  merged = stab_table %>%
    dplyr::select(Genotype, f1, OP, RMSD) %>%
    left_join(ammi_scores, by = 'Genotype')

  list(
    table = merged,
    cor_f1_pc1 = suppressWarnings(cor(merged$f1,   merged$PC1, method = 'spearman')),
    cor_rmsd_pc2 = suppressWarnings(cor(merged$RMSD, abs(merged$PC2), method = 'spearman'))
  )
}

