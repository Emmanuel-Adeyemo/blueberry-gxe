# 10/07

library(ggplot2)
library(dplyr)
library(glue)
library(factoextra)
library(ggcorrplot)


group_cols = c(
  'Above average & stable' = '#1b7837',
  'Above average, unstable' = '#e08214',
  'Below average, stable' = '#4393c3',
  'Below average, unstable' = '#b2182b'
)


plot_corr = function(corr_mat, specs){
  ggcorrplot(corr_mat, 
           hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method = 'square', 
           colors = c('#D95F02', 'white', '#2B5C8F'), 
           title = glue('Genetic Correlation Matrix ({specs})'), 
           ggtheme = theme_minimal())
}


plot_hc = function(corr_mat, specs){
  dist_mat = as.dist(1 - corr_mat)
  hc = hclust(dist_mat, method = 'average')
  min_y = -0.25 * max(hc$height)
  
  fviz_dend(
    hc, 
    main = glue::glue('Envt clustering ({specs})'),
    ggtheme = theme_minimal()
  ) +
    scale_y_continuous(limits = c(min_y, NA)) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = margin(t = 10, r = 10, b = 40, l = 10, unit = 'pt'))

}



# upper left: high performer and stable
plot_op_rmsd = function(stab, label_n = 30, specs) {

  t = stab$table
  grey_lab = grep('^tested in', unique(t$group), value = TRUE)
  cols = c(group_cols, setNames(rep('grey60', length(grey_lab)), grey_lab))
  lab = t %>% slice_max(OP, n = label_n) %>%
    bind_rows(t %>% filter(group == "Above average & stable")) %>% distinct(Genotype, .keep_all = TRUE)

  p = ggplot(t, aes(RMSD, OP, colour = group, size = n_tested)) +
    geom_hline(yintercept = stab$op_cut, linetype = 2, colour = "grey50") +
    geom_vline(xintercept = stab$rmsd_cut, linetype = 2, colour = "grey50") +
    geom_point(alpha = 0.85) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_size_continuous(range = c(1.5, 5), name = "environments\ntested") +
    labs(x = "RMSD (lower = more stable)",
         y = "Overall performance in SD units (higher = better)",
         title = glue("Overall performance vs. stability ({specs})"),
         subtitle = "Dashed lines = (y: mean of genotypes, x: median of genotypes)") +
    theme_minimal()

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_text_repel(data = lab, aes(label = Genotype), size = 3,
                                 show.legend = FALSE, max.overlaps = 30)
  } else {
    p + geom_text(data = lab, aes(label = Genotype), size = 3, vjust = -1,
                  check_overlap = TRUE, show.legend = FALSE)
  }
}

# gxe blups heatmap
plot_blup_heatmap = function(blup_long, stab, tested, specs) {
  env_order = names(sort(stab$pc1_loadings, decreasing = TRUE))
  
  d = blup_long %>%
    left_join(tested %>% mutate(tested = TRUE), by = c("Genotype", "Env")) %>%
    mutate(tested = !is.na(tested),
           Genotype = factor(Genotype, levels = rev(stab$table$Genotype)),  # best OP at top
           Env = factor(Env, levels = env_order))
  
  ggplot(d, aes(Env, Genotype, fill = blup)) +
    geom_tile() +
    geom_point(data = filter(d, tested), size = 0.6, colour = "black") +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0,
                         name = "BLUP\n(SD units)") +
    labs(x = NULL, y = NULL,
         title = glue("Genotype BLUPs by environment ({specs})"),
         subtitle = "Rows ordered by overall performance; columns by PC1 loading; dot = actually tested") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# sensitivity
plot_reaction_norms = function(blup_long, stab, tested, specs, n_highlight = 6) {
  env_order = names(sort(stab$pc1_loadings, decreasing = TRUE))
  top = stab$table %>% filter(group == "Above average & stable") %>%
    slice_max(OP, n = n_highlight) %>% pull(Genotype)

  d = blup_long %>%
    inner_join(tested, by = c("Genotype", "Env")) %>%      # tested cells only
    mutate(Env = factor(Env, levels = env_order),
           hl = ifelse(Genotype %in% top, Genotype, NA))
  
  ggplot(d, aes(Env, blup, group = Genotype)) +
    geom_line(data = filter(d, is.na(hl)), colour = "grey90") +
    geom_line(data = filter(d, !is.na(hl)), aes(colour = hl), linewidth = 1) +
    geom_point(data = filter(d, !is.na(hl)), aes(colour = hl)) +
    labs(x = "Environment (ordered by PC1 loading)", y = "BLUP (SD units)",
         colour = "above av` & stable", title = glue("Reaction norms ({specs})")) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# biplots
plot_fa_biplot = function(pc, specs, pcs = c(1, 2), label = c("extreme", "all", "none"), n_label = 20,
                           arrow_scale = NULL) {
  label = match.arg(label)
  xs = paste0("PC", pcs[1]); ys <- paste0("PC", pcs[2])
  # check loadings
  stopifnot(all(c(xs, ys) %in% colnames(pc$loadings)))
  sc = pc$scores %>% mutate(x = .data[[xs]], y = .data[[ys]])
  ld = tibble(Env = rownames(pc$loadings), x = pc$loadings[, xs], y = pc$loadings[, ys])

  if (is.null(arrow_scale)) # display-only constant so arrows and points share one picture
    arrow_scale = 0.8 * max(abs(c(sc$x, sc$y))) / max(abs(c(ld$x, ld$y)))
  ld = ld %>% mutate(xe = x * arrow_scale, ye = y * arrow_scale)
  pct = round(100 * pc$var_explained[pcs])

  lab = switch(label,
    extreme = sc %>% mutate(d = sqrt(x^2 + y^2)) %>% slice_max(d, n = n_label),
    all = sc, none = sc[0, ])
  
  rep_ok = requireNamespace("ggrepel", quietly = TRUE)
  txt = function(data, ...) if (rep_ok) ggrepel::geom_text_repel(data = data, max.overlaps = 40, ...) else
    geom_text(data = data, vjust = -0.8, check_overlap = TRUE, ...)

  ggplot() +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    # geom_point(data = sc, aes(x, y, size = n_tested, colour = reliable), alpha = 0.85) +
    # scale_colour_manual(values = c("TRUE" = "#2166ac", "FALSE" = "grey60"),
    #                     labels = c("TRUE" = "tested in enough envs", "FALSE" = "tested in few envs"), name = NULL) +
    # scale_size_continuous(range = c(1.5, 5), name = "environments\ntested") +

    geom_point(data = sc, aes(x, y, size = n_tested)) +
  
    geom_segment(data = ld, aes(x = 0, y = 0, xend = xe, yend = ye), colour = "#b2182b",
                 arrow = arrow(length = unit(0.2, "cm"))) +
    txt(ld, aes(xe, ye, label = Env), colour = "#b2182b", fontface = "bold", size = 3.5) +
    txt(lab, aes(x, y, label = Genotype), size = 3, colour = "grey20") +
    coord_fixed() +
    labs(x = paste0(xs, " genotype score  (", pct[1], "% of common variance)"),
         y = paste0(ys, " genotype score  (", pct[2], "% of common variance)"),
         title = paste0("Genotype-environment biplot: ", xs, " vs ", ys),
         subtitle = glue("Points = genotypes, red arrows = environments.\n{specs}")) +
    theme_minimal()
}


plot_uniplot = function(loadings, cov, specs, pcs = c(1, 2)) {
  xs = paste0("PC", pcs[1]); ys <- paste0("PC", pcs[2])
  stopifnot(all(c(xs, ys) %in% colnames(loadings)))
  envs = rownames(loadings)
  stopifnot(all(envs %in% rownames(cov)))
 
  total_var = diag(cov)[envs]
  lam = loadings[, c(xs, ys), drop = FALSE] / sqrt(total_var)
  communality = rowSums(lam^2)
 
  d = tibble::tibble(Env = envs, x = lam[, 1], y = lam[, 2], communality = communality)
  circle = tibble::tibble(t = seq(0, 2 * pi, length.out = 200), x = cos(t), y = sin(t))
 
  rep_ok = requireNamespace("ggrepel", quietly = TRUE)
  txt = function(data, ...) if (rep_ok) ggrepel::geom_text_repel(data = data, max.overlaps = 40, ...) else
    geom_text(data = data, vjust = -0.8, check_overlap = TRUE, ...)
 
  ggplot(d, aes(x, y)) +
    geom_path(data = circle, aes(x, y), colour = "#2166ac", inherit.aes = FALSE) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_segment(aes(x = 0, y = 0, xend = x, yend = y, colour = communality),
                 arrow = arrow(length = unit(0.15, "cm"))) +
    txt(d, aes(x, y, label = Env), size = 3, colour = "grey20") +
    scale_colour_gradient(low = "grey70", high = "#b2182b", limits = c(0, 1), name = "communality") +
    coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
    labs(x = xs, y = ys,
         title = glue("Uniplot: environment genetic correlations ({specs})")) +
    theme_minimal()
}


# testing normal/frost years
plot_group_means = function(blup_long, env_groups, specs, tested = NULL, label_n = 30) {

  stopifnot(length(env_groups) >= 2)
  nm = names(env_groups)[1:2]
  d = blup_long %>% dplyr::select(Genotype, Env, blup)

  if (!is.null(tested)) d <- semi_join(d, tested, by = c("Genotype", "Env"))   # real data only
  
  gm = lapply(nm, function(g) d %>% filter(Env %in% env_groups[[g]]) %>% group_by(Genotype) %>%
                 summarise(m = mean(blup), n = n(), .groups = "drop") %>% rename(!!g := m, !!paste0("n_", g) := n))
  
  w = full_join(gm[[1]], gm[[2]], by = "Genotype")
  miss = w %>% filter(is.na(.data[[nm[1]]]) | is.na(.data[[nm[2]]]))
  if (nrow(miss) > 0) message(nrow(miss), " genotype(s) dropped: no tested environment in one of the groups")
  
  w = w %>% filter(!is.na(.data[[nm[1]]]), !is.na(.data[[nm[2]]])) %>%
    mutate(x = .data[[nm[1]]], y = .data[[nm[2]]], n_total = .data[[paste0("n_", nm[1])]] + .data[[paste0("n_", nm[2])]])
  lab = w %>% mutate(d = sqrt(x^2 + y^2)) %>% slice_max(d, n = label_n)

  p = ggplot(w, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(size = n_total), colour = "#2166ac", alpha = 0.8) +
    scale_size_continuous(range = c(1.5, 5), name = "tested cells\n(both groups)") +
    labs(x = paste0("mean BLUP in '", nm[1], "' environments (SD units)"),
         y = paste0("mean BLUP in '", nm[2], "' environments (SD units)"),
         title = paste0("Performance: ", nm[1], " vs ", nm[2], " environments"),
         subtitle = glue("Dashed line = equal in both; above it = relatively better in the frost envts\n{specs}")) +
    theme_minimal()

  if (requireNamespace("ggrepel", quietly = TRUE)) p + ggrepel::geom_text_repel(data = lab, aes(label = Genotype), size = 3, max.overlaps = 40)
  else p + geom_text(data = lab, aes(label = Genotype), size = 3, vjust = -0.8, check_overlap = TRUE)
}


plot_trio_pair_survival = function(pair_summary, specs) {

  d = pair_summary %>%
    pivot_longer(c(cor_baseline_with_all3, cor_with_third_dropped),
                 names_to = "which", values_to = "r") %>%
    mutate(which = recode(which, cor_baseline_with_all3 = "all 3 in model", cor_with_third_dropped = "3rd dropped"))
  
  ggplot(d, aes(which, r, group = interaction(pair, order), colour = pair)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~order) +
    labs(x = NULL, y = "rg of the remaining pair",
         title = specs)+
    theme_minimal()
}


plot_predictive_ability = function(pred_tab, specs) {
  env_order = pred_tab %>%
    group_by(dropped_env) %>%
    summarise(m = mean(r_p, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>% pull(dropped_env)
 
  d = pred_tab %>%
    mutate(dropped_env = factor(dropped_env, levels = env_order)) %>%
    tidyr::pivot_longer(c(r_p, nrmse), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric, r_p = "predictive ability (r_p)", nrmse = "NRMSE (rmse / env sd)"))
 
  ggplot(d, aes(dropped_env, value, fill = metric)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    geom_hline(data = data.frame(metric = "NRMSE (rmse / env sd)", yint = 1),
               aes(yintercept = yint), linetype = 2, colour = "grey40") +
    facet_grid(metric ~ order, scales = "free_y") +
    labs(x = NULL, y = NULL,
         title = glue("Predictive-ability LOEO ({specs})"),
         subtitle = paste("dashed line = NRMSE of 1 (Normalized RMSE)")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



plot_watch_genotype_tracking = function(watch_summary) {
  ggplot(watch_summary, aes(config, OP, colour = Genotype, group = Genotype)) +
    geom_line() + geom_point(aes(size = n_tested)) +
    facet_wrap(~order) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "NC5289 / MS2162 / SanJoaquin: OP across LOEO configs",
         x = NULL)
}