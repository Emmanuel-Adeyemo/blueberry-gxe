# 09/22

library(dplyr)
library(chillR)

compute_daily_chill_heat = function(hour_temps_one_loc) {

  hourtemps_clean = hour_temps_one_loc %>% dplyr::select(-any_of('Location'))

  results = chillR::daily_chill(
    hourtemps = as.data.frame(hour_temps_one_loc),
    running_mean = 1,
    models = list(Chill_Portions = chillR::Dynamic_Model, GDH = chillR::GDH )
  )
  as_tibble(results$daily_chill)
}


# plot envt covariates
plot_envt_cluster = function(env_cov){

  Env_label = env_cov %>%
    mutate(Env = paste(location, season)) %>%
    pull(Env)

  clust_vars = env_cov %>%
    dplyr::select(Chill_Portions_sum, GDH_sum, frost_score) %>%
    scale()

  rownames(clust_vars) = Env_label

  dist_mat = dist(clust_vars, method = "euclidean")
  hc = hclust(dist_mat, method = "ward.D2")
  min_y = -0.25 * max(hc$height)

  fviz_dend(
    hc, 
    main = "Env clustering: chill + heat + frost",
    ggtheme = theme_minimal()
  ) +
    scale_y_continuous(limits = c(min_y, NA)) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = margin(t = 10, r = 10, b = 40, l = 10, unit = 'pt'))
    
}

