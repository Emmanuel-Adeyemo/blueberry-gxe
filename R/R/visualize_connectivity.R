library(dplyr)
library(ggplot2)
library(ggrepel)
library(igraph)


classify_connectivity_status = function(edge_overlaps, strong_threshold = 8, core_min_fraction = 0.5) {
  long = bind_rows(
    edge_overlaps %>% dplyr::select(environment = "env_a", partner = "env_b", "overlap"),
    edge_overlaps %>% dplyr::select(environment = "env_b", partner = "env_a", "overlap")
  )

  n_partners = long %>%
    distinct(environment, partner) %>%
    count(environment, name = "n_partners")

  long %>%
    group_by(environment) %>%
    summarise(n_strong_ties = sum(overlap >= strong_threshold), .groups = "drop") %>%
    left_join(n_partners, by = "environment") %>%
    mutate(
      strong_tie_fraction = n_strong_ties / n_partners,
      status = case_when(
        n_strong_ties == 0 ~ "isolated / weak link",
        strong_tie_fraction >= core_min_fraction ~ "well-connected core",
        TRUE ~ "partially connected"
      )
    ) %>%
    dplyr::select("environment", "n_strong_ties", "strong_tie_fraction", "status")
}




plot_check_connectivity_network = function(node_counts,
                                             edge_overlaps,
                                             status = NULL,
                                             title = "Check-sharing connectivity between trial environments",
                                             subtitle = NULL,
                                             seed = 42) {
  if (is.null(status)) {
    status = classify_connectivity_status(edge_overlaps)
  }

  node_counts = node_counts %>% left_join(status, by = "environment")

  g = graph_from_data_frame(
    edge_overlaps %>% dplyr::select("env_a", "env_b", weight = "overlap"),
    directed = FALSE,
    vertices = node_counts
  )

  set.seed(seed)
  layout_xy = layout_with_fr(g, weights = E(g)$weight)
  colnames(layout_xy) = c("x", "y")

  node_df = as_tibble(layout_xy) %>%
    mutate(environment = V(g)$name) %>%
    left_join(node_counts, by = "environment")

  edge_df = edge_overlaps %>%
    left_join(node_df %>% dplyr::select("env_a" = "environment", x_a = "x", y_a = "y"), by = "env_a") %>%
    left_join(node_df %>% dplyr::select("env_b" = "environment", x_b = "x", y_b = "y"), by = "env_b")

  status_colors = c(
    "well-connected core" = "#2C7A4B",
    "partially connected" = "#D9A441",
    "isolated / weak link" = "#B23A48"
  )

  ggplot() +
    geom_segment(
      data = edge_df,
      aes(x = x_a, y = y_a, xend = x_b, yend = y_b, linewidth = overlap, alpha = overlap),
      color = "grey40"
    ) +
    scale_linewidth(range = c(0.15, 2.8), guide = "none") +
    scale_alpha(range = c(0.15, 0.85), guide = "none") +
    geom_point(
      data = node_df,
      aes(x = x, y = y, size = total_checks, fill = status),
      shape = 21, color = "white", stroke = 1.2
    ) +
    scale_size(range = c(8, 26), name = "Total checks\nin environment") +
    scale_fill_manual(values = status_colors, name = "Connectivity status") +
    geom_text_repel(
      data = node_df,
      aes(x = x, y = y, label = paste0(environment, "\n(n=", total_checks, ")")),
      size = 3.4, fontface = "bold", segment.color = NA,
      point.padding = 0.6, box.padding = 0.4
    ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0),
      plot.subtitle = element_text(color = "grey30", margin = margin(b = 12)),
      legend.position = "right"
    )
}
