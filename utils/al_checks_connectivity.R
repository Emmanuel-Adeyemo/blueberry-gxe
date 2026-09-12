
source("R/R/visualize_connectivity.R")

al_node_counts = tibble::tribble(
  ~environment, ~total_checks,
  "BW2023", 1,
  "BW2024", 13,
  "BW2025", 13,
  "BW2026", 12,
  "EV2023", 9,
  "EV2024", 21,
  "EV2025", 20,
  "EV2026", 11
)

al_edge_overlaps = tibble::tribble(
  ~env_a, ~env_b, ~overlap,
  "BW2026", "BW2025", 12,
  "BW2026", "EV2026", 10,
  "BW2026", "BW2024", 12,
  "BW2026", "EV2025", 12,
  "BW2026", "BW2023", 1,
  "BW2026", "EV2024", 12,
  "BW2026", "EV2023", 1,
  "BW2025", "EV2026", 10,
  "BW2025", "BW2024", 13,
  "BW2025", "EV2025", 13,
  "BW2025", "BW2023", 1,
  "BW2025", "EV2024", 13,
  "BW2025", "EV2023", 1,
  "EV2026", "BW2024", 10,
  "EV2026", "EV2025", 11,
  "EV2026", "BW2023", 1,
  "EV2026", "EV2024", 11,
  "EV2026", "EV2023", 2,
  "BW2024", "EV2025", 13,
  "BW2024", "BW2023", 1,
  "BW2024", "EV2024", 13,
  "EV2025", "BW2023", 1,
  "EV2025", "EV2024", 20,
  "EV2025", "EV2023", 8,
  "BW2023", "EV2024", 1,
  "BW2023", "EV2023", 1,
  "EV2024", "EV2023", 9
)

p = plot_check_connectivity_network(
  node_counts = al_node_counts,
  edge_overlaps = al_edge_overlaps,
  title = "AL Environments: check-sharing connectivity across 8 trial-years",
  subtitle = "Node size = total checks in that environment  |  Edge width/darkness = number of shared checks"
)

ggplot2::ggsave("outputs/plots/al_checks_connectivity_network.png", p, width = 10, height = 7.5, dpi = 200, bg = "white")
