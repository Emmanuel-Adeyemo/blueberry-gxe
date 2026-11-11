# 08/28/2026

library(targets)

for (f in list.files('R/R', pattern = '\\.R$', full.names = TRUE)) {
  source(f)
}

for (f in list.files('R/pipeline', pattern = 'stage.*\\.R$', full.names = TRUE)){
  source(f)
}

tar_option_set(
  packages = c(
    'dplyr', 'tidyr', 'readr', 'readxl', 'stringr', 'tibble', 'purrr', 'glue', 'asreml', 'sommer', 'ggplot2',
    'ggrepel', 'igraph', 'AGHmatrix', 'stringr', 'chillR', 'readxl', 'statgenGxE', 'factoextra', 'ggcorrplot'
  )
)


c(
  stage01_targets,
  stage02_targets,
  stage03_targets,
  stage04_targets,
  stage05_targets,
  stage05b_targets,
  stage06_targets,
  stage07a_targets,
  stage07b_targets,
  stage07c_targets,
  stage07d_targets,
  stage07e_targets
)


