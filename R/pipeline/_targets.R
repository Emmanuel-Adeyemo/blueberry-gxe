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
    'dplyr', 'tidyr', 'readr', 'readxl', 'stringr', 'tibble', 'purrr', 'glue', 'asreml', 'sommer'
  )
)

c(
  stage00_targets,
  stage01_targets,
  stage03_targets,
  stage04_targets,
  stage05_targets,
  stage06_targets,
  stage07_targets,
  stage07b_targets,
  stage07c_targets
)


