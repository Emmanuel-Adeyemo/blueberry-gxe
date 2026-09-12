
stage03_targets = list(
  tar_target(
    al_bloom_long_clean,
    clean_al_trait_long(al_bloom_long, al_metadata, al_off_type)
  ),
  tar_target(
    al_yield_long_clean,
    clean_al_trait_long(al_yield_long, al_metadata, al_off_type)
  )
)