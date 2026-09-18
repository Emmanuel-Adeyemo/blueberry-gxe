
stage03_targets <- list(
  tar_target(
    al_bloom_prepared,
    exclude_never_bloomed(
      carry_forward_bloom_correction(
        clean_al_trait_long(al_bloom_long, al_metadata, al_off_type)
      )
    )
  ),
  tar_target(al_bloom_long_clean, al_bloom_prepared$bloom_dta),
  tar_target(al_bloom_excluded_dead_plants, al_bloom_prepared$excluded_plants),

  tar_target(
    al_yield_long_clean,
    clean_al_trait_long(al_yield_long, al_metadata, al_off_type)
  )
)