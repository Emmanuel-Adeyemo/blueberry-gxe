
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
  ),

  # nc and fl added 10/03
  # used al dta to get shb only for fl and nc
  tar_target(nc_yield_clean, clean_nc_genotypes(nc_yield, al_yield_long_clean)),
  tar_target(nc_bloom_clean, clean_nc_genotypes(nc_bloom, al_yield_long_clean)),

  tar_target(fl_yield_clean, clean_fl_genotypes(fl_yield, al_yield_long_clean)),
  tar_target(fl_bloom_clean, clean_fl_genotypes(fl_bloom, al_yield_long_clean))
)