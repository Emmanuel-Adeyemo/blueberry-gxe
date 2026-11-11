stage03_targets = list(
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
  tar_target(nc_yield_clean, clean_nc_genotypes(nc_yield, al_yield_long_clean, nc_meta)),
  tar_target(nc_bloom_clean, clean_nc_genotypes(nc_bloom, al_yield_long_clean)),

  tar_target(fl_yield_clean, clean_fl_genotypes(fl_yield, al_yield_long_clean)),
  tar_target(fl_bloom_clean, clean_fl_genotypes(fl_bloom, al_yield_long_clean)),

  # 10/16 quality traits 
  tar_target(al_plant_height_long_clean, clean_al_trait_long(al_plant_height_long, al_metadata, al_off_type)),
  tar_target(al_avg_berry_weight_long_clean, clean_al_trait_long(al_avg_berry_weight_long, al_metadata, al_off_type)),
  tar_target(al_firmness_long_clean, clean_al_trait_long(al_firmness_long, al_metadata, al_off_type)),
  tar_target(al_brix_long_clean, clean_al_trait_long(al_brix_long, al_metadata, al_off_type)),
  tar_target(al_ta_long_clean, clean_al_trait_long(al_ta_long, al_metadata, al_off_type)),

  tar_target(nc_avg_berry_weight_clean, clean_nc_genotypes(nc_avg_berry_weight, al_yield_long_clean)),
  tar_target(nc_firmness_clean, clean_nc_genotypes(nc_firmness, al_yield_long_clean)),
  tar_target(nc_brix_clean, clean_nc_genotypes(nc_brix, al_yield_long_clean)),
  tar_target(nc_ta_clean, clean_nc_genotypes(nc_ta, al_yield_long_clean))
)