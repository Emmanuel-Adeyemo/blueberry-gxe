# 08/28/2026
# AL pheno
# Single rda file with 8 trial dta frames <LOC><YY>.ori + MetaData + OffID


stage00_targets = list(
  
  # pedigree
  # UF
  tar_target(
    uf_pedigree_file,
    'data/raw/pedigree/ped_updated2.rda',
    format = 'file'
  ),
  tar_target(
    uf_pedigree, 
    read_uf_pedigree(uf_pedigree_file)
  ),
  
  # NCSU
  tar_target(
    ncsu_pedigree_file,
    'data/raw/pedigree/NCSU_PedigreeAppJS_spencer.xlsx',
    format = 'file'
  ),
  tar_target(
    ncsu_pedigree,
    read_ncsu_pedigree(ncsu_pedigree_file)
  ),
  
  # UGA
  tar_target(
    uga_pedigree_file,
    'data/raw/pedigree/Pedigree_UGA_2026_Chu.xlsx',
    format = 'file'
  ),
  tar_target(
    uga_crosses_pedigree,
    read_uga_pedigree_crosses(uga_pedigree_file)
  ),
  tar_target(
    uga_ancestors_pedigree,
    read_uga_pedigree_ancestors(uga_pedigree_file)
  ),
  
  # MS USDA
  tar_target(
    ms_usda_pedigree_file,
    'data/raw/pedigree/MS_USDA_pedigree.csv',
    format = 'file'
  ),
  tar_target(
    ms_usda_pedigree,
    read_ms_usda_pedigree(ms_usda_pedigree_file)
  ),
  
  
  # SHB Ebrahiem
  tar_target(
    shb_ebrahiem_pedigree_file,
    'data/raw/pedigree/Pedigree_of_SHB_Ebrahiem.xlsx',
    format = 'file'
  ),
  tar_target(
    shb_ebrahiem_pedigree,
    read_shb_ebrahiem_pedigree(shb_ebrahiem_pedigree_file)
  ),
  
  # loosely combine for stacked view
  tar_target(
    pedigree_all_raw,
    dplyr::bind_rows(
      uf_pedigree,
      ncsu_pedigree,
      dplyr::rename(uga_crosses_pedigree, female_parent_or_pedigree = 'pedigree_raw'),
      dplyr::rename(ms_usda_pedigree, female_parent_or_pedigree = 'pedigree_raw'),
      dplyr::rename(shb_ebrahiem_pedigree, female_parent_or_pedigree = 'pedigree_raw')
    )
  ),
  
  
  # pheno
  # AL
  tar_target(
    al_rda_file,
    'data/raw/al/OriginalData.rda',
    format='file'
  ),
  tar_target(al_rda, read_al_rda(al_rda_file)),
  
  tar_target(al_metadata, read_al_metadata(al_rda)),
  tar_target(al_off_type, read_al_off_type(al_rda)),
  
  tar_target(al_all_dta, ingest_all_al_dta(al_rda)),
  tar_target(al_bloom_long, al_all_dta$bloom_long),
  tar_target(al_berry_long, al_all_dta$berry_percent_long),
  tar_target(al_frost_long, al_all_dta$frost_percent_long),
  tar_target(al_yield_long, al_all_dta$yield_long),
  tar_target(al_health_score_long, al_all_dta$health_score_long),
  tar_target(al_plant_height_long, al_all_dta$plant_height_long),
  tar_target(al_notes_long, al_all_dta$notes_long),
  tar_target(avg_berry_weight_long, al_all_dta$avg_berry_weight_long),
  
  tar_target(firmness_long, al_all_dta$firmness_long),
  tar_target(brix_long, al_all_dta$brix_long),
  tar_target(ta_long, al_all_dta$ta_long),
  tar_target(berry_diameter_long, al_all_dta$berry_diameter_long),

  tar_target(trial_genotypes, get_unique_genotypes(al_yield_long))
  
)

