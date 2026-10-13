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

  # NC and FL added 10/03
  tar_target(
    nc_path25, 'data/raw/nc/2025NCSUData_updated_09282026.xlsx', format = 'file'
  ),

  tar_target(nc_path24, 'data/raw/nc/2024 NCSU Data.xlsx', format = 'file'),

  tar_target(fl_path25, 'data/raw/fl/UF 2025.xlsx', format = 'file'),

  tar_target(fl_path24_yld, 'data/raw/fl/UF 2024 yield.xlsx', format = 'file'),

  tar_target(fl_path24_bloom, 'data/raw/fl/Bloom_UF_2024.xlsx', format = 'file'),

  tar_target(al_rda, read_al_rda(al_rda_file)),
  
  tar_target(al_metadata, read_al_metadata(al_rda)),
  tar_target(al_off_type, read_al_off_type(al_rda)),
  
  tar_target(al_all_dta, ingest_all_al_dta(al_rda)),
  tar_target(al_bloom_long, al_all_dta$bloom_long),
  tar_target(al_yield_long, al_all_dta$yield_long),

  # nc
  tar_target(nc_all_dta, ingest_nc(nc_path25, nc_path24)),
  tar_target(nc_yield, nc_all_dta$yield),
  tar_target(nc_bloom, nc_all_dta$bloom),

  # fl
  tar_target(fl_all_dta, ingest_fl(fl_path25, fl_path24_yld, fl_path24_bloom)),
  tar_target(fl_yield, fl_all_dta$yield),
  tar_target(fl_bloom, fl_all_dta$bloom),

  tar_target(trial_genotypes, get_unique_genotypes(al_metadata))
  
)

