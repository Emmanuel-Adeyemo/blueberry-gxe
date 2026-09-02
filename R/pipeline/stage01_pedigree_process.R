
# all pedigree processing

stage01_targets = list(
  tar_target(
    pedigree_resolved,
    combine_and_resolve_pedigree(
      uf_pedigree = uf_pedigree,
      ncsu_pedigree = ncsu_pedigree,
      uga_crosses_pedigree = uga_crosses_pedigree,
      uga_ancestors_pedigree = uga_ancestors_pedigree,
      ms_usda_pedigree = ms_usda_pedigree,
      shb_ebrahiem_pedigree = shb_ebrahiem_pedigree
    )
  ), 
  tar_target(
    pedigree_processed, pedigree_resolved$pedigree_processed
  ),
  tar_target(
    pedigree_synonym_map, pedigree_resolved$synonym_map
  )
)


