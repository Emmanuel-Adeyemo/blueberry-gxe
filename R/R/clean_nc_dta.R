# 10/02

library(dplyr)
library(stringr)

NC_GENOTYPE_ALIASES = c(
  "O'Neal" = "ONeal",
  'NC5326B' = 'NC5326',
  'New Hanover' = 'NewHanover'
)

normalize_nc_genotype = function(genotype) {
  recoded = unname(NC_GENOTYPE_ALIASES[genotype])
  if_else(!is.na(recoded), recoded, genotype)
}


clean_nc_genotypes = function(dta, al_dta) {
  before = dta$genotype
  dta$genotype = normalize_nc_genotype(before)

  unique_genotypes = unique(al_dta$Genotype)

  changed = before != dta$genotype
  if (any(changed)) {
    changes = unique(paste0(before[changed], " -> ", dta$genotype[changed]))
    message("clean_nc_genotypes(): resolved ", length(changes), " alias(es): ", paste(changes, collapse = ", "))
  }
  dta %>% filter(genotype %in% unique_genotypes)
}
