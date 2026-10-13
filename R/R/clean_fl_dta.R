# 10/2

library(dplyr)
library(stringr) 


FL_GENOTYPE_ALIASES = c(
  'Albus' = '11-51',
  'New Hannover' = 'NewHanover',
  'New Hanover' = 'NewHanover'
)


normalize_fl_genotype = function(genotype) {
  recoded = FL_GENOTYPE_ALIASES[genotype]
  x = if_else(!is.na(recoded), recoded, genotype)
  needs_prefix = str_detect(x, '^[0-9]+-[0-9]+[A-Za-z]?$')
  if_else(needs_prefix, paste0('FL', x), x)
}


clean_fl_genotypes = function(dta, al_dta) {
  before = dta$genotype
  dta$genotype = normalize_fl_genotype(before)

  unique_genotypes = unique(al_dta$Genotype)

  changed = before != dta$genotype
  if (any(changed)) {
    changes <- unique(paste0(before[changed], " -> ", dta$genotype[changed]))
    message("clean_fl_genotypes(): resolved ", length(changes), " names: ", paste(changes, collapse = ", "))
  }
  dta %>% filter(genotype %in% unique_genotypes)
}


