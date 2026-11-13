# 10/02

library(dplyr)
library(stringr)

NC_GENOTYPE_ALIASES = c(
  "O'Neal" = "ONeal",
  'NC5326B' = 'NC5326',
  'New Hanover' = 'NewHanover'
)



#' @description normalizes NC genotype names for matching with AL genotpe set.
#' @param genotype character vector of raw NC genotype names e.g. the `genotype` col of `nc_yield`, output of `ingest_nc()`
#' @output character vector of normalized geno names 
normalize_nc_genotype = function(genotype) {
  recoded = unname(NC_GENOTYPE_ALIASES[genotype])
  if_else(!is.na(recoded), recoded, genotype)
}



#' @description cleans NC genotype names with `normalize_nc_genotype`, drops any geno not in the AL/SHB set
#' @param dta NC trait dta with a `genotype` col e.g. `nc_yield`, output of `ingest_nc()`
#' @param al_dta al dta, already cleaned. used as the reference genotype set
#' @param nc_meta optional. NC trial metadata e.g. `nc_meta`, output of `ingest_nc()`
#' @output df with genotype normalized and filtered to genotypes in the reference set
clean_nc_genotypes = function(dta, al_dta, nc_meta = NULL) {
  before = dta$genotype
  dta$genotype = normalize_nc_genotype(before)

  unique_genotypes = unique(al_dta$Genotype)

  changed = before != dta$genotype
  if (any(changed)) {
    changes = unique(paste0(before[changed], " -> ", dta$genotype[changed]))
    message("clean_nc_genotypes(): resolved ", length(changes), " alias(es): ", paste(changes, collapse = ", "))
  }
  # added subpot filter 10/12
  dta = dta %>% filter(genotype %in% unique_genotypes) %>% filter(!subplot %in% c(4:6))
  if(!is.null(nc_meta)) dta %>% inner_join(nc_meta, by = c('plot', 'genotype')) %>% mutate(Age = season - planting_year) %>% dplyr::select(-planting_year) else dta
}
