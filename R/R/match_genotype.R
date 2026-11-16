# 09/02/2026
# match external genotypes against processed pedigree from stage01, builds matching keys that survives formatting noise
# two fixes: 1. normalize key types of pure formatting noise eg whitespaces, apostrophes. 2. Fix UF's convention that omits UF prefix.
# TODO - in the process_pedigree.R - figure out the fuzzy matching e.g. "Springhigh" vs "Springhieght" or cross-source disagreement (e.g. Gumbo's
# female parent: "Cara'sChoice" (NCSU) vs "G695"/"G 695" (MS_USDA, SHB_Ebrahiem)). Also review with Sushan and Sakshi before deciding

library(dplyr)
library(stringr)
library(tibble)


.SOURCE_PREFIX = c(UF = 'FL')

# fix one - normalize formatting noise
#' @description normalizes a genotype/line name into a matchable key,lowercases it, strips apostrophes, and removes all whitespace
#' @param x character vec of raw genotype/line names 
#' @output a character vec of normalized keys, same length as `x`.
normalize_pedigree_name = function(x) {
  x %>%
    tolower() %>%
    str_remove_all("['\u2019]") %>%
    str_replace_all("\\s+", "")
}


# fix two - add prefix for UF's convention
#' @description builds one or two candidate match keys per name: a plain normalized key and the institution prefix eg FL for UF genos
#' @param names character vec of genotype/line names to build keys for e.g. `trial_genotypes` in `match_trial_genotypes()`, or `processed_pedigree$line_id`
#' @param sources character vec the same length as `names` 
#' @output df
build_match_keys = function(names, sources = NA_character_) {
  if(length(sources) == 1) {
    sources = rep(sources, length(names))
  }
  if (length(sources) != length(names)) {
    stop('build_match_keys(): "source" must be len 1 or len(names).')
  }

  base = tibble = tibble(
    name = names,
    source = sources,
    match_key = normalize_pedigree_name(names),
    key_type = 'normalized'
  )

  has_prefix = sources %in% names(.SOURCE_PREFIX)
  prefixed = tibble(
    name = names[has_prefix],
    source = sources[has_prefix],
    match_key = normalize_pedigree_name(
      paste0(.SOURCE_PREFIX[sources[has_prefix]], names[has_prefix])
    ),
    key_type = 'institution_prefix'
  )
  bind_rows(base, prefixed) %>% 
    distinct()
}


# do matching
#' @description matches trial genotype names against the processed pedigree's `line_id`s by building normalized/institution prefix and inner-joining on `match_key`.
#' @param trial_genotypes character vec of trial genotype names. In the pipeline this is the `trial_genotypes` target, i.e. `get_unique_genotypes(al_metadata)` 
#' @param processed_pedigree a pedigree table with `line_id` and `source` columns. 
#' @output df
match_trial_genotypes = function(trial_genotypes, processed_pedigree) {
  trial_keys = build_match_keys(trial_genotypes) %>%
    dplyr::select('trial_genotype' = 'name', 'match_key')

  pedigree_keys = build_match_keys(processed_pedigree$line_id, processed_pedigree$source) %>%
    dplyr::select('line_id' = 'name', 'source', 'match_key', 'key_type')

  trial_keys %>% 
    inner_join(pedigree_keys, by = 'match_key', relationship = 'many-to-many') %>% 
    distinct(trial_genotype, line_id, source, key_type) %>%
    arrange(trial_genotype, key_type)
}