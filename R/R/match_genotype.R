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
normalize_pedigree_name <- function(x) {
  x %>%
    tolower() %>%
    str_remove_all("['\u2019]") %>%
    str_replace_all("\\s+", "")
}


# fix two - add prefix for UF's convention
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
match_trial_genotypes = function(trail_genotypes, processed_pedigree) {
  trial_keys = build_match_keys(trial_genotypes) %>%
    select('trial_genotype' = 'name', 'match_key')

  pedigree_keys = build_match_keys(processed_pedigree$line_id, processed_pedigree$source) %>%
    select('line_id' = 'name', 'source', 'match_key', 'key_type')

  trial_keys %>% 
    inner_join(pedigree_keys, by = 'match_key', relationship = 'many-to-many') %>% 
    distinct(.data$trial_genotype, .data$line_id, .data$source, .data$key_type) %>%
    arrange(.data$trial_genotype, .data$key_type)
}