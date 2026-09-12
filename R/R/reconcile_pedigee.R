# 09/02/2026

library(dplyr)
library(stringr)
library(tibble)

.pick_identifier_value <- function(values, sources) {
  keep <- !is.na(values) & str_trim(values) != ''
  values <- values[keep]
  sources <- sources[keep]

  if (length(values) == 0) {
    return(tibble(
      canonical_value = NA_character_, n_sources_agreeing = 0L,
      n_sources_total = 0L, confidence = 'no_data', alternatives = NA_character_
    ))
  }

  n_sources_total = length(unique(sources))
  keys = build_match_keys(values, sources)
  
  all_keys <- unique(keys$match_key)
  uf <- .uf_new(all_keys)
  key_groups <- split(keys$match_key, paste(keys$name, keys$source, sep = "\r"))
  for (ks in key_groups) {
    ks <- unique(ks)
    if (length(ks) > 1) {
      for (i in 2:length(ks)) .uf_union(uf, ks[1], ks[i])
    }
  }

  group_root <- vapply(all_keys, function(k) .uf_find(uf, k), character(1))
  keys <- keys %>% mutate(group_id = group_root[match(match_key, all_keys)])

  # a source "supports" a group if ANY of its rows landed in it
  group_support <- keys %>%
    distinct(source, group_id) %>%
    count(group_id, name = "n_sources_agreeing") %>%
    arrange(desc(n_sources_agreeing))

  representatives <- keys %>%
    distinct(name, group_id) %>%
    mutate(
      informative_len = nchar(str_remove_all(name, "\\s+")),
      raw_len = nchar(name)
    ) %>%
    group_by(group_id) %>%
    slice_max(order_by = informative_len, n = 1, with_ties = TRUE) %>%
    slice_min(order_by = raw_len, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select("name", "group_id")

  n_groups <- nrow(group_support)
  winner <- group_support$group_id[1]
  n_agreeing <- group_support$n_sources_agreeing[1]

  confidence = if(n_sources_total == 1) {
    'single_source'
  }else if (n_groups >= 2 && n_agreeing > group_support$n_sources_agreeing[2]) {
    'majority'
  }else{
    'tied'
  }

  identifier_value = representatives$name[representatives$group_id == winner]

  alternatives = representatives %>%
    filter(group_id != winner) %>% 
    pull(name) %>% 
    paste(collapse = '; ')

  tibble(
    identifier_value = identifier_value,
    n_sources_agreeing = n_agreeing,
    n_sources_total = n_sources_total,
    confidence = confidence,
    alternatives = if(identical(alternatives, '') ) NA_character_ else alternatives
  )
}  



# collape output of match_trial_genotype() and output of stage01 pedigree table into one resolved row per genotype.
# joins in female / male parent info per line_id and source. Then resolves disagreements independently for female_parent 
# and male_parent using .pick_identifier_value()

#' @param genotype_matches - output of match_trial_genotype()
#' @param processed_pedigree - output of stage01

reconcile_trial_genotype_pedigree = function(genotype_matches, processed_pedigree) {
  joined = genotype_matches %>%
    distinct(trial_genotype, line_id, source) %>%
    left_join(
      processed_pedigree %>% select('line_id', 'source', 'female_parent', 'male_parent'),
      by = c('line_id', 'source')
    )
  
  joined %>%
    group_by(trial_genotype) %>%
    group_modify(function(rows, ...) {
      female <- .pick_identifier_value(rows$female_parent, rows$source)
      male <- .pick_identifier_value(rows$male_parent, rows$source)
      tibble(
        identifier_female_parent = female$identifier_value,
        female_parent_confidence = female$confidence,
        female_parent_alternatives = female$alternatives,
        identifier_male_parent = male$identifier_value,
        male_parent_confidence = male$confidence,
        male_parent_alternatives = male$alternatives,
        n_sources = length(unique(rows$source)),
        sources = paste(sort(unique(rows$source)), collapse = ", ")
      )
    }) %>%
    ungroup()
}



