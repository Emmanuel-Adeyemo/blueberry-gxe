
# 09/01/2026
# goal here is to untangle pedigree strings and resolve synonyms. I will build a synonym table with 'identifier' to compare against.
# inputs are the six ingested pedigree files from R/R/ingest_pedigree.R from stage00.
# two are already structured, UF and NCSU, with separate male and female parent columns. For that, nothing to really process
# the others are really messy and tangled esp SHB_Ebrahiem. I will use known patterns and flag what I can't resolve for now.

library(dplyr)
library(purrr)
library(stringr)
library(tibble)


# resolve name, MS and UGA has this
.PEDIGREE_OP = c(
  "op", "open pollination", "open-pollination", "open pollinated"
)

# lines with no pedigree info. this are diff from those with open-pollination
.PEDIGREE_UNKNOWN = c('unknown', 'unk', 'na', 'n/a', '.')


# find whitespaces eg 'Un known'
.looks_like_unknown = function(s){
  compact = str_replace_all(tolower(s), "\\s+", "")
  compact %in% .PEDIGREE_UNKNOWN || (nchar(s) <= 25 && str_detect(compact, 'unknown'))
}


.top_level_matches = function(s, pattern){
  locs = str_locate_all(s, pattern)[[1]]
  if(nrow(locs) == 0){
    return(locs)
  }
  chars = str_split(s, "")[[1]]
  depth_after = integer(length(chars))
  d = 0L
  for (i in seq_along(chars)) {
    if (chars[i] %in% c("(", "[")) d = d + 1L
    if (chars[i] %in% c(")", "]")) d = d - 1L
    depth_after[i] = d
  }
  depth_before_pos = function(p) if (p <= 1) 0L else depth_after[p - 1]

  is_top_level = vapply(locs[, "start"], function(p) depth_before_pos(p) == 0L, logical(1))
  locs[is_top_level, , drop = FALSE]
}


.unwrap_if_fully_wrapped = function(s) {
  if (!str_detect(s, "^\\(") || !str_detect(s, "\\)$")) {
    return(s)
  }
  chars = str_split(s, "")[[1]]
  depth = 0
  for (i in seq_along(chars)) {
    if (chars[i] == "(") depth <- depth + 1
    if (chars[i] == ")") depth <- depth - 1
    if (depth == 0) {
      return(if (i == length(chars)) str_trim(str_sub(s, 2, -2)) else s)
    }
  }
  s
}


# GM37 (Jersey x Pioneer)/CU5 problem
# .count_top_level = function(s, pattern) {
#   stripped = str_remove_all(s, "\\([^()]*\\)")
#   str_count(stripped, pattern)
# }

# GM37 (Jersey x Pioneer)/CU5 problem
# strips the (Jersey x Pioneer)
# leaves the clean GM37 and put (Jersey x Pioneer) in the note column
.strip_trailing_parent = function(x) {
  m = str_match(x, "^(.*?)\\s*\\(([^()]*)\\)\\s*$")
  if (!is.na(m[1, 2]) && m[1, 2] != "") {
    list(clean = str_trim(m[1, 2]), note = str_trim(m[1, 3]))
  } else {
    list(clean = x, note = NA_character_)
  }
}

.is_open_pollination = function(x) {
  tolower(x) %in% .PEDIGREE_OP
}


# parse one raw pedigree string into a structured parent pair. Add status to indicate confidence
# thanks messy SHB Ebrahiem!!
.parse_one_pedigree = function(s) {
  blank = function(status) {
    tibble(
      female_parent = NA_character_, male_parent = NA_character_,
      female_parent_note = NA_character_, male_parent_note = NA_character_,
      parse_status = status
    )
  }

  if (is.na(s) || str_trim(s) == "") {
    return(blank("missing"))
  }
  s = str_trim(s)

  # if (tolower(s) %in% .PEDIGREE_UNKOWN) {
  #   return(blank("unknown_pedigree"))
  # }

  if(.looks_like_unknown(s)){
    return(blank('unknown_pedigree'))
  }

  # unwrap srings with parenthesis eg '(FL01-26/Windsor)'. This is important before anyoher pattern matching
  s = .unwrap_if_fully_wrapped(s)
 
   # self pollination
  self_match = str_match(s, regex("^(.+?)\\s+self[- ]?(pollinat\\w*|ed)?$", ignore_case = TRUE))
  if (!is.na(self_match[1, 2])) {
    parent = str_trim(self_match[1, 2])
    return(tibble(
      female_parent = parent, male_parent = parent,
      female_parent_note = NA_character_, male_parent_note = NA_character_,
      parse_status = "self_pollination"
    ))
  }

  # complex hybrids eg V. corymbosum/V. darrowii" in SHB_Ebrahiem.
  # flag it till I figure out what to do
  # TODO: SEE FLAGED
  looks_like_species_name = str_detect(s, regex("\\bhybrid\\b|\\bcomplex\\b", ignore_case = TRUE)) &&
    str_detect(s, "\\b[A-Z]\\.\\s") # eg spp name abbr like V. darrowii
  if (looks_like_species_name) {
    return(blank("unparseable_description"))
  }

  # eg TX52 - prefer pick is '/', to prevent incorrectly spliting T and 52 using x
  # flag as 'ambiguous_unparsed' rather than guess
  x_delim_pattern <- regex("\\s[x\u00d7]\\s", ignore_case = TRUE)

  slash_matches = .top_level_matches(s, fixed("/"))
  x_matches = .top_level_matches(s, x_delim_pattern) # regex("\\sx\\s", ignore_case = TRUE))

  if (nrow(slash_matches) == 1) {
    loc = slash_matches[1, ]
    # parts = str_split(s, fixed("/"), n = 2)[[1]]
    delimiter_status = "parsed_slash"
  } else if (nrow(slash_matches) == 0 && nrow(x_matches) == 1) {
    loc = x_matches[1, ]
    # parts = str_split(s, regex("\\sx\\s", ignore_case = TRUE), n = 2)[[1]]
    delimiter_status = "parsed_x_delimiter"
  } else {
    return(blank("ambiguous_unparsed"))
  }

  parts = c(str_sub(s, 1, loc[['start']] - 1), str_sub(s, loc[['end']] + 1))

  if (any(str_trim(parts) == '')) {
    return(blank('ambiguous_unparsed'))
  }

  # if (length(parts) != 2 || any(str_trim(parts) == "")) {
  #   return(blank("ambiguous_unparsed"))
  # }

  f = .strip_trailing_parent(str_trim(parts[1]))
  m = .strip_trailing_parent(str_trim(parts[2]))

  female_parent = if (.is_open_pollination(f$clean)) NA_character_ else f$clean
  male_parent = if (.is_open_pollination(m$clean)) NA_character_ else m$clean

  status = if (.is_open_pollination(f$clean) || .is_open_pollination(m$clean)) {
    paste0(delimiter_status, "_open_pollination")
  } else if (!is.na(f$note) || !is.na(m$note)) {
    paste0(delimiter_status, "_nested_generation")
  } else {
    delimiter_status
  }

  tibble(
    female_parent = female_parent, male_parent = male_parent,
    female_parent_note = f$note, male_parent_note = m$note,
    parse_status = status
  )
}


# wrapper around .parse_one_ped
#' @param x - character vector of raw pedigree 
#' @output - tibble with one row per x of female parent, male parent, female_parent_note, male_parent_note, and parse status
parse_pedigree_string = function(x) {
  map_dfr(x, .parse_one_pedigree)
}


# SYNONYM RESOLUTION
# captures 'Name=Synonym' pair
.looks_like_plausible_name = function(x){
  !str_detect(x, '[?,]') & str_count(x, fixed('(')) == str_count(x, fixed(')')) & nchar(x) <= 30
}

#' @param x - character vector - may contain NA
#' @output - tibble(name_a, name_b), one per row of only matched element
extract_synonym_pairs = function(x) {
  x = x[!is.na(x)]
  m = str_match(x, "^([^=]+)=([^=]+)$")
  ok = !is.na(m[, 1])
  name_a = str_trim(m[ok, 2])
  name_b = str_trim(m[ok, 3])
  tibble(name_a = name_a, name_b = name_b) %>%
    filter(.looks_like_plausible_name(name_a) & .looks_like_plausible_name(name_b))
}


# line_id is itself ie a 'Name=Synonym' eg in UGA ancestors, 'Rebel = TH-642', replace with Rebel as identifier.
# extract_synonym_pairs() already captured pairs, this gives the row a single identifier to carry till apply_synonym_map()
split_double_names = function(x) {
  m = str_match(x, "^([^=]+)=([^=]+)$")
  name_a = str_trim(m[, 2])
  name_b = str_trim(m[, 3])
  is_plausible_pair = !is.na(m[, 1]) &
    .looks_like_plausible_name(name_a) &
    .looks_like_plausible_name(name_b)
  if_else(is_plausible_pair, name_a, x)
}

# for chains like A=B, B=C ie A=B=C
.uf_new = function(names) {
  env = new.env(parent = emptyenv())
  for (nm in names) assign(nm, nm, envir = env)
  env
}

.uf_find = function(env, x) {
  while (!identical(get(x, envir = env), x)) {
    x = get(x, envir = env)
  }
  x
}

.uf_union = function(env, a, b) {
  ra = .uf_find(env, a)
  rb = .uf_find(env, b)
  if (!identical(ra, rb)) {
    assign(ra, rb, envir = env)
  }
}

# build alias
# TODO: talk to Sushan about the correct name, right now - using a deterministic default (alphabetical)
#' @param pairs - tibble(name_a, name_b) ie output of extract_synonym_pairs()
#' @output - tibble(alias, identifier)
build_synonym_groups = function(pairs) {
  if (nrow(pairs) == 0) {
    return(tibble(alias = character(0), identifier = character(0)))
  }

  all_names = unique(c(pairs$name_a, pairs$name_b))
  uf = .uf_new(all_names)
  for (i in seq_len(nrow(pairs))) {
    .uf_union(uf, pairs$name_a[i], pairs$name_b[i])
  }

  roots = vapply(all_names, function(nm) .uf_find(uf, nm), character(1))

  tibble(alias = all_names, group_root = roots) %>%
    group_by(group_root) %>%
    mutate(identifier = sort(alias)[1]) %>%
    ungroup() %>%
    select("alias", "identifier")
}

# rewrite character vec to identifier name, ie based on synonym_map. values not present in map remain unchanged
#' @param x - character vector 
#' @param synonym_map - tibble(alias, identifier)
#' @output - character vector, same len as x
apply_synonym_map = function(x, synonym_map) {
  idx = match(x, synonym_map$alias)
  if_else(is.na(idx), x, synonym_map$identifier[idx])
}


# Workflow entry - combines all six ped files into one parsed, synonym processed file
# input1: uf_pedigree, ncsu_pedigree - already-structured stage00 tables
#   (source, line_id, female_parent, male_parent, ...)
# input2: uga_crosses_pedigree, uga_ancestors_pedigree, ms_usda_pedigree,
#   shb_ebrahiem_pedigree - unparsed  stage00 tables (source, line_id,
#   pedigree_raw, ...)

combine_and_resolve_pedigree = function(uf_pedigree,
                                          ncsu_pedigree,
                                          uga_crosses_pedigree,
                                          uga_ancestors_pedigree,
                                          ms_usda_pedigree,
                                          shb_ebrahiem_pedigree) {
  
  # these are already structured
  already_structured = bind_rows(
    uf_pedigree %>% select("source", "line_id", "female_parent", "male_parent"),
    ncsu_pedigree %>% select("source", "line_id", "female_parent", "male_parent")
  ) %>%
    mutate(
      pedigree_raw = NA_character_,
      female_parent_note = NA_character_,
      male_parent_note = NA_character_,
      parse_status = "already_structured"
    )

    # for the problem children
  parse_source = function(dta_table) {
    bind_cols(
      dta_table %>% select("source", "line_id", pedigree_raw),
      parse_pedigree_string(dta_table$pedigree_raw)
    )
  }

  parsed_string_sources = bind_rows(
    parse_source(uga_crosses_pedigree),
    parse_source(uga_ancestors_pedigree),
    parse_source(ms_usda_pedigree),
    parse_source(shb_ebrahiem_pedigree)
  )

  combined = bind_rows(already_structured, parsed_string_sources)

  #fix synonym 
  # 'Name=Synonym' sometimes appear as line_id itself especially in UGA ancestors
  synonym_pairs = bind_rows(
    extract_synonym_pairs(combined$line_id),
    extract_synonym_pairs(combined$female_parent),
    extract_synonym_pairs(combined$male_parent)
  ) %>%
    distinct()

  synonym_map = build_synonym_groups(synonym_pairs)

  pedigree_processed = combined %>%
    mutate(
      line_id_raw = line_id,
      female_parent_raw = female_parent,
      male_parent_raw = male_parent,
      
      # happens here
      line_id = split_double_names(line_id),
      female_parent = split_double_names(female_parent),
      male_parent = split_double_names(male_parent)
    ) %>%
    mutate(
      line_id = apply_synonym_map(line_id, synonym_map),
      female_parent = apply_synonym_map(female_parent, synonym_map),
      male_parent = apply_synonym_map(male_parent, synonym_map)
    ) %>%
    select(
      "source", "line_id", "female_parent", "male_parent",
      "female_parent_note", "male_parent_note", "parse_status","pedigree_raw",
      "line_id_raw", "female_parent_raw", "male_parent_raw"
    )

  list(pedigree_processed = pedigree_processed, synonym_map = synonym_map)
}
