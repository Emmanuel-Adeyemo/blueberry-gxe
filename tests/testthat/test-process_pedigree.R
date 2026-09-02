test_that("parse_pedigree_string() parses a clean slash-delimited pair", {
  out <- parse_pedigree_string("FL93-70/FL90-4")
  expect_equal(out$female_parent, "FL93-70")
  expect_equal(out$male_parent, "FL90-4")
  expect_equal(out$parse_status, "parsed_slash")
  expect_true(is.na(out$female_parent_note))
})

test_that("parse_pedigree_string() falls back to ' x ' when there's no slash", {
  out <- parse_pedigree_string("E-66 x Fla 1-3")
  expect_equal(out$female_parent, "E-66")
  expect_equal(out$male_parent, "Fla 1-3")
  expect_equal(out$parse_status, "parsed_x_delimiter")
})

test_that("parse_pedigree_string() resolves both open-pollination spellings to NA male parent", {
  op <- parse_pedigree_string("TH-639/OP")
  expect_equal(op$female_parent, "TH-639")
  expect_true(is.na(op$male_parent))
  expect_equal(op$parse_status, "parsed_slash_open_pollination")

  full_phrase <- parse_pedigree_string("VERNON/Open pollination")
  expect_equal(full_phrase$female_parent, "VERNON")
  expect_true(is.na(full_phrase$male_parent))
  expect_equal(full_phrase$parse_status, "parsed_slash_open_pollination")
})

test_that("parse_pedigree_string() detects self-pollination", {
  out <- parse_pedigree_string("US75 self pollination")
  expect_equal(out$female_parent, "US75")
  expect_equal(out$male_parent, "US75")
  expect_equal(out$parse_status, "self_pollination")
})

test_that("parse_pedigree_string() flags 'Unknown' without guessing parents", {
  out <- parse_pedigree_string("Unknown")
  expect_true(is.na(out$female_parent))
  expect_true(is.na(out$male_parent))
  expect_equal(out$parse_status, "unknown_pedigree")
})

test_that("parse_pedigree_string() flags NA/blank as missing", {
  out <- parse_pedigree_string(c(NA_character_, ""))
  expect_equal(out$parse_status, c("missing", "missing"))
})

test_that("parse_pedigree_string() strips a trailing nested-generation parenthetical but keeps it as a note", {
  out <- parse_pedigree_string("GM37 (Jersey x Pioneer)/CU5 (Stanley x June)")
  expect_equal(out$female_parent, "GM37")
  expect_equal(out$female_parent_note, "Jersey x Pioneer")
  expect_equal(out$male_parent, "CU5")
  expect_equal(out$male_parent_note, "Stanley x June")
  expect_equal(out$parse_status, "parsed_slash_nested_generation")
})

test_that("parse_pedigree_string() flags species-hybrid prose instead of mis-splitting it", {
  out <- parse_pedigree_string("A complex hybrid of V. corymbosum/V. darrowii")
  expect_true(is.na(out$female_parent))
  expect_true(is.na(out$male_parent))
  expect_equal(out$parse_status, "unparseable_description")
})

test_that("parse_pedigree_string() flags ambiguous multi-delimiter strings rather than guessing", {
  out <- parse_pedigree_string("A/B/C")
  expect_true(is.na(out$female_parent))
  expect_true(is.na(out$male_parent))
  expect_equal(out$parse_status, "ambiguous_unparsed")
})

test_that("parse_pedigree_string() is vectorized and preserves order", {
  out <- parse_pedigree_string(c("A/B", "Unknown", "C/D"))
  expect_equal(nrow(out), 3)
  expect_equal(out$female_parent, c("A", NA, "C"))
  expect_equal(out$parse_status, c("parsed_slash", "unknown_pedigree", "parsed_slash"))
})

test_that("extract_synonym_pairs() finds 'Name=Synonym' patterns and ignores everything else", {
  out <- extract_synonym_pairs(c("Rebel=TH-642", "TH-639", NA, "Suziblue=TH-730"))
  expect_equal(nrow(out), 2)
  expect_setequal(out$name_a, c("Rebel", "Suziblue"))
  expect_setequal(out$name_b, c("TH-642", "TH-730"))
})

test_that("split_composite_names() extracts the left-hand name and passes other values through", {
  out <- split_composite_names(c("Rebel=TH-642", "TH-639", NA))
  expect_equal(out, c("Rebel", "TH-639", NA))
})

test_that("build_synonym_groups() resolves transitive chains (A=B, B=C -> one group)", {
  pairs <- tibble::tibble(name_a = c("A", "B"), name_b = c("B", "C"))
  out <- build_synonym_groups(pairs)

  expect_setequal(out$alias, c("A", "B", "C"))
  # alphabetically-first-in-group canonical rule -> "A"
  expect_true(all(out$canonical == "A"))
})

test_that("build_synonym_groups() returns an empty table for no pairs", {
  out <- build_synonym_groups(tibble::tibble(name_a = character(0), name_b = character(0)))
  expect_equal(nrow(out), 0)
})

test_that("apply_synonym_map() rewrites mapped values and passes through everything else, including NA", {
  synonym_map <- tibble::tibble(alias = c("Rebel", "TH-642"), canonical = c("Rebel", "Rebel"))
  out <- apply_synonym_map(c("TH-642", "Rebel", "Unrelated", NA), synonym_map)
  expect_equal(out, c("Rebel", "Rebel", "Unrelated", NA))
})

test_that("combine_and_resolve_pedigree() integrates all 6 sources and resolves a cross-source synonym", {
  pedigree_uf <- tibble::tibble(
    source = "UF", line_id = "Wolcott", female_parent = NA_character_, male_parent = NA_character_
  )
  pedigree_ncsu <- tibble::tibble(
    source = "NCSU", line_id = "GM-37", female_parent = "NorthSedgwick", male_parent = "Coville"
  )
  pedigree_uga_crosses <- tibble::tibble(
    source = "UGA", line_id = "TH-1484", pedigree_raw = "TH-639/Rebel"
  )
  # "Rebel=TH-642" here means THIS row's own identity has two names; its
  # own pedigree is separate from the synonym fact itself.
  pedigree_uga_ancestors <- tibble::tibble(
    source = "UGA_ancestors", block_id = 1L,
    line_id = "Rebel=TH-642", pedigree_raw = "FL.HB-90-04/OP"
  )
  pedigree_ms_usda <- tibble::tibble(
    source = "MS_USDA", line_id = "MS1110", pedigree_raw = "US 841/MS231"
  )
  pedigree_shb_ebrahiem <- tibble::tibble(
    source = "SHB_Ebrahiem", line_id = "Bluecrisp",
    pedigree_raw = "A complex hybrid of V. corymbosum/V. darrowii",
    observation_raw = NA_character_
  )

  out <- combine_and_resolve_pedigree(
    pedigree_uf, pedigree_ncsu, pedigree_uga_crosses,
    pedigree_uga_ancestors, pedigree_ms_usda, pedigree_shb_ebrahiem
  )

  parsed <- out$pedigree_parsed
  expect_equal(nrow(parsed), 6)

  # UF/NCSU pass through untouched, tagged as already_structured
  expect_equal(parsed$parse_status[parsed$source == "UF"], "already_structured")

  # UGA_ancestors' own composite line_id resolves to the canonical form
  # ("Rebel" < "TH-642" alphabetically)
  ancestor_row <- parsed[parsed$source == "UGA_ancestors", ]
  expect_equal(ancestor_row$line_id, "Rebel")
  expect_equal(ancestor_row$line_id_raw, "Rebel=TH-642")

  # ...and the UGA crosses row that references "Rebel" as a male parent
  # stays "Rebel" too (same canonical group, already the canonical form)
  crosses_row <- parsed[parsed$source == "UGA", ]
  expect_equal(crosses_row$male_parent, "Rebel")

  # SHB_Ebrahiem's species-hybrid prose is flagged, not mis-split
  ebrahiem_row <- parsed[parsed$source == "SHB_Ebrahiem", ]
  expect_equal(ebrahiem_row$parse_status, "unparseable_description")
  expect_true(is.na(ebrahiem_row$female_parent))

  # the synonym map itself is exposed separately
  expect_true(any(out$synonym_map$alias == "TH-642" & out$synonym_map$canonical == "Rebel"))
})
