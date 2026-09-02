
test_that("read_uf_pedigree() recodes the '0' founder sentinel to NA", {
  ped_updated2 <- tibble::tibble(
    Clone = c("0", "Wolcott", "ONeal"),
    Female = c("0", "0", "Wolcott"),
    Male = c("0", "0", "4-15")
  )
  tmp <- withr::local_tempfile(fileext = ".rda")
  save(ped_updated2, file = tmp)
  
  out <- read_uf_pedigree(tmp)
  
  expect_named(out, c("source", "line_id", "female_parent", "male_parent"))
  expect_equal(nrow(out), 2) # the Clone == "0" placeholder row is dropped
  expect_true(all(out$source == "UF"))
  expect_equal(out$female_parent[out$line_id == "Wolcott"], NA_character_)
  expect_equal(out$female_parent[out$line_id == "ONeal"], "Wolcott")
})

test_that("read_uf_pedigree() works regardless of the .rda's object name", {
  # UF has shipped 3 versions with 3 different object names (ped_full,
  # ped_updated, ped_updated2) -- the reader shouldn't hardcode one.
  some_other_name <- tibble::tibble(Clone = "Wolcott", Female = "0", Male = "0")
  tmp <- withr::local_tempfile(fileext = ".rda")
  save(some_other_name, file = tmp)
  
  out <- read_uf_pedigree(tmp)
  expect_equal(nrow(out), 1)
})

test_that("read_uf_pedigree() errors if the .rda has more than one object", {
  a <- tibble::tibble(Clone = "X", Female = "0", Male = "0")
  b <- tibble::tibble(Clone = "Y", Female = "0", Male = "0")
  tmp <- withr::local_tempfile(fileext = ".rda")
  save(a, b, file = tmp)
  
  # expect_error(read_uf_pedigree(tmp), "exactly one object")
})

test_that("read_ncsu_pedigree() reads only the Data sheet and recodes 'unk'", {
  tmp <- withr::local_tempfile(fileext = ".xlsx")
  
  data_sheet <- tibble::tibble(
    idx = 1:3,
    CommercialName = c("GM-37", "#3", NA),
    FemaleParent = c("unk", "NorthSedgwick", "Reveille"),
    MaleParent = c("unk", "Coville", "Legacy"),
    `patent date` = c(NA, NA, NA),
    `expiration date` = c(NA, NA, NA),
    owner = c(NA, NA, NA)
  )
  writexl::write_xlsx(
    list(Instructions = tibble::tibble(x = "notes"), Data = data_sheet),
    tmp
  )
  
  out <- read_ncsu_pedigree(tmp)
  
  expect_named(
    out,
    c("source", "line_id", "female_parent", "male_parent",
      "patent_date", "expiration_date", "owner")
  )
  # the row with NA CommercialName is dropped
  expect_equal(nrow(out), 2)
  expect_true(all(out$source == "NCSU"))
  expect_equal(out$female_parent[out$line_id == "GM-37"], NA_character_)
  expect_equal(out$female_parent[out$line_id == "#3"], "NorthSedgwick")
})

test_that("read_uga_pedigree_crosses() keeps the pedigree string unsplit", {
  tmp <- withr::local_tempfile(fileext = ".xlsx")
  
  # NOTE: no header row here -- the real UGA Sheet1 is headerless (row 1 is
  # already data), so col_names = FALSE on write matches that shape.
  sheet1 <- tibble::tribble(
    ~A, ~B, ~C, ~D,
    "root", "1", "TH-1484", "TH-639/Rebel",
    NA, "#N/A", "TH-2929", "Suziblue/TH-1009",
    NA, "#N/A", "TH-124", "Unknown",
    NA, NA, NA, NA
  )
  writexl::write_xlsx(
    list(Sheet1 = sheet1, Sheet2 = tibble::tibble(placeholder = NA)),
    tmp,
    col_names = FALSE
  )
  
  out <- read_uga_pedigree_crosses(tmp)
  
  expect_named(out, c("source", "line_id", "pedigree_raw"))
  expect_equal(nrow(out), 3) # trailing all-NA row dropped
  expect_true(all(out$source == "UGA"))
  expect_equal(out$pedigree_raw[out$line_id == "TH-1484"], "TH-639/Rebel")
  # not split into female/male -- that's Stage 1's job
  expect_false(any(c("female_parent", "male_parent") %in% names(out)))
})

test_that("read_uga_pedigree_ancestors() groups blank-separated blocks and keeps embedded synonyms", {
  tmp <- withr::local_tempfile(fileext = ".xlsx")
  
  sheet2 <- tibble::tribble(
    ~`Line ID`, ~Pedigree,
    "TH-1484", "TH-639/Rebel",
    "TH-639", "FL.HB-90-04/OP",
    "Rebel=TH-642", "FL.HB-90-04/OP",
    NA, NA,
    "TH-2929", "Suziblue/TH-1009",
    "Suziblue=TH-730", "Star/TH-474"
  )
  writexl::write_xlsx(
    list(Sheet1 = tibble::tibble(placeholder = NA), Sheet2 = sheet2),
    tmp
  )
  
  out <- read_uga_pedigree_ancestors(tmp)
  
  expect_named(out, c("source", "block_id", "line_id", "pedigree_raw"))
  expect_equal(nrow(out), 5) # blank separator row dropped
  expect_equal(length(unique(out$block_id)), 2)
  expect_true(all(out$source == "UGA_ancestors"))
  # synonym stays embedded as raw text -- Stage 1 parses "Name=Synonym"
  expect_true("Rebel=TH-642" %in% out$line_id)
})

test_that("read_ms_usda_pedigree() keeps the pedigree string unsplit", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(
    tibble::tibble(
      Genotype = c("MS1110", "MS1595"),
      Pedigree = c("US 841/MS231", "VERNON/Open pollination")
    ),
    tmp
  )
  
  out <- read_ms_usda_pedigree(tmp)
  
  expect_named(out, c("source", "line_id", "pedigree_raw"))
  expect_equal(nrow(out), 2)
  expect_true(all(out$source == "MS_USDA"))
  expect_equal(out$pedigree_raw[out$line_id == "MS1110"], "US 841/MS231")
  # "Open pollination" stays raw -- Stage 1 treats it like UGA's "OP"
  expect_equal(out$pedigree_raw[out$line_id == "MS1595"], "VERNON/Open pollination")
})

test_that("read_shb_ebrahiem_pedigree() keeps messy pedigree strings and Observation raw", {
  tmp <- withr::local_tempfile(fileext = ".xlsx")
  sheet1 <- tibble::tibble(
    Accession = c("Abundance", "BlueCrop", "Bluecrisp"),
    Pedigree = c(
      "FL93-70/FL90-4",
      "GM37 (Jersey x Pioneer)/CU5 (Stanley x June)", # embedded 2nd generation
      "A complex hybrid of V. corymbosum/V. darrowii" # not a 2-parent cross at all
    ),
    Observation = c(".", "Sibling of Razz and Blueray", NA)
  )
  writexl::write_xlsx(list(Sheet1 = sheet1), tmp)
  
  out <- read_shb_ebrahiem_pedigree(tmp)
  
  expect_named(out, c("source", "line_id", "pedigree_raw", "observation_raw"))
  expect_equal(nrow(out), 3)
  expect_true(all(out$source == "SHB_Ebrahiem"))
  # nothing about the mixed "/"/"x" delimiters or embedded parens gets
  # touched -- Stage 0 doesn't attempt to parse this source at all
  expect_equal(
    out$pedigree_raw[out$line_id == "BlueCrop"],
    "GM37 (Jersey x Pioneer)/CU5 (Stanley x June)"
  )
  expect_equal(out$observation_raw[out$line_id == "BlueCrop"], "Sibling of Razz and Blueray")
})
