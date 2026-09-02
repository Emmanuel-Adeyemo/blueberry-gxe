

al_trial_fixture <- function() {
  tibble::tribble(
    ~PlantID, ~BlockNo, ~Ecotype, ~`Bloom%_110122`, ~`Bloom%_110822`,
    ~`Berry%_030623`, ~`Frost%_112123`, ~HealthScore_060323, ~PlantHeight_060223,
    ~Notes_012023, ~Yield_040123, ~Yield_040823, ~AverageBerryWeight_040123,
    ~Brix_052623, ~`Brix_06032023`,
    "P1_G1", 1, "SHB", 0, 20, 15, 5, 8.5, 42.0, "looks healthy", NA, 3.2, 1.1, 12.5, 13.0,
    "P2_G2", 1, "SHB", 10, 60, 30, NA, 7.0, 39.5, NA, 1.5, 2.0, NA, 11.0, NA
  )
}

test_that(".al_parse_file_name() parses '<LOC><YY>.ori' and rejects others", {
  expect_equal(.al_parse_file_name("EV23.ori"), list(location = "EV", year = 2023L))
  expect_equal(.al_parse_file_name("BW26.ori"), list(location = "BW", year = 2026L))
  expect_null(.al_parse_file_name("MetaData"))
  expect_null(.al_parse_file_name("OffID"))
})

test_that("extract_trial_loc_and_year() finds trial objects and warns on unrecognized ones", {
  bundle <- list(
    EV23.ori = al_trial_fixture(),
    BW23.ori = al_trial_fixture(),
    MetaData = tibble::tibble(PlantID = "P1_G1"),
    OffID = "P3_G3",
    SomethingUnexpected = 1
  )
  
  expect_warning(
    specs <- extract_trial_loc_and_year(bundle),
    "SomethingUnexpected"
  )
  
  expect_equal(nrow(specs), 2)
  expect_setequal(specs$file_name, c("EV23.ori", "BW23.ori"))
  expect_setequal(specs$location, c("EV", "BW"))
  expect_true(all(specs$year == 2023L))
})

test_that(".al_parse_date() handles mixed 6- and 8-digit date suffixes", {
  out <- .al_parse_date(c("Brix_052623", "Brix_06032023"), "Brix_")
  expect_equal(out, as.Date(c("2023-05-26", "2023-06-03")))
})

test_that(".al_parse_date() warns on an unexpected digit length", {
  expect_warning(
    out <- .al_parse_date("Brix_2023", "Brix_"),
    'al_parse_date(): unexpected date-string length in columns: '
  )
  expect_true(is.na(out))
})

test_that(".al_resolve_id_cols() keeps known identifiers and warns on unrecognized columns", {
  raw <- tibble::tibble(
    PlantID = "P1_G1", BlockNo = 1, Ecotype = "SHB",
    `Bloom%_110122` = 50,
    Date_1 = "2022-11-01" # unrecognized -- neither AL_ID_COLS nor a trait prefix
  )
  
  expect_warning(id_cols <- .al_resolve_id_cols(raw), "Date_1")
  expect_setequal(id_cols, c("PlantID", "BlockNo", "Ecotype"))
  expect_false("Date_1" %in% id_cols)
})

test_that("ingest_all_al_dta() doesn't error when environments have inconsistent stray columns", {
  # Reproduces the reported bug: one environment has an extra column
  # ("Date_1") that the other doesn't. Before the AL_ID_COLS allowlist fix,
  # this caused a bind_rows() type-mismatch error (logical NA-fill vs
  # character) because id_cols was inferred per-environment instead of
  # from a fixed, known set.
  fixture_with_stray_col <- al_trial_fixture()
  fixture_with_stray_col$Date_1 <- c("2022-11-01", "2022-11-08")
  
  bundle <- list(
    EV23.ori = al_trial_fixture(), # no Date_1
    BW23.ori = fixture_with_stray_col # has Date_1
  )
  
  expect_warning(out <- ingest_all_al_dta(bundle), "Date_1")
  
  expect_equal(nrow(out$bloom_long), 8)
  expect_false("Date_1" %in% names(out$bloom_long))
})

test_that("tidy_single_envt() pivots bloom and the new trait families correctly", {
  raw <- al_trial_fixture()
  out <- tidy_single_envt(raw, year = 2023, location = "EV")
  
  expect_equal(nrow(out$bloom_long), 4) # 2 plants x 2 dates, none NA
  expect_setequal(out$bloom_long$obs_date, as.Date(c("2022-11-01", "2022-11-08")))
  
  expect_equal(nrow(out$berry_percent_long), 2)
  expect_true(all(c("year", "location", "obs_date", "berry_percent") %in% names(out$berry_percent_long)))
  
  # P2's Frost% is NA -- only P1's frost reading survives
  expect_equal(nrow(out$frost_percent_long), 1)
  
  expect_equal(nrow(out$health_score_long), 2)
  expect_equal(nrow(out$plant_height_long), 2)
})

test_that("tidy_single_envt() keeps Notes_* separate with a character value column", {
  raw <- al_trial_fixture()
  out <- tidy_single_envt(raw, year = 2023, location = "EV")
  
  notes <- out$notes_long
  expect_true("notes_text" %in% names(notes))
  expect_true(is.character(notes$notes_text))
  # P2's note is NA and dropped -- only P1's survives
  expect_equal(nrow(notes), 1)
  expect_equal(notes$notes_text, "looks healthy")
  
  # Notes must not leak into id_cols of a numeric trait's long table
  expect_false("Notes_012023" %in% names(out$bloom_long))
})

# test_that("tidy_single_envt() combines DESIGN.md-named out-of-scope traits with a trait column", {
#   raw <- al_trial_fixture()
#   out <- tidy_single_envt(raw, year = 2023, location = "EV")
  
#   other <- out$other_traits_long
#   expect_true(all(c("trait", "value") %in% names(other)))
#   expect_true("avg_berry_weight" %in% other$trait)
#   expect_true("brix" %in% other$trait)
#   # Brix has both a 6- and 8-digit-date reading for P1, plus one for P2
#   expect_equal(sum(other$trait == "brix"), 3)
# })

test_that("ingest_al_all() stacks multiple trial objects from a bundle", {
  bundle <- list(
    EV23.ori = al_trial_fixture(),
    BW23.ori = al_trial_fixture(),
    MetaData = tibble::tibble(PlantID = "P1_G1"),
    OffID = "P3_G3"
  )
  
  out <- ingest_all_al_dta(bundle)
  
  expect_equal(nrow(out$bloom_long), 8) # 4 rows/trial x 2 trials
  expect_setequal(out$bloom_long$location, c("EV", "BW"))
})

test_that("read_al_rda() loads every object an .rda defines", {
  EV23.ori <- al_trial_fixture()
  MetaData <- tibble::tibble(PlantID = "P1_G1")
  tmp <- withr::local_tempfile(fileext = ".rda")
  save(EV23.ori, MetaData, file = tmp)
  
  bundle <- read_al_rda(tmp)
  
  expect_setequal(names(bundle), c("EV23.ori", "MetaData"))
  expect_equal(nrow(bundle[["EV23.ori"]]), 2)
})
