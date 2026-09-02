

test_that("read_al_metadata() extracts MetaData from the bundle unmodified", {
  bundle <- list(
    MetaData = tibble::tibble(
      PlantID = c("P1_G1", "P2_G2"),
      Location = c("EV", "BW"),
      PlantingYear = c(2021, 2022),
      Ecotype = c("SHB", "SHB")
    ),
    OffID = character(0)
  )
  
  out <- read_al_metadata(bundle)
  
  expect_equal(nrow(out), 2)
  expect_true(all(c("PlantID", "Location", "PlantingYear", "Ecotype") %in% names(out)))
  # no filtering/joining -- this is raw pass-through
  expect_equal(out$PlantID, c("P1_G1", "P2_G2"))
})

test_that("read_al_metadata() errors clearly if MetaData is missing from the bundle", {
  expect_error(read_al_metadata(list(OffID = "P1")), "MetaData")
})

test_that("read_al_metadata() wraps a plain ID vector into a tibble", {
  bundle <- list(MetaData = tibble::tibble(), OffID = c("P3_G3", "P7_G7"))
  
  out <- read_al_off_type(bundle)
  
  expect_named(out, "PlantID")
  expect_equal(out$PlantID, c("P3_G3", "P7_G7"))
})

test_that("read_al_off_type() passes a data.frame OffID through as-is", {
  bundle <- list(
    MetaData = tibble::tibble(),
    OffID = tibble::tibble(PlantID = "P3_G3", Reason = "mislabeled")
  )
  
  out <- read_al_off_type(bundle)
  
  expect_named(out, c("PlantID", "Reason"))
  expect_equal(nrow(out), 1)
})

test_that("read_al_off_type() errors clearly if OffID is missing from the bundle", {
  expect_error(read_al_off_type(list(MetaData = tibble::tibble())), "OffID")
})
