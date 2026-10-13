# 10/02

library(dplyr)
library(readxl)
library(tidyr)
library(stringr)


.week_since_jan1_to_date = function(year, week) {
  as.Date(paste0(year, '-01-01')) + (week - 1) * 7
}



# map site details
fl_site_detail_2025_2024 = function(site, block_or_field) {
  is_waldo = !is.na(site) & site == 'Waldo' # %in% not used here since it would silently treat NA site as 'not Waldo'
  year_extracted = str_extract(block_or_field, '(19|20)[0-9]{2}')
  site_detail = dplyr::case_when(
    is.na(site) ~ NA_character_,
    is_waldo & is.na(year_extracted) ~ NA_character_, 
    
    is_waldo ~ paste('Waldo', year_extracted),
    TRUE ~ paste(site, block_or_field) # Citra (and any other site):
  )

  if (any(is_waldo & is.na(year_extracted))) {
    bad = unique(block_or_field[is_waldo & is.na(year_extracted)])
    warning(
      'fl_site_detail_2025_2024(): could not extract a 4-digit year from ',
      'Waldo block/field value(s): ', paste(bad, collapse = ', '),' -- returns NA instead of just  guessing.'
    )
  }
  site_detail
}



# Yield

read_fl_yld_2025_raw = function(path){
  read_excel(path, sheet = 'yield_fq')
}


read_fl_yld_2024_raw = function(path){
  read_excel(path, sheet = 'Sheet 1')
}


process_fl_yld_2025 = function(raw) {
  raw %>%
    filter(!is.na(genotype)) %>%
    group_by(genotype, site, block, bush_plant_number) %>%
    summarize(y = sum(mass, na.rm = TRUE), .groups = 'drop') %>%
    transmute(
      season = 2025L,
      location = 'FL',
      site_detail = fl_site_detail_2025_2024(site, block),
      subplot = as.character(bush_plant_number),
      genotype = as.character(genotype),
      y = y
    )
}



process_fl_yld_2024 = function(raw) {
  raw %>%
    # remv stray note from Mr Paul
    filter(!is.na(genotype), !is.na(site)) %>%
    group_by(genotype, site, field, bush) %>%
    summarise(y = sum(.data$`mass (g)`, na.rm = TRUE), .groups = 'drop') %>%
    transmute(
      season = 2024L,
      location = 'FL',
      site_detail = fl_site_detail_2025_2024(site, field),
      subplot = as.character(bush),
      genotype = as.character(genotype),
      y = y
    )
}


# Bloom

read_fl_bloom_2025_raw = function(path) {
  read_excel(path, sheet = 'info_bloom')
}


read_fl_bloom_2024_raw = function(path) {
  read_excel(path, sheet = 'Sheet1')
}


process_fl_bloom_2025 = function(raw) {
  raw %>%
    transmute(
      season = 2025L,
      location = 'FL',
      site_detail = as.character(Location),
      genotype = as.character(Genotype),
      bloom5 = as.Date(.data$`5% bloom`),
      bloom50 = as.Date(.data$`50% bloom`),
      bloom100 = as.Date(.data$`Full Bloom`)
    )
}


process_fl_bloom_2024 = function(raw) {
  raw %>%
    filter(!is.na(Genotype)) %>%
    transmute(
      season = 2024L,
      location = 'FL',
      site_detail = paste(Location, Block),
      genotype = as.character(Genotype),
      bloom5 = .week_since_jan1_to_date(2024L, .data$`F5 - 5% Bloom`),
      bloom50 = .week_since_jan1_to_date(2024L, .data$`F50 - 50% Bloom`),
      bloom100 = .week_since_jan1_to_date(2024L, .data$`100% Bloom`)
    )
}



# entry point
ingest_fl = function(path_2025, path_2024_yld, path_2024_bloom) {
  yld_2025 = process_fl_yld_2025(read_fl_yld_2025_raw(path_2025))
  yld_2024 = process_fl_yld_2024(read_fl_yld_2024_raw(path_2024_yld))

  bloom_2025 = process_fl_bloom_2025(read_fl_bloom_2025_raw(path_2025))
  bloom_2024 = process_fl_bloom_2024(read_fl_bloom_2024_raw(path_2024_bloom))

  list(
    yield = bind_rows(yld_2025, yld_2024),
    bloom = bind_rows(bloom_2025, bloom_2024)
  )
}