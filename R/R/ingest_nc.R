# 10/2

library(dplyr)
library(readxl)
library(tidyr)
library(stringr)


# helper for h2 and f2 NC 2024
#' @description reads an excel sheet 
#' @param path path to the NC raw excel workbook (e.g. `path2024`/`path2025` as passed into `ingest_nc()`
#' @param sheet name of the sheet to read (e.g. "H2 Yield", "F2 Yield")
#' @param header_row row num that holds the actual column headers
#' @output df
.read_excel_at_header_row = function(path, sheet, header_row) {
  raw = suppressMessages(
    read_excel(path, sheet = sheet, col_names = FALSE, col_types = 'text')
  )
  header_vals = as.character(unlist(raw[header_row, ]))
  header_vals[is.na(header_vals)] = paste0('blank_col_', which(is.na(header_vals)))
  header_vals = make.unique(header_vals, sep = "...")

  df = as.data.frame(raw[-seq_len(header_row), ], stringsAsFactors = FALSE)
  names(df) = header_vals
  df
}



# YIELD
#' @description reads the raw NC 2025 yield summary sheet
#' @param path path to the NC 2025 workbook. In the pipeline it is `nc_path25` file target passed in by `ingest_nc()`
#' @output df
read_nc_yld_2025_raw = function(path){
  read_excel(path, sheet = 'YieldSum')
}


#' @description reads the raw NC 2024 "H2" yield trial sheet 
#' @param path path to the NC 2024 workbook. In the pipeline this is the `nc_path24` file passed in by `ingest_nc()`.
#' @output df
read_nc_yld_2024_h2_raw = function(path){
  .read_excel_at_header_row(path, sheet = 'H2 Yield', header_row = 2)
}


#' @description reads the raw NC 2024 "F2" yield trial sheet 
#' @param path path to the NC workbook whose `F2 Yield` sheet should be read. In the pipeline this is called with `path2024` inside `ingest_nc()`.
#' @output df
read_nc_yld_2025_f2_raw = function(path){
  .read_excel_at_header_row(path, sheet = 'F2 Yield', header_row = 2)
}

# META
#' @description reads the raw NC planting-dates/metadata sheet.
#' @param path path to an NC workbook with a `planting dates` sheet. In the pipeline this is called with `path2025` inside `ingest_nc()`.
#' @output The raw tibble from the `planting dates` sheet
read_nc_metadata = function(path){
  read_excel(path, sheet = 'planting dates')
}


#' @description selects and renames the columns needed from the raw NC planting-dates/metadata sheet.
#' @param raw output of `read_nc_metadata()`.
#' @output df
process_meta = function(raw){
  raw %>% dplyr::select(plot, genotype = name, planting_year)
}


# normalize raw yld into common schema
# 2025
#' @description reshapes the raw NC 2025 yield sheet into the pipeline's common schema for yld
#' @param raw output of `read_nc_yld_2025_raw()`
#' @output df
process_nc_yld_2025 = function(raw){
  raw %>% transmute(
    season = 2025L,
    location = 'NC',
    site = as.character(field),
    plot = as.character(plot),
    subplot = as.character(.data$`Set entry`),
    genotype = as.character(name),
    y = as.numeric(.data$`Total Yield/bush (g)`)

  )
}


# 2024
#' @description reshapes the raw NC 2024 "H2" yield sheet into the pipeline's common yld schema
#' @output df
process_nc_yld_2024_h2 = function(raw) {
  raw %>%
    transmute(
      season = 2024L,
      location = 'NC',
      site = as.character(Field),
      plot = as.character(Plot),
      subplot = as.character(Bush),
      genotype = as.character(Name),
      y = as.numeric(.data$`Total g/Plant`),
      # kept
      yld_note = str_trim(Notes)
    )
}


#' @description reshapes the raw NC 2024 "F2" yield sheet into the pipeline's common yld schema.
#' @param raw output of `read_nc_yld_2025_f2_raw()` 
#' @output df
process_nc_yld_2024_f2 = function(raw) {
  raw %>%
    transmute(
      season = 2024L,
      location = 'NC',
      site = as.character(Field),
      plot = as.character(Plot),
      subplot = as.character(Bush),
      genotype = as.character(Name),
      y = as.numeric(.data$`Total g/Plant`)
    )
}


# Bloom
#' @description reads the raw NC 2025 
#' @param path path to the NC 2025 workbook. In the pipeline this is the `nc_path25` 
#' @output df
read_nc_bloom_2025_raw = function(path) {
  .read_excel_at_header_row(path, sheet = 'Flowering dates', header_row = 2)
}


#' @description reads the raw NC 2024 
#' @param path path to the NC 2024 workbook. In the pipeline this is the `nc_path24` 
#' @output df
read_nc_bloom_2024_raw = function(path) {
  raw = read_excel(path, sheet = '24 flowering w calc')
  names(raw) = str_trim(names(raw))
  raw
}


# build plot - site columns from yield data
#' @description derives a (plot : site) lookup table from the already processed 2025 yld dta
#' @param yield_2025 output of `process_nc_yld_2025()`
#' @output df
derive_nc_plot_site_map_2025 = function(yield_2025) {
  yield_2025 %>% distinct(plot, site)
}



#' @description reshapes the raw NC 2025 bloom sheet into the pipeline's common bloom-timing schema
#' @param raw output of `read_nc_bloom_2025_raw()`
#' @param plot_site_map output of `derive_nc_plot_site_map_2025()`
#' @output df
process_nc_bloom_2025 = function(raw, plot_site_map) {
  .excel_serial_to_date = function(x) as.Date(as.numeric(x), origin = '1899-12-30')

  raw %>%
    transmute(
      season = 2025L,
      location = 'NC',
      plot = as.character(.data$`ID...1`),
      subplot = as.character(.data$`Set entry`),
      genotype = as.character(.data$Name),
      bloom5 = .excel_serial_to_date(.data$`5% bloom date`),
      bloom50 = .excel_serial_to_date(.data$`50% bloom date`),
      bloom100 = .excel_serial_to_date(.data$`100% bloom date`)
    ) %>%
    left_join(plot_site_map, by = 'plot') %>%
    relocate('site', .after = 'location')
}


#' @description reshapes the raw NC 2024 bloom sheet into the pipeline's common bloom-timing schema
#' @param raw output of `read_nc_bloom_2024_raw()`
#' @output df
process_nc_bloom_2024 = function(raw) {
  raw %>%
    transmute(
      season = 2024L,
      location = 'NC',
      site = as.character(Loc),
      plot = as.character(Plot),
      subplot = as.character(Bushes),
      genotype = as.character(Name),
      bloom5 = as.Date(.data$`5% date`),
      bloom50 = as.Date(.data$`50% date`),
      bloom100 = as.Date(.data$`100% date`)
    )
}


# Entry point
#' @description entry point for NC ingestion entry point
#' @param path2025 path to the combined NC 2025 workbook, in the pipeline this is the `nc_path25` file target
#' @param path2024 path to the combined NC 2024 workbook (H2/F2 yield, bloom). In the pipeline this is the `nc_path24` file target
#' @output list
ingest_nc = function(path2025, path2024){

  yld_2025 = process_nc_yld_2025(read_nc_yld_2025_raw(path2025))
  yld_2024_h2 = process_nc_yld_2024_h2(read_nc_yld_2024_h2_raw(path2024))
  yld_2024_f2 = process_nc_yld_2024_f2(read_nc_yld_2025_f2_raw(path2024))

  meta = process_meta(read_nc_metadata(path2025))

  plot_site_map_2025 = derive_nc_plot_site_map_2025(yld_2025) # builds plot - site for bloom data
  bloom_2025 = process_nc_bloom_2025(read_nc_bloom_2025_raw(path2025), plot_site_map_2025)
  bloom_2024 = process_nc_bloom_2024(read_nc_bloom_2024_raw(path2024))

  list(
    yield = bind_rows(
      yld_2025, yld_2024_h2 %>% dplyr::select(-'yld_note'),
      yld_2024_f2
    ),
    bloom = bind_rows(bloom_2025, bloom_2024),
    meta = meta
  )
}