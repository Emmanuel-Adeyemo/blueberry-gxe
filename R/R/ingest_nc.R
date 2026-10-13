# 10/2

library(dplyr)
library(readxl)
library(tidyr)
library(stringr)


# helper for h2 and f2 NC 2024
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

read_nc_yld_2025_raw = function(path){
  read_excel(path, sheet = 'YieldSum')
}


read_nc_yld_2024_h2_raw = function(path){
  .read_excel_at_header_row(path, sheet = 'H2 Yield', header_row = 2)
}


read_nc_yld_2025_f2_raw = function(path){
  .read_excel_at_header_row(path, sheet = 'F2 Yield', header_row = 2)
}


# normalize raw yld into common schema
# 2025
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
process_nc_yld_2024_h2 = function(raw) {
  raw %>%
    transmute(
      season = 2024L,
      location = 'NC',
      site = as.character(Field),
      plot = as.character(Plot),
      subplot = as.character(Plant),
      genotype = as.character(Name),
      y = as.numeric(.data$`Total g/Plant`),
      # kept
      yld_note = str_trim(Notes)
    )
}


process_nc_yld_2024_f2 = function(raw) {
  raw %>%
    transmute(
      season = 2024L,
      location = 'NC',
      site = as.character(Field),
      plot = as.character(Plot),
      subplot = as.character(Plant),
      genotype = as.character(Name),
      y = as.numeric(.data$`Total g/Plant`)
    )
}


# Bloom

read_nc_bloom_2025_raw = function(path) {
  .read_excel_at_header_row(path, sheet = 'Flowering dates', header_row = 2)
}


read_nc_bloom_2024_raw = function(path) {
  raw = read_excel(path, sheet = '24 flowering w calc')
  names(raw) = str_trim(names(raw))
  raw
}


# build plot - site columns from yield data
derive_nc_plot_site_map_2025 = function(yield_2025) {
  yield_2025 %>% distinct(plot, site)
}



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
ingest_nc = function(path2025, path2024){

  yld_2025 = process_nc_yld_2025(read_nc_yld_2025_raw(path2025))
  yld_2024_h2 = process_nc_yld_2024_h2(read_nc_yld_2024_h2_raw(path2024))
  yld_2024_f2 = process_nc_yld_2024_f2(read_nc_yld_2025_f2_raw(path2024))

  plot_site_map_2025 = derive_nc_plot_site_map_2025(yld_2025) # builds plot - site for bloom data
  bloom_2025 = process_nc_bloom_2025(read_nc_bloom_2025_raw(path2025), plot_site_map_2025)
  bloom_2024 = process_nc_bloom_2024(read_nc_bloom_2024_raw(path2024))

  list(
    yield = bind_rows(
      yld_2025, yld_2024_h2 %>% dplyr::select(-'yld_note'),
      yld_2024_f2
    ),
    bloom = bind_rows(bloom_2025, bloom_2024)
  )
}