# 10/16

library(dplyr)
library(readxl)
library(stringr)


# Brix / TA (acidity) -- 2025, sheet 'bapH'
#' @description reads the raw NC 2025 Brix/TA sheet
#' @param path path to the 2025 NCSU workbook. In the pipeline this is called with `path2025` inside `ingest_nc_quality()`.
#' @output the raw tibble from the `bapH` sheet
read_nc_bapH_2025_raw = function(path) {
  read_excel(path, sheet = 'bapH')
}


#' @description reshapes the raw NC 2025 Brix/TA sheet into the pipeline's common quality-trait schema.
#' @param raw output of `read_nc_bapH_2025_raw()`.
#' @output df
process_nc_brix_ta_2025 = function(raw) {
  raw %>%
    transmute(
      season = 2025L,
      location = 'NC',
      site = as.character(Field),
      plot = as.character(Plot),
      subplot = as.character(Bush),
      genotype = as.character(Name),
      brix = as.numeric(Brix),
      ta = as.numeric(Acidity)
    )
}


# avg berry weight - 2025, sheet 'berry', no bush/plant column available.
#' @description reads the raw NC 2025 avg berry weight sheet
#' @param path path to the 2025 NCSU workbook. In the pipeline this is called with `path2025` inside `ingest_nc_quality()`.
#' @output the raw tibble from the `berry` sheet
read_nc_berry_2025_raw = function(path) {
  read_excel(path, sheet = 'berry')
}


#' @description reshapes the raw NC 2025 berry weight sheet into the pipeline's common quality-trait schema.
#' @param raw output of `read_nc_berry_2025_raw()`.
#' @output df
process_nc_avg_berry_weight_2025 = function(raw) {
  raw %>%
    transmute(
      season = 2025L,
      location = 'NC',
      site = as.character(field),
      plot = as.character(plot),
      genotype = as.character(name),
      avg_berry_weight = as.numeric(.data$`Wght/Berry (g)`)
    )
}


# firmness -- 2025, sheet 'TANG'. no bush/plant column available
#' @description reads the raw NC 2025 tang sheet
#' @param path path to the 2025 NCSU workbook. In the pipeline this is called with `path2025` inside `ingest_nc_quality()`.
#' @output the raw tibble from the `TANG` sheet
read_nc_tang_2025_raw = function(path) {
  read_excel(path, sheet = 'TANG')
}


#' @description reshapes the raw NC 2025 firmness sheet into the pipeline's common quality-trait schema
#' @param raw Output of `read_nc_tang_2025_raw()`
#' @output df
process_nc_firmness_2025 = function(raw) {
  raw %>%
    transmute(
      season = 2025L,
      location = 'NC',
      site = as.character(field),
      plot = as.character(plot),
      genotype = as.character(Name),
      firmness = as.numeric(.data$`Caclulated Firmtech Firmness (g/mm)`)
    )
}


# all four quality traits, 2024, sheet is 24 Harvest Data
#' @description reads the raw NC 2024 harvest-data sheet, which holds all four quality traits (berry weight, firmness, Brix, TA) in one place.
#' @param path path to the 2024 NCSU workbook. in the pipeline this is called with `path2024` inside `ingest_nc_quality()`
#' @output df
read_nc_harvest_2024_raw = function(path) {
  read_excel(path, sheet = '24 Harvest Data ')
}

#' @description reshapes the raw NC 2024 harvest-data sheet into the pipeline's common quality-trait schema
#' @param raw output of `read_nc_harvest_2024_raw()`
#' @output df
process_nc_quality_2024 = function(raw) {
  raw %>%
    transmute(
      season = 2024L,
      location = 'NC',
      site = as.character(Field),
      plot = as.character(Orderid),
      subplot = toupper(str_trim(as.character(Bush))),
      genotype = as.character(Name),
      avg_berry_weight = rowMeans(
        cbind(as.numeric(.data$`g/berry 1`), as.numeric(.data$`g/berry 2`)), na.rm = TRUE),
      firmness = as.numeric(.data$`FIRMNESS g`),
      brix = as.numeric(Brix),
      ta = as.numeric(Acid)
    ) %>%
    mutate(avg_berry_weight = ifelse(is.nan(avg_berry_weight), NA_real_, avg_berry_weight))
}


#' @description entry point. ingest NC quality traits (avg berry weight, firmness, Brix, TA) across the two years
#' @param path2025 - path to 2025 NCSU workbook 
#' @param path2024 - path to the 2024 NCSU workbook
#' @output list(avg_berry_weight, firmness, brix, ta)
ingest_nc_quality = function(path2025, path2024) {

  bapH_2025 = process_nc_brix_ta_2025(read_nc_bapH_2025_raw(path2025))
  berry_2025 = process_nc_avg_berry_weight_2025(read_nc_berry_2025_raw(path2025))
  tang_2025 = process_nc_firmness_2025(read_nc_tang_2025_raw(path2025))
  harvest_2024 = process_nc_quality_2024(read_nc_harvest_2024_raw(path2024))

  list(
    avg_berry_weight = bind_rows(
      berry_2025 %>% dplyr::select(season, location, site, plot, genotype, avg_berry_weight),
      harvest_2024 %>% dplyr::select(season, location, site, plot, subplot, genotype, avg_berry_weight)
    ),
    firmness = bind_rows(
      tang_2025 %>% dplyr::select(season, location, site, plot, genotype, firmness),
      harvest_2024 %>% dplyr::select(season, location, site, plot, subplot, genotype, firmness)
    ),
    brix = bind_rows(
      bapH_2025 %>% dplyr::select(season, location, site, plot, subplot, genotype, brix),
      harvest_2024 %>% dplyr::select(season, location, site, plot, subplot, genotype, brix)
    ),
    ta = bind_rows(
      bapH_2025 %>% dplyr::select(season, location, site, plot, subplot, genotype, ta),
      harvest_2024 %>% dplyr::select(season, location, site, plot, subplot, genotype, ta)
    )
  )
}