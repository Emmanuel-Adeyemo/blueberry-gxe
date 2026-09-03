# 08/31/2026

library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(tibble)


read_uf_pedigree = function(path){
  
  env = new.env()
  loaded_names = load(path, envir = env)
  
  # if (!'ped_full' %in% loaded_names){
  #   stop(
  #     'read_uf_pedigree(): expected file name "ped_full" in ', path, 'but found: ', paste(loaded_names, collapse = ', ')
  #   )
  # } 
  
  raw = env[[loaded_names]]
  
  raw %>%
    as_tibble() %>%
    transmute(
      source = 'UF',
      line_id = as.character(.data$Clone),
      female_parent = na_if(as.character(.data$Female), '0'),
      male_parent = na_if(as.character(.data$Male), '0')
                            
    ) %>%
    filter(!(is.na(.data$line_id) | .data$line_id == '0'))
}




read_ncsu_pedigree = function(path) {
  raw = read_excel(path, sheet = "Data")
  
  raw %>%
    transmute(
      source = "NCSU",
      line_id = as.character(.data$CommercialName),
      female_parent = na_if(str_trim(as.character(.data$FemaleParent)), "unk"),
      male_parent = na_if(str_trim(as.character(.data$MaleParent)), "unk"),
      patent_date = .data[["patent date"]],
      expiration_date = .data[["expiration date"]],
      owner = .data$owner
    ) %>%
    filter(!is.na(.data$line_id))
}



read_uga_pedigree_crosses = function(path) {
  raw = read_excel(
    path,
    sheet = "Sheet1",
    col_names = c("root_label", "formula_cache", "line_id", "pedigree_raw"),
    col_types = "text"
  )
  
  raw %>%
    transmute(
      source = "UGA",
      line_id = str_trim(.data$line_id),
      pedigree_raw = str_trim(.data$pedigree_raw)
    ) %>%
    filter(!is.na(.data$line_id) & .data$line_id != "")
}





read_uga_pedigree_ancestors = function(path) {
  # Read using the sheet's own real header row rather than assuming a fixed
  # column count -- the real file carries a few trailing blank columns
  # (up to E), which this approach ignores naturally rather than requiring
  # an exact positional match.
  raw = read_excel(path, sheet = "Sheet2", col_types = "text")
  
  raw = raw %>%
    rename(line_id = "Line ID", pedigree_raw = "Pedigree") %>%
    select("line_id", "pedigree_raw")
  
  is_blank_row = is.na(raw$line_id) & is.na(raw$pedigree_raw)
  
  raw %>%
    mutate(
      source = "UGA_ancestors",
      block_id = cumsum(is_blank_row) + 1L
    ) %>%
    filter(!is_blank_row) %>%
    mutate(
      line_id = str_trim(.data$line_id),
      pedigree_raw = str_trim(.data$pedigree_raw)
    ) %>%
    select("source", "block_id", "line_id", "pedigree_raw")
}




read_ms_usda_pedigree = function(path){
  
  raw = read_csv(path, show_col_types = FALSE)
  
  raw %>%
    transmute(
      source = 'MS_USDA',
      line_id = str_trim(as.character(.data$Genotype)),
      pedigree_raw = str_trim(as.character(.data$Pedigree))
    ) %>% 
    filter(!is.na(.data$line_id) & .data$line_id != '')
}



read_shb_ebrahiem_pedigree <- function(path) {
  raw = read_excel(path, sheet = "Sheet1")
  
  raw %>%
    transmute(
      source = "SHB_Ebrahiem",
      line_id = str_trim(as.character(.data$Accession)),
      pedigree_raw = str_trim(as.character(.data$Pedigree)),
      observation_raw = str_trim(as.character(.data$Observation))
    ) %>%
    filter(!is.na(.data$line_id) & .data$line_id != "")
}
