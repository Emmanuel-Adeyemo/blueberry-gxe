# 09/11
# basically implementing DataPrep() from legacy scripts R

library(dplyr)
library(purrr)
library(tibble)

AL_METADATA_COLS = c(
  'BlockNo', 'SubplotNo', 'RepNo', 'Genotype', 'Ecotype'
)


AL_GENOTYPE_ALIASES = c(
  "MS771S" = "MS771",
  "BLUERIBBON" = "BlueRibbon" # matched after stripping whitespace
)


#' @description builds the per env plant/plot metadata table used to join a raw al-trait-long df. It filters by loc, remvs plants with no planting year or with an off-type flag
#' @param al_metadata_raw raw al plant/plot metadata ie `al_metadata`. output of `read_al_metadata()`
#' @param al_off_type_raw raw al offtype df eg `al_off_type` output of `read_al_off_type()`
#' @param season single season/year to build the table for
#' @param location single loc eg BW or EV to build table for
#' @output a df 
build_table = function(al_metadata_raw, al_off_type_raw, season, location){

  table_dta = al_metadata_raw %>%
    filter(Location == location) %>% 
    filter(!is.na(PlantingYear)) %>%
    filter(!(PlantID %in% al_off_type_raw$PlantID)) %>%
    mutate(Age = season - PlantingYear) %>%
    filter(Age > 0) %>%
    filter(Ecotype == 'SHB')

  # if(season == 2026){
  #   table_dta = table_dta %>% filter(Ecotype == 'SHB')
  # }

  table_dta %>%
    dplyr::select(
      'PlantID', 'BlockNo', 'SubplotNo', 'RepNo', 'Genotype', 'Ecotype', 'PlantingYear', 'ReplantingYear', 'Age'
    )

}


# interpolate_human_error_dta = function(bloom_long, season, location) {

#   season_tmp = bloom_long %>% filter(season == season, location == location)

#   bloom_clean = season_tmp %>% 
#     arrange(PlantID, obs_date) %>%
#     group_by(PlantID) %>%
#     mutate(
#       running_max = cummax(bloom_percent),
#       bloom_percent = if_else(bloom_percent < running_max, running_max, bloom_percent)
#     ) %>% ungroup() %>%
#     select(-running_max)

# }


# found some issues 09/19 when I compared legacy EV25 processed output with dta that goes into trial QC model
# MS771S vs MS771 naming aliases

#' @description normalizes a raw al `PlantID` string used for matching
#' @param plant_id characer vector of raw `PlantID` values e.g PlantID col from `al_bloom_long`
#' @output character vector of normalized `PlantID` values
.normalize_al_plant_id = function(plant_id) {
  x = gsub(' ', '', plant_id, fixed = TRUE)
  for (i in seq_along(AL_GENOTYPE_ALIASES)) {
    x = gsub(names(AL_GENOTYPE_ALIASES)[i], AL_GENOTYPE_ALIASES[i], x, fixed = TRUE)
  }
  x
}



#' @description entry point for AL cleaning - similar to `DataPrep()` from Sushan's scripts
#' @param trait_long raw long format trait data for one or more envts combinations e.g `al_yield_long`
#' @param al_metadata_raw raw al plant/plot metadata e.g. `al_metadata`, output of `read_al_metadata()`
#' @param al_off_tpe_raw raw al offtype df e.g. `al_off_type`, output of `read_al_metadata()`
#' @output triat long df
clean_al_trait_long = function(trait_long, al_metadata_raw, al_off_type_raw){

  trait_long = trait_long %>% mutate(PlantID = .normalize_al_plant_id(PlantID))
  al_metadata_raw = al_metadata_raw %>% mutate(PlantID = .normalize_al_plant_id(PlantID))
  al_off_type_raw = al_off_type_raw %>% mutate(PlantID = .normalize_al_plant_id(PlantID))

  season_location = trait_long %>% distinct(season, location)

  table_dta = pmap_dfr(season_location, function(season, location) {
    build_table(
      al_metadata_raw, al_off_type_raw, season, location
    ) %>%
      mutate(season = season, location = location)
  })

  trait_long %>% dplyr::select(-any_of(AL_METADATA_COLS)) %>%
    inner_join(table_dta, by = c('PlantID', 'season', 'location')) 
}

