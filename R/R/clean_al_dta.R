# 09/11
# basically implementing DataPrep() from legacy scripts R

library(dplyr)
library(purrr)
library(tibble)

AL_METADATA_COLS = c(
  'BlockNo', 'SubplotNo', 'RepNo', 'Genotype', 'Ecotype'
)


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
    select(
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



clean_al_trait_long = function(trait_long, al_metadata_raw, al_off_type_raw){
  season_location = trait_long %>% distinct(season, location)

  table_dta = pmap_dfr(season_location, function(season, location) {
    build_table(
      al_metadata_raw, al_off_type_raw, season, location
    ) %>%
      mutate(season = season, location = location)
  })

  trait_long %>% select(-any_of(AL_METADATA_COLS)) %>%
    inner_join(table_dta, by = c('PlantID', 'season', 'location')) 
}

