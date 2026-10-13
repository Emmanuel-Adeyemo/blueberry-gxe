# 08/28/2026
# blueberry/R/R
# legacy R scripts already processed AL field data. 
# Rerun to try rxn norm for bloom data instead of the 
# interpolation used in the legacy script

# 
# EV23.ori <- read_csv('Data/EV23_071526.csv')
# EV24.ori <- read_csv('Data/EV24_071526.csv')
# EV25.ori <- read_csv('Data/EV25_071526.csv')
# EV26.ori <- read_csv('Data/EV26_071526.csv')
# 
# BW23.ori <- read_csv('Data/BW23_071526.csv')
# BW24.ori <- read_csv('Data/BW24_071526.csv')
# BW25.ori <- read_csv('Data/BW25_071526.csv')
# BW26.ori <- read_csv('Data/BW26_071526.csv')
# 
# OffID <- read_csv('Data/OffType.csv')$PlantID
# MetaData <- read_csv("'Data/MetaData.csv')
# 
# save(EV23.ori, EV24.ori, EV25.ori, EV26.ori,
#      BW23.ori, BW24.ori, BW25.ori, BW26.ori, # typo here. Changed EV26.ori (appeared twice) to BW26.ori
#      OffID, MetaData, file = "OriginalData.rda")
# 
# 
# 
# # AL data is compiled into a rda format
# load('data/raw/al/OriginalData.rda')
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)


AL_TRAIT_PREFIXES = c(
  bloom_percent = 'Bloom%_', 
  yield = 'Yield_'
  
)



AL_ID_COLS = c(
  'PlantID', 'BlockNo', 'SubplotNo', 'RepNo', 'Genotype', 'Ecotype',
  "HealthSummary", 'PlantHeightSummary', 'Harvested"'
)



#' @param path - path to AL rda file
#' @output named list for every file in the rda file
read_al_rda = function(path){
  
  env = new.env()
  load(path, envir = env)
  as.list(env)
}

# a = read_al_rda('OriginalData.rda')



#' @description parse trial file name into loc/season. Should return NULL for MetaData and OffID files
#' @param file_name - a single file name eg EV23.ori
#' @output list(loc,season)

.al_parse_file_name = function(file_name){
  m = stringr::str_match(file_name, '^([A-Za-z]+)([0-9]{2})\\.ori$')
  if (is.na(m[1, 1])){
    return(NULL) # for metadata and offtype
  }
  list(location = m[1, 2], season = 2000L + as.integer(m[1, 3]))
}



#' @description identifies the trial files within the loaded rda file
#' @param rda_file - combination of all csv trials in AL
#' @output tibble(file_name, location, season)
#' 
extract_trial_loc_and_season = function(rda_file){
  
  file_names = names(rda_file)
  parsed = lapply(file_names, .al_parse_file_name)
  is_trial = !vapply(parsed, is.null, logical(1)) # sanity checks for trials, should be 8 and 2 nulls
  # 
  unexpected = setdiff(file_names[!is_trial], c('MetaData', 'OffID'))
  if (length(unexpected) > 0) {
    warning('unexpected file extracted from AL rda file', paste(unexpected, collapse = ', '))
  }
  
  tibble::tibble(
    file_name = file_names[is_trial],
    location = vapply(parsed[is_trial], `[[`, character(1), 'location'),
    season = vapply(parsed[is_trial], `[[`, integer(1), 'season')
  )
  
}

# rm(file_names, parsed, is_trial)
# file_names = a
# parsed = lapply(file_names, .al_parse_file_name)
# is_trial = !vapply(parsed, is.null, logical(1)) 


#' @description parse colname into date eg Yield_MM_DD_YY
#' @param col_name - character vector of column names
#' @param prefix - the trait prefix to strip eg Bloom%_ 

.al_parse_date = function(col_names, prefix) {
  raw_date = str_sub(col_names, nchar(prefix) + 1)
  # in case of / or - in date
  raw_date = str_remove_all(raw_date, '[^0-9]')
  n = nchar(raw_date)

  out = as.Date(rep(NA_character_, length(raw_date)))
  out[n == 6] = as.Date(raw_date[n == 6], format = '%m%d%y')
  out[n == 8] = as.Date(raw_date[n == 8], format = '%m%d%Y')

  # unexpected = !(n %in% c(6, 8))
  # if (any(unexpected)) {
  #   warning(
  #     "al_parse_obs_date(): unexpected date-string length in columns: ",
  #     paste(col_names[unexpected], collapse = ", ")
  #   )
  # }

  out
}




#' @param raw - original wide data ie single envt
#' @param output - character vector of id column names present in raw
.al_resolve_id_cols = function(raw) {

  # quickfix for columnn name misspelling in bw25 and bw26
  if ('SubloptNo' %in% names(raw)){
    raw = raw %>% rename(SubplotNo = SubloptNo)
  }
  
  all_trait_cols = unlist(lapply(
    AL_TRAIT_PREFIXES,
    function(p) names(raw) [stringr::str_starts(names(raw), stringr::fixed(p))]
  ))

  id_cols = intersect(AL_ID_COLS, names(raw))
  
  unexpected = setdiff(names(raw), c(all_trait_cols, id_cols))
  # if(length(unexpected) > 0){
  #   warning(
  #     '.al_resolve_id_cols(): dropping unrecognised column. Check known identifiers AL_TRAIT_PREFIXES and ALL_ID_COLS for details. ',
  #     paste(unexpected, collapse = ', ')
  #   )
  # }
  id_cols
}



#' @param raw - original wide data ie single envt
#' @param prefix - trait column prefix eg Bloom%_
#' @param value_name - new name for column eg Bloom1 or Bloom5
#' @param season - trial season
#' @param location - trial location

.al_trait_wide_to_long = function(raw, id_cols, prefix, value_name, season, location){

  # quickfix for columnn name misspelling in bw25 and bw26
  if ('SubloptNo' %in% names(raw)){
    raw = raw %>% rename(SubplotNo = SubloptNo)
  }
  
  trait_cols = names(raw)[stringr::str_starts(names(raw), stringr::fixed(prefix))]
  
  if(length(trait_cols) == 0){
    return (NULL)
  }
  
  # all_trait_cols = unlist(lapply(
  #   AL_TRAIT_PREFIXES,
  #   function(p) names(raw) [str_starts(names(raw), fixed(p))]
  # ))
  # 
  # id_cols = setdiff(names(raw), all_trait_cols)
#   
  raw %>% 
    dplyr::select(all_of(c(id_cols, trait_cols))) %>%
    pivot_longer(
      cols = all_of(trait_cols),
      names_to = 'obs_date_raw',
      values_to = value_name
    ) %>%
    mutate(
      season = season,
      location = location,
      obs_date = .al_parse_date(obs_date_raw, prefix),
      .keep = 'unused'
    ) %>%
    relocate('season', 'location', 'obs_date', .before = 1) %>%
    filter(!is.na(.data[[value_name]]))
}


#' @description process single trial data in the AL rda file
#' @param raw - single environment
#' @param season - trial season eg 2023
#' @param location - trial location eg BW
#' @output list(bloom_long,..., berry_diameter_long) combined for a single location

tidy_single_envt = function(raw, season, location) {

  # quickfix for columnn name misspelling in bw25 and bw26
  if ('SubloptNo' %in% names(raw)){
    # cat('fixed for: ', location, season)
    raw = raw %>% rename(SubplotNo = SubloptNo)
  }
  
  id_cols = .al_resolve_id_cols(raw)
  
  bloom_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['bloom_percent']], 'bloom_percent', season, location
  )
  
  yield_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['yield']], 'yield', season, location
  )
  
  
  list( 
    bloom_long = bloom_long,
    yield_long = yield_long
    )
  
}



#' @description ingest every trial data in the AL rda file
#' @param rda_file - combination of all csv trials in AL
#' @output list(bloom_long,..., berry_diameter_long) combined for all environments
ingest_all_al_dta = function(rda_file){
  trial_info = extract_trial_loc_and_season(rda_file) # gets name, loc, and season for trial
  
  if(nrow(trial_info) == 0){
    warning(
      'ingest_all_al_data(): no trial object matching "<LOC><YY>.ori" found in the al rda file'
    )
  }
  
  results = purrr::pmap(
    trial_info, function(file_name, location, season){
      tidy_single_envt(rda_file[[file_name]], season = season, location = location)
    }
  )
  
  table_names = c(
    'bloom_long', 'yield_long'
  )
  
  stats::setNames(
    lapply(table_names, function(nm) bind_rows(purrr::map(results, nm))),
    table_names
  )
}

# 09/02/2026
# Adding function to process one of the phenos to get the unique lines for genotypes used in trials
get_unique_genotypes = function(yield) {
  unique(yield$Genotype)
}
# 
# 
# a = read_al_rda('data/raw/al/OriginalData.rda')
# b = ingest_all_al_dta(a)
# 
# c = tidy_single_envt(BW23.ori, season = 2023, location = 'BW')
# d = extract_trial_loc_and_season(a)
# 
# 
# results = purrr::pmap(
#   d, function(file_name, location, season){
#     tidy_single_envt(a[[file_name]], season = season, location = location)
#   }
# )
# 
# table_names = c(
#   'bloom_long', 'berry_long', 'frost_long', 'yield_long', 'health_score_long', 'plant_height_long',
#   'notes_long', 'avg_berry_weight_long', 'firmness_long', 'brix_long', 'ta_long', 'berry_diameter_long'
# )
# 
# stats::setNames(
#   lapply(table_names, function(nm) bind_rows(purrr::map(results, nm))),
#   table_names
# )
# 
# 
# 
# 
# 
# raw = BW23.ori
# raw %>% select(matches('Date')) %>% head
