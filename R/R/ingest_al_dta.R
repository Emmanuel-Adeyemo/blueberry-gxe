# blueberry/R/R
# legacy R scripts already processed AL field data. 
# Rerun to try rxn norm for bloom data instead of the 
# interpolation used in the legacy script

# 
# EV23.ori <- read_csv("Data/EV23_071526.csv")
# EV24.ori <- read_csv("Data/EV24_071526.csv")
# EV25.ori <- read_csv("Data/EV25_071526.csv")
# EV26.ori <- read_csv("Data/EV26_071526.csv")
# 
# BW23.ori <- read_csv("Data/BW23_071526.csv")
# BW24.ori <- read_csv("Data/BW24_071526.csv")
# BW25.ori <- read_csv("Data/BW25_071526.csv")
# BW26.ori <- read_csv("Data/BW26_071526.csv")
# 
# OffID <- read_csv("Data/OffType.csv")$PlantID
# MetaData <- read_csv("Data/MetaData.csv")
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


AL_TRAIT_PREFIXES <- c(
  bloom_percent = "Bloom%_",
  berry_percent = "Berry%_", # fruit-set progression, same shape as bloom
  frost_percent = "Frost%_", # will not use
  yield = "Yield_",
  health_score = "HealthScore_", # dated visit score 
  plant_height = "PlantHeight_", # dated measurement 
  notes = "Notes_", # free text, handled separately in tidy_al_environment()
  avg_berry_weight = "AverageBerryWeight_", # out of scope 
  firmness = "Firmness_", # out of scope
  brix = "Brix_", # out of scope
  ta = "TA_", # out of scope
  berry_diameter = "BerryDiameter_" # out of scope
)



AL_ID_COLS <- c(
  "PlantID", "BlockNo", "SubplotNo", "RepNo", "Genotype", "Ecotype",
  "HealthSummary", "PlantHeightSummary", "Harvested"
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
  list(location = m[1, 2], year = 2000L + as.integer(m[1, 3]))
}



#' @description identifies the trial files within the loaded rda file
#' @param rda_file - combination of all csv trials in AL
#' @output tibble(file_name, location, year)
#' 
extract_trial_loc_and_year = function(rda_file){
  
  file_names = names(rda_file)
  parsed = lapply(file_names, .al_parse_file_name)
  is_trial = !vapply(parsed, is.null, logical(1)) # sanity checks for trials, should be 8 and 2 nulls
  # 
  unexpected = setdiff(file_names[!is_trial], c('MetaData', 'OffID'))
  if (length(unexpected) > 0) {
    warning('unexpected file extracted from AL rda file',
            paste(unexpected, collapse = ', '))
  }
  
  tibble::tibble(
    file_name = file_names[is_trial],
    location = vapply(parsed[is_trial], `[[`, character(1), "location"),
    year = vapply(parsed[is_trial], `[[`, integer(1), "year")
  )
  
}

# rm(file_names, parsed, is_trial)
# file_names = a
# parsed = lapply(file_names, .al_parse_file_name)
# is_trial = !vapply(parsed, is.null, logical(1)) 







#' @description parse colname into date eg Yield_MM_DD_YY
#' @param col_name - character vector of column names
#' @param prefix - the trait prefix to strip eg Bloom%_ 

.al_parse_date = function(col_names, prefix){
  raw_date = stringr::str_sub(col_names, nchar(prefix) + 1) 
  
  # for future, in case date has / or -
  raw_date = stringr::str_remove_all(raw_date, '[^0-9]')
  
  n = nchar(raw_date)
  out = as.Date(rep(NA_character_, length(raw_date)))
  out[n == 6] = as.Date(raw_date[n==6], format = '%m%d%y')
  out[n == 8] = as.Date(raw_date[n==8], format = '%m%d%Y')
  
  unexpected = !(n %in% c(6,8))
  if (any(unexpected)) {
    warning(
      'al_parse_date(): unexpected date-string length in columns: ',
      paste(col_names[unexpected], collaps=', ')
    )
  }
  out
  
}




#' @param raw - original wide data ie single envt
#' @param output - character vector of id column names present in raw
.al_resolve_id_cols = function(raw) {
  
  all_trait_cols = unlist(lapply(
    AL_TRAIT_PREFIXES,
    function(p) names(raw) [stringr::str_starts(names(raw), stringr::fixed(p))]
  ))
  id_cols = intersect(AL_ID_COLS, names(raw))
  
  unexpected = setdiff(names(raw), c(all_trait_cols))
  if(length(unexpected) > 0){
    warning(
      '.al_resolve_id_cols(): dropping unrecognised column. Check known identifiers AL_TRAIT_PREFIXES and ALL_ID_COLS for details. ',
      paste(unexpected, collapse = ', ')
    )
  }
  id_cols
}



#' @param raw - original wide data ie single envt
#' @param prefix - trait column prefix eg Bloom%_
#' @param value_name - new name for column eg Bloom1 or Bloom5
#' @param year - trial year
#' @param location - trial location

.al_trait_wide_to_long = function(raw, id_cols, prefix, value_name, year, location){
  
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
    select(all_of(c(id_cols, trait_cols))) %>%
    pivot_longer(
      cols = all_of(trait_cols),
      names_to = 'obs_date_raw',
      values_to = value_name
    ) %>%
    mutate(
      year = year,
      location = location,
      obs_date = .al_parse_date(.data$obs_date_raw, prefix),
      .keep = 'unused'
    ) %>%
    relocate('year', 'location', 'obs_date', .before = 1) %>%
    filter(!is.na(.data[[value_name]]))
}


#' @description process single trial data in the AL rda file
#' @param raw - single environment
#' @param year - trial year eg 2023
#' @param location - trial location eg BW
#' #' @output list(bloom_long,..., berry_diameter_long) combined for a single location

tidy_single_envt = function(raw, year, location) {
  
  id_cols = .al_resolve_id_cols(raw)
  
  bloom_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['bloom_percent']], 'bloom', year, location
  )
  
  berry_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['berry_percent']], 'berry', year, location
  )
  
  frost_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['frost_percent']], 'frost', year, location
  )
  
  yield_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['yield']], 'yield', year, location
  )
  
  health_score_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['health_score']], 'health_score', year, location
  )
  
  plant_height_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['plant_height']], 'plant_height', year, location
  )
  
  notes_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['notes']], 'notes', year, location
  )
  
  avg_berry_weight_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['avg_berry_weight']], 'avg_berry_weight', year, location
  )
  
  firmness_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['firmness']], 'firmness', year, location
  )
  
  brix_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['brix']], 'brix', year, location
  )
  
  ta_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['ta']], 'ta', year, location
  )
  
  berry_diameter_long = .al_trait_wide_to_long(
    raw, id_cols, AL_TRAIT_PREFIXES[['berry_diameter']], 'berry_diameter', year, location
  )
  
  
  list( 
    bloom_long = bloom_long,
    berry_long = berry_long,
    frost_long = frost_long,
    yield_long = yield_long,
    health_score_long = health_score_long,
    plant_height_long = plant_height_long,
    notes_long = notes_long,
    avg_berry_weight_long = avg_berry_weight_long,
    firmness_long = firmness_long,
    brix_long = brix_long,
    ta_long = ta_long,
    berry_diameter_long = berry_diameter_long
    )
  
}



#' @description ingest every trial data in the AL rda file
#' @param rda_file - combination of all csv trials in AL
#' @output list(bloom_long,..., berry_diameter_long) combined for all environments
ingest_all_al_dta = function(rda_file){
  trial_info = extract_trial_loc_and_year(rda_file) # gets name, loc, and year for trial
  
  if(nrow(trial_info) == 0){
    warning(
      'ingest_all_al_data(): no trial object matching "<LOC><YY>.ori" found in the al rda file'
    )
  }
  
  results = purrr::pmap(
    trial_info, function(file_name, location, year){
      tidy_single_envt(rda_file[[file_name]], year = year, location = location)
    }
  )
  
  table_names = c(
    'bloom_long', 'berry_long', 'frost_long', 'yield_long', 'health_score_long', 'plant_height_long',
    'notes_long', 'avg_berry_weight_long', 'firmness_long', 'brix_long', 'ta_long', 'berry_diameter_long'
  )
  
  stats::setNames(
    lapply(table_names, function(nm) bind_rows(purrr::map(results, nm))),
    table_names
  )
}
# 
# 
# a = read_al_rda('data/raw/al/OriginalData.rda')
# b = ingest_all_al_dta(a)
# 
# c = tidy_single_envt(BW23.ori, year = 2023, location = 'BW')
# d = extract_trial_loc_and_year(a)
# 
# 
# results = purrr::pmap(
#   d, function(file_name, location, year){
#     tidy_single_envt(a[[file_name]], year = year, location = location)
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
