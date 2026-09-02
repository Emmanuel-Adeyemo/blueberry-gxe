library(dplyr)
library(tibble)


#' @param rda_file - combination of all csv trials in AL
read_al_metadata = function(rda_file) {
  
  if (!'MetaData' %in% names(rda_file)) {
    stop("read_al_metadata_raw(): rda_file has no 'MetaData' file")
  }
  as_tibble(rda_file[['MetaData']])
}



#' @param rda_file - combination of all csv trials in AL
read_al_off_type = function(rda_file) {
  if (!'OffID' %in% names(rda_file)){
    stop('read_al_off_type(): rda_file has no "OffID" file' )
  }
  off_id = rda_file[['OffID']]
  
  if (is.data.frame(off_id)) {
    return(as_tibble(off_id))
  }
  tibble(PlantID = off_id)
}