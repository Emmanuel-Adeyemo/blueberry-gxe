# 08/28/2026

library(dplyr)
library(tibble)


#' @description pulls the `MetaData` element out of the combined AL rda obj
#' @param rda_file - combination of all csv trials in AL ie `al_rda`
read_al_metadata = function(rda_file) {
  
  if (!'MetaData' %in% names(rda_file)) {
    stop("read_al_metadata_raw(): rda_file has no 'MetaData' file")
  }
  as_tibble(rda_file[['MetaData']])
}



#' @description pulls the `OffID` element out of the combined AL rda
#' @param rda_file - combination of all csv trials in AL ie `al_rda`
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