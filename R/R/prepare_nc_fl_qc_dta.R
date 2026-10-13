# 10/02

library(dplyr)


prepare_nc_qc_dta = function(nc_dta, trait_col) {
  nc_dta %>%
    transmute(
      genotype = genotype,
      block = subplot,
      y = .data[[trait_col]]
    )
}


prepare_fl_qc_data = function(fl_dta, trait_col) {
  fl_dta %>%
    transmute(
      genotype = genotype,
      block = site_detail,
      y = .data[[trait_col]]
    )
}