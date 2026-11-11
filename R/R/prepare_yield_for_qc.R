# 09/18

library(dplyr)

aggregate_yield_for_qc = function(yield_dta){


  has_age = 'Age' %in% names(yield_dta)

  plant_group_cols = c('season', 'location', 'BlockNo', 'SubplotNo', 'PlantID', 'Genotype')
  if (has_age) plant_group_cols = c(plant_group_cols, 'Age')

  plant_totals = yield_dta %>%
    filter(HealthSummary > 0, Harvested == 'Yes') %>%
    group_by(across(all_of(plant_group_cols))) %>%
    summarise(plant_total_yield = sum(yield, na.rm = TRUE), .groups = 'drop')

  if (!has_age) {
    return(
      plant_totals %>%
        group_by(season, location, BlockNo, SubplotNo, Genotype) %>%
        summarise(TotalYield = mean(plant_total_yield, na.rm = TRUE), .groups = 'drop')
    )
  }

  plant_totals %>%
    group_by(season, location, BlockNo, SubplotNo, Genotype) %>%
    summarise(
      TotalYield = mean(plant_total_yield, na.rm = TRUE),
      Age = {
        n_age = dplyr::n_distinct(Age, na.rm = TRUE)
        if (n_age > 1) {
          stop(
            "aggregate_yield_for_qc(): Age is not constant within subplot (",
            "season=", dplyr::first(season), ", location=", dplyr::first(location),
            ", BlockNo=", dplyr::first(BlockNo), ", SubplotNo=", dplyr::first(SubplotNo),
            ", Genotype=", dplyr::first(Genotype), ") -- likely a replanting event ",
            "put some reps on a different planting cohort than their subplot-mates. ",
            "Resolve before using Age as a fixed effect."
          )
        }
        dplyr::first(Age)
      },
      .groups = 'drop'
    )
}


# 10/21 
#
# compute_yield_log_offset = function(al_yield_long_clean, nc_yield_clean, fl_yield_clean) {

#   al_total = aggregate_yield_for_qc(al_yield_long_clean)$TotalYield
#   nc_total = nc_yield_clean$y
#   fl_total = fl_yield_clean$y

#   all_yield = c(al_total, nc_total, fl_total)
#   all_yield = all_yield[!is.na(all_yield)]

#   min_nonzero = min(all_yield[all_yield > 0])
#   min_nonzero / 2
# }


# 10/16 quality traits 
aggregate_quality_trait_for_qc = function(trait_long, value_col) {

  plant_means = trait_long %>%
    filter(HealthSummary > 0, Harvested == 'Yes') %>%
    group_by(season, location, BlockNo, SubplotNo, PlantID, Genotype) %>%
    summarise(plant_mean = mean(.data[[value_col]], na.rm = TRUE), .groups = 'drop') %>%
    mutate(plant_mean = ifelse(is.nan(plant_mean), NA_real_, plant_mean))

  plant_means %>%
    group_by(season, location, BlockNo, SubplotNo, Genotype) %>%
    summarise(!!value_col := mean(plant_mean, na.rm = TRUE), .groups = 'drop')
}