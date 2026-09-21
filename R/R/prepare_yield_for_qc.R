# 09/18

library(dplyr)

aggregate_yield_for_qc = function(yield_dta){

  plant_totals = yield_dta %>%
    filter(HealthSummary > 0, Harvested == 'Yes') %>%
    group_by(season, location, BlockNo, SubplotNo, PlantID, Genotype) %>%
    summarise(plant_total_yield = sum(yield, na.rm = TRUE), .groups = 'drop')

  plant_totals %>%
    group_by(season, location, BlockNo, SubplotNo, Genotype) %>%
    summarise(TotalYield = mean(plant_total_yield, na.rm = TRUE), .groups = 'drop')
}