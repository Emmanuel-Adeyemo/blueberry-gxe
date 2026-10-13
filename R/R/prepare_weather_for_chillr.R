# 09/22

library(dplyr)
library(lubridate)


prepare_weather_for_chillr = function(weather_raw, location) {
  weather_raw %>%
    mutate(
      Location = location,
      Year = as.integer(.data$YEAR),
      Month = as.integer(.data$MO),
      Day = as.integer(.data$DY),
      JDay = as.integer(format(as.Date(paste(.data$YEAR, .data$MO, .data$DY, sep = "-")), "%j")),
      Hour = as.integer(.data$HR),
      # Hour = as.integer(format(strptime(.data$Time, format = "%I:%M %p"), "%H")),
      # Hour = .parse_12hour_to_24hour(.data$Info),
      Temp = .data$T2M
    ) %>%
    dplyr::select(Location, Year, Month, Day, JDay, Hour, Temp)
}

