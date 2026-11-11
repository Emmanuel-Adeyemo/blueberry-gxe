# 09/22 reworked 10/29 

library(dplyr)
library(readxl)

prepare_station_weather_for_chillr = function(csv_path, location) {
  raw = readr::read_csv(csv_path, show_col_types = FALSE)

  dt_str = sub(" [A-Z]{2,4}$", "", raw[["YYYY-MM-DD hh:mm AM/PM tz"]])
  datetime = as.POSIXct(dt_str, format = "%Y-%m-%d %I:%M %p", tz = "America/Chicago")

  tibble(
    Location = location,
    Year  = as.integer(format(datetime, "%Y")),
    Month = as.integer(format(datetime, "%m")),
    Day   = as.integer(format(datetime, "%d")),
    JDay  = as.integer(format(datetime, "%j")),
    Hour  = as.integer(format(datetime, "%H")),
    Temp  = (raw[["Min Air Temperature"]] - 32) * 5 / 9
  ) %>%
    filter(!is.na(.data$Year)) %>%
    arrange(.data$Year, .data$Month, .data$Day, .data$Hour)
}


# 10/29 -- NC weather prep, mirrors prepare_weather_for_chillr.R's output

prepare_nc_weather_for_chillr = function(xlsx_path, sheet_names = c("2023", "2024", "2025")) {

  read_one_sheet = function(sheet_name) {
    raw = read_excel(xlsx_path, sheet = sheet_name, col_names = FALSE, n_max = 20)
    header_row = which(apply(raw, 1, function(r) any(grepl("Date/Time", r, fixed = TRUE))))
    if (length(header_row) != 1) {
      stop("prepare_nc_weather_for_chillr(): couldn't find a unique 'Date/Time' header row in sheet '", sheet_name, "'.")
    }

    dta = read_excel(xlsx_path, sheet = sheet_name, skip = header_row, col_names = FALSE)
    names(dta)[1:3] = c("datetime", "Temp", "Precip_cm")

    dta %>%
      mutate(
        datetime = as.POSIXct(.data$datetime, tz = "America/New_York"),
        Temp = as.numeric(.data$Temp)
      ) %>%
      filter(!is.na(.data$datetime)) %>%
      arrange(.data$datetime) %>%
      transmute(
        Location = "NC",
        Year  = as.integer(format(.data$datetime, "%Y")),
        Month = as.integer(format(.data$datetime, "%m")),
        Day   = as.integer(format(.data$datetime, "%d")),
        JDay  = as.integer(format(.data$datetime, "%j")),
        Hour  = as.integer(format(.data$datetime, "%H")),
        Temp  = .data$Temp
      )
  }

  bind_rows(lapply(sheet_names, read_one_sheet))
}


# 10/29
#

library(dplyr)

fill_hourly_gaps = function(weather_prepared, location) {
  obs_datetime = as.POSIXct(
    paste(weather_prepared$Year, weather_prepared$Month, weather_prepared$Day, weather_prepared$Hour),
    format = "%Y %m %d %H", tz = "UTC"
  )
  full_grid = tibble(datetime = seq(min(obs_datetime), max(obs_datetime), by = "hour"))

  joined = full_grid %>%
    left_join(
      weather_prepared %>% mutate(datetime = obs_datetime),
      by = "datetime"
    ) %>%
    arrange(.data$datetime) %>%
    mutate(
      Temp_observed = .data$Temp,
      Temp = approx(x = as.numeric(.data$datetime)[!is.na(.data$Temp)],
                     y = .data$Temp[!is.na(.data$Temp)],
                     xout = as.numeric(.data$datetime), rule = 2)$y,
      Temp_imputed = is.na(.data$Temp_observed),
      Location = location,
      Year  = as.integer(format(.data$datetime, "%Y")),
      Month = as.integer(format(.data$datetime, "%m")),
      Day   = as.integer(format(.data$datetime, "%d")),
      JDay  = as.integer(format(.data$datetime, "%j")),
      Hour  = as.integer(format(.data$datetime, "%H"))
    ) %>%
    dplyr::select(.data$Location, .data$Year, .data$Month, .data$Day, .data$JDay, .data$Hour, .data$Temp, .data$Temp_imputed, .data$datetime)

  joined
}



summarize_gap_overlap = function(weather_filled, location, candidate_seasons) {
  obs_date = as.Date(paste(weather_filled$Year, weather_filled$Month, weather_filled$Day, sep = "-"))

  purrr::map_dfr(candidate_seasons, function(season) {
    anchor = season_anchor_date(season)
    days_since_anchor = as.numeric(obs_date - anchor)

    in_chill = days_since_anchor >= 30 & days_since_anchor <= 180
    in_heat  = days_since_anchor >= 122 & days_since_anchor <= 180
    in_frost_scoring = days_since_anchor >= 195 & days_since_anchor <= 211

    tibble(
      location = location,
      season = season,
      imputed_hours_in_chill_window = sum(weather_filled$Temp_imputed[in_chill]),
      imputed_hours_in_heat_window = sum(weather_filled$Temp_imputed[in_heat]),
      imputed_hours_in_frost_scoring_window = sum(weather_filled$Temp_imputed[in_frost_scoring])
    )
  })
}