# 1. LOAD LIBRARIES 
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)

# 2. DEFINE YEARS (Have data downloaded since 2000 but RHV starts in 2012)
years <- 2012:2025

# 3. IMPORT & COMBINE DATA
all_data <- map_dfr(years, function(y) {
  read_csv(
    paste0("raw data/daily_LEAD_", y, ".csv"),
    col_types = cols(`Method Code` = col_character())
  )
})

# 4. FILTERING & DATA CLEANING
filtered_data <- all_data %>%
  filter(grepl("airport", `Local Site Name`, ignore.case = TRUE)) %>%
  filter(`Arithmetic Mean` > 0) %>%
  mutate(`Date Local` = as.Date(`Date Local`))

# 5. VISUALIZATION - Monthly Trends
filtered_data %>%
  mutate(month = format(`Date Local`, "%Y-%m")) %>%
  group_by(`Local Site Name`, month) %>%
  summarise(mean_value = mean(`Arithmetic Mean`, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = as.Date(paste0(month, "-01")), y = mean_value, color = `Local Site Name`)) +
  geom_line(alpha = 0.7) +
  scale_y_log10() +
  labs(
    title = "Monthly Log Arithmetic Mean (Airport Sites, 2000–2025)",
    x = "Month",
    y = "Arithmetic Mean (log scale)"
  ) +
  theme_minimal()

# 6. EXPLORATION
unique(all_data$`County Name`)


####### WEATHER DATA WORK ###########
library(ncdf4)
library(tidyr)
library(lubridate)
library(doParallel)
library(foreach)


process_weather_year <- function(year, data_path = "raw data/") {
  
  # ---- FILE PATHS ----
  uwnd_file <- paste0(data_path, "uwnd.sig995.", year, ".nc")
  vwnd_file <- paste0(data_path, "vwnd.sig995.", year, ".nc")
  air_file  <- paste0(data_path, "air.sig995.", year, ".nc")
  rhum_file <- paste0(data_path, "rhum.sig995.", year, ".nc")
  pres_file <- paste0(data_path, "pres.sfc.", year, ".nc")
  prw_file  <- paste0(data_path, "pr_wtr.eatm.", year, ".nc")
  
  # ---- LOAD ----
  uwnd_nc <- nc_open(uwnd_file)
  vwnd_nc <- nc_open(vwnd_file)
  air_nc  <- nc_open(air_file)
  rhum_nc <- nc_open(rhum_file)
  pres_nc <- nc_open(pres_file)
  prw_nc  <- nc_open(prw_file)
  
  # ---- EXTRACT VARIABLES ----
  uwnd <- ncvar_get(uwnd_nc, "uwnd")
  vwnd <- ncvar_get(vwnd_nc, "vwnd")
  air  <- ncvar_get(air_nc, "air")
  rhum <- ncvar_get(rhum_nc, "rhum")
  prw  <- ncvar_get(prw_nc, "pr_wtr")
  
  time_dim <- ncvar_get(pres_nc, "time")
  n_time <- length(time_dim)
  
  pres_list <- vector("list", n_time)
  
  for (t in 1:n_time) {
    pres_list[[t]] <- ncvar_get(
      pres_nc,
      "pres",
      start = c(1, 1, t),
      count = c(-1, -1, 1)
    )
  }
  
  pres <- simplify2array(pres_list)
  # ---- TIME ----
  time_raw <- ncvar_get(uwnd_nc, "time")
  time_units <- ncatt_get(uwnd_nc, "time", "units")$value
  
  # Extract origin date (e.g. "days since 1800-01-01")
  origin <- sub(".*since ", "", time_units)
  
  # Convert to Date
  time <- as.Date(time_raw, origin = origin)
  
  # ---- LAT/LON ----
  lon <- ncvar_get(uwnd_nc, "lon")
  lat <- ncvar_get(uwnd_nc, "lat")
  
  # ---- CLOSE FILES ----
  nc_close(uwnd_nc); nc_close(vwnd_nc); nc_close(air_nc)
  nc_close(rhum_nc); nc_close(pres_nc); nc_close(prw_nc)
  
  # ---- BUILD GRID ----
  grid <- expand.grid(lon = lon, lat = lat, time = time)
  
  # ---- FLATTEN ARRAYS ----
  grid$uwnd <- as.vector(uwnd)
  grid$vwnd <- as.vector(vwnd)
  grid$air_temp <- as.vector(air) - 273.15   # Kelvin → Celsius
  grid$rhum <- as.vector(rhum)
  grid$pres <- as.vector(pres) / 1000        # Pa → kPa
  grid$pr_wtr <- as.vector(prw)
  
  # ---- DERIVED VARIABLES ----
  grid <- grid %>%
    mutate(
      wind_speed = sqrt(uwnd^2 + vwnd^2),
      date = as.Date(time)
    )
  
  # ---- FIX LONGITUDE (0–360 → -180–180) ----
  grid <- grid %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon))
  
  return(grid)
}
### BUILD FULL WEATHER DATA SET #####
weather_data <- map_dfr(years, process_weather_year)

### PREP AIRPORT MONITOR DATA
airport_data <- filtered_data %>%
  mutate(
    lat = Latitude,
    lon = Longitude,
    date = as.Date(`Date Local`)
  ) %>%
  select(`Local Site Name`, lat, lon, date, `Arithmetic Mean`)

#### MATCH AIPORT TO NEAREST WEATHER

# Get unique weather grid
weather_coords <- weather_data %>%
  distinct(lat, lon)

# Function to find nearest grid point
find_nearest <- function(lat, lon, weather_coords) {
  dists <- (weather_coords$lat - lat)^2 + (weather_coords$lon - lon)^2
  weather_coords[which.min(dists), ]
}

# Apply to airports
airport_with_grid <- airport_data %>%
  rowwise() %>%
  mutate(
    nearest = list(find_nearest(lat, lon, weather_coords))
  ) %>%
  unnest(nearest)

### MERGE WEATHER ONTO AIRPORTS ### 

final_data <- airport_with_grid %>%
  left_join(
    weather_data,
    by = c("lat" = "lat", "lon" = "lon", "date" = "date")
  )



