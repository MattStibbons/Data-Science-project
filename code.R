###########################################################################
# AVIATION GASOLINE LEAD PROJECT
# Base code + weather merge + airport screening
# Time window: 2016-01-01 to 2025-01-01
###########################################################################

# -------------------------------------------------------------------------
# 1. LOAD PACKAGES
# -------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,
  sf,
  lubridate,
  tigris,
  units,
  ncdf4,
  fixest,
  broom,
  stringr
)

options(tigris_use_cache = TRUE)

# -------------------------------------------------------------------------
# 2. SETTINGS
# -------------------------------------------------------------------------
raw_data_folder <- "Raw data"
output_folder <- "outputs"
fig_dir <- file.path(output_folder, "figures")

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

analysis_start <- as.Date("2016-01-01")
analysis_end <- as.Date("2025-01-01")   # exclusive end date

years <- lubridate::year(analysis_start):lubridate::year(analysis_end)

max_distance_km <- 10
policy_date <- as.Date("2022-01-01")

# Airport screening for plots
min_airport_months <- 24

# TWFE sample restrictions
min_pre_months <- 12
min_post_months <- 3
lead_log_constant <- 0.001

# -------------------------------------------------------------------------
# 3. READ EPA DAILY LEAD DATA
# -------------------------------------------------------------------------
lead_files <- file.path(raw_data_folder, paste0("daily_LEAD_", years, ".csv"))

file_check <- tibble(
  year = years,
  file = lead_files,
  exists = file.exists(lead_files)
)
print(file_check)

lead_files_existing <- lead_files[file.exists(lead_files)]

if (length(lead_files_existing) == 0) {
  stop("No daily_LEAD_YYYY.csv files found. Check raw_data_folder and file names.")
}

read_lead_file <- function(file_path) {
  readr::read_csv(
    file_path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    mutate(source_file = basename(file_path))
}

lead_all <- purrr::map_dfr(lead_files_existing, read_lead_file)

# -------------------------------------------------------------------------
# 4. CLEAN LEAD DATA AND CONSTRUCT MONITOR IDS / COUNTY FIPS
# -------------------------------------------------------------------------
lead_clean <- lead_all %>%
  rename(
    date = `Date Local`,
    lead_conc = `Arithmetic Mean`,
    state_code_raw = `State Code`,
    county_code_raw = `County Code`,
    site_num_raw = `Site Num`,
    state_name = `State Name`,
    county_name = `County Name`,
    monitor_lat = Latitude,
    monitor_lon = Longitude
  ) %>%
  mutate(
    date = as.Date(date),
    year = lubridate::year(date),
    month = lubridate::floor_date(date, "month"),
    lead_conc = as.numeric(lead_conc),
    monitor_lat = as.numeric(monitor_lat),
    monitor_lon = as.numeric(monitor_lon),
    
    state_code_raw = str_remove(str_trim(as.character(state_code_raw)), "\\.0$"),
    county_code_raw = str_remove(str_trim(as.character(county_code_raw)), "\\.0$"),
    site_num_raw = str_remove(str_trim(as.character(site_num_raw)), "\\.0$"),
    
    state_code_2 = str_pad(state_code_raw, width = 2, side = "left", pad = "0"),
    county_code_3 = str_pad(county_code_raw, width = 3, side = "left", pad = "0"),
    site_num_4 = str_pad(site_num_raw, width = 4, side = "left", pad = "0"),
    
    county_fips = paste0(state_code_2, county_code_3),
    
    # Keep this consistent throughout the project:
    # state + county + site number.
    monitor_id = paste0(state_code_2, county_code_3, site_num_4)
  ) %>%
  filter(
    date >= analysis_start,
    date < analysis_end,
    !is.na(date),
    !is.na(lead_conc),
    lead_conc >= 0,
    !is.na(monitor_lat),
    !is.na(monitor_lon)
  )

print(
  lead_clean %>%
    summarise(
      n_daily_obs = n(),
      n_monitors = n_distinct(monitor_id),
      n_counties = n_distinct(county_fips),
      min_date = min(date),
      max_date = max(date)
    )
)

# Santa Clara County should be 06085
print(
  lead_clean %>%
    filter(county_fips == "06085") %>%
    summarise(
      n_daily_obs = n(),
      n_monitors = n_distinct(monitor_id),
      min_date = min(date),
      max_date = max(date)
    )
)

# -------------------------------------------------------------------------
# 5. CREATE MONITOR SPATIAL DATA
# -------------------------------------------------------------------------
monitors <- lead_clean %>%
  distinct(
    monitor_id,
    county_fips,
    state_name,
    county_name,
    monitor_lat,
    monitor_lon
  )

monitors_sf <- sf::st_as_sf(
  monitors,
  coords = c("monitor_lon", "monitor_lat"),
  crs = 4326,
  remove = FALSE
)

# -------------------------------------------------------------------------
# 6. READ AIRPORT DATA
# -------------------------------------------------------------------------
airport_file <- file.path(raw_data_folder, "airports_raw.csv")
airport_url <- "https://davidmegginson.github.io/ourairports-data/airports.csv"

if (!file.exists(airport_file)) {
  message("airports_raw.csv not found in Raw data/. Downloading from OurAirports...")
  download.file(airport_url, airport_file, mode = "wb")
}

airports_raw <- readr::read_csv(airport_file, show_col_types = FALSE)

airports_us <- airports_raw %>%
  filter(
    iso_country == "US",
    type %in% c("small_airport", "medium_airport", "large_airport")
  ) %>%
  transmute(
    airport_id = as.character(id),
    airport_ident = ident,
    gps_code = gps_code,
    local_code = local_code,
    airport_name = name,
    airport_type = type,
    airport_lat = latitude_deg,
    airport_lon = longitude_deg
  ) %>%
  filter(
    !is.na(airport_lat),
    !is.na(airport_lon)
  )

airports_sf <- sf::st_as_sf(
  airports_us,
  coords = c("airport_lon", "airport_lat"),
  crs = 4326,
  remove = FALSE
)

# -------------------------------------------------------------------------
# 7. MATCH AIRPORTS TO COUNTY FIPS
# -------------------------------------------------------------------------
counties_sf <- tigris::counties(cb = TRUE, year = 2022, class = "sf") %>%
  sf::st_transform(4326) %>%
  select(
    airport_county_fips = GEOID,
    airport_county_name = NAME
  )

airports_sf <- sf::st_join(
  airports_sf,
  counties_sf,
  join = sf::st_intersects,
  left = TRUE
)

print(
  airports_sf %>%
    sf::st_drop_geometry() %>%
    summarise(
      n_airports = n(),
      n_missing_county = sum(is.na(airport_county_fips))
    )
)

# -------------------------------------------------------------------------
# 8. IDENTIFY REID-HILLVIEW AIRPORT
# -------------------------------------------------------------------------
hillview_airport <- airports_sf %>%
  sf::st_drop_geometry() %>%
  filter(
    str_detect(str_to_lower(airport_name), "reid|hillview") |
      airport_ident %in% c("KRHV", "RHV") |
      gps_code %in% c("KRHV", "RHV") |
      local_code == "RHV"
  )

print(hillview_airport)

if (nrow(hillview_airport) == 0) {
  warning("Reid-Hillview was not identified. Check airport_name / ident / gps_code / local_code.")
}

# -------------------------------------------------------------------------
# 9. AIRPORT-TO-MONITOR DISTANCE MATCHING
# -------------------------------------------------------------------------
# Distance bands:
#   0-1 km
#   1-3 km
#   3-5 km
#   0-4 km, overlapping group used for TWFE
monitors_proj <- sf::st_transform(monitors_sf, 5070) %>%
  mutate(monitor_row_id = row_number())

airports_proj <- sf::st_transform(airports_sf, 5070) %>%
  mutate(airport_row_id = row_number())

within_list <- sf::st_is_within_distance(
  airports_proj,
  monitors_proj,
  dist = units::set_units(max_distance_km, "km")
)

airport_monitor_pairs_base <- tibble(
  airport_row_id = rep(seq_along(within_list), lengths(within_list)),
  monitor_row_id = unlist(within_list)
) %>%
  filter(!is.na(monitor_row_id)) %>%
  mutate(
    distance_m = as.numeric(sf::st_distance(
      airports_proj[airport_row_id, ],
      monitors_proj[monitor_row_id, ],
      by_element = TRUE
    )),
    distance_km = distance_m / 1000
  )

airport_monitor_pairs <- bind_rows(
  airport_monitor_pairs_base %>%
    mutate(
      distance_band = case_when(
        distance_km <= 1 ~ "0-1 km",
        distance_km > 1 & distance_km <= 3 ~ "1-3 km",
        distance_km > 3 & distance_km <= 5 ~ "3-5 km",
        TRUE ~ NA_character_
      )
    ),
  airport_monitor_pairs_base %>%
    filter(distance_km <= 4) %>%
    mutate(distance_band = "0-4 km")
) %>%
  filter(!is.na(distance_band)) %>%
  mutate(
    distance_band = factor(
      distance_band,
      levels = c("0-1 km", "1-3 km", "3-5 km", "0-4 km")
    )
  ) %>%
  left_join(
    airports_proj %>%
      sf::st_drop_geometry() %>%
      select(
        airport_row_id,
        airport_id,
        airport_ident,
        gps_code,
        local_code,
        airport_name,
        airport_type,
        airport_county_fips,
        airport_county_name
      ),
    by = "airport_row_id"
  ) %>%
  left_join(
    monitors_proj %>%
      sf::st_drop_geometry() %>%
      select(
        monitor_row_id,
        monitor_id,
        county_fips,
        state_name,
        county_name,
        monitor_lat,
        monitor_lon
      ),
    by = "monitor_row_id"
  ) %>%
  mutate(
    same_county = county_fips == airport_county_fips,
    is_reid_hillview = str_detect(str_to_lower(airport_name), "reid|hillview") |
      airport_ident %in% c("KRHV", "RHV") |
      gps_code %in% c("KRHV", "RHV") |
      local_code == "RHV"
  )

print(airport_monitor_pairs %>% count(distance_band))

print(
  airport_monitor_pairs %>%
    filter(is_reid_hillview) %>%
    arrange(distance_km)
)

# -------------------------------------------------------------------------
# 10. JOIN DAILY LEAD OBSERVATIONS TO AIRPORT-MONITOR PAIRS
# -------------------------------------------------------------------------
lead_airport_analysis <- lead_clean %>%
  inner_join(
    airport_monitor_pairs %>%
      select(
        airport_id,
        airport_ident,
        gps_code,
        local_code,
        airport_name,
        airport_type,
        airport_county_fips,
        airport_county_name,
        monitor_id,
        distance_km,
        distance_band,
        same_county,
        is_reid_hillview
      ),
    by = "monitor_id"
  )

print(
  lead_airport_analysis %>%
    summarise(
      n_airport_monitor_daily_obs = n(),
      n_airports = n_distinct(airport_id),
      n_monitors = n_distinct(monitor_id),
      min_date = min(date),
      max_date = max(date)
    )
)

print(
  lead_airport_analysis %>%
    group_by(distance_band) %>%
    summarise(
      n_airports = n_distinct(airport_id),
      n_monitors = n_distinct(monitor_id),
      n_daily_obs = n(),
      min_date = min(date),
      max_date = max(date),
      .groups = "drop"
    )
)

# -------------------------------------------------------------------------
# 11. EXTRACT DAILY WEATHER FOR EACH MONITOR-DATE
# -------------------------------------------------------------------------
# Required NCEP files in Raw data/:
#   air.sig995.YYYY.nc
#   uwnd.sig995.YYYY.nc
#   vwnd.sig995.YYYY.nc
#   rhum.sig995.YYYY.nc
#   pres.sfc.YYYY.nc
#
# Important:
# cache name includes the analysis window so old 2000-2025 weather cache
# will not be accidentally reused.
weather_cache_file <- file.path(
  output_folder,
  paste0(
    "weather_monitor_date_cache_",
    format(analysis_start, "%Y%m%d"),
    "_",
    format(analysis_end, "%Y%m%d"),
    ".rds"
  )
)

extract_weather_for_year <- function(obs_year, z, raw_data_folder) {
  if (nrow(obs_year) == 0) return(NULL)
  
  message("Extracting weather for year: ", z)
  
  air_file <- file.path(raw_data_folder, paste0("air.sig995.", z, ".nc"))
  
  if (!file.exists(air_file)) {
    warning("Missing air file for year ", z, ": ", air_file)
    
    return(
      obs_year %>%
        mutate(
          temp_k = NA_real_,
          uwnd = NA_real_,
          vwnd = NA_real_,
          rhum = NA_real_,
          pres = NA_real_,
          air_temp = NA_real_,
          wind_speed = NA_real_,
          wind_dir = NA_real_
        ) %>%
        select(
          monitor_id,
          date,
          air_temp,
          wind_speed,
          wind_dir,
          rhum,
          pres
        )
    )
  }
  
  nc_ref <- ncdf4::nc_open(air_file)
  lon_vals <- ncdf4::ncvar_get(nc_ref, "lon")
  lat_vals <- ncdf4::ncvar_get(nc_ref, "lat")
  ncdf4::nc_close(nc_ref)
  
  obs_year <- obs_year %>%
    mutate(
      lon_360 = if_else(monitor_lon < 0, monitor_lon + 360, monitor_lon),
      day_idx = lubridate::yday(date),
      lat_idx = purrr::map_int(monitor_lat, ~ which.min(abs(lat_vals - .x))),
      lon_idx = purrr::map_int(lon_360, ~ which.min(abs(lon_vals - .x)))
    )
  
  pull_weather <- function(var_name, file_prefix) {
    f_path <- file.path(raw_data_folder, paste0(file_prefix, ".", z, ".nc"))
    
    if (!file.exists(f_path)) {
      warning("Missing weather file: ", f_path)
      return(rep(NA_real_, nrow(obs_year)))
    }
    
    nc <- ncdf4::nc_open(f_path)
    
    vals <- purrr::map_dbl(seq_len(nrow(obs_year)), function(i) {
      ncdf4::ncvar_get(
        nc,
        var_name,
        start = c(obs_year$lon_idx[i], obs_year$lat_idx[i], obs_year$day_idx[i]),
        count = c(1, 1, 1)
      )
    })
    
    ncdf4::nc_close(nc)
    
    vals
  }
  
  obs_year$temp_k <- pull_weather("air", "air.sig995")
  obs_year$uwnd <- pull_weather("uwnd", "uwnd.sig995")
  obs_year$vwnd <- pull_weather("vwnd", "vwnd.sig995")
  obs_year$rhum <- pull_weather("rhum", "rhum.sig995")
  obs_year$pres <- pull_weather("pres", "pres.sfc") / 1000  # Pa to kPa
  
  obs_year %>%
    mutate(
      air_temp = temp_k - 273.15,
      wind_speed = sqrt(uwnd^2 + vwnd^2),
      wind_dir = (atan2(-uwnd, -vwnd) * 180 / pi + 360) %% 360
    ) %>%
    select(
      monitor_id,
      date,
      air_temp,
      wind_speed,
      wind_dir,
      rhum,
      pres
    )
}

if (file.exists(weather_cache_file)) {
  message("Loading cached weather data: ", weather_cache_file)
  weather_daily_clean <- readRDS(weather_cache_file) %>%
    mutate(
      date = as.Date(date),
      month = lubridate::floor_date(date, "month")
    )
} else {
  weather_input <- lead_airport_analysis %>%
    distinct(monitor_id, date, monitor_lat, monitor_lon) %>%
    mutate(year = lubridate::year(date)) %>%
    filter(
      date >= analysis_start,
      date < analysis_end,
      year %in% years
    )
  
  print(
    weather_input %>%
      summarise(
        n_monitor_dates = n(),
        n_monitors = n_distinct(monitor_id),
        min_date = min(date),
        max_date = max(date)
      )
  )
  
  weather_daily_clean <- purrr::map_dfr(years, function(z) {
    obs_year <- weather_input %>% filter(year == z)
    extract_weather_for_year(obs_year, z, raw_data_folder)
  }) %>%
    mutate(
      date = as.Date(date),
      month = lubridate::floor_date(date, "month")
    )
  
  saveRDS(weather_daily_clean, weather_cache_file)
  
  readr::write_csv(
    weather_daily_clean,
    file.path(
      output_folder,
      paste0(
        "weather_monitor_date_cache_",
        format(analysis_start, "%Y%m%d"),
        "_",
        format(analysis_end, "%Y%m%d"),
        ".csv"
      )
    )
  )
}

print(
  weather_daily_clean %>%
    summarise(
      n_obs = n(),
      pct_missing_temp = mean(is.na(air_temp)) * 100,
      pct_missing_wind = mean(is.na(wind_speed)) * 100,
      pct_missing_rhum = mean(is.na(rhum)) * 100,
      pct_missing_pres = mean(is.na(pres)) * 100,
      mean_air_temp = mean(air_temp, na.rm = TRUE),
      mean_wind_speed = mean(wind_speed, na.rm = TRUE)
    )
)

# -------------------------------------------------------------------------
# 12. MERGE WEATHER INTO DAILY LEAD-AIRPORT DATA
# -------------------------------------------------------------------------
lead_airport_weather <- lead_airport_analysis %>%
  left_join(
    weather_daily_clean,
    by = c("monitor_id", "date", "month")
  ) %>%
  mutate(
    wind_dir_rad = wind_dir * pi / 180,
    wind_dir_sin = sin(wind_dir_rad),
    wind_dir_cos = cos(wind_dir_rad)
  )

print(
  lead_airport_weather %>%
    summarise(
      n_obs = n(),
      pct_missing_temp = mean(is.na(air_temp)) * 100,
      pct_missing_wind = mean(is.na(wind_speed)) * 100,
      pct_missing_rhum = mean(is.na(rhum)) * 100,
      pct_missing_pres = mean(is.na(pres)) * 100
    )
)

print(
  lead_airport_weather %>%
    group_by(distance_band) %>%
    summarise(
      n_airports = n_distinct(airport_id),
      n_monitors = n_distinct(monitor_id),
      n_daily_obs = n(),
      pct_missing_temp = mean(is.na(air_temp)) * 100,
      min_date = min(date),
      max_date = max(date),
      .groups = "drop"
    )
)

# -------------------------------------------------------------------------
# 13. AIRPORT-MONTH AGGREGATION WITH WEATHER CONTROLS
# -------------------------------------------------------------------------
airport_monthly_weather <- lead_airport_weather %>%
  group_by(
    distance_band,
    airport_id,
    airport_ident,
    gps_code,
    local_code,
    airport_name,
    airport_type,
    airport_county_fips,
    airport_county_name,
    month
  ) %>%
  summarise(
    mean_lead = mean(lead_conc, na.rm = TRUE),
    n_obs = n(),
    n_monitors = n_distinct(monitor_id),
    is_reid_hillview = any(is_reid_hillview),
    
    air_temp = mean(air_temp, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    rhum = mean(rhum, na.rm = TRUE),
    pres = mean(pres, na.rm = TRUE),
    wind_dir_sin = mean(wind_dir_sin, na.rm = TRUE),
    wind_dir_cos = mean(wind_dir_cos, na.rm = TRUE),
    
    .groups = "drop"
  )

readr::write_csv(
  airport_monitor_pairs,
  file.path(output_folder, "airport_monitor_pairs_within_10km_2016_2024.csv")
)

readr::write_csv(
  lead_airport_weather,
  file.path(output_folder, "lead_airport_daily_with_weather_2016_2024.csv")
)

readr::write_csv(
  airport_monthly_weather,
  file.path(output_folder, "airport_monthly_weather_2016_2024.csv")
)

saveRDS(
  airport_monthly_weather,
  file.path(output_folder, "airport_monthly_weather_2016_2024.rds")
)

print(
  airport_monthly_weather %>%
    distinct(distance_band, airport_id, airport_name) %>%
    count(distance_band)
)

# -------------------------------------------------------------------------
# 14. AIRPORT SCREENING FOR PANEL PLOTS
# -------------------------------------------------------------------------
# This is the airport screening that previously reduced the counts:
# only keep airport x distance_band panels with at least min_airport_months months.
airport_keep <- airport_monthly_weather %>%
  group_by(distance_band, airport_id, airport_name) %>%
  summarise(
    n_months = n_distinct(month),
    avg_lead = mean(mean_lead, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_months >= min_airport_months)

airport_monthly_plot <- airport_monthly_weather %>%
  inner_join(
    airport_keep,
    by = c("distance_band", "airport_id", "airport_name")
  )

print(
  airport_monthly_plot %>%
    distinct(distance_band, airport_id, airport_name) %>%
    count(distance_band)
)

readr::write_csv(
  airport_keep,
  file.path(output_folder, "airport_keep_plot_24_month_filter_2016_2024.csv")
)

# -------------------------------------------------------------------------
# 15. AIRPORT-LEVEL PANEL PLOTS BY DISTANCE BAND
# -------------------------------------------------------------------------
plot_airport_band <- function(band_label) {
  df_band <- airport_monthly_plot %>%
    filter(distance_band == band_label)
  
  ggplot(df_band, aes(x = month, y = mean_lead, group = airport_id)) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = policy_date, linetype = "dashed") +
    facet_wrap(~ airport_name, scales = "free_y", ncol = 3) +
    labs(
      title = "Monthly Ambient Lead Concentration by Airport",
      subtitle = paste0(
        "Distance band: ",
        band_label,
        "; sample: 2016-01-01 to 2025-01-01"
      ),
      x = "Month",
      y = "Mean daily ambient lead concentration"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

all_bands <- levels(airport_monthly_weather$distance_band)

airport_band_plots <- set_names(
  map(all_bands, plot_airport_band),
  all_bands
)

walk2(
  airport_band_plots,
  names(airport_band_plots),
  ~ ggsave(
    filename = file.path(
      fig_dir,
      paste0(
        "airport_monthly_panel_",
        str_replace_all(.y, "[^A-Za-z0-9]+", "_"),
        "_2016_2024.png"
      )
    ),
    plot = .x,
    width = 13,
    height = 9,
    dpi = 300
  )
)

print(airport_band_plots[["0-4 km"]])

# -------------------------------------------------------------------------
# 16. BUILD TWFE SAMPLE FOR 0-4 KM
# -------------------------------------------------------------------------
twfe_data <- airport_monthly_weather %>%
  filter(distance_band == "0-4 km") %>%
  mutate(
    post = month >= policy_date,
    treated = is_reid_hillview,
    treated_post = as.integer(treated & post),
    log_lead = log(mean_lead + lead_log_constant)
  ) %>%
  filter(
    !is.na(log_lead),
    is.finite(log_lead),
    !is.na(air_temp),
    !is.na(wind_speed),
    !is.na(rhum),
    !is.na(pres),
    !is.na(wind_dir_sin),
    !is.na(wind_dir_cos)
  )

coverage_all <- twfe_data %>%
  group_by(airport_id, airport_name, treated) %>%
  summarise(
    n_pre_months = n_distinct(month[month < policy_date]),
    n_post_months = n_distinct(month[month >= policy_date]),
    min_month = min(month),
    max_month = max(month),
    .groups = "drop"
  ) %>%
  arrange(desc(treated), airport_name)

print(coverage_all)

readr::write_csv(
  coverage_all,
  file.path(output_folder, "twfe_airport_coverage_all_2016_2024.csv")
)

# Keep airports with enough pre and post observations.
airport_keep_twfe <- twfe_data %>%
  group_by(airport_id, airport_name) %>%
  summarise(
    n_pre_months = n_distinct(month[month < policy_date]),
    n_post_months = n_distinct(month[month >= policy_date]),
    .groups = "drop"
  ) %>%
  filter(
    n_pre_months >= min_pre_months,
    n_post_months >= min_post_months
  )

twfe_data_clean <- twfe_data %>%
  semi_join(
    airport_keep_twfe,
    by = c("airport_id", "airport_name")
  ) %>%
  # Exclude other Santa Clara County airports from the control group.
  filter(treated | airport_county_fips != "06085")

coverage_clean <- twfe_data_clean %>%
  group_by(airport_id, airport_name, treated) %>%
  summarise(
    n_pre_months = n_distinct(month[month < policy_date]),
    n_post_months = n_distinct(month[month >= policy_date]),
    min_month = min(month),
    max_month = max(month),
    .groups = "drop"
  ) %>%
  arrange(desc(treated), airport_name)

print(coverage_clean)

readr::write_csv(
  coverage_clean,
  file.path(output_folder, "twfe_airport_coverage_clean_2016_2024.csv")
)

readr::write_csv(
  twfe_data_clean,
  file.path(output_folder, "twfe_data_clean_0_4km_weather_2016_2024.csv")
)

# -------------------------------------------------------------------------
# 17. TWFE MODELS: BASELINE AND WEATHER CONTROLS
# -------------------------------------------------------------------------
twfe_model_baseline <- fixest::feols(
  log_lead ~ treated_post | airport_id + month,
  data = twfe_data_clean,
  cluster = ~ airport_id
)

twfe_model_weather <- fixest::feols(
  log_lead ~ treated_post + air_temp + wind_speed + rhum + pres |
    airport_id + month,
  data = twfe_data_clean,
  cluster = ~ airport_id
)

twfe_model_weather_winddir <- fixest::feols(
  log_lead ~ treated_post + air_temp + wind_speed + rhum + pres +
    wind_dir_sin + wind_dir_cos |
    airport_id + month,
  data = twfe_data_clean,
  cluster = ~ airport_id
)

print(summary(twfe_model_baseline))
print(summary(twfe_model_weather))
print(summary(twfe_model_weather_winddir))

fixest::etable(
  twfe_model_baseline,
  twfe_model_weather,
  twfe_model_weather_winddir,
  dict = c(
    treated_post = "RHV x Post-2022",
    air_temp = "Air temperature",
    wind_speed = "Wind speed",
    rhum = "Relative humidity",
    pres = "Surface pressure",
    wind_dir_sin = "Wind direction sine",
    wind_dir_cos = "Wind direction cosine"
  )
)

beta_baseline <- coef(twfe_model_baseline)["treated_post"]
beta_weather <- coef(twfe_model_weather)["treated_post"]
beta_weather_winddir <- coef(twfe_model_weather_winddir)["treated_post"]

percent_effects <- tibble(
  model = c("Baseline", "Weather controls", "Weather + wind direction"),
  beta = c(beta_baseline, beta_weather, beta_weather_winddir),
  approximate_percent_change = (exp(beta) - 1) * 100
)

print(percent_effects)

readr::write_csv(
  percent_effects,
  file.path(output_folder, "twfe_percent_effects_2016_2024.csv")
)

# -------------------------------------------------------------------------
# 18. PLOT RHV VS CONTROL AIRPORTS
# -------------------------------------------------------------------------
plot_data_twfe <- twfe_data_clean %>%
  mutate(group = if_else(treated, "Reid-Hillview", "Control airports")) %>%
  group_by(group, month) %>%
  summarise(
    mean_lead = mean(mean_lead, na.rm = TRUE),
    .groups = "drop"
  )

p_twfe_raw <- ggplot(
  plot_data_twfe,
  aes(x = month, y = mean_lead, color = group)
) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = policy_date, linetype = "dashed") +
  labs(
    title = "Monthly Ambient Lead Concentration: Reid-Hillview vs Control Airports",
    subtitle = "Distance band: 0-4 km; sample: 2016-01-01 to 2025-01-01",
    x = "Month",
    y = "Mean daily ambient lead concentration",
    color = "Group"
  ) +
  theme_minimal(base_size = 13)

print(p_twfe_raw)

ggsave(
  filename = file.path(fig_dir, "twfe_rhv_vs_controls_raw_0_4km_2016_2024.png"),
  plot = p_twfe_raw,
  width = 11,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 19. EVENT STUDY WITH WEATHER CONTROLS, REFERENCE YEAR = 2019
# -------------------------------------------------------------------------
event_data <- twfe_data_clean %>%
  mutate(
    event_year = lubridate::year(month),
    rel_year = event_year - 2022
  ) %>%
  filter(
    rel_year >= -6,   # 2016
    rel_year <= 2     # 2024
  )

# Reference year = 2019, so rel_year = -3.
event_model_ref2019_weather <- fixest::feols(
  log_lead ~ i(rel_year, treated, ref = -3) +
    air_temp + wind_speed + rhum + pres + wind_dir_sin + wind_dir_cos |
    airport_id + event_year,
  data = event_data,
  cluster = ~ airport_id
)

print(summary(event_model_ref2019_weather))

event_results <- broom::tidy(
  event_model_ref2019_weather,
  conf.int = TRUE
) %>%
  filter(str_detect(term, "rel_year::")) %>%
  mutate(
    rel_year = as.numeric(str_extract(term, "-?\\d+"))
  )

p_event <- ggplot(event_results, aes(x = rel_year, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.15
  ) +
  labs(
    title = "Event Study: Reid-Hillview Airport Lead Concentration",
    subtitle = "Distance band: 0-4 km; weather controls included; reference year is 2019",
    x = "Years relative to January 2022",
    y = "Estimated effect on log lead concentration"
  ) +
  theme_minimal(base_size = 13)

print(p_event)

ggsave(
  filename = file.path(fig_dir, "event_study_reid_hillview_0_4km_weather_ref2019_2016_2024.png"),
  plot = p_event,
  width = 9,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 20. WEATHER DIAGNOSTIC PLOTS
# -------------------------------------------------------------------------
p_diag_temp <- lead_airport_weather %>%
  filter(!is.na(air_temp)) %>%
  ggplot(aes(x = monitor_lat, y = air_temp)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "lm") +
  labs(
    title = "Weather Diagnostic: Temperature vs Latitude",
    subtitle = "Expected pattern: colder temperatures at higher latitudes",
    x = "Monitor latitude",
    y = "Air temperature (C)"
  ) +
  theme_minimal(base_size = 13)

print(p_diag_temp)

ggsave(
  file.path(fig_dir, "diagnostic_temperature_vs_latitude_2016_2024.png"),
  p_diag_temp,
  width = 8,
  height = 5,
  dpi = 300
)

p_diag_wind <- lead_airport_weather %>%
  filter(!is.na(wind_speed), lead_conc > 0) %>%
  ggplot(aes(x = wind_speed, y = lead_conc)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  labs(
    title = "Lead Concentration vs Wind Speed",
    subtitle = "Testing whether wind speed is related to measured lead concentration",
    x = "Wind speed (m/s)",
    y = "Lead concentration, log scale"
  ) +
  theme_minimal(base_size = 13)

print(p_diag_wind)

ggsave(
  file.path(fig_dir, "diagnostic_lead_vs_wind_speed_2016_2024.png"),
  p_diag_wind,
  width = 8,
  height = 5,
  dpi = 300
)

cat("\nDONE.\n")
cat("Outputs saved in:", output_folder, "\n")
cat("Figures saved in:", fig_dir, "\n")