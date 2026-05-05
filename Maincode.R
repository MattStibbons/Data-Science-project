###########################################################################
# LEAD & WEATHER SPATIAL-TEMPORAL MERGE (2000-2025)
# Methodology: NCEP Reanalysis 1 Sigma 995 Mapping
###########################################################################

# 1. LOAD LIBRARIES
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, ggplot2, purrr, ncdf4, tibble, lubridate)

# 2. CONFIGURATION
years <- 2000:2025
raw_data_folder <- "raw data/" # Ensure this folder contains your .nc and .csv files

# 3. IMPORT & CLEAN LEAD DATA
all_lead_data <- map_dfr(years, function(y) {
  file_path <- paste0(raw_data_folder, "daily_LEAD_", y, ".csv")
  if (file.exists(file_path)) {
    read_csv(file_path, col_types = cols(`Method Code` = col_character()))
  } else {
    return(NULL)
  }
})

# Filter for all airports (Treatment + Controls)
lead_base <- all_lead_data %>%
  filter(grepl("airport", `Local Site Name`, ignore.case = TRUE)) %>%
  filter(`Arithmetic Mean` > 0) %>%
  mutate(date = as.Date(`Date Local`)) %>%
  rename(latitude = Latitude, longitude = Longitude)

# 4. WEATHER EXTRACTION LOOP
final_list <- list()

for (z in years) {
  cat("\014") # Clear console
  print(paste("EXTRACTING WEATHER FOR YEAR:", z))
  
  # Filter observations for current year
  obs_year <- lead_base %>% filter(format(date, "%Y") == as.character(z))
  if(nrow(obs_year) == 0) next
  
  # Prepare NC coordinate mapping
  # NCEP Longitude is 0-360; Time index is Day of Year
  nc_ref <- nc_open(paste0(raw_data_folder, 'air.sig995.', z, '.nc'))
  lon_vals <- ncvar_get(nc_ref, "lon") 
  lat_vals <- ncvar_get(nc_ref, "lat")
  nc_close(nc_ref)
  
  obs_year <- obs_year %>%
    mutate(
      lon_360 = ifelse(longitude < 0, longitude + 360, longitude),
      day_idx = as.numeric(format(date, "%j"))
    )
  
  # Find nearest grid indices (Saves us from the "North Pole" default error)
  lat_idx <- map_int(obs_year$latitude, ~which.min(abs(lat_vals - .x)))
  lon_idx <- map_int(obs_year$lon_360, ~which.min(abs(lon_vals - .x)))
  
  # Surgical Extraction Function
  pull_weather <- function(var_name, file_prefix) {
    f_path <- paste0(raw_data_folder, file_prefix, ".", z, ".nc")
    if(!file.exists(f_path)) return(rep(NA, nrow(obs_year)))
    
    nc <- nc_open(f_path)
    vals <- map_dbl(1:nrow(obs_year), function(i) {
      ncvar_get(nc, var_name, 
                start = c(lon_idx[i], lat_idx[i], obs_year$day_idx[i]), 
                count = c(1, 1, 1))
    })
    nc_close(nc)
    return(vals)
  }
  
  # Map variables to dataframe
  obs_year$temp_k <- pull_weather("air", "air.sig995")
  obs_year$uwnd   <- pull_weather("uwnd", "uwnd.sig995")
  obs_year$vwnd   <- pull_weather("vwnd", "vwnd.sig995")
  obs_year$rhum   <- pull_weather("rhum", "rhum.sig995")
  obs_year$pres   <- pull_weather("pres", "pres.sfc") / 1000 # Pa to kPa
  
  # Calculate Physics
  obs_year <- obs_year %>%
    mutate(
      air_temp = temp_k - 273.15,
      wind_speed = sqrt(uwnd^2 + vwnd^2),
      wind_dir = (atan2(-uwnd, -vwnd) * 180 / pi + 360) %% 360
    )
  
  final_list[[as.character(z)]] <- obs_year
}

# 5. FINAL MERGE AND SAVE
final_dataset <- bind_rows(final_list)
saveRDS(final_dataset, "final_analysis_data.rds")
write_csv(final_dataset, "airport_lead_weather_2000_2025.csv")

###########################################################################
# DIAGNOSTIC SUITE: CHECK IF IT WORKED
###########################################################################

cat("\n--- RUNNING DIAGNOSTICS ---\n")


temp_summary <- summary(final_dataset$air_temp)
print("Temperature Distribution (Celsius):")
print(temp_summary)

if(temp_summary["Mean"] < 0) {
  warning("CRITICAL: Average temperature is negative. Mapping likely failed!")
} else {
  print("SUCCESS: Average temperature is positive/realistic.")
}

# 2. Visual Check: Temperature vs. Latitude
# Purpose: Hotter airports should be at lower latitudes. 
diag_plot1 <- ggplot(final_dataset, aes(x = latitude, y = air_temp)) +
  geom_point(alpha = 0.2, color = "steelblue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Diagnostic 1: Temperature vs Latitude",
       subtitle = "Expected: Negative slope (Colder as you move North)") +
  theme_minimal()
print(diag_plot1)

# 3. Visual Check: Seasonal Swing (At Reid Hillview)
# Purpose: Check if summer is hotter than winter.
sample_airport <- final_dataset$`Local Site Name`[1]
diag_plot2 <- final_dataset %>%
  filter(`Local Site Name` == "Reid Hillview Airport", year(date) == 2018) %>%
  ggplot(aes(x = date, y = air_temp)) +
  geom_line() +
  labs(title = paste("Diagnostic 2: Seasonal Swing at", sample_airport),
       subtitle = "Expected: A clear 'wave' pattern peaking in July/Aug") +
  theme_minimal()
print(diag_plot2)

# 4. Wind Rose distribution
diag_plot3 <- ggplot(final_dataset, aes(x = wind_dir)) +
  geom_histogram(binwidth = 30, fill = "darkgreen", color = "white") +
  coord_polar(start = 0) +
  labs(title = "Diagnostic 3: Wind Direction (0-360)",
       subtitle = "Ensure directions are spread across the circle") +
  theme_minimal()
print(diag_plot3)

cat("\n--- DIAGNOSTICS COMPLETE ---\n")

###########################################################################
# LEAD TREND VISUALIZATIONS (2000-2025)
###########################################################################

# 1. Monthly Trends - Log Scale
# This helps identify seasonal spikes in lead levels across all monitored airports
final_dataset %>%
  mutate(month = format(date, "%Y-%m")) %>%
  group_by(`Local Site Name`, month) %>%
  summarise(mean_lead = mean(`Arithmetic Mean`, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = as.Date(paste0(month, "-01")), y = mean_lead, color = `Local Site Name`)) +
  geom_line(alpha = 0.6) +
  scale_y_log10() +
  labs(
    title = "Monthly Log Arithmetic Mean Lead (2000–2025)",
    subtitle = "Aggregated across all airport monitor sites",
    x = "Timeline",
    y = "Lead Concentration (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# 2. Yearly Trends - Long-term Comparison
# Useful for the Difference-in-Differences visual baseline
final_dataset %>%
  mutate(year = as.numeric(format(date, "%Y"))) %>%
  group_by(`Local Site Name`, year) %>%
  summarise(mean_lead = mean(`Arithmetic Mean`, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = mean_lead, color = `Local Site Name`)) +
  geom_line(size = 1, alpha = 0.8) +
  scale_y_log10() +
  labs(
    title = "Yearly Log Arithmetic Mean Lead (2000–2025)",
    x = "Year",
    y = "Lead Concentration (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Hidden due to high number of sites

# 3. New Diagnostic: Lead vs Wind Speed
# Since you have the weather data now, let's see if wind actually dilutes lead
final_dataset %>%
  ggplot(aes(x = wind_speed, y = `Arithmetic Mean`)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", color = "red") +
  scale_y_log10() +
  labs(
    title = "Lead Concentration vs. Wind Speed",
    subtitle = "Testing atmospheric dilution hypothesis",
    x = "Wind Speed (m/s)",
    y = "Lead Concentration (log scale)"
  ) +
  theme_minimal()


# Plotting Lead Concentration at US Airports on a map

library(ggplot2)
library(dplyr)
library(maps)

# Get US map polygons
us_map <- map_data("state")

# Aggregate your data (so points aren't overplotted)
plot_data <- final_dataset %>%
  group_by(latitude, longitude) %>%
  summarise(mean_lead = mean(`Arithmetic Mean`, na.rm = TRUE), .groups = "drop")

ggplot() +
  # Draw US map
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "white"
  ) +
  
  # Overlay your data
  geom_point(
    data = plot_data,
    aes(x = longitude, y = latitude, color = mean_lead),
    size = 2,
    alpha = 0.8
  ) +
  
  scale_color_viridis_c(trans = "log10") +
  
  coord_fixed(1.3) +
  
  labs(
    title = "Lead Concentration at US Airports (2000–2025)",
    subtitle = "Mapped across monitoring locations",
    x = "",
    y = "",
    color = "Lead (log scale)"
  ) +
  
  theme_minimal()
