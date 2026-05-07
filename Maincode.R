
###########################################################################
# LEAD STUDY — EXPANDING RADIUS ROBUSTNESS
# Core idea: Does the RHV treatment effect hold as we widen the control
# group from monitors within 1km → 2km → 3km of any US airport?
###########################################################################

# 1. LIBRARIES
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, purrr, ncdf4, lubridate, sf, fixest, ggplot2)

# 2. CONFIG
years          <- 2000:2025
raw_data_folder <- "raw data/"

###########################################################################
# 3. LOAD LEAD DATA
###########################################################################

all_lead_data <- map_dfr(years, function(y) {
  f <- paste0(raw_data_folder, "daily_LEAD_", y, ".csv")
  if (file.exists(f)) read_csv(f, col_types = cols(`Method Code` = col_character()))
})

lead_base <- all_lead_data %>%
  filter(`Arithmetic Mean` > 0) %>%
  mutate(
    date      = as.Date(`Date Local`),
    latitude  = Latitude,
    longitude = Longitude
  )

###########################################################################
# 4. WEATHER MERGE
###########################################################################

final_list <- list()

for (z in years) {
  
  cat("Processing year:", z, "\n")
  obs_year <- lead_base %>% filter(year(date) == z)
  if (nrow(obs_year) == 0) next
  
  nc_ref   <- nc_open(paste0(raw_data_folder, "air.sig995.", z, ".nc"))
  lon_vals <- ncvar_get(nc_ref, "lon")
  lat_vals <- ncvar_get(nc_ref, "lat")
  nc_close(nc_ref)
  
  obs_year <- obs_year %>%
    mutate(
      lon_360 = ifelse(longitude < 0, longitude + 360, longitude),
      day_idx = as.numeric(format(date, "%j"))
    )
  
  lat_idx <- map_int(obs_year$latitude, ~which.min(abs(lat_vals - .x)))
  lon_idx <- map_int(obs_year$lon_360, ~which.min(abs(lon_vals - .x)))
  
  pull_weather <- function(var, prefix) {
    f <- paste0(raw_data_folder, prefix, ".", z, ".nc")
    if (!file.exists(f)) return(rep(NA, nrow(obs_year)))
    nc   <- nc_open(f)
    vals <- map_dbl(1:nrow(obs_year), function(i) {
      ncvar_get(nc, var,
                start = c(lon_idx[i], lat_idx[i], obs_year$day_idx[i]),
                count = c(1, 1, 1))
    })
    nc_close(nc)
    vals
  }
  
  obs_year$temp_k <- pull_weather("air",  "air.sig995")
  obs_year$uwnd   <- pull_weather("uwnd", "uwnd.sig995")
  obs_year$vwnd   <- pull_weather("vwnd", "vwnd.sig995")
  obs_year$rhum   <- pull_weather("rhum", "rhum.sig995")
  obs_year$pres   <- pull_weather("pres", "pres.sfc") / 1000
  obs_year$pr_wtr <- pull_weather("pr_wtr", "pr_wtr.eatm")
  
  obs_year <- obs_year %>%
    mutate(
      air_temp   = temp_k - 273.15,
      wind_speed = sqrt(uwnd^2 + vwnd^2)
    )
  
  final_list[[as.character(z)]] <- obs_year
}

final_dataset <- bind_rows(final_list)

###########################################################################
# 5. AIRPORT DISTANCES
###########################################################################

airport_file <- paste0(raw_data_folder, "airports.csv")
if (!file.exists(airport_file)) {
  download.file(
    "https://davidmegginson.github.io/ourairports-data/airports.csv",
    airport_file, mode = "wb"
  )
}

airports <- read_csv(airport_file, show_col_types = FALSE) %>%
  filter(iso_country == "US",
         type %in% c("small_airport", "medium_airport", "large_airport")) %>%
  st_as_sf(coords = c("longitude_deg", "latitude_deg"), crs = 4326)

monitor_locs <- lead_base %>%
  distinct(latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

monitors_proj <- st_transform(monitor_locs, 5070)
airports_proj <- st_transform(airports,     5070)

dist_matrix          <- st_distance(monitors_proj, airports_proj)
monitor_locs$dist_km <- as.numeric(apply(dist_matrix, 1, min)) / 1000

# Diagnostic — see how many monitors fall in each radius
cat("\n--- Monitor counts by radius ---\n")
cat("Within 1km:", sum(monitor_locs$dist_km <= 1), "\n")
cat("Within 2km:", sum(monitor_locs$dist_km <= 2), "\n")
cat("Within 3km:", sum(monitor_locs$dist_km <= 3), "\n")
cat("Within 4km:", sum(monitor_locs$dist_km <= 4), "\n")
cat("Within 5km:", sum(monitor_locs$dist_km <= 5), "\n")

###########################################################################
# 6. BUILD ANALYSIS DATASET
###########################################################################

analysis_data <- final_dataset %>%
  left_join(
    monitor_locs %>%
      st_drop_geometry() %>%
      select(latitude, longitude, dist_km),
    by = c("latitude", "longitude")
  ) %>%
  mutate(
    year_month        = floor_date(date, "month"),
    year_to_treatment = year(date) - 2022,
    is_rhv            = grepl("Reid Hillview", `Local Site Name`, ignore.case = TRUE),
    is_airport_name   = grepl("airport",       `Local Site Name`, ignore.case = TRUE),
    is_treated        = is_rhv & date >= as.Date("2022-01-01")
  )

# Quick check: confirm RHV monitor distance
cat("\n--- RHV monitor distance to nearest airport ---\n")
analysis_data %>%
  filter(is_rhv) %>%
  distinct(`Local Site Name`, dist_km) %>%
  print()
###########################################################################
# MAIN MODEL: AIRPORT NAME SPECIFICATION
# Monitors with "airport" in their name — most directly comparable to RHV
# This is the headline result; distance specs below are robustness checks
###########################################################################

model_airport_name <- feols(
  log(`Arithmetic Mean`) ~ is_treated + air_temp + wind_speed + rhum + pres |
    `Local Site Name` + year_month,
  data = analysis_data %>% filter(is_airport_name == TRUE)
)

cat("\n========== MAIN MODEL: AIRPORT NAME ==========\n")
summary(model_airport_name)

# Event study — airport name
ev_airport <- analysis_data %>%
  filter(is_airport_name == TRUE,
         year(date) >= 2016,
         year(date) <= 2025)

event_airport_name <- feols(
  log(`Arithmetic Mean`) ~
    i(year_to_treatment, is_rhv, ref = -1) +
    air_temp + wind_speed + rhum + pres |
    `Local Site Name` + year_month,
  data = ev_airport
)

iplot(
  event_airport_name,
  main    = "Event Study: Airport Name Specification (MAIN)",
  xlab    = "Years Relative to 2022 (0 = Ban)",
  ylab    = "Effect on log Lead",
  pt.join = TRUE
)

# Add airport-name to the combined robustness table for comparison
airport_name_row <- tibble(
  radius_km  = NA_real_,
  n_monitors = ev_airport %>% distinct(`Local Site Name`) %>% nrow(),
  n_controls = ev_airport %>% filter(!is_rhv) %>% distinct(`Local Site Name`) %>% nrow(),
  n_obs      = nrow(analysis_data %>% filter(is_airport_name == TRUE)),
  estimate   = coef(model_airport_name)["is_treatedTRUE"],
  std_error  = se(model_airport_name)["is_treatedTRUE"],
  ci_low     = estimate - 1.96 * std_error,
  ci_high    = estimate + 1.96 * std_error,
  label      = "Airport name"
)
###########################################################################
# 7. RADIUS ROBUSTNESS — DiD
#
# For each radius cutoff, run the same DiD spec.
# Sample = all monitors within X km of any US airport.
# Treatment = RHV post-2022. Controls = all other monitors in the sample.
# The coefficient on is_treated tells us:
#   "Did RHV fall more than same-radius control monitors after 2022?"
###########################################################################

radii <- c(1, 2, 3, 4, 5)     # km cutoffs to loop over

did_results <- map_dfr(radii, function(r) {
  
  df <- analysis_data %>% filter(!is.na(dist_km), dist_km <= r)
  
  # Skip if RHV not in sample (shouldn't happen but safety check)
  if (!any(df$is_rhv)) {
    cat("WARNING: No RHV observations within", r, "km — skipping\n")
    return(NULL)
  }
  
  n_monitors  <- df %>% distinct(`Local Site Name`) %>% nrow()
  n_controls  <- df %>% filter(!is_rhv) %>% distinct(`Local Site Name`) %>% nrow()
  n_obs       <- nrow(df)
  
  cat("\n--- Radius:", r, "km | Monitors:", n_monitors,
      "| Controls:", n_controls, "| Obs:", n_obs, "---\n")
  
  mod <- feols(
    log(`Arithmetic Mean`) ~ is_treated + air_temp + wind_speed + rhum + pres |
      `Local Site Name` + year_month,
    data = df
  )
  
  coefs <- coef(mod)
  ses   <- se(mod)
  
  tibble(
    radius_km   = r,
    n_monitors  = n_monitors,
    n_controls  = n_controls,
    n_obs       = n_obs,
    estimate    = coefs["is_treatedTRUE"],
    std_error   = ses["is_treatedTRUE"],
    ci_low      = estimate - 1.96 * std_error,
    ci_high     = estimate + 1.96 * std_error,
    label       = paste0("≤", r, "km radius")   # <-- add this line
  )
})

cat("\n========== DiD ROBUSTNESS TABLE ==========\n")
print(did_results)
# Combine main model + robustness rows into one summary
all_results <- bind_rows(airport_name_row, did_results) %>%
  mutate(label = ifelse(is.na(label), paste0("≤", radius_km, "km radius"), label),
         label = factor(label, levels = c("Airport name",
                                          "≤1km radius", "≤2km radius",
                                          "≤3km radius", "≤5km radius")))

cat("\n========== FULL RESULTS TABLE ==========\n")
print(all_results %>% select(label, n_controls, n_obs, estimate, std_error, ci_low, ci_high))

ggplot(all_results, aes(x = label, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8) +
  geom_point(aes(colour = label == "Airport name"), size = 3.5) +
  geom_text(aes(label = paste0("n=", n_controls)),
            vjust = -1.2, size = 3.2) +
  scale_colour_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                      guide = "none") +
  labs(
    title    = "Treatment Effect Across All Specifications",
    subtitle = "Red = main model (airport name); Blue = distance robustness checks",
    x        = NULL,
    y        = "Estimated effect on log Lead",
    caption  = "95% CIs. All specs: monitor + year-month FEs, weather controls."
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
###########################################################################
# 8. COEFFICIENT STABILITY PLOT — DiD
###########################################################################

ggplot(did_results, aes(x = factor(radius_km), y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8) +
  geom_point(size = 3.5, colour = "steelblue") +
  geom_text(aes(label = paste0("n=", n_controls, " controls")),
            vjust = -1.2, size = 3.2) +
  scale_x_discrete(labels = function(x) paste0("≤", x, " km")) +
  labs(
    title    = "DiD Coefficient by Control Group Radius",
    subtitle = "Treatment: Reid Hillview post-Jan 2022 | Dep. var: log(lead)",
    x        = "Sample radius from nearest airport",
    y        = "Estimated treatment effect",
    caption  = "95% CIs shown. All specs include monitor + year-month FEs."
  ) +
  theme_minimal(base_size = 13)

###########################################################################
# 9. RADIUS ROBUSTNESS — EVENT STUDIES
#
# Same expanding-radius logic but now we show the full year-by-year path.
# This lets you check BOTH coefficient stability AND parallel pre-trends
# across different control groups simultaneously.
###########################################################################

event_models <- map(radii, function(r) {
  
  df <- analysis_data %>%
    filter(!is.na(dist_km),
           dist_km <= r,
           year(date) >= 2016,
           year(date) <= 2025)
  
  if (!any(df$is_rhv)) return(NULL)
  
  feols(
    log(`Arithmetic Mean`) ~
      i(year_to_treatment, is_rhv, ref = -1) +
      air_temp + wind_speed + rhum + pres |
      `Local Site Name` + year_month,
    data = df
  )
})
names(event_models) <- paste0("<=", radii, "km")

# Remove any NULLs
event_models <- Filter(Negate(is.null), event_models)

###########################################################################
# 10. EVENT STUDY PLOTS — ONE PER RADIUS (individual, easy to read)
###########################################################################

for (nm in names(event_models)) {
  iplot(
    event_models[[nm]],
    main    = paste0("Event Study: ", nm, " radius"),
    xlab    = "Years Relative to 2022 (0 = Ban)",
    ylab    = "Effect on log Lead",
    pt.join = TRUE
  )
}

###########################################################################
# 11. COMBINED EVENT STUDY PLOT (all radii overlaid, one clean figure)
#
# Extracts coefficients manually so we can overlay on a single ggplot.
# Easier to see whether the event path shifts as the control group widens.
###########################################################################

event_coefs <- map_dfr(names(event_models), function(nm) {
  mod <- event_models[[nm]]
  
  # Pull only the year_to_treatment:is_rhv interactions
  cf  <- coef(mod)
  se_ <- se(mod)
  
  idx <- grep("year_to_treatment::", names(cf))
  
  tibble(
    radius    = nm,
    term      = names(cf)[idx],
    estimate  = cf[idx],
    std_error = se_[idx],
    ci_low    = estimate - 1.96 * std_error,
    ci_high   = estimate + 1.96 * std_error
  ) %>%
    mutate(
      # Extract the numeric year-relative-to-treatment value from term name
      year_rel = as.numeric(gsub(".*year_to_treatment::([-0-9]+):.*", "\\1", term))
    )
})

# Add the ref year manually (coefficient = 0 by construction)
ref_rows <- tibble(
  radius    = names(event_models),
  term      = "ref",
  estimate  = 0,
  std_error = 0,
  ci_low    = 0,
  ci_high   = 0,
  year_rel  = -1
)

event_coefs <- bind_rows(event_coefs, ref_rows) %>%
  arrange(radius, year_rel)

ggplot(event_coefs, aes(x = year_rel, y = estimate,
                        colour = radius, fill = radius)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey30") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  annotate("text", x = -0.5, y = Inf, label = "Closure →",
           hjust = 1.1, vjust = 1.5, size = 3.5, colour = "grey30") +
  scale_x_continuous(breaks = sort(unique(event_coefs$year_rel))) +
  labs(
    title    = "Event Study: RHV Treatment Effect by Control Group Radius",
    subtitle = "Each line uses a different radius cutoff for sample selection",
    x        = "Years Relative to 2022",
    y        = "Effect on log Lead Concentration",
    colour   = "Sample radius",
    fill     = "Sample radius",
    caption  = "Shaded bands = 95% CI. Ref year = -1. Monitor + year-month FEs."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

