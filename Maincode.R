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
