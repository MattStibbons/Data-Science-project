# 1. LOAD LIBRARIES 
library(readr)
library(dplyr)
library(ggplot2)

# 2. IMPORT DATA 

daily_LEAD_2016 <- read_csv("raw data/daily_LEAD_2016.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2017 <- read_csv("raw data/daily_LEAD_2017.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2018 <- read_csv("raw data/daily_LEAD_2018.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2019 <- read_csv("raw data/daily_LEAD_2019.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2020 <- read_csv("raw data/daily_LEAD_2020.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2021 <- read_csv("raw data/daily_LEAD_2021.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2022 <- read_csv("raw data/daily_LEAD_2022.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2023 <- read_csv("raw data/daily_LEAD_2023.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2024 <- read_csv("raw data/daily_LEAD_2024.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2025 <- read_csv("raw data/daily_LEAD_2025.csv", col_types = cols(`Method Code` = col_character()))

# 3. COMBINE ALL YEARS
all_data <- bind_rows(
  daily_LEAD_2016, daily_LEAD_2017, daily_LEAD_2018, daily_LEAD_2019, 
  daily_LEAD_2020, daily_LEAD_2021, daily_LEAD_2022, daily_LEAD_2023, 
  daily_LEAD_2024, daily_LEAD_2025
)

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
  geom_line() +
  scale_y_log10() +
  labs(
    title = "Monthly Log Arithmetic Mean (Airport Sites)",
    x = "Month",
    y = "Arithmetic Mean (log scale)"
  ) +
  theme_minimal()

# 6. EXPLORATION
unique(all_data$`County Name`)

