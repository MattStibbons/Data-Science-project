unique(daily_LEAD_2022$`County Name`)
unique(daily_LEAD_2022$`State Name`)

colnames(daily_LEAD_2022)
library(dplyr)

filtered_data <- daily_LEAD_2022 %>%
  filter(grepl("airport", `Local Site Name`, ignore.case = TRUE))

library(ggplot2)
library(dplyr)

filtered_data %>%
  mutate(`Date Local` = as.Date(`Date Local`)) %>%
  ggplot(aes(x = `Date Local`, y = `Arithmetic Mean`, color = `Local Site Name`)) +
  geom_line() +
  labs(
    title = "Arithmetic Mean Over Time by Site",
    x = "Date",
    y = "Arithmetic Mean"
  ) +
  theme_minimal()







unique(daily_LEAD_2018$`County Name`)
unique(daily_LEAD_2018$`State Name`)

colnames(daily_LEAD_2018)
library(dplyr)

filtered_data <- daily_LEAD_2018 %>%
  filter(grepl("airport", `Local Site Name`, ignore.case = TRUE))

library(ggplot2)
library(dplyr)

filtered_data %>%
  mutate(`Date Local` = as.Date(`Date Local`)) %>%
  ggplot(aes(x = `Date Local`, y = `Arithmetic Mean`, color = `Local Site Name`)) +
  geom_line() +
  labs(
    title = "Arithmetic Mean Over Time by Site",
    x = "Date",
    y = "Arithmetic Mean"
  ) +
  theme_minimal()



filtered_data %>%
  mutate(`Date Local` = as.Date(`Date Local`)) %>%
  ggplot(aes(x = `Date Local`, y = `Arithmetic Mean`, color = `Local Site Name`)) +
  geom_line() +
  scale_y_log10() +
  labs(
    title = "Log-Scaled Arithmetic Mean Over Time by Site",
    x = "Date",
    y = "Arithmetic Mean (log scale)"
  ) +
  theme_minimal()




daily_LEAD_2016 <- read_csv("daily_LEAD_2016.csv")
daily_LEAD_2017 <- read_csv("daily_LEAD_2017.csv")
daily_LEAD_2018 <- read_csv("daily_LEAD_2018.csv")
daily_LEAD_2019 <- read_csv("daily_LEAD_2019.csv")
daily_LEAD_2020 <- read_csv("daily_LEAD_2020.csv")
daily_LEAD_2021 <- read_csv("daily_LEAD_2021.csv")
daily_LEAD_2022 <- read_csv("daily_LEAD_2022.csv")
daily_LEAD_2023 <- read_csv("daily_LEAD_2023.csv")
daily_LEAD_2024 <- read_csv("daily_LEAD_2024.csv")
daily_LEAD_2025 <- read_csv("daily_LEAD_2025.csv")

library(readr)

daily_LEAD_2016 <- read_csv("daily_LEAD_2016.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2017 <- read_csv("daily_LEAD_2017.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2018 <- read_csv("daily_LEAD_2018.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2019 <- read_csv("daily_LEAD_2019.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2020 <- read_csv("daily_LEAD_2020.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2021 <- read_csv("daily_LEAD_2021.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2022 <- read_csv("daily_LEAD_2022.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2023 <- read_csv("daily_LEAD_2023.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2024 <- read_csv("daily_LEAD_2024.csv", col_types = cols(`Method Code` = col_character()))
daily_LEAD_2025 <- read_csv("daily_LEAD_2025.csv", col_types = cols(`Method Code` = col_character()))

library(dplyr)
library(readr)
library(ggplot2)

all_data <- bind_rows(
  daily_LEAD_2016,
  daily_LEAD_2017,
  daily_LEAD_2018,
  daily_LEAD_2019,
  daily_LEAD_2020,
  daily_LEAD_2021,
  daily_LEAD_2022,
  daily_LEAD_2023,
  daily_LEAD_2024,
  daily_LEAD_2025
)


filtered_data <- all_data %>%
  filter(grepl("airport", `Local Site Name`, ignore.case = TRUE))
filtered_data <- filtered_data %>%
  filter(`Arithmetic Mean` > 0)

filtered_data %>%
  mutate(`Date Local` = as.Date(`Date Local`)) %>%
  ggplot(aes(x = `Date Local`, y = `Arithmetic Mean`, color = `Local Site Name`)) +
  geom_line(alpha = 0.7) +
  scale_y_log10() +
  labs(
    title = "Log-Scaled Arithmetic Mean Over Time (Airport Sites)",
    x = "Date",
    y = "Arithmetic Mean (log scale)"
  ) +
  theme_minimal()


library(dplyr)
library(ggplot2)

filtered_data %>%
  mutate(`Date Local` = as.Date(`Date Local`)) %>%
  filter(`Arithmetic Mean` > 0) %>%
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


unique(all_data$`County Name`)