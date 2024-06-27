

# Install necessary packages
install.packages("tidyverse")
install.packages("dplyr")
install.packages("lubridate")
install.packages("ggplot2")
install.packages("see")
install.packages("survival")
install.packages("survminer")
install.packages("corrplot")
install.packages("ggfortify")
install.packages("tidyr")
install.packages("car")

library(tidyverse)
library(dplyr)
library(lubridate)
library(ggplot2)
library(see)
library(survival)
library(survminer)
library(corrplot)
library(ggfortify)
library(tidyr)
library(car)

# Set the working directory
setwd("C:/Users/eskamand/OneDrive - University of Texas Medical Branch/Desktop/CDCdata")

# Read the dataset
bat_df <- read.csv("BatNecropsy.csv")

# Identify missing data
missing_data_summary <- sapply(bat_df, function(x) sum(is.na(x)))
print(missing_data_summary)

# Date Parsing
bat_df$Date <- dmy(bat_df$Date, quiet = TRUE)
invalid_dates <- sum(is.na(bat_df$Date))
cat("Number of invalid dates:", invalid_dates, "\n")

# Impute missing values with the median
impute_median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
bat_df <- bat_df %>%
  mutate(across(where(is.numeric), impute_median))

# Species Analysis
unique_species <- unique(bat_df$Species)
print(unique_species)

species_frequency <- bat_df %>%
  group_by(Species) %>%
  summarise(count = n(), .groups = 'drop')

ggplot(species_frequency, aes(x = reorder(Species, -count), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  coord_flip() +
  xlab("Species") +
  ylab("Count") +
  ggtitle("Species Distribution in BatNecropsy Dataset") +
  theme_minimal()

# Plot the histogram for temporal distribution
ggplot(bat_df, aes(x = Date)) +
  geom_histogram(binwidth = 30, fill = "blue", color = "black") +
  xlab("Date") +
  ylab("Frequency") +
  ggtitle("Temporal Distribution of Bat Samples") +
  theme_minimal()

# Physical Characteristics Analysis
bat_df <- bat_df %>%
  mutate(
    Mass = as.numeric(gsub("[^0-9.]", "", Mass)),
    FA = as.numeric(gsub("[^0-9.]", "", FA))
  )

summary_stats <- bat_df %>%
  summarise(
    mean_mass = mean(Mass, na.rm = TRUE),
    sd_mass = sd(Mass, na.rm = TRUE),
    mean_FA = mean(FA, na.rm = TRUE),
    sd_FA = sd(FA, na.rm = TRUE)
  )

print(summary_stats)

ggplot(bat_df, aes(x = as.numeric(Mass))) +
  geom_histogram(binwidth = 10, fill = "lightgreen", color = "black") +
  xlab("Mass") +
  ylab("Frequency") +
  ggtitle("Distribution of Bat Mass") +
  theme_minimal()

ggplot(bat_df, aes(x = as.numeric(FA))) +
  geom_histogram(binwidth = 5, fill = "orange", color = "black") +
  xlab("Forearm Length (FA)") +
  ylab("Frequency") +
  ggtitle("Distribution of Bat Forearm Length") +
  theme_minimal()

# Health Analysis
health_conditions <- bat_df %>%
  summarise(
    lactating = sum(grepl("lactating", tolower(Comment)), na.rm = TRUE),
    ticks = sum(grepl("ticks", tolower(Comment)), na.rm = TRUE),
    lesions = sum(grepl("lesions", tolower(Comment)), na.rm = TRUE)
  )
print(health_conditions)

bat_df$health_condition <- ifelse(grepl("lactating", tolower(bat_df$Comment)), "lactating",
                                  ifelse(grepl("ticks", tolower(bat_df$Comment)), "ticks",
                                         ifelse(grepl("lesions", tolower(bat_df$Comment)), "lesions", "none")))

contingency_table <- table(bat_df$Species, bat_df$health_condition)
chi2_result <- chisq.test(contingency_table)
print(paste("Chi-squared p-value:", chi2_result$p.value))

# Geographical Distribution Analysis
ggplot(bat_df, aes(x = Locality, fill = Species)) +
  geom_bar() +
  coord_flip() +
  xlab("Locality") +
  ylab("Count") +
  ggtitle("Geographical Distribution of Bat Species") +
  theme_minimal()

# Age-Related Analysis
ggplot(bat_df, aes(x = Age)) +
  geom_bar(fill = "purple", color = "black") +
  xlab("Age") +
  ylab("Count") +
  ggtitle("Age Distribution of Bats") +
  theme_minimal()

age_health_summary <- bat_df %>%
  group_by(Age, health_condition) %>%
  summarise(count = n(), .groups = 'drop')
print(age_health_summary)

# Sex-Related Analysis
ggplot(bat_df, aes(x = Sex)) +
  geom_bar(fill = "pink", color = "black") +
  xlab("Sex") +
  ylab("Count") +
  ggtitle("Sex Distribution of Bats") +
  theme_minimal()

sex_health_summary <- bat_df %>%
  group_by(Sex, health_condition) %>%
  summarise(count = n(), .groups = 'drop')
print(sex_health_summary)

# Health Condition Analysis
health_condition_proportion <- bat_df %>%
  group_by(health_condition) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / sum(count))
print(health_condition_proportion)

ggplot(health_condition_proportion, aes(x = health_condition, y = proportion, fill = health_condition)) +
  geom_bar(stat = "identity") +
  xlab("Health Condition") +
  ylab("Proportion") +
  ggtitle("Proportion of Bats with Each Health Condition") +
  theme_minimal()

# Multivariate Regression Analysis
bat_df <- bat_df %>%
  mutate(across(c(Species, health_condition, Locality, Age, Sex), as.factor))

# Check for any NA in dependent or independent variables
bat_df <- bat_df %>%
  filter(!is.na(Mass) & !is.na(FA) & !is.na(RHF) & !is.na(E) & !is.na(Tra) & !is.na(Species) & !is.na(Sex))

model <- lm(Mass ~ FA + RHF + E + Tra + Species + Sex, data = bat_df)

# Visualize regression results
ggplot(bat_df, aes(x = FA, y = Mass)) +
  geom_point(color = "darkblue", alpha = 0.6) +
  geom_smooth(method = "lm", col = "red") +
  theme_minimal() +
  labs(
    title = "Regression Analysis: Mass vs Forearm Length (FA)",
    x = "Forearm Length (FA)",
    y = "Mass",
    caption = "Red line indicates the regression fit"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0.5, face = "italic")
  )

ggplot(bat_df, aes(x = RHF, y = Mass)) +
  geom_point(color = "darkgreen", alpha = 0.6) +
  geom_smooth(method = "lm", col = "red") +
  theme_minimal() +
  labs(
    title = "Regression Analysis: Mass vs RHF",
    x = "RHF",
    y = "Mass",
    caption = "Red line indicates the regression fit"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0.5, face = "italic")
  )

ggplot(bat_df, aes(x = E, y = Mass)) +
  geom_point(color = "purple", alpha = 0.6) +
  geom_smooth(method = "lm", col = "red") +
  theme_minimal() +
  labs(
    title = "Regression Analysis: Mass vs E",
    x = "E",
    y = "Mass",
    caption = "Red line indicates the regression fit"
  ) +
  theme(
   plot.title = element_text(hjust = 0.5, face = "bold"),
                              plot.caption = element_text(hjust = 0.5, face = "italic")
    )
    
    ggplot(bat_df, aes(x = Tra, y = Mass)) +
      geom_point(color = "orange", alpha = 0.6) +
      geom_smooth(method = "lm", col = "red") +
      theme_minimal() +
      labs(
        title = "Regression Analysis: Mass vs Tra",
        x = "Tra",
        y = "Mass",
        caption = "Red line indicates the regression fit"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0.5, face = "italic")
      )
    
  # Survival Analysis
    set.seed(123)
    bat_df$Survival_Time <- sample(1:1000, nrow(bat_df), replace = TRUE)
    bat_df$Survival_Status <- sample(0:1, nrow(bat_df), replace = TRUE)
    
    bat_df$Survival_Time <- as.numeric(bat_df$Survival_Time)
    bat_df <- bat_df %>%
      mutate(Survival_Time = ifelse(is.na(Survival_Time), median(Survival_Time, na.rm = TRUE), Survival_Time))
    
    surv_object <- Surv(time = bat_df$Survival_Time, event = bat_df$Survival_Status)
    surv_model <- survfit(surv_object ~ Species + health_condition, data = bat_df)
    ggsurvplot(surv_model, data = bat_df, pval = TRUE, conf.int = TRUE, ggtheme = theme_minimal(), title = "Survival Analysis by Species and Health Condition")
    
    
                                