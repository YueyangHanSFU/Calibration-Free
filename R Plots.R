# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

criteria_labels <- c(
  "Selection of most efficacious dose %",
  "Patients allocation to most efficacious dose %",
  "Selection of non-inferior doses %",
  "Patients allocation to non-inferior doses %",
  "Selection of overdose %",
  "Patients allocation to overdose %"
)

# Create data for Main Simulation (Scenarios 1-5 from Table 1)
main_data_part1 <- data.frame(
  Scenario = rep(1:5, each = 6),
  Criteria = rep(criteria_labels, 5),
  CFH_Design = c(70.6, 39.5, 94.6, 66.0, 0.0, 0.0,    # Scenario 1
                 44.2, 27.3, 80.0, 55.9, 0.0, 0.0,    # Scenario 2
                 27.2, 19.7, 94.4, 77.6, 0.0, 0.0,    # Scenario 3
                 45.4, 27.7, 70.6, 51.9, 16.2, 14.3,  # Scenario 4
                 48.8, 30.6, 78.8, 57.8, 11.4, 11.1), # Scenario 5
  Thompson = c(80.2, 43.9, 99.6, 74.8, 0.0, 0.0,      # Scenario 1
               28.6, 8.8, 83.2, 35.3, 0.0, 0.0,       # Scenario 2
               9.4, 3.0, 98.0, 69.8, 0.0, 0.0,        # Scenario 3
               49.2, 21.8, 73.2, 48.8, 22.4, 7.7,     # Scenario 4
               53.2, 23.7, 81.8, 54.8, 15.0, 5.2),    # Scenario 5
  TwoStage = c(50.8, 12.4, 82.6, 30.1, 0.0, 0.0,      # Scenario 1
               9.4, 1.7, 30.6, 6.5, 0.0, 0.0,         # Scenario 2
               5.0, 1.7, 65.2, 26.8, 0.0, 0.0,        # Scenario 3
               18.0, 4.6, 40.4, 12.6, 6.4, 1.2,       # Scenario 4
               19.6, 4.6, 43.6, 12.6, 4.4, 1.2),      # Scenario 5
  BOINET = c(24.8, 20.0, 96.2, 80.2, 0.0, 0.0,        # Scenario 1
             17.2, 12.1, 60.0, 39.5, 0.0, 0.0,        # Scenario 2
             1.2, 1.7, 65.2, 53.1, 0.0, 0.0,          # Scenario 3
             38.4, 29.0, 73.6, 58.6, 15.6, 9.0,       # Scenario 4
             31.2, 23.4, 86.4, 63.5, 6.4, 3.9)        # Scenario 5
)

# Create data for Main Simulation (Scenarios 6-10 from Table 1 continued)
main_data_part2 <- data.frame(
  Scenario = rep(6:10, each = 6),
  Criteria = rep(criteria_labels, 5),
  CFH_Design = c(92.0, 54.5, 92.0, 54.5, 0.0, 0.0,    # Scenario 6
                 58.8, 35.5, 58.8, 35.5, 0.0, 0.0,    # Scenario 7
                 38.2, 26.7, 87.8, 62.4, 0.0, 0.0,    # Scenario 8
                 52.4, 31.8, 80.2, 61.0, 15.0, 14.9,  # Scenario 9
                 64.4, 37.1, 64.4, 37.1, 11.0, 13.0), # Scenario 10
  Thompson = c(99.6, 72.6, 99.6, 72.6, 0.0, 0.0,      # Scenario 6
               68.0, 29.3, 68.0, 29.3, 0.0, 0.0,      # Scenario 7
               26.6, 12.2, 94.2, 59.1, 0.0, 0.0,      # Scenario 8
               47.2, 30.9, 61.2, 61.2, 38.2, 19.6,    # Scenario 9
               66.2, 43.7, 66.2, 43.7, 26.8, 13.8),   # Scenario 10
  TwoStage = c(80.6, 19.5, 80.6, 19.5, 0.0, 0.0,      # Scenario 6
               34.4, 7.2, 34.4, 7.2, 0.0, 0.0,        # Scenario 7
               20.4, 7.2, 60.6, 20.8, 0.0, 0.0,       # Scenario 8
               32.4, 13.6, 58.8, 33.7, 21.6, 5.5,     # Scenario 9
               43.4, 13.6, 43.4, 13.6, 16.8, 5.5),    # Scenario 10
  BOINET = c(87.2, 63.2, 87.2, 63.2, 0.0, 0.0,        # Scenario 6
             54.4, 36.8, 54.4, 36.8, 0.0, 0.0,        # Scenario 7
             12.6, 10.5, 64.6, 51.7, 0.0, 0.0,        # Scenario 8
             36.4, 45.2, 47.8, 60.5, 45.4, 27.2,      # Scenario 9
             52.2, 45.6, 52.2, 45.6, 11.8, 5.9)       # Scenario 10
)

main_data <- bind_rows(main_data_part1, main_data_part2)

# Create data for Sensitivity Analysis 1 (Reduced number of cohorts to 15; Scenarios 1-5 from Table 2)
sensitivity1_part1 <- data.frame(
  Scenario = rep(1:5, each = 6),
  Criteria = rep(criteria_labels, 5),
  CFH_Design = c(58.8, 30.2, 87.0, 56.1, 0.0, 0.0,    # Scenario 1
                 39.4, 25.1, 72.8, 48.2, 0.0, 0.0,    # Scenario 2
                 23.8, 16.7, 94.2, 74.0, 0.0, 0.0,    # Scenario 3
                 36.6, 23.1, 64.4, 45.6, 16.6, 13.7,  # Scenario 4
                 45.4, 26.1, 81.4, 53.7, 9.8, 10.5),  # Scenario 5
  Thompson = c(77.2, 39.1, 99.4, 72.1, 0.0, 0.0,      # Scenario 1
               14.2, 3.5, 63.0, 23.3, 0.0, 0.0,       # Scenario 2
               2.8, 1.0, 95.6, 61.8, 0.0, 0.0,        # Scenario 3
               41.4, 16.2, 78.0, 43.3, 13.6, 4.0,     # Scenario 4
               41.0, 16.2, 85.0, 46.7, 8.8, 2.5),     # Scenario 5
  TwoStage = c(50.6, 11.5, 84.2, 30.3, 0.0, 0.0,      # Scenario 1
               13.8, 2.6, 31.4, 7.2, 0.0, 0.0,        # Scenario 2
               7.2, 2.6, 66.8, 26.7, 0.0, 0.0,        # Scenario 3
               17.0, 4.6, 37.4, 12.5, 9.8, 2.1,       # Scenario 4
               19.2, 4.6, 43.0, 12.5, 7.2, 2.1),      # Scenario 5
  BOINET = c(23.6, 18.5, 97.6, 77.5, 0.0, 0.0,        # Scenario 1
             13.0, 6.8, 50.8, 26.9, 0.0, 0.0,         # Scenario 2
             1.0, 1.0, 62.0, 45.9, 0.0, 0.0,          # Scenario 3
             35.8, 22.4, 71.2, 50.6, 12.2, 6.2,       # Scenario 4
             24.6, 14.8, 82.4, 54.5, 5.6, 2.7)        # Scenario 5
)

# Create data for Sensitivity Analysis 1 (Scenarios 6-10 from Table 2 continued)
sensitivity1_part2 <- data.frame(
  Scenario = rep(6:10, each = 6),
  Criteria = rep(criteria_labels, 5),
  CFH_Design = c(86.6, 48.0, 86.6, 48.0, 0.0, 0.0,    # Scenario 6
                 57.0, 32.6, 57.0, 32.6, 0.0, 0.0,    # Scenario 7
                 32.4, 22.8, 78.8, 55.2, 0.0, 0.0,    # Scenario 8
                 47.2, 28.5, 77.4, 57.9, 18.0, 16.0,  # Scenario 9
                 58.2, 33.4, 58.2, 33.4, 12.4, 13.4), # Scenario 10
  Thompson = c(97.4, 64.5, 97.4, 64.5, 0.0, 0.0,      # Scenario 6
               49.4, 18.3, 49.4, 18.3, 0.0, 0.0,      # Scenario 7
               23.2, 9.9, 89.0, 48.8, 0.0, 0.0,       # Scenario 8
               47.0, 27.7, 68.8, 59.6, 28.8, 14.0,    # Scenario 9
               69.4, 38.2, 69.4, 38.2, 16.0, 7.7),    # Scenario 10
  TwoStage = c(84.4, 21.2, 84.4, 21.2, 0.0, 0.0,      # Scenario 6
               33.2, 7.0, 33.2, 7.0, 0.0, 0.0,        # Scenario 7
               20.0, 7.0, 61.6, 19.9, 0.0, 0.0,       # Scenario 8
               34.0, 12.9, 62.4, 34.1, 21.0, 5.6,     # Scenario 9
               44.4, 12.9, 44.4, 12.9, 16.4, 5.6),    # Scenario 10
  BOINET = c(86.2, 58.6, 86.2, 58.6, 0.0, 0.0,        # Scenario 6
             48.2, 29.8, 48.2, 29.8, 0.0, 0.0,        # Scenario 7
             13.0, 11.0, 63.8, 46.8, 0.0, 0.0,        # Scenario 8
             38.4, 39.5, 51.4, 58.8, 42.6, 24.5,      # Scenario 9
             49.2, 38.7, 49.2, 38.7, 11.8, 7.3)       # Scenario 10
)

sensitivity1_data <- bind_rows(sensitivity1_part1, sensitivity1_part2)

# Create data for Sensitivity Analysis 2 (Reduced cohort size from 4 to 3; Scenarios 1-5 from Table 3)
sensitivity2_part1 <- data.frame(
  Scenario = rep(1:5, each = 6),
  Criteria = rep(criteria_labels, 5),
  CFH_Design = c(65.6, 35.4, 94.6, 63.1, 0.0, 0.0,    # Scenario 1
                 49.4, 28.4, 82.6, 55.7, 0.0, 0.0,    # Scenario 2
                 26.8, 19.2, 96.4, 78.8, 0.0, 0.0,    # Scenario 3
                 42.8, 26.2, 71.8, 51.0, 19.4, 15.9,  # Scenario 4
                 49.0, 30.3, 75.0, 56.2, 16.4, 13.9), # Scenario 5
  Thompson = c(80.8, 44.4, 98.6, 71.0, 0.0, 0.0,      # Scenario 1
               24.4, 7.4, 73.8, 31.3, 0.0, 0.0,       # Scenario 2
               7.2, 2.8, 96.8, 70.3, 0.0, 0.0,        # Scenario 3
               40.8, 19.9, 71.4, 46.7, 17.0, 6.0,     # Scenario 4
               42.2, 20.5, 73.0, 50.2, 18.4, 5.9),    # Scenario 5
  TwoStage = c(55.2, 13.7, 87.2, 31.3, 0.0, 0.0,      # Scenario 1
               17.8, 2.8, 40.4, 9.1, 0.0, 0.0,        # Scenario 2
               10.4, 2.8, 74.4, 32.6, 0.0, 0.0,       # Scenario 3
               16.6, 6.2, 39.0, 15.8, 17.4, 2.3,      # Scenario 4
               21.0, 6.2, 46.4, 15.8, 11.6, 2.3),     # Scenario 5
  BOINET = c(27.6, 21.4, 94.2, 79.4, 0.0, 0.0,        # Scenario 1
             18.4, 8.2, 57.0, 31.2, 0.0, 0.0,         # Scenario 2
             2.0, 1.0, 67.2, 50.6, 0.0, 0.0,          # Scenario 3
             36.4, 24.3, 73.2, 55.3, 14.6, 6.0,       # Scenario 4
             29.8, 17.8, 85.0, 58.3, 5.8, 2.7)        # Scenario 5
)

# Create data for Sensitivity Analysis 2 (Scenarios 6-10 from Table 3 continued)
sensitivity2_part2 <- data.frame(
  Scenario = rep(6:10, each = 6),
  Criteria = rep(criteria_labels, 5),
  CFH_Design = c(91.2, 51.3, 91.2, 51.3, 0.0, 0.0,    # Scenario 6
                 66.0, 37.1, 66.0, 37.1, 0.0, 0.0,    # Scenario 7
                 35.8, 26.8, 84.4, 61.3, 0.0, 0.0,    # Scenario 8
                 48.0, 30.2, 70.8, 57.2, 23.8, 20.2,  # Scenario 9
                 61.8, 36.1, 61.8, 36.1, 17.0, 16.7), # Scenario 10
  Thompson = c(99.0, 71.7, 99.0, 71.7, 0.0, 0.0,      # Scenario 6
               64.4, 27.7, 64.4, 27.7, 0.0, 0.0,      # Scenario 7
               31.4, 14.3, 92.8, 61.0, 0.0, 0.0,      # Scenario 8
               46.0, 30.6, 64.0, 60.6, 34.6, 18.6,    # Scenario 9
               68.8, 44.8, 68.8, 44.8, 21.6, 12.1),   # Scenario 10
  TwoStage = c(87.4, 20.2, 87.4, 20.2, 0.0, 0.0,      # Scenario 6
               48.0, 9.9, 48.0, 9.9, 0.0, 0.0,        # Scenario 7
               28.6, 9.9, 71.8, 25.4, 0.0, 0.0,       # Scenario 8
               32.0, 15.5, 53.6, 36.3, 31.6, 8.2,     # Scenario 9
               46.2, 15.5, 46.2, 15.5, 23.8, 8.2),    # Scenario 10
  BOINET = c(78.6, 58.7, 78.6, 58.7, 0.0, 0.0,        # Scenario 6
             57.0, 33.6, 57.0, 33.6, 0.0, 0.0,        # Scenario 7
             15.4, 10.8, 68.2, 48.6, 0.0, 0.0,        # Scenario 8
             37.6, 43.6, 52.8, 62.0, 42.8, 23.7,      # Scenario 9
             58.2, 47.3, 58.2, 47.3, 13.2, 5.4)       # Scenario 10
)

sensitivity2_data <- bind_rows(sensitivity2_part1, sensitivity2_part2)

# Function to create sensitivity analysis plots with scenario numbers in ALL facets
create_sensitivity_plot <- function(data, title) {
  
  # Reshape data for ggplot - INCLUDING BOINET
  data_long <- data %>%
    pivot_longer(cols = c(CFH_Design, Thompson, TwoStage, BOINET),
                 names_to = "Design",
                 values_to = "Percentage")
  
  # Clean up design names
  data_long$Design <- gsub("_", " ", data_long$Design)
  data_long$Design <- ifelse(data_long$Design == "TwoStage", "Two-Stage", data_long$Design)
  data_long$Design <- factor(data_long$Design, 
                             levels = c("CFH Design", "Thompson", "Two-Stage", "BOINET"))
  # Create grouping variables
  data_long <- data_long %>%
    mutate(
      Dose_Type = case_when(
        grepl("most efficacious", Criteria) ~ "Most Efficacious",
        grepl("non-inferior", Criteria) ~ "Non-inferior",
        grepl("overdose", Criteria) ~ "Overdose"
      ),
      Metric_Type = ifelse(grepl("Selection", Criteria), "Selection", "Allocation"),
      Facet_Group = paste(Dose_Type, "-", Metric_Type)
    )
  
  # Create factor for facet order
  facet_order <- c(
    "Most Efficacious - Selection",
    "Most Efficacious - Allocation",
    "Non-inferior - Selection", 
    "Non-inferior - Allocation",
    "Overdose - Selection",
    "Overdose - Allocation"
  )
  
  data_long$Facet_Group <- factor(data_long$Facet_Group, levels = facet_order)
  
  # Create facet labels
  facet_labels <- c(
    "Most Efficacious - Selection" = "Selection of Most Efficacious Dose",
    "Most Efficacious - Allocation" = "Patient Allocation to Most Efficacious Dose",
    "Non-inferior - Selection" = "Selection of Non-inferior Doses",
    "Non-inferior - Allocation" = "Patient Allocation to Non-inferior Doses",
    "Overdose - Selection" = "Selection of Overdose",
    "Overdose - Allocation" = "Patient Allocation to Overdose"
  )
  
  # Create the plot
  p <- ggplot(data_long, aes(x = factor(Scenario), y = Percentage, fill = Design)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    facet_wrap(~ Facet_Group, ncol = 2, scales = "free_y",
               labeller = labeller(Facet_Group = facet_labels)) +
    scale_fill_manual(values = c("CFH Design" = "#1f77b4", 
                                 "Thompson" = "#ff7f0e", 
                                 "Two-Stage" = "#2ca02c",
                                 "BOINET" = "#d62728")) +
    labs(title = title,
         x = "Scenario",
         y = "Percentage (%)",
         fill = "Design") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
      strip.placement = "outside",
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = NA),
      panel.spacing = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
  
  return(p)
}

# Create all plots
plot_main <- create_sensitivity_plot(main_data, title = "Main Simulation Study")
plot_sensitivity1 <- create_sensitivity_plot(sensitivity1_data, title = "Sensitivity Analysis: Reduced Number of Cohorts")
plot_sensitivity2 <- create_sensitivity_plot(sensitivity2_data, title = "Sensitivity Analysis: Reduced Cohort Size")

# Display the plots
print(plot_main)
print(plot_sensitivity1)
print(plot_sensitivity2)

# historical data
historical_data_part1 <- data.frame(
  Scenario = rep(1:5, each = 6),
  Criteria = rep(c("Selection of most efficacious dose %",
                   "Patient allocation to most efficacious dose %",
                   "Selection of non-inferior doses %",
                   "Patient allocation to non-inferior doses %",
                   "Selection of overdose %",
                   "Patient allocation to overdose %"), 5),
  No_Historical_Data = c(70.6, 39.5, 94.6, 66.0, 0.0, 0.0,      # Scenario 1
                         44.2, 27.3, 80.0, 55.9, 0.0, 0.0,      # Scenario 2
                         27.2, 19.7, 94.4, 77.6, 0.0, 0.0,      # Scenario 3
                         45.4, 27.7, 70.6, 51.9, 16.2, 14.3,    # Scenario 4
                         48.8, 30.6, 78.8, 57.8, 11.4, 11.1),   # Scenario 5
  Commensurate = c(89.2, 57.4, 99.6, 84.0, 0.0, 0.0,            # Scenario 1
                   72.6, 50.1, 88.4, 82.0, 0.0, 0.0,            # Scenario 2
                   42.6, 32.4, 98.0, 91.1, 0.0, 0.0,            # Scenario 3
                   59.0, 37.3, 75.6, 57.2, 19.6, 29.3,          # Scenario 4
                   76.4, 52.7, 94.8, 85.3, 0.0, 0.0),           # Scenario 5
  Conflicting = c(88.2, 55.5, 99.0, 75.4, 0.0, 0.0,             # Scenario 1
                  69.8, 48.0, 87.8, 79.3, 0.0, 0.0,               # Scenario 2
                  39.4, 32.4, 98.6, 89.6, 0.0, 0.0,             # Scenario 3
                  55.2, 37.6, 71.8, 55.3, 22.4, 32.2,           # Scenario 4
                  59.2, 41.5, 79.4, 73.6, 14.6, 12.7)           # Scenario 5
)

# Create data for historical-data comparison (Scenarios 6-10 from Historical2)
historical_data_part2 <- data.frame(
  Scenario = rep(6:10, each = 6),
  Criteria = rep(c("Selection of most efficacious dose %",
                   "Patient allocation to most efficacious dose %",
                   "Selection of non-inferior doses %",
                   "Patient allocation to non-inferior doses %",
                   "Selection of overdose %",
                   "Patient allocation to overdose %"), 5),
  No_Historical_Data = c(92.0, 54.5, 92.0, 54.5, 0.0, 0.0,      # Scenario 6
                         58.8, 35.5, 58.8, 35.5, 0.0, 0.0,      # Scenario 7
                         38.2, 26.7, 87.8, 62.4, 0.0, 0.0,      # Scenario 8
                         52.4, 31.8, 80.2, 61.0, 15.0, 14.9,    # Scenario 9
                         64.4, 37.1, 64.4, 37.1, 11.0, 13.0),   # Scenario 10
  Commensurate = c(99.8, 74.3, 99.8, 74.3, 0.0, 0.0,            # Scenario 6
                   79.4, 47.9, 79.4, 47.9, 0.0, 0.0,            # Scenario 7
                   55.4, 44.0, 96.4, 86.3, 0.0, 0.0,            # Scenario 8
                   66.6, 40.0, 83.8, 67.3, 14.6, 23.3,          # Scenario 9
                   82.0, 54.0, 82.0, 54.0, 14.6, 23.0),         # Scenario 10
  Conflicting = c(100.0, 70.3, 100.0, 70.3, 0.0, 0.0,           # Scenario 6
                  86.2, 52.7, 86.2, 52.7, 0.0, 0.0,             # Scenario 7
                  49.6, 37.3, 96.6, 84.9, 0.0, 0.0,             # Scenario 8
                  66.8, 41.4, 77.8, 57.9, 20.2, 34.6,           # Scenario 9
                  71.4, 54.5, 71.4, 54.5, 18.8, 20.1)           # Scenario 10
)

# Combine historical-data comparison data
historical_data <- bind_rows(historical_data_part1, historical_data_part2)

# Function to create the historical-data comparison plot
create_historical_plot <- function(data, title) {
  
  data_long <- data %>%
    pivot_longer(cols = c(No_Historical_Data, Commensurate, Conflicting),
                 names_to = "Historical_Setting",
                 values_to = "Percentage") %>%
    mutate(
      Historical_Setting = recode(Historical_Setting,
                                  "No_Historical_Data" = "No Historical Data",
                                  "Commensurate" = "Commensurate",
                                  "Conflicting" = "Conflicting"),
      Historical_Setting = factor(Historical_Setting,
                                  levels = c("No Historical Data",
                                             "Commensurate",
                                             "Conflicting")),
      Dose_Type = case_when(
        grepl("most efficacious", Criteria) ~ "Most Efficacious",
        grepl("non-inferior", Criteria) ~ "Non-inferior",
        grepl("overdose", Criteria) ~ "Overdose"
      ),
      Metric_Type = ifelse(grepl("Selection", Criteria), "Selection", "Allocation"),
      Facet_Group = paste(Dose_Type, "-", Metric_Type)
    )
  
  facet_order <- c(
    "Most Efficacious - Selection",
    "Most Efficacious - Allocation",
    "Non-inferior - Selection",
    "Non-inferior - Allocation",
    "Overdose - Selection",
    "Overdose - Allocation"
  )
  
  facet_labels <- c(
    "Most Efficacious - Selection" = "Selection of Most Efficacious Dose",
    "Most Efficacious - Allocation" = "Patient Allocation to Most Efficacious Dose",
    "Non-inferior - Selection" = "Selection of Non-inferior Doses",
    "Non-inferior - Allocation" = "Patient Allocation to Non-inferior Doses",
    "Overdose - Selection" = "Selection of Overdose",
    "Overdose - Allocation" = "Patient Allocation to Overdose"
  )
  
  data_long$Facet_Group <- factor(data_long$Facet_Group, levels = facet_order)
  
  ggplot(data_long, aes(x = factor(Scenario), y = Percentage, fill = Historical_Setting)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.75) +
    facet_wrap(~ Facet_Group, ncol = 2, scales = "free_y",
               labeller = labeller(Facet_Group = facet_labels)) +
    scale_fill_manual(values = c("No Historical Data" = "#1f77b4",
                                 "Commensurate" = "#ff7f0e",
                                 "Conflicting" = "#d62728")) +
    labs(title = title,
         x = "Scenario",
         y = "Percentage (%)",
         fill = "Historical Data Setting") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
      strip.placement = "outside",
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = NA),
      panel.spacing = grid::unit(1, "lines"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = as.character(1:10))
}

# Create and display the plot
plot_historical <- create_historical_plot(
  historical_data,
  title = "Historical Data Simulation Study"
)

print(plot_historical)