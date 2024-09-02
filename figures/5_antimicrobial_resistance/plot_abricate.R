rm(list = ls())

library("tidyverse")
library(RColorBrewer)
library(patchwork)
library(ggtext)

# read in csv
amr <- read.csv("card_summary.csv")

# pick reps to visualize
reps <- c("Yale-CS00030", "Yale-CS00032", "Yale-CS00034", "Yale-CS00036", "Yale-CS00038",
          "Yale-CS00040", "Yale-CS00174", "Yale-CS00175", "Yale-CS00176")

# create a grid of samples and drug classes
complete_grid <- expand.grid(sample_name = unique(amr$sample_name), 
                             RESISTANCE = unique(amr$RESISTANCE))

# join grid 
amr_complete <- complete_grid %>%
  left_join(amr, by = c("sample_name", "RESISTANCE"))

# make new colname for resistance status
amr_complete <- amr_complete %>%
  mutate(resistance_status = ifelse(is.na(RESISTANCE), "Susceptible", "Resistant"))

unique(amr_complete$resistance_status)


# plot
amr_plot <- amr_complete %>%
  filter(sample_name %in% reps) %>%
  ggplot(aes(sample_name, RESISTANCE, fill = resistance_status)) + 
  geom_tile() + 
  scale_fill_manual(values = c("Resistant" = "purple", "Susceptible" = "lightblue"), 
                    name = "Resistance Status") + 
 # scale_y_discrete(labels = new_drug_class_names) + 
  labs(x = "Sample Name", y = "Drug Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

amr_plot

ggsave("bacseq_spneumo_czid_amr.png",width=350,height=300,units="mm", bg = "white")