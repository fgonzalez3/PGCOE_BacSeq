rm(list = ls())

library("tidyverse")
library(RColorBrewer)
library(patchwork)
library(ggtext)

# read in csv
amr <- read.csv("combined_amr_results.csv") 

# pick reps to visualize
reps <- c("Yale-CS00030", "Yale-CS00047", "Yale-CS00048", "Yale-CS00049",
          "Yale-CS00049", "Yale-CS00050", "Yale-CS00051")

# rename drug classes 
new_drug_class_names <- c(
  "tetracycline antibiotic" = "Tetracycline", 
  "phosphonic acid antibiotic" = "Fosfomycin", 
  "macrolide antibiotic; lincosamide antibiotic" = "Macrolides and Lincsamides",
  "macrolide antibiotic; fluoroquinolone antibiotic; aminoglycoside antibiotic; cephalosporin" = "Aminoglycosides",
  "glycopeptide antibiotic" = "Glycopeptides", 
  "fluoroquinolone antibiotic" = "Fluoroquinolone", 
  "cephalosporin; penam; peptide antibiotic" = "Cephalosporins and Penams",
  "carbapenem" = "Carbapenems"
)

# plot
amr %>%
  filter(sample_name %in% reps) %>%
  drop_na(resistance_mechanism) %>%
  ggplot(aes(sample_name, drug_class, fill = resistance_mechanism)) + 
  geom_tile() + 
  scale_fill_brewer(palette = "Set3", 
                    name = "Resistance Mechanism",
                    labels = c("Antibiotic Efflux", "Antibiotic Inactivation", "Antibiotic Target Alteration")) + 
  scale_y_discrete(labels = new_drug_class_names) + 
  labs(x = "Sample Name", y = "Drug Class") +
  theme_minimal()



