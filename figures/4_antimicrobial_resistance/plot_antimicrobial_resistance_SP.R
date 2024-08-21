rm(list = ls())

library("tidyverse")
library(RColorBrewer)
library(patchwork)
library(ggtext)

# read in csv
amr <- read.csv("combined_amr_results_diffsamples.csv") 

# pick reps to visualize
reps <- c("Yale-CS00030", "Yale-CS00032", "Yale-CS00034", "Yale-CS00036", "Yale-CS00038",
          "Yale-CS00040", "Yale-CS00174", "Yale-CS00175", "Yale-CS00176", "Yale-SP00057",
          "Yale-SP00058", "Yale-SP00059")

# rename drug classes 
new_drug_class_names <- c(
  "tetracycline antibiotic" = "Tetracycline", 
  "phosphonic acid antibiotic" = "Fosfomycin", 
  "macrolide antibiotic; lincosamide antibiotic" = "Macrolides and Lincosamides",
  "macrolide antibiotic; fluoroquinolone antibiotic; aminoglycoside antibiotic; cephalosporin" = "Aminoglycosides",
  "glycopeptide antibiotic" = "Glycopeptides", 
  "macrolide antibiotic" = "Macrolides",
  "fluoroquinolone antibiotic" = "Fluoroquinolone", 
  "cephalosporin; penam; peptide antibiotic" = "Cephalosporins and Penams",
  "carbapenem" = "Carbapenems"
)

# plot
amr_plot <- amr %>%
  slice(-23) %>%
  filter(sample_name %in% reps) %>%
  drop_na(resistance_mechanism) %>%
  ggplot(aes(sample_name, drug_class, fill = resistance_mechanism)) + 
  geom_tile() + 
  scale_fill_brewer(palette = "Set3", 
                    name = "Resistance Mechanism",
                    labels = c("Antibiotic Efflux", "Antibiotic Inactivation", 
                               "Antibiotic Target Alteration", "Antibiotic Target Protection")) + 
  scale_y_discrete(labels = new_drug_class_names) + 
  labs(x = "Sample Name", y = "Drug Class") +
  theme_minimal()

amr_plot

ggsave("bacseq_spneumo_czid_amr.png",width=350,height=300,units="mm", bg = "white")




