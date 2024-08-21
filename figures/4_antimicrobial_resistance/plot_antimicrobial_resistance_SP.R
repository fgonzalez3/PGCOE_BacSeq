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
  "monobactam; carbapenem; cephalosporin; cephamycin; penam; phenicol antibiotic; penem" = "Beta-lactams and Phenicols",
  "mupirocin-like antibiotic" = "Mupirocin",
  "glycopeptide antibiotic" = "Glycopeptides", 
  "fluoroquinolone antibiotic" = "Fluoroquinolone", 
  "cephalosporin; penam; peptide antibiotic" = "Cephalosporins and Penams",
  "carbapenem" = "Carbapenems"
)

# create a grid of samples and drug classes
complete_grid <- expand.grid(sample_name = unique(amr$sample_name), 
                             drug_class = unique(amr$drug_class))

# join grid 
amr_complete <- complete_grid %>%
  left_join(amr, by = c("sample_name", "drug_class"))

# make new colname for resistance status
amr_complete <- amr_complete %>%
  mutate(resistance_status = ifelse(is.na(resistance_mechanism), "Susceptible", "Resistant"))

# filter out repetitive resistant markers 
amr_filtered <- amr_complete %>%
  filter(drug_class != "macrolide antibiotic") %>%
  filter(drug_class != "macrolide antibiotic; fluoroquinolone antibiotic; monobactam; aminoglycoside antibiotic; carbapenem; cephalosporin; cephamycin; penam; tetracycline antibiotic; phenicol antibiotic; penem; disinfecting agents and antiseptics")

# plot
amr_plot <- amr_filtered %>%
  filter(sample_name %in% reps) %>%
  ggplot(aes(sample_name, drug_class, fill = resistance_status)) + 
  geom_tile() + 
  scale_fill_manual(values = c("Resistant" = "purple", "Susceptible" = "lightblue"), 
                    name = "Resistance Status") + 
  scale_y_discrete(labels = new_drug_class_names) + 
  labs(x = "Sample Name", y = "Drug Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

amr_plot

ggsave("bacseq_spneumo_czid_amr.png",width=350,height=300,units="mm", bg = "white")


############# old code ###########

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






