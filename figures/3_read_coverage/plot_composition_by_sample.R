rm(list = ls())

library("tidyverse")
library(RColorBrewer)
library(patchwork)
library(ggtext)

# read counts by sample type & sequencing method -------------------------------

heattab <- read.table("data/BacSeq_amplicon_mNGS_comparison_CZID_heatmap_12AUG2024.csv",sep=",",header=1)
heattab$ID <- gsub("_S.+","",heattab$sample_name,perl=T)

#seqtab <- read.table("GLab_seq_master - Seq opt.txt",sep="\t",header=1)
seqtab <- read.csv("metadata/PGCOE_DataSheet - S. pneumoniae.csv")

brewer.pal(5,name = "Blues")
brewer.pal(5,name = "Greens")


cols = c("Streptococcus pneumoniae"="#006D2C","Streptococcus mitis"="#BDD7E7",
         "Streptococcus oralis"="#6BAED6","Streptococcus sp. oral taxon 061"="#3182BD","Streptococcus sp. oral taxon 431"="#08519C")

# create a vector of grey scale colors
grey_scale <- colorRampPalette(c("grey80", "grey20"))(length(unique(heattab$taxon_name)) - length(cols))

# identify taxa that are not in the cols vector
non_streptococcus_taxa <- setdiff(unique(heattab$taxon_name), names(cols))

# assign grey scale colors to these taxa
grey_cols <- setNames(grey_scale, non_streptococcus_taxa)

# combine the grey scale colors with the existing cols vector
all_cols <- c(cols, grey_cols)

heattab <- heattab %>%
  rename(Seq_ID = sample_name)

heattab <- merge(heattab, seqtab) %>% filter(taxon_name %in% names(all_cols))
heattab$taxon_name <- factor(heattab$taxon_name, levels=rev(names(all_cols)), ordered=T)

# filter out samples where we don't have a match 
reps <- c("Yale-SP00051", "Yale-SP00052", "Yale-SP00053", "Yale-SP00054",
          "Yale-SP00055", "Yale-SP00056", "Yale-SP00057", "Yale-SP00058", "Yale-SP00059",
          "Yale-CS00171", "Yale-CS00172", "Yale-CS00173", "Yale-CS00174", "Yale-CS00175",
          "Yale-CS00176", "Yale-CS00162", "Yale-CS00163", "Yale-CS00164")


# rename IDs
name_mapping <- c("Yale-SP00051" = "A889-Unamplified", 
                  "Yale-SP00052" = "B042-Unamplified",
                  "Yale-SP00053" = "C677-Unamplified", 
                  "Yale-SP00054" = "A889-Unamplified",
                  "Yale-SP00055" = "B042-Unamplified", 
                  "Yale-SP00056" = "C677-Unamplified", 
                  "Yale-SP00057" = "W1527-Unamplified", 
                  "Yale-SP00058" = "W3317-Unamplified", 
                  "Yale-SP00059" = "W4034-Unamplified",
                  "Yale-CS00171" = "A889-Amplified", 
                  "Yale-CS00172" = "B042-Amplified", 
                  "Yale-CS00173" = "C677-Amplified", 
                  "Yale-CS00174" = "W1527-Amplified", 
                  "Yale-CS00175" = "W3317-Amplified",
                  "Yale-CS00176" = "W4034-Amplified", 
                  "Yale-CS00162" = "A889-Amplified", 
                  "Yale-CS00163" = "B042-Amplified", 
                  "Yale-CS00164" = "C677-Amplified")


# rename reps
renamed_reps <- name_mapping[reps]

# filter and rename IDs in the data frame
filtered_heattab <- heattab %>% 
  filter(ID %in% reps) %>%
  mutate(ID = name_mapping[ID])

# calculate total RPM for each ID, including all IDs
id_order <- filtered_heattab %>%
  group_by(ID) %>%
  summarise(total_rpm = sum(NT_rpm), max_rpm = max(NT_rpm)) %>%
  arrange(desc(total_rpm)) %>%
  pull(ID)

# create a complete list of IDs including both -Amp and -NoAmp types
all_ids <- unique(c(id_order, gsub("-Amp", "-NoAmp", id_order), gsub("-NoAmp", "-Amp", id_order)))

# set the factor levels of ID based on the complete list
filtered_heattab$ID <- factor(filtered_heattab$ID, levels = all_ids)

# create a new column for labels based on the Method column
#filtered_heattab <- filtered_heattab %>%
  #mutate(Label = case_when(
   # Method == "Ampllicon" ~ "Amplified",
   # Method == "mNGS" ~ "Unamplified"
 # ))

# create a summary data frame for labels
#label_data <- filtered_heattab %>%
 # group_by(ID) %>%
  #summarise(Label = first(Label), max_NT_rpm = max(NT_rpm))

# plot reads per million by sample ---------------------------------------------

# plot
sppplot <- ggplot(filtered_heattab, aes(x = ID, y = NT_rpm, fill = taxon_name)) + 
  geom_bar(stat = "identity") + 
  #geom_text(data = label_data, aes(x = ID, y = max_NT_rpm + 1, label = Label), vjust = -0.5, size = 3, inherit.aes = FALSE) +  # Add labels above bars
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        axis.title.x = element_blank()) + 
  facet_grid(~Sample_Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = all_cols, 
                    breaks = c("Streptococcus pneumoniae", "Streptococcus mitis",
                               "Streptococcus oralis", "Streptococcus sp. oral taxon 061",
                               "Streptococcus sp. oral taxon 431")) +  
  theme(axis.text.x = ggtext::element_markdown())

sppplot
ggsave("bacseq_spneumo_czid_props.png",width=350,height=300,units="mm")


# plot percentage of reads by sample -------------------------------------------

# rename reps
renamed_reps <- name_mapping[reps]

# filter for the representative sequences we want to visualize
filtered_heattab <- heattab %>% 
  filter(ID %in% reps) %>%
  mutate(ID = name_mapping[ID])

# calculate total amount of reads per sequencing run
total_reads_per_id <- filtered_heattab %>%
  group_by(ID) %>%
  summarise(total_reads = sum(NT_r))

# calculate % of reads that map to each taxon 
filtered_heattab <- filtered_heattab %>%
  left_join(total_reads_per_id, by = "ID") %>%
  mutate(percentage_reads = (NT_r / total_reads) * 100)

# arrange ID order 
id_order <- filtered_heattab %>%
  group_by(ID) %>%
  summarise(max_percentage = max(percentage_reads)) %>%
  arrange(desc(max_percentage)) %>%
  pull(ID)

# make list of ID order by seq method
all_ids <- unique(c(id_order, gsub("-Amp", "-NoAmp", id_order), gsub("-NoAmp", "-Amp", id_order)))

# factor levels of IDs
filtered_heattab$ID <- factor(filtered_heattab$ID, levels = all_ids)

# plot
sppplot <- ggplot(filtered_heattab, aes(x = ID, y = percentage_reads, fill = taxon_name)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        axis.title.x = element_blank()) + 
  facet_grid(~Sample_Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = all_cols, 
                    breaks = c("Streptococcus pneumoniae", "Streptococcus mitis",
                               "Streptococcus oralis", "Streptococcus sp. oral taxon 061",
                               "Streptococcus sp. oral taxon 431")) +  
  theme(axis.text.x = ggtext::element_markdown()) + 
  ylab("% Reads")

sppplot

ggsave("bacseq_spneumo_czid_pc.png",width=350,height=300,units="mm")
