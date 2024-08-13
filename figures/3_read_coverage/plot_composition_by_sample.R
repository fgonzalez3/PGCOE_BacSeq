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

heattab <- heattab %>%
  rename(Seq_ID = sample_name)

heattab <- merge(heattab,seqtab) %>% filter(taxon_name %in% names(cols))
heattab$taxon_name <- factor(heattab$taxon_name,levels=rev(names(cols)),ordered=T)

# filter out samples where we don't have a match 
reps <- c("Yale-SP00051", "Yale-SP00052", "Yale-SP00053", "Yale-SP00054",
          "Yale-SP00055", "Yale-SP00056", "Yale-SP00057", "Yale-SP00058", "Yale-SP00059",
          "Yale-CS00171", "Yale-CS00172", "Yale-CS00173", "Yale-CS00174", "Yale-CS00175",
          "Yale-CS00176", "Yale-CS00162", "Yale-CS00163", "Yale-CS00164")

# rename 
name_mapping <- c("Yale-SP00051" = "A889-NoAMP", 
                  "Yale-SP00052" = "B042-NoAMP",
                  "Yale-SP00053" = "C677-NoAmp", 
                  "Yale-SP00054" = "A889-NoAmp",
                  "Yale-SP00055" = "B042-NoAmp", 
                  "Yale-SP00056" = "C677-NoAmp", 
                  "Yale-SP00057" = "W1527-NoAmp", 
                  "Yale-SP00058" = "W3317-NoAmp", 
                  "Yale-SP00059" = "W4034-NoAmp",
                  "Yale-CS00171" = "A889-Amp", 
                  "Yale-CS00172" = "B042-Amp", 
                  "Yale-CS00173" = "C677-Amp", 
                  "Yale-CS00174" = "W1527-Amp", 
                  "Yale-CS00175" = "W3317-Amp",
                  "Yale-CS00176" = "W4034-Amp", 
                  "Yale-CS00162" = "A889-Amp", 
                  "Yale-CS00163" = "B042-Amp", 
                  "Yale-CS00164" = "C677-Amp")

renamed_reps <- name_mapping[reps]

filtered_heattab <- heattab %>% filter(ID %in% reps)
filtered_heattab$ID <- name_mapping[filtered_heattab$ID]
filtered_heattab$ID <- factor(filtered_heattab$ID, levels = unique(name_mapping[reps]))

sppplot <- ggplot(filtered_heattab, aes(x = ID, y = NT_rpm, fill = taxon_name)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        axis.title.x = element_blank()) + 
  facet_grid(. ~ Sample_Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = ggtext::element_markdown())
sppplot

ggsave("bacseq_spneumo_czid_props.png",width=300,height=140,units="mm")












# ------------------------------------------------------------------------------

covtab <- read.table("coverage_stats.csv",header=1,sep=",") %>% 
  rename("strain"="ID") %>% 
  mutate(ID=paste("Yale-",sample,sep="")) %>% 
  filter(readsample=="5m")
covtab <- merge(covtab,seqtab)

covplot <- ggplot(covtab,aes(x=ID,y=coverage,fill=meandepth)) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position = "bottom",
        axis.title.x = element_blank()) + 
  facet_grid(. ~ Sample.type,scales="free_x",space="free_x")
covplot
ggsave("bacseq_spneumo_czid_cov.png",width=300,height=140,units="mm")


sppplot + covplot + plot_layout(heights=c(3,2))

