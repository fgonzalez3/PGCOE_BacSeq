library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(plyr)
library(wesanderson)

homewd= "/Users/flg9/Desktop/Developer/grubaugh_lab/snakemake_workflows/read_alignment/R_script/"
setwd(paste0(homewd))

# load in mNGS data
mNGS_seq <- read.csv("mNGS_cov_stats_12JUN2024.csv", 
                     header = T, stringsAsFactors = F, sep = "\t")

# load in PGCOE metadata 
pgcoe_metadat <- read.csv("PGCOE_DataSheet_12JUN2024.csv", 
                       header = T, stringsAsFactors = F)

# modify variable names in pgcoe dataset to match those of coverage data
pgcoe_metadat$sample <- str_replace(pgcoe_metadat$sample, ".*\\-", "")

# create new df to merge original PGCOE metadata to new coverage data 
mNGS_seq_merge <- merge(mNGS_seq, pgcoe_metadat, by = "sample", sort = T, no.dups = T)

# sort sample ID order for plot 
unique(mNGS_seq_merge$Original_ID) %>% sort

mNGS_seq_merge$ID <- factor(mNGS_seq_merge$Original_ID,      # Reordering group factor levels
                       levels = c("A889-NP", "A889-S", "B042-NP",  "B042-S",
                                  "C677-NP", "C677-S", "NTC3", "W1527-CI", 
                                  "W3317-CI", "W4034-CI"))

# plot coverage for samples from mNGS run
mNGS_seq_merge %>% ggplot() +
  geom_point(aes(Original_ID, coverage, col = subsample)) + 
  ylab("Coverage (%)") + xlab("Sample ID") 

## save mNGS merge to a csv
write.csv(mNGS_seq_merge, file = "mNGS_seq_13JUN2024.csv")

################## now do the same with untrimmed amplicon data #######################

# load in amplicon data
amp_seq <- read.csv("amp_cov_stats_12JUN2024.csv", 
                     header = T, stringsAsFactors = F, sep = "\t")

# create new df to merge original PGCOE metadata to new coverage data 
amp_seq_merge <- merge(amp_seq, pgcoe_metadat, by = "sample", sort = T, no.dups = T)

# look like two of the NEC samples aren't on the pgcoe metadata so those don't merge
# that's ok since we have plenty of other NECs 
unique(amp_seq$sample)
unique(amp_seq_merge$sample) 
unique(amp_seq_merge$Original_ID)

# plot coverage for samples from ampseq run
amp_seq_merge %>%
  filter(Primer_Conc == "100uM" & subsample == "1.0" & (Sample_Dilution == "1e-03" | !is.na(Sample_Dilution | !is.na(Ct1)))) %>%
  ggplot() +
  geom_point(aes(Ct1, coverage, col = Sample_Type), size = 4) + 
  ylab("Coverage(%)") + xlab ("Ct") + 
  geom_hline(yintercept = 80, linetype = "dashed", color = "salmon") +
  geom_text(x = max(amp_seq_merge$Ct1), y = 80, label = "80%", vjust = 2, hjust = 2, color = "black", size = 5) + 
  scale_y_continuous(limits = c(-15, 105), breaks = c(0, 25, 50, 75, 100)) + 
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) + 
  scale_color_manual(values = wes_palette("Darjeeling2", 8, type = "continuous"))

# save df to a csv 
# first save raw file
write.csv(amp_seq_merge, file = "amplicon_seq_13JUN2024.csv")

# now save filtered file for Chaney to access relevant data easier 
amp_seq_merged_filtered <- amp_seq_merge %>%
  filter(Primer_Conc == "100uM" & subsample == "1.0" & (Sample_Dilution == "1e-03" | is.na(Sample_Dilution | is.na(Ct1))))
  
################### now compare this with the trimmed amplicon data ####################

# load in amplicon data
ivar_amp_seq <- read.csv("ivar_amp_cov_stats.csv", 
                    header = T, stringsAsFactors = F, sep = "\t")


# create new df to merge original PGCOE metadata to new coverage data 
ivar_amp_seq_merge <- merge(ivar_amp_seq, pgcoe_metadat, by = "sample", sort = T, no.dups = T)

## something is going wrong with the trimming so we'll only plot the untrimmed amplicon sequencing 



