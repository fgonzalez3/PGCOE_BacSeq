rm(list = ls())

library(tidyverse)
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)
library(readr)
library(wesanderson)

# file paths
SP_amp_cov  <- "SP_amplicon_seq.csv"
SP_mNGS_cov <- "SP_mNGS_seq.csv"
TB_amp_cov  <- "TB_combined_coverage.tsv"

SP_metadat  <- "PGCOE_DataSheet.csv"
TB_metadat  <- "PGCOE_Datasheet - M.tuberculosis_Seq.csv"

outpng      <- "Fig3_combined_coverage_plot.png"

t <- read.csv(SP_metadat)


t1 <- read_table(SP_amp_cov)
t2 <- read_table(SP_mNGS_cov)

t1 <- t1 %>%
  rename_with(~ str_remove _all(., "^\"|\"$")) %>%
  mutate(
    Seq_ID = str_remove_all(Seq_ID, "^\"|\"$"),
    meanmapq = as.numeric(str_remove_all(meanmapq, "^\"|\"$")),
    subsample = as.character(subsample)  # Ensure subsample is character
  )

t2 <- t2 %>%
  rename_with(~ str_remove_all(., "^\"|\"$")) %>%
  mutate(
    Seq_ID = str_remove_all(Seq_ID, "^\"|\"$"),
    meanmapq = as.numeric(str_remove_all(meanmapq, "^\"|\"$")),
    subsample = as.character(subsample)  # Ensure subsample is character
  )


data_combined_spn <- bind_rows(t,t1, t2) %>%
  mutate(sample = str_c("Yale-", Seq_ID)) %>%
  mutate(Seq_ID = str_remove(Seq_ID, "^Yale-")) %>%
  mutate(pathogen = "SPn")


# load and clean data ------------------------------------------------------------

cleanup_data <- function(SP_amp_cov, SP_mNGS_cov, TB_amp_cov, SP_metadat, TB_metadat) {
  
  # remove leading and trailing quotes
  clean_column_names <- function(df) {
    colnames(df) <- gsub('^"|"$', '', colnames(df))  
    return(df)
  }
  
  # load data 
  coverageinfo_tb <- read_delim(TB_amp_cov, delim = "\t")
  datasheet_tb <- read_delim(TB_metadat, delim = ",") 
  coverageinfo_spn_amp <- read_table(SP_amp_cov)
  coverageinfo_spn_mngs <- read_table(SP_mNGS_cov)
  datasheet_spn <- read_delim(SP_metadat, delim = ",")
  
  # clean colnames
  coverageinfo_tb <- clean_column_names(coverageinfo_tb)
  datasheet_tb <- clean_column_names(datasheet_tb)
  coverageinfo_spn_amp <- clean_column_names(coverageinfo_spn_amp)
  coverageinfo_spn_mngs <- clean_column_names(coverageinfo_spn_mngs)
  datasheet_spn <- clean_column_names(datasheet_spn)
  
  # match seqid across all dfs
  coverageinfo_tb$Seq_ID <- as.character(coverageinfo_tb$Seq_ID)
  datasheet_tb$Seq_ID <- as.character(datasheet_tb$Seq_ID)
  coverageinfo_spn_amp$Seq_ID <- as.character(coverageinfo_spn_amp$Seq_ID)
  coverageinfo_spn_mngs$Seq_ID <- as.character(coverageinfo_spn_mngs$Seq_ID)
  datasheet_spn$Seq_ID <- as.character(datasheet_spn$Seq_ID)
  
  # merge TB data 
  data_combined_tb <- coverageinfo_tb %>%
    left_join(datasheet_tb, by = 'Seq_ID') %>%
    mutate(pathogen = "TB")
  
  # Clean column names and values
  coverageinfo_spn_amp <- coverageinfo_spn_amp %>%
    rename_with(~ str_remove_all(., "^\"|\"$")) %>%
    mutate(
      Seq_ID = str_remove_all(Seq_ID, "^\"|\"$"),
      meanmapq = as.numeric(str_remove_all(meanmapq, "^\"|\"$")), 
      subsample = as.character(subsample) 
    )
  
  coverageinfo_spn_mngs <- coverageinfo_spn_mngs %>%
    rename_with(~ str_remove_all(., "^\"|\"$")) %>%
    mutate(
      Seq_ID = str_remove_all(Seq_ID, "^\"|\"$"),
      meanmapq = as.numeric(str_remove_all(meanmapq, "^\"|\"$")), 
      subsample = as.character(subsample) 
    )
  
  
  # Perform the joins
  data_combined_spn <- coverageinfo_spn_amp %>%
    full_join(coverageinfo_spn_mngs, by = 'Seq_ID') %>%
    mutate(sample = str_c("Yale-", Seq_ID)) %>%
    mutate(Seq_ID = str_remove(Seq_ID, "^Yale-")) %>%
    left_join(datasheet_spn, by = 'Seq_ID') %>%
    mutate(pathogen = "SPn")
  
  list(data_combined_tb = data_combined_tb, data_combined_spn = data_combined_spn)
}

# example usage
result <- cleanup_data(SP_amp_cov, SP_mNGS_cov, TB_amp_cov, SP_metadat, TB_metadat)
data_combined_tb <- result$data_combined_tb
data_combined_spn <- result$data_combined_spn













# visualize data ------------------------------------------------------
visualize <- function(merged_df) {
  
  # visualize coverage data 
  plot <- merged_df %>%
    filter(Primer_Conc == "100uM" & subsample == "1.0" & (Sample_Dilution == "1e-03" | !is.na(Sample_Dilution | !is.na(Ct1)))) %>%
    ggplot() +
    geom_point(aes(Ct1, coverage, col = Sample_Type), size = 4) + 
    ylab("Coverage(%)") + xlab ("Ct") + 
    geom_hline(yintercept = 80, linetype = "dashed", color = "salmon") +
    geom_text(x = max(merged_df$Ct1), y = 80, label = "80%", vjust = 2, hjust = 2, color = "black", size = 5) + 
    scale_y_continuous(limits = c(-15, 105), breaks = c(0, 25, 50, 75, 100)) + 
    scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) + 
    scale_color_manual(values = wes_palette("Darjeeling2", 8, type = "continuous"))
  
  plot
}

coverage_files <- c(SP_amp_cov, SP_mNGS_cov)
metadata_file <- metadat

cleaned_data <- cleanup_data(coverage_files, metadata_file)
merged_df <- cleaned_data$merged_df

# generate and display the plot
coverage_plot <- visualize(merged_df)
print(coverage_plot)
