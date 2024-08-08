rm(list = ls())

library(tidyverse)
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)
library(readr)
library(wesanderson)

SP_amp_cov  <- "data/SP_amplicon_seq.csv"
SP_mNGS_cov <- "data/SP_mNGS_seq.csv"
TB_amp_cov  <- "data/TB_combined_coverage.tsv"

SP_metadat  <- "metadata/PGCOE_DataSheet - S. pneumoniae.csv"
TB_metadat  <- "metadata/PGCOE_Datasheet - M.tuberculosis_Seq.csv"

outpng      <- "Fig3_combined_coverage_plot.png"


# s.pneumo sample type coverage plot -------------------------------------------

process_data <- function(metadata_path, amp_cov_path, mngs_cov_path) {
  
  # read in files
  SP_metadat <- read.csv(metadata_path)
  SP_amp_cov <- read_table(amp_cov_path)
  SP_mngs_cov <- read_table(mngs_cov_path)
  
  # clean up files for easier merging
  SP_amp_cov <- SP_amp_cov %>%
    rename_with(~ str_remove_all(., "^\"|\"$")) %>%
    mutate(
      Seq_ID = str_remove_all(Seq_ID, "^\"|\"$"),
      meanmapq = as.numeric(str_remove_all(meanmapq, "^\"|\"$")),
      Seq_ID = str_replace(Seq_ID, "0000", "000"),
      subsample = as.character(subsample)  
    )
  
  SP_mngs_cov <- SP_mngs_cov %>%
    rename_with(~ str_remove_all(., "^\"|\"$")) %>%
    mutate(
      Seq_ID = str_remove_all(Seq_ID, "^\"|\"$"),
      meanmapq = as.numeric(str_remove_all(meanmapq, "^\"|\"$")),
      subsample = as.character(subsample) 
    )
  
  SP_metadat <- SP_metadat %>%
    mutate(Qubit_Pool2 = str_trim(Qubit_Pool2)) %>%  
    filter(str_to_lower(Qubit_Pool2) != "contaminated") %>%
    mutate(Seq_ID = str_remove(Seq_ID, "^Yale-"))
  
  # merge amp, mngs, and pgcoe datasheet
  SP_merged <- SP_metadat %>%
    full_join(SP_amp_cov, by = "Seq_ID") %>%
    full_join(SP_mngs_cov, by = "Seq_ID") %>%
    rename(
      subsample_x = subsample.x,
      rname_x = rname.x,
      startpos_x = startpos.x,
      endpos_x = endpos.x,
      numreads_x = numreads.x,
      covbases_x = covbases.x,
      coverage_x = coverage.x, 
      meandepth_x = meandepth.x,
      meanbaseq_x = meanbaseq.x,
      meanmapq_x = meanmapq.x,
      subsample_y = subsample.y,
      rname_y = rname.y,
      startpos_y = startpos.y,
      endpos_y = endpos.y,
      numreads_y = numreads.y,
      covbases_y = covbases.y,
      coverage_y = coverage.y, 
      meandepth_y = meandepth.y,
      meanbaseq_y = meanbaseq.y,
      meanmapq_y = meanmapq.y
    ) %>%
    mutate(
      subsample = coalesce(subsample_x, subsample_y),
      rname = coalesce(rname_x, rname_y),
      startpos = coalesce(startpos_x, startpos_y),
      endpos = coalesce(endpos_x, endpos_y),
      numreads = coalesce(numreads_x, numreads_y),
      covbases = coalesce(covbases_x, covbases_y),
      coverage = coalesce(coverage_x, coverage_y),
      meandepth = coalesce(meandepth_x, meandepth_y),
      meanbaseq = coalesce(meanbaseq_x, meanbaseq_y),
      meanmapq = coalesce(meanmapq_x, meanmapq_y)
    ) %>%
    select("Seq_ID", "Original_ID", "Sample_Type", "Sample_Dilution", "PCR_Primer_Scheme",
           "Primer_Conc", "NGS_Prep_Method", "NGS_Run_Date", "Qubit_Pool1", "Qubit_Pool2",      
           "Ct1", "Ct2", "Method", "subsample", "rname", "startpos", "endpos", "numreads", 
           "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq") %>%
    mutate(pathogen = "SPn")
  
  # return
  list(SP_metadat = SP_metadat, SP_amp_cov = SP_amp_cov, SP_mngs_cov = SP_mngs_cov, SP_merged = SP_merged)
}

plot_sample_types <- function(data) {
  
  # Filter the data
  SP_merged_filtered <- data %>%
    filter(Sample_Dilution == "0.001" & 
             Primer_Conc == "100uM" & 
             (subsample == "1" | subsample == "1.0") & 
             !is.na(Ct1) & 
             Seq_ID != "CS00163")
  
  # Plot the data
  ggplot(SP_merged_filtered, aes(x = Ct1, y = coverage, col = Method, group = interaction(Ct1, Sample_Type))) + 
    geom_point(alpha = 1, size = 4) +
    geom_line(color = "black") +
    scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
    facet_wrap(~Sample_Type) 
}  

# assemble combined plots ------------------------------------------------------

SP_coverage_plot <- process_data(SP_metadat, SP_amp_cov, SP_mNGS_cov)
plot_sample_types(SP_coverage_plot$SP_merged)

ggsave(outpng, width = 200, height = 200, dpi = 400, units = "mm")


# m.tuberculosis vs s.pneumo coverage plot ------------------------------------






