# HPC R Script: process_and_save.R
library(tidyr)

# Function to read and save data
read_and_save <- function(file_base) {
  dir <- "MSD_analysis/"
  
  # Read the data
  t <- read.table(paste0(dir, file_base, "_t_values.txt"), header = FALSE)
  msd_total <- read.table(paste0(dir, file_base, "_msd_total.txt"), header = FALSE)
  msd_av <- read.table(paste0(dir, file_base, "_msd_av.txt"), header = FALSE)
  
  # Create a list containing all data
  msd_data <- list(t = t, msd_total = msd_total, msd_av = msd_av)
  
  # Save the list to a .rda file
  save(msd_data, file = paste0(dir, "Pure",file_base, ".rda"))
}

# Process each file
files <- c("MSDtestCOMx", "MSDtestCOMy", "MSDtestCOMz")
lapply(files, read_and_save)
