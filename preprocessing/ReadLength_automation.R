BiocManager::install("ShortRead")

library(ShortRead)

# This script calculates the maximum read length from the FASTQ file. The value is used to determine the --sjdbOverhang parameter during the STAR genome index generation step. 

# Function to determine maximum read length
determine_read_length <- function(fastq_file) {
  # Load the FASTQ file
  fq <- readFastq(fastq_file)
  
  # Extract sequences from the fastq object
  sequences <- sread(fq)
  
  # Get the length of sequences
  lengths <- width(sequences)
  
  # Determine the maximum length
  max_length <- max(lengths)
  
  return(max_length)
}

# Directory containing the FASTQ files
fastq_directory <- "where/fastq/files/are/located"

# Use list.files to get all files that start with 'SRR' and end with '.fastq' (when fastq begins with SRR)
# fastq_files <- list.files(path = fastq_directory, pattern = "^SRR.*\\.fastq$", full.names = TRUE)


# Use list.files to get all files that end with '.fastq' (for others where fastq does not begin with SRR)
# (This also works for SRR-beginnig file when all the files in the directory begins with SRR and ends with fsatq)
fastq_files <- list.files(path = fastq_directory, pattern = "\\.fastq$", full.names = TRUE)


# Loop through each file and calculate the maximum read length
results_for_this <- setNames(lapply(fastq_files, determine_read_length), fastq_files)


# Print the results
print(results_for_this)
