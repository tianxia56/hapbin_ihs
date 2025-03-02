# Install necessary packages if they are not already installed
if (!require("dplyr")) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if (!require("data.table")) install.packages("data.table", repos = "http://cran.us.r-project.org")

# Load necessary libraries
library(dplyr)
library(data.table)

# Define the populations
populations <- c("IND", "MGN", "SASI", "SASP", "PJL", "KOR")

# Create the output directory
output_dir <- "/home/tx56/palmer_scratch/100kga/ihs/AA_bfile"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read all TSV files once and store them in a list
tsv_data_list <- lapply(1:22, function(i) {
  tsv_file <- sprintf("/vast/palmer/pi/reilly/jfa38/datasets_for_annotation/ErinG_ARG_hg19_AncestralAlleles/tidied_AA_tsvs/chr%d_ErinG_AA_hg19.tsv.gz", i)
  tsv_data <- fread(cmd = sprintf("zcat %s", tsv_file), header = TRUE)
  tsv_data <- tsv_data %>%
    mutate(POS = as.integer(POS))
  return(tsv_data)
})

# Loop through each population and chromosome file
for (pop in populations) {
  for (i in 1:22) {
    # Define input and output file names
    bim_file <- sprintf("/home/tx56/palmer_scratch/100kga/ihs/qc_bfile/%s.%d.qc1.bim", pop, i)
    output_file <- sprintf("%s/%s.%d.AA.bim", output_dir, pop, i)
    
    # Read the BIM file
    bim_data <- fread(bim_file, header = FALSE)
    colnames(bim_data) <- c("chr", "rsid", "cm", "pos", "allele1", "allele2")
    
    # Get the corresponding TSV data
    tsv_data <- tsv_data_list[[i]]
    
    # Filter rows by matching pos
    filtered_data <- bim_data %>%
      filter(pos %in% tsv_data$POS)
    
    # Ensure either allele matches ErinG_ARG_AA
    correct_alleles <- tsv_data$ErinG_ARG_AA[match(filtered_data$pos, tsv_data$POS)]
    
    # Create a logical vector to identify rows where allele1 should be swapped
    swap_alleles <- filtered_data$allele2 == correct_alleles
    
    # Swap alleles if necessary
    filtered_data[swap_alleles, c("allele1", "allele2")] <- filtered_data[swap_alleles, c("allele2", "allele1")]
    
    # Remove rows where neither allele matches AA
    filtered_data <- filtered_data %>%
      filter(allele1 %in% correct_alleles)
    
    # Remove duplicate rsid rows
    filtered_data <- filtered_data %>%
      distinct(rsid, .keep_all = TRUE)
    
    # Write the filtered data to the output file
    fwrite(filtered_data, output_file, sep = "\t", col.names = FALSE)
    
    # Print the number of rows from the TSV file, and the number of rows kept and removed
    total_bim_rows <- nrow(bim_data)
    kept_rows <- nrow(filtered_data)
    removed_rows <- total_bim_rows - kept_rows
    cat(sprintf("Population: %s, Chromosome: %d\n", pop, i))
    cat(sprintf("Total rows from TSV: %d\n", nrow(tsv_data)))
    cat(sprintf("Rows kept: %d\n", kept_rows))
    cat(sprintf("Rows removed: %d\n\n", removed_rows))
  }
}
