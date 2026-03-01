# main_analysis.R
source("R/utils_io.R")

# 1. Setup the folder where results go
output_directory <- setup_output_dir()

# 2. Get the input file(s)
# (Using the loop version you requested)
files_to_process <- get_multiple_inputs()

# Process the first file selected
if (length(files_to_process) > 0) {
  current_file <- files_to_process[1]
  df <- read.csv(current_file, stringsAsFactors = FALSE)
  
  # 3. Display Headers
  cat("\nSUCCESS: File Loaded.\n")
  cat("Headers found in this file:\n")
  print(colnames(df))
  
  # 4. Ask which column to query
  target <- readline(prompt = "Which column would you like to query? ")
  
  if (target %in% colnames(df)) {
    # 5. Calculate unique entries and frequency
    counts_table <- as.data.frame(table(df[[target]]))
    colnames(counts_table) <- c("Unique_Entry", "Count")
    
    # Sort by frequency (highest at top)
    counts_table <- counts_table[order(counts_table$Count, decreasing = TRUE), ]
    
    # 6. Show results and Save
    print(counts_table)
    
    out_file <- file.path(output_directory, paste0("query_", target, ".csv"))
    write.csv(counts_table, out_file, row.names = FALSE)
    
    cat("\nSUCCESS: Unique counts saved to:", out_file, "\n")
  } else {
    cat("\nERROR: Column name not found. Please check spelling/capitalization.\n")
  }
}