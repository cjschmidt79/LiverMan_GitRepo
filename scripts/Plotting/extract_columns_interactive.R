# ------------------------------------------------------------
# ============================================================
# Script: extract_columns_interactive.R
#
# Description:
# This interactive R script allows a user to extract specific
# columns from a "column table" (e.g., data matrix) based on
# values listed in the rows of another "row table" (e.g., list
# of variable names like genes or metabolites).
#
# ▶️ Main Functionality:
# - The script takes the contents in the rows of one table
#   (e.g., a column containing gene names) and extracts those
#   columns from the second table into a new data frame.
#
# - It also allows the user to select additional metadata
#   columns from the column table (e.g., Sample ID, Treatment)
#   to be prepended to the output.
#
# 🧭 Interactive Steps:
# 1. Choose the row table CSV (contains values for extraction).
# 2. Choose the column table CSV (contains data to extract from).
# 3. View the first 5 rows of the row table to help identify the
#    correct column for extraction.
# 4. Enter the column name in the row table that contains the
#    values (e.g., "Variable", "Gene", etc.).
# 5. View the first 10 columns of the column table and select
#    any you want to include (e.g., "SampleID").
# 6. Choose the output folder and filename (a default is proposed).
# 7. The result is saved as a CSV and previewed in the console.
#
# ✅ Example use case:
# Extract columns for a specific gene signature from a gene
# expression matrix, and include sample metadata in the result.
#
# Author: Carl Schmidt
# Date:05/14/2025
# ============================================================

# ------------------------------------------------------------

library(dplyr)
library(tools)

if (interactive()) {
    # === FILE SELECTION ========================================
    cat("\n📁 Please select your row extraction CSV file:\n")
    row_file <- file.choose()
    row_df <- read.csv(row_file, stringsAsFactors = FALSE)
    
    cat("\n📁 Please select your column extraction CSV file:\n")
    col_file <- file.choose()
    col_df <- read.csv(col_file, stringsAsFactors = FALSE)
    
    row_basename <- file_path_sans_ext(basename(row_file))
    col_basename <- file_path_sans_ext(basename(col_file))
    
    # === PREVIEW ROW TABLE =====================================
    cat("\n🔍 First 5 rows of the row extraction table:\n")
    print(head(row_df, 5))
    
    # === ASK FOR COLUMN NAME FOR EXTRACTION ====================
    row_col_name <- readline(prompt = "\n✏️  Enter the name of the column in the row table that contains column names to extract from the column table: ")
    
    if (!(row_col_name %in% colnames(row_df))) {
        stop(paste0("❌ Column '", row_col_name, "' not found in the row table."))
    }
    
    # === ASK WHICH col_df COLUMNS TO APPEND ====================
    cat("\n📋 First 10 columns in the column table:\n")
    for (i in seq_len(min(10, ncol(col_df)))) {
        cat(paste0("[", i, "] ", colnames(col_df)[i], "\n"))
    }
    
    append_input <- readline(prompt = paste0(
        "\n➕ Enter column numbers from the column table to append to the output.\n",
        "   👉 Example: 1,3 for columns 1 and 3 | 2:5 for columns 2 to 5\n",
        "   [Press Enter to select column 1 by default]: "
    ))
    
    append_cols <- if (append_input == "") {
        1
    } else {
        tryCatch(eval(parse(text = paste0("c(", append_input, ")"))),
                 error = function(e) stop("❌ Invalid input format for column numbers."))
    }
    
    if (any(append_cols < 1 | append_cols > ncol(col_df))) {
        stop("❌ One or more selected columns to append are out of range.")
    }
    
    col_metadata <- col_df[, append_cols, drop = FALSE]
    
    # === PERFORM EXTRACTION =====================================
    columns_to_extract <- unique(row_df[[row_col_name]])
    columns_to_keep <- intersect(columns_to_extract, colnames(col_df))
    missing_columns <- setdiff(columns_to_extract, colnames(col_df))
    
    if (length(missing_columns) > 0) {
        warning("⚠️ Some requested columns were not found in the column table:\n", paste(missing_columns, collapse = ", "))
    }
    
    if (length(columns_to_keep) == 0) {
        stop("❌ No matching columns found. Nothing to extract.")
    }
    
    extracted_df <- col_df %>% select(all_of(columns_to_keep))
    cat("\n✅ Extracted", length(columns_to_keep), "columns from the column table.\n")
    
    # === COMBINE FINAL OUTPUT ===================================
    final_df <- cbind(col_metadata, extracted_df)
    
    # === OUTPUT DIRECTORY & FILE ================================
    cat("\n📂 Please select output directory (choose any file inside it):\n")
    output_dir <- dirname(file.choose())
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("📁 Output directory created:", output_dir, "\n")
    }
    
    suggested_filename <- paste0(row_basename, "_", col_basename, ".csv")
    user_filename <- readline(prompt = paste0("\n💾 Enter output filename [default: ", suggested_filename, "]: "))
    if (user_filename == "") {
        user_filename <- suggested_filename
    }
    
    output_path <- file.path(output_dir, user_filename)
    write.csv(final_df, output_path, row.names = FALSE)
    cat("\n✅ Output saved to:\n", output_path, "\n")
    
    # === FINAL PREVIEW ==========================================
    cat("\n🔍 Preview of final output:\n")
    print(head(final_df))
}
