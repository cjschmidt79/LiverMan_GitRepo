#!/usr/bin/env Rscript
# ======================================================================
# WIDE -> TIDY (LONG) transcriptome converter (Mac + RStudio-friendly)
#
# INPUT (WIDE):
#   - Column 1: Day
#   - Column 2: SampleID
#   - Columns 3..N: genes (column names = gene symbols), values = abundance
#
# OUTPUT (TIDY / LONG):
#   Day, SampleID, Gene, Abundance
#
# Features:
#   - RStudio-native file chooser for input when available (most reliable on Mac)
#   - Folder chooser for output (RStudio-native when available; tcltk fallback)
#   - Output name: TIDY_<input_basename>.csv
#   - Optional low-expression filtering:
#       * User provides mean-abundance threshold across all rows.
#       * Default = 1; enter 0 for NO FILTER.
#   - Writes a log file documenting parameters and row/gene counts
#
# Notes on "Mac-style directory":
#   - You will be prompted to choose a folder via a Mac/RStudio dialog when possible.
#   - If a dialog is unavailable, you must paste an absolute macOS path like:
#       /Users/yourname/Desktop/output_folder
# ======================================================================

# -----------------------------
# Helpers
# -----------------------------
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

safe_basename_noext <- function(path) {
  b <- basename(path)
  sub("\\.[^.]*$", "", b)
}

stop2 <- function(...) stop(paste0(...), call. = FALSE)

is_rstudio <- function() {
  requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
}

# Ask user for an input file (RStudio picker preferred; else file.choose)
choose_input_file <- function() {
  message("Select the input WIDE table (CSV/TSV)...")
  flush.console()
  
  f <- ""
  if (is_rstudio()) {
    f <- tryCatch(rstudioapi::selectFile(caption = "Select input WIDE table"),
                  error = function(e) "")
  } else {
    f <- tryCatch(file.choose(), error = function(e) "")
  }
  
  f <- trimws(f)
  if (!nzchar(f)) stop2("No input file selected.")
  if (!file.exists(f)) stop2("Selected file does not exist: ", f)
  
  message("Input file: ", f)
  flush.console()
  f
}

# Ask user for an output directory (RStudio folder picker preferred; else tcltk; else paste path)
choose_output_dir <- function() {
  message("Choose output folder for the tidy table and log...")
  flush.console()
  
  out <- ""
  if (is_rstudio()) {
    out <- tryCatch(rstudioapi::selectDirectory(caption = "Choose output folder"),
                    error = function(e) "")
  }
  
  if (!nzchar(out) && requireNamespace("tcltk", quietly = TRUE)) {
    out <- tryCatch(tcltk::tclvalue(tcltk::tk_choose.dir(caption = "Choose output folder")),
                    error = function(e) "")
  }
  
  out <- trimws(out)
  
  if (!nzchar(out)) {
    message("\nCould not open a folder chooser.")
    message("Please paste a macOS folder path to save outputs, for example:")
    message("  /Users/yourname/Desktop/output_folder")
    out <- trimws(readline(prompt = "> "))
  }
  
  if (!nzchar(out)) stop2("No output directory provided.")
  if (!dir.exists(out)) {
    # Create it if it does not exist
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(out)) stop2("Could not create output directory: ", out)
  
  message("Output folder: ", out)
  flush.console()
  out
}

read_table_auto <- function(path) {
  # Auto-detect delimiter by extension; fallback tries both CSV and tab-delimited
  ext <- tolower(sub(".*\\.", "", path))
  if (ext %in% c("tsv", "txt")) {
    df <- utils::read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    df <- tryCatch(
      utils::read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(df)) {
      df <- utils::read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  df
}

prompt_numeric <- function(prompt, default) {
  message(prompt)
  message("Press Enter for default: ", default)
  x <- trimws(readline(prompt = "> "))
  if (!nzchar(x)) return(as.numeric(default))
  val <- suppressWarnings(as.numeric(x))
  if (!is.finite(val)) {
    message("Invalid numeric input; using default: ", default)
    return(as.numeric(default))
  }
  val
}

write_log <- function(log_path, lines) {
  cat(paste0(lines, collapse = "\n"), file = log_path)
}

# -----------------------------
# Main
# -----------------------------
message("=== STARTING WIDE->TIDY SCRIPT: reached Main section ===")
flush.console()

infile <- choose_input_file()
outdir <- choose_output_dir()

base <- safe_basename_noext(infile)
outfile <- file.path(outdir, paste0("TIDY_", base, ".csv"))
logfile <- file.path(outdir, paste0("TIDY_", base, "_LOG_", ts_stamp(), ".txt"))

log_lines <- c(
  "WIDE -> TIDY (LONG) converter",
  paste0("Timestamp: ", Sys.time()),
  paste0("Input file: ", infile),
  paste0("Output dir: ", outdir),
  paste0("Planned output: ", outfile),
  ""
)

# Read
df <- read_table_auto(infile)
if (!is.data.frame(df) || nrow(df) == 0) stop2("Input file read failed or is empty: ", infile)

# Validate columns
if (ncol(df) < 3) stop2("Input table must have at least 3 columns: Day, SampleID, and >=1 gene column.")

colnames_in <- names(df)
if (!identical(colnames_in[1:2], c("Day", "SampleID"))) {
  log_lines <- c(log_lines, "NOTE: First two columns were not exactly 'Day' and 'SampleID'. Renaming column 1 -> Day, column 2 -> SampleID.")
  names(df)[1:2] <- c("Day", "SampleID")
}

# Gene columns begin at the 3rd column
gene_cols <- names(df)[3:ncol(df)]
if (length(gene_cols) < 1) stop2("No gene columns detected (expected columns 3..N).")

if (anyDuplicated(gene_cols)) {
  log_lines <- c(log_lines, "WARNING: Duplicate gene column names detected. Making them unique with suffixes.")
  gene_cols_unique <- make.unique(gene_cols)
  names(df)[3:ncol(df)] <- gene_cols_unique
  gene_cols <- gene_cols_unique
}

# Coerce types
df$Day <- suppressWarnings(as.integer(df$Day))
if (any(!is.finite(df$Day))) {
  log_lines <- c(log_lines, "WARNING: 'Day' column could not be fully coerced to integer. Non-integer values may appear as NA.")
}
df$SampleID <- as.character(df$SampleID)


# Convert to tidy (long)
tidy <- NULL

if (requireNamespace("data.table", quietly = TRUE)) {
  DT <- data.table::as.data.table(df)
  tidy <- data.table::melt(
    DT,
    id.vars = c("Day", "SampleID"),
    measure.vars = gene_cols,
    variable.name = "Gene",
    value.name = "Abundance",
    variable.factor = FALSE
  )
  tidy[, Abundance := suppressWarnings(as.numeric(Abundance))]
  tidy <- tidy[is.finite(Abundance)]
} else {
  stacked <- stack(df[gene_cols])
  tidy <- data.frame(
    Day = rep(df$Day, times = length(gene_cols)),
    SampleID = rep(df$SampleID, times = length(gene_cols)),
    Gene = as.character(stacked$ind),
    Abundance = suppressWarnings(as.numeric(stacked$values)),
    stringsAsFactors = FALSE
  )
  tidy <- tidy[is.finite(tidy$Abundance), ]
}

n_rows_pre_filter <- if (inherits(tidy, "data.table")) nrow(tidy) else nrow(tidy)
n_genes_pre_filter <- length(unique(if (inherits(tidy, "data.table")) tidy$Gene else tidy$Gene))

# --- ADD THIS BEFORE THE FILTER SECTION ---
filter_threshold <- prompt_numeric(
  prompt = "\nLow-expression gene filter:\nEnter a MEAN abundance threshold (computed across all rows) for removing low-expression genes.\n- Default = 1\n- Enter 0 to apply NO FILTER",
  default = 1
)
# ------------------------------------------
# --- FIXED FILTERING SECTION ---

# 1. Identify which columns are actually genes (columns 3 to N)
# We use gene_cols which was defined earlier in your script
gene_data <- df[, gene_cols, drop = FALSE]

# 2. Ensure gene data is numeric (strips out any stray characters)
gene_data[] <- lapply(gene_data, function(x) as.numeric(as.character(x)))

# 3. Calculate the mean for each GENE across all SAMPLES (column means)
col_means <- colMeans(gene_data, na.rm = TRUE)

# 4. Identify genes that pass the threshold
keep_genes_idx <- which(col_means >= filter_threshold)
keep_genes_names <- names(keep_genes_idx)

# 5. Update the data frame to keep Day, SampleID, and only the 'keep' genes
df_filtered <- df[, c("Day", "SampleID", keep_genes_names), drop = FALSE]

# 6. Log the results
genes_removed <- length(gene_cols) - length(keep_genes_names)
log_lines <- c(log_lines, 
               paste0("Filter threshold applied: ", filter_threshold),
               paste0("Genes before filtering: ", length(gene_cols)),
               paste0("Genes removed: ", genes_removed),
               paste0("Genes remaining: ", length(keep_genes_names)))

# 7. Update df and gene_cols for the rest of the script
df <- df_filtered
gene_cols <- keep_genes_names
# -------------------------------

n_rows_post_filter <- if (inherits(tidy, "data.table")) nrow(tidy) else nrow(tidy)
n_genes_post_filter <- length(unique(if (inherits(tidy, "data.table")) tidy$Gene else tidy$Gene))

log_lines <- c(
  log_lines,
  "",
  paste0("Input dimensions: rows=", nrow(df), " cols=", ncol(df)),
  paste0("Gene columns detected: ", length(gene_cols)),
  paste0("Tidy rows (pre-filter): ", n_rows_pre_filter),
  paste0("Unique genes (pre-filter): ", n_genes_pre_filter),
  paste0("Tidy rows (post-filter): ", n_rows_post_filter),
  paste0("Unique genes (post-filter): ", n_genes_post_filter)
)

# Ensure exact output columns and order; write
if (inherits(tidy, "data.table")) {
  tidy_out <- tidy[, .(Day, SampleID, Gene, Abundance)]
  data.table::setorder(tidy_out, Day, SampleID, Gene)
  data.table::fwrite(tidy_out, outfile)
} else {
  tidy_out <- tidy[, c("Day", "SampleID", "Gene", "Abundance")]
  tidy_out <- tidy_out[order(tidy_out$Day, tidy_out$SampleID, tidy_out$Gene), ]
  utils::write.csv(tidy_out, outfile, row.names = FALSE)
}

log_lines <- c(
  log_lines,
  "",
  paste0("Wrote tidy output: ", outfile),
  paste0("Wrote log: ", logfile)
)
write_log(logfile, log_lines)

message("\nDone.")
message("Tidy output: ", outfile)
message("Log file:    ", logfile)
flush.console()
