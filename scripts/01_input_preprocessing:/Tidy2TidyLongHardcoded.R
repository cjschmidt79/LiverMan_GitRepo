############################################################
# Standardize the uploaded "TIDY_Transpose of ImmuneGenes.csv"
# (It is already tidy-long: Day, SampleID, Gene, fpkm)
# Base R only | Safe to run via RStudio "Source"
############################################################

standardize_immune_tidy_long <- function(input_file = NULL, output_dir = NULL,
                                         parse_rep_from_sampleid = FALSE) {
  
  message("=== Standardize Immune Tidy-Long Table ===")
  
  # ---- Select input file ----
  if (is.null(input_file)) input_file <- file.choose()
  if (!file.exists(input_file)) stop("Input file does not exist:\n", input_file)
  
  # ---- Output directory ----
  if (is.null(output_dir)) output_dir <- dirname(input_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Read CSV ----
  df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # ---- Validate expected columns (your uploaded file has these exact names) ----
  required <- c("Day", "SampleID", "Gene", "fpkm")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         "\nFound columns: ", paste(colnames(df), collapse = ", "))
  }
  
  # ---- Coerce types ----
  df$Day <- suppressWarnings(as.integer(df$Day))
  df$fpkm <- suppressWarnings(as.numeric(df$fpkm))
  
  # ---- Standardize column names (canonical long schema) ----
  out <- data.frame(
    Day      = df$Day,
    SampleID = df$SampleID,
    Feature  = df$Gene,
    Value    = df$fpkm,
    stringsAsFactors = FALSE
  )
  
  # ---- Optional: parse replicate from SampleID (best-effort, non-destructive) ----
  # This tries to detect patterns like "14_3" or "..._rep3" etc.
  if (isTRUE(parse_rep_from_sampleid)) {
    # Pattern 1: "..._<rep>" where rep is trailing digits after underscore
    rep1 <- sub("^.*_([0-9]+)$", "\\1", out$SampleID)
    ok1 <- grepl("^[0-9]+$", rep1)
    
    # Pattern 2: "...rep<digits>..."
    rep2 <- sub("^.*[Rr][Ee][Pp]([0-9]+).*$", "\\1", out$SampleID)
    ok2 <- grepl("^[0-9]+$", rep2)
    
    Replicate <- rep(NA_integer_, nrow(out))
    Replicate[ok1] <- as.integer(rep1[ok1])
    Replicate[!ok1 & ok2] <- as.integer(rep2[!ok1 & ok2])
    
    out$Replicate <- Replicate
  }
  
  # ---- Write outputs ----
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base <- tools::file_path_sans_ext(basename(input_file))
  
  out_file <- file.path(output_dir, paste0(base, "_STANDARD_LONG_", ts, ".csv"))
  write.csv(out, out_file, row.names = FALSE)
  
  message("Done.")
  message("Input rows:  ", nrow(df))
  message("Output rows: ", nrow(out))
  message("Wrote: ", out_file)
  
  invisible(list(long_table = out, output_file = out_file))
}

############################################################
# AUTO-RUN WHEN SOURCED
############################################################
if (interactive()) {
  res <- standardize_immune_tidy_long(parse_rep_from_sampleid = FALSE)
  try(View(res$long_table), silent = TRUE)
}
