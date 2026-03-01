#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
in_txt  <- if (length(args) >= 1) args[[1]] else "Systemic_Canalization.txt"
out_csv <- if (length(args) >= 2) args[[2]] else "Systemic_Canalization_FIXED.csv"

if (!file.exists(in_txt)) stop("Input .txt not found: ", in_txt)

lines <- readLines(in_txt, warn = FALSE, encoding = "UTF-8")

# Strip UTF-8 BOM if present on first line
if (length(lines) > 0) lines[1] <- sub("^\ufeff", "", lines[1])

# Records start at "Reference Type:" (allow leading whitespace)
starts <- which(grepl("^\\s*Reference Type:", lines))
if (length(starts) == 0) stop("No records found (no 'Reference Type:' lines).")

ends <- c(starts[-1] - 1, length(lines))

# Tag line detector: optional leading whitespace + "Field Name:" pattern
is_tag_line <- function(x) {
  grepl("^\\s*[A-Za-z][A-Za-z0-9 /\\-\\(\\)\\.]*:\\s*", x)
}

extract_field <- function(rec_lines, key) {
  # Find line that starts this key (allow leading whitespace)
  idx <- which(grepl(paste0("^\\s*", key, ":"), rec_lines))
  if (length(idx) == 0) return(NA_character_)
  i <- idx[1]
  
  # Value begins after "Key:"
  first <- sub(paste0("^\\s*", key, ":\\s*"), "", rec_lines[i])
  
  # Continuation lines: until next tag line
  out <- first
  if (i < length(rec_lines)) {
    for (j in (i + 1):length(rec_lines)) {
      ln <- rec_lines[j]
      if (is_tag_line(ln)) break
      out <- paste(out, ln, sep = "\n")
    }
  }
  
  val <- trimws(out)
  
  # Keep newlines internally, but normalize excessive spaces on each line
  val <- paste(vapply(strsplit(val, "\n", fixed = TRUE)[[1]],
                      function(z) gsub("[ \t]+", " ", trimws(z)),
                      character(1)), collapse = "\n")
  
  if (!nzchar(val)) NA_character_ else val
}

recs <- vector("list", length(starts))

for (k in seq_along(starts)) {
  rec_lines <- lines[starts[k]:ends[k]]
  
  recs[[k]] <- data.frame(
    RecordNumber    = extract_field(rec_lines, "Record Number"),
    ReferenceType   = extract_field(rec_lines, "Reference Type"),
    Title           = extract_field(rec_lines, "Title"),
    Abstract        = extract_field(rec_lines, "Abstract"),
    Year            = extract_field(rec_lines, "Year"),
    Journal         = extract_field(rec_lines, "Journal"),
    DOI             = extract_field(rec_lines, "DOI"),
    Authors         = extract_field(rec_lines, "Author"),
    Keywords        = extract_field(rec_lines, "Keywords"),
    AccessionNumber = extract_field(rec_lines, "Accession Number"),
    URL             = extract_field(rec_lines, "URL"),
    stringsAsFactors = FALSE
  )
}

df <- do.call(rbind, recs)

cat("Parsed records: ", nrow(df), "\n", sep = "")
cat("Titles present: ", sum(!is.na(df$Title) & nzchar(df$Title)), "\n", sep = "")
cat("Abstracts present: ", sum(!is.na(df$Abstract) & nzchar(df$Abstract)), "\n", sep = "")

write.csv(df, out_csv, row.names = FALSE)
cat("Wrote CSV: ", normalizePath(out_csv, winslash = "/", mustWork = FALSE), "\n", sep = "")
