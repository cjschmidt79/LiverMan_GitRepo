#!/usr/bin/env Rscript
# ============================================================
# Z-score Function Scoring + Per-Group Trajectory Plots
# Script: Z-score_function.R
#
# PURPOSE
#   Given a gene-by-sample expression table with metadata columns and expression
#   columns named using a DAY_REP convention (e.g., "14_3"), this script:
#     1) Z-scores expression within each gene across ALL samples (gene-wise)
#     2) Computes per-gene, per-day mean Z and mean |Z| across replicates
#     3) Aggregates these gene-level day summaries within a user-selected
#        grouping column (e.g., Function / Pathway / Category)
#     4) Writes:
#        - Wide group table (all groups)
#        - Wide group table filtered by a user-defined gene-count cutoff N
#        - Long/tidy group table (group x day)
#        - Gene-level wide table (GeneMeanZ_DayX, GeneMeanAbsZ_DayX)
#        - Gene-level long/tidy table (gene x day)
#     5) Produces per-group plots (signed mean Z and mean |Z|) across days (600 dpi PNG)
#     6) Builds a Quarto QMD + HTML report in the SAME output directory.
#
# ASSUMPTIONS
#   - There is a column containing gene identifiers (often "Gene Symbol")
#   - Expression columns follow "DAY_REP" naming with an underscore, and the
#     day is the substring BEFORE the first underscore (e.g., "14" in "14_3")
#   - Expression columns are numeric (or coercible to numeric)
#
# NOTES
#   - No weighting is applied (each gene contributes equally within a group)
#   - n_genes counts UNIQUE gene identifiers per group
# ============================================================

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(utils)
  library(tools)
  library(stats)
  library(ggplot2)
  library(jsonlite)
})

# =========================
# Capture Script Metadata
# =========================
script_info <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    this_path <- rstudioapi::getSourceEditorContext()$path
  } else {
    this_path <- sys.frame(1)$ofile
  }
  if (is.null(this_path) || this_path == "") {
    list(name = "Z-score_function", path = getwd(), full = file.path(getwd(), "Z-score_function.R"))
  } else {
    list(
      name = tools::file_path_sans_ext(basename(this_path)),
      path = normalizePath(dirname(this_path), winslash = "/", mustWork = FALSE),
      full = normalizePath(this_path, winslash = "/", mustWork = FALSE)
    )
  }
}, error = function(e) {
  list(name = "Z-score_function", path = getwd(), full = file.path(getwd(), "Z-score_function.R"))
})

script_name <- script_info$name
script_path <- script_info$path
script_full <- script_info$full

# =========================
# Helper: Mac-friendly file chooser
# =========================
choose_input_file <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::selectFile(caption = "Choose input CSV file", filter = "CSV (*.csv)")
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }
  p <- file.choose()
  normalizePath(p, winslash = "/", mustWork = TRUE)
}

choose_output_parent_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    d <- rstudioapi::selectDirectory(caption = "Choose parent directory for outputs (optional)")
    if (!is.null(d) && nzchar(d)) return(normalizePath(d, winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

# =========================
# Helper: safe filenames
# =========================
safe_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\:]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "Unnamed"
  x
}

# =========================
# Interactive prompts
# =========================
cat("\n====================================\n")
cat("Z-score Group Scoring Pipeline\n")
cat("====================================\n\n")

cat("Choose input CSV file...\n")
input_path <- choose_input_file()
cat("Using input:\n  ", input_path, "\n\n", sep = "")

df <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
cols <- colnames(df)

# Show columns once
cat("Detected columns:\n")
for (i in seq_along(cols)) cat(sprintf("  [%d] %s\n", i, cols[i]))
cat("\n")

# =========================
# Choose gene identifier column (gene_col)
# =========================
default_gene_col <- if ("Gene Symbol" %in% cols) "Gene Symbol" else NA_character_

if (!is.na(default_gene_col)) {
  gene_col <- default_gene_col
  ans <- readline("Gene identifier column detected as 'Gene Symbol'. Use this? (y/n) [default=y]: ")
  if (nzchar(ans) && tolower(trimws(ans)) == "n") {
    idx <- suppressWarnings(as.integer(readline("Enter the column number that contains gene identifiers: ")))
    if (is.na(idx) || idx < 1 || idx > length(cols)) stop("Invalid gene column selection.")
    gene_col <- cols[idx]
  }
} else {
  cat("No default 'Gene Symbol' column detected.\n")
  idx <- suppressWarnings(as.integer(readline("Enter the column number that contains gene identifiers: ")))
  if (is.na(idx) || idx < 1 || idx > length(cols)) stop("Invalid gene column selection.")
  gene_col <- cols[idx]
}
cat("Using gene identifier column: ", gene_col, "\n\n", sep = "")

# =========================
# Choose grouping column (func_col)
# =========================
cat("Select the grouping column (category / pathway / function label).\n")
idx <- suppressWarnings(as.integer(readline("Enter the column number to use for grouping: ")))
if (is.na(idx) || idx < 1 || idx > length(cols)) stop("Invalid grouping column selection.")
func_col <- cols[idx]
cat("Using grouping column: ", func_col, "\n\n", sep = "")

# =========================
# Output directory (ALWAYS under the current working directory)
# =========================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Capture the current working directory ONCE and build outputs beneath it
wd_now <- getwd()
wd_now <- normalizePath(wd_now, winslash = "/", mustWork = TRUE)

cat("Current working directory:\n  ", wd_now, "\n\n", sep = "")

out_base_name <- readline("Enter a name for the output directory (e.g., Zscore_Groups): ")
out_base_name <- trimws(out_base_name)
if (!nzchar(out_base_name)) out_base_name <- "Zscore_Groups"

outputs_root <- file.path(wd_now, "outputs")
if (!dir.exists(outputs_root)) {
  dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)
}

output_dir <- file.path(outputs_root, paste0(out_base_name, "_", timestamp))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

cat("Output directory:\n  ", output_dir, "\n\n", sep = "")

# =========================
# Cutoff N
# =========================
cutoff_n <- suppressWarnings(as.integer(readline("Enter minimum number of UNIQUE gene IDs (N) to include in cutoff outputs: ")))
if (is.na(cutoff_n) || cutoff_n < 1) stop("Cutoff N must be a positive integer.")
cat("Using cutoff N = ", cutoff_n, "\n\n", sep = "")

# =========================
# Identify expression columns + infer days
# =========================
is_expr_name <- grepl("^[0-9]+_.+", cols)
expr_cols <- cols[is_expr_name]

if (length(expr_cols) == 0) {
  stop("No expression columns detected. Expected columns named like '14_3' (DAY_REP).")
}

# Coerce expression columns to numeric (keeping NAs)
for (cname in expr_cols) df[[cname]] <- suppressWarnings(as.numeric(df[[cname]]))

day_from_col <- function(cname) strsplit(cname, "_", fixed = TRUE)[[1]][1]
days_raw <- unique(vapply(expr_cols, day_from_col, character(1)))

days_num <- suppressWarnings(as.integer(days_raw))
if (any(is.na(days_num))) {
  bad_days <- days_raw[is.na(days_num)]
  stop(paste0("Could not parse numeric day from some columns: ", paste(bad_days, collapse = ", ")))
}
days <- as.character(sort(unique(days_num)))

day_cols <- lapply(days, function(d) expr_cols[grepl(paste0("^", d, "_"), expr_cols)])
names(day_cols) <- days

rep_counts <- vapply(day_cols, length, integer(1))
cat("Detected days and replicate column counts:\n")
for (d in days) cat(sprintf("  Day %s: %d columns\n", d, rep_counts[[d]]))
cat("\n")

# =========================
# Gene-wise Z-scoring across ALL samples
# =========================
X <- as.matrix(df[, expr_cols, drop = FALSE])

gene_mean <- rowMeans(X, na.rm = TRUE)
gene_sd   <- apply(X, 1, sd, na.rm = TRUE)
gene_sd[gene_sd == 0] <- NA_real_

Z <- sweep(X, 1, gene_mean, FUN = "-")
Z <- sweep(Z, 1, gene_sd,   FUN = "/")

# =========================
# Gene-level day summaries (wide matrices)
# =========================
gene_day_meanZ <- matrix(
  NA_real_, nrow = nrow(df), ncol = length(days),
  dimnames = list(NULL, paste0("GeneMeanZ_Day", days))
)
gene_day_meanAbsZ <- matrix(
  NA_real_, nrow = nrow(df), ncol = length(days),
  dimnames = list(NULL, paste0("GeneMeanAbsZ_Day", days))
)

for (i in seq_along(days)) {
  d <- days[i]
  cols_d <- day_cols[[d]]
  idx <- match(cols_d, expr_cols)
  gene_day_meanZ[, i]    <- rowMeans(Z[, idx, drop = FALSE], na.rm = TRUE)
  gene_day_meanAbsZ[, i] <- rowMeans(abs(Z[, idx, drop = FALSE]), na.rm = TRUE)
}

gene_level_wide <- data.frame(
  df[, setdiff(cols, expr_cols), drop = FALSE],
  gene_day_meanZ,
  gene_day_meanAbsZ,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================================
# POST-SUMMARY BLOCK:
#   gene_long -> group aggregation -> writes -> plots -> manifest/inventory
#   (Strictly de-duplicated; no reserved column names like `function` or `file`)
# ============================================================

# ---------- Guardrails (fail fast if run out-of-order)
needed_now <- c("df","cols","gene_col","func_col","output_dir","days","gene_day_meanZ","gene_day_meanAbsZ","gene_level_wide","cutoff_n","rep_counts")
missing_now <- needed_now[!vapply(needed_now, exists, logical(1), inherits = TRUE)]
if (length(missing_now) > 0) {
  stop("Internal error: this block was reached before required objects existed: ",
       paste(missing_now, collapse = ", "), call. = FALSE)
}

# =========================
# Gene-level tidy/long output
# =========================
gene_ids_vec    <- as.character(df[[gene_col]])
group_vals_vec  <- as.character(df[[func_col]])

gene_long <- do.call(rbind, lapply(seq_along(days), function(day_index) {
  day_label <- days[day_index]
  data.frame(
    gene_id         = gene_ids_vec,
    group_label     = group_vals_vec,
    day             = as.integer(day_label),
    gene_mean_z     = gene_day_meanZ[, day_index],
    gene_mean_abs_z = gene_day_meanAbsZ[, day_index],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
rownames(gene_long) <- NULL

# Remove rows with missing/blank group or gene id
gene_long_ok <- gene_long[
  !is.na(gene_long$group_label) & nzchar(trimws(gene_long$group_label)) &
    !is.na(gene_long$gene_id) & nzchar(trimws(gene_long$gene_id)),
  , drop = FALSE
]

# =========================
# Group-level aggregation (long)
# =========================
agg_mean_z <- aggregate(
  gene_mean_z ~ group_label + day,
  data = gene_long_ok,
  FUN = function(v) mean(v, na.rm = TRUE)
)
agg_mean_abs <- aggregate(
  gene_mean_abs_z ~ group_label + day,
  data = gene_long_ok,
  FUN = function(v) mean(v, na.rm = TRUE)
)

uniq_gene_count_by_group <- aggregate(
  gene_id ~ group_label,
  data = gene_long_ok,
  FUN = function(v) length(unique(v))
)
colnames(uniq_gene_count_by_group)[colnames(uniq_gene_count_by_group) == "gene_id"] <- "n_genes"

group_long <- merge(
  agg_mean_z,
  agg_mean_abs,
  by = c("group_label", "day"),
  all = TRUE
)
colnames(group_long)[colnames(group_long) == "gene_mean_z"] <- "mean_z"
colnames(group_long)[colnames(group_long) == "gene_mean_abs_z"] <- "mean_abs_z"

group_long <- merge(group_long, uniq_gene_count_by_group, by = "group_label", all.x = TRUE)
group_long$day <- as.integer(group_long$day)
group_long <- group_long[order(group_long$group_label, group_long$day), , drop = FALSE]

group_long_cutoff <- group_long[!is.na(group_long$n_genes) & group_long$n_genes >= cutoff_n, , drop = FALSE]

# =========================
# Group-level wide tables
# =========================
group_levels <- sort(unique(group_long$group_label))
group_wide <- data.frame(group_label = group_levels, stringsAsFactors = FALSE, check.names = FALSE)
group_wide <- merge(group_wide, uniq_gene_count_by_group, by = "group_label", all.x = TRUE)
group_wide$n_genes[is.na(group_wide$n_genes)] <- 0L

fill_day_col <- function(day_value, value_col_name, out_col_prefix) {
  day_value_int <- as.integer(day_value)
  sub_df <- group_long[group_long$day == day_value_int, c("group_label", value_col_name), drop = FALSE]
  out_col <- paste0(out_col_prefix, day_value_int)
  tmp <- merge(group_wide["group_label"], sub_df, by = "group_label", all.x = TRUE)
  group_wide[[out_col]] <<- tmp[[value_col_name]]
}

for (d in days) {
  fill_day_col(d, "mean_z", "mean_z_day")
  fill_day_col(d, "mean_abs_z", "mean_abs_z_day")
}

day_col_order <- c(
  unlist(lapply(days, function(d) paste0("mean_z_day", as.integer(d))), use.names = FALSE),
  unlist(lapply(days, function(d) paste0("mean_abs_z_day", as.integer(d))), use.names = FALSE)
)
day_col_order <- day_col_order[day_col_order %in% colnames(group_wide)]
group_wide <- group_wide[, c("group_label", "n_genes", day_col_order), drop = FALSE]

group_wide_cutoff <- group_wide[group_wide$n_genes >= cutoff_n, , drop = FALSE]

# =========================
# Write outputs (single, de-duplicated)
# =========================
fn_group_wide_all <- file.path(output_dir, "Z-scores_all.csv")
fn_group_wide_cut <- file.path(output_dir, "Z-scores_cutoff.csv")
fn_group_long     <- file.path(output_dir, "Z-scores_group_long.csv")
fn_group_long_cut <- file.path(output_dir, "Z-scores_group_long_cutoff.csv")
fn_gene_wide       <- file.path(output_dir, "Z-scores_gene_wide.csv")
fn_gene_long       <- file.path(output_dir, "Z-scores_gene_long.csv")

write.csv(group_wide,        fn_group_wide_all, row.names = FALSE)
write.csv(group_wide_cutoff, fn_group_wide_cut, row.names = FALSE)
write.csv(group_long,        fn_group_long,     row.names = FALSE)
write.csv(group_long_cutoff, fn_group_long_cut, row.names = FALSE)
write.csv(gene_level_wide,   fn_gene_wide,      row.names = FALSE)
write.csv(gene_long,         fn_gene_long,      row.names = FALSE)

cat("\nWrote tables:\n")
cat("  ", fn_group_wide_all, "\n", sep = "")
cat("  ", fn_group_wide_cut, "\n", sep = "")
cat("  ", fn_group_long, "\n", sep = "")
cat("  ", fn_group_long_cut, "\n", sep = "")
cat("  ", fn_gene_wide, "\n", sep = "")
cat("  ", fn_gene_long, "\n\n", sep = "")

# =========================
# Per-group plots (600 dpi PNG) — cutoff only
# =========================
plot_signed_dir <- file.path(output_dir, "plots_signed_meanZ_cutoff")
plot_abs_dir    <- file.path(output_dir, "plots_magnitude_meanAbsZ_cutoff")
dir.create(plot_signed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_abs_dir, recursive = TRUE, showWarnings = FALSE)

plot_df <- group_long_cutoff
plot_groups <- unique(plot_df$group_label)

cat("Generating per-group plots (cutoff only):\n")
cat("  Groups to plot: ", length(plot_groups), "\n\n", sep = "")

for (g in plot_groups) {
  dfg <- plot_df[plot_df$group_label == g, , drop = FALSE]
  dfg$day <- as.integer(dfg$day)
  
  p1 <- ggplot(dfg, aes(x = day, y = mean_z)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    labs(
      title = g,
      subtitle = sprintf("Signed mean Z across days (n_genes=%d)", unique(dfg$n_genes)),
      x = "Day",
      y = "mean_z"
    ) +
    theme_bw(base_size = 11)
  
  out1 <- file.path(plot_signed_dir, paste0("Group_", safe_filename(g), "_meanZ.png"))
  ggsave(filename = out1, plot = p1, width = 6.5, height = 4.0, dpi = 600, units = "in")
  
  p2 <- ggplot(dfg, aes(x = day, y = mean_abs_z)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    labs(
      title = g,
      subtitle = sprintf("Mean |Z| across days (n_genes=%d)", unique(dfg$n_genes)),
      x = "Day",
      y = "mean_abs_z"
    ) +
    theme_bw(base_size = 11)
  
  out2 <- file.path(plot_abs_dir, paste0("Group_", safe_filename(g), "_meanAbsZ.png"))
  ggsave(filename = out2, plot = p2, width = 6.5, height = 4.0, dpi = 600, units = "in")
}

cat("Plots written to:\n")
cat("  ", plot_signed_dir, "\n", sep = "")
cat("  ", plot_abs_dir, "\n\n", sep = "")

# =========================
# Manifest + inventory
# =========================
manifest_list <- list(
  script_name = "Z-score_function.R",
  detected_script_name = script_name,
  script_path = script_path,
  script_full = script_full,
  input_file = normalizePath(input_path, winslash = "/", mustWork = TRUE),
  output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
  timestamp = timestamp,
  gene_identifier_column = gene_col,
  grouping_column = func_col,
  day_parsing_rule = "Day = substring before first underscore in expression columns (DAY_REP)",
  cutoff_n_genes = cutoff_n,
  days_detected = days,
  replicate_cols_per_day = as.list(rep_counts),
  weighting = "none",
  tables_written = list(
    Z_scores_all = basename(fn_group_wide_all),
    Z_scores_cutoff = basename(fn_group_wide_cut),
    Z_scores_group_long = basename(fn_group_long),
    Z_scores_group_long_cutoff = basename(fn_group_long_cut),
    Z_scores_gene_wide = basename(fn_gene_wide),
    Z_scores_gene_long = basename(fn_gene_long)
  ),
  plot_dirs = list(
    signed_meanZ_cutoff = basename(plot_signed_dir),
    magnitude_meanAbsZ_cutoff = basename(plot_abs_dir)
  )
)

manifest_path <- file.path(output_dir, "Project_Manifest.json")
writeLines(jsonlite::toJSON(manifest_list, pretty = TRUE, auto_unbox = TRUE), manifest_path)

inventory_df <- data.frame(
  path = list.files(output_dir, recursive = TRUE, full.names = TRUE),
  stringsAsFactors = FALSE
)
inventory_df$path <- normalizePath(inventory_df$path, winslash = "/", mustWork = FALSE)
inventory_df$basename <- basename(inventory_df$path)
inventory_df$ext <- tools::file_ext(inventory_df$basename)

inventory_path <- file.path(output_dir, "Project_Manifest_Files.csv")
write.csv(inventory_df, inventory_path, row.names = FALSE)

cat("Manifest written:\n  ", manifest_path, "\n", sep = "")
cat("Inventory written:\n  ", inventory_path, "\n\n", sep = "")
# ============================================================
# QMD/HTML REPORT GENERATION BLOCK (paste-at-end compliant)
# ============================================================
# Design goals:
#   - Write BOTH the .qmd and .html into output_dir
#   - Avoid path/working-directory rendering failures
#   - Minimize redundancy
#   - Avoid using R reserved keywords / problematic shadowing names
# ============================================================

# ---- Guardrails: required objects must exist (created earlier in your script)
required_objects <- c(
  "script_name", "script_full", "input_path", "output_dir",
  "timestamp", "gene_col", "func_col", "cutoff_n", "days"
)
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects) > 0) {
  stop(
    "Report block cannot run because these required objects are missing: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# ---- Absolute paths (capture once)
abs_script_path  <- normalizePath(as.character(script_full), winslash = "/", mustWork = FALSE)
abs_input_file   <- normalizePath(as.character(input_path),  winslash = "/", mustWork = TRUE)
abs_output_path  <- normalizePath(as.character(output_dir),  winslash = "/", mustWork = TRUE)

# ---- Capture commented header (best-effort; avoid huge capture)
header_capture_text <- "# (Header capture skipped: script file not found.)"
if (file.exists(abs_script_path)) {
  script_text_vec <- readLines(abs_script_path, warn = FALSE)
  comment_only_vec <- script_text_vec[grepl("^\\s*#", script_text_vec)]
  if (length(comment_only_vec) > 0) {
    header_capture_text <- paste(utils::head(comment_only_vec, 250), collapse = "\n")
  } else {
    header_capture_text <- "# (No leading comment lines detected.)"
  }
}

# ---- Filenames inside output_dir
qmd_file_name  <- paste0(script_name, "_Report.qmd")
html_file_name <- paste0(tools::file_path_sans_ext(qmd_file_name), ".html")

qmd_full_path  <- file.path(abs_output_path, qmd_file_name)
html_full_path <- file.path(abs_output_path, html_file_name)

# ---- Dependencies list (documentary only)
dependency_pkg_vec <- c("utils", "tools", "stats", "ggplot2", "jsonlite", "quarto", "knitr")

# ---- Inventory outputs (relative paths from output_dir)
output_inventory_rel <- list.files(abs_output_path, recursive = TRUE, full.names = FALSE)

# ---- Build the QMD text
math_section_vec <- c(
  "",
  "### Scoring Logic (Formulas)",
  "",
  "Gene-wise Z-score across all samples:",
  "",
  "$$ Z_{g,s} = \\\\frac{X_{g,s} - \\\\mu_g}{\\\\sigma_g} $$",
  "",
  "Per-gene per-day summaries (across replicates within a day):",
  "",
  "$$ \\\\overline{Z}_{g,d} = \\\\frac{1}{R_d}\\\\sum_{r \\\\in d} Z_{g,r} \\\\quad\\\\quad \\\\overline{|Z|}_{g,d} = \\\\frac{1}{R_d}\\\\sum_{r \\\\in d} |Z_{g,r}| $$",
  "",
  "Function-level aggregation (unweighted mean across genes in a group):",
  "",
  "$$ \\\\text{mean\\\\_z}_{f,d} = \\\\frac{1}{N_f}\\\\sum_{g \\\\in f} \\\\overline{Z}_{g,d} \\\\quad\\\\quad \\\\text{mean\\\\_abs\\\\_z}_{f,d} = \\\\frac{1}{N_f}\\\\sum_{g \\\\in f} \\\\overline{|Z|}_{g,d} $$",
  ""
)

report_text_vec <- c(
  "---",
  paste0("title: \"", script_name, " Report\""),
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 3",
  "    theme: flatly",
  "    number-sections: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "## 1. Summary",
  "",
  paste0(
    "This report documents the outputs from **", script_name, "**. ",
    "The script performs gene-wise Z-scoring across all expression columns, computes gene-level day means ",
    "(signed Z and |Z|), aggregates values within a user-selected functional grouping column, writes wide and ",
    "tidy/long tables, generates trajectory plots (cutoff-filtered), and records a manifest."
  ),
  "",
  "### Input data format and structure",
  "",
  "- A CSV table containing a gene identifier column (selected interactively).",
  "- A functional grouping column (selected interactively).",
  "- Expression columns named using `DAY_REP` convention (e.g., `14_3`).",
  "- Days are inferred from the substring before the first underscore.",
  "",
  "### Interactivity / user prompts",
  "",
  "- Choose input CSV (Mac-style chooser when available).",
  "- Select gene symbol column and functional grouping column.",
  "- Provide cutoff N (minimum unique gene symbols per functional group).",
  "",
  "## 2. Script Header (captured comments)",
  "",
  "```text",
  header_capture_text,
  "```",
  "",
  "## 3. Metadata",
  "",
  paste0("- **Script name:** `", script_name, "`"),
  paste0("- **Detected script file:** `", abs_script_path, "`"),
  paste0("- **Input file:** `", abs_input_file, "`"),
  paste0("- **Output directory:** `", abs_output_path, "`"),
  paste0("- **Timestamp:** `", timestamp, "`"),
  paste0("- **Gene symbol column:** `", gene_col, "`"),
  paste0("- **Functional grouping column:** `", func_col, "`"),
  paste0("- **Cutoff N (unique genes):** `", cutoff_n, "`"),
  paste0("- **Days detected:** `", paste(days, collapse = ", "), "`"),
  "",
  "## 4. Dependencies",
  "",
  paste0("- ", paste(dependency_pkg_vec, collapse = "\n- ")),
  ""
)

report_text_vec <- c(report_text_vec, math_section_vec)

report_text_vec <- c(
  report_text_vec,
  "## 5. Generated Output Files",
  "",
  "```text",
  if (length(output_inventory_rel) > 0) paste(output_inventory_rel, collapse = "\n") else "(No files detected in output directory.)",
  "```",
  "",
  "### Key output tables (brief descriptions)",
  "",
  "- `Z-scores_all.csv`: Wide function table for all groups.",
  "- `Z-scores_cutoff.csv`: Wide function table filtered to groups with `n_genes >= N`.",
  "- `Z-scores_function_long.csv`: Tidy/long function table (function × day).",
  "- `Z-scores_function_long_cutoff.csv`: Tidy/long function table filtered by cutoff.",
  "- `Z-scores_gene_wide.csv`: Gene-level wide table with `GeneMeanZ_DayX` and `GeneMeanAbsZ_DayX`.",
  "- `Z-scores_gene_long.csv`: Gene-level tidy/long table (gene × day).",
  "- `Project_Manifest.json`: Provenance manifest.",
  "- `Project_Manifest_Files.csv`: Inventory of outputs.",
  "",
  "## 6. Figures / Plots",
  "",
  "The following plots (PNG/PDF) found in the output directory are embedded below, if any exist.",
  "",
  "```{r}",
  "plot_path_vec <- list.files('.', pattern = '\\\\.(png|pdf)$', recursive = TRUE, full.names = TRUE)",
  "if (length(plot_path_vec) == 0) {",
  "  cat('No plot files found.')",
  "} else {",
  "  for (plot_item in plot_path_vec) {",
  "    cat('\\n### ', plot_item, '\\n\\n', sep = '')",
  "    knitr::include_graphics(plot_item)",
  "    cat('\\n\\n')",
  "  }",
  "}",
  "```",
  "",
  "## 7. Interpretation of Results",
  "",
  "- **mean_z (signed):** Directional coordination of the functional group on a given day relative to each gene’s baseline.",
  "- **mean_abs_z (magnitude):** Strength/engagement of the functional group on a given day irrespective of direction.",
  "- **n_genes:** Number of **unique gene symbols** contributing to the functional group score.",
  "",
  "## 8. Session Info",
  "",
  "```{r}",
  "sessionInfo()",
  "```",
  ""
)

# ---- Write + Render: force Quarto to render in output_dir
prev_working_dir <- getwd()
on.exit(setwd(prev_working_dir), add = TRUE)
setwd(abs_output_path)

writeLines(report_text_vec, qmd_file_name)
cat("\nQuarto .qmd report written to:\n  ", qmd_full_path, "\n", sep = "")

cat("\nRendering HTML report...\n")
render_ok <- TRUE
render_msg <- NULL

tryCatch(
  {
    # Render in the current working directory (output_dir), and name the output explicitly.
    quarto::quarto_render(
      input         = qmd_file_name,
      output_format = "html",
      output_file   = html_file_name
    )
  },
  error = function(err_obj) {
    render_ok <<- FALSE
    render_msg <<- conditionMessage(err_obj)
  }
)

if (!render_ok) {
  cat("\nFAILED: Quarto render error:\n  ", render_msg, "\n", sep = "")
} else {
  cat("\nSUCCESS: HTML report saved to:\n  ", html_full_path, "\n", sep = "")
  if (interactive() && file.exists(html_full_path)) {
    browseURL(html_full_path)
  }
}

cat("\nDone.\n")
