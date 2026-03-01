#!/usr/bin/env Rscript
# ======================================================================
# Sequential Jaccard + Directional Overlap Decomposition (Up/Down + sign-switch)
#
# PURPOSE
#   Given a CSV where:
#     - Column 1 is Day1
#     - Column 2 is Day2
#     - Columns 3..N are gene-list columns (semicolon/comma-delimited strings)
#   This script:
#     1) Computes sequential Jaccard similarity for each gene-list column (long + wide).
#     2) Computes sequential Jaccard for a Combined Set per transition (union across gene-list columns).
#     3) If the input contains both an "Up" and "Down" gene-list column (auto-detected),
#        it additionally computes overlap categories across sequential transitions:
#           - Up∩Up (same-direction Up recurrence)
#           - Down∩Down (same-direction Down recurrence)
#           - Up∩Down and Down∩Up (sign-switch recurrence)
#        and writes gene lists for each overlap category.
#
# OUTPUTS (under outputs/<run_id>/)
#   - Jaccard_Long_ByColumn.csv
#   - Jaccard_Wide_ByColumn.csv
#   - Jaccard_CombinedSet.csv
#   - Jaccard_CombinedSet.png
#   - Jaccard_ByColumn_Heatmap.png
#   - Directional_Overlap_Decomposition.csv              (if Up/Down detected)
#   - Directional_Overlap_Genes__UpUp.csv                (if Up/Down detected)
#   - Directional_Overlap_Genes__DownDown.csv            (if Up/Down detected)
#   - Directional_Overlap_Genes__UpDown.csv              (if Up/Down detected)
#   - Directional_Overlap_Genes__DownUp.csv              (if Up/Down detected)
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - <script_name>_Report.qmd + rendered HTML
# ======================================================================

# -----------------------------
# Internal feature flags (MANDATORY)
# -----------------------------
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

# -----------------------------
# Runtime tracking (MANDATORY if enabled)
# -----------------------------
start_time <- Sys.time()

# -----------------------------
# Output policy (MANDATORY): all outputs under outputs/<run_id>/
# -----------------------------
analysis_name <- "Jaccard_AnyGeneListColumns_withDirectionalOverlaps"
source_label  <- "CSV"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0(analysis_name, "_", source_label, "_", timestamp)

outputs_root <- normalizePath(file.path(getwd(), "outputs"), winslash = "/", mustWork = FALSE)
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

output_dir <- normalizePath(file.path(outputs_root, run_id), winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Script identity capture (MANDATORY block)
# -----------------------------
resolve_script_path <- function() {
  p <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      rstudioapi::getSourceEditorContext()$path
    } else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  
  NA_character_
}

known_script_filename <- "Jaccard_DirectionalOverlap_Decomposition.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()
if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

cat("\n==================== SCRIPT IDENTITY ====================\n")
cat("script_name: ", script_name, "\n", sep = "")
cat("script_path: ", script_path, "\n", sep = "")
cat("script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
if (is.na(script_full)) cat("NOTE: script_full path detection failed; using fallback known_script_filename.\n")
cat("run_id:      ", run_id, "\n", sep = "")
cat("output_dir:  ", output_dir, "\n", sep = "")
cat("========================================================\n\n")

# -----------------------------
# Dependency handling (MANDATORY)
# -----------------------------
req_pkgs <- c("jsonlite", "knitr", "rmarkdown", "ggplot2")
if (enable_plotly) req_pkgs <- c(req_pkgs, "plotly")

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}
install_if_missing(req_pkgs)

suppressPackageStartupMessages({
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
  library(ggplot2)
  if (enable_plotly) library(plotly)
})

get_pkg_versions <- function(pkgs) {
  do.call(rbind, lapply(pkgs, function(p) {
    data.frame(package = p, version = as.character(utils::packageVersion(p)), stringsAsFactors = FALSE)
  }))
}
dependencies_df <- get_pkg_versions(req_pkgs)

# -----------------------------
# Input selection + STRICT FORMAT VALIDATION (MANDATORY)
# -----------------------------
pick_file <- function(prompt = "Select input CSV") {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::selectFile(caption = prompt), error = function(e) "")
    if (nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  cat(prompt, "\n")
  cat("Enter path to CSV file: ")
  p <- trimws(readLines(con = stdin(), n = 1))
  if (!nzchar(p)) stop("No input file provided.")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

validate_jaccard_wide <- function(dat) {
  # 1) Minimum columns
  if (ncol(dat) < 3) {
    return(list(ok = FALSE, msg = "File has <3 columns. Expected: Day1, Day2, and ≥1 gene-list column."))
  }
  
  # 2) Require Day columns by NAME *or* accept by position but warn
  cn <- colnames(dat)
  has_named <- all(c("Day1","Day2") %in% cn)
  
  # If named, use those. If not named, still allow but REQUIRE numeric-like content.
  if (!has_named) {
    # Position-based fallback validation
    d1 <- suppressWarnings(as.numeric(dat[[1]]))
    d2 <- suppressWarnings(as.numeric(dat[[2]]))
    if (any(is.na(d1)) || any(is.na(d2))) {
      return(list(
        ok = FALSE,
        msg = paste0(
          "First two columns are not named Day1/Day2 and are not clean numeric days.\n",
          "Fix by naming columns exactly 'Day1' and 'Day2', OR ensure column 1/2 are numeric day codes.\n",
          "Column 1 name: '", cn[1], "'; Column 2 name: '", cn[2], "'."
        )
      ))
    }
  }
  
  # 3) Must have at least one gene-list column with delimiter-like content in some rows
  gene_cols <- cn[3:ncol(dat)]
  nonempty <- sapply(gene_cols, function(x) {
    v <- dat[[x]]
    any(!is.na(v) & nzchar(trimws(as.character(v))))
  })
  if (!any(nonempty)) {
    return(list(
      ok = FALSE,
      msg = paste0(
        "None of columns 3..N contain any non-empty gene-list strings.\n",
        "Expected semicolon/comma-delimited gene lists in at least one column after Day1/Day2."
      )
    ))
  }
  
  # 4) Strong reject: looks like an expression matrix (many numeric columns, no delimiters)
  #    (common wrong pick: wide expression for DTW)
  looks_expr <- {
    # if >50% of columns 3..N are numeric and almost none contain ';' or ','
    is_num <- sapply(gene_cols, function(x) is.numeric(dat[[x]]) || all(grepl("^\\s*[-0-9.]+\\s*$", na.omit(as.character(dat[[x]])))))
    has_delim <- sapply(gene_cols, function(x) any(grepl("[;,|]", na.omit(as.character(dat[[x]])))))
    mean(is_num) > 0.5 && mean(has_delim) < 0.1
  }
  if (looks_expr) {
    return(list(
      ok = FALSE,
      msg = paste0(
        "This file looks like a WIDE EXPRESSION MATRIX (numeric columns), not a Jaccard gene-list file.\n",
        "This script requires gene lists like 'STAT1; IRF7; ...' in columns 3..N."
      )
    ))
  }
  
  list(ok = TRUE, msg = "OK")
}

args <- commandArgs(trailingOnly = TRUE)

prompt_txt <- paste0(
  "Select the **Jaccard-input WIDE** CSV produced by your PCA script:\n",
  "  Required layout:\n",
  "    - Column 1 = Day1 (or named 'Day1')\n",
  "    - Column 2 = Day2 (or named 'Day2')\n",
  "    - Columns 3..N = gene-list columns (strings like 'GENE1; GENE2; ...')\n",
  "  NOT an expression matrix.\n"
)

if (length(args) >= 1 && nzchar(args[1])) {
  input_csv <- normalizePath(args[1], winslash = "/", mustWork = FALSE)
} else {
  input_csv <- pick_file(prompt_txt)
}

if (!file.exists(input_csv)) stop("Input CSV does not exist: ", input_csv)

dat <- utils::read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

v <- validate_jaccard_wide(dat)
if (!isTRUE(v$ok)) {
  stop(
    "❌ Wrong input file format.\n\n",
    v$msg, "\n\n",
    "You selected: ", input_csv, "\n",
    "Columns found: ", paste(colnames(dat), collapse = ", "), "\n"
  )
}

cat("✅ Input format validated as Jaccard gene-list WIDE file.\n\n")

# Enforce by position
day1_colname <- colnames(dat)[1]
day2_colname <- colnames(dat)[2]
gene_cols <- colnames(dat)[3:ncol(dat)]

cat("Assuming Day1 = column 1: ", day1_colname, "\n", sep = "")
cat("Assuming Day2 = column 2: ", day2_colname, "\n", sep = "")
cat("Gene-list columns (3..N): ", paste(gene_cols, collapse = ", "), "\n\n", sep = "")

# -----------------------------
# Helpers
# -----------------------------
guess_split_regex <- function(x) {
  if (grepl(";", x, fixed = TRUE)) return("\\s*;\\s*")
  if (grepl(",", x, fixed = TRUE)) return("\\s*,\\s*")
  if (grepl("\\|", x)) return("\\s*\\|\\s*")
  "\\s+"
}

parse_gene_list <- function(x) {
  if (is.na(x)) return(character(0))
  x <- trimws(as.character(x))
  if (!nzchar(x)) return(character(0))
  rx <- guess_split_regex(x)
  g <- unlist(strsplit(x, rx, perl = TRUE), use.names = FALSE)
  g <- trimws(g)
  g <- g[nzchar(g)]
  unique(g)
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

collapse_genes <- function(x) {
  x <- unique(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return("")
  paste(sort(x), collapse = "; ")
}

# -----------------------------
# Sort transitions
# -----------------------------
d1n <- suppressWarnings(as.numeric(dat[[1]]))
d2n <- suppressWarnings(as.numeric(dat[[2]]))
if (all(!is.na(d1n)) && all(!is.na(d2n))) {
  dat <- dat[order(d1n, d2n), , drop = FALSE]
} else {
  dat <- dat[order(as.character(dat[[1]]), as.character(dat[[2]])), , drop = FALSE]
}

dat$transition <- paste0(dat[[1]], "-", dat[[2]])
if (nrow(dat) < 2) stop("Need at least two rows (transitions) to compute sequential overlaps.")

# Sanity metric: fraction of non-empty cells per gene-list column
nonempty_fraction <- sapply(gene_cols, function(cn) {
  x <- dat[[cn]]
  mean(!is.na(x) & nzchar(trimws(as.character(x))))
})

# -----------------------------
# Pre-parse sets for each gene-list column
# -----------------------------
sets_by_col <- lapply(gene_cols, function(cn) lapply(dat[[cn]], parse_gene_list))
names(sets_by_col) <- gene_cols

# Combined set per transition
combined_set <- vector("list", nrow(dat))
for (i in seq_len(nrow(dat))) {
  s <- character(0)
  for (cn in gene_cols) s <- union(s, sets_by_col[[cn]][[i]])
  combined_set[[i]] <- s
}

# -----------------------------
# Sequential Jaccard per gene-list column (long)
# -----------------------------
long_rows <- list()
k <- 1
for (i in 1:(nrow(dat) - 1)) {
  t1 <- dat$transition[i]
  t2 <- dat$transition[i + 1]
  for (cn in gene_cols) {
    A <- sets_by_col[[cn]][[i]]
    B <- sets_by_col[[cn]][[i + 1]]
    long_rows[[k]] <- data.frame(
      column = cn,
      transition_1 = t1,
      transition_2 = t2,
      pair = paste0(t1, " → ", t2),
      jaccard = jaccard(A, B),
      shared_n = length(intersect(A, B)),
      union_n  = length(union(A, B)),
      size_A   = length(A),
      size_B   = length(B),
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }
}
jaccard_long <- do.call(rbind, long_rows)

# Wide (pair order forced)
pair_order <- paste0(dat$transition[1:(nrow(dat)-1)], " → ", dat$transition[2:nrow(dat)])
jaccard_wide <- data.frame(pair = pair_order, stringsAsFactors = FALSE)
for (cn in gene_cols) {
  tmp <- jaccard_long[jaccard_long$column == cn, ]
  jaccard_wide[[cn]] <- tmp$jaccard[match(jaccard_wide$pair, tmp$pair)]
}

# Combined set sequential Jaccard
comb_rows <- list()
for (i in 1:(nrow(dat) - 1)) {
  t1 <- dat$transition[i]
  t2 <- dat$transition[i + 1]
  A <- combined_set[[i]]
  B <- combined_set[[i + 1]]
  comb_rows[[i]] <- data.frame(
    transition_1 = t1,
    transition_2 = t2,
    pair = paste0(t1, " → ", t2),
    jaccard_combined = jaccard(A, B),
    shared_n = length(intersect(A, B)),
    union_n  = length(union(A, B)),
    size_A   = length(A),
    size_B   = length(B),
    stringsAsFactors = FALSE
  )
}
jaccard_combined <- do.call(rbind, comb_rows)

# -----------------------------
# Directional overlap decomposition (Up/Down + sign-switch)
# Auto-detect Up/Down columns among gene_cols.
#   - "Up" column: name contains "up" and NOT "down"
#   - "Down" column: name contains "down"
# If multiple candidates exist, choose the first by appearance.
# -----------------------------
to_lower <- function(x) tolower(gsub("\\s+", "", x))

up_candidates <- gene_cols[grepl("up", to_lower(gene_cols)) & !grepl("down", to_lower(gene_cols))]
down_candidates <- gene_cols[grepl("down", to_lower(gene_cols))]

up_col <- if (length(up_candidates) >= 1) up_candidates[1] else NA_character_
down_col <- if (length(down_candidates) >= 1) down_candidates[1] else NA_character_

directional_ok <- !is.na(up_col) && !is.na(down_col)

if (directional_ok) {
  cat("Directional overlap detected.\n")
  cat("Using Up column:   ", up_col, "\n", sep = "")
  cat("Using Down column: ", down_col, "\n\n", sep = "")
} else {
  cat("Directional overlap NOT computed (could not auto-detect both Up and Down columns).\n")
  if (length(up_candidates) == 0) cat(" - No Up-like column found.\n")
  if (length(down_candidates) == 0) cat(" - No Down-like column found.\n")
  cat("\n")
}

directional_decomp <- NULL
genes_upup <- NULL
genes_downdown <- NULL
genes_updown <- NULL
genes_downup <- NULL

if (directional_ok) {
  rows <- list()
  g_upup <- list()
  g_dd <- list()
  g_ud <- list()
  g_du <- list()
  
  for (i in 1:(nrow(dat) - 1)) {
    t1 <- dat$transition[i]
    t2 <- dat$transition[i + 1]
    pair <- paste0(t1, " → ", t2)
    
    Up1 <- sets_by_col[[up_col]][[i]]
    Down1 <- sets_by_col[[down_col]][[i]]
    Up2 <- sets_by_col[[up_col]][[i + 1]]
    Down2 <- sets_by_col[[down_col]][[i + 1]]
    
    UpUp <- intersect(Up1, Up2)
    DownDown <- intersect(Down1, Down2)
    UpDown <- intersect(Up1, Down2)
    DownUp <- intersect(Down1, Up2)
    
    rows[[i]] <- data.frame(
      transition_1 = t1,
      transition_2 = t2,
      pair = pair,
      
      upup_n = length(UpUp),
      downdown_n = length(DownDown),
      updown_n = length(UpDown),
      downup_n = length(DownUp),
      
      upup_genes = collapse_genes(UpUp),
      downdown_genes = collapse_genes(DownDown),
      updown_genes = collapse_genes(UpDown),
      downup_genes = collapse_genes(DownUp),
      
      stringsAsFactors = FALSE
    )
    
    # Also store long-form gene lists for separate files (one gene per row)
    if (length(UpUp) > 0) g_upup[[length(g_upup)+1]] <- data.frame(pair = pair, gene = sort(UpUp), stringsAsFactors = FALSE)
    if (length(DownDown) > 0) g_dd[[length(g_dd)+1]] <- data.frame(pair = pair, gene = sort(DownDown), stringsAsFactors = FALSE)
    if (length(UpDown) > 0) g_ud[[length(g_ud)+1]] <- data.frame(pair = pair, gene = sort(UpDown), stringsAsFactors = FALSE)
    if (length(DownUp) > 0) g_du[[length(g_du)+1]] <- data.frame(pair = pair, gene = sort(DownUp), stringsAsFactors = FALSE)
  }
  
  directional_decomp <- do.call(rbind, rows)
  genes_upup <- if (length(g_upup) > 0) do.call(rbind, g_upup) else data.frame(pair = character(0), gene = character(0), stringsAsFactors = FALSE)
  genes_downdown <- if (length(g_dd) > 0) do.call(rbind, g_dd) else data.frame(pair = character(0), gene = character(0), stringsAsFactors = FALSE)
  genes_updown <- if (length(g_ud) > 0) do.call(rbind, g_ud) else data.frame(pair = character(0), gene = character(0), stringsAsFactors = FALSE)
  genes_downup <- if (length(g_du) > 0) do.call(rbind, g_du) else data.frame(pair = character(0), gene = character(0), stringsAsFactors = FALSE)
}

# -----------------------------
# Write outputs
# -----------------------------
csv_long <- file.path(output_dir, "Jaccard_Long_ByColumn.csv")
csv_wide <- file.path(output_dir, "Jaccard_Wide_ByColumn.csv")
csv_comb <- file.path(output_dir, "Jaccard_CombinedSet.csv")

utils::write.csv(jaccard_long, csv_long, row.names = FALSE)
utils::write.csv(jaccard_wide, csv_wide, row.names = FALSE)
utils::write.csv(jaccard_combined, csv_comb, row.names = FALSE)

# Directional outputs (if available)
csv_dir_decomp <- NA_character_
csv_upup_genes <- NA_character_
csv_dd_genes <- NA_character_
csv_ud_genes <- NA_character_
csv_du_genes <- NA_character_

if (directional_ok) {
  csv_dir_decomp <- file.path(output_dir, "Directional_Overlap_Decomposition.csv")
  csv_upup_genes <- file.path(output_dir, "Directional_Overlap_Genes__UpUp.csv")
  csv_dd_genes <- file.path(output_dir, "Directional_Overlap_Genes__DownDown.csv")
  csv_ud_genes <- file.path(output_dir, "Directional_Overlap_Genes__UpDown.csv")
  csv_du_genes <- file.path(output_dir, "Directional_Overlap_Genes__DownUp.csv")
  
  utils::write.csv(directional_decomp, csv_dir_decomp, row.names = FALSE)
  utils::write.csv(genes_upup, csv_upup_genes, row.names = FALSE)
  utils::write.csv(genes_downdown, csv_dd_genes, row.names = FALSE)
  utils::write.csv(genes_updown, csv_ud_genes, row.names = FALSE)
  utils::write.csv(genes_downup, csv_du_genes, row.names = FALSE)
}

# -----------------------------
# Plots (static always-on)
# -----------------------------
p_comb <- ggplot(jaccard_combined, aes(x = pair, y = jaccard_combined)) +
  geom_col() +
  labs(
    title = "Sequential Jaccard Similarity (Combined Gene Set Across All Gene-List Columns)",
    x = "Sequential transition pairs",
    y = "Jaccard"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png_comb <- file.path(output_dir, "Jaccard_CombinedSet.png")
ggsave(png_comb, p_comb, width = 11, height = 6, dpi = 220)

heat_df <- do.call(rbind, lapply(gene_cols, function(cn) {
  data.frame(pair = jaccard_wide$pair, column = cn, jaccard = jaccard_wide[[cn]], stringsAsFactors = FALSE)
}))

p_heat <- ggplot(heat_df, aes(x = pair, y = column, fill = jaccard)) +
  geom_tile() +
  labs(
    title = "Sequential Jaccard by Column (each gene-list column treated independently)",
    x = "Sequential transition pairs",
    y = "Gene-list column"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png_heat <- file.path(output_dir, "Jaccard_ByColumn_Heatmap.png")
ggsave(png_heat, p_heat, width = 12, height = max(4, 0.35 * length(gene_cols) + 2), dpi = 220)

# If directional overlap exists, make one simple plot of counts per pair
png_dir_counts <- NA_character_
if (directional_ok) {
  dfc <- directional_decomp
  dfc_long <- do.call(rbind, list(
    data.frame(pair = dfc$pair, category = "Up∩Up", n = dfc$upup_n, stringsAsFactors = FALSE),
    data.frame(pair = dfc$pair, category = "Down∩Down", n = dfc$downdown_n, stringsAsFactors = FALSE),
    data.frame(pair = dfc$pair, category = "Up∩Down", n = dfc$updown_n, stringsAsFactors = FALSE),
    data.frame(pair = dfc$pair, category = "Down∩Up", n = dfc$downup_n, stringsAsFactors = FALSE)
  ))
  
  p_dir <- ggplot(dfc_long, aes(x = pair, y = n)) +
    geom_col() +
    facet_wrap(~category, ncol = 2, scales = "free_y") +
    labs(
      title = "Directional overlap counts across sequential transition pairs",
      x = "Sequential transition pairs",
      y = "Number of recurring genes"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png_dir_counts <- file.path(output_dir, "Directional_Overlap_Counts.png")
  ggsave(png_dir_counts, p_dir, width = 12, height = 7, dpi = 220)
}

# -----------------------------
# Manifest + inventory (MANDATORY)
# -----------------------------
manifest_json <- file.path(output_dir, "Project_Manifest.json")
manifest_files_csv <- file.path(output_dir, "Project_Manifest_Files.csv")

list_files_rel <- function(dir) {
  files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    return(data.frame(file = character(0), full_path = character(0), size_bytes = numeric(0),
                      modified_time = character(0), stringsAsFactors = FALSE))
  }
  data.frame(
    file = gsub(paste0("^", gsub("([\\^\\$\\.|\\*\\+\\?\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", dir), "/?"), "", files),
    full_path = normalizePath(files, winslash = "/", mustWork = FALSE),
    size_bytes = file.info(files)$size,
    modified_time = as.character(file.info(files)$mtime),
    stringsAsFactors = FALSE
  )
}

runtime_seconds <- NA_real_
if (enable_runtime_tracking) runtime_seconds <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)

deps_list <- lapply(seq_len(nrow(dependencies_df)), function(i) {
  list(package = dependencies_df$package[i], version = dependencies_df$version[i])
})

manifest <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script = list(name = script_name, path = script_path, full_path = ifelse(is.na(script_full), NA_character_, script_full)),
  input = list(input_csv = normalizePath(input_csv, winslash = "/", mustWork = FALSE), n_rows = nrow(dat), n_cols = ncol(dat)),
  parameters = list(
    analysis_name = analysis_name,
    source_label = source_label,
    enable_plotly = enable_plotly,
    enable_runtime_tracking = enable_runtime_tracking,
    runtime_seconds = runtime_seconds,
    day_columns_by_position = list(col1 = day1_colname, col2 = day2_colname),
    gene_list_columns = gene_cols,
    gene_list_nonempty_fraction = as.list(nonempty_fraction),
    parsing_rule = "Split by detected delimiter per cell (prefers ';', then ',', then '|', then whitespace).",
    comparisons = "Sequential transitions in sorted row order (numeric Day1/Day2 if possible).",
    directional_overlap = list(
      computed = directional_ok,
      up_column = ifelse(directional_ok, up_col, NA_character_),
      down_column = ifelse(directional_ok, down_col, NA_character_)
    )
    
  ),
  dependencies = deps_list,
  outputs = list(outputs_root = outputs_root, output_dir = output_dir)
)

writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"), con = manifest_json)
utils::write.csv(list_files_rel(output_dir), manifest_files_csv, row.names = FALSE)

# -----------------------------
# Build Quarto QMD + render to HTML (MANDATORY final step)
# -----------------------------
get_header_block <- function() {
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
    take_n <- min(length(lines), 220)
    head_lines <- lines[seq_len(take_n)]
    keep <- c()
    started <- FALSE
    for (ln in head_lines) {
      if (!started) {
        if (grepl("^\\s*#|^\\s*$", ln)) { keep <- c(keep, ln); started <- TRUE } else break
      } else {
        if (grepl("^\\s*#|^\\s*$", ln)) keep <- c(keep, ln) else break
      }
    }
    paste(keep, collapse = "\n")
  } else {
    paste0("# Header block unavailable (script_full not detected). Fallback: ", known_script_filename)
  }
}

header_block <- get_header_block()

# Pre-normalize paths as plain strings (prevents winslash quoting bugs)
norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

input_csv_norm <- norm(input_csv)
output_dir_norm <- norm(output_dir)
csv_long_norm <- norm(csv_long)
csv_wide_norm <- norm(csv_wide)
csv_comb_norm <- norm(csv_comb)
png_comb_norm <- norm(png_comb)
png_heat_norm <- norm(png_heat)
manifest_files_csv_norm <- norm(manifest_files_csv)
png_dir_counts_norm <- if (!is.na(png_dir_counts)) norm(png_dir_counts) else ""
csv_dir_decomp_norm <- if (!is.na(csv_dir_decomp)) norm(csv_dir_decomp) else ""
csv_upup_genes_norm <- if (!is.na(csv_upup_genes)) norm(csv_upup_genes) else ""
csv_dd_genes_norm <- if (!is.na(csv_dd_genes)) norm(csv_dd_genes) else ""
csv_ud_genes_norm <- if (!is.na(csv_ud_genes)) norm(csv_ud_genes) else ""
csv_du_genes_norm <- if (!is.na(csv_du_genes)) norm(csv_du_genes) else ""

# YAML-safe quoting helper
q <- function(x) {
  x <- gsub("\\\\", "/", x)
  x <- gsub("\"", "\\\\\"", x)
  paste0("\"", x, "\"")
}

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_out <- file.path(output_dir, paste0(script_name, "_Report.html"))

qmd <- paste0(
  "---\n",
  "title: \"Sequential Jaccard + Directional Overlap Decomposition\"\n",
  "format:\n",
  "  html:\n",
  "    toc: true\n",
  "execute:\n",
  "  echo: false\n",
  "  warning: false\n",
  "  message: false\n",
  "params:\n",
  "  run_id: ", q(run_id), "\n",
  "  output_dir: ", q(output_dir_norm), "\n",
  "  input_csv: ", q(input_csv_norm), "\n",
  "  csv_long: ", q(csv_long_norm), "\n",
  "  csv_wide: ", q(csv_wide_norm), "\n",
  "  csv_comb: ", q(csv_comb_norm), "\n",
  "  png_comb: ", q(png_comb_norm), "\n",
  "  png_heat: ", q(png_heat_norm), "\n",
  "  manifest_files_csv: ", q(manifest_files_csv_norm), "\n",
  "  enable_plotly: ", ifelse(enable_plotly, "true", "false"), "\n",
  "  runtime_seconds: ", ifelse(is.na(runtime_seconds), "null", as.character(runtime_seconds)), "\n",
  "  directional_ok: ", ifelse(directional_ok, "true", "false"), "\n",
  "  dir_decomp_csv: ", q(csv_dir_decomp_norm), "\n",
  "  dir_genes_upup_csv: ", q(csv_upup_genes_norm), "\n",
  "  dir_genes_downdown_csv: ", q(csv_dd_genes_norm), "\n",
  "  dir_genes_updown_csv: ", q(csv_ud_genes_norm), "\n",
  "  dir_genes_downup_csv: ", q(csv_du_genes_norm), "\n",
  "  dir_counts_png: ", q(png_dir_counts_norm), "\n",
  "---\n\n",
  
  "# Script header\n\n",
  "```\n", header_block, "\n```\n\n",
  
  "# Metadata\n\n",
  "- **run_id:** `", run_id, "`\n",
  "- **output_dir:** `", output_dir_norm, "`\n",
  "- **input_csv:** `", input_csv_norm, "`\n\n",
  
  "# Runtime\n\n",
  "Runtime (seconds): `r params$runtime_seconds`\n\n",
  
  "# Dependencies\n\n",
  "```{r}\n",
  "suppressPackageStartupMessages(library(knitr))\n",
  "deps <- ", "data.frame(", 
  "package = c(", paste(sprintf("\"%s\"", dependencies_df$package), collapse = ","), "), ",
  "version = c(", paste(sprintf("\"%s\"", dependencies_df$version), collapse = ","), "), ",
  "stringsAsFactors = FALSE)\n",
  "knitr::kable(deps)\n",
  "```\n\n",
  
  "# Analytical logic\n\n",
  "**Column conventions:** Column 1 is Day1, Column 2 is Day2, Columns 3..N are treated as gene-list columns.\n\n",
  "**Jaccard similarity:** For each gene-list column and for the combined set (union across gene-list columns), we compute sequential similarity across adjacent transition rows:\n\n",
  "\\[ J(A,B) = \\frac{|A \\cap B|}{|A \\cup B|} \\]\n\n",
  "**Directional overlap (if Up/Down columns detected):** For sequential transition pairs, we decompose recurrence into:\n\n",
  "- Up∩Up (same-direction Up recurrence)\n",
  "- Down∩Down (same-direction Down recurrence)\n",
  "- Up∩Down and Down∩Up (sign-switch recurrence)\n\n",
  
  "# Results tables\n\n",
  "## Combined-set sequential Jaccard\n\n",
  "```{r}\n",
  "knitr::kable(read.csv(params$csv_comb, stringsAsFactors = FALSE))\n",
  "```\n\n",
  "## Per-column sequential Jaccard (long)\n\n",
  "```{r}\n",
  "knitr::kable(read.csv(params$csv_long, stringsAsFactors = FALSE))\n",
  "```\n\n",
  "## Per-column sequential Jaccard (wide)\n\n",
  "```{r}\n",
  "knitr::kable(read.csv(params$csv_wide, stringsAsFactors = FALSE))\n",
  "```\n\n",
  
  "# Figures\n\n",
  "```{r}\n",
  "knitr::include_graphics(params$png_comb)\n",
  "```\n\n",
  "```{r}\n",
  "knitr::include_graphics(params$png_heat)\n",
  "```\n\n",
  
  "# Directional overlap decomposition\n\n",
  "```{r}\n",
  "if (isTRUE(params$directional_ok)) {\n",
  "  knitr::kable(read.csv(params$dir_decomp_csv, stringsAsFactors = FALSE))\n",
  "} else {\n",
  "  cat('Directional overlap tables were not generated because Up/Down columns could not be auto-detected.')\n",
  "}\n",
  "```\n\n",
  
  "## Directional overlap gene lists\n\n",
  "### Up∩Up genes\n\n",
  "```{r}\n",
  "if (isTRUE(params$directional_ok)) {\n",
  "  knitr::kable(read.csv(params$dir_genes_upup_csv, stringsAsFactors = FALSE))\n",
  "}\n",
  "```\n\n",
  "### Down∩Down genes\n\n",
  "```{r}\n",
  "if (isTRUE(params$directional_ok)) {\n",
  "  knitr::kable(read.csv(params$dir_genes_downdown_csv, stringsAsFactors = FALSE))\n",
  "}\n",
  "```\n\n",
  "### Up∩Down genes (sign-switch)\n\n",
  "```{r}\n",
  "if (isTRUE(params$directional_ok)) {\n",
  "  knitr::kable(read.csv(params$dir_genes_updown_csv, stringsAsFactors = FALSE))\n",
  "}\n",
  "```\n\n",
  "### Down∩Up genes (sign-switch)\n\n",
  "```{r}\n",
  "if (isTRUE(params$directional_ok)) {\n",
  "  knitr::kable(read.csv(params$dir_genes_downup_csv, stringsAsFactors = FALSE))\n",
  "}\n",
  "```\n\n",
  
  "```{r}\n",
  "if (isTRUE(params$directional_ok) && nzchar(params$dir_counts_png)) {\n",
  "  knitr::include_graphics(params$dir_counts_png)\n",
  "}\n",
  "```\n\n",
  
  "# Interactive plot (plotly; if enabled)\n\n",
  "```{r}\n",
  "if (isTRUE(params$enable_plotly)) {\n",
  "  suppressPackageStartupMessages(library(ggplot2))\n",
  "  suppressPackageStartupMessages(library(plotly))\n",
  "  comb <- read.csv(params$csv_comb, stringsAsFactors = FALSE)\n",
  "  gp <- ggplot(comb, aes(x = pair, y = jaccard_combined)) +\n",
  "    geom_col() +\n",
  "    labs(title = \"Interactive: Combined-set sequential Jaccard\", x = \"Sequential transition pairs\", y = \"Jaccard\") +\n",
  "    theme_bw() +\n",
  "    theme(axis.text.x = element_text(angle = 45, hjust = 1))\n",
  "  tryCatch(plotly::ggplotly(gp), error = function(e) gp)\n",
  "}\n",
  "```\n\n",
  
  "# Interpretation\n\n",
  "- Higher combined-set Jaccard indicates reuse of the same genes as transition drivers across adjacent transitions.\n",
  "- If Up∩Up and Down∩Down are low but sign-switch overlaps (Up∩Down / Down∩Up) are high, the same genes recur but change direction.\n",
  "- If all overlap categories are low, the transition is driven by a new gene set (turnover).\n\n",
  
  "# Reproducibility\n\n",
  "```{r}\n",
  "sessionInfo()\n",
  "```\n\n",
  
  "# Generated files inventory\n\n",
  "```{r}\n",
  "knitr::kable(read.csv(params$manifest_files_csv, stringsAsFactors = FALSE))\n",
  "```\n"
)

writeLines(qmd, con = qmd_path)

# Render (final step)
render_ok <- TRUE
render_msg <- NULL

if (enable_runtime_tracking) {
  runtime_seconds <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  manifest$parameters$runtime_seconds <- runtime_seconds
  writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"), con = manifest_json)
}

tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_file = basename(html_out),
    output_dir = output_dir,
    params = list(runtime_seconds = runtime_seconds, enable_plotly = enable_plotly),
    quiet = TRUE
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# Refresh inventory after render (MANDATORY)
utils::write.csv(list_files_rel(output_dir), manifest_files_csv, row.names = FALSE)

if (render_ok && file.exists(html_out)) {
  cat("\nHTML report created: ", normalizePath(html_out, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("\nERROR: HTML report render failed.\n")
  if (!is.null(render_msg)) cat("Reason: ", render_msg, "\n", sep = "")
  cat("QMD written at: ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
}

cat("\n==================== OUTPUT SUMMARY ====================\n")
cat("Long Jaccard CSV:   ", normalizePath(csv_long, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Wide Jaccard CSV:   ", normalizePath(csv_wide, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Combined Jaccard:   ", normalizePath(csv_comb, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Plot (combined):    ", normalizePath(png_comb, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Plot (heatmap):     ", normalizePath(png_heat, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (directional_ok) {
  cat("Directional decomp: ", normalizePath(csv_dir_decomp, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("Up∩Up genes:        ", normalizePath(csv_upup_genes, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("Down∩Down genes:    ", normalizePath(csv_dd_genes, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("Up∩Down genes:      ", normalizePath(csv_ud_genes, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("Down∩Up genes:      ", normalizePath(csv_du_genes, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (!is.na(png_dir_counts)) cat("Directional plot:   ", normalizePath(png_dir_counts, winslash = "/", mustWork = FALSE), "\n", sep = "")
}
cat("Manifest JSON:      ", normalizePath(manifest_json, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Files CSV:          ", normalizePath(manifest_files_csv, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (enable_runtime_tracking) cat("Runtime (s):        ", runtime_seconds, "\n", sep = "")
cat("========================================================\n\n")
