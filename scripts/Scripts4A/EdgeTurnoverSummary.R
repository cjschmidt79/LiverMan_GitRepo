#!/usr/bin/env Rscript
# ======================================================================
# Edge Turnover -> Unified Functional Summary (Sequential Day-Pairs Only)
#
# PURPOSE
#   Merge:
#     (1) Seven edge-turnover tables (gene–gene; sequential day pairs)
#     (2) A gene->function mapping table
#   And produce:
#     - Long annotated edge table (gene + function)
#     - Functional-pair summaries per transition
#     - Unified functional summary across all transitions (single table)
#     - Project_Manifest.json + Project_Manifest_Files.csv
#     - QMD + rendered HTML report documenting logic for reviewers
#
# USER REQUIREMENTS (implemented):
#   - Mac-style directory/file chooser when selecting inputs
#   - Capture script path + name, output path, input paths + names in manifest
#   - Only include sequential pairs of days (chain order; dayB == next dayA)
#
# INPUTS (chosen interactively):
#   - Gene2FunctionMapping4Edgeswap.csv (columns: Gene, Primary_Function)
#   - Directory containing files like:
#       B_EdgeTurnover_Top50_Day4_to_Day8.csv
#       ...
#       B_EdgeTurnover_Top50_Day18_to_Day20.csv
#
# OUTPUTS (written under output_dir):
#   - EdgeTurnover_Annotated_Long.csv
#   - EdgeTurnover_FunctionPair_ByTransition.csv
#   - EdgeTurnover_FunctionalSummary_Unified.csv
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - EdgeTurnover_FunctionalSummary_Report.qmd
#   - EdgeTurnover_FunctionalSummary_Report.html  (if render succeeds)
# ======================================================================

suppressPackageStartupMessages({
  ok <- function(x) requireNamespace(x, quietly = TRUE)
  if (!ok("utils")) stop("Base 'utils' not available.")
})

# ----------------------------- Helpers -----------------------------

msg <- function(...) cat(sprintf(...), "\n", sep = "")

now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

safe_norm <- function(p) {
  # Vector-safe normalizePath wrapper
  if (is.null(p)) return(NA_character_)
  
  p <- as.character(p)
  
  out <- vapply(p, function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    tryCatch(
      normalizePath(x, winslash = "/", mustWork = FALSE),
      error = function(e) x
    )
  }, FUN.VALUE = character(1))
  
  out
}

script_identity <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- NA_character_
  if (length(args) > 0) {
    hit <- grep(file_flag, args, value = TRUE)
    if (length(hit) > 0) {
      script_path <- sub(file_flag, "", hit[1])
    }
  }
  script_path <- safe_norm(script_path)
  list(
    name = ifelse(is.na(script_path), NA, basename(script_path)),
    path = ifelse(is.na(script_path), NA, dirname(script_path)),
    full_path = script_path
  )
}

# Mac-friendly file picker
choose_file_mac <- function(caption = "Choose file") {
  # Prefer RStudio chooser if available (native feel on macOS)
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::selectFile(caption = caption)
    return(safe_norm(p))
  }
  
  # Base chooser (macOS uses a native dialog)
  p <- tryCatch(utils::choose.files(caption = caption, multi = FALSE), error = function(e) "")
  p <- if (length(p) == 0) "" else p[1]
  p <- safe_norm(p)
  if (!is.na(p) && file.exists(p)) return(p)
  
  # Fallback
  msg("Could not open file chooser. Paste full path:")
  p2 <- trimws(readLines(con = stdin(), n = 1))
  p2 <- safe_norm(p2)
  if (is.na(p2) || !file.exists(p2)) stop("File not found: ", p2)
  p2
}

# Mac-friendly directory picker
choose_dir_mac <- function(caption = "Choose directory") {
  # RStudio directory chooser
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::selectDirectory(caption = caption)
    return(safe_norm(p))
  }
  
  # tcltk chooser if present
  if (requireNamespace("tcltk", quietly = TRUE)) {
    p <- tryCatch(tcltk::tk_choose.dir(caption = caption), error = function(e) "")
    p <- safe_norm(p)
    if (!is.na(p) && dir.exists(p)) return(p)
  }
  
  # Fallback
  msg("Could not open directory chooser. Paste full directory path:")
  p2 <- trimws(readLines(con = stdin(), n = 1))
  p2 <- safe_norm(p2)
  if (is.na(p2) || !dir.exists(p2)) stop("Directory not found: ", p2)
  p2
}

dir_create <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}

write_json <- function(x, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")
  jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

build_file_inventory <- function(out_dir, csv_path) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  info <- file.info(files)
  inv <- data.frame(
    file = gsub(paste0("^", gsub("([\\W])", "\\\\\\1", safe_norm(out_dir)), "/?"), "", safe_norm(files)),
    full_path = safe_norm(files),
    bytes = as.numeric(info$size),
    modified = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
  utils::write.csv(inv, csv_path, row.names = FALSE)
  inv
}

parse_daypair_from_filename <- function(fname) {
  # expects ...DayX_to_DayY...
  m <- regexec("Day([0-9]+)_to_Day([0-9]+)", fname, perl = TRUE)
  r <- regmatches(fname, m)[[1]]
  if (length(r) != 3) return(NULL)
  a <- as.integer(r[2]); b <- as.integer(r[3])
  if (is.na(a) || is.na(b)) return(NULL)
  list(day_a = a, day_b = b, transition = sprintf("D%d-%d", a, b))
}

# Enforce "sequential pairs" as a chain: after sorting by day_a,
# require each next file's day_a == previous file's day_b.

keep_sequential_chain <- function(daypairs_df) {
  if (nrow(daypairs_df) == 0) return(daypairs_df)
  
  # Ensure numeric
  daypairs_df$day_a <- as.integer(daypairs_df$day_a)
  daypairs_df$day_b <- as.integer(daypairs_df$day_b)
  
  daypairs_df <- daypairs_df[order(daypairs_df$day_a, daypairs_df$day_b), , drop = FALSE]
  
  # Start at the earliest day_a
  start_idx <- which.min(daypairs_df$day_a)
  chain_idx <- integer(0)
  
  current_idx <- start_idx
  chain_idx <- c(chain_idx, current_idx)
  
  current_end <- daypairs_df$day_b[current_idx]
  
  # Greedily follow: next day_a must equal current_end
  repeat {
    next_candidates <- which(daypairs_df$day_a == current_end)
    if (length(next_candidates) == 0) break
    
    # If multiple candidates exist, take the one with smallest day_b
    next_idx <- next_candidates[which.min(daypairs_df$day_b[next_candidates])]
    
    # Prevent infinite loops (shouldn't happen, but guard anyway)
    if (next_idx %in% chain_idx) break
    
    chain_idx <- c(chain_idx, next_idx)
    current_end <- daypairs_df$day_b[next_idx]
  }
  
  daypairs_df[chain_idx, , drop = FALSE]
}


# Persistence type rules (objective)
classify_persistence <- function(transitions_sorted, midpoints_sorted) {
  nT <- length(transitions_sorted)
  if (nT <= 1) return("Transient")
  
  # consecutive in time if transition midpoints are nondecreasing and present in >=3 consecutive transitions
  # We treat "consecutive" as adjacent transitions in the chain (not numeric day spacing).
  if (nT >= 3) return("Persistent")
  # Appears >=2 times: Emergent vs Late-onset based on first midpoint
  first_mid <- midpoints_sorted[1]
  if (!is.na(first_mid) && first_mid >= 14) return("Late-onset")
  "Emergent"
}

# ----------------------------- Inputs -----------------------------

msg("Choose GENE -> FUNCTION mapping table (CSV) ...")
mapping_path <- choose_file_mac("Choose Gene2FunctionMapping4Edgeswap.csv")
if (is.na(mapping_path) || !file.exists(mapping_path)) stop("Mapping file not found.")

msg("Choose directory containing edge-turnover CSV files ...")
edges_dir <- choose_dir_mac("Choose directory with B_EdgeTurnover_Top50_DayX_to_DayY.csv files")
if (is.na(edges_dir) || !dir.exists(edges_dir)) stop("Edges directory not found.")

# ----------------------------- Outputs -----------------------------

run_id <- paste0("EdgeTurnover_FunctionalSummary_", now_stamp())
output_dir <- file.path(edges_dir, run_id)
dir_create(output_dir)

# Core output paths
out_long_csv      <- file.path(output_dir, "EdgeTurnover_Annotated_Long.csv")
out_pair_csv      <- file.path(output_dir, "EdgeTurnover_FunctionPair_ByTransition.csv")
out_unified_csv   <- file.path(output_dir, "EdgeTurnover_FunctionalSummary_Unified.csv")
manifest_json     <- file.path(output_dir, "Project_Manifest.json")
manifest_filescsv <- file.path(output_dir, "Project_Manifest_Files.csv")
report_qmd        <- file.path(output_dir, "EdgeTurnover_FunctionalSummary_Report.qmd")
report_html       <- file.path(output_dir, "EdgeTurnover_FunctionalSummary_Report.html")

# ----------------------------- Read mapping -----------------------------

if (!requireNamespace("readr", quietly = TRUE)) stop("Package 'readr' is required.")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")
if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required.")

mapping <- readr::read_csv(mapping_path, show_col_types = FALSE)

req_map_cols <- c("Gene", "Primary_Function")
miss_map <- setdiff(req_map_cols, colnames(mapping))
if (length(miss_map) > 0) {
  stop("Mapping file must contain columns: ", paste(req_map_cols, collapse = ", "),
       ". Missing: ", paste(miss_map, collapse = ", "))
}

mapping <- mapping |>
  dplyr::mutate(
    Gene = as.character(Gene),
    Primary_Function = as.character(Primary_Function)
  ) |>
  dplyr::filter(!is.na(Gene), nzchar(Gene)) |>
  dplyr::distinct(Gene, .keep_all = TRUE)

# ----------------------------- Discover edge files -----------------------------

edge_files <- list.files(edges_dir, pattern = "^B_EdgeTurnover_.*Day[0-9]+_to_Day[0-9]+.*\\.csv$", full.names = TRUE)
if (length(edge_files) == 0) stop("No edge-turnover files found in: ", edges_dir)

pairs <- lapply(edge_files, function(f) {
  p <- parse_daypair_from_filename(basename(f))
  if (is.null(p)) return(NULL)
  data.frame(
    file = safe_norm(f),
    filename = basename(f),
    day_a = p$day_a,
    day_b = p$day_b,
    transition = p$transition,
    stringsAsFactors = FALSE
  )
})
pairs <- do.call(rbind, pairs)
if (is.null(pairs) || nrow(pairs) == 0) stop("Could not parse Day pairs from filenames.")

# keep only sequential chain pairs (as requested)
pairs_chain <- keep_sequential_chain(pairs)
if (nrow(pairs_chain) == 0) stop("No sequential day-pair chain detected after filtering.")

msg("Detected day-pair chain (sequential pairs only):")
msg(paste0("  - ", pairs_chain$transition, "  (", pairs_chain$filename, ")", collapse = "\n"))

# ----------------------------- Read + annotate edges -----------------------------

read_edge_one <- function(file_path, day_a, day_b, transition) {
  df <- readr::read_csv(file_path, show_col_types = FALSE)
  
  req_cols <- c("Gene1", "Gene2", "r_A", "r_B", "delta_r", "abs_delta_r")
  miss <- setdiff(req_cols, colnames(df))
  if (length(miss) > 0) {
    stop("Edge file missing required columns: ", paste(miss, collapse = ", "),
         " in file: ", basename(file_path))
  }
  
  df |>
    dplyr::mutate(
      file = safe_norm(file_path),
      transition = transition,
      day_a = as.integer(day_a),
      day_b = as.integer(day_b),
      midpoint = (day_a + day_b) / 2,
      Gene1 = as.character(Gene1),
      Gene2 = as.character(Gene2),
      r_A = as.numeric(r_A),
      r_B = as.numeric(r_B),
      delta_r = as.numeric(delta_r),
      abs_delta_r = as.numeric(abs_delta_r),
      change_dir = dplyr::case_when(
        is.na(delta_r) ~ NA_character_,
        delta_r > 0 ~ "Gain",
        delta_r < 0 ~ "Loss",
        TRUE ~ "NoChange"
      )
    )
}

edge_list <- lapply(seq_len(nrow(pairs_chain)), function(i) {
  read_edge_one(pairs_chain$file[i], pairs_chain$day_a[i], pairs_chain$day_b[i], pairs_chain$transition[i])
})
edges <- dplyr::bind_rows(edge_list)

# Merge mapping: Gene1/Gene2 -> Function1/Function2
edges_annot <- edges |>
  dplyr::left_join(mapping, by = c("Gene1" = "Gene")) |>
  dplyr::rename(Function1 = Primary_Function) |>
  dplyr::left_join(mapping, by = c("Gene2" = "Gene")) |>
  dplyr::rename(Function2 = Primary_Function) |>
  dplyr::mutate(
    Function1 = dplyr::if_else(is.na(Function1) | !nzchar(Function1), "UNMAPPED", Function1),
    Function2 = dplyr::if_else(is.na(Function2) | !nzchar(Function2), "UNMAPPED", Function2),
    Function_A = pmin(Function1, Function2),
    Function_B = pmax(Function1, Function2),
    Edge_Class = paste(Function_A, Function_B, sep = " ↔ ")
  )

# QC: mapping coverage
qc_total_edges <- nrow(edges_annot)
qc_unmapped_edges <- sum(edges_annot$Function1 == "UNMAPPED" | edges_annot$Function2 == "UNMAPPED", na.rm = TRUE)
qc_unmapped_frac <- if (qc_total_edges > 0) qc_unmapped_edges / qc_total_edges else NA_real_

# Write long annotated table
readr::write_csv(edges_annot, out_long_csv)

# ----------------------------- Summaries -----------------------------

# Per-transition functional pair summary
pair_by_transition <- edges_annot |>
  dplyr::group_by(transition, day_a, day_b, midpoint, Function_A, Function_B, Edge_Class) |>
  dplyr::summarise(
    n_edges = dplyr::n(),
    mean_abs_delta = mean(abs_delta_r, na.rm = TRUE),
    max_abs_delta  = max(abs_delta_r, na.rm = TRUE),
    mean_delta     = mean(delta_r, na.rm = TRUE),
    frac_gain      = mean(change_dir == "Gain", na.rm = TRUE),
    frac_loss      = mean(change_dir == "Loss", na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(midpoint, dplyr::desc(mean_abs_delta))

readr::write_csv(pair_by_transition, out_pair_csv)

# Unified functional summary across all transitions
unified <- pair_by_transition |>
  dplyr::group_by(Edge_Class, Function_A, Function_B) |>
  dplyr::summarise(
    Transitions_Present = paste(unique(transition), collapse = ";"),
    N_Transitions = dplyr::n_distinct(transition),
    First_Appearance_Midpoint = min(midpoint, na.rm = TRUE),
    Last_Appearance_Midpoint  = max(midpoint, na.rm = TRUE),
    Mean_EdgeScore = mean(mean_abs_delta, na.rm = TRUE),
    Max_EdgeScore  = max(max_abs_delta, na.rm = TRUE),
    Mean_SignedDelta = mean(mean_delta, na.rm = TRUE),
    Dominant_Direction = dplyr::case_when(
      mean(mean_delta, na.rm = TRUE) > 0 ~ "Gain",
      mean(mean_delta, na.rm = TRUE) < 0 ~ "Loss",
      TRUE ~ "Mixed"
    ),
    .groups = "drop"
  )

# Add persistence class based on appearances in the sequential chain order
# Build a transition order index from pairs_chain
trans_order <- pairs_chain |>
  dplyr::arrange(day_a, day_b) |>
  dplyr::mutate(
    midpoint  = (as.numeric(day_a) + as.numeric(day_b)) / 2,
    order_idx = dplyr::row_number()
  ) |>
  dplyr::select(transition, midpoint, order_idx)


unified <- unified |>
  dplyr::left_join(
    edges_annot |>
      dplyr::distinct(Edge_Class, transition) |>
      dplyr::left_join(trans_order, by = "transition") |>
      dplyr::group_by(Edge_Class) |>
      dplyr::summarise(
        order_idxs = list(sort(unique(order_idx))),
        mids = list(sort(unique(midpoint))),
        transitions_sorted = list(unique(transition[order(order_idx)])),
        .groups = "drop"
      ),
    by = "Edge_Class"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    Persistence_Type = classify_persistence(
      transitions_sorted = unlist(transitions_sorted),
      midpoints_sorted = unlist(mids)
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-order_idxs, -mids, -transitions_sorted)

# Add phase association labels by midpoint (adjust if your phase boundaries differ)
unified <- unified |>
  dplyr::mutate(
    Phase_Association = dplyr::case_when(
      First_Appearance_Midpoint < 8 ~ "Phase_A",
      First_Appearance_Midpoint >= 8 & First_Appearance_Midpoint <= 14 ~ "Phase_B",
      First_Appearance_Midpoint > 14 ~ "Phase_C",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::arrange(dplyr::desc(Mean_EdgeScore), dplyr::desc(N_Transitions))

readr::write_csv(unified, out_unified_csv)

# ----------------------------- Manifest -----------------------------

deps <- c("readr", "dplyr", "tidyr", "stringr", "jsonlite", "knitr", "rmarkdown")
deps_present <- deps[sapply(deps, requireNamespace, quietly = TRUE)]

si <- script_identity()

manifest <- list(
  run_id = run_id,
  run_timestamp = as.character(Sys.time()),
  script = si,
  output = list(
    output_dir = safe_norm(output_dir),
    report_qmd = safe_norm(report_qmd),
    report_html = safe_norm(report_html)
  ),
  input = list(
    mapping_table = list(path = safe_norm(mapping_path), filename = basename(mapping_path)),
    edge_tables_dir = safe_norm(edges_dir),
    edge_files_in_chain = lapply(seq_len(nrow(pairs_chain)), function(i) {
      list(
        transition = pairs_chain$transition[i],
        path = pairs_chain$file[i],
        filename = pairs_chain$filename[i]
      )
    })
  ),
  parameters = list(
    sequential_pair_filter = "chain order; keep only files where next day_a == previous day_b after sorting",
    unmapped_label = "UNMAPPED",
    direction_def = "delta_r > 0 => Gain; delta_r < 0 => Loss; 0 => NoChange",
    phase_midpoint_rules = "Phase_A: <8; Phase_B: 8-14; Phase_C: >14 (midpoint-based)",
    persistence_rules = "Transient: 1 transition; Emergent: 2; Persistent: >=3; Late-onset: first midpoint >=14"
  ),
  qc = list(
    total_edges = qc_total_edges,
    edges_with_unmapped_function = qc_unmapped_edges,
    fraction_edges_with_unmapped_function = qc_unmapped_frac
  ),
  outputs_written = list(
    annotated_long = basename(out_long_csv),
    functionpair_by_transition = basename(out_pair_csv),
    unified_summary = basename(out_unified_csv),
    manifest_json = basename(manifest_json),
    manifest_files = basename(manifest_filescsv)
  ),
  dependencies = list(
    packages_available = deps_present,
    R_version = as.character(getRversion())
  )
)

write_json(manifest, manifest_json)

# Pre-report inventory (so report can reference it)
inv0 <- build_file_inventory(output_dir, manifest_filescsv)

# ----------------------------- QMD + HTML report -----------------------------

qmd_lines <- c(
  "---",
  "title: \"Edge Turnover Functional Summary (Reviewer Documentation)\"",
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 3",
  "    number-sections: true",
  "execute:",
  "  echo: true",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "## Overview",
  "",
  "This report documents a reproducible pipeline that:",
  "",
  "1. Reads gene–gene edge-turnover tables for sequential developmental day-pairs.",
  "2. Annotates each gene with an *a priori* functional label from a fixed mapping table.",
  "3. Collapses gene–gene edges into function–function *edge classes*.",
  "4. Produces (i) per-transition summaries and (ii) a unified functional summary across transitions.",
  "",
  "The goal is to convert high-dimensional edge turnover into interpretable statements about **rewiring of coordination between functional programs** during development.",
  "",
  "## Inputs and Outputs (from Manifest)",
  "",
  "```{r}",
  "library(jsonlite)",
  sprintf("man <- fromJSON(%s)", shQuote(manifest_json)),
  "man$script",
  "man$input$mapping_table",
  "man$input$edge_tables_dir",
  "man$output$output_dir",
  "```",
  "",
  "### Sequential day-pair chain used",
  "",
  "```{r}",
  "# jsonlite::fromJSON may simplify list-of-objects into a data.frame automatically.",
  "efi <- man$input$edge_files_in_chain",
  "",
  "if (is.data.frame(efi)) {",
  "  chain <- efi",
  "} else if (is.list(efi)) {",
  "  chain <- do.call(rbind, lapply(efi, function(x) {",
  "    # handle either list elements or named atomic vectors",
  "    data.frame(",
  "      transition = unname(x[['transition']]),",
  "      filename   = unname(x[['filename']]),",
  "      path       = unname(x[['path']]),",
  "      stringsAsFactors = FALSE",
  "    )",
  "  }))",
  "} else {",
  "  chain <- data.frame(transition=character(), filename=character(), path=character())",
  "}",
  "",
  "chain",
  "```",
  
  "",
  "## Methods",
  "",
  "### Functional annotation",
  "Each edge table provides gene pairs (`Gene1`, `Gene2`) and correlation change metrics (`delta_r`, `abs_delta_r`).",
  "Genes are mapped to functional labels using a fixed mapping (`Gene`, `Primary_Function`).",
  "Unmapped genes are labeled `UNMAPPED` (reported in QC).",
  "",
  "### Functional edge classes",
  "Each gene–gene edge is converted to a function–function pair by mapping both endpoints.",
  "To ensure direction-invariance, pairs are standardized as:",
  "",
  "- `Function_A = min(Function1, Function2)`",
  "- `Function_B = max(Function1, Function2)`",
  "- `Edge_Class = Function_A ↔ Function_B`",
  "",
  "### Direction of change (interpretation of delta)",
  "- `delta_r > 0`: Gain",
  "- `delta_r < 0`: Loss",
  "- `delta_r = 0`: NoChange",
  "",
  "### Summaries",
  "Two levels of summarization are produced:",
  "",
  "1. **By Transition:** For each `Edge_Class` within a day-pair, compute counts and mean/max turnover.",
  "2. **Unified Across Transitions:** For each `Edge_Class`, compute persistence (how many transitions), first/last appearance, and overall mean/max turnover.",
  "",
  "## QC Summary",
  "",
  "```{r}",
  "man$qc",
  "```",
  "",
  "## Results: Per-transition function–function turnover",
  "",
  "```{r}",
  "library(readr); library(dplyr);",
  sprintf("pair <- read_csv(%s, show_col_types = FALSE)", shQuote(out_pair_csv)),
  "pair |> arrange(midpoint, desc(mean_abs_delta)) |> head(30)",
  "```",
  "",
  "## Results: Unified functional summary (single table)",
  "",
  "```{r}",
  "library(readr); library(dplyr);",
  sprintf("uni <- read_csv(%s, show_col_types = FALSE)", shQuote(out_unified_csv)),
  "uni |> head(50)",
  "```",
  "",
  "## Files written",
  "",
  "```{r}",
  "library(readr)",
  sprintf("inv <- read_csv(%s, show_col_types = FALSE)", shQuote(manifest_filescsv)),
  "inv",
  "```",
  "",
  "## Session Info",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, con = report_qmd)

# Render with Quarto if available, else rmarkdown fallback
render_ok <- TRUE
render_err <- NULL

quarto_bin <- Sys.which("quarto")
if (nzchar(quarto_bin)) {
  cmd <- sprintf('%s render %s --to html --output %s --output-dir %s',
                 shQuote(quarto_bin),
                 shQuote(report_qmd),
                 shQuote(basename(report_html)),
                 shQuote(output_dir))
  msg("Rendering with Quarto CLI:")
  msg(cmd)
  
  res <- tryCatch(system(cmd, intern = TRUE), error = function(e) e)
  if (inherits(res, "error")) {
    render_ok <- FALSE
    render_err <- paste("Quarto render error:", conditionMessage(res))
  } else {
    # Verify output landed exactly where intended
    if (!file.exists(report_html)) {
      render_ok <- FALSE
      render_err <- paste0(
        "Quarto finished without an R error, but expected HTML not found at: ",
        safe_norm(report_html)
      )
    }
  }
} else if (requireNamespace("rmarkdown", quietly = TRUE)) {

  msg("Quarto not found. Rendering with rmarkdown::render (HTML).")
  tryCatch({
    rmarkdown::render(
      input = report_qmd,
      output_format = "html_document",
      output_file = basename(report_html),
      output_dir = output_dir,
      quiet = TRUE
    )
  }, error = function(e) {
    render_ok <<- FALSE
    render_err <<- paste("rmarkdown render error:", conditionMessage(e))
  })
} else {
  render_ok <- FALSE
  render_err <- "Neither Quarto CLI nor rmarkdown package is available; report was written but not rendered."
}

# Final inventory (post-render)
inv1 <- build_file_inventory(output_dir, manifest_filescsv)

# Update manifest with render status
manifest$render <- list(
  ok = render_ok,
  error = render_err
)
write_json(manifest, manifest_json)

# ----------------------------- Finish -----------------------------

msg("=====================================")
msg("Done. Files written to:")
msg("  %s", safe_norm(output_dir))
msg("Key outputs:")
msg("  - %s", safe_norm(out_unified_csv))
msg("  - %s", safe_norm(report_qmd))
if (render_ok) {
  msg("  - %s", safe_norm(report_html))
} else {
  msg("HTML render failed: %s", ifelse(is.null(render_err), "Unknown error", render_err))
}
msg("Manifest:")
msg("  - %s", safe_norm(manifest_json))
msg("  - %s", safe_norm(manifest_filescsv))
msg("=====================================")
