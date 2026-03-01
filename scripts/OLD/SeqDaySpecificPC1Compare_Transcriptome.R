#!/usr/bin/env Rscript
# ======================================================================
# Sequential Day-Pair DAY-SPECIFIC PC1 Comparison (Transcriptome)
#
# PURPOSE (what this script is ACTUALLY useful for)
#
# This script is designed to answer the exact question:
#   “How does the dominant WITHIN-DAY coordination axis (PC1) change
#    as the system moves through sequential days?”
#
# It does this by running PCA SEPARATELY within each day, then producing
# a gene-level comparison table for each sequential day-pair (Day1 -> Day2):
#
#   For the shared genes, it outputs:
#     - PC1 loading on Day1 (from Day1-only PCA)
#     - PC1 loading on Day2 (from Day2-only PCA), sign-aligned to Day1 PC1
#     - PC1 loading change (Day2 - Day1)
#     - Δmean expression (mean Day2 - mean Day1) for directional context
#
# KEY INTERPRETATION
#   - PC1 loadings describe a gene’s participation in the dominant coordinated
#     expression program *within that day* (i.e., within-day organization).
#   - Comparing Day1-PC1 vs Day2-PC1 (after sign alignment) tells you whether
#     the dominant coordination structure is conserved (high similarity) or
#     reorganized/rotated (low similarity).
#   - PCA itself is NOT temporally directional; Δmean provides time direction.
#
# WHAT THIS SCRIPT IS NOT
#   - It is NOT a JOINT PCA on pooled days (that answers a different “contrast axis” question).
#   - It is NOT a sample-score table. Outputs are gene-level loading comparisons.
#
# INPUT FORMAT
#   One tidy-long transcriptome table with columns:
#     Day, SampleID, Gene, Abundance
#
# OUTPUTS (STRICT)
#   All outputs are written ONLY under:
#     outputs/<run_id>/
#
# KEY OUTPUT FILES (human-interpretable names)
#   For each sequential pair DayXX_to_DayYY:
#     - GENE_PC1COMPARE__DayXX_to_DayYY__SharedGenes.csv
#     - DAYPCA_Loadings__DayXX__PC1_PC2.csv
#     - DAYPCA_Loadings__DayYY__PC1_PC2.csv
#     - PLOT_PC1LoadingsScatter__DayXX_to_DayYY.png
#     - DIR_DeltaMean__DayYY_minus_DayXX.csv
#
#   Across all pairs:
#     - SUMMARY_SequentialPairs__DaySpecificPC1Similarity.csv
#
#   Provenance / reproducibility:
#     - Project_Manifest.json
#     - Project_Manifest_Files.csv
#     - <script_name>_Report.qmd
#     - <script_name>_Report.html (rendered as final step)
#
# ======================================================================

# --------------------------- Dependency handling ---------------------------

quiet_install <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message(sprintf("Installing missing package: %s", p))
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

required_pkgs <- c(
  "jsonlite", "knitr", "rmarkdown",
  "dplyr", "tidyr", "readr", "ggplot2"
)

quiet_install(required_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
})

stop2 <- function(...) stop(paste0(...), call. = FALSE)

# --------------------------- Script identity capture (MANDATORY) ---------------------------

resolve_script_path <- function() {
  # 1) RStudio editor context
  p <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      rstudioapi::getSourceEditorContext()$path
    } else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) {
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  # 2) Rscript --file=
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) {
      return(normalizePath(p2, winslash = "/", mustWork = FALSE))
    }
  }
  
  # 3) source() fallback (sometimes present)
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) {
    return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  }
  
  # 4) Unknown
  NA_character_
}

known_script_filename <- "SeqDaySpecificPC1Compare_Transcriptome.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  message("NOTE: Script path detection failed; using fallback filename: ", known_script_filename)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

message("---- Script Identity ----")
message("script_name: ", script_name)
message("script_path: ", script_path)
message("script_full: ", ifelse(is.na(script_full), "NA", script_full))
message("-------------------------")

# --------------------------- Output directory policy (MANDATORY) ---------------------------

analysis_name <- "SeqDaySpecificPC1Compare"
source_label  <- "Transcriptome"

timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0(analysis_name, "_", source_label, "_", timestamp_tag)

outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- Input handling ---------------------------

# Optional command-line args:
#   --input=/path/to/file.csv
#   --scale=TRUE/FALSE
#   --center=TRUE/FALSE
#   --anchor=THRSP
args <- commandArgs(trailingOnly = FALSE)

get_arg_value <- function(prefix) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", prefix, "="), "", hit[1])
}

input_csv <- get_arg_value("--input")
center_flag <- get_arg_value("--center")
scale_flag  <- get_arg_value("--scale")
anchor_gene_user <- get_arg_value("--anchor")

center <- if (!is.na(center_flag)) tolower(center_flag) %in% c("true", "t", "1", "yes", "y") else TRUE
scale_ <- if (!is.na(scale_flag))  tolower(scale_flag)  %in% c("true", "t", "1", "yes", "y") else TRUE

if (is.na(input_csv) || !nzchar(input_csv) || !file.exists(input_csv)) {
  message("Choose input tidy-long transcriptome CSV (must contain: Day, SampleID, Gene, Abundance)")
  input_csv <- tryCatch(file.choose(), error = function(e) NA_character_)
  if (is.na(input_csv) || !file.exists(input_csv)) stop2("No valid input file selected. Aborting.")
}
input_csv <- normalizePath(input_csv, winslash = "/", mustWork = FALSE)

# Anchor gene preference order
anchor_preference <- c()
if (!is.na(anchor_gene_user) && nzchar(anchor_gene_user)) {
  anchor_preference <- c(anchor_preference, anchor_gene_user)
}
anchor_preference <- c(anchor_preference, "THRSP", "HMGCS2")

# --------------------------- Read and validate input ---------------------------

df <- suppressMessages(readr::read_csv(input_csv, show_col_types = FALSE))

needed_cols <- c("Day", "SampleID", "Gene", "Abundance")
missing_cols <- setdiff(needed_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop2("Input file is missing required columns: ", paste(missing_cols, collapse = ", "),
        "\nExpected columns: ", paste(needed_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    Day = as.integer(Day),
    SampleID = as.character(SampleID),
    Gene = as.character(Gene),
    Abundance = as.numeric(Abundance)
  )

if (anyNA(df$Day)) stop2("Day contains NA after coercion to integer.")
if (anyNA(df$Abundance)) {
  message("NOTE: Abundance contains NA values; matrices will drop genes with incomplete cases per day.")
}

days <- sort(unique(df$Day))
if (length(days) < 2) stop2("Need at least 2 unique Day values to run sequential comparisons.")

day_pairs <- data.frame(
  Day1 = days[-length(days)],
  Day2 = days[-1],
  stringsAsFactors = FALSE
)

# --------------------------- Utilities ---------------------------

get_pkg_versions <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    data.frame(package = p, version = v, stringsAsFactors = FALSE)
  })
  bind_rows(out)
}

choose_anchor_gene <- function(gene_vec, anchor_pref) {
  for (g in anchor_pref) {
    if (g %in% gene_vec) return(g)
  }
  gene_vec[1]
}

cosine_similarity <- function(a, b) {
  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(a * b) / denom
}

write_plot_safe <- function(p, path, width = 7.5, height = 5.5, dpi = 300) {
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi, units = "in")
  if (!file.exists(path)) stop2("Failed to write plot: ", path)
}

# Build sample-by-gene matrix for a single day
# - conservative filtering: complete cases by gene across samples
# - drop zero-variance genes when scaling
build_sample_by_gene_matrix <- function(df_day) {
  mat_gene_sample <- df_day %>%
    select(Gene, SampleID, Abundance) %>%
    pivot_wider(names_from = SampleID, values_from = Abundance)
  
  if (nrow(mat_gene_sample) == 0) return(NULL)
  
  genes <- mat_gene_sample$Gene
  mat <- as.data.frame(mat_gene_sample[, -1, drop = FALSE])
  rownames(mat) <- genes
  
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 10) return(NULL)
  
  X <- t(as.matrix(mat))  # samples x genes
  
  # drop non-finite and zero-variance columns (required if scale.=TRUE)
  col_sd <- apply(X, 2, sd, na.rm = TRUE)
  keep <- is.finite(col_sd) & (col_sd > 0)
  X <- X[, keep, drop = FALSE]
  
  if (ncol(X) < 10 || nrow(X) < 2) return(NULL)
  X
}

pca_var_expl <- function(pca) {
  v <- (pca$sdev^2) / sum(pca$sdev^2)
  list(ve1 = v[1], ve2 = if (length(v) >= 2) v[2] else NA_real_)
}

# Anchor PC1 sign within a day using an anchor gene (if present)
anchor_day_loadings <- function(load_df, anchor_pref) {
  # load_df must have Gene, PC1, PC2 (PC2 may be NA)
  anchor_used <- choose_anchor_gene(load_df$Gene, anchor_pref)
  flipped <- FALSE
  v <- load_df$PC1[match(anchor_used, load_df$Gene)]
  if (!is.na(v) && v < 0) {
    load_df$PC1 <- -load_df$PC1
    if ("PC2" %in% names(load_df)) load_df$PC2 <- -load_df$PC2
    flipped <- TRUE
  }
  list(loadings = load_df, anchor_used = anchor_used, flipped = flipped)
}

# Align Day2 PC1 to Day1 PC1 using correlation on shared genes (post-anchor)
align_day2_to_day1 <- function(v_day1, v_day2) {
  r0 <- suppressWarnings(stats::cor(v_day1, v_day2, use = "pairwise.complete.obs"))
  flipped <- FALSE
  v2 <- v_day2
  if (!is.na(r0) && r0 < 0) {
    v2 <- -v2
    flipped <- TRUE
  }
  list(v2 = v2, flipped = flipped, cor_before = r0)
}

# --------------------------- Run sequential comparisons ---------------------------

pair_summaries <- list()

for (i in seq_len(nrow(day_pairs))) {
  d1 <- day_pairs$Day1[i]
  d2 <- day_pairs$Day2[i]
  pair_tag <- sprintf("Day%02d_to_Day%02d", d1, d2)
  
  message("---- Running DAY-SPECIFIC PCA PC1 comparison for ", pair_tag, " ----")
  
  df_d1 <- df %>% filter(Day == d1)
  df_d2 <- df %>% filter(Day == d2)
  
  X1 <- build_sample_by_gene_matrix(df_d1)
  X2 <- build_sample_by_gene_matrix(df_d2)
  
  if (is.null(X1) || is.null(X2)) {
    message("Skipping ", pair_tag, ": insufficient data after filtering for one or both days.")
    next
  }
  
  pca1 <- prcomp(X1, center = center, scale. = scale_)
  pca2 <- prcomp(X2, center = center, scale. = scale_)
  
  ve1 <- pca_var_expl(pca1)
  ve2 <- pca_var_expl(pca2)
  
  message(sprintf(
    "Variance explained by PC1: Day %02d = %.2f%% | Day %02d = %.2f%%",
    d1, 100 * ve1$ve1,
    d2, 100 * ve2$ve1
  ))
  
  # Loadings tables (genes)
  L1 <- as.data.frame(pca1$rotation)
  L1$Gene <- rownames(L1)
  rownames(L1) <- NULL
  L1 <- L1 %>% transmute(Gene, PC1 = .data$PC1, PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_)
  
  L2 <- as.data.frame(pca2$rotation)
  L2$Gene <- rownames(L2)
  rownames(L2) <- NULL
  L2 <- L2 %>% transmute(Gene, PC1 = .data$PC1, PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_)
  
  # Anchor within-day sign (so “PC1 positive pole” is consistent across runs)
  anch1 <- anchor_day_loadings(L1, anchor_preference)
  anch2 <- anchor_day_loadings(L2, anchor_preference)
  L1a <- anch1$loadings
  L2a <- anch2$loadings
  
  # Write day-specific loadings (human readable)
  day1_load_path <- file.path(output_dir, sprintf("DAYPCA_Loadings__Day%02d__PC1_PC2.csv", d1))
  day2_load_path <- file.path(output_dir, sprintf("DAYPCA_Loadings__Day%02d__PC1_PC2.csv", d2))
  # These files may be overwritten across pairs if days repeat; avoid overwrite by including run-wide uniqueness:
  # Add a suffix noting it was generated in context of this pair.
  day1_load_path <- file.path(output_dir, sprintf("DAYPCA_Loadings__Day%02d__PC1_PC2__Context_%s.csv", d1, pair_tag))
  day2_load_path <- file.path(output_dir, sprintf("DAYPCA_Loadings__Day%02d__PC1_PC2__Context_%s.csv", d2, pair_tag))
  readr::write_csv(L1a %>% arrange(desc(abs(PC1))), day1_load_path)
  readr::write_csv(L2a %>% arrange(desc(abs(PC1))), day2_load_path)
  
  # Shared genes and sign-aligned comparison
  shared <- intersect(L1a$Gene, L2a$Gene)
  
  if (length(shared) < 50) {
    message("Skipping gene-level PC1 compare for ", pair_tag, ": too few shared genes (", length(shared), ").")
    next
  }
  
  v1 <- L1a$PC1[match(shared, L1a$Gene)]
  v2 <- L2a$PC1[match(shared, L2a$Gene)]
  
  aligned <- align_day2_to_day1(v1, v2)
  v2a <- aligned$v2
  flipped_day2_to_day1 <- aligned$flipped
  
  cor_pc1 <- suppressWarnings(stats::cor(v1, v2a, use = "pairwise.complete.obs"))
  cos_pc1 <- cosine_similarity(v1, v2a)
  angle_deg <- NA_real_
  if (is.finite(cos_pc1)) {
    cc <- max(min(cos_pc1, 1), -1)
    angle_deg <- acos(cc) * 180 / pi
  }
  
  # Directional context: Δmean (Day2 - Day1)
  mean_tbl <- df %>%
    filter(Day %in% c(d1, d2)) %>%
    group_by(Day, Gene) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Day, values_from = mean_abundance)
  
  col_d1 <- as.character(d1)
  col_d2 <- as.character(d2)
  if (!(col_d1 %in% names(mean_tbl))) mean_tbl[[col_d1]] <- NA_real_
  if (!(col_d2 %in% names(mean_tbl))) mean_tbl[[col_d2]] <- NA_real_
  
  mean_tbl <- mean_tbl %>%
    mutate(
      delta_mean = .data[[col_d2]] - .data[[col_d1]]
    ) %>%
    rename(
      mean_Day1 = !!col_d1,
      mean_Day2 = !!col_d2
    ) %>%
    select(Gene, mean_Day1, mean_Day2, delta_mean)
  
  delta_path <- file.path(output_dir, sprintf("DIR_DeltaMean__Day%02d_minus_Day%02d.csv", d2, d1))
  readr::write_csv(mean_tbl %>% arrange(desc(abs(delta_mean))), delta_path)
  
  # Build the exact “what you want” table
  comp_tbl <- data.frame(
    Gene = shared,
    PC1_Day1 = v1,
    PC1_Day2_aligned = v2a,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      PC1_change_Day2_minus_Day1 = PC1_Day2_aligned - PC1_Day1
    ) %>%
    left_join(mean_tbl, by = "Gene") %>%
    arrange(desc(abs(PC1_change_Day2_minus_Day1)))
  
  gene_comp_path <- file.path(output_dir, sprintf("GENE_PC1COMPARE__Day%02d_to_Day%02d__SharedGenes.csv", d1, d2))
  readr::write_csv(comp_tbl, gene_comp_path)
  
  # Plot: PC1 loading scatter (Day1 vs Day2)
  p_scatter <- ggplot(comp_tbl, aes(x = PC1_Day1, y = PC1_Day2_aligned)) +
    geom_point(alpha = 0.45) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    labs(
      title = sprintf("Day-specific PC1 loadings: Day %d vs Day %d (shared genes)", d1, d2),
      subtitle = paste0(
        "Shared genes=", nrow(comp_tbl),
        " | cor=", sprintf("%.3f", cor_pc1),
        " | cos=", sprintf("%.3f", cos_pc1),
        " | angle=", sprintf("%.1f°", angle_deg),
        " | Anchor Day1=", anch1$anchor_used,
        " | Anchor Day2=", anch2$anchor_used,
        if (flipped_day2_to_day1) " | Day2 PC1 flipped to align with Day1" else ""
      ),
      x = paste0("PC1 loading (Day ", d1, ")"),
      y = paste0("PC1 loading (Day ", d2, ", aligned)")
    )
  
  scatter_path <- file.path(output_dir, sprintf("PLOT_PC1LoadingsScatter__Day%02d_to_Day%02d.png", d1, d2))
  write_plot_safe(p_scatter, scatter_path)
  
  # Pair summary (human interpretable)
  pair_summaries[[pair_tag]] <- data.frame(
    Pair = pair_tag,
    Day1 = d1,
    Day2 = d2,
    center = center,
    scale = scale_,
    n_shared_genes = nrow(comp_tbl),
    Day1_anchor_gene = anch1$anchor_used,
    Day1_anchor_flipped = anch1$flipped,
    Day2_anchor_gene = anch2$anchor_used,
    Day2_anchor_flipped = anch2$flipped,
    Day2_flipped_to_align_with_Day1 = flipped_day2_to_day1,
    PC1_similarity_correlation = cor_pc1,
    PC1_similarity_cosine = cos_pc1,
    PC1_similarity_angle_deg = angle_deg,
    var_expl_Day1_PC1 = ve1$ve1,
    var_expl_Day2_PC1 = ve2$ve1,
    var_expl_Day1_PC1_pct = 100 * ve1$ve1,
    var_expl_Day2_PC1_pct = 100 * ve2$ve1,
    gene_compare_table = basename(gene_comp_path),
    delta_mean_table = basename(delta_path),
    scatter_plot = basename(scatter_path),
    stringsAsFactors = FALSE
  )
}

if (length(pair_summaries) == 0) {
  stop2("No sequential pairs produced outputs. Check input data completeness and filtering.")
}

summary_tbl <- bind_rows(pair_summaries)
summary_path <- file.path(output_dir, "SUMMARY_SequentialPairs__DaySpecificPC1Similarity.csv")
readr::write_csv(summary_tbl, summary_path)

# --------------------------- Manifest (MANDATORY) ---------------------------

deps_tbl <- get_pkg_versions(required_pkgs)

manifest_path_json <- file.path(output_dir, "Project_Manifest.json")
manifest_path_csv  <- file.path(output_dir, "Project_Manifest_Files.csv")

# Script header embedding for QMD
header_text <- NULL
if (!is.na(script_full) && file.exists(script_full)) {
  lines <- readLines(script_full, warn = FALSE)
  keep_idx <- which(grepl("^#!", lines) | grepl("^\\s*$", lines) | grepl("^\\s*#", lines))
  if (length(keep_idx) > 0) {
    last <- 0
    for (k in keep_idx) {
      if (k == last + 1) last <- k else break
    }
    header_text <- paste(lines[seq_len(last)], collapse = "\n")
  }
}
if (is.null(header_text) || !nzchar(header_text)) {
  header_text <- paste(
    "# Sequential Day-Pair DAY-SPECIFIC PC1 Comparison (Transcriptome)",
    "# (Header extraction fallback: script path not available or unreadable.)",
    sep = "\n"
  )
}

manifest <- list(
  run_id = run_id,
  run_timestamp = as.character(Sys.time()),
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA_character_, script_full)
  ),
  input = list(
    transcriptome_csv = input_csv
  ),
  parameters = list(
    analysis_name = analysis_name,
    source_label = source_label,
    days_found = days,
    sequential_pairs = day_pairs,
    center = center,
    scale = scale_,
    anchor_preference = anchor_preference,
    outputs_description = list(
      gene_pc1_compare = "Shared-gene table with Day1 PC1 and Day2 PC1 (aligned), PC1 change, and Δmean",
      day_loadings = "Day-specific PCA loadings (PC1/PC2) for context within each pair",
      delta_mean = "Directional context: mean(Day2) - mean(Day1) per gene",
      similarity_summary = "Across-pair similarity metrics (cor/cos/angle) for Day-specific PC1"
    )
  ),
  dependencies = lapply(seq_len(nrow(deps_tbl)), function(i) {
    list(package = deps_tbl$package[i], version = deps_tbl$version[i])
  }),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
)

writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"), manifest_path_json)

# --------------------------- File inventory (MANDATORY) ---------------------------

write_inventory_csv <- function(out_dir, out_csv) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  rel <- sub(
    paste0("^", gsub("\\\\", "/", normalizePath(out_dir, winslash = "/", mustWork = FALSE)), "/?"),
    "",
    gsub("\\\\", "/", files)
  )
  info <- file.info(files)
  inv <- data.frame(
    file = rel,
    size_bytes = as.numeric(info$size),
    mtime = as.character(info$mtime),
    stringsAsFactors = FALSE
  ) %>% arrange(file)
  readr::write_csv(inv, out_csv)
  inv
}

inv0 <- write_inventory_csv(output_dir, manifest_path_csv)

# --------------------------- Quarto QMD (MANDATORY) ---------------------------

qmd_path  <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))

qmd_text <- paste0(
  "---\n",
  "title: \"", script_name, " Report\"\n",
  "format:\n",
  "  html:\n",
  "    toc: true\n",
  "execute:\n",
  "  echo: false\n",
  "  warning: false\n",
  "  message: false\n",
  "---\n\n",
  
  "# Run metadata\n\n",
  "```{r}\n",
  "library(jsonlite)\n",
  "library(knitr)\n",
  "library(dplyr)\n",
  "manifest <- jsonlite::fromJSON('Project_Manifest.json', simplifyVector = FALSE)\n",
  "inv <- read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE)\n",
  "meta_tbl <- data.frame(\n",
  "  Field = c('run_id','run_timestamp','script_name','script_path','script_full_path','outputs_root','output_dir','input_csv','center','scale','anchor_preference'),\n",
  "  Value = c(\n",
  "    manifest[['run_id']],\n",
  "    manifest[['run_timestamp']],\n",
  "    manifest[['script']][['name']],\n",
  "    manifest[['script']][['path']],\n",
  "    ifelse(is.null(manifest[['script']][['full_path']]) || is.na(manifest[['script']][['full_path']]), 'NA', manifest[['script']][['full_path']]),\n",
  "    manifest[['outputs']][['outputs_root']],\n",
  "    manifest[['outputs']][['output_dir']],\n",
  "    manifest[['input']][['transcriptome_csv']],\n",
  "    as.character(manifest[['parameters']][['center']]),\n",
  "    as.character(manifest[['parameters']][['scale']]),\n",
  "    paste(manifest[['parameters']][['anchor_preference']], collapse=', ')\n",
  "  ),\n",
  "  stringsAsFactors = FALSE\n",
  ")\n",
  "knitr::kable(meta_tbl)\n",
  "```\n\n",
  
  "# Script header\n\n",
  "```r\n",
  gsub("```", "``\\`", header_text),
  "\n```\n\n",
  
  "# What this analysis does\n\n",
  "This report documents **day-specific PCA** comparisons across sequential days.\n\n",
  "## Concept: day-specific PC1 as within-day organization\n\n",
  "For each day, PCA is computed using only samples from that day. The **PC1 loading vector** describes how genes participate in the dominant coordinated expression program within that day.\n\n",
  "## Comparing Day1 PC1 to Day2 PC1\n\n",
  "For each sequential day-pair, we restrict to **shared genes** and compare the PC1 loading vectors.\n\n",
  "- Sign ambiguity is handled by:\n",
  "  1) anchoring each day’s PC1 so an anchor gene has positive loading (if present), then\n",
  "  2) flipping Day2 PC1 if needed so its correlation with Day1 PC1 is positive.\n\n",
  "Similarity is summarized using correlation, cosine similarity, and the implied angle:\n\n",
  "\\[\n",
  "\\cos(\\theta) = \\frac{v_1 \\cdot v_2}{\\lVert v_1 \\rVert \\lVert v_2 \\rVert}, \\quad \\theta = \\arccos(\\cos(\\theta)).\n",
  "\\]\n\n",
  "## Directionality (Δmean)\n\n",
  "Because PCA is not intrinsically directional in time, we add gene-level mean change:\n\n",
  "\\[\n",
  "\\Delta\\text{mean}(g) = \\text{mean}_{Day2}(g) - \\text{mean}_{Day1}(g).\n",
  "\\]\n\n",
  
  "# Dependencies\n\n",
  "```{r}\n",
  "deps <- manifest[['dependencies']]\n",
  "deps_tbl <- do.call(rbind, lapply(deps, function(x) data.frame(package=x$package, version=x$version, stringsAsFactors=FALSE)))\n",
  "knitr::kable(deps_tbl)\n",
  "```\n\n",
  
  "# Generated files\n\n",
  "```{r}\n",
  "knitr::kable(inv)\n",
  "```\n\n",
  
  "# Summary across sequential pairs\n\n",
  "```{r}\n",
  "summary_path <- 'SUMMARY_SequentialPairs__DaySpecificPC1Similarity.csv'\n",
  "if (file.exists(summary_path)) {\n",
  "  s <- read.csv(summary_path, stringsAsFactors = FALSE)\n",
  "  knitr::kable(s)\n",
  "} else {\n",
  "  cat('Summary table not found.')\n",
  "}\n",
  "```\n\n",
  
  "# Example plot\n\n",
  "```{r}\n",
  "pngs <- inv$file[grepl('^PLOT_PC1LoadingsScatter__.*\\\\.png$', inv$file)]\n",
  "if (length(pngs) > 0) {\n",
  "  knitr::include_graphics(pngs[1])\n",
  "} else {\n",
  "  cat('No scatter plots found to display.')\n",
  "}\n",
  "```\n\n",
  
  "# Interpretation guide\n\n",
  "- **High similarity (high cor/cos, small angle):** dominant within-day coordination is conserved from Day1 to Day2.\n",
  "- **Low similarity (low cor/cos, large angle):** dominant coordination structure reorganized/rotated between days.\n",
  "- **Δmean:** indicates whether the gene increased/decreased in average abundance, which is distinct from coordination structure.\n\n",
  
  "# Reproducibility\n\n",
  "```{r}\n",
  "sessionInfo()\n",
  "```\n"
)

writeLines(qmd_text, qmd_path)

# --------------------------- Render report as final step (MANDATORY) ---------------------------

render_ok <- TRUE
render_msg <- NULL

tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_format = "html_document",
    output_file = basename(html_path),
    output_dir = output_dir,
    quiet = TRUE
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# Refresh inventory AFTER rendering (MANDATORY)
inv1 <- write_inventory_csv(output_dir, manifest_path_csv)

if (render_ok && file.exists(html_path)) {
  message("HTML report created: ", normalizePath(html_path, winslash = "/", mustWork = FALSE))
} else {
  message("ERROR: HTML report render failed.")
  if (!is.null(render_msg)) message("Reason: ", render_msg)
  message("QMD written at: ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE))
  message("Check outputs in: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
}

message("Run complete. Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
