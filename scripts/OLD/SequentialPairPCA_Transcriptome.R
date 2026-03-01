#!/usr/bin/env Rscript
# ======================================================================
# Sequential Day-Pair PCA Toolkit (Transcriptome): Joint PCA + Day-Specific PCA + Subspace Similarity + Directional Context
#
# PURPOSE (What this script is ACTUALLY useful for)
#
# This script performs FOUR complementary analyses for each sequential day-pair (Day1, Day2):
#
# (A) JOINT PCA (Day1 + Day2 pooled; “contrast axis”)
#     - What it answers:
#         • Which coordinated gene program most strongly DISTINGUISHES Day1 from Day2?
#         • Which genes have the largest |PC1| loadings on the pooled data and therefore
#           define the dominant BETWEEN-day contrast?
#     - What it is NOT:
#         • It is not “PC1 of Day1 compared to PC1 of Day2.”
#         • It does not provide intrinsic temporal directionality (PCA is symmetric).
#     - When to use it:
#         • Identify transition markers / discriminating programs between adjacent stages.
#
# (B) DAY-SPECIFIC PCA + PC1 COMPARISON (separate PCA within each day; “organization axis”)
#     - What it answers:
#         • How does the INTERNAL coordination structure change from Day1 to Day2?
#           (i.e., does the dominant within-day axis rotate / reorganize?)
#         • Are Day1-PC1 and Day2-PC1 similar (continuity) or dissimilar (reorganization)?
#     - How we quantify similarity:
#         • Cosine similarity and angle between Day1 PC1 and Day2 PC1 loadings on shared genes
#         • Correlation between Day1 and Day2 PC1 loadings (after sign alignment)
#     - Caveat:
#         • PC1 alone can be unstable if eigenvalues are close (axes can swap or rotate).
#
# (C) DAY-SPECIFIC SUBSPACE SIMILARITY (PC1–PCK; “dominant manifold”)
#     - What it answers:
#         • Is the dominant low-dimensional coordination structure conserved across days,
#           even if individual PCs reorder or rotate?
#     - Why it matters:
#         • During reorganization, Day1-PC1 may become Day2-PC2 (or be distributed across PC1–PC3).
#           Comparing ONLY PC1 can therefore overstate reorganization.
#         • Comparing the K-dimensional subspace spanned by PC1..PCK is often more stable.
#     - How we quantify:
#         • Principal angles between the two K-dimensional subspaces, derived from SVD.
#         • Summaries: mean angle, max angle, and mean cos^2(angle) similarity.
#
# (D) DIRECTIONAL CONTEXT (adds “time arrow” not provided by PCA)
#     - PCA itself is non-directional (sign and axis orientation are arbitrary).
#       To add time ordering, we compute gene-level mean change:
#           Δmean(gene) = mean_Day2(gene) - mean_Day1(gene)
#     - We then report:
#         • correlation( joint_PC1_loading , Δmean ) across genes
#         • top genes with largest positive/negative Δmean
#
# INPUT FORMAT
#   One tidy-long transcriptome table with columns:
#     Day, SampleID, Gene, Abundance
#
# OUTPUTS (STRICT)
#   All outputs are written ONLY under:
#     outputs/<run_id>/
#
# HUMAN-INTERPRETABLE FILE NAMING STANDARD
#   - A_JointPCA__... : pooled Day1+Day2 “contrast axis” outputs
#   - B_DayPCA__...   : within-day PCA outputs (Day1 and Day2 separately) + PC1 comparison
#   - C_Subspace__... : PC1..PCk subspace similarity (principal angles)
#   - D_Direction__...: Δmean tables and directional summaries
#   - Z_Manifest__... : manifests / inventories / README mapping files to meaning
#   - Z_Summary__...  : cross-pair summary tables
#
# KEY OUTPUTS PER PAIR
#   - (A) Joint PCA:
#       • A_JointPCA__Pair_Day04_vs_Day06__GeneLoadings_PC1_PC2.csv
#       • A_JointPCA__Pair_Day04_vs_Day06__SampleScores_PC1_PC2_withDayLabel.csv
#       • A_JointPCA__Pair_Day04_vs_Day06__Plot_PC1ScoresByDay_BoxJitter.png
#       • A_JointPCA__Pair_Day04_vs_Day06__Plot_Top25GeneLoadingsByAbsPC1_Bar.png
#   - (B) Day-specific PCA:
#       • B_DayPCA__Pair_Day04_vs_Day06__Within_Day04__GeneLoadings_PC1_PC2.csv
#       • B_DayPCA__Pair_Day04_vs_Day06__Within_Day06__GeneLoadings_PC1_PC2.csv
#       • B_DayPCA__Pair_Day04_vs_Day06__Plot_Day04_PC1vsDay06_PC1_LoadingScatter_Aligned.png
#   - (C) Subspace similarity:
#       • C_Subspace__Pair_Day04_vs_Day06__PrincipalAngles_Subspace_PC1_to_PC5.csv
#   - (D) Directional context:
#       • D_Direction__Pair_Day04_vs_Day06__GeneDeltaMean_Day06_minus_Day04.csv
#
# PROJECT-LEVEL OUTPUTS
#   - Z_Summary__AllPairs__Metrics_JointPCA_DayPCA_Subspace_Direction.csv
#   - Z_Manifest__Project_Manifest.json
#   - Z_Manifest__FileInventory.csv
#   - Z_Manifest__README_FileDescriptions.csv
#   - <script_name>_Report.qmd and <script_name>_Report.html
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

known_script_filename <- "SequentialPairPCA_Transcriptome.R"
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

# Console summary (MANDATORY)
message("---- Script Identity ----")
message("script_name: ", script_name)
message("script_path: ", script_path)
message("script_full: ", ifelse(is.na(script_full), "NA", script_full))
message("-------------------------")

# --------------------------- Output directory policy (MANDATORY) ---------------------------

analysis_name <- "SeqPairPCA_Toolkit"
source_label  <- "Transcriptome"

timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0(analysis_name, "_", source_label, "_", timestamp_tag)

outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("---- Output Policy ----")
message("All outputs will be written under: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message("-----------------------")

# --------------------------- Input handling ---------------------------

# Optional command-line args:
#   --input=/path/to/file.csv
#   --scale=TRUE/FALSE
#   --center=TRUE/FALSE
#   --anchor=THRSP
#   --k=5   (for subspace similarity PC1..PCK)
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
k_flag <- get_arg_value("--k")

center <- if (!is.na(center_flag)) tolower(center_flag) %in% c("true", "t", "1", "yes", "y") else TRUE
scale_ <- if (!is.na(scale_flag))  tolower(scale_flag)  %in% c("true", "t", "1", "yes", "y") else TRUE

K_subspace <- 5L
if (!is.na(k_flag) && nzchar(k_flag)) {
  K_subspace <- suppressWarnings(as.integer(k_flag))
  if (!is.finite(K_subspace) || is.na(K_subspace) || K_subspace < 1) {
    stop2("Invalid --k value: ", k_flag, " (must be integer >= 1)")
  }
}

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

message("---- Run Parameters ----")
message("input_csv: ", input_csv)
message("center: ", center)
message("scale: ", scale_)
message("K_subspace (requested): ", K_subspace)
message("anchor_preference: ", paste(anchor_preference, collapse = " -> "))
message("------------------------")

# --------------------------- Read and validate input ---------------------------

df <- suppressMessages(readr::read_csv(input_csv, show_col_types = FALSE))

needed_cols <- c("Day", "SampleID", "Gene", "Abundance")
missing_cols <- setdiff(needed_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop2(
    "Input file is missing required columns: ", paste(missing_cols, collapse = ", "),
    "\nExpected columns: ", paste(needed_cols, collapse = ", ")
  )
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
  message("NOTE: Abundance contains NA values; matrices will drop genes with incomplete cases per analysis.")
}

days <- sort(unique(df$Day))
if (length(days) < 2) stop2("Need at least 2 unique Day values to run sequential pair analyses.")

day_pairs <- data.frame(
  Day1 = days[-length(days)],
  Day2 = days[-1],
  stringsAsFactors = FALSE
)

message("Days found: ", paste(days, collapse = ", "))
message("Sequential day-pairs: ", paste(sprintf("%s-%s", day_pairs$Day1, day_pairs$Day2), collapse = ", "))

# --------------------------- PCA utilities ---------------------------

choose_anchor_gene <- function(gene_vec, anchor_pref) {
  for (g in anchor_pref) {
    if (g %in% gene_vec) return(g)
  }
  gene_vec[1]
}

anchor_vector_sign_by_gene <- function(loadings_df, anchor_pref) {
  anchor_used <- choose_anchor_gene(loadings_df$Gene, anchor_pref)
  v <- loadings_df$PC1[match(anchor_used, loadings_df$Gene)]
  flipped <- FALSE
  if (!is.na(v) && v < 0) {
    loadings_df$PC1 <- -loadings_df$PC1
    if ("PC2" %in% names(loadings_df)) loadings_df$PC2 <- -loadings_df$PC2
    flipped <- TRUE
  }
  list(loadings = loadings_df, anchor_used = anchor_used, flipped = flipped)
}

align_sign_to_reference <- function(v_ref, v_target) {
  if (length(v_ref) != length(v_target)) stop2("Vector length mismatch in align_sign_to_reference().")
  r <- suppressWarnings(stats::cor(v_ref, v_target, use = "pairwise.complete.obs"))
  flipped <- FALSE
  v2 <- v_target
  if (!is.na(r) && r < 0) {
    v2 <- -v2
    flipped <- TRUE
  }
  list(v = v2, flipped = flipped, cor_before = r)
}

cosine_similarity <- function(a, b) {
  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(a * b) / denom
}

write_plot_safe <- function(p, path, width = 8, height = 5, dpi = 300) {
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi, units = "in")
  if (!file.exists(path)) stop2("Failed to write plot: ", path)
}

build_sample_by_gene_matrix <- function(df_in) {
  # Returns: sample x gene matrix X
  # - Drops genes with any NA across samples
  # - Drops zero-variance genes (required when scale.=TRUE; safe even if scale.=FALSE)
  mat_gene_sample <- df_in %>%
    select(Gene, SampleID, Abundance) %>%
    pivot_wider(names_from = SampleID, values_from = Abundance)
  
  if (nrow(mat_gene_sample) == 0) return(NULL)
  
  genes <- mat_gene_sample$Gene
  mat <- as.data.frame(mat_gene_sample[, -1, drop = FALSE])
  rownames(mat) <- genes
  
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)
  
  X <- t(as.matrix(mat))  # sample x gene
  
  col_sd <- apply(X, 2, sd, na.rm = TRUE)
  keep_var <- is.finite(col_sd) & (col_sd > 0)
  X <- X[, keep_var, drop = FALSE]
  
  if (ncol(X) < 10 || nrow(X) < 2) return(NULL)
  X
}

pca_from_matrix <- function(X, center = TRUE, scale_ = TRUE) {
  prcomp(X, center = center, scale. = scale_)
}

pca_var_expl <- function(pca) {
  v <- (pca$sdev^2) / sum(pca$sdev^2)
  list(ve1 = v[1], ve2 = if (length(v) >= 2) v[2] else NA_real_)
}

principal_angles_subspace <- function(A, B) {
  # A, B: gene x k matrices of loadings on shared genes.
  # Restricting to shared genes can break orthonormality; re-orthonormalize via QR before SVD.
  if (nrow(A) != nrow(B)) stop2("principal_angles_subspace: row mismatch.")
  if (ncol(A) != ncol(B)) stop2("principal_angles_subspace: col mismatch.")
  k <- ncol(A)
  
  Qa <- qr.Q(qr(A))
  Qb <- qr.Q(qr(B))
  
  M <- t(Qa) %*% Qb
  s <- svd(M)$d
  
  s <- pmin(1, pmax(0, s))
  ang <- acos(s)
  
  data.frame(
    component = paste0("theta", seq_len(k)),
    cos = s,
    cos2 = s^2,
    angle_rad = ang,
    angle_deg = ang * 180 / pi,
    stringsAsFactors = FALSE
  )
}

# --------------------------- Run analyses for sequential pairs ---------------------------

pair_summaries <- list()

for (i in seq_len(nrow(day_pairs))) {
  d1 <- day_pairs$Day1[i]
  d2 <- day_pairs$Day2[i]
  
  pair_lbl <- sprintf("Pair_Day%02d_vs_Day%02d", d1, d2)
  pair_tag <- sprintf("Day%02d_Day%02d", d1, d2)
  
  message(sprintf("---- %s ----", pair_lbl))
  
  df_pair <- df %>% filter(Day %in% c(d1, d2))
  samp_day <- df_pair %>% distinct(SampleID, Day)
  
  # ========= (A) JOINT PCA (pooled Day1+Day2) =========
  X_joint <- build_sample_by_gene_matrix(df_pair)
  if (is.null(X_joint)) {
    message("Skipping pair: joint PCA matrix could not be built after filtering (too many NA or too few genes).")
    next
  }
  
  pca_joint <- pca_from_matrix(X_joint, center = center, scale_ = scale_)
  ve_joint <- pca_var_expl(pca_joint)
  
  load_joint <- as.data.frame(pca_joint$rotation)
  load_joint$Gene <- rownames(load_joint)
  rownames(load_joint) <- NULL
  
  load_joint <- load_joint %>%
    transmute(
      Gene,
      PC1 = .data$PC1,
      PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_
    )
  
  scores_joint <- as.data.frame(pca_joint$x)
  scores_joint$SampleID <- rownames(scores_joint)
  rownames(scores_joint) <- NULL
  
  scores_joint <- scores_joint %>%
    transmute(
      SampleID,
      PC1 = .data$PC1,
      PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_
    ) %>%
    left_join(samp_day, by = "SampleID")
  
  anchored_joint <- anchor_vector_sign_by_gene(load_joint, anchor_preference)
  load_joint2 <- anchored_joint$loadings
  anchor_used_joint <- anchored_joint$anchor_used
  flipped_joint <- anchored_joint$flipped
  
  if (flipped_joint) {
    scores_joint$PC1 <- -scores_joint$PC1
    scores_joint$PC2 <- -scores_joint$PC2
  }
  
  top_k_plot <- 25
  
  loadings_joint_path <- file.path(output_dir, paste0(
    "A_JointPCA__", pair_lbl, "__GeneLoadings_PC1_PC2.csv"
  ))
  scores_joint_path <- file.path(output_dir, paste0(
    "A_JointPCA__", pair_lbl, "__SampleScores_PC1_PC2_withDayLabel.csv"
  ))
  joint_scores_plot_path <- file.path(output_dir, paste0(
    "A_JointPCA__", pair_lbl, "__Plot_PC1ScoresByDay_BoxJitter.png"
  ))
  joint_load_plot_path <- file.path(output_dir, paste0(
    "A_JointPCA__", pair_lbl, "__Plot_Top", top_k_plot, "GeneLoadingsByAbsPC1_Bar.png"
  ))
  
  readr::write_csv(load_joint2 %>% arrange(desc(abs(PC1))), loadings_joint_path)
  readr::write_csv(scores_joint %>% arrange(Day, SampleID), scores_joint_path)
  
  p_joint_scores <- ggplot(scores_joint, aes(x = factor(Day), y = PC1)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    labs(
      title = paste0("JOINT PCA (pooled) PC1 sample scores: Day ", d1, " vs Day ", d2),
      subtitle = paste0(
        "PC1 variance=", sprintf("%.2f%%", 100 * ve_joint$ve1),
        " | Anchor gene=", anchor_used_joint,
        if (flipped_joint) " | PC1 flipped by anchor rule" else ""
      ),
      x = "Day",
      y = "PC1 score (pooled PCA)"
    )
  write_plot_safe(p_joint_scores, joint_scores_plot_path)
  
  top_joint <- load_joint2 %>%
    arrange(desc(abs(PC1))) %>%
    slice_head(n = top_k_plot) %>%
    mutate(Gene = factor(Gene, levels = rev(Gene)))
  
  p_joint_load <- ggplot(top_joint, aes(x = Gene, y = PC1)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0("JOINT PCA: Top ", top_k_plot, " genes by |PC1 loading| (pooled Day ", d1, "+", d2, ")"),
      subtitle = paste0(
        "Variance explained: PC1=", sprintf("%.2f%%", 100 * ve_joint$ve1),
        ifelse(is.na(ve_joint$ve2), "", paste0(", PC2=", sprintf("%.2f%%", 100 * ve_joint$ve2))),
        " | Anchor gene=", anchor_used_joint
      ),
      x = "Gene",
      y = "PC1 loading (pooled PCA)"
    )
  write_plot_safe(p_joint_load, joint_load_plot_path)
  
  # ========= (D) Directional context: Δmean per gene =========
  mean_tbl <- df_pair %>%
    group_by(Day, Gene) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Day, values_from = mean_abundance)
  
  col_d1 <- as.character(d1)
  col_d2 <- as.character(d2)
  
  if (!(col_d1 %in% names(mean_tbl)) || !(col_d2 %in% names(mean_tbl))) {
    mean_tbl$delta_mean <- NA_real_
  } else {
    mean_tbl$delta_mean <- mean_tbl[[col_d2]] - mean_tbl[[col_d1]]
  }
  
  delta_tbl <- mean_tbl %>% select(Gene, delta_mean)
  
  delta_path <- file.path(output_dir, paste0(
    "D_Direction__", pair_lbl, "__GeneDeltaMean_Day", sprintf("%02d", d2),
    "_minus_Day", sprintf("%02d", d1), ".csv"
  ))
  
  readr::write_csv(delta_tbl %>% arrange(desc(abs(delta_mean))), delta_path)
  
  joint_delta_merge <- load_joint2 %>% inner_join(delta_tbl, by = "Gene")
  cor_joint_loading_delta <- suppressWarnings(stats::cor(
    joint_delta_merge$PC1, joint_delta_merge$delta_mean, use = "pairwise.complete.obs"
  ))
  
  # ========= (B) Day-specific PCA (Day1 and Day2 separately) + PC1 comparison =========
  df_d1 <- df %>% filter(Day == d1)
  df_d2 <- df %>% filter(Day == d2)
  
  X_d1 <- build_sample_by_gene_matrix(df_d1)
  X_d2 <- build_sample_by_gene_matrix(df_d2)
  
  day_pca_ok <- !(is.null(X_d1) || is.null(X_d2))
  
  cos_sim_pc1 <- NA_real_
  angle_deg_pc1 <- NA_real_
  cor_pc1 <- NA_real_
  anchor_used_day <- NA_character_
  flipped_d2_to_d1 <- FALSE
  
  load_d1_path <- NA_character_
  load_d2_path <- NA_character_
  comp_scatter_path <- NA_character_
  
  # ========= (C) Subspace similarity PC1..PCk =========
  subspace_ok <- FALSE
  k_used <- NA_integer_
  mean_cos2 <- NA_real_
  mean_angle_deg <- NA_real_
  max_angle_deg <- NA_real_
  subspace_angles_path <- NA_character_
  
  if (day_pca_ok) {
    pca_d1 <- pca_from_matrix(X_d1, center = center, scale_ = scale_)
    pca_d2 <- pca_from_matrix(X_d2, center = center, scale_ = scale_)
    
    load_d1_full <- as.data.frame(pca_d1$rotation)
    load_d1_full$Gene <- rownames(load_d1_full)
    rownames(load_d1_full) <- NULL
    
    load_d2_full <- as.data.frame(pca_d2$rotation)
    load_d2_full$Gene <- rownames(load_d2_full)
    rownames(load_d2_full) <- NULL
    
    # (B) PC1/PC2 loadings for interpretability and PC1 comparison
    load_d1_pc12 <- load_d1_full %>%
      transmute(Gene, PC1 = .data$PC1, PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_)
    load_d2_pc12 <- load_d2_full %>%
      transmute(Gene, PC1 = .data$PC1, PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_)
    
    # Anchor signs within each day (so “positive pole” is consistent-ish)
    anch1 <- anchor_vector_sign_by_gene(load_d1_pc12, anchor_preference)
    load_d1a <- anch1$loadings
    anchor_used_day <- anch1$anchor_used
    
    anch2 <- anchor_vector_sign_by_gene(load_d2_pc12, anchor_preference)
    load_d2a <- anch2$loadings
    
    load_d1_path <- file.path(output_dir, paste0(
      "B_DayPCA__", pair_lbl, "__Within_Day", sprintf("%02d", d1), "__GeneLoadings_PC1_PC2.csv"
    ))
    load_d2_path <- file.path(output_dir, paste0(
      "B_DayPCA__", pair_lbl, "__Within_Day", sprintf("%02d", d2), "__GeneLoadings_PC1_PC2.csv"
    ))
    
    readr::write_csv(load_d1a %>% arrange(desc(abs(PC1))), load_d1_path)
    readr::write_csv(load_d2a %>% arrange(desc(abs(PC1))), load_d2_path)
    
    # PC1 comparison on shared genes
    shared_pc1 <- intersect(load_d1a$Gene, load_d2a$Gene)
    
    if (length(shared_pc1) >= 50) {
      v1 <- load_d1a$PC1[match(shared_pc1, load_d1a$Gene)]
      v2 <- load_d2a$PC1[match(shared_pc1, load_d2a$Gene)]
      
      aligned <- align_sign_to_reference(v1, v2)
      v2a <- aligned$v
      flipped_d2_to_d1 <- aligned$flipped
      
      cos_sim_pc1 <- cosine_similarity(v1, v2a)
      if (is.finite(cos_sim_pc1)) {
        cos_clip <- max(min(cos_sim_pc1, 1), -1)
        angle_deg_pc1 <- acos(cos_clip) * 180 / pi
      }
      cor_pc1 <- suppressWarnings(stats::cor(v1, v2a, use = "pairwise.complete.obs"))
      
      comp_df <- data.frame(
        Gene = shared_pc1,
        PC1_loading_Day1 = v1,
        PC1_loading_Day2 = v2a,
        stringsAsFactors = FALSE
      )
      
      p_comp <- ggplot(comp_df, aes(x = PC1_loading_Day1, y = PC1_loading_Day2)) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = 0, linewidth = 0.2) +
        geom_vline(xintercept = 0, linewidth = 0.2) +
        labs(
          title = paste0("DAY-SPECIFIC PC1 loading comparison (aligned): Day ", d1, " vs Day ", d2),
          subtitle = paste0(
            "Shared genes=", length(shared_pc1),
            " | cosine=", sprintf("%.3f", cos_sim_pc1),
            " | angle=", sprintf("%.1f°", angle_deg_pc1),
            " | cor=", sprintf("%.3f", cor_pc1),
            " | anchor=", anchor_used_day,
            if (flipped_d2_to_d1) " | Day2 PC1 sign flipped to match Day1" else ""
          ),
          x = paste0("PC1 loading within Day ", d1),
          y = paste0("PC1 loading within Day ", d2)
        )
      
      comp_scatter_path <- file.path(output_dir, paste0(
        "B_DayPCA__", pair_lbl, "__Plot_Day", sprintf("%02d", d1),
        "_PC1vsDay", sprintf("%02d", d2), "_PC1_LoadingScatter_Aligned.png"
      ))
      write_plot_safe(p_comp, comp_scatter_path)
    } else {
      message("NOTE: Too few shared genes for day-specific PC1 comparison in ", pair_lbl, " (shared=", length(shared_pc1), ").")
    }
    
    # (C) Subspace similarity (PC1..PCk) using full rotation matrices
    n_pc_d1 <- sum(grepl("^PC\\d+$", names(load_d1_full)))
    n_pc_d2 <- sum(grepl("^PC\\d+$", names(load_d2_full)))
    k_max <- min(K_subspace, n_pc_d1, n_pc_d2)
    
    if (k_max >= 2) {
      shared_sub <- intersect(load_d1_full$Gene, load_d2_full$Gene)
      
      if (length(shared_sub) >= 200) {
        pc_cols <- paste0("PC", seq_len(k_max))
        if (!all(pc_cols %in% names(load_d1_full)) || !all(pc_cols %in% names(load_d2_full))) {
          message("NOTE: Missing required PC columns for subspace comparison in ", pair_lbl, "; skipping subspace.")
        } else {
          A <- as.matrix(load_d1_full[match(shared_sub, load_d1_full$Gene), pc_cols, drop = FALSE])
          B <- as.matrix(load_d2_full[match(shared_sub, load_d2_full$Gene), pc_cols, drop = FALSE])
          
          ok_rows <- apply(A, 1, function(x) all(is.finite(x))) & apply(B, 1, function(x) all(is.finite(x)))
          A <- A[ok_rows, , drop = FALSE]
          B <- B[ok_rows, , drop = FALSE]
          
          if (nrow(A) >= 200) {
            k_used <- k_max
            ang_tbl <- principal_angles_subspace(A, B)
            subspace_ok <- TRUE
            
            mean_cos2 <- mean(ang_tbl$cos2, na.rm = TRUE)
            mean_angle_deg <- mean(ang_tbl$angle_deg, na.rm = TRUE)
            max_angle_deg <- max(ang_tbl$angle_deg, na.rm = TRUE)
            
            subspace_angles_path <- file.path(output_dir, paste0(
              "C_Subspace__", pair_lbl, "__PrincipalAngles_Subspace_PC1_to_PC", k_used, ".csv"
            ))
            readr::write_csv(ang_tbl, subspace_angles_path)
          } else {
            message("NOTE: Too few finite shared genes for subspace angles in ", pair_lbl, ".")
          }
        }
      } else {
        message("NOTE: Too few shared genes for subspace comparison in ", pair_lbl, " (shared=", length(shared_sub), ").")
      }
    } else {
      message("NOTE: Too few PCs available for subspace comparison in ", pair_lbl, " (k_max=", k_max, ").")
    }
  } else {
    message("NOTE: Day-specific PCA skipped for ", pair_lbl, " (insufficient data after filtering).")
  }
  
  # Summaries: Joint PCA poles + Directional poles
  top_pos_joint <- load_joint2 %>% arrange(desc(PC1)) %>% slice_head(n = 10) %>% pull(Gene)
  top_neg_joint <- load_joint2 %>% arrange(PC1) %>% slice_head(n = 10) %>% pull(Gene)
  
  top_up <- delta_tbl %>% arrange(desc(delta_mean)) %>% slice_head(n = 10) %>% pull(Gene)
  top_down <- delta_tbl %>% arrange(delta_mean) %>% slice_head(n = 10) %>% pull(Gene)
  
  pair_summaries[[pair_lbl]] <- data.frame(
    Pair = pair_lbl,
    Day1 = d1,
    Day2 = d2,
    center = center,
    scale = scale_,
    
    # Joint PCA
    joint_anchor_gene = anchor_used_joint,
    joint_pc1_flipped_by_anchor = flipped_joint,
    joint_var_expl_PC1 = ve_joint$ve1,
    joint_var_expl_PC2 = ve_joint$ve2,
    joint_Top10_Pos_PC1 = paste(top_pos_joint, collapse = "; "),
    joint_Top10_Neg_PC1 = paste(top_neg_joint, collapse = "; "),
    
    # Directional context
    corr_jointPC1_vs_deltaMean = cor_joint_loading_delta,
    deltaMean_Top10_Up = paste(top_up, collapse = "; "),
    deltaMean_Top10_Down = paste(top_down, collapse = "; "),
    
    # Day-specific PC1 comparison
    dayPC1_anchor_gene = anchor_used_day,
    cos_sim_dayPC1 = cos_sim_pc1,
    angle_deg_dayPC1 = angle_deg_pc1,
    cor_dayPC1 = cor_pc1,
    day2_pc1_flipped_to_match_day1 = flipped_d2_to_d1,
    
    # Subspace similarity PC1..PCK
    K_subspace_requested = K_subspace,
    K_subspace_used = k_used,
    subspace_ok = subspace_ok,
    subspace_mean_cos2 = mean_cos2,
    subspace_mean_angle_deg = mean_angle_deg,
    subspace_max_angle_deg = max_angle_deg,
    
    stringsAsFactors = FALSE
  )
}

if (length(pair_summaries) == 0) {
  stop2("No day-pairs produced valid outputs. Check input data for missingness or filtering effects.")
}

pair_summary_tbl <- bind_rows(pair_summaries)

pair_summary_path <- file.path(output_dir, paste0(
  "Z_Summary__AllPairs__Metrics_JointPCA_DayPCA_Subspace_Direction.csv"
))
readr::write_csv(pair_summary_tbl, pair_summary_path)

# --------------------------- Manifest (MANDATORY) ---------------------------

get_pkg_versions <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    data.frame(package = p, version = v, stringsAsFactors = FALSE)
  })
  bind_rows(out)
}

deps_tbl <- get_pkg_versions(required_pkgs)

manifest_path_json <- file.path(output_dir, "Z_Manifest__Project_Manifest.json")
manifest_path_csv  <- file.path(output_dir, "Z_Manifest__FileInventory.csv")

# Script header embedding
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
    "# Sequential Day-Pair PCA Toolkit (Transcriptome)",
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
    day_pairs = day_pairs,
    center = center,
    scale = scale_,
    anchor_preference = anchor_preference,
    K_subspace_requested = K_subspace,
    analyses = list(
      joint_pca = list(purpose = "Between-day contrast axis on pooled samples (Day1+Day2)"),
      day_specific_pca = list(purpose = "Within-day organization axes and PC1 similarity across days"),
      subspace_similarity = list(purpose = "PC1..PCk subspace similarity (principal angles) for stability under PC rotation/reordering"),
      directional_context = list(purpose = "Δmean (Day2-Day1) adds temporal direction; PCA itself is non-directional")
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

writeLines(
  jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"),
  manifest_path_json
)

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

# --------------------------- Human README mapping files -> meaning ---------------------------

readme_path <- file.path(output_dir, "Z_Manifest__README_FileDescriptions.csv")

describe_file <- function(f) {
  if (grepl("^A_JointPCA__", f) && grepl("GeneLoadings", f)) {
    return("JOINT PCA (pooled Day1+Day2): gene loadings for PC1/PC2. Rows are genes; columns are PC loadings. Use to interpret the between-day contrast axis.")
  }
  if (grepl("^A_JointPCA__", f) && grepl("SampleScores", f)) {
    return("JOINT PCA (pooled Day1+Day2): sample scores for PC1/PC2 with Day label. Rows are samples; use PC1 to see separation between days.")
  }
  if (grepl("^A_JointPCA__", f) && grepl("Plot_PC1ScoresByDay", f)) {
    return("JOINT PCA plot: distribution of pooled-PC1 scores by day (box+jitter). Shows whether pooled PC1 separates Day1 and Day2.")
  }
  if (grepl("^A_JointPCA__", f) && grepl("Top.*GeneLoadingsByAbsPC1", f)) {
    return("JOINT PCA plot: top genes by absolute PC1 loading (pooled). Identifies strongest contributors to between-day contrast axis.")
  }
  if (grepl("^B_DayPCA__", f) && grepl("Within_.*GeneLoadings", f)) {
    return("DAY-SPECIFIC PCA: gene loadings for PC1/PC2 computed within a single day. Use to understand within-day coordination structure.")
  }
  if (grepl("^B_DayPCA__", f) && grepl("LoadingScatter_Aligned", f)) {
    return("DAY-SPECIFIC PC1 comparison plot: scatter of Day1 PC1 loadings vs Day2 PC1 loadings (sign-aligned). Tight diagonal = conserved axis; diffuse/rotated = reorganization.")
  }
  if (grepl("^C_Subspace__", f) && grepl("PrincipalAngles", f)) {
    return("SUBSPACE similarity table: principal angles between Day1 and Day2 PC1..PCk subspaces. Smaller angles / higher cos^2 = conserved dominant manifold.")
  }
  if (grepl("^D_Direction__", f) && grepl("GeneDeltaMean", f)) {
    return("DIRECTIONAL context table: gene-level Δmean = mean(Day2) - mean(Day1). Adds temporal direction not provided by PCA.")
  }
  if (grepl("^Z_Summary__", f)) {
    return("RUN SUMMARY table: one row per day-pair, includes variance explained, joint poles, Δmean correlation, day-PC1 similarity, and subspace similarity metrics.")
  }
  if (grepl("^Z_Manifest__", f) && grepl("Project_Manifest\\.json$", f)) {
    return("PROJECT MANIFEST (JSON): machine-readable run metadata (run_id, inputs, parameters, dependencies, outputs).")
  }
  if (grepl("^Z_Manifest__", f) && grepl("FileInventory\\.csv$", f)) {
    return("FILE INVENTORY (CSV): all files generated in output directory with sizes and timestamps.")
  }
  if (grepl("^Z_Manifest__", f) && grepl("README_FileDescriptions\\.csv$", f)) {
    return("README mapping file: links each output filename to a plain-English description.")
  }
  if (grepl("_Report\\.qmd$", f)) {
    return("Quarto source report (QMD): inputs, logic, plots, outputs inventory, and reproducibility section.")
  }
  if (grepl("_Report\\.html$", f)) {
    return("Rendered HTML report: human-readable summary of run, including plots, tables, and reproducibility info.")
  }
  "Other output file."
}

file_desc <- data.frame(
  file = inv0$file,
  description = vapply(inv0$file, describe_file, character(1)),
  stringsAsFactors = FALSE
)
readr::write_csv(file_desc, readme_path)

# --------------------------- Refresh inventory before report ---------------------------

inv_pre_report <- write_inventory_csv(output_dir, manifest_path_csv)

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
  "manifest <- jsonlite::fromJSON('", basename(manifest_path_json), "', simplifyVector = FALSE)\n",
  "inv <- read.csv('", basename(manifest_path_csv), "', stringsAsFactors = FALSE)\n",
  "readme <- read.csv('Z_Manifest__README_FileDescriptions.csv', stringsAsFactors = FALSE)\n",
  "meta_tbl <- data.frame(\n",
  "  Field = c('run_id','run_timestamp','script_name','script_path','script_full_path','outputs_root','output_dir','input_csv','center','scale','K_subspace'),\n",
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
  "    as.character(manifest[['parameters']][['K_subspace_requested']])\n",
  "  ),\n",
  "  stringsAsFactors = FALSE\n",
  ")\n",
  "knitr::kable(meta_tbl)\n",
  "```\n\n",
  
  "# Script header\n\n",
  "```r\n",
  gsub("```", "``\\`", header_text),
  "\n```\n\n",
  
  "# What each analysis means (and does not mean)\n\n",
  "## (A) Joint PCA (pooled Day1 + Day2): a **between-day contrast axis**\n\n",
  "- **Useful for:** identifying coordinated gene programs that most strongly distinguish Day1 from Day2.\n",
  "- **Not useful for:** comparing “PC1 of Day1” to “PC1 of Day2.” Joint PCA produces a single axis for the pooled samples.\n",
  "- **Directionality note:** PCA is symmetric; temporal direction is not intrinsic and must be added using time ordering and mean changes.\n\n",
  
  "## (B) Day-specific PCA: a **within-day organization axis**\n\n",
  "- **Useful for:** quantifying whether the dominant within-day coordination axis (PC1) is preserved or reorganized between Day1 and Day2.\n",
  "- **Quantified by:** cosine similarity, angle, and correlation between Day1 and Day2 PC1 loadings on shared genes.\n",
  "- **Important caveat:** PC1 alone can be unstable if the top eigenvalues are close. That is why we also compute subspace similarity.\n\n",
  
  "## (C) Subspace similarity (PC1..PCk): a **dominant coordination manifold**\n\n",
  "- **Useful for:** detecting continuity of the dominant low-dimensional structure when individual PCs reorder or rotate.\n",
  "- **How it works:** compare the K-dimensional subspaces spanned by PC1..PCk for Day1 and Day2 using **principal angles**.\n\n",
  "If \\(Q_1\\) and \\(Q_2\\) are orthonormal bases for the two subspaces, compute singular values \\(\\sigma_i\\) of \\(Q_1^T Q_2\\).\n",
  "Principal angles are:\n",
  "\\[\n",
  "\\theta_i = \\arccos(\\sigma_i)\n",
  "\\]\n",
  "We report \\(\\cos^2(\\theta_i)=\\sigma_i^2\\) (similarity) and summarize by mean and max angles.\n\n",
  
  "## (D) Directional context (Δmean)\n\n",
  "For each gene, \\(\\Delta\\text{mean} = \\text{mean}_{Day2} - \\text{mean}_{Day1}\\). This adds temporal context to the otherwise non-directional PCA contrast.\n\n",
  
  "# Dependencies\n\n",
  "```{r}\n",
  "deps <- manifest[['dependencies']]\n",
  "deps_tbl <- do.call(rbind, lapply(deps, function(x) data.frame(package=x$package, version=x$version, stringsAsFactors=FALSE)))\n",
  "knitr::kable(deps_tbl)\n",
  "```\n\n",
  
  "# Outputs: file inventory (machine) + human README\n\n",
  "## File inventory\n\n",
  "```{r}\n",
  "knitr::kable(inv)\n",
  "```\n\n",
  
  "## Human-readable file descriptions\n\n",
  "```{r}\n",
  "knitr::kable(readme)\n",
  "```\n\n",
  
  "# Key plot (auto-selected)\n\n",
  "```{r}\n",
  "pngs <- inv$file[grepl('\\\\.png$', inv$file)]\n",
  "# Prefer a day-specific PC1 scatter plot; if absent, prefer joint PC1 scores plot; else any PNG.\n",
  "pref1 <- pngs[grepl('^B_DayPCA__.*LoadingScatter_Aligned\\\\.png$', pngs)]\n",
  "pref2 <- pngs[grepl('^A_JointPCA__.*Plot_PC1ScoresByDay_BoxJitter\\\\.png$', pngs)]\n",
  "if (length(pref1) > 0) {\n",
  "  knitr::include_graphics(pref1[1])\n",
  "} else if (length(pref2) > 0) {\n",
  "  knitr::include_graphics(pref2[1])\n",
  "} else if (length(pngs) > 0) {\n",
  "  knitr::include_graphics(pngs[1])\n",
  "} else {\n",
  "  cat('No PNG plots found to display.')\n",
  "}\n",
  "```\n\n",
  
  "# Interpretation guide (how to read the outputs)\n\n",
  "- **A_JointPCA GeneLoadings:** genes defining the dominant pooled contrast between the two days.\n",
  "- **A_JointPCA SampleScores:** sample positions along the pooled contrast axis; separation by day indicates a coordinated shift.\n",
  "- **B_DayPCA Within_Day Loadings:** within-day coordination axes (organization structure) for each day.\n",
  "- **B_DayPCA LoadingScatter (aligned):** continuity vs reorganization of Day1-PC1 vs Day2-PC1.\n",
  "- **C_Subspace PrincipalAngles:** stability of the dominant PC1..PCk manifold under PC swapping/rotation.\n",
  "- **D_Direction Δmean:** temporal direction of mean expression change; not the same as coordination structure.\n\n",
  
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

# --------------------------- Refresh inventory AFTER rendering (MANDATORY) ---------------------------

inv_post_report <- write_inventory_csv(output_dir, manifest_path_csv)

# Also refresh README after render (so the HTML/QMD files are described)
file_desc2 <- data.frame(
  file = inv_post_report$file,
  description = vapply(inv_post_report$file, describe_file, character(1)),
  stringsAsFactors = FALSE
)
readr::write_csv(file_desc2, readme_path)

# Final inventory refresh (README updated)
inv_final <- write_inventory_csv(output_dir, manifest_path_csv)

# Console confirmation
if (render_ok && file.exists(html_path)) {
  message("HTML report created: ", normalizePath(html_path, winslash = "/", mustWork = FALSE))
} else {
  message("ERROR: HTML report render failed.")
  if (!is.null(render_msg)) message("Reason: ", render_msg)
  message("QMD written at: ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE))
  message("Check outputs in: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
}

message("Run complete. Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
