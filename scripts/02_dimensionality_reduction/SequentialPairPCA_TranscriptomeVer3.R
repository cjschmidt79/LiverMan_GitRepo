#!/usr/bin/env Rscript
# ======================================================================
# Sequential Day-Pair PCA Toolkit (Transcriptome): Joint PCA + Day-Specific PCA + Subspace Similarity + Directional Context
#
# UPDATE (2026-01-14):
#   - Adds user-selectable transform:
#       (0) raw
#       (1) log2(Abundance + pseudocount)
#   - Transform is applied consistently to:
#       • PCA matrices (joint and day-specific)
#       • Δmean calculations
#       • Outputs/manifest/QMD report annotations
#
# INPUT FORMAT
#   One tidy-long transcriptome table with columns:
#     Day, SampleID, Gene, Abundance
#
# OUTPUTS (STRICT)
#   All outputs are written ONLY under:
#     outputs/<run_id>/
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
  "jsonlite", "knitr",
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

# --------------------------- Argument parsing ---------------------------

args <- commandArgs(trailingOnly = FALSE)

get_arg_value <- function(prefix) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", prefix, "="), "", hit[1])
}

# Inputs / knobs
input_csv <- get_arg_value("--input")
center_flag <- get_arg_value("--center")
scale_flag  <- get_arg_value("--scale")
anchor_gene_user <- get_arg_value("--anchor")
k_flag <- get_arg_value("--k")

# NEW: transform controls
transform_flag <- get_arg_value("--transform")     # raw | log2p
pseudocount_flag <- get_arg_value("--pseudocount") # numeric

center <- if (!is.na(center_flag)) tolower(center_flag) %in% c("true", "t", "1", "yes", "y") else TRUE
scale_ <- if (!is.na(scale_flag))  tolower(scale_flag)  %in% c("true", "t", "1", "yes", "y") else TRUE

K_subspace <- 5L
if (!is.na(k_flag) && nzchar(k_flag)) {
  K_subspace <- suppressWarnings(as.integer(k_flag))
  if (!is.finite(K_subspace) || is.na(K_subspace) || K_subspace < 1) {
    stop2("Invalid --k value: ", k_flag, " (must be integer >= 1)")
  }
}

# Input file
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

# --------------------------- Transform choice (NEW) ---------------------------

normalize_transform_choice <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  y <- tolower(trimws(x))
  if (y %in% c("0", "raw")) return("raw")
  if (y %in% c("1", "log2p", "log2+pseudocount", "log2pseudocount", "log2")) return("log2p")
  NA_character_
}

transform_mode <- normalize_transform_choice(transform_flag)

pseudocount <- 1
if (!is.na(pseudocount_flag) && nzchar(pseudocount_flag)) {
  pc <- suppressWarnings(as.numeric(pseudocount_flag))
  if (!is.finite(pc) || is.na(pc) || pc < 0) stop2("Invalid --pseudocount value: ", pseudocount_flag, " (must be numeric >= 0)")
  pseudocount <- pc
}

# If transform not specified by CLI, prompt interactively.
# (We do not ask extra questions beyond the minimal choice.)
if (is.na(transform_mode)) {
  # If not interactive (e.g., running under strict non-interactive batch), default to raw.
  is_interactive_like <- interactive() || (!is.null(Sys.getenv("RSTUDIO", unset = NA)) && nzchar(Sys.getenv("RSTUDIO")))
  if (is_interactive_like) {
    message("---- Transform Choice ----")
    message("Choose data scale for analysis:")
    message("  (0) raw Abundance")
    message("  (1) log2(Abundance + pseudocount)")
    ans <- readline(prompt = "Enter 0 or 1 [default 0]: ")
    ans2 <- trimws(ans)
    if (!nzchar(ans2)) ans2 <- "0"
    transform_mode <- normalize_transform_choice(ans2)
    if (is.na(transform_mode)) {
      message("Invalid choice; defaulting to raw (0).")
      transform_mode <- "raw"
    }
  } else {
    message("NOTE: Non-interactive run and --transform not provided; defaulting to raw.")
    transform_mode <- "raw"
  }
}

transform_label <- if (transform_mode == "raw") "raw" else paste0("log2p_pc", pseudocount)

message("---- Run Parameters ----")
message("input_csv: ", input_csv)
message("center: ", center)
message("scale: ", scale_)
message("K_subspace (requested): ", K_subspace)
message("anchor_preference: ", paste(anchor_preference, collapse = " -> "))
message("transform_mode: ", transform_mode)
message("pseudocount: ", pseudocount)
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

# Apply transform (NEW): Includes CPM Normalization and Log2 VST-Lite
apply_transform_cpm <- function(df_in, mode, pc) {
  if (mode == "raw") {
    return(df_in %>% mutate(Abundance_use = Abundance))
  }
  df_in %>%
    group_by(SampleID) %>%
    mutate(lib_size = sum(Abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      CPM = (Abundance / lib_size) * 1e6,
      Abundance_use = log2(CPM + pc)
    ) %>%
    select(-lib_size, -CPM)
}

df <- apply_transform_cpm(df, transform_mode, pseudocount)

if (transform_mode == "log2p") {
  bad_n <- sum(!is.finite(df$Abundance_use) & !is.na(df$Abundance_use))
  if (bad_n > 0) {
    message("NOTE: ", bad_n, " values became non-finite after log2 transform. These will be imputed to the pseudocount floor.")
  }
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
  ggplot2::ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi, units = "in")
  if (!file.exists(path)) stop2("Failed to write plot: ", path)
}
write_plot_rds_safe <- function(p, path) {
  saveRDS(p, file = path)
  if (!file.exists(path)) stop2("Failed to write plot RDS: ", path)
}

build_sample_by_gene_matrix <- function(df_in, mode_in, pc_val) {
  # Returns: sample x gene matrix X (using Abundance_use)
  mat_gene_sample <- df_in %>%
    select(Gene, SampleID, Abundance_use) %>%
    pivot_wider(names_from = SampleID, values_from = Abundance_use)
  
  if (nrow(mat_gene_sample) == 0) return(NULL)
  
  genes <- mat_gene_sample$Gene
  mat <- as.matrix(mat_gene_sample[, -1, drop = FALSE])
  rownames(mat) <- genes
  
  # Impute missing/non-finite with the scale floor (log2 pseudocount or 0)
  fill_val <- if(mode_in == "log2p") log2(pc_val) else 0
  mat[!is.finite(mat) | is.na(mat)] <- fill_val
  
  X <- t(mat) # sample x gene
  
  if (ncol(X) < 10 || nrow(X) < 2) return(NULL)
  
  col_sd <- apply(X, 2, sd, na.rm = TRUE)
  keep_var <- is.finite(col_sd) & (col_sd > 0)
  X <- X[, keep_var, drop = FALSE]
  
  if (ncol(X) < 10 || nrow(X) < 2) return(NULL)
  X
}
pca_from_matrix <- function(X, center = TRUE, scale_ = TRUE) {
  stats::prcomp(X, center = center, scale. = scale_)
}

pca_var_expl <- function(pca) {
  v <- (pca$sdev^2) / sum(pca$sdev^2)
  list(ve1 = v[1], ve2 = if (length(v) >= 2) v[2] else NA_real_)
}

principal_angles_subspace <- function(A, B) {
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
  
  message(sprintf("---- %s ----", pair_lbl))
  
  df_pair <- df %>% filter(Day %in% c(d1, d2))
  samp_day <- df_pair %>% distinct(SampleID, Day)
  
  # ========= (A) JOINT PCA (pooled Day1+Day2) =========
  # ========= (A) JOINT PCA (pooled Day1+Day2) =========
  X_joint <- build_sample_by_gene_matrix(df_pair, transform_mode, pseudocount)
  if (is.null(X_joint)) {
    message("Skipping pair: joint PCA matrix could not be built.")
    next
  }
  
  pca_joint <- pca_from_matrix(X_joint, center = center, scale_ = scale_)
  ve_joint <- pca_var_expl(pca_joint)
  
  # SCREE PLOT
  var_df <- data.frame(
    PC = seq_along(pca_joint$sdev),
    VarExpl = (pca_joint$sdev^2) / sum(pca_joint$sdev^2)
  )
  p_scree <- ggplot(var_df[1:min(10, nrow(var_df)), ], aes(x = factor(PC), y = VarExpl)) +
    geom_col(fill = "steelblue") +
    geom_line(aes(x = PC), group = 1, color = "red") +
    geom_point(aes(x = PC), color = "red") +
    labs(title = paste("Scree Plot:", pair_lbl), x = "Principal Component", y = "Prop. Variance Explained") +
    theme_minimal()
  write_plot_safe(p_scree, file.path(output_dir, paste0("A_ScreePlot__", pair_lbl, ".png")))
  
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
    "A_JointPCA__", pair_lbl, "__", transform_label, "__GeneLoadings_PC1_PC2.csv"
  ))
  scores_joint_path <- file.path(output_dir, paste0(
    "A_JointPCA__", pair_lbl, "__", transform_label, "__SampleScores_PC1_PC2_withDayLabel.csv"
  ))
  joint_scores_plot_path <- file.path(output_dir, paste0(
    "A_JointPCA__", pair_lbl, "__", transform_label, "__Plot_PC1ScoresByDay_BoxJitter.png"
  ))
  joint_load_plot_path <- file.path(output_dir, paste0(
    "A_JointPCA__", pair_lbl, "__", transform_label, "__Plot_Top", top_k_plot, "GeneLoadingsByAbsPC1_Bar.png"
  ))
  
  readr::write_csv(load_joint2 %>% arrange(desc(abs(PC1))), loadings_joint_path)
  readr::write_csv(scores_joint %>% arrange(Day, SampleID), scores_joint_path)
  
  p_joint_scores <- ggplot(scores_joint, aes(x = factor(Day), y = PC1)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    labs(
      title = paste0("JOINT PCA (pooled) PC1 sample scores: Day ", d1, " vs Day ", d2),
      subtitle = paste0(
        "Transform=", transform_label,
        " | PC1 variance=", sprintf("%.2f%%", 100 * ve_joint$ve1),
        " | Anchor gene=", anchor_used_joint,
        if (flipped_joint) " | PC1 flipped by anchor rule" else ""
      ),
      x = "Day",
      y = "PC1 score (pooled PCA)"
    )
  write_plot_safe(p_joint_scores, joint_scores_plot_path)
  joint_scores_plot_rds_path <- sub("\\.png$", ".rds", joint_scores_plot_path)
  write_plot_rds_safe(p_joint_scores, joint_scores_plot_rds_path)
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
        "Transform=", transform_label,
        " | PC1=", sprintf("%.2f%%", 100 * ve_joint$ve1),
        ifelse(is.na(ve_joint$ve2), "", paste0(", PC2=", sprintf("%.2f%%", 100 * ve_joint$ve2))),
        " | Anchor gene=", anchor_used_joint
      ),
      x = "Gene",
      y = "PC1 loading (pooled PCA)"
    )
  write_plot_safe(p_joint_load, joint_load_plot_path)
  joint_load_plot_rds_path <- sub("\\.png$", ".rds", joint_load_plot_path)
  write_plot_rds_safe(p_joint_load, joint_load_plot_rds_path)
  
  
  # ========= (D) Directional context: Δmean per gene (USING Abundance_use) =========
  mean_tbl <- df_pair %>%
    group_by(Day, Gene) %>%
    summarize(mean_abundance = mean(Abundance_use, na.rm = TRUE), .groups = "drop") %>%
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
    "D_Direction__", pair_lbl, "__", transform_label,
    "__GeneDeltaMean_Day", sprintf("%02d", d2),
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
  
  X_d1 <- build_sample_by_gene_matrix(df_d1, transform_mode, pseudocount)
  X_d2 <- build_sample_by_gene_matrix(df_d2, transform_mode, pseudocount)
  day_pca_ok <- !(is.null(X_d1) || is.null(X_d2))
  
  cos_sim_pc1 <- NA_real_
  angle_deg_pc1 <- NA_real_
  cor_pc1 <- NA_real_
  anchor_used_day <- NA_character_
  flipped_d2_to_d1 <- FALSE
  
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
    
    load_d1_pc12 <- load_d1_full %>%
      transmute(Gene, PC1 = .data$PC1, PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_)
    load_d2_pc12 <- load_d2_full %>%
      transmute(Gene, PC1 = .data$PC1, PC2 = if ("PC2" %in% names(.)) .data$PC2 else NA_real_)
    
    anch1 <- anchor_vector_sign_by_gene(load_d1_pc12, anchor_preference)
    load_d1a <- anch1$loadings
    anchor_used_day <- anch1$anchor_used
    
    anch2 <- anchor_vector_sign_by_gene(load_d2_pc12, anchor_preference)
    load_d2a <- anch2$loadings
    
    load_d1_path <- file.path(output_dir, paste0(
      "B_DayPCA__", pair_lbl, "__", transform_label, "__Within_Day", sprintf("%02d", d1), "__GeneLoadings_PC1_PC2.csv"
    ))
    load_d2_path <- file.path(output_dir, paste0(
      "B_DayPCA__", pair_lbl, "__", transform_label, "__Within_Day", sprintf("%02d", d2), "__GeneLoadings_PC1_PC2.csv"
    ))
    
    readr::write_csv(load_d1a %>% arrange(desc(abs(PC1))), load_d1_path)
    readr::write_csv(load_d2a %>% arrange(desc(abs(PC1))), load_d2_path)
    
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
            "Transform=", transform_label,
            " | Shared genes=", length(shared_pc1),
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
        "B_DayPCA__", pair_lbl, "__", transform_label, "__Plot_Day", sprintf("%02d", d1),
        "_PC1vsDay", sprintf("%02d", d2), "_PC1_LoadingScatter_Aligned.png"
      ))
      write_plot_safe(p_comp, comp_scatter_path)
      comp_scatter_rds_path <- sub("\\.png$", ".rds", comp_scatter_path)
      write_plot_rds_safe(p_comp, comp_scatter_rds_path)
    } else {
      message("NOTE: Too few shared genes for day-specific PC1 comparison in ", pair_lbl, " (shared=", length(shared_pc1), ").")
    }
    
    # (C) Subspace similarity (PC1..PCk) using full rotation matrices
    n_pc_d1 <- sum(grepl("^PC\\d+$", names(load_d1_full)))
    n_pc_d2 <- sum(grepl("^PC\\d+$", names(load_d2_full)))
    k_max <- min(K_subspace, n_pc_d1, n_pc_d2)
    
    if (k_max >= 2) {
      shared_sub <- intersect(load_d1_full$Gene, load_d2_full$Gene)
      
      if (length(shared_sub) >= 100) {
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
              "C_Subspace__", pair_lbl, "__", transform_label, "__PrincipalAngles_Subspace_PC1_to_PC", k_used, ".csv"
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
  
  # ---- Top-N lists for downstream Jaccard (rank-threshold sensitivity control) ----
  topN_small <- 10
  topN_large <- 25
  
  top_pos_joint_10 <- load_joint2 %>% arrange(desc(PC1)) %>% slice_head(n = topN_small) %>% pull(Gene)
  top_neg_joint_10 <- load_joint2 %>% arrange(PC1)      %>% slice_head(n = topN_small) %>% pull(Gene)
  
  top_pos_joint_25 <- load_joint2 %>% arrange(desc(PC1)) %>% slice_head(n = topN_large) %>% pull(Gene)
  top_neg_joint_25 <- load_joint2 %>% arrange(PC1)      %>% slice_head(n = topN_large) %>% pull(Gene)
  
  
  top_up <- delta_tbl %>% arrange(desc(delta_mean)) %>% slice_head(n = 10) %>% pull(Gene)
  top_down <- delta_tbl %>% arrange(delta_mean) %>% slice_head(n = 10) %>% pull(Gene)
  
  pair_summaries[[pair_lbl]] <- data.frame(
    Pair = pair_lbl,
    Day1 = d1,
    Day2 = d2,
    center = center,
    scale = scale_,
    Transform = transform_mode,
    Pseudocount = pseudocount,
    
    joint_anchor_gene = anchor_used_joint,
    joint_pc1_flipped_by_anchor = flipped_joint,
    joint_var_expl_PC1 = ve_joint$ve1,
    joint_var_expl_PC2 = ve_joint$ve2,
    joint_Top10_Pos_PC1 = paste(top_pos_joint_10, collapse = "; "),
    joint_Top10_Neg_PC1 = paste(top_neg_joint_10, collapse = "; "),
    joint_Top25_Up_PC1 = paste(top_pos_joint_25, collapse = "; "),
    joint_Top25_Down_PC1 = paste(top_neg_joint_25, collapse = "; "),
    
    
    corr_jointPC1_vs_deltaMean = cor_joint_loading_delta,
    deltaMean_Top10_Up = paste(top_up, collapse = "; "),
    deltaMean_Top10_Down = paste(top_down, collapse = "; "),
    
    dayPC1_anchor_gene = anchor_used_day,
    cos_sim_dayPC1 = cos_sim_pc1,
    angle_deg_dayPC1 = angle_deg_pc1,
    cor_dayPC1 = cor_pc1,
    day2_pc1_flipped_to_match_day1 = flipped_d2_to_d1,
    
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

# ==========================================================
# JACCARD INPUT FILES (WIDE = required for decomposition; LONG = optional)
# ==========================================================

# ---- (1) WIDE Jaccard input (THIS is what the decomposition script expects) ----
jaccard_wide <- pair_summary_tbl %>%
  transmute(
    Day1,
    Day2,
    joint_Top25_Up_PC1   = ifelse(is.na(joint_Top25_Up_PC1),   "", joint_Top25_Up_PC1),
    joint_Top25_Down_PC1 = ifelse(is.na(joint_Top25_Down_PC1), "", joint_Top25_Down_PC1)
  )

jaccard_wide_path <- file.path(output_dir, paste0(
  "JaccardInput__JointPCA_Top25_PC1__", transform_label, "__Wide.csv"
))
readr::write_csv(jaccard_wide, jaccard_wide_path)
message("Wrote Jaccard WIDE input: ", normalizePath(jaccard_wide_path, winslash = "/", mustWork = FALSE))

# ---- (2) LONG (optional; useful for debugging / plotting your own summaries) ----
split_gene_list <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  trimws(unlist(strsplit(x, ";", fixed = TRUE)))
}

make_long_set <- function(df, set_col, set_type, pc = "PC1", rank_n = 25L) {
  out <- lapply(seq_len(nrow(df)), function(i) {
    x <- df[[set_col]][i]
    genes <- split_gene_list(x)
    if (length(genes) == 0) return(NULL)
    data.frame(
      Pair = df$Pair[i],
      Day1 = df$Day1[i],
      Day2 = df$Day2[i],
      set_type = set_type,       # "UP" or "DOWN"
      pc = pc,
      rank_n = rank_n,
      rank = seq_along(genes),
      Gene = genes,
      stringsAsFactors = FALSE
    )
  })
  if (all(vapply(out, is.null, logical(1)))) return(data.frame())
  do.call(rbind, out)
}

long_up25   <- make_long_set(pair_summary_tbl, set_col = "joint_Top25_Up_PC1",   set_type = "UP",   rank_n = 25L)
long_down25 <- make_long_set(pair_summary_tbl, set_col = "joint_Top25_Down_PC1", set_type = "DOWN", rank_n = 25L)

jaccard_long <- rbind(long_up25, long_down25)

jaccard_long_path <- file.path(output_dir, paste0(
  "JaccardReady__JointPCA_Top25_PC1__", transform_label, "__Long.csv"
))
readr::write_csv(jaccard_long, jaccard_long_path)
message("Wrote Jaccard LONG table (optional): ", normalizePath(jaccard_long_path, winslash = "/", mustWork = FALSE))


message("Wrote Jaccard-ready long table: ", normalizePath(jaccard_long_path, winslash = "/", mustWork = FALSE))



pair_summary_path <- file.path(output_dir, paste0(
  "Z_Summary__AllPairs__", transform_label, "__Metrics_JointPCA_DayPCA_Subspace_Direction.csv"
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

# Script header extraction (robust)
header_text <- NULL
if (!is.na(script_full) && file.exists(script_full)) {
  lines <- readLines(script_full, warn = FALSE)
  keep <- c()
  started <- FALSE
  for (ln in lines[seq_len(min(250, length(lines)))]) {
    if (!started) {
      if (grepl("^\\s*#!", ln) || grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
        keep <- c(keep, ln)
        started <- TRUE
      } else {
        break
      }
    } else {
      if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
        keep <- c(keep, ln)
      } else {
        break
      }
    }
  }
  if (length(keep) > 0) header_text <- paste(keep, collapse = "\n")
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
    transform = list(
      mode = transform_mode,
      pseudocount = pseudocount,
      transform_label = transform_label,
      applied_to = c("PCA matrices", "Δmean calculations")
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
    return("JOINT PCA (pooled Day1+Day2): gene loadings for PC1/PC2 on Abundance_use (raw or log2+pseudocount).")
  }
  if (grepl("^A_JointPCA__", f) && grepl("SampleScores", f)) {
    return("JOINT PCA (pooled Day1+Day2): sample scores for PC1/PC2 with Day label (Abundance_use).")
  }
  if (grepl("^A_JointPCA__", f) && grepl("Plot_PC1ScoresByDay", f)) {
    return("JOINT PCA plot: pooled-PC1 scores by day (Abundance_use).")
  }
  if (grepl("^A_JointPCA__", f) && grepl("GeneLoadingsByAbsPC1", f)) {
    return("JOINT PCA plot: top genes by |PC1 loading| (Abundance_use).")
  }
  if (grepl("^B_DayPCA__", f) && grepl("Within_.*GeneLoadings", f)) {
    return("DAY-SPECIFIC PCA: within-day gene loadings for PC1/PC2 (Abundance_use).")
  }
  if (grepl("^B_DayPCA__", f) && grepl("LoadingScatter_Aligned", f)) {
    return("DAY-SPECIFIC PC1 comparison: Day1 vs Day2 PC1 loadings (sign-aligned), using Abundance_use.")
  }
  if (grepl("^C_Subspace__", f) && grepl("PrincipalAngles", f)) {
    return("SUBSPACE similarity: principal angles between Day1 and Day2 PC1..PCk subspaces (Abundance_use).")
  }
  if (grepl("^D_Direction__", f) && grepl("GeneDeltaMean", f)) {
    return("DIRECTIONAL context: gene-level Δmean = mean(Day2) - mean(Day1) computed on Abundance_use.")
  }
  if (grepl("^Z_Summary__", f)) {
    return("RUN SUMMARY: per-pair metrics including transform choice.")
  }
  if (grepl("^Z_Manifest__", f) && grepl("Project_Manifest\\.json$", f)) {
    return("PROJECT MANIFEST (JSON): run metadata (inputs, parameters, transform, dependencies, outputs).")
  }
  if (grepl("^Z_Manifest__", f) && grepl("FileInventory\\.csv$", f)) {
    return("FILE INVENTORY (CSV): all files generated with sizes and timestamps.")
  }
  if (grepl("^Z_Manifest__", f) && grepl("README_FileDescriptions\\.csv$", f)) {
    return("README mapping: output filename to description.")
  }
  if (grepl("_Report\\.qmd$", f)) {
    return("Quarto source report (QMD).")
  }
  if (grepl("_Report\\.html$", f)) {
    return("Rendered HTML report (Quarto).")
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

# Safer header embedding: print as verbatim via cat() inside an R chunk
header_lines_escaped <- gsub("\\\\", "\\\\\\\\", header_text)
header_lines_escaped <- gsub("\"", "\\\\\"", header_lines_escaped)

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
  "  Field = c('run_id','run_timestamp','script_name','script_path','script_full_path','outputs_root','output_dir','input_csv','center','scale','K_subspace','transform_mode','pseudocount','transform_label'),\n",
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
  "    as.character(manifest[['parameters']][['K_subspace_requested']]),\n",
  "    manifest[['parameters']][['transform']][['mode']],\n",
  "    as.character(manifest[['parameters']][['transform']][['pseudocount']]),\n",
  "    manifest[['parameters']][['transform']][['transform_label']]\n",
  "  ),\n",
  "  stringsAsFactors = FALSE\n",
  ")\n",
  "knitr::kable(meta_tbl)\n",
  "```\n\n",
  
  "# Transform definition\n\n",
  "The analysis used **Abundance_use**, defined as:\n\n",
  "- If `transform_mode = raw`:  \\(Abundance\\_use = Abundance\\)\n",
  "- If `transform_mode = log2p`: \\(Abundance\\_use = \\log_2(Abundance + pseudocount)\\)\n\n",
  
  "# Script header\n\n",
  "```{r}\n",
  "cat(\"", gsub("\n", "\\\\n", header_lines_escaped), "\")\n",
  "```\n\n",
  
  "# Dependencies\n\n",
  "```{r}\n",
  "deps <- manifest[['dependencies']]\n",
  "deps_tbl <- do.call(rbind, lapply(deps, function(x) data.frame(package=x$package, version=x$version, stringsAsFactors=FALSE)))\n",
  "knitr::kable(deps_tbl)\n",
  "```\n\n",
  
  "# Outputs\n\n",
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
  
  "# All generated plots (static PNG)\n\n",
  "```{r}\n",
  "pngs <- inv$file[grepl('\\\\.png$', inv$file)]\n",
  "if (length(pngs) > 0) {\n",
  "  for (p in pngs) {\n",
  "    knitr::include_graphics(p)\n",
  "    cat('\\n\\n')\n",
  "  }\n",
  "} else {\n",
  "  cat('No PNG plots found to display.')\n",
  "}\n",
  "```\n\n",
  
  "# Interactive plots (plotly; all RDS plot objects)\n\n",
  "```{r}\n",
  "ok_plotly <- requireNamespace('plotly', quietly = TRUE)\n",
  "ok_htmltools <- requireNamespace('htmltools', quietly = TRUE)\n",
  "if (!ok_plotly) {\n",
  "  cat('plotly is not available in this R session; interactive plots cannot be rendered.')\n",
  "} else if (!ok_htmltools) {\n",
  "  cat('htmltools is not available in this R session; interactive plots cannot be emitted reliably.')\n",
  "} else {\n",
  "  rds_files <- inv$file[grepl('\\\\.rds$', inv$file)]\n",
  "  if (length(rds_files) == 0) {\n",
  "    cat('No RDS plot objects found (no *.rds files in inventory). If you want interactive plots, save ggplot objects to .rds alongside the PNGs.')\n",
  "  } else {\n",
  "    widgets <- list()\n",
  "    for (f in rds_files) {\n",
  "      obj <- tryCatch(readRDS(f), error = function(e) NULL)\n",
  "      if (is.null(obj)) {\n",
  "        widgets[[length(widgets) + 1]] <- htmltools::tags$p(paste0('Could not read RDS: ', f))\n",
  "      } else {\n",
  "        w <- tryCatch({\n",
  "          if (inherits(obj, 'plotly') || inherits(obj, 'htmlwidget')) {\n",
  "            obj\n",
  "          } else {\n",
  "            plotly::ggplotly(obj)\n",
  "          }\n",
  "        }, error = function(e) NULL)\n",
  "        if (is.null(w)) {\n",
  "          widgets[[length(widgets) + 1]] <- htmltools::tags$p(paste0('plotly rendering failed for: ', f))\n",
  "        } else {\n",
  "          widgets[[length(widgets) + 1]] <- htmltools::tagList(htmltools::tags$h4(basename(f)), w)\n",
  "        }\n",
  "      }\n",
  "    }\n",
  "    htmltools::tagList(widgets)\n",
  "  }\n",
  "}\n",
  "```\n\n",
  
  "# Reproducibility\n\n",
  "```{r}\n",
  "sessionInfo()\n",
  "```\n"
)

writeLines(qmd_text, qmd_path)

# --------------------------- Render report (Quarto-first; MANDATORY) ---------------------------

render_ok <- TRUE
render_msg <- NULL

render_with_quarto <- function(qmd_file, out_dir) {
  qbin <- Sys.which("quarto")
  if (!nzchar(qbin)) return(FALSE)
  
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(dirname(qmd_file))
  
  cmd <- c("render", basename(qmd_file), "--to", "html", "--output-dir", out_dir, "--quiet")
  
  # Capture output and check status; do NOT pretend success.
  res <- tryCatch(
    system2(qbin, args = cmd, stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    stop2("Quarto CLI render failed to start: ", conditionMessage(res))
  }
  
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    stop2("Quarto CLI render failed. Output:\n", paste(res, collapse = "\n"))
  }
  
  TRUE
}

tryCatch({
  ok <- render_with_quarto(qmd_path, output_dir)
  if (!ok) {
    stop2(
      "Quarto rendering not available. Install Quarto (CLI) or the R package 'quarto'.\n",
      "Attempted: quarto::quarto_render() then system 'quarto'."
    )
  }
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# --------------------------- Refresh inventory AFTER rendering (MANDATORY) ---------------------------

inv_post_report <- write_inventory_csv(output_dir, manifest_path_csv)

# Refresh README after render
file_desc2 <- data.frame(
  file = inv_post_report$file,
  description = vapply(inv_post_report$file, describe_file, character(1)),
  stringsAsFactors = FALSE
)
readr::write_csv(file_desc2, readme_path)

# Final inventory refresh
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
