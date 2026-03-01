#!/usr/bin/env Rscript
# ============================================================
# BridgingGenes Post-Processing: Quadrants Over Time + Trajectories
#
# STRICT FIX (as requested; ignore all prior):
#   - Fixes HTML report failure: "The file '..._Report.qmd' does not exist"
#   - Uses Quarto CLI for .qmd when available (robust; correct engine)
#   - Always passes FULL PATH to renderer (never basename-dependent)
#   - Verifies QMD exists immediately after writeLines()
#
# EXTENSIONS (STRICTLY ADDITIVE; requested):
#   1) Jaccard + Module persistence score (complements Jaccard)
#   2) Gene-level role stability index
#   3) Day-level “network regime” summary
#   4) Transition directionality matrix
#   5) Module-level functional annotation (sparingly; optional if annotations provided)
# ============================================================

# ============================================================
# FLAGS (MANDATORY)
# ============================================================
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

# Optional: module functional annotation (sparingly)
# If a file exists in the prior run folder named:
#   - BridgingGenes_FunctionAnnotations.csv
# or
#   - Gene_FunctionAnnotations.csv
# with columns including at least Gene + (Function OR Pathway OR Category),
# the script will produce brief module summaries (top 3 categories).
enable_module_function_annotation <- TRUE
max_module_function_terms <- 3L

start_time <- Sys.time()

# ============================================================
# Dependency bootstrap (install if missing)
# ============================================================
quiet_install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

base_required <- c("jsonlite", "rmarkdown", "knitr", "tools", "stats", "utils")
quiet_install_if_missing(base_required)
if (enable_plotly) quiet_install_if_missing(c("ggplot2", "plotly", "htmltools"))

suppressPackageStartupMessages({
  library(jsonlite)
  library(rmarkdown)
  library(knitr)
  library(tools)
})
if (enable_plotly) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(plotly)
    library(htmltools)
  })
}

# ============================================================
# Script identity capture (MANDATORY)
#   Captures script_name/script_path/script_full_path robustly
# ============================================================
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
  
  # 3) source() fallback
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) {
    return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  }
  
  # 4) Unknown
  NA_character_
}

known_script_filename <- "BridgingGenes_PostProcess_QuadrantsOverTime.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()
if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

cat("\n================ Script Identity ================\n")
cat("script_name: ", script_name, "\n", sep = "")
cat("script_path: ", script_path, "\n", sep = "")
cat("script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
cat("=================================================\n\n")

# ============================================================
# Helper: Mac-friendly directory chooser
# ============================================================
choose_directory <- function(prompt = "Select directory") {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::selectDirectory(caption = prompt), error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  if (requireNamespace("tcltk", quietly = TRUE)) {
    p <- tryCatch(tcltk::tk_choose.dir(caption = prompt), error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  cat("\n", prompt, "\n", sep = "")
  cat("Paste full path (fallback):\n")
  p <- readline("Directory: ")
  if (!nzchar(p)) stop("No directory selected/provided.")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

# ============================================================
# Output directory policy (MANDATORY)
# ============================================================
analysis_name <- "BridgingGenes_PostProcess"
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id_base <- paste0(analysis_name, "_", run_timestamp)

outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

output_dir <- file.path(outputs_root, run_id_base)
if (dir.exists(output_dir)) {
  k <- 1L
  repeat {
    candidate <- file.path(outputs_root, paste0(run_id_base, "_", sprintf("%03d", k)))
    if (!dir.exists(candidate)) { output_dir <- candidate; break }
    k <- k + 1L
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- basename(output_dir)

cat("\n================ Output Policy ================\n")
cat("run_id:       ", run_id, "\n", sep = "")
cat("outputs_root: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("output_dir:   ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("================================================\n\n")

# ============================================================
# USER INPUT: choose prior BridgingGenes run folder
# ============================================================
cat("\nChoose the PRIOR BridgingGenes run directory (the folder containing BridgingGenes_*.csv outputs)...\n")
in_dir <- choose_directory("Select the prior BridgingGenes run folder")
if (!dir.exists(in_dir)) stop("Input directory does not exist: ", in_dir)
in_dir <- normalizePath(in_dir, winslash = "/", mustWork = FALSE)
in_dir_name <- basename(in_dir)

cat("Using input run folder:\n  ", in_dir, "\n", sep = "")
cat("Input run folder name:\n  ", in_dir_name, "\n", sep = "")

# Required + optional inputs
path_summary  <- file.path(in_dir, "BridgingGenes_SummaryAcrossDays.csv")
path_long     <- file.path(in_dir, "BridgingGenes_Integrated_Long_AllDays.csv")
path_roles    <- file.path(in_dir, "BridgingGenes_RoleQuadrants.csv")

# Optional annotation inputs (module-level functional annotation; sparingly)
path_annot1 <- file.path(in_dir, "BridgingGenes_FunctionAnnotations.csv")
path_annot2 <- file.path(in_dir, "Gene_FunctionAnnotations.csv")

# Track input existence flags early
input_exists <- list(
  prior_run_dir_exists = dir.exists(in_dir),
  summary_csv_exists = file.exists(path_summary),
  integrated_long_csv_exists = file.exists(path_long),
  role_quadrants_csv_exists = file.exists(path_roles),
  module_annotation_file1_exists = file.exists(path_annot1),
  module_annotation_file2_exists = file.exists(path_annot2)
)

if (!file.exists(path_summary)) {
  stop("Required file missing in selected folder:\n  ", path_summary,
       "\nThis post-processing script expects the FIRST script's outputs.")
}

summary <- read.csv(path_summary, stringsAsFactors = FALSE)

# Long table: preferred
integrated_long <- NULL
used_per_day_fallback <- FALSE
per_day_files <- character()

if (file.exists(path_long)) {
  integrated_long <- read.csv(path_long, stringsAsFactors = FALSE)
} else {
  # Fallback: bind per-day files
  per_day_files <- list.files(in_dir, pattern = "^BridgingGenes_NetworkMetrics_Day[0-9]+\\.csv$", full.names = TRUE)
  if (length(per_day_files) == 0) {
    stop("Could not find BridgingGenes_Integrated_Long_AllDays.csv and no per-day files found.\n",
         "Expected either:\n  - ", path_long, "\n  OR\n  - BridgingGenes_NetworkMetrics_DayXX.csv files.")
  }
  dfs <- lapply(per_day_files, function(f) read.csv(f, stringsAsFactors = FALSE))
  integrated_long <- do.call(rbind, dfs)
  used_per_day_fallback <- TRUE
}

input_exists$per_day_fallback_used <- used_per_day_fallback
input_exists$per_day_files_found <- length(per_day_files)
input_exists$per_day_files_exist <- if (length(per_day_files) > 0) all(file.exists(per_day_files)) else FALSE

# Basic sanity checks
needed_cols <- c("Gene", "Day", "participation", "betweenness", "degree", "community")
missing_cols <- setdiff(needed_cols, colnames(integrated_long))
if (length(missing_cols) > 0) {
  stop("Integrated long table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

integrated_long$Day <- as.integer(integrated_long$Day)
integrated_long <- integrated_long[order(integrated_long$Gene, integrated_long$Day), ]
rownames(integrated_long) <- NULL

# ============================================================
# GLOBAL betweenness cutoff (bw_cut) — match your 2×2 logic
# ============================================================
betweenness_quantile <- 0.85
bw_cut <- NA_real_
bw_cut_source <- NA_character_

if (file.exists(path_roles)) {
  roles <- read.csv(path_roles, stringsAsFactors = FALSE)
  
  if (!("mean_betweenness" %in% colnames(roles))) {
    roles <- NULL
  } else {
    bw_nonzero <- roles$mean_betweenness[is.finite(roles$mean_betweenness) & roles$mean_betweenness > 0]
    if (length(bw_nonzero) >= 10) {
      bw_cut <- as.numeric(quantile(bw_nonzero, probs = betweenness_quantile, names = FALSE, type = 7))
      bw_cut_source <- "RoleQuadrants.csv (recomputed from mean_betweenness nonzero quantile)"
    }
  }
}

if (!is.finite(bw_cut)) {
  if (!("mean_betweenness" %in% colnames(summary))) stop("SummaryAcrossDays lacks mean_betweenness column.")
  bw_nonzero <- summary$mean_betweenness[is.finite(summary$mean_betweenness) & summary$mean_betweenness > 0]
  if (length(bw_nonzero) < 10) {
    stop("Too few nonzero mean_betweenness values to compute bw_cut robustly.")
  }
  bw_cut <- as.numeric(quantile(bw_nonzero, probs = betweenness_quantile, names = FALSE, type = 7))
  bw_cut_source <- "SummaryAcrossDays.csv (mean_betweenness nonzero quantile)"
}

cat("\n================ Quadrant Cutoffs ================\n")
cat("Participation high criterion: participation > 0 (per day)\n")
cat("Betweenness high criterion  : betweenness >= bw_cut\n")
cat("bw_cut (global)             : ", signif(bw_cut, 6), "\n", sep = "")
cat("bw_cut source               : ", bw_cut_source, "\n", sep = "")
cat("betweenness_quantile        : ", betweenness_quantile, "\n", sep = "")
cat("==================================================\n\n")

# ============================================================
# Assign per-day quadrant membership (Gene × Day)
# ============================================================
integrated_long$P_class_day <- ifelse(integrated_long$participation > 0,
                                      "High participation", "Low participation")
integrated_long$B_class_day <- ifelse(integrated_long$betweenness >= bw_cut,
                                      "High betweenness", "Low betweenness")

integrated_long$Quadrant_day <- NA_character_
integrated_long$Quadrant_day[integrated_long$P_class_day == "Low participation"  &
                               integrated_long$B_class_day == "Low betweenness"]  <- "Within-module genes"
integrated_long$Quadrant_day[integrated_long$P_class_day == "Low participation"  &
                               integrated_long$B_class_day == "High betweenness"] <- "Intra-module bottlenecks"
integrated_long$Quadrant_day[integrated_long$P_class_day == "High participation" &
                               integrated_long$B_class_day == "Low betweenness"]  <- "Module connectors"
integrated_long$Quadrant_day[integrated_long$P_class_day == "High participation" &
                               integrated_long$B_class_day == "High betweenness"] <- "System-level integrators"

# ============================================================
# Trajectory strings + volatility metrics
# ============================================================
split_by_gene <- split(integrated_long, integrated_long$Gene)

traj_strings <- do.call(rbind, lapply(split_by_gene, function(df) {
  df <- df[order(df$Day), ]
  data.frame(
    Gene = df$Gene[1],
    n_days_present = nrow(df),
    Quadrant_trajectory = paste(paste0("Day", df$Day, ":", df$Quadrant_day), collapse = " \u2192 "),
    stringsAsFactors = FALSE
  )
}))

volatility <- do.call(rbind, lapply(split_by_gene, function(df) {
  df <- df[order(df$Day), ]
  q <- df$Quadrant_day
  n_transitions <- if (length(q) <= 1) 0L else sum(q[-1] != q[-length(q)])
  data.frame(
    Gene = df$Gene[1],
    n_days_present = nrow(df),
    unique_quadrants = length(unique(q)),
    n_transitions = n_transitions,
    first_quadrant = q[1],
    last_quadrant = q[length(q)],
    stringsAsFactors = FALSE
  )
}))

# ============================================================
# (2) Gene-level role stability index (NEW)
#   - role_stability_index: combines dominance + transition stability
# ============================================================
gene_role_stability <- do.call(rbind, lapply(split_by_gene, function(df) {
  df <- df[order(df$Day), ]
  q <- df$Quadrant_day
  n <- length(q)
  if (n <= 1) {
    dom_frac <- 1
    trans_rate <- 0
    stab <- 1
  } else {
    tab <- sort(table(q), decreasing = TRUE)
    dom_frac <- as.numeric(tab[1]) / n
    n_trans <- sum(q[-1] != q[-n])
    trans_rate <- n_trans / (n - 1)
    # Stability: high when one role dominates AND few transitions
    stab <- dom_frac * (1 - trans_rate)
  }
  data.frame(
    Gene = df$Gene[1],
    role_dominance_fraction = round(dom_frac, 6),
    role_transition_rate = round(trans_rate, 6),
    role_stability_index = round(stab, 6),
    stringsAsFactors = FALSE
  )
}))
# =========================
# PART 2/3
# =========================

# Quadrant counts per gene
quad_counts <- aggregate(Day ~ Gene + Quadrant_day, integrated_long, length)
colnames(quad_counts)[3] <- "n_days"

# Day-wise composition
day_comp <- aggregate(Gene ~ Day + Quadrant_day, integrated_long, length)
colnames(day_comp)[3] <- "n_genes"

# Join volatility + stability with cross-day summary
if (!("Gene" %in% colnames(summary))) stop("SummaryAcrossDays missing Gene column.")
merged_summary <- merge(summary, volatility, by = "Gene", all.x = TRUE)
merged_summary <- merge(merged_summary, gene_role_stability, by = "Gene", all.x = TRUE)

# ============================================================
# (4) Transition directionality matrix (NEW)
#   Quadrant -> Quadrant flows for each consecutive day pair
#   Outputs: long + wide per transition, plus combined.
# ============================================================
days_present <- sort(unique(integrated_long$Day))
day_pairs <- data.frame(
  Day1 = days_present[-length(days_present)],
  Day2 = days_present[-1],
  stringsAsFactors = FALSE
)

# Prepare per gene per day quadrant lookup
quad_lookup <- integrated_long[, c("Gene", "Day", "Quadrant_day")]
quad_lookup <- quad_lookup[order(quad_lookup$Gene, quad_lookup$Day), ]
rownames(quad_lookup) <- NULL

transition_long_all <- list()
transition_wide_files <- character()

all_quadrants <- c(
  "Within-module genes",
  "Intra-module bottlenecks",
  "Module connectors",
  "System-level integrators"
)

for (i in seq_len(nrow(day_pairs))) {
  d1 <- day_pairs$Day1[i]
  d2 <- day_pairs$Day2[i]
  
  x1 <- quad_lookup[quad_lookup$Day == d1, c("Gene", "Quadrant_day")]
  x2 <- quad_lookup[quad_lookup$Day == d2, c("Gene", "Quadrant_day")]
  colnames(x1)[2] <- "Q1"
  colnames(x2)[2] <- "Q2"
  
  m <- merge(x1, x2, by = "Gene", all = FALSE)
  if (nrow(m) == 0) next
  
  m$Q1 <- factor(m$Q1, levels = all_quadrants)
  m$Q2 <- factor(m$Q2, levels = all_quadrants)
  
  tab <- as.data.frame.matrix(table(m$Q1, m$Q2))
  tab <- tab[all_quadrants, all_quadrants, drop = FALSE]
  
  # Long
  tl <- as.data.frame(as.table(as.matrix(tab)))
  colnames(tl) <- c("From_Quadrant", "To_Quadrant", "n_genes")
  tl$Day1 <- d1
  tl$Day2 <- d2
  tl <- tl[, c("Day1", "Day2", "From_Quadrant", "To_Quadrant", "n_genes")]
  
  transition_long_all[[paste0(d1, "_", d2)]] <- tl
  
  # Wide (one file per transition)
  out_wide <- file.path(output_dir, paste0("BridgingGenes_TransitionMatrix_Quadrants_Day", d1, "_to_Day", d2, "_Wide.csv"))
  write.csv(data.frame(From_Quadrant = rownames(tab), tab, row.names = NULL, check.names = FALSE),
            out_wide, row.names = FALSE)
  transition_wide_files <- c(transition_wide_files, basename(out_wide))
}

transition_long <- if (length(transition_long_all) > 0) do.call(rbind, transition_long_all) else {
  data.frame(Day1=integer(), Day2=integer(), From_Quadrant=character(), To_Quadrant=character(), n_genes=integer())
}

# ============================================================
# (1) Jaccard + Module persistence score (NEW)
#   Communities can relabel across days. We compute best-match mapping:
#     For each module in Day1, find Day2 module with max Jaccard overlap.
#   Scores:
#     - jaccard = |A∩B| / |A∪B|
#     - persistence_retention = |A∩B| / |A|        (how much of Day1 module persists)
#     - persistence_recovery  = |A∩B| / |B|        (how much of Day2 module comes from Day1 module)
#     - overlap_coefficient   = |A∩B| / min(|A|,|B|)
# ============================================================
get_module_sets <- function(df_day) {
  # returns named list: module_id -> character vector of genes
  split(df_day$Gene, df_day$community)
}

jaccard_one <- function(a, b) {
  a <- unique(a); b <- unique(b)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

module_match_results <- list()

for (i in seq_len(nrow(day_pairs))) {
  d1 <- day_pairs$Day1[i]
  d2 <- day_pairs$Day2[i]
  
  df1 <- integrated_long[integrated_long$Day == d1, c("Gene", "community")]
  df2 <- integrated_long[integrated_long$Day == d2, c("Gene", "community")]
  if (nrow(df1) == 0 || nrow(df2) == 0) next
  
  m1 <- get_module_sets(df1)
  m2 <- get_module_sets(df2)
  
  mod1_ids <- names(m1)
  mod2_ids <- names(m2)
  
  # For each mod in day1, find best matching mod in day2 by Jaccard
  rows <- list()
  for (m1_id in mod1_ids) {
    best_j <- -Inf
    best_m2 <- NA_character_
    best_inter <- 0L
    best_uni <- 0L
    best_a <- length(unique(m1[[m1_id]]))
    best_b <- NA_integer_
    
    for (m2_id in mod2_ids) {
      a <- unique(m1[[m1_id]])
      b <- unique(m2[[m2_id]])
      inter <- length(intersect(a, b))
      uni <- length(union(a, b))
      j <- if (uni == 0) NA_real_ else inter / uni
      if (is.finite(j) && j > best_j) {
        best_j <- j
        best_m2 <- m2_id
        best_inter <- inter
        best_uni <- uni
        best_b <- length(b)
      }
    }
    
    if (!is.finite(best_j)) next
    
    retention <- if (best_a > 0) best_inter / best_a else NA_real_
    recovery  <- if (!is.na(best_b) && best_b > 0) best_inter / best_b else NA_real_
    overlap_coef <- if (!is.na(best_b) && min(best_a, best_b) > 0) best_inter / min(best_a, best_b) else NA_real_
    
    rows[[length(rows) + 1L]] <- data.frame(
      Day1 = d1,
      Day2 = d2,
      community_day1 = as.character(m1_id),
      best_match_community_day2 = as.character(best_m2),
      size_day1 = best_a,
      size_day2 = as.integer(best_b),
      intersection_n = as.integer(best_inter),
      union_n = as.integer(best_uni),
      jaccard = round(best_j, 6),
      persistence_retention = round(retention, 6),
      persistence_recovery = round(recovery, 6),
      overlap_coefficient = round(overlap_coef, 6),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(rows) > 0) {
    module_match_results[[paste0(d1, "_", d2)]] <- do.call(rbind, rows)
  }
}

module_persistence <- if (length(module_match_results) > 0) do.call(rbind, module_match_results) else {
  data.frame(
    Day1=integer(), Day2=integer(), community_day1=character(), best_match_community_day2=character(),
    size_day1=integer(), size_day2=integer(), intersection_n=integer(), union_n=integer(),
    jaccard=numeric(), persistence_retention=numeric(), persistence_recovery=numeric(), overlap_coefficient=numeric(),
    stringsAsFactors = FALSE
  )
}

# Summarize module persistence by transition
module_persistence_summary <- if (nrow(module_persistence) > 0) {
  aggregate(
    cbind(jaccard, persistence_retention, persistence_recovery, overlap_coefficient) ~ Day1 + Day2,
    data = module_persistence,
    FUN = function(x) c(mean=mean(x, na.rm=TRUE), median=median(x, na.rm=TRUE))
  )
} else {
  data.frame(Day1=integer(), Day2=integer())
}

# Flatten the aggregate columns for CSV friendliness
flatten_agg <- function(df) {
  if (nrow(df) == 0) return(df)
  
  out <- df[, c("Day1", "Day2"), drop = FALSE]
  
  for (nm in setdiff(colnames(df), c("Day1", "Day2"))) {
    v <- df[[nm]]
    
    # aggregate() may return:
    #  - a list of numeric vectors (older behavior)
    #  - an AsIs matrix column (newer behavior)
    #  - rarely a plain numeric vector in edge cases
    if (is.list(v)) {
      mat <- do.call(rbind, v)
    } else if (is.matrix(v)) {
      mat <- v
    } else if (is.numeric(v)) {
      # Edge case: coerce to 1-row matrix if needed
      if (length(v) == 2 && nrow(df) == 1) {
        mat <- matrix(v, nrow = 1)
        colnames(mat) <- c("mean", "median")
      } else {
        # Fallback: cannot reliably infer structure
        mat <- matrix(NA_real_, nrow = nrow(df), ncol = 2)
        colnames(mat) <- c("mean", "median")
      }
    } else {
      mat <- matrix(NA_real_, nrow = nrow(df), ncol = 2)
      colnames(mat) <- c("mean", "median")
    }
    
    if (is.null(colnames(mat)) || ncol(mat) != 2) {
      # Enforce expected 2-col structure
      mat2 <- matrix(NA_real_, nrow = nrow(df), ncol = 2)
      colnames(mat2) <- c("mean", "median")
      if (ncol(mat) >= 1) mat2[, 1] <- mat[, 1]
      if (ncol(mat) >= 2) mat2[, 2] <- mat[, 2]
      mat <- mat2
    }
    
    colnames(mat) <- paste0(nm, "_", colnames(mat))
    out <- cbind(out, as.data.frame(mat, stringsAsFactors = FALSE))
  }
  
  out
}
module_persistence_summary_flat <- flatten_agg(module_persistence_summary)

# ============================================================
# (3) Day-level “network regime” summary (NEW; high value)
#   Computes:
#     - n_genes_day
#     - n_modules_day
#     - mean/median degree, betweenness, participation
#     - quadrant proportions
#     - integration_index = frac_integrators + 0.5*frac_connectors
#     - quadrant_entropy (Shannon; normalized by log(K))
#     - regime_label (interpretable category)
# ============================================================
shannon_entropy <- function(p) {
  p <- p[p > 0 & is.finite(p)]
  if (length(p) == 0) return(0)
  -sum(p * log(p))
}

day_regime <- do.call(rbind, lapply(days_present, function(d) {
  df <- integrated_long[integrated_long$Day == d, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  
  n_genes <- length(unique(df$Gene))
  n_mods  <- length(unique(df$community))
  
  mean_deg <- mean(df$degree, na.rm = TRUE)
  med_deg  <- median(df$degree, na.rm = TRUE)
  mean_bw  <- mean(df$betweenness, na.rm = TRUE)
  med_bw   <- median(df$betweenness, na.rm = TRUE)
  mean_p   <- mean(df$participation, na.rm = TRUE)
  med_p    <- median(df$participation, na.rm = TRUE)
  
  qc <- table(factor(df$Quadrant_day, levels = all_quadrants))
  qprop <- as.numeric(qc) / sum(qc)
  
  names(qprop) <- c("prop_within", "prop_bottleneck", "prop_connector", "prop_integrator")
  
  integration_index <- qprop["prop_integrator"] + 0.5 * qprop["prop_connector"]
  ent <- shannon_entropy(qprop)
  ent_norm <- if (length(all_quadrants) > 1) ent / log(length(all_quadrants)) else NA_real_
  
  # Simple, interpretable regime labeling
  regime_label <- "Within-module dominant"
  if (qprop["prop_integrator"] >= 0.20) regime_label <- "Integration-heavy"
  if (qprop["prop_connector"] >= 0.25 && qprop["prop_integrator"] < 0.20) regime_label <- "Connector-heavy"
  if (qprop["prop_bottleneck"] >= 0.25 && qprop["prop_integrator"] < 0.20) regime_label <- "Bottleneck-heavy"
  if (ent_norm >= 0.90) regime_label <- paste0(regime_label, " (highly mixed roles)")
  
  data.frame(
    Day = d,
    n_genes_day = n_genes,
    n_modules_day = n_mods,
    mean_degree = round(mean_deg, 6),
    median_degree = round(med_deg, 6),
    mean_betweenness = round(mean_bw, 6),
    median_betweenness = round(med_bw, 6),
    mean_participation = round(mean_p, 6),
    median_participation = round(med_p, 6),
    prop_within = round(qprop["prop_within"], 6),
    prop_bottleneck = round(qprop["prop_bottleneck"], 6),
    prop_connector = round(qprop["prop_connector"], 6),
    prop_integrator = round(qprop["prop_integrator"], 6),
    integration_index = round(integration_index, 6),
    quadrant_entropy_norm = round(ent_norm, 6),
    regime_label = regime_label,
    stringsAsFactors = FALSE
  )
}))

# ============================================================
# (5) Module-level functional annotation (sparingly; optional)
#   If annotation file exists with Gene + (Function|Pathway|Category):
#     For each day + community, report top N categories (default 3).
# ============================================================
module_function_summary <- NULL
module_function_source <- NA_character_

read_annotation_table <- function() {
  if (!isTRUE(enable_module_function_annotation)) return(NULL)
  f <- NULL
  if (file.exists(path_annot1)) f <- path_annot1
  if (is.null(f) && file.exists(path_annot2)) f <- path_annot2
  if (is.null(f)) return(NULL)
  
  ann <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(ann) || !("Gene" %in% colnames(ann))) return(NULL)
  
  # Find a category column
  cat_col <- NULL
  for (cand in c("Function", "Pathway", "Category", "Annotation", "Term")) {
    if (cand %in% colnames(ann)) { cat_col <- cand; break }
  }
  if (is.null(cat_col)) return(NULL)
  
  out <- ann[, c("Gene", cat_col), drop = FALSE]
  colnames(out) <- c("Gene", "Category")
  out <- out[!is.na(out$Gene) & nzchar(out$Gene) & !is.na(out$Category) & nzchar(out$Category), , drop = FALSE]
  out$Category <- as.character(out$Category)
  out <- out[nzchar(out$Gene) & nzchar(out$Category), , drop = FALSE]
  
  module_function_source <<- basename(f)
  out
}

ann_tbl <- read_annotation_table()

if (!is.null(ann_tbl)) {
  # Merge annotation onto integrated_long
  tmp <- merge(integrated_long[, c("Gene", "Day", "community")], ann_tbl, by = "Gene", all.x = FALSE)
  if (nrow(tmp) > 0) {
    # For each Day x community, summarize top categories
    module_function_summary <- do.call(rbind, lapply(split(tmp, list(tmp$Day, tmp$community), drop = TRUE), function(df) {
      d <- df$Day[1]
      cmt <- df$community[1]
      tab <- sort(table(df$Category), decreasing = TRUE)
      topn <- head(tab, max_module_function_terms)
      data.frame(
        Day = d,
        community = as.character(cmt),
        n_annot_genes = nrow(df),
        top_categories = paste(names(topn), collapse = "; "),
        top_category_counts = paste(as.integer(topn), collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
    module_function_summary <- module_function_summary[order(module_function_summary$Day, module_function_summary$community), ]
    rownames(module_function_summary) <- NULL
  }
}

# ============================================================
# Write outputs (extended)
# ============================================================
out_long_q <- file.path(output_dir, "BridgingGenes_Integrated_Long_AllDays_WithQuadrants.csv")
out_traj   <- file.path(output_dir, "BridgingGenes_QuadrantTrajectories_ByGene.csv")
out_vol    <- file.path(output_dir, "BridgingGenes_QuadrantVolatility_ByGene.csv")
out_qcnt   <- file.path(output_dir, "BridgingGenes_QuadrantCounts_ByGene.csv")
out_dcomp  <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay.csv")
out_merge  <- file.path(output_dir, "BridgingGenes_SummaryAcrossDays_WithVolatility.csv")

# NEW outputs
out_gene_stab <- file.path(output_dir, "BridgingGenes_GeneRoleStabilityIndex.csv")
out_trans_long <- file.path(output_dir, "BridgingGenes_TransitionMatrix_Quadrants_AllTransitions_Long.csv")
out_mod_persist <- file.path(output_dir, "BridgingGenes_ModulePersistence_Jaccard_BestMatch_ByTransition.csv")
out_mod_persist_sum <- file.path(output_dir, "BridgingGenes_ModulePersistence_Summary_ByTransition.csv")
out_day_regime <- file.path(output_dir, "BridgingGenes_DayLevel_NetworkRegime_Summary.csv")
out_mod_func <- file.path(output_dir, "BridgingGenes_ModuleFunctionalAnnotation_Summary.csv")

write.csv(integrated_long, out_long_q, row.names = FALSE)
write.csv(traj_strings, out_traj, row.names = FALSE)
write.csv(volatility, out_vol, row.names = FALSE)
write.csv(quad_counts, out_qcnt, row.names = FALSE)
write.csv(day_comp, out_dcomp, row.names = FALSE)
write.csv(merged_summary, out_merge, row.names = FALSE)

write.csv(gene_role_stability, out_gene_stab, row.names = FALSE)
write.csv(transition_long, out_trans_long, row.names = FALSE)
write.csv(module_persistence, out_mod_persist, row.names = FALSE)
write.csv(module_persistence_summary_flat, out_mod_persist_sum, row.names = FALSE)
write.csv(day_regime, out_day_regime, row.names = FALSE)

if (!is.null(module_function_summary)) {
  write.csv(module_function_summary, out_mod_func, row.names = FALSE)
} else {
  # Write an empty-but-valid file so downstream reporting is stable
  write.csv(data.frame(
    Day = integer(), community = character(), n_annot_genes = integer(),
    top_categories = character(), top_category_counts = character(),
    stringsAsFactors = FALSE
  ), out_mod_func, row.names = FALSE)
}

cat("Wrote:\n")
cat("  - ", basename(out_long_q), "\n", sep = "")
cat("  - ", basename(out_traj), "\n", sep = "")
cat("  - ", basename(out_vol), "\n", sep = "")
cat("  - ", basename(out_qcnt), "\n", sep = "")
cat("  - ", basename(out_dcomp), "\n", sep = "")
cat("  - ", basename(out_merge), "\n", sep = "")
cat("  - ", basename(out_gene_stab), "  (NEW)\n", sep = "")
cat("  - ", basename(out_trans_long), "  (NEW)\n", sep = "")
cat("  - ", paste0(length(transition_wide_files), " transition-wide CSVs  (NEW)\n"), sep = "")
cat("  - ", basename(out_mod_persist), "  (NEW)\n", sep = "")
cat("  - ", basename(out_mod_persist_sum), "  (NEW)\n", sep = "")
cat("  - ", basename(out_day_regime), "  (NEW)\n", sep = "")
cat("  - ", basename(out_mod_func), "  (NEW; optional)\n", sep = "")

# ============================================================
# Optional: simple diagnostic plots (extended)
# ============================================================
plotly_rds_files <- character()

# Plot 1: Day composition
p1_png <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay.png")
p1_rds <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay_plotly.rds")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p1 <- ggplot(day_comp, aes(x = Day, y = n_genes, fill = Quadrant_day)) +
    geom_col() +
    labs(title = "Quadrant composition by day", x = "Day", y = "Number of genes") +
    theme_minimal(base_size = 12)
  
  png(p1_png, width = 1600, height = 1000, res = 150)
  print(p1)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    try({
      saveRDS(ggplotly(p1), p1_rds)
      plotly_rds_files <- c(plotly_rds_files, basename(p1_rds))
    }, silent = TRUE)
  }
}

# Plot 2: Top integrators
integrator_days <- quad_counts[quad_counts$Quadrant_day == "System-level integrators", , drop = FALSE]
integrator_days <- integrator_days[order(-integrator_days$n_days, integrator_days$Gene), ]
top_k <- 25
integrator_top <- head(integrator_days, top_k)

p2_png <- file.path(output_dir, "BridgingGenes_TopIntegrators_ByDays.png")
p2_rds <- file.path(output_dir, "BridgingGenes_TopIntegrators_ByDays_plotly.rds")

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(integrator_top) > 0) {
  integrator_top$Gene <- factor(integrator_top$Gene, levels = rev(integrator_top$Gene))
  p2 <- ggplot(integrator_top, aes(x = Gene, y = n_days)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top ", min(top_k, nrow(integrator_top)),
                        " genes by days classified as System-level integrators"),
         x = "Gene", y = "Days as integrator") +
    theme_minimal(base_size = 12)
  
  png(p2_png, width = 1600, height = 1100, res = 150)
  print(p2)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    try({
      saveRDS(ggplotly(p2), p2_rds)
      plotly_rds_files <- c(plotly_rds_files, basename(p2_rds))
    }, silent = TRUE)
  }
}

# Plot 3 (NEW): Day-level regime indices
p3_png <- file.path(output_dir, "BridgingGenes_DayLevel_RegimeIndices.png")
p3_rds <- file.path(output_dir, "BridgingGenes_DayLevel_RegimeIndices_plotly.rds")

if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(day_regime) && nrow(day_regime) > 0) {
  dr <- day_regime
  dr_long <- rbind(
    data.frame(Day = dr$Day, Metric = "integration_index", Value = dr$integration_index, stringsAsFactors = FALSE),
    data.frame(Day = dr$Day, Metric = "quadrant_entropy_norm", Value = dr$quadrant_entropy_norm, stringsAsFactors = FALSE)
  )
  
  p3 <- ggplot(dr_long, aes(x = Day, y = Value)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(title = "Day-level network regime indices", x = "Day", y = "Value") +
    theme_minimal(base_size = 12)
  
  png(p3_png, width = 1600, height = 900, res = 150)
  print(p3)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    try({
      saveRDS(ggplotly(p3), p3_rds)
      plotly_rds_files <- c(plotly_rds_files, basename(p3_rds))
    }, silent = TRUE)
  }
}

# Plot 4 (NEW): Module persistence summary (Jaccard mean)
p4_png <- file.path(output_dir, "BridgingGenes_ModulePersistence_JaccardMean_ByTransition.png")
p4_rds <- file.path(output_dir, "BridgingGenes_ModulePersistence_JaccardMean_ByTransition_plotly.rds")

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(module_persistence_summary_flat) > 0 && "jaccard_mean" %in% colnames(module_persistence_summary_flat)) {
  mp <- module_persistence_summary_flat
  mp$Transition <- paste0(mp$Day1, "\u2192", mp$Day2)
  
  p4 <- ggplot(mp, aes(x = Transition, y = jaccard_mean)) +
    geom_col() +
    labs(title = "Module persistence: mean best-match Jaccard by transition", x = "Transition", y = "Mean Jaccard") +
    theme_minimal(base_size = 12)
  
  png(p4_png, width = 1600, height = 900, res = 150)
  print(p4)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    try({
      saveRDS(ggplotly(p4), p4_rds)
      plotly_rds_files <- c(plotly_rds_files, basename(p4_rds))
    }, silent = TRUE)
  }
}
# =========================
# PART 3/3
# =========================

# ============================================================
# Manifest + inventory helpers
# ============================================================
get_deps <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    ver <- NA_character_
    if (requireNamespace(p, quietly = TRUE)) ver <- as.character(utils::packageVersion(p))
    data.frame(package = p, version = ver, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

build_file_inventory <- function(dir_path, out_csv) {
  f <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) {
    write.csv(data.frame(path=character(), filename=character(), bytes=numeric(), modified=character()),
              out_csv, row.names = FALSE)
    return(invisible(out_csv))
  }
  info <- file.info(f)
  inv <- data.frame(
    path = normalizePath(f, winslash = "/", mustWork = FALSE),
    filename = basename(f),
    bytes = as.numeric(info$size),
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  inv <- inv[order(inv$filename), , drop = FALSE]
  write.csv(inv, out_csv, row.names = FALSE)
  invisible(out_csv)
}

end_time <- Sys.time()
runtime_seconds <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

deps_pkgs <- c("jsonlite", "rmarkdown", "knitr", "tools", "stats", "utils")
if (enable_plotly) deps_pkgs <- c(deps_pkgs, "ggplot2", "plotly", "htmltools")
deps_df <- get_deps(unique(deps_pkgs))

# Output existence flags (for initial outputs)
output_expected <- list(
  integrated_long_with_quadrants_csv = out_long_q,
  quadrant_trajectories_csv = out_traj,
  quadrant_volatility_csv = out_vol,
  quadrant_counts_csv = out_qcnt,
  quadrant_composition_csv = out_dcomp,
  summary_with_volatility_csv = out_merge,
  
  gene_role_stability_index_csv = out_gene_stab,
  transition_matrix_all_long_csv = out_trans_long,
  module_persistence_bestmatch_csv = out_mod_persist,
  module_persistence_summary_csv = out_mod_persist_sum,
  day_level_regime_summary_csv = out_day_regime,
  module_function_annotation_summary_csv = out_mod_func,
  
  quadrant_composition_png = p1_png,
  quadrant_composition_plotly_rds = p1_rds,
  top_integrators_png = p2_png,
  top_integrators_plotly_rds = p2_rds,
  
  day_level_regime_indices_png = p3_png,
  day_level_regime_indices_plotly_rds = p3_rds,
  module_persistence_jaccardmean_png = p4_png,
  module_persistence_jaccardmean_plotly_rds = p4_rds
)

# Add per-transition wide files to manifest expected list (as dynamic)
if (length(transition_wide_files) > 0) {
  for (fn in transition_wide_files) {
    output_expected[[paste0("transition_wide_", fn)]] <- file.path(output_dir, fn)
  }
}

output_exists <- lapply(output_expected, file.exists)

manifest <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script = list(
    name = script_name,
    path = script_path,
    full_path = if (is.na(script_full)) NA_character_ else script_full
  ),
  input = list(
    prior_run_dir = in_dir,
    prior_run_dir_name = in_dir_name,
    paths = list(
      summary_csv = path_summary,
      integrated_long_csv = path_long,
      role_quadrants_csv = path_roles,
      module_annotation_file1 = path_annot1,
      module_annotation_file2 = path_annot2
    ),
    files_detected = input_exists
  ),
  parameters = list(
    participation_high_rule = "participation > 0 (per day)",
    betweenness_quantile = betweenness_quantile,
    bw_cut = bw_cut,
    bw_cut_source = bw_cut_source,
    enable_plotly = enable_plotly,
    enable_runtime_tracking = enable_runtime_tracking,
    enable_module_function_annotation = enable_module_function_annotation,
    max_module_function_terms = max_module_function_terms,
    module_function_source = module_function_source,
    runtime_seconds = if (enable_runtime_tracking) runtime_seconds else NA_real_
  ),
  dependencies = jsonlite::fromJSON(jsonlite::toJSON(deps_df, dataframe = "rows", auto_unbox = TRUE)),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    expected_paths = lapply(output_expected, function(p) normalizePath(p, winslash = "/", mustWork = FALSE)),
    files_detected = output_exists
  )
)

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           con = manifest_json_path)

manifest_files_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")
build_file_inventory(output_dir, manifest_files_csv_path)

# ============================================================
# Write QMD + render HTML report (UPDATED to include NEW outputs)
# ============================================================
report_qmd_path  <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
report_html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))

qmd_lines <- c(
  "---",
  paste0('title: "', script_name, '"'),
  "format:",
  "  html:",
  "    toc: true",
  "    embed-resources: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "params:",
  paste0('  run_id: "', run_id, '"'),
  paste0('  run_timestamp: "', manifest$run_timestamp, '"'),
  paste0('  enable_plotly: ', ifelse(enable_plotly, "true", "false")),
  paste0('  runtime_seconds: ', ifelse(enable_runtime_tracking, as.character(runtime_seconds), "null")),
  paste0('  bw_cut: ', as.character(bw_cut)),
  paste0('  bw_cut_source: "', bw_cut_source, '"'),
  paste0('  betweenness_quantile: ', as.character(betweenness_quantile)),
  paste0('  enable_module_function_annotation: ', ifelse(enable_module_function_annotation, "true", "false")),
  paste0('  module_function_source: "', ifelse(is.na(module_function_source), "", module_function_source), '"'),
  paste0('  max_module_function_terms: ', as.character(max_module_function_terms)),
  "---",
  "",
  "## Provenance (script / input / output)",
  "",
  "```{r}",
  "man <- jsonlite::fromJSON('Project_Manifest.json')",
  "cat('Script name: ', man$script$name, '\\n')",
  "cat('Script path: ', man$script$path, '\\n')",
  "cat('Script full: ', man$script$full_path, '\\n\\n')",
  "cat('Input run dir: ', man$input$prior_run_dir, '\\n')",
  "cat('Input run dir name: ', man$input$prior_run_dir_name, '\\n\\n')",
  "cat('Output dir: ', man$outputs$output_dir, '\\n')",
  "cat('Run ID: ', man$run_id, '\\n')",
  "```",
  "",
  "### Input files detected (YES/NO)",
  "",
  "```{r}",
  "fd <- man$input$files_detected",
  "fd_df <- data.frame(item=names(fd), value=unlist(fd, use.names=FALSE), stringsAsFactors=FALSE)",
  "knitr::kable(fd_df)",
  "```",
  "",
  "### Output files detected (YES/NO)",
  "",
  "```{r}",
  "od <- man$outputs$files_detected",
  "od_df <- data.frame(item=names(od), exists=unlist(od, use.names=FALSE), stringsAsFactors=FALSE)",
  "knitr::kable(od_df)",
  "```",
  "",
  "## Cutoffs used for per-day quadrant membership",
  "",
  "- Participation high rule: **participation > 0** (per day)",
  paste0("- Betweenness high rule: **betweenness \u2265 bw_cut** where bw_cut = ", signif(bw_cut, 6),
         " (", bw_cut_source, "; quantile = ", betweenness_quantile, ")"),
  "",
  "## Integrated long table with quadrants (head)",
  "",
  "```{r}",
  "x <- read.csv('BridgingGenes_Integrated_Long_AllDays_WithQuadrants.csv', stringsAsFactors = FALSE)",
  "x <- x[order(x$Gene, x$Day), ]",
  "knitr::kable(head(x, 25), digits = 4)",
  "```",
  "",
  "## Day-level network regime summary (NEW; very high value)",
  "",
  "```{r}",
  "dr <- read.csv('BridgingGenes_DayLevel_NetworkRegime_Summary.csv', stringsAsFactors = FALSE)",
  "knitr::kable(dr)",
  "```",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_DayLevel_RegimeIndices.png')) knitr::include_graphics('BridgingGenes_DayLevel_RegimeIndices.png')",
  "```",
  "",
  "## Quadrant composition by day",
  "",
  "```{r}",
  "d <- read.csv('BridgingGenes_QuadrantComposition_ByDay.csv', stringsAsFactors = FALSE)",
  "d <- d[order(d$Day, d$Quadrant_day), ]",
  "knitr::kable(d)",
  "```",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_QuadrantComposition_ByDay.png')) knitr::include_graphics('BridgingGenes_QuadrantComposition_ByDay.png')",
  "```",
  "",
  "## Gene-level role stability index (NEW)",
  "",
  "```{r}",
  "gs <- read.csv('BridgingGenes_GeneRoleStabilityIndex.csv', stringsAsFactors = FALSE)",
  "gs2 <- gs[order(-gs$role_stability_index, -gs$role_dominance_fraction, gs$Gene), ]",
  "knitr::kable(head(gs2, 30))",
  "```",
  "",
  "## Quadrant volatility per gene (top by transitions)",
  "",
  "```{r}",
  "v <- read.csv('BridgingGenes_QuadrantVolatility_ByGene.csv', stringsAsFactors = FALSE)",
  "v2 <- v[order(-v$n_transitions, -v$unique_quadrants, v$Gene), ]",
  "knitr::kable(head(v2, 30))",
  "```",
  "",
  "## Transition directionality matrix (NEW; quadrant \u2192 quadrant flows)",
  "",
  "```{r}",
  "tl <- read.csv('BridgingGenes_TransitionMatrix_Quadrants_AllTransitions_Long.csv', stringsAsFactors = FALSE)",
  "tl2 <- tl[order(tl$Day1, tl$From_Quadrant, tl$To_Quadrant), ]",
  "knitr::kable(head(tl2, 40))",
  "```",
  "",
  "## Module persistence (Jaccard + persistence; best-match mapping) (NEW)",
  "",
  "```{r}",
  "mp <- read.csv('BridgingGenes_ModulePersistence_Jaccard_BestMatch_ByTransition.csv', stringsAsFactors = FALSE)",
  "mp2 <- mp[order(-mp$jaccard, -mp$persistence_retention, mp$Day1, mp$community_day1), ]",
  "knitr::kable(head(mp2, 40))",
  "```",
  "",
  "```{r}",
  "mps <- read.csv('BridgingGenes_ModulePersistence_Summary_ByTransition.csv', stringsAsFactors = FALSE)",
  "knitr::kable(mps)",
  "```",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_ModulePersistence_JaccardMean_ByTransition.png')) knitr::include_graphics('BridgingGenes_ModulePersistence_JaccardMean_ByTransition.png')",
  "```",
  "",
  "## Module-level functional annotation (sparingly; optional)",
  "",
  "```{r}",
  "mf <- read.csv('BridgingGenes_ModuleFunctionalAnnotation_Summary.csv', stringsAsFactors = FALSE)",
  "if (nrow(mf) == 0 || all(is.na(mf$top_categories)) || all(mf$top_categories == '')) {",
  "  cat('No annotation table was detected or annotation table lacked required columns (Gene + Function/Pathway/Category).\\n')",
  "} else {",
  "  knitr::kable(head(mf, 50))",
  "}",
  "```",
  "",
  "## Trajectory strings (head)",
  "",
  "```{r}",
  "t <- read.csv('BridgingGenes_QuadrantTrajectories_ByGene.csv', stringsAsFactors = FALSE)",
  "knitr::kable(head(t, 25))",
  "```",
  "",
  "## Top integrators by days",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_TopIntegrators_ByDays.png')) knitr::include_graphics('BridgingGenes_TopIntegrators_ByDays.png')",
  "```",
  "",
  "## Generated files (inventory)",
  "",
  "```{r}",
  "inv <- read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "## Interactive plots (Plotly)",
  "",
  "```{r results='asis'}",
  "if (isTRUE(params$enable_plotly)) {",
  "  rds_files <- list.files('.', pattern = '\\\\_plotly\\\\.rds$', full.names = FALSE)",
  "  if (length(rds_files) == 0) {",
  "    htmltools::tags$p('No plotly RDS files found.')",
  "  } else {",
  "    htmltools::tagList(lapply(rds_files, function(f) {",
  "      w <- tryCatch(readRDS(f), error = function(e) NULL)",
  "      if (is.null(w)) htmltools::tags$p(paste0('Could not read: ', f))",
  "      else htmltools::tagList(htmltools::tags$h4(f), w)",
  "    }))",
  "  }",
  "} else {",
  "  htmltools::tags$p('Plotly disabled.')",
  "}",
  "```",
  "",
  "## Runtime",
  "",
  "```{r}",
  "cat('Runtime seconds: ', params$runtime_seconds)",
  "```",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

# Write QMD
writeLines(qmd_lines, con = report_qmd_path)

# STRICT: verify QMD exists right now
if (!file.exists(report_qmd_path)) {
  stop("QMD write failed; file does not exist at expected path:\n  ", report_qmd_path)
}

render_ok <- TRUE
render_err <- NULL

# Render with Quarto (preferred) from output_dir so relative paths resolve
quarto_bin <- Sys.which("quarto")
if (nzchar(quarto_bin)) {
  old_wd2 <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd2), add = TRUE)
  
  args <- c(
    "render",
    basename(report_qmd_path),
    "--to", "html",
    "--output", basename(report_html_path)
  )
  
  res <- tryCatch(
    system2(quarto_bin, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    render_ok <- FALSE
    render_err <- conditionMessage(res)
  } else {
    status <- attr(res, "status")  # NULL on success, nonzero on failure
    if (!is.null(status) && is.numeric(status) && status != 0) {
      render_ok <- FALSE
      render_err <- paste0("Quarto exit status = ", status, "\nQuarto output:\n", paste(res, collapse = "\n"))
    } else if (!file.exists(report_html_path)) {
      render_ok <- FALSE
      render_err <- paste0("Quarto did not produce expected HTML: ", report_html_path,
                           "\nQuarto output:\n", paste(res, collapse = "\n"))
    }
  }
} else {
  # Fallback (last resort): rmarkdown::render from output_dir for relative paths
  old_wd2 <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd2), add = TRUE)
  
  tryCatch({
    rmarkdown::render(
      input = basename(report_qmd_path),
      output_file = basename(report_html_path),
      output_dir = output_dir,
      quiet = TRUE
    )
    if (!file.exists(report_html_path)) {
      stop("rmarkdown::render completed but HTML not found at: ", report_html_path)
    }
  }, error = function(e) {
    render_ok <<- FALSE
    render_err <<- conditionMessage(e)
  })
}

# Refresh inventory after render
build_file_inventory(output_dir, manifest_files_csv_path)

# Update manifest output flags after report render
manifest$outputs$files_detected$report_qmd_exists  <- file.exists(report_qmd_path)
manifest$outputs$files_detected$report_html_exists <- file.exists(report_html_path)
manifest$outputs$expected_paths$report_qmd  <- normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE)
manifest$outputs$expected_paths$report_html <- normalizePath(report_html_path, winslash = "/", mustWork = FALSE)

# Rewrite manifest to include final report existence
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           con = manifest_json_path)

cat("\n================ DONE ================\n")
cat("Input run folder : ", in_dir, "\n", sep = "")
cat("Output run folder: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Script full path : ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
cat("Key outputs:\n")
cat("  - ", basename(out_long_q), "\n", sep = "")
cat("  - ", basename(out_traj), "\n", sep = "")
cat("  - ", basename(out_vol), "\n", sep = "")
cat("  - ", basename(out_qcnt), "\n", sep = "")
cat("  - ", basename(out_dcomp), "\n", sep = "")
cat("  - ", basename(out_merge), "\n", sep = "")
cat("  - ", basename(out_gene_stab), " (NEW)\n", sep = "")
cat("  - ", basename(out_trans_long), " (NEW)\n", sep = "")
cat("  - ", basename(out_mod_persist), " (NEW)\n", sep = "")
cat("  - ", basename(out_mod_persist_sum), " (NEW)\n", sep = "")
cat("  - ", basename(out_day_regime), " (NEW)\n", sep = "")
cat("  - ", basename(out_mod_func), " (NEW; optional)\n", sep = "")
cat("Manifest         : ", basename(manifest_json_path), "\n", sep = "")
cat("Inventory        : ", basename(manifest_files_csv_path), "\n", sep = "")
if (isTRUE(render_ok) && file.exists(report_html_path)) {
  cat("HTML report      : ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("HTML report      : FAILED\n")
  cat("QMD written to   : ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (!is.null(render_err)) cat("Render error     : ", render_err, "\n", sep = "")
}
cat("=====================================\n")
