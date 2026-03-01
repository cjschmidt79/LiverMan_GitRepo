#TripletPermutationsVer13
#!/usr/bin/env Rscript
# ========================================================================
# 🔬 Triplet Interaction Pipeline (Interactive / RStudio)
#    - Genes in rows, samples in columns
#    - A ~ log(B/C) * group  with PERMUTATION P-VALUE and bootstrapping
#    - Every gene tested as A, B, and C (3 orientations / triplet)
# ==============================================================================
# 🧬 Triplet Regulatory Interaction Analysis (v5.6 - Q-value REMOVED, Syntax FIXED)
# ==============================================================================
#
# FUNCTION:
# This script performs a global scan for **3-way gene regulatory interactions**
# using the 'rheostat' model: log(A) ~ log(B/C) * Group. It tests every
# combination of three genes (A, B, C) in three possible orientations.
#
# The goal is to identify triplets where the relationship (slope) between
# log(A) and log(B/C) significantly changes between two defined sample groups.
#
# Key Steps:
# 1. Filtering genes by minimum expression and variance.
# 2. Linear Model fitting, testing the A ~ log(B/C) * Group interaction term.
# 3. Permutation testing for a non-parametric p-value estimate.
# 4. **REMOVED:** False Discovery Rate (FDR) correction (q-value).
# 5. Bootstrapping on the top N significant triplets (filtered by perm_p_interaction).
#
# ------------------------------------------------------------------------------
# 💾 INPUT FILES:
# ------------------------------------------------------------------------------
# - **User-Selected CSV:** A matrix of raw expression counts.
#    - Rows must be **Genes/Labels**.
#    - Columns must be **Samples**.
#    - Column headers must contain a **prefix** (e.g., '1', '2', '4', '8') that
#      allows for assignment into two user-defined groups (e.g., 'early', 'late').
#
# ------------------------------------------------------------------------------
# 📊 OUTPUT FILES (Saved to a time-stamped directory):
# ------------------------------------------------------------------------------
# 1. `all_triplets.csv`
# 2. `significant_triplets.csv`
# 3. `bootstrap_summary_top.csv`
# 4. `final_triplet_plots/*.png`
# 5. `perm_pval_distribution.png`
# 6. `bootstrap_log.csv` (if bootstraps run)
# 7. `target_gene_filter_diagnostics.csv` (if target file used)
# 8. `README.qmd`
# 9. `OUTPUT_INVENTORY.txt`
#
# ------------------------------------------------------------------------------
# NOTE:
# - qvalue/fdrtool are NOT used in this Q-value-free version.
# - No dplyr/tidyverse is used. If dplyr is loaded, the script stops.
# ========================================================================

# ----------------------------- NO dplyr policy -----------------------------
if ("package:dplyr" %in% search()) {
  stop("This script must not be run with dplyr loaded. Please restart R and run again without loading dplyr/tidyverse.",
       call. = FALSE)
}

# ----------------------------- Minimal utilities (biologist-friendly) -----------------------------
has_rstudioapi <- requireNamespace("rstudioapi", quietly = TRUE)

is_interactive_rstudio <- function() {
  interactive() && has_rstudioapi && rstudioapi::isAvailable()
}

pick_file <- function(prompt = "Select a file") {
  cat("\n", prompt, "\n", sep = "")
  if (is_interactive_rstudio()) {
    p <- rstudioapi::selectFile(caption = prompt)
    if (is.null(p) || !nzchar(p)) stop("No file selected.", call. = FALSE)
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  } else {
    p <- file.choose()
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }
}

safe_slug <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- substr(x, 1, 80)
  if (!nzchar(x)) x <- "run"
  x
}

get_repoish_root <- function() {
  out <- tryCatch(system("git rev-parse --show-toplevel", intern = TRUE),
                  error = function(e) character(0))
  if (length(out) == 0 || !nzchar(out[1])) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  normalizePath(out[1], winslash = "/", mustWork = TRUE)
}

capture_script_identity <- function(verbose = TRUE) {
  script_full <- NA_character_
  
  if (is_interactive_rstudio()) {
    ctx <- rstudioapi::getActiveDocumentContext()
    if (!is.null(ctx$path) && nzchar(ctx$path) && file.exists(ctx$path)) {
      script_full <- normalizePath(ctx$path, winslash = "/", mustWork = TRUE)
    }
  }
  
  if (is.na(script_full)) {
    ca <- commandArgs(trailingOnly = FALSE)
    hit <- grep("--file=", ca, value = TRUE)
    if (length(hit) == 1) {
      cand <- sub("--file=", "", hit, fixed = TRUE)
      if (file.exists(cand)) script_full <- normalizePath(cand, winslash = "/", mustWork = TRUE)
    }
  }
  
  if (is.na(script_full)) script_full <- "UNKNOWN_SCRIPT_PATH"
  
  script_path <- if (script_full == "UNKNOWN_SCRIPT_PATH") "UNKNOWN_DIR" else dirname(script_full)
  script_name <- if (script_full == "UNKNOWN_SCRIPT_PATH") "UNKNOWN_SCRIPT" else basename(script_full)
  
  if (verbose) {
    cat("\n[script_name] ", script_name, "\n", sep = "")
    cat("[script_path] ", script_path, "\n", sep = "")
    cat("[script_full] ", script_full, "\n", sep = "")
  }
  
  list(script_name = script_name, script_path = script_path, script_full = script_full)
}

git_safe <- function(cmd) {
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  if (length(out) == 0 || !nzchar(out[1])) return(NA_character_)
  out[1]
}

git_meta <- list(
  repo_root  = git_safe("git rev-parse --show-toplevel"),
  git_branch = git_safe("git rev-parse --abbrev-ref HEAD"),
  git_commit = git_safe("git rev-parse HEAD"),
  git_remote = git_safe("git config --get remote.origin.url")
)

write_output_inventory <- function(dir_path) {
  inv <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  writeLines(inv, con = file.path(dir_path, "OUTPUT_INVENTORY.txt"))
  invisible(TRUE)
}

# Capture provenance early
id <- capture_script_identity(verbose = TRUE)
script_name <- id$script_name
script_path <- id$script_path
script_full <- id$script_full

cat("\n[git repo_root]  ", git_meta$repo_root,  "\n", sep = "")
cat("[git branch]    ", git_meta$git_branch, "\n", sep = "")
cat("[git commit]    ", git_meta$git_commit, "\n", sep = "")
cat("[git remote]    ", git_meta$git_remote, "\n", sep = "")

# ----------------------------- Packages -----------------------------
need <- c("ggplot2", "gridExtra", "furrr", "future", "progressr")
for (p in need) {
  if (!suppressWarnings(require(p, character.only = TRUE, quietly = TRUE))) {
    install.packages(p, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
}
progressr::handlers(global = TRUE)
op_contr <- options(contrasts = c("contr.treatment","contr.poly"))
on.exit(options(op_contr), add = TRUE)

# ----------------------------- I/O setup (UPDATED: outputs/<run_name>_<timestamp>/ under repo root) -----------------------------
cat("\n💾 Select input CSV (rows=genes, cols=samples):\n")
file_path <- pick_file("Select input CSV (rows=genes, cols=samples)")

input_paths <- normalizePath(c(file_path), winslash = "/", mustWork = TRUE)

input_name <- tools::file_path_sans_ext(basename(file_path))
timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")

root_dir <- if (!is.na(git_meta$repo_root)) normalizePath(git_meta$repo_root, winslash = "/", mustWork = TRUE) else get_repoish_root()
outputs_base <- file.path(root_dir, "outputs")
if (!dir.exists(outputs_base)) dir.create(outputs_base, recursive = TRUE, showWarnings = FALSE)

cat("\nEnter run name (e.g., TripletInteractions_D4vsD20)\n")
run_name <- safe_slug(readline("run_name: "))
run_folder <- paste0(run_name, "_", timestamp)

output_dir <- file.path(outputs_base, run_folder)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir   <- file.path(output_dir, "final_triplet_plots")
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)

cat("✅ Outputs will go to:\n  ", normalizePath(output_dir, winslash="/", mustWork=FALSE), "\n", sep = "")

# ----------------------------- Load data -----------------------------
cat("\n🧬 Loading matrix...\n")
data_raw <- read.csv(file_path, row.names = 1, check.names = FALSE)
if (nrow(data_raw) < 3 || ncol(data_raw) < 3) stop("Matrix too small.", call. = FALSE)

sample_names <- colnames(data_raw)

# ----------------------------- Grouping -----------------------------
cat("\n🔠 Define group NUMBER ranges.\n",
    "  Example: early = 1-7, late = 8-15\n", sep="")

extract_number <- function(x) {
  m <- regexpr("^\\d+", x)
  if (m == -1) return(NA_integer_)
  as.integer(regmatches(x, m))
}

sample_numbers <- sapply(sample_names, extract_number)

if (any(is.na(sample_numbers))) {
  warning("Some sample names do not begin with a number and were set to NA.")
}

cat("Sample numbers extracted:\n")
print(sample_numbers)

early_range <- readline("Enter EARLY range (e.g., 1-7): ")
late_range  <- readline("Enter LATE  range (e.g., 8-15): ")

parse_range <- function(txt) {
  parts <- as.integer(strsplit(txt, "-")[[1]])
  if (length(parts) != 2 || any(is.na(parts))) stop("Invalid range: ", txt, call. = FALSE)
  seq(parts[1], parts[2])
}

early_nums <- parse_range(early_range)
late_nums  <- parse_range(late_range)

group_assignments <- ifelse(
  sample_numbers %in% early_nums, "early",
  ifelse(sample_numbers %in% late_nums, "late", NA)
)

valid <- !is.na(group_assignments)
if (!any(valid)) stop("❌ No columns matched the provided ranges.", call. = FALSE)

data_raw     <- data_raw[, valid, drop = FALSE]
sample_names <- colnames(data_raw)
lvl_order    <- unique(group_assignments[valid])
group_factor <- factor(group_assignments[valid], levels = lvl_order)

cat("\n✅ Grouping:\n")
print(table(group_factor))

# ----------------------------- Parameters (interactive) -----------------------------
read_num <- function(prompt, default = NA_real_) {
  x <- readline(prompt)
  if (!nzchar(x) && !is.na(default)) return(default)
  as.numeric(x)
}

min_expr <- read_num("Minimum expression (default 10): ", 10)
if (!is.finite(min_expr)) min_expr <- 10

pseudocount_delta <- read_num("Pseudocount delta for logs (default 1e-4): ", 1e-4)
if (!is.finite(pseudocount_delta) || pseudocount_delta <= 0) pseudocount_delta <- 1e-4

num_permutations <- read_num("Number of permutations (default 1000): ", 1000)
if (!is.finite(num_permutations) || num_permutations < 1) num_permutations <- 1000

num_bootstrap <- read_num("Number of bootstraps per triplet (default 1000): ", 1000)
if (!is.finite(num_bootstrap) || num_bootstrap < 100) num_bootstrap <- 1000

p_val_thresh <- read_num("Permutation P-value threshold (default 0.05): ", 0.05)
if (!is.finite(p_val_thresh) || p_val_thresh <= 0 || p_val_thresh >= 1) p_val_thresh <- 0.05

min_r2_in <- readline("Optional minimum R^2 filter (Enter to skip; e.g., 0.30): ")
min_r2    <- suppressWarnings(as.numeric(min_r2_in))
use_r2    <- is.finite(min_r2) && min_r2 > 0 && min_r2 < 1

# ----------------------------- Variance inspection & cutoff -----------------------------
cat("\n📊 Inspecting gene variance distribution...\n")
var_vec <- apply(data_raw, 1, var)

suppressWarnings({
  hist(log10(var_vec + 1), breaks = 100,
       main = "Log10 Gene Variance Distribution",
       xlab  = "log10(variance)")
})
suggested_cut <- as.numeric(quantile(var_vec, 0.90, na.rm = TRUE))
abline(v = log10(suggested_cut), lty = 2, lwd = 2)
cat("Suggested 90th percentile variance cutoff:", format(suggested_cut, scientific=TRUE), "\n")

xwin <- readline("Limit x-axis? Enter two numbers like '6 9' for 1e6–1e9, or Enter to skip: ")
if (nzchar(xwin)) {
  parts <- suppressWarnings(as.numeric(strsplit(xwin, "\\s+")[[1]]))
  if (length(parts) == 2 && all(is.finite(parts))) {
    try({
      if (.Platform$OS.type == "windows") {
        windows()
      } else {
        if (capabilities("aqua")) quartz() else x11()
      }
    }, silent = TRUE)
    suppressWarnings({
      hist(log10(var_vec + 1), breaks = 100,
           main = "Log10 Gene Variance Distribution (windowed)",
           xlab = "log10(variance)", xlim = range(parts))
    })
    abline(v = log10(suggested_cut), lty = 2, lwd = 2)
  }
}

cat("\nChoose variance cutoff method:\n",
    "  1 = Absolute (e.g., 1e6)\n",
    "  2 = Percentile (e.g., 0.90 for 90th)\n", sep="")
mth <- readline("Method [1/2, default 2]: ")
if (mth != "1") mth <- "2"

if (mth == "1") {
  abs_txt <- readline("Enter absolute variance cutoff (e.g., 1e6) [default 1e6]: ")
  if (!nzchar(abs_txt)) abs_txt <- "1e6"
  min_var_cutoff <- suppressWarnings(as.numeric(abs_txt))
  if (!is.finite(min_var_cutoff)) min_var_cutoff <- 1e6
  cutoff_desc <- paste0("absolute >= ", format(min_var_cutoff, scientific = TRUE))
} else {
  pct_txt <- readline("Enter percentile in (0,1), e.g. 0.90 [default 0.90]: ")
  if (!nzchar(pct_txt)) pct_txt <- "0.90"
  pct <- suppressWarnings(as.numeric(pct_txt))
  if (!is.finite(pct) || pct <= 0 || pct >= 1) pct <- 0.90
  min_var_cutoff <- as.numeric(quantile(var_vec, probs = pct, na.rm = TRUE))
  cutoff_desc <- paste0("percentile ", pct, " => >= ", format(min_var_cutoff, scientific = TRUE))
}
cat("\n✅ Variance cutoff chosen: ", cutoff_desc, "\n", sep="")

# ----------------------------- Filtering -----------------------------
n_samples <- ncol(data_raw)
min_samples_expr <- max(3, ceiling(0.25 * n_samples))

expr_filter <- rowSums(data_raw >= min_expr) >= min_samples_expr
var_filter  <- var_vec >= min_var_cutoff
keep        <- expr_filter & var_filter

cat("\n🔎 Genes passing filters:\n",
    "  By expression (>= min_expr in ≥", min_samples_expr, " samples): ", sum(expr_filter), "\n",
    "  By variance (", cutoff_desc, "): ", sum(var_filter), "\n",
    "  Passing BOTH: ", sum(keep), "\n\n", sep="")
if (sum(keep) < 3) stop("❌ Too few genes after filtering (need ≥3). Relax cutoff and retry.", call. = FALSE)

data_filtered <- data_raw[keep, , drop = FALSE]
cat("✅ Retained ", nrow(data_filtered), " genes after filtering\n", sep = "")

# ----------------------------- Optional target restriction -----------------------------
genes     <- rownames(data_filtered)
all_genes <- rownames(data_raw)

cat("\n📄 OPTIONAL: Select a file containing target genes (one per line, gene names in the FIRST column).\n")
cat("  👉 Press 'Cancel' to skip and use ALL genes.\n")

target_file  <- tryCatch(file.choose(), error = function(e) NA)
target_genes <- character(0)
target_raw   <- character(0)

if (!is.na(target_file)) {
  cat("📥 Selected target gene file: ", target_file, "\n", sep = "")
  
  tmp <- tryCatch(
    read.csv(target_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(tmp) || nrow(tmp) == 0) {
    cat("⚠️ Target file contained no readable entries.\n")
  } else {
    target_raw <- tmp[[1]]
    target_raw <- unique(trimws(as.character(target_raw)))
    target_raw <- target_raw[nzchar(target_raw)]
    
    cat("📌 Loaded ", length(target_raw),
        " unique raw target entries from column 1.\n", sep = "")
    
    if (length(target_raw) > 0) {
      in_matrix     <- target_raw[target_raw %in% all_genes]
      not_in_matrix <- setdiff(target_raw, all_genes)
      
      if (length(in_matrix) > 0) {
        diag_df <- data.frame(
          gene      = in_matrix,
          in_matrix = TRUE,
          pass_expr = expr_filter[in_matrix],
          pass_var  = var_filter[in_matrix],
          pass_both = keep[in_matrix],
          stringsAsFactors = FALSE
        )
        
        diag_df$reason <- ifelse(
          diag_df$pass_both, "kept_after_filters",
          ifelse(!diag_df$pass_expr & !diag_df$pass_var, "failed_expr_and_var",
                 ifelse(!diag_df$pass_expr, "failed_expr_only",
                        ifelse(!diag_df$pass_var, "failed_var_only", "other")))
        )
        
        n_in_file    <- length(target_raw)
        n_in_matrix  <- nrow(diag_df)
        n_not_matrix <- length(not_in_matrix)
        n_kept       <- sum(diag_df$pass_both, na.rm = TRUE)
        n_fail_expr  <- sum(!diag_df$pass_expr &  diag_df$pass_var, na.rm = TRUE)
        n_fail_var   <- sum( diag_df$pass_expr & !diag_df$pass_var, na.rm = TRUE)
        n_fail_both  <- sum(!diag_df$pass_expr & !diag_df$pass_var, na.rm = TRUE)
        
        cat("\n📊 Target gene diagnostics:\n",
            "  Total entries in file (non-empty): ", n_in_file, "\n",
            "  Found in expression matrix:        ", n_in_matrix, "\n",
            "  Not found in matrix:               ", n_not_matrix, "\n",
            "  -------- Among those in matrix --------\n",
            "  Kept (pass expr & var filters):     ", n_kept, "\n",
            "  Failed expr only:                   ", n_fail_expr, "\n",
            "  Failed variance only:               ", n_fail_var, "\n",
            "  Failed both expr & var:             ", n_fail_both, "\n\n",
            sep = "")
        
        diag_out <- diag_df
        if (length(not_in_matrix) > 0) {
          diag_out <- rbind(
            diag_out,
            data.frame(
              gene      = not_in_matrix,
              in_matrix = FALSE,
              pass_expr = NA,
              pass_var  = NA,
              pass_both = NA,
              reason    = "not_in_expression_matrix",
              stringsAsFactors = FALSE
            )
          )
        }
        
        write.csv(diag_out, file.path(output_dir, "target_gene_filter_diagnostics.csv"), row.names = FALSE)
        cat("💾 Wrote target gene diagnostics to target_gene_filter_diagnostics.csv\n")
        
        target_genes <- diag_df$gene[diag_df$pass_both]
        
        if (length(target_genes) == 0) {
          cat("⚠️ None of the target genes passed BOTH expression and variance filters.\n",
              "   ➜ No target restriction will be applied to triplets.\n", sep = "")
        } else {
          cat("✅ ", length(target_genes),
              " target gene(s) survived filters and can be used in triplets.\n", sep = "")
        }
      } else {
        cat("⚠️ None of the genes in the target file match any rownames in the expression matrix.\n")
        if (length(not_in_matrix) > 0) {
          write.csv(
            data.frame(gene = not_in_matrix, reason = "not_in_expression_matrix", stringsAsFactors = FALSE),
            file.path(output_dir, "target_gene_filter_diagnostics.csv"),
            row.names = FALSE
          )
          cat("💾 Wrote list of unmatched target genes to target_gene_filter_diagnostics.csv\n")
        }
      }
    } else {
      cat("⚠️ Target file had only blank/whitespace entries in column 1.\n")
    }
  }
} else {
  cat("ℹ️ No target file selected (Cancel pressed). Using ALL genes.\n")
}

if (length(target_genes) > 0) {
  genes_for_triplets <- target_genes
  cat("🔬 Using ", length(genes_for_triplets),
      " target gene(s) that passed filters to build triplets.\n", sep = "")
} else {
  genes_for_triplets <- genes
  cat("ℹ️ No target gene restriction after filters; using all ",
      length(genes_for_triplets), " filtered genes to build triplets.\n", sep = "")
}

# ----------------------------- Build all triplets -----------------------------
triplet_list <- combn(genes_for_triplets, 3, simplify = FALSE)
if (length(triplet_list) == 0) stop("❌ No triplets remain after applying target gene / count restrictions.", call. = FALSE)

# ----------------------------- Triplet analysis function (3 orientations) -----------------------------
analyze_triplet_3orient <- function(tri) {
  do_one <- function(A, B, C) {
    A_vals <- as.numeric(data_filtered[A, ]) + pseudocount_delta
    B_vals <- as.numeric(data_filtered[B, ]) + pseudocount_delta
    C_vals <- as.numeric(data_filtered[C, ]) + pseudocount_delta
    
    df <- data.frame(
      A     = log(A_vals),
      BDivC = log(B_vals / C_vals),
      group = group_factor
    )
    if (sd(df$BDivC) < 1e-8) return(NULL)
    
    fit <- tryCatch(lm(A ~ BDivC * group, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    
    s  <- summary(fit)
    cf <- coef(fit)
    
    g_levels <- levels(df$group)
    if (length(g_levels) != 2) return(NULL)
    baseline <- g_levels[1]
    other    <- g_levels[2]
    iname    <- paste0("BDivC:group", other)
    
    smc <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    p_int <- if (!is.null(smc) && (iname %in% rownames(smc))) smc[iname, 4] else NA_real_
    
    slope_base  <- if ("BDivC" %in% names(cf)) cf["BDivC"] else NA_real_
    slope_other <- if (!is.na(p_int) && (iname %in% names(cf))) slope_base + cf[iname] else slope_base
    
    perm_fun <- function() {
      dfp <- df
      dfp$group <- sample(dfp$group)
      m <- tryCatch(lm(A ~ BDivC * group, data = dfp), error = function(e) NULL)
      if (is.null(m)) return(NA_real_)
      sm <- tryCatch(coef(summary(m)), error = function(e) NULL)
      if (is.null(sm) || !(iname %in% rownames(sm))) return(NA_real_)
      sm[iname, 4]
    }
    
    perm_p <- NA_real_
    if (is.finite(p_int)) {
      pvals <- replicate(num_permutations, perm_fun())
      perm_p_val_count <- sum(pvals <= p_int, na.rm = TRUE)
      perm_p <- (perm_p_val_count + 1) / (num_permutations + 1)
    }
    
    data.frame(
      A = A, B = B, C = C,
      baseline_level = baseline,
      other_level    = other,
      slope_baseline = as.numeric(round(slope_base, 3)),
      slope_other    = as.numeric(round(slope_other, 3)),
      abs_slope      = as.numeric(round(abs(slope_base), 3)),
      r_squared      = as.numeric(round(s$r.squared, 3)),
      p_interaction  = p_int,
      perm_p_interaction = perm_p,
      stringsAsFactors = FALSE
    )
  }
  
  rbind(
    do_one(tri[1], tri[2], tri[3]),
    do_one(tri[2], tri[1], tri[3]),
    do_one(tri[3], tri[1], tri[2])
  )
}

# ----------------------------- Parallel plan -----------------------------
cores <- max(1, parallel::detectCores() - 1)
future::plan(multisession, workers = cores)
cat("\n🔧 Parallel plan set: multisession with ", cores, " workers\n", sep = "")

cat("\n⏳ Fitting triplets (3 orientations each)...\n")
results_list <- progressr::with_progress({
  furrr::future_map(
    triplet_list,
    analyze_triplet_3orient,
    .progress = TRUE,
    .options  = furrr::furrr_options(seed = TRUE)
  )
})

results_df <- do.call(rbind, Filter(Negate(is.null), results_list))
if (is.null(results_df) || nrow(results_df) == 0) stop("❌ No valid triplets fit.", call. = FALSE)

# =========================================================
# Triplet-level multiple testing control
# Collapse 3 orientations -> 1 p per unordered triplet, then BH
# Uses orientation-level perm_p_interaction already computed.
# =========================================================

# Helper: stable unordered triplet key
make_triplet_key <- function(A, B, C) {
  paste(sort(c(A, B, C)), collapse = "|")
}

# Attach unordered triplet key to each orientation row
results_df$triplet_key <- mapply(make_triplet_key, results_df$A, results_df$B, results_df$C)

# For each unordered triplet, combine the 3 orientation p-values into one:
# p_comb = 1 - (1 - p_min)^3   (Sidak correction for min-p across 3 tests)
combine_orient_p <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
  if (length(pvals) == 0) return(NA_real_)
  pmin <- min(pvals)
  1 - (1 - pmin)^3
}

# Compute triplet-level combined p and record which orientation "won" (smallest p)
triplet_keys <- unique(results_df$triplet_key)

triplet_summary <- lapply(triplet_keys, function(k) {
  idx <- which(results_df$triplet_key == k)
  pv  <- results_df$perm_p_interaction[idx]
  
  p_comb <- combine_orient_p(pv)
  
  # identify best orientation (min perm p) for reporting
  best_i <- idx[which.min(ifelse(is.finite(pv), pv, Inf))]
  best_p <- if (length(best_i) == 1 && is.finite(results_df$perm_p_interaction[best_i])) results_df$perm_p_interaction[best_i] else NA_real_
  
  # store genes (unordered) + best orientation (ordered A,B,C from that row)
  genes_unordered <- strsplit(k, "\\|", fixed = FALSE)[[1]]
  
  data.frame(
    triplet_key = k,
    g1 = genes_unordered[1],
    g2 = genes_unordered[2],
    g3 = genes_unordered[3],
    p_min_orient = best_p,
    p_triplet_combined = p_comb,
    best_A = results_df$A[best_i],
    best_B = results_df$B[best_i],
    best_C = results_df$C[best_i],
    best_r_squared = results_df$r_squared[best_i],
    stringsAsFactors = FALSE
  )
})

triplet_df <- do.call(rbind, triplet_summary)

# BH across triplets (one hypothesis per unordered triplet)
triplet_df$q_triplet_BH <- NA_real_
ok <- is.finite(triplet_df$p_triplet_combined)
if (sum(ok) > 0) {
  triplet_df$q_triplet_BH[ok] <- p.adjust(triplet_df$p_triplet_combined[ok], method = "BH")
}

# Join triplet-level p/q back onto every orientation row (so you can still bootstrap per orientation if you want)
results_df$p_triplet_combined <- triplet_df$p_triplet_combined[match(results_df$triplet_key, triplet_df$triplet_key)]
results_df$q_triplet_BH       <- triplet_df$q_triplet_BH[match(results_df$triplet_key, triplet_df$triplet_key)]

# Save triplet-level multiple-testing table
write.csv(triplet_df, file.path(output_dir, "triplet_level_results.csv"), row.names = FALSE)

cat("\n✅ Triplet-level collapse complete:\n",
    "  Unordered triplets: ", nrow(triplet_df), "\n",
    "  Finite combined p:  ", sum(is.finite(triplet_df$p_triplet_combined)), "\n",
    "  BH q < 0.10:        ", sum(triplet_df$q_triplet_BH < 0.10, na.rm = TRUE), "\n\n",
    sep = "")

cat("\n🔎 Filtering by triplet-level BH q-value...\n")

if (use_r2) {
  sig_results <- subset(results_df,
                        is.finite(q_triplet_BH) & q_triplet_BH < 0.10 &
                          is.finite(r_squared) & r_squared >= min_r2)
} else {
  sig_results <- subset(results_df,
                        is.finite(q_triplet_BH) & q_triplet_BH < 0.10)
}

# Optional: if you only want ONE row per unordered triplet in significant file,
# keep only the best orientation row per triplet:
sig_results <- sig_results[order(sig_results$triplet_key, sig_results$perm_p_interaction), ]
sig_results <- sig_results[!duplicated(sig_results$triplet_key), ]
# ----------------------------- Save primary outputs -----------------------------
write.csv(results_df,  file.path(output_dir, "all_triplets.csv"), row.names = FALSE)
write.csv(sig_results, file.path(output_dir, "significant_triplets.csv"), row.names = FALSE)

cat("\n📊 Summary:\n",
    "  Total tested rows (3x orientations): ", nrow(results_df), "\n",
    "  Permutation P-value < ", p_val_thresh, ": ",
    sum(is.finite(results_df$perm_p_interaction) & results_df$perm_p_interaction < p_val_thresh), "\n",
    if (use_r2) paste0("  Perm P < ", p_val_thresh, " & R^2 >= ", min_r2, ": ", nrow(sig_results), "\n") else "",
    sep="")

# ----------------------------- Permutation p-value histogram -----------------------------
valid_p <- results_df$perm_p_interaction[is.finite(results_df$perm_p_interaction)]
if (length(valid_p) > 0) {
  p_hist <- ggplot(data.frame(p = valid_p), aes(x = p)) +
    geom_histogram(binwidth = 0.01, color = "black", fill = "lightblue") +
    geom_vline(xintercept = p_val_thresh, linetype = "dashed", linewidth = 1, color = "red") +
    labs(title = "Permutation P-value Distribution", x = "Permutation P-value", y = "Count") +
    theme_minimal(base_size = 14)
  ggsave(file.path(output_dir, "perm_pval_distribution.png"), p_hist, width = 6, height = 4, dpi = 150)
}

# ----------------------------- Bootstrapping -----------------------------
if (nrow(results_df) > 0) {
  n_top <- suppressWarnings(as.integer(
    readline("How many top triplets (by permutation p-value) to bootstrap & plot? [e.g., 50]: ")
  ))
  if (!is.finite(n_top) || n_top < 1) n_top <- min(50L, nrow(results_df))
  
  ord_ok <- which(is.finite(results_df$perm_p_interaction))
  top_boot <- head(results_df[ord_ok, , drop = FALSE][order(results_df$perm_p_interaction[ord_ok]), , drop = FALSE], n_top)
  
  if (nrow(top_boot) == 0) {
    cat("❌ No triplets with finite permutation p-values for bootstrapping.\n")
  } else {
    
    # ---- UPDATED: removed purrr dependency (base R only in this function) ----
    bootstrap_one_triplet <- function(triplet_info, num_bootstrap, data_filtered, group_factor, pseudocount_delta, p_val_thresh_local, top_boot_df) {
      A <- triplet_info[["A"]]; B <- triplet_info[["B"]]; C <- triplet_info[["C"]]
      triplet_name <- paste(A, B, C, sep = "_")
      
      A_vals <- as.numeric(data_filtered[A, ]) + pseudocount_delta
      B_vals <- as.numeric(data_filtered[B, ]) + pseudocount_delta
      C_vals <- as.numeric(data_filtered[C, ]) + pseudocount_delta
      df <- data.frame(
        A     = log(A_vals),
        BDivC = log(B_vals / C_vals),
        group = droplevels(group_factor)
      )
      
      if (length(levels(df$group)) != 2) {
        return(list(
          name = triplet_name, slopes = NULL,
          summary = data.frame(
            A = A, B = B, C = C, baseline_level = NA, other_level = NA,
            bootstrap_diff = NA, ci_lower = NA, ci_upper = NA,
            agree_with_p_val = "no_two_groups", stringsAsFactors = FALSE
          ),
          log = data.frame(triplet = triplet_name, reason = "not_two_groups", stringsAsFactors = FALSE)
        ))
      }
      
      g1 <- levels(df$group)[1]; g2 <- levels(df$group)[2]
      iname <- paste0("BDivC:group", g2)
      
      if (sum(df$group == g1) < 2 || sum(df$group == g2) < 2) {
        return(list(
          name = triplet_name, slopes = NULL,
          summary = data.frame(
            A = A, B = B, C = C, baseline_level = g1, other_level = g2,
            bootstrap_diff = NA, ci_lower = NA, ci_upper = NA,
            agree_with_p_val = "too_few_per_group", stringsAsFactors = FALSE
          ),
          log = data.frame(triplet = triplet_name, reason = "too_few_per_group", stringsAsFactors = FALSE)
        ))
      }
      
      boot_one <- function(j) {
        idx1 <- sample(which(df$group == g1), replace = TRUE)
        idx2 <- sample(which(df$group == g2), replace = TRUE)
        d <- rbind(df[idx1, , drop = FALSE], df[idx2, , drop = FALSE])
        
        m <- tryCatch(lm(A ~ BDivC * group, data = d), error = function(e) NULL)
        if (is.null(m)) return(c(slope_1 = NA_real_, slope_2 = NA_real_))
        cf <- coef(m)
        
        if (!("BDivC" %in% names(cf))) return(c(slope_1 = NA_real_, slope_2 = NA_real_))
        s1 <- cf["BDivC"]
        s2 <- if (iname %in% names(cf)) s1 + cf[iname] else s1
        c(slope_1 = as.numeric(s1), slope_2 = as.numeric(s2))
      }
      
      boot_mat <- do.call(rbind, lapply(seq_len(num_bootstrap), boot_one))
      boot_df <- data.frame(slope_1 = boot_mat[, "slope_1"], slope_2 = boot_mat[, "slope_2"])
      boot_df <- boot_df[is.finite(boot_df$slope_1) & is.finite(boot_df$slope_2), , drop = FALSE]
      
      if (nrow(boot_df) < 10) {
        return(list(
          name = triplet_name, slopes = NULL,
          summary = data.frame(
            A = A, B = B, C = C, baseline_level = g1, other_level = g2,
            bootstrap_diff = NA, ci_lower = NA, ci_upper = NA,
            agree_with_p_val = "insufficient_valid_boots", stringsAsFactors = FALSE
          ),
          log = data.frame(triplet = triplet_name, reason = "insufficient_valid_boots", stringsAsFactors = FALSE)
        ))
      }
      
      slope_diff <- boot_df$slope_2 - boot_df$slope_1
      ci <- stats::quantile(slope_diff, c(0.025, 0.975), na.rm = TRUE)
      
      hit <- which(top_boot_df$A == A & top_boot_df$B == B & top_boot_df$C == C)
      p_here <- if (length(hit) >= 1) top_boot_df$perm_p_interaction[hit[1]] else NA_real_
      
      if (!is.finite(p_here)) {
        agree <- "no_p_val"
      } else {
        strong_boot <- (ci[1] > 0 || ci[2] < 0)
        null_boot   <- (ci[1] <= 0 && ci[2] >= 0)
        
        if (p_here < p_val_thresh_local && strong_boot) {
          agree <- "yes"
        } else if (p_here >= p_val_thresh_local && null_boot) {
          agree <- "yes"
        } else {
          agree <- "no"
        }
      }
      
      list(
        name   = triplet_name,
        slopes = list(early = boot_df$slope_1, late = boot_df$slope_2, diff = slope_diff),
        summary = data.frame(
          A = A, B = B, C = C,
          baseline_level = g1, other_level = g2,
          bootstrap_diff = mean(slope_diff, na.rm = TRUE),
          ci_lower = unname(ci[1]), ci_upper = unname(ci[2]),
          agree_with_p_val = agree,
          stringsAsFactors = FALSE
        ),
        log = data.frame(triplet = triplet_name, reason = "ok", stringsAsFactors = FALSE)
      )
    }
    
    cat("\n🔁 Bootstrapping ", nrow(top_boot), " triplets using ", future::nbrOfWorkers(), " workers…\n", sep="")
    cat("🔍 Before bootstrap future_map:\n")
    print(future::plan())
    cat("🔢 Workers (bootstrap): ", future::nbrOfWorkers(), "\n\n", sep = "")
    
    boot_res <- progressr::with_progress({
      furrr::future_map(
        seq_len(nrow(top_boot)),
        ~ bootstrap_one_triplet(
          list(A = top_boot$A[.x], B = top_boot$B[.x], C = top_boot$C[.x]),
          num_bootstrap, data_filtered, group_factor, pseudocount_delta, p_val_thresh, top_boot
        ),
        .progress = TRUE,
        .options = furrr::furrr_options(seed = TRUE)
      )
    })
    
    boot_logs <- do.call(rbind, lapply(boot_res, `[[`, "log"))
    write.csv(boot_logs, file.path(output_dir, "bootstrap_log.csv"), row.names = FALSE)
    
    ok_res <- Filter(function(x) !is.null(x$slopes), boot_res)
    
    if (length(ok_res) == 0) {
      cat("⚠️ No bootstrap sets reached ≥10 valid resamples. See bootstrap_log.csv for reasons.\n")
      write.csv(data.frame(), file.path(output_dir, "bootstrap_summary_top.csv"), row.names = FALSE)
    } else {
      all_boots <- setNames(lapply(ok_res, `[[`, "slopes"),
                            vapply(ok_res, `[[`, "name", FUN.VALUE = character(1)))
      saveRDS(all_boots, file = file.path(output_dir, "bootstrap_results_all.rds"))
      
      boot_sum <- do.call(rbind, lapply(ok_res, `[[`, "summary"))
      write.csv(boot_sum, file.path(output_dir, "bootstrap_summary_top.csv"), row.names = FALSE)
      
      cat("🎨 Saving scatter + density plots…\n")
      for (nm in names(all_boots)) {
        A_B_C <- unlist(strsplit(nm, "_", fixed = TRUE))
        A <- A_B_C[1]; B <- A_B_C[2]; C <- A_B_C[3]
        if (!(A %in% rownames(data_filtered) &&
              B %in% rownames(data_filtered) &&
              C %in% rownames(data_filtered))) next
        
        A_vals <- as.numeric(data_filtered[A, ]) + pseudocount_delta
        B_vals <- as.numeric(data_filtered[B, ]) + pseudocount_delta
        C_vals <- as.numeric(data_filtered[C, ]) + pseudocount_delta
        dfp <- data.frame(
          A     = log(A_vals),
          BDivC = log(B_vals / C_vals),
          group = group_factor
        )
        
        boot <- all_boots[[nm]]
        p1 <- ggplot(dfp, aes(BDivC, A, color = group)) +
          geom_point(size = 2, alpha = 0.7) +
          geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
          labs(title = paste0("Triplet: ", A, " ~ log(", B, "/", C, ")"),
               x = paste0("log(", B, "/", C, ")"),
               y = paste0("log(", A, ")")) +
          theme_minimal(base_size = 14)
        
        g_levels <- levels(group_factor)
        dens_df <- data.frame(
          Slope = c(boot$early, boot$late),
          Group = factor(c(rep(g_levels[1], length(boot$early)),
                           rep(g_levels[2], length(boot$late))),
                         levels = g_levels)
        )
        p2 <- ggplot(dens_df, aes(Slope, fill = Group)) +
          geom_density(alpha = 0.5) +
          labs(title = "Bootstrapped slope distributions") +
          theme_minimal(base_size = 14)
        
        g <- gridExtra::grid.arrange(p1, p2, ncol = 2)
        ggsave(file.path(plot_dir, paste0("triplet_", nm, ".png")), g, width = 12, height = 5, dpi = 150)
        grid::grid.newpage()
      }
    }
  }
}

# ----------------------------- README.qmd (UPDATED: add script/git provenance lines) -----------------------------
readme <- file.path(output_dir, "README.qmd")
sink(readme)
cat(
  "---\n",
  paste0("title: \"Triplet Interaction Pipeline Summary (", input_name, ")\""), "\n",
  "format: html\n",
  "---\n\n",
  "## 🧾 Run details\n",
  "- Input file: ", normalizePath(file_path, winslash="/", mustWork=FALSE), "\n",
  "- Output dir: ", normalizePath(output_dir, winslash="/", mustWork=FALSE), "\n",
  "- Date: ", as.character(Sys.Date()), "\n",
  "- Script name: ", script_name, "\n",
  "- Script path: ", script_path, "\n",
  "- Script full: ", script_full, "\n",
  "- repo_root: ", git_meta$repo_root, "\n",
  "- git_branch: ", git_meta$git_branch, "\n",
  "- git_commit: ", git_meta$git_commit, "\n",
  "- git_remote: ", git_meta$git_remote, "\n\n",
  "## ⚙️ Parameters\n",
  "- min_expr: ", min_expr, "\n",
  "- pseudocount_delta: ", pseudocount_delta, "\n",
  "- num_permutations: ", num_permutations, "\n",
  "- num_bootstrap: ", num_bootstrap, "\n",
  "- permutation_p_val_threshold: ", p_val_thresh, "\n",
  if (use_r2) paste0("- min_r2: ", min_r2, "\n") else "- min_r2: (none)\n",
  "\n",
  "## 📁 Outputs\n",
  "- all_triplets.csv\n",
  "- significant_triplets.csv\n",
  "- perm_pval_distribution.png\n",
  "- bootstrap_results_all.rds (if bootstraps run)\n",
  "- bootstrap_summary_top.csv (if bootstraps run)\n",
  "- bootstrap_log.csv (if bootstraps run)\n",
  "- final_triplet_plots/*.png (if bootstraps run)\n",
  "- OUTPUT_INVENTORY.txt\n",
  "\n",
  "## 🧪 Session Info\n",
  "```{r}\n",
  "sessionInfo()\n",
  "```\n"
)
sink()
cat("📄 Wrote README.qmd\n")

# ----------------------------- OUTPUT_INVENTORY.txt (MANDATORY) -----------------------------
write_output_inventory(output_dir)

cat("\n✅ Done.\n")
cat("Output directory:\n  ", normalizePath(output_dir, winslash="/", mustWork=FALSE), "\n", sep = "")