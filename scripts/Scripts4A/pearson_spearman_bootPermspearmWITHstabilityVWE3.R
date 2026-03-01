#!/usr/bin/env Rscript
###############################################################################
# Correlation Analysis with Bootstrap Confidence Intervals, FDR Correction,
# and Bootstrap-Based Stability Metrics
#
# UPDATE (requested):
#   1) Add option to run RAW (0) or LOG2 + pseudocount (1) and execute all
#      downstream analyses on the chosen scale.
#   2) Write day-wise correlation matrices (Pearson + Spearman) into separate
#      subdirectories (plus PearsonZ kept, also separated).
#input format WIDE_Liver5SDTranscriptome4PEARSON.csv
###############################################################################

suppressPackageStartupMessages({
  library(parallel)
  library(ggplot2)
  library(pbapply)
})

###############################################################################
# 🛠 Helper: script path + sanitation
###############################################################################
.get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg)))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
      if (nzchar(p)) return(normalizePath(p))
    }
  }
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile))
  NA_character_
}

.sanitize <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

.select_output_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
      d <- tryCatch(rstudioapi::selectDirectory(caption = "Select output directory"), error = function(e) NULL)
      if (!is.null(d) && nzchar(d)) return(normalizePath(d))
    }
  }
  cat("\nRStudio directory chooser not available. Please pick any file INSIDE the folder you want outputs written to.\n")
  dirname(file.choose())
}

###############################################################################
# 🛠 User Inputs
###############################################################################
cat("\nChoose analysis scale:\n")
cat("  0 = RAW\n")
cat("  1 = LOG2(x + pseudocount)\n")
cat("Enter choice [default = 0]: ")
scale_choice <- suppressWarnings(as.integer(readLines(con = stdin(), n = 1)))
if (is.na(scale_choice) || !scale_choice %in% c(0, 1)) scale_choice <- 0

pseudocount <- 1
if (scale_choice == 1) {
  cat("\nEnter pseudocount to use in LOG2(x + pseudocount) [default = 1]: ")
  pc_in <- suppressWarnings(as.numeric(readLines(con = stdin(), n = 1)))
  if (!is.na(pc_in) && is.finite(pc_in) && pc_in > 0) pseudocount <- pc_in
}

analysis_scale_tag <- if (scale_choice == 0) "RAW" else paste0("LOG2P", .sanitize(as.character(pseudocount)))
cat("\n✅ Using scale: ", analysis_scale_tag, "\n", sep = "")

cat("\nEnter number of bootstrap replicates [default = 1000]: ")
boot_R <- suppressWarnings(as.integer(readLines(con = stdin(), n = 1)))
if (is.na(boot_R) || boot_R <= 0) boot_R <- 1000

max_cores <- parallel::detectCores()
cat("\nEnter number of processors to use [default = ", floor(max_cores / 2), "]: ", sep = "")
use_cores <- suppressWarnings(as.integer(readLines(con = stdin(), n = 1)))
if (is.na(use_cores) || use_cores <= 0) use_cores <- floor(max_cores / 2)
use_cores <- max(1L, min(use_cores, max_cores))

cat("\nGuidance for choosing |r| threshold based on per-day sample size:\n\n")
cat("  n < 10        -> |r| >= 0.40  (conservative; reduce noise)\n")
cat("  10 <= n <= 25 -> |r| >= 0.30  (recommended default)\n")
cat("  n > 25        -> |r| >= 0.20-0.30 (higher power; optional sensitivity check)\n\n")

cat("\nEnter |r| threshold for 'magnitude stability' [default = 0.30]: ")
stab_r_thr <- suppressWarnings(as.numeric(readLines(con = stdin(), n = 1)))
if (is.na(stab_r_thr) || stab_r_thr <= 0 || stab_r_thr >= 1) stab_r_thr <- 0.30
cat("✅ Using magnitude stability threshold |r| >= ", stab_r_thr, "\n", sep = "")

cat("\nEnter minimum pairwise complete observations per pair (min_n) [default = 4]: ")
min_n <- suppressWarnings(as.integer(readLines(con = stdin(), n = 1)))
if (is.na(min_n) || min_n < 3) min_n <- 4
cat("✅ Using min_n = ", min_n, "\n", sep = "")

###############################################################################
# 📂 Load Data + Output Directory Setup
###############################################################################
if (interactive()) {
  cat("\nSelect input CSV file:\n")
  input_file <- file.choose()
  if (!file.exists(input_file)) stop("Input file not found: ", input_file)
  
  mydata <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!("Day" %in% names(mydata))) {
    stop("Input must contain a column named 'Day'. Found columns: ", paste(names(mydata), collapse = ", "))
  }
  mydata$Day <- as.factor(mydata$Day)
  
  cat("\nSelect output directory:\n")
  base_dir <- .select_output_dir()
  
  cat("\nEnter a name prefix for outputs (or press Enter to use input file name): ")
  prefix_input <- readLines(con = stdin(), n = 1)
  
  input_name <- tools::file_path_sans_ext(basename(input_file))
  if (nzchar(prefix_input)) input_name <- prefix_input
  input_name <- .sanitize(input_name)
  
  timestamp       <- format(Sys.time(), "%Y%m%d_%H%M%S")
  timestamp_human <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Include scale tag to prevent overwriting and to make outputs self-identifying.
  output_dir <- file.path(base_dir, paste0(input_name, "_", analysis_scale_tag, "_", timestamp))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("✅ Output will be saved in: ", normalizePath(output_dir), "\n", sep = "")
  
} else {
  if (!exists("input_file") || is.null(input_file) || !nzchar(input_file)) {
    stop("No interactive session and no 'input_file' provided. Run interactively or set input_file before sourcing.")
  }
  if (!file.exists(input_file)) stop("input_file not found: ", input_file)
  
  mydata <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("Day" %in% names(mydata))) {
    stop("Input must contain a column named 'Day'. Found columns: ", paste(names(mydata), collapse = ", "))
  }
  mydata$Day <- as.factor(mydata$Day)
  
  input_name <- tools::file_path_sans_ext(basename(input_file))
  input_name <- .sanitize(input_name)
  
  timestamp       <- format(Sys.time(), "%Y%m%d_%H%M%S")
  timestamp_human <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  output_dir <- file.path(getwd(), paste0(input_name, "_", analysis_scale_tag, "_", timestamp))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("✅ Output will be saved in: ", normalizePath(output_dir), "\n", sep = "")
}

# Directories (create ONCE, early, for both modes)
expr_dir       <- file.path(output_dir, "expression_levels_by_day")
cor_by_day_dir <- file.path(output_dir, "correlations_by_day")

# Requested: separate subdirectories for Pearson and Spearman day-wise correlation matrices
corr_mat_root      <- file.path(output_dir, "correlation_matrices_by_day")
corr_mat_pearson   <- file.path(corr_mat_root, "Pearson")
corr_mat_spearman  <- file.path(corr_mat_root, "Spearman")
corr_mat_pearson_z <- file.path(corr_mat_root, "PearsonZ")

dir.create(expr_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cor_by_day_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(corr_mat_root, showWarnings = FALSE, recursive = TRUE)
dir.create(corr_mat_pearson, showWarnings = FALSE, recursive = TRUE)
dir.create(corr_mat_spearman, showWarnings = FALSE, recursive = TRUE)
dir.create(corr_mat_pearson_z, showWarnings = FALSE, recursive = TRUE)

script_path <- .get_script_path()
script_name <- if (!is.na(script_path)) tools::file_path_sans_ext(basename(script_path)) else "Rscript"
script_name <- .sanitize(script_name)

assign("script_path", script_path, envir = .GlobalEnv)
assign("script_name", script_name, envir = .GlobalEnv)
assign("input_file",  input_file,  envir = .GlobalEnv)
assign("input_name",  input_name,  envir = .GlobalEnv)
assign("output_dir",  output_dir,  envir = .GlobalEnv)
assign("timestamp",   timestamp,   envir = .GlobalEnv)
assign("timestamp_human", timestamp_human, envir = .GlobalEnv)

assign("analysis_scale_tag", analysis_scale_tag, envir = .GlobalEnv)
assign("scale_choice", scale_choice, envir = .GlobalEnv)
assign("pseudocount", pseudocount, envir = .GlobalEnv)

assign("corr_mat_root", corr_mat_root, envir = .GlobalEnv)
assign("corr_mat_pearson", corr_mat_pearson, envir = .GlobalEnv)
assign("corr_mat_spearman", corr_mat_spearman, envir = .GlobalEnv)
assign("corr_mat_pearson_z", corr_mat_pearson_z, envir = .GlobalEnv)

###############################################################################
# Force all columns except Day to be numeric
###############################################################################
num_vars <- setdiff(names(mydata), "Day")
if (!length(num_vars)) stop("No analyte columns found (only 'Day' present).")

mydata[num_vars] <- lapply(mydata[num_vars], function(x) suppressWarnings(as.numeric(x)))

###############################################################################
# Apply chosen scale (RAW vs LOG2 + pseudocount)
###############################################################################
if (scale_choice == 1) {
  # Guardrails: log2(x + pseudocount) requires x + pseudocount > 0 for finite values.
  min_val <- suppressWarnings(min(unlist(mydata[num_vars]), na.rm = TRUE))
  if (is.finite(min_val) && (min_val + pseudocount) <= 0) {
    stop(
      "LOG2 transform invalid: min(x) + pseudocount <= 0.\n",
      "min(x)=", min_val, ", pseudocount=", pseudocount, ".\n",
      "Use RAW mode or choose a larger pseudocount (or fix negative values upstream)."
    )
  }
  mydata[num_vars] <- lapply(mydata[num_vars], function(x) {
    ifelse(is.na(x), NA_real_, log2(x + pseudocount))
  })
}

cat("\n📌 Unique non-NA values per variable by Day (diagnostic; post-scale):\n")
tab <- by(mydata, mydata$Day, function(df) {
  sapply(df[, setdiff(colnames(df), "Day"), drop = FALSE],
         function(col) length(unique(col[!is.na(col)])))
})
print(tab)

###############################################################################
# Variance filtering and plots (post-scale)
###############################################################################
var_vec <- sapply(mydata[num_vars], function(x) var(x, na.rm = TRUE))

default_cutoff    <- 1e-6
percentile_cutoff <- suppressWarnings(quantile(var_vec, 0.10, na.rm = TRUE))

grid_cutoffs <- 10^(seq(floor(log10(default_cutoff)), 1, by = 1))
cutoffs_df <- rbind(
  data.frame(Cutoff = default_cutoff,                Source = "Default"),
  data.frame(Cutoff = as.numeric(percentile_cutoff), Source = "10th percentile"),
  data.frame(Cutoff = grid_cutoffs,                  Source = "Grid")
)

cutoffs_df <- cutoffs_df[order(cutoffs_df$Cutoff), , drop = FALSE]
cutoffs_df <- cutoffs_df[!duplicated(round(cutoffs_df$Cutoff, 12)), , drop = FALSE]

total_vars <- sum(!is.na(var_vec))
cutoffs_df$Kept         <- sapply(cutoffs_df$Cutoff, function(th) sum(var_vec > th, na.rm = TRUE))
cutoffs_df$Dropped      <- total_vars - cutoffs_df$Kept
cutoffs_df$Prop_Kept    <- round(cutoffs_df$Kept / total_vars, 3)
cutoffs_df$Cutoff_Label <- format(cutoffs_df$Cutoff, scientific = TRUE)

cat("\n📌 Variance cutoff sweep:\n")
print(cutoffs_df)

out_sweep <- file.path(output_dir, "Variance_Filter_Sweep.csv")
write.csv(cutoffs_df, out_sweep, row.names = FALSE)
cat("Wrote variance sweep table to: ", out_sweep, "\n", sep = "")

p_sweep <- ggplot(cutoffs_df, aes(x = Cutoff, y = Kept, color = Source)) +
  geom_line() + geom_point() +
  scale_x_log10() +
  theme_minimal() +
  labs(title = paste0("Variables kept vs variance cutoff (", analysis_scale_tag, ")"),
       x = "Variance cutoff (log10 scale)", y = "# variables kept")
ggsave(file.path(output_dir, "Variance_Filter_Sweep.png"), p_sweep, width = 7, height = 4, dpi = 200)

png(file.path(output_dir, "Variance_Histogram.png"), width = 900, height = 650)
if (sum(var_vec > 0, na.rm = TRUE) > 0) {
  log_var_vec <- log10(var_vec[var_vec > 0])
  hist(log_var_vec, breaks = 100,
       main = paste0("log10(Feature Variance Distribution) (", analysis_scale_tag, ")"),
       xlab = "log10(Variance)", col = "lightblue", border = "white")
  abline(v = log10(default_cutoff),    col = "red",       lwd = 2, lty = 2)
  abline(v = log10(percentile_cutoff), col = "darkgreen", lwd = 2, lty = 2)
  legend("topright",
         legend = c(paste("Default (red) = ", signif(default_cutoff, 3)),
                    paste("10th percentile (green) = ", signif(percentile_cutoff, 3))),
         col = c("red", "darkgreen"), lty = 2, lwd = 2, bty = "n")
} else {
  hist(var_vec, breaks = 100, main = paste0("Feature Variance Distribution (", analysis_scale_tag, ")"),
       xlab = "Variance", col = "lightblue", border = "white")
}
dev.off()

cat("\nVariance cutoff options (choose by row):\n")
tbl_print <- transform(cutoffs_df, Row = seq_len(nrow(cutoffs_df)))
tbl_print <- tbl_print[, c("Row","Source","Cutoff","Cutoff_Label","Kept","Dropped","Prop_Kept")]
print(tbl_print, row.names = FALSE)

cat("\nEnter the row number of your chosen cutoff [default = 1]: ")
row_choice <- suppressWarnings(as.integer(readLines(con = stdin(), n = 1)))
if (is.na(row_choice) || row_choice < 1 || row_choice > nrow(cutoffs_df)) row_choice <- 1

# ------------------------------------------------------------
# SAFETY CONFIRMATION: feature retention check
# ------------------------------------------------------------
chosen_cutoff <- cutoffs_df$Cutoff[row_choice]
kept_n        <- cutoffs_df$Kept[row_choice]
total_n       <- total_vars

cat("\n⚠️  You selected:\n")
cat("   Row: ", row_choice, "\n", sep = "")
cat("   Source: ", cutoffs_df$Source[row_choice], "\n", sep = "")
cat("   Variance cutoff: ", signif(chosen_cutoff, 6), "\n", sep = "")
cat("   Features kept: ", kept_n, " / ", total_n, "\n\n", sep = "")

cat("Proceed with this cutoff? (y/n) [default = n]: ")
confirm <- tolower(trimws(readLines(con = stdin(), n = 1)))

if (!nzchar(confirm) || confirm != "y") {
  cat("\n❌ Cutoff not confirmed. Please re-run and choose a different row.\n")
  stop("Aborted by user at variance-cutoff confirmation.")
}


variance_threshold <- cutoffs_df$Cutoff[row_choice]
cat(sprintf("✅ Using variance cutoff = %.12f (row %d, Source: %s)\n",
            variance_threshold, row_choice, cutoffs_df$Source[row_choice]))

keep_vars <- names(var_vec)[var_vec > variance_threshold]
cat("Keeping ", length(keep_vars), " of ", length(num_vars),
    " variables after filtering low-variance features (variance > ", variance_threshold, ").\n", sep = "")

if (length(keep_vars) < 2) stop("After variance filtering, fewer than 2 variables remain. Lower the cutoff.")

mydata <- mydata[, c("Day", keep_vars), drop = FALSE]

###############################################################################
# Choose FDR settings
###############################################################################
cat("\nEnter FDR correction method (BH, BY, holm, bonferroni) [default = BH]: ")
user_fdr_method <- toupper(readLines(con = stdin(), n = 1))
if (!user_fdr_method %in% c("BH","BY","HOLM","BONFERRONI")) user_fdr_method <- "BH"

cat("\nEnter FDR significance cutoff [default = 0.05]: ")
user_fdr_cutoff <- suppressWarnings(as.numeric(readLines(con = stdin(), n = 1)))
if (is.na(user_fdr_cutoff) || user_fdr_cutoff <= 0) user_fdr_cutoff <- 0.05

assign("user_fdr_method", user_fdr_method, envir = .GlobalEnv)
assign("user_fdr_cutoff", user_fdr_cutoff, envir = .GlobalEnv)
cat("✅ Using FDR method: ", user_fdr_method, " with cutoff: ", user_fdr_cutoff, "\n", sep = "")

###############################################################################
# Save per-day expression (post-filtering; post-scale)
###############################################################################
split_groups <- split(mydata, mydata$Day)
for (group_name in names(split_groups)) {
  safe_group <- gsub("[^A-Za-z0-9_]", "_", group_name)
  out_file <- file.path(expr_dir, paste0("Group_", safe_group, ".csv"))
  write.csv(split_groups[[group_name]], out_file, row.names = FALSE)
  cat("Wrote expression data: ", out_file, "\n", sep = "")
}

###############################################################################
# 📊 Helper Functions
###############################################################################
correlation_matrix <- function(mat, method = "pearson") {
  cor(mat, method = method, use = "pairwise.complete.obs")
}

pairwise_n_matrix <- function(mat) {
  ok <- !is.na(mat)
  nmat <- t(ok) %*% ok
  nmat <- as.matrix(nmat)
  dimnames(nmat) <- list(colnames(mat), colnames(mat))
  nmat
}

p_from_r_matrix <- function(rmat, nmat) {
  p <- matrix(NA_real_, nrow = nrow(rmat), ncol = ncol(rmat), dimnames = dimnames(rmat))
  valid <- !is.na(rmat) & !is.na(nmat) & (nmat >= 3) & (abs(rmat) <= 1)
  if (any(valid)) {
    r <- rmat[valid]
    n <- nmat[valid]
    t_val <- r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
    p[valid] <- 2 * stats::pt(-abs(t_val), df = n - 2)
  }
  diag(p) <- NA_real_
  p
}

mask_by_min_n <- function(xmat, nmat, min_n) {
  xmat[nmat < min_n] <- NA_real_
  xmat
}

bootstrap_matrix <- function(mat, R = 1000, method = "pearson", n_cores = 2, r_threshold = 0.30) {
  n <- nrow(mat)
  if (n < 4) stop("bootstrap_matrix requires n >= 4 rows in mat.")
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(cl, varlist = c("mat", "method"), envir = environment())
  parallel::clusterEvalQ(cl, library(stats))
  
  boot_results <- parallel::parLapply(cl, 1:R, function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    cor(mat[idx, , drop = FALSE], method = method, use = "pairwise.complete.obs")
  })
  
  boot_array <- simplify2array(boot_results)
  
  lower <- apply(boot_array, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot_array, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
  med   <- apply(boot_array, c(1, 2), median,   na.rm = TRUE)
  width <- upper - lower
  
  prob_pos <- apply(boot_array > 0, c(1, 2), mean, na.rm = TRUE)
  prob_neg <- apply(boot_array < 0, c(1, 2), mean, na.rm = TRUE)
  
  directional_stability <- abs(prob_pos - prob_neg)
  magnitude_stability   <- apply(abs(boot_array) >= r_threshold, c(1, 2), mean, na.rm = TRUE)
  boot_p_sign           <- 2 * pmin(prob_pos, prob_neg)
  
  diag(lower) <- diag(upper) <- diag(med) <- diag(width) <- NA_real_
  diag(directional_stability) <- diag(magnitude_stability) <- diag(boot_p_sign) <- NA_real_
  
  list(lower = lower, upper = upper,
       median = med, width = width,
       dir_stab = directional_stability,
       mag_stab = magnitude_stability,
       boot_p_sign = boot_p_sign)
}

###############################################################################
# 🚀 Per-Day Analysis
###############################################################################
analyze_day <- function(day_val) {
  sub <- mydata[mydata$Day == day_val, , drop = FALSE]
  vars <- setdiff(colnames(sub), "Day")
  mat <- as.matrix(sub[, vars, drop = FALSE])
  
  if (nrow(mat) < 4) {
    cat("Skipping Day=", as.character(day_val), " (n<4)\n", sep = "")
    return(NULL)
  }
  
  safe_day <- gsub("[^A-Za-z0-9_]", "_", as.character(day_val))
  
  nmat <- pairwise_n_matrix(mat)
  
  pearson_cor  <- correlation_matrix(mat, "pearson")
  spearman_cor <- correlation_matrix(mat, "spearman")
  
  pearson_cor  <- mask_by_min_n(pearson_cor,  nmat, min_n)
  spearman_cor <- mask_by_min_n(spearman_cor, nmat, min_n)
  
  pearson_p  <- p_from_r_matrix(pearson_cor,  nmat)
  spearman_p <- p_from_r_matrix(spearman_cor, nmat)
  
  diag(pearson_cor)  <- 1
  diag(spearman_cor) <- 1
  
  # Requested: separate subdirectories for Pearson and Spearman matrices
  write.csv(
    pearson_cor,
    file.path(corr_mat_pearson, paste0("CorrMatrix_Pearson_Day", safe_day, ".csv")),
    row.names = TRUE, quote = FALSE
  )
  write.csv(
    spearman_cor,
    file.path(corr_mat_spearman, paste0("CorrMatrix_Spearman_Day", safe_day, ".csv")),
    row.names = TRUE, quote = FALSE
  )
  
  # Keep PearsonZ (z-scored features) as before, but separate directory
  mat_z <- scale(mat)
  pearson_cor_z <- correlation_matrix(mat_z, "pearson")
  pearson_cor_z <- mask_by_min_n(pearson_cor_z, nmat, min_n)
  diag(pearson_cor_z) <- 1
  
  write.csv(
    pearson_cor_z,
    file.path(corr_mat_pearson_z, paste0("CorrMatrix_PearsonZ_Day", safe_day, ".csv")),
    row.names = TRUE, quote = FALSE
  )
  
  boot_res <- bootstrap_matrix(mat, R = boot_R, method = "pearson",
                               n_cores = use_cores, r_threshold = stab_r_thr)
  
  out <- vector("list", length = 0)
  combn(vars, 2, function(v) {
    i <- match(v[1], vars)
    j <- match(v[2], vars)
    
    nij <- nmat[i, j]
    if (is.na(nij) || nij < min_n) {
      out[[length(out) + 1]] <<- data.frame(
        Pair = paste(v, collapse = " vs "),
        Metab1 = v[1], Metab2 = v[2],
        Day = as.numeric(as.character(day_val)),
        N_pair = as.integer(nij),
        Pearson_Correlation = NA_real_, Pearson_P_value = NA_real_,
        Spearman_Correlation = NA_real_, Spearman_P_value = NA_real_,
        Boot_LowerCI = NA_real_, Boot_UpperCI = NA_real_,
        Boot_Median = NA_real_, Boot_Width = NA_real_,
        Stability_Directional = NA_real_, Stability_Magnitude = NA_real_,
        Boot_P_Sign = NA_real_, CI_ExcludesZero = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      low <- boot_res$lower[i, j]
      up  <- boot_res$upper[i, j]
      out[[length(out) + 1]] <<- data.frame(
        Pair = paste(v, collapse = " vs "),
        Metab1 = v[1], Metab2 = v[2],
        Day = as.numeric(as.character(day_val)),
        N_pair = as.integer(nij),
        
        Pearson_Correlation  = pearson_cor[i, j],
        Pearson_P_value      = pearson_p[i, j],
        Spearman_Correlation = spearman_cor[i, j],
        Spearman_P_value     = spearman_p[i, j],
        
        Boot_LowerCI = low,
        Boot_UpperCI = up,
        Boot_Median  = boot_res$median[i, j],
        Boot_Width   = boot_res$width[i, j],
        
        Stability_Directional = boot_res$dir_stab[i, j],
        Stability_Magnitude   = boot_res$mag_stab[i, j],
        Boot_P_Sign           = boot_res$boot_p_sign[i, j],
        CI_ExcludesZero       = as.integer(!is.na(low) && !is.na(up) && (low * up > 0)),
        
        stringsAsFactors = FALSE
      )
    }
  })
  
  do.call(rbind, out)
}

days <- unique(mydata$Day)
cat("\n⏳ Running correlation analysis by Day...\n")

results_list <- pblapply(days, analyze_day)
results_list <- results_list[!sapply(results_list, is.null)]
if (length(results_list) == 0) stop("No results for any Day group. Check your data.")
cor_df <- do.call(rbind, results_list)

###############################################################################
# 🧮 FDR Corrections (WITHIN Day)
###############################################################################
fdr_methods <- c("BH", "BY", "holm", "bonferroni")
for (m in fdr_methods) {
  cor_df[[paste0("Pearson_FDR_", m)]] <- ave(
    cor_df$Pearson_P_value, cor_df$Day,
    FUN = function(p) if (all(is.na(p))) rep(NA_real_, length(p)) else p.adjust(p, method = m)
  )
  cor_df[[paste0("Spearman_FDR_", m)]] <- ave(
    cor_df$Spearman_P_value, cor_df$Day,
    FUN = function(p) if (all(is.na(p))) rep(NA_real_, length(p)) else p.adjust(p, method = m)
  )
}

stab_dir_min <- 0.70
stab_mag_min <- 0.70
fdr_col_suffix <- if (user_fdr_method %in% c("BH","BY")) user_fdr_method else tolower(user_fdr_method)
fdr_col_name   <- paste0("Pearson_FDR_", fdr_col_suffix)

if (!fdr_col_name %in% names(cor_df)) {
  stop("FDR column not found: ", fdr_col_name,
       ". Available columns: ", paste(names(cor_df), collapse = ", "))
}

fdr_vec <- cor_df[[fdr_col_name]]
cor_df$StableEdge <-
  (cor_df$CI_ExcludesZero == 1) &
  (cor_df$Stability_Directional >= stab_dir_min) &
  (cor_df$Stability_Magnitude   >= stab_mag_min) &
  !is.na(fdr_vec) & (fdr_vec < user_fdr_cutoff)

combined_file <- file.path(output_dir, "Correlation_by_Day_with_FDR.csv")
write.csv(cor_df, combined_file, row.names = FALSE)
cat("Wrote combined per-Day correlation data to: ", combined_file, "\n", sep = "")

###############################################################################
# 📊 Save by-Day Results
###############################################################################
cor_split <- split(cor_df, cor_df$Day)
for (day_name in names(cor_split)) {
  safe_day <- gsub("[^A-Za-z0-9_]", "_", day_name)
  out_file <- file.path(cor_by_day_dir, paste0("Correlation_by_Day_", safe_day, ".csv"))
  write.csv(cor_split[[day_name]], out_file, row.names = FALSE)
  cat("Wrote correlation data: ", out_file, "\n", sep = "")
}

###############################################################################
# 📈 Plot Top Changing Pairs
###############################################################################
sd_values <- tapply(cor_df$Pearson_Correlation, cor_df$Pair, function(x) sd(x, na.rm = TRUE))
sd_values <- sd_values[!is.na(sd_values)]
if (length(sd_values) == 0) stop("No valid SDs for correlations (all NA).")

top_pairs <- names(sort(sd_values, decreasing = TRUE))[1:min(10, length(sd_values))]
top_df <- cor_df[cor_df$Pair %in% top_pairs, , drop = FALSE]

if (nrow(top_df) > 0) {
  plot <- ggplot(top_df, aes(x = Day, y = Pearson_Correlation)) +
    geom_line(aes(group = Pair), linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Boot_LowerCI, ymax = Boot_UpperCI),
                  width = 0.2, color = "gray50") +
    facet_wrap(~ Pair, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(title = paste0("Top Changing Correlations (Pearson, 95% CI) [", analysis_scale_tag, "]"),
         x = "Day", y = "Pearson Correlation")
  ggsave(file.path(output_dir, "Top_Correlations_Faceted.png"), plot, width = 12, height = 7, dpi = 300)
} else {
  cat("No top pairs to plot.\n")
}

###############################################################################
# 🌍 Global Correlations
###############################################################################
vars <- setdiff(colnames(mydata), "Day")
mat_all <- as.matrix(mydata[, vars, drop = FALSE])

nmat_all <- pairwise_n_matrix(mat_all)
global_cor_pearson  <- correlation_matrix(mat_all, "pearson")
global_cor_spearman <- correlation_matrix(mat_all, "spearman")

global_cor_pearson  <- mask_by_min_n(global_cor_pearson,  nmat_all, min_n)
global_cor_spearman <- mask_by_min_n(global_cor_spearman, nmat_all, min_n)

global_p_pearson  <- p_from_r_matrix(global_cor_pearson,  nmat_all)
global_p_spearman <- p_from_r_matrix(global_cor_spearman, nmat_all)

global_results <- list()
combn(vars, 2, function(v) {
  i <- match(v[1], vars)
  j <- match(v[2], vars)
  
  pearson_r  <- global_cor_pearson[i, j]
  spearman_r <- global_cor_spearman[i, j]
  abs_diff   <- if (is.na(pearson_r) || is.na(spearman_r)) NA_real_ else abs(spearman_r - pearson_r)
  
  global_results[[length(global_results) + 1]] <<- data.frame(
    Pair       = paste(v, collapse = " vs "),
    Metab1     = v[1],
    Metab2     = v[2],
    N_pair     = as.integer(nmat_all[i, j]),
    Pearson    = pearson_r,
    Spearman   = spearman_r,
    AbsDiff    = abs_diff,
    P_Pearson  = global_p_pearson[i, j],
    P_Spearman = global_p_spearman[i, j],
    stringsAsFactors = FALSE
  )
})

global_results <- do.call(rbind, global_results)

for (m in fdr_methods) {
  global_results[[paste0("Pearson_FDR_", m)]]  <- p.adjust(global_results$P_Pearson,  method = m)
  global_results[[paste0("Spearman_FDR_", m)]] <- p.adjust(global_results$P_Spearman, method = m)
}

global_out_path <- file.path(output_dir, "Global_Correlations_with_FDR.csv")
write.csv(global_results, global_out_path, row.names = FALSE)
cat("Wrote Global Correlations to: ", global_out_path, "\n", sep = "")

cat("\n✅ Analysis complete. Results in: ", output_dir, "\n", sep = "")

###############################################################################
# =========================
# Auto-build Quarto Report
# =========================
###############################################################################
{
  timestamp   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  script_path_q <- .get_script_path()
  script_name_q <- if (!is.na(script_path_q)) basename(script_path_q) else "UNKNOWN_SCRIPT.R"
  script_dir_q  <- if (!is.na(script_path_q)) dirname(script_path_q)  else getwd()
  
  input_file_path  <- if (exists("input_file")) normalizePath(input_file, winslash = "/") else NA_character_
  output_dir_path  <- if (exists("output_dir")) normalizePath(output_dir, winslash = "/") else getwd()
  
  combined_out_csv <- file.path(output_dir_path, "Correlation_by_Day_with_FDR.csv")
  global_out_csv   <- file.path(output_dir_path, "Global_Correlations_with_FDR.csv")
  top_plot_png     <- file.path(output_dir_path, "Top_Correlations_Faceted.png")
  
  cor_by_day_dir_q  <- file.path(output_dir_path, "correlations_by_day")
  expr_by_day_dir_q <- file.path(output_dir_path, "expression_levels_by_day")
  
  header_text <- character(0)
  if (!is.na(script_path_q) && file.exists(script_path_q)) {
    sl <- readLines(script_path_q, warn = FALSE)
    if (length(sl)) {
      is_comment <- grepl("^\\s*#", sl)
      if (length(is_comment) && is_comment[1]) {
        end_idx <- which(!is_comment)
        end_idx <- if (length(end_idx)) end_idx[1] - 1 else length(sl)
        block <- sl[seq_len(end_idx)]
        header_text <- gsub("^\\s*#\\s?", "", block)
      }
    }
  }
  
  deps <- c("parallel", "stats", "ggplot2", "pbapply")
  
  yaml <- c(
    "---",
    paste0("title: \"", tools::file_path_sans_ext(script_name_q), " - Correlation Report\""),
    paste0("author: \"", Sys.info()[["user"]], "\""),
    paste0("date: \"", timestamp, "\""),
    "format:",
    "  html:",
    "    toc: true",
    "    toc-depth: 3",
    "    number-sections: false",
    "    theme: cosmo",
    "    df-print: paged",
    "execute:",
    "  echo: true",
    "  warning: false",
    "  message: false",
    "editor: visual",
    "---",
    ""
  )
  
  run_settings <- c(
    "## Run Settings",
    "",
    paste0("- **Scale:** ", analysis_scale_tag),
    if (scale_choice == 1) paste0("- **Pseudocount:** ", pseudocount) else NULL,
    paste0("- **Bootstrap replicates:** ", boot_R),
    paste0("- **Cores used:** ", use_cores),
    paste0("- **Stability |r| threshold:** ", stab_r_thr),
    paste0("- **min_n:** ", min_n),
    paste0("- **FDR method:** ", user_fdr_method),
    paste0("- **FDR cutoff:** ", user_fdr_cutoff),
    ""
  )
  run_settings <- run_settings[!is.null(run_settings)]
  
  header_section <- if (length(header_text)) c(
    "## Script Header",
    "",
    "```",
    header_text,
    "```",
    ""
  ) else character(0)
  
  deps_section <- c(
    "## Dependencies",
    "",
    paste0("- ", deps),
    ""
  )
  
  metadata <- c(
    "## Key Metadata",
    "",
    paste0("- **Script name:** ", script_name_q),
    paste0("- **Script directory:** ", script_dir_q),
    paste0("- **Input file:** ", if (is.na(input_file_path)) "N/A (provided at runtime)" else input_file_path),
    paste0("- **Output directory:** ", output_dir_path),
    paste0("- **Timestamp:** ", timestamp),
    ""
  )
  
  outputs <- c(
    "## Generated Outputs",
    "",
    paste0("- **Combined per-Day correlations:** `", basename(combined_out_csv), "`"),
    paste0("- **Global correlations:** `", basename(global_out_csv), "`"),
    "- **Per-Day correlation tables:** `correlations_by_day/Correlation_by_Day_<Day>.csv`",
    "- **Per-Day expression tables:** `expression_levels_by_day/Group_<Day>.csv`",
    "- **Per-Day correlation matrices:**",
    "  - `correlation_matrices_by_day/Pearson/CorrMatrix_Pearson_Day<Day>.csv`",
    "  - `correlation_matrices_by_day/Spearman/CorrMatrix_Spearman_Day<Day>.csv`",
    "  - `correlation_matrices_by_day/PearsonZ/CorrMatrix_PearsonZ_Day<Day>.csv`",
    "- **Top-changing pairs plot:** `Top_Correlations_Faceted.png`",
    ""
  )
  
  checks_chunk <- c(
    "## Quick Previews",
    "",
    "```{r}",
    "qmd_path <- knitr::current_input()",
    "out_dir <- if (!is.null(qmd_path) && nzchar(qmd_path)) dirname(normalizePath(qmd_path, winslash = \"/\")) else getwd()",
    "cat(\"Using output directory:\", out_dir, \"\\n\\n\")",
    "",
    "combined_csv <- file.path(out_dir, \"Correlation_by_Day_with_FDR.csv\")",
    "global_csv   <- file.path(out_dir, \"Global_Correlations_with_FDR.csv\")",
    "plot_png     <- file.path(out_dir, \"Top_Correlations_Faceted.png\")",
    "",
    "cat(\"Contents of out_dir (first 30):\\n\")",
    "print(utils::head(list.files(out_dir, full.names = TRUE), 30))",
    "cat(\"\\n\\n\")",
    "",
    "if (file.exists(combined_csv)) {",
    "  cat(\"\\n### Combined Results (head)\\n\\n\")",
    "  cor_df <- read.csv(combined_csv, nrows = 2000, check.names = FALSE)",
    "  print(utils::head(cor_df, 10))",
    "} else {",
    "  cat(\"\\n(Combined results CSV not found at render time)\\n\\n\")",
    "}",
    "",
    "if (file.exists(global_csv)) {",
    "  gdf <- read.csv(global_csv, check.names = FALSE)",
    "  cat(\"\\n### Global Summary Statistics\\n\")",
    "  cat(\"Total pairs:\", nrow(gdf), \"\\n\")",
    "  sig_pearson <- sum(gdf$Pearson_FDR_BH < 0.05, na.rm = TRUE)",
    "  sig_spearman <- sum(gdf$Spearman_FDR_BH < 0.05, na.rm = TRUE)",
    "  cat(\"Significant Pearson (BH FDR < 0.05):\", sig_pearson, \"\\n\")",
    "  cat(\"Significant Spearman (BH FDR < 0.05):\", sig_spearman, \"\\n\")",
    "  cat(\"\\nTop 5 |rho - r| differences:\\n\")",
    "  print(utils::head(gdf[order(-gdf$AbsDiff), c('Pair','Pearson','Spearman','AbsDiff')], 5))",
    "  library(ggplot2)",
    "  print(ggplot(gdf, aes(x = Pearson, y = Spearman, color = AbsDiff)) +",
    "        geom_point(alpha = 0.6) +",
    "        theme_minimal() +",
    "        labs(title = 'Global Pearson vs Spearman',",
    "             x = 'Pearson r', y = 'Spearman rho', color = '|rho - r|'))",
    "  print(ggplot(gdf, aes(x = AbsDiff)) +",
    "        geom_histogram(bins = 30, fill = 'skyblue', color = 'white') +",
    "        theme_minimal() +",
    "        labs(title = '|rho - r| distribution', x = '|rho - r|', y = 'Count'))",
    "} else {",
    "  cat(\"\\n(Global results CSV not found at render time)\\n\\n\")",
    "}",
    "",
    "if (file.exists(plot_png)) {",
    "  cat(\"\\n### Top Pairs Plot\\n\\n\")",
    "  knitr::include_graphics(plot_png)",
    "} else {",
    "  cat(\"\\n(Top plot image not found at render time)\\n\\n\")",
    "}",
    "```",
    ""
  )
  
  qmd_lines <- c(
    yaml,
    run_settings,
    header_section,
    deps_section,
    metadata,
    outputs,
    checks_chunk
  )
  
  report_base <- paste0(tools::file_path_sans_ext(script_name_q), "_Report_", timestamp_tag)
  qmd_dir     <- if (dir.exists(output_dir_path)) output_dir_path else script_dir_q
  qmd_path    <- file.path(qmd_dir, paste0(report_base, ".qmd"))
  
  writeLines(qmd_lines, con = qmd_path)
  message("Wrote QMD: ", qmd_path)
  
  if (requireNamespace("quarto", quietly = TRUE)) {
    message("Rendering QMD via quarto R package...")
    tryCatch({
      quarto::quarto_render(qmd_path)
      message("Rendered: ", sub("\\.qmd$", ".html", qmd_path))
    }, error = function(e) {
      message("Quarto render failed: ", conditionMessage(e),
              "\nOpen the QMD in RStudio and click 'Render', or install Quarto.")
    })
  } else {
    message("quarto R package not found; trying Quarto CLI if available...")
    ok <- FALSE
    tryCatch({
      res <- suppressWarnings(system2("quarto", c("render", shQuote(qmd_path)), stdout = TRUE, stderr = TRUE))
      ok <- TRUE
      message(paste(res, collapse = "\n"))
    }, error = function(e) {
      ok <- FALSE
    })
    if (!ok) {
      message("Quarto not available. To render HTML, install Quarto (https://quarto.org) or install the R package 'quarto'.")
    }
  }
}
