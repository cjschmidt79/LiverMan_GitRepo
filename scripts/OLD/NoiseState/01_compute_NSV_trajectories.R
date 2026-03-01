#!/usr/bin/env Rscript
# ============================================================
# NSV Starter Script (robust to unknown batch; no conditions)
# Input: HeidiLiiverCompleteTranscriptome.csv
# Output: NSV tables + plots + robustness summaries
#
# UPDATES (this revision)
#   1) log2(FPKM + 1) handling (toggleable)
#   2) interval-aware structural change:
#        - delta_day and S_per_day = S_adj / delta_day
#        - QC warning for unequal time gaps (e.g., missing Day 6)
#   3) future-based parallel subsampling (future.apply)
#   4) short QC summary written to QC_SUMMARY.txt
#
# Core idea: compute NSV trajectories under multiple "nuisance" models:
#   (1) RAW
#   (2) Residualize Day means (RESID_DAY)
#   (3) Residualize Day means + remove top k latent PCs (RESID_DAY_PCk)
#
# This does NOT "prove" batch doesn't exist; it shows whether
# NSV conclusions are robust to plausible hidden technical structure.
# ============================================================

suppressPackageStartupMessages({
  library(stats)
  library(utils)
})
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

# -------------------------- User settings --------------------------
infile <- "WIDE_LiverMetabolome.csv"   # change if needed
outdir <- paste0("NSV_START_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

time_col <- "Day"
id_col   <- "SampleID"

# FPKM handling
apply_log2_fpkm1 <- TRUE  # set FALSE if already log-scale

# Robustness knobs
set.seed(1)
n_boot <- 250              # (reserved) bootstrap replicates
#latent_k_grid <- c(0, 2, 5)  # remove 0/2/5 latent PCs from residuals
use_shrinkage_if_available <- TRUE

# How reviewers will read this (important)
# 
# Including k = 2 looks careful and reasonable.
# 
# Including k = 5 looks defensive and conservative, but can raise the question:
#   
#   “Are you removing biology?”
# 
# So:
#   
#   Exploratory / main analysis: c(0, 2)

#

# Parallelization
use_parallel <- TRUE


# Feature filters
min_sd <- 1e-8             # drop constant genes
top_var_features <- 2000   # cap features for correlation-based steps for speed

n_subsample <- 250
latent_k_grid <- c(0, 2)
parallel_workers <- 8
top_var_features <- 2000

# -------------------------- Helpers --------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

mad_robust <- function(x) {
  stats::mad(x, center = stats::median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
}

safe_scale <- function(X) {
  # scale columns; drop zero-sd columns
  s <- apply(X, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(s) & (s > min_sd)
  X <- X[, keep, drop = FALSE]
  scale(X, center = TRUE, scale = TRUE)
}

# Shrinkage correlation if available; otherwise standard cor
cor_robust <- function(Z) {
  if (use_shrinkage_if_available && requireNamespace("corpcor", quietly = TRUE)) {
    corpcor::cor.shrink(Z, verbose = FALSE)
  } else {
    stats::cor(Z, use = "pairwise.complete.obs", method = "pearson")
  }
}

effective_dimensionality <- function(R) {
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[is.finite(ev) & (ev > 0)]
  if (length(ev) < 2) return(NA_real_)
  (sum(ev)^2) / sum(ev^2)
}

frobenius_cosine_similarity <- function(A, B) {
  num <- sum(A * B, na.rm = TRUE)
  den <- sqrt(sum(A^2, na.rm = TRUE) * sum(B^2, na.rm = TRUE))
  if (!is.finite(den) || den == 0) return(NA_real_)
  num / den
}

# Compute D(t): median per-gene MAD across replicates
compute_D <- function(Xt) {
  per_gene <- apply(Xt, 2, mad_robust)
  stats::median(per_gene, na.rm = TRUE)
}

# Compute O(t): fraction of genes with any replicate beyond robust z threshold
compute_O <- function(Xt, z_thresh = 3) {
  Z <- apply(Xt, 2, function(v) {
    m <- stats::median(v, na.rm = TRUE)
    s <- mad_robust(v)
    if (!is.finite(s) || s <= min_sd) return(rep(0, length(v)))
    (v - m) / s
  })
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
  extreme <- apply(abs(Z), 2, function(v) any(v > z_thresh, na.rm = TRUE))
  mean(extreme, na.rm = TRUE)
}

# Residualize day means gene-by-gene: X ~ Day
residualize_day <- function(X, day) {
  day_f <- as.factor(day)
  M <- model.matrix(~ day_f)
  Xt <- as.matrix(X)
  beta <- solve(t(M) %*% M, t(M) %*% Xt)
  resid <- Xt - M %*% beta
  resid
}

# Remove top k PCs from a matrix (samples x genes)
remove_top_pcs <- function(X, k) {
  if (k <= 0) return(X)
  Z <- scale(X, center = TRUE, scale = FALSE) # center only
  sv <- svd(Z, nu = k, nv = 0)
  scores <- sv$u %*% diag(sv$d[1:k], nrow = k)
  B <- solve(t(scores) %*% scores, t(scores) %*% Z)
  nuisance <- scores %*% B
  Z - nuisance
}

# Choose top-variable genes (to make correlation feasible and stable)
select_top_var <- function(X, top_n) {
  if (ncol(X) <= top_n) return(X)
  v <- apply(X, 2, stats::var, na.rm = TRUE)
  keep <- order(v, decreasing = TRUE)[seq_len(top_n)]
  X[, keep, drop = FALSE]
}

# Compute NSV per day for a given matrix X (samples x genes) and day labels
compute_nsv_table <- function(X, day) {
  days <- sort(unique(day))
  nsv <- data.frame(
    Day = days,
    n = as.integer(NA),
    D = NA_real_,
    ED = NA_real_,
    O = NA_real_
  )
  
  R_by_day <- vector("list", length(days))
  names(R_by_day) <- as.character(days)
  
  for (i in seq_along(days)) {
    d <- days[i]
    idx <- which(day == d)
    Xt <- X[idx, , drop = FALSE]
    nsv$n[i] <- nrow(Xt)
    nsv$D[i]  <- compute_D(Xt)
    nsv$O[i]  <- compute_O(Xt)
    
    Xt2 <- select_top_var(Xt, top_var_features)
    Z <- safe_scale(Xt2)
    R <- cor_robust(Z)
    nsv$ED[i] <- effective_dimensionality(R)
    R_by_day[[i]] <- R
  }
  
  # Interval-aware structural change
  nsv$delta_day <- c(NA_real_, diff(nsv$Day))
  
  # Structural change S_adj: 1 - similarity between consecutive day correlation matrices
  nsv$S_adj <- NA_real_
  for (i in 2:nrow(nsv)) {
    Rprev <- R_by_day[[i - 1]]
    Rcur  <- R_by_day[[i]]
    
    if (!is.null(rownames(Rprev)) && !is.null(rownames(Rcur))) {
      common <- intersect(rownames(Rprev), rownames(Rcur))
      Rprev <- Rprev[common, common, drop = FALSE]
      Rcur  <- Rcur[common, common, drop = FALSE]
    }
    sim <- frobenius_cosine_similarity(Rprev, Rcur)
    nsv$S_adj[i] <- if (is.finite(sim)) 1 - sim else NA_real_
  }
  
  # Reconfiguration rate (per day) to handle missing day(s)
  nsv$S_per_day <- nsv$S_adj / nsv$delta_day
  
  # Curvature in (D, ED, O) space
  nsv$curvature <- NA_real_
  if (nrow(nsv) >= 3) {
    for (i in 2:(nrow(nsv) - 1)) {
      v1 <- as.numeric(nsv[i, c("D", "ED", "O")]) - as.numeric(nsv[i - 1, c("D", "ED", "O")])
      v2 <- as.numeric(nsv[i + 1, c("D", "ED", "O")]) - as.numeric(nsv[i, c("D", "ED", "O")])
      n1 <- sqrt(sum(v1^2, na.rm = TRUE))
      n2 <- sqrt(sum(v2^2, na.rm = TRUE))
      if (n1 > 0 && n2 > 0) {
        cosang <- sum(v1 * v2, na.rm = TRUE) / (n1 * n2)
        cosang <- max(-1, min(1, cosang))
        nsv$curvature[i] <- 1 - cosang
      }
    }
  }
  
  nsv
}

# Equalize sample sizes per day by subsampling to n_min
subsample_equal_n <- function(X, day, n_min) {
  days <- sort(unique(day))
  idx_keep <- integer(0)
  for (d in days) {
    idx <- which(day == d)
    if (length(idx) < n_min) stop("A day has fewer than n_min samples; cannot subsample equally.")
    idx_keep <- c(idx_keep, sample(idx, n_min, replace = FALSE))
  }
  X2 <- X[idx_keep, , drop = FALSE]
  day2 <- day[idx_keep]
  list(X = X2, day = day2)
}

# -------------------------- Load data --------------------------
cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading:", infile, "\n")
df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(time_col %in% colnames(df), id_col %in% colnames(df))

day <- df[[time_col]]
sample_id <- df[[id_col]]

# Expression matrix
X <- as.matrix(df[, !(colnames(df) %in% c(time_col, id_col)), drop = FALSE])
storage.mode(X) <- "double"

# log2(FPKM+1) handling
if (apply_log2_fpkm1) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] Applying log2(FPKM + 1) transform\n")
  X <- log2(X + 1)
}

# -------------------------- QC summary --------------------------
qc_lines <- character(0)
qc_lines <- c(qc_lines, paste("QC SUMMARY:", Sys.time()))
qc_lines <- c(qc_lines, paste("Input:", infile))
qc_lines <- c(qc_lines, paste("Samples:", nrow(X), "Features:", ncol(X)))

tab_day <- table(day)
qc_lines <- c(qc_lines, "Per-day n:")
qc_lines <- c(qc_lines, paste(capture.output(print(tab_day)), collapse = "\n"))

days_sorted <- sort(unique(day))
deltas <- diff(days_sorted)
qc_lines <- c(qc_lines, paste("Observed days:", paste(days_sorted, collapse = ", ")))
qc_lines <- c(qc_lines, paste("Inter-day gaps:", paste(deltas, collapse = ", ")))

if (length(unique(deltas)) > 1) {
  qc_lines <- c(qc_lines, "WARNING: Unequal time gaps detected. Interpret S_adj as total change over the interval.")
  qc_lines <- c(qc_lines, "         Use S_per_day for approximate rate normalization.")
}
if (any(deltas > min(deltas))) {
  # simple inferred missing-day hint: any gap > min gap
  qc_lines <- c(qc_lines, paste("NOTE: Largest gap =", max(deltas), "days; possible missing intermediate timepoint(s)."))
}

writeLines(qc_lines, file.path(outdir, "QC_SUMMARY.txt"))
cat(paste(qc_lines, collapse = "\n"), "\n\n")

# -------------------------- Basic QC (console) --------------------------
cat("Per-day n:\n")
print(tab_day)

# Drop constant genes globally
sd_all <- apply(X, 2, stats::sd, na.rm = TRUE)
keep <- is.finite(sd_all) & (sd_all > min_sd)
X <- X[, keep, drop = FALSE]
cat("Features after removing constant:", ncol(X), "\n")

n_min <- min(as.integer(tab_day))
cat("n_min for equal-n subsampling:", n_min, "\n")

# -------------------------- Parallel setup (future.apply) --------------------------
if (use_parallel) {
  if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
  if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
  suppressPackageStartupMessages({
    library(future)
    library(future.apply)
  })
  
  # safest cross-platform plan
  future::plan(future::multisession, workers = parallel_workers)
  cat("[", format(Sys.time(), "%H:%M:%S"), "] Parallel subsampling enabled with",
      parallel_workers, "workers\n")
} else {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] Parallel subsampling disabled\n")
}

# -------------------------- Compute NSV: modes --------------------------
modes <- list()

# RAW
modes[["RAW"]] <- list(X = X)

# RESID_DAY
X_resid <- residualize_day(X, day)
modes[["RESID_DAY"]] <- list(X = X_resid)

# RESID_DAY_PCk
for (k in latent_k_grid[latent_k_grid > 0]) {
  Xk <- remove_top_pcs(X_resid, k = k)
  modes[[paste0("RESID_DAY_PC", k)]] <- list(X = Xk)
}

# Main result tables + robustness by subsampling
all_results <- list()
all_subsample <- list()

for (nm in names(modes)) {
  cat("\n[", format(Sys.time(), "%H:%M:%S"), "] Mode:", nm, "\n")
  Xm <- modes[[nm]]$X
  
  # Full-data NSV
  nsv <- compute_nsv_table(Xm, day)
  all_results[[nm]] <- nsv
  write.csv(nsv, file.path(outdir, paste0("NSV_", nm, "_full.csv")), row.names = FALSE)
  
  # Subsampling stability (equal n per day) -- parallel if enabled
  if (use_parallel) {
    sub_tbls <- future.apply::future_lapply(seq_len(n_subsample), function(b) {
      sub <- subsample_equal_n(Xm, day, n_min)
      compute_nsv_table(sub$X, sub$day)
    })
  } else {
    sub_tbls <- vector("list", n_subsample)
    for (b in seq_len(n_subsample)) {
      sub <- subsample_equal_n(Xm, day, n_min)
      sub_tbls[[b]] <- compute_nsv_table(sub$X, sub$day)
    }
  }
  
  # Summarize subsampling distribution per Day
  days <- sort(unique(day))
  summ <- data.frame(Day = days)
  
  # include S_per_day in summary too
  for (col in c("D", "ED", "O", "S_adj", "S_per_day", "curvature")) {
    mat <- sapply(sub_tbls, function(tb) tb[[col]])
    summ[[paste0(col, "_med")]] <- apply(mat, 1, stats::median, na.rm = TRUE)
    summ[[paste0(col, "_lo")]]  <- apply(mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
    summ[[paste0(col, "_hi")]]  <- apply(mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  }
  
  all_subsample[[nm]] <- summ
  write.csv(summ, file.path(outdir, paste0("NSV_", nm, "_subsample_CI.csv")), row.names = FALSE)
}

# -------------------------- Plotting (base R) --------------------------
png(file.path(outdir, "NSV_panels.png"), width = 1600, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot_panel <- function(mode_name, y_col, ylab) {
  s <- all_subsample[[mode_name]]
  x <- s$Day
  y <- s[[paste0(y_col, "_med")]]
  lo <- s[[paste0(y_col, "_lo")]]
  hi <- s[[paste0(y_col, "_hi")]]
  
  plot(x, y, type = "b", pch = 16, xlab = "Day", ylab = ylab,
       main = paste(mode_name, "-", y_col))
  arrows(x0 = x, y0 = lo, x1 = x, y1 = hi, angle = 90, code = 3, length = 0.03)
}

mode_to_plot <- "RAW"
plot_panel(mode_to_plot, "D", "Dispersion D(t) (median gene MAD) [log2(FPKM+1) if enabled]")
plot_panel(mode_to_plot, "ED", "Effective dimensionality ED(t)")
plot_panel(mode_to_plot, "S_adj", "Structural change S_adj (1 - similarity)")
plot_panel(mode_to_plot, "S_per_day", "Structural change rate S_per_day (per day)")
dev.off()

# Phase portrait
png(file.path(outdir, "NSV_phase_D_vs_ED.png"), width = 1000, height = 800)
nsv_raw <- all_results[["RAW"]]
plot(nsv_raw$D, nsv_raw$ED, type = "b", pch = 16,
     xlab = "D(t) dispersion", ylab = "ED(t) effective dimensionality",
     main = "Phase portrait: D vs ED (RAW)")
text(nsv_raw$D, nsv_raw$ED, labels = nsv_raw$Day, pos = 3, cex = 0.9)
dev.off()

# Save a short README
readme <- c(
  "NSV Starter Run",
  paste("Run time:", Sys.time()),
  paste("Input:", infile),
  paste("Log transform:", ifelse(apply_log2_fpkm1, "log2(FPKM+1) applied", "no log transform")),
  "",
  "QC:",
  "- QC_SUMMARY.txt includes time gaps; unequal gaps trigger warnings",
  "",
  "Outputs:",
  "- NSV_<MODE>_full.csv: per-day NSV values on full data",
  "- NSV_<MODE>_subsample_CI.csv: subsampling-based median and 95% CI per day",
  "- NSV_panels.png: D, ED, S_adj, S_per_day panels for RAW",
  "- NSV_phase_D_vs_ED.png: phase portrait for RAW",
  "",
  "Modes:",
  "RAW: no correction",
  "RESID_DAY: regress out day means (removes mean trajectory)",
  "RESID_DAY_PCk: remove top k latent PCs after day-mean residualization (proxy for unknown batch)",
  "",
  "Interval-aware structure:",
  "- S_adj is total reconfiguration over the observed interval",
  "- S_per_day = S_adj / delta_day (approximate rate normalization)"
)
writeLines(readme, file.path(outdir, "README.txt"))

cat("\nDone. Results written to:", outdir, "\n")
