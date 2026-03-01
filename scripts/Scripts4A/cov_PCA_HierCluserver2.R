#!/usr/bin/env Rscript
# ======================================================================
# NAME: Omics_Tidy_ByDay_MeanCoV_MAD_MADresid_PCA_Clustering.R
#
# PURPOSE
#   From a TIDY omics table (biological replicates per day), compute:
#     (A) By-day Mean / Variance / SD / CoV / MAD per Source×Feature×Day
#     (B) Mean-adjusted robust dispersion: MAD_resid per Source×Feature×Day
#         - For each Source×Day, fit log(MAD+eps) ~ smooth(log(|Mean|+eps)) across features
#         - Define MAD_resid as residuals (mean-independent by construction within that day)
#     (C) Feature-by-day matrices for:
#         1) Mean trajectories      (Features × Days)
#         2) CoV trajectories       (Features × Days)
#         3) MAD trajectories       (Features × Days)
#         4) MAD_resid trajectories (Features × Days)  [mean-adjusted robust dispersion]
#     (D) PCA + hierarchical clustering (Ward.D2) on all matrices
#     (E) Multi-Source operation:
#         - Runs per Source independently (separate PCA/HC per Source)
#         - Also optional "ALL sources combined" PCA (rows = Source||Metabolite)
#     (F) Writes outputs to a NEW run directory and generates an auto-executing QMD report
#         including diagnostics:
#           correlation(Mean, MAD) vs correlation(Mean, MAD_resid), per day (per source)
#
# INPUT (TIDY; EXACT 5 columns)
#   Source, Day, SampleID, Metabolite, Abundance
#
# USER OPTIONS
#   - Transform:
#       0 = raw input Abundance
#       1 = log2(Abundance + pseudocount)
#   - Pseudocount: numeric (default 1)
#
# REPLICATE FILTER (requested)
#   - Dispersion metrics (Variance/SD/CoV/MAD) require N >= 3 replicates per day.
#
# PCA SCALING (explicit, requested)
#   - Mean PCA:      center=TRUE, scale.=TRUE  (shape-based trajectories)
#   - CoV PCA:       center=TRUE, scale.=FALSE (preserve amplitude)
#   - MAD PCA:       center=TRUE, scale.=FALSE (preserve amplitude)
#   - MAD_resid PCA: center=TRUE, scale.=FALSE (preserve amplitude)
#
# NOTES
#   - MAD is robust to outliers but not inherently mean-independent.
#   - MAD_resid is the mean-adjusted robust dispersion, computed per Source×Day.
#   - Residualization is performed on the analysis scale (raw or log2+pc), with log links inside the fit.
# ======================================================================

# -----------------------------
# Minimal packages
# -----------------------------
need_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
need_pkg("tools")
need_pkg("knitr")
have_quarto <- requireNamespace("quarto", quietly = TRUE)

library(tools)
library(knitr)

# -----------------------------
# User settings (defaults)
# -----------------------------
set_seed <- 123
coverage_threshold <- 0.75
cv_eps <- 1e-6
hc_method <- "ward.D2"
k_clusters_mean <- 5
k_clusters_cov  <- 5
k_clusters_mad  <- 5
k_clusters_madresid <- 5

min_reps_dispersion <- 3

do_all_sources_combined <- TRUE

# Mean–MAD residualization settings
trend_eps <- 1e-8         # protects logs for MAD/Mean when very small
loess_span <- 0.75        # smoothing span for loess
min_features_for_trend <- 10  # require at least this many features to fit a trend per Source×Day

# -----------------------------
# Utilities
# -----------------------------
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
mk <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
norm_abs <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p) || !nzchar(p)) return(NA_character_)
  tryCatch(normalizePath(p, winslash = "/", mustWork = FALSE),
           error = function(e) as.character(p))
}

stop_schema <- function(required_cols, df) {
  missing_cols <- setdiff(required_cols, names(df))
  extra_cols   <- setdiff(names(df), required_cols)
  if (length(missing_cols) > 0) {
    stop("Input must be tidy with EXACT columns:\n  ",
         paste(required_cols, collapse = ", "),
         "\nMissing:\n  ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (length(extra_cols) > 0) {
    stop("Input must contain ONLY these columns:\n  ",
         paste(required_cols, collapse = ", "),
         "\nExtra columns detected:\n  ", paste(extra_cols, collapse = ", "), call. = FALSE)
  }
}

ask_choice <- function(prompt, allowed = c("0","1"), default = "0") {
  repeat {
    ans <- readline(prompt)
    ans <- trimws(ans)
    if (ans == "" && !is.null(default)) ans <- default
    if (ans %in% allowed) return(ans)
    cat("Invalid input. Allowed values: ", paste(allowed, collapse=", "),
        ". Default is ", default, ".\n", sep = "")
  }
}

ask_numeric <- function(prompt, default = 1) {
  repeat {
    ans <- readline(prompt)
    ans <- trimws(ans)
    if (ans == "") return(as.numeric(default))
    val <- suppressWarnings(as.numeric(ans))
    if (is.finite(val)) return(val)
    cat("Invalid numeric input. Please enter a finite number.\n")
  }
}

# -----------------------------
# I/O + Transform selection
# -----------------------------
cat("\nSelect input CSV (TIDY; EXACT 5 columns):\n")
cat("  Source, Day, SampleID, Metabolite, Abundance\n")
infile <- file.choose()

cat("\nChoose abundance transform:\n")
cat("  0 = raw input Abundance\n")
cat("  1 = log2(Abundance + pseudocount)\n")
transform_choice <- ask_choice("Enter 0 or 1 [default 0]: ", allowed = c("0","1"), default = "0")
transform_choice <- as.integer(transform_choice)

pseudocount <- 1
if (transform_choice == 1) {
  pseudocount <- ask_numeric("Enter pseudocount for log2(Abundance + pseudocount) [default 1]: ", default = 1)
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    stop("Pseudocount must be a finite positive number.", call. = FALSE)
  }
}

transform_label <- if (transform_choice == 0) "raw" else paste0("log2_pc_", pseudocount)

cat("\nSelect a base output directory (pick ANY file inside the desired folder):\n")
base_out_dir <- dirname(file.choose())
if (is.na(base_out_dir) || base_out_dir == "") base_out_dir <- dirname(infile)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
script_name <- "Omics_Tidy_ByDay_MeanCoV_MAD_MADresid_PCA_Clustering.R"
run_dir <- file.path(base_out_dir, paste0("MeanDispersion_PCAClust_Run_", timestamp, "_", transform_label))
tables_dir <- file.path(run_dir, "tables")
plots_dir  <- file.path(run_dir, "plots")
report_dir <- file.path(run_dir, "report")
mk(tables_dir); mk(plots_dir); mk(report_dir)

cat("\nRun directory:\n  ", run_dir, "\n", sep = "")
cat("Transform:\n  ", transform_label, "\n", sep = "")

# -----------------------------
# Read + enforce schema
# -----------------------------
required_cols <- c("Source","Day","SampleID","Metabolite","Abundance")
df_raw <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)
stop_schema(required_cols, df_raw)

# Coerce types safely
df_raw$Source     <- trimws(as.character(df_raw$Source))
df_raw$SampleID   <- trimws(as.character(df_raw$SampleID))
df_raw$Metabolite <- trimws(as.character(df_raw$Metabolite))
df_raw$Day        <- suppressWarnings(as.numeric(as.character(df_raw$Day)))
df_raw$Abundance  <- suppressWarnings(as.numeric(as.character(df_raw$Abundance)))

# Basic checks
if (any(!is.finite(df_raw$Day))) stop("Day contains non-numeric/NA after coercion.", call. = FALSE)
if (any(!is.finite(df_raw$Abundance))) stop("Abundance contains non-numeric/NA after coercion.", call. = FALSE)
if (any(df_raw$Source == "")) stop("Source has empty strings after trimws().", call. = FALSE)
if (any(df_raw$Metabolite == "")) stop("Metabolite has empty strings after trimws().", call. = FALSE)

# SampleID policy
df_raw$SampleID[df_raw$SampleID %in% c("NA","N/A","na","n/a","NULL","null",".","-","")] <- NA_character_
bad_id <- which(is.na(df_raw$SampleID))
if (length(bad_id) > 0) {
  warning("SampleID missing in ", length(bad_id), " rows. Filling placeholders deterministically.")
  df_raw$SampleID[bad_id] <- paste0("MISSING_SAMPLEID_ROW_", bad_id)
}

set.seed(set_seed)

cat("\nTidy input accepted.\n")
cat("Rows:", nrow(df_raw),
    " | Unique Sources:", length(unique(df_raw$Source)),
    " | Unique Days:", length(unique(df_raw$Day)),
    " | Unique Features:", length(unique(df_raw$Metabolite)), "\n\n", sep = "")

# -----------------------------
# Apply transform
# -----------------------------
tidy <- df_raw
if (transform_choice == 1) {
  if (any(tidy$Abundance < 0, na.rm = TRUE)) {
    stop("Negative Abundance values found; log2(Abundance + pseudocount) is invalid.", call. = FALSE)
  }
  tidy$Abundance <- log2(tidy$Abundance + pseudocount)
}

# ======================================================================
# (A) By-day mean/variance/SD/CoV/MAD per Source×Metabolite×Day
# ======================================================================
sources <- sort(unique(tidy$Source))

byday_list <- list()
ii <- 1

for (src in sources) {
  dsrc <- tidy[tidy$Source == src, , drop = FALSE]
  days <- sort(unique(dsrc$Day))
  mets <- sort(unique(dsrc$Metabolite))
  
  for (met in mets) {
    dmet <- dsrc[dsrc$Metabolite == met, , drop = FALSE]
    for (dd in days) {
      v <- dmet$Abundance[dmet$Day == dd]
      v <- v[is.finite(v)]
      n <- length(v)
      
      mu <- if (n > 0) mean(v) else NA_real_
      
      if (n >= min_reps_dispersion) {
        va <- var(v)
        sdv <- sd(v)
        cv <- if (is.finite(mu) && abs(mu) > cv_eps) sdv / abs(mu) else NA_real_
        madv <- stats::mad(v, center = stats::median(v), constant = 1.4826, na.rm = TRUE)
      } else {
        va <- NA_real_
        sdv <- NA_real_
        cv <- NA_real_
        madv <- NA_real_
      }
      
      byday_list[[ii]] <- data.frame(
        Source = src,
        Metabolite = met,
        Day = dd,
        N = n,
        Mean = mu,
        Variance = va,
        SD = sdv,
        CoV = cv,
        MAD = madv,
        MAD_resid = NA_real_,   # filled below
        Transform = transform_label,
        stringsAsFactors = FALSE
      )
      ii <- ii + 1
    }
  }
}

byday_df <- do.call(rbind, byday_list)

# ======================================================================
# (B) Compute MAD_resid per Source×Day (mean-adjusted robust dispersion)
#   Fit across metabolites within each Source×Day:
#     y = log(MAD + eps)
#     x = log(|Mean| + eps)
#   Use LOESS smoothing; MAD_resid = y - yhat
# ======================================================================
cat("Computing MAD_resid (mean-adjusted robust dispersion) per Source×Day...\n")

for (src in sources) {
  dsrc_idx <- which(byday_df$Source == src)
  if (length(dsrc_idx) == 0) next
  days <- sort(unique(byday_df$Day[dsrc_idx]))
  
  for (dd in days) {
    idx <- which(byday_df$Source == src & byday_df$Day == dd)
    if (length(idx) == 0) next
    
    sub <- byday_df[idx, , drop = FALSE]
    ok <- is.finite(sub$Mean) & is.finite(sub$MAD) & (sub$MAD >= 0)
    
    # Need enough features to fit
    if (sum(ok) < min_features_for_trend) {
      byday_df$MAD_resid[idx] <- NA_real_
      next
    }
    
    x <- log(abs(sub$Mean[ok]) + trend_eps)
    y <- log(sub$MAD[ok] + trend_eps)
    
    # LOESS fit; protect with tryCatch
    fit <- tryCatch({
      stats::loess(y ~ x, span = loess_span, degree = 2,
                   control = stats::loess.control(surface = "direct"))
    }, error = function(e) NULL)
    
    if (is.null(fit)) {
      byday_df$MAD_resid[idx] <- NA_real_
      next
    }
    
    yhat <- tryCatch(stats::predict(fit, newdata = data.frame(x = x)),
                     error = function(e) rep(NA_real_, length(x)))
    
    resid <- y - yhat
    
    # Fill back
    out <- rep(NA_real_, nrow(sub))
    out[ok] <- resid
    byday_df$MAD_resid[idx] <- out
  }
}

byday_path <- file.path(tables_dir, "ByDay_MeanVarSD_CoV_MAD_MADresid.csv")
write.csv(byday_df, byday_path, row.names = FALSE)

# ======================================================================
# Diagnostics: correlation(Mean, MAD) vs correlation(Mean, MAD_resid) per Source×Day
#   computed across metabolites within each day
# ======================================================================
cat("Computing per-day diagnostics correlations...\n")
diag_list <- list()
di <- 1

for (src in sources) {
  dsrc <- byday_df[byday_df$Source == src, , drop = FALSE]
  days <- sort(unique(dsrc$Day))
  
  for (dd in days) {
    sub <- dsrc[dsrc$Day == dd, , drop = FALSE]
    
    ok_mad <- is.finite(sub$Mean) & is.finite(sub$MAD)
    ok_res <- is.finite(sub$Mean) & is.finite(sub$MAD_resid)
    
    # Pearson + Spearman (helps when distributions are heavy-tailed)
    cor_mean_mad_p <- if (sum(ok_mad) >= 3) suppressWarnings(stats::cor(sub$Mean[ok_mad], sub$MAD[ok_mad], method="pearson")) else NA_real_
    cor_mean_res_p <- if (sum(ok_res) >= 3) suppressWarnings(stats::cor(sub$Mean[ok_res], sub$MAD_resid[ok_res], method="pearson")) else NA_real_
    
    cor_mean_mad_s <- if (sum(ok_mad) >= 3) suppressWarnings(stats::cor(sub$Mean[ok_mad], sub$MAD[ok_mad], method="spearman")) else NA_real_
    cor_mean_res_s <- if (sum(ok_res) >= 3) suppressWarnings(stats::cor(sub$Mean[ok_res], sub$MAD_resid[ok_res], method="spearman")) else NA_real_
    
    diag_list[[di]] <- data.frame(
      Source = src,
      Day = dd,
      Transform = transform_label,
      N_features_MeanMAD = sum(ok_mad),
      N_features_MeanMADresid = sum(ok_res),
      Cor_Pearson_Mean_vs_MAD = cor_mean_mad_p,
      Cor_Pearson_Mean_vs_MADresid = cor_mean_res_p,
      Cor_Spearman_Mean_vs_MAD = cor_mean_mad_s,
      Cor_Spearman_Mean_vs_MADresid = cor_mean_res_s,
      stringsAsFactors = FALSE
    )
    di <- di + 1
  }
}

diag_df <- if (length(diag_list) > 0) do.call(rbind, diag_list) else data.frame()
diag_path <- file.path(tables_dir, "Diagnostics_Cor_Mean_vs_MAD_vs_MADresid_ByDay.csv")
write.csv(diag_df, diag_path, row.names = FALSE)

# ======================================================================
# Helpers to build matrices and run PCA/HC
# ======================================================================
build_feature_by_day_matrix <- function(df_byday,
                                        value_col = c("Mean","CoV","MAD","MAD_resid"),
                                        coverage_threshold = 0.75) {
  value_col <- match.arg(value_col)
  df <- df_byday[, c("Metabolite","Day", value_col)]
  df <- df[is.finite(df$Day), , drop = FALSE]
  df$DayName <- paste0("Day", df$Day)
  
  mets <- sort(unique(df$Metabolite))
  days <- sort(unique(df$DayName))
  
  mat <- matrix(NA_real_, nrow = length(mets), ncol = length(days),
                dimnames = list(mets, days))
  
  for (r in seq_len(nrow(df))) {
    mat[df$Metabolite[r], df$DayName[r]] <- df[[value_col]][r]
  }
  
  keep <- rowMeans(is.finite(mat)) >= coverage_threshold
  mat <- mat[keep, , drop = FALSE]
  
  mat
}

run_pca_and_cluster <- function(mat,
                                center = TRUE,
                                scale. = TRUE,
                                hc_method = "ward.D2",
                                k = 5) {
  keep_cols <- colSums(is.finite(mat)) > 0
  mat2 <- mat[, keep_cols, drop = FALSE]
  
  cc <- complete.cases(mat2)
  mat3 <- mat2[cc, , drop = FALSE]
  
  if (nrow(mat3) < 3 || ncol(mat3) < 2) {
    return(list(
      mat_used = mat3,
      pca = NULL,
      hc = NULL,
      clusters = NULL,
      dropped_rows = sum(!cc),
      dropped_cols_allNA = sum(!keep_cols)
    ))
  }
  
  pca <- prcomp(mat3, center = center, scale. = scale.)
  hc  <- hclust(dist(mat3), method = hc_method)
  clusters <- cutree(hc, k = k)
  
  list(
    mat_used = mat3,
    pca = pca,
    hc = hc,
    clusters = clusters,
    dropped_rows = sum(!cc),
    dropped_cols_allNA = sum(!keep_cols)
  )
}

save_pca_tables <- function(pca_obj, prefix, outdir) {
  scores <- as.data.frame(pca_obj$x)
  scores$Metabolite <- rownames(scores)
  
  loadings <- as.data.frame(pca_obj$rotation)
  loadings$Day <- rownames(loadings)
  
  write.csv(scores,  file.path(outdir, paste0(prefix, "_Scores.csv")), row.names = FALSE)
  write.csv(loadings,file.path(outdir, paste0(prefix, "_Loadings.csv")), row.names = FALSE)
  
  ve <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  ve_df <- data.frame(PC = paste0("PC", seq_along(ve)), VarExplained = ve, stringsAsFactors = FALSE)
  write.csv(ve_df, file.path(outdir, paste0(prefix, "_VarianceExplained.csv")), row.names = FALSE)
  
  invisible(TRUE)
}

plot_pca2 <- function(pca_obj, title, outfile, label_top = 0) {
  ve <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  xlab <- paste0("PC1 (", round(100*ve[1], 1), "%)")
  ylab <- paste0("PC2 (", round(100*ve[2], 1), "%)")
  
  scores <- as.data.frame(pca_obj$x[, 1:2, drop = FALSE])
  scores$Metabolite <- rownames(scores)
  
  png(outfile, width = 1400, height = 1100, res = 150)
  plot(scores[,1], scores[,2], pch = 16,
       xlab = xlab, ylab = ylab, main = title)
  if (is.numeric(label_top) && label_top > 0) {
    r <- sqrt(scores[,1]^2 + scores[,2]^2)
    ord <- order(r, decreasing = TRUE)
    pick <- head(ord, min(label_top, length(ord)))
    text(scores[pick,1], scores[pick,2], labels = scores$Metabolite[pick], pos = 4, cex = 0.7)
  }
  dev.off()
}

plot_dendro <- function(hc_obj, title, outfile) {
  png(outfile, width = 1400, height = 1100, res = 150)
  plot(hc_obj, main = title, xlab = "", sub = "", cex = 0.6)
  dev.off()
}

# ======================================================================
# (C) Per-source: build matrices, run PCA/HC, write outputs
# ======================================================================
results_index <- list()
ri <- 1

for (src in sources) {
  dsrc <- byday_df[byday_df$Source == src, , drop = FALSE]
  if (nrow(dsrc) == 0) next
  
  MeanMat <- build_feature_by_day_matrix(dsrc, value_col = "Mean", coverage_threshold = coverage_threshold)
  CoVMat  <- build_feature_by_day_matrix(dsrc, value_col = "CoV",  coverage_threshold = coverage_threshold)
  MADMat  <- build_feature_by_day_matrix(dsrc, value_col = "MAD",  coverage_threshold = coverage_threshold)
  MRMat   <- build_feature_by_day_matrix(dsrc, value_col = "MAD_resid", coverage_threshold = coverage_threshold)
  
  write.csv(MeanMat, file.path(tables_dir, paste0("MeanMatrix_", safe_name(src), ".csv")))
  write.csv(CoVMat,  file.path(tables_dir, paste0("CoVMatrix_",  safe_name(src), ".csv")))
  write.csv(MADMat,  file.path(tables_dir, paste0("MADMatrix_",  safe_name(src), ".csv")))
  write.csv(MRMat,   file.path(tables_dir, paste0("MADresidMatrix_", safe_name(src), ".csv")))
  
  # PCA + HC (explicit scaling per matrix)
  mean_res <- run_pca_and_cluster(MeanMat, center = TRUE, scale. = TRUE,  hc_method = hc_method, k = k_clusters_mean)
  cov_res  <- run_pca_and_cluster(CoVMat,  center = TRUE, scale. = FALSE, hc_method = hc_method, k = k_clusters_cov)
  mad_res  <- run_pca_and_cluster(MADMat,  center = TRUE, scale. = FALSE, hc_method = hc_method, k = k_clusters_mad)
  mr_res   <- run_pca_and_cluster(MRMat,   center = TRUE, scale. = FALSE, hc_method = hc_method, k = k_clusters_madresid)
  
  # Clusters
  if (!is.null(mean_res$clusters)) {
    write.csv(data.frame(Metabolite = names(mean_res$clusters), MeanCluster = as.integer(mean_res$clusters)),
              file.path(tables_dir, paste0("Clusters_Mean_", safe_name(src), ".csv")), row.names = FALSE)
  } else write.csv(data.frame(), file.path(tables_dir, paste0("Clusters_Mean_", safe_name(src), ".csv")), row.names = FALSE)
  
  if (!is.null(cov_res$clusters)) {
    write.csv(data.frame(Metabolite = names(cov_res$clusters), CoVCluster = as.integer(cov_res$clusters)),
              file.path(tables_dir, paste0("Clusters_CoV_", safe_name(src), ".csv")), row.names = FALSE)
  } else write.csv(data.frame(), file.path(tables_dir, paste0("Clusters_CoV_", safe_name(src), ".csv")), row.names = FALSE)
  
  if (!is.null(mad_res$clusters)) {
    write.csv(data.frame(Metabolite = names(mad_res$clusters), MADCluster = as.integer(mad_res$clusters)),
              file.path(tables_dir, paste0("Clusters_MAD_", safe_name(src), ".csv")), row.names = FALSE)
  } else write.csv(data.frame(), file.path(tables_dir, paste0("Clusters_MAD_", safe_name(src), ".csv")), row.names = FALSE)
  
  if (!is.null(mr_res$clusters)) {
    write.csv(data.frame(Metabolite = names(mr_res$clusters), MADresidCluster = as.integer(mr_res$clusters)),
              file.path(tables_dir, paste0("Clusters_MADresid_", safe_name(src), ".csv")), row.names = FALSE)
  } else write.csv(data.frame(), file.path(tables_dir, paste0("Clusters_MADresid_", safe_name(src), ".csv")), row.names = FALSE)
  
  # PCA tables + plots
  if (!is.null(mean_res$pca)) {
    save_pca_tables(mean_res$pca, paste0("PCA_Mean_", safe_name(src)), tables_dir)
    plot_pca2(mean_res$pca, paste0("PCA (Means; centered+scaled) — ", src),
              file.path(plots_dir, paste0("PCA_Mean_", safe_name(src), ".png")), label_top = 10)
    if (!is.null(mean_res$hc)) plot_dendro(mean_res$hc, paste0("HC (Means) — ", src),
                                           file.path(plots_dir, paste0("Dendrogram_Mean_", safe_name(src), ".png")))
  }
  
  if (!is.null(cov_res$pca)) {
    save_pca_tables(cov_res$pca, paste0("PCA_CoV_", safe_name(src)), tables_dir)
    plot_pca2(cov_res$pca, paste0("PCA (CoV; centered, NOT scaled) — ", src),
              file.path(plots_dir, paste0("PCA_CoV_", safe_name(src), ".png")), label_top = 10)
    if (!is.null(cov_res$hc)) plot_dendro(cov_res$hc, paste0("HC (CoV) — ", src),
                                          file.path(plots_dir, paste0("Dendrogram_CoV_", safe_name(src), ".png")))
  }
  
  if (!is.null(mad_res$pca)) {
    save_pca_tables(mad_res$pca, paste0("PCA_MAD_", safe_name(src)), tables_dir)
    plot_pca2(mad_res$pca, paste0("PCA (MAD; centered, NOT scaled) — ", src),
              file.path(plots_dir, paste0("PCA_MAD_", safe_name(src), ".png")), label_top = 10)
    if (!is.null(mad_res$hc)) plot_dendro(mad_res$hc, paste0("HC (MAD) — ", src),
                                          file.path(plots_dir, paste0("Dendrogram_MAD_", safe_name(src), ".png")))
  }
  
  if (!is.null(mr_res$pca)) {
    save_pca_tables(mr_res$pca, paste0("PCA_MADresid_", safe_name(src)), tables_dir)
    plot_pca2(mr_res$pca, paste0("PCA (MAD_resid; centered, NOT scaled) — ", src),
              file.path(plots_dir, paste0("PCA_MADresid_", safe_name(src), ".png")), label_top = 10)
    if (!is.null(mr_res$hc)) plot_dendro(mr_res$hc, paste0("HC (MAD_resid) — ", src),
                                         file.path(plots_dir, paste0("Dendrogram_MADresid_", safe_name(src), ".png")))
  }
  
  results_index[[ri]] <- data.frame(
    Source = src,
    Transform = transform_label,
    MinRepsDispersion = min_reps_dispersion,
    N_features_mean_matrix = nrow(MeanMat),
    N_features_cov_matrix  = nrow(CoVMat),
    N_features_mad_matrix  = nrow(MADMat),
    N_features_madresid_matrix = nrow(MRMat),
    DroppedRows_Mean_PCA_NA = mean_res$dropped_rows,
    DroppedRows_CoV_PCA_NA  = cov_res$dropped_rows,
    DroppedRows_MAD_PCA_NA  = mad_res$dropped_rows,
    DroppedRows_MADresid_PCA_NA = mr_res$dropped_rows,
    DroppedColsAllNA_Mean = mean_res$dropped_cols_allNA,
    DroppedColsAllNA_CoV  = cov_res$dropped_cols_allNA,
    DroppedColsAllNA_MAD  = mad_res$dropped_cols_allNA,
    DroppedColsAllNA_MADresid = mr_res$dropped_cols_allNA,
    stringsAsFactors = FALSE
  )
  ri <- ri + 1
}

index_df <- if (length(results_index) > 0) do.call(rbind, results_index) else data.frame()
write.csv(index_df, file.path(tables_dir, "Run_Index_BySource.csv"), row.names = FALSE)

# ======================================================================
# (D) Optional combined PCA across sources (rows = Source||Metabolite)
# ======================================================================
if (isTRUE(do_all_sources_combined) && length(sources) > 1) {
  df_all <- byday_df
  df_all$RowKey <- paste(df_all$Source, df_all$Metabolite, sep = "||")
  df_all$DayName <- paste0("Day", df_all$Day)
  
  rowkeys <- sort(unique(df_all$RowKey))
  daynames <- sort(unique(df_all$DayName))
  
  MeanAll <- matrix(NA_real_, nrow = length(rowkeys), ncol = length(daynames),
                    dimnames = list(rowkeys, daynames))
  CoVAll  <- MeanAll
  MADAll  <- MeanAll
  MRAll   <- MeanAll
  
  for (r in seq_len(nrow(df_all))) {
    rk <- df_all$RowKey[r]
    dk <- df_all$DayName[r]
    MeanAll[rk, dk] <- df_all$Mean[r]
    CoVAll[rk, dk]  <- df_all$CoV[r]
    MADAll[rk, dk]  <- df_all$MAD[r]
    MRAll[rk, dk]   <- df_all$MAD_resid[r]
  }
  
  keep_m <- rowMeans(is.finite(MeanAll)) >= coverage_threshold
  keep_c <- rowMeans(is.finite(CoVAll))  >= coverage_threshold
  keep_a <- rowMeans(is.finite(MADAll))  >= coverage_threshold
  keep_r <- rowMeans(is.finite(MRAll))   >= coverage_threshold
  
  MeanAll <- MeanAll[keep_m, , drop = FALSE]
  CoVAll  <- CoVAll[keep_c, , drop = FALSE]
  MADAll  <- MADAll[keep_a, , drop = FALSE]
  MRAll   <- MRAll[keep_r, , drop = FALSE]
  
  write.csv(MeanAll, file.path(tables_dir, "MeanMatrix_ALL_Sources.csv"))
  write.csv(CoVAll,  file.path(tables_dir, "CoVMatrix_ALL_Sources.csv"))
  write.csv(MADAll,  file.path(tables_dir, "MADMatrix_ALL_Sources.csv"))
  write.csv(MRAll,   file.path(tables_dir, "MADresidMatrix_ALL_Sources.csv"))
  
  mean_all <- run_pca_and_cluster(MeanAll, center = TRUE, scale. = TRUE,  hc_method = hc_method, k = 6)
  cov_all  <- run_pca_and_cluster(CoVAll,  center = TRUE, scale. = FALSE, hc_method = hc_method, k = 6)
  mad_all  <- run_pca_and_cluster(MADAll,  center = TRUE, scale. = FALSE, hc_method = hc_method, k = 6)
  mr_all   <- run_pca_and_cluster(MRAll,   center = TRUE, scale. = FALSE, hc_method = hc_method, k = 6)
  
  if (!is.null(mean_all$pca)) {
    save_pca_tables(mean_all$pca, "PCA_Mean_ALL_Sources", tables_dir)
    plot_pca2(mean_all$pca, "PCA (Means; centered+scaled) — ALL Sources (rows=Source||Metabolite)",
              file.path(plots_dir, "PCA_Mean_ALL_Sources.png"), label_top = 0)
  }
  if (!is.null(cov_all$pca)) {
    save_pca_tables(cov_all$pca, "PCA_CoV_ALL_Sources", tables_dir)
    plot_pca2(cov_all$pca, "PCA (CoV; centered, NOT scaled) — ALL Sources (rows=Source||Metabolite)",
              file.path(plots_dir, "PCA_CoV_ALL_Sources.png"), label_top = 0)
  }
  if (!is.null(mad_all$pca)) {
    save_pca_tables(mad_all$pca, "PCA_MAD_ALL_Sources", tables_dir)
    plot_pca2(mad_all$pca, "PCA (MAD; centered, NOT scaled) — ALL Sources (rows=Source||Metabolite)",
              file.path(plots_dir, "PCA_MAD_ALL_Sources.png"), label_top = 0)
  }
  if (!is.null(mr_all$pca)) {
    save_pca_tables(mr_all$pca, "PCA_MADresid_ALL_Sources", tables_dir)
    plot_pca2(mr_all$pca, "PCA (MAD_resid; centered, NOT scaled) — ALL Sources (rows=Source||Metabolite)",
              file.path(plots_dir, "PCA_MADresid_ALL_Sources.png"), label_top = 0)
  }
}

# ======================================================================
# (E) Write QMD + render (optional)
# ======================================================================
qmd_path <- file.path(report_dir, "run_report.qmd")
html_path <- file.path(report_dir, "run_report.html")

q <- function(x) paste0('"', gsub('"', '\\"', x), '"')

qmd_lines <- c(
  "---",
  "title: \"Mean vs Dispersion (CoV, MAD, MAD_resid) PCA/Clustering — Run Report\"",
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 3",
  "    number-sections: false",
  "    code-fold: true",
  "    code-tools: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "## Run metadata",
  "",
  paste0("- **Script name:** ", script_name),
  paste0("- **Timestamp:** ", timestamp),
  paste0("- **Input file:** `", norm_abs(infile), "`"),
  paste0("- **Run directory:** `", norm_abs(run_dir), "`"),
  paste0("- **Transform:** ", transform_label),
  paste0("- **Min replicates for dispersion metrics (Variance/SD/CoV/MAD):** ", min_reps_dispersion),
  paste0("- **Coverage threshold (row keep):** ", coverage_threshold),
  paste0("- **MAD_resid method:** per Source×Day LOESS on log(MAD+eps) ~ log(|Mean|+eps) across features; residuals retained"),
  "",
  "## Output tables",
  "",
  "```{r}",
  "library(knitr)",
  "tbl_dir <- ", q(norm_abs(tables_dir)),
  "files <- list.files(tbl_dir, full.names = TRUE)",
  "df <- data.frame(File = basename(files), Path = files, stringsAsFactors = FALSE)",
  "kable(df)",
  "```",
  "",
  "## Per-source summary",
  "",
  "```{r}",
  "idx <- read.csv(file.path(", q(norm_abs(tables_dir)), ", 'Run_Index_BySource.csv'), stringsAsFactors=FALSE)",
  "kable(idx)",
  "```",
  "",
  "## Diagnostics: mean–dispersion coupling by day",
  "",
  "This table compares, per day, the across-metabolite correlation between **Mean vs raw MAD** and **Mean vs MAD_resid**.",
  "",
  "```{r}",
  "diag <- read.csv(file.path(", q(norm_abs(tables_dir)), ", 'Diagnostics_Cor_Mean_vs_MAD_vs_MADresid_ByDay.csv'), stringsAsFactors=FALSE)",
  "kable(diag)",
  "```",
  "",
  "```{r}",
  "# Simple diagnostic plot: per Source, scatter of Pearson correlations by day",
  "diag <- diag[order(diag$Source, diag$Day), ]",
  "sources <- unique(diag$Source)",
  "png(file.path(", q(norm_abs(plots_dir)), ", 'Diagnostics_CorrelationShift_Pearson.png'), width=1500, height=900, res=150)",
  "par(mfrow=c(1,1))",
  "x <- diag$Cor_Pearson_Mean_vs_MAD",
  "y <- diag$Cor_Pearson_Mean_vs_MADresid",
  "plot(x, y, pch=16, xlab='cor(Mean, MAD) Pearson', ylab='cor(Mean, MAD_resid) Pearson',",
  "     main='Mean–Dispersion Coupling Reduction by Residualization (per Source×Day)')",
  "abline(h=0, v=0, lty=2)",
  "abline(a=0, b=1, lty=3)  # y=x reference",
  "text(x, y, labels=paste0(diag$Source, ':D', diag$Day), pos=4, cex=0.6)",
  "dev.off()",
  "cat('Wrote diagnostic plot: Diagnostics_CorrelationShift_Pearson.png')",
  "```",
  "",
  "## Key plots",
  "",
  "```{r}",
  "plot_dir <- ", q(norm_abs(plots_dir)),
  "pngs <- list.files(plot_dir, pattern='\\\\.png$', full.names=TRUE)",
  "if (length(pngs) == 0) {",
  "  cat('No plots found.')",
  "} else {",
  "  for (p in pngs) {",
  "    cat('\\n### ', basename(p), '\\n\\n', sep='')",
  "    cat('![](', p, '){width=95%}\\n\\n', sep='')",
  "  }",
  "}",
  "```",
  "",
  "## Interpretation guide",
  "",
  "- **Mean PCA** (centered+scaled) groups features by *trajectory shape* of abundance across days.",
  "- **CoV PCA** (centered, not scaled) groups features by *trajectory of relative dispersion* across days.",
  "- **MAD PCA** (centered, not scaled) groups features by *trajectory of robust absolute dispersion* across days.",
  "- **MAD_resid PCA** (centered, not scaled) groups features by *trajectory of robust dispersion after removing mean–dispersion trend* per day."
)

writeLines(qmd_lines, qmd_path)

if (have_quarto) {
  old <- getwd()
  setwd(report_dir)
  on.exit(setwd(old), add = TRUE)
  tryCatch({
    quarto::quarto_render(input = qmd_path, output_file = basename(html_path))
  }, error = function(e) {
    warning("Quarto render failed: ", e$message)
  })
} else {
  warning("Quarto not installed/available. QMD written but HTML not rendered.")
}

# ======================================================================
# Done
# ======================================================================
cat("\nDone.\n")
cat("Run directory:\n  ", run_dir, "\n", sep = "")
cat("Transform:\n  ", transform_label, "\n", sep = "")
cat("By-day table:\n  ", byday_path, "\n", sep = "")
cat("Diagnostics table:\n  ", diag_path, "\n", sep = "")
cat("Tables directory:\n  ", tables_dir, "\n", sep = "")
cat("Plots directory:\n  ", plots_dir, "\n", sep = "")
cat("Report QMD:\n  ", qmd_path, "\n", sep = "")
if (file.exists(html_path)) cat("Report HTML:\n  ", html_path, "\n", sep = "")
cat("\n")
