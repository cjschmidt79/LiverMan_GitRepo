# ===============================================================
# 📜 TIME-SERIES GENE EXPRESSION: limma DE  ➜  DTW CLUSTERING
# Biologist-first, SOURCE-to-run, base-R (no dplyr/tidyverse)
# ===============================================================
#
# What this does (in plain terms):
# 1) Reads a WIDE expression CSV with a Day column (rows = samples/replicates).
# 2) Runs limma to find genes that change over time (using replicates if present).
# 3) Builds a per-Day mean time-course for significant genes (for DTW clustering).
# 4) Computes DTW distances, hierarchical clustering, and partitional clustering (k=2..10 or forced k).
# 5) Saves: DE tables, significant genes, per-day mean matrix, DTW distance, dendrogram, silhouette scores, clusters, heatmap.
#
# Expected input format (WIDE):
# - Column "Day" (numeric or integer-like)
# - All other columns = genes (numeric expression)
# - Each row = one biological replicate/sample at that Day (replicates are GOOD)
#
# If your file already contains one row per Day (no replicates),
# limma will not have replication to estimate variance well; the script will still run,
# but DE results will be less trustworthy.
## ===============================================================
# 🔧 USER-ADJUSTABLE PARAMETERS (What you can safely change)
# ===============================================================
#
# Parameter: variance_threshold
# Default : 2
# Meaning : Minimum variance (across biological samples) required
#           for a gene to be considered dynamic.
# Biology : Removes genes that are essentially flat / invariant
#           across development.
# Change if:
#   - You want to be more permissive → LOWER (e.g. 0.5–1)
#   - You want only strongly dynamic genes → RAISE (e.g. 3–5)
#
# ---------------------------------------------------------------
# Parameter: fdr_threshold
# Default : 0.05
# Meaning : Adjusted p-value cutoff for limma differential
#           expression across days.
# Biology : Controls how confidently a gene must change over time
#           to be included.
# Change if:
#   - Exploratory analysis → 0.1
#   - Very stringent / conservative → 0.01
#
# ---------------------------------------------------------------
# Parameter: min_max_expr
# Default : 10   (recommended for this dataset)
# Meaning : Minimum absolute expression a gene must reach at
#           ANY time point to be kept.
# Biology : Removes near-background transcripts that can have
#           inflated relative variance but little biological impact.
# Change if:
#   - Raw counts / high scale data → RAISE (e.g. 50–100)
#   - Normalized / log-like data → KEEP LOW (5–20)
#
# ---------------------------------------------------------------
# Parameter: force_k
# Default : NULL
# Meaning : Forces the number of DTW clusters.
# Biology : Overrides automatic model selection.
# Usage :
#   - NULL → choose best k by silhouette (recommended)
#   - Integer (e.g. 2, 3, 5) → force that many clusters
# ⚠️ NOTE : If set to a number, silhouette is IGNORED.
#
# ---------------------------------------------------------------
# Parameter: k_min
# Default : 2
# Meaning : Smallest number of clusters tested.
# Biology : Sets the minimum number of distinct temporal programs
#           the model is allowed to consider.
#
# ---------------------------------------------------------------
# Parameter: k_max
# Default : 10
# Meaning : Largest number of clusters tested.
# Biology : Upper bound on how many distinct temporal programs
#           the data can be split into.
# Change if:
#   - Small gene set → LOWER (e.g. 5–6)
#   - Large, heterogeneous gene set → RAISE (e.g. 12–15)
#
# ---------------------------------------------------------------
# Parameter: dtw_distance_method
# Default : "dtw_basic"
# Meaning : Distance metric used for Dynamic Time Warping.
# Biology : Controls how strictly trajectories must align in time.
# Notes   : "dtw_basic" is robust and recommended.
#
# ---------------------------------------------------------------
# Parameter: use_parallel
# Default : FALSE
# Meaning : Enables parallel computation for DTW distance matrix.
# Biology : No effect on results; only affects runtime.
# Change if:
#   - Large gene sets AND many CPU cores available → TRUE
#
# ---------------------------------------------------------------
# Parameter: make_quarto_report
# Default : FALSE
# Meaning : Writes a Quarto (.qmd) report template.
# Biology : No effect on analysis; documentation only.
#
# ---------------------------------------------------------------
# INTERNAL (DO NOT CHANGE UNLESS YOU KNOW WHY):
# - Z-scoring of per-gene trajectories before DTW ensures
#   clustering is based on TEMPORAL SHAPE, not expression magnitude.
# ===============================================================

# ===============================================================

# ---------------------------
# USER SETTINGS (edit if desired)
# ---------------------------
variance_threshold <- 2       # keep genes with variance > this (computed across SAMPLES, before DE)
fdr_threshold      <- 0.05    # adjusted p-value cutoff for DE genes
force_k            <- NULL       # set to NULL for auto (best silhouette), or a number like 5
k_min              <- 2
k_max              <- 10
# ---------------------------
# Permutation test (DTW cluster validity)
# ---------------------------
do_perm_test <- TRUE
n_perm       <- 200          # 200 quick / 500 decent / 1000 strong
perm_seed    <- 123
perm_stat    <- "R"          # "R" (between/within ratio) is recommended
make_perm_plot <- TRUE


# DTW distance choice for clustering
# dtwclust supports several; "dtw_basic" is a safe default.
dtw_distance_method <- "dtw_basic"

# Parallel (optional). Default OFF for biologist-friendly stability.
use_parallel <- TRUE

# Optional Quarto report (OFF by default). If ON, writes a .qmd and tries to render HTML.
make_quarto_report <- TRUE

# ---------------------------
# AUTO-RUN WHEN SOURCED (RStudio)
# ---------------------------
AUTO_RUN_WHEN_SOURCED <- TRUE

# ---------------------------
# Safe package loader/installer
# ---------------------------
need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

need_pkg("limma")
need_pkg("dtw")
need_pkg("dtwclust")
need_pkg("pheatmap")
need_pkg("clusterCrit")
need_pkg("jsonlite")

# Optional for parallel
if (use_parallel) {
  need_pkg("future")
  need_pkg("future.apply")
}

# ---------------------------
# Helper: message + stop formatting
# ---------------------------
say <- function(...) cat(paste0(...), "\n")
die <- function(...) stop(paste0(...), call. = FALSE)

# ---------------------------
# Helper: pick file (mac-friendly)
# ---------------------------
pick_csv_file <- function(prompt_text = "Select your WIDE expression CSV file") {
  say("\n📂 ", prompt_text)
  path <- file.choose()
  if (!nzchar(path)) die("No file selected.")
  path
}

# ---------------------------
# Helper: safe numeric coercion for expression columns
# ---------------------------
coerce_numeric_matrix <- function(df, gene_cols) {
  X <- df[, gene_cols, drop = FALSE]
  for (nm in gene_cols) {
    if (!is.numeric(X[[nm]])) {
      X[[nm]] <- suppressWarnings(as.numeric(X[[nm]]))
    }
  }
  X
}

# ---------------------------
# Core runner
# ---------------------------
run_dtw_timeseries <- function(file_path = NULL) {
  
  # ---- choose file interactively if needed
  if (is.null(file_path) || !nzchar(file_path)) {
    file_path <- pick_csv_file()
  }
  if (!file.exists(file_path)) die("Input file not found: ", file_path)
  
  # ---- output directory next to input file
  base_name  <- sub("\\.csv$", "", basename(file_path), ignore.case = TRUE)
  timestamp  <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  out_root   <- file.path(dirname(file_path), paste0("DTW_TimeSeriesClustering_", base_name, "_", timestamp))
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  # ---- write a tiny README and manifest (provenance)
  manifest <- list(
    input_file = normalizePath(file_path, winslash = "/"),
    output_dir = normalizePath(out_root, winslash = "/"),
    timestamp  = timestamp,
    params     = list(
      variance_threshold = variance_threshold,
      fdr_threshold      = fdr_threshold,
      force_k            = force_k,
      k_min              = k_min,
      k_max              = k_max,
      dtw_distance_method= dtw_distance_method,
      use_parallel       = use_parallel,
      make_quarto_report = make_quarto_report
    )
  )
  jsonlite::write_json(manifest, file.path(out_root, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(
    c(
      "DTW Time-Series Clustering Outputs",
      "--------------------------------",
      paste0("Input: ", manifest$input_file),
      paste0("Run time: ", timestamp),
      "",
      "Key files:",
      "- DE_full_limma.csv",
      "- DE_significant_genes.csv",
      "- sig_expr_byDay_mean.csv",
      "- DTW_distance_matrix.csv",
      "- DTW_hclust_dendrogram.png",
      "- DTW_partitional_silhouette_scores.csv",
      "- DTW_partitional_cluster_assignments.csv",
      "- DTW_clustered_heatmap.png",
      "",
      "Notes:",
      "- limma uses sample-level replicates if present (recommended).",
      "- DTW clustering uses per-Day means of significant genes."
    ),
    con = file.path(out_root, "README_outputs.txt")
  )
  
  say("📁 Output directory: ", out_root)
  
  # ---- load CSV (base R)
  say("📥 Reading CSV: ", file_path)
  df <- tryCatch(read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) die("Failed to read CSV: ", e$message))
  
  if (!("Day" %in% names(df))) die("Your CSV must contain a column named 'Day'.")
  
  # ---- clean Day
  df$Day <- suppressWarnings(as.integer(as.character(df$Day)))
  if (any(!is.finite(df$Day))) die("'Day' column contains non-numeric values that could not be converted.")
  
  # ---- identify gene columns
  gene_cols <- setdiff(names(df), "Day")
  if (length(gene_cols) < 2) die("I only found ", length(gene_cols), " non-Day columns. Need gene columns to proceed.")
  
  # ---- coerce expression columns to numeric
  expr_df <- coerce_numeric_matrix(df, gene_cols)
  bad_frac <- mean(!is.finite(as.matrix(expr_df)))
  if (bad_frac > 0) {
    say("⚠️ Note: ", round(100 * bad_frac, 2), "% of expression values are NA/Inf after numeric coercion.")
  }
  
  # ---- drop rows that are entirely missing across genes
  row_ok <- rowSums(is.finite(as.matrix(expr_df))) > 0
  if (!all(row_ok)) {
    say("⚠️ Dropping ", sum(!row_ok), " rows with no finite gene values.")
    df      <- df[row_ok, , drop = FALSE]
    expr_df <- expr_df[row_ok, , drop = FALSE]
  }
  
  # ============================================================
  # 1) VARIANCE FILTER (across SAMPLES)
  # ============================================================
  say("🧹 Variance filtering genes (threshold = ", variance_threshold, ") ...")
  gene_var <- apply(expr_df, 2, var, na.rm = TRUE)
  keep_genes <- names(gene_var)[is.finite(gene_var) & gene_var > variance_threshold]
  
  say("✅ Genes retained: ", length(keep_genes), " / ", ncol(expr_df))
  if (length(keep_genes) < 2) die("Too few genes passed variance filtering. Lower variance_threshold.")
  
  expr_filt <- expr_df[, keep_genes, drop = FALSE]
  
#============================================================
  # 2) limma DE ACROSS DAYS (uses replicates if present)
  # ============================================================
  expr_mat_gene_by_sample <- t(as.matrix(expr_filt))
  
  day_factor <- factor(df$Day)
  reps_per_day <- table(day_factor)
  say("🔁 Replicates detected per day: ", paste(names(reps_per_day), reps_per_day, sep="=", collapse=", "))
  
  if (all(reps_per_day == 1)) {
    die("No biological replicates detected (each Day appears once). ",
        "You likely supplied a Day-aggregated file. Use the sample-level file instead.")
  }
  
  design <- stats::model.matrix(~ 0 + day_factor)
  colnames(design) <- levels(day_factor)
  
  reps_per_day <- table(day_factor)
  if (any(reps_per_day < 2)) {
    say("⚠️ Some days have <2 samples (replicates). limma will run, but DE inference may be unstable.")
    say("   Replicates per day: ", paste(names(reps_per_day), reps_per_day, sep = "=", collapse = ", "))
  }
  
  say("🧪 Running limma (time-series / overall F-test across days) ...")
  fit <- limma::lmFit(expr_mat_gene_by_sample, design)
  fit <- limma::eBayes(fit)
  
  # IMPORTANT:
  # - topTableF() is deprecated and its API changed.
  # - Use topTable() sorted by overall F statistic.
  tt <- limma::topTable(
    fit,
    number = Inf,
    adjust.method = "BH",
    sort.by = "F"
  )
  
  write.csv(tt, file.path(out_root, "DE_full_limma.csv"), row.names = TRUE)
  say("💾 Saved: DE_full_limma.csv")
  
  sig_tt <- tt[is.finite(tt$adj.P.Val) & tt$adj.P.Val < fdr_threshold, , drop = FALSE]
  write.csv(sig_tt, file.path(out_root, "DE_significant_genes.csv"), row.names = TRUE)
  say("✅ Significant genes (FDR < ", fdr_threshold, "): ", nrow(sig_tt))
  
  if (nrow(sig_tt) < 2) die("Fewer than 2 significant genes. DTW clustering skipped.")
  
  sig_genes <- rownames(sig_tt)
  
  
  # ============================================================
  # 3) BUILD PER-DAY MEAN TIME SERIES (for DTW)
  # ============================================================
  # We build a gene x day matrix using the ORIGINAL sample-level data, averaging replicates within Day.
  say("🧾 Building per-Day mean profiles for significant genes (used for DTW) ...")
  
  unique_days <- sort(unique(df$Day))
  sig_expr_byDay <- matrix(NA_real_, nrow = length(sig_genes), ncol = length(unique_days))
  rownames(sig_expr_byDay) <- sig_genes
  colnames(sig_expr_byDay) <- as.character(unique_days)
  
  for (d in unique_days) {
    idx <- which(df$Day == d)
    # average across samples for each gene
    sig_expr_byDay[, as.character(d)] <- colMeans(expr_filt[idx, sig_genes, drop = FALSE], na.rm = TRUE)
      }
  
  # save for inspection
  write.csv(sig_expr_byDay, file.path(out_root, "sig_expr_byDay_mean.csv"), row.names = TRUE)
  say("💾 Saved: sig_expr_byDay_mean.csv")
  
  # ============================================================
  # FILTER GENES WITH NEGLIGIBLE ABSOLUTE EXPRESSION
  # ============================================================
  say("🧹 Filtering genes with negligible absolute expression...")
  
  min_max_expr <- 20   # adjust if needed; 50–200 typical for RNA-seq scale
  gene_max <- apply(sig_expr_byDay, 1, max, na.rm = TRUE)
  
  keep_expr <- gene_max >= min_max_expr
  
  say("   Removing ", sum(!keep_expr),
      " genes with max expression < ", min_max_expr)
  
  sig_expr_byDay <- sig_expr_byDay[keep_expr, , drop = FALSE]
  sig_genes      <- sig_genes[keep_expr]
  
  
  # ============================================================
  # SAVE LIST OF FILTERED-OUT LOW-EXPRESSION GENES
  # ============================================================
  filtered_out_genes <- names(gene_max)[!keep_expr]
  
  write.csv(
    data.frame(
      Gene = filtered_out_genes,
      MaxExpression = gene_max[!keep_expr]
    ),
    file.path(out_root, "Filtered_LowExpression_Genes.csv"),
    row.names = FALSE
  )
  
  say("📄 Saved: Filtered_LowExpression_Genes.csv")
  
  
  # ============================================================
  # 4) DTW DISTANCE MATRIX (gene-gene)
  # ============================================================
  say("🔄 Computing DTW distance matrix (this can be slow if many genes) ...")
  n_genes <- nrow(sig_expr_byDay)
  if (n_genes > 2000) {
    say("⚠️ You have ", n_genes, " significant genes. DTW all-pairs may be very slow/heavy.")
    say("   Consider using a stricter FDR or higher variance_threshold.")
  }
  
  # pairwise DTW distance (symmetric)
  # We'll compute upper triangle and fill.
  dtw_dist <- matrix(0, n_genes, n_genes)
  rownames(dtw_dist) <- rownames(sig_expr_byDay)
  colnames(dtw_dist) <- rownames(sig_expr_byDay)
  
  # optional parallel
  if (use_parallel) {
    future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
    on.exit(future::plan(future::sequential), add = TRUE)
    say("🚀 Parallel ON (future multisession).")
  } else {
    say("🧠 Parallel OFF (stable default).")
  }
  
  # build index pairs
  pairs <- combn(seq_len(n_genes), 2)
  pair_list <- split(pairs, rep(seq_len(ncol(pairs)), each = nrow(pairs)))
  
  dtw_one <- function(ii_jj) {
    i <- ii_jj[1]; j <- ii_jj[2]
    a <- sig_expr_byDay[i, ]
    b <- sig_expr_byDay[j, ]
    # dtw::dtw returns alignment; distance.only keeps it light
    dtw::dtw(a, b, distance.only = TRUE)$distance
  }
  
  if (use_parallel) {
    dists <- future.apply::future_sapply(pair_list, dtw_one)
  } else {
    dists <- sapply(pair_list, dtw_one)
  }
  
  # fill matrix
  k <- 1
  for (col in seq_len(ncol(pairs))) {
    i <- pairs[1, col]; j <- pairs[2, col]
    dtw_dist[i, j] <- dists[k]
    dtw_dist[j, i] <- dists[k]
    k <- k + 1
  }
  
  write.csv(dtw_dist, file.path(out_root, "DTW_distance_matrix.csv"), row.names = TRUE)
  say("💾 Saved: DTW_distance_matrix.csv")
  
  # ============================================================
  # 5) HIERARCHICAL CLUSTERING + DENDROGRAM
  # ============================================================
  say("🌳 Hierarchical clustering (average linkage) ...")
  hc <- hclust(as.dist(dtw_dist), method = "average")
  
  png(file.path(out_root, "DTW_hclust_dendrogram.png"), width = 1200, height = 800)
  plot(hc, main = "Hierarchical Clustering (DTW distance)", xlab = "", sub = "", cex = 0.4)
  dev.off()
  say("🖼️ Saved: DTW_hclust_dendrogram.png")
  
  # ============================================================
  # 6) PARTITIONAL CLUSTERING (dtwclust) + SILHOUETTE SELECTION
  # ============================================================
  say("🧩 Partitional clustering with dtwclust; testing k = ", k_min, " .. ", k_max, " ...")
  
  ks <- k_min:k_max
  sil_scores <- rep(NA_real_, length(ks))
  names(sil_scores) <- as.character(ks)
  models <- vector("list", length(ks))
  names(models) <- as.character(ks)
  
  # dtwclust expects series as rows (objects) and columns (time)
  # We already have gene x day => perfect.
  # ============================================================
  # SHAPE-ONLY NORMALIZATION FOR DTW (RECOMMENDED)
  # ============================================================
  say("🔄 Z-scoring gene trajectories for shape-only DTW clustering...")
  
  X <- t(scale(t(sig_expr_byDay)))
  X[!is.finite(X)] <- 0
  
  
  for (ii in seq_along(ks)) {
    k_try <- ks[ii]
    say("   • k = ", k_try)
    
    mod <- tryCatch(
      dtwclust::tsclust(
        X,
        type     = "partitional",
        k        = k_try,
        distance = dtw_distance_method,
        seed     = 123,
        trace    = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(mod)) {
      sil_scores[ii] <- NA_real_
      models[[ii]] <- NULL
      next
    }
    
    cl <- as.integer(mod@cluster)
    
    # silhouette from clusterCrit requires matrix (objects x features) and cluster labels
    # Our objects=genes, features=days.
    ok <- length(unique(cl)) > 1 && all(table(cl) >= 2)
    if (!ok) {
      sil_scores[ii] <- NA_real_
      models[[ii]] <- mod
      next
    }
    
    sil <- tryCatch(
      clusterCrit::intCriteria(as.matrix(X), cl, c("Silhouette"))$silhouette,
      error = function(e) NA_real_
    )
    
    sil_scores[ii] <- sil
    models[[ii]] <- mod
  }
  
  sil_df <- data.frame(
    k = ks,
    silhouette = as.numeric(sil_scores)
  )
  write.csv(sil_df, file.path(out_root, "DTW_partitional_silhouette_scores.csv"), row.names = FALSE)
  say("💾 Saved: DTW_partitional_silhouette_scores.csv")
  
  # silhouette plot
  png(file.path(out_root, "DTW_partitional_silhouette_scores.png"), width = 900, height = 600)
  plot(sil_df$k, sil_df$silhouette, type = "b",
       xlab = "k (number of clusters)", ylab = "Silhouette (higher is better)",
       main = "DTW Partitional Clustering: Silhouette by k")
  dev.off()
  say("🖼️ Saved: DTW_partitional_silhouette_scores.png")
  
  # choose best k (or force)
  if (!is.null(force_k)) {
    best_k <- force_k
    say("⚠️ force_k set: using k = ", best_k)
  } else {
    if (all(!is.finite(sil_df$silhouette))) die("All silhouette scores are NA. Try different k range or check data.")
    best_k <- sil_df$k[which.max(sil_df$silhouette)]
    say("🏆 Best k by silhouette: ", best_k)
  }
  
  best_model <- models[[as.character(best_k)]]
  if (is.null(best_model)) die("Best model is NULL (k=", best_k, "). Try force_k to a different value.")
  
  cluster_assignments <- data.frame(
    Gene = rownames(X),
    Cluster = as.integer(best_model@cluster),
    stringsAsFactors = FALSE
  )
  write.csv(cluster_assignments, file.path(out_root, "DTW_partitional_cluster_assignments.csv"), row.names = FALSE)
  say("💾 Saved: DTW_partitional_cluster_assignments.csv")
  
  # ============================================================
  # 6b) PERMUTATION TEST: Are DTW clusters better than chance?
  #     (Label permutation; preserves cluster sizes)
  # ============================================================
  if (isTRUE(do_perm_test)) {
    say("🧪 Permutation test (label randomization) for DTW clustering...")
    
    # ---- sanity checks
    if (!exists("dtw_dist")) die("Permutation test needs dtw_dist in memory.")
    if (!is.matrix(dtw_dist) || nrow(dtw_dist) != ncol(dtw_dist)) die("dtw_dist must be a square matrix.")
    if (nrow(cluster_assignments) < 4) die("Too few genes for permutation test.")
    
    # Align labels to dtw_dist order
    genes <- rownames(dtw_dist)
    lab_df <- cluster_assignments
    if (!all(c("Gene", "Cluster") %in% names(lab_df))) die("cluster_assignments must have columns Gene and Cluster.")
    
    m <- match(genes, lab_df$Gene)
    if (any(!is.finite(m))) {
      missing <- genes[!is.finite(m)][1:min(10, sum(!is.finite(m)))]
      die("Permutation test: some genes in dtw_dist missing from cluster_assignments. Examples: ",
          paste(missing, collapse = ", "))
    }
    cl_obs <- as.integer(lab_df$Cluster[m])
    
    # This permutation test is most interpretable for k=2
    ucl <- sort(unique(cl_obs))
    if (length(ucl) != 2) {
      say("⚠️ Permutation test block currently targets k=2. You have k=", length(ucl), ".")
      say("   (You can still run it, but interpretation changes.)")
    }
    
    # Helper to compute within/between and ratio, using dtw_dist
    # Uses upper triangle pairs only (efficient, no double counting)
    pair_stats <- function(labels, D) {
      n <- length(labels)
      if (n < 4) return(c(W = NA_real_, B = NA_real_, R = NA_real_))
      
      # Upper triangle indices
      ut <- which(upper.tri(D), arr.ind = TRUE)
      i <- ut[, 1]; j <- ut[, 2]
      dij <- D[ut]
      
      same <- labels[i] == labels[j]
      W <- mean(dij[same], na.rm = TRUE)
      B <- mean(dij[!same], na.rm = TRUE)
      
      # Ratio (within / between); lower = better separation
      R <- W / B
      c(W = W, B = B, R = R)
    }
    
    # Observed statistics
    obs <- pair_stats(cl_obs, dtw_dist)
    
    # ============================================================
    # NEW: Save OBSERVED W/B/R to disk (this was missing)
    # ============================================================
    obs_df <- data.frame(
      perm_id = 0,
      W = as.numeric(obs["W"]),
      B = as.numeric(obs["B"]),
      R = as.numeric(obs["R"])
    )
    
    write.csv(
      obs_df,
      file.path(out_root, "DTW_permtest_observed.csv"),
      row.names = FALSE
    )
    
    say("✅ Observed stats: W=", signif(obs["W"], 6),
        "  B=", signif(obs["B"], 6),
        "  R=", signif(obs["R"], 6))
    
    
    # Permutations: shuffle labels while preserving cluster sizes
    set.seed(perm_seed)
    
    # Preserve exact counts
    tab <- table(cl_obs)
    labs <- rep(as.integer(names(tab)), times = as.integer(tab))  # e.g., c(1,1,1,2,2,...)
    if (length(labs) != length(cl_obs)) die("Internal error building size-preserving labels.")
    
    perm_W <- numeric(n_perm)
    perm_B <- numeric(n_perm)
    perm_R <- numeric(n_perm)
    
    for (pp in seq_len(n_perm)) {
      cl_p <- sample(labs, size = length(labs), replace = FALSE)
      st <- pair_stats(cl_p, dtw_dist)
      perm_W[pp] <- st["W"]
      perm_B[pp] <- st["B"]
      perm_R[pp] <- st["R"]
      if (pp %% 50 == 0) say("   Perm ", pp, " / ", n_perm)
    }
    
    perm_df <- data.frame(
      perm_id = seq_len(n_perm),
      W = perm_W,
      B = perm_B,
      R = perm_R
    )
    # ============================================================
    # NEW: prepend observed row to permutation table (single self-contained CSV)
    # ============================================================
    perm_df_with_obs <- rbind(obs_df, perm_df)
    
    # Empirical p-values (one-sided: "as good or better than observed")
    # For separation, higher R is better; lower W is better.
    p_R <- (1 + sum(perm_df$R <= obs["R"], na.rm = TRUE)) / (n_perm + 1)
    p_W <- (1 + sum(perm_df$W <= obs["W"], na.rm = TRUE)) / (n_perm + 1)
    
    # Save results
    out_csv <- file.path(out_root, "DTW_permtest_label_based.csv")
    write.csv(perm_df_with_obs, out_csv, row.names = FALSE)
    
    out_txt <- file.path(out_root, "DTW_permtest_summary.txt")
    writeLines(c(
      "DTW permutation test (label randomization; size-preserving)",
      "----------------------------------------------------------",
      paste0("n_perm: ", n_perm),
      paste0("seed: ", perm_seed),
      "",
      "Observed statistics:",
      paste0("  Within mean distance (W): ", signif(obs["W"], 6)),
      paste0("  Between mean distance (B): ", signif(obs["B"], 6)),
      paste0("  Separation ratio (R=W/B):  ", signif(obs["R"], 6)),
      "",
      "Empirical p-values (one-sided):",
      paste0("  p_R (perm R <= obs R): ", signif(p_R, 4)),
      paste0("  p_W (perm W <= obs W): ", signif(p_W, 4)),
      "",
      "Interpretation:",
      "- Small p_R suggests the observed clusters are more separated than random labelings.",
      "- Large p_R suggests the 2-cluster partition is not more separated than chance (in DTW space).",
      "- p_W is an alternative view focusing on within-cluster tightness."
    ), con = out_txt)
    
    say("💾 Saved: DTW_permtest_label_based.csv")
    say("📝 Saved: DTW_permtest_summary.txt")
    
    # Plot null distribution with observed line
    if (isTRUE(make_perm_plot)) {
      png(file.path(out_root, "DTW_permtest_R_hist.png"), width = 900, height = 600)
      hist(
        perm_df$R,
        breaks = 40,
        main = "Permutation null: DTW ratio R = (within / between)",
        xlab = "R = W/B (lower = better separation)"
      )
      abline(v = obs["R"], lwd = 3)
      legend("topright",
             legend = c(paste0("Observed R = ", signif(obs["R"], 4)),
                        paste0("Empirical p_R = ", signif(p_R, 3))),
             bty = "n")
      dev.off()
      say("🖼️ Saved: DTW_permtest_R_hist.png")
    }
  }
  
  
  X <- read.csv(file.path(out_root, "sig_expr_byDay_mean.csv"), row.names = 1, check.names = FALSE)
  cl <- read.csv(file.path(out_root, "DTW_partitional_cluster_assignments.csv"), stringsAsFactors = FALSE)
  
  # per-gene expression summaries across days
  gene_mean <- rowMeans(X, na.rm = TRUE)
  gene_max  <- apply(X, 1, max, na.rm = TRUE)
  
  tmp <- data.frame(Gene = rownames(X), mean = gene_mean, max = gene_max)
  tmp <- merge(tmp, cl, by = "Gene")
  
  # show lowest-expression genes overall
  tmp[order(tmp$max), ][1:20, ]
  
  # show cluster medians (this tells you if an entire cluster is low)
  tapply(tmp$max, tmp$Cluster, median, na.rm = TRUE)
  tapply(tmp$mean, tmp$Cluster, median, na.rm = TRUE)
  
  # ============================================================
  # SAVE DTW CLUSTER CENTROIDS (per-day trajectories)
  # ============================================================
  say("📈 Saving DTW cluster centroids...")
  
  centroid_list <- best_model@centroids
  
  centroid_mat <- do.call(rbind, centroid_list)
  rownames(centroid_mat) <- paste0("Cluster_", seq_len(nrow(centroid_mat)))
  colnames(centroid_mat) <- colnames(sig_expr_byDay)
  
  centroid_df <- as.data.frame(centroid_mat)
  centroid_df$Cluster <- rownames(centroid_df)
  
  # reorder columns: Cluster first
  centroid_df <- centroid_df[, c("Cluster", colnames(sig_expr_byDay))]
  
  write.csv(
    centroid_df,
    file.path(out_root, "DTW_cluster_centroids_byDay.csv"),
    row.names = FALSE
  )
  
  say("💾 Saved: DTW_cluster_centroids_byDay.csv")
  # ============================================================
  # PLOT DTW CLUSTER CENTROIDS (colored by cluster)
  # ============================================================
  say("🖼️ Plotting DTW cluster centroids...")
  
  # centroid_mat is clusters x days
  # ensure days are numeric and ordered correctly
  day_vals <- suppressWarnings(as.numeric(colnames(centroid_mat)))
  if (any(!is.finite(day_vals))) {
    # fallback: use index positions if day labels aren't numeric
    day_vals <- seq_len(ncol(centroid_mat))
  }
  
  # order by day (in case columns are out of order)
  ord <- order(day_vals)
  day_vals <- day_vals[ord]
  centroid_mat_plot <- centroid_mat[, ord, drop = FALSE]
  
  nC <- nrow(centroid_mat_plot)
  
  # simple distinct palette (base R)
  cols <- grDevices::hcl.colors(nC, palette = "Dark 3")
  
  png(file.path(out_root, "DTW_cluster_centroids_byDay.png"), width = 1000, height = 700)
  
  # set y-range across all centroids
  yl <- range(centroid_mat_plot, finite = TRUE)
  
  # plot first centroid to initialize canvas
  plot(
    day_vals,
    centroid_mat_plot[1, ],
    type = "l",
    lwd = 3,
    col = cols[1],
    xlab = "Day",
    ylab = "Centroid expression (DTW prototype)",
    main = "DTW Cluster Centroids (per-day trajectories)",
    ylim = yl
  )
  
  # add remaining centroids
  if (nC > 1) {
    for (i in 2:nC) {
      lines(day_vals, centroid_mat_plot[i, ], lwd = 3, col = cols[i])
    }
  }
  
  legend(
    "topleft",
    legend = rownames(centroid_mat_plot),
    col = cols,
    lwd = 3,
    bty = "n",
    cex = 0.9
  )
  
  dev.off()
  
  say("💾 Saved: DTW_cluster_centroids_byDay.png")
  
  # ============================================================
  # 7) HEATMAP (clustered genes ordered by cluster)
  # ============================================================
  say("🎨 Heatmap (genes ordered by cluster) ...")
  ord <- order(cluster_assignments$Cluster)
  X_ord <- X[ord, , drop = FALSE]
  ann <- data.frame(Cluster = factor(cluster_assignments$Cluster[ord]))
  rownames(ann) <- rownames(X_ord)
  
  pheatmap::pheatmap(
    X_ord,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    annotation_row = ann,
    show_rownames = FALSE,
    main = "Significant genes (per-Day mean), ordered by DTW clusters",
    filename = file.path(out_root, "DTW_clustered_heatmap.png"),
    width = 9,
    height = 10
  )
  say("🖼️ Saved: DTW_clustered_heatmap.png")
  
  # ============================================================
  # 8) OPTIONAL: Quarto report (.qmd + optional HTML render)
  # ============================================================
  if (make_quarto_report) {
    
    # write QMD
    qmd <- file.path(out_root, "DTW_TimeSeries_Report.qmd")
    writeLines(c(
      "---",
      "title: \"DTW Time-Series Clustering Report\"",
      "format: html",
      "---",
      "",
      "## Run summary",
      "",
      paste0("- Input: `", basename(file_path), "`"),
      paste0("- Output folder: `", basename(out_root), "`"),
      paste0("- variance_threshold: ", variance_threshold),
      paste0("- fdr_threshold: ", fdr_threshold),
      paste0("- min_max_expr: ", ifelse(exists("min_max_expr", inherits = TRUE),
                                        min_max_expr, "NA")),
      paste0("- force_k: ", ifelse(is.null(force_k), "NULL", force_k)),
      paste0("- Selected k: ", best_k),
      "",
      "## Files produced",
      "",
      "- `DE_full_limma.csv`",
      "- `DE_significant_genes.csv`",
      "- `sig_expr_byDay_mean.csv`",
      "- `DTW_distance_matrix.csv`",
      "- `DTW_hclust_dendrogram.png`",
      "- `DTW_partitional_silhouette_scores.png`",
      "- `DTW_partitional_cluster_assignments.csv`",
      "- `DTW_cluster_centroids_byDay.csv`",
      "- `DTW_cluster_centroids_byDay.png`",
      "- `DTW_clustered_heatmap.png`",
      "- `Filtered_LowExpression_Genes.csv` (if enabled)",
      "",
      "## Key figures",
      "",
      "![](DTW_cluster_centroids_byDay.png)",
      "",
      "![](DTW_partitional_silhouette_scores.png)"
      ), qmd)
    
    say("📝 Wrote Quarto stub: ", qmd)
    
    # optional HTML render
    if (nzchar(Sys.which("quarto"))) {
      
      # R package 'quarto' is needed for quarto::quarto_render()
      need_pkg("quarto")
      
      say("🧾 Rendering Quarto HTML report...")
      try(
        quarto::quarto_render(
          input = qmd,
          output_format = "html"
        ),
        silent = TRUE
      )
      
      # tell user where the HTML should be
      html_out <- sub("\\.qmd$", ".html", qmd, ignore.case = TRUE)
      if (file.exists(html_out)) {
        say("✅ Rendered HTML: ", html_out)
      } else {
        say("⚠️ Quarto render attempted but HTML was not found. Check Quarto install.")
      }
      
    } else {
      say("⚠️ Quarto not found on system PATH — skipping HTML render.")
      say("   Install Quarto (system app), then re-run with make_quarto_report <- TRUE.")
    }
  }
  
  
  # ---- inventory file list
  inv <- list.files(out_root, full.names = FALSE)
  writeLines(inv, con = file.path(out_root, "file_inventory.txt"))
  
  # ---- final summary
  say("\n📋 FINAL SUMMARY")
  say("• Input file: ", file_path)
  say("• Output dir: ", out_root)
  say("• Genes after variance filter: ", length(keep_genes))
  say("• Significant genes (FDR<", fdr_threshold, "): ", length(sig_genes))
  say("• Clusters (k): ", best_k)
  
  invisible(list(
    out_dir = out_root,
    keep_genes = keep_genes,
    sig_genes = sig_genes,
    best_k = best_k,
    cluster_assignments = cluster_assignments
  ))
}

# ---------------------------
# RUN WHEN SOURCED (if interactive)
# ---------------------------
if (interactive() && isTRUE(AUTO_RUN_WHEN_SOURCED)) {
  say("Running DTW time-series clustering interactively...")
  # If you want to hardcode a path, set it here:
  # result <- run_dtw_timeseries(file_path = "/path/to/your.csv")
  result <- run_dtw_timeseries(file_path = NULL)
  say("Done. Output written to:\n", result$out_dir)
  # Open folder in Finder (mac) or Explorer (Windows) if you want:
  # mac: system(paste("open", shQuote(result$out_dir)))
  # win: shell.exec(result$out_dir)
}
