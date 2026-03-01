#!/usr/bin/env Rscript
# ======================================================================
# RV_MDS_ANALYSIS.R
#Day–Day Correlation-Structure Similarity: NON-HEATMAP Publication Plots
# + Integrated "Optimized" RV Bootstrap (row-resampling; memory-conscious)
#
# PRIMARY INPUT (CSV; long pairwise): Day1, Day2, matrix_corr, frobenius, RV
# OPTIONAL INPUT (for RV bootstrap only; tidy expression): Day, SampleID, Gene, Abundance
#
# OUTPUT (timestamped folder):
#   - similarity_long_clean.csv
#   - similarity_wide_matrix_corr.csv
#   - similarity_wide_RV.csv
#   - similarity_wide_frobenius.csv
#   - RefDay_Profile_<metric>_data.csv
#   - Plot_RefDay_Profile_<metric>.png/.pdf
#   - Plot_MDS_<metric>.png/.pdf
#   - Plot_DayNetwork_<metric>.png/.pdf
#   - mds_coords_<metric>.csv
#   - network_edges_<metric>_topk.csv
#   - RefDay_Profile_RV_bootstrap_CI.csv              (if metric=RV and bootstrap runs)
#   - RV_bootstrap_draws_profile_wide.csv             (optional; if --save_boot_draws=TRUE)
#
# USAGE (pairwise only):
#   Rscript day_similarity_pretty.R --in=A_GlobalSimilarity_ByDayPairs.csv --metric=RV --ref_day=14 --top_k=4
#
# USAGE (with bootstrap):
#   Rscript day_similarity_pretty.R --in=A_GlobalSimilarity_ByDayPairs.csv --metric=RV --ref_day=14 --top_k=4 \
#          --expr=Transcriptome_Tidy_Long.csv --nboot=200 --seed=1 --ci=0.95 --transform=log1p --max_genes=0
#
# NOTES ON BOOTSTRAP MODE:
#   This integrates an "optimized" bootstrap that resamples ROWS of each day's
#   Sample×Gene matrix (after preprocessing). This is fast and memory-conscious.
#   It is a pragmatic robustness check, but it is NOT identical to resampling
#   biological replicates at the SampleID level if samples differ in gene missingness.
#
# IMPORTANT:
#   All if/else logic is written in BLOCK FORM to avoid R's "unexpected else" errors.
# ======================================================================

# ------------------------------- Helpers --------------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2    <- function(...) stop(paste0(...), call. = FALSE)
msg      <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

ensure_pkgs <- function(pkgs) {
  missing <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(missing) > 0) {
    msg("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

arg_val <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

pick_file <- function(prompt = "Choose file") {
  if (interactive() && exists("file.choose")) {
    msg(prompt)
    return(file.choose())
  }
  stop2(prompt, " (non-interactive). Provide a path via command line args.")
}

ask_int <- function(prompt, default) {
  if (!interactive()) return(as.integer(default))
  ans <- readline(paste0(prompt, " [default ", default, "]: "))
  if (nchar(ans) == 0) return(as.integer(default))
  out <- suppressWarnings(as.integer(ans))
  if (is.na(out) || out < 1) return(as.integer(default))
  out
}

ask_yesno <- function(prompt, default_yes = TRUE) {
  if (!interactive()) return(default_yes)
  def <- if (default_yes) "Y/n" else "y/N"
  ans <- readline(paste0(prompt, " [", def, "]: "))
  if (nchar(ans) == 0) return(default_yes)
  ans <- tolower(ans)
  if (ans %in% c("y", "yes")) return(TRUE)
  if (ans %in% c("n", "no")) return(FALSE)
  default_yes
}

is_true <- function(x) {
  if (is.null(x) || is.na(x)) return(FALSE)
  x <- tolower(as.character(x))
  x %in% c("1", "true", "t", "yes", "y")
}

# ------------------------------- Params ---------------------------------

in_path <- arg_val("--in", default = NA_character_)
metric  <- arg_val("--metric", default = "RV")     # RV | matrix_corr | frobenius
ref_day <- suppressWarnings(as.numeric(arg_val("--ref_day", default = "14")))
top_k   <- suppressWarnings(as.integer(arg_val("--top_k", default = "4")))
out_dir <- arg_val("--outdir", default = file.path(getwd(), paste0("DaySimilarityPlots_", ts_stamp())))

# Bootstrap params (RV only)
expr_path         <- arg_val("--expr", default = NA_character_)
nboot_arg         <- arg_val("--nboot", default = NA_character_)
seed              <- suppressWarnings(as.integer(arg_val("--seed", default = "1")))
ci_level          <- suppressWarnings(as.numeric(arg_val("--ci", default = "0.95")))
transform         <- arg_val("--transform", default = "log1p")  # none | log1p
max_genes         <- suppressWarnings(as.integer(arg_val("--max_genes", default = "0"))) # 0 = no cap
run_boot_arg      <- arg_val("--run_boot", default = NA_character_) # TRUE/FALSE; if absent, prompt
save_boot_draws   <- is_true(arg_val("--save_boot_draws", default = "FALSE")) # save wide draws CSV

# Debug toggles
debug            <- is_true(arg_val("--debug", default = "TRUE"))
save_debug_files <- is_true(arg_val("--save_debug_files", default = "TRUE"))

# Resolve inputs
if (is.na(in_path) || nchar(in_path) == 0) in_path <- pick_file("Choose pairwise Day1–Day2 similarity CSV (Day1, Day2, matrix_corr, frobenius, RV)")
if (!file.exists(in_path)) stop2("Input file not found: ", in_path)

if (tolower(metric) == "rv") metric <- "RV"
if (!metric %in% c("matrix_corr", "RV", "frobenius")) stop2("--metric must be one of: RV, matrix_corr, frobenius")

if (ci_level <= 0 || ci_level >= 1) stop2("--ci must be between 0 and 1 (e.g., 0.95).")
alpha <- (1 - ci_level) / 2

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

msg("Input:   ", in_path)
msg("Metric:  ", metric)
msg("Ref day: ", ref_day)
msg("Top_k:   ", top_k)
msg("Output:  ", out_dir)
msg("Debug:   ", debug, " | save_debug_files=", save_debug_files)

# ----------------------------- Dependencies -----------------------------

ensure_pkgs(c("ggplot2", "igraph"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(igraph)
})

theme_pub <- theme_classic(base_size = 13) +
  theme(
    axis.title    = element_text(face = "bold"),
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.title  = element_text(face = "bold"),
    panel.grid    = element_blank()
  )

# ------------------------------ Read + QC -------------------------------

df <- read.csv(in_path, stringsAsFactors = FALSE, check.names = FALSE)

req <- c("Day1", "Day2", "matrix_corr", "frobenius", "RV")
missing_cols <- setdiff(req, colnames(df))
if (length(missing_cols) > 0) stop2("Missing required columns: ", paste(missing_cols, collapse = ", "))

df$Day1 <- suppressWarnings(as.numeric(df$Day1))
df$Day2 <- suppressWarnings(as.numeric(df$Day2))
if (any(is.na(df$Day1)) || any(is.na(df$Day2))) stop2("Day1/Day2 contain non-numeric values after coercion.")

for (m in c("matrix_corr", "frobenius", "RV")) {
  df[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  if (any(is.na(df[[m]]))) stop2("Metric '", m, "' has non-numeric values or NA.")
}

# Deduplicate by averaging (safe)
key <- paste(df$Day1, df$Day2, sep = "_")
if (any(duplicated(key))) {
  msg("Duplicate Day1-Day2 pairs found; averaging duplicates.")
  df <- aggregate(df[, c("matrix_corr", "frobenius", "RV")],
                  by = list(Day1 = df$Day1, Day2 = df$Day2),
                  FUN = mean)
}

# Symmetrize
df_rev <- df
df_rev$Day1 <- df$Day2
df_rev$Day2 <- df$Day1
df_sym <- rbind(df, df_rev)

days <- sort(unique(c(df_sym$Day1, df_sym$Day2)))

# Add diagonals
diag_df <- data.frame(Day1 = days, Day2 = days, matrix_corr = 1, RV = 1, frobenius = 0)
df_sym <- rbind(df_sym, diag_df)

# Export cleaned long
write.csv(df_sym, file.path(out_dir, "similarity_long_clean.csv"), row.names = FALSE)

# ------------------------------ Wide matrices ---------------------------

to_wide <- function(d, value_col, diag_value, days_vec) {
  mat <- matrix(NA_real_,
                nrow = length(days_vec),
                ncol = length(days_vec),
                dimnames = list(as.character(days_vec), as.character(days_vec)))
  for (i in seq_len(nrow(d))) {
    r <- as.character(d$Day1[i])
    c <- as.character(d$Day2[i])
    mat[r, c] <- d[[value_col]][i]
  }
  diag(mat) <- diag_value
  mat
}

mat_mc <- to_wide(df_sym, "matrix_corr", 1, days)
mat_rv <- to_wide(df_sym, "RV",         1, days)
mat_fr <- to_wide(df_sym, "frobenius",  0, days)

write.csv(mat_mc, file.path(out_dir, "similarity_wide_matrix_corr.csv"))
write.csv(mat_rv, file.path(out_dir, "similarity_wide_RV.csv"))
write.csv(mat_fr, file.path(out_dir, "similarity_wide_frobenius.csv"))

# Choose matrix by metric (BLOCK if/else)
if (metric == "matrix_corr") {
  M <- mat_mc
} else if (metric == "RV") {
  M <- mat_rv
} else {
  M <- mat_fr
}

# ------------------------------ Plot 1: Reference-day profile -----------

ref_profile_path <- file.path(out_dir, paste0("RefDay_Profile_", metric, "_data.csv"))

if (!is.na(ref_day) && (ref_day %in% days)) {
  
  vals <- as.numeric(M[as.character(days), as.character(ref_day)])
  df_ref <- data.frame(Day = days, Value = vals, Metric = metric, RefDay = ref_day)
  write.csv(df_ref, ref_profile_path, row.names = FALSE)
  
  if (metric == "frobenius") {
    ylab <- "Frobenius distance (lower = more similar)"
  } else {
    ylab <- paste0(metric, " similarity to Day ", ref_day)
  }
  
  p_ref <- ggplot(df_ref, aes(x = Day, y = Value)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_x_continuous(breaks = days) +
    labs(
      title = paste0("Global transcriptome structure relative to Day ", ref_day),
      subtitle = paste0("Metric: ", metric),
      x = "Day",
      y = ylab
    ) +
    theme_pub
  
  ggsave(file.path(out_dir, paste0("Plot_RefDay_Profile_", metric, ".png")), p_ref,
         width = 8.2, height = 4.6, dpi = 300)
  ggsave(file.path(out_dir, paste0("Plot_RefDay_Profile_", metric, ".pdf")), p_ref,
         width = 8.2, height = 4.6)
  
} else {
  msg("Reference day ", ref_day, " not present; skipping reference-day profile plot.")
}

# ------------------------------ Plot 2: MDS embedding -------------------

# Convert similarity -> distance for embedding
if (metric == "frobenius") {
  D <- M
} else {
  S <- pmin(pmax(M, 0), 1)  # clamp to [0,1]
  D <- 1 - S
}

if (any(is.na(D))) stop2("Distance matrix contains NA; missing day pairs in input table.")

mds <- cmdscale(as.dist(D), k = 2, eig = TRUE)
coords <- data.frame(Day = days, Dim1 = mds$points[, 1], Dim2 = mds$points[, 2])
write.csv(coords, file.path(out_dir, paste0("mds_coords_", metric, ".csv")), row.names = FALSE)

p_mds <- ggplot(coords, aes(x = Dim1, y = Dim2, label = Day)) +
  geom_point(size = 2.6) +
  geom_text(vjust = -0.7, size = 4) +
  labs(
    title = paste0("2D embedding of days from global structure (", metric, ")"),
    x = "MDS dimension 1",
    y = "MDS dimension 2"
  ) +
  theme_pub

ggsave(file.path(out_dir, paste0("Plot_MDS_", metric, ".png")), p_mds,
       width = 6.6, height = 5.4, dpi = 300)
ggsave(file.path(out_dir, paste0("Plot_MDS_", metric, ".pdf")), p_mds,
       width = 6.6, height = 5.4)

# ------------------------------ Plot 3: Day network ---------------------

edges <- data.frame()

for (d in days) {
  vals <- M[as.character(d), as.character(days)]
  names(vals) <- as.character(days)
  vals <- vals[names(vals) != as.character(d)]
  
  if (metric == "frobenius") {
    ord <- order(vals, decreasing = FALSE)  # smaller distance = closer
  } else {
    ord <- order(vals, decreasing = TRUE)   # larger similarity = closer
  }
  
  keep <- head(ord, top_k)
  to_days <- as.numeric(names(vals)[keep])
  
  sub <- data.frame(from = d, to = to_days, weight = as.numeric(vals[keep]))
  edges <- rbind(edges, sub)
}

if (metric == "frobenius") {
  edges <- edges[order(edges$weight, decreasing = FALSE), ]
} else {
  edges <- edges[order(edges$weight, decreasing = TRUE), ]
}

# Deduplicate undirected edges (keep best tie)
edges$pair <- apply(edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "_"))
edges <- edges[!duplicated(edges$pair), ]
edges$pair <- NULL

stopifnot(nrow(edges) > 0)

write.csv(edges, file.path(out_dir, paste0("network_edges_", metric, "_topk.csv")), row.names = FALSE)

g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = as.character(days)))

set.seed(1)
lay <- layout_with_fr(g)

png(file.path(out_dir, paste0("Plot_DayNetwork_", metric, ".png")), width = 1400, height = 1000, res = 150)
plot(g, layout = lay,
     vertex.size = 30,
     vertex.label.cex = 1.1,
     edge.width = 1.2,
     main = paste0("Day network from global structure (", metric, "), top_k=", top_k))
dev.off()

pdf(file.path(out_dir, paste0("Plot_DayNetwork_", metric, ".pdf")), width = 9, height = 6.5)
plot(g, layout = lay,
     vertex.size = 22,
     vertex.label.cex = 1.0,
     edge.width = 1.0,
     main = paste0("Day network from global structure (", metric, "), top_k=", top_k))
dev.off()

# ======================================================================
# RV BOOTSTRAP (integrated optimized row-resampling approach)
# ======================================================================

# Build Sample×Gene matrix for a given day (fast-ish; avoids per-cell loops where possible)
build_day_matrix <- function(expr_df, day, trans = "log1p", max_g = 0) {
  d <- expr_df[expr_df$Day == day, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  
  samps <- unique(d$SampleID)
  genes <- unique(d$Gene)
  
  if (max_g > 0 && length(genes) > max_g) {
    genes <- genes[seq_len(max_g)]
    d <- d[d$Gene %in% genes, , drop = FALSE]
  }
  
  # Map to indices; drop rows that do not map cleanly (should not happen, but safe)
  i <- match(d$SampleID, samps)
  j <- match(d$Gene, genes)
  ok <- !(is.na(i) | is.na(j) | is.na(d$Abundance))
  d  <- d[ok, , drop = FALSE]
  i  <- i[ok]
  j  <- j[ok]
  
  mat <- matrix(NA_real_, nrow = length(samps), ncol = length(genes),
                dimnames = list(samps, genes))
  mat[cbind(i, j)] <- d$Abundance
  
  # Impute missing per gene (median) + transform
  for (col in seq_len(ncol(mat))) {
    v <- mat[, col]
    if (all(is.na(v))) {
      # leave all-NA; will be removed by keep-filter later
      next
    }
    med <- median(v, na.rm = TRUE)
    v[is.na(v)] <- med
    if (trans == "log1p") {
      v <- log1p(v)
    } else if (trans == "none") {
      # no-op
    } else {
      stop2("Unknown --transform: ", trans, " (use none|log1p)")
    }
    mat[, col] <- v
  }
  
  # Scale & drop bad columns (NA/Inf or zero variance)
  mat <- scale(mat, center = TRUE, scale = TRUE)
  
  keep <- apply(mat, 2, function(x) {
    if (!all(is.finite(x))) return(FALSE)
    if (stats::sd(x) <= 0) return(FALSE)
    TRUE
  })
  
  mat <- mat[, keep, drop = FALSE]
  if (is.null(dim(mat)) || ncol(mat) < 2 || nrow(mat) < 2) return(NULL)
  mat
}

# RV between two day matrices via gene–gene correlation matrices (DIMENSION SAFE)
calc_rv <- function(matA, matB) {
  if (is.null(matA) || is.null(matB)) return(NA_real_)
  
  # Ensure overlapping gene set and identical column order
  g <- intersect(colnames(matA), colnames(matB))
  if (length(g) < 2) return(NA_real_)
  g <- sort(g)
  
  matA <- matA[, g, drop = FALSE]
  matB <- matB[, g, drop = FALSE]
  
  c1 <- stats::cor(matA, use = "pairwise.complete.obs")
  c2 <- stats::cor(matB, use = "pairwise.complete.obs")
  
  c1[is.na(c1)] <- 0
  c2[is.na(c2)] <- 0
  
  # final sanity
  if (!identical(dim(c1), dim(c2))) return(NA_real_)
  
  num <- sum(c1 * c2)
  den <- sqrt(sum(c1^2) * sum(c2^2))
  
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

maybe_run_bootstrap <- FALSE
if (metric == "RV") {
  if (!is.na(run_boot_arg) && nchar(run_boot_arg) > 0) {
    maybe_run_bootstrap <- is_true(run_boot_arg)
  } else {
    maybe_run_bootstrap <- ask_yesno("Run optimized RV bootstrap (row-resampling)?", default_yes = TRUE)
  }
}

if (metric == "RV" && maybe_run_bootstrap) {
  
  if (is.na(expr_path) || nchar(expr_path) == 0) {
    expr_path <- pick_file("Choose tidy expression CSV for bootstrap (Day, SampleID, Gene, Abundance)")
  }
  if (!file.exists(expr_path)) stop2("Expression file not found: ", expr_path)
  
  nboot <- NA_integer_
  if (is.na(nboot_arg) || nchar(nboot_arg) == 0) {
    nboot <- ask_int("How many bootstrap replicates?", default = 200)
  } else {
    nboot <- suppressWarnings(as.integer(nboot_arg))
    if (is.na(nboot) || nboot < 1) stop2("Invalid --nboot (must be positive integer).")
  }
  
  msg("Bootstrap settings: expr=", expr_path,
      " | nboot=", nboot,
      " | seed=", seed,
      " | ci=", ci_level,
      " | transform=", transform,
      " | max_genes=", max_genes)
  
  # Read expression
  expr_data <- read.csv(expr_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  req_expr <- c("Day", "SampleID", "Gene", "Abundance")
  missing_expr <- setdiff(req_expr, colnames(expr_data))
  if (length(missing_expr) > 0) stop2("Expr file missing required columns: ", paste(missing_expr, collapse = ", "))
  
  expr_data$Day <- suppressWarnings(as.numeric(expr_data$Day))
  expr_data$Abundance <- suppressWarnings(as.numeric(expr_data$Abundance))
  if (any(is.na(expr_data$Day))) stop2("Expr Day has non-numeric values.")
  if (any(is.na(expr_data$Abundance))) stop2("Expr Abundance has non-numeric values (or NA).")
  
  if (!(ref_day %in% unique(expr_data$Day))) stop2("Reference day ", ref_day, " not found in expression data.")
  
  msg("Pre-computing per-day matrices (this can take a while for large gene sets)...")
  master_mats <- vector("list", length(days))
  names(master_mats) <- as.character(days)
  
  debug_day_summary <- data.frame(
    Day = days,
    n_samples = NA_integer_,
    n_genes_pre = NA_integer_,
    n_genes_post = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (idx in seq_along(days)) {
    d <- days[idx]
    
    # quick pre-counts from tidy input (debug)
    if (debug) {
      dd <- expr_data[expr_data$Day == d, , drop = FALSE]
      debug_day_summary$n_samples[idx]   <- length(unique(dd$SampleID))
      debug_day_summary$n_genes_pre[idx] <- length(unique(dd$Gene))
    }
    
    master_mats[[as.character(d)]] <- build_day_matrix(expr_data, d, trans = transform, max_g = max_genes)
    
    if (is.null(master_mats[[as.character(d)]])) {
      stop2("Day ", d, " produced NULL matrix after preprocessing (too few genes/samples). ",
            "Try --max_genes=1000 for troubleshooting, and verify each day has >=2 samples.")
    }
    
    if (debug) {
      msg("Day ", d, ": ", nrow(master_mats[[as.character(d)]]), " samples × ",
          ncol(master_mats[[as.character(d)]]), " genes AFTER preprocessing")
      debug_day_summary$n_genes_post[idx] <- ncol(master_mats[[as.character(d)]])
    }
    
    if (d %% 2 == 0) gc()
  }
  
  if (debug && save_debug_files) {
    write.csv(debug_day_summary, file.path(out_dir, "DEBUG_day_matrix_summary.csv"), row.names = FALSE)
    msg("Wrote debug day summary: ", file.path(out_dir, "DEBUG_day_matrix_summary.csv"))
  }
  
  # ------------------------------------------------------------
  # Harmonize gene set across days (critical for RV conformance)
  # ------------------------------------------------------------
  msg("Harmonizing gene set across days (intersection of columns)...")
  gene_sets <- lapply(master_mats, colnames)
  common_genes <- Reduce(intersect, gene_sets)
  
  if (length(common_genes) < 2) {
    stop2(
      "After preprocessing, fewer than 2 genes are shared across all days (",
      length(common_genes), "). RV cannot be computed.\n",
      "Likely causes: aggressive filtering from missingness/zero-variance; days differ in gene coverage.\n",
      "Try: (i) ensure consistent gene coverage across days; (ii) set --max_genes to a smaller consistent subset; ",
      "(iii) relax preprocessing filters."
    )
  }
  
  common_genes <- sort(common_genes)
  for (nm in names(master_mats)) {
    master_mats[[nm]] <- master_mats[[nm]][, common_genes, drop = FALSE]
  }
  msg("Common genes retained across all days: ", length(common_genes))
  
  if (debug && save_debug_files) {
    # Save common gene list
    write.csv(data.frame(Gene = common_genes), file.path(out_dir, "DEBUG_common_genes_across_days.csv"), row.names = FALSE)
    msg("Wrote common gene list: ", file.path(out_dir, "DEBUG_common_genes_across_days.csv"))
  }
  
  ref_mat <- master_mats[[as.character(ref_day)]]
  if (is.null(ref_mat)) stop2("Reference day matrix is NULL after preprocessing.")
  
  # Observed RV profile
  msg("Computing observed RV-to-Day", ref_day, " profile...")
  rv_obs <- sapply(days, function(d) calc_rv(master_mats[[as.character(d)]], ref_mat))
  
  if (debug) {
    bad <- which(!is.finite(rv_obs))
    if (length(bad) > 0) {
      msg("DEBUG: Non-finite RV observed for days: ", paste(days[bad], collapse = ", "))
    }
  }
  
  # Bootstrap draws (row-resampling)
  set.seed(seed)
  boot_results <- matrix(NA_real_, nrow = nboot, ncol = length(days),
                         dimnames = list(NULL, as.character(days)))
  
  msg("Starting bootstrap loop...")
  for (b in seq_len(nboot)) {
    
    ref_b <- ref_mat[sample.int(nrow(ref_mat), size = nrow(ref_mat), replace = TRUE), , drop = FALSE]
    
    for (d in as.character(days)) {
      target_mat <- master_mats[[d]]
      target_b <- target_mat[sample.int(nrow(target_mat), size = nrow(target_mat), replace = TRUE), , drop = FALSE]
      boot_results[b, d] <- calc_rv(target_b, ref_b)
    }
    
    if (debug && b == 1 && save_debug_files) {
      # Save first draw as a quick check
      write.csv(data.frame(Day = days, RV_b1 = as.numeric(boot_results[b, ])),
                file.path(out_dir, "DEBUG_bootstrap_draw1_profile.csv"),
                row.names = FALSE)
      msg("Wrote debug first bootstrap draw: ", file.path(out_dir, "DEBUG_bootstrap_draw1_profile.csv"))
    }
    
    if (b %% max(1, floor(nboot / 10)) == 0) {
      msg("  bootstrap ", b, "/", nboot)
      gc()
    }
  }
  
  ci_low  <- apply(boot_results, 2, quantile, probs = alpha,     na.rm = TRUE)
  ci_high <- apply(boot_results, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
  
  out_df <- data.frame(
    Day = days,
    RV_observed = as.numeric(rv_obs),
    CI_low  = as.numeric(ci_low[as.character(days)]),
    CI_high = as.numeric(ci_high[as.character(days)]),
    RefDay = ref_day,
    nboot = nboot,
    ci_level = ci_level,
    seed = seed,
    transform = transform,
    max_genes = max_genes,
    n_common_genes = length(common_genes)
  )
  
  write.csv(out_df, file.path(out_dir, "RefDay_Profile_RV_bootstrap_CI.csv"), row.names = FALSE)
  msg("Bootstrap complete. Wrote: ", file.path(out_dir, "RefDay_Profile_RV_bootstrap_CI.csv"))
  
  if (save_boot_draws) {
    write.csv(boot_results, file.path(out_dir, "RV_bootstrap_draws_profile_wide.csv"), row.names = FALSE)
    msg("Saved bootstrap draws (wide): ", file.path(out_dir, "RV_bootstrap_draws_profile_wide.csv"))
  }
  
} else {
  if (metric == "RV") {
    msg("RV bootstrap not run (skipped or disabled).")
  }
}

# ------------------------------ Finish ----------------------------------

msg("Done. Files written to: ", out_dir)
msg("Key outputs:")
msg(" - similarity_long_clean.csv")
msg(" - similarity_wide_*.csv")
msg(" - RefDay_Profile_", metric, "_data.csv")
msg(" - Plot_RefDay_Profile_", metric, ".png/.pdf")
msg(" - Plot_MDS_", metric, ".png/.pdf")
msg(" - Plot_DayNetwork_", metric, ".png/.pdf")
msg(" - mds_coords_", metric, ".csv")
msg(" - network_edges_", metric, "_topk.csv")
if (metric == "RV") {
  msg(" - RefDay_Profile_RV_bootstrap_CI.csv (if bootstrap ran)")
  msg(" - DEBUG_day_matrix_summary.csv, DEBUG_common_genes_across_days.csv (if --debug=TRUE)")
}
