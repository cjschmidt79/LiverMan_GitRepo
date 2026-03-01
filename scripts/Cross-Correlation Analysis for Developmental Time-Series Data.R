#!/usr/bin/env Rscript
# ======================================================================
# Δ–Δ Cross-Correlation Analysis for Developmental Time-Series Data
#
# PURPOSE
#   This script quantifies coordinated temporal changes between genes across
#   developmental time by computing cross-correlations on first differences
#   (Δ–Δ correlations).
#
#   The analysis:
#     1) Optionally transforms expression (raw OR log2(fpkm + pseudocount)).
#     2) Aggregates biological replicates within each time point
#        (per-day medians or means).
#     3) Computes first differences across successive time points.
#     4) Calculates cross-correlations between Δ-series for a target gene
#        and query genes across user-defined lags.
#     5) Reports Pearson and Spearman correlations, permutation-based p-values,
#        and BH-FDR-adjusted q-values.
#     6) Generates a single Quarto HTML report in the same output_dir as results,
#        using a WD-safe rendering pattern.
#
# EXPECTED INPUT FORMAT (LONG TABLE)
#   A CSV/TSV with at least:
#     Day      : numeric/integer developmental time point
#     SampleID : biological replicate identifier (e.g., "4_1" = bird 1 at day 4)
#     Source   : gene identifier (e.g., gene symbol)
#     fpkm     : numeric expression value
#
# IMPORTANT ASSUMPTIONS
#   • Samples within a Day are biological replicates.
#   • Replicates are NOT tracked longitudinally across days.
#   • First differences are computed on per-day summary statistics (not per-sample).
#
# OUTPUTS
#   • DailySummary_ByDayGene.csv
#   • DeltaDelta_XCorr_AllTests.csv  (Pearson+Spearman, perm p, BH q)
#   • DeltaDelta_XCorr_Lag0.csv      (lag=0 subset)
#   • DeltaDelta_XCorr_BestByAbsR.csv
#   • Figures/*.png                 (minimal QC plots)
#   • <script_name>_Report.qmd + .html in output_dir
#
# NO dplyr
# ======================================================================

suppressWarnings(suppressMessages({
  if (!requireNamespace("future", quietly = TRUE)) stop("Package 'future' is required.")
  if (!requireNamespace("furrr", quietly = TRUE)) stop("Package 'furrr' is required.")
  if (!requireNamespace("quarto", quietly = TRUE)) stop("Package 'quarto' is required for reporting.")
  if (!requireNamespace("knitr", quietly = TRUE)) stop("Package 'knitr' is required for report embedding.")
}))

# ---------------------------
# Utilities (no dplyr)
# ---------------------------
timestamp_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

normalize_slash <- function(x, mustWork = FALSE) {
  normalizePath(x, winslash = "/", mustWork = mustWork)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

rstudio_select_file <- function(caption = "Select input file") {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(rstudioapi::selectFile(caption = caption))
  }
  return(NULL)
}

rstudio_select_dir <- function(caption = "Select output directory") {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(rstudioapi::selectDirectory(caption = caption))
  }
  return(NULL)
}

get_script_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- grep("^--file=", cmdArgs, value = TRUE)
  if (length(fileArg) == 1) return(sub("^--file=", "", fileArg))
  if (!is.null(sys.frames()[[1]]$ofile)) return(sys.frames()[[1]]$ofile)
  return(NA_character_)
}

read_table_auto <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  if (length(first) == 0) stop("Input file is empty: ", path)
  sep <- if (grepl("\t", first)) "\t" else ","
  df <- tryCatch(
    read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read input file: ", e$message)
  )
  df
}

safe_cor_stats <- function(x, y, method = c("pearson","spearman")) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3) return(list(r = NA_real_, n = n, t = NA_real_, p_param = NA_real_))
  r <- suppressWarnings(cor(x, y, method = method))
  if (!is.finite(r)) r <- NA_real_
  
  if (is.na(r) || abs(r) >= 1) {
    tval <- if (!is.na(r) && abs(r) == 1) Inf else NA_real_
    pval <- if (!is.na(r) && abs(r) == 1) 0 else NA_real_
  } else {
    tval <- r * sqrt((n - 2) / (1 - r^2))
    pval <- 2 * pt(-abs(tval), df = n - 2)
  }
  list(r = r, n = n, t = tval, p_param = pval)
}

perm_p_cor <- function(x, y, method = c("pearson","spearman"), B = 2000L, seed = 1L) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4) return(NA_real_)
  r_obs <- suppressWarnings(cor(x, y, method = method))
  if (!is.finite(r_obs)) return(NA_real_)
  set.seed(seed)
  r_null <- numeric(B)
  for (b in seq_len(B)) {
    r_null[b] <- suppressWarnings(cor(x, sample(y, size = n, replace = FALSE), method = method))
  }
  (sum(abs(r_null) >= abs(r_obs)) + 1) / (B + 1)
}

align_lag <- function(dx, dy, day_d, k) {
  if (length(dx) != length(dy)) stop("dx and dy lengths differ unexpectedly.")
  L <- length(dx)
  if (k >= 0) {
    if (L - k <= 1) return(list(x = numeric(0), y = numeric(0), days = numeric(0)))
    x <- dx[seq_len(L - k)]
    y <- dy[(1 + k):L]
    days <- day_d[(1 + k):length(day_d)]
  } else {
    kk <- abs(k)
    if (L - kk <= 1) return(list(x = numeric(0), y = numeric(0), days = numeric(0)))
    x <- dx[(1 + kk):L]
    y <- dy[seq_len(L - kk)]
    days <- day_d[seq_len(length(day_d) - kk)]
  }
  list(x = x, y = y, days = days)
}

# ---------------------------
# NEW: Parallel worker suggestion + optional auto-benchmark
# ---------------------------
suggest_workers <- function(n_tasks,
                            perm_B,
                            avail = future::availableCores()) {
  avail <- max(1L, as.integer(avail))
  n_tasks <- max(1L, as.integer(n_tasks))
  
  frac <- if (perm_B >= 10000L) 0.90 else if (perm_B >= 5000L) 0.80 else 0.65
  if (interactive()) frac <- min(frac, 0.70)
  
  w <- as.integer(floor(avail * frac))
  w <- max(1L, w)
  w <- min(w, n_tasks)
  
  if (perm_B < 2000L) w <- min(w, 8L)
  
  w
}

benchmark_workers <- function(gene_subset,
                              compute_gene_tests_fun,
                              candidates,
                              seed = 1L) {
  candidates <- unique(as.integer(candidates))
  candidates <- candidates[candidates >= 1L]
  out <- data.frame(workers = candidates, seconds = NA_real_, stringsAsFactors = FALSE)
  
  for (i in seq_along(candidates)) {
    w <- candidates[i]
    future::plan(future::multisession, workers = w)
    t0 <- proc.time()[["elapsed"]]
    invisible(furrr::future_map(gene_subset, compute_gene_tests_fun,
                                .options = furrr::furrr_options(seed = TRUE)))
    out$seconds[i] <- proc.time()[["elapsed"]] - t0
  }
  
  out[order(out$seconds), , drop = FALSE]
}

# ---------------------------
# Interactive inputs
# ---------------------------
cat("\n=== Δ–Δ Cross-Correlation (Developmental Time-Series) ===\n")

in_path <- rstudio_select_file("Select input CSV/TSV (Day, SampleID, Source, fpkm)")
if (is.null(in_path) || is.na(in_path) || !nzchar(in_path)) {
  cat("No RStudio file picker selection detected.\n")
  in_path <- file.choose()
}
in_path <- normalize_slash(in_path, mustWork = TRUE)

cat("\nInput file:\n  ", in_path, "\n", sep = "")

run_name <- readline("\nEnter run name (e.g., DeltaDelta_GRHPR) [default: DeltaDelta]: ")
if (!nzchar(run_name)) run_name <- "DeltaDelta"

base_out <- rstudio_select_dir("Select base output directory (or Cancel to use ./outputs)")
if (is.null(base_out) || is.na(base_out) || !nzchar(base_out)) {
  base_out <- file.path(getwd(), "outputs")
}
safe_dir_create(base_out)

ts <- timestamp_now()
output_dir <- file.path(base_out, paste0(run_name, "_", ts))
safe_dir_create(output_dir)
output_dir <- normalize_slash(output_dir, mustWork = TRUE)

cat("\nOutput directory:\n  ", output_dir, "\n", sep = "")

# Transform selection as requested:
# 1 = raw (default)
# 2 = log2 + pseudocount
cat("\nExpression transform:\n  1 = raw (default)\n  2 = log2(fpkm + pseudocount)\n")
transform_choice <- readline("Choose 1 or 2 [default: 1]: ")
transform_choice <- trimws(transform_choice)
if (!nzchar(transform_choice)) transform_choice <- "1"
transform_choice <- suppressWarnings(as.integer(transform_choice))
if (!is.finite(transform_choice) || is.na(transform_choice) || !transform_choice %in% c(1L, 2L)) transform_choice <- 1L

pseudocount <- 1
transform_label <- "raw"
if (transform_choice == 2L) {
  pc_in <- readline("Pseudocount for log2(fpkm + pseudocount) [default: 1]: ")
  pc_in <- suppressWarnings(as.numeric(trimws(pc_in)))
  if (is.finite(pc_in) && !is.na(pc_in) && pc_in > 0) pseudocount <- pc_in else pseudocount <- 1
  transform_label <- paste0("log2p(pc=", pseudocount, ")")
}

daily_stat <- readline("\nDaily summary statistic (median | mean) [default: median]: ")
daily_stat <- tolower(trimws(daily_stat))
if (!daily_stat %in% c("median","mean")) daily_stat <- "median"

target_gene <- readline("\nTarget gene (Source) [default: GRHPR]: ")
if (!nzchar(target_gene)) target_gene <- "GRHPR"

max_abs_lag <- readline("\nMaximum |lag| in steps (e.g., 1 gives -1:1) [default: 1]: ")
max_abs_lag <- suppressWarnings(as.integer(trimws(max_abs_lag)))
if (!is.finite(max_abs_lag) || is.na(max_abs_lag) || max_abs_lag < 0) max_abs_lag <- 1L
lags_steps <- seq.int(-max_abs_lag, max_abs_lag, by = 1L)

sampling_interval_days <- readline("\nSampling interval in days (for lag_days reporting) [default: auto-infer]: ")
sampling_interval_days <- suppressWarnings(as.numeric(trimws(sampling_interval_days)))
if (!is.finite(sampling_interval_days) || is.na(sampling_interval_days) || sampling_interval_days <= 0) {
  sampling_interval_days <- NA_real_
}

# (Already present in your script; kept intact)
perm_B <- readline("\nPermutation iterations for p-values [default: 2000]: ")
perm_B <- suppressWarnings(as.integer(trimws(perm_B)))
if (!is.finite(perm_B) || is.na(perm_B) || perm_B < 200) perm_B <- 2000L

set_seed <- readline("\nRandom seed [default: 1]: ")
set_seed <- suppressWarnings(as.integer(trimws(set_seed)))
if (!is.finite(set_seed) || is.na(set_seed)) set_seed <- 1L

# Query genes selection
cat("\nQuery gene selection:\n",
    "  1) Analyze ALL genes except target (can be slow)\n",
    "  2) Provide a comma-separated list of genes\n", sep = "")
mode_sel <- readline("Choose 1 or 2 [default: 1]: ")
mode_sel <- trimws(mode_sel); if (!nzchar(mode_sel)) mode_sel <- "1"

query_genes_user <- NULL
if (mode_sel == "2") {
  gl <- readline("Enter genes (comma-separated): ")
  gl <- trimws(gl)
  if (nzchar(gl)) {
    query_genes_user <- unique(trimws(strsplit(gl, ",", fixed = TRUE)[[1]]))
    query_genes_user <- query_genes_user[nzchar(query_genes_user)]
  }
}

# ---------------------------
# Read + validate input
# ---------------------------
df <- read_table_auto(in_path)

required_cols <- c("Day","SampleID","Source","fpkm")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Basic type coercions
df$Day <- suppressWarnings(as.numeric(df$Day))
df$fpkm <- suppressWarnings(as.numeric(df$fpkm))
df$SampleID <- as.character(df$SampleID)
df$Source <- as.character(df$Source)

if (any(!is.finite(df$Day))) stop("Non-numeric Day values detected after coercion.")
if (any(!is.finite(df$fpkm))) stop("Non-numeric fpkm values detected after coercion.")

# Unique days and inferred interval if not provided
days <- sort(unique(df$Day))
if (length(days) < 4) stop("Need at least 4 distinct time points (days) to compute Δ–Δ reliably.")

if (is.na(sampling_interval_days)) {
  diffs <- diff(days)
  # Use median spacing as "reporting" interval; analysis itself uses steps
  sampling_interval_days <- stats::median(diffs)
}

# ---------------------------
# Optional expression transform (applied BEFORE aggregation)
# ---------------------------
if (transform_choice == 2L) {
  if (any(df$fpkm < 0, na.rm = TRUE)) stop("Negative fpkm values detected; cannot log-transform.")
  df$fpkm_transformed <- log2(df$fpkm + pseudocount)
} else {
  df$fpkm_transformed <- df$fpkm
}

# ---------------------------
# Aggregate within Day x Source across biological replicates
# ---------------------------
agg_fun <- if (daily_stat == "median") stats::median else base::mean

# aggregate() returns data.frame with columns: Day, Source, fpkm_transformed
daily <- aggregate(df$fpkm_transformed,
                   by = list(Day = df$Day, Source = df$Source),
                   FUN = agg_fun, na.rm = TRUE)
names(daily)[names(daily) == "x"] <- "daily_value"

# Save daily summary (long)
daily_out <- file.path(output_dir, "DailySummary_ByDayGene.csv")
write.csv(daily, daily_out, row.names = FALSE)

# Wide matrix: rows = Day, cols = Source
daily_wide <- reshape(daily,
                      idvar = "Day",
                      timevar = "Source",
                      direction = "wide")
# Clean column names: daily_value.GENE -> GENE
colnames(daily_wide) <- sub("^daily_value\\.", "", colnames(daily_wide))
daily_wide <- daily_wide[order(daily_wide$Day), , drop = FALSE]

if (!target_gene %in% colnames(daily_wide)) {
  stop("Target gene not found in daily summary: ", target_gene)
}

# Determine query genes
all_genes <- setdiff(colnames(daily_wide), "Day")
all_genes <- unique(all_genes)
query_genes <- if (is.null(query_genes_user)) setdiff(all_genes, target_gene) else setdiff(intersect(all_genes, query_genes_user), target_gene)

if (length(query_genes) == 0) stop("No query genes available after filtering.")

cat("\nGenes in analysis:\n")
cat("  Target:", target_gene, "\n")
cat("  #Query:", length(query_genes), "\n")

# ---------------------------
# Build Δ-series for target
# ---------------------------
Y <- daily_wide[[target_gene]]
day_vec <- daily_wide$Day
dY <- diff(Y)
dY_day <- day_vec[-1]  # Δ at t corresponds to day t (excluding first)

# ---------------------------
# Per-gene computation (function defined BEFORE parallel setup)
# ---------------------------
compute_gene_tests <- function(gene_name) {
  X <- daily_wide[[gene_name]]
  dX <- diff(X)
  
  rows <- vector("list", length(lags_steps))
  for (i in seq_along(lags_steps)) {
    k <- lags_steps[i]
    al <- align_lag(dX, dY, dY_day, k)
    x_al <- al$x; y_al <- al$y
    
    pear <- safe_cor_stats(x_al, y_al, method = "pearson")
    spear <- safe_cor_stats(x_al, y_al, method = "spearman")
    
    # Permutation p-values (primary)
    p_perm_pear <- if (pear$n >= 4 && is.finite(pear$r)) perm_p_cor(x_al, y_al, method = "pearson", B = perm_B, seed = set_seed) else NA_real_
    p_perm_spear <- if (spear$n >= 4 && is.finite(spear$r)) perm_p_cor(x_al, y_al, method = "spearman", B = perm_B, seed = set_seed) else NA_real_
    
    rows[[i]] <- data.frame(
      gene = gene_name,
      target = target_gene,
      lag_steps = k,
      lag_days = k * sampling_interval_days,
      n_pairs = pear$n,
      r_pearson = pear$r,
      t_pearson = pear$t,
      p_param_pearson = pear$p_param,
      p_perm_pearson = p_perm_pear,
      r_spearman = spear$r,
      t_spearman = spear$t,
      p_param_spearman = spear$p_param,
      p_perm_spearman = p_perm_spear,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---------------------------
# NEW: Parallel worker selection (suggested / manual / auto-benchmark)
# ---------------------------
avail <- future::availableCores()
workers_suggested <- suggest_workers(n_tasks = length(query_genes), perm_B = perm_B, avail = avail)

cat("\nParallel configuration:\n")
cat("  Available cores: ", avail, "\n", sep = "")
cat("  Suggested workers: ", workers_suggested, "\n", sep = "")
cat("  Worker selection modes:\n")
cat("    1) Use suggested workers (default)\n")
cat("    2) Manual workers\n")
cat("    3) Auto-benchmark (pilot subset)\n")

w_mode <- readline("Choose 1, 2, or 3 [default: 1]: ")
w_mode <- trimws(w_mode)
if (!nzchar(w_mode)) w_mode <- "1"
w_mode <- suppressWarnings(as.integer(w_mode))
if (!is.finite(w_mode) || is.na(w_mode) || !w_mode %in% c(1L, 2L, 3L)) w_mode <- 1L

workers_in <- workers_suggested

if (w_mode == 2L) {
  w_in <- readline(paste0("Number of workers [default: ", workers_suggested, "]: "))
  w_in <- suppressWarnings(as.integer(trimws(w_in)))
  if (is.finite(w_in) && !is.na(w_in) && w_in >= 1L) {
    workers_in <- min(as.integer(w_in), avail, length(query_genes))
  } else {
    workers_in <- workers_suggested
  }
} else if (w_mode == 3L) {
  set.seed(set_seed)
  pilot_n <- min(50L, length(query_genes))
  pilot_genes <- sample(query_genes, pilot_n)
  
  # Candidate grid (bounded by avail and task count)
  cands <- unique(c(1L, 2L, 4L, 6L, 8L, 10L, 12L, 14L, 16L))
  cands <- cands[cands <= avail]
  cands <- cands[cands <= length(query_genes)]
  if (length(cands) == 0) cands <- 1L
  
  cat("\nAuto-benchmarking workers on pilot set (n=", pilot_n, ")...\n", sep = "")
  bench <- benchmark_workers(pilot_genes, compute_gene_tests, cands, seed = set_seed)
  print(bench)
  
  workers_in <- bench$workers[1]
  cat("Auto-selected workers: ", workers_in, "\n", sep = "")
}

cat("Using workers: ", workers_in, "\n", sep = "")

# ---------------------------
# Parallel setup (unchanged behavior except workers chosen above)
# ---------------------------
future::plan(future::multisession, workers = workers_in)
furrr::furrr_options(seed = TRUE)

# ---------------------------
# Per-gene computation (parallel)
# ---------------------------
cat("\nRunning Δ–Δ cross-correlations in parallel...\n")
res_list <- furrr::future_map(query_genes, compute_gene_tests, .options = furrr::furrr_options(seed = TRUE))
res <- do.call(rbind, res_list)

# BH-FDR across all tests within each method (family = all gene x lag for this run)
res$q_bh_pearson  <- p.adjust(res$p_perm_pearson,  method = "BH")
res$q_bh_spearman <- p.adjust(res$p_perm_spearman, method = "BH")

# Save full table
out_all <- file.path(output_dir, "DeltaDelta_XCorr_AllTests.csv")
write.csv(res, out_all, row.names = FALSE)

# Lag 0 subset
res_lag0 <- res[res$lag_steps == 0, , drop = FALSE]
out_lag0 <- file.path(output_dir, "DeltaDelta_XCorr_Lag0.csv")
write.csv(res_lag0, out_lag0, row.names = FALSE)

# Best lag per gene by |Pearson r| (tie-break by smaller |lag|, then larger n_pairs)
best_by_absr <- do.call(rbind, lapply(split(res, res$gene), function(d) {
  d <- d[order(-abs(d$r_pearson), abs(d$lag_steps), -d$n_pairs), , drop = FALSE]
  d[1, , drop = FALSE]
}))
out_best <- file.path(output_dir, "DeltaDelta_XCorr_BestByAbsR.csv")
write.csv(best_by_absr, out_best, row.names = FALSE)

# ---------------------------
# Minimal figures (base R) to help QC
# ---------------------------
fig_dir <- file.path(output_dir, "Figures")
safe_dir_create(fig_dir)

png(file.path(fig_dir, paste0("Target_", target_gene, "_Daily_", daily_stat, "_", transform_label, ".png")), width = 1600, height = 900, res = 200)
plot(day_vec, Y, type = "b", xlab = "Day", ylab = paste0(target_gene, " (daily ", daily_stat, ", ", transform_label, ")"),
     main = paste0("Target trajectory: ", target_gene))
dev.off()

png(file.path(fig_dir, "Lag0_PearsonR_Hist.png"), width = 1600, height = 900, res = 200)
hist(res_lag0$r_pearson, breaks = 40, main = "Lag 0 Δ–Δ Pearson r (all query genes)", xlab = "r")
dev.off()

png(file.path(fig_dir, "Lag0_SpearmanR_Hist.png"), width = 1600, height = 900, res = 200)
hist(res_lag0$r_spearman, breaks = 40, main = "Lag 0 Δ–Δ Spearman rho (all query genes)", xlab = "rho")
dev.off()

# ---------------------------
# Quarto report generation (WD-safe canonical pattern)
# ---------------------------
script_path <- get_script_path()
script_path_norm <- if (!is.na(script_path) && nzchar(script_path)) normalize_slash(script_path, mustWork = FALSE) else NA_character_
script_name <- if (!is.na(script_path_norm)) basename(script_path_norm) else "DeltaDelta_XCorr.R"
script_full <- tryCatch(if (!is.na(script_path_norm) && file.exists(script_path_norm)) paste(readLines(script_path_norm, warn = FALSE), collapse = "\n") else NA_character_,
                        error = function(e) NA_character_)

# Extract header comment block (first ~200 lines until first non-comment, loosely)
header_block <- NA_character_
if (!is.na(script_path_norm) && file.exists(script_path_norm)) {
  L <- readLines(script_path_norm, warn = FALSE)
  Lh <- L[seq_len(min(250, length(L)))]
  # keep initial contiguous comment lines + shebang
  keep <- logical(length(Lh))
  for (i in seq_along(Lh)) {
    line <- Lh[i]
    if (i == 1 && grepl("^#!", line)) { keep[i] <- TRUE; next }
    if (grepl("^\\s*#", line) || grepl("^\\s*$", line)) keep[i] <- TRUE else break
  }
  header_block <- paste(Lh[keep], collapse = "\n")
}

# Build QMD in a single character vector (no fragments)
qmd_title <- paste0("Δ–Δ Cross-Correlation Report: ", run_name)
qmd_file <- file.path(output_dir, paste0(sub("\\.R$", "", script_name), "_Report.qmd"))
qmd_base <- basename(qmd_file)

report_lines <- c(
  "---",
  paste0('title: "', gsub('"', '\\"', qmd_title), '"'),
  "format:",
  "  html:",
  "    toc: true",
  "    code-fold: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "## Summary",
  "",
  "This report summarizes a Δ–Δ (first-difference) cross-correlation analysis computed on per-day aggregated expression across biological replicates.",
  "",
  "## Script header",
  "",
  "```",
  if (!is.na(header_block)) header_block else "(Header block unavailable.)",
  "```",
  "",
  "## Metadata",
  "",
  "```",
  paste0("Timestamp: ", ts),
  paste0("Script name: ", script_name),
  paste0("Script path: ", if (!is.na(script_path_norm)) script_path_norm else "NA"),
  paste0("Input path: ", in_path),
  paste0("Output dir: ", output_dir),
  paste0("Target gene: ", target_gene),
  paste0("Expression transform: ", transform_label),
  paste0("Daily summary: ", daily_stat),
  paste0("Max |lag| steps: ", max_abs_lag),
  paste0("Sampling interval days (for lag_days): ", sampling_interval_days),
  paste0("Permutation B: ", perm_B),
  paste0("Seed: ", set_seed),
  paste0("Workers: ", workers_in),
  paste0("#Query genes: ", length(query_genes)),
  "```",
  "",
  "## Dependencies",
  "",
  "```{r}",
  "pkgs <- c('future','furrr','quarto','knitr')",
  "ip <- installed.packages()[, c('Package','Version')]",
  "ip <- ip[ip[, 'Package'] %in% pkgs, , drop = FALSE]",
  "print(ip)",
  "```",
  "",
  "## Analytical logic",
  "",
  "Let $X_t$ be the per-day summary (median or mean) for a query gene and $Y_t$ for the target gene (after optional transformation). First differences are computed across successive sampled days:",
  "",
  "$$\\Delta X_t = X_t - X_{t-1}, \\quad \\Delta Y_t = Y_t - Y_{t-1}.$$",
  "",
  "For lag $k$ (in sampling steps), the Δ–Δ cross-correlation is computed as:",
  "",
  "$$\\mathrm{cor}(\\Delta X_t, \\Delta Y_{t+k}).$$",
  "",
  "Both Pearson and Spearman correlations are reported, with permutation p-values and BH-FDR-adjusted q-values (computed across all gene × lag tests separately for each correlation type).",
  "",
  "## Generated outputs",
  "",
  "```{r}",
  "out_dir <- normalizePath('.', winslash='/', mustWork=TRUE)",
  "files <- list.files(out_dir, recursive=TRUE)",
  "cat(paste(files, collapse='\\n'))",
  "```",
  "",
  "## Key tables",
  "",
  paste0("- `", basename(daily_out), "`: per-day aggregated expression (long)"),
  paste0("- `", basename(out_all), "`: all Δ–Δ tests (lags × genes)"),
  paste0("- `", basename(out_lag0), "`: lag 0 subset"),
  paste0("- `", basename(out_best), "`: best lag by |Pearson r| per gene"),
  "",
  "## Figures",
  "",
  "```{r}",
  "fig_dir <- file.path(normalizePath('.', winslash='/', mustWork=TRUE), 'Figures')",
  "imgs <- list.files(fig_dir, pattern='\\\\.(png|pdf)$', full.names=TRUE, ignore.case=TRUE)",
  "if (length(imgs) == 0) {",
  "  cat('No figures found in Figures/.')",
  "} else {",
  "  for (p in imgs) {",
  "    cat('### ', basename(p), '\\n')",
  "    knitr::include_graphics(p)",
  "    cat('\\n\\n')",
  "  }",
  "}",
  "```",
  "",
  "## Session info",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

# WD-safe render pattern
old_wd <- getwd()
setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)

writeLines(report_lines, qmd_base)

# Render HTML (no output_dir argument)
quarto::quarto_render(qmd_base)

html_path <- sub("\\.qmd$", ".html", qmd_file)

# ---------------------------
# Final console output
# ---------------------------
cat("\n\n=== DONE ===\n")
cat("Daily summary (long):\n  ", daily_out, "\n", sep = "")
cat("All tests:\n  ", out_all, "\n", sep = "")
cat("Lag 0 subset:\n  ", out_lag0, "\n", sep = "")
cat("Best-by-|r|:\n  ", out_best, "\n", sep = "")
cat("Figures dir:\n  ", normalize_slash(fig_dir, mustWork = TRUE), "\n", sep = "")
cat("QMD report:\n  ", qmd_file, "\n", sep = "")
cat("HTML report:\n  ", html_path, "\n", sep = "")

if (interactive()) {
  cat("\nOpening HTML report in browser...\n")
  utils::browseURL(html_path)
}
