#!/usr/bin/env Rscript
# ======================================================================
# Δ–Δ Cross-Correlation Analysis for Developmental Time-Series Data
#
# PURPOSE
#   Compute Δ–Δ cross-correlations on per-day aggregated expression.
#   Supports:
#     • Single-target vs all genes
#     • All-targets (all-by-all) via Top-K shortlisting per target
#     • Pearson + Spearman correlations
#     • Permutation p-values + BH FDR (run on shortlisted pairs in all-targets mode)
#     • Optional transform: 1=raw, 2=log2(fpkm + pseudocount) (default pseudocount=1)
#     • future/furrr parallel processing with suggested workers + optional benchmark
#     • Single Quarto HTML report in output_dir using WD-safe rendering
#
# EXPECTED INPUT FORMAT (LONG TABLE)
#   Required columns:
#     Day      : numeric/integer time point
#     SampleID : replicate id (e.g., "4_1" bird 1 day 4) (replicates within day)
#     Source   : gene id (e.g., symbol)
#     fpkm     : numeric expression value
#
# NOTES
#   Replicates are within-day biological replicates, not longitudinal tracks.
#   Δ is computed on per-day aggregate (median or mean).
#
# OUTPUTS (key)
#   • DailySummary_ByDayGene.csv
#   • Shortlist_TopK_ByTarget.csv                (all-targets mode)
#   • DeltaDelta_XCorr_AllTests.csv              (full inference table)
#   • DeltaDelta_XCorr_Lag0.csv                  (lag=0 subset)
#   • DeltaDelta_XCorr_BestByAbsR.csv            (best lag per (target,gene) by |Pearson r|)
#   • Figures/*.png
#   • <script_name>_Report.qmd + .html
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

# Deterministic per-test seed so permutations don't reuse identical RNG streams
seed_from_key <- function(base_seed, target, gene, lag, method) {
  key <- paste(base_seed, target, gene, lag, method, sep = "|")
  v <- utf8ToInt(key)
  s <- base_seed + sum(v * seq_along(v)) %% 1000000000L
  as.integer(s)
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
# Parallel worker suggestion + optional benchmark
# ---------------------------
suggest_workers <- function(n_tasks, perm_B, avail = future::availableCores()) {
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

benchmark_workers <- function(task_subset, task_fun, candidates) {
  candidates <- unique(as.integer(candidates))
  candidates <- candidates[candidates >= 1L]
  out <- data.frame(workers = candidates, seconds = NA_real_, stringsAsFactors = FALSE)
  
  for (i in seq_along(candidates)) {
    w <- candidates[i]
    future::plan(future::multisession, workers = w)
    t0 <- proc.time()[["elapsed"]]
    invisible(furrr::future_map(task_subset, task_fun, .options = furrr::furrr_options(seed = TRUE)))
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

run_name <- readline("\nEnter run name [default: DeltaDelta]: ")
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

# Run mode
cat("\nRun mode:\n  1) Single target vs all genes\n  2) ALL targets (all genes as targets; uses shortlist)\n")
run_mode <- readline("Choose 1 or 2 [default: 2]: ")
run_mode <- trimws(run_mode)
if (!nzchar(run_mode)) run_mode <- "2"
run_mode <- suppressWarnings(as.integer(run_mode))
if (!is.finite(run_mode) || is.na(run_mode) || !run_mode %in% c(1L, 2L)) run_mode <- 2L

target_gene <- NA_character_
if (run_mode == 1L) {
  target_gene <- readline("\nTarget gene (Source) [default: GRHPR]: ")
  if (!nzchar(target_gene)) target_gene <- "GRHPR"
}

max_abs_lag <- readline("\nMaximum |lag| in steps (e.g., 1 gives -1:1) [default: 1]: ")
max_abs_lag <- suppressWarnings(as.integer(trimws(max_abs_lag)))
if (!is.finite(max_abs_lag) || is.na(max_abs_lag) || max_abs_lag < 0) max_abs_lag <- 1L
lags_steps <- seq.int(-max_abs_lag, max_abs_lag, by = 1L)

sampling_interval_days <- readline("\nSampling interval in days (for lag_days reporting) [default: auto-infer]: ")
sampling_interval_days <- suppressWarnings(as.numeric(trimws(sampling_interval_days)))
if (!is.finite(sampling_interval_days) || is.na(sampling_interval_days) || sampling_interval_days <= 0) {
  sampling_interval_days <- NA_real_
}

# Permutations
perm_B <- readline("\nPermutation iterations for p-values [default: 2000]: ")
perm_B <- suppressWarnings(as.integer(trimws(perm_B)))
if (!is.finite(perm_B) || is.na(perm_B) || perm_B < 200) perm_B <- 2000L

set_seed <- readline("\nRandom seed [default: 1]: ")
set_seed <- suppressWarnings(as.integer(trimws(set_seed)))
if (!is.finite(set_seed) || is.na(set_seed)) set_seed <- 1L

# Shortlist configuration (applies in all-targets mode; optional in single-target mode)
shortlist_enable <- TRUE
shortlist_method <- "topk"
shortlist_k <- 25L
shortlist_score_rule <- "max_over_lags"  # or "lag0_only"

if (run_mode == 2L) {
  cat("\nShortlist options (applies to ALL targets mode):\n")
  cat("  1) Top-K per target (default K=25)\n")
  cat("  2) No shortlist (NOT recommended unless gene count is small)\n")
  sm <- readline("Choose 1 or 2 [default: 1]: ")
  sm <- trimws(sm); if (!nzchar(sm)) sm <- "1"
  sm <- suppressWarnings(as.integer(sm))
  if (!is.finite(sm) || is.na(sm) || !sm %in% c(1L, 2L)) sm <- 1L
  
  if (sm == 2L) {
    shortlist_enable <- FALSE
  } else {
    shortlist_enable <- TRUE
    shortlist_method <- "topk"
    k_in <- readline("\nTop-K per target [default: 25]: ")
    k_in <- suppressWarnings(as.integer(trimws(k_in)))
    if (is.finite(k_in) && !is.na(k_in) && k_in >= 1L) shortlist_k <- as.integer(k_in) else shortlist_k <- 25L
    
    cat("\nShortlist scoring rule:\n  1) lag0 only\n  2) max over lags (default)\n")
    sr <- readline("Choose 1 or 2 [default: 2]: ")
    sr <- trimws(sr); if (!nzchar(sr)) sr <- "2"
    sr <- suppressWarnings(as.integer(sr))
    if (!is.finite(sr) || is.na(sr) || !sr %in% c(1L, 2L)) sr <- 2L
    shortlist_score_rule <- if (sr == 1L) "lag0_only" else "max_over_lags"
  }
}

# ---------------------------
# Read + validate input
# ---------------------------
df <- read_table_auto(in_path)

required_cols <- c("Day","SampleID","Source","fpkm")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

df$Day <- suppressWarnings(as.numeric(df$Day))
df$fpkm <- suppressWarnings(as.numeric(df$fpkm))
df$SampleID <- as.character(df$SampleID)
df$Source <- as.character(df$Source)

if (any(!is.finite(df$Day))) stop("Non-numeric Day values detected after coercion.")
if (any(!is.finite(df$fpkm))) stop("Non-numeric fpkm values detected after coercion.")

days <- sort(unique(df$Day))
if (length(days) < 4) stop("Need at least 4 distinct time points (days) to compute Δ–Δ reliably.")

if (is.na(sampling_interval_days)) {
  sampling_interval_days <- stats::median(diff(days))
}

# Transform before aggregation
if (transform_choice == 2L) {
  if (any(df$fpkm < 0, na.rm = TRUE)) stop("Negative fpkm values detected; cannot log-transform.")
  df$fpkm_transformed <- log2(df$fpkm + pseudocount)
} else {
  df$fpkm_transformed <- df$fpkm
}

# Aggregate within Day x Source across replicates
agg_fun <- if (daily_stat == "median") stats::median else base::mean

daily <- aggregate(df$fpkm_transformed,
                   by = list(Day = df$Day, Source = df$Source),
                   FUN = agg_fun, na.rm = TRUE)
names(daily)[names(daily) == "x"] <- "daily_value"

daily_out <- file.path(output_dir, "DailySummary_ByDayGene.csv")
write.csv(daily, daily_out, row.names = FALSE)

daily_wide <- reshape(daily, idvar = "Day", timevar = "Source", direction = "wide")
colnames(daily_wide) <- sub("^daily_value\\.", "", colnames(daily_wide))
daily_wide <- daily_wide[order(daily_wide$Day), , drop = FALSE]

all_genes <- unique(setdiff(colnames(daily_wide), "Day"))
if (length(all_genes) < 2) stop("Need at least 2 genes in wide matrix.")

# Δ time index (same for all genes)
day_vec <- daily_wide$Day
d_day <- day_vec[-1]  # Δ corresponds to this day index

# ---------------------------
# Shortlisting (Top-K per target) and inference computation
# ---------------------------

# For score computation, we use Pearson r on aligned Δ-series.
# Score = either |r(lag0)| or max_k |r(k)| across lags.
score_pair <- function(dX, dY, score_rule, lags_steps) {
  if (score_rule == "lag0_only") {
    al <- align_lag(dX, dY, d_day, 0L)
    st <- safe_cor_stats(al$x, al$y, method = "pearson")
    return(list(score = abs(st$r), best_lag = 0L, r_best = st$r, n_pairs = st$n))
  }
  # max over lags
  best <- list(score = -Inf, best_lag = NA_integer_, r_best = NA_real_, n_pairs = 0L)
  for (k in lags_steps) {
    al <- align_lag(dX, dY, d_day, k)
    st <- safe_cor_stats(al$x, al$y, method = "pearson")
    sc <- abs(st$r)
    if (is.finite(sc) && sc > best$score) {
      best$score <- sc
      best$best_lag <- as.integer(k)
      best$r_best <- st$r
      best$n_pairs <- st$n
    }
  }
  if (!is.finite(best$score)) best$score <- NA_real_
  best
}

# Compute full stats (Pearson + Spearman, param p, perm p) for one pair across all lags
compute_pair_all_lags <- function(target, gene, dY, dX, lags_steps, sampling_interval_days, perm_B, set_seed) {
  rows <- vector("list", length(lags_steps))
  for (i in seq_along(lags_steps)) {
    k <- lags_steps[i]
    al <- align_lag(dX, dY, d_day, k)
    x_al <- al$x; y_al <- al$y
    
    pear <- safe_cor_stats(x_al, y_al, method = "pearson")
    spear <- safe_cor_stats(x_al, y_al, method = "spearman")
    
    seed_p <- seed_from_key(set_seed, target, gene, k, "pearson")
    seed_s <- seed_from_key(set_seed, target, gene, k, "spearman")
    
    p_perm_pear <- if (pear$n >= 4 && is.finite(pear$r)) perm_p_cor(x_al, y_al, method = "pearson", B = perm_B, seed = seed_p) else NA_real_
    p_perm_spear <- if (spear$n >= 4 && is.finite(spear$r)) perm_p_cor(x_al, y_al, method = "spearman", B = perm_B, seed = seed_s) else NA_real_
    
    rows[[i]] <- data.frame(
      target = target,
      gene   = gene,
      lag_steps = as.integer(k),
      lag_days  = as.numeric(k) * sampling_interval_days,
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
# Build task list depending on run mode
# ---------------------------
if (run_mode == 1L) {
  if (!target_gene %in% all_genes) stop("Target gene not found: ", target_gene)
  target_genes <- target_gene
} else {
  target_genes <- all_genes
}

# ---------------------------
# Worker selection (suggest/manual/benchmark)
#   Tasks are targets (single target or all targets)
# ---------------------------
avail <- future::availableCores()
tasks_n <- length(target_genes)
workers_suggested <- suggest_workers(n_tasks = tasks_n, perm_B = perm_B, avail = avail)

cat("\nParallel configuration:\n")
cat("  Available cores: ", avail, "\n", sep = "")
cat("  Targets (tasks): ", tasks_n, "\n", sep = "")
cat("  Suggested workers: ", workers_suggested, "\n", sep = "")
cat("  Worker selection modes:\n")
cat("    1) Use suggested workers (default)\n")
cat("    2) Manual workers\n")
cat("    3) Auto-benchmark (pilot subset of targets)\n")

w_mode <- readline("Choose 1, 2, or 3 [default: 1]: ")
w_mode <- trimws(w_mode)
if (!nzchar(w_mode)) w_mode <- "1"
w_mode <- suppressWarnings(as.integer(w_mode))
if (!is.finite(w_mode) || is.na(w_mode) || !w_mode %in% c(1L, 2L, 3L)) w_mode <- 1L

workers_in <- workers_suggested

# Define the per-target task function now (used by benchmark too)
# This function returns:
#  - In single-target mode: full results for all genes
#  - In all-targets mode: shortlist + full results for shortlisted partners
process_one_target <- function(tgt) {
  Y <- daily_wide[[tgt]]
  dY <- diff(Y)
  
  query_genes <- setdiff(all_genes, tgt)
  if (length(query_genes) == 0) return(NULL)
  
  # --- PASS 1 (shortlist) ---
  shortlist_df <- NULL
  if (run_mode == 2L && shortlist_enable) {
    scores <- data.frame(
      target = character(0),
      gene = character(0),
      score = numeric(0),
      best_lag_steps = integer(0),
      r_best = numeric(0),
      n_pairs_best = integer(0),
      stringsAsFactors = FALSE
    )
    
    # score each candidate gene
    # (No permutations here)
    for (g in query_genes) {
      X <- daily_wide[[g]]
      dX <- diff(X)
      sc <- score_pair(dX, dY, shortlist_score_rule, lags_steps)
      scores <- rbind(scores, data.frame(
        target = tgt,
        gene = g,
        score = sc$score,
        best_lag_steps = sc$best_lag,
        r_best = sc$r_best,
        n_pairs_best = sc$n_pairs,
        stringsAsFactors = FALSE
      ))
    }
    
    # order and keep topK
    scores <- scores[order(-scores$score, abs(scores$best_lag_steps), -scores$n_pairs_best), , drop = FALSE]
    if (nrow(scores) > shortlist_k) scores <- scores[seq_len(shortlist_k), , drop = FALSE]
    shortlist_df <- scores
    
    # shortlisted genes become query set for PASS 2
    query_genes <- shortlist_df$gene
  }
  
  # --- PASS 2 (inference on chosen pairs) ---
  # In single-target mode, query_genes is all genes except target (no shortlisting by default)
  # In all-targets mode, query_genes is shortlisted set (unless shortlist disabled)
  out_list <- vector("list", length(query_genes))
  for (j in seq_along(query_genes)) {
    g <- query_genes[j]
    X <- daily_wide[[g]]
    dX <- diff(X)
    out_list[[j]] <- compute_pair_all_lags(tgt, g, dY, dX, lags_steps, sampling_interval_days, perm_B, set_seed)
  }
  res_df <- do.call(rbind, out_list)
  
  # Attach shortlist score metadata (if available) to each row (helps downstream)
  if (!is.null(shortlist_df) && nrow(shortlist_df) > 0) {
    # map score fields onto results by gene
    m <- match(res_df$gene, shortlist_df$gene)
    res_df$shortlist_score <- shortlist_df$score[m]
    res_df$shortlist_best_lag_steps <- shortlist_df$best_lag_steps[m]
  } else {
    res_df$shortlist_score <- NA_real_
    res_df$shortlist_best_lag_steps <- NA_integer_
  }
  
  list(shortlist = shortlist_df, results = res_df)
}

if (w_mode == 2L) {
  w_in <- readline(paste0("Number of workers [default: ", workers_suggested, "]: "))
  w_in <- suppressWarnings(as.integer(trimws(w_in)))
  if (is.finite(w_in) && !is.na(w_in) && w_in >= 1L) {
    workers_in <- min(as.integer(w_in), avail, tasks_n)
  } else {
    workers_in <- workers_suggested
  }
} else if (w_mode == 3L) {
  set.seed(set_seed)
  pilot_n <- min(10L, tasks_n)
  pilot_targets <- sample(target_genes, pilot_n)
  
  cands <- unique(c(1L, 2L, 4L, 6L, 8L, 10L, 12L, 14L, 16L))
  cands <- cands[cands <= avail]
  cands <- cands[cands <= tasks_n]
  if (length(cands) == 0) cands <- 1L
  
  cat("\nAuto-benchmarking workers on pilot targets (n=", pilot_n, ")...\n", sep = "")
  bench <- benchmark_workers(pilot_targets, function(tt) process_one_target(tt), cands)
  print(bench)
  workers_in <- bench$workers[1]
  cat("Auto-selected workers: ", workers_in, "\n", sep = "")
}

cat("Using workers: ", workers_in, "\n", sep = "")
future::plan(future::multisession, workers = workers_in)
furrr::furrr_options(seed = TRUE)

# ---------------------------
# Run in parallel across targets
# ---------------------------
cat("\nRunning Δ–Δ cross-correlations...\n")
out_list <- furrr::future_map(target_genes, function(tt) process_one_target(tt),
                              .options = furrr::furrr_options(seed = TRUE))

# Combine outputs
shortlists <- lapply(out_list, function(x) if (!is.null(x)) x$shortlist else NULL)
shortlists <- shortlists[!vapply(shortlists, is.null, logical(1))]

results_list <- lapply(out_list, function(x) if (!is.null(x)) x$results else NULL)
results_list <- results_list[!vapply(results_list, is.null, logical(1))]

if (length(results_list) == 0) stop("No results were produced (unexpected).")

res <- do.call(rbind, results_list)

# Save shortlist table (all-targets mode)
shortlist_out <- NA_character_
if (run_mode == 2L && shortlist_enable && length(shortlists) > 0) {
  shortlist_df_all <- do.call(rbind, shortlists)
  shortlist_out <- file.path(output_dir, "Shortlist_TopK_ByTarget.csv")
  write.csv(shortlist_df_all, shortlist_out, row.names = FALSE)
}

# ---------------------------
# FDR control
# Recommended: adjust within each (method × lag) family across all tested pairs
# ---------------------------
res$q_bh_pearson <- ave(res$p_perm_pearson, res$lag_steps,
                        FUN = function(p) p.adjust(p, method = "BH"))
res$q_bh_spearman <- ave(res$p_perm_spearman, res$lag_steps,
                         FUN = function(p) p.adjust(p, method = "BH"))

# ---------------------------
# Write outputs
# ---------------------------
out_all <- file.path(output_dir, "DeltaDelta_XCorr_AllTests.csv")
write.csv(res, out_all, row.names = FALSE)

res_lag0 <- res[res$lag_steps == 0, , drop = FALSE]
out_lag0 <- file.path(output_dir, "DeltaDelta_XCorr_Lag0.csv")
write.csv(res_lag0, out_lag0, row.names = FALSE)

# Best lag per (target,gene) by |Pearson r|
best_by_absr <- do.call(rbind, lapply(split(res, paste(res$target, res$gene, sep = "||")), function(d) {
  d <- d[order(-abs(d$r_pearson), abs(d$lag_steps), -d$n_pairs), , drop = FALSE]
  d[1, , drop = FALSE]
}))
out_best <- file.path(output_dir, "DeltaDelta_XCorr_BestByAbsR.csv")
write.csv(best_by_absr, out_best, row.names = FALSE)

# ---------------------------
# Minimal figures (base R)
# ---------------------------
fig_dir <- file.path(output_dir, "Figures")
safe_dir_create(fig_dir)

# If single-target mode, plot target trajectory as before
if (run_mode == 1L && !is.na(target_gene) && target_gene %in% all_genes) {
  Y <- daily_wide[[target_gene]]
  png(file.path(fig_dir, paste0("Target_", target_gene, "_Daily_", daily_stat, "_", transform_label, ".png")),
      width = 1600, height = 900, res = 200)
  plot(day_vec, Y, type = "b",
       xlab = "Day",
       ylab = paste0(target_gene, " (daily ", daily_stat, ", ", transform_label, ")"),
       main = paste0("Target trajectory: ", target_gene))
  dev.off()
}

png(file.path(fig_dir, "Lag0_PearsonR_Hist.png"), width = 1600, height = 900, res = 200)
hist(res_lag0$r_pearson, breaks = 40,
     main = "Lag 0 Δ–Δ Pearson r (tested pairs)", xlab = "r")
dev.off()

png(file.path(fig_dir, "Lag0_SpearmanR_Hist.png"), width = 1600, height = 900, res = 200)
hist(res_lag0$r_spearman, breaks = 40,
     main = "Lag 0 Δ–Δ Spearman rho (tested pairs)", xlab = "rho")
dev.off()

# ---------------------------
# Quarto report generation (WD-safe canonical pattern)
# ---------------------------
script_path <- get_script_path()
script_path_norm <- if (!is.na(script_path) && nzchar(script_path)) normalize_slash(script_path, mustWork = FALSE) else NA_character_
script_name <- if (!is.na(script_path_norm)) basename(script_path_norm) else "DeltaDelta_XCorr.R"

# Extract header comment block
header_block <- NA_character_
if (!is.na(script_path_norm) && file.exists(script_path_norm)) {
  L <- readLines(script_path_norm, warn = FALSE)
  Lh <- L[seq_len(min(250, length(L)))]
  keep <- logical(length(Lh))
  for (i in seq_along(Lh)) {
    line <- Lh[i]
    if (i == 1 && grepl("^#!", line)) { keep[i] <- TRUE; next }
    if (grepl("^\\s*#", line) || grepl("^\\s*$", line)) keep[i] <- TRUE else break
  }
  header_block <- paste(Lh[keep], collapse = "\n")
}

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
  paste0("Run mode: ", if (run_mode == 1L) "Single target" else "ALL targets"),
  paste0("Target gene: ", if (run_mode == 1L) target_gene else "NA (all targets)"),
  paste0("Expression transform: ", transform_label),
  paste0("Daily summary: ", daily_stat),
  paste0("Max |lag| steps: ", max_abs_lag),
  paste0("Sampling interval days (for lag_days): ", sampling_interval_days),
  paste0("Permutation B: ", perm_B),
  paste0("Seed: ", set_seed),
  paste0("Workers: ", workers_in),
  paste0("Shortlist enabled: ", if (run_mode == 2L) shortlist_enable else "NA"),
  paste0("Shortlist K: ", if (run_mode == 2L && shortlist_enable) shortlist_k else "NA"),
  paste0("Shortlist score rule: ", if (run_mode == 2L && shortlist_enable) shortlist_score_rule else "NA"),
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
  "Let $X_t$ be the per-day summary (median or mean) for a gene and $Y_t$ for a target gene (after optional transformation). First differences are computed across successive sampled days:",
  "",
  "$$\\Delta X_t = X_t - X_{t-1}, \\quad \\Delta Y_t = Y_t - Y_{t-1}.$$",
  "",
  "For lag $k$ (in sampling steps), the Δ–Δ cross-correlation is computed as:",
  "",
  "$$\\mathrm{cor}(\\Delta X_t, \\Delta Y_{t+k}).$$",
  "",
  "Pearson and Spearman correlations are reported with permutation p-values and BH-adjusted q-values (BH applied within each lag for each method).",
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
  if (!is.na(shortlist_out)) paste0("- `", basename(shortlist_out), "`: shortlist (Top-K by target)") else "- (No shortlist file written)",
  paste0("- `", basename(out_all), "`: all Δ–Δ tests (pairs × lags)"),
  paste0("- `", basename(out_lag0), "`: lag 0 subset"),
  paste0("- `", basename(out_best), "`: best lag per (target,gene) by |Pearson r|"),
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

old_wd <- getwd()
setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)

writeLines(report_lines, qmd_base)
quarto::quarto_render(qmd_base)

html_path <- sub("\\.qmd$", ".html", qmd_file)

# ---------------------------
# Final console output
# ---------------------------
cat("\n\n=== DONE ===\n")
cat("Daily summary (long):\n  ", daily_out, "\n", sep = "")
if (!is.na(shortlist_out)) cat("Shortlist:\n  ", shortlist_out, "\n", sep = "")
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
