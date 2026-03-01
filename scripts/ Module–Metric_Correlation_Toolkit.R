###############################################################################
# Module–Metric Correlation Runner (SOURCE-AND-RUN, Biologist-first, base R)
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES (general, not hardwired):
#   1) Prompts you to select an EXPRESSION matrix file (wide; genes x samples OR
#      samples x genes).
#   2) Prompts you to select a MODULE definition file (tidy: Module,Gene).
#   3) Optionally prompts you to select a METRICS file (tidy: SampleID, Day, metrics…).
#   4) Computes module scores per sample (mean-z by gene; optional PC1).
#   5) Correlates module scores vs selected metrics overall + within windows.
#   6) Saves tables + plots + log + manifest to an outputs folder.
#
# ----------------------------- INPUT FORMATS ---------------------------------
# (A) EXPRESSION MATRIX (CSV/TSV)  [REQUIRED]
#   Accepted formats:
#
#   Format A1 (genes in rows; samples in columns)  <-- like your IFN table
#     Gene,4_1,4_2,...,20_12
#     STAT1,9.9,43.8,...,15.1
#     ...
#
#   Format A2 (samples in rows; genes in columns)
#     SampleID,STAT1,STAT2,IRF7,...
#     4_1,9.9,2.4,6.1,...
#     ...
#
#   Requirements:
#     - First column is gene IDs (A1) or SampleID (A2), OR you will be prompted.
#     - Sample IDs may be "Day_BR" (e.g. "14_6") so Day can be parsed.
#
# (B) MODULE DEFINITIONS (CSV/TSV)  [REQUIRED]
#   Tidy mapping with at least:
#     Module,Gene
#     IFN,STAT1
#     Liver5SD,ALB
#     ...
#
# (C) METRICS TABLE (CSV/TSV)  [OPTIONAL]
#   Tidy mapping with at least:
#     SampleID,Day,ED,Neff,median_CV,...
#   Day optional (can be parsed from SampleID if missing).
#
# OUTPUTS:
#   ./outputs/ModuleMetricCorr_<timestamp>/
#     - ModuleScores_PerSample.csv
#     - ModuleCoverage.csv
#     - Correlations_Overall.csv
#     - Correlations_ByWindow.csv (if windows defined)
#     - Plots/*.png
#     - LOG_<timestamp>.txt
#     - manifest.json
#
# ----------------------------- HOW TO RUN ------------------------------------
#   Click "Source" in RStudio. The script will prompt for files and run.
###############################################################################

# =========================
# USER SETTINGS (edit if desired)
# =========================
RUN_ON_SOURCE <- TRUE                 # <--- keep TRUE for "click Source = run"
OUT_DIR_BASE  <- "./outputs"

# If you want defaults (no picker), you can set paths here; otherwise leave NULL:
DEFAULT_EXPR_FILE    <- NULL          # e.g. "BulkLiver_Expression.csv"
DEFAULT_MODULE_FILE  <- NULL          # e.g. "Modules_IFN_Liver5SD.csv"
DEFAULT_METRICS_FILE <- NULL          # e.g. "Constraint_Metrics.csv" or NULL

# Parsing sample IDs like "Day_BR"
PARSE_DAY_FROM_SAMPLEID <- TRUE
DAY_SEP <- "_"                        # "14_6" -> day 14

# Analysis choices
SCORE_METHOD <- "mean_z"              # "mean_z" or "pc1"
CORR_METHOD  <- "spearman"            # "spearman" or "pearson"
MIN_GENES_PRESENT <- 5
MAKE_PLOTS <- TRUE

# Optional: define phase windows (edit or leave empty list())
WINDOWS <- list(
  # Canalization = c(12,14),
  # Early = c(4,8,10),
  # Late = c(16,18,20)
)

# If metrics_file is provided and you want to pick metrics interactively, set TRUE
PICK_METRIC_COLS_INTERACTIVE <- TRUE

# =========================
# Helper functions (robust; base R)
# =========================
.now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

.guess_sep <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) return("\t")
  if (grepl("\\.csv$", path, ignore.case = TRUE)) return(",")
  txt <- readLines(path, n = 5, warn = FALSE)
  txt <- txt[nzchar(txt)]
  if (length(txt) == 0) return(",")
  c_count <- length(strsplit(txt[1], ",", fixed = TRUE)[[1]])
  t_count <- length(strsplit(txt[1], "\t", fixed = TRUE)[[1]])
  if (t_count > c_count) "\t" else ","
}

.read_table <- function(path, sep = NULL) {
  if (is.null(sep)) sep <- .guess_sep(path)
  df <- tryCatch(
    read.table(path, sep = sep, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, quote = "\"", comment.char = ""),
    error = function(e) e
  )
  if (inherits(df, "error")) stop("Failed to read file: ", path, "\n", df$message, call. = FALSE)
  df
}

.pick_file <- function(path = NULL, prompt = "Select a file") {
  if (!is.null(path) && nzchar(path) && file.exists(path)) return(normalizePath(path))
  message(prompt)
  if (!interactive()) stop("Non-interactive session: please set file paths in DEFAULT_* variables.", call. = FALSE)
  p <- file.choose()
  if (!file.exists(p)) stop("No valid file selected.", call. = FALSE)
  normalizePath(p)
}

.safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

.parse_day <- function(sample_ids, sep = "_") {
  parts <- strsplit(sample_ids, sep, fixed = TRUE)
  suppressWarnings(as.numeric(vapply(parts, function(x) x[1], character(1))))
}

# Convert genes-in-rows or samples-in-rows to internal standard:
#   X: genes x samples numeric matrix
#   gene_ids, sample_ids
.to_gene_by_sample_matrix <- function(expr_df) {
  
  has_Gene <- "Gene" %in% names(expr_df)
  has_SampleID <- "SampleID" %in% names(expr_df)
  
  orientation <- NULL
  id_col <- NULL
  
  if (has_Gene) { orientation <- "genes_in_rows"; id_col <- "Gene" }
  if (is.null(orientation) && has_SampleID) { orientation <- "samples_in_rows"; id_col <- "SampleID" }
  
  if (is.null(orientation)) {
    # Prompt user to pick which column is the ID
    message("I can't see a 'Gene' or 'SampleID' column.")
    message("Columns detected:\n  ", paste(names(expr_df), collapse = ", "))
    if (!interactive()) stop("Non-interactive: add 'Gene' or 'SampleID' column, or set it explicitly.", call. = FALSE)
    idx <- as.integer(readline("Enter the column number that contains IDs (Gene or SampleID): "))
    if (!is.finite(idx) || idx < 1 || idx > ncol(expr_df)) stop("Invalid column number.", call. = FALSE)
    id_col <- names(expr_df)[idx]
    # Heuristic: if remaining columns look like many samples, assume genes-in-rows
    orientation <- "genes_in_rows"
    warning("Assuming genes-in-rows with ID column: ", id_col, call. = FALSE)
  }
  
  if (orientation == "genes_in_rows") {
    gene_ids <- expr_df[[id_col]]
    mat_df <- expr_df[, setdiff(names(expr_df), id_col), drop = FALSE]
    for (j in seq_along(mat_df)) mat_df[[j]] <- suppressWarnings(as.numeric(mat_df[[j]]))
    X <- as.matrix(mat_df)
    rownames(X) <- gene_ids
    sample_ids <- colnames(X)
  } else {
    sample_ids <- expr_df[[id_col]]
    mat_df <- expr_df[, setdiff(names(expr_df), id_col), drop = FALSE]
    for (j in seq_along(mat_df)) mat_df[[j]] <- suppressWarnings(as.numeric(mat_df[[j]]))
    X <- t(as.matrix(mat_df))
    rownames(X) <- colnames(mat_df)
    colnames(X) <- sample_ids
    gene_ids <- rownames(X)
  }
  
  # Cleanup: non-finite -> NA
  X[!is.finite(X)] <- NA_real_
  
  list(X = X, gene_ids = gene_ids, sample_ids = sample_ids, orientation = orientation, id_col = id_col)
}

# Compute module scores
.compute_module_scores <- function(X, modules, score_method = "mean_z", min_genes_present = 5) {
  gene_universe <- rownames(X)
  
  # z-score per gene across samples (mean_z scoring uses this)
  gene_means <- rowMeans(X, na.rm = TRUE)
  gene_sds   <- apply(X, 1, sd, na.rm = TRUE)
  gene_sds[gene_sds == 0 | is.na(gene_sds)] <- NA_real_
  Xz <- sweep(X, 1, gene_means, "-")
  Xz <- sweep(Xz, 1, gene_sds, "/")
  
  scores <- data.frame(SampleID = colnames(X), stringsAsFactors = FALSE)
  cov <- data.frame(Module = names(modules),
                    Genes_in_module = vapply(modules, length, integer(1)),
                    Genes_present = vapply(modules, function(g) sum(unique(g) %in% gene_universe), integer(1)),
                    stringsAsFactors = FALSE)
  
  for (m in names(modules)) {
    g <- unique(modules[[m]])
    gp <- g[g %in% gene_universe]
    colname <- paste0(m, "_score")
    
    if (length(gp) < min_genes_present) {
      scores[[colname]] <- NA_real_
      next
    }
    
    if (score_method == "mean_z") {
      scores[[colname]] <- colMeans(Xz[gp, , drop = FALSE], na.rm = TRUE)
    } else {
      M <- t(X[gp, , drop = FALSE]) # samples x genes
      good <- apply(M, 2, function(v) sum(is.finite(v)) >= 3 && sd(v, na.rm = TRUE) > 0)
      M <- M[, good, drop = FALSE]
      if (ncol(M) < 2) {
        scores[[colname]] <- NA_real_
      } else {
        pc <- prcomp(M, center = TRUE, scale. = TRUE)
        scores[[colname]] <- pc$x[,1]
      }
    }
  }
  
  list(scores = scores, coverage = cov)
}

# Correlation helper
.corr_one <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 6) return(list(r = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  list(r = unname(ct$estimate), p = unname(ct$p.value), n = sum(ok))
}
.compute_state_metrics_by_day <- function(X, sample_ids, day_sep = "_",
                                          min_reps_per_day = 4,
                                          eps = 1e-8) {
  # X is genes x samples numeric matrix
  Day <- suppressWarnings(as.numeric(vapply(strsplit(sample_ids, day_sep, fixed = TRUE),
                                            function(x) x[1], character(1))))
  out <- list()
  
  uniq_days <- sort(unique(Day[is.finite(Day)]))
  if (length(uniq_days) == 0) stop("No valid Day parsed from sample IDs.", call. = FALSE)
  
  for (d in uniq_days) {
    idx <- which(Day == d)
    nrep <- length(idx)
    
    # Need replicates for within-day dispersion and eigen spectrum
    if (nrep < min_reps_per_day) {
      out[[length(out) + 1]] <- data.frame(
        Day = d, n_rep = nrep,
        median_CV = NA_real_, median_IQR = NA_real_, median_SD = NA_real_,
        Neff = NA_real_,
        median_resid_log_var = NA_real_, iqr_resid_log_var = NA_real_, sd_resid_log_var = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    Xd <- X[, idx, drop = FALSE]  # genes x reps
    
    # per-gene mean / sd / var within day
    mu  <- rowMeans(Xd, na.rm = TRUE)
    sdg <- apply(Xd, 1, sd, na.rm = TRUE)
    varg <- sdg^2
    
    # robust gene-level dispersion summaries
    iqr_g <- apply(Xd, 1, IQR, na.rm = TRUE)
    cv_g  <- sdg / (abs(mu) + eps)
    
    median_CV  <- median(cv_g, na.rm = TRUE)
    median_IQR <- median(iqr_g, na.rm = TRUE)
    median_SD  <- median(sdg,  na.rm = TRUE)
    
    # ----- ED-like residual dispersion: remove mean–variance dependence -----
    # log(var) ~ log(mean)
    log_mu  <- log(abs(mu) + eps)
    log_var <- log(varg + eps)
    
    ok <- is.finite(log_mu) & is.finite(log_var)
    if (sum(ok) >= 20) {
      fit <- lm(log_var[ok] ~ log_mu[ok])
      resid <- rep(NA_real_, length(log_var))
      resid[ok] <- resid(fit)
      median_resid <- median(resid[ok], na.rm = TRUE)
      iqr_resid    <- IQR(resid[ok], na.rm = TRUE)
      sd_resid     <- sd(resid[ok], na.rm = TRUE)
    } else {
      median_resid <- NA_real_
      iqr_resid    <- NA_real_
      sd_resid     <- NA_real_
    }
    
    # ----- Neff from eigenvalue spectrum (gene–gene correlation) -----
    # Observations = replicates; variables = genes
    # Use correlation across genes within day (may be large; handle with care)
    Neff <- NA_real_
    # Remove genes with zero variance or too many NA within day
    good_gene <- apply(Xd, 1, function(v) sum(is.finite(v)) >= 3 && sd(v, na.rm = TRUE) > 0)
    Xd2 <- Xd[good_gene, , drop = FALSE]
    
    if (nrow(Xd2) >= 5) {
      # genes x reps -> reps x genes
      M <- t(Xd2)
      # correlation among genes
      C <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
      # ensure finite
      if (all(is.finite(C))) {
        ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
        ev <- ev[ev > 0 & is.finite(ev)]
        if (length(ev) >= 2) {
          # participation ratio: (sum ev)^2 / sum(ev^2)
          Neff <- (sum(ev)^2) / sum(ev^2)
        }
      }
    }
    
    out[[length(out) + 1]] <- data.frame(
      Day = d, n_rep = nrep,
      median_CV = median_CV,
      median_IQR = median_IQR,
      median_SD = median_SD,
      Neff = Neff,
      median_resid_log_var = median_resid,
      iqr_resid_log_var = iqr_resid,
      sd_resid_log_var = sd_resid,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

# =========================
# MAIN RUNNER
# =========================
if (isTRUE(RUN_ON_SOURCE)) {
  
  stamp <- .now_stamp()
  out_dir <- file.path(OUT_DIR_BASE, paste0("ModuleMetricCorr_", stamp))
  plot_dir <- file.path(out_dir, "Plots")
  .safe_dir(out_dir)
  if (MAKE_PLOTS) .safe_dir(plot_dir)
  
  log_path <- file.path(out_dir, paste0("LOG_", stamp, ".txt"))
  log_con <- file(log_path, open = "wt")
  on.exit(try(close(log_con), silent = TRUE), add = TRUE)
  
  .log <- function(...) {
    msg <- paste0(..., collapse = "")
    writeLines(msg, log_con)
    cat(msg, "\n", sep = "")
  }
  
  .log("=== Module–Metric Correlation (SOURCE-AND-RUN) ===")
  .log("Timestamp: ", stamp)
  .log("Working directory: ", getwd())
  .log("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  
  # ---- File selection prompts ----
  expr_file   <- .pick_file(DEFAULT_EXPR_FILE,   "Select EXPRESSION matrix (CSV/TSV; genes×samples or samples×genes)")
  module_file <- .pick_file(DEFAULT_MODULE_FILE, "Select MODULE definition file (CSV/TSV; columns: Module,Gene)")
  
  metrics_file <- DEFAULT_METRICS_FILE
  if (interactive()) {
    ans <- readline("Do you have a METRICS table to include (ED/Neff/etc.)? [y/n]: ")
    if (tolower(substr(ans, 1, 1)) == "y") {
      metrics_file <- .pick_file(metrics_file, "Select METRICS file (CSV/TSV; SampleID, Day, metrics...)")
    } else {
      metrics_file <- NULL
    }
  }
  
  .log("Expression file: ", expr_file)
  .log("Module file: ", module_file)
  .log("Metrics file: ", ifelse(is.null(metrics_file), "<none>", metrics_file))
  
  # ---- Read inputs ----
  expr_df <- .read_table(expr_file)
  mod_df  <- .read_table(module_file)
  
  if (!all(c("Module","Gene") %in% names(mod_df))) {
    stop("Module file must contain columns: Module, Gene", call. = FALSE)
  }
  
  mod_df$Module <- trimws(mod_df$Module)
  mod_df$Gene   <- trimws(mod_df$Gene)
  mod_df <- mod_df[nzchar(mod_df$Module) & nzchar(mod_df$Gene), ]
  mod_df <- mod_df[!duplicated(mod_df[, c("Module","Gene")]), ]
  modules <- split(mod_df$Gene, mod_df$Module)
  
  conv <- .to_gene_by_sample_matrix(expr_df)
  X <- conv$X
  .log("Detected expression orientation: ", conv$orientation, " (ID col: ", conv$id_col, ")")
  .log("Expression matrix: ", nrow(X), " genes x ", ncol(X), " samples")
  
  # ---- Day parsing ----
  sample_ids <- colnames(X)
  Day <- rep(NA_real_, length(sample_ids))
  if (isTRUE(PARSE_DAY_FROM_SAMPLEID)) {
    Day <- .parse_day(sample_ids, sep = DAY_SEP)
    if (all(is.na(Day))) .log("WARNING: Could not parse Day from SampleID; Day will be NA.")
  }
  sample_meta <- data.frame(SampleID = sample_ids, Day = Day, stringsAsFactors = FALSE)
  
  # ---- OPTIONAL: Compute state metrics directly from expression (per Day) ----
  CALCULATE_STATE_METRICS <- TRUE
  
  state_metrics <- NULL
  if (isTRUE(CALCULATE_STATE_METRICS)) {
    .log("Computing per-Day state metrics from expression...")
    state_metrics <- .compute_state_metrics_by_day(
      X = X,
      sample_ids = colnames(X),
      day_sep = DAY_SEP,
      min_reps_per_day = 4
    )
    write.csv(state_metrics, file.path(out_dir, "StateMetrics_ByDay.csv"), row.names = FALSE)
    .log("State metrics saved: StateMetrics_ByDay.csv")
  }
  
  # ---- Metrics table ----
  metrics_df <- sample_meta
  metric_cols <- character()
  
  if (!is.null(metrics_file)) {
    met_df <- .read_table(metrics_file)
    if (!("SampleID" %in% names(met_df))) stop("Metrics file must contain SampleID column.", call. = FALSE)
    if (!("Day" %in% names(met_df))) {
      met_df$Day <- if (PARSE_DAY_FROM_SAMPLEID) .parse_day(met_df$SampleID, sep = DAY_SEP) else NA_real_
      .log("WARNING: Metrics file had no Day column; attempted parsing from SampleID.")
    }
    metrics_df <- merge(sample_meta, met_df, by = "SampleID", all.x = TRUE)
    
    # Choose metric columns
    cand <- setdiff(names(metrics_df), c("SampleID","Day"))
    # coerce candidates to numeric where possible
    for (nm in cand) metrics_df[[nm]] <- suppressWarnings(as.numeric(metrics_df[[nm]]))
    num_cand <- cand[vapply(metrics_df[, cand, drop = FALSE], function(v) sum(is.finite(v)) > 0, logical(1))]
    
    if (length(num_cand) == 0) {
      .log("WARNING: No numeric metric columns detected. You can still compute module scores and plots won't be made for metrics.")
    } else if (isTRUE(PICK_METRIC_COLS_INTERACTIVE) && interactive()) {
      .log("Metric columns detected: ", paste(num_cand, collapse = ", "))
      ans <- readline("Enter metric columns to analyze (comma-separated), or press Enter to use ALL detected: ")
      if (nzchar(ans)) {
        sel <- trimws(strsplit(ans, ",")[[1]])
        metric_cols <- intersect(sel, num_cand)
      } else {
        metric_cols <- num_cand
      }
    } else {
      metric_cols <- num_cand
    }
  }
  
  .log("Metric columns used: ", ifelse(length(metric_cols)==0, "<none>", paste(metric_cols, collapse = ", ")))
  
  # ---- Compute module scores ----
  sc <- .compute_module_scores(X, modules, score_method = SCORE_METHOD, min_genes_present = MIN_GENES_PRESENT)
  scores <- sc$scores
  coverage <- sc$coverage
  
  .log("Module coverage:")
  for (i in seq_len(nrow(coverage))) {
    .log("  - ", coverage$Module[i], ": ", coverage$Genes_present[i], "/", coverage$Genes_in_module[i])
  }
  
  full_df <- merge(metrics_df, scores, by = "SampleID", all.x = TRUE)
  
  # ---- Correlations ----
  score_cols <- grep("_score$", names(full_df), value = TRUE)
  # Attach day-level state metrics to each sample (so correlations/plots work)
  if (!is.null(state_metrics)) {
    full_df <- merge(full_df, state_metrics, by = "Day", all.x = TRUE)
    .log("Merged state metrics into sample table by Day.")
  }
  if (!is.null(state_metrics)) {
    auto_state_cols <- c("Neff","median_CV","median_IQR",
                         "median_resid_log_var","iqr_resid_log_var","sd_resid_log_var")
    metric_cols <- intersect(auto_state_cols, names(full_df))
    .log("Metric columns used (auto from state metrics): ", paste(metric_cols, collapse = ", "))
  }
  
  corr_overall <- NULL
  if (length(metric_cols) > 0) {
    rows <- list()
    for (scn in score_cols) {
      for (mc in metric_cols) {
        res <- .corr_one(full_df[[scn]], full_df[[mc]], method = CORR_METHOD)
        rows[[length(rows)+1]] <- data.frame(
          ModuleScore = scn, Metric = mc, Method = CORR_METHOD,
          R = res$r, P = res$p, N = res$n, stringsAsFactors = FALSE
        )
      }
    }
    corr_overall <- do.call(rbind, rows)
    corr_overall$P_adj_BH <- p.adjust(corr_overall$P, method = "BH")
  }
  
  corr_windows <- NULL
  if (length(metric_cols) > 0 && length(WINDOWS) > 0) {
    wrows <- list()
    for (wname in names(WINDOWS)) {
      days <- WINDOWS[[wname]]
      dfw <- full_df[full_df$Day %in% days, , drop = FALSE]
      for (scn in score_cols) {
        for (mc in metric_cols) {
          res <- .corr_one(dfw[[scn]], dfw[[mc]], method = CORR_METHOD)
          wrows[[length(wrows)+1]] <- data.frame(
            Window = wname, Days = paste(days, collapse = ","),
            ModuleScore = scn, Metric = mc, Method = CORR_METHOD,
            R = res$r, P = res$p, N = res$n, stringsAsFactors = FALSE
          )
        }
      }
    }
    corr_windows <- do.call(rbind, wrows)
    corr_windows$P_adj_BH <- p.adjust(corr_windows$P, method = "BH")
  }
  
  # ---- Plots ----
  if (MAKE_PLOTS && length(metric_cols) > 0) {
    day_vals <- full_df$Day
    day_levels <- sort(unique(day_vals[is.finite(day_vals)]))
    pal <- grDevices::rainbow(max(1, length(day_levels)))
    day_to_col <- setNames(pal, as.character(day_levels))
    pt_col <- ifelse(is.finite(day_vals), day_to_col[as.character(day_vals)], "grey70")
    
    for (scn in score_cols) {
      for (mc in metric_cols) {
        png(file.path(plot_dir, paste0(scn, "_vs_", mc, ".png")), width = 1100, height = 850)
        x <- full_df[[scn]]; y <- full_df[[mc]]
        ok <- is.finite(x) & is.finite(y)
        plot(x, y, pch = 19, col = pt_col, xlab = scn, ylab = mc,
             main = paste0(scn, " vs ", mc, " (", CORR_METHOD, ")"))
        if (sum(ok) >= 6) {
          fit <- tryCatch(stats::loess(y ~ x, subset = ok), error = function(e) NULL)
          if (!is.null(fit)) {
            xs <- seq(min(x[ok]), max(x[ok]), length.out = 200)
            ys <- predict(fit, newdata = data.frame(x = xs))
            lines(xs, ys, lwd = 3)
          }
          cr <- .corr_one(x, y, method = CORR_METHOD)
          legend("topleft", bty = "n",
                 legend = c(paste0("R = ", signif(cr$r, 3)),
                            paste0("P = ", signif(cr$p, 3)),
                            paste0("N = ", cr$n)))
        }
        if (length(day_levels) <= 12) {
          legend("bottomright", bty = "n",
                 legend = paste0("Day ", day_levels), col = pal, pch = 19, cex = 0.9)
        }
        dev.off()
      }
    }
  }
  
  # ---- Write outputs ----
  write.csv(scores,   file.path(out_dir, "ModuleScores_PerSample.csv"), row.names = FALSE)
  write.csv(coverage, file.path(out_dir, "ModuleCoverage.csv"), row.names = FALSE)
  if (!is.null(corr_overall)) write.csv(corr_overall, file.path(out_dir, "Correlations_Overall.csv"), row.names = FALSE)
  if (!is.null(corr_windows)) write.csv(corr_windows, file.path(out_dir, "Correlations_ByWindow.csv"), row.names = FALSE)
  
  # ---- Manifest (simple) ----
  manifest <- list(
    timestamp = stamp,
    expr_file = expr_file,
    module_file = module_file,
    metrics_file = ifelse(is.null(metrics_file), NA, metrics_file),
    score_method = SCORE_METHOD,
    corr_method = CORR_METHOD,
    min_genes_present = MIN_GENES_PRESENT,
    windows = WINDOWS,
    metric_cols = metric_cols,
    outputs = list(out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
                   log = basename(log_path))
  )
  # minimal JSON writer
  .to_json <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      nm <- names(x)
      if (is.null(nm)) {
        paste0("[", paste(vapply(x, .to_json, character(1)), collapse = ","), "]")
      } else {
        inner <- paste(vapply(seq_along(x), function(i) {
          paste0("\"", nm[i], "\":", .to_json(x[[i]]))
        }, character(1)), collapse = ",")
        paste0("{", inner, "}")
      }
    } else if (is.character(x)) {
      if (length(x) == 1) {
        if (is.na(x)) "null" else paste0("\"", gsub("\"", "\\\\\"", x), "\"")
      } else paste0("[", paste(vapply(x, function(s) if (is.na(s)) "null" else paste0("\"", gsub("\"","\\\\\"",s), "\""), character(1)), collapse=","), "]")
    } else if (is.numeric(x)) {
      if (length(x) == 1) if (is.na(x)) "null" else as.character(x)
      else paste0("[", paste(ifelse(is.na(x),"null",as.character(x)), collapse=","), "]")
    } else if (is.logical(x)) {
      if (length(x) == 1) if (is.na(x)) "null" else ifelse(x,"true","false")
      else paste0("[", paste(ifelse(is.na(x),"null",ifelse(x,"true","false")), collapse=","), "]")
    } else {
      paste0("\"", gsub("\"","\\\\\"", as.character(x)), "\"")
    }
  }
  writeLines(.to_json(manifest), con = file.path(out_dir, "manifest.json"))
  
  .log("DONE. Outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}
