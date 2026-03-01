#!/usr/bin/env Rscript
# ======================================================================
# Transcriptome Temporal Structure Pipeline (Small-T safe; descriptive, not over-claiming)
#
# PURPOSE
#   This script summarizes how transcriptome-wide variability and coordination change
#   across a small number of developmental timepoints ("small T", e.g., ~9 days).
#   It outputs day-level metrics and several complementary, explicitly-defined
#   “important day” candidates, each tied to a distinct statistical concept:
#
#     A) Canalization day (variance minimum)
#        - Defined as the day where a global variability metric is smallest.
#        - Operationally: argmin across days of the *median* gene-level CV or IQR.
#
#     B) Transition / inflection day (best two-slope piecewise linear breakpoint)
#        - Defined as the day that best separates two linear regimes (before/after),
#          with continuity enforced at the breakpoint.
#        - Operationally: the day minimizing SSE over interior candidate breakpoints.
#
#     C) Excursion day (largest one-step jump)
#        - Defined as the day following the largest absolute one-step change in a
#          metric (captures impulse-like shifts).
#        - Operationally: argmax |Δ metric|, assigned to the later day of the jump.
#
#     D) AMOC two-regime split day (mean/variance step split; reported cautiously)
#        - A single change-point day from changepoint::cpt.meanvar(method="AMOC").
#        - This is reported as a *two-regime split candidate*, not a definitive
#          “critical day”, because small T makes change-point inference fragile.
#
# INPUT FORMAT (tidy / long)
#   Required columns (case-sensitive):
#     - Day        : integer-like day label (e.g., 4,6,8,...)
#     - SampleID   : biological replicate identifier within Day
#     - Gene       : feature identifier (gene symbol/ID)
#     - Abundance  : numeric expression value (e.g., FPKM/TPM/count-derived)
#
#   Each row is one Gene measurement for one SampleID at one Day.
#
# KEY COMPUTATIONS
#   1) Per-gene × day summaries (computed from replicate values):
#        - mean
#        - CV  = sd(x) / mean(x)      (only when mean>0 and sufficient non-zeros)
#        - IQR = IQR(x)
#
#   2) Global day-level metrics (robust across genes):
#        - median_cv  : median across genes of gene-level CV at that day
#        - median_iqr : median across genes of gene-level IQR at that day
#        - eff_dim    : "effective dimensionality" (participation ratio) of the
#                      within-day sample×gene standardized matrix.
#                      Computed from eigenvalues of the covariance implied by SVD:
#                        Neff = ( (Σ λ)^2 ) / (Σ λ^2 )
#                      where λ are eigenvalues of the standardized data covariance.
#
#   3) Uncertainty quantification (PRIMARY INFERENCE)
#        - Replicate-stratified bootstrap: resample SampleID *within each day*
#          with replacement; recompute day metrics and important-day calls.
#        - Output support for each candidate day is reported as frequency / B.
#
# USER CONTROLS (interactive prompts)
#   - B_TOTAL      : number of bootstrap replicates (default 200)
#   - USE_LOG1P    : apply log1p(Abundance) transform (default TRUE)
#   - MIN_NONZERO  : minimum non-zero replicate count required for CV (default 2)
#   - AMOC penalty : MBIC/BIC/AIC/MANUAL (default MBIC)
#   - MINSEGLEN    : AMOC minimum segment length (default 2)
#
# HOW TO RUN
#   Interactive (will prompt for input file and output directory):
#     Rscript TemporalStructure.R
#
#   Non-interactive (recommended for reproducibility):
#     Rscript TemporalStructure.R infile=/path/to/input.csv outdir=/path/to/output_parent
#
# OUTPUTS (created under a timestamped folder: TemporalStructure_YYYYMMDD_HHMMSS/)
#
#   Core tables:
#     1) Table_1_Daily_Metrics.csv
#        Columns: Day, median_cv, median_iqr, eff_dim
#
#     2) Table_1b_Daily_Metrics_BootstrapCI.csv
#        Columns: Day, Metric, lo, mid, hi
#        (bootstrap 2.5%, 50%, 97.5% quantiles per day/metric)
#
#     3) Bootstrap_Raw_AllStats.csv
#        Columns: b, stat, metric, day
#
#   “Important day” support tables:
#     4) Table_2_Canalization_Minima_Bootstrap.csv
#     5) Table_3_Piecewise_Inflection_Bootstrap.csv
#     6) Table_4_Excursion_Day_Bootstrap.csv
#     7) Table_5_AMOC_TwoRegimeSplit_Bootstrap.csv
#
#   Driver-gene tables (only written when relevant day is defined):
#     8) Table_6_DriverGenes_Canalization_D14.csv
#     9) Table_7_DriverGenes_Transition_D10.csv
#
#   Plots:
#     10) Plots/Metric_Trends_Annotated.png
#         IMPORTANT (as implemented below):
#           - Pointranges + connecting line: bootstrap median (mid) ± 95% CI (lo, hi)
#           - Observed series is NOT plotted (removed), per user request.
#           - Vertical lines: observed important-day calls (linetype-coded).
#
#   Log:
#     11) LOG_YYYYMMDD_HHMMSS.txt
#
# IMPORTANT CAVEATS (Small-T framing)
#   - All “important day” calls are operational definitions, not mechanistic claims.
#   - Interpret robustness primarily via bootstrap support (Tables 2–5) and CI widths (Table 1b).
# ======================================================================

suppressPackageStartupMessages({
  pkgs <- c("data.table","matrixStats","changepoint","ggplot2","quarto","knitr")
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")
  library(data.table)
  library(matrixStats)
  library(changepoint)
  library(ggplot2)
})

# -----------------------------
# Helpers
# -----------------------------
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

arg_val <- function(flag, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

is_nonempty_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && isTRUE(nzchar(x))
}

stop2 <- function(...) stop(paste0(...), call. = FALSE)

prompt_int <- function(msg, default) {
  cat(msg, " [default=", default, "]: ", sep="")
  x <- trimws(readline())
  if (!is_nonempty_string(x)) return(as.integer(default))
  v <- suppressWarnings(as.integer(x))
  if (!is.finite(v) || v < 0) return(as.integer(default))
  v
}

prompt_yesno <- function(msg, default_yes = TRUE) {
  def <- if (default_yes) "Y" else "N"
  cat(msg, " (Y/N) [default=", def, "]: ", sep="")
  x <- toupper(trimws(readline()))
  if (!is_nonempty_string(x)) return(default_yes)
  x %in% c("Y","YES")
}

choose_output_dir <- function() {
  outdir_arg <- arg_val("outdir", NA_character_)
  if (is_nonempty_string(outdir_arg)) return(outdir_arg)
  
  # 1) RStudio directory chooser (if available)
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    d <- tryCatch(rstudioapi::selectDirectory(caption = "Choose output folder"),
                  error = function(e) "")
    if (is_nonempty_string(d)) return(normalizePath(d, winslash="/", mustWork=FALSE))
  }
  
  # 2) macOS native Finder chooser (works well for Rscript in Terminal)
  if (Sys.info()[["sysname"]] == "Darwin") {
    d <- tryCatch({
      cmd <- c(
        "-e",
        'try',
        "-e",
        'set p to POSIX path of (choose folder with prompt "Choose output folder")',
        "-e",
        'on error',
        "-e",
        'set p to ""',
        "-e",
        'end try',
        "-e",
        'return p'
      )
      x <- system2("osascript", cmd, stdout = TRUE, stderr = TRUE)
      trimws(paste(x, collapse = "\n"))
    }, error = function(e) "")
    
    if (is_nonempty_string(d)) return(normalizePath(d, winslash="/", mustWork=FALSE))
  }
  
  # 3) Tk chooser (may feel less native; can fail in some headless contexts)
  if (requireNamespace("tcltk", quietly = TRUE)) {
    d <- tryCatch(tcltk::tclvalue(tcltk::tk_choose.dir(caption = "Choose output folder")),
                  error = function(e) "")
    if (is_nonempty_string(d)) return(normalizePath(d, winslash="/", mustWork=FALSE))
  }
  
  # 4) Last resort: manual path
  cat("\nPaste an output directory path, e.g. /Users/you/Desktop/OUT :\n> ")
  d <- trimws(readline())
  if (!is_nonempty_string(d)) stop2("No output directory provided.")
  normalizePath(d, winslash="/", mustWork=FALSE)
}


make_logger <- function(path) {
  force(path)
  function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- paste0("[", stamp, "] ", msg)
    cat(line, "\n", sep="")
    cat(line, "\n", sep="", file=path, append=TRUE)
  }
}

# -----------------------------
# Provenance / Reproducibility Tracking (NEW)
#   - Script path/name capture
#   - Input/output file registries
#   - QMD + HTML report generation in output directory
# -----------------------------
get_script_path <- function() {
  # 1) Most reliable for: Rscript /path/to/script.R
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = FALSE))
  }
  
  # 2) If called via source(".../script.R"), this is usually available inside the sourced file
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE))
  }
  
  # 3) RStudio interactive: active document path (if saved)
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  NA_character_
}

read_script_header <- function(script_path, n = 180L) {
  if (is.na(script_path) || !file.exists(script_path)) {
    return("Script path not available; header cannot be embedded.")
  }
  x <- readLines(script_path, warn = FALSE)
  paste(head(x, n), collapse = "\n")
}

# Registries
.input_files  <- character()
.output_files <- character()

register_input <- function(path) {
  .input_files <<- unique(c(.input_files, normalizePath(path, winslash = "/", mustWork = FALSE)))
  invisible(path)
}

register_output <- function(path) {
  .output_files <<- unique(c(.output_files, normalizePath(path, winslash = "/", mustWork = FALSE)))
  invisible(path)
}

# QMD writer (in output directory)
write_qmd_report <- function(qmd_path,
                             title,
                             script_path,
                             script_header_text,
                             input_files,
                             output_files,
                             output_dir,
                             dependencies) {
  md_list <- function(x) {
    if (!length(x)) return("_None registered_")
    paste0("- `", x, "`", collapse = "\n")
  }
  
  qmd <- c(
    "---",
    paste0("title: \"", title, "\""),
    "format:",
    "  html:",
    "    toc: true",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "## Metadata",
    paste0("- **Run timestamp:** `", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "`"),
    paste0("- **Script:** `", if (!is.na(script_path)) script_path else "NA", "`"),
    paste0("- **Output directory:** `", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "`"),
    "",
    "## Script header (embedded)",
    "```",
    script_header_text,
    "```",
    "",
    "## Inputs",
    md_list(input_files),
    "",
    "## Outputs / Generated files",
    md_list(output_files),
    "",
    "## Dependencies",
    paste0("- `", dependencies, "`", collapse = "\n"),
    "",
    "## Analytical logic",
    "This report summarizes the outputs produced by the pipeline and documents the operational definitions used to identify candidate “important days” under a small-T design.",
    "",
    "### Metrics",
    "- **Gene-level CV (within day):** `sd(x) / mean(x)` (computed when mean>0, sufficient non-zeros, and at least 2 replicates).",
    "- **Gene-level IQR (within day):** `IQR(x)`.",
    "- **Day-level median CV / IQR:** median across genes of gene-level CV/IQR for that day.",
    "- **Effective dimensionality (participation ratio):** computed from eigenvalues λ of the standardized within-day covariance implied by SVD:",
    "",
    "$$",
    "N_{\\mathrm{eff}} = \\frac{\\left(\\sum_i \\lambda_i\\right)^2}{\\sum_i \\lambda_i^2}",
    "$$",
    "",
    "### Statistical / inferential logic (small-T safe)",
    "- Primary uncertainty quantification uses replicate-stratified bootstrap resampling of SampleID within each day.",
    "- “Important day” candidates are operational definitions (argmin, piecewise breakpoint, max one-step jump, AMOC split) and should be interpreted via bootstrap support and CI widths.",
    "",
    "## Figures",
    "```{r}",
    "out_files <- c(",
    if (length(output_files)) paste0("  ", paste(sprintf("\"%s\"", gsub("\\\\", "/", output_files)), collapse = ",\n  ")) else "  character()",
    ")",
    "img_files <- out_files[file.exists(out_files) & grepl(\"\\\\.(png|jpg|jpeg|pdf)$\", out_files, ignore.case=TRUE)]",
    "if (length(img_files)) knitr::include_graphics(img_files)",
    "```",
    "",
    "## Interpretation of results",
    "Interpret day-level changes by triangulating: (i) shifts in bootstrap median trajectories, (ii) narrowing or widening of bootstrap CIs, and (iii) bootstrap support for candidate “important day” calls (Tables 2–5).",
    "",
    "## Reproducibility",
    "```{r}",
    "sessionInfo()",
    "```"
  )
  
  writeLines(qmd, con = qmd_path)
  invisible(qmd_path)
}

# Render QMD -> HTML into same directory
render_qmd_to_html <- function(qmd_path, output_dir) {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop2("Package 'quarto' is required to render .qmd. Install with install.packages('quarto').")
  }
  quarto::quarto_render(input = qmd_path, output_dir = output_dir)
  out_html <- file.path(output_dir, sub("\\.qmd$", ".html", basename(qmd_path)))
  invisible(out_html)
}

# -----------------------------
# Script identity (robust capture + fallback to known name)
# -----------------------------
resolve_script_path <- function() {
  # 1) RStudio editor context (most reliable when running interactively)
  p <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      rstudioapi::getSourceEditorContext()$path
    } else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) {
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  # 2) Rscript --file=...
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) {
      return(normalizePath(p2, winslash = "/", mustWork = FALSE))
    }
  }
  
  # 3) source() fallback (sometimes available)
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) {
    return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  }
  
  NA_character_
}

# Your actual script name (fallback for execution contexts where path cannot be detected)
known_script_filename <- "Change_Point_Directionv4ErrorBarsV2wthHTML.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_path_full <- resolve_script_path()

if (is.na(script_path_full)) {
  script_path <- NA_character_
  script_name <- known_script_stem
  # Header cannot be read without a path; keep a clear sentinel
  script_header_text <- "Script path not available; header cannot be embedded."
} else {
  script_path <- script_path_full
  script_name <- tools::file_path_sans_ext(basename(script_path_full))
  script_header_text <- read_script_header(script_path_full, n = 180L)
}

# AMOC step-change (reported as "two-regime split", not "critical day")
detect_amoc_split_day <- function(y, days, penalty="MBIC", minseglen=2L) {
  ok <- is.finite(y) & is.finite(days)
  y2 <- y[ok]; d2 <- days[ok]
  if (length(y2) < 5) return(NA_integer_)
  fit <- tryCatch(
    changepoint::cpt.meanvar(y2, method="AMOC", penalty=penalty, minseglen=minseglen, class=TRUE),
    error=function(e) NULL
  )
  if (is.null(fit)) return(NA_integer_)
  cp <- changepoint::cpts(fit)
  if (length(cp) < 1 || !is.finite(cp[1])) return(NA_integer_)
  as.integer(d2[cp[1]])
}

# Effective dimensionality (SVD-based) from samples x genes
calc_neff_from_matrix <- function(X) {
  if (is.null(X) || !is.matrix(X) || nrow(X) < 3 || ncol(X) < 2) return(NA_real_)
  X <- scale(X, center=TRUE, scale=TRUE)
  X[!is.finite(X)] <- 0
  s <- tryCatch(svd(X, nu=0, nv=0)$d, error=function(e) NULL)
  if (is.null(s) || length(s) == 0) return(NA_real_)
  ev <- (s^2) / (nrow(X) - 1)
  ev[ev < 1e-12] <- 0
  if (sum(ev) <= 0) return(NA_real_)
  (sum(ev)^2) / sum(ev^2)
}

make_day_matrix <- function(DT_day, genes_keep=NULL) {
  if (is.null(DT_day) || nrow(DT_day) == 0) return(NULL)
  dd <- DT_day
  if (!is.null(genes_keep)) dd <- dd[Gene %in% genes_keep]
  if (nrow(dd) == 0) return(NULL)
  W <- dcast(dd, SampleID ~ Gene, value.var="Abundance", fill=0)
  M <- as.matrix(W[, -1, with=FALSE])
  rownames(M) <- W$SampleID
  M
}

# --- Importance metrics (small-T safe) ---
pick_argmin_day <- function(y, days) {
  ok <- is.finite(y) & is.finite(days)
  if (!any(ok)) return(NA_integer_)
  days[ok][which.min(y[ok])]
}

piecewise_break_day <- function(y, days, min_points_each_side=3L) {
  ok <- is.finite(y) & is.finite(days)
  y <- y[ok]; days <- days[ok]
  n <- length(y)
  if (n < (2L*min_points_each_side)) return(NA_integer_)
  cand_idx <- seq(min_points_each_side, n - min_points_each_side)
  sse_best <- Inf
  best_idx <- NA_integer_
  
  for (k in cand_idx) {
    tb <- days[k]
    t <- days - tb
    Ileft <- as.numeric(t <= 0)
    Iright <- 1 - Ileft
    X <- cbind(1, t * Ileft, t * Iright)
    
    fit <- tryCatch(qr.solve(X, y), error=function(e) NULL)
    if (is.null(fit)) next
    yhat <- as.numeric(X %*% fit)
    sse <- sum((y - yhat)^2, na.rm=TRUE)
    if (is.finite(sse) && sse < sse_best) {
      sse_best <- sse
      best_idx <- k
    }
  }
  if (!is.finite(sse_best) || is.na(best_idx)) return(NA_integer_)
  as.integer(days[best_idx])
}

excursion_day <- function(y, days) {
  ok <- is.finite(y) & is.finite(days)
  y <- y[ok]; days <- days[ok]
  if (length(y) < 3) return(NA_integer_)
  dd <- abs(diff(y))
  as.integer(days[which.max(dd) + 1L])
}

# Bootstrap: replicate-stratified resampling within day
resample_DT_within_day <- function(DT, days_vec) {
  DT_by_day <- split(DT, DT$Day)
  out <- vector("list", length(days_vec))
  for (i in seq_along(days_vec)) {
    d <- days_vec[i]
    DD <- DT_by_day[[as.character(d)]]
    if (is.null(DD) || nrow(DD) == 0) {
      out[[i]] <- DT[0]
      next
    }
    sids <- unique(DD$SampleID)
    if (length(sids) == 0) {
      out[[i]] <- DT[0]
      next
    }
    draw <- sample(sids, size=length(sids), replace=TRUE)
    map <- data.table(SampleID = draw,
                      BootSampleID = paste0(draw, "__B", d, "_", seq_along(draw)))
    tmp <- DD[map, on="SampleID", allow.cartesian=TRUE]
    tmp[, SampleID := BootSampleID]
    tmp[, BootSampleID := NULL]
    out[[i]] <- tmp
  }
  rbindlist(out, use.names=TRUE)
}

# Compute per-day metrics from DT (median CV, median IQR, eff_dim)
compute_day_metrics <- function(DT, days_vec, MIN_NONZERO=2L, genes_univ=NULL, use_log1p=TRUE) {
  DT0 <- copy(DT)
  DT0[, Day := as.integer(Day)]
  DT0[, Abundance := as.numeric(Abundance)]
  DT0 <- DT0[is.finite(Day) & nzchar(SampleID) & nzchar(Gene) & is.finite(Abundance)]
  if (use_log1p) DT0[, Abundance := log1p(Abundance)]
  
  gd <- DT0[, {
    x <- Abundance
    nz <- sum(x != 0, na.rm=TRUE)
    m <- mean(x, na.rm=TRUE)
    s <- sd(x, na.rm=TRUE)
    list(
      cv  = if (is.finite(m) && m > 0 && nz >= MIN_NONZERO && .N >= 2) s/m else NA_real_,
      iqr = if (.N >= 2) IQR(x, na.rm=TRUE) else NA_real_
    )
  }, by=.(Day, Gene)]
  
  keep_days_chr <- as.character(days_vec)
  cvw  <- dcast(gd, Gene ~ Day, value.var="cv")
  iqrw <- dcast(gd, Gene ~ Day, value.var="iqr")
  
  for (dc in keep_days_chr) if (!dc %in% names(cvw))  cvw[, (dc) := NA_real_]
  for (dc in keep_days_chr) if (!dc %in% names(iqrw)) iqrw[, (dc) := NA_real_]
  setcolorder(cvw,  c("Gene", keep_days_chr))
  setcolorder(iqrw, c("Gene", keep_days_chr))
  
  cvm  <- as.matrix(cvw[, -1, with=FALSE]); rownames(cvm)  <- cvw$Gene;  colnames(cvm)  <- keep_days_chr
  iqrm <- as.matrix(iqrw[, -1, with=FALSE]); rownames(iqrm) <- iqrw$Gene; colnames(iqrm) <- keep_days_chr
  
  if (is.null(genes_univ)) genes_univ <- rownames(cvm)
  
  if (!all(genes_univ %in% rownames(cvm))) {
    miss <- setdiff(genes_univ, rownames(cvm))
    cvm <- rbind(cvm, matrix(NA_real_, nrow=length(miss), ncol=ncol(cvm),
                             dimnames=list(miss, colnames(cvm))))
  }
  if (!all(genes_univ %in% rownames(iqrm))) {
    miss <- setdiff(genes_univ, rownames(iqrm))
    iqrm <- rbind(iqrm, matrix(NA_real_, nrow=length(miss), ncol=ncol(iqrm),
                               dimnames=list(miss, colnames(iqrm))))
  }
  cvm  <- cvm[genes_univ, , drop=FALSE]
  iqrm <- iqrm[genes_univ, , drop=FALSE]
  
  m_cv  <- colMedians(cvm,  na.rm=TRUE)
  m_iqr <- colMedians(iqrm, na.rm=TRUE)
  
  DT_by_day <- split(DT0, DT0$Day)
  eff <- rep(NA_real_, length(days_vec))
  for (i in seq_along(days_vec)) {
    d <- days_vec[i]
    M <- make_day_matrix(DT_by_day[[as.character(d)]], genes_keep=genes_univ)
    eff[i] <- calc_neff_from_matrix(M)
  }
  
  list(
    day_metrics = data.table(Day=days_vec,
                             median_cv=as.numeric(m_cv),
                             median_iqr=as.numeric(m_iqr),
                             eff_dim=as.numeric(eff)),
    cv_mat = cvm,
    iqr_mat = iqrm,
    genes_univ = genes_univ
  )
}

support_table <- function(vec, name) {
  v <- vec[is.finite(vec)]
  if (length(v) == 0) return(data.table(stat=name, day=NA_integer_, support=0, n=0))
  tab <- sort(table(v), decreasing=TRUE)
  data.table(stat=name,
             day=as.integer(names(tab)),
             support=as.numeric(tab)/length(vec),
             n=as.integer(tab))
}

# -----------------------------
# Prompts
# -----------------------------
cat("\n=== Transcriptome Temporal Structure Pipeline (Small-T safe) ===\n")
cat("Input must be tidy/long: Day, SampleID, Gene, Abundance\n\n")

B_TOTAL <- prompt_int("Replicate-stratified bootstraps (within-day SampleID resampling)", 200L)
CORES   <- prompt_int("Cores to use (this rewrite uses base loop; keep modest)", 1L)
USE_LOG1P <- prompt_yesno("Apply log1p transform to Abundance?", TRUE)
MIN_NONZERO <- prompt_int("Min non-zero samples per Day×Gene for CV (recommended 2)", 2L)

cat("AMOC penalty (MBIC / BIC / AIC / Manual): [default=MBIC] > ")
PENALTY <- toupper(trimws(readline()))
if (!is_nonempty_string(PENALTY)) PENALTY <- "MBIC"
if (!PENALTY %in% c("MBIC","BIC","AIC","MANUAL")) PENALTY <- "MBIC"

MINSEGLEN <- prompt_int("AMOC minseglen (recommend 2 to avoid endpoint artifacts)", 2L)

# -----------------------------
# IO
# -----------------------------
infile_arg <- arg_val("infile", NA_character_)
infile <- if (is_nonempty_string(infile_arg)) infile_arg else {
  if (interactive()) file.choose() else stop2("Provide infile=/path/to/file.csv")
}
# Track input
register_input(infile)

out_parent <- choose_output_dir()
if (!dir.exists(out_parent)) dir.create(out_parent, recursive=TRUE, showWarnings=FALSE)

out_root <- file.path(out_parent, paste0("TemporalStructure_", ts_stamp()))
dir.create(out_root, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out_root, "Plots"), recursive=TRUE, showWarnings=FALSE)

logfile <- file.path(out_root, paste0("LOG_", ts_stamp(), ".txt"))
# Track log output
register_output(logfile)

logf <- make_logger(logfile)

logf("START")
logf(paste0("Script path: ", if (!is.na(script_path)) script_path else "NA"))
logf(paste0("Input file: ", infile))
logf(paste0("Output dir: ", out_root))
logf(paste0("Bootstraps: ", B_TOTAL, " | log1p: ", USE_LOG1P, " | MIN_NONZERO: ", MIN_NONZERO))
logf(paste0("AMOC penalty: ", PENALTY, " | minseglen: ", MINSEGLEN))

# -----------------------------
# Read + validate
# -----------------------------
DT <- fread(infile)
req <- c("Day","SampleID","Gene","Abundance")
miss <- setdiff(req, names(DT))
if (length(miss)) stop2("Missing required columns: ", paste(miss, collapse=", "))

DT[, Day := suppressWarnings(as.integer(gsub("[^0-9]+", "", as.character(Day))))]
DT[, SampleID := as.character(SampleID)]
DT[, Gene := as.character(Gene)]
DT[, Abundance := suppressWarnings(as.numeric(Abundance))]
DT <- DT[is.finite(Day) & nzchar(SampleID) & nzchar(Gene) & is.finite(Abundance)]
if (nrow(DT) == 0) stop2("No valid rows after cleaning.")

days_vec <- sort(unique(DT$Day))
if (length(days_vec) < 6) stop2("Too few unique days for reliable structure inference: ", length(days_vec))
logf(paste0("Unique days: ", paste(days_vec, collapse=", ")))

rep_counts <- DT[, uniqueN(SampleID), by=Day][order(Day)]
logf("Replicates per day (unique SampleID):")
for (i in seq_len(nrow(rep_counts))) logf(paste0("  Day ", rep_counts$Day[i], ": ", rep_counts$V1[i]))
if (any(rep_counts$V1 < 2)) logf("WARNING: days with <2 replicates will destabilize CV/eff_dim.")

# -----------------------------
# Observed metrics
# -----------------------------
logf("Computing observed daily metrics...")
obs <- compute_day_metrics(DT, days_vec, MIN_NONZERO=MIN_NONZERO, genes_univ=NULL, use_log1p=USE_LOG1P)
day_metrics <- obs$day_metrics

out_table1 <- file.path(out_root, "Table_1_Daily_Metrics.csv")
fwrite(day_metrics, out_table1)
register_output(out_table1)
logf("Wrote Table_1_Daily_Metrics.csv")

obs_canal_cv  <- pick_argmin_day(day_metrics$median_cv,  day_metrics$Day)
obs_canal_iqr <- pick_argmin_day(day_metrics$median_iqr, day_metrics$Day)

obs_infl_cv   <- piecewise_break_day(day_metrics$median_cv,  day_metrics$Day, min_points_each_side=3L)
obs_infl_iqr  <- piecewise_break_day(day_metrics$median_iqr, day_metrics$Day, min_points_each_side=3L)
obs_infl_eff  <- piecewise_break_day(day_metrics$eff_dim,    day_metrics$Day, min_points_each_side=3L)

obs_exc_cv    <- excursion_day(day_metrics$median_cv,  day_metrics$Day)
obs_exc_iqr   <- excursion_day(day_metrics$median_iqr, day_metrics$Day)
obs_exc_eff   <- excursion_day(day_metrics$eff_dim,    day_metrics$Day)

obs_amoc_cv   <- detect_amoc_split_day(day_metrics$median_cv,  day_metrics$Day, penalty=PENALTY, minseglen=MINSEGLEN)
obs_amoc_iqr  <- detect_amoc_split_day(day_metrics$median_iqr, day_metrics$Day, penalty=PENALTY, minseglen=MINSEGLEN)
obs_amoc_eff  <- detect_amoc_split_day(day_metrics$eff_dim,    day_metrics$Day, penalty=PENALTY, minseglen=MINSEGLEN)

logf(paste0("Observed canalization (argmin) CV=", obs_canal_cv, " | IQR=", obs_canal_iqr))
logf(paste0("Observed inflection (piecewise) CV=", obs_infl_cv, " | IQR=", obs_infl_iqr, " | eff_dim=", obs_infl_eff))
logf(paste0("Observed excursion (max|Δ|) CV=", obs_exc_cv, " | IQR=", obs_exc_iqr, " | eff_dim=", obs_exc_eff))
logf(paste0("Observed AMOC two-regime split CV=", obs_amoc_cv, " | IQR=", obs_amoc_iqr, " | eff_dim=", obs_amoc_eff))

# -----------------------------
# Bootstrap (replicate-stratified) — primary small-T inference
# -----------------------------
logf(paste0("Running replicate-stratified bootstrap: B=", B_TOTAL))
genes_univ <- obs$genes_univ

D <- length(days_vec)
boot_vals <- list(
  median_cv  = matrix(NA_real_, nrow=B_TOTAL, ncol=D),
  median_iqr = matrix(NA_real_, nrow=B_TOTAL, ncol=D),
  eff_dim    = matrix(NA_real_, nrow=B_TOTAL, ncol=D)
)
colnames(boot_vals$median_cv)  <- as.character(days_vec)
colnames(boot_vals$median_iqr) <- as.character(days_vec)
colnames(boot_vals$eff_dim)    <- as.character(days_vec)

boot <- data.table(
  b = integer(),
  stat = character(),
  metric = character(),
  day = integer()
)

for (b in seq_len(B_TOTAL)) {
  if (b %% 25 == 0) logf(paste0("  bootstrap ", b, "/", B_TOTAL))
  
  DTb <- resample_DT_within_day(DT, days_vec)
  bm <- compute_day_metrics(DTb, days_vec, MIN_NONZERO=MIN_NONZERO, genes_univ=genes_univ, use_log1p=USE_LOG1P)$day_metrics
  
  boot_vals$median_cv[b, ]  <- bm$median_cv
  boot_vals$median_iqr[b, ] <- bm$median_iqr
  boot_vals$eff_dim[b, ]    <- bm$eff_dim
  
  boot <- rbind(boot, data.table(b=b, stat="canalization", metric="median_cv",  day=pick_argmin_day(bm$median_cv, bm$Day)))
  boot <- rbind(boot, data.table(b=b, stat="canalization", metric="median_iqr", day=pick_argmin_day(bm$median_iqr, bm$Day)))
  
  boot <- rbind(boot, data.table(b=b, stat="inflection", metric="median_cv",  day=piecewise_break_day(bm$median_cv, bm$Day, 3L)))
  boot <- rbind(boot, data.table(b=b, stat="inflection", metric="median_iqr", day=piecewise_break_day(bm$median_iqr, bm$Day, 3L)))
  boot <- rbind(boot, data.table(b=b, stat="inflection", metric="eff_dim",    day=piecewise_break_day(bm$eff_dim, bm$Day, 3L)))
  
  boot <- rbind(boot, data.table(b=b, stat="excursion", metric="median_cv",  day=excursion_day(bm$median_cv, bm$Day)))
  boot <- rbind(boot, data.table(b=b, stat="excursion", metric="median_iqr", day=excursion_day(bm$median_iqr, bm$Day)))
  boot <- rbind(boot, data.table(b=b, stat="excursion", metric="eff_dim",    day=excursion_day(bm$eff_dim, bm$Day)))
  
  boot <- rbind(boot, data.table(b=b, stat="amoc_split", metric="median_cv",  day=detect_amoc_split_day(bm$median_cv,  bm$Day, penalty=PENALTY, minseglen=MINSEGLEN)))
  boot <- rbind(boot, data.table(b=b, stat="amoc_split", metric="median_iqr", day=detect_amoc_split_day(bm$median_iqr, bm$Day, penalty=PENALTY, minseglen=MINSEGLEN)))
  boot <- rbind(boot, data.table(b=b, stat="amoc_split", metric="eff_dim",    day=detect_amoc_split_day(bm$eff_dim,    bm$Day, penalty=PENALTY, minseglen=MINSEGLEN)))
}

out_boot_raw <- file.path(out_root, "Bootstrap_Raw_AllStats.csv")
fwrite(boot, out_boot_raw)
register_output(out_boot_raw)
logf("Wrote Bootstrap_Raw_AllStats.csv")

ci_dt <- rbindlist(lapply(names(boot_vals), function(met) {
  M <- boot_vals[[met]]
  data.table(
    Day = days_vec,
    Metric = met,
    lo = as.numeric(apply(M, 2, quantile, probs=0.025, na.rm=TRUE)),
    mid = as.numeric(apply(M, 2, quantile, probs=0.50,  na.rm=TRUE)),
    hi = as.numeric(apply(M, 2, quantile, probs=0.975, na.rm=TRUE))
  )
}), use.names=TRUE)

out_ci <- file.path(out_root, "Table_1b_Daily_Metrics_BootstrapCI.csv")
fwrite(ci_dt, out_ci)
register_output(out_ci)
logf("Wrote Table_1b_Daily_Metrics_BootstrapCI.csv")

summ <- boot[is.finite(day), .N, by=.(stat, metric, day)]
summ[, support := N / B_TOTAL]
setorder(summ, stat, metric, -support)

out_t2 <- file.path(out_root, "Table_2_Canalization_Minima_Bootstrap.csv")
out_t3 <- file.path(out_root, "Table_3_Piecewise_Inflection_Bootstrap.csv")
out_t4 <- file.path(out_root, "Table_4_Excursion_Day_Bootstrap.csv")
out_t5 <- file.path(out_root, "Table_5_AMOC_TwoRegimeSplit_Bootstrap.csv")

fwrite(summ[stat=="canalization"], out_t2); register_output(out_t2)
fwrite(summ[stat=="inflection"],   out_t3); register_output(out_t3)
fwrite(summ[stat=="excursion"],    out_t4); register_output(out_t4)
fwrite(summ[stat=="amoc_split"],   out_t5); register_output(out_t5)

logf("Wrote bootstrap support tables: Table_2..Table_5")

# -----------------------------
# Driver genes aligned to the RIGHT concepts (not “cp_day”)
# -----------------------------
logf("Computing driver genes aligned to canalization/transition (if days identified)...")

DT2 <- copy(DT)
if (USE_LOG1P) DT2[, Abundance := log1p(Abundance)]
gd_full <- DT2[, {
  x <- Abundance
  nz <- sum(x != 0, na.rm=TRUE)
  m <- mean(x, na.rm=TRUE)
  s <- sd(x, na.rm=TRUE)
  list(
    mean = m,
    cv   = if (is.finite(m) && m > 0 && nz >= MIN_NONZERO && .N >= 2) s/m else NA_real_,
    iqr  = if (.N >= 2) IQR(x, na.rm=TRUE) else NA_real_
  )
}, by=.(Day, Gene)]

wide_mat <- function(var) {
  W <- dcast(gd_full, Gene ~ Day, value.var=var)
  for (dc in as.character(days_vec)) if (!dc %in% names(W)) W[, (dc) := NA_real_]
  setcolorder(W, c("Gene", as.character(days_vec)))
  M <- as.matrix(W[, -1, with=FALSE])
  rownames(M) <- W$Gene
  colnames(M) <- as.character(days_vec)
  M
}

cvM   <- wide_mat("cv")
meanM <- wide_mat("mean")

canal_day <- obs_canal_cv
if (is.finite(canal_day)) {
  idx <- match(as.character(canal_day), colnames(cvM))
  if (!is.na(idx) && idx > 1) {
    drop_in <- cvM[, idx-1] - cvM[, idx]
    stays <- rep(TRUE, nrow(cvM))
    if (idx < ncol(cvM)) stays <- is.finite(cvM[, idx]) & is.finite(cvM[, idx+1]) & (cvM[, idx+1] <= cvM[, idx] + 1e-9)
    score <- ifelse(stays, drop_in, NA_real_)
    DTd <- data.table(Gene=rownames(cvM), canal_day=canal_day, drop_cv_into=as.numeric(score))
    DTd <- DTd[is.finite(drop_cv_into)]
    setorder(DTd, -drop_cv_into)
    DTd <- DTd[1:min(200L, .N)]
    out_t6 <- file.path(out_root, "Table_6_DriverGenes_Canalization_D14.csv")
    fwrite(DTd, out_t6)
    register_output(out_t6)
    logf(paste0("Wrote Table_6_DriverGenes_Canalization_D14.csv (canal_day=", canal_day, ")"))
  }
}

cand_infl <- c(obs_infl_cv, obs_infl_iqr)
infl_day <- suppressWarnings(as.integer(round(median(cand_infl[is.finite(cand_infl)]))))
if (is.finite(infl_day)) {
  idx <- match(as.character(infl_day), colnames(meanM))
  if (!is.na(idx) && idx > 1) {
    delta_mean <- meanM[, idx] - meanM[, idx-1]
    DTt <- data.table(Gene=rownames(meanM), infl_day=infl_day, delta_mean=as.numeric(delta_mean))
    DTt <- DTt[is.finite(delta_mean)]
    DTt[, abs_delta := abs(delta_mean)]
    setorder(DTt, -abs_delta)
    DTt <- DTt[1:min(200L, .N)]
    DTt[, abs_delta := NULL]
    out_t7 <- file.path(out_root, "Table_7_DriverGenes_Transition_D10.csv")
    fwrite(DTt, out_t7)
    register_output(out_t7)
    logf(paste0("Wrote Table_7_DriverGenes_Transition_D10.csv (infl_day=", infl_day, ")"))
  }
}

# -----------------------------
# Plots (annotated with the right concepts) + bootstrap uncertainty
#   USER REQUEST IMPLEMENTED HERE:
#     - Remove observed series entirely.
#     - Draw the line THROUGH the bootstrap medians (ci_dt_plot$mid).
#     - Keep 95% intervals as pointranges.
# -----------------------------
logf("Generating annotated plots (bootstrap median ± 95% CI; line through bootstrap medians)...")

ci_dt_plot <- copy(ci_dt)
ci_dt_plot[, Day := as.numeric(Day)]

anno <- data.table(
  Metric = c("median_cv","median_iqr","eff_dim",
             "median_cv","median_iqr","eff_dim",
             "median_cv","median_iqr","eff_dim",
             "median_cv","median_iqr","eff_dim"),
  stat = c(rep("canalization",3), rep("inflection",3), rep("excursion",3), rep("amoc_split",3)),
  day  = c(obs_canal_cv, obs_canal_iqr, NA_integer_,
           obs_infl_cv,  obs_infl_iqr,  obs_infl_eff,
           obs_exc_cv,   obs_exc_iqr,   obs_exc_eff,
           obs_amoc_cv,  obs_amoc_iqr,  obs_amoc_eff)
)
anno <- anno[is.finite(day)]
anno[, lab := stat]
anno[, day := as.numeric(day)]

p1 <- ggplot() +
  geom_pointrange(
    data = ci_dt_plot,
    aes(x = Day, y = mid, ymin = lo, ymax = hi),
    linewidth = 0.7
  ) +
  geom_line(
    data = ci_dt_plot,
    aes(x = Day, y = mid),
    linewidth = 1
  ) +
  geom_point(
    data = ci_dt_plot,
    aes(x = Day, y = mid),
    size = 2
  ) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  geom_vline(
    data = anno,
    aes(xintercept = day, linetype = lab),
    linewidth = 0.8
  ) +
  scale_linetype_manual(values = c(
    canalization = "solid",
    inflection   = "dotted",
    excursion    = "twodash",
    amoc_split   = "dashed"
  )) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Daily transcriptome structure metrics (bootstrap median ± 95% CI)",
    subtitle = "Points/line: bootstrap median (SampleID resampling within day). Pointranges: 95% bootstrap percentile interval.",
    x = "Day",
    y = NULL
  )

out_plot <- file.path(out_root, "Plots", "Metric_Trends_Annotated.png")
ggsave(
  filename = out_plot,
  plot = p1,
  width = 9,
  height = 11,
  dpi = 600
)
register_output(out_plot)
logf("Saved: Plots/Metric_Trends_Annotated.png")

# -----------------------------
# QMD + HTML Report (NEW): write + render into out_root
#   Checklist implemented:
#     - Write QMD in output_dir (out_root)
#     - YAML execute: echo: false
#     - Embed script header text (read from file when available)
#     - Add metadata (script, input, output)
#     - List dependencies
#     - List generated files
#     - Describe logic + formulas
#     - Insert plots with knitr::include_graphics()
#     - Interpretation section stub
#     - Reproducibility via sessionInfo()
#     - Render to HTML and confirm path
# -----------------------------
logf("Writing QMD + rendering HTML report...")

qmd_path <- file.path(out_root, "TemporalStructure_Report.qmd")
html_path <- file.path(out_root, "TemporalStructure_Report.html")

# Ensure output registries include the report artifacts
register_output(qmd_path)
register_output(html_path)

deps <- c("data.table","matrixStats","changepoint","ggplot2","quarto","knitr")

write_qmd_report(
  qmd_path = qmd_path,
  title = "Transcriptome Temporal Structure Pipeline Report",
  script_path = script_path,
  script_header_text = script_header_text,
  input_files = .input_files,
  output_files = .output_files,
  output_dir = out_root,
  dependencies = deps
)

# Render; Quarto names HTML based on input unless overridden—ensure we end with the expected name
out_html_actual <- tryCatch(render_qmd_to_html(qmd_path, out_root), error=function(e) {
  logf(paste0("ERROR rendering QMD: ", conditionMessage(e)))
  NA_character_
})

# If Quarto produced a differently-named HTML, track it and also copy to the expected filename
if (is.character(out_html_actual) && length(out_html_actual) == 1L && is_nonempty_string(out_html_actual) && file.exists(out_html_actual)) {
  if (normalizePath(out_html_actual, winslash="/", mustWork=FALSE) != normalizePath(html_path, winslash="/", mustWork=FALSE)) {
    ok_copy <- tryCatch(file.copy(out_html_actual, html_path, overwrite=TRUE), error=function(e) FALSE)
    if (isTRUE(ok_copy)) {
      logf(paste0("Copied HTML to: ", html_path))
    } else {
      logf("WARNING: Could not copy HTML to expected filename.")
    }
  }
  register_output(out_html_actual)
  logf(paste0("HTML report written to: ", normalizePath(html_path, winslash="/", mustWork=FALSE)))
} else {
  logf("WARNING: HTML report path could not be confirmed (render may have failed or named output unexpectedly).")
}

# -----------------------------
# Finish
# -----------------------------
logf("DONE")
cat("\n==============================================\n")
cat("COMPLETE\n")
cat("Results saved in:\n  ", out_root, "\n", sep="")
cat("Key outputs:\n")
cat("  - Table_1_Daily_Metrics.csv\n")
cat("  - Table_2_Canalization_Minima_Bootstrap.csv\n")
cat("  - Table_3_Piecewise_Inflection_Bootstrap.csv\n")
cat("  - Table_4_Excursion_Day_Bootstrap.csv\n")
cat("  - Table_5_AMOC_TwoRegimeSplit_Bootstrap.csv\n")
cat("Plots:\n")
cat("  - Plots/Metric_Trends_Annotated.png\n")
cat("Report:\n")
cat("  - TemporalStructure_Report.qmd\n")
cat("  - TemporalStructure_Report.html\n")
cat("Log:\n")
cat("  - ", basename(logfile), "\n", sep="")
cat("==============================================\n\n")
