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
#   Interactive (will prompt for input file):
#     Rscript TemporalStructure.R
#
#   Non-interactive (recommended for reproducibility):
#     Rscript TemporalStructure.R infile=/path/to/input.csv
#
# OUTPUT POLICY (MANDATORY)
#   All outputs are written under:
#     outputs/<run_id>/
#   where:
#     run_id = <analysis_name>_<Source>_<YYYYMMDD_HHMMSS>
#
# IMPORTANT CAVEATS (Small-T framing)
#   - All “important day” calls are operational definitions, not mechanistic claims.
#   - Interpret robustness primarily via bootstrap support (Tables 2–5) and CI widths (Table 1b).
# ======================================================================

# -----------------------------
# Packages (analysis + retrofit layer)
# -----------------------------
suppressPackageStartupMessages({
  pkgs <- c(
    "data.table","matrixStats","changepoint","ggplot2",
    # retrofit layer requirements
    "jsonlite","knitr","rmarkdown"
  )
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing)) {
    install.packages(missing, repos="https://cloud.r-project.org", quiet=TRUE)
  }
  library(data.table)
  library(matrixStats)
  library(changepoint)
  library(ggplot2)
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
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

prompt_string <- function(msg, default="") {
  cat(msg, " [default=", default, "]: ", sep="")
  x <- trimws(readline())
  if (!is_nonempty_string(x)) return(default)
  x
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
# Script identity capture (MANDATORY; provided block + required behavior)
# -----------------------------
resolve_script_path <- function() {
  # 1) RStudio editor context
  p <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      rstudioapi::getSourceEditorContext()$path
    } else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) {
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  # 2) Rscript --file=
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) {
      return(normalizePath(p2, winslash = "/", mustWork = FALSE))
    }
  }
  
  # 3) source() fallback (sometimes present)
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) {
    return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  }
  
  # 4) Unknown
  NA_character_
}

# If the script filename is known, set it here as a fallback:
known_script_filename <- "Change_Point_Directionv4ErrorBarsV2wthHTML.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  # Path cannot be detected in this execution mode; still record a valid script_name.
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

read_script_header <- function(script_full_path, n=220L) {
  if (is.na(script_full_path) || !file.exists(script_full_path)) {
    return("Script path not available; header cannot be embedded (script_full is NA or file not found).")
  }
  x <- readLines(script_full_path, warn = FALSE)
  paste(head(x, n), collapse = "\n")
}

script_header_text <- read_script_header(script_full, n=220L)

# -----------------------------
# Output policy (MANDATORY): outputs/<run_id> under current working directory
# -----------------------------
outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

# Registries (for manifest)
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

# Dependencies (packages + versions)
deps_table <- function(pkgs) {
  ok <- pkgs[pkgs %in% rownames(installed.packages())]
  data.table(
    package = ok,
    version = vapply(ok, function(p) as.character(utils::packageVersion(p)), character(1))
  )
}

write_manifest_json <- function(manifest, path) {
  jsonlite::write_json(manifest, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(path)
}

write_files_inventory <- function(output_dir, inventory_csv_path) {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, include.dirs = FALSE)
  files <- files[file.exists(files)]
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(output_dir, winslash="/", mustWork=FALSE)), "/?"), "", gsub("\\\\", "/", files))
  info <- file.info(files)
  dt <- data.table(
    relative_path = rel,
    full_path = gsub("\\\\", "/", files),
    bytes = as.numeric(info$size),
    modified = as.character(info$mtime)
  )
  data.table::fwrite(dt, inventory_csv_path)
  invisible(inventory_csv_path)
}

# -----------------------------
# Analysis functions (UNCHANGED)
# -----------------------------
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
  W <- data.table::dcast(dd, SampleID ~ Gene, value.var="Abundance", fill=0)
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
  data.table::rbindlist(out, use.names=TRUE)
}

# Compute per-day metrics from DT (median CV, median IQR, eff_dim)
compute_day_metrics <- function(DT, days_vec, MIN_NONZERO=2L, genes_univ=NULL, use_log1p=TRUE) {
  DT0 <- data.table::copy(DT)
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
  cvw  <- data.table::dcast(gd, Gene ~ Day, value.var="cv")
  iqrw <- data.table::dcast(gd, Gene ~ Day, value.var="iqr")
  
  for (dc in keep_days_chr) if (!dc %in% names(cvw))  cvw[, (dc) := NA_real_]
  for (dc in keep_days_chr) if (!dc %in% names(iqrw)) iqrw[, (dc) := NA_real_]
  data.table::setcolorder(cvw,  c("Gene", keep_days_chr))
  data.table::setcolorder(iqrw, c("Gene", keep_days_chr))
  
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
  
  m_cv  <- matrixStats::colMedians(cvm,  na.rm=TRUE)
  m_iqr <- matrixStats::colMedians(iqrm, na.rm=TRUE)
  
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

# -----------------------------
# QMD writer (MANDATORY) + render via rmarkdown::render (MANDATORY)
# -----------------------------
write_qmd_report <- function(qmd_path,
                             script_name,
                             run_id,
                             run_timestamp,
                             script_name_val,
                             script_path_val,
                             script_full_val,
                             input_files,
                             output_dir,
                             inventory_csv,
                             dependencies_dt,
                             plot_relpaths = character(0)) {
  
  deps_lines <- if (nrow(dependencies_dt)) {
    paste0("- `", dependencies_dt$package, "` (", dependencies_dt$version, ")", collapse = "\n")
  } else {
    "_No dependencies recorded_"
  }
  
  qmd <- c(
    "---",
    paste0("title: \"", script_name, " Report\""),
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
    paste0("- **Run ID:** `", run_id, "`"),
    paste0("- **Run timestamp:** `", run_timestamp, "`"),
    paste0("- **Script name:** `", script_name_val, "`"),
    paste0("- **Script path:** `", script_path_val, "`"),
    paste0("- **Script full path:** `", ifelse(is.na(script_full_val), "NA", script_full_val), "`"),
    paste0("- **Output directory:** `", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "`"),
    "",
    "## Inputs",
    if (length(input_files)) paste0("- `", input_files, "`", collapse = "\n") else "_None registered_",
    "",
    "## Script header (embedded)",
    "```",
    script_header_text,
    "```",
    "",
    "## Dependencies (packages + versions)",
    deps_lines,
    "",
    "## Generated files inventory",
    "```{r}",
    "inv <- data.table::fread(\"Project_Manifest_Files.csv\")",
    "knitr::kable(inv, format = \"html\")",
    "```",
    "",
    "## Analytical logic + formulas",
    "- **Gene-level CV (within day):** `sd(x) / mean(x)` (computed when mean>0, sufficient non-zeros, and at least 2 replicates).",
    "- **Gene-level IQR (within day):** `IQR(x)`.",
    "- **Day-level median CV / IQR:** median across genes of gene-level CV/IQR for that day.",
    "- **Effective dimensionality (participation ratio):** computed from eigenvalues λ of the standardized within-day covariance implied by SVD:",
    "",
    "$$",
    "N_{\\mathrm{eff}} = \\frac{\\left(\\sum_i \\lambda_i\\right)^2}{\\sum_i \\lambda_i^2}",
    "$$",
    "",
    "- **Primary uncertainty quantification:** replicate-stratified bootstrap resampling of SampleID within each day (small-T safe).",
    "",
    "## Figures",
    "```{r}",
    "img_files <- c(",
    if (length(plot_relpaths)) paste0("  ", paste(sprintf("\"%s\"", plot_relpaths), collapse = ",\n  ")) else "  character()",
    ")",
    "img_files <- img_files[file.exists(img_files)]",
    "if (length(img_files)) knitr::include_graphics(img_files)",
    "```",
    "",
    "## Interpretation of results",
    "Interpret day-level changes by triangulating: (i) shifts in bootstrap median trajectories, (ii) narrowing or widening of bootstrap CIs, and (iii) bootstrap support for candidate “important day” calls (Tables 2–5). Under small-T, treat all candidate days as operational definitions rather than mechanistic claims; emphasize robustness (support and CI widths) over point estimates.",
    "",
    "## Reproducibility",
    "```{r}",
    "sessionInfo()",
    "```"
  )
  
  writeLines(qmd, con = qmd_path)
  invisible(qmd_path)
}

render_qmd_to_html <- function(qmd_path, output_html_path) {
  # rmarkdown::render() will use Quarto when rendering .qmd if available.
  # This is a mandatory requirement in the packet.
  res <- tryCatch({
    rmarkdown::render(
      input = qmd_path,
      output_file = basename(output_html_path),
      output_dir = dirname(output_html_path),
      quiet = TRUE
    )
  }, error = function(e) e)
  
  res
}

# -----------------------------
# Console header
# -----------------------------
cat("\n=== Transcriptome Temporal Structure Pipeline (Small-T safe) ===\n")
cat("Input must be tidy/long: Day, SampleID, Gene, Abundance\n\n")

# Mandatory: note if script_full detection failed
if (is.na(script_full)) {
  cat("NOTE: Script full-path detection failed in this execution mode; using fallback filename.\n")
}

# -----------------------------
# Prompts (UNCHANGED, plus Source prompt required for deterministic run_id)
# -----------------------------
B_TOTAL <- prompt_int("Replicate-stratified bootstraps (within-day SampleID resampling)", 200L)
CORES   <- prompt_int("Cores to use (this rewrite uses base loop; keep modest)", 1L) # retained, not used
USE_LOG1P <- prompt_yesno("Apply log1p transform to Abundance?", TRUE)
MIN_NONZERO <- prompt_int("Min non-zero samples per Day×Gene for CV (recommended 2)", 2L)

cat("AMOC penalty (MBIC / BIC / AIC / Manual): [default=MBIC] > ")
PENALTY <- toupper(trimws(readline()))
if (!is_nonempty_string(PENALTY)) PENALTY <- "MBIC"
if (!PENALTY %in% c("MBIC","BIC","AIC","MANUAL")) PENALTY <- "MBIC"

MINSEGLEN <- prompt_int("AMOC minseglen (recommend 2 to avoid endpoint artifacts)", 2L)

# New (for run_id policy): Source string used in run_id
SOURCE <- prompt_string("Source label for run_id (e.g., liver, plasma, lineA)", "UnknownSource")
# sanitize source for filesystem
SOURCE <- gsub("[^A-Za-z0-9._-]+", "_", SOURCE)

# analysis_name string used in run_id (infer; can be overridden by arg analysis_name=...)
analysis_name_arg <- arg_val("analysis_name", NA_character_)
analysis_name <- if (is_nonempty_string(analysis_name_arg)) analysis_name_arg else "TemporalStructure"
analysis_name <- gsub("[^A-Za-z0-9._-]+", "_", analysis_name)

run_id <- paste0(analysis_name, "_", SOURCE, "_", ts_stamp())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Plots"), recursive = TRUE, showWarnings = FALSE)

# Manifest paths (MANDATORY)
manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
manifest_files_csv <- file.path(output_dir, "Project_Manifest_Files.csv")

# Register manifest placeholders early
register_output(manifest_json_path)
register_output(manifest_files_csv)

# Logging
logfile <- file.path(output_dir, paste0("LOG_", ts_stamp(), ".txt"))
register_output(logfile)
logf <- make_logger(logfile)

logf("START")
logf(paste0("run_id: ", run_id))
logf(paste0("outputs_root: ", normalizePath(outputs_root, winslash="/", mustWork=FALSE)))
logf(paste0("output_dir: ", normalizePath(output_dir, winslash="/", mustWork=FALSE)))

logf(paste0("script_name: ", script_name))
logf(paste0("script_path: ", script_path))
logf(paste0("script_full: ", ifelse(is.na(script_full), "NA", script_full)))

logf(paste0("Bootstraps: ", B_TOTAL, " | log1p: ", USE_LOG1P, " | MIN_NONZERO: ", MIN_NONZERO))
logf(paste0("AMOC penalty: ", PENALTY, " | minseglen: ", MINSEGLEN))
logf(paste0("Source: ", SOURCE))
logf(paste0("analysis_name: ", analysis_name))

# -----------------------------
# IO: input file (prompt if not provided)
# -----------------------------
infile_arg <- arg_val("infile", NA_character_)

if (is_nonempty_string(infile_arg)) {
  infile <- infile_arg
} else {
  if (interactive()) {
    cat("\nSelect the INPUT CSV (tidy/long: Day, SampleID, Gene, Abundance)\n")
    infile <- file.choose()
  } else {
    stop2("Provide infile=/path/to/file.csv (non-interactive mode cannot prompt).")
  }
}

register_input(infile)
logf(paste0("Input file: ", infile))

# -----------------------------
# Read + validate (UNCHANGED)
# -----------------------------
DT <- data.table::fread(infile)
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
# Observed metrics (UNCHANGED)
# -----------------------------
logf("Computing observed daily metrics...")
obs <- compute_day_metrics(DT, days_vec, MIN_NONZERO=MIN_NONZERO, genes_univ=NULL, use_log1p=USE_LOG1P)
day_metrics <- obs$day_metrics

out_table1 <- file.path(output_dir, "Table_1_Daily_Metrics.csv")
data.table::fwrite(day_metrics, out_table1)
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
# Bootstrap (replicate-stratified) — primary small-T inference (UNCHANGED)
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

out_boot_raw <- file.path(output_dir, "Bootstrap_Raw_AllStats.csv")
data.table::fwrite(boot, out_boot_raw)
register_output(out_boot_raw)
logf("Wrote Bootstrap_Raw_AllStats.csv")

ci_dt <- data.table::rbindlist(lapply(names(boot_vals), function(met) {
  M <- boot_vals[[met]]
  data.table(
    Day = days_vec,
    Metric = met,
    lo = as.numeric(apply(M, 2, quantile, probs=0.025, na.rm=TRUE)),
    mid = as.numeric(apply(M, 2, quantile, probs=0.50,  na.rm=TRUE)),
    hi = as.numeric(apply(M, 2, quantile, probs=0.975, na.rm=TRUE))
  )
}), use.names=TRUE)

out_ci <- file.path(output_dir, "Table_1b_Daily_Metrics_BootstrapCI.csv")
data.table::fwrite(ci_dt, out_ci)
register_output(out_ci)
logf("Wrote Table_1b_Daily_Metrics_BootstrapCI.csv")

summ <- boot[is.finite(day), .N, by=.(stat, metric, day)]
summ[, support := N / B_TOTAL]
data.table::setorder(summ, stat, metric, -support)

out_t2 <- file.path(output_dir, "Table_2_Canalization_Minima_Bootstrap.csv")
out_t3 <- file.path(output_dir, "Table_3_Piecewise_Inflection_Bootstrap.csv")
out_t4 <- file.path(output_dir, "Table_4_Excursion_Day_Bootstrap.csv")
out_t5 <- file.path(output_dir, "Table_5_AMOC_TwoRegimeSplit_Bootstrap.csv")

data.table::fwrite(summ[stat=="canalization"], out_t2); register_output(out_t2)
data.table::fwrite(summ[stat=="inflection"],   out_t3); register_output(out_t3)
data.table::fwrite(summ[stat=="excursion"],    out_t4); register_output(out_t4)
data.table::fwrite(summ[stat=="amoc_split"],   out_t5); register_output(out_t5)

logf("Wrote bootstrap support tables: Table_2..Table_5")

# -----------------------------
# Driver genes aligned to the RIGHT concepts (UNCHANGED)
# -----------------------------
logf("Computing driver genes aligned to canalization/transition (if days identified)...")

DT2 <- data.table::copy(DT)
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
  W <- data.table::dcast(gd_full, Gene ~ Day, value.var=var)
  for (dc in as.character(days_vec)) if (!dc %in% names(W)) W[, (dc) := NA_real_]
  data.table::setcolorder(W, c("Gene", as.character(days_vec)))
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
    data.table::setorder(DTd, -drop_cv_into)
    DTd <- DTd[1:min(200L, .N)]
    out_t6 <- file.path(output_dir, sprintf("Table_6_DriverGenes_Canalization_D%s.csv", canal_day))
    data.table::fwrite(DTd, out_t6)
    register_output(out_t6)
    logf(paste0("Wrote Table_6_DriverGenes_Canalization_D", canal_day, ".csv (canal_day=", canal_day, ")"))
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
    data.table::setorder(DTt, -abs_delta)
    DTt <- DTt[1:min(200L, .N)]
    DTt[, abs_delta := NULL]
    out_t7 <- file.path(output_dir, sprintf("Table_7_DriverGenes_Transition_D%s.csv", infl_day))
    data.table::fwrite(DTt, out_t7)
    register_output(out_t7)
    logf(paste0("Wrote Table_7_DriverGenes_Transition_D", infl_day, ".csv (infl_day=", infl_day, ")"))
  }
}

# -----------------------------
# Plots (UNCHANGED scientific content; routed into outputs/<run_id>/Plots)
#   USER REQUEST (already implemented): no observed series, line through bootstrap medians.
# -----------------------------
logf("Generating annotated plots (bootstrap median ± 95% CI; line through bootstrap medians)...")

ci_dt_plot <- data.table::copy(ci_dt)
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

out_plot <- file.path(output_dir, "Plots", "Metric_Trends_Annotated.png")
ggplot2::ggsave(
  filename = out_plot,
  plot = p1,
  width = 9,
  height = 11,
  dpi = 600
)
register_output(out_plot)
logf("Saved: Plots/Metric_Trends_Annotated.png")

# -----------------------------
# Manifest (MANDATORY): write initial Project_Manifest.json (before render),
# then write inventory after render, then refresh inventory per rule.
# -----------------------------
deps_needed <- c("data.table","matrixStats","changepoint","ggplot2","jsonlite","knitr","rmarkdown")
deps_dt <- deps_table(deps_needed)

run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

manifest <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA, script_full)
  ),
  input = list(
    infile = normalizePath(infile, winslash="/", mustWork=FALSE),
    source = SOURCE,
    days = as.integer(days_vec)
  ),
  parameters = list(
    analysis_name = analysis_name,
    B_TOTAL = as.integer(B_TOTAL),
    USE_LOG1P = isTRUE(USE_LOG1P),
    MIN_NONZERO = as.integer(MIN_NONZERO),
    AMOC_penalty = PENALTY,
    MINSEGLEN = as.integer(MINSEGLEN)
  ),
  dependencies = split(deps_dt, seq_len(nrow(deps_dt))),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash="/", mustWork=FALSE),
    output_dir = normalizePath(output_dir, winslash="/", mustWork=FALSE)
  ),
  generated_files = list() # populated after inventory write
)

write_manifest_json(manifest, manifest_json_path)
logf("Wrote Project_Manifest.json (initial)")

# -----------------------------
# Quarto report (MANDATORY): write <script_name>_Report.qmd and render to HTML as FINAL step
# - Use rmarkdown::render(...)
# -----------------------------
logf("Writing QMD report...")

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))
register_output(qmd_path)
register_output(html_path)

# First inventory write (pre-render) so report can read a file inventory even if render fails
write_files_inventory(output_dir, manifest_files_csv)
register_output(manifest_files_csv)
logf("Wrote Project_Manifest_Files.csv (pre-render)")

# Write QMD (uses inventory file and includes graphics)
# Provide plot relative paths from within output_dir (so include_graphics resolves)
plot_rel <- c(file.path("Plots", "Metric_Trends_Annotated.png"))

write_qmd_report(
  qmd_path = qmd_path,
  script_name = script_name,
  run_id = run_id,
  run_timestamp = run_timestamp,
  script_name_val = script_name,
  script_path_val = script_path,
  script_full_val = script_full,
  input_files = .input_files,
  output_dir = output_dir,
  inventory_csv = manifest_files_csv,
  dependencies_dt = deps_dt,
  plot_relpaths = plot_rel
)
logf(paste0("Wrote QMD: ", qmd_path))

logf("Rendering QMD to HTML via rmarkdown::render (final step)...")
render_res <- render_qmd_to_html(qmd_path, html_path)

render_ok <- is.character(render_res) && length(render_res) == 1L && file.exists(html_path)
if (render_ok) {
  logf(paste0("HTML report created: ", normalizePath(html_path, winslash="/", mustWork=FALSE)))
} else {
  # render_res may be an error object or a path; produce a clear message
  if (inherits(render_res, "error")) {
    logf(paste0("ERROR rendering report: ", conditionMessage(render_res)))
  } else {
    logf("WARNING: Report render did not confirm expected HTML output.")
  }
  logf(paste0("QMD remains at: ", normalizePath(qmd_path, winslash="/", mustWork=FALSE)))
}

# -----------------------------
# File inventory refresh after rendering (MANDATORY) + manifest update
# -----------------------------
write_files_inventory(output_dir, manifest_files_csv)
logf("Refreshed Project_Manifest_Files.csv (post-render)")

inv_dt <- data.table::fread(manifest_files_csv)
manifest$generated_files <- split(inv_dt[, .(relative_path, full_path, bytes, modified)], seq_len(nrow(inv_dt)))

write_manifest_json(manifest, manifest_json_path)
logf("Updated Project_Manifest.json (with generated_files)")

# -----------------------------
# Finish: console summary (MANDATORY includes script_name, script_path, script_full)
# -----------------------------
logf("DONE")

cat("\n==============================================\n")
cat("COMPLETE\n")
cat("Run ID:\n  ", run_id, "\n", sep="")
cat("Outputs written under:\n  ", normalizePath(output_dir, winslash="/", mustWork=FALSE), "\n", sep="")

cat("\nScript identity:\n")
cat("  script_name: ", script_name, "\n", sep="")
cat("  script_path: ", script_path, "\n", sep="")
cat("  script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep="")
if (is.na(script_full)) {
  cat("  NOTE: script_full could not be detected; fallback filename was used.\n")
}

cat("\nKey outputs:\n")
cat("  - Table_1_Daily_Metrics.csv\n")
cat("  - Table_1b_Daily_Metrics_BootstrapCI.csv\n")
cat("  - Table_2_Canalization_Minima_Bootstrap.csv\n")
cat("  - Table_3_Piecewise_Inflection_Bootstrap.csv\n")
cat("  - Table_4_Excursion_Day_Bootstrap.csv\n")
cat("  - Table_5_AMOC_TwoRegimeSplit_Bootstrap.csv\n")
cat("  - Project_Manifest.json\n")
cat("  - Project_Manifest_Files.csv\n")
cat("Plots:\n")
cat("  - Plots/Metric_Trends_Annotated.png\n")
cat("Log:\n")
cat("  - ", normalizePath(logfile, winslash="/", mustWork=FALSE), "\n", sep="")

if (file.exists(html_path)) {
  cat("\nHTML report created:\n  ", normalizePath(html_path, winslash="/", mustWork=FALSE), "\n", sep="")
} else {
  cat("\nHTML report creation FAILED.\n")
  cat("QMD is available at:\n  ", normalizePath(qmd_path, winslash="/", mustWork=FALSE), "\n", sep="")
}
cat("==============================================\n\n")
