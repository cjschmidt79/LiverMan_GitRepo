#!/usr/bin/env Rscript
# ============================================================
# Mean + Dispersion + Canalization Highlight Figure Builder
#
# INPUT (tidy long CSV; one row per measurement):
#   Day | SampleID | Feature | Abundance
# Optional:
#   Source
#
# OUTPUT (always under ./outputs/<run_id>/):
#   - Multi-panel figure PNG (+ optional PDF)
#   - Daily summary CSV
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - <script_name>_Report.qmd + rendered HTML
#
# PANELS (default):
#   A: System mean across features (per Day) with error bars across SampleID replicate-means
#   B: System dispersion across features (per Day) with interval (bootstrap CI or MAD)
#   C: Landmark day of maximal dispersion (from Panel B)
#   D: CV vs Dispersion “quadrant” plot (per Day), highlights scaling vs absolute variability
#   E: Mean–dispersion decoupling index Δ(d) = Z_mean(d) - Z_disp(d) (single curve)
#
# USAGE (CLI):
#   Rscript make_mean_dispersion_canalization_figure.R --in="path/to/tidy.csv"
#
# KEY OPTIONS:
#   --log_abundance=0/1
#   --split_by_source=0/1 (default auto if Source has >1 level)
#   --include_pooled=0/1
#
# Panel A interval:
#   --a_interval=mad|ci|se|none   (default: mad)
#   --a_center=mean|median        (default: mean)
#   --mad_k=1                     (MAD scaling; use 1.4826 for SD-equivalent)
#
# Panel B interval:
#   --b_interval=ci|mad|none      (default: ci)
#   --disp_metric=sd|cv|iqr|mad   (per feature within day; default: sd)
#   --disp_agg=median|mean        (aggregate across features; default: median)
#   --ci=0/1                      (default: 1)
#   --ci_level=0.95
#   --ci_mode=sample|feature
#   --n_boot=500
#   --ask_n_boot=0/1              (force interactive prompt when CI is on)
#
# Panels toggle:
#   --panel_C=0/1  --panel_D=0/1  --panel_E=0/1
#
# ============================================================

# --------------------------
# Script identity capture (MANDATORY BLOCK; do not modify structure)
# --------------------------
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
known_script_filename <- "make_mean_dispersion_canalization_figure.R"
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

# --------------------------
# Dependency handling (install if missing; do not uninstall)
# --------------------------
quiet_install_if_missing <- function(pkgs, repos = getOption("repos")) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("[deps] Installing missing package: ", p)
      install.packages(p, repos = repos, quiet = TRUE)
    }
  }
}

quiet_install_if_missing(c("jsonlite", "knitr", "rmarkdown"))
quiet_install_if_missing(c("data.table", "ggplot2", "patchwork"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

stop2 <- function(...) stop(paste0(...), call. = FALSE)
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

safe_slug <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "NA", x)
}

# --------------------------
# Lightweight CLI parsing
# --------------------------
args <- commandArgs(trailingOnly = TRUE)

arg_val <- function(flag, default = NULL) {
  hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", flag, "="), "", hit[1])
}
arg_flag <- function(flag, default = FALSE) {
  v <- arg_val(flag, NULL)
  if (is.null(v)) return(default)
  v <- tolower(trimws(v))
  v %in% c("1", "true", "t", "yes", "y")
}

# --------------------------
# Console identity summary (mandatory)
# --------------------------
message("------------------------------------------------------------")
message("Script identity")
message("  script_name: ", script_name)
message("  script_path: ", script_path)
message("  script_full: ", ifelse(is.na(script_full), "NA", script_full))
if (is.na(script_full)) {
  message("  NOTE: script path detection failed in this execution mode; using fallback filename: ", known_script_filename)
}
message("------------------------------------------------------------")

# --------------------------
# Inputs + parameters (outputs policy enforced)
# --------------------------
in_path <- arg_val("in", NULL)

analysis_name <- arg_val("analysis_name", "MeanDispersionCanalization")
run_source_tag <- arg_val("run_source", "")

log_abundance <- arg_flag("log_abundance", FALSE)
write_pdf <- arg_flag("pdf", FALSE)

# Source handling
split_by_source_raw <- arg_val("split_by_source", NULL) # NULL => auto
include_pooled <- arg_flag("include_pooled", FALSE)
source_regex <- arg_val("source_regex", "")

# Panel toggles
panel_C <- arg_flag("panel_C", TRUE)
panel_D <- arg_flag("panel_D", TRUE)
panel_E <- arg_flag("panel_E", TRUE)

# Panel A options
a_interval <- tolower(arg_val("a_interval", "mad"))   # mad|ci|se|none
a_center   <- tolower(arg_val("a_center", "mean"))    # mean|median
mad_k      <- as.numeric(arg_val("mad_k", "1"))

# Panel B options
disp_metric <- tolower(arg_val("disp_metric", "sd"))  # sd|cv|iqr|mad
disp_agg    <- tolower(arg_val("disp_agg", "median")) # median|mean
b_interval  <- tolower(arg_val("b_interval", "ci"))   # ci|mad|none

# CI options (Panel B; also usable for Panel A if a_interval=ci)
ci_on      <- arg_flag("ci", TRUE)
ci_level   <- as.numeric(arg_val("ci_level", "0.95"))
ci_mode    <- tolower(arg_val("ci_mode", "sample"))   # sample|feature
seed       <- as.integer(arg_val("seed", "1"))
ask_n_boot <- arg_flag("ask_n_boot", FALSE)

n_boot_raw <- arg_val("n_boot", NULL)
suggested_default_boot <- 500L

# Column overrides
col_day <- arg_val("col_day", NULL)
col_sample <- arg_val("col_sample", NULL)
col_feature <- arg_val("col_feature", NULL)
col_abundance <- arg_val("col_abundance", NULL)
col_source <- arg_val("col_source", NULL)

# Validate
if (!(ci_level > 0 && ci_level < 1)) stop2("--ci_level must be between 0 and 1.")
if (!(ci_mode %in% c("sample", "feature"))) stop2("--ci_mode must be 'sample' or 'feature'.")
if (!(a_interval %in% c("mad", "ci", "se", "none"))) stop2("--a_interval must be mad|ci|se|none.")
if (!(a_center %in% c("mean", "median"))) stop2("--a_center must be mean|median.")
if (!(b_interval %in% c("ci", "mad", "none"))) stop2("--b_interval must be ci|mad|none.")
if (!(disp_metric %in% c("sd", "cv", "iqr", "mad"))) stop2("--disp_metric must be sd|cv|iqr|mad.")
if (!(disp_agg %in% c("median", "mean"))) stop2("--disp_agg must be median|mean.")
if (!is.finite(mad_k) || is.na(mad_k) || mad_k <= 0) stop2("--mad_k must be a positive number.")

set.seed(seed)

# --------------------------
# Resolve n_boot (prompt logic)
# --------------------------
resolve_n_boot <- function(current_value) {
  if (!interactive()) return(current_value)
  prompt <- paste0("How many bootstrap replicates? [default ", suggested_default_boot,
                   if (!is.null(current_value)) paste0(", current ", current_value) else "",
                   "]: ")
  ans <- trimws(readline(prompt))
  if (nchar(ans) == 0) {
    if (!is.null(current_value)) return(current_value)
    return(suggested_default_boot)
  }
  nb <- suppressWarnings(as.integer(ans))
  if (!is.finite(nb) || is.na(nb) || nb < 10) {
    message("Invalid entry; using ", suggested_default_boot)
    return(suggested_default_boot)
  }
  nb
}

if (ci_on && (b_interval == "ci" || a_interval == "ci")) {
  if (ask_n_boot && interactive()) {
    n_boot <- resolve_n_boot(if (!is.null(n_boot_raw) && nchar(n_boot_raw) > 0) as.integer(n_boot_raw) else NULL)
  } else if (is.null(n_boot_raw) || nchar(n_boot_raw) == 0) {
    if (interactive()) {
      n_boot <- resolve_n_boot(NULL)
    } else {
      n_boot <- suggested_default_boot
      message("n_boot not provided; using default n_boot = ", n_boot)
    }
  } else {
    n_boot <- suppressWarnings(as.integer(n_boot_raw))
    if (!is.finite(n_boot) || is.na(n_boot) || n_boot < 10) stop2("--n_boot must be an integer >= 10.")
  }
} else {
  n_boot <- 0L
}

# --------------------------
# Input selection (interactive-friendly)
# --------------------------
if (is.null(in_path) || nchar(in_path) == 0) {
  if (interactive() && exists("file.choose")) {
    message("\nSelect a tidy-long CSV with columns like:")
    message("  Day | SampleID | Feature | Abundance   (optional: Source)")
    message("Tip: One row = one Feature measurement in one Sample on one Day.\n")
    in_path <- file.choose()
  } else {
    stop2("No input provided. Run:\n  Rscript ", known_script_filename, " --in=\"path/to/tidy.csv\"")
  }
}

message("Input:  ", in_path)
message("Panels: A + B",
        if (panel_C) " + C" else "",
        if (panel_D) " + D" else "",
        if (panel_E) " + E" else "")
message("Panel A interval: ", a_interval, " | center=", a_center)
message("Panel B interval: ", b_interval, " | disp_metric=", disp_metric, " | disp_agg=", disp_agg)
message("CI: ", if (ci_on) paste0("ON (", round(ci_level * 100), "%, mode=", ci_mode, ", n_boot=", n_boot, ")") else "OFF")

# --------------------------
# Output directory policy (MANDATORY)
# --------------------------
outputs_root <- file.path(getwd(), "outputs")
dir.create(outputs_root, showWarnings = FALSE, recursive = TRUE)

# --------------------------
# Read data
# --------------------------
dt0 <- fread(in_path, showProgress = FALSE)
if (nrow(dt0) == 0) stop2("Input file has 0 rows.")
setDT(dt0)

# --------------------------
# Column detection helpers
# --------------------------
nm0 <- names(dt0)
nm <- tolower(nm0)

pick_col <- function(preferred, patterns, label) {
  if (!is.null(preferred)) {
    if (!(preferred %in% nm0)) stop2("Provided ", label, " override not found: ", preferred)
    return(preferred)
  }
  hits <- which(Reduce(`|`, lapply(patterns, function(p) grepl(p, nm))))
  if (length(hits) == 1) return(nm0[hits])
  if (length(hits) == 0) return(NULL)
  nm0[hits]
}
is_multi <- function(x) is.character(x) && length(x) > 1

cand_day <- pick_col(col_day, c("^day$", "day", "timepoint", "age"), "Day column")
cand_sample <- pick_col(col_sample, c("^sampleid$", "^sample$", "sample", "replicate", "bird", "id$"), "SampleID column")
cand_feature <- pick_col(col_feature, c("^feature$", "^gene$", "^metabolite$", "gene", "metab", "compound", "feature"), "Feature column")
cand_abund <- pick_col(col_abundance, c("^abundance$", "^value$", "^expr$", "abund", "count", "tpm", "fpkm", "intensity", "area", "value"), "Abundance column")
cand_source <- pick_col(col_source, c("^source$", "tissue", "compartment"), "Source column")

if (is_multi(cand_day) || is_multi(cand_sample) || is_multi(cand_feature) || is_multi(cand_abund)) {
  message("\nColumn auto-detection is ambiguous. Candidates:")
  if (is_multi(cand_day))     message("  Day candidates:       ", paste(cand_day, collapse = ", "))
  if (is_multi(cand_sample))  message("  Sample candidates:    ", paste(cand_sample, collapse = ", "))
  if (is_multi(cand_feature)) message("  Feature candidates:   ", paste(cand_feature, collapse = ", "))
  if (is_multi(cand_abund))   message("  Abundance candidates: ", paste(cand_abund, collapse = ", "))
  message("\nRe-run with explicit overrides, e.g.:")
  message("  Rscript ", known_script_filename, " --in=\"...\" --col_day=\"Day\" --col_sample=\"SampleID\" --col_feature=\"Feature\" --col_abundance=\"Abundance\"")
  stop2("Stopping due to ambiguous column detection.")
}

if (is.null(cand_day) || is.null(cand_sample) || is.null(cand_feature) || is.null(cand_abund)) {
  message("\nCould not auto-detect required columns.")
  message("Available columns: ", paste(nm0, collapse = ", "))
  message("\nRe-run with overrides:")
  message("  --col_day= --col_sample= --col_feature= --col_abundance=")
  stop2("Stopping due to missing required columns.")
}

col_day <- cand_day
col_sample <- cand_sample
col_feature <- cand_feature
col_abundance <- cand_abund
if (!is_multi(cand_source) && !is.null(cand_source)) col_source <- cand_source

message("\nUsing columns:")
message("  Day:       ", col_day)
message("  SampleID:  ", col_sample)
message("  Feature:   ", col_feature)
message("  Abundance: ", col_abundance)
if (!is.null(col_source)) message("  Source:    ", col_source)

# --------------------------
# Standardize columns
# --------------------------
dt0[, Day := suppressWarnings(as.numeric(get(col_day)))]
dt0[, SampleID := as.character(get(col_sample))]
dt0[, Feature := as.character(get(col_feature))]
dt0[, Abundance := suppressWarnings(as.numeric(get(col_abundance)))]

if (!is.null(col_source)) {
  dt0[, Source := as.character(get(col_source))]
} else {
  dt0[, Source := NA_character_]
}

dt0 <- dt0[is.finite(Day) & is.finite(Abundance) & !is.na(SampleID) & !is.na(Feature)]
if (nrow(dt0) == 0) stop2("After filtering non-finite values, 0 rows remain.")

if (log_abundance) {
  dt0[, Abundance := log1p(Abundance)]
  message("Applied log1p transform to Abundance.")
}

if (!is.null(col_source) && nchar(source_regex) > 0) {
  keep <- grepl(source_regex, dt0$Source, ignore.case = TRUE)
  if (!any(keep)) stop2("source_regex='", source_regex, "' matched 0 rows.")
  dt0 <- dt0[keep]
  message("Applied source_regex filter; rows now: ", nrow(dt0))
}

has_source <- !is.null(col_source) && !all(is.na(dt0$Source))
n_sources <- if (has_source) length(unique(dt0[!is.na(Source), Source])) else 0L

if (is.null(split_by_source_raw)) {
  split_by_source <- has_source && n_sources > 1
} else {
  split_by_source <- tolower(trimws(split_by_source_raw)) %in% c("1", "true", "t", "yes", "y")
}

message("Source handling: ",
        if (!has_source) "no Source column detected"
        else paste0(n_sources, " source level(s) detected; split_by_source=", split_by_source,
                    if (include_pooled) " (also pooled)" else ""))

# --------------------------
# Run folder creation (MANDATORY)
#   run_id = <analysis_name>_<Source>_<YYYYMMDD_HHMMSS>
# --------------------------
infer_run_source <- function() {
  if (nzchar(run_source_tag)) return(run_source_tag)
  if (has_source && split_by_source && n_sources > 1) return("MultiSource")
  if (has_source && n_sources == 1) {
    s <- unique(dt0[!is.na(Source), Source])[1]
    return(safe_slug(s))
  }
  "POOLED"
}

run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_id <- paste0(safe_slug(analysis_name), "_", safe_slug(infer_run_source()), "_", ts_stamp())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Outputs root: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE))
message("Run folder:   ", normalizePath(output_dir,  winslash = "/", mustWork = FALSE))

# --------------------------
# Manifest scaffolding (mandatory)
# --------------------------
get_deps <- function(pkgs) {
  lapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) {
      list(package = p, version = as.character(utils::packageVersion(p)))
    } else {
      list(package = p, version = NA_character_)
    }
  })
}

manifest_path_json <- file.path(output_dir, "Project_Manifest.json")
manifest_path_files <- file.path(output_dir, "Project_Manifest_Files.csv")

write_manifest_json <- function(manifest) {
  jsonlite::write_json(manifest, manifest_path_json, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

refresh_file_inventory <- function() {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (length(files) == 0) {
    inv <- data.table(file = character(), relative_path = character(), bytes = numeric(), modified_time = character())
  } else {
    out_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
    files_norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
    rel <- sub(paste0("^", gsub("([\\^\\$\\.|\\+\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", out_norm), "/?"), "", files_norm)
    finfo <- file.info(files)
    inv <- data.table(
      file = basename(files),
      relative_path = rel,
      bytes = as.numeric(finfo$size),
      modified_time = format(finfo$mtime, "%Y-%m-%d %H:%M:%S")
    )[order(relative_path)]
  }
  fwrite(inv, manifest_path_files)
  inv
}

required_pkgs <- c("jsonlite", "knitr", "rmarkdown", "data.table", "ggplot2", "patchwork")

manifest <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA_character_, script_full)
  ),
  input = list(
    in_path = normalizePath(in_path, winslash = "/", mustWork = FALSE),
    detected_columns = list(
      Day = col_day,
      SampleID = col_sample,
      Feature = col_feature,
      Abundance = col_abundance,
      Source = if (!is.null(col_source)) col_source else NA_character_
    )
  ),
  parameters = list(
    analysis_name = analysis_name,
    log_abundance = log_abundance,
    split_by_source = split_by_source,
    include_pooled = include_pooled,
    source_regex = source_regex,
    panels = list(C = panel_C, D = panel_D, E = panel_E),
    a_interval = a_interval,
    a_center = a_center,
    b_interval = b_interval,
    mad_k = mad_k,
    disp_metric = disp_metric,
    disp_agg = disp_agg,
    ci = ci_on,
    ci_level = ci_level,
    ci_mode = ci_mode,
    n_boot = n_boot,
    ask_n_boot = ask_n_boot,
    seed = seed,
    pdf = write_pdf
  ),
  dependencies = get_deps(required_pkgs),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  ),
  generated_files = list()
)

write_manifest_json(manifest)
refresh_file_inventory()

# --------------------------
# Core computation helpers
# --------------------------
disp_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  if (disp_metric == "sd")  return(sd(x))
  if (disp_metric == "cv")  { m <- mean(x); if (!is.finite(m) || m == 0) return(NA_real_); return(sd(x) / abs(m)) }
  if (disp_metric == "iqr") return(IQR(x))
  if (disp_metric == "mad") return(mad(x, constant = 1))
  stop2("Unknown --disp_metric: ", disp_metric)
}

agg_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (disp_agg == "median") return(median(x))
  if (disp_agg == "mean")   return(mean(x))
  stop2("Unknown --disp_agg: ", disp_agg)
}

robust_z <- function(x) {
  x <- as.numeric(x)
  m <- median(x, na.rm = TRUE)
  s <- mad(x, constant = 1, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

# Extract initial comment/header block (for QMD 4607).
extract_initial_comment_block <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(lines) || length(lines) == 0) return(NULL)
  
  out <- character()
  i <- 1L
  if (grepl("^#!", lines[1])) {
    out <- c(out, lines[1])
    i <- 2L
  }
  while (i <= length(lines)) {
    ln <- lines[i]
    if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
      out <- c(out, ln)
      i <- i + 1L
    } else {
      break
    }
  }
  if (length(out) == 0) return(NULL)
  out
}

# --------------------------
# Bootstrap CI helpers
# --------------------------
bootstrap_ci_over_samples <- function(dt_day, stat_fun, n_boot, q_lo, q_hi) {
  # dt_day has columns: SampleID, Feature, Abundance (one Day)
  sids <- unique(dt_day$SampleID)
  if (length(sids) < 2) return(c(NA_real_, NA_real_))
  boot_vals <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    sids_b <- sample(sids, size = length(sids), replace = TRUE)
    dboot <- rbindlist(lapply(sids_b, function(id) dt_day[SampleID == id]), use.names = TRUE)
    boot_vals[b] <- stat_fun(dboot)
  }
  c(
    as.numeric(quantile(boot_vals, q_lo, na.rm = TRUE)),
    as.numeric(quantile(boot_vals, q_hi, na.rm = TRUE))
  )
}

bootstrap_ci_over_features <- function(gene_disp_day, agg_fun_local, n_boot, q_lo, q_hi) {
  # gene_disp_day has column feat_disp (one Day)
  vals <- gene_disp_day$feat_disp
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2) return(c(NA_real_, NA_real_))
  nF <- length(vals)
  boot_vals <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    samp <- sample(vals, size = nF, replace = TRUE)
    boot_vals[b] <- agg_fun_local(samp)
  }
  c(
    as.numeric(quantile(boot_vals, q_lo, na.rm = TRUE)),
    as.numeric(quantile(boot_vals, q_hi, na.rm = TRUE))
  )
}

# --------------------------
# Main per-source computation
# --------------------------
compute_panels <- function(dt, source_label, out_dir, file_tag) {
  # dt must have: Day, SampleID, Feature, Abundance
  
  # ---- Panel A: replicate means across features ----
  rep_means <- dt[, .(rep_mean = mean(Abundance, na.rm = TRUE)), by = .(Day, SampleID)]
  
  if (a_center == "mean") {
    A_df <- rep_means[, .(
      center = mean(rep_mean, na.rm = TRUE),
      n_rep = .N
    ), by = Day][order(Day)]
  } else {
    A_df <- rep_means[, .(
      center = median(rep_mean, na.rm = TRUE),
      n_rep = .N
    ), by = Day][order(Day)]
  }
  
  A_df[, `:=`(lo = NA_real_, hi = NA_real_, interval_type = "none")]
  
  if (a_interval == "se") {
    A_se <- rep_means[, .(
      se = sd(rep_mean, na.rm = TRUE) / sqrt(.N)
    ), by = Day]
    setkey(A_df, Day); setkey(A_se, Day)
    A_df[A_se, `:=`(lo = center - i.se, hi = center + i.se, interval_type = "SE")]
  } else if (a_interval == "mad") {
    A_mad <- rep_means[, .(
      m = mad(rep_mean, constant = 1, na.rm = TRUE)
    ), by = Day]
    setkey(A_df, Day); setkey(A_mad, Day)
    A_df[A_mad, `:=`(
      lo = center - mad_k * i.m,
      hi = center + mad_k * i.m,
      interval_type = paste0("MAD (k=", mad_k, ")")
    )]
  } else if (a_interval == "ci" && ci_on) {
    alpha <- (1 - ci_level) / 2
    q_lo <- alpha
    q_hi <- 1 - alpha
    days <- sort(unique(rep_means$Day))
    ci_list <- vector("list", length(days))
    
    # Bootstrap over SampleID within day on rep_means (stat = mean or median)
    for (i in seq_along(days)) {
      d <- days[i]
      x <- rep_means[Day == d]$rep_mean
      x <- x[is.finite(x)]
      if (length(x) < 2) {
        ci_list[[i]] <- data.table(Day = d, lo = NA_real_, hi = NA_real_)
        next
      }
      boot_vals <- numeric(n_boot)
      nS <- length(x)
      for (b in seq_len(n_boot)) {
        samp <- sample(x, size = nS, replace = TRUE)
        boot_vals[b] <- if (a_center == "mean") mean(samp) else median(samp)
      }
      ci_list[[i]] <- data.table(
        Day = d,
        lo = as.numeric(quantile(boot_vals, q_lo, na.rm = TRUE)),
        hi = as.numeric(quantile(boot_vals, q_hi, na.rm = TRUE))
      )
    }
    A_ci <- rbindlist(ci_list)
    setkey(A_df, Day); setkey(A_ci, Day)
    A_df[A_ci, `:=`(lo = i.lo, hi = i.hi, interval_type = paste0(round(ci_level * 100), "% bootstrap CI (n_boot=", n_boot, ")"))]
  }
  
  # ---- Panel B: feature-wise dispersion within day, aggregated across features ----
  gene_disp <- dt[, .(feat_disp = disp_fun(Abundance)), by = .(Day, Feature)]
  gene_disp <- gene_disp[is.finite(feat_disp)]
  if (nrow(gene_disp) == 0) stop2("No feature-level dispersion values computed for Source=", source_label)
  
  # Also compute system-level CV (median feature-wise CV) for Panel D
  gene_cv <- dt[, .(feat_cv = {
    x <- Abundance[is.finite(Abundance)]
    if (length(x) < 2) NA_real_ else {
      m <- mean(x)
      if (!is.finite(m) || m == 0) NA_real_ else sd(x) / abs(m)
    }
  }), by = .(Day, Feature)]
  gene_cv <- gene_cv[is.finite(feat_cv)]
  
  B_df <- gene_disp[, .(
    disp = agg_fun(feat_disp),
    n_feat = .N,
    mad_featdisp = mad(feat_disp, constant = 1, na.rm = TRUE)
  ), by = Day][order(Day)]
  
  B_df[, `:=`(lo = NA_real_, hi = NA_real_, interval_type = "none")]
  
  if (b_interval == "mad") {
    B_df[, `:=`(
      lo = disp - mad_k * mad_featdisp,
      hi = disp + mad_k * mad_featdisp,
      interval_type = paste0("MAD (k=", mad_k, ")")
    )]
  } else if (b_interval == "ci" && ci_on) {
    alpha <- (1 - ci_level) / 2
    q_lo <- alpha
    q_hi <- 1 - alpha
    days <- sort(unique(dt$Day))
    ci_list <- vector("list", length(days))
    
    if (ci_mode == "feature") {
      for (i in seq_along(days)) {
        d <- days[i]
        gdd <- gene_disp[Day == d]
        ci <- bootstrap_ci_over_features(gdd, agg_fun, n_boot, q_lo, q_hi)
        ci_list[[i]] <- data.table(Day = d, lo = ci[1], hi = ci[2])
      }
    } else {
      # ci_mode == sample
      for (i in seq_along(days)) {
        d <- days[i]
        ddt <- dt[Day == d]
        stat_fun <- function(dboot) {
          gd_b <- dboot[, .(feat_disp = disp_fun(Abundance)), by = Feature]
          agg_fun(gd_b$feat_disp)
        }
        ci <- bootstrap_ci_over_samples(ddt, stat_fun, n_boot, q_lo, q_hi)
        ci_list[[i]] <- data.table(Day = d, lo = ci[1], hi = ci[2])
      }
    }
    
    B_ci <- rbindlist(ci_list)
    setkey(B_df, Day); setkey(B_ci, Day)
    B_df[B_ci, `:=`(lo = i.lo, hi = i.hi,
                    interval_type = paste0(round(ci_level * 100), "% bootstrap CI (mode=", ci_mode, ", n_boot=", n_boot, ")"))]
  }
  
  # ---- Panel C landmark (max dispersion day) ----
  landmark_day <- B_df[which.max(disp), Day]
  landmark_disp <- B_df[which.max(disp), disp]
  
  # ---- Panel D: CV vs Dispersion quadrant ----
  CV_df <- gene_cv[, .(sys_cv = median(feat_cv, na.rm = TRUE)), by = Day][order(Day)]
  setkey(CV_df, Day); setkey(B_df, Day)
  D_df <- merge(B_df[, .(Day, disp)], CV_df, by = "Day", all = TRUE)
  
  # ---- Panel E: Mean–dispersion decoupling index Δ(d) = Z_mean(d) - Z_disp(d) ----
  # Use Panel A mean center as the "mean" curve; use Panel B disp as "disp" curve.
  # Robust z-scoring across days within this source.
  dec_df <- merge(A_df[, .(Day, mean_center = center)], B_df[, .(Day, disp)], by = "Day", all = TRUE)
  dec_df[, z_mean := robust_z(mean_center)]
  dec_df[, z_disp := robust_z(disp)]
  dec_df[, delta := z_mean - z_disp]
  
  # --------------------------
  # Plot A: mean with error bars (NO ribbon)
  # --------------------------
  A_interval_on <- any(is.finite(A_df$lo)) && any(is.finite(A_df$hi)) && a_interval != "none"
  pA <- ggplot(A_df, aes(x = Day, y = center)) +
    { if (A_interval_on) geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) } +
    geom_line() +
    geom_point() +
    labs(
      title = "A",
      x = "Day",
      y = "System mean (across features)",
      subtitle = paste0("Interval: ", if (A_interval_on) unique(A_df$interval_type)[1] else "none",
                        " | Source=", source_label)
    ) +
    theme_classic()
  
  # --------------------------
  # Plot B: dispersion with error bars (NO ribbon)
  # --------------------------
  B_interval_on <- any(is.finite(B_df$lo)) && any(is.finite(B_df$hi)) && b_interval != "none"
  pB <- ggplot(B_df, aes(x = Day, y = disp)) +
    { if (B_interval_on) geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) } +
    geom_line() +
    geom_point() +
    labs(
      title = "B",
      x = "Day",
      y = paste0(disp_agg, " feature-wise ", disp_metric, " (within day)"),
      subtitle = paste0("Interval: ", if (B_interval_on) unique(B_df$interval_type)[1] else "none",
                        " | Source=", source_label)
    ) +
    theme_classic()
  
  # --------------------------
  # Plot C: landmark day of maximal dispersion
  # --------------------------
  pC <- ggplot(B_df, aes(x = Day, y = disp)) +
    geom_line() +
    geom_point() +
    geom_point(data = B_df[Day == landmark_day], size = 3) +
    labs(
      title = "C",
      x = "Day",
      y = "Dispersion (Panel B)",
      subtitle = paste0("Landmark day of maximal dispersion: Day ", landmark_day,
                        " (disp=", signif(landmark_disp, 4), ") | Source=", source_label)
    ) +
    theme_classic()
  
  # --------------------------
  # Plot D: CV vs Dispersion quadrant
  # --------------------------
  pD <- ggplot(D_df, aes(x = sys_cv, y = disp)) +
    geom_point() +
    geom_text(aes(label = Day), vjust = -0.8, check_overlap = TRUE) +
    geom_point(data = D_df[Day == landmark_day], size = 3) +
    labs(
      title = "D",
      x = "System CV (median feature-wise CV)",
      y = "System dispersion (Panel B)",
      subtitle = paste0("Quadrant view: separates absolute variability (disp) from scaling effects (CV). ",
                        "Highlighted point = Day ", landmark_day, " | Source=", source_label)
    ) +
    theme_classic()
  
  # --------------------------
  # Plot E: mean–dispersion decoupling Δ(d)
  # --------------------------
  pE <- ggplot(dec_df, aes(x = Day, y = delta)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line() +
    geom_point() +
    geom_point(data = dec_df[Day == landmark_day], size = 3) +
    labs(
      title = "E",
      x = "Day",
      y = expression(Delta(d) == Z[mean](d) - Z[disp](d)),
      subtitle = paste0("Decoupling index: positive => mean high relative to dispersion (canalization-like regime). ",
                        "Highlighted day = Day ", landmark_day, " | Source=", source_label)
    ) +
    theme_classic()
  
  # --------------------------
  # Assemble figure
  # --------------------------
  fig <- pA / pB
  if (panel_C) fig <- fig / pC
  if (panel_D) fig <- fig / pD
  if (panel_E) fig <- fig / pE
  
  # ---- Write outputs ----
  base <- paste0("Figure_mean_dispersion_canalization_", safe_slug(file_tag), "_", ts_stamp())
  out_png <- file.path(out_dir, paste0(base, ".png"))
  out_pdf <- file.path(out_dir, paste0(base, ".pdf"))
  out_csv <- file.path(out_dir, paste0(base, "_daily_summary.csv"))
  
  daily_summary <- merge(
    merge(A_df, B_df, by = "Day", all = TRUE, suffixes = c(".A", ".B")),
    D_df, by = "Day", all = TRUE
  )
  daily_summary[, `:=`(
    Source = source_label,
    Landmark_Day = landmark_day,
    Decoupling_Delta = dec_df$delta[match(Day, dec_df$Day)]
  )]
  
  fwrite(daily_summary, out_csv)
  ggsave(out_png, fig,
         width = 7.5,
         height = 7.5 + 2.5 * (panel_C + panel_D + panel_E),
         dpi = 300)
  message("Wrote: ", out_png)
  message("Wrote: ", out_csv)
  
  if (write_pdf) {
    ggsave(out_pdf, fig,
           width = 7.5,
           height = 7.5 + 2.5 * (panel_C + panel_D + panel_E))
    message("Wrote: ", out_pdf)
  }
  
  invisible(list(
    A_df = A_df, B_df = B_df, D_df = D_df, dec_df = dec_df,
    png = out_png, csv = out_csv, pdf = if (write_pdf) out_pdf else NA_character_,
    landmark_day = landmark_day
  ))
}

# --------------------------
# Run analyses
# --------------------------
results <- list()

if (split_by_source && has_source) {
  src_vals <- sort(unique(dt0[!is.na(Source), Source]))
  for (s in src_vals) {
    dts <- dt0[Source == s]
    if (nrow(dts) == 0) next
    message("\n--- Analyzing Source: ", s, " (rows=", nrow(dts), ") ---")
    results[[as.character(s)]] <- compute_panels(
      dts,
      source_label = s,
      out_dir = output_dir,
      file_tag = paste0("Source_", s)
    )
  }
  
  if (include_pooled) {
    message("\n--- Analyzing POOLED across sources (rows=", nrow(dt0), ") ---")
    results[["POOLED"]] <- compute_panels(
      dt0,
      source_label = "POOLED",
      out_dir = output_dir,
      file_tag = "POOLED"
    )
  }
} else {
  source_label <- if (has_source && n_sources == 1) unique(dt0[!is.na(Source), Source])[1] else "POOLED"
  message("\n--- Analyzing dataset (split_by_source=FALSE); Source label: ", source_label, " ---")
  results[[source_label]] <- compute_panels(
    dt0,
    source_label = source_label,
    out_dir = output_dir,
    file_tag = source_label
  )
}

# --------------------------
# Refresh manifest + inventory after analysis outputs exist
# --------------------------
inv_mid <- refresh_file_inventory()
manifest$generated_files <- list(
  inventory_csv = normalizePath(manifest_path_files, winslash = "/", mustWork = FALSE),
  n_files = nrow(inv_mid)
)
write_manifest_json(manifest)

# --------------------------
# Quarto QMD report (MANDATORY) + render to HTML as FINAL step
# --------------------------
header_lines <- extract_initial_comment_block(script_full)
if (is.null(header_lines)) {
  header_lines <- c(
    "# (Header block unavailable: script_full could not be resolved in this execution mode.)",
    "# See console identity summary and manifest for script identity and parameters."
  )
}
SCRIPT_HEADER_TEXT <- paste(header_lines, collapse = "\n")

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_name <- paste0(script_name, "_Report.html")
html_path <- file.path(output_dir, html_name)

# Choose a representative always-on plot: first PNG in output_dir
png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
plot_to_include <- if (length(png_files) > 0) normalizePath(png_files[1], winslash = "/", mustWork = FALSE) else NA_character_

manifest$generated_files$preferred_plot <- if (!is.na(plot_to_include)) plot_to_include else NA_character_
write_manifest_json(manifest)

qmd_lines <- c(
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
  "## Run metadata",
  "",
  paste0("- **run_id:** ", run_id),
  paste0("- **run_timestamp:** ", run_timestamp),
  paste0("- **output_dir:** ", normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
  "",
  "### Script identity",
  "",
  paste0("- **script_name:** ", script_name),
  paste0("- **script_path:** ", script_path),
  paste0("- **script_full:** ", ifelse(is.na(script_full), "NA", script_full)),
  "",
  "### Input",
  "",
  paste0("- **input_csv:** ", normalizePath(in_path, winslash = "/", mustWork = FALSE)),
  "",
  "## Key parameters",
  "",
  "```{r}",
  "library(jsonlite)",
  "library(knitr)",
  "man <- jsonlite::fromJSON(file.path('.', 'Project_Manifest.json'), simplifyVector = TRUE)",
  "params <- man$parameters",
  "params_flat <- unlist(params, recursive = TRUE, use.names = TRUE)",
  "params_df <- data.frame(Parameter = names(params_flat), Value = as.character(params_flat), stringsAsFactors = FALSE)",
  "knitr::kable(params_df, align = 'l')",
  "```",
  "",
  "## Script header",
  "",
  "```",
  SCRIPT_HEADER_TEXT,
  "```",
  "",
  "## Dependencies",
  "",
  "```{r}",
  "deps <- man$dependencies",
  "deps_df <- if (length(deps) == 0) data.frame() else do.call(rbind, lapply(deps, as.data.frame))",
  "knitr::kable(deps_df, align = 'l')",
  "```",
  "",
  "## Generated files",
  "",
  "```{r}",
  "library(data.table)",
  "inv <- data.table::fread(file.path('.', 'Project_Manifest_Files.csv'))",
  "knitr::kable(inv, align = 'l')",
  "```",
  "",
  "## Analytical logic and formulas",
  "",
  "This analysis produces a multi-panel, reader-facing summary from tidy-long data (Day, SampleID, Feature, Abundance).",
  "",
  "### Panel A: system mean across features with uncertainty across SampleIDs",
  "",
  "For each Day \\(d\\) and SampleID \\(s\\), compute the replicate mean across features:",
  "",
  "\\[ \\bar{x}_{d,s} = \\frac{1}{F} \\sum_{f=1}^{F} x_{d,s,f} \\]",
  "",
  "Then compute the Day-level center (mean or median across \\(\\bar{x}_{d,s}\\)) and an interval across SampleIDs using one of:",
  "",
  "- **SE**: \\( \\mathrm{SE}(\\bar{x}_{d,\\cdot}) = \\mathrm{sd}(\\bar{x}_{d,\\cdot})/\\sqrt{n} \\)",
  "- **MAD**: center \\(\\pm k\\cdot \\mathrm{MAD}(\\bar{x}_{d,\\cdot})\\)",
  "- **Bootstrap CI** (optional): resample SampleIDs within Day and recompute the center.",
  "",
  "### Panel B: system dispersion across features (within-day)",
  "",
  "For each Day \\(d\\) and Feature \\(f\\), compute a within-day dispersion metric over SampleIDs:",
  "",
  "- SD: \\(\\mathrm{sd}(x_{d,\\cdot,f})\\)",
  "- CV: \\(\\mathrm{sd}(x_{d,\\cdot,f}) / |\\mathrm{mean}(x_{d,\\cdot,f})|\\)",
  "- IQR: \\(\\mathrm{IQR}(x_{d,\\cdot,f})\\)",
  "- MAD: \\(\\mathrm{mad}(x_{d,\\cdot,f};\\, c=1)\\)",
  "",
  "Aggregate feature-wise dispersions within each Day using the median (default) or mean:",
  "",
  "\\[ D_d = \\mathrm{median}_f\\{\\delta_{d,f}\\} \\quad \\text{or} \\quad D_d = \\mathrm{mean}_f\\{\\delta_{d,f}\\} \\]",
  "",
  "Intervals are shown either as a descriptive MAD spread across features (\\(D_d \\pm k\\cdot \\mathrm{MAD}_f(\\delta_{d,f})\\)) or as a bootstrap CI.",
  "",
  "### Panel D: CV vs Dispersion quadrant",
  "",
  "This panel separates **absolute variability** (dispersion) from **scaling effects** (CV). Days with low dispersion but high CV are immediately interpretable as low absolute variability at low baseline abundance.",
  "",
  "### Panel E: Mean–dispersion decoupling index",
  "",
  "Let \\(\\mu(d)\\) be the Day-level mean center from Panel A and \\(D(d)\\) be the dispersion from Panel B.",
  "Compute robust z-scores across days within source:",
  "",
  "\\[ Z_{\\mu}(d) = \\frac{\\mu(d) - \\mathrm{median}(\\mu)}{\\mathrm{MAD}(\\mu)} \\quad,\\quad",
  "   Z_{D}(d) = \\frac{D(d) - \\mathrm{median}(D)}{\\mathrm{MAD}(D)} \\]",
  "",
  "Then the decoupling index is:",
  "",
  "\\[ \\Delta(d) = Z_{\\mu}(d) - Z_{D}(d) \\]",
  "",
  "Positive values indicate mean is high relative to dispersion (a canalization-consistent regime).",
  "",
  "## Plots",
  "",
  "```{r}",
  "plot_path <- man$generated_files$preferred_plot",
  "if (is.null(plot_path) || is.na(plot_path) || !file.exists(plot_path)) {",
  "  pngs <- list.files('.', pattern = '\\\\.png$', full.names = TRUE)",
  "  if (length(pngs) > 0) plot_path <- pngs[1] else plot_path <- NA_character_",
  "}",
  "if (!is.na(plot_path) && file.exists(plot_path)) knitr::include_graphics(plot_path) else cat('No PNG plot found.')",
  "```",
  "",
  "## Interpretation of results",
  "",
  "- **Panel A**: system mean trajectory with replicate-level uncertainty across SampleIDs.",
  "- **Panel B**: system-wide within-day heterogeneity aggregated across features; interval indicates uncertainty (CI) or descriptive spread (MAD).",
  "- **Panel C**: makes timing of maximal dispersion unambiguous.",
  "- **Panel D**: distinguishes absolute dispersion from CV scaling.",
  "- **Panel E**: single-curve summary of mean–dispersion decoupling relevant to canalization framing.",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, qmd_path)

render_ok <- FALSE
render_err <- NULL
tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_file = html_name,
    output_dir = output_dir,
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
  render_ok <- file.exists(html_path)
}, error = function(e) {
  render_err <<- e$message
  render_ok <- FALSE
})

# After rendering, refresh inventory (mandatory)
inv_final <- refresh_file_inventory()

manifest$generated_files <- list(
  inventory_csv = normalizePath(manifest_path_files, winslash = "/", mustWork = FALSE),
  n_files = nrow(inv_final),
  preferred_plot = manifest$generated_files$preferred_plot,
  qmd = normalizePath(qmd_path, winslash = "/", mustWork = FALSE),
  html = if (render_ok) normalizePath(html_path, winslash = "/", mustWork = FALSE) else NA_character_
)
write_manifest_json(manifest)

if (render_ok) {
  message("HTML report created: ", normalizePath(html_path, winslash = "/", mustWork = FALSE))
} else {
  message("FAILED to render HTML report.")
  message("QMD written at: ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE))
  if (!is.null(render_err)) message("Render error: ", render_err)
}

message("\nDone.")
