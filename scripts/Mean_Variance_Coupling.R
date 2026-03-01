#!/usr/bin/env Rscript
# =============================================================================
# Mean–Variance Decoupling Control (Gene Level)
#
# PURPOSE
#   Provide a reviewer-proof control demonstrating that observed variance
#   restructuring across developmental time is not a trivial artifact of
#   mean-expression shifts (mean–variance dependence / heteroscedasticity).
#
# INPUT (broadly compatible "tidy long" table)
#   A delimited text file (CSV/TSV) in long format with (at minimum):
#     - a time/day column (e.g., Day)
#     - a feature column (e.g., Gene)
#     - a replicate/sample column (e.g., Bird, SampleID)
#     - a value column (e.g., Expression; log2-scale recommended)
#
# OUTPUTS (written ONLY under outputs/<run_id>/)
#   - MeanVariance_ByDay_Summary.csv
#   - MeanVariance_GeneLevel_ByDay.csv (optional; can be large)
#   - Diagnostic plots (PNG):
#       * Plot_ResidualLogVar_ByDay.png
#       * Plot_MeanVariance_Scatter_Day<mid>.png
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - <script_name>_Report.qmd + rendered HTML
#
# METHODS (core idea)
#   For each day d and feature g:
#     mu_{g,d}  = mean(value across replicates)
#     var_{g,d} = variance(value across replicates)
#   Fit mean–variance trend within each day:
#     log(var_{g,d}) ~ f(mu_{g,d}) using LOESS (or fallback to linear)
#   Residual variance:
#     resid_{g,d} = log(var_{g,d}) - fhat(mu_{g,d})
#   Summarize per day:
#     median(resid_{g,d}), IQR, N features
#
# NOTES
#   - This script does NOT hardcode the input file.
#   - The script ALWAYS asks whether to apply a transform:
#       1 = no transform
#       2 = log2(Value + pseudocount)
# =============================================================================

# ----------------------------- Session Start ---------------------------------
start_time <- Sys.time()

# ----------------------------- User Toggles ----------------------------------
analysis_name         <- "MeanVarianceDecoupling"
enable_plotly         <- TRUE     # interactive versions in report when possible
runtime_tracking      <- TRUE
write_gene_level_csv  <- TRUE     # can be large for big feature sets

# LOESS controls (robust defaults)
loess_span            <- 0.75
min_genes_per_day     <- 50       # below this, fallback to linear (more stable)
min_reps_per_gene_day <- 3        # variance needs >=2; 3 is safer for stability

# -------------------------- Dependency Handling ------------------------------
quiet_install <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

base_pkgs   <- c("data.table", "ggplot2", "jsonlite", "rmarkdown", "knitr")
plotly_pkgs <- c("plotly", "htmlwidgets")
rstudio_pkg <- c("rstudioapi") # optional, only for script path + nicer UX

quiet_install(base_pkgs)
quiet_install(rstudio_pkg)
if (isTRUE(enable_plotly)) quiet_install(plotly_pkgs)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(rmarkdown)
  library(knitr)
})

# ---------------------- Script identity capture (MANDATORY) -------------------
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
known_script_filename <- "MeanVarianceDecoupling_Control.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  # Path cannot be detected in this execution mode; still record a valid script_name.
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  message("NOTE: Script path detection failed; using fallback filename: ", known_script_filename)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

# ----------------------------- Console Summary --------------------------------
cat("\n==================== Script Identity ====================\n")
cat("script_name: ", script_name, "\n", sep = "")
cat("script_path: ", script_path, "\n", sep = "")
cat("script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
cat("=========================================================\n\n")

# ------------------------------ Output Policy --------------------------------
if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0(analysis_name, "_Source_", timestamp)  # updated after input selection
output_dir <- file.path("outputs", run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: write file inventory CSV
write_inventory <- function(out_dir, inventory_path) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    inv <- data.frame(file = character(), bytes = numeric(), modified = character(), stringsAsFactors = FALSE)
  } else {
    info <- file.info(files)
    inv <- data.frame(
      file = normalizePath(files, winslash = "/", mustWork = FALSE),
      bytes = info$size,
      modified = as.character(info$mtime),
      stringsAsFactors = FALSE
    )
  }
  data.table::fwrite(inv, inventory_path)
}

# ------------------------------ Input Loading --------------------------------
pick_input_file <- function() {
  # Prefer interactive chooser when possible
  p <- tryCatch({
    if (interactive()) file.choose() else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) return(p)
  
  # Non-interactive: allow commandArgs --input=...
  ca <- commandArgs(trailingOnly = FALSE)
  in_arg <- grep("^--input=", ca, value = TRUE)
  if (length(in_arg) > 0) {
    p2 <- sub("^--input=", "", in_arg[1])
    if (nzchar(p2) && file.exists(p2)) return(p2)
  }
  
  stop("No input file selected/found. Run interactively to use file chooser, or pass --input=/path/to/file.csv")
}

input_path <- pick_input_file()
input_path_norm <- normalizePath(input_path, winslash = "/", mustWork = FALSE)
source_label <- tools::file_path_sans_ext(basename(input_path_norm))

# Update run_id with source label (deterministic-ish, informative)
run_id <- paste0(analysis_name, "_", source_label, "_", timestamp)
output_dir <- file.path("outputs", run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read CSV/TSV robustly
read_delim_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    return(data.table::fread(path, sep = "\t", data.table = TRUE))
  }
  data.table::fread(path, data.table = TRUE)
}

dt <- read_delim_any(input_path_norm)
if (nrow(dt) == 0) stop("Input appears empty: ", input_path_norm)

# ------------------------- Column Mapping (broad) -----------------------------
guess_col <- function(cands, names_vec) {
  hit <- intersect(tolower(cands), tolower(names_vec))
  if (length(hit) == 0) return(NA_character_)
  names_vec[match(hit[1], tolower(names_vec))]
}

cn <- names(dt)

day_col  <- guess_col(c("day", "time", "age", "d", "timepoint"), cn)
feat_col <- guess_col(c("gene", "feature", "analyte", "metabolite", "id"), cn)
rep_col  <- guess_col(c("bird", "sample", "sampleid", "replicate", "individual"), cn)
val_col  <- guess_col(c("value", "expr", "expression", "abundance", "intensity", "log2", "log2expr"), cn)

choose_col <- function(prompt, cn) {
  if (!interactive()) stop("Cannot resolve columns non-interactively. Provide standardized column names or run interactively.")
  idx <- menu(cn, title = prompt)
  if (idx < 1) stop("No column selected for: ", prompt)
  cn[idx]
}

if (any(is.na(c(day_col, feat_col, rep_col, val_col)))) {
  message("Column auto-detection incomplete. Please map columns.")
  if (is.na(day_col))  day_col  <- choose_col("Select the Day/Time column", cn)
  if (is.na(feat_col)) feat_col <- choose_col("Select the Feature (Gene/Metabolite) column", cn)
  if (is.na(rep_col))  rep_col  <- choose_col("Select the Replicate/Sample column", cn)
  if (is.na(val_col))  val_col  <- choose_col("Select the Value/Expression column", cn)
}

# Standardize working columns
dt2 <- dt[, .(
  Day       = suppressWarnings(as.integer(get(day_col))),
  Feature   = as.character(get(feat_col)),
  Replicate = as.character(get(rep_col)),
  Value     = suppressWarnings(as.numeric(get(val_col)))
)]

dt2 <- dt2[is.finite(Day) & nzchar(Feature) & nzchar(Replicate) & is.finite(Value)]
if (nrow(dt2) == 0) stop("After cleaning, no valid rows remain. Check Day/Feature/Replicate/Value columns.")

# ---------------------- ALWAYS ASK: Transform Mode ----------------------------
transform_mode <- 1L
pseudocount <- 1.0

ask_transform <- function() {
  cat("\n================ Value Transform =================\n")
  cat("Choose transformation for Value prior to analysis:\n")
  cat("  1 = No transform\n")
  cat("  2 = log2(Value + pseudocount)\n")
  cat("==================================================\n")
  ans <- readline("Enter 1 or 2 (default 1): ")
  if (!nzchar(ans)) return(list(mode = 1L, pseudocount = 1.0))
  mode <- suppressWarnings(as.integer(ans))
  if (is.na(mode) || !mode %in% c(1L, 2L)) mode <- 1L
  
  pc <- 1.0
  if (mode == 2L) {
    ans2 <- readline("Enter pseudocount (default 1.0; must be >= 0): ")
    if (nzchar(ans2)) pc <- suppressWarnings(as.numeric(ans2))
    if (!is.finite(pc) || pc < 0) pc <- 1.0
  }
  list(mode = mode, pseudocount = pc)
}

if (interactive()) {
  tr <- ask_transform()
  transform_mode <- tr$mode
  pseudocount <- tr$pseudocount
} else {
  # Non-interactive fallback: allow --transform=1/2 and --pseudocount=...
  ca <- commandArgs(trailingOnly = FALSE)
  t_arg <- grep("^--transform=", ca, value = TRUE)
  if (length(t_arg) > 0) {
    transform_mode <- suppressWarnings(as.integer(sub("^--transform=", "", t_arg[1])))
    if (is.na(transform_mode) || !transform_mode %in% c(1L, 2L)) transform_mode <- 1L
  }
  pc_arg <- grep("^--pseudocount=", ca, value = TRUE)
  if (length(pc_arg) > 0) {
    pseudocount <- suppressWarnings(as.numeric(sub("^--pseudocount=", "", pc_arg[1])))
    if (!is.finite(pseudocount) || pseudocount < 0) pseudocount <- 1.0
  }
  cat("NOTE: Non-interactive run; using transform_mode=", transform_mode,
      ", pseudocount=", pseudocount, "\n", sep = "")
}

# Apply transform
if (transform_mode == 2L) {
  if (!is.finite(pseudocount) || pseudocount < 0) {
    stop("pseudocount must be a finite, non-negative number.")
  }
  if (any(dt2$Value + pseudocount <= 0, na.rm = TRUE)) {
    stop("log2(Value + pseudocount) invalid: some Value + pseudocount <= 0. ",
         "Increase pseudocount or verify input scale.")
  }
  dt2[, Value := log2(Value + pseudocount)]
}

# ----------------------- Compute mean/variance per gene/day -------------------
rep_counts <- dt2[, .N, by = .(Day, Feature)]
keep_pairs <- rep_counts[N >= min_reps_per_gene_day, .(Day, Feature)]
setkey(keep_pairs, Day, Feature)
setkey(dt2, Day, Feature)
dt2k <- dt2[keep_pairs]

gd <- dt2k[, .(
  mu    = mean(Value, na.rm = TRUE),
  var   = var(Value, na.rm = TRUE),
  n_rep = .N
), by = .(Day, Feature)]

# log variance; guard against zero
eps <- 1e-12
gd[, log_var := log(pmax(var, eps))]

# ---------------------- Residualize mean–variance trend -----------------------
fit_residuals_one_day <- function(df_day, span = loess_span, min_genes = min_genes_per_day) {
  df_day <- df_day[is.finite(mu) & is.finite(log_var)]
  if (nrow(df_day) < 10) {
    df_day[, trend := NA_real_]
    df_day[, resid_log_var := NA_real_]
    df_day[, fit_method := "insufficient_genes"]
    return(df_day)
  }
  
  if (nrow(df_day) < min_genes) {
    fit <- lm(log_var ~ mu, data = df_day)
    df_day[, trend := predict(fit, newdata = df_day)]
    df_day[, resid_log_var := log_var - trend]
    df_day[, fit_method := "linear_fallback"]
    return(df_day)
  }
  
  ok <- TRUE
  fit <- tryCatch({
    loess(log_var ~ mu, data = df_day, span = span, degree = 2,
          control = loess.control(surface = "direct"))
  }, error = function(e) { ok <<- FALSE; NULL })
  
  if (!ok || is.null(fit)) {
    fit2 <- lm(log_var ~ mu, data = df_day)
    df_day[, trend := predict(fit2, newdata = df_day)]
    df_day[, resid_log_var := log_var - trend]
    df_day[, fit_method := "linear_after_loess_fail"]
    return(df_day)
  }
  
  df_day[, trend := predict(fit, newdata = df_day)]
  df_day[, resid_log_var := log_var - trend]
  df_day[, fit_method := "loess"]
  df_day
}

gd_res <- rbindlist(lapply(split(gd, gd$Day), fit_residuals_one_day), use.names = TRUE, fill = TRUE)

# ------------------------------ Summaries -------------------------------------
by_day <- gd_res[, .(
  n_features              = sum(is.finite(resid_log_var)),
  median_mu               = median(mu, na.rm = TRUE),
  median_log_var          = median(log_var, na.rm = TRUE),
  median_resid_log_var    = median(resid_log_var, na.rm = TRUE),
  iqr_resid_log_var       = IQR(resid_log_var, na.rm = TRUE),
  frac_linear_fallback    = mean(fit_method %in% c("linear_fallback", "linear_after_loess_fail"), na.rm = TRUE)
), by = .(Day)][order(Day)]

# ------------------------------ Write outputs --------------------------------
summary_csv <- file.path(output_dir, "MeanVariance_ByDay_Summary.csv")
data.table::fwrite(by_day, summary_csv)

gene_level_csv <- file.path(output_dir, "MeanVariance_GeneLevel_ByDay.csv")
if (isTRUE(write_gene_level_csv)) {
  data.table::fwrite(gd_res, gene_level_csv)
}

# ------------------------------ Plotting -------------------------------------
plot1_path <- file.path(output_dir, "Plot_ResidualLogVar_ByDay.png")
p_resid <- ggplot(by_day, aes(x = Day, y = median_resid_log_var)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  labs(
    title = "Mean–Variance Decoupling Control",
    subtitle = "Median residual log-variance per day (variance after controlling for mean)",
    x = "Day",
    y = "Median residual log(variance)"
  ) +
  theme_bw(base_size = 11)

ggsave(plot1_path, p_resid, width = 7.5, height = 4.5, dpi = 300)

# Mean–variance scatter for a representative day (median day)
days <- sort(unique(gd_res$Day))
pick_day <- days[ceiling(length(days) / 2)]
df_sc <- gd_res[Day == pick_day & is.finite(mu) & is.finite(log_var)]
plot2_path <- file.path(output_dir, paste0("Plot_MeanVariance_Scatter_Day", pick_day, ".png"))

p_sc <- ggplot(df_sc, aes(x = mu, y = log_var)) +
  geom_point(alpha = 0.6, size = 1.1) +
  labs(
    title = paste0("Mean–Variance Relationship (Day ", pick_day, ")"),
    subtitle = "Each point is a feature; log(variance) across replicates vs mean value",
    x = "Mean value (mu)",
    y = "log(variance)"
  ) +
  theme_bw(base_size = 11)

ggsave(plot2_path, p_sc, width = 7.0, height = 5.0, dpi = 300)

# ------------------------------ Manifest -------------------------------------
get_pkg_versions <- function(pkgs) {
  lapply(pkgs, function(p) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    list(package = p, version = v)
  })
}

deps <- unique(c(base_pkgs, if (isTRUE(enable_plotly)) plotly_pkgs else character(), rstudio_pkg))
manifest_path  <- file.path(output_dir, "Project_Manifest.json")
inventory_path <- file.path(output_dir, "Project_Manifest_Files.csv")

# runtime placeholder; filled later if enabled
runtime_seconds <- NA_real_

manifest <- list(
  run_id = run_id,
  run_timestamp = format(start_time, "%Y-%m-%d %H:%M:%S"),
  script = list(
    name = script_name,
    path = script_path,
    full_path = if (is.na(script_full)) NULL else script_full
  ),
  input = list(
    input_path = input_path_norm,
    source_label = source_label,
    detected_columns = list(day_col = day_col, feature_col = feat_col, replicate_col = rep_col, value_col = val_col)
  ),
  parameters = list(
    analysis_name = analysis_name,
    enable_plotly = enable_plotly,
    runtime_tracking = runtime_tracking,
    write_gene_level_csv = write_gene_level_csv,
    loess_span = loess_span,
    min_genes_per_day = min_genes_per_day,
    min_reps_per_gene_day = min_reps_per_gene_day,
    transform_mode = transform_mode,
    pseudocount = pseudocount,
    runtime_seconds = runtime_seconds
  ),
  dependencies = get_pkg_versions(deps),
  outputs = list(
    outputs_root = normalizePath("outputs", winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  ),
  generated_files = list()
)

write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")

# ------------------------------ QMD Report -----------------------------------
read_header_block <- function(script_full, max_lines = 140) {
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
  } else {
    lines <- c("# Script header unavailable (script_full not detected).")
  }
  if (length(lines) == 0) return("Header not available.")
  out <- character()
  for (i in seq_len(min(length(lines), max_lines))) {
    if (grepl("^\\s*#", lines[i]) || i == 1) {
      out <- c(out, lines[i])
    } else {
      break
    }
  }
  paste(out, collapse = "\n")
}

header_txt <- read_header_block(script_full)

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_out <- file.path(output_dir, paste0(script_name, "_Report.html"))

qmd <- c(
  "---",
  paste0('title: "', script_name, ' Report"'),
  "format:",
  "  html:",
  "    toc: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "params:",
  paste0('  run_id: "', run_id, '"'),
  paste0('  input_path: "', gsub('"', '\\"', input_path_norm), '"'),
  paste0('  output_dir: "', gsub('"', '\\"', normalizePath(output_dir, winslash="/", mustWork=FALSE)), '"'),
  paste0('  enable_plotly: ', ifelse(enable_plotly, "true", "false")),
  paste0('  transform_mode: ', transform_mode),
  paste0('  pseudocount: ', pseudocount),
  paste0('  runtime_seconds: ', ifelse(is.na(runtime_seconds), "null", runtime_seconds)),
  "---",
  "",
  "# Overview",
  "This report documents a mean–variance decoupling control to test whether variance restructuring",
  "across days can be explained by mean-expression dependence.",
  "",
  "## Metadata",
  "```{r}",
  "cat('Run ID: ', params$run_id, '\\n', sep='')",
  "cat('Input: ', params$input_path, '\\n', sep='')",
  "cat('Output dir: ', params$output_dir, '\\n', sep='')",
  "cat('Interactive plots enabled: ', params$enable_plotly, '\\n', sep='')",
  "cat('Transform mode: ', params$transform_mode, ' (1=no transform, 2=log2+pc)\\n', sep='')",
  "cat('Pseudocount: ', params$pseudocount, '\\n', sep='')",
  "```",
  "",
  "## Runtime",
  "```{r}",
  "if (!is.null(params$runtime_seconds)) {",
  "  cat('Total runtime (seconds): ', params$runtime_seconds, '\\n', sep='')",
  "} else {",
  "  cat('Runtime tracking was not recorded.\\n')",
  "}",
  "```",
  "",
  "# Script Header (embedded)",
  "```",
  header_txt,
  "```",
  "",
  "# Dependencies",
  "```{r}",
  "library(jsonlite)",
  "m <- fromJSON(file.path(params$output_dir, 'Project_Manifest.json'))",
  "deps <- as.data.frame(m$dependencies)",
  "knitr::kable(deps)",
  "```",
  "",
  "# Generated Files",
  "```{r}",
  "library(data.table)",
  "inv <- fread(file.path(params$output_dir, 'Project_Manifest_Files.csv'))",
  "knitr::kable(inv)",
  "```",
  "",
  "# Analytical logic and formulas",
  "We compute, for each day and feature (gene/metabolite):",
  "",
  "- Mean across replicates:  \\(\\mu_{g,d}\\)",
  "- Variance across replicates: \\(s^2_{g,d}\\)",
  "- Mean–variance trend within day: \\(\\log s^2_{g,d} \\sim f(\\mu_{g,d})\\) via LOESS (fallback linear)",
  "- Residual log-variance: \\(r_{g,d} = \\log s^2_{g,d} - \\hat f(\\mu_{g,d})\\)",
  "",
  "We then summarize \\(r_{g,d}\\) across features within day (median, IQR).",
  "",
  "# Results",
  "## Required diagnostic plot (static)",
  "```{r}",
  "knitr::include_graphics(file.path(params$output_dir, 'Plot_ResidualLogVar_ByDay.png'))",
  "```",
  "",
  "## Mean–variance scatter (static)",
  "```{r}",
  "files <- list.files(params$output_dir, pattern='Plot_MeanVariance_Scatter_', full.names=TRUE)",
  "if (length(files) > 0) knitr::include_graphics(files[1])",
  "```",
  "",
  "## Interactive plot (if enabled)",
  "```{r}",
  "ok <- FALSE",
  "if (isTRUE(params$enable_plotly)) {",
  "  ok <- tryCatch({",
  "    library(plotly)",
  "    library(data.table)",
  "    s <- fread(file.path(params$output_dir, 'MeanVariance_ByDay_Summary.csv'))",
  "    g <- ggplot(s, aes(x=Day, y=median_resid_log_var)) + geom_line() + geom_point() +",
  "      labs(title='Interactive: Median residual log-variance by day', x='Day', y='Median residual log(variance)')",
  "    print(plotly::ggplotly(g))",
  "    TRUE",
  "  }, error = function(e) {",
  "    FALSE",
  "  })",
  "}",
  "if (!ok) {",
  "  cat('Interactive plot unavailable (plotly disabled or failed). Static plots above provide diagnostics.\\n')",
  "}",
  "```",
  "",
  "# Interpretation",
  "If the day/phase pattern in dispersion persists in the residual log-variance summary,",
  "then variance restructuring cannot be explained by mean-expression shifts alone.",
  "",
  "# Reproducibility",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd, qmd_path)

# ------------------------- End time + runtime (MANDATORY) ---------------------
end_time <- Sys.time()
if (isTRUE(runtime_tracking)) {
  runtime_seconds <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  manifest$parameters$runtime_seconds <- runtime_seconds
  write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  qmd_lines <- readLines(qmd_path, warn = FALSE)
  qmd_lines <- sub("^  runtime_seconds: null$", paste0("  runtime_seconds: ", runtime_seconds), qmd_lines)
  writeLines(qmd_lines, qmd_path)
}

# ------------------------------ Render Report --------------------------------
render_ok <- TRUE
render_err <- NULL

tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_file = basename(html_out),
    output_dir = output_dir,
    quiet = TRUE,
    params = list(
      run_id = run_id,
      input_path = input_path_norm,
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      enable_plotly = enable_plotly,
      transform_mode = transform_mode,
      pseudocount = pseudocount,
      runtime_seconds = ifelse(isTRUE(runtime_tracking), runtime_seconds, NULL)
    )
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- conditionMessage(e)
})

# ------------------------------ Inventory Refresh -----------------------------
write_inventory(output_dir, inventory_path)

inv_dt <- tryCatch(data.table::fread(inventory_path), error = function(e) NULL)
if (!is.null(inv_dt)) {
  manifest$generated_files <- inv_dt$file
  write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
}

# ------------------------------ Console Footer --------------------------------
if (render_ok && file.exists(html_out)) {
  cat("\nHTML report created: ", normalizePath(html_out, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("\nREPORT RENDER FAILED.\n")
  if (!is.null(render_err)) cat("Error: ", render_err, "\n", sep = "")
  cat("QMD saved at: ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
}

cat("\nOutput directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Manifest: ", normalizePath(manifest_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("File inventory: ", normalizePath(inventory_path, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
