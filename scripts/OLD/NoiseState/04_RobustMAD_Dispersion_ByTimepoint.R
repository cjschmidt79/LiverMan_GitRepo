#!/usr/bin/env Rscript
# ======================================================================
# Robust MAD Dispersion (MAD_resid) Analysis (Input-Agnostic; Manifest + Report)
#
# PURPOSE
#   Compute statistically orthodox, robust within-timepoint dispersion using MAD:
#     - Feature-level MAD across replicates within each timepoint
#     - Timepoint summary: D_t = median_i( MAD_{i,t} )
#   Optionally compute MAD_resid by removing a global mean–variance trend:
#     - Fit loess on log(MAD) ~ log(mean) across all (feature, timepoint) points
#     - Residualize: logMAD_resid = logMAD - loess_fit(logMean)
#     - Summarize per timepoint: D_t_resid = median_i( logMAD_resid_{i,t} )
#
# INPUTS (either format is accepted)
#   A) "Tidy long" table:
#        - Columns resembling: Day, SampleID, Feature (Gene/Metabolite), Abundance
#        - Example: Transcriptome_Tidy_Long_FIXED4columnOnly.csv
#
#   B) "Wide" matrix-like table:
#        - One row per sample (replicate), one column Day, optional SampleID,
#          and remaining columns are features.
#        - Example: WIDE_LiverMetabolome.csv
#
# OUTPUTS
#   All outputs are written to: outputs/<run_id>/
#     - MAD_ByFeatureDay.csv
#     - MAD_ByDay_Summary.csv
#     - Plots (PNG)
#     - Project_Manifest.json
#     - Project_Manifest_Files.csv
#     - <script_name>_Report.qmd and rendered HTML report
#
# NOTES
#   - Designed to be robust to outliers and mean dependence (via residualization).
#   - Keeps computations orthodox and transparent; avoids overclaiming.
# ======================================================================

# ----------------------------- Dependency bootstrap -----------------------------
quiet_install <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

required_pkgs <- c(
  "jsonlite", "rmarkdown", "knitr",
  "ggplot2", "dplyr", "tidyr", "readr", "tibble", "stringr"
)
quiet_install(required_pkgs)

suppressPackageStartupMessages({
  library(jsonlite)
  library(rmarkdown)
  library(knitr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
})

# ----------------------------- Script identity capture (MANDATORY) -----------------------------
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
known_script_filename <- "RobustMAD_Analysis.R"
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

# Console summary (MANDATORY)
cat("\n==================== Script Identity ====================\n")
cat("script_name:", script_name, "\n")
cat("script_path:", script_path, "\n")
cat("script_full:", ifelse(is.na(script_full), "NA", script_full), "\n")
if (is.na(script_full)) {
  cat("NOTE: Script full path detection failed; using fallback filename:", known_script_filename, "\n")
}
cat("=========================================================\n\n")

# ----------------------------- Output directory policy (MANDATORY) -----------------------------
outputs_root <- normalizePath(file.path(getwd(), "outputs"), winslash = "/", mustWork = FALSE)
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Utilities -----------------------------
stop2 <- function(...) stop(paste0(...), call. = FALSE)

now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

get_deps <- function(pkgs) {
  deps <- lapply(pkgs, function(p) {
    ver <- NA_character_
    if (requireNamespace(p, quietly = TRUE)) ver <- as.character(utils::packageVersion(p))
    data.frame(package = p, version = ver, stringsAsFactors = FALSE)
  })
  do.call(rbind, deps)
}

# Robust MAD (unscaled) and a scaled version (sigma-approx) if desired
mad_unscaled <- function(x) stats::mad(x, center = stats::median(x, na.rm = TRUE),
                                       constant = 1, na.rm = TRUE)
mad_scaled <- function(x) stats::mad(x, center = stats::median(x, na.rm = TRUE),
                                     constant = 1.4826, na.rm = TRUE)

coerce_day_numeric <- function(x) {
  # Accept numeric, integer, or strings like "Day14" / "D14"
  if (is.numeric(x)) return(as.numeric(x))
  s <- as.character(x)
  s2 <- gsub("[^0-9\\.\\-]+", "", s)
  suppressWarnings(as.numeric(s2))
}

# Identify columns case-insensitively
find_col <- function(nms, candidates) {
  nms_low <- tolower(nms)
  cand_low <- tolower(candidates)
  hit <- which(nms_low %in% cand_low)
  if (length(hit) > 0) return(nms[hit[1]])
  # try contains matches
  hit2 <- which(vapply(cand_low, function(cc) any(grepl(cc, nms_low, fixed = TRUE)), logical(1)))
  if (length(hit2) > 0) {
    cc <- cand_low[hit2[1]]
    idx <- which(grepl(cc, nms_low, fixed = TRUE))[1]
    return(nms[idx])
  }
  NA_character_
}

# ----------------------------- Input selection (interactive + CLI) -----------------------------
args <- commandArgs(trailingOnly = TRUE)

# Supported CLI:
#   Rscript RobustMAD_Analysis.R --input=/path/to/file.csv --analysis_name=MADDispersion --residualize=TRUE
parse_arg <- function(flag) {
  # expects --flag=value
  hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^--", flag, "="), "", hit[1])
}

input_path_cli <- parse_arg("input")
analysis_name_cli <- parse_arg("analysis_name")
residualize_cli <- parse_arg("residualize")
mad_scale_cli <- parse_arg("mad_scale") # "unscaled" or "scaled"

analysis_name <- if (!is.na(analysis_name_cli) && nzchar(analysis_name_cli)) analysis_name_cli else "MADDispersion"
do_residualize <- if (!is.na(residualize_cli) && nzchar(residualize_cli)) {
  tolower(residualize_cli) %in% c("true", "t", "1", "yes", "y")
} else TRUE

mad_scale_mode <- if (!is.na(mad_scale_cli) && nzchar(mad_scale_cli)) tolower(mad_scale_cli) else "unscaled"
if (!mad_scale_mode %in% c("unscaled", "scaled")) mad_scale_mode <- "unscaled"

mad_fun <- if (mad_scale_mode == "scaled") mad_scaled else mad_unscaled

choose_input <- function() {
  if (!is.na(input_path_cli) && nzchar(input_path_cli)) {
    p <- input_path_cli
    if (!file.exists(p)) stop2("Input file does not exist: ", p)
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  if (interactive()) {
    cat("[Input] Select a CSV file (either tidy-long or wide format).\n")
    p <- file.choose()
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  stop2("Non-interactive mode requires --input=/path/to/file.csv")
}

input_path <- choose_input()

# ----------------------------- Read input (agnostic) -----------------------------
df_raw <- suppressWarnings(readr::read_csv(input_path, show_col_types = FALSE, progress = FALSE))

if (nrow(df_raw) == 0) stop2("Input file appears empty: ", input_path)

# Determine format
nms <- names(df_raw)

day_col <- find_col(nms, c("day", "time", "timepoint", "tp"))
sample_col <- find_col(nms, c("sampleid", "sample", "replicate", "bird", "id", "subject"))
feature_col <- find_col(nms, c("feature", "gene", "metabolite", "analyte", "compound"))
abundance_col <- find_col(nms, c("abundance", "value", "count", "counts", "intensity", "area", "signal"))

has_tidy_signature <- !is.na(day_col) && !is.na(feature_col) && !is.na(abundance_col)

# Heuristic: if feature_col & abundance_col exist, treat as tidy-long; else treat as wide.
input_format <- if (has_tidy_signature) "TidyLong" else "Wide"

# Derive Source label for run_id
source_label <- input_format

# run_id (deterministic + informative)
run_id <- paste0(analysis_name, "_", source_label, "_", now_stamp())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

cat("\n==================== Run Configuration ===================\n")
cat("run_id:", run_id, "\n")
cat("input_path:", input_path, "\n")
cat("input_format:", input_format, "\n")
cat("analysis_name:", analysis_name, "\n")
cat("do_residualize:", do_residualize, "\n")
cat("mad_scale_mode:", mad_scale_mode, "\n")
cat("output_dir:", output_dir, "\n")
cat("=========================================================\n\n")

# ----------------------------- Canonicalize to tidy long: Day, SampleID, Feature, Value -----------------------------
if (input_format == "TidyLong") {
  # Ensure sample column; if missing, synthesize replicate index within day
  if (is.na(sample_col)) {
    sample_col <- ".SampleID"
    df_raw[[sample_col]] <- ave(seq_len(nrow(df_raw)), df_raw[[day_col]], FUN = seq_along)
  }
  
  df_long <- df_raw %>%
    transmute(
      Day = coerce_day_numeric(.data[[day_col]]),
      SampleID = as.character(.data[[sample_col]]),
      Feature = as.character(.data[[feature_col]]),
      Value = as.numeric(.data[[abundance_col]])
    ) %>%
    filter(!is.na(Day), !is.na(Feature), !is.na(SampleID)) %>%
    filter(!is.na(Value))
  
} else {
  # Wide: requires Day; SampleID optional; all other non-Day/sample columns are features
  if (is.na(day_col)) {
    stop2("Wide format detected but could not find a Day/Time column. Columns are: ",
          paste(names(df_raw), collapse = ", "))
  }
  
  if (is.na(sample_col)) {
    # Create a stable per-row SampleID
    sample_col <- ".SampleID"
    df_raw[[sample_col]] <- paste0("S", seq_len(nrow(df_raw)))
  }
  
  feature_cols <- setdiff(names(df_raw), c(day_col, sample_col))
  if (length(feature_cols) < 2) {
    stop2("Wide format requires at least 2 feature columns besides Day/SampleID.")
  }
  
  df_long <- df_raw %>%
    mutate(
      Day = coerce_day_numeric(.data[[day_col]]),
      SampleID = as.character(.data[[sample_col]])
    ) %>%
    select(Day, SampleID, all_of(feature_cols)) %>%
    pivot_longer(cols = all_of(feature_cols), names_to = "Feature", values_to = "Value") %>%
    mutate(
      Feature = as.character(Feature),
      Value = as.numeric(Value)
    ) %>%
    filter(!is.na(Day), !is.na(Value), !is.na(Feature), !is.na(SampleID))
}

if (nrow(df_long) == 0) stop2("No valid (Day, SampleID, Feature, Value) rows after parsing.")

# Basic guardrails
n_days <- n_distinct(df_long$Day)
n_features <- n_distinct(df_long$Feature)
if (n_days < 2) stop2("Need at least 2 unique timepoints (Day). Found: ", n_days)
if (n_features < 10) {
  cat("NOTE: Fewer than 10 features detected (", n_features, "). Analysis will run but interpret cautiously.\n", sep = "")
}

# ----------------------------- Compute feature-by-day statistics -----------------------------
# Feature-by-Day mean and MAD across replicates
by_fd <- df_long %>%
  group_by(Day, Feature) %>%
  summarise(
    n_reps = n_distinct(SampleID),
    mean_value = mean(Value, na.rm = TRUE),
    mad_value = mad_fun(Value),
    .groups = "drop"
  ) %>%
  filter(is.finite(mean_value), is.finite(mad_value))

# Some features may have n_reps == 1 at some days; MAD will be 0; keep but flag
by_fd <- by_fd %>%
  mutate(
    logMean = log1p(abs(mean_value)),
    logMAD  = log1p(abs(mad_value))
  )

# ----------------------------- Residualize mean–variance trend (optional) -----------------------------
if (do_residualize) {
  # Fit loess on all points; robust family for stability if available
  lo <- tryCatch(
    stats::loess(logMAD ~ logMean, data = by_fd, span = 0.75, degree = 2,
                 control = stats::loess.control(surface = "direct")),
    error = function(e) NULL
  )
  
  if (is.null(lo)) {
    cat("WARNING: loess fit failed; proceeding without residualization.\n")
    by_fd$logMAD_resid <- by_fd$logMAD
    residualize_used <- FALSE
  } else {
    pred <- tryCatch(stats::predict(lo, newdata = by_fd), error = function(e) rep(NA_real_, nrow(by_fd)))
    if (all(is.na(pred))) {
      cat("WARNING: loess prediction failed; proceeding without residualization.\n")
      by_fd$logMAD_resid <- by_fd$logMAD
      residualize_used <- FALSE
    } else {
      by_fd$logMAD_resid <- by_fd$logMAD - pred
      residualize_used <- TRUE
    }
  }
} else {
  by_fd$logMAD_resid <- by_fd$logMAD
  residualize_used <- FALSE
}

# ----------------------------- Timepoint summaries: D_t and D_t_resid -----------------------------
by_day <- by_fd %>%
  group_by(Day) %>%
  summarise(
    n_features = n_distinct(Feature),
    med_n_reps = median(n_reps, na.rm = TRUE),
    D_t_mad = median(mad_value, na.rm = TRUE),
    D_t_logMAD = median(logMAD, na.rm = TRUE),
    D_t_logMAD_resid = median(logMAD_resid, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Day)

# ----------------------------- Write primary outputs -----------------------------
out_by_fd <- file.path(output_dir, "MAD_ByFeatureDay.csv")
out_by_day <- file.path(output_dir, "MAD_ByDay_Summary.csv")

readr::write_csv(by_fd, out_by_fd)
readr::write_csv(by_day, out_by_day)

# ----------------------------- Plots -----------------------------
plot1_path <- file.path(output_dir, "MAD_ByDay_Dt.png")
p1 <- ggplot(by_day, aes(x = Day, y = D_t_logMAD_resid)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Robust Dispersion by Timepoint (D_t from logMAD_resid)",
    x = "Timepoint (Day)",
    y = "D_t = median_i(logMAD_resid)"
  ) +
  theme_classic(base_size = 12)
ggsave(plot1_path, p1, width = 8, height = 4.5, dpi = 300)

plot2_path <- file.path(output_dir, "MAD_ByFeatureDay_LogMeanVsLogMAD.png")
p2 <- ggplot(by_fd, aes(x = logMean, y = logMAD)) +
  geom_point(alpha = 0.35) +
  labs(
    title = "Mean–Variance Structure (log1p|mean| vs log1p|MAD|)",
    x = "log1p(|mean|)",
    y = "log1p(|MAD|)"
  ) +
  theme_classic(base_size = 12)
ggsave(plot2_path, p2, width = 7, height = 5, dpi = 300)

plot3_path <- file.path(output_dir, "MAD_ByFeatureDay_Resid_Distribution.png")
p3 <- ggplot(by_fd, aes(x = logMAD_resid)) +
  geom_histogram(bins = 60) +
  labs(
    title = "Distribution of logMAD_resid Across Feature×Day",
    x = "logMAD_resid",
    y = "Count"
  ) +
  theme_classic(base_size = 12)
ggsave(plot3_path, p3, width = 7, height = 5, dpi = 300)

# ----------------------------- Manifest construction (MANDATORY) -----------------------------
run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
deps_tbl <- get_deps(required_pkgs)

manifest <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA, script_full)
  ),
  input = list(
    input_path = input_path,
    input_format = input_format
  ),
  parameters = list(
    analysis_name = analysis_name,
    do_residualize = do_residualize,
    residualize_used = residualize_used,
    mad_scale_mode = mad_scale_mode,
    mad_definition = if (mad_scale_mode == "scaled")
      "stats::mad(constant=1.4826): sigma-approx scaled MAD"
    else
      "stats::mad(constant=1): unscaled MAD"
  ),
  dependencies = split(deps_tbl, seq_len(nrow(deps_tbl))),
  outputs = list(
    outputs_root = outputs_root,
    output_dir = output_dir
  )
)

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"),
           con = manifest_json_path)

# File inventory CSV (initial; will refresh after render)
inventory_path <- file.path(output_dir, "Project_Manifest_Files.csv")

write_inventory <- function(dir_path, out_csv) {
  files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    inv <- data.frame(
      file = character(0),
      rel_path = character(0),
      bytes = numeric(0),
      modified = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    info <- file.info(files)
    inv <- data.frame(
      file = normalizePath(files, winslash = "/", mustWork = FALSE),
      rel_path = gsub(paste0("^", normalizePath(dir_path, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(files, winslash = "/", mustWork = FALSE)),
      bytes = as.numeric(info$size),
      modified = as.character(info$mtime),
      stringsAsFactors = FALSE
    ) %>% arrange(rel_path)
  }
  readr::write_csv(inv, out_csv)
  inv
}

write_inventory(output_dir, inventory_path)

# ----------------------------- QMD report (MANDATORY) -----------------------------
# Extract header text (initial comment block) to embed in report
extract_header_block <- function(script_full_path) {
  # If script path is unknown, embed this script's header from the literal comments above.
  # Prefer reading the executing file if available.
  if (!is.na(script_full_path) && file.exists(script_full_path)) {
    lines <- readLines(script_full_path, warn = FALSE)
    # Grab from top until first non-comment, non-blank line after initial header
    # We'll stop when we hit the first line that doesn't start with '#' after having started.
    out <- character(0)
    started <- FALSE
    for (ln in lines) {
      if (!started) {
        if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
          started <- TRUE
          out <- c(out, ln)
        } else {
          # if first line is code, no header
          break
        }
      } else {
        if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
          out <- c(out, ln)
        } else {
          break
        }
      }
    }
    return(paste(out, collapse = "\n"))
  } else {
    # Fallback minimal header
    return(paste(
      "# Robust MAD Dispersion (MAD_resid) Analysis",
      "# (Header unavailable: script_full path not detected in this mode.)",
      sep = "\n"
    ))
  }
}

header_text <- extract_header_block(script_full)

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_out <- file.path(output_dir, paste0(script_name, "_Report.html"))

# Build QMD content
qmd <- c(
  "---",
  paste0("title: \"", script_name, " — Robust MAD Dispersion Report\""),
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
  "```{r}",
  "library(jsonlite)",
  "library(readr)",
  "library(dplyr)",
  "man <- jsonlite::fromJSON('Project_Manifest.json', simplifyVector = FALSE)",
  "inv <- readr::read_csv('Project_Manifest_Files.csv', show_col_types = FALSE)",
  "```",
  "",
  "```{r}",
  "meta_tbl <- tibble::tibble(",
  "  Field = c(",
  "    'run_id','run_timestamp',",
  "    'script_name','script_path','script_full_path',",
  "    'input_path','input_format',",
  "    'outputs_root','output_dir',",
  "    'do_residualize','residualize_used','mad_scale_mode'",
  "  ),",
  "  Value = c(",
  "    man$run_id, man$run_timestamp,",
  "    man$script$name, man$script$path, as.character(man$script$full_path),",
  "    man$input$input_path, man$input$input_format,",
  "    man$outputs$outputs_root, man$outputs$output_dir,",
  "    as.character(man$parameters$do_residualize),",
  "    as.character(man$parameters$residualize_used),",
  "    man$parameters$mad_scale_mode",
  "  )",
  ")",
  "knitr::kable(meta_tbl)",
  "```",
  "",
  "## Script header",
  "",
  "```text",
  header_text,
  "```",
  "",
  "## Dependencies",
  "",
  "```{r}",
  "deps <- dplyr::bind_rows(man$dependencies)",
  "knitr::kable(deps)",
  "```",
  "",
  "## Analytical logic",
  "",
  "This report summarizes robust within-timepoint dispersion using **MAD** across biological replicates, aggregated across features.",
  "",
  "**Feature-level dispersion (within timepoint):**",
  "",
  "$$\\mathrm{MAD}_{i,t} = \\mathrm{median}_r\\left(\\left|x_{i,r,t} - \\mathrm{median}_r(x_{i,r,t})\\right|\\right)$$",
  "",
  "**Timepoint summary:**",
  "",
  "$$D_t = \\mathrm{median}_i\\left(\\mathrm{MAD}_{i,t}\\right)$$",
  "",
  "If residualization is enabled, we remove the global mean–variance trend by fitting a smooth model on all feature×timepoint points:",
  "",
  "$$\\log(1+\\mathrm{MAD}_{i,t}) = f\\left(\\log(1+|\\mu_{i,t}|)\\right) + \\varepsilon_{i,t}$$",
  "",
  "and define:",
  "",
  "$$\\mathrm{logMAD\\_resid}_{i,t} = \\log(1+\\mathrm{MAD}_{i,t}) - \\hat{f}(\\log(1+|\\mu_{i,t}|))$$",
  "",
  "The primary reported trajectory is:",
  "",
  "$$D^{(resid)}_t = \\mathrm{median}_i\\left(\\mathrm{logMAD\\_resid}_{i,t}\\right)$$",
  "",
  "## Results overview",
  "",
  "```{r}",
  "by_day <- readr::read_csv('MAD_ByDay_Summary.csv', show_col_types = FALSE)",
  "knitr::kable(by_day)",
  "```",
  "",
  "## Plots",
  "",
  "```{r}",
  "knitr::include_graphics('MAD_ByDay_Dt.png')",
  "```",
  "",
  "```{r}",
  "knitr::include_graphics('MAD_ByFeatureDay_LogMeanVsLogMAD.png')",
  "```",
  "",
  "```{r}",
  "knitr::include_graphics('MAD_ByFeatureDay_Resid_Distribution.png')",
  "```",
  "",
  "## Interpretation guidance",
  "",
  "- **Higher** $D^{(resid)}_t$ indicates a broader accessible state space (more between-individual variability at that timepoint), after accounting for mean-dependent effects.",
  "- **Lower** $D^{(resid)}_t$ indicates tighter constraint (a narrower set of tolerated molecular configurations).",
  "- Interpret in conjunction with coordination/structure metrics (e.g., correlation architecture, effective dimensionality) if available; dispersion alone does not identify mechanism.",
  "",
  "## Generated files",
  "",
  "```{r}",
  "inv <- readr::read_csv('Project_Manifest_Files.csv', show_col_types = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd, con = qmd_path)

# ----------------------------- Render report (MANDATORY final step) -----------------------------
render_ok <- FALSE
render_msg <- NULL

# Render into output_dir; set wd to output_dir so relative paths work
old_wd <- getwd()
setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)

tryCatch({
  rmarkdown::render(
    input = basename(qmd_path),
    output_file = basename(html_out),
    output_dir = output_dir,
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
  render_ok <- file.exists(html_out)
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# Refresh inventory AFTER rendering (MANDATORY)
write_inventory(output_dir, inventory_path)

# Console confirmation (MANDATORY)
cat("\n==================== Report Status =======================\n")
if (render_ok) {
  cat("HTML report created:", normalizePath(html_out, winslash = "/", mustWork = FALSE), "\n")
} else {
  cat("FAILED to render HTML report.\n")
  cat("QMD path:", normalizePath(qmd_path, winslash = "/", mustWork = FALSE), "\n")
  if (!is.null(render_msg)) cat("Error:", render_msg, "\n")
}
cat("=========================================================\n\n")

cat("Done.\n")
