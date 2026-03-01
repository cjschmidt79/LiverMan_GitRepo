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
#   All outputs are written to: outputs/<run_name>_<timestamp>/
#     - MAD_ByFeatureDay.csv
#     - MAD_ByDay_Summary.csv
#     - Plots (PNG)  [base R]
#     - Project_Manifest.json
#     - Project_Manifest_Files.csv
#     - OUTPUT_INVENTORY.txt
#     - <script_name>_Report.qmd and rendered HTML report (if quarto available)
#
# NOTES
#   - Designed to be robust to outliers and mean dependence (via residualization).
#   - Keeps computations orthodox and transparent; avoids overclaiming.
#
# CONTRACT COMPLIANCE (minimal-change update)
#   - Base R only; NO dplyr/tidyverse. If dplyr loaded -> stop.
#   - Optional packages: rstudioapi (pickers), quarto (render). Script runs without them.
#   - Output dir policy: ./outputs/<run_name>_<timestamp>/ under git repo root if available, else getwd().
#   - Capture provenance early: script_name/path/full + git metadata.
#   - Build QMD as ONE character vector + ONE writeLines().
#   - Quarto render via setwd(output_dir) + quarto_render(basename(qmd_path)); never pass output_dir=.
# ======================================================================

# ----------------------------- Hard stop if dplyr is loaded -----------------------------
if ("package:dplyr" %in% search()) {
  stop("This script must not be run with dplyr loaded.", call. = FALSE)
}

# ----------------------------- Optional packages (allowed) -----------------------------
has_rstudioapi <- requireNamespace("rstudioapi", quietly = TRUE)
has_quarto     <- requireNamespace("quarto", quietly = TRUE)

# ----------------------------- Utilities -----------------------------
stop2 <- function(...) stop(paste0(...), call. = FALSE)
now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

safe_slug <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)   # non-alnum -> _
  x <- gsub("^_+|_+$", "", x)          # trim underscores
  x <- substr(x, 1, 80)               # keep reasonable length
  if (!nzchar(x)) x <- "run"
  x
}

is_interactive_rstudio <- function() {
  interactive() && has_rstudioapi && rstudioapi::isAvailable()
}

pick_file <- function(prompt = "Select a file") {
  cat("\n", prompt, "\n", sep = "")
  if (is_interactive_rstudio()) {
    p <- rstudioapi::selectFile(caption = prompt)
    if (is.null(p) || !nzchar(p)) stop2("No file selected.")
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  } else if (interactive()) {
    p <- file.choose()
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  } else {
    stop2("Non-interactive mode requires --input=/path/to/file.csv")
  }
}

coerce_day_numeric <- function(x) {
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
  # try "contains" matches
  hit2 <- which(vapply(cand_low, function(cc) any(grepl(cc, nms_low, fixed = TRUE)), logical(1)))
  if (length(hit2) > 0) {
    cc <- cand_low[hit2[1]]
    idx <- which(grepl(cc, nms_low, fixed = TRUE))[1]
    return(nms[idx])
  }
  NA_character_
}

# Robust MAD (unscaled) and scaled (sigma-approx)
mad_unscaled <- function(x) stats::mad(x, center = stats::median(x, na.rm = TRUE),
                                       constant = 1, na.rm = TRUE)
mad_scaled <- function(x) stats::mad(x, center = stats::median(x, na.rm = TRUE),
                                     constant = 1.4826, na.rm = TRUE)

# Git helpers (base R only)
git_safe <- function(cmd) {
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  if (length(out) == 0 || !nzchar(out[1])) return(NA_character_)
  out[1]
}

get_repoish_root <- function() {
  out <- tryCatch(system("git rev-parse --show-toplevel", intern = TRUE), error = function(e) character(0))
  if (length(out) == 0 || !nzchar(out[1])) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  normalizePath(out[1], winslash = "/", mustWork = TRUE)
}

# ----------------------------- Script identity capture (MANDATORY) -----------------------------
resolve_script_path <- function() {
  # 1) RStudio active document
  p <- tryCatch({
    if (has_rstudioapi && rstudioapi::isAvailable()) {
      ctx <- rstudioapi::getActiveDocumentContext()
      ctx$path
    } else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  
  # 2) Rscript --file=
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  
  # 3) source() fallback (sometimes present)
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  
  NA_character_
}

# If the script filename is known, set it here as a fallback:
known_script_filename <- "RobustMAD_Analysis.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

git_meta <- list(
  repo_root  = git_safe("git rev-parse --show-toplevel"),
  git_branch = git_safe("git rev-parse --abbrev-ref HEAD"),
  git_commit = git_safe("git rev-parse HEAD"),
  git_remote = git_safe("git config --get remote.origin.url")
)

cat("\n==================== Script Identity ====================\n")
cat("script_name:", script_name, "\n")
cat("script_path:", script_path, "\n")
cat("script_full:", ifelse(is.na(script_full), "NA", script_full), "\n")
if (is.na(script_full)) cat("NOTE: Script full path detection failed; using fallback filename:", known_script_filename, "\n")
cat("repo_root:",  git_meta$repo_root,  "\n")
cat("git_branch:", git_meta$git_branch, "\n")
cat("git_commit:", git_meta$git_commit, "\n")
cat("git_remote:", git_meta$git_remote, "\n")
cat("=========================================================\n\n")

# ----------------------------- CLI args -----------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag) {
  # expects --flag=value
  hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^--", flag, "="), "", hit[1])
}

input_path_cli    <- parse_arg("input")
analysis_name_cli <- parse_arg("analysis_name")
residualize_cli   <- parse_arg("residualize")
mad_scale_cli     <- parse_arg("mad_scale") # "unscaled" or "scaled"
run_name_cli      <- parse_arg("run_name")  # optional convenience

analysis_name <- if (!is.na(analysis_name_cli) && nzchar(analysis_name_cli)) analysis_name_cli else "MADDispersion"
do_residualize <- if (!is.na(residualize_cli) && nzchar(residualize_cli)) {
  tolower(residualize_cli) %in% c("true", "t", "1", "yes", "y")
} else TRUE

mad_scale_mode <- if (!is.na(mad_scale_cli) && nzchar(mad_scale_cli)) tolower(mad_scale_cli) else "unscaled"
if (!mad_scale_mode %in% c("unscaled", "scaled")) mad_scale_mode <- "unscaled"
mad_fun <- if (mad_scale_mode == "scaled") mad_scaled else mad_unscaled

# ----------------------------- Output directory policy (MANDATORY) -----------------------------
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
timestamp_tag <- now_stamp()

root_dir <- get_repoish_root()  # git root if available; else getwd()
outputs_root <- normalizePath(file.path(root_dir, "outputs"), winslash = "/", mustWork = FALSE)
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

# Ask user for run name (contract), with CLI override
default_run <- paste0(analysis_name, "_", timestamp_tag)
if (!is.na(run_name_cli) && nzchar(run_name_cli)) {
  run_name <- safe_slug(run_name_cli)
} else {
  if (interactive()) {
    cat("\nEnter a run name (e.g., MAD_Dispersion_Transcriptome, Day14_Check, Rev1_ExtFig3, etc.)\n")
    rn <- readline(paste0("run_name [default ", default_run, "]: "))
    rn <- trimws(rn)
    if (!nzchar(rn)) rn <- default_run
    run_name <- safe_slug(rn)
  } else {
    # Non-interactive: use deterministic default
    run_name <- safe_slug(default_run)
  }
}

run_folder <- paste0(run_name, "_", timestamp_tag)
output_dir <- file.path(outputs_root, run_folder)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

cat("\n==================== Run Configuration ===================\n")
cat("run_name:", run_name, "\n")
cat("timestamp:", timestamp, "\n")
cat("do_residualize:", do_residualize, "\n")
cat("mad_scale_mode:", mad_scale_mode, "\n")
cat("output_dir:", output_dir, "\n")
cat("=========================================================\n\n")

# ----------------------------- Input selection (interactive + CLI) -----------------------------
choose_input <- function() {
  if (!is.na(input_path_cli) && nzchar(input_path_cli)) {
    p <- input_path_cli
    if (!file.exists(p)) stop2("Input file does not exist: ", p)
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }
  pick_file("[Input] Select a CSV file (either tidy-long or wide format).")
}

# input_paths contract variable
input_paths <- character(0)

input_path <- choose_input()
input_paths <- c(input_paths, input_path)

# ----------------------------- Read input (agnostic, base R) -----------------------------
df_raw <- tryCatch(read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE),
                   error = function(e) NULL)
if (is.null(df_raw) || nrow(df_raw) == 0) stop2("Input file appears empty or unreadable as CSV: ", input_path)

nms <- names(df_raw)

day_col <- find_col(nms, c("day", "time", "timepoint", "tp"))
sample_col <- find_col(nms, c("sampleid", "sample", "replicate", "bird", "id", "subject"))
feature_col <- find_col(nms, c("feature", "gene", "metabolite", "analyte", "compound"))
abundance_col <- find_col(nms, c("abundance", "value", "count", "counts", "intensity", "area", "signal"))

has_tidy_signature <- !is.na(day_col) && !is.na(feature_col) && !is.na(abundance_col)
input_format <- if (has_tidy_signature) "TidyLong" else "Wide"
source_label <- input_format

cat("input_path:", input_path, "\n")
cat("input_format:", input_format, "\n\n")

# ----------------------------- Canonicalize to tidy long: Day, SampleID, Feature, Value -----------------------------
if (input_format == "TidyLong") {
  
  if (is.na(sample_col)) {
    sample_col <- ".SampleID"
    # replicate index within day (stable)
    df_raw[[sample_col]] <- ave(seq_len(nrow(df_raw)), df_raw[[day_col]], FUN = seq_along)
  }
  
  Day      <- coerce_day_numeric(df_raw[[day_col]])
  SampleID <- as.character(df_raw[[sample_col]])
  Feature  <- as.character(df_raw[[feature_col]])
  Value    <- suppressWarnings(as.numeric(df_raw[[abundance_col]]))
  
  keep <- !is.na(Day) & !is.na(SampleID) & nzchar(SampleID) & !is.na(Feature) & nzchar(Feature) & !is.na(Value)
  df_long <- data.frame(Day = Day[keep], SampleID = SampleID[keep], Feature = Feature[keep], Value = Value[keep],
                        stringsAsFactors = FALSE)
  
} else {
  
  if (is.na(day_col)) {
    stop2("Wide format detected but could not find a Day/Time column. Columns are: ",
          paste(names(df_raw), collapse = ", "))
  }
  
  if (is.na(sample_col)) {
    sample_col <- ".SampleID"
    df_raw[[sample_col]] <- paste0("S", seq_len(nrow(df_raw)))
  }
  
  feature_cols <- setdiff(names(df_raw), c(day_col, sample_col))
  if (length(feature_cols) < 2) stop2("Wide format requires at least 2 feature columns besides Day/SampleID.")
  
  Day      <- coerce_day_numeric(df_raw[[day_col]])
  SampleID <- as.character(df_raw[[sample_col]])
  
  # Build long by stacking columns (base R)
  n_feat <- length(feature_cols)
  n_row  <- nrow(df_raw)
  
  Day_rep      <- rep(Day, times = n_feat)
  SampleID_rep <- rep(SampleID, times = n_feat)
  Feature_rep  <- rep(feature_cols, each = n_row)
  
  # c() column-wise values in feature_cols order
  vals <- suppressWarnings(as.numeric(unlist(df_raw[feature_cols], use.names = FALSE)))
  
  keep <- !is.na(Day_rep) & !is.na(SampleID_rep) & nzchar(SampleID_rep) & !is.na(Feature_rep) & nzchar(Feature_rep) & !is.na(vals)
  
  df_long <- data.frame(
    Day = Day_rep[keep],
    SampleID = SampleID_rep[keep],
    Feature = Feature_rep[keep],
    Value = vals[keep],
    stringsAsFactors = FALSE
  )
}

if (nrow(df_long) == 0) stop2("No valid (Day, SampleID, Feature, Value) rows after parsing.")

# Contract: keep input_paths normalized absolute
input_paths <- normalizePath(input_paths, winslash = "/", mustWork = TRUE)

# Basic guardrails
n_days <- length(unique(df_long$Day))
n_features <- length(unique(df_long$Feature))
if (n_days < 2) stop2("Need at least 2 unique timepoints (Day). Found: ", n_days)
if (n_features < 10) cat("NOTE: Fewer than 10 features detected (", n_features, "). Analysis will run but interpret cautiously.\n", sep = "")

# ----------------------------- Compute feature-by-day statistics -----------------------------
# Aggregate by Day+Feature
key <- paste(df_long$Day, df_long$Feature, sep = "||")
groups <- split(seq_len(nrow(df_long)), key)

by_fd_list <- lapply(groups, function(idx) {
  d <- df_long$Day[idx][1]
  f <- df_long$Feature[idx][1]
  reps <- df_long$SampleID[idx]
  x <- df_long$Value[idx]
  
  n_reps <- length(unique(reps))
  mean_value <- mean(x, na.rm = TRUE)
  mad_value  <- mad_fun(x)
  
  data.frame(
    Day = d,
    Feature = f,
    n_reps = n_reps,
    mean_value = mean_value,
    mad_value = mad_value,
    stringsAsFactors = FALSE
  )
})

by_fd <- do.call(rbind, by_fd_list)
if (is.null(by_fd) || nrow(by_fd) == 0) stop2("Failed to compute feature-by-day table.")

# Filter finite
keep_fd <- is.finite(by_fd$mean_value) & is.finite(by_fd$mad_value)
by_fd <- by_fd[keep_fd, , drop = FALSE]

by_fd$logMean <- log1p(abs(by_fd$mean_value))
by_fd$logMAD  <- log1p(abs(by_fd$mad_value))

# ----------------------------- Residualize mean–variance trend (optional) -----------------------------
residualize_used <- FALSE
by_fd$logMAD_resid <- by_fd$logMAD

if (do_residualize) {
  lo <- tryCatch(
    stats::loess(logMAD ~ logMean, data = by_fd, span = 0.75, degree = 2,
                 control = stats::loess.control(surface = "direct")),
    error = function(e) NULL
  )
  if (is.null(lo)) {
    cat("WARNING: loess fit failed; proceeding without residualization.\n")
    residualize_used <- FALSE
  } else {
    pred <- tryCatch(stats::predict(lo, newdata = by_fd),
                     error = function(e) rep(NA_real_, nrow(by_fd)))
    if (all(is.na(pred))) {
      cat("WARNING: loess prediction failed; proceeding without residualization.\n")
      residualize_used <- FALSE
    } else {
      by_fd$logMAD_resid <- by_fd$logMAD - pred
      residualize_used <- TRUE
    }
  }
}

# ----------------------------- Timepoint summaries: D_t and D_t_resid -----------------------------
days <- sort(unique(by_fd$Day))
by_day <- data.frame(
  Day = days,
  n_features = NA_integer_,
  med_n_reps = NA_real_,
  D_t_mad = NA_real_,
  D_t_logMAD = NA_real_,
  D_t_logMAD_resid = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(days)) {
  d <- days[i]
  sub <- by_fd[by_fd$Day == d, , drop = FALSE]
  by_day$n_features[i] <- length(unique(sub$Feature))
  by_day$med_n_reps[i] <- stats::median(sub$n_reps, na.rm = TRUE)
  by_day$D_t_mad[i] <- stats::median(sub$mad_value, na.rm = TRUE)
  by_day$D_t_logMAD[i] <- stats::median(sub$logMAD, na.rm = TRUE)
  by_day$D_t_logMAD_resid[i] <- stats::median(sub$logMAD_resid, na.rm = TRUE)
}

# ----------------------------- Write primary outputs -----------------------------
out_by_fd  <- file.path(output_dir, "MAD_ByFeatureDay.csv")
out_by_day <- file.path(output_dir, "MAD_ByDay_Summary.csv")
write.csv(by_fd,  out_by_fd,  row.names = FALSE)
write.csv(by_day, out_by_day, row.names = FALSE)

# ----------------------------- Plots (base R; PNG) -----------------------------
save_png <- function(path, width = 2000, height = 1400, res = 250, expr) {
  png(filename = path, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  expr
}

plot1_path <- file.path(output_dir, "MAD_ByDay_Dt.png")
save_png(plot1_path, width = 2000, height = 1200, res = 250, expr = {
  plot(by_day$Day, by_day$D_t_logMAD_resid, type = "b",
       xlab = "Timepoint (Day)",
       ylab = "D_t = median_i(logMAD_resid)",
       main = "Robust Dispersion by Timepoint (D_t from logMAD_resid)")
  grid()
})

plot2_path <- file.path(output_dir, "MAD_ByFeatureDay_LogMeanVsLogMAD.png")
save_png(plot2_path, width = 1800, height = 1400, res = 250, expr = {
  plot(by_fd$logMean, by_fd$logMAD, pch = 16, cex = 0.5,
       xlab = "log1p(|mean|)",
       ylab = "log1p(|MAD|)",
       main = "Mean–Variance Structure (log1p|mean| vs log1p|MAD|)")
  grid()
})

plot3_path <- file.path(output_dir, "MAD_ByFeatureDay_Resid_Distribution.png")
save_png(plot3_path, width = 1800, height = 1400, res = 250, expr = {
  hist(by_fd$logMAD_resid, breaks = 60,
       main = "Distribution of logMAD_resid Across Feature×Day",
       xlab = "logMAD_resid")
  grid()
})

# ----------------------------- Manifest construction (kept; base R JSON writer) -----------------------------
# Minimal JSON writer (sufficient for this manifest; avoids jsonlite)
json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\n", "\\\\n", x)
  x
}

json_kv <- function(k, v, indent = 2) {
  pad <- paste(rep(" ", indent), collapse = "")
  if (is.na(v)) return(paste0(pad, "\"", k, "\": null"))
  if (is.logical(v)) return(paste0(pad, "\"", k, "\": ", ifelse(v, "true", "false")))
  if (is.numeric(v)) return(paste0(pad, "\"", k, "\": ", v))
  paste0(pad, "\"", k, "\": \"", json_escape(as.character(v)), "\"")
}

json_obj <- function(lines, indent = 0) {
  pad <- paste(rep(" ", indent), collapse = "")
  c(paste0(pad, "{"),
    paste0(lines, collapse = ",\n"),
    paste0(pad, "}"))
}

run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

manifest <- c(
  paste0("  \"run_id\": \"", json_escape(paste0(run_name, "_", source_label, "_", timestamp_tag)), "\","),
  paste0("  \"run_timestamp\": \"", json_escape(run_timestamp), "\","),
  "  \"script\": {",
  paste0("    \"name\": \"", json_escape(script_name), "\","),
  paste0("    \"path\": \"", json_escape(script_path), "\","),
  paste0("    \"full_path\": ", if (is.na(script_full)) "null" else paste0("\"", json_escape(script_full), "\"")),
  "  },",
  "  \"git\": {",
  paste0("    \"repo_root\": ", if (is.na(git_meta$repo_root)) "null" else paste0("\"", json_escape(git_meta$repo_root), "\""), ","),
  paste0("    \"git_branch\": ", if (is.na(git_meta$git_branch)) "null" else paste0("\"", json_escape(git_meta$git_branch), "\""), ","),
  paste0("    \"git_commit\": ", if (is.na(git_meta$git_commit)) "null" else paste0("\"", json_escape(git_meta$git_commit), "\""), ","),
  paste0("    \"git_remote\": ", if (is.na(git_meta$git_remote)) "null" else paste0("\"", json_escape(git_meta$git_remote), "\"")),
  "  },",
  "  \"input\": {",
  paste0("    \"input_path\": \"", json_escape(input_path), "\","),
  paste0("    \"input_format\": \"", json_escape(input_format), "\""),
  "  },",
  "  \"parameters\": {",
  paste0("    \"analysis_name\": \"", json_escape(analysis_name), "\","),
  paste0("    \"do_residualize\": ", ifelse(do_residualize, "true", "false"), ","),
  paste0("    \"residualize_used\": ", ifelse(residualize_used, "true", "false"), ","),
  paste0("    \"mad_scale_mode\": \"", json_escape(mad_scale_mode), "\","),
  paste0("    \"mad_definition\": \"", json_escape(if (mad_scale_mode == "scaled")
    "stats::mad(constant=1.4826): sigma-approx scaled MAD"
    else
      "stats::mad(constant=1): unscaled MAD"), "\""),
  "  },",
  "  \"outputs\": {",
  paste0("    \"outputs_root\": \"", json_escape(outputs_root), "\","),
  paste0("    \"output_dir\": \"", json_escape(output_dir), "\""),
  "  }"
)

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
writeLines(c("{", manifest, "}"), con = manifest_json_path)

# ----------------------------- File inventory CSV + OUTPUT_INVENTORY.txt (MANDATORY) -----------------------------
inventory_path <- file.path(output_dir, "Project_Manifest_Files.csv")
output_inventory_txt <- file.path(output_dir, "OUTPUT_INVENTORY.txt")

write_inventory <- function(dir_path, out_csv, out_txt) {
  files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  
  # TXT inventory (contract)
  writeLines(files, con = out_txt)
  
  # CSV inventory (kept)
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
    dir_norm <- normalizePath(dir_path, winslash = "/", mustWork = FALSE)
    files_norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
    rel <- gsub(paste0("^", dir_norm, "/?"), "", files_norm)
    
    inv <- data.frame(
      file = files_norm,
      rel_path = rel,
      bytes = as.numeric(info$size),
      modified = as.character(info$mtime),
      stringsAsFactors = FALSE
    )
    inv <- inv[order(inv$rel_path), , drop = FALSE]
  }
  
  write.csv(inv, out_csv, row.names = FALSE)
  inv
}

write_inventory(output_dir, inventory_path, output_inventory_txt)

# ----------------------------- QMD report (MANDATORY structure; Quarto if available) -----------------------------
extract_header_block <- function(script_full_path) {
  if (!is.na(script_full_path) && file.exists(script_full_path)) {
    lines <- readLines(script_full_path, warn = FALSE)
    out <- character(0)
    for (ln in lines) {
      if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) out <- c(out, ln) else break
    }
    return(paste(out, collapse = "\n"))
  }
  paste(
    "# Robust MAD Dispersion (MAD_resid) Analysis",
    "# (Header unavailable: script_full path not detected in this mode.)",
    sep = "\n"
  )
}

header_text <- extract_header_block(script_full)

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))

# Build QMD content as ONE character vector + ONE writeLines() (contract)
qmd <- c(
  "---",
  paste0("title: \"", script_name, " — Robust MAD Dispersion Report\""),
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 2",
  "execute:",
  "  echo: true",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "## Summary",
  "This report summarizes robust within-timepoint dispersion using MAD across replicates, aggregated across features.",
  "",
  "## Script header (verbatim)",
  "```",
  header_text,
  "```",
  "",
  "## Metadata",
  paste0("- **timestamp:** ", timestamp),
  paste0("- **script_name:** ", script_name),
  paste0("- **script_path:** ", script_path),
  paste0("- **script_full:** ", ifelse(is.na(script_full), "NA", script_full)),
  paste0("- **output_dir:** ", output_dir),
  paste0("- **repo_root:** ", git_meta$repo_root),
  paste0("- **git_branch:** ", git_meta$git_branch),
  paste0("- **git_commit:** ", git_meta$git_commit),
  paste0("- **git_remote:** ", git_meta$git_remote),
  "",
  "### Input files (absolute paths)",
  "```",
  paste(input_paths, collapse = "\n"),
  "```",
  "",
  "### Parameters",
  "```",
  paste(
    c(
      paste0("analysis_name: ", analysis_name),
      paste0("input_format: ", input_format),
      paste0("do_residualize: ", do_residualize),
      paste0("residualize_used: ", residualize_used),
      paste0("mad_scale_mode: ", mad_scale_mode)
    ),
    collapse = "\n"
  ),
  "```",
  "",
  "## Analytical logic",
  "",
  "Feature-level dispersion (within timepoint):",
  "",
  "$$\\mathrm{MAD}_{i,t} = \\mathrm{median}_r\\left(\\left|x_{i,r,t} - \\mathrm{median}_r(x_{i,r,t})\\right|\\right)$$",
  "",
  "Timepoint summary:",
  "",
  "$$D_t = \\mathrm{median}_i\\left(\\mathrm{MAD}_{i,t}\\right)$$",
  "",
  "If residualization is enabled, we fit a smooth mean–variance trend across all feature×timepoint points:",
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
  "## Results table",
  "",
  "```{r}",
  "by_day <- read.csv('MAD_ByDay_Summary.csv', check.names = FALSE, stringsAsFactors = FALSE)",
  "by_day",
  "```",
  "",
  "## Plots",
  "",
  "### D_t trajectory",
  paste0("![](", basename(plot1_path), ")"),
  "",
  "### Mean–variance structure",
  paste0("![](", basename(plot2_path), ")"),
  "",
  "### Residual distribution",
  paste0("![](", basename(plot3_path), ")"),
  "",
  "## Generated outputs",
  "```{r}",
  "list.files('.', recursive = TRUE, full.names = TRUE)",
  "```",
  "",
  "## Session info",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd, con = qmd_path)

# ----------------------------- Render report (Quarto contract) -----------------------------
render_ok <- FALSE
render_msg <- NULL
html_path <- sub("\\.qmd$", ".html", qmd_path)

if (!has_quarto) {
  cat("\nNOTE: quarto package not installed. Skipping HTML render.\n")
  cat("QMD written to: ", qmd_path, "\n", sep = "")
} else {
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  tryCatch({
    quarto::quarto_render(basename(qmd_path))
    render_ok <- file.exists(html_path)
  }, error = function(e) {
    render_ok <<- FALSE
    render_msg <<- conditionMessage(e)
  })
}

# Refresh inventories AFTER rendering (MANDATORY)
write_inventory(output_dir, inventory_path, output_inventory_txt)

cat("\n==================== Report Status =======================\n")
if (has_quarto && render_ok) {
  cat("HTML report created:", normalizePath(html_path, winslash = "/", mustWork = FALSE), "\n")
} else if (has_quarto && !render_ok) {
  cat("FAILED to render HTML report.\n")
  cat("QMD path:", normalizePath(qmd_path, winslash = "/", mustWork = FALSE), "\n")
  if (!is.null(render_msg)) cat("Error:", render_msg, "\n")
} else {
  cat("HTML render skipped (quarto not available).\n")
}
cat("OUTPUT_INVENTORY.txt:", normalizePath(output_inventory_txt, winslash = "/", mustWork = FALSE), "\n")
cat("=========================================================\n\n")

cat("Done.\n")
cat("Output directory:\n", output_dir, "\n", sep = "")

# ----------------------------- Auto-run when sourced (contract) -----------------------------
# This script is already "mainline" (not wrapped). The following is a no-op placeholder
# for future refactors; leaving it minimal to avoid changing working behavior.
if (interactive()) {
  # Optional: View key output table for biologists
  # (Will not error if View is unavailable in non-RStudio contexts)
  try(View(by_day), silent = TRUE)
}