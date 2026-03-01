#!/usr/bin/env Rscript
# ============================================================
# Two-panel Figure (Mean + Dispersion) from Tidy Long Data
#   - NO hard-coded column names
#   - Auto-detect Day / SampleID / Feature / Abundance / Source (optional)
#   - If Source exists and split_by_source=1: analyze each Source independently
#
# Panel A: system mean with replicate dispersion band (across SampleID within Day)
# Panel B: system-level dispersion across features (within Day) + optional 95% CI
#
# Expected data format (generic tidy long):
#   Day | SampleID | Feature (Gene/Metabolite/etc) | Abundance
# Optional:
#   Source
#
# Usage (CLI):
#   Rscript make_two_panel_figure.R --in="path/to/file.csv"
#
# Optional flags:
#   --band="iqr"                 # iqr or quantiles
#   --q_lo=0.25 --q_hi=0.75      # if band=quantiles
#   --disp_metric="sd"           # sd, cv, iqr, mad  (per feature within day)
#   --disp_agg="median"          # median or mean aggregation across features
#
#   TRANSFORM (NEW):
#   --transform="raw"            # raw | log2p | log1p
#   --pseudocount=1              # numeric; used for log2p (default 1)
#
#   Backward-compatibility:
#   --log_abundance=0            # legacy: if TRUE and --transform not provided => log1p
#
#   --pdf=0                      # 1 to also write PDF
#
# Source handling:
#   --split_by_source=1          # if Source exists: analyze each Source separately (default auto)
#   --include_pooled=0           # also run pooled analysis across all sources (default 0)
#   --source_regex=""            # optional regex filter on Source values (case-insensitive)
#
# Confidence interval options (Panel B):
#   --ci=1                       # 1 to compute CI for Panel B (default: 1)
#   --ci_level=0.95              # CI level (default 0.95)
#   --ci_mode="sample"           # "sample": resample SampleID within Day (within Source)
#                                # "feature": resample feature dispersions within Day (within Source)
#   --n_boot=500                 # number of bootstrap replicates (if omitted, script will ask interactively)
#   --seed=1                     # random seed
#
# Column overrides (only if auto-detection fails):
#   --col_day= --col_sample= --col_feature= --col_abundance= --col_source=
#
# ============================================================
# RETROFIT LAYER (Standardized outputs + manifest + self-report)
# - All outputs go to: outputs/<run_id>/
# - Writes:
#     Project_Manifest.json
#     Project_Manifest_Files.csv
#     <script_name>_Report.qmd
#     <script_name>_Report.html (attempted)
# - Optional features (internal flags):
#     enable_plotly <- TRUE
#     enable_runtime_tracking <- TRUE
# ============================================================

# --------------------------
# Internal flags (MANDATORY)
# --------------------------
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

# Runtime tracking start (MANDATORY if enabled)
start_time <- Sys.time()

# --------------------------
# Dependency handling (quiet install)
# --------------------------
quiet_install_if_missing <- function(pkgs) {
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) {
    suppressWarnings(utils::install.packages(to_install, repos = "https://cloud.r-project.org", quiet = TRUE))
  }
  invisible(TRUE)
}

base_pkgs <- c("data.table", "ggplot2", "patchwork", "rmarkdown", "knitr", "jsonlite")
if (enable_plotly) base_pkgs <- c(base_pkgs, "plotly")
quiet_install_if_missing(base_pkgs)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --------------------------
# Script identity capture (MANDATORY block)
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
known_script_filename <- "make_two_panel_figure.R"
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
  v %in% c("1","true","t","yes","y")
}
stop2 <- function(...) stop(paste0(...), call. = FALSE)

safe_slug <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "NA", x)
}

# --------------------------
# Inputs (analysis settings)
# --------------------------
in_path <- arg_val("in", NULL)

band_mode <- tolower(arg_val("band", "iqr"))                 # iqr | quantiles
q_lo_band <- as.numeric(arg_val("q_lo", "0.25"))
q_hi_band <- as.numeric(arg_val("q_hi", "0.75"))

disp_metric <- tolower(arg_val("disp_metric", "sd"))         # sd | cv | iqr | mad
disp_agg <- tolower(arg_val("disp_agg", "median"))           # median | mean

# Legacy transform flag (kept for backwards compatibility)
log_abundance <- arg_flag("log_abundance", FALSE)            # legacy: log1p
write_pdf <- arg_flag("pdf", FALSE)

# NEW transform flags
transform_raw <- arg_val("transform", NULL)                  # raw | log2p | log1p (NULL => prompt or fallback)
pseudocount <- suppressWarnings(as.numeric(arg_val("pseudocount", "1")))

# Source handling
split_by_source_raw <- arg_val("split_by_source", NULL)      # may be NULL -> auto
include_pooled <- arg_flag("include_pooled", FALSE)
source_regex <- arg_val("source_regex", "")

# CI options
ci_on    <- arg_flag("ci", TRUE)
ci_level <- as.numeric(arg_val("ci_level", "0.95"))
ci_mode  <- tolower(arg_val("ci_mode", "sample"))            # sample | feature
seed     <- as.integer(arg_val("seed", "1"))

# n_boot: if missing, ask the user (interactive) with default suggestion
n_boot_raw <- arg_val("n_boot", NULL)
suggested_default_boot <- 500L

# Optional explicit column overrides (only needed if detection fails)
col_day <- arg_val("col_day", NULL)
col_sample <- arg_val("col_sample", NULL)
col_feature <- arg_val("col_feature", NULL)
col_abundance <- arg_val("col_abundance", NULL)
col_source <- arg_val("col_source", NULL)

# Analysis name for run_id (inferred; can be changed here without touching CLI)
analysis_name <- "TwoPanelMeanDispersion"

# Validate CI inputs
if (!(ci_level > 0 && ci_level < 1)) stop2("--ci_level must be between 0 and 1.")
if (!(ci_mode %in% c("sample","feature"))) stop2("--ci_mode must be 'sample' or 'feature'.")

# Validate pseudocount if supplied/used
if (!is.finite(pseudocount) || is.na(pseudocount) || pseudocount <= 0) {
  pseudocount <- 1
}

# Resolve n_boot with prompt if needed
if (ci_on) {
  if (is.null(n_boot_raw) || nchar(n_boot_raw) == 0) {
    if (interactive()) {
      prompt <- paste0("How many bootstrap replicates for Panel B CI? [default ", suggested_default_boot, "]: ")
      ans <- trimws(readline(prompt))
      if (nchar(ans) == 0) {
        n_boot <- suggested_default_boot
      } else {
        n_boot <- suppressWarnings(as.integer(ans))
        if (!is.finite(n_boot) || is.na(n_boot) || n_boot < 10) {
          message("Invalid entry. Using default n_boot = ", suggested_default_boot)
          n_boot <- suggested_default_boot
        }
      }
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

set.seed(seed)

# --------------------------
# Resolve transform choice (NEW)
# --------------------------
normalize_transform <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  x <- tolower(trimws(x))
  if (x %in% c("raw","none","direct")) return("raw")
  if (x %in% c("log2p","log2pc","log2+pseudocount","log2_pseudocount","log2")) return("log2p")
  if (x %in% c("log1p","log1","ln1p")) return("log1p")
  x
}

transform_choice <- normalize_transform(transform_raw)

if (is.null(transform_choice)) {
  # If legacy flag is TRUE and user didn't specify new transform, preserve old behavior.
  if (isTRUE(log_abundance)) {
    transform_choice <- "log1p"
  } else if (interactive()) {
    message("\nTransform option:")
    message("  [0] Analyze Abundance directly (raw)")
    message("  [1] log2(Abundance + pseudocount)")
    ans <- trimws(readline("Choose 0 or 1 [default 0]: "))
    if (nchar(ans) == 0 || ans == "0") {
      transform_choice <- "raw"
    } else if (ans == "1") {
      transform_choice <- "log2p"
      ans_pc <- trimws(readline(paste0("Pseudocount for log2(Abundance + pseudocount)? [default ", pseudocount, "]: ")))
      if (nchar(ans_pc) > 0) {
        pc2 <- suppressWarnings(as.numeric(ans_pc))
        if (is.finite(pc2) && !is.na(pc2) && pc2 > 0) pseudocount <- pc2
      }
    } else {
      message("Unrecognized choice; defaulting to raw.")
      transform_choice <- "raw"
    }
  } else {
    transform_choice <- "raw"
  }
} else {
  if (!(transform_choice %in% c("raw","log2p","log1p"))) {
    stop2("--transform must be one of: raw | log2p | log1p (got: ", transform_raw, ")")
  }
}

transform_label <- switch(
  transform_choice,
  raw   = "raw",
  log1p = "log1p",
  log2p = paste0("log2p_pc", safe_slug(pseudocount)),
  "raw"
)

message("\nTransform selection:")
message("  transform: ", transform_choice)
if (transform_choice == "log2p") message("  pseudocount: ", pseudocount)

# --------------------------
# Input path
# --------------------------
if (is.null(in_path) || nchar(in_path) == 0) {
  if (interactive() && exists("file.choose")) {
    in_path <- file.choose()
  } else {
    stop2("No input provided. Run:\n  Rscript make_two_panel_figure.R --in=\"path/to/file.csv\"")
  }
}

# --------------------------
# Output directory policy (MANDATORY)
# --------------------------
outputs_root <- file.path(getwd(), "outputs")
dir.create(outputs_root, showWarnings = FALSE, recursive = TRUE)

run_timestamp <- Sys.time()
run_stamp <- format(run_timestamp, "%Y%m%d_%H%M%S")

run_id <- paste0(analysis_name, "_TEMP_", run_stamp)
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Console identity summary (MANDATORY)
message("Script identity:")
message("  script_name: ", script_name)
message("  script_path: ", script_path)
message("  script_full: ", ifelse(is.na(script_full), "NA", script_full))
if (is.na(script_full)) {
  message("NOTE: script path detection failed; using fallback known_script_filename = ", known_script_filename)
}

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
  return(nm0[hits])  # candidates if ambiguous
}
is_multi <- function(x) is.character(x) && length(x) > 1

cand_day <- pick_col(col_day, c("^day$", "day", "timepoint", "age"), "Day column")
cand_sample <- pick_col(col_sample, c("^sampleid$", "^sample$", "sample", "replicate", "bird", "id$"), "SampleID column")
cand_feature <- pick_col(col_feature, c("^gene$", "^metabolite$", "^feature$", "gene", "metab", "compound", "feature"), "Feature column")
cand_abund <- pick_col(col_abundance, c("^abundance$", "^value$", "^expr$", "abund", "count", "fpkm", "tpm", "intensity", "area", "value"), "Abundance column")
cand_source <- pick_col(col_source, c("^source$", "tissue", "compartment"), "Source column")

if (is_multi(cand_day) || is_multi(cand_sample) || is_multi(cand_feature) || is_multi(cand_abund)) {
  message("\nColumn auto-detection is ambiguous. Candidates:")
  if (is_multi(cand_day))     message("  Day candidates:       ", paste(cand_day, collapse = ", "))
  if (is_multi(cand_sample))  message("  Sample candidates:    ", paste(cand_sample, collapse = ", "))
  if (is_multi(cand_feature)) message("  Feature candidates:   ", paste(cand_feature, collapse = ", "))
  if (is_multi(cand_abund))   message("  Abundance candidates: ", paste(cand_abund, collapse = ", "))
  message("\nRe-run with explicit overrides, e.g.:")
  message("  Rscript make_two_panel_figure.R --in=\"...\" --col_day=\"Day\" --col_sample=\"SampleID\" --col_feature=\"Gene\" --col_abundance=\"Abundance\"")
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
# Standardize columns into a working dt
# --------------------------
dt0[, Day := as.numeric(get(col_day))]
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

# Optional filter on Source values
if (!is.null(col_source) && nchar(source_regex) > 0) {
  keep <- grepl(source_regex, dt0$Source, ignore.case = TRUE)
  if (!any(keep)) stop2("source_regex='", source_regex, "' matched 0 rows.")
  dt0 <- dt0[keep]
  message("Applied source_regex filter; rows now: ", nrow(dt0))
}

# Determine split_by_source default: if Source exists and has >1 unique non-NA value -> TRUE
has_source <- !is.null(col_source) && !all(is.na(dt0$Source))
n_sources <- if (has_source) length(unique(dt0[!is.na(Source), Source])) else 0L

if (is.null(split_by_source_raw)) {
  split_by_source <- has_source && n_sources > 1
} else {
  split_by_source <- tolower(trimws(split_by_source_raw)) %in% c("1","true","t","yes","y")
}

message("Source handling: ",
        if (!has_source) "no Source column detected"
        else paste0(n_sources, " source level(s) detected; split_by_source=", split_by_source,
                    if (include_pooled) " (also pooled)" else ""))

# --------------------------
# Finalize run_id now that Source mode is known (MANDATORY)
# --------------------------
source_tag_for_run <- if (!has_source) {
  "NOSOURCE"
} else if (split_by_source && n_sources > 1) {
  "MULTI"
} else if (n_sources == 1) {
  unique(dt0[!is.na(Source), Source])[1]
} else {
  "POOLED"
}

run_id_final <- paste0(
  analysis_name, "_",
  safe_slug(source_tag_for_run), "_",
  safe_slug(transform_label), "_",
  run_stamp
)

output_dir_final <- file.path(outputs_root, run_id_final)
dir.create(output_dir_final, showWarnings = FALSE, recursive = TRUE)
output_dir <- output_dir_final
run_id <- run_id_final

message("\nRun outputs:")
message("  outputs_root: ", outputs_root)
message("  run_id:       ", run_id)
message("  output_dir:   ", output_dir)

# --------------------------
# Transform helper (NEW)
# --------------------------
apply_transform <- function(dt, transform_choice, pseudocount) {
  out <- copy(dt)
  if (transform_choice == "raw") {
    return(out)
  }
  if (transform_choice == "log1p") {
    out[, Abundance := log1p(Abundance)]
    return(out)
  }
  if (transform_choice == "log2p") {
    min_val <- suppressWarnings(min(out$Abundance, na.rm = TRUE))
    if (!is.finite(min_val)) stop2("Abundance contains no finite values; cannot transform.")
    if (any((out$Abundance + pseudocount) <= 0, na.rm = TRUE)) {
      stop2("Cannot apply log2(Abundance + pseudocount): some values are <= 0 after adding pseudocount. ",
            "Check Abundance scale or increase --pseudocount. Current pseudocount=", pseudocount, ".")
    }
    out[, Abundance := log2(Abundance + pseudocount)]
    return(out)
  }
  stop2("Internal error: unknown transform_choice=", transform_choice)
}

# --------------------------
# Core functions for one dataset (one Source subset)
# --------------------------
disp_fun <- function(x, disp_metric_local = disp_metric) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  if (disp_metric_local == "sd")  return(sd(x))
  if (disp_metric_local == "cv")  { m <- mean(x); if (m == 0) return(NA_real_); return(sd(x) / abs(m)) }
  if (disp_metric_local == "iqr") return(IQR(x))
  if (disp_metric_local == "mad") return(mad(x, constant = 1))
  stop2("Unknown --disp_metric: ", disp_metric_local, " (sd, cv, iqr, mad)")
}
agg_fun <- function(x, disp_agg_local = disp_agg) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (disp_agg_local == "median") return(median(x))
  if (disp_agg_local == "mean")   return(mean(x))
  stop2("Unknown --disp_agg: ", disp_agg_local, " (median or mean)")
}

compute_two_panel <- function(dt, source_label, output_dir, file_tag,
                              transform_choice, transform_label, pseudocount) {
  
  dt_use <- apply_transform(dt, transform_choice = transform_choice, pseudocount = pseudocount)
  
  # ---- Panel A ----
  rep_means <- dt_use[, .(rep_mean = mean(Abundance, na.rm = TRUE)), by = .(Day, SampleID)]
  
  if (band_mode == "iqr") {
    A_df <- rep_means[, .(
      mean_val = mean(rep_mean, na.rm = TRUE),
      band_lo = as.numeric(quantile(rep_mean, 0.25, na.rm = TRUE)),
      band_hi = as.numeric(quantile(rep_mean, 0.75, na.rm = TRUE)),
      n_rep = .N
    ), by = Day][order(Day)]
    band_label <- "Replicate band (IQR)"
  } else if (band_mode == "quantiles") {
    if (!(q_lo_band < q_hi_band)) stop2("q_lo must be < q_hi.")
    A_df <- rep_means[, .(
      mean_val = mean(rep_mean, na.rm = TRUE),
      band_lo = as.numeric(quantile(rep_mean, q_lo_band, na.rm = TRUE)),
      band_hi = as.numeric(quantile(rep_mean, q_hi_band, na.rm = TRUE)),
      n_rep = .N
    ), by = Day][order(Day)]
    band_label <- paste0("Replicate band (q", q_lo_band, "–q", q_hi_band, ")")
  } else {
    stop2("Unknown --band mode: ", band_mode, " (use iqr or quantiles)")
  }
  
  # ---- Panel B observed ----
  feat_disp <- dt_use[, .(feat_disp = disp_fun(Abundance)), by = .(Day, Feature)]
  feat_disp <- feat_disp[is.finite(feat_disp)]
  if (nrow(feat_disp) == 0) stop2("No feature-level dispersion values computed for Source=", source_label)
  
  B_df <- feat_disp[, .(
    disp   = agg_fun(feat_disp),
    n_feat = sum(is.finite(feat_disp))
  ), by = Day][order(Day)]
  
  B_df[, `:=`(disp_lo = NA_real_, disp_hi = NA_real_, ci_mode = if (ci_on) ci_mode else NA_character_)]
  
  # ---- Panel B CI ----
  if (ci_on) {
    alpha <- (1 - ci_level) / 2
    q_lo_ci <- alpha
    q_hi_ci <- 1 - alpha
    days <- sort(unique(dt_use$Day))
    
    if (ci_mode == "feature") {
      fd_by_day <- split(feat_disp, feat_disp$Day)
      ci_list <- vector("list", length(days))
      
      for (i in seq_along(days)) {
        d <- days[i]
        fdd <- fd_by_day[[as.character(d)]]
        if (is.null(fdd) || nrow(fdd) < 2) {
          ci_list[[i]] <- data.table(Day = d, disp_lo = NA_real_, disp_hi = NA_real_)
          next
        }
        vals <- fdd$feat_disp
        vals <- vals[is.finite(vals)]
        if (length(vals) < 2) {
          ci_list[[i]] <- data.table(Day = d, disp_lo = NA_real_, disp_hi = NA_real_)
          next
        }
        
        boot_vals <- numeric(n_boot)
        nF <- length(vals)
        for (b in seq_len(n_boot)) {
          samp <- sample(vals, size = nF, replace = TRUE)
          boot_vals[b] <- agg_fun(samp)
        }
        
        ci_list[[i]] <- data.table(
          Day = d,
          disp_lo = as.numeric(quantile(boot_vals, q_lo_ci, na.rm = TRUE)),
          disp_hi = as.numeric(quantile(boot_vals, q_hi_ci, na.rm = TRUE))
        )
      }
      
      B_ci <- rbindlist(ci_list)
      setkey(B_df, Day); setkey(B_ci, Day)
      B_df[B_ci, c("disp_lo", "disp_hi") := .(i.disp_lo, i.disp_hi)]
      
    } else if (ci_mode == "sample") {
      dt_by_day <- split(dt_use, dt_use$Day)
      ci_list <- vector("list", length(days))
      
      for (i in seq_along(days)) {
        d <- days[i]
        ddt <- dt_by_day[[as.character(d)]]
        
        sample_list <- split(ddt, ddt$SampleID)
        snames <- names(sample_list)
        nS <- length(snames)
        
        if (nS < 2) {
          ci_list[[i]] <- data.table(Day = d, disp_lo = NA_real_, disp_hi = NA_real_)
          next
        }
        
        boot_vals <- numeric(n_boot)
        for (b in seq_len(n_boot)) {
          idx_b <- sample.int(nS, size = nS, replace = TRUE)
          dboot <- rbindlist(sample_list[idx_b], use.names = TRUE)
          fd_b <- dboot[, .(feat_disp = disp_fun(Abundance)), by = Feature]
          boot_vals[b] <- agg_fun(fd_b$feat_disp)
        }
        
        ci_list[[i]] <- data.table(
          Day = d,
          disp_lo = as.numeric(quantile(boot_vals, q_lo_ci, na.rm = TRUE)),
          disp_hi = as.numeric(quantile(boot_vals, q_hi_ci, na.rm = TRUE))
        )
      }
      B_ci <- rbindlist(ci_list)
      setkey(B_df, Day); setkey(B_ci, Day)
      B_df[B_ci, c("disp_lo", "disp_hi") := .(i.disp_lo, i.disp_hi)]
    }
    
    setorder(B_df, Day)
  }
  
  # ---- Plot ----
  transform_note <- paste0(
    "Transform=", transform_label,
    if (transform_choice == "log2p") paste0(" (pc=", pseudocount, ")") else ""
  )
  
  pA <- ggplot(A_df, aes(x = Day)) +
    geom_ribbon(aes(ymin = band_lo, ymax = band_hi), alpha = 0.25) +
    geom_line(aes(y = mean_val)) +
    geom_point(aes(y = mean_val)) +
    labs(
      title = "A",
      x = "Day",
      y = "System mean (across features)",
      subtitle = paste0(
        band_label,
        if (!is.na(source_label)) paste0(" | Source=", source_label) else "",
        " | ", transform_note
      )
    ) +
    theme_classic()
  
  pB <- ggplot(B_df, aes(x = Day)) +
    {if (ci_on) geom_ribbon(aes(ymin = disp_lo, ymax = disp_hi), alpha = 0.25)} +
    geom_line(aes(y = disp)) +
    geom_point(aes(y = disp)) +
    labs(
      title = "B",
      x = "Day",
      y = paste0(disp_agg, " feature-wise ", disp_metric, " (within day)"),
      subtitle = paste0(
        if (ci_on) paste0(round(ci_level * 100), "% bootstrap CI (mode=", ci_mode, ", n_boot=", n_boot, ")") else "CI off",
        if (!is.na(source_label)) paste0(" | Source=", source_label) else "",
        " | ", transform_note
      )
    ) +
    theme_classic()
  
  fig <- pA / pB
  
  # ---- Write outputs ----
  base <- paste0("FigureX_mean_dispersion_", safe_slug(file_tag), "_", safe_slug(transform_label))
  out_png <- file.path(output_dir, paste0(base, ".png"))
  out_csv <- file.path(output_dir, paste0(base, "_daily_summary.csv"))
  out_pdf <- file.path(output_dir, paste0(base, ".pdf"))
  out_rds_pA <- file.path(output_dir, paste0(base, "_pA.rds"))
  out_rds_pB <- file.path(output_dir, paste0(base, "_pB.rds"))
  
  daily_summary <- merge(A_df, B_df, by = "Day", all = TRUE, suffixes = c(".A", ".B"))
  daily_summary[, Source := source_label]
  daily_summary[, Transform := transform_label]
  daily_summary[, Pseudocount := ifelse(transform_choice == "log2p", pseudocount, NA_real_)]
  
  fwrite(daily_summary, out_csv)
  ggsave(out_png, fig, width = 7.5, height = 7.5, dpi = 300)
  saveRDS(pA, out_rds_pA)
  saveRDS(pB, out_rds_pB)
  
  message("Wrote: ", out_png)
  message("Wrote: ", out_csv)
  message("Wrote: ", out_rds_pA)
  message("Wrote: ", out_rds_pB)
  
  if (write_pdf) {
    ggsave(out_pdf, fig, width = 7.5, height = 7.5)
    message("Wrote: ", out_pdf)
  }
  
  invisible(list(A_df = A_df, B_df = B_df, daily = daily_summary,
                 png = out_png, csv = out_csv, pA_rds = out_rds_pA, pB_rds = out_rds_pB))
}

# --------------------------
# Run analyses
# --------------------------
results <- list()
all_plot_pngs <- character(0)
all_daily_csvs <- character(0)
all_plot_rds <- character(0)

if (split_by_source && has_source) {
  src_vals <- sort(unique(dt0[!is.na(Source), Source]))
  for (s in src_vals) {
    dts <- dt0[Source == s]
    if (nrow(dts) == 0) next
    message("\n--- Analyzing Source: ", s, " (rows=", nrow(dts), ") ---")
    res <- compute_two_panel(
      dts,
      source_label = s,
      output_dir = output_dir,
      file_tag = paste0("Source_", s),
      transform_choice = transform_choice,
      transform_label = transform_label,
      pseudocount = pseudocount
    )
    results[[as.character(s)]] <- res
    all_plot_pngs <- c(all_plot_pngs, res$png)
    all_daily_csvs <- c(all_daily_csvs, res$csv)
    all_plot_rds <- c(all_plot_rds, res$pA_rds, res$pB_rds)
  }
  
  if (include_pooled) {
    message("\n--- Analyzing POOLED across sources (rows=", nrow(dt0), ") ---")
    res <- compute_two_panel(
      dt0,
      source_label = "POOLED",
      output_dir = output_dir,
      file_tag = "POOLED",
      transform_choice = transform_choice,
      transform_label = transform_label,
      pseudocount = pseudocount
    )
    results[["POOLED"]] <- res
    all_plot_pngs <- c(all_plot_pngs, res$png)
    all_daily_csvs <- c(all_daily_csvs, res$csv)
    all_plot_rds <- c(all_plot_rds, res$pA_rds, res$pB_rds)
  }
  
} else {
  source_label <- if (has_source && n_sources == 1) unique(dt0[!is.na(Source), Source])[1] else "POOLED"
  message("\n--- Analyzing dataset (split_by_source=FALSE); Source label: ", source_label, " ---")
  res <- compute_two_panel(
    dt0,
    source_label = source_label,
    output_dir = output_dir,
    file_tag = source_label,
    transform_choice = transform_choice,
    transform_label = transform_label,
    pseudocount = pseudocount
  )
  results[[source_label]] <- res
  all_plot_pngs <- c(all_plot_pngs, res$png)
  all_daily_csvs <- c(all_daily_csvs, res$csv)
  all_plot_rds <- c(all_plot_rds, res$pA_rds, res$pB_rds)
}

# --------------------------
# Manifest helpers (MANDATORY)
# --------------------------
pkg_versions <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    data.frame(package = p, version = v, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

write_file_inventory <- function(output_dir, inventory_path) {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  if (length(files) == 0) {
    inv <- data.table(
      file = character(0),
      rel_path = character(0),
      bytes = numeric(0),
      mtime = character(0)
    )
  } else {
    info <- file.info(files)
    inv <- data.table(
      file = files,
      rel_path = gsub(paste0("^", gsub("\\\\", "/", normalizePath(output_dir, winslash = "/", mustWork = FALSE)), "/?"), "", gsub("\\\\", "/", files)),
      bytes = as.numeric(info$size),
      mtime = as.character(info$mtime)
    )
  }
  fwrite(inv, inventory_path)
  inv
}

# --------------------------
# Runtime tracking end (MANDATORY if enabled)
# --------------------------
end_time <- Sys.time()
runtime_seconds <- NA_real_
if (enable_runtime_tracking) {
  runtime_seconds <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
}

# --------------------------
# Write initial manifest (before report render; updated again after)
# --------------------------
manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
manifest_files_csv <- file.path(output_dir, "Project_Manifest_Files.csv")

deps_df <- pkg_versions(unique(base_pkgs))
deps_list <- lapply(seq_len(nrow(deps_df)), function(i) list(package = deps_df$package[i], version = deps_df$version[i]))

manifest <- list(
  run_id = run_id,
  run_timestamp = as.character(run_timestamp),
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA_character_, script_full)
  ),
  input = list(
    input_path = normalizePath(in_path, winslash = "/", mustWork = FALSE),
    detected_columns = list(
      Day = col_day,
      SampleID = col_sample,
      Feature = col_feature,
      Abundance = col_abundance,
      Source = ifelse(is.null(col_source), NA_character_, col_source)
    )
  ),
  parameters = list(
    analysis_name = analysis_name,
    band_mode = band_mode,
    q_lo_band = q_lo_band,
    q_hi_band = q_hi_band,
    disp_metric = disp_metric,
    disp_agg = disp_agg,
    transform = transform_choice,
    transform_label = transform_label,
    pseudocount = ifelse(transform_choice == "log2p", pseudocount, NA_real_),
    legacy_log_abundance_flag = log_abundance,
    write_pdf = write_pdf,
    split_by_source = split_by_source,
    include_pooled = include_pooled,
    source_regex = source_regex,
    ci_on = ci_on,
    ci_level = ci_level,
    ci_mode = ci_mode,
    n_boot = n_boot,
    seed = seed,
    enable_plotly = enable_plotly,
    enable_runtime_tracking = enable_runtime_tracking,
    runtime_seconds = runtime_seconds
  ),
  dependencies = deps_list,
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
)

jsonlite::write_json(manifest, manifest_json_path, auto_unbox = TRUE, pretty = TRUE)
write_file_inventory(output_dir, manifest_files_csv)

# --------------------------
# QMD generation (MANDATORY)
# --------------------------
get_header_block <- function() {
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
    take_n <- min(length(lines), 220)
    head_lines <- lines[seq_len(take_n)]
    keep <- c()
    started <- FALSE
    for (ln in head_lines) {
      if (!started) {
        if (grepl("^\\s*#!/", ln)) {
          keep <- c(keep, ln)
          next
        }
        if (grepl("^\\s*#|^\\s*$", ln)) {
          keep <- c(keep, ln)
          started <- TRUE
        } else {
          break
        }
      } else {
        if (grepl("^\\s*#|^\\s*$", ln)) {
          keep <- c(keep, ln)
        } else {
          break
        }
      }
    }
    paste(keep, collapse = "\n")
  } else {
    paste0("# Header block unavailable (script_full not detected). Fallback filename: ", known_script_filename)
  }
}

header_text <- get_header_block()

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))

static_png_for_report <- if (length(all_plot_pngs) > 0) all_plot_pngs[1] else NA_character_
log_note_lines <- if (transform_choice %in% c("log2p","log1p")) {
  c(
    "",
    "## Statistical note (log-scale mean)",
    "",
    "When the analysis uses a log transform, the 'System mean' shown in Panel A is the **arithmetic mean on the log scale**.",
    "This corresponds to a **geometric mean on the original scale** (after reversing the transform).",
    "This is standard in omics workflows because it reduces leverage from extreme values and emphasizes fold-change logic."
  )
} else {
  character(0)
}

qmd <- c(
  "---",
  paste0('title: "', script_name, " report\""),
  "format:",
  "  html:",
  "    toc: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "params:",
  paste0("  run_id: \"", run_id, "\""),
  paste0("  run_timestamp: \"", as.character(run_timestamp), "\""),
  paste0("  script_name: \"", script_name, "\""),
  paste0("  script_path: \"", gsub('\"', '\\\\\"', script_path), "\""),
  paste0("  script_full: \"", gsub('\"', '\\\\\"', ifelse(is.na(script_full), "NA", script_full)), "\""),
  paste0("  input_path: \"", gsub('\"', '\\\\\"', normalizePath(in_path, winslash = "/", mustWork = FALSE)), "\""),
  paste0("  output_dir: \"", gsub('\"', '\\\\\"', normalizePath(output_dir, winslash = "/", mustWork = FALSE)), "\""),
  paste0("  runtime_seconds: \"", ifelse(is.na(runtime_seconds), "NA", as.character(runtime_seconds)), "\""),
  paste0("  enable_plotly: ", ifelse(enable_plotly, "true", "false")),
  paste0("  transform: \"", transform_choice, "\""),
  paste0("  transform_label: \"", transform_label, "\""),
  paste0("  pseudocount: \"", ifelse(transform_choice == "log2p", as.character(pseudocount), "NA"), "\""),
  "---",
  "",
  "## Overview",
  "",
  "This report summarizes a two-panel time-series visualization built from tidy long data.",
  "",
  "### Transform applied to Abundance",
  "",
  paste0("- **transform:** `", transform_choice, "`"),
  paste0("- **transform_label:** `", transform_label, "`"),
  paste0("- **pseudocount (log2p only):** `", ifelse(transform_choice == "log2p", as.character(pseudocount), "NA"), "`"),
  "",
  "Panel A (Mean + replicate band):",
  "- For each Day and SampleID, compute the within-sample mean across features.",
  "- Then summarize across SampleID within Day: the mean (line/points) and a replicate dispersion band (IQR or user-specified quantiles).",
  "",
  "Panel B (System-level dispersion across features):",
  "- For each Day and Feature, compute a within-day dispersion metric across SampleID (SD, CV, IQR, or MAD).",
  "- Aggregate feature dispersions within Day using median or mean.",
  "- Optionally compute a bootstrap confidence interval either by resampling features or resampling SampleIDs within day.",
  "",
  "## Metadata",
  "",
  paste0("- **Run ID:** `", run_id, "`"),
  paste0("- **Run timestamp:** `", as.character(run_timestamp), "`"),
  paste0("- **Input:** `", normalizePath(in_path, winslash = "/", mustWork = FALSE), "`"),
  paste0("- **Output directory:** `", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "`"),
  "",
  "### Script identity",
  "",
  paste0("- **script_name:** `", script_name, "`"),
  paste0("- **script_path:** `", script_path, "`"),
  paste0("- **script_full:** `", ifelse(is.na(script_full), "NA", script_full), "`"),
  "",
  "## Script header",
  "",
  "```",
  header_text,
  "```",
  "",
  "## Dependencies (package versions)",
  "",
  "```{r}",
  "manifest_path <- file.path(params$output_dir, 'Project_Manifest.json')",
  "deps_out <- NULL",
  "if (file.exists(manifest_path)) {",
  "  man <- tryCatch(jsonlite::fromJSON(manifest_path), error = function(e) NULL)",
  "  deps <- tryCatch(man$dependencies, error = function(e) NULL)",
  "",
  "  # Normalize deps into a data.frame(package, version) robustly",
  "  if (is.data.frame(deps)) {",
  "    # Common case: deps already a data.frame with columns package/version",
  "    if (all(c('package','version') %in% names(deps))) {",
  "      deps_out <- deps[, c('package','version')]",
  "    } else if (ncol(deps) >= 2) {",
  "      deps_out <- data.frame(package = deps[[1]], version = deps[[2]], stringsAsFactors = FALSE)",
  "    }",
  "  } else if (is.list(deps)) {",
  "    # List-of-lists case",
  "    deps_out <- tryCatch({",
  "      do.call(rbind, lapply(deps, function(x) {",
  "        if (is.list(x) && !is.null(x$package)) {",
  "          data.frame(package = x$package, version = x$version, stringsAsFactors = FALSE)",
  "        } else if (is.character(x) && length(x) >= 2) {",
  "          data.frame(package = x[[1]], version = x[[2]], stringsAsFactors = FALSE)",
  "        } else {",
  "          NULL",
  "        }",
  "      }))",
  "    }, error = function(e) NULL)",
  "  } else if (is.atomic(deps) && length(deps) > 0) {",
  "    # Fallback: atomic vector; cannot reliably parse versions",
  "    deps_out <- data.frame(package = as.character(deps), version = NA_character_, stringsAsFactors = FALSE)",
  "  }",
  "}",
  "",
  "if (!is.null(deps_out) && nrow(deps_out) > 0) {",
  "  knitr::kable(deps_out)",
  "} else {",
  "  cat('Dependency list not available.')",
  "}",
  "```",
  "",
  "## Generated files",
  "",
  "```{r}",
  "inv <- data.table::fread(file.path(params$output_dir, 'Project_Manifest_Files.csv'))",
  "knitr::kable(inv)",
  "```",
  "",
  "## Key plots",
  "",
  "### Static figure (always-on)",
  "",
  "```{r}",
  "png_path <- NA_character_",
  if (!is.na(static_png_for_report)) paste0("png_path <- '", gsub("\\\\", "/", static_png_for_report), "'") else "png_path <- NA_character_",
  "if (!is.na(png_path) && nzchar(png_path) && file.exists(png_path)) {",
  "  knitr::include_graphics(png_path)",
  "} else {",
  "  cat('No PNG figure found to display.')",
  "}",
  "```",
  "",
  "### Interactive plots (plotly; if enabled)",
  "",
  "```{r}",
  "if (isTRUE(params$enable_plotly)) {",
  "  ok <- requireNamespace('plotly', quietly = TRUE)",
  "  if (!ok) {",
  "    cat('plotly is not available in this R session; falling back to static plot above.')",
  "  } else {",
  "    rds_files <- list.files(params$output_dir, pattern = '\\\\_pA\\\\.rds$|\\\\_pB\\\\.rds$', full.names = TRUE)",
  "    if (length(rds_files) == 0) {",
  "      cat('No plot RDS files found for interactive rendering.')",
  "    } else {",
  "      p <- tryCatch(readRDS(rds_files[1]), error = function(e) NULL)",
  "      if (is.null(p)) {",
  "        cat('Could not read plot RDS; falling back to static plot above.')",
  "      } else {",
  "        tryCatch({",
  "          plotly::ggplotly(p)",
  "        }, error = function(e) {",
  "          cat('plotly rendering failed; falling back to static plot above.')",
  "        })",
  "      }",
  "    }",
  "  }",
  "} else {",
  "  cat('Interactive plots disabled (enable_plotly=FALSE).')",
  "}",
  "```",
  "",
  "## Interpretation guide",
  "",
  "- **Panel A** reflects how the system-wide mean signal changes with Day, with the band showing between-sample dispersion of per-sample means.",
  "- **Panel B** reflects how heterogeneous features are within each day (after computing within-day, within-feature dispersion across replicates), summarized to a single system-level dispersion value per day.",
  "- A rise in Panel B can occur without a rise in Panel A, consistent with greater within-day variability or diminished constraint even if mean levels are stable.",
  log_note_lines,
  "",
  "## Runtime",
  "",
  "```{r}",
  "if (!is.na(params$runtime_seconds) && params$runtime_seconds != 'NA') {",
  "  cat('Total runtime (seconds): ', params$runtime_seconds)",
  "} else {",
  "  cat('Runtime tracking disabled or unavailable.')",
  "}",
  "```",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```",
  ""
)

writeLines(qmd, qmd_path)
message("\nQMD written at: ", qmd_path)

# --------------------------
# Render report to HTML (MANDATORY final step)
# --------------------------
render_ok <- FALSE
render_err <- NULL

tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_format = "html_document",
    output_file = basename(html_path),
    output_dir = output_dir,
    quiet = TRUE,
    params = list(
      run_id = run_id,
      run_timestamp = as.character(run_timestamp),
      script_name = script_name,
      script_path = script_path,
      script_full = ifelse(is.na(script_full), "NA", script_full),
      input_path = normalizePath(in_path, winslash = "/", mustWork = FALSE),
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      runtime_seconds = ifelse(is.na(runtime_seconds), "NA", as.character(runtime_seconds)),
      enable_plotly = enable_plotly,
      transform = transform_choice,
      transform_label = transform_label,
      pseudocount = ifelse(transform_choice == "log2p", as.character(pseudocount), "NA")
    )
  )
  render_ok <- TRUE
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- conditionMessage(e)
})

if (render_ok && file.exists(html_path)) {
  message("HTML report created: ", html_path)
} else {
  message("ERROR: HTML report render failed.")
  if (!is.null(render_err)) message("Reason: ", render_err)
  message("QMD available at: ", qmd_path)
}

# --------------------------
# Refresh file inventory after render (MANDATORY)
# --------------------------
inv <- write_file_inventory(output_dir, manifest_files_csv)

manifest$generated_files <- as.list(inv$rel_path)
jsonlite::write_json(manifest, manifest_json_path, auto_unbox = TRUE, pretty = TRUE)

message("\nDone.")
