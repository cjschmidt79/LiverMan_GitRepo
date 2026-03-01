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
#   --log_abundance=0            # 1 to log1p-transform Abundance
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
# Retrofits (MANDATORY PACKET):
#   - All outputs forced under: ./outputs/<run_id>/
#   - Project manifest JSON + file inventory CSV
#   - Script identity capture (RStudio + Rscript + fallback)
#   - Quarto QMD report written + rendered to HTML as FINAL step
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

# Minimum packages for the retrofit layer (mandatory)
quiet_install_if_missing(c("jsonlite", "knitr", "rmarkdown"))

# Analysis packages (existing script)
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
# Inputs + parameters (existing behavior preserved; outputs policy enforced)
# --------------------------
in_path <- arg_val("in", NULL)

# Output routing policy (MANDATORY): always under ./outputs/<run_id>/
outputs_root <- file.path(getwd(), "outputs")
dir.create(outputs_root, showWarnings = FALSE, recursive = TRUE)

# The script previously accepted --out_dir; now ignored by design to satisfy policy.
out_dir_ignored <- arg_val("out_dir", NULL)
if (!is.null(out_dir_ignored) && nzchar(out_dir_ignored)) {
  message("[policy] NOTE: --out_dir was provided but is ignored. All outputs are written under: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE))
}

band_mode <- tolower(arg_val("band", "iqr"))                 # iqr | quantiles
q_lo_band <- as.numeric(arg_val("q_lo", "0.25"))
q_hi_band <- as.numeric(arg_val("q_hi", "0.75"))

disp_metric <- tolower(arg_val("disp_metric", "sd"))         # sd | cv | iqr | mad
disp_agg <- tolower(arg_val("disp_agg", "median"))           # median | mean

log_abundance <- arg_flag("log_abundance", FALSE)
write_pdf <- arg_flag("pdf", FALSE)

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

# Optional run naming parameters
analysis_name <- arg_val("analysis_name", "TwoPanelFigure")  # used in run_id; user may override
run_source_tag <- arg_val("run_source", "")                 # optional: user override for run_id's Source component

# Validate CI inputs
if (!(ci_level > 0 && ci_level < 1)) stop2("--ci_level must be between 0 and 1.")
if (!(ci_mode %in% c("sample", "feature"))) stop2("--ci_mode must be 'sample' or 'feature'.")

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
# Input path (existing behavior preserved)
# --------------------------
if (is.null(in_path) || nchar(in_path) == 0) {
  if (interactive() && exists("file.choose")) {
    in_path <- file.choose()
  } else {
    stop2("No input provided. Run:\n  Rscript make_two_panel_figure.R --in=\"path/to/file.csv\"")
  }
}

message("Input:  ", in_path)
message("Panel B CI: ", if (ci_on) paste0("ON (", round(ci_level * 100), "%, mode=", ci_mode, ", n_boot=", n_boot, ")") else "OFF")

# --------------------------
# Read data (unchanged analysis logic)
# --------------------------
dt0 <- fread(in_path, showProgress = FALSE)
if (nrow(dt0) == 0) stop2("Input file has 0 rows.")
setDT(dt0)

# --------------------------
# Column detection helpers (unchanged)
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

# Detect columns unless explicitly provided
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
# Standardize columns into a working dt (unchanged)
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

if (log_abundance) {
  dt0[, Abundance := log1p(Abundance)]
  message("Applied log1p transform to Abundance.")
}

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
  split_by_source <- tolower(trimws(split_by_source_raw)) %in% c("1", "true", "t", "yes", "y")
}

message("Source handling: ",
        if (!has_source) "no Source column detected"
        else paste0(n_sources, " source level(s) detected; split_by_source=", split_by_source,
                    if (include_pooled) " (also pooled)" else ""))

# --------------------------
# Run folder creation (MANDATORY: outputs/<run_id>/)
#   run_id should be <analysis_name>_<Source>_<YYYYMMDD_HHMMSS>
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
  deps <- lapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) {
      list(package = p, version = as.character(utils::packageVersion(p)))
    } else {
      list(package = p, version = NA_character_)
    }
  })
  deps
}

manifest_path_json <- file.path(output_dir, "Project_Manifest.json")
manifest_path_files <- file.path(output_dir, "Project_Manifest_Files.csv")

write_manifest_json <- function(manifest) {
  jsonlite::write_json(manifest, manifest_path_json, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

refresh_file_inventory <- function() {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (length(files) == 0) {
    inv <- data.table(
      file = character(),
      relative_path = character(),
      bytes = numeric(),
      modified_time = character()
    )
  } else {
    rel <- sub(paste0("^", gsub("([\\^\\$\\.|\\+\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", normalizePath(output_dir, winslash = "/", mustWork = FALSE)), "/?"), "", normalizePath(files, winslash = "/", mustWork = FALSE))
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

# Initial manifest write (before outputs exist; will be refreshed later)
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
    band = band_mode,
    q_lo = q_lo_band,
    q_hi = q_hi_band,
    disp_metric = disp_metric,
    disp_agg = disp_agg,
    log_abundance = log_abundance,
    pdf = write_pdf,
    split_by_source = split_by_source,
    include_pooled = include_pooled,
    source_regex = source_regex,
    ci = ci_on,
    ci_level = ci_level,
    ci_mode = ci_mode,
    n_boot = n_boot,
    seed = seed
  ),
  dependencies = get_deps(required_pkgs),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  ),
  generated_files = list() # populated via CSV inventory; refreshed at end
)
write_manifest_json(manifest)
refresh_file_inventory()

# --------------------------
# Core functions for one dataset (one Source subset)
# --------------------------
disp_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  if (disp_metric == "sd")  return(sd(x))
  if (disp_metric == "cv")  { m <- mean(x); if (m == 0) return(NA_real_); return(sd(x) / abs(m)) }
  if (disp_metric == "iqr") return(IQR(x))
  if (disp_metric == "mad") return(mad(x, constant = 1))
  stop2("Unknown --disp_metric: ", disp_metric, " (sd, cv, iqr, mad)")
}
agg_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (disp_agg == "median") return(median(x))
  if (disp_agg == "mean")   return(mean(x))
  stop2("Unknown --disp_agg: ", disp_agg, " (median or mean)")
}

compute_two_panel <- function(dt, source_label, out_dir, file_tag) {
  # dt must have: Day, SampleID, Feature, Abundance
  # source_label used only for labeling outputs
  # file_tag used in filenames
  
  # ---- Panel A ----
  rep_means <- dt[, .(rep_mean = mean(Abundance, na.rm = TRUE)), by = .(Day, SampleID)]
  
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
  gene_disp <- dt[, .(feat_disp = disp_fun(Abundance)), by = .(Day, Feature)]
  gene_disp <- gene_disp[is.finite(feat_disp)]
  if (nrow(gene_disp) == 0) stop2("No feature-level dispersion values computed for Source=", source_label)
  
  B_df <- gene_disp[, .(
    disp   = agg_fun(feat_disp),
    n_feat = sum(is.finite(feat_disp))
  ), by = Day][order(Day)]
  
  B_df[, `:=`(disp_lo = NA_real_, disp_hi = NA_real_, ci_mode = if (ci_on) ci_mode else NA_character_)]
  
  # ---- Panel B CI ----
  if (ci_on) {
    alpha <- (1 - ci_level) / 2
    q_lo_ci <- alpha
    q_hi_ci <- 1 - alpha
    days <- sort(unique(dt$Day))
    
    if (ci_mode == "feature") {
      gd_by_day <- split(gene_disp, gene_disp$Day)
      ci_list <- vector("list", length(days))
      
      for (i in seq_along(days)) {
        d <- days[i]
        gdd <- gd_by_day[[as.character(d)]]
        if (is.null(gdd) || nrow(gdd) < 2) {
          ci_list[[i]] <- data.table(Day = d, disp_lo = NA_real_, disp_hi = NA_real_)
          next
        }
        vals <- gdd$feat_disp
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
      dt_by_day <- split(dt, dt$Day)
      ci_list <- vector("list", length(days))
      
      for (i in seq_along(days)) {
        d <- days[i]
        ddt <- dt_by_day[[as.character(d)]]
        sids <- unique(ddt$SampleID)
        if (length(sids) < 2) {
          ci_list[[i]] <- data.table(Day = d, disp_lo = NA_real_, disp_hi = NA_real_)
          next
        }
        
        boot_vals <- numeric(n_boot)
        for (b in seq_len(n_boot)) {
          sids_b <- sample(sids, size = length(sids), replace = TRUE)
          dboot <- rbindlist(lapply(sids_b, function(id) ddt[SampleID == id]), use.names = TRUE)
          gd_b <- dboot[, .(feat_disp = disp_fun(Abundance)), by = Feature]
          boot_vals[b] <- agg_fun(gd_b$feat_disp)
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
  pA <- ggplot(A_df, aes(x = Day)) +
    geom_ribbon(aes(ymin = band_lo, ymax = band_hi), alpha = 0.25) +
    geom_line(aes(y = mean_val)) +
    geom_point(aes(y = mean_val)) +
    labs(
      title = "A",
      x = "Day",
      y = "System mean (across features)",
      subtitle = paste0(band_label, if (!is.na(source_label)) paste0(" | Source=", source_label) else "")
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
      subtitle = if (ci_on) paste0(round(ci_level * 100), "% bootstrap CI (mode=", ci_mode, ", n_boot=", n_boot, ")",
                                   if (!is.na(source_label)) paste0(" | Source=", source_label) else "") else NULL
    ) +
    theme_classic()
  
  fig <- pA / pB
  
  # ---- Write outputs ----
  base <- paste0("FigureX_mean_dispersion_", safe_slug(file_tag), "_", ts_stamp())
  out_png <- file.path(out_dir, paste0(base, ".png"))
  out_csv <- file.path(out_dir, paste0(base, "_daily_summary.csv"))
  out_pdf <- file.path(out_dir, paste0(base, ".pdf"))
  
  daily_summary <- merge(A_df, B_df, by = "Day", all = TRUE, suffixes = c(".A", ".B"))
  daily_summary[, Source := source_label]
  
  fwrite(daily_summary, out_csv)
  ggsave(out_png, fig, width = 7.5, height = 7.5, dpi = 300)
  
  message("Wrote: ", out_png)
  message("Wrote: ", out_csv)
  
  if (write_pdf) {
    ggsave(out_pdf, fig, width = 7.5, height = 7.5)
    message("Wrote: ", out_pdf)
  }
  
  invisible(list(A_df = A_df, B_df = B_df, daily = daily_summary, png = out_png, csv = out_csv, pdf = if (write_pdf) out_pdf else NA_character_))
}

# --------------------------
# Run analyses (unchanged; only output_dir routing)
# --------------------------
results <- list()

if (split_by_source && has_source) {
  src_vals <- sort(unique(dt0[!is.na(Source), Source]))
  for (s in src_vals) {
    dts <- dt0[Source == s]
    if (nrow(dts) == 0) next
    message("\n--- Analyzing Source: ", s, " (rows=", nrow(dts), ") ---")
    results[[as.character(s)]] <- compute_two_panel(
      dts,
      source_label = s,
      out_dir = output_dir,
      file_tag = paste0("Source_", s)
    )
  }
  
  if (include_pooled) {
    message("\n--- Analyzing POOLED across sources (rows=", nrow(dt0), ") ---")
    results[["POOLED"]] <- compute_two_panel(
      dt0,
      source_label = "POOLED",
      out_dir = output_dir,
      file_tag = "POOLED"
    )
  }
  
} else {
  # Single run (pooled or single-source data)
  source_label <- if (has_source && n_sources == 1) unique(dt0[!is.na(Source), Source])[1] else "POOLED"
  message("\n--- Analyzing dataset (split_by_source=FALSE); Source label: ", source_label, " ---")
  results[[source_label]] <- compute_two_panel(
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
# Embed the script's header text (initial comment block) as a fenced code block.
SCRIPT_HEADER_TEXT <- paste0(
  "#!/usr/bin/env Rscript\n",
  "# ============================================================\n",
  "# Two-panel Figure (Mean + Dispersion) from Tidy Long Data\n",
  "#   - NO hard-coded column names\n",
  "#   - Auto-detect Day / SampleID / Feature / Abundance / Source (optional)\n",
  "#   - If Source exists and split_by_source=1: analyze each Source independently\n",
  "#\n",
  "# Panel A: system mean with replicate dispersion band (across SampleID within Day)\n",
  "# Panel B: system-level dispersion across features (within Day) + optional 95% CI\n",
  "#\n",
  "# Expected data format (generic tidy long):\n",
  "#   Day | SampleID | Feature (Gene/Metabolite/etc) | Abundance\n",
  "# Optional:\n",
  "#   Source\n",
  "#\n",
  "# Usage (CLI):\n",
  "#   Rscript make_two_panel_figure.R --in=\"path/to/file.csv\"\n",
  "#\n",
  "# Optional flags:\n",
  "#   --band=\"iqr\"                 # iqr or quantiles\n",
  "#   --q_lo=0.25 --q_hi=0.75      # if band=quantiles\n",
  "#   --disp_metric=\"sd\"           # sd, cv, iqr, mad  (per feature within day)\n",
  "#   --disp_agg=\"median\"          # median or mean aggregation across features\n",
  "#   --log_abundance=0            # 1 to log1p-transform Abundance\n",
  "#   --pdf=0                      # 1 to also write PDF\n",
  "#\n",
  "# Source handling:\n",
  "#   --split_by_source=1          # if Source exists: analyze each Source separately (default auto)\n",
  "#   --include_pooled=0           # also run pooled analysis across all sources (default 0)\n",
  "#   --source_regex=\"\"            # optional regex filter on Source values (case-insensitive)\n",
  "#\n",
  "# Confidence interval options (Panel B):\n",
  "#   --ci=1                       # 1 to compute CI for Panel B (default: 1)\n",
  "#   --ci_level=0.95              # CI level (default 0.95)\n",
  "#   --ci_mode=\"sample\"           # \"sample\": resample SampleID within Day (within Source)\n",
  "#                                # \"feature\": resample feature dispersions within Day (within Source)\n",
  "#   --n_boot=500                 # number of bootstrap replicates (if omitted, script will ask interactively)\n",
  "#   --seed=1                     # random seed\n",
  "#\n",
  "# Column overrides (only if auto-detection fails):\n",
  "#   --col_day= --col_sample= --col_feature= --col_abundance= --col_source=\n",
  "# ============================================================\n"
)

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_name <- paste0(script_name, "_Report.html")
html_path <- file.path(output_dir, html_name)

# Choose a representative always-on plot: first PNG in output_dir (if present)
png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
plot_to_include <- if (length(png_files) > 0) normalizePath(png_files[1], winslash = "/", mustWork = FALSE) else NA_character_

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
  "## Script header",
  "",
  "```",
  SCRIPT_HEADER_TEXT,
  "```",
  "",
  "## Dependencies",
  "",
  "```{r}",
  "library(jsonlite)",
  "man <- jsonlite::fromJSON(file.path(\".\", \"Project_Manifest.json\"), simplifyVector = TRUE)",
  "deps <- man$dependencies",
  "deps_df <- if (length(deps) == 0) data.frame() else do.call(rbind, lapply(deps, as.data.frame))",
  "knitr::kable(deps_df, align = 'l')",
  "```",
  "",
  "## Generated files",
  "",
  "```{r}",
  "inv <- data.table::fread(file.path(\".\", \"Project_Manifest_Files.csv\"))",
  "knitr::kable(inv, align = 'l')",
  "```",
  "",
  "## Analytical logic and formulas",
  "",
  "This analysis produces a two-panel summary across post-hoc timepoints (Day) from tidy long data.",
  "",
  "- **Panel A (System mean):** For each Day and SampleID, compute the replicate-level mean across features:",
  "  ",
  "  \\[ \\bar{x}_{d,s} = \\frac{1}{F} \\sum_{f=1}^{F} x_{d,s,f} \\]",
  "  ",
  "  Then compute the Day-level mean across SampleIDs and a replicate dispersion band (IQR or user-specified quantiles) over \\(\\bar{x}_{d,s}\\).",
  "",
  "- **Panel B (System dispersion across features):** For each Day and Feature, compute a within-day dispersion metric over SampleIDs:",
  "  ",
  "  - SD: \\(\\mathrm{sd}(x_{d,\\cdot,f})\\)",
  "  - CV: \\(\\mathrm{sd}(x_{d,\\cdot,f}) / |\\mathrm{mean}(x_{d,\\cdot,f})|\\)",
  "  - IQR: \\(\\mathrm{IQR}(x_{d,\\cdot,f})\\)",
  "  - MAD: \\(\\mathrm{mad}(x_{d,\\cdot,f};\\, c=1)\\)",
  "  ",
  "  Aggregate feature-wise dispersions within each Day using the median (default) or mean:",
  "  ",
  "  \\[ D_d = \\mathrm{median}_f\\{\\delta_{d,f}\\} \\quad \\text{or} \\quad D_d = \\mathrm{mean}_f\\{\\delta_{d,f}\\} \\]",
  "",
  "- **Bootstrap CI (optional):**",
  "  - `ci_mode = sample`: resample SampleIDs within Day and recompute \\(D_d\\).",
  "  - `ci_mode = feature`: resample feature dispersions \\(\\delta_{d,f}\\) within Day and recompute \\(D_d\\).",
  "",
  "## Plots",
  "",
  "```{r}",
  "plot_path <- man$generated_files$preferred_plot",
  "if (is.null(plot_path) || is.na(plot_path) || !file.exists(plot_path)) {",
  "  # fallback: choose first png in the directory",
  "  pngs <- list.files(\".\", pattern = \"\\\\.png$\", full.names = TRUE)",
  "  if (length(pngs) > 0) plot_path <- pngs[1] else plot_path <- NA_character_",
  "}",
  "if (!is.na(plot_path) && file.exists(plot_path)) {",
  "  knitr::include_graphics(plot_path)",
  "} else {",
  "  cat(\"No PNG plot was found to include in the report.\")",
  "}",
  "```",
  "",
  "## Interpretation of results",
  "",
  "Panel A summarizes how the system-wide mean abundance changes across Day, with the band capturing between-replicate dispersion of replicate means within each Day. ",
  "Panel B summarizes how within-day variability (computed feature-wise across replicates) behaves across Day, optionally with a bootstrap confidence interval to indicate uncertainty under resampling assumptions.",
  "",
  "When `split_by_source` is enabled and multiple Source levels exist, the above logic is applied independently per Source (and optionally to the pooled dataset), producing separate output figures and daily summary tables.",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

# Persist a preferred plot path into manifest for report convenience
manifest$generated_files$preferred_plot <- if (!is.na(plot_to_include)) plot_to_include else NA_character_
write_manifest_json(manifest)

writeLines(qmd_lines, qmd_path)

render_ok <- FALSE
render_err <- NULL
tryCatch({
  # rmarkdown can render .qmd; it will use Quarto if available; otherwise it will still attempt via rmarkdown machinery.
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

# Mandatory: after rendering (success or failure), refresh inventory
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
