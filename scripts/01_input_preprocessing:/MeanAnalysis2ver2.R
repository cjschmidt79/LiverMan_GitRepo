#!/usr/bin/env Rscript
# ======================================================================
# Script for per-feature mean and variability analysis (log2(x + 1) transformed abundance)
#
#
#Mean Expression / Abundance by Day (TIDY-ONLY)
#   - Computes per-Source × Feature × Day mean + SD + SE + 95% CI
#   - Computes per-Source “average feature mean” by Day (with 95% CI)
#   - OPTION A (requested): Global ALL-sources plot uses BIOLOGICAL REPLICATES:
#       * For each Source×Day×SampleID, compute "system mean" = mean across features.
#       * Across SampleIDs on each Day (pooling Sources), compute:
#            - grey band = IQR or P10–P90 of system means (inter-individual dispersion)
#       # black line = mean or median of system means
#            - CI bars = t-based 95% CI of the MEAN of system means
#       Grey band is dispersion, not uncertainty.
#   - Writes CSV tables + saves plots (base R)
## Usage (RStudio or command line):
#   Rscript MeanAnalysis2.R
#   # or source("MeanAnalysis2.R") in RStudio
# Inputs:
#   - User-selected CSV with EXACT 5 columns: Source, Day, SampleID, Metabolite, Abundance

# TIDY INPUT REQUIRED (5 columns EXACT):
#   Source, Day, SampleID, Metabolite, Abundance
# ======================================================================

# ======================================================================
# Retrofit layer: standardized outputs + manifest + self-report (QMD->HTML)
#   - All outputs forced under: outputs/<run_id>/
#   - Captures script identity (RStudio + Rscript + source() fallback)
#   - Writes:
#       * Project_Manifest.json
#       * Project_Manifest_Files.csv
#       * <script_name>_Report.qmd
#       * <script_name>_Report.html (rendered final step)
#   - Optional:
#       * plotly interactive plots in report
#       * runtime tracking
# ======================================================================

# -----------------------------
# Internal feature flags (mandatory)
# -----------------------------
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

# -----------------------------
# Runtime tracking (mandatory if enabled)
# -----------------------------
start_time <- Sys.time()

# -----------------------------
# Dependency bootstrap (install if missing)
# -----------------------------
quiet_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  requireNamespace(pkg, quietly = TRUE)
}

# Minimum retrofit packages
quiet_require("jsonlite")
quiet_require("rmarkdown")
quiet_require("knitr")

# Optional (plotly layer)
if (isTRUE(enable_plotly)) {
  quiet_require("ggplot2")
  quiet_require("plotly")
} else {
  # Still allow ggplot2 static key plot if present; install only if missing later when needed
  quiet_require("ggplot2")
}

# -----------------------------
# Script identity capture (mandatory)
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

# Fallback filename (fill-in required by packet; inferred here)
known_script_filename <- "MeanByDay_TidyOnly_Retrofit.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  # Path cannot be detected in this execution mode; still record a valid script_name.
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  cat("[NOTE] Script path detection failed; using fallback known_script_filename.\n")
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

# -----------------------------
# User settings
# -----------------------------
set_seed <- 123

log_transform <- TRUE
log_pseudocount <- 1

make_plots <- TRUE
top_n_features_per_source <- 24
min_n_per_day_for_stats <- 2

plot_width <- 1600
plot_height <- 1000
plot_res <- 150

sampleid_policy <- "fill"   # "error" | "drop" | "fill"

# Band options for ALL Sources plot (replicate-based band)
add_global_band <- TRUE
global_band_type <- "IQR"   # "IQR" or "P10P90"
global_center <- "mean"     # "mean" or "median"
band_fill_col <- "gray85"

do_levene_perm <- TRUE
n_perm_levene  <- 1000  # Set to 10000 for publication

do_system_cv_perm <- TRUE
n_perm_system_cv  <- 5000
focus_peak_day    <- NA   # e.g., 8 (set NA to skip day-specific test)
focus_trough_day  <- NA   # e.g., 14

use_parallel <- TRUE
n_cores <- max(1, parallel::detectCores() - 2)

# Analysis name (for run_id; inferred)
analysis_name <- "MeanByDay"

# -----------------------------
# Helpers
# -----------------------------
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

mean_ci95 <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) return(c(mean = NA_real_, sd = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n = 0))
  mu <- mean(x)
  if (n < 2) return(c(mean = mu, sd = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n = n))
  s <- sd(x)
  se <- s / sqrt(n)
  tcrit <- qt(0.975, df = n - 1)
  c(mean = mu, sd = s, se = se, ci_low = mu - tcrit * se, ci_high = mu + tcrit * se, n = n)
}

plot_feature <- function(dd, main_title, ylab_txt, out_file, plot_width, plot_height, plot_res) {
  dd <- dd[order(dd$Day), , drop = FALSE]
  okp <- is.finite(dd$Day) & is.finite(dd$Mean)
  if (sum(okp) < 2) return(invisible(FALSE))
  dd <- dd[okp, , drop = FALSE]
  
  png(out_file, width = plot_width, height = plot_height, res = plot_res)
  
  ylim_vals <- dd$Mean
  okci <- is.finite(dd$CI_Low) & is.finite(dd$CI_High)
  if (any(okci)) ylim_vals <- c(ylim_vals, dd$CI_Low[okci], dd$CI_High[okci])
  ylim_vals <- ylim_vals[is.finite(ylim_vals)]
  if (length(ylim_vals) < 2) { dev.off(); return(invisible(FALSE)) }
  
  plot(dd$Day, dd$Mean, type = "n",
       main = main_title, xlab = "Day", ylab = ylab_txt,
       ylim = range(ylim_vals))
  
  if (any(okci)) segments(dd$Day[okci], dd$CI_Low[okci], dd$Day[okci], dd$CI_High[okci])
  lines(dd$Day, dd$Mean, lwd = 2)
  points(dd$Day, dd$Mean, pch = 16)
  
  dev.off()
  invisible(TRUE)
}

get_header_block <- function() {
  # Extract the initial comment block from the script file if available
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
    take_n <- min(length(lines), 200)
    head_lines <- lines[seq_len(take_n)]
    keep <- c()
    started <- FALSE
    for (ln in head_lines) {
      if (!started) {
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
    return(paste(keep, collapse = "\n"))
  }
  "# Header block unavailable (script_full not detected)."
}

pkg_versions <- function(pkgs) {
  out <- data.frame(package = character(0), version = character(0), stringsAsFactors = FALSE)
  for (p in pkgs) {
    v <- NA_character_
    if (requireNamespace(p, quietly = TRUE)) {
      v <- as.character(utils::packageVersion(p))
    }
    out <- rbind(out, data.frame(package = p, version = v, stringsAsFactors = FALSE))
  }
  out
}

write_file_inventory <- function(output_dir, csv_path) {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  if (length(files) == 0) {
    inv <- data.frame(
      rel_path = character(0),
      size_bytes = numeric(0),
      modified_time = character(0),
      stringsAsFactors = FALSE
    )
    utils::write.csv(inv, csv_path, row.names = FALSE)
    return(inv)
  }
  rel <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", normalizePath(output_dir, winslash = "/")), "/?"), "", normalizePath(files, winslash = "/"))
  info <- file.info(files)
  inv <- data.frame(
    rel_path = rel,
    size_bytes = as.numeric(info$size),
    modified_time = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
  inv <- inv[order(inv$rel_path), , drop = FALSE]
  utils::write.csv(inv, csv_path, row.names = FALSE)
  inv
}

# -----------------------------
# I/O (input selection retained; outputs forced to outputs/<run_id>/)
# -----------------------------
# I/O (input selection via CLI or interactive fallback)
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
input_arg <- grep("^--input_csv=", args, value = TRUE)
if (length(input_arg) > 0) {
  infile <- sub("^--input_csv=", "", input_arg)
  if (!file.exists(infile)) {
    stop("Specified input file does not exist: ", infile)
  }
  infile <- normalizePath(infile, winslash = "/", mustWork = FALSE)
  cat("[INFO] Using input CSV from command line:\n  ", infile, "\n")
} else {
  cat("\nSelect input CSV (TIDY-ONLY; EXACT 5 columns):\n")
  cat("  Source, Day, SampleID, Metabolite, Abundance\n")
  infile <- file.choose()
  infile <- normalizePath(infile, winslash = "/", mustWork = FALSE)
}

# -----------------------------
# Read data (schema enforced)
# -----------------------------
df_raw <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("Source","Day","SampleID","Metabolite","Abundance")
missing_cols <- setdiff(required_cols, names(df_raw))
extra_cols   <- setdiff(names(df_raw), required_cols)

if (length(missing_cols) > 0) {
  stop("Input must be tidy with EXACT columns:\n  ", paste(required_cols, collapse = ", "),
       "\nMissing:\n  ", paste(missing_cols, collapse = ", "))
}
if (length(extra_cols) > 0) {
  stop("Input must contain ONLY these columns:\n  ", paste(required_cols, collapse = ", "),
       "\nExtra columns detected:\n  ", paste(extra_cols, collapse = ", "))
}

df_raw$Source     <- trimws(as.character(df_raw$Source))
df_raw$SampleID   <- trimws(as.character(df_raw$SampleID))
df_raw$Metabolite <- trimws(as.character(df_raw$Metabolite))

df_raw$SampleID[df_raw$SampleID %in% c("NA","N/A","na","n/a","NULL","null",".","-")] <- NA_character_
df_raw$SampleID <- trimws(df_raw$SampleID)

bad_id <- which(is.na(df_raw$SampleID) | df_raw$SampleID == "")
if (length(bad_id) > 0) {
  diag <- df_raw[bad_id, c("Source","Day","Metabolite","Abundance","SampleID")]
  diag <- diag[seq_len(min(10, nrow(diag))), , drop = FALSE]
  msg <- paste0(
    "SampleID is empty/NA after trimws() in ", length(bad_id), " row(s).\n",
    "Example offending rows (up to 10):\n",
    paste(capture.output(print(diag, row.names = TRUE)), collapse = "\n")
  )
  
  if (identical(sampleid_policy, "error")) {
    stop(msg)
  } else if (identical(sampleid_policy, "drop")) {
    warning(msg, "\nDropping offending rows and continuing.")
    df_raw <- df_raw[-bad_id, , drop = FALSE]
  } else if (identical(sampleid_policy, "fill")) {
    warning(msg, "\nFilling missing SampleID with deterministic placeholders and continuing.")
    df_raw$SampleID[bad_id] <- paste0("MISSING_SAMPLEID_ROW_", bad_id)
  } else {
    stop("sampleid_policy must be one of: 'error', 'drop', 'fill'.")
  }
}

df_raw$Day <- suppressWarnings(as.numeric(as.character(df_raw$Day)))
if (any(!is.finite(df_raw$Day))) {
  bad <- which(!is.finite(df_raw$Day))[1:min(10, sum(!is.finite(df_raw$Day)))]
  stop("Day contains non-numeric / NA after coercion. Example rows: ", paste(bad, collapse = ", "))
}

df_raw$Abundance <- suppressWarnings(as.numeric(as.character(df_raw$Abundance)))
if (any(!is.finite(df_raw$Abundance))) {
  bad <- which(!is.finite(df_raw$Abundance))[1:min(10, sum(!is.finite(df_raw$Abundance)))]
  stop("Abundance contains non-numeric / NA after coercion. Example rows: ", paste(bad, collapse = ", "))
}

if (any(df_raw$Source == "")) stop("Source has empty strings after trimws().")
if (any(df_raw$Metabolite == "")) stop("Metabolite has empty strings after trimws().")

tidy <- df_raw

# -----------------------------
# Optional transform
# -----------------------------
tidy$Abundance_use <- tidy$Abundance
transform_mode <- "raw"
if (isTRUE(log_transform)) {
  cat("[INFO] Applying log2(x +", log_pseudocount, ") transform to Abundance.\n")
  tidy$Abundance_use <- log2(tidy$Abundance_use + log_pseudocount)
  transform_mode <- "log2p1"
}

cat("\nTidy input accepted.\n")
cat("Rows:", nrow(tidy),
    " | Unique Sources:", length(unique(tidy$Source)),
    " | Unique Days:", length(unique(tidy$Day)),
    " | Unique Features:", length(unique(tidy$Metabolite)), "\n\n")

sources  <- sort(unique(tidy$Source))
features <- sort(unique(tidy$Metabolite))

# -----------------------------
# Output directory policy (mandatory)
# -----------------------------
outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

source_tag <- if (length(sources) == 1) safe_name(sources[1]) else "MultiSource"
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0(safe_name(analysis_name), "_", source_tag, "_", run_timestamp)

output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Standardized prefix inside output_dir
out_prefix <- file.path(output_dir, paste0("MeanByDay_", run_timestamp))
set.seed(set_seed)

# Console summary (mandatory)
cat("------------------------------------------------------------\n")
cat("Script identity:\n")
cat("  script_name: ", script_name, "\n", sep = "")
cat("  script_path: ", script_path, "\n", sep = "")
cat("  script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
cat("Run:\n")
cat("  run_id:      ", run_id, "\n", sep = "")
cat("  output_dir:  ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Input:\n")
cat("  infile:      ", infile, "\n", sep = "")
cat("------------------------------------------------------------\n\n")

# ======================================================================
# Levene + permutation hooks (Source × Metabolite) using replicate-level tidy
# ======================================================================
levene_perm_df <- data.frame()

levene_F_median <- function(values, groups) {
  ok <- is.finite(values) & is.finite(groups)
  values <- values[ok]
  groups <- groups[ok]
  groups <- as.factor(groups)
  
  k <- nlevels(groups)
  if (k < 2) return(list(F = NA_real_, df1 = NA_integer_, df2 = NA_integer_, p = NA_real_))
  
  meds <- tapply(values, groups, median, na.rm = TRUE)
  z <- abs(values - meds[as.character(groups)])
  
  a <- aov(z ~ groups)
  s <- summary(a)[[1]]
  F <- as.numeric(s["groups", "F value"])
  df1 <- as.integer(s["groups", "Df"])
  df2 <- as.integer(s["Residuals", "Df"])
  p <- as.numeric(s["groups", "Pr(>F)"])
  
  list(F = F, df1 = df1, df2 = df2, p = p)
}

cat("[INFO] Starting Levene permutations...\n")

if (isTRUE(do_levene_perm) && nrow(tidy) > 0) {
  
  key_vec <- paste(tidy$Source, tidy$Metabolite, sep = "||")
  idx_list <- split(seq_len(nrow(tidy)), key_vec)
  keys <- names(idx_list)
  
  # Reproducible parallel RNG
  RNGkind("L'Ecuyer-CMRG")
  set.seed(set_seed)
  
  worker_fun <- function(ii) {
    idx <- idx_list[[ii]]
    
    src <- tidy$Source[idx[1]]
    met <- tidy$Metabolite[idx[1]]
    
    v <- tidy$Abundance_use[idx]
    g <- tidy$Day[idx]
    
    ok <- is.finite(v) & is.finite(g)
    v <- v[ok]; g <- g[ok]
    
    if (length(v) < 3) return(NULL)
    if (length(unique(g)) < 2) return(NULL)
    
    obs <- levene_F_median(v, g)
    
    perm_p <- NA_real_
    if (is.finite(obs$F) && n_perm_levene >= 100) {
      Fnull <- numeric(n_perm_levene)
      for (b in seq_len(n_perm_levene)) {
        gp <- sample(g, length(g), replace = FALSE)
        tmp <- levene_F_median(v, gp)
        Fnull[b] <- tmp$F
      }
      perm_p <- (sum(Fnull >= obs$F, na.rm = TRUE) + 1) / (sum(is.finite(Fnull)) + 1)
    }
    
    data.frame(
      Source = src,
      Metabolite = met,
      N_total = length(v),
      N_days = length(unique(g)),
      Levene_F = obs$F,
      Levene_df1 = obs$df1,
      Levene_df2 = obs$df2,
      Levene_p = obs$p,
      Perm_p = perm_p,
      N_perm = if (is.finite(perm_p)) n_perm_levene else 0L,
      Levene_Definition = "Median-centered Levene across days; replicate-level Abundance_use",
      Transform = transform_mode,
      stringsAsFactors = FALSE
    )
  }
  
  if (isTRUE(use_parallel) && n_cores > 1) {
    res_list <- parallel::mclapply(seq_along(keys), worker_fun, mc.cores = n_cores)
  } else {
    res_list <- lapply(seq_along(keys), worker_fun)
  }
  
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  levene_perm_df <- if (length(res_list) > 0) do.call(rbind, res_list) else data.frame()
}

out_levene_perm <- paste0(out_prefix, "_LevenePerm_ByFeature.csv")
write.csv(levene_perm_df, out_levene_perm, row.names = FALSE)
cat("[INFO] Finished Levene permutations.\n")

# ======================================================================
# 1) Per-Source × Feature × Day summaries
# ======================================================================
sum_rows <- list()
si <- 1

for (src in sources) {
  sub_src <- tidy[tidy$Source == src, , drop = FALSE]
  days <- sort(unique(sub_src$Day))
  
  for (met in features) {
    sub_met <- sub_src[sub_src$Metabolite == met, , drop = FALSE]
    if (nrow(sub_met) == 0) next
    
    for (d in days) {
      v <- sub_met$Abundance_use[sub_met$Day == d]
      v <- v[is.finite(v)]
      if (length(v) == 0) next
      
      stats <- mean_ci95(v)
      
      sum_rows[[si]] <- data.frame(
        Source = src,
        Metabolite = met,
        Day = d,
        N = as.integer(stats["n"]),
        Mean = as.numeric(stats["mean"]),
        SD = as.numeric(stats["sd"]),
        SE = as.numeric(stats["se"]),
        CI_Low = as.numeric(stats["ci_low"]),
        CI_High = as.numeric(stats["ci_high"]),
        Transform = transform_mode,
        stringsAsFactors = FALSE
      )
      si <- si + 1
    }
  }
}

mean_byday_df <- if (length(sum_rows) > 0) do.call(rbind, sum_rows) else data.frame()
out_mean_byday <- paste0(out_prefix, "_Feature_ByDay_Mean_WithError.csv")
write.csv(mean_byday_df, out_mean_byday, row.names = FALSE)

# ======================================================================
# CV-only table: Source × Metabolite × Day (biological replicate CV)
# ======================================================================
cv_only_df <- data.frame()
if (nrow(mean_byday_df) > 0) {
  cv_only_df <- mean_byday_df[, c("Source","Metabolite","Day","N","Mean","SD","Transform")]
  cv_only_df$CV <- ifelse(is.finite(cv_only_df$SD) & is.finite(cv_only_df$Mean) & cv_only_df$Mean != 0,
                          cv_only_df$SD / abs(cv_only_df$Mean),
                          NA_real_)
  cv_only_df$CV_Definition <- "CV = SD/|Mean| across biological replicates within day"
}

out_cv_only <- paste0(out_prefix, "_CVONLY_FeatureByDay_BioReplicates.csv")
write.csv(cv_only_df, out_cv_only, row.names = FALSE)

# ======================================================================
# 1B) Group-mean CV by Day
# ======================================================================
group_mean_cv_df <- data.frame()

if (exists("cv_only_df") && nrow(cv_only_df) > 0) {
  
  cv_ok <- is.finite(cv_only_df$CV) & is.finite(cv_only_df$Day) &
    is.finite(cv_only_df$N) & cv_only_df$N >= min_n_per_day_for_stats
  
  tmp <- cv_only_df[cv_ok, , drop = FALSE]
  
  if (nrow(tmp) > 0) {
    
    g_list <- list()
    gi <- 1
    
    for (src in sort(unique(tmp$Source))) {
      sub <- tmp[tmp$Source == src, , drop = FALSE]
      days <- sort(unique(sub$Day))
      
      for (d in days) {
        x <- sub$CV[sub$Day == d]
        x <- x[is.finite(x)]
        
        if (length(x) < min_n_per_day_for_stats) next
        
        st <- mean_ci95(x)
        
        g_list[[gi]] <- data.frame(
          Source = src,
          Day = d,
          N_features_with_CV = as.integer(st["n"]),
          Mean_CV = as.numeric(st["mean"]),
          SD_CV = as.numeric(st["sd"]),
          SE_CV = as.numeric(st["se"]),
          CI_Low_CV = as.numeric(st["ci_low"]),
          CI_High_CV = as.numeric(st["ci_high"]),
          CV_Definition = "Mean(CV across metabolites) per day; CV from bio-reps within day",
          Transform = transform_mode,
          stringsAsFactors = FALSE
        )
        gi <- gi + 1
      }
    }
    
    if (length(g_list) > 0) group_mean_cv_df <- do.call(rbind, g_list)
  }
}

out_group_mean_cv <- paste0(out_prefix, "_SourceByDay_GroupMeanCV_WithCI.csv")
write.csv(group_mean_cv_df, out_group_mean_cv, row.names = FALSE)

# ======================================================================
# 1B2) Dispersion band for Mean_CV vs Day (across metabolites)
# ======================================================================
group_cv_band_df <- data.frame()

if (exists("cv_only_df") && nrow(cv_only_df) > 0) {
  
  cv_ok <- is.finite(cv_only_df$CV) & is.finite(cv_only_df$Day) &
    is.finite(cv_only_df$N) & cv_only_df$N >= min_n_per_day_for_stats
  
  tmp <- cv_only_df[cv_ok, c("Source","Day","CV"), drop = FALSE]
  tmp$Day <- as.numeric(tmp$Day)
  
  if (nrow(tmp) > 0) {
    
    b_list <- list()
    bi <- 1
    
    for (src in sort(unique(tmp$Source))) {
      sub <- tmp[tmp$Source == src, , drop = FALSE]
      days <- sort(unique(sub$Day))
      
      for (d in days) {
        x <- sub$CV[sub$Day == d]
        x <- x[is.finite(x)]
        if (length(x) < min_n_per_day_for_stats) next
        
        q25 <- as.numeric(quantile(x, probs = 0.25, names = FALSE, type = 7, na.rm = TRUE))
        q75 <- as.numeric(quantile(x, probs = 0.75, names = FALSE, type = 7, na.rm = TRUE))
        q10 <- as.numeric(quantile(x, probs = 0.10, names = FALSE, type = 7, na.rm = TRUE))
        q90 <- as.numeric(quantile(x, probs = 0.90, names = FALSE, type = 7, na.rm = TRUE))
        
        b_list[[bi]] <- data.frame(
          Source = src,
          Day = d,
          N_features_with_CV = length(x),
          Q10_CV = q10,
          Q25_CV = q25,
          Q75_CV = q75,
          Q90_CV = q90,
          Band_Definition = "Quantile band across metabolites of feature-level CV within day",
          CV_Definition = "Feature CV = SD/|Mean| across biological replicates within day",
          Transform = transform_mode,
          stringsAsFactors = FALSE
        )
        bi <- bi + 1
      }
    }
    
    if (length(b_list) > 0) group_cv_band_df <- do.call(rbind, b_list)
  }
}

out_group_cv_band <- paste0(out_prefix, "_SourceByDay_GroupMeanCV_DispersionBand.csv")
write.csv(group_cv_band_df, out_group_cv_band, row.names = FALSE)

# ======================================================================
# System-wide permutation test for CV peaks/troughs (Source-level)
# ======================================================================
system_cv_perm_results <- data.frame()

if (exists("cv_only_df") && nrow(cv_only_df) > 0 &&
    isTRUE(exists("do_system_cv_perm")) && isTRUE(do_system_cv_perm)) {
  
  cv_ok <- is.finite(cv_only_df$CV) & is.finite(cv_only_df$Day) &
    is.finite(cv_only_df$N) & cv_only_df$N >= min_n_per_day_for_stats
  
  base <- cv_only_df[cv_ok, c("Source","Metabolite","Day","CV"), drop = FALSE]
  base$Day <- as.numeric(base$Day)
  
  if (nrow(base) > 0) {
    
    res_list <- list()
    rj <- 1
    
    for (src in sort(unique(base$Source))) {
      
      sub <- base[base$Source == src, , drop = FALSE]
      if (nrow(sub) < 5) next
      
      days_src <- sort(unique(sub$Day))
      
      obs_mean_byday <- tapply(sub$CV, sub$Day, mean, na.rm = TRUE)
      obs_mean_byday <- obs_mean_byday[as.character(days_src)]
      obs_peak  <- max(obs_mean_byday, na.rm = TRUE)
      obs_trough <- min(obs_mean_byday, na.rm = TRUE)
      
      obs_peak_focus <- NA_real_
      obs_trough_focus <- NA_real_
      if (is.finite(focus_peak_day) && focus_peak_day %in% days_src) {
        obs_peak_focus <- obs_mean_byday[as.character(focus_peak_day)]
      }
      if (is.finite(focus_trough_day) && focus_trough_day %in% days_src) {
        obs_trough_focus <- obs_mean_byday[as.character(focus_trough_day)]
      }
      
      keys <- unique(sub$Metabolite)
      
      peak_null <- numeric(n_perm_system_cv)
      trough_null <- numeric(n_perm_system_cv)
      peak_focus_null <- if (is.finite(obs_peak_focus)) numeric(n_perm_system_cv) else NULL
      trough_focus_null <- if (is.finite(obs_trough_focus)) numeric(n_perm_system_cv) else NULL
      
      for (b in seq_len(n_perm_system_cv)) {
        
        perm_day <- sub$Day
        
        for (met in keys) {
          idx <- which(sub$Metabolite == met)
          if (length(idx) < 2) next
          perm_day[idx] <- sample(sub$Day[idx], length(idx), replace = FALSE)
        }
        
        perm_mean <- tapply(sub$CV, perm_day, mean, na.rm = TRUE)
        perm_mean <- perm_mean[as.character(days_src)]
        peak_null[b] <- max(perm_mean, na.rm = TRUE)
        trough_null[b] <- min(perm_mean, na.rm = TRUE)
        
        if (!is.null(peak_focus_null)) {
          peak_focus_null[b] <- perm_mean[as.character(focus_peak_day)]
        }
        if (!is.null(trough_focus_null)) {
          trough_focus_null[b] <- perm_mean[as.character(focus_trough_day)]
        }
      }
      
      p_peak_global <- (sum(peak_null >= obs_peak, na.rm = TRUE) + 1) / (sum(is.finite(peak_null)) + 1)
      p_trough_global <- (sum(trough_null <= obs_trough, na.rm = TRUE) + 1) / (sum(is.finite(trough_null)) + 1)
      
      p_peak_focus <- NA_real_
      p_trough_focus <- NA_real_
      if (!is.null(peak_focus_null)) {
        p_peak_focus <- (sum(peak_focus_null >= obs_peak_focus, na.rm = TRUE) + 1) / (sum(is.finite(peak_focus_null)) + 1)
      }
      if (!is.null(trough_focus_null)) {
        p_trough_focus <- (sum(trough_focus_null <= obs_trough_focus, na.rm = TRUE) + 1) / (sum(is.finite(trough_focus_null)) + 1)
      }
      
      res_list[[rj]] <- data.frame(
        Source = src,
        N_perm = n_perm_system_cv,
        Observed_GlobalPeak_MeanCV = obs_peak,
        P_GlobalPeak = p_peak_global,
        Observed_GlobalTrough_MeanCV = obs_trough,
        P_GlobalTrough = p_trough_global,
        FocusPeakDay = focus_peak_day,
        Observed_FocusPeakDay_MeanCV = obs_peak_focus,
        P_FocusPeakDay = p_peak_focus,
        FocusTroughDay = focus_trough_day,
        Observed_FocusTroughDay_MeanCV = obs_trough_focus,
        P_FocusTroughDay = p_trough_focus,
        Null = "Permute Day labels within Source×Metabolite CV trajectories",
        Transform = transform_mode,
        stringsAsFactors = FALSE
      )
      rj <- rj + 1
    }
    
    if (length(res_list) > 0) system_cv_perm_results <- do.call(rbind, res_list)
  }
}

out_system_perm <- paste0(out_prefix, "_SystemWide_CVPeakTrough_PermutationTest.csv")
write.csv(system_cv_perm_results, out_system_perm, row.names = FALSE)

# ======================================================================
# 2) Per-Source “average feature mean” by Day
# ======================================================================
avg_rows <- list()
ai <- 1

if (nrow(mean_byday_df) > 0) {
  for (src in sources) {
    sub <- mean_byday_df[mean_byday_df$Source == src, , drop = FALSE]
    if (nrow(sub) == 0) next
    days <- sort(unique(sub$Day))
    
    for (d in days) {
      x <- sub$Mean[sub$Day == d]
      x <- x[is.finite(x)]
      if (length(x) == 0) next
      
      st <- mean_ci95(x)
      
      avg_rows[[ai]] <- data.frame(
        Source = src,
        Day = d,
        N_features = as.integer(st["n"]),
        Mean = as.numeric(st["mean"]),
        SD = as.numeric(st["sd"]),
        SE = as.numeric(st["se"]),
        CI_Low = as.numeric(st["ci_low"]),
        CI_High = as.numeric(st["ci_high"]),
        Transform = transform_mode,
        stringsAsFactors = FALSE
      )
      ai <- ai + 1
    }
  }
}

avg_df <- if (length(avg_rows) > 0) do.call(rbind, avg_rows) else data.frame()
out_avg_byday <- paste0(out_prefix, "_SourceByDay_AverageFeatureMean_WithError.csv")
write.csv(avg_df, out_avg_byday, row.names = FALSE)

# ======================================================================
# 3) Replicate-level "system mean" per Source×Day×SampleID
# ======================================================================
rep_sys <- data.frame()
if (nrow(tidy) > 0) {
  keys <- unique(paste(tidy$Source, tidy$Day, tidy$SampleID, sep = "||"))
  rep_list <- vector("list", length(keys))
  ri <- 1
  
  for (k in keys) {
    parts <- strsplit(k, "\\|\\|", fixed = FALSE)[[1]]
    src <- parts[1]; d <- as.numeric(parts[2]); sid <- parts[3]
    
    sub <- tidy[tidy$Source == src & tidy$Day == d & tidy$SampleID == sid, , drop = FALSE]
    v <- sub$Abundance_use
    v <- v[is.finite(v)]
    if (length(v) == 0) next
    
    rep_list[[ri]] <- data.frame(
      Source = src,
      Day = d,
      SampleID = sid,
      SystemMean = mean(v),
      N_features_in_sample = length(v),
      Transform = transform_mode,
      stringsAsFactors = FALSE
    )
    ri <- ri + 1
  }
  
  rep_list <- rep_list[!vapply(rep_list, is.null, logical(1))]
  if (length(rep_list) > 0) rep_sys <- do.call(rbind, rep_list)
}

out_rep_sys <- paste0(out_prefix, "_ReplicateSystemMean_BySourceDaySample.csv")
write.csv(rep_sys, out_rep_sys, row.names = FALSE)

# ======================================================================
# 3B) CV of system mean across biological replicates (Source × Day)
# ======================================================================
sysmean_cv_df <- data.frame()

if (nrow(rep_sys) > 0) {
  
  sm_list <- list()
  si2 <- 1
  
  for (src in sort(unique(rep_sys$Source))) {
    sub <- rep_sys[rep_sys$Source == src, , drop = FALSE]
    days <- sort(unique(sub$Day))
    
    for (d in days) {
      x <- sub$SystemMean[sub$Day == d]
      x <- x[is.finite(x)]
      if (length(x) < min_n_per_day_for_stats) next
      
      st <- mean_ci95(x)
      mu <- as.numeric(st["mean"])
      s  <- as.numeric(st["sd"])
      
      cv <- if (is.finite(s) && is.finite(mu) && mu != 0) s / abs(mu) else NA_real_
      
      sm_list[[si2]] <- data.frame(
        Source = src,
        Day = d,
        N_reps = as.integer(st["n"]),
        Mean_SystemMean = mu,
        SD_SystemMean = s,
        CV_SystemMean = cv,
        SE_SystemMean = as.numeric(st["se"]),
        CI_Low_SystemMean = as.numeric(st["ci_low"]),
        CI_High_SystemMean = as.numeric(st["ci_high"]),
        CV_Definition = "CV = SD/|Mean| of replicate SystemMean within Source×Day",
        Transform = transform_mode,
        stringsAsFactors = FALSE
      )
      si2 <- si2 + 1
    }
  }
  
  if (length(sm_list) > 0) sysmean_cv_df <- do.call(rbind, sm_list)
}

out_sysmean_cv <- paste0(out_prefix, "_SourceByDay_SystemMean_CV_AcrossReplicates.csv")
write.csv(sysmean_cv_df, out_sysmean_cv, row.names = FALSE)

# ======================================================================
# 4) Global ALL-sources replicate system-mean by Day + band
# ======================================================================
glob_rep_df <- data.frame()
band_rep_df <- data.frame()

if (nrow(rep_sys) > 0) {
  days_all <- sort(unique(rep_sys$Day))
  
  if (identical(global_band_type, "P10P90")) { probs_lo <- 0.10; probs_hi <- 0.90 } else { probs_lo <- 0.25; probs_hi <- 0.75 }
  
  gr_list <- list()
  br_list <- list()
  gi <- 1
  bi <- 1
  
  for (d in days_all) {
    x <- rep_sys$SystemMean[rep_sys$Day == d]
    x <- x[is.finite(x)]
    if (length(x) == 0) next
    
    center_val <- if (identical(global_center, "median")) median(x) else mean(x)
    st <- mean_ci95(x)
    
    gr_list[[gi]] <- data.frame(
      Day = d,
      N_reps = as.integer(st["n"]),
      Center = center_val,
      Mean = as.numeric(st["mean"]),
      SD = as.numeric(st["sd"]),
      SE = as.numeric(st["se"]),
      CI_Low = as.numeric(st["ci_low"]),
      CI_High = as.numeric(st["ci_high"]),
      Transform = transform_mode,
      stringsAsFactors = FALSE
    )
    gi <- gi + 1
    
    if (length(x) >= 2) {
      qlo <- as.numeric(quantile(x, probs = probs_lo, names = FALSE, type = 7, na.rm = TRUE))
      qhi <- as.numeric(quantile(x, probs = probs_hi, names = FALSE, type = 7, na.rm = TRUE))
    } else {
      qlo <- NA_real_; qhi <- NA_real_
    }
    
    br_list[[bi]] <- data.frame(
      Day = d,
      Band_Low = qlo,
      Band_High = qhi,
      N_reps = length(x),
      Band_Type = global_band_type,
      stringsAsFactors = FALSE
    )
    bi <- bi + 1
  }
  
  if (length(gr_list) > 0) glob_rep_df <- do.call(rbind, gr_list)
  if (length(br_list) > 0) band_rep_df <- do.call(rbind, br_list)
}

out_global <- paste0(out_prefix, "_ALLSources_ReplicateSystemMean_ByDay_WithError.csv")
write.csv(glob_rep_df, out_global, row.names = FALSE)

out_band <- paste0(out_prefix, "_ALLSources_ReplicateBand_ByDay.csv")
write.csv(band_rep_df, out_band, row.names = FALSE)

# ======================================================================
# Always-on static key plot (for report include_graphics requirement)
#   - Generates a ggplot2 PNG even if base plotting is disabled
# ======================================================================
keyplot_path <- file.path(output_dir, "KeyPlot_GlobalSystemMean.png")
try({
  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(glob_rep_df) >= 2) {
    gdf <- glob_rep_df[is.finite(glob_rep_df$Day) & is.finite(glob_rep_df$Center), , drop = FALSE]
    bdf <- band_rep_df
    p <- ggplot2::ggplot(gdf, ggplot2::aes(x = Day, y = Center)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
        title = "Global replicate system mean by day (ALL Sources)",
        x = "Day",
        y = "Replicate system mean"
      )
    
    # Optional band
    if (isTRUE(add_global_band) && nrow(bdf) > 0) {
      okb <- is.finite(bdf$Day) & is.finite(bdf$Band_Low) & is.finite(bdf$Band_High)
      b2 <- bdf[okb, , drop = FALSE]
      if (nrow(b2) >= 2) {
        p <- p + ggplot2::geom_ribbon(
          data = b2,
          ggplot2::aes(x = Day, ymin = Band_Low, ymax = Band_High),
          inherit.aes = FALSE,
          alpha = 0.3
        )
      }
    }
    
    ggplot2::ggsave(filename = keyplot_path, plot = p, width = 9, height = 5, dpi = 150)
  } else {
    # Create a minimal placeholder image so include_graphics always has a target
    png(keyplot_path, width = 1200, height = 700, res = 150)
    plot.new()
    text(0.5, 0.5, "Key plot unavailable (insufficient data).")
    dev.off()
  }
}, silent = TRUE)

# ======================================================================
# PLOTTING (base R; preserved logic)
# ======================================================================
if (isTRUE(make_plots)) {
  
  mk <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE) }
  
  plot_root <- file.path(output_dir, "plots")
  mk(plot_root)
  
  dir_feat <- file.path(plot_root, "01_FeatureMean_Trajectories")
  dir_avg  <- file.path(plot_root, "02_AverageFeatureMean")
  mk(dir_feat); mk(dir_avg)
  
  # ---- 2E) Mean_CV vs Day with metabolite-dispersion band (IQR) ----
  if (exists("group_mean_cv_df") && nrow(group_mean_cv_df) > 0 &&
      exists("group_cv_band_df") && nrow(group_cv_band_df) > 0) {
    
    for (src in sort(unique(group_mean_cv_df$Source))) {
      
      m <- group_mean_cv_df[group_mean_cv_df$Source == src, , drop = FALSE]
      b <- group_cv_band_df[group_cv_band_df$Source == src, , drop = FALSE]
      
      m <- m[is.finite(m$Day) & is.finite(m$Mean_CV), , drop = FALSE]
      b <- b[is.finite(b$Day) & is.finite(b$Q25_CV) & is.finite(b$Q75_CV), , drop = FALSE]
      
      if (nrow(m) < 2 || nrow(b) < 2) next
      
      m <- m[order(m$Day), , drop = FALSE]
      b <- b[order(b$Day), , drop = FALSE]
      
      out_file <- file.path(dir_avg, paste0("MeanCV_ByDay_WithIQRBand_", safe_name(src), ".png"))
      png(out_file, width = plot_width, height = plot_height, res = plot_res)
      
      ylim_vals <- c(m$Mean_CV, b$Q25_CV, b$Q75_CV)
      ylim_vals <- ylim_vals[is.finite(ylim_vals)]
      if (length(ylim_vals) < 2) { dev.off(); next }
      
      plot(m$Day, m$Mean_CV, type = "n",
           main = paste0("Mean_CV vs Day with IQR band (", src, ")"),
           xlab = "Day",
           ylab = "Mean(CV across metabolites); band = IQR across metabolites",
           ylim = range(ylim_vals))
      
      ord <- order(b$Day)
      xd <- b$Day[ord]
      lo <- b$Q25_CV[ord]
      hi <- b$Q75_CV[ord]
      polygon(c(xd, rev(xd)), c(lo, rev(hi)), col = "gray85", border = NA)
      
      lines(m$Day, m$Mean_CV, lwd = 2)
      points(m$Day, m$Mean_CV, pch = 16)
      
      dev.off()
    }
  }
  
  ylab_feat <- if (identical(transform_mode, "log2p1")) "Mean (log2(x+pseudocount))" else "Mean (Abundance)"
  
  # ---- 2F) SystemMean CV across replicates by Source ----
  if (exists("sysmean_cv_df") && nrow(sysmean_cv_df) > 0) {
    
    out_file <- file.path(dir_avg, "SystemMeanCV_AcrossReplicates_BySource.png")
    png(out_file, width = plot_width, height = plot_height, res = plot_res)
    
    dd <- sysmean_cv_df[
      is.finite(sysmean_cv_df$Day) & is.finite(sysmean_cv_df$CV_SystemMean),
      , drop = FALSE
    ]
    dd <- dd[order(dd$Source, dd$Day), , drop = FALSE]
    
    if (nrow(dd) >= 2) {
      
      ylim_vals <- dd$CV_SystemMean
      ylim_vals <- ylim_vals[is.finite(ylim_vals)]
      
      plot(NA,
           xlim = range(dd$Day),
           ylim = range(ylim_vals),
           xlab = "Day",
           ylab = "CV of replicate system mean (SD/|Mean| across SampleIDs)",
           main = "Between-individual heterogeneity: CV(SystemMean) by day (by Source)")
      
      srcs <- sort(unique(dd$Source))
      cols <- seq_along(srcs)
      
      for (i in seq_along(srcs)) {
        s <- dd[dd$Source == srcs[i], , drop = FALSE]
        s <- s[order(s$Day), , drop = FALSE]
        if (nrow(s) < 2) next
        
        lines(s$Day, s$CV_SystemMean, lwd = 2, col = cols[i])
        points(s$Day, s$CV_SystemMean, pch = 16, col = cols[i])
      }
      
      legend("topright", legend = srcs, col = cols, lwd = 2, pch = 16, bty = "n")
    }
    
    dev.off()
  }
  
  # ---- 1) Feature mean trajectories (top N by overall mean within each Source) ----
  if (nrow(mean_byday_df) > 0) {
    for (src in sources) {
      sdf <- mean_byday_df[mean_byday_df$Source == src, , drop = FALSE]
      if (nrow(sdf) == 0) next
      
      feat_means <- tapply(sdf$Mean, sdf$Metabolite, function(z) mean(z[is.finite(z)], na.rm = TRUE))
      feat_means <- feat_means[is.finite(feat_means)]
      if (length(feat_means) == 0) next
      
      keep_feats <- head(names(sort(feat_means, decreasing = TRUE)), top_n_features_per_source)
      
      sdir <- file.path(dir_feat, safe_name(src))
      mk(sdir)
      
      for (met in keep_feats) {
        dd <- sdf[sdf$Metabolite == met, , drop = FALSE]
        dd <- dd[order(dd$Day), , drop = FALSE]
        okp <- is.finite(dd$Day) & is.finite(dd$Mean)
        if (sum(okp) < 2) next
        
        out_file <- file.path(sdir, paste0("MeanTrajectory_", safe_name(met), ".png"))
        plot_feature(dd, paste0("Mean trajectory: ", met, " (", src, ")"),
                     ylab_feat, out_file, plot_width, plot_height, plot_res)
      }
    }
  }
  
  # ---- 2) Per-Source average feature mean by day (with CI) ----
  if (nrow(avg_df) > 0) {
    for (src in sources) {
      s <- avg_df[avg_df$Source == src, , drop = FALSE]
      s <- s[order(s$Day), , drop = FALSE]
      okp <- is.finite(s$Day) & is.finite(s$Mean)
      if (sum(okp) < 2) next
      s <- s[okp, , drop = FALSE]
      
      out_file <- file.path(dir_avg, paste0("AverageFeatureTrajectory_", safe_name(src), ".png"))
      png(out_file, width = plot_width, height = plot_height, res = plot_res)
      
      ylim_vals <- s$Mean
      okci <- is.finite(s$CI_Low) & is.finite(s$CI_High)
      if (any(okci)) ylim_vals <- c(ylim_vals, s$CI_Low[okci], s$CI_High[okci])
      ylim_vals <- ylim_vals[is.finite(ylim_vals)]
      if (length(ylim_vals) < 2) { dev.off(); next }
      
      plot(s$Day, s$Mean, type = "n",
           main = paste0("Average feature mean by day (", src, ")"),
           xlab = "Day",
           ylab = if (identical(transform_mode, "log2p1")) "Mean of feature means (log2 scale)" else "Mean of feature means",
           ylim = range(ylim_vals))
      
      if (any(okci)) segments(s$Day[okci], s$CI_Low[okci], s$Day[okci], s$CI_High[okci])
      lines(s$Day, s$Mean, lwd = 2)
      points(s$Day, s$Mean, pch = 16)
      
      dev.off()
    }
  }
  
  # ---- 3) Global ALL Sources plot (grey band + CI + line) ----
  if (nrow(glob_rep_df) > 0) {
    
    gplot <- glob_rep_df[order(glob_rep_df$Day), , drop = FALSE]
    bplot <- band_rep_df[order(band_rep_df$Day), , drop = FALSE]
    
    okp <- is.finite(gplot$Day) & is.finite(gplot$Center)
    if (sum(okp) < 2) {
      warning("Global ALL-sources plot skipped: fewer than 2 finite day points.")
    } else {
      
      out_file <- file.path(dir_avg, "AverageFeatureTrajectory_ALLSources.png")
      png(out_file, width = plot_width, height = plot_height, res = plot_res)
      
      okband_idx <- rep(FALSE, nrow(bplot))
      if (isTRUE(add_global_band) && nrow(bplot) > 0) {
        okband_idx <- is.finite(bplot$Day) & is.finite(bplot$Band_Low) & is.finite(bplot$Band_High)
      }
      okband_any <- any(okband_idx)
      
      okci <- is.finite(gplot$CI_Low) & is.finite(gplot$CI_High) & is.finite(gplot$Day) & is.finite(gplot$Mean)
      
      ylim_vals <- gplot$Center[okp]
      if (okband_any) ylim_vals <- c(ylim_vals, bplot$Band_Low[okband_idx], bplot$Band_High[okband_idx])
      if (any(okci))  ylim_vals <- c(ylim_vals, gplot$CI_Low[okci], gplot$CI_High[okci])
      ylim_vals <- ylim_vals[is.finite(ylim_vals)]
      
      if (length(ylim_vals) < 2) {
        dev.off()
      } else {
        
        plot(gplot$Day[okp], gplot$Center[okp], type = "n",
             main = "Average feature mean by day (ALL Sources; replicate band)",
             xlab = "Day",
             ylab = if (identical(transform_mode, "log2p1")) "Replicate system mean (log2 scale)" else "Replicate system mean",
             ylim = range(ylim_vals))
        
        if (okband_any) {
          ord <- order(bplot$Day[okband_idx])
          xd <- bplot$Day[okband_idx][ord]
          lo <- bplot$Band_Low[okband_idx][ord]
          hi <- bplot$Band_High[okband_idx][ord]
          polygon(c(xd, rev(xd)), c(lo, rev(hi)), col = band_fill_col, border = NA)
        }
        
        if (any(okci)) segments(gplot$Day[okci], gplot$CI_Low[okci], gplot$Day[okci], gplot$CI_High[okci])
        
        lines(gplot$Day[okp], gplot$Center[okp], lwd = 2)
        points(gplot$Day[okp], gplot$Center[okp], pch = 16)
        
        dev.off()
      }
    }
  }
  
  cat("Plots saved to:\n  ", normalizePath(plot_root, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
}

# ======================================================================
# Manifest (enhanced): write initial manifest + file inventory (pre-render)
# ======================================================================
manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
manifest_files_csv <- file.path(output_dir, "Project_Manifest_Files.csv")

deps_list <- pkg_versions(c(
  "jsonlite", "rmarkdown", "knitr", "ggplot2",
  if (isTRUE(enable_plotly)) "plotly" else NULL,
  "parallel"
))

execution_mode <- if (interactive()) "interactive" else "script"
r_info <- R.version

manifest <- list(
  run_id = run_id,
  run_timestamp = as.character(Sys.time()),
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), "NA", script_full),
    execution_mode = execution_mode,
    command_args = commandArgs(trailingOnly = TRUE)
  ),
  system = list(
    R_version = r_info$version.string,
    platform = r_info$platform
  ),
  input = list(
    input_csv = infile,
    required_columns = required_cols,
    sampleid_policy = sampleid_policy
  ),
  parameters = list(
    analysis_name = analysis_name,
    set_seed = set_seed,
    log_transform = log_transform,
    log_pseudocount = log_pseudocount,
    make_plots = make_plots,
    top_n_features_per_source = top_n_features_per_source,
    min_n_per_day_for_stats = min_n_per_day_for_stats,
    add_global_band = add_global_band,
    global_band_type = global_band_type,
    global_center = global_center,
    band_fill_col = band_fill_col,
    do_levene_perm = do_levene_perm,
    n_perm_levene = n_perm_levene,
    do_system_cv_perm = do_system_cv_perm,
    n_perm_system_cv = n_perm_system_cv,
    focus_peak_day = focus_peak_day,
    focus_trough_day = focus_trough_day,
    use_parallel = use_parallel,
    n_cores = n_cores,
    enable_plotly = enable_plotly,
    enable_runtime_tracking = enable_runtime_tracking,
    runtime_seconds = NA_real_
  ),
  dependencies = deps_list,
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
)

jsonlite::write_json(manifest, manifest_json_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
write_file_inventory(output_dir, manifest_files_csv)

# ======================================================================
# Quarto QMD: Generate enhanced report scaffold
# ======================================================================
qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_name <- paste0(script_name, "_Report.html")
html_full <- file.path(output_dir, html_name)

header_block <- get_header_block()

qmd_lines <- c(
  "---",
  paste0('title: "', script_name, ' Report"'),
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 3",
  "    number-sections: true",
  "    theme: flatly",
  "    code-fold: true",
  "    df-print: paged",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "params:",
  paste0('  run_id: "', run_id, '"'),
  paste0('  run_timestamp: "', manifest$run_timestamp, '"'),
  paste0('  script_name: "', script_name, '"'),
  paste0('  script_path: "', script_path, '"'),
  paste0('  script_full: "', ifelse(is.na(script_full), "NA", script_full), '"'),
  paste0('  input_csv: "', infile, '"'),
  paste0('  output_dir: "', normalizePath(output_dir, winslash = "/", mustWork = FALSE), '"'),
  paste0('  runtime_seconds: NA'),
  paste0('  enable_plotly: ', ifelse(isTRUE(enable_plotly), "true", "false")),
  paste0('  enable_runtime_tracking: ', ifelse(isTRUE(enable_runtime_tracking), "true", "false")),
  "---",
  "",
  "# Script Header",
  "::: {.panel}",
  "```",
  header_block,
  "```",
  ":::",
  "",
  "# Metadata Overview",
  "",
  "```{r}",
  "meta <- data.frame(",
  "  key = c('run_id','run_timestamp','script_name','script_path','script_full','execution_mode','input_csv','output_dir','runtime_seconds','sampleid_policy','enable_plotly','enable_runtime_tracking'),",
  "  value = c(params$run_id, params$run_timestamp, params$script_name, params$script_path, params$script_full, ",
  "            ifelse(interactive(), 'interactive', 'script'), params$input_csv, params$output_dir, ",
  "            params$runtime_seconds, man$input$sampleid_policy, params$enable_plotly, params$enable_runtime_tracking),",
  "  stringsAsFactors = FALSE",
  ")",
  "knitr::kable(meta)",
  "```",
  "",
  "# System Information",
  "```{r}",
  "man <- jsonlite::fromJSON(file.path(params$output_dir, 'Project_Manifest.json'))",
  "sys_df <- as.data.frame(man$system, stringsAsFactors = FALSE)",
  "knitr::kable(t(sys_df), caption = 'R session and platform info')",
  "```",
  "",
  "# Dependencies",
  "```{r}",
  "dep_df <- as.data.frame(man$dependencies, stringsAsFactors = FALSE)",
  "knitr::kable(dep_df, caption = 'R package dependencies')",
  "```",
  "",
  "# Outputs",
  "```{r}",
  "inv <- read.csv(file.path(params$output_dir, 'Project_Manifest_Files.csv'), check.names = FALSE)",
  "knitr::kable(inv, caption = 'Files generated during run')",
  "```",
  "",
  "# Method Summary",
  "::: {.panel}",
  "This report summarizes trajectory and dispersion metrics from tidy longitudinal data using log2(x + pseudocount) transformation, biological replicate error propagation, and optional permutation tests for variance (Levene) and trajectory shape.",
  ":::",
  "",
  "# Key Static Diagnostic Plot",
  "```{r}",
  "kp <- file.path(params$output_dir, 'KeyPlot_GlobalSystemMean.png')",
  "if (file.exists(kp)) knitr::include_graphics(kp) else cat('Key plot not available.')",
  "```",
  "",
  "# Optional Interactive Plot",
  "```{r}",
  "if (isTRUE(params$enable_plotly)) {",
  "  library(plotly); library(ggplot2)",
  "  f <- list.files(params$output_dir, pattern = 'ALLSources_ReplicateSystemMean_ByDay_WithError.csv$', full.names = TRUE)",
  "  if (length(f) >= 1) {",
  "    gdf <- read.csv(f[1])",
  "    p <- ggplot(gdf, aes(x = Day, y = Center)) + geom_point() + geom_line() + labs(title = 'Global Replicate System Mean')",
  "    plotly::ggplotly(p)",
  "  } else cat('Interactive plot unavailable.')",
  "} else cat('Plotly disabled.')",
  "```",
  "",
  "# Session Info",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, qmd_path)

