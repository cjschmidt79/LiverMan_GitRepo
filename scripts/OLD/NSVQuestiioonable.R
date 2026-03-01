#!/usr/bin/env Rscript
# ============================================================
# 01_compute_NSV_trajectories.R
# NSV Starter Script (robust to unknown batch; no conditions)
#
# PURPOSE
#   Compute NSV trajectories under multiple nuisance models:
#     (1) RAW
#     (2) Residualize Day means (RESID_DAY)
#     (3) Residualize Day means + remove top k latent PCs (RESID_DAY_PCk)
#
# INPUT
#   Wide CSV with:
#     - Day column
#     - SampleID column
#     - Feature columns (genes/metabolites)
#
# TRANSFORM (NEW)
#   User choice:
#     0 = raw
#     1 = log2(TPM + pseudocount)
#
# OUTPUT POLICY (MANDATORY)
#   All outputs under:
#     outputs/<run_id>/
#   Including:
#     - NSV tables + plots
#     - Project_Manifest.json
#     - Project_Manifest_Files.csv
#     - <script_name>_Report.qmd + rendered HTML
#
# OPTIONAL (FLAGS; MANDATORY TO EXIST)
#   enable_plotly <- TRUE
#   enable_runtime_tracking <- TRUE
#
# NOTES
#   This script does not "prove" batch doesn't exist; it probes robustness to
#   plausible hidden technical structure (residualization + latent PC removal).
# ============================================================

suppressPackageStartupMessages({
  library(stats)
  library(utils)
})

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

# -------------------------- Internal feature flags (MANDATORY) --------------------------
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

# Runtime tracking start (MANDATORY if enabled)
start_time <- if (isTRUE(enable_runtime_tracking)) Sys.time() else NA

# ------------------------------- Helpers --------------------------------
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2    <- function(...) stop(paste0(...), call. = FALSE)
msg      <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# Diagnostics toggles (added)
enable_self_parse_diagnostic <- TRUE
write_warnings_file <- TRUE
debug_matrix_sizes <- TRUE

ensure_pkgs <- function(pkgs) {
  ip <- rownames(installed.packages())
  missing <- pkgs[!(pkgs %in% ip)]
  if (length(missing) > 0) {
    msg("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

arg_val <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

is_true <- function(x) {
  if (is.null(x) || is.na(x)) return(FALSE)
  x <- tolower(as.character(x))
  x %in% c("1", "true", "t", "yes", "y")
}

pick_file <- function(prompt = "Choose file") {
  if (interactive() && exists("file.choose")) {
    msg(prompt)
    return(file.choose())
  }
  stop2(prompt, " (non-interactive). Provide a path via command line args.")
}

ask_choice01 <- function(prompt, default = 0) {
  if (!interactive()) return(as.integer(default))
  ans <- readline(paste0(prompt, " [0/1; default ", default, "]: "))
  if (nchar(ans) == 0) return(as.integer(default))
  out <- suppressWarnings(as.integer(ans))
  if (is.na(out) || !(out %in% c(0, 1))) return(as.integer(default))
  out
}

ask_num <- function(prompt, default = 1) {
  if (!interactive()) return(as.numeric(default))
  ans <- readline(paste0(prompt, " [default ", default, "]: "))
  if (nchar(ans) == 0) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(ans))
  if (is.na(out) || out < 0) return(as.numeric(default))
  out
}

deps_snapshot <- function(pkgs) {
  pkgs <- unique(pkgs)
  out <- vector("list", length(pkgs))
  for (i in seq_along(pkgs)) {
    p <- pkgs[i]
    ver <- NA_character_
    if (requireNamespace(p, quietly = TRUE)) {
      ver <- as.character(utils::packageVersion(p))
    }
    out[[i]] <- list(package = p, version = ver)
  }
  out
}

sanitize_for_json <- function(x) {
  if (is.list(x)) {
    for (nm in names(x)) x[[nm]] <- sanitize_for_json(x[[nm]])
    return(x)
  }
  if (is.atomic(x)) {
    if (is.numeric(x)) x[!is.finite(x)] <- NA_real_
    if (is.character(x)) x <- enc2utf8(x)
    return(x)
  }
  x
}

write_json_atomic <- function(obj, path) {
  ensure_pkgs(c("jsonlite"))
  tmp <- paste0(path, ".tmp_", ts_stamp())
  jsonlite::write_json(obj, tmp, pretty = TRUE, auto_unbox = TRUE, na = "null")
  ok <- TRUE
  tryCatch({
    jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  }, error = function(e) {
    ok <<- FALSE
  })
  if (!ok) stop2("Manifest JSON failed sanity parse after write: ", tmp)
  file.rename(tmp, path)
  invisible(TRUE)
}

build_file_inventory <- function(output_dir, inventory_csv_path) {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, include.dirs = FALSE)
  files <- files[file.exists(files)]
  out_dir_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  files_norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
  
  esc <- gsub("([\\^\\$\\.|\\+\\*\\?\\(\\)\\[\\]\\{\\}\\\\\\|])", "\\\\\\1", out_dir_norm)
  rel <- sub(paste0("^", esc, "/?"), "", files_norm)
  
  info <- file.info(files_norm)
  inv <- data.frame(
    file = rel,
    full_path = files_norm,
    size_bytes = as.numeric(info$size),
    mtime = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
  inv <- inv[order(inv$file), , drop = FALSE]
  utils::write.csv(inv, inventory_csv_path, row.names = FALSE)
  inv
}
# ---------------------- Script identity capture (MANDATORY) ----------------------
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
known_script_filename <- "01_compute_NSV_trajectories.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  script_full_path <- NA_character_
  script_path_detection_note <- "Script path detection failed; using fallback filename."
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
  script_full_path <- script_full
  script_path_detection_note <- ""
}

msg("Script identity:")
msg("  script_name: ", script_name)
msg("  script_path: ", script_path)
msg("  script_full: ", ifelse(is.na(script_full_path), "NA", script_full_path))
if (nzchar(script_path_detection_note)) msg("NOTE: ", script_path_detection_note)

# ---------------------- Self-parse diagnostic (NEW) ----------------------
if (isTRUE(enable_self_parse_diagnostic) && !is.na(script_full_path) && file.exists(script_full_path)) {
  parse_ok <- TRUE
  parse_err <- NULL
  tryCatch({
    parse(file = script_full_path, keep.source = TRUE)
  }, error = function(e) {
    parse_ok <<- FALSE
    parse_err <<- conditionMessage(e)
  })
  if (!parse_ok) {
    stop2(
      "SELF-PARSE FAILED: the script file on disk has a syntax issue.\n",
      "Parser message: ", parse_err, "\n",
      "Fix the on-disk script and re-run."
    )
  }
}

# ------------------------------- Params ---------------------------------
in_path <- arg_val("--in", default = "WIDE_LiverMetabolome.csv")

time_col <- arg_val("--time_col", default = "Day")
id_col   <- arg_val("--id_col",   default = "SampleID")

transform_arg <- arg_val("--transform", default = NA_character_)
pseudocount_arg <- arg_val("--pseudocount", default = NA_character_)

transform_choice <- suppressWarnings(as.integer(transform_arg))
if (is.na(transform_choice) || !(transform_choice %in% c(0, 1))) transform_choice <- NA_integer_

pseudocount <- suppressWarnings(as.numeric(pseudocount_arg))
if (is.na(pseudocount) || pseudocount < 0) pseudocount <- NA_real_

min_sd <- 1e-8
top_var_features <- suppressWarnings(as.integer(arg_val("--top_var_features", default = "2000")))
n_subsample <- suppressWarnings(as.integer(arg_val("--n_subsample", default = "250")))
parallel_workers <- suppressWarnings(as.integer(arg_val("--parallel_workers", default = "8")))
latent_k_grid <- c(0, 2)
use_parallel <- TRUE
use_shrinkage_if_available <- TRUE

enable_plotly <- if (!is.na(arg_val("--enable_plotly", default = NA_character_))) is_true(arg_val("--enable_plotly")) else enable_plotly
enable_runtime_tracking <- if (!is.na(arg_val("--enable_runtime_tracking", default = NA_character_))) is_true(arg_val("--enable_runtime_tracking")) else enable_runtime_tracking

if (!file.exists(in_path)) {
  if (interactive()) {
    msg("Default input not found: ", in_path)
    in_path <- pick_file("Choose wide input CSV (must contain Day, SampleID, and feature columns)")
  } else {
    stop2("Input file not found: ", in_path, "\nProvide via: --in=/path/to/file.csv")
  }
}
in_path <- normalizePath(in_path, winslash = "/", mustWork = FALSE)

if (is.na(transform_choice)) {
  transform_choice <- ask_choice01("Choose transform: raw (0) or log2(TPM + pseudocount) (1)", default = 0)
}
if (transform_choice == 1 && is.na(pseudocount)) {
  pseudocount <- ask_num("Pseudocount for log2(TPM + pseudocount)", default = 1)
}
if (transform_choice == 0) pseudocount <- NA_real_

transform_label <- if (transform_choice == 0) "raw" else paste0("log2(TPM + ", pseudocount, ")")

# ---------------------- Output directory policy (MANDATORY) ----------------------
outputs_root <- file.path(getwd(), "outputs")
dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

analysis_name <- "NSV_TRAJECTORIES"
source_stem <- tools::file_path_sans_ext(basename(in_path))
run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_id <- paste0(analysis_name, "_", source_stem, "_", ts_stamp())

out_dir <- file.path(outputs_root, run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest_json_path <- file.path(out_dir, "Project_Manifest.json")
manifest_files_csv_path <- file.path(out_dir, "Project_Manifest_Files.csv")

# Initial inventory
build_file_inventory(out_dir, manifest_files_csv_path)

msg("Run summary:")
msg("  run_id:    ", run_id)
msg("  input:     ", in_path)
msg("  output:    ", out_dir)
msg("  transform: ", transform_label)
msg("  plotly:    ", enable_plotly, " | runtime_tracking: ", enable_runtime_tracking)

# Ensure warnings get captured to a file inside output_dir (NEW)
warnings_path <- file.path(out_dir, "WARNINGS.txt")
if (isTRUE(write_warnings_file)) {
  options(warn = 1)
}
# ----------------------------- Dependencies -----------------------------
base_deps <- c("jsonlite", "knitr", "rmarkdown", "ggplot2")
if (isTRUE(enable_plotly)) base_deps <- c(base_deps, "plotly", "htmltools")

ensure_pkgs(base_deps)
if (use_parallel) ensure_pkgs(c("future", "future.apply"))

if (use_shrinkage_if_available && !requireNamespace("corpcor", quietly = TRUE)) {
  msg("NOTE: corpcor not installed; using standard cor(). (Optional shrinkage)")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

# -------------------------- Scientific functions --------------------------
mad_robust <- function(x) {
  stats::mad(x, center = stats::median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
}

safe_scale <- function(X) {
  if (is.null(dim(X)) || ncol(X) < 2 || nrow(X) < 2) return(NULL)
  s <- apply(X, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(s) & (s > min_sd)
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 2 || nrow(X) < 2) return(NULL)
  scale(X, center = TRUE, scale = TRUE)
}

cor_robust <- function(Z) {
  if (is.null(Z) || ncol(Z) < 2 || nrow(Z) < 2) return(NULL)
  if (use_shrinkage_if_available && requireNamespace("corpcor", quietly = TRUE)) {
    corpcor::cor.shrink(Z, verbose = FALSE)
  } else {
    stats::cor(Z, use = "pairwise.complete.obs")
  }
}

effective_dimensionality <- function(R) {
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[is.finite(ev) & ev > 0]
  if (length(ev) < 2) return(NA_real_)
  (sum(ev)^2) / sum(ev^2)
}

frobenius_cosine_similarity <- function(A, B) {
  if (!identical(dim(A), dim(B))) return(NA_real_)
  num <- sum(A * B, na.rm = TRUE)
  den <- sqrt(sum(A^2, na.rm = TRUE) * sum(B^2, na.rm = TRUE))
  if (!is.finite(den) || den == 0) return(NA_real_)
  num / den
}

compute_D <- function(Xt) {
  per_gene <- apply(Xt, 2, mad_robust)
  stats::median(per_gene, na.rm = TRUE)
}

compute_O <- function(Xt, z_thresh = 3) {
  Z <- apply(Xt, 2, function(v) {
    m <- stats::median(v, na.rm = TRUE)
    s <- mad_robust(v)
    if (!is.finite(s) || s <= min_sd) return(rep(0, length(v)))
    (v - m) / s
  })
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
  extreme <- apply(abs(Z), 2, function(v) any(v > z_thresh, na.rm = TRUE))
  mean(extreme, na.rm = TRUE)
}

select_top_var <- function(X, top_n) {
  if (is.null(dim(X)) || ncol(X) <= top_n) return(X)
  v <- apply(X, 2, stats::var, na.rm = TRUE)
  keep <- order(v, decreasing = TRUE)[seq_len(top_n)]
  X[, keep, drop = FALSE]
}

# ✅ FIXED: Single correct definition
compute_nsv_table <- function(X, day) {
  days <- sort(unique(day))  # <-- FIXED HERE
  nsv <- data.frame(
    Day = days,
    n = as.integer(NA),
    D = NA_real_,
    ED = NA_real_,
    O = NA_real_
  )
  
  R_by_day <- vector("list", length(days))
  names(R_by_day) <- as.character(days)
  
  for (i in seq_along(days)) {
    d <- days[i]
    idx <- which(day == d)
    Xt <- X[idx, , drop = FALSE]
    
    nsv$n[i] <- nrow(Xt)
    nsv$D[i]  <- compute_D(Xt)
    nsv$O[i]  <- compute_O(Xt)
    
    Xt2 <- select_top_var(Xt, top_var_features)
    Z <- safe_scale(Xt2)
    
    if (is.null(Z)) {
      nsv$ED[i] <- NA_real_
      R_by_day[[i]] <- NULL
    } else {
      R <- cor_robust(Z)
      if (is.null(R)) {
        nsv$ED[i] <- NA_real_
        R_by_day[[i]] <- NULL
      } else {
        nsv$ED[i] <- effective_dimensionality(R)
        R_by_day[[i]] <- R
      }
    }
  }
  
  nsv$delta_day <- c(NA_real_, diff(nsv$Day))
  nsv$S_adj <- NA_real_
  
  for (i in 2:nrow(nsv)) {
    Rprev <- R_by_day[[i - 1]]
    Rcur  <- R_by_day[[i]]
    
    if (is.null(Rprev) || is.null(Rcur)) {
      nsv$S_adj[i] <- NA_real_
      next
    }
    
    if (!is.null(rownames(Rprev)) && !is.null(rownames(Rcur))) {
      common <- intersect(rownames(Rprev), rownames(Rcur))
      if (length(common) >= 2) {
        Rprev <- Rprev[common, common, drop = FALSE]
        Rcur  <- Rcur[common, common, drop = FALSE]
      } else {
        nsv$S_adj[i] <- NA_real_
        next
      }
    }
    
    sim <- frobenius_cosine_similarity(Rprev, Rcur)
    nsv$S_adj[i] <- if (is.finite(sim)) 1 - sim else NA_real_
  }
  
  nsv$S_per_day <- nsv$S_adj / nsv$delta_day
  
  nsv$curvature <- NA_real_
  if (nrow(nsv) >= 3) {
    for (i in 2:(nrow(nsv) - 1)) {
      v1 <- as.numeric(nsv[i, c("D", "ED", "O")]) - as.numeric(nsv[i - 1, c("D", "ED", "O")])
      v2 <- as.numeric(nsv[i + 1, c("D", "ED", "O")]) - as.numeric(nsv[i, c("D", "ED", "O")])
      n1 <- sqrt(sum(v1^2, na.rm = TRUE))
      n2 <- sqrt(sum(v2^2, na.rm = TRUE))
      if (is.finite(n1) && is.finite(n2) && n1 > 0 && n2 > 0) {
        cosang <- sum(v1 * v2, na.rm = TRUE) / (n1 * n2)
        cosang <- max(-1, min(1, cosang))
        nsv$curvature[i] <- 1 - cosang
      }
    }
  }
  
  nsv
}

residualize_day <- function(X, day) {
  day_f <- as.factor(day)
  M <- model.matrix(~ day_f)
  Xt <- as.matrix(X)
  beta <- solve(t(M) %*% M, t(M) %*% Xt)
  resid <- Xt - M %*% beta
  resid
}

remove_top_pcs <- function(X, k) {
  if (k <= 0 || ncol(X) < 2 || nrow(X) < 2) return(X)
  Z <- scale(X, center = TRUE, scale = FALSE)
  keep <- apply(Z, 2, function(v) sd(v, na.rm = TRUE) > min_sd)
  Z <- Z[, keep, drop = FALSE]
  if (ncol(Z) < 2) return(X)
  kept_cols <- colnames(Z)
  k_eff <- min(k, nrow(Z) - 1, ncol(Z) - 1)
  sv <- svd(Z, nu = k_eff, nv = 0)
  scores <- sv$u %*% diag(sv$d[1:k_eff], nrow = k_eff)
  XtX <- crossprod(scores)
  if (qr(XtX)$rank < ncol(scores)) return(X)
  B <- solve(XtX, crossprod(scores, Z))
  nuisance <- scores %*% B
  Z2 <- Z - nuisance
  X_out <- X
  X_out[, kept_cols] <- Z2
  X_out
}
# -------------------------- Load data --------------------------
msg("Loading: ", in_path)
df <- read.csv(in_path, stringsAsFactors = FALSE, check.names = FALSE)

# Case-insensitive column resolution (prevents Day vs day failures)
cn <- colnames(df)
if (!(time_col %in% cn)) {
  hit <- which(tolower(cn) == tolower(time_col))
  if (length(hit) == 1) time_col <- cn[hit]
}
if (!(id_col %in% cn)) {
  hit <- which(tolower(cn) == tolower(id_col))
  if (length(hit) == 1) id_col <- cn[hit]
}

stopifnot(time_col %in% colnames(df), id_col %in% colnames(df))

day <- df[[time_col]]
sample_id <- df[[id_col]]

# FIX: canonical day vector used throughout (prevents "object 'days' not found")
days <- sort(unique(day))

X <- as.matrix(df[, !(colnames(df) %in% c(time_col, id_col)), drop = FALSE])
storage.mode(X) <- "double"

# -------------------------- Transform (NEW; user choice) --------------------------
if (transform_choice == 1) {
  msg("Applying ", transform_label, " transform")
  if (any(X < 0, na.rm = TRUE)) stop2("Found negative values; cannot apply log2(TPM + pseudocount).")
  X <- log2(X + pseudocount)
} else {
  msg("Using raw values (no transform)")
}

# -------------------------- QC summary --------------------------
qc_lines <- character(0)
qc_lines <- c(qc_lines, paste("QC SUMMARY:", Sys.time()))
qc_lines <- c(qc_lines, paste("Input:", in_path))
qc_lines <- c(qc_lines, paste("Transform:", transform_label))
qc_lines <- c(qc_lines, paste("Samples:", nrow(X), "Features:", ncol(X)))

tab_day <- table(day)
qc_lines <- c(qc_lines, "Per-day n:")
qc_lines <- c(qc_lines, paste(capture.output(print(tab_day)), collapse = "\n"))

days_sorted <- sort(unique(day))
days <- days_sorted   # FIX: canonical days vector used throughout the script
deltas <- diff(days_sorted)

qc_lines <- c(qc_lines, paste("Observed days:", paste(days_sorted, collapse = ", ")))
qc_lines <- c(qc_lines, paste("Inter-day gaps:", paste(deltas, collapse = ", ")))

if (length(unique(deltas)) > 1) {
  qc_lines <- c(qc_lines, "WARNING: Unequal time gaps detected. Interpret S_adj as total change over the interval.")
  qc_lines <- c(qc_lines, "         Use S_per_day for approximate rate normalization.")
}
if (length(deltas) > 0 && any(deltas > min(deltas))) {
  qc_lines <- c(qc_lines, paste("NOTE: Largest gap =", max(deltas), "days; possible missing intermediate timepoint(s)."))
}

writeLines(qc_lines, file.path(out_dir, "QC_SUMMARY.txt"))
msg("Wrote QC summary: QC_SUMMARY.txt")

# -------------------------- Basic QC --------------------------
msg("Per-day n:")
print(tab_day)

sd_all <- apply(X, 2, stats::sd, na.rm = TRUE)
msg("Range of feature SDs before filtering: ", paste0("[", min(sd_all, na.rm = TRUE), ", ", max(sd_all, na.rm = TRUE), "]"))

keep <- is.finite(sd_all) & (sd_all >= min_sd)
if (sum(keep) < 2) {
  msg("WARNING: Too few variable features (", sum(keep), "); keeping top 2 by SD anyway")
  keep <- order(sd_all, decreasing = TRUE)[1:2]
}
X <- X[, keep, drop = FALSE]
msg("Features after removing low-variance features: ", ncol(X))


n_min <- min(as.integer(tab_day))
msg("n_min for equal-n subsampling: ", n_min)

# Diagnostics (NEW): feature counts per day + mode will be recorded later
debug_day_feature_counts <- data.frame()

# -------------------------- Parallel setup (future.apply) --------------------------
if (use_parallel) {
  suppressPackageStartupMessages({
    library(future)
    library(future.apply)
  })
  future::plan(future::multisession, workers = parallel_workers)
  msg("Parallel subsampling enabled with ", parallel_workers, " workers")
} else {
  msg("Parallel subsampling disabled")
}

# -------------------------- Compute NSV: modes --------------------------
modes <- list()
modes[["RAW"]] <- list(X = X)

X_resid <- residualize_day(X, day)
modes[["RESID_DAY"]] <- list(X = X_resid)

if (is.null(dim(X_resid)) || nrow(X_resid) < 2 || ncol(X_resid) < 2) {
  msg("WARNING: RESID_DAY matrix is <2x2; skipping RESID_DAY_PCk modes.")
} else {
  for (k in latent_k_grid[latent_k_grid > 0]) {
    Xk <- remove_top_pcs(X_resid, k = k)
    modes[[paste0("RESID_DAY_PC", k)]] <- list(X = Xk)
  }
}

all_results <- list()
all_subsample <- list()

for (nm in names(modes)) {
  msg("Mode: ", nm)
  Xm <- modes[[nm]]$X
  
  # Diagnostics (NEW): per-day dimension snapshot
  if (isTRUE(debug_matrix_sizes)) {
    days_u <- sort(unique(day))
    for (d in days_u) {
      idx <- which(day == d)
      Xt <- Xm[idx, , drop = FALSE]
      debug_day_feature_counts <- rbind(
        debug_day_feature_counts,
        data.frame(
          Mode = nm,
          Day = d,
          n_samples = nrow(Xt),
          n_features = ncol(Xt),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  nsv <- compute_nsv_table(Xm, day)
  all_results[[nm]] <- nsv
  write.csv(nsv, file.path(out_dir, paste0("NSV_", nm, "_full.csv")), row.names = FALSE)
  
  if (use_parallel) {
    sub_tbls <- future.apply::future_lapply(seq_len(n_subsample), function(b) {
      sub <- subsample_equal_n(Xm, day, n_min)
      compute_nsv_table(sub$X, sub$day)
    })
  } else {
    sub_tbls <- vector("list", n_subsample)
    for (b in seq_len(n_subsample)) {
      sub <- subsample_equal_n(Xm, day, n_min)
      sub_tbls[[b]] <- compute_nsv_table(sub$X, sub$day)
    }
  }
  
  days_u <- sort(unique(day))
  summ <- data.frame(Day = days_u)
  
  for (col in c("D", "ED", "O", "S_adj", "S_per_day", "curvature")) {
    mat <- sapply(sub_tbls, function(tb) tb[[col]])
    summ[[paste0(col, "_med")]] <- apply(mat, 1, stats::median, na.rm = TRUE)
    summ[[paste0(col, "_lo")]]  <- apply(mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
    summ[[paste0(col, "_hi")]]  <- apply(mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  }
  
  all_subsample[[nm]] <- summ
  write.csv(summ, file.path(out_dir, paste0("NSV_", nm, "_subsample_CI.csv")), row.names = FALSE)
}

if (isTRUE(debug_matrix_sizes) && nrow(debug_day_feature_counts) > 0) {
  write.csv(debug_day_feature_counts, file.path(out_dir, "DEBUG_day_sample_feature_counts_by_mode.csv"), row.names = FALSE)
  msg("Wrote diagnostics: DEBUG_day_sample_feature_counts_by_mode.csv")
}
# -------------------------- Plots (write into output_dir) --------------------------
png(file.path(out_dir, "NSV_panels.png"), width = 1600, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot_panel <- function(mode_name, y_col, ylab) {
  s <- all_subsample[[mode_name]]
  x <- s$Day
  y <- s[[paste0(y_col, "_med")]]
  lo <- s[[paste0(y_col, "_lo")]]
  hi <- s[[paste0(y_col, "_hi")]]
  
  plot(x, y, type = "b", pch = 16, xlab = "Day", ylab = ylab,
       main = paste(mode_name, "-", y_col))
  arrows(x0 = x, y0 = lo, x1 = x, y1 = hi, angle = 90, code = 3, length = 0.03)
}

mode_to_plot <- "RAW"
plot_panel(mode_to_plot, "D", paste0("Dispersion D(t) (median gene MAD) [", transform_label, "]"))
plot_panel(mode_to_plot, "ED", "Effective dimensionality ED(t)")
plot_panel(mode_to_plot, "S_adj", "Structural change S_adj (1 - similarity)")
plot_panel(mode_to_plot, "S_per_day", "Structural change rate S_per_day (per day)")
dev.off()

png(file.path(out_dir, "NSV_phase_D_vs_ED.png"), width = 1000, height = 800)
nsv_raw <- all_results[["RAW"]]
plot(nsv_raw$D, nsv_raw$ED, type = "b", pch = 16,
     xlab = "D(t) dispersion", ylab = "ED(t) effective dimensionality",
     main = paste0("Phase portrait: D vs ED (RAW; ", transform_label, ")"))
text(nsv_raw$D, nsv_raw$ED, labels = nsv_raw$Day, pos = 3, cex = 0.9)
dev.off()

# Save ggplot objects for report
df_pan <- all_subsample[["RAW"]]
df_pan_long <- data.frame(
  Day = rep(df_pan$Day, 4),
  Metric = rep(c("D", "ED", "S_adj", "S_per_day"), each = nrow(df_pan)),
  med = c(df_pan$D_med, df_pan$ED_med, df_pan$S_adj_med, df_pan$S_per_day_med),
  lo  = c(df_pan$D_lo,  df_pan$ED_lo,  df_pan$S_adj_lo,  df_pan$S_per_day_lo),
  hi  = c(df_pan$D_hi,  df_pan$ED_hi,  df_pan$S_adj_hi,  df_pan$S_per_day_hi),
  stringsAsFactors = FALSE
)

pA <- ggplot(df_pan_long, aes(x = Day, y = med)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.0) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "NSV panels (RAW; subsampling median ± 95% CI)",
    subtitle = paste0("Transform: ", transform_label),
    x = "Day",
    y = "Value"
  )

pB <- ggplot(nsv_raw, aes(x = D, y = ED, label = Day)) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_text(vjust = -0.7, size = 3.6) +
  theme_classic(base_size = 12) +
  labs(
    title = "Phase portrait: D vs ED (RAW)",
    subtitle = paste0("Transform: ", transform_label),
    x = "Dispersion D(t)",
    y = "Effective dimensionality ED(t)"
  )

saveRDS(pA, file.path(out_dir, paste0(script_name, "_pA.rds")))
saveRDS(pB, file.path(out_dir, paste0(script_name, "_pB.rds")))

# -------------------------- README --------------------------
readme <- c(
  "NSV Starter Run",
  paste("Run time:", Sys.time()),
  paste("Input:", in_path),
  paste("Transform:", transform_label),
  "",
  "QC:",
  "- QC_SUMMARY.txt includes time gaps; unequal gaps trigger warnings",
  "",
  "Outputs:",
  "- NSV_<MODE>_full.csv: per-day NSV values on full data",
  "- NSV_<MODE>_subsample_CI.csv: subsampling-based median and 95% CI per day",
  "- NSV_panels.png: D, ED, S_adj, S_per_day panels for RAW",
  "- NSV_phase_D_vs_ED.png: phase portrait for RAW",
  "- <script_name>_pA.rds / <script_name>_pB.rds: ggplot objects for report/plotly",
  "",
  "Modes:",
  "RAW: no correction",
  "RESID_DAY: regress out day means (removes mean trajectory)",
  "RESID_DAY_PCk: remove top k latent PCs after day-mean residualization (proxy for unknown batch)",
  "",
  "Interval-aware structure:",
  "- S_adj is total reconfiguration over the observed interval",
  "- S_per_day = S_adj / delta_day (approximate rate normalization)"
)
writeLines(readme, file.path(out_dir, "README.txt"))

# -------------------------- Manifest Finalize --------------------------
runtime_seconds <- if (isTRUE(enable_runtime_tracking)) round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2) else NA_real_

manifest_obj <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(name = script_name, path = script_path, full_path = script_full_path),
  input = list(wide_csv = in_path),
  parameters = list(
    time_col = time_col,
    id_col = id_col,
    transform_choice = transform_choice,
    pseudocount = if (transform_choice == 1) pseudocount else NA_real_,
    transform_label = transform_label,
    top_var_features = top_var_features,
    n_subsample = n_subsample,
    latent_k_grid = latent_k_grid,
    parallel_workers = parallel_workers,
    enable_plotly = enable_plotly,
    enable_runtime_tracking = enable_runtime_tracking,
    runtime_seconds = runtime_seconds
  ),
  dependencies = deps_snapshot(base_deps),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  )
)
manifest_obj <- sanitize_for_json(manifest_obj)
write_json_atomic(manifest_obj, manifest_json_path)

# Final file inventory
build_file_inventory(out_dir, manifest_files_csv_path)

# -------------------------- Finish --------------------------
msg("Done. Files written to: ", out_dir)
msg("Manifest:")
msg(" - ", manifest_json_path)
msg(" - ", manifest_files_csv_path)

# -------------------------- Quarto QMD + Render --------------------------

report_qmd_path  <- file.path(out_dir, paste0(script_name, "_Report.qmd"))
report_html_path <- file.path(out_dir, paste0(script_name, "_Report.html"))

# Write QMD file (qmd_lines is assumed to be defined earlier — as in your original script)
writeLines(qmd_lines, con = report_qmd_path)

# Render the QMD file to HTML
render_ok <- TRUE
render_err <- NULL

tryCatch({
  suppressPackageStartupMessages(library(rmarkdown))
  rmarkdown::render(
    input = report_qmd_path,
    output_file = basename(report_html_path),
    output_dir = out_dir,
    quiet = TRUE,
    knit_root_dir = out_dir,
    params = list(
      output_dir = ".",
      runtime_seconds = runtime_seconds,
      enable_plotly = isTRUE(enable_plotly)
    )
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- conditionMessage(e)
})

# Capture warnings to file
if (isTRUE(write_warnings_file)) {
  w <- warnings()
  if (!is.null(w)) {
    writeLines(c("WARNINGS CAPTURED:", capture.output(print(w))), warnings_path)
    msg("Wrote warnings file: ", warnings_path)
  }
}

# Refresh inventory AFTER render
build_file_inventory(out_dir, manifest_files_csv_path)

# Update manifest with report info
manifest_obj$notes$report_qmd_path <- normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE)
manifest_obj$notes$report_html_path <- normalizePath(report_html_path, winslash = "/", mustWork = FALSE)
manifest_obj$notes$report_render_success <- render_ok
manifest_obj$notes$report_render_error <- if (!render_ok) render_err else NA_character_
manifest_obj <- sanitize_for_json(manifest_obj)
write_json_atomic(manifest_obj, manifest_json_path)

# Final message about report status
if (render_ok && file.exists(report_html_path)) {
  msg("HTML report created: ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE))
} else {
  msg("REPORT RENDER FAILED. QMD written at: ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE))
  if (!render_ok && !is.null(render_err)) msg("Render error: ", render_err)
}

