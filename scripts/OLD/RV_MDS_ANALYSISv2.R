#!/usr/bin/env Rscript
# ======================================================================
# RV_MDS_ANALYSIS.R
# Day–Day Correlation-Structure Similarity: NON-HEATMAP Publication Plots
# + Integrated "Optimized" RV Bootstrap (row-resampling; memory-conscious)
#
# PRIMARY INPUT (CSV; long pairwise): Day1, Day2, matrix_corr, frobenius, RV
# OPTIONAL INPUT (for RV bootstrap only; tidy expression): Day, SampleID, Gene, Abundance
#
# OUTPUT (standardized under outputs/<run_id>/):
#   - similarity_long_clean.csv
#   - similarity_wide_matrix_corr.csv
#   - similarity_wide_RV.csv
#   - similarity_wide_frobenius.csv
#   - RefDay_Profile_<metric>_data.csv
#   - Plot_RefDay_Profile_<metric>.png/.pdf
#   - Plot_MDS_<metric>.png/.pdf
#   - Plot_DayNetwork_<metric>.png/.pdf
#   - mds_coords_<metric>.csv
#   - network_edges_<metric>_topk.csv
#   - RefDay_Profile_RV_bootstrap_CI.csv              (if metric=RV and bootstrap runs)
#   - RV_bootstrap_draws_profile_wide.csv             (optional; if --save_boot_draws=TRUE)
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - <script_name>_Report.qmd + rendered HTML
#
# USAGE (pairwise only):
#   Rscript RV_MDS_ANALYSIS.R --in=A_GlobalSimilarity_ByDayPairs.csv --metric=RV --ref_day=14 --top_k=4
#
# USAGE (with bootstrap):
#   Rscript RV_MDS_ANALYSIS.R --in=A_GlobalSimilarity_ByDayPairs.csv --metric=RV --ref_day=14 --top_k=4 \
#          --expr=Transcriptome_Tidy_Long.csv --nboot=200 --seed=1 --ci=0.95 --transform=log1p --max_genes=0
#
# NOTES ON BOOTSTRAP MODE:
#   This integrates an "optimized" bootstrap that resamples ROWS of each day's
#   Sample×Gene matrix (after preprocessing). This is fast and memory-conscious.
#   It is a pragmatic robustness check, but it is NOT identical to resampling
#   biological replicates at the SampleID level if samples differ in gene missingness.
#
# IMPORTANT:
#   All if/else logic is written in BLOCK FORM to avoid R's "unexpected else" errors.
# ======================================================================

# ------------------------------- Helpers --------------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2    <- function(...) stop(paste0(...), call. = FALSE)
msg      <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

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

pick_file <- function(prompt = "Choose file") {
  if (interactive() && exists("file.choose")) {
    msg(prompt)
    return(file.choose())
  }
  stop2(prompt, " (non-interactive). Provide a path via command line args.")
}

ask_int <- function(prompt, default) {
  if (!interactive()) return(as.integer(default))
  ans <- readline(paste0(prompt, " [default ", default, "]: "))
  if (nchar(ans) == 0) return(as.integer(default))
  out <- suppressWarnings(as.integer(ans))
  if (is.na(out) || out < 1) return(as.integer(default))
  out
}

ask_yesno <- function(prompt, default_yes = TRUE) {
  if (!interactive()) return(default_yes)
  def <- if (default_yes) "Y/n" else "y/N"
  ans <- readline(paste0(prompt, " [", def, "]: "))
  if (nchar(ans) == 0) return(default_yes)
  ans <- tolower(ans)
  if (ans %in% c("y", "yes")) return(TRUE)
  if (ans %in% c("n", "no")) return(FALSE)
  default_yes
}

is_true <- function(x) {
  if (is.null(x) || is.na(x)) return(FALSE)
  x <- tolower(as.character(x))
  x %in% c("1", "true", "t", "yes", "y")
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
known_script_filename <- "RV_MDS_ANALYSIS.R"
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

# ---------------------- Manifest + inventory utilities (MANDATORY) ----------------------

build_file_inventory <- function(output_dir, inventory_csv_path) {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, include.dirs = FALSE)
  files <- files[file.exists(files)]
  out_dir_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  files_norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
  
  # relative path
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

# Sanitize a nested list so jsonlite always writes valid JSON
sanitize_for_json <- function(x) {
  if (is.list(x)) {
    for (nm in names(x)) x[[nm]] <- sanitize_for_json(x[[nm]])
    return(x)
  }
  if (is.atomic(x)) {
    # convert non-finite numerics to NA (prevents NaN/Inf tokens)
    if (is.numeric(x)) {
      x[!is.finite(x)] <- NA_real_
    }
    # ensure UTF-8 for characters (prevents weird encoding artifacts)
    if (is.character(x)) {
      x <- enc2utf8(x)
    }
    return(x)
  }
  x
}

write_json_atomic <- function(obj, path) {
  ensure_pkgs(c("jsonlite"))
  tmp <- paste0(path, ".tmp_", ts_stamp())
  jsonlite::write_json(obj, tmp, pretty = TRUE, auto_unbox = TRUE, na = "null")
  # basic sanity check that it parses
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

# ------------------------------- Params ---------------------------------

in_path <- arg_val("--in", default = NA_character_)
metric  <- arg_val("--metric", default = "RV")     # RV | matrix_corr | frobenius
ref_day <- suppressWarnings(as.numeric(arg_val("--ref_day", default = "14")))
top_k   <- suppressWarnings(as.integer(arg_val("--top_k", default = "4")))

# Bootstrap params (RV only)
expr_path         <- arg_val("--expr", default = NA_character_)
nboot_arg         <- arg_val("--nboot", default = NA_character_)
seed              <- suppressWarnings(as.integer(arg_val("--seed", default = "1")))
ci_level          <- suppressWarnings(as.numeric(arg_val("--ci", default = "0.95")))
transform         <- arg_val("--transform", default = "log1p")  # none | log1p
max_genes         <- suppressWarnings(as.integer(arg_val("--max_genes", default = "0"))) # 0 = no cap
run_boot_arg      <- arg_val("--run_boot", default = NA_character_) # TRUE/FALSE; if absent, prompt
save_boot_draws   <- is_true(arg_val("--save_boot_draws", default = "FALSE"))

# Debug toggles
debug            <- is_true(arg_val("--debug", default = "TRUE"))
save_debug_files <- is_true(arg_val("--save_debug_files", default = "TRUE"))

msg("NOTE: This script expects a file named 'A_GlobalSimilarity_ByDayPairs.csv'.")
msg("      This file is typically generated by 'Top50_TxCorrStructure_ByDay_Liver.R'.")
Sys.sleep(2)
# Resolve inputs
if (is.na(in_path) || nchar(in_path) == 0) {
  in_path <- pick_file("Choose pairwise Day1–Day2 similarity CSV (Day1, Day2, matrix_corr, frobenius, RV)")
}
if (!file.exists(in_path)) stop2("Input file not found: ", in_path)


if (tolower(metric) == "rv") metric <- "RV"
if (!metric %in% c("matrix_corr", "RV", "frobenius")) stop2("--metric must be one of: RV, matrix_corr, frobenius")
if (ci_level <= 0 || ci_level >= 1) stop2("--ci must be between 0 and 1 (e.g., 0.95).")
alpha <- (1 - ci_level) / 2

# ---------------------- Output directory policy (MANDATORY) ----------------------

outputs_root <- file.path(getwd(), "outputs")
dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

analysis_name <- "RV_MDS_ANALYSIS"
source_stem <- tools::file_path_sans_ext(basename(in_path))
run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_id <- paste0(analysis_name, "_", source_stem, "_", ts_stamp())

out_dir <- file.path(outputs_root, run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest_json_path <- file.path(out_dir, "Project_Manifest.json")
manifest_files_csv_path <- file.path(out_dir, "Project_Manifest_Files.csv")

# initial inventory
build_file_inventory(out_dir, manifest_files_csv_path)

msg("Run summary:")
msg("  run_id:    ", run_id)
msg("  input:     ", in_path)
msg("  metric:    ", metric)
msg("  ref_day:   ", ref_day)
msg("  top_k:     ", top_k)
msg("  output:    ", out_dir)
msg("  debug:     ", debug, " | save_debug_files=", save_debug_files)

# ----------------------------- Dependencies -----------------------------

ensure_pkgs(c("ggplot2", "igraph", "jsonlite", "knitr", "rmarkdown"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(igraph)
})

theme_pub <- theme_classic(base_size = 13) +
  theme(
    axis.title    = element_text(face = "bold"),
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.title  = element_text(face = "bold"),
    panel.grid    = element_blank()
  )

# ------------------------------ Read + QC -------------------------------

df <- read.csv(in_path, stringsAsFactors = FALSE, check.names = FALSE)

req <- c("Day1", "Day2", "matrix_corr", "frobenius", "RV")
missing_cols <- setdiff(req, colnames(df))
if (length(missing_cols) > 0) stop2("Missing required columns: ", paste(missing_cols, collapse = ", "))

df$Day1 <- suppressWarnings(as.numeric(df$Day1))
df$Day2 <- suppressWarnings(as.numeric(df$Day2))
if (any(is.na(df$Day1)) || any(is.na(df$Day2))) stop2("Day1/Day2 contain non-numeric values after coercion.")

for (m in c("matrix_corr", "frobenius", "RV")) {
  df[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  if (any(is.na(df[[m]]))) stop2("Metric '", m, "' has non-numeric values or NA.")
}

# Deduplicate by averaging
key <- paste(df$Day1, df$Day2, sep = "_")
if (any(duplicated(key))) {
  msg("Duplicate Day1-Day2 pairs found; averaging duplicates.")
  df <- aggregate(df[, c("matrix_corr", "frobenius", "RV")],
                  by = list(Day1 = df$Day1, Day2 = df$Day2),
                  FUN = mean)
}

# Symmetrize
df_rev <- df
df_rev$Day1 <- df$Day2
df_rev$Day2 <- df$Day1
df_sym <- rbind(df, df_rev)

days <- sort(unique(c(df_sym$Day1, df_sym$Day2)))

# Add diagonals
diag_df <- data.frame(Day1 = days, Day2 = days, matrix_corr = 1, RV = 1, frobenius = 0)
df_sym <- rbind(df_sym, diag_df)

write.csv(df_sym, file.path(out_dir, "similarity_long_clean.csv"), row.names = FALSE)

# ------------------------------ Wide matrices ---------------------------

to_wide <- function(d, value_col, diag_value, days_vec) {
  mat <- matrix(NA_real_,
                nrow = length(days_vec),
                ncol = length(days_vec),
                dimnames = list(as.character(days_vec), as.character(days_vec)))
  for (i in seq_len(nrow(d))) {
    r <- as.character(d$Day1[i])
    c <- as.character(d$Day2[i])
    mat[r, c] <- d[[value_col]][i]
  }
  diag(mat) <- diag_value
  mat
}

mat_mc <- to_wide(df_sym, "matrix_corr", 1, days)
mat_rv <- to_wide(df_sym, "RV",         1, days)
mat_fr <- to_wide(df_sym, "frobenius",  0, days)

write.csv(mat_mc, file.path(out_dir, "similarity_wide_matrix_corr.csv"))
write.csv(mat_rv, file.path(out_dir, "similarity_wide_RV.csv"))
write.csv(mat_fr, file.path(out_dir, "similarity_wide_frobenius.csv"))

if (metric == "matrix_corr") {
  M <- mat_mc
} else if (metric == "RV") {
  M <- mat_rv
} else {
  M <- mat_fr
}

# ------------------------------ Plot 1: Reference-day profile -----------

ref_profile_path <- file.path(out_dir, paste0("RefDay_Profile_", metric, "_data.csv"))

if (!is.na(ref_day) && (ref_day %in% days)) {
  
  vals <- as.numeric(M[as.character(days), as.character(ref_day)])
  df_ref <- data.frame(Day = days, Value = vals, Metric = metric, RefDay = ref_day, stringsAsFactors = FALSE)
  write.csv(df_ref, ref_profile_path, row.names = FALSE)
  
  if (metric == "frobenius") {
    ylab <- "Frobenius distance (lower = more similar)"
  } else {
    ylab <- paste0(metric, " similarity to Day ", ref_day)
  }
  
  p_ref <- ggplot(df_ref, aes(x = Day, y = Value)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_x_continuous(breaks = days) +
    labs(
      title = paste0("Global transcriptome structure relative to Day ", ref_day),
      subtitle = paste0("Metric: ", metric),
      x = "Day",
      y = ylab
    ) +
    theme_pub
  
  ggsave(file.path(out_dir, paste0("Plot_RefDay_Profile_", metric, ".png")), p_ref,
         width = 8.2, height = 4.6, dpi = 300)
  ggsave(file.path(out_dir, paste0("Plot_RefDay_Profile_", metric, ".pdf")), p_ref,
         width = 8.2, height = 4.6)
  
} else {
  msg("Reference day ", ref_day, " not present; skipping reference-day profile plot.")
}

# ------------------------------ Plot 2: MDS embedding -------------------

if (metric == "frobenius") {
  D <- M
} else {
  S <- pmin(pmax(M, 0), 1)
  D <- 1 - S
}

if (any(is.na(D))) stop2("Distance matrix contains NA; missing day pairs in input table.")

mds <- cmdscale(as.dist(D), k = 2, eig = TRUE)
coords <- data.frame(Day = days, Dim1 = mds$points[, 1], Dim2 = mds$points[, 2], stringsAsFactors = FALSE)
write.csv(coords, file.path(out_dir, paste0("mds_coords_", metric, ".csv")), row.names = FALSE)

p_mds <- ggplot(coords, aes(x = Dim1, y = Dim2, label = Day)) +
  geom_point(size = 2.6) +
  geom_text(vjust = -0.7, size = 4) +
  labs(
    title = paste0("2D embedding of days from global structure (", metric, ")"),
    x = "MDS dimension 1",
    y = "MDS dimension 2"
  ) +
  theme_pub

ggsave(file.path(out_dir, paste0("Plot_MDS_", metric, ".png")), p_mds,
       width = 6.6, height = 5.4, dpi = 300)
ggsave(file.path(out_dir, paste0("Plot_MDS_", metric, ".pdf")), p_mds,
       width = 6.6, height = 5.4)

# ------------------------------ Plot 3: Day network ---------------------

edges <- data.frame(from = numeric(0), to = numeric(0), weight = numeric(0), stringsAsFactors = FALSE)

for (d in days) {
  vals <- M[as.character(d), as.character(days)]
  names(vals) <- as.character(days)
  vals <- vals[names(vals) != as.character(d)]
  
  if (metric == "frobenius") {
    ord <- order(vals, decreasing = FALSE)  # smaller = closer
  } else {
    ord <- order(vals, decreasing = TRUE)   # larger = closer
  }
  
  keep <- head(ord, top_k)
  to_days <- as.numeric(names(vals)[keep])
  
  sub <- data.frame(from = d, to = to_days, weight = as.numeric(vals[keep]), stringsAsFactors = FALSE)
  edges <- rbind(edges, sub)
}

if (metric == "frobenius") {
  edges <- edges[order(edges$weight, decreasing = FALSE), ]
} else {
  edges <- edges[order(edges$weight, decreasing = TRUE), ]
}

edges$pair <- apply(edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "_"))
edges <- edges[!duplicated(edges$pair), ]
edges$pair <- NULL

if (nrow(edges) < 1) stop2("No edges were constructed (top_k too small or matrix invalid).")

write.csv(edges, file.path(out_dir, paste0("network_edges_", metric, "_topk.csv")), row.names = FALSE)

g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = as.character(days)))

set.seed(1)
lay <- layout_with_fr(g)

png(file.path(out_dir, paste0("Plot_DayNetwork_", metric, ".png")), width = 1400, height = 1000, res = 150)
plot(g, layout = lay,
     vertex.size = 30,
     vertex.label.cex = 1.1,
     edge.width = 1.2,
     main = paste0("Day network from global structure (", metric, "), top_k=", top_k))
dev.off()

pdf(file.path(out_dir, paste0("Plot_DayNetwork_", metric, ".pdf")), width = 9, height = 6.5)
plot(g, layout = lay,
     vertex.size = 22,
     vertex.label.cex = 1.0,
     edge.width = 1.0,
     main = paste0("Day network from global structure (", metric, "), top_k=", top_k))
dev.off()

# ======================================================================
# RV BOOTSTRAP (integrated optimized row-resampling approach)
# ======================================================================

build_day_matrix <- function(expr_df, day, trans = "log1p", max_g = 0) {
  d <- expr_df[expr_df$Day == day, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  
  samps <- unique(d$SampleID)
  genes <- unique(d$Gene)
  
  if (max_g > 0 && length(genes) > max_g) {
    genes <- genes[seq_len(max_g)]
    d <- d[d$Gene %in% genes, , drop = FALSE]
  }
  
  i <- match(d$SampleID, samps)
  j <- match(d$Gene, genes)
  ok <- !(is.na(i) | is.na(j) | is.na(d$Abundance))
  d  <- d[ok, , drop = FALSE]
  i  <- i[ok]
  j  <- j[ok]
  
  mat <- matrix(NA_real_, nrow = length(samps), ncol = length(genes),
                dimnames = list(samps, genes))
  mat[cbind(i, j)] <- d$Abundance
  
  for (col in seq_len(ncol(mat))) {
    v <- mat[, col]
    if (all(is.na(v))) next
    med <- median(v, na.rm = TRUE)
    v[is.na(v)] <- med
    if (trans == "log1p") {
      v <- log1p(v)
    } else if (trans == "none") {
      # no-op
    } else {
      stop2("Unknown --transform: ", trans, " (use none|log1p)")
    }
    mat[, col] <- v
  }
  
  mat <- scale(mat, center = TRUE, scale = TRUE)
  
  keep <- apply(mat, 2, function(x) {
    if (!all(is.finite(x))) return(FALSE)
    if (stats::sd(x) <= 0) return(FALSE)
    TRUE
  })
  
  mat <- mat[, keep, drop = FALSE]
  if (is.null(dim(mat)) || ncol(mat) < 2 || nrow(mat) < 2) return(NULL)
  mat
}

calc_rv <- function(matA, matB) {
  if (is.null(matA) || is.null(matB)) return(NA_real_)
  g <- intersect(colnames(matA), colnames(matB))
  if (length(g) < 2) return(NA_real_)
  g <- sort(g)
  
  matA <- matA[, g, drop = FALSE]
  matB <- matB[, g, drop = FALSE]
  
  c1 <- stats::cor(matA, use = "pairwise.complete.obs")
  c2 <- stats::cor(matB, use = "pairwise.complete.obs")
  
  c1[is.na(c1)] <- 0
  c2[is.na(c2)] <- 0
  
  if (!identical(dim(c1), dim(c2))) return(NA_real_)
  
  num <- sum(c1 * c2)
  den <- sqrt(sum(c1^2) * sum(c2^2))
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

maybe_run_bootstrap <- FALSE
if (metric == "RV") {
  if (!is.na(run_boot_arg) && nchar(run_boot_arg) > 0) {
    maybe_run_bootstrap <- is_true(run_boot_arg)
  } else {
    maybe_run_bootstrap <- ask_yesno("Run optimized RV bootstrap (row-resampling)?", default_yes = TRUE)
  }
}

if (metric == "RV" && maybe_run_bootstrap) {
  
  if (is.na(expr_path) || nchar(expr_path) == 0) {
    expr_path <- pick_file("Choose tidy expression CSV for bootstrap (Day, SampleID, Gene, Abundance)")
  }
  if (!file.exists(expr_path)) stop2("Expression file not found: ", expr_path)
  
  nboot <- NA_integer_
  if (is.na(nboot_arg) || nchar(nboot_arg) == 0) {
    nboot <- ask_int("How many bootstrap replicates?", default = 200)
  } else {
    nboot <- suppressWarnings(as.integer(nboot_arg))
    if (is.na(nboot) || nboot < 1) stop2("Invalid --nboot (must be positive integer).")
  }
  
  msg("Bootstrap settings: expr=", expr_path,
      " | nboot=", nboot,
      " | seed=", seed,
      " | ci=", ci_level,
      " | transform=", transform,
      " | max_genes=", max_genes)
  
  expr_data <- read.csv(expr_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  req_expr <- c("Day", "SampleID", "Gene", "Abundance")
  missing_expr <- setdiff(req_expr, colnames(expr_data))
  if (length(missing_expr) > 0) stop2("Expr file missing required columns: ", paste(missing_expr, collapse = ", "))
  
  expr_data$Day <- suppressWarnings(as.numeric(expr_data$Day))
  expr_data$Abundance <- suppressWarnings(as.numeric(expr_data$Abundance))
  if (any(is.na(expr_data$Day))) stop2("Expr Day has non-numeric values.")
  if (any(is.na(expr_data$Abundance))) stop2("Expr Abundance has non-numeric values (or NA).")
  
  if (!(ref_day %in% unique(expr_data$Day))) stop2("Reference day ", ref_day, " not found in expression data.")
  
  msg("Pre-computing per-day matrices (this can take a while for large gene sets)...")
  master_mats <- vector("list", length(days))
  names(master_mats) <- as.character(days)
  
  debug_day_summary <- data.frame(
    Day = days,
    n_samples = NA_integer_,
    n_genes_pre = NA_integer_,
    n_genes_post = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (idx in seq_along(days)) {
    d <- days[idx]
    
    if (debug) {
      dd <- expr_data[expr_data$Day == d, , drop = FALSE]
      debug_day_summary$n_samples[idx]   <- length(unique(dd$SampleID))
      debug_day_summary$n_genes_pre[idx] <- length(unique(dd$Gene))
    }
    
    master_mats[[as.character(d)]] <- build_day_matrix(expr_data, d, trans = transform, max_g = max_genes)
    
    if (is.null(master_mats[[as.character(d)]])) {
      stop2("Day ", d, " produced NULL matrix after preprocessing (too few genes/samples). ",
            "Try --max_genes=1000 for troubleshooting, and verify each day has >=2 samples.")
    }
    
    if (debug) {
      msg("Day ", d, ": ", nrow(master_mats[[as.character(d)]]), " samples × ",
          ncol(master_mats[[as.character(d)]]), " genes AFTER preprocessing")
      debug_day_summary$n_genes_post[idx] <- ncol(master_mats[[as.character(d)]])
    }
    
    if (d %% 2 == 0) gc()
  }
  
  if (debug && save_debug_files) {
    write.csv(debug_day_summary, file.path(out_dir, "DEBUG_day_matrix_summary.csv"), row.names = FALSE)
    msg("Wrote debug day summary: ", file.path(out_dir, "DEBUG_day_matrix_summary.csv"))
  }
  
  msg("Harmonizing gene set across days (intersection of columns)...")
  gene_sets <- lapply(master_mats, colnames)
  common_genes <- Reduce(intersect, gene_sets)
  
  if (length(common_genes) < 2) {
    stop2(
      "After preprocessing, fewer than 2 genes are shared across all days (",
      length(common_genes), "). RV cannot be computed.\n",
      "Likely causes: aggressive filtering from missingness/zero-variance; days differ in gene coverage.\n",
      "Try: (i) ensure consistent gene coverage across days; (ii) set --max_genes to a smaller consistent subset; ",
      "(iii) relax preprocessing filters."
    )
  }
  
  common_genes <- sort(common_genes)
  for (nm in names(master_mats)) {
    master_mats[[nm]] <- master_mats[[nm]][, common_genes, drop = FALSE]
  }
  msg("Common genes retained across all days: ", length(common_genes))
  
  if (debug && save_debug_files) {
    write.csv(data.frame(Gene = common_genes), file.path(out_dir, "DEBUG_common_genes_across_days.csv"), row.names = FALSE)
    msg("Wrote common gene list: ", file.path(out_dir, "DEBUG_common_genes_across_days.csv"))
  }
  
  ref_mat <- master_mats[[as.character(ref_day)]]
  if (is.null(ref_mat)) stop2("Reference day matrix is NULL after preprocessing.")
  
  msg("Computing observed RV-to-Day", ref_day, " profile...")
  rv_obs <- sapply(days, function(d) calc_rv(master_mats[[as.character(d)]], ref_mat))
  
  set.seed(seed)
  boot_results <- matrix(NA_real_, nrow = nboot, ncol = length(days),
                         dimnames = list(NULL, as.character(days)))
  
  msg("Starting bootstrap loop...")
  for (b in seq_len(nboot)) {
    
    ref_b <- ref_mat[sample.int(nrow(ref_mat), size = nrow(ref_mat), replace = TRUE), , drop = FALSE]
    
    for (d in as.character(days)) {
      target_mat <- master_mats[[d]]
      target_b <- target_mat[sample.int(nrow(target_mat), size = nrow(target_mat), replace = TRUE), , drop = FALSE]
      boot_results[b, d] <- calc_rv(target_b, ref_b)
    }
    
    if (debug && b == 1 && save_debug_files) {
      write.csv(data.frame(Day = days, RV_b1 = as.numeric(boot_results[b, ])),
                file.path(out_dir, "DEBUG_bootstrap_draw1_profile.csv"),
                row.names = FALSE)
      msg("Wrote debug first bootstrap draw: ", file.path(out_dir, "DEBUG_bootstrap_draw1_profile.csv"))
    }
    
    if (b %% max(1, floor(nboot / 10)) == 0) {
      msg("  bootstrap ", b, "/", nboot)
      gc()
    }
  }
  
  ci_low  <- apply(boot_results, 2, quantile, probs = alpha,     na.rm = TRUE)
  ci_high <- apply(boot_results, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
  
  out_df <- data.frame(
    Day = days,
    RV_observed = as.numeric(rv_obs),
    CI_low  = as.numeric(ci_low[as.character(days)]),
    CI_high = as.numeric(ci_high[as.character(days)]),
    RefDay = ref_day,
    nboot = nboot,
    ci_level = ci_level,
    seed = seed,
    transform = transform,
    max_genes = max_genes,
    n_common_genes = length(common_genes),
    stringsAsFactors = FALSE
  )
  
  write.csv(out_df, file.path(out_dir, "RefDay_Profile_RV_bootstrap_CI.csv"), row.names = FALSE)
  msg("Bootstrap complete. Wrote: ", file.path(out_dir, "RefDay_Profile_RV_bootstrap_CI.csv"))
  
  if (save_boot_draws) {
    write.csv(boot_results, file.path(out_dir, "RV_bootstrap_draws_profile_wide.csv"), row.names = FALSE)
    msg("Saved bootstrap draws (wide): ", file.path(out_dir, "RV_bootstrap_draws_profile_wide.csv"))
  }
  
} else {
  if (metric == "RV") msg("RV bootstrap not run (skipped or disabled).")
}

# ------------------------------ Write manifest (JSON + inventory CSV) ------------------------------

dep_pkgs <- unique(c("jsonlite", "knitr", "rmarkdown", "ggplot2", "igraph", "rstudioapi"))
dependencies <- deps_snapshot(dep_pkgs)

manifest_obj <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(
    name = script_name,
    path = script_path,
    full_path = script_full_path
  ),
  input = list(
    pairwise_similarity_csv = normalizePath(in_path, winslash = "/", mustWork = FALSE),
    expression_tidy_csv = if (!is.na(expr_path) && nzchar(expr_path)) normalizePath(expr_path, winslash = "/", mustWork = FALSE) else NA_character_
  ),
  parameters = list(
    metric = metric,
    ref_day = ref_day,
    top_k = top_k,
    bootstrap = list(
      ran = (metric == "RV" && exists("maybe_run_bootstrap") && isTRUE(maybe_run_bootstrap) &&
               !is.na(expr_path) && nzchar(expr_path) && file.exists(expr_path)),
      nboot = if (exists("nboot") && !is.na(nboot)) nboot else NA_integer_,
      seed = seed,
      ci_level = ci_level,
      transform = transform,
      max_genes = max_genes,
      save_boot_draws = save_boot_draws
    ),
    debug = list(
      debug = debug,
      save_debug_files = save_debug_files
    )
  ),
  dependencies = dependencies,
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  ),
  generated_files = list(
    inventory_csv = "Project_Manifest_Files.csv"
  ),
  notes = list(
    script_path_detection_note = ifelse(nzchar(script_path_detection_note), script_path_detection_note, NA_character_)
  )
)

manifest_obj <- sanitize_for_json(manifest_obj)

# atomic write + sanity parse
write_json_atomic(manifest_obj, manifest_json_path)

# update inventory after manifest exists
build_file_inventory(out_dir, manifest_files_csv_path)

# ------------------------------ Quarto QMD + Render (no dplyr) ------------------------------

report_qmd_path  <- file.path(out_dir, paste0(script_name, "_Report.qmd"))
report_html_path <- file.path(out_dir, paste0(script_name, "_Report.html"))

# Embed the header block (verbatim-ish)
header_text <- paste(readLines(textConnection(c(
  "# ======================================================================",
  "# RV_MDS_ANALYSIS.R",
  "# Day–Day Correlation-Structure Similarity: NON-HEATMAP Publication Plots",
  "# + Integrated \"Optimized\" RV Bootstrap (row-resampling; memory-conscious)",
  "#",
  "# PRIMARY INPUT (CSV; long pairwise): Day1, Day2, matrix_corr, frobenius, RV",
  "# OPTIONAL INPUT (for RV bootstrap only; tidy expression): Day, SampleID, Gene, Abundance",
  "#",
  "# OUTPUT (standardized under outputs/<run_id>/):",
  "#   - similarity_long_clean.csv",
  "#   - similarity_wide_matrix_corr.csv",
  "#   - similarity_wide_RV.csv",
  "#   - similarity_wide_frobenius.csv",
  "#   - RefDay_Profile_<metric>_data.csv",
  "#   - Plot_RefDay_Profile_<metric>.png/.pdf",
  "#   - Plot_MDS_<metric>.png/.pdf",
  "#   - Plot_DayNetwork_<metric>.png/.pdf",
  "#   - mds_coords_<metric>.csv",
  "#   - network_edges_<metric>_topk.csv",
  "#   - RefDay_Profile_RV_bootstrap_CI.csv (if metric=RV and bootstrap runs)",
  "#   - RV_bootstrap_draws_profile_wide.csv (optional)",
  "#   - Project_Manifest.json",
  "#   - Project_Manifest_Files.csv",
  "#   - <script_name>_Report.qmd + rendered HTML",
  "#",
  "# IMPORTANT:",
  "#   All if/else logic is written in BLOCK FORM to avoid R's \"unexpected else\" errors.",
  "# ======================================================================"
))), collapse = "\n")

qmd_lines <- c(
  "---",
  paste0("title: \"", script_name, " Report\""),
  paste0("date: \"", run_timestamp, "\""),
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
  "library(knitr)",
  "manifest <- jsonlite::fromJSON('Project_Manifest.json', simplifyVector = FALSE)",
  "inv <- tryCatch(read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE), error = function(e) NULL)",
  "```",
  "",
  "```{r}",
  "meta_tbl <- data.frame(",
  "  Field = c('run_id','run_timestamp','script_name','script_path','script_full_path','outputs_root','output_dir','metric','ref_day','top_k'),",
  "  Value = c(",
  "    manifest[['run_id']],",
  "    manifest[['run_timestamp']],",
  "    manifest[['script']][['name']],",
  "    manifest[['script']][['path']],",
  "    ifelse(is.null(manifest[['script']][['full_path']]) || is.na(manifest[['script']][['full_path']]), 'NA', manifest[['script']][['full_path']]),",
  "    manifest[['outputs']][['outputs_root']],",
  "    manifest[['outputs']][['output_dir']],",
  "    manifest[['parameters']][['metric']],",
  "    manifest[['parameters']][['ref_day']],",
  "    manifest[['parameters']][['top_k']]",
  "  ),",
  "  stringsAsFactors = FALSE",
  ")",
  "knitr::kable(meta_tbl, align = c('l','l'))",
  "```",
  "",
  "## Script header",
  "",
  "```text",
  header_text,
  "```",
  "",
  "## Inputs",
  "",
  "```{r}",
  "in_tbl <- data.frame(",
  "  Input = c('pairwise_similarity_csv','expression_tidy_csv'),",
  "  Path  = c(manifest[['input']][['pairwise_similarity_csv']], manifest[['input']][['expression_tidy_csv']]),",
  "  stringsAsFactors = FALSE",
  ")",
  "knitr::kable(in_tbl, align = c('l','l'))",
  "```",
  "",
  "## Dependencies (packages + versions)",
  "",
  "```{r}",
  "dep <- manifest[['dependencies']]",
  "if (is.null(dep)) {",
  "  cat('No dependency list found in manifest.')",
  "} else {",
  "  pkgs <- character(0); vers <- character(0)",
  "  for (i in seq_along(dep)) {",
  "    item <- dep[[i]]",
  "    pkgs <- c(pkgs, if (!is.null(item[['package']])) as.character(item[['package']]) else NA_character_)",
  "    vers <- c(vers, if (!is.null(item[['version']])) as.character(item[['version']]) else NA_character_)",
  "  }",
  "  dep_tbl <- data.frame(package = pkgs, version = vers, stringsAsFactors = FALSE)",
  "  knitr::kable(dep_tbl, align = c('l','l'))",
  "}",
  "```",
  "",
  "## Generated files",
  "",
  "```{r}",
  "if (is.null(inv)) {",
  "  cat('File inventory could not be read.')",
  "} else {",
  "  knitr::kable(inv[, c('file','size_bytes','mtime')], align = c('l','r','l'))",
  "}",
  "```",
  "",
  "## Analytical logic and formulas",
  "",
  "This analysis consumes a day-by-day similarity table (Day1, Day2) with three metrics:",
  "",
  "- **matrix_corr**: similarity in correlation structure (assumed on a [0,1] scale here).",
  "- **RV**: Robert–Escoufier RV coefficient (a normalized inner product between two matrices).",
  "- **frobenius**: Frobenius distance between matrices (lower means more similar).",
  "",
  "### Symmetrization and diagonals",
  "",
  "Pairs are symmetrized by adding (Day2, Day1) for each (Day1, Day2). Diagonals are enforced as:",
  "",
  "- matrix_corr(Day,Day) = 1",
  "- RV(Day,Day) = 1",
  "- frobenius(Day,Day) = 0",
  "",
  "### Distance transform for MDS",
  "",
  "For similarity metrics (matrix_corr, RV), the distance matrix is computed as:",
  "",
  "$$D = 1 - S$$",
  "",
  "with clamping of S to [0,1]. For Frobenius, the distance matrix is the Frobenius matrix itself.",
  "",
  "MDS is computed via classical scaling (cmdscale) on $D$ to obtain a 2D embedding.",
  "",
  "### Reference-day profile",
  "",
  "Given a reference day $r$, the profile is the column $S_{:,r}$ (or $D_{:,r}$ for Frobenius), plotted over days.",
  "",
  "### Day network (top-k neighbors)",
  "",
  "For each day, the top_k most similar (or smallest distance for Frobenius) other days are selected to create edges; undirected duplicates are removed by keeping the first instance after sorting by edge strength.",
  "",
  "### RV bootstrap (row-resampling; optional; RV metric only)",
  "",
  "If enabled, the script constructs per-day Sample×Gene matrices, applies per-gene median imputation, optional log1p transform, scaling, and then computes RV between gene–gene correlation matrices:",
  "",
  "$$RV(C_1,C_2)=\\frac{\\sum_{ij} C_{1,ij} C_{2,ij}}{\\sqrt{\\left(\\sum_{ij} C_{1,ij}^2\\right)\\left(\\sum_{ij} C_{2,ij}^2\\right)}}$$",
  "",
  "Bootstrap replicates resample rows (samples) within each day matrix with replacement to generate confidence intervals for the RV-to-reference-day profile.",
  "",
  "## Key plots",
  "",
  "```{r}",
  "metric <- manifest[['parameters']][['metric']]",
  "p1 <- paste0('Plot_MDS_', metric, '.png')",
  "p2 <- paste0('Plot_DayNetwork_', metric, '.png')",
  "plot_path <- if (file.exists(p1)) p1 else if (file.exists(p2)) p2 else NA_character_",
  "if (!is.na(plot_path)) {",
  "  knitr::include_graphics(plot_path)",
  "} else {",
  "  cat('No plot image found to embed.')",
  "}",
  "```",
  "",
  "## Interpretation guidance",
  "",
  "- **High RV / matrix_corr to a reference day** indicates that the global gene–gene coordination structure for that day resembles the reference day.",
  "- **Low Frobenius distance to a reference day** indicates the same, expressed as a distance.",
  "- In the MDS map, nearby days share similar global correlation structure, while distant days indicate reconfiguration of coordination architecture.",
  "- The day-network visualization is a local neighborhood summary highlighting each day’s strongest structural affinities under the chosen metric and top_k.",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, con = report_qmd_path)

if (!file.exists(report_qmd_path)) stop2("QMD write failed; file not found at: ", report_qmd_path)

render_ok <- TRUE
render_err <- NULL

# Render with knit_root_dir so relative reads (Project_Manifest.json, plots) work
tryCatch({
  rmarkdown::render(
    input = report_qmd_path,
    output_file = basename(report_html_path),
    output_dir = out_dir,
    quiet = TRUE,
    knit_root_dir = out_dir
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- conditionMessage(e)
})

# Refresh inventory AFTER render (MANDATORY)
build_file_inventory(out_dir, manifest_files_csv_path)

# Update manifest with report status/paths (atomic rewrite)
manifest_obj$notes$report_qmd_path <- normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE)
manifest_obj$notes$report_html_path <- normalizePath(report_html_path, winslash = "/", mustWork = FALSE)
manifest_obj$notes$report_render_success <- render_ok
manifest_obj$notes$report_render_error <- if (!render_ok) render_err else NA_character_
manifest_obj <- sanitize_for_json(manifest_obj)
write_json_atomic(manifest_obj, manifest_json_path)

# ------------------------------ Finish ----------------------------------

msg("Done. Files written to: ", out_dir)
msg("Manifest:")
msg(" - ", manifest_json_path)
msg(" - ", manifest_files_csv_path)

if (render_ok && file.exists(report_html_path)) {
  msg("HTML report created: ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE))
} else {
  msg("REPORT RENDER FAILED. QMD written at: ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE))
  if (!render_ok && !is.null(render_err)) msg("Render error: ", render_err)
}
