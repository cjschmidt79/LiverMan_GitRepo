#!/usr/bin/env Rscript
# ======================================================================
# PC1 Maintainer Finder (Within-day PCA loadings across Days 4–20)
#
# PURPOSE
#   Given a wide gene-loading table with per-day PC1 columns (e.g., PC1_D4, PC1_D8, ...),
#   this script identifies genes that "maintain" PC1 across days using a user-selected
#   cutoff definition (Top-N, Top-Percent, or Absolute |loading| threshold), computes
#   persistence across days, and ranks genes by:
#     - mean(|PC1|)
#     - SD(|PC1|)
#     - IQR(|PC1|)
#     - SD(rank(|PC1|) across days)
#     - IQR(rank(|PC1|) across days)
#   Then classifies "core" vs "secondary" maintainers based on persistence + rank stability.
#
# OUTPUT POLICY (MANDATORY)
#   All outputs are written under: outputs/<run_id>/
#   Includes:
#     - Project_Manifest.json
#     - Project_Manifest_Files.csv
#     - tables + plots
#     - <script_name>_Report.qmd + rendered HTML
#
# FLAGS
#   enable_plotly <- TRUE
#   enable_runtime_tracking <- TRUE
# ======================================================================

# -------------------------- Flags (MANDATORY) --------------------------
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

start_time <- if (enable_runtime_tracking) Sys.time() else NULL

# -------------------------- Robust script identity (MANDATORY) --------------------------
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
known_script_filename <- "PC1_Maintainer_Finder.R"
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

cat("\n================ Script Identity ================\n")
cat("script_name : ", script_name, "\n", sep = "")
cat("script_path : ", script_path, "\n", sep = "")
cat("script_full : ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
if (is.na(script_full)) {
  cat("NOTE: Script path detection failed; using fallback known_script_filename.\n")
}
cat("=================================================\n\n")

# -------------------------- Dependencies (install if missing) --------------------------
install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

base_pkgs <- c("jsonlite", "knitr", "rmarkdown")
analysis_pkgs <- c("data.table", "ggplot2")

if (enable_plotly) {
  analysis_pkgs <- c(analysis_pkgs, "plotly", "htmlwidgets", "htmltools")
}

install_if_missing(unique(c(base_pkgs, analysis_pkgs)))

suppressPackageStartupMessages({
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
  library(data.table)
  library(ggplot2)
})

if (enable_plotly) {
  suppressPackageStartupMessages({
    library(plotly)
    library(htmlwidgets)
    library(htmltools)
  })
}

# -------------------------- Output directory policy (MANDATORY) --------------------------
analysis_name <- "PC1_Maintainers"
outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE)

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id <- paste0(analysis_name, "_", run_timestamp)
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE)

log_path <- file.path(output_dir, "execution_log.txt")
writeLines(paste(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "Log initialized. Output folder created:" , output_dir),
           con = log_path)

log_msg <- function(...) {
  msg <- paste0(...)
  ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S] ")
  cat(ts, msg, "\n", sep = "")
  cat(ts, msg, "\n", sep = "", file = log_path, append = TRUE)
}

log_msg("Output directory: ", output_dir)

# -------------------------- Input selection (RStudio picker or CLI) --------------------------
pick_input_csv <- function() {
  # Prefer RStudio picker when available
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    cat("\n>>> ACTION: Select the INPUT CSV (gene loadings table)\n")
    cat("The file picker will open in 2 seconds...\n")
    Sys.sleep(2)
    p <- tryCatch(rstudioapi::selectFile(caption = "Select Input CSV", label = "Open"), error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(p)
  }
  # CLI fallback: use --input= or prompt
  ca <- commandArgs(trailingOnly = FALSE)
  inp_arg <- grep("^--input=", ca, value = TRUE)
  if (length(inp_arg) > 0) {
    p2 <- sub("^--input=", "", inp_arg[1])
    if (nzchar(p2) && file.exists(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  p3 <- readline("Enter full path to input CSV: ")
  if (!nzchar(p3) || !file.exists(p3)) stop("Input CSV not found: ", p3)
  normalizePath(p3, winslash = "/", mustWork = FALSE)
}

input_csv <- pick_input_csv()
input_csv <- normalizePath(input_csv, winslash = "/", mustWork = FALSE)
log_msg("Selected input: ", input_csv)

# -------------------------- Load + validate input --------------------------
dt <- data.table::fread(input_csv, data.table = TRUE)

# Identify gene column
gene_col_candidates <- c("Gene", "gene", "GENE", "Feature", "feature", "ID", "id")
gene_col <- gene_col_candidates[gene_col_candidates %in% names(dt)]
if (length(gene_col) == 0) {
  # fallback: first column
  gene_col <- names(dt)[1]
  log_msg("WARNING: No standard gene column found; using first column as gene id: ", gene_col)
} else {
  gene_col <- gene_col[1]
}
setnames(dt, gene_col, "Gene")

# Identify PC1 day columns: PC1_D4, PC1_4, PC1.Day4, etc.
pc1_cols <- grep("^PC1(_|\\.|\\s)?(D)?\\d+$", names(dt), value = TRUE, ignore.case = TRUE)
if (length(pc1_cols) < 2) {
  # Allow alternative: columns like "PC1_Day04" or "PC1_04"
  pc1_cols <- grep("^PC1.*(4|8|10|12|14|16|18|20)$", names(dt), value = TRUE, ignore.case = TRUE)
}
if (length(pc1_cols) < 2) {
  stop("Could not detect PC1 day columns. Expected columns like PC1_D4, PC1_D8, ... Found: ",
       paste(names(dt), collapse = ", "))
}

# Parse day from column names (extract digits)
extract_day <- function(x) as.integer(gsub(".*?(\\d+)$", "\\1", x))
days <- extract_day(pc1_cols)
ord <- order(days)
pc1_cols <- pc1_cols[ord]
days <- days[ord]

log_msg("Detected PC1 columns: ", paste(pc1_cols, collapse = ", "))
log_msg("Detected days: ", paste(days, collapse = ", "))

# Build matrix of absolute loadings
X <- as.matrix(dt[, ..pc1_cols])
storage.mode(X) <- "numeric"
absX <- abs(X)

# -------------------------- User-selected cutoff definition --------------------------
cat("\n================ Cutoff Selection ================\n")
cat("Choose how to define 'strong PC1 contributor' within each day:\n")
cat("  1) Top N genes by |PC1| per day\n")
cat("  2) Top P percent genes by |PC1| per day\n")
cat("  3) Absolute |PC1| >= threshold per day\n")
choice <- suppressWarnings(as.integer(readline("Enter 1, 2, or 3: ")))
if (is.na(choice) || !(choice %in% c(1,2,3))) stop("Invalid choice.")

cutoff_params <- list()
if (choice == 1) {
  N <- suppressWarnings(as.integer(readline("Top N per day (e.g., 20): ")))
  if (is.na(N) || N < 1) stop("Top N must be >= 1.")
  cutoff_params$type <- "topN"
  cutoff_params$N <- N
} else if (choice == 2) {
  P <- suppressWarnings(as.numeric(readline("Top percent per day (e.g., 10 for top 10%): ")))
  if (is.na(P) || P <= 0 || P >= 100) stop("Percent must be between 0 and 100.")
  cutoff_params$type <- "topPercent"
  cutoff_params$P <- P
} else {
  thr <- suppressWarnings(as.numeric(readline("Absolute loading threshold (e.g., 0.05): ")))
  if (is.na(thr) || thr <= 0) stop("Threshold must be > 0.")
  cutoff_params$type <- "absThreshold"
  cutoff_params$thr <- thr
}

# Persistence threshold
cat("\nHow many days must a gene meet the within-day cutoff to be called a maintainer?\n")
cat("You have ", length(days), " days total.\n", sep = "")
persist_k <- suppressWarnings(as.integer(readline(paste0("Enter required days (e.g., 5): "))))
if (is.na(persist_k) || persist_k < 1 || persist_k > length(days)) stop("Invalid persistence requirement.")

# Core/Secondary split based on rank stability
cat("\nCore vs Secondary classification uses rank stability across days.\n")
rank_sd_cut <- suppressWarnings(as.numeric(readline("Enter SD(rank(|PC1|)) cutoff for CORE (e.g., 200): ")))
if (is.na(rank_sd_cut) || rank_sd_cut <= 0) stop("rank_sd_cut must be > 0.")
rank_iqr_cut <- suppressWarnings(as.numeric(readline("Enter IQR(rank(|PC1|)) cutoff for CORE (e.g., 250): ")))
if (is.na(rank_iqr_cut) || rank_iqr_cut <= 0) stop("rank_iqr_cut must be > 0.")

cat("==================================================\n\n")

log_msg("Cutoff type: ", cutoff_params$type)
log_msg("Cutoff params: ", toJSON(cutoff_params, auto_unbox = TRUE))
log_msg("Persistence requirement (days): ", persist_k)
log_msg("Core rank SD cutoff: ", rank_sd_cut, " | Core rank IQR cutoff: ", rank_iqr_cut)

# -------------------------- Compute within-day cutoff hits --------------------------
n_genes <- nrow(dt)
hit <- matrix(FALSE, nrow = n_genes, ncol = length(days))
colnames(hit) <- paste0("D", days)

for (j in seq_along(days)) {
  v <- absX[, j]
  if (cutoff_params$type == "topN") {
    N <- min(cutoff_params$N, length(v))
    idx <- order(v, decreasing = TRUE)[seq_len(N)]
    hit[idx, j] <- TRUE
  } else if (cutoff_params$type == "topPercent") {
    P <- cutoff_params$P
    N <- max(1, floor(length(v) * (P / 100)))
    idx <- order(v, decreasing = TRUE)[seq_len(N)]
    hit[idx, j] <- TRUE
  } else {
    thr <- cutoff_params$thr
    hit[, j] <- (v >= thr)
  }
}

persistence_days <- rowSums(hit)
is_maintainer <- persistence_days >= persist_k

# -------------------------- Rank stability + loading stability --------------------------
# For each day: rank genes by |PC1| (1 = strongest)
rank_mat <- apply(absX, 2, function(v) rank(-v, ties.method = "average"))
colnames(rank_mat) <- paste0("rank_D", days)

# Per-gene summaries
iqr_vec <- function(x) as.numeric(stats::quantile(x, 0.75, na.rm = TRUE) - stats::quantile(x, 0.25, na.rm = TRUE))

gene_metrics <- data.table(
  Gene = dt$Gene,
  persistence_days = persistence_days,
  maintainer = is_maintainer,
  mean_abs_PC1 = rowMeans(absX, na.rm = TRUE),
  sd_abs_PC1 = apply(absX, 1, stats::sd, na.rm = TRUE),
  iqr_abs_PC1 = apply(absX, 1, iqr_vec),
  mean_rank_abs_PC1 = rowMeans(rank_mat, na.rm = TRUE),
  sd_rank_abs_PC1 = apply(rank_mat, 1, stats::sd, na.rm = TRUE),
  iqr_rank_abs_PC1 = apply(rank_mat, 1, iqr_vec)
)

# Core vs Secondary among maintainers
gene_metrics[, class := "non_maintainer"]
gene_metrics[maintainer == TRUE, class := "secondary"]
gene_metrics[maintainer == TRUE &
               sd_rank_abs_PC1 <= rank_sd_cut &
               iqr_rank_abs_PC1 <= rank_iqr_cut,
             class := "core"]

# Rank outputs (requested)
gene_metrics[, rank_mean_abs_PC1 := frank(-mean_abs_PC1, ties.method = "average")]
gene_metrics[, rank_sd_abs_PC1   := frank(sd_abs_PC1, ties.method = "average")]  # higher sd => higher rank number (can invert if desired)
gene_metrics[, rank_iqr_abs_PC1  := frank(iqr_abs_PC1, ties.method = "average")]

# -------------------------- Write tables --------------------------
out_metrics_csv <- file.path(output_dir, "Gene_PC1_Maintainer_Metrics.csv")
out_core_csv <- file.path(output_dir, "Core_PC1_Maintainers.csv")
out_secondary_csv <- file.path(output_dir, "Secondary_PC1_Maintainers.csv")

data.table::fwrite(gene_metrics[order(rank_mean_abs_PC1)], out_metrics_csv)
data.table::fwrite(gene_metrics[class == "core"][order(rank_mean_abs_PC1)], out_core_csv)
data.table::fwrite(gene_metrics[class == "secondary"][order(rank_mean_abs_PC1)], out_secondary_csv)

log_msg("Wrote: ", out_metrics_csv)
log_msg("Wrote: ", out_core_csv)
log_msg("Wrote: ", out_secondary_csv)

# Also output daywise top contributors (by day) for transparency
topN_for_export <- if (cutoff_params$type == "topN") cutoff_params$N else 20
day_top_list <- list()
for (j in seq_along(days)) {
  v <- absX[, j]
  idx <- order(v, decreasing = TRUE)[seq_len(min(topN_for_export, length(v)))]
  day_top_list[[j]] <- data.table(
    Day = days[j],
    Gene = dt$Gene[idx],
    abs_PC1 = v[idx],
    rank = seq_along(idx)
  )
}
day_top_dt <- rbindlist(day_top_list)
out_daytop_csv <- file.path(output_dir, paste0("TopGenes_PerDay_absPC1_Top", topN_for_export, ".csv"))
data.table::fwrite(day_top_dt, out_daytop_csv)
log_msg("Wrote: ", out_daytop_csv)

# -------------------------- Plots (static always-on; interactive optional) --------------------------
p1 <- ggplot(gene_metrics[maintainer == TRUE], aes(x = persistence_days)) +
  geom_histogram(binwidth = 1) +
  labs(title = "Persistence of strong |PC1| contribution", x = "Days meeting cutoff", y = "Gene count")

p2 <- ggplot(gene_metrics[maintainer == TRUE], aes(x = mean_abs_PC1, y = sd_rank_abs_PC1, shape = class)) +
  geom_point(alpha = 0.7) +
  labs(title = "PC1 maintainers: mean(|PC1|) vs SD(rank(|PC1|))", x = "mean(|PC1|) across days", y = "SD(rank(|PC1|)) across days")

png1 <- file.path(output_dir, "Plot_Persistence_Hist.png")
png2 <- file.path(output_dir, "Plot_MeanAbs_vs_SDRank.png")
ggsave(png1, p1, width = 7, height = 4, dpi = 150)
ggsave(png2, p2, width = 7, height = 4, dpi = 150)
log_msg("Wrote: ", png1)
log_msg("Wrote: ", png2)

# Save RDS for QMD (flat-folder logic)
rds_p1 <- file.path(output_dir, "Plot_Persistence_Hist.rds")
rds_p2 <- file.path(output_dir, "Plot_MeanAbs_vs_SDRank.rds")
saveRDS(p1, rds_p1)
saveRDS(p2, rds_p2)

# Optional Plotly widget saved to HTML (for QMD inclusion / verification)
plotly_html <- NULL
if (enable_plotly) {
  w_ok <- TRUE
  w <- tryCatch({
    ggplotly(p2)
  }, error = function(e) {
    w_ok <<- FALSE
    NULL
  })
  if (w_ok && !is.null(w)) {
    plotly_html <- file.path(output_dir, "Plotly_MeanAbs_vs_SDRank.html")
    tryCatch({
      htmlwidgets::saveWidget(w, file = plotly_html, selfcontained = TRUE)
      log_msg("Wrote: ", plotly_html)
    }, error = function(e) {
      log_msg("WARNING: Could not save plotly widget: ", conditionMessage(e))
      plotly_html <<- NULL
    })
  } else {
    log_msg("WARNING: plotly::ggplotly failed; interactive plot will fall back to static in report.")
  }
}

# -------------------------- Manifest (mandatory) --------------------------
get_pkg_versions <- function(pkgs) {
  ip <- installed.packages()
  out <- lapply(pkgs, function(p) {
    if (p %in% rownames(ip)) {
      list(package = p, version = as.character(ip[p, "Version"]))
    } else {
      list(package = p, version = NA_character_)
    }
  })
  out
}

deps <- unique(c(base_pkgs, analysis_pkgs))
dependencies <- get_pkg_versions(deps)

parameters <- list(
  analysis_name = analysis_name,
  cutoff_type = cutoff_params$type,
  cutoff_params = cutoff_params,
  persistence_required_days = persist_k,
  core_rank_sd_cutoff = rank_sd_cut,
  core_rank_iqr_cutoff = rank_iqr_cut,
  enable_plotly = enable_plotly,
  enable_runtime_tracking = enable_runtime_tracking
)

runtime_seconds <- NULL
if (enable_runtime_tracking) {
  end_time <- Sys.time()
  runtime_seconds <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  parameters$runtime_seconds <- runtime_seconds
}

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
manifest_files_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")

manifest <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA_character_, script_full)
  ),
  input = list(
    input_csv = input_csv,
    detected_gene_column = "Gene",
    detected_pc1_columns = pc1_cols,
    detected_days = as.list(days)
  ),
  parameters = parameters,
  dependencies = dependencies,
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
)

write_json(manifest, manifest_json_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
log_msg("Wrote: ", manifest_json_path)

# -------------------------- File inventory (mandatory; refreshed after render too) --------------------------
build_file_inventory <- function(dir_path, out_csv) {
  files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    data.table::fwrite(data.table(path = character(), filename = character(), bytes = numeric()), out_csv)
    return(invisible(NULL))
  }
  info <- file.info(files)
  inv <- data.table(
    path = normalizePath(files, winslash = "/", mustWork = FALSE),
    filename = basename(files),
    bytes = as.numeric(info$size),
    mtime = as.character(info$mtime)
  )
  data.table::fwrite(inv[order(filename)], out_csv)
}

build_file_inventory(output_dir, manifest_files_csv_path)
log_msg("Wrote: ", manifest_files_csv_path)

# -------------------------- QMD report (mandatory; flat-folder logic) --------------------------
# Embed script header text (initial comment block) if script_full available
read_script_header <- function(path_full) {
  if (is.na(path_full) || !file.exists(path_full)) return("Header unavailable (script path not detectable in this execution mode).")
  lines <- readLines(path_full, warn = FALSE)
  # take until first non-comment non-empty line after header block
  hdr <- character()
  for (ln in lines) {
    if (grepl("^\\s*#", ln) || !nzchar(trimws(ln))) {
      hdr <- c(hdr, ln)
    } else {
      break
    }
  }
  paste(hdr, collapse = "\n")
}

script_header_text <- read_script_header(script_full)

report_qmd_name <- paste0(script_name, "_Report.qmd")
report_html_name <- paste0(script_name, "_Report.html")
report_qmd_path <- file.path(output_dir, report_qmd_name)
report_html_path <- file.path(output_dir, report_html_name)

# QMD uses only filenames; render sets knit_root_dir = output_dir
qmd_lines <- c(
  "---",
  paste0('title: "', script_name, " Report\""),
  "format:",
  "  html:",
  "    toc: true",
  "    embed-resources: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "params:",
  paste0('  run_id: "', run_id, '"'),
  paste0('  output_dir: "', normalizePath(output_dir, winslash = "/", mustWork = FALSE), '"'),
  paste0('  input_csv: "', normalizePath(input_csv, winslash = "/", mustWork = FALSE), '"'),
  paste0('  enable_plotly: ', ifelse(enable_plotly, "true", "false")),
  paste0('  runtime_seconds: ', ifelse(is.null(runtime_seconds), "null", runtime_seconds)),
  "---",
  "",
  "```{r setup}",
  "knitr::opts_knit$set(root.dir = params$output_dir)",
  "```",
  "",
  "## Script header",
  "```",
  script_header_text,
  "```",
  "",
  "## Run metadata",
  "```{r}",
  "meta <- list(",
  "  run_id = params$run_id,",
  "  input_csv = params$input_csv,",
  "  output_dir = params$output_dir,",
  "  enable_plotly = params$enable_plotly,",
  "  runtime_seconds = params$runtime_seconds",
  ")",
  "meta",
  "```",
  "",
  "## Dependencies",
  "```{r}",
  "deps <- jsonlite::fromJSON('Project_Manifest.json')$dependencies",
  "deps",
  "```",
  "",
  "## Key parameters",
  "```{r}",
  "params_list <- jsonlite::fromJSON('Project_Manifest.json')$parameters",
  "params_list",
  "```",
  "",
  "## Analytical logic (definitions)",
  "- **Within-day strong contributor cutoff** is applied to **|PC1 loading|** per day (PCA sign is arbitrary).",
  "- **Persistence** = number of days a gene meets the within-day cutoff.",
  "- **Maintainers** = genes with persistence ≥ user-defined threshold.",
  "- **Rank stability** = SD and IQR of within-day ranks of |PC1| across days.",
  "- **Core maintainers** = maintainers whose rank-stability metrics are below user-defined cutoffs.",
  "",
  "## Outputs produced",
  "```{r}",
  "inv <- data.table::fread('Project_Manifest_Files.csv')",
  "inv",
  "```",
  "",
  "## Tables",
  "```{r}",
  "metrics <- data.table::fread('Gene_PC1_Maintainer_Metrics.csv')",
  "core <- data.table::fread('Core_PC1_Maintainers.csv')",
  "secondary <- data.table::fread('Secondary_PC1_Maintainers.csv')",
  "list(",
  "  n_total_genes = nrow(metrics),",
  "  n_maintainers = sum(metrics$maintainer),",
  "  n_core = nrow(core),",
  "  n_secondary = nrow(secondary)",
  ")",
  "```",
  "",
  "### Core maintainers (top 25 by mean(|PC1|))",
  "```{r}",
  "core[order(rank_mean_abs_PC1)][1:min(25, .N)]",
  "```",
  "",
  "## Plots (static)",
  "```{r}",
  "knitr::include_graphics('Plot_Persistence_Hist.png')",
  "```",
  "",
  "```{r}",
  "knitr::include_graphics('Plot_MeanAbs_vs_SDRank.png')",
  "```",
  "",
  "## Plots (interactive if enabled)",
  "```{r, results='asis'}",
  "if (isTRUE(params$enable_plotly)) {",
  "  rds_files <- list.files(params$output_dir, pattern = '\\\\.(rds)$', full.names = FALSE)",
  "  # Only render a curated set of plots interactively (avoid accidental RDS of non-plots)",
  "  keep <- c('Plot_MeanAbs_vs_SDRank.rds')",
  "  keep <- keep[keep %in% rds_files]",
  "  widgets <- list()",
  "  for (f in keep) {",
  "    p <- tryCatch(readRDS(f), error = function(e) NULL)",
  "    if (is.null(p)) {",
  "      widgets[[length(widgets) + 1]] <- htmltools::tags$p(paste0('Could not read plot RDS: ', f))",
  "    } else {",
  "      w <- tryCatch({",
  "        if (inherits(p, 'plotly') || inherits(p, 'htmlwidget')) p else plotly::ggplotly(p)",
  "      }, error = function(e) NULL)",
  "      if (is.null(w)) {",
  "        widgets[[length(widgets) + 1]] <- htmltools::tags$p(paste0('plotly rendering failed for: ', f))",
  "      } else {",
  "        widgets[[length(widgets) + 1]] <- htmltools::tagList(",
  "          htmltools::tags$h4(f),",
  "          w",
  "        )",
  "      }",
  "    }",
  "  }",
  "  htmltools::tagList(widgets)",
  "} else {",
  "  knitr::asis_output('Plotly disabled.')",
  "}",
  "```",
  "",
  "## Interpretation (what the outputs mean)",
  "- Genes with **high persistence** contribute strongly to the **dominant within-day covariance axis (PC1)** across multiple days.",
  "- Among maintainers, **core** genes show **stable ranks** (low SD/IQR of rank), indicating consistent dominance on PC1.",
  "- **Secondary** maintainers meet persistence but have higher rank variability, indicating context-dependent reweighting.",
  "",
  "## Runtime",
  "```{r}",
  "params$runtime_seconds",
  "```",
  "",
  "## Reproducibility",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, con = report_qmd_path)
log_msg("Wrote QMD: ", report_qmd_path)

# -------------------------- Render QMD to HTML (final step) --------------------------
render_ok <- TRUE
render_err <- NULL

if (enable_runtime_tracking) end_time <- Sys.time()
if (enable_runtime_tracking) {
  runtime_seconds <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  # Update manifest parameter before render
  manifest$parameters$runtime_seconds <- runtime_seconds
  write_json(manifest, manifest_json_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

tryCatch({
  suppressPackageStartupMessages(library(rmarkdown))
  rmarkdown::render(
    input = report_qmd_path,
    output_file = report_html_name,
    output_dir = output_dir,
    quiet = TRUE,
    params = list(
      run_id = run_id,
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      input_csv = normalizePath(input_csv, winslash = "/", mustWork = FALSE),
      enable_plotly = enable_plotly,
      runtime_seconds = runtime_seconds
    ),
    knit_root_dir = output_dir
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- conditionMessage(e)
})

# Refresh inventory after render (mandatory)
build_file_inventory(output_dir, manifest_files_csv_path)

# Final console summary
cat("\n===================== DONE =====================\n")
cat("Outputs directory:\n  ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Manifest:\n  ", normalizePath(manifest_json_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("File inventory:\n  ", normalizePath(manifest_files_csv_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (render_ok && file.exists(report_html_path)) {
  cat("HTML report created:\n  ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("HTML report rendering FAILED.\n")
  cat("QMD path:\n  ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("Error:\n  ", render_err, "\n", sep = "")
}
if (enable_runtime_tracking) {
  cat("Runtime (seconds): ", runtime_seconds, "\n", sep = "")
}
cat("===============================================\n\n")
