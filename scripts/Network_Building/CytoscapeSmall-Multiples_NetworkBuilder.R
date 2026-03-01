#!/usr/bin/env Rscript
# ============================================================
# Cytoscape Small-Multiples Network Builder (FULL + BRIDGING)
#
# PURPOSE
#   Build two sets of Cytoscape networks (one per day):
#     (A) FULL (unrestricted) correlation networks
#     (B) BRIDGING (restricted to BridgingGenes list)
#
# INPUTS
#   - Directory containing: CorrMatrix_Pearson_Day4.csv, Day8, ..., Day20
#     (square correlation matrices with gene names as row/col headers)
#   - BridgingGenes_SummaryAcrossDays_WithVolatility.csv (optional but recommended)
#
# OUTPUTS
#   outputs/<run_id>/
#     cytoscape_full/edges_dayXX.csv
#     cytoscape_bridging/edges_dayXX.csv
#     nodes/nodes_full_union.csv
#     nodes/nodes_bridging.csv
#     Project_Manifest.json
#     Project_Manifest_Files.csv
#     <script_name>_Report.qmd + rendered HTML
#
# NOTES
#   - Networks are undirected; edges are unique (upper triangle).
#   - Thresholding and/or top-N cap is applied separately for FULL and BRIDGING.
#   - QMD/HTML includes provenance (script/input/output) and inventories.
# ============================================================

# ============================================================
# FLAGS
# ============================================================
enable_plotly <- FALSE   # keep FALSE unless you truly want plotly in report
enable_runtime_tracking <- TRUE

start_time <- Sys.time()

# ============================================================
# Dependency bootstrap
# ============================================================
quiet_install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

base_required <- c("jsonlite", "knitr", "tools")
quiet_install_if_missing(base_required)

suppressPackageStartupMessages({
  library(jsonlite)
  library(knitr)
  library(tools)
})

# ============================================================
# Script identity capture (robust)
# ============================================================
resolve_script_path <- function() {
  p <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      rstudioapi::getSourceEditorContext()$path
    } else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  
  NA_character_
}

known_script_filename <- "Make_Cytoscape_SmallMultiples_FULL_and_BRIDGING.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()
if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

cat("\n================ Script Identity ================\n")
cat("script_name: ", script_name, "\n", sep = "")
cat("script_path: ", script_path, "\n", sep = "")
cat("script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
cat("=================================================\n\n")

# ============================================================
# Helpers: choose directory/file (Mac/RStudio-friendly)
# ============================================================
choose_directory <- function(prompt = "Select directory") {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::selectDirectory(caption = prompt), error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  if (requireNamespace("tcltk", quietly = TRUE)) {
    p <- tryCatch(tcltk::tk_choose.dir(caption = prompt), error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  cat("\n", prompt, "\n", sep = "")
  p <- readline("Directory (full path): ")
  if (!nzchar(p)) stop("No directory provided.")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

choose_file_optional <- function(prompt = "Select file (optional; press Enter to skip)") {
  cat("\n", prompt, "\n", sep = "")
  p <- readline("File path (or Enter to skip): ")
  if (!nzchar(p)) return(NA_character_)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}
choose_file <- function(prompt = "Select file") {
  # 1) RStudio file chooser (preferred)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::selectFile(caption = prompt), error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  
  # 2) Mac native fallback (tk)
  if (requireNamespace("tcltk", quietly = TRUE)) {
    p <- tryCatch(tcltk::tk_getOpenFile(), error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  
  # 3) Final fallback: manual path
  cat("\n", prompt, "\n", sep = "")
  p <- readline("File path: ")
  if (!nzchar(p)) stop("No file selected/provided.")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

# ============================================================
# Output policy (MANDATORY): outputs/<run_id>/
# ============================================================
analysis_name <- "Cytoscape_SmallMultiples_FULL_and_BRIDGING"
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id_base <- paste0(analysis_name, "_", run_timestamp)

outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

output_dir <- file.path(outputs_root, run_id_base)
if (dir.exists(output_dir)) {
  k <- 1L
  repeat {
    candidate <- file.path(outputs_root, paste0(run_id_base, "_", sprintf("%03d", k)))
    if (!dir.exists(candidate)) { output_dir <- candidate; break }
    k <- k + 1L
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- basename(output_dir)

cat("\n================ Output Policy ================\n")
cat("run_id:       ", run_id, "\n", sep = "")
cat("outputs_root: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("output_dir:   ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("================================================\n\n")

# ============================================================
# PARAMETERS (edit defaults as desired)
#   FULL network should be more stringent to avoid hairballs.
# ============================================================
# FULL network selection:
full_abs_r_threshold <- 0.85
full_top_edges_per_day <- 2500  # cap; set NA to disable

# BRIDGING network selection:
bridging_abs_r_threshold <- 0.70
bridging_top_edges_per_day <- 1500  # cap; set NA to disable

use_upper_triangle_only <- TRUE

# ============================================================
# USER INPUTS
# ============================================================
cat("Choose directory containing CorrMatrix_Pearson_Day*.csv files...\n")
corr_dir <- choose_directory("Select directory with CorrMatrix_Pearson_Day*.csv")

pattern <- "^CorrMatrix_Pearson_Day[0-9]+\\.csv$"
mat_files <- list.files(corr_dir, pattern = pattern, full.names = TRUE)
if (length(mat_files) == 0) {
  stop("No files matched pattern ", pattern, " in:\n  ", corr_dir)
}

infer_day_from_filename <- function(fn) {
  m <- regexec("Day([0-9]+)\\.", basename(fn))
  reg <- regmatches(basename(fn), m)[[1]]
  if (length(reg) < 2) return(NA_integer_)
  as.integer(reg[2])
}
days <- vapply(mat_files, infer_day_from_filename, integer(1))
if (any(is.na(days))) {
  bad <- mat_files[is.na(days)]
  stop("Could not infer day from one or more filenames:\n  ", paste(basename(bad), collapse = "\n  "))
}
ord <- order(days)
mat_files <- mat_files[ord]
days <- days[ord]

cat("\nFound matrices for days: ", paste(days, collapse = ", "), "\n", sep = "")
cat("Matrix directory: ", normalizePath(corr_dir, winslash="/", mustWork=FALSE), "\n", sep = "")

cat("\nSelect BridgingGenes_SummaryAcrossDays_WithVolatility.csv (required for BRIDGING networks).\n")
bg_path <- choose_file("Select BridgingGenes summary CSV")

bridging_genes <- NULL
bg_nodes <- NULL

if (!file.exists(bg_path)) {
  stop("Selected BridgingGenes file does not exist:\n  ", bg_path)
}

bg_nodes <- read.csv(bg_path, stringsAsFactors = FALSE)
if (!("Gene" %in% colnames(bg_nodes))) {
  stop("BridgingGenes summary missing required column: 'Gene'")
}

bridging_genes <- unique(bg_nodes$Gene)
cat("Loaded bridging genes: n = ", length(bridging_genes), "\n", sep = "")


if (!is.na(bg_path) && file.exists(bg_path)) {
  bg_nodes <- read.csv(bg_path, stringsAsFactors = FALSE)
  if (!("Gene" %in% colnames(bg_nodes))) stop("BridgingGenes summary missing 'Gene' column.")
  bridging_genes <- unique(bg_nodes$Gene)
  cat("Loaded bridging genes: n = ", length(bridging_genes), "\n", sep = "")
} else {
  cat("No bridging gene file provided. BRIDGING networks will not be produced.\n")
}

# Output subdirs
out_full_dir <- file.path(output_dir, "cytoscape_full")
out_brid_dir <- file.path(output_dir, "cytoscape_bridging")
out_nodes_dir <- file.path(output_dir, "nodes")
dir.create(out_full_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_nodes_dir, recursive = TRUE, showWarnings = FALSE)
if (!is.null(bridging_genes)) dir.create(out_brid_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Matrix parsing: handles "rowname column + numeric cols" format
# ============================================================
read_corr_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Expect first column is rownames (gene IDs); remaining columns are numeric correlations
  rn <- df[[1]]
  if (anyDuplicated(rn)) stop("Row identifier column has duplicates in: ", basename(path))
  mat_df <- df[, -1, drop = FALSE]
  
  # coerce to numeric matrix
  mat <- as.matrix(sapply(mat_df, function(x) as.numeric(x)))
  rownames(mat) <- rn
  colnames(mat) <- colnames(mat_df)
  
  mat
}

# ============================================================
# Convert matrix to Cytoscape edge list
# ============================================================
matrix_to_edges <- function(mat, abs_thr, top_n, upper_only = TRUE) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  
  rn <- rownames(mat); cn <- colnames(mat)
  if (is.null(rn) || is.null(cn)) stop("Matrix must have rownames and colnames.")
  if (!identical(rn, cn)) {
    common <- intersect(rn, cn)
    if (length(common) < 3) stop("Row/col names mismatch with too small overlap.")
    mat <- mat[common, common, drop = FALSE]
    rn <- rownames(mat); cn <- colnames(mat)
  }
  
  diag(mat) <- NA_real_
  if (upper_only) mat[lower.tri(mat, diag = FALSE)] <- NA_real_
  
  keep <- which(is.finite(mat) & abs(mat) >= abs_thr, arr.ind = TRUE)
  if (nrow(keep) == 0) {
    return(data.frame(source=character(), target=character(), r=numeric(), abs_r=numeric(), sign=character(),
                      stringsAsFactors = FALSE))
  }
  
  edges <- data.frame(
    source = rn[keep[, 1]],
    target = cn[keep[, 2]],
    r      = mat[keep],
    stringsAsFactors = FALSE
  )
  edges$abs_r <- abs(edges$r)
  edges$sign  <- ifelse(edges$r >= 0, "positive", "negative")
  
  if (!is.na(top_n) && is.finite(top_n) && nrow(edges) > top_n) {
    edges <- edges[order(-edges$abs_r), , drop = FALSE]
    edges <- edges[seq_len(top_n), , drop = FALSE]
  }
  
  edges
}

# ============================================================
# Build two network sets by day
# ============================================================
all_full_genes <- character()
all_brid_genes <- character()

edges_written_full <- character()
edges_written_brid <- character()

for (i in seq_along(mat_files)) {
  f <- mat_files[i]
  day <- days[i]
  cat("Processing Day ", day, ": ", basename(f), "\n", sep = "")
  
  mat <- read_corr_matrix(f)
  
  # ---------------- FULL edges ----------------
  edges_full <- matrix_to_edges(
    mat = mat,
    abs_thr = full_abs_r_threshold,
    top_n = full_top_edges_per_day,
    upper_only = use_upper_triangle_only
  )
  out_full <- file.path(out_full_dir, paste0("edges_day", day, ".csv"))
  write.csv(edges_full, out_full, row.names = FALSE)
  edges_written_full <- c(edges_written_full, out_full)
  all_full_genes <- union(all_full_genes, unique(c(edges_full$source, edges_full$target)))
  
  # ---------------- BRIDGING edges ----------------
  if (!is.null(bridging_genes)) {
    common <- intersect(bridging_genes, intersect(rownames(mat), colnames(mat)))
    if (length(common) < 3) {
      edges_brid <- data.frame(source=character(), target=character(), r=numeric(), abs_r=numeric(), sign=character(),
                               stringsAsFactors = FALSE)
    } else {
      mat_b <- mat[common, common, drop = FALSE]
      edges_brid <- matrix_to_edges(
        mat = mat_b,
        abs_thr = bridging_abs_r_threshold,
        top_n = bridging_top_edges_per_day,
        upper_only = use_upper_triangle_only
      )
    }
    
    out_brid <- file.path(out_brid_dir, paste0("edges_day", day, ".csv"))
    write.csv(edges_brid, out_brid, row.names = FALSE)
    edges_written_brid <- c(edges_written_brid, out_brid)
    all_brid_genes <- union(all_brid_genes, unique(c(edges_brid$source, edges_brid$target)))
  }
}

# ============================================================
# Node tables (for Cytoscape styling)
# ============================================================
# Nodes: FULL union of all genes that appear in FULL edge lists
nodes_full <- data.frame(Gene = sort(unique(all_full_genes)), stringsAsFactors = FALSE)

# Nodes: BRIDGING list (if provided), plus any attributes
nodes_brid <- NULL
if (!is.null(bridging_genes)) {
  # Prefer to export the full bridging list for consistent node set across days
  nodes_brid <- data.frame(Gene = sort(unique(bridging_genes)), stringsAsFactors = FALSE)
  
  # Merge selected columns if present
  keep_cols <- intersect(colnames(bg_nodes), c(
    "Gene","mean_participation","max_participation","mean_betweenness","max_betweenness",
    "mean_degree","n_days_present","unique_quadrants","n_transitions","first_quadrant","last_quadrant"
  ))
  if (length(keep_cols) > 1) {
    bg_sub <- bg_nodes[, keep_cols, drop = FALSE]
    nodes_brid <- merge(nodes_brid, bg_sub, by = "Gene", all.x = TRUE)
  }
}

out_nodes_full <- file.path(out_nodes_dir, "nodes_full_union.csv")
write.csv(nodes_full, out_nodes_full, row.names = FALSE)

out_nodes_brid <- NA_character_
if (!is.null(nodes_brid)) {
  out_nodes_brid <- file.path(out_nodes_dir, "nodes_bridging.csv")
  write.csv(nodes_brid, out_nodes_brid, row.names = FALSE)
}

# ============================================================
# Manifest + inventory helpers
# ============================================================
get_deps <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    ver <- NA_character_
    if (requireNamespace(p, quietly = TRUE)) ver <- as.character(utils::packageVersion(p))
    data.frame(package = p, version = ver, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

build_file_inventory <- function(dir_path, out_csv) {
  f <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) {
    write.csv(data.frame(path=character(), filename=character(), bytes=numeric(), modified=character()),
              out_csv, row.names = FALSE)
    return(invisible(out_csv))
  }
  info <- file.info(f)
  inv <- data.frame(
    path = normalizePath(f, winslash = "/", mustWork = FALSE),
    filename = basename(f),
    bytes = as.numeric(info$size),
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  inv <- inv[order(inv$filename), , drop = FALSE]
  write.csv(inv, out_csv, row.names = FALSE)
  invisible(out_csv)
}

end_time <- Sys.time()
runtime_seconds <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

deps_pkgs <- c("jsonlite", "knitr", "tools")
deps_df <- get_deps(unique(deps_pkgs))

# Track input existence
input_exists <- list(
  corr_dir_exists = dir.exists(corr_dir),
  corr_files_found = length(mat_files),
  bridging_gene_file_provided = !is.null(bridging_genes),
  bridging_gene_file_exists = ifelse(is.na(bg_path), FALSE, file.exists(bg_path))
)

# Track outputs
output_expected <- list(
  nodes_full_union_csv = out_nodes_full,
  nodes_bridging_csv = out_nodes_brid,
  cytoscape_full_dir = out_full_dir,
  cytoscape_bridging_dir = if (is.null(bridging_genes)) NA_character_ else out_brid_dir
)
output_exists <- lapply(output_expected, function(p) {
  if (is.na(p) || is.null(p)) return(FALSE)
  file.exists(p) || dir.exists(p)
})

manifest <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script = list(
    name = script_name,
    path = script_path,
    full_path = if (is.na(script_full)) NA_character_ else script_full
  ),
  input = list(
    corr_dir = normalizePath(corr_dir, winslash = "/", mustWork = FALSE),
    corr_files = lapply(mat_files, function(p) normalizePath(p, winslash = "/", mustWork = FALSE)),
    corr_days = as.integer(days),
    bridging_gene_file = ifelse(is.na(bg_path), NA_character_, normalizePath(bg_path, winslash = "/", mustWork = FALSE)),
    files_detected = input_exists
  ),
  parameters = list(
    full_abs_r_threshold = full_abs_r_threshold,
    full_top_edges_per_day = full_top_edges_per_day,
    bridging_abs_r_threshold = bridging_abs_r_threshold,
    bridging_top_edges_per_day = bridging_top_edges_per_day,
    upper_triangle_only = use_upper_triangle_only,
    enable_runtime_tracking = enable_runtime_tracking,
    runtime_seconds = if (enable_runtime_tracking) runtime_seconds else NA_real_
  ),
  dependencies = jsonlite::fromJSON(jsonlite::toJSON(deps_df, dataframe = "rows", auto_unbox = TRUE)),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    expected_paths = lapply(output_expected, function(p) if (is.na(p)) NA_character_ else normalizePath(p, winslash="/", mustWork=FALSE)),
    files_detected = output_exists
  )
)

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           con = manifest_json_path)

manifest_files_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")
build_file_inventory(output_dir, manifest_files_csv_path)

# ============================================================
# QMD + HTML report (Quarto preferred; robust tables)
# ============================================================
report_qmd_path  <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
report_html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))

qmd_lines <- c(
  "---",
  paste0('title: "', script_name, '"'),
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
  paste0('  run_timestamp: "', manifest$run_timestamp, '"'),
  paste0('  runtime_seconds: ', ifelse(enable_runtime_tracking, as.character(runtime_seconds), "null")),
  "---",
  "",
  "## Provenance (script / input / output)",
  "",
  "```{r}",
  "man <- jsonlite::fromJSON('Project_Manifest.json')",
  "cat('Script name: ', man$script$name, '\\n')",
  "cat('Script path: ', man$script$path, '\\n')",
  "cat('Script full: ', man$script$full_path, '\\n\\n')",
  "cat('Output dir: ', man$outputs$output_dir, '\\n')",
  "cat('Run ID: ', man$run_id, '\\n')",
  "```",
  "",
  "### Input files (matrices)",
  "",
  "```{r}",
  "man <- jsonlite::fromJSON('Project_Manifest.json')",
  "df_in <- data.frame(day = man$input$corr_days, file = unlist(man$input$corr_files), stringsAsFactors = FALSE)",
  "knitr::kable(df_in)",
  "```",
  "",
  "### Parameters",
  "",
  "```{r}",
  "p <- man$parameters",
  "p_df <- data.frame(parameter = names(p), value = unlist(p, use.names = FALSE), stringsAsFactors = FALSE)",
  "knitr::kable(p_df)",
  "```",
  "",
  "## Outputs overview",
  "",
  "### Nodes tables",
  "",
  "```{r}",
  "files <- c('nodes/nodes_full_union.csv', 'nodes/nodes_bridging.csv')",
  "exists <- file.exists(files)",
  "knitr::kable(data.frame(file = files, exists = exists, stringsAsFactors = FALSE))",
  "```",
  "",
  "### Edge lists per day",
  "",
  "```{r}",
  "full_edges <- list.files('cytoscape_full', pattern = '^edges_day[0-9]+\\\\.csv$', full.names = TRUE)",
  "brid_edges <- if (dir.exists('cytoscape_bridging')) list.files('cytoscape_bridging', pattern = '^edges_day[0-9]+\\\\.csv$', full.names = TRUE) else character()",
  "knitr::kable(data.frame(set = c(rep('FULL', length(full_edges)), rep('BRIDGING', length(brid_edges))),",
  "                 file = c(basename(full_edges), basename(brid_edges)), stringsAsFactors = FALSE))",
  "```",
  "",
  "## Cytoscape import recipe (Small multiples by day)",
  "",
  "- Import a day network: **File → Import → Network from File** using `cytoscape_full/edges_dayXX.csv` or `cytoscape_bridging/edges_dayXX.csv`.",
  "- Import node table once: **File → Import → Table from File** using `nodes/nodes_full_union.csv` or `nodes/nodes_bridging.csv`.",
  "- Compute layout on Day 4 network once, then **Copy Layout** and **Paste Layout** onto other day networks to keep node positions fixed.",
  "",
  "## File inventory",
  "",
  "```{r}",
  "inv <- read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "## Runtime",
  "",
  "```{r}",
  "cat('Runtime seconds: ', params$runtime_seconds)",
  "```",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, con = report_qmd_path)
if (!file.exists(report_qmd_path)) stop("QMD write failed: ", report_qmd_path)

render_ok <- TRUE
render_err <- NULL

quarto_bin <- Sys.which("quarto")
if (nzchar(quarto_bin)) {
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  args <- c("render", basename(report_qmd_path), "--to", "html", "--output", basename(report_html_path))
  res <- tryCatch(system2(quarto_bin, args = args, stdout = TRUE, stderr = TRUE),
                  error = function(e) e)
  
  if (inherits(res, "error")) {
    render_ok <- FALSE
    render_err <- conditionMessage(res)
  } else {
    status <- attr(res, "status")
    if (!is.null(status) && is.numeric(status) && status != 0) {
      render_ok <- FALSE
      render_err <- paste0("Quarto exit status = ", status, "\nQuarto output:\n", paste(res, collapse = "\n"))
    } else if (!file.exists(report_html_path)) {
      render_ok <- FALSE
      render_err <- paste0("Quarto did not produce expected HTML: ", report_html_path,
                           "\nQuarto output:\n", paste(res, collapse = "\n"))
    }
  }
} else {
  # Fallback: rmarkdown (may fail for qmd depending on setup)
  quiet_install_if_missing(c("rmarkdown"))
  suppressPackageStartupMessages(library(rmarkdown))
  
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  tryCatch({
    rmarkdown::render(
      input = basename(report_qmd_path),
      output_file = basename(report_html_path),
      output_dir = output_dir,
      quiet = TRUE
    )
    if (!file.exists(report_html_path)) stop("HTML not found after rmarkdown::render: ", report_html_path)
  }, error = function(e) {
    render_ok <<- FALSE
    render_err <<- conditionMessage(e)
  })
}

# Refresh inventory after render
build_file_inventory(output_dir, manifest_files_csv_path)

# Update manifest with report presence
manifest$outputs$files_detected$report_qmd_exists <- file.exists(report_qmd_path)
manifest$outputs$files_detected$report_html_exists <- file.exists(report_html_path)
manifest$outputs$expected_paths$report_qmd <- normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE)
manifest$outputs$expected_paths$report_html <- normalizePath(report_html_path, winslash = "/", mustWork = FALSE)

writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           con = manifest_json_path)

cat("\n================ DONE ================\n")
cat("Output run folder: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("FULL edges dir   : ", normalizePath(out_full_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (!is.null(bridging_genes)) {
  cat("BRIDGING edges   : ", normalizePath(out_brid_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("BRIDGING edges   : (not produced; no bridging gene file)\n")
}
cat("Nodes table FULL : ", normalizePath(out_nodes_full, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (!is.null(bridging_genes)) {
  cat("Nodes table BRID : ", normalizePath(out_nodes_brid, winslash = "/", mustWork = FALSE), "\n", sep = "")
}
cat("Manifest         : ", normalizePath(manifest_json_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Inventory        : ", normalizePath(manifest_files_csv_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (isTRUE(render_ok) && file.exists(report_html_path)) {
  cat("HTML report      : ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("HTML report      : FAILED\n")
  cat("QMD written to   : ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (!is.null(render_err)) cat("Render error     : ", render_err, "\n", sep = "")
}
cat("=====================================\n")
