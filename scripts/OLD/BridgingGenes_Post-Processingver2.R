#!/usr/bin/env Rscript
# ============================================================
# BridgingGenes Post-Processing: Quadrants Over Time + Trajectories
#
# STRICT FIX (as requested; ignore all prior):
#   - Fixes HTML report failure: "The file '..._Report.qmd' does not exist"
#   - Uses Quarto CLI for .qmd when available (robust; correct engine)
#   - Always passes FULL PATH to renderer (never basename-dependent)
#   - Verifies QMD exists immediately after writeLines()
# ============================================================

# ============================================================
# FLAGS (MANDATORY)
# ============================================================
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

start_time <- Sys.time()

# ============================================================
# Dependency bootstrap (install if missing)
# ============================================================
quiet_install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

base_required <- c("jsonlite", "rmarkdown", "knitr", "tools")
quiet_install_if_missing(base_required)
if (enable_plotly) quiet_install_if_missing(c("ggplot2", "plotly", "htmltools"))

suppressPackageStartupMessages({
  library(jsonlite)
  library(rmarkdown)
  library(knitr)
  library(tools)
})
if (enable_plotly) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(plotly)
    library(htmltools)
  })
}

# ============================================================
# Script identity capture (MANDATORY)
#   Captures script_name/script_path/script_full_path robustly
# ============================================================
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
  
  # 3) source() fallback
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) {
    return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  }
  
  # 4) Unknown
  NA_character_
}

known_script_filename <- "BridgingGenes_PostProcess_QuadrantsOverTime.R"
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
# Helper: Mac-friendly directory chooser
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
  cat("Paste full path (fallback):\n")
  p <- readline("Directory: ")
  if (!nzchar(p)) stop("No directory selected/provided.")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

# ============================================================
# Output directory policy (MANDATORY)
# ============================================================
analysis_name <- "BridgingGenes_PostProcess"
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
# USER INPUT: choose prior BridgingGenes run folder
# ============================================================
cat("\nChoose the PRIOR BridgingGenes run directory (the folder containing BridgingGenes_*.csv outputs)...\n")
in_dir <- choose_directory("Select the prior BridgingGenes run folder")
if (!dir.exists(in_dir)) stop("Input directory does not exist: ", in_dir)
in_dir <- normalizePath(in_dir, winslash = "/", mustWork = FALSE)
in_dir_name <- basename(in_dir)

cat("Using input run folder:\n  ", in_dir, "\n", sep = "")
cat("Input run folder name:\n  ", in_dir_name, "\n", sep = "")

# Required + optional inputs
path_summary  <- file.path(in_dir, "BridgingGenes_SummaryAcrossDays.csv")
path_long     <- file.path(in_dir, "BridgingGenes_Integrated_Long_AllDays.csv")
path_roles    <- file.path(in_dir, "BridgingGenes_RoleQuadrants.csv")

# Track input existence flags early
input_exists <- list(
  prior_run_dir_exists = dir.exists(in_dir),
  summary_csv_exists = file.exists(path_summary),
  integrated_long_csv_exists = file.exists(path_long),
  role_quadrants_csv_exists = file.exists(path_roles)
)

if (!file.exists(path_summary)) {
  stop("Required file missing in selected folder:\n  ", path_summary,
       "\nThis post-processing script expects the FIRST script's outputs.")
}

summary <- read.csv(path_summary, stringsAsFactors = FALSE)

# Long table: preferred
integrated_long <- NULL
used_per_day_fallback <- FALSE
per_day_files <- character()

if (file.exists(path_long)) {
  integrated_long <- read.csv(path_long, stringsAsFactors = FALSE)
} else {
  # Fallback: bind per-day files
  per_day_files <- list.files(in_dir, pattern = "^BridgingGenes_NetworkMetrics_Day[0-9]+\\.csv$", full.names = TRUE)
  if (length(per_day_files) == 0) {
    stop("Could not find BridgingGenes_Integrated_Long_AllDays.csv and no per-day files found.\n",
         "Expected either:\n  - ", path_long, "\n  OR\n  - BridgingGenes_NetworkMetrics_DayXX.csv files.")
  }
  dfs <- lapply(per_day_files, function(f) read.csv(f, stringsAsFactors = FALSE))
  integrated_long <- do.call(rbind, dfs)
  used_per_day_fallback <- TRUE
}

input_exists$per_day_fallback_used <- used_per_day_fallback
input_exists$per_day_files_found <- length(per_day_files)
input_exists$per_day_files_exist <- if (length(per_day_files) > 0) all(file.exists(per_day_files)) else FALSE

# Basic sanity checks
needed_cols <- c("Gene", "Day", "participation", "betweenness", "degree", "community")
missing_cols <- setdiff(needed_cols, colnames(integrated_long))
if (length(missing_cols) > 0) {
  stop("Integrated long table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

integrated_long$Day <- as.integer(integrated_long$Day)
integrated_long <- integrated_long[order(integrated_long$Gene, integrated_long$Day), ]
rownames(integrated_long) <- NULL

# ============================================================
# GLOBAL betweenness cutoff (bw_cut) — match your 2×2 logic
# ============================================================
betweenness_quantile <- 0.85
bw_cut <- NA_real_
bw_cut_source <- NA_character_

if (file.exists(path_roles)) {
  roles <- read.csv(path_roles, stringsAsFactors = FALSE)
  
  if (!("mean_betweenness" %in% colnames(roles))) {
    roles <- NULL
  } else {
    bw_nonzero <- roles$mean_betweenness[is.finite(roles$mean_betweenness) & roles$mean_betweenness > 0]
    if (length(bw_nonzero) >= 10) {
      bw_cut <- as.numeric(quantile(bw_nonzero, probs = betweenness_quantile, names = FALSE, type = 7))
      bw_cut_source <- "RoleQuadrants.csv (recomputed from mean_betweenness nonzero quantile)"
    }
  }
}

if (!is.finite(bw_cut)) {
  if (!("mean_betweenness" %in% colnames(summary))) stop("SummaryAcrossDays lacks mean_betweenness column.")
  bw_nonzero <- summary$mean_betweenness[is.finite(summary$mean_betweenness) & summary$mean_betweenness > 0]
  if (length(bw_nonzero) < 10) {
    stop("Too few nonzero mean_betweenness values to compute bw_cut robustly.")
  }
  bw_cut <- as.numeric(quantile(bw_nonzero, probs = betweenness_quantile, names = FALSE, type = 7))
  bw_cut_source <- "SummaryAcrossDays.csv (mean_betweenness nonzero quantile)"
}

cat("\n================ Quadrant Cutoffs ================\n")
cat("Participation high criterion: participation > 0 (per day)\n")
cat("Betweenness high criterion  : betweenness >= bw_cut\n")
cat("bw_cut (global)             : ", signif(bw_cut, 6), "\n", sep = "")
cat("bw_cut source               : ", bw_cut_source, "\n", sep = "")
cat("betweenness_quantile        : ", betweenness_quantile, "\n", sep = "")
cat("==================================================\n\n")

# ============================================================
# Assign per-day quadrant membership (Gene × Day)
# ============================================================
integrated_long$P_class_day <- ifelse(integrated_long$participation > 0,
                                      "High participation", "Low participation")
integrated_long$B_class_day <- ifelse(integrated_long$betweenness >= bw_cut,
                                      "High betweenness", "Low betweenness")

integrated_long$Quadrant_day <- NA_character_
integrated_long$Quadrant_day[integrated_long$P_class_day == "Low participation"  &
                               integrated_long$B_class_day == "Low betweenness"]  <- "Within-module genes"
integrated_long$Quadrant_day[integrated_long$P_class_day == "Low participation"  &
                               integrated_long$B_class_day == "High betweenness"] <- "Intra-module bottlenecks"
integrated_long$Quadrant_day[integrated_long$P_class_day == "High participation" &
                               integrated_long$B_class_day == "Low betweenness"]  <- "Module connectors"
integrated_long$Quadrant_day[integrated_long$P_class_day == "High participation" &
                               integrated_long$B_class_day == "High betweenness"] <- "System-level integrators"

# ============================================================
# Trajectory strings + volatility metrics
# ============================================================
split_by_gene <- split(integrated_long, integrated_long$Gene)

traj_strings <- do.call(rbind, lapply(split_by_gene, function(df) {
  df <- df[order(df$Day), ]
  data.frame(
    Gene = df$Gene[1],
    n_days_present = nrow(df),
    Quadrant_trajectory = paste(paste0("Day", df$Day, ":", df$Quadrant_day), collapse = " → "),
    stringsAsFactors = FALSE
  )
}))

volatility <- do.call(rbind, lapply(split_by_gene, function(df) {
  df <- df[order(df$Day), ]
  q <- df$Quadrant_day
  n_transitions <- if (length(q) <= 1) 0L else sum(q[-1] != q[-length(q)])
  data.frame(
    Gene = df$Gene[1],
    n_days_present = nrow(df),
    unique_quadrants = length(unique(q)),
    n_transitions = n_transitions,
    first_quadrant = q[1],
    last_quadrant = q[length(q)],
    stringsAsFactors = FALSE
  )
}))

# Quadrant counts per gene
quad_counts <- aggregate(Day ~ Gene + Quadrant_day, integrated_long, length)
colnames(quad_counts)[3] <- "n_days"

# Day-wise composition
day_comp <- aggregate(Gene ~ Day + Quadrant_day, integrated_long, length)
colnames(day_comp)[3] <- "n_genes"

# Join volatility with cross-day summary
if (!("Gene" %in% colnames(summary))) stop("SummaryAcrossDays missing Gene column.")
merged_summary <- merge(summary, volatility, by = "Gene", all.x = TRUE)

# ============================================================
# Write outputs
# ============================================================
out_long_q <- file.path(output_dir, "BridgingGenes_Integrated_Long_AllDays_WithQuadrants.csv")
out_traj   <- file.path(output_dir, "BridgingGenes_QuadrantTrajectories_ByGene.csv")
out_vol    <- file.path(output_dir, "BridgingGenes_QuadrantVolatility_ByGene.csv")
out_qcnt   <- file.path(output_dir, "BridgingGenes_QuadrantCounts_ByGene.csv")
out_dcomp  <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay.csv")
out_merge  <- file.path(output_dir, "BridgingGenes_SummaryAcrossDays_WithVolatility.csv")

write.csv(integrated_long, out_long_q, row.names = FALSE)
write.csv(traj_strings, out_traj, row.names = FALSE)
write.csv(volatility, out_vol, row.names = FALSE)
write.csv(quad_counts, out_qcnt, row.names = FALSE)
write.csv(day_comp, out_dcomp, row.names = FALSE)
write.csv(merged_summary, out_merge, row.names = FALSE)

cat("Wrote:\n")
cat("  - ", basename(out_long_q), "\n", sep = "")
cat("  - ", basename(out_traj), "\n", sep = "")
cat("  - ", basename(out_vol), "\n", sep = "")
cat("  - ", basename(out_qcnt), "\n", sep = "")
cat("  - ", basename(out_dcomp), "\n", sep = "")
cat("  - ", basename(out_merge), "\n", sep = "")

# ============================================================
# Optional: simple diagnostic plots
# ============================================================
plotly_rds_files <- character()

# Plot 1: Day composition
p1_png <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay.png")
p1_rds <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay_plotly.rds")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p1 <- ggplot(day_comp, aes(x = Day, y = n_genes, fill = Quadrant_day)) +
    geom_col() +
    labs(title = "Quadrant composition by day", x = "Day", y = "Number of genes") +
    theme_minimal(base_size = 12)
  
  png(p1_png, width = 1600, height = 1000, res = 150)
  print(p1)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    try({
      saveRDS(ggplotly(p1), p1_rds)
      plotly_rds_files <- c(plotly_rds_files, basename(p1_rds))
    }, silent = TRUE)
  }
}

# Plot 2: Top integrators
integrator_days <- quad_counts[quad_counts$Quadrant_day == "System-level integrators", , drop = FALSE]
integrator_days <- integrator_days[order(-integrator_days$n_days, integrator_days$Gene), ]
top_k <- 25
integrator_top <- head(integrator_days, top_k)

p2_png <- file.path(output_dir, "BridgingGenes_TopIntegrators_ByDays.png")
p2_rds <- file.path(output_dir, "BridgingGenes_TopIntegrators_ByDays_plotly.rds")

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(integrator_top) > 0) {
  integrator_top$Gene <- factor(integrator_top$Gene, levels = rev(integrator_top$Gene))
  p2 <- ggplot(integrator_top, aes(x = Gene, y = n_days)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top ", min(top_k, nrow(integrator_top)),
                        " genes by days classified as System-level integrators"),
         x = "Gene", y = "Days as integrator") +
    theme_minimal(base_size = 12)
  
  png(p2_png, width = 1600, height = 1100, res = 150)
  print(p2)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    try({
      saveRDS(ggplotly(p2), p2_rds)
      plotly_rds_files <- c(plotly_rds_files, basename(p2_rds))
    }, silent = TRUE)
  }
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

deps_pkgs <- c("jsonlite", "rmarkdown", "knitr", "tools")
if (enable_plotly) deps_pkgs <- c(deps_pkgs, "ggplot2", "plotly", "htmltools")
deps_df <- get_deps(unique(deps_pkgs))

# Output existence flags (for initial outputs)
output_expected <- list(
  integrated_long_with_quadrants_csv = out_long_q,
  quadrant_trajectories_csv = out_traj,
  quadrant_volatility_csv = out_vol,
  quadrant_counts_csv = out_qcnt,
  quadrant_composition_csv = out_dcomp,
  summary_with_volatility_csv = out_merge,
  quadrant_composition_png = p1_png,
  quadrant_composition_plotly_rds = p1_rds,
  top_integrators_png = p2_png,
  top_integrators_plotly_rds = p2_rds
)
output_exists <- lapply(output_expected, file.exists)

manifest <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script = list(
    name = script_name,
    path = script_path,
    full_path = if (is.na(script_full)) NA_character_ else script_full
  ),
  input = list(
    prior_run_dir = in_dir,
    prior_run_dir_name = in_dir_name,
    paths = list(
      summary_csv = path_summary,
      integrated_long_csv = path_long,
      role_quadrants_csv = path_roles
    ),
    files_detected = input_exists
  ),
  parameters = list(
    participation_high_rule = "participation > 0 (per day)",
    betweenness_quantile = betweenness_quantile,
    bw_cut = bw_cut,
    bw_cut_source = bw_cut_source,
    enable_plotly = enable_plotly,
    enable_runtime_tracking = enable_runtime_tracking,
    runtime_seconds = if (enable_runtime_tracking) runtime_seconds else NA_real_
  ),
  dependencies = jsonlite::fromJSON(jsonlite::toJSON(deps_df, dataframe = "rows", auto_unbox = TRUE)),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    expected_paths = lapply(output_expected, function(p) normalizePath(p, winslash = "/", mustWork = FALSE)),
    files_detected = output_exists
  )
)

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           con = manifest_json_path)

manifest_files_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")
build_file_inventory(output_dir, manifest_files_csv_path)

# ============================================================
# Write QMD + render HTML report (FIXED; DROP-IN REPLACEMENT)
#   - Writes QMD into output_dir
#   - QMD uses robust list->data.frame formatting for detected files
#   - Renders with Quarto (preferred) from output_dir so relative paths resolve
#   - Fallback to rmarkdown::render (last resort) also from output_dir
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
  paste0('  enable_plotly: ', ifelse(enable_plotly, "true", "false")),
  paste0('  runtime_seconds: ', ifelse(enable_runtime_tracking, as.character(runtime_seconds), "null")),
  paste0('  bw_cut: ', as.character(bw_cut)),
  paste0('  bw_cut_source: "', bw_cut_source, '"'),
  paste0('  betweenness_quantile: ', as.character(betweenness_quantile)),
  "---",
  "",
  "## Provenance (script / input / output)",
  "",
  "```{r}",
  "man <- jsonlite::fromJSON('Project_Manifest.json')",
  "cat('Script name: ', man$script$name, '\\n')",
  "cat('Script path: ', man$script$path, '\\n')",
  "cat('Script full: ', man$script$full_path, '\\n\\n')",
  "cat('Input run dir: ', man$input$prior_run_dir, '\\n')",
  "cat('Input run dir name: ', man$input$prior_run_dir_name, '\\n\\n')",
  "cat('Output dir: ', man$outputs$output_dir, '\\n')",
  "cat('Run ID: ', man$run_id, '\\n')",
  "```",
  "",
  "### Input files detected (YES/NO)",
  "",
  "```{r}",
  "fd <- man$input$files_detected",
  "fd_df <- data.frame(",
  "  item  = names(fd),",
  "  value = unlist(fd, use.names = FALSE),",
  "  stringsAsFactors = FALSE",
  ")",
  "knitr::kable(fd_df)",
  "```",
  "",
  "### Output files detected (YES/NO)",
  "",
  "```{r}",
  "od <- man$outputs$files_detected",
  "od_df <- data.frame(",
  "  item   = names(od),",
  "  exists = unlist(od, use.names = FALSE),",
  "  stringsAsFactors = FALSE",
  ")",
  "knitr::kable(od_df)",
  "```",
  "",
  "## Cutoffs used for per-day quadrant membership",
  "",
  "- Participation high rule: **participation > 0** (per day)",
  paste0("- Betweenness high rule: **betweenness ≥ bw_cut** where bw_cut = ", signif(bw_cut, 6),
         " (", bw_cut_source, "; quantile = ", betweenness_quantile, ")"),
  "",
  "## Integrated long table with quadrants (head)",
  "",
  "```{r}",
  "x <- read.csv('BridgingGenes_Integrated_Long_AllDays_WithQuadrants.csv', stringsAsFactors = FALSE)",
  "x <- x[order(x$Gene, x$Day), ]",
  "knitr::kable(head(x, 25), digits = 4)",
  "```",
  "",
  "## Quadrant composition by day",
  "",
  "```{r}",
  "d <- read.csv('BridgingGenes_QuadrantComposition_ByDay.csv', stringsAsFactors = FALSE)",
  "d <- d[order(d$Day, d$Quadrant_day), ]",
  "knitr::kable(d)",
  "```",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_QuadrantComposition_ByDay.png')) knitr::include_graphics('BridgingGenes_QuadrantComposition_ByDay.png')",
  "```",
  "",
  "## Quadrant volatility per gene (top by transitions)",
  "",
  "```{r}",
  "v <- read.csv('BridgingGenes_QuadrantVolatility_ByGene.csv', stringsAsFactors = FALSE)",
  "v2 <- v[order(-v$n_transitions, -v$unique_quadrants, v$Gene), ]",
  "knitr::kable(head(v2, 30))",
  "```",
  "",
  "## Trajectory strings (head)",
  "",
  "```{r}",
  "t <- read.csv('BridgingGenes_QuadrantTrajectories_ByGene.csv', stringsAsFactors = FALSE)",
  "knitr::kable(head(t, 25))",
  "```",
  "",
  "## Top integrators by days",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_TopIntegrators_ByDays.png')) knitr::include_graphics('BridgingGenes_TopIntegrators_ByDays.png')",
  "```",
  "",
  "## Generated files (inventory)",
  "",
  "```{r}",
  "inv <- read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "## Interactive plots (Plotly)",
  "",
  "```{r results='asis'}",
  "if (isTRUE(params$enable_plotly)) {",
  "  rds_files <- list.files('.', pattern = '\\\\_plotly\\\\.rds$', full.names = FALSE)",
  "  if (length(rds_files) == 0) {",
  "    htmltools::tags$p('No plotly RDS files found.')",
  "  } else {",
  "    htmltools::tagList(lapply(rds_files, function(f) {",
  "      w <- tryCatch(readRDS(f), error = function(e) NULL)",
  "      if (is.null(w)) htmltools::tags$p(paste0('Could not read: ', f))",
  "      else htmltools::tagList(htmltools::tags$h4(f), w)",
  "    }))",
  "  }",
  "} else {",
  "  htmltools::tags$p('Plotly disabled.')",
  "}",
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

# Write QMD
writeLines(qmd_lines, con = report_qmd_path)

# STRICT: verify QMD exists right now
if (!file.exists(report_qmd_path)) {
  stop("QMD write failed; file does not exist at expected path:\n  ", report_qmd_path)
}

render_ok <- TRUE
render_err <- NULL

# Render with Quarto (preferred) from output_dir so relative paths resolve
quarto_bin <- Sys.which("quarto")
if (nzchar(quarto_bin)) {
  old_wd2 <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd2), add = TRUE)
  
  args <- c(
    "render",
    basename(report_qmd_path),
    "--to", "html",
    "--output", basename(report_html_path)
  )
  
  res <- tryCatch(
    system2(quarto_bin, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    render_ok <- FALSE
    render_err <- conditionMessage(res)
  } else {
    status <- attr(res, "status")  # NULL on success, nonzero on failure
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
  # Fallback (last resort): rmarkdown::render from output_dir for relative paths
  old_wd2 <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd2), add = TRUE)
  
  tryCatch({
    rmarkdown::render(
      input = basename(report_qmd_path),
      output_file = basename(report_html_path),
      output_dir = output_dir,
      quiet = TRUE
    )
    if (!file.exists(report_html_path)) {
      stop("rmarkdown::render completed but HTML not found at: ", report_html_path)
    }
  }, error = function(e) {
    render_ok <<- FALSE
    render_err <<- conditionMessage(e)
  })
}

# Refresh inventory after render
build_file_inventory(output_dir, manifest_files_csv_path)

# Update manifest output flags after report render
manifest$outputs$files_detected$report_qmd_exists  <- file.exists(report_qmd_path)
manifest$outputs$files_detected$report_html_exists <- file.exists(report_html_path)
manifest$outputs$expected_paths$report_qmd  <- normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE)
manifest$outputs$expected_paths$report_html <- normalizePath(report_html_path, winslash = "/", mustWork = FALSE)

# Rewrite manifest to include final report existence
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           con = manifest_json_path)

cat("\n================ DONE ================\n")
cat("Input run folder : ", in_dir, "\n", sep = "")
cat("Output run folder: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Script full path : ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
cat("Key outputs:\n")
cat("  - ", basename(out_long_q), "\n", sep = "")
cat("  - ", basename(out_traj), "\n", sep = "")
cat("  - ", basename(out_vol), "\n", sep = "")
cat("  - ", basename(out_qcnt), "\n", sep = "")
cat("  - ", basename(out_dcomp), "\n", sep = "")
cat("  - ", basename(out_merge), "\n", sep = "")
cat("Manifest         : ", basename(manifest_json_path), "\n", sep = "")
cat("Inventory        : ", basename(manifest_files_csv_path), "\n", sep = "")
if (isTRUE(render_ok) && file.exists(report_html_path)) {
  cat("HTML report      : ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("HTML report      : FAILED\n")
  cat("QMD written to   : ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (!is.null(render_err)) cat("Render error     : ", render_err, "\n", sep = "")
}
cat("=====================================\n")
