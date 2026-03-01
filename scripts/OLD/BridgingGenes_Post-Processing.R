#!/usr/bin/env Rscript
# ============================================================
# BridgingGenes Post-Processing: Quadrants Over Time + Trajectories
#
# PURPOSE
#   Consume outputs from Identify_Bridging_Genes_FullNetworkver2.R
#   (a single run directory), and produce:
#
#   1) Per-(Gene, Day) quadrant membership using global betweenness cutoff
#   2) Gene trajectory strings (Day4:Q -> Day8:Q -> ...)
#   3) Transition/volatility metrics per gene (n_transitions, unique_quadrants)
#   4) Summary tables suitable for Extended Data:
#        - Quadrant counts per gene
#        - Day-wise quadrant composition
#        - Top genes by time spent as integrator/connector/bottleneck
#   5) Optional diagnostic plots (static + plotly RDS)
#   6) QMD + rendered HTML report
#   7) Project_Manifest.json + Project_Manifest_Files.csv
#
# INPUT
#   User selects a prior BridgingGenes run directory that contains at least:
#     - BridgingGenes_SummaryAcrossDays.csv
#   And preferably also:
#     - BridgingGenes_Integrated_Long_AllDays.csv
#   Optionally:
#     - BridgingGenes_RoleQuadrants.csv (contains the bw_cut logic already applied)
#
# OUTPUT POLICY
#   Writes a NEW outputs/<run_id>/ directory under current working dir.
#   Does NOT modify the original run folder.
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
cat("run_id:     ", run_id, "\n", sep = "")
cat("output_dir: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("================================================\n\n")

# ============================================================
# USER INPUT: choose prior BridgingGenes run folder
# ============================================================
cat("\nChoose the PRIOR BridgingGenes run directory (the folder containing BridgingGenes_*.csv outputs)...\n")
in_dir <- choose_directory("Select the prior BridgingGenes run folder")
if (!dir.exists(in_dir)) stop("Input directory does not exist: ", in_dir)
in_dir <- normalizePath(in_dir, winslash = "/", mustWork = FALSE)

cat("Using input run folder:\n  ", in_dir, "\n", sep = "")

# Required + optional inputs
path_summary  <- file.path(in_dir, "BridgingGenes_SummaryAcrossDays.csv")
path_long     <- file.path(in_dir, "BridgingGenes_Integrated_Long_AllDays.csv")
path_roles    <- file.path(in_dir, "BridgingGenes_RoleQuadrants.csv")

if (!file.exists(path_summary)) {
  stop("Required file missing in selected folder:\n  ", path_summary,
       "\nThis post-processing script expects the FIRST script's outputs.")
}

summary <- read.csv(path_summary, stringsAsFactors = FALSE)

# Long table: preferred
integrated_long <- NULL
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
}

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
# Preferred: infer bw_cut from RoleQuadrants.csv if it exists, because it was
# computed from mean_betweenness nonzero quantile (0.85) in your first script.
# If missing, recompute from SummaryAcrossDays using same rule.

betweenness_quantile <- 0.85
bw_cut <- NA_real_
bw_cut_source <- NA_character_

if (file.exists(path_roles)) {
  roles <- read.csv(path_roles, stringsAsFactors = FALSE)
  
  # Recompute from the SAME underlying summary rule rather than trying to parse text:
  # top 15% of nonzero mean_betweenness
  if (!("mean_betweenness" %in% colnames(roles))) {
    # roles should have mean_betweenness; if not, fallback
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
quad_counts <- aggregate(
  Day ~ Gene + Quadrant_day,
  integrated_long,
  length
)
colnames(quad_counts)[3] <- "n_days"

# Day-wise composition
day_comp <- aggregate(
  Gene ~ Day + Quadrant_day,
  integrated_long,
  length
)
colnames(day_comp)[3] <- "n_genes"

# Join volatility with cross-day summary (if genes overlap)
summary_keep <- summary
if (!("Gene" %in% colnames(summary_keep))) stop("SummaryAcrossDays missing Gene column.")
merged_summary <- merge(summary_keep, volatility, by = "Gene", all.x = TRUE)

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
#   (1) Day-wise quadrant composition
#   (2) Top genes by time spent as System-level integrators
# ============================================================
plotly_rds_files <- character()

# Plot 1: Day composition (stacked bars)
p1_png <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay.png")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p1 <- ggplot(day_comp, aes(x = Day, y = n_genes, fill = Quadrant_day)) +
    geom_col() +
    labs(title = "Quadrant composition by day", x = "Day", y = "Number of genes") +
    theme_minimal(base_size = 12)
  
  png(p1_png, width = 1600, height = 1000, res = 150)
  print(p1)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    p1_rds <- file.path(output_dir, "BridgingGenes_QuadrantComposition_ByDay_plotly.rds")
    try({ saveRDS(ggplotly(p1), p1_rds); plotly_rds_files <- c(plotly_rds_files, basename(p1_rds)) }, silent = TRUE)
  }
}

# Plot 2: Top genes by days as integrator
integrator_days <- quad_counts[quad_counts$Quadrant_day == "System-level integrators", , drop = FALSE]
integrator_days <- integrator_days[order(-integrator_days$n_days, integrator_days$Gene), ]
top_k <- 25
integrator_top <- head(integrator_days, top_k)

p2_png <- file.path(output_dir, "BridgingGenes_TopIntegrators_ByDays.png")
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(integrator_top) > 0) {
  integrator_top$Gene <- factor(integrator_top$Gene, levels = rev(integrator_top$Gene))
  p2 <- ggplot(integrator_top, aes(x = Gene, y = n_days)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top ", min(top_k, nrow(integrator_top)), " genes by days classified as System-level integrators"),
         x = "Gene", y = "Days as integrator") +
    theme_minimal(base_size = 12)
  
  png(p2_png, width = 1600, height = 1100, res = 150)
  print(p2)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    p2_rds <- file.path(output_dir, "BridgingGenes_TopIntegrators_ByDays_plotly.rds")
    try({ saveRDS(ggplotly(p2), p2_rds); plotly_rds_files <- c(plotly_rds_files, basename(p2_rds)) }, silent = TRUE)
  }
}

# ============================================================
# Manifest + inventory helpers (same style)
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

manifest <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input = list(
    prior_run_dir = in_dir,
    files_detected = list(
      summary = file.exists(path_summary),
      integrated_long = file.exists(path_long),
      role_quadrants = file.exists(path_roles)
    )
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
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
)

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
           con = manifest_json_path)

manifest_files_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")
build_file_inventory(output_dir, manifest_files_csv_path)

# ============================================================
# Write QMD + render HTML report
# ============================================================
script_name <- "BridgingGenes_PostProcess_QuadrantsOverTime"
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
  paste0('  prior_run_dir: "', in_dir, '"'),
  paste0('  bw_cut: ', as.character(bw_cut)),
  paste0('  bw_cut_source: "', bw_cut_source, '"'),
  paste0('  betweenness_quantile: ', as.character(betweenness_quantile)),
  paste0('  enable_plotly: ', ifelse(enable_plotly, "true", "false")),
  paste0('  runtime_seconds: ', ifelse(enable_runtime_tracking, as.character(runtime_seconds), "null")),
  "---",
  "",
  "## Inputs",
  "",
  "```{r}",
  "man <- jsonlite::fromJSON('Project_Manifest.json')",
  "print(man$input)",
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
  "knitr::kable(head(x, 20))",
  "```",
  "",
  "## Quadrant composition by day",
  "",
  "```{r}",
  "d <- read.csv('BridgingGenes_QuadrantComposition_ByDay.csv', stringsAsFactors = FALSE)",
  "knitr::kable(d[order(d$Day, d$Quadrant_day), ], digits = 3)",
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
  "knitr::kable(head(v2, 25))",
  "```",
  "",
  "## Trajectory strings (selected head)",
  "",
  "```{r}",
  "t <- read.csv('BridgingGenes_QuadrantTrajectories_ByGene.csv', stringsAsFactors = FALSE)",
  "knitr::kable(head(t, 20))",
  "```",
  "",
  "## Top integrators by days",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_TopIntegrators_ByDays.png')) knitr::include_graphics('BridgingGenes_TopIntegrators_ByDays.png')",
  "```",
  "",
  "## Generated files",
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
  "```"
)

writeLines(qmd_lines, con = report_qmd_path)

render_ok <- TRUE
render_err <- NULL
tryCatch({
  rmarkdown::render(
    input = report_qmd_path,
    output_file = basename(report_html_path),
    output_dir = output_dir,
    quiet = TRUE
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- conditionMessage(e)
})

# Refresh inventory after render
build_file_inventory(output_dir, manifest_files_csv_path)

cat("\n================ DONE ================\n")
cat("Input run folder : ", in_dir, "\n", sep = "")
cat("Output run folder: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Key outputs:\n")
cat("  - ", basename(out_long_q), "\n", sep = "")
cat("  - ", basename(out_traj), "\n", sep = "")
cat("  - ", basename(out_vol), "\n", sep = "")
cat("  - ", basename(out_qcnt), "\n", sep = "")
cat("  - ", basename(out_dcomp), "\n", sep = "")
cat("  - ", basename(out_merge), "\n", sep = "")
if (file.exists(report_html_path) && render_ok) {
  cat("HTML report      : ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("HTML report      : FAILED\n")
  cat("QMD written to    : ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (!is.null(render_err)) cat("Render error      : ", render_err, "\n", sep = "")
}
cat("Manifest         : ", basename(manifest_json_path), "\n", sep = "")
cat("=====================================\n")
