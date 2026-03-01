#!/usr/bin/env Rscript
# ============================================================
# Identify Bridging Genes from Correlation Matrices
#
# PURPOSE
#   - Automatically detect correlation-matrix CSVs
#   - Infer developmental day from filenames
#   - Build fixed-density gene–gene networks
#   - Compute:
#       * Participation coefficient
#       * Betweenness centrality
#       * Degree (NEW: integrated cross-day tables requested)
#   - Summarize bridging behavior across days
#   - Draw 2×2 schematic:
#       Participation (low/high) vs Betweenness (low/high)
#
# ASSUMPTIONS
#   - Each CSV is a symmetric correlation matrix
#   - Row names = column names = gene symbols
#   - Filenames contain "DayXX" (e.g., Day14, Day16)
#
# OUTPUT (written into outputs/<run_id>/)
#   - BridgingGenes_NetworkMetrics_DayXX.csv (per-day)
#   - BridgingGenes_ExcludedGenes_ByDay.csv (QC; genes excluded per day)
#   - BridgingGenes_SummaryAcrossDays.csv (cross-day summary)
#   - BridgingGenes_RoleQuadrants.csv (summary + quadrant labels)
#   - BridgingGenes_RoleQuadrants_2x2.png / .pdf (schematic)
#   - NEW:
#       * BridgingGenes_Integrated_Long_AllDays.csv
#       * BridgingGenes_Integrated_Wide_Participation.csv
#       * BridgingGenes_Integrated_Wide_Betweenness.csv
#       * BridgingGenes_Integrated_Wide_Degree.csv
#       * BridgingGenes_Trajectories_TopByMeanParticipation.png
#       * (optional) BridgingGenes_Trajectories_TopByMeanParticipation_plotly.rds
#   - <script_name>_Report.qmd + <script_name>_Report.html
#   - Project_Manifest.json + Project_Manifest_Files.csv
# ============================================================

# ============================================================
# FLAGS (MANDATORY)
# ============================================================
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

# Runtime tracking start (MANDATORY if enabled)
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

base_required <- c("igraph", "jsonlite", "rmarkdown", "knitr", "tools")
quiet_install_if_missing(base_required)
if (enable_plotly) quiet_install_if_missing(c("ggplot2", "plotly", "htmltools"))

suppressPackageStartupMessages({
  library(igraph)
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
# Script identity capture (MANDATORY BLOCK; do not edit)
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
  
  # 3) source() fallback (sometimes present)
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) {
    return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  }
  
  # 4) Unknown
  NA_character_
}

# If the script filename is known, set it here as a fallback:
known_script_filename <- "Identify_Bridging_Genes_FullNetworkver2.R"
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

# Console summary
cat("\n================ Script Identity ================\n")
cat("script_name:  ", script_name, "\n", sep = "")
cat("script_path:  ", script_path, "\n", sep = "")
cat("script_full:  ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
if (is.na(script_full)) {
  cat("NOTE: Could not detect script path in this execution mode; using fallback known_script_filename.\n")
}
cat("=================================================\n\n")

# ============================================================
# Output directory policy (MANDATORY)
#   - Always write under ./outputs/<run_id>/
#   - Ensure run_id is unique even if launched repeatedly within same second
# ============================================================
analysis_name <- "BridgingGenes"
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_id_base <- paste0(analysis_name, "_CorrMatrices_", run_timestamp)

outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

# Make output_dir unique (no overwrite)
output_dir <- file.path(outputs_root, run_id_base)
if (dir.exists(output_dir)) {
  # If launched multiple times in the same second (or folder already exists), append a numeric suffix.
  k <- 1L
  repeat {
    candidate <- file.path(outputs_root, paste0(run_id_base, "_", sprintf("%03d", k)))
    if (!dir.exists(candidate)) {
      output_dir <- candidate
      break
    }
    k <- k + 1L
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Final run_id must match the folder name actually used
run_id <- basename(output_dir)

cat("\n================ Output Policy ================\n")
cat("outputs_root: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("run_id:       ", run_id, "\n", sep = "")
cat("output_dir:   ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("================================================\n\n")

# ============================================================
# Helper: Mac-friendly directory chooser (INPUT only)
# ============================================================
choose_directory <- function(prompt = "Select directory containing correlation matrix CSV files") {
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
# USER INPUT
# ============================================================
cat("\nChoose directory containing correlation matrix CSV files...\n")
mat_dir <- choose_directory("Select the folder that contains correlation-matrix CSVs")

if (!dir.exists(mat_dir)) stop("Directory does not exist: ", mat_dir)
cat("Using input directory:\n  ", mat_dir, "\n", sep = "")

cat("\nOptional filename pattern [default = '.csv']:\n")
pat <- readline("Pattern: ")
if (!nzchar(pat)) pat <- "\\.csv$"

files <- list.files(mat_dir, pattern = pat, full.names = TRUE)
if (length(files) == 0) stop("No files found matching pattern in ", mat_dir, " using pattern: ", pat)

cat("\nFound correlation matrices:\n")
for (f in files) cat("  - ", basename(f), "\n", sep = "")

# ----------------------------
# Infer day from filename
# ----------------------------
infer_day <- function(fname) {
  m <- regmatches(fname, regexpr("Day[0-9]+", fname))
  if (length(m) == 0 || m == "") return(NA_integer_)
  as.integer(sub("Day", "", m))
}

days <- sapply(basename(files), infer_day)
if (any(is.na(days))) {
  bad <- basename(files)[is.na(days)]
  stop("Could not infer Day from one or more filenames.\n",
       "Expected pattern like 'Day14'. Problem files:\n  - ",
       paste(bad, collapse = "\n  - "))
}

ord_days <- order(days)
files <- files[ord_days]
days  <- days[ord_days]

# Input identity fields (captured once; used in manifest + QMD)
mat_dir_norm <- normalizePath(mat_dir, winslash = "/", mustWork = FALSE)
mat_dir_name <- basename(mat_dir_norm)
files_full   <- normalizePath(files, winslash = "/", mustWork = FALSE)
files_base   <- basename(files)
files_by_day <- stats::setNames(files_base, as.character(days))

cat("\n================ Input Summary ================\n")
cat("correlation_matrix_dir:      ", mat_dir_norm, "\n", sep = "")
cat("correlation_matrix_dir_name: ", mat_dir_name, "\n", sep = "")
cat("n_files:                     ", length(files_full), "\n", sep = "")
cat("files (basename):\n")
for (nm in files_base) cat("  - ", nm, "\n", sep = "")
cat("================================================\n\n")

# ============================================================
# Parameters (keep analysis logic unchanged)
# ============================================================
edge_density <- 0.01
use_weighted_paths <- FALSE
na_frac_thr <- 0.05
betweenness_quantile <- 0.85
n_exemplars_per_quadrant <- 6

# NEW (integrated tables + plot controls)
top_n_for_trajectory_plot <- 25  # number of genes to show in trajectory plot (ranked by mean participation)
trajectory_plot_png_width <- 1600
trajectory_plot_png_height <- 1000
trajectory_plot_png_res <- 150

# ============================================================
# Analysis helpers
# ============================================================
build_graph_fixed_density <- function(R, density) {
  genes <- colnames(R)
  iu <- which(upper.tri(R, diag = FALSE), arr.ind = TRUE)
  rvals <- abs(R[iu])
  k <- max(1, floor(density * length(rvals)))
  ord <- order(rvals, decreasing = TRUE)[seq_len(k)]
  edges <- iu[ord, , drop = FALSE]
  
  g <- graph_from_edgelist(
    cbind(genes[edges[, 1]], genes[edges[, 2]]),
    directed = FALSE
  )
  E(g)$weight <- rvals[ord]
  g
}

participation_coeff <- function(g, membership_vec) {
  vnames <- V(g)$name
  if (is.null(vnames) || any(!nzchar(vnames))) stop("Graph vertices must have non-empty names.")
  
  if (is.null(names(membership_vec))) stop("membership_vec must be a named vector (names = vertex names).")
  if (!all(vnames %in% names(membership_vec))) {
    missing <- setdiff(vnames, names(membership_vec))
    stop("membership_vec missing membership for vertices: ", paste(head(missing, 10), collapse = ", "),
         if (length(missing) > 10) " ...")
  }
  
  ki <- degree(g)
  if (is.null(names(ki))) names(ki) <- vnames
  
  P <- numeric(length(vnames))
  names(P) <- vnames
  
  for (vid in seq_along(vnames)) {
    vname <- vnames[vid]
    deg <- ki[vname]
    if (is.na(deg) || deg == 0) { P[vname] <- 0; next }
    
    nbrs <- neighbors(g, vid)
    nbr_names <- V(g)[nbrs]$name
    if (length(nbr_names) == 0) { P[vname] <- 0; next }
    
    comms <- membership_vec[nbr_names]
    tab <- table(comms)
    frac <- as.numeric(tab) / deg
    P[vname] <- 1 - sum(frac^2)
  }
  
  P
}

plot_2x2 <- function(device_fun, out_path, width = 10, height = 7,
                     ex1, ex2, ex3, ex4, bw_cut, betweenness_quantile) {
  device_fun(out_path, width = width, height = height)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  
  par(mar = c(1.8, 1.8, 3.2, 1.2))
  plot.new()
  plot.window(xlim = c(0, 2), ylim = c(0, 2))
  
  rect(0, 0, 2, 2, border = "black", lwd = 1.3)
  segments(1, 0, 1, 2, lwd = 1.0)
  segments(0, 1, 2, 1, lwd = 1.0)
  
  title(main = "Gene network roles (2×2): Participation vs Betweenness", cex.main = 1.25)
  
  text(1, 2.06, "High betweenness", cex = 0.95)
  text(1, -0.08, "Low betweenness", cex = 0.95)
  text(-0.08, 1, "Low participation", srt = 90, cex = 0.95)
  text(2.08, 1, "High participation", srt = 90, cex = 0.95)
  
  panel_text <- function(cx, cy, heading, genes) {
    text(cx, cy + 0.40, heading, font = 2, cex = 0.95)
    if (length(genes) == 0) {
      text(cx, cy + 0.10, "(none)", cex = 0.85)
    } else {
      text(cx, cy + 0.05, paste(genes, collapse = "\n"), cex = 0.85)
    }
  }
  
  panel_text(0.5, 1.5, "Intra-module bottlenecks", ex2)
  panel_text(1.5, 1.5, "System-level integrators", ex4)
  panel_text(0.5, 0.5, "Within-module genes", ex1)
  panel_text(1.5, 0.5, "Module connectors", ex3)
  
  mtext(
    paste0(
      "Cutoffs: participation > 0 = High; betweenness High = top ",
      round((1 - betweenness_quantile) * 100), "% of nonzero mean_betweenness (cut = ",
      signif(bw_cut, 3), ")."
    ),
    side = 1, line = 0.4, cex = 0.82
  )
}

# ============================================================
# NEW: integrated table builders
# ============================================================
make_integrated_long <- function(per_day_results) {
  # Bind per-day results into a single long table.
  # Enforce column set and ordering.
  out <- do.call(rbind, per_day_results)
  out$Day <- as.integer(out$Day)
  out <- out[order(out$Gene, out$Day), ]
  out <- out[, c("Gene", "Day", "participation", "betweenness", "degree", "community"), drop = FALSE]
  rownames(out) <- NULL
  out
}

make_wide_metric <- function(long_df, metric_col, day_levels) {
  # Base-R reshape to wide: Gene rows, Day columns.
  df <- long_df[, c("Gene", "Day", metric_col), drop = FALSE]
  colnames(df)[3] <- "value"
  
  df$Day <- as.integer(df$Day)
  df <- df[!is.na(df$value), , drop = FALSE]
  
  wide <- reshape(df, idvar = "Gene", timevar = "Day", direction = "wide")
  
  # Column naming: value.4 -> Day4
  cn <- colnames(wide)
  cn <- sub("^value\\.", "Day", cn)
  colnames(wide) <- cn
  
  # Ensure all day columns exist, in desired order
  desired <- paste0("Day", day_levels)
  for (d in desired) {
    if (!d %in% colnames(wide)) wide[[d]] <- NA_real_
  }
  wide <- wide[, c("Gene", desired), drop = FALSE]
  wide
}

# ============================================================
# Containers
# ============================================================
per_day_results <- list()
all_genes <- character()
excluded_genes <- list()

# ============================================================
# Main loop
# ============================================================
for (i in seq_along(files)) {
  
  f <- files[i]
  day <- days[i]
  cat("\nProcessing Day ", day, " (", basename(f), ")\n", sep = "")
  
  R <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  
  if (nrow(R) < 2 || ncol(R) < 2) {
    warning("Matrix too small in file: ", basename(f), " (", nrow(R), "x", ncol(R), "). Skipping.")
    next
  }
  
  if (!all(colnames(R) == rownames(R))) stop("Matrix is not symmetric with matching row/column names: ", f)
  
  suppressWarnings({ mode(R) <- "numeric" })
  diag(R) <- 1
  
  if (any(!is.finite(R))) {
    na_frac <- rowMeans(!is.finite(R))
    bad_genes <- names(na_frac)[na_frac > na_frac_thr]
    
    cat("  Non-finite entries present.\n", sep = "")
    cat("  Genes with NA-fraction >", na_frac_thr, ": ", length(bad_genes), "\n", sep = "")
    if (length(bad_genes) > 0) cat("  Example excluded: ", paste(head(bad_genes, 8), collapse = ", "), "\n", sep = "")
    
    excluded_genes[[as.character(day)]] <- data.frame(
      Day = day,
      Gene = bad_genes,
      NA_fraction = unname(na_frac[bad_genes]),
      Reason = paste0("Non-finite correlation fraction > ", na_frac_thr),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
    keep <- setdiff(colnames(R), bad_genes)
    if (length(keep) < 2) {
      warning("After excluding invalid genes (NA_fraction>", na_frac_thr,
              "), fewer than 2 genes remain for Day ", day, ". Skipping.")
      next
    }
    
    R <- R[keep, keep, drop = FALSE]
    diag(R) <- 1
  }
  
  g <- build_graph_fixed_density(R, edge_density)
  
  if (vcount(g) < 2 || ecount(g) < 1) {
    warning("Graph too small after thresholding for file: ", basename(f),
            " (vcount=", vcount(g), ", ecount=", ecount(g), "). Skipping day ", day, ".")
    next
  }
  
  comm <- cluster_louvain(g, weights = E(g)$weight)
  memb <- membership(comm)
  
  vnames <- V(g)$name
  if (is.null(names(memb))) names(memb) <- vnames
  
  part <- participation_coeff(g, memb)
  
  if (use_weighted_paths) {
    betw <- betweenness(g, weights = 1 / E(g)$weight, normalized = TRUE)
  } else {
    betw <- betweenness(g, normalized = TRUE)
  }
  
  deg <- degree(g)
  
  if (is.null(names(betw))) names(betw) <- vnames
  if (is.null(names(deg)))  names(deg)  <- vnames
  
  df <- data.frame(
    Gene = names(part),
    Day = day,
    participation = unname(as.numeric(part)),
    betweenness = unname(betw[names(part)]),
    degree = unname(deg[names(part)]),
    community = unname(memb[names(part)]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  per_day_results[[as.character(day)]] <- df
  all_genes <- union(all_genes, df$Gene)
  
  out_csv <- file.path(output_dir, paste0("BridgingGenes_NetworkMetrics_Day", day, ".csv"))
  write.csv(df, out_csv, row.names = FALSE)
  cat("  Wrote: ", basename(out_csv), "\n", sep = "")
}

if (length(per_day_results) == 0) stop("No per-day results were produced (all days skipped or failed).")

# ----------------------------
# Write excluded gene audit table
# ----------------------------
excluded_csv <- NULL
if (length(excluded_genes) > 0) {
  excl_df <- do.call(rbind, excluded_genes)
  excl_df <- excl_df[order(excl_df$Day, -excl_df$NA_fraction, excl_df$Gene), ]
  excluded_csv <- file.path(output_dir, "BridgingGenes_ExcludedGenes_ByDay.csv")
  write.csv(excl_df, excluded_csv, row.names = FALSE)
  cat("  Wrote: ", basename(excluded_csv), "\n", sep = "")
}

# ============================================================
# NEW: Integrated cross-day tables (LONG + WIDE)
# ============================================================
cat("\nBuilding integrated cross-day tables (participation, betweenness, degree)\n")

# Long table
integrated_long <- make_integrated_long(per_day_results)
out_long <- file.path(output_dir, "BridgingGenes_Integrated_Long_AllDays.csv")
write.csv(integrated_long, out_long, row.names = FALSE)
cat("  Wrote: ", basename(out_long), "\n", sep = "")

# Wide tables (per metric)
day_levels <- sort(unique(integrated_long$Day))

wide_part <- make_wide_metric(integrated_long, "participation", day_levels)
wide_betw <- make_wide_metric(integrated_long, "betweenness", day_levels)
wide_deg  <- make_wide_metric(integrated_long, "degree", day_levels)

out_wide_part <- file.path(output_dir, "BridgingGenes_Integrated_Wide_Participation.csv")
out_wide_betw <- file.path(output_dir, "BridgingGenes_Integrated_Wide_Betweenness.csv")
out_wide_deg  <- file.path(output_dir, "BridgingGenes_Integrated_Wide_Degree.csv")

write.csv(wide_part, out_wide_part, row.names = FALSE)
write.csv(wide_betw, out_wide_betw, row.names = FALSE)
write.csv(wide_deg,  out_wide_deg,  row.names = FALSE)

cat("  Wrote: ", basename(out_wide_part), "\n", sep = "")
cat("  Wrote: ", basename(out_wide_betw), "\n", sep = "")
cat("  Wrote: ", basename(out_wide_deg), "\n", sep = "")

# ============================================================
# Cross-day summary (unchanged logic)
# ============================================================
cat("\nBuilding cross-day summary\n")

summary <- data.frame(
  Gene = all_genes,
  mean_participation = NA_real_,
  max_participation  = NA_real_,
  mean_betweenness   = NA_real_,
  max_betweenness    = NA_real_,
  mean_degree        = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(all_genes)) {
  gsym <- all_genes[i]
  
  pvals <- c()
  bvals <- c()
  dvals <- c()
  
  for (df in per_day_results) {
    if (gsym %in% df$Gene) {
      row <- df[df$Gene == gsym, , drop = FALSE]
      pvals <- c(pvals, row$participation)
      bvals <- c(bvals, row$betweenness)
      dvals <- c(dvals, row$degree)
    }
  }
  
  summary$mean_participation[i] <- mean(pvals, na.rm = TRUE)
  summary$max_participation[i]  <- max(pvals, na.rm = TRUE)
  summary$mean_betweenness[i]   <- mean(bvals, na.rm = TRUE)
  summary$max_betweenness[i]    <- max(bvals, na.rm = TRUE)
  summary$mean_degree[i]        <- mean(dvals, na.rm = TRUE)
}

summary <- summary[order(-summary$mean_participation, -summary$mean_betweenness), ]

out_sum <- file.path(output_dir, "BridgingGenes_SummaryAcrossDays.csv")
write.csv(summary, out_sum, row.names = FALSE)
cat("  Wrote: ", basename(out_sum), "\n", sep = "")

# ============================================================
# NEW: Trajectory plot (top genes by mean participation)
# ============================================================
cat("\nBuilding participation trajectory plot for top genes by mean participation\n")

traj_png <- file.path(output_dir, "BridgingGenes_Trajectories_TopByMeanParticipation.png")
traj_plotly_rds <- NULL

# Pick top genes by mean_participation (computed above)
top_n <- min(top_n_for_trajectory_plot, nrow(summary))
top_genes <- summary$Gene[seq_len(top_n)]

traj_df <- integrated_long[integrated_long$Gene %in% top_genes, ]
traj_df$Gene <- factor(traj_df$Gene, levels = top_genes)  # keep ordering
traj_df <- traj_df[order(traj_df$Gene, traj_df$Day), ]

# Static plot using ggplot2 if available, otherwise base
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_traj <- ggplot(traj_df, aes(x = Day, y = participation, group = Gene, color = Gene)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    geom_point(size = 1.6, alpha = 0.9) +
    labs(
      title = paste0("Participation trajectories (top ", top_n, " genes by mean participation)"),
      x = "Day",
      y = "Participation coefficient"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  png(traj_png, width = trajectory_plot_png_width, height = trajectory_plot_png_height, res = trajectory_plot_png_res)
  print(p_traj)
  dev.off()
  
  if (isTRUE(enable_plotly)) {
    traj_plotly_rds <- file.path(output_dir, "BridgingGenes_Trajectories_TopByMeanParticipation_plotly.rds")
    try({
      w <- ggplotly(p_traj)
      saveRDS(w, traj_plotly_rds)
    }, silent = TRUE)
  }
} else {
  # Base fallback (no colors specified)
  png(traj_png, width = trajectory_plot_png_width, height = trajectory_plot_png_height, res = trajectory_plot_png_res)
  par(mar = c(4.5, 4.5, 2.8, 1.2))
  plot(NA, xlim = range(traj_df$Day), ylim = range(traj_df$participation, na.rm = TRUE),
       xlab = "Day", ylab = "Participation coefficient",
       main = paste0("Participation trajectories (top ", top_n, " genes by mean participation)"))
  for (g in top_genes) {
    dd <- traj_df[traj_df$Gene == g, ]
    lines(dd$Day, dd$participation, lwd = 1)
    points(dd$Day, dd$participation, pch = 16, cex = 0.6)
  }
  dev.off()
}

cat("  Wrote: ", basename(traj_png), "\n", sep = "")
if (!is.null(traj_plotly_rds)) cat("  Wrote: ", basename(traj_plotly_rds), "\n", sep = "")

# ============================================================
# 2×2 quadrant classification + schematic figure
# ============================================================
cat("\nBuilding 2×2 role quadrants (Participation vs Betweenness)\n")

bw_nonzero <- summary$mean_betweenness[is.finite(summary$mean_betweenness) & summary$mean_betweenness > 0]

out_roles <- NULL
out_png <- NULL
out_pdf <- NULL

if (length(bw_nonzero) < 10) {
  warning("Too few nonzero mean_betweenness values to compute robust quantile cutoff. Skipping 2×2 schematic.")
} else {
  
  bw_cut <- as.numeric(quantile(bw_nonzero, probs = betweenness_quantile, names = FALSE, type = 7))
  
  summary$P_class <- ifelse(summary$mean_participation > 0, "High participation", "Low participation")
  summary$B_class <- ifelse(summary$mean_betweenness >= bw_cut, "High betweenness", "Low betweenness")
  
  summary$Role <- NA_character_
  summary$Role[summary$P_class == "Low participation"  & summary$B_class == "Low betweenness"]  <- "Within-module genes"
  summary$Role[summary$P_class == "Low participation"  & summary$B_class == "High betweenness"] <- "Intra-module bottlenecks"
  summary$Role[summary$P_class == "High participation" & summary$B_class == "Low betweenness"]  <- "Module connectors"
  summary$Role[summary$P_class == "High participation" & summary$B_class == "High betweenness"] <- "System-level integrators"
  
  out_roles <- file.path(output_dir, "BridgingGenes_RoleQuadrants.csv")
  write.csv(summary, out_roles, row.names = FALSE)
  cat("  Wrote: ", basename(out_roles), "\n", sep = "")
  
  pick_exemplars <- function(dfq, n) {
    if (nrow(dfq) == 0) return(character())
    dfq <- dfq[order(-dfq$mean_betweenness, -dfq$mean_participation, -dfq$mean_degree), , drop = FALSE]
    head(dfq$Gene, n)
  }
  
  q1 <- summary[summary$Role == "Within-module genes", , drop = FALSE]
  q2 <- summary[summary$Role == "Intra-module bottlenecks", , drop = FALSE]
  q3 <- summary[summary$Role == "Module connectors", , drop = FALSE]
  q4 <- summary[summary$Role == "System-level integrators", , drop = FALSE]
  
  ex1 <- pick_exemplars(q1, n_exemplars_per_quadrant)
  ex2 <- pick_exemplars(q2, n_exemplars_per_quadrant)
  ex3 <- pick_exemplars(q3, n_exemplars_per_quadrant)
  ex4 <- pick_exemplars(q4, n_exemplars_per_quadrant)
  
  out_png <- file.path(output_dir, "BridgingGenes_RoleQuadrants_2x2.png")
  out_pdf <- file.path(output_dir, "BridgingGenes_RoleQuadrants_2x2.pdf")
  
  plot_2x2(
    device_fun = function(path, width, height) png(path, width = 1600, height = 1000, res = 150),
    out_path = out_png,
    ex1 = ex1, ex2 = ex2, ex3 = ex3, ex4 = ex4,
    bw_cut = bw_cut, betweenness_quantile = betweenness_quantile
  )
  
  plot_2x2(
    device_fun = function(path, width, height) pdf(path, width = 10, height = 7),
    out_path = out_pdf,
    ex1 = ex1, ex2 = ex2, ex3 = ex3, ex4 = ex4,
    bw_cut = bw_cut, betweenness_quantile = betweenness_quantile
  )
  
  cat("  Wrote: ", basename(out_png), "\n", sep = "")
  cat("  Wrote: ", basename(out_pdf), "\n", sep = "")
}

# ============================================================
# Plot objects for report (static + optional plotly)
# ============================================================
# Always-on static diagnostic: excluded gene counts by day
excluded_counts_png <- file.path(output_dir, "ExcludedGenes_CountsByDay.png")
try({
  excl_counts <- NULL
  if (!is.null(excluded_csv) && file.exists(excluded_csv)) {
    excl_df <- read.csv(excluded_csv, stringsAsFactors = FALSE)
    excl_counts <- aggregate(Gene ~ Day, excl_df, length)
    colnames(excl_counts) <- c("Day", "ExcludedGenes_n")
  } else {
    excl_counts <- data.frame(Day = unique(days), ExcludedGenes_n = 0)
    excl_counts <- excl_counts[order(excl_counts$Day), , drop = FALSE]
  }
  
  png(excluded_counts_png, width = 1400, height = 900, res = 150)
  par(mar = c(4.5, 4.5, 2.5, 1.2))
  plot(excl_counts$Day, excl_counts$ExcludedGenes_n, type = "b",
       xlab = "Day", ylab = "Excluded genes (count)",
       main = "Excluded genes per day (non-finite correlation QC)")
  dev.off()
}, silent = TRUE)

# Optional plotly artifact saved as RDS (for QMD flat-folder logic)
plotly_rds <- NULL
if (enable_plotly) {
  plotly_rds <- file.path(output_dir, "ExcludedGenes_CountsByDay_plotly.rds")
  try({
    if (exists("excl_counts")) {
      p <- ggplot(excl_counts, aes(x = Day, y = ExcludedGenes_n)) +
        geom_line() + geom_point() +
        labs(title = "Excluded genes per day (QC)", x = "Day", y = "Excluded genes (count)")
      w <- ggplotly(p)
      saveRDS(w, plotly_rds)
    }
  }, silent = TRUE)
}

# ============================================================
# Manifest + inventory helpers
# ============================================================
get_deps <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    ver <- NA_character_
    if (requireNamespace(p, quietly = TRUE)) {
      ver <- as.character(utils::packageVersion(p))
    }
    data.frame(package = p, version = ver, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

build_file_inventory <- function(dir_path, out_csv) {
  f <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) {
    write.csv(
      data.frame(path = character(), filename = character(), bytes = numeric(), modified = character()),
      out_csv, row.names = FALSE
    )
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

# ============================================================
# Runtime tracking end (capture before render; mandatory if enabled)
# ============================================================
end_time <- Sys.time()
runtime_seconds <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

# ============================================================
# Build manifest (pre-render)
# ============================================================
deps_pkgs <- c("igraph", "jsonlite", "rmarkdown", "knitr", "tools")
if (enable_plotly) deps_pkgs <- c(deps_pkgs, "ggplot2", "plotly", "htmltools")
deps_df <- get_deps(unique(deps_pkgs))

manifest <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script = list(
    name = script_name,
    path = script_path,
    full_path = if (is.na(script_full)) NA_character_ else script_full
  ),
  input = list(
    correlation_matrix_dir = mat_dir_norm,
    correlation_matrix_dir_name = mat_dir_name,
    filename_pattern = pat,
    n_files = length(files_full),
    
    # Full paths (authoritative)
    files_full = files_full,
    
    # Filenames only (human-readable)
    files_basename = files_base,
    
    # Map day -> filename (basename)
    files_by_day = files_by_day
  ),
  parameters = list(
    analysis_name = analysis_name,
    edge_density = edge_density,
    use_weighted_paths = use_weighted_paths,
    na_frac_thr = na_frac_thr,
    betweenness_quantile = betweenness_quantile,
    n_exemplars_per_quadrant = n_exemplars_per_quadrant,
    top_n_for_trajectory_plot = top_n_for_trajectory_plot,
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
writeLines(
  jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
  con = manifest_json_path
)

manifest_files_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")
build_file_inventory(output_dir, manifest_files_csv_path)

# ============================================================
# Write QMD (flat-folder logic; no paths; embed-resources: true)
# ============================================================
report_qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
report_html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))

# Extract this script header block (best-effort). If not readable, leave placeholder.
script_header_text <- NULL
try({
  if (!is.na(script_full) && file.exists(script_full)) {
    txt <- readLines(script_full, warn = FALSE)
    keep <- txt[1:min(length(txt), 140)]
    script_header_text <- paste(keep, collapse = "\n")
  }
}, silent = TRUE)
if (is.null(script_header_text)) {
  script_header_text <- "<Script header not available in this execution mode.>"
}

qmd_lines <- c(
  "---",
  paste0('title: "', script_name, " report\""),
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
  paste0('  script_name: "', script_name, '"'),
  paste0('  script_path: "', script_path, '"'),
  paste0('  script_full: "', ifelse(is.na(script_full), "NA", script_full), '"'),
  paste0('  output_dir: "', normalizePath(output_dir, winslash = "/", mustWork = FALSE), '"'),
  paste0('  enable_plotly: ', ifelse(enable_plotly, "true", "false")),
  paste0('  enable_runtime_tracking: ', ifelse(enable_runtime_tracking, "true", "false")),
  paste0('  runtime_seconds: ', ifelse(enable_runtime_tracking, as.character(runtime_seconds), "null")),
  "---",
  "",
  "## Script header",
  "```text",
  script_header_text,
  "```",
  "",
  "## Metadata",
  "",
  "- **Run ID:** `r params$run_id`",
  "- **Run timestamp:** `r params$run_timestamp`",
  "- **Script name:** `r params$script_name`",
  "- **Script path:** `r params$script_path`",
  "- **Script full path:** `r params$script_full`",
  "- **Output directory:** `r params$output_dir`",
  "",
  "## Inputs and parameters",
  "",
  "```{r}",
  "man <- jsonlite::fromJSON('Project_Manifest.json')",
  "cat('Input directory:', man$input$correlation_matrix_dir, '\\n')",
  "cat('Input directory name:', man$input$correlation_matrix_dir_name, '\\n')",
  "cat('Filename pattern:', man$input$filename_pattern, '\\n')",
  "cat('Number of files:', man$input$n_files, '\\n\\n')",
  "cat('Files (basename):\\n')",
  "print(man$input$files_basename)",
  "cat('\\nFiles by day (basename):\\n')",
  "print(man$input$files_by_day)",
  "cat('\\nParameters:\\n')",
  "print(man$parameters)",
  "```",
  "",
  "## Dependencies",
  "",
  "```{r}",
  "deps <- man$dependencies",
  "knitr::kable(deps)",
  "```",
  "",
  "## Generated files",
  "",
  "```{r}",
  "inv <- read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "## Integrated cross-day tables (requested)",
  "",
  "### Long format (Gene × Day)",
  "```{r}",
  "if (file.exists('BridgingGenes_Integrated_Long_AllDays.csv')) {",
  "  x <- read.csv('BridgingGenes_Integrated_Long_AllDays.csv', stringsAsFactors = FALSE)",
  "  knitr::kable(head(x, 20))",
  "} else {",
  "  cat('Missing: BridgingGenes_Integrated_Long_AllDays.csv')",
  "}",
  "```",
  "",
  "### Wide format (participation)",
  "```{r}",
  "if (file.exists('BridgingGenes_Integrated_Wide_Participation.csv')) {",
  "  x <- read.csv('BridgingGenes_Integrated_Wide_Participation.csv', stringsAsFactors = FALSE)",
  "  knitr::kable(head(x, 15))",
  "} else {",
  "  cat('Missing: BridgingGenes_Integrated_Wide_Participation.csv')",
  "}",
  "```",
  "",
  "### Wide format (betweenness)",
  "```{r}",
  "if (file.exists('BridgingGenes_Integrated_Wide_Betweenness.csv')) {",
  "  x <- read.csv('BridgingGenes_Integrated_Wide_Betweenness.csv', stringsAsFactors = FALSE)",
  "  knitr::kable(head(x, 15))",
  "} else {",
  "  cat('Missing: BridgingGenes_Integrated_Wide_Betweenness.csv')",
  "}",
  "```",
  "",
  "### Wide format (degree)",
  "```{r}",
  "if (file.exists('BridgingGenes_Integrated_Wide_Degree.csv')) {",
  "  x <- read.csv('BridgingGenes_Integrated_Wide_Degree.csv', stringsAsFactors = FALSE)",
  "  knitr::kable(head(x, 15))",
  "} else {",
  "  cat('Missing: BridgingGenes_Integrated_Wide_Degree.csv')",
  "}",
  "```",
  "",
  "## Cross-day summary",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_SummaryAcrossDays.csv')) {",
  "  s <- read.csv('BridgingGenes_SummaryAcrossDays.csv', stringsAsFactors = FALSE)",
  "  knitr::kable(head(s, 25))",
  "} else {",
  "  cat('Missing: BridgingGenes_SummaryAcrossDays.csv')",
  "}",
  "```",
  "",
  "## Trajectory plot (participation)",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_Trajectories_TopByMeanParticipation.png')) {",
  "  knitr::include_graphics('BridgingGenes_Trajectories_TopByMeanParticipation.png')",
  "} else {",
  "  cat('Missing: BridgingGenes_Trajectories_TopByMeanParticipation.png')",
  "}",
  "```",
  "",
  "## Analytical logic",
  "",
  "This analysis constructs fixed-density gene–gene networks from per-day correlation matrices by retaining the top |r| edges at a target density.",
  "",
  "**Edge selection:** Let $r_{ij}$ be the correlation between genes $i$ and $j$. We rank all off-diagonal pairs by $|r_{ij}|$ and retain the top $k$ edges where $k = \\lfloor density \\cdot \\binom{p}{2} \\rfloor$ for $p$ genes.",
  "",
  "**Community detection:** Louvain clustering is applied on the retained network (weights = |r|).",
  "",
  "**Participation coefficient:**",
  "$$P_i = 1 - \\sum_c \\left(\\frac{k_{i,c}}{k_i}\\right)^2$$",
  "where $k_i$ is the degree of node $i$ and $k_{i,c}$ counts edges from $i$ to community $c$.",
  "",
  "**Betweenness centrality:** computed on the unweighted graph (unless weighted paths are enabled in parameters).",
  "",
  "## Required static diagnostic figure",
  "",
  "```{r}",
  "if (file.exists('ExcludedGenes_CountsByDay.png')) {",
  "  knitr::include_graphics('ExcludedGenes_CountsByDay.png')",
  "} else {",
  "  cat('Static diagnostic plot not found: ExcludedGenes_CountsByDay.png')",
  "}",
  "```",
  "",
  "## 2×2 role schematic",
  "",
  "```{r}",
  "if (file.exists('BridgingGenes_RoleQuadrants_2x2.png')) {",
  "  knitr::include_graphics('BridgingGenes_RoleQuadrants_2x2.png')",
  "} else {",
  "  cat('Schematic not found: BridgingGenes_RoleQuadrants_2x2.png')",
  "}",
  "```",
  "",
  "## Interactive plots (Plotly)",
  "",
  "```{r results='asis'}",
  "if (isTRUE(params$enable_plotly)) {",
  "  rds_files <- list.files(params$output_dir, pattern = '\\\\_plotly\\\\.rds$', full.names = FALSE)",
  "  widgets <- list()",
  "  if (length(rds_files) == 0) {",
  "    widgets[[1]] <- htmltools::tags$p('No plotly RDS files found.')",
  "  } else {",
  "    for (f in rds_files) {",
  "      w <- tryCatch(readRDS(f), error = function(e) NULL)",
  "      if (is.null(w)) {",
  "        widgets[[length(widgets) + 1]] <- htmltools::tags$p(paste0('Could not read plot RDS: ', f))",
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
  "  cat('Plotly disabled (enable_plotly = FALSE).')",
  "}",
  "```",
  "",
  "## Runtime",
  "",
  "```{r}",
  "if (isTRUE(params$enable_runtime_tracking)) {",
  "  cat('Total runtime (seconds): ', params$runtime_seconds)",
  "} else {",
  "  cat('Runtime tracking disabled.')",
  "}",
  "```",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, con = report_qmd_path)

# ============================================================
# Render (final step): render into same output_dir, using params
# ============================================================
render_ok <- TRUE
render_err <- NULL

tryCatch({
  rmarkdown::render(
    input = report_qmd_path,
    output_file = basename(report_html_path),
    output_dir = output_dir,
    params = list(
      run_id = run_id,
      run_timestamp = manifest$run_timestamp,
      script_name = script_name,
      script_path = script_path,
      script_full = ifelse(is.na(script_full), "NA", script_full),
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      enable_plotly = enable_plotly,
      enable_runtime_tracking = enable_runtime_tracking,
      runtime_seconds = if (enable_runtime_tracking) runtime_seconds else NULL
    ),
    quiet = TRUE
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- conditionMessage(e)
})

# Refresh inventory AFTER render (mandatory)
build_file_inventory(output_dir, manifest_files_csv_path)

if (render_ok && file.exists(report_html_path)) {
  cat("\nHTML report created: ", normalizePath(report_html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("\nHTML report render FAILED.\n")
  cat("QMD written to: ", normalizePath(report_qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (!is.null(render_err)) cat("Render error: ", render_err, "\n", sep = "")
}

cat("\nDONE\n")
cat("Outputs root: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Run dir     : ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Run ID      : ", run_id, "\n", sep = "")
cat("Manifest    : ", basename(manifest_json_path), "\n", sep = "")
cat("Inventory   : ", basename(manifest_files_csv_path), "\n", sep = "")
cat("Summary CSV : ", basename(out_sum), "\n", sep = "")
cat("Long table  : ", basename(out_long), "\n", sep = "")
cat("Wide (P)    : ", basename(out_wide_part), "\n", sep = "")
cat("Wide (B)    : ", basename(out_wide_betw), "\n", sep = "")
cat("Wide (Deg)  : ", basename(out_wide_deg), "\n", sep = "")
cat("Traj plot   : ", basename(traj_png), "\n", sep = "")
if (!is.null(out_roles)) cat("Role CSV    : ", basename(out_roles), "\n", sep = "")
if (!is.null(out_png))   cat("Schematic   : ", basename(out_png), "\n", sep = "")
if (!is.null(out_pdf))   cat("Schematic   : ", basename(out_pdf), "\n", sep = "")
cat("\n")
