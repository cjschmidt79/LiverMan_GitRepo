#!/usr/bin/env Rscript
# ============================================================
# Regulatory Landscape + DTW Functional Breakdown (Source-to-run)
#   - Mac/RStudio-friendly interactive pickers
#   - Builds integrated table (optional), summary tables, and plots
#   - Generates Quarto HTML report (minimal required sections)
#
# REQUIRED CONCEPTS:
#  - PCA_strength = max(|PC1_max|, |PC2_max|)
#  - Variance restructuring = Levene_F
#  - DTW clusters = Cluster
#  - Driver genes = high PCA_strength AND high variance AND high network centrality
#
# OUTPUTS:
#  - Liver_Regulatory_Summary_Table.csv (if rebuilt)
#  - DTW_cluster_functional_frequency_ANNOTATED_ONLY.csv
#  - DTW_cluster_functional_frequency_CORRECTED.csv (Unknown filled)
#  - CandidateDrivers_top10.csv, CandidateDrivers_top20.csv
#  - Regulatory_Landscape_Day14Boundary.png
#  - Regulatory_Landscape_DTWcolored.png
#  - Regulatory_Landscape_DTWcolored_NUMBERED.png + .pdf
#  - Regulatory_Landscape_PointID_Key.csv
#  - Quarto report HTML in output_dir
# ============================================================

# -----------------------------
# Small utilities (base R)
# -----------------------------
is_interactive_rstudio <- function() {
  interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
}

pick_file <- function(prompt = "Select a file") {
  cat("\n", prompt, "\n", sep = "")
  if (is_interactive_rstudio()) {
    p <- rstudioapi::selectFile(caption = prompt)
    if (is.null(p) || !nzchar(p)) stop("No file selected.")
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  } else {
    p <- file.choose()
    return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }
}

pick_dir <- function(prompt = "Select an output directory") {
  cat("\n", prompt, "\n", sep = "")
  if (is_interactive_rstudio()) {
    d <- rstudioapi::selectDirectory(caption = prompt)
    if (is.null(d) || !nzchar(d)) stop("No directory selected.")
    return(normalizePath(d, winslash = "/", mustWork = TRUE))
  } else {
    cat("RStudio not available; please type/paste the output directory path:\n")
    d <- readline("output_dir: ")
    if (!nzchar(d)) stop("No directory provided.")
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(d, winslash = "/", mustWork = TRUE))
  }
}

require_cols <- function(df, required, df_name = "data.frame") {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(paste0(df_name, " is missing required columns: ", paste(missing, collapse = ", ")))
  }
  TRUE
}

# -----------------------------
# Script identity capture (works when sourced in RStudio)
# -----------------------------
capture_script_identity <- function(verbose = TRUE) {
  script_full <- NA_character_
  
  # 1) RStudio active document path
  if (is_interactive_rstudio()) {
    ctx <- rstudioapi::getActiveDocumentContext()
    if (!is.null(ctx$path) && nzchar(ctx$path) && file.exists(ctx$path)) {
      script_full <- normalizePath(ctx$path, winslash = "/", mustWork = TRUE)
    }
  }
  
  # 2) Rscript invocation
  if (is.na(script_full)) {
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- "--file="
    hit <- grep(file_arg, ca, value = TRUE)
    if (length(hit) == 1) {
      cand <- sub(file_arg, "", hit, fixed = TRUE)
      if (file.exists(cand)) script_full <- normalizePath(cand, winslash = "/", mustWork = TRUE)
    }
  }
  
  # 3) Fallback (unknown)
  if (is.na(script_full)) script_full <- "UNKNOWN_SCRIPT_PATH"
  
  script_path <- if (script_full == "UNKNOWN_SCRIPT_PATH") "UNKNOWN_DIR" else dirname(script_full)
  script_name <- if (script_full == "UNKNOWN_SCRIPT_PATH") "UNKNOWN_SCRIPT" else basename(script_full)
  
  if (verbose) {
    cat("\n[script_name] ", script_name, "\n", sep = "")
    cat("[script_path] ", script_path, "\n", sep = "")
    cat("[script_full] ", script_full, "\n", sep = "")
  }
  
  list(script_name = script_name, script_path = script_path, script_full = script_full)
}

# -----------------------------
# Read CSV safely
# -----------------------------
read_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# -----------------------------
# Plot helper (base R)
# -----------------------------
save_png <- function(path, width = 2000, height = 1400, res = 250, expr) {
  png(filename = path, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  expr
}

save_pdf <- function(path, width = 10, height = 7, expr) {
  pdf(file = path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  expr
}

# ============================================================
# Step A — Capture absolute paths (once, early)
# ============================================================
meta <- capture_script_identity(verbose = TRUE)
script_name <- meta$script_name
script_path <- meta$script_path
script_full <- meta$script_full

timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Output directory (normalized, exists)
output_dir <- pick_dir("Select OUTPUT directory (all files + HTML will be written here)")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

# Ask if integrated summary table exists
cat("\nDo you already have the integrated summary table (Liver_Regulatory_Summary_Table.csv)?\n")
cat("Enter Y or N\n")
has_integrated <- toupper(trimws(readline("Have integrated summary table? [Y/N]: ")))
if (!has_integrated %in% c("Y", "N")) stop("Please enter Y or N.")

input_paths <- character(0)

# ============================================================
# Inputs (sequential prompts) + column requirements
# ============================================================
integrated_path <- NA_character_
dtw_assign_path <- NA_character_
pca_path <- NA_character_
levene_path <- NA_character_
func_map_path <- NA_character_
edges_path <- NA_character_

if (has_integrated == "Y") {
  integrated_path <- pick_file("Select integrated summary table CSV (Liver_Regulatory_Summary_Table.csv)")
  input_paths <- c(input_paths, integrated_path)
  
} else {
  cat("\nOK — we will REBUILD the integrated table.\n")
  cat("\nYou will be prompted sequentially for the required component tables.\n")
  
  cat("\n(1) DTW cluster assignments\n")
  cat("Required columns: Gene, Cluster\n")
  dtw_assign_path <- pick_file("Select DTWclust_BestSweep_Cluster_Assignments.csv")
  input_paths <- c(input_paths, dtw_assign_path)
  
  cat("\n(2) PCA loadings (gene-level)\n")
  cat("Required columns: Gene, PC1, PC2 (plus an interval column like day_paiir/day_pair is OK)\n")
  pca_path <- pick_file("Select gene_Pc1_PC2Loadings_all days.csv")
  input_paths <- c(input_paths, pca_path)
  
  cat("\n(3) Variance restructuring table (Levene)\n")
  cat("Required columns: Gene OR Metabolite (gene names), Levene_F, Levene_p, Perm_p\n")
  levene_path <- pick_file("Select MeanByDay_*_LevenePerm_ByFeature.csv")
  input_paths <- c(input_paths, levene_path)
  
  cat("\n(4) Functional mapping (pathway)\n")
  cat("Required columns: gene (or Gene), pathway\n")
  func_map_path <- pick_file("Select Gene2FunctiionMappiing4Edgeswap.csv")
  input_paths <- c(input_paths, func_map_path)
  
  cat("\n(5) Network metrics (optional but recommended)\n")
  cat("If provided, required columns: name (or Gene), Degree, BetweennessCentrality\n")
  cat("Enter to skip, or select file.\n")
  maybe <- trimws(readline("Provide Cytoscape network metrics? [Y/N]: "))
  maybe <- toupper(maybe)
  if (maybe == "Y") {
    edges_path <- pick_file("Select GeneEdges_cytoscape.csv")
    input_paths <- c(input_paths, edges_path)
  } else {
    edges_path <- NA_character_
  }
}

# Normalize input paths vector
input_paths <- normalizePath(input_paths, winslash = "/", mustWork = TRUE)

# ============================================================
# Load / Build integrated table
# ============================================================
integrated <- NULL

if (has_integrated == "Y") {
  integrated <- read_csv(integrated_path)
  
  # Normalize gene column name if needed
  gene_col <- NULL
  if ("Gene" %in% names(integrated)) gene_col <- "Gene"
  if (is.null(gene_col)) {
    # your earlier file had "Gene of Liver_Regulatory_Summary_Table"
    alt <- "Gene of Liver_Regulatory_Summary_Table"
    if (alt %in% names(integrated)) {
      names(integrated)[names(integrated) == alt] <- "Gene"
      gene_col <- "Gene"
    }
  }
  if (is.null(gene_col)) stop("Integrated table must contain a Gene column (Gene or 'Gene of Liver_Regulatory_Summary_Table').")
  
  # Minimal required columns for plotting/tables
  require_cols(integrated, c("Gene", "Cluster", "PC1_max", "PC2_max", "Levene_F", "Levene_p", "Perm_p"), "Integrated summary table")
  
  # Ensure pathway exists (if missing, fill Unknown)
  if (!("pathway" %in% names(integrated))) integrated$pathway <- "Unknown_or_Poorly_Characterized"
  
} else {
  # Build from components with validation
  dtw <- read_csv(dtw_assign_path); require_cols(dtw, c("Gene", "Cluster"), "DTW assignments")
  
  pca <- read_csv(pca_path); require_cols(pca, c("Gene", "PC1", "PC2"), "PCA loadings")
  pca_sum <- aggregate(cbind(PC1, PC2) ~ Gene, data = pca, FUN = function(x) c(mean = mean(x), max = max(x)))
  # pca_sum has matrix columns; expand
  PC1_mean <- pca_sum$PC1[, "mean"]; PC1_max <- pca_sum$PC1[, "max"]
  PC2_mean <- pca_sum$PC2[, "mean"]; PC2_max <- pca_sum$PC2[, "max"]
  pca_summary <- data.frame(Gene = pca_sum$Gene, PC1_mean = PC1_mean, PC1_max = PC1_max,
                            PC2_mean = PC2_mean, PC2_max = PC2_max, stringsAsFactors = FALSE)
  
  lev <- read_csv(levene_path)
  # accept Gene or Metabolite as gene field
  if (!("Gene" %in% names(lev))) {
    if ("Metabolite" %in% names(lev)) {
      names(lev)[names(lev) == "Metabolite"] <- "Gene"
    } else {
      stop("Levene table must contain Gene or Metabolite column.")
    }
  }
  require_cols(lev, c("Gene", "Levene_F", "Levene_p", "Perm_p"), "Levene table")
  lev_sel <- lev[, c("Gene", "Levene_F", "Levene_p", "Perm_p")]
  
  fmap <- read_csv(func_map_path)
  if (!("Gene" %in% names(fmap))) {
    if ("gene" %in% names(fmap)) names(fmap)[names(fmap) == "gene"] <- "Gene"
  }
  require_cols(fmap, c("Gene", "pathway"), "Functional mapping table")
  fmap_sel <- fmap[, c("Gene", "pathway")]
  
  # Optional network
  net_sel <- NULL
  if (!is.na(edges_path)) {
    net <- read_csv(edges_path)
    if (!("Gene" %in% names(net))) {
      if ("name" %in% names(net)) names(net)[names(net) == "name"] <- "Gene"
    }
    # Degree + Betweenness required if provided
    require_cols(net, c("Gene", "Degree", "BetweennessCentrality"), "Network metrics table")
    keep <- intersect(c("Gene", "Degree", "BetweennessCentrality", "ClosenessCentrality", "ClusteringCoefficient"), names(net))
    net_sel <- net[, keep]
  } else {
    net_sel <- data.frame(Gene = dtw$Gene, Degree = NA_real_, BetweennessCentrality = NA_real_, stringsAsFactors = FALSE)
  }
  
  # Merge
  integrated <- merge(dtw, pca_summary, by = "Gene", all.x = TRUE)
  integrated <- merge(integrated, fmap_sel, by = "Gene", all.x = TRUE)
  integrated <- merge(integrated, net_sel, by = "Gene", all.x = TRUE)
  integrated <- merge(integrated, lev_sel, by = "Gene", all.x = TRUE)
  
  # Fill pathway missing
  integrated$pathway[is.na(integrated$pathway) | !nzchar(integrated$pathway)] <- "Unknown_or_Poorly_Characterized"
  
  # Save integrated table
  integrated_out <- file.path(output_dir, "Liver_Regulatory_Summary_Table.csv")
  write.csv(integrated, integrated_out, row.names = FALSE)
  integrated_path <- integrated_out
}

# Ensure core derived variable
integrated$PCA_strength <- pmax(abs(integrated$PC1_max), abs(integrated$PC2_max), na.rm = TRUE)

# ============================================================
# Generate tables
# ============================================================

# A) Functional frequency (ANNOTATED ONLY)
annot_only <- integrated[integrated$pathway != "Unknown_or_Poorly_Characterized", ]
freq_annot <- as.data.frame(table(Cluster = annot_only$Cluster, FunctionalCategory = annot_only$pathway), stringsAsFactors = FALSE)
names(freq_annot)[names(freq_annot) == "Freq"] <- "GeneCount"
# Percent within cluster (annotated denominator)
tot_annot <- aggregate(GeneCount ~ Cluster, data = freq_annot, sum)
freq_annot <- merge(freq_annot, tot_annot, by = "Cluster", suffixes = c("", "_ClusterTotal"))
freq_annot$PercentWithinCluster <- round(100 * freq_annot$GeneCount / freq_annot$GeneCount_ClusterTotal, 2)
freq_annot <- freq_annot[order(freq_annot$Cluster, -freq_annot$GeneCount), ]

freq_annot_path <- file.path(output_dir, "DTW_cluster_functional_frequency_ANNOTATED_ONLY.csv")
write.csv(freq_annot, freq_annot_path, row.names = FALSE)

# B) Functional frequency (CORRECTED, includes Unknown)
freq_corr <- as.data.frame(table(Cluster = integrated$Cluster, FunctionalCategory = integrated$pathway), stringsAsFactors = FALSE)
names(freq_corr)[names(freq_corr) == "Freq"] <- "GeneCount"
cluster_sizes <- as.data.frame(table(Cluster = integrated$Cluster), stringsAsFactors = FALSE)
names(cluster_sizes)[names(cluster_sizes) == "Freq"] <- "ClusterSize"
freq_corr <- merge(freq_corr, cluster_sizes, by = "Cluster", all.x = TRUE)
freq_corr$PercentWithinCluster <- round(100 * freq_corr$GeneCount / freq_corr$ClusterSize, 2)
freq_corr <- freq_corr[order(freq_corr$Cluster, -freq_corr$GeneCount), ]

freq_corr_path <- file.path(output_dir, "DTW_cluster_functional_frequency_CORRECTED.csv")
write.csv(freq_corr, freq_corr_path, row.names = FALSE)

# C) Candidate driver lists (top 10% and top 20%)
calc_driver_list <- function(df, q = 0.9) {
  # High PCA
  pca_thr <- quantile(df$PCA_strength, probs = q, na.rm = TRUE)
  
  # High variance change
  lev_thr <- quantile(df$Levene_F, probs = q, na.rm = TRUE)
  
  # High network centrality (Degree or Betweenness)
  # If network missing, these may be NA; handle by requiring at least one non-NA
  deg_thr <- quantile(df$Degree, probs = q, na.rm = TRUE)
  bet_thr <- quantile(df$BetweennessCentrality, probs = q, na.rm = TRUE)
  
  high_pca <- df$PCA_strength >= pca_thr
  high_var <- df$Levene_F >= lev_thr
  
  high_net <- rep(FALSE, nrow(df))
  if (all(is.na(df$Degree)) && all(is.na(df$BetweennessCentrality))) {
    # No network info: relax this filter (but warn in report)
    high_net <- rep(TRUE, nrow(df))
  } else {
    high_net <- (df$Degree >= deg_thr) | (df$BetweennessCentrality >= bet_thr)
    high_net[is.na(high_net)] <- FALSE
  }
  
  out <- df[high_pca & high_var & high_net, ]
  out <- out[order(-out$Levene_F), ]
  out
}

drivers10 <- calc_driver_list(integrated, q = 0.9)
drivers20 <- calc_driver_list(integrated, q = 0.8)

drivers10_path <- file.path(output_dir, "CandidateDrivers_top10.csv")
drivers20_path <- file.path(output_dir, "CandidateDrivers_top20.csv")

write.csv(drivers10, drivers10_path, row.names = FALSE)
write.csv(drivers20, drivers20_path, row.names = FALSE)

# ============================================================
# Generate plots
# ============================================================

# Thresholds for "Day14 boundary" visualization (use top 20% like our plots)
pc_thresh <- as.numeric(quantile(integrated$PCA_strength, probs = 0.8, na.rm = TRUE))
lev_thresh <- as.numeric(quantile(integrated$Levene_F, probs = 0.8, na.rm = TRUE))

# Drivers used for labeling in plots (top20 intersection)
drivers_for_plot <- drivers20

# 1) Day14 boundary plot (uncolored)
p_day14 <- file.path(output_dir, "Regulatory_Landscape_Day14Boundary.png")
save_png(p_day14, expr = {
  plot(integrated$PCA_strength, integrated$Levene_F,
       xlab = "PCA Coordination Strength (max |PC1_max| or |PC2_max|)",
       ylab = "Variance Restructuring (Levene_F)",
       main = "Regulatory Landscape with Canalization Boundary (Day ~14)")
  abline(v = pc_thresh)
  abline(h = lev_thresh)
  points(drivers_for_plot$PCA_strength, drivers_for_plot$Levene_F, pch = 23)
  # label only drivers to keep readable
  if (nrow(drivers_for_plot) > 0) {
    text(drivers_for_plot$PCA_strength, drivers_for_plot$Levene_F, labels = drivers_for_plot$Gene,
         pos = 4, cex = 0.7)
  }
})

# 2) DTW-colored plot
p_dtw <- file.path(output_dir, "Regulatory_Landscape_DTWcolored.png")
save_png(p_dtw, expr = {
  cl <- sort(unique(integrated$Cluster))
  plot(integrated$PCA_strength, integrated$Levene_F, type = "n",
       xlab = "PCA Coordination Strength (max |PC1_max| or |PC2_max|)",
       ylab = "Variance Restructuring (Levene_F)",
       main = "Regulatory Landscape Colored by DTW Trajectory Cluster")
  # Use default palette without hardcoding specific colors
  pal <- seq_along(cl)
  for (i in seq_along(cl)) {
    idx <- integrated$Cluster == cl[i]
    points(integrated$PCA_strength[idx], integrated$Levene_F[idx], pch = 16)
  }
  abline(v = pc_thresh)
  abline(h = lev_thresh)
  points(drivers_for_plot$PCA_strength, drivers_for_plot$Levene_F, pch = 23, cex = 1.1)
  if (nrow(drivers_for_plot) > 0) {
    text(drivers_for_plot$PCA_strength, drivers_for_plot$Levene_F, labels = drivers_for_plot$Gene,
         pos = 4, cex = 0.7)
  }
  legend("topleft", legend = paste0("DTW Cluster ", cl),
         pch = 16, bty = "n")
})

# 3) Numbered plot + key
# Stable numbering by alphabetical Gene
num_df <- integrated[order(integrated$Gene), ]
num_df$PointID <- seq_len(nrow(num_df))

key_path <- file.path(output_dir, "Regulatory_Landscape_PointID_Key.csv")
key_cols <- intersect(c("PointID","Gene","pathway","Cluster","PC1_max","PC2_max","PCA_strength","Degree","BetweennessCentrality","Levene_F","Levene_p","Perm_p"),
                      names(num_df))
write.csv(num_df[, key_cols], key_path, row.names = FALSE)

p_num_png <- file.path(output_dir, "Regulatory_Landscape_DTWcolored_NUMBERED.png")
p_num_pdf <- file.path(output_dir, "Regulatory_Landscape_DTWcolored_NUMBERED.pdf")

make_numbered_plot <- function() {
  plot(num_df$PCA_strength, num_df$Levene_F, type = "n",
       xlab = "PCA Coordination Strength (max |PC1_max| or |PC2_max|)",
       ylab = "Variance Restructuring (Levene_F)",
       main = "Regulatory Landscape (Numbered Points) — see Key CSV")
  cl <- sort(unique(num_df$Cluster))
  for (i in seq_along(cl)) {
    idx <- num_df$Cluster == cl[i]
    points(num_df$PCA_strength[idx], num_df$Levene_F[idx], pch = 16)
  }
  abline(v = pc_thresh)
  abline(h = lev_thresh)
  # Annotate every point
  ok <- is.finite(num_df$PCA_strength) & is.finite(num_df$Levene_F)
  text(num_df$PCA_strength[ok], num_df$Levene_F[ok], labels = num_df$PointID[ok], cex = 0.55)
  # Highlight drivers (top20 list)
  points(drivers_for_plot$PCA_strength, drivers_for_plot$Levene_F, pch = 23, cex = 1.2)
  legend("topleft", legend = c(paste0("DTW Cluster ", cl), "High PCA + High variance"),
         pch = c(rep(16, length(cl)), 23), bty = "n")
}

save_png(p_num_png, width = 2600, height = 1800, res = 300, expr = { make_numbered_plot() })
save_pdf(p_num_pdf, width = 11, height = 8.5, expr = { make_numbered_plot() })

# ============================================================
# Build QMD as a single character vector (Step C/D)
# ============================================================

# Capture script header (verbatim comment block)
script_header <- ""
if (script_full != "UNKNOWN_SCRIPT_PATH" && file.exists(script_full)) {
  lines <- readLines(script_full, warn = FALSE)
  # take contiguous leading comments up to first non-comment/non-blank line
  keep <- character(0)
  for (ln in lines) {
    if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
      keep <- c(keep, ln)
    } else {
      break
    }
  }
  script_header <- paste(keep, collapse = "\n")
} else {
  script_header <- "# Script header unavailable (UNKNOWN_SCRIPT_PATH)"
}

# Quarto dependencies
deps <- c("base R", "quarto (for rendering report)")
if (is_interactive_rstudio()) deps <- c(deps, "rstudioapi (for file pickers)")

# Parameters used
params_list <- list(
  has_integrated = has_integrated,
  threshold_PCA_quantile = 0.8,
  threshold_Levene_quantile = 0.8,
  drivers10_quantile = 0.9,
  drivers20_quantile = 0.8
)

# Write QMD to output_dir
qmd_name <- paste0("Regulatory_Landscape_Report_", timestamp_tag, ".qmd")
qmd_path <- file.path(output_dir, qmd_name)

# Minimal required QMD sections ONLY (Step D)
qmd <- c(
  # 4445. YAML header
  "---",
  "title: \"Regulatory Landscape + DTW Functional Breakdown\"",
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 2",
  "execute:",
  "  echo: true",
  "  warning: false",
  "  message: false",
  "---",
  "",
  # 4446. Summary
  "## Summary",
  "This report summarizes gene-level regulatory structure by integrating:",
  "- DTW trajectory cluster membership",
  "- PCA coordination strength (max |PC1_max| or |PC2_max|)",
  "- variance restructuring across development (Levene_F)",
  "- (optional) network centrality metrics",
  "",
  # 4447. Script header (verbatim comment block)
  "## Script header (verbatim)",
  "```",
  script_header,
  "```",
  "",
  # 4448. Metadata (paths, timestamp, parameters)
  "## Metadata",
  paste0("- **timestamp:** ", timestamp),
  paste0("- **script_name:** ", script_name),
  paste0("- **script_path:** ", script_path),
  paste0("- **script_full:** ", script_full),
  paste0("- **output_dir:** ", output_dir),
  "",
  "### Input files (absolute paths)",
  "```",
  paste(input_paths, collapse = "\n"),
  "```",
  "",
  "### Parameters",
  "```",
  paste(capture.output(str(params_list)), collapse = "\n"),
  "```",
  "",
  # 4449. Dependencies
  "## Dependencies",
  paste0("- ", paste(deps, collapse = "\n- ")),
  "",
  # 4450. Analytical logic (Δ–Δ definition, equations)
  "## Analytical logic",
  "We define:",
  "",
  "- **PCA_strength** = max(|PC1_max|, |PC2_max|), capturing the strongest participation of a gene in either dominant coordination axis.",
  "- **Variance restructuring** is quantified by **Levene_F**, capturing developmental changes in inter-individual variance.",
  "- **Candidate regulatory drivers** are genes satisfying simultaneous high PCA_strength, high variance restructuring, and (if available) high network centrality.",
  "",
  "Thresholding used for visualization in the landscape plots:",
  "",
  "- PCA_strength threshold = top 20% quantile",
  "- Levene_F threshold = top 20% quantile",
  "",
  # 4451. Generated outputs (list.files(output_dir))
  "## Generated outputs",
  "```{r}",
  "list.files('.', full.names = TRUE)",
  "```",
  "",
  # 4452. Figures (auto-embed .png/.pdf)
  "## Figures",
  "### Regulatory landscape with Day-14 boundary",
  paste0("![](", basename(p_day14), ")"),
  "",
  "### Regulatory landscape colored by DTW cluster",
  paste0("![](", basename(p_dtw), ")"),
  "",
  "### Numbered regulatory landscape (PNG)",
  paste0("![](", basename(p_num_png), ")"),
  "",
  "### Numbered regulatory landscape (PDF for zooming)",
  paste0("- File: `", basename(p_num_pdf), "`"),
  "",
  # 4453. Session info
  "## Session info",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd, con = qmd_path)

# ============================================================
# Step B/E/F — Render with temporary WD switch (NO output_dir arg)
# ============================================================
if (!requireNamespace("quarto", quietly = TRUE)) {
  cat("\nNOTE: quarto R package not installed. Install with install.packages('quarto')\n")
  cat("QMD was written to: ", qmd_path, "\n", sep = "")
} else {
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  # Render using basename only, no output_dir parameter
  quarto::quarto_render(basename(qmd_path))
  
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  cat("\nReport written to:\n", html_path, "\n", sep = "")
  
  # Optional open in RStudio
  if (interactive()) {
    ans <- toupper(trimws(readline("Open HTML now? [Y/N]: ")))
    if (ans == "Y" && file.exists(html_path)) {
      browseURL(html_path)
    }
  }
}

cat("\nDONE.\n")
cat("Output directory:\n", output_dir, "\n", sep = "")