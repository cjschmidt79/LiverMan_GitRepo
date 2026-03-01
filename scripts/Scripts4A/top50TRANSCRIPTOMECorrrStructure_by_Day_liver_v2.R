#!/usr/bin/env Rscript

# ============================================================
# Correlation-Structure Comparisons (Transcriptome; FPKM; TIDY INPUT)
#
# PURPOSE
#   Apply the same correlation-structure analyses you ran for metabolomics to
#   transcriptome FPKM data, with transcriptomics-appropriate preprocessing:
#     - log2(FPKM + pseudocount)
#     - expression filtering (remove low/flat genes that destabilize correlations)
#     - within-day z-scoring per gene (across biological replicates)
#   Then compute:
#     A) Global similarity (RV, matrix correlation, Frobenius distance)
#     B) Edge turnover (top |Δr| gene-gene edges between day pairs)
#     C) Low-rank structure comparison (PC1 cosine similarity + angle)
#     D) Network modularity/integration at fixed edge density (igraph)
#
# INPUT (TIDY/LONG) REQUIRED columns (case-insensitive accepted):
#   Source, Day, SampleID, Gene, FPKM
#
# OUTPUT
#   outputs/<run_folder>/ with CSVs, optional heatmaps, summary plots,
#   a machine-readable manifest, and an auto-rendered HTML report.
#
# NOTES (important)
#   - Do NOT run on raw counts with this script.
#   - Filtering is essential for transcriptome correlation stability.
# ============================================================

# ----------------------------
# Script identity (robust; uses known filename fallback)
# ----------------------------
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

# User-provided current script filename (used if path cannot be detected)
known_script_filename <- "top50TRANSCRIPTOMECorrrStructure_by_Day_liver_v2.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  cat("\n[NOTE] Script full path could not be detected. Using:\n")
  cat("       script_name = ", script_name, "\n", sep = "")
  cat("       script_path = ", script_path, "\n", sep = "")
  cat("       script_full = NA\n")
  cat("       Tip: In RStudio, save the file and run via Source.\n")
  cat("            In terminal, run: Rscript ", known_script_filename, "\n\n", sep = "")
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

timestamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

# Always write into ./outputs
outputs_root <- file.path(getwd(), "outputs")
dir.create(outputs_root, showWarnings = FALSE, recursive = TRUE)

# Extract the script header comment block for embedding in the report
extract_header_text <- function(path, max_lines = 160) {
  if (is.na(path) || !file.exists(path)) {
    return(paste0(
      "Header text unavailable (script path not detected).\n",
      "If you want the header embedded, run via:\n",
      "  - RStudio: save file, then Source\n",
      "  - terminal: Rscript ", known_script_filename, "\n"
    ))
  }
  x <- readLines(path, warn = FALSE)
  x <- x[seq_len(min(length(x), max_lines))]
  keep <- character(0)
  for (ln in x) {
    if (grepl("^\\s*#!", ln) || grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
      keep <- c(keep, ln)
    } else {
      break
    }
  }
  paste(keep, collapse = "\n")
}

script_header_text <- extract_header_text(script_full)

# ----------------------------
# Packages
# ----------------------------
suppressWarnings(suppressMessages({
  req_pkgs <- c("igraph", "rmarkdown", "knitr", "jsonlite")
  for (p in req_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  library(igraph)
}))

# ----------------------------
# Helpers
# ----------------------------
pick_file_interactive <- function() {
  cat("\nSelect tidy transcriptome CSV file...\n")
  path <- tryCatch(file.choose(), error = function(e) NULL)
  if (is.null(path) || !nzchar(path)) {
    path <- readline("Enter full path to the tidy CSV: ")
  }
  if (!file.exists(path)) stop("File not found: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

safe_read_csv <- function(path) {
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read CSV: ", conditionMessage(e))
  )
}

# Case-insensitive mapping to required transcriptome columns
normalize_colnames_tx <- function(df) {
  nms <- names(df)
  lower <- tolower(nms)
  req <- c("source", "day", "sampleid", "gene", "fpkm")
  idx <- match(req, lower)
  if (any(is.na(idx))) {
    missing <- req[is.na(idx)]
    stop(
      "Missing required columns: ", paste(missing, collapse = ", "),
      "\nFound columns: ", paste(names(df), collapse = ", ")
    )
  }
  out <- df[, idx]
  names(out) <- c("Source", "Day", "SampleID", "Gene", "FPKM")
  out
}

as_numeric_day <- function(x) suppressWarnings(as.integer(as.character(x)))

# Build SampleID x Gene matrix (duplicates averaged)
wide_day_matrix <- function(df_day) {
  df_day$FPKM <- suppressWarnings(as.numeric(df_day$FPKM))
  df_day <- df_day[is.finite(df_day$FPKM), , drop = FALSE]
  if (nrow(df_day) == 0) stop("No finite FPKM values for this day after filtering.")
  
  key <- paste(df_day$SampleID, df_day$Gene, sep = "\t")
  if (any(duplicated(key))) {
    df_day <- aggregate(FPKM ~ SampleID + Gene, df_day, mean, na.rm = TRUE)
  }
  
  samples <- unique(df_day$SampleID)
  genes <- unique(df_day$Gene)
  
  mat <- matrix(NA_real_, nrow = length(samples), ncol = length(genes),
                dimnames = list(samples, genes))
  for (i in seq_len(nrow(df_day))) {
    mat[df_day$SampleID[i], df_day$Gene[i]] <- df_day$FPKM[i]
  }
  mat
}

# log2 transform
log2_fpkm <- function(mat, pseudocount = 0.1) log2(mat + pseudocount)

# Filter genes within a day using mean + SD on log2(FPKM+p)
filter_genes_day <- function(logmat, min_non_na = 3, min_mean = 0.0, min_sd = 0.25) {
  keep <- rep(TRUE, ncol(logmat))
  for (j in seq_len(ncol(logmat))) {
    x <- logmat[, j]
    nn <- sum(is.finite(x))
    if (nn < min_non_na) { keep[j] <- FALSE; next }
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(m) || !is.finite(s)) { keep[j] <- FALSE; next }
    if (m < min_mean) keep[j] <- FALSE
    if (s < min_sd) keep[j] <- FALSE
  }
  logmat[, keep, drop = FALSE]
}

# Z-score each gene within day across samples
zscore_cols <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]
    mu <- mean(x, na.rm = TRUE)
    s  <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      mat[, j] <- NA_real_
    } else {
      mat[, j] <- (x - mu) / s
    }
  }
  mat
}

correlation_matrix <- function(zmat) {
  if (ncol(zmat) < 10) stop("Too few genes after filtering (need >= 10 for stable structure).")
  R <- suppressWarnings(stats::cor(zmat, use = "pairwise.complete.obs", method = "pearson"))
  diag(R) <- 1
  R
}

upper_vec <- function(M) M[upper.tri(M, diag = FALSE)]

matrix_correlation <- function(A, B) {
  a <- upper_vec(A); b <- upper_vec(B)
  stats::cor(a, b, use = "complete.obs")
}

frobenius_distance <- function(A, B) sqrt(sum((A - B)^2, na.rm = TRUE))

rv_coefficient <- function(A, B) {
  # RV(A,B) = tr(A B) / sqrt(tr(A A) tr(B B)), for symmetric matrices.
  num <- sum(diag(A %*% B))
  den <- sqrt(sum(diag(A %*% A)) * sum(diag(B %*% B)))
  if (!is.finite(den) || den == 0) return(NA_real_)
  num / den
}

pc1_similarity <- function(A, B) {
  ea <- eigen(A, symmetric = TRUE)
  eb <- eigen(B, symmetric = TRUE)
  va <- ea$vectors[, 1]
  vb <- eb$vectors[, 1]
  cos <- sum(va * vb) / sqrt(sum(va^2) * sum(vb^2))
  cos_abs <- abs(cos)
  angle_deg <- acos(pmin(1, pmax(-1, cos_abs))) * (180 / pi)
  list(
    cos = cos,
    cos_abs = cos_abs,
    angle_deg = angle_deg,
    lambda1_over_sum_A = ea$values[1] / sum(ea$values),
    lambda1_over_sum_B = eb$values[1] / sum(eb$values)
  )
}

edge_turnover <- function(A, B, top_n = 50) {
  genes <- colnames(A)
  if (!identical(genes, colnames(B))) stop("Edge turnover requires matching gene sets.")
  iu <- which(upper.tri(A, diag = FALSE), arr.ind = TRUE)
  rA <- A[iu]
  rB <- B[iu]
  dr <- rB - rA
  out <- data.frame(
    Gene1 = genes[iu[, 1]],
    Gene2 = genes[iu[, 2]],
    r_A = rA,
    r_B = rB,
    delta_r = dr,
    abs_delta_r = abs(dr),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$abs_delta_r), , drop = FALSE]
  head(out, top_n)
}

graph_metrics_fixed_density <- function(R, density = 0.01) {
  genes <- colnames(R)
  p <- length(genes)
  if (p < 30) stop("Too few genes for network metrics (need >= 30).")
  
  n_possible <- p * (p - 1) / 2
  k <- max(1, floor(density * n_possible))
  
  iu <- which(upper.tri(R, diag = FALSE), arr.ind = TRUE)
  rvals <- R[iu]
  ord <- order(abs(rvals), decreasing = TRUE)
  take <- ord[seq_len(min(k, length(ord)))]
  edges <- iu[take, , drop = FALSE]
  w <- abs(rvals[take])
  
  g <- igraph::graph_from_edgelist(cbind(genes[edges[, 1]], genes[edges[, 2]]), directed = FALSE)
  igraph::E(g)$weight <- w
  
  comps <- igraph::components(g)
  giant_size <- if (length(comps$csize) > 0) max(comps$csize) else 0L
  giant_frac <- giant_size / igraph::vcount(g)
  
  clust <- tryCatch(igraph::transitivity(g, type = "globalundirected"), error = function(e) NA_real_)
  
  mod <- NA_real_
  n_comm <- NA_integer_
  if (igraph::ecount(g) > 0 && igraph::vcount(g) > 0) {
    comm <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$weight),
                     error = function(e) NULL)
    if (!is.null(comm)) {
      mod <- igraph::modularity(comm)
      n_comm <- length(unique(igraph::membership(comm)))
    }
  }
  
  data.frame(
    genes = p,
    edges_kept = igraph::ecount(g),
    density_target = density,
    clustering_global = clust,
    modularity_louvain = mod,
    n_communities = n_comm,
    giant_component_frac = giant_frac,
    stringsAsFactors = FALSE
  )
}

save_base_heatmap <- function(R, out_png, main = "") {
  png(out_png, width = 1400, height = 1200, res = 140)
  op <- par(mar = c(5, 5, 4, 2))
  image(
    1:ncol(R), 1:nrow(R),
    t(R[nrow(R):1, , drop = FALSE]),
    axes = FALSE,
    xlab = "", ylab = "",
    main = main
  )
  box()
  par(op)
  dev.off()
}

# A couple always-on summary plots (so the report can reliably include graphics)
plot_global_similarity <- function(sim_df, out_png) {
  png(out_png, width = 1300, height = 900, res = 140)
  op <- par(mar = c(10, 5, 4, 2))
  labs <- paste(sim_df$Day1, sim_df$Day2, sep = "→")
  barplot(sim_df$matrix_corr, names.arg = labs, las = 2,
          ylab = "Matrix correlation (upper triangle)", main = "Global similarity across day pairs")
  par(op)
  dev.off()
}

plot_network_metrics <- function(net_df, out_png) {
  png(out_png, width = 1300, height = 900, res = 140)
  op <- par(mar = c(5, 5, 4, 2))
  plot(net_df$Day, net_df$giant_component_frac, type = "b",
       xlab = "Day", ylab = "Giant component fraction",
       main = "Network integration at fixed edge density")
  par(op)
  dev.off()
}

# ----------------------------
# Main
# ----------------------------
cat("\n=== Day10/14/16 Correlation-Structure Analysis (Transcriptome; FPKM) ===\n")

in_path <- pick_file_interactive()
raw <- safe_read_csv(in_path)
dat <- normalize_colnames_tx(raw)

dat$Day <- as_numeric_day(dat$Day)
dat <- dat[!is.na(dat$Day), , drop = FALSE]

dat$FPKM <- suppressWarnings(as.numeric(dat$FPKM))
dat <- dat[is.finite(dat$FPKM), , drop = FALSE]

dat$Source   <- as.character(dat$Source)
dat$SampleID <- as.character(dat$SampleID)
dat$Gene     <- as.character(dat$Gene)

if (nrow(dat) == 0) stop("No usable rows after cleaning.")

# Choose Source
sources <- sort(unique(dat$Source))
cat("\nSources found:\n")
for (i in seq_along(sources)) cat(sprintf("  [%d] %s\n", i, sources[i]))
src_choice <- readline(sprintf("Choose Source index [default 1 = %s]: ", sources[1]))
if (!nzchar(src_choice)) src_choice <- "1"
src_i <- suppressWarnings(as.integer(src_choice))
if (!is.finite(src_i) || src_i < 1 || src_i > length(sources)) stop("Invalid Source index.")
SRC <- sources[src_i]
cat("Using Source: ", SRC, "\n", sep = "")

# Days
dflt_days <- c(10, 14, 16)
days_in <- sort(unique(dat$Day[dat$Source == SRC]))
cat("\nDays in this Source:\n  ", paste(days_in, collapse = ", "), "\n", sep = "")
day_str <- readline("Enter days to analyze (comma-separated) [default: 10,14,16]: ")
if (!nzchar(day_str)) {
  days <- dflt_days
} else {
  days <- suppressWarnings(as.integer(trimws(strsplit(day_str, ",")[[1]])))
  days <- days[is.finite(days)]
}
days <- sort(unique(days))
if (length(days) < 2) stop("Need at least 2 days for comparisons.")
missing_days <- setdiff(days, days_in)
if (length(missing_days) > 0) stop("Requested days not present for this Source: ", paste(missing_days, collapse = ", "))

# Transcriptome-specific preprocessing parameters
pc_str <- readline("log2 pseudocount added to FPKM [default 0.1]: ")
pseudocount <- if (!nzchar(pc_str)) 0.1 else suppressWarnings(as.numeric(pc_str))
if (!is.finite(pseudocount) || pseudocount < 0) pseudocount <- 0.1

mm_str <- readline("Gene filter: minimum mean log2(FPKM+p) [default 0.0]: ")
min_mean <- if (!nzchar(mm_str)) 0.0 else suppressWarnings(as.numeric(mm_str))
if (!is.finite(min_mean)) min_mean <- 0.0

msd_str <- readline("Gene filter: minimum SD log2(FPKM+p) within-day [default 0.25]: ")
min_sd <- if (!nzchar(msd_str)) 0.25 else suppressWarnings(as.numeric(msd_str))
if (!is.finite(min_sd) || min_sd < 0) min_sd <- 0.25

# Network density (lower default than metabolome)
dens_str <- readline("Network edge density (fraction) [default 0.01]: ")
density <- if (!nzchar(dens_str)) 0.01 else suppressWarnings(as.numeric(dens_str))
if (!is.finite(density) || density <= 0 || density >= 1) stop("Invalid density; use e.g. 0.005–0.02")

topn_str <- readline("Top-N edges for turnover tables [default 50]: ")
topN <- if (!nzchar(topn_str)) 50L else suppressWarnings(as.integer(topn_str))
if (!is.finite(topN) || topN < 5) topN <- 50L

make_heatmaps <- readline("Save correlation heatmaps? (y/n) [default n]: ")
save_heat <- (nzchar(make_heatmaps) && tolower(substr(make_heatmaps, 1, 1)) == "y")

# Output directory (ALWAYS under ./outputs)
run_id <- paste0("TxCorrStructure_ByDay_", gsub("\\s+", "_", SRC), "_", timestamp())
out_dir <- file.path(outputs_root, run_id)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("\nOutput directory:\n  ", out_dir, "\n", sep = "")

# Build correlation matrices per day
R_list <- list()
qc <- data.frame()

for (d in days) {
  cat("Processing Day ", d, "...\n", sep = "")
  df_d <- dat[dat$Source == SRC & dat$Day == d, c("SampleID","Gene","FPKM"), drop = FALSE]
  if (nrow(df_d) == 0) stop("No rows for Day ", d)
  
  mat <- wide_day_matrix(df_d)                         # samples x genes (FPKM)
  logmat <- log2_fpkm(mat, pseudocount = pseudocount)  # log2(FPKM + p)
  logmat <- filter_genes_day(
    logmat,
    min_non_na = 3,
    min_mean = min_mean,
    min_sd = min_sd
  )
  
  # QC
  qc <- rbind(qc, data.frame(
    Day = d,
    n_samples = nrow(logmat),
    n_genes_used = ncol(logmat),
    stringsAsFactors = FALSE
  ))
  
  if (ncol(logmat) < 10) {
    stop("After filtering at Day ", d, ", too few genes remain (", ncol(logmat),
         "). Relax min_mean/min_sd or verify data.")
  }
  
  zmat <- zscore_cols(logmat)
  R <- correlation_matrix(zmat)
  
  R_list[[as.character(d)]] <- R
  
  if (isTRUE(save_heat)) {
    save_base_heatmap(R,
                      file.path(out_dir, sprintf("CorrHeatmap_Day%d.png", d)),
                      main = sprintf("Gene–gene correlation (log2 FPKM, z-scored within day) - %s - Day %d", SRC, d)
    )
  }
}

write.csv(qc, file.path(out_dir, "QC_Samples_Genes_ByDay.csv"), row.names = FALSE)

# Harmonize to common gene intersection
common_genes <- Reduce(intersect, lapply(R_list, colnames))
if (length(common_genes) < 50) {
  stop("Too few common genes across selected days (", length(common_genes),
       "). Consider relaxing filters or using a pre-defined expressed-gene set.")
}
cat("\nCommon genes across days: ", length(common_genes), "\n", sep = "")

for (nm in names(R_list)) {
  R_list[[nm]] <- R_list[[nm]][common_genes, common_genes, drop = FALSE]
}

# Day pairs
pairs <- combn(as.character(days), 2, simplify = FALSE)

# ----------------------------
# A) Global similarity
# ----------------------------
sim_rows <- list()
for (pp in pairs) {
  d1 <- pp[1]; d2 <- pp[2]
  A <- R_list[[d1]]; B <- R_list[[d2]]
  sim_rows[[paste(d1, d2, sep = "_")]] <- data.frame(
    Day1 = as.integer(d1),
    Day2 = as.integer(d2),
    matrix_corr = matrix_correlation(A, B),
    frobenius = frobenius_distance(A, B),
    RV = rv_coefficient(A, B),
    stringsAsFactors = FALSE
  )
}
sim_df <- do.call(rbind, sim_rows)
write.csv(sim_df, file.path(out_dir, "A_GlobalSimilarity_ByDayPairs.csv"), row.names = FALSE)

# ----------------------------
# B) Edge turnover
# ----------------------------
for (pp in pairs) {
  d1 <- pp[1]; d2 <- pp[2]
  A <- R_list[[d1]]; B <- R_list[[d2]]
  et <- edge_turnover(A, B, top_n = topN)
  out <- file.path(out_dir, sprintf("B_EdgeTurnover_Top%d_Day%s_to_Day%s.csv", topN, d1, d2))
  write.csv(et, out, row.names = FALSE)
}

# ----------------------------
# C) PC1 similarity + PC1 loadings
# ----------------------------
pc_rows <- list()
for (pp in pairs) {
  d1 <- pp[1]; d2 <- pp[2]
  A <- R_list[[d1]]; B <- R_list[[d2]]
  pc <- pc1_similarity(A, B)
  pc_rows[[paste(d1, d2, sep = "_")]] <- data.frame(
    Day1 = as.integer(d1),
    Day2 = as.integer(d2),
    abs_cos_PC1 = pc$cos_abs,
    angle_deg = pc$angle_deg,
    lambda1_over_sum_day1 = pc$lambda1_over_sum_A,
    lambda1_over_sum_day2 = pc$lambda1_over_sum_B,
    stringsAsFactors = FALSE
  )
}
pc_df <- do.call(rbind, pc_rows)
write.csv(pc_df, file.path(out_dir, "C_PC1Similarity_ByDayPairs.csv"), row.names = FALSE)

# Save PC1 loadings per day
for (d in as.character(days)) {
  R <- R_list[[d]]
  ee <- eigen(R, symmetric = TRUE)
  pc1 <- ee$vectors[, 1]
  # Stabilize sign by forcing largest-magnitude loading positive
  j <- which.max(abs(pc1))
  if (pc1[j] < 0) pc1 <- -pc1
  pc1_df <- data.frame(Gene = common_genes, PC1_Loading = pc1, stringsAsFactors = FALSE)
  pc1_df <- pc1_df[order(-abs(pc1_df$PC1_Loading)), , drop = FALSE]
  write.csv(pc1_df, file.path(out_dir, sprintf("C_PC1Loadings_Day%s.csv", d)), row.names = FALSE)
}

# ----------------------------
# D) Network metrics at fixed density
# ----------------------------
net_rows <- list()
for (d in as.character(days)) {
  R <- R_list[[d]]
  met <- graph_metrics_fixed_density(R, density = density)
  met$Day <- as.integer(d)
  net_rows[[d]] <- met
}
net_df <- do.call(rbind, net_rows)
net_df <- net_df[, c("Day","genes","edges_kept","density_target","clustering_global",
                     "modularity_louvain","n_communities","giant_component_frac")]
write.csv(net_df, file.path(out_dir, "D_NetworkMetrics_FixedDensity.csv"), row.names = FALSE)

# ----------------------------
# Always-on summary plots for the report
# ----------------------------
plot_sim_png <- file.path(out_dir, "Plot_GlobalSimilarity_MatrixCorr.png")
plot_net_png <- file.path(out_dir, "Plot_Network_GiantComponentFrac.png")
plot_global_similarity(sim_df, plot_sim_png)
plot_network_metrics(net_df, plot_net_png)

# ----------------------------
# Project manifest (machine-readable)
# ----------------------------
deps_versions <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    ver <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
    list(package = p, version = ver)
  })
  out
}

generated_files <- function(dirp) {
  sort(normalizePath(list.files(dirp, full.names = TRUE, recursive = FALSE), winslash = "/", mustWork = FALSE))
}

manifest <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  script = list(
    name = script_name,
    path = script_path,
    full_path = script_full
  ),
  input = list(
    input_csv = in_path,
    source = SRC,
    days = as.integer(days)
  ),
  parameters = list(
    pseudocount = pseudocount,
    min_mean_log2 = min_mean,
    min_sd_log2 = min_sd,
    network_density = density,
    edge_turnover_topN = as.integer(topN),
    save_heatmaps = isTRUE(save_heat)
  ),
  outputs = list(
    outputs_root = outputs_root,
    output_dir = out_dir
  ),
  dependencies = deps_versions(c("igraph","rmarkdown","knitr","jsonlite")),
  key_output_files_expected = c(
    "QC_Samples_Genes_ByDay.csv",
    "A_GlobalSimilarity_ByDayPairs.csv",
    "B_EdgeTurnover_*.csv",
    "C_PC1Similarity_ByDayPairs.csv",
    "C_PC1Loadings_Day*.csv",
    "D_NetworkMetrics_FixedDensity.csv",
    "Plot_GlobalSimilarity_MatrixCorr.png",
    "Plot_Network_GiantComponentFrac.png"
  )
)

manifest_json_path <- file.path(out_dir, "Project_Manifest.json")
jsonlite::write_json(manifest, manifest_json_path, pretty = TRUE, auto_unbox = TRUE)

manifest_csv_path <- file.path(out_dir, "Project_Manifest_Files.csv")
mf <- data.frame(
  file = basename(generated_files(out_dir)),
  full_path = generated_files(out_dir),
  stringsAsFactors = FALSE
)
write.csv(mf, manifest_csv_path, row.names = FALSE)

# ----------------------------
# QMD + HTML report (auto-rendered)
# ----------------------------
qmd_path <- file.path(out_dir, paste0(script_name, "_Report.qmd"))
html_path <- file.path(out_dir, paste0(script_name, "_Report.html"))

heat_pngs <- sort(list.files(out_dir, pattern = "^CorrHeatmap_Day\\d+\\.png$", full.names = TRUE))
summary_pngs <- c(plot_sim_png, plot_net_png)

dep_tbl <- deps_versions(c("igraph","rmarkdown","knitr","jsonlite"))
dep_md <- paste0(
  "| Package | Version |\n|---|---|\n",
  paste(vapply(dep_tbl, function(x) paste0("| ", x$package, " | ", x$version, " |"), character(1)), collapse = "\n"),
  "\n"
)

qmd_lines <- c(
  "---",
  paste0('title: "', script_name, ' — Run Report"'),
  'format:',
  '  html:',
  '    toc: true',
  '    toc-depth: 3',
  'execute:',
  '  echo: false',
  '  warning: false',
  '  message: false',
  "params:",
  paste0("  out_dir: \"", gsub("\\\\", "/", out_dir), "\""),
  "---",
  "",
  "# Overview",
  "",
  "This report was generated automatically as the final step of the analysis script.",
  "",
  "## Script identity",
  "",
  paste0("- **Script name:** `", script_name, "`"),
  paste0("- **Script path:** `", script_path, "`"),
  paste0("- **Script full path:** `", ifelse(is.na(script_full), "NA", script_full), "`"),
  "",
  "## Run metadata",
  "",
  paste0("- **Run ID:** `", run_id, "`"),
  paste0("- **Run timestamp:** `", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "`"),
  paste0("- **Input CSV:** `", in_path, "`"),
  paste0("- **Source:** `", SRC, "`"),
  paste0("- **Days analyzed:** `", paste(days, collapse = ", "), "`"),
  paste0("- **Output directory:** `", out_dir, "`"),
  "",
  "# Embedded script header",
  "",
  "```",
  script_header_text,
  "```",
  "",
  "# Analytical logic and formulas",
  "",
  "## Preprocessing",
  "",
  "- Transform: **log2(FPKM + pseudocount)**.",
  "- Within each day, filter genes using mean and SD on log2-scale:",
  "  - Keep gene *g* if it has at least 3 finite observations, mean ≥ min_mean, and SD ≥ min_sd.",
  "- Z-score each gene within day across biological replicates:",
  "  - For gene *g* in day *d*:  \n    \\( z_{ig} = (x_{ig} - \\bar{x}_{g}) / s_g \\).",
  "",
  "## A) Global similarity between correlation matrices",
  "",
  "- **Matrix correlation**: Pearson correlation between upper-triangular entries of \\(R_A\\) and \\(R_B\\).",
  "- **Frobenius distance**: \\( \\|R_A - R_B\\|_F = \\sqrt{\\sum_{i,j} (R_{A,ij} - R_{B,ij})^2 } \\).",
  "- **RV coefficient**: \\( RV(R_A, R_B) = \\frac{\\mathrm{tr}(R_A R_B)}{\\sqrt{\\mathrm{tr}(R_A R_A)\\,\\mathrm{tr}(R_B R_B)}} \\).",
  "",
  "## B) Edge turnover",
  "",
  "- For each day pair, compute \\(\\Delta r_{ij} = r^{(B)}_{ij} - r^{(A)}_{ij}\\) for \\((i<j)\\). Rank by \\(|\\Delta r_{ij}|\\).",
  "",
  "## C) Low-rank structure comparison (PC1)",
  "",
  "- Compare leading eigenvectors using \\(|\\cos(\\theta)|\\) and \\(\\theta\\) (degrees). Report \\(\\lambda_1/\\sum\\lambda\\) per day.",
  "",
  "## D) Network metrics at fixed edge density",
  "",
  "- Keep top edges by \\(|r_{ij}|\\) at density \\(\\rho\\); compute clustering, Louvain modularity, community count, and giant component fraction.",
  "",
  "# Dependencies",
  "",
  dep_md,
  "",
  "# Key outputs (tables)",
  "",
  "```{r}",
  "sim_df <- read.csv(file.path(params$out_dir, 'A_GlobalSimilarity_ByDayPairs.csv'))",
  "net_df <- read.csv(file.path(params$out_dir, 'D_NetworkMetrics_FixedDensity.csv'))",
  "knitr::kable(sim_df, digits = 4, caption = 'A) Global similarity across day pairs')",
  "```",
  "",
  "```{r}",
  "knitr::kable(net_df, digits = 4, caption = 'D) Network metrics at fixed density')",
  "```",
  "",
  "# Plots",
  "",
  "## Summary plots (always generated)",
  "",
  "```{r}",
  "knitr::include_graphics(c(",
  paste0("  '", gsub("\\\\", "/", summary_pngs), "'", collapse = ",\n"),
  "))",
  "```",
  "",
  "## Correlation heatmaps (only if enabled)",
  "",
  "```{r}",
  "heat_pngs <- list.files(params$out_dir, pattern='^CorrHeatmap_Day\\\\d+\\\\.png$', full.names=TRUE)",
  "if (length(heat_pngs) == 0) {",
  "  cat('No heatmaps were generated (heatmap option disabled).')",
  "} else {",
  "  knitr::include_graphics(heat_pngs)",
  "}",
  "```",
  "",
  "# Interpretation guidance",
  "",
  "- High matrix correlation + low Frobenius distance: global structure is similar across days.",
  "- Large top-|Δr| edges: localized rewiring even when global similarity is moderate.",
  "- High \\(\\lambda_1/\\sum\\lambda\\): stronger low-rank (single-axis) organization.",
  "- Giant component fraction and community structure summarize integration vs partitioning at fixed density.",
  "",
  "# Reproducibility",
  "",
  "## Generated files",
  "",
  "```{r}",
  "mf <- read.csv(file.path(params$out_dir, 'Project_Manifest_Files.csv'))",
  "knitr::kable(mf, caption = 'Generated files in the run directory')",
  "```",
  "",
  "## Session information",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, qmd_path)

render_ok <- TRUE
render_msg <- NULL
tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_format = "html_document",
    output_file = basename(html_path),
    output_dir = out_dir,
    quiet = TRUE
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# Refresh file list after render attempt
mf2 <- data.frame(
  file = basename(generated_files(out_dir)),
  full_path = generated_files(out_dir),
  stringsAsFactors = FALSE
)
write.csv(mf2, manifest_csv_path, row.names = FALSE)

# ----------------------------
# Console summary
# ----------------------------
cat("\n=== DONE ===\n")
cat("Script identity:\n")
cat("  script_name: ", script_name, "\n", sep = "")
cat("  script_path: ", script_path, "\n", sep = "")
cat("  script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
cat("Input:\n  ", in_path, "\n", sep = "")
cat("Source:\n  ", SRC, "\n", sep = "")
cat("Days:\n  ", paste(days, collapse = ", "), "\n", sep = "")
cat("Common genes used:\n  ", length(common_genes), "\n", sep = "")
cat("\nPreprocessing:\n")
cat("  log2(FPKM + ", pseudocount, ")\n", sep = "")
cat("  filter: min_mean(log2) = ", min_mean, "; min_sd(log2) = ", min_sd, "\n", sep = "")
cat("  network density = ", density, "\n", sep = "")
cat("\nKey outputs:\n")
cat("  QC_Samples_Genes_ByDay.csv\n")
cat("  A_GlobalSimilarity_ByDayPairs.csv\n")
cat("  B_EdgeTurnover_*.csv\n")
cat("  C_PC1Similarity_ByDayPairs.csv\n")
cat("  C_PC1Loadings_Day*.csv\n")
cat("  D_NetworkMetrics_FixedDensity.csv\n")
cat("  Plot_GlobalSimilarity_MatrixCorr.png\n")
cat("  Plot_Network_GiantComponentFrac.png\n")
if (isTRUE(save_heat)) cat("  CorrHeatmap_Day*.png\n")
cat("  Project_Manifest.json\n")
cat("  Project_Manifest_Files.csv\n")
cat("  ", basename(qmd_path), "\n", sep = "")
if (render_ok && file.exists(html_path)) {
  cat("  ", basename(html_path), "\n", sep = "")
  cat("\nHTML report created:\n  ", html_path, "\n", sep = "")
} else {
  cat("\nHTML report render FAILED.\n")
  cat("QMD written to:\n  ", qmd_path, "\n", sep = "")
  if (!is.null(render_msg)) cat("Render error:\n  ", render_msg, "\n", sep = "")
}
cat("\nOutput directory:\n  ", out_dir, "\n", sep = "")

# End
