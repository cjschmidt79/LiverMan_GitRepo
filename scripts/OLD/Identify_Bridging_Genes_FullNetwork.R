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
#   - Summarize bridging behavior across days
#
# ASSUMPTIONS
#   - Each CSV is a symmetric correlation matrix
#   - Row names = column names = gene symbols
#   - Filenames contain "DayXX" (e.g., Day14, Day16)
#
# OUTPUT
#   - Per-day network metrics
#   - Cross-day gene-level summary
# ============================================================

suppressPackageStartupMessages({
  library(igraph)
})

# ----------------------------
# USER INPUT (minimal)
# ----------------------------
cat("\nEnter directory containing correlation matrix CSV files:\n")
cat("Example: outputs/TxCorrStructure_ByDay_*/\n")
mat_dir <- readline("Directory: ")

if (!dir.exists(mat_dir)) {
  stop("Directory does not exist: ", mat_dir)
}

cat("\nOptional filename pattern [default = '.csv']:\n")
pat <- readline("Pattern: ")
if (!nzchar(pat)) pat <- "\\.csv$"

files <- list.files(mat_dir, pattern = pat, full.names = TRUE)

if (length(files) == 0) {
  stop("No files found matching pattern in ", mat_dir)
}

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
  stop("Could not infer Day from one or more filenames. ",
       "Expected pattern like 'Day14'.")
}

names(days) <- basename(files)

# Optional: process in day order (keeps output tidy)
ord_days <- order(days)
files <- files[ord_days]
days  <- days[ord_days]

# ----------------------------
# Parameters
# ----------------------------
edge_density <- 0.01
use_weighted_paths <- FALSE

# ----------------------------
# Helper: build fixed-density network
# ----------------------------
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

# ----------------------------
# Participation coefficient (FIXED)
#   - Avoids v$name; iterates over vertex indices/names
# ----------------------------
participation_coeff <- function(g, membership) {
  vnames <- V(g)$name
  ki <- degree(g)
  if (is.null(names(ki))) names(ki) <- vnames
  
  P <- numeric(length(vnames))
  names(P) <- vnames
  
  for (vid in seq_along(vnames)) {
    vname <- vnames[vid]
    deg <- ki[vname]
    
    if (is.na(deg) || deg == 0) {
      P[vname] <- 0
      next
    }
    
    # neighbors() can take a vertex id (vid) or name; id is faster/safer here
    nbrs <- neighbors(g, vid)
    nbr_names <- as_ids(nbrs)
    
    if (length(nbr_names) == 0) {
      P[vname] <- 0
      next
    }
    
    comms <- membership[nbr_names]
    tab <- table(comms)
    frac <- as.numeric(tab) / deg
    P[vname] <- 1 - sum(frac^2)
  }
  
  P
}

# ----------------------------
# Containers
# ----------------------------
per_day_results <- list()
all_genes <- character()

# ----------------------------
# Main loop
# ----------------------------
for (i in seq_along(files)) {
  
  f <- files[i]
  day <- days[i]
  
  cat("\nProcessing Day ", day, " (", basename(f), ")\n", sep = "")
  
  R <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  
  if (!all(colnames(R) == rownames(R))) {
    stop("Matrix is not symmetric with matching row/column names: ", f)
  }
  
  g <- build_graph_fixed_density(R, edge_density)
  
  # Guard: if graph has too few vertices/edges, metrics may fail
  if (vcount(g) < 2 || ecount(g) < 1) {
    warning("Graph too small after thresholding for file: ", basename(f),
            " (vcount=", vcount(g), ", ecount=", ecount(g), "). Skipping day ", day, ".")
    next
  }
  
  comm <- cluster_louvain(g, weights = E(g)$weight)
  memb <- membership(comm)
  
  part <- participation_coeff(g, memb)
  
  if (use_weighted_paths) {
    # betweenness weights are interpreted as "cost"; use inverse-correlation as cost
    betw <- betweenness(g, weights = 1 / E(g)$weight, normalized = TRUE)
  } else {
    betw <- betweenness(g, normalized = TRUE)
  }
  
  deg <- degree(g)
  
  # Ensure named vectors align to vertex names
  vnames <- V(g)$name
  if (is.null(names(betw))) names(betw) <- vnames
  if (is.null(names(deg)))  names(deg)  <- vnames
  if (is.null(names(memb))) names(memb) <- vnames
  
  df <- data.frame(
    Gene = names(part),
    Day = day,
    participation = as.numeric(part),
    betweenness = betw[names(part)],
    degree = deg[names(part)],
    community = memb[names(part)],
    stringsAsFactors = FALSE
  )
  
  per_day_results[[as.character(day)]] <- df
  all_genes <- union(all_genes, df$Gene)
  
  out_csv <- file.path(
    mat_dir,
    paste0("BridgingGenes_NetworkMetrics_Day", day, ".csv")
  )
  write.csv(df, out_csv, row.names = FALSE)
  
  cat("  Wrote: ", basename(out_csv), "\n", sep = "")
}

if (length(per_day_results) == 0) {
  stop("No per-day results were produced (all days skipped or failed).")
}

# ----------------------------
# Cross-day summary
# ----------------------------
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

summary <- summary[order(-summary$mean_participation,
                         -summary$mean_betweenness), ]

out_sum <- file.path(mat_dir, "BridgingGenes_SummaryAcrossDays.csv")
write.csv(summary, out_sum, row.names = FALSE)

cat("\nDONE\n")
cat("Per-day files: BridgingGenes_NetworkMetrics_DayXX.csv\n")
cat("Summary file : BridgingGenes_SummaryAcrossDays.csv\n")
