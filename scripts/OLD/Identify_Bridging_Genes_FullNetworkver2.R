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
#   - Draw 2×2 schematic:
#       Participation (low/high) vs Betweenness (low/high)
#
# ASSUMPTIONS
#   - Each CSV is a symmetric correlation matrix
#   - Row names = column names = gene symbols
#   - Filenames contain "DayXX" (e.g., Day14, Day16)
#
# OUTPUT (written into the selected directory)
#   - BridgingGenes_NetworkMetrics_DayXX.csv (per-day)
#   - BridgingGenes_ExcludedGenes_ByDay.csv (QC; genes excluded per day)
#   - BridgingGenes_SummaryAcrossDays.csv (cross-day summary)
#   - BridgingGenes_RoleQuadrants.csv (summary + quadrant labels)
#   - BridgingGenes_RoleQuadrants_2x2.png / .pdf (schematic)
# ============================================================

suppressPackageStartupMessages({
  library(igraph)
})

# ----------------------------
# Mac-friendly directory chooser
# ----------------------------
choose_directory <- function(prompt = "Select directory containing correlation matrix CSV files") {
  # Prefer RStudio's native picker when available (best on macOS; no Tcl/Tk).
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (rstudioapi::isAvailable()) {
      p <- tryCatch({
        rstudioapi::selectDirectory(caption = prompt)
      }, error = function(e) NULL)
      if (!is.null(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  
  # Fallback: Tcl/Tk chooser (can fail on some mac setups; only try if available).
  if (requireNamespace("tcltk", quietly = TRUE)) {
    p <- tryCatch({
      tcltk::tk_choose.dir(caption = prompt)
    }, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  # Last resort: manual paste.
  cat("\n", prompt, "\n", sep = "")
  cat("Paste full path (fallback):\n")
  p <- readline("Directory: ")
  if (!nzchar(p)) stop("No directory selected/provided.")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

# ----------------------------
# USER INPUT
# ----------------------------
cat("\nChoose directory containing correlation matrix CSV files...\n")
mat_dir <- choose_directory("Select the folder that contains correlation-matrix CSVs")

if (!dir.exists(mat_dir)) {
  stop("Directory does not exist: ", mat_dir)
}
cat("Using directory:\n  ", mat_dir, "\n", sep = "")

cat("\nOptional filename pattern [default = '.csv']:\n")
pat <- readline("Pattern: ")
if (!nzchar(pat)) pat <- "\\.csv$"

files <- list.files(mat_dir, pattern = pat, full.names = TRUE)
if (length(files) == 0) {
  stop("No files found matching pattern in ", mat_dir, " using pattern: ", pat)
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
  bad <- basename(files)[is.na(days)]
  stop("Could not infer Day from one or more filenames.\n",
       "Expected pattern like 'Day14'. Problem files:\n  - ",
       paste(bad, collapse = "\n  - "))
}

# Process in day order
ord_days <- order(days)
files <- files[ord_days]
days  <- days[ord_days]

# ----------------------------
# Parameters
# ----------------------------
edge_density <- 0.01
use_weighted_paths <- FALSE

# NA filtering rule (prevents single NA-gene from wiping the day)
na_frac_thr <- 0.05  # exclude genes with >5% non-finite correlations across partners

# 2x2 schematic thresholds
betweenness_quantile <- 0.85  # top 15% of nonzero mean_betweenness
n_exemplars_per_quadrant <- 6

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
# Participation coefficient
#   P_i = 1 - sum_c (k_i,c / k_i)^2
# ----------------------------
participation_coeff <- function(g, membership_vec) {
  vnames <- V(g)$name
  if (is.null(vnames) || any(!nzchar(vnames))) stop("Graph vertices must have non-empty names.")
  
  # Ensure membership is named by vertex name
  if (is.null(names(membership_vec))) {
    stop("membership_vec must be a named vector (names = vertex names).")
  }
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
    
    if (is.na(deg) || deg == 0) {
      P[vname] <- 0
      next
    }
    
    nbrs <- neighbors(g, vid)
    nbr_names <- V(g)[nbrs]$name
    
    if (length(nbr_names) == 0) {
      P[vname] <- 0
      next
    }
    
    comms <- membership_vec[nbr_names]
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
excluded_genes <- list()  # list of data.frames (Day, Gene, NA_fraction, Reason)

# ----------------------------
# Main loop
# ----------------------------
for (i in seq_along(files)) {
  
  f <- files[i]
  day <- days[i]
  
  cat("\nProcessing Day ", day, " (", basename(f), ")\n", sep = "")
  
  R <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  
  if (nrow(R) < 2 || ncol(R) < 2) {
    warning("Matrix too small in file: ", basename(f), " (", nrow(R), "x", ncol(R), "). Skipping.")
    next
  }
  
  if (!all(colnames(R) == rownames(R))) {
    stop("Matrix is not symmetric with matching row/column names: ", f)
  }
  
  suppressWarnings({ mode(R) <- "numeric" })
  
  # Ensure diagonal is sane for downstream use
  diag(R) <- 1
  
  # Handle non-finite correlations without wiping entire day
  if (any(!is.finite(R))) {
    
    na_frac <- rowMeans(!is.finite(R))
    bad_genes <- names(na_frac)[na_frac > na_frac_thr]
    
    cat("  Non-finite entries present.\n", sep = "")
    cat("  Genes with NA-fraction >", na_frac_thr, ": ", length(bad_genes), "\n", sep = "")
    if (length(bad_genes) > 0) {
      cat("  Example excluded: ", paste(head(bad_genes, 8), collapse = ", "), "\n", sep = "")
    }
    
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
  
  # Build graph
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
  
  # Avoid rowname warnings: unname() + row.names=NULL
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
  
  out_csv <- file.path(mat_dir, paste0("BridgingGenes_NetworkMetrics_Day", day, ".csv"))
  write.csv(df, out_csv, row.names = FALSE)
  cat("  Wrote: ", basename(out_csv), "\n", sep = "")
}

if (length(per_day_results) == 0) {
  stop("No per-day results were produced (all days skipped or failed).")
}

# ----------------------------
# Write excluded gene audit table
# ----------------------------
if (length(excluded_genes) > 0) {
  excl_df <- do.call(rbind, excluded_genes)
  excl_df <- excl_df[order(excl_df$Day, -excl_df$NA_fraction, excl_df$Gene), ]
  out_excl <- file.path(mat_dir, "BridgingGenes_ExcludedGenes_ByDay.csv")
  write.csv(excl_df, out_excl, row.names = FALSE)
  cat("  Wrote: ", basename(out_excl), "\n", sep = "")
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

summary <- summary[order(-summary$mean_participation, -summary$mean_betweenness), ]

out_sum <- file.path(mat_dir, "BridgingGenes_SummaryAcrossDays.csv")
write.csv(summary, out_sum, row.names = FALSE)
cat("  Wrote: ", basename(out_sum), "\n", sep = "")

# ============================================================
# 2×2 quadrant classification + schematic figure
# ============================================================
cat("\nBuilding 2×2 role quadrants (Participation vs Betweenness)\n")

bw_nonzero <- summary$mean_betweenness[is.finite(summary$mean_betweenness) & summary$mean_betweenness > 0]

if (length(bw_nonzero) < 10) {
  warning("Too few nonzero mean_betweenness values to compute robust quantile cutoff. Skipping 2×2 schematic.")
} else {
  
  bw_cut <- as.numeric(quantile(bw_nonzero, probs = betweenness_quantile, names = FALSE, type = 7))
  
  # Class assignments
  summary$P_class <- ifelse(summary$mean_participation > 0, "High participation", "Low participation")
  summary$B_class <- ifelse(summary$mean_betweenness >= bw_cut, "High betweenness", "Low betweenness")
  
  summary$Role <- NA_character_
  summary$Role[summary$P_class == "Low participation"  & summary$B_class == "Low betweenness"]  <- "Within-module genes"
  summary$Role[summary$P_class == "Low participation"  & summary$B_class == "High betweenness"] <- "Intra-module bottlenecks"
  summary$Role[summary$P_class == "High participation" & summary$B_class == "Low betweenness"]  <- "Module connectors"
  summary$Role[summary$P_class == "High participation" & summary$B_class == "High betweenness"] <- "System-level integrators"
  
  # Write role table
  out_roles <- file.path(mat_dir, "BridgingGenes_RoleQuadrants.csv")
  write.csv(summary, out_roles, row.names = FALSE)
  cat("  Wrote: ", basename(out_roles), "\n", sep = "")
  
  # Exemplar selection helper
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
  
  # Draw schematic (base R): all plotting inside function, no par(op) restoration
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
    
    panel_text(0.5, 1.5, "Intra-module bottlenecks", ex2)     # Low P / High B
    panel_text(1.5, 1.5, "System-level integrators", ex4)     # High P / High B
    panel_text(0.5, 0.5, "Within-module genes", ex1)          # Low P / Low B
    panel_text(1.5, 0.5, "Module connectors", ex3)            # High P / Low B
    
    mtext(
      paste0(
        "Cutoffs: participation > 0 = High; betweenness High = top ",
        round((1 - betweenness_quantile) * 100), "% of nonzero mean_betweenness (cut = ",
        signif(bw_cut, 3), ")."
      ),
      side = 1, line = 0.4, cex = 0.82
    )
  }
  
  out_png <- file.path(mat_dir, "BridgingGenes_RoleQuadrants_2x2.png")
  out_pdf <- file.path(mat_dir, "BridgingGenes_RoleQuadrants_2x2.pdf")
  
  # Use explicit PNG pixels to avoid "figure margins too large"
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

cat("\nDONE\n")
cat("Per-day files: BridgingGenes_NetworkMetrics_DayXX.csv\n")
cat("Excluded QC  : BridgingGenes_ExcludedGenes_ByDay.csv (if any)\n")
cat("Summary file : BridgingGenes_SummaryAcrossDays.csv\n")
cat("Role table   : BridgingGenes_RoleQuadrants.csv\n")
cat("Schematic    : BridgingGenes_RoleQuadrants_2x2.(png|pdf)\n")
