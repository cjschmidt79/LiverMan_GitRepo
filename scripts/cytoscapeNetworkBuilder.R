#!/usr/bin/env Rscript
# ============================================================
# Cytoscape Network Export Builder (Gene–Gene Correlations)
# ------------------------------------------------------------
# PURPOSE (manuscript-grade, reproducible):
#   Build Cytoscape-compatible edge and node tables from a long-format
#   gene–gene correlation table with per-slice (e.g., Day) edges and
#   optional statistical fields (FDR, bootstrap CI, stable-edge flags).
#
#   The script is intentionally "biologist-first":
#     - SOURCE-to-run in RStudio (interactive pickers)
#     - Minimal dependencies (NO dplyr)
#     - Explicit prompts and sensible defaults
#     - Standardized outputs under ./outputs/<run>_<timestamp>/
#     - Manifest + inventory for reproducibility
#     - Optional Quarto QMD/HTML report (safe contract)
#
# INPUTS (files; selected interactively):
#   REQUIRED:
#     1) Correlation edge table (long format): one row per gene pair per slice
#        Must include: source gene, target gene, slice (Day/Time/Condition), correlation value
#        Optional columns: p-values, FDR, CI_low/CI_high OR CI_excludes_zero flag, StableEdge flag
#
#   OPTIONAL (node attributes; merged by Gene ID):
#     2) GeneRole_Integrated_Table.csv   (Gene, Functional_Category, VarRatio..., Max_PC1_Window, etc.)
#     3) JointPCA_Agroup.csv             (Gene + columns like 4_8PC1 ... 18_20PC2) -> derive Axis_Class
#     4) DTW cluster table               (Gene, DTW_cluster)
#     5) Gene2Function mapping table     (Gene, Functional_Category) (optional override/augment)
#
# OUTPUTS (CSV; Cytoscape-ready):
#   edges/
#     - edges_all_combined.csv
#     - edges_filtered_combined.csv
#     - edges_slice_<slice>.csv                (filtered)
#     - edges_slice_<slice>_all.csv            (optional)
#   nodes/
#     - nodes_master.csv                       (one row per gene; merged annotations)
#     - unmatched_genes_in_edges.csv
#     - unmatched_genes_in_annotations.csv
#   metrics/ (optional, if chosen)
#     - node_metrics_by_slice.csv              (long: gene × slice)
#     - bridge_calls_by_slice.csv              (long: gene × slice; is_bridge)
#   report/ (optional)
#     - cytoscape_export_report.qmd
#     - cytoscape_export_report.html
#   provenance/
#     - manifest.csv
#     - inventory.csv
#
# NOTES ON INTERPRETATION:
#   - "Bridging" is computed per slice network (time-local), not globally.
#   - "Ever_Bridging" is a summary label; avoid plotting it as if simultaneous.
#
# ============================================================

# -----------------------------
# Libraries (minimal)
# -----------------------------
suppressWarnings(suppressMessages({
  ok_igraph <- requireNamespace("igraph", quietly = TRUE)
  ok_parallel <- requireNamespace("parallel", quietly = TRUE)
  ok_rstudioapi <- requireNamespace("rstudioapi", quietly = TRUE)
  ok_quarto <- requireNamespace("quarto", quietly = TRUE)
}))

if (!ok_igraph) {
  stop("Package 'igraph' is required. Install with: install.packages('igraph')", call. = FALSE)
}

# -----------------------------
# Utility: friendly file/directory pickers
# -----------------------------
pick_file <- function(prompt, must_exist = TRUE) {
  cat("\n", prompt, "\n", sep = "")
  path <- NA_character_
  
  if (interactive() && ok_rstudioapi) {
    path <- tryCatch(rstudioapi::selectFile(caption = prompt), error = function(e) NA_character_)
  }
  if (is.na(path) || !nzchar(path)) {
    path <- tryCatch(file.choose(), error = function(e) NA_character_)
  }
  if (must_exist && (is.na(path) || !file.exists(path))) {
    stop("File not found / not selected. Aborting.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = must_exist)
}

pick_dir <- function(prompt) {
  cat("\n", prompt, "\n", sep = "")
  path <- NA_character_
  
  if (interactive() && ok_rstudioapi) {
    path <- tryCatch(rstudioapi::selectDirectory(caption = prompt), error = function(e) NA_character_)
  }
  if (is.na(path) || !nzchar(path)) {
    path <- readline("Enter output parent directory (blank for current working dir): ")
    if (!nzchar(path)) path <- getwd()
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# -----------------------------
# Utility: safe prompts with defaults
# -----------------------------
prompt_choice <- function(title, choices, default_idx = 1) {
  cat("\n", title, "\n", sep = "")
  for (i in seq_along(choices)) cat(sprintf("  [%d] %s\n", i, choices[i]))
  ans <- readline(sprintf("Select [default %d]: ", default_idx))
  if (!nzchar(ans)) return(default_idx)
  idx <- suppressWarnings(as.integer(ans))
  if (is.na(idx) || idx < 1 || idx > length(choices)) {
    cat("Invalid choice; using default.\n")
    return(default_idx)
  }
  idx
}

prompt_yn <- function(question, default_yes = TRUE) {
  def <- if (default_yes) "Y" else "N"
  ans <- readline(sprintf("%s [Y/N, default %s]: ", question, def))
  if (!nzchar(ans)) return(default_yes)
  ans <- toupper(substr(trimws(ans), 1, 1))
  if (ans == "Y") return(TRUE)
  if (ans == "N") return(FALSE)
  default_yes
}

prompt_number <- function(question, default_val, min_val = -Inf, max_val = Inf) {
  ans <- readline(sprintf("%s [default %s]: ", question, as.character(default_val)))
  if (!nzchar(ans)) return(default_val)
  x <- suppressWarnings(as.numeric(ans))
  if (is.na(x) || x < min_val || x > max_val) {
    cat("Invalid value; using default.\n")
    return(default_val)
  }
  x
}

# -----------------------------
# Utility: script identity (best-effort)
# -----------------------------
get_script_identity <- function() {
  script_full <- NA_character_
  script_path <- NA_character_
  script_name <- NA_character_
  
  # Best-effort: works when running via Rscript or RStudio Source
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) > 0) {
    script_full <- sub(file_arg, "", hit[1])
    if (file.exists(script_full)) script_full <- normalizePath(script_full, winslash = "/", mustWork = TRUE)
  } else if (interactive()) {
    # RStudio: try to get active document path
    if (ok_rstudioapi) {
      script_full <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NA_character_)
      if (!is.na(script_full) && nzchar(script_full) && file.exists(script_full)) {
        script_full <- normalizePath(script_full, winslash = "/", mustWork = TRUE)
      } else {
        script_full <- NA_character_
      }
    }
  }
  
  if (!is.na(script_full)) {
    script_path <- dirname(script_full)
    script_name <- basename(script_full)
  } else {
    script_path <- getwd()
    script_name <- "Cytoscape_Export_Builder.R"
  }
  
  list(script_name = script_name, script_path = script_path, script_full = script_full)
}

# -----------------------------
# Utility: read CSV with basic checks
# -----------------------------
read_csv_checked <- function(path, required_cols = NULL) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  df <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) stop("Failed to read CSV: ", path, "\n", e$message, call. = FALSE))
  if (!is.null(required_cols)) {
    missing <- setdiff(required_cols, names(df))
    if (length(missing) > 0) {
      stop("Input missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }
  df
}

# -----------------------------
# Utility: column mapping for a dataframe
# -----------------------------
choose_column <- function(df, role, hint_patterns = NULL, allow_none = FALSE) {
  cols <- names(df)
  # Suggest default based on patterns
  default_idx <- 1
  if (!is.null(hint_patterns)) {
    hits <- integer()
    for (p in hint_patterns) {
      h <- grep(p, cols, ignore.case = TRUE)
      hits <- unique(c(hits, h))
    }
    if (length(hits) > 0) default_idx <- hits[1]
  }
  
  choices <- cols
  if (allow_none) choices <- c("<NONE>", choices)
  
  idx <- prompt_choice(
    title = sprintf("Select column for: %s", role),
    choices = choices,
    default_idx = if (allow_none) default_idx + 1 else default_idx
  )
  
  if (allow_none && idx == 1) return(NA_character_)
  if (allow_none) return(choices[idx])
  choices[idx]
}

# -----------------------------
# Utility: clean gene IDs consistently
# -----------------------------
clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

# -----------------------------
# Utility: make outputs folder
# -----------------------------
make_output_dir <- function(parent_dir, run_name) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe_run <- gsub("[^A-Za-z0-9_\\-]+", "_", run_name)
  out_dir <- file.path(parent_dir, "outputs", paste0(safe_run, "_", ts))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  subdirs <- c("edges", "nodes", "metrics", "report", "provenance", "figures", "tables")
  for (sd in subdirs) {
    d <- file.path(out_dir, sd)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(out_dir, winslash = "/", mustWork = TRUE)
}

# -----------------------------
# Utility: write manifest + inventory (CSV)
# -----------------------------
write_manifest_csv <- function(out_dir, manifest_list) {
  prov_dir <- file.path(out_dir, "provenance")
  if (!dir.exists(prov_dir)) dir.create(prov_dir, recursive = TRUE, showWarnings = FALSE)
  mpath <- file.path(prov_dir, "manifest.csv")
  
  # Flatten list -> 2-column key/value CSV
  keys <- names(manifest_list)
  vals <- vapply(manifest_list, function(v) {
    if (length(v) == 0) return("")
    paste(as.character(v), collapse = " | ")
  }, character(1))
  
  mdf <- data.frame(key = keys, value = vals, stringsAsFactors = FALSE)
  write.csv(mdf, mpath, row.names = FALSE)
  normalizePath(mpath, winslash = "/", mustWork = TRUE)
}

write_inventory_csv <- function(out_dir) {
  prov_dir <- file.path(out_dir, "provenance")
  ipath <- file.path(prov_dir, "inventory.csv")
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  files <- normalizePath(files, winslash = "/", mustWork = FALSE)
  rel <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", out_dir), "/?"), "", files)
  
  idf <- data.frame(
    rel_path = rel,
    abs_path = files,
    size_bytes = file.info(files)$size,
    stringsAsFactors = FALSE
  )
  write.csv(idf, ipath, row.names = FALSE)
  normalizePath(ipath, winslash = "/", mustWork = TRUE)
}

# -----------------------------
# Utility: CV + Axis_Class derivation from JointPCA loadings table
# -----------------------------
derive_axis_class <- function(jointpca_csv_path) {
  df <- read_csv_checked(jointpca_csv_path)
  
  # Detect gene column: first column often contains gene names in your file
  gene_col <- names(df)[1]
  # If there is a column literally named Gene, prefer it
  if ("Gene" %in% names(df)) gene_col <- "Gene"
  if ("gene" %in% names(df)) gene_col <- "gene"
  
  df$Gene <- clean_gene_id(df[[gene_col]])
  
  pc1_cols <- grep("PC1$", names(df), value = TRUE)
  pc2_cols <- grep("PC2$", names(df), value = TRUE)
  
  if (length(pc1_cols) < 2 || length(pc2_cols) < 2) {
    stop("JointPCA file does not appear to contain multiple *PC1 and *PC2 columns (e.g., 4_8PC1...).", call. = FALSE)
  }
  
  cv_fun <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) < 3) return(NA_real_)
    sd(x) / (mean(abs(x)) + 1e-6)
  }
  
  CV_PC1 <- apply(df[, pc1_cols, drop = FALSE], 1, cv_fun)
  CV_PC2 <- apply(df[, pc2_cols, drop = FALSE], 1, cv_fun)
  
  out <- data.frame(Gene = df$Gene, CV_PC1 = CV_PC1, CV_PC2 = CV_PC2, stringsAsFactors = FALSE)
  out <- out[!is.na(out$Gene), , drop = FALSE]
  
  # Quantile-based classes (stable and interpretable; avoids hardcoded magic numbers)
  q1 <- as.numeric(quantile(out$CV_PC1, probs = 0.33, na.rm = TRUE))
  q2 <- as.numeric(quantile(out$CV_PC1, probs = 0.67, na.rm = TRUE))
  q2_pc2 <- as.numeric(quantile(out$CV_PC2, probs = 0.67, na.rm = TRUE))
  
  axis_class <- rep("Intermediate", nrow(out))
  axis_class[out$CV_PC1 <= q1] <- "PC1-stable"
  axis_class[out$CV_PC1 >= q2] <- "PC1-volatile"
  axis_class[out$CV_PC2 >= q2_pc2 & out$CV_PC1 < q2] <- "PC2-dominant"
  
  out$Axis_Class <- axis_class
  out
}

# -----------------------------
# Utility: suggest cores (avoid thrashing)
# -----------------------------
suggest_cores <- function(n_slices, edges_per_slice) {
  avail <- if (ok_parallel) parallel::detectCores(logical = TRUE) else 1
  if (is.na(avail) || avail < 1) avail <- 1
  
  # Start with common-sense upper bound:
  #  - no more than slices
  #  - leave one core free
  base <- min(max(1, avail - 1), max(1, n_slices))
  
  # Data-driven scaling: if graphs are small, parallel overhead dominates
  med_e <- median(edges_per_slice, na.rm = TRUE)
  max_e <- max(edges_per_slice, na.rm = TRUE)
  
  # Conservative caps to reduce memory thrash on large graphs
  # (betweenness can be expensive; too many parallel graphs can blow RAM)
  if (med_e < 20000) {
    rec <- min(base, 4)
  } else if (med_e < 80000) {
    rec <- min(base, 8)
  } else if (med_e < 200000) {
    rec <- min(base, 12)
  } else {
    # Very large graphs: fewer parallel workers is often better
    rec <- min(base, 10)
  }
  
  list(avail = avail, recommended = max(1, rec), med_edges = med_e, max_edges = max_e)
}

# -----------------------------
# Core function
# -----------------------------
cytoscape_export_builder <- function() {
  
  # ---- Script identity + run name
  ident <- get_script_identity()
  cat("\n============================================================\n")
  cat("Cytoscape Network Export Builder\n")
  cat("Script: ", ident$script_name, "\n", sep = "")
  if (!is.na(ident$script_full)) cat("Full path: ", ident$script_full, "\n", sep = "")
  cat("============================================================\n")
  
  cat("\nThis script will create Cytoscape-ready CSVs:\n")
  cat("  - Edge lists (combined + per-slice)\n")
  cat("  - Node table with merged annotations\n")
  cat("  - Optional: per-slice betweenness (bridging) metrics\n")
  cat("  - Provenance manifest + inventory\n")
  cat("  - Optional: Quarto HTML report\n")
  
  run_name <- readline("\nEnter a short run name (e.g., liver_RV_networks): ")
  if (!nzchar(run_name)) run_name <- "cytoscape_export"
  
  parent_dir <- pick_dir("Select (or type) the parent directory where ./outputs/ will be created:")
  out_dir <- make_output_dir(parent_dir, run_name)
  cat("\nOutput directory:\n  ", out_dir, "\n", sep = "")
  
  # ---- Required: correlation edge table
  cat("\nREQUIRED FILE 1/?: Correlation edge table (long format)\n")
  cat("Expected: one row per gene pair per slice (e.g., Day)\n")
  cat("Must contain: source gene, target gene, slice variable, correlation/weight\n")
  edge_path <- pick_file("Select your correlation edge CSV (gene–gene correlations by slice):")
  edges_raw <- read_csv_checked(edge_path)
  
  # ---- Map edge columns (maximum flexibility)
  cat("\n--- Column mapping for the edge table ---\n")
  source_col <- choose_column(edges_raw, "SOURCE gene (node A)", hint_patterns = c("^Gene1$", "^Metab1$", "source", "from"))
  target_col <- choose_column(edges_raw, "TARGET gene (node B)", hint_patterns = c("^Gene2$", "^Metab2$", "target", "to"))
  slice_col  <- choose_column(edges_raw, "SLICE variable (e.g., Day/Time/Condition)", hint_patterns = c("^Day$", "time", "stage", "window", "slice", "group"))
  weight_col <- choose_column(edges_raw, "EDGE WEIGHT (correlation/association)", hint_patterns = c("Pearson", "Spearman", "cor", "rho", "weight"))
  
  # Optional statistical columns (allow NONE)
  fdr_col <- choose_column(edges_raw, "FDR column (optional; used for filtering)", hint_patterns = c("FDR", "BH", "qval"), allow_none = TRUE)
  ci_ex_col <- choose_column(edges_raw, "CI excludes zero flag (optional; 0/1 or TRUE/FALSE)", hint_patterns = c("CI_ExcludesZero", "excludes"), allow_none = TRUE)
  ci_low_col <- choose_column(edges_raw, "Bootstrap CI LOWER (optional)", hint_patterns = c("LowerCI", "CI_low", "Boot_Lower"), allow_none = TRUE)
  ci_high_col <- choose_column(edges_raw, "Bootstrap CI UPPER (optional)", hint_patterns = c("UpperCI", "CI_high", "Boot_Upper"), allow_none = TRUE)
  stable_col <- choose_column(edges_raw, "Stable edge flag (optional; TRUE/FALSE)", hint_patterns = c("StableEdge", "stable"), allow_none = TRUE)
  
  # ---- Bug-detection stanza: basic edge table sanity
  cat("\n--- Validating edge table ---\n")
  edges_raw[[source_col]] <- clean_gene_id(edges_raw[[source_col]])
  edges_raw[[target_col]] <- clean_gene_id(edges_raw[[target_col]])
  edges_raw[[slice_col]]  <- edges_raw[[slice_col]]
  
  # Convert weight to numeric
  w <- suppressWarnings(as.numeric(edges_raw[[weight_col]]))
  if (all(!is.finite(w))) stop("Edge weight column could not be parsed as numeric.", call. = FALSE)
  edges_raw[[weight_col]] <- w
  
  # Drop rows with missing essentials
  keep <- !is.na(edges_raw[[source_col]]) & !is.na(edges_raw[[target_col]]) &
    !is.na(edges_raw[[slice_col]]) & is.finite(edges_raw[[weight_col]])
  dropped <- sum(!keep)
  if (dropped > 0) cat("Dropping ", dropped, " rows with missing source/target/slice/weight.\n", sep = "")
  edges_raw <- edges_raw[keep, , drop = FALSE]
  
  if (nrow(edges_raw) < 10) stop("Too few edges after basic cleaning (<10). Check inputs/column mapping.", call. = FALSE)
  
  # Remove self-loops (optional but recommended)
  self_loops <- edges_raw[[source_col]] == edges_raw[[target_col]]
  if (any(self_loops, na.rm = TRUE)) {
    nsl <- sum(self_loops, na.rm = TRUE)
    cat("Removing ", nsl, " self-loops (source==target).\n", sep = "")
    edges_raw <- edges_raw[!self_loops, , drop = FALSE]
  }
  
  # ---- Decide filtering strategy (sensible defaults + explanation)
  cat("\n--- Edge filtering strategy ---\n")
  cat("Filtering controls network density. Dense networks can obscure structure and may be slow in Cytoscape.\n")
  cat("You can export both unfiltered and filtered edges (recommended).\n")
  
  has_stable <- !is.na(stable_col) && stable_col %in% names(edges_raw)
  has_fdr <- !is.na(fdr_col) && fdr_col %in% names(edges_raw)
  has_ci_flag <- !is.na(ci_ex_col) && ci_ex_col %in% names(edges_raw)
  has_ci_bounds <- (!is.na(ci_low_col) && !is.na(ci_high_col) &&
                      ci_low_col %in% names(edges_raw) && ci_high_col %in% names(edges_raw))
  
  # Default filter mode priority: StableEdge > CI flag/bounds > FDR > top-N > |r|
  filter_choices <- c(
    "No filtering (export all edges)",
    "Stable-edge flag (e.g., StableEdge == TRUE)",
    "Bootstrap CI excludes zero (CI flag or CI bounds)",
    "FDR cutoff (choose column + threshold)",
    "Absolute correlation threshold (|weight| >= cutoff)",
    "Top-N edges per slice (by |weight|)"
  )
  
  default_filter <- 6
  if (has_stable) default_filter <- 2
  else if (has_ci_flag || has_ci_bounds) default_filter <- 3
  else if (has_fdr) default_filter <- 4
  
  f_idx <- prompt_choice(
    title = "Choose how to filter edges (briefly: increases interpretability + Cytoscape performance):",
    choices = filter_choices,
    default_idx = default_filter
  )
  
  # Prepare working edge table with standardized columns
  edges <- data.frame(
    source = edges_raw[[source_col]],
    target = edges_raw[[target_col]],
    slice  = edges_raw[[slice_col]],
    weight = edges_raw[[weight_col]],
    stringsAsFactors = FALSE
  )
  edges$abs_weight <- abs(edges$weight)
  edges$sign <- ifelse(edges$weight >= 0, "pos", "neg")
  
  # Attach optional stats columns if present (for user reference in Cytoscape)
  add_if_present <- function(out_df, raw_df, colname, newname) {
    if (!is.na(colname) && colname %in% names(raw_df)) {
      out_df[[newname]] <- raw_df[[colname]]
    }
    out_df
  }
  edges <- add_if_present(edges, edges_raw, fdr_col, "FDR")
  edges <- add_if_present(edges, edges_raw, ci_ex_col, "CI_ExcludesZero")
  edges <- add_if_present(edges, edges_raw, ci_low_col, "CI_Low")
  edges <- add_if_present(edges, edges_raw, ci_high_col, "CI_High")
  edges <- add_if_present(edges, edges_raw, stable_col, "StableEdge")
  
  # ---- Build filtered edge set
  edges_f <- edges
  filter_desc <- "none"
  filter_params <- ""
  
  if (f_idx == 1) {
    # no filtering
    filter_desc <- "none"
  } else if (f_idx == 2) {
    if (!("StableEdge" %in% names(edges_f))) stop("StableEdge column not available.", call. = FALSE)
    cat("Using StableEdge == TRUE.\n")
    filter_desc <- "StableEdge==TRUE"
    edges_f <- edges_f[as.logical(edges_f$StableEdge) %in% TRUE, , drop = FALSE]
  } else if (f_idx == 3) {
    # CI excludes zero
    if ("CI_ExcludesZero" %in% names(edges_f)) {
      cat("Using CI_ExcludesZero == 1/TRUE.\n")
      filter_desc <- "CI_ExcludesZero"
      keep_ci <- edges_f$CI_ExcludesZero
      keep_ci <- if (is.logical(keep_ci)) keep_ci else (as.numeric(keep_ci) == 1)
      edges_f <- edges_f[keep_ci %in% TRUE, , drop = FALSE]
    } else if (all(c("CI_Low", "CI_High") %in% names(edges_f))) {
      cat("Using CI bounds: keep edges where CI_Low > 0 OR CI_High < 0.\n")
      filter_desc <- "CI_bounds_exclude_zero"
      lo <- suppressWarnings(as.numeric(edges_f$CI_Low))
      hi <- suppressWarnings(as.numeric(edges_f$CI_High))
      keep_ci <- is.finite(lo) & is.finite(hi) & (lo > 0 | hi < 0)
      edges_f <- edges_f[keep_ci, , drop = FALSE]
    } else {
      stop("No CI flag or CI bounds available.", call. = FALSE)
    }
  } else if (f_idx == 4) {
    if (!("FDR" %in% names(edges_f))) stop("FDR column not available.", call. = FALSE)
    cat("FDR column should be numeric in [0,1].\n")
    thr <- prompt_number("Enter FDR threshold (commonly 0.05 or 0.10)", default_val = 0.05, min_val = 0, max_val = 1)
    filter_desc <- "FDR<thr"
    filter_params <- paste0("thr=", thr)
    fdrv <- suppressWarnings(as.numeric(edges_f$FDR))
    keep_fdr <- is.finite(fdrv) & fdrv <= thr
    edges_f <- edges_f[keep_fdr, , drop = FALSE]
  } else if (f_idx == 5) {
    cat("This keeps edges with strong correlations regardless of p-values.\n")
    thr <- prompt_number("Enter absolute weight cutoff (e.g., 0.5, 0.7)", default_val = 0.6, min_val = 0, max_val = 1)
    filter_desc <- "|weight|>=thr"
    filter_params <- paste0("thr=", thr)
    edges_f <- edges_f[edges_f$abs_weight >= thr, , drop = FALSE]
  } else if (f_idx == 6) {
    cat("This keeps the strongest edges per slice. Best for Cytoscape performance.\n")
    # Default N based on dataset size (sensible and safe)
    slices <- unique(edges_f$slice)
    n_slices <- length(slices)
    # Heuristic default: 1000 edges/slice (adjustable)
    default_n <- 1000
    N <- as.integer(prompt_number("Enter N edges per slice (ranked by |weight|)", default_val = default_n, min_val = 10, max_val = 1e7))
    filter_desc <- "TopN_per_slice"
    filter_params <- paste0("N=", N)
    
    keep_rows <- logical(nrow(edges_f))
    for (s in slices) {
      idx <- which(edges_f$slice == s)
      if (length(idx) <= N) {
        keep_rows[idx] <- TRUE
      } else {
        ord <- order(edges_f$abs_weight[idx], decreasing = TRUE)
        keep_rows[idx[ord[seq_len(N)]]] <- TRUE
      }
    }
    edges_f <- edges_f[keep_rows, , drop = FALSE]
  }
  
  # ---- Bug-detection stanza: post-filter sanity
  if (nrow(edges_f) < 10) {
    stop("Too few edges after filtering (<10). Loosen filter or choose a different strategy.", call. = FALSE)
  }
  
  # ---- Write edges (combined + per-slice)
  edges_dir <- file.path(out_dir, "edges")
  all_combined_path <- file.path(edges_dir, "edges_all_combined.csv")
  filt_combined_path <- file.path(edges_dir, "edges_filtered_combined.csv")
  
  write.csv(edges, all_combined_path, row.names = FALSE)
  write.csv(edges_f, filt_combined_path, row.names = FALSE)
  
  # Per-slice exports
  slices <- sort(unique(edges$slice))
  for (s in slices) {
    safe_s <- gsub("[^A-Za-z0-9_\\-]+", "_", as.character(s))
    # filtered
    p_f <- file.path(edges_dir, paste0("edges_slice_", safe_s, ".csv"))
    write.csv(edges_f[edges_f$slice == s, , drop = FALSE], p_f, row.names = FALSE)
    # unfiltered (optional; can be huge)
    if (prompt_yn(paste0("Also export UNFILTERED edges for slice '", s, "'? (can be large)"), default_yes = FALSE)) {
      p_a <- file.path(edges_dir, paste0("edges_slice_", safe_s, "_all.csv"))
      write.csv(edges[edges$slice == s, , drop = FALSE], p_a, row.names = FALSE)
    }
  }
  
  cat("\nEdges written:\n")
  cat("  ", normalizePath(all_combined_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("  ", normalizePath(filt_combined_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  
  # ---- Build node table backbone from unique genes in FILTERED edges (Cytoscape practicality)
  nodes <- data.frame(id = sort(unique(c(edges_f$source, edges_f$target))), stringsAsFactors = FALSE)
  nodes$id <- clean_gene_id(nodes$id)
  
  # ---- Optional: merge node annotations
  cat("\n--- Optional node annotations ---\n")
  cat("If you provide annotation tables, they will be merged by Gene ID.\n")
  
  # 1) GeneRole integrated table
  gene_role_path <- NA_character_
  if (prompt_yn("Do you have a GeneRole_Integrated_Table CSV to merge (Functional_Category, VarRatio, Max_PC1_Window, etc.)?", TRUE)) {
    cat("\nOPTIONAL FILE: GeneRole_Integrated_Table.csv\n")
    cat("Expected: a 'Gene' column plus annotation columns.\n")
    gene_role_path <- pick_file("Select GeneRole_Integrated_Table CSV:")
    gr <- read_csv_checked(gene_role_path)
    gene_col <- if ("Gene" %in% names(gr)) "Gene" else names(gr)[1]
    gr$Gene <- clean_gene_id(gr[[gene_col]])
    # Merge selected columns (keep all other columns too; Cytoscape can ignore extras)
    nodes <- merge(nodes, gr, by.x = "id", by.y = "Gene", all.x = TRUE, sort = FALSE)
  }
  
  # 2) JointPCA -> Axis_Class
  jointpca_path <- NA_character_
  if (prompt_yn("Do you have a JointPCA loadings CSV to derive Axis_Class (PC1-stable/volatile/PC2-dominant)?", TRUE)) {
    cat("\nOPTIONAL FILE: JointPCA_Agroup.csv\n")
    cat("Expected: gene column + multiple *PC1 and *PC2 columns (e.g., 4_8PC1, 4_8PC2, ...)\n")
    jointpca_path <- pick_file("Select JointPCA loadings CSV:")
    ax <- derive_axis_class(jointpca_path)
    nodes <- merge(nodes, ax, by.x = "id", by.y = "Gene", all.x = TRUE, sort = FALSE)
  }
  
  # 3) DTW cluster table
  dtw_path <- NA_character_
  if (prompt_yn("Do you have a DTW cluster assignment CSV (Gene, DTW_cluster) to merge?", FALSE)) {
    cat("\nOPTIONAL FILE: DTW cluster table\n")
    cat("Expected columns: Gene, DTW_cluster\n")
    dtw_path <- pick_file("Select DTW cluster CSV:")
    dtw <- read_csv_checked(dtw_path)
    gene_col <- if ("Gene" %in% names(dtw)) "Gene" else names(dtw)[1]
    dtw$Gene <- clean_gene_id(dtw[[gene_col]])
    # Try to find a DTW cluster column
    dtw_col <- if ("DTW_cluster" %in% names(dtw)) "DTW_cluster" else choose_column(dtw, "DTW cluster column", hint_patterns = c("DTW", "cluster"))
    dtw2 <- data.frame(Gene = dtw$Gene, DTW_cluster = dtw[[dtw_col]], stringsAsFactors = FALSE)
    nodes <- merge(nodes, dtw2, by.x = "id", by.y = "Gene", all.x = TRUE, sort = FALSE)
  }
  
  # 4) Function mapping table (optional override/augment)
  funcmap_path <- NA_character_
  if (prompt_yn("Do you have a separate gene-to-function mapping CSV to merge/override Functional_Category?", FALSE)) {
    cat("\nOPTIONAL FILE: Gene-to-function mapping table\n")
    cat("Expected: Gene + Functional_Category (or similar)\n")
    funcmap_path <- pick_file("Select gene-to-function mapping CSV:")
    fm <- read_csv_checked(funcmap_path)
    gene_col <- if ("Gene" %in% names(fm)) "Gene" else names(fm)[1]
    fm$Gene <- clean_gene_id(fm[[gene_col]])
    func_col <- if ("Functional_Category" %in% names(fm)) "Functional_Category" else choose_column(fm, "Functional category column", hint_patterns = c("Function", "Category"))
    fm2 <- data.frame(Gene = fm$Gene, Functional_Category_map = fm[[func_col]], stringsAsFactors = FALSE)
    nodes <- merge(nodes, fm2, by.x = "id", by.y = "Gene", all.x = TRUE, sort = FALSE)
    
    # If nodes already has Functional_Category from GeneRole, offer override
    if ("Functional_Category" %in% names(nodes)) {
      if (prompt_yn("Override existing Functional_Category with mapping table where available?", FALSE)) {
        has_map <- !is.na(nodes$Functional_Category_map) & nzchar(as.character(nodes$Functional_Category_map))
        nodes$Functional_Category[has_map] <- nodes$Functional_Category_map[has_map]
      }
    } else {
      # Promote mapped column
      nodes$Functional_Category <- nodes$Functional_Category_map
    }
  }
  
  # ---- Diagnostics: unmatched annotation
  ann_present <- setdiff(names(nodes), "id")
  nodes_dir <- file.path(out_dir, "nodes")
  nodes_master_path <- file.path(nodes_dir, "nodes_master.csv")
  
  write.csv(nodes, nodes_master_path, row.names = FALSE)
  
  # Unmatched lists: genes in edges not present in annotations are already in nodes table (with NA attributes)
  # But we can report genes that had *any* annotation matched vs none.
  if (length(ann_present) > 0) {
    any_ann <- apply(nodes[, ann_present, drop = FALSE], 1, function(r) any(!is.na(r) & nzchar(as.character(r))))
    unmatched_edges <- nodes$id[!any_ann]
    write.csv(data.frame(Gene = unmatched_edges, stringsAsFactors = FALSE),
              file.path(nodes_dir, "unmatched_genes_in_annotations.csv"), row.names = FALSE)
  }
  
  # Also list genes in annotation tables not present in filtered edges (if GeneRole supplied)
  if (!is.na(gene_role_path)) {
    gr2 <- read_csv_checked(gene_role_path)
    gene_col <- if ("Gene" %in% names(gr2)) "Gene" else names(gr2)[1]
    gset <- sort(unique(clean_gene_id(gr2[[gene_col]])))
    missing_in_edges <- setdiff(gset, nodes$id)
    write.csv(data.frame(Gene = missing_in_edges, stringsAsFactors = FALSE),
              file.path(nodes_dir, "unmatched_genes_in_edges.csv"), row.names = FALSE)
  }
  
  cat("\nNodes written:\n")
  cat("  ", normalizePath(nodes_master_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  
  # ---- Optional: compute per-slice metrics (betweenness, degree) and bridge calls
  metrics_dir <- file.path(out_dir, "metrics")
  do_metrics <- prompt_yn("\nCompute per-slice network metrics (degree, betweenness) to depict transient bridging?", TRUE)
  
  node_metrics <- NULL
  bridge_calls <- NULL
  
  if (do_metrics) {
    cat("\nMetrics help depict bridging as a time-local role.\n")
    cat("We compute metrics separately for each slice network built from FILTERED edges.\n")
    
    # Estimate compute size per slice
    slices_f <- sort(unique(edges_f$slice))
    edges_per_slice <- vapply(slices_f, function(s) sum(edges_f$slice == s), numeric(1))
    
    # Core availability and recommendation
    avail <- if (ok_parallel) parallel::detectCores(logical = TRUE) else 1
    rec <- suggest_cores(length(slices_f), edges_per_slice)
    
    cat("\nCores detected: ", rec$avail, "\n", sep = "")
    cat("Filtered edges per slice: median ", format(rec$med_edges, big.mark = ","), ", max ", format(rec$max_edges, big.mark = ","), "\n", sep = "")
    cat("Recommended cores for metric computation (speed vs memory thrash): ", rec$recommended, "\n", sep = "")
    
    use_parallel <- prompt_yn("Use parallel computation for per-slice metrics?", default_yes = (rec$recommended > 1))
    n_cores <- 1L
    if (use_parallel && ok_parallel) {
      n_cores <- as.integer(prompt_number("How many cores to use?", default_val = rec$recommended, min_val = 1, max_val = rec$avail))
    } else {
      use_parallel <- FALSE
      n_cores <- 1L
    }
    
    # Bridge calling rule
    cat("\nBridge calling rule:\n")
    cat("This flags the most integrative nodes within each slice (does NOT imply permanent hubs).\n")
    bc_choices <- c(
      "Top K betweenness per slice (simple, interpretable)",
      "Top quantile of betweenness per slice (scale-aware)"
    )
    bc_idx <- prompt_choice("Choose a bridge-calling rule:", bc_choices, default_idx = 1)
    top_k <- 10L
    top_q <- 0.05
    if (bc_idx == 1) {
      top_k <- as.integer(prompt_number("Enter K (top K nodes per slice)", default_val = 10, min_val = 1, max_val = 100000))
    } else {
      top_q <- prompt_number("Enter quantile (e.g., 0.05 = top 5%)", default_val = 0.05, min_val = 0.001, max_val = 0.5)
    }
    
    # Worker function: compute metrics for one slice
    compute_one_slice <- function(s) {
      sub <- edges_f[edges_f$slice == s, c("source", "target", "weight"), drop = FALSE]
      # Build igraph graph (undirected correlation network)
      g <- igraph::graph_from_data_frame(sub[, c("source", "target")], directed = FALSE)
      
      # Degree
      deg <- igraph::degree(g)
      
      # Betweenness: unweighted by default for interpretability.
      # (Weighted betweenness exists but changes meaning; keep unweighted unless you explicitly need it.)
      btw <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)
      
      data.frame(
        Gene = names(deg),
        slice = s,
        degree = as.numeric(deg),
        betweenness = as.numeric(btw),
        stringsAsFactors = FALSE
      )
    }
    
    # Compute metrics across slices (parallel optional)
    if (use_parallel && n_cores > 1L && ok_parallel) {
      # mclapply works on macOS/Linux; if not supported, fallback to lapply
      metrics_list <- tryCatch(
        parallel::mclapply(slices_f, compute_one_slice, mc.cores = n_cores),
        error = function(e) {
          cat("Parallel execution failed; falling back to single-core.\n")
          lapply(slices_f, compute_one_slice)
        }
      )
    } else {
      metrics_list <- lapply(slices_f, compute_one_slice)
    }
    
    node_metrics <- do.call(rbind, metrics_list)
    
    # Bridge calls per slice
    node_metrics$is_bridge <- FALSE
    for (s in slices_f) {
      idx <- which(node_metrics$slice == s)
      b <- node_metrics$betweenness[idx]
      if (length(b) == 0) next
      if (bc_idx == 1) {
        # Top K
        ord <- order(b, decreasing = TRUE)
        k <- min(top_k, length(ord))
        node_metrics$is_bridge[idx[ord[seq_len(k)]]] <- TRUE
      } else {
        # Top quantile
        thr <- as.numeric(stats::quantile(b, probs = 1 - top_q, na.rm = TRUE))
        node_metrics$is_bridge[idx[b >= thr]] <- TRUE
      }
    }
    
    node_metrics_path <- file.path(metrics_dir, "node_metrics_by_slice.csv")
    bridge_calls_path <- file.path(metrics_dir, "bridge_calls_by_slice.csv")
    
    write.csv(node_metrics, node_metrics_path, row.names = FALSE)
    write.csv(node_metrics[node_metrics$is_bridge, c("Gene", "slice", "betweenness", "degree", "is_bridge")],
              bridge_calls_path, row.names = FALSE)
    
    cat("\nMetrics written:\n")
    cat("  ", normalizePath(node_metrics_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
    cat("  ", normalizePath(bridge_calls_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  }
  
  # ---- Provenance: manifest + inventory
  manifest <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script_name = ident$script_name,
    script_path = ident$script_path,
    script_full = ifelse(is.na(ident$script_full), "", ident$script_full),
    out_dir = out_dir,
    input_edges = edge_path,
    input_gene_role = ifelse(is.na(gene_role_path), "", gene_role_path),
    input_jointpca = ifelse(is.na(jointpca_path), "", jointpca_path),
    input_dtw = ifelse(is.na(dtw_path), "", dtw_path),
    input_funcmap = ifelse(is.na(funcmap_path), "", funcmap_path),
    edge_source_col = source_col,
    edge_target_col = target_col,
    edge_slice_col  = slice_col,
    edge_weight_col = weight_col,
    filter_mode = filter_desc,
    filter_params = filter_params,
    metrics_computed = do_metrics
  )
  manifest_path <- write_manifest_csv(out_dir, manifest)
  inventory_path <- write_inventory_csv(out_dir)
  
  cat("\nProvenance written:\n")
  cat("  ", manifest_path, "\n", sep = "")
  cat("  ", inventory_path, "\n", sep = "")
  
  # ---- Optional Quarto report (your contract; safe, minimal)
  do_report <- prompt_yn("\nGenerate Quarto QMD + HTML report summarizing inputs/parameters/outputs?", TRUE)
  qmd_path <- NA_character_
  html_path <- NA_character_
  
  if (do_report) {
    report_dir <- file.path(out_dir, "report")
    qmd_path <- file.path(report_dir, "cytoscape_export_report.qmd")
    
    # Capture script header block (first comment block)
    header_text <- ""
    if (!is.na(ident$script_full) && file.exists(ident$script_full)) {
      lines <- readLines(ident$script_full, warn = FALSE)
      take_n <- min(length(lines), 200)
      head_lines <- lines[seq_len(take_n)]
      keep <- character()
      for (ln in head_lines) {
        if (grepl("^\\s*#|^\\s*$", ln)) keep <- c(keep, ln) else break
      }
      header_text <- paste(keep, collapse = "\n")
    } else {
      header_text <- "# Script header unavailable (script path not detected in this session)."
    }
    
    # Build QMD as a single character vector
    qmd <- c(
      "---",
      "title: \"Cytoscape Export Report\"",
      "format: html",
      "execute:",
      "  echo: false",
      "  warning: false",
      "  message: false",
      "---",
      "",
      "## Summary",
      "This report documents the creation of Cytoscape-compatible network files from a long-format gene–gene correlation table.",
      "",
      "## Script header (verbatim)",
      "```",
      header_text,
      "```",
      "",
      "## Metadata",
      paste0("- Timestamp: ", manifest$timestamp),
      paste0("- Output directory: ", manifest$out_dir),
      paste0("- Edge file: ", manifest$input_edges),
      paste0("- GeneRole table: ", ifelse(nzchar(manifest$input_gene_role), manifest$input_gene_role, "(not provided)")),
      paste0("- JointPCA table: ", ifelse(nzchar(manifest$input_jointpca), manifest$input_jointpca, "(not provided)")),
      paste0("- DTW table: ", ifelse(nzchar(manifest$input_dtw), manifest$input_dtw, "(not provided)")),
      paste0("- Function map: ", ifelse(nzchar(manifest$input_funcmap), manifest$input_funcmap, "(not provided)")),
      "",
      "## Parameters",
      paste0("- Source column: ", manifest$edge_source_col),
      paste0("- Target column: ", manifest$edge_target_col),
      paste0("- Slice column: ", manifest$edge_slice_col),
      paste0("- Weight column: ", manifest$edge_weight_col),
      paste0("- Filter mode: ", manifest$filter_mode),
      paste0("- Filter params: ", ifelse(nzchar(manifest$filter_params), manifest$filter_params, "(none)")),
      paste0("- Metrics computed: ", manifest$metrics_computed),
      "",
      "## Generated outputs",
      "```{r}",
      "out_dir <- '.'",
      "files <- list.files(out_dir, recursive = TRUE)",
      "files <- files[order(files)]",
      "print(files)",
      "```",
      "",
      "## Session info",
      "```{r}",
      "sessionInfo()",
      "```"
    )
    
    writeLines(qmd, qmd_path)
    
    if (ok_quarto) {
      # Canonical contract: Quarto writes to current WD.
      old_wd <- getwd()
      setwd(report_dir)
      on.exit(setwd(old_wd), add = TRUE)
      
      # Render in report_dir; QMD uses relative paths to itself
      tryCatch(
        quarto::quarto_render(basename(qmd_path)),
        error = function(e) {
          cat("Quarto render failed: ", e$message, "\n", sep = "")
        }
      )
      html_path <- sub("\\.qmd$", ".html", qmd_path)
      if (!file.exists(html_path)) html_path <- NA_character_
    } else {
      cat("Package 'quarto' not available; wrote QMD only.\n")
    }
    
    cat("\nReport written:\n")
    cat("  ", normalizePath(qmd_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
    if (!is.na(html_path)) cat("  ", normalizePath(html_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  }
  
  # ---- Return (for interactive use)
  list(
    out_dir = out_dir,
    edges_all = normalizePath(all_combined_path, winslash = "/", mustWork = TRUE),
    edges_filtered = normalizePath(filt_combined_path, winslash = "/", mustWork = TRUE),
    nodes_master = normalizePath(nodes_master_path, winslash = "/", mustWork = TRUE),
    manifest = manifest_path,
    inventory = inventory_path,
    qmd = ifelse(is.na(qmd_path), NA_character_, normalizePath(qmd_path, winslash = "/", mustWork = TRUE)),
    html = ifelse(is.na(html_path), NA_character_, normalizePath(html_path, winslash = "/", mustWork = TRUE))
  )
}

# ============================================================
# AUTO-RUN WHEN SOURCED (RStudio / interactive use)
# ============================================================
if (interactive()) {
  message("Running cytoscape_export_builder() interactively...")
  res <- cytoscape_export_builder()
  message("Done. Results written to:")
  message(res$out_dir)
  
  # Helpful: open node table in RStudio if available
  if (!is.null(res$nodes_master) && file.exists(res$nodes_master)) {
    tab <- tryCatch(read.csv(res$nodes_master, check.names = FALSE), error = function(e) NULL)
    if (!is.null(tab)) try(View(tab), silent = TRUE)
  }
  
  message("Key files:")
  message(paste0("  edges (filtered): ", res$edges_filtered))
  message(paste0("  nodes:            ", res$nodes_master))
  message(paste0("  manifest:         ", res$manifest))
  message(paste0("  inventory:        ", res$inventory))
  if (!is.na(res$html)) message(paste0("  report (html):    ", res$html))
}