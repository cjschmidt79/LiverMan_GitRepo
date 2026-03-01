#!/usr/bin/env Rscript
# ======================================================================
# Function–Function Edge Turnover Enrichment (Option 1: Gene-label permutation null)
#
# PURPOSE
#   Given an edge-turnover table with:
#     - a transition identifier (e.g., "Days" column like "8-10")
#     - a signed delta correlation (delta_r)
#     - two endpoints that are functional classes (e.g., GeneA/GeneB now contain classes)
#     - OPTIONAL original gene columns saved earlier (e.g., GeneA_orig, GeneB_orig)
#
#   This script:
#     1) Detects/asks for required columns
#     2) For each transition:
#         - Builds observed function–function summaries:
#             * n_edges per function-pair
#             * sum_abs_delta = Σ|delta_r|
#             * sign balance: n_pos, n_neg, net_signed = Σdelta_r
#         - Tests enrichment vs a null where functional labels are permuted across genes
#           (Option 1; preserves gene pairs + delta_r values, tests functional specificity)
#         - Produces per-pair empirical p-values and z-scores for:
#             * n_edges
#             * sum_abs_delta
#         - Produces per-transition global metrics:
#             * concentration / entropy of function-pair distribution
#             * deviation-from-null (z) of entropy (optional)
#     3) Writes all outputs to a timestamped subdirectory of the working directory
#
# OUTPUTS (CSV)
#   - Transition_Summary.csv
#   - Pair_Enrichment_nEdges.csv
#   - Pair_Enrichment_sumAbsDelta.csv
#   - Pair_SignBalance.csv
#   - Diagnostics_RunInfo.txt
#
# NOTES
#   - Recommended for your "top-50 edges per transition" setup.
#   - Requires original gene columns (GeneA_orig/GeneB_orig) to define the permutation null properly.
#     If you do NOT have original gene columns, the script can still run a weaker null by shuffling
#     endpoint classes directly (not recommended). This script will warn you and offer the fallback.
# ======================================================================
source("R/project_utils.R")
THIS_SCRIPT <- "scripts/edge_turnover_function_enrichment_permnull.R"  # <-- change to your actual filename

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)

# ----------------------------- Manifest helpers -----------------------------

get_script_info <- function(explicit_path = getOption("this_script_path", NA_character_)) {
  
  # 0) Explicit path wins (robust for RStudio Source)
  script_path <- explicit_path
  if (!is.na(script_path) && nzchar(script_path) && file.exists(script_path)) {
    script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
    return(list(name = basename(script_path), path = script_path))
  }
  
  # 1) Rscript case (--file=)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_path <- tryCatch(normalizePath(script_path, winslash = "/", mustWork = TRUE),
                            error = function(e) script_path)
    return(list(name = basename(script_path), path = script_path))
  }
  
  # 2) RStudio active document fallback
  if (interactive() &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NA_character_)
    if (!is.na(p) && nzchar(p) && file.exists(p)) {
      script_path <- normalizePath(p, winslash = "/", mustWork = TRUE)
      return(list(name = basename(script_path), path = script_path))
    }
  }
  
  # 3) Give up
  list(name = "interactive_session", path = NA_character_)
}


build_in_meta <- function(path, df = NULL) {
  info <- file.info(path)
  list(
    path = path,
    name = basename(path),
    size_bytes = if (!is.na(info$size)) as.numeric(info$size) else NA_real_,
    modified_time = if (!is.na(info$mtime)) as.character(info$mtime) else NA_character_,
    n_rows = if (!is.null(df)) nrow(df) else NA_integer_,
    n_cols = if (!is.null(df)) ncol(df) else NA_integer_,
    colnames = if (!is.null(df)) names(df) else NULL
  )
}

write_manifest <- function(manifest, out_dir, filename = "Project_Manifest.json") {
  json_path <- file.path(out_dir, filename)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(manifest, json_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    # Fallback if jsonlite is not installed: write an RDS
    saveRDS(manifest, file.path(out_dir, sub("\\.json$", ".rds", filename)))
    json_path <- file.path(out_dir, sub("\\.json$", ".rds", filename))
  }
  json_path
}

read_csv_safely <- function(path) {
  if (!file.exists(path)) stop2("File not found: ", path)
  df <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop2("Failed to read CSV: ", path, "\n", e$message)
  )
  if (nrow(df) == 0) stop2("CSV has 0 rows: ", path)
  df
}

pick_file <- function(prompt_txt) {
  cat("\n", prompt_txt, "\n", sep = "")
  p <- tryCatch(file.choose(), error = function(e) NA_character_)
  if (is.na(p) || !nzchar(p)) stop2("No file selected.")
  normalizePath(p, winslash = "/", mustWork = TRUE)
}

choose_from_list <- function(options, prompt_txt) {
  if (length(options) == 0) stop2("No options available to choose from.")
  cat("\n", prompt_txt, "\n", sep = "")
  for (i in seq_along(options)) cat(sprintf("  [%d] %s\n", i, options[i]))
  ans <- readline("Enter selection number: ")
  idx <- suppressWarnings(as.integer(ans))
  if (is.na(idx) || idx < 1 || idx > length(options)) stop2("Invalid selection.")
  options[idx]
}

guess_days_col <- function(nms) {
  cand <- c("Days","days","Transition","transition","DayPair","daypair","A_B","A-B","AtoB","A_to_B")
  hit <- intersect(cand, nms)
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

guess_delta_col <- function(nms) {
  cand <- c("delta_r","Delta_r","deltaR","DeltaR","dR","dr","DeltaCorr","delta_corr")
  hit <- intersect(cand, nms)
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

guess_endpoints_cols <- function(nms) {
  pairs <- list(
    c("GeneA","GeneB"),
    c("geneA","geneB"),
    c("Gene1","Gene2"),
    c("gene1","gene2"),
    c("From","To"),
    c("from","to"),
    c("source","target"),
    c("Source","Target"),
    c("node1","node2"),
    c("Node1","Node2")
  )
  for (p in pairs) if (all(p %in% nms)) return(p)
  c(NA_character_, NA_character_)
}

guess_orig_cols <- function(nms, ep1, ep2) {
  c1 <- paste0(ep1, "_orig")
  c2 <- paste0(ep2, "_orig")
  if (c1 %in% nms && c2 %in% nms) return(c(c1, c2))
  
  # common fallbacks
  cand_pairs <- list(
    c("GeneA_orig","GeneB_orig"),
    c("geneA_orig","geneB_orig"),
    c("Gene1_orig","Gene2_orig"),
    c("gene1_orig","gene2_orig"),
    c("From_orig","To_orig"),
    c("from_orig","to_orig")
  )
  for (p in cand_pairs) if (all(p %in% nms)) return(p)
  
  c(NA_character_, NA_character_)
}

canon_pair <- function(a, b, sep="__") {
  # undirected canonical pair label
  ifelse(a <= b, paste0(a, sep, b), paste0(b, sep, a))
}

entropy_shannon <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

# Build per-transition observed summaries
summarize_transition <- function(df_t, ep1, ep2, delta_col) {
  a <- trimws(as.character(df_t[[ep1]]))
  b <- trimws(as.character(df_t[[ep2]]))
  d <- as.numeric(df_t[[delta_col]])
  if (anyNA(d)) stop2("delta_r contains NA within a transition; clean data first.")
  
  pair <- canon_pair(a, b)
  absd <- abs(d)
  pos <- d > 0
  neg <- d < 0
  
  tab_n <- tapply(rep(1, length(pair)), pair, sum)
  tab_abs <- tapply(absd, pair, sum)
  tab_signed <- tapply(d, pair, sum)
  tab_pos <- tapply(as.integer(pos), pair, sum)
  tab_neg <- tapply(as.integer(neg), pair, sum)
  
  # Align
  pairs <- sort(unique(pair))
  out <- data.frame(
    pair = pairs,
    n_edges = as.integer(tab_n[pairs]),
    sum_abs_delta = as.numeric(tab_abs[pairs]),
    sum_signed_delta = as.numeric(tab_signed[pairs]),
    n_pos = as.integer(tab_pos[pairs]),
    n_neg = as.integer(tab_neg[pairs]),
    stringsAsFactors = FALSE
  )
  out$frac_pos <- ifelse(out$n_edges > 0, out$n_pos / out$n_edges, NA_real_)
  out$frac_neg <- ifelse(out$n_edges > 0, out$n_neg / out$n_edges, NA_real_)
  out
}

# Permutation null: shuffle function labels across genes, preserve gene pairs and delta_r
permute_once <- function(genes1, genes2, gene_to_func_perm, delta_vec) {
  f1 <- unname(gene_to_func_perm[genes1])
  f2 <- unname(gene_to_func_perm[genes2])
  pair <- canon_pair(f1, f2)
  absd <- abs(delta_vec)
  
  # Return two named vectors:
  #  - counts per pair
  #  - sum_abs per pair
  list(
    n = tapply(rep(1, length(pair)), pair, sum),
    abs = tapply(absd, pair, sum)
  )
}

# Empirical p-value helper (two-sided by default)
emp_p_two_sided <- function(null_vec, obs_val) {
  # two-sided using absolute deviation from null mean
  mu <- mean(null_vec)
  dev_obs <- abs(obs_val - mu)
  dev_null <- abs(null_vec - mu)
  (sum(dev_null >= dev_obs) + 1) / (length(null_vec) + 1)
}

# ------------------------------ MAIN ------------------------------

cat("=== Function–Function Edge Turnover Enrichment (Option 1) ===\n")
cat("Working directory:\n  ", getwd(), "\n", sep = "")

in_path <- pick_file("Select the EDGE TURNOVER CSV to process (function-swapped table).")
df <- read_csv_safely(in_path)
nms <- names(df)

# Identify columns
days_col <- guess_days_col(nms)
delta_col <- guess_delta_col(nms)
eps <- guess_endpoints_cols(nms)

if (is.na(days_col))  days_col  <- choose_from_list(nms, "Which column indicates the transition (e.g., '8-10')?")
if (is.na(delta_col)) delta_col <- choose_from_list(nms, "Which column contains delta_r (signed change in correlation)?")
if (any(is.na(eps))) {
  ep1 <- choose_from_list(nms, "Select endpoint column 1 (currently functional class labels).")
  ep2 <- choose_from_list(nms, "Select endpoint column 2 (currently functional class labels).")
} else {
  ep1 <- eps[1]; ep2 <- eps[2]
  cat("\nDetected endpoint columns:\n  ", ep1, " and ", ep2, "\n", sep = "")
}

origs <- guess_orig_cols(nms, ep1, ep2)
orig1 <- origs[1]; orig2 <- origs[2]

use_fallback <- FALSE
if (any(is.na(c(orig1, orig2)))) {
  cat("\nWARNING: Original gene columns (e.g., ", ep1, "_orig / ", ep2, "_orig) were not found.\n", sep = "")
  cat("Option 1 (recommended) requires original genes to permute function labels across genes.\n")
  cat("Fallback option: shuffle endpoint functional labels directly within each transition (WEAKER null).\n")
  ans <- tolower(trimws(readline("Proceed with fallback? (yes/no): ")))
  if (ans %in% c("y","yes")) {
    use_fallback <- TRUE
  } else {
    stop2("Stopping because original gene columns are required for Option 1.")
  }
} else {
  cat("\nDetected original gene columns:\n  ", orig1, " and ", orig2, "\n", sep = "")
}

# Ask permutations
n_perm <- suppressWarnings(as.integer(trimws(readline("\nNumber of permutations per transition? [default 2000]: "))))
if (is.na(n_perm) || n_perm <= 0) n_perm <- 2000L
cat("Using n_perm = ", n_perm, "\n", sep = "")

set_seed_txt <- trimws(readline("Random seed? [default 1]: "))
seed <- suppressWarnings(as.integer(set_seed_txt))
if (is.na(seed)) seed <- 1L
set.seed(seed)

proj <- init_project_output(
  subdir = paste0("FunctionEnrichment_EdgeTurnover_", ts_stamp())
)

out_dir <- proj$output_root
# ----------------------------- Manifest assembly (MUST happen once) -----------------------------

out_meta <- list(
  path = out_dir,
  script_info = get_script_info(explicit_path = THIS_SCRIPT),
  timestamp = ts_stamp(),
  r_version = R.version.string
)

in_meta <- list(
  edge_turnover_csv = build_in_meta(in_path, df)
)

project_manifest <- list(
  script_name = out_meta$script_info$name,
  script_path = out_meta$script_info$path,
  output_dir  = out_meta$path,
  input_files = in_meta,
  run_parameters = list(
    n_perm = n_perm,
    seed = seed,
    days_col = days_col,
    delta_col = delta_col,
    endpoint_cols = c(ep1, ep2),
    orig_gene_cols = if (use_fallback) NA_character_ else c(orig1, orig2),
    null_model = if (use_fallback) "Fallback (shuffle endpoint classes)" else "Option 1 (permute function labels across genes)"
  )
)


cat("\nOutput directory:\n  ", out_dir, "\n", sep = "")


# Clean required columns
df[[days_col]] <- trimws(as.character(df[[days_col]]))
df[[ep1]] <- trimws(as.character(df[[ep1]]))
df[[ep2]] <- trimws(as.character(df[[ep2]]))
df[[delta_col]] <- suppressWarnings(as.numeric(df[[delta_col]]))
if (anyNA(df[[delta_col]])) stop2("delta_r column contains non-numeric values or NA after coercion.")

# Per-transition loop
transitions <- sort(unique(df[[days_col]]))
cat("\nFound ", length(transitions), " transitions.\n", sep = "")

pair_enrich_n_all <- list()
pair_enrich_abs_all <- list()
pair_sign_all <- list()
transition_summary <- list()

for (tr in transitions) {
  cat("\n--- Transition: ", tr, " ---\n", sep = "")
  df_t <- df[df[[days_col]] == tr, , drop = FALSE]
  if (nrow(df_t) == 0) next
  
  # Observed per-pair summary
  obs_pairs <- summarize_transition(df_t, ep1, ep2, delta_col)
  obs_pairs$transition <- tr
  
  # Global concentration metrics (based on n_edges and sum_abs)
  p_n <- obs_pairs$n_edges / sum(obs_pairs$n_edges)
  H_n <- entropy_shannon(p_n)
  p_abs <- obs_pairs$sum_abs_delta / sum(obs_pairs$sum_abs_delta)
  H_abs <- entropy_shannon(p_abs)
  
  # Build null distributions for each observed pair
  # We'll store null values per pair for:
  #   - n_edges
  #   - sum_abs_delta
  obs_pair_names <- obs_pairs$pair
  null_n <- matrix(0, nrow = length(obs_pair_names), ncol = n_perm,
                   dimnames = list(obs_pair_names, NULL))
  null_abs <- matrix(0, nrow = length(obs_pair_names), ncol = n_perm,
                     dimnames = list(obs_pair_names, NULL))
  
  if (!use_fallback) {
    # Option 1: permute function labels across genes (within this transition's involved genes)
    g1 <- trimws(as.character(df_t[[orig1]]))
    g2 <- trimws(as.character(df_t[[orig2]]))
    dvec <- as.numeric(df_t[[delta_col]])
    
    # Map observed gene->function using current swapped endpoints
    # (This assumes swapped endpoints are consistent with the mapping used earlier)
    # Create gene->function from both sides
    f1_obs <- trimws(as.character(df_t[[ep1]]))
    f2_obs <- trimws(as.character(df_t[[ep2]]))
    
    gene_to_func <- c(setNames(f1_obs, g1), setNames(f2_obs, g2))
    # Resolve conflicts (same gene mapped to multiple classes) by first occurrence; warn
    dup_genes <- names(gene_to_func)[duplicated(names(gene_to_func))]
    if (length(dup_genes) > 0) {
      cat("WARNING: Some genes map to multiple functional classes within transition (keeping first). Example(s): ",
          paste(unique(head(dup_genes, 8)), collapse = ", "), "\n", sep = "")
      gene_to_func <- gene_to_func[!duplicated(names(gene_to_func))]
    }
    
    genes_all <- names(gene_to_func)
    funcs_all <- unname(gene_to_func)
    
    for (b in seq_len(n_perm)) {
      perm_funcs <- sample(funcs_all, size = length(funcs_all), replace = FALSE)
      gene_to_func_perm <- setNames(perm_funcs, genes_all)
      
      res <- permute_once(g1, g2, gene_to_func_perm, dvec)
      
      # Fill only for observed pairs (others treated as 0)
      n_vec <- rep(0, length(obs_pair_names)); names(n_vec) <- obs_pair_names
      abs_vec <- rep(0, length(obs_pair_names)); names(abs_vec) <- obs_pair_names
      
      if (!is.null(res$n)) {
        common <- intersect(names(res$n), obs_pair_names)
        n_vec[common] <- as.numeric(res$n[common])
      }
      if (!is.null(res$abs)) {
        common <- intersect(names(res$abs), obs_pair_names)
        abs_vec[common] <- as.numeric(res$abs[common])
      }
      
      null_n[, b] <- n_vec
      null_abs[, b] <- abs_vec
    }
    
  } else {
    # Fallback: shuffle endpoint classes directly (weaker null; does not preserve gene->function structure)
    f1 <- trimws(as.character(df_t[[ep1]]))
    f2 <- trimws(as.character(df_t[[ep2]]))
    dvec <- as.numeric(df_t[[delta_col]])
    
    for (b in seq_len(n_perm)) {
      f1p <- sample(f1, length(f1), replace = FALSE)
      f2p <- sample(f2, length(f2), replace = FALSE)
      pairp <- canon_pair(f1p, f2p)
      absd <- abs(dvec)
      
      res_n <- tapply(rep(1, length(pairp)), pairp, sum)
      res_abs <- tapply(absd, pairp, sum)
      
      n_vec <- rep(0, length(obs_pair_names)); names(n_vec) <- obs_pair_names
      abs_vec <- rep(0, length(obs_pair_names)); names(abs_vec) <- obs_pair_names
      
      if (!is.null(res_n)) {
        common <- intersect(names(res_n), obs_pair_names)
        n_vec[common] <- as.numeric(res_n[common])
      }
      if (!is.null(res_abs)) {
        common <- intersect(names(res_abs), obs_pair_names)
        abs_vec[common] <- as.numeric(res_abs[common])
      }
      
      null_n[, b] <- n_vec
      null_abs[, b] <- abs_vec
    }
  }
  
  # Per-pair enrichment stats
  obs_n <- obs_pairs$n_edges; names(obs_n) <- obs_pairs$pair
  obs_abs <- obs_pairs$sum_abs_delta; names(obs_abs) <- obs_pairs$pair
  
  mu_n <- rowMeans(null_n); sd_n <- apply(null_n, 1, sd)
  z_n <- (obs_n - mu_n) / ifelse(sd_n == 0, NA_real_, sd_n)
  p_n_emp <- mapply(function(pn, ob) emp_p_two_sided(pn, ob),
                    split(null_n, row(null_n)), obs_n, SIMPLIFY = TRUE)
  
  mu_abs <- rowMeans(null_abs); sd_abs <- apply(null_abs, 1, sd)
  z_abs <- (obs_abs - mu_abs) / ifelse(sd_abs == 0, NA_real_, sd_abs)
  p_abs_emp <- mapply(function(pa, ob) emp_p_two_sided(pa, ob),
                      split(null_abs, row(null_abs)), obs_abs, SIMPLIFY = TRUE)
  
  enrich_n <- data.frame(
    transition = tr,
    pair = obs_pairs$pair,
    obs_n_edges = obs_pairs$n_edges,
    null_mean = as.numeric(mu_n[obs_pairs$pair]),
    null_sd = as.numeric(sd_n[obs_pairs$pair]),
    z = as.numeric(z_n[obs_pairs$pair]),
    p_empirical_2sided = as.numeric(p_n_emp),
    stringsAsFactors = FALSE
  )
  enrich_abs <- data.frame(
    transition = tr,
    pair = obs_pairs$pair,
    obs_sum_abs_delta = obs_pairs$sum_abs_delta,
    null_mean = as.numeric(mu_abs[obs_pairs$pair]),
    null_sd = as.numeric(sd_abs[obs_pairs$pair]),
    z = as.numeric(z_abs[obs_pairs$pair]),
    p_empirical_2sided = as.numeric(p_abs_emp),
    stringsAsFactors = FALSE
  )
  
  # Sign-balance (observed only)
  sign_tab <- obs_pairs[, c("transition","pair","n_edges","n_pos","n_neg","frac_pos","frac_neg","sum_signed_delta")]
  sign_tab$direction_bias <- ifelse(sign_tab$n_edges > 0, (sign_tab$n_pos - sign_tab$n_neg) / sign_tab$n_edges, NA_real_)
  
  # Transition summary
  trans_sum <- data.frame(
    transition = tr,
    n_edges = nrow(df_t),
    n_unique_pairs = nrow(obs_pairs),
    entropy_nEdges = H_n,
    entropy_sumAbsDelta = H_abs,
    top_pair_by_n = obs_pairs$pair[which.max(obs_pairs$n_edges)],
    top_pair_by_abs = obs_pairs$pair[which.max(obs_pairs$sum_abs_delta)],
    stringsAsFactors = FALSE
  )
  
  pair_enrich_n_all[[tr]] <- enrich_n
  pair_enrich_abs_all[[tr]] <- enrich_abs
  pair_sign_all[[tr]] <- sign_tab
  transition_summary[[tr]] <- trans_sum
}

# Bind outputs
Pair_Enrichment_nEdges <- do.call(rbind, pair_enrich_n_all)
Pair_Enrichment_sumAbsDelta <- do.call(rbind, pair_enrich_abs_all)
Pair_SignBalance <- do.call(rbind, pair_sign_all)
Transition_Summary <- do.call(rbind, transition_summary)

# Write outputs
write.csv(Transition_Summary, file.path(out_dir, "Transition_Summary.csv"), row.names = FALSE)
write.csv(Pair_Enrichment_nEdges, file.path(out_dir, "Pair_Enrichment_nEdges.csv"), row.names = FALSE)
write.csv(Pair_Enrichment_sumAbsDelta, file.path(out_dir, "Pair_Enrichment_sumAbsDelta.csv"), row.names = FALSE)
write.csv(Pair_SignBalance, file.path(out_dir, "Pair_SignBalance.csv"), row.names = FALSE)

# Diagnostics
diag_path <- file.path(out_dir, "Diagnostics_RunInfo.txt")
con <- file(diag_path, open = "wt")
writeLines(c(
  "Function–Function Edge Turnover Enrichment (Option 1) - Run Info",
  paste0("Timestamp: ", ts_stamp()),
  paste0("Input file: ", in_path),
  paste0("Working dir: ", getwd()),
  paste0("Output dir: ", out_dir),
  paste0("Transitions: ", length(transitions)),
  paste0("Days column: ", days_col),
  paste0("Delta column: ", delta_col),
  paste0("Endpoint columns: ", ep1, ", ", ep2),
  paste0("Original gene columns: ", ifelse(use_fallback, "NOT USED (fallback)", paste0(orig1, ", ", orig2))),
  paste0("Permutation count: ", n_perm),
  paste0("Seed: ", seed),
  paste0("Null model: ", ifelse(use_fallback, "Fallback (shuffle endpoint classes)", "Option 1 (permute function labels across genes)"))
), con)
close(con)


cat("\nWrote outputs:\n")
cat("  - ", file.path(out_dir, "Transition_Summary.csv"), "\n", sep = "")
cat("  - ", file.path(out_dir, "Pair_Enrichment_nEdges.csv"), "\n", sep = "")
cat("  - ", file.path(out_dir, "Pair_Enrichment_sumAbsDelta.csv"), "\n", sep = "")
cat("  - ", file.path(out_dir, "Pair_SignBalance.csv"), "\n", sep = "")
cat("  - ", diag_path, "\n", sep = "")

project_manifest$output_files <- list(
  Transition_Summary = file.path(out_dir, "Transition_Summary.csv"),
  Pair_Enrichment_nEdges = file.path(out_dir, "Pair_Enrichment_nEdges.csv"),
  Pair_Enrichment_sumAbsDelta = file.path(out_dir, "Pair_Enrichment_sumAbsDelta.csv"),
  Pair_SignBalance = file.path(out_dir, "Pair_SignBalance.csv"),
  Diagnostics = diag_path
)
manifest_path <- write_manifest(project_manifest, out_dir)

# ======================================================================
# Quarto report generation (self-contained; output_dir-safe)
# ======================================================================

if (requireNamespace("quarto", quietly = TRUE)) {
  
  # -----------------------------
  # 1. Capture absolute paths FIRST
  # -----------------------------
  out_dir_abs <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  
  paths <- list(
    Transition_Summary = file.path(out_dir_abs, "Transition_Summary.csv"),
    Pair_nEdges = file.path(out_dir_abs, "Pair_Enrichment_nEdges.csv"),
    Pair_sumAbs = file.path(out_dir_abs, "Pair_Enrichment_sumAbsDelta.csv"),
    Pair_Sign = file.path(out_dir_abs, "Pair_SignBalance.csv"),
    Diagnostics = file.path(out_dir_abs, "Diagnostics_RunInfo.txt"),
    Manifest = normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  )
  
  # -----------------------------
  # 2. Construct QMD content
  # -----------------------------
  qmd_lines <- c(
    "---",
    "title: \"Function–Function Edge Turnover Enrichment Report\"",
    "format:",
    "  html:",
    "    toc: true",
    "execute:",
    "  echo: false",
    "  message: false",
    "  warning: false",
    "---",
    "",
    "## Overview",
    "",
    "```{r}",
    sprintf("cat(readLines('%s'), sep = '\\n')", paths$Diagnostics),
    "```",
    "",
    "## Executed Script",
    "",
    "```{r}",
    sprintf("mp <- '%s'", paths$Manifest),
    "if (grepl('\\\\.rds$', mp, ignore.case = TRUE)) {",
    "  manifest <- readRDS(mp)",
    "} else {",
    "  manifest <- jsonlite::read_json(mp, simplifyVector = TRUE)",
    "}",
    "cat('Script name:', manifest$script_name, '\\n')",
    "cat('Script path:', manifest$script_path, '\\n')",
    "```",
    
    "",
    
    "## Generated Files",
    "",
    "```{r}",
    sprintf("print(readLines('%s'))", paths$Manifest),
    "```",
    "",
    "## Transition Summary",
    "",
    "```{r}",
    sprintf("knitr::kable(read.csv('%s'))", paths$Transition_Summary),
    "```",
    "",
    "## Function–Function Enrichment (Edge Counts)",
    "",
    "```{r}",
    sprintf("knitr::kable(head(read.csv('%s'), 20))", paths$Pair_nEdges),
    "```",
    "",
    "## Function–Function Enrichment (Σ|Δr|)",
    "",
    "```{r}",
    sprintf("knitr::kable(head(read.csv('%s'), 20))", paths$Pair_sumAbs),
    "```",
    "",
    "## Sign Balance",
    "",
    "```{r}",
    sprintf("knitr::kable(head(read.csv('%s'), 20))", paths$Pair_Sign),
    "```",
    "",
    "## Reproducibility",
    "",
    "```{r}",
    "sessionInfo()",
    "```"
  )
  
  # -----------------------------
  # 3. Write QMD into out_dir
  # -----------------------------
  qmd_path <- file.path(out_dir_abs, "report.qmd")
  writeLines(qmd_lines, qmd_path)
  
  # -----------------------------
  # 4. Render by temporarily switching wd
  # -----------------------------
  old_wd <- getwd()
  setwd(out_dir_abs)
  
  quarto::quarto_render(
    input = "report.qmd",
    output_file = "report.html"
  )
  
  setwd(old_wd)
  
  # -----------------------------
  # 5. Report final path
  # -----------------------------
  cat(
    "\nQuarto report written to:\n  ",
    normalizePath(file.path(out_dir_abs, "report.html"), winslash = "/"),
    "\n",
    sep = ""
  )
  
} else {
  cat("\nQuarto not installed; skipping report generation.\n")
}

cat("\nDone.\n")
