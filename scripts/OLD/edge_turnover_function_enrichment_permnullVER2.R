#!/usr/bin/env Rscript
# ======================================================================
# Function–Function Edge Turnover Enrichment (Option 1: Gene-label permutation null)
#Are the metabolome and transcriptome organized similarly, and at what lag (0, 1 step, 2 steps)?”
#Any manuscript claim about cross-layer coordination or timing offset.
#example input:
#EdgeTable_FunctionSwapped_Edgetunover_concatenatted_20260103_003124.csv
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
#     3) Writes all outputs under outputs/<run_id>/ and produces a Quarto-style HTML report
#
# OUTPUTS (CSV)
#   - Transition_Summary.csv
#   - Pair_Enrichment_nEdges.csv
#   - Pair_Enrichment_sumAbsDelta.csv
#   - Pair_SignBalance.csv
#   - Diagnostics_RunInfo.txt
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - <script_name>_Report.qmd
#   - <script_name>_Report.html
#
# NOTES
#   - Recommended for your "top-50 edges per transition" setup.
#   - Requires original gene columns (GeneA_orig/GeneB_orig) to define the permutation null properly.
#     If you do NOT have original gene columns, the script can still run a weaker null by shuffling
#     endpoint classes directly (not recommended). This script will warn you and offer the fallback.
# ======================================================================

# ----------------------------- Flags (MANDATORY) -----------------------------
enable_plotly <- TRUE
enable_runtime_tracking <- TRUE

start_time <- if (enable_runtime_tracking) Sys.time() else NULL

# ----------------------------- Utilities -----------------------------
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)

# ----------------------------- Script identity (MANDATORY block) -----------------------------
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
known_script_filename <- "edge_turnover_function_enrichment_permnull.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  .script_path_note <- "NOTE: Script full path could not be detected in this execution mode; using fallback filename + getwd()."
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
  .script_path_note <- NULL
}

# ----------------------------- Package management (MANDATORY) -----------------------------
ensure_pkg <- function(pkgs) {
  pkgs <- unique(pkgs)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  invisible(TRUE)
}

required_pkgs <- c("jsonlite", "rmarkdown", "knitr", "ggplot2")
if (enable_plotly) required_pkgs <- c(required_pkgs, "plotly")
ensure_pkg(required_pkgs)

pkg_versions <- function(pkgs) {
  data.frame(
    package = pkgs,
    version = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)),
    stringsAsFactors = FALSE
  )
}

# ----------------------------- Output directory policy (MANDATORY) -----------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

analysis_name <- "FunctionEnrichment_EdgeTurnover"
run_id <- paste0(analysis_name, "_", ts_stamp())

outputs_root <- ensure_dir(file.path(getwd(), "outputs"))
output_dir <- ensure_dir(file.path(outputs_root, run_id))

# ----------------------------- Manifest + inventory (MANDATORY) -----------------------------
write_manifest_json <- function(manifest, out_dir, filename = "Project_Manifest.json") {
  path <- file.path(out_dir, filename)
  jsonlite::write_json(manifest, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  path
}

write_file_inventory <- function(out_dir, filename = "Project_Manifest_Files.csv") {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, include.dirs = FALSE)
  if (length(files) == 0) {
    inv <- data.frame(
      rel_path = character(0),
      abs_path = character(0),
      size_bytes = numeric(0),
      modified_time = character(0),
      stringsAsFactors = FALSE
    )
    utils::write.csv(inv, file.path(out_dir, filename), row.names = FALSE)
    return(file.path(out_dir, filename))
  }
  info <- file.info(files)
  inv <- data.frame(
    rel_path = sub(paste0("^", gsub("\\\\", "/", normalizePath(out_dir, winslash = "/")), "/?"), "", gsub("\\\\", "/", files)),
    abs_path = gsub("\\\\", "/", files),
    size_bytes = as.numeric(info$size),
    modified_time = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
  inv <- inv[order(inv$rel_path), , drop = FALSE]
  utils::write.csv(inv, file.path(out_dir, filename), row.names = FALSE)
  file.path(out_dir, filename)
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
  ifelse(a <= b, paste0(a, sep, b), paste0(b, sep, a))
}

entropy_shannon <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

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

permute_once <- function(genes1, genes2, gene_to_func_perm, delta_vec) {
  f1 <- unname(gene_to_func_perm[genes1])
  f2 <- unname(gene_to_func_perm[genes2])
  pair <- canon_pair(f1, f2)
  absd <- abs(delta_vec)
  
  list(
    n = tapply(rep(1, length(pair)), pair, sum),
    abs = tapply(absd, pair, sum)
  )
}

emp_p_two_sided <- function(null_vec, obs_val) {
  mu <- mean(null_vec)
  dev_obs <- abs(obs_val - mu)
  dev_null <- abs(null_vec - mu)
  (sum(dev_null >= dev_obs) + 1) / (length(null_vec) + 1)
}

# ----------------------------- Header block capture for report (MANDATORY) -----------------------------
get_header_block <- function() {
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
    take_n <- min(length(lines), 220)
    head_lines <- lines[seq_len(take_n)]
    keep <- character(0)
    started <- FALSE
    for (ln in head_lines) {
      if (!started) {
        # start at first comment block line OR blank after shebang
        if (grepl("^\\s*#|^\\s*$", ln) || grepl("^#!/", ln)) {
          keep <- c(keep, ln)
          started <- TRUE
        } else {
          break
        }
      } else {
        if (grepl("^\\s*#|^\\s*$", ln) || grepl("^#!/", ln)) {
          keep <- c(keep, ln)
        } else {
          break
        }
      }
    }
    paste(keep, collapse = "\n")
  } else {
    paste0("# Header block unavailable (script_full not detected). Fallback: ", known_script_filename)
  }
}

# ------------------------------ MAIN ------------------------------
cat("=== Function–Function Edge Turnover Enrichment (Option 1) ===\n")
cat("Working directory:\n  ", normalizePath(getwd(), winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Outputs root:\n  ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Run output_dir:\n  ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Script identity:\n")
cat("  script_name: ", script_name, "\n", sep = "")
cat("  script_path: ", script_path, "\n", sep = "")
cat("  script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
if (!is.null(.script_path_note)) cat(.script_path_note, "\n")

# Input
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
  
  obs_pairs <- summarize_transition(df_t, ep1, ep2, delta_col)
  obs_pairs$transition <- tr
  
  p_n <- obs_pairs$n_edges / sum(obs_pairs$n_edges)
  H_n <- entropy_shannon(p_n)
  p_abs <- obs_pairs$sum_abs_delta / sum(obs_pairs$sum_abs_delta)
  H_abs <- entropy_shannon(p_abs)
  
  obs_pair_names <- obs_pairs$pair
  null_n <- matrix(0, nrow = length(obs_pair_names), ncol = n_perm,
                   dimnames = list(obs_pair_names, NULL))
  null_abs <- matrix(0, nrow = length(obs_pair_names), ncol = n_perm,
                     dimnames = list(obs_pair_names, NULL))
  
  if (!use_fallback) {
    g1 <- trimws(as.character(df_t[[orig1]]))
    g2 <- trimws(as.character(df_t[[orig2]]))
    dvec <- as.numeric(df_t[[delta_col]])
    
    f1_obs <- trimws(as.character(df_t[[ep1]]))
    f2_obs <- trimws(as.character(df_t[[ep2]]))
    
    gene_to_func <- c(setNames(f1_obs, g1), setNames(f2_obs, g2))
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
  
  sign_tab <- obs_pairs[, c("transition","pair","n_edges","n_pos","n_neg","frac_pos","frac_neg","sum_signed_delta")]
  sign_tab$direction_bias <- ifelse(sign_tab$n_edges > 0, (sign_tab$n_pos - sign_tab$n_neg) / sign_tab$n_edges, NA_real_)
  
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

Pair_Enrichment_nEdges <- do.call(rbind, pair_enrich_n_all)
Pair_Enrichment_sumAbsDelta <- do.call(rbind, pair_enrich_abs_all)
Pair_SignBalance <- do.call(rbind, pair_sign_all)
Transition_Summary <- do.call(rbind, transition_summary)

# Write outputs (ONLY within output_dir)
paths <- list(
  Transition_Summary = file.path(output_dir, "Transition_Summary.csv"),
  Pair_nEdges = file.path(output_dir, "Pair_Enrichment_nEdges.csv"),
  Pair_sumAbs = file.path(output_dir, "Pair_Enrichment_sumAbsDelta.csv"),
  Pair_Sign = file.path(output_dir, "Pair_SignBalance.csv"),
  Diagnostics = file.path(output_dir, "Diagnostics_RunInfo.txt")
)

utils::write.csv(Transition_Summary, paths$Transition_Summary, row.names = FALSE)
utils::write.csv(Pair_Enrichment_nEdges, paths$Pair_nEdges, row.names = FALSE)
utils::write.csv(Pair_Enrichment_sumAbsDelta, paths$Pair_sumAbs, row.names = FALSE)
utils::write.csv(Pair_SignBalance, paths$Pair_Sign, row.names = FALSE)

# Diagnostics
con <- file(paths$Diagnostics, open = "wt")
writeLines(c(
  "Function–Function Edge Turnover Enrichment (Option 1) - Run Info",
  paste0("Run ID: ", run_id),
  paste0("Timestamp: ", ts_stamp()),
  paste0("Input file: ", in_path),
  paste0("Working dir: ", normalizePath(getwd(), winslash = "/", mustWork = FALSE)),
  paste0("Output dir: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
  paste0("Transitions: ", length(transitions)),
  paste0("Days column: ", days_col),
  paste0("Delta column: ", delta_col),
  paste0("Endpoint columns: ", ep1, ", ", ep2),
  paste0("Original gene columns: ", ifelse(use_fallback, "NOT USED (fallback)", paste0(orig1, ", ", orig2))),
  paste0("Permutation count: ", n_perm),
  paste0("Seed: ", seed),
  paste0("Null model: ", ifelse(use_fallback, "Fallback (shuffle endpoint classes)", "Option 1 (permute function labels across genes)")),
  paste0("enable_plotly: ", enable_plotly),
  paste0("enable_runtime_tracking: ", enable_runtime_tracking)
), con)
close(con)

# ----------------------------- Always-on diagnostic plot (MANDATORY) -----------------------------
plot_counts_path <- file.path(output_dir, "Plot_EdgesByTransition.png")
try({
  dplot <- Transition_Summary
  if (!is.null(dplot) && nrow(dplot) > 0) {
    p <- ggplot2::ggplot(dplot, ggplot2::aes(x = transition, y = n_edges)) +
      ggplot2::geom_col() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Edge counts by transition", x = "Transition", y = "Number of edges") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(plot_counts_path, plot = p, width = 8.5, height = 4.8, dpi = 150)
  }
}, silent = TRUE)

# ----------------------------- Manifest (MANDATORY sections) -----------------------------
deps_df <- pkg_versions(required_pkgs)

manifest <- list(
  run_id = run_id,
  run_timestamp = as.character(Sys.time()),
  script = list(
    name = script_name,
    path = script_path,
    full_path = if (is.na(script_full)) NA_character_ else script_full
  ),
  input = list(
    edge_turnover_csv = in_path,
    columns = list(
      days_col = days_col,
      delta_col = delta_col,
      endpoint_cols = c(ep1, ep2),
      orig_gene_cols = if (use_fallback) NA_character_ else c(orig1, orig2)
    )
  ),
  parameters = list(
    n_perm = n_perm,
    seed = seed,
    null_model = if (use_fallback) "Fallback (shuffle endpoint classes)" else "Option 1 (permute function labels across genes)",
    enable_plotly = enable_plotly,
    enable_runtime_tracking = enable_runtime_tracking,
    runtime_seconds = NA_real_
  ),
  dependencies = deps_df,
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  ),
  output_files = list(
    Transition_Summary = paths$Transition_Summary,
    Pair_Enrichment_nEdges = paths$Pair_nEdges,
    Pair_Enrichment_sumAbsDelta = paths$Pair_sumAbs,
    Pair_SignBalance = paths$Pair_Sign,
    Diagnostics = paths$Diagnostics,
    Plot_EdgesByTransition = if (file.exists(plot_counts_path)) plot_counts_path else NA_character_
  )
)

manifest_path <- write_manifest_json(manifest, output_dir, filename = "Project_Manifest.json")

# Initial inventory (pre-report)
inventory_path <- write_file_inventory(output_dir, filename = "Project_Manifest_Files.csv")

# ----------------------------- Runtime tracking (MANDATORY behavior when enabled) -----------------------------
end_time <- if (enable_runtime_tracking) Sys.time() else NULL
runtime_seconds <- if (enable_runtime_tracking) round(as.numeric(difftime(end_time, start_time, units = "secs")), 2) else NA_real_

manifest$parameters$runtime_seconds <- runtime_seconds
manifest_path <- write_manifest_json(manifest, output_dir, filename = "Project_Manifest.json")

# ----------------------------- QMD + HTML report (MANDATORY) -----------------------------
# Report filenames
qmd_name <- paste0(script_name, "_Report.qmd")
html_name <- paste0(script_name, "_Report.html")
qmd_path <- file.path(output_dir, qmd_name)
html_path <- file.path(output_dir, html_name)

header_text <- get_header_block()

# QMD content with params (runtime_seconds included)
qmd_lines <- c(
  "---",
  paste0("title: \"Function–Function Edge Turnover Enrichment Report\""),
  "format:",
  "  html:",
  "    toc: true",
  "execute:",
  "  echo: false",
  "  message: false",
  "  warning: false",
  "params:",
  paste0("  run_id: \"", run_id, "\""),
  paste0("  run_timestamp: \"", as.character(Sys.time()), "\""),
  paste0("  script_name: \"", script_name, "\""),
  paste0("  script_path: \"", gsub("\"", "\\\"", script_path), "\""),
  paste0("  script_full: \"", gsub("\"", "\\\"", ifelse(is.na(script_full), "NA", script_full)), "\""),
  paste0("  input_file: \"", gsub("\"", "\\\"", normalizePath(in_path, winslash = "/", mustWork = TRUE)), "\""),
  paste0("  runtime_seconds: ", ifelse(is.na(runtime_seconds), "null", runtime_seconds)),
  paste0("  enable_plotly: ", ifelse(enable_plotly, "true", "false")),
  "---",
  "",
  "## Metadata",
  "```{r}",
  "cat('Run ID: ', params$run_id, '\\n', sep='')",
  "cat('Run timestamp: ', params$run_timestamp, '\\n', sep='')",
  "cat('Script name: ', params$script_name, '\\n', sep='')",
  "cat('Script path: ', params$script_path, '\\n', sep='')",
  "cat('Script full: ', params$script_full, '\\n', sep='')",
  "cat('Input file: ', params$input_file, '\\n', sep='')",
  "```",
  "",
  "## Script header",
  "```",
  header_text,
  "```",
  "",
  "## Analytical logic and formulas",
  "",
  "- **Observed summaries per function-pair** (within each transition):",
  "  - Edge count:  \\( n_{pair} = \\sum 1 \\)",
  "  - Magnitude:   \\( S_{pair} = \\sum |\\Delta r| \\)",
  "  - Sign balance: \\( n_{+}, n_{-}, \\sum \\Delta r \\)",
  "- **Null model (preferred)**: permute function labels across genes while preserving the original gene–gene pairs and \\(\\Delta r\\) values.",
  "- **Empirical two-sided p-value**: compares the observed deviation from the null mean to the null distribution of deviations.",
  "",
  "## Dependencies (packages + versions)",
  "```{r}",
  paste0("deps <- read.csv('", gsub("\\\\", "/", file.path(output_dir, "Project_Manifest.json")), "', stringsAsFactors = FALSE)"),
  "```",
  "```{r}",
  paste0("mp <- '", gsub("\\\\", "/", manifest_path), "'"),
  "manifest <- jsonlite::read_json(mp, simplifyVector = TRUE)",
  "knitr::kable(manifest$dependencies)",
  "```",
  "",
  "## Generated files (inventory)",
  "```{r}",
  paste0("inv <- read.csv('", gsub("\\\\", "/", inventory_path), "', stringsAsFactors = FALSE)"),
  "knitr::kable(inv)",
  "```",
  "",
  "## Runtime",
  "```{r}",
  "if (!is.null(params$runtime_seconds) && !is.na(params$runtime_seconds)) {",
  "  cat('Total session duration (seconds): ', params$runtime_seconds, '\\n', sep='')",
  "} else {",
  "  cat('Runtime tracking disabled or unavailable.\\n')",
  "}",
  "```",
  "",
  "## Diagnostics",
  "```{r}",
  paste0("cat(readLines('", gsub("\\\\", "/", paths$Diagnostics), "'), sep='\\n')"),
  "```",
  "",
  "## Required static plot (always-on)",
  "```{r}",
  paste0("if (file.exists('", gsub("\\\\", "/", plot_counts_path), "')) knitr::include_graphics('", gsub("\\\\", "/", plot_counts_path), "')"),
  "```",
  "",
  "## Optional interactive plot (plotly)",
  "```{r}",
  "ok <- FALSE",
  "if (isTRUE(params$enable_plotly)) {",
  "  ok <- tryCatch({",
  "    dfp <- read.csv(manifest$output_files$Transition_Summary, stringsAsFactors = FALSE)",
  "    p <- ggplot2::ggplot(dfp, ggplot2::aes(x = transition, y = n_edges)) + ggplot2::geom_col() + ggplot2::theme_minimal() +",
  "      ggplot2::labs(title='Edge counts by transition (interactive)', x='Transition', y='Number of edges') +",
  "      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))",
  "    plotly::ggplotly(p)",
  "    TRUE",
  "  }, error = function(e) {",
  "    cat('Plotly rendering failed; falling back to static plot.\\n')",
  "    FALSE",
  "  })",
  "}",
  "if (!ok) {",
  "  if (file.exists(manifest$output_files$Plot_EdgesByTransition)) knitr::include_graphics(manifest$output_files$Plot_EdgesByTransition)",
  "}",
  "```",
  "",
  "## Results tables (first 30 rows each)",
  "",
  "### Transition Summary",
  "```{r}",
  paste0("knitr::kable(utils::head(read.csv('", gsub("\\\\", "/", paths$Transition_Summary), "', stringsAsFactors = FALSE), 30))"),
  "```",
  "",
  "### Function–Function Enrichment (Edge Counts)",
  "```{r}",
  paste0("knitr::kable(utils::head(read.csv('", gsub("\\\\", "/", paths$Pair_nEdges), "', stringsAsFactors = FALSE), 30))"),
  "```",
  "",
  "### Function–Function Enrichment (Σ|Δr|)",
  "```{r}",
  paste0("knitr::kable(utils::head(read.csv('", gsub("\\\\", "/", paths$Pair_sumAbs), "', stringsAsFactors = FALSE), 30))"),
  "```",
  "",
  "### Sign Balance",
  "```{r}",
  paste0("knitr::kable(utils::head(read.csv('", gsub("\\\\", "/", paths$Pair_Sign), "', stringsAsFactors = FALSE), 30))"),
  "```",
  "",
  "## Interpretation notes",
  "",
  "- **High positive z for edge counts** indicates a function–function pair occurs more often than expected under the permutation null, consistent with transition-specific functional coupling.",
  "- **High positive z for Σ|Δr|** indicates not only more edges but also a larger aggregate magnitude of correlation reorganization concentrated in that function–function interface.",
  "- **Sign balance** (fraction positive/negative; direction_bias) helps distinguish coordinated strengthening vs weakening of associations in that interface.",
  "- **Entropy summaries** quantify concentration of reorganization across pairs; lower entropy implies that turnover is focused into fewer functional interfaces.",
  "",
  "## Reproducibility",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, qmd_path)

# Render behavior (MANDATORY): use rmarkdown::render as final step
render_ok <- TRUE
render_err <- NULL
tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_format = "html_document",
    output_file = basename(html_path),
    output_dir = output_dir,
    quiet = TRUE,
    params = list(
      run_id = run_id,
      run_timestamp = as.character(Sys.time()),
      script_name = script_name,
      script_path = script_path,
      script_full = ifelse(is.na(script_full), "NA", script_full),
      input_file = normalizePath(in_path, winslash = "/", mustWork = TRUE),
      runtime_seconds = runtime_seconds,
      enable_plotly = enable_plotly
    )
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_err <<- e$message
})

# Refresh inventory AFTER rendering (MANDATORY)
inventory_path <- write_file_inventory(output_dir, filename = "Project_Manifest_Files.csv")

# Update manifest with final report paths + refreshed inventory reference
manifest$report <- list(
  qmd = normalizePath(qmd_path, winslash = "/", mustWork = FALSE),
  html = normalizePath(html_path, winslash = "/", mustWork = FALSE),
  render_ok = render_ok,
  render_error = if (!render_ok) render_err else NA_character_
)
manifest$generated_files_inventory <- normalizePath(inventory_path, winslash = "/", mustWork = FALSE)
manifest$parameters$runtime_seconds <- runtime_seconds
manifest_path <- write_manifest_json(manifest, output_dir, filename = "Project_Manifest.json")

# Console summary (MANDATORY)
cat("\nWrote outputs under:\n  ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Key outputs:\n")
cat("  - ", paths$Transition_Summary, "\n", sep = "")
cat("  - ", paths$Pair_nEdges, "\n", sep = "")
cat("  - ", paths$Pair_sumAbs, "\n", sep = "")
cat("  - ", paths$Pair_Sign, "\n", sep = "")
cat("  - ", paths$Diagnostics, "\n", sep = "")
cat("  - ", manifest_path, "\n", sep = "")
cat("  - ", inventory_path, "\n", sep = "")
cat("  - ", qmd_path, "\n", sep = "")
if (render_ok) {
  cat("\nHTML report created:\n  ", normalizePath(html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("\nERROR: HTML report render failed.\nReason: ", render_err, "\nQMD written at:\n  ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
}

cat("\nDone.\n")
