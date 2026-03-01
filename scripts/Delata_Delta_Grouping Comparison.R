#!/usr/bin/env Rscript
# ============================================================
# Δ–Δ Grouping Comparison (SOURCE-button-safe, Base R)
#   - No dplyr/tidyr/optparse
#   - Mac-friendly pickers (rstudioapi) with fallbacks
#   - Standardized output: ./outputs/<run_name>_<timestamp>/
#   - Manifest + inventory + Quarto report (WD-safe, fixed render path)
#
# REQUIRED INPUT:
#   A CSV with columns representing:
#     - Day (numeric/integer)
#     - Gene ID (e.g., Source)
#     - daily_value (per-day aggregated expression)
#
# OPTIONAL INPUT:
#   Gene2Pathway mapping CSV with columns:
#     gene, pathway
# ============================================================

# ------------------------------
# 0) Logging
# ------------------------------
msg <- function(...) { cat(paste0(..., collapse=""), "\n"); flush.console() }

msg("[START] Δ–Δ Grouping Comparison script is running.")
msg("[INFO] interactive() = ", interactive())

# ------------------------------
# 1) Optional packages (no hard dependencies)
# ------------------------------
have_rstudioapi <- requireNamespace("rstudioapi", quietly = TRUE) && interactive() && rstudioapi::isAvailable()
have_igraph     <- requireNamespace("igraph", quietly = TRUE)
have_quarto     <- requireNamespace("quarto", quietly = TRUE)

# ------------------------------
# 2) Reproducibility seed
# ------------------------------
SET_SEED <- 1
set.seed(SET_SEED)

# ------------------------------
# 3) Small utilities
# ------------------------------
timestamp_string <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
norm_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
safe_dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

pick_file <- function(caption = "Select a file") {
  if (have_rstudioapi) {
    p <- rstudioapi::selectFile(caption = caption)
    if (!is.null(p) && nzchar(p)) return(p)
  }
  p <- tryCatch(file.choose(), error = function(e) "")
  if (nzchar(p)) return(p)
  readline(paste0(caption, " (paste full path): "))
}

choose_from_list <- function(options, prompt, default_idx = 1) {
  msg("")
  msg(prompt)
  for (i in seq_along(options)) msg("  [", i, "] ", options[i])
  ans <- readline(paste0("Choose [", default_idx, "]: "))
  if (!nzchar(ans)) return(options[default_idx])
  if (grepl("^[0-9]+$", ans)) {
    k <- as.integer(ans)
    if (is.finite(k) && k >= 1 && k <= length(options)) return(options[k])
  }
  if (ans %in% options) return(ans)
  stop("Invalid selection.")
}

choose_column <- function(df, prompt, default = NULL) {
  cols <- colnames(df)
  msg("")
  msg(prompt)
  for (i in seq_along(cols)) msg("  [", i, "] ", cols[i])
  if (!is.null(default) && default %in% cols) msg("Press Enter for default: ", default)
  ans <- readline("Choose column by number or name: ")
  if (!nzchar(ans) && !is.null(default) && default %in% cols) return(default)
  if (grepl("^[0-9]+$", ans)) {
    idx <- as.integer(ans)
    if (is.finite(idx) && idx >= 1 && idx <= length(cols)) return(cols[idx])
  }
  if (ans %in% cols) return(ans)
  stop("Invalid column choice.")
}

# Script identity (best-effort)
get_script_identity <- function() {
  script_full <- ""
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) script_full <- sub("^--file=", "", file_arg)
  script_full <- if (nzchar(script_full)) norm_path(script_full) else ""
  list(
    script_name = if (nzchar(script_full)) basename(script_full) else "",
    script_path = if (nzchar(script_full)) dirname(script_full) else "",
    script_full = script_full
  )
}

# Simple JSON writer (no jsonlite)
to_json_safe <- function(x, indent = 0) {
  sp <- paste(rep("  ", indent), collapse = "")
  if (is.null(x)) return("null")
  if (is.logical(x)) return(ifelse(isTRUE(x), "true", "false"))
  if (is.numeric(x)) return(ifelse(is.finite(x), format(x, scientific = FALSE), "null"))
  if (is.character(x)) {
    s <- gsub("\\\\", "\\\\\\\\", x)
    s <- gsub("\"", "\\\\\"", s)
    return(paste0("\"", s, "\""))
  }
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms) && any(nzchar(nms))) {
      parts <- character(0)
      for (nm in nms) {
        parts <- c(parts, paste0(sp, "  ", to_json_safe(nm), ": ", to_json_safe(x[[nm]], indent + 1)))
      }
      return(paste0("{\n", paste(parts, collapse = ",\n"), "\n", sp, "}"))
    } else {
      parts <- vapply(x, function(xx) to_json_safe(xx, indent + 1), character(1))
      return(paste0("[", paste(parts, collapse = ", "), "]"))
    }
  }
  "\"\""
}

# ------------------------------
# 4) Output directory standardization
# ------------------------------
setup_output_dir <- function() {
  safe_dir_create("outputs")
  run_name <- readline("Enter a run name [DeltaDelta_Grouping]: ")
  if (!nzchar(run_name)) run_name <- "DeltaDelta_Grouping"
  ts <- timestamp_string()
  out_dir <- file.path("outputs", paste0(run_name, "_", ts))
  safe_dir_create(out_dir)
  safe_dir_create(file.path(out_dir, "Tables"))
  safe_dir_create(file.path(out_dir, "Figures"))
  safe_dir_create(file.path(out_dir, "Manifest"))
  safe_dir_create(file.path(out_dir, "Report"))
  list(out_dir = norm_path(out_dir), run_name = run_name, timestamp = ts)
}

# ------------------------------
# 5) Read input and build Δ profiles (base R)
# ------------------------------
read_csv_strict <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

build_delta_profiles <- function(df, day_col, gene_col, value_col) {
  day  <- as.integer(df[[day_col]])
  gene <- as.character(df[[gene_col]])
  val  <- as.numeric(df[[value_col]])
  
  if (any(!is.finite(day))) stop("Day column has non-numeric values after coercion.")
  if (any(!nzchar(gene))) stop("Gene column contains empty values.")
  if (all(is.na(val))) stop("Value column is all NA after coercion.")
  
  days  <- sort(unique(day))
  genes <- sort(unique(gene))
  
  X <- matrix(NA_real_, nrow = length(genes), ncol = length(days),
              dimnames = list(genes, as.character(days)))
  
  # if multiple rows per (gene, day), take mean
  key <- paste(gene, day, sep = "||")
  agg <- tapply(val, key, mean, na.rm = TRUE)
  keys <- names(agg)
  spl <- strsplit(keys, "\\|\\|")
  
  for (i in seq_along(spl)) {
    g <- spl[[i]][1]
    d <- spl[[i]][2]
    X[g, d] <- as.numeric(agg[[i]])
  }
  
  if (ncol(X) < 2) stop("Need at least 2 unique days to compute Δ.")
  
  D <- t(apply(X, 1, function(v) diff(v)))
  rownames(D) <- rownames(X)
  colnames(D) <- paste0("d", days[-1], "_minus_d", days[-length(days)])
  
  list(days = days, X = X, D = D)
}

filter_delta_matrix <- function(D, min_transitions = 3) {
  ok <- apply(D, 1, function(v) {
    v2 <- v[is.finite(v)]
    length(v2) >= min_transitions && stats::sd(v2) > 0
  })
  D[ok, , drop = FALSE]
}

# ------------------------------
# 6) Methods
# ------------------------------
# A) Change-point timing groups = transition with max |Δ|
method_timing_groups <- function(D) {
  idx <- apply(D, 1, function(v) {
    v2 <- v; v2[!is.finite(v2)] <- NA_real_
    if (all(is.na(v2))) return(NA_integer_)
    which.max(abs(v2))
  })
  timing_transition <- rep(NA_character_, length(idx))
  timing_transition[!is.na(idx)] <- colnames(D)[idx[!is.na(idx)]]
  
  max_abs_delta <- apply(D, 1, function(v) {
    v2 <- v; v2[!is.finite(v2)] <- NA_real_
    if (all(is.na(v2))) return(NA_real_)
    max(abs(v2), na.rm = TRUE)
  })
  
  data.frame(
    Source = rownames(D),
    timing_transition = timing_transition,
    max_abs_delta = max_abs_delta,
    stringsAsFactors = FALSE
  )
}

# B) Sign-pattern groups (+/−/0)
method_sign_patterns <- function(D, eps_flat = 0) {
  sign_code <- function(v) {
    v2 <- v; v2[!is.finite(v2)] <- NA_real_
    s <- rep("N", length(v2))
    s[is.finite(v2) & v2 >  eps_flat] <- "+"
    s[is.finite(v2) & v2 < -eps_flat] <- "-"
    s[is.finite(v2) & abs(v2) <= eps_flat] <- "0"
    paste(s, collapse = "")
  }
  pat <- apply(D, 1, sign_code)
  data.frame(Source = rownames(D), sign_pattern = pat, stringsAsFactors = FALSE)
}

safe_cor <- function(D, method = "spearman") {
  C <- suppressWarnings(stats::cor(t(D), method = method, use = "pairwise.complete.obs"))
  C[!is.finite(C)] <- 0
  C <- pmax(pmin(C, 1), -1)
  C
}

# C) Network communities (igraph optional)
method_network_communities <- function(D, cor_method = "spearman",
                                       build = c("topk", "thr"),
                                       topk = 10,
                                       abs_thr = 0.85,
                                       use_abs = TRUE) {
  build <- match.arg(build)
  assignments <- data.frame(Source = rownames(D), community = NA_integer_, stringsAsFactors = FALSE)
  edges_out   <- data.frame(from = character(0), to = character(0), cor = numeric(0), stringsAsFactors = FALSE)
  
  if (!have_igraph) {
    return(list(assignments = assignments, edges = edges_out, note = "igraph not installed; skipped."))
  }
  
  C <- safe_cor(D, method = cor_method)
  diag(C) <- 0
  genes <- rownames(D)
  n <- nrow(C)
  
  edges <- list()
  
  if (build == "topk") {
    for (i in seq_len(n)) {
      r <- C[i, ]
      rr <- if (use_abs) abs(r) else r
      ord <- order(rr, decreasing = TRUE, na.last = NA)
      ord <- ord[ord != i]
      k <- min(topk, length(ord))
      if (k < 1) next
      nbrs <- ord[seq_len(k)]
      for (j in nbrs) edges[[length(edges) + 1]] <- c(genes[i], genes[j], r[j])
    }
  } else {
    for (i in seq_len(n)) {
      if (i == n) next
      for (j in (i + 1):n) {
        r <- C[i, j]
        rr <- if (use_abs) abs(r) else r
        if (is.finite(rr) && rr >= abs_thr) edges[[length(edges) + 1]] <- c(genes[i], genes[j], r)
      }
    }
  }
  
  if (length(edges) == 0) {
    return(list(assignments = assignments, edges = edges_out,
                note = "No edges under current settings; increase TOPK or lower threshold."))
  }
  
  E <- do.call(rbind, edges)
  edges_out <- data.frame(from = E[,1], to = E[,2], cor = as.numeric(E[,3]), stringsAsFactors = FALSE)
  
  g <- igraph::graph_from_data_frame(
    d = data.frame(from = edges_out$from, to = edges_out$to, weight = abs(edges_out$cor),
                   stringsAsFactors = FALSE),
    directed = FALSE,
    vertices = data.frame(name = genes, stringsAsFactors = FALSE)
  )
  
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  memb <- igraph::membership(cl)
  
  assignments <- data.frame(Source = names(memb), community = as.integer(memb), stringsAsFactors = FALSE)
  list(assignments = assignments, edges = edges_out, note = "OK")
}

mean_upper_tri <- function(M) {
  if (nrow(M) < 2) return(NA_real_)
  idx <- upper.tri(M, diag = FALSE)
  mean(M[idx], na.rm = TRUE)
}

# D) Pathway coherence (requires gene,pathway mapping)
method_pathway_coherence <- function(D, map_df, cor_method = "spearman",
                                     min_size = 5, max_size = 500,
                                     null_B = 2000, seed = 1) {
  if (!all(c("gene", "pathway") %in% colnames(map_df))) stop("Mapping must have columns: gene, pathway")
  
  set.seed(seed)
  map_df$gene <- as.character(map_df$gene)
  map_df$pathway <- as.character(map_df$pathway)
  
  all_genes <- rownames(D)
  map_df <- map_df[map_df$gene %in% all_genes, , drop = FALSE]
  pathways <- sort(unique(map_df$pathway))
  
  res <- list()
  
  for (pw in pathways) {
    genes <- unique(map_df$gene[map_df$pathway == pw])
    n <- length(genes)
    if (n < min_size || n > max_size) next
    
    Cobs <- safe_cor(D[genes, , drop = FALSE], method = cor_method)
    diag(Cobs) <- NA
    obs <- mean_upper_tri(abs(Cobs))
    
    null <- numeric(null_B)
    for (b in seq_len(null_B)) {
      gset <- sample(all_genes, size = n, replace = FALSE)
      Cn <- safe_cor(D[gset, , drop = FALSE], method = cor_method)
      diag(Cn) <- NA
      null[b] <- mean_upper_tri(abs(Cn))
    }
    
    p_hi <- mean(null >= obs, na.rm = TRUE)
    z <- (obs - mean(null, na.rm = TRUE)) / stats::sd(null, na.rm = TRUE)
    
    res[[length(res) + 1]] <- data.frame(
      pathway = pw,
      n_genes = n,
      coherence_obs = obs,
      null_mean = mean(null, na.rm = TRUE),
      null_sd = stats::sd(null, na.rm = TRUE),
      z = z,
      p_empirical_hi = p_hi,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(res) == 0) return(data.frame())
  out <- do.call(rbind, res)
  out$q_empirical_hi <- p.adjust(out$p_empirical_hi, method = "BH")
  out <- out[order(out$q_empirical_hi, -out$z), , drop = FALSE]
  out
}

write_contingency <- function(df, colA, colB, outfile) {
  tab <- table(df[[colA]], df[[colB]], useNA = "ifany")
  m <- as.matrix(tab)
  out <- cbind(GroupA = rownames(m), as.data.frame(m, check.names = FALSE))
  write.csv(out, outfile, row.names = FALSE)
}

# ------------------------------
# 7) Quarto report (WD-safe, FIXED render path)
# ------------------------------
build_and_render_report <- function(output_dir,
                                    report_dir,
                                    run_name,
                                    timestamp,
                                    input_file,
                                    pathway_file,
                                    manifest_json,
                                    deps,
                                    script_header = character(0)) {
  
  if (!dir.exists(output_dir)) stop("output_dir does not exist: ", output_dir)
  if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  
  qmd_name <- paste0("DeltaDelta_Grouping_Report_", timestamp, ".qmd")
  qmd_path <- file.path(report_dir, qmd_name)
  
  esc <- function(x) gsub("\"", "\\\\\"", x)
  
  qmd_lines <- c(
    "---",
    paste0("title: \"Δ–Δ Grouping Comparison: ", esc(run_name), "\""),
    paste0("date: \"", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\""),
    "format: html",
    "---",
    "",
    "# Summary",
    "",
    "This report compares four ways of grouping Δ–Δ (day-to-day change) profiles across development:",
    "",
    "- **Change-point timing groups**: transition with maximal |Δ| per gene",
    "- **Sign-pattern groups**: +/−/0 strings across transitions",
    "- **Network communities**: community detection on a Δ–Δ similarity graph (if `igraph` available)",
    "- **Pathway-first coherence**: coherence within predefined pathways vs random sets (if mapping provided)",
    "",
    "# Script header (verbatim; best-effort)",
    "",
    "```",
    if (length(script_header) > 0) paste(script_header, collapse = "\n") else "(script header not available; likely sourced via RStudio)",
    "```",
    "",
    "# Metadata",
    "",
    "```{r}",
    paste0("output_dir <- \"", esc(norm_path(output_dir)), "\""),
    paste0("input_file <- \"", esc(norm_path(input_file)), "\""),
    paste0("pathway_file <- \"", if (nzchar(pathway_file) && file.exists(pathway_file)) esc(norm_path(pathway_file)) else "", "\""),
    paste0("manifest_json <- \"", esc(norm_path(manifest_json)), "\""),
    paste0("timestamp <- \"", esc(timestamp), "\""),
    "cat(paste0('Output dir: ', output_dir, '\\n'))",
    "cat(paste0('Input file: ', input_file, '\\n'))",
    "cat(paste0('Pathway file: ', pathway_file, '\\n'))",
    "cat(paste0('Manifest: ', manifest_json, '\\n'))",
    "cat(paste0('Timestamp: ', timestamp, '\\n'))",
    "```",
    "",
    "# Parameters (Manifest.json)",
    "",
    "```{r}",
    "if (file.exists(manifest_json)) {",
    "  cat(readLines(manifest_json), sep='\\n')",
    "} else {",
    "  cat('Manifest.json not found.\\n')",
    "}",
    "```",
    "",
    "# Dependencies",
    "",
    paste0("- ", paste(deps, collapse = "\n- ")),
    "",
    "# Analytical logic",
    "",
    "We define the Δ profile for gene *g* as:",
    "",
    "$$\\Delta_g(t_i) = x_g(t_{i+1}) - x_g(t_i)$$",
    "",
    "where $x_g(t)$ is the per-day aggregated expression and $t_i$ are the sampled days.",
    "",
    "# Generated outputs",
    "",
    "```{r}",
    "files <- list.files(output_dir, recursive = TRUE)",
    "cat(paste(files, collapse='\\n'))",
    "```",
    "",
    "# Figures",
    "",
    "```{r}",
    "fig_dir <- file.path(output_dir, 'Figures')",
    "imgs <- list.files(fig_dir, pattern='\\\\.(png|pdf)$', full.names=TRUE, ignore.case=TRUE)",
    "if (length(imgs)==0) {",
    "  cat('No figures found.\\n')",
    "} else {",
    "  for (im in imgs) cat('![](', im, ')\\n\\n', sep='')",
    "}",
    "```",
    "",
    "# Session info",
    "",
    "```{r}",
    "sessionInfo()",
    "```"
  )
  
  writeLines(qmd_lines, con = qmd_path, useBytes = TRUE)
  
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  # Critical: render by relative path from output_dir
  qmd_rel <- file.path("Report", basename(qmd_path))
  
  msg("[OK] QMD written: ", norm_path(qmd_path))
  msg("[OK] Rendering via relative path from output_dir: ", qmd_rel)
  
  if (!file.exists(qmd_rel)) stop("QMD missing at expected location: ", file.path(output_dir, qmd_rel))
  
  if (!have_quarto) {
    msg("[WARN] quarto package not available; skipping report render.")
    return(list(rendered = FALSE, qmd_path = norm_path(qmd_path), html_path = norm_path(html_path),
                note = "quarto not available"))
  }
  
  ok <- TRUE
  note <- "OK"
  tryCatch({
    quarto::quarto_render(qmd_rel)
  }, error = function(e) {
    ok <<- FALSE
    note <<- paste0("Quarto render failed: ", conditionMessage(e))
  })
  
  # Fallback HTML if needed
  if (!ok || !file.exists(html_path)) {
    msg("[WARN] Quarto did not produce HTML; writing fallback HTML.")
    fallback <- c(
      "<html><head><meta charset='utf-8'><title>Δ–Δ Grouping Report (Fallback)</title></head><body>",
      paste0("<h1>Δ–Δ Grouping Comparison: ", run_name, "</h1>"),
      paste0("<p><b>Timestamp:</b> ", timestamp, "</p>"),
      paste0("<p><b>Output dir:</b> ", norm_path(output_dir), "</p>"),
      "<h2>Status</h2>",
      paste0("<pre>", note, "</pre>"),
      "<h2>Generated outputs</h2>",
      "<pre>",
      paste(list.files(output_dir, recursive = TRUE), collapse = "\n"),
      "</pre>",
      "</body></html>"
    )
    writeLines(fallback, con = html_path, useBytes = TRUE)
    return(list(rendered = FALSE, qmd_path = norm_path(qmd_path), html_path = norm_path(html_path),
                note = note))
  }
  
  list(rendered = TRUE, qmd_path = norm_path(qmd_path), html_path = norm_path(html_path), note = "OK")
}

# ------------------------------
# 8) MAIN workflow (prompts)
# ------------------------------
msg("[STEP] Select input file (DailySummary_ByDayGene.csv)...")
input_file <- norm_path(pick_file("Select DailySummary_ByDayGene.csv"))
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

df <- read_csv_strict(input_file)

day_col  <- choose_column(df, "Choose the Day column:", default = "Day")
gene_col <- choose_column(df, "Choose the gene ID column:", default = "Source")
val_col  <- choose_column(df, "Choose the value column:", default = "daily_value")

use_pathway <- tolower(trimws(readline("Run pathway-first coherence? (y/n) [n]: ")))
if (!nzchar(use_pathway)) use_pathway <- "n"

pathway_file <- ""
if (use_pathway == "y") {
  pathway_file <- norm_path(pick_file("Select Gene2Pathway.csv (columns: gene,pathway)"))
  if (!file.exists(pathway_file)) {
    msg("[WARN] Mapping file not found; skipping pathway coherence.")
    pathway_file <- ""
  }
}

eps_flat <- readline("Sign-pattern flat threshold EPS (abs(Δ) <= EPS => 0) [0]: ")
eps_flat <- if (nzchar(eps_flat)) as.numeric(eps_flat) else 0
if (!is.finite(eps_flat)) eps_flat <- 0

min_trans <- readline("Minimum finite transitions per gene to include [3]: ")
min_trans <- if (nzchar(min_trans)) as.integer(min_trans) else 3
if (!is.finite(min_trans) || min_trans < 2) min_trans <- 3

cor_method <- tolower(trimws(readline("Correlation method: spearman or pearson [spearman]: ")))
if (!nzchar(cor_method)) cor_method <- "spearman"
if (!cor_method %in% c("spearman", "pearson")) cor_method <- "spearman"

net_build <- tolower(trimws(readline("Network build: topk or thr [topk]: ")))
if (!nzchar(net_build)) net_build <- "topk"
if (!net_build %in% c("topk", "thr")) net_build <- "topk"

net_topk <- readline("Network TOPK neighbors (if topk) [10]: ")
net_topk <- if (nzchar(net_topk)) as.integer(net_topk) else 10
if (!is.finite(net_topk) || net_topk < 1) net_topk <- 10

net_thr <- readline("Network abs correlation threshold (if thr) [0.85]: ")
net_thr <- if (nzchar(net_thr)) as.numeric(net_thr) else 0.85
if (!is.finite(net_thr)) net_thr <- 0.85

path_nullB <- readline("Pathway coherence null permutations per pathway [2000]: ")
path_nullB <- if (nzchar(path_nullB)) as.integer(path_nullB) else 2000
if (!is.finite(path_nullB) || path_nullB < 100) path_nullB <- 2000

# Output setup
od <- setup_output_dir()
output_dir <- od$out_dir
tables_dir <- file.path(output_dir, "Tables")
fig_dir    <- file.path(output_dir, "Figures")
man_dir    <- file.path(output_dir, "Manifest")
rep_dir    <- file.path(output_dir, "Report")

msg("[OK] Output directory: ", output_dir)

# ------------------------------
# Build Δ profiles
# ------------------------------
msg("[STEP] Building Δ profiles...")
built <- build_delta_profiles(df, day_col, gene_col, val_col)
D <- built$D
D <- filter_delta_matrix(D, min_transitions = min_trans)
msg("[OK] Genes retained: ", nrow(D), " | Transitions: ", ncol(D))

delta_csv <- file.path(tables_dir, "DeltaProfiles.csv")
write.csv(data.frame(Source = rownames(D), D, check.names = FALSE), delta_csv, row.names = FALSE)
msg("[OK] Wrote: ", norm_path(delta_csv))

# ------------------------------
# Run methods
# ------------------------------
msg("[STEP] Method A: timing groups...")
timing_df <- method_timing_groups(D)
timing_out <- file.path(tables_dir, "Method_ChangePoint_TimingGroups.csv")
write.csv(timing_df, timing_out, row.names = FALSE)

msg("[STEP] Method B: sign patterns...")
sign_df <- method_sign_patterns(D, eps_flat = eps_flat)
sign_out <- file.path(tables_dir, "Method_SignPattern_Groups.csv")
write.csv(sign_df, sign_out, row.names = FALSE)

msg("[STEP] Method C: network communities...")
net <- method_network_communities(D, cor_method = cor_method, build = net_build,
                                  topk = net_topk, abs_thr = net_thr, use_abs = TRUE)
net_out <- file.path(tables_dir, "Method_NetworkCommunities.csv")
write.csv(net$assignments, net_out, row.names = FALSE)
edge_out <- file.path(tables_dir, "Network_Edges.csv")
write.csv(net$edges, edge_out, row.names = FALSE)
msg("[OK] Network: ", net$note)

msg("[STEP] Method D: pathway coherence (optional)...")
path_out <- file.path(tables_dir, "Method_PathwayCoherence.csv")
path_df <- data.frame()
path_note <- "SKIPPED"
if (nzchar(pathway_file) && file.exists(pathway_file)) {
  map_df <- read_csv_strict(pathway_file)
  path_df <- method_pathway_coherence(D, map_df, cor_method = cor_method,
                                      min_size = 5, max_size = 500,
                                      null_B = path_nullB, seed = SET_SEED)
  path_note <- "OK"
}
write.csv(path_df, path_out, row.names = FALSE)
msg("[OK] Pathway coherence: ", path_note)

# ------------------------------
# Merge assignments + contingency comparisons
# ------------------------------
msg("[STEP] Comparing methods (assignments + contingency tables)...")
assign <- data.frame(Source = rownames(D), stringsAsFactors = FALSE)

m1 <- match(assign$Source, timing_df$Source)
assign$timing_group <- timing_df$timing_transition[m1]

m2 <- match(assign$Source, sign_df$Source)
assign$sign_group <- sign_df$sign_pattern[m2]

m3 <- match(assign$Source, net$assignments$Source)
assign$network_community <- net$assignments$community[m3]

assign_out <- file.path(tables_dir, "AllMethods_Assignments.csv")
write.csv(assign, assign_out, row.names = FALSE)

write_contingency(assign, "timing_group", "sign_group",
                  file.path(tables_dir, "Compare_Contingency_Timing_vs_Sign.csv"))
write_contingency(assign, "timing_group", "network_community",
                  file.path(tables_dir, "Compare_Contingency_Timing_vs_Network.csv"))
write_contingency(assign, "sign_group", "network_community",
                  file.path(tables_dir, "Compare_Contingency_Sign_vs_Network.csv"))

# ------------------------------
# Figures (simple, deterministic)
# ------------------------------
msg("[STEP] Writing summary figures...")
tim_tab <- sort(table(timing_df$timing_transition), decreasing = TRUE)
png(file.path(fig_dir, "TimingGroup_Counts.png"), width = 1600, height = 1000, res = 150)
par(mar=c(10,5,3,2))
barplot(tim_tab, las=2, main="Change-point timing groups (max |Δ| transition)", ylab="Number of genes")
dev.off()

pat_tab <- sort(table(sign_df$sign_pattern), decreasing = TRUE)
topN <- min(25, length(pat_tab))
png(file.path(fig_dir, "SignPattern_Top25.png"), width = 1600, height = 1000, res = 150)
par(mar=c(10,5,3,2))
barplot(pat_tab[seq_len(topN)], las=2, main="Top sign-pattern groups (top 25)", ylab="Number of genes")
dev.off()

com_tab <- sort(table(net$assignments$community), decreasing = TRUE)
png(file.path(fig_dir, "NetworkCommunity_Sizes.png"), width = 1600, height = 1000, res = 150)
par(mar=c(10,5,3,2))
barplot(com_tab, las=2, main="Network community sizes (Δ similarity graph)", ylab="Number of genes")
dev.off()

if (nrow(path_df) > 0) {
  showN <- min(25, nrow(path_df))
  png(file.path(fig_dir, "PathwayCoherence_Top25.png"), width = 1600, height = 1000, res = 150)
  par(mar=c(10,5,3,2))
  barplot(path_df$z[seq_len(showN)], names.arg = path_df$pathway[seq_len(showN)],
          las=2, main="Top pathway coherence (Z vs random sets)", ylab="Z score")
  dev.off()
}

# ------------------------------
# Manifest + inventory
# ------------------------------
msg("[STEP] Writing manifest + inventory...")
ident <- get_script_identity()

inventory_path <- file.path(man_dir, "Output_Inventory.txt")
inv <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
writeLines(norm_path(inv), con = inventory_path)

params <- list(
  seed = SET_SEED,
  day_col = day_col,
  gene_col = gene_col,
  value_col = val_col,
  eps_flat = eps_flat,
  min_transitions = min_trans,
  cor_method = cor_method,
  net_build = net_build,
  net_topk = net_topk,
  net_thr = net_thr,
  pathway_nullB = path_nullB,
  pathway_mapping_used = (nzchar(pathway_file) && file.exists(pathway_file)),
  pathway_status = path_note,
  igraph_available = have_igraph,
  quarto_available = have_quarto
)

manifest <- list(
  script_name = ident$script_name,
  script_path = ident$script_path,
  script_full = ident$script_full,
  input_paths = c(norm_path(input_file),
                  if (nzchar(pathway_file) && file.exists(pathway_file)) norm_path(pathway_file) else ""),
  output_dir = output_dir,
  timestamp = od$timestamp,
  run_name = od$run_name,
  parameters = params
)

manifest_json <- file.path(man_dir, "Manifest.json")
writeLines(to_json_safe(manifest, indent = 0), con = manifest_json)

msg("[OK] Manifest: ", norm_path(manifest_json))
msg("[OK] Inventory: ", norm_path(inventory_path))

# Dependencies list for report
deps <- c(
  "base R",
  if (have_rstudioapi) "rstudioapi (available)" else "rstudioapi (not available)",
  if (have_igraph) "igraph (available)" else "igraph (not available; network communities may be empty)",
  if (have_quarto) "quarto (available)" else "quarto (not available; HTML report will fall back)"
)

# Script header capture (best-effort; sourced scripts often can't get full header)
script_header <- c(
  "# Δ–Δ Grouping Comparison (Base R)",
  paste0("# Run: ", od$run_name, " | Timestamp: ", od$timestamp),
  paste0("# Input: ", input_file)
)

# ------------------------------
# Quarto report
# ------------------------------
msg("[STEP] Generating Quarto report (if quarto available)...")
qr <- build_and_render_report(
  output_dir = output_dir,
  report_dir = rep_dir,
  run_name = od$run_name,
  timestamp = od$timestamp,
  input_file = input_file,
  pathway_file = pathway_file,
  manifest_json = manifest_json,
  deps = deps,
  script_header = script_header
)

# ------------------------------
# Final summary
# ------------------------------
msg("")
msg("[DONE] Completed.")
msg("Output directory: ", output_dir)
msg("Tables:          ", norm_path(tables_dir))
msg("Figures:         ", norm_path(fig_dir))
msg("Manifest:        ", norm_path(man_dir))
if (isTRUE(qr$rendered)) {
  msg("Report QMD:      ", qr$qmd_path)
  msg("Report HTML:     ", qr$html_path)
} else {
  msg("Report HTML:     ", qr$html_path)
  msg("Report note:     ", qr$note)
}
