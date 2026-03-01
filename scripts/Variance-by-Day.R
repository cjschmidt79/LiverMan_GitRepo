###############################################################################
# Variance-by-Day + Cluster “Variance-of-Variance” Boxplots (WIDE input)
# SOURCE-to-run (Biologist-first, no dplyr)
# ---------------------------------------------------------------------------
# INPUT:
#   - WIDE expression table:
#       rows = biological replicates (samples)
#       columns = genes (numeric expression)
#   - Sample IDs encode Day and Rep as "#_#" (e.g., "14_3")
#       Sample ID can be:
#         (A) rownames, OR
#         (B) first column (non-numeric)
#   - OPTIONAL gene->cluster table:
#       columns: gene (gene name), cluster (numeric ID)
#
# OUTPUT:
#   - Per-gene per-day variance stats (tables): variance, SD, mean, CV
#   - Optional mean-corrected variability: residual log2-variance vs log2-mean per day
#   - Added columns for plotting:
#       log2Var = log2(Variance + pseudocount)
#       log2ResVar (alias of residuals if enabled)
#   - Per-cluster figures:
#       (i) boxplots by day of per-gene variability (raw or log2Var)
#       (ii) optional boxplots by day of residual variability (log2 residuals)
#       (iii) mean/median trajectory plots across genes in cluster vs day
#       (iv) top-N gene variance trajectories (default top 10) per cluster
#   - NEW: Explicit "variance-of-variance" summary columns per cluster/day:
#       var_metric = var(metric across genes within cluster at that day)
#       sd_metric  = sd(metric across genes within cluster at that day)
#     plus ALL-cluster trajectory plots of var_metric and sd_metric.
#
# BIologist-first requirements met:
#   - RStudio file pickers when available
#   - outputs/<run_name>_<input_stub>_<timestamp>/ standardized
#   - provenance manifest + inventory
#   - Quarto HTML report via canonical contract (setwd(output_dir) during render)
#   - AUTO-RUN when sourced in interactive RStudio
###############################################################################

# =========================
# Libraries (minimal)
# =========================
suppressWarnings({
  suppressMessages({
    has_rstudioapi <- requireNamespace("rstudioapi", quietly = TRUE)
    has_ggplot2    <- requireNamespace("ggplot2", quietly = TRUE)
    has_quarto     <- requireNamespace("quarto", quietly = TRUE)
    has_knitr      <- requireNamespace("knitr", quietly = TRUE)
  })
})

# =========================
# Helpers: pickers, paths
# =========================
pick_file <- function(prompt = "Select a file") {
  if (has_rstudioapi && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::selectFile(caption = prompt), error = function(e) NA_character_)
    if (!is.na(p) && nzchar(p)) return(p)
  }
  file.choose()
}

norm_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

safe_mkdir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}

timestamp_string <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

# =========================
# Script identity capture (best-effort)
# =========================
capture_script_identity <- function() {
  script_full <- NA_character_
  script_path <- NA_character_
  script_name <- NA_character_
  
  if (has_rstudioapi && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) script_full <- ctx$path
  }
  
  if (is.na(script_full) || !nzchar(script_full)) {
    script_full <- tryCatch(norm_path(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
  }
  
  if (!is.na(script_full) && nzchar(script_full)) {
    script_path <- dirname(script_full)
    script_name <- basename(script_full)
  } else {
    script_name <- "VarianceByDay_Cluster_Boxplots.R"
    script_path <- norm_path(getwd())
    script_full <- NA_character_
  }
  
  list(script_name = script_name, script_path = script_path, script_full = script_full)
}

# =========================
# Header block capture (for Quarto)
# =========================
get_header_block <- function(script_full, fallback_name = "Script") {
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
    take_n <- min(length(lines), 220)
    head_lines <- lines[seq_len(take_n)]
    keep <- character()
    for (ln in head_lines) {
      if (grepl("^\\s*#|^\\s*$", ln)) keep <- c(keep, ln) else break
    }
    paste(keep, collapse = "\n")
  } else {
    paste0("# Header block unavailable (script_full not detected). Fallback: ", fallback_name)
  }
}

# =========================
# Robust CSV/TSV reader (base)
# =========================
read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# =========================
# Parse sample IDs: "#_#" -> Day, Rep
# Accepts: "14_3" or "Day14_3" or "D14_3"
# =========================
parse_day_rep <- function(sample_ids) {
  m <- regexec("^\\D*(\\d+)_(\\d+)\\D*$", sample_ids)
  reg <- regmatches(sample_ids, m)
  
  day <- rep(NA_integer_, length(sample_ids))
  repn <- rep(NA_integer_, length(sample_ids))
  
  for (i in seq_along(reg)) {
    if (length(reg[[i]]) == 3) {
      day[i]  <- as.integer(reg[[i]][2])
      repn[i] <- as.integer(reg[[i]][3])
    }
  }
  
  if (any(is.na(day)) || any(is.na(repn))) {
    bad <- sample_ids[is.na(day) | is.na(repn)]
    stop(
      "SampleID parsing failed for some rows.\n",
      "Expected '#_#' (e.g., '14_3'), optionally with prefix like 'Day14_3'.\n",
      "Examples of problematic IDs:\n  - ",
      paste(head(bad, 10), collapse = "\n  - ")
    )
  }
  
  data.frame(SampleID = sample_ids, Day = day, Rep = repn, stringsAsFactors = FALSE)
}

# =========================
# Compute per-day per-gene variance across replicates
# =========================
compute_gene_day_stats <- function(expr_mat, day_vec, pseudocount = 1e-8) {
  genes <- colnames(expr_mat)
  days  <- sort(unique(day_vec))
  
  out_list <- vector("list", length(days))
  names(out_list) <- as.character(days)
  
  for (d in days) {
    idx <- which(day_vec == d)
    sub <- expr_mat[idx, , drop = FALSE]
    
    g_mean <- colMeans(sub, na.rm = TRUE)
    g_var  <- apply(sub, 2, var, na.rm = TRUE)
    g_sd   <- sqrt(g_var)
    g_cv   <- g_sd / pmax(abs(g_mean), pseudocount)
    
    out_list[[as.character(d)]] <- data.frame(
      Day = d,
      Gene = genes,
      Mean = as.numeric(g_mean),
      Variance = as.numeric(g_var),
      SD = as.numeric(g_sd),
      CV = as.numeric(g_cv),
      n_samples_day = length(idx),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out_list)
}

# =========================
# Mean-corrected variability: residual log2(Var) ~ log2(Mean) per day
# Adds: log2Mean, log2Var, ResLog2Var, log2ResVar (alias)
# =========================
add_residual_log2var <- function(gene_day_long, pseudocount = 1e-8) {
  gene_day_long$log2Mean <- log2(gene_day_long$Mean + pseudocount)
  gene_day_long$log2Var  <- log2(gene_day_long$Variance + pseudocount)
  gene_day_long$ResLog2Var <- NA_real_
  gene_day_long$log2ResVar <- NA_real_
  
  days <- sort(unique(gene_day_long$Day))
  for (d in days) {
    idx <- which(gene_day_long$Day == d)
    df  <- gene_day_long[idx, , drop = FALSE]
    
    ok <- is.finite(df$log2Mean) & is.finite(df$log2Var)
    if (sum(ok) < 5) next
    
    fit <- lm(log2Var ~ log2Mean, data = df[ok, , drop = FALSE])
    r <- resid(fit)
    gene_day_long$ResLog2Var[idx][ok] <- r
    gene_day_long$log2ResVar[idx][ok] <- r
  }
  
  gene_day_long
}

# =========================
# Summaries for clusters by day (mean/median/IQR + NEW var/sd of metric across genes)
# =========================
summarize_cluster_day <- function(gene_day_long, metric_col = "Variance") {
  stopifnot(metric_col %in% colnames(gene_day_long))
  days <- sort(unique(gene_day_long$Day))
  clus <- sort(unique(gene_day_long$Cluster))
  
  rows <- list()
  k <- 1
  for (c in clus) {
    for (d in days) {
      idx <- which(gene_day_long$Cluster == c & gene_day_long$Day == d)
      vals <- gene_day_long[idx, metric_col]
      vals <- vals[is.finite(vals)]
      
      rows[[k]] <- data.frame(
        Cluster = c,
        Day = d,
        n_genes = sum(gene_day_long$Cluster == c & gene_day_long$Day == d),
        metric = metric_col,
        mean_metric = if (length(vals)) mean(vals) else NA_real_,
        median_metric = if (length(vals)) median(vals) else NA_real_,
        iqr_metric = if (length(vals)) IQR(vals) else NA_real_,
        # NEW: explicit "variance-of-variance" (dispersion across genes of per-gene variance metric)
        var_metric = if (length(vals) > 1) var(vals) else NA_real_,
        sd_metric  = if (length(vals) > 1) sd(vals) else NA_real_,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  do.call(rbind, rows)
}

# =========================
# Plot helpers
# =========================
save_boxplot_by_cluster <- function(df, out_png, metric_col, title, dpi = 600, w = 8, h = 5) {
  if (has_ggplot2) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(Day), y = .data[[metric_col]])) +
      ggplot2::geom_boxplot(outlier.size = 0.8) +
      ggplot2::labs(title = title, x = "Day", y = metric_col) +
      ggplot2::theme_bw()
    ggplot2::ggsave(out_png, p, dpi = dpi, width = w, height = h, units = "in")
  } else {
    png(out_png, width = w, height = h, units = "in", res = dpi)
    boxplot(df[[metric_col]] ~ factor(df$Day),
            main = title, xlab = "Day", ylab = metric_col, outline = TRUE)
    dev.off()
  }
}

save_trajplot_cluster_summary <- function(summary_df, out_png, metric_label, dpi = 600, w = 8, h = 5) {
  # summary_df: Cluster, Day, mean_metric, median_metric
  if (has_ggplot2) {
    p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = Day)) +
      ggplot2::geom_line(ggplot2::aes(y = mean_metric), linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = mean_metric), size = 1.6) +
      ggplot2::geom_line(ggplot2::aes(y = median_metric), linewidth = 0.8, linetype = 2) +
      ggplot2::geom_point(ggplot2::aes(y = median_metric), size = 1.6) +
      ggplot2::labs(
        title = paste0("Cluster ", unique(summary_df$Cluster), ": mean/median across genes vs Day"),
        x = "Day",
        y = metric_label
      ) +
      ggplot2::theme_bw()
    ggplot2::ggsave(out_png, p, dpi = dpi, width = w, height = h, units = "in")
  } else {
    png(out_png, width = w, height = h, units = "in", res = dpi)
    plot(summary_df$Day, summary_df$mean_metric, type = "b", pch = 16,
         xlab = "Day", ylab = metric_label,
         main = paste0("Cluster ", unique(summary_df$Cluster), ": mean/median vs Day"))
    lines(summary_df$Day, summary_df$median_metric, type = "b", pch = 1, lty = 2)
    legend("topright", legend = c("Mean", "Median"), lty = c(1, 2), pch = c(16, 1), bty = "n")
    dev.off()
  }
}

save_topN_gene_trajectories <- function(df_cluster, out_png, metric_col, top_genes, dpi = 600, w = 9, h = 6) {
  sub <- df_cluster[df_cluster$Gene %in% top_genes, , drop = FALSE]
  sub <- sub[order(sub$Gene, sub$Day), , drop = FALSE]
  
  if (has_ggplot2) {
    p <- ggplot2::ggplot(sub, ggplot2::aes(x = Day, y = .data[[metric_col]], group = Gene)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::labs(
        title = paste0("Top ", length(top_genes), " genes: ", metric_col, " vs Day (Cluster ", unique(sub$Cluster), ")"),
        x = "Day", y = metric_col
      ) +
      ggplot2::theme_bw()
    ggplot2::ggsave(out_png, p, dpi = dpi, width = w, height = h, units = "in")
  } else {
    png(out_png, width = w, height = h, units = "in", res = dpi)
    days <- sort(unique(sub$Day))
    mat <- matrix(NA_real_, nrow = length(days), ncol = length(top_genes),
                  dimnames = list(paste0("Day", days), top_genes))
    for (g in top_genes) {
      gg <- sub[sub$Gene == g, , drop = FALSE]
      mat[paste0("Day", gg$Day), g] <- gg[[metric_col]]
    }
    matplot(days, mat, type = "l", lty = 1, xlab = "Day", ylab = metric_col,
            main = paste0("Top ", length(top_genes), " genes: ", metric_col, " trajectories (Cluster ", unique(sub$Cluster), ")"))
    legend("topright", legend = top_genes, cex = 0.6, bty = "n")
    dev.off()
  }
}

# NEW: plot variance-of-variance trajectories (var_metric / sd_metric) across clusters
save_varofvar_trajectories_all <- function(summ_all, out_png, ycol = "var_metric",
                                           title = "Variance-of-variance across genes vs Day (by cluster)",
                                           dpi = 600, w = 10, h = 6) {
  if (!has_ggplot2) return(invisible(FALSE))
  if (!(ycol %in% colnames(summ_all))) return(invisible(FALSE))
  df <- summ_all[is.finite(summ_all[[ycol]]), , drop = FALSE]
  if (!nrow(df)) return(invisible(FALSE))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Day, y = .data[[ycol]])) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~Cluster, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = "Day", y = ycol)
  
  ggplot2::ggsave(out_png, p, dpi = dpi, width = w, height = h, units = "in")
  invisible(TRUE)
}

# =========================
# Manifest + inventory
# =========================
write_manifest <- function(out_dir, manifest_list) {
  json_path <- file.path(out_dir, "Project_Manifest.json")
  csv_path  <- file.path(out_dir, "Project_Manifest_Files.csv")
  
  to_json_scalar <- function(x) {
    if (length(x) == 0 || all(is.na(x))) return("null")
    if (is.numeric(x)) return(as.character(x))
    s <- gsub("\\\\", "\\\\\\\\", as.character(x))
    s <- gsub("\"", "\\\\\"", s)
    paste0("\"", s, "\"")
  }
  to_json_value <- function(x) {
    if (length(x) > 1) paste0("[", paste(vapply(x, to_json_scalar, character(1)), collapse = ", "), "]")
    else to_json_scalar(x)
  }
  
  keys <- names(manifest_list)
  json_lines <- c("{", paste0("  \"", keys, "\": ", vapply(manifest_list, to_json_value, character(1)), collapse = ",\n"), "\n}")
  writeLines(json_lines, con = json_path)
  
  inv <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  inv_df <- data.frame(
    file = basename(inv),
    rel_path = sub(paste0("^", gsub("([\\^\\$\\.|\\(\\)\\[\\]\\*\\+\\?\\\\])", "\\\\\\1", norm_path(out_dir)), "/?"), "", norm_path(inv)),
    abs_path = norm_path(inv),
    stringsAsFactors = FALSE
  )
  write.csv(inv_df, csv_path, row.names = FALSE)
  
  list(json_path = json_path, files_csv = csv_path)
}

# =========================
# Quarto report generation (canonical contract)
# =========================
render_quarto_report <- function(output_dir, header_text, params_list) {
  if (!has_quarto) {
    message("Quarto not available (package 'quarto' not installed). Skipping report.")
    return(list(qmd_path = NA_character_, html_path = NA_character_))
  }
  
  qmd_path <- file.path(output_dir, "Variance_ByDay_Cluster_Report.qmd")
  
  qmd <- c(
    "---",
    "title: \"Variance-by-Day and Cluster Variance-of-Variance\"",
    "format: html",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "## Summary",
    "This report computes per-gene inter-individual variance across biological replicates within each day, and summarizes the distribution of these per-gene variances across gene clusters. Boxplots show the *variance-of-variance*: for each cluster and day, the distribution of per-gene variability across genes in that cluster. Additional plots summarize mean/median variability trajectories, explicit var/sd across genes (var-of-var), and top-gene trajectories.",
    "",
    "## Script header (verbatim)",
    "```",
    header_text,
    "```",
    "",
    "## Metadata",
    "```",
    paste0("timestamp: ", params_list$timestamp),
    paste0("output_dir: ", params_list$output_dir),
    paste0("expression_file: ", params_list$expression_file),
    paste0("cluster_file: ", ifelse(is.na(params_list$cluster_file), "NA", params_list$cluster_file)),
    paste0("pseudocount: ", params_list$pseudocount),
    paste0("plot_scale_option: ", params_list$plot_scale_option),
    paste0("do_residuals: ", params_list$do_residuals),
    paste0("topN_per_cluster: ", params_list$topN),
    "```",
    "",
    "## Generated outputs",
    "```{r}",
    "files <- list.files(getwd(), recursive = TRUE)",
    "print(files)",
    "```",
    "",
    "## Figures",
    "```{r}",
    "fig_dir <- file.path(getwd(), \"figures\")",
    "figs <- list.files(fig_dir, pattern = \"\\\\.png$\", full.names = TRUE)",
    "if (length(figs) == 0) {",
    "  cat('No figures found.\\n')",
    "} else {",
    "  for (f in figs) {",
    "    cat('\\n### ', basename(f), '\\n', sep='')",
    if (has_knitr) "    knitr::include_graphics(f)" else "    cat('(knitr not installed; cannot embed images)\\n')",
    "  }",
    "}",
    "```",
    "",
    "## Session info",
    "```{r}",
    "sessionInfo()",
    "```"
  )
  
  writeLines(qmd, con = qmd_path)
  
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  quarto::quarto_render(basename(qmd_path))
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  list(qmd_path = qmd_path, html_path = html_path)
}

# =========================
# Main function
# =========================
variance_by_day_cluster <- function() {
  sid <- capture_script_identity()
  script_name <- sid$script_name
  script_path <- sid$script_path
  script_full <- sid$script_full
  
  message("\nINPUT EXPECTED (IMPORTANT):")
  message(" - WIDE expression table:")
  message("     * rows   = biological replicates (samples)")
  message("     * cols   = genes (numeric expression)")
  message(" - Sample IDs MUST be the row identifiers, not gene column headers.")
  message("     Acceptable formats:")
  message("       (A) Sample IDs in the FIRST COLUMN (non-numeric), e.g., '14_3'")
  message("       (B) Sample IDs as ROWNAMES (row names), e.g., '14_3'")
  message(" - Gene names MUST be the COLUMN HEADERS (all other columns must be numeric).")
  message(" - Sample IDs must encode Day and Rep as '#_#' (e.g., '14_3').")
  message("")
  
  expr_file <- norm_path(pick_file("Select WIDE expression CSV/TSV (rows=samples w/ SampleID in rownames or first column; cols=genes)"))
  # NEW: use input filename stub in output naming
  input_stub <- tools::file_path_sans_ext(basename(expr_file))
  input_stub <- gsub("[^A-Za-z0-9_\\-]+", "_", input_stub)
  
  has_cluster <- tolower(trimws(readline("\nDo you have a gene→cluster table (columns: gene, cluster)? (y/n): ")))
  cluster_file <- NA_character_
  if (has_cluster %in% c("y", "yes")) {
    cluster_file <- norm_path(pick_file("Select gene→cluster table CSV/TSV (columns: gene, cluster)"))
  }
  
  message("\nPLOT SCALE OPTION:")
  message("  1 = raw variance")
  message("  2 = log2(variance + pseudocount)  [recommended for boxplots]")
  opt <- trimws(readline("Choose 1 or 2 (default 2): "))
  plot_scale_option <- suppressWarnings(as.integer(opt))
  if (!is.finite(plot_scale_option) || !(plot_scale_option %in% c(1, 2))) plot_scale_option <- 2
  
  pc_in <- trimws(readline("\nPseudocount for log2(variance + pseudocount) (default 1e-8): "))
  pseudocount <- suppressWarnings(as.numeric(pc_in))
  if (!is.finite(pseudocount)) pseudocount <- 1e-8
  
  do_resid <- tolower(trimws(readline("Compute mean-corrected residual log2-variance (recommended)? (y/n, default y): ")))
  do_residuals <- !(do_resid %in% c("n", "no"))
  
  topN_in <- trimws(readline("\nTop N genes per cluster to plot trajectories (default 10; set 0 to skip): "))
  topN <- suppressWarnings(as.integer(topN_in))
  if (!is.finite(topN) || topN < 0) topN <- 10
  
  safe_mkdir("outputs")
  run_name <- gsub("[^A-Za-z0-9_\\-]+", "_", trimws(readline("\nEnter a run name (short, no spaces preferred): ")))
  if (!nzchar(run_name)) run_name <- "VarianceByDayCluster"
  ts <- timestamp_string()
  
  # NEW: include input_stub in folder name
  out_dir <- norm_path(file.path("outputs", paste0(run_name, "_", input_stub, "_", ts)))
  safe_mkdir(file.path(out_dir, "tables"))
  safe_mkdir(file.path(out_dir, "figures"))
  message("\nOutput directory (absolute):\n", out_dir)
  
  # Load expression table
  expr_df <- read_table_auto(expr_file)
  # ---- Quick transpose sanity check: if MANY column names look like '#_#', user likely provided samples as columns ----
  cn <- colnames(expr_df)
  looks_like_sample_ids <- sum(grepl("^\\D*\\d+_\\d+\\D*$", cn))
  if (looks_like_sample_ids >= max(3, round(0.3 * length(cn)))) {
    stop(
      "It looks like your SAMPLE IDs are in the COLUMN HEADERS (file may be transposed).\n",
      "This script expects SAMPLE IDs as rownames or the FIRST COLUMN, and GENE names as column headers.\n",
      "Fix: transpose the table so rows=samples and cols=genes, then re-run."
    )
  }
  # Determine SampleID
  sample_ids <- NULL
  if (ncol(expr_df) >= 2) {
    first_col <- expr_df[[1]]
    if (!is.numeric(first_col) && any(grepl("_", first_col))) {
      sample_ids <- as.character(first_col)
      expr_df <- expr_df[, -1, drop = FALSE]
    }
  }
  if (is.null(sample_ids)) {
    if (!is.null(rownames(expr_df)) && all(nzchar(rownames(expr_df)))) sample_ids <- rownames(expr_df)
    else stop("Could not determine SampleID. Provide sample IDs as first column or rownames.")
  }
  
  # Expression matrix
  expr_mat <- as.matrix(expr_df)
  suppressWarnings(storage.mode(expr_mat) <- "numeric")
  colnames(expr_mat) <- colnames(expr_df)
  
  meta <- parse_day_rep(sample_ids)
  
  # Per-gene per-day stats
  gene_day <- compute_gene_day_stats(expr_mat, meta$Day, pseudocount = pseudocount)
  
  # Add log2Var column always
  gene_day$log2Var <- log2(gene_day$Variance + pseudocount)
  
  # Residuals (log2)
  if (do_residuals) {
    gene_day <- add_residual_log2var(gene_day, pseudocount = pseudocount)
    if (!("log2ResVar" %in% colnames(gene_day))) gene_day$log2ResVar <- gene_day$ResLog2Var
  } else {
    gene_day$log2Mean <- log2(gene_day$Mean + pseudocount)
    gene_day$ResLog2Var <- NA_real_
    gene_day$log2ResVar <- NA_real_
  }
  
  # Save core tables (long + wide)
  gene_day_long_path <- file.path(out_dir, "tables", "GeneByDay_Variance_Long.csv")
  write.csv(gene_day, gene_day_long_path, row.names = FALSE)
  
  days <- sort(unique(gene_day$Day))
  genes <- unique(gene_day$Gene)
  
  make_wide <- function(metric) {
    wide <- matrix(NA_real_, nrow = length(genes), ncol = length(days),
                   dimnames = list(genes, paste0("Day", days)))
    for (d in days) {
      idx <- which(gene_day$Day == d)
      g <- gene_day$Gene[idx]
      v <- gene_day[[metric]][idx]
      wide[g, paste0("Day", d)] <- v
    }
    data.frame(Gene = rownames(wide), wide, row.names = NULL, check.names = FALSE)
  }
  
  write.csv(make_wide("Variance"), file.path(out_dir, "tables", "GeneByDay_Variance_Wide.csv"), row.names = FALSE)
  write.csv(make_wide("CV"),       file.path(out_dir, "tables", "GeneByDay_CV_Wide.csv"), row.names = FALSE)
  write.csv(make_wide("log2Var"),  file.path(out_dir, "tables", "GeneByDay_log2Var_Wide.csv"), row.names = FALSE)
  if (do_residuals) {
    write.csv(make_wide("ResLog2Var"), file.path(out_dir, "tables", "GeneByDay_ResLog2Var_Wide.csv"), row.names = FALSE)
  }
  
  # Cluster work (optional)
  gene_day_with_cluster <- gene_day
  fig_dir <- file.path(out_dir, "figures")
  
  if (!is.na(cluster_file)) {
    cl_df <- read_table_auto(cluster_file)
    message("\nCluster table columns detected:\n  ", paste(colnames(cl_df), collapse = ", "))
    message("\nUsing cluster assignment table:\n  ", cluster_file)
    message("\nIMPORTANT:")
    message("The following prompts refer to the GENE→CLUSTER ASSIGNMENT TABLE,")
    message("NOT the expression matrix you selected earlier.")
    message("")
    message("This table should contain:")
    message("  - one column with GENE NAMES")
    message("  - one column with CLUSTER IDs (numeric)")
    message("")
    message("You are now being asked to specify column names FROM THE CLUSTER TABLE.")
    message("")
    
    gene_col <- trimws(readline("Enter gene column name (default 'gene'): "))
    if (!nzchar(gene_col)) gene_col <- "gene"
    cluster_col <- trimws(readline("Enter cluster column name (default 'cluster'): "))
    if (!nzchar(cluster_col)) cluster_col <- "cluster"
    
    if (!(gene_col %in% colnames(cl_df))) stop("Gene column not found: ", gene_col)
    if (!(cluster_col %in% colnames(cl_df))) stop("Cluster column not found: ", cluster_col)
    
    map <- data.frame(
      Gene = as.character(cl_df[[gene_col]]),
      Cluster = suppressWarnings(as.integer(cl_df[[cluster_col]])),
      stringsAsFactors = FALSE
    )
    if (any(is.na(map$Cluster))) stop("Cluster column must be numeric/integer-like. Found non-numeric values.")
    
    gene_day_with_cluster <- merge(gene_day, map, by = "Gene", all.x = FALSE, all.y = FALSE)
    if (nrow(gene_day_with_cluster) == 0) stop("No overlap between expression genes and cluster genes.")
    
    merged_long_path <- file.path(out_dir, "tables", "GeneByDay_Variance_Long_WITH_Cluster.csv")
    write.csv(gene_day_with_cluster, merged_long_path, row.names = FALSE)
    
    metric_main <- if (plot_scale_option == 1) "Variance" else "log2Var"
    metric_label <- if (plot_scale_option == 1) "Variance (raw)" else "log2(Variance + pseudocount)"
    
    clusters <- sort(unique(gene_day_with_cluster$Cluster))
    topgenes_rows <- list()
    kk <- 1
    
    for (c in clusters) {
      dfc <- gene_day_with_cluster[gene_day_with_cluster$Cluster == c, , drop = FALSE]
      dfc <- dfc[is.finite(dfc[[metric_main]]), , drop = FALSE]
      if (nrow(dfc) < 5) next
      
      out_box <- file.path(fig_dir, paste0("Boxplot_", metric_main, "_Cluster_", c, ".png"))
      save_boxplot_by_cluster(
        df = dfc, out_png = out_box, metric_col = metric_main,
        title = paste0("Cluster ", c, ": per-gene variability across birds by Day (", metric_label, ")")
      )
      
      if (do_residuals) {
        dfc2 <- dfc[is.finite(dfc$ResLog2Var), , drop = FALSE]
        if (nrow(dfc2) >= 5) {
          out_box_res <- file.path(fig_dir, paste0("Boxplot_ResLog2Var_Cluster_", c, ".png"))
          save_boxplot_by_cluster(
            df = dfc2, out_png = out_box_res, metric_col = "ResLog2Var",
            title = paste0("Cluster ", c, ": mean-corrected variability by Day (Res log2 Var ~ log2 Mean)")
          )
        }
      }
      
      summ <- summarize_cluster_day(dfc, metric_col = metric_main)
      
      out_traj <- file.path(fig_dir, paste0("Trajectory_MeanMedian_", metric_main, "_Cluster_", c, ".png"))
      save_trajplot_cluster_summary(summ[summ$Cluster == c, , drop = FALSE], out_traj, metric_label)
      
      summ_path <- file.path(out_dir, "tables", paste0("ClusterByDay_Summary_", metric_main, "_Cluster_", c, ".csv"))
      write.csv(summ[summ$Cluster == c, , drop = FALSE], summ_path, row.names = FALSE)
      
      if (topN > 0) {
        genes_c <- unique(dfc$Gene)
        gene_score <- rep(NA_real_, length(genes_c))
        for (i in seq_along(genes_c)) {
          g <- genes_c[i]
          vv <- dfc[dfc$Gene == g, metric_main]
          vv <- vv[is.finite(vv)]
          gene_score[i] <- if (length(vv)) mean(vv) else NA_real_
        }
        ord <- order(gene_score, decreasing = TRUE, na.last = NA)
        top_genes <- genes_c[ord][seq_len(min(topN, length(ord)))]
        
        for (i in seq_along(top_genes)) {
          topgenes_rows[[kk]] <- data.frame(
            Cluster = c,
            Rank = i,
            Gene = top_genes[i],
            mean_metric_across_days = gene_score[match(top_genes[i], genes_c)],
            metric = metric_main,
            stringsAsFactors = FALSE
          )
          kk <- kk + 1
        }
        
        out_top <- file.path(fig_dir, paste0("Top", length(top_genes), "_Trajectories_", metric_main, "_Cluster_", c, ".png"))
        save_topN_gene_trajectories(dfc, out_top, metric_main, top_genes)
      }
    }
    
    if (length(topgenes_rows)) {
      topgenes_df <- do.call(rbind, topgenes_rows)
      topgenes_path <- file.path(out_dir, "tables", paste0("TopGenes_ByCluster_", if (plot_scale_option == 1) "VarianceRaw" else "log2Var", ".csv"))
      write.csv(topgenes_df, topgenes_path, row.names = FALSE)
    }
    
    # Global cluster summaries (ALL) including NEW var/sd across genes
    summ_all <- summarize_cluster_day(gene_day_with_cluster, metric_col = metric_main)
    summ_all_path <- file.path(out_dir, "tables", paste0("ClusterByDay_Summary_", metric_main, "_ALL.csv"))
    write.csv(summ_all, summ_all_path, row.names = FALSE)
    
    # NEW: var-of-var trajectory plots across all clusters (faceted)
    save_varofvar_trajectories_all(
      summ_all,
      out_png = file.path(fig_dir, "VarOfVar_Trajectory_ALL.png"),
      ycol = "var_metric",
      title = paste0("Var across genes of ", metric_main, " (variance-of-variance) vs Day (by cluster)")
    )
    save_varofvar_trajectories_all(
      summ_all,
      out_png = file.path(fig_dir, "SDofVar_Trajectory_ALL.png"),
      ycol = "sd_metric",
      title = paste0("SD across genes of ", metric_main, " vs Day (by cluster)")
    )
    
    if (do_residuals) {
      summ_res_all <- summarize_cluster_day(gene_day_with_cluster, metric_col = "ResLog2Var")
      summ_res_path <- file.path(out_dir, "tables", "ClusterByDay_Summary_ResLog2Var_ALL.csv")
      write.csv(summ_res_all, summ_res_path, row.names = FALSE)
      
      save_varofvar_trajectories_all(
        summ_res_all,
        out_png = file.path(fig_dir, "VarOfVar_Trajectory_ResLog2Var_ALL.png"),
        ycol = "var_metric",
        title = "Var across genes of ResLog2Var vs Day (mean-corrected variance-of-variance; by cluster)"
      )
      save_varofvar_trajectories_all(
        summ_res_all,
        out_png = file.path(fig_dir, "SDofVar_Trajectory_ResLog2Var_ALL.png"),
        ycol = "sd_metric",
        title = "SD across genes of ResLog2Var vs Day (mean-corrected; by cluster)"
      )
    }
    
  } else {
    message("\nNo cluster table provided.")
    message("You still have per-gene by-day variance tables.")
    message("If you want gene-set boxplots, provide a gene→cluster table where each gene set is a cluster ID.")
  }
  
  # README
  readme_path <- file.path(out_dir, "README_outputs.txt")
  readme <- c(
    "Variance-by-Day + Cluster Variance-of-Variance Outputs",
    "-----------------------------------------------------",
    paste0("Output directory: ", out_dir),
    "",
    "Core tables:",
    " - tables/GeneByDay_Variance_Long.csv",
    " - tables/GeneByDay_Variance_Wide.csv",
    " - tables/GeneByDay_CV_Wide.csv",
    " - tables/GeneByDay_log2Var_Wide.csv",
    if (do_residuals) " - tables/GeneByDay_ResLog2Var_Wide.csv" else NULL,
    "",
    "If a cluster table was provided:",
    " - tables/GeneByDay_Variance_Long_WITH_Cluster.csv",
    " - tables/ClusterByDay_Summary_<metric>_ALL.csv  (includes iqr_metric, var_metric, sd_metric)",
    " - figures/Boxplot_<metric>_Cluster_<k>.png",
    " - figures/Trajectory_MeanMedian_<metric>_Cluster_<k>.png",
    " - figures/VarOfVar_Trajectory_ALL.png  (var_metric vs Day, faceted by cluster)",
    " - figures/SDofVar_Trajectory_ALL.png   (sd_metric vs Day, faceted by cluster)",
    if (do_residuals) " - figures/Boxplot_ResLog2Var_Cluster_<k>.png" else NULL,
    if (do_residuals) " - figures/VarOfVar_Trajectory_ResLog2Var_ALL.png" else NULL,
    if (do_residuals) " - figures/SDofVar_Trajectory_ResLog2Var_ALL.png" else NULL,
    " - figures/TopN_Trajectories_<metric>_Cluster_<k>.png (if enabled)",
    "",
    "Metric choices:",
    " - Variance: raw per-gene variance across birds (within day)",
    " - log2Var: log2(Variance + pseudocount) (recommended for boxplots)",
    if (do_residuals) " - ResLog2Var: residuals from log2Var ~ log2Mean (mean-corrected variability)" else NULL,
    "",
    "Interpretation note:",
    "Boxplots show the distribution of per-gene variability values across genes within a cluster at each day.",
    "This directly visualizes 'variance of variance' for a gene set.",
    "In the cluster/day summary tables:",
    " - iqr_metric is a robust spread across genes of the per-gene variability metric at that day.",
    " - var_metric is the literal variance across genes of that metric (explicit variance-of-variance).",
    " - sd_metric is the SD across genes of that metric."
  )
  writeLines(readme[!vapply(readme, is.null, logical(1))], con = readme_path)
  
  # Manifest + inventory
  manifest_list <- list(
    script_name = script_name,
    script_path = script_path,
    script_full = script_full,
    timestamp = ts,
    input_paths = c(expr_file, if (!is.na(cluster_file)) cluster_file else NULL),
    output_dir = out_dir,
    pseudocount = pseudocount,
    plot_scale_option = plot_scale_option,
    do_residuals = do_residuals,
    topN = topN
  )
  man_paths <- write_manifest(out_dir, manifest_list)
  
  # Quarto report
  header_text <- get_header_block(script_full, fallback_name = script_name)
  params_list <- list(
    timestamp = ts,
    output_dir = out_dir,
    expression_file = expr_file,
    cluster_file = cluster_file,
    pseudocount = pseudocount,
    plot_scale_option = plot_scale_option,
    do_residuals = do_residuals,
    topN = topN
  )
  qr <- render_quarto_report(out_dir, header_text, params_list)
  
  # Return
  list(
    out_dir = out_dir,
    tables = list.files(file.path(out_dir, "tables"), full.names = TRUE),
    figures = list.files(file.path(out_dir, "figures"), full.names = TRUE),
    readme = readme_path,
    manifest_json = man_paths$json_path,
    manifest_files_csv = man_paths$files_csv,
    quarto_qmd = qr$qmd_path,
    quarto_html = qr$html_path,
    preview_table = if (!is.na(cluster_file)) gene_day_with_cluster else gene_day
  )
}

# =========================
# AUTO-RUN WHEN SOURCED (RStudio / interactive use)
# =========================
if (interactive()) {
  message("Running variance_by_day_cluster() interactively...")
  res <- variance_by_day_cluster()
  message("Done. Results written to:\n", res$out_dir)
  if (!is.null(res$preview_table)) {
    try(View(res$preview_table), silent = TRUE)
  }
  if (!is.na(res$quarto_html) && file.exists(res$quarto_html)) {
    message("Quarto HTML report:\n", res$quarto_html)
  }
}
