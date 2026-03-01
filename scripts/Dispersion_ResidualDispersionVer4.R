#!/usr/bin/env Rscript
# ======================================================================
# Publication Plotter (Canonical, Interactive, Provenance-First)
#
# HOW TO USE (RStudio; biologist-friendly)
#   1) Open this script in RStudio.
#   2) Set your working directory to your project root (so ./outputs/ is created there).
#      (Session > Set Working Directory > Choose Directory...)
#   3) Click Source.
#
# WHAT YOU WILL BE ASKED
#   A) Input format:
#      [1] WIDE  = one row per FEATURE (gene/metabolite), many sample columns named DAY_REP (e.g., 14_1)
#      [2] TIDY  = one row per observation, columns like SampleID, Day (optional), Value, Feature (optional), Group (optional)
#
#   B) Roles (TIDY mode; no second-guessing):
#      - FEATURE column (gene/metabolite ID): optional (REQUIRED for PCA/cutoffs)
#      - SAMPLE ID column: required (may be DAY_REP like 14_3)
#      - GROUP column: optional stratifier (Source/Tissue/Pathway/Line/etc.)
#      - VALUE column: required numeric (fpkm/tpm/count/log2p1/etc.)
#      - DAY column: optional; choose None to parse Day from SAMPLE ID
#
# OUTPUTS (always under ./outputs/)
#   outputs/<run_name>_<timestamp>/
#     manifest.json (or manifest.txt)
#     inventory.txt
#     tables/
#       long_tidy_all.csv
#       long_tidy_cutoff.csv (optional)
#       system_per_sample_all.csv
#       system_per_day_all.csv
#       wide_summary_all.csv
#       wide_summary_cutoff.csv (optional; feature-level)
#     plots/
#       png/      (600 dpi)
#       pdf/
#       plotly/   (self-contained HTML)
#     plots_report.qmd
#     plots_report.html   (if Quarto installed) OR plots_index.html fallback
#
# NOTES
#   - No silent overrides: if you select a column for a role, the script uses it.
#   - “GROUP” is a biological stratifier; SAMPLE IDs are NOT “grouping” (they are identifiers).
#   - Quarto contract is enforced: setwd(output_dir) before render, then restore.
# ======================================================================

# ---------------------------- Utilities ----------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (nchar(x) == 0) "unnamed" else x
}

norm_abs <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

as_numeric_strict <- function(x, name = "value") {
  y <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(y))) stop2("Column '", name, "' must be numeric (found non-numeric or NA after coercion).")
  y
}

n_unique <- function(x) length(unique(x[!is.na(x)]))

choose_yes_no <- function(prompt, default = "y") {
  ans <- readline(paste0(prompt, " (y/n) [", default, "]: "))
  if (ans == "") ans <- default
  tolower(ans) == "y"
}

# ---------------------------- Package bootstrap ----------------------------

needed <- c("ggplot2", "plotly", "htmlwidgets")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

has_rstudioapi <- requireNamespace("rstudioapi", quietly = TRUE)
has_jsonlite  <- requireNamespace("jsonlite", quietly = TRUE)
has_ragg      <- requireNamespace("ragg", quietly = TRUE)
has_quarto_r  <- requireNamespace("quarto", quietly = TRUE) # optional; we can also use system quarto

suppressPackageStartupMessages({
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
})

# ---------------------------- Script identity (best-effort, early) ----------------------------

get_script_identity <- function() {
  script_full <- NA_character_
  script_path <- NA_character_
  script_name <- NA_character_
  
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) script_full <- sub("^--file=", "", file_arg)
  
  if (is.na(script_full) && has_rstudioapi) {
    script_full <- tryCatch({
      p <- rstudioapi::getSourceEditorContext()$path
      if (is.character(p) && nzchar(p)) p else NA_character_
    }, error = function(e) NA_character_)
  }
  
  if (!is.na(script_full) && nzchar(script_full)) {
    script_full <- norm_abs(script_full)
    script_path <- dirname(script_full)
    script_name <- basename(script_full)
  }
  
  list(
    script_name = script_name,
    script_path = script_path,
    script_full = script_full
  )
}

script_id <- get_script_identity()

# ---------------------------- Mac-friendly pickers + fallbacks ----------------------------

pick_file_csv <- function() {
  if (interactive() && has_rstudioapi) {
    p <- tryCatch(
      rstudioapi::selectFile(caption = "Select input CSV", filter = "CSV (*.csv)"),
      error = function(e) NULL
    )
    if (!is.null(p) && nzchar(p)) return(p)
  }
  if (interactive()) return(file.choose())
  stop2("Non-interactive session. Run interactively or hardcode input_csv.")
}

# ---------------------------- Column selection UI ----------------------------

show_cols <- function(df) {
  cols <- colnames(df)
  for (i in seq_along(cols)) cat(sprintf("  [%d] %s\n", i, cols[i]))
  invisible(cols)
}

choose_column <- function(df, prompt, allow_none = FALSE) {
  cat("\n", prompt, "\n", sep = "")
  if (allow_none) cat("  [0] None\n")
  show_cols(df)
  idx <- suppressWarnings(as.integer(readline("Enter column number: ")))
  cols <- colnames(df)
  if (allow_none && identical(idx, 0L)) return(NULL)
  if (is.na(idx) || idx < 1 || idx > length(cols)) stop2("Invalid column selection.")
  cols[idx]
}

# ---------------------------- Deterministic sample name parsing DAY_REP ----------------------------

parse_sample_names_day_rep <- function(sample_names) {
  has_us <- grepl("_", sample_names, fixed = TRUE)
  if (!all(has_us)) {
    bad <- sample_names[!has_us]
    stop2("Sample naming rule violated: sample IDs missing underscore '_': ",
          paste(head(bad, 15), collapse = ", "),
          if (length(bad) > 15) " ..." else "")
  }
  
  day <- sub("_.*$", "", sample_names)
  rep <- sub("^[^_]+_", "", sample_names)
  
  if (any(day == "" | rep == "")) {
    bad <- sample_names[day == "" | rep == ""]
    stop2("Sample naming rule violated: empty Day or Rep detected in: ",
          paste(head(bad, 15), collapse = ", "),
          if (length(bad) > 15) " ..." else "")
  }
  
  day_num <- suppressWarnings(as.numeric(day))
  day_is_num <- all(!is.na(day_num))
  
  data.frame(
    sample  = sample_names,
    Day     = day,
    Day_num = if (day_is_num) day_num else NA_real_,
    Rep     = rep,
    stringsAsFactors = FALSE
  )
}

# ---------------------------- WIDE -> LONG (base R) ----------------------------

wide_to_long <- function(df, feature_col, group_col, sample_cols) {
  feature <- as.character(df[[feature_col]])
  group <- if (!is.null(group_col)) as.character(df[[group_col]]) else rep(NA_character_, nrow(df))
  
  X <- as.matrix(df[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(X) <- "numeric")
  
  n_feat <- nrow(X)
  n_samp <- ncol(X)
  
  data.frame(
    gene   = rep(feature, times = n_samp),
    group  = rep(group, times = n_samp),
    sample = rep(colnames(X), each = n_feat),
    value  = as.vector(X),
    stringsAsFactors = FALSE
  )
}

# ---------------------------- Summaries (system-level + feature-level) ----------------------------

# system_per_sample: within each sample, aggregate across genes (mean/median)
compute_system_per_sample <- function(long_df) {
  # long_df must have: sample, Day, Day_num, value; may have group (stratifier)
  if (!all(c("sample","Day","value") %in% names(long_df))) stop2("Internal: long_df missing required columns.")
  
  # If group exists and has non-NA, compute per (group, sample); else per sample
  has_group <- ("group" %in% names(long_df)) && any(!is.na(long_df$group))
  key <- if (has_group) paste(long_df$group, long_df$sample, sep = "||") else long_df$sample
  
  idxs <- split(seq_len(nrow(long_df)), key)
  rows <- lapply(names(idxs), function(k) {
    ii <- idxs[[k]]
    v <- long_df$value[ii]
    out <- data.frame(
      sample = long_df$sample[ii][1],
      Day    = long_df$Day[ii][1],
      Day_num = if ("Day_num" %in% names(long_df)) long_df$Day_num[ii][1] else NA_real_,
      n_features = n_unique(long_df$gene[ii]),
      n_obs = sum(is.finite(v)),
      sample_mean   = mean(v, na.rm = TRUE),
      sample_median = median(v, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (has_group) out$group <- long_df$group[ii][1]
    out
  })
  out <- do.call(rbind, rows)
  
  # order by numeric day if available
  if ("Day_num" %in% names(out) && all(is.finite(out$Day_num))) {
    out <- out[order(out$Day_num, out$sample), , drop = FALSE]
  }
  out
}

# bootstrap CI for a vector
boot_ci <- function(x, stat = c("mean","median"), R = 2000, conf = 0.95) {
  stat <- match.arg(stat)
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  fun <- if (stat == "mean") mean else median
  n <- length(x)
  sims <- replicate(R, fun(sample(x, size = n, replace = TRUE)))
  alpha <- (1 - conf) / 2
  as.numeric(quantile(sims, probs = c(alpha, 1 - alpha), na.rm = TRUE))
}

# system_per_day: across samples within day, aggregate sample_mean or sample_median with CI
compute_system_per_day <- function(sys_sample_df, which_stat = c("sample_mean","sample_median"),
                                   ci_R = 2000, ci_conf = 0.95) {
  which_stat <- match.arg(which_stat)
  if (!all(c("Day","sample") %in% names(sys_sample_df))) stop2("Internal: sys_sample_df missing Day/sample.")
  if (!which_stat %in% names(sys_sample_df)) stop2("Internal: requested stat missing in sys_sample_df.")
  
  has_group <- ("group" %in% names(sys_sample_df)) && any(!is.na(sys_sample_df$group))
  key <- if (has_group) paste(sys_sample_df$group, sys_sample_df$Day, sep = "||") else sys_sample_df$Day
  idxs <- split(seq_len(nrow(sys_sample_df)), key)
  
  rows <- lapply(names(idxs), function(k) {
    ii <- idxs[[k]]
    x <- sys_sample_df[[which_stat]][ii]
    ci <- boot_ci(x, stat = if (which_stat == "sample_mean") "mean" else "median", R = ci_R, conf = ci_conf)
    out <- data.frame(
      Day = sys_sample_df$Day[ii][1],
      Day_num = if ("Day_num" %in% names(sys_sample_df)) sys_sample_df$Day_num[ii][1] else NA_real_,
      n_samples = n_unique(sys_sample_df$sample[ii]),
      center = if (which_stat == "sample_mean") mean(x, na.rm = TRUE) else median(x, na.rm = TRUE),
      ci_low = ci[1],
      ci_high = ci[2],
      stringsAsFactors = FALSE
    )
    if (has_group) out$group <- sys_sample_df$group[ii][1]
    out
  })
  out <- do.call(rbind, rows)
  
  if ("Day_num" %in% names(out) && all(is.finite(out$Day_num))) {
    ord_cols <- c("Day_num")
    if ("group" %in% names(out)) out <- out[order(out$group, out$Day_num), , drop = FALSE] else out <- out[order(out$Day_num), , drop = FALSE]
  }
  out
}

# feature-level summary by day (+group)
summarise_feature_values_by_day <- function(long_df) {
  # long_df must have gene, Day, value; may have group
  if (!all(c("gene","Day","value") %in% names(long_df))) stop2("Internal: long_df missing required columns.")
  has_group <- ("group" %in% names(long_df)) && any(!is.na(long_df$group))
  
  # group keys
  key <- if (has_group) paste(long_df$group, long_df$Day, sep = "||") else long_df$Day
  idxs <- split(seq_len(nrow(long_df)), key)
  
  q025 <- function(x) as.numeric(quantile(x, 0.025, na.rm = TRUE))
  q975 <- function(x) as.numeric(quantile(x, 0.975, na.rm = TRUE))
  
  rows <- lapply(names(idxs), function(k) {
    ii <- idxs[[k]]
    v <- long_df$value[ii]
    out <- data.frame(
      Day = long_df$Day[ii][1],
      Day_num = if ("Day_num" %in% names(long_df)) long_df$Day_num[ii][1] else NA_real_,
      n_features = n_unique(long_df$gene[ii]),
      n_obs = sum(is.finite(v)),
      mean = mean(v, na.rm = TRUE),
      median = median(v, na.rm = TRUE),
      sd = sd(v, na.rm = TRUE),
      q025 = q025(v),
      q975 = q975(v),
      stringsAsFactors = FALSE
    )
    if (has_group) out$group <- long_df$group[ii][1]
    out
  })
  
  out <- do.call(rbind, rows)
  if ("Day_num" %in% names(out) && all(is.finite(out$Day_num))) {
    if ("group" %in% names(out)) out <- out[order(out$group, out$Day_num), , drop = FALSE] else out <- out[order(out$Day_num), , drop = FALSE]
  }
  out
}

# ---------------------------- Plot theme + constructors ----------------------------

theme_pub <- function(base_size = 11, base_family = "Helvetica", legend_pos = "right") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = grid::unit(2.2, "mm"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(size = ggplot2::rel(0.95)),
      plot.margin = ggplot2::margin(6, 10, 6, 6),
      legend.position = legend_pos,
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank()
    )
}

plot_timecourse_system <- function(sys_day_df, legend_pos = "right",
                                   title = "System time course",
                                   xlab = "Day", ylab = "Center (per sample)",
                                   gray_ci = TRUE) {
  has_group <- ("group" %in% names(sys_day_df)) && any(!is.na(sys_day_df$group))
  xcol <- if ("Day_num" %in% names(sys_day_df) && all(is.finite(sys_day_df$Day_num))) "Day_num" else "Day"
  
  if (has_group) {
    p <- ggplot(sys_day_df, aes(x = .data[[xcol]], y = center, color = group, group = group)) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), linewidth = 0.6, width = 0) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.0, alpha = 0.9) +
      labs(title = title, x = xlab, y = ylab) +
      theme_pub(legend_pos = legend_pos)
  } else {
    p <- ggplot(sys_day_df, aes(x = .data[[xcol]], y = center, group = 1)) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                    linewidth = 0.6, width = 0,
                    color = if (gray_ci) "grey30" else NULL) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.0, alpha = 0.9) +
      labs(title = title, x = xlab, y = ylab) +
      theme_pub(legend_pos = legend_pos)
  }
  p
}

plot_hist_values <- function(long_df, legend_pos = "right", title = "Histogram", xlab = "Value") {
  has_group <- ("group" %in% names(long_df)) && any(!is.na(long_df$group))
  if (has_group) {
    p <- ggplot(long_df, aes(x = value, fill = group)) +
      geom_histogram(bins = 30, alpha = 0.55, color = "white", linewidth = 0.2, position = "identity")
  } else {
    p <- ggplot(long_df, aes(x = value)) +
      geom_histogram(bins = 30, alpha = 0.9, color = "white", linewidth = 0.2)
  }
  p + labs(title = title, x = xlab, y = "Count") + theme_pub(legend_pos = legend_pos)
}

plot_timecourse_feature_summary <- function(wide_summary_df, legend_pos = "right",
                                            title = "Feature-value summary time course",
                                            xlab = "Day", ylab = "Mean across all values",
                                            show_ribbon = TRUE, gray_ribbon = TRUE) {
  has_group <- ("group" %in% names(wide_summary_df)) && any(!is.na(wide_summary_df$group))
  xcol <- if ("Day_num" %in% names(wide_summary_df) && all(is.finite(wide_summary_df$Day_num))) "Day_num" else "Day"
  
  if (has_group) {
    p <- ggplot(wide_summary_df, aes(x = .data[[xcol]], y = mean, color = group, fill = group, group = group))
    if (show_ribbon) p <- p + geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.18, color = NA)
    p <- p + geom_line(linewidth = 0.9) + geom_point(size = 2.0, alpha = 0.9)
  } else {
    p <- ggplot(wide_summary_df, aes(x = .data[[xcol]], y = mean, group = 1))
    if (show_ribbon) {
      p <- p + geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.18,
                           fill = if (gray_ribbon) "grey70" else NA, color = NA)
    }
    p <- p + geom_line(linewidth = 0.9) + geom_point(size = 2.0, alpha = 0.9)
  }
  
  p + labs(title = title, x = xlab, y = ylab) + theme_pub(legend_pos = legend_pos)
}

plot_pca_samples <- function(sample_by_feature_mat, sample_meta_df = NULL,
                             color_col = NULL, title = "PCA (samples)",
                             legend_pos = "right") {
  # sample_by_feature_mat: rows = samples, cols = features
  X <- as.matrix(sample_by_feature_mat)
  ok <- apply(X, 2, function(z) all(is.finite(z)))
  if (!all(ok)) X <- X[, ok, drop = FALSE]
  if (ncol(X) < 2) stop2("PCA requires >= 2 numeric feature columns after filtering non-finite.")
  
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  colnames(scores) <- c("PC1", "PC2")
  scores$sample <- rownames(X)
  
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  xlab <- paste0("PC1 (", round(100 * var_expl[1], 1), "%)")
  ylab <- paste0("PC2 (", round(100 * var_expl[2], 1), "%)")
  
  if (!is.null(sample_meta_df) && "sample" %in% names(sample_meta_df)) {
    scores <- merge(scores, sample_meta_df, by.x = "sample", by.y = "sample", all.x = TRUE, sort = FALSE)
  }
  
  if (!is.null(color_col) && color_col %in% names(scores)) {
    p <- ggplot(scores, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
      geom_point(size = 2.0, alpha = 0.9)
  } else {
    p <- ggplot(scores, aes(x = PC1, y = PC2)) +
      geom_point(size = 2.0, alpha = 0.9)
  }
  
  p + labs(title = title, x = xlab, y = ylab) + theme_pub(legend_pos = legend_pos)
}

plot_xy_scatter <- function(df, x_col, y_col, group_col = NULL,
                            title = "XY scatter", xlab = NULL, ylab = NULL,
                            legend_pos = "right", add_identity = FALSE, add_lm = FALSE) {
  if (!x_col %in% names(df) || !y_col %in% names(df)) stop2("XY scatter: selected columns not found.")
  if (!is.numeric(df[[x_col]]) || !is.numeric(df[[y_col]])) stop2("XY scatter: x and y must be numeric.")
  
  if (!is.null(group_col) && group_col %in% names(df)) {
    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[group_col]])) +
      geom_point(alpha = 0.75, size = 2.0)
  } else {
    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
      geom_point(alpha = 0.75, size = 2.0, color = "black")
  }
  
  if (add_identity) p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7)
  if (add_lm) p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7)
  
  p + labs(title = title, x = xlab %||% x_col, y = ylab %||% y_col) + theme_pub(legend_pos = legend_pos)
}

# ---------------------------- Export functions (pub + plotly) ----------------------------

save_pub <- function(p, png_dir, pdf_dir, stem, width = 6.5, height = 4.0) {
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf_path <- file.path(pdf_dir, paste0(stem, ".pdf"))
  png_path <- file.path(png_dir, paste0(stem, ".png"))
  
  if (capabilities("cairo")) {
    ggsave(pdf_path, p, width = width, height = height, units = "in", device = cairo_pdf)
  } else {
    ggsave(pdf_path, p, width = width, height = height, units = "in")
  }
  
  if (has_ragg) {
    ggsave(png_path, p, width = width, height = height, units = "in", dpi = 600, device = ragg::agg_png)
  } else {
    ggsave(png_path, p, width = width, height = height, units = "in", dpi = 600)
  }
  
  list(pdf = pdf_path, png = png_path)
}

save_plotly_html <- function(p, plotly_dir, stem, width_in = 6.5, height_in = 4.0) {
  dir.create(plotly_dir, recursive = TRUE, showWarnings = FALSE)
  html_path <- file.path(plotly_dir, paste0(stem, ".plotly.html"))
  
  w_px <- as.integer(width_in * 96)
  h_px <- as.integer(height_in * 96)
  
  gg <- plotly::ggplotly(p, width = w_px, height = h_px)
  gg <- plotly::config(gg, responsive = TRUE)
  
  htmlwidgets::saveWidget(gg, file = html_path, selfcontained = TRUE)
  html_path
}

# ---------------------------- Manifest + inventory ----------------------------

write_manifest <- function(out_dir, manifest) {
  if (!has_jsonlite) {
    p <- file.path(out_dir, "manifest.txt")
    lines <- unlist(Map(function(k, v) paste0(k, ": ", paste(v, collapse = "; ")),
                        names(manifest), manifest))
    writeLines(lines, p)
    return(p)
  }
  p <- file.path(out_dir, "manifest.json")
  jsonlite::write_json(manifest, path = p, pretty = TRUE, auto_unbox = TRUE)
  p
}

write_inventory <- function(out_dir) {
  inv <- file.path(out_dir, "inventory.txt")
  files <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
  writeLines(files, inv)
  inv
}

# ---------------------------- Quarto report (canonical contract) ----------------------------

write_quarto_report <- function(output_dir, title, script_header_lines, manifest_path, inventory_path, plots_dir, tables_dir) {
  qmd_path <- file.path(output_dir, "plots_report.qmd")
  
  pngs <- list.files(file.path(plots_dir, "png"), pattern = "\\.png$", full.names = FALSE)
  pdfs <- list.files(file.path(plots_dir, "pdf"), pattern = "\\.pdf$", full.names = FALSE)
  plotly_htmls <- list.files(file.path(plots_dir, "plotly"), pattern = "\\.plotly\\.html$", full.names = FALSE)
  tables <- list.files(tables_dir, pattern = "\\.csv$", full.names = FALSE)
  
  blocks <- c(
    "---",
    paste0('title: "', gsub('"', '\\"', title), '"'),
    "format: html",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "## Summary",
    "",
    "This report documents the outputs generated by the Publication Plotter.",
    "",
    "## Script header",
    "",
    "```",
    script_header_lines,
    "```",
    "",
    "## Metadata",
    "",
    paste0("- Manifest: `", basename(manifest_path), "`"),
    paste0("- Inventory: `", basename(inventory_path), "`"),
    ""
  )
  
  blocks <- c(blocks, "## Tables", "")
  if (length(tables) == 0) {
    blocks <- c(blocks, "_No CSV tables were written._", "")
  } else {
    for (f in tables) blocks <- c(blocks, paste0("- `tables/", f, "`"))
    blocks <- c(blocks, "")
  }
  
  blocks <- c(blocks, "## Interactive plots (Plotly)", "")
  if (length(plotly_htmls) == 0) {
    blocks <- c(blocks, "_No interactive plots were written._", "")
  } else {
    for (f in plotly_htmls) {
      blocks <- c(
        blocks,
        paste0("### ", f),
        "",
        paste0("```{=html}\n<iframe src=\"plots/plotly/", f,
               "\" width=\"100%\" height=\"650\" style=\"border:0;\"></iframe>\n```"),
        ""
      )
    }
  }
  
  blocks <- c(blocks, "## Figures (PNG)", "")
  if (length(pngs) == 0) {
    blocks <- c(blocks, "_No PNG plots were written._", "")
  } else {
    for (f in pngs) blocks <- c(blocks, paste0("![](plots/png/", f, ")"), "")
  }
  
  blocks <- c(blocks, "## Figures (PDF)", "")
  if (length(pdfs) == 0) {
    blocks <- c(blocks, "_No PDF plots were written._", "")
  } else {
    for (f in pdfs) blocks <- c(blocks, paste0("- `plots/pdf/", f, "`"))
    blocks <- c(blocks, "")
  }
  
  blocks <- c(
    blocks,
    "## Session info",
    "",
    "```{r}",
    "sessionInfo()",
    "```",
    ""
  )
  
  writeLines(blocks, qmd_path)
  
  # Render per canonical contract
  qbin <- Sys.which("quarto")
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  
  if (nzchar(qbin)) {
    old_wd <- getwd()
    setwd(output_dir)
    on.exit(setwd(old_wd), add = TRUE)
    system2(qbin, args = c("render", basename(qmd_path), "--to", "html"), stdout = TRUE, stderr = TRUE)
    if (file.exists(html_path)) return(list(qmd = qmd_path, html = html_path, mode = "quarto"))
  }
  
  # fallback index
  idx_path <- file.path(output_dir, "plots_index.html")
  body <- c(
    "<!doctype html>",
    "<html><head><meta charset='utf-8'>",
    paste0("<title>", title, "</title>"),
    "<style>body{font-family: Arial, sans-serif; margin: 24px;} iframe{margin: 12px 0;}</style>",
    "</head><body>",
    paste0("<h1>", title, "</h1>"),
    paste0("<p>Manifest: <code>", basename(manifest_path), "</code></p>"),
    paste0("<p>Inventory: <code>", basename(inventory_path), "</code></p>")
  )
  
  if (length(plotly_htmls) > 0) {
    for (f in plotly_htmls) {
      body <- c(body,
                paste0("<h2>", f, "</h2>"),
                paste0("<iframe src='plots/plotly/", f, "' width='100%' height='650' style='border:0;'></iframe>"))
    }
  } else {
    body <- c(body, "<p><em>No interactive plots were generated.</em></p>")
  }
  
  body <- c(body, "</body></html>")
  writeLines(body, idx_path)
  list(qmd = qmd_path, html = idx_path, mode = "fallback")
}

# ---------------------------- Helpers: build sample x feature matrix for PCA ----------------------------

# From LONG: returns matrix rows=samples cols=features (genes/metabolites), with NA filled as 0 (or removed)
long_to_sample_feature_matrix <- function(long_df, fill_value = 0) {
  # requires gene, sample, value
  if (any(is.na(long_df$gene))) stop2("PCA requires FEATURE IDs (gene/metabolite). Your long table has missing FEATURE IDs.")
  if (n_unique(long_df$gene) < 2) stop2("PCA requires >=2 unique features. Found: ", n_unique(long_df$gene))
  if (n_unique(long_df$sample) < 2) stop2("PCA requires >=2 unique samples. Found: ", n_unique(long_df$sample))
  
  # aggregate duplicates (gene,sample) by mean
  key <- paste(long_df$sample, long_df$gene, sep = "||")
  idxs <- split(seq_len(nrow(long_df)), key)
  samp <- character(length(idxs))
  gene <- character(length(idxs))
  val  <- numeric(length(idxs))
  
  i <- 0
  for (k in names(idxs)) {
    i <- i + 1
    ii <- idxs[[k]]
    samp[i] <- long_df$sample[ii][1]
    gene[i] <- long_df$gene[ii][1]
    val[i]  <- mean(long_df$value[ii], na.rm = TRUE)
  }
  
  samp_levels <- unique(samp)
  gene_levels <- unique(gene)
  
  M <- matrix(fill_value, nrow = length(samp_levels), ncol = length(gene_levels),
              dimnames = list(samp_levels, gene_levels))
  # place values
  r <- match(samp, samp_levels)
  c <- match(gene, gene_levels)
  M[cbind(r, c)] <- val
  
  M
}

# ---------------------------- MAIN ----------------------------

cat("\n=== Publication Plotter (Canonical) ===\n")

# Step 1: input CSV picker
input_csv <- pick_file_csv()
input_csv <- norm_abs(input_csv)
cat("Input CSV:\n  ", input_csv, "\n", sep = "")

df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
cat("\nLoaded table: ", nrow(df), " rows x ", ncol(df), " cols\n", sep = "")

# Step 1b: format choice
cat("\nIs your input WIDE (feature x sample columns) or TIDY (one row per feature-sample)?\n")
cat("  [1] WIDE  = sample columns named DAY_REP like 4_1, 4_2, 6_1...\n")
cat("  [2] TIDY  = columns like SampleID, Day (optional), Value, Feature (optional)\n")
fmt <- suppressWarnings(as.integer(readline("Enter number [2]: ")))
if (is.na(fmt) || !(fmt %in% 1:2)) fmt <- 2
input_mode <- if (fmt == 1) "wide" else "tidy"
cat("Input mode: ", input_mode, "\n", sep = "")

# Placeholders
long_all <- NULL
sample_meta_df <- NULL
feature_col <- NULL
group_col <- NULL
sample_id_col <- NULL
day_col <- NULL
value_col <- NULL
sample_cols <- NULL

# ---------------------------- Build canonical long_all ----------------------------

if (input_mode == "wide") {
  
  cat("\nSelect roles (WIDE mode):\n")
  cat("  FEATURE = gene/metabolite ID column (required)\n")
  cat("  GROUP   = optional stratifier (optional)\n")
  cat("  SAMPLE columns = all remaining columns; must be named DAY_REP\n")
  
  feature_col <- choose_column(df, "Select FEATURE column (gene/metabolite ID) [required]:", allow_none = FALSE)
  group_col   <- choose_column(df, "Select GROUP column (optional stratifier) [optional]:", allow_none = TRUE)
  
  exclude <- c(feature_col, group_col)
  exclude <- exclude[!is.na(exclude)]
  sample_cols <- setdiff(colnames(df), exclude)
  if (length(sample_cols) < 2) stop2("WIDE mode: need >=2 sample columns after excluding FEATURE/GROUP.")
  
  cat("\nDetected sample columns: ", length(sample_cols), "\n", sep = "")
  cat("Validating sample naming rule DAY_REP...\n")
  sample_meta_df <- parse_sample_names_day_rep(sample_cols)
  cat("Sample naming validated.\n")
  
  long_all <- wide_to_long(df, feature_col = feature_col, group_col = group_col, sample_cols = sample_cols)
  
  # merge Day/Rep onto long
  long_all <- merge(long_all, sample_meta_df, by = "sample", all.x = TRUE, sort = FALSE)
  if (any(is.na(long_all$Day))) stop2("WIDE mode: internal error mapping sample->Day failed.")
  
} else {
  
  cat("\nSelect roles (TIDY mode):\n")
  cat("  FEATURE = gene/metabolite ID (optional; REQUIRED for PCA/cutoffs)\n")
  cat("  SAMPLE  = sample/replicate ID (required; may be DAY_REP like 14_3)\n")
  cat("  GROUP   = optional stratifier (Source/Tissue/Pathway/Line/etc.)\n")
  cat("  VALUE   = numeric expression/abundance (required)\n")
  cat("  DAY     = optional; choose None to parse from SAMPLE\n")
  
  feature_col   <- choose_column(df, "Select FEATURE column (gene/metabolite ID) [optional]:", allow_none = TRUE)
  sample_id_col <- choose_column(df, "Select SAMPLE ID column (replicate ID; may be DAY_REP) [required]:", allow_none = FALSE)
  group_col     <- choose_column(df, "Select GROUP column (optional stratifier) [optional]:", allow_none = TRUE)
  value_col     <- choose_column(df, "Select numeric VALUE column (fpkm/tpm/count/log2p1) [required]:", allow_none = FALSE)
  day_col       <- choose_column(df, "Select DAY column if present, or choose [0] None to parse Day from SAMPLE ID:", allow_none = TRUE)
  
  long_all <- data.frame(
    sample = as.character(df[[sample_id_col]]),
    value  = as_numeric_strict(df[[value_col]], name = value_col),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(feature_col)) long_all$gene <- as.character(df[[feature_col]]) else long_all$gene <- NA_character_
  if (!is.null(group_col))   long_all$group <- as.character(df[[group_col]]) else long_all$group <- NA_character_
  
  if (!is.null(day_col)) {
    long_all$Day <- as.character(df[[day_col]])
    long_all$Day_num <- suppressWarnings(as.numeric(long_all$Day))
    # Rep not required if Day exists; but if you want, you can still parse Rep later.
    sample_meta_df <- data.frame(
      sample = unique(long_all$sample),
      stringsAsFactors = FALSE
    )
    # optionally parse Day_num order if day is numeric
  } else {
    # parse Day/Rep from sample IDs
    sample_names <- unique(long_all$sample)
    sample_meta_df <- parse_sample_names_day_rep(sample_names)
    long_all <- merge(long_all, sample_meta_df, by = "sample", all.x = TRUE, sort = FALSE)
    if (any(is.na(long_all$Day))) stop2("TIDY mode: Day parsing failed for some samples.")
  }
  
  # Ensure Day_num exists (even if Day column was selected)
  if (!"Day_num" %in% names(long_all)) long_all$Day_num <- suppressWarnings(as.numeric(long_all$Day))
}

cat("\nCanonical long table ready.\n")
cat("  Rows: ", nrow(long_all), "\n", sep = "")
cat("  Unique samples: ", n_unique(long_all$sample), "\n", sep = "")
cat("  Unique days: ", n_unique(long_all$Day), "\n", sep = "")
cat("  FEATURE provided: ", if (all(is.na(long_all$gene))) "NO" else "YES", "\n", sep = "")
cat("  GROUP provided: ", if (("group" %in% names(long_all)) && any(!is.na(long_all$group))) "YES" else "NO", "\n", sep = "")

# ---------------------------- Output directory standardization ----------------------------

project_root <- getwd()
outputs_root <- file.path(project_root, "outputs")
dir.create(outputs_root, showWarnings = FALSE, recursive = TRUE)

run_name <- readline("\nEnter run name (short label; letters/numbers/underscore recommended) [plotrun]: ")
if (run_name == "") run_name <- "plotrun"
run_name <- safe_stem(run_name)

timestamp <- ts_stamp()
output_dir <- file.path(outputs_root, paste0(run_name, "_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tables_dir <- file.path(output_dir, "tables")
plots_dir  <- file.path(output_dir, "plots")
png_dir    <- file.path(plots_dir, "png")
pdf_dir    <- file.path(plots_dir, "pdf")
plotly_dir <- file.path(plots_dir, "plotly")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotly_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nOutput directory:\n  ", norm_abs(output_dir), "\n", sep = "")

# ---------------------------- Plot selection + options ----------------------------

cat("\nWhich plots to generate? Enter comma-separated list.\n")
cat("  [1] System time course (per-sample center by day + bootstrap CI)\n")
cat("  [2] Histogram of all values\n")
cat("  [3] Feature-value summary time course (mean + 95% interval across all values)\n")
cat("  [4] PCA (samples)  -- REQUIRES FEATURE column\n")
cat("  [5] XY scatter (choose any two numeric columns from tables we write)\n")
sel <- readline("Enter plot numbers [1,2]: ")
if (sel == "") sel <- "1,2"
want <- suppressWarnings(as.integer(strsplit(gsub("\\s+", "", sel), ",")[[1]]))
want <- want[want %in% 1:5]
if (length(want) == 0) stop2("No valid plot choices selected.")

cat("\nLegend position:\n  [1] right\n  [2] bottom\n  [3] none\n")
lp <- suppressWarnings(as.integer(readline("Enter number [1]: ")))
if (is.na(lp) || !(lp %in% 1:3)) lp <- 1
legend_pos <- c("right", "bottom", "none")[lp]

# system statistic choice for Plot [1]
sys_stat <- "sample_mean"
if (1 %in% want) {
  cat("\nSystem statistic (within each sample, across features/values):\n")
  cat("  [1] mean\n")
  cat("  [2] median\n")
  sidx <- suppressWarnings(as.integer(readline("Enter number [1]: ")))
  if (is.na(sidx) || !(sidx %in% 1:2)) sidx <- 1
  sys_stat <- if (sidx == 1) "sample_mean" else "sample_median"
}

# cutoff option (requires FEATURE and is meaningful mainly for WIDE or TIDY-with-feature)
do_cutoff <- FALSE
cutoff_N <- NA_integer_
if (!all(is.na(long_all$gene))) {
  do_cutoff <- choose_yes_no("\nApply cutoff filter by global mean (top N features)?", default = "n")
  if (do_cutoff) {
    cutoff_N <- suppressWarnings(as.integer(readline("Enter cutoff N (e.g., 100): ")))
    if (!is.finite(cutoff_N) || cutoff_N < 1) stop2("Invalid cutoff N.")
  }
} else {
  cat("\nNOTE: FEATURE column not provided; cutoff and PCA are disabled.\n")
}

# ---------------------------- Tables: write long + summaries ----------------------------

long_all_path <- file.path(tables_dir, "long_tidy_all.csv")
write.csv(long_all, long_all_path, row.names = FALSE)

# System summaries
sys_sample_all <- compute_system_per_sample(long_all)
sys_sample_all_path <- file.path(tables_dir, "system_per_sample_all.csv")
write.csv(sys_sample_all, sys_sample_all_path, row.names = FALSE)

sys_day_all <- compute_system_per_day(sys_sample_all, which_stat = sys_stat, ci_R = 2000, ci_conf = 0.95)
sys_day_all_path <- file.path(tables_dir, paste0("system_per_day_all_", sys_stat, ".csv"))
write.csv(sys_day_all, sys_day_all_path, row.names = FALSE)

# Feature-value summary by day (always possible, even if gene is NA; it summarizes values)
wide_all <- summarise_feature_values_by_day(long_all)
wide_all_path <- file.path(tables_dir, "wide_summary_all.csv")
write.csv(wide_all, wide_all_path, row.names = FALSE)

# Optional cutoff
long_cut <- NULL
sys_sample_cut <- NULL
sys_day_cut <- NULL
wide_cut <- NULL

if (do_cutoff) {
  # compute gene means across all samples
  # (Use long table: mean value per gene across all rows)
  if (all(is.na(long_all$gene))) stop2("Cutoff requested but FEATURE is missing.")
  key <- long_all$gene
  idxs <- split(seq_len(nrow(long_all)), key)
  gene_means <- sapply(idxs, function(ii) mean(long_all$value[ii], na.rm = TRUE))
  ord <- names(sort(gene_means, decreasing = TRUE))
  keep <- head(ord, cutoff_N)
  
  long_cut <- long_all[long_all$gene %in% keep, , drop = FALSE]
  long_cut_path <- file.path(tables_dir, "long_tidy_cutoff.csv")
  write.csv(long_cut, long_cut_path, row.names = FALSE)
  
  sys_sample_cut <- compute_system_per_sample(long_cut)
  sys_sample_cut_path <- file.path(tables_dir, "system_per_sample_cutoff.csv")
  write.csv(sys_sample_cut, sys_sample_cut_path, row.names = FALSE)
  
  sys_day_cut <- compute_system_per_day(sys_sample_cut, which_stat = sys_stat, ci_R = 2000, ci_conf = 0.95)
  sys_day_cut_path <- file.path(tables_dir, paste0("system_per_day_cutoff_", sys_stat, ".csv"))
  write.csv(sys_day_cut, sys_day_cut_path, row.names = FALSE)
  
  wide_cut <- summarise_feature_values_by_day(long_cut)
  wide_cut_path <- file.path(tables_dir, "wide_summary_cutoff.csv")
  write.csv(wide_cut, wide_cut_path, row.names = FALSE)
}

cat("\nTables written:\n")
cat("  ", long_all_path, "\n", sep = "")
cat("  ", sys_sample_all_path, "\n", sep = "")
cat("  ", sys_day_all_path, "\n", sep = "")
cat("  ", wide_all_path, "\n", sep = "")
if (do_cutoff) cat("  ", long_cut_path, "\n  ", sys_sample_cut_path, "\n  ", sys_day_cut_path, "\n  ", wide_cut_path, "\n", sep = "")

# ---------------------------- Plots ----------------------------

plotly_files <- character(0)

# [1] System time course
if (1 %in% want) {
  p1 <- plot_timecourse_system(
    sys_day_all, legend_pos = legend_pos,
    title = paste0("System time course (", gsub("sample_", "", sys_stat), " across features)"),
    xlab = "Day",
    ylab = paste0("Per-day center of per-sample ", gsub("sample_", "", sys_stat)),
    gray_ci = TRUE
  )
  stem <- paste0("system_timecourse_all_", sys_stat)
  save_pub(p1, png_dir, pdf_dir, stem, width = 6.8, height = 4.2)
  save_plotly_html(p1, plotly_dir, stem, width_in = 6.8, height_in = 4.2)
  plotly_files <- c(plotly_files, paste0(stem, ".plotly.html"))
  
  if (do_cutoff && !is.null(sys_day_cut)) {
    p1c <- plot_timecourse_system(
      sys_day_cut, legend_pos = legend_pos,
      title = paste0("System time course (cutoff top ", cutoff_N, " features)"),
      xlab = "Day",
      ylab = paste0("Per-day center of per-sample ", gsub("sample_", "", sys_stat)),
      gray_ci = TRUE
    )
    stemc <- paste0("system_timecourse_cutoff_top", cutoff_N, "_", sys_stat)
    save_pub(p1c, png_dir, pdf_dir, stemc, width = 6.8, height = 4.2)
    save_plotly_html(p1c, plotly_dir, stemc, width_in = 6.8, height_in = 4.2)
    plotly_files <- c(plotly_files, paste0(stemc, ".plotly.html"))
  }
}

# [2] Histogram
if (2 %in% want) {
  p2 <- plot_hist_values(long_all, legend_pos = legend_pos, title = "Histogram (all values)", xlab = "Value")
  stem <- "hist_all_values"
  save_pub(p2, png_dir, pdf_dir, stem, width = 6.8, height = 4.2)
  save_plotly_html(p2, plotly_dir, stem, width_in = 6.8, height_in = 4.2)
  plotly_files <- c(plotly_files, paste0(stem, ".plotly.html"))
}

# [3] Feature-value summary time course
if (3 %in% want) {
  p3 <- plot_timecourse_feature_summary(
    wide_all, legend_pos = legend_pos,
    title = "Feature-value summary time course (mean + 95% interval)",
    xlab = "Day", ylab = "Mean across all values",
    show_ribbon = TRUE, gray_ribbon = TRUE
  )
  stem <- "feature_value_summary_timecourse_all"
  save_pub(p3, png_dir, pdf_dir, stem, width = 6.8, height = 4.2)
  save_plotly_html(p3, plotly_dir, stem, width_in = 6.8, height_in = 4.2)
  plotly_files <- c(plotly_files, paste0(stem, ".plotly.html"))
  
  if (do_cutoff && !is.null(wide_cut)) {
    p3c <- plot_timecourse_feature_summary(
      wide_cut, legend_pos = legend_pos,
      title = paste0("Feature-value summary (cutoff top ", cutoff_N, " features)"),
      xlab = "Day", ylab = "Mean across all values",
      show_ribbon = TRUE, gray_ribbon = TRUE
    )
    stemc <- paste0("feature_value_summary_timecourse_cutoff_top", cutoff_N)
    save_pub(p3c, png_dir, pdf_dir, stemc, width = 6.8, height = 4.2)
    save_plotly_html(p3c, plotly_dir, stemc, width_in = 6.8, height_in = 4.2)
    plotly_files <- c(plotly_files, paste0(stemc, ".plotly.html"))
  }
}

# [4] PCA
if (4 %in% want) {
  if (all(is.na(long_all$gene))) {
    stop2("Plot [4] PCA requested, but FEATURE column is missing. Re-run and select a FEATURE column.")
  }
  # Build sample x feature matrix
  M <- long_to_sample_feature_matrix(long_all, fill_value = 0)
  
  # Attach sample meta if we have it
  sm <- NULL
  if (!is.null(sample_meta_df) && "sample" %in% names(sample_meta_df)) {
    sm <- sample_meta_df
  } else {
    sm <- data.frame(sample = rownames(M), stringsAsFactors = FALSE)
  }
  
  # Choose PCA color column from sample metadata if available
  color_col <- NULL
  cat("\nPCA coloring (optional):\n")
  cat("  [0] None\n")
  cols_sm <- setdiff(names(sm), "sample")
  if (length(cols_sm) > 0) {
    for (i in seq_along(cols_sm)) cat(sprintf("  [%d] %s\n", i, cols_sm[i]))
    ci <- suppressWarnings(as.integer(readline("Choose column number for PCA color [0]: ")))
    if (!is.na(ci) && ci >= 1 && ci <= length(cols_sm)) color_col <- cols_sm[ci]
  } else {
    cat("  (No sample metadata columns available)\n")
  }
  
  p4 <- plot_pca_samples(M, sample_meta_df = sm, color_col = color_col,
                         title = "PCA (samples)", legend_pos = legend_pos)
  stem <- "pca_samples"
  save_pub(p4, png_dir, pdf_dir, stem, width = 6.8, height = 5.2)
  save_plotly_html(p4, plotly_dir, stem, width_in = 6.8, height_in = 5.2)
  plotly_files <- c(plotly_files, paste0(stem, ".plotly.html"))
}

# [5] XY scatter
if (5 %in% want) {
  # Choose which table to use
  cat("\nXY scatter: choose a table to plot from:\n")
  cat("  [1] system_per_day_all\n")
  cat("  [2] system_per_sample_all\n")
  cat("  [3] wide_summary_all\n")
  ti <- suppressWarnings(as.integer(readline("Enter number [1]: ")))
  if (is.na(ti) || !(ti %in% 1:3)) ti <- 1
  
  dxy <- if (ti == 1) sys_day_all else if (ti == 2) sys_sample_all else wide_all
  
  # Show numeric columns
  num_cols <- names(dxy)[vapply(dxy, is.numeric, logical(1))]
  if (length(num_cols) < 2) stop2("Selected table has <2 numeric columns for XY scatter.")
  cat("\nNumeric columns:\n")
  for (i in seq_along(num_cols)) cat(sprintf("  [%d] %s\n", i, num_cols[i]))
  
  xi <- suppressWarnings(as.integer(readline("Select X column number: ")))
  yi <- suppressWarnings(as.integer(readline("Select Y column number: ")))
  if (is.na(xi) || xi < 1 || xi > length(num_cols)) stop2("Invalid X selection.")
  if (is.na(yi) || yi < 1 || yi > length(num_cols)) stop2("Invalid Y selection.")
  
  x_col <- num_cols[xi]
  y_col <- num_cols[yi]
  
  # Optional grouping/color column
  group_xy <- NULL
  if ("group" %in% names(dxy) && any(!is.na(dxy$group))) {
    if (choose_yes_no("Color by GROUP?", default = "y")) group_xy <- "group"
  }
  
  add_id <- choose_yes_no("Add y=x line?", default = "n")
  add_lm <- choose_yes_no("Add linear regression line?", default = "n")
  
  p5 <- plot_xy_scatter(dxy, x_col, y_col, group_col = group_xy,
                        title = paste0("XY: ", x_col, " vs ", y_col),
                        legend_pos = legend_pos,
                        add_identity = add_id, add_lm = add_lm)
  
  stem <- paste0("xy_", safe_stem(x_col), "_vs_", safe_stem(y_col))
  save_pub(p5, png_dir, pdf_dir, stem, width = 6.8, height = 4.8)
  save_plotly_html(p5, plotly_dir, stem, width_in = 6.8, height_in = 4.8)
  plotly_files <- c(plotly_files, paste0(stem, ".plotly.html"))
}

# ---------------------------- Provenance (manifest + inventory) ----------------------------

manifest <- list(
  script_name = script_id$script_name,
  script_path = script_id$script_path,
  script_full = script_id$script_full,
  input_paths = c(input_csv),
  input_mode = input_mode,
  output_dir = norm_abs(output_dir),
  timestamp = timestamp,
  feature_col = feature_col %||% NA_character_,
  group_col = group_col %||% NA_character_,
  sample_id_col = sample_id_col %||% NA_character_,
  day_col = day_col %||% NA_character_,
  value_col = value_col %||% NA_character_,
  sample_naming_rule = "DAY_REP: Day = substring before first underscore; Rep = substring after first underscore",
  cutoff_enabled = do_cutoff,
  cutoff_N = cutoff_N %||% NA_integer_,
  plots_selected = paste(want, collapse = ","),
  legend_pos = legend_pos
)

manifest_path <- write_manifest(output_dir, manifest)
inventory_path <- write_inventory(output_dir)

# capture script header (best effort)
script_header_lines <- c("Script header not available (script path not detected).")
if (!is.na(script_id$script_full) && file.exists(script_id$script_full)) {
  hdr <- readLines(script_id$script_full, warn = FALSE)
  script_header_lines <- hdr[seq_len(min(80, length(hdr)))]
}

rep <- write_quarto_report(
  output_dir = output_dir,
  title = paste0("Publication Plotter Report: ", run_name),
  script_header_lines = script_header_lines,
  manifest_path = manifest_path,
  inventory_path = inventory_path,
  plots_dir = plots_dir,
  tables_dir = tables_dir
)

cat("\nReport written:\n  ", rep$html, "\nMode: ", rep$mode, "\n", sep = "")
cat("\nManifest:\n  ", manifest_path, "\n", sep = "")
cat("Inventory:\n  ", inventory_path, "\n", sep = "")

if (interactive() && file.exists(rep$html) && choose_yes_no("Open report now in browser?", default = "y")) {
  browseURL(rep$html)
}

cat("\nDone.\n")
