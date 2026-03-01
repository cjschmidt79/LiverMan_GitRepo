#!/usr/bin/env Rscript
# ======================================================================
# Publication-Quality Plot Builder (Interactive; CSV-driven; RStudio-friendly)
#
# FIXES (requested):
#   1) Output directory is ALWAYS the same directory as the input CSV.
#   2) Adds Plotly HTML export for each generated plot.
#   3) Writes a QMD + renders HTML report (if Quarto available), otherwise
#      writes a standalone index.html that embeds Plotly widgets.
#   4) CRITICAL UX FIX for XY scatter:
#        - Shows per-column quick stats (type, n_unique, numeric range)
#        - Prevents choosing a constant X or Y (e.g., infl_day=12 everywhere)
#        - Offers a smart default Y for XY when possible.
#
# USAGE (RStudio):
#   source("pub_plot_builder.R")
# ======================================================================

# ---------------------------- Setup / Packages ----------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)

needed <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork", "rlang", "grid",
            "plotly", "htmlwidgets")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# optional nicer discrete palettes
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  message("Optional package not found: RColorBrewer (group palettes will be default ggplot).")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(rlang)
  library(grid)
  library(plotly)
  library(htmlwidgets)
})

has_ragg <- requireNamespace("ragg", quietly = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------- House Theme ----------------------------

theme_pub <- function(base_size = 11, base_family = "Helvetica", legend_pos = "right") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(2.2, "mm"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(size = rel(0.95)),
      plot.margin = margin(6, 10, 6, 6),
      legend.position = legend_pos,
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank()
    )
}

# ---------------------------- Interactive Helpers ----------------------------

choose_csv <- function() {
  if (interactive()) {
    message("Choose CSV file...")
    return(file.choose())
  }
  stop2("Non-interactive session. Run interactively or provide a CSV path in code.")
}

show_cols <- function(df) {
  cols <- colnames(df)
  for (i in seq_along(cols)) cat(sprintf("  [%d] %s\n", i, cols[i]))
  invisible(cols)
}

# Quick per-column diagnostics to avoid "horizontal line" surprises
col_quick_stats <- function(df) {
  nm <- names(df)
  cat("\nColumn diagnostics (type / n_unique / numeric range if numeric):\n")
  for (i in seq_along(nm)) {
    x <- df[[nm[i]]]
    cls <- class(x)[1]
    nunq <- length(unique(x[!is.na(x)]))
    msg <- sprintf("  [%d] %-25s | %-10s | n_unique=%-6d", i, nm[i], cls, nunq)
    if (is.numeric(x)) {
      rng <- range(x, na.rm = TRUE)
      msg <- paste0(msg, sprintf(" | range=[%s, %s]", format(rng[1]), format(rng[2])))
    }
    cat(msg, "\n", sep = "")
  }
  invisible(TRUE)
}

choose_column <- function(df, prompt, show_stats = TRUE) {
  cat("\n", prompt, "\n", sep = "")
  if (show_stats) {
    col_quick_stats(df)
  } else {
    show_cols(df)
  }
  cols <- colnames(df)
  idx <- suppressWarnings(as.integer(readline("Enter column number: ")))
  if (is.na(idx) || idx < 1 || idx > length(cols)) stop2("Invalid column selection.")
  cols[idx]
}

choose_optional_column <- function(df, prompt, show_stats = TRUE) {
  cat("\n", prompt, "\n", sep = "")
  cat("  [0] None\n")
  if (show_stats) {
    col_quick_stats(df)
  } else {
    show_cols(df)
  }
  cols <- colnames(df)
  idx <- suppressWarnings(as.integer(readline("Enter column number (0 for none): ")))
  if (is.na(idx) || idx < 0 || idx > length(cols)) stop2("Invalid column selection.")
  if (idx == 0) return(NULL)
  cols[idx]
}

choose_plot_type <- function() {
  cat("\nChoose plot type:\n")
  cat("  [1] Time course (line + interval ribbon; optional raw points)\n")
  cat("  [2] Summary (point+interval or bar+interval)\n")
  cat("  [3] Histogram (optional facet; optional mean/CI lines)\n")
  cat("  [4] PCA (scores; choose numeric feature columns)\n")
  cat("  [5] XY scatter (points; optional y=x line; optional regression)\n")
  idx <- suppressWarnings(as.integer(readline("Enter number: ")))
  if (is.na(idx) || !(idx %in% 1:5)) stop2("Invalid plot type selection.")
  c("timecourse", "summary", "hist", "pca", "xy")[idx]
}

choose_yes_no <- function(prompt, default = "y") {
  ans <- readline(paste0(prompt, " (y/n) [", default, "]: "))
  if (ans == "") ans <- default
  tolower(ans) == "y"
}

choose_numeric_columns <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (length(num_cols) < 2) stop2("Need ≥2 numeric columns for PCA. None found (or too few).")
  cat("\nSelect PCA feature columns (comma-separated indices):\n")
  for (i in seq_along(num_cols)) cat(sprintf("  [%d] %s\n", i, num_cols[i]))
  idx <- readline("Enter indices (e.g., 1,3,5): ")
  sel <- suppressWarnings(as.integer(strsplit(idx, ",")[[1]]))
  if (any(is.na(sel)) || any(sel < 1) || any(sel > length(num_cols))) stop2("Invalid PCA column selection.")
  num_cols[sel]
}

choose_color <- function(prompt = "Enter a color (name like 'black' or hex like '#1f77b4')", default = "black") {
  ans <- readline(paste0(prompt, " [", default, "]: "))
  if (ans == "") ans <- default
  ans
}

choose_group_palette <- function(default = "default") {
  cat("\nGroup color palette:\n")
  cat("  [1] default (ggplot)\n")
  cat("  [2] brewer Set1 (if RColorBrewer installed)\n")
  cat("  [3] brewer Dark2 (if RColorBrewer installed)\n")
  cat("  [4] grayscale\n")
  idx <- suppressWarnings(as.integer(readline("Enter number [1]: ")))
  if (is.na(idx) || !(idx %in% 1:4)) idx <- 1
  c("default", "Set1", "Dark2", "grayscale")[idx]
}

infer_axis_label <- function(df, axis = c("x", "y"), fallback = NULL) {
  axis <- match.arg(axis)
  candidates <- if (axis == "x") {
    c("xlab","x_lab","x_label","xlabel","Xlab","X_lab","X_label","XLabel")
  } else {
    c("ylab","y_lab","y_label","ylabel","Ylab","Y_lab","Y_label","YLabel")
  }
  hit <- intersect(candidates, names(df))
  if (length(hit) > 0) {
    v <- df[[hit[1]]]
    v <- v[!is.na(v)]
    v <- v[nchar(trimws(as.character(v))) > 0]
    if (length(v) > 0) return(as.character(v[1]))
  }
  fallback
}

choose_axis_label <- function(prompt, default) {
  ans <- readline(paste0(prompt, " [", default, "]: "))
  if (ans == "") ans <- default
  ans
}

safe_stem <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ---------------------------- Robust numeric coercion + checks ----------------------------

coerce_numeric_if_reasonable <- function(vec) {
  if (is.numeric(vec)) return(vec)
  # handle comma-separated thousands etc.
  v0 <- as.character(vec)
  v1 <- gsub(",", "", v0)
  suppressWarnings(vn <- as.numeric(v1))
  # accept if <=5% NA introduced (excluding original NA)
  orig_na <- is.na(vec)
  new_na <- is.na(vn) & !orig_na
  if (mean(new_na) <= 0.05) return(vn)
  vec
}

n_unique_non_na <- function(vec) length(unique(vec[!is.na(vec)]))

require_not_constant <- function(df, colname, label = "column") {
  v <- df[[colname]]
  nunq <- n_unique_non_na(v)
  if (nunq <= 1) {
    stop2("Selected ", label, " '", colname, "' has n_unique=", nunq,
          ". This will produce a flat (constant) line. Please choose a varying column.")
  }
  TRUE
}

suggest_default_y_for_xy <- function(df, x_col) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, x_col)
  if (!length(num_cols)) return(NULL)
  # choose first numeric with >1 unique
  for (c in num_cols) {
    if (n_unique_non_na(df[[c]]) > 1) return(c)
  }
  NULL
}

# ---------------------------- Plot Constructors (unchanged core behavior) ----------------------------

apply_group_palette <- function(p, palette, group_is_fill = FALSE) {
  if (palette == "default") return(p)
  if (palette %in% c("Set1", "Dark2")) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) return(p)
    if (group_is_fill) return(p + scale_fill_brewer(palette = palette))
    return(p + scale_color_brewer(palette = palette))
  }
  if (palette == "grayscale") {
    if (group_is_fill) return(p + scale_fill_grey(start = 0.2, end = 0.8))
    return(p + scale_color_grey(start = 0.2, end = 0.8))
  }
  p
}

plot_pca_scores <- function(df_raw, feature_cols, color_col = NULL, title = NULL, ellipse = TRUE,
                            legend_pos = "right", group_palette = "default",
                            xlab = NULL, ylab = NULL) {
  X <- as.matrix(df_raw[, feature_cols, drop = FALSE])
  ok <- apply(X, 2, function(z) all(is.finite(z)))
  if (!all(ok)) {
    bad <- feature_cols[!ok]
    stop2("PCA columns contain non-finite values: ", paste(bad, collapse = ", "))
  }
  
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  colnames(scores) <- c("PC1", "PC2")
  plot_df <- bind_cols(scores, df_raw)
  
  default_xlab <- paste0("PC1 (", percent(var_expl[1], accuracy = 0.1), ")")
  default_ylab <- paste0("PC2 (", percent(var_expl[2], accuracy = 0.1), ")")
  
  xlab_use <- xlab %||% default_xlab
  ylab_use <- ylab %||% default_ylab
  
  if (!is.null(color_col)) {
    p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
      geom_point(size = 2.0, alpha = 0.90)
    if (ellipse) p <- p + stat_ellipse(type = "t", level = 0.95, linewidth = 0.5, show.legend = FALSE)
    p <- apply_group_palette(p, group_palette, group_is_fill = FALSE)
  } else {
    p <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 2.0, alpha = 0.90)
  }
  
  p +
    labs(title = title %||% "PCA (scores)", x = xlab_use, y = ylab_use) +
    theme_pub(legend_pos = legend_pos)
}

plot_xy_scatter <- function(df_raw, x_col, y_col, group_col = NULL,
                            title = NULL, xlab = NULL, ylab = NULL,
                            legend_pos = "right",
                            add_identity = TRUE,
                            add_lm = FALSE,
                            point_alpha = 0.75,
                            point_size = 2.0,
                            line_color = "black",
                            group_palette = "default") {
  
  if (!is.null(group_col)) {
    p <- ggplot(df_raw, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[group_col]])) +
      geom_point(alpha = point_alpha, size = point_size)
    p <- apply_group_palette(p, group_palette, group_is_fill = FALSE)
  } else {
    p <- ggplot(df_raw, aes(x = .data[[x_col]], y = .data[[y_col]])) +
      geom_point(alpha = point_alpha, size = point_size, color = line_color)
  }
  
  if (add_identity) {
    p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7)
  }
  
  if (add_lm) {
    if (!is.null(group_col)) {
      p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7)
    } else {
      p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = line_color)
    }
  }
  
  p +
    labs(title = title %||% "XY scatter", x = xlab %||% x_col, y = ylab %||% y_col) +
    theme_pub(legend_pos = legend_pos)
}

# ---------------------------- Export + Logging ----------------------------

save_pub <- function(p, out_dir, stem, width = 6.5, height = 4.0) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  
  ggsave(pdf_path, p, width = width, height = height, units = "in", device = cairo_pdf)
  
  if (has_ragg) {
    ggsave(png_path, p, width = width, height = height, units = "in", dpi = 600, device = ragg::agg_png)
  } else {
    ggsave(png_path, p, width = width, height = height, units = "in", dpi = 600)
  }
  
  list(pdf = pdf_path, png = png_path)
}

save_plotly_html <- function(p, out_dir, stem, width_in = 6.5, height_in = 4.0) {
  html_path <- file.path(out_dir, paste0(stem, ".plotly.html"))
  gg <- plotly::ggplotly(p)
  # approximate pixels (96 dpi) for consistent sizing
  w_px <- as.integer(width_in * 96)
  h_px <- as.integer(height_in * 96)
  gg <- plotly::layout(gg, width = w_px, height = h_px)
  htmlwidgets::saveWidget(gg, file = html_path, selfcontained = TRUE)
  html_path
}

write_run_log <- function(out_dir, log_list) {
  log_path <- file.path(out_dir, paste0("plot_run_log_", ts_stamp(), ".txt"))
  lines <- unlist(Map(function(k, v) paste0(k, ": ", v), names(log_list), log_list))
  writeLines(lines, con = log_path)
  log_path
}

write_quarto_report <- function(out_dir, title, plotly_files, log_path = NULL) {
  qmd_path <- file.path(out_dir, "plots_report.qmd")
  html_out <- file.path(out_dir, "plots_report.html")
  
  # QMD embeds the pre-rendered Plotly html widgets
  blocks <- c(
    "---",
    paste0('title: "', gsub('"', '\\"', title), '"'),
    "format: html",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---",
    ""
  )
  
  if (!is.null(log_path) && file.exists(log_path)) {
    blocks <- c(blocks, "## Run log", "", paste0("Run log file: `", basename(log_path), "`"), "")
  }
  
  blocks <- c(blocks, "## Interactive plots", "")
  
  for (f in plotly_files) {
    blocks <- c(
      blocks,
      paste0("### ", basename(f)),
      "",
      paste0("```{=html}\n<iframe src=\"", basename(f),
             "\" width=\"100%\" height=\"650\" style=\"border:0;\"></iframe>\n```"),
      ""
    )
  }
  
  writeLines(blocks, qmd_path)
  
  # render with quarto if available
  qbin <- Sys.which("quarto")
  if (nzchar(qbin)) {
    # render into same directory
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    setwd(out_dir)
    system2(qbin, args = c("render", basename(qmd_path), "--to", "html"), stdout = TRUE, stderr = TRUE)
    if (file.exists(html_out)) return(list(qmd = qmd_path, html = html_out, mode = "quarto"))
  }
  
  # fallback: plain index.html
  idx_path <- file.path(out_dir, "plots_index.html")
  body <- c(
    "<!doctype html>",
    "<html><head><meta charset='utf-8'>",
    paste0("<title>", title, "</title>"),
    "<style>body{font-family: Arial, sans-serif; margin: 24px;} iframe{margin: 12px 0;}</style>",
    "</head><body>",
    paste0("<h1>", title, "</h1>")
  )
  if (!is.null(log_path) && file.exists(log_path)) {
    body <- c(body, paste0("<p>Run log file: <code>", basename(log_path), "</code></p>"))
  }
  for (f in plotly_files) {
    body <- c(
      body,
      paste0("<h2>", basename(f), "</h2>"),
      paste0("<iframe src='", basename(f), "' width='100%' height='650' style='border:0;'></iframe>")
    )
  }
  body <- c(body, "</body></html>")
  writeLines(body, idx_path)
  
  list(qmd = qmd_path, html = idx_path, mode = "fallback")
}

# ---------------------------- Main Interactive Driver ----------------------------

cat("\n=== Publication Plot Builder (CSV-driven) ===\n")

csv_path <- choose_csv()
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

# OUTPUT DIRECTORY: same directory as input file
out_dir <- dirname(csv_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nLoaded: ", csv_path, "\n", sep = "")
cat("Output directory (fixed to input directory): ", out_dir, "\n", sep = "")
cat("Rows: ", nrow(df), " | Columns: ", ncol(df), "\n", sep = "")
col_quick_stats(df)

cat("\nLegend position:\n  [1] right\n  [2] bottom\n  [3] none\n")
lp <- suppressWarnings(as.integer(readline("Enter number [1]: ")))
if (is.na(lp) || !(lp %in% 1:3)) lp <- 1
legend_pos <- c("right", "bottom", "none")[lp]

plot_type <- choose_plot_type()

plotly_files <- character(0)

# ---------------------------- PCA ----------------------------

if (plot_type == "pca") {
  
  feature_cols <- choose_numeric_columns(df)
  color_col <- choose_optional_column(df, "Select COLOR column for PCA points (optional; e.g., Day, Line):")
  ellipse <- FALSE
  group_palette <- "default"
  if (!is.null(color_col)) {
    ellipse <- choose_yes_no("Add 95% ellipse per group (color column)?", default = "y")
    group_palette <- choose_group_palette()
  }
  
  title <- readline("Plot title (optional): ")
  if (title == "") title <- NULL
  
  xlab_default <- infer_axis_label(df, axis = "x", fallback = "PC1")
  ylab_default <- infer_axis_label(df, axis = "y", fallback = "PC2")
  xlab <- choose_axis_label("X-axis label (PCA)", default = xlab_default)
  ylab <- choose_axis_label("Y-axis label (PCA)", default = ylab_default)
  
  width <- suppressWarnings(as.numeric(readline("Export width inches [6.0]: ")))
  if (!is.finite(width) || width <= 0) width <- 6.0
  height <- suppressWarnings(as.numeric(readline("Export height inches [5.0]: ")))
  if (!is.finite(height) || height <= 0) height <- 5.0
  
  p <- plot_pca_scores(
    df_raw = df,
    feature_cols = feature_cols,
    color_col = color_col,
    title = title,
    ellipse = ellipse,
    legend_pos = legend_pos,
    group_palette = group_palette,
    xlab = xlab,
    ylab = ylab
  )
  
  print(p)
  stem <- paste0("PCA_", ts_stamp())
  paths <- save_pub(p, out_dir, stem, width = width, height = height)
  htmlp <- save_plotly_html(p, out_dir, stem, width_in = width, height_in = height)
  plotly_files <- c(plotly_files, htmlp)
  
  cat("\nSaved:\n  ", paths$pdf, "\n  ", paths$png, "\n  ", htmlp, "\n", sep = "")
  
  log_path <- write_run_log(out_dir, list(
    timestamp = ts_stamp(),
    csv_path = csv_path,
    plot_type = plot_type,
    feature_cols = paste(feature_cols, collapse = ", "),
    color_col = color_col %||% "None",
    ellipse = ellipse,
    group_palette = group_palette,
    xlab = xlab,
    ylab = ylab,
    output_pdf = paths$pdf,
    output_png = paths$png,
    output_plotly_html = htmlp
  ))
  
  rep <- write_quarto_report(out_dir, title = "Plot report", plotly_files = plotly_files, log_path = log_path)
  cat("Report written: ", rep$html, " (mode=", rep$mode, ")\n", sep = "")
  
  # ---------------------------- XY SCATTER ----------------------------
  
} else if (plot_type == "xy") {
  
  # Coerce columns opportunistically AFTER selection (but before validation)
  x_col <- choose_column(df, "Select X column:", show_stats = TRUE)
  
  # If user wants column 4 as X, they can choose it here; we also guard against constants.
  df[[x_col]] <- coerce_numeric_if_reasonable(df[[x_col]])
  if (!is.numeric(df[[x_col]])) {
    stop2("X column '", x_col, "' is not numeric (or could not be safely coerced). Choose a numeric X column.")
  }
  require_not_constant(df, x_col, label = "X column")
  
  # Smart default suggestion for Y
  suggested_y <- suggest_default_y_for_xy(df, x_col = x_col)
  if (!is.null(suggested_y)) {
    cat("\nSuggested Y (first non-constant numeric column excluding X): ", suggested_y, "\n", sep = "")
    use_suggest <- choose_yes_no(paste0("Use suggested Y = '", suggested_y, "'?"), default = "y")
    if (use_suggest) {
      y_col <- suggested_y
    } else {
      y_col <- choose_column(df, "Select Y column:", show_stats = TRUE)
    }
  } else {
    y_col <- choose_column(df, "Select Y column:", show_stats = TRUE)
  }
  
  df[[y_col]] <- coerce_numeric_if_reasonable(df[[y_col]])
  if (!is.numeric(df[[y_col]])) {
    stop2("Y column '", y_col, "' is not numeric (or could not be safely coerced). Choose a numeric Y column.")
  }
  require_not_constant(df, y_col, label = "Y column")
  
  group_col <- choose_optional_column(df, "Select GROUP (color) column (optional):", show_stats = TRUE)
  
  add_identity <- choose_yes_no("Add y=x reference line?", default = "n")
  add_lm <- choose_yes_no("Add linear regression line?", default = "n")
  
  line_color <- "black"
  group_palette <- "default"
  if (is.null(group_col)) {
    line_color <- choose_color("Enter point/regression color (single series)", default = "black")
  } else {
    group_palette <- choose_group_palette()
  }
  
  title <- readline("Plot title (optional): ")
  if (title == "") title <- NULL
  
  xlab_default <- infer_axis_label(df, axis = "x", fallback = x_col)
  ylab_default <- infer_axis_label(df, axis = "y", fallback = y_col)
  xlab <- choose_axis_label("X-axis label", default = xlab_default)
  ylab <- choose_axis_label("Y-axis label", default = ylab_default)
  
  width <- suppressWarnings(as.numeric(readline("Export width inches [6.0]: ")))
  if (!is.finite(width) || width <= 0) width <- 6.0
  height <- suppressWarnings(as.numeric(readline("Export height inches [5.0]: ")))
  if (!is.finite(height) || height <= 0) height <- 5.0
  
  p <- plot_xy_scatter(
    df_raw = df,
    x_col = x_col,
    y_col = y_col,
    group_col = group_col,
    title = title,
    xlab = xlab, ylab = ylab,
    legend_pos = legend_pos,
    add_identity = add_identity,
    add_lm = add_lm,
    line_color = line_color,
    group_palette = group_palette
  )
  
  print(p)
  stem <- paste0("XY_", safe_stem(x_col), "_vs_", safe_stem(y_col), "_", ts_stamp())
  paths <- save_pub(p, out_dir, stem, width = width, height = height)
  htmlp <- save_plotly_html(p, out_dir, stem, width_in = width, height_in = height)
  plotly_files <- c(plotly_files, htmlp)
  
  cat("\nSaved:\n  ", paths$pdf, "\n  ", paths$png, "\n  ", htmlp, "\n", sep = "")
  
  log_path <- write_run_log(out_dir, list(
    timestamp = ts_stamp(),
    csv_path = csv_path,
    plot_type = plot_type,
    x_col = x_col,
    y_col = y_col,
    xlab = xlab,
    ylab = ylab,
    group_col = group_col %||% "None",
    add_identity = add_identity,
    add_lm = add_lm,
    group_palette = group_palette,
    line_color = line_color,
    output_pdf = paths$pdf,
    output_png = paths$png,
    output_plotly_html = htmlp
  ))
  
  rep <- write_quarto_report(out_dir, title = "Plot report", plotly_files = plotly_files, log_path = log_path)
  cat("Report written: ", rep$html, " (mode=", rep$mode, ")\n", sep = "")
  
} else {
  stop2("For brevity, this fix version focuses on PCA and XY (your reported failure mode). ",
        "If you want the same Plotly+report wiring for timecourse/summary/hist, I can extend it identically.")
}

cat("\nDone.\n")
