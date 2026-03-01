#!/usr/bin/env Rscript
# ======================================================================
# Publication-Quality Plot Builder (Interactive; CSV-driven; RStudio-friendly)
#
# FIX (this issue)
#   The previous UX was confusing for tables like:
#     Day | Metric | lo | mid | hi
#   because you were asked to pick Y columns from (lo, mid, hi) even though
#   lo/hi are interval bounds, not the plotted value.
#
#   This version:
#     1) Detects common "precomputed CI" patterns (lo/mid/hi or lo/mean/hi).
#     2) Offers a simple guided choice:
#          "Use mid as Y and lo/hi as bounds automatically? (y/n)"
#        If yes: it sets uncertainty mode to precomputed and DOES NOT make you
#        manually map mean/lo/hi repeatedly.
#     3) Adds an explicit "Metric plotting mode":
#          - One combined plot with all metrics (colored)
#          - Separate plot per metric (one file per metric)
#
# USAGE (RStudio):
#   source("pub_plot_builder.R")
# ======================================================================

# ---------------------------- Setup / Packages ----------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)

needed <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork", "rlang", "grid")
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
})

has_ragg <- requireNamespace("ragg", quietly = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------- House Theme (Improved) ----------------------------

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

choose_out_dir <- function(default_dir = getwd()) {
  if (!interactive()) return(default_dir)
  cat("\nOutput directory:\n")
  cat("  [1] Use working directory: ", default_dir, "\n", sep = "")
  cat("  [2] Choose another directory (type path)\n")
  idx <- readline("Enter choice (1/2): ")
  if (idx == "2") {
    p <- readline("Enter output directory path: ")
    if (!dir.exists(p)) {
      ans <- readline("Directory does not exist. Create it? (y/n): ")
      if (tolower(ans) != "y") stop2("Output directory not created; aborting.")
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
    return(p)
  }
  default_dir
}

show_cols <- function(df) {
  cols <- colnames(df)
  for (i in seq_along(cols)) cat(sprintf("  [%d] %s\n", i, cols[i]))
  invisible(cols)
}

choose_column <- function(df, prompt) {
  cat("\n", prompt, "\n", sep = "")
  cols <- colnames(df)
  show_cols(df)
  idx <- suppressWarnings(as.integer(readline("Enter column number: ")))
  if (is.na(idx) || idx < 1 || idx > length(cols)) stop2("Invalid column selection.")
  cols[idx]
}

choose_optional_column <- function(df, prompt) {
  cat("\n", prompt, "\n", sep = "")
  cat("  [0] None\n")
  cols <- colnames(df)
  show_cols(df)
  idx <- suppressWarnings(as.integer(readline("Enter column number (0 for none): ")))
  if (is.na(idx) || idx < 0 || idx > length(cols)) stop2("Invalid column selection.")
  if (idx == 0) return(NULL)
  cols[idx]
}

choose_multi_columns <- function(df, prompt, default_first = TRUE) {
  cat("\n", prompt, "\n", sep = "")
  cols <- colnames(df)
  show_cols(df)
  cat("\nType ONE index (e.g., 5) or a comma list (e.g., 5,6,9) then press Enter.\n")
  raw <- trimws(readline("Indices: "))
  if (!nzchar(raw)) {
    if (default_first) return(cols[1])
    stop2("No selection provided.")
  }
  parts <- unlist(strsplit(raw, ","))
  sel <- suppressWarnings(as.integer(trimws(parts)))
  sel <- sel[!is.na(sel)]
  if (!length(sel)) stop2("No valid indices provided.")
  if (any(sel < 1) || any(sel > length(cols))) stop2("Index out of range.")
  unique(cols[sel])
}

choose_plot_type <- function() {
  cat("\nChoose plot type:\n")
  cat("  [1] Time course (line + interval ribbon/error bars; optional raw points)\n")
  cat("  [2] Summary (point+interval or bar+interval)\n")
  cat("  [3] Histogram (optional facet; optional mean/CI lines)\n")
  cat("  [4] PCA (scores; choose numeric feature columns)\n")
  cat("  [5] XY scatter (points; optional y=x line; optional regression)\n")
  idx <- suppressWarnings(as.integer(readline("Enter number: ")))
  if (is.na(idx) || !(idx %in% 1:5)) stop2("Invalid plot type selection.")
  c("timecourse", "summary", "hist", "pca", "xy")[idx]
}

choose_uncertainty_mode <- function() {
  cat("\nUncertainty / error bars:\n")
  cat("  [0] None\n")
  cat("  [1] Standard error (SE)\n")
  cat("  [2] 95% CI (t-based)\n")
  cat("  [3] 95% CI (bootstrap)\n")
  cat("  [4] Use precomputed columns (you will pick mean/lo/hi, OR auto-detect)\n")
  idx <- suppressWarnings(as.integer(readline("Enter number: ")))
  if (is.na(idx) || !(idx %in% 0:4)) stop2("Invalid uncertainty selection.")
  idx
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

choose_precomputed_ci_cols <- function(df) {
  mean_col <- choose_column(df, "Pick column for mean/central estimate (this is the plotted Y):")
  lo_col   <- choose_column(df, "Pick column for lower bound (lo):")
  hi_col   <- choose_column(df, "Pick column for upper bound (hi):")
  list(mean = mean_col, lo = lo_col, hi = hi_col)
}

choose_color <- function(prompt = "Enter a line color (name like 'black' or hex like '#1f77b4')", default = "black") {
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

choose_output_mode <- function(context = "timecourse") {
  cat("\nOutput mode (", context, "):\n", sep = "")
  cat("  [1] Composite only (patchwork)\n")
  cat("  [2] Separate plots only (one file per panel)\n")
  cat("  [3] Both composite + separate\n")
  idx <- suppressWarnings(as.integer(readline("Enter number [2]: ")))
  if (is.na(idx) || !(idx %in% 1:3)) idx <- 2
  c("composite", "separate", "both")[idx]
}

# NEW: for a GROUP like Metric, ask if user wants one plot or separate per metric
choose_group_plot_mode <- function(group_col) {
  if (is.null(group_col)) return("together")
  cat("\nGroup plotting for ", group_col, ":\n", sep = "")
  cat("  [1] One plot with all group levels (colored lines/points)\n")
  cat("  [2] Separate plot for each ", group_col, " level (one file per level)\n", sep = "")
  idx <- suppressWarnings(as.integer(readline("Enter number [2]: ")))
  if (is.na(idx) || !(idx %in% 1:2)) idx <- 2
  c("together", "separate")[idx]
}

# ---- Axis label helpers ----

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

# ---- Precomputed CI auto-detection (KEY FIX) ----

norm_names <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

detect_precomputed_ci_triplet <- function(df) {
  n <- names(df)
  nn <- norm_names(n)
  
  # candidate patterns: lo / mid / hi, or lo / mean / hi, or lower/upper, etc.
  find_one <- function(opts) {
    hits <- which(nn %in% opts)
    if (length(hits) >= 1) n[hits[1]] else NULL
  }
  
  lo <- find_one(c("lo","lower","low","lcl","cilow","cilower","lowerci","ci_lower","ci_lo"))
  hi <- find_one(c("hi","upper","high","ucl","cihigh","ciupper","upperci","ci_upper","ci_hi"))
  mid <- find_one(c("mid","med","mean","mu","center","central","estimate","est","value"))
  
  
  # common exact: lo/mid/hi
  # if mid missing but mean exists, above captures it
  list(lo = lo, mid = mid, hi = hi)
}

safe_stem <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ---------------------------- Interval / Summary Functions ----------------------------

se_mean <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(mean = mean(x), se = NA_real_, n = n))
  m <- mean(x)
  se <- stats::sd(x) / sqrt(n)
  c(mean = m, se = se, n = n)
}

ci_t <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(mean = mean(x), lo = NA_real_, hi = NA_real_, n = n))
  m <- mean(x)
  se <- stats::sd(x) / sqrt(n)
  alpha <- 1 - conf
  tcrit <- stats::qt(1 - alpha/2, df = n - 1)
  c(mean = m, lo = m - tcrit * se, hi = m + tcrit * se, n = n)
}

ci_boot <- function(x, conf = 0.95, R = 2000, seed = 1) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(mean = mean(x), lo = NA_real_, hi = NA_real_, n = n))
  set.seed(seed)
  boots <- replicate(R, mean(sample(x, size = n, replace = TRUE)))
  alpha <- 1 - conf
  qs <- stats::quantile(boots, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
  c(mean = mean(x), lo = qs[1], hi = qs[2], n = n)
}

summarise_with_uncertainty <- function(df, x_col, y_col, group_col = NULL,
                                       mode = 2, boot_R = 2000, boot_seed = 1) {
  stopifnot(x_col %in% names(df), y_col %in% names(df))
  group_vars <- c(x_col, if (!is.null(group_col)) group_col)
  
  out <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      stats = list({
        v <- .data[[y_col]]
        if (mode == 1) se_mean(v)
        else if (mode == 2) ci_t(v, conf = 0.95)
        else if (mode == 3) ci_boot(v, conf = 0.95, R = boot_R, seed = boot_seed)
        else c(mean = mean(v, na.rm = TRUE), lo = NA_real_, hi = NA_real_, n = sum(is.finite(v)))
      }),
      .groups = "drop"
    )
  
  out$mean <- vapply(out$stats, function(z) unname(z["mean"]), numeric(1))
  out$n    <- vapply(out$stats, function(z) unname(z["n"]), numeric(1))
  
  if (mode == 1) {
    out$se <- vapply(out$stats, function(z) unname(z["se"]), numeric(1))
  } else if (mode %in% c(2, 3)) {
    out$lo <- vapply(out$stats, function(z) unname(z["lo"]), numeric(1))
    out$hi <- vapply(out$stats, function(z) unname(z["hi"]), numeric(1))
  } else {
    out$lo <- NA_real_
    out$hi <- NA_real_
  }
  
  out$stats <- NULL
  out
}

# ---------------------------- Plot Constructors ----------------------------

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

plot_timecourse <- function(df_raw, df_sum, x_col, y_col, group_col = NULL,
                            mode = 2, show_points = TRUE, title = NULL,
                            xlab = NULL, ylab = NULL,
                            legend_pos = "right",
                            line_color = "black",
                            group_palette = "default") {
  
  aes_base <- aes(x = .data[[x_col]], y = mean)
  if (!is.null(group_col)) aes_base <- aes(x = .data[[x_col]], y = mean, group = .data[[group_col]])
  
  p <- ggplot(df_sum, aes_base)
  
  if (mode %in% c(2, 3, 4)) {
    if (!is.null(group_col)) {
      p <- p + geom_ribbon(aes(ymin = lo, ymax = hi, fill = .data[[group_col]]),
                           alpha = 0.18, color = NA, show.legend = FALSE)
    } else {
      p <- p + geom_ribbon(aes(ymin = lo, ymax = hi),
                           alpha = 0.20, color = NA)
    }
  }
  
  if (!is.null(group_col)) {
    p <- p +
      geom_line(aes(color = .data[[group_col]]), linewidth = 0.75) +
      geom_point(aes(color = .data[[group_col]]), size = 1.8)
    p <- apply_group_palette(p, group_palette, group_is_fill = FALSE)
  } else {
    p <- p +
      geom_line(linewidth = 0.8, color = line_color) +
      geom_point(size = 2.0, color = line_color)
  }
  
  if (mode == 1 && "se" %in% names(df_sum)) {
    p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                           width = 0.12, linewidth = 0.45)
  }
  
  if (show_points) {
    if (!is.null(group_col)) {
      p <- p + geom_point(
        data = df_raw,
        aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[group_col]]),
        inherit.aes = FALSE, alpha = 0.25, size = 1.1
      )
      p <- apply_group_palette(p, group_palette, group_is_fill = FALSE)
    } else {
      p <- p + geom_point(
        data = df_raw,
        aes(x = .data[[x_col]], y = .data[[y_col]]),
        inherit.aes = FALSE, alpha = 0.25, size = 1.1, color = line_color
      )
    }
  }
  
  x_vals <- df_sum[[x_col]]
  if (is.numeric(x_vals)) {
    breaks <- sort(unique(x_vals))
    if (length(breaks) <= 25) {
      p <- p + scale_x_continuous(breaks = breaks, minor_breaks = NULL)
    }
  }
  
  p +
    labs(title = title %||% "Time course",
         x = xlab %||% x_col,
         y = ylab %||% y_col,
         color = group_col %||% NULL) +
    theme_pub(legend_pos = legend_pos)
}

plot_summary <- function(df_raw, df_sum, x_col, y_col, group_col = NULL,
                         mode = 2, as_bar = FALSE, title = NULL,
                         xlab = NULL, ylab = NULL,
                         legend_pos = "right",
                         group_palette = "default") {
  
  dodge <- position_dodge(width = 0.6)
  
  if (!is.null(group_col)) {
    if (as_bar) {
      p <- ggplot(df_sum, aes(x = .data[[x_col]], y = mean, fill = .data[[group_col]])) +
        geom_col(position = dodge, width = 0.72)
      if (mode %in% c(2, 3, 4)) {
        p <- p + geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.14, linewidth = 0.45)
      } else if (mode == 1 && "se" %in% names(df_sum)) {
        p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = dodge, width = 0.14, linewidth = 0.45)
      }
      p <- apply_group_palette(p, group_palette, group_is_fill = TRUE)
    } else {
      p <- ggplot(df_sum, aes(x = .data[[x_col]], y = mean, color = .data[[group_col]])) +
        geom_point(position = dodge, size = 2.1)
      if (mode %in% c(2, 3, 4)) {
        p <- p + geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.14, linewidth = 0.45)
      } else if (mode == 1 && "se" %in% names(df_sum)) {
        p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = dodge, width = 0.14, linewidth = 0.45)
      }
      p <- apply_group_palette(p, group_palette, group_is_fill = FALSE)
    }
  } else {
    if (as_bar) {
      p <- ggplot(df_sum, aes(x = .data[[x_col]], y = mean)) +
        geom_col(width = 0.72)
      if (mode %in% c(2, 3, 4)) {
        p <- p + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.45)
      } else if (mode == 1 && "se" %in% names(df_sum)) {
        p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.14, linewidth = 0.45)
      }
    } else {
      p <- ggplot(df_sum, aes(x = .data[[x_col]], y = mean)) +
        geom_point(size = 2.1)
      if (mode %in% c(2, 3, 4)) {
        p <- p + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.45)
      } else if (mode == 1 && "se" %in% names(df_sum)) {
        p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.14, linewidth = 0.45)
      }
    }
  }
  
  x_vals <- df_sum[[x_col]]
  if (is.numeric(x_vals)) {
    breaks <- sort(unique(x_vals))
    if (length(breaks) <= 25) {
      p <- p + scale_x_continuous(breaks = breaks, minor_breaks = NULL)
    }
  }
  
  p +
    labs(title = title %||% "Summary",
         x = xlab %||% x_col,
         y = ylab %||% y_col) +
    theme_pub(legend_pos = legend_pos) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_hist <- function(df_raw, y_col, facet_col = NULL, show_ci_lines = TRUE,
                      ci_mode = 2, boot_R = 2000, boot_seed = 1,
                      bins = 30, title = NULL, xlab = NULL,
                      legend_pos = "right") {
  
  p <- ggplot(df_raw, aes(x = .data[[y_col]])) +
    geom_histogram(bins = bins)
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(vars(.data[[facet_col]]), scales = "free_y")
  } else if (show_ci_lines) {
    v <- df_raw[[y_col]]
    stats <- if (ci_mode == 1) se_mean(v)
    else if (ci_mode == 2) ci_t(v, conf = 0.95)
    else if (ci_mode == 3) ci_boot(v, conf = 0.95, R = boot_R, seed = boot_seed)
    else NULL
    
    if (!is.null(stats)) {
      p <- p + geom_vline(xintercept = unname(stats["mean"]), linewidth = 0.7)
      if (ci_mode == 1 && "se" %in% names(stats)) {
        p <- p + geom_vline(
          xintercept = unname(stats["mean"]) + c(-1, 1) * unname(stats["se"]),
          linetype = "dashed", linewidth = 0.55
        )
      }
      if (ci_mode %in% c(2, 3) && all(c("lo", "hi") %in% names(stats))) {
        p <- p + geom_vline(
          xintercept = unname(stats[c("lo", "hi")]),
          linetype = "dashed", linewidth = 0.55
        )
      }
    }
  }
  
  p +
    labs(title = title %||% "Histogram", x = xlab %||% y_col, y = "Count") +
    theme_pub(legend_pos = legend_pos)
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

write_run_log <- function(out_dir, log_list) {
  log_path <- file.path(out_dir, paste0("plot_run_log_", ts_stamp(), ".txt"))
  lines <- unlist(Map(function(k, v) paste0(k, ": ", v), names(log_list), log_list))
  writeLines(lines, con = log_path)
  log_path
}

# ---------------------------- Main Interactive Driver ----------------------------

cat("\n=== Publication Plot Builder (CSV-driven) ===\n")

csv_path <- choose_csv()
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

# ---- NSV wide-CI helper: detect columns like D_med/D_lo/D_hi, ED_med/ED_lo/ED_hi ----
is_nsv_wide_ci <- function(df) {
  n <- names(df)
  any(grepl("_(med|mid|mean|lo|hi)$", tolower(n))) && any(grepl("_(lo|hi)$", tolower(n)))
}

pivot_nsv_wide_to_long <- function(df, day_candidates = c("Day","day","time","Time")) {
  # pick Day column if present
  day_col <- intersect(day_candidates, names(df))
  if (length(day_col) == 0) stop("Could not find a Day/time column for NSV pivot.")
  day_col <- day_col[1]
  
  # keep only columns with suffix _med/_mid/_mean/_lo/_hi plus day
  keep <- c(day_col, names(df)[grepl("_(med|mid|mean|lo|hi)$", tolower(names(df)))])
  df2 <- df[, keep, drop = FALSE]
  
  # pivot longer: Metric + stat
  df_long <- df2 |>
    tidyr::pivot_longer(
      cols = -all_of(day_col),
      names_to = c("Metric","stat"),
      names_pattern = "^(.*)_(med|mid|mean|lo|hi)$",
      values_to = "value"
    ) |>
    tidyr::pivot_wider(names_from = stat, values_from = value)
  
  # normalize central column name to 'mid'
  if (!"mid" %in% names(df_long) && "med" %in% names(df_long)) df_long$mid <- df_long$med
  if (!"mid" %in% names(df_long) && "mean" %in% names(df_long)) df_long$mid <- df_long$mean
  
  # standardize Day column name
  names(df_long)[names(df_long) == day_col] <- "Day"
  
  # final columns: Day, Metric, mid, lo, hi
  df_long <- df_long[, intersect(c("Day","Metric","mid","lo","hi"), names(df_long)), drop = FALSE]
  df_long
}

# Offer pivot if wide CI looks like NSV
if (is_nsv_wide_ci(df)) {
  cat("\nDetected wide precomputed CI columns (e.g., D_med/D_lo/D_hi). This looks like an NSV summary.\n")
  if (choose_yes_no("Pivot to long format Day|Metric|mid|lo|hi for easier plotting?", default = "y")) {
    df <- pivot_nsv_wide_to_long(df)
    cat("Pivoted to long format. Columns now:\n")
    show_cols(df)
  }
}


cat("\nLoaded: ", csv_path, "\n", sep = "")
cat("Rows: ", nrow(df), " | Columns: ", ncol(df), "\n", sep = "")
cat("Column preview:\n")
show_cols(df)

cat("\nLegend position:\n  [1] right\n  [2] bottom\n  [3] none\n")
lp <- suppressWarnings(as.integer(readline("Enter number [1]: ")))
if (is.na(lp) || !(lp %in% 1:3)) lp <- 1
legend_pos <- c("right", "bottom", "none")[lp]

out_dir <- choose_out_dir(getwd())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_type <- choose_plot_type()

x_col <- NULL; group_col <- NULL; facet_col <- NULL
unc_mode <- NULL
group_palette <- "default"

# ------------- TIMECOURSE / SUMMARY -------------

if (plot_type %in% c("timecourse", "summary")) {
  
  x_col <- choose_column(df, "Select X column (e.g., Day, Time):")
  
  # Coerce X to numeric if it looks numeric
  if (!is.numeric(df[[x_col]])) {
    suppressWarnings({ x_num <- as.numeric(df[[x_col]]) })
    if (sum(is.na(x_num)) <= 0.05 * length(x_num)) {
      df[[x_col]] <- x_num
    }
  }
  
  group_col <- choose_optional_column(df, "Select GROUP column (optional; e.g., Metric, Line, Treatment):")
  
  # Ask: one plot vs separate plot per group level
  group_plot_mode <- "together"
  if (!is.null(group_col)) {
    group_plot_mode <- choose_group_plot_mode(group_col)
    if (group_plot_mode == "together") {
      group_palette <- choose_group_palette()
    }
  }
  
  # Detect precomputed CI pattern early (so we can guide Y selection)
  trip <- detect_precomputed_ci_triplet(df)
  
  # Choose uncertainty mode (with a guided option)
  unc_mode <- choose_uncertainty_mode()
  
  boot_R <- 2000
  boot_seed <- 1
  if (unc_mode == 3) {
    boot_R <- suppressWarnings(as.integer(readline("Bootstrap iterations R [2000]: ")))
    if (is.na(boot_R) || boot_R < 200) boot_R <- 2000
    boot_seed <- suppressWarnings(as.integer(readline("Bootstrap seed [1]: ")))
    if (is.na(boot_seed)) boot_seed <- 1
  }
  
  use_precomputed <- (unc_mode == 4)
  
  # GUIDED FIX: if precomputed chosen and we detect lo/mid/hi (or lo/mean/hi), offer auto-map
  auto_ci_map <- FALSE
  auto_ci <- list(mean = NULL, lo = NULL, hi = NULL)
  
  if (use_precomputed && !is.null(trip$lo) && !is.null(trip$hi) && !is.null(trip$mid)) {
    cat("\nDetected likely precomputed CI columns:\n")
    cat("  mean/central: ", trip$mid, "\n", sep = "")
    cat("  lo:           ", trip$lo, "\n", sep = "")
    cat("  hi:           ", trip$hi, "\n", sep = "")
    auto_ci_map <- choose_yes_no("Use these automatically? (recommended for Day|Metric|lo|mid|hi tables)", default = "y")
    if (auto_ci_map) {
      auto_ci <- list(mean = trip$mid, lo = trip$lo, hi = trip$hi)
    }
  }
  
  # Y selection:
  # - If precomputed auto-map is ON: choose plotted value column(s) defaults to detected mid/mean.
  # - Otherwise: user chooses one or more VALUE columns to plot (explicitly warns not to choose lo/hi).
  if (use_precomputed && auto_ci_map) {
    cat("\nVALUE columns to plot:\n")
    cat("  Using detected central estimate column: ", auto_ci$mean, "\n", sep = "")
    y_cols <- auto_ci$mean
  } else {
    cat("\nSelect VALUE column(s) to PLOT.\n")
    cat("Important: Do NOT select bound columns like 'lo' or 'hi' as Y.\n")
    y_cols <- choose_multi_columns(df, "Select VALUE column(s) (ONE or multiple):")
    
    # Guardrail: if user picked lo/hi as Y in precomputed mode, warn and allow reselect
    if (use_precomputed) {
      bad <- y_cols[norm_names(y_cols) %in% c("lo","lower","low","hi","upper","high")]
      if (length(bad) > 0) {
        cat("\nYou selected: ", paste(bad, collapse = ", "), "\n", sep = "")
        cat("Those look like interval bounds, not the plotted value column.\n")
        redo <- choose_yes_no("Reselect VALUE column(s)?", default = "y")
        if (redo) {
          y_cols <- choose_multi_columns(df, "Reselect VALUE column(s) to PLOT:")
        }
      }
    }
  }
  
  show_points <- TRUE
  if (plot_type == "timecourse") {
    show_points <- choose_yes_no("Overlay raw points on summary? (timecourse only)", default = "y")
  }
  
  # Color logic
  line_color <- "black"
  if (is.null(group_col) || group_plot_mode == "separate") {
    line_color <- choose_color("Enter line/point color (single series plots)", default = "black")
  }
  
  # Output mode:
  # Only meaningful when there are multiple Y panels (composite stack).
  # Separate-by-group already produces separate files per metric; we only ask about composite stacks if multi-Y.
  out_mode <- "separate"
  if (length(y_cols) > 1) {
    out_mode <- choose_output_mode(context = plot_type)
  } else {
    out_mode <- "separate"
  }
  
  title_global <- readline("Plot title (optional; used for composite annotation): ")
  if (title_global == "") title_global <- NULL
  
  xlab_default <- infer_axis_label(df, axis = "x", fallback = x_col)
  xlab <- choose_axis_label("X-axis label", default = xlab_default)
  
  width <- suppressWarnings(as.numeric(readline("Export width inches [6.5]: ")))
  if (!is.finite(width) || width <= 0) width <- 6.5
  height <- suppressWarnings(as.numeric(readline("Export height inches [4.0]: ")))
  if (!is.finite(height) || height <= 0) height <- 4.0
  
  per_panel_height <- height
  if (length(y_cols) > 1 && out_mode %in% c("composite", "both")) {
    per_panel_height <- suppressWarnings(as.numeric(readline(paste0(
      "Per-panel height inches for composite [", height, "]: "
    ))))
    if (!is.finite(per_panel_height) || per_panel_height <= 0) per_panel_height <- height
  }
  
  # Precomputed CI mapping:
  # - If auto_ci_map: reuse lo/hi for all Y, and mean is the Y itself (so map is per Y).
  # - Else: prompt once per Y for mean/lo/hi.
  ci_cols_map <- list()
  if (use_precomputed) {
    if (auto_ci_map) {
      for (yy in y_cols) {
        ci_cols_map[[yy]] <- list(mean = yy, lo = auto_ci$lo, hi = auto_ci$hi)
      }
    } else {
      for (yy in y_cols) {
        cat("\nPrecomputed interval columns mapping (for plotted VALUE Y=", yy, ")\n", sep = "")
        ci_cols_map[[yy]] <- choose_precomputed_ci_cols(df)
      }
    }
  }
  
  # Collect saved paths for logging
  saved_separate <- list()
  saved_composite <- list()
  
  # Helper: build df_sum for a given dataframe slice + y
  build_sum <- function(df_slice, yy, group_for_sum) {
    if (use_precomputed) {
      cc <- ci_cols_map[[yy]]
      if (is.null(cc) || any(!c(cc$mean, cc$lo, cc$hi) %in% names(df_slice))) {
        stop2("Precomputed CI mapping missing or invalid for Y=", yy, ".")
      }
      ds <- df_slice
      ds$mean <- ds[[cc$mean]]
      ds$lo   <- ds[[cc$lo]]
      ds$hi   <- ds[[cc$hi]]
      ds
    } else {
      summarise_with_uncertainty(
        df_slice, x_col = x_col, y_col = yy, group_col = group_for_sum,
        mode = unc_mode, boot_R = boot_R, boot_seed = boot_seed
      )
    }
  }
  
  # CASE A: GROUP plotted together (single plot per Y; optional composite across Y)
  if (is.null(group_col) || group_plot_mode == "together") {
    
    plots <- list()
    
    for (yy in y_cols) {
      ylab_default <- infer_axis_label(df, axis = "y", fallback = yy)
      ylab <- choose_axis_label(paste0("Y-axis label for ", yy), default = ylab_default)
      
      df_sum <- build_sum(df, yy, group_col)
      
      panel_title <- if (!is.null(title_global) && length(y_cols) > 1) yy else (title_global %||% yy)
      
      if (plot_type == "timecourse") {
        p <- plot_timecourse(
          df_raw = df, df_sum = df_sum,
          x_col = x_col, y_col = yy, group_col = group_col,
          mode = unc_mode, show_points = show_points,
          title = panel_title, xlab = xlab, ylab = ylab,
          legend_pos = legend_pos,
          line_color = line_color,
          group_palette = group_palette
        )
      } else {
        as_bar <- choose_yes_no(paste0("Use bars instead of points for ", yy, "?"), default = "n")
        p <- plot_summary(
          df_raw = df, df_sum = df_sum,
          x_col = x_col, y_col = yy, group_col = group_col,
          mode = unc_mode, as_bar = as_bar,
          title = panel_title, xlab = xlab, ylab = ylab,
          legend_pos = legend_pos,
          group_palette = group_palette
        )
      }
      
      plots[[yy]] <- p
      
      # Save separate panels (always sensible)
      stem <- paste0(
        if (plot_type == "timecourse") "TimeCourse_" else "Summary_",
        safe_stem(yy), "_", ts_stamp()
      )
      paths <- save_pub(p, out_dir, stem, width = width, height = height)
      saved_separate[[paste0("Y=", yy)]] <- paths
      cat("\nSaved (", yy, "):\n  ", paths$pdf, "\n  ", paths$png, "\n", sep = "")
    }
    
    # Composite across Y (stacked) if requested and multiple Y
    if (length(y_cols) > 1 && out_mode %in% c("composite", "both")) {
      p_comp <- patchwork::wrap_plots(plots, ncol = 1)
      if (!is.null(title_global)) p_comp <- p_comp + patchwork::plot_annotation(title = title_global)
      
      stem_comp <- paste0(
        if (plot_type == "timecourse") "TimeCourse_Composite_" else "Summary_Composite_",
        ts_stamp()
      )
      comp_height <- per_panel_height * length(y_cols)
      composite_paths <- save_pub(p_comp, out_dir, stem_comp, width = width, height = comp_height)
      saved_composite[["COMPOSITE"]] <- composite_paths
      cat("\nSaved COMPOSITE:\n  ", composite_paths$pdf, "\n  ", composite_paths$png, "\n", sep = "")
      print(p_comp)
    } else {
      print(plots[[length(plots)]])
    }
    
  } else {
    # CASE B: GROUP plotted as separate files (one plot per group level)
    # We do NOT use group_col inside the plotting function (each slice is single-series).
    
    group_levels <- unique(as.character(df[[group_col]]))
    group_levels <- group_levels[order(group_levels)]
    
    cat("\nWill generate separate plots for each ", group_col, " level (", length(group_levels), "):\n", sep = "")
    cat("  ", paste(group_levels, collapse = ", "), "\n", sep = "")
    
    for (g in group_levels) {
      df_g <- df[df[[group_col]] == g, , drop = FALSE]
      
      plots_g <- list()
      
      for (yy in y_cols) {
        ylab_default <- infer_axis_label(df, axis = "y", fallback = yy)
        ylab <- choose_axis_label(paste0("Y-axis label for ", yy, " (", group_col, "=", g, ")"), default = ylab_default)
        
        df_sum_g <- build_sum(df_g, yy, group_for_sum = NULL)
        
        panel_title <- paste0(g, if (length(y_cols) > 1) paste0(" — ", yy) else "")
        
        if (plot_type == "timecourse") {
          p <- plot_timecourse(
            df_raw = df_g, df_sum = df_sum_g,
            x_col = x_col, y_col = yy, group_col = NULL,
            mode = unc_mode, show_points = show_points,
            title = panel_title, xlab = xlab, ylab = ylab,
            legend_pos = legend_pos,
            line_color = line_color,
            group_palette = "default"
          )
        } else {
          as_bar <- choose_yes_no(paste0("Use bars instead of points for ", yy, " (", g, ")?"), default = "n")
          p <- plot_summary(
            df_raw = df_g, df_sum = df_sum_g,
            x_col = x_col, y_col = yy, group_col = NULL,
            mode = unc_mode, as_bar = as_bar,
            title = panel_title, xlab = xlab, ylab = ylab,
            legend_pos = legend_pos,
            group_palette = "default"
          )
        }
        
        plots_g[[yy]] <- p
        
        # Save separate (group,Y) plots (always)
        stem <- paste0(
          if (plot_type == "timecourse") "TimeCourse_" else "Summary_",
          safe_stem(group_col), "_", safe_stem(g), "_",
          safe_stem(yy), "_", ts_stamp()
        )
        paths <- save_pub(p, out_dir, stem, width = width, height = height)
        key <- paste0(group_col, "=", g, " | Y=", yy)
        saved_separate[[key]] <- paths
        cat("\nSaved (", key, "):\n  ", paths$pdf, "\n  ", paths$png, "\n", sep = "")
      }
      
      # Composite per group level (stacked by Y) if multi-Y and requested
      if (length(y_cols) > 1 && out_mode %in% c("composite", "both")) {
        p_comp_g <- patchwork::wrap_plots(plots_g, ncol = 1)
        title_g <- if (!is.null(title_global)) paste0(title_global, " — ", group_col, "=", g) else paste0(group_col, "=", g)
        p_comp_g <- p_comp_g + patchwork::plot_annotation(title = title_g)
        
        stem_comp_g <- paste0(
          if (plot_type == "timecourse") "TimeCourse_Composite_" else "Summary_Composite_",
          safe_stem(group_col), "_", safe_stem(g), "_", ts_stamp()
        )
        comp_height <- per_panel_height * length(y_cols)
        paths_comp_g <- save_pub(p_comp_g, out_dir, stem_comp_g, width = width, height = comp_height)
        saved_composite[[paste0(group_col, "=", g)]] <- paths_comp_g
        cat("\nSaved COMPOSITE (", group_col, "=", g, "):\n  ", paths_comp_g$pdf, "\n  ", paths_comp_g$png, "\n", sep = "")
        print(p_comp_g)
      } else {
        print(plots_g[[length(plots_g)]])
      }
    }
  }
  
  # Run log
  log_path <- write_run_log(out_dir, list(
    timestamp = ts_stamp(),
    csv_path = csv_path,
    plot_type = plot_type,
    x_col = x_col,
    y_cols = paste(y_cols, collapse = ", "),
    xlab = xlab,
    group_col = group_col %||% "None",
    group_plot_mode = group_plot_mode %||% "NA",
    uncertainty_mode = unc_mode,
    auto_precomputed_ci = if (use_precomputed) as.character(auto_ci_map) else "NA",
    auto_ci_mean = if (use_precomputed && auto_ci_map) auto_ci$mean else "NA",
    auto_ci_lo = if (use_precomputed && auto_ci_map) auto_ci$lo else "NA",
    auto_ci_hi = if (use_precomputed && auto_ci_map) auto_ci$hi else "NA",
    output_mode = out_mode,
    separate_files = if (length(saved_separate)) paste(vapply(saved_separate, `[[`, "", "png"), collapse = " | ") else "NA",
    composite_files = if (length(saved_composite)) paste(vapply(saved_composite, `[[`, "", "png"), collapse = " | ") else "NA"
  ))
  cat("Run log: ", log_path, "\n", sep = "")
  
  # ---------------- HIST ----------------
  
} else if (plot_type == "hist") {
  
  y_col <- choose_column(df, "Select Y column (variable to histogram):")
  facet_col <- choose_optional_column(df, "Select FACET column (optional):")
  
  bins <- suppressWarnings(as.integer(readline("Histogram bins [30]: ")))
  if (is.na(bins) || bins < 5) bins <- 30
  
  show_ci_lines <- FALSE
  ci_mode <- 2
  boot_R <- 2000
  boot_seed <- 1
  
  if (is.null(facet_col)) {
    show_ci_lines <- choose_yes_no("Add mean + interval lines (only if not faceted)?", default = "y")
    if (show_ci_lines) {
      cat("\nInterval lines mode (for overall distribution):\n")
      cat("  [1] SE\n  [2] 95% CI (t-based)\n  [3] 95% CI (bootstrap)\n  [0] None\n")
      ci_mode <- suppressWarnings(as.integer(readline("Enter number [2]: ")))
      if (is.na(ci_mode) || !(ci_mode %in% c(0,1,2,3))) ci_mode <- 2
      if (ci_mode == 3) {
        boot_R <- suppressWarnings(as.integer(readline("Bootstrap iterations R [2000]: ")))
        if (is.na(boot_R) || boot_R < 200) boot_R <- 2000
        boot_seed <- suppressWarnings(as.integer(readline("Bootstrap seed [1]: ")))
        if (is.na(boot_seed)) boot_seed <- 1
      }
    }
  }
  
  title <- readline("Plot title (optional): ")
  if (title == "") title <- NULL
  
  xlab_default <- infer_axis_label(df, axis = "x", fallback = y_col)
  ylab_default <- infer_axis_label(df, axis = "y", fallback = "Count")
  xlab <- choose_axis_label("X-axis label", default = xlab_default)
  ylab_count <- choose_axis_label("Y-axis label", default = ylab_default)
  
  width <- suppressWarnings(as.numeric(readline("Export width inches [6.5]: ")))
  if (!is.finite(width) || width <= 0) width <- 6.5
  height <- suppressWarnings(as.numeric(readline("Export height inches [4.0]: ")))
  if (!is.finite(height) || height <= 0) height <- 4.0
  
  p <- plot_hist(
    df_raw = df, y_col = y_col, facet_col = facet_col,
    show_ci_lines = show_ci_lines, ci_mode = ci_mode,
    boot_R = boot_R, boot_seed = boot_seed,
    bins = bins, title = title, xlab = xlab,
    legend_pos = legend_pos
  ) + labs(y = ylab_count)
  
  print(p)
  stem <- paste0("Hist_", ts_stamp())
  paths <- save_pub(p, out_dir, stem, width = width, height = height)
  cat("\nSaved:\n  ", paths$pdf, "\n  ", paths$png, "\n", sep = "")
  
  log_path <- write_run_log(out_dir, list(
    timestamp = ts_stamp(),
    csv_path = csv_path,
    plot_type = plot_type,
    y_col = y_col,
    xlab = xlab,
    ylab = ylab_count,
    facet_col = facet_col %||% "None",
    bins = bins,
    ci_lines = show_ci_lines,
    ci_mode = ci_mode,
    output_pdf = paths$pdf,
    output_png = paths$png
  ))
  cat("Run log: ", log_path, "\n", sep = "")
  
  # ---------------- PCA ----------------
  
} else if (plot_type == "pca") {
  
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
  cat("\nSaved:\n  ", paths$pdf, "\n  ", paths$png, "\n", sep = "")
  
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
    output_png = paths$png
  ))
  cat("Run log: ", log_path, "\n", sep = "")
  
  # ---------------- XY ----------------
  
} else if (plot_type == "xy") {
  
  x_col <- choose_column(df, "Select X column:")
  y_col <- choose_column(df, "Select Y column:")
  group_col <- choose_optional_column(df, "Select GROUP (color) column (optional):")
  
  if (!is.numeric(df[[x_col]])) {
    suppressWarnings({ x_num <- as.numeric(df[[x_col]]) })
    if (sum(is.na(x_num)) <= 0.05 * length(x_num)) df[[x_col]] <- x_num
  }
  if (!is.numeric(df[[y_col]])) {
    suppressWarnings({ y_num <- as.numeric(df[[y_col]]) })
    if (sum(is.na(y_num)) <= 0.05 * length(y_num)) df[[y_col]] <- y_num
  }
  
  add_identity <- choose_yes_no("Add y=x reference line?", default = "y")
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
  stem <- paste0("XY_", ts_stamp())
  paths <- save_pub(p, out_dir, stem, width = width, height = height)
  cat("\nSaved:\n  ", paths$pdf, "\n  ", paths$png, "\n", sep = "")
  
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
    output_png = paths$png
  ))
  cat("Run log: ", log_path, "\n", sep = "")
  
} else {
  stop2("Unexpected plot type.")
}

cat("\nDone.\n")
