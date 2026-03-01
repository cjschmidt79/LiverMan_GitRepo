###############################################################################
# Dispersion + Residual Dispersion by Day (biologist-friendly, SOURCE-and-run)
# Base R only (ggplot2 optional)
###############################################################################

analyze_gene_dispersion <- function(input_file = NULL,
                                    run_name = NULL,
                                    analysis_mode = c("raw", "log2p1"),
                                    pseudocount = 1,
                                    make_quarto = TRUE,
                                    make_plots = TRUE) {
  
  # ------------------------------- helpers -----------------------------------
  .safe_cat <- function(...) cat(..., "\n", sep = "")
  .now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
  
  .is_rstudio <- function() {
    if (!requireNamespace("rstudioapi", quietly = TRUE)) return(FALSE)
    isTRUE(rstudioapi::isAvailable())
  }
  
  .pick_file <- function(caption = "Select file") {
    if (.is_rstudio()) {
      p <- rstudioapi::selectFile(caption = caption)
      if (is.null(p) || !nzchar(p)) return(NA_character_)
      return(p)
    } else {
      return(file.choose())
    }
  }
  
  .readline_default <- function(prompt, default) {
    ans <- readline(prompt = paste0(prompt, " [default: ", default, "]: "))
    ans <- trimws(ans)
    if (!nzchar(ans)) default else ans
  }
  
  .select_from <- function(choices, prompt, default = 1L) {
    .safe_cat("")
    .safe_cat(prompt)
    for (i in seq_along(choices)) .safe_cat("  ", i, ") ", choices[i])
    ans <- suppressWarnings(as.integer(readline(paste0("Enter 1-", length(choices),
                                                       " (default ", default, "): "))))
    if (is.na(ans) || ans < 1 || ans > length(choices)) ans <- default
    choices[ans]
  }
  
  .sanitize_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) x <- "run"
    x
  }
  
  .safe_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
    x <- gsub("_+", "_", x)
    x
  }
  
  .write_csv <- function(df, path) {
    utils::write.csv(df, file = path, row.names = FALSE)
    path
  }
  
  .mad1 <- function(x) stats::mad(x, constant = 1, na.rm = TRUE)
  
  .bootstrap_ci <- function(x, nboot = 1000, conf = 0.95, fun = stats::median) {
    x <- x[is.finite(x)]
    if (length(x) < 3) return(c(NA_real_, NA_real_))
    n <- length(x)
    b <- replicate(nboot, fun(sample(x, size = n, replace = TRUE)))
    alpha <- (1 - conf) / 2
    stats::quantile(b, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
  }
  
  .json_escape <- function(s) {
    s <- gsub("\\\\", "\\\\\\\\", s)
    s <- gsub("\"", "\\\\\"", s)
    s <- gsub("\n", "\\\\n", s)
    s
  }
  
  .write_manifest <- function(path, lst) {
    keys <- names(lst)
    vals <- lst
    lines <- c("{")
    for (i in seq_along(keys)) {
      k <- keys[i]
      v <- vals[[i]]
      if (is.character(v)) {
        if (length(v) == 1) {
          vv <- paste0("\"", .json_escape(v), "\"")
        } else {
          vv <- paste0("[", paste0("\"", .json_escape(v), "\"", collapse = ", "), "]")
        }
      } else if (is.numeric(v) || is.logical(v)) {
        if (length(v) == 1) {
          vv <- if (is.logical(v)) tolower(as.character(v)) else as.character(v)
        } else {
          vv <- paste0("[", paste0(as.character(v), collapse = ", "), "]")
        }
      } else {
        vv <- paste0("\"", .json_escape(paste(v, collapse = "; ")), "\"")
      }
      comma <- if (i < length(keys)) "," else ""
      lines <- c(lines, paste0("  \"", .json_escape(k), "\": ", vv, comma))
    }
    lines <- c(lines, "}")
    writeLines(lines, con = path)
    path
  }
  
  # --------- robust script identity capture (with interactive fallback) -------
  .get_script_identity_strict <- function() {
    script_full <- NA_character_
    
    # 1) Best: RStudio source editor context
    if (.is_rstudio()) {
      try({
        ctx <- rstudioapi::getSourceEditorContext()
        if (!is.null(ctx$path) && nzchar(ctx$path)) script_full <- ctx$path
      }, silent = TRUE)
    }
    
    # 2) If run via Rscript --file
    if (is.na(script_full)) {
      try({
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        idx <- grep("^--file=", cmdArgs)
        if (length(idx) > 0) script_full <- sub("^--file=", "", cmdArgs[idx[1]])
      }, silent = TRUE)
    }
    
    # 3) If still NA and interactive: prompt user to pick the script file ONCE
    if (is.na(script_full) && interactive()) {
      .safe_cat("")
      .safe_cat("Could not automatically detect script path.")
      .safe_cat("Please select THIS script file (Dispersion_ResidualDispersion.R) so it can be recorded in the manifest.")
      picked <- .pick_file("Select the script file currently being sourced")
      if (is.character(picked) && nzchar(picked) && file.exists(picked)) {
        script_full <- picked
      }
    }
    
    if (!is.na(script_full)) {
      script_full <- normalizePath(script_full, winslash = "/", mustWork = FALSE)
      return(list(
        script_full = script_full,
        script_name = basename(script_full),
        script_path = dirname(script_full)
      ))
    } else {
      return(list(
        script_full = NA_character_,
        script_name = NA_character_,
        script_path = NA_character_
      ))
    }
  }
  
  # ----------------------------- inputs / UI --------------------------------
  timestamp <- .now_stamp()
  sid <- .get_script_identity_strict()
  
  .safe_cat("")
  .safe_cat("Script identity:")
  .safe_cat("  script_full: ", ifelse(is.na(sid$script_full), "NA", sid$script_full))
  .safe_cat("  script_name: ", ifelse(is.na(sid$script_name), "NA", sid$script_name))
  .safe_cat("  script_path: ", ifelse(is.na(sid$script_path), "NA", sid$script_path))
  
  if (is.null(input_file) || !nzchar(input_file)) {
    .safe_cat("")
    .safe_cat("Select the input CSV file (tidy long table).")
    input_file <- .pick_file("Select input CSV (tidy long table)")
  }
  if (is.null(input_file) || !nzchar(input_file) || !file.exists(input_file)) {
    stop("Input file not found. Aborting.")
  }
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  
  if (is.null(run_name) || !nzchar(run_name)) {
    run_name <- .readline_default("Enter a short run name", "GeneDispersion")
  }
  run_name <- .sanitize_name(run_name)
  
  analysis_mode <- match.arg(analysis_mode)
  
  if (interactive()) {
    analysis_mode <- .select_from(c("raw", "log2p1"),
                                  "Choose analysis scale:",
                                  default = ifelse(analysis_mode == "raw", 1L, 2L))
    if (analysis_mode == "log2p1") {
      pc_in <- .readline_default("Enter pseudocount PC for log2(TPM+PC)", as.character(pseudocount))
      pseudocount <- suppressWarnings(as.numeric(pc_in))
      if (!is.finite(pseudocount) || pseudocount < 0) {
        pseudocount <- 1
        .safe_cat("Invalid PC; using 1.")
      }
    }
  }
  
  # --------------------------- output directory ------------------------------
  proj_root <- getwd()
  outputs_root <- file.path(proj_root, "outputs")
  if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)
  
  out_dir <- file.path(outputs_root, paste0(run_name, "_", timestamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  
  fig_dir <- file.path(out_dir, "Figures"); dir.create(fig_dir, showWarnings = FALSE)
  tab_dir <- file.path(out_dir, "Tables");  dir.create(tab_dir, showWarnings = FALSE)
  
  .safe_cat("")
  .safe_cat("Output directory:")
  .safe_cat(out_dir)
  
  # ------------------------------- load data --------------------------------
  dat <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(dat) < 4) stop("Input must have at least 4 columns (tidy long).")
  
  cn <- colnames(dat)
  
  day_col    <- .select_from(cn, "Select the DAY column (numeric/integer-like):", default = 1L)
  sample_col <- .select_from(cn, "Select the SAMPLE ID column (replicate identifier):", default = 2L)
  gene_col   <- .select_from(cn, "Select the GENE SYMBOL/ID column:", default = 3L)
  value_col  <- .select_from(cn, "Select the ABUNDANCE/EXPRESSION column (numeric):", default = 4L)
  
  Day    <- dat[[day_col]]
  Sample <- dat[[sample_col]]
  Gene   <- dat[[gene_col]]
  Xraw   <- dat[[value_col]]
  
  Day_num <- suppressWarnings(as.integer(as.character(Day)))
  if (anyNA(Day_num)) {
    .safe_cat("")
    .safe_cat("WARNING: Day column contains non-integers. Attempting Day from SampleID as DAY_REP (before first underscore).")
    Day_str <- as.character(Sample)
    Day_num2 <- suppressWarnings(as.integer(sub("_.*$", "", Day_str)))
    if (anyNA(Day_num2)) stop("Could not parse Day. Provide a valid Day column or DAY_REP sample naming.")
    Day_num <- Day_num2
  }
  
  Sample_chr <- as.character(Sample)
  Gene_chr   <- as.character(Gene)
  
  Xraw_num <- suppressWarnings(as.numeric(as.character(Xraw)))
  if (any(!is.finite(Xraw_num))) {
    bad <- sum(!is.finite(Xraw_num))
    stop(paste0("Abundance column has ", bad, " non-numeric/non-finite values. Fix input."))
  }
  
  if (analysis_mode == "raw") {
    X <- Xraw_num
    scale_label <- "raw"
    transform_label <- "none"
  } else {
    X <- log2(Xraw_num + pseudocount)
    scale_label <- paste0("log2p", pseudocount)
    transform_label <- paste0("log2(TPM+PC), PC=", pseudocount)
  }
  
  df <- data.frame(Day = Day_num, SampleID = Sample_chr, Gene = Gene_chr, Value = X,
                   stringsAsFactors = FALSE)
  
  days <- sort(unique(df$Day))
  n_genes <- length(unique(df$Gene))
  n_samples_by_day <- tapply(df$SampleID, df$Day, function(z) length(unique(z)))
  
  .safe_cat("")
  .safe_cat("Days: ", paste(days, collapse = ", "))
  .safe_cat("Unique genes: ", n_genes)
  .safe_cat("Unique samples per day:")
  for (d in names(n_samples_by_day)) .safe_cat("  Day ", d, ": ", n_samples_by_day[[d]])
  
  # ------------------- A) within-gene dispersion per day ---------------------
  key <- paste(df$Day, df$Gene, sep = "||")
  split_idx <- split(seq_len(nrow(df)), key)
  
  per_gene_day <- lapply(names(split_idx), function(k) {
    ii <- split_idx[[k]]
    dd <- df$Day[ii[1]]
    gg <- df$Gene[ii[1]]
    x  <- df$Value[ii]
    
    mu <- mean(x)
    va <- if (length(x) >= 2) stats::var(x) else NA_real_
    sd <- sqrt(va)
    iq <- stats::IQR(x, na.rm = TRUE)
    md <- .mad1(x)
    cv <- if (is.finite(mu) && abs(mu) > 1e-12) sd / abs(mu) else NA_real_
    
    data.frame(Day = dd, Gene = gg,
               n_reps = length(x),
               mean = mu, var = va, sd = sd,
               iqr = iq, mad = md, cv = cv,
               stringsAsFactors = FALSE)
  })
  per_gene_day <- do.call(rbind, per_gene_day)
  
  summarize_day <- function(metric) {
    out <- lapply(days, function(d) {
      v <- per_gene_day[per_gene_day$Day == d, metric]
      v <- v[is.finite(v)]
      med <- if (length(v) > 0) stats::median(v) else NA_real_
      ci <- .bootstrap_ci(v, nboot = 1000, conf = 0.95, fun = stats::median)
      data.frame(Day = d,
                 metric = metric,
                 n_genes = length(unique(per_gene_day$Gene[per_gene_day$Day == d])),
                 median = med,
                 ci_low = ci[1], ci_high = ci[2],
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }
  
  day_disp_all <- rbind(
    summarize_day("sd"),
    summarize_day("iqr"),
    summarize_day("mad"),
    summarize_day("cv")
  )
  day_disp_all$analysis_scale <- scale_label
  day_disp_all$transform <- transform_label
  
  # ------------------ B) residual dispersion (mean–variance) -----------------
  eps_var <- 1e-12
  per_gene_day$resid_log_var <- NA_real_
  per_gene_day$fitted_log_var <- NA_real_
  
  for (d in days) {
    sub <- per_gene_day[per_gene_day$Day == d, ]
    ok <- is.finite(sub$mean) & is.finite(sub$var) & sub$var >= 0
    sub_ok <- sub[ok, , drop = FALSE]
    if (nrow(sub_ok) < 10) next
    
    y <- log(sub_ok$var + eps_var)
    
    if (analysis_mode == "raw") {
      x <- log(sub_ok$mean + max(pseudocount, 1e-6))
    } else {
      x <- sub_ok$mean
    }
    
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    genes_keep <- sub_ok$Gene[keep]
    if (length(x) < 10) next
    
    fit <- try(stats::loess(y ~ x, span = 0.75, degree = 2), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    yhat <- stats::predict(fit, newdata = data.frame(x = x))
    resid <- y - yhat
    
    for (i in seq_along(genes_keep)) {
      idx <- which(per_gene_day$Day == d & per_gene_day$Gene == genes_keep[i])
      if (length(idx) == 1) {
        per_gene_day$fitted_log_var[idx] <- yhat[i]
        per_gene_day$resid_log_var[idx]  <- resid[i]
      }
    }
  }
  
  day_resid <- lapply(days, function(d) {
    v <- per_gene_day$resid_log_var[per_gene_day$Day == d]
    v <- v[is.finite(v)]
    med <- if (length(v) > 0) stats::median(v) else NA_real_
    ci <- .bootstrap_ci(v, nboot = 1000, conf = 0.95, fun = stats::median)
    data.frame(Day = d,
               metric = "resid_log_var",
               n_genes = length(unique(per_gene_day$Gene[per_gene_day$Day == d])),
               median = med,
               ci_low = ci[1], ci_high = ci[2],
               analysis_scale = scale_label,
               transform = transform_label,
               stringsAsFactors = FALSE)
  })
  day_resid <- do.call(rbind, day_resid)
  
  day_resid_mad <- lapply(days, function(d) {
    v <- per_gene_day$resid_log_var[per_gene_day$Day == d]
    v <- v[is.finite(v)]
    data.frame(Day = d,
               metric = "mad_resid_log_var",
               n_genes = length(unique(per_gene_day$Gene[per_gene_day$Day == d])),
               median = .mad1(v),
               ci_low = NA_real_, ci_high = NA_real_,
               analysis_scale = scale_label,
               transform = transform_label,
               stringsAsFactors = FALSE)
  })
  day_resid_mad <- do.call(rbind, day_resid_mad)
  
  # ------------------------------- write tables ------------------------------
  out_per_gene_day <- .write_csv(per_gene_day,
                                 file.path(tab_dir, .safe_filename(paste0("PerGene_ByDay_Stats_", scale_label, ".csv"))))
  
  out_day_disp <- .write_csv(day_disp_all,
                             file.path(tab_dir, .safe_filename(paste0("PerDay_Dispersion_Summary_", scale_label, ".csv"))))
  
  out_day_resid <- .write_csv(rbind(day_resid, day_resid_mad),
                              file.path(tab_dir, .safe_filename(paste0("PerDay_ResidualDispersion_Summary_", scale_label, ".csv"))))
  
  # Wide summary
  wide_day <- data.frame(Day = days, stringsAsFactors = FALSE)
  
  .merge_metric <- function(metric_name) {
    tmp <- day_disp_all[day_disp_all$metric == metric_name, c("Day", "median")]
    colnames(tmp)[2] <- paste0("median_", metric_name)
    tmp
  }
  
  wide_day <- merge(wide_day, .merge_metric("sd"),  by = "Day", all.x = TRUE)
  wide_day <- merge(wide_day, .merge_metric("iqr"), by = "Day", all.x = TRUE)
  wide_day <- merge(wide_day, .merge_metric("mad"), by = "Day", all.x = TRUE)
  wide_day <- merge(wide_day, .merge_metric("cv"),  by = "Day", all.x = TRUE)
  
  tmp_res <- day_resid[, c("Day", "median")]
  colnames(tmp_res)[2] <- "median_resid_log_var"
  wide_day <- merge(wide_day, tmp_res, by = "Day", all.x = TRUE)
  
  wide_day$analysis_scale <- scale_label
  wide_day$transform <- transform_label
  
  out_wide_day <- .write_csv(wide_day,
                             file.path(tab_dir, .safe_filename(paste0("PerDay_WideSummary_", scale_label, ".csv"))))
  
  # ------------------------------- plots ------------------------------------
  plot_paths <- list()
  
  if (isTRUE(make_plots)) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      dd_plot <- day_disp_all[day_disp_all$metric %in% c("sd", "iqr", "mad"), ]
      
      p1 <- ggplot2::ggplot(dd_plot, ggplot2::aes(x = Day, y = median, group = metric)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, ymax = ci_high, fill = metric),
                             alpha = 0.15, colour = NA) +
        ggplot2::labs(title = paste0("Within-gene dispersion across replicates (", scale_label, ")"),
                      subtitle = "Medians across genes; 95% bootstrap CI",
                      x = "Day", y = "Dispersion (median across genes)") +
        ggplot2::theme_bw()
      
      f1 <- file.path(fig_dir, .safe_filename(paste0("Dispersion_ByDay_", scale_label, ".png")))
      ggplot2::ggsave(filename = f1, plot = p1, width = 7, height = 4.5, dpi = 600)
      plot_paths$dispersion <- f1
      
      p2 <- ggplot2::ggplot(day_resid, ggplot2::aes(x = Day, y = median)) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, ymax = ci_high),
                             alpha = 0.15, colour = NA) +
        ggplot2::labs(title = paste0("Residual dispersion (mean–variance corrected) (", scale_label, ")"),
                      subtitle = "Residual = log(var) - fitted trend; median across genes; 95% bootstrap CI",
                      x = "Day", y = "Median residual log-variance") +
        ggplot2::theme_bw()
      
      f2 <- file.path(fig_dir, .safe_filename(paste0("ResidualDispersion_ByDay_", scale_label, ".png")))
      ggplot2::ggsave(filename = f2, plot = p2, width = 7, height = 4.5, dpi = 600)
      plot_paths$residual_dispersion <- f2
    }
  }
  
  # ------------------------------ provenance --------------------------------
  manifest <- list(
    timestamp = timestamp,
    run_name = run_name,
    output_dir = out_dir,
    input_file = input_file,
    analysis_mode = analysis_mode,
    pseudocount = pseudocount,
    transform = transform_label,
    day_column = day_col,
    sample_column = sample_col,
    gene_column = gene_col,
    value_column = value_col,
    script_name = sid$script_name,
    script_path = sid$script_path,
    script_full = sid$script_full,
    n_rows = nrow(df),
    n_genes = n_genes,
    days = as.character(days)
  )
  
  manifest_path <- file.path(out_dir, "manifest.json")
  .write_manifest(manifest_path, manifest)
  
  inventory_path <- file.path(out_dir, "inventory_files.txt")
  inv <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  writeLines(inv, con = inventory_path)
  
  # ------------------------------ Quarto report ------------------------------
  qmd_path <- file.path(out_dir, .safe_filename(paste0("Dispersion_Report_", scale_label, ".qmd")))
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  
  if (isTRUE(make_quarto) && requireNamespace("quarto", quietly = TRUE)) {
    qmd <- c(
      "---",
      paste0("title: \"Gene dispersion and residual dispersion (", scale_label, ")\""),
      paste0("date: \"", format(Sys.Date()), "\""),
      "format: html",
      "---",
      "",
      "## Summary",
      "",
      "Within-gene dispersion across biological replicates per day, and mean–variance-corrected residual dispersion.",
      "",
      "## Metadata",
      "",
      "```{r}",
      "params <- list(",
      paste0("  script_full = \"", .json_escape(ifelse(is.na(sid$script_full), "", sid$script_full)), "\","),
      paste0("  script_name = \"", .json_escape(ifelse(is.na(sid$script_name), "", sid$script_name)), "\","),
      paste0("  script_path = \"", .json_escape(ifelse(is.na(sid$script_path), "", sid$script_path)), "\","),
      paste0("  input_file = \"", .json_escape(input_file), "\","),
      paste0("  output_dir = \"", .json_escape(out_dir), "\","),
      paste0("  analysis_mode = \"", .json_escape(analysis_mode), "\","),
      paste0("  pseudocount = ", pseudocount),
      ")",
      "print(params)",
      "```",
      "",
      "## Key outputs",
      "",
      "```{r}",
      "print(list.files(\"Tables\", full.names = TRUE))",
      "```",
      "",
      "```{r, results='asis', echo=FALSE}",
      "imgs <- list.files(\"Figures\", pattern = \"\\\\.png$\", full.names = TRUE)",
      "for (p in imgs) cat(\"![](\", p, \")\\n\\n\", sep = \"\")",
      "```",
      "",
      "## Session info",
      "",
      "```{r}",
      "sessionInfo()",
      "```"
    )
    
    writeLines(qmd, con = qmd_path)
    
    old_wd <- getwd()
    setwd(out_dir)
    on.exit(setwd(old_wd), add = TRUE)
    
    try(quarto::quarto_render(basename(qmd_path)), silent = TRUE)
  }
  
  # ------------------------------- return ------------------------------------
  result <- list(
    out_dir = out_dir,
    tables = list(
      per_gene_by_day = out_per_gene_day,
      per_day_dispersion = out_day_disp,
      per_day_residual = out_day_resid,
      per_day_wide = out_wide_day
    ),
    figures = plot_paths,
    manifest = manifest_path,
    inventory = inventory_path,
    quarto = list(
      attempted = isTRUE(make_quarto),
      qmd = if (file.exists(qmd_path)) qmd_path else NA_character_,
      html = if (file.exists(html_path)) html_path else NA_character_
    ),
    params = manifest
  )
  
  if (interactive()) {
    .safe_cat("")
    .safe_cat("Done. Results written to:")
    .safe_cat(out_dir)
    if (file.exists(out_wide_day)) {
      try(utils::View(utils::read.csv(out_wide_day, stringsAsFactors = FALSE)), silent = TRUE)
    }
  }
  
  result
}

###############################################################################
# AUTO-RUN WHEN SOURCED (RStudio / interactive use)
###############################################################################
if (interactive()) {
  message("Running analyze_gene_dispersion() interactively...")
  res <- analyze_gene_dispersion()
  message("Done. Results written to:")
  message(res$out_dir)
  if (!is.null(res$tables$per_day_wide) && file.exists(res$tables$per_day_wide)) {
    try(View(read.csv(res$tables$per_day_wide, stringsAsFactors = FALSE)), silent = TRUE)
  }
}
