#!/usr/bin/env Rscript
# ======================================================================
# Test 1: Lagged cross-modal structure coupling (Metabolome NSV vs Transcriptome NSV)
#
# PURPOSE
#   Use existing metabolome NSV tables (S_per_day) and compute transcriptome NSV
#   (S_per_day) from a tidy long transcriptome table, then evaluate lead–lag
#   coupling via lagged rank correlations.
#
# INPUTS (defaults; override via CLI or interactive file picker)
#   --met_nsv=/path/to/NSV_RAW_full.csv
#   --tx_tidy=/path/to/Transcriptome_Tidy_Long_FIXED4columnOnly.csv
#
# Transcriptome tidy expected columns:
#   Day, SampleID, Gene, Abundance
#
# DEFINITIONS
#   - Within each day: Sample x Feature matrix, log2(x+1)
#   - Feature-feature correlation matrix R_day (Pearson, pairwise complete)
#   - RV similarity between adjacent days
#   - S_adj = 1 - RV
#   - S_per_day = S_adj / delta_day
#
# LAG CONVENTION
#   corr( metabolome(t), transcriptome(t + lag_days) )
#   => positive lag_days suggests metabolome leads transcriptome
#
# OUTPUT POLICY (MANDATORY)
#   All outputs are written under: outputs/<run_id>/
#   A Quarto QMD report is written and rendered to HTML as the final step.
# ======================================================================

# ---------------------------- Dependency handling ----------------------------
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "jsonlite", "knitr", "rmarkdown")
install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
})

# ---------------------------- Script identity capture (MANDATORY) ----------------------------
resolve_script_path <- function() {
  p <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      rstudioapi::getSourceEditorContext()$path
    } else ""
  }, error = function(e) "")
  
  if (nzchar(p) && file.exists(p)) {
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) {
      return(normalizePath(p2, winslash = "/", mustWork = FALSE))
    }
  }
  
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) {
    return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  }
  
  NA_character_
}

known_script_filename <- "03_lagged_structure_coupling.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  script_full_note <- "NOTE: script_full path detection failed; used fallback filename."
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
  script_full_note <- NULL
}

# ---------------------------- Output directory policy (MANDATORY) ----------------------------
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

analysis_name <- "LaggedStructureCoupling"
source_tag <- "MetTx"

outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

run_id <- paste0(analysis_name, "_", source_tag, "_", ts_stamp())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- Safe I/O helpers ----------------------------
stop2 <- function(...) stop(paste0(...), call. = FALSE)

safe_write_lines <- function(text, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(text, con = con, sep = "\n", useBytes = TRUE)
  invisible(TRUE)
}

safe_write_csv <- function(df, path, row.names = FALSE) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(df, con, row.names = row.names)
  invisible(TRUE)
}

pkg_versions_tbl <- function(pkgs) {
  ip <- installed.packages()
  keep <- pkgs[pkgs %in% rownames(ip)]
  data.frame(
    package = keep,
    version = unname(ip[keep, "Version"]),
    stringsAsFactors = FALSE
  )
}

write_file_inventory <- function(output_dir, inventory_path) {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  files <- files[file.info(files)$isdir %in% c(FALSE, NA)]
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(output_dir, winslash = "/", mustWork = FALSE)), "/?"),
             "", gsub("\\\\", "/", files))
  info <- file.info(files)
  inv <- data.frame(
    file = rel,
    size_bytes = as.numeric(info$size),
    mtime = as.character(info$mtime),
    stringsAsFactors = FALSE
  ) %>% arrange(file)
  safe_write_csv(inv, inventory_path, row.names = FALSE)
  inv
}

# ---------------------------- Core math (UNCHANGED) ----------------------------
rv_coeff <- function(A, B) {
  # Robert–Escoufier RV coefficient for two matrices of same dimension
  if (!all(dim(A) == dim(B))) stop2("RV: A and B must have identical dimensions.")
  
  # NA-safe sums (NA in correlation matrices is common when a feature is constant)
  num <- sum(A * B, na.rm = TRUE)
  sA  <- sum(A * A, na.rm = TRUE)
  sB  <- sum(B * B, na.rm = TRUE)
  
  den <- sqrt(sA * sB)
  
  # Guard NA/Inf/zero
  if (!is.finite(den) || den <= 0) return(NA_real_)
  
  num / den
}


compute_day_cor <- function(X) {
  stats::cor(X, use = "pairwise.complete.obs", method = "pearson")
}

compute_tx_nsv_from_tidy <- function(tidy_df) {
  req <- c("Day", "SampleID", "Gene", "Abundance")
  if (!all(req %in% names(tidy_df))) {
    stop2("Transcriptome tidy file must contain columns: ", paste(req, collapse = ", "))
  }
  
  wide <- tidy_df %>%
    mutate(Day = as.numeric(Day), Abundance = as.numeric(Abundance)) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = Abundance) %>%
    arrange(Day, SampleID)
  
  days <- sort(unique(wide$Day))
  
  cor_list <- vector("list", length(days))
  names(cor_list) <- as.character(days)
  
  for (d in days) {
    sub <- wide %>% filter(Day == d)
    X <- as.matrix(sub %>% select(-Day, -SampleID))
    X <- log2(X + 1)
    cor_list[[as.character(d)]] <- compute_day_cor(X)
  }
  
  out <- data.frame(
    Day = days, RV = NA_real_, S_adj = NA_real_, delta_day = NA_real_, S_per_day = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(days)) {
    if (i == 1) next
    d_prev <- days[i - 1]
    d_cur  <- days[i]
    delta  <- d_cur - d_prev
    
    R_prev <- cor_list[[as.character(d_prev)]]
    R_cur  <- cor_list[[as.character(d_cur)]]
    # diagnostic: how many NA entries in the correlation matrices?
    # (helps confirm the root cause if needed)
    # cat("[DIAG] Day ", d_prev, " NA in R_prev: ", sum(is.na(R_prev)),
    #     " | Day ", d_cur, " NA in R_cur: ", sum(is.na(R_cur)), "\n", sep = "")
    
    rv <- rv_coeff(R_prev, R_cur)
    s_adj <- ifelse(is.na(rv), NA_real_, 1 - rv)
    s_pd  <- ifelse(is.na(s_adj) || is.na(delta) || delta == 0, NA_real_, s_adj / delta)
    
    out$RV[i] <- rv
    out$S_adj[i] <- s_adj
    out$delta_day[i] <- delta
    out$S_per_day[i] <- s_pd
  }
  
  out
}

lagged_rank_coupling <- function(met_tbl, tx_tbl, lags_days = seq(-6, 6, by = 2),
                                 met_col = "S_per_day", tx_col = "S_per_day",
                                 method = "spearman") {
  if (!("Day" %in% names(met_tbl))) stop2("Metabolome table must contain column Day.")
  if (!("Day" %in% names(tx_tbl)))  stop2("Transcriptome table must contain column Day.")
  if (!(met_col %in% names(met_tbl))) stop2("Metabolome table missing column: ", met_col)
  if (!(tx_col %in% names(tx_tbl)))   stop2("Transcriptome table missing column: ", tx_col)
  
  met_tbl <- met_tbl %>% mutate(Day = as.numeric(Day))
  tx_tbl  <- tx_tbl  %>% mutate(Day = as.numeric(Day))
  
  res <- lapply(lags_days, function(lag) {
    joined <- met_tbl %>%
      transmute(Day = Day, met = .data[[met_col]]) %>%
      inner_join(
        tx_tbl %>% transmute(Day = Day - lag, tx = .data[[tx_col]]),
        by = "Day"
      ) %>%
      filter(is.finite(met), is.finite(tx))
    
    n <- nrow(joined)
    if (n < 4) return(data.frame(lag_days = lag, n = n, rho = NA_real_, p = NA_real_))
    
    ct <- suppressWarnings(stats::cor.test(joined$met, joined$tx, method = method, exact = FALSE))
    data.frame(lag_days = lag, n = n, rho = unname(ct$estimate), p = unname(ct$p.value))
  })
  
  bind_rows(res) %>% arrange(lag_days)
}

# ---------------------------- Inputs (defaults + CLI overrides + interactive fallback) ----------------------------
met_nsv_path <- "/mnt/data/NSV_RAW_full.csv"
tx_tidy_path <- "/mnt/data/Transcriptome_Tidy_Long_FIXED4columnOnly.csv"

ca <- commandArgs(trailingOnly = FALSE)
met_arg <- grep("^--met_nsv=", ca, value = TRUE)
tx_arg  <- grep("^--tx_tidy=", ca, value = TRUE)
if (length(met_arg) > 0) met_nsv_path <- sub("^--met_nsv=", "", met_arg[1])
if (length(tx_arg) > 0)  tx_tidy_path <- sub("^--tx_tidy=", "", tx_arg[1])

# If defaults don't exist locally, prompt interactively (or fail fast non-interactively)
if (!file.exists(met_nsv_path)) {
  if (interactive()) {
    message("Metabolome NSV path not found: ", met_nsv_path)
    message("Choose metabolome NSV CSV now...")
    met_nsv_path <- file.choose()
  } else {
    stop2("Metabolome NSV path not found: ", met_nsv_path, "\n",
          "Provide via CLI: --met_nsv=/path/to/NSV_RAW_full.csv")
  }
}
if (!file.exists(tx_tidy_path)) {
  if (interactive()) {
    message("Transcriptome tidy path not found: ", tx_tidy_path)
    message("Choose transcriptome tidy CSV now...")
    tx_tidy_path <- file.choose()
  } else {
    stop2("Transcriptome tidy path not found: ", tx_tidy_path, "\n",
          "Provide via CLI: --tx_tidy=/path/to/Transcriptome_Tidy_Long_FIXED4columnOnly.csv")
  }
}

# ---------------------------- Console summary (MANDATORY) ----------------------------
cat("\n==================== Script identity ====================\n")
cat("script_name: ", script_name, "\n", sep = "")
cat("script_path: ", script_path, "\n", sep = "")
cat("script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
if (!is.null(script_full_note)) cat(script_full_note, "\n")
cat("=========================================================\n\n")

cat("==================== Outputs ====================\n")
cat("outputs_root: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("output_dir:   ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("run_id:       ", run_id, "\n", sep = "")
cat("===============================================\n\n")

cat("==================== Inputs ====================\n")
cat("met_nsv_path: ", normalizePath(met_nsv_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("tx_tidy_path: ", normalizePath(tx_tidy_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("==============================================\n\n")

# ---------------------------- Run analysis (guarded) ----------------------------
run_ok <- FALSE
analysis_err <- NA_character_

# outputs
out_met_copy <- file.path(output_dir, "Input_Metabolome_NSV.csv")
out_tx_copy  <- file.path(output_dir, "Input_Transcriptome_Tidy.csv")
out_tx_nsv   <- file.path(output_dir, "Transcriptome_NSV_RAW_computed.csv")
out_lag_csv  <- file.path(output_dir, "Lagged_Coupling_Spearman_SperDay.csv")
out_best_txt <- file.path(output_dir, "Best_Lag_Summary.txt")
out_plot_png <- file.path(output_dir, "Lagged_Coupling_rho_vs_lag.png")

tryCatch({
  met <- read.csv(met_nsv_path, stringsAsFactors = FALSE)
  if (!all(c("Day", "S_per_day") %in% names(met))) {
    stop2("Metabolome NSV table must include at least columns: Day, S_per_day. Found: ",
          paste(names(met), collapse = ", "))
  }
  
  cat("[CHECK] Metabolome NSV columns: ", paste(names(met), collapse = ", "), "\n", sep = "")
  if (!all(c("Day", "S_per_day") %in% names(met))) {
    stop2(
      "Selected metabolome file is NOT an NSV-per-day table.\n",
      "It must contain columns Day and S_per_day.\n",
      "You selected: ", met_nsv_path, "\n",
      "Columns found: ", paste(names(met), collapse = ", ")
    )
  }
  
  
  met_use <- met %>% transmute(Day = as.numeric(Day), S_per_day = as.numeric(S_per_day))
  
  tx_tidy <- read.csv(tx_tidy_path, stringsAsFactors = FALSE)
  cat("[CHECK] Transcriptome tidy columns: ", paste(names(tx_tidy), collapse = ", "), "\n", sep = "")
  req <- c("Day", "SampleID", "Gene", "Abundance")
  if (!all(req %in% names(tx_tidy))) {
    stop2(
      "Selected transcriptome file is NOT the expected tidy long format.\n",
      "It must contain: ", paste(req, collapse = ", "), "\n",
      "You selected: ", tx_tidy_path, "\n",
      "Columns found: ", paste(names(tx_tidy), collapse = ", ")
    )
  }
  
  tx_nsv <- compute_tx_nsv_from_tidy(tx_tidy)
  
  lags <- seq(-6, 6, by = 2)
  lag_tbl <- lagged_rank_coupling(met_use, tx_nsv, lags_days = lags, met_col = "S_per_day", tx_col = "S_per_day")
  
  best <- lag_tbl %>%
    filter(is.finite(rho)) %>%
    mutate(abs_rho = abs(rho)) %>%
    arrange(desc(abs_rho), p) %>%
    slice(1)
  
  # copy inputs into output_dir for provenance
  tryCatch(file.copy(met_nsv_path, out_met_copy, overwrite = TRUE), error = function(e) NULL)
  tryCatch(file.copy(tx_tidy_path, out_tx_copy, overwrite = TRUE), error = function(e) NULL)
  
  safe_write_csv(tx_nsv, out_tx_nsv, row.names = FALSE)
  safe_write_csv(lag_tbl, out_lag_csv, row.names = FALSE)
  
  best_lines <- c(
    paste0("run_id: ", run_id),
    paste0("met_nsv_path: ", met_nsv_path),
    paste0("tx_tidy_path: ", tx_tidy_path),
    "metric: S_per_day (Spearman)",
    "lag convention: corr(metabolome(t), transcriptome(t + lag_days))",
    if (nrow(best) == 1) {
      paste0("best_lag_days: ", best$lag_days,
             " | n: ", best$n,
             " | rho: ", sprintf("%.4f", best$rho),
             " | p: ", sprintf("%.4g", best$p))
    } else {
      "best_lag_days: NA (no usable overlap / insufficient n)"
    }
  )
  safe_write_lines(best_lines, out_best_txt)
  
  plt <- ggplot(lag_tbl, aes(x = lag_days, y = rho)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    labs(
      title = "Lagged coupling: metabolome S_per_day vs transcriptome S_per_day",
      x = "lag_days (positive = metabolome leads)",
      y = "Spearman rho"
    )
  ggsave(out_plot_png, plt, width = 7.5, height = 4.5, dpi = 300)
  
  # print to console
  cat("\n==================== Test 1 results ====================\n")
  print(lag_tbl, row.names = FALSE)
  if (nrow(best) == 1) {
    cat("Best (by |rho|): lag_days = ", best$lag_days,
        " | n = ", best$n,
        " | rho = ", sprintf("%.3f", best$rho),
        " | p = ", sprintf("%.4g", best$p),
        "\n", sep = "")
  } else {
    cat("No usable lag had enough overlap to compute correlation.\n")
  }
  cat("========================================================\n\n")
  
  run_ok <- TRUE
}, error = function(e) {
  analysis_err <<- e
  run_ok <<- FALSE
})

# ---------------------------- Manifest + Report (only if run_ok) ----------------------------
manifest_path <- file.path(output_dir, "Project_Manifest.json")
inventory_path <- file.path(output_dir, "Project_Manifest_Files.csv")

deps_tbl <- pkg_versions_tbl(c(required_pkgs, "rstudioapi"))

manifest <- list(
  run_id = run_id,
  run_timestamp = as.character(Sys.time()),
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA_character_, script_full)
  ),
  input = list(
    metabolome_nsv_path = met_nsv_path,
    transcriptome_tidy_path = tx_tidy_path
  ),
  parameters = list(
    analysis_name = analysis_name,
    source_tag = source_tag,
    metric = "S_per_day",
    correlation_method = "spearman",
    lags_days = as.list(seq(-6, 6, by = 2)),
    lag_convention = "corr(metabolome(t), transcriptome(t + lag_days))",
    transcriptome_transform = "log2(x+1)",
    within_day_correlation = "pearson (pairwise complete)"
  ),
  dependencies = lapply(seq_len(nrow(deps_tbl)), function(i) {
    list(package = deps_tbl$package[i], version = deps_tbl$version[i])
  }),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  ),
  status = list(
    run_ok = run_ok,
    analysis_error = ifelse(is.null(analysis_err), NA_character_, analysis_err)
  )
)

# Write manifest + initial inventory regardless (so failures are recorded)
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
inv0 <- write_file_inventory(output_dir, inventory_path)

# Prepare header text for report (best effort)
header_text <- NULL
if (!is.na(script_full) && file.exists(script_full)) {
  lines <- readLines(script_full, warn = FALSE)
  cap <- min(length(lines), 120)
  header_text <- paste(lines[1:cap], collapse = "\n")
} else {
  header_text <- "Header text unavailable (script_full could not be resolved in this execution mode)."
}

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_name <- paste0(script_name, "_Report.html")
html_path <- file.path(output_dir, html_name)

qmd <- c(
  "---",
  paste0('title: "', script_name, ' Report"'),
  "format:",
  "  html:",
  "    toc: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "## Run metadata",
  "```{r}",
  "library(jsonlite)",
  "library(knitr)",
  sprintf("manifest <- jsonlite::fromJSON('%s', simplifyVector = FALSE)", basename(manifest_path)),
  sprintf("inv <- read.csv('%s', stringsAsFactors = FALSE)", basename(inventory_path)),
  "meta_tbl <- data.frame(",
  "  Field = c('run_id','run_timestamp','run_ok','analysis_error','script_name','script_path','script_full_path','outputs_root','output_dir','met_nsv_path','tx_tidy_path'),",
  "  Value = c(",
  "    manifest[['run_id']],",
  "    manifest[['run_timestamp']],",
  "    as.character(manifest[['status']][['run_ok']]),",
  "    ifelse(is.null(manifest[['status']][['analysis_error']]) || is.na(manifest[['status']][['analysis_error']]), 'NA', manifest[['status']][['analysis_error']]),",
  "    manifest[['script']][['name']],",
  "    manifest[['script']][['path']],",
  "    ifelse(is.null(manifest[['script']][['full_path']]) || is.na(manifest[['script']][['full_path']]), 'NA', manifest[['script']][['full_path']]),",
  "    manifest[['outputs']][['outputs_root']],",
  "    manifest[['outputs']][['output_dir']],",
  "    manifest[['input']][['metabolome_nsv_path']],",
  "    manifest[['input']][['transcriptome_tidy_path']]",
  "  ), stringsAsFactors = FALSE)",
  "knitr::kable(meta_tbl)",
  "```",
  "",
  "## Script header (as executed)",
  "```",
  header_text,
  "```",
  "",
  "## Dependencies",
  "```{r}",
  "deps <- do.call(rbind, lapply(manifest[['dependencies']], as.data.frame))",
  "knitr::kable(deps)",
  "```",
  "",
  "## Analytical logic + formulas",
  "- Compute within-day correlation matrices `R_day`.",
  "- Similarity between adjacent days uses RV.",
  "- Structural change: `S_adj(t) = 1 - RV(R_t, R_{t-1})`.",
  "- Rate: `S_per_day(t) = S_adj(t) / delta_day`.",
  "- Lagged coupling: `corr(metabolome(t), transcriptome(t + lag_days))` (Spearman).",
  "",
  "## Key plot (if available)",
  "```{r}",
  sprintf("p <- '%s'; if (file.exists(p)) knitr::include_graphics(p) else cat('Plot not available (analysis failed before plot creation).')",
          basename(out_plot_png)),
  "```",
  "",
  "## Lag table (if available)",
  "```{r}",
  sprintf("p <- '%s'; if (file.exists(p)) knitr::kable(read.csv(p)) else cat('Lag table not available (analysis failed before table creation).')",
          basename(out_lag_csv)),
  "```",
  "",
  "## Generated files",
  "```{r}",
  "inv <- read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "## Interpretation guidance",
  "- Positive best lag suggests metabolome tends to lead transcriptome by that lag; negative suggests transcriptome leads metabolome.",
  "- With small numbers of days, effect size and stability across preprocessing modes are more informative than p-values alone.",
  "",
  "## Reproducibility",
  "```{r}",
  "sessionInfo()",
  "```"
)

safe_write_lines(qmd, qmd_path)

# Render report as final step (MANDATORY) – but do not crash if analysis failed
render_ok <- TRUE
render_msg <- NULL
tryCatch({
  rmarkdown::render(
    input = qmd_path,
    output_format = "html_document",
    output_file = html_name,
    output_dir = output_dir,
    quiet = TRUE
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# Refresh inventory after render (MANDATORY)
inv1 <- write_file_inventory(output_dir, inventory_path)

# ---------------------------- Final console report ----------------------------
cat("\n==================== Wrapup ====================\n")

if (!isTRUE(run_ok)) {
  msg <- tryCatch({
    m <- conditionMessage(analysis_err)
    if (!nzchar(m)) "Unknown error (empty conditionMessage)." else m
  }, error = function(e) {
    m2 <- tryCatch(paste(analysis_err), error = function(e2) "")
    if (!nzchar(m2)) "Unknown error (analysis_err not printable)." else m2
  })
  
  cat("ANALYSIS FAILED: ", msg, "\n", sep = "")
  
} else {
  cat("ANALYSIS OK.\n")
}

if (isTRUE(render_ok) && isTRUE(file.exists(html_path))) {
  cat("HTML report created: ", normalizePath(html_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("HTML report render FAILED.\n")
  cat("QMD written to:      ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (!is.null(render_msg)) cat("Render error:        ", render_msg, "\n", sep = "")
}
cat("Manifest:            ", normalizePath(manifest_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("File inventory:      ", normalizePath(inventory_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Output directory:    ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("===============================================\n\n")
