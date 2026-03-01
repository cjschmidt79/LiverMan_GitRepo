#!/usr/bin/env Rscript
# ======================================================================
# EffDim Bootstrap Analysis Pipeline
#
# PURPOSE
#   Consume Bootstrap_EffDim_ByDay.csv (long-form: b, Day, eff_dim) and compute:
#     A) Per-day summaries: median, IQR, 95% CI, SD, N, skewness (optional)
#     B) Between-day contrasts:
#         - P(A > B) dominance probability
#         - median difference (median_A - median_B)
#         - optional density overlap coefficient
#     C) Phase-level summaries (user-defined phases):
#         - phase median, IQR, 95% CI
#         - phase dominance P(phase1 > phase2) + median diff
#
# INPUT
#   A CSV with at least columns: Day, eff_dim
#   Optional: b (bootstrap id). If missing, script will create b within each Day.
# eff dim table 1 from  Change_Point_Directionv4ErrorBarsV4wthHTML.R
# OUTPUT POLICY (MANDATORY)
#   All outputs go under: outputs/<run_id>/
#   plus a log file and (optional) HTML report.
#
# NO HARDCODING
#   - You will be prompted for:
#       * input file
#       * output labels (analysis_name, source label)
#       * per-day statistics toggles (skewness, overlap)
#       * which day-pair contrasts to compute (custom list OR "all adjacent")
#       * phases (names + day sets)
#       * optional HTML report
# ======================================================================

suppressPackageStartupMessages({
  pkgs <- c("data.table","jsonlite","knitr","rmarkdown","ggplot2")
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org", quiet=TRUE)
  library(data.table)
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
  library(ggplot2)
})

# -----------------------------
# Helpers (prompts + logging)
# -----------------------------
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

is_nonempty_string <- function(x) is.character(x) && length(x)==1L && !is.na(x) && nzchar(x)

prompt_string <- function(msg, default="") {
  cat(msg, " [default=", default, "]: ", sep="")
  x <- trimws(readline())
  if (!is_nonempty_string(x)) return(default)
  x
}

prompt_yesno <- function(msg, default_yes=TRUE) {
  def <- if (default_yes) "Y" else "N"
  cat(msg, " (Y/N) [default=", def, "]: ", sep="")
  x <- toupper(trimws(readline()))
  if (!is_nonempty_string(x)) return(default_yes)
  x %in% c("Y","YES")
}

prompt_int <- function(msg, default) {
  cat(msg, " [default=", default, "]: ", sep="")
  x <- trimws(readline())
  if (!is_nonempty_string(x)) return(as.integer(default))
  v <- suppressWarnings(as.integer(x))
  if (!is.finite(v)) return(as.integer(default))
  v
}

stop2 <- function(...) stop(paste0(...), call.=FALSE)

make_logger <- function(path) {
  force(path)
  function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- paste0("[", stamp, "] ", msg)
    cat(line, "\n", sep="")
    cat(line, "\n", sep="", file=path, append=TRUE)
  }
}

# -----------------------------
# Output policy
# -----------------------------
outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive=TRUE, showWarnings=FALSE)

# registries for manifest
.input_files  <- character()
.output_files <- character()

register_input <- function(path) {
  .input_files <<- unique(c(.input_files, normalizePath(path, winslash="/", mustWork=FALSE)))
  invisible(path)
}
register_output <- function(path) {
  .output_files <<- unique(c(.output_files, normalizePath(path, winslash="/", mustWork=FALSE)))
  invisible(path)
}

deps_table <- function(pkgs) {
  ok <- pkgs[pkgs %in% rownames(installed.packages())]
  data.table(
    package = ok,
    version = vapply(ok, function(p) as.character(utils::packageVersion(p)), character(1))
  )
}

write_manifest_json <- function(manifest, path) {
  jsonlite::write_json(manifest, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(path)
}

write_files_inventory <- function(output_dir, inventory_csv_path) {
  files <- list.files(output_dir, recursive=TRUE, full.names=TRUE, all.files=FALSE, include.dirs=FALSE)
  files <- files[file.exists(files)]
  rel <- sub(paste0("^", gsub("\\\\","/", normalizePath(output_dir, winslash="/", mustWork=FALSE)), "/?"),
             "", gsub("\\\\","/", files))
  info <- file.info(files)
  dt <- data.table(
    relative_path = rel,
    full_path = gsub("\\\\","/", files),
    bytes = as.numeric(info$size),
    modified = as.character(info$mtime)
  )
  data.table::fwrite(dt, inventory_csv_path)
  invisible(inventory_csv_path)
}

# -----------------------------
# Core stats
# -----------------------------
skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x); s <- sd(x)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  mean((x - m)^3) / (s^3)
}

# Overlap coefficient between two distributions using KDE on a common grid
overlap_coef <- function(x, y, n=2048L) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 10 || length(y) < 10) return(NA_real_)
  rx <- range(x); ry <- range(y)
  lo <- min(rx[1], ry[1]); hi <- max(rx[2], ry[2])
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(NA_real_)
  dx <- density(x, from=lo, to=hi, n=n)
  dy <- density(y, from=lo, to=hi, n=n)
  # densities already on same grid (from/to/n)
  sum(pmin(dx$y, dy$y)) * (dx$x[2] - dx$x[1])
}

# Dominance probability: P(x > y)
p_dominance <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  # Efficient exact computation for bootstrap-aligned reps:
  # if b is aligned across days, compare by b; otherwise fallback to Monte Carlo pairing.
  # We'll do aligned if lengths equal; else MC.
  if (length(x) == length(y)) return(mean(x > y))
  # MC pairing
  B <- min(20000L, length(x) * length(y))
  ix <- sample.int(length(x), B, replace=TRUE)
  iy <- sample.int(length(y), B, replace=TRUE)
  mean(x[ix] > y[iy])
}

# -----------------------------
# User inputs (NO HARDCODING)
# -----------------------------
cat("\n=== EffDim Bootstrap Analysis Pipeline ===\n")
cat("This script expects a CSV with columns including Day and eff_dim.\n\n")

analysis_name <- prompt_string("analysis_name (used in run_id)", "EffDimBootstrapAnalysis")
analysis_name <- gsub("[^A-Za-z0-9._-]+", "_", analysis_name)

SOURCE <- prompt_string("source label (used in run_id; e.g., liver_tx)", "UnknownSource")
SOURCE <- gsub("[^A-Za-z0-9._-]+", "_", SOURCE)

# Input file
infile <- NA_character_
if (interactive()) {
  cat("\nSelect the INPUT CSV (Bootstrap_EffDim_ByDay.csv; long-form recommended)\n")
  infile <- file.choose()
} else {
  stop2("Non-interactive mode: please run in R interactive OR adapt script to pass infile path.")
}
register_input(infile)

# Toggles
do_skew    <- prompt_yesno("Compute skewness per day?", TRUE)
do_overlap <- prompt_yesno("Compute overlap coefficient for day contrasts (KDE)?", FALSE)
do_html    <- prompt_yesno("Write QMD + render HTML report?", TRUE)

# Contrast selection
contrast_mode <- prompt_string(
  "Day-contrast mode: 'adjacent' (all adjacent pairs) OR 'custom' (you type pairs)",
  "adjacent"
)
contrast_mode <- tolower(trimws(contrast_mode))
if (!contrast_mode %in% c("adjacent","custom")) contrast_mode <- "adjacent"

custom_pairs_text <- ""
if (contrast_mode == "custom") {
  cat("\nEnter day pairs as 'd1-d2' separated by commas.\n")
  cat("Example: 4-10, 10-14, 14-16, 18-20\n")
  custom_pairs_text <- prompt_string("Day pairs", "")
}

# Phase definitions
do_phases <- prompt_yesno("Compute phase-level summaries (you define phases)?", TRUE)
phase_text <- ""
if (do_phases) {
  cat("\nDefine phases as NAME:day,day,day; NAME2:day,day\n")
  cat("Example: Early:4,8,10; Transition:12,14; Late:16,18,20\n")
  phase_text <- prompt_string("Phases", "Early:4,8,10; Transition:12,14; Late:16,18,20")
}

# -----------------------------
# Create run directory, logger, manifest paths
# -----------------------------
run_id <- paste0(analysis_name, "_", SOURCE, "_", ts_stamp())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_dir, "Plots"), recursive=TRUE, showWarnings=FALSE)

logfile <- file.path(output_dir, paste0("LOG_", ts_stamp(), ".txt"))
register_output(logfile)
logf <- make_logger(logfile)

manifest_json_path <- file.path(output_dir, "Project_Manifest.json")
manifest_files_csv <- file.path(output_dir, "Project_Manifest_Files.csv")
register_output(manifest_json_path)
register_output(manifest_files_csv)

logf("START")
logf(paste0("run_id: ", run_id))
logf(paste0("Input file: ", infile))

# -----------------------------
# Read + validate
# -----------------------------
DT <- data.table::fread(infile)
needed <- c("Day","eff_dim")
miss <- setdiff(needed, names(DT))
if (length(miss)) stop2("Missing required columns: ", paste(miss, collapse=", "))

DT[, Day := suppressWarnings(as.integer(as.character(Day)))]
DT[, eff_dim := suppressWarnings(as.numeric(eff_dim))]
DT <- DT[is.finite(Day) & is.finite(eff_dim)]

if (!"b" %in% names(DT)) {
  logf("Column 'b' missing: constructing bootstrap id within each Day.")
  DT[, b := seq_len(.N), by=Day]
} else {
  DT[, b := suppressWarnings(as.integer(b))]
  if (any(!is.finite(DT$b))) {
    logf("WARNING: non-finite b values detected; reconstructing b within each Day.")
    DT[, b := seq_len(.N), by=Day]
  }
}

days_vec <- sort(unique(DT$Day))
if (length(days_vec) < 2) stop2("Need at least 2 days in data.")
logf(paste0("Days found: ", paste(days_vec, collapse=", ")))

# -----------------------------
# A) Per-day summaries
# -----------------------------
logf("Computing per-day summaries...")
day_summ <- DT[, .(
  n = .N,
  median_eff_dim = median(eff_dim),
  mean_eff_dim   = mean(eff_dim),
  sd_eff_dim     = sd(eff_dim),
  iqr_eff_dim    = IQR(eff_dim),
  ci_lo_2p5      = as.numeric(quantile(eff_dim, 0.025)),
  ci_hi_97p5     = as.numeric(quantile(eff_dim, 0.975)),
  skew_eff_dim   = if (isTRUE(do_skew)) skewness(eff_dim) else NA_real_
), by=Day][order(Day)]

out_day <- file.path(output_dir, "EffDim_ByDay_Summary.csv")
data.table::fwrite(day_summ, out_day)
register_output(out_day)
logf("Wrote EffDim_ByDay_Summary.csv")

# -----------------------------
# B) Between-day contrasts
# -----------------------------
logf("Computing between-day contrasts...")

parse_pairs <- function(text) {
  text <- gsub("\\s+", "", text)
  if (!nzchar(text)) return(data.table(d1=integer(), d2=integer()))
  parts <- unlist(strsplit(text, ",", fixed=TRUE))
  out <- lapply(parts, function(p) {
    qq <- unlist(strsplit(p, "-", fixed=TRUE))
    if (length(qq) != 2) return(NULL)
    d1 <- suppressWarnings(as.integer(qq[1]))
    d2 <- suppressWarnings(as.integer(qq[2]))
    if (!is.finite(d1) || !is.finite(d2)) return(NULL)
    data.table(d1=d1, d2=d2)
  })
  rbindlist(out, fill=TRUE)
}

pairs_dt <- NULL
if (contrast_mode == "adjacent") {
  pairs_dt <- data.table(
    d1 = days_vec[-length(days_vec)],
    d2 = days_vec[-1]
  )
} else {
  pairs_dt <- parse_pairs(custom_pairs_text)
}

# keep only pairs present
pairs_dt <- unique(pairs_dt[is.finite(d1) & is.finite(d2)])
pairs_dt <- pairs_dt[d1 %in% days_vec & d2 %in% days_vec]
if (nrow(pairs_dt) == 0) {
  logf("WARNING: no valid day pairs for contrasts; skipping contrast table.")
  contrasts <- data.table()
} else {
  get_x <- function(d) DT[Day==d, eff_dim]
  contrasts <- pairs_dt[, {
    x <- get_x(d1)
    y <- get_x(d2)
    list(
      n1 = length(x), n2 = length(y),
      median1 = median(x), median2 = median(y),
      median_diff = median(x) - median(y),
      P_d1_gt_d2  = p_dominance(x, y),
      overlap     = if (isTRUE(do_overlap)) overlap_coef(x, y) else NA_real_
    )
  }, by=.(d1, d2)]
  
  # Add an interpretability column (direction label)
  contrasts[, direction := fifelse(median_diff > 0, "d1 higher-dim", fifelse(median_diff < 0, "d2 higher-dim", "equal"))]
}

out_con <- file.path(output_dir, "EffDim_DayContrasts.csv")
data.table::fwrite(contrasts, out_con)
register_output(out_con)
logf("Wrote EffDim_DayContrasts.csv")

# -----------------------------
# C) Phase-level summaries
# -----------------------------
phase_summ <- data.table()
phase_contrasts <- data.table()

parse_phases <- function(text) {
  # NAME:4,8,10; NAME2:12,14
  parts <- unlist(strsplit(text, ";", fixed=TRUE))
  out <- lapply(parts, function(p) {
    p <- trimws(p)
    if (!nzchar(p)) return(NULL)
    kv <- unlist(strsplit(p, ":", fixed=TRUE))
    if (length(kv) != 2) return(NULL)
    nm <- trimws(kv[1])
    ds <- unlist(strsplit(gsub("\\s+","", kv[2]), ",", fixed=TRUE))
    dv <- suppressWarnings(as.integer(ds))
    dv <- dv[is.finite(dv)]
    if (!nzchar(nm) || length(dv) == 0) return(NULL)
    data.table(phase=nm, Day=dv)
  })
  rbindlist(out, fill=TRUE)
}

if (isTRUE(do_phases)) {
  logf("Computing phase summaries...")
  ph_map <- parse_phases(phase_text)
  if (nrow(ph_map) == 0) {
    logf("WARNING: could not parse phases; skipping phase summaries.")
  } else {
    # Keep only days present in data
    ph_map <- unique(ph_map[Day %in% days_vec])
    if (nrow(ph_map) == 0) {
      logf("WARNING: phase days not found in dataset; skipping phase summaries.")
    } else {
      DTp <- DT[ph_map, on="Day", nomatch=0][, .(phase, Day, b, eff_dim)]
      phase_summ <- DTp[, .(
        n = .N,
        median_eff_dim = median(eff_dim),
        mean_eff_dim   = mean(eff_dim),
        sd_eff_dim     = sd(eff_dim),
        iqr_eff_dim    = IQR(eff_dim),
        ci_lo_2p5      = as.numeric(quantile(eff_dim, 0.025)),
        ci_hi_97p5     = as.numeric(quantile(eff_dim, 0.975))
      ), by=phase][order(phase)]
      
      out_ph <- file.path(output_dir, "EffDim_Phase_Summary.csv")
      data.table::fwrite(phase_summ, out_ph)
      register_output(out_ph)
      logf("Wrote EffDim_Phase_Summary.csv")
      
      # Phase vs phase contrasts (all pairs of phases)
      phs <- sort(unique(DTp$phase))
      if (length(phs) >= 2) {
        ph_pairs <- CJ(p1=phs, p2=phs)[p1 < p2]
        phase_contrasts <- ph_pairs[, {
          x <- DTp[phase==p1, eff_dim]
          y <- DTp[phase==p2, eff_dim]
          list(
            median_diff = median(x) - median(y),
            P_p1_gt_p2  = p_dominance(x, y),
            overlap     = if (isTRUE(do_overlap)) overlap_coef(x, y) else NA_real_
          )
        }, by=.(p1, p2)]
        out_phc <- file.path(output_dir, "EffDim_Phase_Contrasts.csv")
        data.table::fwrite(phase_contrasts, out_phc)
        register_output(out_phc)
        logf("Wrote EffDim_Phase_Contrasts.csv")
      }
    }
  }
}

# -----------------------------
# Plots (minimal; optional)
# -----------------------------
logf("Writing minimal summary plot (median ± 95% CI by day)...")
p <- ggplot(day_summ, aes(x=Day, y=median_eff_dim)) +
  geom_pointrange(aes(ymin=ci_lo_2p5, ymax=ci_hi_97p5)) +
  geom_line() +
  geom_point(size=2) +
  theme_minimal(base_size=14) +
  labs(title="Effective dimensionality by day (bootstrap median ± 95% CI)",
       x="Day", y="eff_dim (Neff)")

out_plot <- file.path(output_dir, "Plots", "EffDim_ByDay_Median_CI.png")
ggsave(out_plot, plot=p, width=9, height=5, dpi=300)
register_output(out_plot)
logf("Saved plot: Plots/EffDim_ByDay_Median_CI.png")

# -----------------------------
# Manifest (before HTML)
# -----------------------------
deps_dt <- deps_table(c("data.table","jsonlite","knitr","rmarkdown","ggplot2"))
run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

manifest <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(
    name = basename(commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))] |> sub("^--file=","", x=_) %||% "EffDim_Bootstrap_Analysis.R"),
    path = normalizePath(getwd(), winslash="/", mustWork=FALSE),
    full_path = NA
  ),
  input = list(
    infile = normalizePath(infile, winslash="/", mustWork=FALSE),
    days = as.integer(days_vec)
  ),
  parameters = list(
    analysis_name = analysis_name,
    source = SOURCE,
    do_skew = isTRUE(do_skew),
    do_overlap = isTRUE(do_overlap),
    contrast_mode = contrast_mode,
    custom_pairs = custom_pairs_text,
    do_phases = isTRUE(do_phases),
    phases = phase_text,
    do_html = isTRUE(do_html)
  ),
  dependencies = split(deps_dt, seq_len(nrow(deps_dt))),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash="/", mustWork=FALSE),
    output_dir = normalizePath(output_dir, winslash="/", mustWork=FALSE)
  ),
  generated_files = list()
)

write_manifest_json(manifest, manifest_json_path)
logf("Wrote Project_Manifest.json (initial)")

# Inventory pre-render
write_files_inventory(output_dir, manifest_files_csv)
register_output(manifest_files_csv)
logf("Wrote Project_Manifest_Files.csv (pre-render)")

# -----------------------------
# Optional QMD + HTML report
# -----------------------------
write_qmd_report <- function(qmd_path, script_label="EffDim Bootstrap Analysis", output_dir) {
  qmd <- c(
    "---",
    paste0("title: \"", script_label, " Report\""),
    "format:",
    "  html:",
    "    toc: true",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "params:",
    "  run_id: \"\"",
    "  run_timestamp: \"\"",
    "  output_dir: \"\"",
    "---",
    "",
    "```{r}",
    "suppressPackageStartupMessages({",
    "  library(data.table)",
    "  library(knitr)",
    "  library(ggplot2)",
    "})",
    "```",
    "",
    "## Metadata",
    "- **Run ID:** `r params$run_id`",
    "- **Run timestamp:** `r params$run_timestamp`",
    "- **Output directory:** `r params$output_dir`",
    "",
    "## Day-level summary",
    "```{r}",
    "f <- file.path(params$output_dir, 'EffDim_ByDay_Summary.csv')",
    "dt <- fread(f)",
    "kable(dt, format='html')",
    "```",
    "",
    "## Day contrasts",
    "```{r}",
    "f <- file.path(params$output_dir, 'EffDim_DayContrasts.csv')",
    "if (file.exists(f)) kable(fread(f), format='html') else cat('No contrast table found.\\n')",
    "```",
    "",
    "## Phase summaries (if requested)",
    "```{r}",
    "f1 <- file.path(params$output_dir, 'EffDim_Phase_Summary.csv')",
    "f2 <- file.path(params$output_dir, 'EffDim_Phase_Contrasts.csv')",
    "if (file.exists(f1)) {",
    "  cat('### Phase summary\\n')",
    "  kable(fread(f1), format='html')",
    "} else {",
    "  cat('No phase summary produced.\\n')",
    "}",
    "if (file.exists(f2)) {",
    "  cat('\\n\\n### Phase contrasts\\n')",
    "  kable(fread(f2), format='html')",
    "}",
    "```",
    "",
    "## Figure",
    "```{r}",
    "img <- file.path(params$output_dir, 'Plots', 'EffDim_ByDay_Median_CI.png')",
    "if (file.exists(img)) knitr::include_graphics(img)",
    "```",
    "",
    "## Files inventory",
    "```{r}",
    "inv <- file.path(params$output_dir, 'Project_Manifest_Files.csv')",
    "if (file.exists(inv)) kable(fread(inv), format='html') else cat('Inventory missing.\\n')",
    "```"
  )
  writeLines(qmd, qmd_path)
  invisible(qmd_path)
}

render_qmd_to_html <- function(qmd_path, html_path, params_list) {
  tryCatch({
    rmarkdown::render(
      input = qmd_path,
      output_file = basename(html_path),
      output_dir = dirname(html_path),
      quiet = TRUE,
      params = params_list
    )
  }, error=function(e) e)
}

if (isTRUE(do_html)) {
  logf("Writing QMD + rendering HTML...")
  qmd_path  <- file.path(output_dir, "EffDim_Bootstrap_Analysis_Report.qmd")
  html_path <- file.path(output_dir, "EffDim_Bootstrap_Analysis_Report.html")
  register_output(qmd_path); register_output(html_path)
  
  write_qmd_report(qmd_path, script_label="EffDim Bootstrap Analysis", output_dir=output_dir)
  res <- render_qmd_to_html(
    qmd_path, html_path,
    params_list = list(run_id=run_id, run_timestamp=run_timestamp, output_dir=normalizePath(output_dir, winslash="/", mustWork=FALSE))
  )
  if (inherits(res, "error") || !file.exists(html_path)) {
    logf(paste0("WARNING: HTML render failed: ", if (inherits(res,"error")) conditionMessage(res) else "unknown"))
    logf(paste0("QMD retained at: ", normalizePath(qmd_path, winslash="/", mustWork=FALSE)))
  } else {
    logf(paste0("HTML report created: ", normalizePath(html_path, winslash="/", mustWork=FALSE)))
  }
} else {
  logf("HTML report disabled by user choice.")
}

# -----------------------------
# Final inventory + manifest update
# -----------------------------
write_files_inventory(output_dir, manifest_files_csv)
logf("Refreshed Project_Manifest_Files.csv (post-run)")

inv_dt <- data.table::fread(manifest_files_csv)
manifest$generated_files <- split(inv_dt[, .(relative_path, full_path, bytes, modified)], seq_len(nrow(inv_dt)))
write_manifest_json(manifest, manifest_json_path)
logf("Updated Project_Manifest.json (with generated_files)")

logf("DONE")

cat("\n==============================================\n")
cat("COMPLETE\n")
cat("Run ID:\n  ", run_id, "\n", sep="")
cat("Outputs written under:\n  ", normalizePath(output_dir, winslash="/", mustWork=FALSE), "\n", sep="")
cat("\nKey outputs:\n")
cat("  - EffDim_ByDay_Summary.csv\n")
cat("  - EffDim_DayContrasts.csv\n")
cat("  - (optional) EffDim_Phase_Summary.csv\n")
cat("  - (optional) EffDim_Phase_Contrasts.csv\n")
cat("  - Plots/EffDim_ByDay_Median_CI.png\n")
cat("  - Project_Manifest.json\n")
cat("  - Project_Manifest_Files.csv\n")
cat("Log:\n  - ", normalizePath(logfile, winslash="/", mustWork=FALSE), "\n", sep="")
cat("==============================================\n\n")

# small helper for manifest script field (safe fallback)
`%||%` <- function(a,b) if (!is.null(a) && length(a) && is_nonempty_string(a)) a else b
