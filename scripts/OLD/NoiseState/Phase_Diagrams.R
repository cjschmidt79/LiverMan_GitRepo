
#!/usr/bin/env Rscript
# ======================================================================
# Constructed Phase Diagrams (NSV day-level summaries)
#
# WHAT THIS SCRIPT IS GOOD FOR
#   This script constructs a descriptive “state-space / phase-style” trajectory from
#   day-level summary metrics such as D(t) and ED(t), optionally with uncertainty
#   envelopes (lo/hi intervals), and with delta diagnostics:
#
#     1) Phase trajectory:         (D(t), ED(t)) ordered by Day, with arrows.
#     2) Uncertainty envelopes:    rectangles using (D_lo, D_hi) × (ED_lo, ED_hi).
#     3) ΔX vs ΔY panel:           scatter of ΔD vs ΔED between consecutive days.
#     4) Δ-vector arrows panel:    vectors from origin (0,0) to (ΔD, ΔED).
#     5) Faceting across modes:    RAW / RESID_DAY / RESID_DAY_PCk (if provided).
#
# IMPORTANT SCIENTIFIC NOTE
#   This is a CONSTRUCTED visualization. It does not assume continuous dynamics.
#   Arrows indicate the ordering of observed timepoints only.
#
# WHAT YOU WILL PROVIDE (INPUT TABLE REQUIREMENTS)
#   You will select one or more CSV files containing day-level summaries.
#
#   Minimum required columns (per file):
#     - Day (numeric or coercible to numeric)
#     - D central estimate (e.g., D_med or D or D_mid or D_mean)
#     - ED central estimate (e.g., ED_med or ED or ED_mid or ED_mean)
#
#   Optional columns for uncertainty rectangles:
#     - D_lo and D_hi
#     - ED_lo and ED_hi
#
#   Optional column for faceting if present in the file:
#     - Mode (e.g., RAW / RESID_DAY / RESID_DAY_PC2 ...)
#
#   If Mode is not present:
#     - You will be prompted once for a mode label per file (e.g., "RAW").
#
# OUTPUTS (ALL WRITTEN UNDER outputs/<run_id>/)
#   - Phase_D_vs_ED_faceted.png/.pdf
#   - Phase_D_vs_ED_faceted_with_uncertainty.png/.pdf (if intervals available + enabled)
#   - Delta_D_vs_ED_faceted.png/.pdf
#   - DeltaVectors_DeltaD_DeltaED_faceted.png/.pdf
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - <script_name>_Report.qmd  (rendered to HTML as final step)
#   - <script_name>_Report.html (if render succeeds)
#
# ======================================================================

# ============================ SCRIPT IDENTITY (MANDATORY BLOCK) ============================

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
known_script_filename <- "PhaseDiagram_Constructed_WithDeltas.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  # Path cannot be detected in this execution mode; still record a valid script_name.
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  script_full_note <- "Script path detection failed; fallback filename was used."
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
  script_full_note <- ""
}

# Special-case: RStudio sometimes runs via ~/.active-rstudio-document
fallback_used_for_reporting <- FALSE
if (!is.na(script_full) && grepl("\\.active-rstudio-document$", script_full)) {
  fallback_used_for_reporting <- TRUE
}

cat("\n================ Script Identity ================\n")
cat("Script name :", script_name, "\n")
cat("Script path :", script_path, "\n")
cat("Script full :", script_full, "\n")
if (nzchar(script_full_note)) cat("NOTE:", script_full_note, "\n")
if (fallback_used_for_reporting) {
  cat("NOTE: Running from .active-rstudio-document; manifest will record this and also record known_script_filename.\n")
}
cat("=================================================\n\n")

# ============================ UTILITIES ============================

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_stem <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

norm_name <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

# ============================ OUTPUT DIRECTORY POLICY (MANDATORY) ============================

analysis_name <- "Phase"
outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

run_id <- paste0(analysis_name, "_NSV_", ts_stamp())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================ PACKAGES (MANDATORY) ============================

needed <- c("ggplot2", "dplyr", "tidyr", "jsonlite", "knitr", "rmarkdown")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!has_patchwork) {
  message("Optional package not found: patchwork. Multi-panel composites will be skipped.")
}

# FIXED: packageVersion() returns a package_version object; coerce to character explicitly.
deps_df <- data.frame(
  package = needed,
  version = vapply(needed, function(p) as.character(utils::packageVersion(p)), character(1)),
  stringsAsFactors = FALSE
)

# ============================ MANIFEST HELPERS ============================

write_inventory_csv <- function(out_dir, csv_path) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(out_dir, winslash = "/")), "/?"), "", gsub("\\\\", "/", files))
  info <- file.info(files)
  inv <- data.frame(
    file = rel,
    size_bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(inv, csv_path, row.names = FALSE)
  inv
}

write_manifest_json <- function(path, manifest_list) {
  jsonlite::write_json(manifest_list, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

# ============================ INPUT INGESTION ============================

choose_csv_one <- function(prompt = "Choose a day-level summary CSV") {
  if (!interactive()) stop2("Run interactively (RStudio) so file.choose() can be used, or add CLI args.")
  message(prompt)
  file.choose()
}

read_one_phase_csv <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df
}

detect_col <- function(df, candidates) {
  nn <- norm_name(names(df))
  cand_norm <- norm_name(candidates)
  hit <- which(nn %in% cand_norm)
  if (length(hit) >= 1) names(df)[hit[1]] else NULL
}

detect_metric_triplet <- function(df, metric_prefix) {
  # e.g., metric_prefix="D" expects D_med/D_mid/D_mean + D_lo + D_hi
  nms <- names(df)
  nn <- norm_name(nms)
  
  pref <- norm_name(metric_prefix)
  is_pref <- function(x) startsWith(x, pref)
  
  mid_candidates <- c(
    paste0(metric_prefix, "_med"),
    paste0(metric_prefix, "_mid"),
    paste0(metric_prefix, "_mean"),
    metric_prefix
  )
  lo_candidates <- c(paste0(metric_prefix, "_lo"), paste0(metric_prefix, "_lower"))
  hi_candidates <- c(paste0(metric_prefix, "_hi"), paste0(metric_prefix, "_upper"))
  
  mid <- detect_col(df, mid_candidates)
  lo  <- detect_col(df, lo_candidates)
  hi  <- detect_col(df, hi_candidates)
  
  list(mid = mid, lo = lo, hi = hi)
}

prompt_for_colname <- function(df, prompt_text) {
  cat("\n", prompt_text, "\n", sep = "")
  cat("Columns:\n")
  for (i in seq_along(names(df))) cat(sprintf("  [%d] %s\n", i, names(df)[i]))
  raw <- readline("Type the column NAME exactly (or an index number): ")
  raw <- trimws(raw)
  if (!nzchar(raw)) stop2("No selection provided.")
  idx <- suppressWarnings(as.integer(raw))
  if (!is.na(idx)) {
    if (idx < 1 || idx > ncol(df)) stop2("Index out of range.")
    return(names(df)[idx])
  }
  if (!(raw %in% names(df))) stop2("Column name not found: ", raw)
  raw
}

# Load one or more files, bind into a single long table with Mode
load_inputs <- function() {
  cat("=== Input selection ===\n")
  cat("You will select one or more day-level CSV(s).\n")
  cat("If a CSV does not contain a Mode column, you will be asked for a mode label.\n\n")
  
  dfs <- list()
  paths <- character(0)
  
  repeat {
    p <- choose_csv_one("Select a day-level summary CSV (cancel to stop selection is NOT supported here).")
    if (!nzchar(p)) stop2("No file selected.")
    paths <- c(paths, p)
    df <- read_one_phase_csv(p)
    
    # Mode
    mode_col <- detect_col(df, c("Mode", "mode"))
    if (is.null(mode_col)) {
      mode_label <- readline(paste0("Mode label for this file (e.g., RAW / RESID_DAY / RESID_DAY_PC2) [RAW]: "))
      if (!nzchar(mode_label)) mode_label <- "RAW"
      df$Mode <- mode_label
    } else {
      # standardize
      df$Mode <- as.character(df[[mode_col]])
    }
    
    df$.__source_file <- normalizePath(p, winslash = "/", mustWork = FALSE)
    
    dfs[[length(dfs) + 1]] <- df
    
    add_more <- readline("Add another CSV? (y/n) [n]: ")
    if (!nzchar(add_more)) add_more <- "n"
    if (tolower(add_more) != "y") break
  }
  
  bind <- dplyr::bind_rows(dfs)
  list(df = bind, paths = paths)
}

# ============================ PLOTTING THEME ============================

theme_pub <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8)
    )
}

# ============================ MAIN ANALYSIS ============================

cat("\n=== Phase Diagrams with Δ Diagnostics ===\n")
cat("Outputs will be written to:\n  ", output_dir, "\n\n", sep = "")

inputs <- load_inputs()
raw_df <- inputs$df
input_paths <- inputs$paths

# Detect Day / D / ED columns (auto; prompt if missing)
day_col <- detect_col(raw_df, c("Day", "day", "time", "Time", "t"))
if (is.null(day_col)) day_col <- prompt_for_colname(raw_df, "Could not auto-detect Day column. Select Day:")

d_trip <- detect_metric_triplet(raw_df, "D")
ed_trip <- detect_metric_triplet(raw_df, "ED")

if (is.null(d_trip$mid))  d_trip$mid  <- prompt_for_colname(raw_df, "Could not auto-detect D central column. Select D central (plotted X):")
if (is.null(ed_trip$mid)) ed_trip$mid <- prompt_for_colname(raw_df, "Could not auto-detect ED central column. Select ED central (plotted Y):")

# Coerce numeric
to_num <- function(x) suppressWarnings(as.numeric(x))
raw_df[[day_col]] <- to_num(raw_df[[day_col]])
raw_df[[d_trip$mid]] <- to_num(raw_df[[d_trip$mid]])
raw_df[[ed_trip$mid]] <- to_num(raw_df[[ed_trip$mid]])

if (!is.null(d_trip$lo)) raw_df[[d_trip$lo]] <- to_num(raw_df[[d_trip$lo]])
if (!is.null(d_trip$hi)) raw_df[[d_trip$hi]] <- to_num(raw_df[[d_trip$hi]])
if (!is.null(ed_trip$lo)) raw_df[[ed_trip$lo]] <- to_num(raw_df[[ed_trip$lo]])
if (!is.null(ed_trip$hi)) raw_df[[ed_trip$hi]] <- to_num(raw_df[[ed_trip$hi]])

# Clean / distinct per Mode-Day
df <- raw_df %>%
  filter(
    is.finite(.data[[day_col]]),
    is.finite(.data[[d_trip$mid]]),
    is.finite(.data[[ed_trip$mid]])
  ) %>%
  mutate(
    Day = .data[[day_col]],
    D_mid = .data[[d_trip$mid]],
    ED_mid = .data[[ed_trip$mid]]
  ) %>%
  distinct(Mode, Day, .keep_all = TRUE) %>%
  arrange(Mode, Day)

if (nrow(df) < 3) stop2("Need at least 3 (Mode,Day) rows in total to construct trajectories.")

# Optional uncertainty
has_uncertainty <- !is.null(d_trip$lo) && !is.null(d_trip$hi) && !is.null(ed_trip$lo) && !is.null(ed_trip$hi)
use_uncertainty <- FALSE
if (has_uncertainty) {
  ans <- readline("Uncertainty rectangles available (D_lo/hi and ED_lo/hi). Include them? (y/n) [y]: ")
  if (!nzchar(ans)) ans <- "y"
  use_uncertainty <- (tolower(ans) == "y")
  df <- df %>%
    mutate(
      D_lo = .data[[d_trip$lo]],
      D_hi = .data[[d_trip$hi]],
      ED_lo = .data[[ed_trip$lo]],
      ED_hi = .data[[ed_trip$hi]]
    )
} else {
  cat("NOTE: Uncertainty rectangles not enabled (missing one or more of D_lo/D_hi/ED_lo/ED_hi).\n")
}

# Build consecutive-day segments (per Mode)
seg <- df %>%
  group_by(Mode) %>%
  arrange(Day, .by_group = TRUE) %>%
  mutate(
    Day_next = dplyr::lead(Day),
    D_next   = dplyr::lead(D_mid),
    ED_next  = dplyr::lead(ED_mid)
  ) %>%
  ungroup() %>%
  filter(is.finite(D_next), is.finite(ED_next), is.finite(Day_next))

# Compute deltas per Mode (Δ relative to previous day)
df_delta <- df %>%
  group_by(Mode) %>%
  arrange(Day, .by_group = TRUE) %>%
  mutate(
    DeltaD  = D_mid  - dplyr::lag(D_mid),
    DeltaED = ED_mid - dplyr::lag(ED_mid),
    DeltaDay = Day - dplyr::lag(Day)
  ) %>%
  ungroup()

# ============================ PLOTS ============================

plot_title <- "Constructed state-space trajectory: (D(t), ED(t)) ordered by Day"
plot_subtitle <- paste0(
  "X = ", d_trip$mid, " (central) | Y = ", ed_trip$mid, " (central) | arrows connect consecutive days"
)

# 1) Phase plot (faceted)
p_phase <- ggplot(df, aes(x = D_mid, y = ED_mid)) +
  geom_segment(
    data = seg,
    aes(x = D_mid, y = ED_mid, xend = D_next, yend = ED_next),
    arrow = arrow(type = "closed", length = grid::unit(3, "mm")),
    linewidth = 0.8,
    lineend = "round"
  ) +
  geom_point(size = 2.6) +
  {
    if (has_ggrepel) {
      ggrepel::geom_text_repel(aes(label = Day), size = 4, max.overlaps = Inf)
    } else {
      geom_text(aes(label = Day), vjust = -0.8, size = 4)
    }
  } +
  facet_wrap(~ Mode, scales = "free") +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "D(t) (central estimate)",
    y = "ED(t) (central estimate)"
  ) +
  theme_pub()

# 1b) Phase plot with uncertainty rectangles (optional)
p_phase_unc <- NULL
if (use_uncertainty) {
  p_phase_unc <- ggplot(df, aes(x = D_mid, y = ED_mid)) +
    geom_rect(
      aes(xmin = D_lo, xmax = D_hi, ymin = ED_lo, ymax = ED_hi),
      fill = "grey80", alpha = 0.55, color = NA
    ) +
    geom_segment(
      data = seg,
      aes(x = D_mid, y = ED_mid, xend = D_next, yend = ED_next),
      arrow = arrow(type = "closed", length = grid::unit(3, "mm")),
      linewidth = 0.8,
      lineend = "round"
    ) +
    geom_point(size = 2.6) +
    {
      if (has_ggrepel) {
        ggrepel::geom_text_repel(aes(label = Day), size = 4, max.overlaps = Inf)
      } else {
        geom_text(aes(label = Day), vjust = -0.8, size = 4)
      }
    } +
    facet_wrap(~ Mode, scales = "free") +
    labs(
      title = "Constructed phase diagram with uncertainty envelopes",
      subtitle = "Rectangles are marginal intervals (D_lo..D_hi × ED_lo..ED_hi). Arrows connect consecutive days.",
      x = "D(t) (central estimate)",
      y = "ED(t) (central estimate)"
    ) +
    theme_pub()
}

# 2) ΔD vs ΔED panel
p_delta <- ggplot(df_delta %>% filter(is.finite(DeltaD), is.finite(DeltaED)), aes(x = DeltaD, y = DeltaED)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_point(size = 2.4) +
  {
    if (has_ggrepel) {
      ggrepel::geom_text_repel(aes(label = Day), size = 4, max.overlaps = Inf)
    } else {
      geom_text(aes(label = Day), vjust = -0.8, size = 4)
    }
  } +
  facet_wrap(~ Mode, scales = "free") +
  labs(
    title = "Δ panel: ΔD vs ΔED (consecutive days)",
    subtitle = "Points labeled by destination Day (current day).",
    x = "ΔD = D(t) - D(t-1)",
    y = "ΔED = ED(t) - ED(t-1)"
  ) +
  theme_pub()

# 3) Δ-vector arrows from origin
df_vec <- df_delta %>% filter(is.finite(DeltaD), is.finite(DeltaED))
p_vec <- ggplot(df_vec, aes(x = 0, y = 0)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_segment(
    aes(xend = DeltaD, yend = DeltaED),
    arrow = arrow(type = "closed", length = grid::unit(3, "mm")),
    linewidth = 0.8,
    lineend = "round"
  ) +
  geom_point(aes(x = DeltaD, y = DeltaED), size = 2.4) +
  {
    if (has_ggrepel) {
      ggrepel::geom_text_repel(aes(x = DeltaD, y = DeltaED, label = Day), size = 4, max.overlaps = Inf)
    } else {
      geom_text(aes(x = DeltaD, y = DeltaED, label = Day), vjust = -0.8, size = 4)
    }
  } +
  facet_wrap(~ Mode, scales = "free") +
  labs(
    title = "Δ-vector arrows: (0,0) → (ΔD, ΔED)",
    subtitle = "Vectors represent the stepwise change between consecutive days; labeled by destination Day.",
    x = "ΔD",
    y = "ΔED"
  ) +
  theme_pub()

# ============================ SAVE PLOTS ============================

save_plot <- function(p, stem, w = 7.2, h = 6.0) {
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, p, width = w, height = h, dpi = 600)
  ggplot2::ggsave(pdf_path, p, width = w, height = h)
  list(png = png_path, pdf = pdf_path)
}

files_out <- list()

files_out$phase <- save_plot(p_phase, "Phase_D_vs_ED_faceted", w = 8.2, h = 6.2)
if (!is.null(p_phase_unc)) {
  files_out$phase_uncertainty <- save_plot(p_phase_unc, "Phase_D_vs_ED_faceted_with_uncertainty", w = 8.2, h = 6.2)
}
files_out$delta <- save_plot(p_delta, "Delta_D_vs_ED_faceted", w = 8.2, h = 6.2)
files_out$delta_vectors <- save_plot(p_vec, "DeltaVectors_DeltaD_DeltaED_faceted", w = 8.2, h = 6.2)

cat("\nSaved plots:\n")
for (nm in names(files_out)) {
  cat(" - ", nm, ":\n", sep = "")
  cat("     ", files_out[[nm]]$png, "\n", sep = "")
  cat("     ", files_out[[nm]]$pdf, "\n", sep = "")
}

# ============================ MANIFEST (MANDATORY) ============================

manifest_path <- file.path(output_dir, "Project_Manifest.json")
inventory_path <- file.path(output_dir, "Project_Manifest_Files.csv")

manifest_list <- list(
  run_id = run_id,
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  script = list(
    name = script_name,
    path = script_path,
    full_path = if (is.na(script_full)) NA_character_ else script_full,
    known_script_filename = known_script_filename,
    note = if (fallback_used_for_reporting) "Executed via .active-rstudio-document; known_script_filename recorded for stable identity." else script_full_note
  ),
  input = list(
    csv_paths = vapply(input_paths, function(p) normalizePath(p, winslash = "/", mustWork = FALSE), character(1)),
    n_rows_combined = nrow(raw_df),
    n_rows_distinct_mode_day = nrow(df)
  ),
  parameters = list(
    analysis_name = analysis_name,
    day_col_detected = day_col,
    d_central_col = d_trip$mid,
    d_lo_col = d_trip$lo %||% NA_character_,
    d_hi_col = d_trip$hi %||% NA_character_,
    ed_central_col = ed_trip$mid,
    ed_lo_col = ed_trip$lo %||% NA_character_,
    ed_hi_col = ed_trip$hi %||% NA_character_,
    uncertainty_rectangles_included = use_uncertainty,
    faceting = "Mode"
  ),
  dependencies = lapply(seq_len(nrow(deps_df)), function(i) list(package = deps_df$package[i], version = deps_df$version[i])),
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
)

write_manifest_json(manifest_path, manifest_list)

# ============================ QUARTO REPORT (QMD + RENDER) ============================

# Embed the initial comment block (hard-coded from this script header for robustness)
script_header_text <- paste(
  readLines(textConnection(c(
    "# ======================================================================",
    "# Constructed Phase Diagrams (NSV day-level summaries)",
    "#",
    "# WHAT THIS SCRIPT IS GOOD FOR",
    "#   This script constructs a descriptive “state-space / phase-style” trajectory from",
    "#   day-level summary metrics such as D(t) and ED(t), optionally with uncertainty",
    "#   envelopes (lo/hi intervals), and with delta diagnostics:",
    "#",
    "#     1) Phase trajectory:         (D(t), ED(t)) ordered by Day, with arrows.",
    "#     2) Uncertainty envelopes:    rectangles using (D_lo, D_hi) × (ED_lo, ED_hi).",
    "#     3) ΔX vs ΔY panel:           scatter of ΔD vs ΔED between consecutive days.",
    "#     4) Δ-vector arrows panel:    vectors from origin (0,0) to (ΔD, ΔED).",
    "#     5) Faceting across modes:    RAW / RESID_DAY / RESID_DAY_PCk (if provided).",
    "#",
    "# IMPORTANT SCIENTIFIC NOTE",
    "#   This is a CONSTRUCTED visualization. It does not assume continuous dynamics.",
    "#   Arrows indicate the ordering of observed timepoints only.",
    "#",
    "# WHAT YOU WILL PROVIDE (INPUT TABLE REQUIREMENTS)",
    "#   You will select one or more CSV files containing day-level summaries.",
    "#",
    "#   Minimum required columns (per file):",
    "#     - Day (numeric or coercible to numeric)",
    "#     - D central estimate (e.g., D_med or D or D_mid or D_mean)",
    "#     - ED central estimate (e.g., ED_med or ED or ED_mid or ED_mean)",
    "#",
    "#   Optional columns for uncertainty rectangles:",
    "#     - D_lo and D_hi",
    "#     - ED_lo and ED_hi",
    "#",
    "#   Optional column for faceting if present in the file:",
    "#     - Mode (e.g., RAW / RESID_DAY / RESID_DAY_PC2 ...)",
    "#",
    "# ======================================================================"
  ))),
  collapse = "\n"
)

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_out <- file.path(output_dir, paste0(script_name, "_Report.html"))

# Pick a plot to always include
main_plot_path <- files_out$phase$png

qmd_lines <- c(
  "---",
  paste0('title: "', script_name, ' Report"'),
  "format: html",
  "toc: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "## Script header",
  "",
  "```",
  script_header_text,
  "```",
  "",
  "## Run metadata",
  "",
  "```{r}",
  "library(jsonlite)",
  "library(knitr)",
  paste0("manifest <- jsonlite::fromJSON('", basename(manifest_path), "', simplifyVector = FALSE)"),
  "meta_tbl <- data.frame(",
  "  Field = c('run_id','run_timestamp','script_name','script_path','script_full_path','outputs_root','output_dir'),",
  "  Value = c(",
  "    manifest[['run_id']],",
  "    manifest[['run_timestamp']],",
  "    manifest[['script']][['name']],",
  "    manifest[['script']][['path']],",
  "    ifelse(is.null(manifest[['script']][['full_path']]) || is.na(manifest[['script']][['full_path']]), 'NA', manifest[['script']][['full_path']]),",
  "    manifest[['outputs']][['outputs_root']],",
  "    manifest[['outputs']][['output_dir']]",
  "  ),",
  "  stringsAsFactors = FALSE",
  ")",
  "knitr::kable(meta_tbl)",
  "```",
  "",
  "## Inputs and parameters",
  "",
  "```{r}",
  "inp <- manifest[['input']]",
  "par <- manifest[['parameters']]",
  "inp_tbl <- data.frame(Field = names(inp), Value = vapply(inp, function(x) paste(x, collapse = '; '), character(1)), stringsAsFactors = FALSE)",
  "par_tbl <- data.frame(Field = names(par), Value = vapply(par, function(x) as.character(x), character(1)), stringsAsFactors = FALSE)",
  "knitr::kable(inp_tbl)",
  "knitr::kable(par_tbl)",
  "```",
  "",
  "## Dependencies",
  "",
  "```{r}",
  "deps <- do.call(rbind, lapply(manifest[['dependencies']], function(x) data.frame(package = x$package, version = x$version, stringsAsFactors = FALSE)))",
  "knitr::kable(deps)",
  "```",
  "",
  "## Analytical logic",
  "",
  "- Let **D(t)** be a day-level dispersion summary (central estimate) and **ED(t)** be day-level effective dimensionality.",
  "- The constructed state-space trajectory plots points **(D(t), ED(t))** and draws arrows from consecutive days.",
  "- Stepwise changes are summarized as:",
  "  - **ΔD = D(t) − D(t−1)**",
  "  - **ΔED = ED(t) − ED(t−1)**",
  "",
  "If interval columns exist (lo/hi), uncertainty is shown as marginal rectangles:",
  "- x-range: [D_lo, D_hi]",
  "- y-range: [ED_lo, ED_hi]",
  "",
  "## Key plots",
  "",
  "```{r}",
  paste0("knitr::include_graphics('", basename(main_plot_path), "')"),
  "```",
  "",
  "## Interpretation guidance (descriptive)",
  "",
  "- Movement along the x-axis indicates changes in D(t) (dispersion).",
  "- Movement along the y-axis indicates changes in ED(t) (effective dimensionality).",
  "- Clusters of nearby points imply similar joint (D, ED) state across adjacent days.",
  "- Large Δ-vectors indicate abrupt stepwise changes between consecutive timepoints.",
  "",
  "## Generated files",
  "",
  "```{r}",
  paste0("inv <- read.csv('", basename(inventory_path), "', stringsAsFactors = FALSE)"),
  "knitr::kable(inv)",
  "```",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd_lines, qmd_path)

# Render QMD to HTML as final step (MANDATORY)
render_ok <- TRUE
render_msg <- NULL
tryCatch({
  # Render into output_dir with explicit output_file name
  rmarkdown::render(
    input = qmd_path,
    output_file = basename(html_out),
    output_dir = output_dir,
    quiet = TRUE
  )
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# Refresh inventory AFTER render attempt (MANDATORY)
inv_df <- write_inventory_csv(output_dir, inventory_path)

# Update manifest with generated_files reference (and rewrite JSON)
manifest_list$generated_files <- list(
  inventory_csv = basename(inventory_path),
  n_files = nrow(inv_df)
)
write_manifest_json(manifest_path, manifest_list)

cat("\n================ FINAL STATUS ================\n")
cat("Manifest JSON: ", manifest_path, "\n", sep = "")
cat("File inventory: ", inventory_path, "\n", sep = "")
cat("QMD written: ", qmd_path, "\n", sep = "")

if (render_ok && file.exists(html_out)) {
  cat("HTML report created: ", normalizePath(html_out, winslash = "/", mustWork = FALSE), "\n", sep = "")
} else {
  cat("HTML report render FAILED.\n")
  if (!is.null(render_msg)) cat("Render error: ", render_msg, "\n", sep = "")
  cat("You can try rendering manually from R:\n")
  cat("  rmarkdown::render('", qmd_path, "', output_dir='", output_dir, "')\n", sep = "")
}
cat("Outputs directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("=============================================\n\n")

cat("Done.\n")
