#!/usr/bin/env Rscript
# ======================================================================
# Cross-layer State-Vector Similarity (CLSVS) Coupling Analysis
#
# PURPOSE
#   Compute a simple, explicit coupling metric between metabolome and
#   transcriptome organizational dynamics across days.
#
# CORE IDEA
#   Each layer is represented at each day by a 2D "organization state vector":
#     s(d) = [ z(D(d)), z(-ED(d)) ]
#   where:
#     - D(d)  = dispersion metric at day d (e.g., median CV / IQR / robust Dt)
#     - ED(d) = effective dimensionality metric at day d
#     - z(.)  = z-score across days within the same layer
#
#   Coupling at metabolome day d with lag L is cosine similarity:
#     C(d;L) = cos_sim( s_met(d), s_tx(d+L) ) in [-1, 1]
#
# INPUTS (provided interactively)
#   1) Metabolome day-metrics CSV with columns:
#        Day, D, ED
#   2) Transcriptome day-metrics CSV with columns:
#        Day, D, ED
#
# OUTPUTS (always under outputs/<run_id>/)
#   - Joint_CLSVS_ByDay.csv
#   - CLSVS_Lag_Summary.csv
#   - CLSVS_Coupling_ByDay.png
#   - CLSVS_StateSpace_Scatter.png
#   - Project_Manifest.json
#   - Project_Manifest_Files.csv
#   - <script_name>_Report.qmd and rendered HTML
#
# NOTES
#   - This script does not claim causality. It quantifies state similarity
#     across layers with explicit lags.
# ======================================================================

# ---------------------------
# 0) Dependency handling
# ---------------------------
required_pkgs <- c("jsonlite", "knitr", "rmarkdown", "ggplot2", "dplyr", "readr", "tibble", "tools")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

stop2 <- function(...) stop(paste0(...), call. = FALSE)

# ---------------------------
# 1) Script identity capture (MANDATORY block)
# ---------------------------
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
known_script_filename <- "CrossLayerCoupling_CLSVS.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()

if (is.na(script_full)) {
  # Path cannot be detected in this execution mode; still record a valid script_name.
  script_name <- known_script_stem
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_name <- tools::file_path_sans_ext(basename(script_full))
  script_path <- normalizePath(dirname(script_full), winslash = "/", mustWork = FALSE)
}

# Print script identity in console summary (MANDATORY)
message("Script identity:")
message("  script_name: ", script_name)
message("  script_path: ", script_path)
message("  script_full: ", ifelse(is.na(script_full), "NA", script_full))
if (is.na(script_full)) {
  message("NOTE: Script path detection failed; using fallback filename: ", known_script_filename)
}

# ---------------------------
# 2) Output directory policy (MANDATORY)
# ---------------------------
analysis_name <- "CrossLayerCoupling"
source_tag <- "Metab_Tx"
run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_id <- paste0(analysis_name, "_", source_tag, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))

outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Outputs:")
message("  outputs_root: ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE))
message("  output_dir:   ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
setwd(output_dir)  # enforce “outputs only” writing

# ---------------------------
# 3) Inputs (interactive file selection)
# ---------------------------
message("\nChoose METABOLOME day-level metrics CSV (must include columns: Day, D, ED)")
met_path <- file.choose()
message("MET input: ", normalizePath(met_path, winslash = "/", mustWork = FALSE))

message("\nChoose TRANSCRIPTOME day-level metrics CSV (must include columns: Day, D, ED)")
tx_path <- file.choose()
message("TX input:  ", normalizePath(tx_path, winslash = "/", mustWork = FALSE))

# ---------------------------
# 4) Load inputs
# ---------------------------
met <- readr::read_csv(met_path, show_col_types = FALSE)
tx  <- readr::read_csv(tx_path,  show_col_types = FALSE)

needed_cols <- c("Day", "D", "ED")
if (!all(needed_cols %in% names(met))) stop2("Metabolome CSV must contain columns: Day, D, ED")
if (!all(needed_cols %in% names(tx)))  stop2("Transcriptome CSV must contain columns: Day, D, ED")

met <- met %>% mutate(Day = as.integer(Day)) %>% arrange(Day)
tx  <- tx  %>% mutate(Day = as.integer(Day)) %>% arrange(Day)

# ---------------------------
# 5) CLSVS computation (do not change scientific logic)
# ---------------------------
z <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

state_vec_tbl <- function(df, layer_label) {
  df %>%
    transmute(
      Day = as.integer(Day),
      layer = layer_label,
      D = as.numeric(D),
      ED = as.numeric(ED),
      zD = z(D),
      zNegED = z(-ED)
    )
}

S_met <- state_vec_tbl(met, "metabolome")
S_tx  <- state_vec_tbl(tx,  "transcriptome")

cos_sim <- function(a1, a2, b1, b2) {
  # cosine similarity between (a1, a2) and (b1, b2)
  num <- a1*b1 + a2*b2
  den <- sqrt(a1^2 + a2^2) * sqrt(b1^2 + b2^2)
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

lags_to_test <- c(0, 2, 4)

# Join by explicit lag: compare met at Day d to tx at Day d+L
coupling_by_lag <- lapply(lags_to_test, function(L) {
  tx_shift <- S_tx %>% transmute(Day = Day - L, zD_tx = zD, zNegED_tx = zNegED)  # so Day aligns with met Day
  S_met %>%
    left_join(tx_shift, by = "Day") %>%
    mutate(
      lag_days = L,
      C = cos_sim(zD, zNegED, zD_tx, zNegED_tx),
      pair_day_tx = Day + L
    ) %>%
    select(Day_met = Day, Day_tx = pair_day_tx, lag_days, C, zD_met = zD, zNegED_met = zNegED, zD_tx, zNegED_tx)
})

coupling_tbl <- bind_rows(coupling_by_lag) %>% arrange(lag_days, Day_met)

# Summary by lag
lag_summary <- coupling_tbl %>%
  group_by(lag_days) %>%
  summarize(
    n_pairs = sum(!is.na(C)),
    mean_C = mean(C, na.rm = TRUE),
    median_C = median(C, na.rm = TRUE),
    .groups = "drop"
  )

# Write core tables
write_csv(coupling_tbl, "Joint_CLSVS_ByDay.csv")
write_csv(lag_summary, "CLSVS_Lag_Summary.csv")

# ---------------------------
# 6) Plots (written into output_dir)
# ---------------------------
# Plot 1: Coupling vs metabolome day for each lag
p1 <- ggplot(coupling_tbl, aes(x = Day_met, y = C, group = factor(lag_days))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ lag_days, ncol = 1, scales = "fixed") +
  scale_x_continuous(breaks = sort(unique(coupling_tbl$Day_met))) +
  labs(
    title = "Cross-layer State-Vector Similarity (CLSVS) by day",
    subtitle = "Cosine similarity between metabolome state at Day d and transcriptome state at Day d + lag",
    x = "Metabolome day (d)",
    y = "CLSVS coupling C(d;lag)"
  ) +
  theme_classic(base_size = 12)

ggsave("CLSVS_Coupling_ByDay.png", plot = p1, width = 7.5, height = 8.5, dpi = 300)

# Plot 2: State-space scatter (zD vs z(-ED)) for each layer (no coupling, just geometry)
state_all <- bind_rows(S_met, S_tx)
p2 <- ggplot(state_all, aes(x = zD, y = zNegED)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_path(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_text(aes(label = Day), vjust = -0.8, size = 3, check_overlap = TRUE) +
  facet_wrap(~ layer, ncol = 1) +
  labs(
    title = "Organization state-space by layer",
    subtitle = "State vector components: z(D) and z(-ED) computed within each layer across days",
    x = "z(D)",
    y = "z(-ED)"
  ) +
  theme_classic(base_size = 12)

ggsave("CLSVS_StateSpace_Scatter.png", plot = p2, width = 7.5, height = 8.5, dpi = 300)

# ---------------------------
# 7) Manifest + file inventory (MANDATORY)
# ---------------------------
deps <- lapply(required_pkgs, function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  list(package = p, version = v)
})

# Inventory function
inventory_outputs <- function(dir_path) {
  files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  tibble(
    file = gsub(paste0("^", gsub("\\\\", "/", normalizePath(dir_path, winslash = "/", mustWork = FALSE)), "/?"), "", gsub("\\\\", "/", files)),
    full_path = gsub("\\\\", "/", files),
    bytes = file.info(files)$size,
    modified = as.character(file.info(files)$mtime)
  ) %>% arrange(file)
}

# Initial inventory (pre-report)
inv_pre <- inventory_outputs(output_dir)
write_csv(inv_pre, "Project_Manifest_Files.csv")

manifest <- list(
  run_id = run_id,
  run_timestamp = run_timestamp,
  script = list(
    name = script_name,
    path = script_path,
    full_path = ifelse(is.na(script_full), NA, script_full)
  ),
  input = list(
    metabolome_csv = normalizePath(met_path, winslash = "/", mustWork = FALSE),
    transcriptome_csv = normalizePath(tx_path, winslash = "/", mustWork = FALSE)
  ),
  parameters = list(
    analysis_name = analysis_name,
    source_tag = source_tag,
    lags_tested_days = lags_to_test,
    state_definition = "s(d) = [ z(D(d)), z(-ED(d)) ]",
    coupling_definition = "C(d;L) = cosine_similarity( s_met(d), s_tx(d+L) )"
  ),
  dependencies = deps,
  outputs = list(
    outputs_root = normalizePath(outputs_root, winslash = "/", mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )
)

write_json(manifest, "Project_Manifest.json", pretty = TRUE, auto_unbox = TRUE, na = "null")

# ---------------------------
# 8) Quarto QMD (MANDATORY) + render to HTML (MANDATORY final step)
# ---------------------------
# Embed the script header text as a fenced code block.
# Prefer reading the actual script if possible; otherwise use a fallback header string.
extract_header_block <- function(path) {
  if (!is.na(path) && file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    # Grab the first contiguous comment/header region from the top
    # (stop when a non-comment, non-empty line is encountered after header begins)
    header <- character(0)
    started <- FALSE
    for (ln in lines) {
      if (!started) {
        if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
          header <- c(header, ln)
          if (grepl("^\\s*#", ln)) started <- TRUE
        } else {
          break
        }
      } else {
        if (grepl("^\\s*#", ln) || grepl("^\\s*$", ln)) {
          header <- c(header, ln)
        } else {
          break
        }
      }
    }
    return(paste(header, collapse = "\n"))
  }
  return(paste0(
    "# ======================================================================\n",
    "# Cross-layer State-Vector Similarity (CLSVS) Coupling Analysis\n",
    "# (Header unavailable because script path was not detected.)\n",
    "# ======================================================================"
  ))
}

header_text <- extract_header_block(script_full)

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))

qmd <- c(
  "---",
  paste0("title: \"", script_name, " Report\""),
  "format: html",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "toc: true",
  "---",
  "",
  "## Run metadata",
  "",
  paste0("- **run_id:** ", run_id),
  paste0("- **run_timestamp:** ", run_timestamp),
  paste0("- **script_name:** ", script_name),
  paste0("- **script_path:** ", script_path),
  paste0("- **script_full_path:** ", ifelse(is.na(script_full), "NA", script_full)),
  paste0("- **outputs_root:** ", normalizePath(outputs_root, winslash = "/", mustWork = FALSE)),
  paste0("- **output_dir:** ", normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
  "",
  "## Inputs",
  "",
  paste0("- **Metabolome CSV:** ", normalizePath(met_path, winslash = "/", mustWork = FALSE)),
  paste0("- **Transcriptome CSV:** ", normalizePath(tx_path, winslash = "/", mustWork = FALSE)),
  "",
  "## Script header",
  "",
  "```",
  header_text,
  "```",
  "",
  "## Analytical logic and formulas",
  "",
  "This report summarizes a coupling analysis between metabolome and transcriptome organizational dynamics.",
  "",
  "- **State vector per layer and day:**",
  "  - \\( s(d) = [ z(D(d)), z(-ED(d)) ] \\)",
  "  - \\(D\\) is a dispersion metric; \\(ED\\) is effective dimensionality; \\(z(\\cdot)\\) is a z-score computed across days *within a layer*.",
  "",
  "- **Coupling with lag \\(L\\):**",
  "  - \\( C(d;L) = \\frac{s_{met}(d) \\cdot s_{tx}(d+L)}{\\|s_{met}(d)\\|\\,\\|s_{tx}(d+L)\\|} \\in [-1,1] \\)",
  "",
  "Positive values indicate aligned organization regimes (after applying the lag). Negative values indicate opposing regimes.",
  "",
  "## Key outputs",
  "",
  "### Coupling by day",
  "",
  "```{r}",
  "knitr::include_graphics('CLSVS_Coupling_ByDay.png')",
  "```",
  "",
  "### Organization state-space by layer",
  "",
  "```{r}",
  "knitr::include_graphics('CLSVS_StateSpace_Scatter.png')",
  "```",
  "",
  "## Generated files inventory",
  "",
  "```{r}",
  "inv <- read.csv('Project_Manifest_Files.csv', stringsAsFactors = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "## Dependencies (packages + versions)",
  "",
  "```{r}",
  "deps <- jsonlite::fromJSON('Project_Manifest.json')$dependencies",
  "knitr::kable(as.data.frame(deps))",
  "```",
  "",
  "## Interpretation (facts and reasonable hypotheses only)",
  "",
  "- **Fact:** CLSVS provides a day-wise, lag-specific similarity score between metabolome and transcriptome organizational state vectors.",
  "- **Fact:** Peaks in CLSVS at a given lag indicate that the metabolome state at day \\(d\\) is most similar to the transcriptome state at day \\(d+L\\) under the defined state representation.",
  "- **Reasonable hypothesis:** If CLSVS is consistently higher for positive lags (e.g., 2–4 days) than for lag 0, this is consistent with temporal offset between layers; this alone does not establish causality.",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(qmd, qmd_path)

render_ok <- FALSE
render_msg <- NULL

tryCatch({
  out_html <- rmarkdown::render(
    input = qmd_path,
    output_file = basename(html_path),
    output_dir = output_dir,
    quiet = TRUE
  )
  render_ok <- TRUE
  render_msg <- paste0("HTML report created: ", normalizePath(out_html, winslash = "/", mustWork = FALSE))
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- paste0("HTML render failed. QMD written at: ", normalizePath(qmd_path, winslash = "/", mustWork = FALSE),
                        "\nError: ", conditionMessage(e))
})

# ---------------------------
# 9) Refresh file inventory after rendering (MANDATORY)
# ---------------------------
inv_post <- inventory_outputs(output_dir)
write_csv(inv_post, "Project_Manifest_Files.csv")

# Update manifest generated_files reference (optional embedded list)
manifest$generated_files <- inv_post$file
write_json(manifest, "Project_Manifest.json", pretty = TRUE, auto_unbox = TRUE, na = "null")

# ---------------------------
# 10) Console confirmation (MANDATORY)
# ---------------------------
message("\nRun complete.")
message("  run_id: ", run_id)
message("  output_dir: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message("  Manifest JSON: ", normalizePath(file.path(output_dir, "Project_Manifest.json"), winslash = "/", mustWork = FALSE))
message("  File inventory: ", normalizePath(file.path(output_dir, "Project_Manifest_Files.csv"), winslash = "/", mustWork = FALSE))
message(render_msg)
