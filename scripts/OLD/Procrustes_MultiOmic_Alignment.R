#!/usr/bin/env Rscript
# ==============================================================================
# Multi-Omic Procrustes: Phase-Specific Coordination & Stability (Omic-aware)
#
# INPUT REQUIREMENTS (Tidy / Long):
#   Columns required:
#     - Day (numeric/integer)
#     - SampleID (character)
#     - Feature (character)   # formerly "Gene" in older tables
#     - Abundance (numeric)
#     - Omic (character)      # must contain at least Transcriptome and Metabolome rows
#
# PHASES:
#   A: Days 4,6,8
#   B: Days 10,12,14
#   C: Days 16,18,20
#
# OUTPUTS:
#   - Missing_Data_Inventory.csv (within-omic full grid audit)
#   - Phase_Coordination_Summary.csv
#   - Procrustes_Plotting_Data.csv
#   - Stability_Trend.png
#   - Trajectory_Map.png
#   - <script>_Report.html (Quarto-like HTML via rmarkdown render of .Rmd written by script)
#
# NOTES:
#   - Procrustes requires shared SampleIDs across Omics within each subset (global / phase).
#   - Script writes to: ./outputs/<run_id>/
# ==============================================================================

# ==============================================================================
# 1. SCRIPT IDENTITY & PATH RESOLUTION
# ==============================================================================
resolve_script_path <- function() {
  p <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      rstudioapi::getSourceEditorContext()$path
    } else ""
  }, error = function(e) "")
  if (nzchar(p) && file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    p2 <- sub("^--file=", "", file_arg[1])
    if (nzchar(p2) && file.exists(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  
  p3 <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(p3) && nzchar(p3) && file.exists(p3)) return(normalizePath(p3, winslash = "/", mustWork = FALSE))
  
  NA_character_
}

known_script_filename <- "Procrustes_Phase_Analysis_OmicAware.R"
script_full <- resolve_script_path()
script_name <- if (!is.na(script_full)) tools::file_path_sans_ext(basename(script_full)) else tools::file_path_sans_ext(known_script_filename)

# ==============================================================================
# 2. PACKAGE SETUP
# ==============================================================================
stop2 <- function(...) stop(paste0(...), call. = FALSE)

needed <- c("tidyverse", "vegan", "jsonlite", "rmarkdown", "knitr", "grid")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop2("Missing packages: ", paste(missing, collapse = ", "),
                               "\nInstall them, e.g.: install.packages(c(", paste0("'", missing, "'", collapse = ", "), "))")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(jsonlite)
  library(rmarkdown)
  library(knitr)
  library(grid)
})

# ==============================================================================
# 3. INTERACTIVE INPUT & OUTPUT DIRECTORIES
# ==============================================================================
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  input_file <- rstudioapi::selectFile(caption = "Select Tidy CSV (requires Omic column)", filter = "CSV Files (*.csv)")
} else {
  input_file <- file.choose()
}

if (is.null(input_file) || input_file == "") stop2("No file selected.")
input_file_abs <- normalizePath(input_file, winslash = "/", mustWork = FALSE)

run_id <- paste0("Procrustes_Analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
output_dir <- file.path(getwd(), "outputs", run_id)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 4. DATA LOAD + VALIDATION + STANDARDIZATION
# ==============================================================================
df <- read.csv(input_file_abs, stringsAsFactors = FALSE, check.names = FALSE)

# Backward compat: if "Gene" exists and "Feature" does not, rename it.
if (!"Feature" %in% names(df) && "Gene" %in% names(df)) {
  df <- df %>% rename(Feature = Gene)
}

required_cols <- c("Day", "SampleID", "Feature", "Abundance", "Omic")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop2("Input file is missing required columns: ", paste(missing_cols, collapse = ", "),
        "\nRequired columns: ", paste(required_cols, collapse = ", "))
}

# Basic types
df <- df %>%
  mutate(
    Day = as.integer(Day),
    SampleID = as.character(SampleID),
    Feature = as.character(Feature),
    Abundance = suppressWarnings(as.numeric(Abundance)),
    Omic = as.character(Omic)
  )

if (any(is.na(df$Abundance))) {
  n_bad <- sum(is.na(df$Abundance))
  warning("Found ", n_bad, " NA Abundance values after numeric coercion. These will propagate to matrices unless handled.")
}

# Phase assignment
df <- df %>%
  mutate(
    Phase = dplyr::case_when(
      Day %in% c(4, 6, 8) ~ "A",
      Day %in% c(10, 12, 14) ~ "B",
      Day %in% c(16, 18, 20) ~ "C",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(df$Phase))) {
  bad_days <- sort(unique(df$Day[is.na(df$Phase)]))
  warning("Some days are not in Phase A/B/C definition and were assigned Phase=NA: ",
          paste(bad_days, collapse = ", "),
          "\nThey will be included in GLOBAL analysis but excluded from PHASE-SPECIFIC splits (A/B/C).")
}

# Validate omics present
omics_present <- sort(unique(df$Omic))
if (!("Transcriptome" %in% omics_present) || !("Metabolome" %in% omics_present)) {
  stop2("Input must contain Omic rows for BOTH Transcriptome and Metabolome.\n",
        "Omics present: ", paste(omics_present, collapse = ", "), "\n",
        "Ensure Omic column includes exactly (or at least) values 'Transcriptome' and 'Metabolome'.")
}

# ==============================================================================
# 5. MISSING DATA AUDIT (WITHIN-OMIC FULL GRID)
#    This checks whether every SampleID has every Feature *within each Omic*.
# ==============================================================================
audit_missing_within_omic <- function(dat) {
  dat %>%
    distinct(Omic, SampleID, Feature) %>%
    group_by(Omic) %>%
    group_modify(~{
      full_grid <- tidyr::expand_grid(
        SampleID = unique(.x$SampleID),
        Feature  = unique(.x$Feature)
      )
      missing <- full_grid %>%
        anti_join(.x, by = c("SampleID", "Feature")) %>%
        mutate(Omic = unique(.x$Omic))
      missing
    }) %>%
    ungroup()
}

missing_data <- audit_missing_within_omic(df)
write.csv(missing_data, file.path(output_dir, "Missing_Data_Inventory.csv"), row.names = FALSE)

# ==============================================================================
# 6. PROCRUSTES CORE
# ==============================================================================
run_procrustes <- function(sub_df, permutations = 999, k = 2) {
  
  # Split by Omic explicitly
  tx_sub <- sub_df %>% filter(Omic == "Transcriptome")
  me_sub <- sub_df %>% filter(Omic == "Metabolome")
  
  if (nrow(tx_sub) == 0) stop2("No Transcriptome rows in this subset. Procrustes requires both layers.")
  if (nrow(me_sub) == 0) stop2("No Metabolome rows in this subset. Procrustes requires both layers.")
  
  tx_wide <- tx_sub %>%
    select(SampleID, Day, Phase, Feature, Abundance) %>%
    tidyr::pivot_wider(names_from = Feature, values_from = Abundance, values_fill = 0)
  
  me_wide <- me_sub %>%
    select(SampleID, Day, Phase, Feature, Abundance) %>%
    tidyr::pivot_wider(names_from = Feature, values_from = Abundance, values_fill = 0)
  
  common <- intersect(tx_wide$SampleID, me_wide$SampleID)
  if (length(common) < 3) {
    stop2("Too few shared SampleID values between layers for Procrustes (need >= 3). Shared = ", length(common))
  }
  
  tx_wide_c <- tx_wide %>% filter(SampleID %in% common) %>% arrange(SampleID)
  me_wide_c <- me_wide %>% filter(SampleID %in% common) %>% arrange(SampleID)
  
  # Matrices
  t_mat <- tx_wide_c %>% select(-SampleID, -Day, -Phase) %>% as.matrix()
  m_mat <- me_wide_c %>% select(-SampleID, -Day, -Phase) %>% as.matrix()
  
  if (nrow(t_mat) < 3 || nrow(m_mat) < 3) stop2("Need >= 3 samples for ordination/procrustes.")
  if (ncol(t_mat) < 2) stop2("Transcriptome matrix has <2 features after widening; cannot ordinate.")
  if (ncol(m_mat) < 2) stop2("Metabolome matrix has <2 features after widening; cannot ordinate.")
  
  # Ordination (MDS on Euclidean distances of log1p)
  t_mds <- cmdscale(vegdist(log1p(t_mat), "euclidean"), k = k)
  m_mds <- cmdscale(vegdist(log1p(m_mat), "euclidean"), k = k)
  
  # Procrustes + permutation test
  test <- protest(t_mds, m_mds, permutations = permutations)
  fit  <- procrustes(t_mds, m_mds, symmetric = TRUE)
  
  # Day labels derived from transcriptome wide (they should match metabolome, but transcriptome is fine)
  days <- tx_wide_c$Day
  phases <- tx_wide_c$Phase
  
  list(
    r = unname(test$t0),
    p = unname(test$signif),
    residuals = residuals(fit),
    tx = fit$X,
    mx = fit$Yrot,
    samples = tx_wide_c$SampleID,
    days = days,
    phases = phases,
    n_tx_features = ncol(t_mat),
    n_me_features = ncol(m_mat),
    n_samples = length(common)
  )
}

# ==============================================================================
# 7. GLOBAL + PHASE-SPECIFIC RUNS
# ==============================================================================
global_res <- run_procrustes(df)

phase_summary <- df %>%
  filter(!is.na(Phase)) %>%
  group_by(Phase) %>%
  group_split() %>%
  purrr::map_dfr(~{
    res <- run_procrustes(.x)
    data.frame(
      Phase = unique(na.omit(.x$Phase))[1],
      Correlation_R = res$r,
      P_Value = res$p,
      N_Samples = res$n_samples,
      N_TX_Features = res$n_tx_features,
      N_ME_Features = res$n_me_features
    )
  })

write.csv(phase_summary, file.path(output_dir, "Phase_Coordination_Summary.csv"), row.names = FALSE)

# ==============================================================================
# 8. EXPORT: PLOTTING TABLE
# ==============================================================================
sample_data <- data.frame(
  SampleID = global_res$samples,
  Day = global_res$days,
  Phase = global_res$phases,
  Residual = global_res$residuals
)

stability_stats <- sample_data %>%
  group_by(Day) %>%
  summarise(
    Mean_Residual = mean(Residual, na.rm = TRUE),
    SEM_Residual = sd(Residual, na.rm = TRUE) / sqrt(sum(!is.na(Residual))),
    .groups = "drop"
  )

centroid_data <- data.frame(
  Day = global_res$days,
  TX_MDS1 = global_res$tx[, 1],
  TX_MDS2 = global_res$tx[, 2],
  ME_MDS1 = global_res$mx[, 1],
  ME_MDS2 = global_res$mx[, 2]
) %>%
  group_by(Day) %>%
  summarise(
    TX_MDS1 = mean(TX_MDS1), TX_MDS2 = mean(TX_MDS2),
    ME_MDS1 = mean(ME_MDS1), ME_MDS2 = mean(ME_MDS2),
    .groups = "drop"
  ) %>%
  mutate(Vector_Length = sqrt((TX_MDS1 - ME_MDS1)^2 + (TX_MDS2 - ME_MDS2)^2))

plotting_table <- left_join(stability_stats, centroid_data, by = "Day")
write.csv(plotting_table, file.path(output_dir, "Procrustes_Plotting_Data.csv"), row.names = FALSE)

cat("\nPlotting table exported to:", file.path(output_dir, "Procrustes_Plotting_Data.csv"), "\n")

# ==============================================================================
# 9. VISUALIZATION OUTPUTS
# ==============================================================================
# Stability bar chart
plot_data <- data.frame(
  SampleID = global_res$samples,
  Day = as.factor(global_res$days),
  residual = global_res$residuals
)

summary_stats <- plot_data %>%
  group_by(Day) %>%
  summarise(
    mean_res = mean(residual, na.rm = TRUE),
    se_res = sd(residual, na.rm = TRUE) / sqrt(sum(!is.na(residual))),
    .groups = "drop"
  )

p_bar <- ggplot(summary_stats, aes(x = Day, y = mean_res, fill = Day)) +
  geom_col(color = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_res - se_res, ymax = mean_res + se_res), width = 0.2) +
  theme_minimal() +
  labs(
    title = "Multi-Omic Decoupling (Stability Trend)",
    y = "Average Procrustes Residual (Higher = Less Coordinated)",
    x = "Day"
  )

# Trajectory map (centroids)
centroids_for_plot <- data.frame(
  Day = as.factor(global_res$days),
  tx1 = global_res$tx[, 1], tx2 = global_res$tx[, 2],
  mx1 = global_res$mx[, 1], mx2 = global_res$mx[, 2]
) %>%
  group_by(Day) %>%
  summarise(
    tx1 = mean(tx1), tx2 = mean(tx2),
    mx1 = mean(mx1), mx2 = mean(mx2),
    .groups = "drop"
  )

p_traj <- ggplot(centroids_for_plot) +
  geom_segment(
    aes(x = tx1, y = tx2, xend = mx1, yend = mx2, color = Day),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 1.2
  ) +
  geom_point(aes(x = tx1, y = tx2, fill = Day), shape = 21, size = 5) +
  theme_minimal() +
  labs(
    title = "Developmental Trajectory: Transcriptome to Metabolome",
    x = "MDS1",
    y = "MDS2"
  )

ggsave(file.path(output_dir, "Stability_Trend.png"), p_bar, width = 7, height = 5)
ggsave(file.path(output_dir, "Trajectory_Map.png"), p_traj, width = 7, height = 5)

# ==============================================================================
# 10. REPORT (Rmarkdown .Rmd written + rendered to HTML)
#     You asked for: "provide 30 rows of all csv tables in the html output".
#     We therefore:
#       - Read every .csv in output_dir
#       - Print head(30) via knitr::kable
# ==============================================================================
rmd_path <- file.path(output_dir, paste0(script_name, "_Report.Rmd"))
html_out <- file.path(output_dir, paste0(script_name, "_Report.html"))

rmd_content <- c(
  "---",
  "title: \"Phase-Specific Procrustes Analysis (Omic-aware)\"",
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_depth: 3",
  "    number_sections: true",
  "params:",
  "  output_dir: NULL",
  "  input_file: NULL",
  "  run_id: NULL",
  "  run_timestamp: NULL",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "```",
  "",
  "## Run Summary",
  "",
  "```{r}",
  "cat('Run ID: ', params$run_id, '\\n')",
  "cat('Run Timestamp: ', params$run_timestamp, '\\n')",
  "cat('Input File: ', params$input_file, '\\n')",
  "cat('Output Directory: ', params$output_dir, '\\n')",
  "```",
  "",
  "## Phase-Specific Coordination",
  "",
  "```{r}",
  "phase_path <- file.path(params$output_dir, 'Phase_Coordination_Summary.csv')",
  "if (file.exists(phase_path)) {",
  "  knitr::kable(utils::head(read.csv(phase_path), 30), caption = 'Phase_Coordination_Summary.csv (first 30 rows)')",
  "} else {",
  "  cat('Phase_Coordination_Summary.csv not found.')",
  "}",
  "```",
  "",
  "## Multi-Omic Stability",
  "",
  "```{r}",
  "img1 <- file.path(params$output_dir, 'Stability_Trend.png')",
  "if (file.exists(img1)) knitr::include_graphics(img1) else cat('Stability_Trend.png not found.')",
  "```",
  "",
  "## Developmental Trajectory",
  "",
  "```{r}",
  "img2 <- file.path(params$output_dir, 'Trajectory_Map.png')",
  "if (file.exists(img2)) knitr::include_graphics(img2) else cat('Trajectory_Map.png not found.')",
  "```",
  "",
  "## All Output CSV Tables (first 30 rows each)",
  "",
  "```{r}",
  "csvs <- list.files(params$output_dir, pattern = '\\\\.[cC][sS][vV]$', full.names = TRUE)",
  "if (length(csvs) == 0) {",
  "  cat('No CSV files found in output directory.')",
  "} else {",
  "  for (f in csvs) {",
  "    cat('\\n\\n### ', basename(f), '\\n')",
  "    tab <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)",
  "    if (is.null(tab)) {",
  "      cat('Could not read: ', basename(f))",
  "    } else {",
  "      knitr::kable(utils::head(tab, 30))",
  "    }",
  "  }",
  "}",
  "```",
  "",
  "## Reproducibility",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

writeLines(rmd_content, rmd_path)

# Render with knit root = output_dir so relative file reads are stable
tryCatch({
  rmarkdown::render(
    input = rmd_path,
    output_file = basename(html_out),
    output_dir = output_dir,
    quiet = TRUE,
    params = list(
      output_dir = output_dir,
      input_file = input_file_abs,
      run_id = run_id,
      run_timestamp = as.character(Sys.time())
    ),
    envir = new.env(parent = globalenv())
  )
}, error = function(e) {
  stop2("HTML report render failed.\nReason: ", conditionMessage(e), "\nRmd written at: ", rmd_path)
})

cat("\nAnalysis Complete. HTML Report created at:\n", html_out, "\n")
cat("\nOutputs directory:\n", output_dir, "\n")
