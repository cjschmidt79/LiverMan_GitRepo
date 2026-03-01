#!/usr/bin/env Rscript
# ======================================================================
# Crosstalk Multi-Phase Statistics (Interactive; RStudio/source-friendly)
#
# PURPOSE
#   Applies 5 system-level analyses across ALL pathway functions and ALL
#   timeframes using existing Crosstalk output tables:
#     (1) Global temporal clustering of events (permutation test)
#     (2) Phase enrichment across functions (Fisher/chi-square-style contingency)
#     (3) Directional coherence across phases (Kruskal + permutation)
#     (4) Focality vs diffuseness across phases (Kruskal + permutation)
#     (5) Change-point alignment (argmax transition per pair; permutation test)
#
# INPUTS (user-selected interactively; prompted by REQUIRED filenames)
#   - Pair_SignBalance.csv
#   - Pair_Enrichment_sumAbsDelta.csv
#   - Pair_Enrichment_nEdges.csv
#   - Transition_Summary.csv
#
# OUTPUTS (MANDATORY POLICY)
#   All outputs are written only to:
#     outputs/<run_id>/
#   including:
#     - Project_Manifest.json
#     - Project_Manifest_Files.csv
#     - UNIFIED_pair_transition_table.csv
#     - ANALYSIS*.csv outputs + diagnostic plots
#     - <script_name>_Report.qmd and rendered HTML
#
# NOTE
#   This script is intended to be run interactively in RStudio via source(),
#   or via Rscript. In non-interactive mode, it will stop.
# ======================================================================

# -----------------------------
# Dependency handling (install if missing, quietly; never uninstall)
# -----------------------------
quiet_install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

min_pkgs <- c("jsonlite", "knitr", "rmarkdown")
analysis_pkgs <- c("dplyr", "tidyr", "stringr", "readr", "purrr", "ggplot2", "tibble")
quiet_install_if_missing(c(min_pkgs, analysis_pkgs))

suppressPackageStartupMessages({
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(tibble)
})

# -----------------------------
# Robust script identity capture (MANDATORY BLOCK; do not edit)
# -----------------------------
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
known_script_filename <- "Crosstalk_Stats_AllPhases.R"
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

# -----------------------------
# Utilities
# -----------------------------
stop2 <- function(...) stop(paste0(...), call. = FALSE)

timestamp_tag <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

normalize_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop2("File not found: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

parse_transition <- function(tr) {
  m <- stringr::str_match(tr, "^\\s*(\\d+)\\s*-\\s*(\\d+)\\s*$")
  if (any(is.na(m[,2])) || any(is.na(m[,3]))) {
    return(tibble(start_day = NA_integer_, end_day = NA_integer_, mid_day = NA_real_))
  }
  s <- as.integer(m[,2]); e <- as.integer(m[,3])
  tibble(start_day = s, end_day = e, mid_day = (s + e) / 2)
}

assign_phase <- function(start_day, end_day) {
  # Phase assignment by transition for your sampled days:
  # Day grid: 4, 8, 10, 12, 14, 16, 18, 20
  if (is.na(start_day) || is.na(end_day)) return(NA_character_)
  
  if (start_day == 4  && end_day == 8)  return("Permissive")
  if (start_day == 8  && end_day == 10) return("Tuning")
  if (start_day == 10 && end_day == 12) return("Tuning")
  if (start_day == 12 && end_day == 14) return("Constraint")
  if (start_day >= 14)                  return("Buffering")
  
  # Fallback based on midpoint (should rarely be used)
  mid <- (start_day + end_day) / 2
  if (mid <= 8) return("Permissive")
  if (mid <= 12) return("Tuning")
  if (mid <= 14) return("Constraint")
  "Buffering"
}

entropy_shannon <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log(p))
}

list_deps <- function(pkgs) {
  tibble(
    package = pkgs,
    version = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
  )
}

# -----------------------------
# Enforce interactive use (as requested)
# -----------------------------
if (!interactive()) {
  stop2(
    "This script is configured for interactive RStudio/source() usage.\n",
    "Run it interactively (e.g., source('", known_script_filename, "'))."
  )
}
# -----------------------------
# Execution mode guard: ensure script is run top-to-bottom
# -----------------------------
if (interactive()) {
  # If you are in RStudio, strongly encourage source() to avoid partial runs.
  # This is a warning (not fatal) because some users intentionally run in segments.
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    cat("[INFO] Recommended: run via source() to ensure inputs/manifest/report execute.\n")
  }
}

# -----------------------------
# Output directory policy (MANDATORY)
# -----------------------------
outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

analysis_name <- "CrosstalkStats"
source_label  <- "CrosstalkTables"
run_id <- paste0(analysis_name, "_", source_label, "_", timestamp_tag())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Console identity summary (MANDATORY)
# -----------------------------
cat("\n==================== Script Identity ====================\n")
cat("script_name: ", script_name, "\n", sep = "")
cat("script_path: ", script_path, "\n", sep = "")
cat("script_full: ", ifelse(is.na(script_full), "NA", script_full), "\n", sep = "")
if (is.na(script_full)) {
  cat("NOTE: script path detection failed; fallback filename used: ", known_script_filename, "\n", sep = "")
}
cat("outputs_root: ", normalize_path(outputs_root), "\n", sep = "")
cat("run_id:       ", run_id, "\n", sep = "")
cat("output_dir:   ", normalize_path(output_dir), "\n", sep = "")
cat("=========================================================\n\n")
cat("\n==================== Input Selection ====================\n")
cat("You will be prompted to select exactly these four files (by exact filename):\n")
cat("  1) Pair_SignBalance.csv\n")
cat("  2) Pair_Enrichment_sumAbsDelta.csv\n")
cat("  3) Pair_Enrichment_nEdges.csv\n")
cat("  4) Transition_Summary.csv\n")
cat("=========================================================\n\n")

# -----------------------------
# Interactive input selection (explicit prompts by REQUIRED filenames)
# -----------------------------
pick_required_file <- function(required_name) {
  repeat {
    cat("\nSelect required input file: ", required_name, "\n", sep = "")
    cat("A file chooser dialog will open. Please select: ", required_name, "\n", sep = "")
    p <- file.choose()
    p <- normalize_path(p)
    
    if (!file.exists(p)) {
      cat("  ERROR: File does not exist. Try again.\n")
      next
    }
    if (basename(p) != required_name) {
      cat("  ERROR: You selected: ", basename(p), "\n", sep = "")
      cat("         Expected:     ", required_name, "\n", sep = "")
      cat("         Please select the correct file.\n")
      next
    }
    return(p)
  }
}

input_sign_path   <- pick_required_file("Pair_SignBalance.csv")
input_sumabs_path <- pick_required_file("Pair_Enrichment_sumAbsDelta.csv")
input_nedges_path <- pick_required_file("Pair_Enrichment_nEdges.csv")
input_trans_path  <- pick_required_file("Transition_Summary.csv")

# -----------------------------
# Run parameters (editable defaults)
# -----------------------------
alpha <- 0.01        # event definition threshold for p_sumabs
B <- 5000            # permutations
min_events_per_function <- 5

# -----------------------------
# Load data
# -----------------------------
sign_df   <- safe_read_csv(input_sign_path)
sumabs_df <- safe_read_csv(input_sumabs_path)
nedges_df <- safe_read_csv(input_nedges_path)
trans_df  <- safe_read_csv(input_trans_path)

# -----------------------------
# HARD VALIDATION: inputs must exist and be data frames before analysis
# -----------------------------
required_objects <- c("sign_df", "sumabs_df", "nedges_df", "trans_df")
missing_objs <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objs) > 0) {
  stop2("Input tables missing before analysis: ", paste(missing_objs, collapse = ", "),
        "\nEnsure all four input files were selected and read successfully.")
}

# Verify they are data frames/tibbles and non-empty
is_df <- function(x) inherits(x, c("data.frame", "tbl_df", "tbl"))
tbls <- list(sign_df = sign_df, sumabs_df = sumabs_df, nedges_df = nedges_df, trans_df = trans_df)

bad_types <- names(tbls)[!vapply(tbls, is_df, logical(1))]
if (length(bad_types) > 0) stop2("These inputs are not data frames: ", paste(bad_types, collapse = ", "))

empty_tbls <- names(tbls)[vapply(tbls, nrow, integer(1)) == 0]
if (length(empty_tbls) > 0) stop2("These inputs have 0 rows: ", paste(empty_tbls, collapse = ", "))


# Required columns check (do not change scientific computations)
required_cols <- list(
  sign   = c("transition","pair","n_edges","n_pos","n_neg","frac_pos","frac_neg","sum_signed_delta","direction_bias"),
  sumabs = c("transition","pair","obs_sum_abs_delta","z","p_empirical_2sided"),
  nedges = c("transition","pair","obs_n_edges","z","p_empirical_2sided"),
  trans  = c("transition","n_edges","n_unique_pairs","entropy_nEdges","entropy_sumAbsDelta","top_pair_by_n","top_pair_by_abs")
)

check_cols <- function(df, need, label) {
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) stop2("Missing required columns in ", label, ": ", paste(miss, collapse = ", "))
}
check_cols(sign_df,   required_cols$sign,   "Pair_SignBalance.csv")
check_cols(sumabs_df, required_cols$sumabs, "Pair_Enrichment_sumAbsDelta.csv")
check_cols(nedges_df, required_cols$nedges, "Pair_Enrichment_nEdges.csv")
check_cols(trans_df,  required_cols$trans,  "Transition_Summary.csv")

# -----------------------------
# Build unified pair-transition table
# -----------------------------
pair_tbl <- sumabs_df %>%
  select(transition, pair, obs_sum_abs_delta, z_sumabs = z, p_sumabs = p_empirical_2sided) %>%
  inner_join(
    nedges_df %>% select(transition, pair, obs_n_edges, z_nedges = z, p_nedges = p_empirical_2sided),
    by = c("transition","pair")
  ) %>%
  inner_join(
    sign_df %>% select(transition, pair, n_edges, n_pos, n_neg, frac_pos, frac_neg, sum_signed_delta, direction_bias),
    by = c("transition","pair")
  ) %>%
  mutate(
    start_day = parse_transition(transition)$start_day,
    end_day   = parse_transition(transition)$end_day,
    mid_day   = parse_transition(transition)$mid_day,
    phase     = mapply(assign_phase, start_day, end_day),
    focality  = ifelse(obs_n_edges > 0, obs_sum_abs_delta / obs_n_edges, NA_real_),
    sign_coherence = abs(direction_bias),
    is_event  = (p_sumabs <= alpha)
  )

unified_path <- file.path(output_dir, "UNIFIED_pair_transition_table.csv")
write_csv(pair_tbl, unified_path)

# -----------------------------
# HARD VALIDATION: pair_tbl must exist before downstream analysis
# -----------------------------
if (!exists("pair_tbl")) {
  stop2("pair_tbl was not created. Upstream join or mutation failed.")
}
if (!inherits(pair_tbl, c("data.frame", "tbl_df", "tbl"))) {
  stop2("pair_tbl is not a data frame.")
}
if (nrow(pair_tbl) == 0) {
  stop2("pair_tbl has 0 rows. Check joins: transitions or pairs may not match across input tables.")
}
required_pair_tbl_cols <- c(
  "transition", "pair", "phase", "is_event",
  "obs_sum_abs_delta", "obs_n_edges",
  "direction_bias", "p_sumabs"
)
missing_cols <- setdiff(required_pair_tbl_cols, names(pair_tbl))
if (length(missing_cols) > 0) {
  stop2("pair_tbl missing required columns: ", paste(missing_cols, collapse = ", "))
}

events_path <- file.path(output_dir, paste0("EVENTS_alpha", gsub("\\.", "", as.character(alpha)), "_flagged.csv"))
write_csv(pair_tbl %>% select(-start_day,-end_day,-mid_day), events_path)

# -----------------------------
# Analysis 1: Global temporal clustering of events
# -----------------------------
set.seed(1)

event_counts <- pair_tbl %>%
  filter(is_event) %>%
  count(transition, name = "n_events") %>%
  right_join(pair_tbl %>% distinct(transition), by = "transition") %>%
  mutate(n_events = tidyr::replace_na(n_events, 0)) %>%
  arrange(transition)

obs_entropy1 <- entropy_shannon(event_counts$n_events)

transitions_all <- pair_tbl %>% distinct(transition) %>% pull(transition)

perm_entropy1 <- replicate(B, {
  if (sum(event_counts$n_events) == 0) return(NA_real_)
  ev_trans <- rep(event_counts$transition, event_counts$n_events)
  ev_trans_perm <- sample(transitions_all, size = length(ev_trans), replace = TRUE)
  perm_counts <- table(factor(ev_trans_perm, levels = transitions_all))
  entropy_shannon(as.numeric(perm_counts))
})

p_temporal_cluster <- mean(perm_entropy1 <= obs_entropy1, na.rm = TRUE)

analysis1 <- tibble(
  metric = "Event temporal clustering entropy (lower = more clustered)",
  alpha = alpha,
  obs_entropy = obs_entropy1,
  p_perm_lower = p_temporal_cluster,
  B = B
)
analysis1_path <- file.path(output_dir, "ANALYSIS1_temporal_clustering.csv")
write_csv(analysis1, analysis1_path)

# -----------------------------
# Analysis 2: Phase enrichment across functions (events only)
# -----------------------------
pair_funcs <- pair_tbl %>%
  mutate(
    funcA = str_split_fixed(pair, "__", 2)[, 1],
    funcB = str_split_fixed(pair, "__", 2)[, 2]
  ) %>%
  select(transition, phase, is_event, funcA, funcB) %>%
  pivot_longer(
    cols = c(funcA, funcB),
    names_to = "side",
    values_to = "function_name"
  ) %>%
  filter(!is.na(function_name) & nzchar(function_name))

func_phase_counts <- pair_funcs %>%
  filter(is_event) %>%
  count(function_name, phase, name = "n_event_mentions") %>%
  pivot_wider(names_from = phase, values_from = n_event_mentions, values_fill = 0)

analysis2_counts_path <- file.path(output_dir, "ANALYSIS2_function_by_phase_event_counts.csv")
write_csv(func_phase_counts, analysis2_counts_path)

eligible <- func_phase_counts %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  filter(total >= min_events_per_function)

phase_cols <- setdiff(names(eligible), c("function_name", "total"))
fisher_p <- NA_real_
if (length(phase_cols) >= 2 && nrow(eligible) >= 2) {
  mat <- as.matrix(eligible[, phase_cols, drop = FALSE])
  rownames(mat) <- eligible$function_name
  fisher_p <- suppressWarnings(fisher.test(mat, simulate.p.value = TRUE, B = 5000)$p.value)
}

analysis2 <- tibble(
  metric = "Function x Phase enrichment (events only; Fisher simulated p)",
  alpha = alpha,
  min_events_per_function = min_events_per_function,
  fisher_sim_p = fisher_p
)
analysis2_path <- file.path(output_dir, "ANALYSIS2_function_phase_fisher.csv")
write_csv(analysis2, analysis2_path)

# -----------------------------
# Analysis 3: Directional coherence across phases (system-wide)
# -----------------------------
set.seed(2)
coh_df <- pair_tbl %>% filter(!is.na(sign_coherence), !is.na(phase))
obs_kw3 <- suppressWarnings(kruskal.test(sign_coherence ~ phase, data = coh_df))
obs_kw3_p <- obs_kw3$p.value
obs_kw3_stat <- as.numeric(obs_kw3$statistic)

perm_kw3 <- replicate(B, {
  ph_perm <- sample(coh_df$phase)
  as.numeric(suppressWarnings(kruskal.test(coh_df$sign_coherence ~ ph_perm)$statistic))
})
p_coherence_perm <- mean(perm_kw3 >= obs_kw3_stat, na.rm = TRUE)

analysis3 <- tibble(
  metric = "Directional coherence differs by phase",
  kruskal_p = obs_kw3_p,
  obs_kruskal_stat = obs_kw3_stat,
  p_perm = p_coherence_perm,
  B = B
)
analysis3_path <- file.path(output_dir, "ANALYSIS3_sign_coherence_by_phase.csv")
write_csv(analysis3, analysis3_path)

# -----------------------------
# Analysis 4: Focality vs diffuseness across phases
# -----------------------------
set.seed(3)
foc_df <- pair_tbl %>% filter(is.finite(focality), !is.na(phase))
obs_kw4 <- suppressWarnings(kruskal.test(focality ~ phase, data = foc_df))
obs_kw4_p <- obs_kw4$p.value
obs_kw4_stat <- as.numeric(obs_kw4$statistic)

perm_kw4 <- replicate(B, {
  ph_perm <- sample(foc_df$phase)
  as.numeric(suppressWarnings(kruskal.test(foc_df$focality ~ ph_perm)$statistic))
})
p_focality_perm <- mean(perm_kw4 >= obs_kw4_stat, na.rm = TRUE)

analysis4 <- tibble(
  metric = "Focality differs by phase",
  kruskal_p = obs_kw4_p,
  obs_kruskal_stat = obs_kw4_stat,
  p_perm = p_focality_perm,
  B = B
)
analysis4_path <- file.path(output_dir, "ANALYSIS4_focality_by_phase.csv")
write_csv(analysis4, analysis4_path)

# -----------------------------
# Analysis 5: Change-point alignment (argmax transition per pair)
# -----------------------------
set.seed(4)
pair_argmax <- pair_tbl %>%
  group_by(pair) %>%
  slice_max(order_by = obs_sum_abs_delta, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(transition, name = "n_pair_peaks") %>%
  right_join(pair_tbl %>% distinct(transition), by = "transition") %>%
  mutate(n_pair_peaks = replace_na(n_pair_peaks, 0)) %>%
  arrange(transition)

obs_entropy5 <- entropy_shannon(pair_argmax$n_pair_peaks)

pair_split <- split(pair_tbl, pair_tbl$pair)

perm_entropy5 <- replicate(B, {
  peaks <- map_chr(pair_split, function(df) {
    tr_perm <- sample(df$transition)
    dfp <- df
    dfp$transition <- tr_perm
    dfp$transition[which.max(dfp$obs_sum_abs_delta)]
  })
  perm_counts <- table(factor(peaks, levels = transitions_all))
  entropy_shannon(as.numeric(perm_counts))
})

p_change_align <- mean(perm_entropy5 <= obs_entropy5, na.rm = TRUE)

analysis5 <- tibble(
  metric = "Change-point alignment (argmax transition entropy; lower = more clustered)",
  obs_entropy = obs_entropy5,
  p_perm_lower = p_change_align,
  B = B
)
analysis5_path <- file.path(output_dir, "ANALYSIS5_changepoint_alignment.csv")
write_csv(analysis5, analysis5_path)

# -----------------------------
# Diagnostic plots (always-on: at least one plot)
# -----------------------------
plot1_path <- file.path(output_dir, "PLOT_events_by_transition.png")
p1 <- event_counts %>%
  ggplot(aes(x = transition, y = n_events)) +
  geom_col() +
  theme_bw() +
  labs(title = paste0("Event counts by transition (alpha=", alpha, ")"),
       x = "Transition", y = "Number of significant events")
ggsave(plot1_path, p1, width = 8, height = 4, dpi = 200)

plot2_path <- file.path(output_dir, "PLOT_sign_coherence_by_phase.png")
p2 <- coh_df %>%
  ggplot(aes(x = phase, y = sign_coherence)) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_bw() +
  labs(title = "Directional coherence by phase", x = "Phase", y = "|direction_bias|")
ggsave(plot2_path, p2, width = 6, height = 4, dpi = 200)

plot3_path <- file.path(output_dir, "PLOT_focality_by_phase.png")
p3 <- foc_df %>%
  ggplot(aes(x = phase, y = focality)) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_bw() +
  labs(title = "Focality by phase", x = "Phase", y = "obs_sum_abs_delta / obs_n_edges")
ggsave(plot3_path, p3, width = 6, height = 4, dpi = 200)

# -----------------------------
# Manifest (machine-readable + human-readable inventory)
# -----------------------------
deps_tbl <- list_deps(c(min_pkgs, analysis_pkgs))

write_manifest <- function(generated_files_tbl) {
  manifest <- list(
    run_id = run_id,
    run_timestamp = as.character(Sys.time()),
    script = list(
      name = script_name,
      path = script_path,
      full_path = ifelse(is.na(script_full), NA, script_full)
    ),
    input = list(
      Pair_SignBalance = input_sign_path,
      Pair_Enrichment_sumAbsDelta = input_sumabs_path,
      Pair_Enrichment_nEdges = input_nedges_path,
      Transition_Summary = input_trans_path
    ),
    parameters = list(
      alpha = alpha,
      B = B,
      min_events_per_function = min_events_per_function,
      phase_definition = list(
        Permissive = "4-8",
        Tuning = "8-10 and 10-12",
        Constraint = "12-14",
        Buffering = ">=14 (14-16,16-18,18-20)"
      )
    ),
    dependencies = split(deps_tbl, seq_len(nrow(deps_tbl))) |> lapply(as.list),
    outputs = list(
      outputs_root = normalize_path(outputs_root),
      output_dir = normalize_path(output_dir)
    ),
    generated_files = generated_files_tbl
  )
  
  json_path <- file.path(output_dir, "Project_Manifest.json")
  jsonlite::write_json(manifest, json_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  json_path
}

inventory_files <- function() {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  tibble(
    file = normalizePath(files, winslash = "/", mustWork = FALSE),
    rel_path = sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(files, winslash = "/", mustWork = FALSE)),
    bytes = file.info(files)$size,
    modified = as.character(file.info(files)$mtime)
  ) %>% arrange(rel_path)
}

# Initial inventory (pre-report)
inv0 <- inventory_files()
inv_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")
write_csv(inv0, inv_csv_path)

# Write manifest (pre-report)
manifest_json_path <- write_manifest(split(inv0, seq_len(nrow(inv0))) |> lapply(as.list))

# -----------------------------
# Build Quarto QMD (MANDATORY) + render to HTML (final step)
# -----------------------------
# Capture script header block (initial comment block) for embedding
get_header_block <- function() {
  # If script_full is available, read from it; otherwise embed a placeholder.
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
    # Capture from top until first non-comment, non-empty line after header block.
    # We'll take first ~120 lines or until a line that doesn't start with "#".
    take_n <- min(length(lines), 160)
    head_lines <- lines[seq_len(take_n)]
    # keep consecutive comment/blank lines from top
    keep <- c()
    started <- FALSE
    for (ln in head_lines) {
      if (!started) {
        if (grepl("^\\s*#|^\\s*$", ln)) {
          keep <- c(keep, ln)
          started <- TRUE
        } else {
          # no header comment at top
          break
        }
      } else {
        if (grepl("^\\s*#|^\\s*$", ln)) {
          keep <- c(keep, ln)
        } else {
          break
        }
      }
    }
    paste(keep, collapse = "\n")
  } else {
    paste0("# Header block unavailable (script_full not detected). Fallback: ", known_script_filename)
  }
}

header_text <- get_header_block()

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
html_path <- file.path(output_dir, paste0(script_name, "_Report.html"))

qmd_text <- c(
  "---",
  "---",
  paste0("title: \"", script_name, " Report\""),
  "format:",
  "  html:",
  "    toc: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "params:",
  "  run_id: null",
  "  run_timestamp: null",
  "  script_name: null",
  "  script_path: null",
  "  script_full: null",
  "  in_sign: null",
  "  in_sumabs: null",
  "  in_nedges: null",
  "  in_trans: null",
  "  in_selected: null",
  "  outputs_root: null",
  "  output_dir: null",
  "  alpha: null",
  "---",
  
  "",
  "# Script header",
  "```",
  header_text,
  "```",
  "",
  "# Run metadata",
  "",
  "**Run ID:** `r params$run_id`  ",
  "**Timestamp:** `r params$run_timestamp`  ",
  "",
  "## Script identity",
  "- **script_name:** `r params$script_name`",
  "- **script_path:** `r params$script_path`",
  "- **script_full:** `r params$script_full`",
  "",
  "## Inputs",
  "- Pair_SignBalance.csv: `r params$in_sign`",
  "- Pair_Enrichment_sumAbsDelta.csv: `r params$in_sumabs`",
  "- Pair_Enrichment_nEdges.csv: `r params$in_nedges`",
  "- Transition_Summary.csv: `r params$in_trans`",
  "",
  "## Outputs",
  "- outputs_root: `r params$outputs_root`",
  "- output_dir: `r params$output_dir`",
  "",
  "# Analytical logic and formulas",
  "",
  "## Event definition",
  "We define an event for a pathway-pair transition when:",
  "",
  "$$",
  "\\mathrm{event} \\iff p_{\\mathrm{sumAbsDelta}} \\le \\alpha",
  "$$",
  "",
  "with $\\alpha = `r params$alpha`.",
  "",
  "## Phase definitions (by transition)",
  "- **Permissive:** 4–8",
  "- **Tuning:** 8–10 and 10–12",
  "- **Constraint:** 12–14",
  "- **Buffering:** 14–16, 16–18, 18–20",
  "",
  "## Analysis overview",
  "1. **Temporal clustering** of events across transitions using Shannon entropy and permutation tests.",
  "2. **Function × Phase enrichment** of event mentions using contingency tables and Fisher’s exact test (simulated p).",
  "3. **Directional coherence** across phases using $|\\mathrm{direction\\_bias}|$ and Kruskal–Wallis + permutation.",
  "4. **Focality** across phases using $\\mathrm{sumAbsDelta}/\\mathrm{nEdges}$ and Kruskal–Wallis + permutation.",
  "5. **Change-point alignment** by argmax transition per pair and entropy-based permutation tests.",
  "",
  "# Key results (tables)",
  "",
  "## Analysis 1: Temporal clustering",
  "```{r}",
  "a1 <- readr::read_csv(file.path(params$output_dir, 'ANALYSIS1_temporal_clustering.csv'), show_col_types = FALSE)",
  "knitr::kable(a1)",
  "```",
  "",
  "## Analysis 2: Function × Phase enrichment",
  "```{r}",
  "a2 <- readr::read_csv(file.path(params$output_dir, 'ANALYSIS2_function_phase_fisher.csv'), show_col_types = FALSE)",
  "knitr::kable(a2)",
  "```",
  "",
  "## Analysis 3: Directional coherence",
  "```{r}",
  "a3 <- readr::read_csv(file.path(params$output_dir, 'ANALYSIS3_sign_coherence_by_phase.csv'), show_col_types = FALSE)",
  "knitr::kable(a3)",
  "```",
  "",
  "## Analysis 4: Focality",
  "```{r}",
  "a4 <- readr::read_csv(file.path(params$output_dir, 'ANALYSIS4_focality_by_phase.csv'), show_col_types = FALSE)",
  "knitr::kable(a4)",
  "```",
  "",
  "## Analysis 5: Change-point alignment",
  "```{r}",
  "a5 <- readr::read_csv(file.path(params$output_dir, 'ANALYSIS5_changepoint_alignment.csv'), show_col_types = FALSE)",
  "knitr::kable(a5)",
  "```",
  "",
  "# Plots",
  "",
  "## Event counts by transition (always-on)",
  "```{r}",
  "knitr::include_graphics(file.path(params$output_dir, 'PLOT_events_by_transition.png'))",
  "```",
  "",
  "## Directional coherence by phase",
  "```{r}",
  "knitr::include_graphics(file.path(params$output_dir, 'PLOT_sign_coherence_by_phase.png'))",
  "```",
  "",
  "## Focality by phase",
  "```{r}",
  "knitr::include_graphics(file.path(params$output_dir, 'PLOT_focality_by_phase.png'))",
  "```",
  "",
  "# Generated files",
  "```{r}",
  "inv <- readr::read_csv(file.path(params$output_dir, 'Project_Manifest_Files.csv'), show_col_types = FALSE)",
  "knitr::kable(inv)",
  "```",
  "",
  "# Dependencies",
  "```{r}",
  "deps <- readr::read_csv(file.path(params$output_dir, 'Dependencies.csv'), show_col_types = FALSE)",
  "knitr::kable(deps)",
  "```",
  "",
  "# Interpretation guide",
  "",
  "## How to read the statistics",
  "- **Lower entropy (Analyses 1 & 5)** indicates stronger temporal clustering (more discrete phase behavior).",
  "- **Directional coherence (Analysis 3)** increases when interaction signs become more one-sided across pathway pairs, consistent with constrained coordination.",
  "- **Higher focality (Analysis 4)** indicates concentrated reorganization (large change per edge) rather than diffuse remodeling.",
  "",
  "## Biological interpretation (high-level)",
  "Evidence for discrete, phase-structured reorganization supports a staged control model of development in which the system shifts from permissive organization (early), to tuning (mid), to a constrained regime (constraint landmark), followed by buffering/maintenance.",
  "",
  "# Reproducibility",
  "```{r}",
  "sessionInfo()",
  "```"
)

# Write deps CSV (for QMD to read)
deps_csv_path <- file.path(output_dir, "Dependencies.csv")
write_csv(deps_tbl, deps_csv_path)

writeLines(qmd_text, qmd_path)

# Render as final step (MANDATORY)
render_ok <- TRUE
render_msg <- NULL
tryCatch({
  # -------------------------------------------------------------
  # FIX: Prepare input file metadata BEFORE rendering (IMPORTANT)
  # -------------------------------------------------------------
  input_files <- tibble::tibble(
    input_label = c("Pair_SignBalance", "Pair_Enrichment_sumAbsDelta", "Pair_Enrichment_nEdges", "Transition_Summary"),
    required_filename = c("Pair_SignBalance.csv", "Pair_Enrichment_sumAbsDelta.csv", "Pair_Enrichment_nEdges.csv", "Transition_Summary.csv"),
    selected_path = c(input_sign_path, input_sumabs_path, input_nedges_path, input_trans_path),
    selected_filename = basename(c(input_sign_path, input_sumabs_path, input_nedges_path, input_trans_path))
  )
  
  # Save the table to CSV
  write_csv(input_files, file.path(output_dir, "Input_Files_Selected.csv"))
  
  # -------------------------------------------------------------
  # Render the Quarto report with valid params list
  # -------------------------------------------------------------
  render_ok <- TRUE
  render_msg <- NULL
  
  tryCatch({
    rmarkdown::render(
      input = qmd_path,
      output_format = "html_document",
      output_file = basename(html_path),
      output_dir = output_dir,
      quiet = TRUE,
      params = list(
        run_id = run_id,
        run_timestamp = as.character(Sys.time()),
        script_name = script_name,
        script_path = script_path,
        script_full = ifelse(is.na(script_full), NA, script_full),
        in_sign = input_sign_path,
        in_sumabs = input_sumabs_path,
        in_nedges = input_nedges_path,
        in_trans = input_trans_path,
        outputs_root = normalize_path(outputs_root),
        output_dir = normalize_path(output_dir),
        alpha = alpha
      )
    )
  }, error = function(e) {
    render_ok <<- FALSE
    render_msg <<- conditionMessage(e)
  })
  
}, error = function(e) {
  render_ok <<- FALSE
  render_msg <<- conditionMessage(e)
})

# Refresh inventory AFTER rendering (MANDATORY)
inv_final <- inventory_files()
write_csv(inv_final, inv_csv_path)

# Refresh manifest AFTER rendering (MANDATORY)
manifest_json_path <- write_manifest(split(inv_final, seq_len(nrow(inv_final))) |> lapply(as.list))

# Console confirmation (MANDATORY)
if (render_ok && file.exists(html_path)) {
  cat("\nHTML report created: ", normalize_path(html_path), "\n", sep = "")
} else {
  cat("\nERROR: HTML report render failed.\n")
  if (!is.null(render_msg)) cat("Reason: ", render_msg, "\n", sep = "")
  cat("QMD written at: ", normalize_path(qmd_path), "\n", sep = "")
}

cat("\nRun folder: ", normalize_path(output_dir), "\n", sep = "")
cat("Manifest JSON: ", normalize_path(manifest_json_path), "\n", sep = "")
cat("File inventory CSV: ", normalize_path(inv_csv_path), "\n\n", sep = "")
