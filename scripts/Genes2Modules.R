#!/usr/bin/env Rscript
# ======================================================================
# Constraint Modules Over Time (Genes + Pathway-Event Network)
# SOURCE-to-run / Biologist-first (RStudio friendly) — NO dplyr/tidyverse
#
# WHAT THIS DOES
#   A) Gene-level "constraint-signature modules"
#      - Groups genes by which sampled days they are constrained (binary signature).
#      - Outputs module membership + summary tables.
#
#   B) Pathway/function interaction modules per transition
#      - Uses UNIFIED_pair_transition_table.csv (pathway-pair × transition)
#      - Defines event edges as p_sumabs <= alpha (or is_event == TRUE)
#      - Builds a graph per transition (nodes = functions, edges = TRUE events)
#      - Community detection per transition (Louvain)
#      - Outputs module membership tables and module stability (Jaccard overlaps)
#
#   C) Optional projection: genes → function-modules by transition
#      - If you provide a Gene↔Function mapping table, the script will project
#        gene membership into each transition’s function-modules.
#
# INPUTS (interactive pickers)
#   1) Gene constraint table CSV  (must contain gene + days_constrained)
#   2) UNIFIED_pair_transition_table.csv
#   3) (Optional) Gene↔Function mapping CSV (gene + function_name)
#
# OUTPUTS (MANDATORY POLICY)
#   Writes ONLY to: ./outputs/<run_name>_<timestamp>/
#   Includes:
#     - Project_Manifest.json
#     - Project_Manifest_Files.csv
#     - Dependencies.csv
#     - Tables in ./tables/
#     - Plots in ./plots/ (600 dpi PNG)
#     - Quarto report: <script_name>_Report.qmd and HTML
#
# NOTE
#   Designed for RStudio: open this script and press SOURCE.
# ======================================================================

# -----------------------------
# Dependency handling (install if missing; never uninstall)
# -----------------------------
quiet_install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}

min_pkgs <- c("jsonlite", "data.table", "knitr")
analysis_pkgs <- c("igraph")

quiet_install_if_missing(c(min_pkgs, analysis_pkgs))

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(igraph)
})

# -----------------------------
# Robust script identity capture (MANDATORY BLOCK; do not edit)
# -----------------------------
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

known_script_filename <- "Constraint_Modules_OverTime.R"
known_script_stem <- tools::file_path_sans_ext(known_script_filename)

script_full <- resolve_script_path()
if (is.na(script_full)) {
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
say <- function(...) cat(paste0(...), "\n")

safe_fread <- function(path) {
  if (!file.exists(path)) stop2("File not found: ", path)
  data.table::fread(path)
}

pick_file <- function(prompt_text) {
  p <- ""
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::selectFile(caption = prompt_text), error = function(e) "")
  }
  if (!nzchar(p)) {
    say(prompt_text)
    p <- file.choose()
  }
  normalize_path(p)
}

pick_yes_no <- function(question) {
  ans <- readline(paste0(question, " [y/n]: "))
  ans <- tolower(trimws(ans))
  ans %in% c("y", "yes")
}

pick_column <- function(dt, prompt_text, default_guess = NULL) {
  cn <- names(dt)
  say("\n", prompt_text)
  for (i in seq_along(cn)) say(sprintf("  %d) %s", i, cn[i]))
  if (!is.null(default_guess) && default_guess %in% cn) {
    say(sprintf("Press ENTER for default: %s", default_guess))
  }
  idx <- readline("Select a column number: ")
  idx <- trimws(idx)
  if (idx == "" && !is.null(default_guess) && default_guess %in% cn) return(default_guess)
  idx_num <- suppressWarnings(as.integer(idx))
  if (is.na(idx_num) || idx_num < 1 || idx_num > length(cn)) stop2("Invalid selection.")
  cn[idx_num]
}

# -----------------------------
# Enforce interactive SOURCE-to-run
# -----------------------------
if (!interactive()) {
  stop2("This script is intended for interactive RStudio SOURCE-to-run usage.\n",
        "Open it in RStudio and press SOURCE.")
}

# -----------------------------
# Output directory policy (MANDATORY)
# -----------------------------
outputs_root <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_root)) dir.create(outputs_root, recursive = TRUE, showWarnings = FALSE)

say("\n==================== Run Setup ====================")
run_name <- readline("Enter a short run name (e.g., Liver5SD_ConstraintModules): ")
run_name <- gsub("[^A-Za-z0-9_\\-]+", "_", trimws(run_name))
if (!nzchar(run_name)) run_name <- "ConstraintModules"

run_id <- paste0(run_name, "_", timestamp_tag())
output_dir <- file.path(outputs_root, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tables_dir <- file.path(output_dir, "tables"); dir.create(tables_dir, showWarnings = FALSE)
plots_dir  <- file.path(output_dir, "plots");  dir.create(plots_dir, showWarnings = FALSE)

say("output_dir: ", normalize_path(output_dir))
say("===================================================\n")

# -----------------------------
# Console identity summary (MANDATORY)
# -----------------------------
say("==================== Script Identity ====================")
say("script_name: ", script_name)
say("script_path: ", script_path)
say("script_full: ", ifelse(is.na(script_full), "NA", script_full))
if (is.na(script_full)) say("NOTE: script path detection failed; fallback filename used: ", known_script_filename)
say("outputs_root: ", normalize_path(outputs_root))
say("run_id:       ", run_id)
say("output_dir:   ", normalize_path(output_dir))
say("=========================================================\n")

# -----------------------------
# Input selection
# -----------------------------
say("==================== Input Selection ====================")
gene_table_path <- pick_file("Select GENE constraint table CSV (must include gene + days_constrained).")
unified_path    <- pick_file("Select UNIFIED_pair_transition_table.csv")
use_map <- pick_yes_no("Do you want to provide a Gene↔Function mapping CSV (for projecting genes into function-modules)?")
map_path <- NA_character_
if (use_map) {
  map_path <- pick_file("Select Gene↔Function mapping CSV (must include gene + function_name).")
}
say("Gene table: ", gene_table_path)
say("Unified:    ", unified_path)
if (use_map) say("Mapping:    ", map_path)
say("=========================================================\n")

# -----------------------------
# Parameters
# -----------------------------
alpha <- 0.01
say("Event threshold alpha = ", alpha)

# -----------------------------
# Load data
# -----------------------------
gene_dt    <- safe_fread(gene_table_path)
unified_dt <- safe_fread(unified_path)
map_dt <- NULL
if (use_map) map_dt <- safe_fread(map_path)

# -----------------------------
# Column selection (biologist-friendly)
# -----------------------------
gene_col <- pick_column(gene_dt, "Select the GENE SYMBOL column in the gene constraint table:", default_guess = "gene")
days_col <- pick_column(gene_dt, "Select the DAYS_CONSTRAINED column (comma-separated list):", default_guess = "days_constrained")

need_unified <- c("transition", "pair")
miss_unified <- setdiff(need_unified, names(unified_dt))
if (length(miss_unified) > 0) stop2("UNIFIED table missing required columns: ", paste(miss_unified, collapse = ", "))

has_is_event <- "is_event" %in% names(unified_dt)
has_p_sumabs <- "p_sumabs" %in% names(unified_dt)
if (!has_is_event && !has_p_sumabs) stop2("UNIFIED table must have either 'is_event' or 'p_sumabs'.")

if (has_is_event) {
  if (is.character(unified_dt$is_event)) {
    unified_dt[, is_event := tolower(trimws(is_event)) %in% c("true", "t", "1", "yes", "y")]
  } else if (is.numeric(unified_dt$is_event)) {
    unified_dt[, is_event := is_event != 0]
  } else {
    unified_dt[, is_event := as.logical(is_event)]
  }
} else {
  unified_dt[, is_event := is.finite(p_sumabs) & (p_sumabs <= alpha)]
}

if (!("obs_sum_abs_delta" %in% names(unified_dt))) unified_dt[, obs_sum_abs_delta := NA_real_]
if (!("direction_bias" %in% names(unified_dt))) unified_dt[, direction_bias := NA_real_]
if (!("focality" %in% names(unified_dt))) {
  if (all(c("obs_n_edges", "obs_sum_abs_delta") %in% names(unified_dt))) {
    unified_dt[, focality := ifelse(obs_n_edges > 0, obs_sum_abs_delta / obs_n_edges, NA_real_)]
  } else {
    unified_dt[, focality := NA_real_]
  }
}

# -----------------------------
# A) Gene constraint-signature modules
# -----------------------------
parse_days <- function(s) {
  if (is.na(s)) return(integer())
  s <- gsub("\\s+", "", as.character(s))
  if (!nzchar(s)) return(integer())
  parts <- unlist(strsplit(s, ",", fixed = TRUE))
  suppressWarnings(as.integer(parts[is.finite(as.integer(parts))]))
}

all_days <- sort(unique(unlist(lapply(gene_dt[[days_col]], parse_days))))
if (length(all_days) == 0) stop2("Could not parse any day values from column: ", days_col)

make_signature <- function(days_vec, day_grid) {
  present <- day_grid %in% days_vec
  paste0(as.integer(present), collapse = "")
}

gene_dt[, gene_symbol := as.character(get(gene_col))]
gene_dt[, days_vec := lapply(get(days_col), parse_days)]
gene_dt[, signature := vapply(days_vec, make_signature, character(1), day_grid = all_days)]
gene_dt[, n_days_sig := nchar(gsub("0", "", signature))]
gene_dt[, module_id_sig := paste0("SIG_", signature)]

gene_sig_membership <- gene_dt[, .(
  gene = gene_symbol,
  n_days_constrained = n_days_sig,
  days_constrained = vapply(days_vec, function(v) paste(v, collapse = ","), character(1)),
  signature = signature,
  module_id = module_id_sig
)]

sig_summary <- gene_sig_membership[, .(
  n_genes = .N,
  n_days = unique(n_days_constrained),
  day_grid = paste(all_days, collapse = ",")
), by = .(module_id, signature)]
setorder(sig_summary, -n_genes, -n_days)

sig_membership_path <- file.path(tables_dir, "GeneModules_ByConstraintSignature_membership.csv")
sig_summary_path    <- file.path(tables_dir, "GeneModules_ByConstraintSignature_summary.csv")
fwrite(gene_sig_membership, sig_membership_path)
fwrite(sig_summary, sig_summary_path)

# -----------------------------
# B) Function/pathway modules per transition from TRUE events
# -----------------------------
split_pair <- function(p) {
  p <- as.character(p)
  parts <- unlist(strsplit(p, "__", fixed = TRUE))
  if (length(parts) != 2) return(c(NA_character_, NA_character_))
  parts
}
pair_parts <- t(vapply(unified_dt$pair, split_pair, character(2)))
unified_dt[, funcA := pair_parts[, 1]]
unified_dt[, funcB := pair_parts[, 2]]

event_dt <- unified_dt[is_event == TRUE & !is.na(funcA) & !is.na(funcB)]
if (nrow(event_dt) == 0) stop2("No TRUE events found in UNIFIED table using alpha = ", alpha)

if (!("weight" %in% names(event_dt))) {
  event_dt[, weight := as.numeric(obs_sum_abs_delta)]
  event_dt[!is.finite(weight), weight := 1]
}

transitions <- sort(unique(event_dt$transition))

membership_list <- list()
module_summary_list <- list()

for (tr in transitions) {
  dt_tr <- event_dt[transition == tr, .(funcA, funcB, weight)]
  if (nrow(dt_tr) == 0) next
  
  g <- graph_from_data_frame(
    data.frame(from = dt_tr$funcA, to = dt_tr$funcB, weight = dt_tr$weight),
    directed = FALSE
  )
  if (ecount(g) == 0 || vcount(g) == 0) next
  
  cl <- tryCatch(cluster_louvain(g, weights = E(g)$weight), error = function(e) NULL)
  if (is.null(cl)) next
  
  mem <- membership(cl)
  mem_dt <- data.table(
    transition = tr,
    function_name = names(mem),
    func_module_id = paste0("TR_", gsub("[^A-Za-z0-9\\-]+", "_", tr), "_M", as.integer(mem))
  )
  membership_list[[tr]] <- mem_dt
  
  mod_sizes <- mem_dt[, .(
    n_functions = .N,
    functions = paste(sort(function_name), collapse = ";")
  ), by = .(transition, func_module_id)]
  module_summary_list[[tr]] <- mod_sizes
}

func_membership_dt <- rbindlist(membership_list, use.names = TRUE, fill = TRUE)
func_modules_dt    <- rbindlist(module_summary_list, use.names = TRUE, fill = TRUE)
if (nrow(func_membership_dt) == 0) stop2("Community detection produced no modules. Check igraph/event content.")

func_membership_path <- file.path(tables_dir, "FunctionModules_ByTransition_membership_long.csv")
func_modules_path    <- file.path(tables_dir, "FunctionModules_ByTransition_summary.csv")
fwrite(func_membership_dt, func_membership_path)
fwrite(func_modules_dt, func_modules_path)

get_set <- function(dt, tr, mid) dt[transition == tr & func_module_id == mid, unique(function_name)]
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

stability_rows <- list()
if (length(transitions) >= 2) {
  for (i in seq_len(length(transitions) - 1)) {
    tr1 <- transitions[i]
    tr2 <- transitions[i + 1]
    mods1 <- unique(func_membership_dt[transition == tr1, func_module_id])
    mods2 <- unique(func_membership_dt[transition == tr2, func_module_id])
    for (m1 in mods1) {
      s1 <- get_set(func_membership_dt, tr1, m1)
      for (m2 in mods2) {
        s2 <- get_set(func_membership_dt, tr2, m2)
        stability_rows[[length(stability_rows) + 1]] <- data.table(
          transition_1 = tr1,
          module_1 = m1,
          transition_2 = tr2,
          module_2 = m2,
          jaccard = jaccard(s1, s2),
          n1 = length(unique(s1)),
          n2 = length(unique(s2))
        )
      }
    }
  }
}
stability_dt <- if (length(stability_rows) > 0) rbindlist(stability_rows) else data.table()
stability_path <- file.path(tables_dir, "FunctionModules_Stability_Jaccard_adjacentTransitions.csv")
fwrite(stability_dt, stability_path)

# -----------------------------
# C) Optional projection: genes → function-modules by transition
# -----------------------------
gene_module_projection_path <- NA_character_
if (use_map) {
  map_gene_col <- pick_column(map_dt, "Select the GENE SYMBOL column in the Gene↔Function mapping table:", default_guess = "gene")
  map_func_col <- pick_column(map_dt, "Select the FUNCTION name column (must match unified function names):", default_guess = "function_name")
  
  map_dt[, gene := as.character(get(map_gene_col))]
  map_dt[, function_name := as.character(get(map_func_col))]
  
  proj_dt <- merge(map_dt[, .(gene, function_name)], func_membership_dt, by = "function_name", allow.cartesian = TRUE)
  
  proj_sum <- proj_dt[, .(
    n_genes = uniqueN(gene),
    genes = paste(sort(unique(gene)), collapse = ";")
  ), by = .(transition, func_module_id)]
  
  gene_module_projection_path <- file.path(tables_dir, "GeneContent_PerFunctionModule_ByTransition.csv")
  fwrite(proj_sum, gene_module_projection_path)
  
  proj_long <- unique(proj_dt[, .(gene, transition, func_module_id)])
  proj_long_path <- file.path(tables_dir, "Gene_to_FunctionModule_membership_long.csv")
  fwrite(proj_long, proj_long_path)
}

# -----------------------------
# Plots (base R, 600 dpi PNG)
# -----------------------------
png600 <- function(path, w = 3600, h = 2400) {
  png(filename = path, width = w, height = h, res = 600)
}

plot1_path <- file.path(plots_dir, "PLOT_events_by_transition.png")
event_counts <- unified_dt[, .(n_events = sum(is_event, na.rm = TRUE)), by = .(transition)]
setorder(event_counts, transition)

png600(plot1_path, w = 4800, h = 2400)
par(mar = c(9, 5, 3, 1))
barplot(event_counts$n_events, names.arg = event_counts$transition, las = 2,
        main = sprintf("Event counts by transition (alpha = %.3f)", alpha),
        ylab = "Number of TRUE events")
dev.off()

plot2_path <- file.path(plots_dir, "PLOT_sign_coherence_by_phase.png")
tmp2 <- unified_dt[is_event == TRUE & is.finite(direction_bias)]
tmp2[, sign_coherence := abs(direction_bias)]
if (!("phase" %in% names(tmp2))) tmp2[, phase := "NA"]

png600(plot2_path, w = 3000, h = 2400)
par(mar = c(7, 5, 3, 1))
boxplot(sign_coherence ~ phase, data = as.data.frame(tmp2),
        main = "Directional coherence by phase (TRUE events)",
        ylab = "|direction_bias|", xlab = "Phase", las = 2, outline = TRUE)
dev.off()

plot3_path <- file.path(plots_dir, "PLOT_focality_by_phase.png")
tmp3 <- unified_dt[is_event == TRUE & is.finite(focality)]
if (!("phase" %in% names(tmp3))) tmp3[, phase := "NA"]

png600(plot3_path, w = 3000, h = 2400)
par(mar = c(7, 5, 3, 1))
boxplot(focality ~ phase, data = as.data.frame(tmp3),
        main = "Focality by phase (TRUE events)",
        ylab = "obs_sum_abs_delta / obs_n_edges", xlab = "Phase", las = 2, outline = TRUE)
dev.off()

# -----------------------------
# Manifest + inventory (MANDATORY)
# -----------------------------
list_deps <- function(pkgs) {
  data.table(
    package = pkgs,
    version = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
  )
}
deps_tbl <- list_deps(c(min_pkgs, analysis_pkgs))
deps_path <- file.path(output_dir, "Dependencies.csv")
fwrite(deps_tbl, deps_path)

inventory_files <- function() {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  data.table(
    file = normalizePath(files, winslash = "/", mustWork = FALSE),
    rel_path = sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/?"), "",
                   normalizePath(files, winslash = "/", mustWork = FALSE)),
    bytes = file.info(files)$size,
    modified = as.character(file.info(files)$mtime)
  )[order(rel_path)]
}

write_manifest <- function(generated_files_dt) {
  manifest <- list(
    run_id = run_id,
    run_timestamp = as.character(Sys.time()),
    script = list(
      name = script_name,
      path = script_path,
      full_path = ifelse(is.na(script_full), NA, script_full)
    ),
    input = list(
      gene_constraint_table = gene_table_path,
      unified_pair_transition_table = unified_path,
      gene_function_mapping = ifelse(use_map, map_path, NA)
    ),
    parameters = list(
      alpha = alpha,
      day_grid = as.list(all_days)
    ),
    dependencies = lapply(seq_len(nrow(deps_tbl)), function(i) as.list(deps_tbl[i])),
    outputs = list(
      outputs_root = normalize_path(outputs_root),
      output_dir = normalize_path(output_dir)
    ),
    generated_files = lapply(seq_len(nrow(generated_files_dt)), function(i) as.list(generated_files_dt[i]))
  )
  json_path <- file.path(output_dir, "Project_Manifest.json")
  jsonlite::write_json(manifest, json_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  json_path
}

inv0 <- inventory_files()
inv_csv_path <- file.path(output_dir, "Project_Manifest_Files.csv")
fwrite(inv0, inv_csv_path)
manifest_json_path <- write_manifest(inv0)

# -----------------------------
# Quarto report (WORKING; NO execute-params; RELATIVE FILE READS)
#   Key idea: render FROM output_dir so relative reads work.
#   This avoids every params-related failure mode you've hit.
# -----------------------------
get_header_block <- function() {
  if (!is.na(script_full) && file.exists(script_full)) {
    lines <- readLines(script_full, warn = FALSE)
    take_n <- min(length(lines), 180)
    head_lines <- lines[seq_len(take_n)]
    keep <- character()
    for (ln in head_lines) {
      if (grepl("^\\s*#|^\\s*$", ln)) keep <- c(keep, ln) else break
    }
    paste(keep, collapse = "\n")
  } else {
    paste0("# Header block unavailable (script_full not detected). Fallback: ", known_script_filename)
  }
}
header_text <- get_header_block()

qmd_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
qmd_basename <- basename(qmd_path)
html_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(qmd_basename), ".html"))

qmd_text <- c(
  "---",
  paste0("title: \"", script_name, " Report\""),
  "format:",
  "  html:",
  "    toc: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "# Summary",
  "This report defines module contents over the developmental time course in two complementary ways:",
  "",
  "- **Gene constraint-signature modules:** genes grouped by which sampled days they are constrained.",
  "- **Function interaction modules per transition:** communities in the TRUE-event function network per transition (Louvain).",
  "",
  "# Script header",
  "```",
  header_text,
  "```",
  "",
  "# Run metadata",
  "",
  paste0("**Run ID:** ", run_id, "  "),
  paste0("**Timestamp:** ", as.character(Sys.time()), "  "),
  paste0("**alpha:** ", alpha, "  "),
  paste0("**Day grid:** ", paste(all_days, collapse = ", "), "  "),
  "",
  "## Inputs",
  paste0("- gene constraint table: ", gene_table_path),
  paste0("- unified table: ", unified_path),
  paste0("- gene↔function mapping (optional): ", ifelse(use_map, map_path, "NA")),
  "",
  "# Dependencies",
  "```{r}",
  "deps <- data.table::fread('Dependencies.csv')",
  "knitr::kable(deps)",
  "```",
  "",
  "# Key outputs",
  "## Gene constraint-signature modules",
  "```{r}",
  "summ <- data.table::fread(file.path('tables', 'GeneModules_ByConstraintSignature_summary.csv'))",
  "knitr::kable(head(summ, 25))",
  "```",
  "",
  "## Function modules by transition",
  "```{r}",
  "fmods <- data.table::fread(file.path('tables', 'FunctionModules_ByTransition_summary.csv'))",
  "knitr::kable(head(fmods, 25))",
  "```",
  "",
  "## Module stability (adjacent transitions)",
  "```{r}",
  "stab_path <- file.path('tables', 'FunctionModules_Stability_Jaccard_adjacentTransitions.csv')",
  "if (file.exists(stab_path)) {",
  "  stab <- data.table::fread(stab_path)",
  "  if ('jaccard' %in% names(stab)) {",
  "    stab <- stab[is.finite(jaccard)]",
  "    stab <- stab[order(-jaccard)]",
  "  }",
  "  knitr::kable(head(stab, 25))",
  "} else {",
  "  cat('No stability table found.')",
  "}",
  "```",
  "",
  "# Figures",
  "```{r}",
  "knitr::include_graphics(file.path('plots', 'PLOT_events_by_transition.png'))",
  "knitr::include_graphics(file.path('plots', 'PLOT_sign_coherence_by_phase.png'))",
  "knitr::include_graphics(file.path('plots', 'PLOT_focality_by_phase.png'))",
  "```",
  "",
  "# Generated files",
  "```{r}",
  "inv <- data.table::fread('Project_Manifest_Files.csv')",
  "knitr::kable(inv)",
  "```",
  "",
  "# Session info",
  "```{r}",
  "sessionInfo()",
  "```"
)
writeLines(qmd_text, qmd_path)

# ---- Render with Quarto CLI (THIS IS THE ONLY RENDER PATH)
render_ok <- TRUE
render_msg <- NULL
out <- character()

quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin)) {
  quarto_bin <- "/Applications/RStudio 2.app/Contents/Resources/app/quarto/bin/quarto"
}

if (!file.exists(quarto_bin)) {
  render_ok <- FALSE
  render_msg <- paste0("Quarto binary not found at: ", quarto_bin)
} else {
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  cmd <- c("render", qmd_basename)
  
  out <- tryCatch(
    system2(quarto_bin, args = cmd, stdout = TRUE, stderr = TRUE),
    error = function(e) { render_ok <<- FALSE; render_msg <<- conditionMessage(e); character() }
  )
  
  if (render_ok && !file.exists(html_path)) {
    render_ok <- FALSE
    render_msg <- paste0(
      "Quarto finished but HTML not found at: ", html_path,
      "\n\nQuarto output:\n", paste(out, collapse = "\n")
    )
  }
  if (!render_ok && is.null(render_msg)) render_msg <- paste(out, collapse = "\n")
}

# Refresh inventory + manifest AFTER rendering (MANDATORY)
inv_final <- inventory_files()
fwrite(inv_final, inv_csv_path)
manifest_json_path <- write_manifest(inv_final)

# -----------------------------
# Return object + "AUTO-RUN WHEN SOURCED" block
# -----------------------------
Constraint_Modules_OverTime <- function() {
  list(
    out_dir = normalize_path(output_dir),
    tables_dir = normalize_path(tables_dir),
    plots_dir = normalize_path(plots_dir),
    gene_signature_membership = normalize_path(sig_membership_path),
    gene_signature_summary = normalize_path(sig_summary_path),
    function_module_membership = normalize_path(func_membership_path),
    function_module_summary = normalize_path(func_modules_path),
    stability_adjacent = normalize_path(stability_path),
    gene_projection = ifelse(is.na(gene_module_projection_path), NULL, normalize_path(gene_module_projection_path)),
    manifest_json = normalize_path(manifest_json_path),
    manifest_files_csv = normalize_path(inv_csv_path),
    report_qmd = normalize_path(qmd_path),
    report_html = if (file.exists(html_path)) normalize_path(html_path) else NULL,
    render_ok = render_ok,
    render_msg = render_msg
  )
}

say("\n==================== DONE ====================")
if (render_ok && file.exists(html_path)) {
  say("HTML report created: ", normalize_path(html_path))
} else {
  say("WARNING: HTML report render failed.")
  if (!is.null(render_msg)) say("Reason: ", render_msg)
  say("QMD written at: ", normalize_path(qmd_path))
}
say("Run folder:       ", normalize_path(output_dir))
say("Manifest JSON:    ", normalize_path(manifest_json_path))
say("File inventory:   ", normalize_path(inv_csv_path))
say("=============================================\n")

# AUTO-RUN WHEN SOURCED (RStudio / interactive use)
if (interactive()) {
  message("Running Constraint_Modules_OverTime() interactively...")
  res <- Constraint_Modules_OverTime()
  message("Done. Results written to:")
  message(res$out_dir)
  if (!is.null(res$gene_signature_summary) && file.exists(res$gene_signature_summary)) {
    dt_view <- tryCatch(data.table::fread(res$gene_signature_summary), error = function(e) NULL)
    if (!is.null(dt_view)) try(View(dt_view), silent = TRUE)
  }
}
