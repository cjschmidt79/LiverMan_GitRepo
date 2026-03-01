# reactome_camera_tidy_module.R
# ======================================================================
# Reactome (msigdbr) + limma::camera enrichment module (TIDY input)
# ======================================================================

# ------------------------- Utilities -------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
sanitize_prefix <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)

ensure_pkgs <- function(pkgs, bioc = character()) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    for (p in bioc) {
      if (!requireNamespace(p, quietly = TRUE)) {
        BiocManager::install(p, ask = FALSE, update = FALSE)
      }
    }
  }
  invisible(TRUE)
}

pick_input_csv <- function() {
  if (!interactive()) stop("Requires interactive() for file.choose().")
  cat("\n📥 Select tidy input CSV (tidy/long counts)\n")
  file.choose()
}

pick_output_dir_with_prefix <- function(prefix_default = "Reactome_Results") {
  if (!interactive()) stop("Requires interactive() for file.choose().")
  cat("\n📂 Select any file inside your desired output folder\n")
  base_dir <- dirname(file.choose())
  pref <- readline(prompt = paste0("Prefix [", prefix_default, "]: "))
  if (pref == "") pref <- prefix_default
  pref <- sanitize_prefix(pref)
  out <- file.path(base_dir, paste0(pref, "_", ts_stamp()))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  out
}

# ------------------------- Column detection + confirm -------------------------

.norm_key <- function(x) {
  x2 <- tolower(as.character(x))
  gsub("[^a-z0-9]", "", x2)
}

.default_synonyms <- list(
  Day = c("day", "time", "timepoint", "tp", "age", "d", "dayposthatch", "posthatchday", "dph"),
  SampleID = c("sampleid", "sample", "bird", "birdid", "animal", "animalid", "individual", "replicate", "rep", "id"),
  Gene = c("gene", "genesymbol", "symbol", "feature", "featurename", "geneid", "ensemblgene", "ensembl", "ensgene", "transcript"),
  Abundance = c("abundance", "count", "counts", "rawcount", "rawcounts", "readcount", "readcounts", "expression", "expr", "value", "n")
)

.detect_column <- function(df_names, default, synonyms) {
  nm_raw <- df_names
  nm_key <- .norm_key(nm_raw)
  default_key <- .norm_key(default)
  syn_keys <- unique(.norm_key(synonyms))
  
  idx <- which(nm_key == default_key)
  if (length(idx) >= 1) return(nm_raw[idx[1]])
  
  idx2 <- which(nm_key %in% syn_keys)
  if (length(idx2) >= 1) return(nm_raw[idx2[1]])
  
  hits <- logical(length(nm_key))
  for (sk in syn_keys) {
    if (nchar(sk) >= 3) hits <- hits | grepl(sk, nm_key, fixed = TRUE)
  }
  idx3 <- which(hits)
  if (length(idx3) == 1) return(nm_raw[idx3])
  
  NA_character_
}

.confirm_or_override_column <- function(df, label, proposed, fallback_default) {
  if (!interactive()) {
    if (!is.na(proposed) && proposed %in% names(df)) return(proposed)
    return(fallback_default)
  }
  
  cat("\n------------------------------\n")
  cat("🔎 Column mapping: ", label, "\n", sep = "")
  cat("   Available columns:\n")
  cat("   - ", paste(names(df), collapse = "\n   - "), "\n", sep = "")
  
  suggested <- proposed
  if (is.na(suggested) || !(suggested %in% names(df))) suggested <- fallback_default
  
  cat("\n✅ Suggested column for ", label, ": '", suggested, "'\n", sep = "")
  cat("   Press Enter to accept, or type the correct column name:\n")
  ans <- readline(prompt = paste0("   ", label, " column [", suggested, "]: "))
  
  chosen <- if (ans == "") suggested else ans
  
  if (!(chosen %in% names(df))) {
    cat("\n❌ Column '", chosen, "' not found.\n", sep = "")
    chosen2 <- readline(prompt = paste0("   ", label, " column (required): "))
    if (!(chosen2 %in% names(df))) stop("Column not found: ", chosen2)
    chosen <- chosen2
  }
  chosen
}

auto_map_columns <- function(df, day_col_default = "day", sample_col_default = "sampleID", 
                             gene_col_default = "gene", abundance_col_default = "Abundance") {
  nms <- names(df)
  syns <- .default_synonyms
  
  prop_day <- .detect_column(nms, day_col_default, syns$Day)
  prop_sample <- .detect_column(nms, sample_col_default, syns$SampleID)
  prop_gene <- .detect_column(nms, gene_col_default, syns$Gene)
  prop_abund <- .detect_column(nms, abundance_col_default, syns$Abundance)
  
  list(
    day_col = .confirm_or_override_column(df, "Day", prop_day, day_col_default),
    sample_col = .confirm_or_override_column(df, "SampleID", prop_sample, sample_col_default),
    gene_col = .confirm_or_override_column(df, "Gene", prop_gene, gene_col_default),
    abundance_col = .confirm_or_override_column(df, "Abundance", prop_abund, abundance_col_default)
  )
}

# ------------------------- Data Processing -------------------------

validate_tidy_input <- function(df, day_col, sample_col, gene_col, abundance_col) {
  req <- c(day_col, sample_col, gene_col, abundance_col)
  missing <- setdiff(req, names(df))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

tidy_to_counts_matrix <- function(df, day_col, sample_col, gene_col, abundance_col, agg_fun = sum) {
  ensure_pkgs(pkgs = "data.table")
  validate_tidy_input(df, day_col, sample_col, gene_col, abundance_col)
  
  dt <- data.table::as.data.table(df)
  data.table::setnames(dt, old = c(day_col, sample_col, gene_col, abundance_col),
                       new = c("Day", "SampleID", "Gene", "Abundance"))
  
  dt[, `:=`(Day = as.character(Day), SampleID = as.character(SampleID), 
            Gene = as.character(Gene), Abundance = as.numeric(Abundance))]
  
  dt_agg <- dt[, .(Abundance = agg_fun(Abundance, na.rm = TRUE)), by = .(Day, SampleID, Gene)]
  
  day_per_sample <- dt_agg[, .(V1 = data.table::uniqueN(Day)), by = SampleID]
  
  if (any(day_per_sample$V1 != 1)) {
    bad <- day_per_sample[V1 != 1, SampleID]
    stop("Each SampleID must map to exactly one Day. Problem SampleIDs: ", paste(head(bad), collapse = ", "))
  }
  
  sample_day <- dt_agg[, .(Day = unique(Day)), by = SampleID]
  sample_day <- sample_day[order(sample_day$Day, sample_day$SampleID)]
  
  wide <- data.table::dcast(dt_agg, Gene ~ SampleID, value.var = "Abundance", fill = 0)
  genes <- wide$Gene
  mat <- as.matrix(wide[, -"Gene"])
  rownames(mat) <- genes
  mat <- mat[, sample_day$SampleID, drop = FALSE]
  
  list(counts = mat, day = factor(sample_day$Day, levels = unique(sample_day$Day)), sample_table = sample_day)
}

# ------------------------- Analysis Engine -------------------------

build_voom_object <- function(counts, day_factor) {
  ensure_pkgs(pkgs = c("limma", "edgeR"), bioc = "edgeR")
  design <- model.matrix(~0 + day_factor)
  colnames(design) <- paste0("Day", levels(day_factor))
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)
  list(v = limma::voom(dge, design), design = design)
}

get_reactome_genesets_all <- function(species = "Homo sapiens", min_genes = 2, universe_genes = NULL) {
  ensure_pkgs(pkgs = "msigdbr")
  msig <- msigdbr::msigdbr(species = species, collection = "C2", subcollection = "REACTOME")
  gs_list <- split(msig$gene_symbol, msig$gs_name)
  if (!is.null(universe_genes)) {
    gs_list <- lapply(gs_list, function(s) intersect(s, universe_genes))
  }
  Filter(function(s) length(s) >= min_genes, gs_list)
}

run_camera_pairwise <- function(v, design, gene_sets, output_dir, file_prefix = "Reactome") {
  day_levels <- sub("^Day", "", colnames(design))
  combos <- combn(day_levels, 2, simplify = FALSE)
  res_list <- list()
  
  for (pair in combos) {
    contrast_label <- paste0("Day", pair[2], " - Day", pair[1])
    message("Running camera for: ", contrast_label)
    cm <- limma::makeContrasts(contrasts = contrast_label, levels = design)
    cres <- limma::camera(v, index = gene_sets, design = design, contrast = cm)
    cres <- cres[order(cres$PValue), ]
    
    write.csv(cres, file.path(output_dir, paste0(file_prefix, "_", pair[2], "_vs_", pair[1], ".csv")))
    res_list[[paste0(pair[2], "_vs_", pair[1])]] <- cres
  }
  invisible(res_list)
}

# ------------------------- Report Function -------------------------

write_quarto_reactome_report <- function(output_dir, manifest = NULL, render = TRUE) {
  qmd_path <- file.path(output_dir, "Reactome_Analysis_Report.qmd")
  log_path <- file.path(output_dir, "execution_log.txt")
  
  # 1. Prepare Input File List text
  input_text <- "No input files recorded."
  if (!is.null(manifest$input_files)) {
    input_names <- sapply(manifest$input_files, function(f) f$name)
    input_text <- paste("- ", input_names, collapse = "\n")
  }
  
  # 2. Read the log if it exists
  log_content <- "No log file found."
  if (file.exists(log_path)) {
    log_content <- readLines(log_path)
  }
  
  report_lines <- c(
    "---",
    "title: \"Reactome Pathway Analysis Report\"",
    paste0("date: \"", Sys.Date(), "\""),
    "format:",
    "  html:",
    "    toc: true",
    "    theme: cosmo",
    "---",
    "",
    "# 1. Provenance (Metadata)",
    "This section tracks the files and scripts used for this specific run.",
    "",
    "### Execution Details",
    paste0("- **Processing Script:** `", manifest$script_name, "`"),
    paste0("- **Script Path:** `", manifest$script_path, "`"),
    paste0("- **Output Directory:** `", output_dir, "`"),
    "",
    "### Input Files Used",
    input_text,
    "",
    "---",
    "",
    "# 2. Analysis Summary",
    paste0("Analysis automatically generated on: ", Sys.time()),
    "",
    "## Methods",
    "1. **Normalization:** Log-CPM transformation via `limma::voom`.",
    "2. **Gene Sets:** All Reactome pathways from `msigdbr` (C2:REACTOME).",
    "3. **Statistics:** Competitive gene set testing via `limma::camera`.",
    "",
    "---",
    "",
    "# 3. Detailed Execution Log",
    "This is the step-by-step history of the analysis run recorded during execution.",
    "",
    "```text",
    log_content,
    "```"
  )
  
  writeLines(report_lines, qmd_path)
  
  if (render && requireNamespace("quarto", quietly = TRUE)) {
    message("🧪 Rendering HTML report...")
    quarto::quarto_render(qmd_path, quiet = TRUE)
  }
}

# ------------------------- Wrapper -------------------------

run_reactome_camera_pipeline_tidy <- function(input_csv = NULL, output_dir = NULL, manifest = NULL, make_report = TRUE) {
  # Fallbacks if called without the launcher
  if (is.null(input_csv)) input_csv <- pick_input_csv()
  if (is.null(output_dir)) output_dir <- pick_output_dir_with_prefix()
  
  df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
  maps <- auto_map_columns(df)
  
  tc <- tidy_to_counts_matrix(df, maps$day_col, maps$sample_col, maps$gene_col, maps$abundance_col)
  vv <- build_voom_object(tc$counts, tc$day)
  gs <- get_reactome_genesets_all(universe_genes = rownames(vv$v))
  
  results <- run_camera_pairwise(vv$v, vv$design, gs, output_dir)
  
  # Use the manifest for the report
  if (make_report) {
    write_quarto_reactome_report(output_dir, manifest = manifest)
  }
  
  message("✅ Analysis complete. Results saved to: ", output_dir)
  return(results)
}