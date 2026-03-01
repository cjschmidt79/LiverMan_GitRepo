# reactome_camera_tidy_module.R
# ======================================================================
# Reactome (msigdbr) + limma::camera enrichment module (TIDY input)
#
# PURPOSE
#   - Read tidy/long RNA-seq counts with metadata
#   - Auto-detect likely column names (case-insensitive; punctuation-insensitive)
#     and confirm interactively
#   - Build voom-normalized model matrix by Day
#   - Run limma::camera competitive gene set testing for ALL Reactome pathways
#     (no "METABOLISM" filtering)
#   - Write one CSV per pairwise Day comparison + optional Quarto report
#
# INPUT (tidy/long) required columns (names can vary; auto-detected):
#   Day, SampleID, Gene, Abundance
#
# DEFAULT expected headers (used if detection fails or non-interactive):
#   day_col       = "day"
#   sample_col    = "sampleID"
#   gene_col      = "gene"
#   abundance_col = "Abundance"
#
# NOTES
#   - Reactome sets retrieved via msigdbr (MSigDB C2:REACTOME) for Homo sapiens.
#   - Your Gene identifiers must match the gene_symbol space used by msigdbr for
#     the selected species (typically human HGNC symbols for Homo sapiens).
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
    cat("   Expected concepts: day, sampleID, gene, abundance (names can vary)\n\n")
    file.choose()
}

pick_output_dir_with_prefix <- function(prefix_default = "Reactome_Results") {
    if (!interactive()) stop("Requires interactive() for file.choose().")
    cat("\n📂 Select any file inside your desired output folder\n")
    base_dir <- dirname(file.choose())
    cat("\n📝 Enter a prefix for the output directory name (default shown):\n")
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

# You can extend these synonyms as needed (they are normalized internally).
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
    
    # 1) exact match to default
    idx <- which(nm_key == default_key)
    if (length(idx) >= 1) return(nm_raw[idx[1]])
    
    # 2) exact match to any synonym
    idx2 <- which(nm_key %in% syn_keys)
    if (length(idx2) >= 1) return(nm_raw[idx2[1]])
    
    # 3) unique contains-match to any synonym key (e.g., "raw_counts" contains "counts")
    hits <- logical(length(nm_key))
    for (sk in syn_keys) {
        if (nchar(sk) >= 3) hits <- hits | grepl(sk, nm_key, fixed = TRUE)
    }
    idx3 <- which(hits)
    if (length(idx3) == 1) return(nm_raw[idx3])
    
    NA_character_
}

.confirm_or_override_column <- function(df, label, proposed, fallback_default) {
    # Non-interactive: use proposed if valid else fallback_default
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
        cat("   Please type one of the available column names exactly.\n")
        chosen2 <- readline(prompt = paste0("   ", label, " column (required): "))
        if (!(chosen2 %in% names(df))) stop("Column not found after retry: ", chosen2)
        chosen <- chosen2
    }
    
    chosen
}

auto_map_columns <- function(df,
                             day_col_default = "day",
                             sample_col_default = "sampleID",
                             gene_col_default = "gene",
                             abundance_col_default = "Abundance",
                             synonyms = .default_synonyms) {
    nms <- names(df)
    
    prop_day <- .detect_column(nms, day_col_default, synonyms$Day)
    prop_sample <- .detect_column(nms, sample_col_default, synonyms$SampleID)
    prop_gene <- .detect_column(nms, gene_col_default, synonyms$Gene)
    prop_abund <- .detect_column(nms, abundance_col_default, synonyms$Abundance)
    
    day_col <- .confirm_or_override_column(df, "Day", prop_day, day_col_default)
    sample_col <- .confirm_or_override_column(df, "SampleID", prop_sample, sample_col_default)
    gene_col <- .confirm_or_override_column(df, "Gene", prop_gene, gene_col_default)
    abundance_col <- .confirm_or_override_column(df, "Abundance", prop_abund, abundance_col_default)
    
    list(
        day_col = day_col,
        sample_col = sample_col,
        gene_col = gene_col,
        abundance_col = abundance_col,
        proposed = list(Day = prop_day, SampleID = prop_sample, Gene = prop_gene, Abundance = prop_abund)
    )
}

# ------------------------- Data validation + tidy → matrix -------------------------

validate_tidy_input <- function(df, day_col, sample_col, gene_col, abundance_col) {
    req <- c(day_col, sample_col, gene_col, abundance_col)
    missing <- setdiff(req, names(df))
    if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
    
    if (anyNA(df[[day_col]])) stop("NA values found in Day column.")
    if (anyNA(df[[sample_col]])) stop("NA values found in SampleID column.")
    if (anyNA(df[[gene_col]])) stop("NA values found in Gene column.")
    if (anyNA(df[[abundance_col]])) stop("NA values found in Abundance column.")
    
    invisible(TRUE)
}

tidy_to_counts_matrix <- function(df,
                                  day_col,
                                  sample_col,
                                  gene_col,
                                  abundance_col,
                                  agg_fun = sum) {
    ensure_pkgs(pkgs = c("data.table"))
    validate_tidy_input(df, day_col, sample_col, gene_col, abundance_col)
    
    dt <- data.table::as.data.table(df)
    data.table::setnames(
        dt,
        old = c(day_col, sample_col, gene_col, abundance_col),
        new = c("Day", "SampleID", "Gene", "Abundance"),
        skip_absent = FALSE
    )
    
    dt[, Day := as.character(Day)]
    dt[, SampleID := as.character(SampleID)]
    dt[, Gene := as.character(Gene)]
    dt[, Abundance := as.numeric(Abundance)]
    
    if (any(!is.finite(dt$Abundance))) stop("Non-finite Abundance values found in abundance column.")
    
    # Aggregate duplicates: multiple rows per (Day, SampleID, Gene)
    dt_agg <- dt[, .(Abundance = agg_fun(Abundance, na.rm = TRUE)), by = .(Day, SampleID, Gene)]
    
    # Enforce one Day per SampleID
    day_per_sample <- dt_agg[, uniqueN(Day), by = SampleID]
    if (any(day_per_sample$V1 != 1)) {
        bad <- day_per_sample[V1 != 1, SampleID]
        stop(
            "Each SampleID must map to exactly one Day. Problem SampleIDs (first 20): ",
            paste(head(bad, 20), collapse = ", "),
            if (length(bad) > 20) " ..." else ""
        )
    }
    
    sample_day <- dt_agg[, .(Day = unique(Day)), by = SampleID]
    sample_day <- sample_day[order(sample_day$Day, sample_day$SampleID)]
    
    # Wide cast: rows = Gene, cols = SampleID; fill missing with 0
    wide <- data.table::dcast(dt_agg, Gene ~ SampleID, value.var = "Abundance", fill = 0)
    
    genes <- wide$Gene
    mat <- as.matrix(wide[, -"Gene"])
    rownames(mat) <- genes
    
    # Align columns to sample_day order
    mat <- mat[, sample_day$SampleID, drop = FALSE]
    
    day_factor <- factor(sample_day$Day, levels = unique(sample_day$Day))
    
    list(counts = mat, day = day_factor, sample_table = sample_day)
}

# ------------------------- Voom + design -------------------------

build_voom_object <- function(counts, day_factor) {
    ensure_pkgs(pkgs = c("limma", "edgeR"), bioc = "edgeR")
    suppressPackageStartupMessages({
        library(limma)
        library(edgeR)
    })
    
    design <- model.matrix(~0 + day_factor)
    colnames(design) <- paste0("Day", levels(day_factor))
    
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge)
    v <- limma::voom(dge, design)
    
    list(v = v, design = design)
}

# ------------------------- Reactome gene sets (ALL Reactome) -------------------------

get_reactome_genesets_all <- function(
        species = "Homo sapiens",
        min_genes = 2,
        universe_genes = NULL
) {
    ensure_pkgs(pkgs = c("msigdbr"))
    suppressPackageStartupMessages(library(msigdbr))
    
    msig <- msigdbr::msigdbr(species = species, collection = "C2", subcollection = "REACTOME")
    if (nrow(msig) == 0) stop("msigdbr returned 0 rows; check species/collection/subcollection.")
    
    geneSetList <- split(msig$gene_symbol, msig$gs_name)
    
    if (!is.null(universe_genes)) {
        universe_genes <- unique(universe_genes)
        geneSetList <- lapply(geneSetList, function(set) intersect(set, universe_genes))
    }
    
    geneSetList <- Filter(function(set) length(set) >= min_genes, geneSetList)
    if (length(geneSetList) == 0) stop("All gene sets filtered out (min_genes/universe too strict or gene IDs mismatch).")
    
    geneSetList
}

# ------------------------- camera pairwise -------------------------

run_camera_pairwise <- function(
        v,
        design,
        gene_sets,
        day_levels = NULL,
        output_dir = NULL,
        file_prefix = "Reactome",
        sort_by = c("PValue", "FDR"),
        verbose = TRUE
) {
    ensure_pkgs(pkgs = c("limma"))
    suppressPackageStartupMessages(library(limma))
    
    sort_by <- match.arg(sort_by)
    
    if (is.null(day_levels)) day_levels <- sub("^Day", "", colnames(design))
    if (length(day_levels) < 2) stop("Need >=2 day levels to do pairwise contrasts.")
    
    combos <- combn(day_levels, 2, simplify = FALSE)
    res_list <- list()
    
    if (!is.null(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (pair in combos) {
        contrast_label <- paste0("Day", pair[2], " - Day", pair[1])
        if (verbose) cat("Running camera for:", contrast_label, "\n")
        
        cm <- makeContrasts(contrasts = contrast_label, levels = design)
        cres <- camera(v, index = gene_sets, design = design, contrast = cm)
        
        if (sort_by %in% names(cres)) cres <- cres[order(cres[[sort_by]]), , drop = FALSE]
        
        key <- paste0(pair[2], "_vs_", pair[1])
        res_list[[key]] <- cres
        
        if (!is.null(output_dir)) {
            outname <- paste0(file_prefix, "_", key, ".csv")
            write.csv(cres, file = file.path(output_dir, outname))
        }
    }
    
    invisible(res_list)
}

# ------------------------- Optional Quarto report -------------------------

write_quarto_reactome_report <- function(
        output_dir,
        input_path = NULL,
        script_name = "Reactome_Tidy_Module_Run",
        script_path = getwd(),
        title_prefix = "Reactome Pathway Analysis Report",
        render = TRUE
) {
    qmd_path <- file.path(output_dir, paste0(script_name, "_Reactome_Report.qmd"))
    html_path <- sub("\\.qmd$", ".html", qmd_path)
    output_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
    
    voom_formula <- "$$\\text{logCPM}_{ij} = \\log_2\\left(\\frac{Y_{ij}}{L_j} \\times 10^6 + 0.5\\right)$$"
    camera_description <- "Competitive gene set testing using linear models with inter-gene correlation adjustment (`limma::camera`)."
    
    report_lines <- c(
        "---",
        paste0("title: \"", title_prefix, " — ", script_name, "\""),
        "format: html",
        "editor: visual",
        "---",
        "",
        "# Reactome Enrichment Report (ALL Reactome Pathways)",
        "",
        "## Input",
        if (!is.null(input_path)) paste0("- Tidy CSV: `", basename(input_path), "`") else "- Tidy CSV: (not recorded)",
        if (!is.null(input_path)) paste0("- Full path: `", input_path, "`") else NULL,
        "",
        "## Methods",
        "### Voom transformation",
        voom_formula,
        "",
        "### Camera",
        camera_description,
        "",
        "## Outputs",
        paste0("- Output directory: `", output_dir, "`"),
        "- CSV result files:",
        if (length(output_files)) paste0("  - ", paste(basename(output_files), collapse = "\n  - ")) else "  - (none found)",
        ""
    )
    
    writeLines(report_lines[!vapply(report_lines, is.null, logical(1))], qmd_path)
    cat("\n📝 Quarto .qmd written to:\n", qmd_path, "\n")
    
    if (render) {
        if (!requireNamespace("quarto", quietly = TRUE)) {
            warning("quarto R package not available; wrote .qmd but did not render HTML.")
            return(invisible(list(qmd = qmd_path, html = html_path)))
        }
        cat("\n🧪 Rendering HTML report...\n")
        quarto::quarto_render(qmd_path)
        cat("✅ HTML report saved to:\n", html_path, "\n")
        if (interactive()) browseURL(html_path)
    }
    
    invisible(list(qmd = qmd_path, html = html_path))
}

# ------------------------- High-level wrapper (tidy pipeline; ALL Reactome) -------------------------

run_reactome_camera_pipeline_tidy <- function(
        input_csv = NULL,
        output_dir = NULL,
        
        # Defaults you requested
        day_col = "day",
        sample_col = "sampleID",
        gene_col = "gene",
        abundance_col = "Abundance",
        
        # Reactome settings (ALL pathways; no regex filter)
        species = "Homo sapiens",
        min_genes = 2,
        
        file_prefix = "Reactome",
        make_report = TRUE,
        verbose = TRUE
) {
    ensure_pkgs(pkgs = c("data.table", "limma", "msigdbr", "edgeR"), bioc = "edgeR")
    
    if (is.null(input_csv)) input_csv <- pick_input_csv()
    if (is.null(output_dir)) output_dir <- pick_output_dir_with_prefix("Reactome_Results")
    
    df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
    
    # Auto-detect + confirm (interactive) without changing the interface
    mapped <- auto_map_columns(
        df,
        day_col_default = day_col,
        sample_col_default = sample_col,
        gene_col_default = gene_col,
        abundance_col_default = abundance_col,
        synonyms = .default_synonyms
    )
    
    day_col <- mapped$day_col
    sample_col <- mapped$sample_col
    gene_col <- mapped$gene_col
    abundance_col <- mapped$abundance_col
    
    if (verbose) {
        cat("\n✅ Final column mapping:\n")
        cat("   Day       :", day_col, "\n")
        cat("   SampleID  :", sample_col, "\n")
        cat("   Gene      :", gene_col, "\n")
        cat("   Abundance :", abundance_col, "\n")
    }
    
    # Convert tidy -> counts + day
    tc <- tidy_to_counts_matrix(
        df = df,
        day_col = day_col,
        sample_col = sample_col,
        gene_col = gene_col,
        abundance_col = abundance_col,
        agg_fun = sum
    )
    
    # Voom
    vv <- build_voom_object(tc$counts, tc$day)
    
    # Reactome sets restricted to genes in voom object (ALL Reactome)
    gene_sets <- get_reactome_genesets_all(
        species = species,
        min_genes = min_genes,
        universe_genes = rownames(vv$v)
    )
    
    # camera
    res <- run_camera_pairwise(
        v = vv$v,
        design = vv$design,
        gene_sets = gene_sets,
        output_dir = output_dir,
        file_prefix = file_prefix,
        verbose = verbose
    )
    
    # Report
    if (make_report) {
        script_info <- tryCatch({
            if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
                this_path <- rstudioapi::getSourceEditorContext()$path
            } else {
                this_path <- sys.frame(1)$ofile
            }
            if (is.null(this_path) || this_path == "") {
                list(name = "Reactome_Tidy_Module_Run", path = getwd())
            } else {
                list(
                    name = tools::file_path_sans_ext(basename(this_path)),
                    path = normalizePath(dirname(this_path))
                )
            }
        }, error = function(e) list(name = "Reactome_Tidy_Module_Run", path = getwd()))
        
        write_quarto_reactome_report(
            output_dir = output_dir,
            input_path = input_csv,
            script_name = script_info$name,
            script_path = script_info$path,
            render = TRUE
        )
    }
    
    cat("\n✅ Reactome camera() complete (ALL Reactome).\n")
    cat("   Output dir: ", output_dir, "\n", sep = "")
    cat("   Days: ", paste(levels(tc$day), collapse = ", "), "\n", sep = "")
    cat("   Gene sets used: ", length(gene_sets), "\n", sep = "")
    
    invisible(list(
        input_csv = input_csv,
        output_dir = output_dir,
        sample_table = tc$sample_table,
        design = vv$design,
        n_gene_sets = length(gene_sets),
        results = res,
        column_mapping = list(
            day_col = day_col,
            sample_col = sample_col,
            gene_col = gene_col,
            abundance_col = abundance_col,
            proposed = mapped$proposed
        )
    ))
}

# ======================================================================
# QUICK START (example run script)
#
# source("reactome_camera_tidy_module.R")
# out <- run_reactome_camera_pipeline_tidy(
#   input_csv = NULL,   # interactive file picker
#   output_dir = NULL,  # interactive output dir picker
#   make_report = TRUE
# )
# ======================================================================
