#!/usr/bin/env Rscript
# ============================================================
# PC1 Driver Genes: Mean-PCA vs CV-PCA (Mac-native dialogs; NO Tcl/Tk)
#
# What it does:
#   - Select Mean and CV PCA gene-score CSVs (Finder / RStudio dialogs)
#   - Rank PC1 driver genes by |PC1|
#   - Quantify overlap (Jaccard, overlap counts)
#   - Compute global PC1 similarity across genes (cosine + Pearson)
#   - Write clean output tables + plots
#   - Biological annotation via:
#       A) GO enrichment (clusterProfiler + org.Gg.eg.db) if available
#       B) GMT enrichment (fgsea) if you provide gmt=/path/to/genesets.gmt
#       C) Optional merge of user-provided annotation CSV
#
# Inputs (required):
#   - PCA_Mean_Liver_Scores.csv (rows=genes; columns include Gene-ish ID + PC1..)
#   - PCA_CoV_Liver_Scores.csv  (rows=genes; columns include Gene-ish ID + PC1..)
#
# Usage (interactive, Mac-native dialogs):
#   Rscript PC1_Drivers_Overlap_Enrichment_MacNative.R
#
# Usage (non-interactive, pass paths explicitly):
#   Rscript PC1_Drivers_Overlap_Enrichment_MacNative.R \
#     mean_scores="/path/to/PCA_Mean_Liver_Scores.csv" \
#     cov_scores="/path/to/PCA_CoV_Liver_Scores.csv" \
#     outdir="/path/to/output" \
#     topN=50 enrich=AUTO ont=BP
#
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos="https://cloud.r-project.org")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="https://cloud.r-project.org")
  library(data.table)
  library(ggplot2)
})

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)

arg_val <- function(flag, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

is_nonempty <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)

# ---- Script path capture (best-effort) ----
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  p <- sub("^--file=", "", hit[1])
  tryCatch(normalizePath(p, winslash="/", mustWork=FALSE), error=function(e) p)
}

# ============================================================
# Mac-native dialogs (NO Tcl/Tk)
#   Priority:
#     1) RStudio dialogs (selectFile/selectDirectory) if available
#     2) macOS Finder dialogs via AppleScript (osascript)
#   Non-interactive: must pass paths via args
# ============================================================

is_rstudio <- function() {
  Sys.getenv("RSTUDIO") == "1" || Sys.getenv("TERM_PROGRAM") == "RStudio"
}

choose_file_mac <- function(prompt = "Choose a file") {
  if (!interactive()) stop2("Non-interactive run: provide file paths via mean_scores= and cov_scores= arguments.")
  
  # RStudio-native
  if (is_rstudio() && requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- rstudioapi::selectFile(caption = prompt)
    if (is.character(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  # Finder via AppleScript
  if (Sys.info()[["sysname"]] == "Darwin") {
    script <- paste0(
      'set theFile to choose file with prompt "', gsub('"', '\\"', prompt), '"\n',
      'POSIX path of theFile'
    )
    p <- tryCatch(system2("osascript", c("-e", script), stdout = TRUE, stderr = TRUE),
                  error = function(e) character())
    if (length(p) > 0) {
      p <- trimws(p[1])
      if (nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  
  stop2("Could not open a Mac-native file chooser. If you are in RStudio, install rstudioapi; otherwise run with explicit arguments.")
}

choose_dir_mac <- function(prompt = "Choose a folder") {
  if (!interactive()) stop2("Non-interactive run: provide outdir= argument.")
  
  # RStudio-native
  if (is_rstudio() && requireNamespace("rstudioapi", quietly = TRUE)) {
    d <- rstudioapi::selectDirectory(caption = prompt)
    if (is.character(d) && nzchar(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  
  # Finder via AppleScript
  if (Sys.info()[["sysname"]] == "Darwin") {
    script <- paste0(
      'set theFolder to choose folder with prompt "', gsub('"', '\\"', prompt), '"\n',
      'POSIX path of theFolder'
    )
    d <- tryCatch(system2("osascript", c("-e", script), stdout = TRUE, stderr = TRUE),
                  error = function(e) character())
    if (length(d) > 0) {
      d <- trimws(d[1])
      if (nzchar(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  
  stop2("Could not open a Mac-native folder chooser. If you are in RStudio, install rstudioapi; otherwise run with outdir= argument.")
}

# ---- Math helpers ----
cos_sim <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(x * y) / den
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

# ---- Column inference ----
infer_gene_col <- function(DT) {
  candidates <- c("Gene","gene","SYMBOL","Symbol","symbol","ID","id","Feature","feature","Ensembl","ensembl")
  hit <- intersect(candidates, names(DT))
  if (length(hit) > 0) return(hit[1])
  
  # otherwise choose first non-PC column
  non_pc <- names(DT)[!grepl("^PC[0-9]+$", names(DT))]
  if (length(non_pc) == 0) stop2("Could not infer gene ID column: only PC columns found.")
  non_pc[1]
}

require_pc1 <- function(DT, label="table") {
  if (!("PC1" %in% names(DT))) stop2("PC1 column not found in ", label, ". Columns: ", paste(names(DT), collapse=", "))
}

# ---- Optional annotation merge (user provided) ----
read_annotation_table <- function(path) {
  if (!is_nonempty(path) || !file.exists(path)) return(NULL)
  A <- fread(path)
  gc <- infer_gene_col(A)
  setnames(A, gc, "Gene")
  A
}

# ---- Enrichment options ----
have_clusterprof <- requireNamespace("clusterProfiler", quietly = TRUE)
have_orggg <- requireNamespace("org.Gg.eg.db", quietly = TRUE) # Gallus gallus
have_fgsea <- requireNamespace("fgsea", quietly = TRUE)

maybe_install <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")
}

run_go_enrichment_gallus <- function(gene_symbols, universe_symbols, out_csv, ont="BP") {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Gg.eg.db)
  })
  
  map <- suppressWarnings(clusterProfiler::bitr(
    unique(c(gene_symbols, universe_symbols)),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Gg.eg.db::org.Gg.eg.db
  ))
  
  set_e <- unique(map[SYMBOL %in% gene_symbols, ENTREZID])
  uni_e <- unique(map[SYMBOL %in% universe_symbols, ENTREZID])
  
  if (length(set_e) < 10) {
    fwrite(data.table(note="Too few mapped genes for GO enrichment", n_set=length(set_e)), out_csv)
    return(invisible(NULL))
  }
  
  ego <- clusterProfiler::enrichGO(
    gene          = set_e,
    universe      = uni_e,
    OrgDb         = org.Gg.eg.db::org.Gg.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
  
  res <- as.data.table(ego@result)
  if (nrow(res) == 0) {
    fwrite(data.table(note="No significant GO terms", ont=ont), out_csv)
    return(invisible(NULL))
  }
  fwrite(res, out_csv)
  invisible(res)
}

read_gmt <- function(gmt_path) {
  if (!is_nonempty(gmt_path) || !file.exists(gmt_path)) return(NULL)
  lines <- readLines(gmt_path, warn = FALSE)
  sets <- lapply(lines, function(s) {
    parts <- strsplit(s, "\t")[[1]]
    if (length(parts) < 3) return(NULL)
    list(name = parts[1], genes = unique(parts[3:length(parts)]))
  })
  sets <- sets[!vapply(sets, is.null, logical(1))]
  setNames(lapply(sets, `[[`, "genes"), vapply(sets, `[[`, character(1), "name"))
}

run_gmt_enrichment_fgsea <- function(ranks_named, gmt_sets, out_csv, minSize=10, maxSize=500) {
  maybe_install(c("fgsea"))
  suppressPackageStartupMessages(library(fgsea))
  
  if (is.null(gmt_sets) || length(gmt_sets) == 0) {
    fwrite(data.table(note="No GMT sets loaded"), out_csv)
    return(invisible(NULL))
  }
  
  ranks_named <- ranks_named[is.finite(ranks_named)]
  ranks_named <- sort(ranks_named, decreasing = TRUE)
  if (length(ranks_named) < 50) {
    fwrite(data.table(note="Too few ranked genes for fgsea", n=length(ranks_named)), out_csv)
    return(invisible(NULL))
  }
  
  fg <- fgsea::fgsea(
    pathways = gmt_sets,
    stats    = ranks_named,
    minSize  = minSize,
    maxSize  = maxSize
  )
  res <- as.data.table(fg)[order(padj, -abs(NES))]
  fwrite(res, out_csv)
  invisible(res)
}

# ============================================================
# MAIN
# ============================================================

cat("\n=== PC1 Driver Comparison: Mean-PCA vs CV-PCA (Mac-native; NO Tcl/Tk) ===\n")

mean_path <- arg_val("mean_scores", NA_character_)
cov_path  <- arg_val("cov_scores",  NA_character_)
outdir    <- arg_val("outdir",      NA_character_)

topN <- suppressWarnings(as.integer(arg_val("topN", "50")))
if (!is.finite(topN) || topN < 10) topN <- 50

enrich_mode <- toupper(arg_val("enrich", "AUTO"))  # AUTO | GO | GMT | NONE
ont <- toupper(arg_val("ont", "BP"))               # BP | MF | CC
if (!ont %in% c("BP","MF","CC")) ont <- "BP"

annotation_csv <- arg_val("annotation_csv", NA_character_)
gmt_file <- arg_val("gmt", NA_character_)

# Mac-native selection if not provided
if (!is_nonempty(mean_path)) mean_path <- choose_file_mac("Select PCA_Mean_Liver_Scores.csv")
if (!is_nonempty(cov_path))  cov_path  <- choose_file_mac("Select PCA_CoV_Liver_Scores.csv")
if (!is_nonempty(outdir))    outdir    <- choose_dir_mac("Select output directory")

if (!file.exists(mean_path)) stop2("Mean scores file not found: ", mean_path)
if (!file.exists(cov_path))  stop2("CV scores file not found: ", cov_path)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Plots"), recursive = TRUE, showWarnings = FALSE)

cat("\nInputs:\n")
cat("  Mean scores: ", mean_path, "\n", sep="")
cat("  CV scores  : ", cov_path, "\n", sep="")
cat("Output:\n")
cat("  outdir     : ", outdir, "\n", sep="")
cat("Settings:\n")
cat("  topN       : ", topN, "\n", sep="")
cat("  enrich     : ", enrich_mode, " (ont=", ont, ")\n", sep="")
if (is_nonempty(annotation_csv)) cat("  annotation_csv: ", annotation_csv, "\n", sep="")
if (is_nonempty(gmt_file)) cat("  gmt file      : ", gmt_file, "\n", sep="")

# Read
Mean <- fread(mean_path)
CoV  <- fread(cov_path)

require_pc1(Mean, "Mean scores")
require_pc1(CoV,  "CV scores")

gcol_mean <- infer_gene_col(Mean)
gcol_cov  <- infer_gene_col(CoV)

setnames(Mean, gcol_mean, "Gene")
setnames(CoV,  gcol_cov,  "Gene")

Mean <- Mean[!is.na(Gene) & nzchar(as.character(Gene))]
CoV  <- CoV[!is.na(Gene) & nzchar(as.character(Gene))]

Mean[, PC1 := as.numeric(PC1)]
CoV[,  PC1 := as.numeric(PC1)]

# Align by gene
M <- merge(
  Mean[, .(Gene, PC1_mean = PC1)],
  CoV[,  .(Gene, PC1_cv   = PC1)],
  by = "Gene", all = FALSE
)

if (nrow(M) < 50) warning("Only ", nrow(M), " matched genes. Similarity/overlap may be unstable.")

# Similarity across all genes
sim <- data.table(
  matched_genes = nrow(M),
  cosine_PC1_scores = cos_sim(M$PC1_mean, M$PC1_cv),
  pearson_PC1_scores = suppressWarnings(cor(M$PC1_mean, M$PC1_cv, use="pairwise.complete.obs"))
)
fwrite(sim, file.path(outdir, "PC1_Score_Similarity.csv"))

# Rank drivers
M[, abs_PC1_mean := abs(PC1_mean)]
M[, abs_PC1_cv   := abs(PC1_cv)]

top_mean <- M[order(-abs_PC1_mean)][1:min(topN, .N)]
top_cv   <- M[order(-abs_PC1_cv)][1:min(topN, .N)]

overlap_genes <- intersect(top_mean$Gene, top_cv$Gene)
overlap_tbl <- M[Gene %in% overlap_genes][
  order(-pmax(abs_PC1_mean, abs_PC1_cv)),
  .(Gene, PC1_mean, PC1_cv, abs_PC1_mean, abs_PC1_cv)
]

fwrite(top_mean[, .(Gene, PC1_mean, abs_PC1_mean)], file.path(outdir, "TopGenes_Mean_PC1.csv"))
fwrite(top_cv[,   .(Gene, PC1_cv,   abs_PC1_cv)],   file.path(outdir, "TopGenes_CV_PC1.csv"))
fwrite(overlap_tbl, file.path(outdir, "TopGenes_Overlap_PC1.csv"))

jac <- jaccard(top_mean$Gene, top_cv$Gene)

summary_lines <- c(
  "PC1 Driver Overlap Summary",
  paste0("Timestamp: ", ts_stamp()),
  paste0("Script: ", get_script_path()),
  paste0("Mean scores: ", mean_path),
  paste0("CV scores: ", cov_path),
  paste0("Output dir: ", outdir),
  paste0("Matched genes: ", nrow(M)),
  paste0("TopN: ", nrow(top_mean)),
  "",
  "PC1 score-vector similarity across all matched genes:",
  paste0("  cosine  = ", sprintf("%.6f", sim$cosine_PC1_scores)),
  paste0("  pearson = ", sprintf("%.6f", sim$pearson_PC1_scores)),
  "",
  "TopN overlap:",
  paste0("  overlap_n = ", length(overlap_genes)),
  paste0("  Jaccard   = ", sprintf("%.6f", jac)),
  "",
  "Files written:",
  "  - PC1_Score_Similarity.csv",
  "  - TopGenes_Mean_PC1.csv",
  "  - TopGenes_CV_PC1.csv",
  "  - TopGenes_Overlap_PC1.csv"
)
writeLines(summary_lines, con = file.path(outdir, "Overlap_Summary_PC1.txt"))

# Optional annotation merge
Anno <- read_annotation_table(annotation_csv)
if (!is.null(Anno)) {
  top_mean_anno <- merge(top_mean, Anno, by="Gene", all.x=TRUE)
  top_cv_anno   <- merge(top_cv,   Anno, by="Gene", all.x=TRUE)
  fwrite(top_mean_anno, file.path(outdir, "TopGenes_Mean_PC1_Annotated.csv"))
  fwrite(top_cv_anno,   file.path(outdir, "TopGenes_CV_PC1_Annotated.csv"))
}

# Plots
p_bar <- rbind(
  top_mean[, .(Gene, absPC1 = abs_PC1_mean, Set = "Mean_PC1")],
  top_cv[,   .(Gene, absPC1 = abs_PC1_cv,   Set = "CV_PC1")]
)

p1 <- ggplot(p_bar, aes(x = reorder(Gene, absPC1), y = absPC1)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Set, scales="free_y") +
  labs(title = paste0("Top ", topN, " genes by |PC1|"),
       x = NULL, y = "|PC1 score|") +
  theme_classic(base_size = 12)

ggsave(file.path(outdir, "Plots", "PC1_TopGeneScores.png"), p1, width=9, height=8, dpi=300)

p2 <- ggplot(M, aes(x = PC1_mean, y = PC1_cv)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth=0.3) +
  geom_vline(xintercept = 0, linewidth=0.3) +
  labs(title = "PC1 gene scores: Mean-PCA vs CV-PCA",
       subtitle = paste0("cosine=", sprintf("%.3f", sim$cosine_PC1_scores),
                         " | pearson=", sprintf("%.3f", sim$pearson_PC1_scores),
                         " | Jaccard(topN)=", sprintf("%.3f", jac)),
       x = "PC1 (Mean-PCA scores)", y = "PC1 (CV-PCA scores)") +
  theme_classic(base_size = 12)

ggsave(file.path(outdir, "Plots", "PC1_ScoreScatter.png"), p2, width=7, height=6, dpi=300)

# Enrichment / biological annotation
do_go  <- (enrich_mode %in% c("AUTO","GO"))  && have_clusterprof && have_orggg
do_gmt <- (enrich_mode %in% c("AUTO","GMT")) && is_nonempty(gmt_file) && file.exists(gmt_file)

universe <- M$Gene

if (enrich_mode == "NONE") {
  message("\nEnrichment disabled (enrich=NONE).")
} else if (do_go) {
  message("\nRunning GO enrichment (Gallus; clusterProfiler; ont=", ont, ") ...")
  run_go_enrichment_gallus(
    gene_symbols = top_mean$Gene,
    universe_symbols = universe,
    out_csv = file.path(outdir, paste0("Enrichment_Mean_PC1_GO_", ont, ".csv")),
    ont = ont
  )
  run_go_enrichment_gallus(
    gene_symbols = top_cv$Gene,
    universe_symbols = universe,
    out_csv = file.path(outdir, paste0("Enrichment_CV_PC1_GO_", ont, ".csv")),
    ont = ont
  )
} else if (do_gmt) {
  message("\nRunning GMT enrichment (fgsea) ...")
  gsets <- read_gmt(gmt_file)
  ranks_mean <- setNames(M$PC1_mean, M$Gene)
  ranks_cv   <- setNames(M$PC1_cv,   M$Gene)
  
  run_gmt_enrichment_fgsea(
    ranks_named = ranks_mean,
    gmt_sets = gsets,
    out_csv = file.path(outdir, "Enrichment_Mean_PC1_GMT_fgsea.csv")
  )
  run_gmt_enrichment_fgsea(
    ranks_named = ranks_cv,
    gmt_sets = gsets,
    out_csv = file.path(outdir, "Enrichment_CV_PC1_GMT_fgsea.csv")
  )
} else {
  message("\nEnrichment not run.")
  message("To enable GO enrichment for Gallus gallus, install (Bioconductor):")
  message("  install.packages('BiocManager'); BiocManager::install(c('clusterProfiler','org.Gg.eg.db'))")
  message("Or provide a GMT file for fgsea: gmt=/path/to/genesets.gmt and set enrich=GMT or AUTO.")
}

cat("\nDONE\n")
cat("Outputs in: ", outdir, "\n", sep="")
