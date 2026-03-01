#!/usr/bin/env Rscript
################################################################################
# Network_Volatility_Engine_no_dplyr.R
#
# NO dplyr / NO tidyr / NO purrr required
# - data.table + base R for all wrangling
# - sequential file read (stable)
# - permutation p-values computed via streaming exceedance counts (memory-safe)
# - optional LIMMA merge (data.table)
# - outputs: full table, FDR<0.05 table, ego-edge CSVs, plots, Quarto report

# Genes were stratified into tiers separating inferentially significant network 
# pivots (FDR < 0.05) from discovery-prioritized structural coordinators defined 
# by temporal persistence and nominal statistical support. A final heuristic tier 
# captured high-persistence network anchors that define stable hepatic identity 
# but do not exhibit punctate rewiring detectable by permutation testing.
################################################################################

# =========================
# 0) Packages (minimal)
# =========================
safe_require <- function(pkgs) {
  miss <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
  invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

safe_require(c("data.table", "ggplot2", "ggrepel", "quarto", "knitr"))

# =========================
# 1) Small helpers
# =========================
is_rstudio <- function() {
  Sys.getenv("RSTUDIO") == "1" &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()
}

pick_dir <- function(caption) {
  if (is_rstudio()) {
    d <- rstudioapi::selectDirectory(caption = caption)
    if (!is.null(d) && d != "") return(d)
  }
  d <- readline(paste0(caption, " (paste full path): "))
  if (d == "") return(NULL)
  d
}

ask_yes_no <- function(prompt) {
  ans <- tolower(trimws(readline(paste0(prompt, " [y/n]: "))))
  ans %in% c("y", "yes")
}

stop_if_cancel <- function(x, msg) {
  if (is.null(x) || is.na(x) || x == "") stop(msg, call. = FALSE)
  x
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# =========================
# 2) Inputs
# =========================
cat("\n--- Input Configuration ---\n")

turnover_dir <- pick_dir("Select directory containing Turnover CSVs")
turnover_dir <- stop_if_cancel(turnover_dir, "User cancelled: Turnover directory is required.")
turnover_dir <- normalizePath(turnover_dir, winslash = "/", mustWork = TRUE)

do_limma <- ask_yes_no("Do you have LIMMA CSVs to integrate?")
limma_dir <- NULL
if (do_limma) {
  limma_dir <- pick_dir("Select LIMMA Results directory")
  limma_dir <- stop_if_cancel(limma_dir, "User cancelled: LIMMA directory required.")
  limma_dir <- normalizePath(limma_dir, winslash = "/", mustWork = TRUE)
}

turnover_files <- list.files(turnover_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(turnover_files)) stop("No CSV files found in turnover directory.", call. = FALSE)

# Sort files: by first number in filename if present, else alphabetical
fname <- basename(turnover_files)
ord_num <- suppressWarnings(as.numeric(sub(".*?(\\d+).*", "\\1", fname)))
ord_num[is.na(ord_num) | ord_num == fname] <- Inf
o <- order(ord_num, fname)
turnover_files <- turnover_files[o]

cat("\nPreview (first file): ", basename(turnover_files[1]), "\n", sep = "")
preview_dt <- data.table::fread(turnover_files[1], nrows = 8)
print(preview_dt)

gene1_col <- readline("Enter column name for Gene1: ")
gene2_col <- readline("Enter column name for Gene2: ")
delta_col <- readline("Enter column name for Delta R (or abs delta): ")
if (gene1_col == "" || gene2_col == "" || delta_col == "") stop("Column names cannot be blank.", call. = FALSE)

run_name <- readline("Enter a run name (e.g., Yolk_to_Feed): ")
if (run_name == "") run_name <- "Analysis"

# Outputs
output_root <- file.path(getwd(), "outputs")
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
output_dir <- file.path(output_root, paste0(run_name, "_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "cytoscape_sif"), recursive = TRUE, showWarnings = FALSE)

# =========================
# 3) Read all edges (stable, sequential)
# =========================
cat("\n--- Reading turnover files (sequential) ---\n")

edge_list <- vector("list", length(turnover_files))
for (i in seq_along(turnover_files)) {
  f <- turnover_files[i]
  dt <- data.table::fread(f, showProgress = FALSE)
  
  needed <- c(gene1_col, gene2_col, delta_col)
  miss <- setdiff(needed, names(dt))
  if (length(miss)) {
    stop(paste0("Missing columns in file ", basename(f), ": ", paste(miss, collapse = ", ")), call. = FALSE)
  }
  
  # Standardize schema (Gene1, Gene2, delta)
  out <- data.table::data.table(
    Gene1 = as.character(dt[[gene1_col]]),
    Gene2 = as.character(dt[[gene2_col]]),
    delta = suppressWarnings(as.numeric(dt[[delta_col]])),
    Source_File = basename(f),
    Source_Path = f
  )
  out <- out[!is.na(Gene1) & !is.na(Gene2) & !is.na(delta)]
  edge_list[[i]] <- out
  
  cat("Read ", basename(f), " (", nrow(out), " edges)\n", sep = "")
}

edge_long <- data.table::rbindlist(edge_list, use.names = TRUE, fill = TRUE)
if (!nrow(edge_long)) stop("No edges after cleaning. Check column mappings.", call. = FALSE)

edge_long[, abs_delta := abs(delta)]

# =========================
# 4) Per-file per-gene metrics using data.table
#    Breadth = unique partners per file
#    Intensity = mean |Δr| per file
# =========================
cat("\n--- Computing Intensity/Breadth/Duration ---\n")

dt1 <- edge_long[, .(Source_File, Gene = Gene1, Partner = Gene2, abs_delta)]
dt2 <- edge_long[, .(Source_File, Gene = Gene2, Partner = Gene1, abs_delta)]
gene_partner <- data.table::rbindlist(list(dt1, dt2))

gene_metrics_by_file <- gene_partner[
  , .(
    Intensity = mean(abs_delta, na.rm = TRUE),
    Breadth = data.table::uniqueN(Partner)
  ),
  by = .(Source_File, Gene)
]

pivot_summary <- gene_metrics_by_file[
  , .(
    Avg_Intensity = mean(Intensity, na.rm = TRUE),
    Max_Intensity = max(Intensity, na.rm = TRUE),
    Total_Breadth = sum(Breadth, na.rm = TRUE),
    Duration = .N
  ),
  by = Gene
]

pivot_summary[, Volatility_Score := Avg_Intensity * Total_Breadth * Duration]
data.table::setorder(pivot_summary, -Volatility_Score)

# =========================
# 5) Permutation p-values (memory-safe streaming)
#    Null: shuffle endpoint labels within each Source_File
# =========================
cat("\n--- Permutation p-values (streaming exceedance counts) ---\n")
n_perm <- suppressWarnings(as.integer(readline("Number of permutations? [Default 1000; use 10000 for final]: ")))
if (is.na(n_perm) || n_perm < 100) n_perm <- 1000

# Split edges by file
edges_split <- split(edge_long, edge_long$Source_File)
files_vec <- names(edges_split)

observed_genes <- pivot_summary$Gene
obs_scores <- pivot_summary$Volatility_Score
names(obs_scores) <- observed_genes

# Exceedance counts: how often null >= observed
exceed <- numeric(length(observed_genes))
names(exceed) <- observed_genes

# Perm loop
for (b in seq_len(n_perm)) {
  # Accumulate per-file per-gene metrics for this permutation
  perm_metrics_all <- NULL
  
  for (fn in files_vec) {
    ed <- edges_split[[fn]]
    nE <- nrow(ed)
    
    endpoints <- c(ed$Gene1, ed$Gene2)
    endpoints_perm <- sample(endpoints, length(endpoints), replace = FALSE)
    
    g1p <- endpoints_perm[1:nE]
    g2p <- endpoints_perm[(nE + 1):(2 * nE)]
    
    # Build permuted partner table (both directions)
    p1 <- data.table::data.table(Source_File = fn, Gene = g1p, Partner = g2p, abs_delta = ed$abs_delta)
    p2 <- data.table::data.table(Source_File = fn, Gene = g2p, Partner = g1p, abs_delta = ed$abs_delta)
    gp <- data.table::rbindlist(list(p1, p2))
    
    pm <- gp[
      , .(
        Intensity = mean(abs_delta, na.rm = TRUE),
        Breadth = data.table::uniqueN(Partner)
      ),
      by = .(Source_File, Gene)
    ]
    
    perm_metrics_all <- if (is.null(perm_metrics_all)) pm else data.table::rbindlist(list(perm_metrics_all, pm))
  }
  
  # Aggregate to perm volatility per gene
  perm_summary <- perm_metrics_all[
    , .(
      Avg_Intensity = mean(Intensity, na.rm = TRUE),
      Total_Breadth = sum(Breadth, na.rm = TRUE),
      Duration = .N
    ),
    by = Gene
  ]
  perm_summary[, Volatility_Score := Avg_Intensity * Total_Breadth * Duration]
  perm_summary <- perm_summary[, .(Gene, Volatility_Score)]
  
  # Compare to observed: only for observed_genes
  # Genes absent in perm_summary imply perm score = 0 (conservative)
  tmp <- data.table::data.table(Gene = observed_genes, Obs = obs_scores)
  tmp <- merge(tmp, perm_summary, by = "Gene", all.x = TRUE)
  tmp[is.na(Volatility_Score), Volatility_Score := 0]
  
  exceed <- exceed + as.numeric(tmp$Volatility_Score >= tmp$Obs)
  
  if (b %% 50 == 0) cat("perm ", b, "/", n_perm, "\n", sep = "")
}

# Empirical p-values with +1 correction
p_value <- (exceed + 1) / (n_perm + 1)
pivot_summary[, p_value := p_value[Gene]]

# BH FDR
pivot_summary[, q_value := p.adjust(p_value, method = "BH")]

# Tiering
pivot_summary[, Status := "Background"]
pivot_summary[q_value < 0.05, Status := "Tier 1: FDR Pivot"]
pivot_summary[q_value >= 0.05 & p_value < 0.05 & Duration > 15, Status := "Tier 2: Structural Coordinator (Discovery)"]
pivot_summary[q_value >= 0.05 & p_value >= 0.05 & p_value < 0.1 & Duration > 10, Status := "Tier 3: Candidate Driver (Discovery)"]
pivot_summary[Duration > 18 & Status == "Background", Status := "Tier 4: High-Persistence Anchor (Heuristic)"]

cat("\nTier counts:\n")
print(table(pivot_summary$Status))

# =========================
# 6) LIMMA integration (optional, data.table)
# =========================
if (do_limma) {
  cat("\n--- Integrating LIMMA ---\n")
  limma_files <- list.files(limma_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(limma_files)) {
    warning("No LIMMA CSVs found; skipping.")
  } else {
    de_all <- rbindlist(lapply(limma_files, function(f) {
      dt <- fread(f, showProgress = FALSE)
      dt[, Source := basename(f)]
      dt
    }), fill = TRUE)
    
    gene_col_guess <- if ("V1" %in% names(de_all)) "V1" else names(de_all)[1]
    setnames(de_all, gene_col_guess, "Limma_Gene")
    
    if (all(c("logFC", "adj.P.Val") %in% names(de_all))) {
      de_sum <- de_all[
        , .(
          Max_abs_logFC = max(abs(as.numeric(logFC)), na.rm = TRUE),
          Min_adjP = min(as.numeric(adj.P.Val), na.rm = TRUE)
        ),
        by = Limma_Gene
      ]
      pivot_summary <- merge(pivot_summary, de_sum, by.x = "Gene", by.y = "Limma_Gene", all.x = TRUE)
    } else {
      warning("LIMMA missing logFC/adj.P.Val; joining only by gene id.")
      de_unique <- unique(de_all[, .(Limma_Gene)])
      pivot_summary <- merge(pivot_summary, de_unique, by.x = "Gene", by.y = "Limma_Gene", all.x = TRUE)
    }
  }
}

# =========================
# 7) Cytoscape ego-edge exports
# =========================
cat("\n--- Exporting ego-edge lists ---\n")

# top pivots: prefer FDR pivots; fallback to top 10 by volatility
top_pivots <- pivot_summary[q_value < 0.05][1:min(10, .N)]
if (nrow(top_pivots) == 0) top_pivots <- pivot_summary[1:min(10, .N)]

for (g in top_pivots$Gene) {
  ego <- edge_long[Gene1 == g | Gene2 == g]
  if (nrow(ego)) {
    ego_out <- data.table::data.table(
      Source_File = ego$Source_File,
      Source = ifelse(ego$Gene1 == g, ego$Gene1, ego$Gene2),
      Target = ifelse(ego$Gene1 == g, ego$Gene2, ego$Gene1),
      Weight = ego$abs_delta
    )
    fwrite(ego_out, file.path(output_dir, "cytoscape_sif", paste0(g, "_ego_edges.csv")))
  }
}

# =========================
# 8) Write tables
# =========================
fwrite(pivot_summary, file.path(output_dir, "Master_Pivot_Analysis_Full.csv"))

filtered <- pivot_summary[q_value < 0.05]
fwrite(filtered, file.path(output_dir, "Master_Pivot_Analysis_Filtered_FDR05.csv"))

# =========================
# 9) Plots (ggplot2 only)
# =========================
cat("\n--- Plotting ---\n")
plot_dt <- as.data.frame(pivot_summary)

# Plot 1: Volatility vs Duration
p1 <- ggplot(plot_dt, aes(x = Duration, y = log10(Volatility_Score + 1), color = Status)) +
  geom_point(alpha = 0.6, size = 3) +
  ggrepel::geom_text_repel(data = head(plot_dt, 10), aes(label = Gene), box.padding = 0.4, max.overlaps = 50) +
  theme_minimal() +
  labs(
    title = "Network Volatility vs Temporal Persistence",
    x = "Duration (number of transition files)",
    y = "log10(Volatility Score + 1)"
  )
ggsave(file.path(output_dir, "plots", "volatility_vs_duration.png"), p1, width = 8, height = 6, dpi = 300)

# Plot 2: Intensity vs Breadth
p2 <- ggplot(plot_dt, aes(x = Total_Breadth, y = Avg_Intensity, color = Status)) +
  geom_point(alpha = 0.6) +
  ggrepel::geom_text_repel(data = head(plot_dt, 8), aes(label = Gene), max.overlaps = 50) +
  theme_minimal() +
  labs(
    title = "Gene Rewiring Role Space",
    x = "Total Breadth (sum unique partners per transition)",
    y = "Average Intensity (mean |Δr|)"
  )
ggsave(file.path(output_dir, "plots", "intensity_vs_breadth.png"), p2, width = 8, height = 6, dpi = 300)

# Plot 3: top gene timeline
top_gene <- plot_dt$Gene[1]
timeline <- gene_metrics_by_file[Gene == top_gene]
timeline_df <- as.data.frame(timeline)
p3 <- ggplot(timeline_df, aes(x = Source_File, y = Intensity, group = 1)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = paste("Activity Profile:", top_gene),
    x = "Transition file",
    y = "Intensity (mean |Δr|)"
  )
ggsave(file.path(output_dir, "plots", "top_pivot_timeline.png"), p3, width = 10, height = 5, dpi = 300)

# --- 10. QUARTO REPORTING (Canonical + Full Provenance + Figures + Tables) ---
cat("\n--- Rendering Quarto Report ---\n")
old_wd <- getwd()
setwd(output_dir)

# Make sure key outputs exist before report tries to read them
full_csv <- "Master_Pivot_Analysis_Full.csv"
fdr_csv  <- "Master_Pivot_Analysis_Filtered_FDR05.csv"

if (!file.exists(full_csv)) {
  warning("Missing Master_Pivot_Analysis_Full.csv; report will be limited.")
}
if (!file.exists(fdr_csv)) {
  # Create empty placeholder if not present (keeps report from failing)
  write.csv(data.frame(), fdr_csv, row.names = FALSE)
}

# Build input file list now (absolute paths)
turnover_files_abs <- list.files(turnover_dir, pattern = "\\.csv$", full.names = TRUE)
turnover_files_abs <- normalizePath(turnover_files_abs, winslash = "/", mustWork = FALSE)

# Output inventory (absolute paths)
out_files_rel <- list.files(output_dir, recursive = TRUE)
out_files_abs <- normalizePath(file.path(output_dir, out_files_rel), winslash = "/", mustWork = FALSE)

# Build QMD content
qmd_filename <- paste0(run_name, "_Report.qmd")

qmd_content <- c(
  "---",
  paste0("title: \"Network Volatility Report: ", run_name, "\""),
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 3",
  "    code-fold: true",
  "    theme: cosmo",
  "    embed-resources: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "# Summary",
  "This report identifies candidate **pivot genes** driving network rewiring across transition windows, using turnover edge lists (Top-50) and optional LIMMA differential expression summaries.",
  "",
  "# Provenance",
  paste0("- **Script path**: `", normalizePath(script_full, winslash = "/", mustWork = FALSE), "`"),
  paste0("- **Turnover input directory**: `", normalizePath(turnover_dir, winslash = "/", mustWork = FALSE), "`"),
  paste0("- **LIMMA integration**: ", ifelse(do_limma, "Enabled", "Disabled")),
  ifelse(do_limma, paste0("- **LIMMA directory**: `", normalizePath(limma_dir, winslash = "/", mustWork = FALSE), "`"), "- **LIMMA directory**: (not used)"),
  paste0("- **Output directory**: `", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "`"),
  paste0("- **Timestamp**: ", timestamp),
  "",
  "## Input files",
  "```{r}",
  "in_files <- readLines('turnover_inputs_manifest.txt')",
  "if(length(in_files)==0) { cat('No input CSVs listed.') } else {",
  "  knitr::kable(data.frame(Input_CSV = in_files))",
  "}",
  "```",
  "",
  "## Output files",
  "```{r}",
  "out_files <- readLines('output_files_manifest.txt')",
  "knitr::kable(data.frame(Output_File = out_files))",
  "```",
  "",
  "# Statistical definitions",
  "Volatility score:",
  "",
  "$$Volatility = \\overline{Intensity} \\times \\sum(Breadth) \\times Duration$$",
  "",
  "Empirical p-value (with +1 correction):",
  "",
  "$$p(\\mathrm{gene}) = \\frac{\\#\\{b : V^{(b)}_{\\mathrm{null}} \\ge V_{\\mathrm{obs}}\\} + 1}{B + 1}$$",
  "",
  "FDR control: Benjamini–Hochberg adjustment of empirical p-values.",
  "",
  "# Key figures",
  "",
  "## Volatility vs duration",
  "```{r}",
  "p <- 'plots/volatility_vs_duration.png'",
  "if (file.exists(p)) knitr::include_graphics(p) else cat('Missing: ', p)",
  "```",
  "",
  "## Intensity vs breadth",
  "```{r}",
  "p <- 'plots/intensity_vs_breadth.png'",
  "if (file.exists(p)) knitr::include_graphics(p) else cat('Missing: ', p)",
  "```",
  "",
  "## Top pivot timeline",
  "```{r}",
  "p <- 'plots/top_pivot_timeline.png'",
  "if (file.exists(p)) knitr::include_graphics(p) else cat('Missing: ', p)",
  "```",
  "",
  "# Top pivot genes",
  "```{r}",
  "full <- if (file.exists('Master_Pivot_Analysis_Full.csv')) read.csv('Master_Pivot_Analysis_Full.csv') else NULL",
  "fdr  <- if (file.exists('Master_Pivot_Analysis_Filtered_FDR05.csv')) read.csv('Master_Pivot_Analysis_Filtered_FDR05.csv') else NULL",
  "if (!is.null(fdr) && nrow(fdr) > 0) {",
  "  cat('Showing FDR < 0.05 (top 15)\\n\\n')",
  "  knitr::kable(head(fdr, 15))",
  "} else if (!is.null(full) && nrow(full) > 0) {",
  "  cat('No FDR < 0.05 hits. Showing top 15 by Volatility_Score.\\n\\n')",
  "  # try to sort if column exists",
  "  if ('Volatility_Score' %in% names(full)) full <- full[order(-full$Volatility_Score), , drop=FALSE]",
  "  knitr::kable(head(full, 15))",
  "} else {",
  "  cat('No results tables found.')",
  "}",
  "```",
  "",
  "# Extended tables",
  "```{r}",
  "csvs <- list.files('.', pattern='\\\\.csv$', full.names=FALSE)",
  "main <- c('Master_Pivot_Analysis_Full.csv','Master_Pivot_Analysis_Filtered_FDR05.csv')",
  "extra <- setdiff(csvs, main)",
  "if (length(extra)==0) {",
  "  cat('No additional CSV tables found in output directory.')",
  "} else {",
  "  for (f in extra) {",
  "    cat('\\n\\n## ', f, '\\n', sep='')",
  "    df <- tryCatch(read.csv(f), error=function(e) NULL)",
  "    if (is.null(df)) { cat('Could not read: ', f) } else { knitr::kable(head(df, 25)) }",
  "  }",
  "}",
  "```",
  "",
  "# Output inventory (relative)",
  "```{r}",
  "knitr::kable(data.frame(File = list.files('.', recursive = TRUE)))",
  "```",
  "",
  "# Session information",
  "```{r}",
  "sessionInfo()",
  "```"
)

# Write manifests used by the QMD (so it can read absolute paths)
writeLines(turnover_files_abs, "turnover_inputs_manifest.txt")
writeLines(out_files_abs, "output_files_manifest.txt")

# Write and render QMD
writeLines(qmd_content, qmd_filename)
quarto::quarto_render(qmd_filename)

# Restore WD
setwd(old_wd)

html_path <- file.path(output_dir, paste0(run_name, "_Report.html"))
cat("\n--- ANALYSIS COMPLETE ---\n")
cat("Absolute Output Path: ", output_dir, "\n")
cat("Report Path: ", html_path, "\n")
if (interactive() && file.exists(html_path)) browseURL(html_path)
