#!/usr/bin/env Rscript
# ============================================================
# ABC_Gene_Column_Enrichment_Analysis.R  (UPDATED for fixed-A design)
#
# Key update: When A is balanced by construction, "A enrichment" is not meaningful.
# We instead interrogate the triplet universe via:
#  - Predictor polarity (B vs C) per gene (primary result)
#  - Conditional predictor selection given A: P(B|A), P(C|A)
#  - Pair structure: P(B,C) (and optionally per A)
#
# Biologist-first SOURCE-to-run (base R; optional rstudioapi/jsonlite)
# - Interactive file pickers
# - Standardized outputs/<run_name>_<timestamp>/
# - Captures script identity + git metadata (if available)
# - Writes provenance QMD/HTML + manifests + output inventory
# ============================================================

# -----------------------------
# Small utilities (minimal deps)
# -----------------------------
safe_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

quiet_require <- function(pkg) {
  suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))
}

ensure_pkg <- function(pkg) {
  if (!quiet_require(pkg)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    if (!quiet_require(pkg)) stop("Package install/load failed: ", pkg)
  }
  TRUE
}

pick_file <- function(prompt = "Select a file") {
  cat("\n", prompt, "\n", sep = "")
  if (quiet_require("rstudioapi") && rstudioapi::isAvailable()) {
    return(rstudioapi::selectFile(caption = prompt))
  }
  return(file.choose())
}

capture_script_identity <- function(verbose = TRUE) {
  script_path <- NA_character_
  
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) > 0) script_path <- sub(file_arg, "", hit[1])
  
  if (is.na(script_path) || !nzchar(script_path)) {
    if (!is.null(sys.frames()[[1]]$ofile)) script_path <- sys.frames()[[1]]$ofile
  }
  
  if ((is.na(script_path) || !nzchar(script_path)) && quiet_require("rstudioapi") && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) script_path <- ctx$path
  }
  
  script_path <- tryCatch(normalizePath(script_path, winslash = "/", mustWork = FALSE), error = function(e) script_path)
  script_name <- if (!is.na(script_path) && nzchar(script_path)) basename(script_path) else "UNKNOWN_SCRIPT.R"
  
  out <- list(
    script_name = script_name,
    script_path = script_path,
    timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    wd          = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  )
  
  if (isTRUE(verbose)) {
    cat("\n--- Script identity ---\n")
    cat("script_name: ", out$script_name, "\n", sep = "")
    cat("script_path: ", out$script_path, "\n", sep = "")
    cat("timestamp  : ", out$timestamp, "\n", sep = "")
    cat("wd         : ", out$wd, "\n", sep = "")
  }
  out
}

get_git_metadata <- function() {
  run <- function(cmd) tryCatch(system(cmd, intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
  list(
    git_branch = paste(run("git rev-parse --abbrev-ref HEAD"), collapse = " "),
    git_commit = paste(run("git rev-parse HEAD"), collapse = " "),
    git_status = paste(run("git status --porcelain"), collapse = "\n")
  )
}

write_text <- function(path, lines) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
}

write_csv <- function(df, path) utils::write.csv(df, path, row.names = FALSE, quote = TRUE)

inventory_outputs <- function(out_dir) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir == FALSE]
  info <- file.info(files)
  data.frame(
    file = normalizePath(files, winslash = "/", mustWork = FALSE),
    bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

write_manifest <- function(out_dir, run_name, script_id, git_meta, inputs, params) {
  manifest <- list(
    run_name = run_name,
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    output_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    script = script_id,
    git = git_meta,
    inputs = inputs,
    parameters = params
  )
  
  manifest_json_path <- file.path(out_dir, "Project_Manifest.json")
  
  if (quiet_require("jsonlite")) {
    jsonlite::write_json(manifest, manifest_json_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    esc <- function(x) gsub('"', '\\"', x)
    lines <- c(
      "{",
      paste0('  "run_name": "', esc(run_name), '",'),
      paste0('  "created": "', esc(manifest$created), '",'),
      paste0('  "output_dir": "', esc(manifest$output_dir), '",'),
      paste0('  "script_name": "', esc(script_id$script_name), '",'),
      paste0('  "script_path": "', esc(as.character(script_id$script_path)), '",'),
      paste0('  "git_branch": "', esc(as.character(git_meta$git_branch)), '",'),
      paste0('  "git_commit": "', esc(as.character(git_meta$git_commit)), '"'),
      "}"
    )
    write_text(manifest_json_path, lines)
  }
  
  inv <- inventory_outputs(out_dir)
  write_csv(inv, file.path(out_dir, "Project_Manifest_Files.csv"))
  invisible(list(manifest_json_path = manifest_json_path))
}

write_qmd_and_render <- function(out_dir, title, script_id, git_meta, key_outputs, notes = NULL) {
  qmd_path <- file.path(out_dir, "Provenance.qmd")
  html_path <- file.path(out_dir, "Provenance.html")
  
  body <- c(
    "---",
    paste0('title: "', title, '"'),
    "format: html",
    "execute:",
    "  echo: true",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "## Run provenance",
    "",
    paste0("- **Timestamp:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("- **Working directory:** ", script_id$wd),
    paste0("- **Script name:** ", script_id$script_name),
    paste0("- **Script path:** ", script_id$script_path),
    "",
    "## Git metadata (best-effort)",
    "",
    paste0("- **Branch:** ", git_meta$git_branch),
    paste0("- **Commit:** ", git_meta$git_commit),
    "",
    "## Outputs",
    ""
  )
  
  for (p in key_outputs) body <- c(body, paste0("- ", basename(p)))
  
  if (!is.null(notes) && length(notes) > 0) body <- c(body, "", "## Notes", "", notes)
  
  body <- c(body, "", "## Session info", "", "```{r}", "sessionInfo()", "```")
  
  write_text(qmd_path, body)
  
  quarto_ok <- nzchar(Sys.which("quarto"))
  if (quarto_ok) {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(out_dir)
    cmd <- paste0('quarto render "', basename(qmd_path), '"')
    suppressWarnings(suppressMessages(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)))
  } else {
    write_text(html_path, c(
      "<html><body>",
      "<h2>Quarto not found</h2>",
      "<p>Provenance.qmd was created, but Quarto is not available on PATH to render HTML.</p>",
      "</body></html>"
    ))
  }
  
  invisible(list(qmd_path = qmd_path, html_path = html_path))
}

# -----------------------------
# Helpers for ABC universe interrogation
# -----------------------------
must_char <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & nzchar(x)]
}

col_counts <- function(v) as.data.frame(table(v), stringsAsFactors = FALSE, responseName = "Observed", row.names = NULL)

design_diagnostics <- function(A, B, C) {
  tabA <- table(A); tabB <- table(B); tabC <- table(C)
  # Balanced A means all nonzero counts are identical
  A_vals <- as.integer(tabA)
  A_balanced <- (length(unique(A_vals)) == 1)
  
  out <- data.frame(
    metric = c("n_rows", "n_unique_A", "n_unique_B", "n_unique_C",
               "A_min_count", "A_max_count", "A_balanced"),
    value = c(length(A), length(tabA), length(tabB), length(tabC),
              min(A_vals), max(A_vals), A_balanced),
    stringsAsFactors = FALSE
  )
  out
}

predictor_polarity <- function(B, C) {
  tabB <- table(B); tabC <- table(C)
  genes <- sort(unique(c(names(tabB), names(tabC))))
  b <- as.integer(tabB[genes]); b[is.na(b)] <- 0L
  c <- as.integer(tabC[genes]); c[is.na(c)] <- 0L
  tot <- b + c
  
  # Polarity in [-1,1]
  pol <- ifelse(tot > 0, (b - c) / tot, NA_real_)
  
  # Exact binomial test under H0: P(B)=P(C)=0.5 conditional on appearing as predictor at all
  p_binom <- rep(NA_real_, length(genes))
  for (i in seq_along(genes)) {
    if (tot[i] > 0) p_binom[i] <- binom.test(b[i], tot[i], p = 0.5)$p.value
  }
  
  role <- ifelse(pol >= 0.2, "B_polar",
                 ifelse(pol <= -0.2, "C_polar", "neutral"))
  
  data.frame(
    Gene = genes,
    B_count = b,
    C_count = c,
    BC_total = tot,
    Polarity = pol,
    Binom_p = p_binom,
    Role = role,
    stringsAsFactors = FALSE
  )[order(-abs(pol), -tot), ]
}

conditional_matrix_long <- function(A, X, label = "B") {
  # Returns long table: A_gene, X_gene, count, prop_within_A, enrichment_vs_global
  tabAX <- table(A, X)
  A_genes <- rownames(tabAX)
  X_genes <- colnames(tabAX)
  
  # Global frequency of X across all rows
  tabX <- colSums(tabAX)
  p_global <- tabX / sum(tabX)
  
  out <- data.frame(
    A_gene = character(),
    X_gene = character(),
    Count = integer(),
    Prop_within_A = numeric(),
    Enrich_vs_global = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (a in A_genes) {
    counts <- tabAX[a, ]
    n <- sum(counts)
    if (n == 0) next
    props <- as.numeric(counts) / n
    enrich <- ifelse(p_global > 0, props / p_global, NA_real_)
    tmp <- data.frame(
      A_gene = rep(a, length(X_genes)),
      X_gene = X_genes,
      Count = as.integer(counts),
      Prop_within_A = props,
      Enrich_vs_global = as.numeric(enrich),
      stringsAsFactors = FALSE
    )
    out <- rbind(out, tmp)
  }
  
  # Keep label-consistent column name
  names(out)[names(out) == "X_gene"] <- paste0(label, "_gene")
  out
}

bc_pair_counts <- function(B, C) {
  bc <- paste0(B, "||", C)
  tab <- sort(table(bc), decreasing = TRUE)
  parts <- strsplit(names(tab), "\\|\\|", fixed = FALSE)
  B_gene <- vapply(parts, `[`, "", 1)
  C_gene <- vapply(parts, `[`, "", 2)
  data.frame(
    B_gene = B_gene,
    C_gene = C_gene,
    Count = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# MAIN
# -----------------------------
cat("\n============================================================\n")
cat("ABC Triplet Universe Diagnostics + Polarity + Conditional Selection\n")
cat("============================================================\n")

script_id <- capture_script_identity(verbose = TRUE)
git_meta  <- get_git_metadata()

abc_csv <- pick_file("Select the triplet CSV with columns A, B, C")
if (is.null(abc_csv) || !nzchar(abc_csv)) stop("No input file selected.")

df <- read.csv(abc_csv, check.names = FALSE, stringsAsFactors = FALSE)
req <- c("A", "B", "C")
missing <- setdiff(req, names(df))
if (length(missing) > 0) stop("Input CSV is missing required columns: ", paste(missing, collapse = ", "))

# Optional integration tables
cat("\nOptional: merge with annotation tables (must include a Gene column).\n")
cat("Cancel to skip.\n")

integration_paths <- character()
if (quiet_require("rstudioapi") && rstudioapi::isAvailable()) {
  integration_paths <- tryCatch(rstudioapi::selectFile(caption = "Select 0+ annotation CSVs to merge (Cancel to skip)", multi = TRUE),
                                error = function(e) character())
  if (is.null(integration_paths)) integration_paths <- character()
} else {
  cat("\nNon-RStudio mode: optionally select ONE integration CSV. Cancel to skip.\n")
  one <- tryCatch(file.choose(), error = function(e) "")
  if (nzchar(one)) integration_paths <- one
}

# Output dir
run_base <- tools::file_path_sans_ext(basename(abc_csv))
run_name <- safe_slug(paste0("ABC_universe_", run_base))
out_root <- file.path(getwd(), "outputs")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
out_dir <- file.path(out_root, paste0(run_name, "_", ts_stamp()))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
file.copy(abc_csv, file.path(out_dir, basename(abc_csv)), overwrite = TRUE)

cat("\nOutput folder:\n  ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")

# Extract columns (clean)
A <- must_char(df[["A"]])
B <- must_char(df[["B"]])
C <- must_char(df[["C"]])

# Ensure same length after filtering
n <- min(length(A), length(B), length(C))
A <- A[seq_len(n)]; B <- B[seq_len(n)]; C <- C[seq_len(n)]

# -----------------------------
# 1) Design diagnostics
# -----------------------------
diag <- design_diagnostics(A, B, C)
diag_path <- file.path(out_dir, "ABC_design_diagnostics.csv")
write_csv(diag, diag_path)

A_balanced <- as.logical(diag$value[diag$metric == "A_balanced"])

if (A_balanced) {
  cat("\nNOTE: Column A is BALANCED by construction (all genes appear equally often as outcomes).\n")
  cat("      Therefore, 'A enrichment' is not a meaningful statistic for this dataset.\n")
}

# -----------------------------
# 2) Predictor polarity (PRIMARY)
# -----------------------------
pol <- predictor_polarity(B, C)
pol_path <- file.path(out_dir, "ABC_predictor_polarity.csv")
write_csv(pol, pol_path)

# -----------------------------
# 3) Conditional selection given A: P(B|A) and P(C|A)
# -----------------------------
condB <- conditional_matrix_long(A, B, label = "B")
condC <- conditional_matrix_long(A, C, label = "C")

condB_path <- file.path(out_dir, "ABC_conditional_B_given_A.csv")
condC_path <- file.path(out_dir, "ABC_conditional_C_given_A.csv")
write_csv(condB[order(condB$A_gene, -condB$Enrich_vs_global, -condB$Count), ], condB_path)
write_csv(condC[order(condC$A_gene, -condC$Enrich_vs_global, -condC$Count), ], condC_path)

# -----------------------------
# 4) BC pair structure
# -----------------------------
bc_pairs <- bc_pair_counts(B, C)
bc_path <- file.path(out_dir, "ABC_BC_pair_counts.csv")
write_csv(bc_pairs, bc_path)

# -----------------------------
# 5) Role index (adapted)
#    Here A frequency is documented but not treated as an enrichment signal.
# -----------------------------
tabA <- table(A); tabB <- table(B); tabC <- table(C)
genes <- sort(unique(c(names(tabA), names(tabB), names(tabC))))
Obs_A <- as.integer(tabA[genes]); Obs_A[is.na(Obs_A)] <- 0L
Obs_B <- as.integer(tabB[genes]); Obs_B[is.na(Obs_B)] <- 0L
Obs_C <- as.integer(tabC[genes]); Obs_C[is.na(Obs_C)] <- 0L

role <- data.frame(Gene = genes, Obs_A = Obs_A, Obs_B = Obs_B, Obs_C = Obs_C, stringsAsFactors = FALSE)
role$Obs_total <- role$Obs_A + role$Obs_B + role$Obs_C

role$pA <- ifelse(role$Obs_total > 0, role$Obs_A / role$Obs_total, NA_real_)
role$pB <- ifelse(role$Obs_total > 0, role$Obs_B / role$Obs_total, NA_real_)
role$pC <- ifelse(role$Obs_total > 0, role$Obs_C / role$Obs_total, NA_real_)

# Polarity duplicated for convenience
pol_match <- match(role$Gene, pol$Gene)
role$B_count <- pol$B_count[pol_match]
role$C_count <- pol$C_count[pol_match]
role$Polarity <- pol$Polarity[pol_match]
role$Binom_p <- pol$Binom_p[pol_match]
role$PredictorRole <- pol$Role[pol_match]

role$DominantRole <- apply(role[, c("Obs_A","Obs_B","Obs_C")], 1, function(v) {
  lab <- c("A_outcome","B_numerator","C_denominator")
  lab[which.max(v)]
})

# Design-aware note column
role$A_is_balanced_design <- A_balanced

role <- role[order(-abs(role$Polarity), -role$Obs_total), ]

role_path <- file.path(out_dir, "ABC_role_index.csv")
write_csv(role, role_path)

# -----------------------------
# 6) OPTIONAL integrated table merge
# -----------------------------
integrated_path <- NULL
integrated <- NULL

normalize_gene_col <- function(d) {
  nms <- names(d)
  if ("Gene" %in% nms) return(d)
  hit <- which(tolower(nms) == "gene")
  if (length(hit) == 1) { names(d)[hit] <- "Gene"; return(d) }
  stop("Integration table lacks a 'Gene' column (case-insensitive). Columns are: ", paste(names(d), collapse = ", "))
}

if (length(integration_paths) > 0) {
  cat("\nMerging integration tables into ABC integrated gene table...\n")
  integrated <- role
  
  for (p in integration_paths) {
    if (!nzchar(p) || is.na(p)) next
    d2 <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
    d2 <- normalize_gene_col(d2)
    
    dup <- intersect(setdiff(names(d2), "Gene"), names(integrated))
    if (length(dup) > 0) {
      prefix <- safe_slug(tools::file_path_sans_ext(basename(p)))
      names(d2)[match(dup, names(d2))] <- paste0(prefix, "__", dup)
    }
    
    integrated <- merge(integrated, d2, by = "Gene", all = TRUE)
  }
  
  integrated <- integrated[order(-abs(integrated$Polarity), -integrated$Obs_total), ]
  integrated_path <- file.path(out_dir, "ABC_integrated_gene_table.csv")
  write_csv(integrated, integrated_path)
  
  for (p in integration_paths) {
    if (nzchar(p) && file.exists(p)) file.copy(p, file.path(out_dir, basename(p)), overwrite = TRUE)
  }
} else {
  cat("\nNo integration tables selected. Skipping integrated gene table.\n")
}

# -----------------------------
# README
# -----------------------------
readme_lines <- c(
  "ABC Triplet Universe Interrogation - Outputs",
  "-------------------------------------------",
  paste0("Run name   : ", run_name),
  paste0("Timestamp  : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Output dir : ", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
  "",
  "Core outputs:",
  paste0(" - ", basename(diag_path), " (Design diagnostics; detects balanced A)"),
  paste0(" - ", basename(pol_path), " (PRIMARY: predictor polarity B vs C per gene; binomial p-values)"),
  paste0(" - ", basename(condB_path), " (Conditional selection: P(B|A) + enrichment vs global B)"),
  paste0(" - ", basename(condC_path), " (Conditional selection: P(C|A) + enrichment vs global C)"),
  paste0(" - ", basename(bc_path), " (BC pair counts across the universe)"),
  paste0(" - ", basename(role_path), " (Convenience summary of counts + polarity + design flags)"),
  if (!is.null(integrated_path)) paste0(" - ", basename(integrated_path), " (Merged with your annotation tables)") else NULL,
  "",
  "Interpretation tips:",
  " - If A_balanced=TRUE, A enrichment is not meaningful because the outcome column was fixed by construction.",
  " - Use ABC_predictor_polarity.csv to classify genes as numerator-biased (B_polar) vs denominator-biased (C_polar).",
  " - Use conditional tables to see which predictors are preferentially selected for each outcome gene (P(B|A), P(C|A)).",
  " - Use BC pair counts to identify dominant regulatory ratio pairs in the full triplet universe.",
  "",
  "Provenance:",
  " - Provenance.qmd (+ Provenance.html if Quarto is available)",
  " - Project_Manifest.json + Project_Manifest_Files.csv"
)

readme_path <- file.path(out_dir, "README_outputs.txt")
write_text(readme_path, readme_lines)

# -----------------------------
# Provenance + Manifest
# -----------------------------
key_outputs <- c(diag_path, pol_path, condB_path, condC_path, bc_path, role_path, readme_path)
if (!is.null(integrated_path)) key_outputs <- c(key_outputs, integrated_path)

write_qmd_and_render(
  out_dir = out_dir,
  title = "ABC Triplet Universe — Provenance Report",
  script_id = script_id,
  git_meta  = git_meta,
  key_outputs = basename(key_outputs),
  notes = c(
    paste0("Input triplet file: ", basename(abc_csv)),
    paste0("A_balanced_design: ", A_balanced),
    if (length(integration_paths) > 0) paste0("Integration files: ", paste(basename(integration_paths), collapse = ", ")) else "Integration files: (none)"
  )
)

params <- list()
inputs <- list(
  triplet_csv = normalizePath(abc_csv, winslash = "/", mustWork = FALSE),
  integration_csvs = if (length(integration_paths) > 0) normalizePath(integration_paths, winslash = "/", mustWork = FALSE) else character()
)
write_manifest(out_dir, run_name, script_id, git_meta, inputs, params)

cat("\nDONE.\nKey outputs:\n")
for (p in key_outputs) cat(" - ", basename(p), "\n", sep = "")
cat("\n")