#!/usr/bin/env Rscript
# ============================================================
# Gene–Gene Correlation Edge Tables (per Day) from TIDY-LONG data
#   - INPUT (tidy long): Day, Gene, Expression
#   - Replicates inferred as row index within each Day×Gene (1..n)
#   - OUTPUT: long edge list for Cytoscape
#       Day, Gene1, Gene2, Cor, AbsCor, n_reps
#   - Base R only (NO dplyr)
#   - SOURCE-to-run in RStudio
# ============================================================

gene_gene_cor_by_day_tidy <- function() {
  
  has_rstudioapi <- requireNamespace("rstudioapi", quietly = TRUE)
  has_quarto    <- requireNamespace("quarto", quietly = TRUE)
  has_jsonlite  <- requireNamespace("jsonlite", quietly = TRUE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # ============================================================
  # Script identity capture (robust)
  # Captures script_name, script_path, script_full, timestamp
  # Works with SOURCE button, source(), and Rscript
  # ============================================================
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  get_script_path <- function() {
    
    script_full <- NA_character_
    
    # ---- Method 1: Rscript execution
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    fileArg <- grep("^--file=", cmdArgs, value = TRUE)
    
    if (length(fileArg) == 1) {
      script_full <- sub("^--file=", "", fileArg)
    }
    
    # ---- Method 2: RStudio editor
    if (is.na(script_full) && requireNamespace("rstudioapi", quietly = TRUE)) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
      if (nzchar(p)) script_full <- p
    }
    
    # ---- Method 3: source() call stack
    if (is.na(script_full)) {
      calls <- sys.calls()
      
      for (i in seq_along(calls)) {
        
        cl <- calls[[i]]
        
        if (length(cl) >= 2) {
          
          fn <- as.character(cl[[1]])
          
          if (fn %in% c("source", "sourceUTF8")) {
            
            if ("file" %in% names(as.list(cl))) {
              
              script_full <- tryCatch(
                eval(cl[["file"]], envir = parent.frame()),
                error = function(e) NA_character_
              )
              
            } else {
              
              script_full <- tryCatch(
                eval(cl[[2]], envir = parent.frame()),
                error = function(e) NA_character_
              )
              
            }
            
            if (!is.na(script_full)) break
          }
        }
      }
    }
    
    if (!is.na(script_full) && nzchar(script_full)) {
      script_full <- normalizePath(script_full, winslash = "/", mustWork = FALSE)
    }
    
    script_full
  }
  
  script_full <- get_script_path()
  script_path <- if (!is.na(script_full)) dirname(script_full) else NA_character_
  script_name <- if (!is.na(script_full)) basename(script_full) else NA_character_
  

  
  
  # ---- input
  message("\nINPUT expected (TIDY-LONG): Day/Group, Gene, Expression/Abundance")
  input_file <- NULL
  if (interactive() && has_rstudioapi) {
    input_file <- tryCatch(
      rstudioapi::selectFile(caption = "Select tidy LONG expression CSV", filter = "CSV (*.csv)"),
      error = function(e) NULL
    )
  }
  if (is.null(input_file) || !nzchar(input_file)) input_file <- tryCatch(file.choose(), error = function(e) "")
  if (!nzchar(input_file) || !file.exists(input_file)) stop("No valid input file selected.")
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  
  # ---- outputs
  if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
  run_name <- if (interactive()) readline("Enter a short run name (e.g., Liv5SD_corrNet): ") else ""
  if (!nzchar(run_name)) run_name <- "GeneGeneCor"
  
  safe_slug <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) x <- "run"
    x
  }
  
  output_dir <- file.path("outputs", paste0(safe_slug(run_name), "_", timestamp))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  message("\nOutput directory:\n", output_dir)
  
  # ---- read
  df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) < 3) stop("Input table too small.")
  
  message("\nDetected columns:")
  print(names(df))
  message("\nPreview (first 12 rows):")
  print(utils::head(df, 12))
  
  pick_one_col <- function(prompt, colnames_vec) {
    message("\n", prompt)
    for (i in seq_along(colnames_vec)) message(sprintf("  [%d] %s", i, colnames_vec[i]))
    ans <- if (interactive()) readline("Type the number: ") else ""
    idx <- suppressWarnings(as.integer(ans))
    if (is.na(idx) || idx < 1 || idx > length(colnames_vec)) stop("Invalid selection.")
    colnames_vec[idx]
  }
  
  cn <- names(df)
  day_col  <- pick_one_col("Which column is DAY / GROUP (e.g., 4,6,8,...)?", cn)
  gene_col <- pick_one_col("Which column is GENE SYMBOL?", cn)
  expr_col <- pick_one_col("Which column is EXPRESSION / ABUNDANCE (numeric)?", cn)
  
  dat <- df[, c(day_col, gene_col, expr_col)]
  names(dat) <- c("Day", "Gene", "Expr")
  
  dat$Day  <- as.character(dat$Day)
  dat$Gene <- as.character(dat$Gene)
  dat$Expr <- suppressWarnings(as.numeric(dat$Expr))
  
  keep <- !is.na(dat$Expr) & nzchar(dat$Day) & nzchar(dat$Gene)
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " row(s) with missing Day/Gene/Expr.")
    dat <- dat[keep, , drop = FALSE]
  }
  if (nrow(dat) < 3) stop("Too few valid rows after cleaning.")
  
  # ---- settings
  method <- "pearson"
  if (interactive()) {
    m <- tolower(trimws(readline("Choose method [pearson]: ")))
    if (nzchar(m)) method <- m
  }
  if (!method %in% c("pearson", "spearman")) stop("Method must be pearson or spearman.")
  
  abs_cutoff <- NA_real_
  if (interactive()) {
    message("\nOptional: filter edges by absolute correlation (|cor| >= cutoff).")
    ctf <- trimws(readline("Enter cutoff (blank for no filter): "))
    if (nzchar(ctf)) abs_cutoff <- suppressWarnings(as.numeric(ctf))
    if (!is.na(abs_cutoff) && (abs_cutoff < 0 || abs_cutoff > 1)) stop("Cutoff must be between 0 and 1.")
  }
  
  min_reps <- 3L
  if (interactive()) {
    mr <- trimws(readline("Minimum replicates required per Day [3]: "))
    if (nzchar(mr)) min_reps <- as.integer(mr)
  }
  if (is.na(min_reps) || min_reps < 3) min_reps <- 3L
  
  # ---- infer replicate index within each Day×Gene
  dat <- dat[order(dat$Day, dat$Gene), , drop = FALSE]
  dat$RepIndex <- ave(dat$Expr, dat$Day, dat$Gene, FUN = seq_along)
  
  # ---- helpers
  cor_to_edges <- function(C, day_label, n_reps) {
    g <- colnames(C)
    if (length(g) < 2) return(NULL)
    idx <- which(upper.tri(C), arr.ind = TRUE)
    data.frame(
      Day = day_label,
      Gene1 = g[idx[, 1]],
      Gene2 = g[idx[, 2]],
      Cor = C[idx],
      n_reps = n_reps,
      stringsAsFactors = FALSE
    )
  }
  
  # Build replicate×gene matrix using base R indexing
  build_matrix_for_day <- function(subd) {
    genes <- sort(unique(subd$Gene))
    reps  <- sort(unique(subd$RepIndex))
    mat <- matrix(NA_real_, nrow = length(reps), ncol = length(genes),
                  dimnames = list(paste0("rep", reps), genes))
    
    # Fill, averaging if duplicates exist for same (RepIndex, Gene)
    key <- paste(subd$RepIndex, subd$Gene, sep = "||")
    if (any(duplicated(key))) {
      agg <- tapply(subd$Expr, key, mean, na.rm = TRUE)
      ks <- names(agg)
      r  <- as.integer(sub("\\|\\|.*$", "", ks))
      g  <- sub("^.*\\|\\|", "", ks)
      for (i in seq_along(agg)) mat[paste0("rep", r[i]), g[i]] <- agg[[i]]
    } else {
      for (i in seq_len(nrow(subd))) {
        mat[paste0("rep", subd$RepIndex[i]), subd$Gene[i]] <- subd$Expr[i]
      }
    }
    mat
  }
  
  # ---- compute per day
  days <- sort(unique(dat$Day))
  corrmat_dir <- file.path(output_dir, "corr_matrices")
  dir.create(corrmat_dir, showWarnings = FALSE, recursive = TRUE)
  
  edges_all <- list()
  
  for (d in days) {
    subd <- dat[dat$Day == d, , drop = FALSE]
    reps <- sort(unique(subd$RepIndex))
    n_reps <- length(reps)
    
    if (n_reps < min_reps) {
      warning("Day ", d, " has only ", n_reps, " inferred replicates; skipping (need >= ", min_reps, ").")
      next
    }
    
    mat <- build_matrix_for_day(subd)
    
    # Drop genes with zero variance
    sds <- apply(mat, 2, sd, na.rm = TRUE)
    keep <- which(!is.na(sds) & sds > 0)
    if (length(keep) < 2) {
      warning("Day ", d, ": fewer than 2 genes with nonzero variance; skipping.")
      next
    }
    mat2 <- mat[, keep, drop = FALSE]
    
    C <- suppressWarnings(stats::cor(mat2, use = "pairwise.complete.obs", method = method))
    write.csv(C, file.path(corrmat_dir, paste0("CorrMatrix_Day", safe_slug(d), "_", method, ".csv")),
              row.names = TRUE)
    
    edges <- cor_to_edges(C, d, n_reps)
    if (!is.null(edges)) {
      if (!is.na(abs_cutoff)) edges <- edges[abs(edges$Cor) >= abs_cutoff, , drop = FALSE]
      edges_all[[d]] <- edges
    }
  }
  
  if (length(edges_all) == 0) stop("No edge tables generated. Check replicate counts and missing values.")
  
  edges_long <- do.call(rbind, edges_all)
  edges_long$AbsCor <- abs(edges_long$Cor)
  edges_long <- edges_long[order(edges_long$Day, -edges_long$AbsCor), , drop = FALSE]
  
  out_edges <- file.path(output_dir, paste0("GeneGene_Edges_LONG_byDay_", method, ".csv"))
  write.csv(edges_long, out_edges, row.names = FALSE)
  
  out_edges_filt <- NULL
  if (!is.na(abs_cutoff)) {
    out_edges_filt <- file.path(output_dir, paste0("GeneGene_Edges_LONG_byDay_", method, "_AbsCorGE_", abs_cutoff, ".csv"))
    write.csv(edges_long, out_edges_filt, row.names = FALSE)
  }
  
  # ---- manifest + inventory
  manifest <- list(
    script_name = script_name,
    script_path = script_path,
    script_full = script_full,
    timestamp = timestamp,
    input_file = input_file,
    output_dir = output_dir,
    chosen_columns = list(day_col = day_col, gene_col = gene_col, expr_col = expr_col),
    params = list(method = method, abs_cutoff = abs_cutoff, min_reps = min_reps),
    replicate_inference = "RepIndex = 1..n within each Day×Gene (order by Day,Gene)"
  )
  
  manifest_path <- file.path(output_dir, "Project_Manifest.json")
  if (has_jsonlite) {
    writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), manifest_path)
  } else {
    writeLines(capture.output(str(manifest)), manifest_path)
  }
  
  inventory_path <- file.path(output_dir, "Project_Manifest_Files.csv")
  inv <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  inv_df <- data.frame(
    file = basename(inv),
    rel_path = sub(paste0("^", gsub("\\\\", "/", output_dir), "/?"), "", gsub("\\\\", "/", inv)),
    abs_path = gsub("\\\\", "/", inv),
    stringsAsFactors = FALSE
  )
  write.csv(inv_df, inventory_path, row.names = FALSE)
  
  # ---- Quarto report (canonical contract; minimal)
  qmd_path <- file.path(output_dir, "GeneGene_Correlation_Report.qmd")
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  
  qmd <- c(
    "---",
    'title: "Gene–Gene Correlation Networks by Day (Tidy Long Input)"',
    "format:",
    "  html:",
    "    toc: true",
    "    toc-depth: 3",
    "execute:",
    "  echo: true",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "## Summary",
    "This analysis computes gene–gene correlation networks separately for each day.",
    "Input is tidy-long with Day, Gene, and Expression columns.",
    "Replicates are inferred as repeated measurements per Day×Gene (RepIndex = 1..n).",
    "",
    "## Metadata",
    "```",
    paste0("input_file: ", input_file),
    paste0("output_dir: ", output_dir),
    paste0("timestamp: ", timestamp),
    paste0("day_col: ", day_col),
    paste0("gene_col: ", gene_col),
    paste0("expr_col: ", expr_col),
    paste0("cor_method: ", method),
    paste0("abs_cutoff: ", ifelse(is.na(abs_cutoff), "none", abs_cutoff)),
    paste0("min_reps: ", min_reps),
    "```",
    "",
    "## Outputs",
    "```{r}",
    "list.files('.', recursive = TRUE)",
    "```",
    "",
    "## Edge table preview",
    "```{r}",
    paste0("edges <- read.csv('", basename(out_edges), "', stringsAsFactors = FALSE)"),
    "head(edges, 20)",
    "```",
    "",
    "## Session info",
    "```{r}",
    "sessionInfo()",
    "```"
  )
  writeLines(qmd, qmd_path)
  
  if (has_quarto) {
    old_wd <- getwd()
    setwd(output_dir)
    on.exit(setwd(old_wd), add = TRUE)
    tryCatch(quarto::quarto_render(basename(qmd_path)),
             error = function(e) warning("Quarto render failed: ", conditionMessage(e)))
  }
  
  res <- list(
    out_dir = output_dir,
    edges_long_csv = out_edges,
    edges_filtered_csv = out_edges_filt,
    corr_matrix_dir = corrmat_dir,
    manifest_json = manifest_path,
    inventory_csv = inventory_path,
    report_qmd = qmd_path,
    report_html = if (file.exists(html_path)) html_path else NA_character_
  )
  return(res)
}

# ------------------------------------------------------------
# AUTO-RUN WHEN SOURCED (RStudio / interactive use)
# ------------------------------------------------------------
if (interactive()) {
  message("Running gene_gene_cor_by_day_tidy() interactively...")
  res <- gene_gene_cor_by_day_tidy()
  message("Done. Results written to:")
  message(res$out_dir)
  message("Edge table:")
  message(res$edges_long_csv)
  if (!is.null(res$edges_filtered_csv)) message(res$edges_filtered_csv)
  if (!is.na(res$report_html)) message(paste("HTML report:", res$report_html))
  try(View(read.csv(res$edges_long_csv, stringsAsFactors = FALSE)), silent = TRUE)
}