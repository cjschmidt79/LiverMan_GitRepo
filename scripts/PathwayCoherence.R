# ======================================================================
# Pathway Coherence + Timing + Dispersion (Biologist-friendly, SOURCE to run)
# ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Small utilities (no fancy R)
# ----------------------------

.safe_message <- function(...) {
  msg <- paste0(...)
  message(msg)
}

.has_rstudio <- function() {
  isTRUE(Sys.getenv("RSTUDIO") == "1") || ("rstudioapi" %in% rownames(installed.packages()))
}

.pick_file <- function(prompt = "Select a file", filters = NULL) {
  if (.has_rstudio() && requireNamespace("rstudioapi", quietly = TRUE)) {
    f <- tryCatch(rstudioapi::selectFile(caption = prompt, filter = filters), error = function(e) NA_character_)
    if (!is.na(f) && nzchar(f)) return(f)
  }
  .safe_message(prompt)
  return(file.choose())
}

.pick_dir <- function(prompt = "Select a directory") {
  if (.has_rstudio() && requireNamespace("rstudioapi", quietly = TRUE)) {
    d <- tryCatch(rstudioapi::selectDirectory(caption = prompt), error = function(e) NA_character_)
    if (!is.na(d) && nzchar(d)) return(d)
  }
  .safe_message(prompt, " (fallback: using current working directory)")
  return(getwd())
}

.ask_text <- function(prompt, default = "") {
  ans <- readline(paste0(prompt, if (nzchar(default)) paste0(" [", default, "]") else "", ": "))
  if (!nzchar(ans)) ans <- default
  ans
}

.ask_yesno <- function(prompt, default = "n") {
  default <- tolower(default)
  ans <- tolower(readline(paste0(prompt, " [", default, "] (y/n): ")))
  if (!nzchar(ans)) ans <- default
  ans %in% c("y", "yes")
}

.choose_column <- function(dt, prompt, default_guess = NULL) {
  nms <- names(dt)
  default <- ""
  if (!is.null(default_guess) && default_guess %in% nms) default <- default_guess
  .safe_message(prompt)
  .safe_message("Available columns: ", paste(nms, collapse = ", "))
  col <- .ask_text("Type column name exactly", default = default)
  if (!col %in% nms) stop("Column not found: ", col)
  col
}

.normalize_abs <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.get_script_identity <- function() {
  script_full <- NA_character_
  if (.has_rstudio() && requireNamespace("rstudioapi", quietly = TRUE)) {
    script_full <- tryCatch({
      ctx <- rstudioapi::getSourceEditorContext()
      if (!is.null(ctx$path) && nzchar(ctx$path)) ctx$path else NA_character_
    }, error = function(e) NA_character_)
  }
  if (is.na(script_full)) {
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- "--file="
    hit <- grep(file_arg, ca, value = TRUE)
    if (length(hit) > 0) script_full <- sub(file_arg, "", hit[1])
  }
  script_full <- .normalize_abs(script_full)
  script_name <- if (!is.na(script_full) && nzchar(script_full)) basename(script_full) else NA_character_
  script_path <- if (!is.na(script_full) && nzchar(script_full)) dirname(script_full) else NA_character_
  list(script_name = script_name, script_path = script_path, script_full = script_full)
}

# ----------------------------
# Stats helpers
# ----------------------------

.safe_cor_stat <- function(mat, method = "pearson") {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (nrow(mat) < 2) return(NA_real_)
  
  v <- apply(mat, 1, function(x) stats::var(x, na.rm = TRUE))
  keep <- which(is.finite(v) & v > 0)
  if (length(keep) < 2) return(NA_real_)
  mat2 <- mat[keep, , drop = FALSE]
  
  C <- suppressWarnings(stats::cor(t(mat2), method = method, use = "pairwise.complete.obs"))
  if (!is.matrix(C) || any(dim(C) == 0)) return(NA_real_)
  
  vals <- abs(C[upper.tri(C)])
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

.z_against_random <- function(obs, null_vals) {
  mu <- mean(null_vals, na.rm = TRUE)
  sdv <- stats::sd(null_vals, na.rm = TRUE)
  if (!is.finite(obs) || !is.finite(mu) || !is.finite(sdv) || sdv == 0) {
    return(list(z = NA_real_, p = NA_real_, null_mean = mu, null_sd = sdv))
  }
  z <- (obs - mu) / sdv
  p <- mean(abs(null_vals - mu) >= abs(obs - mu), na.rm = TRUE)
  list(z = z, p = p, null_mean = mu, null_sd = sdv)
}

.bh_q <- function(p) {
  p2 <- p
  p2[!is.finite(p2)] <- NA_real_
  q <- rep(NA_real_, length(p2))
  ok <- which(is.finite(p2))
  if (length(ok) > 0) q[ok] <- p.adjust(p2[ok], method = "BH")
  q
}

.safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

# ----------------------------
# NEW: input expectations + validation
# ----------------------------

.print_expected_headers <- function(kind = c("expr", "map")) {
  kind <- match.arg(kind)
  if (kind == "expr") {
    .safe_message("\nEXPECTED EXPRESSION INPUT (tidy long):")
    .safe_message("  Required conceptual fields (any headers ok; you will map them):")
    .safe_message("    Day (integer day)")
    .safe_message("    SampleID (replicate/bird ID; must distinguish replicates within day)")
    .safe_message("    Gene (gene ID / symbol)")
    .safe_message("    Expr (numeric expression; TPM/FPKM/counts ok)")
    .safe_message("  Example headers: Day, SampleID, Gene, Expr\n")
  } else {
    .safe_message("\nEXPECTED MAPPING INPUT (Gene -> Pathway):")
    .safe_message("  Required conceptual fields (any headers ok; you will map them):")
    .safe_message("    Gene (must match expression Gene IDs)")
    .safe_message("    Pathway (function/program label)")
    .safe_message("  Example headers: Gene, Pathway\n")
  }
}

.validate_after_mapping <- function(dt, mp) {
  # dt must have Day, SampleID, Gene, Expr, Expr_log
  need_dt <- c("Day", "SampleID", "Gene", "Expr", "Expr_log")
  miss_dt <- setdiff(need_dt, names(dt))
  if (length(miss_dt) > 0) stop("Expression table is missing required internal columns: ", paste(miss_dt, collapse = ", "))
  
  need_mp <- c("Gene", "Pathway")
  miss_mp <- setdiff(need_mp, names(mp))
  if (length(miss_mp) > 0) stop("Mapping table is missing required internal columns: ", paste(miss_mp, collapse = ", "))
  
  if (!is.integer(dt$Day)) dt[, Day := as.integer(Day)]
  
  # Basic sanity: long tidy should have multiple rows per gene (across days and/or samples)
  if (nrow(dt) < 100) {
    .safe_message("WARNING: expression table has very few rows (", nrow(dt), "). Are you sure this is long tidy gene expression?")
  }
  
  # Overlap diagnostics
  g_expr <- unique(dt$Gene)
  g_map  <- unique(mp$Gene)
  overlap <- intersect(g_expr, g_map)
  
  .safe_message("\nGene overlap diagnostics:")
  .safe_message("  Unique genes in expression: ", length(g_expr))
  .safe_message("  Unique genes in mapping   : ", length(g_map))
  .safe_message("  Overlap genes             : ", length(overlap))
  
  if (length(overlap) == 0) {
    stop(
      "ZERO overlap between expression Gene IDs and mapping Gene IDs.\n",
      "This is the #1 reason coherence results are empty.\n",
      "Likely wrong mapping file OR Gene identifiers differ (e.g., Ensembl vs symbols)."
    )
  }
  invisible(TRUE)
}

# ----------------------------
# Main runner
# ----------------------------

run_pathway_coherence_timing <- function(
    expr_file = NULL,
    map_file  = NULL,
    base_out_dir = NULL,
    run_name = NULL,
    pseudocount = 1,
    corr_method = "pearson",
    day_window_k = 3,
    trans_window_k = 3,
    n_perm = 2000,
    seed = 1,
    min_genes_per_pathway = 5,
    top_n_paths_plot = 8,
    render_report = TRUE,
    open_report = FALSE
) {
  
  # ---------- Choose inputs ----------
  .print_expected_headers("expr")
  if (is.null(expr_file)) {
    expr_file <- .pick_file("Select tidy long expression CSV (Day, SampleID, Gene, Expr)")
  }
  .print_expected_headers("map")
  if (is.null(map_file)) {
    map_file <- .pick_file("Select Gene -> Pathway/Function mapping CSV")
  }
  expr_file <- .normalize_abs(expr_file)
  map_file  <- .normalize_abs(map_file)
  
  # ---------- Choose output directory + run naming ----------
  if (is.null(base_out_dir)) {
    base_out_dir <- file.path(getwd(), "outputs")
    if (.ask_yesno("Select a different base output directory?", default = "n")) {
      base_out_dir <- .pick_dir("Select base output directory (will create a subfolder inside)")
    }
  }
  base_out_dir <- .normalize_abs(base_out_dir)
  if (!dir.exists(base_out_dir)) dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(run_name)) {
    run_name <- .ask_text("Enter a short run name", default = "CoherenceTiming")
  }
  run_name <- .safe_filename(run_name)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  out_dir <- file.path(base_out_dir, paste0(run_name, "_", timestamp))
  fig_dir <- file.path(out_dir, "Figures")
  tab_dir <- file.path(out_dir, "Tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_dir_abs <- .normalize_abs(out_dir)
  .safe_message("Output directory:\n  ", out_dir_abs)
  writeLines(out_dir_abs, con = file.path(out_dir, "OUTPUT_DIR.txt"))
  
  # ---------- Provenance ----------
  sid <- .get_script_identity()
  
  # ---------- Load data ----------
  .safe_message("Reading expression file...")
  dt <- fread(expr_file)
  
  .safe_message("Reading mapping file...")
  mp <- fread(map_file)
  
  # ---------- Choose columns (interactive) ----------
  .safe_message("\nNow select columns for the EXPRESSION table.")
  day_col    <- .choose_column(dt, "Which column is Day?", default_guess = if ("Day" %in% names(dt)) "Day" else NULL)
  sample_col <- .choose_column(dt, "Which column is Sample ID / replicate?", default_guess = if ("SampleID" %in% names(dt)) "SampleID" else NULL)
  gene_col   <- .choose_column(dt, "Which column is Gene ID / symbol?", default_guess = if ("Gene" %in% names(dt)) "Gene" else NULL)
  expr_col   <- .choose_column(dt, "Which column is expression value (TPM/FPKM/etc)?", default_guess = if ("Expr" %in% names(dt)) "Expr" else NULL)
  
  .safe_message("\nNow select columns for the MAPPING table (Gene -> Pathway).")
  mgene_col  <- .choose_column(mp, "Which column is Gene ID (must match expression genes)?", default_guess = if ("Gene" %in% names(mp)) "Gene" else NULL)
  mpath_col  <- .choose_column(mp, "Which column is Pathway/Function group?", default_guess = if ("Pathway" %in% names(mp)) "Pathway" else NULL)
  
  # Standardize column names internally
  setnames(dt, c(day_col, sample_col, gene_col, expr_col), c("Day", "SampleID", "Gene", "Expr"))
  setnames(mp, c(mgene_col, mpath_col), c("Gene", "Pathway"))
  
  # Coerce types
  dt[, Day := as.integer(Day)]
  dt[, SampleID := as.character(SampleID)]
  dt[, Gene := as.character(Gene)]
  dt[, Expr := suppressWarnings(as.numeric(Expr))]
  
  mp[, Gene := as.character(Gene)]
  mp[, Pathway := as.character(Pathway)]
  
  # Basic filtering
  dt <- dt[is.finite(Day) & is.finite(Expr) & !is.na(Gene) & nzchar(Gene) & !is.na(SampleID) & nzchar(SampleID)]
  mp <- mp[!is.na(Gene) & nzchar(Gene) & !is.na(Pathway) & nzchar(Pathway)]
  
  # Transform
  dt[, Expr_log := log2(Expr + pseudocount)]
  
  # Validate overlap + structure (HARD STOP if mismatch)
  .validate_after_mapping(dt, mp)
  
  # ---------- Prepare gene universe ----------
  mean_dt <- dt[, .(Expr_log = mean(Expr_log, na.rm = TRUE)), by = .(Gene, Day)]
  days <- sort(unique(mean_dt$Day))
  if (length(days) < day_window_k) stop("Not enough unique days for day_window_k = ", day_window_k)
  
  mean_wide <- dcast(mean_dt, Gene ~ Day, value.var = "Expr_log")
  gene_universe <- mean_wide$Gene
  X <- as.matrix(mean_wide[, -1, with = FALSE])
  rownames(X) <- mean_wide$Gene
  colnames(X) <- as.character(sort(as.integer(colnames(X))))
  
  # Restrict mapping to genes observed
  mp <- mp[Gene %in% gene_universe]
  pathways <- unique(mp$Pathway)
  
  # Membership summary
  mem_sum <- mp[, .(n_genes = uniqueN(Gene), n_rows = .N), by = .(Pathway)][order(-n_genes)]
  fwrite(mem_sum, file.path(tab_dir, "PathwayMembership_Summary.csv"))
  
  # If nothing can pass the threshold, STOP EARLY with an informative message
  pass_paths <- mem_sum[n_genes >= min_genes_per_pathway]$Pathway
  if (length(pass_paths) == 0) {
    stop(
      "No pathways have >= min_genes_per_pathway after overlap filtering.\n",
      "min_genes_per_pathway = ", min_genes_per_pathway, "\n",
      "Check: (1) mapping gene IDs match expression, (2) mapping pathway sizes."
    )
  }
  
  # ---------- Day windows ----------
  day_windows <- lapply(seq_len(length(days) - day_window_k + 1), function(i) days[i:(i + day_window_k - 1)])
  names(day_windows) <- sapply(day_windows, function(w) paste0("d", w[1], "_to_d", w[length(w)]))
  
  # ---------- Option 1: Sliding-window coherence ----------
  set.seed(seed)
  .safe_message("\nRunning Option 1: sliding-window coherence (this can take a bit)...")
  win_results <- list()
  
  for (wname in names(day_windows)) {
    wdays <- day_windows[[wname]]
    cols <- match(as.character(wdays), colnames(X))
    Xw <- X[, cols, drop = FALSE]
    
    for (pw in pass_paths) {
      genes <- unique(mp[Pathway == pw, Gene])
      genes <- genes[genes %in% rownames(Xw)]
      if (length(genes) < min_genes_per_pathway) next
      
      obs <- .safe_cor_stat(Xw[genes, , drop = FALSE], method = corr_method)
      
      null_vals <- numeric(n_perm)
      for (b in seq_len(n_perm)) {
        g0 <- sample(rownames(Xw), size = length(genes), replace = FALSE)
        null_vals[b] <- .safe_cor_stat(Xw[g0, , drop = FALSE], method = corr_method)
      }
      
      zp <- .z_against_random(obs, null_vals)
      
      win_results[[length(win_results) + 1]] <- data.table(
        window = wname,
        start_day = wdays[1],
        end_day = wdays[length(wdays)],
        window_mid = (wdays[1] + wdays[length(wdays)]) / 2,
        Pathway = pw,
        n_genes = length(genes),
        stat_mean_abs_cor = obs,
        z = zp$z,
        p_emp = zp$p,
        null_mean = zp$null_mean,
        null_sd = zp$null_sd
      )
    }
    .safe_message("  finished window: ", wname)
  }
  
  # IMPORTANT FIX: handle empty results robustly
  if (length(win_results) == 0) {
    win_tab <- data.table()
    .safe_message("WARNING: No window coherence results were produced. Check mapping overlap and min_genes_per_pathway.")
  } else {
    win_tab <- rbindlist(win_results, fill = TRUE)
    if ("p_emp" %in% names(win_tab)) win_tab[, q := .bh_q(p_emp)] else win_tab[, q := NA_real_]
  }
  fwrite(win_tab, file.path(tab_dir, "WindowCoherence_ByPathway.csv"))
  
  # Plot heatmap/trajectories only if results exist
  if (nrow(win_tab) > 0) {
    hm <- copy(win_tab)
    hm[, window := factor(window, levels = names(day_windows))]
    ord <- hm[, .(maxz = max(z, na.rm = TRUE)), by = Pathway][order(-maxz)]$Pathway
    hm[, Pathway := factor(Pathway, levels = ord)]
    
    p_hm <- ggplot(hm, aes(x = window, y = Pathway, fill = z)) +
      geom_tile(color = "white") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(
        title = "Sliding-window pathway coherence (Z vs random gene sets)",
        subtitle = paste0("Coherence = mean |cor| (", corr_method, ") on log2(Expr + ", pseudocount, ")"),
        x = "Day window",
        y = "Pathway",
        fill = "Z"
      )
    
    ggsave(file.path(fig_dir, "WindowCoherence_Heatmap_Z.png"), p_hm, width = 12, height = 6, dpi = 600)
    ggsave(file.path(fig_dir, "WindowCoherence_Heatmap_Z.pdf"), p_hm, width = 12, height = 6)
    
    top_paths_win <- hm[, .(maxz = max(z, na.rm = TRUE)), by = Pathway][order(-maxz)][1:min(top_n_paths_plot, .N)]$Pathway
    tr <- hm[Pathway %in% top_paths_win]
    
    p_tr <- ggplot(tr, aes(x = window_mid, y = z, group = Pathway, color = Pathway)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      theme_bw(base_size = 12) +
      labs(
        title = "Top pathways: sliding-window coherence across development",
        x = "Developmental time (window midpoint day)",
        y = "Coherence Z",
        color = "Pathway"
      )
    
    ggsave(file.path(fig_dir, "WindowCoherence_TopTrajectories.png"), p_tr, width = 10, height = 6, dpi = 600)
    ggsave(file.path(fig_dir, "WindowCoherence_TopTrajectories.pdf"), p_tr, width = 10, height = 6)
  }
  
  # ---------- Option 2: Transition-window coherence ondelta profiles ----------
  .safe_message("\nRunning Option 2: transition-window coherence on delta profiles...")
  
  day_pairs <- data.table(
    d1 = days[-length(days)],
    d2 = days[-1]
  )
  trans_names <- paste0("d", day_pairs$d2, "_minus_d", day_pairs$d1)
  
  D <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(day_pairs))
  rownames(D) <- rownames(X)
  colnames(D) <- trans_names
  
  for (i in seq_len(nrow(day_pairs))) {
    c1 <- match(as.character(day_pairs$d1[i]), colnames(X))
    c2 <- match(as.character(day_pairs$d2[i]), colnames(X))
    D[, i] <- X[, c2] - X[, c1]
  }
  
  if (ncol(D) < trans_window_k) stop("Not enough transitions for trans_window_k = ", trans_window_k)
  
  t_windows <- lapply(seq_len(ncol(D) - trans_window_k + 1), function(i) seq(i, i + trans_window_k - 1))
  names(t_windows) <- sapply(t_windows, function(idx) paste0(trans_names[idx[1]], "_to_", trans_names[idx[length(idx)]]))
  
  tr_results <- list()
  
  for (tw in names(t_windows)) {
    idx <- t_windows[[tw]]
    Dw <- D[, idx, drop = FALSE]
    
    center_idx <- idx[ceiling(length(idx) / 2)]
    center_trans <- trans_names[center_idx]
    to_day <- suppressWarnings(as.integer(sub("^d([0-9]+)_minus_d[0-9]+$", "\\1", center_trans)))
    if (!is.finite(to_day)) to_day <- NA_integer_
    
    for (pw in pass_paths) {
      genes <- unique(mp[Pathway == pw, Gene])
      genes <- genes[genes %in% rownames(Dw)]
      if (length(genes) < min_genes_per_pathway) next
      
      obs <- .safe_cor_stat(Dw[genes, , drop = FALSE], method = corr_method)
      
      null_vals <- numeric(n_perm)
      for (b in seq_len(n_perm)) {
        g0 <- sample(rownames(Dw), size = length(genes), replace = FALSE)
        null_vals[b] <- .safe_cor_stat(Dw[g0, , drop = FALSE], method = corr_method)
      }
      
      zp <- .z_against_random(obs, null_vals)
      
      tr_results[[length(tr_results) + 1]] <- data.table(
        trans_window = tw,
        center_transition = center_trans,
        trans_time = to_day,
        Pathway = pw,
        n_genes = length(genes),
        stat_mean_abs_cor = obs,
        z = zp$z,
        p_emp = zp$p,
        null_mean = zp$null_mean,
        null_sd = zp$null_sd
      )
    }
    .safe_message("  finished transition-window: ", tw)
  }
  
  # IMPORTANT FIX: handle empty results robustly
  if (length(tr_results) == 0) {
    tr_tab <- data.table()
    .safe_message("WARNING: No transition-window coherence results were produced. Check mapping overlap and min_genes_per_pathway.")
  } else {
    tr_tab <- rbindlist(tr_results, fill = TRUE)
    if ("p_emp" %in% names(tr_tab)) tr_tab[, q := .bh_q(p_emp)] else tr_tab[, q := NA_real_]
  }
  fwrite(tr_tab, file.path(tab_dir, "TransitionWindowCoherence_ByPathway.csv"))
  
  if (nrow(tr_tab) > 0) {
    hm2 <- copy(tr_tab)
    hm2[, center_transition := factor(center_transition, levels = trans_names)]
    ord2 <- hm2[, .(maxz = max(z, na.rm = TRUE)), by = Pathway][order(-maxz)]$Pathway
    hm2[, Pathway := factor(Pathway, levels = ord2)]
    
    p_hm2 <- ggplot(hm2, aes(x = center_transition, y = Pathway, fill = z)) +
      geom_tile(color = "white") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(
        title = "Transition-window coherence on Delta profiles (Z vs random gene sets)",
        subtitle = "Each tile reflects coordinated switching across a small window of adjacent transitions",
        x = "Center transition",
        y = "Pathway",
        fill = "Z"
      )
    
    ggsave(file.path(fig_dir, "TransitionCoherence_Heatmap_Z.png"), p_hm2, width = 12, height = 6, dpi = 600)
    ggsave(file.path(fig_dir, "TransitionCoherence_Heatmap_Z.pdf"), p_hm2, width = 12, height = 6)
  }
  
  # ---------- Option 3: Per-day dispersion within pathways ----------
  .safe_message("\nRunning Option 3: per-day dispersion (MAD across replicates)...")
  
  disp_gene <- dt[, .(
    mad = mad(Expr_log, constant = 1, na.rm = TRUE),
    iqr = IQR(Expr_log, na.rm = TRUE),
    n_rep = sum(is.finite(Expr_log))
  ), by = .(Gene, Day)]
  
  disp_pw <- merge(disp_gene, mp, by = "Gene", allow.cartesian = TRUE)
  
  disp_tab <- disp_pw[, .(
    n_genes     = as.integer(uniqueN(Gene)),
    median_mad  = as.numeric(median(as.numeric(mad), na.rm = TRUE)),
    median_iqr  = as.numeric(median(as.numeric(iqr), na.rm = TRUE)),
    median_nrep = as.numeric(median(as.numeric(n_rep), na.rm = TRUE))
  ), by = .(Pathway = Pathway, Day)][order(Pathway, Day)]
  
  
  fwrite(disp_tab, file.path(tab_dir, "Dispersion_ByPathwayDay.csv"))
  
  if (nrow(disp_tab) > 0) {
    top_paths_by_size <- mem_sum[n_genes >= min_genes_per_pathway][1:min(top_n_paths_plot, .N)]$Pathway
    dplt <- disp_tab[Pathway %in% top_paths_by_size]
    
    p_disp <- ggplot(dplt, aes(x = Day, y = median_mad, group = Pathway, color = Pathway)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      theme_bw(base_size = 12) +
      labs(
        title = "Per-day dispersion within pathways (median MAD on log2 scale)",
        subtitle = "Lower MAD = tighter buffering (less bird-to-bird variability)",
        x = "Day",
        y = "Median MAD (log2 scale)",
        color = "Pathway"
      )
    
    ggsave(file.path(fig_dir, "Dispersion_Trajectories.png"), p_disp, width = 10, height = 6, dpi = 600)
    ggsave(file.path(fig_dir, "Dispersion_Trajectories.pdf"), p_disp, width = 10, height = 6)
  }
  
  # ---------- Aligned biology-forward panel ----------
  .safe_message("\nCreating aligned biology-forward panel...")
  
  top_by_win <- character()
  if (nrow(win_tab) > 0) {
    top_by_win <- win_tab[, .(maxz = max(z, na.rm = TRUE)), by = Pathway][order(-maxz)][1:min(top_n_paths_plot, .N)]$Pathway
  }
  top_by_tr <- character()
  if (nrow(tr_tab) > 0) {
    top_by_tr <- tr_tab[, .(maxz = max(z, na.rm = TRUE)), by = Pathway][order(-maxz)][1:min(top_n_paths_plot, .N)]$Pathway
  }
  plot_paths <- unique(c(top_by_win, top_by_tr))
  if (length(plot_paths) == 0) plot_paths <- unique(mp$Pathway)
  
  aligned_list <- list()
  
  if (nrow(win_tab) > 0) {
    aligned_list[[length(aligned_list) + 1]] <- win_tab[Pathway %in% plot_paths, .(
      Pathway, Time = window_mid, Metric = "Window coherence (Z)", Value = z
    )]
  }
  if (nrow(tr_tab) > 0) {
    aligned_list[[length(aligned_list) + 1]] <- tr_tab[Pathway %in% plot_paths & is.finite(trans_time), .(
      Pathway, Time = as.numeric(trans_time), Metric = "Transition coherence (Z)", Value = z
    )]
  }
  if (nrow(disp_tab) > 0) {
    aligned_list[[length(aligned_list) + 1]] <- disp_tab[Pathway %in% plot_paths, .(
      Pathway, Time = as.numeric(Day), Metric = "Dispersion (MAD)", Value = median_mad
    )]
  }
  
  aligned_long <- if (length(aligned_list) > 0) rbindlist(aligned_list, fill = TRUE) else data.table()
  
  if (nrow(aligned_long) > 0) {
    aligned_long[, Metric := factor(Metric, levels = c("Window coherence (Z)", "Transition coherence (Z)", "Dispersion (MAD)"))]
    
    ord_paths <- plot_paths
    if (nrow(win_tab) > 0) {
      ord_paths <- win_tab[Pathway %in% plot_paths, .(maxz = max(z, na.rm = TRUE)), by = Pathway][order(-maxz)]$Pathway
    } else if (nrow(mem_sum) > 0) {
      ord_paths <- mem_sum[Pathway %in% plot_paths][order(-n_genes)]$Pathway
    }
    aligned_long[, Pathway := factor(Pathway, levels = ord_paths)]
    
    fwrite(aligned_long, file.path(tab_dir, "AlignedPanel_Long.csv"))
    
    p_aligned <- ggplot(aligned_long, aes(x = Time, y = Value, group = Pathway, color = Pathway)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      facet_wrap(~ Metric, ncol = 1, scales = "free_y") +
      theme_bw(base_size = 12) +
      labs(
        title = "Aligned program structure across development",
        subtitle = "Same pathways: coherence (state) + coherence (switching) + dispersion (buffering) on one timeline",
        x = "Developmental time (day)",
        y = NULL,
        color = "Pathway"
      )
    
    ggsave(file.path(fig_dir, "AlignedPanel_CoherenceTransitionDispersion.png"), p_aligned, width = 11, height = 9, dpi = 600)
    ggsave(file.path(fig_dir, "AlignedPanel_CoherenceTransitionDispersion.pdf"), p_aligned, width = 11, height = 9)
  }
  
  # ---------- Manifest + Inventory ----------
  # ---------- Manifest (CSV, human-readable) ----------
  manifest_dt <- data.table(
    param = c(
      "timestamp",
      "output_dir",
      "expr_file",
      "map_file",
      "script_name",
      "script_path",
      "script_full",
      "pseudocount",
      "transform",
      "corr_method",
      "day_window_k",
      "trans_window_k",
      "n_perm",
      "seed",
      "min_genes_per_pathway",
      "R_version"
    ),
    value = c(
      timestamp,
      out_dir_abs,
      expr_file,
      map_file,
      sid$script_name,
      sid$script_path,
      sid$script_full,
      as.character(pseudocount),
      paste0("log2(Expr + ", pseudocount, ")"),
      corr_method,
      as.character(day_window_k),
      as.character(trans_window_k),
      as.character(n_perm),
      as.character(seed),
      as.character(min_genes_per_pathway),
      R.version.string
    )
  )
  fwrite(manifest_dt, file.path(out_dir, "Manifest.csv"))
  
  inv <- data.table(file = list.files(out_dir, recursive = TRUE, full.names = FALSE))
  fwrite(inv, file.path(out_dir, "Inventory.csv"))
  
  # ---------- Optional Quarto report (HTML) ----------
  if (isTRUE(render_report)) {
    .safe_message("\nWriting + rendering Quarto HTML report...")
    
    if (!requireNamespace("quarto", quietly = TRUE)) {
      .safe_message("Quarto not available (R package 'quarto' not installed). Skipping report.")
      .safe_message("Install with: install.packages('quarto')")
    } else if (!requireNamespace("knitr", quietly = TRUE)) {
      .safe_message("Package 'knitr' not installed (needed for include_graphics). Skipping report.")
      .safe_message("Install with: install.packages('knitr')")
    } else {
      
      qmd_path <- file.path(out_dir, "Pathway_Coherence_Timing_Report.qmd")
      
      qmd <- c(
        "---",
        "title: \"Pathway Coherence + Timing + Dispersion\"",
        "format:",
        "  html:",
        "    toc: true",
        "    toc-depth: 2",
        "---",
        "",
        "## What this report contains",
        "",
        "- Window coherence: pathway coordination within sliding day windows (state-level coherence).",
        "- Transition coherence: pathway coordination on Delta profiles across adjacent transitions (switching coherence).",
        "- Dispersion: within-pathway bird-to-bird variability per day (median MAD on log2 scale).",
        "",
        "## Output directory",
        "",
        "```{r}",
        "cat(readLines('OUTPUT_DIR.txt'), sep='\\n')",
        "```",
        "",
        "## Manifest",
        "",
        "```{r}",
        "meta <- read.csv('Manifest.csv', stringsAsFactors = FALSE)",
        "print(meta, row.names = FALSE)",
        "```",
        "",
        "## Key tables present",
        "",
        "```{r}",
        "tabs <- list.files('Tables', pattern='\\\\.csv$', full.names=TRUE)",
        "if (length(tabs) == 0) {",
        "  cat('No tables found in Tables/.\\n')",
        "} else {",
        "  cat(paste(basename(tabs), collapse='\\n'), '\\n')",
        "}",
        "```",
        "",
        "## Figures present",
        "",
        "```{r}",
        "figs <- list.files('Figures', pattern='\\\\.(png|pdf)$', full.names=TRUE, ignore.case=TRUE)",
        "if (length(figs) == 0) {",
        "  cat('No figures found in Figures/.\\n')",
        "} else {",
        "  for (f in figs) {",
        "    cat('\\n### ', basename(f), '\\n', sep='')",
        "    if (grepl('\\\\.png$', f, ignore.case=TRUE)) {",
        "      knitr::include_graphics(f)",
        "    } else {",
        "      cat('PDF saved: ', basename(f), '\\n', sep='')",
        "    }",
        "  }",
        "}",
        "```",
        "",
        "## Snapshot of core results (if available)",
        "",
        "### Window coherence (first 10 rows)",
        "",
        "```{r}",
        "f <- file.path('Tables', 'WindowCoherence_ByPathway.csv')",
        "if (file.exists(f)) {",
        "  x <- read.csv(f, stringsAsFactors=FALSE)",
        "  print(utils::head(x, 10), row.names = FALSE)",
        "} else {",
        "  cat('WindowCoherence_ByPathway.csv not found.\\n')",
        "}",
        "```",
        "",
        "### Transition coherence (first 10 rows)",
        "",
        "```{r}",
        "f <- file.path('Tables', 'TransitionWindowCoherence_ByPathway.csv')",
        "if (file.exists(f)) {",
        "  x <- read.csv(f, stringsAsFactors=FALSE)",
        "  print(utils::head(x, 10), row.names = FALSE)",
        "} else {",
        "  cat('TransitionWindowCoherence_ByPathway.csv not found.\\n')",
        "}",
        "```",
        "",
        "### Dispersion (first 10 rows)",
        "",
        "```{r}",
        "f <- file.path('Tables', 'Dispersion_ByPathwayDay.csv')",
        "if (file.exists(f)) {",
        "  x <- read.csv(f, stringsAsFactors=FALSE)",
        "  print(utils::head(x, 10), row.names = FALSE)",
        "} else {",
        "  cat('Dispersion_ByPathwayDay.csv not found.\\n')",
        "}",
        "```",
        "",
        "## Session info",
        "",
        "```{r}",
        "sessionInfo()",
        "```"
      )
      
      writeLines(qmd, con = qmd_path)
      
      # WD-safe render contract
      old_wd <- getwd()
      setwd(out_dir)
      on.exit(setwd(old_wd), add = TRUE)
      
      .safe_message("Rendering Quarto report (HTML)...")
      quarto::quarto_render(basename(qmd_path))
      
      html_path <- sub("\\.qmd$", ".html", qmd_path)
      .safe_message("Report written: ", .normalize_abs(html_path))
      
      if (isTRUE(open_report) && interactive()) {
        try(utils::browseURL(.normalize_abs(html_path)), silent = TRUE)
      }
    }
  }
  
  res <- list(
    out_dir = out_dir_abs,
    tables = list(
      window = win_tab,
      transition = tr_tab,
      dispersion = disp_tab,
      aligned_long = aligned_long,
      membership = mem_sum
    ),
    figures_dir = .normalize_abs(fig_dir),
    tables_dir = .normalize_abs(tab_dir)
  )
  return(res)
}

# ======================================================================
# AUTO-RUN WHEN SOURCED (RStudio / interactive use)
# ======================================================================
if (interactive()) {
  message("Running run_pathway_coherence_timing() interactively...")
  res <- run_pathway_coherence_timing()
  message("Done. Results written to:")
  message(res$out_dir)
  if (!is.null(res$tables$aligned_long) && nrow(res$tables$aligned_long) > 0) {
    try(View(res$tables$aligned_long), silent = TRUE)
  }
}
