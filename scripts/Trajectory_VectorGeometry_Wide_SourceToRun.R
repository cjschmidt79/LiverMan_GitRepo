###############################################################################
# Trajectory Vector Geometry (Wide CSV) — SOURCE-to-run (Biologist-first)
# ---------------------------------------------------------------------------
# What this does:
# 1) Reads a WIDE gene expression CSV (rows = birds/replicates; columns = genes).
# 2) Parses Day (and Rep) from SampleID = "DAY_REP" if Day column not provided.
# 3) Computes per-day state vectors in TWO spaces:
#    (A) Mean-expression space (engine trajectory)
#    (B) Residual-dispersion space (constraint/sensors trajectory)
#        - within-day gene variance across birds
#        - log-variance regressed on log-mean to remove mean–variance dependence
# 4) Evaluates trajectory structure across the ENTIRE experiment (Day 4–20 or any set):
#    Layer 1 (PRIMARY): Global axis projection (first day → last day)
#    Layer 2: Tangent field (local direction continuity) for EVERY adjacent interval
#    Layer 3: Piecewise regime vectors (user-defined boundaries; default 4|8|14|20)
# 5) Bias controls:
#    - Break ranking across ALL adjacent intervals (no privileged days)
#    - PCA replication space (K=20) with fixed basis fit on early window (default ≤14)
# 6) Uncertainty:
#    - Bird-within-day bootstrap CIs (primary)
# 7) Outputs:
#    - outputs/<run_name>_<timestamp>/
#      tables/, figures/, report/, logs/
#    - manifest (JSON if jsonlite available; else txt fallback)
#    - inventory of generated files
#    - Quarto QMD + HTML report (if Quarto found)
#
# Requirements:
# - CSV input only
# - Base R only (NO dplyr/tidyverse). Optional: rstudioapi, quarto, jsonlite
#
###############################################################################

trajectory_vector_geometry_wide <- function(input_file = NULL,
                                            run_name = NULL,
                                            make_quarto = TRUE,
                                            make_plots = TRUE,
                                            transform = c("log2p1","raw"),
                                            pseudocount = 1,
                                            pca_k = 20,
                                            pca_train_max_day = 14,
                                            bootstrap_B = 500,
                                            permute_P = 0,      # set >0 to compute permutation p-values for cosines
                                            eps_var = 1e-8,
                                            verbose = TRUE) {
  # ------------------------------- helpers -----------------------------------
  say <- function(...) if (isTRUE(verbose)) message(...)
  now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
  
  is_rstudio <- function() {
    if (!requireNamespace("rstudioapi", quietly = TRUE)) return(FALSE)
    isTRUE(rstudioapi::isAvailable())
  }
  
  pick_file <- function() {
    if (is_rstudio()) {
      f <- try(rstudioapi::selectFile(caption = "Select WIDE expression CSV"), silent = TRUE)
      if (!inherits(f, "try-error") && nzchar(f)) return(f)
    }
    return(file.choose())
  }
  
  pick_dir <- function() {
    if (is_rstudio()) {
      d <- try(rstudioapi::selectDirectory(caption = "Select output parent directory (or Cancel to use ./outputs)"),
               silent = TRUE)
      if (!inherits(d, "try-error") && nzchar(d)) return(d)
    }
    return("") # caller will handle fallback
  }
  
  menu_pick <- function(prompt, choices, allow_none = FALSE) {
    say(prompt)
    if (allow_none) choices2 <- c("NONE", choices) else choices2 <- choices
    idx <- utils::menu(choices2, title = prompt)
    if (idx <= 0) return(NA_character_)
    return(choices2[idx])
  }
  
  norm_path <- function(p) {
    tryCatch(normalizePath(p, winslash = "/", mustWork = FALSE), error = function(e) p)
  }
  
  safe_mkdir <- function(p) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  
  # script identity (best-effort)
  get_script_identity <- function() {
    script_full <- NA_character_
    script_path <- NA_character_
    script_name <- NA_character_
    
    # common rstudio pattern
    if (is_rstudio()) {
      sp <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
      if (!inherits(sp, "try-error") && nzchar(sp)) script_full <- sp
    }
    
    # commandArgs fallback
    if (!nzchar(script_full) || is.na(script_full)) {
      ca <- commandArgs(trailingOnly = FALSE)
      file_arg <- "--file="
      hit <- grep(file_arg, ca, value = TRUE)
      if (length(hit) == 1) {
        script_full <- sub(file_arg, "", hit, fixed = TRUE)
      }
    }
    
    if (nzchar(script_full) && !is.na(script_full)) {
      script_full <- norm_path(script_full)
      script_path <- dirname(script_full)
      script_name <- basename(script_full)
    }
    
    list(script_name = script_name, script_path = script_path, script_full = script_full)
  }
  
  # vector math
  dot <- function(a,b) sum(a*b)
  vnorm <- function(a) sqrt(sum(a*a))
  cosine <- function(a,b) {
    na <- vnorm(a); nb <- vnorm(b)
    if (!is.finite(na) || !is.finite(nb) || na < 1e-12 || nb < 1e-12) return(NA_real_)
    dot(a,b)/(na*nb)
  }
  
  # write manifest JSON if jsonlite available; else simple key=value text
  write_manifest <- function(path, manifest_list) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      txt <- jsonlite::toJSON(manifest_list, pretty = TRUE, auto_unbox = TRUE, null = "null")
      writeLines(txt, con = path)
    } else {
      # fallback
      con <- file(path, open = "wt")
      on.exit(close(con), add = TRUE)
      writeLines("MANIFEST_FALLBACK_TXT (jsonlite not installed)\n", con)
      for (nm in names(manifest_list)) {
        val <- manifest_list[[nm]]
        if (is.null(val)) val <- "NULL"
        if (length(val) > 1) val <- paste(val, collapse = "; ")
        writeLines(paste0(nm, " = ", val), con)
      }
    }
  }
  
  # -------------------------- interactive inputs -----------------------------
  if (is.null(input_file) || !nzchar(input_file)) {
    say("📄 Select input WIDE CSV...")
    input_file <- pick_file()
  }
  if (!file.exists(input_file)) stop("Input file not found: ", input_file)
  input_file_abs <- norm_path(input_file)
  
  transform <- match.arg(transform)
  
  # output base directory
  say("📁 Select OUTPUT parent folder (the script will create ./outputs/<run>_<timestamp>/ under it).")
  say("   Tip: Cancel the picker to use the default: ./outputs/ under your current working directory.")
  out_parent <- pick_dir()
  
  if (!nzchar(out_parent)) {
    out_parent <- file.path(getwd(), "outputs")
    say("📁 No folder selected — using default output parent: ", norm_path(out_parent))
  } else {
    say("📁 Output parent selected: ", norm_path(out_parent))
  }
  
  out_parent <- norm_path(out_parent)
  safe_mkdir(out_parent)

  
  if (is.null(run_name) || !nzchar(run_name)) {
    run_name <- readline("Enter run name (short, no spaces; e.g., VecGeom_Liver5SD): ")
    if (!nzchar(run_name)) run_name <- "VecGeom"
  }
  
  stamp <- now_stamp()
  out_dir <- norm_path(file.path(out_parent, paste0(run_name, "_", stamp)))
  safe_mkdir(out_dir)
  
  # subfolders
  dir_tables <- file.path(out_dir, "tables")
  dir_figs   <- file.path(out_dir, "figures")
  dir_report <- file.path(out_dir, "report")
  dir_logs   <- file.path(out_dir, "logs")
  safe_mkdir(dir_tables); safe_mkdir(dir_figs); safe_mkdir(dir_report); safe_mkdir(dir_logs)
  
  say("✅ Output directory:")
  say(out_dir)
  
  # logging
  log_file <- file.path(dir_logs, paste0("LOG_", stamp, ".txt"))
  sink(log_file, split = TRUE)
  on.exit({
    try(sink(), silent = TRUE)
  }, add = TRUE)
  
  say("===== Trajectory Vector Geometry START =====")
  say("Input: ", input_file_abs)
  say("Transform: ", transform)
  
  # ------------------------------- read data ---------------------------------
  say("📥 Reading CSV...")
  df <- try(utils::read.csv(input_file_abs, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(df, "try-error")) stop("Failed to read CSV. Is it a valid comma-delimited CSV?")
  
  if (nrow(df) < 3) stop("Input has too few rows to analyze.")
  if (ncol(df) < 5) stop("Input has too few columns to be a wide expression table.")
  
  cn <- colnames(df)
  
  # select columns
  sample_col <- menu_pick("Select SAMPLE ID column (required; e.g., 14_3):", cn, allow_none = FALSE)
  if (is.na(sample_col)) stop("No SAMPLE ID column selected.")
  
  day_col <- menu_pick("Select DAY column (optional). Choose NONE to parse Day from SampleID:", cn, allow_none = TRUE)
  if (identical(day_col, "NONE") || is.na(day_col)) day_col <- NA_character_
  
  group_col <- menu_pick("Select optional GROUP/SOURCE column (optional). Choose NONE if not applicable:", cn, allow_none = TRUE)
  if (identical(group_col, "NONE") || is.na(group_col)) group_col <- NA_character_
  
  # parse Day/Rep from SampleID
  sample_id <- as.character(df[[sample_col]])
  if (any(!nzchar(sample_id))) stop("SampleID column contains empty values.")
  
  parse_day_rep <- function(x) {
    # deterministic: Day = substring before first underscore
    u <- strsplit(x, "_", fixed = TRUE)
    day <- vapply(u, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1))
    rep <- vapply(u, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
    list(day = day, rep = rep)
  }
  
  parsed <- parse_day_rep(sample_id)
  day_parsed <- suppressWarnings(as.numeric(parsed$day))
  rep_parsed <- suppressWarnings(as.integer(parsed$rep))
  
  if (any(!is.finite(day_parsed))) {
    stop("Failed to parse numeric Day from SampleID. Expected DAY_REP like 14_3.")
  }
  if (any(is.na(rep_parsed))) {
    say("⚠️ Replicate index (after underscore) could not be parsed for some rows. This is OK unless you need Rep explicitly.")
  }
  
  if (!is.na(day_col)) {
    day_given <- suppressWarnings(as.numeric(df[[day_col]]))
    if (any(!is.finite(day_given))) stop("Provided Day column has non-numeric values.")
    # validate match
    mismatch <- which(day_given != day_parsed)
    if (length(mismatch) > 0) {
      stop("Day column does not match Day parsed from SampleID for at least one row. Fix naming or choose NONE for Day column.")
    }
    day <- day_given
  } else {
    day <- day_parsed
  }
  
  df$.__DAY__ <- day
  df$.__SAMPLE__ <- sample_id
  if (!is.na(group_col)) df$.__GROUP__ <- as.character(df[[group_col]]) else df$.__GROUP__ <- "ALL"
  
  # choose group handling
  groups <- sort(unique(df$.__GROUP__))
  run_mode <- "pooled"
  if (length(groups) > 1) {
    ans <- readline(paste0("Multiple groups detected (", paste(groups, collapse = ", "),
                           "). Run (1) pooled ALL or (2) each group separately? [default 1]: "))
    if (nzchar(ans) && trimws(ans) == "2") run_mode <- "per_group"
  }
  
  # detect gene columns: numeric columns excluding id cols
  exclude_cols <- c(sample_col, day_col, group_col, ".__DAY__", ".__SAMPLE__", ".__GROUP__")
  keep_idx <- !(cn %in% exclude_cols)
  cand_cols <- cn[keep_idx]
  
  # numeric check
  is_num <- vapply(df[cand_cols], function(x) is.numeric(x) || is.integer(x), logical(1))
  gene_cols <- cand_cols[is_num]
  nonnum <- cand_cols[!is_num]
  
  if (length(gene_cols) < 10) stop("Fewer than 10 numeric gene columns detected. Check input format.")
  if (length(nonnum) > 0) {
    say("⚠️ Dropping non-numeric candidate gene columns (first 20 shown): ",
        paste(head(nonnum, 20), collapse = ", "))
  }
  
  # extract expression matrix
  X <- as.matrix(df[gene_cols])
  storage.mode(X) <- "double"
  
  # missing values handling
  na_frac <- sum(!is.finite(X)) / length(X)
  if (na_frac > 0) say(sprintf("⚠️ Non-finite values in expression matrix: %.3f%%", 100*na_frac))
  if (na_frac > 0.10) say("⚠️ >10% non-finite values. Downstream results may be unreliable.")
  
  # transform
  if (transform == "log2p1") {
    say("🔁 Applying log2(x + pseudocount) transform with pseudocount=", pseudocount)
    X <- log2(X + pseudocount)
  } else {
    say("🔁 Using raw expression values (no transform).")
  }
  
  # --------------------- core computation functions --------------------------
  compute_day_vectors <- function(Xmat, day_vec, group_vec, group_value,
                                  eps_var = 1e-8) {
    # returns:
    # - days_sorted
    # - mean_by_day (list named by day: numeric vector length G)
    # - disp_resid_by_day (list named by day: numeric vector length G) (may be NA if insufficient reps)
    # - mu_by_day, logvar_by_day (optional debug)
    idxg <- which(group_vec == group_value)
    if (length(idxg) < 3) stop("Group has too few samples: ", group_value)
    
    days <- sort(unique(day_vec[idxg]))
    days <- days[is.finite(days)]
    if (length(days) < 3) stop("Need at least 3 days to analyze for group: ", group_value)
    
    # report missing expected day grid if looks like 4..20 even spacing
    # (non-fatal, just informative)
    day_report <- paste(days, collapse = ", ")
    say("Days in analysis: ", day_report)
    
    mean_list <- list()
    disp_list <- list()
    mu_list <- list()
    logvar_list <- list()
    nrep <- integer(length(days))
    names(nrep) <- as.character(days)
    
    for (d in days) {
      idxd <- idxg[day_vec[idxg] == d]
      nrep[as.character(d)] <- length(idxd)
      if (length(idxd) < 1) next
      Xd <- Xmat[idxd, , drop = FALSE]
      # day mean
      m <- colMeans(Xd, na.rm = TRUE)
      mean_list[[as.character(d)]] <- m
      
      # within-day variance for dispersion
      if (nrow(Xd) >= 2) {
        mu <- m
        vv <- apply(Xd, 2, stats::var, na.rm = TRUE)
        vv[!is.finite(vv)] <- NA_real_
        logvar <- log(vv + eps_var)
        logmu <- log(abs(mu) + eps_var)
        
        # residualize logvar ~ logmu across genes (within day)
        ok <- is.finite(logvar) & is.finite(logmu)
        resid <- rep(NA_real_, length(logvar))
        if (sum(ok) >= 10) {
          fit <- stats::lm(logvar[ok] ~ logmu[ok])
          resid[ok] <- stats::residuals(fit)
        }
        disp_list[[as.character(d)]] <- resid
        mu_list[[as.character(d)]] <- mu
        logvar_list[[as.character(d)]] <- logvar
      } else {
        # insufficient reps for variance
        disp_list[[as.character(d)]] <- rep(NA_real_, ncol(Xmat))
        mu_list[[as.character(d)]] <- m
        logvar_list[[as.character(d)]] <- rep(NA_real_, ncol(Xmat))
      }
    }
    
    list(days = days,
         nrep = nrep,
         mean_by_day = mean_list,
         disp_resid_by_day = disp_list,
         mu_by_day = mu_list,
         logvar_by_day = logvar_list)
  }
  
  eval_layers <- function(days, x_by_day, label, pca_k = 20,
                          pca_train_max_day = 14,
                          bootstrap_CI = NULL,
                          permute_P = 0,
                          out_prefix = "MEAN",
                          dir_tables = ".", dir_figs = ".") {
    # days: numeric sorted
    # x_by_day: list named by day, each numeric vector length G
    # Returns list with tables, plus optional PCA replication results.
    
    # build day matrix (D x G)
    D <- length(days)
    day_names <- as.character(days)
    G <- length(x_by_day[[day_names[1]]])
    
    Xday <- matrix(NA_real_, nrow = D, ncol = G)
    for (i in seq_along(days)) {
      Xday[i, ] <- x_by_day[[as.character(days[i])]]
    }
    
    # bug detector: too many NAs?
    na_frac <- sum(!is.finite(Xday)) / length(Xday)
    if (na_frac > 0.10) {
      say(sprintf("⚠️ %s day-matrix has %.2f%% non-finite entries.", label, 100*na_frac))
    }
    
    # adjacent intervals
    intervals <- data.frame(
      day_from = days[-length(days)],
      day_to   = days[-1],
      stringsAsFactors = FALSE
    )
    
    # global axis (first -> last)
    v_global <- Xday[D, ] - Xday[1, ]
    vgn <- vnorm(v_global)
    if (!is.finite(vgn) || vgn < 1e-10) stop(label, ": Global axis norm is near zero (no measurable change).")
    vhat <- v_global / vgn
    
    # compute per-interval metrics
    n_int <- nrow(intervals)
    cos_to_global <- rep(NA_real_, n_int)
    along_s <- rep(NA_real_, n_int)
    orth_norm <- rep(NA_real_, n_int)
    
    for (i in seq_len(n_int)) {
      d1 <- intervals$day_from[i]
      d2 <- intervals$day_to[i]
      x1 <- x_by_day[[as.character(d1)]]
      x2 <- x_by_day[[as.character(d2)]]
      delta <- x2 - x1
      cval <- cosine(delta, v_global)
      s <- dot(delta, vhat)
      perp <- delta - s*vhat
      cos_to_global[i] <- cval
      along_s[i] <- s
      orth_norm[i] <- vnorm(perp)
    }
    
    # tangent continuity across adjacent boundaries
    # tangent at index i uses i-1 and i+1: v_i = x_{i+1} - x_{i-1}
    tangents <- matrix(NA_real_, nrow = D, ncol = G)
    for (i in 2:(D-1)) {
      tangents[i, ] <- Xday[i+1, ] - Xday[i-1, ]
    }
    tangent_continuity <- rep(NA_real_, n_int)
    # continuity for boundary i (between day i and i+1) uses tangents at i and i+1
    # valid for i=2..D-2
    for (i in 2:(D-2)) {
      tangent_continuity[i] <- cosine(tangents[i, ], tangents[i+1, ])
    }
    
    # permutation p-values for cosines to global (gene-label permutation)
    perm_p <- rep(NA_real_, n_int)
    if (permute_P > 0) {
      say("🔀 Permutation test (gene-label) for cosine alignment to global axis: P=", permute_P)
      for (i in seq_len(n_int)) {
        d1 <- intervals$day_from[i]; d2 <- intervals$day_to[i]
        delta <- x_by_day[[as.character(d2)]] - x_by_day[[as.character(d1)]]
        obs <- cosine(delta, v_global)
        if (!is.finite(obs)) next
        ge <- 0
        for (p in seq_len(permute_P)) {
          perm <- sample.int(length(delta))
          cperm <- cosine(delta[perm], v_global)
          if (is.finite(cperm) && cperm >= obs) ge <- ge + 1
        }
        perm_p[i] <- (ge + 1) / (permute_P + 1)
      }
    }
    
    tab <- data.frame(
      day_from = intervals$day_from,
      day_to = intervals$day_to,
      cosine_to_global = cos_to_global,
      along_s = along_s,
      orth_norm = orth_norm,
      tangent_continuity = tangent_continuity,
      perm_p_cosine_ge = perm_p,
      stringsAsFactors = FALSE
    )
    
    # Rank "break-likeness" across all adjacent intervals
    # Break score = high orth_norm + low tangent_continuity (where defined)
    # We compute ranks/percentiles for transparency (anti-bias)
    # For tangent continuity, lower is more "break"; for orth_norm, higher is more "break"
    orth_rank <- rank(-tab$orth_norm, ties.method = "average")  # descending orth_norm
    orth_pct  <- (orth_rank - 1) / (length(orth_rank) - 1)
    # tangent only for defined indices; others NA
    tang_rank <- rep(NA_real_, nrow(tab))
    tang_pct <- rep(NA_real_, nrow(tab))
    ok_t <- is.finite(tab$tangent_continuity)
    if (sum(ok_t) >= 2) {
      r <- rank(tab$tangent_continuity[ok_t], ties.method = "average") # ascending continuity (low = break)
      # convert to "breakiness" percentile: low continuity => high break percentile
      # so break_pct = 1 - ((rank-1)/(n-1))
      tang_break_pct <- 1 - ((r - 1) / (length(r) - 1))
      tang_rank[ok_t] <- r
      tang_pct[ok_t] <- tang_break_pct
    }
    
    tab$orth_break_percentile <- orth_pct
    tab$tangent_break_percentile <- tang_pct
    
    # PCA replication (K=20) with fixed basis fit on days <= train_max_day
    # (if insufficient early days, fallback to all days)
    pca_res <- NULL
    scores_all <- NULL
    try({
      train_idx <- which(days <= pca_train_max_day)
      if (length(train_idx) < 3) {
        say("⚠️ Not enough days <= pca_train_max_day for PCA training; using all days.")
        train_idx <- seq_along(days)
      }
      Xtrain <- Xday[train_idx, , drop = FALSE]
      # center genes; do not scale by default
      pca <- stats::prcomp(Xtrain, center = TRUE, scale. = FALSE)
      k_eff <- min(pca_k, ncol(pca$x), ncol(pca$rotation))
      # project ALL days onto training basis
      # Use training center:
      center_vec <- pca$center
      Xc <- sweep(Xday, 2, center_vec, FUN = "-")
      scores <- Xc %*% pca$rotation[, 1:k_eff, drop = FALSE]
      scores_all <- scores
      
      # compute same interval metrics in PCA space
      vG <- scores[D, ] - scores[1, ]
      vgn2 <- vnorm(vG)
      if (is.finite(vgn2) && vgn2 > 1e-12) {
        vhat2 <- vG / vgn2
        cos2 <- along2 <- orth2 <- rep(NA_real_, n_int)
        for (i in seq_len(n_int)) {
          di <- i
          delta2 <- scores[di+1, ] - scores[di, ]
          cos2[i] <- cosine(delta2, vG)
          s2 <- dot(delta2, vhat2)
          perp2 <- delta2 - s2*vhat2
          along2[i] <- s2
          orth2[i] <- vnorm(perp2)
        }
        pca_tab <- data.frame(
          day_from = intervals$day_from,
          day_to = intervals$day_to,
          cosine_to_global_PCA = cos2,
          along_s_PCA = along2,
          orth_norm_PCA = orth2,
          stringsAsFactors = FALSE
        )
        pca_res <- list(k_eff = k_eff, train_max_day = pca_train_max_day, table = pca_tab, scores = scores)
      }
    }, silent = TRUE)
    
    # merge bootstrap CIs if provided
    if (!is.null(bootstrap_CI)) {
      # bootstrap_CI is a list with ci tables keyed by metric name
      # expects columns: day_from, day_to, metric, ci_low, ci_high
      tab <- merge(tab, bootstrap_CI, by = c("day_from","day_to"), all.x = TRUE, sort = FALSE)
    }
    
    # write tables
    out_csv <- file.path(dir_tables, paste0("AdjacentIntervals_VectorMetrics_", out_prefix, ".csv"))
    utils::write.csv(tab, out_csv, row.names = FALSE)
    
    # PCA table
    out_pca_csv <- NA_character_
    if (!is.null(pca_res) && !is.null(pca_res$table)) {
      out_pca_csv <- file.path(dir_tables, paste0("AdjacentIntervals_VectorMetrics_", out_prefix, "_PCAK", pca_res$k_eff, ".csv"))
      utils::write.csv(pca_res$table, out_pca_csv, row.names = FALSE)
    }
    
    # plots
    plot_files <- list()
    if (isTRUE(make_plots)) {
      # 1) Along-axis vs Orthogonal by interval
      f1 <- file.path(dir_figs, paste0("Intervals_AlongVsOrth_", out_prefix, ".png"))
      grDevices::png(f1, width = 1600, height = 900, res = 200)
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mar = c(7,5,4,2))
      xlab <- paste0(tab$day_from, "→", tab$day_to)
      plot(tab$along_s, type = "b", pch = 16, xaxt = "n",
           xlab = "", ylab = "Value", main = paste0(label, ": Along-axis progress vs Orthogonal excursion"))
      axis(1, at = seq_along(xlab), labels = xlab, las = 2, cex.axis = 0.8)
      lines(tab$orth_norm, type = "b", pch = 1)
      legend("topright", legend = c("along_s (projection)", "orth_norm (perp magnitude)"),
             bty = "n", pch = c(16,1))
      grDevices::dev.off()
      plot_files$along_orth <- f1
      
      # 2) Tangent continuity
      f2 <- file.path(dir_figs, paste0("Intervals_TangentContinuity_", out_prefix, ".png"))
      grDevices::png(f2, width = 1600, height = 900, res = 200)
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mar = c(7,5,4,2))
      y <- tab$tangent_continuity
      plot(y, type = "b", pch = 16, xaxt = "n", ylim = c(-1,1),
           xlab = "", ylab = "cosine(tangent_i, tangent_{i+1})",
           main = paste0(label, ": Tangent continuity across adjacent intervals"))
      axis(1, at = seq_along(xlab), labels = xlab, las = 2, cex.axis = 0.8)
      abline(h = 0, lty = 2)
      grDevices::dev.off()
      plot_files$tangent <- f2
      
      # 3) PCA arrow plot (PC1 vs PC2)
      if (!is.null(pca_res) && !is.null(pca_res$scores) && ncol(pca_res$scores) >= 2) {
        f3 <- file.path(dir_figs, paste0("PCA_Arrows_", out_prefix, "_K", pca_res$k_eff, ".png"))
        grDevices::png(f3, width = 1400, height = 1000, res = 200)
        op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
        par(mar = c(5,5,4,2))
        sc <- pca_res$scores
        plot(sc[,1], sc[,2], pch = 16, xlab = "PC1 (fixed-basis)", ylab = "PC2 (fixed-basis)",
             main = paste0(label, ": Day centroids in PCA space (training ≤ ", pca_res$train_max_day, ")"))
        text(sc[,1], sc[,2], labels = day_names, pos = 3, cex = 0.9)
        for (i in 1:(nrow(sc)-1)) {
          arrows(sc[i,1], sc[i,2], sc[i+1,1], sc[i+1,2], length = 0.08)
        }
        grDevices::dev.off()
        plot_files$pca_arrows <- f3
      }
    }
    
    list(table = tab,
         table_path = out_csv,
         pca_table_path = out_pca_csv,
         pca = pca_res,
         plot_files = plot_files)
  }
  
  bootstrap_interval_CI <- function(days, day_vec, group_vec, group_value, Xmat,
                                    mode = c("mean","disp"),
                                    B = 500, eps_var = 1e-8, pca_train_max_day = 14, pca_k = 20) {
    mode <- match.arg(mode)
    # Precompute original day vectors for dimension sizes
    base <- compute_day_vectors(Xmat, day_vec, group_vec, group_value, eps_var = eps_var)
    d <- base$days
    D <- length(d)
    intervals <- data.frame(day_from = d[-D], day_to = d[-1], stringsAsFactors = FALSE)
    n_int <- nrow(intervals)
    
    # storage: along_s, orth_norm, cosine_to_global
    along_mat <- matrix(NA_real_, nrow = B, ncol = n_int)
    orth_mat  <- matrix(NA_real_, nrow = B, ncol = n_int)
    cos_mat   <- matrix(NA_real_, nrow = B, ncol = n_int)
    tang_mat  <- matrix(NA_real_, nrow = B, ncol = n_int)
    
    # identify rows per day
    idxg <- which(group_vec == group_value)
    rows_by_day <- lapply(d, function(dd) idxg[day_vec[idxg] == dd])
    names(rows_by_day) <- as.character(d)
    
    for (b in seq_len(B)) {
      # resample birds within each day
      sampled_idx <- integer(0)
      for (dd in d) {
        rows <- rows_by_day[[as.character(dd)]]
        if (length(rows) < 1) next
        sampled_idx <- c(sampled_idx, sample(rows, size = length(rows), replace = TRUE))
      }
      Xb <- Xmat[sampled_idx, , drop = FALSE]
      dayb <- day_vec[sampled_idx]
      groupb <- group_vec[sampled_idx]
      
      vb <- compute_day_vectors(Xb, dayb, groupb, rep(group_value, length(groupb)), eps_var = eps_var)
      
      xlist <- if (mode == "mean") vb$mean_by_day else vb$disp_resid_by_day
      # build day matrix
      Xday <- do.call(rbind, lapply(as.character(d), function(nm) xlist[[nm]]))
      storage.mode(Xday) <- "double"
      
      # global axis
      vG <- Xday[D,] - Xday[1,]
      if (!is.finite(vnorm(vG)) || vnorm(vG) < 1e-12) next
      vhat <- vG / vnorm(vG)
      
      # interval metrics
      for (i in seq_len(n_int)) {
        delta <- Xday[i+1,] - Xday[i,]
        cos_mat[b,i] <- cosine(delta, vG)
        s <- dot(delta, vhat)
        perp <- delta - s*vhat
        along_mat[b,i] <- s
        orth_mat[b,i] <- vnorm(perp)
      }
      
      # tangent continuity
      tangents <- matrix(NA_real_, nrow = D, ncol = ncol(Xday))
      for (i in 2:(D-1)) tangents[i,] <- Xday[i+1,] - Xday[i-1,]
      for (i in 2:(D-2)) tang_mat[b,i] <- cosine(tangents[i,], tangents[i+1,])
    }
    
    # CI helper
    ci <- function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 20) return(c(NA_real_, NA_real_))
      stats::quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    }
    
    # build CI table wide
    ci_cos <- t(apply(cos_mat, 2, ci))
    ci_along <- t(apply(along_mat, 2, ci))
    ci_orth <- t(apply(orth_mat, 2, ci))
    ci_tang <- t(apply(tang_mat, 2, ci))
    
    out <- intervals
    out$cos_ci_low <- ci_cos[,1]; out$cos_ci_high <- ci_cos[,2]
    out$along_ci_low <- ci_along[,1]; out$along_ci_high <- ci_along[,2]
    out$orth_ci_low <- ci_orth[,1]; out$orth_ci_high <- ci_orth[,2]
    out$tangent_ci_low <- ci_tang[,1]; out$tangent_ci_high <- ci_tang[,2]
    out
  }
  
  # piecewise regime vectors (user boundaries)
  piecewise_regimes <- function(days, x_by_day, boundaries) {
    # boundaries is numeric vector including endpoints, e.g., c(4,8,14,20)
    # we map to closest available days in 'days'
    closest_day <- function(target) {
      days[which.min(abs(days - target))]
    }
    b2 <- unique(vapply(boundaries, closest_day, numeric(1)))
    b2 <- sort(b2)
    # ensure endpoints included
    if (b2[1] != days[1]) b2 <- c(days[1], b2)
    if (b2[length(b2)] != days[length(days)]) b2 <- c(b2, days[length(days)])
    b2 <- unique(b2)
    b2 <- sort(b2)
    
    segs <- list()
    for (i in 1:(length(b2)-1)) {
      a <- b2[i]; b <- b2[i+1]
      seg_days <- days[days >= a & days <= b]
      if (length(seg_days) < 2) next
      v <- x_by_day[[as.character(seg_days[length(seg_days)])]] - x_by_day[[as.character(seg_days[1])]]
      segs[[paste0(a, "_to_", b)]] <- list(start = seg_days[1], end = seg_days[length(seg_days)],
                                           days = seg_days, v = v)
    }
    segs
  }
  
  segment_coherence_table <- function(days, x_by_day, segs) {
    # within-segment alignment: cosine(step, segment_vector) for steps in segment
    if (length(segs) == 0) return(NULL)
    rows <- list()
    for (nm in names(segs)) {
      seg <- segs[[nm]]
      v <- seg$v
      if (!is.finite(vnorm(v)) || vnorm(v) < 1e-12) next
      seg_days <- seg$days
      if (length(seg_days) < 2) next
      cs <- c()
      for (i in 1:(length(seg_days)-1)) {
        d1 <- seg_days[i]; d2 <- seg_days[i+1]
        delta <- x_by_day[[as.character(d2)]] - x_by_day[[as.character(d1)]]
        cs <- c(cs, cosine(delta, v))
      }
      rows[[nm]] <- data.frame(
        segment = nm,
        start = seg$start,
        end = seg$end,
        n_steps = length(seg_days)-1,
        mean_within_alignment = if (length(cs)>0) mean(cs, na.rm = TRUE) else NA_real_,
        median_within_alignment = if (length(cs)>0) stats::median(cs, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
  
  # ------------------------------ run analysis -------------------------------
  sid <- get_script_identity()
  
  # manifest scaffolding (will fill after outputs exist)
  manifest <- list(
    script_name = sid$script_name,
    script_path = sid$script_path,
    script_full = sid$script_full,
    input_file = input_file_abs,
    output_dir = out_dir,
    timestamp = stamp,
    transform = transform,
    pseudocount = pseudocount,
    pca_k = pca_k,
    pca_train_max_day = pca_train_max_day,
    bootstrap_B = bootstrap_B,
    permute_P = permute_P,
    sample_col = sample_col,
    day_col = if (is.na(day_col)) "PARSED_FROM_SAMPLEID" else day_col,
    group_col = if (is.na(group_col)) "NONE" else group_col,
    run_mode = run_mode,
    n_rows = nrow(df),
    n_gene_cols = length(gene_cols)
  )
  
  results_by_group <- list()
  
  groups_to_run <- if (run_mode == "per_group") groups else "ALL"
  if (run_mode == "per_group") say("Running per group: ", paste(groups_to_run, collapse = ", "))
  
  # regime boundaries prompt (defaults)
  default_bounds <- c(4, 8, 14, 20)
  b_in <- readline(paste0("Piecewise regime boundaries (comma-separated; default ",
                          paste(default_bounds, collapse = ","), "): "))
  if (!nzchar(b_in)) {
    bounds <- default_bounds
  } else {
    bounds <- suppressWarnings(as.numeric(strsplit(gsub(" ", "", b_in), ",", fixed = TRUE)[[1]]))
    bounds <- bounds[is.finite(bounds)]
    if (length(bounds) < 2) bounds <- default_bounds
  }
  manifest$piecewise_boundaries <- bounds
  
  for (gval in groups_to_run) {
    say("----- Group: ", gval, " -----")
    vecs <- compute_day_vectors(X, df$.__DAY__, df$.__GROUP__, gval, eps_var = eps_var)
    
    # replicate count diagnostics
    nrep <- vecs$nrep
    bad_days <- names(nrep)[nrep < 2]
    if (length(bad_days) > 0) {
      say("⚠️ Days with <2 birds in group '", gval, "': ", paste(bad_days, collapse = ", "),
          " — dispersion/residual variance will be NA for those days.")
    }
    
    days <- vecs$days
    
    # ---------------- Bootstrap CIs (bird-within-day) ----------------
    say("🧪 Bootstrap (bird-within-day) B=", bootstrap_B, " [MEAN]")
    ci_mean <- bootstrap_interval_CI(days, df$.__DAY__, df$.__GROUP__, gval, X,
                                     mode = "mean", B = bootstrap_B, eps_var = eps_var)
    say("🧪 Bootstrap (bird-within-day) B=", bootstrap_B, " [DISP]")
    ci_disp <- bootstrap_interval_CI(days, df$.__DAY__, df$.__GROUP__, gval, X,
                                     mode = "disp", B = bootstrap_B, eps_var = eps_var)
    
    # ---------------- Layer evaluation: MEAN ----------------
    prefix_mean <- if (gval == "ALL") "MEAN" else paste0("MEAN_", gval)
    res_mean <- eval_layers(
      days = days,
      x_by_day = vecs$mean_by_day,
      label = paste0("Mean-expression space (", gval, ")"),
      pca_k = pca_k,
      pca_train_max_day = pca_train_max_day,
      bootstrap_CI = ci_mean,
      permute_P = permute_P,
      out_prefix = prefix_mean,
      dir_tables = dir_tables,
      dir_figs = dir_figs
    )
    
    # ---------------- Layer evaluation: DISP (residual dispersion) ----------------
    prefix_disp <- if (gval == "ALL") "DISP" else paste0("DISP_", gval)
    res_disp <- eval_layers(
      days = days,
      x_by_day = vecs$disp_resid_by_day,
      label = paste0("Residual-dispersion space (", gval, ")"),
      pca_k = pca_k,
      pca_train_max_day = pca_train_max_day,
      bootstrap_CI = ci_disp,
      permute_P = permute_P,
      out_prefix = prefix_disp,
      dir_tables = dir_tables,
      dir_figs = dir_figs
    )
    
    # ---------------- Piecewise regimes (interpretation layer) ----------------
    segs_mean <- piecewise_regimes(days, vecs$mean_by_day, bounds)
    segs_disp <- piecewise_regimes(days, vecs$disp_resid_by_day, bounds)
    
    segtab_mean <- segment_coherence_table(days, vecs$mean_by_day, segs_mean)
    segtab_disp <- segment_coherence_table(days, vecs$disp_resid_by_day, segs_disp)
    
    seg_mean_path <- NA_character_
    seg_disp_path <- NA_character_
    if (!is.null(segtab_mean)) {
      seg_mean_path <- file.path(dir_tables, paste0("SegmentCoherence_", prefix_mean, ".csv"))
      utils::write.csv(segtab_mean, seg_mean_path, row.names = FALSE)
    }
    if (!is.null(segtab_disp)) {
      seg_disp_path <- file.path(dir_tables, paste0("SegmentCoherence_", prefix_disp, ".csv"))
      utils::write.csv(segtab_disp, seg_disp_path, row.names = FALSE)
    }
    
    results_by_group[[gval]] <- list(
      days = days,
      nrep = nrep,
      mean = res_mean,
      disp = res_disp,
      segment_mean_path = seg_mean_path,
      segment_disp_path = seg_disp_path
    )
  }
  
  # --------------------------- inventory + manifest ---------------------------
  inv_file <- file.path(out_dir, "INVENTORY_files.txt")
  inv <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  writeLines(inv, con = inv_file)
  
  manifest$inventory_file <- inv_file
  manifest$log_file <- log_file
  
  manifest_path <- file.path(out_dir, "MANIFEST.json")
  write_manifest(manifest_path, manifest)
  
  # ----------------------------- Quarto report --------------------------------
  qmd_path <- file.path(dir_report, paste0("Trajectory_VectorGeometry_Report_", stamp, ".qmd"))
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  
  if (isTRUE(make_quarto)) {
    have_quarto <- FALSE
    if (requireNamespace("quarto", quietly = TRUE)) {
      # check if quarto executable exists
      qbin <- Sys.which("quarto")
      # quarto pkg can still render using bundled quarto in some installs, but PATH issues are common.
      # We'll accept either PATH or quarto:::quarto_path() if available.
      have_quarto <- nzchar(qbin) || TRUE
    }
    if (!requireNamespace("quarto", quietly = TRUE)) {
      say("⚠️ quarto R package not installed; skipping HTML render.")
      have_quarto <- FALSE
    }
    
    # Build QMD as a single character vector (canonical)
    fig_pngs <- list.files(dir_figs, pattern = "\\.png$", full.names = FALSE)
    fig_rel <- file.path("figures", fig_pngs)
    
    out_files_rel <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
    
    qmd <- c(
      "---",
      "title: \"Trajectory Vector Geometry (Wide CSV)\"",
      "format:",
      "  html:",
      "    toc: true",
      "    toc-depth: 3",
      "    number-sections: true",
      "execute:",
      "  echo: false",
      "  warning: false",
      "  message: false",
      "---",
      "",
      "# Summary",
      "",
      "This report quantifies developmental trajectory geometry across all adjacent day intervals using:",
      "",
      "- **Mean-expression state vectors** (engine trajectory)",
      "- **Residual-dispersion state vectors** (constraint/sensors trajectory; mean–variance residualized)",
      "",
      "Primary inference uses **gene space**. Replication uses **fixed-basis PCA (K=20)** trained on early days (default ≤14).",
      "",
      "# Script header (verbatim)",
      "",
      "The script is designed for **SOURCE-to-run** use in RStudio (base R; no dplyr/tidyverse).",
      "",
      "# Metadata",
      "",
      "## Paths",
      paste0("- **Input CSV:** ", input_file_abs),
      paste0("- **Output dir:** ", out_dir),
      paste0("- **Manifest:** ", manifest_path),
      paste0("- **Inventory:** ", inv_file),
      paste0("- **Log:** ", log_file),
      "",
      "## Parameters",
      paste0("- transform: ", transform),
      paste0("- pseudocount: ", pseudocount),
      paste0("- bootstrap_B (birds within day): ", bootstrap_B),
      paste0("- permute_P: ", permute_P),
      paste0("- PCA K: ", pca_k),
      paste0("- PCA training max day: ", pca_train_max_day),
      paste0("- Piecewise boundaries: ", paste(bounds, collapse = ", ")),
      "",
      "# Dependencies",
      "",
      "- Base R: stats, utils, grDevices",
      "- Optional: rstudioapi (file pickers), jsonlite (manifest JSON), quarto (HTML render)",
      "",
      "# Analytical logic",
      "",
      "## Definitions",
      "",
      "Let $x_t$ be the per-day state vector (either mean-expression or residual-dispersion).",
      "",
      "- Global axis: $v_{global} = x_{last} - x_{first}$",
      "- Adjacent step: $\\Delta_t = x_{t+1} - x_t$",
      "- Cosine alignment: $\\cos(\\theta) = \\frac{\\Delta_t \\cdot v_{global}}{\\|\\Delta_t\\|\\|v_{global}\\|}$",
      "- Along-axis progress: $s_t = \\Delta_t \\cdot \\hat v_{global}$",
      "- Orthogonal excursion: $\\|\\Delta_{t,\\perp}\\| = \\|\\Delta_t - s_t\\hat v_{global}\\|$",
      "",
      "## Tangent continuity (local direction)",
      "",
      "Tangent at day index $i$ uses symmetric difference $v_i = x_{i+1} - x_{i-1}$.",
      "Continuity across adjacent boundaries is computed as $\\cos(v_i, v_{i+1})$ where defined.",
      "",
      "## Residual dispersion (mean–variance adjustment)",
      "",
      "For each day, per-gene within-day variance is computed across birds, then:",
      "",
      "- $\\log(Var_g)$ is regressed on $\\log(|Mean_g|)$ across genes within the day",
      "- Residuals define the residual-dispersion state vector",
      "",
      "# Key outputs",
      "",
      "Tables are in `tables/` and figures in `figures/`.",
      "",
      "## Figures",
      if (length(fig_rel) > 0) "" else "(No figures found.)",
      unlist(lapply(fig_rel, function(f) c(paste0("![](", f, ")"), ""))),
      "",
      "## Generated outputs (inventory)",
      "",
      "```",
      paste(out_files_rel, collapse = "\n"),
      "```",
      "",
      "# Session info",
      "",
      "```",
      capture.output(sessionInfo()),
      "```"
    )
    
    writeLines(qmd, con = qmd_path)
    
    if (have_quarto) {
      say("🧾 Rendering Quarto HTML report (canonical WD switch)...")
      old_wd <- getwd()
      setwd(out_dir)
      on.exit(setwd(old_wd), add = TRUE)
      
      # IMPORTANT: render using relative path from output_dir
      qmd_rel <- file.path("report", basename(qmd_path))
      try(quarto::quarto_render(qmd_rel), silent = TRUE)
    }
  } else {
    say("Skipping Quarto report generation (make_quarto=FALSE).")
  }
  
  # final object
  res <- list(
    out_dir = out_dir,
    tables_dir = dir_tables,
    figures_dir = dir_figs,
    report_dir = dir_report,
    log_file = log_file,
    manifest_path = manifest_path,
    inventory_path = inv_file,
    qmd_path = qmd_path,
    html_path = html_path,
    results_by_group = results_by_group
  )
  
  say("===== Trajectory Vector Geometry DONE =====")
  say("Output dir: ", out_dir)
  say("Manifest: ", manifest_path)
  say("Inventory: ", inv_file)
  say("Log: ", log_file)
  if (file.exists(html_path)) say("HTML report: ", html_path) else say("HTML report not found (Quarto may be missing).")
  
  return(res)
}

# ---------------------------- AUTO-RUN WHEN SOURCED --------------------------
if (interactive()) {
  message("Running trajectory_vector_geometry_wide() interactively...")
  res <- trajectory_vector_geometry_wide()
  message("Done. Results written to:")
  message(res$out_dir)
  if (!is.null(res$results_by_group)) {
    # try to open a main table (MEAN pooled if present)
    tab_try <- file.path(res$tables_dir, "AdjacentIntervals_VectorMetrics_MEAN.csv")
    if (file.exists(tab_try)) {
      tt <- try(read.csv(tab_try, stringsAsFactors = FALSE), silent = TRUE)
      if (!inherits(tt, "try-error")) {
        try(View(tt), silent = TRUE)
      }
    }
  }
}