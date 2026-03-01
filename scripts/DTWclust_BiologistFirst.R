
###############################################################################
# DTWclust_BiologistFirst.R  (SOURCE-to-run in RStudio; biologist-first)
#
# Purpose:
# - Read WIDE expression CSV with a required Day column.
# - Collapse replicates by Day (mean of numeric feature columns).
# - Optional filter to top N most variable features.
# - DTW clustering (dtwclust) with:
#     (1) k sweep for silhouette + stability (Jaccard) across k=2..10
#     (2) optional override best_k
#     (3) sensitivity sweep across type × distance
#     (4) repeated runs for consensus seed profiles
# - Writes ALL outputs into: ./outputs/<run_name>_<timestamp>/
# - Captures log, manifest, inventory, and optional Quarto report.
#
# Dependencies (CRAN):
#   dtwclust, dtw, cluster, ggplot2
# Optional:
#   rstudioapi (for Mac-friendly pickers), jsonlite (manifest), quarto (report),
#   future.apply (parallel repeats; optional; script will run serial if missing)
#
# NO dplyr / NO tidyverse.
###############################################################################

# ----------------------------- user-facing knobs -----------------------------
DEFAULT_TOP_N <- 200
DEFAULT_KS <- 2:10
DEFAULT_N_REPEATS <- 10
DEFAULT_METHODS <- c("partitional", "hierarchical")
DEFAULT_DISTANCES <- c("dtw_basic", "gak")   # dtwclust distances
DEFAULT_REPEAT_METHOD <- "partitional"
DEFAULT_REPEAT_DISTANCE <- "gak"
DEFAULT_SEED <- 123

# ------------------------------- helpers -------------------------------------
say <- function(...) cat(..., "\n", sep = "")
now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

is_rstudio <- function() {
  if (!requireNamespace("rstudioapi", quietly = TRUE)) return(FALSE)
  isTRUE(rstudioapi::isAvailable())
}

pick_file <- function(caption = "Select input CSV") {
  if (is_rstudio()) {
    f <- try(rstudioapi::selectFile(caption = caption), silent = TRUE)
    if (!inherits(f, "try-error") && nzchar(f)) return(f)
  }
  file.choose()
}

pick_dir <- function(caption = "Select output parent directory (optional)") {
  if (is_rstudio()) {
    d <- try(rstudioapi::selectDirectory(caption = caption), silent = TRUE)
    if (!inherits(d, "try-error") && nzchar(d)) return(d)
  }
  return("")
}

safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

script_identity <- function() {
  # best-effort: where is THIS script?
  script_path <- NA_character_
  script_name <- NA_character_
  script_full <- NA_character_
  
  # Common RStudio / source() patterns
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  hit <- grep(fileArg, cmdArgs, value = TRUE)
  if (length(hit) == 1) {
    script_full <- sub(fileArg, "", hit)
  } else {
    # RStudio API fallback
    if (is_rstudio()) {
      ctx <- try(rstudioapi::getActiveDocumentContext(), silent = TRUE)
      if (!inherits(ctx, "try-error")) script_full <- ctx$path
    }
  }
  
  if (!is.na(script_full) && nzchar(script_full) && file.exists(script_full)) {
    script_full <- normalizePath(script_full, winslash = "/", mustWork = TRUE)
    script_path <- dirname(script_full)
    script_name <- basename(script_full)
  }
  
  list(script_name = script_name, script_path = script_path, script_full = script_full)
}

write_manifest <- function(path, manifest_list) {
  txt <- NULL
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    txt <- jsonlite::toJSON(manifest_list, pretty = TRUE, auto_unbox = TRUE, null = "null")
    writeLines(txt, con = path)
  } else {
    # fallback: simple dput
    dput(manifest_list, file = path)
  }
}

write_inventory <- function(out_dir, file = "INVENTORY_files.txt") {
  inv <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  inv_rel <- sub(paste0("^", gsub("([\\^\\$\\.|\\+\\*\\?\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", out_dir), "/?"),
                 "", inv)
  writeLines(inv_rel, con = file.path(out_dir, file))
  invisible(inv_rel)
}

# Collapse replicates by Day (mean)
collapse_by_day <- function(df) {
  if (!("Day" %in% names(df))) stop("❌ Required column 'Day' not found.")
  day <- df$Day
  num_cols <- sapply(df, is.numeric)
  # exclude Day itself from features
  num_cols["Day"] <- FALSE
  feat_names <- names(df)[num_cols]
  if (length(feat_names) < 2) stop("❌ Not enough numeric feature columns found.")
  
  days <- sort(unique(day))
  out <- data.frame(Day = days, stringsAsFactors = FALSE)
  
  for (fn in feat_names) {
    out[[fn]] <- vapply(days, function(d) {
      vals <- df[day == d, fn]
      mean(vals, na.rm = TRUE)
    }, numeric(1))
  }
  out
}

top_var_filter <- function(expr_mat, top_n = 200) {
  v <- apply(expr_mat, 1, var, na.rm = TRUE)
  ord <- order(v, decreasing = TRUE)
  keep <- ord[seq_len(min(top_n, length(ord)))]
  expr_mat[keep, , drop = FALSE]
}

# Silhouette from precomputed distance object/matrix
silhouette_from_dist <- function(dist_obj, labels) {
  # labels must be integer/factor of length n
  if (!requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
  labs <- as.integer(as.factor(labels))
  if (length(unique(labs)) < 2) return(NA_real_)
  if (any(table(labs) < 2)) return(NA_real_)
  sil <- try(cluster::silhouette(labs, dist_obj), silent = TRUE)
  if (inherits(sil, "try-error")) return(NA_real_)
  mean(sil[, "sil_width"])
}

# Jaccard between two sets of feature names
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

# Match clusters between two labelings by maximizing Jaccard (Hungarian-lite greedy)
# returns mean Jaccard across matched clusters
mean_jaccard_matched <- function(labels1, labels2, feature_names) {
  labs1 <- as.integer(as.factor(labels1))
  labs2 <- as.integer(as.factor(labels2))
  k1 <- length(unique(labs1)); k2 <- length(unique(labs2))
  k <- min(k1, k2)
  
  sets1 <- lapply(seq_len(k1), function(c) feature_names[labs1 == c])
  sets2 <- lapply(seq_len(k2), function(c) feature_names[labs2 == c])
  
  # Build similarity matrix
  sim <- matrix(0, nrow = k1, ncol = k2)
  for (i in seq_len(k1)) for (j in seq_len(k2)) sim[i, j] <- jaccard(sets1[[i]], sets2[[j]])
  
  # Greedy match: repeatedly take max remaining cell
  used_i <- rep(FALSE, k1)
  used_j <- rep(FALSE, k2)
  chosen <- numeric(0)
  for (step in seq_len(k)) {
    sim2 <- sim
    sim2[used_i, ] <- -1
    sim2[, used_j] <- -1
    mx <- max(sim2)
    if (mx < 0) break
    idx <- which(sim2 == mx, arr.ind = TRUE)[1, ]
    used_i[idx[1]] <- TRUE
    used_j[idx[2]] <- TRUE
    chosen <- c(chosen, mx)
  }
  if (length(chosen) == 0) return(NA_real_)
  mean(chosen)
}

# Compute stability for a given k via repeated runs: mean matched Jaccard vs run1
stability_for_k <- function(expr_mat, k, n_repeats, type = "partitional", distance = "dtw_basic",
                            seed0 = 1000) {
  if (!requireNamespace("dtwclust", quietly = TRUE)) stop("dtwclust not installed.")
  feat <- rownames(expr_mat)
  
  labels_list <- vector("list", n_repeats)
  for (r in seq_len(n_repeats)) {
    set.seed(seed0 + r)
    model <- try(dtwclust::tsclust(expr_mat,
                                   type = type,
                                   k = k,
                                   distance = distance,
                                   seed = seed0 + r,
                                   trace = FALSE),
                 silent = TRUE)
    if (inherits(model, "try-error")) return(NA_real_)
    labels_list[[r]] <- model@cluster
  }
  
  base <- labels_list[[1]]
  vals <- numeric(0)
  for (r in 2:n_repeats) {
    vals <- c(vals, mean_jaccard_matched(base, labels_list[[r]], feat))
  }
  mean(vals, na.rm = TRUE)
}

# Choose seed profiles: closest to centroid within each cluster
pick_seed_profiles <- function(expr_mat, labels, centroids) {
  if (!requireNamespace("dtw", quietly = TRUE)) stop("dtw package needed for seed profiles.")
  feat <- rownames(expr_mat)
  k <- length(unique(labels))
  out <- data.frame(Cluster = seq_len(k), SeedFeature = NA_character_, stringsAsFactors = FALSE)
  
  for (c in seq_len(k)) {
    members <- feat[labels == c]
    if (length(members) == 0) next
    centroid <- centroids[[c]]
    # compute DTW distance to centroid for each member
    dists <- vapply(members, function(f) {
      dtw::dtw(as.numeric(expr_mat[f, ]), as.numeric(centroid), distance.only = TRUE)$distance
    }, numeric(1))
    out$SeedFeature[c] <- members[which.min(dists)]
  }
  out
}

# Simple long-format without tidyverse
to_long <- function(df_wide, id_col = "k") {
  # df_wide: columns: id_col + metrics...
  metrics <- setdiff(names(df_wide), id_col)
  out <- data.frame(k = integer(0), Metric = character(0), Value = numeric(0), stringsAsFactors = FALSE)
  for (m in metrics) {
    out <- rbind(out, data.frame(k = df_wide[[id_col]], Metric = m, Value = df_wide[[m]],
                                 stringsAsFactors = FALSE))
  }
  out
}

# ----------------------------- main analysis ---------------------------------
run_dtwclust_biologist_first <- function(make_quarto = TRUE) {
  
  
  # ---- load packages (fail cleanly) ----
  needed <- c("dtwclust", "ggplot2", "cluster", "proxy")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) stop("❌ Missing packages: ", paste(missing, collapse = ", "),
                                "\nInstall them first (install.packages).")
  
  # dtw only needed for seed profiles; warn later if missing
  if (!requireNamespace("dtw", quietly = TRUE)) {
    say("⚠️ Package 'dtw' not installed: seed-profile selection will be skipped unless installed.")
  }
  
  ts <- now_stamp()
  ids <- script_identity()
  
  # ---- interactive inputs ----
  say("\n📂 Select input CSV (WIDE). Must contain column: Day")
  input_file <- pick_file("Select input CSV (must contain Day column)")
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  
  run_name <- readline("📝 Enter run name (e.g., DTW_Liver5SD): ")
  if (!nzchar(run_name)) run_name <- "DTW_run"
  
  say("\n📁 Output location")
  say("   - Default is ./outputs/")
  say("   - Optionally choose a different parent directory (Enter to skip)")
  parent_choice <- ""
  if (interactive()) {
    resp <- readline("Choose output parent directory with picker? (y/N): ")
    if (tolower(trimws(resp)) == "y") parent_choice <- pick_dir("Select output parent directory")
  }
  
  parent_dir <- if (nzchar(parent_choice)) parent_choice else file.path(getwd(), "outputs")
  parent_dir <- safe_mkdir(parent_dir)
  
  out_dir <- safe_mkdir(file.path(parent_dir, paste0(run_name, "_", ts)))
  plots_dir <- safe_mkdir(file.path(out_dir, "plots"))
  tables_dir <- safe_mkdir(file.path(out_dir, "tables"))
  
  # ---- capture log ----
  log_file <- file.path(out_dir, paste0("LOG_", ts, ".txt"))
  sink(log_file, split = TRUE)
  on.exit({
    sink()
  }, add = TRUE)
  
  say("✅ Output directory: ", out_dir)
  say("🧾 Log file: ", log_file)
  
  # ---- read data ----
  say("\n🔎 Reading: ", input_file)
  data_raw <- read.csv(input_file, check.names = FALSE)
  if (!("Day" %in% names(data_raw))) stop("❌ 'Day' column not found in input file.")
  # ensure Day numeric-ish but keep ordering
  # do not coerce if it breaks; just carry
  say("   Rows: ", nrow(data_raw), "  Cols: ", ncol(data_raw))
  
  # ---- collapse by Day ----
  say("\n🧮 Collapsing replicates by Day (mean of numeric feature columns)...")
  collapsed <- collapse_by_day(data_raw)
  write.csv(collapsed, file.path(tables_dir, "CollapsedByDay_Mean.csv"), row.names = FALSE)
  
  timepoints <- collapsed$Day
  # expr_matrix: features x timepoints
  expr_matrix <- t(as.matrix(collapsed[, setdiff(names(collapsed), "Day"), drop = FALSE]))
  rownames(expr_matrix) <- setdiff(names(collapsed), "Day")
  colnames(expr_matrix) <- as.character(timepoints)
  
  say("✅ Collapsed matrix dimensions (features x timepoints): ",
      nrow(expr_matrix), " x ", ncol(expr_matrix))
  
  # ---- optional filter ----
  top_n_in <- readline(paste0("🔎 Top N most-variable features to keep [default ", DEFAULT_TOP_N, "]: "))
  top_n <- suppressWarnings(as.integer(top_n_in))
  if (is.na(top_n) || top_n < 10) top_n <- DEFAULT_TOP_N
  expr_filt <- top_var_filter(expr_matrix, top_n = top_n)
  say("✅ After filter: ", nrow(expr_filt), " features retained.")
  write.csv(data.frame(Feature = rownames(expr_filt)), file.path(tables_dir, "Features_Kept.csv"),
            row.names = FALSE)
  
  # ---- precompute DTW distance matrix once (dtw_basic) for silhouette ----
  say("\n🔄 Precomputing DTW distances for silhouette (dtw_basic) via proxy::dist()...")
  
  # dtwclust registers "dtw_basic" with proxy in most installations, but be defensive:
  dist_obj <- try(proxy::dist(expr_filt, method = "dtw_basic"), silent = TRUE)
  
  if (inherits(dist_obj, "try-error")) {
    say("⚠️ proxy::dist(method='dtw_basic') failed. Falling back to proxy::dist(method='DTW') ...")
    dist_obj <- try(proxy::dist(expr_filt, method = "DTW"), silent = TRUE)
  }
  
  if (inherits(dist_obj, "try-error")) {
    stop("❌ Could not compute DTW distances via proxy. ",
         "This usually means dtwclust did not register DTW distances with proxy in your install.\n",
         "Try: install.packages(c('dtwclust','proxy','dtw')) and restart R.")
  }
  
  saveRDS(dist_obj, file.path(out_dir, "DTW_basic_dist.rds"))
  
  # Optional CSV export (can be large); keep but warn
  dist_mat <- as.matrix(dist_obj)
  write.csv(dist_mat, file.path(tables_dir, "DTW_Distance_Matrix_dtw_basic.csv"))
  
  # ---- k sweep: silhouette + stability ----
  say("\n📊 Evaluating silhouette + stability across k...")
  ks <- DEFAULT_KS
  sil <- rep(NA_real_, length(ks))
  jac <- rep(NA_real_, length(ks))
  
  for (i in seq_along(ks)) {
    k <- ks[i]
    say("   🔹 k = ", k)
    
    # cluster once for silhouette (partitional dtw_basic)
    model <- try(dtwclust::tsclust(expr_filt, type = "partitional", k = k,
                                   distance = "dtw_basic", seed = DEFAULT_SEED, trace = FALSE),
                 silent = TRUE)
    if (!inherits(model, "try-error")) {
      labels <- model@cluster
      sil[i] <- silhouette_from_dist(dist_obj, labels)
    }
    
    # stability via repeats (mean matched Jaccard)
    jac[i] <- stability_for_k(expr_filt, k = k, n_repeats = DEFAULT_N_REPEATS,
                              type = "partitional", distance = "dtw_basic", seed0 = 100)
    say("      silhouette=", round(sil[i], 4), "  stability(Jaccard)=", round(jac[i], 4))
  }
  
  combo <- data.frame(k = ks, Silhouette = sil, JaccardStability = jac)
  write.csv(combo, file.path(tables_dir, "DTWclust_Silhouette_Jaccard_Combo.csv"), row.names = FALSE)
  
  # plot combo
  combo_long <- to_long(combo, id_col = "k")
  p <- ggplot2::ggplot(combo_long, ggplot2::aes(x = k, y = Value, color = Metric)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(title = "Silhouette Score & Jaccard Stability vs Number of Clusters",
                  x = "Number of Clusters (k)", y = "Score")
  ggplot2::ggsave(filename = file.path(plots_dir, "DTWclust_Silhouette_Jaccard_ComboPlot.png"),
                  plot = p, width = 8, height = 10, dpi = 600)
  
  # choose auto k by silhouette (you can change this rule)
  auto_k <- combo$k[which.max(combo$Silhouette)]
  say("\n🏆 Auto-selected best k = ", auto_k,
      " (max silhouette = ", round(max(combo$Silhouette, na.rm = TRUE), 4), ")")
  
  user_input <- readline(paste0("🔒 Enter k to override (Enter to accept ", auto_k, "): "))
  if (nzchar(user_input)) {
    best_k <- suppressWarnings(as.integer(user_input))
    if (is.na(best_k) || best_k < 2) {
      say("⚠️ Invalid input. Using auto_k.")
      best_k <- auto_k
    } else {
      say("✅ Using override best_k = ", best_k)
    }
  } else {
    best_k <- auto_k
    say("✅ Proceeding with auto-selected best_k = ", best_k)
  }
  
  # ---- sensitivity sweep: type × distance (uses best_k) ----
  say("\n🧪 Sensitivity sweep (type × distance) using best_k = ", best_k, " ...")
  methods <- DEFAULT_METHODS
  distances <- DEFAULT_DISTANCES
  
  sens <- data.frame(Method = character(0), Distance = character(0), Config = character(0),
                     Silhouette = numeric(0), ClusterFile = character(0),
                     stringsAsFactors = FALSE)
  
  for (type in methods) {
    for (dist_name in distances) {
      config_id <- paste(type, dist_name, sep = "_")
      say("   ⚙️ Trying ", config_id)
      
      sil_here <- NA_real_
      cluster_file <- "NA"
      
      model <- try(dtwclust::tsclust(expr_filt, type = type, k = best_k,
                                     distance = dist_name, seed = DEFAULT_SEED, trace = FALSE),
                   silent = TRUE)
      if (!inherits(model, "try-error")) {
        labels <- model@cluster
        cluster_file <- file.path(tables_dir, paste0("Cluster_Labels_", config_id, ".csv"))
        write.csv(data.frame(Feature = rownames(expr_filt), Cluster = labels),
                  cluster_file, row.names = FALSE)
        
        # silhouette only valid when we have a distance consistent with dist_obj
        # For comparability, we report silhouette computed on dtw_basic dist_obj regardless of dist_name.
        sil_here <- silhouette_from_dist(dist_obj, labels)
      } else {
        say("      ❌ Failed: ", as.character(model))
      }
      
      sens <- rbind(sens, data.frame(Method = type, Distance = dist_name, Config = config_id,
                                     Silhouette = sil_here,
                                     ClusterFile = basename(cluster_file),
                                     stringsAsFactors = FALSE))
    }
  }
  
  write.csv(sens, file.path(tables_dir, "DTWclust_Sensitivity_Full_Table.csv"), row.names = FALSE)
  
  # pick best by silhouette (computed on dtw_basic dist)
  ok <- which(!is.na(sens$Silhouette))
  if (length(ok) == 0) stop("❌ No valid clustering results in sensitivity sweep.")
  best_idx <- ok[which.max(sens$Silhouette[ok])]
  best_cfg <- sens[best_idx, , drop = FALSE]
  say("\n✅ Best sweep config: ", best_cfg$Config, "  silhouette=", round(best_cfg$Silhouette, 4))
  
  best_assign_path <- file.path(tables_dir, best_cfg$ClusterFile)
  best_assign <- read.csv(best_assign_path, stringsAsFactors = FALSE)
  write.csv(best_assign, file.path(tables_dir, "DTWclust_BestSweep_Cluster_Assignments.csv"),
            row.names = FALSE)
  
  # ---- repeats for consensus seeds (optional) ----
  say("\n🔁 Repeating clustering for consensus seeds...")
  repeat_method <- DEFAULT_REPEAT_METHOD
  repeat_distance <- DEFAULT_REPEAT_DISTANCE
  n_repeats <- DEFAULT_N_REPEATS
  
  do_parallel <- requireNamespace("future.apply", quietly = TRUE)
  if (do_parallel) {
    # best effort parallel; if it fails, fall back
    try({
      future::plan(future::multisession)
    }, silent = TRUE)
  }
  
  repeat_fun <- function(i) {
    set.seed(100 + i)
    model <- dtwclust::tsclust(expr_filt, type = repeat_method, k = best_k,
                               distance = repeat_distance, seed = 100 + i, trace = FALSE)
    labels <- model@cluster
    
    if (!requireNamespace("dtw", quietly = TRUE)) {
      return(list(labels = labels, seeds = NULL))
    }
    
    # centroids: dtwclust stores as list-like structure for partitional
    cents <- model@centroids
    # normalize centroid container
    # Some versions store a matrix; some store list; we handle both
    if (is.matrix(cents)) {
      centroids <- lapply(seq_len(nrow(cents)), function(r) cents[r, ])
    } else if (is.list(cents)) {
      centroids <- cents
    } else {
      centroids <- NULL
    }
    
    seeds_df <- NULL
    if (!is.null(centroids)) {
      seeds_df <- pick_seed_profiles(expr_filt, labels, centroids)
    }
    list(labels = labels, seeds = seeds_df)
  }
  
  if (do_parallel) {
    reps <- try(future.apply::future_lapply(seq_len(n_repeats), repeat_fun, future.seed = TRUE),
                silent = TRUE)
    if (inherits(reps, "try-error")) {
      say("⚠️ Parallel repeats failed; running serially.")
      reps <- lapply(seq_len(n_repeats), repeat_fun)
    }
  } else {
    reps <- lapply(seq_len(n_repeats), repeat_fun)
  }
  
  # save membership matrix
  membership <- do.call(cbind, lapply(reps, function(x) x$labels))
  rownames(membership) <- rownames(expr_filt)
  colnames(membership) <- paste0("Run", seq_len(n_repeats))
  write.csv(membership, file.path(tables_dir, "ClusterMembership_AllRuns.csv"))
  
  # consensus seeds if available
  if (requireNamespace("dtw", quietly = TRUE)) {
    all_seeds <- do.call(rbind, lapply(seq_along(reps), function(i) {
      sd <- reps[[i]]$seeds
      if (is.null(sd)) return(NULL)
      sd$Run <- i
      sd
    }))
    
    if (!is.null(all_seeds) && nrow(all_seeds) > 0) {
      # frequency table: (Cluster, SeedFeature)
      key <- paste(all_seeds$Cluster, all_seeds$SeedFeature, sep = "||")
      freq <- table(key)
      parts <- strsplit(names(freq), "\\|\\|")
      consensus <- data.frame(
        Cluster = as.integer(vapply(parts, `[`, "", 1)),
        SeedFeature = vapply(parts, `[`, "", 2),
        Frequency = as.integer(freq),
        stringsAsFactors = FALSE
      )
      
      # pick most frequent per cluster
      consensus <- consensus[order(consensus$Cluster, -consensus$Frequency), ]
      keep <- !duplicated(consensus$Cluster)
      consensus_best <- consensus[keep, ]
      write.csv(consensus_best, file.path(tables_dir, "DTWclust_Consensus_SeedProfiles.csv"),
                row.names = FALSE)
      
      # plot seed frequencies
      pp <- ggplot2::ggplot(consensus, ggplot2::aes(x = reorder(SeedFeature, -Frequency),
                                                    y = Frequency)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~Cluster, scales = "free_y") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(title = "Seed Feature Frequency Across Repeats",
                      x = "Seed Feature", y = "Frequency")
      ggplot2::ggsave(file.path(plots_dir, "SeedFrequency_PerCluster.png"),
                      plot = pp, width = 10, height = 6, dpi = 300)
    } else {
      say("⚠️ Seed profiles were not generated (centroids unavailable or dtw missing).")
    }
  }
  
  # ---- stability summary per cluster (Jaccard across repeats) ----
  say("\n📊 Cluster-wise stability summary (mean pairwise Jaccard across runs)...")
  # Use run1 clusters as reference ordering; compute pairwise Jaccard per cluster label index
  # Note: label switching handled by matched method; for per-cluster stats we keep it simple:
  # report overall matched stability at best_k (already computed in combo), plus a per-cluster
  # "run1 membership reproducibility" via matched mapping to each run.
  feat <- rownames(expr_filt)
  base_labels <- membership[, 1]
  
  per_cluster_stats <- data.frame(Cluster = seq_len(best_k),
                                  MeanJaccard = NA_real_,
                                  MedianJaccard = NA_real_,
                                  SDJaccard = NA_real_,
                                  stringsAsFactors = FALSE)
  
  # build base cluster sets
  base_sets <- lapply(seq_len(best_k), function(c) feat[base_labels == c])
  
  # for each run, match its clusters to base clusters greedily and record matched jaccards
  jacc_by_cluster <- vector("list", best_k)
  for (c in seq_len(best_k)) jacc_by_cluster[[c]] <- numeric(0)
  
  for (r in 2:n_repeats) {
    labs_r <- membership[, r]
    # build sets for run r
    k_r <- length(unique(labs_r))
    sets_r <- lapply(seq_len(k_r), function(c) feat[labs_r == c])
    
    # similarity matrix base_k x k_r
    sim <- matrix(0, nrow = best_k, ncol = k_r)
    for (i in seq_len(best_k)) for (j in seq_len(k_r)) sim[i, j] <- jaccard(base_sets[[i]], sets_r[[j]])
    
    used_j <- rep(FALSE, k_r)
    for (i in seq_len(best_k)) {
      # pick best available j for base cluster i
      vals <- sim[i, ]
      vals[used_j] <- -1
      jbest <- which.max(vals)
      if (length(jbest) == 1 && vals[jbest] >= 0) {
        used_j[jbest] <- TRUE
        jacc_by_cluster[[i]] <- c(jacc_by_cluster[[i]], vals[jbest])
      }
    }
  }
  
  for (c in seq_len(best_k)) {
    v <- jacc_by_cluster[[c]]
    if (length(v) > 0) {
      per_cluster_stats$MeanJaccard[c] <- mean(v)
      per_cluster_stats$MedianJaccard[c] <- median(v)
      per_cluster_stats$SDJaccard[c] <- sd(v)
    }
  }
  write.csv(per_cluster_stats, file.path(tables_dir, "AllClusters_Jaccard_SummaryStats.csv"),
            row.names = FALSE)
  
  # ---- provenance: manifest + inventory + session info ----
  manifest <- list(
    script_name = ids$script_name,
    script_path = ids$script_path,
    script_full = ids$script_full,
    input_file = input_file,
    output_dir = out_dir,
    timestamp = ts,
    parameters = list(
      top_n = top_n,
      ks = as.integer(ks),
      n_repeats = n_repeats,
      auto_k = auto_k,
      best_k = best_k,
      sensitivity_methods = methods,
      sensitivity_distances = distances,
      repeat_method = repeat_method,
      repeat_distance = repeat_distance
    )
  )
  write_manifest(file.path(out_dir, "MANIFEST.json"), manifest)
  write_inventory(out_dir)
  capture.output(sessionInfo(), file = file.path(out_dir, "SESSIONINFO.txt"))
  
  # ---- optional Quarto report (canonical contract) ----
  qmd_path <- file.path(out_dir, paste0("DTWclust_Report_", ts, ".qmd"))
  html_path <- sub("\\.qmd$", ".html", qmd_path)
  
  if (isTRUE(make_quarto) && requireNamespace("quarto", quietly = TRUE)) {
    say("\n🧾 Building Quarto report...")
    header_lines <- c(
      "###############################################################################",
      "# Script header block (verbatim):",
      "# DTWclust_BiologistFirst.R",
      "###############################################################################"
    )
    
    qmd <- c(
      "---",
      "title: \"DTWclust Report\"",
      "format: html",
      "---",
      "",
      "## Summary",
      "- Collapsed replicates by Day and clustered feature trajectories using DTW.",
      "- Selected k using silhouette and assessed stability using matched Jaccard across repeats.",
      "",
      "## Script header",
      "```",
      header_lines,
      "```",
      "",
      "## Metadata",
      paste0("- **Timestamp:** ", ts),
      paste0("- **Input file:** ", input_file),
      paste0("- **Output dir:** ", out_dir),
      paste0("- **top_n:** ", top_n),
      paste0("- **auto_k:** ", auto_k),
      paste0("- **best_k:** ", best_k),
      "",
      "## Dependencies",
      "Packages: dtwclust, dtw (optional), cluster, ggplot2, jsonlite (optional), quarto (optional), future.apply (optional)",
      "",
      "## Outputs generated",
      "```{r}",
      paste0("list.files('", out_dir, "', recursive = TRUE)"),
      "```",
      "",
      "## Key figures",
      paste0("![](", file.path("plots", "DTWclust_Silhouette_Jaccard_ComboPlot.png"), ")"),
      "",
      "## Session info",
      "```{r}",
      "sessionInfo()",
      "```"
    )
    
    writeLines(qmd, con = qmd_path)
    
    # Canonical contract: setwd(output_dir) before render, restore immediately after
    old_wd <- getwd()
    setwd(out_dir)
    on.exit(setwd(old_wd), add = TRUE)
    
    say("🧾 Rendering Quarto HTML report (working dir = output_dir)...")
    try(quarto::quarto_render(basename(qmd_path)), silent = TRUE)
    
    # restore wd happens via on.exit
  } else {
    say("\n⚠️ Quarto not available or make_quarto=FALSE — skipping report.")
  }
  
  say("\n✅ DONE.")
  say("Output directory: ", out_dir)
  say("Key files:")
  say(" - ", file.path(tables_dir, "DTWclust_Silhouette_Jaccard_Combo.csv"))
  say(" - ", file.path(plots_dir, "DTWclust_Silhouette_Jaccard_ComboPlot.png"))
  say(" - ", file.path(tables_dir, "DTWclust_Sensitivity_Full_Table.csv"))
  say(" - ", file.path(tables_dir, "AllClusters_Jaccard_SummaryStats.csv"))
  if (file.exists(html_path)) say(" - ", html_path)
  
  # return object for interactive use
  list(
    out_dir = out_dir,
    tables_dir = tables_dir,
    plots_dir = plots_dir,
    combo = combo,
    best_k = best_k,
    dist_rds = file.path(out_dir, "DTW_basic_dist.rds"),
    report_html = if (file.exists(html_path)) html_path else NA_character_
  )
}

# ----------------------- AUTO-RUN WHEN SOURCED (RStudio) ----------------------
if (interactive()) {
  say("Running run_dtwclust_biologist_first() interactively...")
  res <- run_dtwclust_biologist_first()
  say("Done. Results written to:")
  say(res$out_dir)
  # capture "output file": show log path + open View of combo table if possible
  log_guess <- list.files(res$out_dir, pattern = "^LOG_.*\\.txt$", full.names = TRUE)
  if (length(log_guess) >= 1) say("Log file: ", log_guess[1])
  if (!is.null(res$combo)) {
    try(View(res$combo), silent = TRUE)
  }
}
