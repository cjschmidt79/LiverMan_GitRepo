#!/usr/bin/env Rscript
# ============================================================
# PlotFaithful_EarlyLate_Biology_WITH_Canalization.R
# ============================================================
# Biologist-first SOURCE-to-run script (base R)
# Optional: rstudioapi (for pickers), jsonlite (manifest), quarto (render)
#
# INPUT (required):
#   Expression CSV: rows = samples/replicates, cols = genes
#   - Sample IDs must be either rownames or the first column.
#   - Sample IDs must contain a day integer somewhere (e.g., "14_3", "D14_3", "Day14_rep3").
#
# INPUT (optional):
#   Mapping CSV: must contain columns Gene and Program
#   - Used for program-level summaries and annotated outputs.
#
# OUTPUTS (outputs/<run>_<timestamp>/):
#   - gene_summary.csv
#   - gene_day_variance_long.csv
#   - canalization_gene_scores.csv
#   - gene_integrated_master.csv
#   - gene_regime_calls.csv
#   - top_genes_early_up.csv
#   - top_genes_late_up.csv
#   - top_genes_transition.csv
#   - top_canalization_genes.csv
#   - early_vs_late_volcanoish.png
#   - top_transition_gene_trajectories.png
#   - top_canalization_variance.png
#   - OPTIONAL: program_summary.csv, program_effect_scatter.png, program_trajectories.png
#   - README_outputs.txt
#   - Provenance.qmd (+ Provenance.html if quarto exists)
#   - Project_Manifest.json + Project_Manifest_Files.csv
#
# Canalization definition used here:
#   - "variance checkpoint": replicate variance is minimized at breakpoint day bp
#     collapse_ratio = var(bp) / median(var(other days))
#   - "stabilization": |slope_pre| - |slope_post|  (positive means damped after bp)
#   - Score uses collapse as primary axis; stabilization is a secondary support term
# ============================================================

# -----------------------------
# utilities
# -----------------------------
safe_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}
ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
quiet_require <- function(pkg) suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))

pick_file <- function(prompt = "Select a file") {
  cat("\n", prompt, "\n", sep = "")
  if (quiet_require("rstudioapi") && rstudioapi::isAvailable()) {
    p <- rstudioapi::selectFile(caption = prompt)
    if (!is.null(p) && nzchar(p)) return(p)
    return("")
  }
  return(file.choose())
}

capture_script_identity <- function(verbose = TRUE) {
  script_path <- NA_character_
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", args, value = TRUE)
  if (length(hit) > 0) script_path <- sub("--file=", "", hit[1])
  if (is.na(script_path) || !nzchar(script_path)) {
    if (!is.null(sys.frames()[[1]]$ofile)) script_path <- sys.frames()[[1]]$ofile
  }
  if ((is.na(script_path) || !nzchar(script_path)) && quiet_require("rstudioapi") && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) script_path <- ctx$path
  }
  script_path <- tryCatch(normalizePath(script_path, winslash = "/", mustWork = FALSE),
                          error = function(e) script_path)
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
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
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
  out_json <- file.path(out_dir, "Project_Manifest.json")
  if (quiet_require("jsonlite")) {
    jsonlite::write_json(manifest, out_json, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    esc <- function(x) gsub('"', '\\"', x)
    write_text(out_json, c(
      "{",
      paste0('  "run_name": "', esc(run_name), '",'),
      paste0('  "created": "', esc(manifest$created), '",'),
      paste0('  "output_dir": "', esc(manifest$output_dir), '",'),
      paste0('  "script_name": "', esc(script_id$script_name), '",'),
      paste0('  "script_path": "', esc(as.character(script_id$script_path)), '"'),
      "}"
    ))
  }
  write_csv(inventory_outputs(out_dir), file.path(out_dir, "Project_Manifest_Files.csv"))
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
  
  if (nzchar(Sys.which("quarto"))) {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(out_dir)
    suppressWarnings(suppressMessages(system(paste0('quarto render "', basename(qmd_path), '"'),
                                             ignore.stdout = TRUE, ignore.stderr = TRUE)))
  } else {
    write_text(html_path, c(
      "<html><body>",
      "<h2>Quarto not found</h2>",
      "<p>Provenance.qmd was created, but Quarto is not available on PATH to render HTML.</p>",
      "</body></html>"
    ))
  }
}

# -----------------------------
# parsing + stats
# -----------------------------
parse_day <- function(sample_id) {
  x <- as.character(sample_id)
  if (any(is.na(x) | !nzchar(x))) {
    bad <- which(is.na(x) | !nzchar(x))
    stop("Some sample IDs are NA/empty. First few bad rows: ", paste(head(bad, 10), collapse = ", "))
  }
  # FIRST integer anywhere in the ID
  day_str <- sub("^.*?([0-9]+).*$", "\\1", x, perl = TRUE)
  day <- suppressWarnings(as.integer(day_str))
  if (any(is.na(day))) {
    bad <- which(is.na(day))
    stop("Could not parse day from some sample IDs (no digits found). Examples:\n  ",
         paste(head(x[bad], 10), collapse = "\n  "))
  }
  day
}

cohen_d <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) return(NA_real_)
  sx <- stats::var(x); sy <- stats::var(y)
  sp <- ((nx - 1) * sx + (ny - 1) * sy) / (nx + ny - 2)
  if (!is.finite(sp) || sp <= 0) return(NA_real_)
  (mean(x) - mean(y)) / sqrt(sp)
}

lin_slope <- function(day, expr) {
  ok <- is.finite(day) & is.finite(expr)
  if (sum(ok) < 3) return(c(slope = NA_real_, r2 = NA_real_))
  fit <- lm(expr[ok] ~ day[ok])
  c(slope = unname(coef(fit)[2]), r2 = summary(fit)$r.squared)
}

piecewise_slopes <- function(day, expr, bp) {
  ok <- is.finite(day) & is.finite(expr)
  pre <- ok & (day <= bp)
  post <- ok & (day > bp)
  out <- c(slope_pre = NA_real_, r2_pre = NA_real_, slope_post = NA_real_, r2_post = NA_real_)
  if (sum(pre) >= 3) {
    f1 <- lm(expr[pre] ~ day[pre])
    out["slope_pre"] <- unname(coef(f1)[2])
    out["r2_pre"] <- summary(f1)$r.squared
  }
  if (sum(post) >= 3) {
    f2 <- lm(expr[post] ~ day[post])
    out["slope_post"] <- unname(coef(f2)[2])
    out["r2_post"] <- summary(f2)$r.squared
  }
  out
}

# -----------------------------
# MAIN
# -----------------------------
cat("\n============================================================\n")
cat("Plot-faithful Early/Late Biology Extractor — FULL REWRITE + CANALIZATION\n")
cat("============================================================\n")

script_id <- capture_script_identity(TRUE)
git_meta <- get_git_metadata()

# ---- pick expression
expr_csv <- pick_file("Select expression CSV (rows=replicates, cols=genes; sample IDs like '14_3' or 'D14_3')")
if (is.null(expr_csv) || !nzchar(expr_csv)) stop("No expression file selected.")

expr <- read.csv(expr_csv, check.names = FALSE, stringsAsFactors = FALSE)

# promote first column to rownames if it looks like IDs
if (ncol(expr) >= 2 && !is.numeric(expr[[1]])) {
  rownames(expr) <- as.character(expr[[1]])
  expr <- expr[, -1, drop = FALSE]
}
if (!any(rownames(expr) != "")) stop("No rownames detected. Provide rownames or first column with sample IDs.")

# drop non-numeric columns robustly
num_ok <- vapply(expr, function(z) {
  zz <- suppressWarnings(as.numeric(as.character(z)))
  all(is.na(zz) | is.finite(zz))
}, logical(1))

if (!all(num_ok)) {
  cat("\nWARNING: Dropping non-numeric columns:\n  ",
      paste(names(expr)[!num_ok], collapse = ", "), "\n", sep = "")
  expr <- expr[, num_ok, drop = FALSE]
}

X <- as.matrix(expr)
mode(X) <- "numeric"
gene_names <- colnames(X)
if (is.null(gene_names) || length(gene_names) < 2) stop("Gene columns not detected after numeric cleaning.")

# parse day
day <- parse_day(rownames(X))
days_u <- sort(unique(day))
cat("\nDay counts (sanity check):\n")
print(sort(table(day)))

# ---- early/late cut + breakpoint
cat("\nDefine EARLY vs LATE by day_cut.\n")
cat("EARLY <= day_cut; LATE > day_cut. Default day_cut = 14\n")
day_cut <- suppressWarnings(as.integer(readline("day_cut [14]: ")))
if (!is.finite(day_cut)) day_cut <- 14L

cat("\nBreakpoint day for piecewise slopes + canalization. Default bp = 14\n")
bp <- suppressWarnings(as.integer(readline("bp [14]: ")))
if (!is.finite(bp)) bp <- 14L

early_idx <- which(day <= day_cut)
late_idx <- which(day > day_cut)

# if the chosen cut yields too few samples, switch to a quantile split (automatic, deterministic)
if (length(early_idx) < 3 || length(late_idx) < 3) {
  cat("\nNOTE: day_cut produced too few samples in EARLY or LATE.\n")
  cat("      Switching to quantile split: bottom 60% days = EARLY, top 40% = LATE.\n")
  thr <- stats::quantile(day, probs = 0.60, na.rm = TRUE)
  early_idx <- which(day <= thr)
  late_idx <- which(day > thr)
}
if (length(early_idx) < 3 || length(late_idx) < 3) stop("Not enough samples in early or late group after fallback split.")

# ---- output folder
run_base <- tools::file_path_sans_ext(basename(expr_csv))
run_name <- safe_slug(paste0("PlotFaithful_", run_base, "_cut", day_cut, "_bp", bp))
out_dir <- file.path(getwd(), "outputs", paste0(run_name, "_", ts_stamp()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(expr_csv, file.path(out_dir, basename(expr_csv)), overwrite = TRUE)

cat("\nOutput folder:\n  ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")

# ---- optional mapping
cat("\nOptional: select Gene->Program mapping CSV with columns Gene and Program.\n")
cat("Cancel to skip.\n")
map_path <- ""
if (quiet_require("rstudioapi") && rstudioapi::isAvailable()) {
  tmp <- tryCatch(rstudioapi::selectFile(caption = "Select mapping CSV (Cancel to skip)"), error = function(e) "")
  if (!is.null(tmp) && nzchar(tmp)) map_path <- tmp
} else {
  ans <- trimws(readline("Mapping CSV path (or press Enter to skip): "))
  if (nzchar(ans) && file.exists(ans)) map_path <- ans
}

# guard: if user accidentally selects expression file as mapping, ignore mapping
if (nzchar(map_path)) {
  if (normalizePath(map_path, winslash = "/", mustWork = FALSE) ==
      normalizePath(expr_csv, winslash = "/", mustWork = FALSE)) {
    cat("\nWARNING: Mapping file equals expression file. Ignoring mapping.\n")
    map_path <- ""
  }
}

map <- NULL
if (nzchar(map_path)) {
  map <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE)
  # normalize Gene column
  if (!("Gene" %in% names(map))) {
    hit <- which(tolower(names(map)) == "gene")
    if (length(hit) == 1) names(map)[hit] <- "Gene"
  }
  # normalize Program column (accept common synonyms)
  if (!("Program" %in% names(map))) {
    hit <- which(tolower(names(map)) %in% c("program","functional_program","module","category"))
    if (length(hit) == 1) names(map)[hit] <- "Program"
  }
  if (!all(c("Gene","Program") %in% names(map))) {
    cat("\nWARNING: mapping file lacks Gene+Program columns; ignoring mapping.\n")
    map <- NULL
    map_path <- ""
  } else {
    map <- map[map$Gene %in% gene_names, , drop = FALSE]
    if (nrow(map) == 0) {
      cat("\nWARNING: mapping genes do not match expression genes; ignoring mapping.\n")
      map <- NULL
      map_path <- ""
    } else {
      file.copy(map_path, file.path(out_dir, basename(map_path)), overwrite = TRUE)
    }
  }
}

# ============================================================
# 1) gene_summary: early/late + evidence + slopes + transition
# ============================================================
pc <- 1e-6
ng <- length(gene_names)

gene_summary <- data.frame(
  Gene = gene_names,
  mean_early = NA_real_,
  mean_late = NA_real_,
  log2FC_late_over_early = NA_real_,
  cohen_d_late_minus_early = NA_real_,
  wilcox_p = NA_real_,
  q_BH = NA_real_,
  slope_all = NA_real_,
  r2_all = NA_real_,
  slope_pre = NA_real_,
  r2_pre = NA_real_,
  slope_post = NA_real_,
  r2_post = NA_real_,
  transition_delta_slope = NA_real_,
  transition_score = NA_real_,
  rank_transition = NA_integer_,
  rank_late_up = NA_integer_,
  rank_early_up = NA_integer_,
  stringsAsFactors = FALSE
)

for (j in seq_len(ng)) {
  g <- X[, j]
  e <- g[early_idx]
  l <- g[late_idx]
  
  gene_summary$mean_early[j] <- mean(e, na.rm = TRUE)
  gene_summary$mean_late[j]  <- mean(l, na.rm = TRUE)
  gene_summary$log2FC_late_over_early[j] <- log2((gene_summary$mean_late[j] + pc) / (gene_summary$mean_early[j] + pc))
  
  gene_summary$cohen_d_late_minus_early[j] <- cohen_d(l, e)
  
  if (sum(is.finite(e)) >= 3 && sum(is.finite(l)) >= 3) {
    gene_summary$wilcox_p[j] <- suppressWarnings(wilcox.test(l, e)$p.value)
  }
  
  sl <- lin_slope(day, g)
  gene_summary$slope_all[j] <- sl["slope"]
  gene_summary$r2_all[j]    <- sl["r2"]
  
  pw <- piecewise_slopes(day, g, bp)
  gene_summary$slope_pre[j]  <- pw["slope_pre"]
  gene_summary$r2_pre[j]     <- pw["r2_pre"]
  gene_summary$slope_post[j] <- pw["slope_post"]
  gene_summary$r2_post[j]    <- pw["r2_post"]
  
  ds <- abs(gene_summary$slope_post[j] - gene_summary$slope_pre[j])
  gene_summary$transition_delta_slope[j] <- ds
  r2w <- mean(c(gene_summary$r2_pre[j], gene_summary$r2_post[j]), na.rm = TRUE)
  if (!is.finite(r2w)) r2w <- 0
  gene_summary$transition_score[j] <- ds * (0.25 + r2w)
}

gene_summary$q_BH <- p.adjust(gene_summary$wilcox_p, method = "BH")
gene_summary$rank_transition <- rank(-gene_summary$transition_score, ties.method = "min")
gene_summary$rank_late_up <- rank(-gene_summary$cohen_d_late_minus_early, ties.method = "min")
gene_summary$rank_early_up <- rank(gene_summary$cohen_d_late_minus_early, ties.method = "min")

gene_summary_path <- file.path(out_dir, "gene_summary.csv")
write_csv(gene_summary[order(gene_summary$rank_transition), ], gene_summary_path)

# convenience: top lists
topN <- 25
top_early_path <- file.path(out_dir, "top_genes_early_up.csv")
top_late_path  <- file.path(out_dir, "top_genes_late_up.csv")
top_trans_path <- file.path(out_dir, "top_genes_transition.csv")

write_csv(gene_summary[order(gene_summary$rank_early_up), ][1:min(topN, ng), ], top_early_path)
write_csv(gene_summary[order(gene_summary$rank_late_up),  ][1:min(topN, ng), ], top_late_path)
write_csv(gene_summary[order(gene_summary$rank_transition),][1:min(topN, ng), ], top_trans_path)

# ============================================================
# 2) canalization: per-gene per-day replicate variance + collapse
# ============================================================
eps <- 1e-12
var_long <- data.frame(Gene=character(), Day=integer(), Var=numeric(), N=integer(), stringsAsFactors=FALSE)

for (gname in gene_names) {
  yy <- X[, gname]
  for (d0 in days_u) {
    idx <- which(day == d0)
    vals <- yy[idx]
    vals <- vals[is.finite(vals)]
    v <- if (length(vals) >= 2) stats::var(vals) else NA_real_
    var_long <- rbind(var_long, data.frame(Gene=gname, Day=as.integer(d0), Var=v, N=length(vals), stringsAsFactors=FALSE))
  }
}

var_long_path <- file.path(out_dir, "gene_day_variance_long.csv")
write_csv(var_long, var_long_path)

can <- data.frame(
  Gene = gene_names,
  var_at_bp = NA_real_,
  median_var_other = NA_real_,
  collapse_ratio = NA_real_,
  collapse_strength = NA_real_,      # -log10(collapse_ratio)
  slope_pre = gene_summary$slope_pre,
  slope_post = gene_summary$slope_post,
  stabilization_gain = NA_real_,     # |pre| - |post|
  stab_term = NA_real_,              # compressed stabilization term (log1p)
  canalization_score = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(gene_names)) {
  g <- gene_names[i]
  vv <- var_long[var_long$Gene == g, ]
  
  v_bp <- vv$Var[vv$Day == bp]
  v_bp <- if (length(v_bp) == 1) v_bp else NA_real_
  
  v_other <- vv$Var[vv$Day != bp]
  v_other <- v_other[is.finite(v_other)]
  med_other <- if (length(v_other) >= 1) stats::median(v_other) else NA_real_
  
  can$var_at_bp[i] <- v_bp
  can$median_var_other[i] <- med_other
  
  if (is.finite(v_bp) && is.finite(med_other)) {
    cr <- (v_bp + eps)/(med_other + eps)
    can$collapse_ratio[i] <- cr
    can$collapse_strength[i] <- -log10(cr)
  }
  
  pre <- can$slope_pre[i]; post <- can$slope_post[i]
  if (is.finite(pre) && is.finite(post)) {
    can$stabilization_gain[i] <- abs(pre) - abs(post)
    can$stab_term[i] <- sign(can$stabilization_gain[i]) * log1p(abs(can$stabilization_gain[i]))
  }
  
  # Score: collapse is primary; stabilization is secondary (compressed)
  can$canalization_score[i] <- can$collapse_strength[i] + 0.5 * can$stab_term[i]
}

can <- can[order(-can$canalization_score), ]
can_path <- file.path(out_dir, "canalization_gene_scores.csv")
write_csv(can, can_path)

top_canal_n <- 20
top_can <- can[1:min(top_canal_n, nrow(can)), ]
top_can_path <- file.path(out_dir, "top_canalization_genes.csv")
write_csv(top_can, top_can_path)

# ============================================================
# 3) integrated master table + explicit regime calls
# ============================================================
master <- merge(gene_summary, can, by = "Gene", all.x = TRUE)
master$mlog10p <- -log10(master$wilcox_p + 1e-300)

# add program annotation if available
if (!is.null(map)) {
  master <- merge(master, unique(map[, c("Gene","Program")]), by = "Gene", all.x = TRUE)
}

master_path <- file.path(out_dir, "gene_integrated_master.csv")
write_csv(master[order(-master$transition_score), ], master_path)

# regime calls: data-driven thresholds (quartiles) + clear semantics
q_collapse <- stats::quantile(master$collapse_ratio, probs = 0.25, na.rm = TRUE)   # strong variance collapse
q_trans    <- stats::quantile(master$transition_score, probs = 0.80, na.rm = TRUE) # top 20% transition

# effect thresholds aligned to your volcano thinking
d_cut <- 0.8
p_cut <- 2  # -log10(p) >= 2 -> p <= 0.01

reg <- master
reg$Regime <- "other"

# Base calls
reg$Regime[ reg$collapse_ratio <= q_collapse ] <- "canalized_checkpoint"
reg$Regime[ (reg$collapse_ratio <= q_collapse) & (reg$stabilization_gain > 0) ] <- "canalized_locked"

reg$Regime[ reg$transition_score >= q_trans ] <- "transition_reallocation"

reg$Regime[ (reg$cohen_d_late_minus_early >= d_cut) & (reg$mlog10p >= p_cut) ] <- "late_deployment"
reg$Regime[ (reg$cohen_d_late_minus_early <= -d_cut) & (reg$mlog10p >= p_cut) ] <- "early_program"

# Priority order (enforce): locked > checkpoint > transition > deployment > early
priority <- c("canalized_locked","canalized_checkpoint","transition_reallocation","late_deployment","early_program","other")
reg$Regime <- factor(reg$Regime, levels = priority)
reg <- reg[order(reg$Regime, reg$collapse_ratio, -reg$transition_score), ]
reg$Regime <- as.character(reg$Regime)

reg_path <- file.path(out_dir, "gene_regime_calls.csv")
write_csv(reg, reg_path)

# ============================================================
# 4) figures
# ============================================================
# volcano-ish: effect size vs evidence (label all genes; you can edit cex if cluttered)
volcano_path <- file.path(out_dir, "early_vs_late_volcanoish.png")
png(volcano_path, width = 1100, height = 900, res = 140)
par(mar = c(5,5,3,2))
x <- gene_summary$cohen_d_late_minus_early
y <- -log10(gene_summary$wilcox_p + 1e-300)
plot(x, y, pch = 16,
     xlab = "Cohen's d (late - early)",
     ylab = "-log10(Wilcoxon p)",
     main = "Early vs Late: effect size vs evidence")
text(x, y, labels = gene_summary$Gene, cex = 0.7, pos = 3)
abline(v = c(-d_cut, d_cut), lty = 2)
abline(h = p_cut, lty = 2)
dev.off()

# transition trajectories
traj_path <- file.path(out_dir, "top_transition_gene_trajectories.png")
topG <- gene_summary$Gene[order(gene_summary$rank_transition)][1:min(10, ng)]

png(traj_path, width = 1200, height = 900, res = 140)
par(mfrow = c(2,5), mar = c(4,4,2,1))
for (gname in topG) {
  yy <- X[, gname]
  plot(day, yy, pch = 16, xlab = "Day", ylab = "Expression", main = gname)
  abline(v = bp, lty = 2)
  ud <- sort(unique(day))
  m <- vapply(ud, function(d0) mean(yy[day == d0], na.rm = TRUE), numeric(1))
  lines(ud, m, lwd = 2)
}
dev.off()

# canalization variance trajectories for top canal genes
can_fig_path <- file.path(out_dir, "top_canalization_variance.png")
topG_can <- top_can$Gene[1:min(10, nrow(top_can))]

png(can_fig_path, width = 1200, height = 900, res = 140)
par(mfrow = c(2,5), mar = c(4,4,2,1))
for (gname in topG_can) {
  vv <- var_long[var_long$Gene == gname, ]
  ord <- order(vv$Day)
  plot(vv$Day, vv$Var, pch = 16, xlab = "Day", ylab = "Replicate variance", main = gname)
  abline(v = bp, lty = 2)
  lines(vv$Day[ord], vv$Var[ord], lwd = 2)
}
dev.off()

# ============================================================
# 5) optional program-level summaries (if mapping provided)
# ============================================================
program_summary_path <- NULL
program_scatter_path <- NULL
program_traj_path <- NULL

if (!is.null(map)) {
  programs <- sort(unique(map$Program))
  prog_sum <- data.frame(
    Program = programs,
    n_genes = 0L,
    mean_d = NA_real_,
    mean_transition = NA_real_,
    mean_collapse_strength = NA_real_,
    mean_canalization_score = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # program trajectories (mean expression across genes)
  program_traj_path <- file.path(out_dir, "program_trajectories.png")
  png(program_traj_path, width = 1200, height = 900, res = 140)
  par(mar = c(5,5,3,2))
  plot(NA, xlim = range(day), ylim = range(X, na.rm = TRUE),
       xlab = "Day", ylab = "Program mean expression",
       main = "Program trajectories (mean across genes)")
  abline(v = bp, lty = 2)
  
  ud <- sort(unique(day))
  
  for (i in seq_along(programs)) {
    pr <- programs[i]
    genes_pr <- map$Gene[map$Program == pr]
    jj <- match(genes_pr, gene_names)
    jj <- jj[is.finite(jj)]
    if (length(jj) < 1) next
    
    M <- rowMeans(X[, jj, drop = FALSE], na.rm = TRUE)
    m <- vapply(ud, function(d0) mean(M[day == d0], na.rm = TRUE), numeric(1))
    lines(ud, m, lwd = 2)
    text(ud[length(ud)], m[length(m)], labels = pr, cex = 0.8, pos = 4)
    
    ss <- master[master$Gene %in% genes_pr, ]
    prog_sum$n_genes[i] <- nrow(ss)
    prog_sum$mean_d[i] <- mean(ss$cohen_d_late_minus_early, na.rm = TRUE)
    prog_sum$mean_transition[i] <- mean(ss$transition_score, na.rm = TRUE)
    prog_sum$mean_collapse_strength[i] <- mean(ss$collapse_strength, na.rm = TRUE)
    prog_sum$mean_canalization_score[i] <- mean(ss$canalization_score, na.rm = TRUE)
  }
  dev.off()
  
  prog_sum <- prog_sum[order(-prog_sum$mean_transition), ]
  program_summary_path <- file.path(out_dir, "program_summary.csv")
  write_csv(prog_sum, program_summary_path)
  
  # program scatter: mean_d vs mean_transition
  program_scatter_path <- file.path(out_dir, "program_effect_scatter.png")
  png(program_scatter_path, width = 1100, height = 900, res = 140)
  par(mar = c(5,5,3,2))
  plot(prog_sum$mean_d, prog_sum$mean_transition, pch = 16,
       xlab = "Mean Cohen's d (late - early)",
       ylab = "Mean transition_score",
       main = "Program-level: deployment vs checkpoint reallocation")
  text(prog_sum$mean_d, prog_sum$mean_transition, labels = prog_sum$Program, cex = 0.8, pos = 3)
  abline(v = 0, lty = 2)
  dev.off()
}

# ============================================================
# 6) README + provenance + manifest
# ============================================================
readme_lines <- c(
  "Plot-faithful Early/Late Biology Extractor — FULL REWRITE + CANALIZATION",
  "-----------------------------------------------------------------------",
  paste0("Run name     : ", run_name),
  paste0("Timestamp    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Output dir   : ", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
  paste0("Expression   : ", basename(expr_csv)),
  paste0("Day cut      : ", day_cut, "  (EARLY <= cut; LATE > cut; fallback quantile split if needed)"),
  paste0("Breakpoint   : ", bp, "  (transition slopes + canalization variance checkpoint evaluated at bp)"),
  if (nzchar(map_path)) paste0("Mapping      : ", basename(map_path)) else "Mapping      : (none)",
  "",
  "Core outputs:",
  paste0(" - ", basename(gene_summary_path), " (early/late, slopes, transition metrics)"),
  paste0(" - ", basename(var_long_path), " (replicate variance per gene per day)"),
  paste0(" - ", basename(can_path), " (canalization metrics: collapse + stabilization + score)"),
  paste0(" - ", basename(master_path), " (integrated master table: summary + canalization + Program if given)"),
  paste0(" - ", basename(reg_path), " (explicit regime calls)"),
  paste0(" - ", basename(top_early_path)),
  paste0(" - ", basename(top_late_path)),
  paste0(" - ", basename(top_trans_path)),
  paste0(" - ", basename(top_can_path)),
  paste0(" - ", basename(volcano_path)),
  paste0(" - ", basename(traj_path)),
  paste0(" - ", basename(can_fig_path)),
  if (!is.null(program_summary_path)) paste0(" - ", basename(program_summary_path)) else NULL,
  if (!is.null(program_scatter_path)) paste0(" - ", basename(program_scatter_path)) else NULL,
  if (!is.null(program_traj_path)) paste0(" - ", basename(program_traj_path)) else NULL,
  "",
  "Interpretation keys:",
  " - cohen_d_late_minus_early > 0 => higher in late; < 0 => higher in early.",
  " - transition_score: large => strong reallocation around bp (slope change + fit support).",
  " - collapse_ratio = var(bp)/median(var(other days)): smaller => stronger variance collapse at bp (canalization).",
  " - canalization_score = collapse_strength + 0.5*compressed(stabilization_gain). Collapse dominates by design.",
  "",
  "Regime calls (gene_regime_calls.csv):",
  " - canalized_checkpoint: collapse_ratio in lowest quartile.",
  " - canalized_locked: checkpoint + stabilization_gain > 0.",
  " - transition_reallocation: transition_score in top 20%.",
  " - late_deployment / early_program: |d|>=0.8 and p<=0.01 (by -log10(p)>=2).",
  "",
  "Provenance:",
  " - Provenance.qmd (+ Provenance.html if Quarto is available)",
  " - Project_Manifest.json + Project_Manifest_Files.csv"
)

readme_path <- file.path(out_dir, "README_outputs.txt")
write_text(readme_path, readme_lines)

key_outputs <- c(
  gene_summary_path, var_long_path, can_path, master_path, reg_path,
  top_early_path, top_late_path, top_trans_path, top_can_path,
  volcano_path, traj_path, can_fig_path,
  readme_path
)
if (!is.null(program_summary_path)) key_outputs <- c(key_outputs, program_summary_path)
if (!is.null(program_scatter_path)) key_outputs <- c(key_outputs, program_scatter_path)
if (!is.null(program_traj_path)) key_outputs <- c(key_outputs, program_traj_path)

write_qmd_and_render(
  out_dir = out_dir,
  title = "Plot-faithful Early/Late Biology Extractor — Provenance",
  script_id = script_id,
  git_meta = git_meta,
  key_outputs = basename(key_outputs),
  notes = c(
    paste0("Input expression: ", basename(expr_csv)),
    paste0("Day cut: ", day_cut),
    paste0("Breakpoint: ", bp),
    if (nzchar(map_path)) paste0("Mapping: ", basename(map_path)) else "Mapping: (none)"
  )
)

inputs <- list(
  expression_csv = normalizePath(expr_csv, winslash = "/", mustWork = FALSE),
  mapping_csv = if (nzchar(map_path)) normalizePath(map_path, winslash = "/", mustWork = FALSE) else ""
)
params <- list(day_cut = day_cut, bp = bp, topN = topN, top_canal_n = top_canal_n, d_cut = d_cut, p_cut = p_cut)
write_manifest(out_dir, run_name, script_id, git_meta, inputs, params)

cat("\nDONE.\nKey outputs:\n")
for (p in key_outputs) cat(" - ", basename(p), "\n", sep = "")
cat("\n")