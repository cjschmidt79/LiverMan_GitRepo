#!/usr/bin/env Rscript
# ============================================================
# Triplet Report Assembler (HTML)
#  - Filters to top stable triplets (e.g., stability_score > 2)
#  - Embeds triplet stats + existing permutation plot PNGs
#  - ADDED: Provenance metadata block embedded in the MAIN HTML
# Base R. SOURCE-to-run in RStudio.
# Outputs: outputs/<run_name>_<timestamp>/
# ============================================================

# -----------------------------
# Utilities
# -----------------------------
safe_slug <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[nchar(x) == 0] <- "x"
  x
}

timestamp_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

capture_script_identity <- function() {
  path <- NA_character_
  ca <- commandArgs(trailingOnly = FALSE)
  fi <- grep("^--file=", ca, value = TRUE)
  if (length(fi) == 1) path <- sub("^--file=", "", fi)
  
  if (is.na(path)) {
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)
      if (!is.null(f$ofile)) { path <- f$ofile; break }
    }
  }
  
  if (is.na(path) && requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- try(rstudioapi::getSourceEditorContext(), silent = TRUE)
    if (!inherits(ctx, "try-error")) {
      if (!is.null(ctx$path) && nzchar(ctx$path)) path <- ctx$path
    }
  }
  
  list(
    script_path = path,
    script_name = if (!is.na(path)) basename(path) else NA_character_,
    wd = getwd(),
    timestamp = as.character(Sys.time()),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform
  )
}

write_inventory <- function(out_dir) {
  files <- list.files(out_dir, full.names = TRUE, recursive = TRUE)
  info <- file.info(files)
  inv <- data.frame(
    path = files,
    size_bytes = info$size,
    mtime = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
  inv_csv <- file.path(out_dir, "Project_Manifest_Files.csv")
  write.csv(inv, inv_csv, row.names = FALSE)
  inv_csv
}

# --- NEW: simple command runner (for git) ---
run_cmd <- function(cmd) {
  out <- tryCatch(system(cmd, intern = TRUE), warning = function(w) NA_character_, error = function(e) NA_character_)
  if (length(out) == 0 || all(is.na(out))) return(NA_character_)
  paste(out, collapse = "\n")
}

# --- NEW: git provenance (best-effort; safe if not in repo) ---
git_provenance <- function(start_dir = getwd()) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  
  ok <- TRUE
  tryCatch(setwd(start_dir), error = function(e) ok <<- FALSE)
  if (!ok) {
    return(list(
      git_repo_root = NA_character_,
      git_branch = NA_character_,
      git_commit = NA_character_,
      git_describe = NA_character_,
      git_is_dirty = NA_character_
    ))
  }
  
  root <- run_cmd("git rev-parse --show-toplevel 2> /dev/null")
  if (is.na(root) || !nzchar(root)) {
    return(list(
      git_repo_root = NA_character_,
      git_branch = NA_character_,
      git_commit = NA_character_,
      git_describe = NA_character_,
      git_is_dirty = NA_character_
    ))
  }
  
  setwd(root)
  branch <- run_cmd("git rev-parse --abbrev-ref HEAD 2> /dev/null")
  commit <- run_cmd("git rev-parse HEAD 2> /dev/null")
  describe <- run_cmd("git describe --always --dirty --tags 2> /dev/null")
  dirty <- run_cmd("git status --porcelain 2> /dev/null")
  is_dirty <- if (!is.na(dirty) && nzchar(dirty)) "YES" else "NO"
  
  list(
    git_repo_root = root,
    git_branch = branch,
    git_commit = commit,
    git_describe = describe,
    git_is_dirty = is_dirty
  )
}

# -----------------------------
# Interactive inputs
# -----------------------------
cat("\nSelect triplet_stability_ranking.csv\n")
csv_path <- file.choose()
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

cat("\nSelect the folder that contains the permutation plot PNGs\n")
perm_dir <- if (requireNamespace("rstudioapi", quietly = TRUE)) {
  tryCatch(rstudioapi::selectDirectory(caption = "Choose permutation plots folder"),
           error = function(e) NULL)
} else NULL
if (is.null(perm_dir) || is.na(perm_dir) || !nzchar(perm_dir)) {
  cat("RStudio picker unavailable or canceled. Please paste permutation plots folder path:\n")
  perm_dir <- normalizePath(readLines("stdin", n = 1), mustWork = TRUE)
} else {
  perm_dir <- normalizePath(perm_dir, mustWork = TRUE)
}

# -----------------------------
# Knobs (edit)
# -----------------------------
run_name <- "TripletPermutation_HTML"
stability_min <- 2
topN <- 40              # set to Inf to include all passing stability_min
sort_by <- "stability"  # "stability" or "abs_effect"

# How to match plot files:
#   "strict"  = require all three genes in filename (recommended if your filenames include A,B,C)
#   "relaxed" = accept any two genes match if strict finds nothing
match_mode <- "strict"

# Optional: include a "global permutation plot" if present
# If your original script outputs a single overall perm plot, set a regex here.
global_perm_regex <- "(global|overall|null|permutation).*\\.png$"

# -----------------------------
# Validate columns
# -----------------------------
need <- c("A","B","C","stability_score","med_diff","ci_lower","ci_upper")
miss <- setdiff(need, names(df))
if (length(miss)) stop("Missing columns in CSV: ", paste(miss, collapse = ", "))

# -----------------------------
# Filter + rank
# -----------------------------
keep <- df[df$stability_score > stability_min, , drop = FALSE]
if (nrow(keep) == 0) stop("No triplets pass stability_score > ", stability_min)

keep$triplet <- paste0(keep$A, " | ", keep$B, " | ", keep$C)
keep$triplet_slug <- safe_slug(paste0(keep$A, "_", keep$B, "_", keep$C))
keep$triplet_slug <- make.unique(keep$triplet_slug, sep = "_")

if (tolower(sort_by) == "abs_effect") {
  keep <- keep[order(abs(keep$med_diff), decreasing = TRUE), , drop = FALSE]
} else {
  keep <- keep[order(keep$stability_score, decreasing = TRUE), , drop = FALSE]
}
if (is.finite(topN) && nrow(keep) > topN) keep <- keep[1:topN, , drop = FALSE]

# -----------------------------
# Find permutation plot files
# -----------------------------
perm_files <- list.files(perm_dir, pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE)
if (length(perm_files) == 0) stop("No PNG files found in permutation plot folder: ", perm_dir)

# helper: does filename contain token (case-insensitive)?
has_tok <- function(fname, tok) grepl(tok, fname, ignore.case = TRUE)

find_triplet_plots <- function(A, B, C, files, mode = "strict") {
  # Strict: require all 3 gene symbols appear somewhere in filename
  strict_hits <- files[has_tok(basename(files), A) & has_tok(basename(files), B) & has_tok(basename(files), C)]
  if (length(strict_hits) > 0) return(strict_hits)
  
  if (tolower(mode) == "relaxed") {
    # Relaxed: accept any 2-of-3 (useful if filenames abbreviate)
    hits2 <- files[
      (has_tok(basename(files), A) & has_tok(basename(files), B)) |
        (has_tok(basename(files), A) & has_tok(basename(files), C)) |
        (has_tok(basename(files), B) & has_tok(basename(files), C))
    ]
    return(hits2)
  }
  
  character(0)
}

# Map each triplet to matching plots
plot_map <- vector("list", nrow(keep))
for (i in seq_len(nrow(keep))) {
  A <- keep$A[i]; B <- keep$B[i]; C <- keep$C[i]
  hits <- find_triplet_plots(A, B, C, perm_files, mode = match_mode)
  plot_map[[i]] <- hits
}

# Global plot (optional)
global_hits <- perm_files[grepl(global_perm_regex, basename(perm_files), ignore.case = TRUE)]
global_hit <- if (length(global_hits) >= 1) global_hits[1] else NA_character_

# -----------------------------
# Output directory + provenance
# -----------------------------
ts <- timestamp_now()
out_dir <- file.path(getwd(), "outputs", paste0(safe_slug(run_name), "_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prov <- capture_script_identity()

# --- NEW: git provenance collected now (best-effort) ---
# Prefer script directory (most likely inside repo); fall back to wd.
git_start <- if (!is.na(prov$script_path) && nzchar(prov$script_path)) dirname(prov$script_path) else prov$wd
git <- git_provenance(start_dir = git_start)

# Copy input CSV for provenance
file.copy(csv_path, file.path(out_dir, basename(csv_path)), overwrite = TRUE)

# Copy plots into output (so HTML is portable)
plots_out_dir <- file.path(out_dir, "plots")
dir.create(plots_out_dir, showWarnings = FALSE)

copy_if_exists <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) return(character(0))
  out_paths <- file.path(plots_out_dir, basename(ok))
  file.copy(ok, out_paths, overwrite = TRUE)
  out_paths
}

global_out <- if (!is.na(global_hit)) copy_if_exists(global_hit) else character(0)

triplet_plot_out <- vector("list", nrow(keep))
for (i in seq_len(nrow(keep))) {
  triplet_plot_out[[i]] <- copy_if_exists(plot_map[[i]])
}

# Save filtered table
keep_csv <- file.path(out_dir, "top_scored_triplets.csv")
write.csv(keep, keep_csv, row.names = FALSE)

# Provenance text
prov_path <- file.path(out_dir, "provenance.txt")
writeLines(c(
  "=== Triplet Permutation HTML Report provenance ===",
  paste0("timestamp: ", prov$timestamp),
  paste0("script_path: ", prov$script_path),
  paste0("script_name: ", prov$script_name),
  paste0("working_dir: ", prov$wd),
  paste0("R_version: ", prov$r_version),
  paste0("platform: ", prov$platform),
  paste0("input_csv: ", csv_path),
  paste0("perm_dir: ", perm_dir),
  paste0("stability_min: ", stability_min),
  paste0("topN: ", topN),
  paste0("sort_by: ", sort_by),
  paste0("match_mode: ", match_mode),
  paste0("global_perm_regex: ", global_perm_regex),
  "",
  "=== Git provenance (best-effort) ===",
  paste0("git_repo_root: ", git$git_repo_root),
  paste0("git_branch: ", git$git_branch),
  paste0("git_commit: ", git$git_commit),
  paste0("git_describe: ", git$git_describe),
  paste0("git_is_dirty: ", git$git_is_dirty),
  ""
), prov_path)

# -----------------------------
# Build HTML
# -----------------------------
html_path <- file.path(out_dir, "Triplet_Permutation_Report.html")

# Table rows
table_rows <- apply(keep, 1, function(r) {
  cols <- c("A","B","C","stability_score","med_diff","ci_lower","ci_upper")
  tds <- vapply(cols, function(cc) paste0("<td>", html_escape(r[[cc]]), "</td>"), character(1))
  paste0("<tr>", paste(tds, collapse = ""), "</tr>")
})

# Per-triplet plot sections
triplet_sections <- character(0)
for (i in seq_len(nrow(keep))) {
  trip <- keep$triplet[i]
  slug <- keep$triplet_slug[i]
  imgs <- triplet_plot_out[[i]]
  
  if (length(imgs) == 0) {
    img_html <- "<p class='warn'>No permutation plot PNG matched this triplet (based on filename matching).</p>"
  } else {
    img_tags <- vapply(imgs, function(p) {
      rp <- file.path("plots", basename(p))
      paste0("<img src='", rp, "' alt='", html_escape(trip), " plot'>")
    }, character(1))
    img_html <- paste(img_tags, collapse = "\n")
  }
  
  stats_line <- sprintf(
    "<p class='meta'><b>stability_score</b>: %s &nbsp; | &nbsp; <b>med_diff</b>: %s &nbsp; | &nbsp; <b>CI</b>: [%s, %s]</p>",
    html_escape(keep$stability_score[i]),
    html_escape(keep$med_diff[i]),
    html_escape(keep$ci_lower[i]),
    html_escape(keep$ci_upper[i])
  )
  
  triplet_sections <- c(
    triplet_sections,
    sprintf("<div class='card' id='%s'><h3>%s</h3>%s%s</div>", slug, html_escape(trip), stats_line, img_html)
  )
}

# Global plot section
global_section <- if (length(global_out) == 1) {
  paste0(
    "<div class='card'><h2>Global permutation plot</h2>",
    "<img src='plots/", basename(global_out), "' alt='global permutation plot'>",
    "</div>"
  )
} else {
  "<div class='card'><h2>Global permutation plot</h2><p class='meta'>No global permutation PNG detected (regex didn’t match any file).</p></div>"
}

# --- NEW: Provenance block embedded into main HTML (details/summary) ---
provenance_block <- c(
  "<div class='card'>",
  "<details open>",
  "<summary><b>Provenance / Run Metadata (click to collapse)</b></summary>",
  "<div style='margin-top:10px'>",
  
  "<h3 style='margin:10px 0 6px 0'>Run</h3>",
  "<table><tbody>",
  sprintf("<tr><th>Generated (timestamp)</th><td>%s</td></tr>", html_escape(prov$timestamp)),
  sprintf("<tr><th>Script name</th><td>%s</td></tr>", html_escape(prov$script_name)),
  sprintf("<tr><th>Script path</th><td>%s</td></tr>", html_escape(prov$script_path)),
  sprintf("<tr><th>Working directory</th><td>%s</td></tr>", html_escape(prov$wd)),
  sprintf("<tr><th>R version</th><td>%s</td></tr>", html_escape(prov$r_version)),
  sprintf("<tr><th>Platform</th><td>%s</td></tr>", html_escape(prov$platform)),
  "</tbody></table>",
  
  "<h3 style='margin:14px 0 6px 0'>Inputs</h3>",
  "<table><tbody>",
  sprintf("<tr><th>triplet_stability_ranking.csv</th><td>%s</td></tr>", html_escape(csv_path)),
  sprintf("<tr><th>Permutation plots folder</th><td>%s</td></tr>", html_escape(perm_dir)),
  "</tbody></table>",
  
  "<h3 style='margin:14px 0 6px 0'>Parameters</h3>",
  "<table><tbody>",
  sprintf("<tr><th>stability_min</th><td>%s</td></tr>", html_escape(stability_min)),
  sprintf("<tr><th>topN</th><td>%s</td></tr>", html_escape(topN)),
  sprintf("<tr><th>sort_by</th><td>%s</td></tr>", html_escape(sort_by)),
  sprintf("<tr><th>match_mode</th><td>%s</td></tr>", html_escape(match_mode)),
  sprintf("<tr><th>global_perm_regex</th><td>%s</td></tr>", html_escape(global_perm_regex)),
  "</tbody></table>",
  
  "<h3 style='margin:14px 0 6px 0'>Outputs</h3>",
  "<table><tbody>",
  sprintf("<tr><th>Output directory</th><td>%s</td></tr>", html_escape(out_dir)),
  sprintf("<tr><th>Main HTML report</th><td>%s</td></tr>", html_escape(html_path)),
  sprintf("<tr><th>Filtered table</th><td>%s</td></tr>", html_escape(keep_csv)),
  sprintf("<tr><th>Provenance text file</th><td>%s</td></tr>", html_escape(prov_path)),
  sprintf("<tr><th>Plots directory</th><td>%s</td></tr>", html_escape(plots_out_dir)),
  "</tbody></table>",
  
  "<h3 style='margin:14px 0 6px 0'>Git (best-effort)</h3>",
  "<table><tbody>",
  sprintf("<tr><th>Repo root</th><td>%s</td></tr>", html_escape(git$git_repo_root)),
  sprintf("<tr><th>Branch</th><td>%s</td></tr>", html_escape(git$git_branch)),
  sprintf("<tr><th>Commit</th><td>%s</td></tr>", html_escape(git$git_commit)),
  sprintf("<tr><th>Describe</th><td>%s</td></tr>", html_escape(git$git_describe)),
  sprintf("<tr><th>Working tree dirty?</th><td>%s</td></tr>", html_escape(git$git_is_dirty)),
  "</tbody></table>",
  
  "<p class='meta'>A plain-text copy is also written to <code>provenance.txt</code>.</p>",
  "</div>",
  "</details>",
  "</div>"
)

html <- c(
  "<!doctype html>",
  "<html><head><meta charset='utf-8'>",
  "<title>Triplet Permutation Report</title>",
  "<style>",
  "body{font-family:Arial,Helvetica,sans-serif;margin:24px;max-width:1200px}",
  "h1{margin-bottom:0.2em}",
  ".meta{color:#444;font-size:0.95em;margin-top:0.2em}",
  ".warn{color:#a33}",
  ".grid{display:grid;grid-template-columns:1fr;gap:14px}",
  ".card{border:1px solid #e6e6e6;border-radius:10px;padding:14px;background:#fff}",
  "table{border-collapse:collapse;width:100%;margin:10px 0}",
  "th,td{border:1px solid #ddd;padding:6px 8px;font-size:0.95em;vertical-align:top}",
  "th{background:#f5f5f5;text-align:left}",
  "img{max-width:100%;border:1px solid #eee;margin:10px 0;border-radius:8px}",
  "code{background:#f7f7f7;padding:2px 4px;border-radius:4px}",
  "summary{cursor:pointer}",
  "</style></head><body>",
  
  "<h1>Triplet Permutation Report</h1>",
  sprintf("<p class='meta'><b>Generated:</b> %s<br><b>Input CSV:</b> %s<br><b>Permutation plots folder:</b> %s<br><b>Filter:</b> stability_score &gt; %s<br><b>Showing:</b> %s triplets (sorted by %s)</p>",
          html_escape(prov$timestamp),
          html_escape(basename(csv_path)),
          html_escape(perm_dir),
          html_escape(as.character(stability_min)),
          html_escape(as.character(nrow(keep))),
          html_escape(sort_by)),
  
  provenance_block,
  
  global_section,
  
  "<div class='card'>",
  "<h2>Top triplets table</h2>",
  "<p class='meta'>Columns: <code>A</code>, <code>B</code>, <code>C</code>, <code>stability_score</code>, <code>med_diff</code>, <code>ci_lower</code>, <code>ci_upper</code></p>",
  "<div style='max-height:520px;overflow:auto;border:1px solid #eee'>",
  "<table><thead><tr>",
  "<th>A</th><th>B</th><th>C</th><th>stability_score</th><th>med_diff</th><th>ci_lower</th><th>ci_upper</th>",
  "</tr></thead><tbody>",
  table_rows,
  "</tbody></table></div>",
  "</div>",
  
  "<div class='card'>",
  "<h2>Triplet permutation plots</h2>",
  "<p class='meta'>Below: each triplet’s stats + any matched permutation plot PNGs copied into <code>plots/</code>.</p>",
  "</div>",
  
  "<div class='grid'>",
  triplet_sections,
  "</div>",
  
  "<div class='card'><h2>Files</h2><ul>",
  sprintf("<li><b>HTML report:</b> %s</li>", html_escape(basename(html_path))),
  sprintf("<li><b>Filtered triplets table:</b> %s</li>", html_escape(basename(keep_csv))),
  sprintf("<li><b>Provenance:</b> %s</li>", html_escape(basename(prov_path))),
  "</ul></div>",
  
  "</body></html>"
)

writeLines(unlist(html), html_path)

# Inventory + README
inv_csv <- write_inventory(out_dir)

readme_path <- file.path(out_dir, "README_outputs.txt")
writeLines(c(
  "Triplet Permutation HTML Report — outputs",
  "----------------------------------------",
  paste0("Output dir: ", out_dir),
  "",
  "Key files:",
  "- Triplet_Permutation_Report.html  (open in browser)",
  "- top_scored_triplets.csv          (filtered/ranked table)",
  "- plots/                           (copied permutation plot PNGs)",
  "- provenance.txt                   (script+run metadata)",
  "- Project_Manifest_Files.csv       (inventory of outputs)",
  "",
  "NOTE: Plot matching is filename-based. If some triplets show 'No permutation plot matched',",
  "rename plot PNGs to include the gene symbols (A, B, C) or switch match_mode to 'relaxed'.",
  ""
), readme_path)

cat("\nDONE.\nOutput directory:\n", out_dir, "\n\nOpen:\n", html_path, "\n", sep = "")
cat("\nInventory:\n", inv_csv, "\n", sep = "")