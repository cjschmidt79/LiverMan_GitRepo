#!/usr/bin/env Rscript
# ============================================================
# Triplet Stability Report (HTML) — Top-scored / stable triplets
# Base R, SOURCE-to-run friendly
# Outputs: outputs/<run_name>_<timestamp>/
# ============================================================

# -----------------------------
# Small utilities
# -----------------------------
safe_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (nchar(x) == 0) "run" else x
}

timestamp_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

capture_script_identity <- function() {
  # Best-effort: works for Rscript, source(), and RStudio
  path <- NA_character_
  # 1) commandArgs (Rscript)
  ca <- commandArgs(trailingOnly = FALSE)
  fi <- grep("^--file=", ca, value = TRUE)
  if (length(fi) == 1) path <- sub("^--file=", "", fi)
  
  # 2) sys.frames trick (source)
  if (is.na(path)) {
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)
      if (!is.null(f$ofile)) { path <- f$ofile; break }
    }
  }
  
  # 3) RStudio API if available
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

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
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

# -----------------------------
# User inputs (interactive)
# -----------------------------
cat("\nSelect triplet_stability_ranking.csv\n")
csv_path <- file.choose()
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

# ---- knobs (edit as desired) ----
run_name <- "TripletStabilityReport"
stability_min <- 2
topN <- 40              # set to Inf to include all passing stability_min
sort_by <- "stability"  # "stability" or "abs_effect"

# -----------------------------
# Validate
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

if (tolower(sort_by) == "abs_effect") {
  keep <- keep[order(abs(keep$med_diff), decreasing = TRUE), , drop = FALSE]
} else {
  keep <- keep[order(keep$stability_score, decreasing = TRUE), , drop = FALSE]
}

if (is.finite(topN) && nrow(keep) > topN) keep <- keep[1:topN, , drop = FALSE]

# -----------------------------
# Output directory
# -----------------------------
ts <- timestamp_now()
out_dir <- file.path(getwd(), "outputs", paste0(safe_slug(run_name), "_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Copy input for provenance
file.copy(csv_path, file.path(out_dir, basename(csv_path)), overwrite = TRUE)

# Capture identity/provenance
prov <- capture_script_identity()
prov_path <- file.path(out_dir, "provenance.txt")
writeLines(c(
  "=== Triplet Stability Report provenance ===",
  paste0("timestamp: ", prov$timestamp),
  paste0("script_path: ", prov$script_path),
  paste0("script_name: ", prov$script_name),
  paste0("working_dir: ", prov$wd),
  paste0("R_version: ", prov$r_version),
  paste0("platform: ", prov$platform),
  paste0("input_csv: ", csv_path),
  paste0("stability_min: ", stability_min),
  paste0("topN: ", topN),
  paste0("sort_by: ", sort_by),
  ""
), prov_path)

# -----------------------------
# Save filtered table
# -----------------------------
keep_csv <- file.path(out_dir, "top_scored_triplets.csv")
write.csv(keep, keep_csv, row.names = FALSE)

# -----------------------------
# Plots (base R)
# -----------------------------
# 1) Stability bar (topN)
bar_png <- file.path(out_dir, "plot_stability_bar.png")
png(bar_png, width = 1400, height = max(700, 22 * nrow(keep)), res = 150)
op <- par(no.readonly = TRUE)
par(mar = c(5, 22, 3, 2))
y <- rev(seq_len(nrow(keep)))
plot(keep$stability_score, y,
     yaxt = "n", ylab = "", xlab = "stability_score",
     pch = 16)
axis(2, at = y, labels = rev(keep$triplet), las = 2, cex.axis = 0.75)
mtext("Top scored triplets (ranked)", side = 3, line = 1, adj = 0, cex = 0.9)
par(op)
dev.off()

# 2) Forest plot of effect (med_diff + CI)
forest_png <- file.path(out_dir, "plot_forest_med_diff_ci.png")
png(forest_png, width = 1400, height = max(700, 22 * nrow(keep)), res = 150)
op <- par(no.readonly = TRUE)
par(mar = c(5, 22, 3, 2))
y <- rev(seq_len(nrow(keep)))
xlim <- range(c(keep$ci_lower, keep$ci_upper, 0), finite = TRUE)
plot(keep$med_diff, y,
     xlim = xlim,
     yaxt = "n", ylab = "", xlab = "med_diff with 95% CI",
     pch = 16)
axis(2, at = y, labels = rev(keep$triplet), las = 2, cex.axis = 0.75)
segments(keep$ci_lower, y, keep$ci_upper, y, lwd = 2)
abline(v = 0, lty = 2)
mtext(sprintf("Filtered: stability_score > %s | showing %d triplets", stability_min, nrow(keep)),
      side = 3, line = 1, adj = 0, cex = 0.9)
par(op)
dev.off()

# -----------------------------
# Simple self-contained HTML report (no Quarto required)
# -----------------------------
html_path <- file.path(out_dir, "Triplet_Stability_Report.html")

# Build an HTML table
table_rows <- apply(keep, 1, function(r) {
  # r is character; use names to pick columns
  cols <- c("A","B","C","stability_score","med_diff","ci_lower","ci_upper")
  tds <- vapply(cols, function(cc) {
    val <- r[[cc]]
    paste0("<td>", html_escape(val), "</td>")
  }, character(1))
  paste0("<tr>", paste(tds, collapse = ""), "</tr>")
})

html <- c(
  "<!doctype html>",
  "<html><head><meta charset='utf-8'>",
  "<title>Triplet Stability Report</title>",
  "<style>",
  "body{font-family:Arial,Helvetica,sans-serif;margin:24px;max-width:1100px}",
  "h1,h2{margin-bottom:0.2em}",
  ".meta{color:#444;font-size:0.95em;margin-top:0}",
  "table{border-collapse:collapse;width:100%;margin:14px 0}",
  "th,td{border:1px solid #ddd;padding:6px 8px;font-size:0.95em}",
  "th{background:#f5f5f5;text-align:left;position:sticky;top:0}",
  "img{max-width:100%;border:1px solid #eee;margin:10px 0}",
  "code{background:#f7f7f7;padding:2px 4px;border-radius:4px}",
  "</style></head><body>",
  "<h1>Triplet Stability Report</h1>",
  sprintf("<p class='meta'><b>Generated:</b> %s<br><b>Input:</b> %s<br><b>Filter:</b> stability_score &gt; %s<br><b>Shown:</b> top %s (after filter), sorted by %s</p>",
          html_escape(prov$timestamp),
          html_escape(basename(csv_path)),
          html_escape(as.character(stability_min)),
          html_escape(as.character(topN)),
          html_escape(sort_by)),
  "<h2>Plots</h2>",
  sprintf("<h3>Stability scores</h3><img src='%s' alt='stability bar'>", basename(bar_png)),
  sprintf("<h3>Effect size (med_diff) with 95%% CI</h3><img src='%s' alt='forest plot'>", basename(forest_png)),
  "<h2>Top triplets table</h2>",
  "<p>Columns: <code>A</code>, <code>B</code>, <code>C</code>, <code>stability_score</code>, <code>med_diff</code>, <code>ci_lower</code>, <code>ci_upper</code></p>",
  "<div style='max-height:520px;overflow:auto;border:1px solid #eee'>",
  "<table><thead><tr>",
  "<th>A</th><th>B</th><th>C</th><th>stability_score</th><th>med_diff</th><th>ci_lower</th><th>ci_upper</th>",
  "</tr></thead><tbody>",
  table_rows,
  "</tbody></table></div>",
  "<h2>Files</h2>",
  "<ul>",
  sprintf("<li><b>Filtered table:</b> %s</li>", html_escape(basename(keep_csv))),
  sprintf("<li><b>Provenance:</b> %s</li>", html_escape(basename(prov_path))),
  "</ul>",
  "</body></html>"
)

writeLines(unlist(html), html_path)

# Inventory
inv_csv <- write_inventory(out_dir)

# README
readme_path <- file.path(out_dir, "README_outputs.txt")
writeLines(c(
  "Triplet Stability Report — outputs",
  "---------------------------------",
  paste0("Output dir: ", out_dir),
  "",
  "Key files:",
  paste0("- Triplet_Stability_Report.html  (open in browser)"),
  paste0("- top_scored_triplets.csv        (filtered/ranked table)"),
  paste0("- plot_stability_bar.png         (stability plot)"),
  paste0("- plot_forest_med_diff_ci.png    (forest plot)"),
  paste0("- provenance.txt                 (script+run metadata)"),
  paste0("- Project_Manifest_Files.csv     (inventory of files)"),
  ""
), readme_path)

cat("\nDONE.\nOutput directory:\n", out_dir, "\n\nOpen:\n", html_path, "\n", sep = "")
cat("\nInventory:\n", inv_csv, "\n", sep = "")