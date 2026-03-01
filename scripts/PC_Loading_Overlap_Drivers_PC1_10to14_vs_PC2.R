# ==========================================================
# Identify genes driving PC1(10-14) -> PC2 similarity (variance PCA)
# + cumulative contribution + first-derivative diagnostics
# Base R only
# ==========================================================

cat("\nSelect your PCA loadings CSV (e.g., JointPCA_Agroup.csv)\n")
csv_path <- file.choose()

df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)

# ---- Detect gene column ----
gene_col <- NULL
if ("Gene" %in% names(df)) gene_col <- "Gene"
if (is.null(gene_col) && "gene" %in% names(df)) gene_col <- "gene"
if (is.null(gene_col)) {
  nonnum <- names(df)[!sapply(df, is.numeric)]
  if (length(nonnum) >= 1) gene_col <- nonnum[1]
}
if (is.null(gene_col)) stop("Could not detect a gene name column. Please rename a column to 'Gene'.")

genes <- df[[gene_col]]

# ---- Helper to find columns flexibly ----
find_cols <- function(patterns) {
  hits <- rep(FALSE, ncol(df))
  for (p in patterns) hits <- hits | grepl(p, names(df), perl = TRUE)
  names(df)[hits]
}

# ---- Locate PC1 columns for 10-12 and 12-14 ----
pc1_10_12_cols <- find_cols(c("10[_\\- ]12.*PC1", "PC1.*10[_\\- ]12"))
pc1_12_14_cols <- find_cols(c("12[_\\- ]14.*PC1", "PC1.*12[_\\- ]14"))

if (length(pc1_10_12_cols) < 1 || length(pc1_12_14_cols) < 1) {
  cat("\nAvailable columns:\n")
  print(names(df))
  stop("\nCouldn't find PC1 columns for 10-12 and/or 12-14. Please check your column names.")
}

pc1_10_12 <- df[[pc1_10_12_cols[1]]]
pc1_12_14 <- df[[pc1_12_14_cols[1]]]
pc1_10_14 <- (pc1_10_12 + pc1_12_14) / 2

# ---- PC2 reference vector ----
pc2_10_12_cols <- find_cols(c("10[_\\- ]12.*PC2", "PC2.*10[_\\- ]12"))
pc2_12_14_cols <- find_cols(c("12[_\\- ]14.*PC2", "PC2.*12[_\\- ]14"))

use_pc2_centroid <- TRUE
pc2_ref <- NULL
pc2_ref_label <- NULL

if (length(pc2_10_12_cols) >= 1 && length(pc2_12_14_cols) >= 1) {
  pc2_10_12 <- df[[pc2_10_12_cols[1]]]
  pc2_12_14 <- df[[pc2_12_14_cols[1]]]
  pc2_ref <- (pc2_10_12 + pc2_12_14) / 2
  pc2_ref_label <- "PC2_10_14(avg of 10-12 and 12-14)"
  use_pc2_centroid <- FALSE
}

if (use_pc2_centroid) {
  pc2_cols <- find_cols(c("PC2"))
  if (length(pc2_cols) < 1) stop("Couldn't find any PC2 columns to build a centroid.")
  pc2_mat <- as.matrix(df[pc2_cols])
  pc2_ref <- rowMeans(pc2_mat, na.rm = TRUE)
  pc2_ref_label <- paste0("PC2_centroid(mean of ", length(pc2_cols), " PC2 columns)")
}

# ---- Choose signed vs absolute scoring ----
cat("\nHow was your clustering done?\n")
cat("1: Signed loadings (direction matters)\n")
cat("2: Absolute loadings (magnitude only)\n")
mode <- readline("Enter 1 or 2: ")
mode <- trimws(mode)
if (!(mode %in% c("1","2"))) mode <- "1"

if (mode == "1") {
  score <- pc1_10_14 * pc2_ref
  sim_cor <- suppressWarnings(cor(pc1_10_14, pc2_ref, use = "pairwise.complete.obs"))
  score_label <- "signed: PC1_10_14 * PC2_ref"
} else {
  score <- abs(pc1_10_14) * abs(pc2_ref)
  sim_cor <- suppressWarnings(cor(abs(pc1_10_14), abs(pc2_ref), use = "pairwise.complete.obs"))
  score_label <- "absolute: |PC1_10_14| * |PC2_ref|"
}

abs_score <- abs(score)

# ---- Rank genes by impact (largest AbsScore first) ----
ord <- order(abs_score, decreasing = TRUE, na.last = NA)

genes_ord    <- genes[ord]
pc1_10_12_o  <- pc1_10_12[ord]
pc1_12_14_o  <- pc1_12_14[ord]
pc1_10_14_o  <- pc1_10_14[ord]
pc2_ref_o    <- pc2_ref[ord]
score_o      <- score[ord]
abs_o        <- abs_score[ord]

# ---- Cumulative contribution ----
total_abs <- sum(abs_o, na.rm = TRUE)
cum_abs   <- cumsum(abs_o)
cum_frac  <- if (total_abs > 0) cum_abs / total_abs else rep(NA_real_, length(cum_abs))

# ---- Smooth curve + first derivative (discrete) ----
# x = rank index
x <- seq_along(abs_o)

# Smooth AbsScore vs rank using loess (robust, base R)
# span controls smoothness; 0.25–0.4 usually works well for N~150
span_use <- 0.30
lo <- loess(abs_o ~ x, span = span_use, degree = 2, family = "gaussian")
yhat <- predict(lo, x)

# First derivative dy/dx via finite differences
# derivative at midpoints; we’ll store as length N with NA at first element
dy <- c(NA_real_, diff(yhat))  # since dx=1 per rank

# Identify “elbow candidates” = places where slope magnitude drops sharply
# Using second difference of smoothed curve as a curvature proxy
d2y <- c(NA_real_, NA_real_, diff(yhat, differences = 2))

# Rank elbow candidates by largest positive curvature magnitude (|d2y|)
# (Typically elbow corresponds to strong curvature early in the curve)
cand_idx <- which(!is.na(d2y))
# take top 10 curvature points by absolute curvature
top_cand <- cand_idx[order(abs(d2y[cand_idx]), decreasing = TRUE)][1:min(10, length(cand_idx))]

# ---- Build output table ----
out <- data.frame(
  Gene = genes_ord,
  PC1_10_12 = pc1_10_12_o,
  PC1_12_14 = pc1_12_14_o,
  PC1_10_14 = pc1_10_14_o,
  PC2_ref = pc2_ref_o,
  Score = score_o,
  AbsScore = abs_o,
  CumAbsScore = cum_abs,
  CumFrac = cum_frac,
  SmoothedAbsScore = yhat,
  dAbsScore_dx = dy,
  Curvature_d2 = d2y,
  Rank = x,
  stringsAsFactors = FALSE
)

cat("\nSimilarity check:\n")
cat("PC2 reference:", pc2_ref_label, "\n")
cat("Scoring:", score_label, "\n")
cat("Correlation used as a sanity check =", signif(sim_cor, 4), "\n")
cat("LOESS span used =", span_use, "\n\n")

# ---- Write outputs ----
out_dir <- file.path(dirname(csv_path), "PC1_10_14_into_PC2_drivers")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

top_n <- 50
write.csv(out, file.path(out_dir, "AllGenes_ranked_PC1_10_14_into_PC2.csv"), row.names = FALSE)
write.csv(out[1:min(top_n, nrow(out)), ],
          file.path(out_dir, paste0("Top", top_n, "_drivers_PC1_10_14_into_PC2.csv")),
          row.names = FALSE)

# ---- Plot 1: AbsScore scree-like plot + LOESS ----
png(file.path(out_dir, "AbsScore_scree_with_elbow_candidates.png"), width = 2400, height = 1600, res = 300)
par(mar = c(5, 5, 4, 2))
plot(x, abs_o, pch = 16, cex = 0.7,
     main = "AbsScore ranked (drivers of PC1(10–14) similarity to PC2)",
     xlab = "Rank (1 = strongest driver)",
     ylab = "AbsScore = |PC1_10_14 * PC2_ref|")
lines(x, yhat, lwd = 3)
# Mark elbow candidates
points(top_cand, abs_o[top_cand], pch = 19, cex = 1.0)
text(top_cand, abs_o[top_cand], labels = top_cand, pos = 3, cex = 0.7)
mtext(paste0("PC2 ref = ", pc2_ref_label, " | scoring = ", score_label, " | loess span=", span_use),
      side = 3, line = 0.5, cex = 0.85)
legend("topright", legend = c("AbsScore", "Loess smooth", "Elbow candidates (by curvature)"),
       pch = c(16, NA, 19), lty = c(NA, 1, NA), lwd = c(NA, 3, NA), bty = "n")
dev.off()

# ---- Plot 2: Cumulative contribution ----
png(file.path(out_dir, "AbsScore_cumulative.png"), width = 2400, height = 1600, res = 300)
par(mar = c(5, 5, 4, 2))
plot(x, cum_frac, type = "l", lwd = 3,
     main = "Cumulative contribution of top drivers",
     xlab = "Rank (1 = strongest driver)",
     ylab = "Cumulative fraction of total AbsScore")
abline(h = c(0.5, 0.7, 0.8, 0.9), lty = 3)
mtext("Dashed lines at 50%, 70%, 80%, 90%", side = 3, line = 0.5, cex = 0.85)
dev.off()

# ---- Plot 3: First derivative of smoothed curve ----
png(file.path(out_dir, "AbsScore_derivative.png"), width = 2400, height = 1600, res = 300)
par(mar = c(5, 5, 4, 2))
plot(x, dy, type = "l", lwd = 3,
     main = "First derivative of smoothed AbsScore (slope vs rank)",
     xlab = "Rank",
     ylab = "d(AbsScore_smooth)/dx")
abline(h = 0, lty = 2)
# Mark elbow candidates on derivative plot too
points(top_cand, dy[top_cand], pch = 19, cex = 1.0)
text(top_cand, dy[top_cand], labels = top_cand, pos = 3, cex = 0.7)
mtext("Large slope changes early often indicate a natural cutoff region", side = 3, line = 0.5, cex = 0.85)
dev.off()

# ---- Write elbow candidates text summary ----
elbow_path <- file.path(out_dir, "Elbow_candidates.txt")
con <- file(elbow_path, open = "wt")
writeLines(c(
  "Elbow candidates identified by largest curvature (|second difference|) of LOESS-smoothed AbsScore.",
  paste0("PC2 reference: ", pc2_ref_label),
  paste0("Scoring: ", score_label),
  paste0("LOESS span: ", span_use),
  "",
  "Top curvature-ranked candidate ranks (Rank : Gene : AbsScore : CumFrac):"
), con)
for (i in top_cand) {
  line <- paste0(
    i, " : ", out$Gene[i],
    " : AbsScore=", signif(out$AbsScore[i], 5),
    " : CumFrac=", signif(out$CumFrac[i], 4)
  )
  writeLines(line, con)
}
close(con)

cat("Done.\nOutputs written to:\n", out_dir, "\n", sep = "")
