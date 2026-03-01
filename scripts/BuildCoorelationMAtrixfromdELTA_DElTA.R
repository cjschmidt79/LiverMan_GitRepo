# Build correlation-distance matrix from Δ (day-to-day change) profiles
# Input: DailySummary_ByDayGene.csv with columns: Day, Source, daily_value
# must be in same directory as infile:DailySummary_ByDayGene.csv
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

infile <- "DailySummary_ByDayGene.csv"   # adjust path if needed
method <- "spearman"                    # "spearman" or "pearson"

df <- read.csv(infile, stringsAsFactors = FALSE)

# Wide matrix: genes x days
mat <- df %>%
  mutate(Day = as.integer(Day),
         Source = as.character(Source)) %>%
  select(Source, Day, daily_value) %>%
  pivot_wider(names_from = Day, values_from = daily_value)

gene_ids <- mat$Source
X <- as.matrix(mat[, -1, drop = FALSE])

# Ensure days are in increasing order
day_cols <- as.integer(colnames(X))
ord <- order(day_cols)
X <- X[, ord, drop = FALSE]
day_cols <- day_cols[ord]

# Δ profiles: genes x transitions (diff across consecutive days)
D <- t(apply(X, 1, function(v) diff(v)))

colnames(D) <- paste0("d", day_cols[-1], "_minus_d", day_cols[-length(day_cols)])
rownames(D) <- gene_ids

# Drop genes with undefined Δ profiles (all NA or zero variance across transitions)
ok <- apply(D, 1, function(v) {
  v <- v[is.finite(v)]
  length(v) >= 3 && stats::sd(v) > 0
})
D2 <- D[ok, , drop = FALSE]

# Correlation similarity (genes x genes), then correlation distance
C <- cor(t(D2), method = method, use = "pairwise.complete.obs")
dist_mat <- as.dist(1 - C)

# Outputs
saveRDS(dist_mat, file = paste0("CorrDistance_DeltaProfiles_", method, ".rds"))
write.csv(as.matrix(dist_mat),
          file = paste0("CorrDistance_DeltaProfiles_", method, ".csv"),
          row.names = TRUE)

cat("Wrote distance matrix for", nrow(D2), "genes using", method, "correlation.\n")
