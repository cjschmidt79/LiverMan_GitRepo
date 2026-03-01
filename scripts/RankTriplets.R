# Rank triplets by effect-size stability (NOT p-values)
# Input: bootstrap_results_all.rds created by your triplet pipeline
# Output: triplet_stability_ranking.csv in the same folder as the RDS

in_rds <- file.choose()  # pick bootstrap_results_all.rds
boots <- readRDS(in_rds)

mad1 <- function(x) mad(x, constant = 1, na.rm = TRUE)

rank_one <- function(name, obj) {
  d <- obj$diff
  d <- d[is.finite(d)]
  if (length(d) < 50) return(NULL)
  
  med <- median(d, na.rm = TRUE)
  ci  <- as.numeric(quantile(d, c(0.025, 0.975), na.rm = TRUE))
  ci_w <- ci[2] - ci[1]
  
  p_pos <- mean(d > 0)
  p_neg <- mean(d < 0)
  sign_consistency <- max(p_pos, p_neg)
  
  # robust standardized effect: median / MAD
  z_stable <- if (mad1(d) > 0) med / mad1(d) else NA_real_
  
  # CI excludes zero?
  ci_excludes_0 <- (ci[1] > 0) || (ci[2] < 0)
  
  # stability score: big robust effect, consistent sign, narrow CI
  score <- if (is.finite(z_stable)) abs(z_stable) * sign_consistency / (1 + ci_w) else NA_real_
  
  # parse A_B_C from name
  parts <- strsplit(name, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  A <- parts[1]; B <- parts[2]; C <- parts[3]
  
  data.frame(
    A=A, B=B, C=C,
    n_boot=length(d),
    med_diff=med,
    ci_lower=ci[1], ci_upper=ci[2],
    ci_width=ci_w,
    sign_consistency=sign_consistency,
    z_stable=z_stable,
    ci_excludes_0=ci_excludes_0,
    stability_score=score,
    stringsAsFactors = FALSE
  )
}

rows <- mapply(rank_one, names(boots), boots, SIMPLIFY = FALSE)
rank_df <- do.call(rbind, Filter(Negate(is.null), rows))
if (is.null(rank_df) || nrow(rank_df) == 0) stop("No rankable triplets found (maybe too few valid bootstraps).")

rank_df <- rank_df[order(-rank_df$stability_score), ]

# ---- 2-line HUB annotation (frequency across bootstrapped triplets) ----
hub_thr <- as.numeric(quantile(table(c(rank_df$A, rank_df$B, rank_df$C)), 0.90))
rank_df$is_hub_triplet <- apply(rank_df[, c("A","B","C")], 1, function(x) any(table(c(rank_df$A, rank_df$B, rank_df$C))[x] >= hub_thr))

out_csv <- file.path(dirname(in_rds), "triplet_stability_ranking.csv")
write.csv(rank_df, out_csv, row.names = FALSE)

cat("Wrote:", out_csv, "\n")
print(head(rank_df, 20))