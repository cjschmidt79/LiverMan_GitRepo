# ==========================================================
# GeneRole_Integrated_Table.R
# Build ONE row per gene with:
#   - DTW cluster
#   - Max PC1 loading window (from WIDE JointPCA_Agroup.csv)
#   - Ever a bridging gene? (YES/NO)
#   - Pre- vs post-14 variance ratio (post>=16 / pre<=12)
#   - Functional category (broad)
#
# Base R only | SOURCE-to-run | Interactive file pickers
# ==========================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# Helper: stop with columns list
# -----------------------------
stop_cols <- function(df, msg) {
  stop(paste0(msg, "\nColumns are:\n  ", paste(names(df), collapse = ", ")))
}

# -----------------------------
# Helper: detect a column by candidates
# -----------------------------
pick_col <- function(df, candidates) {
  for (cc in candidates) if (cc %in% names(df)) return(cc)
  return(NULL)
}

# -----------------------------
# 1) DTW cluster assignments
# -----------------------------
cat("\nSelect DTW cluster assignment file (e.g., DTW_partitional_cluster_assignments.csv)\n")
dtw_path <- file.choose()
dtw <- read.csv(dtw_path, check.names = FALSE)

# Detect columns
dtw_gene_col <- pick_col(dtw, c("Gene", "gene", "GENE"))
dtw_cluster_col <- pick_col(dtw, c("Cluster", "cluster", "cluster_id", "ClusterID", "DTW_cluster", "DTWcluster"))
if (is.null(dtw_gene_col) || is.null(dtw_cluster_col)) {
  stop_cols(dtw, "DTW file must contain gene + cluster columns (e.g., Gene, Cluster).")
}

dtw_df <- dtw[, c(dtw_gene_col, dtw_cluster_col)]
names(dtw_df) <- c("Gene", "DTW_cluster")
dtw_df$Gene <- as.character(dtw_df$Gene)

cat("\nDTW mapping:\n  Gene column: ", dtw_gene_col, "\n  Cluster column: ", dtw_cluster_col, "\n", sep = "")

# -----------------------------
# 2) PCA loadings (WIDE format)
#    JointPCA_Agroup.csv: first col is gene; PC1 columns end with 'PC1'
# -----------------------------
cat("\nSelect PCA loadings file (e.g., JointPCA_Agroup.csv)\n")
pca_path <- file.choose()
pca <- read.csv(pca_path, check.names = FALSE)

# First column is gene in your JointPCA_Agroup.csv
names(pca)[1] <- "Gene"
pca$Gene <- as.character(pca$Gene)

pc1_cols <- grep("PC1$", names(pca), value = TRUE)
if (length(pc1_cols) == 0) stop_cols(pca, "PCA file has no columns ending in 'PC1'.")

# Abs loadings matrix
abs_mat <- abs(as.matrix(pca[, pc1_cols, drop = FALSE]))

# Handle potential non-numeric PC columns
if (!is.numeric(abs_mat)) {
  suppressWarnings(abs_mat <- abs(apply(pca[, pc1_cols, drop = FALSE], 2, as.numeric)))
}
if (anyNA(abs_mat)) {
  # still OK, but warn
  cat("\nWARNING: Some PC1 loadings are NA after numeric conversion.\n")
}

imax <- max.col(abs_mat, ties.method = "first")
Max_PC1_Window <- sub("PC1$", "", pc1_cols[imax])

pc1_df <- data.frame(
  Gene = pca$Gene,
  Max_PC1_Window = Max_PC1_Window,
  stringsAsFactors = FALSE
)

cat("\nPCA mapping:\n  PC1 columns detected: ", length(pc1_cols), "\n", sep = "")

# -----------------------------
# 3) Bridging genes summary (ever/never)
# -----------------------------
cat("\nSelect BridgingGenes_SummaryAcrossDays.csv\n")
bridge_path <- file.choose()
bridge <- read.csv(bridge_path, check.names = FALSE)

bridge_gene_col <- pick_col(bridge, c("Gene", "gene", "GENE"))
if (is.null(bridge_gene_col)) stop_cols(bridge, "Bridging summary must contain a Gene column.")

bridge_df <- data.frame(
  Gene = unique(as.character(bridge[[bridge_gene_col]])),
  Ever_Bridging = "YES",
  stringsAsFactors = FALSE
)

cat("\nBridging mapping:\n  Gene column: ", bridge_gene_col, "\n", sep = "")

# -----------------------------
# 4) Variance table -> Pre/Post ratio
#    Input: MeanVariance_GeneLevel_ByDay.csv (column names may vary)
#    We standardize to: var$Gene, var$Day, var$VarianceValue
# -----------------------------
cat("\nSelect MeanVariance_GeneLevel_ByDay.csv\n")
var_path <- file.choose()
var <- read.csv(var_path, check.names = FALSE)

# Detect required columns
var_gene_col <- pick_col(var, c("Gene", "gene", "GENE"))
var_day_col  <- pick_col(var, c("Day", "day", "DAY", "Time", "time"))
# Pick the variance-like column you actually want:
# Put preferred column names FIRST.
var_val_col  <- pick_col(var, c("Variance", "var", "Var", "GeneVar", "gene_var", "ResVar", "ResLog2Var", "log2Var", "Log2Var"))

if (is.null(var_gene_col)) stop_cols(var, "Variance table missing Gene column.")
if (is.null(var_day_col))  stop_cols(var, "Variance table missing Day/Time column.")
if (is.null(var_val_col))  stop_cols(var, "Variance table missing a variance-like column (e.g., Variance, ResLog2Var, log2Var).")

# Standardize
var$Gene <- as.character(var[[var_gene_col]])
var$Day  <- suppressWarnings(as.numeric(var[[var_day_col]]))
var$VarianceValue <- suppressWarnings(as.numeric(var[[var_val_col]]))

if (anyNA(var$Day)) cat("\nWARNING: Some Day values are NA after numeric conversion.\n")
if (anyNA(var$VarianceValue)) cat("\nWARNING: Some VarianceValue values are NA after numeric conversion.\n")

cat("\nVariance mapping:\n  Gene column: ", var_gene_col,
    "\n  Day column: ", var_day_col,
    "\n  Variance column used: ", var_val_col, "\n", sep = "")

# Define pre/post sets (as specified)
pre  <- var[var$Day <= 12, ]
post <- var[var$Day >= 16, ]

# Aggregate mean variance value by gene in pre and post
pre_var  <- aggregate(VarianceValue ~ Gene, data = pre,  FUN = mean, na.rm = TRUE)
post_var <- aggregate(VarianceValue ~ Gene, data = post, FUN = mean, na.rm = TRUE)

var_ratio <- merge(pre_var, post_var, by = "Gene", all = TRUE)
names(var_ratio) <- c("Gene", "Var_Pre14", "Var_Post14")
var_ratio$VarRatio_Post14_vs_Pre14 <- var_ratio$Var_Post14 / var_ratio$Var_Pre14

# Keep only requested ratio column
var_ratio <- var_ratio[, c("Gene", "VarRatio_Post14_vs_Pre14")]

# -----------------------------
# 5) Functional category mapping
#    Gene2FunctiionMappiing4Edgeswap.csv columns: gene + pathway (or function_name)
# -----------------------------
cat("\nSelect Gene2FunctionMapping file (e.g., Gene2FunctiionMappiing4Edgeswap.csv)\n")
func_path <- file.choose()
func <- read.csv(func_path, check.names = FALSE)

func_gene_col <- pick_col(func, c("Gene", "gene", "GENE"))
func_cat_col  <- pick_col(func, c("pathway", "function_name", "Functional_Category", "Category", "Function"))
if (is.null(func_gene_col) || is.null(func_cat_col)) {
  stop_cols(func, "Function mapping must contain gene + category columns (e.g., gene, pathway).")
}

func_df <- unique(data.frame(
  Gene = as.character(func[[func_gene_col]]),
  Functional_Category = as.character(func[[func_cat_col]]),
  stringsAsFactors = FALSE
))

cat("\nFunction mapping:\n  Gene column: ", func_gene_col, "\n  Category column: ", func_cat_col, "\n", sep = "")

# -----------------------------
# 6) Merge everything into one gene table
# -----------------------------
gene_table <- Reduce(
  function(x, y) merge(x, y, by = "Gene", all = TRUE),
  list(dtw_df, pc1_df, bridge_df, var_ratio, func_df)
)

# Fill bridging missing as NO
gene_table$Ever_Bridging[is.na(gene_table$Ever_Bridging)] <- "NO"

# -----------------------------
# 7) Save outputs
# -----------------------------
out_dir <- dirname(dtw_path)
out_csv <- file.path(out_dir, "GeneRole_Integrated_Table.csv")
write.csv(gene_table, out_csv, row.names = FALSE)

cat("\nDONE.\nSaved:\n  ", out_csv, "\n", sep = "")

# Optional quick diagnostics
cat("\nQuick diagnostics:\n")
cat("  Rows (genes): ", nrow(gene_table), "\n", sep = "")
cat("  DTW_cluster missing: ", sum(is.na(gene_table$DTW_cluster)), "\n", sep = "")
cat("  Max_PC1_Window missing: ", sum(is.na(gene_table$Max_PC1_Window)), "\n", sep = "")
cat("  VarRatio missing: ", sum(is.na(gene_table$VarRatio_Post14_vs_Pre14)), "\n", sep = "")
cat("  Functional_Category missing: ", sum(is.na(gene_table$Functional_Category)), "\n", sep = "")
cat("  Ever_Bridging YES count: ", sum(gene_table$Ever_Bridging == "YES"), "\n", sep = "")