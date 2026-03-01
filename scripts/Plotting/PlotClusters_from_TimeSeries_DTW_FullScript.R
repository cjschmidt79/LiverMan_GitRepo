
# ============================================
# 📜 Script: CI_Cluster_Trajectories_Analysis_Final_Fixed.R
# Description: Handles biological replicates with proxy names and computes 95% CI ribbons.
# ============================================

# -------------------------------
# Load Required Libraries
# -------------------------------
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("data.table")) install.packages("data.table")
if (!require("fs")) install.packages("fs")

library(tidyverse)
library(data.table)
library(fs)

# -------------------------------
# 📁 File Selection (Interactive Mac-Friendly)
# -------------------------------
if (interactive()) {
  cat("\nPlease select your replicate-level expression CSV file (e.g. Liver5SD.csv):\n")
  input_expr <- file.choose()

  cat("\nPlease select your cluster assignment CSV file (e.g. DTWclust_Optimal_Cluster_Assignments.csv):\n")
  input_cluster <- file.choose()

  cat("\nPlease choose your base output directory:\n")
  base_output <- dirname(file.choose())
}

# -------------------------------
# 📂 Create Output Directory
# -------------------------------
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_dir <- file.path(base_output, paste0("CI_Trajectory_Output_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 🧪 Load and Process Expression Data
# -------------------------------
expr_df <- fread(input_expr)
setnames(expr_df, 1, "Day")
expr_df[, (2:ncol(expr_df)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(expr_df)]

# Transpose expression values
expr_data <- expr_df[, -"Day"]
transposed_expr <- as.data.table(t(expr_data))

# Generate proxy column names after transpose to preserve order
proxy_ids_fixed <- paste0("D", rep(expr_df$Day, each = 1), "_rep", sequence(rle(expr_df$Day)$lengths))
setnames(transposed_expr, proxy_ids_fixed)

# Coerce to numeric
transposed_expr[, (proxy_ids_fixed) := lapply(.SD, as.numeric), .SDcols = proxy_ids_fixed]

# Add gene names
transposed_expr[, Gene := colnames(expr_data)]
setcolorder(transposed_expr, c("Gene", setdiff(names(transposed_expr), "Gene")))

# Generate and save metadata
proxy_metadata <- data.table(Proxy = proxy_ids_fixed, Day = expr_df$Day)
fwrite(proxy_metadata, file.path(output_dir, "Proxy_Metadata.csv"))

# -------------------------------
# 🧩 Load Cluster Assignments
# -------------------------------
cluster_df <- fread(input_cluster)
if (!"Gene" %in% colnames(cluster_df)) {
  setnames(cluster_df, 1, "Gene")
}

# Merge and reshape
merged_df <- merge(transposed_expr, cluster_df, by = "Gene")
# Ensure all columns except Gene and Cluster are numeric
expr_cols <- setdiff(names(merged_df), c("Gene", "Cluster"))
merged_df[, (expr_cols) := lapply(.SD, as.numeric), .SDcols = expr_cols]

long_df <- melt(merged_df, id.vars = c("Gene", "Cluster"),
                variable.name = "Proxy", value.name = "Expression")

# Map back to real day using metadata
proxy_metadata <- unique(proxy_metadata)
long_df <- merge(long_df, proxy_metadata, by = "Proxy")
long_df[, Day := as.numeric(as.character(Day))]
long_df[, Expression := as.numeric(Expression)]

# -------------------------------
# 🧮 Compute 95% CI via Bootstrapping
# -------------------------------
set.seed(123)
boot_ci <- function(x, reps = 1000, conf = 0.95) {
  bt <- replicate(reps, mean(sample(x, replace = TRUE), na.rm = TRUE))
  ci <- quantile(bt, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE)
  return(c(mean = mean(x, na.rm = TRUE), lower = ci[1], upper = ci[2]))
}

agg_ci <- long_df[, {
  ci <- boot_ci(Expression)
  .(mean = ci[1], lower = ci[2], upper = ci[3])
}, by = .(Cluster, Day)]

# -------------------------------
# 📊 Plot with Shaded CI
# -------------------------------
plot_path <- file.path(output_dir, "Cluster_CI_Trajectories.png")
p <- ggplot(agg_ci, aes(x = Day, y = mean)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
  facet_wrap(~ Cluster, scales = "free_y") +
  labs(title = "Cluster Trajectories with 95% Confidence Intervals",
       x = "Day", y = "Mean Expression") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggsave(plot_path, plot = p, width = 10, height = 6, dpi = 300)

# -------------------------------
# 📝 Generate QMD Summary
# -------------------------------
qmd_path <- file.path(output_dir, "Cluster_Trajectory_Summary.qmd")
script_name <- tryCatch(basename(sys.frame(1)$ofile), error = function(e) "CI_Cluster_Trajectories_Analysis_Final_Fixed.R")

tree_output <- paste(capture.output(fs::dir_tree(output_dir, recurse = 2)), collapse = "\n")
input_expr_tree <- paste(capture.output(fs::dir_tree(dirname(input_expr), recurse = 1)), collapse = "\n")
input_cluster_tree <- paste(capture.output(fs::dir_tree(dirname(input_cluster), recurse = 1)), collapse = "\n")

writeLines(c(
  "---",
  "title: \"Cluster Trajectories with 95% Confidence Intervals\"",
  "format: html",
  "---",
  "",
  paste0("**Run Timestamp:** `", timestamp, "`"),
  paste0("**Script Name:** `", script_name, "`"),
  "",
  "### 📂 Input Files",
  "",
  "```plaintext",
  "Expression File:",
  input_expr_tree,
  "",
  "Cluster Assignment File:",
  input_cluster_tree,
  "```",
  "",
  "### 📁 Output Directory Structure",
  "```plaintext",
  tree_output,
  "```",
  "",
  "### 🔍 Biological Significance",
  "",
  "- Expression replicates were assigned proxy names to preserve day grouping.",
  "- 95% CIs were calculated per cluster per day based on replicates.",
  "- This highlights both trends and the certainty of observed patterns.",
  "",
  "### 📊 Figure Preview",
  "",
  paste0("![](", basename(plot_path), ")")
), qmd_path)

cat("\n✅ All outputs saved in:", output_dir, "\n")

# -------------------------------
# 📄 Auto-generate and Render Quarto Report
# -------------------------------

if (!require("quarto")) install.packages("quarto")
if (!require("knitr")) install.packages("knitr")

# Required variables assumed from earlier in the script
# - output_dir (already created)
# - plot_path (main plot file)
# - input_expr, input_cluster (input files)

qmd_path <- file.path(output_dir, "Cluster_Trajectory_Summary.qmd")
html_path <- file.path(output_dir, "Cluster_Trajectory_Summary.html")

# Extract script name and location
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "CI_Cluster_Trajectories_Analysis_Final_Fixed.R")
script_dir <- dirname(script_path)
script_name <- basename(script_path)

# Get file listings
output_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
input_files <- c(input_expr, input_cluster)

# Create QMD report content
report_lines <- c(
  "---",
  "title: \"Cluster Trajectories with 95% Confidence Intervals\"",
  "format: html",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "# 📊 Summary Report",
  "",
  "**Script:**", paste0("`", script_name, "`"),
  "",
  "**Directory:**", paste0("`", script_dir, "`"),
  "",
  "**Generated On:**", paste0("`", format(Sys.time(), '%Y-%m-%d %H:%M:%S'), "`"),
  "",
  "## 🧠 Script Overview",
  "",
  "- This R script analyzes replicate-level gene expression data across timepoints.",
  "- Cluster assignments are provided separately, and each gene is matched to its cluster.",
  "- For each cluster and time point, the script computes a **bootstrapped 95% confidence interval** of expression levels.",
  "- The final output is a faceted line plot per cluster with shaded CI ribbons.",
  "- Output files are stored in a timestamped folder for reproducibility.",
  "",
  "## 📥 Inputs",
  "",
  "**Expression File:**",
  paste0("- `", normalizePath(input_expr), "`"),
  "",
  "**Cluster Assignment File:**",
  paste0("- `", normalizePath(input_cluster), "`"),
  "",
  "## 📤 Outputs",
  "",
  "**Output Directory:**",
  paste0("`", normalizePath(output_dir), "`"),
  "",
  "**Files Generated:**",
  paste0("```plaintext\n", paste(basename(output_files), collapse = "\n"), "\n```"),
  "",
  "## 📦 Dependencies",
  "",
  "- tidyverse",
  "- data.table",
  "- ggplot2",
  "- fs",
  "- knitr",
  "- quarto",
  "",
  "## 🧮 Statistical Logic",
  "",
  "We calculate the mean and 95% confidence interval for each cluster's expression values per day using **non-parametric bootstrapping**:",
  "",
  "$$ CI = [Q_{\\alpha/2}, Q_{1-\\alpha/2}] \\text{ from resampled means} $$",
  "",
  "Where $Q$ denotes the quantile from bootstrapped resamples.",
  "",
  "The code uses the following function:",
  "```r",
  "boot_ci <- function(x, reps = 1000, conf = 0.95) {",
  "  x <- x[!is.na(x)]",
  "  if (length(x) == 0) return(c(mean = NA, lower = NA, upper = NA))",
  "  bt <- replicate(reps, mean(sample(x, replace = TRUE)))",
  "  ci <- quantile(bt, probs = c((1 - conf)/2, 1 - (1 - conf)/2))",
  "  return(c(mean = mean(x), lower = ci[1], upper = ci[2]))",
  "}",
  "```",
  "",
  "## 📈 Cluster Trajectory Plot",
  "",
  "```{r echo=FALSE, out.width='100%'}",
  paste0("knitr::include_graphics(\"", normalizePath(plot_path), "\")"),
  "```",
  "",
  "## 🧾 Session Info",
  "",
  "```{r}",
  "sessionInfo()",
  "```"
)

# Write QMD file
writeLines(report_lines, qmd_path)
cat("\n📝 Quarto .qmd report written to:\n", qmd_path, "\n")

# Render HTML
cat("\n🧪 Rendering HTML report...\n")
quarto::quarto_render(qmd_path)

# Confirm output
cat("\n✅ HTML report saved to:\n", html_path, "\n")

# Open in browser if interactive
if (interactive()) browseURL(html_path)
