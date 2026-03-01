# ==============================================================================
# 📊 EXPRESSION PLOTTER WITH PARALLELIZED OUTPUT AND FDR-CORRECTED T-TESTS
# ==============================================================================
# DESCRIPTION:
# This R script generates expression plots (PNG, 300 dpi) for multiple genes across
# biological conditions (e.g., time, treatment), where:
#   - Columns = genes
#   - Rows = samples
#   - One column contains the grouping/category label (e.g., treatment group)

# The script supports:
#   ✅ Interactive file input and output directory selection
#   ✅ User-defined grouping column and ignored columns
#   ✅ Parallelized plot generation with a progress bar
#   ✅ Bar plots with p-values & FDR for 2 groups
#   ✅ Line or bar plots for >2 groups (user decides)
#   ✅ Choice of uncertainty display: Standard Error or 95% CI (shadow ribbon)
#   ✅ Full ANOVA + Tukey HSD results exported when >2 groups

# ==============================================================================
# 🧪 HOW TO USE THE SCRIPT
# ==============================================================================

# 1️⃣ SELECT INPUT FILE — CSV with samples in rows, genes in columns
# 2️⃣ SELECT OUTPUT DIRECTORY — a timestamped subfolder will be created
# 3️⃣ CHOOSE GROUPING COLUMN — used to define biological replicates
# 4️⃣ IGNORE COLUMNS — e.g., SampleID or notes columns
# 5️⃣ CHOOSE PLOT TYPE — bar or line plot (for >2 groups)
# 6️⃣ CHOOSE UNCERTAINTY TYPE — SE bars or 95% CI (shadow ribbon for lines)
# 7️⃣ CHOOSE NUMBER OF CORES — for fast parallel processing

# OUTPUT:
#   📁 Individual gene plots (PNG, 300 dpi)
#   📄 (If >2 groups): CSV with ANOVA and Tukey HSD results
# ANOVA and Tukey HSD Results Explanation:
# - diff: Mean difference between group pairs (from Tukey HSD)
# - lwr, upr: Lower and upper bounds of the 95% confidence interval (from Tukey HSD)
# - padj: Adjusted p-value for multiple comparisons (Tukey HSD)
# - Comparison: Group pair being compared (e.g., "14-12")
# - Gene: Name of the gene tested
# - FDR: False Discovery Rate (may be additional correction across all genes)
# Note: ANOVA (aov) is used to detect overall significance; Tukey is used post-ANOVA
#       to identify which specific group pairs differ.

# ==============================================================================
# 🛠️ DEPENDENCIES:
# install.packages(c("ggplot2", "tidyr", "readr", "tools",
#                    "future", "furrr", "progressr", "multcomp"))
# ==============================================================================

# ---------------------- Load Required Libraries ----------------------
library(ggplot2)
library(tidyr)
library(readr)
library(tools)
library(future)
library(furrr)
library(progressr)
library(multcomp)

# ---------------------- Script Metadata ----------------------
script_info <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    this_path <- rstudioapi::getSourceEditorContext()$path
  } else {
    this_path <- sys.frame(1)$ofile
  }
  if (is.null(this_path) || this_path == "") {
    list(name = "Unnamed_Script", path = getwd())
  } else {
    list(
      name = tools::file_path_sans_ext(basename(this_path)),
      path = normalizePath(dirname(this_path))
    )
  }
}, error = function(e) list(name = "Unnamed_Script", path = getwd()))

script_name <- script_info$name
script_path <- script_info$path

# ---------------------- File Selection ----------------------
if (interactive()) {
  cat("\n📂 Please select your input CSV file:\n")
  file_path <- file.choose()
} else {
  stop("❌ This script must be run in interactive mode to select files.")
}

cat("\n📂 Select any file in your desired output directory:\n")
base_dir <- dirname(file.choose())

input_name <- file_path_sans_ext(basename(file_path))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_dir, paste0(input_name, "_", timestamp))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("✅ Output will be saved in:", output_dir, "\n")

# ---------------------- Read Data ----------------------
data <- read_csv(file_path, show_col_types = FALSE)

# ---------------------- Ask Which Column Is Grouping Column ----------------------
cat("\n📋 First 10 column headers:\n")
header_preview <- head(colnames(data), 10)
for (i in seq_along(header_preview)) {
  cat(paste0(i, ": ", header_preview[i], "\n"))
}

repeat {
  group_col_index <- as.integer(readline(prompt = "🔍 Enter the number of the grouping column: "))
  if (!is.na(group_col_index) &&
      group_col_index >= 1 &&
      group_col_index <= length(colnames(data))) break
  cat("⚠️ Invalid selection. Try again.\n")
}

group_col <- colnames(data)[group_col_index]

# ---------------------- Ask Which Columns to Ignore ----------------------
cat("\n❓ Would you like to ignore any columns (e.g., sample ID)?\n")
for (i in seq_along(colnames(data))) {
  cat(paste0(i, ": ", colnames(data)[i], "\n"))
}
ignore_input  <- readline(prompt = "📛 Enter comma-separated numbers of columns to ignore (or press Enter to skip): ")
ignore_indices <- as.integer(strsplit(ignore_input, ",")[[1]])
ignore_indices <- ignore_indices[!is.na(ignore_indices)]
ignore_cols   <- unique(c(group_col, colnames(data)[ignore_indices]))
gene_cols     <- setdiff(colnames(data), ignore_cols)

# ---------------------- Reshape Data ----------------------
long_df <- tidyr::pivot_longer(
  data,
  cols      = all_of(gene_cols),
  names_to  = "Gene",
  values_to = "Expression"
)

colnames(long_df)[which(colnames(long_df) == group_col)] <- "Group"
long_df$Group <- as.character(long_df$Group)

# ---------------------- Mean, SE, and 95% CI (Base R) ----------------------
mean_se_list <- split(long_df, list(long_df$Gene, long_df$Group), drop = TRUE)

summary_rows <- lapply(mean_se_list, function(sub) {
  gene  <- unique(sub$Gene)
  group <- unique(sub$Group)
  vals  <- na.omit(sub$Expression)
  n     <- length(vals)
  
  if (n == 0) {
    mean_val <- NA
    se_val   <- NA
    ci_lower <- NA
    ci_upper <- NA
  } else if (n == 1) {
    mean_val <- vals[1]
    se_val   <- NA
    ci_lower <- NA
    ci_upper <- NA
  } else {
    mean_val <- mean(vals)
    se_val   <- sd(vals) / sqrt(n)
    t_crit   <- qt(0.975, df = n - 1)
    ci_lower <- mean_val - t_crit * se_val
    ci_upper <- mean_val + t_crit * se_val
  }
  
  data.frame(
    Gene     = gene,
    Group    = group,
    N        = n,
    Mean     = mean_val,
    SE       = se_val,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)

# ---------------------- Ask User How Many Cores to Use ----------------------
total_cores <- future::availableCores()
cat(paste0("\n🧠 You have ", total_cores, " available logical cores.\n"))

max_cores <- max(1, total_cores - 2)
repeat {
  user_cores <- as.integer(readline(prompt = paste0("🔧 How many cores would you like to use? (1–", max_cores, "): ")))
  if (!is.na(user_cores) && user_cores >= 1 && user_cores <= max_cores) break
  cat("⚠️  Invalid input. Please enter a number between 1 and ", max_cores, ".\n")
}

# ---------------------- Set Up Parallel Plan ----------------------
plan(multisession, workers = user_cores)
handlers(global = TRUE)

# ---------------------- Get Group Levels and Ask Plot Type ----------------------
group_levels_raw <- sort(unique(summary_df$Group), na.last = NA)

suppressWarnings({
  numeric_possible <- !any(is.na(as.numeric(group_levels_raw)))
})

if (numeric_possible) {
  group_levels <- sort(as.numeric(group_levels_raw))
  group_levels <- as.character(group_levels)
} else {
  group_levels <- sort(group_levels_raw)
}

plot_type <- "line"

if (length(group_levels) > 2) {
  repeat {
    user_plot <- readline(prompt = "📊 You have more than 2 groups. Choose plot type: [1] Bar [2] Line: ")
    if (user_plot %in% c("1", "2")) break
    cat("⚠️ Please enter 1 or 2.\n")
  }
  plot_type <- if (user_plot == "1") "bar" else "line"
} else {
  plot_type <- "bar"
}

# ---------------------- Ask How to Display Uncertainty ----------------------
cat("\n📈 How should uncertainty be displayed on plots?\n")
cat("  [1] Standard error bars (Mean ± SE)\n")
cat("  [2] 95% confidence interval (shadow ribbon for line plots; CI error bars for bars)\n")

repeat {
  user_interval <- readline(prompt = "👉 Enter 1 for SE or 2 for 95% CI: ")
  if (user_interval %in% c("1", "2")) break
  cat("⚠️ Please enter 1 or 2.\n")
}

interval_type <- if (user_interval == "1") "se" else "ci"

# ---------------------- Precompute t-tests or ANOVA ----------------------
# ---------------------- Precompute t-tests or ANOVA ----------------------
if (length(group_levels) == 2) {
  # ----- Two-group comparison -----
  ttest_list <- list()
  group1_name <- group_levels[1]
  group2_name <- group_levels[2]
  
  for (g in unique(long_df$Gene)) {
    raw_data <- long_df[long_df$Gene == g, ]
    group1_vals <- raw_data$Expression[raw_data$Group == group1_name]
    group2_vals <- raw_data$Expression[raw_data$Group == group2_name]
    
    # 🚿 Remove NA values BEFORE checking lengths
    group1_vals <- na.omit(group1_vals)
    group2_vals <- na.omit(group2_vals)
    
    if (length(group1_vals) > 1 && length(group2_vals) > 1) {
      # Running t-test: Group2 vs Group1
      ttest <- tryCatch(t.test(group2_vals, group1_vals), error = function(e) NULL)
      
      if (!is.null(ttest)) {
        # Fetch means from the pre-computed summary table
        mean1 <- summary_df$Mean[summary_df$Gene == g & summary_df$Group == group1_name]
        mean2 <- summary_df$Mean[summary_df$Gene == g & summary_df$Group == group2_name]
        
        ttest_list[[g]] <- data.frame(
          Gene = g,
          Comparison = paste0(group2_name, " vs. ", group1_name),
          Mean_Group1 = mean1,
          Mean_Group2 = mean2,
          Mean_Difference = ttest$estimate[1] - ttest$estimate[2], # Mean2 - Mean1
          Raw_P_Value = ttest$p.value,
          T_Statistic = ttest$statistic,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  ttest_df <- do.call(rbind, ttest_list)
  ttest_df$Adjusted_P_Value <- p.adjust(ttest_df$Raw_P_Value, method = "fdr")
  
  # New file output for 2-group results
  write_csv(ttest_df, file.path(output_dir, "ttest_fdr_results.csv"))
  cat("📄 T-test + FDR results saved to ttest_fdr_results.csv\n")
  
  # Re-define fdr_vals as a named vector for the plotting function later (CRITICAL)
  fdr_vals <- setNames(ttest_df$Adjusted_P_Value, ttest_df$Gene)
  
} else {
  # ----- Multi-group ANOVA + Tukey -----
  # ... (Keep existing multi-group code here)
}
  # ----- Multi-group ANOVA + Tukey -----
  tukey_list <- list()
  
  for (g in unique(long_df$Gene)) {
    raw_data <- long_df[long_df$Gene == g, ]
    raw_data <- na.omit(raw_data[, c("Expression", "Group")])
    
    model <- tryCatch(aov(Expression ~ Group, data = raw_data), error = function(e) NULL)
    
    if (!is.null(model)) {
      tukey <- tryCatch(TukeyHSD(model), error = function(e) NULL)
      if (!is.null(tukey)) {
        tdf <- as.data.frame(tukey$Group)
        tdf$Comparison <- rownames(tdf)
        tdf$Gene       <- g
        tukey_list[[g]] <- tdf
      }
    }
  }
  
  tukey_df      <- do.call(rbind, tukey_list)
  tukey_df$FDR  <- p.adjust(tukey_df$`p adj`, method = "fdr")
  write_csv(tukey_df, file.path(output_dir, "anova_tukey_results.csv"))
  cat("📄 ANOVA + Tukey HSD results saved to anova_tukey_results.csv\n")


# ---------------------- Define Plotting Function ----------------------
plot_gene <- function(g) {
  plot_data <- summary_df[summary_df$Gene == g, ]
  raw_data  <- long_df[long_df$Gene == g, ]
  
  plot_data$Group <- factor(plot_data$Group, levels = group_levels)
  
  # Precompute y-range for error bars if needed
  if (interval_type == "se") {
    plot_data$ymin <- plot_data$Mean - plot_data$SE
    plot_data$ymax <- plot_data$Mean + plot_data$SE
  } else {
    plot_data$ymin <- plot_data$CI_lower
    plot_data$ymax <- plot_data$CI_upper
  }
  
  if (length(group_levels) == 2) {
    # ----- Two-group: bar + stats in subtitle -----
    group1_vals <- raw_data$Expression[raw_data$Group == group_levels[1]]
    group2_vals <- raw_data$Expression[raw_data$Group == group_levels[2]]
    
    p_val  <- tryCatch(t.test(group1_vals, group2_vals)$p.value, error = function(e) NA)
    fdr_val <- fdr_vals[g]
    
    label_text <- paste0("p = ", signif(p_val, 3),
                         if (!is.na(fdr_val)) paste0(", FDR = ", signif(fdr_val, 3)) else "")
    
    p <- ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "black") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3) +
      labs(
        title    = paste("Expression of", g),
        subtitle = label_text,
        x        = group_col,
        y        = if (interval_type == "se") "Mean Expression ± SE" else "Mean Expression (95% CI)"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
    
  } else if (plot_type == "bar") {
    # ----- Multi-group: bar plot -----
    p <- ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "black") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3) +
      labs(
        title = paste("Expression of", g),
        x     = group_col,
        y     = if (interval_type == "se") "Mean Expression ± SE" else "Mean Expression (95% CI)"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
    
  } else {
    # ----- Multi-group: line plot -----
    if (interval_type == "se") {
      # Line + SE error bars
      p <- ggplot(plot_data, aes(x = Group, y = Mean, group = 1)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3) +
        labs(
          title = paste("Expression of", g),
          x     = group_col,
          y     = "Mean Expression ± SE"
        ) +
        theme_minimal(base_size = 14)
    } else {
      # Line + shadowed 95% CI ribbon
      plot_data$x_num <- as.numeric(plot_data$Group)
      
      p <- ggplot(plot_data, aes(x = x_num, y = Mean, group = 1)) +
        geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        scale_x_continuous(
          breaks = seq_along(levels(plot_data$Group)),
          labels = levels(plot_data$Group)
        ) +
        labs(
          title = paste("Expression of", g),
          x     = group_col,
          y     = "Mean Expression (95% CI)"
        ) +
        theme_minimal(base_size = 14)
    }
  }
  
  ggsave(
    filename = file.path(output_dir, paste0(g, "_expression.png")),
    plot     = p,
    dpi      = 600,
    width    = 8,
    height   = 6,
    units    = "in",
    type     = "cairo"
  )
}

# ---------------------- Run Parallel Plotting with Progress Bar ----------------------
genes <- unique(summary_df$Gene)
cat("\n⏳ Generating plots in parallel...\n")

with_progress({
  future_map(genes, plot_gene, .progress = TRUE)
})

cat("✅ All plots saved to:", output_dir, "\n")

# ---------------------- Auto-generate .qmd Report ----------------------
interval_description <- if (interval_type == "se") {
  "Standard error bars (Mean ± SE)"
} else {
  "95% confidence intervals (shadow ribbon for line plots; CI error bars for bar plots)"
}

qmd_lines <- c(
  "---",
  paste0("title: \"", script_name, " Report\""),
  "author: \"Auto-generated\"",
  paste0("date: \"", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\""),
  "format: html",
  "---",
  "",
  "# 📊 Script Summary",
  "",
  "This Quarto report documents the results and configuration of the R script:",
  "",
  paste0("- **Script name**: `", script_name, ".R`"),
  paste0("- **Script path**: `", script_path, "`"),
  paste0("- **Input file**: `", normalizePath(file_path), "`"),
  paste0("- **Output directory**: `", normalizePath(output_dir), "`"),
  paste0("- **Uncertainty display**: ", interval_description),
  paste0("- **Generated on**: `", format(Sys.time(), '%Y-%m-%d %H:%M:%S'), "`"),
  "",
  "# 📁 What This Script Does",
  "",
  "This script generates publication-quality expression plots (PNG format, 300 dpi) for multiple genes across different biological conditions.",
  "",
  "## Workflow Overview",
  "",
  "- Takes a CSV input file with samples as rows and genes as columns.",
  "- Allows the user to select a grouping column and columns to ignore.",
  "- Uses parallel processing to generate plots quickly.",
  "- Applies:",
  "  - **t-tests** and FDR correction if 2 groups are present.",
  "  - **ANOVA and Tukey HSD** if >2 groups are detected.",
  "- Lets the user choose whether to display uncertainty as **standard error bars** or **95% confidence intervals**.",
  "- Outputs include plots and statistical test results.",
  "",
  "# 📦 Dependencies",
  "",
  "The following R packages are required:",
  "",
  "- `ggplot2`, `tidyr`, `readr`, `tools`, `future`, `furrr`, `progressr`, `multcomp`",
  "",
  "# 🔍 Statistical Output Documentation",
  "",
  "ANOVA and Tukey HSD Results Explanation:",
  "",
  "- `diff`: Mean difference between group pairs (from Tukey HSD)",
  "- `lwr`, `upr`: Lower and upper bounds of the 95% confidence interval (from Tukey HSD)",
  "- `padj`: Adjusted p-value for multiple comparisons (Tukey HSD)",
  "- `Comparison`: Group pair being compared (e.g., \"14-12\")",
  "- `Gene`: Name of the gene tested",
  "- `FDR`: False Discovery Rate (may be additional correction across all genes)",
  "- Note: ANOVA (aov) is used to detect overall significance; Tukey is used post-ANOVA to identify which specific group pairs differ.",
  "",
  "# 🗂️ Output Files",
  "",
  "- PNG plots: One per gene, saved at 300 dpi.",
  if (length(group_levels) > 2) {
    "- `anova_tukey_results.csv`: Contains full ANOVA and Tukey HSD test results with FDR correction."
  } else {
    "- No ANOVA output file (only t-tests used across genes)."
  },
  "",
  "# ✅ Run Complete",
  "",
  "This report was generated automatically upon completion of the script."
)

report_path <- file.path(output_dir, paste0(script_name, "_Report.qmd"))
writeLines(qmd_lines, con = report_path)
cat("📝 Quarto report written to:", report_path, "\n")

# ---------------------- Render the .qmd Report Automatically ----------------------
if (requireNamespace("quarto", quietly = TRUE)) {
  tryCatch({
    quarto::quarto_render(report_path)
    cat("✅ Quarto report successfully rendered.\n")
  }, error = function(e) {
    cat("⚠️ Rendering failed: ", e$message, "\n")
  })
} else {
  cat("⚠️ 'quarto' package is not installed. Skipping rendering.\n")
}
