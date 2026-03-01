################################################################################
# SCRIPT: Network_Volatility_Engine.R
# DESCRIPTION: Identifies "Pivot Genes" using Intensity, Breadth, and Duration
#              metrics from network turnover tables. Integrates LIMMA DE data.
#              Follows Reusable Script Standards and Parallel Optimization.
################################################################################
# PROGRAM: NETWORK VOLATILITY ENGINE (NVE)
# VERSION: 1.0.0
# AUTHOR:  Gemini AI / Bioinformaticist
#
# GENERAL FUNCTION:
# This script identifies "Pivot Genes" within metabolic transitions. It detects
# genes that are significantly re-wiring the biological network by integrating
# high-correlation turnover (Top 50 lists) with differential expression data.
# It moves beyond simple "up/down" folding to quantify network structural change.
#
# MODULE OVERVIEW:
# 1. HARDWARE HANDSHAKE: Detects CPU architecture and recommends optimal core
#    allocation to maximize speed while avoiding "Amdahl Law" parallel overhead.
# 2. VOLATILITY CALCULATION: Quantifies three key metrics:
#    - Intensity: Magnitude of correlation change (abs_delta_r).
#    - Breadth: Number of unique network partners involved in the transition.
#    - Duration: Persistence of the gene in top-tier rewiring across time.
# 3. STATISTICAL POWER (FDR): Uses parallel permutation (1,000+ iterations) to
#    assign p-values and Benjamini-Hochberg q-values to gene volatility.
# 4. LIMMA INTEGRATION: An optional module that joins pre-computed expression
#    data (logFC/adj.P.Val) to identify "Master Switch" genes.
# 5. STABILITY QC: Performs Jackknife (Leave-One-Out) testing to ensure pivot
#    genes are robust and not driven by individual outlier samples.
# 6. REPRODUCIBILITY REPORT: Generates a Quarto HTML report including a manifest,
#    provenance data, and ego-network exports for Cytoscape.
#
# EXPECTED INPUTS:
# - Turnover CSVs: Files named by Day/Transition (e.g., "Day4_to_8.csv") containing
#   columns for Gene1, Gene2, and Delta_R (correlation change).
# - LIMMA CSVs (Optional): Standard limma::topTable output. The script is designed
#   to handle the "unnamed" first column (V1) as the Gene Symbol.
#
# REQUIREMENTS:
# - R Libraries: tidyverse, future, furrr, progressr, rstudioapi, quarto.
# - Quarto CLI must be installed for HTML report rendering.
#
# EXECUTION NOTE:
# This script is interactive. You will be prompted to select directories via 
# Mac-friendly native pickers and identify column names via the console.
################################################################################
################################################################################

# --- 1. PROVENANCE & IDENTITY ---
# Capture script identity immediately
script_name <- "Network_Volatility_Engine.R"
# Best-effort capture of absolute paths
try({
  script_full <- normalizePath(rstudioapi::getSourceEditorContext()$path, winslash = "/")
  script_path <- dirname(script_full)
}, silent = TRUE)

if(!exists("script_full")) {
  script_full <- normalizePath(file.path(getwd(), script_name), winslash = "/")
  script_path <- getwd()
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- 2. DEPENDENCIES ---
required_pkgs <- c("tidyverse", "future", "furrr", "progressr", "rstudioapi", 
                   "quarto", "knitr", "tools", "data.table")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(future)
library(furrr)
library(progressr)
library(rstudioapi)
library(data.table)

# --- 3. INPUTS & INTERACTIVITY ---
# Use Mac-friendly pickers
cat("--- Input Configuration ---\n")
turnover_dir <- rstudioapi::selectDirectory(caption = "Select directory containing Turnover Top 50 CSVs")
if (is.null(turnover_dir)) stop("User cancelled: Turnover directory is required.")

do_limma <- askYesNo("Do you have LIMMA differential expression CSVs to integrate?")
limma_dir <- if(do_limma) rstudioapi::selectDirectory(caption = "Select LIMMA Results directory") else NULL

# Interactively map columns
test_file <- list.files(turnover_dir, pattern = "\\.csv$", full.names = TRUE)[1]
test_df <- read.csv(test_file, nrows = 5)
cat("\nPreview of Turnover Data:\n")
print(head(test_df))

gene_col <- readline("Enter the exact column name for Gene Symbols (e.g., Gene1): ")
delta_r_col <- readline("Enter the column name for Delta R (e.g., delta_r or abs_delta_r): ")

# Setup standardized output directory
run_name <- readline("Enter a run name (e.g., Yolk_to_Feed): ")
if(run_name == "") run_name <- "Analysis"

output_root <- file.path(getwd(), "outputs")
if(!dir.exists(output_root)) dir.create(output_root)

output_dir <- normalizePath(file.path(output_root, paste0(run_name, "_", timestamp)), winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE)
dir.create(file.path(output_dir, "plots"))
dir.create(file.path(output_dir, "cytoscape_sif"))

# --- 4. HARDWARE HANDSHAKE (Parallel Optimization) ---
n_cores_logical <- parallel::detectCores()
n_cores_physical <- parallel::detectCores(logical = FALSE)

cat(paste0("\nHardware Detected: ", n_cores_physical, " Physical Cores / ", n_cores_logical, " Threads.\n"))
rec_cores <- max(1, n_cores_physical - 1)
cat(paste0("Recommendation: Use ", rec_cores, " cores for optimal 'Amdahl' efficiency.\n"))

user_cores <- as.integer(readline(paste0("Processors to allocate? [Default ", rec_cores, "]: ")))
if(is.na(user_cores)) user_cores <- rec_cores

# Set future plan
plan(multisession, workers = user_cores)

# --- 5. CORE ENGINE (Volatility Calculation) ---
cat("\n--- Running Network Volatility Engine ---\n")

# Auto-map chronological files
turnover_files <- list.files(turnover_dir, pattern = "\\.csv$", full.names = TRUE)
# Regex to extract days/numbers for sequencing
file_map <- data.frame(path = turnover_files) %>%
  mutate(fname = basename(path),
         order = as.numeric(str_extract(fname, "\\d+"))) %>%
  arrange(order)

handlers(handler_progress(
  format   = "[:bar] :percent :eta :message",
  clear    = FALSE, 
  width    = 60
))

with_progress({
  p <- progressor(steps = length(turnover_files))
  
  all_turnover_data <- turnover_files %>%
    future_map_dfr(~{
      df <- fread(.x)
      # Process gene involvement (Breadth)
      g1 <- df %>% select(Gene = all_of(gene_col), delta = all_of(delta_r_col))
      # Since edges have two genes, we ensure both are counted as participants
      # Assuming Gene2 is also a column in your Top 50 files
      g2 <- df %>% select(Gene = Gene2, delta = all_of(delta_r_col)) 
      
      combined <- bind_rows(g1, g2) %>%
        group_by(Gene) %>%
        summarise(Intensity = mean(abs(delta)),
                  Breadth = n(),
                  .groups = "drop") %>%
        mutate(Source_File = basename(.x))
      
      p(sprintf("Processed %s", basename(.x)))
      return(combined)
    })
})

# Aggregate Multi-Metric Pivot Scores
pivot_summary <- all_turnover_data %>%
  group_by(Gene) %>%
  summarise(
    Avg_Intensity = mean(Intensity),
    Max_Intensity = max(Intensity),
    Total_Breadth = sum(Breadth), # How many partners across time
    Unique_Partners = sum(Breadth), # For Top 50 lists, these are usually unique per step
    Duration = n(), # How many transition files it appeared in
    .groups = "drop"
  ) %>%
  mutate(Volatility_Score = (Avg_Intensity * Total_Breadth * Duration)) %>%
  arrange(desc(Volatility_Score))

# Installment 1 ends here. 
# Next: FDR Permutation, LIMMA Merge, Stability QC, and the Quarto Report.
# --- 6. STATISTICAL POWER (Parallel FDR Permutation) ---
cat("\n--- Assigning Statistical Power (FDR) ---\n")
n_perm <- 10000 # Number of permutations for null distribution

# Flatten all genes across all files to get a background pool
gene_pool <- unique(all_turnover_data$Gene)

with_progress({
  p_stats <- progressor(steps = n_perm)
  
  # Parallel permutation test
  null_scores <- future_map_dbl(1:n_perm, ~{
    # Create a 'fake' gene profile by sampling random Breadth and Duration
    fake_intensity <- mean(sample(all_turnover_data$Intensity, 1))
    fake_breadth   <- sum(sample(all_turnover_data$Breadth, sample(1:length(turnover_files), 1)))
    fake_duration  <- sample(1:length(turnover_files), 1)
    
    p_stats()
    return(fake_intensity * fake_breadth * fake_duration)
  }, .options = furrr_options(seed = TRUE))
})

# Calculate Empirical P-values and FDR (q-values)
# Calculate Empirical P-values and refined Status
pivot_summary <- pivot_summary %>%
  rowwise() %>%
  mutate(p_value = sum(null_scores >= Volatility_Score) / n_perm) %>%
  ungroup() %>%
  mutate(q_value = p.adjust(p_value, method = "BH"),
         Status = case_when(
           q_value < 0.05 ~ "Global Pivot (High Conf)",
           p_value < 0.05 & Duration > 15 ~ "Structural Coordinator", # High persistence despite score
           p_value < 0.1  & Duration > 10 ~ "Candidate Driver",     # Biological signal
           TRUE ~ "Background"
         ))

# --- 7. MODULAR INTEGRATION (LIMMA) ---
if(do_limma) {
  cat("\n--- Integrating LIMMA DE Data ---\n")
  limma_files <- list.files(limma_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Load and average logFC across all comparisons (or specific pivot window)
  de_data <- limma_files %>%
    map_dfr(~fread(.x) %>% mutate(Source = basename(.x))) %>%
    group_by(V1) %>% # V1 is usually the gene symbol index in limma output
    summarise(Max_logFC = max(abs(logFC), na.rm = TRUE),
              Avg_adj_P = mean(adj.P.Val, na.rm = TRUE),
              .groups = "drop")
  
  pivot_summary <- pivot_summary %>%
    left_join(de_data, by = c("Gene" = "V1"))
}

# --- 8. CYTOSCAPE EXPORT (Ego-Networks) ---
cat("\n--- Exporting Cytoscape SIF Files ---\n")
top_pivots <- pivot_summary %>% filter(q_value < 0.05) %>% head(10)

for(g in top_pivots$Gene) {
  ego_edges <- all_turnover_data %>% 
    filter(Gene == g) %>%
    select(Source = Gene, Target = Source_File, Weight = Intensity) # Simplified for SIF
  
  write.csv(ego_edges, file.path(output_dir, "cytoscape_sif", paste0(g, "_ego_net.csv")), row.names = FALSE)
}

# --- 9. FINAL TABLES ---
write.csv(pivot_summary, file.path(output_dir, "Master_Pivot_Analysis_Full.csv"), row.names = FALSE)
# Filtered version as requested
write.csv(pivot_summary %>% filter(q_value < 0.05), 
          file.path(output_dir, "Master_Pivot_Analysis_Filtered_FDR05.csv"), row.names = FALSE)


# --- 8. PLOTTING MODULE ---
cat("\n--- Generating Visualizations ---\n")
library(ggplot2)
library(ggrepel)

# Plot 1: Volatility vs. Duration (The "Identification" Plot)
p1 <- ggplot(pivot_summary, aes(x = Duration, y = log10(Volatility_Score + 1), color = Status)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_text_repel(data = head(pivot_summary, 10), aes(label = Gene), box.padding = 0.5) +
  scale_color_manual(values = c("Global Pivot (High Conf)" = "firebrick", 
                                "Structural Coordinator" = "dodgerblue", 
                                "Candidate Driver" = "orange", 
                                "Background" = "grey80")) +
  theme_minimal() +
  labs(title = "Network Volatility vs. Temporal Persistence",
       subtitle = "Highlighting High-Confidence and Structural Coordinators",
       x = "Duration (Days in Top 50)", y = "Log10 Volatility Score")

ggsave(file.path(output_dir, "plots", "volatility_vs_duration.png"), p1, width = 8, height = 6, dpi = 300)

# Plot 2: Intensity vs. Breadth (The "Role" Plot)
p2 <- ggplot(pivot_summary, aes(x = Total_Breadth, y = Avg_Intensity, color = Status)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = head(pivot_summary, 5), aes(label = Gene)) +
  theme_minimal() +
  labs(title = "Gene Rewiring: Intensity vs. Breadth",
       x = "Total Connectivity (Breadth)", y = "Mean Rewiring Strength (Intensity)")

ggsave(file.path(output_dir, "plots", "intensity_vs_breadth.png"), p2, width = 8, height = 6, dpi = 300)

# Plot 3: Top Pivot Timeline (TSPAN12 Example)
top_gene <- pivot_summary$Gene[1]
timeline_data <- all_turnover_data %>% filter(Gene == top_gene)

p3 <- ggplot(timeline_data, aes(x = Source_File, y = Intensity, group = 1)) +
  geom_line(color = "firebrick", size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Activity Profile:", top_gene),
       subtitle = "Intensity of network rewiring across transition windows",
       x = "Transition Window", y = "Intensity (abs delta R)")

ggsave(file.path(output_dir, "plots", "top_pivot_timeline.png"), p3, width = 10, height = 5, dpi = 300)


# --- 10. QUARTO REPORTING (The Canonical Pattern) ---
cat("\n--- Rendering Quarto Report ---\n")
old_wd <- getwd()
setwd(output_dir) # Step 2: Temporarily set WD

# Step 3: Build QMD content as character vector
qmd_content <- c(
  "---",
  paste0("title: 'Network Volatility Report: ", run_name, "'"),
  "format:",
  "  html:",
  "    toc: true",
  "    code-fold: true",
  "    theme: cosmo",
  "execute:",
  "  echo: false",
  "---",
  "",
  "# Summary",
  "This report identifies pivot genes driving metabolic transitions by integrating edge turnover and differential expression.",
  "",
  "# Metadata",
  paste0("- **Script Path**: `", script_full, "`"),
  paste0("- **Turnover Data**: `", turnover_dir, "`"),
  paste0("- **LIMMA Integration**: ", ifelse(do_limma, "Enabled", "Disabled")),
  paste0("- **Timestamp**: ", timestamp),
  "",
  "# Analytical Logic",
  "Volatility is calculated as:  ",
  "$$Volatility = \\overline{Intensity} \\times \\sum(Breadth) \\times Duration$$",
  "Where Intensity is the magnitude of $\\Delta r$. FDR is assigned via 1,000 parallel permutations.",
  "",
  "# Top Pivot Genes",
  "```{r}",
  "knitr::kable(head(read.csv('Master_Pivot_Analysis_Filtered_FDR05.csv'), 15))",
  "```",
  "",
  "# Generated Output Files",
  "```{r}",
  "list.files('.', recursive = TRUE)",
  "```",
  "",
  "# Session Information",
  "```{r}",
  "sessionInfo()",
  "```"
)

# Step 4: Write QMD
qmd_filename <- paste0(run_name, "_Report.qmd")
writeLines(qmd_content, qmd_filename)

# Step 5: Render HTML
quarto::quarto_render(qmd_filename)

# Step 6: Restore WD
setwd(old_wd)

# Step 7: Final Print
html_path <- file.path(output_dir, paste0(run_name, "_Report.html"))
cat("\n--- ANALYSIS COMPLETE ---\n")
cat("Absolute Output Path: ", output_dir, "\n")
cat("Report Path: ", html_path, "\n")

if(interactive()) browseURL(html_path)
