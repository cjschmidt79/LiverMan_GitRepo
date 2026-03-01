###############################################################################
# Script Name: CompositeScores_NetworkAnalysis.R
#
# Purpose:
#   This script reads a Cytoscape "Network Analysis" description table and
#   computes integrated network scores for each node (e.g., metabolite, protein,
#   gene). These scores combine multiple standard network metrics into six
#   biologically meaningful categories: HUB, Bottleneck, InfluenceScore,
#   VulnerabilityScore, Modularity, and SignalIntegration.
#
# Expected Input:
#   - CSV file exported from Cytoscape’s "Network Analyzer" (or similar tool).
#   - File should have:
#       * First column: node identifier (e.g., "name" or "Metabolite").
#       * Other columns: network metrics such as:
#           Degree, NeighborhoodConnectivity, ClosenessCentrality, Radiality,
#           BetweennessCentrality, Stress, AverageShortestPathLength, Eccentricity,
#           ClusteringCoefficient, TopologicalCoefficient
#   - Extra columns will be ignored. Missing metrics are allowed (script adapts).
#
# Output:
#   All results are saved in a timestamped output directory:
#       1) ZScores_<timestamp>.csv
#       2) CompositeScores_<timestamp>.csv
#       3) CompositeScores_Categorized_<timestamp>.csv
#       4) Figures (PNGs + interactive HTMLs)

###############################################################################
# Interpretation Guide for Scatter Plots (HUB vs Bottleneck, HUB vs Modularity)
#
# These plots show composite category means with:
#   - X-axis = HUB (Mean z-score)
#   - Y-axis = Bottleneck or Modularity (Mean z-score)
#   - Color  = SignalIntegration (Mean z-score; blue=low, gray≈0, red=high)
#   - Dashed line = regression fit (trend line summarizing average relationship)
#   - Solid axes = zero lines (split quadrants at network average)
#
# Quadrants (using axis zero lines):
#   Q1 (HUB > 0, Y > 0) → "Super nodes"
#       - High hubness and high bottleneck/modularity
#       - Likely anchors or critical routers
#   Q2 (HUB < 0, Y > 0) → "Bridges/gatekeepers"
#       - Not hubs, but still embedded (Modularity) or chokepoints (Bottleneck)
#       - Can disconnect modules without being globally central
#   Q3 (HUB < 0, Y < 0) → "Peripheral/background"
#       - Weakly connected, low flow control, low embedding
#       - Typically low biological leverage
#   Q4 (HUB > 0, Y < 0) → "Cross-module hubs"
#       - High connectivity, but not strongly modular or flow-controlling
#       - Potential global connectors between clusters
#
# Trend line interpretation (dashed LM fit):
#   - Nodes ON the line → follow expected relationship between HUB and Y
#   - Nodes ABOVE line → stronger Y-role than expected given their hubness
#       * Example: at HUB ≈ –2 and Bottleneck ≈ –2, nodes are still "above line"
#         → weakly connected, weak flow controllers in absolute terms,
#           but slightly MORE bottlenecking than expected.
#         → interpreted as "peripheral weak links" or "local connectors"
#   - Nodes BELOW line → weaker Y-role than expected
#       * Example: high HUB but low Bottleneck → connectors that are not chokepoints
#
# Color interpretation (SignalIntegration):
#   - Red points → cross-module integrators (bridges, global gatekeepers)
#   - Blue points → local/module-specific roles
#   - Gray points → neutral, average integration
#
# Biological summary:
#   - High HUB + Modularity → module anchors (intra-module leaders)
#   - High HUB + Bottleneck → super routers (globally critical nodes)
#   - Moderate HUB + high Bottleneck (above line) → bridges/gatekeepers
#   - Low/low scores (Q3) → background nodes; only local weak links if above line
###############################################################################

#
###############################################################################

###############################################################################
# Script Name: CompositeScores_NetworkAnalysis.R (base R version, no dplyr/tidyr)
###############################################################################

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(scales)
  library(GGally)
  library(plotly)
  library(htmlwidgets)
})

# -------------------------
# Helper functions
# -------------------------

.safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

# Embed helpers
.embed_image <- function(file, width = "80%") {
  paste0("![](", file, "){width=", width, "}")
}
.embed_iframe <- function(file, height = "650px") {
  paste0("<iframe src=\"", file, "\" width=\"100%\" height=\"", height, "\" frameborder=\"0\"></iframe>")
}

# -------------------------
# 1. Select input CSV
# -------------------------
if (interactive()) {
  cat("\n📂 Please select your Cytoscape description table (CSV):\n")
  file_path <- file.choose()
} else {
  stop("This script must be run interactively to select a file.")
}

# -------------------------
# 2. Select output directory
# -------------------------
cat("\n📂 Select output directory:\n")
base_dir <- dirname(file.choose())

input_name <- tools::file_path_sans_ext(basename(file_path))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_dir, paste0(input_name, "_", timestamp))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("✅ Output will be saved in:", output_dir, "\n")

# -------------------------
# 3. Category definitions
# -------------------------
metric_pool <- c("Degree","NeighborhoodConnectivity","ClosenessCentrality","Radiality",
                 "BetweennessCentrality","Stress","AverageShortestPathLength","Eccentricity",
                 "ClusteringCoefficient","TopologicalCoefficient")

categories <- list(
  HUB = c("Degree","NeighborhoodConnectivity","ClosenessCentrality","Radiality"),
  Bottleneck = c("BetweennessCentrality","Stress","AverageShortestPathLength","Eccentricity"),
  InfluenceScore = c("ClosenessCentrality","Radiality","Stress","BetweennessCentrality"),
  VulnerabilityScore = c("BetweennessCentrality","Stress","Eccentricity","ClusteringCoefficient"),
  Modularity = c("Degree","NeighborhoodConnectivity","TopologicalCoefficient","ClusteringCoefficient"),
  SignalIntegration = c("Degree","TopologicalCoefficient","NeighborhoodConnectivity","BetweennessCentrality")
)

# -------------------------
# 4. Load + clean input
# -------------------------
raw <- readr::read_csv(file_path, show_col_types = FALSE)

name_candidates <- c("Metabolite","name","Name","Node","ID","id")
name_col <- intersect(name_candidates, colnames(raw))
name_col <- if (length(name_col) == 0) colnames(raw)[1] else name_col[1]

present_metrics <- intersect(metric_pool, colnames(raw))
if (length(present_metrics) == 0) stop("No approved metric columns found in the input file.")

DATA <- raw[, c(name_col, present_metrics), drop = FALSE]
colnames(DATA)[1] <- "Metabolite"
DATA$Metabolite <- as.factor(DATA$Metabolite)

# Drop zero-variance metrics
sd_vals <- sapply(DATA[, -1, drop = FALSE], function(x) sd(as.numeric(x), na.rm = TRUE))
keep_cols <- names(sd_vals)[sd_vals > 0]
DATA <- DATA[, c("Metabolite", keep_cols), drop = FALSE]

# -------------------------
# 5. Z-scores + flips
# -------------------------
numeric_data <- as.data.frame(lapply(DATA[, -1, drop = FALSE], as.numeric))
z_data <- as.data.frame(scale(numeric_data))
z_data <- cbind(Metabolite = DATA$Metabolite, z_data)

if ("AverageShortestPathLength" %in% colnames(z_data)) z_data$AverageShortestPathLength <- -z_data$AverageShortestPathLength
if ("Eccentricity" %in% colnames(z_data)) z_data$Eccentricity <- -z_data$Eccentricity
if ("ClusteringCoefficient" %in% colnames(z_data)) z_data$ClusteringCoefficient <- -z_data$ClusteringCoefficient

# -------------------------
# 6. Composite calculator
# -------------------------
calculate_composites <- function(df, metrics, name) {
  metrics <- intersect(metrics, colnames(df))
  scores <- df[, metrics, drop = FALSE]
  out <- data.frame(
    Metabolite = df$Metabolite,
    Sum = rowSums(scores, na.rm = TRUE),
    Mean = rowMeans(scores, na.rm = TRUE),
    Product = apply(scores, 1, function(x) prod(replace(x, x == 0, 1e-10), na.rm = TRUE))
  )
  colnames(out)[2:4] <- paste0(name, "_", c("Sum","Mean","Product"))
  out
}

composite_scores <- z_data["Metabolite"]
for (cat in names(categories)) {
  comp <- calculate_composites(z_data, categories[[cat]], cat)
  composite_scores <- merge(composite_scores, comp, by = "Metabolite")
}

# -------------------------
# 7. Classification  (SAFE / robust to tied quantiles)
# -------------------------
metric_counts <- sapply(categories, function(m) length(intersect(m, colnames(z_data))))

classify_score <- function(score, type, category) {
    if (type == "Mean") {
        if (score > 1) "High" else if (score < -1) "Low" else "Moderate"
    } else if (type == "Sum") {
        n <- metric_counts[[category]]
        if (score > n) "High" else if (score < -n) "Low" else "Moderate"
    } else {
        score
    }
}

# Robust tertile labeling that never produces duplicate breaks
safe_tertiles <- function(x) {
    v <- x[is.finite(x)]
    if (!length(v) || length(unique(v)) <= 1L) {
        # all equal or empty -> everything Moderate
        out <- rep("Moderate", length(x))
        out[!is.finite(x)] <- NA
        return(out)
    }
    qs <- as.numeric(quantile(v, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE, type = 7))
    if (!is.finite(qs[1]) || !is.finite(qs[2]) || qs[1] >= qs[2]) {
        # fallback: median split with center band at exact median
        med <- median(v, na.rm = TRUE)
        out <- ifelse(x < med, "Low", ifelse(x > med, "High", "Moderate"))
        out[!is.finite(x)] <- NA
        return(out)
    }
    br <- c(-Inf, qs, Inf)  # strictly increasing by construction
    out <- cut(x, breaks = br, labels = c("Low", "Moderate", "High"),
               include.lowest = TRUE, right = TRUE)
    as.character(out)
}

class_df <- data.frame(Metabolite = composite_scores$Metabolite)
for (cat in names(categories)) {
    class_df[[paste0(cat, "_Mean_Class")]] <- sapply(
        composite_scores[[paste0(cat,"_Mean")]], classify_score, type="Mean", category=cat
    )
    class_df[[paste0(cat, "_Sum_Class")]]  <- sapply(
        composite_scores[[paste0(cat,"_Sum")]], classify_score, type="Sum", category=cat
    )
    prod_vals <- composite_scores[[paste0(cat,"_Product")]]
    class_df[[paste0(cat, "_Product_Class")]] <- safe_tertiles(prod_vals)
}


# -------------------------
# 8. Save outputs
# -------------------------
write.csv(z_data, file.path(output_dir, paste0("ZScores_", timestamp, ".csv")), row.names = FALSE)
write.csv(composite_scores, file.path(output_dir, paste0("CompositeScores_", timestamp, ".csv")), row.names = FALSE)
write.csv(class_df, file.path(output_dir, paste0("CompositeScores_Categorized_", timestamp, ".csv")), row.names = FALSE)

# -------------------------
# 9. Visualization (with tracking)
# -------------------------
generated_pngs <- c()
generated_htmls <- c()

# Top-N barplots
mean_cols <- grep("_Mean$", names(composite_scores), value = TRUE)
cs_mean_long <- reshape(composite_scores[, c("Metabolite", mean_cols)],
                        varying = mean_cols, v.names = "Score",
                        timevar = "Category", times = gsub("_Mean","",mean_cols),
                        new.row.names = 1:(nrow(composite_scores)*length(mean_cols)),
                        direction = "long")

TOP_N <- 15L
topN <- do.call(rbind, lapply(split(cs_mean_long, cs_mean_long$Category), function(df) {
  df[order(-df$Score), ][1:min(TOP_N, nrow(df)), ]
}))
topN$Metabolite <- factor(topN$Metabolite, levels = unique(topN$Metabolite[order(topN$Score)]))

p_top <- ggplot(topN, aes(x=Metabolite, y=Score, fill=Category)) +
  geom_col(show.legend=FALSE) + coord_flip() +
  facet_wrap(~Category, scales="free_y", ncol=2) +
  labs(title=paste0("Top ", TOP_N, " by Category (Mean composite)"),
       x=NULL,y="Mean Z-score") + theme_minimal()
out_top <- file.path(output_dir, paste0("Top", TOP_N, "_CategoryMeans_", timestamp, ".png"))
ggsave(out_top, p_top, width=10, height=10, dpi=300, bg="white")
generated_pngs <- c(generated_pngs, basename(out_top))

# Violin plot
p_violin <- ggplot(cs_mean_long, aes(x = Category, y = Score, fill = Category)) +
  geom_violin(trim = TRUE, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.9, show.legend = FALSE) +
  scale_y_continuous(labels = number_format(accuracy = 0.1)) +
  labs(title = "Distribution of Category Mean Scores", x = NULL, y = "Mean Z-score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
out_violin <- file.path(output_dir, paste0("CategoryMeans_Distributions_", timestamp, ".png"))
ggsave(out_violin, p_violin, width = 10, height = 6, dpi = 300, bg = "white")
generated_pngs <- c(generated_pngs, basename(out_violin))

# Wide table of means
wide_means <- reshape(cs_mean_long,
                      idvar = "Metabolite",
                      timevar = "Category",
                      direction = "wide")
colnames(wide_means) <- sub("Score.", "", colnames(wide_means))

plot_scatter <- function(df, xcat, ycat, color_cat = "SignalIntegration") {
  has_color <- color_cat %in% names(df)
  p <- ggplot(df, aes_string(x = xcat, y = ycat)) +
    geom_hline(yintercept = 0, colour = "gray60", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "gray60", linewidth = 0.3) +
    (if (has_color) geom_point(aes_string(color = color_cat), size = 2, alpha = 0.85)
     else geom_point(color = "gray50", size = 2, alpha = 0.85)) +
    geom_smooth(method = "lm", se = FALSE, color = "gray40", linetype = "dashed") +
    (if (has_color) scale_color_gradient2(
      name = paste0(color_cat, " (Mean)\nblue=low • gray≈0 • red=high"),
      low = "#2c7bb6", mid = "gray80", high = "#d7191c", midpoint = 0
    ) else NULL) +
    labs(title = paste(ycat, "vs", xcat, "(Mean)"),
         x = paste0(xcat, " (Mean)"),
         y = paste0(ycat, " (Mean)")) +
    theme_minimal(base_size = 11)
  
  out <- file.path(output_dir, paste0("Scatter_", .safe_name(xcat), "_vs_", .safe_name(ycat), "_", timestamp, ".png"))
  ggsave(out, p, width = 7, height = 6, dpi = 300, bg = "white")
  generated_pngs <<- c(generated_pngs, basename(out))
}

plot_scatter(wide_means, "HUB", "Modularity")
plot_scatter(wide_means, "HUB", "Bottleneck")
if (requireNamespace("GGally", quietly = TRUE)) {
  K <- 30L
  avg_means <- rowMeans(wide_means[, setdiff(names(wide_means), "Metabolite")], na.rm = TRUE)
  top_idx <- order(avg_means, decreasing = TRUE)[1:min(K, length(avg_means))]
  topK <- wide_means[top_idx, ]
  
  p_pc <- GGally::ggparcoord(
    data = topK,
    columns = which(names(topK) != "Metabolite"),
    groupColumn = NULL,
    showPoints = FALSE,
    scale = "globalminmax",
    alphaLines = 0.6
  ) +
    theme_minimal(base_size = 11) +
    labs(title = paste0("Parallel Coordinates (Top ", K, " by Avg Mean)"),
         x = "Category", y = "Mean Z-score") +
    theme(legend.position = "none")
  
  out_pc <- file.path(output_dir, paste0("ParallelCoords_Top", K, "_", timestamp, ".png"))
  ggsave(out_pc, p_pc, width = 10, height = 7, dpi = 300, bg = "white")
  generated_pngs <- c(generated_pngs, basename(out_pc))
}
if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
  make_interactive_scatter <- function(df, x, y, color = "SignalIntegration", file_stub) {
    tooltip <- paste0("Metabolite: ", df$Metabolite,
                      "<br>", x, ": ", round(df[[x]], 3),
                      "<br>", y, ": ", round(df[[y]], 3),
                      if (color %in% names(df)) paste0("<br>", color, ": ", round(df[[color]], 3)) else "")
    
    p <- plotly::plot_ly(
      x = df[[x]], y = df[[y]],
      type = "scatter", mode = "markers",
      color = if (color %in% names(df)) df[[color]] else NULL,
      colors = c("#2c7bb6", "gray80", "#d7191c"),
      text = tooltip,
      hovertemplate = "%{text}<extra></extra>",
      marker = list(size = 8, opacity = 0.85, line = list(width = 1, color = "rgba(0,0,0,0.3)"))
    ) %>% plotly::layout(
      title = paste0(y, " vs ", x, " (Mean)"),
      xaxis = list(title = paste0(x, " (Mean)")),
      yaxis = list(title = paste0(y, " (Mean)")),
      dragmode = "lasso"
    )
    
    out <- file.path(output_dir, paste0(file_stub, "_", timestamp, ".html"))
    htmlwidgets::saveWidget(p, file = out, selfcontained = TRUE)
    generated_htmls <<- c(generated_htmls, basename(out))
  }
  
  make_interactive_scatter(wide_means, "HUB", "Modularity", "SignalIntegration",
                           "Interactive_Scatter_Modularity_vs_HUB")
  make_interactive_scatter(wide_means, "HUB", "Bottleneck", "SignalIntegration",
                           "Interactive_Scatter_Bottleneck_vs_HUB")
}
if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
  if (all(c("HUB","Bottleneck","Modularity","SignalIntegration") %in% names(wide_means))) {
    p <- plotly::plot_ly(
      data = wide_means,
      x = ~HUB, y = ~Bottleneck,
      type = "scatter", mode = "markers",
      color = ~SignalIntegration, colors = c("#2c7bb6", "gray80", "#d7191c"),
      size = ~pmax(0.1, abs(Modularity)) + 2, sizes = c(6, 24),
      text = ~paste0("Metabolite: ", Metabolite,
                     "<br>HUB: ", round(HUB, 3),
                     "<br>Bottleneck: ", round(Bottleneck, 3),
                     "<br>Modularity: ", round(Modularity, 3),
                     "<br>SignalIntegration: ", round(SignalIntegration, 3)),
      hoverinfo = "text"
    ) %>% plotly::layout(
      title = "HUB vs Bottleneck (bubble size = |Modularity|, color = SignalIntegration)",
      xaxis = list(title = "HUB (Mean)"),
      yaxis = list(title = "Bottleneck (Mean)")
    )
    
    out_bubble <- file.path(output_dir, paste0("Bubble_HUB_vs_Bottleneck_", timestamp, ".html"))
    htmlwidgets::saveWidget(p, file = out_bubble, selfcontained = TRUE)
    generated_htmls <- c(generated_htmls, basename(out_bubble))
  }
}


# Export tracking
assign("generated_pngs", generated_pngs, envir = .GlobalEnv)
assign("generated_htmls", generated_htmls, envir = .GlobalEnv)

message("✅ Analysis complete. Files tracked for report builder.")

# =========================
# QMD builder + auto-render
# =========================

# --- guards ---
if (!exists("output_dir")) stop("output_dir not found")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- small helpers (no regex) ---
embed_image  <- function(file, width = "90%") sprintf("![](%s){width=%s}", file, width)
embed_iframe <- function(file, height = "650px") sprintf('<iframe src="%s" width="100%%" height="%s" frameborder="0"></iframe>', file, height)

get_script_path <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    f  <- sub("^--file=", "", ca[grepl("^--file=", ca)])
    if (length(f) == 1 && file.exists(f)) return(normalizePath(f, winslash = "/", mustWork = FALSE))
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
        p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
        if (nzchar(p) && file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
    NA_character_
}

capture_header <- function(path) {
    if (is.na(path) || !file.exists(path)) return(character(0))
    ln <- readLines(path, warn = FALSE)
    run <- character(0)
    for (x in ln) {
        if (grepl("^\\s*#", x)) run <- c(run, sub("^\\s*#\\s?", "", x)) else break
    }
    run
}

# --- detect assets already written by the script ---
ls_out <- function(pat) {
    ff <- list.files(output_dir, pattern = pat, full.names = FALSE)
    ff[order(ff)]
}

all_files       <- list.files(output_dir, full.names = FALSE)
z_csvs          <- ls_out("^ZScores_.*\\.csv$")
comp_csvs       <- ls_out("^CompositeScores_\\d{8}_\\d{6}\\.csv$")
class_csvs      <- ls_out("^CompositeScores_Categorized_.*\\.csv$")
top_plots       <- ls_out("^Top\\d+_CategoryMeans_.*\\.png$")
violin_plots    <- ls_out("^CategoryMeans_Distributions_.*\\.png$")
scatter_plots   <- ls_out("^Scatter_.*\\.png$")
parallel_plots  <- ls_out("^ParallelCoords_.*\\.png$")
bubble_htmls    <- ls_out("^Bubble_HUB_vs_Bottleneck_.*\\.html$")
scatter_htmls   <- ls_out("^Interactive_Scatter_.*\\.html$")

# --- metadata + header ---
script_path <- get_script_path()
script_name <- if (!is.na(script_path)) basename(script_path) else "CompositeScores_NetworkAnalysis.R"
input_stub  <- {
    # try to infer original input name from folder like "<input>_YYYYMMDD_HHMMSS"
    b <- basename(output_dir)
    sub("_[0-9]{8}_[0-9]{6}$", "", b)
}
header_lines <- capture_header(script_path)

# --- build QMD content ---
yaml <- c(
    "---",
    sprintf('title: "%s — Report"', tools::file_path_sans_ext(script_name)),
    "format:",
    "  html:",
    "    toc: true",
    "    toc-depth: 3",
    "    number-sections: false",
    "    code-fold: show",
    "    theme: cosmo",
    "editor: visual",
    "---",
    ""
)

sec_979 <- c(
    "# 979. Summary", "",
    "- This report compiles outputs from the network composite scoring workflow: inputs, methods, generated files, visual summaries (static + interactive), a shortlist of multi-role high nodes, and session details.", ""
)

sec_980 <- c(
    "# 980. Script Header", "",
    if (length(header_lines)) c("```text", header_lines, "```", "") else "_No commented header detected._", ""
)

sec_981 <- c(
    "# 981. Metadata", "",
    sprintf("- **Script name:** `%s`", script_name),
    sprintf("- **Script path:** `%s`", if (!is.na(script_path)) script_path else "(unknown)"),
    sprintf("- **Output directory:** `%s`", normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
    sprintf("- **Input (inferred):** `%s`", input_stub),
    if (exists("file_path")) sprintf("- **Input path:** `%s`", file_path) else NULL,
    if (exists("timestamp")) sprintf("- **Timestamp:** `%s`", timestamp) else NULL,
    ""
)

sec_982 <- c(
    "# 982. Dependencies", "",
    "- readr", "- ggplot2", "- scales", "- GGally", "- plotly", "- htmlwidgets", "- knitr", ""
)

sec_983 <- c(
    "# 983. Generated Output Files", "",
    "```",
    if (length(all_files)) all_files else "(no files found in output_dir)",
    "```",
    ""
)

sec_984 <- c(
    "# 984. Statistical / Analytical Logic", "",
    "**Z-scoring & sign alignment**: cost-like metrics (AverageShortestPathLength, Eccentricity, ClusteringCoefficient where treated as costs) are sign-flipped so “higher is better” before aggregation.", "",
    "**Composites per category C (metrics m ∈ C, node i):**", "",
    "- **Sum:** $S_i^{(C)} = \\sum_{m\\in C} z_{i,m}$",
    "- **Mean:** $\\bar{S}_i^{(C)} = \\frac{1}{|C|}\\sum_{m\\in C} z_{i,m}$",
    "- **Product (stabilized):** $P_i^{(C)} = \\prod_{m\\in C} \\max(z_{i,m}, \\varepsilon)$ with $\\varepsilon=10^{-10}$.", "",
    "**Classification**: Mean High(>1)/Low(<−1); Sum High(>|C|)/Low(<−|C|); Product by tertiles (25%, 75%).", ""
)

# figures / plots
sec_985 <- c("# 985. Figures / Plots", "")
plots <- c(top_plots, violin_plots, scatter_plots, parallel_plots)
if (length(plots)) {
    for (f in plots) sec_985 <- c(sec_985, paste0("### ", tools::file_path_sans_ext(f)), embed_image(f), "")
} else {
    sec_985 <- c(sec_985, "_No static plots detected._", "")
}
sec_985 <- c(sec_985, "## Interactive Visualizations", "")
htmls <- c(scatter_htmls, bubble_htmls)
if (length(htmls)) {
    for (f in htmls) sec_985 <- c(sec_985, paste0("### ", tools::file_path_sans_ext(f)), embed_iframe(f), "")
} else {
    sec_985 <- c(sec_985, "_No interactive HTML plots detected._", "")
}

sec_986 <- c(
    "# 986. Interpretation of Results", "",
    "- **HUB**: high connectivity/reach; anchors.",
    "- **Bottleneck**: flow controllers/chokepoints.",
    "- **InfluenceScore**: reach × control.",
    "- **VulnerabilityScore**: sensitivity to node loss.",
    "- **Modularity**: embedding in cohesive modules.",
    "- **SignalIntegration**: cross-module bridging.", "",
    "**Quadrants (HUB vs Bottleneck/Modularity)**: Q1 super nodes; Q2 bridges/gatekeepers; Q3 peripheral; Q4 cross-module hubs.", ""
)

# multi-role shortlist (reads categorized CSV if present)
sec_multi <- c(
    "# Multi-role High Nodes (Prioritization Shortlist)", "",
    "```{r}",
    "suppressPackageStartupMessages(library(knitr))",
    "cats <- c('HUB','Bottleneck','InfluenceScore','VulnerabilityScore','Modularity','SignalIntegration')",
    "pat  <- '^CompositeScores_Categorized_.*\\\\.csv$'",

    "cand <- list.files('.', pattern = pat, full.names = TRUE)",
    "if (length(cand)) {",
    "  x <- tryCatch(read.csv(cand[order(file.info(cand)$mtime, decreasing = TRUE)[1]], check.names = FALSE), error = function(e) NULL)",
    "  if (!is.null(x)) {",
    "    mcols <- intersect(paste0(cats, '_Mean_Class'), names(x))",
    "    if (length(mcols)) {",
    "      x$High_Count <- rowSums(x[mcols] == 'High', na.rm = TRUE)",
    "      y <- subset(x, High_Count >= 3, select = c('Metabolite', mcols, 'High_Count'))",
    "      if (nrow(y)) knitr::kable(y, caption = 'Nodes scoring High in ≥3 categories (Mean composites)') else cat('No nodes met the ≥3 High criteria.')",
    "    } else cat('Expected columns not found in categorized CSV.')",
    "  } else cat('Could not read categorized composite CSV.')",
    "} else cat('No categorized composite CSV found in output_dir.')",
    "```",
    ""
)

sec_987 <- c("# 987. Session Info", "", "```{r}", "sessionInfo()", "```", "")

# --- write QMD next to outputs ---
qmd_name <- "CompositeScores_NetworkAnalysis_Report.qmd"
qmd_path <- file.path(output_dir, qmd_name)
qmd <- c(
    yaml,
    sec_979,
    sec_980,
    sec_981,
    sec_982,
    sec_983,
    sec_984,
    sec_985,
    sec_multi,
    sec_987
)
writeLines(qmd, qmd_path)
message("✓ QMD written: ", qmd_path)

# --- render HTML in the same folder (works with older quarto pkgs) ---
render_ok <- FALSE
if (requireNamespace("quarto", quietly = TRUE)) {
    ow <- getwd(); on.exit(setwd(ow), add = TRUE)
    setwd(output_dir)
    try({
        quarto::quarto_render(input = basename(qmd_path), quiet = TRUE)
        render_ok <- file.exists(sub("\\.qmd$", ".html", qmd_path))
    }, silent = TRUE)
}
if (!render_ok && nzchar(Sys.which("quarto"))) {
    # fallback to CLI
    try(system2("quarto", c("render", shQuote(qmd_path))), silent = TRUE)
    render_ok <- file.exists(sub("\\.qmd$", ".html", qmd_path))
}

if (render_ok) {
    message("✓ HTML rendered: ", sub("\\.qmd$", ".html", qmd_path))
} else {
    message("⚠ Could not render automatically. Try: quarto render \"", qmd_path, "\"")
}
