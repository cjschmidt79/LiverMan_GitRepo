#!/usr/bin/env Rscript
# ============================================================
# Generic Violin Plot Utility
#
# PURPOSE:
#   Create violin plots from any CSV with user-selected columns.
#
# USAGE (interactive):
#   source("Generic_ViolinPlot.R")
#
# USAGE (CLI):
#   Rscript Generic_ViolinPlot.R \
#     --input=data.csv \
#     --x=Day \
#     --y=Value \
#     --fill=Group \
#     --out=violin.png
#
# PARAMETERS:
#   --input   CSV file path
#   --x       Column name for x-axis (factor)
#   --y       Column name for y-axis (numeric)
#   --fill    OPTIONAL column for coloring violins
#   --out     Output PNG filename
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

# ----------------------------- CLI parsing -----------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag) {
  hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^--", flag, "="), "", hit[1])
}

infile <- parse_arg("input")
xcol   <- parse_arg("x")
ycol   <- parse_arg("y")
fillcol <- parse_arg("fill")
outfile <- parse_arg("out")

# Interactive fallback
if (interactive()) {
  if (is.na(infile)) {
    cat("Select input CSV\n")
    infile <- file.choose()
  }
  if (is.na(xcol) || is.na(ycol)) {
    stop("Interactive mode still requires --x and --y")
  }
}

# Defaults
if (is.na(outfile)) outfile <- "violin_plot.png"

# ----------------------------- Validation -----------------------------
if (!file.exists(infile)) stop("Input file not found: ", infile)

df <- read_csv(infile, show_col_types = FALSE)

if (!xcol %in% names(df)) stop("x column not found: ", xcol)
if (!ycol %in% names(df)) stop("y column not found: ", ycol)
if (!is.na(fillcol) && !fillcol %in% names(df)) {
  stop("fill column not found: ", fillcol)
}

df <- df[is.finite(df[[ycol]]), ]

df[[xcol]] <- as.factor(df[[xcol]])
if (!is.na(fillcol)) df[[fillcol]] <- as.factor(df[[fillcol]])

# ----------------------------- Plot construction -----------------------------
aes_map <- aes_string(x = xcol, y = ycol)
if (!is.na(fillcol)) aes_map$fill <- as.name(fillcol)

y_min <- min(df[[ycol]], na.rm = TRUE)
y_max <- max(df[[ycol]], na.rm = TRUE)
y_pad <- 0.05 * (y_max - y_min)

p <- ggplot(df, aes_map) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "grey30",
    linewidth = 0.3
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad)) +
  labs(
    x = xcol,
    y = ycol,
    fill = fillcol
  ) +
  theme_classic(base_size = 12)

# ----------------------------- Save -----------------------------
ggsave(outfile, p, width = 8, height = 5, dpi = 300)
print(p)

cat("Saved violin plot to:", outfile, "\n")