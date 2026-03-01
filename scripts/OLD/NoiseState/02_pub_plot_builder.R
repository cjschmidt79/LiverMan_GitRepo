#!/usr/bin/env Rscript
# ======================================================================
# Constructed Phase Diagram with Uncertainty Envelopes
#
# X-axis: D(t) median
# Y-axis: ED(t) median
# Rectangles: marginal CI for D and ED at each day
# Arrows: connect consecutive days (ordered by Day)
#
# This is a DESCRIPTIVE, CONSTRUCTED state-space trajectory.
# No continuous dynamics are assumed.
# ======================================================================

# ----------------------------- Utilities -----------------------------

ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
stop2 <- function(...) stop(paste0(...), call. = FALSE)

# ----------------------------- Packages ------------------------------

needed <- c("ggplot2", "dplyr", "grid")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# ----------------------------- Helpers -------------------------------

choose_csv <- function() {
  if (!interactive()) stop2("Run interactively so file.choose() can be used.")
  message("Select DAY-LEVEL summary CSV (with med / lo / hi columns)")
  file.choose()
}

show_cols <- function(df) {
  for (i in seq_along(names(df))) {
    cat(sprintf("  [%d] %s\n", i, names(df)[i]))
  }
}

choose_column <- function(df, prompt) {
  cat("\n", prompt, "\n", sep = "")
  show_cols(df)
  idx <- suppressWarnings(as.integer(readline("Enter column number: ")))
  if (is.na(idx) || idx < 1 || idx > ncol(df)) stop2("Invalid column selection.")
  names(df)[idx]
}

# ----------------------------- Main ----------------------------------

cat("\n=== Constructed Phase Diagram (with Uncertainty Envelopes) ===\n")

csv_path <- choose_csv()
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

cat("\nLoaded: ", csv_path, "\n", sep = "")
cat("Rows: ", nrow(df), " | Columns: ", ncol(df), "\n", sep = "")

# ---- Explicit column selection (NO AUTO-DETECTION) ----

day_col   <- choose_column(df, "Select DAY column:")
d_col     <- choose_column(df, "Select D MEDIAN column:")
d_lo_col  <- choose_column(df, "Select D LOWER CI column:")
d_hi_col  <- choose_column(df, "Select D UPPER CI column:")
ed_col    <- choose_column(df, "Select ED MEDIAN column:")
ed_lo_col <- choose_column(df, "Select ED LOWER CI column:")
ed_hi_col <- choose_column(df, "Select ED UPPER CI column:")

# ---- Coerce numeric ----

num_cols <- c(day_col, d_col, d_lo_col, d_hi_col, ed_col, ed_lo_col, ed_hi_col)
for (cc in num_cols) {
  df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
}

# ---- Clean + order ----

df2 <- df %>%
  filter(
    is.finite(.data[[day_col]]),
    is.finite(.data[[d_col]]),
    is.finite(.data[[ed_col]])
  ) %>%
  distinct(.data[[day_col]], .keep_all = TRUE) %>%
  arrange(.data[[day_col]])

if (nrow(df2) < 3) stop2("Need at least 3 days for a trajectory.")

# ---- Build arrow segments ----

seg <- df2 %>%
  mutate(
    d_next  = dplyr::lead(.data[[d_col]]),
    ed_next = dplyr::lead(.data[[ed_col]])
  ) %>%
  filter(is.finite(d_next), is.finite(ed_next))

# ----------------------------- Plot ----------------------------------

title_txt <- "Constructed Phase Diagram (State-Space Trajectory)"
subtitle_txt <- paste0(
  "X = ", d_col, " | Y = ", ed_col,
  " | arrows connect consecutive days | rectangles = marginal CI\n",
  "Input: ", basename(csv_path)
)

p <- ggplot(df2, aes(x = .data[[d_col]], y = .data[[ed_col]])) +
  
  # Uncertainty envelopes (RECTANGLES)
  geom_rect(
    aes(
      xmin = .data[[d_lo_col]],
      xmax = .data[[d_hi_col]],
      ymin = .data[[ed_lo_col]],
      ymax = .data[[ed_hi_col]]
    ),
    fill = "grey80",
    alpha = 0.6,
    color = NA
  ) +
  
  # Arrowed trajectory
  geom_segment(
    data = seg,
    aes(
      xend = d_next,
      yend = ed_next
    ),
    arrow = arrow(type = "closed", length = unit(3, "mm")),
    linewidth = 0.8,
    lineend = "round"
  ) +
  
  geom_point(size = 2.8) +
  
  {
    if (has_ggrepel) {
      ggrepel::geom_text_repel(aes(label = .data[[day_col]]), size = 4)
    } else {
      geom_text(aes(label = .data[[day_col]]), vjust = -0.8, size = 4)
    }
  } +
  
  labs(
    title = title_txt,
    subtitle = subtitle_txt,
    x = "D (dispersion)",
    y = "ED (effective dimensionality)"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

print(p)

# ----------------------------- Export --------------------------------

out_dir <- getwd()
stem <- paste0("Phase_D_vs_ED_with_Uncertainty_", ts_stamp())

png_path <- file.path(out_dir, paste0(stem, ".png"))
pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))

ggsave(png_path, p, width = 6.8, height = 6.0, dpi = 600)
ggsave(pdf_path, p, width = 6.8, height = 6.0)

cat("\nSaved:\n  ", png_path, "\n  ", pdf_path, "\n", sep = "")
cat("Done.\n")
