# ==========================================================
# 🖼️ COMBINE PNG IMAGES INTO A FACETED GRID
# ==========================================================

# Required packages
library(grid)
library(gridExtra)
library(png)
library(tools)
library(magick)

# -----------------------------------------
# Ask user to choose input mode
# -----------------------------------------
cat("\U0001F4C1 Select PNG import mode:\n")
cat("1: Select PNGs manually (one by one)\n")
cat("2: Load all PNGs from a folder\n")
mode <- readline("Enter 1 or 2: ")

# -----------------------------------------
# Collect PNG files
# -----------------------------------------
png_files <- character()

if (mode == "1") {
    repeat {
        cat("\U0001F4C4 Select a PNG file (or hit ENTER to stop):\n")
        file <- tryCatch(file.choose(), error = function(e) "")
        if (file == "") break
        if (!grepl("\\.png$", file, ignore.case = TRUE)) {
            cat("\u26A0\uFE0F Not a PNG file. Please try again.\n")
            next
        }
        png_files <- c(png_files, file)
    }
} else if (mode == "2") {
    cat("\U0001F4C2 Select ANY file inside your target directory:\n")
    any_file <- file.choose()
    dir_path <- dirname(any_file)
    png_files <- list.files(dir_path, pattern = "\\.png$", full.names = TRUE)
    cat("\u2705 Found", length(png_files), "PNG files in directory.\n")
} else {
    stop("\u274C Invalid mode. Please enter 1 or 2.")
}

if (length(png_files) == 0) stop("\u274C No PNG files selected.")

# -----------------------------------------
# Ask for output directory
# -----------------------------------------
cat("\n\U0001F4C2 Select ANY file inside your desired output directory:\n")
base_dir <- dirname(file.choose())
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- file.path(base_dir, paste0("combined_plot_", timestamp, ".png"))

# -----------------------------------------
# Ask layout preferences
# -----------------------------------------
cat("\n\U0001F9E9 Layout: Combining your images into a grid of rows and columns.\n")
cat("   - You control the number of columns (images per row).\n")
cat("   - The script will automatically calculate how many rows are needed.\n")
cat("   - For example:\n")
cat("       3 columns = 3 images per row (good for wide layout)\n")
cat("       1 column  = single column (stacked vertically)\n")
cat("       Leave blank to auto-select based on number of images.\n")
cat("\n\U0001F522 Enter number of columns for layout (default = auto): ")
col_input <- readline()

if (col_input == "") {
    n_col <- ceiling(sqrt(length(png_files)))
} else {
    n_col <- as.integer(col_input)
}
n_row <- ceiling(length(png_files) / n_col)

# -----------------------------------------
# Ask if user wants custom titles
# -----------------------------------------
cat("\n\U0001F4DD Do you want to add custom titles above each image? (yes/no): ")
add_titles <- tolower(trimws(readline()))

titles <- NULL
if (add_titles %in% c("yes", "y")) {
    titles <- character(length(png_files))
    cat("\u270F\uFE0F Enter a short title for each image:\n")
    for (i in seq_along(png_files)) {
        default_title <- file_path_sans_ext(basename(png_files[i]))
        prompt <- paste0("  • Title for '", basename(png_files[i]), "' [default: ", default_title, "]: ")
        title <- readline(prompt)
        if (title == "") title <- default_title
        titles[i] <- title
    }
}

# -----------------------------------------
# Save combined PNGs and PDFs with paging
# -----------------------------------------
n_row_per_page <- 3
images_per_page <- n_col * n_row_per_page
page_chunks <- split(png_files, ceiling(seq_along(png_files) / images_per_page))

for (page_idx in seq_along(page_chunks)) {
    this_page <- page_chunks[[page_idx]]
    this_titles <- if (!is.null(titles)) titles[seq_along(this_page) + (page_idx - 1) * images_per_page] else NULL
    
    grobs <- lapply(seq_along(this_page), function(i) {
        trimmed_file <- image_read(this_page[i]) %>%
            image_trim(fuzz = 5) %>%
            image_write(tempfile(fileext = ".png"))
        img <- readPNG(trimmed_file)
        raster <- rasterGrob(
            img,
            interpolate = TRUE,
            width  = unit(1, "npc"),   # <-- fills cell width
            height = unit(1, "npc")    # <-- fills cell height (kills vertical gaps)
        )
        if (!is.null(this_titles)) {
            title_grob <- textGrob(this_titles[i], gp = gpar(fontsize = 9, fontface = "bold"))
            return(arrangeGrob(
                title_grob,
                raster,
                ncol = 1,
                heights = unit(c(3, 1), c("mm", "null"))
            )
            )
        } else {
            return(raster)
        }
    })
    
    width_px <- 1800 * n_col
    height_px <- 1800 * n_row_per_page
    suffix <- paste0("_page", page_idx)
    output_file_page <- sub("\\.png$", paste0(suffix, ".png"), output_file)
    pdf_output_file_page <- sub("\\.png$", paste0(suffix, ".pdf"), output_file)
    
    cat("\U0001F9E9 Creating layout for page", page_idx, "with", length(grobs), "images...\n")
    
    # Build the layout ONCE (no divergence possible)
    layout_grob <- arrangeGrob(grobs = grobs, ncol = n_col, padding = unit(0, "mm"))  ## <-- change this Adjust 0 to e.g. 0.2, 0.5, 1:
    
    # --- Save PNG ---
    png(output_file_page, width = width_px, height = height_px, units = "px", res = 600)
    grid.draw(layout_grob)
    dev.off()
    cat("\u2705 Saved PNG:", output_file_page, "\n")
    
    # --- Save PDF (same layout_grob) ---
    pdf(pdf_output_file_page, width = width_px / 600, height = height_px / 600)
    grid.draw(layout_grob)
    dev.off()
    cat("\u2705 Saved PDF:", pdf_output_file_page, "\n")
    
}