
#removes the quotes from the edge table
#make sure the node import table has: name as the common ID 

cat("\nSelect the folder containing edge CSVs (edges/)\n")
edge_dir <- if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  rstudioapi::selectDirectory("Select edges/ folder")
} else {
  readline("Enter full path to edges folder: ")
}

files <- list.files(edge_dir, pattern = "\\.csv$", full.names = TRUE)
files <- files[grepl("edges_slice_|edges_filtered_combined|edges_all_combined", basename(files))]

strip_quotes <- function(x) {
  x <- as.character(x)
  x <- gsub('^"+|"+$', "", x)
  x <- gsub("^'+|'+$", "", x)
  trimws(x)
}

for (f in files) {
  edges <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Only fix if source/target exist
  if (!all(c("source", "target") %in% names(edges))) next
  
  edges$source <- strip_quotes(edges$source)
  edges$target <- strip_quotes(edges$target)
  
  out_f <- file.path(edge_dir, paste0(sub("\\.csv$", "", basename(f)), "_FIXED.csv"))
  write.csv(edges, out_f, row.names = FALSE, quote = FALSE)
  cat("Fixed: ", basename(out_f), "\n", sep = "")
}

cat("\nDone. Import the *_FIXED.csv files into Cytoscape.\n")