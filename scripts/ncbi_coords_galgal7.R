# ===================================================================
# GRCg7b (galGal7) Local Gene Coordinator - ROBUST COLUMN VERSION
#Must have a GFF3 file from NCBI: galGal7_annotation.gff
# ===================================================================

# ---- 1. Dependencies & Setup ----
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("rtracklayer")
}
if (!requireNamespace("quarto", quietly = TRUE)) install.packages("quarto")

library(rtracklayer)

# Provenance Setup
script_full <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) "Unknown")
script_name <- basename(script_full)
timestamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- 2. Interactive Inputs ----
message("\n====================================================")
message(">>> STEP 1: SELECT YOUR GENE LIST")
message("====================================================\n")

gene_file <- rstudioapi::selectFile(caption = "Select Gene CSV/TSV")
if (is.null(gene_file)) stop("Cancelled.")

raw_input <- if(tools::file_ext(gene_file) == "csv") read.csv(gene_file) else read.delim(gene_file)
message("Columns found: ", paste(colnames(raw_input), collapse = ", "))
gene_col <- readline("Type the EXACT name of the column with GENE SYMBOLS: ")

message("\n====================================================")
message(">>> STEP 2: SELECT THE galGal7 GFF3")
message("====================================================\n")

gff_file <- rstudioapi::selectFile(caption = "Select genomic.gff.gz")
if (is.null(gff_file)) stop("Cancelled.")

run_name <- readline("Enter a name for this run: ")
if(run_name == "") run_name <- "canalization"

# ---- 3. Output Directory ----
output_dir <- normalizePath(file.path(getwd(), "outputs", paste0(run_name, "_", timestamp)), mustWork = FALSE)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 4. Robust Extraction Logic ----
message("\n>>> Importing GFF3 (this can take ~60s)...")
gff_raw <- import(gff_file)
gff_df  <- as.data.frame(gff_raw)

# --- STEP 4a: Robust Chromosome Mapping ---
message(">>> Mapping RefSeq IDs to Chromosome names...")

# Look for any column that is named 'chromosome' (case-insensitive)
col_idx <- grep("^chromosome$", colnames(gff_df), ignore.case = TRUE)

if(length(col_idx) > 0) {
  chr_col_name <- colnames(gff_df)[col_idx[1]]
  # Map RefSeq (seqnames) to the found chromosome column
  chr_map <- unique(gff_df[!is.na(gff_df[[chr_col_name]]), c("seqnames", chr_col_name)])
  colnames(chr_map) <- c("seqnames", "CleanChr")
} else {
  message("!!! Warning: No 'chromosome' attribute column found. Using Accession IDs.")
  chr_map <- data.frame(seqnames = unique(gff_df$seqnames), CleanChr = unique(gff_df$seqnames))
}

# --- STEP 4b: Extract Genes ---
target_list <- unique(as.character(raw_input[[gene_col]]))
# Base R filtering
coords <- gff_df[gff_df$type == "gene" & (gff_df$gene %in% target_list | gff_df$Name %in% target_list), ]

# --- STEP 4c: Merge and Cleanup ---
if (nrow(coords) == 0) stop("No genes from your list were found in this GFF3 file.")

# Force merge the CleanChr name onto our coordinate table
final_coords <- merge(coords, chr_map, by = "seqnames", all.x = TRUE)

# Ensure the column exists before selecting to avoid 'undefined columns' error
final_table <- data.frame(
  Chromosome = if("CleanChr" %in% colnames(final_coords)) final_coords$CleanChr else final_coords$seqnames,
  Start      = final_coords$start,
  End        = final_coords$end,
  Strand     = final_coords$strand,
  Symbol     = if("gene" %in% colnames(final_coords)) final_coords$gene else final_coords$Name
)

# Replace NAs in Chromosome with seqnames (Accessions)
final_table$Chromosome[is.na(final_table$Chromosome)] <- as.character(final_coords$seqnames[is.na(final_table$Chromosome)])

# ---- 5. Save Outputs ----
write.table(final_table, file.path(output_dir, "resolved_coords.tsv"), sep="\t", quote=F, row.names=F)

bed_out <- final_table[, c("Chromosome", "Start", "End", "Symbol")]
bed_out$Start <- bed_out$Start - 1
write.table(bed_out, file.path(output_dir, "targets.bed"), sep="\t", quote=F, row.names=F, col.names=F)

# ---- 6. Quarto Report ----
qmd_content <- c(
  "---",
  paste0("title: 'Canalization Mapping: ", run_name, "'"),
  "format: html",
  "---",
  "## Top 20 Genes Found",
  "```{r echo=FALSE}",
  "df <- read.delim('resolved_coords.tsv')",
  "knitr::kable(head(df, 20))",
  "```"
)
writeLines(qmd_content, file.path(output_dir, "report.qmd"))

old_wd <- getwd()
setwd(output_dir)
tryCatch({ quarto::quarto_render("report.qmd") }, finally = { setwd(old_wd) })

if (interactive()) {
  message(">>> SUCCESS. Results in: ", output_dir)
  View(final_table)
  browseURL(file.path(output_dir, "report.html"))
}