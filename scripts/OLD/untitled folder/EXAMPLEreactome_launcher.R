# ======================================================================
# reactome_launcher.R (MANUAL SELECTION + LOGGING)
# ======================================================================

cat("\n======================================================================\n")
cat("STEP 1: Locate your Logic Scripts\n")
cat("======================================================================\n")

# --- Ask for utils_io.R ---
cat("\n>>> ACTION: Select 'utils_io.R'\n")
utils_path <- file.choose()
source(utils_path)
message("✅ Loaded: utils_io.R")

# --- Ask for reactome_camera_tidy_module.R ---
cat("\n>>> ACTION: Select 'reactome_camera_tidy_module.R'\n")
module_path <- file.choose()
source(module_path)
message("✅ Loaded: reactome_camera_tidy_module.R")


cat("\n======================================================================\n")
cat("STEP 2: Define Project Paths\n")
cat("======================================================================\n")

# 1. Setup Output & Get Log Path
out_meta <- setup_output_dir()
log_f <- get_log_path(out_meta$path) # Connects to the log file

# 2. Get Inputs
in_meta <- get_multiple_inputs()

if (length(in_meta) == 0) {
  log_msg("ERROR: No input files selected. Aborting.", log_f)
  stop("❌ No input files selected. Analysis aborted.")
}

# 3. Log the Project Start
log_msg(paste("User selected input file:", in_meta[[1]]$name), log_f)
log_msg(paste("Output directory set to:", out_meta$path), log_f)


cat("\n======================================================================\n")
cat("STEP 3: Run the Analysis\n")
cat("======================================================================\n")

# Assemble the manifest for the report
project_manifest <- list(
  script_name = out_meta$script_info$name,
  script_path = out_meta$script_info$path,
  output_dir  = out_meta$path,
  input_files = in_meta
)

# Execute the camera pipeline
log_msg("Starting Camera Pipeline math...", log_f)

run_reactome_camera_pipeline_tidy(
  input_csv  = project_manifest$input_files[[1]]$path, 
  output_dir = project_manifest$output_dir,
  manifest   = project_manifest,
  make_report = TRUE
)

log_msg("Pipeline run complete. Report generated.", log_f)

cat("\n======================================================================\n")
cat("🎉 ALL TASKS COMPLETE\n")
cat("Check your output folder for CSVs, the log, and the Quarto report.\n")
cat("======================================================================\n")