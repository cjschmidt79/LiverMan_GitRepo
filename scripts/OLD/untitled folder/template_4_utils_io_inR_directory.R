#The Universal R Launcher TemplateR# ======================================================================
# [PROJECT NAME] LAUNCHER
# Description: [Briefly describe what this analysis does]
# Date: 2026-XX-XX
# ======================================================================

cat("\n======================================================================\n")
cat("STEP 1: Load Utility and Logic Modules\n")
cat("======================================================================\n")

# --- 1. Load Universal Utils (usually stays the same) ---
cat("\n>>> ACTION: Select 'utils_io.R'\n")
# If utils_io.R is always in the same place, you can hardcode the path here
utils_path <- file.choose() 
source(utils_path)
message("✅ Loaded: utils_io.R")

# --- 2. Load Specific Logic Module (CHANGE THIS for new scripts) ---
cat("\n>>> ACTION: Select your specific Logic/Analysis Script\n")
# TODO: Rename this variable to match your analysis type (e.g., gsea_module)
module_path <- file.choose() 
source(module_path)
message("✅ Loaded Logic Module")


cat("\n======================================================================\n")
cat("STEP 2: Project Environment & Logging\n")
cat("======================================================================\n")

# 1. Setup Output & Get Log Path
# (Assumes setup_output_dir() is defined in your utils_io.R)
out_meta <- setup_output_dir()
log_f <- get_log_path(out_meta$path) 

# 2. Get Input Files
in_meta <- get_multiple_inputs()

if (length(in_meta) == 0) {
  log_msg("ERROR: No input files selected. Aborting.", log_f)
  stop("❌ No input files selected. Analysis aborted.")
}

# 3. Log Initial Details
log_msg(paste("Project:", "[INSERT PROJECT NAME HERE]"), log_f)
log_msg(paste("User selected input file:", in_meta[[1]]$name), log_f)
log_msg(paste("Output directory set to:", out_meta$path), log_f)


cat("\n======================================================================\n")
cat("STEP 3: Run Analysis\n")
cat("======================================================================\n")

# TODO: Call the main function from your sourced module here
# Example: run_main_analysis(in_meta, out_meta, log_f)

log_msg("SUCCESS: Analysis script completed.", log_f)
cat("\n✅ Process Complete. Check log for details.\n")
#🛠 What to Change for Each New ScriptTo make this template specific to a new project, follow this checklist:SectionTargetActionHeader[PROJECT NAME]Update the title so you know which launcher you are looking at.Step 1module_pathUpdate the prompt message to tell you which specific logic file to pick.Step 2log_msgUpdate the hardcoded string "Project: ..." to match your current work.Step 3The Function CallThis is the most important part. You 
#need to call the actual function defined inside your "Logic Module."