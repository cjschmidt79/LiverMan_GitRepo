# R/utils_io.R
# ======================================================================
# Utilities for structured Input/Output and Metadata Tracking
# ======================================================================

get_multiple_inputs <- function() {
  all_files_meta <- list()
  add_more <- TRUE
  
  while (add_more) {
    cat("\n>>> ACTION: Select an INPUT FILE (CSV)\n")
    cat("The file picker will open in 3 seconds...\n")
    Sys.sleep(3)
    
    selected_path <- rstudioapi::selectFile(
      caption = "Select an Input File",
      label = "Open"
    )
    
    if (is.null(selected_path)) {
      if (length(all_files_meta) == 0) message("No files selected.")
      break
    }
    
    # Store both the path and the name for the QMD report
    file_info <- list(
      path = selected_path,
      name = basename(selected_path)
    )
    
    all_files_meta[[length(all_files_meta) + 1]] <- file_info
    message("✅ Added: ", file_info$name)
    
    response <- utils::askYesNo("Would you like to add another input file?")
    if (is.na(response) || response == FALSE) {
      add_more <- FALSE
    }
  }
  return(all_files_meta)
}

setup_output_dir <- function() {
  cat("\n>>> ACTION: Select the MAIN FOLDER where you want to save your results\n")
  cat("The directory picker will open in 3 seconds...\n")
  Sys.sleep(3)
  
  base_path <- rstudioapi::selectDirectory(caption = "Select Folder for All Outputs")
  
  if (is.null(base_path)) stop("No folder selected. Script stopped.")
  
  # --- Capture Execution Script Info ---
  script_context <- rstudioapi::getActiveDocumentContext()
  script_path <- script_context$path
  script_name <- basename(script_path)
  
  if (script_name == "") {
    script_name <- "Console_Run"
    script_path <- "Interactive Session"
  }
  
  # --- Create Folder ---
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  folder_name <- paste0(gsub("\\.R$", "", script_name), "_", timestamp)
  full_output_path <- file.path(base_path, folder_name)
  
  dir.create(full_output_path, recursive = TRUE)
  
  # --- CRITICAL FIX: Initialize Log File Immediately ---
  log_path <- file.path(full_output_path, "execution_log.txt")
  init_msg <- paste(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "Log initialized. Folder created.")
  writeLines(init_msg, log_path)
  
  message("📂 Created output folder: ", full_output_path)
  
  # --- Return structured metadata ---
  return(list(
    path = full_output_path,
    log_path = log_path,  # Added this to the return for easy access
    timestamp = timestamp,
    script_info = list(
      name = script_name,
      path = script_path
    )
  ))
}

# --- Logging Functions ---

log_msg <- function(message, log_path) {
  # Creates a timestamped entry: [2024-05-20 10:00:00] Your Message
  ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  formatted_msg <- paste(ts, message)
  
  # Append to the file (creates it if it doesn't exist)
  cat(formatted_msg, file = log_path, sep = "\n", append = TRUE)
  
  # Also print to console so you see it live
  cat(formatted_msg, "\n")
}

# Helper to define the log file path
get_log_path <- function(output_dir) {
  file.path(output_dir, "execution_log.txt")
}