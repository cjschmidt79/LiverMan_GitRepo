# ============================================================
# Contract utilities
# ============================================================

write_script_contract_readme <- function(output_dir,
                                         contract,
                                         filename = "README_outputs.txt") {
  stopifnot(is.character(output_dir), length(output_dir) == 1)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!is.list(contract) || is.null(names(contract)) || any(names(contract) == "")) {
    stop("`contract` must be a *named list* with non-empty names.")
  }
  
  norm_val <- function(v) {
    if (is.null(v)) return(character(0))
    if (length(v) == 1 && is.na(v)) return(character(0))
    if (is.character(v)) return(v)
    as.character(v)
  }
  
  canonical <- c(
    "SCRIPT NAME",
    "BIOLOGICAL QUESTION",
    "BIOLOGICAL CLAIM SUPPORTED",
    "UPSTREAM DEPENDENCIES",
    "DOWNSTREAM CONSUMERS",
    "FIGURES FED",
    "MANUSCRIPT PLACEMENT",
    "CANALIZATION DEFINITION SUPPORTED",
    "KEY OUTPUT FILES"
  )
  
  ordered_keys <- c(intersect(canonical, names(contract)),
                    setdiff(names(contract), canonical))
  
  lines <- character(0)
  for (k in ordered_keys) {
    lines <- c(lines, paste0(k, ":"), norm_val(contract[[k]]), "")
  }
  
  out_path <- file.path(output_dir, filename)
  writeLines(lines, out_path)
  invisible(out_path)
}
