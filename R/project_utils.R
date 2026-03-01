# ======================================================================
# Project root + output directory enforcement
# ======================================================================

find_project_root <- function() normalizePath(getwd(), winslash = "/", mustWork = TRUE)


init_project_output <- function(subdir = NULL) {
  root <- find_project_root()
  out <- file.path(root, "outputs")
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  
  if (!is.null(subdir)) {
    out <- file.path(out, subdir)
    if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  }
  
  list(
    project_root = root,
    output_root = out
  )
}
