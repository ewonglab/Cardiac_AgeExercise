# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Resolve configured R workflow paths from named environment variables.
# Steps: 1) Read a configured root. 2) Fail clearly when it is unset. 3) Join optional path components reproducibly.
# =======================================
publication_path <- function(variable, ...) {
  root <- Sys.getenv(variable, unset = "")
  if (!nzchar(root)) stop(sprintf("Set %s in config/.env before running this workflow.", variable), call. = FALSE)
  components <- list(...)
  if (!length(components)) return(root)
  do.call(file.path, c(list(root), components))
}
