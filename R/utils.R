## Internal utilities shared across all tool files.

# Package-level cache environment (avoids repeated GO.db full-table scans)
.pkg_env <- new.env(parent = emptyenv())

#' @importFrom stats na.omit setNames
#' @importFrom utils URLencode
NULL

# Null-coalescing operator (internal only)
`%||%` <- function(a, b) {
  if (length(a) == 1L && !is.na(a)) a else b
}

# Accept a gene symbol OR Entrez ID; always returns a character Entrez ID.
.lookup_entrez <- function(id) {
  id <- trimws(as.character(id))
  if (grepl("^[0-9]+$", id)) return(id)
  hits <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = toupper(id),
    columns = "ENTREZID",
    keytype = "SYMBOL"
  ))
  valid <- hits$ENTREZID[!is.na(hits$ENTREZID)]
  if (length(valid) == 0L)
    stop(sprintf("Symbol '%s' not found in org.Hs.eg.db.", id))
  valid[[1L]]
}
