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

# Look up GO IDs in GO.db and return a data.frame with go_id, term, ontology.
# Returns an empty data.frame with those columns when ids is length 0.
.lookup_go_terms <- function(ids) {
  if (length(ids) == 0L)
    return(data.frame(go_id = character(), term = character(),
                      ontology = character(), stringsAsFactors = FALSE))
  t <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db, keys = ids,
    columns = c("GOID", "TERM", "ONTOLOGY"), keytype = "GOID"
  ))
  data.frame(go_id = t$GOID, term = t$TERM, ontology = t$ONTOLOGY,
             stringsAsFactors = FALSE)
}


#' Sanitize a scalar argument that an LLM may have JSON-array-wrapped.
#' e.g. '["GO:0010464"]' -> "GO:0010464"
#'      "[50]"           -> "50"
.sanitize_scalar <- function(x) {
  if (!is.character(x) || length(x) != 1) return(x)
  x <- trimws(x)
  parsed <- tryCatch(jsonlite::fromJSON(x), error = function(e) x)
  if (length(parsed) == 1) return(as.character(parsed))
  x  # return original if parse yields unexpected length
}

#' Validate that a string is a syntactically and ontologically valid GO ID.
#' Returns TRUE only if the ID matches GO:XXXXXXX and exists in GO.db.
.is_valid_go_id <- function(go_id) {
  if (!grepl("^GO:\\d{7}$", go_id)) return(FALSE)
  go_id %in% AnnotationDbi::keys(GO.db, keytype = "GOID")
}

#' Uniform negative-result list returned when inputs are nonsensical.
.na_result <- function(message) {
  list(result = NA_character_, message = message)
}


