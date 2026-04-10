## MCP tools: search_go_terms (local GO.db), search_go_terms_ols (EBI OLS4)

# ── search_go_terms (local, cached) ──────────────────────────────────────────

get_go_term_table <- function() {
  if (!is.null(.pkg_env$go_term_cache)) return(.pkg_env$go_term_cache)
  tbl <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys    = AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"),
    columns = c("GOID", "TERM", "DEFINITION", "ONTOLOGY"),
    keytype = "GOID"
  ))
  .pkg_env$go_term_cache <- tbl
  tbl
}

fn_search_go_terms <- function(query, ontology = "ALL", max_results = 20L) {
  ont_arg     <- toupper(trimws(ontology))
  max_results <- as.integer(max_results)
  if (!ont_arg %in% c("ALL", "BP", "MF", "CC"))
    stop("ontology must be one of: ALL, BP, MF, CC")

  tbl <- get_go_term_table()
  if (ont_arg != "ALL")
    tbl <- tbl[!is.na(tbl$ONTOLOGY) & tbl$ONTOLOGY == ont_arg, ]

  term_hit <- grepl(query, tbl$TERM,       ignore.case = TRUE, perl = TRUE)
  defn_hit <- grepl(query, tbl$DEFINITION, ignore.case = TRUE, perl = TRUE)
  score    <- term_hit * 2L + defn_hit

  hits       <- tbl[score > 0L, ]
  hits_score <- score[score > 0L]
  total      <- nrow(hits)

  if (total == 0L)
    return(list(query = query, ontology = ont_arg, total = 0L, results = list()))

  ord  <- order(-hits_score, hits$TERM)
  hits <- hits[ord, ]
  if (nrow(hits) > max_results) hits <- hits[seq_len(max_results), ]

  results <- lapply(seq_len(nrow(hits)), function(i) list(
    go_id      = hits$GOID[[i]],
    ontology   = hits$ONTOLOGY[[i]],
    term       = hits$TERM[[i]],
    definition = if (is.na(hits$DEFINITION[[i]])) NA_character_
                 else substr(hits$DEFINITION[[i]], 1L, 200L)
  ))

  list(
    query    = query,
    ontology = ont_arg,
    total    = total,
    results  = results
  )
}

#' MCP tool: search GO terms by keyword (local, offline)
#'
#' Searches `TERM` and `DEFINITION` fields in GO.db. Results are cached after
#' the first call. Term-name matches are ranked above definition-only matches.
#'
#' @format An object of class `ToolDef`.
#' @export
tool_search_go_terms <- ellmer::tool(
  fn_search_go_terms,
  description = "Search for Gene Ontology terms matching a keyword or phrase by grepping
   term names and definitions in the local GO.db database. Works offline.
   Returns a structured object with query, ontology, total match count, and a
   'results' array. Each result has go_id, ontology, term, and definition fields.
   Term-name matches are ranked above definition-only matches.",
  arguments = list(
    query       = ellmer::type_string(
      "A keyword or phrase, e.g. 'apoptosis', 'kinase activity', 'DNA damage response'."
    ),
    ontology    = ellmer::type_string(
      "Restrict to: 'BP', 'MF', 'CC', or 'ALL' (default)."
    ),
    max_results = ellmer::type_integer(
      "Maximum GO terms to return. Defaults to 20."
    )
  )
)

# ── search_go_terms_ols (EBI OLS4 REST API) ───────────────────────────────────

fn_search_go_terms_ols <- function(query, ontology = "ALL", max_results = 20L) {
  ont_arg     <- toupper(trimws(ontology))
  max_results <- as.integer(max_results)
  if (!ont_arg %in% c("ALL", "BP", "MF", "CC"))
    stop("ontology must be one of: ALL, BP, MF, CC")

  # Over-fetch to allow post-hoc namespace filtering
  fetch_rows <- if (ont_arg == "ALL") max_results else min(max_results * 4L, 200L)

  full_url <- paste0(
    "https://www.ebi.ac.uk/ols4/api/search",
    "?q=",          utils::URLencode(query, reserved = TRUE),
    "&ontology=go",
    "&type=class",
    "&rows=",       fetch_rows,
    "&fieldList=obo_id,label,description,ontology_prefix,short_form"
  )

  result <- tryCatch(
    jsonlite::fromJSON(full_url, simplifyDataFrame = TRUE),
    error = function(e) stop(sprintf(
      "OLS4 API request failed for '%s': %s\nTip: use search_go_terms() for offline search.",
      query, e$message
    ))
  )

  docs <- result$response$docs
  if (is.null(docs) || nrow(docs) == 0L)
    return(list(query = query, ontology = ont_arg, total = 0L, results = list()))

  # Keep only well-formed GO IDs
  valid <- !is.na(docs$obo_id) & grepl("^GO:[0-9]{7}$", docs$obo_id)
  docs  <- docs[valid, , drop = FALSE]
  if (nrow(docs) == 0L)
    return(list(query = query, ontology = ont_arg, total = 0L, results = list()))

  # Enrich with namespace from local GO.db (authoritative, fast)
  ns_tbl  <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys    = docs$obo_id,
    columns = c("GOID", "ONTOLOGY"),
    keytype = "GOID"
  ))
  ns_map  <- stats::setNames(ns_tbl$ONTOLOGY, ns_tbl$GOID)
  docs$ns <- ns_map[docs$obo_id]

  if (ont_arg != "ALL") {
    docs <- docs[!is.na(docs$ns) & docs$ns == ont_arg, , drop = FALSE]
    if (nrow(docs) == 0L)
      return(list(query = query, ontology = ont_arg, total = 0L, results = list()))
  }

  total_hits <- result$response$numFound
  if (nrow(docs) > max_results) docs <- docs[seq_len(max_results), ]

  fmt_desc <- function(d) {
    if (is.null(d) || (length(d) == 1L && is.na(d))) return(NA_character_)
    if (is.list(d)) d <- unlist(d)
    substr(paste(d, collapse = " "), 1L, 200L)
  }

  results <- lapply(seq_len(nrow(docs)), function(i) {
    row  <- docs[i, ]
    list(
      go_id       = row$obo_id,
      ontology    = if (!is.na(row$ns)) row$ns else NA_character_,
      term        = row$label %||% NA_character_,
      description = if ("description" %in% names(docs)) fmt_desc(docs$description[[i]])
                    else NA_character_
    )
  })

  list(
    query    = query,
    ontology = ont_arg,
    total    = total_hits,
    results  = results
  )
}

#' MCP tool: search GO terms via EBI OLS4 API (network, fuzzy-ranked)
#'
#' Uses BM25 full-text search with synonym expansion. Requires network access.
#' Namespace (BP/MF/CC) is resolved locally from GO.db after retrieval.
#'
#' @format An object of class `ToolDef`.
#' @export
tool_search_go_terms_ols <- ellmer::tool(
  fn_search_go_terms_ols,
  description = "Search the EBI Ontology Lookup Service (OLS4) for GO terms matching a
   natural language description. Uses BM25 ranking with synonym expansion.
   Returns a structured object with query, ontology, total match count, and a
   'results' array. Each result has go_id, ontology, term, and description
   fields. Requires network access.",
  arguments = list(
    query       = ellmer::type_string(
      "Natural language phrase, e.g. 'programmed cell death', 'chromatin remodeling',
       'response to oxidative stress', 'protein folding in ER'."
    ),
    ontology    = ellmer::type_string(
      "Restrict results to: 'BP', 'MF', 'CC', or 'ALL' (default)."
    ),
    max_results = ellmer::type_integer(
      "Maximum GO terms to return. Defaults to 20."
    )
  )
)
