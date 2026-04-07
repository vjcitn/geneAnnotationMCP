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
    return(sprintf(
      "No GO terms matched '%s' (ontology=%s). Try broader keywords or use search_go_terms_ols().",
      query, ont_arg
    ))

  ord  <- order(-hits_score, hits$TERM)
  hits <- hits[ord, ]
  if (nrow(hits) > max_results) hits <- hits[seq_len(max_results), ]

  lines <- sprintf(
    "[%s] %s (%s)\n  %s",
    hits$ONTOLOGY, hits$GOID, hits$TERM,
    ifelse(is.na(hits$DEFINITION), "(no definition)",
           substr(hits$DEFINITION, 1L, 140L))
  )
  paste0(
    sprintf("GO term search for '%s' (ontology=%s) — %d of %d results:\n\n",
            query, ont_arg, nrow(hits), total),
    paste(lines, collapse = "\n\n")
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
  "Search for Gene Ontology terms matching a keyword or phrase by grepping
   term names and definitions in the local GO.db database. Works offline.
   Term-name matches are ranked above definition-only matches. Use
   search_go_terms_ols() for fuzzy or synonym-aware search.",
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
    return(sprintf(
      "OLS4 returned no GO matches for '%s'. Try different keywords or use search_go_terms().",
      query
    ))

  # Keep only well-formed GO IDs
  valid <- !is.na(docs$obo_id) & grepl("^GO:[0-9]{7}$", docs$obo_id)
  docs  <- docs[valid, , drop = FALSE]
  if (nrow(docs) == 0L)
    return(sprintf("OLS4 results for '%s' contained no valid GO IDs.", query))

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
      return(sprintf(
        "OLS4 found GO results for '%s' but none were in the '%s' namespace.",
        query, ont_arg
      ))
  }

  total_hits <- result$response$numFound
  if (nrow(docs) > max_results) docs <- docs[seq_len(max_results), ]

  fmt_desc <- function(d) {
    if (is.null(d) || (length(d) == 1L && is.na(d))) return("(no definition)")
    if (is.list(d)) d <- unlist(d)
    substr(paste(d, collapse = " "), 1L, 140L)
  }

  lines <- vapply(seq_len(nrow(docs)), function(i) {
    row  <- docs[i, ]
    desc <- if ("description" %in% names(docs)) fmt_desc(docs$description[[i]])
            else "(no definition)"
    sprintf("[%s] %s (%s)\n  %s",
            if (!is.na(row$ns)) row$ns else "?",
            row$obo_id,
            row$label %||% "unknown",
            desc)
  }, character(1L))

  paste0(
    sprintf("OLS4 GO search for '%s' (ontology=%s) — showing %d of %d total hits:\n\n",
            query, ont_arg, nrow(docs), total_hits),
    paste(lines, collapse = "\n\n")
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
  "Search the EBI Ontology Lookup Service (OLS4) for GO terms matching a
   natural language description. Uses BM25 ranking with synonym expansion,
   handling paraphrased queries better than the local grep-based search.
   Requires network access. Use search_go_terms() for offline keyword search.",
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
