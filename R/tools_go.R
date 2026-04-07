## MCP tools: gene_go_terms, go_term_info, go_to_genes, compare_gene_go

# ── gene_go_terms ─────────────────────────────────────────────────────────────

fn_gene_go_terms <- function(gene_id, ontology = "ALL") {
  entrez  <- .lookup_entrez(gene_id)
  ont_arg <- toupper(trimws(ontology))
  if (!ont_arg %in% c("ALL", "BP", "MF", "CC"))
    stop("ontology must be one of: ALL, BP, MF, CC")

  res <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = entrez,
    columns = c("GO", "ONTOLOGY", "EVIDENCE"),
    keytype = "ENTREZID"
  ))
  res <- res[!is.na(res$GO), , drop = FALSE]
  if (nrow(res) == 0L)
    return(sprintf("No GO annotations found for '%s'.", gene_id))
  if (ont_arg != "ALL") {
    res <- res[!is.na(res$ONTOLOGY) & res$ONTOLOGY == ont_arg, , drop = FALSE]
    if (nrow(res) == 0L)
      return(sprintf("No %s GO annotations found for '%s'.", ont_arg, gene_id))
  }

  go_ids <- unique(res$GO)
  terms  <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys    = go_ids,
    columns = c("GOID", "TERM", "ONTOLOGY"),
    keytype = "GOID"
  ))
  merged <- merge(res, terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
  merged <- merged[!duplicated(merged$GO), , drop = FALSE]
  lines  <- sprintf(
    "[%s] %s: %s (evidence: %s)",
    ifelse(is.na(merged$ONTOLOGY.y), merged$ONTOLOGY.x, merged$ONTOLOGY.y),
    merged$GO,
    ifelse(is.na(merged$TERM), "unknown", merged$TERM),
    ifelse(is.na(merged$EVIDENCE), "?", merged$EVIDENCE)
  )
  paste0(
    sprintf("GO annotations for '%s' (Entrez %s), ontology=%s:\n",
            gene_id, entrez, ont_arg),
    paste(sort(lines), collapse = "\n")
  )
}

#' MCP tool: retrieve GO annotations for a human gene
#' @format An object of class `ToolDef`.
#' @export
tool_gene_go_terms <- ellmer::tool(
  fn_gene_go_terms,
  "Retrieve Gene Ontology (GO) annotations for a human gene from org.Hs.eg.db,
   including GO ID, term name, ontology namespace (BP/MF/CC), and evidence code.
   Filter with the ontology argument.",
  gene_id  = ellmer::type_string(
    "Gene symbol (e.g. 'TP53') or Entrez Gene ID (e.g. '7157')."
  ),
  ontology = ellmer::type_string(
    "Ontology filter: 'ALL' (default), 'BP' (Biological Process),
     'MF' (Molecular Function), or 'CC' (Cellular Component)."
  )
)

# ── go_term_info ──────────────────────────────────────────────────────────────

fn_go_term_info <- function(go_id) {
  go_id <- trimws(go_id)
  if (!grepl("^GO:[0-9]{7}$", go_id))
    stop(sprintf("'%s' is not a valid GO ID. Expected format: GO:0000000", go_id))
  res <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys    = go_id,
    columns = c("GOID", "TERM", "DEFINITION", "ONTOLOGY"),
    keytype = "GOID"
  ))
  if (nrow(res) == 0L) stop(sprintf("GO ID %s not found in GO.db.", go_id))
  r <- res[1L, ]
  paste(c(
    sprintf("GO ID:      %s", r$GOID),
    sprintf("Ontology:   %s", r$ONTOLOGY   %||% "NA"),
    sprintf("Term:       %s", r$TERM       %||% "NA"),
    sprintf("Definition: %s", r$DEFINITION %||% "NA")
  ), collapse = "\n")
}

#' MCP tool: look up GO term name and definition
#' @format An object of class `ToolDef`.
#' @export
tool_go_term_info <- ellmer::tool(
  fn_go_term_info,
  "Look up the name, definition, and ontology namespace (BP/MF/CC) for a
   Gene Ontology term by its GO ID.",
  go_id = ellmer::type_string(
    "A GO term identifier in the format GO:XXXXXXX, e.g. 'GO:0006915'."
  )
)

# ── go_to_genes ───────────────────────────────────────────────────────────────

fn_go_to_genes <- function(go_id, max_genes = 50L) {
  go_id     <- trimws(go_id)
  max_genes <- as.integer(max_genes)
  if (!grepl("^GO:[0-9]{7}$", go_id))
    stop(sprintf("'%s' is not a valid GO ID.", go_id))
  res <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = go_id,
    columns = c("ENTREZID", "SYMBOL", "EVIDENCE"),
    keytype = "GOALL"
  ))
  res <- res[!is.na(res$ENTREZID) & !duplicated(res$ENTREZID), , drop = FALSE]
  if (nrow(res) == 0L) return(sprintf("No human genes annotated to %s.", go_id))
  res   <- res[order(res$SYMBOL), ]
  total <- nrow(res)
  if (total > max_genes) res <- res[seq_len(max_genes), ]
  term_res  <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db, keys = go_id, columns = "TERM", keytype = "GOID"
  ))
  term_name <- if (nrow(term_res) > 0L) term_res$TERM[[1L]] else "unknown"
  paste0(
    sprintf("Genes annotated to %s (%s) — showing %d of %d:\n",
            go_id, term_name, nrow(res), total),
    paste(sprintf("%s (Entrez: %s) [%s]", res$SYMBOL, res$ENTREZID, res$EVIDENCE),
          collapse = "\n")
  )
}

#' MCP tool: find human genes annotated to a GO term
#' @format An object of class `ToolDef`.
#' @export
tool_go_to_genes <- ellmer::tool(
  fn_go_to_genes,
  "Find human genes annotated to a GO term (including annotations inherited
   through child terms via the transitive closure) using org.Hs.eg.db. Returns
   gene symbol, Entrez ID, and evidence code.",
  go_id     = ellmer::type_string(
    "A GO ID, e.g. 'GO:0006915' (apoptotic process)."
  ),
  max_genes = ellmer::type_integer(
    "Maximum number of genes to return. Defaults to 50."
  )
)

# ── compare_gene_go ───────────────────────────────────────────────────────────

fn_compare_gene_go <- function(gene_id_1, gene_id_2, ontology = "BP") {
  ont_arg <- toupper(trimws(ontology))
  if (!ont_arg %in% c("BP", "MF", "CC", "ALL"))
    stop("ontology must be one of: BP, MF, CC, ALL")

  get_go <- function(gid) {
    entrez <- .lookup_entrez(gid)
    res    <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = entrez,
      columns = c("GO", "ONTOLOGY"),
      keytype = "ENTREZID"
    ))
    if (ont_arg != "ALL")
      res <- res[!is.na(res$ONTOLOGY) & res$ONTOLOGY == ont_arg, , drop = FALSE]
    unique(stats::na.omit(res$GO))
  }

  fmt_terms <- function(ids) {
    if (length(ids) == 0L) return("  (none)")
    t <- suppressMessages(AnnotationDbi::select(
      GO.db::GO.db, keys = ids,
      columns = c("GOID", "TERM"), keytype = "GOID"
    ))
    paste(sprintf("  %s: %s", t$GOID, t$TERM), collapse = "\n")
  }

  go1    <- get_go(gene_id_1)
  go2    <- get_go(gene_id_2)
  shared <- intersect(go1, go2)
  only1  <- setdiff(go1, go2)
  only2  <- setdiff(go2, go1)

  paste0(
    sprintf("GO comparison: %s vs %s  (ontology: %s)\n",
            gene_id_1, gene_id_2, ont_arg),
    sprintf("\nShared terms (%d):\n",          length(shared)), fmt_terms(shared), "\n",
    sprintf("\nOnly in %s (%d):\n", gene_id_1, length(only1)),  fmt_terms(only1),  "\n",
    sprintf("\nOnly in %s (%d):\n", gene_id_2, length(only2)),  fmt_terms(only2)
  )
}

#' MCP tool: compare GO annotations between two genes
#' @format An object of class `ToolDef`.
#' @export
tool_compare_gene_go <- ellmer::tool(
  fn_compare_gene_go,
  "Compare GO annotations of two human genes to find shared and unique GO
   terms. Useful for assessing functional similarity. Accepts gene symbols or
   Entrez IDs.",
  gene_id_1 = ellmer::type_string("First gene — symbol or Entrez ID."),
  gene_id_2 = ellmer::type_string("Second gene — symbol or Entrez ID."),
  ontology  = ellmer::type_string(
    "Ontology: 'BP' (default), 'MF', 'CC', or 'ALL'."
  )
)
