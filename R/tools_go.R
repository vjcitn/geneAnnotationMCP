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
  if (ont_arg != "ALL") {
    res <- res[!is.na(res$ONTOLOGY) & res$ONTOLOGY == ont_arg, , drop = FALSE]
  }

  if (nrow(res) == 0L)
    return(list(gene_id = gene_id, entrez_id = entrez, ontology = ont_arg,
                annotations = list()))

  go_ids <- unique(res$GO)
  terms  <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys    = go_ids,
    columns = c("GOID", "TERM", "ONTOLOGY"),
    keytype = "GOID"
  ))
  merged <- merge(res, terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
  merged <- merged[!duplicated(merged$GO), , drop = FALSE]
  merged <- merged[order(merged$GO), , drop = FALSE]

  annotations <- lapply(seq_len(nrow(merged)), function(i) list(
    go_id    = merged$GO[[i]],
    ontology = if (is.na(merged$ONTOLOGY.y[[i]])) merged$ONTOLOGY.x[[i]]
               else merged$ONTOLOGY.y[[i]],
    term     = if (is.na(merged$TERM[[i]])) NA_character_ else merged$TERM[[i]],
    evidence = if (is.na(merged$EVIDENCE[[i]])) NA_character_ else merged$EVIDENCE[[i]]
  ))

  list(
    gene_id     = gene_id,
    entrez_id   = entrez,
    ontology    = ont_arg,
    annotations = annotations
  )
}

#' MCP tool: retrieve GO annotations for a human gene
#' @format An object of class `ToolDef`.
#' @export
tool_gene_go_terms <- ellmer::tool(
  fn_gene_go_terms,
  description = "Retrieve Gene Ontology (GO) annotations for a human gene from org.Hs.eg.db.
   Returns a structured object with fields gene_id, entrez_id, ontology, and an
   'annotations' array. Each annotation entry has go_id, ontology, term, and
   evidence fields. Filter by ontology namespace with the ontology argument.",
  arguments = list(
    gene_id  = ellmer::type_string(
      "Gene symbol (e.g. 'TP53') or Entrez Gene ID (e.g. '7157')."
    ),
    ontology = ellmer::type_string(
      "Ontology filter: 'ALL' (default), 'BP' (Biological Process),
       'MF' (Molecular Function), or 'CC' (Cellular Component)."
    )
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
  list(
    go_id      = r$GOID,
    ontology   = r$ONTOLOGY   %||% NA_character_,
    term       = r$TERM       %||% NA_character_,
    definition = r$DEFINITION %||% NA_character_
  )
}

#' MCP tool: look up GO term name and definition
#' @format An object of class `ToolDef`.
#' @export
tool_go_term_info <- ellmer::tool(
  fn_go_term_info,
  description = "Look up the name, definition, and ontology namespace (BP/MF/CC) for a
   Gene Ontology term by its GO ID. Returns a structured object with fields
   go_id, ontology, term, and definition.",
  arguments = list(
    go_id = ellmer::type_string(
      "A GO term identifier in the format GO:XXXXXXX, e.g. 'GO:0006915'."
    )
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
  term_res  <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db, keys = go_id, columns = "TERM", keytype = "GOID"
  ))
  term_name <- if (nrow(term_res) > 0L) term_res$TERM[[1L]] else NA_character_
  if (nrow(res) == 0L)
    return(list(go_id = go_id, term = term_name, total_genes = 0L, genes = list()))
  res   <- res[order(res$SYMBOL), ]
  total <- nrow(res)
  if (total > max_genes) res <- res[seq_len(max_genes), ]
  genes <- lapply(seq_len(nrow(res)), function(i) list(
    symbol    = if (is.na(res$SYMBOL[[i]]))   NA_character_ else res$SYMBOL[[i]],
    entrez_id = res$ENTREZID[[i]],
    evidence  = if (is.na(res$EVIDENCE[[i]])) NA_character_ else res$EVIDENCE[[i]]
  ))
  list(
    go_id       = go_id,
    term        = term_name,
    total_genes = total,
    genes       = genes
  )
}

#' MCP tool: find human genes annotated to a GO term
#' @format An object of class `ToolDef`.
#' @export
tool_go_to_genes <- ellmer::tool(
  fn_go_to_genes,
  description = "Find human genes annotated to a GO term (including annotations inherited
   through child terms via the transitive closure) using org.Hs.eg.db. Returns
   a structured object with go_id, term, total_genes count, and a 'genes' array.
   Each gene entry has symbol, entrez_id, and evidence fields.",
  arguments = list(
    go_id     = ellmer::type_string(
      "A GO ID, e.g. 'GO:0006915' (apoptotic process)."
    ),
    max_genes = ellmer::type_integer(
      "Maximum number of genes to return. Defaults to 50."
    )
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

  lookup_terms <- function(ids) {
    if (length(ids) == 0L) return(list())
    t <- suppressMessages(AnnotationDbi::select(
      GO.db::GO.db, keys = ids,
      columns = c("GOID", "TERM"), keytype = "GOID"
    ))
    lapply(seq_len(nrow(t)), function(i) list(
      go_id = t$GOID[[i]],
      term  = if (is.na(t$TERM[[i]])) NA_character_ else t$TERM[[i]]
    ))
  }

  go1    <- get_go(gene_id_1)
  go2    <- get_go(gene_id_2)
  shared <- intersect(go1, go2)
  only1  <- setdiff(go1, go2)
  only2  <- setdiff(go2, go1)

  list(
    gene_id_1         = gene_id_1,
    gene_id_2         = gene_id_2,
    ontology          = ont_arg,
    shared            = lookup_terms(shared),
    only_in_gene_id_1 = lookup_terms(only1),
    only_in_gene_id_2 = lookup_terms(only2)
  )
}

#' MCP tool: compare GO annotations between two genes
#' @format An object of class `ToolDef`.
#' @export
tool_compare_gene_go <- ellmer::tool(
  fn_compare_gene_go,
  description = "Compare GO annotations of two human genes to find shared and unique GO
   terms. Returns a structured object with gene_id_1, gene_id_2, ontology, and
   three arrays: 'shared', 'only_in_gene_id_1', and 'only_in_gene_id_2'. Each
   array entry has go_id and term fields. Accepts gene symbols or Entrez IDs.",
  arguments = list(
    gene_id_1 = ellmer::type_string("First gene — symbol or Entrez ID."),
    gene_id_2 = ellmer::type_string("Second gene — symbol or Entrez ID."),
    ontology  = ellmer::type_string(
      "Ontology: 'BP' (default), 'MF', 'CC', or 'ALL'."
    )
  )
)
