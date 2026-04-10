## MCP tools: gene_info, symbol_to_entrez, entrez_to_symbol

# ── gene_info ─────────────────────────────────────────────────────────────────

fn_gene_info <- function(gene_id) {
  entrez <- .lookup_entrez(gene_id)
  avail  <- AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db)
  cols   <- intersect(c("SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL", "CHR", "MAP"), avail)
  res    <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = entrez,
    columns = cols,
    keytype = "ENTREZID"
  ))
  if (nrow(res) == 0L) stop(sprintf("No annotation found for Entrez ID %s.", entrez))
  ensembl_ids <- paste(unique(stats::na.omit(res$ENSEMBL)), collapse = ", ")
  if (!nzchar(ensembl_ids)) ensembl_ids <- NA_character_
  r <- res[1L, ]
  list(
    symbol      = r$SYMBOL   %||% NA_character_,
    entrez_id   = r$ENTREZID %||% NA_character_,
    full_name   = r$GENENAME %||% NA_character_,
    chromosome  = r$CHR      %||% NA_character_,
    cytoband    = r$MAP      %||% NA_character_,
    ensembl_ids = ensembl_ids
  )
}

#' MCP tool: look up basic annotation for a human gene
#'
#' Returns a pre-built [ellmer::tool()] object. Pass this (or the full list
#' from [gene_annotation_tools]) to [gene_annotation_mcp_server()] or
#' directly to `mcptools::mcp_server(tools = ...)`.
#'
#' @format An object of class `ToolDef` (from the ellmer package).
#' @export
tool_gene_info <- ellmer::tool(
  fn_gene_info,
  description = "Return basic annotation for a human gene. gene_id can be a gene symbol
   (e.g. 'TP53') or a numeric Entrez ID (e.g. '7157'). Returns a structured
   object with fields: symbol, entrez_id, full_name, chromosome, cytoband,
   ensembl_ids.",
  arguments = list(
    gene_id = ellmer::type_string(
      "Gene symbol such as 'TP53' or 'BRCA1', or a numeric Entrez ID such as '7157'."
    )
  )
)

# ── symbol_to_entrez ──────────────────────────────────────────────────────────

fn_symbol_to_entrez <- function(symbols) {
  sv  <- unique(trimws(strsplit(symbols, "[,;[:space:]]+")[[1L]]))
  sv  <- sv[nchar(sv) > 0L]
  res <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = toupper(sv),
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "SYMBOL"
  ))
  if (nrow(res) == 0L)
    return(list(mappings = list()))
  list(
    mappings = lapply(seq_len(nrow(res)), function(i) list(
      symbol    = res$SYMBOL[[i]],
      entrez_id = if (is.na(res$ENTREZID[[i]])) NA_character_ else res$ENTREZID[[i]]
    ))
  )
}

#' MCP tool: convert gene symbols to Entrez IDs
#' @format An object of class `ToolDef`.
#' @export
tool_symbol_to_entrez <- ellmer::tool(
  fn_symbol_to_entrez,
  description = "Convert one or more human gene symbols to NCBI Entrez Gene IDs using
   org.Hs.eg.db. Returns a structured object with a 'mappings' array; each
   entry has fields 'symbol' and 'entrez_id'. Accepts a comma- or
   semicolon-separated list of symbols.",
  arguments = list(
    symbols = ellmer::type_string(
      "One or more gene symbols separated by commas or semicolons, e.g. 'TP53, BRCA1, EGFR'."
    )
  )
)

# ── entrez_to_symbol ──────────────────────────────────────────────────────────

fn_entrez_to_symbol <- function(entrez_ids) {
  iv  <- unique(trimws(strsplit(entrez_ids, "[,;[:space:]]+")[[1L]]))
  iv  <- iv[nchar(iv) > 0L]
  res <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = iv,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
  ))
  if (nrow(res) == 0L)
    return(list(mappings = list()))
  list(
    mappings = lapply(seq_len(nrow(res)), function(i) list(
      entrez_id = res$ENTREZID[[i]],
      symbol    = if (is.na(res$SYMBOL[[i]]))   NA_character_ else res$SYMBOL[[i]],
      gene_name = if (is.na(res$GENENAME[[i]])) NA_character_ else res$GENENAME[[i]]
    ))
  )
}

#' MCP tool: convert Entrez IDs to gene symbols
#' @format An object of class `ToolDef`.
#' @export
tool_entrez_to_symbol <- ellmer::tool(
  fn_entrez_to_symbol,
  description = "Convert one or more NCBI Entrez Gene IDs to official gene symbols and
   full gene names using org.Hs.eg.db. Returns a structured object with a
   'mappings' array; each entry has fields 'entrez_id', 'symbol', and
   'gene_name'.",
  arguments = list(
    entrez_ids = ellmer::type_string(
      "One or more Entrez Gene IDs separated by commas or semicolons, e.g. '7157, 672, 1956'."
    )
  )
)
