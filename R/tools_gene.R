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
  if (!nzchar(ensembl_ids)) ensembl_ids <- "NA"
  r <- res[1L, ]
  paste(c(
    sprintf("Symbol:      %s", r$SYMBOL    %||% "NA"),
    sprintf("Entrez ID:   %s", r$ENTREZID  %||% "NA"),
    sprintf("Full name:   %s", r$GENENAME  %||% "NA"),
    sprintf("Chromosome:  %s", r$CHR       %||% "NA"),
    sprintf("Cytoband:    %s", r$MAP       %||% "NA"),
    sprintf("Ensembl IDs: %s", ensembl_ids)
  ), collapse = "\n")
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
  "Return basic annotation (symbol, full name, Entrez ID, Ensembl IDs,
   chromosome, cytoband) for a human gene. gene_id can be a gene symbol
   (e.g. 'TP53') or a numeric Entrez ID (e.g. '7157').",
  gene_id = ellmer::type_string(
    "Gene symbol such as 'TP53' or 'BRCA1', or a numeric Entrez ID such as '7157'."
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
  if (nrow(res) == 0L) return("No matches found.")
  paste(
    sprintf("%s -> %s", res$SYMBOL,
            ifelse(is.na(res$ENTREZID), "(not found)", res$ENTREZID)),
    collapse = "\n"
  )
}

#' MCP tool: convert gene symbols to Entrez IDs
#' @format An object of class `ToolDef`.
#' @export
tool_symbol_to_entrez <- ellmer::tool(
  fn_symbol_to_entrez,
  "Convert one or more human gene symbols to NCBI Entrez Gene IDs using
   org.Hs.eg.db. Accepts a comma- or semicolon-separated list of symbols.",
  symbols = ellmer::type_string(
    "One or more gene symbols separated by commas or semicolons, e.g. 'TP53, BRCA1, EGFR'."
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
  if (nrow(res) == 0L) return("No matches found.")
  paste(
    sprintf("%s -> %s (%s)", res$ENTREZID,
            ifelse(is.na(res$SYMBOL),   "(not found)", res$SYMBOL),
            ifelse(is.na(res$GENENAME), "",            res$GENENAME)),
    collapse = "\n"
  )
}

#' MCP tool: convert Entrez IDs to gene symbols
#' @format An object of class `ToolDef`.
#' @export
tool_entrez_to_symbol <- ellmer::tool(
  fn_entrez_to_symbol,
  "Convert one or more NCBI Entrez Gene IDs to official gene symbols and
   full gene names using org.Hs.eg.db.",
  entrez_ids = ellmer::type_string(
    "One or more Entrez Gene IDs separated by commas or semicolons, e.g. '7157, 672, 1956'."
  )
)
