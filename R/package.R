#' geneAnnotationMCP: MCP Tools for Human Gene Annotation and GO Queries
#'
#' Provides nine Model Context Protocol (MCP) tools for interrogating human
#' gene function and Gene Ontology annotation via `org.Hs.eg.db` and `GO.db`.
#'
#' ## Tools
#'
#' | Exported object | What it does |
#' |---|---|
#' | [tool_gene_info] | Symbol, name, chromosome, Ensembl IDs |
#' | [tool_symbol_to_entrez] | Batch symbol -> Entrez ID |
#' | [tool_entrez_to_symbol] | Batch Entrez ID -> symbol + name |
#' | [tool_gene_go_terms] | GO annotations for a gene (filterable by namespace) |
#' | [tool_go_term_info] | Name, definition, namespace for a GO ID |
#' | [tool_go_to_genes] | Genes annotated to a GO term (transitive) |
#' | [tool_compare_gene_go] | Shared vs. unique GO terms between two genes |
#' | [tool_search_go_terms] | Keyword search in local GO.db (offline) |
#' | [tool_search_go_terms_ols] | Fuzzy search via EBI OLS4 REST API |
#'
#' ## Quick start
#'
#' ```r
#' # In your .Rprofile (makes session available to Claude Desktop):
#' mcptools::mcp_session()
#'
#' # In claude_desktop_config.json:
#' # {
#' #   "mcpServers": {
#' #     "gene-annotation": {
#' #       "command": "Rscript",
#' #       "args": ["-e", "geneAnnotationMCP::gene_annotation_mcp_server()"]
#' #     }
#' #   }
#' # }
#' ```
#'
#' @keywords internal
"_PACKAGE"

#' All nine gene annotation MCP tools as a list
#'
#' A named list of all [ellmer::tool()] objects defined in this package,
#' ready to pass to `mcptools::mcp_server(tools = gene_annotation_tools)`.
#'
#' @export
gene_annotation_tools <- list(
  tool_gene_info,
  tool_symbol_to_entrez,
  tool_entrez_to_symbol,
  tool_gene_go_terms,
  tool_go_term_info,
  tool_go_to_genes,
  tool_compare_gene_go,
  tool_search_go_terms,
  tool_search_go_terms_ols
)

#' Start the gene annotation MCP server
#'
#' A thin wrapper around [mcptools::mcp_server()] that pre-loads the full
#' set of gene annotation tools. Blocks the calling process (intended for
#' non-interactive use via `Rscript`).
#'
#' @param ... Additional arguments forwarded to [mcptools::mcp_server()],
#'   e.g. `type = "http"`, `port = 9090`.
#'
#' @return Called for its side effect (blocks indefinitely).
#' @export
#'
#' @examples
#' \dontrun{
#' # From a terminal:
#' # Rscript -e "geneAnnotationMCP::gene_annotation_mcp_server()"
#' gene_annotation_mcp_server()
#' }
gene_annotation_mcp_server <- function(...) {
  mcptools::mcp_server(tools = gene_annotation_tools, ...)
}
