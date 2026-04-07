# geneAnnotationMCP

An R package providing nine **Model Context Protocol (MCP)** tools for
conversational querying of human gene function and Gene Ontology (GO) annotation,
backed by Bioconductor's `org.Hs.eg.db` and `GO.db`.

## Tools

| Tool | Description |
|---|---|
| `tool_gene_info` | Symbol, full name, Entrez ID, Ensembl IDs, chromosome, cytoband |
| `tool_symbol_to_entrez` | Batch gene symbol → Entrez ID conversion |
| `tool_entrez_to_symbol` | Batch Entrez ID → gene symbol + full name |
| `tool_gene_go_terms` | GO annotations for a gene, filterable by BP / MF / CC |
| `tool_go_term_info` | Name, definition, and namespace for a GO ID |
| `tool_go_to_genes` | Genes annotated to a GO term (transitive, via GOALL) |
| `tool_compare_gene_go` | Shared vs. unique GO terms between two genes |
| `tool_search_go_terms` | Keyword search in local GO.db (offline, cached) |
| `tool_search_go_terms_ols` | Fuzzy search via EBI OLS4 REST API |

## Installation

```r
# Install Bioconductor dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("org.Hs.eg.db", "GO.db", "AnnotationDbi"))

# Install CRAN dependencies
install.packages(c("ellmer", "mcptools", "jsonlite"))

# Install this package
devtools::install("path/to/geneAnnotationMCP")
```

## Claude Desktop setup

Add to `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "gene-annotation": {
      "command": "Rscript",
      "args": ["-e", "geneAnnotationMCP::gene_annotation_mcp_server()"]
    }
  }
}
```

Add to your `~/.Rprofile` to make your interactive RStudio session available:

```r
mcptools::mcp_session()
```

## Development

```r
devtools::document()   # regenerate man/ and NAMESPACE
devtools::test()       # run testthat suite
devtools::check()      # full R CMD check
```
