#!/usr/bin/env Rscript
# Standalone entry point for Claude Desktop / Claude Code.
#
# Preferred usage (after package installation):
#   Rscript -e "geneAnnotationMCP::gene_annotation_mcp_server()"
#
# Alternative direct usage (without installing):
#   Rscript inst/scripts/serve.R
library(geneAnnotationMCP)
gene_annotation_mcp_server()
