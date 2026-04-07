test_that("gene_info returns expected fields for TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_info("TP53")
  expect_match(out, "Symbol:.*TP53")
  expect_match(out, "Entrez ID:.*7157")
  expect_match(out, "Chromosome:")
})

test_that("gene_info accepts Entrez ID directly", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_info("7157")
  expect_match(out, "TP53")
})

test_that("gene_info errors on unknown symbol", {
  skip_if_no_orgdb()
  expect_error(geneAnnotationMCP:::fn_gene_info("NOTAREALGENEXYZ"))#, "not found")
})

test_that("symbol_to_entrez converts TP53 correctly", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_symbol_to_entrez("TP53")
  expect_match(out, "TP53 -> 7157")
})

test_that("symbol_to_entrez handles comma-separated list", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_symbol_to_entrez("TP53, BRCA1")
  expect_match(out, "TP53")
  expect_match(out, "BRCA1")
})

test_that("entrez_to_symbol converts 7157 to TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_entrez_to_symbol("7157")
  expect_match(out, "7157 -> TP53")
})
