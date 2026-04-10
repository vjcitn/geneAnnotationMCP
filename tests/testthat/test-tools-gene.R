test_that("gene_info returns named list with expected fields for TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_info("TP53")
  expect_true(is.list(out))
  expect_named(out, c("symbol", "entrez_id", "full_name", "chromosome", "cytoband", "ensembl_ids"))
  expect_equal(out$symbol, "TP53")
  expect_equal(out$entrez_id, "7157")
  expect_false(is.na(out$chromosome))
})

test_that("gene_info accepts Entrez ID directly", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_info("7157")
  expect_equal(out$symbol, "TP53")
})

test_that("gene_info errors on unknown symbol", {
  skip_if_no_orgdb()
  expect_error(geneAnnotationMCP:::fn_gene_info("NOTAREALGENEXYZ"))
})

test_that("symbol_to_entrez converts TP53 correctly", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_symbol_to_entrez("TP53")
  expect_true(is.list(out))
  expect_named(out, "mappings")
  expect_equal(out$mappings[[1]]$symbol, "TP53")
  expect_equal(out$mappings[[1]]$entrez_id, "7157")
})

test_that("symbol_to_entrez handles comma-separated list", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_symbol_to_entrez("TP53, BRCA1")
  symbols <- vapply(out$mappings, `[[`, character(1), "symbol")
  expect_true("TP53" %in% symbols)
  expect_true("BRCA1" %in% symbols)
})

test_that("entrez_to_symbol converts 7157 to TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_entrez_to_symbol("7157")
  expect_true(is.list(out))
  expect_named(out, "mappings")
  expect_equal(out$mappings[[1]]$entrez_id, "7157")
  expect_equal(out$mappings[[1]]$symbol, "TP53")
})
