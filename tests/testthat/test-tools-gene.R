test_that("gene_info returns a data.frame with expected fields for TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_info("TP53")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_true(all(c("symbol", "entrez_id", "full_name", "chromosome",
                    "cytoband", "ensembl_ids") %in% names(out)))
  expect_equal(out$symbol, "TP53")
  expect_equal(out$entrez_id, "7157")
  expect_false(is.na(out$chromosome))
})

test_that("gene_info accepts Entrez ID directly", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_info("7157")
  expect_s3_class(out, "data.frame")
  expect_equal(out$symbol, "TP53")
})

test_that("gene_info errors on unknown symbol", {
  skip_if_no_orgdb()
  expect_error(geneAnnotationMCP:::fn_gene_info("NOTAREALGENEXYZ"))
})

test_that("symbol_to_entrez returns a data.frame with correct columns", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_symbol_to_entrez("TP53")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("symbol", "entrez_id") %in% names(out)))
  expect_equal(out$symbol[out$symbol == "TP53"], "TP53")
  expect_equal(out$entrez_id[out$symbol == "TP53"], "7157")
})

test_that("symbol_to_entrez handles comma-separated list", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_symbol_to_entrez("TP53, BRCA1")
  expect_s3_class(out, "data.frame")
  expect_true("TP53" %in% out$symbol)
  expect_true("BRCA1" %in% out$symbol)
})

test_that("entrez_to_symbol returns a data.frame with correct columns", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_entrez_to_symbol("7157")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("entrez_id", "symbol", "gene_name") %in% names(out)))
  expect_equal(out$entrez_id, "7157")
  expect_equal(out$symbol, "TP53")
})
