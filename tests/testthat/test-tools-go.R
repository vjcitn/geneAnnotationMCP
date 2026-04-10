test_that("gene_go_terms returns a data.frame with BP terms for TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_go_terms("TP53", "BP")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("go_id", "ontology", "term", "evidence") %in% names(out)))
  expect_true(nrow(out) > 0L)
  expect_true(all(out$ontology == "BP"))
  expect_true(all(grepl("^GO:[0-9]{7}$", out$go_id)))
})

test_that("gene_go_terms rejects invalid ontology argument", {
  skip_if_no_orgdb()
  expect_error(geneAnnotationMCP:::fn_gene_go_terms("TP53", "XX"), "ontology must be")
})

test_that("go_term_info returns a data.frame with correct fields for GO:0006915", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_go_term_info("GO:0006915")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_true(all(c("go_id", "ontology", "term", "definition") %in% names(out)))
  expect_equal(out$go_id, "GO:0006915")
  expect_equal(out$ontology, "BP")
  expect_true(grepl("apoptotic process", out$term, ignore.case = TRUE))
})

test_that("go_term_info errors on malformed GO ID", {
  expect_error(geneAnnotationMCP:::fn_go_term_info("GO:123"), "not a valid GO ID")
})

test_that("go_to_genes returns a data.frame with genes for GO:0006915", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_go_to_genes("GO:0006915", max_genes = 5L)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("go_id", "symbol", "entrez_id", "evidence") %in% names(out)))
  expect_true(nrow(out) > 0L)
  expect_true(all(out$go_id == "GO:0006915"))
})

test_that("compare_gene_go returns a data.frame for TP53 and MDM2", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_compare_gene_go("TP53", "MDM2", "BP")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("go_id", "term", "ontology", "category") %in% names(out)))
  expect_true("shared" %in% out$category)
  expect_true("TP53" %in% out$category)
  expect_true("MDM2" %in% out$category)
})
