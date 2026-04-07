test_that("gene_go_terms returns BP terms for TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_go_terms("TP53", "BP")
  expect_match(out, "\\[BP\\]")
  expect_match(out, "GO:[0-9]{7}")
})

test_that("gene_go_terms rejects invalid ontology argument", {
  skip_if_no_orgdb()
  expect_error(geneAnnotationMCP:::fn_gene_go_terms("TP53", "XX"), "ontology must be")
})

test_that("go_term_info returns correct fields for GO:0006915", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_go_term_info("GO:0006915")
  expect_match(out, "GO:0006915")
  expect_match(out, "Ontology:.*BP")
  expect_match(out, "apoptotic process", ignore.case = TRUE)
})

test_that("go_term_info errors on malformed GO ID", {
  expect_error(geneAnnotationMCP:::fn_go_term_info("GO:123"), "not a valid GO ID")
})

test_that("go_to_genes returns genes for GO:0006915", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_go_to_genes("GO:0006915", max_genes = 5L)
  expect_match(out, "apoptotic process", ignore.case = TRUE)
  expect_match(out, "Entrez:")
})

test_that("compare_gene_go finds shared BP terms for TP53 and MDM2", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_compare_gene_go("TP53", "MDM2", "BP")
  expect_match(out, "Shared terms")
  expect_match(out, "Only in TP53")
  expect_match(out, "Only in MDM2")
})
