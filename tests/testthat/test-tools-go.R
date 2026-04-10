test_that("gene_go_terms returns named list with BP annotations for TP53", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_gene_go_terms("TP53", "BP")
  expect_true(is.list(out))
  expect_named(out, c("gene_id", "entrez_id", "ontology", "annotations"))
  expect_equal(out$gene_id, "TP53")
  expect_equal(out$ontology, "BP")
  expect_true(length(out$annotations) > 0)
  first <- out$annotations[[1]]
  expect_true(grepl("^GO:[0-9]{7}$", first$go_id))
  expect_equal(first$ontology, "BP")
})

test_that("gene_go_terms rejects invalid ontology argument", {
  skip_if_no_orgdb()
  expect_error(geneAnnotationMCP:::fn_gene_go_terms("TP53", "XX"), "ontology must be")
})

test_that("go_term_info returns named list with correct fields for GO:0006915", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_go_term_info("GO:0006915")
  expect_true(is.list(out))
  expect_named(out, c("go_id", "ontology", "term", "definition"))
  expect_equal(out$go_id, "GO:0006915")
  expect_equal(out$ontology, "BP")
  expect_true(grepl("apoptotic process", out$term, ignore.case = TRUE))
})

test_that("go_term_info errors on malformed GO ID", {
  expect_error(geneAnnotationMCP:::fn_go_term_info("GO:123"), "not a valid GO ID")
})

test_that("go_to_genes returns named list for GO:0006915", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_go_to_genes("GO:0006915", max_genes = 5L)
  expect_true(is.list(out))
  expect_named(out, c("go_id", "term", "total_genes", "genes"))
  expect_equal(out$go_id, "GO:0006915")
  expect_true(grepl("apoptotic process", out$term, ignore.case = TRUE))
  expect_true(out$total_genes > 0)
  expect_true(length(out$genes) <= 5L)
  first <- out$genes[[1]]
  expect_true(!is.null(first$symbol))
  expect_true(!is.null(first$entrez_id))
})

test_that("compare_gene_go returns named list with shared/unique arrays for TP53 and MDM2", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_compare_gene_go("TP53", "MDM2", "BP")
  expect_true(is.list(out))
  expect_named(out, c("gene_id_1", "gene_id_2", "ontology",
                       "shared", "only_in_gene_id_1", "only_in_gene_id_2"))
  expect_equal(out$gene_id_1, "TP53")
  expect_equal(out$gene_id_2, "MDM2")
  expect_true(is.list(out$shared))
  expect_true(length(out$shared) > 0)
})
