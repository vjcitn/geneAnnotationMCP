test_that("search_go_terms returns a data.frame with apoptosis terms", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_search_go_terms("apoptosis", "BP", 10L)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("ontology", "go_id", "term", "definition") %in% names(out)))
  expect_true(nrow(out) > 0L)
  expect_true(all(grepl("^GO:[0-9]{7}$", out$go_id)))
  expect_true(all(out$ontology == "BP"))
})

test_that("search_go_terms returns empty data.frame for no match", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_search_go_terms("xyzzy_no_such_term_12345", "ALL", 10L)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  expect_true(all(c("ontology", "go_id", "term", "definition") %in% names(out)))
})

test_that("search_go_terms_ols queries OLS4 and returns a data.frame with GO IDs", {
  skip_if_no_orgdb()
  skip_if_offline()
  out <- geneAnnotationMCP:::fn_search_go_terms_ols("programmed cell death", "BP", 5L)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("ontology", "go_id", "term", "description") %in% names(out)))
  expect_true(all(grepl("^GO:[0-9]{7}$", out$go_id)))
})
