test_that("search_go_terms returns named list with results for apoptosis", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_search_go_terms("apoptosis", "BP", 10L)
  expect_true(is.list(out))
  expect_named(out, c("query", "ontology", "total", "results"))
  expect_equal(out$query, "apoptosis")
  expect_equal(out$ontology, "BP")
  expect_true(out$total > 0)
  first <- out$results[[1]]
  expect_true(grepl("^GO:[0-9]{7}$", first$go_id))
  expect_equal(first$ontology, "BP")
})

test_that("search_go_terms returns empty results for no match", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_search_go_terms("xyzzy_no_such_term_12345", "ALL", 10L)
  expect_true(is.list(out))
  expect_equal(out$total, 0L)
  expect_equal(length(out$results), 0L)
})

test_that("search_go_terms_ols queries OLS4 and returns named list with GO IDs", {
  skip_if_no_orgdb()
  skip_if_offline()
  out <- geneAnnotationMCP:::fn_search_go_terms_ols("programmed cell death", "BP", 5L)
  expect_true(is.list(out))
  expect_named(out, c("query", "ontology", "total", "results"))
  expect_true(length(out$results) > 0)
  expect_true(grepl("^GO:[0-9]{7}$", out$results[[1]]$go_id))
})
