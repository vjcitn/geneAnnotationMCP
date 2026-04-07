test_that("search_go_terms finds apoptosis terms", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_search_go_terms("apoptosis", "BP", 10L)
  expect_match(out, "GO:[0-9]{7}")
  expect_match(out, "\\[BP\\]")
})

test_that("search_go_terms returns informative message for no match", {
  skip_if_no_orgdb()
  out <- geneAnnotationMCP:::fn_search_go_terms("xyzzy_no_such_term_12345", "ALL", 10L)
  expect_match(out, "No GO terms matched")
})

test_that("search_go_terms_ols queries OLS4 and returns GO IDs", {
  skip_if_no_orgdb()
  skip_if_offline()
  out <- geneAnnotationMCP:::fn_search_go_terms_ols("programmed cell death", "BP", 5L)
  expect_match(out, "GO:[0-9]{7}")
})
