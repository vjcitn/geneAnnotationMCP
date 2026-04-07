# Shared test helpers

skip_if_no_orgdb <- function() {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GO.db")
}

skip_if_offline <- function() {
  testthat::skip_if(
    tryCatch({
      con <- url("https://www.ebi.ac.uk", "r")
      close(con)
      FALSE
    }, error = function(e) TRUE),
    "Network not available"
  )
}
