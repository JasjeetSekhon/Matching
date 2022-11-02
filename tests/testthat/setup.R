is_mac <- function() {
  return(grepl('darwin', R.Version()$platform, ignore.case = TRUE) && grepl('x86', R.Version()$platform, ignore.case = TRUE))
}

skip_if_not_mac <- function() {
  if(is_mac()) {
    return(invisible(TRUE))
  }

  skip('Not run on non-mac platforms due to non-deterministic randomness in implementation of std::uniform_int_distribution across different compilers. Test cases were written only for macOS (clang). We have not had time to write tests for other compilers like gcc yet.')
}
