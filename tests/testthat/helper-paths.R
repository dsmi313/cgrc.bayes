# Locate data/pacutes.csv regardless of the working directory the tests run in.
# testthat's test_dir()/test_check() set the working directory to
# tests/testthat/, whereas an interactive run sits at the project root. Without
# this, the real-data reproduction tests silently skip on file.exists().
cgrc_data_path <- function() {
  cands <- c("data/pacutes.csv", "../../data/pacutes.csv")
  for (p in cands) if (file.exists(p)) return(p)
  "data/pacutes.csv"   # fall back to the root-relative path for the skip message
}
