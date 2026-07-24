library(testthat)
# Treat a local run like devtools::test() rather than a CRAN check, so the
# KDE-convergence test (skip_on_cran) actually executes here.
if (Sys.getenv("NOT_CRAN") == "") Sys.setenv(NOT_CRAN = "true")

# Dual-mode so BOTH invocations work:
#  * `Rscript tests/testthat.R` from the repo root  -> source R/ directly (no
#    install needed; this is the path documented in the README).
#  * `R CMD check` / `devtools::check()`             -> the package is installed
#    and the working directory is <pkg>.Rcheck/tests/, where R/ does not exist,
#    so run the standard test_check() against the installed namespace.
if (dir.exists("R") && length(list.files("R", pattern = "[.]R$"))) {
  for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
  test_dir("tests/testthat")
} else {
  library(cgrc.bayes)
  test_check("cgrc.bayes")
}
