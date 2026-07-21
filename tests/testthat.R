library(testthat)
# Treat a local run like devtools::test() rather than a CRAN check, so the
# KDE-convergence test (skip_on_cran) actually executes here.
if (Sys.getenv("NOT_CRAN") == "") Sys.setenv(NOT_CRAN = "true")
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
test_dir("tests/testthat")
