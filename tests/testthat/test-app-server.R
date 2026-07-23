# End-to-end test of the Shiny app's server logic via shiny::testServer (no
# browser). Exercises Panel B: CSV upload -> column mapping -> coding
# normalisation -> cgrc() -> ROPE, and asserts the observed-CGR identity holds
# on the uploaded data (the check the brief requires).

test_that("Panel B analyses an uploaded trial and the identity holds", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  skip_on_cran()
  app_dir <- system.file("app", package = "cgrc.bayes")
  skip_if(app_dir == "", "app not installed")

  set.seed(1)
  d <- sim_aeb(200, p_cg = 0.7, dte_on = TRUE)
  # deliberately non-standard codings, to exercise normalisation
  csv <- data.frame(arm = ifelse(d$condition == "AC", "drug", "placebo"),
                    guess = ifelse(d$guess == "AC", "drug", "placebo"),
                    outcome = d$value)
  f <- tempfile(fileext = ".csv"); write.csv(csv, f, row.names = FALSE)

  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "t.csv"))
    session$setInputs(col_cond = "arm", col_guess = "guess", col_value = "outcome",
                      direction = "1", rope = 0.1, analyse = 1)
    ff <- safe_fit()
    expect_false(inherits(ff, "cgrc_err"))
    expect_equal(nrow(ff$fit$summary), 2)
    # ROPE regions exhaustive
    tot <- ff$rope$p_harm + ff$rope$p_negligible + ff$rope$p_benefit
    expect_true(all(abs(tot - 1) < 1e-9))
    # the no-op identity: curve at observed CGR == raw arm-mean difference
    z <- cgr_reference_line_test(ff$trial, ff$fit$observed_cgr)
    expect_lt(abs(z$D_at_obs - z$raw_mean_diff), 1e-10)
  })
})

test_that("Panel B fails gracefully when an uploaded trial has an empty stratum", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  app_dir <- system.file("app", package = "cgrc.bayes")
  skip_if(app_dir == "", "app not installed")
  # a high-CGR trial where nobody on placebo guessed active -> PLAC is empty and
  # the estimand is undefined. Must be caught, not crash the app.
  csv <- data.frame(
    arm     = c(rep("drug", 20), rep("placebo", 20)),
    guess   = c(rep("drug", 20), rep("placebo", 20)),   # zero discordant guesses
    outcome = rnorm(40))
  f <- tempfile(fileext = ".csv"); write.csv(csv, f, row.names = FALSE)
  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "t.csv"),
                      col_cond = "arm", col_guess = "guess", col_value = "outcome",
                      direction = "1", rope = 0.1, analyse = 1)
    ff <- safe_fit()
    expect_true(inherits(ff, "cgrc_err"))
    expect_match(as.character(ff), "empty stratum")
  })
})

test_that("Panel B errors clearly on an unmappable coding", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  app_dir <- system.file("app", package = "cgrc.bayes")
  skip_if(app_dir == "", "app not installed")
  csv <- data.frame(arm = c("banana", "drug"), guess = c("drug", "placebo"),
                    outcome = c(1, 2))
  f <- tempfile(fileext = ".csv"); write.csv(csv, f, row.names = FALSE)
  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "t.csv"))
    session$setInputs(col_cond = "arm", col_guess = "guess", col_value = "outcome",
                      direction = "1", rope = 0.1, analyse = 1)
    ff <- safe_fit()
    expect_true(inherits(ff, "cgrc_err"))       # caught, not crashed
    expect_match(as.character(ff), "banana")
  })
})
