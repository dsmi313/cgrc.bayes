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

test_that("Panel B reproduces the Rmd's PANAS analysis on the real data", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  skip_on_cran()
  app_dir <- system.file("app", package = "cgrc.bayes")
  skip_if(app_dir == "", "app not installed")
  skip_if_not(file.exists(cgrc_data_path()), "pacutes.csv not present")

  raw <- read.csv(cgrc_data_path(), stringsAsFactors = FALSE)
  d <- raw[raw$test_name == "PANAS" & raw$tp == "w1s1", ]
  csv <- data.frame(condition = d$condition, guess = d$guess, value = d$value) # MD/PL
  f <- tempfile(fileext = ".csv"); write.csv(csv, f, row.names = FALSE)

  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "panas.csv"),
                      col_cond = "condition", col_guess = "guess",
                      col_value = "value", direction = "1", rope = 0.1, analyse = 1)
    set.seed(2)
    ff <- safe_fit()
    expect_false(inherits(ff, "cgrc_err"))
    st <- cgr_strata(ff$trial)
    # deterministic quantities must be EXACT
    expect_equal(unname(lengths(st)[c("ACAC","ACPL","PLAC","PLPL")]),
                 c(48L, 43L, 39L, 102L))
    expect_equal(ff$fit$observed_cgr, 0.6466, tolerance = 1e-3)
    z <- cgr_reference_line_test(ff$trial, ff$fit$observed_cgr)
    expect_lt(abs(z$D_at_obs - z$raw_mean_diff), 1e-10)   # no-op identity, exact
    # posterior summaries match the Rmd to Monte Carlo precision
    s <- ff$fit$summary
    expect_equal(s$post_mean[1], 3.157, tolerance = 0.10)  # unadjusted
    expect_equal(s$post_mean[2], 1.080, tolerance = 0.15)  # adjusted at 0.50
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

# ---- UNKNOWN-preserving extension in the app --------------------------------

# a CSV with all six strata, two genuinely-missing rows, and mixed UNKNOWN tokens
unknown_csv <- function() {
  set.seed(42)
  mk <- function(cond, g, k) data.frame(
    condition = cond, guess = g,
    value = rnorm(k, 10 + (cond == "active") * 2, 3), stringsAsFactors = FALSE)
  d <- rbind(
    mk("active",  "active",  26), mk("active",  "placebo", 6),
    mk("active",  "unsure",  12), mk("placebo", "active",  8),
    mk("placebo", "placebo", 22), mk("placebo", "I don't know", 14),
    data.frame(condition = c("active", ""), guess = c("", "placebo"),
               value = c(5, 6), stringsAsFactors = FALSE))          # 2 missing rows
  d <- d[sample(nrow(d)), ]
  f <- tempfile(fileext = ".csv"); write.csv(d, f, row.names = FALSE); f
}

test_that("Panel B preserves UNKNOWN by default, shows six strata, drops nothing silently", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  app_dir <- system.file("app", package = "cgrc.bayes"); skip_if(app_dir == "", "app not installed")
  f <- unknown_csv()
  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "u.csv"))
    session$setInputs(col_cond = "condition", col_guess = "guess", col_value = "value",
                      direction = "1", threshold_mode = "sd", rope = 0.1,
                      seed_b = 1, unknown_level = "UNKNOWN", analyse = 1)
    expect_true(unknown_present())                    # detected
    ff <- safe_fit()
    expect_false(inherits(ff, "cgrc_err"))
    expect_identical(ff$mode, "unknown")              # preserve is the default
    expect_s3_class(ff$ufit, "cgrc_unknown")
    st <- cgr_unknown_strata(ff$trial)
    expect_equal(length(st), 6L)                      # six strata
    expect_true(all(lengths(st)[UNKNOWN_STRATA] > 0)) # all populated
    expect_equal(ff$ufit$n_unknown, 26L)              # 12 + 14 UNKNOWN retained
    # the two missing rows were excluded; the UNKNOWN rows were NOT
    expect_equal(unname(ff$audit$summary[["missing_condition"]]), 1L)
    expect_equal(unname(ff$audit$summary[["missing_guess"]]), 1L)
    expect_equal(unname(ff$audit$summary[["observed_unknown"]]), 26L)
    expect_equal(nrow(ff$trial), 88L)                 # 90 - 2 missing
    # identity holds on the preserved analysis
    z <- cgr_unknown_reference_line_test(ff$trial)
    expect_lt(abs(z$D_at_obs - z$raw_mean_diff), 1e-9)
    # a Markdown report is produced
    rep <- cgrc_build_report(ff)
    expect_true(any(grepl("UNKNOWN-preserving", rep)))
    expect_true(any(grepl("Identity check", rep)))
  })
})

test_that("Panel B complete-case mode excludes UNKNOWN and reports the count", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  app_dir <- system.file("app", package = "cgrc.bayes"); skip_if(app_dir == "", "app not installed")
  f <- unknown_csv()
  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "u.csv"),
                      col_cond = "condition", col_guess = "guess", col_value = "value",
                      direction = "1", threshold_mode = "sd", rope = 0.1, seed_b = 1,
                      unknown_level = "UNKNOWN", unknown_mode = "completecase", analyse = 1)
    ff <- safe_fit()
    expect_false(inherits(ff, "cgrc_err"))
    expect_identical(ff$mode, "binary")
    expect_equal(ff$n_excl_unknown, 26L)              # UNKNOWN explicitly excluded
    expect_true(all(ff$trial$guess %in% c("AC", "PL")))
    expect_equal(nrow(ff$trial), 62L)                 # 88 analysable - 26 UNKNOWN
  })
})

test_that("Panel B does not offer the binary design bridge for an UNKNOWN analysis", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  app_dir <- system.file("app", package = "cgrc.bayes"); skip_if(app_dir == "", "app not installed")
  f <- unknown_csv()
  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "u.csv"),
                      col_cond = "condition", col_guess = "guess", col_value = "value",
                      direction = "1", rope = 0.1, seed_b = 1, unknown_level = "UNKNOWN",
                      analyse = 1)
    session$setInputs(n = 200, pcg = 0.6)
    session$setInputs(do_bridge = 1)                  # must be a no-op in UNKNOWN mode
    expect_identical(safe_fit()$mode, "unknown")
    expect_equal(input$n, 200)                        # design sliders unchanged
    expect_equal(input$pcg, 0.6)
  })
})

test_that("Panel A switches to the UNKNOWN design lookup when u > 0", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  app_dir <- system.file("app", package = "cgrc.bayes"); skip_if(app_dir == "", "app not installed")
  skip_if(is.null(cgrc_unknown_lookup()), "UNKNOWN design lookup not built")
  shiny::testServer(app_dir, {
    session$setInputs(n = 200, pcg = 0.7, eff = 3, mu_aeb = 7.7, u_rate = 0)
    expect_false(grepl("UNKNOWN-preserving", output$verdict$html))   # binary at u = 0
    session$setInputs(u_rate = 0.3)
    expect_match(output$verdict$html, "UNKNOWN-preserving")          # switches at u > 0
    expect_match(output$verdict$html, "UNKNOWN rate 30")
    expect_match(output$feasibility$html, "six strata")
    session$setInputs(n_trials = 30, run_exact = 1)
    op <- exact_rv()
    expect_equal(nrow(op), 4L)
    expect_true("u" %in% names(op))                                  # UNKNOWN model was used
  })
})

test_that("Panel B runs an on-demand UNKNOWN design check (not the binary lookup)", {
  skip_if_not(requireNamespace("shiny", quietly = TRUE), "shiny not installed")
  app_dir <- system.file("app", package = "cgrc.bayes"); skip_if(app_dir == "", "app not installed")
  f <- unknown_csv()
  shiny::testServer(app_dir, {
    session$setInputs(csv = list(datapath = f, name = "u.csv"),
                      col_cond = "condition", col_guess = "guess", col_value = "value",
                      direction = "1", rope = 0.1, seed_b = 1, unknown_level = "UNKNOWN",
                      eff = 3, mu_aeb = 7.7, n_trials = 40, analyse = 1)
    expect_identical(safe_fit()$mode, "unknown")
    session$setInputs(run_unknown_design = 1)
    op <- udesign_rv()
    expect_equal(nrow(op), 4L)                          # 4 DTE x AEB scenarios
    expect_equal(attr(op, "n_trials"), 40)              # uses the trial slider
    expect_true(all(op$coverage95 > 0.8 & op$coverage95 <= 1))
    expect_true("empty_stratum_rate" %in% names(op))
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
