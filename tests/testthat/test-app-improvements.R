# Tests for the app-readability/honesty pass (see inst/app/app.R and R/09_app.R).
# These exercise the acceptance criteria from the improvement brief:
#   - direction-matched frequentist criterion
#   - separate P>0.95 and P>0.975 outputs
#   - pure-expectancy relabelling ("false treatment attribution")
#   - thin-stratum / fragility flagging
#   - arm-specific correct-guess rates
#   - reproducibility under a fixed seed
#   - the analysis-summary and HTML report builders
#   - the app parses/launches without error
# Helper-level tests need no shiny; the server-level test skips when the package
# is not installed (it runs under devtools::test()/R CMD check, which install it).

# Locate inst/app/app.R whether tests run from tests/testthat/ or the repo root.
app_src_path <- function() {
  for (p in c("../../inst/app/app.R", "inst/app/app.R")) if (file.exists(p)) return(p)
  "inst/app/app.R"
}

test_that("cgr_operating uses a DIRECTION-MATCHED frequentist criterion", {
  # A strong real effect in the NEGATIVE direction: highly significant, but not in
  # the +favourable direction. With direction = +1 it must NOT be counted as a
  # favourable frequentist event; with direction = -1 it must be counted.
  op_pos <- cgr_operating(n_trials = 120, n = 300, p_cg = 0.6, mu_dte = -6,
                          mu_aeb = 7.7, seed = 1, direction = 1)
  op_neg <- cgr_operating(n_trials = 120, n = 300, p_cg = 0.6, mu_dte = -6,
                          mu_aeb = 7.7, seed = 1, direction = -1)
  r_pos <- op_pos[op_pos$DTE == 1 & op_pos$AEB == 0, ]
  r_neg <- op_neg[op_neg$DTE == 1 & op_neg$AEB == 0, ]
  expect_lt(r_pos$freq_sig, 0.05)   # wrong-direction significance is excluded
  expect_gt(r_neg$freq_sig, 0.80)   # same trials, favourable direction, counted
})

test_that("both P>0.95 and P>0.975 flags are reported and correctly ordered", {
  op <- cgr_operating(n_trials = 120, n = 200, p_cg = 0.7, mu_dte = 3, seed = 1)
  expect_true(all(c("p_fav_gt_95", "p_fav_gt_975") %in% names(op)))
  # the stricter level can never exceed the looser one
  expect_true(all(op$p_fav_gt_975 <= op$p_fav_gt_95 + 1e-9))
  # cgrc_op_at surfaces both from a lookup (skip if the lookup is not built)
  lut <- tryCatch(cgrc_lookup(), error = function(e) NULL)
  skip_if(is.null(lut), "lookup table not built")
  r <- cgrc_op_at(lut, 200, 0.7, 3, 1, 0)
  expect_true(all(c("p_fav_gt_95", "p_fav_gt_975") %in% names(r)))
})

test_that("cgr_unknown_operating is also direction-matched", {
  op1 <- cgr_unknown_operating(n_trials = 80, n = 300, p_cg = 0.6, u = 0.2,
                               mu_dte = -6, seed = 1, direction = 1)
  opm1 <- cgr_unknown_operating(n_trials = 80, n = 300, p_cg = 0.6, u = 0.2,
                                mu_dte = -6, seed = 1, direction = -1)
  r1  <- op1[op1$DTE == 1 & op1$AEB == 0, ]
  rm1 <- opm1[opm1$DTE == 1 & opm1$AEB == 0, ]
  expect_lt(r1$freq_sig, 0.10)
  expect_gt(rm1$freq_sig, 0.60)
})

test_that("app relabels the pure-expectancy outcome (no 'false positive')", {
  txt <- paste(readLines(app_src_path(), warn = FALSE), collapse = "\n")
  expect_true(grepl("false treatment attribution", txt, fixed = TRUE))
  expect_false(grepl("false positive", txt, fixed = TRUE))
  # the misleading "Bayesian equivalent" claim is gone
  expect_false(grepl("Bayesian equivalent", txt, fixed = TRUE))
})

test_that("cgrc_reliability gives four non-binary categories, not safe/unsafe", {
  expect_equal(cgrc_reliability(60,  0.00)$category, "Reliable under simulated conditions")
  expect_equal(cgrc_reliability(20,  0.00)$category, "Use with caution")
  expect_equal(cgrc_reliability(10,  0.00)$category, "Fragile design")
  expect_equal(cgrc_reliability(60,  0.20)$category, "Adjustment undefined in many simulated trials")
  cats <- vapply(list(c(60,0), c(20,0), c(10,0), c(60,0.2)),
                 function(z) cgrc_reliability(z[1], z[2])$category, character(1))
  expect_false(any(grepl("safe|unsafe", cats, ignore.case = TRUE)))
})

test_that("cgr_expected_strata shows all four strata and sums to n", {
  es <- cgr_expected_strata(120, 0.85)
  expect_named(es, c("ACAC", "ACPL", "PLAC", "PLPL"))
  expect_equal(sum(es), 120)
  expect_equal(min(es), cgr_min_stratum(120, 0.85))
})

test_that("cgr_guess_rates reports correct ARM-SPECIFIC guess rates", {
  # 3 active (2 guess active), 4 placebo (1 guess placebo): overall hides the gap
  trial <- data.frame(
    condition = c("AC","AC","AC","PL","PL","PL","PL"),
    guess     = c("AC","AC","PL","PL","AC","AC","AC"),
    value     = seq_len(7))
  gr <- cgr_guess_rates(trial)
  expect_equal(gr$active,  2/3)
  expect_equal(gr$placebo, 1/4)
  expect_equal(gr$overall, 3/7)
  expect_equal(gr$n_active, 3L); expect_equal(gr$n_placebo, 4L)
  expect_equal(gr$guess_active, 5L); expect_equal(gr$guess_placebo, 2L)
})

test_that("a fixed seed makes cgrc reproducible; a different seed changes it", {
  set.seed(1); d <- sim_aeb(160, p_cg = 0.7, dte_on = TRUE)
  a <- cgrc(d, n_draws = 3000, seed = 42)$summary$post_mean
  b <- cgrc(d, n_draws = 3000, seed = 42)$summary$post_mean
  expect_equal(a, b)                                   # identical under same seed
  cc <- cgrc(d, n_draws = 3000, seed = 99)$summary$post_mean
  expect_false(isTRUE(all.equal(a, cc)))               # different seed differs
})

test_that("analysis summary and HTML report build for an uploaded trial", {
  set.seed(1); d <- sim_aeb(200, p_cg = 0.7, dte_on = TRUE)
  aud   <- cgrc_input_audit(d$condition, d$guess, d$value)
  clean <- aud$clean
  delta <- 0.5 * stats::sd(clean$value)
  f <- list(mode = "binary", trial = clean, dir = 1, audit = aud, delta = delta,
            seed = 1, n_excl_unknown = 0L,
            fit  = cgrc(clean, n_draws = 2000, direction = 1, seed = 1),
            head = cgrc_headline(clean, direction = 1, delta = delta,
                                 n_draws = 2000, seed = 1))
  summ <- cgrc_analysis_summary(f)
  expect_true(all(c("quantity", "value") %in% names(summ)))
  expect_true(any(grepl("counterfactual", summ$quantity)))     # caveat present
  expect_true(any(grepl("Smallest stratum size", summ$quantity)))
  html <- cgrc_build_html_report(f)
  expect_true(grepl("<!doctype html>", html, fixed = TRUE))
  expect_true(grepl("counterfactual", html))                   # not the true effect
  expect_true(grepl("cgrc.bayes", html))                       # package version line
})

test_that("app.R parses without error (launch smoke test)", {
  expect_silent(parse(app_src_path()))
})
