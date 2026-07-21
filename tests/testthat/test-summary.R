# Tests for the cgr_summary_table bug fixes (CH-13 grid-snap guard, CH-14
# pct_attenuation suppression).

make_cur <- function(obs_row, obs_cgr = 0.65) {
  # minimal curve with a row at 0.5 and a row at obs_cgr
  rbind(
    data.frame(cgr = 0.5,     est = 1.0,        lo = -1.5, hi = 3.6,
               p_fav = 0.79, stringsAsFactors = FALSE),
    data.frame(cgr = obs_cgr, est = obs_row$est, lo = obs_row$lo,
               hi = obs_row$hi, p_fav = 0.5, stringsAsFactors = FALSE))
}

test_that("pct_attenuation is suppressed when the unadjusted CrI includes zero", {
  # unadjusted estimate not distinguishable from zero -> ratio is meaningless
  cur <- make_cur(data.frame(est = -0.011, lo = -0.169, hi = 0.151))
  tab <- cgr_summary_table(cur, obs_cgr = 0.65, "CPS")
  expect_true(is.na(tab$pct_attenuation[2]))
  # abs_attenuation is still reported
  expect_false(is.na(tab$abs_attenuation[2]))
})

test_that("pct_attenuation is reported when the unadjusted CrI excludes zero", {
  cur <- make_cur(data.frame(est = 3.157, lo = 0.70, hi = 5.73))
  tab <- cgr_summary_table(cur, obs_cgr = 0.65, "PANAS")
  expect_false(is.na(tab$pct_attenuation[2]))
})

test_that("cgrc() returns a summarised, printable adjusted analysis", {
  set.seed(1)
  d <- sim_aeb(300, dte_on = TRUE)
  res <- cgrc(d, n_draws = 4000)
  expect_s3_class(res, "cgrc")
  expect_named(res, c("curve", "summary", "observed_cgr"))
  # the curve includes the EXACT observed CGR, so no grid-snapping
  expect_true(any(abs(res$curve$cgr - res$observed_cgr) < 1e-12))
  expect_equal(nrow(res$summary), 2)              # observed + perfect blinding
  expect_output(print(res), "CGRC-adjusted analysis")
})

test_that("cgr_summary_table warns when the grid does not contain the target", {
  cur <- make_cur(data.frame(est = 3.2, lo = 0.7, hi = 5.7), obs_cgr = 0.65)
  # asking for 0.6466 when the grid only has 0.65 must warn (grid-snapping)
  expect_warning(cgr_summary_table(cur, obs_cgr = 0.6466, "PANAS"),
                 "grid-snapping")
  # and no warning when the exact target is present
  cur2 <- make_cur(data.frame(est = 3.16, lo = 0.7, hi = 5.7), obs_cgr = 0.6466)
  expect_silent(cgr_summary_table(cur2, obs_cgr = 0.6466, "PANAS"))
})
