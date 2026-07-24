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
  expect_named(res, c("curve", "summary", "observed_cgr", "seed"))
  # the curve includes the EXACT observed CGR, so no grid-snapping
  expect_true(any(abs(res$curve$cgr - res$observed_cgr) < 1e-12))
  expect_equal(nrow(res$summary), 2)              # observed + perfect blinding
  expect_output(print(res), "CGRC-adjusted analysis")
})

test_that("cgrc_headline gives two probabilities before and after the correction", {
  set.seed(2)
  # AEB on: the raw signal is inflated by expectancy, so the correction to
  # perfect blinding should PULL the probabilities DOWN.
  d   <- sim_aeb(400, p_cg = 0.85, dte_on = TRUE, aeb_on = TRUE, noise = "all")
  h   <- cgrc_headline(d, direction = 1, n_draws = 6000)
  expect_s3_class(h, "cgrc_headline")
  probs <- c(h$p_dir_obs, h$p_dir_blind, h$p_meaningful_obs, h$p_meaningful_blind)
  expect_true(all(probs >= 0 & probs <= 1))
  # P(direction) is never below P(meaningful): a stricter bar is less probable
  expect_gte(h$p_dir_obs,   h$p_meaningful_obs)
  expect_gte(h$p_dir_blind, h$p_meaningful_blind)
  # expectancy inflation: correcting to perfect blinding lowers both probs
  expect_gte(h$p_dir_obs,        h$p_dir_blind)
  expect_gte(h$p_meaningful_obs, h$p_meaningful_blind)
  # the adjusted CrI brackets the point estimate, and delta is positive
  expect_true(h$adj_lo <= h$adj_est && h$adj_est <= h$adj_hi)
  expect_gt(h$delta, 0)
  # the plain-language line states both probabilities and never a verdict
  expect_output(print(h), "favourable effect")
  expect_output(print(h), "meaningful")
  expect_false(grepl("significant|safe|unsafe", h$text))
})

test_that("cgrc_headline pluralises the delta unit correctly", {
  set.seed(3)
  d <- sim_aeb(200, dte_on = TRUE)
  # delta of exactly 1 reads "1 point"; a non-unit delta reads "points"
  h1 <- cgrc_headline(d, delta = 1,   n_draws = 2000)
  h2 <- cgrc_headline(d, delta = 2.5, n_draws = 2000)
  expect_match(h1$text, "beyond 1 point\\b")
  expect_false(grepl("beyond 1 points", h1$text))
  expect_match(h2$text, "beyond 2.5 points")
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
