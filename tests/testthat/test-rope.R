# Tests for the two additions ported from the self-contained demonstration:
# the ROPE decomposition and the reference-line diagnostic.

test_that("cgr_rope regions are exhaustive and sum to 1 at every CGR", {
  set.seed(1)
  d <- sim_aeb(300, dte_on = TRUE)
  z <- cgr_rope(d, grid = seq(0, 1, length.out = 11), n_draws = 4000)
  tot <- z$p_harm + z$p_negligible + z$p_benefit
  expect_true(all(abs(tot - 1) < 1e-12))
  expect_true(all(z$p_harm >= 0 & z$p_benefit >= 0 & z$p_negligible >= 0))
})

test_that("cgr_rope delta defaults to 0.1 * pooled SD", {
  set.seed(2)
  d <- sim_aeb(200)
  z <- cgr_rope(d, grid = 0.5, n_draws = 2000)
  expect_equal(z$delta, 0.1 * sd(d$value), tolerance = 1e-12)
  z2 <- cgr_rope(d, grid = 0.5, n_draws = 2000, delta = 3)
  expect_equal(z2$delta, 3, tolerance = 1e-12)
})

test_that("wider ROPE never lowers P(negligible)", {
  set.seed(3)
  d <- sim_aeb(300, dte_on = TRUE)
  s <- cgr_rope_sensitivity(d, at_cgr = 0.5,
                            fracs = c(0.05, 0.1, 0.2, 0.4), n_draws = 6000)
  expect_true(all(diff(s$p_negligible) >= -0.02))   # monotone up to MC noise
  expect_true(all(s$delta_in_SD == c(0.05, 0.1, 0.2, 0.4)))
})

test_that("direction flips harm and benefit but not the negligible mass", {
  set.seed(4)
  d <- sim_aeb(300, dte_on = TRUE)
  up   <- cgr_rope(d, grid = 0.5, n_draws = 8000, direction =  1)
  down <- cgr_rope(d, grid = 0.5, n_draws = 8000, direction = -1)
  expect_equal(up$p_benefit, down$p_harm,    tolerance = 0.02)
  expect_equal(up$p_harm,    down$p_benefit, tolerance = 0.02)
  expect_equal(up$p_negligible, down$p_negligible, tolerance = 0.02)
})

test_that("reference-line test: at the observed CGR it reduces to the identity", {
  set.seed(5)
  d <- sim_aeb(400, dte_on = TRUE, aeb_on = TRUE)
  o <- cgr_observed(cgr_strata(d))
  z <- cgr_reference_line_test(d, orig_cgr = o, published_unadj = 5)
  # at c = c_obs the curve equals the raw mean difference (the identity)
  expect_equal(z$D_at_obs, z$raw_mean_diff, tolerance = 1e-10)
  # and evaluating at orig_cgr = observed CGR must match D_at_obs exactly
  expect_equal(z$D_at_orig_cgr, z$D_at_obs, tolerance = 1e-10)
  expect_equal(z$err_at_orig, z$err_at_obs, tolerance = 1e-10)
})

test_that("reference-line test flags a misplaced line on the real PANAS data", {
  skip_if_not(file.exists(cgrc_data_path()))
  raw <- read.csv(cgrc_data_path(), stringsAsFactors = FALSE)
  x <- raw[raw$test_name == "PANAS" & raw$tp == "w1s1", ]
  d <- data.frame(condition = ifelse(x$condition == "MD", "AC", "PL"),
                  guess = ifelse(x$guess == "MD", "AC", "PL"),
                  value = x$value, stringsAsFactors = FALSE)
  z <- cgr_reference_line_test(d, orig_cgr = 0.72, published_unadj = 3.2)
  expect_lt(abs(z$err_at_obs), 0.15)     # lands on the published unadjusted
  expect_gt(abs(z$err_at_orig), 0.8)     # overshoots badly at 0.72
})
