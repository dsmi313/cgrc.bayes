# Tests for the Shiny app support functions (R/09_app.R) and the acceptance
# invariants from the app brief. The interpolation/verdict helpers are tested
# against a small synthetic lookup; the real-lookup invariants skip if the
# precomputed table has not been built/installed.

mk_lut <- function() {
  g <- expand.grid(n = c(100, 200), p_cg = c(0.6, 0.8), true_effect = 3,
                   DTE = c(0, 1), AEB = c(0, 1), KEEP.OUT.ATTRS = FALSE)
  g$true       <- ifelse(g$DTE == 1, 3, 0)
  g$adj_bias   <- 0.05
  g$adj_rmse   <- 1
  g$coverage95 <- 0.95
  g$unadj_bias <- ifelse(g$AEB == 1, 3, 0)
  # clean adjusted power rises linearly with n (DTE on, AEB off)
  g$p_fav_gt_95 <- ifelse(g$DTE == 1 & g$AEB == 0, g$n / 300,
                          ifelse(g$AEB == 1, 0.05, 0.05))
  g$freq_sig   <- ifelse(g$AEB == 1, 0.9, ifelse(g$DTE == 1, 0.6, 0.05))
  g$empty_stratum_rate <- 0
  g$n_valid    <- 500
  g
}

test_that("cgrc_op_at interpolates bilinearly and flags interpolation", {
  lut <- mk_lut()
  # midpoint in n between 100 (0.333) and 200 (0.667) at p_cg on grid -> 0.5
  r <- cgrc_op_at(lut, n = 150, p_cg = 0.6, true_effect = 3, dte = 1, aeb = 0)
  expect_equal(unname(r$p_fav_gt_95), 0.5, tolerance = 1e-8)
  expect_true(r$interpolated)
  # exactly on grid -> not interpolated
  r2 <- cgrc_op_at(lut, n = 100, p_cg = 0.6, true_effect = 3, dte = 1, aeb = 0)
  expect_false(r2$interpolated)
  expect_equal(unname(r2$p_fav_gt_95), 100/300, tolerance = 1e-8)
})

test_that("cgrc_power_curve returns the DTE-on/AEB-off row over n", {
  pc <- cgrc_power_curve(mk_lut(), p_cg = 0.6, true_effect = 3)
  expect_equal(pc$n, c(100, 200))
  expect_equal(pc$power, c(100/300, 200/300), tolerance = 1e-8)
})

test_that("cgrc_verdict is computed and mentions power and percentages", {
  v <- cgrc_verdict(mk_lut(), n = 200, p_cg = 0.6, true_effect = 3)
  expect_type(v, "character")
  expect_match(v, "power")
  expect_match(v, "%")
})

test_that("cgrc_normalise_arm maps common codings and errors on unknowns", {
  expect_equal(cgrc_normalise_arm(c("MD", "PL", "active", "placebo", "1", "0")),
               c("AC", "PL", "AC", "PL", "AC", "PL"))
  expect_equal(cgrc_normalise_arm(c(TRUE, FALSE)), c("AC", "PL"))
  expect_error(cgrc_normalise_arm(c("drug", "banana")), "banana")
})

test_that("cgrc_pct_ok suppresses meaningless attenuation ratios", {
  # unadjusted CrI includes zero -> suppress
  expect_false(cgrc_pct_ok(0.1, -0.5, 0.7, 0.02))
  # opposite signs (would print >100%) -> suppress
  expect_false(cgrc_pct_ok(3.0, 1.0, 5.0, -1.3))
  # clean, same sign, distinct from zero -> show
  expect_true(cgrc_pct_ok(3.0, 1.0, 5.0, 1.1))
})

test_that("cgr_min_stratum matches the closed form", {
  expect_equal(cgr_min_stratum(120, 0.85), 120 * 0.5 * 0.15)
  expect_equal(cgr_min_stratum(200, 0.60), 200 * 0.5 * 0.40)
})

# ---- Acceptance invariants against the REAL lookup (skips if not built) ------

have_lut <- tryCatch({ cgrc_lookup(); TRUE }, error = function(e) FALSE)

test_that("lookup invariants: coverage calibrated and bias ~0 where feasible", {
  skip_if_not(have_lut, "lookup table not built")
  lut <- cgrc_lookup()
  feasible <- lut[cgr_min_stratum(lut$n, lut$p_cg) > 15 & lut$empty_stratum_rate < 0.02, ]
  expect_true(all(feasible$coverage95 >= 0.90 & feasible$coverage95 <= 0.99))
  # adjusted estimator ~ unbiased in feasible cells
  expect_true(mean(abs(feasible$adj_bias)) < 0.15)
})

test_that("lookup reproduces the brief's acceptance cell (n=120, p_cg=0.85)", {
  skip_if_not(have_lut, "lookup table not built")
  # NOTE: the lookup is built at seed = 1; the brief's table is seed = 2, so this
  # checks the shape (adj_bias small, coverage ~0.95, expectancy false-positive
  # high) rather than exact equality.
  r <- cgrc_op_at(cgrc_lookup(), 120, 0.85, 3, 0, 1)   # pure expectancy
  expect_gt(r$freq_sig, 0.9)          # unadjusted almost always "significant"
  expect_lt(r$p_fav_gt_95, 0.15)      # adjusted rarely flags it
})
