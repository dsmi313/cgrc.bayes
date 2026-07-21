# Requirement 7: automated tests for weights, denominators, endpoints,
# invalid codes, empty strata, tiny strata, and seed reproducibility.

test_that("weights sum to 1 and split c / 1-c by guess class", {
  for (cc in c(0, 0.25, 0.5, 0.646552, 1)) {
    for (rr in c(0.1, 0.5, 0.9)) for (ss in c(0.1, 0.5, 0.9)) {
      w <- cgr_weights(cc, rr, ss)
      expect_equal(sum(w), 1, tolerance = 1e-12)
      expect_equal(unname(w[["ACAC"]] + w[["PLPL"]]), cc, tolerance = 1e-12)
      expect_equal(unname(w[["ACPL"]] + w[["PLAC"]]), 1 - cc, tolerance = 1e-12)
    }
  }
})

test_that("treatment-arm denominators stay positive on the interior", {
  for (cc in seq(0.001, 0.999, length.out = 50)) {
    w <- cgr_weights(cc, 0.68, 0.524)
    expect_gt(w[["ACAC"]] + w[["ACPL"]], 0)
    expect_gt(w[["PLPL"]] + w[["PLAC"]], 0)
  }
})

test_that("endpoints reduce to single-stratum contrasts", {
  mu <- list(ACAC = 12, ACPL = 10, PLAC = 11, PLPL = 9)
  expect_equal(cgr_delta(1, mu, 0.4, 0.6), mu$ACAC - mu$PLPL)
  expect_equal(cgr_delta(0, mu, 0.4, 0.6), mu$ACPL - mu$PLAC)
})

test_that("near-endpoint behaviour is continuous, not explosive", {
  mu <- list(ACAC = 12, ACPL = 10, PLAC = 11, PLPL = 9)
  near0 <- cgr_delta(1e-8, mu, 0.4, 0.6); near1 <- cgr_delta(1 - 1e-8, mu, 0.4, 0.6)
  expect_lt(abs(near0 - (mu$ACPL - mu$PLAC)), 1e-4)
  expect_lt(abs(near1 - (mu$ACAC - mu$PLPL)), 1e-4)
})

test_that("degenerate r or s errors instead of returning NaN", {
  mu <- list(ACAC = 12, ACPL = 10, PLAC = 11, PLPL = 9)
  expect_error(cgr_delta(1, mu, r = 1, s = 0.5), "degenerate")
})

test_that("invalid treatment/guess codes are rejected", {
  bad <- data.frame(condition = c("AC", "banana"), guess = c("AC", "PL"),
                    value = c(1, 2))
  expect_error(cgr_strata(bad), "AC/PL")
})

test_that("empty strata error rather than silently dropping", {
  d <- data.frame(condition = c("AC", "AC", "PL"),
                  guess = c("AC", "AC", "PL"), value = 1:3)
  expect_error(cgr_strata(d), "empty stratum")
})

test_that("very small strata still produce a finite, wider posterior", {
  set.seed(99)
  d <- data.frame(
    condition = c(rep("AC", 6), rep("PL", 6)),
    guess     = c(rep("AC", 3), rep("PL", 3), rep("AC", 3), rep("PL", 3)),
    value     = rnorm(12, 10, 3))
  z <- cgr_conjugate(d, seq(0.1, 0.9, length.out = 9), n_draws = 2000)
  expect_true(all(is.finite(z$est)))
  expect_true(all(z$hi > z$lo))
  big <- cgr_conjugate(sim_aeb(400), seq(0.1, 0.9, length.out = 9),
                       n_draws = 2000)
  expect_gt(mean(z$hi - z$lo), mean(big$hi - big$lo))
})

test_that("results are reproducible under a fixed seed", {
  d <- sim_aeb(200)
  set.seed(7); a <- cgr_conjugate(d, seq(0, 1, length.out = 11), n_draws = 3000)
  set.seed(7); b <- cgr_conjugate(d, seq(0, 1, length.out = 11), n_draws = 3000)
  expect_identical(a, b)
})

test_that("supplementary stratum allocation reproduces 54/58/42/46", {
  ex <- data.frame(
    condition = c(rep("PL", 65), rep("AC", 35), rep("PL", 25), rep("AC", 75)),
    guess     = c(rep("PL", 65), rep("PL", 35), rep("AC", 25), rep("AC", 75)),
    value = 0)
  got <- cgr_sizes(cgr_strata(ex), 0.5)
  expect_identical(got[c("ACAC","ACPL","PLAC","PLPL")],
                   c(ACAC = 54, ACPL = 58, PLAC = 42, PLPL = 46))
})

test_that("observed-CGR identity holds across many configurations", {
  set.seed(3)
  for (p in c(0.5, 0.7, 0.9)) for (n in c(120, 400)) {
    for (z in list(c(TRUE,TRUE), c(FALSE,TRUE), c(TRUE,FALSE))) {
      d  <- sim_aeb(n, p, z[1], z[2])
      st <- cgr_strata(d); rat <- cgr_ratios(st)
      got  <- cgr_delta(cgr_observed(st), lapply(st, mean), rat$r, rat$s)
      want <- mean(d$value[d$condition == "AC"]) -
              mean(d$value[d$condition == "PL"])
      expect_equal(got, want, tolerance = 1e-10)
    }
  }
})

test_that("nig_draws can return paired mu and sigma2", {
  set.seed(1)
  y <- rnorm(50, 10, 3)
  o <- nig_draws(y, n_draws = 500, return_sigma2 = TRUE)
  expect_named(o, c("mu", "sigma2", "hyper"))
  expect_length(o$mu, 500); expect_length(o$sigma2, 500)
  expect_true(all(o$sigma2 > 0))
  # with a vague prior the posterior mean sits on the sample mean
  expect_equal(unname(o$hyper[["mn"]]), mean(y), tolerance = 1e-5)
  expect_equal(unname(o$hyper[["ss"]]), sum((y - mean(y))^2), tolerance = 1e-8)
})

test_that("mu and sigma2 are dependent WITHIN a draw", {
  # larger sigma2 must pair with more dispersed mu; that dependence is the
  # mechanism by which variance uncertainty reaches the mean
  set.seed(2)
  o <- nig_draws(rnorm(40, 10, 5), n_draws = 20000, return_sigma2 = TRUE)
  dev <- abs(o$mu - o$hyper[["mn"]])
  hi <- o$sigma2 > median(o$sigma2)
  expect_gt(mean(dev[hi]), mean(dev[!hi]))
})

test_that("direction argument flips which tail counts as favourable", {
  set.seed(4)
  d <- sim_aeb(300, dte_on = TRUE)
  up   <- cgr_conjugate(d, 0.5, n_draws = 4000, direction =  1)
  down <- cgr_conjugate(d, 0.5, n_draws = 4000, direction = -1)
  expect_equal(up$p_fav + down$p_fav, 1, tolerance = 0.02)
})
