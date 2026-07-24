# Tests for the UNKNOWN-preserving CGRC extension (R/10_unknown.R, six-stratum
# JAGS in R/03_jags.R, cgrc_normalise_guess in R/09_app.R). The binary method is
# tested unchanged in test-estimand.R / test-rope.R; here we test only the
# extension and, crucially, that it REDUCES to the binary method at u = 0.

# Build a six-stratum trial from a named count vector (labels in UNKNOWN_STRATA).
mk_unknown <- function(counts, seed = 1, arm_shift = 2, sd = 3) {
  set.seed(seed)
  do.call(rbind, Map(function(nm, k) {
    if (k == 0) return(NULL)
    cond <- substr(nm, 1, 2)
    g <- substring(nm, 3); g <- ifelse(g == "U", "UNKNOWN", g)
    data.frame(condition = cond, guess = g,
               value = rnorm(k, 10 + (cond == "AC") * arm_shift, sd),
               stringsAsFactors = FALSE)
  }, names(counts), counts))
}

## ============================ A. Weight algebra =============================
test_that("six weights sum to 1 and split mass by class", {
  for (u in c(0, 0.25, 0.5, 0.8)) for (cc in c(0, 0.3, 0.5686, 1)) {
    for (rr in c(0.1, 0.5, 0.9)) for (ss in c(0.2, 0.7)) for (tt in c(0.3, 0.6)) {
      w <- cgr_unknown_weights(cc, u, rr, ss, tt)
      expect_equal(sum(w), 1, tolerance = 1e-12)
      expect_equal(unname(w[["ACAC"]] + w[["PLPL"]]), (1 - u) * cc, tolerance = 1e-12)
      expect_equal(unname(w[["ACPL"]] + w[["PLAC"]]), (1 - u) * (1 - cc), tolerance = 1e-12)
      expect_equal(unname(w[["ACU"]] + w[["PLU"]]), u, tolerance = 1e-12)
      expect_true(all(w >= 0))
    }
  }
})

test_that("arm denominators are positive on the interior of a full design", {
  for (cc in seq(0.001, 0.999, length.out = 40)) for (u in c(0.1, 0.4)) {
    w <- cgr_unknown_weights(cc, u, 0.4, 0.6, 0.55)
    expect_gt(w[["ACAC"]] + w[["ACPL"]] + w[["ACU"]], 0)
    expect_gt(w[["PLAC"]] + w[["PLPL"]] + w[["PLU"]], 0)
  }
})

## ============================ B. Observed identity ==========================
test_that("Delta(c_obs, u_obs) equals the raw arm-mean difference", {
  set.seed(11)
  for (rep in 1:25) {
    counts <- setNames(sample(3:30, 6, replace = TRUE), UNKNOWN_STRATA)
    d  <- mk_unknown(counts, seed = rep, arm_shift = runif(1, -3, 3))
    st <- cgr_unknown_strata(d); o <- cgr_unknown_observed(st)
    rat <- cgr_unknown_ratios(st)
    mu <- setNames(lapply(UNKNOWN_STRATA,
      function(nm) if (length(st[[nm]])) mean(st[[nm]]) else NA_real_), UNKNOWN_STRATA)
    got <- cgr_unknown_delta(o$c_obs, o$u_obs, mu, rat$r, rat$s, rat$t)
    want <- mean(d$value[d$condition == "AC"]) - mean(d$value[d$condition == "PL"])
    expect_equal(got, want, tolerance = 1e-10)
  }
})

## ======================= C. Reduction to the binary method ==================
test_that("at u = 0 the weights equal the four-stratum weights", {
  for (cc in c(0, 0.3, 0.5, 0.7, 1)) {
    wu <- cgr_unknown_weights(cc, 0, r = 0.4, s = 0.6, t = 0.5)
    wb <- cgr_weights(cc, r = 0.4, s = 0.6)
    expect_equal(unname(wu[c("ACAC", "ACPL", "PLAC", "PLPL")]),
                 unname(wb[c("ACAC", "ACPL", "PLAC", "PLPL")]), tolerance = 1e-12)
    expect_equal(unname(wu[["ACU"]] + wu[["PLU"]]), 0)
  }
})

test_that("at u = 0 the point estimate equals cgr_delta() exactly", {
  set.seed(3)
  d <- sim_aeb(300, 0.7, dte_on = TRUE)                 # binary, no UNKNOWN
  stb <- cgr_strata(d); ratb <- cgr_ratios(stb); mub <- lapply(stb, mean)
  stu <- cgr_unknown_strata(d); ratu <- cgr_unknown_ratios(stu)
  muu <- setNames(lapply(UNKNOWN_STRATA,
    function(nm) if (length(stu[[nm]])) mean(stu[[nm]]) else NA_real_), UNKNOWN_STRATA)
  for (cc in seq(0, 1, length.out = 11)) {
    expect_equal(cgr_unknown_delta(cc, 0, muu, ratu$r, ratu$s, ratu$t),
                 cgr_delta(cc, mub, ratb$r, ratb$s), tolerance = 1e-12)
  }
})

test_that("at u = 0 the posterior summaries agree with cgr_conjugate within MC error", {
  set.seed(5)
  d <- sim_aeb(400, 0.7, dte_on = TRUE)
  grid <- seq(0.1, 0.9, length.out = 9)
  set.seed(9); a <- cgr_conjugate(d, grid, n_draws = 20000)
  set.seed(9); b <- cgr_unknown_conjugate(d, grid, u_target = 0, n_draws = 20000)
  mcse <- sqrt(a$mcse^2 + b$mcse^2)
  expect_true(max(abs(a$est - b$est) / mcse) < 5)
  expect_lt(max(abs(a$p_fav - b$p_fav)), 0.02)
})

## ======================= D. Exact count calculations ========================
test_that("Santana-Penin count table reproduces exactly", {
  d <- mk_unknown(c(ACAC = 26, ACPL = 1, ACU = 12, PLAC = 21, PLPL = 3, PLU = 14))
  o <- cgr_unknown_observed(cgr_unknown_strata(d))
  expect_equal(o$n_total, 77L)
  expect_equal(o$n_unknown, 26L)
  expect_equal(o$n_directional, 51L)
  expect_equal(o$u_obs, 26 / 77, tolerance = 1e-12)
  expect_equal(o$c_obs, 29 / 51, tolerance = 1e-12)      # correct = 26 + 3
  expect_true(all(o$counts > 0))
})

test_that("ketamine count table reproduces exactly", {
  d <- mk_unknown(c(ACAC = 9, ACPL = 5, ACU = 5, PLAC = 8, PLPL = 5, PLU = 6))
  o <- cgr_unknown_observed(cgr_unknown_strata(d))
  expect_equal(o$n_total, 38L)
  expect_equal(o$n_unknown, 11L)
  expect_equal(o$n_directional, 27L)
  expect_equal(o$u_obs, 11 / 38, tolerance = 1e-12)
  expect_equal(o$c_obs, 14 / 27, tolerance = 1e-12)      # correct = 9 + 5
  expect_true(all(o$counts > 0))
})

## ============================ E. Normalisation ==============================
test_that("UNKNOWN synonyms map correctly and casing/whitespace is ignored", {
  x <- c("drug", "placebo", "UNKNOWN", "unsure", " Uncertain ", "not sure",
         "do not know", "don't know", "I don't know", "DK", "idk")
  out <- cgrc_normalise_guess(x, allow_unknown = TRUE)
  expect_equal(out, c("AC", "PL", rep("UNKNOWN", 9)))
})

test_that("a curly apostrophe in don't-know is recognised", {
  expect_equal(cgrc_normalise_guess("I don’t know", allow_unknown = TRUE), "UNKNOWN")
})

test_that("a custom unknown_labels token is recognised", {
  expect_equal(cgrc_normalise_guess(c("weiss nicht", "AC"), allow_unknown = TRUE,
                                    unknown_labels = "weiss nicht"),
               c("UNKNOWN", "AC"))
})

test_that("blank and NA stay missing, not UNKNOWN and not an error", {
  out <- cgrc_normalise_guess(c("", NA, "  ", "AC"), allow_unknown = TRUE)
  expect_equal(out, c(NA, NA, NA, "AC"))
})

test_that("unrecognised values error and name the offending labels", {
  expect_error(cgrc_normalise_guess(c("AC", "banana", "kiwi"), allow_unknown = TRUE),
               "banana")
})

test_that("allow_unknown = FALSE rejects UNKNOWN explicitly", {
  expect_error(cgrc_normalise_guess(c("AC", "unsure"), allow_unknown = FALSE), "unsure")
})

## ======================= F. Sparse and empty cells ==========================
test_that("no UNKNOWN responses reduces to the binary strata and is defined", {
  d <- sim_aeb(200, 0.7, dte_on = TRUE)
  o <- cgr_unknown_observed(cgr_unknown_strata(d))
  expect_equal(o$n_unknown, 0L)
  expect_equal(o$u_obs, 0)
  expect_silent(cgr_unknown_estimable(cgr_unknown_strata(d), 0,
                                      seq(0, 1, length.out = 11), warn = TRUE))
})

test_that("one empty UNKNOWN arm cell gets structurally zero weight (no pseudo-counts)", {
  # PLU empty: t = ACU/(ACU+PLU) = 1, so w_PLU = u*(1-t) = 0 for all u
  d <- mk_unknown(c(ACAC = 20, ACPL = 8, ACU = 10, PLAC = 9, PLPL = 18, PLU = 0))
  st <- cgr_unknown_strata(d)
  expect_equal(length(st[["PLU"]]), 0L)
  rat <- cgr_unknown_ratios(st)
  w <- cgr_unknown_weights(0.5, cgr_unknown_observed(st)$u_obs, rat$r, rat$s, rat$t)
  expect_equal(unname(w[["PLU"]]), 0)
  # the estimate is still finite and the identity still holds
  set.seed(1); fit <- cgr_unknown_conjugate(d, c(0.5, cgr_unknown_observed(st)$c_obs),
                                            n_draws = 4000)
  expect_true(all(is.finite(fit$est)))
})

test_that("an empty correct or incorrect directional class is a clear undefined error", {
  # correct class empty (ACAC = PLPL = 0)
  d1 <- mk_unknown(c(ACAC = 0, ACPL = 10, ACU = 5, PLAC = 10, PLPL = 0, PLU = 5))
  expect_error(cgr_unknown_estimable(cgr_unknown_strata(d1), 0.2, c(0.5), warn = FALSE),
               "correct-guess class")
  # incorrect class empty (ACPL = PLAC = 0), interior c requested
  d2 <- mk_unknown(c(ACAC = 12, ACPL = 0, ACU = 5, PLAC = 0, PLPL = 12, PLU = 5))
  expect_error(cgr_unknown_estimable(cgr_unknown_strata(d2), 0.2, c(0.5), warn = FALSE),
               "incorrect-guess class")
})

test_that("u > 0 with no UNKNOWN observations is a clear undefined error", {
  d <- sim_aeb(120, 0.7, dte_on = TRUE)                  # no UNKNOWN
  expect_error(cgr_unknown_estimable(cgr_unknown_strata(d), 0.3, c(0.5), warn = FALSE),
               "no UNKNOWN responses")
})

test_that("thin cells warn but do not change the estimate", {
  d <- mk_unknown(c(ACAC = 20, ACPL = 2, ACU = 10, PLAC = 12, PLPL = 18, PLU = 9))
  st <- cgr_unknown_strata(d); o <- cgr_unknown_observed(st); rat <- cgr_unknown_ratios(st)
  expect_warning(cgr_unknown_estimable(st, o$u_obs, c(0.5), warn = TRUE), "sparse")
  mu <- setNames(lapply(UNKNOWN_STRATA,
    function(nm) if (length(st[[nm]])) mean(st[[nm]]) else NA_real_), UNKNOWN_STRATA)
  # the point estimate is a deterministic function of the cell means, unaffected
  d_obs <- cgr_unknown_delta(o$c_obs, o$u_obs, mu, rat$r, rat$s, rat$t)
  raw <- mean(d$value[d$condition == "AC"]) - mean(d$value[d$condition == "PL"])
  expect_equal(d_obs, raw, tolerance = 1e-10)
})

## ============================ G. Posterior output ===========================
test_that("posterior credible interval is ordered and p_fav lies in [0,1]", {
  d <- mk_unknown(c(ACAC = 25, ACPL = 12, ACU = 15, PLAC = 14, PLPL = 22, PLU = 13), seed = 4)
  set.seed(2); cur <- cgr_unknown_conjugate(d, seq(0, 1, length.out = 11), n_draws = 6000)
  expect_true(all(cur$hi >= cur$lo))
  expect_true(all(cur$p_fav >= 0 & cur$p_fav <= 1))
  expect_true(all(is.finite(cur$mcse)))
})

test_that("UNKNOWN ROPE regions are exhaustive and sum to 1 at every CGR", {
  d <- mk_unknown(c(ACAC = 25, ACPL = 12, ACU = 15, PLAC = 14, PLPL = 22, PLU = 13), seed = 6)
  set.seed(1); z <- cgr_unknown_rope(d, grid = seq(0, 1, length.out = 11), n_draws = 4000)
  tot <- z$p_harm + z$p_negligible + z$p_benefit
  expect_true(all(abs(tot - 1) < 1e-12))
  expect_true(all(z$p_harm >= 0 & z$p_negligible >= 0 & z$p_benefit >= 0))
})

test_that("UNKNOWN posterior is reproducible under a fixed seed", {
  d <- mk_unknown(c(ACAC = 20, ACPL = 10, ACU = 12, PLAC = 11, PLPL = 18, PLU = 10), seed = 8)
  set.seed(7); a <- cgr_unknown_conjugate(d, seq(0, 1, length.out = 11), n_draws = 3000)
  set.seed(7); b <- cgr_unknown_conjugate(d, seq(0, 1, length.out = 11), n_draws = 3000)
  expect_identical(a, b)
})

test_that("cgrc_unknown grid includes the observed CGR and 0.50 exactly", {
  d <- mk_unknown(c(ACAC = 26, ACPL = 1, ACU = 12, PLAC = 21, PLPL = 3, PLU = 14))
  set.seed(1); fit <- suppressWarnings(cgrc_unknown(d, n_draws = 3000))
  expect_true(any(abs(fit$curve$cgr - 0.5) < 1e-12))
  expect_true(any(abs(fit$curve$cgr - fit$observed_directional_cgr) < 1e-9))
  expect_equal(fit$target_unknown_rate, fit$observed_unknown_rate)
  expect_s3_class(fit, "cgrc_unknown")
})

test_that("direction flips which tail counts as favourable", {
  d <- mk_unknown(c(ACAC = 25, ACPL = 12, ACU = 15, PLAC = 14, PLPL = 22, PLU = 13), seed = 5)
  set.seed(3); up   <- cgr_unknown_conjugate(d, 0.5, n_draws = 6000, direction =  1)
  set.seed(3); down <- cgr_unknown_conjugate(d, 0.5, n_draws = 6000, direction = -1)
  expect_equal(up$p_fav + down$p_fav, 1, tolerance = 0.02)
})

## ================= Optional sensitivities (Section 17E, 17F) =================
test_that("ratio-uncertainty is off by default and widens intervals when on", {
  d <- mk_unknown(c(ACAC = 25, ACPL = 12, ACU = 15, PLAC = 14, PLPL = 22, PLU = 13), seed = 7)
  set.seed(1); base <- cgr_unknown_conjugate(d, seq(0.1, 0.9, length.out = 9), n_draws = 8000)
  set.seed(1); base2 <- cgr_unknown_conjugate(d, seq(0.1, 0.9, length.out = 9), n_draws = 8000,
                                              ratio_uncertainty = FALSE)
  expect_identical(base, base2)                       # default path unchanged
  set.seed(1); ru <- cgr_unknown_conjugate(d, seq(0.1, 0.9, length.out = 9), n_draws = 8000,
                                           ratio_uncertainty = TRUE)
  # propagating ratio uncertainty cannot, on average, narrow the intervals
  expect_gt(mean(ru$hi - ru$lo), mean(base$hi - base$lo) - 1e-6)
  expect_true(all(is.finite(ru$est)))
})

test_that("ratio draws hold a structural (empty-cell) share fixed", {
  d <- mk_unknown(c(ACAC = 20, ACPL = 8, ACU = 10, PLAC = 9, PLPL = 18, PLU = 0))
  rd <- cgr_unknown_ratio_draws(cgr_unknown_strata(d), n_draws = 500)
  expect_length(rd$t, 1L)          # PLU empty -> t fixed at structural 1
  expect_equal(rd$t, 1)
  expect_length(rd$r, 500L)        # correct class has both cells -> randomised
})

test_that("independent shared-guess estimand runs and defaults to the observed marginal", {
  d <- mk_unknown(c(ACAC = 25, ACPL = 12, ACU = 15, PLAC = 14, PLPL = 22, PLU = 13), seed = 3)
  z <- cgr_unknown_independent(d, n_draws = 8000, seed = 1)
  expect_s3_class(z, "cgr_unknown_independent")
  expect_equal(sum(z$q), 1, tolerance = 1e-12)
  expect_true(z$p_favourable >= 0 && z$p_favourable <= 1)
  expect_true(z$hi >= z$lo)
  # seed is recorded and reproducible
  z2 <- cgr_unknown_independent(d, n_draws = 8000, seed = 1)
  expect_equal(z$est, z2$est, tolerance = 1e-12)
})

test_that("independent estimand errors when a weighted class has an empty arm cell", {
  d <- mk_unknown(c(ACAC = 20, ACPL = 10, ACU = 0, PLAC = 12, PLPL = 18, PLU = 8))
  # default q gives the UNKNOWN class positive weight (PLU nonempty) but ACU empty
  expect_error(cgr_unknown_independent(d, n_draws = 2000), "empty arm cell")
})

## ============================ H. Backend agreement ==========================
test_that("conjugate and normal-JAGS UNKNOWN backends agree within MC error", {
  skip_if_not(requireNamespace("rjags", quietly = TRUE), "rjags/JAGS not installed")
  d <- mk_unknown(c(ACAC = 40, ACPL = 22, ACU = 20, PLAC = 25, PLPL = 35, PLU = 18), seed = 7)
  set.seed(1)
  res <- cgr_unknown_check_backends(d, grid = seq(0, 1, length.out = 11),
                                    n_draws = 30000,
                                    jags_args = list(n_iter = 6000, n_burn = 1500, n_chains = 4))
  expect_lt(res$max_abs_z, 5)
  expect_lt(res$max_rhat, 1.05)
  expect_gt(res$min_ess, 1000)
  expect_match(res$verdict, "PASS")
})

test_that("the hierarchical (partial-pooling) UNKNOWN sensitivity runs", {
  skip_if_not(requireNamespace("rjags", quietly = TRUE), "rjags/JAGS not installed")
  d <- mk_unknown(c(ACAC = 30, ACPL = 16, ACU = 14, PLAC = 18, PLPL = 26, PLU = 12), seed = 11)
  b <- cgr_unknown_jags(d, grid = c(0.3, 0.5, 0.7), pooling = "partial",
                        n_iter = 4000, n_burn = 1000, n_chains = 2)
  expect_true(all(is.finite(b$est)))
  expect_match(b$method[1], "pooled")
  expect_identical(attr(b, "pooling"), "partial")
  expect_true(all(is.finite(b$rhat)))
})

test_that("the Student-t UNKNOWN path runs and returns diagnostics", {
  skip_if_not(requireNamespace("rjags", quietly = TRUE), "rjags/JAGS not installed")
  d <- mk_unknown(c(ACAC = 30, ACPL = 16, ACU = 14, PLAC = 18, PLPL = 26, PLU = 12), seed = 9)
  b <- cgr_unknown_jags(d, grid = c(0.3, 0.5, 0.7), likelihood = "t",
                        n_iter = 3000, n_burn = 1000, n_chains = 2)
  expect_true(all(is.finite(b$est)))
  expect_true(all(is.finite(b$ess)))
  expect_true(is.finite(attr(b, "nu")))
})
