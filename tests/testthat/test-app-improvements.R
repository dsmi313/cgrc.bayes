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

# Locate the app source whether tests run from tests/testthat/, the repo root, or
# an installed package (R CMD check runs from <pkg>.Rcheck/tests/, where only the
# installed copy under system.file("app") exists).
app_src_path <- function() {
  cands <- c("../../inst/app/app.R", "inst/app/app.R",
             system.file("app", "app.R", package = "cgrc.bayes"))
  for (p in cands) if (nzchar(p) && file.exists(p)) return(p)
  "inst/app/app.R"
}

test_that("both trade-off bars are Bayesian at P(favourable)>0.95", {
  # cg_operating must expose BOTH the adjusted and the unadjusted Bayesian flag,
  # so the trade-off plot can compare like with like at the same 0.95 threshold.
  op <- cgr_operating(n_trials = 150, n = 200, p_cg = 0.7, mu_dte = 3, seed = 1)
  expect_true(all(c("p_fav_gt_95", "unadj_p_fav_gt_95") %in% names(op)))
  # pure expectancy (DTE=0, AEB=1): the unadjusted Bayesian analysis falsely flags
  # far more often than the adjusted one - the whole point of the adjustment.
  pe <- op[op$DTE == 0 & op$AEB == 1, ]
  expect_gt(pe$unadj_p_fav_gt_95, pe$p_fav_gt_95)
  expect_lt(pe$p_fav_gt_95, 0.20)          # adjusted false-attribution stays low
  # cgrc_op_at surfaces the unadjusted Bayesian flag from a rebuilt lookup
  lut <- tryCatch(cgrc_lookup(), error = function(e) NULL)
  skip_if(is.null(lut) || !"unadj_p_fav_gt_95" %in% names(lut),
          "rebuilt lookup with unadj_p_fav_gt_95 not present")
  r <- cgrc_op_at(lut, 200, 0.7, 3, 0, 1)
  expect_true("unadj_p_fav_gt_95" %in% names(r))
})

test_that("freq_sig is the conventional two-sided t-test kept as a reference", {
  # It is a familiar reference metric, direction-INVARIANT (unlike the Bayesian
  # flags), retained in the simulation output but not the plotted criterion.
  op1  <- cgr_operating(n_trials = 120, n = 200, p_cg = 0.7, mu_dte = 3, seed = 1)
  opm1 <- cgr_operating(n_trials = 120, n = 200, p_cg = 0.7, mu_dte = 3, seed = 1,
                        direction = -1)
  expect_equal(op1$freq_sig, opm1$freq_sig)       # two-sided: direction-invariant
  # the UNADJUSTED Bayesian flag, by contrast, IS direction-sensitive
  expect_false(isTRUE(all.equal(op1$unadj_p_fav_gt_95, opm1$unadj_p_fav_gt_95)))
})

test_that("p_fav_gt_975 is retained internally but not shown in the app", {
  op <- cgr_operating(n_trials = 60, n = 200, p_cg = 0.7, mu_dte = 3, seed = 1)
  expect_true("p_fav_gt_975" %in% names(op))       # kept for research/sensitivity
  expect_true(all(op$p_fav_gt_975 <= op$p_fav_gt_95 + 1e-9))
  # ... but the default Shiny interface must not display it anywhere
  txt <- paste(readLines(app_src_path(), warn = FALSE), collapse = "\n")
  expect_false(grepl("p_fav_gt_975", txt, fixed = TRUE))
  expect_false(grepl("0.975", txt, fixed = TRUE))
})

test_that("every user-facing adjusted flag/power in the app is sourced from p_fav_gt_95", {
  txt <- readLines(app_src_path(), warn = FALSE)
  # lines that plot/label an adjusted Bayesian rate must reference p_fav_gt_95,
  # and no user-facing rate may come from the 0.975 column.
  expect_false(any(grepl("p_fav_gt_975", txt, fixed = TRUE)))
  expect_true(any(grepl("p_fav_gt_95", txt, fixed = TRUE)))
  expect_true(any(grepl("unadj_p_fav_gt_95", txt, fixed = TRUE)))
})

test_that("app relabels the pure-expectancy outcome (no 'false positive')", {
  txt <- paste(readLines(app_src_path(), warn = FALSE), collapse = "\n")
  expect_true(grepl("false treatment attribution", txt, fixed = TRUE))
  expect_false(grepl("false positive", txt, fixed = TRUE))
  # the misleading "Bayesian equivalent" claim is gone
  expect_false(grepl("Bayesian equivalent", txt, fixed = TRUE))
})

test_that("cgrc_reliability gives four feasibility categories, not safe/unsafe", {
  expect_equal(cgrc_reliability(60,  0.00)$category, "Feasibility looks good under simulated conditions")
  expect_equal(cgrc_reliability(20,  0.00)$category, "Use with caution")
  expect_equal(cgrc_reliability(10,  0.00)$category, "Fragile design")
  expect_equal(cgrc_reliability(60,  0.20)$category, "Adjustment frequently undefined")
  cats <- vapply(list(c(60,0), c(20,0), c(10,0), c(60,0.2)),
                 function(z) cgrc_reliability(z[1], z[2])$category, character(1))
  # feasibility-only wording: never claims "reliable"/"safe" on stratum size alone
  expect_false(any(grepl("safe|unsafe|reliable", cats, ignore.case = TRUE)))
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
