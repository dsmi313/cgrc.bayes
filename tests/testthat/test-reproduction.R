# Published-target regression tests. Values verified against an independent
# reimplementation. Real-data tests need data/pacutes.csv present.

test_that("Table 1: adjustment recovers the true DTE in all four cells", {
  set.seed(101)
  for (z in list(c(0,0,0), c(1,0,3), c(0,1,0), c(1,1,3))) {
    d  <- sim_aeb(20000, 0.7, as.logical(z[1]), as.logical(z[2]), noise = "all")
    st <- cgr_strata(d); rt <- cgr_ratios(st); mu <- lapply(st, mean)
    adj   <- cgr_delta(0.5, mu, rt$r, rt$s)
    unadj <- cgr_delta(cgr_observed(st), mu, rt$r, rt$s)
    expect_equal(adj, z[3], tolerance = 0.35)
    if (z[2] == 1) expect_gt(unadj - z[3], 1.5) else
      expect_equal(unadj, z[3], tolerance = 0.35)
  }
})

test_that("Table 2: real data reproduces the published effects", {
  skip_if_not(file.exists(cgrc_data_path()))
  raw <- read.csv(cgrc_data_path(), stringsAsFactors = FALSE)
  tg <- list(PANAS = c(3.2, 1.1, 232), mood = c(6.4, 2.7, 232),
             energy = c(11.5, 6.8, 232), CPS = c(0.0, 0.0, 186))
  for (sc in names(tg)) {
    x <- raw[raw$test_name == sc & raw$tp == "w1s1", ]
    d <- data.frame(condition = ifelse(x$condition == "MD", "AC", "PL"),
                    guess = ifelse(x$guess == "MD", "AC", "PL"),
                    value = x$value, stringsAsFactors = FALSE)
    expect_equal(nrow(d), tg[[sc]][3], info = sc)
    st <- cgr_strata(d); rt <- cgr_ratios(st); mu <- lapply(st, mean)
    obs <- cgr_observed(st)
    expect_equal(cgr_delta(obs, mu, rt$r, rt$s), tg[[sc]][1], tolerance = 0.15)
    expect_equal(cgr_delta(0.5, mu, rt$r, rt$s), tg[[sc]][2], tolerance = 0.35)
  }
})

test_that("the observed CGR is 0.647, NOT the quoted 0.72", {
  skip_if_not(file.exists(cgrc_data_path()))
  raw <- read.csv(cgrc_data_path(), stringsAsFactors = FALSE)
  x <- raw[raw$test_name == "PANAS" & raw$tp == "w1s1", ]
  d <- data.frame(condition = ifelse(x$condition == "MD", "AC", "PL"),
                  guess = ifelse(x$guess == "MD", "AC", "PL"),
                  value = x$value, stringsAsFactors = FALSE)
  obs <- cgr_observed(cgr_strata(d))
  expect_equal(obs, 0.6466, tolerance = 0.001)
  expect_false(isTRUE(all.equal(obs, 0.72, tolerance = 0.005)))

  # the 0.72 hypothesis: correct-guess rate WITHIN the placebo arm
  pl <- d[d$condition == "PL", ]
  expect_equal(mean(pl$condition == pl$guess), 0.7234, tolerance = 0.001)
})

test_that("KDE resampling converges to the analytic value", {
  skip_on_cran()
  set.seed(12)
  d <- sim_aeb(400, dte_on = TRUE, aeb_on = TRUE)
  st <- cgr_strata(d); rt <- cgr_ratios(st)
  analytic <- cgr_delta(0.5, lapply(st, mean), rt$r, rt$s)
  k <- cgr_kde(d, 0.5, n_rep = 4000)
  expect_lt(abs(k$est - analytic), 4 * k$est_mcse + 0.05)
})
