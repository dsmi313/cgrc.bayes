# Faithful reimplementation of the ORIGINAL Szigeti procedure, for comparison.
#
# Paper's description: split into four strata; fit KDE per stratum (scikit-learn
# KernelDensity, "all parameters left at default value", i.e. Gaussian kernel
# with bandwidth = 1.0); draw samples so the combined pseudo-dataset has the
# target CGR; fit outcome ~ treatment; repeat 100 times; average the estimates
# and p-values.
#
# NOTE ON BANDWIDTH: sklearn's default bandwidth is a FIXED 1.0, not a
# data-adaptive rule. On a 0-100 VAS that is a very narrow kernel; on PANAS it
# is narrow too. Because the estimand only uses stratum MEANS, and KDE
# smoothing is mean-preserving, this turns out not to matter for the point
# estimate - see reports/CHANGELOG.md F4. It would matter for a quantile or
# tail-based estimand.
#
# Averaging p-values across resamples has no accepted inferential
# interpretation. It is reproduced here to characterise the original method,
# not endorsed.

rkde <- function(y, m, bw = 1.0) {
  # sample from a Gaussian KDE: pick a datum, add N(0, bw^2)
  y[sample.int(length(y), m, replace = TRUE)] + stats::rnorm(m, 0, bw)
}

cgr_kde <- function(df, cgr, n_rep = 100, bw = 1.0) {
  st <- cgr_strata(df)
  nk <- cgr_sizes(st, cgr)

  est <- numeric(n_rep); pval <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    ys <- list(); trt <- list()
    for (k in STRATA) {
      m <- nk[[k]]
      if (m <= 0) next
      ys[[k]]  <- rkde(st[[k]], m, bw)
      trt[[k]] <- rep(substr(k, 1, 2) == "AC", m)
    }
    y <- unlist(ys, use.names = FALSE)
    t <- unlist(trt, use.names = FALSE)
    tt <- stats::t.test(y[t], y[!t], var.equal = TRUE)
    est[i]  <- unname(diff(rev(tt$estimate)))
    pval[i] <- tt$p.value
  }
  list(est = mean(est), est_mcse = stats::sd(est) / sqrt(n_rep),
       p = mean(pval), n_rep = n_rep, draws = est)
}

# Separates the three sources of difference between original and Bayesian:
#   (a) Monte Carlo noise from using only 100 resamples
#   (b) KDE vs a Gaussian stratum model
#   (c) averaging p-values vs summarising a posterior
cgr_kde_ladder <- function(df, cgr, reps = c(100, 1000, 10000), bw = 1.0) {
  st <- cgr_strata(df); rat <- cgr_ratios(st)
  analytic <- cgr_delta(cgr, lapply(st, mean), rat$r, rat$s)
  out <- do.call(rbind, lapply(reps, function(m) {
    k <- cgr_kde(df, cgr, n_rep = m, bw = bw)
    data.frame(n_rep = m, kde_est = k$est, kde_mcse = k$est_mcse,
               mean_p = k$p, analytic = analytic,
               diff_from_analytic = k$est - analytic)
  }))
  rownames(out) <- NULL
  out
}
