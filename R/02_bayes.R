# Normal-Inverse-Gamma posterior for each stratum mean.
#
#   y_ij | mu_j, sigma2_j ~ N(mu_j, sigma2_j)
#   mu_j | sigma2_j       ~ N(m0, sigma2_j / k0)      <- NOTE: depends on sigma2
#   sigma2_j              ~ InverseGamma(a0, b0)
#
# Posterior update:
#   k_n = k0 + n
#   m_n = (k0*m0 + n*ybar) / k_n
#   a_n = a0 + n/2
#   b_n = b0 + SS/2 + k0*n*(ybar - m0)^2 / (2*k_n)
#
# k0 is a prior sample size in pseudo-observations: with k0 = 1e-6 against
# n = 48 real observations, m_n differs from ybar by ~4e-07.
#
# We draw sigma2 FIRST and then mu | sigma2 because that is the factorisation
# of the joint posterior: p(mu, sigma2 | y) = p(sigma2 | y) p(mu | sigma2, y).
# The marginal p(sigma2 | y) is Inverse-Gamma in closed form; the conditional
# p(mu | sigma2, y) is Normal. Drawing in that order gives exact joint draws.
#
# Iterations are INDEPENDENT of one another. Within a single iteration mu and
# sigma2 are DEPENDENT: a draw with large sigma2 comes paired with a more
# dispersed mu. Averaging over iterations therefore propagates variance
# uncertainty into the mean, which is the whole point.

nig_draws <- function(y, n_draws = 20000, m0 = 0, k0 = 1e-6,
                      a0 = 1e-3, b0 = 1e-3, return_sigma2 = FALSE) {
  n <- length(y); ybar <- mean(y); ss <- sum((y - ybar)^2)

  kn <- k0 + n
  mn <- (k0 * m0 + n * ybar) / kn
  an <- a0 + n / 2
  bn <- b0 + 0.5 * ss + (k0 * n * (ybar - m0)^2) / (2 * kn)

  sigma2 <- bn / stats::rgamma(n_draws, shape = an, rate = 1)  # InvGamma(an, bn)
  mu     <- mn + stats::rnorm(n_draws) * sqrt(sigma2 / kn)

  if (return_sigma2) {
    return(list(mu = mu, sigma2 = sigma2,
                hyper = c(n = n, ybar = ybar, ss = ss,
                          kn = kn, mn = mn, an = an, bn = bn)))
  }
  mu
}

# n_draws controls Monte Carlo error in the posterior SUMMARIES only.
# It is not a sample size and carries no information about participants.
cgr_conjugate <- function(df, grid = seq(0, 1, length.out = 101),
                          n_draws = 20000, prior = list(),
                          direction = 1) {
  st  <- cgr_strata(df); rat <- cgr_ratios(st)
  mu  <- lapply(st, function(y)
    do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior)))

  d <- lapply(grid, cgr_delta, mu = mu, r = rat$r, s = rat$s)
  sdv <- vapply(d, stats::sd, numeric(1))

  out <- data.frame(
    cgr = grid, method = "conjugate",
    est = vapply(d, mean, numeric(1)), sd = sdv,
    lo = vapply(d, function(x) unname(stats::quantile(x, 0.025)), numeric(1)),
    hi = vapply(d, function(x) unname(stats::quantile(x, 0.975)), numeric(1)),
    # P(favourable): direction = +1 if higher is better, -1 if lower is better
    p_fav = vapply(d, function(x) mean(direction * x > 0), numeric(1)),
    stringsAsFactors = FALSE
  )
  out$mcse <- sdv / sqrt(n_draws)   # iid draws: no autocorrelation correction
  out
}
