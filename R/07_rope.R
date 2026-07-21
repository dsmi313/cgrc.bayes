# Region of practical equivalence (ROPE) decomposition of the CGR-adjusted
# contrast - a more robust inferential summary than a single tail probability.
#
# For each grid point the posterior of Delta(c) is split into three exhaustive,
# mutually exclusive regions relative to a band [-delta, +delta]:
#
#   p_harm       = P(direction * Delta < -delta | y)
#   p_negligible = P(|Delta| <= delta         | y)
#   p_benefit    = P(direction * Delta >  delta | y)
#
# The three sum to 1 at every CGR. `direction` orients favourability: +1 if
# higher scores are better, -1 if lower are better.
#
# Why this instead of the two-sided tail probability. The tail statistic is
# magnitude-blind (a tiny precise effect and a large uncertain one can score
# alike) and it conflates "the posterior is concentrated near zero, so the
# effect is negligible" with "the posterior is diffuse, so we have learned
# nothing". The middle region above separates exactly those two states. The
# cost is that delta must be declared - which is a feature: it forces the
# practical-significance question into the open. Default is 0.1 * the outcome's
# pooled SD (a conventional choice); cgr_rope_sensitivity() shows how the
# conclusion moves as delta widens.
#
# A Bayes factor against the point null Delta = 0 is deliberately NOT offered:
# with the vague default prior (k0 = 1e-6) it depends on the prior width without
# limit (the Jeffreys-Lindley paradox), so it would measure the prior, not the
# data. A region-based summary has no such failure mode.
cgr_rope <- function(df, grid = seq(0, 1, length.out = 101),
                     n_draws = 20000, delta = NULL, delta_sd_frac = 0.1,
                     direction = 1, prior = list()) {
  st  <- cgr_strata(df); rat <- cgr_ratios(st)
  if (is.null(delta)) delta <- delta_sd_frac * stats::sd(df$value)
  mu  <- lapply(st, function(y)
    do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior)))
  d   <- lapply(grid, cgr_delta, mu = mu, r = rat$r, s = rat$s)
  q   <- function(x, p) unname(stats::quantile(x, p))

  data.frame(
    cgr = grid, delta = delta,
    est   = vapply(d, mean, numeric(1)),
    lo95  = vapply(d, function(x) q(x, 0.025), numeric(1)),
    hi95  = vapply(d, function(x) q(x, 0.975), numeric(1)),
    lo50  = vapply(d, function(x) q(x, 0.250), numeric(1)),
    hi50  = vapply(d, function(x) q(x, 0.750), numeric(1)),
    # region probabilities are oriented by `direction` so + always = benefit
    p_harm      = vapply(d, function(x) mean(direction * x < -delta), numeric(1)),
    p_negligible= vapply(d, function(x) mean(abs(x) <= delta),        numeric(1)),
    p_benefit   = vapply(d, function(x) mean(direction * x >  delta), numeric(1)),
    stringsAsFactors = FALSE)
}

# How the ROPE conclusion at a single CGR moves as the width delta is varied.
# A region-based conclusion is only as robust as the region; this is the check
# that keeps it honest. Returns P(negligible) and P(benefit) at `at_cgr` for a
# ladder of delta = frac * sd(outcome).
cgr_rope_sensitivity <- function(df, at_cgr = 0.5,
                                 fracs = c(0.05, 0.10, 0.20, 0.30),
                                 n_draws = 12000, direction = 1) {
  sdy <- stats::sd(df$value)
  do.call(rbind, lapply(fracs, function(f) {
    z <- cgr_rope(df, grid = at_cgr, n_draws = n_draws,
                  delta_sd_frac = f, direction = direction)
    data.frame(delta_in_SD = f, delta = f * sdy,
               p_negligible = z$p_negligible, p_benefit = z$p_benefit,
               stringsAsFactors = FALSE)
  }))
}
