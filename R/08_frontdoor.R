# Front door. Two entry points that map onto what a researcher actually wants:
# adjust my trial, and (as the original paper's limitations section asked) first
# check whether the adjustment is trustworthy for my design. The second is
# cgr_operating() in R/05_sim.R; this file adds the one-call adjuster.

# cgrc(df): the CGR-adjusted analysis of one trial in a single call. `df` needs
# columns condition (AC/PL), guess (AC/PL) and value. direction = +1 if higher
# scores are better, -1 if lower are better. Returns the posterior curve, a
# summary at the observed CGR and at a target CGR of 0.50 (guessing at chance,
# not proof of perfect blinding), and the observed
# CGR. The grid always includes the exact observed CGR, so the unadjusted row is
# never read off a snapped grid point.
# `seed` (optional): when given, set.seed(seed) is called before drawing and the
# value is recorded in the returned object for reproducible reports. Default NULL
# preserves the previous behaviour (reproducibility via an external set.seed).
cgrc <- function(df, n_draws = 20000, direction = 1, prior = list(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  o    <- cgr_observed(cgr_strata(df))
  grid <- sort(unique(c(seq(0, 1, length.out = 101), o)))
  cur  <- cgr_conjugate(df, grid, n_draws = n_draws,
                        direction = direction, prior = prior)
  structure(list(curve = cur,
                 summary = cgr_summary_table(cur, o),
                 observed_cgr = o, seed = seed),
            class = "cgrc")
}

print.cgrc <- function(x, ...) {
  cat(sprintf("CGRC-adjusted analysis  (observed CGR = %.4f)\n", x$observed_cgr))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

# plot(cgrc(df)) draws the effect curve with its 95% credible band and the
# P(favourable) panel. Needs ggplot2.
plot.cgrc <- function(x, ...) {
  cgr_plot(x$curve, obs_cgr = x$observed_cgr, ...)
}

# cgrc_headline(df): the most interpretable summary of a CGR-adjusted analysis -
# TWO plain probabilities, before and after the blinding correction, because a
# trialist has two questions the p-value conflates:
#   * "is there an effect?"       -> P(favourable) = P(direction * Delta > 0)
#   * "is it big enough to care?"  -> P(meaningful) = P(direction * Delta > delta)
# Each is reported at the observed CGR (raw) and at a target CGR 0.50 (guessing
# at chance),
# with the adjusted point estimate and 95% CrI. No bright-line threshold is
# imposed - these are continuous probabilities, deliberately not turned back into
# a "significant/not" verdict. delta defaults to 0.5 * outcome SD: half a standard
# deviation is the minimum important difference Norman (2003) argues for and the
# one Szigeti's own 2024 escitalopram trial adopts, so "meaningful" here means
# clinically meaningful rather than merely non-zero. Narrow delta_sd_frac for a
# stricter view. direction = +1 if higher scores are better, -1 if lower are better.
cgrc_headline <- function(df, direction = 1, delta = NULL, delta_sd_frac = 0.5,
                          n_draws = 20000, prior = list(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  st  <- cgr_strata(df); rat <- cgr_ratios(st)
  if (is.null(delta)) delta <- delta_sd_frac * stats::sd(df$value)
  mu  <- lapply(st, function(y)
    do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior)))
  o   <- cgr_observed(st)
  d_obs <- cgr_delta(o,   mu, rat$r, rat$s)
  d_bl  <- cgr_delta(0.5, mu, rat$r, rat$s)
  eff_obs <- direction * d_obs; eff_bl <- direction * d_bl   # + = favourable
  q <- function(x, p) unname(stats::quantile(x, p))
  res <- list(
    observed_cgr = o, delta = delta, direction = direction, seed = seed,
    adj_est = mean(d_bl), adj_lo = q(d_bl, 0.025), adj_hi = q(d_bl, 0.975),
    p_dir_obs        = mean(eff_obs > 0),     p_dir_blind        = mean(eff_bl > 0),
    p_meaningful_obs = mean(eff_obs > delta), p_meaningful_blind = mean(eff_bl > delta))
  dtxt <- trimws(sprintf("%.2g", delta))
  unit <- if (identical(dtxt, "1")) "point" else "points"
  res$text <- sprintf(paste0(
    "Raw, this trial shows a %.0f%% probability of a favourable effect and %.0f%% ",
    "that it is meaningful (beyond %s %s). Reweighted to a correct-guess rate of ",
    "0.50 (guessing at chance), those become %.0f%% and %.0f%% (adjusted effect ",
    "%.2f, 95%% CrI %.2f to %.2f)."),
    100 * res$p_dir_obs, 100 * res$p_meaningful_obs, dtxt, unit,
    100 * res$p_dir_blind, 100 * res$p_meaningful_blind,
    res$adj_est, res$adj_lo, res$adj_hi)
  structure(res, class = "cgrc_headline")
}

print.cgrc_headline <- function(x, ...) { cat(x$text, "\n"); invisible(x) }
