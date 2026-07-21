# Front door. Two entry points that map onto what a researcher actually wants:
# adjust my trial, and (as the original paper's limitations section asked) first
# check whether the adjustment is trustworthy for my design. The second is
# cgr_operating() in R/05_sim.R; this file adds the one-call adjuster.

# cgrc(df): the CGR-adjusted analysis of one trial in a single call. `df` needs
# columns condition (AC/PL), guess (AC/PL) and value. direction = +1 if higher
# scores are better, -1 if lower are better. Returns the posterior curve, a
# summary at the observed CGR and at perfect blinding (0.50), and the observed
# CGR. The grid always includes the exact observed CGR, so the unadjusted row is
# never read off a snapped grid point.
cgrc <- function(df, n_draws = 20000, direction = 1, prior = list()) {
  o    <- cgr_observed(cgr_strata(df))
  grid <- sort(unique(c(seq(0, 1, length.out = 101), o)))
  cur  <- cgr_conjugate(df, grid, n_draws = n_draws,
                        direction = direction, prior = prior)
  structure(list(curve = cur,
                 summary = cgr_summary_table(cur, o),
                 observed_cgr = o),
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
