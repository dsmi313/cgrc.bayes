# Plotting. The primary Bayesian figure shows posterior mean, 95% credible
# interval, and P(favourable) - NOT a p-value-shaped curve. There is no
# horizontal line at 0.05 anywhere in this file.

cgr_plot <- function(cur, obs_cgr, title = NULL,
                     show_95_line = TRUE, direction_label = "positive") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("needs ggplot2", call. = FALSE)
  }
  lev <- c("Treatment effect (posterior mean, 95% CrI)",
           sprintf("Posterior probability the effect is %s", direction_label))

  eff <- data.frame(cgr = cur$cgr, method = cur$method, est = cur$est,
                    lo = cur$lo, hi = cur$hi, quantity = lev[1],
                    stringsAsFactors = FALSE)
  pp  <- data.frame(cgr = cur$cgr, method = cur$method, est = cur$p_fav,
                    lo = NA_real_, hi = NA_real_, quantity = lev[2],
                    stringsAsFactors = FALSE)
  dat <- rbind(eff, pp); dat$quantity <- factor(dat$quantity, levels = lev)

  href <- data.frame(quantity = factor(lev, levels = lev),
                     y = c(0, if (show_95_line) 0.95 else NA_real_))
  href <- href[!is.na(href$y), ]

  pal <- c(analytic = "grey30", conjugate = "#C0392B",
           jags = "#2471A3", `jags-t` = "#117A65")

  ggplot2::ggplot() +
    ggplot2::geom_hline(data = href, ggplot2::aes(yintercept = y),
                        linetype = "dotted", colour = "grey40") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        colour = "black") +
    ggplot2::geom_vline(xintercept = obs_cgr, linetype = "dashed",
                        colour = "darkgreen") +
    ggplot2::geom_ribbon(data = dat[!is.na(dat$lo), ],
                         ggplot2::aes(x = cgr, ymin = lo, ymax = hi,
                                      fill = method),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(data = dat,
                       ggplot2::aes(x = cgr, y = est, colour = method),
                       linewidth = 0.7) +
    ggplot2::facet_wrap(~ quantity, ncol = 1, scales = "free_y",
                        strip.position = "left") +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::labs(x = "Correct guess rate (CGR)", y = NULL,
                  colour = NULL, fill = NULL, title = title,
                  subtitle = sprintf(paste("black dashed = perfect blinding",
                                           "(0.50); green dashed = observed",
                                           "CGR (%.3f)"), obs_cgr)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(strip.placement = "outside",
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top")
}

# Summary at the two CGRs that matter, with attenuation reported both ways.
#
# The observed-CGR row must be read off the curve at EXACTLY the observed CGR,
# not the nearest grid point. On a 101-point grid c_obs = 0.6466 would snap to
# 0.65, and Delta(0.65) overstates Delta(0.6466) by ~0.05 - the same class of
# grid error this document criticises elsewhere. `at()` therefore warns if the
# grid does not contain the target within `tol`; pass a grid that includes
# c_obs (e.g. sort(unique(c(grid, c_obs)))) so no snapping occurs.
cgr_summary_table <- function(cur, obs_cgr, label = "", tol = 1e-6) {
  at <- function(t) {
    i <- which.min(abs(cur$cgr - t))
    if (abs(cur$cgr[i] - t) > tol) {
      warning(sprintf(paste0(
        "cgr_summary_table: nearest grid point %.4f is %.2e from the requested ",
        "%.4f; include the target in the grid to avoid grid-snapping."),
        cur$cgr[i], abs(cur$cgr[i] - t), t), call. = FALSE)
    }
    cur[i, ]
  }
  a <- at(obs_cgr); h <- at(0.5)
  # Percentage attenuation is a ratio to the unadjusted estimate. Suppress it
  # when that estimate is not distinguishable from zero (its 95% CrI includes
  # zero): dividing by a near-zero, sign-ambiguous denominator is meaningless
  # (e.g. cognitive performance, unadj -0.011 with CrI crossing 0).
  unadj_distinct <- !(a$lo <= 0 && a$hi >= 0)
  pct <- if (unadj_distinct) 100 * (a$est - h$est) / a$est else NA_real_
  data.frame(
    outcome = label,
    cgr = c(round(obs_cgr, 4), 0.5),
    what = c("observed (unadjusted)", "perfect blinding (adjusted)"),
    post_mean = round(c(a$est, h$est), 3),
    cri_lo = round(c(a$lo, h$lo), 3),
    cri_hi = round(c(a$hi, h$hi), 3),
    p_favourable = round(c(a$p_fav, h$p_fav), 3),
    abs_attenuation = round(c(NA, a$est - h$est), 3),
    pct_attenuation = round(c(NA, pct), 1),
    stringsAsFactors = FALSE
  )
}
