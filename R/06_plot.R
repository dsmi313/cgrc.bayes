# Plotting. The primary Bayesian figure shows posterior mean, 95% credible
# interval, and P(favourable) - NOT a p-value-shaped curve. There is no
# horizontal line at 0.05 anywhere in this file.

cgr_plot <- function(cur, obs_cgr, title = NULL,
                     show_95_line = TRUE, direction_label = "positive") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("needs ggplot2", call. = FALSE)
  }
  lev <- c("Treatment effect\n(95% CrI)",
           sprintf("P(effect %s)", direction_label))

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
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(strip.placement = "outside",
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top")
}

# Reproduce Szigeti's published twin-axis figure in its OWN visual grammar so a
# reader can lay it beside the paper: blue averaged p-value on the left axis, red
# treatment estimate on the right, a magenta 0.05 threshold, a green "Original
# CGR" line and a black "True blind CGR" (0.5) line. `z` is the output of
# cgr_kde_curve(). This deliberately plots the p-value presentation the rest of
# the package argues against - because that is what reproduction means.
szigeti_panel <- function(z, title, orig_cgr, blind_cgr = 0.5, legend = FALSE) {
  lo_p <- -0.05; hi_p <- 0.75                       # left axis, as published
  elo <- min(z$est - z$est_sd, na.rm = TRUE)
  ehi <- max(z$est + z$est_sd, na.rm = TRUE)
  pad <- 0.12 * (ehi - elo); elo <- elo - pad; ehi <- ehi + pad
  sc  <- function(v) (v - elo) / (ehi - elo) * (hi_p - lo_p) + lo_p

  op <- graphics::par(mar = c(4.2, 4.2, 3, 4.6), bg = "white")
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(NA, xlim = c(0, 1), ylim = c(lo_p, hi_p), xaxs = "i",
                 xlab = "Correct guess rate (CGR)", ylab = "", axes = FALSE,
                 main = title, font.main = 2, cex.main = 1.1)
  graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3],
                 graphics::par("usr")[2], graphics::par("usr")[4],
                 col = "grey92", border = NA)
  graphics::grid(col = "white", lty = 1, lwd = 1.1)

  graphics::polygon(c(z$cgr, rev(z$cgr)),
                    c(sc(z$est - z$est_sd), rev(sc(z$est + z$est_sd))),
                    col = grDevices::adjustcolor("red", 0.22), border = NA)
  graphics::lines(z$cgr, sc(z$est), col = "red", lwd = 2)

  graphics::polygon(c(z$cgr, rev(z$cgr)),
                    c(pmax(z$p - z$p_sd, lo_p), rev(pmin(z$p + z$p_sd, hi_p))),
                    col = grDevices::adjustcolor("blue", 0.22), border = NA)
  graphics::lines(z$cgr, z$p, col = "blue", lwd = 2)

  graphics::abline(h = 0.05, col = "magenta", lty = 2, lwd = 1.4)
  graphics::abline(v = blind_cgr, col = "black",     lty = 2, lwd = 1.4)
  graphics::abline(v = orig_cgr,  col = "darkgreen", lty = 2, lwd = 1.4)

  graphics::axis(1); graphics::box()
  graphics::axis(2, at = seq(-0.05, 0.75, by = 0.10), las = 1, col.axis = "blue")
  graphics::mtext("Treatment p-value", side = 2, line = 2.9, col = "blue")
  tk <- pretty(c(elo, ehi), 8); tk <- tk[tk >= elo & tk <= ehi]
  graphics::axis(4, at = sc(tk), labels = format(tk, trim = TRUE), las = 1,
                 col.axis = "red")
  graphics::mtext("Treatment estimate", side = 4, line = 3.2, col = "red")

  if (legend) graphics::legend("topleft", bty = "n", cex = 0.72,
    legend = c("Treatment p-value", "Treatment estimate",
               "Original CGR", "True blind CGR", "Sig. threshold"),
    col = c("blue", "red", "darkgreen", "black", "magenta"),
    lty = c(1, 1, 2, 2, 2), lwd = c(2, 2, 1.4, 1.4, 1.4))
  invisible(NULL)
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
