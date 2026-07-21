# AEB generative model (Szigeti Fig 1, Eq 1-4) and the repeated-simulation
# study used to establish operating characteristics.
#
#   Eq 1  p(B_TRT = 1) = 0.5
#   Eq 2  p(B_PT = B_TRT) = p_CG
#   Eq 3  B_TE = B_PT
#   Eq 4  OUT = N_NH + B_TRT*N_DTE + B_TE*N_AEB
#
# `noise` resolves the ambiguity in Eq 4. Default "all" is EVIDENCE-BASED, not
# a convenience choice; see reports/CHANGELOG.md F1:
#   noise="all": Hedges g = 0.4011 (published 0.40); Table 1 unadjusted
#                significance 0.056/0.876/0.826/0.992 (published
#                0.05/0.86/0.78/0.99)
#   noise="arm": Hedges g = 0.5023; rates 0.056/0.970/0.930/0.998
# The author's source code was NOT located, so this is convergent empirical
# evidence rather than verification against the original implementation.

sim_aeb <- function(n = 230, p_cg = 0.7, dte_on = FALSE, aeb_on = FALSE,
                    noise = c("all", "arm"),
                    mu_nh = 10, sd_nh = 4, mu_dte = 3, sd_dte = 6.2,
                    mu_aeb = 7.7, sd_aeb = 6.2) {
  noise <- match.arg(noise)
  trt     <- stats::runif(n) < 0.5
  correct <- stats::runif(n) < p_cg
  pt      <- ifelse(correct, trt, !trt)
  y <- stats::rnorm(n, mu_nh, sd_nh)
  if (dte_on) y <- y + if (noise == "arm") trt * stats::rnorm(n, mu_dte, sd_dte)
                       else trt * mu_dte + stats::rnorm(n, 0, sd_dte)
  if (aeb_on) y <- y + if (noise == "arm") pt * stats::rnorm(n, mu_aeb, sd_aeb)
                       else pt * mu_aeb + stats::rnorm(n, 0, sd_aeb)
  data.frame(condition = ifelse(trt, "AC", "PL"),
             guess = ifelse(pt, "AC", "PL"), value = y,
             stringsAsFactors = FALSE)
}

# Operating characteristics over many INDEPENDENT trials. A single simulated
# trial is an illustration, not validation.
#
# True Delta(0.5) equals the direct treatment effect: at perfect blinding the
# guess distribution is identical in both arms, so the AEB term contributes
# equally to each and cancels.
cgr_operating <- function(n_trials = 500, n = 230, p_cg = 0.7,
                          noise = "all", n_draws = 4000, seed = 1) {
  set.seed(seed)
  cfg <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  rows <- lapply(cfg, function(z) {
    dte <- as.logical(z[1]); aeb <- as.logical(z[2])
    truth <- if (dte) 3 else 0
    adj <- unadj <- numeric(n_trials)
    cov <- fav <- sig <- logical(n_trials)
    for (i in seq_len(n_trials)) {
      d  <- sim_aeb(n, p_cg, dte, aeb, noise)
      st <- cgr_strata(d); rat <- cgr_ratios(st)
      mu <- lapply(st, nig_draws, n_draws = n_draws)
      dd <- cgr_delta(0.5, mu, rat$r, rat$s)
      q  <- stats::quantile(dd, c(0.025, 0.975))
      adj[i] <- mean(dd)
      cov[i] <- q[1] <= truth && truth <= q[2]
      fav[i] <- mean(dd > 0) > 0.95
      a <- d$value[d$condition == "AC"]; b <- d$value[d$condition == "PL"]
      unadj[i] <- mean(a) - mean(b)
      sig[i]   <- stats::t.test(a, b, var.equal = TRUE)$p.value < 0.05
    }
    data.frame(DTE = z[1], AEB = z[2], true = truth,
               unadj_mean = mean(unadj), unadj_bias = mean(unadj) - truth,
               adj_mean = mean(adj), adj_bias = mean(adj) - truth,
               adj_rmse = sqrt(mean((adj - truth)^2)),
               coverage95 = mean(cov), p_fav_gt_95 = mean(fav),
               freq_sig = mean(sig))
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}
