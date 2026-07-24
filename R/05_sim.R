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
# `direction` orients favourability (+1 higher-is-better, -1 lower-is-better) and
# is applied identically to the Bayesian flags and the frequentist criterion, so
# the two are compared on the SAME tail. At the default +1 every column is
# byte-identical to the older two-sided build EXCEPT freq_sig (see below).
cgr_operating <- function(n_trials = 500, n = 230, p_cg = 0.7,
                          noise = "all", n_draws = 4000, seed = 1,
                          mu_dte = 3, mu_aeb = 7.7, direction = 1) {
  set.seed(seed)
  cfg <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  rows <- lapply(cfg, function(z) {
    dte <- as.logical(z[1]); aeb <- as.logical(z[2])
    truth <- if (dte) mu_dte else 0
    adj <- unadj <- numeric(n_trials)
    cov <- fav <- fav975 <- sig <- valid <- logical(n_trials)
    for (i in seq_len(n_trials)) {
      d  <- sim_aeb(n, p_cg, dte, aeb, noise, mu_dte = mu_dte, mu_aeb = mu_aeb)
      # At a high correct-guess rate with small n, a wrong-guess stratum can come
      # up empty, and the estimand is then undefined. Skip that trial rather than
      # crash - the RATE of such trials is itself a warning that CGR adjustment
      # is fragile at these parameters, and it is reported below.
      if (!all(STRATA %in% paste0(d$condition, d$guess))) next
      valid[i] <- TRUE
      st <- cgr_strata(d); rat <- cgr_ratios(st)
      mu <- lapply(st, nig_draws, n_draws = n_draws)
      dd <- cgr_delta(0.5, mu, rat$r, rat$s)
      q  <- stats::quantile(dd, c(0.025, 0.975))
      pd <- mean(direction * dd > 0)        # posterior prob. of a FAVOURABLE effect
      adj[i] <- mean(dd)
      cov[i] <- q[1] <= truth && truth <= q[2]
      fav[i]    <- pd > 0.95                # one-sided flag (standard)
      fav975[i] <- pd > 0.975               # approx. matched to a two-sided p<0.05
      a <- d$value[d$condition == "AC"]; b <- d$value[d$condition == "PL"]
      raw_effect <- mean(a) - mean(b)
      unadj[i] <- raw_effect
      # Direction-MATCHED frequentist event: significant AND in the favourable
      # direction. The Bayesian flags count only the favourable tail, so counting
      # two-sided significance here (as older builds did) compared unlike tails and
      # overstated the match. freq_sig is the ONLY column this changes at dir = +1.
      sig[i] <- stats::t.test(a, b, var.equal = TRUE)$p.value < 0.05 &&
                (direction * raw_effect) > 0
    }
    v <- which(valid); nv <- length(v)
    m <- function(x) if (nv) mean(x[v]) else NA_real_
    data.frame(DTE = z[1], AEB = z[2], true = truth,
               unadj_mean = m(unadj), unadj_bias = m(unadj) - truth,
               adj_mean = m(adj), adj_bias = m(adj) - truth,
               adj_rmse = if (nv) sqrt(mean((adj[v] - truth)^2)) else NA_real_,
               coverage95 = m(cov), p_fav_gt_95 = m(fav),
               p_fav_gt_975 = m(fav975), freq_sig = m(sig),
               empty_stratum_rate = round(mean(!valid), 4), n_valid = nv)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Expected size of the smallest stratum for a design of size n at correct-guess
# rate p_cg. Under balanced allocation the four strata have expected shares
# 0.5*p_cg (the two concordant/correct strata) and 0.5*(1 - p_cg) (the two
# discordant/wrong strata), so the smallest is n * 0.5 * min(p_cg, 1 - p_cg).
# Closed form, no simulation - the single best early warning that CGR adjustment
# may be infeasible: when this drops below ~15 the discordant strata are thin and
# simulated trials start coming up empty (see cgr_operating()'s empty_stratum_rate).
cgr_min_stratum <- function(n, p_cg) n * 0.5 * pmin(p_cg, 1 - p_cg)

# Expected inflation of an UNADJUSTED treatment estimate caused by activated
# expectancy bias, in outcome points. In the AEB model the expectancy term adds
# mu_aeb * E[B_TE | arm] to each arm's mean, and B_TE (perceived treatment) has
# E[B_TE | active] = p_cg and E[B_TE | placebo] = 1 - p_cg, so an unadjusted
# analysis picks up mu_aeb * (p_cg - (1 - p_cg)) = mu_aeb * (2 * p_cg - 1).
# Closed form, no simulation (matches simulation to < 0.02 points). Turns the
# opaque expectancy magnitude into a number a user can weigh against their own
# effect size: "at CGR 0.85, mu_aeb = 7.7 inflates an unadjusted estimate by
# 5.4 points".
cgr_aeb_inflation <- function(mu_aeb, p_cg) mu_aeb * (2 * p_cg - 1)

# ---- UNKNOWN-aware generative model and operating characteristics ------------
#
# Extends the AEB model (sim_aeb) with an observed "I do not know" response, to
# characterise the UNKNOWN-preserving estimator (R/10_unknown.R) the way
# cgr_operating() characterises the binary one. This is a NEW generative model
# with EXPLICIT assumptions, developed deliberately (not invented): see
# reports/UNRESOLVED.md U10.
#
# Assumptions (both chosen deliberately; alternatives are stress tests, not the
# default):
#   A1  Each participant answers UNKNOWN with probability `u`, INDEPENDENT of arm.
#   A2  An UNKNOWN responder carries NO expectancy (they state no directional
#       belief, so the perceived-treatment -> expectancy path is absent: B_TE = 0).
# Among directional responders the guess is correct with probability p_cg (the
# directional CGR the estimator conditions on), and expectancy follows the
# perceived treatment exactly as in sim_aeb. `noise` matches sim_aeb.
#
# Consequence used by the operating-characteristics study: at directional
# CGR 0.50 with the UNKNOWN rate held fixed, the directional expectancy is
# balanced across arms and the UNKNOWN mass is inert, so Delta(0.50, u_obs)
# targets the direct treatment effect mu_dte - which coverage/bias then check.
sim_aeb_unknown <- function(n = 230, p_cg = 0.7, u = 0.2,
                            dte_on = FALSE, aeb_on = FALSE, noise = c("all", "arm"),
                            mu_nh = 10, sd_nh = 4, mu_dte = 3, sd_dte = 6.2,
                            mu_aeb = 7.7, sd_aeb = 6.2) {
  noise <- match.arg(noise)
  trt     <- stats::runif(n) < 0.5
  is_unk  <- stats::runif(n) < u                 # A1: arm-independent UNKNOWN rate
  correct <- stats::runif(n) < p_cg
  pt      <- ifelse(correct, trt, !trt)          # perceived treatment (directional)
  b_te    <- ifelse(is_unk, 0, as.numeric(pt))   # A2: UNKNOWN -> no expectancy
  y <- stats::rnorm(n, mu_nh, sd_nh)
  if (dte_on) y <- y + if (noise == "arm") trt * stats::rnorm(n, mu_dte, sd_dte)
                       else trt * mu_dte + stats::rnorm(n, 0, sd_dte)
  if (aeb_on) y <- y + if (noise == "arm") b_te * stats::rnorm(n, mu_aeb, sd_aeb)
                       else b_te * mu_aeb + stats::rnorm(n, 0, sd_aeb)
  data.frame(condition = ifelse(trt, "AC", "PL"),
             guess = ifelse(is_unk, "UNKNOWN", ifelse(pt, "AC", "PL")),
             value = y, stringsAsFactors = FALSE)
}

# Operating characteristics of the UNKNOWN-PRESERVING estimator over many
# independent trials, the six-stratum analogue of cgr_operating(). At each trial
# it fits the six-stratum conjugate posterior and evaluates Delta(0.50, u_obs).
# Reports adjusted bias, RMSE, 95% coverage, favourable-flag rates (one-sided
# > 0.95 and matched > 0.975), the unadjusted t-test significance rate, and the
# empty-stratum rate. The six strata are thinner than four, so `empty_stratum_rate`
# here is the honest early warning that the extension is fragile at small n / high
# guess or UNKNOWN rates - a trial with an empty directional class is skipped
# (the estimand is undefined) and counted, never fabricated.
cgr_unknown_operating <- function(n_trials = 500, n = 230, p_cg = 0.7, u = 0.2,
                                  noise = "all", n_draws = 4000, seed = 1,
                                  mu_dte = 3, mu_aeb = 7.7, direction = 1) {
  set.seed(seed)
  cfg <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  rows <- lapply(cfg, function(z) {
    dte <- as.logical(z[1]); aeb <- as.logical(z[2])
    truth <- if (dte) mu_dte else 0
    adj <- unadj <- numeric(n_trials)
    cov <- fav <- fav975 <- sig <- valid <- logical(n_trials)
    for (i in seq_len(n_trials)) {
      d  <- sim_aeb_unknown(n, p_cg, u, dte, aeb, noise, mu_dte = mu_dte, mu_aeb = mu_aeb)
      st <- cgr_unknown_strata(d); o <- cgr_unknown_observed(st)
      # estimand defined at c = 0.5, u = u_obs: both directional classes present,
      # and the UNKNOWN class present if u_obs > 0. Otherwise skip and count it.
      if (o$n_correct == 0 || o$n_incorrect == 0 ||
          (o$u_obs > 0 && o$n_unknown == 0)) next
      valid[i] <- TRUE
      rat <- cgr_unknown_ratios(st)
      mu  <- stats::setNames(lapply(UNKNOWN_STRATA, function(nm) {
        y <- st[[nm]]; if (length(y)) nig_draws(y, n_draws = n_draws) else NA_real_
      }), UNKNOWN_STRATA)
      dd <- cgr_unknown_delta(0.5, o$u_obs, mu, rat$r, rat$s, rat$t)
      q  <- stats::quantile(dd, c(0.025, 0.975)); pd <- mean(direction * dd > 0)
      adj[i] <- mean(dd)
      cov[i] <- q[1] <= truth && truth <= q[2]
      fav[i] <- pd > 0.95; fav975[i] <- pd > 0.975   # matched flag, favourable tail
      a <- d$value[d$condition == "AC"]; b <- d$value[d$condition == "PL"]
      raw_effect <- mean(a) - mean(b)
      unadj[i] <- raw_effect
      # Direction-matched frequentist event (see cgr_operating): favourable tail
      # only, so it compares like-for-like with the Bayesian flags.
      sig[i] <- stats::t.test(a, b, var.equal = TRUE)$p.value < 0.05 &&
                (direction * raw_effect) > 0
    }
    v <- which(valid); nv <- length(v)
    m <- function(x) if (nv) mean(x[v]) else NA_real_
    data.frame(DTE = z[1], AEB = z[2], true = truth, u = u,
               unadj_mean = m(unadj), unadj_bias = m(unadj) - truth,
               adj_mean = m(adj), adj_bias = m(adj) - truth,
               adj_rmse = if (nv) sqrt(mean((adj[v] - truth)^2)) else NA_real_,
               coverage95 = m(cov), p_fav_gt_95 = m(fav),
               p_fav_gt_975 = m(fav975), freq_sig = m(sig),
               empty_stratum_rate = round(mean(!valid), 4), n_valid = nv)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Expected size of the smallest of the SIX strata for an UNKNOWN-aware design.
# With arm-independent UNKNOWN rate u, the directional mass is n(1-u), split into
# correct 0.5*p_cg and wrong 0.5*(1-p_cg) shares per arm; the UNKNOWN mass is n*u
# split ~0.5 per arm. The smallest expected cell is the thinnest of these - the
# feasibility early warning for the six-stratum estimand (cf. cgr_min_stratum).
cgr_unknown_min_stratum <- function(n, p_cg, u) {
  # thinnest OCCUPIED cell: a structurally-empty class (u = 0 or u = 1) is not
  # counted, so at u = 0 this reduces to the binary cgr_min_stratum().
  cands <- c(0.5 * (1 - u) * p_cg, 0.5 * (1 - u) * (1 - p_cg), 0.5 * u)
  n * min(cands[cands > 0])
}
