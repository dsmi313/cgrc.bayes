# UNKNOWN-preserving CGRC extension.
#
# This is an EXTENSION implemented by cgrc.bayes. It is NOT part of the original
# CGRC formulation of Szigeti et al. and is not claimed by them. The binary
# four-stratum estimand in R/01_estimand.R is unchanged; nothing here alters it.
#
# Real trials often let a participant answer "I do not know" to the guess
# question. Those participants must not be dropped, counted as wrong, counted as
# placebo, or split into AC/PL without an explicit assumption. This file keeps
# UNKNOWN as an observed third response category via six strata:
#
#   ACAC received active,  guessed active
#   ACPL received active,  guessed placebo
#   ACU  received active,  answered UNKNOWN
#   PLAC received placebo, guessed active
#   PLPL received placebo, guessed placebo
#   PLU  received placebo, answered UNKNOWN
#
# The estimand holds the observed UNKNOWN-response rate u fixed and varies the
# DIRECTIONAL correct-guess rate c (the correct-guess rate AMONG participants who
# gave an AC/PL guess). At u = 0 it reduces exactly to the binary CGRC; see the
# tests in tests/testthat/test-unknown.R.
#
# Reweighting to c = 0.50 is NOT "perfect blinding" without qualification: it is
# "directional guessing at chance while holding the observed UNKNOWN-response
# rate fixed". The wording throughout is deliberate.

UNKNOWN_STRATA <- c("ACAC", "ACPL", "ACU", "PLAC", "PLPL", "PLU")

# Split a normalised trial (condition in {AC,PL}, guess in {AC,PL,UNKNOWN}) into
# the six strata. Unlike cgr_strata(), this does NOT error on an empty cell: some
# cells can be legitimately empty and are then structurally zero-weighted (their
# preserved within-class arm share is 0). Estimability is checked separately by
# cgr_unknown_estimable(). Returns a length-6 list in UNKNOWN_STRATA order, with
# empty cells present as zero-length vectors.
cgr_unknown_strata <- function(df) {
  g <- as.character(df$guess)
  gcode <- ifelse(g == "UNKNOWN", "U", g)          # AC/PL/UNKNOWN -> AC/PL/U
  key <- factor(paste0(df$condition, gcode), levels = UNKNOWN_STRATA)
  if (anyNA(key)) {
    bad <- unique(paste0(df$condition, "/", g)[is.na(key)])
    stop("condition must be AC/PL and guess must be AC/PL/UNKNOWN; offending ",
         "condition/guess: ", paste(shQuote(bad), collapse = ", "), call. = FALSE)
  }
  split(df$value, key)                              # keeps empty levels (drop=FALSE)
}

# Observed counts and rates. c_obs is the DIRECTIONAL correct-guess rate (correct
# among AC/PL responders); u_obs is the UNKNOWN-response rate over all responders.
cgr_unknown_observed <- function(st) {
  n <- lengths(st)
  n_total    <- sum(n)
  n_unknown  <- n[["ACU"]] + n[["PLU"]]
  n_direct   <- n_total - n_unknown
  n_correct  <- n[["ACAC"]] + n[["PLPL"]]
  n_incorrect<- n[["ACPL"]] + n[["PLAC"]]
  list(
    counts = n[UNKNOWN_STRATA],
    n_total = n_total, n_unknown = n_unknown, n_directional = n_direct,
    n_correct = n_correct, n_incorrect = n_incorrect,
    u_obs = if (n_total > 0) n_unknown / n_total else NA_real_,
    c_obs = if (n_direct > 0) n_correct / n_direct else NA_real_)
}

# Preserved within-class arm shares. These are what reweighting holds fixed, so
# it cannot manufacture an arm imbalance the trial never had.
#   r = placebo's share of the CORRECT-guess class     PLPL / (PLPL + ACAC)
#   s = active's  share of the INCORRECT-guess class   ACPL / (ACPL + PLAC)
#   t = active's  share of the UNKNOWN-response class   ACU  / (ACU  + PLU)
# A share is NA when its class is empty (0/0); the estimability check turns the
# NAs that actually matter into a clear error.
cgr_unknown_ratios <- function(st) {
  n <- lengths(st)
  share <- function(a, b) if ((a + b) > 0) a / (a + b) else NA_real_
  list(r = share(n[["PLPL"]], n[["ACAC"]]),
       s = share(n[["ACPL"]], n[["PLAC"]]),
       t = share(n[["ACU"]],  n[["PLU"]]))
}

# The six weights at target directional CGR c and target UNKNOWN rate u.
#   correct  : w_ACAC + w_PLPL = (1 - u) * c
#   incorrect: w_ACPL + w_PLAC = (1 - u) * (1 - c)
#   unknown  : w_ACU  + w_PLU  = u
# and all six sum to 1. When a class is empty its share is NA; but that class's
# total target mass is then also structurally forced to zero here only if the
# caller has confirmed estimability. To keep 0 * NA = 0 (an empty cell is never
# weighted), NA shares are treated as 0 when their class mass is 0.
cgr_unknown_weights <- function(c, u, r, s, t) {
  z <- function(v) if (is.na(v)) 0 else v          # empty-class share -> 0
  r <- z(r); s <- z(s); t <- z(t)
  c(ACAC = (1 - u) * c       * (1 - r),
    ACPL = (1 - u) * (1 - c) * s,
    ACU  = u * t,
    PLAC = (1 - u) * (1 - c) * (1 - s),
    PLPL = (1 - u) * c       * r,
    PLU  = u * (1 - t))[UNKNOWN_STRATA]
}

# Delta(c, u): reweighted active mean minus reweighted placebo mean. `mu` is a
# named list over UNKNOWN_STRATA of either six scalars (point estimate) or six
# equal-length draw vectors. A structurally zero-weight cell contributes exactly
# nothing and its mu is never touched (so an empty cell may carry NA mu safely).
cgr_unknown_delta <- function(c, u, mu, r, s, t) {
  w <- cgr_unknown_weights(c, u, r, s, t)
  term <- function(nm) if (w[[nm]] == 0) 0 else w[[nm]] * mu[[nm]]
  den_ac <- w[["ACAC"]] + w[["ACPL"]] + w[["ACU"]]
  den_pl <- w[["PLAC"]] + w[["PLPL"]] + w[["PLU"]]
  if (den_ac <= 0 || den_pl <= 0) {
    stop(sprintf("degenerate weights at c = %g, u = %g (an arm has no mass)", c, u),
         call. = FALSE)
  }
  (term("ACAC") + term("ACPL") + term("ACU")) / den_ac -
    (term("PLAC") + term("PLPL") + term("PLU")) / den_pl
}

# Estimability of the UNKNOWN-preserving estimand at target UNKNOWN rate u over a
# grid of directional CGR values. Distinguishes three states, per the brief:
#   - undefined  (error): a required directional class is empty, or u > 0 with an
#                         empty UNKNOWN class;
#   - fragile    (warn) : defined but the smallest weighted cell is thin;
#   - well populated    : otherwise.
# Returns an invisible classification; raises an error for the undefined case.
cgr_unknown_estimable <- function(st, u_target, grid, thin = 5L, warn = TRUE) {
  o <- cgr_unknown_observed(st)
  n <- o$counts
  interior <- any(grid > 0 & grid < 1) || any(abs(grid - 0.5) < 1e-12)
  if (o$n_correct == 0)
    stop("undefined estimand: the correct-guess class (ACAC + PLPL) is empty.",
         call. = FALSE)
  if (interior && o$n_incorrect == 0)
    stop("undefined estimand: the incorrect-guess class (ACPL + PLAC) is empty, ",
         "so Delta is undefined at interior c (including c = 0.50).", call. = FALSE)
  if (u_target > 0 && o$n_unknown == 0)
    stop("undefined estimand: target UNKNOWN rate u = ", signif(u_target, 4),
         " > 0 but no UNKNOWN responses were observed.", call. = FALSE)
  # cells that can carry nonzero weight somewhere on the grid
  used <- names(n)[n > 0]
  thin_cells <- used[n[used] < thin]
  state <- if (length(thin_cells)) "fragile" else "well_populated"
  if (warn && length(thin_cells))
    warning("sparse UNKNOWN-preserving strata (n < ", thin, "): ",
            paste(sprintf("%s=%d", thin_cells, n[thin_cells]), collapse = ", "),
            "; the reweighted estimate is defined but fragile.", call. = FALSE)
  invisible(list(state = state, thin_cells = thin_cells,
                 min_stratum = min(n[used])))
}

# Conjugate Normal-Inverse-Gamma posterior for the six-stratum contrast. Reuses
# nig_draws() (no duplicate posterior code). Draws are generated only for cells
# with observations; a structurally zero-weight (empty) cell is given NA mu and
# is never referenced. n_draws controls Monte Carlo precision only - it is not a
# sample size.
cgr_unknown_conjugate <- function(df, grid = seq(0, 1, length.out = 101),
                                  u_target = NULL, n_draws = 20000,
                                  prior = list(), direction = 1) {
  st <- cgr_unknown_strata(df)
  o  <- cgr_unknown_observed(st)
  u  <- if (is.null(u_target)) o$u_obs else u_target
  cgr_unknown_estimable(st, u, grid)
  rat <- cgr_unknown_ratios(st)
  mu <- stats::setNames(lapply(UNKNOWN_STRATA, function(nm) {
    y <- st[[nm]]
    if (length(y) == 0) NA_real_
    else do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior))
  }), UNKNOWN_STRATA)

  d   <- lapply(grid, cgr_unknown_delta, u = u, mu = mu,
                r = rat$r, s = rat$s, t = rat$t)
  sdv <- vapply(d, stats::sd, numeric(1))
  out <- data.frame(
    cgr = grid, method = "conjugate",
    est = vapply(d, mean, numeric(1)), sd = sdv,
    lo = vapply(d, function(x) unname(stats::quantile(x, 0.025)), numeric(1)),
    hi = vapply(d, function(x) unname(stats::quantile(x, 0.975)), numeric(1)),
    p_fav = vapply(d, function(x) mean(direction * x > 0), numeric(1)),
    stringsAsFactors = FALSE)
  out$mcse <- sdv / sqrt(n_draws)
  attr(out, "u") <- u
  out
}

# Summary at the two directional CGRs that matter: the observed directional CGR
# (the reweighting is a no-op there) and 0.50 (directional guessing at chance,
# UNKNOWN rate held fixed). Deliberately mirrors cgr_summary_table() but names
# the quantities for the UNKNOWN extension.
cgr_unknown_summary_table <- function(cur, obs_cgr, u, tol = 1e-6) {
  at <- function(target) {
    i <- which.min(abs(cur$cgr - target))
    if (abs(cur$cgr[i] - target) > tol)
      warning(sprintf(paste0("cgr_unknown_summary_table: nearest grid point %.4f ",
        "is %.2e from %.4f; include the target in the grid."),
        cur$cgr[i], abs(cur$cgr[i] - target), target), call. = FALSE)
    cur[i, ]
  }
  a <- at(obs_cgr); h <- at(0.5)
  unadj_distinct <- !(a$lo <= 0 && a$hi >= 0)
  pct <- if (unadj_distinct) 100 * (a$est - h$est) / a$est else NA_real_
  data.frame(
    directional_cgr = c(round(obs_cgr, 4), 0.5),
    unknown_rate = round(c(u, u), 4),
    what = c("observed (unadjusted)",
             "directional CGR 0.50 (UNKNOWN rate held fixed)"),
    post_mean = round(c(a$est, h$est), 3),
    cri_lo = round(c(a$lo, h$lo), 3),
    cri_hi = round(c(a$hi, h$hi), 3),
    p_favourable = round(c(a$p_fav, h$p_fav), 3),
    abs_attenuation = round(c(NA, a$est - h$est), 3),
    pct_attenuation = round(c(NA, pct), 1),
    stringsAsFactors = FALSE)
}

# The reference-line diagnostic for the UNKNOWN extension. At the observed
# directional CGR and observed UNKNOWN rate the reweighting is a no-op, so Delta
# there equals the raw active-minus-placebo mean difference. That makes the
# observed-value identity checkable exactly as in the binary case.
cgr_unknown_reference_line_test <- function(df, orig_cgr = NULL,
                                            published_unadj = NA_real_) {
  st  <- cgr_unknown_strata(df)
  o   <- cgr_unknown_observed(st)
  cgr_unknown_estimable(st, o$u_obs, c(o$c_obs, 0.5), warn = FALSE)
  rat <- cgr_unknown_ratios(st)
  mu  <- stats::setNames(lapply(UNKNOWN_STRATA, function(nm)
    if (length(st[[nm]])) mean(st[[nm]]) else NA_real_), UNKNOWN_STRATA)
  if (is.null(orig_cgr)) orig_cgr <- o$c_obs
  raw <- mean(df$value[df$condition == "AC"]) -
         mean(df$value[df$condition == "PL"])
  d_obs  <- cgr_unknown_delta(o$c_obs,  o$u_obs, mu, rat$r, rat$s, rat$t)
  d_orig <- cgr_unknown_delta(orig_cgr, o$u_obs, mu, rat$r, rat$s, rat$t)
  data.frame(
    computed_obs_directional_cgr = o$c_obs,
    observed_unknown_rate = o$u_obs,
    D_at_obs = d_obs, raw_mean_diff = raw, published_unadj = published_unadj,
    orig_cgr = orig_cgr, D_at_orig_cgr = d_orig,
    err_at_obs = d_obs - published_unadj, err_at_orig = d_orig - published_unadj,
    stringsAsFactors = FALSE)
}

# ROPE decomposition of the UNKNOWN-preserving contrast, using the same
# Delta(c, u) posterior draws. Regions relative to the band [-delta, +delta] are
# exhaustive and sum to 1 at every directional CGR.
cgr_unknown_rope <- function(df, grid = seq(0, 1, length.out = 101),
                             u_target = NULL, n_draws = 20000, delta = NULL,
                             delta_sd_frac = 0.1, direction = 1, prior = list()) {
  st <- cgr_unknown_strata(df)
  o  <- cgr_unknown_observed(st)
  u  <- if (is.null(u_target)) o$u_obs else u_target
  cgr_unknown_estimable(st, u, grid)
  rat <- cgr_unknown_ratios(st)
  if (is.null(delta)) delta <- delta_sd_frac * stats::sd(df$value)
  mu <- stats::setNames(lapply(UNKNOWN_STRATA, function(nm) {
    y <- st[[nm]]; if (length(y) == 0) NA_real_
    else do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior))
  }), UNKNOWN_STRATA)
  d <- lapply(grid, cgr_unknown_delta, u = u, mu = mu, r = rat$r, s = rat$s, t = rat$t)
  q <- function(x, p) unname(stats::quantile(x, p))
  data.frame(
    cgr = grid, delta = delta,
    est  = vapply(d, mean, numeric(1)),
    lo95 = vapply(d, function(x) q(x, 0.025), numeric(1)),
    hi95 = vapply(d, function(x) q(x, 0.975), numeric(1)),
    lo50 = vapply(d, function(x) q(x, 0.250), numeric(1)),
    hi50 = vapply(d, function(x) q(x, 0.750), numeric(1)),
    p_harm       = vapply(d, function(x) mean(direction * x < -delta), numeric(1)),
    p_negligible = vapply(d, function(x) mean(abs(x) <= delta),        numeric(1)),
    p_benefit    = vapply(d, function(x) mean(direction * x >  delta), numeric(1)),
    stringsAsFactors = FALSE)
}

# How the UNKNOWN-preserving ROPE conclusion at a single directional CGR moves as
# the band width delta is varied.
cgr_unknown_rope_sensitivity <- function(df, at_cgr = 0.5, u_target = NULL,
                                         fracs = c(0.05, 0.10, 0.20, 0.30),
                                         n_draws = 12000, direction = 1) {
  sdy <- stats::sd(df$value)
  do.call(rbind, lapply(fracs, function(f) {
    z <- cgr_unknown_rope(df, grid = at_cgr, u_target = u_target,
                          n_draws = n_draws, delta_sd_frac = f, direction = direction)
    data.frame(delta_in_SD = f, delta = f * sdy,
               p_negligible = z$p_negligible, p_benefit = z$p_benefit,
               stringsAsFactors = FALSE)
  }))
}

# cgrc_unknown(df): the one-call UNKNOWN-preserving analysis. `df` needs columns
# condition, guess and value. condition normalises to AC/PL; guess normalises to
# AC/PL/UNKNOWN (unknown_level names the raw token that means UNKNOWN, in
# addition to the recognised synonyms). unknown_rate = NULL holds u at the
# observed UNKNOWN-response rate; a numeric in [0, 1) is a sensitivity target.
# direction = +1 if higher is better, -1 if lower is better. Returns an S3
# "cgrc_unknown" object (NOT inheriting "cgrc": the binary plot/print methods
# would mislabel the directional x-axis).
cgrc_unknown <- function(df, unknown_level = "UNKNOWN", unknown_rate = NULL,
                         n_draws = 20000, direction = 1, prior = list()) {
  if (!is.null(unknown_rate) && (unknown_rate < 0 || unknown_rate >= 1))
    stop("unknown_rate must be NULL or in [0, 1).", call. = FALSE)
  trial <- data.frame(
    condition = cgrc_normalise_arm(df$condition, "treatment received"),
    guess     = cgrc_normalise_guess(df$guess, allow_unknown = TRUE,
                                     unknown_labels = unknown_level),
    value     = as.numeric(df$value), stringsAsFactors = FALSE)
  # genuinely missing rows (NA condition/guess/value) are dropped as incomplete;
  # an observed UNKNOWN is NOT missing and is retained.
  trial <- trial[stats::complete.cases(trial), ]

  st <- cgr_unknown_strata(trial)
  o  <- cgr_unknown_observed(st)
  u  <- if (is.null(unknown_rate)) o$u_obs else unknown_rate
  grid <- sort(unique(c(seq(0, 1, length.out = 101), o$c_obs, 0.5)))
  cur  <- cgr_unknown_conjugate(trial, grid, u_target = u, n_draws = n_draws,
                                direction = direction, prior = prior)
  structure(list(
    curve = cur,
    summary = cgr_unknown_summary_table(cur, o$c_obs, u),
    observed_directional_cgr = o$c_obs,
    observed_unknown_rate = o$u_obs,
    target_unknown_rate = u,
    counts = o$counts,
    n_total = o$n_total, n_directional = o$n_directional, n_unknown = o$n_unknown,
    direction = direction,
    method = "UNKNOWN-preserving CGRC extension"),
    class = "cgrc_unknown")
}

print.cgrc_unknown <- function(x, ...) {
  cat(sprintf("%s\n", x$method))
  cat(sprintf("observed directional CGR = %.4f; UNKNOWN rate = %.1f%% held at %.1f%%\n",
              x$observed_directional_cgr, 100 * x$observed_unknown_rate,
              100 * x$target_unknown_rate))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

plot.cgrc_unknown <- function(x, ...) {
  cgr_unknown_plot(x$curve, obs_cgr = x$observed_directional_cgr,
                   u = x$target_unknown_rate, ...)
}

# cgrc_unknown_headline(df): the interpretable summary for the UNKNOWN extension.
# Reports the two plain probabilities before/after reweighting directional guesses
# to 0.50 (UNKNOWN rate held fixed), and states plainly that this is an extension,
# not the original Szigeti estimand, using reweighting rather than causal wording.
cgrc_unknown_headline <- function(df, unknown_level = "UNKNOWN", unknown_rate = NULL,
                                  direction = 1, delta = NULL, delta_sd_frac = 0.5,
                                  n_draws = 20000, prior = list()) {
  trial <- data.frame(
    condition = cgrc_normalise_arm(df$condition, "treatment received"),
    guess     = cgrc_normalise_guess(df$guess, allow_unknown = TRUE,
                                     unknown_labels = unknown_level),
    value     = as.numeric(df$value), stringsAsFactors = FALSE)
  trial <- trial[stats::complete.cases(trial), ]
  st  <- cgr_unknown_strata(trial)
  o   <- cgr_unknown_observed(st)
  u   <- if (is.null(unknown_rate)) o$u_obs else unknown_rate
  cgr_unknown_estimable(st, u, c(o$c_obs, 0.5))
  rat <- cgr_unknown_ratios(st)
  if (is.null(delta)) delta <- delta_sd_frac * stats::sd(trial$value)
  mu <- stats::setNames(lapply(UNKNOWN_STRATA, function(nm) {
    y <- st[[nm]]; if (length(y) == 0) NA_real_
    else do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior))
  }), UNKNOWN_STRATA)
  d_obs <- cgr_unknown_delta(o$c_obs, u, mu, rat$r, rat$s, rat$t)
  d_bl  <- cgr_unknown_delta(0.5,     u, mu, rat$r, rat$s, rat$t)
  eff_obs <- direction * d_obs; eff_bl <- direction * d_bl
  q <- function(x, p) unname(stats::quantile(x, p))
  res <- list(
    observed_directional_cgr = o$c_obs, observed_unknown_rate = o$u_obs,
    target_unknown_rate = u, delta = delta, direction = direction,
    n_total = o$n_total, n_directional = o$n_directional, n_unknown = o$n_unknown,
    counts = o$counts,
    adj_est = mean(d_bl), adj_lo = q(d_bl, 0.025), adj_hi = q(d_bl, 0.975),
    p_dir_obs = mean(eff_obs > 0), p_dir_blind = mean(eff_bl > 0),
    p_meaningful_obs = mean(eff_obs > delta), p_meaningful_blind = mean(eff_bl > delta))
  dtxt <- trimws(sprintf("%.2g", delta))
  unit <- if (identical(dtxt, "1")) "point" else "points"
  res$text <- sprintf(paste0(
    "UNKNOWN-preserving extension (not the original Szigeti estimand). Raw, the ",
    "directional correct-guess rate was %.1f%%, and %.1f%% of participants ",
    "answered UNKNOWN. Reweighting directional guesses to 50%% while holding the ",
    "UNKNOWN-response rate at %.1f%% changed the estimated effect from %.2f to ",
    "%.2f (95%% CrI %.2f to %.2f). The probability of a favourable effect went ",
    "from %.0f%% to %.0f%%, and of a meaningful effect (beyond %s %s) from %.0f%% ",
    "to %.0f%%. This is a reweighted estimate under the specified guess ",
    "distribution, not an expectancy-removed causal effect."),
    100 * o$c_obs, 100 * o$u_obs, 100 * u, mean(d_obs), mean(d_bl),
    res$adj_lo, res$adj_hi, 100 * res$p_dir_obs, 100 * res$p_dir_blind,
    dtxt, unit, 100 * res$p_meaningful_obs, 100 * res$p_meaningful_blind)
  structure(res, class = "cgrc_unknown_headline")
}

print.cgrc_unknown_headline <- function(x, ...) { cat(x$text, "\n"); invisible(x) }
