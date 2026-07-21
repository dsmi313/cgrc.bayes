# Core CGRC estimand. Base R only: dplyr::filter and stats::filter collide,
# and stats::filter has a `method` argument, so a masked call fails with the
# baffling error "object 'method' not found".

STRATA <- c("ACAC", "ACPL", "PLAC", "PLPL")

cgr_strata <- function(df) {
  key <- factor(paste0(df$condition, df$guess), levels = STRATA)
  if (anyNA(key)) stop("condition/guess must be coded AC/PL", call. = FALSE)
  out <- split(df$value, key)
  empty <- STRATA[lengths(out) == 0L]
  if (length(empty)) {
    stop("empty stratum: ", paste(empty, collapse = ", "),
         "; the estimand is undefined", call. = FALSE)
  }
  out
}

cgr_observed <- function(st) {
  (length(st$ACAC) + length(st$PLPL)) / sum(lengths(st))
}

# r = placebo's share of the CORRECT-guess mass
# s = active's  share of the INCORRECT-guess mass
# Reweighting preserves these, so it cannot manufacture an arm imbalance the
# trial never had.
cgr_ratios <- function(st, legacy_round = FALSE) {
  n <- sum(lengths(st)); rho <- lengths(st) / n
  if (legacy_round) rho <- round(rho, 2)   # reproduces Szigeti's shipped code
  list(r = unname(rho[["PLPL"]] / (rho[["PLPL"]] + rho[["ACAC"]])),
       s = unname(rho[["ACPL"]] / (rho[["ACPL"]] + rho[["PLAC"]])))
}

# w_ACAC = c(1-r) ; w_PLPL = cr ; w_ACPL = (1-c)s ; w_PLAC = (1-c)(1-s)
# Correct-guess weights sum to c; incorrect-guess weights sum to 1-c.
cgr_weights <- function(c, r, s) {
  c(ACAC = c * (1 - r), ACPL = (1 - c) * s,
    PLAC = (1 - c) * (1 - s), PLPL = c * r)
}

# Delta(c): reweighted active mean minus reweighted placebo mean.
# mu may be four scalars (point estimate) or four equal-length draw vectors.
cgr_delta <- function(c, mu, r, s) {
  w <- cgr_weights(c, r, s)
  den_ac <- w[["ACAC"]] + w[["ACPL"]]
  den_pl <- w[["PLPL"]] + w[["PLAC"]]
  if (den_ac <= 0 || den_pl <= 0) {
    stop("degenerate weights at cgr = ", c, call. = FALSE)
  }
  (w[["ACAC"]] * mu$ACAC + w[["ACPL"]] * mu$ACPL) / den_ac -
    (w[["PLPL"]] * mu$PLPL + w[["PLAC"]] * mu$PLAC) / den_pl
}

# Stratum sizes at a target CGR. Not used by the estimand; validates against
# Szigeti's supplementary worked example (65/35/25/75 -> 54/58/42/46 at 0.5).
cgr_sizes <- function(st, cgr, legacy_round = FALSE) {
  n <- sum(lengths(st)); rat <- cgr_ratios(st, legacy_round)
  n_corr <- round(cgr * n); n_inc <- n - n_corr
  plpl <- round(n_corr * rat$r); acpl <- round(n_inc * rat$s)
  c(ACAC = n_corr - plpl, ACPL = acpl, PLAC = n_inc - acpl, PLPL = plpl)
}

cgr_analytic <- function(df, grid = seq(0, 1, length.out = 101)) {
  st <- cgr_strata(df); rat <- cgr_ratios(st); mu <- lapply(st, mean)
  data.frame(cgr = grid, method = "analytic",
             est = vapply(grid, cgr_delta, numeric(1),
                          mu = mu, r = rat$r, s = rat$s),
             stringsAsFactors = FALSE)
}

# Optional EXTENSION (not part of the original estimand, off by default).
# The original CGRC conditions on the observed within-class ratios. Treating
# them as uncertain propagates binomial sampling error in the stratum
# composition, which is a different (larger) uncertainty statement.
cgr_ratio_draws <- function(st, n_draws, alpha = 1) {
  n <- lengths(st)
  list(r = stats::rbeta(n_draws, n[["PLPL"]] + alpha, n[["ACAC"]] + alpha),
       s = stats::rbeta(n_draws, n[["ACPL"]] + alpha, n[["PLAC"]] + alpha))
}
