# CORRECTED JAGS model.
#
# The previous version used independent priors:
#     mu[j]  ~ dnorm(0, 1.0E-6)
#     tau[j] ~ dgamma(1.0E-3, 1.0E-3)
# which is NOT the Normal-Inverse-Gamma model, because NIG has
# Var(mu | sigma2) = sigma2 / k0, i.e. precision of mu equal to k0 * tau.
# Fixing the prior precision at 1e-6 regardless of tau is a different model,
# so the earlier claim that the two backends target an identical posterior
# was false as written.
#
# Corrected below. m0, k0, a0, b0 come through the data list so the two
# implementations cannot drift apart.
#
# STATUS: NOT VERIFIED. No JAGS was available when this was written. Until
# cgr_check_backends() is run and passes, the identity claim is unsupported.

jags_model_string <- function(likelihood = c("normal", "t")) {
  likelihood <- match.arg(likelihood)
  lik <- if (likelihood == "normal") {
    "  for (i in 1:N) { y[i] ~ dnorm(mu[k[i]], tau[k[i]]) }\n"
  } else {
    paste0("  for (i in 1:N) { y[i] ~ dt(mu[k[i]], tau[k[i]], nu) }\n",
           "  nu ~ dexp(0.1) T(2, 100)\n")
  }
  paste0(
    "model {\n", lik,
    "  for (j in 1:4) {\n",
    "    tau[j] ~ dgamma(a0, b0)\n",
    "    mu[j]  ~ dnorm(m0, k0 * tau[j])\n",   # <- conditional precision
    "  }\n",
    "  for (m in 1:M) {\n",
    "    w_acac[m] <- cgr[m] * (1 - r)\n",
    "    w_acpl[m] <- (1 - cgr[m]) * s\n",
    "    w_plpl[m] <- cgr[m] * r\n",
    "    w_plac[m] <- (1 - cgr[m]) * (1 - s)\n",
    "    delta[m] <- (w_acac[m]*mu[1] + w_acpl[m]*mu[2]) / (w_acac[m]+w_acpl[m]) -\n",
    "                (w_plpl[m]*mu[4] + w_plac[m]*mu[3]) / (w_plpl[m]+w_plac[m])\n",
    "  }\n}\n")
}

cgr_jags <- function(df, grid = seq(0, 1, length.out = 101),
                     likelihood = c("normal", "t"),
                     n_iter = 10000, n_burn = 2000, n_chains = 4, seed = 1,
                     prior = list(m0 = 0, k0 = 1e-6, a0 = 1e-3, b0 = 1e-3),
                     direction = 1) {

  likelihood <- match.arg(likelihood)
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("needs rjags + a JAGS install; see https://mcmc-jags.sourceforge.io",
         call. = FALSE)
  }

  st <- cgr_strata(df); rat <- cgr_ratios(st)
  eps <- 1e-6
  grid_safe <- pmin(pmax(grid, eps), 1 - eps)   # c in {0,1} can zero a weight

  y <- unlist(st, use.names = FALSE)
  k <- rep(seq_along(st), lengths(st))

  dat <- c(list(y = y, k = k, N = length(y),
                cgr = grid_safe, M = length(grid_safe),
                r = rat$r, s = rat$s), prior)

  inits <- lapply(seq_len(n_chains), function(i)
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + i))

  con <- textConnection(jags_model_string(likelihood))
  on.exit(close(con), add = TRUE)
  m <- rjags::jags.model(con, data = dat, inits = inits,
                         n.chains = n_chains, quiet = TRUE)
  stats::update(m, n_burn)
  mon <- if (likelihood == "t") c("delta", "nu") else "delta"
  samp <- rjags::coda.samples(m, variable.names = mon, n.iter = n_iter)

  all <- as.matrix(samp)
  dc <- grep("^delta\\[", colnames(all), value = TRUE)
  dc <- dc[order(as.integer(sub("^delta\\[(\\d+)\\]$", "\\1", dc)))]
  dr <- all[, dc, drop = FALSE]

  ess <- rep(NA_real_, length(dc)); rh <- ess
  if (requireNamespace("coda", quietly = TRUE)) {
    e <- try(coda::effectiveSize(samp), silent = TRUE)
    if (!inherits(e, "try-error")) ess <- unname(e[dc])
    if (n_chains > 1) {
      g <- try(coda::gelman.diag(samp, multivariate = FALSE)$psrf[, 1],
               silent = TRUE)
      if (!inherits(g, "try-error")) rh <- unname(g[dc])
    }
  }

  out <- data.frame(
    cgr = grid, method = if (likelihood == "t") "jags-t" else "jags",
    est = unname(colMeans(dr)),
    sd  = unname(apply(dr, 2, stats::sd)),
    lo  = unname(apply(dr, 2, function(x) unname(stats::quantile(x, .025)))),
    hi  = unname(apply(dr, 2, function(x) unname(stats::quantile(x, .975)))),
    p_fav = unname(apply(dr, 2, function(x) mean(direction * x > 0))),
    ess = ess, rhat = rh, stringsAsFactors = FALSE
  )
  out$mcse <- out$sd / sqrt(out$ess)   # autocorrelation-corrected
  out
}

# DECISION-1 verification. Reports max abs differences in posterior means,
# both credible limits, and posterior probabilities, as required.
cgr_check_backends <- function(df, grid = seq(0, 1, length.out = 101),
                               n_draws = 40000, jags_args = list()) {
  a <- cgr_conjugate(df, grid, n_draws = n_draws)
  b <- do.call(cgr_jags, c(list(df = df, grid = grid), jags_args))
  a <- a[order(a$cgr), ]; b <- b[order(b$cgr), ]
  mcse <- sqrt(a$mcse^2 + b$mcse^2)
  res <- list(
    max_abs_diff_mean = max(abs(a$est   - b$est)),
    max_abs_diff_lo   = max(abs(a$lo    - b$lo)),
    max_abs_diff_hi   = max(abs(a$hi    - b$hi)),
    max_abs_diff_pfav = max(abs(a$p_fav - b$p_fav)),
    max_abs_z         = max(abs((a$est - b$est) / mcse)),
    max_rhat          = max(b$rhat, na.rm = TRUE),
    min_ess           = min(b$ess,  na.rm = TRUE)
  )
  res$verdict <- if (is.finite(res$max_abs_z) && res$max_abs_z < 5) {
    "PASS - differences within Monte Carlo error; identity claim supported"
  } else {
    "FAIL - remove the identity claim; describe as comparable weak-prior models"
  }
  res
}
