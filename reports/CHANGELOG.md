# Change log

Every substantive modification, with rationale. Nothing here was applied
silently. Each item can be reverted independently.

## Statistical corrections

### CH-01  JAGS prior corrected  [BREAKING - affects results if verification runs]
**Was:** `mu[j] ~ dnorm(0, 1.0E-6)`; `tau[j] ~ dgamma(1.0E-3, 1.0E-3)`
**Now:** `tau[j] ~ dgamma(a0, b0)`; `mu[j] ~ dnorm(m0, k0 * tau[j])`
**Why:** Normal-Inverse-Gamma requires Var(mu | sigma2) = sigma2/k0, i.e.
precision k0*tau. The old code fixed the prior precision at 1e-6 regardless of
tau, which is a different model. The prior claim that both backends target an
identical posterior was therefore FALSE AS WRITTEN.
**Status:** VERIFIED 2026-07-21 (JAGS 4.3.2). `cgr_check_backends()` ran over
the full grid and PASSED: posterior means agree to <= 0.006, max |z| = 0.63,
max Rhat = 1.0001, min ESS = 39664. The identity claim is now supported. The
standing instruction is unchanged for any future regression: if it ever fails,
delete the claim - do not weaken the conjugate prior to force agreement.

### CH-02  Two-sided posterior tail curve removed  [BREAKING - changes figures]
**Was:** `p_two = 2 * pmin(p_pos, 1 - p_pos)` plotted with a magenta line at 0.05.
**Now:** `p_fav = P(direction * Delta > 0 | y)`, labelled "Posterior probability
the effect is [favourable]". Optional line at 0.95, explicitly described as
descriptive rather than a universal cutoff.
**Why:** the old curve invited reading a posterior as a p-value. It was
constructed to look like one and placed on a 0.05 threshold.
**Added:** `direction` argument so outcomes where lower is better (QIDS, STAIT)
declare favourability rather than assuming positive = good.

### CH-03  `S` renamed `n_draws`
**Why:** `S` was ambiguous against participant sample size. Text now states
that draws reduce Monte Carlo error in posterior summaries only and carry no
information about participants.

### CH-04  AEB `noise` ambiguity resolved to `"all"` on evidence
**Was:** unresolved; document claimed neither variant reproduced both published
quantities.
**Now:** `"all"` is the documented default.
**Why:** that earlier claim was WRONG, and wrong because it compared a posterior
probability against a frequentist p-value. Proper simulation gives:
  - Hedges g: "all" 0.4011 vs published 0.40; "arm" 0.5023
  - Table 1 unadjusted significance rates:
      published  0.05 / 0.86 / 0.78 / 0.99
      "all"      0.056/ 0.876/ 0.826/ 0.992
      "arm"      0.056/ 0.970/ 0.930/ 0.998
**Caveat:** author's source code NOT located (repo 404s on all source paths).
Convergent empirical evidence, not verification. Figure 3 is labelled a
reproduction *consistent with* the published operating characteristics.

### CH-05  Single-trial illustration replaced by 500-trial operating characteristics
**Why:** one simulated trial is an illustration, not validation. Now reports
bias, RMSE, 95% coverage, P(posterior > 0.95), and frequentist significance
rate for all four AEB scenarios. Coverage is 0.936-0.960 against nominal 0.95 -
the first direct evidence the intervals are calibrated.

### CH-06  Original KDE procedure implemented and run at 100/1000/10000
**Why:** required to separate three confounded sources of difference. Result:
KDE converges to the analytic value at 10 000 resamples in every cell, so the
KDE-vs-Gaussian choice contributes essentially nothing; but 100 resamples
carries 0.12-0.29 points of Monte Carlo SE. The averaged p-values reproduce the
published ones (PANAS 0.41 vs 0.43, Energy 0.043 vs 0.04), confirming the port
is faithful.
**Consequence for framing:** the contribution is computational stability plus a
change of inferential summary - NOT a change of estimand, and not better
inference about the trial.

### CH-07  Comparative language removed
"Better", "superior", and "improves on" do not appear. No comparative claim
appears before the simulation study. Section 13 states narrowly what the
approach does and does not add.

### CH-08  Optional Beta uncertainty in r and s  [OFF BY DEFAULT]
`cgr_ratio_draws()`. Off by default because the original estimand CONDITIONS on
the observed within-class ratios. Labelled an extension, not a reproduction.

### CH-09  Data pinned locally with provenance
`data-raw/download_data.R` records URL, access date, SHA-256, size, and record
count; the Rmd verifies the hash at render time and stops on mismatch.
No network access during rendering.

### CH-10  Implementation moved to `R/`
Six modules sourced by the Rmd, which is now a vignette rather than the only
copy of the code.

### CH-11  All dplyr verbs namespaced / removed
`R/` is base R only. Cause: `dplyr::filter` and `stats::filter` collide, and
`stats::filter` has a `method` argument, so a masked
`filter(x, method == "conjugate")` fails with "object 'method' not found".
This was an actual render failure, not a hypothetical.

## Findings that changed the document's claims

### F-A  The 0.72 CGR is probably the placebo-arm rate
Overall week-1 CGR is 0.6466. Pooling timepoints does not reach 0.72 (w1 0.647,
w2 0.679, w3 0.620, w4 0.651, pooled 0.649; participant-level mean 0.6468).
But the correct-guess rate WITHIN the placebo arm is 0.7234, and within the
microdose arm only 0.5275. Hypothesis with strong numerical support; NOT
confirmed. For the author.

### F-B  n = 232 vs published n = 233 remains unexplained
No scale at any timepoint has n = 233. Eight scales at w1s1 have exactly 232.
No duplicate trial_ids, no missing outcome values. Off by one record. For the
author.

## Preserved unchanged

- The estimand: strata, r, s, weights, Delta(c).
- The observed-CGR identity check.
- Default vague priors (m0=0, k0=1e-6, a0=b0=1e-3).
- Week-1-only filter and use of raw value as outcome.
- The supplementary stratum-allocation validation.
- The published empirical results: PANAS 3.16 -> 1.08 etc. are unchanged.
