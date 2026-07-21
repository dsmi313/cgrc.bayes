# Unresolved discrepancies

Ranked by how much they would change a conclusion. Do NOT close any of these by
choosing the convenient answer.

## U1  Corrected JAGS prior   [RESOLVED 2026-07-21 - PASS]
No JAGS in the authoring environment, so `cgr_check_backends()` had never run.
It has now been run in this environment (JAGS 4.3.2, rjags) over the full
101-point grid: conjugate backend at 40 000 draws vs the corrected 4-chain JAGS
model at 10 000 iterations/chain.
**Result:** PASS. Posterior means agree to <= 0.006 across the grid; 95%
credible limits to <= 0.026 (lower) and <= 0.079 (upper); posterior
probabilities to <= 0.005; the largest mean difference is 0.63 Monte Carlo SE
(threshold 5). JAGS convergence is clean (max Rhat 1.0001, min ESS 39664). The
two backends target the same posterior, so Section 8's identity claim is
supported.
**Caveat:** this requires a JAGS install. A render without rjags skips the check
and the claim reverts to unverified (Section 8 prints the NOT RUN note in that
case). The instruction stands for any future failure: if it ever fails, delete
the identity claim rather than weakening the conjugate prior to force agreement.

## U2  n = 232 versus published n = 233   [needs the author]
No filtering rule over the public file reproduces 233. Eight w1s1 scales have
exactly 232; no duplicate IDs; no missing values. One record short.
**Candidates:** a record withheld from the public release; a typo; a different
inclusion rule not described in the paper.

## U3  The quoted CGR of 0.72   [needs the author, strong candidate answer]
Overall week-1 CGR is 0.647. The placebo-arm correct-guess rate is 0.7234.
**Hypothesis:** 0.72 is conditional on receiving placebo. Numerically compelling
but unconfirmed. Ask directly.

## U4  Szigeti's generative source code was not located
`szb37/mcrds_public` serves `data/pacutes.csv` but returns 404 for `README.md`
and every source path probed; the GitHub API was rate-limited. The `noise`
resolution (CH-04) rests on convergent empirical evidence, not on reading the
author's implementation.
**Resolution:** ask the author for the AEB simulation code and re-check.

## U5  No head-to-head operating-characteristic comparison
Section 9 characterises the adjusted ESTIMAND under the AEB model. It does not
compare the KDE procedure and the Bayesian procedure on bias, RMSE, coverage,
and power. Section 8 shows they converge to the same point estimate, so any
difference would be in interval behaviour and Monte Carlo stability.
**Do not claim the Bayesian method performs better until this is run.**

## U6  Benign versus malicious unblinding is not testable here
If people guess correctly BECAUSE they improved, CGR adjustment removes a real
effect. The paper argues for malicious unblinding from self-reported cues (55%
bodily sensations vs 23% mental benefits) and from effect size relative to
day-to-day variability. That is an argument, not a demonstration, and no
analysis of these data can settle it.

## U7  KDE bandwidth was never a deliberate choice
The original uses sklearn's default fixed bandwidth of 1.0, not a data-adaptive
rule. On a 0-100 VAS that is very narrow. It does not matter here because the
estimand uses only means and KDE smoothing preserves means - but it would
matter for any quantile or tail-based estimand.

## U8  Stratum normality is assumed, not established
Skew and excess kurtosis are reported per stratum in Section 11, but no formal
check or robust-likelihood comparison has been run. JAGS is now available
(see U1), and `cgr_jags(likelihood = "t")` implements the Student-t sensitivity,
but the Rmd does not exercise it at render time, so the robust-likelihood
comparison remains outstanding (TESTS.md N4).
