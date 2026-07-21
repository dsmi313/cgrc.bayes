# Unresolved discrepancies

Ranked by how much they would change a conclusion. Do NOT close any of these by
choosing the convenient answer.

> **Update (2026-07-21).** The three Szigeti source papers were obtained and
> read (see `reports/SOURCES.md`). U1 is resolved (JAGS ran, PASS). U3, U4 and
> U6 are now grounded in the papers rather than in recollection, and U2's
> methodology is confirmed from the paper. Where an item is genuinely closed it
> says so; where the paper narrows but does not settle it, the residual is
> stated plainly. Nothing was closed by picking a convenient answer.

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

## U2  n = 232 versus published n = 233   [methodology confirmed; off-by-one open]
**Source check (Szigeti 2023, Sci Rep 13:12107 — see SOURCES.md S1).** The paper
states the analysis rule verbatim: "we only used data from the first week of the
experiment, thus, each datapoint is independent... In the current analysis
**n = 233** datapoints were included." That is exactly the `tp == "w1s1"` filter
used here, so the analysis approach is confirmed identical - this is not a
filtering mistake on our side. The public file nonetheless yields 232 on every
acute scale (no duplicate IDs, no missing values). The discrepancy is a genuine
one-record gap between the published n and the public release.
**Still open:** the papers do not account for the extra record. Most likely a
record withheld from the public file or a typo in the paper. Confirm with the
author.

## U3  The quoted CGR of 0.72   [explained; provenance worth one confirmation]
**Source check (Szigeti 2023 Fig 4 caption — S1; Szigeti & Heifets 2024 — S3).**
The 2023 paper does quote it: the Fig 4 caption reads "vertical green dashed line
corresponds to the trial's original CGR (= 0.72)", and the line is drawn at
~0.72. But 0.72 cannot be the analyzed data's overall CGR:
- The public week-1 data gives overall CGR **0.647** and reproduces the paper's
  own Table 2 estimates exactly (both are independent of where the CGR line sits,
  so the analyzed data's observed CGR is ~0.647, not 0.72).
- The **placebo-arm** correct-guess rate is **0.7234 ~= 0.72**.
- The same author's 2024 review states the microdose correct-guess rate is
  "only ~65% to 70%" - again not 0.72.
**Conclusion:** 0.72 is the placebo-conditional correct-guess rate, not the
overall CGR; the Fig 4 reference line is therefore misplaced by ~0.07 on the CGR
axis - the exact class of error the observed-CGR identity check (Section 4)
exists to catch. Numerically compelling and now backed by the sources; the only
thing left is a one-line confirmation from the author that 0.72 was the
placebo-conditional figure.

## U4  Szigeti's generative source code was not located
`szb37/mcrds_public` serves `data/pacutes.csv` but returns 404 for `README.md`
and every source path probed; the GitHub API was rate-limited. The `noise`
resolution (CH-04) rests on convergent empirical evidence, not on reading the
author's implementation.
**Source check (Szigeti 2023 Table 1 — S1).** The published Table 1 rates
(0.05 / 0.86 / 0.78 / 0.99) and the DTE effect size (Hedges g = 0.4), plus the
stated regime "n ~ 200, CGR ~ 0.7, treatment effect ~ 0.4 Hedges' g", are now
confirmed against the actual paper and all match `noise = "all"` (0.4011), not
`noise = "arm"` (0.50). This upgrades the resolution from "recollection" to
"matches the published operating characteristics".
**Still open:** the exact Eq. 4 SD scope lives in the paper's Supplementary
Table 1, which was not in the supplied PDFs, so this is not yet a line-by-line
code check. Ask the author for the AEB simulation code to close it fully.

## U5  No head-to-head operating-characteristic comparison
Section 9 characterises the adjusted ESTIMAND under the AEB model. It does not
compare the KDE procedure and the Bayesian procedure on bias, RMSE, coverage,
and power. Section 8 shows they converge to the same point estimate, so any
difference would be in interval behaviour and Monte Carlo stability.
**Do not claim the Bayesian method performs better until this is run.**

## U6  Benign versus malicious unblinding is not testable here
If people guess correctly BECAUSE they improved, CGR adjustment removes a real
effect.
**Source check (Szigeti 2023 — S1; Szigeti & Heifets 2024 — S3).** The paper's
argument for malicious unblinding is now sourced exactly: 55% of participants
named "body/perceptual sensations" as the primary cue (muscle tension 58%,
stomach discomfort 27%) versus only 23% naming "mental/psychological benefits";
and the placebo-microdose PANAS positive/negative difference (2.1/0.8) is
~500-750% smaller than the natural within-subject day-to-day variability
(~10/~6). The 2024 review adds that microdose effects fall below the 0.5-SMD
minimally-important difference - "too small to be noticeable".
**Still open by nature:** that is an argument from evidence, not a demonstration.
No analysis of these data can settle whether unblinding is benign or malicious;
this item stays open on principle, not for lack of a number.

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
