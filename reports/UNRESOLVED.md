# Unresolved discrepancies

Ranked by how much they would change a conclusion. Do NOT close any of these by
choosing the convenient answer.

> **Update (2026-07-21).** The three Szigeti source papers *and the author's
> analysis source code* (`szb37/CorrectGuessRateCurve`) were obtained and read
> (see `reports/SOURCES.md`). Now resolved: U1 (JAGS ran, PASS), U3 (0.72 is a
> hardcoded constant in the code), U4 (source located; estimand independently
> confirmed). U6 is grounded in the papers; U2 remains open. A new item U9
> records a guess-rate pattern that contradicts the 2024 review. Where an item
> is genuinely closed it says so; residuals are stated plainly. Nothing was
> closed by picking a convenient answer.

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

## U3  The quoted CGR of 0.72   [RESOLVED 2026-07-21 - hardcoded constant]
**Settled from the author's source code (SOURCES.md S4).**
`szb37/CorrectGuessRateCurve/src/config.py` contains
`trial_cgrs = {'sbmd': 0.72}` - a **hardcoded constant**, drawn as the Figure 4
reference line. It is not computed from the data: the same code computes the
trial CGR as `(n_plpl + n_acac) / n` = 0.647 and does not use it for the line.
So 0.72 is a fixed annotation, not the data's correct guess rate; it coincides
with the **placebo-arm** correct-guess rate (0.7234). The reference-line test
(Rmd Section 8) shows the curve equals the reported unadjusted values at 0.647,
not at 0.72. This is no longer a hypothesis.
**Residual (cosmetic):** whether the author intended 0.72 as the placebo-arm
figure or simply mis-set the constant is a question only they can answer, but it
does not affect any estimate here.

## U4  Szigeti's generative source code   [RESOLVED 2026-07-21 - located]
The code was named in the 2023 paper's data-availability statement all along, at
`github.com/szb37/CorrectGuessRateCurve`; earlier turns probed the *data* mirror
`szb37/mcrds_public` instead and wrongly reported it not located.
**Confirmed from source (SOURCES.md S4):** `get_strata_ratio` /
`get_strata_sample_sizes` form `r = PLPL/(PLPL+ACAC)` and `s = ACPL/(ACPL+PLAC)`,
identical to this implementation - independent confirmation of the estimand.
The default `strata_sampling = 'all_prop'` confirms the `noise = "all"` reading
(CH-04), matching published Table 1 (0.05/0.86/0.78/0.99) and Hedges g = 0.4.
Two documentation corrections followed (both in CHANGELOG CH-16): the resample
count is 32/13-grid, not 100; and `get_strata_ratio` does `round(x, 2)`, so
`legacy_round = TRUE` is the faithful reproduction path.
**Residual:** the exact AEB-simulation SD scope (Eq. 4) was not re-derived
line-by-line from the simulation module; the operating-characteristic match is
strong enough that this is low priority.

## U9  Guess rates are reversed vs the author's stated pattern   [for the author]
Szigeti & Heifets (2024, S3) state correct guess rates are "generally higher in
the active arms". In this dataset the pattern is the **reverse**: the placebo
arm guesses correctly at 0.723 and the active (microdose) arm at only 0.528
(Rmd Section 10). It is driven by a response bias toward guessing "placebo"
(9508 placebo vs 6115 microdose guesses across the dataset), so participants on
placebo are "right" far more often than those on microdose. This does not break
the estimand - the CGRC conditions on the observed strata through r and s - but
it contradicts the generalisation in the review and is worth raising with the
author, because it changes the intuition for which arm expectancy inflates.

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

## U10  UNKNOWN-preserving extension: assumptions   [operating characteristics now simulated; assumptions remain]
The six-stratum UNKNOWN extension (CH-17) is an addition by this package, not the
Szigeti estimand. Two interpretive assumptions remain open by nature, and the
one empirical gap that was open has now been addressed (CH-23):
- **The preserved UNKNOWN arm share `t` is an assumption.** Reweighting holds
  `t = ACU/(ACU+PLU)` fixed as it scales the UNKNOWN mass, exactly as the CGRC
  holds r and s fixed. This preserves the observed within-UNKNOWN arm ratio; it
  does NOT prove that treatment assignment is independent of the UNKNOWN
  response. The independent-guess-distribution estimand (`cgr_unknown_independent`,
  CH-20) targets a different assumption and is offered as a contrast, not a
  resolution.
- **`c = 0.50` is not "perfect blinding".** It is directional guessing at chance
  while the UNKNOWN rate is held fixed. Whether an UNKNOWN response is itself a
  sign of good blinding, or of disengagement, is not decidable from the counts.

**Operating characteristics — now simulated [ADDRESSED 2026-07-24, CH-23].** A
purpose-built UNKNOWN-aware generative model (`sim_aeb_unknown`) and an operating-
characteristics study (`cgr_unknown_operating`) now exist. Under the model's two
explicit assumptions — (A1) the UNKNOWN-response rate is equal in both arms, and
(A2) an UNKNOWN responder carries no expectancy (no directional belief -> no
PT->TE path) — the six-stratum estimator of `Delta(0.50, u_obs)` is essentially
unbiased for the direct effect (|bias| < 0.06 at n=300) with 95% coverage
(~0.95-0.96 across all four DTE x AEB scenarios), and under pure expectancy the
adjusted false-favourable rate stays ~0.05 while the naive t-test flags ~0.73.
The app's Panel B runs this on demand at the uploaded trial's n, directional CGR
and observed UNKNOWN rate.
A precomputed UNKNOWN design lookup now exists too (CH-24,
`inst/extdata/cgrc_unknown_lookup.rds`, over n x p_cg x true_effect x u at
mu_aeb = 7.7), and Panel A's UNKNOWN-response-rate slider reads it, so the design
tool covers UNKNOWN designs interactively - not only the on-demand Panel B run.
**Still conditional:** this is validation *under A1 and A2*. It is not a claim
that the extension outperforms the binary method in general, nor that A1/A2 hold
in any real trial. Differential UNKNOWN rates (A1 relaxed) and UNKNOWN responders
who do carry expectancy (A2 relaxed) are not yet characterised; the six strata
are also thinner than four, so `empty_stratum_rate` climbs faster at small n /
high guess or UNKNOWN rates. The UNKNOWN lookup fixes mu_aeb at 7.7 (the adjusted
estimator is ~insensitive to it under A2; the unadjusted comparator is not), so
the mu_aeb control does not move the UNKNOWN design panel.

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
