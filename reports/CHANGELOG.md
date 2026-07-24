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

## Bug fixes (2026-07-21)

### CH-13  Grid-snapping in the empirical Results table  [BREAKING - corrects a number]
The results table reported the unadjusted PANAS effect at the nearest grid point
(0.65) instead of the exact observed CGR (0.6466): Delta(0.65) = 3.207 vs
Delta(0.6466) = 3.157, an 0.05 overstatement (Mood 0.09, Energy 0.10). Awkward,
because the document flags a 1.03 discrepancy in the paper's 0.72 line while
carrying its own error from the same class of mistake.
**Fix:** the empirical chunk evaluates on `sort(unique(c(GRID, c_obs)))` so the
posterior is read at exactly c_obs; `cgr_summary_table()` now warns if the grid
does not contain the requested target within tol, so snapping cannot recur
silently.

### CH-14  pct_attenuation blew up for cognitive performance
`pct_attenuation` was 151.5% for CPS - a ratio to an unadjusted estimate of
-0.011 whose 95% CrI is [-0.169, 0.151], i.e. not distinguishable from zero.
**Fix:** `cgr_summary_table()` suppresses `pct_attenuation` (NA) when the
unadjusted estimate's 95% CrI includes zero. abs_attenuation is retained.

### CH-15  Student-t robustness sensitivity now runs
The Rmd mentioned `cgr_jags(likelihood = "t")` but never executed it (zero
occurrences of `jags-t` in any render). With JAGS working there was no blocker,
and the skew table (PANAS/ACAC skew -0.94, excess kurtosis 1.42) makes it worth
running. **Added:** a Section 11 chunk that fits the normal and Student-t
likelihoods and reports the estimated nu. Result: nu ~ 18, so the t collapses
toward the normal - the Gaussian conclusion is robust. `cgr_jags()` now exposes
the posterior-mean nu via `attr(out, "nu")`.

### CH-17  cgr_operating() crashed on empty strata  [BREAKING - was fatal]
`cgr_operating(n = 120, p_cg = 0.85, n_trials = 500)` - the exact "is CGR safe
for my trial?" call advertised in the README and Section 15 - could abort with
"empty stratum: PLAC; the estimand is undefined". At a high correct guess rate
with small n, a wrong-guess stratum occasionally comes up empty in a simulated
trial; `cgr_strata()` then throws and killed the whole run.
**Fix:** the simulation loop skips degenerate trials instead of crashing and
reports two new columns - `empty_stratum_rate` (fraction of trials with an empty
stratum) and `n_valid` (trials actually used). The rate is 0 at n=230/p_cg=0.7
but can exceed 0.5 for small, badly-unblinded designs, which is itself the
"CGR adjustment is fragile here" signal the tool exists to give. Regression test
added in tests/testthat/test-reproduction.R.

## Findings from the author's source code (2026-07-21)

### CH-16  Source code located: szb37/CorrectGuessRateCurve
Named in the 2023 paper's data-availability statement; read at last (see
SOURCES.md S4). Consequences:
- **0.72 is hardcoded** (`trial_cgrs = {'sbmd': 0.72}`), not computed - settles
  U3. The code separately computes the data CGR as 0.647 and does not use it for
  the reference line.
- **Estimand confirmed independently**: the code's r and s formulas are
  identical to this implementation - settles U4.
- **"100 times" qualified**: the repo's Figure-4 config `cgrC_low` specifies
  `n_cgrc_trials = 32` over `np.linspace(0, 1, 13)` (options 32/64/96, never
  100). Stated as the repo's configuration, not proof of the final figure - the
  repo tracks the preprint and the Figure-4 block is behind `if False:`. If the
  count is 32 the Monte Carlo error is ~1.8x the 100-resample value. Section 2
  and the KDE ladder updated (32 added).
- **legacy_round is the faithful reproduction**: `get_strata_ratio` does
  `round(x, 2)`, so `legacy_round = TRUE` reproduces Szigeti's numbers
  (+0.010 PANAS, +0.019 Energy). Documented in Section 8; the exact ratios
  remain the default for the estimand/Bayesian sections.
- **U9 raised**: guess rates are reversed vs the 2024 review's "higher in the
  active arms" (placebo 0.723 vs microdose 0.528 here).

## Extensions (additive; the binary method is unchanged)

### CH-17  UNKNOWN-preserving CGRC extension  [ADDITIVE - no existing result changes]
A six-stratum extension (ACAC/ACPL/ACU/PLAC/PLPL/PLU, R/10_unknown.R) that keeps
an observed "I do not know" guess as a third response category instead of
dropping it, counting it as wrong, counting it as placebo, or splitting it
across arms. It holds the observed UNKNOWN-response rate `u` fixed and varies the
DIRECTIONAL correct-guess rate `c` (correct among AC/PL responders), reweighting
class mass while preserving the observed within-class arm shares r, s and a third
share t = ACU/(ACU+PLU).
**Verified properties (tests/testthat/test-unknown.R):** six weights sum to 1
and split mass by class; `Delta(c_obs, u_obs)` equals the raw arm-mean difference
exactly; at `u = 0` every formula reduces exactly to the four-stratum
`cgr_delta()`; an empty stratum is always forced to structurally zero weight
through r/s/t and is never estimated; exact Santana-Penin (u = 26/77, c = 29/51)
and ketamine (u = 11/38, c = 14/27) count tables reproduce.
**Backend agreement:** `cgr_unknown_check_backends()` ran (JAGS 4.3.2) - PASS
(max |z| = 1.13, max Rhat = 1.0, min ESS ~31700). Student-t path runs.
**Explicitly NOT claimed:** that this proves independence of assignment and all
three guess categories; that `c = 0.50` is "perfect blinding"; that the
extension outperforms anything (no operating-characteristic simulation with
UNKNOWN exists - see UNRESOLVED U10). It is an extension by cgrc.bayes, not the
Szigeti estimand.

### CH-18  Guess normalisation, input audit, and Markdown report
`cgrc_normalise_guess()` recognises explicit UNKNOWN synonyms and keeps blank/NA
as MISSING (missing data and an observed UNKNOWN are different). `cgrc_input_audit()`
classifies every uploaded row (missing condition/guess/outcome, non-numeric
outcome, observed UNKNOWN) instead of a silent `complete.cases()`.
`cgrc_build_report()` emits a self-contained Markdown report.

### CH-19  Shiny Panel B: UNKNOWN mode, unit threshold, seed, audit downloads
Panel B detects UNKNOWN responses and offers preserve (default) vs binary
complete-case (which reports the exact number excluded); shows six-cell strata
and directional-CGR labels; lets the meaningful-difference threshold be a
fraction of SD OR outcome units (no 0.5 cap); takes an explicit seed; offers
cleaned-data, exclusion-log and report downloads; and disables the binary design
lookup bridge for an UNKNOWN-preserving analysis (the lookup does not model
UNKNOWN).

### CH-20  Optional sensitivities (all default-off)
`cgr_unknown_conjugate(ratio_uncertainty=TRUE)` propagates r/s/t sampling error;
`cgr_unknown_independent()` is a separate EXPERIMENTAL shared-guess-distribution
estimand; `cgr_unknown_jags(pooling="partial")` is an assumption-dependent
hierarchical partial-pooling sensitivity. The independent-stratum, conditional
model remains the default everywhere.

### CH-21  Terminology: "perfect blinding" qualified
Unqualified "perfect blinding" replaced with "target CGR 0.50 / guessing at
chance" in the CGR-curve subtitle, the summary-table label, and the front-door
docs, and reworded in the binary headline. CGR 0.50 is a target correct-guess
rate, not proof that assignment and guessing are independent.

### CH-23  UNKNOWN-aware generative model and operating characteristics  [addresses U10]
`sim_aeb_unknown()` extends the AEB model with an observed "I do not know"
response under two explicit, deliberately-chosen assumptions: (A1) the
UNKNOWN-response rate is equal in both arms, and (A2) an UNKNOWN responder carries
no expectancy. `cgr_unknown_operating()` is the six-stratum analogue of
`cgr_operating()` and reports adjusted bias, RMSE, 95% coverage, favourable-flag
rates, the naive-t significance rate, and the (higher, six-stratum)
empty-stratum rate; `cgr_unknown_min_stratum()` is the matching feasibility
early-warning. **Result under A1/A2:** the estimator of `Delta(0.50, u_obs)` is
essentially unbiased for the direct effect with ~0.95 coverage, and controls the
pure-expectancy false-favourable rate (~0.05 adjusted vs ~0.73 for the naive
t-test). This addresses the empirical half of U10 (the interpretive assumptions
remain). The Shiny Panel B now runs this on demand at the uploaded trial's n,
directional CGR and observed UNKNOWN rate, replacing the previously-disabled
design bridge; the binary Panel A lookup is untouched and still not reused for
UNKNOWN designs. `cgrc_normalise_guess()` also gained apostrophe-free and a few
extra UNKNOWN synonyms ("dont know", "no idea", "cant tell", ...).

### CH-22  Seed arguments on the posterior front doors
`cgrc()`, `cgrc_headline()`, `cgrc_unknown()`, `cgrc_unknown_headline()` and
`cgr_unknown_independent()` gained an optional `seed`, recorded in the returned
object. Default `NULL` preserves prior behaviour (an external `set.seed`).

## Preserved unchanged

- The estimand: strata, r, s, weights, Delta(c).
- The observed-CGR identity check (exact-ratio path).
- Default vague priors (m0=0, k0=1e-6, a0=b0=1e-3).
- Week-1-only filter and use of raw value as outcome.
- The supplementary stratum-allocation validation.
- The published empirical results: PANAS 3.16 -> 1.08 etc. are unchanged.
- The binary API: cgrc(), cgr_strata(), cgr_delta(), cgr_conjugate(), cgr_jags()
  and every published four-stratum result are byte-for-byte unchanged; the
  extension lives in new functions.
