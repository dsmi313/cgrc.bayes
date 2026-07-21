# CGRC refactor — brief, corrections, and findings

Read this before touching anything else.

> **Execution status (2026-07-21).** The "first four actions" in Part 5 have now
> been run in an environment with R 4.3.3, pandoc 3.1.3, and JAGS 4.3.2:
> 1. Data fetched and checksum-verified (SHA-256 matches).
> 2. `CGRC_bayes.Rmd` renders end to end with **zero chunk errors**.
> 3. `cgr_check_backends()` ran with JAGS and **PASSES** — the conjugate and
>    corrected-JAGS posteriors agree to <= 0.006 in the mean across the grid
>    (max |z| = 0.63, max R-hat = 1.0001, min ESS = 39664), so DECISION-1's
>    identity claim (U1) is now **supported**, not merely asserted.
> 4. `tests/testthat.R` is green: **297 pass / 0 fail / 0 skip**.
>
> The only formerly-blocked item still outstanding is the Student-t robustness
> sensitivity (TESTS.md N4), which the Rmd does not exercise at render time.
> The DECISION flags below remain live for the human to review; nothing about
> the three decisions was changed by this execution. See `reports/TESTS.md`,
> `reports/UNRESOLVED.md` (U1), and `reports/CHANGELOG.md` (CH-01) for details.

## Status of the approval gate

The source prompt asked for an outline and a list of statistical corrections
**before** rewriting. The instruction accompanying it asked for the finished
work. Both are here: this file is the outline and correction list, and the
rewritten `CGRC_bayes.Rmd` follows it. Nothing below was silently applied — every
substantive change is in `reports/CHANGELOG.md` with a rationale, and every
change can be reverted independently.

Three decisions were made without approval because the work could not proceed
otherwise. They are the first three items in the correction list and are flagged
`DECISION` so they can be overruled cheaply.

---

## Part 1 — Empirical findings established during this refactor

These were computed, not asserted. Each is reproducible from the scripts in
`R/` and `data-raw/`.

### F1. The AEB `noise` ambiguity is RESOLVED in favour of `"all"`

This overturns a claim made in earlier drafts of this document, which said
neither variant reproduced both published quantities and left the choice open.
That claim was wrong, and it was wrong because it compared incomparable things
(a posterior probability against a frequentist p-value).

Two independent lines of evidence, both from 500–2000 simulated trials at
n = 230, p_CG = 0.7:

**Evidence 1 — published effect size.** Szigeti states the DTE corresponds to
Hedges' g = 0.40.

| variant | Hedges' g | error |
|---|---|---|
| `noise = "all"` | **0.4011** | 0.0011 |
| `noise = "arm"` | 0.5023 | 0.1023 |

**Evidence 2 — published unadjusted significance rates (Table 1).**

| DTE | AEB | published | `"all"` | `"arm"` |
|---|---|---|---|---|
| off | off | 0.05 | **0.056** | 0.056 |
| on | off | 0.86 | **0.876** | 0.970 |
| off | on | 0.78 | **0.826** | 0.930 |
| on | on | 0.99 | **0.992** | 0.998 |

`"all"` matches on every row; `"arm"` is badly off on rows 2 and 3.

**Caveat, stated plainly.** Szigeti's source code was NOT located. The public
repo `szb37/mcrds_public` serves `data/pacutes.csv` but returns 404 for
`README.md` and for every source path probed. This is strong convergent
empirical evidence, not a reading of the author's code. Figure 3 should be
described as a **reproduction consistent with the published operating
characteristics**, not a verified line-by-line reproduction. If the author can
supply the generative code, this should be re-checked.

**Action:** `noise = "all"` is now the documented default with the evidence
recorded inline. The `"arm"` option is retained so the comparison can be re-run.

### F2. The 0.72 vs 0.647 CGR discrepancy has a strong candidate explanation

The paper and associated materials quote a trial CGR of 0.72. The public week-1
data gives 0.647. Pooling all timepoints does **not** explain it:

| timepoint | n | CGR |
|---|---|---|
| w1s1 | 232 | 0.6466 |
| w2s1 | 218 | 0.6789 |
| w3s1 | 205 | 0.6195 |
| w4s1 | 192 | 0.6510 |
| all pooled | 847 | 0.6494 |
| every row, every scale | 15623 | 0.6370 |
| participant-level mean | 239 | 0.6468 |

But splitting week 1 by **actual allocation**:

| actual arm | n | correct | rate |
|---|---|---|---|
| microdose | 91 | 48 | 0.5275 |
| **placebo** | **141** | **102** | **0.7234** |

**0.7234 ≈ 0.72.** The most likely explanation is that the quoted 0.72 is the
correct-guess rate *conditional on having received placebo*, not the overall
correct-guess rate. This is a hypothesis with strong numerical support, not a
confirmed fact — it should be put to the author directly.

It matters substantively, not just cosmetically. The arms are badly unbalanced
(141 placebo vs 91 microdose in week 1) and guessing is heavily skewed toward
"placebo" (9508 placebo guesses vs 6115 microdose guesses across the dataset).
Participants who received the microdose identified it barely above chance
(0.528). A single scalar "CGR" hides that asymmetry, and the CGRC estimand
conditions on it through `s`.

**Action:** the reference line is computed, never hardcoded, and the audit table
reports the per-arm rates alongside the overall rate.

### F3. n = 232 vs the published n = 233 is unexplained

No scale at any timepoint has n = 233. Eight scales at w1s1 have exactly 232
(PANAS, positive, negative, energy, mood, creativity, focus, temper);
`intensity` has 231; cognitive scales have 186–191. There are no duplicate
`trial_id` values and no blank or `NA` outcome values at w1s1.

The discrepancy is one record and cannot be reproduced from the public file. It
is most likely a record withheld from the public release or a typo. **Not
resolved.** Flagged for the author.

### F4. KDE resampling is an inefficient Monte Carlo estimator of the analytic value

A faithful port of the original procedure (sklearn `KernelDensity` with default
parameters, i.e. Gaussian kernel and bandwidth 1.0; resample to target CGR; fit
`outcome ~ treatment`; average estimates and p-values) gives:

| scale | target | 100 reps | 1 000 reps | 10 000 reps | analytic |
|---|---|---|---|---|---|
| PANAS | obs | 2.856 | 3.227 | 3.153 | **3.157** |
| PANAS | 0.50 | 0.895 | 1.120 | 1.090 | **1.080** |
| Mood VAS | obs | 5.889 | 6.315 | 6.329 | **6.337** |
| Mood VAS | 0.50 | 2.432 | 2.487 | 2.553 | **2.517** |
| Energy VAS | obs | 11.678 | 11.270 | 11.348 | **11.377** |
| Energy VAS | 0.50 | 7.433 | 6.940 | 7.103 | **7.104** |
| CPS | obs | 0.014 | 0.001 | −0.009 | **−0.011** |
| CPS | 0.50 | −0.000 | 0.003 | 0.007 | **0.006** |

Two conclusions, and they are the core of the contribution:

1. **KDE vs a Gaussian stratum model contributes essentially nothing to the
   point estimate.** At 10 000 resamples the KDE average converges on the
   analytic value in every cell. The KDE is unbiased for the stratum mean, and
   the estimand only uses means.
2. **100 resamples is materially noisy.** Monte Carlo SE at 100 reps is
   0.12–0.29 points on the VAS/PANAS scales. PANAS at CGR 0.50 moves from 0.895
   to 1.090 — roughly 20% — purely from resampling noise.

The averaged p-values also reproduce the published ones closely (PANAS 0.409 vs
published 0.43; Mood 0.384 vs 0.42; Energy 0.043 vs 0.04; CPS 0.504 vs 0.52),
which is what confirms the port is faithful rather than merely plausible.

**Framing consequence.** The Bayesian version should NOT be sold as better
inference about the trial. Its defensible claims are narrower and stronger: it
removes a Monte Carlo error term that has no scientific content, and it replaces
an average of p-values — a quantity with no accepted interpretation — with a
posterior distribution that has one.

### F5. Operating characteristics (500 trials per scenario, n = 230, `noise="all"`)

True Δ(0.5) equals the direct treatment effect, because AEB cancels at perfect
blinding.

| DTE | AEB | true | unadj bias | **adj bias** | RMSE | 95% coverage | P(post > 0.95) | freq p < .05 |
|---|---|---|---|---|---|---|---|---|
| off | off | 0 | 0.02 | **0.05** | 0.59 | 0.938 | 0.066 | 0.056 |
| on | off | 3 | 0.06 | **0.05** | 1.09 | 0.936 | 0.892 | 0.876 |
| off | on | 0 | **3.13** | **0.02** | 1.09 | 0.940 | 0.056 | 0.826 |
| on | on | 3 | **3.12** | **0.05** | 1.33 | 0.960 | 0.704 | 0.992 |

Findings:

- **The adjusted estimator is unbiased in all four scenarios** (|bias| ≤ 0.05),
  including both scenarios where AEB inflates the unadjusted estimate by ~3.1
  points.
- **Credible interval coverage is nominal** (0.936–0.960 against 0.95). This is
  the first evidence in this project that the intervals mean what they say.
- **The false-favourable rate is controlled** — 0.056 and 0.066 against a 0.05
  target — while the unadjusted frequentist test fires at 0.826 in the pure-AEB
  scenario.
- **Adjustment costs power**: 0.892 → 0.704 in the partial-mediation row. This
  is a real cost and must be reported, not buried.

These are operating characteristics of the *adjusted estimand*, computed under
the AEB data-generating model. They are not evidence that the Bayesian
implementation beats the KDE implementation; F4 already shows both target the
same quantity. A head-to-head comparison is listed as unresolved.

---

## Part 2 — Statistical corrections

### DECISION-1 — JAGS prior corrected to match the conjugate model

Was (independent priors, so **not** the same posterior as the conjugate model):

```
mu[j]  ~ dnorm(0, 1.0E-6)
tau[j] ~ dgamma(1.0E-3, 1.0E-3)
```

Now (conditional precision, matching Normal-Inverse-Gamma exactly):

```
mu[j]  ~ dnorm(m0, k0 * tau[j])
tau[j] ~ dgamma(a0, b0)
```

with `m0`, `k0`, `a0`, `b0` passed through the data list. Under NIG,
Var(μ | σ²) = σ²/k₀, so the precision of μ is k₀·τ. The old code fixed the prior
precision at 1e-6 regardless of τ, which is a different model.

**This has not been executed.** No JAGS in the authoring environment. The claim
that the two backends target an identical posterior is now
*conditional on a verification that has not been run*. The Rmd states this and
`reports/TESTS.md` lists it as NOT RUN. If verification fails, delete the
identity claim rather than weakening the prior to match.

### DECISION-2 — the p-value-shaped curve is removed

The two-sided posterior tail probability with a line at 0.05 is gone. It
invited reading a posterior as a p-value. Replaced with
**P(Δ(c) favourable | y)**, labelled in full, with an optional descriptive line
at 0.95 explicitly marked as not a universal cutoff. A `direction` argument
handles outcomes where lower is better (QIDS, STAIT), so favourability is
declared rather than assumed.

### DECISION-3 — `S` renamed `n_draws` throughout

`S` was ambiguous against participant sample size. The Rmd now states in the
text that posterior draws reduce Monte Carlo error in posterior summaries only,
and carry no information about participants.

### C4 — reproduction and extension separated

Sections 2–4 reproduce the published estimand and curves. Sections 5–7 develop
the Bayesian implementation. No comparative claim appears before the simulation
study in Section 9. The word "better" does not appear in the document.

### C5 — worked example added

`Δ(c_obs)` and `Δ(0.5)` for PANAS, every intermediate number printed, hand
checkable. See Part 3 below.

### C6 — `nig_draws()` can return paired (μ, σ²)

`nig_draws(..., return_sigma2 = TRUE)` returns both so the conditional
dependence within a draw can be inspected. The Rmd states explicitly that draws
are independent *across* iterations but μ and σ² are dependent *within* a draw.

### C7 — original KDE procedure implemented for comparison

`R/cgrc_kde.R`, run at 100/1000/10000 resamples. See F4.

### C8 — repeated-simulation validation replaces single-trial illustration

`R/cgrc_sim.R`, 500 trials per scenario. See F5.

### C9 — data pinned locally with provenance

`data-raw/download_data.R` records URL, access date, and SHA-256. The Rmd reads
the local copy and verifies the hash. No network access at render time.

- URL: `https://raw.githubusercontent.com/szb37/mcrds_public/master/data/pacutes.csv`
- Accessed: 2026-07-21
- SHA-256: `86aa784528ee045c61fadf3eacfd3e1897d16aae9839cee7cb4bfe839a7cc4e3`
- Size: 977 449 bytes, 15 623 records

### C10 — optional uncertainty in r and s

`Beta(n_PLPL + α, n_ACAC + α)` and `Beta(n_ACPL + α, n_PLAC + α)`, **off by
default**, because the original estimand conditions on the observed ratios.
Turning it on is an extension, not a reproduction, and the Rmd says so.

### C11 — reusable code moved to `R/`

The Rmd sources `R/` and no longer holds the implementation.

---

## Part 3 — The worked example (hand-checkable)

PANAS, week 1 acute (`tp == "w1s1"`), n = 232.

**Stratum counts and means**

| stratum | n | mean | sd | sum |
|---|---|---|---|---|
| ACAC | 48 | 19.5625 | 8.4272 | 939 |
| ACPL | 43 | 12.3256 | 10.7276 | 530 |
| PLAC | 39 | 18.3590 | 7.2855 | 716 |
| PLPL | 102 | 10.9314 | 9.7962 | 1115 |

**Observed CGR** = (48 + 102) / 232 = 150/232 = **0.646552**

**Ratios**
- r = 102 / (102 + 48) = **0.680000**
- s = 43 / (43 + 39) = **0.524390**

**At the observed CGR (c = 0.646552)**

| weight | formula | value |
|---|---|---|
| w_ACAC | c(1−r) | 0.206897 |
| w_ACPL | (1−c)s | 0.185345 |
| w_PLAC | (1−c)(1−s) | 0.168103 |
| w_PLPL | cr | 0.439655 |

Checks: 0.206897 + 0.439655 = 0.646552 = c ✓  0.185345 + 0.168103 = 0.353448 = 1−c ✓

- active arm = 6.331897 / 0.392241 = 16.142857
- placebo arm = 7.892241 / 0.607759 = 12.985816
- **Δ(0.646552) = 3.157042**

**At CGR = 0.50**

| weight | value |
|---|---|
| w_ACAC | 0.160000 |
| w_ACPL | 0.262195 |
| w_PLAC | 0.237805 |
| w_PLPL | 0.340000 |

- active arm = 6.361707 / 0.422195 = 15.068169
- placebo arm = 8.082520 / 0.577805 = 13.988321
- **Δ(0.50) = 1.079847**

**Identity check.** Raw active-minus-placebo sample mean difference =
3.1570415400. Δ(c_obs) = 3.1570415400. Absolute error **0.00e+00** (exact).

**NIG posterior, stratum ACAC**, with m₀ = 0, k₀ = 1e−6, a₀ = b₀ = 1e−3:

n = 48, ȳ = 19.562500, SS = 3337.8125
- kₙ = 1e−6 + 48 = 48.000001
- mₙ = 48 × 19.5625 / 48.000001 = 19.562500 (differs from ȳ by −4.08e−07)
- aₙ = 0.001 + 24 = 24.0010
- bₙ = 0.001 + 1668.90625 + 3.4e−05 = 1668.9074
- E[σ²] = bₙ/(aₙ−1) = 72.5580, against sample variance 71.0173

The gap between mₙ and ȳ is 4e−07, which is what "vague prior" means
operationally: the prior contributes k₀ = 1e−6 pseudo-observations against 48
real ones.

---

## Part 4 — Proposed document structure

1. Research question and executive summary
2. The original CGRC procedure, described accurately
3. The four-stratum estimand
4. Worked calculation by hand
5. Bayesian probability model
6. Posterior derivation
7. R implementation
8. Validation against the original method
9. Simulation study
10. Empirical microdosing results
11. Sensitivity analyses
12. Assumptions and limitations
13. What the Bayesian approach adds
14. Unresolved questions
15. Session information and reproducibility

---

## Part 5 — What Claude Code should do first

1. Render `CGRC_bayes.Rmd`. It has **never been rendered**; see
   `reports/TESTS.md` for what is verified vs assumed.
2. Install JAGS and run the DECISION-1 verification. Until it passes, the
   identity claim in Section 8 is unsupported. If it fails, delete the claim.
3. Run `testthat::test_dir("tests/testthat")`.
4. Do not resolve anything in `reports/UNRESOLVED.md` by choosing the
   convenient answer. Two items need the original author.
