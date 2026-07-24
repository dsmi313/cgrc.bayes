# CGRC — a Bayesian implementation of the Correct Guess Rate Curve

A reproduction of the Correct Guess Rate Curve (CGRC) method of
**Szigeti et al. (2023)** and an alternative implementation that expresses the
uncertainty around the CGR-adjusted estimate as a **posterior distribution**
rather than through resampling. The estimand is unchanged throughout — this is a
different way of computing uncertainty around the published quantity, not a
different quantity.

The CGRC asks a counterfactual of an imperfectly blinded trial: *what would the
treatment effect have been if the correct guess rate had been 50% — i.e. if
blinding had held?* It reweights the four treatment × guess strata to a target
guess rate while holding the within-class arm ratios fixed, so only *how much
guessing happened* changes, never *who got what*.

## Use it on your own trial

Szigeti et al. (2023) ask, in their limitations, that "researchers wishing to use
CGR adjustment should first run simulations to determine whether CGR produces
acceptable error rates for the parameters of their data." This package is that
tool.

**Try it right now** — no data file needed; `sim_aeb()` makes a fake unblinded
trial you can adjust immediately after `R CMD INSTALL .`:

```r
library(cgrc.bayes)
demo <- sim_aeb(n = 230, p_cg = 0.7, dte_on = TRUE)   # simulated trial
cgrc(demo)                                            # adjusted estimate + curve
plot(cgrc(demo))                                      # the effect + P(favourable) figure
```

**On your own data.** `cgrc()` needs one row per participant and exactly three
columns — `condition` (AC/PL), `guess` (AC/PL) and `value`. Build that data
frame however you like, e.g. from a CSV:

```r
d <- read.csv("my_trial.csv")
my_trial <- data.frame(
  condition = ifelse(d$arm   == "drug",  "AC", "PL"),   # active / placebo
  guess     = ifelse(d$guess == "drug",  "AC", "PL"),   # guessed active / placebo
  value     = d$outcome
)
cgrc(my_trial)                        # for lower-is-better outcomes: cgrc(my_trial, direction = -1)
```

**Before trusting it** — the pre-flight simulation the paper asks for, at *your*
sample size, guess rate and effect size:

```r
cgr_operating(n = 120, p_cg = 0.85, n_trials = 500)   # bias, RMSE, 95% coverage, error rates
```

**Or use the interactive app** — "Is CGR adjustment safe for my trial?":

```r
install.packages("shiny")
cgrc_app()
```

Panel A (Design) reads a precomputed simulation grid instantly: a plain-language
verdict, a power-vs-n curve, the false-positive/power trade-off with and without
adjustment, a feasibility readout (smallest expected stratum, % of simulated
trials with an empty stratum), and a button to run the exact simulation at your
settings. Panel B (Analyse) takes a CSV of your own trial, adjusts it with
`cgrc()`, and shows the CGR curve, the ROPE decomposition, and the observed-CGR
identity check — then offers to run the design check at your trial's own n and
observed guess rate.

`cgrc()` returns the CGR-adjusted estimate (at a target correct-guess rate of
0.50 — guessing at chance, not a claim of proven perfect blinding), a
95% credible interval, and the posterior probability the effect is favourable.
`cgr_operating()` tells you whether that adjustment is trustworthy for a trial of
that size and blinding quality in the first place — including an
`empty_stratum_rate` column: at a high guess rate with small `n`, a wrong-guess
stratum can come up empty and the estimand is undefined, so a high rate there is
a warning that CGR adjustment is fragile for your design.

## Function reference

| Function | Purpose |
|---|---|
| `cgrc_app()` | launch the "is CGR safe for my trial?" Shiny app |
| `cgrc(df)` | one-call adjuster: curve + estimate at CGR 0.50 + P(favourable) |
| `cgr_operating(n, p_cg, ...)` | pre-flight simulation: bias, RMSE, 95% coverage, error rates (+ `empty_stratum_rate`) |
| `cgr_min_stratum(n, p_cg)` | expected smallest stratum — the feasibility early warning |
| `cgr_conjugate(df, grid)` | full Normal-Inverse-Gamma posterior curve over a CGR grid |
| `cgr_jags(df, likelihood)` | JAGS backend; `"normal"` or robust `"t"` likelihood |
| `cgr_check_backends(df)` | verify the conjugate and JAGS posteriors agree |
| `cgr_rope(df)` | region-of-practical-equivalence decomposition (harm/negligible/benefit) |
| `cgr_kde(df, cgr)`, `cgr_kde_curve(df)` | faithful port of the original KDE resampling |
| `cgr_reference_line_test(df, orig_cgr)` | check a reference line against the observed-CGR identity |
| `szigeti_panel(cgr_kde_curve(df), ...)` | reproduce the published twin-axis figure |
| `cgr_strata`, `cgr_ratios`, `cgr_weights`, `cgr_delta` | the estimand primitives |
| `sim_aeb(...)` | the activated-expectancy-bias generative model |
| `cgrc_unknown(df)` | UNKNOWN-preserving extension: six-stratum adjuster keeping "I do not know" |
| `cgrc_unknown_headline(df)` | two-probability plain-language summary of the UNKNOWN extension |
| `cgr_unknown_jags(df, pooling)` | six-stratum JAGS backend (`"normal"`/`"t"`; `"none"`/`"partial"` pooling) |
| `cgr_unknown_check_backends(df)` | verify the UNKNOWN conjugate and JAGS posteriors agree |
| `cgr_unknown_independent(df)` | experimental shared-guess-distribution estimand (distinct from CGRC) |
| `cgrc_normalise_guess(x)` | map a guess column to AC/PL/UNKNOWN; blank/NA stay missing |
| `cgrc_input_audit(...)`, `cgrc_build_report(f)` | per-row input audit and a Markdown analysis report |

## Extension for UNKNOWN treatment guesses

Real trials often let a participant answer **"I do not know"** to the guess
question. Those responses must not be silently dropped, counted as incorrect
guesses, counted as placebo, or split across arms without an explicit assumption.
`cgrc.bayes` adds an **UNKNOWN-preserving extension** that keeps them as a third
response category.

> **This is an extension implemented by cgrc.bayes. It is not part of the
> original CGRC formulation of Szigeti et al.** The binary four-stratum method
> and every published result above are unchanged.

**Six strata** instead of four — `ACAC, ACPL, ACU, PLAC, PLPL, PLU` (arm ×
guess, where `U` = UNKNOWN). Let

```
u_obs = n_unknown / n_total                          (observed UNKNOWN-response rate)
c_obs = n_correct / n_directional                    (correct among AC/PL responders)
r = PLPL/(PLPL+ACAC)   s = ACPL/(ACPL+PLAC)   t = ACU/(ACU+PLU)
```

be the observed UNKNOWN rate, the **directional** correct-guess rate, and the
three preserved within-class arm shares. For a target directional CGR `c` and
target UNKNOWN rate `u` the six weights are

```
w_ACAC=(1-u)c(1-r)   w_PLPL=(1-u)c·r        (correct class, mass (1-u)c)
w_ACPL=(1-u)(1-c)s   w_PLAC=(1-u)(1-c)(1-s) (incorrect class, mass (1-u)(1-c))
w_ACU =u·t           w_PLU =u(1-t)          (UNKNOWN class, mass u)
```

and the contrast is `Δ(c,u) = μ_AC(c,u) − μ_PL(c,u)` with each arm mean the
weight-averaged stratum mean within that arm. Key properties:

- **Observed-value identity.** `Δ(c_obs, u_obs)` equals the raw active−placebo
  mean difference exactly (each within-arm weight reduces to `n_stratum/n_total`).
- **Reduction.** At `u = 0` every formula collapses **exactly** to the binary
  `cgr_delta()` — the extension is a strict generalisation.
- **Empty cells are safe.** An empty stratum forces its own weight to zero
  through `r/s/t` (e.g. empty `ACU` ⇒ `t = 0` ⇒ `w_ACU = 0`), so a
  structurally-zero cell is never estimated. An empty *correct* or *incorrect*
  directional class is a clear "undefined estimand" error, not a silent guess.

**Why UNKNOWN ≠ an incorrect guess.** A participant who says "I don't know" has
made no directional claim; scoring them as wrong (or as placebo) invents
information the trial did not collect. Holding `u` fixed and reweighting only the
directional guesses keeps the UNKNOWN mass where it was observed.

**Interpretation and limits.** `c = 0.50` means *directional guessing at chance
with the UNKNOWN rate held fixed* — **not** "perfect blinding", and not proof
that assignment and all three guess categories are independent. The preserved
share `t` is an assumption (the observed within-UNKNOWN arm ratio), exactly as
`r, s` are. There is **no** operating-characteristic simulation with UNKNOWN yet,
so the extension is not calibration-validated and must not be claimed to
outperform the binary method (see `reports/UNRESOLVED.md` U10).

```r
u_trial <- data.frame(condition = ..., guess = ..., value = ...)   # guess may be UNKNOWN
cgrc_unknown(u_trial)                 # six-stratum adjusted analysis
cgrc_unknown_headline(u_trial)        # the two-probability plain-language summary
cgr_unknown_check_backends(u_trial)   # conjugate vs JAGS agreement (needs JAGS)
```

Optional, default-off sensitivities: `cgr_unknown_conjugate(ratio_uncertainty=TRUE)`
(propagate r/s/t uncertainty), `cgr_unknown_jags(pooling="partial")` (hierarchical
partial pooling), and `cgr_unknown_independent()` (a *separate*, experimental
shared-guess-distribution estimand). In the Shiny app, Panel B detects UNKNOWN
responses and offers "preserve" (default) or "binary complete-case (exclude)"
with a full input audit and downloadable cleaned data, exclusion log, and report.

## What's here

| Path | Contents |
|---|---|
| `CGRC_bayes.Rmd` | the writeup: reproduction, worked example, Bayesian model, validation, simulation study, empirical results, sensitivity, limitations |
| `R/01_estimand.R` | strata, ratios, weights, Δ(c), analytic curve (base R only) |
| `R/02_bayes.R` | `nig_draws()` — Normal-Inverse-Gamma posterior; conjugate CGR curve |
| `R/03_jags.R` | corrected JAGS backend and `cgr_check_backends()` |
| `R/04_kde.R` | faithful port of the original KDE resampling procedure |
| `R/05_sim.R` | AEB generative model and operating-characteristics study |
| `R/06_plot.R` | figures and summary tables |
| `R/07_rope.R` | `cgr_rope()` region-of-practical-equivalence decomposition |
| `R/08_frontdoor.R` | `cgrc()` one-call adjuster and `cgrc_headline()` two-probability plain-language summary |
| `R/09_app.R` | Shiny app support: lookup accessor, interpolation, verdict, input normalisation, audit, report |
| `R/10_unknown.R` | UNKNOWN-preserving six-stratum extension (estimand, posterior, ROPE, headline, sensitivities) |
| `inst/app/app.R` | the Shiny app (Design + Analyse panels) |
| `data-raw/build_lookup.R` | precompute the operating-characteristics grid (run once, ~1h) |
| `inst/extdata/cgrc_lookup.rds` | the precomputed grid, shipped as package data |
| `DESCRIPTION`, `NAMESPACE` | R package metadata (`cgrc.bayes`) |
| `tests/testthat/` | estimand + reproduction + rope + app test suites |
| `data-raw/download_data.R` | fetch + checksum-verify the public data |
| `data/pacutes.csv` | pinned public dataset (SHA-256 in `PROVENANCE.txt`) |
| `00_BRIEF.md` | the design brief, corrections list, and empirical findings |
| `reports/` | changelog, unresolved items, test status, audit, plain-language, sources |

## Install and reproduce

The code is an installable R package, `cgrc.bayes`; the writeup is a vignette
that consumes it via `library(cgrc.bayes)`.

```sh
# system: R (>= 4.3), pandoc, and JAGS (>= 4.3) for the backend check
R CMD INSTALL .                           # install the cgrc.bayes package
Rscript data-raw/download_data.R          # fetch + verify data/pacutes.csv
Rscript -e 'rmarkdown::render("CGRC_bayes.Rmd")'
Rscript tests/testthat.R                  # runs the testthat suite (source, no install needed)
```

R package dependencies: `ggplot2`, `digest`, `testthat`, `knitr`, `rmarkdown`,
and — for the backend check — `rjags` + `coda` (needs a system JAGS install).
Without JAGS the render still completes; the conjugate-vs-JAGS check is skipped
and reported as not run. The test suite sources `R/` directly, so it does not
require the package to be installed first.

## Status of this reproduction

- The estimand reproduces the published Table 2 (PANAS 3.2 → 1.1, Energy VAS
  11.5 → 6.8, only Energy surviving adjustment). The observed-CGR identity is
  exact.
- The conjugate posterior and the corrected 4-chain JAGS model target the **same
  posterior** (agreement to ≤ 0.006 in the mean across the grid; verified with
  JAGS 4.3.2).
- Over 500 simulated trials per scenario the adjusted estimator is unbiased
  (|bias| ≤ 0.05) with nominal 95% credible-interval coverage (0.936–0.964).

The Bayesian version is **not** claimed to be better inference about the trial:
it targets the same estimand, removes a Monte Carlo error term with no
scientific content, and replaces an average of p-values with a posterior
probability. Open questions and the source-paper checks behind them are in
`reports/UNRESOLVED.md`; every substantive change is logged in
`reports/CHANGELOG.md`.

## Sources

Full details and per-finding citations are in `reports/SOURCES.md`.

- Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D. "The difference between
  'placebo group' and 'placebo control': a case study in psychedelic
  microdosing." *Scientific Reports* 13:12107 (2023).
  <https://doi.org/10.1038/s41598-023-34938-7> — **the method paper reproduced
  here.**
- Szigeti B, Weiss B, Rosas FE, Erritzoe D, Nutt D, Carhart-Harris R. "Assessing
  expectancy and suggestibility in a trial of escitalopram v. psilocybin for
  depression." *Psychological Medicine* 54, 1717–1724 (2024).
  <https://doi.org/10.1017/S0033291723003653>
- Szigeti B, Heifets BD. "Expectancy Effects in Psychedelic Trials."
  *Biological Psychiatry: Cognitive Neuroscience and Neuroimaging* 9:512–521
  (2024). <https://doi.org/10.1016/j.bpsc.2024.02.004>
- Public dataset: Szigeti B, et al. "Self-blinding citizen science to explore
  psychedelic microdosing." *eLife* 10:e62878 (2021).
  <https://doi.org/10.7554/eLife.62878>

This project reproduces and re-analyses published work; it is independent of and
not endorsed by the original authors.
