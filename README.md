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
| `tests/testthat/` | estimand + reproduction test suites |
| `data-raw/download_data.R` | fetch + checksum-verify the public data |
| `data/pacutes.csv` | pinned public dataset (SHA-256 in `PROVENANCE.txt`) |
| `00_BRIEF.md` | the design brief, corrections list, and empirical findings |
| `reports/` | changelog, unresolved items, test status, audit, plain-language, sources |

## Reproduce it

```sh
# system: R (>= 4.3), pandoc, and JAGS (>= 4.3) for the backend check
Rscript data-raw/download_data.R          # fetch + verify data/pacutes.csv
Rscript -e 'rmarkdown::render("CGRC_bayes.Rmd")'
Rscript tests/testthat.R                  # 297 pass / 0 fail / 0 skip
```

R package dependencies: `ggplot2`, `digest`, `testthat`, `knitr`, `rmarkdown`,
and — for the backend check — `rjags` + `coda` (needs a system JAGS install).
Without JAGS the render still completes; the conjugate-vs-JAGS check is skipped
and reported as not run.

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
