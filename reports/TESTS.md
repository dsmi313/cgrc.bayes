# Test status

Honest three-way split. "VERIFIED" means it was actually computed. "NOT RUN"
means no R or no JAGS was available and nobody should assume it passes.

**Update (2026-07-24, this environment):** the suite has grown well beyond the
July-21 state as the UNKNOWN-preserving extension and the Shiny readability pass
landed. On R 4.3.3 the full testthat suite is green — **1711 pass / 0 fail**,
with **3 skips** that are JAGS-only (the `cgr_unknown_jags` checks; JAGS was not
installed in this run). The Student-t robustness sensitivity (N4) DOES run at
render time via the `sens-student-t` chunk (`eval = has_jags`) and passes when
JAGS is present (nu ~ 18; see below and CHANGELOG CH-15) — it is not outstanding.

**Earlier update (2026-07-21):** R 4.3.3 + pandoc 3.1.3 + JAGS 4.3.2 were
installed and everything previously blocked was run; the Rmd rendered end to end
with zero chunk errors, the then-297-test suite was green (0 fail / 0 skip), and
the load-bearing backend check PASSED.

## VERIFIED - computed against an independent reimplementation

| # | Test | Result |
|---|---|---|
| T1 | Supplementary stratum allocation -> 54/58/42/46 | PASS |
| T2 | Weights sum to 1; correct-guess weights sum to c; incorrect to 1-c | PASS |
| T3 | Endpoint behaviour at c=0 and c=1 reduces to single-stratum contrasts | PASS |
| T4 | Observed-CGR identity, 36 simulated configurations | PASS, worst err 7.1e-15 |
| T5 | Observed-CGR identity on all four real outcomes | PASS, err 0.00e+00 |
| T6 | Szigeti Table 1 recovery, n=20000, all four cells | PASS |
| T7 | Szigeti Table 2 reproduction, real data | PASS (unadj <=0.12, adj <=0.30) |
| T8 | KDE ladder converges to analytic at 10 000 resamples | PASS, all 8 cells |
| T9 | KDE averaged p-values reproduce published p-values | PASS |
| T10 | Hedges g under noise="all" equals published 0.40 | PASS (0.4011) |
| T11 | 500-trial coverage of 95% CrI | PASS (0.936-0.960) |
| T12 | 500-trial bias of adjusted estimator, all 4 scenarios | PASS (<=0.05) |
| T13 | Observed CGR is NOT 0.72 on the public week-1 data | PASS (0.6466) |

## NOW RUN - executed 2026-07-21 with R 4.3.3 + JAGS 4.3.2

| # | Test | Result |
|---|---|---|
| N1 | The Rmd renders end to end | PASS, zero chunk errors -> CGRC_bayes.html |
| N2 | `cgr_check_backends()` conjugate vs corrected JAGS | PASS, verdict "identity claim supported" |
| N3 | JAGS convergence (Rhat, ESS) | PASS, max Rhat 1.0001, min ESS 39664 |
| N5 | testthat suite executes | PASS, 1711 pass / 0 fail (2026-07-24; 3 JAGS-only skips without JAGS. Was 297 pass / 0 fail / 0 skip on 2026-07-21 with JAGS) |
| N6 | Prior sensitivity table renders | PASS |
| N7 | Small-strata degradation table renders | PASS |

**N2 was load-bearing and now RESOLVES to PASS.** Full-grid (101-point)
comparison of the conjugate backend (40 000 draws) against the corrected 4-chain
JAGS model (10 000 iterations/chain): posterior means agree to <= 0.006 across
the grid, both 95% credible limits to <= 0.026 (lower) / 0.079 (upper),
posterior probabilities to <= 0.005, and the largest mean difference is 0.63
Monte Carlo standard errors. Section 8's claim that the two backends target an
identical posterior is therefore supported by this render. Reproducing it
requires a JAGS install; a render without rjags falls back to the NOT RUN note
and the claim reverts to unverified.

## NOW RUN (continued)

| # | Test | Result |
|---|---|---|
| N4 | Student-t robustness sensitivity | PASS - `cgr_jags(likelihood = "t")` runs; estimated nu ~ 18 (PANAS), so the t collapses toward the normal and the Gaussian conclusion is robust despite the ACAC skew (-0.94) |

## STILL NOT RUN

*(none — every previously-blocked item has now been executed)*

## KNOWN FAILURE, now fixed

| # | Failure | Fix |
|---|---|---|
| X1 | Render halted: "object 'method' not found" in chunk `fig3-table` | `stats::filter` was masking `dplyr::filter`; `stats::filter` has a `method` argument, so `method == "conjugate"` was evaluated eagerly as its `filter` argument. `R/` is now base R only. |

## Automated tests added

`tests/testthat/`: weights, denominators, endpoints, invalid codes, empty
strata, tiny strata, seed reproducibility, identity, allocation, published
targets.
