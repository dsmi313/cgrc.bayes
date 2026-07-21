# Statistical audit: every assumption

Grouped by where it enters. "Testable here" means testable with these data.

## A. Assumptions of the CGRC estimand itself (inherited, not introduced)

| # | Assumption | Consequence if false | Testable here |
|---|---|---|---|
| A1 | Unblinding is malicious: PT -> TE -> OUT. Guessing causes expectancy, which causes outcome. | If benign (OUT -> PT), adjustment removes a REAL effect and manufactures a false negative. The paper warns about this explicitly. | No |
| A2 | The four observed strata are exchangeable with the corresponding strata in a hypothetical perfectly blinded trial. | Reweighting transports the wrong population. Someone who guesses right under poor blinding may differ from someone who guesses right by luck under good blinding. | No |
| A3 | Within-class arm ratios r and s are the correct thing to hold fixed. | Reweighting distorts arm balance while changing guess rate. | Partly (Sec 11) |
| A4 | r and s are known constants, not estimates. | Understated uncertainty. Quantified in Sec 11: Beta-uncertain widens the interval. | Yes (Sec 11) |
| A5 | Only stratum MEANS matter; the estimand is linear in them. | Any tail- or quantile-based question needs a different estimator. | By construction |
| A6 | Treatment belief is adequately captured by one binary bit. | A confident correct guess and a lucky one are treated identically. Confidence data would separate them. | No |
| A7 | Expectancy acts additively with the direct treatment effect. | The AEB model assumes linear addition; the paper flags that effects may not be additive. | No |
| A8 | CGR = 0.50 is the right counterfactual target. | Perfect blinding is defined as chance-level guessing. Defensible, but a definition rather than a fact. | No |

## B. Assumptions introduced by the Bayesian implementation

| # | Assumption | Consequence if false | Testable here |
|---|---|---|---|
| B1 | Outcomes within a stratum are normal. | Posterior means are robust; intervals less so. Skew/kurtosis reported in Sec 11. | Partly |
| B2 | Observations are independent within stratum. | Week-1-only filter is what buys this. Pooling weeks would break it via repeated measures. | By design |
| B3 | Strata have independent parameters (no pooling). | No borrowing of strength. Conservative for small strata; a hierarchical model would shrink. | Yes, not run |
| B4 | Priors are vague enough not to matter. | Sec 11 shows the vague and weakly-informative priors agree; the strongly informative one shifts the estimate, as it should. | Yes (Sec 11) |
| B5 | Each stratum has its own variance. | No homoscedasticity imposed. Strictly weaker than the pooled-variance t-test in the original. | By construction |
| B6 | n_draws is large enough that Monte Carlo error in summaries is negligible. | MCSE reported alongside every estimate. | Yes |
| B7 | Conjugate and corrected-JAGS posteriors coincide. | **UNVERIFIED.** See UNRESOLVED U1. | Yes, NOT RUN |

## C. Assumptions about the data

| # | Assumption | Status |
|---|---|---|
| C1 | Week 1 only, so each datapoint is independent. | Follows the paper. |
| C2 | Raw value is the outcome; no baseline adjustment. | Forced: no baseline exists in the public data. Reproduces the published curves. |
| C3 | condition MD/PL maps to active/placebo. | Verified in the audit table. |
| C4 | guess MD/PL maps to guessed-active/guessed-placebo. | Verified. |
| C5 | The public file matches the analysed file. | **QUESTIONABLE**: n = 232 vs published 233. See UNRESOLVED U2. |
| C6 | The observed CGR is 0.647, not 0.72. | Verified from data; 0.72 is a hardcoded constant in the author's code (`trial_cgrs`), matching the placebo-arm rate. See U3, SOURCES.md S4. |
| C7 | No missing outcome values at w1s1. | Verified: zero blanks on all four scales. |
| C8 | No duplicate participants at w1s1. | Verified: 232 rows, 232 unique trial_id. |

## D. Assumptions in the simulation study

| # | Assumption | Status |
|---|---|---|
| D1 | The AEB model is the right data-generating process. | It is the paper's model. Operating characteristics are conditional on it and would differ under a misspecified DGP. |
| D2 | `noise = "all"` is the intended Eq 4 reading. | Confirmed from source: default `strata_sampling = 'all_prop'` in CorrectGuessRateCurve. See U4, SOURCES.md S4. |
| D3 | True Delta(0.5) equals the direct treatment effect. | Follows from the model: at CGR 0.5 the guess distribution is identical across arms, so AEB cancels. |
| D4 | 500 trials is enough. | Monte Carlo SE on coverage is about sqrt(.05*.95/500) = 0.010, so 0.936 vs 0.95 is within noise. |

## E. Claims deliberately NOT made

- That the Bayesian method is better. Sec 8 shows both converge to the same estimand.
- That the energy VAS result is a direct drug effect.
- That attenuation measures the amount of effect caused by expectancy.
- That reproducing the published curve validates any causal assumption.
- That the energy VAS result survives every robustness axis untested (the
  Student-t check now runs; nu ~ 18, so it is robust for PANAS).
- That the AEB simulation was reproduced line-by-line from the author's
  simulation module (the estimand's r/s formulas are confirmed from source, but
  the Eq. 4 SD scope was not re-derived line-by-line).
