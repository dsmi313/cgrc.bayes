# Putting the Correct Guess Rate Curve into practice

### A Bayesian formulation with explicit uncertainty, extended to ‘I don’t know’ responses

*2026-08-01*
# Summary

The Correct Guess Rate Curve (Szigeti et al. 2023) is reimplemented here with the estimand unchanged, inside a Bayesian framework that reports posterior uncertainty directly instead of averaging resampled p-values. It is then extended to the answer the four-cell design cannot hold, “I don’t know,” under a six-cell scheme that reduces exactly to the original whenever no such answers occur. Both run behind a point-and-click front end, so a trial researcher can apply the adjustment, and check whether it is trustworthy for their design, without writing code.

**The reproduction holds.** The observed-CGR identity holds to machine precision for every outcome tested. The CGR-adjusted point estimates match the published Table 2 to within about 0.3 points across all four scales. The unadjusted significance rates reproduce across all four AEB scenarios, and a 20,000-trial paired simulation shows the original and Bayesian procedures recovering the same point estimate, correlation 0.995 or better in every scenario. The published *CGR-adjusted* rejection rates are not reproduced, and source rounding does not explain the gap; the most likely explanation is a mis-specification in my implementation of the AEB generative model, and I would welcome a pointer to where.

**Uncertainty is reported directly.** Stratum means get conjugate Normal-Inverse-Gamma posteriors and are pushed through the same reweighting, giving a credible interval and posterior probabilities for the adjusted contrast in place of an average of resampled p-values. No MCMC.

**“I don’t know” is preserved rather than discarded.** Six-cell reweighting keeps the third response as its own category instead of recoding it as a wrong guess or dropping the participant. It satisfies the same observed-value identity, reduces exactly to the original estimator when no such responses occur, and comes with its own operating characteristics. It is also sensitive to unequal UNKNOWN rates across arms, so the package reports that split as a diagnostic.

**A feasibility warning.** The adjustment needs participants in all four cells, and the wrong-guess cells fill only from the minority who guess incorrectly. At a high guess rate with a small sample they empty out and the recalculation is undefined rather than imprecise. `cgr_operating()` runs that check for an arbitrary design, following the paper’s own recommendation that researchers simulate before applying.

Two things the numbers do not support. Shrinkage under adjustment is not a measure of how much of an effect was expectancy: participants guess after treatment, and what they guess is partly caused by the drug, so reweighting shuffles who counts how much rather than tracing where the effect came from. Loewinger et al. (2026) give the formal version of this. The CGRC is treated throughout as a sensitivity analysis for unblinding and expectancy, not a causal decomposition.

What the shrinkage does tell you is how much weight a result can bear. The PANAS effect does not survive the assumption that guessing drove it; the Energy VAS effect does. That is a statement about fragility, and it is the one these numbers support.

# 1 The estimand

Crossing assigned arm with guessed arm gives four strata, written arm first: `ACAC`, `ACPL`, `PLAC`, `PLPL`. Write $\bar y_k$ for a stratum’s mean outcome and $\rho_k$ for its share of participants, so the observed correct-guess rate is $c_{obs}=\rho_{ACAC}+\rho_{PLPL}$. Two ratios record how correct and incorrect guessers split between arms: 

$$r=\frac{\rho_{PLPL}}{\rho_{PLPL}+\rho_{ACAC}},\qquad s=\frac{\rho_{ACPL}}{\rho_{ACPL}+\rho_{PLAC}}.$$

 To impose a target guess rate $c$, reweight the four strata to $w_{ACAC}=c(1-r)$, $w_{PLPL}=cr$, $w_{ACPL}=(1-c)s$, $w_{PLAC}=(1-c)(1-s)$, preserving those arm splits within the correct and incorrect classes and renormalizing each arm: 

$$\Delta(c)=\frac{w_{ACAC}\bar y_{ACAC}+w_{ACPL}\bar y_{ACPL}}{w_{ACAC}+w_{ACPL}} -\frac{w_{PLPL}\bar y_{PLPL}+w_{PLAC}\bar y_{PLAC}}{w_{PLPL}+w_{PLAC}}.$$

 The primary target is $\Delta(0.5)$.

Because $r$ and $s$ are read off the data, setting $c=c_{obs}$ leaves the weighting unchanged and returns the plain arm contrast, $\Delta(c_{obs})=\bar y_{AC}-\bar y_{PL}$. This **observed-CGR identity** is a necessary condition only: any estimator deriving its weights from the data’s own $r$ and $s$ satisfies it automatically, including one built on a mistaken estimand. It catches implementation errors, not conceptual ones, and is used that way throughout.

```r
st  <- cgr_strata(panas); rat <- cgr_ratios(st); o <- cgr_observed(st)
idn <- cgr_reference_line_test(panas, orig_cgr = o)
knitr::kable(
  data.frame(stratum = STRATA, n = lengths(st)[STRATA],
             mean = round(vapply(st[STRATA], mean, numeric(1)), 2)),
  caption = sprintf("PANAS week 1: four strata; observed CGR = %.4f, r = %.3f, s = %.3f",
                    o, rat$r, rat$s))
```

|      | stratum |   n |  mean |
|:-----|:--------|----:|------:|
| ACAC | ACAC    |  48 | 19.56 |
| ACPL | ACPL    |  43 | 12.33 |
| PLAC | PLAC    |  39 | 18.36 |
| PLPL | PLPL    | 102 | 10.93 |

PANAS week 1: four strata; observed CGR = 0.6466, r = 0.680, s = 0.524

The identity holds to machine precision: $\Delta(c_{obs})=$ 3.157 against a raw difference of 3.157, discrepancy 0.0e+00. The adjusted estimate is $\Delta(0.5)=$ 1.08, so the PANAS contrast falls from roughly 3.2 points to about 1 once correct guessing is imagined away. (Reproducing the shipped code exactly requires a rounded-ratio variant, retained as `legacy_round = TRUE`; the estimand uses the exact ratios above.)

# 2 A Bayesian formulation

The original summarizes inference by averaging p-values across resampled pseudo-datasets. An average of p-values is not a conventional inferential quantity, so this implementation propagates uncertainty through the estimand directly instead.

Each stratum mean and variance gets a weak Normal-Inverse-Gamma prior, which is conjugate, so the posterior is closed form and can be sampled exactly with no MCMC. The draws exist because $\Delta(c)$ is a nonlinear function of four stratum means: each mean has a tractable posterior, their ratio of weighted sums does not. Drawing the stratum means jointly and pushing each complete draw through the same $\Delta(c)$ yields a posterior sample of the adjusted contrast, from which the credible interval and probabilities below are read. The draws are a transformation device, not an inference engine, and being independent their Monte Carlo error falls as $1/\sqrt{n}$.

Two probabilities carry the practical message and answer different questions: $P(\text{favorable})$, the posterior probability the adjusted contrast is on the beneficial side of zero, and $P(\text{meaningful})$, the probability it exceeds a pre-declared smallest difference worth caring about. Where a trial declares a minimal important difference that value is used; otherwise the fallback is half the outcome’s standard deviation, following Norman et al. (2003), whose review across health-related quality-of-life measures found minimal important differences clustering near that value. Szigeti et al. (2024) adopt the same 0.5 SMD bound for equivalence testing, and Szigeti and Heifets (2024) use it as the benchmark against which microdose effect sizes are judged. It remains a cross-condition generalization rather than an outcome-specific threshold, so a result reported against it should be read as “beyond half an SD” rather than as validated clinical meaningfulness.

```r
fit_panas <- cgrc(panas, n_draws = N_DRAWS, seed = 2)
# cgrc_headline() supplies the meaningfulness probabilities against the declared
# delta; its two values correspond to the two rows of fit_panas$summary.
h_panas <- cgrc_headline(panas, n_draws = N_DRAWS, seed = 2)
tab_panas <- fit_panas$summary[, c("what", "post_mean", "cri_lo", "cri_hi", "p_favourable")]
tab_panas$p_meaningful <- c(h_panas$p_meaningful_obs, h_panas$p_meaningful_blind)
knitr::kable(tab_panas, digits = 3, row.names = FALSE,
             caption = sprintf(paste0("PANAS: CGR-adjusted posterior (Normal-Inverse-Gamma ",
                                      "conjugate). p_meaningful is measured against ",
                                      "delta = %.2f points, half the outcome SD."),
                               h_panas$delta))
```

| what                              | post_mean | cri_lo | cri_hi | p_favourable | p_meaningful |
|:----------------------------------|----------:|-------:|-------:|-------------:|-------------:|
| observed (unadjusted)             |     3.151 |  0.659 |  5.700 |        0.993 |        0.074 |
| reweighted to CGR 0.50 (adjusted) |     1.074 | -1.594 |  3.785 |        0.784 |        0.002 |

PANAS: CGR-adjusted posterior (Normal-Inverse-Gamma conjugate). p_meaningful is measured against delta = 5.01 points, half the outcome SD.

# 3 Checking it against the original

**Reproducing the published behavior.** Szigeti et al. characterized the method by simulating trials under their activated expectancy bias (AEB) model. Re-running that simulation at the paper’s settings recovers the same pattern of significance rates across all four scenarios, though not the same numbers. The largest absolute discrepancy is 0.020. Absolute size is the wrong yardstick, though, because a rate near 0.99 is estimated far more precisely from 500 trials than a rate near 0.8. Scaled by the Monte Carlo standard error of each published rate, no scenario differs by more than 1.3 standard errors. That comparison is conservative, since the published rates carry the same sampling error from their own 500 trials and the original used a different seed and a different resampling procedure. Agreement at this level is what a faithful reimplementation should produce; it is not an exact match and should not be read as one. The same 500-trial run is reused as the main stress test below, so nothing is cherry-picked.

```r
knitr::kable(data.frame(scenario = scen,
  `unadjusted p<0.05 (2-sided)` = round(op$freq_sig, 3),
  `published Table 1` = pub_t1, check.names = FALSE),
  caption = "AEB simulation, n = 230, CGR = 0.7, 500 trials, against Szigeti et al. (2023) Table 1")
```

| scenario                  | unadjusted p\<0.05 (2-sided) | published Table 1 |
|:--------------------------|-----------------------------:|------------------:|
| no effect / no expectancy |                        0.050 |              0.05 |
| effect / no expectancy    |                        0.870 |              0.86 |
| no effect / expectancy    |                        0.800 |              0.78 |
| effect / expectancy       |                        0.996 |              0.99 |

AEB simulation, n = 230, CGR = 0.7, 500 trials, against Szigeti et al. (2023) Table 1

**The identity holds for every outcome.** At the actual guess rate the adjusted curve returns the plain unadjusted effect, up to floating-point rounding.

```r
id_tab <- do.call(rbind, lapply(names(sb), function(nm) {
  z <- cgr_reference_line_test(sb[[nm]], orig_cgr = cgr_observed(cgr_strata(sb[[nm]])))
  data.frame(outcome = nm, observed_CGR = round(z$computed_obs_cgr, 4),
             `Delta(c_obs)` = round(z$D_at_obs, 3), `raw AC-PL` = round(z$raw_mean_diff, 3),
             `abs error` = sprintf("%.1e", abs(z$D_at_obs - z$raw_mean_diff)),
             check.names = FALSE)
}))
id_max <- max(vapply(sb, function(d) {
  z <- cgr_reference_line_test(d, orig_cgr = cgr_observed(cgr_strata(d)))
  abs(z$D_at_obs - z$raw_mean_diff)
}, numeric(1)))
knitr::kable(id_tab, row.names = FALSE,
             caption = "Observed-CGR identity across the four outcomes")
```

| outcome               | observed_CGR | Delta(c_obs) | raw AC-PL | abs error |
|:----------------------|-------------:|-------------:|----------:|:----------|
| PANAS                 |       0.6466 |        3.157 |     3.157 | 0.0e+00   |
| Mood VAS              |       0.6466 |        6.337 |     6.337 | 1.4e-14   |
| Energy VAS            |       0.6466 |       11.377 |    11.377 | 0.0e+00   |
| Cognitive performance |       0.6290 |       -0.011 |    -0.011 | 0.0e+00   |

Observed-CGR identity across the four outcomes

The largest residual across the four outcomes is 1.4e-14, which is the scale of double-precision arithmetic on these numbers rather than a substantive difference. This is a check on the implementation, not on the estimand.

One further check lives in the package rather than here: an exact, unrounded KDE ladder confirming that the original’s resampling converges on the analytic value. A 20,000-trial paired comparison of the two inference procedures is reported below, under operating characteristics.

# 4 Worked example: the microdosing trial

Following the paper, only each participant’s first-week single-timepoint record is used. One point I could not resolve: the public data give an overall guess rate of 0.647 whereas the published figure annotates 0.72, a value hard-coded in the plotting script. I could not determine from the public materials which sample or timepoint 0.72 refers to; the identity check places the unadjusted effect at the data-derived value, so that is what is used here.

```r
emp <- do.call(rbind, lapply(names(sb), function(nm) {
  d <- sb[[nm]]; f <- cgrc(d, n_draws = N_DRAWS, seed = 2)
  h <- cgrc_headline(d, n_draws = N_DRAWS, seed = 2); s <- f$summary
  data.frame(outcome = nm, n = nrow(d), observed_CGR = round(f$observed_cgr, 3),
             `raw` = round(s$post_mean[1], 2), `adjusted (0.50)` = round(s$post_mean[2], 2),
             `95% CrI` = sprintf("[%.2f, %.2f]", s$cri_lo[2], s$cri_hi[2]),
             `P(fav)` = sprintf("%.0f%%", 100 * s$p_favourable[2]),
             `P(meaningful)` = sprintf("%.0f%%", 100 * h$p_meaningful_blind),
             check.names = FALSE)
}))
# Published CGR-adjusted point estimates, Szigeti et al. (2023) Table 2, in the
# order PANAS / Mood VAS / Energy VAS / CPS.
emp$`published (Table 2)` <- c(1.1, 2.7, 6.8, 0.0)
knitr::kable(emp, caption = "CGR-adjusted results, four outcomes (adjusted = target CGR 0.50; meaningful = beyond half an outcome SD, per Norman et al. 2003)")
```

| outcome               |   n | observed_CGR |   raw | adjusted (0.50) | 95% CrI         | P(fav) | P(meaningful) | published (Table 2) |
|:----------------------|----:|-------------:|------:|----------------:|:----------------|:-------|:--------------|--------------------:|
| PANAS                 | 232 |        0.647 |  3.15 |            1.07 | \[-1.59, 3.79\] | 78%    | 0%            |                 1.1 |
| Mood VAS              | 232 |        0.647 |  6.33 |            2.51 | \[-3.15, 8.28\] | 80%    | 0%            |                 2.7 |
| Energy VAS            | 232 |        0.647 | 11.36 |            7.09 | \[2.18, 12.10\] | 100%   | 10%           |                 6.8 |
| Cognitive performance | 186 |        0.629 | -0.01 |            0.01 | \[-0.16, 0.17\] | 53%    | 0%            |                 0.0 |

CGR-adjusted results, four outcomes (adjusted = target CGR 0.50; meaningful = beyond half an outcome SD, per Norman et al. 2003)

```r
plot(fit_panas, title = "PANAS, week 1")
```

![](https://raw.githubusercontent.com/dsmi313/cgrc.bayes/main/report-files/summary/fig01.png)

The final column is the published CGR-adjusted estimate for the same outcome. The posterior means agree with it to within about 0.3 points across all four scales, which is the closest thing here to a direct reproduction check on empirical data rather than on simulation.

For PANAS the raw effect of about 3.2 falls to 1.1 at the target rate, its adjusted interval (\[-1.6, 3.8\]) straddles zero, and $P(\text{favorable})$ falls to about 78%. Mood VAS behaves the same way. Energy VAS is the survivor, its adjusted interval sitting entirely above zero. Cognitive performance is near zero throughout, so the adjustment changes little, a reassuring null. The adjustment is therefore not a blunt instrument that shrinks everything: it separates results fragile to the blinding assumption from those that are not.

# 5 Operating characteristics

Szigeti et al. recommended that researchers simulate how the adjustment behaves for their own design before applying it. Reusing the same 500-trials-per-scenario run, across no effect, a treatment component, an expectancy component, and both. In each row, “flag” means the fraction of simulated trials in which the posterior probability of a favorable effect exceeded 0.95, the Bayesian equivalent of a positive test result; “adjusted flag” applies that bar to the CGR-adjusted contrast, “unadjusted flag” to the raw one:

```r
n_eff <- op$n_valid    # trials where the adjustment is defined (returned by cgr_operating)
knitr::kable(data.frame(
  scenario = scen,
  `true Delta` = round(op$true, 1),
  `adj bias` = round(op$adj_bias, 3),
  RMSE = round(op$adj_rmse, 2),
  `95% coverage` = round(op$coverage95, 3),
  `coverage MCSE` = round(mcse_prop(op$coverage95, n_eff), 3),
  `flagged, adjusted (P(fav.) > 0.95)` = round(op$p_fav_gt_95, 3),
  `flagged, unadjusted (P(fav.) > 0.95)` = round(op$unadj_p_fav_gt_95, 3),
  check.names = FALSE),
  caption = "Operating characteristics, n = 230, CGR = 0.7, 500 trials/scenario")
```

| scenario                  | true Delta | adj bias | RMSE | 95% coverage | coverage MCSE | flagged, adjusted (P(fav.) \> 0.95) | flagged, unadjusted (P(fav.) \> 0.95) |
|:--------------------------|-----------:|---------:|-----:|-------------:|--------------:|------------------------------------:|--------------------------------------:|
| no effect / no expectancy |          0 |    0.025 | 0.59 |        0.944 |         0.010 |                               0.068 |                                 0.056 |
| effect / no expectancy    |          3 |    0.015 | 1.06 |        0.964 |         0.008 |                               0.874 |                                 0.942 |
| no effect / expectancy    |          0 |   -0.031 | 1.11 |        0.944 |         0.010 |                               0.048 |                                 0.868 |
| effect / expectancy       |          3 |    0.003 | 1.42 |        0.936 |         0.011 |                               0.694 |                                 0.998 |

Operating characteristics, n = 230, CGR = 0.7, 500 trials/scenario

The adjusted estimate is essentially unbiased throughout, and its interval covers at a rate consistent with nominal. Monte Carlo standard error on a coverage estimate near 0.95 with 500 replicates is about 0.010, so differences of a point or two in that column are noise rather than scenario-specific behavior.

The headline result is row 3: where the effect is pure expectancy, the raw analysis crosses the P(favorable) \> 0.95 bar almost every time (0.868) while the adjusted analysis stays quiet (0.048). The price appears in rows 2 and 4, where the adjusted analysis crosses that bar less often than the unadjusted one. This holds even in row 2, where no expectancy effect exists at all: re-estimating under a counterfactual guess distribution relies on the treatment by guess strata rather than the simple randomized comparison, so the estimate is more variable. The correction is not free.

The same machinery answers a feasibility question. The adjustment needs participants in each of the four cells, and the two wrong-guess cells are filled from the $(1-\text{CGR})$ minority. At a high guess rate with a small sample they empty out and the recalculation is undefined. A deliberately punishing design shows the failure:

```r
sp  <- cgr_operating(n = 50, p_cg = 0.95, mu_dte = 3, mu_aeb = 7.7, n_trials = 300, seed = 1)
r10 <- sp[sp$DTE == 1 & sp$AEB == 0, ]
knitr::kable(data.frame(design = "n = 50, CGR = 0.95",
  `empty-stratum rate` = sprintf("%.0f%%", 100 * r10$empty_stratum_rate),
  `adjusted bias (feasible trials)` = round(r10$adj_bias, 2),
  `95% coverage` = round(r10$coverage95, 3),
  `coverage MCSE` = round(mcse_prop(r10$coverage95, r10$n_valid), 3),
  check.names = FALSE),
  caption = "Sparse-stratum diagnostic: a thin discordant stratum leaves many trials undefined")
```

| design             | empty-stratum rate | adjusted bias (feasible trials) | 95% coverage | coverage MCSE |
|:-------------------|:-------------------|--------------------------------:|-------------:|--------------:|
| n = 50, CGR = 0.95 | 52%                |                           -0.01 |        0.793 |         0.034 |

Sparse-stratum diagnostic: a thin discordant stratum leaves many trials undefined

About half those trials have no defined answer, and those that do have degraded coverage. `cgr_operating(n, p_cg, n_trials)` runs this for an arbitrary design, which is the check a researcher should run before committing to a CGRC adjustment.

## 5.1 The two inference procedures compared directly

A separate 20,000-trial experiment holds the estimand and the data fixed and varies only the inferential procedure. Each replicate generated one dataset with `sim_aeb()` at `n = 230`, `p_cg = 0.70`, `mu_dte = 3`, `mu_aeb = 7.7`, and analyzed it both ways: the original-style KDE resampling (exact, unrounded) and the conjugate Bayesian posterior, paired at the trial level so differences cannot come from the two methods seeing different data.

```r
pair_sum    <- read_ext("paired_cgrc_simulation_summary.csv")
pair_trials <- read_ext("paired_cgrc_simulation_trials.csv")
scen_order  <- c("no effect / no expectancy", "effect / no expectancy",
                 "no effect / expectancy",    "effect / expectancy")
pair_sum    <- pair_sum[match(scen_order, pair_sum$scenario), ]

knitr::kable(data.frame(
  scenario = scen_order,
  `true Delta` = pair_sum$true_delta,
  `KDE bias` = round(pair_sum$kde_bias, 3),
  `Bayes bias` = round(pair_sum$bayes_bias, 3),
  `KDE RMSE` = round(pair_sum$kde_rmse, 3),
  `Bayes RMSE` = round(pair_sum$bayes_rmse, 3),
  `estimate correlation` = round(pair_sum$estimate_correlation, 3),
  check.names = FALSE), row.names = FALSE,
  caption = "Paired point-estimate comparison, 5,000 trials per scenario")
```

| scenario                  | true Delta | KDE bias | Bayes bias | KDE RMSE | Bayes RMSE | estimate correlation |
|:--------------------------|-----------:|---------:|-----------:|---------:|-----------:|---------------------:|
| no effect / no expectancy |          0 |    0.007 |      0.008 |    0.577 |      0.575 |                0.995 |
| effect / no expectancy    |          3 |   -0.010 |     -0.011 |    1.076 |      1.068 |                0.996 |
| no effect / expectancy    |          0 |    0.007 |      0.010 |    1.072 |      1.068 |                0.996 |
| effect / expectancy       |          3 |   -0.012 |     -0.009 |    1.406 |      1.401 |                0.996 |

Paired point-estimate comparison, 5,000 trials per scenario

**The two procedures compute the same quantity.** Estimate correlation is 0.995 or above in every scenario, biases agree to within 0.003 outcome units, and RMSEs to within 0.008. Whatever separates them, it is not the target.

The decision rules differ. Averaging p-values across resampled pseudo-datasets pulls the summary toward the middle of its distribution, diluting extremes at both ends, so the KDE rule rejects in the null scenarios at roughly a fifth of its nominal level and is correspondingly less sensitive when an effect is present. Bayesian interval coverage is close to nominal throughout.

One discrepancy is unresolved. The KDE detection rates obtained here, roughly 0.61 with a treatment effect alone and 0.33 with treatment and expectancy together, are well below the published CGR-adjusted 0.84 and 0.82. A separate sensitivity run confirmed source-code rounding does not account for it: `source_faithful = TRUE` and `FALSE` on 1,000 identical DTE + AEB datasets gave the same flag proportion to three decimal places. The most likely place to look is my implementation of the AEB generative model, in particular my reading of the Eq. 4 variance terms, which I could not pin down from the published description. Clarification would settle it quickly. Until then, this experiment should be read as a comparison of two inferential procedures under this package’s implementation of the AEB model, not as a characterization of the published procedure. Full decision-rate tables and the paired 2x2 breakdown are in the repository.

# 6 Handling “I don’t know” answers

Real blinding questionnaires often offer a third response. The four-cell method has no place for it, forcing a choice between discarding randomized data and inventing a belief that was never reported. **The six-stratum extension below is my proposal, not part of the original method.**

Crossing assignment with three responses gives $ACAC$, $ACPL$, $ACU$, $PLAC$, $PLPL$, $PLU$. Let $u_{obs}=\rho_{ACU}+\rho_{PLU}$ be the observed UNKNOWN rate, and let the directional correct-guess rate be computed among directional guessers only, $c_{obs}=(\rho_{ACAC}+\rho_{PLPL})/(1-u_{obs})$. An UNKNOWN response is not an incorrect guess. Three ratios preserve the arm split within the correct, incorrect and UNKNOWN classes: 

$$r=\frac{\rho_{PLPL}}{\rho_{PLPL}+\rho_{ACAC}},\qquad s=\frac{\rho_{ACPL}}{\rho_{ACPL}+\rho_{PLAC}},\qquad t=\frac{\rho_{ACU}}{\rho_{ACU}+\rho_{PLU}}.$$

 Targets $c$ and $u$ then give 

$$\begin{aligned} w_{ACAC}&amp;=(1-u)c(1-r), &amp; w_{PLPL}&amp;=(1-u)cr,\\ w_{ACPL}&amp;=(1-u)(1-c)s, &amp; w_{PLAC}&amp;=(1-u)(1-c)(1-s),\\ w_{ACU}&amp;=ut, &amp; w_{PLU}&amp;=u(1-t), \end{aligned}$$

 and $\Delta(c,u)$ is the same active-minus-placebo comparison of weighted means over the three cells per arm. For the primary adjustment I hold $u=u_{obs}$ and set $c=0.5$: keep the observed proportion who said “I don’t know,” and reweight those who did guess so their directional accuracy is chance.

Two algebraic checks follow. At $c_{obs}$ and $u_{obs}$ each weight reduces to its observed stratum proportion, so $\Delta(c_{obs},u_{obs})=\bar y_{AC}-\bar y_{PL}$. And $u=0$ zeroes the UNKNOWN weights and returns the remaining four to the binary expressions, so $\Delta(c,0)=\Delta_{\text{binary}}(c)$: the extension is strict. Both are asserted in the package’s tests and checked numerically below against weights written out directly from the algebra rather than taken from the package.

```r
# Delta(c, u) from the six target weights, implemented from the algebra above.
delta_cu <- function(ybar, r, s, t, c, u) {
  w <- c(ACAC = (1 - u) * c * (1 - r), ACPL = (1 - u) * (1 - c) * s,
         ACU  = u * t,
         PLPL = (1 - u) * c * r,       PLAC = (1 - u) * (1 - c) * (1 - s),
         PLU  = u * (1 - t))
  arm <- function(k) { k <- k[w[k] > 0]; sum(w[k] * ybar[k]) / sum(w[k]) }
  arm(c("ACAC", "ACPL", "ACU")) - arm(c("PLPL", "PLAC", "PLU"))
}

# Six-stratum counts, means and ratios; UNKNOWN cells padded for a binary trial.
six_summary <- function(d, has_unknown) {
  st <- if (has_unknown) cgr_unknown_strata(d) else cgr_strata(d)
  nm <- if (has_unknown) UNKNOWN_STRATA else STRATA
  n  <- lengths(st)[nm]; names(n) <- nm
  m  <- vapply(nm, function(k) if (n[[k]] > 0) mean(st[[k]]) else 0, numeric(1))
  if (!has_unknown) { n <- c(n, ACU = 0L, PLU = 0L); m <- c(m, ACU = 0, PLU = 0) }
  rho <- n / sum(n)
  u   <- rho[["ACU"]] + rho[["PLU"]]
  list(ybar = m, u = u,
       c = (rho[["ACAC"]] + rho[["PLPL"]]) / (1 - u),
       r = rho[["PLPL"]] / (rho[["PLPL"]] + rho[["ACAC"]]),
       s = rho[["ACPL"]] / (rho[["ACPL"]] + rho[["PLAC"]]),
       t = if (u > 0) rho[["ACU"]] / (rho[["ACU"]] + rho[["PLU"]]) else 0)
}

raw_diff <- function(d) mean(d$value[d$condition == "AC"]) -
                        mean(d$value[d$condition == "PL"])

set.seed(20260724)
sim_u <- sim_aeb_unknown(n = 240, p_cg = 0.70, u = 0.25, dte_on = TRUE)

rows <- list()

# Identity 1: at the observed c and u, the curve returns the raw arm contrast.
for (nm in c("simulated 6-stratum", "PANAS (no UNKNOWN)")) {
  d  <- if (nm == "simulated 6-stratum") sim_u else panas
  hu <- nm == "simulated 6-stratum"
  g  <- six_summary(d, hu)
  lhs <- delta_cu(g$ybar, g$r, g$s, g$t, g$c, g$u)
  rhs <- raw_diff(d)
  rows[[length(rows) + 1]] <- data.frame(
    check = "Delta(c_obs, u_obs) = raw AC - PL", dataset = nm,
    computed = round(lhs, 6), reference = round(rhs, 6),
    `abs error` = sprintf("%.1e", abs(lhs - rhs)), check.names = FALSE)
}

# Identity 2: at u = 0 the six-stratum contrast equals the binary estimand.
gb <- six_summary(panas, FALSE)
bst <- cgr_strata(panas); brt <- cgr_ratios(bst); bmu <- lapply(bst, mean)
for (cc in c(0.30, 0.50, 0.70)) {
  lhs <- delta_cu(gb$ybar, gb$r, gb$s, 0, cc, 0)
  rhs <- cgr_delta(cc, bmu, brt$r, brt$s)
  rows[[length(rows) + 1]] <- data.frame(
    check = sprintf("Delta(c, 0) = Delta_binary(c), c = %.2f", cc),
    dataset = "PANAS", computed = round(lhs, 6), reference = round(rhs, 6),
    `abs error` = sprintf("%.1e", abs(lhs - rhs)), check.names = FALSE)
}

knitr::kable(do.call(rbind, rows), row.names = FALSE,
  caption = "Both UNKNOWN-extension identities, checked against independent references")
```

| check                                   | dataset             |  computed | reference | abs error |
|:----------------------------------------|:--------------------|----------:|----------:|:----------|
| Delta(c_obs, u_obs) = raw AC - PL       | simulated 6-stratum |  2.724509 |  2.724509 | 0.0e+00   |
| Delta(c_obs, u_obs) = raw AC - PL       | PANAS (no UNKNOWN)  |  3.157042 |  3.157042 | 0.0e+00   |
| Delta(c, 0) = Delta_binary(c), c = 0.30 | PANAS               | -1.711059 | -1.711059 | 0.0e+00   |
| Delta(c, 0) = Delta_binary(c), c = 0.50 | PANAS               |  1.079847 |  1.079847 | 0.0e+00   |
| Delta(c, 0) = Delta_binary(c), c = 0.70 | PANAS               |  3.932466 |  3.932466 | 0.0e+00   |

Both UNKNOWN-extension identities, checked against independent references

## 6.1 What holding $u$ fixed does and does not fix

Pinning $u=u_{obs}$ holds the total mass on UNKNOWN responses fixed but not their within-arm influence. The active-arm denominator is $(1-u)[c(1-r)+(1-c)s]+ut$, whose derivative in $c$ is $(1-u)[(1-r)-s]$. Unless $(1-r)=s$ it changes as $c$ is swept, so the share of the active-arm mean contributed by UNKNOWN responders moves along the curve even with $u$ constant. This is a property of the reweighting, the direct analogue of the arm renormalization already present in the binary case, but the counterfactual should be described precisely: it holds the population share of UNKNOWN responders fixed, not their leverage on either arm mean.

Holding $u=u_{obs}$ and preserving $t$ also assume the UNKNOWN rate and its arm composition are untouched by the reweighting. As with $r$ and $s$ in the original, these define the counterfactual calculation rather than establishing a mechanism.

In small samples, an observed zero cell within an otherwise nonempty response class drives the corresponding ratio to 0 or 1 and gives that cell zero weight. The estimand stays defined but leans entirely on a boundary proportion, so the implementation flags it. An entirely absent response class leaves the reweighting undefined rather than silently invented.

```r
set.seed(3)
du  <- sim_aeb_unknown(n = 240, p_cg = 0.7, u = 0.25, dte_on = TRUE)
uf  <- cgrc_unknown(du, n_draws = 8000, direction = 1, seed = 1)
cnt <- lengths(cgr_unknown_strata(du))[UNKNOWN_STRATA]
zu  <- cgr_unknown_reference_line_test(du)

knitr::kable(
  data.frame(t(cnt), check.names = FALSE),
  caption = sprintf(
    "Six strata (simulated): observed directional CGR = %.3f, UNKNOWN rate = %.0f%%; identity error = %.1e",
    uf$observed_directional_cgr,
    100 * uf$observed_unknown_rate,
    abs(zu$D_at_obs - zu$raw_mean_diff)
)
)
```

| ACAC | ACPL | ACU | PLAC | PLPL | PLU |
|-----:|-----:|----:|-----:|-----:|----:|
|   64 |   32 |  27 |   28 |   58 |  31 |

Six strata (simulated): observed directional CGR = 0.670, UNKNOWN rate = 24%; identity error = 1.8e-15

The six-stratum estimator gets its own operating characteristics, under two assumptions made for the simulation only: the probability of answering UNKNOWN is equal across arms (A1), and an UNKNOWN responder carries no directional expectancy (A2). Under those, the target equals the simulated direct treatment component. Neither assumption is required by the algebra: $\Delta(c,u)$ is computable, and the identities above hold, regardless of whether UNKNOWN reporting is balanced across arms. What depends on A1 and A2 is whether the adjustment recovers the *right* answer, and that turns out to matter more than the identities suggest.

```r
uop <- cgr_unknown_operating(
  n = 230, p_cg = 0.70, u = 0.25,
  mu_dte = 3, mu_aeb = 7.7,
  n_trials = 500, seed = 1
)

u_scen <- c(
  "no effect / no expectancy",
  "effect / no expectancy",
  "no effect / expectancy",
  "effect / expectancy"
)

u_expectancy_only <- uop[uop$DTE == 0 & uop$AEB == 1,, drop = FALSE]
u_n_eff <- uop$n_valid

knitr::kable(
  data.frame(
    scenario = u_scen,
    `true Delta` = round(uop$true, 1),
    `adj bias` = round(uop$adj_bias, 3),
    RMSE = round(uop$adj_rmse, 2),
    `95% coverage` = round(uop$coverage95, 3),
    `coverage MCSE` = round(mcse_prop(uop$coverage95, u_n_eff), 3),
    `flagged, adjusted (P(fav.) > 0.95)` = round(uop$p_fav_gt_95, 3),
    `flagged, unadjusted (P(fav.) > 0.95)` = round(uop$unadj_p_fav_gt_95, 3),
    `empty-stratum rate` = round(uop$empty_stratum_rate, 3),
    check.names = FALSE
),
  caption = "UNKNOWN-preserving operating characteristics, n = 230, directional CGR = 0.70, UNKNOWN rate = 0.25, 500 trials/scenario"
)
```

| scenario                  | true Delta | adj bias | RMSE | 95% coverage | coverage MCSE | flagged, adjusted (P(fav.) \> 0.95) | flagged, unadjusted (P(fav.) \> 0.95) | empty-stratum rate |
|:--------------------------|-----------:|---------:|-----:|-------------:|--------------:|------------------------------------:|--------------------------------------:|-------------------:|
| no effect / no expectancy |          0 |    0.034 | 0.56 |        0.956 |         0.009 |                               0.064 |                                 0.064 |                  0 |
| effect / no expectancy    |          3 |   -0.015 | 1.05 |        0.942 |         0.010 |                               0.878 |                                 0.912 |                  0 |
| no effect / expectancy    |          0 |   -0.042 | 1.06 |        0.948 |         0.010 |                               0.048 |                                 0.650 |                  0 |
| effect / expectancy       |          3 |   -0.081 | 1.43 |        0.934 |         0.011 |                               0.672 |                                 0.980 |                  0 |

UNKNOWN-preserving operating characteristics, n = 230, directional CGR = 0.70, UNKNOWN rate = 0.25, 500 trials/scenario

The largest absolute bias is 0.08 outcome units and coverage ranges from 0.934 to 0.956, a spread comparable to Monte Carlo error. In the expectancy-only scenario the unadjusted analysis flags in 65.0% of valid trials against 4.8% after adjustment. These are operating characteristics under the stated UNKNOWN generative assumptions. They show the estimator behaves coherently when those assumptions hold, not that UNKNOWN responses are arm-independent or expectancy-free in a real trial. Spreading the same sample across six cells also raises sparsity risk, so an `empty_stratum_rate` is reported and `cgr_unknown_min_stratum()` gives the expected-smallest-cell check.

## 6.2 What happens when A1 fails

`cgr_unknown_operating()` cannot test A1 directly, because its own generator holds the UNKNOWN rate arm-independent by construction. Two purpose-built generators check what happens when it is not: one breaks only A1, giving the active and placebo arms different UNKNOWN rates (0.35 versus 0.15) while keeping A2 intact; the other breaks A1 and A2 together, additionally giving UNKNOWN responders their own outcome shift, the case where declining to guess is itself correlated with the outcome. Both hold `n = 230`, `p_cg = 0.70`, `mu_dte = 3`, `mu_aeb = 7.7`, matching the table above, and run the package’s own `cgrc_unknown()` at `c = 0.5`, `u = u_obs` a thousand times against the known true effect of 3.

| Condition                                                                | UNKNOWN rate (AC, PL) | Bias  | 95% coverage |
|--------------------------------------------------------------------------|-----------------------|-------|--------------|
| A1 and A2 both hold (baseline)                                           | 0.25, 0.25            | -0.05 | 0.953        |
| A1 fails, mild imbalance                                                 | 0.30, 0.20            | -0.36 | 0.949        |
| A1 fails, larger imbalance (A2 still holds)                              | 0.35, 0.15            | -0.82 | 0.901        |
| A1 and A2 both fail (larger imbalance, plus an outcome-correlated shift) | 0.35, 0.15            | -2.02 | 0.681        |

Halving the arm gap in UNKNOWN rate, from a 20-point gap (0.35 vs 0.15) to a 10-point gap (0.30 vs 0.20), roughly halves the bias (-0.82 to -0.36) and returns coverage nearly to nominal (0.901 to 0.949). That is close to linear in the size of the imbalance, which is informative: this is not a threshold effect where the adjustment is fine until some tipping point and then fails abruptly. It degrades gracefully and in proportion to how unbalanced the UNKNOWN rates actually are. Adding a real dependence between answering UNKNOWN and the outcome, on top of the larger imbalance, roughly doubles the damage again: bias exceeds two thirds of the true effect and the interval covers the truth only about two trials in three.

The likely mechanism: holding a single pooled `u` and a single pooled `t` (the arm split among UNKNOWN responders) and moving only `c` to 0.5 implicitly treats the UNKNOWN summary as applying symmetrically to both arms. When the arms actually lose different fractions of their AEB-affected participants to UNKNOWN, the remaining directional strata in each arm are purified by different amounts, and reweighting toward `c = 0.5` on that basis does not cleanly recover the direct effect. The near-linear scaling with the arm gap is consistent with this account, though it has only been checked at these three imbalance levels and at one sample size, so it should still be read as a well-supported working explanation rather than a proof.

The practical implication is concrete: **if UNKNOWN response rates differ visibly between arms, the default adjustment should not be trusted without checking this.** `t`, the arm split among UNKNOWN responders, is reported by the app and is the diagnostic to inspect; a `t` far from the trial’s overall arm balance is the signal that this simulation’s assumptions may not hold.

# 7 Three independent trials

Each trial below runs through the same pipeline: audit, detect whether guesses are binary or include UNKNOWN, then apply the four- or six-cell analysis. They were chosen for the regime each illustrates.

| Trial                       | Design                        | Guess structure       | Regime                             |
|-----------------------------|-------------------------------|-----------------------|------------------------------------|
| Cavanna et al. (2022)       | psilocybin microdosing        | binary, reconstructed | large effect, guessing near chance |
| Santana-Penín et al. (2023) | sham dental therapy           | AC / PL / UNKNOWN     | a third answered “I don’t know”    |
| Lii et al. (2023)           | ketamine masked by anesthesia | AC / PL / UNKNOWN     | directional guessing at chance     |

Two measurement caveats matter more than the results. For **Cavanna** the guess is not observed: the public data record only whether the participant was correct, so the guessed arm is reconstructed as the received arm when that indicator is 1 and the opposite when 0. This is exact only under a forced active-versus-placebo choice with no third option; if participants could decline, the reconstruction converts those into directional guesses, the failure the UNKNOWN extension exists to avoid. For **Santana** the guess used is the one collected at six months, alongside the outcome, and for **Lii** belief was recorded at Closeout and applied to an outcome averaged over Days 1 to 3. In both, the reweighting conditions on a variable measured after, and plausibly caused by, the outcome it adjusts. That is the causal objection from the opening section in concrete form.

Outcomes were also selected to illustrate a regime rather than being registered primary endpoints. That matters most for Cavanna, where the outcome is a dosing-day sensation VAS derived from the raw daily ratings, chosen as a deliberately high-signal measure.

```r
cav <- read_ext("cavanna_cgrc.csv")
san <- read_ext("santana_cgrc_unk.csv")
lii <- read_ext("ketamine_cgrc_unk.csv")

# The shipped `vas_dose_days` column is the four-day mean (Wed/Thu/Fri/Sat), not
# the dosing-day mean, so it is deliberately not used. Derive the dosing-day
# average (Wednesday and Friday) from the raw daily ratings so the column name
# and the text describing it agree.
cav$vas_mean_dose_days <- rowMeans(
  cav[, c("acute_effect_vas_wed", "acute_effect_vas_fri")]
)

A <- list(
  Cavanna  = app_analyze(cav, "treatment_received", "treatment_guessed", "vas_mean_dose_days",      direction =  1),
  Santana  = app_analyze(san, "treatment_received", "treatment_guessed", "pain_improvement_6m_itt", direction =  1, delta = 1.5),
  Lii      = app_analyze(lii, "treatment_received", "treatment_guessed", "madrs_mean_day1_3",       direction = -1))
row_of <- function(nm, a) {
  s <- a$fit$summary
  ocg <- if (a$mode == "unknown") a$fit$observed_directional_cgr else a$fit$observed_cgr
  pct <- function(p) sprintf("%.0f%%", 100 * p)
  data.frame(trial = nm, mode = a$mode, n = nrow(a$trial),
             `obs CGR` = round(ocg, 3),
             raw = round(s$post_mean[1], 2), `adjusted (0.50)` = round(s$post_mean[2], 2),
             `adj 95% CrI` = sprintf("[%.2f, %.2f]", s$cri_lo[2], s$cri_hi[2]),
             `P(fav) raw >>> adj` = paste(pct(s$p_favourable[1]), ">>>", pct(s$p_favourable[2])),
             `P(meaningful) raw >>> adj` = paste(pct(a$head$p_meaningful_obs), ">>>", pct(a$head$p_meaningful_blind)),
             check.names = FALSE)
}
knitr::kable(do.call(rbind, Map(row_of, names(A), A)), row.names = FALSE,
  caption = "Three independent trials through the app's analysis pipeline (Cavanna outcome = dosing-day VAS, a selected high-signal outcome)")
```

| trial   | mode    |   n | obs CGR |   raw | adjusted (0.50) | adj 95% CrI     | P(fav) raw \>\>\> adj | P(meaningful) raw \>\>\> adj |
|:--------|:--------|----:|--------:|------:|----------------:|:----------------|:----------------------|:-----------------------------|
| Cavanna | binary  |  34 |   0.559 | 21.55 |           19.93 | \[0.38, 39.22\] | 99% \>\>\> 98%        | 81% \>\>\> 75%               |
| Santana | unknown |  77 |   0.569 |  1.58 |            1.43 | \[0.49, 2.41\]  | 100% \>\>\> 100%      | 57% \>\>\> 45%               |
| Lii     | unknown |  38 |   0.519 | -3.02 |           -2.80 | \[-8.36, 2.63\] | 87% \>\>\> 86%        | 28% \>\>\> 25%               |

Three independent trials through the app’s analysis pipeline (Cavanna outcome = dosing-day VAS, a selected high-signal outcome)

**Cavanna** guesses only a little above chance (0.559), so the adjustment barely moves the estimate; the interval is wide with 34 participants but stays above zero. The trial’s own conclusion, that acute VAS effects were larger for active than placebo only among correct guessers, is informally the four-stratum split this method formalizes; the CGRC adds a continuous reweighting and a posterior in place of a single subgroup comparison.

**Santana** stays strongly positive after adjustment, with $P(\text{favorable})$ at 100% and an interval entirely above zero, but clinical meaningfulness is the open question: against the trial’s 1.5-point threshold the probability slips from 57% to 45% because the adjusted estimate (1.43) lands just under the bar. It is also the live example of the feasibility caution: keeping the 34% who answered “I don’t know” is why all 77 remain in the analysis, but the discordant cells are thin, one holding 1 participant. That number is defined but fragile.

*What the trial claims, and what this adds.* Santana-Penín et al. reported an adjusted six-month pain reduction of 1.54 points relative to sham, just over their own 1.5-point threshold for clinical importance, and concluded that the active therapy significantly reduced pain and improved mouth opening, with no qualification tied to blinding or expectancy. The published analysis never used the guess data for this outcome, despite collecting assignment-awareness responses at both post-therapy and six months. The CGRC analysis is the first use of that data here. It leaves the direction of benefit intact but moves the clinical-importance question from just above the trial’s own bar to just below it once directional guessing is reweighted to chance, and it flags that the reweighting itself is fragile in this trial: the `ACPL` stratum, guessed placebo but received the real therapy, holds a single participant.

**Lii** has directional guessing essentially at chance (0.519), so the adjustment is nearly inert and all six cells fill comfortably. The trial itself reran its MADRS results with participants regrouped by guess alone, discarding randomized assignment; the six-stratum CGRC instead crosses guess with assignment and preserves the randomized contrast. With guessing this close to chance the two land in the same place, which is what should happen when guesses carry little information.

```r
app_curve(A$Santana, "favorable") +   # UNKNOWN preserved; effect at the clinical MID
  ggplot2::labs(title = "Santana-Pen\u00edn et al. (2023): dental therapy") +
  ggplot2::theme(text = ggplot2::element_text(size = 15))
```

![Santana-Penín et al. (2023), sham dental therapy: pain improvement at 6 months, UNKNOWN preserved](https://raw.githubusercontent.com/dsmi313/cgrc.bayes/main/report-files/summary/fig02.png)

*Santana-Penín et al. (2023), sham dental therapy: pain improvement at 6 months, UNKNOWN preserved*

Santana-Penín et al. (2023), sham dental therapy: pain improvement at 6 months, UNKNOWN preserved

```r
app_curve(A$Lii, "favorable") +       # near-chance guessing; correction is small
  ggplot2::labs(title = "Lii et al. (2023): ketamine masked by anesthesia") +
  ggplot2::theme(text = ggplot2::element_text(size = 15))
```

![Lii et al. (2023), ketamine masked by anesthesia: MADRS depression score, directional guessing near chance](https://raw.githubusercontent.com/dsmi313/cgrc.bayes/main/report-files/summary/fig03.png)

*Lii et al. (2023), ketamine masked by anesthesia: MADRS depression score, directional guessing near chance*

Lii et al. (2023), ketamine masked by anesthesia: MADRS depression score, directional guessing near chance

Participant counts, derived variables, reconstruction assumptions and stratum breakdowns for all three datasets are recorded in [PROVENANCE.txt](https://github.com/dsmi313/cgrc.bayes/blob/main/inst/extdata/external_studies/PROVENANCE.txt).

# 8 What this does and does not claim

The estimand is unchanged from the original, as is the self-test identity. The Bayesian formulation adds a directly interpretable uncertainty summary in place of an averaged p-value. The simulations add a way to see, before applying it, whether the adjustment is approximately unbiased, approximately calibrated, and even computable for a given design. The UNKNOWN extension adds a fully specified six-cell reweighting that keeps genuine “don’t know” answers, with its own identities, operating characteristics and sparsity diagnostics.

What none of it establishes:

- Shrinkage under adjustment is not a measure of how much of an effect was expectancy. Randomization identifies the total effect of assignment; guessing happens afterward, so the difficulty is mediation rather than confounding, and reweighting on guesses does not isolate a pharmacological pathway. Loewinger et al. (2026, Appendix E.4.1) derive the causal quantity the CGRC weights correspond to and show it need not equal a controlled direct effect, since it is built from belief-stratified means and so inherits their exposure to collider bias. Szigeti et al. (2023) anticipate a version of this in their own limitations; the Loewinger analysis shows the exposure is more general than that conditional, and gives a counterexample in which the adjusted contrast carries the opposite sign from the controlled direct effect.
- Reproducing the original validates the computation, not its causal assumptions. The identity is a bug check, not a correctness proof.
- The simulations validate behavior within their own model.
- The UNKNOWN extension fixes a data-handling problem, not the causal limits.

Read as a fragility check rather than a decomposition, the CGRC does something no other available tool does: it puts a number and an uncertainty band on how much a trial result depends on participants having been blind. The package makes that check runnable on a researcher’s own data, with a feasibility diagnostic beforehand and honest uncertainty after, via `cgrc.bayes::cgrc_app()`.

# 9 Sources

- Szigeti, B., Nutt, D., Carhart-Harris, R., & Erritzoe, D. (2023). The difference between ‘placebo group’ and ‘placebo control’: A case study in psychedelic microdosing. *Scientific Reports*, *13*, 12107. <https://doi.org/10.1038/s41598-023-34938-7>
- Loewinger, G., Stensrud, M. J., Nayak, S. M., Yaden, D., & Levis, A. W. (2026). Causal inference in studies with functional unmasking: Psychedelics and beyond \[Preprint\]. medRxiv. <https://doi.org/10.64898/2025.12.05.25341713>
- Szigeti, B., Kartner, L., Blemings, A., Rosas, F., Feilding, A., Nutt, D. J., Carhart-Harris, R. L., & Erritzoe, D. (2021). Self-blinding citizen science to explore psychedelic microdosing. *eLife*, *10*, e62878. <https://doi.org/10.7554/eLife.62878>
- Cavanna, F., Muller, S., de la Fuente, L. A., Zamberlan, F., Palmucci, M., Janeckova, L., Kuchar, M., Pallavicini, C., & Tagliazucchi, E. (2022). Microdosing with psilocybin mushrooms: A double-blind placebo-controlled study. *Translational Psychiatry*, *12*, 307. <https://doi.org/10.1038/s41398-022-02039-0> Data: <https://doi.org/10.5281/zenodo.5745892>
- Santana-Penín, U., Santana-Mora, U., López-Solache, A., Mora, M. J., Collier, T., Pocock, S. J., Lorenzo-Franco, F., Varela-Centelles, P., & López-Cedrún, J. L. (2023). Remodeling dental anatomy vs sham therapy for chronic temporomandibular disorders. A placebo-controlled randomized clinical trial. *Annals of Anatomy*, *250*, 152117. <https://doi.org/10.1016/j.aanat.2023.152117> Data: <https://doi.org/10.5061/dryad.zkh189370>
- Lii, T. R., Smith, A. E., Flohr, J. R., Okada, R. L., Nyongesa, C. A., Cianfichi, L. J., Hack, L. M., Schatzberg, A. F., & Heifets, B. D. (2023). Randomized trial of ketamine masked by surgical anesthesia in patients with depression. *Nature Mental Health*, *1*(11), 876-886. <https://doi.org/10.1038/s44220-023-00140-x> Data: <https://doi.org/10.17605/OSF.IO/ZDKR8>
- Szigeti, B., Weiss, B., Rosas, F. E., Erritzoe, D., Nutt, D., & Carhart-Harris, R. (2024). Assessing expectancy and suggestibility in a trial of escitalopram v. psilocybin for depression. *Psychological Medicine*, *54*, 1717-1724. <https://doi.org/10.1017/S0033291723003653>
- Szigeti, B., & Heifets, B. D. (2024). Expectancy effects in psychedelic trials. *Biological Psychiatry: Cognitive Neuroscience and Neuroimaging*, *9*, 512-521. <https://doi.org/10.1016/j.bpsc.2024.02.004>
- Norman, G. R., Sloan, J. A., & Wyrwich, K. W. (2003). Interpretation of changes in health-related quality of life: The remarkable universality of half a standard deviation. *Medical Care*, *41*(5), 582-592.
- Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian data analysis* (3rd ed.). CRC Press. Chapter 3 covers conjugate inference for a normal mean and variance.

Package, code and test suite: <https://github.com/dsmi313/cgrc.bayes>. Per-finding citations: <https://github.com/dsmi313/cgrc.bayes/blob/main/reports/SOURCES.md>. Per-dataset provenance and data-preparation notes: <https://github.com/dsmi313/cgrc.bayes/blob/main/inst/extdata/external_studies/PROVENANCE.txt>.

This project reproduces and re-analyzes published work; it is independent of and not endorsed by the original authors.
