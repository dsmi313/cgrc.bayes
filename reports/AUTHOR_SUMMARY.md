# Note for Balazs Szigeti

I have been working through the CGRC method in detail, partly to understand it
properly and partly to see whether the uncertainty around the adjusted estimate
could be expressed as a posterior distribution rather than through resampling.
The estimand is unchanged throughout - this is an alternative way of computing
uncertainty around your quantity, not a different quantity.

Three things came out of it that seem worth sending you.

**The reproduction works.** Using the public week-1 data I get 3.16 unadjusted
and 1.08 at CGR 0.50 for PANAS, against your published 3.2 and 1.1; energy VAS
gives 11.38 and 7.10 against 11.5 and 6.8. I also reimplemented the KDE
resampling procedure directly, and when I run it at 10 000 resamples instead of
100 it converges on the closed-form value in every cell I tried. The averaged
p-values come out close to your published ones too (PANAS 0.41 vs 0.43, energy
0.043 vs 0.04), which is what convinced me the port was faithful rather than
just plausible.

That comparison also suggests the KDE step is doing less work than it might
appear. Because the estimand only uses stratum means and kernel smoothing
preserves means, the KDE and a Gaussian stratum model give the same answer; the
visible difference between 100 resamples and the closed form is Monte Carlo
noise, worth roughly 0.12-0.29 points on these scales. So the honest summary of
what a Bayesian version adds is fairly narrow: it removes that noise, and it
gives a posterior probability instead of an average of p-values. It is not
better inference about the trial, and I have tried to keep the writeup from
implying otherwise.

**The 0.72, resolved from your code.** I eventually found
`CorrectGuessRateCurve` (it was in your data-availability statement - I had been
looking in `mcrds_public`). The finding is four lines from your own repo:

```python
# core.py:181  — feeds the comparison table
trial_cgr = (n_plpl+n_acac)/tmp_trial_data.shape[0]   # = 0.647, rounds to 0.65

# figures.py:493 — draws the green line (cgr <- config.trial_cgrs['sbmd'] = 0.72)
ax1.axvline(x=cgr, color='green', ls='--')  # Vertical line, empirical CGR
```

The same run computes the empirical CGR as 0.647 for the table, but the green
line is drawn at the hardcoded constant 0.72 - and the comment on that line even
calls it "empirical CGR". So the annotation and the computed value disagree.
0.72 happens to equal the placebo-arm correct-guess rate (0.7234; the microdose
arm is only 0.5275), which is probably where the constant came from. It does not
affect any estimate, but Figure 4's reference line is not the trial's guess rate,
so it is worth a correction or a footnote.

**One thing I still could not resolve.** You report n = 233 for week 1; the
public file gives me 232 on every acute scale, with no duplicate trial IDs and
no missing outcome values. I could not find a filtering rule that produces 233 -
one record short.

**A pattern that runs opposite to your 2024 review.** The review says correct
guess rates are generally higher in the active arms. Here it is the other way:
placebo guesses correctly 72% of the time and microdose only 53%, because guesses
skew heavily toward "placebo" (9508 vs 6115 across the set). It does not change
the method - the curve conditions on the strata through r and s - but it flips
the intuition about which arm expectancy inflates, so I wanted to flag it.

**One ambiguity in the AEB equations.** Equation 4 does not say whether the SD
of the DTE and AEB terms applies to the whole sample or only to the gated
subgroup. I simulated both. The whole-sample reading gives Hedges' g = 0.4011
against your stated 0.40, and unadjusted significance rates of
0.056/0.876/0.826/0.992 against your Table 1's 0.05/0.86/0.78/0.99; the
subgroup reading gives g = 0.50 and rates that miss on two rows. So I have gone
with the whole-sample reading, but I could not find the simulation code in the
public repo, so this is inference from your published numbers rather than
confirmation. If the code is available anywhere I would like to check it
properly.

I also ran the adjusted estimator over 500 simulated trials per scenario. It
came out unbiased in all four AEB configurations and the 95% intervals covered
at 0.936-0.960, which was reassuring. The cost is power: in the
partial-mediation scenario the adjusted analysis flags a favourable effect 70%
of the time against 89% for the unadjusted one. That seems like a reasonable
trade for the false-positive rate dropping from 0.83 to 0.06 when expectancy is
the only thing driving the result, but it is a real cost and I have reported it
as one.

Everything is in a repo with the data pinned by checksum and the code split out
of the writeup. Happy to send it over, and happy to be told I have
misunderstood something - the causal assumption about malicious versus benign
unblinding in particular is one I have restated in your terms rather than
tested, since I do not think these data can settle it.
