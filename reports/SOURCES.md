# Sources

The primary and secondary literature this project reproduces and cites. The PDFs
themselves are copyrighted journal articles and are **not** committed to the
repo; cite them by DOI. Read 2026-07-21.

## S1 — the CGRC method paper (primary)

Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D. **"The difference between
'placebo group' and 'placebo control': a case study in psychedelic
microdosing."** *Scientific Reports* 13:12107 (2023).
<https://doi.org/10.1038/s41598-023-34938-7>

This is the source of everything reproduced here: the activated expectancy bias
(AEB) model (Fig 1, Eq 1-4), the Correct Guess Rate Curve (Fig 2), Table 1 (the
simulation operating characteristics), Table 2 (the microdose re-analysis), and
Fig 4 (the empirical CGR curves). Facts confirmed against it during this work:

- **Analysis rule (U2).** "We only used data from the first week... each
  datapoint is independent... In the current analysis n = 233 datapoints were
  included." Confirms the `tp == "w1s1"` filter; the public file gives 232.
- **The 0.72 CGR (U3).** Fig 4 caption: "vertical green dashed line corresponds
  to the trial's original CGR (= 0.72)". Shown here to be the placebo-arm rate
  (0.7234), not the overall CGR (0.647).
- **`noise = "all"` (U4).** Table 1 rates 0.05 / 0.86 / 0.78 / 0.99 and DTE
  Hedges g = 0.4; stated regime "n ~ 200, CGR ~ 0.7, effect ~ 0.4 Hedges' g".
  All match `noise = "all"`. The exact Eq. 4 SD scope is in the paper's
  Supplementary Table 1 (not in the supplied PDFs).
- **Malicious-unblinding argument (U6).** 55% body/perceptual cues (muscle
  tension 58%, stomach discomfort 27%) vs 23% mental/psychological benefits;
  placebo-microdose PANAS diff 2.1/0.8 vs ~10/~6 natural within-subject
  variability.
- **Table 2 reproduction.** PANAS 3.2 -> 1.1, Mood VAS 6.4 -> 2.7, Energy VAS
  11.5 -> 6.8 (g 0.58 -> 0.34), only Energy VAS surviving CGR adjustment - all
  reproduced from the public data.

## S2 — companion expectancy paper (context, different trial)

Szigeti B, Weiss B, Rosas FE, Erritzoe D, Nutt D, Carhart-Harris R. **"Assessing
expectancy and suggestibility in a trial of escitalopram v. psilocybin for
depression."** *Psychological Medicine* 54, 1717-1724 (2024).
<https://doi.org/10.1017/S0033291723003653>

A different trial (escitalopram vs psilocybin for depression), not the microdose
CGRC. Included because it develops the same expectancy/blinding-failure theme
(e.g. it notes 94% correct guessing in another psilocybin trial). It does not
bear directly on any item in UNRESOLVED.md.

## S3 — review (secondary)

Szigeti B, Heifets BD. **"Expectancy Effects in Psychedelic Trials."**
*Biological Psychiatry: Cognitive Neuroscience and Neuroimaging* 9:512-521
(May 2024). <https://doi.org/10.1016/j.bpsc.2024.02.004>

Same-author review. Used here for two facts that bear on U3 and U6:

- The microdose correct-guess rate is "only ~65% to 70%" (a truly blind trial is
  ~50%) - consistent with the public data's 0.647, not with the 0.72 in the 2023
  Fig 4 caption. Supports the reading that 0.72 is the placebo-conditional rate.
- Microdose effects fall below the 0.5-SMD minimally-important difference,
  "too small to be noticeable" - reinforces the U6 magnitude argument.

## The public dataset

`data/pacutes.csv`, from the self-blinding microdose trial (Szigeti et al.,
*eLife* 10:e62878, 2021), mirrored at `szb37/mcrds_public`. Provenance and
SHA-256 in `data/PROVENANCE.txt`.
