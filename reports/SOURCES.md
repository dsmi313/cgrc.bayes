# Sources

The literature this project reproduces and cites. The PDFs themselves are
copyrighted journal articles and are **not** committed to the repo; cite them by
DOI. Read 2026-07-21.

**Source roles at a glance** (only the first two underpin the package):

1. **Primary methodological source** — Szigeti et al. 2023, *Sci Rep* (S1 below):
   the CGRC method reproduced here.
2. **Empirical dataset** — Szigeti et al. 2021, *eLife* ("The empirical dataset"
   below): the source study for `data/pacutes.csv`, reanalysed here.
3. **Supporting literature** — Szigeti & Heifets 2024, *Biol Psychiatry CNNI*
   (S3 below): contextual interpretation, not a CGRC methodological source.
4. **Related expectancy literature** — Szigeti et al. 2024, *Psychological
   Medicine* (S2 below): context only; not the CGRC method or the data used here.

The author's source code (S4) and the pinned dataset are described after the
papers. Section labels S1–S4 are stable identifiers referenced from
`reports/UNRESOLVED.md`; they are ordering tags, not a ranking of importance.

## S1 — the CGRC method paper (primary methodological source)

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
  to the trial's original CGR (= 0.72)". The source code (S4) shows this is a
  **hardcoded constant** (`trial_cgrs = {'sbmd': 0.72}`), not the data's CGR
  (0.647); its value most likely originates from the placebo-arm rate (0.7234),
  which is the one part still a hypothesis.
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

## S2 — related expectancy literature (context only, different trial)

Szigeti B, Weiss B, Rosas FE, Erritzoe D, Nutt D, Carhart-Harris R. **"Assessing
expectancy and suggestibility in a trial of escitalopram v. psilocybin for
depression."** *Psychological Medicine* 54, 1717-1724 (2024).
<https://doi.org/10.1017/S0033291723003653>

A different trial (escitalopram vs psilocybin for depression), **not** the
microdose CGRC and **not** a source of the CGRC method or the data used here.
Included only because it develops the same expectancy/blinding-failure theme
(e.g. it notes 94% correct guessing in another psilocybin trial). It does not
bear directly on any item in UNRESOLVED.md.

## S3 — supporting literature (contextual, not a CGRC methodological source)

Szigeti B, Heifets BD. **"Expectancy Effects in Psychedelic Trials."**
*Biological Psychiatry: Cognitive Neuroscience and Neuroimaging* 9:512-521
(May 2024). <https://doi.org/10.1016/j.bpsc.2024.02.004>

Same-author review; supports *interpretation* of expectancy effects, blinding
failure and the approximate correct-guess rates, but does not provide the CGRC
method or estimand. Used here for two facts that bear on U3 and U6:

- The microdose correct-guess rate is "only ~65% to 70%" (a truly blind trial is
  ~50%) - consistent with the public data's 0.647, not with the hardcoded 0.72
  in the 2023 Fig 4 line. (The code, S4, is the primary evidence that 0.72 is a
  fixed constant; this review just corroborates that the real rate is ~65-70%.)
- Microdose effects fall below the 0.5-SMD minimally-important difference,
  "too small to be noticeable" - reinforces the U6 magnitude argument.

## S4 — the author's source code (located 2026-07-21)

Szigeti B. **CorrectGuessRateCurve** (analysis code for S1). MIT licensed,
DrugNerdsLab, 2022. <https://github.com/szb37/CorrectGuessRateCurve>

Named in S1's data-availability statement; earlier turns of this project probed
`szb37/mcrds_public` (the *data* mirror) and reported the source "not located".
It was in the paper the whole time. Reading it closes two open items and
corrects two document claims:

- **The 0.72 is hardcoded (closes U3).** `src/config.py`:
  `trial_cgrs = {'sbmd': 0.72}` — a fixed constant, drawn as the Figure 4
  reference line. The same code computes the trial CGR from data as
  `(n_plpl + n_acac) / n` = 0.647 (`src/cgrc/core.py` line ~181), but that
  computed value is not what the plotted line uses. So 0.72 is an annotation,
  not the data's CGR; it coincides with the placebo-arm rate 0.7234.
- **The estimand is independently confirmed (closes U4).**
  `get_strata_ratio` / `get_strata_sample_sizes` in `src/cgrc/core.py` form
  `r = PLPL / (PLPL + ACAC)` and `s = ACPL / (ACPL + PLAC)` — identical to this
  implementation — and the default `strata_sampling = 'all_prop'` confirms the
  `noise = "all"` reading (CH-04).
- **Resample count (qualifies "100 times").** The repo's Figure-4 config
  `cgrC_low` (`config.py`) specifies `n_cgrc_trials = 32` over
  `np.linspace(0, 1, 13)` (options 32/64/96, never 100). Two caveats: the repo's
  README says it reproduces the **preprint** ("On the fallibility of the placebo
  control...", a different title, then under review), so the config may be
  preprint-era; and the Figure-4 block in `run.py` sits inside `if False:`. So
  this is the repository's stated configuration, not proof of what generated the
  final published figure. If the real count is 32, the Monte Carlo error is
  ~1.8x a 100-resample assumption.
- **Rounding (documents the reproduction default).** `get_strata_ratio` does
  `round(x, 2)` on every stratum proportion, so `legacy_round = TRUE` is the
  faithful reproduction path (+0.010 PANAS, +0.019 Energy vs the exact ratios).

## The empirical dataset (reanalysed here)

The source study for the public data reanalysed in the empirical reproduction:
Szigeti B, et al. **"Self-blinding citizen science to explore psychedelic
microdosing."** *eLife* 10:e62878 (2021).
<https://doi.org/10.7554/eLife.62878>

The package pins a checksum-verified public copy as `data/pacutes.csv` (mirrored
at `szb37/mcrds_public`) for reproducibility. Provenance and SHA-256 in
`data/PROVENANCE.txt`.
