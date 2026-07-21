# The model in plain language

Read this and then try to explain it back without looking. The checkpoints at
the end are the parts people usually get wrong.

## The problem

You run a trial. Half get the drug, half get placebo. At the end you ask
everyone to guess which they got - and 65% guess right. That is a problem,
because people who think they got the real drug tend to feel better *because*
they think that. Since more of the drug group correctly believe they got the
drug, the "feeling better because you believe" bonus is spread unevenly. Some
of your apparent drug effect is really that unevenness.

## The idea

Ask a hypothetical: what would this trial have shown if only 50% had guessed
right - which is what you would see if the blinding had actually worked?

You cannot rerun the trial. But you can reweight the people you already have.

## The four boxes

Everyone lands in one of four boxes: what they GOT crossed with what they
GUESSED.

|  | Guessed drug | Guessed placebo |
|---|---|---|
| **Got drug** | ACAC (right) | ACPL (wrong) |
| **Got placebo** | PLAC (wrong) | PLPL (right) |

Two boxes are correct guesses, two are wrong guesses. The correct guess rate is
just the share of people in the two "right" boxes.

## The reweighting

To simulate 50% correct guessing, put half the total weight on the two "right"
boxes and half on the two "wrong" boxes.

But there is a trap. If you are careless, you also change the drug-to-placebo
balance - and then you would be comparing a different trial entirely. So you
hold two proportions fixed:

- **r**: among people who guessed right, what fraction were on placebo?
- **s**: among people who guessed wrong, what fraction were on the drug?

Keeping these fixed means the reweighting changes only *how much guessing
happened*, never *who got what*. That is the whole trick.

Then you compute the drug group's average and the placebo group's average using
the new weights, and take the difference. That is Delta(0.5).

## The check that proves it works

Set the target to the ACTUAL observed guess rate rather than 0.50. Then the
reweighting should do nothing at all - and the answer should be exactly the
ordinary difference in averages you would have got without any of this.

It is. Exactly, to fifteen decimal places. If it were not, something would be
broken.

## Where the Bayesian part comes in

So far this uses each box's *average*. But an average from 40 people is
uncertain. The Bayesian step replaces each box's single average with a **cloud
of plausible averages**, given the data.

Here is the part worth getting right. For each box we do not just ask "what is
the average" - we ask about the average AND the spread together, because how
uncertain you are about an average depends on how spread out the data are.

So each draw works in two steps:

1. Draw a plausible spread for this box.
2. Given that spread, draw a plausible average.

Do that 20 000 times per box. Then for each of the 20 000 sets, run the
reweighting arithmetic. You end up with 20 000 plausible values of Delta(0.5) -
a distribution, not a number.

From that distribution: the middle is your estimate, the 2.5th to 97.5th
percentiles are your 95% credible interval, and the fraction above zero is the
probability the effect is positive.

## What the number at the bottom of the plot is

It is **the probability that the treatment effect is positive**, given the
model and the data.

It is NOT a p-value. A p-value answers "if there were no effect, how surprising
is this data?" This answers "given this data, how likely is it there IS an
effect?" Those are different questions and the second is the one people
usually think they are asking.

## Three checkpoints

**1. Does running 20 000 draws instead of 1 000 give you more information about
the patients?**

No. There are 232 people either way. More draws only makes the *summary of the
distribution* more precise - it smooths out the computer's own randomness. It
tells you nothing new about anyone in the trial.

**2. If the adjusted effect is smaller than the unadjusted one, does that
measure how much was caused by expectancy?**

No, and this is the easiest mistake to make. It measures what happens to the
estimate under a reweighting assumption. Calling the gap "the expectancy
effect" assumes the whole causal story is right.

**3. What if people guessed correctly BECAUSE they got better?**

Then the entire method backfires. The logic assumes people guess right because
of side effects, and then feeling better follows from the belief. If instead
they felt better first and inferred the drug from that, then adjusting away the
guessing also adjusts away a genuine effect - and you would produce a false
negative rather than removing a false positive.

Nothing in this data can tell you which of those is happening. The original
paper says so, and that warning matters more than any of the arithmetic.
