<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# ALDEx2: behavioural specification

This is the specification the native `aldex2` is written from. It describes
**what** the reported numbers are, as measurable behaviour, so that an
implementation can be written from it alone.

## Provenance

The method is from Fernandes et al. 2013 (PLoS ONE 8:e67019) and 2014
(Microbiome 2:15). The papers do not pin down two details that decide the
reported numbers:

- how the per-instance p-values are combined into an expected p-value;
- how the random pairs behind the effect size are drawn.

Bioconductor ALDEx2's public documentation describes the first only loosely.
NEWS for 1.33 says "all p-values calculated are now posterior p-values with
consistent sign", and the vignette calls `we.ep` "a posterior predictive
p-value".

On 2026-10-02 the native and Bioconductor 1.42.0 results disagreed beyond
Monte-Carlo noise. To find out why, the coding agent (Claude, session
https://claude.ai/code/session_019Y6e175rh1LbrsEhQvCrJ1) read the reference
package's R source (GPL ≥ 3): `aldex.ttest`, `t.fast`, `wilcox.fast` and
`aldex.effect`. No code was copied. The behaviour below was then confirmed
numerically on the reference package's own Monte-Carlo instances: on the mock,
gut and soil datasets, the p-value rules reproduce its `we.ep`, `we.eBH` and
`wi.ep` to within 1e-14.

The owner ruled that the native code is written from this document by someone
who has **not** read the reference source. That is what happened: a separate
agent that was barred from the reference source, R and web search wrote
`src/Methods/ALDEx2.jl` from this document (spec committed alone as `87a6316`,
before the code).

## Notation

- Input: integer counts with `n` samples and `d` features, and a two-level
  condition per sample.
- `L1 < L2` are the two condition levels in sorted order.
- `inst[i, j, k]` is the CLR value of sample `i`, feature `j`, in Monte-Carlo
  instance `k = 1..mc`. These are produced by the existing
  `dirichlet_clr_instances`, which is unchanged by this document.
- **Test group 1** is the condition of the *first sample in input order*.
  **Test group 2** is the other condition. With `n1` and `n2` samples
  respectively, `n1 + n2 = n`.

**Calling convention.** The function is `aldex2_from_instances(inst, i1, i2; rng, first_is_i1 = true)`.

- `i1` and `i2` are the sample rows of `L1` and `L2`.
- `first_is_i1` says whether the first input sample belongs to `L1`. If it does,
  `i1` is test group 1; if not, `i2` is.
- `aldex2` computes `first_is_i1` from `groups[1]`.
- The property test P6 calls the function without this flag and gets the
  default.

Test group 1 only matters for the exact Wilcoxon p-value, because that
distribution is discrete. The rule is part of the reference behaviour: on the
soil dataset, whose first sample is in `L2`, using `L1` instead moves `wi_pvalue`
by up to 0.0097.

## Per-instance one-sided p-values

For each instance `k` and feature `j`, let `a` be the values of test group 1 and
`b` the values of test group 2.

- **Welch.** `t = (mean(a) − mean(b)) / sqrt(s_a²/n1 + s_b²/n2)`, using sample
  variances. The degrees of freedom are Welch–Satterthwaite. The upper-tail
  p-value is `pW = P(T_df ≥ t)`.
- **Wilcoxon rank sum.** `U = #{(x, y) : x ∈ a, y ∈ b, x > y}`.
  - If `n1 < 50` **and** `n2 < 50`: `pR = P(U' ≥ U)` under the exact null
    distribution of the Mann–Whitney statistic for sizes `(n1, n2)`.
  - Otherwise use the normal approximation **without continuity correction**:
    `pR = 1 − Φ((U − n1·n2/2) / sqrt(n1·n2·(n1+n2+1)/12))`.
  - CLR instances are continuous, so ties have probability zero. If any do
    occur, count a tie as ½ in `U`; the result in that case is not
    oracle-checked.

## Consistent-sign expected p-values

For a per-instance upper-tail p-value `p` (either `pW` or `pR`), define the two
doubled directional values:

```
g = min(1, 2p)          # "group 1 higher"
l = min(1, 2(1 − p))    # "group 1 lower"   (note: 1 − p, not the lower tail)
```

Then, per test:

- **Expected p:** `min( mean_k g_k , mean_k l_k )`. This is reported as
  `pvalue` for Welch and `wi_pvalue` for Wilcoxon.
- **Expected BH:** within each instance, BH-adjust the vector of `g` over
  features, and separately the vector of `l` over features, using the existing
  `bh_adjust`. Average each over instances, then take the elementwise minimum of
  the two averages. This is reported as `qvalue` and `wi_qvalue`.

This is not the mean of two-sided p-values: by Jensen, the mean of the minimum
is at most the minimum of the means. The two agree when a feature's direction is
the same in every instance.

## Effect size

The effect is computed on the **sorted** levels, so its sign is `L2 − L1`. The
test-group rule does not apply here.

For each feature `j`:

- `x1` pools all instance values of `L1` samples, `n_L1·mc` values in total.
  `x2` pools those of `L2` samples.
- `m = min(length(x1), length(x2))`.
- **Within-group spread.** For each group `g`, take two independent uniformly
  random permutations of `x_g` and keep the first `m` elements of each. `win_g`
  is their elementwise absolute difference.
- **Between-group difference.** Take one further independent random
  permutation of `x1` and one of `x2`, each truncated to `m`. Then
  `btw = perm(x2) − perm(x1)`.
  - These draws are **independent** of the within-group draws. Reusing the same
    elements in `btw` and `win` biases the effect: on soil, the gap was 0.18
    against noise of 0.03.
- `winmax = max.(win_1, win_2)`.
- `ratio = btw ./ winmax`. Any `0/0` becomes `0`; a nonzero value divided by 0
  stays `±Inf`.
- The reported values are:
  - `effect = median(ratio)`
  - `diff_btw = median(btw)`
  - `diff_win = median(winmax)`
  - `rab_all = median` of every value of feature `j` over all samples and
    instances.

**Deliberate difference from the reference, required by P6.** The permutations
are drawn as **index** permutations once and shared by every feature. Relabelling
the features then permutes every output column exactly. The reference draws them
per feature. The marginal distribution of each feature's output is the same
under both; only the correlation across features differs. Sharing gives exact
permutation-equivariance.

## Determinism

All randomness comes from one `StableRNG(seed)`. The draws happen in this order:

1. the Dirichlet instances;
2. the effect-size index permutations.

## What is checked, and how

- **Exact:** given the same instances, the p-value columns must match the
  reference to floating-point precision. A diagnostic run on the reference's own
  instances reproduced this.
- **Statistical:** with independent Monte-Carlo draws, the oracle requires
  `|native − R| ≤ 2 × |R(seed 1) − R(seed 2)| + floor` per column
  (`test/oracle/runoracle.jl`).
- **Not proved:** calibration; see `PROOF-STATUS.md` and the README.
