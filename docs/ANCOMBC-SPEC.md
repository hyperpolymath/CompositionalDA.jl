<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->

# ANCOM-BC: the implemented model and its limits

`ancombc` implements the **original two-group model** of Lin and Peddada,
*Analysis of compositions of microbiomes with bias correction*, Nature
Communications 11, 3514 (2020), doi:10.1038/s41467-020-17041-7. It was derived
from the published equations. The R estimator source was subsequently inspected
to diagnose reference discrepancies and document its initialization, stopping
and variance conventions; this is not a clean-room conformance exercise.
It is not ANCOM-BC2. Covariates, interactions, repeated measures, global tests
and pairwise testing of more than two groups are outside this API.

## Data and estimates

Rows are independent samples; columns are taxa. Counts must be nonnegative
integers, sample totals positive, identifiers unique and groups complete.
Two samples per group is the computational minimum, not an assurance of
asymptotic calibration. All-zero taxa and taxa below `prevalence` are removed
and named in provenance. At least four tested taxa must remain; this is a
computational refusal rule, not a sufficient statistical sample size.

Both native and reference APIs default to `structural_zeros=false`, matching
the reference method. With `structural_zeros=true`, structural absence is
determined separately in each group from zero
prevalence; `negative_lower_bound=true` additionally uses the published
normal lower bound. Those decisions have explicit statuses and **no numerical
p-value**, effect or q-value. They do not enter BH or the bias fit.

For retained taxa, `y = log(count + pseudocount)`, with explicit default
pseudocount 1. Within each group, the sample offset is its across-taxon mean
minus the group grand mean. Taxon means are group means of `y`. The contrast
is sorted group 2 minus group 1. Residual sum of squares divided by the square
of the group size gives the variance of each group mean (paper equations 6
and 14); the two variances add.

The contrast distribution is fitted as three normal components. The middle
component has known taxon-specific variance `v[i]`; the two alternatives
have variances `v[i] + extra[k]`. Ordered component means are estimated by
weighted isotonic projection. Extra variances are nonnegative. An EM/ECM
iteration evaluates densities in the log domain, verifies nondecreasing
likelihood and requires both likelihood and parameter convergence.
Both stopping checks use `tolerance`; the parameter check is scaled by the
largest absolute parameter. A square-root tolerance would stop too early on
the soil fixture's slowly converging fit.
Initial means are the 0.1/0.5/0.9 quantiles, weights 0.1/0.8/0.1, and
alternative extra variances the observed contrast variance.

The variance-coordinate search examines a fixed logarithmic grid, refines
local brackets by 80 golden-section steps, and retains an update only if it
improves the conditional objective. This finite numerical search does **not**
certify a global mixture-likelihood maximum. Component collapse, non-finite
variance, decreasing likelihood or iteration exhaustion raises
`CompositionError`; no other method is substituted.

The effect subtracts the fitted middle mean, independently of `conserve`,
matching the pinned reference API's distinction between coefficient and
variance estimation. With `conserve=true`, the standard error is
`sqrt(v[i]) + sqrt(bias_variance)`, using the WLS variance approximation in
equation 31 and the conservative bound in equation 35. The paper's equation
35 also substitutes the WLS point estimate in its numerator; this API retains
the EM coefficient when changing the uncertainty option. The WLS point
estimate remains in provenance for comparison, not as the reported effect.
For this WLS approximation, taxa below the fitted lower-component proportion's
quantile form the lower class; taxa at or above the upper-tail quantile form
the upper class. The middle interval is lower-inclusive and upper-exclusive,
matching the reference's tie conventions. This is recorded as
`bias_classification = "mixture-weight-quantiles"` in provenance.
Two-sided normal Wald p-values are BH-adjusted over tested taxa. Natural-log
effects and all fit settings are recorded in provenance.

## What certainty is available

`BiasModel.agda` proves observational shift invariance, exact bias
cancellation and the expression for residual estimation error over any
commutative ring. It does not prove that the estimated bias is the true bias.
The conditional model can shift abundance and sampling bias in opposite
directions without changing observations. Identification therefore needs
assumptions about non-differential taxa; counts cannot establish them alone.

The finite data, independence, weak taxon-dependence, asymptotic normality and
identification assumptions are not formalised. Neither floating-point
execution, statistical calibration, FDR control nor biological truth is
proved. Structural-zero classification is also a statistical decision.

`test/method_ancombc.jl` checks published moment equations against a hand
calculation, nuisance-shift invariance, planted effects, feature permutation,
structural zeros and refusals. `test/oracle/ancombc_oracle.jl` is a separate
conformance gate against ANCOMBC 2.12.0's original `ancombc`. Its thresholds
were fixed before observing results. **A paper-derived implementation is not
automatically equivalent to that R release; the gate must pass before such a
claim or downstream integration is released.**

The reference wrapper allows 20,000 iterations by default. The pinned R
implementation returns its current parameters at its iteration cap without a
convergence error. The oracle records its internal iteration count and residual
and requires actual convergence. Both original fits converged below 1,000
iterations: increasing the reference budget did not resolve their disagreement.

## Reference defect and the qualified comparison

Unmodified ANCOMBC 2.12.0 fails the soil effect, p-value and q-value thresholds.
Its variance optimizer can stop at zero despite an improving positive direction.
For a single squared residual 1 and known variance 0.99, the exact conditional
optimum is extra variance 0.01. The bounded Nelder-Mead call used by the reference
returns zero when started there. This defect also occurs in the latest upstream
original-ANCOM-BC estimator. It is addressed in
[ANCOMBC PR 344](https://github.com/FrederickHuangLin/ANCOMBC/pull/344).

For fixed residuals r_i, known positive variances v_i and nonnegative weights
w_i, the conditional negative log likelihood, apart from a constant, is
Q(t) = (1/2) sum_i w_i [log(v_i+t) + r_i^2/(v_i+t)]. Its derivative is
(1/2) sum_i w_i [(v_i+t)-r_i^2]/(v_i+t)^2. It is nonnegative once
t >= max_i r_i^2, supplying a finite search interval. In the equal-variance
case the optimum is max(sum_i w_i r_i^2 / sum_i w_i - v, 0).
These are calculus derivations, not additional Agda theorems.

The positive conformance gate loads only the bias estimator and its scalar
variance helper from signed upstream-repair commit
`6d094731ed85e9b0458053ad3510e9bdc238e152` into the otherwise pinned ANCOMBC
2.12.0 reference environment. The native Julia optimizer and the R repair use
different searches (grid/golden-section versus Brent). All 27 checks passed at
native revision `480b10bd405c2766f08408105435b33d5dd05050`, including convergence
checks, with the original thresholds unchanged. The largest soil discrepancies
were 3.57e-6 in effect, 4.21e-5 in SE, 1.10e-4 in p and 1.54e-4 in q.

A separately named pristine-reference control reproduces the analytic boundary
failure and marks exactly those three soil comparisons as expected failures.
An unexpected pass fails that control. Without that explicit control mode,
the unmodified-reference oracle still fails normally. A green combined CI is
**not a claim of equivalence to unmodified ANCOMBC 2.12.0**: qualification is
to the stated estimator repair, pending upstream review. Changing that pin
requires a fresh comparison and an updated receipt.
