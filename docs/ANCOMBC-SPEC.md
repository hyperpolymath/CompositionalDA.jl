<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->

# ANCOM-BC: the implemented model and its limits

`ancombc` implements the **original two-group model** of Lin and Peddada,
*Analysis of compositions of microbiomes with bias correction*, Nature
Communications 11, 3514 (2020), doi:10.1038/s41467-020-17041-7. It was derived
from the published equations, not copied or translated from the GPL R source.
It is not ANCOM-BC2. Covariates, interactions, repeated measures, global tests
and pairwise testing of more than two groups are outside this API.

## Data and estimates

Rows are independent samples; columns are taxa. Counts must be nonnegative
integers, sample totals positive, identifiers unique and groups complete.
Two samples per group is the computational minimum, not an assurance of
asymptotic calibration. All-zero taxa and taxa below `prevalence` are removed
and named in provenance. At least four tested taxa must remain; this is a
computational refusal rule, not a sufficient statistical sample size.

Structural absence is determined separately in each group from zero
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
Initial means are the 0.1/0.5/0.9 quantiles, weights 0.1/0.8/0.1, and
alternative extra variances the observed contrast variance.

The variance-coordinate search examines a fixed logarithmic grid, refines
local brackets by 80 golden-section steps, and retains an update only if it
improves the conditional objective. This finite numerical search does **not**
certify a global mixture-likelihood maximum. Component collapse, non-finite
variance, decreasing likelihood or iteration exhaustion raises
`CompositionError`; no other method is substituted.

The ordinary effect subtracts the fitted middle mean. With `conserve=true`,
the effect uses the hard-classification weighted bias estimate (equation 22)
and the standard error is `sqrt(v[i]) + sqrt(bias_variance)` (equation 35).
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
