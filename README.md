<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# CompositionalDA.jl

Native-Julia compositional differential abundance for microbiome count tables.
Methods: ALDEx2 (`src/Methods/ALDEx2.jl`) and the original two-group ANCOM-BC
model (`src/Methods/ANCOMBC.jl`). Multinomial and Dirichlet-multinomial regression
live in a separate package,
[`CompositionalCounts.jl`](https://github.com/hyperpolymath/CompositionalCounts.jl).
It needs no Python and no TensorFlow. R is needed only for the optional reference
backend.

> **Status: v0.3 development.** Native implementations and reference adapters
> exist. Release and integration are gated by `docs/ISSUE-23-ACCEPTANCE.md`;
> implementation is not a declaration of reference conformance. The ALDEx2
> gate now checks all statistics and shared Monte-Carlo instances. ANCOM-BC's
> precise model and limitations are described in `docs/ANCOMBC-SPEC.md`.

## Input contract

Raw integer counts, as a samples × features matrix. Normalised or transformed
input is refused by type. Zeros are never handled silently: every log-ratio path
takes an explicit `ZeroPolicy` (`RefuseZeros`, `Pseudocount`,
`MultiplicativeReplacement`).

```julia
using CompositionalDA
counts = [10 0 5; 3 7 2]
z = clr(MultiplicativeReplacement(), counts)   # rows sum to 0
q = bh_adjust([0.01, 0.04, 0.2])                # equals R's p.adjust(, "BH")
```

## Proved and tested are kept apart

`proofs/agda` holds machine-checked proofs (Debian 13's Agda 2.6.4.3 with agda-stdlib 2.1, `--safe --without-K`,
no postulates) of the algebraic properties the methods rely on. Each CLR theorem is
mirrored by a Julia property test. Each property test is in turn shown to reject
a planted mutant (`test/mutants.jl`). The Benjamini–Hochberg step-up (P5: q ≥ p,
q ≤ 1, order-preserving, at the original indices) is proved over exact
rationals, assuming the input was sorted correctly; the sort itself is not
verified, and nothing about floating point is proved. Statistical properties such as calibration,
consistency and FDR control are **not** proved. R reference comparisons are
numerical checks, not proofs of those properties. PosteriorP and BiasModel
add finite aggregation and bias identities. The full table is in
[`proofs/PROOF-STATUS.md`](proofs/PROOF-STATUS.md).

## The R bridge is temporary by design

The optional reference backend calls
Bioconductor ALDEx2/ANCOMBC through RCall, shipped as a Julia package extension.
It exists for two purposes: to serve as
the test oracle, and to help researchers move off R. It has a written sunset
criterion; see [`docs/R-SUNSET.md`](docs/R-SUNSET.md).

## Running the checks

```sh
julia --project -e 'using Pkg; Pkg.test()'
proofs/check.sh                                   # Debian 13: apt install agda agda-stdlib
```

## Licence

Code: MPL-2.0. Documentation: CC-BY-SA-4.0. Per-file SPDX headers are
authoritative.

## Applications in practice

The concrete target is MetaManifold's Differential Abundance tab
([JoshuaJewell/MetaManifold-WebUI](https://github.com/JoshuaJewell/MetaManifold-WebUI),
proposed as [issue #23](https://github.com/JoshuaJewell/MetaManifold-WebUI/issues/23)),
which today fits each taxon with a negative-binomial GLM or a Welch t-test on
CLR values through R. The real-world value of doing the compositional methods
natively is not speed: it is that `proofs/PROOF-STATUS.md` lets a review answer
*"is this safe and sound, or can it produce funny results?"* claim by claim --
each CLR property is either proved under `--safe` or marked tested, so a
pipeline maintainer adopts exactly as much certainty as the method carries,
instead of borrowing all of it from a paper. CI includes separate
Bioconductor comparisons for ALDEx2 and the original two-group ANCOM-BC.
The tests' existence is not evidence that they pass: the release receipt
must identify the source revision, reference versions and actual outcomes.
See `docs/ISSUE-23-ACCEPTANCE.md` for the independent release gates.
