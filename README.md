<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# CompositionalDA.jl

Native-Julia compositional differential abundance for microbiome count tables.
Methods: ALDEx2 and ANCOM-BC. Multinomial and Dirichlet-multinomial regression
live in a separate package,
[`CompositionalCounts.jl`](https://github.com/hyperpolymath/CompositionalCounts.jl).
It needs no Python and no TensorFlow. R is needed only for the optional reference
backend.

> **Status: v0.2.0.** The core layer (closure, the centred log-ratio, zero
> policies, Benjamini–Hochberg, result/provenance types) and two-group ALDEx2
> (`aldex2`) are implemented. ALDEx2 is defined by `docs/ALDEx2-SPEC.md` and
> agrees with Bioconductor ALDEx2 1.42.0 within Monte-Carlo noise on the mock,
> gut and soil datasets. ANCOM-BC lands in v0.3.
> Nothing here claims a method that is not yet in `src/`.

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

`proofs/agda` holds machine-checked proofs (Agda 2.7.0.1, `--safe --without-K`,
no postulates) of the algebraic properties the methods rely on. Each CLR theorem is
mirrored by a Julia property test. Each property test is in turn shown to reject
a planted mutant (`test/mutants.jl`). The Benjamini–Hochberg step-up (P5: q ≥ p,
q ≤ 1, order-preserving, at the original indices) is proved over exact
rationals, assuming the input was sorted correctly; the sort itself is not
verified, and nothing about floating point is proved. Statistical properties such as calibration,
consistency and FDR control are **not** proved. They are to be tested against
the R reference implementations; that oracle is not built yet. The full table is in
[`proofs/PROOF-STATUS.md`](proofs/PROOF-STATUS.md).

## The R bridge is temporary by design

**Planned, not yet in the tree.** The optional reference backend will call
Bioconductor ALDEx2/ANCOMBC through RCall, shipped as a Julia package extension.
It will exist for two purposes: to serve as
the test oracle, and to help researchers move off R. It has a written sunset
criterion; see [`docs/R-SUNSET.md`](docs/R-SUNSET.md).

## Running the checks

```sh
julia --project -e 'using Pkg; Pkg.test()'
proofs/bootstrap.sh && proofs/tests/gate-selftest.sh   # needs Agda 2.7.0.1
```

## Licence

Code: MPL-2.0. Documentation: CC-BY-SA-4.0. Per-file SPDX headers are
authoritative.
