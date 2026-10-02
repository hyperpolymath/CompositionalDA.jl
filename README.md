<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) -->

# CompositionalDA.jl

Native-Julia compositional differential abundance for microbiome count tables.
Methods: ALDEx2, ANCOM-BC, and a Songbird-style multinomial regression.
It needs no Python and no TensorFlow. R is needed only for the optional reference
backend.

> **Status: v0.1.0, core layer only.** Closure, the centred log-ratio, zero
> policies, Benjamini–Hochberg and the result/provenance types are implemented.
> The methods land one per release: ALDEx2 (v0.2), ANCOM-BC (v0.3), multinomial
> (v0.4). Nothing here claims a method that is not yet in `src/`.

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
no postulates) of the algebraic properties the methods rely on. Each theorem is
mirrored by a Julia property test. Each property test is in turn shown to reject
a planted mutant (`test/mutants.jl`). Statistical properties such as calibration,
consistency and FDR control are **not** proved. They are tested against the R
reference implementations. The full table is in
[`proofs/PROOF-STATUS.md`](proofs/PROOF-STATUS.md).

## The R bridge is temporary by design

The optional reference backend calls Bioconductor ALDEx2/ANCOMBC through RCall.
It ships as a Julia package extension and exists for two purposes: to serve as
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
