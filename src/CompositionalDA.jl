# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

"""
    CompositionalDA

Compositional differential-abundance methods for count tables, implemented
natively in Julia: ALDEx2 (`aldex2`) and original two-group ANCOM-BC (`ancombc`).
See the method specifications and acceptance gates. Input is raw integer counts with
samples as rows and features (taxa) as columns.

The core layer is exported: closure, the centred log-ratio, zero policies,
Benjamini–Hochberg and the result/provenance types. Its algebraic properties are
proved in `proofs/agda` and mirrored by property tests in `test/`. See
`proofs/PROOF-STATUS.md` for what is proved and what is only tested.
"""
module CompositionalDA

using DataFrames
using Distributions: Gamma, Normal, ccdf
using HypothesisTests: UnequalVarianceTTest, pvalue
using LinearAlgebra
using Random
using StableRNGs: StableRNG
using Statistics

include("Core/errors.jl")
include("Core/zeros.jl")
include("Core/composition.jl")
include("Core/bh.jl")
include("Core/result.jl")
include("Methods/ALDEx2.jl")
include("Methods/ANCOMBC.jl")

export CompositionError
export ZeroPolicy, RefuseZeros, Pseudocount, MultiplicativeReplacement, replace_zeros
export validate_counts, closure, clr
export bh_adjust
export Provenance, DAResult
export aldex2, aldex2_reference, dirichlet_clr_instances
export ancombc, ancombc_reference

end # module
