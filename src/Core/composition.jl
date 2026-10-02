# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

"""
    validate_counts(counts) -> counts

Check that `counts` is a non-empty samples × features matrix of non-negative
integers in which every sample has a positive total. Return it unchanged, or
throw a `CompositionError`. Only integer counts are accepted, so normalised or
already-transformed input is refused by type.
"""
function validate_counts(counts::AbstractMatrix{<:Integer})
    n, d = size(counts)
    n > 0 || throw(CompositionError("no samples"))
    d > 1 || throw(CompositionError("a composition needs at least 2 features, got $d"))
    any(<(0), counts) && throw(CompositionError("negative counts present"))
    for i in 1:n
        sum(@view counts[i, :]) > 0 || throw(CompositionError("sample $i has zero total count"))
    end
    return counts
end

validate_counts(::AbstractMatrix) =
    throw(CompositionError("counts must be an integer matrix (raw counts, not normalised values)"))

"""
    closure(x) -> Matrix{Float64}

Divide each row of a strictly positive samples × features matrix by its row sum,
so that each row sums to 1.
"""
function closure(x::AbstractMatrix{<:Real})
    all(v -> isfinite(v) && v > 0, x) || throw(CompositionError("closure needs strictly positive finite parts"))
    return x ./ sum(x; dims = 2)
end

"""
    clr(x) -> Matrix{Float64}

Centred log-ratio of each row of a strictly positive samples × features matrix:
`clr(x)ᵢⱼ = log xᵢⱼ - mean_k log xᵢₖ`.

Proved properties (`proofs/PROOF-STATUS.md`, P1):

- each row sums to zero;
- each row is invariant to multiplying that row by a positive constant;
- relabelling features permutes the columns in the same way.
"""
function clr(x::AbstractMatrix{<:Real})
    all(v -> isfinite(v) && v > 0, x) || throw(CompositionError("clr needs strictly positive finite parts; apply a ZeroPolicy first"))
    l = log.(x)
    return l .- mean(l; dims = 2)
end

"""
    clr(policy, counts) -> Matrix{Float64}

Validate raw counts, apply the zero policy, then take the centred log-ratio.
"""
clr(policy::ZeroPolicy, counts::AbstractMatrix{<:Integer}) =
    clr(replace_zeros(policy, validate_counts(counts)))
