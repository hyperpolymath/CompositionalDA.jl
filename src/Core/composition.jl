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
        # Do not sum machine integers: a valid library can overflow Int/UInt.
        any(>(0), @view counts[i, :]) || throw(CompositionError("sample $i has zero total count"))
    end
    return counts
end

validate_counts(::AbstractMatrix) =
    throw(CompositionError("counts must be an integer matrix (raw counts, not normalised values)"))

# Both native methods and their reference adapters share this boundary. In
# particular, converting missing labels to the string "missing" is forbidden.
function validate_two_groups(counts, features, groups)
    validate_counts(counts)
    n, d = size(counts)
    length(features) == d || throw(CompositionError("feature names do not match columns"))
    length(groups) == n || throw(CompositionError("group labels do not match samples"))
    any(ismissing, features) && throw(CompositionError("feature names cannot be missing"))
    any(ismissing, groups) && throw(CompositionError("group labels cannot be missing"))
    f = string.(features)
    g = string.(groups)
    any(isempty ∘ strip, f) && throw(CompositionError("feature names cannot be empty"))
    any(isempty ∘ strip, g) && throw(CompositionError("group labels cannot be empty"))
    allunique(f) || throw(CompositionError("feature names must be unique"))
    original_levels = try
        sort(unique(groups))
    catch
        throw(CompositionError("group labels must be mutually sortable"))
    end
    levels = string.(original_levels)
    allunique(levels) || throw(CompositionError("distinct group labels must have distinct display names"))
    length(levels) == 2 || throw(CompositionError("exactly two groups are required"))
    i1, i2 = (findall(==(level), groups) for level in original_levels)
    min(length(i1), length(i2)) >= 2 || throw(CompositionError("each group needs at least two samples"))
    return (; features = f, groups = g, levels, i1, i2)
end

function validate_mc(mc_samples, seed)
    mc_samples >= 2 || throw(ArgumentError("mc_samples must be at least 2"))
    0 <= seed <= typemax(UInt64) || throw(ArgumentError("seed must fit UInt64"))
    return nothing
end

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
