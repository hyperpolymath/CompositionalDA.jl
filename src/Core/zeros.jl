# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

"""
    ZeroPolicy

How zeros are handled before a log-ratio transform. Logs of zero are undefined,
so every log-ratio path must choose one of these policies explicitly:

- [`RefuseZeros`](@ref): throw a `CompositionError` if any zero is present.
- [`Pseudocount`](@ref): add a constant to every cell.
- [`MultiplicativeReplacement`](@ref): replace zeros by `δ` and shrink the
  non-zero parts so that each row keeps its total (Martín-Fernández et al. 2003).
"""
abstract type ZeroPolicy end

"""
    RefuseZeros()

Zero policy that refuses any zero cell.
"""
struct RefuseZeros <: ZeroPolicy end

"""
    Pseudocount(α = 0.5)

Zero policy that adds `α > 0` to every cell. The 0.5 default matches the
Dirichlet prior that ALDEx2 uses.
"""
struct Pseudocount <: ZeroPolicy
    α::Float64
    function Pseudocount(α::Real = 0.5)
        (isfinite(α) && α > 0) || throw(ArgumentError("pseudocount must be finite and > 0, got $α"))
        stored = Float64(α)
        isfinite(stored) && stored > 0 || throw(ArgumentError("pseudocount must be representable as a positive finite Float64"))
        new(stored)
    end
end

"""
    MultiplicativeReplacement(δ = 0.65)

Zero policy for multiplicative replacement on closed data. Within each row,
zeros become `δ × (1 / row total)` of the closed composition, and the non-zero
parts are scaled by `1 - (sum of replacements)`. Ratios among the non-zero parts
are therefore preserved. `δ` is a fraction of one count, in (0, 1).
"""
struct MultiplicativeReplacement <: ZeroPolicy
    δ::Float64
    function MultiplicativeReplacement(δ::Real = 0.65)
        (isfinite(δ) && 0 < δ < 1) || throw(ArgumentError("δ must be in (0, 1), got $δ"))
        stored = Float64(δ)
        0 < stored < 1 || throw(ArgumentError("δ must be representable as a Float64 strictly between 0 and 1"))
        new(stored)
    end
end

"""
    replace_zeros(policy, counts) -> Matrix{Float64}

Validate a samples × features count matrix, apply a zero policy and return a
strictly positive matrix. For `MultiplicativeReplacement` the result is closed:
each row sums to 1. For the other policies, rows keep their count scale.
"""
function replace_zeros end

function replace_zeros(::RefuseZeros, counts::AbstractMatrix{<:Integer})
    out = float_count_matrix(counts)
    any(iszero, out) && throw(CompositionError("zeros present and the zero policy is RefuseZeros"))
    return out
end

function replace_zeros(p::Pseudocount, counts::AbstractMatrix{<:Integer})
    out = float_count_matrix(counts) .+ p.α
    all(v -> isfinite(v) && v > 0, out) || throw(CompositionError("pseudocount addition is not representable as positive finite Float64 values"))
    return out
end

function replace_zeros(p::MultiplicativeReplacement, counts::AbstractMatrix{<:Integer})
    validate_counts(counts)
    out = Matrix{Float64}(undef, size(counts))
    for i in axes(counts, 1)
        row = @view counts[i, :]
        scale = maximum(row)
        scaled = Float64.(row ./ scale)
        scaled_total = sum(scaled)
        nz = count(iszero, row)
        r = nz == 0 ? 0.0 : Float64(p.δ / scale / scaled_total)
        nz == 0 || r > 0 || throw(CompositionError("sample $i: zero replacement underflows Float64"))
        shrink = 1 - nz * r
        shrink > 0 || throw(CompositionError("sample $i: replacement δ too large for $nz zeros"))
        for j in axes(counts, 2)
            out[i, j] = iszero(row[j]) ? r : shrink * scaled[j] / scaled_total
        end
    end
    all(v -> isfinite(v) && v > 0, out) || throw(CompositionError("replacement proportions are not representable as positive finite Float64 values"))
    return out
end
