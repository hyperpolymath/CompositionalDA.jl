# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Property predicates mirroring the Agda theorems in proofs/PROOF-STATUS.md.
# Each takes the implementation under test as an argument, so the same predicate
# is run against the real function (must hold) and a planted mutant (must fail).

using StableRNGs, Random

"""Random strictly positive samples × features matrices from a seeded RNG."""
function positive_matrices(seed; n = 50)
    map(1:n) do k
        r = StableRNG(seed + k)
        exp.(3 .* randn(r, rand(r, 1:12), rand(r, 2:20)))
    end
end

"""P1a: every row of `f(x)` sums to zero."""
prop_clr_sums_to_zero(f; seed = 1) =
    all(x -> all(abs.(sum(f(x); dims = 2)) .< 1e-9), positive_matrices(seed))

"""P1b: scaling each row by a positive constant leaves `f` unchanged."""
prop_clr_scale_invariant(f; seed = 2) = all(positive_matrices(seed)) do x
    c = exp.(3 .* randn(StableRNG(seed), size(x, 1)))
    isapprox(f(c .* x), f(x); atol = 1e-9)
end

"""P1c: permuting the feature columns permutes the output columns identically."""
prop_clr_permutation_equivariant(f; seed = 3) = all(positive_matrices(seed)) do x
    σ = randperm(StableRNG(seed + size(x, 2)), size(x, 2))
    isapprox(f(x[:, σ]), f(x)[:, σ]; atol = 1e-12)
end

"""Random p-value vectors, some with ties and NaNs."""
pvectors(seed; n = 200) = (begin
        r = StableRNG(seed + k)
        p = round.(rand(r, rand(r, 1:40)); digits = rand(r, 1:3))
        p[rand(r, 1:length(p))] = rand(r) < 0.2 ? NaN : p[1]
        p
    end for k in 1:n)

"""P5: adjusted values are ≥ raw, ≤ 1, and order-preserving (p_i ≤ p_j ⇒ q_i ≤ q_j)."""
prop_bh_monotone(f; seed = 5) = all(pvectors(seed)) do p
    q = f(p)
    ok = findall(!isnan, p)
    all(isnan, q[findall(isnan, p)]) &&
        all(q[ok] .>= p[ok] .- 1e-15) && all(q[ok] .<= 1) &&
        all(q[i] <= q[j] + 1e-15 for i in ok, j in ok if p[i] <= p[j])
end

"""Random count matrices (samples × features) with two groups of ≥ 3 samples."""
function count_tables(seed; n = 8)
    map(1:n) do k
        r = StableRNG(seed + k)
        s = rand(r, 6:10); d = rand(r, 3:12)
        x = rand(r, 0:200, s, d); x[:, 1] .+= 1
        g = vcat(fill("a", s ÷ 2), fill("b", s - s ÷ 2))
        (x, g)
    end
end

"""P3: every Monte-Carlo instance from generator `f(rng, counts)` is centred
like a CLR: each sample's row of each instance sums to zero."""
prop_instances_are_clr(f; seed = 7) = all(count_tables(seed)) do (x, _)
    inst = f(StableRNG(seed), x)
    all(abs.(sum(inst; dims = 2)) .< 1e-9)
end

"""P6: given fixed instances, relabelling the features permutes every output
column of estimator `f(inst, i1, i2; rng)` identically."""
prop_aldex2_equivariant(f; seed = 8) = all(count_tables(seed)) do (x, g)
    inst = CompositionalDA.dirichlet_clr_instances(StableRNG(seed), x; mc_samples = 16)
    i1 = findall(==("a"), g); i2 = findall(==("b"), g)
    σ = randperm(StableRNG(seed + size(x, 2)), size(x, 2))
    t = f(inst, i1, i2; rng = StableRNG(seed))
    tσ = f(inst[:, σ, :], i1, i2; rng = StableRNG(seed))
    all(isapprox(Matrix(tσ), Matrix(t)[σ, :]; atol = 1e-12))
end
