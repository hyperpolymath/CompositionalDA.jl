# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "Exhaustive small-sample rank distributions" begin
    for n1 in 2:4, n2 in 2:4
        n = n1+n2
        allocations = [[i for i in 1:n if !iszero(mask & (1 << (i-1)))]
            for mask in 0:((1 << n)-1) if count_ones(mask) == n1]
        statistics = [sum(a .> b' ) for a in allocations for b in [setdiff(1:n, a)]]
        tail = CompositionalDA.mann_whitney_upper_tail(n1, n2)
        reported = Float64[]
        exact_p = Rational{Int}[]
        for (a, u) in zip(allocations, statistics)
            b = setdiff(1:n, a)
            greater, less = CompositionalDA.wilcoxon_tails(a, b, tail)
            expected_greater = count(>=(u), statistics) // length(statistics)
            expected_less = count(<=(u), statistics) // length(statistics)
            @test greater ≈ Float64(expected_greater) atol=2e-15
            @test less ≈ Float64(expected_less) atol=2e-15
            reversed = CompositionalDA.wilcoxon_tails(b, a,
                CompositionalDA.mann_whitney_upper_tail(n2, n1))
            @test greater ≈ reversed[2] atol=2e-15
            @test less ≈ reversed[1] atol=2e-15
            push!(reported, min(1, 2greater, 2less))
            push!(exact_p, min(1, 2expected_greater, 2expected_less))
        end
        # Every equiprobable allocation is enumerated: this is an exact finite
        # null check, not a Monte Carlo estimate of type-I error.
        for alpha in unique(exact_p)
            @test count(p -> p <= Float64(alpha)+1e-14, reported) / length(reported) <= Float64(alpha)+1e-14
        end
    end
    tails = CompositionalDA.wilcoxon_tails([1., 2], [3., 4],
        CompositionalDA.mann_whitney_upper_tail(2, 2))
    @test min(1, 2tails[1], 2tails[2]) ≈ 1/3
end

@testset "Directional tails, ties and numerical refusals" begin
    greater, less = CompositionalDA.welch_tails([1e8, 1e8+1, 1e8+2, 1e8+3], [1., 2, 3, 4])
    @test 0 < greater < 1e-20
    opposite = CompositionalDA.welch_tails([1., 2, 3, 4], [1e8, 1e8+1, 1e8+2, 1e8+3])
    @test opposite[2] == greater
    @test less == opposite[1]
    @test_throws CompositionError CompositionalDA.welch_tails(ones(3), ones(3))
    @test_throws CompositionError CompositionalDA.wilcoxon_tails(ones(3), ones(3), nothing)
    tied = CompositionalDA.wilcoxon_tails([1., 1, 2], [2., 3, 3],
        CompositionalDA.mann_whitney_upper_tail(3, 3))
    # Midrank U = 1/2 and tie-adjusted variance = 24/5, without continuity correction.
    @test tied[2] ≈ CompositionalDA.ccdf(CompositionalDA.Normal(), 4 / sqrt(24/5))
    @test sum(tied) ≈ 1
    x = [2 4; 3 5; 4 2; 5 3]
    @test_throws ArgumentError dirichlet_clr_instances(StableRNG(1), x; prior=BigFloat("1e-400"))
    @test_throws ArgumentError dirichlet_clr_instances(StableRNG(1), x; prior=BigFloat("1e400"))
    @test_throws CompositionError aldex2(fill(big(10)^300, 4, 2), ["a", "b"], [1, 1, 2, 2]; mc_samples=2)
    @test_throws ArgumentError ancombc(x, ["a", "b"], [1, 1, 2, 2]; pseudocount=BigFloat("1e-400"))
    @test_throws ArgumentError ancombc(x, ["a", "b"], [1, 1, 2, 2]; tolerance=BigFloat("1e-400"))
end

@testset "Seeded ALDEx2 input and instance fuzz" begin
    for seed in (23, 20261107, 0xcafe), case in 1:16
        @testset "seed=$seed case=$case" begin
            rng = StableRNG(UInt64(seed)+UInt64(case)*1009)
            n1, n2, d = rand(rng, 2:6), rand(rng, 2:6), rand(rng, 2:10)
            x = rand(rng, 1:100, n1+n2, d)
            features = ["f$j" for j in 1:d]
            groups = vcat(fill(2, n1), fill(10, n2))
            first = aldex2(x, features, groups; seed=case, mc_samples=8)
            @test isequal(first.table, aldex2(x, features, groups; seed=case, mc_samples=8).table)
            @test first.provenance.parameters[:contrast] == "10 vs 2"
            for (p, q) in ((:pvalue, :qvalue), (:wi_pvalue, :wi_qvalue))
                @test all(isfinite, first.table[!, p])
                @test all(0 .<= first.table[!, p] .<= first.table[!, q] .<= 1)
            end
            inst = dirichlet_clr_instances(StableRNG(case), x; mc_samples=8)
            rows1, rows2 = collect(1:n1), collect(n1+1:n1+n2)
            direct = CompositionalDA.aldex2_from_instances(inst, rows1, rows2; rng=StableRNG(case))
            reverse = CompositionalDA.aldex2_from_instances(inst, rows1, rows2;
                rng=StableRNG(case), first_is_i1=false)
            for column in (:pvalue, :qvalue, :wi_pvalue, :wi_qvalue)
                @test direct[!, column] ≈ reverse[!, column] atol=2e-14
            end
            permutation = randperm(rng, d)
            permuted = CompositionalDA.aldex2_from_instances(inst[:, permutation, :], rows1, rows2;
                rng=StableRNG(case))
            for column in propertynames(direct)
                @test isequal(permuted[!, column], direct[permutation, column])
            end
        end
    end
end
