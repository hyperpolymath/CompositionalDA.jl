# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "ANCOM-BC original two-group model" begin
    @testset "equations 6 and 14, known answer and nuisance offsets" begin
        y = [1.0 3 5 7; 3 1 7 5; 2 5 4 9; 4 3 8 5]
        moments = CompositionalDA.ancombc_moments(y, [1, 2], [3, 4])
        @test moments.difference == [1, 2, 0, 1]
        @test moments.variance == [1, 1, 2.5, 2.5]
        @test moments.offsets == zeros(4)
        shifted = CompositionalDA.ancombc_moments(y .+ [10, -4, 12, -2], [1, 2], [3, 4])
        @test shifted.difference == moments.difference .+ 2
        @test shifted.variance == moments.variance
        @test shifted.offsets == [7, -7, 7, -7]
    end
    rng = StableRNG(2026)
    n, d = 20, 40
    base = collect(range(1000, 5000; length = d))
    x = round.(Int, permutedims(base) .* exp.(0.25 .* randn(rng, n, d)))
    x[11:20, 1:3] .*= 8
    x[11:20, 4:5] .= max.(1, x[11:20, 4:5] .÷ 4)
    features = ["f$j" for j in 1:d]
    groups = vcat(fill("a", 10), fill("b", 10))
    result = ancombc(x, features, groups)
    @test all(result.table.qvalue[1:5] .< 0.05)
    @test all(result.table.effect[1:3] .> 0)
    @test all(result.table.effect[4:5] .< 0)
    @test result.provenance.parameters[:iterations] < 1000
    @test result.table == ancombc(x, features, groups).table
    ordinary = ancombc(x, features, groups; conserve=false)
    @test ordinary.table.effect == result.table.effect
    @test all(result.table.se .>= ordinary.table.se)

    perm = reverse(1:d)
    permuted = ancombc(x[:, perm], features[perm], groups)
    @test permuted.table.effect ≈ result.table.effect[perm] atol = 1e-8
    @test permuted.table.qvalue ≈ result.table.qvalue[perm] atol = 1e-8

    # Structural absence is a separate decision, not a fabricated p-value.
    absent = vcat(zeros(Int, 10), fill(200, 10))
    withzero = ancombc(hcat(x, absent, zeros(Int, n)), vcat(features, ["absent_a", "absent_all"]), groups; structural_zeros=true)
    row = withzero.table[end, :]
    @test row.status == "structural_zero"
    @test row.structural_zero_group1 && !row.structural_zero_group2
    @test isnan(row.pvalue) && isnan(row.qvalue) && isnan(row.effect)
    @test withzero.provenance.parameters[:excluded_features] == ["absent_all"]
    @test withzero.table.effect[1:d] == result.table.effect

    @test_throws CompositionError ancombc(x, features, groups; max_iterations = 1)
    @test_throws ArgumentError ancombc(x, features, groups; pseudocount = NaN)
    @test_throws ArgumentError ancombc(x, features, groups; tolerance = Inf)
    @test_throws CompositionError ancombc(fill(1, 4, 4), ["a", "b", "c", "d"], [1, 1, 2, 2])
    @test_throws CompositionError ancombc(x[:, 1:3], features[1:3], groups)
end
