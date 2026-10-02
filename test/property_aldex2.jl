# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

using StableRNGs

@testset "ALDEx2" begin
    @testset "P3 every Monte-Carlo instance is a CLR" begin
        gen(rng, x) = dirichlet_clr_instances(rng, x; mc_samples = 8)
        @test prop_instances_are_clr(gen)
        @test !prop_instances_are_clr(mutant_instances_uncentred)
    end
    @testset "P6 equivariant under feature relabelling (given instances)" begin
        @test prop_aldex2_equivariant(CompositionalDA.aldex2_from_instances)
        @test !prop_aldex2_equivariant(mutant_aldex2_per_feature_pairs)
    end
    @testset "deterministic for a fixed seed" begin
        x, g = count_tables(11)[1]
        f = ["f$j" for j in 1:size(x, 2)]
        @test aldex2(x, f, g; mc_samples = 16, seed = 3).table ==
              aldex2(x, f, g; mc_samples = 16, seed = 3).table
    end
    @testset "planted signal is found; null features are not" begin
        r = StableRNG(21)
        base = [400, 300, 200, 100, 80, 60, 40, 30, 20, 10]
        rows = [rand(r, 0.8:0.01:1.2, 10) .* base for _ in 1:16]
        for i in 9:16; rows[i][1] *= 8; end           # feature 1 up 8× in group b
        x = round.(Int, permutedims(reduce(hcat, rows)) .* 5)
        g = vcat(fill("a", 8), fill("b", 8))
        t = aldex2(x, ["f$j" for j in 1:10], g; mc_samples = 64).table
        @test t.qvalue[1] < 0.05
        @test t.effect[1] > 1
        # The other features shift only through the closure (by the CLR of the
        # change), so their effects must be small and negative, not significant.
        @test all(abs.(t.effect[2:end]) .< abs(t.effect[1]))
    end
    @testset "input validation" begin
        x = [1 2; 3 4; 5 6; 7 8]
        @test_throws CompositionError aldex2(x, ["a"], ["u", "u", "v", "v"])
        @test_throws CompositionError aldex2(x, ["a", "b"], ["u", "u", "u", "u"])
        @test_throws CompositionError aldex2(x, ["a", "b"], ["u", "v", "v", "v"])
        @test_throws CompositionError aldex2(-x, ["a", "b"], ["u", "u", "v", "v"])
    end
end
