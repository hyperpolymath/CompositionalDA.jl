# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "Core units" begin
    @testset "validate_counts refuses non-compositions" begin
        @test validate_counts([1 2; 3 4]) == [1 2; 3 4]
        @test_throws CompositionError validate_counts([1.0 2.0; 3.0 4.0])
        @test_throws CompositionError validate_counts([1 -2; 3 4])
        @test_throws CompositionError validate_counts([0 0; 3 4])
        @test_throws CompositionError validate_counts(reshape([1, 2], 2, 1))
    end
    @testset "zero policies" begin
        c = [0 5 5; 2 2 6]
        @test_throws CompositionError replace_zeros(RefuseZeros(), c)
        @test replace_zeros(Pseudocount(0.5), c) == Float64.(c) .+ 0.5
        m = replace_zeros(MultiplicativeReplacement(0.65), c)
        @test all(m .> 0)
        @test all(isapprox.(sum(m; dims = 2), 1))
        @test m[1, 2] / m[1, 3] ≈ 1.0          # non-zero ratios preserved
        @test m[2, :] ≈ [0.2, 0.2, 0.6]        # zero-free row is just closed
        @test_throws ArgumentError Pseudocount(0)
        @test_throws ArgumentError MultiplicativeReplacement(1.0)
    end
    @testset "clr known answer" begin
        x = [1.0 2.0 4.0]
        @test clr(x) ≈ [-log(2) 0.0 log(2)]
        @test_throws CompositionError clr([0.0 1.0])
        @test clr(Pseudocount(1), [0 1 3]) ≈ clr([1.0 2.0 4.0])
    end
    @testset "bh_adjust matches R p.adjust(\"BH\")" begin
        # Reference: R 4.5 `p.adjust(c(0.01,0.04,0.03,0.005,0.5,0.04,NA,0.2), "BH")`
        p = [0.01, 0.04, 0.03, 0.005, 0.5, 0.04, NaN, 0.2]
        ref = [0.035, 0.056, 0.056, 0.035, 0.5, 0.056, NaN, 0.23333333333333336]
        q = bh_adjust(p)
        @test isnan(q[7])
        @test q[[1:6; 8]] ≈ ref[[1:6; 8]] atol = 1e-15
        @test all(isnan, bh_adjust([NaN, NaN]))
        @test_throws ArgumentError bh_adjust([0.5, 1.5])
    end
    @testset "DAResult requires its columns" begin
        prov = Provenance(:test, :native; seed = 7, mc_samples = 3)
        @test prov.seed == UInt64(7) && prov.parameters[:mc_samples] == 3
        t = DataFrame(feature = ["a"], effect = [0.1], pvalue = [0.2], qvalue = [0.2])
        @test DAResult(t, prov).table === t
        @test_throws ArgumentError DAResult(select(t, Not(:qvalue)), prov)
    end
end
