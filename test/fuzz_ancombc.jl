# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "Seeded ANCOM-BC metamorphic and likelihood checks" begin
    for seed in (23, 20261107, 0xcafe)
        rng = StableRNG(seed)
        n, d = 20, 40
        x = round.(Int, permutedims(range(1000,5000; length=d)) .* exp.(0.25 .* randn(rng,n,d)))
        x[11:20,1:3] .*= 8
        x[11:20,4:5] .= max.(1,x[11:20,4:5] .÷ 4)
        f = string.(1:d); g = repeat([1,2]; inner=10)
        result = ancombc(x,f,g)
        @test all(isfinite,result.table.effect)
        @test all(0 .<= result.table.pvalue .<= result.table.qvalue .<= 1)
        @test all(result.table.se .> 0)
        @test result.table == ancombc(x,f,g).table
        perm = randperm(rng,n)
        reordered = ancombc(x[perm,:],f,g[perm])
        reversed = ancombc(x,f,3 .- g)
        # Scaling counts AND the pseudocount adds a constant to every log
        # count. The estimator removes it exactly in real arithmetic.
        scaled = ancombc(7 .* x,f,g; pseudocount=7)
        for col in (:effect,:se,:pvalue,:qvalue)
            @test reordered.table[!,col] ≈ result.table[!,col] atol=1e-7
            @test scaled.table[!,col] ≈ result.table[!,col] atol=1e-7
            @test reversed.table[!,col] ≈ (col == :effect ? -1 : 1) .* result.table[!,col] atol=1e-7
        end
        y = log.(x .+ 1)
        moments = CompositionalDA.ancombc_moments(y,1:10,11:20)
        fit = CompositionalDA.ancombc_bias(moments.difference,moments.variance)
        responsibilities, likelihood = CompositionalDA.mixture_responsibilities(
            moments.difference,moments.variance,fit.centres,fit.extra,fit.mixing)
        @test vec(sum(responsibilities;dims=2)) ≈ ones(d)
        @test likelihood ≈ fit.likelihood
        @test issorted(fit.centres)
        @test sum(fit.mixing) ≈ 1
        @test all(fit.extra .>= 0) && fit.extra[2] == 0
        # Homoscedastic alternative component has an analytic optimum.
        for _ in 1:64
            values = randn(rng,8); variance = rand(rng)+0.01
            weights = rand(rng,8); centre = randn(rng)
            expected = max(sum(weights .* (values .- centre).^2)/sum(weights)-variance,0)
            actual = CompositionalDA.variance_component(values,fill(variance,8),weights,centre)
            @test actual ≈ expected atol=2e-7 rtol=2e-7
        end
        sparse = hcat(x,vcat(ones(Int,1),zeros(Int,19)))
        bounded = ancombc(sparse,vcat(f,"rare"),g; prevalence=0,structural_zeros=true,negative_lower_bound=true)
        @test bounded.table.status[end] == "structural_zero_both"
        @test isnan(bounded.table.pvalue[end]) && isnan(bounded.table.qvalue[end])
        @test_throws ArgumentError ancombc(x,f,g; negative_lower_bound=true)
        @test_throws CompositionError ancombc(fill(big(10)^400,4,4),string.(1:4),[1,1,2,2])
    end
end
