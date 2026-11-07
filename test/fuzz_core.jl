# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

# A policy that could hide invalid counts if the CLR boundary skipped validation.
struct AbsoluteCountsPolicy <: ZeroPolicy end
CompositionalDA.replace_zeros(::AbsoluteCountsPolicy, counts::AbstractMatrix{<:Integer}) =
    Float64.(abs.(counts)) .+ 1

@testset "Adversarial core boundaries" begin
    @test closure(fill(typemax(Int), 1, 2)) ≈ fill(0.5, 1, 2)
    @test closure(fill(floatmax(Float64), 1, 2)) ≈ fill(0.5, 1, 2)
    @test closure(fill(big(10)^1000, 1, 2)) ≈ fill(0.5, 1, 2)
    for T in (Float16, Float32)
        narrow = [nextfloat(zero(T)) floatmax(T)]
        expected = Float64.(BigFloat.(narrow) ./ sum(BigFloat.(narrow)))
        @test all(>(0), closure(narrow))
        @test closure(narrow) ≈ expected rtol=2e-15
        @test closure(Matrix{Real}(narrow)) ≈ expected rtol=2e-15
    end
    @test_throws CompositionError clr(AbsoluteCountsPolicy(), [-1 2])
    @test clr(AbsoluteCountsPolicy(), [1 2]) ≈ clr([2.0 3.0])
    near_one = replace_zeros(MultiplicativeReplacement(prevfloat(1.0)), [3 2 0 0 0 0 0])
    @test all(>(0), near_one)
    @test sum(near_one) ≈ 1.0
    @test near_one[1, 1] / near_one[1, 2] ≈ 3 / 2
    @test near_one[1, 1] ≈ 3 / 5 * (1 - prevfloat(1.0))
    many_zeros = reshape(vcat([1025, 1025], zeros(Int, 2050)), 1, :)
    @test all(>(0), replace_zeros(MultiplicativeReplacement(prevfloat(1.0)), many_zeros))
    @test_throws CompositionError replace_zeros(MultiplicativeReplacement(prevfloat(1.0)), [1 0 0])
    for empty in (zeros(0, 2), zeros(1, 0))
        @test_throws CompositionError closure(empty)
        @test_throws CompositionError clr(empty)
    end
    @test_throws CompositionError closure([BigFloat("1e-1000") BigFloat(1)])
    @test replace_zeros(MultiplicativeReplacement(), fill(typemax(UInt64), 1, 2)) ≈ fill(0.5, 1, 2)
    @test replace_zeros(MultiplicativeReplacement(), fill(big(10)^1000, 1, 2)) ≈ fill(0.5, 1, 2)
    @test_throws CompositionError replace_zeros(MultiplicativeReplacement(), [big(0) big(10)^1000])
    for policy in (RefuseZeros(), Pseudocount(), MultiplicativeReplacement())
        @test_throws CompositionError replace_zeros(policy, [-1 2])
        @test_throws CompositionError replace_zeros(policy, zeros(Int, 1, 2))
    end
    @test_throws CompositionError replace_zeros(Pseudocount(1e308), fill(big(10)^308, 1, 2))
    @test_throws CompositionError replace_zeros(RefuseZeros(), fill(big(10)^1000, 1, 2))
    @test_throws ArgumentError Pseudocount(BigFloat("1e-400"))
    @test_throws ArgumentError Pseudocount(BigFloat("1e400"))
    @test_throws ArgumentError MultiplicativeReplacement(BigFloat(1) - big(2.0)^(-200))
    @test_throws ArgumentError bh_adjust([BigFloat(1) + big(2.0)^(-200)])
    @test_throws ArgumentError bh_adjust([-BigFloat("1e-400")])

    x = [2 4; 3 5; 4 2; 5 3]
    signed = CompositionalDA.validate_two_groups(x, ["a", "b"], [-0.0, -0.0, 0.0, 0.0])
    @test signed.i1 == [1, 2] && signed.i2 == [3, 4]
    mixed = CompositionalDA.validate_two_groups(x, ["a", "b"], Any[2, 2.0, 10, 10.0])
    @test mixed.groups == ["2", "2", "10", "10"]
    for missing_label in (missing, nothing, NaN, Inf, -Inf)
        @test_throws CompositionError CompositionalDA.validate_two_groups(x, ["a", "b"],
            Any[missing_label, missing_label, 1, 1])
    end
    @test_throws CompositionError CompositionalDA.validate_two_groups(x, ["a", nothing], [1, 1, 2, 2])
end

@testset "BH preserves pointwise bounds without multiply-divide cancellation" begin
    @test all(bh_adjust([0.1,0.3,0.7]) .>= [0.1,0.3,0.7])
    for seed in (23,20261107,0xcafe)
        rng=StableRNG(seed)
        for _ in 1:256
            probabilities=rand(rng,rand(rng,1:100))
            @test all(probabilities .<= bh_adjust(probabilities) .<= 1)
        end
    end
end

@testset "Seeded count and composition fuzz" begin
    for seed in (23, 20261107, 0xcafe), case in 1:128
        @testset "seed=$seed case=$case" begin
            rng = StableRNG(UInt64(seed) + UInt64(case) * 1009)
            T = (Int8, Int16, Int64, UInt64, BigInt)[mod1(case, 5)]
            n, d = rand(rng, 1:8), rand(rng, 2:12)
            x = T.(rand(rng, 0:100, n, d))
            if T != BigInt && case % 3 == 0
                x[rand(rng, 1:n), rand(rng, 1:d)] = typemax(T)
            end
            for i in 1:n
                any(>(0), x[i, :]) || (x[i, 1] = one(T))
            end
            @test validate_counts(x) === x
            positive = max.(x, one(T))
            exact = BigFloat.(positive)
            reference = Float64.(exact ./ sum(exact; dims=2))
            @test closure(positive) ≈ reference rtol=2e-15
            scale = exp2(rand(rng, -300:300))
            @test closure(Float64.(positive) .* scale) ≈ reference rtol=2e-15
            @test clr(Float64.(positive) .* scale) ≈ clr(Float64.(positive)) atol=2e-12

            policy = MultiplicativeReplacement(rand(rng, (0.1, 0.5, 0.65, 0.9)))
            totals = sum(BigFloat.(x); dims=2)
            replacements = BigFloat(policy.δ) ./ totals
            shrinks = 1 .- sum(iszero.(x); dims=2) .* replacements
            if any(<=(0), shrinks)
                @test_throws CompositionError replace_zeros(policy, x)
            else
                expected = [Float64(iszero(x[i,j]) ? replacements[i] :
                    shrinks[i] * BigFloat(x[i,j]) / totals[i]) for i in 1:n, j in 1:d]
                result = replace_zeros(policy, x)
                @test result ≈ expected rtol=3e-15
                @test all(>(0), result)
                @test all(isapprox.(sum(result; dims=2), 1; atol=3e-15))
            end
            @test replace_zeros(Pseudocount(0.5), x) == Float64.(x) .+ 0.5
            if any(iszero, x)
                @test_throws CompositionError replace_zeros(RefuseZeros(), x)
            else
                @test replace_zeros(RefuseZeros(), x) == Float64.(x)
            end
        end
    end
end

@testset "Seeded BH fuzz against a high-precision threshold definition" begin
    for seed in (23, 20261107, 0xcafe), case in 1:256
        rng = StableRNG(UInt64(seed) + UInt64(case) * 1009)
        p = round.(rand(rng, rand(rng, 0:32)); digits=rand(rng, 1:5))
        if !isempty(p)
            p[rand(rng, eachindex(p))] = rand(rng, (0.0, 1.0, NaN))
        end
        q = bh_adjust(p)
        finite = findall(!isnan, p)
        bigp = BigFloat.(p[finite])
        for (k, i) in enumerate(finite)
            # Independent definition by all candidate rejection thresholds;
            # no sorted cumulative-minimum loop is reused from the implementation.
            expected = minimum(vcat(BigFloat[1],
                [length(bigp) * t / count(<=(t), bigp) for t in bigp if t >= bigp[k]]))
            @test q[i] ≈ Float64(expected) atol=2e-15
        end
        @test all(isnan, q[findall(isnan, p)])
        permutation = randperm(rng, length(p))
        @test isequal(bh_adjust(p[permutation]), q[permutation])
    end
end
