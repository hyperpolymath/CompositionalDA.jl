# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "ALDEx2 data boundary and units" begin
    x = [2 4 8; 3 6 9; 5 8 3; 6 9 4]
    f = ["x", "y", "z"]
    g = ["a", "a", "b", "b"]
    run(x, f, g) = aldex2(x, f, g; mc_samples = 8)
    @test_throws CompositionError run(x, ["x", "x", "z"], g)
    @test_throws CompositionError run(x, ["x", missing, "z"], g)
    @test_throws CompositionError run(x, f, ["a", missing, "b", "b"])
    @test_throws CompositionError run(x, f, ["a", "a", "b"])
    @test_throws CompositionError run(x, f, ["a", "a", " ", " "])
    @test_throws ArgumentError aldex2(x, f, g; mc_samples = 1)
    @test_throws ArgumentError aldex2(x, f, g; seed = -1)
    @test_throws ArgumentError dirichlet_clr_instances(StableRNG(1), x; prior = Inf)
    @test_throws CompositionError dirichlet_clr_instances(StableRNG(1), -x)
    @test_throws CompositionError dirichlet_clr_instances(StableRNG(1), zeros(Int, 4, 3))
    @test validate_counts(fill(typemax(Int), 2, 2)) == fill(typemax(Int), 2, 2)

    baseline = run(x, f, g)
    withzero = run(hcat(zeros(Int, 4), x), vcat(["absent"], f), g)
    @test withzero.table == baseline.table
    @test withzero.provenance.parameters[:excluded_all_zero] == ["absent"]
    @test_throws CompositionError run(hcat(x[:, 1], zeros(Int, 4, 2)), f, g)

    # Independent reconstruction of the first draw catches natural-log output:
    # effect and test statistics alone cannot detect the wrong log base.
    r = StableRNG(12)
    gamma = [rand(r, CompositionalDA.Gamma(x[1, j] + 0.5, 1.0)) for j in 1:3]
    expected = log2.(gamma) .- mean(log2.(gamma))
    inst = dirichlet_clr_instances(StableRNG(12), x; mc_samples = 2)
    @test inst[1, :, 1] ≈ expected atol = 1e-14
    @test baseline.provenance.parameters[:log_base] == 2
end
