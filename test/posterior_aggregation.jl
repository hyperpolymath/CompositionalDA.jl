# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

@testset "Consistent-sign aggregation (exact rational mirror)" begin
    # Opposite posterior directions are the load-bearing counterexample:
    # averaging each instance's two-sided p-value is not ALDEx2's rule.
    g = [1//10, 1//1]
    l = [1//1, 1//10]
    @test min(sum(g), sum(l)) / 2 == 11//20
    @test sum(min.(g, l)) / 2 == 1//10
    @test sum(min.(g, l)) < min(sum(g), sum(l))
    rng = StableRNG(98)
    for _ in 1:100
        left = rand(rng, -20:20, 12) .// 20
        right = rand(rng, -20:20, 12) .// 20
        @test sum(min.(left, right)) <= min(sum(left), sum(right))
    end
end
