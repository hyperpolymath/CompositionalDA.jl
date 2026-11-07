# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
@testset "Reference adapter validation precedes the R runner" begin
    x = [10 20 30 40; 11 22 29 41; 30 20 10 42; 31 21 11 43]
    features = string.(1:4); groups = [1,1,2,2]
    calls = Ref(0)
    runner = f -> (calls[] += 1; error("invalid input entered R runner"))
    for method in (aldex2_reference,ancombc_reference)
        for bad in (fill(-1,4,4),zeros(Int,4,4),fill(big(10)^400,4,4))
            @test_throws CompositionError method(bad,features,groups;runner)
        end
        for labels in ([1,1,2,missing],[NaN,NaN,2,2],[1,1,1,2],[1,1,1,1])
            @test_throws CompositionError method(x,features,labels;runner)
        end
        @test_throws CompositionError method(x,["a","a","b","c"],groups;runner)
    end
    @test_throws ArgumentError aldex2_reference(x,features,groups;seed=big(2)^40,runner)
    @test_throws ArgumentError aldex2_reference(x,features,groups;mc_samples=1,runner)
    @test_throws ArgumentError ancombc_reference(x,features,groups;negative_lower_bound=true,runner)
    @test_throws ArgumentError ancombc_reference(x,features,groups;pseudocount=big"1e-1000",runner)
    @test_throws ArgumentError ancombc_reference(x,features,groups;tolerance=big"1e-1000",runner)
    @test calls[] == 0
end
