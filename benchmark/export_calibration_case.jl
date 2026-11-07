# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath)
include("calibration_pilot.jl")
destination=only(ARGS)
mkpath(destination)
for replicate in 1:32
    x,features,groups,truth=simulate(REGIMES[3],20261107+30000+replicate)
    open(joinpath(destination,"replicate-$replicate.tsv"),"w") do io
        println(io,join(vcat("feature",["s$i" for i in axes(x,1)]),'\t'))
        for j in axes(x,2)
            println(io,join(vcat(features[j],string.(x[:,j])),'\t'))
        end
    end
    try
        fitted=ancombc(x,features,groups;prevalence=.1,structural_zeros=false,conserve=true)
        open(joinpath(destination,"native-$replicate.csv"),"w") do io
            println(io,"feature,pvalue,qvalue,se,effect")
            for row in eachrow(fitted.table)
                println(io,join((row.feature,row.pvalue,row.qvalue,row.se,row.effect),','))
            end
        end
    catch error
        error isa CompositionError || rethrow()
        write(joinpath(destination,"native-refused-$replicate.txt"),sprint(showerror,error))
    end
end
