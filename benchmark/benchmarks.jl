# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath)
using BenchmarkTools, CompositionalDA, Statistics, TOML, Dates
include("calibration_pilot.jl")
destination=only(ARGS)
mkpath(destination)
provenance=source_provenance(pkgdir(CompositionalDA);require_clean=true)
save_run_environment(destination)
rows=NamedTuple[]
for (n,d) in ((20,40),(20,400),(100,400))
    regime=(name="benchmark",n1=n÷2,n2=n÷2,depth1=10000,depth2=10000,sd=.35,sparse=false,signals=4)
    x,f,g,_=simulate(regime,20261107+n+d;d)
    for method in (:closure,:clr,:bh,:aldex2,:ancombc)
        p=collect(range(.001,.999;length=d))
        call = method==:closure ? (() -> closure(x)) :
            method==:clr ? (() -> clr(Pseudocount(1),x)) :
            method==:bh ? (() -> bh_adjust(p)) :
            method==:aldex2 ? (() -> aldex2(x,f,g;mc_samples=128,seed=123)) :
            (() -> ancombc(x,f,g;structural_zeros=false))
        try
            call() # compilation/warmup is excluded and labelled as such
            trial=@benchmark $call() samples=10 evals=1 seconds=15
            push!(rows,(method=String(method),samples=n,features=d,status="measured",
                trials=length(trial.times),median_ns=median(trial.times),p95_ns=quantile(trial.times,.95),
                bytes=trial.memory,allocations=trial.allocs))
        catch error
            error isa CompositionError || rethrow()
            push!(rows,(method=String(method),samples=n,features=d,status="refused",
                trials=0,median_ns=NaN,p95_ns=NaN,bytes=0,allocations=0))
            println("REFUSED ",method," ",n,"x",d,": ",sprint(showerror,error))
        end
    end
end
open(joinpath(destination,"benchmarks.csv"),"w") do io
    println(io,join(keys(first(rows)),','))
    foreach(row->println(io,join(values(row),',')),rows)
end
metadata=Dict("julia"=>string(VERSION),"threads"=>Threads.nthreads(),"cpu"=>Sys.CPU_NAME,
    "machine"=>Sys.MACHINE,"word_size"=>Sys.WORD_SIZE,"date_utc"=>string(now(UTC)),
    "benchmarktools"=>string(pkgversion(BenchmarkTools)),"warmup_excluded"=>true,
    "claim"=>"exploratory warm-process baseline; no reference speedup or acceptance budget")
merge!(metadata,provenance)
open(io->TOML.print(io,metadata),joinpath(destination,"benchmarks-metadata.toml"),"w")
