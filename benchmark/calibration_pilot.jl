# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# Run with the package project; results are written only to the supplied directory.
using CompositionalDA, Random, Statistics, Dates, Printf, SHA, TOML
using CompositionalDA: StableRNG
using Distributions: Multinomial

const REGIMES = [
    (name="balanced_null",n1=10,n2=10,depth1=10000,depth2=10000,sd=0.35,sparse=false,signals=0),
    (name="depth_imbalance_null",n1=10,n2=10,depth1=1000,depth2=10000,sd=0.35,sparse=false,signals=0),
    (name="small_unbalanced_null",n1=4,n2=16,depth1=1000,depth2=1000,sd=0.35,sparse=false,signals=0),
    (name="sparse_null",n1=6,n2=6,depth1=200,depth2=200,sd=0.7,sparse=true,signals=0),
    (name="sparse_signal",n1=10,n2=10,depth1=10000,depth2=10000,sd=0.35,sparse=false,signals=4),
    (name="dense_balanced_signal",n1=10,n2=10,depth1=10000,depth2=10000,sd=0.35,sparse=false,signals=20),
]

function simulate(regime,seed;d=40)
    rng=StableRNG(seed)
    baseline=regime.sparse ? exp.(range(-4,0;length=d)) : collect(range(1000,5000;length=d))
    logeffect=zeros(d)
    half=regime.signals÷2
    logeffect[1:half].=log(4)
    logeffect[half+1:2half].=-log(4)
    x=Matrix{Int}(undef,regime.n1+regime.n2,d)
    for i in axes(x,1)
        condition=i>regime.n1
        weights=baseline .* exp.(regime.sd .* randn(rng,d) .+ condition .* logeffect)
        x[i,:]=rand(rng,Multinomial(condition ? regime.depth2 : regime.depth1,weights/sum(weights)))
    end
    return x,string.(1:d),vcat(fill(0,regime.n1),fill(1,regime.n2)),logeffect
end

function wilson(events,total)
    total==0 && return (NaN,NaN)
    z=1.959963984540054; p=events/total; scale=1+z^2/total
    centre=(p+z^2/(2total))/scale
    radius=z*sqrt(p*(1-p)/total+z^2/(4total^2))/scale
    return (max(0,centre-radius),min(1,centre+radius))
end

function run_pilot(destination)
    mkpath(destination)
    protocol=joinpath(@__DIR__,"..","docs","CALIBRATION-PILOT-PROTOCOL.adoc")
    metadata=Dict("protocol_sha256"=>bytes2hex(sha256(read(protocol))),
        "julia"=>string(VERSION),"package_version"=>string(pkgversion(CompositionalDA)),
        "source_head"=>strip(read(`git -C $(@__DIR__) rev-parse HEAD`,String)),
        "started_utc"=>string(now(UTC)),"replicates_per_regime"=>32,
        "mc_samples"=>128,"seed_formula"=>"20261107 + 10000*regime_index + replicate",
        "status"=>"exploratory pilot; not a calibration certificate")
    records=NamedTuple[]
    failures=NamedTuple[]
    for (regime_index,regime) in enumerate(REGIMES), replicate in 1:32
        seed=20261107+10000regime_index+replicate
        x,features,groups,truth=simulate(regime,seed)
        for method in (:aldex2,:ancombc)
            start=time_ns()
            fitted=try
                method==:aldex2 ? aldex2(x,features,groups;mc_samples=128,seed=seed+1_000_000) :
                    ancombc(x,features,groups;prevalence=0.1,structural_zeros=false,conserve=true)
            catch error
                if error isa CompositionError
                    push!(failures,(regime=regime.name,replicate,method=string(method),seed,
                        reason=sprint(showerror,error),seconds=(time_ns()-start)/1e9))
                    continue
                end
                rethrow()
            end
            elapsed=(time_ns()-start)/1e9
            table=fitted.table
            @assert all(0 .<= table.pvalue .<= table.qvalue .<= 1)
            ids=parse.(Int,table.feature)
            for alpha in (0.05,0.10)
                discovered=table.qvalue .<= alpha
                false_calls=count(discovered .& iszero.(truth[ids]))
                true_calls=count(discovered .& .!iszero.(truth[ids]))
                push!(records,(regime=regime.name,replicate,method=string(method),seed,alpha,
                    tested=length(ids),excluded=40-length(ids),false_calls,true_calls,
                    fdp=false_calls/max(false_calls+true_calls,1),
                    power=regime.signals==0 ? NaN : true_calls/regime.signals,seconds=elapsed))
            end
        end
        replicate==32 && println("Completed ",regime.name)
        flush(stdout)
    end
    open(joinpath(destination,"calibration-pilot.csv"),"w") do io
        println(io,join(string.(keys(first(records))),','))
        for record in records
            println(io,join(values(record),','))
        end
    end
    metadata["finished_utc"]=string(now(UTC))
    metadata["successful_fits"]=length(records)÷2
    metadata["refused_fits"]=length(failures)
    metadata["refusals"]=[Dict(string(k)=>v for (k,v) in pairs(f)) for f in failures]
    open(io->TOML.print(io,metadata),joinpath(destination,"calibration-pilot-metadata.toml"),"w")
    open(joinpath(destination,"calibration-pilot-summary.md"),"w") do io
        println(io,"# Exploratory calibration pilot\n\n32 independent datasets per regime. Metrics condition on successful fits; refusals are separate. No calibration certification.\n")
        println(io,"| Regime | Method | q threshold | Fits/32 | Mean FDP | Mean power | Any false discovery: fraction [95% Wilson interval] |")
        println(io,"| --- | --- | --- | --- | --- | --- | --- |")
        for regime in REGIMES,method in ("aldex2","ancombc"),alpha in (0.05,0.10)
            rows=filter(r->r.regime==regime.name && r.method==method && r.alpha==alpha,records)
            isempty(rows) && (println(io,"| $(regime.name) | $method | $alpha | 0/32 | — | — | — |");continue)
            events=count(r->r.false_calls>0,rows)
            lo,hi=wilson(events,length(rows))
            power=regime.signals==0 ? "N/A" : @sprintf("%.3f",mean(r.power for r in rows))
            @printf(io,"| %s | %s | %.2f | %d/32 | %.3f | %s | %.3f [%.3f, %.3f] |\n",
                regime.name,method,alpha,length(rows),mean(r.fdp for r in rows),power,events/length(rows),lo,hi)
        end
        println(io,"\nFor complete-null regimes, any false discovery equals the false-discovery-rate event. For signal regimes that event is FWER, not mean FDP. The interval covers the event probability among successful fits, not arbitrary populations or the mean FDP. See the fixed protocol and refusal records.")
    end
    println("Successful fits: ",metadata["successful_fits"],"; refusals: ",length(failures))
end

if abspath(PROGRAM_FILE)==(@__FILE__)
    run_pilot(isempty(ARGS) ? error("Supply an output directory") : only(ARGS))
end
