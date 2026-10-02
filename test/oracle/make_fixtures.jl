# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Rebuilds test/oracle/fixtures/{gut,soil}_*.tsv from the Nearing et al. (2022)
# figshare deposit (CC BY 4.0). The fixtures are committed; this script exists so
# anyone can check that they are what the deposit and the filter below produce.
#
# Run: julia --project=test/oracle test/oracle/make_fixtures.jl [path/to/tarball]
# With no argument the 253 MB tarball is downloaded to a temporary directory.

using DelimitedFiles, Downloads, SHA

const URL = "https://ndownloader.figshare.com/files/27857901"
const SHA256 = "53337642e26735fc138ed17089f6a467135719413bafbd1f9e5bff048c164b84"
const STUDIES = (gut = ("par_scheperjans", "par_scheperjans_metadata.tsv"),
                 soil = ("Blueberry", "Blueberry_metadata.tsv"))
const MIN_PREVALENCE = 0.25
const MAX_FEATURES = 300
const OUT = joinpath(@__DIR__, "fixtures")

"""
    fetch_tarball(arg) -> path

Return the path to the deposit tarball, downloading it when `arg` is `nothing`,
and refuse it unless its SHA-256 matches the recorded digest.
"""
function fetch_tarball(arg)
    path = arg === nothing ? Downloads.download(URL, joinpath(mktempdir(), "da.tar.gz")) : arg
    got = bytes2hex(open(sha256, path))
    got == SHA256 || error("checksum mismatch for $path: $got")
    return path
end

"""
    read_study(dir, study, metafile) -> (counts, features, samples, groups)

Read one study's raw ASV table (features × samples) and its two-column
metadata, keeping only samples present in both, in table order.
"""
function read_study(dir, study, metafile)
    t = readdlm(joinpath(dir, study, "$(study)_ASVs_table.tsv"), '\t', String; comments = false)
    m = readdlm(joinpath(dir, study, metafile), '\t', String; comments = false)
    group = Dict(m[i, 1] => m[i, 2] for i in 2:size(m, 1))
    keep = [j for j in 2:size(t, 2) if haskey(group, t[1, j])]
    samples = t[1, keep]
    features = t[2:end, 1]
    raw = parse.(Float64, t[2:end, keep])   # Blueberry writes counts as "474.0"
    all(isinteger, raw) || error("$study: non-integer counts; not a raw count table")
    counts = Int.(raw)
    return counts, features, samples, [group[s] for s in samples]
end

"""
    filter_features(counts, features) -> (counts, features)

Keep features present in at least `MIN_PREVALENCE` of samples, then at most the
`MAX_FEATURES` most prevalent (ties by total count, then source order). Source
row order is preserved in the output.
"""
function filter_features(counts, features)
    prev = vec(sum(counts .> 0; dims = 2)) ./ size(counts, 2)
    tot = vec(sum(counts; dims = 2))
    cand = findall(>=(MIN_PREVALENCE), prev)
    ranked = sort(cand; by = i -> (-prev[i], -tot[i], i))
    rows = sort(ranked[1:min(MAX_FEATURES, length(ranked))])
    return counts[rows, :], features[rows]
end

"""
    write_fixture(name, counts, features, samples, groups)

Write `<name>_counts.tsv` (features × samples, header row of sample ids) and
`<name>_groups.tsv` (sample id, group) under `OUT`.
"""
function write_fixture(name, counts, features, samples, groups)
    all(>(0), sum(counts; dims = 1)) || error("$name: a sample has no reads after filtering")
    mkpath(OUT)
    open(joinpath(OUT, "$(name)_counts.tsv"), "w") do io
        println(io, join(vcat("feature", samples), '\t'))
        for i in eachindex(features)
            println(io, join(vcat(features[i], string.(counts[i, :])), '\t'))
        end
    end
    open(joinpath(OUT, "$(name)_groups.tsv"), "w") do io
        println(io, "sample\tgroup")
        foreach(k -> println(io, samples[k], '\t', groups[k]), eachindex(samples))
    end
    println(name, ": ", length(features), " features × ", length(samples), " samples; groups ",
            Dict(g => count(==(g), groups) for g in unique(groups)))
end

"""Extract the needed studies from the tarball and write both fixtures."""
function main(args)
    tarball = fetch_tarball(isempty(args) ? nothing : args[1])
    dir = mktempdir()
    members = ["Hackathon/Studies/$(s[1])" for s in values(STUDIES)]
    run(`tar -xzf $tarball -C $dir $members`)
    for (name, (study, metafile)) in pairs(STUDIES)
        c, f, s, g = read_study(joinpath(dir, "Hackathon", "Studies"), study, metafile)
        c, f = filter_features(c, f)
        write_fixture(String(name), c, f, s, g)
    end
end

main(ARGS)
