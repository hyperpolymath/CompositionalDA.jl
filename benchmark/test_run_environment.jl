# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath)
using Test, TOML
include("run_environment.jl")

@testset "Run provenance" begin
    mktempdir() do source
        run(`git -C $source init -q`)
        write(joinpath(source,"source.jl"),"original\n")
        run(`git -C $source add source.jl`)
        run(`git -C $source -c user.name=Test -c user.email=test@example.invalid commit -qm baseline`)
        head=strip(read(`git -C $source rev-parse HEAD`,String))
        @test source_provenance(source;require_clean=true)==Dict("source_head"=>head,"source_tree"=>"clean")
        write(joinpath(source,"source.jl"),"modified\n")
        @test source_provenance(source)["source_head"]==head
        @test_throws ErrorException source_provenance(source;require_clean=true)
        run(`git -C $source checkout -- source.jl`)
        write(joinpath(source,"untracked.jl"),"new source\n")
        @test_throws ErrorException source_provenance(source;require_clean=true)
        withenv("PATH"=>"") do
            provenance=@test_logs (:warn,r"Source revision could not be captured") source_provenance(source;require_clean=true)
            @test provenance==Dict("source_head"=>"unknown","source_tree"=>"unknown")
        end
    end
    mktempdir() do source
        provenance=@test_logs (:warn,r"Source revision could not be captured") source_provenance(source)
        @test provenance["source_head"]=="unknown"
    end
end

@testset "Resolved environment snapshot" begin
    mktempdir() do destination
        save_run_environment(destination)
        environment=Pkg.Types.EnvCache()
        @test read(joinpath(destination,"Project.toml"))==read(environment.project_file)
        @test read(joinpath(destination,"Manifest.toml"))==read(environment.manifest_file)
        manifest=TOML.parsefile(joinpath(destination,"Manifest.toml"))
        @test haskey(only(manifest["deps"]["StatsFuns"]),"version")
    end
    mktempdir() do project
        Pkg.activate(project) do
            write(joinpath(project,"Project.toml"),"[deps]\n")
            @test_throws ErrorException save_run_environment(project)
        end
    end
end
