# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath)
using Pkg

function source_provenance(source; require_clean=false)
    head, state = try
        head = strip(read(`git -C $source rev-parse HEAD`, String))
        status = read(`git -C $source status --porcelain --untracked-files=all`, String)
        (head, isempty(status) ? "clean" : "dirty")
    catch error
        @warn "Source revision could not be captured; Git provenance is unavailable" exception=(error, catch_backtrace())
        ("unknown", "unknown")
    end
    require_clean && state == "dirty" &&
        error("Refusing to benchmark a dirty source tree; commit source changes and keep outputs outside the repository")
    return Dict("source_head" => head, "source_tree" => state)
end

function save_run_environment(destination)
    # EnvCache selects the manifest used by Pkg, including version-specific names.
    environment = Pkg.Types.EnvCache()
    isfile(environment.manifest_file) ||
        error("A resolved manifest is required; instantiate the active project before running")
    cp(environment.project_file, joinpath(destination, "Project.toml"); force=true)
    cp(environment.manifest_file, joinpath(destination, "Manifest.toml"); force=true)
end
