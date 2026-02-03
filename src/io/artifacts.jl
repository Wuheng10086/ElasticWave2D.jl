using JSON

struct ArtifactsManifest
    paths::Dict{Symbol,String}
end

ArtifactsManifest() = ArtifactsManifest(Dict{Symbol,String}())

function record_artifact!(m::ArtifactsManifest, key::Symbol, path::AbstractString)
    m.paths[key] = String(path)
    return m
end

get_artifact(m::ArtifactsManifest, key::Symbol, default=nothing) = get(m.paths, key, default)

function write_manifest(o::OutputPaths, m::ArtifactsManifest; filename::AbstractString="manifest.json")
    ensure_output_dirs(o)
    path = resolve_output_path(o, :logs, filename)
    open(path, "w") do io
        JSON.print(io, m.paths, 2)
    end
    return path
end
