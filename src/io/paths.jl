struct OutputPaths
    base_dir::String
end

OutputPaths(; base_dir::AbstractString="outputs") = OutputPaths(String(base_dir))

OutputPaths(base_dir::AbstractString) = OutputPaths(String(base_dir))

results_dir(o::OutputPaths) = o.base_dir
videos_dir(o::OutputPaths) = o.base_dir
figures_dir(o::OutputPaths) = o.base_dir
logs_dir(o::OutputPaths) = o.base_dir

function ensure_output_dirs(o::OutputPaths)
    mkpath(o.base_dir)
    return o
end

function resolve_output_path(o::OutputPaths, kind::Symbol, filename::AbstractString)
    kind === :results || kind === :videos || kind === :figures || kind === :logs || error("Unknown output kind: $kind")
    joinpath(o.base_dir, String(filename))
end

default_result_filename(; shot_id=nothing) =
    shot_id === nothing ? "result.jld2" : "shot_$(shot_id).jld2"

default_video_filename(field::Symbol; ext::AbstractString="mp4", shot_id=nothing) =
    shot_id === nothing ? "wavefield_$(field).$(ext)" : "shot_$(shot_id)_wavefield_$(field).$(ext)"
