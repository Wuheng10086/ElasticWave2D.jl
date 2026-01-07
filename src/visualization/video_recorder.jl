# ==============================================================================
# visualization/video_recorder.jl
#
# Wavefield recording callback (no plotting dependencies)
# Video generation is in the CairoMakie extension
# ==============================================================================

using Printf

# ==============================================================================
# Video Configuration
# ==============================================================================

"""
    VideoConfig

Configuration for video recording.

# Fields
- `fields`: Vector of field symbols to record (:p, :vel, :vx, :vz)
- `output_dir`: Directory for output files
- `prefix`: Filename prefix
- `time_scale`: Playback speed relative to real time
- `skip`: Record every N steps
- `downsample`: Spatial downsampling factor
- `clim`: Color limits (nothing for auto)
- `colormap`: Colormap to use

# Example
```julia
config = VideoConfig(
    fields = [:p, :vx],
    time_scale = 0.01,
    skip = 10
)
```
"""
Base.@kwdef struct VideoConfig
    fields::Vector{Symbol} = [:p]
    output_dir::String = "videos"
    prefix::String = "wave"
    time_scale::Float64 = 0.01
    skip::Int = 10
    downsample::Int = 1
    clim::Union{Nothing, Tuple{Float64, Float64}} = nothing
    colormap::Symbol = :RdBu
end

# ==============================================================================
# Field Recorder (saves data, no plotting)
# ==============================================================================

"""
    FieldRecorder

Records wavefield snapshots to memory for later video generation.
Use with CairoMakie extension to generate videos.

# Example
```julia
recorder = FieldRecorder(nx, nz, config)

# Use as callback
for it in 1:nt
    # ... simulation step ...
    record!(recorder, wavefield, it, dt)
end

# Generate video (requires CairoMakie)
using CairoMakie
generate_video(recorder, "output.mp4")
```
"""
mutable struct FieldRecorder
    config::VideoConfig
    nx::Int
    nz::Int
    frames::Vector{Dict{Symbol, Matrix{Float32}}}
    times::Vector{Float32}
    frame_count::Int
end

function FieldRecorder(nx::Int, nz::Int, config::VideoConfig)
    # Adjust for downsampling
    nx_out = nx ÷ config.downsample
    nz_out = nz ÷ config.downsample
    
    FieldRecorder(config, nx_out, nz_out, [], Float32[], 0)
end

"""
Record a frame if it's time (based on skip setting).
"""
function record!(recorder::FieldRecorder, wavefield::Wavefield, it::Int, dt::Float32)
    if mod(it - 1, recorder.config.skip) != 0
        return
    end
    
    ds = recorder.config.downsample
    frame = Dict{Symbol, Matrix{Float32}}()
    
    for field in recorder.config.fields
        data = if field == :p
            # Pressure = -(txx + tzz) / 2
            p = @. -(wavefield.txx + wavefield.tzz) / 2
            Array(p)
        elseif field == :vel
            # Velocity magnitude
            vel = @. sqrt(wavefield.vx^2 + wavefield.vz^2)
            Array(vel)
        elseif field == :vx
            Array(wavefield.vx)
        elseif field == :vz
            Array(wavefield.vz)
        else
            error("Unknown field: $field")
        end
        
        # Downsample
        if ds > 1
            data = data[1:ds:end, 1:ds:end]
        end
        
        frame[field] = Float32.(data)
    end
    
    push!(recorder.frames, frame)
    push!(recorder.times, (it - 1) * dt)
    recorder.frame_count += 1
end

"""
Get total number of recorded frames.
"""
n_frames(recorder::FieldRecorder) = recorder.frame_count

"""
Clear all recorded frames.
"""
function clear!(recorder::FieldRecorder)
    empty!(recorder.frames)
    empty!(recorder.times)
    recorder.frame_count = 0
end

# ==============================================================================
# Multi-Field Recorder (callback interface)
# ==============================================================================

"""
    MultiFieldRecorder

Callback-compatible recorder for use with run_shots!().

# Example
```julia
config = VideoConfig(fields=[:p], skip=10)
recorder = MultiFieldRecorder(nx, nz, dt, config)

results = run_shots!(backend, wavefield, medium, habc, fd_coeffs,
                     rec, shots, params;
                     on_step = recorder)

# Access recorded data
frames = recorder.recorder.frames
```
"""
mutable struct MultiFieldRecorder
    recorder::FieldRecorder
    dt::Float32
end

function MultiFieldRecorder(nx::Int, nz::Int, dt::Float32, config::VideoConfig)
    MultiFieldRecorder(FieldRecorder(nx, nz, config), dt)
end

# Make it callable as a callback
function (mfr::MultiFieldRecorder)(wavefield::Wavefield, it::Int)
    record!(mfr.recorder, wavefield, it, mfr.dt)
end

# Cleanup (no-op for memory recorder)
function Base.close(mfr::MultiFieldRecorder)
    @info "Recording complete" frames=mfr.recorder.frame_count
end

# ==============================================================================
# Video Generation
# ==============================================================================

"""
    generate_video(recorder::FieldRecorder, output::String; fps=30, colormap=:RdBu)

Generate MP4 video from recorded frames.

# Arguments
- `recorder`: FieldRecorder with recorded frames
- `output`: Output filename (e.g., "wavefield.mp4")
- `fps`: Frames per second (default: 30)
- `colormap`: Colormap to use (default: :RdBu)

# Example
```julia
recorder = FieldRecorder(nx, nz, VideoConfig(fields=[:p], skip=5))
# ... run simulation with recorder callback ...
generate_video(recorder, "output.mp4"; fps=30)
```
"""
function generate_video(recorder::FieldRecorder, output::String; 
                        fps::Int=30, colormap::Symbol=:RdBu)
    
    frames = recorder.frames
    if isempty(frames)
        @warn "No frames recorded!"
        return nothing
    end
    
    # Get dimensions from first frame
    first_frame = first(values(frames))[1]
    nz_ds, nx_ds = size(first_frame)
    
    # Create figure
    fig = Figure(size=(800, 600))
    
    # Determine color limits
    all_data = vcat([vec(f) for f in first(values(frames))]...)
    clim_val = maximum(abs.(all_data)) * 0.8
    if clim_val == 0
        clim_val = 1.0
    end
    
    # Create observable for animation
    field_name = first(keys(frames))
    data = Observable(frames[field_name][1])
    time_text = Observable("t = 0.000 s")
    
    ax = Axis(fig[1, 1], 
              title=time_text,
              xlabel="X grid", 
              ylabel="Z grid",
              aspect=DataAspect())
    
    hm = heatmap!(ax, data, 
                  colormap=colormap, 
                  colorrange=(-clim_val, clim_val))
    ax.yreversed = true
    
    Colorbar(fig[1, 2], hm, label=string(field_name))
    
    # Record video
    n_frames = length(frames[field_name])
    @info "Generating video" frames=n_frames output=output
    
    record(fig, output, 1:n_frames; framerate=fps) do i
        data[] = frames[field_name][i]
        time_text[] = @sprintf("t = %.4f s", recorder.times[i])
    end
    
    @info "Video saved" path=output
    return output
end

# Legacy compatibility - VideoRecorder is now an alias
const VideoRecorder = FieldRecorder
