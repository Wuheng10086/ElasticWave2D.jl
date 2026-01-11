# ==============================================================================
# simulation/api.jl
#
# High-level API for simplified simulation workflow
# ==============================================================================

# ==============================================================================
# Configuration Structs
# ==============================================================================

"""
    SimulationConfig

Configuration for regular surface simulation.

# Fields
- `nbc`: Absorbing boundary layers (default: 50)
- `fd_order`: Finite difference order (default: 8)
- `dt`: Time step, auto-computed if nothing (default: nothing)
- `nt`: Number of time steps (default: 3000)
- `cfl`: CFL number for auto dt computation (default: 0.4)
- `f0`: Source dominant frequency in Hz (default: 15.0)
- `free_surface`: Enable free surface at top (default: true)
- `output_dir`: Output directory (default: "outputs")
- `save_gather`: Save gather to binary file (default: true)
- `show_progress`: Show progress bar (default: true)
- `plot_gather`: Save gather plot as PNG (default: true)
"""
Base.@kwdef struct SimulationConfig
    nbc::Int = 50
    fd_order::Int = 8
    dt::Union{Float32,Nothing} = nothing
    nt::Int = 3000
    cfl::Float32 = 0.4f0
    f0::Float32 = 15.0f0
    free_surface::Bool = true
    output_dir::String = "outputs"
    save_gather::Bool = true
    show_progress::Bool = true
    plot_gather::Bool = true
end

"""
    IrregularSurfaceConfig

Configuration for irregular surface simulation with IBM.

# Fields
All fields from SimulationConfig, plus:
- `ibm_method`: :direct_zero (stable) or :mirror (accurate) (default: :direct_zero)
- `ibm_iterations`: Number of IBM iterations (default: 3)
- `src_depth`: Source depth below surface in meters (default: 30.0)
- `rec_depth`: Receiver depth below surface in meters (default: 0.0)
"""
Base.@kwdef struct IrregularSurfaceConfig
    nbc::Int = 50
    fd_order::Int = 8
    dt::Union{Float32,Nothing} = nothing
    nt::Int = 3000
    cfl::Float32 = 0.5f0
    f0::Float32 = 15.0f0
    ibm_method::Symbol = :direct_zero
    ibm_iterations::Int = 3
    src_depth::Float32 = 30.0f0
    rec_depth::Float32 = 0.0f0
    output_dir::String = "outputs"
    save_gather::Bool = true
    show_progress::Bool = true
    plot_gather::Bool = true
    plot_model::Bool = true
end

"""
    SimulationResult

Container for simulation results.

# Fields
- `gather`: Recorded seismogram [nt × n_rec]
- `dt`: Time step used
- `nt`: Number of time steps
- `video_file`: Path to generated video (or nothing)
- `gather_file`: Path to saved gather (or nothing)
- `gather_plot`: Path to gather plot (or nothing)
"""
struct SimulationResult
    gather::Matrix{Float32}
    dt::Float32
    nt::Int
    video_file::Union{String,Nothing}
    gather_file::Union{String,Nothing}
    gather_plot::Union{String,Nothing}
end

# ==============================================================================
# Regular Surface API
# ==============================================================================

"""
    simulate!(model, src_x, src_z, rec_x, rec_z; config, video_config=nothing)

Run simulation with regular (flat) free surface.

# Arguments
- `model::VelocityModel`: Velocity model
- `src_x, src_z`: Source position in meters
- `rec_x, rec_z`: Receiver position arrays in meters
- `config`: SimulationConfig
- `video_config`: VideoConfig for wavefield recording (optional)

# Returns
- `SimulationResult`

# Example
```julia
model = VelocityModel(vp, vs, rho, 10.0, 10.0)

# Without video
result = simulate!(model, 1500.0, 50.0, rec_x, rec_z;
    config = SimulationConfig(nt=3000))

# With video
result = simulate!(model, 1500.0, 50.0, rec_x, rec_z;
    config = SimulationConfig(nt=3000),
    video_config = VideoConfig(fields=[:vz], skip=5)
)
```
"""
function simulate!(model::VelocityModel,
    src_x::Real, src_z::Real,
    rec_x::Vector{<:Real}, rec_z::Vector{<:Real};
    config::SimulationConfig=SimulationConfig(),
    video_config::Union{VideoConfig,Nothing}=nothing)

    mkpath(config.output_dir)

    be = is_cuda_available() ? backend(:cuda) : backend(:cpu)
    @info "Simulation started" backend = typeof(be) model = model.name

    vp_max = maximum(model.vp)
    dt = config.dt === nothing ? config.cfl * min(model.dx, model.dz) / vp_max : config.dt
    dt = Float32(dt)

    @info "Parameters" dt_ms = round(dt * 1000, digits=3) nt = config.nt

    medium = init_medium(model, config.nbc, config.fd_order, be; free_surface=config.free_surface)
    habc = init_habc(medium.nx, medium.nz, config.nbc, dt, model.dx, model.dz, vp_max, be)
    fd_coeffs = to_device(get_fd_coefficients(config.fd_order), be)
    wavefield = Wavefield(medium.nx, medium.nz, be)
    params = SimParams(dt, config.nt, model.dx, model.dz, config.fd_order)

    src_i = round(Int, src_x / model.dx) + medium.pad + 1
    src_j = round(Int, src_z / model.dz) + medium.pad + 1
    wavelet = ricker_wavelet(config.f0, dt, config.nt)
    src = Source(src_i, src_j, to_device(wavelet, be))

    n_rec = length(rec_x)
    rec_i = [round(Int, x / model.dx) + medium.pad + 1 for x in rec_x]
    rec_j = [round(Int, z / model.dz) + medium.pad + 1 for z in rec_z]
    rec = Receivers(
        to_device(rec_i, be),
        to_device(rec_j, be),
        to_device(zeros(Float32, config.nt, n_rec), be),
        :vz
    )

    @info "Geometry" source = (src_x, src_z) n_receivers = n_rec

    recorder = nothing
    if video_config !== nothing
        recorder = MultiFieldRecorder(medium.nx, medium.nz, dt, video_config)
    end

    reset!(be, wavefield)

    run_time_loop!(be, wavefield, medium, habc, fd_coeffs, src, rec, params;
        progress=config.show_progress,
        on_step=recorder === nothing ? nothing : (W, info) -> begin
            Fomo.record!(recorder.recorder, W, info.k, dt)
            return true
        end
    )

    @info "Simulation complete"

    gather = be isa CUDABackend ? Array(rec.data) : copy(rec.data)

    video_file, gather_file, gather_plot = _save_outputs(
        gather, dt, config, recorder, video_config
    )

    return SimulationResult(gather, dt, config.nt, video_file, gather_file, gather_plot)
end

# ==============================================================================
# Irregular Surface API
# ==============================================================================

"""
    simulate_irregular!(model, z_surface, src_x, rec_x; config, video_config=nothing)

Run simulation with irregular free surface using IBM.

# Arguments
- `model::VelocityModel`: Velocity model
- `z_surface::Vector{Float32}`: Surface elevation at each x grid point (meters)
- `src_x`: Source x position in meters
- `rec_x`: Receiver x positions in meters
- `config`: IrregularSurfaceConfig
- `video_config`: VideoConfig for wavefield recording (optional)

# Returns
- `SimulationResult`

# Example
```julia
model = VelocityModel(vp, vs, rho, 10.0, 10.0)
z_surface = sinusoidal_surface(nx, dx; amplitude=30, wavelength=1000)

# Without video
result = simulate_irregular!(model, z_surface, 2000.0, rec_x;
    config = IrregularSurfaceConfig(nt=3000))

# With video  
result = simulate_irregular!(model, z_surface, 2000.0, rec_x;
    config = IrregularSurfaceConfig(nt=3000),
    video_config = VideoConfig(fields=[:vz], skip=10)
)
```
"""
function simulate_irregular!(model::VelocityModel,
    z_surface::Vector{Float32},
    src_x::Real,
    rec_x::Vector{<:Real};
    config::IrregularSurfaceConfig=IrregularSurfaceConfig(),
    video_config::Union{VideoConfig,Nothing}=nothing)

    mkpath(config.output_dir)

    be = is_cuda_available() ? backend(:cuda) : backend(:cpu)
    @info "Irregular surface simulation" backend = typeof(be) ibm_method = config.ibm_method

    vp_max = maximum(model.vp)
    dt = config.dt === nothing ? config.cfl * min(model.dx, model.dz) / vp_max : config.dt
    dt = Float32(dt)

    @info "Parameters" dt_ms = round(dt * 1000, digits=3) nt = config.nt

    medium = init_medium(model, config.nbc, config.fd_order, be; free_surface=false)

    surface_cpu = init_irregular_surface(z_surface, medium;
        n_iter=config.ibm_iterations,
        method=config.ibm_method,
        backend=CPU_BACKEND)
    surface = be isa CUDABackend ? to_gpu(surface_cpu) : surface_cpu

    @info "Irregular surface" ghost_points = surface_cpu.n_ghost method = config.ibm_method

    habc = init_habc(medium.nx, medium.nz, config.nbc, dt, model.dx, model.dz, vp_max, be)
    fd_coeffs = to_device(get_fd_coefficients(config.fd_order), be)
    wavefield = Wavefield(medium.nx, medium.nz, be)
    params = SimParams(dt, config.nt, model.dx, model.dz, config.fd_order)

    wavelet = ricker_wavelet(config.f0, dt, config.nt)
    src = setup_irregular_source(Float32(src_x), config.src_depth, wavelet,
        surface_cpu, medium; backend=be)

    n_rec = length(rec_x)
    rec_depths = fill(config.rec_depth, n_rec)
    rec = setup_irregular_receivers(Float32.(rec_x), rec_depths, surface_cpu, medium, config.nt;
        type=:vz, backend=be)

    @info "Geometry" source_x = src_x src_depth = config.src_depth n_receivers = n_rec

    recorder = nothing
    if video_config !== nothing
        recorder = MultiFieldRecorder(medium.nx, medium.nz, dt, video_config)
    end

    if config.plot_model
        _plot_irregular_setup(model, z_surface, src_x, config.src_depth,
            Float32.(rec_x), fill(config.rec_depth, n_rec),
            config.ibm_method,
            joinpath(config.output_dir, "model_setup.png"))
    end

    reset!(be, wavefield)

    if config.show_progress
        prog = Progress(config.nt, desc="Simulating: ")
    end

    for k in 1:config.nt
        time_step_irregular!(be, wavefield, medium, habc, fd_coeffs,
            src, rec, k, params, surface)

        if recorder !== nothing
            Fomo.record!(recorder.recorder, wavefield, k, dt)
        end

        if config.show_progress
            next!(prog)
        end
    end

    synchronize(be)
    @info "Simulation complete"

    gather = be isa CUDABackend ? Array(rec.data) : copy(rec.data)

    # Save surface elevation
    open(joinpath(config.output_dir, "surface_elevation.txt"), "w") do io
        for (i, z) in enumerate(z_surface)
            println(io, "$((i-1) * model.dx) $z")
        end
    end

    video_file, gather_file, gather_plot = _save_outputs(
        gather, dt, config, recorder, video_config
    )

    return SimulationResult(gather, dt, config.nt, video_file, gather_file, gather_plot)
end

# ==============================================================================
# Surface Shape Helper Functions
# ==============================================================================

"""
    flat_surface(nx, dx, depth)

Create a flat surface at constant depth.
"""
flat_surface(nx::Int, dx::Real, depth::Real) = fill(Float32(depth), nx)

"""
    sinusoidal_surface(nx, dx; base_depth=50, amplitude=20, wavelength=1000)

Create a sinusoidal surface.
"""
function sinusoidal_surface(nx::Int, dx::Real;
    base_depth::Real=50.0,
    amplitude::Real=20.0,
    wavelength::Real=1000.0)
    x = Float32.((0:nx-1) .* dx)
    return Float32.(base_depth .+ amplitude .* sin.(2π .* x ./ wavelength))
end

"""
    gaussian_valley(nx, dx; base_depth=50, valley_depth=30, center=nothing, width=200)

Create a surface with a Gaussian valley (depression).
"""
function gaussian_valley(nx::Int, dx::Real;
    base_depth::Real=50.0,
    valley_depth::Real=30.0,
    center::Union{Real,Nothing}=nothing,
    width::Real=200.0)
    x = Float32.((0:nx-1) .* dx)
    center = center === nothing ? nx * dx / 2 : Float32(center)
    return Float32.(base_depth .+ valley_depth .* exp.(-(x .- center) .^ 2 ./ (2 * width^2)))
end

"""
    gaussian_hill(nx, dx; base_depth=80, hill_height=30, center=nothing, width=200)

Create a surface with a Gaussian hill (elevation).
"""
function gaussian_hill(nx::Int, dx::Real;
    base_depth::Real=80.0,
    hill_height::Real=30.0,
    center::Union{Real,Nothing}=nothing,
    width::Real=200.0)
    x = Float32.((0:nx-1) .* dx)
    center = center === nothing ? nx * dx / 2 : Float32(center)
    return Float32.(base_depth .- hill_height .* exp.(-(x .- center) .^ 2 ./ (2 * width^2)))
end

"""
    tilted_surface(nx, dx; depth_left=30, depth_right=70)

Create a linearly tilted surface.
"""
function tilted_surface(nx::Int, dx::Real;
    depth_left::Real=30.0,
    depth_right::Real=70.0)
    return Float32.(range(depth_left, depth_right, length=nx))
end

"""
    step_surface(nx, dx; depth_left=30, depth_right=70, step_position=nothing)

Create a surface with a step (cliff).
"""
function step_surface(nx::Int, dx::Real;
    depth_left::Real=30.0,
    depth_right::Real=70.0,
    step_position::Union{Real,Nothing}=nothing)
    step_idx = step_position === nothing ? nx ÷ 2 : round(Int, step_position / dx)
    z = zeros(Float32, nx)
    z[1:step_idx] .= Float32(depth_left)
    z[step_idx+1:end] .= Float32(depth_right)
    return z
end

"""
    random_surface(nx, dx; base_depth=50, amplitude=10, smoothness=5)

Create a random rough surface.
"""
function random_surface(nx::Int, dx::Real;
    base_depth::Real=50.0,
    amplitude::Real=10.0,
    smoothness::Int=5)
    noise = randn(Float32, nx) .* Float32(amplitude)
    z = zeros(Float32, nx)
    for i in 1:nx
        i_start = max(1, i - smoothness)
        i_end = min(nx, i + smoothness)
        z[i] = sum(noise[i_start:i_end]) / (i_end - i_start + 1)
    end
    return Float32.(base_depth .+ z)
end

"""
    combine_surfaces(surfaces...; method=:add)

Combine multiple surface perturbations. Methods: :add, :min, :max
"""
function combine_surfaces(surfaces...; method::Symbol=:add)
    nx = length(surfaces[1])
    if method == :add
        return Float32.(sum(surfaces))
    elseif method == :min
        result = copy(surfaces[1])
        for s in surfaces[2:end]
            result .= min.(result, s)
        end
        return Float32.(result)
    elseif method == :max
        result = copy(surfaces[1])
        for s in surfaces[2:end]
            result .= max.(result, s)
        end
        return Float32.(result)
    end
end

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

function _save_outputs(gather, dt, config, recorder, video_config)
    video_file = nothing
    gather_file = nothing
    gather_plot = nothing

    n_rec = size(gather, 2)

    if config.save_gather
        gather_file = joinpath(config.output_dir, "gather.bin")
        open(gather_file, "w") do io
            write(io, Int32(config.nt))
            write(io, Int32(n_rec))
            write(io, gather)
        end
        @info "Saved" file = gather_file
    end

    if recorder !== nothing && video_config !== nothing
        # Generate video for each recorded field
        for field in video_config.fields
            video_file = joinpath(config.output_dir, "wavefield_$(field).mp4")
            generate_video(recorder.recorder, video_file;
                fps=video_config.fps, colormap=video_config.colormap)
            @info "Saved" file = video_file
        end
    end

    if config.plot_gather
        gather_plot = joinpath(config.output_dir, "gather.png")
        _plot_gather_simple(gather, dt, gather_plot)
        @info "Saved" file = gather_plot
    end

    return video_file, gather_file, gather_plot
end

function _plot_gather_simple(gather::Matrix{Float32}, dt::Float32, output::String)
    nt, n_rec = size(gather)
    t_axis = (0:nt-1) .* dt

    fig = CairoMakie.Figure(size=(900, 700))
    ax = CairoMakie.Axis(fig[1, 1],
        xlabel="Receiver", ylabel="Time (s)", title="Shot Gather")

    gmax = maximum(abs.(gather))
    data = gmax > 0 ? gather ./ gmax : gather

    hm = CairoMakie.heatmap!(ax, 1:n_rec, t_axis, data',
        colormap=:seismic, colorrange=(-0.3, 0.3))
    ax.yreversed = true
    CairoMakie.Colorbar(fig[1, 2], hm, label="Normalized Amplitude")
    CairoMakie.save(output, fig)
end

function _plot_irregular_setup(model::VelocityModel, z_surface::Vector{Float32},
    src_x::Real, src_depth::Real,
    rec_x::Vector{Float32}, rec_depth::Vector{Float32},
    ibm_method::Symbol, output::String)

    fig = CairoMakie.Figure(size=(1000, 700))

    x_axis = range(0, (model.nx - 1) * model.dx, length=model.nx)
    z_axis = range(0, (model.nz - 1) * model.dz, length=model.nz)

    ax = CairoMakie.Axis(fig[1, 1],
        xlabel="X (m)", ylabel="Z (m)",
        title="Model with Irregular Surface (IBM: $ibm_method)",
        aspect=CairoMakie.DataAspect())

    hm = CairoMakie.heatmap!(ax, x_axis, z_axis, model.vp', colormap=:viridis)

    surf_x = Float32.((0:length(z_surface)-1) .* model.dx)
    CairoMakie.lines!(ax, surf_x, z_surface, color=:white, linewidth=3, label="Free surface")

    src_idx = clamp(round(Int, src_x / model.dx) + 1, 1, length(z_surface))
    src_z = z_surface[src_idx] + src_depth
    CairoMakie.scatter!(ax, [Float32(src_x)], [Float32(src_z)],
        marker=:star5, markersize=20, color=:red, label="Source")

    rec_z = Float32[]
    for i in 1:length(rec_x)
        idx = clamp(round(Int, rec_x[i] / model.dx) + 1, 1, length(z_surface))
        push!(rec_z, z_surface[idx] + rec_depth[i])
    end
    step = max(1, length(rec_x) ÷ 20)
    CairoMakie.scatter!(ax, rec_x[1:step:end], rec_z[1:step:end],
        marker=:dtriangle, markersize=8, color=:cyan, label="Receivers")

    ax.yreversed = true
    CairoMakie.Colorbar(fig[1, 2], hm, label="Vp (m/s)")
    CairoMakie.axislegend(ax, position=:rb)
    CairoMakie.save(output, fig)
end