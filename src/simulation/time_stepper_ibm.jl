# ==============================================================================
# simulation/time_stepper_irregular.jl
#
# Time stepper with irregular free surface support
# ==============================================================================

"""
    time_step_irregular!(backend, W, M, H, a, src, rec, k, params, surface)

Execute single time step (with irregular free surface support)

# Time Step Sequence
1. Backup boundary values (for HABC)
2. Inject source
3. Update velocity + apply HABC
4. Update stress + apply HABC  
5. Apply irregular free surface (IBM) ← new
6. Record receivers
"""
function time_step_irregular!(backend::AbstractBackend, W::Wavefield, M::Medium, H::HABCConfig,
                               a, src::Source, rec::Receivers, k::Int, params::SimParams,
                               surface::Union{IrregularSurface, IrregularSurfaceGPU})
    
    # 1. Backup boundary (for HABC)
    backup_boundary!(backend, W, H, M)
    
    # 2. Source injection
    inject_source!(backend, W, src, k, params.dt)
    
    # 3. Velocity update + HABC
    update_velocity!(backend, W, M, a, params)
    apply_habc_velocity!(backend, W, H, M)
    
    # 4. Stress update + HABC
    update_stress!(backend, W, M, a, params)
    apply_habc_stress!(backend, W, H, M)
    
    # 5. Apply irregular free surface (IBM)
    # This replaces regular apply_free_surface!
    if is_irregular_surface_enabled(surface)
        apply_irregular_free_surface!(backend, W, surface)
    else
        # If irregular surface not enabled, use regular flat free surface
        apply_free_surface!(backend, W, M)
    end
    
    # 6. Record receivers
    record_receivers!(backend, W, rec, k)
    
    return nothing
end

"""
    run_time_loop_irregular!(backend, W, M, H, a, src, rec, params, surface; 
                             progress=true, on_step=nothing)

Run complete time loop (with irregular surface support)
"""
function run_time_loop_irregular!(backend::AbstractBackend, W::Wavefield, M::Medium, H::HABCConfig,
                                   a, src::Source, rec::Receivers, params::SimParams,
                                   surface::Union{IrregularSurface, IrregularSurfaceGPU};
                                   progress::Bool=true, on_step=nothing)
    
    nt = params.nt
    dt = params.dt
    
    if progress
        p = Progress(nt; desc="Simulating (IBM): ", color=:cyan)
    end
    
    for k in 1:nt
        time_step_irregular!(backend, W, M, H, a, src, rec, k, params, surface)
        
        # Callback
        if on_step !== nothing
            info = TimeStepInfo(k, k * dt, nt)
            cont = on_step(W, info)
            if cont === false
                @info "Simulation stopped early at step $k"
                break
            end
        end
        
        if progress
            next!(p)
        end
    end
    
    synchronize(backend)
    return nothing
end

# ==============================================================================
# Helper functions
# ==============================================================================

"""
Create source for this shot(from MultiShotConfig)
"""
function _create_source_for_shot(shot_config::MultiShotConfig, shot_id::Int, M::Medium, backend)
    shot = shot_config.shots[shot_id]
    i = round(Int, shot.source_x / M.dx) + M.pad + 1
    j = round(Int, shot.source_z / M.dz) + M.pad + 1
    wavelet = to_device(shot_config.wavelet, backend)
    return Source(i, j, wavelet)
end

"""
Create receiver data buffer
"""
function _create_receiver_buffer(template::Receivers, nt::Int, backend::CPUBackend)
    data = zeros(Float32, nt, length(template.i))
    return Receivers(copy(template.i), copy(template.j), data, template.type)
end

function _create_receiver_buffer(template::Receivers, nt::Int, backend::CUDABackend)
    i_cpu = template.i isa Array ? template.i : Array(template.i)
    j_cpu = template.j isa Array ? template.j : Array(template.j)
    i_gpu = CuArray(Int32.(i_cpu))
    j_gpu = CuArray(Int32.(j_cpu))
    data = CUDA.zeros(Float32, nt, length(i_gpu))
    return Receivers(i_gpu, j_gpu, data, template.type)
end

# ==============================================================================
# Complete simulation function(irregular surface version)
# ==============================================================================

"""
    run_shot_irregular!(backend, wavefield, medium, habc, fd_coeffs, 
                        rec, src, params, surface; kwargs...)

Run single shot simulation(with irregular surface support)
"""
function run_shot_irregular!(backend::AbstractBackend, W::Wavefield, M::Medium, H::HABCConfig,
                              a, rec_template::Receivers, src::Source, params::SimParams,
                              surface::Union{IrregularSurface, IrregularSurfaceGPU};
                              shot_id::Int=1, show_progress::Bool=true, on_step=nothing)
    
    # Reset wavefield
    reset!(backend, W)
    
    # Create receiver data buffer for this shot
    rec = _create_receiver_buffer(rec_template, params.nt, backend)
    
    # Run simulation
    run_time_loop_irregular!(backend, W, M, H, a, src, rec, params, surface;
                              progress=show_progress, on_step=on_step)
    
    # Collect results
    gather = Array(rec.data)
    
    return ShotResult(
        gather, shot_id,
        src.i, src.j,
        Array(rec.i), Array(rec.j)
    )
end

"""
    run_shots_irregular!(backend, wavefield, medium, habc, fd_coeffs,
                         rec_template, shot_config, params, surface; kwargs...)

Run multi-shot simulation(with irregular surface support)
"""
function run_shots_irregular!(backend::AbstractBackend, W::Wavefield, M::Medium, H::HABCConfig,
                               a, rec_template::Receivers, shot_config::MultiShotConfig,
                               params::SimParams, surface::Union{IrregularSurface, IrregularSurfaceGPU};
                               show_progress::Bool=true,
                               on_step=nothing, on_shot_complete=nothing)
    
    n_shots = length(shot_config.shots)
    results = Vector{ShotResult}(undef, n_shots)
    
    @info "Running $n_shots shots with irregular free surface"
    
    for shot_id in 1:n_shots
        # Create source for this shot
        src = _create_source_for_shot(shot_config, shot_id, M, backend)
        
        # Run simulation
        result = run_shot_irregular!(backend, W, M, H, a, rec_template, src, params, surface;
                                      shot_id=shot_id, show_progress=show_progress, on_step=on_step)
        
        results[shot_id] = result
        
        # Callback
        if on_shot_complete !== nothing
            on_shot_complete(result)
        end
        
        @info "Shot $shot_id/$n_shots completed"
    end
    
    return results
end

# ==============================================================================
# Automated run functions
# ==============================================================================

"""
    run_shots_auto_irregular!(model, z_surface, rec_x, rec_depth, 
                              src_x, src_depth, wavelet, params; kwargs...)

Automated multi-shot simulation(irregular surface version)

# Arguments
- `model`: VelocityModel
- `z_surface`: Surface elevation array
- `rec_x`, `rec_depth`: Receiver positions(x coordinate and depth relative to surface)
- `src_x`, `src_depth`: Source position array
- `wavelet`: Source wavelet
- `params`: SimParams

# Keyword Arguments
- `nbc`: absorbing boundary layers (default: 50)
- `fd_order`: FD order (default: 8)
- `n_iter`: IBMnumber of iterations (default: 20)
- `output_dir`: Output directory
"""
function run_shots_auto_irregular!(model::VelocityModel,
                                    z_surface::Vector{<:Real},
                                    rec_x::Vector{<:Real}, rec_depth::Vector{<:Real},
                                    src_x::Vector{<:Real}, src_depth::Vector{<:Real},
                                    wavelet::Vector{Float32}, params::SimParams;
                                    nbc::Int=50, fd_order::Int=8, n_iter::Int=20,
                                    output_dir::String="outputs/",
                                    on_shot_complete=nothing)
    
    # Select backend
    be = is_cuda_available() ? backend(:cuda) : backend(:cpu)
    @info "Using backend: $(typeof(be))"
    
    # Initialize medium
    medium = init_medium(model, nbc, fd_order, be; free_surface=false)  # do not use regular free surface
    
    # Initializeirregular surface
    surface_cpu = init_irregular_surface(Float32.(z_surface), medium; n_iter=n_iter, backend=CPU_BACKEND)
    
    # If GPU, convert to GPU format
    surface = be isa CUDABackend ? to_gpu(surface_cpu) : surface_cpu
    
    # Initialize other components
    habc = init_habc(medium.nx, medium.nz, nbc, params.dt, model.dx, model.dz, 
                     sum(model.vp)/length(model.vp), be)
    fd_coeffs = to_device(get_fd_coefficients(fd_order), be)
    wavefield = Wavefield(medium.nx, medium.nz, be)
    
    # Setup receivers(using irregular surface version)
    rec_template = setup_irregular_receivers(Float32.(rec_x), Float32.(rec_depth),
                                              surface_cpu, medium, params.nt;
                                              type=:vz, backend=be)
    
    # Create output directory
    mkpath(output_dir)
    
    # Run all shots
    results = Vector{ShotResult}()
    n_shots = length(src_x)
    
    for shot_id in 1:n_shots
        # Setup source for this shot
        src = setup_irregular_source(Float32(src_x[shot_id]), Float32(src_depth[shot_id]),
                                      wavelet, surface_cpu, medium; backend=be)
        
        # Run simulation
        result = run_shot_irregular!(be, wavefield, medium, habc, fd_coeffs,
                                      rec_template, src, params, surface;
                                      shot_id=shot_id, show_progress=true)
        
        push!(results, result)
        
        # Save results
        if on_shot_complete !== nothing
            on_shot_complete(result)
        end
        
        # Default save
        save_gather(result, joinpath(output_dir, "shot_$(shot_id).bin"))
        
        @info "Shot $shot_id/$n_shots completed"
    end
    
    return results
end
