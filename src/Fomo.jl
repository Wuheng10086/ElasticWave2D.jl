# ==============================================================================
# Fomo.jl - Forward Modeling
#
# High-Performance 2D Elastic Wave Simulation Framework
#
# Features:
# - Backend-dispatched kernels (CPU/CUDA)
# - High-order staggered-grid finite-difference
# - Hybrid Absorbing Boundary Condition (HABC)
# - Free surface modeling
# - Multi-GPU parallel execution
#
# Author: zswh
# License: MIT
# ==============================================================================

module Fomo

# ==============================================================================
# Dependencies
# ==============================================================================

using LoopVectorization
using ProgressMeter
using Printf
using Statistics
using CairoMakie
using JLD2
using JSON

# CUDA support (conditional loading)
const CUDA_AVAILABLE = Ref(false)
const CUDA_FUNCTIONAL = Ref(false)

function __init__()
    # Try to check if CUDA is functional
    try
        @eval using CUDA
        if CUDA.functional()
            CUDA_AVAILABLE[] = true
            CUDA_FUNCTIONAL[] = true
            @info "CUDA available: GPU acceleration enabled"
        else
            CUDA_AVAILABLE[] = true
            CUDA_FUNCTIONAL[] = false
            @debug "CUDA loaded but not functional"
        end
    catch e
        CUDA_AVAILABLE[] = false
        CUDA_FUNCTIONAL[] = false
        @debug "CUDA not available" exception=e
    end
end

# ==============================================================================
# Exports
# ==============================================================================

# Backend system
export AbstractBackend, CPUBackend, CUDABackend
export CPU_BACKEND, CUDA_BACKEND
export backend, to_device, synchronize

# Data structures
export Wavefield, Medium, HABCConfig
export Source, Receivers, SimParams
export ShotConfig, MultiShotConfig, ShotResult

# Initialization
export init_medium, init_habc, setup_receivers
export get_fd_coefficients, ricker_wavelet

# Kernels (for advanced users)
export update_velocity!, update_stress!
export apply_habc!, apply_habc_velocity!, apply_habc_stress!
export backup_boundary!, apply_free_surface!
export inject_source!, record_receivers!, reset!

# Simulation interface
export TimeStepInfo
export time_step!, run_time_loop!
export run_shot!, run_shots!

# Parallel execution
export get_gpu_info, print_hardware_info
export run_shots_multi_gpu!, run_shots_auto!

# Visualization
export VideoConfig, VideoRecorder, MultiFieldRecorder
export plot_setup, generate_video

# IO
export save_gather, load_gather

# Model IO
export VelocityModel, load_model, load_model_files, save_model, convert_model, model_info

# Geometry IO (for migration)
export SurveyGeometry, MultiShotGeometry
export create_geometry, save_geometry, load_geometry

# Utilities
export is_cuda_available, is_cuda_functional
is_cuda_available() = CUDA_AVAILABLE[]
is_cuda_functional() = CUDA_FUNCTIONAL[]

# ==============================================================================
# Include Files (order matters!)
# ==============================================================================

# Backend abstraction (must be first)
include("backends/backend.jl")

# Core structures
include("core/structures.jl")

# Kernels
include("kernels/velocity.jl")
include("kernels/stress.jl")
include("kernels/boundary.jl")
include("kernels/source_receiver.jl")

# IO (must be before utils/init.jl which uses VelocityModel)
include("io/output.jl")
include("io/model_loader.jl")

# Utilities (uses VelocityModel from model_loader.jl)
include("utils/init.jl")

# Geometry IO (uses ShotResult from shot_manager.jl)
include("io/geometry_io.jl")

# Simulation
include("simulation/time_stepper.jl")
include("simulation/shot_manager.jl")
include("simulation/parallel_shots.jl")

# Visualization
include("visualization/video_recorder.jl")
include("visualization/setup_check.jl")

end # module Fomo
