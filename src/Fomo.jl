# ==============================================================================
# Fomo.jl - Forward Modeling
#
# High-Performance 2D Elastic Wave Simulation Framework
#
# Project Structure:
# ==================
# src/
# ├── Fomo.jl                 # Main module
# ├── backends/               # Hardware abstraction
# │   └── backend.jl          # CPU/CUDA backend
# ├── types/                  # Data structures
# │   ├── structures.jl       # Core types (Wavefield, Medium, Source, etc.)
# │   └── model.jl            # VelocityModel and related
# ├── kernels/                # Compute kernels (hot loops)
# │   ├── velocity.jl         # Velocity update
# │   ├── stress.jl           # Stress update  
# │   ├── boundary.jl         # HABC, free surface
# │   ├── source_receiver.jl  # Source injection, receiver recording
# │   └── ibm.jl              # Immersed Boundary Method
# ├── surface/                # Irregular surface handling
# │   └── irregular.jl        # Surface initialization and helpers
# ├── simulation/             # Simulation control
# │   ├── init.jl             # Medium/wavefield initialization
# │   ├── time_stepper.jl     # Time stepping (regular)
# │   ├── time_stepper_ibm.jl # Time stepping (irregular surface)
# │   ├── shots.jl            # Shot management
# │   └── parallel.jl         # Multi-GPU parallel execution
# ├── io/                     # Input/Output
# │   ├── model_io.jl         # Model load/save
# │   ├── gather_io.jl        # Gather/results save/load
# │   └── geometry_io.jl      # Survey geometry
# └── visualization/          # Plotting and video
#     ├── video.jl            # Wavefield video recording
#     └── plots.jl            # Static plots
#
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
using CUDA

# ==============================================================================
# CUDA Support
# ==============================================================================

const CUDA_AVAILABLE = Ref(false)

function __init__()
    if CUDA.functional()
        CUDA_AVAILABLE[] = true
        @info "Fomo: CUDA functional, GPU acceleration enabled"
    else
        @info "Fomo: CUDA not functional (no GPU), using CPU mode"
    end
end

is_cuda_available() = CUDA_AVAILABLE[]
is_cuda_functional() = CUDA_AVAILABLE[]

# ==============================================================================
# Exports
# ==============================================================================

# --- Backend ---
export AbstractBackend, CPUBackend, CUDABackend
export CPU_BACKEND, CUDA_BACKEND
export backend, to_device, synchronize
export is_cuda_available, is_cuda_functional

# --- Types ---
export Wavefield, Medium, HABCConfig
export Source, Receivers, SimParams
export ShotConfig, MultiShotConfig, ShotResult
export VelocityModel

# --- Irregular Surface ---
export IrregularSurface, IrregularSurfaceGPU
export SurfacePoint, GhostPoint
export init_irregular_surface, init_flat_surface
export setup_irregular_source, setup_irregular_receivers, setup_surface_receivers
export is_irregular_surface_enabled, to_gpu

# --- Initialization ---
export init_medium, init_habc, setup_receivers
export get_fd_coefficients, ricker_wavelet

# --- Kernels (advanced) ---
export update_velocity!, update_stress!
export apply_habc!, apply_habc_velocity!, apply_habc_stress!
export backup_boundary!, apply_free_surface!
export apply_irregular_free_surface!
export inject_source!, record_receivers!, reset!

# --- Simulation ---
export TimeStepInfo
export time_step!, run_time_loop!
export time_step_irregular!, run_time_loop_irregular!
export run_shot!, run_shots!
export run_shot_irregular!, run_shots_irregular!
export run_shots_auto_irregular!

# --- High-level API ---
export SimulationConfig, SimulationResult, simulate!
export IrregularSurfaceConfig, simulate_irregular!

# --- Surface shape helpers ---
export flat_surface, sinusoidal_surface, gaussian_valley, gaussian_hill
export tilted_surface, step_surface, random_surface, combine_surfaces

# --- Parallel ---
export get_gpu_info, print_hardware_info
export run_shots_multi_gpu!, run_shots_auto!

# --- Visualization ---
export VideoConfig, FieldRecorder, MultiFieldRecorder
export generate_video
export plot_setup, plot_irregular_setup
export plot_surface_comparison, plot_gather_with_surface

# --- IO ---
export save_gather, load_gather
export load_model, load_model_files, save_model
export convert_model, model_info, resample_model, suggest_grid_spacing
export SurveyGeometry, MultiShotGeometry
export create_geometry, save_geometry, load_geometry

# ==============================================================================
# Include Files
# ==============================================================================

# 1. Backend (must be first - defines AbstractBackend)
include("backends/backend.jl")

# 2. Types (data structures)
include("types/structures.jl")
include("types/model.jl")

# 3. Surface (MUST be before ibm.jl - defines IrregularSurface)
include("surface/irregular.jl")

# 4. Kernels (compute-intensive code)
include("kernels/velocity.jl")
include("kernels/stress.jl")
include("kernels/boundary.jl")
include("kernels/source_receiver.jl")
include("kernels/ibm.jl")  # Uses IrregularSurface from irregular.jl

# 5. Visualization (MUST be before api.jl - defines VideoConfig)
include("visualization/video.jl")
include("visualization/plots.jl")

# 6. Simulation - ORDER IS CRITICAL HERE!
include("simulation/init.jl")
include("simulation/time_stepper.jl")      # Defines run_time_loop!
include("simulation/shots.jl")             # Uses run_time_loop!, defines MultiShotConfig
include("simulation/time_stepper_ibm.jl")  # Uses MultiShotConfig
include("simulation/parallel.jl")
include("simulation/api.jl")               # High-level API (uses VideoConfig)

# 7. IO (file operations)
include("io/model_io.jl")
include("io/gather_io.jl")
include("io/geometry_io.jl")

end # module Fomo