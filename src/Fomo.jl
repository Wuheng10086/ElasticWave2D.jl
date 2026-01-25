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
# │   └── vacuum.jl           # Vacuum formulation for irregular surface
# ├── simulation/             # Simulation control
# │   ├── init.jl             # Medium/wavefield initialization (flat surface)
# │   ├── init_vacuum.jl      # Medium initialization with vacuum formulation
# │   ├── time_stepper.jl     # Time stepping
# │   ├── shots.jl            # Shot management
# │   ├── batch.jl            # High-performance batch simulation
# │   ├── parallel.jl         # Multi-GPU parallel execution
# │   └── api.jl              # High-level API
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
export BoundaryConfig

# --- Initialization (Flat Surface) ---
export init_medium, init_wavefield, init_habc, setup_receivers
export get_fd_coefficients, ricker_wavelet

# --- Vacuum Formulation (Irregular Surface) ---
export init_medium_vacuum
export setup_vacuum_formulation!, compute_staggered_params_vacuum
export apply_vacuum_mask!, compute_surface_indices
export setup_receivers_on_surface, setup_source_on_surface
# Surface shape generators
export flat_surface, sinusoidal_surface
export gaussian_valley, gaussian_hill, tilted_surface, step_surface
export random_surface, combine_surfaces
export validate_surface_elevation

# --- Kernels (advanced) ---
export update_velocity!, update_stress!
export apply_habc!, apply_habc_velocity!, apply_habc_stress!
export backup_boundary!, apply_free_surface!
export inject_source!, record_receivers!, reset!

# --- Simulation ---
export TimeStepInfo
export time_step!, run_time_loop!
export run_shot!, run_shots!, run_shots_fast!

# --- High-level API ---
export SimulationConfig, SimulationResult, simulate!
export SimpleConfig, TopographyConfig
export create_homogeneous_model, create_layered_model, create_gradient_model
export SourceConfig, ReceiverConfig
export seismic_survey

# --- Batch Simulation (High Performance) ---
export BatchSimulator, simulate_shot!, simulate_shots!, benchmark_shots

# --- Parallel ---
export get_gpu_info, print_hardware_info
export run_shots_multi_gpu!, run_shots_auto!

# --- Visualization ---
export VideoConfig, FieldRecorder, MultiFieldRecorder
export generate_video
export plot_setup, plot_gather

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
include("types/boundary_config.jl")

# 3. Kernels (compute-intensive code)
include("kernels/velocity.jl")
include("kernels/stress.jl")
include("kernels/boundary.jl")
include("kernels/source_receiver.jl")
include("kernels/vacuum.jl")  # Vacuum formulation for irregular surface

# 4. Visualization (MUST be before api.jl - defines VideoConfig)
include("visualization/video.jl")
include("visualization/plots.jl")

# 5. Simulation - ORDER IS CRITICAL HERE!
include("simulation/init.jl")           # Flat surface initialization
include("simulation/init_vacuum.jl")    # Vacuum formulation initialization
include("simulation/time_stepper.jl")   # Defines run_time_loop!
include("simulation/shots.jl")          # Uses run_time_loop!, defines MultiShotConfig
include("simulation/parallel.jl")
include("simulation/api.jl")            # High-level API (defines SimulationConfig)
include("simulation/batch.jl")          # High-performance batch (uses SimulationConfig)
include("simulation/simple_api.jl")     # Simplified high-level API

# 6. IO (file operations)
include("io/model_io.jl")
include("io/gather_io.jl")
include("io/geometry_io.jl")

end # module Fomo