# Fomo.jl - Documentation

## Overview
Fomo.jl is a high-performance 2D isotropic elastic wave simulator written in Julia. It provides a flexible and efficient framework for seismic forward modeling with support for CPU and GPU computation.

## Key Features
- Backend-dispatched architecture (same code runs on CPU or GPU)
- High-order staggered-grid finite differences (2nd to 10th order)
- Hybrid absorbing boundary conditions (HABC)
- Free surface modeling for accurate Rayleigh wave simulation
- Multi-GPU parallel processing with automatic load balancing
- Support for multiple data formats (SEG-Y, Binary, MAT, NPY, HDF5, JLD2)
- Real-time wavefield visualization and video recording

## Installation
```julia
using Pkg
Pkg.add(url="https://github.com/Wuheng10086/Fomo.jl")
```

## Core Components

### 1. Model Creation
Three main types of velocity models can be created:

#### Homogeneous Model
```julia
model = create_homogeneous_model(vp, vs, rho, (nz, nx), dx; name="model_name")
```

#### Layered Model
```julia
layers = [
    (thickness=80.0, vp=2500.0, vs=1500.0, rho=2000.0),  # Top layer
    (thickness=150.0, vp=3500.0, vs=2100.0, rho=2300.0), # Middle layer
    (thickness=200.0, vp=4500.0, vs=2600.0, rho=2500.0)  # Bottom layer
]
model = create_layered_model(layers, dx; name="layered_model")
```

#### Gradient Model
```julia
vp_func(z, x) = 2000.0 + 1.2*z + 0.1*x
vs_func(z, x) = 1200.0 + 0.7*z + 0.05*x
rho_func(z, x) = 2000.0 + 0.25*z

model = create_gradient_model(vp_func, vs_func, rho_func, (nz, nx), dx; name="gradient_model")
```

### 2. Simulation Configuration
The main simulation parameters are configured using the `SimulationConfig` struct:

```julia
config = SimulationConfig(
    nt = 2000,              # Number of time steps
    f0 = 15.0f0,            # Source frequency (Hz)
    nbc = 50,               # Absorbing boundary layers
    fd_order = 8,           # Finite difference order
    free_surface = true,    # Enable free surface
    output_dir = "outputs", # Output directory
    save_gather = true,     # Save seismic gather
    plot_gather = true      # Create gather plot
)
```

### 3. Video Recording
To record wavefield videos, use the `VideoConfig` struct:

```julia
video_config = VideoConfig(
    fields = [:vz],         # Fields to record (:vx, :vz, :vel, :p)
    skip = 10,              # Record every N steps
    fps = 30,               # Video frame rate
    colormap = :seismic     # Color scheme
)
```

### 4. Running Simulations
The main simulation function is `simulate!`:

```julia
result = simulate!(
    model,                  # Velocity model
    src_x, src_z,          # Source coordinates
    rec_x, rec_z;          # Receiver coordinates
    config = config,       # Simulation configuration
    video_config = video_config  # Optional video configuration
)
```

## Surface Shape Generators
Fomo.jl provides several functions to create surface topographies:

- `flat_surface(nx, dx, depth)` - Flat surface at constant depth
- `sinusoidal_surface(nx, dx; base_depth=50, amplitude=25, wavelength=1000)` - Sinusoidal surface
- `gaussian_valley(nx, dx; base_depth=50, valley_depth=30, width=200)` - Gaussian depression
- `gaussian_hill(nx, dx; base_depth=50, hill_height=30, width=200)` - Gaussian elevation
- `tilted_surface(nx, dx; depth_left=50, depth_right=80)` - Linear tilt
- `step_surface(nx, dx; depth_left=50, depth_right=80)` - Step/cliff
- `random_surface(nx, dx; amplitude=20, smoothness=2.0)` - Random rough surface
- `combine_surfaces(s1, s2, ...)` - Combine multiple shapes

## Performance Tips
1. Use GPU backend when available for significant speedup
2. Choose appropriate finite difference order (8th order is typically optimal)
3. Balance model size with available memory
4. Use appropriate time step based on stability criteria
5. Consider using multi-GPU parallelization for large models

## File Organization
```
Fomo.jl/
├── src/
│   ├── Fomo.jl                 # Main module
│   ├── backends/               # CPU/CUDA abstraction
│   │   └── backend.jl
│   ├── types/                  # Data structures
│   │   ├── model.jl            # Velocity model definition
│   │   └── structures.jl       # Simulation structures
│   ├── kernels/                # Finite difference kernels
│   │   ├── velocity.jl         # Velocity update kernels
│   │   ├── stress.jl           # Stress update kernels
│   │   ├── boundary.jl         # Boundary condition kernels
│   │   ├── source_receiver.jl  # Source and receiver kernels
│   │   └── vacuum.jl           # Vacuum formulation kernels
│   ├── simulation/             # Simulation logic
│   │   ├── api.jl              # High-level API
│   │   ├── simple_api.jl       # Simplified API
│   │   ├── init.jl             # Initialization routines
│   │   ├── init_vacuum.jl      # Vacuum initialization
│   │   ├── time_stepper.jl     # Time stepping routines
│   │   ├── shots.jl            # Shot processing
│   │   ├── batch.jl            # Batch processing
│   │   └── parallel.jl         # Parallel processing
│   ├── io/                     # Input/output operations
│   │   ├── model_io.jl         # Model I/O
│   │   ├── gather_io.jl        # Gather data I/O
│   │   └── geometry_io.jl      # Geometry I/O
│   └── visualization/          # Plotting & video
│       ├── plots.jl            # Static plots
│       └── video.jl            # Video generation
├── examples/
│   ├── run_demo.jl
│   ├── run_irregular_with_video.jl
│   └── run_vacuum_topography.jl
├── README.md
├── README_zh.md
└── Project.toml
```

## Examples
Several example scripts are provided in the `examples/` directory:
- `quick_start_demo.jl` - Basic usage example
- `comprehensive_demo.jl` - Complete feature demonstration
- `run_demo.jl` - Standard demonstration
- `run_irregular_with_video.jl` - Irregular surface with video
- `run_vacuum_topography.jl` - Vacuum topography modeling

## Troubleshooting
1. If CUDA is not available, the code will automatically fall back to CPU
2. Ensure sufficient memory for large models
3. Check that source and receiver coordinates are within model bounds
4. Verify that time step satisfies stability criteria
5. For video recording, ensure required video encoding libraries are installed