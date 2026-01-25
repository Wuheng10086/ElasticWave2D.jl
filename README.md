# Fomo.jl

**Fomo.jl** (Forward Modeling) is a high-performance 2D elastic wave equation solver written in Julia. It is designed for seismic exploration and geophysics research, featuring GPU acceleration (CUDA), flexible boundary conditions, and a vacuum formulation for handling complex topography and internal voids.

## Key Features

*   **High Performance**:
    *   Optimized finite-difference kernels (up to 8th order spatial accuracy).
    *   Seamless GPU acceleration via `CUDA.jl`.
    *   Staggered grid formulation (Virieux, 1986).
*   **Advanced Boundary Handling**:
    *   **HABC**: Higdon Absorbing Boundary Conditions for efficient wave absorption.
    *   **Free Surface**: Explicit free surface implementation or Vacuum formulation.
    *   **Vacuum Formulation**: Handle arbitrary topography and internal cavities (cracks, tunnels) by setting density to zero.
*   **User-Friendly API**:
    *   Simple `seismic_survey` interface for generating shot gathers.
    *   `simulate!` for detailed wavefield modeling.
    *   Built-in visualization tools (Makie-based) for models, geometry, and wavefield movies.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/Fomo.jl")
```

Or for development:
```bash
git clone https://github.com/yourusername/Fomo.jl
cd Fomo.jl
julia --project=.
```

## Quick Start

### 1. Run a Demo
The best way to get started is to run the included examples:

```bash
# High-resolution elastic wave propagation in a two-layer medium
julia -t 4 examples/elastic_wave_demo.jl

# Generate seismic shot gathers with different boundary conditions
julia -t 4 examples/seismic_survey_demo.jl

# High-performance batch simulation of multiple shots
julia -t 4 examples/batch_simulation_demo.jl
```

### 2. Simple Simulation Script

```julia
using Fomo

# 1. Define Model (200x100 grid, 10m spacing)
dx, dz = 10.0f0, 10.0f0
nx, nz = 200, 100
vp = fill(3000.0f0, nz, nx)
vs = fill(1800.0f0, nz, nx)
rho = fill(2200.0f0, nz, nx)

model = VelocityModel(vp, vs, rho, dx, dz; name="simple_model")

# 2. Define Source and Receivers
src_x, src_z = 1000.0f0, 20.0f0
rec_x = Float32.(collect(100:10:1900))
rec_z = fill(2.0f0, length(rec_x))

# 3. Run Simulation
seismic_survey(model, [(src_x, src_z)], (rec_x, rec_z);
    simulate_surface_waves = true, # Enable free surface
    config = SimulationConfig(
        nt = 1000,
        f0 = 20.0f0,
        output_dir = "outputs/my_test"
    )
)
```

## Documentation

For detailed usage instructions, API reference, and physics explanation, please refer to the [User Manual](docs/user_manual.md).

## Project Structure

*   `src/`: Source code.
    *   `kernels/`: Low-level finite difference and physics kernels (CPU/GPU).
    *   `simulation/`: High-level simulation logic and time stepping.
    *   `visualization/`: Plotting and video generation tools.
*   `examples/`: Ready-to-run demonstration scripts.
*   `outputs/`: Default directory for simulation results.
