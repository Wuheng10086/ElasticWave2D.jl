# Fomo.jl User Manual

This manual provides detailed instructions on how to use Fomo.jl for elastic wave simulation, covering everything from basic API usage to advanced boundary conditions and performance optimization.

## Table of Contents
1. [Core Concepts](#core-concepts)
2. [Simulation Workflow](#simulation-workflow)
3. [API Reference](#api-reference)
4. [Boundary Conditions](#boundary-conditions)
5. [Demos & Examples](#demos--examples)

---

## Core Concepts

### Staggered Grid
Fomo.jl uses a standard Virieux (1986) staggered grid for high stability and accuracy.
- **Stress Nodes ($T_{xx}, T_{zz}$)**: Located at integer coordinates $(i, j)$.
- **Velocity Nodes ($V_x$)**: Located at $(i+1/2, j)$.
- **Velocity Nodes ($V_z$)**: Located at $(i, j+1/2)$.
- **Shear Stress ($T_{xz}$)**: Located at $(i+1/2, j+1/2)$.

### Units
- **Distance**: Meters (m)
- **Time**: Seconds (s)
- **Velocity**: Meters per second (m/s)
- **Density**: Kilograms per cubic meter (kg/m³)
- **Frequency**: Hertz (Hz)

---

## Simulation Workflow

A typical simulation consists of three steps:
1.  **Model Definition**: Create velocity and density models.
2.  **Geometry Setup**: Define source and receiver locations.
3.  **Execution**: Run the simulation using `seismic_survey` or `simulate!`.

### 1. Model Definition
Use the `VelocityModel` struct to define your medium.
```julia
# Create 200x100 grid with 10m spacing
vp = fill(3000.0f0, 100, 200)  # Note: [nz, nx] dimensions!
vs = fill(1800.0f0, 100, 200)
rho = fill(2200.0f0, 100, 200)
model = VelocityModel(vp, vs, rho, 10.0f0, 10.0f0; name="my_model")
```

### 2. Geometry Setup
Sources and receivers are defined by their coordinates.
```julia
# Single source at x=1000m, z=20m
sources = [(1000.0f0, 20.0f0)]

# Line of receivers at surface
rec_x = Float32.(collect(0:10:2000))
rec_z = fill(0.0f0, length(rec_x))
```

### 3. Execution
Run the survey.
```julia
seismic_survey(model, sources, (rec_x, rec_z);
    simulate_surface_waves = true,
    config = SimulationConfig(nt=1000, f0=25.0f0)
)
```

---

## API Reference

### `VelocityModel`
Stores physical properties of the medium.
- `vp`: P-wave velocity matrix `[nz, nx]`
- `vs`: S-wave velocity matrix `[nz, nx]`
- `rho`: Density matrix `[nz, nx]`
- `dx`, `dz`: Grid spacing

### `SimulationConfig`
Controls simulation parameters.
- `nt`: Number of time steps (Integer).
- `f0`: Source center frequency (Hz).
- `dt`: Time step (Optional, auto-calculated from CFL if `nothing`).
- `output_dir`: Directory to save results.
- `save_gather`: Save seismic data to binary/image (Bool).
- `plot_gather`: Generate plots of shot gathers (Bool).
- `show_progress`: Display progress bar (Bool).

### `VideoConfig`
Controls wavefield animation recording.
- `fields`: Array of symbols to record (e.g., `[:vz, :vx, :p, :vel]`).
- `fps`: Frames per second for the output video.
- `skip`: Number of time steps to skip between frames.
- `colormap`: Symbol for Makie colormap (e.g., `:seismic`, `:inferno`).

---

## Boundary Conditions

Fomo.jl supports three primary boundary behaviors:

1.  **Absorbing (HABC)**:
    -   Used when `simulate_surface_waves = false`.
    -   Simulates an infinite medium by absorbing outgoing waves at all boundaries (Top, Bottom, Left, Right).
    -   **Note**: The code automatically pads the model upwards to ensure sources are not too close to the absorbing boundary.

2.  **Explicit Free Surface**:
    -   Used when `simulate_surface_waves = true` (and `free_surface=true` in config).
    -   Simulates the top boundary as a stress-free surface ($T_{zz}=0, T_{xz}=0$).
    -   Generates strong Rayleigh waves and surface reflections.

3.  **Vacuum Formulation (Topography/Cavities)**:
    -   Achieved by setting `rho = 0`, `vp = 0`, `vs = 0` in the model.
    -   The solver automatically handles the vacuum-solid interface as a free surface.
    -   Allows for arbitrary topography (mountains, valleys) and internal voids (caves, tunnels).

---

## Demos & Examples

### 1. Elastic Wave Demo (`examples/elastic_wave_demo.jl`)
**Purpose**: Visualize detailed wave propagation physics.
-   **Scenario**: A two-layer medium with a high-velocity bottom layer.
-   **Features**:
    -   High-resolution grid (0.5m spacing).
    -   Records Vz wavefield movies.
    -   Demonstrates reflection, transmission, and P-to-S mode conversion.

### 2. Seismic Survey Demo (`examples/seismic_survey_demo.jl`)
**Purpose**: Compare different acquisition scenarios.
-   **Scenarios**:
    1.  **Absorbing**: No surface waves, only body waves.
    2.  **Free Surface**: Strong ground roll (Rayleigh waves).
    3.  **Vacuum Flat**: Using vacuum layer to simulate surface.
    4.  **Vacuum Topography**: Using vacuum to simulate irregular terrain.

### 3. Batch Simulation Demo (`examples/batch_simulation_demo.jl`)
**Purpose**: High-performance production run.
-   **Features**:
    -   Simulates multiple shots sequentially.
    -   Optimized for throughput (minimal I/O).
    -   Pre-allocates memory for efficiency.
