# Fomo.jl - API Reference

## Main Module Exports

### Core Types
- `VelocityModel` - Represents the velocity model with VP, VS, and density
- `SimulationResult` - Contains simulation results including gather data

### Configuration Structures
- `SimulationConfig` - Main simulation configuration
- `VideoConfig` - Video recording configuration

### Model Creation Functions
- `create_homogeneous_model(vp, vs, rho, (nz, nx), dx; name="")` - Creates a homogeneous velocity model
- `create_layered_model(layers, dx; name="")` - Creates a layered velocity model
- `create_gradient_model(vp_func, vs_func, rho_func, (nz, nx), dx; name="")` - Creates a gradient velocity model

### Simulation Functions
- `simulate!(model, src_x, src_z, rec_x, rec_z; config, video_config=nothing)` - Runs the simulation

### Surface Shape Generators
- `flat_surface(nx, dx, depth)` - Creates a flat surface
- `sinusoidal_surface(nx, dx; base_depth=50, amplitude=25, wavelength=1000)` - Creates a sinusoidal surface
- `gaussian_valley(nx, dx; base_depth=50, valley_depth=30, width=200)` - Creates a gaussian valley
- `gaussian_hill(nx, dx; base_depth=50, hill_height=30, width=200)` - Creates a gaussian hill
- `tilted_surface(nx, dx; depth_left=50, depth_right=80)` - Creates a tilted surface
- `step_surface(nx, dx; depth_left=50, depth_right=80)` - Creates a step surface
- `random_surface(nx, dx; amplitude=20, smoothness=2.0)` - Creates a random surface
- `combine_surfaces(s1, s2, ...)` - Combines multiple surface shapes

### Utility Functions
- `model_info(model)` - Prints information about a velocity model

---

## Detailed API Documentation

### SimulationConfig
Configuration for the main simulation:

```julia
SimulationConfig(;
    nt = 1000,              # Number of time steps
    f0 = 15.0f0,            # Source frequency (Hz)
    dt = nothing,            # Time step (computed automatically if nothing)
    cfl = 0.5,              # CFL number
    nbc = 50,               # Absorbing boundary layers
    fd_order = 8,           # Finite difference order
    free_surface = true,    # Enable free surface condition (controls Rayleigh waves)
    output_dir = "outputs", # Output directory
    save_gather = true,     # Save seismic gather
    plot_gather = false,    # Create gather plot
    show_progress = true    # Show progress bar
)
```

The `free_surface` parameter controls whether surface waves (like Rayleigh waves) are generated:
- `free_surface = true`: Enables free surface condition, allowing surface waves to propagate
- `free_surface = false`: Disables free surface condition, suppressing surface waves
```

### VideoConfig
Configuration for video recording:

```julia
VideoConfig(;
    fields = [:vz],         # Fields to record (:vx, :vz, :vel, :p, :txx, :tzz, :txz)
                            # Can specify multiple fields for multi-field visualization
    skip = 10,              # Record every N steps
    downsample = 1,         # Spatial downsampling factor
    colormap = :seismic,    # Color scheme
    fps = 30,               # Video frame rate
    show_boundary = true    # Show physical domain boundary
)
```

### VelocityModel
Represents a velocity model with P-wave velocity, S-wave velocity, and density:

```julia
VelocityModel(vp, vs, rho, dx, dz; name="")
```

Where:
- `vp` - P-wave velocity matrix (nz × nx)
- `vs` - S-wave velocity matrix (nz × nx)
- `rho` - Density matrix (nz × nx)
- `dx`, `dz` - Spatial sampling intervals
- `name` - Optional model name

### Simulation Result
The `simulate!` function returns a `SimulationResult` object with the following fields:
- `gather` - Seismic gather data (nt × nreceivers)
- `dt` - Time sampling interval
- `nt` - Number of time samples
- `video_file` - Path to video file (if recorded)
- `gather_file` - Path to gather data file
- `gather_plot` - Path to gather plot (if created)

---

## Example Usage Patterns

### Basic Simulation
```julia
using Fomo

# Create model
model = create_homogeneous_model(3000.0, 1800.0, 2200.0, (100, 200), 10.0)

# Set up geometry
src_x, src_z = 1000.0f0, 50.0f0
rec_x = Float32.(collect(100.0:20.0:1900.0))
rec_z = fill(10.0f0, length(rec_x))

# Configure simulation
config = SimulationConfig(
    nt = 1000,
    f0 = 15.0f0,
    output_dir = "outputs/basic"
)

# Run simulation
result = simulate!(model, src_x, src_z, rec_x, rec_z; config=config)
```

### Simulation with Video
```julia
using Fomo

# ... model and geometry setup ...

# Configure simulation and video
sim_config = SimulationConfig(
    nt = 2000,
    f0 = 15.0f0,
    output_dir = "outputs/video"
)

video_config = VideoConfig(
    fields = [:vz],
    skip = 10,
    fps = 30
)

# Run simulation with video
result = simulate!(model, src_x, src_z, rec_x, rec_z;
    config = sim_config,
    video_config = video_config)
```

### Layered Model
```julia
using Fomo

# Define layers
layers = [
    (thickness=80.0, vp=2500.0, vs=1500.0, rho=2000.0),  # Sediment
    (thickness=150.0, vp=3500.0, vs=2100.0, rho=2300.0), # Basement
    (thickness=200.0, vp=4500.0, vs=2600.0, rho=2500.0)  # Bedrock
]

model = create_layered_model(layers, 10.0; name="three_layer")
```