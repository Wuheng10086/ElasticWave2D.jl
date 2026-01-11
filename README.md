# Fomo.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.9%20|%201.10%20|%201.11-blue)](https://julialang.org/)

[中文文档](README_zh.md) | [English](README.md)

**Fomo** - **Fo**rward **Mo**deling: High-performance 2D isotropic elastic wave simulator in Julia.

```ibm_method=:mirror```  is **not stable** yet, please wait for the next release.

## ✨ Features

- 🚀 **Backend-dispatched architecture** - Same code runs on CPU or GPU
- 📐 **High-order staggered-grid FD** - 2nd to 10th order spatial accuracy
- 🛡️ **Hybrid absorbing boundary (HABC)** - Effective boundary reflection suppression
- 🌊 **Free surface modeling** - Accurate Rayleigh wave simulation
- 🏔️ **Irregular topography (IBM)** - Immersed Boundary Method for complex surfaces
- ⚡ **Multi-GPU parallel** - Automatic load balancing
- 📁 **Multiple formats** - SEG-Y, Binary, MAT, NPY, HDF5, JLD2
- 🎬 **Video recording** - Real-time wavefield visualization

## 📋 Requirements

- **Julia 1.9, 1.10, or 1.11** (1.12 not yet supported due to CairoMakie compatibility)
- CUDA GPU (optional, for GPU acceleration)

## 🔧 Installation

### From GitHub

```julia
using Pkg
Pkg.add(url="https://github.com/Wuheng10086/Fomo.jl")
```

### Local Development

```bash
git clone https://github.com/Wuheng10086/Fomo.jl.git
cd Fomo.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## 🚀 Quick Start

### High-level API (Recommended)

```julia
using Fomo

# Create velocity model
nx, nz = 400, 200
dx = 10.0f0

vp = fill(3000.0f0, nz, nx)
vs = fill(1800.0f0, nz, nx)
rho = fill(2200.0f0, nz, nx)

vp[100:end, :] .= 4000.0f0
vs[100:end, :] .= 2400.0f0

model = VelocityModel(vp, vs, rho, dx, dx; name="Two-layer model")

# Run simulation (without video)
result = simulate!(
    model,
    2000.0f0, 50.0f0,                    # source (x, z) in meters
    Float32.(100:20:3900),               # receiver x positions
    fill(10.0f0, 190);                   # receiver z positions
    config = SimulationConfig(nt=3000, f0=15.0f0, output_dir="outputs")
)

# Run simulation (with video) - VideoConfig is a separate parameter
result = simulate!(
    model,
    2000.0f0, 50.0f0,
    Float32.(100:20:3900),
    fill(10.0f0, 190);
    config = SimulationConfig(nt=3000, f0=15.0f0, output_dir="outputs"),
    video_config = VideoConfig(fields=[:vz], skip=5, fps=30)
)
```

### Irregular Free Surface

```julia
using Fomo

model = VelocityModel(vp, vs, rho, dx, dx)

# Define surface shape using helper functions
z_surface = sinusoidal_surface(nx, dx; base_depth=50, amplitude=30, wavelength=1000)

# Or combine multiple shapes
z_surface = combine_surfaces(
    sinusoidal_surface(nx, dx; amplitude=20),
    gaussian_valley(nx, dx; valley_depth=25, width=300)
)

# Run simulation
result = simulate_irregular!(
    model,
    z_surface,                           # your surface shape
    2000.0f0,                            # source x position
    Float32.(100:20:3900);               # receiver x positions
    config = IrregularSurfaceConfig(
        nt = 3000,
        ibm_method = :direct_zero,       # or :mirror for higher accuracy
        src_depth = 30.0f0,              # depth below surface
        output_dir = "outputs_irregular"
    ),
    video_config = VideoConfig(fields=[:vz], skip=10)
)
```

### Surface Shape Helpers

| Function | Description |
|----------|-------------|
| `flat_surface(nx, dx, depth)` | Flat surface at constant depth |
| `sinusoidal_surface(nx, dx; amplitude, wavelength)` | Sinusoidal surface |
| `gaussian_valley(nx, dx; valley_depth, width)` | Gaussian depression |
| `gaussian_hill(nx, dx; hill_height, width)` | Gaussian elevation |
| `tilted_surface(nx, dx; depth_left, depth_right)` | Linear tilt |
| `step_surface(nx, dx; depth_left, depth_right)` | Step/cliff |
| `random_surface(nx, dx; amplitude, smoothness)` | Random rough surface |
| `combine_surfaces(s1, s2, ...)` | Combine multiple shapes |

## 📂 Examples

### `examples/run_demo.jl` - Comprehensive Demo

Demonstrates three scenarios:

```julia
using Fomo

# Demo 1: Quick test (homogeneous model, no video)
result1 = simulate!(model, src_x, src_z, rec_x, rec_z;
    config = SimulationConfig(nt=1000, output_dir="outputs/demo1"))

# Demo 2: Surface waves visualization (with video)
result2 = simulate!(model2, src_x, src_z, rec_x, rec_z;
    config = SimulationConfig(nt=4000, f0=20.0f0, output_dir="outputs/demo2"),
    video_config = VideoConfig(fields=[:vz], skip=5, fps=30))

# Demo 3: Irregular surface with IBM
z_surface = combine_surfaces(
    sinusoidal_surface(nx, dx; base_depth=50, amplitude=25),
    gaussian_valley(nx, dx; valley_depth=20, width=250)
)
result3 = simulate_irregular!(model3, z_surface, src_x, rec_x;
    config = IrregularSurfaceConfig(nt=3000, ibm_method=:direct_zero),
    video_config = VideoConfig(fields=[:vz], skip=10))
```

Run:
```bash
julia --project=. examples/run_demo.jl
```

### `examples/run_regular_surface_video.jl` - Surface Wave Visualization

Two-layer model showing P-wave, S-wave, and Rayleigh wave:

```julia
using Fomo

# Two-layer model
vp[1:160, :] .= 2500.0f0    # Upper layer
vp[161:end, :] .= 4000.0f0  # Lower layer

model = VelocityModel(vp, vs, rho, 5.0f0, 5.0f0)

video_cfg = VideoConfig(fields=[:vz], skip=5, fps=30, colormap=:seismic)

result = simulate!(model, nx*dx/2, 50.0f0, rec_x, rec_z;
    config = SimulationConfig(nt=4000, f0=20.0f0, output_dir="outputs_regular"),
    video_config = video_cfg)
```

Expected waves in video:
- **P-wave**: ~2500 m/s (fastest)
- **S-wave**: ~1500 m/s
- **Rayleigh wave**: ~1380 m/s (along surface)

Run:
```bash
julia --project=. examples/run_regular_surface_video.jl
```

### `examples/run_irregular_with_video.jl` - Irregular Topography

Demonstrates custom surface shapes:

```julia
using Fomo

# Example 1: Using helper functions
z_surface = sinusoidal_surface(nx, dx; base_depth=50, amplitude=25, wavelength=1500)

# Example 2: Fully custom shape
x = Float32.((0:nx-1) .* dx)
z_custom = Float32.(60.0 .+ 20.0 .* sin.(2π .* x ./ 1500.0) .+
                    10.0 .* sin.(2π .* x ./ 300.0))

# Add a canyon
for i in 1:nx
    xi = (i-1) * dx
    if 1300 < xi < 1700
        z_custom[i] += 30.0f0 * (1 - abs(xi - 1500) / 200)
    end
end

result = simulate_irregular!(model, z_custom, src_x, rec_x;
    config = IrregularSurfaceConfig(nt=3000, src_depth=40.0f0),
    video_config = VideoConfig(fields=[:vz], skip=10))
```

Run:
```bash
julia --project=. examples/run_irregular_with_video.jl
```

## 🏔️ Immersed Boundary Method (IBM)

The IBM enables accurate modeling of complex topography without fine grid refinement:

| Method | Grid Size | Time Steps | Total Cost |
|--------|-----------|------------|------------|
| Fine grid + staircase | 4N | 2T | **8×** |
| **IBM (this package)** | N | T | **~1.08×** |

Two IBM methods available:
- `:direct_zero` - Stable, recommended for most cases
- `:mirror` - Higher accuracy, may require smaller time step

## 📂 Project Structure

```
Fomo.jl/
├── src/
│   ├── Fomo.jl                 # Main module
│   ├── backends/               # CPU/CUDA abstraction
│   ├── types/                  # Data structures
│   ├── kernels/                # FD kernels
│   │   ├── velocity.jl
│   │   ├── stress.jl
│   │   ├── boundary.jl
│   │   └── ibm.jl
│   ├── surface/                # Irregular surface
│   ├── simulation/             # Time stepping
│   │   ├── api.jl              # High-level API
│   │   ├── time_stepper.jl
│   │   └── time_stepper_ibm.jl
│   ├── io/                     # Model/gather I/O
│   └── visualization/          # Plotting & video
├── examples/
│   ├── run_demo.jl
│   ├── run_regular_surface_video.jl
│   └── run_irregular_with_video.jl
└── README.md
```

## 📚 API Reference

### Configuration Structs

```julia
SimulationConfig(
    nt = 3000,              # Number of time steps
    f0 = 15.0f0,            # Source frequency (Hz)
    nbc = 50,               # Absorbing boundary layers
    fd_order = 8,           # Finite difference order
    free_surface = true,    # Enable free surface
    output_dir = "outputs"  # Output directory
)

IrregularSurfaceConfig(
    nt = 3000,
    f0 = 15.0f0,
    ibm_method = :direct_zero,  # :direct_zero or :mirror
    ibm_iterations = 3,
    src_depth = 30.0f0,         # Source depth below surface
    rec_depth = 0.0f0,          # Receiver depth (0 = on surface)
    output_dir = "outputs"
)

VideoConfig(
    fields = [:vz],         # Fields to record (:vx, :vz, :vel, :p)
    skip = 10,              # Record every N steps
    fps = 30,               # Video frame rate
    colormap = :seismic     # Color scheme
)
```

### High-level Functions

| Function | Description |
|----------|-------------|
| `simulate!(model, src_x, src_z, rec_x, rec_z; config, video_config)` | Regular surface simulation |
| `simulate_irregular!(model, z_surface, src_x, rec_x; config, video_config)` | Irregular surface simulation |

### Surface Helpers

| Function | Description |
|----------|-------------|
| `flat_surface`, `sinusoidal_surface`, `gaussian_valley`, `gaussian_hill` | Basic shapes |
| `tilted_surface`, `step_surface`, `random_surface` | More shapes |
| `combine_surfaces(s1, s2, ...; method=:add)` | Combine shapes |

## 📖 References

1. Luo, Y., & Schuster, G. (1990). Parsimonious staggered grid finite-differencing of the wave equation. *Geophysical Research Letters*, 17(2), 155-158.

2. Ren, Z., & Liu, Y. (2014). Numerical modeling of the first-order elastic equations with the hybrid absorbing boundary condition. *Chinese Journal of Geophysics*, 57(2), 595-606. doi:10.6038/cjg20140223

3. Li, X., Yao, G., Niu, F., Wu, D., & Liu, N. (2023). Waveform inversion of seismic first arrivals acquired on irregular surface. *Geophysics*, 88(3), R289-R302.

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

## 👤 Author

Wuheng - 2025
