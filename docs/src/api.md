# API Reference

## Backends

```@docs
backend
CPUBackend
CUDABackend
is_cuda_available
```

## Data Structures

```@docs
VelocityModel
Wavefield
Medium
SimParams
```

## Initialization

```@docs
init_medium
init_habc
get_fd_coefficients
ricker_wavelet
setup_receivers
```

## Simulation

```@docs
run_shots!
run_shots_auto!
run_shots_multi_gpu!
```

## I/O

```@docs
load_model
load_model_files
save_model
save_gather
load_gather
save_geometry
load_geometry
```

## Visualization

```@docs
plot_setup
```
