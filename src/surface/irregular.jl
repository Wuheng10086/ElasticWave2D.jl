# ==============================================================================
# core/irregular_surface.jl
#
# Irregular Free Surface using Immersed Boundary Method (IBM)
# 
# Supports two methods:
# - :direct_zero - Simple and stable, lower accuracy
# - :mirror      - Higher accuracy mirror interpolation (with safety checks)
#
# Implementation principles:
# 1. "Immerse" irregular free surface into regular Cartesian grid
# 2. Points above surface are ghost points, boundary conditions set via mirroring
# 3. Use iterative interpolation to compute wavefield values at mirror points
# 
# For elastic waves, free surface condition: tzz = 0, txz = 0
# ==============================================================================

"""
    SurfacePoint

Information for a single point on the surface.
"""
struct SurfacePoint
    x::Float32           # X coordinate (physical)
    z::Float32           # Z coordinate (physical, surface elevation)
    normal_x::Float32    # Normal vector x component
    normal_z::Float32    # Normal vector z component
end

"""
    GhostPoint

Complete information for a single ghost point.
"""
struct GhostPoint
    # Ghost point grid indices
    gi::Int
    gj::Int
    
    # Mirror point physical coordinates (usually on fractional grid)
    mirror_x::Float32
    mirror_z::Float32
    
    # Interpolation stencil: indices and weights of surrounding grid points
    # Using 4x4 = 16 point bicubic interpolation
    interp_i::NTuple{16, Int}      # X direction indices
    interp_j::NTuple{16, Int}      # Z direction indices  
    interp_w::NTuple{16, Float32}  # Interpolation weights
end

"""
    IrregularSurface{T}

Irregular free surface configuration.

# Fields
- `z_surface`: Surface elevation array, z_surface[i] corresponds to elevation at x = (i-1)*dx
- `ghost_points`: Information for all ghost points
- `n_ghost`: Total number of ghost points
- `n_iter`: Number of iterations (default 20)
- `enabled`: Whether irregular surface is enabled
- `method`: IBM method (:direct_zero or :mirror)
"""
struct IrregularSurface{T<:AbstractVector{Float32}}
    z_surface::T              # Surface elevation z(x)
    ghost_points::Vector{GhostPoint}
    n_ghost::Int
    n_iter::Int
    enabled::Bool
    method::Symbol            # :direct_zero or :mirror (NEW)
    
    # Grid parameters (for coordinate conversion)
    dx::Float32
    dz::Float32
    pad::Int
end

"""
    IrregularSource

Source definition for irregular surface.

# Fields
- `x`: Source x coordinate (physical)
- `depth`: Source depth (relative to surface, positive = underground)
- `wavelet`: Source wavelet
"""
struct IrregularSource{V<:AbstractVector{Float32}}
    x::Float32
    depth::Float32      # Depth relative to surface
    wavelet::V
    
    # Computed absolute grid positions (internal use)
    i::Int              # Grid x index
    j::Int              # Grid z index
end

"""
    IrregularReceivers

Receiver definition for irregular surface.

# Fields
- `x`: Receiver x coordinates array
- `depth`: Receiver depth array (relative to surface, 0 = on surface)
- `type`: Recording type (:vz, :vx, :p)
"""
struct IrregularReceivers{T<:AbstractMatrix{Float32}, V<:AbstractVector{Float32}}
    x::Vector{Float32}
    depth::Vector{Float32}
    type::Symbol
    data::T
    
    # Computed physical coordinates and interpolation info
    abs_x::V                    # Absolute x coordinates
    abs_z::V                    # Absolute z coordinates
    interp_i::Vector{NTuple{16, Int}}
    interp_j::Vector{NTuple{16, Int}}
    interp_w::Vector{NTuple{16, Float32}}
end

# ==============================================================================
# Initialization Functions
# ==============================================================================

"""
    init_irregular_surface(z_surface, medium; n_iter=20, method=:direct_zero)

Initialize irregular free surface.

# Arguments
- `z_surface`: Surface elevation array, length should equal model nx (without padding)
              z_surface[i] is the surface height at x = (i-1)*dx
- `medium`: Medium struct
- `n_iter`: Number of iterations for mirror method (default 20)
- `method`: IBM method, either :direct_zero (default, stable) or :mirror (higher accuracy)

# Returns
- `IrregularSurface` struct

# Example
```julia
# Create a simple hill topography
nx = 200
x = range(0, 2000, length=nx)
z_surface = 100 .* sin.(2pi .* x ./ 2000) .+ 50  # Undulating terrain

# Use direct zero method (stable)
surface = init_irregular_surface(z_surface, medium; method=:direct_zero)

# Or use mirror method (higher accuracy, may need tuning)
surface = init_irregular_surface(z_surface, medium; method=:mirror, n_iter=3)
```
"""
function init_irregular_surface(z_surface::Vector{<:Real}, M::Medium; 
                                 n_iter::Int=20, method::Symbol=:direct_zero,
                                 backend::AbstractBackend=CPU_BACKEND)
    
    # Validate method
    if !(method in (:direct_zero, :mirror))
        error("Unknown IBM method: $method. Use :direct_zero or :mirror")
    end
    
    dx, dz = M.dx, M.dz
    pad = M.pad
    nx_inner = M.nx - 2 * pad
    nz_inner = M.nz - 2 * pad
    fd_half = M.M  # FD half-stencil width
    
    # Validate input
    @assert length(z_surface) == nx_inner "z_surface length must equal model nx (without padding)"
    
    # Extend z_surface to full grid including padding
    z_full = zeros(Float32, M.nx)
    z_full[pad+1:pad+nx_inner] .= Float32.(z_surface)
    # Boundary extrapolation
    z_full[1:pad] .= z_full[pad+1]
    z_full[pad+nx_inner+1:end] .= z_full[pad+nx_inner]
    
    # Collect all ghost points
    ghost_list = GhostPoint[]
    
    for i in 1:M.nx
        # Surface height at this x position converted to j index
        z_surf_here = z_full[i]
        j_surf = round(Int, z_surf_here / dz) + pad + 1
        
        # Ghost points: fd_half layers above surface
        for layer in 1:fd_half
            gj = j_surf - layer
            if gj >= 1
                # Compute mirror point coordinates
                # Ghost point physical coordinates
                gx = (i - pad - 1) * dx
                gz = (gj - pad - 1) * dz
                
                # Mirror point: symmetric about surface
                # Simplified: direct vertical mirroring
                mz = 2 * z_surf_here - gz
                mx = gx  # x unchanged for vertical mirroring
                
                # Compute interpolation weights (4x4 bicubic)
                interp_i, interp_j, interp_w = _compute_interp_weights(mx, mz, dx, dz, pad, M.nx, M.nz)
                
                push!(ghost_list, GhostPoint(
                    i, gj, 
                    Float32(mx), Float32(mz),
                    interp_i, interp_j, interp_w
                ))
            end
        end
    end
    
    n_ghost = length(ghost_list)
    @info "Irregular surface initialized" n_ghost_points=n_ghost n_iterations=n_iter method=method
    
    # Move to device
    z_device = to_device(z_full, backend)
    
    return IrregularSurface(
        z_device,
        ghost_list,
        n_ghost,
        n_iter,
        true,
        method,
        Float32(dx), Float32(dz), pad
    )
end

"""
Create a disabled (flat) free surface for interface compatibility.
"""
function init_flat_surface(M::Medium)
    return IrregularSurface(
        zeros(Float32, M.nx),
        GhostPoint[],
        0,
        0,
        false,
        :direct_zero,
        M.dx, M.dz, M.pad
    )
end

"""
Compute bicubic interpolation weights (4x4 stencil).
"""
function _compute_interp_weights(x::Real, z::Real, dx::Real, dz::Real, 
                                  pad::Int, nx::Int, nz::Int)
    # Find grid position enclosing target point
    # Physical coordinates to grid index
    fi = x / dx + pad + 1  # Float index
    fj = z / dz + pad + 1
    
    # Center grid point (nearest)
    i0 = floor(Int, fi)
    j0 = floor(Int, fj)
    
    # Local coordinates (between 0 and 1)
    fx = fi - i0
    fz = fj - j0
    
    # 4x4 stencil indices (centered upper-left)
    i_indices = ntuple(k -> clamp(i0 - 1 + k, 1, nx), 4)
    j_indices = ntuple(k -> clamp(j0 - 1 + k, 1, nz), 4)
    
    # Compute bicubic interpolation weights
    wx = _cubic_weights(fx)
    wz = _cubic_weights(fz)
    
    # Flatten to 16 points
    interp_i = NTuple{16, Int}(ntuple(k -> i_indices[(k-1) % 4 + 1], 16))
    interp_j = NTuple{16, Int}(ntuple(k -> j_indices[(k-1) ÷ 4 + 1], 16))
    
    # Weights are outer product of wx and wz
    weights = Float32[]
    for jj in 1:4
        for ii in 1:4
            push!(weights, wx[ii] * wz[jj])
        end
    end
    interp_w = NTuple{16, Float32}(Tuple(weights))
    
    return interp_i, interp_j, interp_w
end

"""
Cubic interpolation kernel function (Catmull-Rom).
"""
function _cubic_weights(t::Real)
    t = Float32(t)
    t2 = t * t
    t3 = t2 * t
    
    # Catmull-Rom spline weights
    w0 = -0.5f0*t3 + t2 - 0.5f0*t
    w1 = 1.5f0*t3 - 2.5f0*t2 + 1.0f0
    w2 = -1.5f0*t3 + 2.0f0*t2 + 0.5f0*t
    w3 = 0.5f0*t3 - 0.5f0*t2
    
    return (w0, w1, w2, w3)
end

# ==============================================================================
# Source and Receiver Initialization
# ==============================================================================

"""
    setup_irregular_source(x, depth, wavelet, surface, medium; backend=CPU_BACKEND)

Setup source below irregular surface.

# Arguments
- `x`: Source x coordinate
- `depth`: Source depth (relative to surface at that location)
- `wavelet`: Source wavelet
- `surface`: IrregularSurface
- `medium`: Medium

# Returns
- `Source` struct (compatible with regular Source)
"""
function setup_irregular_source(x::Real, depth::Real, wavelet::Vector{Float32},
                                 surface::IrregularSurface, M::Medium;
                                 backend::AbstractBackend=CPU_BACKEND)
    dx, dz = M.dx, M.dz
    pad = M.pad
    
    # Compute surface height at this x location
    i_approx = round(Int, x / dx) + pad + 1
    i_approx = clamp(i_approx, 1, M.nx)
    
    # Get surface height (z_surface may be on GPU, need conversion)
    z_surf_array = surface.z_surface isa Array ? surface.z_surface : Array(surface.z_surface)
    z_surface_here = surface.enabled ? z_surf_array[i_approx] : 0.0f0
    
    # Absolute z coordinate = surface height + depth
    z_abs = z_surface_here + depth
    
    # Convert to grid indices
    i = round(Int, x / dx) + pad + 1
    j = round(Int, z_abs / dz) + pad + 1
    
    # Boundary check
    i = clamp(i, 1, M.nx)
    j = clamp(j, 1, M.nz)
    
    @info "Irregular source setup" x=x depth=depth z_surface=z_surface_here z_absolute=z_abs grid_i=i grid_j=j
    
    # Return standard Source struct
    return Source(i, j, to_device(wavelet, backend))
end

"""
    setup_irregular_receivers(x_positions, depths, surface, medium, nt; 
                              type=:vz, backend=CPU_BACKEND)

Setup receivers on/below irregular surface.

# Arguments
- `x_positions`: Receiver x coordinates array
- `depths`: Receiver depth array (relative to surface at each position, 0 = on surface)
- `surface`: IrregularSurface
- `medium`: Medium
- `nt`: Number of time steps
- `type`: Recording type (:vz, :vx, :p)

# Returns
- `Receivers` struct (compatible with regular Receivers)
"""
function setup_irregular_receivers(x::Vector{<:Real}, depth::Vector{<:Real},
                                    surface::IrregularSurface, M::Medium, nt::Int;
                                    type::Symbol=:vz, backend::AbstractBackend=CPU_BACKEND)
    
    @assert length(x) == length(depth) "x and depth must have same length"
    
    n_rec = length(x)
    dx, dz = M.dx, M.dz
    pad = M.pad
    
    i_rec = Vector{Int}(undef, n_rec)
    j_rec = Vector{Int}(undef, n_rec)
    
    # Get surface height array (ensure on CPU)
    z_surf_array = surface.z_surface isa Array ? surface.z_surface : Array(surface.z_surface)
    
    for r in 1:n_rec
        # Surface height at this x position
        i_approx = round(Int, x[r] / dx) + pad + 1
        i_approx = clamp(i_approx, 1, M.nx)
        
        z_surface_here = surface.enabled ? z_surf_array[i_approx] : 0.0f0
        
        # Absolute z coordinate
        z_abs = z_surface_here + depth[r]
        
        # Grid indices
        i_rec[r] = clamp(round(Int, x[r] / dx) + pad + 1, 1, M.nx)
        j_rec[r] = clamp(round(Int, z_abs / dz) + pad + 1, 1, M.nz)
    end
    
    # Create data buffer
    data = zeros(Float32, nt, n_rec)
    
    return Receivers(
        to_device(i_rec, backend),
        to_device(j_rec, backend),
        to_device(data, backend),
        type
    )
end

"""
Convenience function: distribute receivers uniformly on surface.
"""
function setup_surface_receivers(x_positions::Vector{<:Real}, surface::IrregularSurface,
                                  M::Medium, nt::Int; type::Symbol=:vz,
                                  backend::AbstractBackend=CPU_BACKEND)
    # All receivers at depth 0 (on surface)
    depths = zeros(Float32, length(x_positions))
    return setup_irregular_receivers(x_positions, depths, surface, M, nt; 
                                      type=type, backend=backend)
end
