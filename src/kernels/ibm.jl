# ==============================================================================
# kernels/immersed_boundary.jl
#
# Immersed Boundary Method (IBM) kernels for irregular free surface
#
# Supports two methods:
# 1. :direct_zero - Simple and stable, sets stress to zero above surface
# 2. :mirror      - Higher accuracy mirror interpolation (with safety checks)
#
# Elastic wave free surface boundary condition:
# - τzz = 0 (normal stress is zero)
# - τxz = 0 (shear stress is zero)
# ==============================================================================

"""
Available IBM methods:
- `:direct_zero` - Simple direct zero method (stable, lower accuracy)
- `:mirror`      - Mirror interpolation method (higher accuracy, with safety checks)
"""
const IBM_METHODS = (:direct_zero, :mirror)

"""
    apply_irregular_free_surface!(backend, W, surface)

Apply irregular free surface boundary condition.

The method used depends on the surface configuration.
"""
function apply_irregular_free_surface! end

# ==============================================================================
# GPU Structure
# ==============================================================================

"""
    IrregularSurfaceGPU

Irregular surface data on GPU.
"""
struct IrregularSurfaceGPU
    z_surface::CuVector{Float32}  # Surface elevation array
    gi::CuVector{Int32}           # Ghost point i indices
    gj::CuVector{Int32}           # Ghost point j indices
    interp_i::CuMatrix{Int32}     # Interpolation i indices (16 × n_ghost)
    interp_j::CuMatrix{Int32}     # Interpolation j indices
    interp_w::CuMatrix{Float32}   # Interpolation weights
    safe_mask::CuVector{Bool}     # Mask for safe interpolation points
    n_ghost::Int
    n_iter::Int
    enabled::Bool
    method::Symbol                # :direct_zero or :mirror
    dz::Float32
    pad::Int32
end

"""
    to_gpu(surface::IrregularSurface)

Convert IrregularSurface to GPU-optimized format.
"""
function to_gpu(surface::IrregularSurface)
    if !surface.enabled || surface.n_ghost == 0
        return IrregularSurfaceGPU(
            CUDA.zeros(Float32, 0),
            CUDA.zeros(Int32, 0),
            CUDA.zeros(Int32, 0),
            CUDA.zeros(Int32, 0, 0),
            CUDA.zeros(Int32, 0, 0),
            CUDA.zeros(Float32, 0, 0),
            CUDA.zeros(Bool, 0),
            0, 0, false, :direct_zero, 0.0f0, Int32(0)
        )
    end
    
    n_ghost = surface.n_ghost
    
    # Get z_surface on GPU
    z_surf_cpu = surface.z_surface isa Array ? surface.z_surface : Array(surface.z_surface)
    z_surf_gpu = CuArray(Float32.(z_surf_cpu))
    
    gi = Int32[gp.gi for gp in surface.ghost_points]
    gj = Int32[gp.gj for gp in surface.ghost_points]
    
    # Build interpolation matrices and safety mask
    interp_i_host = zeros(Int32, 16, n_ghost)
    interp_j_host = zeros(Int32, 16, n_ghost)
    interp_w_host = zeros(Float32, 16, n_ghost)
    safe_mask_host = zeros(Bool, n_ghost)
    
    for (idx, gp) in enumerate(surface.ghost_points)
        # Check if all interpolation points are below the ghost point (safe)
        # A point is safe if its interpolation stencil doesn't include ghost points
        min_interp_j = minimum(gp.interp_j)
        is_safe = min_interp_j > gp.gj + 1  # Stencil must be clearly below ghost layer
        safe_mask_host[idx] = is_safe
        
        for k in 1:16
            interp_i_host[k, idx] = gp.interp_i[k]
            interp_j_host[k, idx] = gp.interp_j[k]
            interp_w_host[k, idx] = gp.interp_w[k]
        end
    end
    
    n_safe = sum(safe_mask_host)
    @debug "IBM GPU: $n_safe / $n_ghost ghost points are safe for mirror interpolation"
    
    return IrregularSurfaceGPU(
        z_surf_gpu,
        CuArray(gi), CuArray(gj),
        CuArray(interp_i_host),
        CuArray(interp_j_host),
        CuArray(interp_w_host),
        CuArray(safe_mask_host),
        n_ghost,
        surface.n_iter,
        true,
        surface.method,
        Float32(surface.dz),
        Int32(surface.pad)
    )
end

# ==============================================================================
# CPU Implementation - Direct Zero Method
# ==============================================================================

function _apply_direct_zero_cpu!(W::Wavefield, z_surf::Vector{Float32}, dz::Float32, pad::Int)
    nx = length(z_surf)
    nz = size(W.tzz, 2)
    
    @inbounds for i in 1:nx
        z_here = z_surf[i]
        j_surf = round(Int, z_here / dz) + pad + 1
        j_surf = min(j_surf, nz)
        
        for j in 1:j_surf
            W.tzz[i, j] = 0.0f0
            W.txz[i, j] = 0.0f0
        end
    end
end

# ==============================================================================
# CPU Implementation - Mirror Method (with safety checks)
# ==============================================================================

"""
Interpolate field value using precomputed weights.
"""
@inline function _interpolate_field(field::AbstractMatrix, gp::GhostPoint)
    val = 0.0f0
    @inbounds for k in 1:16
        val += gp.interp_w[k] * field[gp.interp_i[k], gp.interp_j[k]]
    end
    return val
end

"""
Check if ghost point interpolation is safe (stencil doesn't include ghost points).
"""
@inline function _is_safe_interpolation(gp::GhostPoint)
    min_j = minimum(gp.interp_j)
    return min_j > gp.gj + 1
end

function _apply_mirror_cpu!(W::Wavefield, surface::IrregularSurface)
    ghost_points = surface.ghost_points
    n_iter = surface.n_iter
    
    # Only iterate once for stability (higher iterations can cause feedback)
    actual_iter = min(n_iter, 3)
    
    for iter in 1:actual_iter
        for gp in ghost_points
            if _is_safe_interpolation(gp)
                # Safe: use mirror interpolation
                tzz_mirror = _interpolate_field(W.tzz, gp)
                txz_mirror = _interpolate_field(W.txz, gp)
                
                @inbounds W.tzz[gp.gi, gp.gj] = -tzz_mirror
                @inbounds W.txz[gp.gi, gp.gj] = -txz_mirror
            else
                # Unsafe: use direct zero
                @inbounds W.tzz[gp.gi, gp.gj] = 0.0f0
                @inbounds W.txz[gp.gi, gp.gj] = 0.0f0
            end
        end
    end
end

# ==============================================================================
# CPU Dispatch
# ==============================================================================

function apply_irregular_free_surface!(::CPUBackend, W::Wavefield, surface::IrregularSurface)
    if !surface.enabled
        return nothing
    end
    
    z_surf = surface.z_surface isa Array ? Vector{Float32}(surface.z_surface) : 
             Vector{Float32}(Array(surface.z_surface))
    
    if surface.method == :direct_zero
        _apply_direct_zero_cpu!(W, z_surf, surface.dz, surface.pad)
    elseif surface.method == :mirror
        _apply_mirror_cpu!(W, surface)
    else
        error("Unknown IBM method: $(surface.method)")
    end
    
    return nothing
end

# ==============================================================================
# CUDA Kernels
# ==============================================================================

# Direct zero kernel
function _direct_zero_kernel!(tzz, txz, z_surface, dz, pad, nx, nz)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    
    if i <= nx
        @inbounds z_here = z_surface[i]
        j_surf = round(Int32, z_here / dz) + pad + Int32(1)
        j_surf = min(j_surf, nz)
        
        @inbounds for j in Int32(1):j_surf
            tzz[i, j] = 0.0f0
            txz[i, j] = 0.0f0
        end
    end
    
    return nothing
end

# Mirror kernel with safety mask
function _mirror_kernel!(tzz, txz, gi, gj, interp_i, interp_j, interp_w, safe_mask, n_ghost)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    
    if idx <= n_ghost
        @inbounds begin
            i = gi[idx]
            j = gj[idx]
            
            if safe_mask[idx]
                # Safe: use mirror interpolation
                tzz_mirror = 0.0f0
                txz_mirror = 0.0f0
                
                for k in 1:16
                    ii = interp_i[k, idx]
                    jj = interp_j[k, idx]
                    ww = interp_w[k, idx]
                    
                    tzz_mirror += ww * tzz[ii, jj]
                    txz_mirror += ww * txz[ii, jj]
                end
                
                tzz[i, j] = -tzz_mirror
                txz[i, j] = -txz_mirror
            else
                # Unsafe: direct zero
                tzz[i, j] = 0.0f0
                txz[i, j] = 0.0f0
            end
        end
    end
    
    return nothing
end

# ==============================================================================
# CUDA Dispatch - IrregularSurface
# ==============================================================================

function apply_irregular_free_surface!(::CUDABackend, W::Wavefield, surface::IrregularSurface)
    if !surface.enabled
        return nothing
    end
    
    z_surf = surface.z_surface isa CuArray ? surface.z_surface : 
             CuArray(Float32.(surface.z_surface isa Array ? surface.z_surface : Array(surface.z_surface)))
    
    if surface.method == :direct_zero
        nx = length(z_surf)
        nz = size(W.tzz, 2)
        
        threads = 256
        blocks = cld(nx, threads)
        
        @cuda threads=threads blocks=blocks _direct_zero_kernel!(
            W.tzz, W.txz, z_surf,
            Float32(surface.dz), Int32(surface.pad), Int32(nx), Int32(nz)
        )
    elseif surface.method == :mirror
        # For non-GPU-optimized surface, fall back to direct zero on GPU
        # (mirror requires precomputed safe_mask which is only in IrregularSurfaceGPU)
        nx = length(z_surf)
        nz = size(W.tzz, 2)
        
        threads = 256
        blocks = cld(nx, threads)
        
        @cuda threads=threads blocks=blocks _direct_zero_kernel!(
            W.tzz, W.txz, z_surf,
            Float32(surface.dz), Int32(surface.pad), Int32(nx), Int32(nz)
        )
    else
        error("Unknown IBM method: $(surface.method)")
    end
    
    return nothing
end

# ==============================================================================
# CUDA Dispatch - IrregularSurfaceGPU
# ==============================================================================

function apply_irregular_free_surface!(::CUDABackend, W::Wavefield, surface::IrregularSurfaceGPU)
    if !surface.enabled
        return nothing
    end
    
    if surface.method == :direct_zero
        nx = length(surface.z_surface)
        nz = size(W.tzz, 2)
        
        threads = 256
        blocks = cld(nx, threads)
        
        @cuda threads=threads blocks=blocks _direct_zero_kernel!(
            W.tzz, W.txz, surface.z_surface,
            surface.dz, surface.pad, Int32(nx), Int32(nz)
        )
    elseif surface.method == :mirror
        n_ghost = surface.n_ghost
        n_iter = min(surface.n_iter, 3)  # Limit iterations for stability
        
        threads = 256
        blocks = cld(n_ghost, threads)
        
        for iter in 1:n_iter
            @cuda threads=threads blocks=blocks _mirror_kernel!(
                W.tzz, W.txz,
                surface.gi, surface.gj,
                surface.interp_i, surface.interp_j, surface.interp_w,
                surface.safe_mask, n_ghost
            )
        end
    else
        error("Unknown IBM method: $(surface.method)")
    end
    
    return nothing
end

# ==============================================================================
# Helper functions
# ==============================================================================

"""
Check if irregular surface is enabled.
"""
is_irregular_surface_enabled(s::IrregularSurface) = s.enabled
is_irregular_surface_enabled(s::IrregularSurfaceGPU) = s.enabled

"""
Get IBM method being used.
"""
get_ibm_method(s::IrregularSurface) = s.method
get_ibm_method(s::IrregularSurfaceGPU) = s.method
