# ==============================================================================
# kernels/boundary.jl
#
# Boundary condition kernels - HABC and Free Surface
# Directly copied from original Kernels.jl with minimal changes
# ==============================================================================

using LoopVectorization

# ==============================================================================
# Boundary Strip Backup (Exact copy from original)
# ==============================================================================

"""
    copy_boundary_strip!(old, new, nbc, nx, nz, is_free_surface)

Backup boundary field values into `old` for HABC extrapolation.
Copied exactly from original Kernels.jl.
"""
function copy_boundary_strip!(old, new, nbc, nx, nz, is_free_surface)
    j_top = is_free_surface ? nbc + 1 : 1

    # Vertical strips (Left/Right)
    @inbounds for i in 1:nbc+2
        @views old[i, j_top:nz] .= new[i, j_top:nz]
    end
    @inbounds for i in (nx-nbc-1):nx
        @views old[i, j_top:nz] .= new[i, j_top:nz]
    end

    # Horizontal strips (Top/Bottom)
    @inbounds for j in j_top:nbc+2
        @views old[nbc+3:nx-nbc-2, j] .= new[nbc+3:nx-nbc-2, j]
    end
    @inbounds for j in (nz-nbc-1):nz
        @views old[nbc+3:nx-nbc-2, j] .= new[nbc+3:nx-nbc-2, j]
    end
end

"""
    backup_boundary!(backend, W, H, M)

Backup all wavefield boundary values.
"""
function backup_boundary!(::CPUBackend, W::Wavefield, H::HABCConfig, M::Medium)
    nx, nz = M.nx, M.nz
    nbc = H.nbc
    is_fs = M.is_free_surface
    
    copy_boundary_strip!(W.vx_old, W.vx, nbc, nx, nz, is_fs)
    copy_boundary_strip!(W.vz_old, W.vz, nbc, nx, nz, is_fs)
    copy_boundary_strip!(W.txx_old, W.txx, nbc, nx, nz, is_fs)
    copy_boundary_strip!(W.tzz_old, W.tzz, nbc, nx, nz, is_fs)
    copy_boundary_strip!(W.txz_old, W.txz, nbc, nx, nz, is_fs)
    return nothing
end

# ==============================================================================
# HABC Application (Exact copy from original Kernels.jl)
# ==============================================================================

"""
    apply_habc!(f, f_old, H, weights, nx, nz, is_free_surface)

Apply Higdon Absorbing Boundary Conditions (HABC) to field `f`.
Copied exactly from original Kernels.jl.
"""
function apply_habc!(f, f_old, H, weights, nx, nz, is_free_surface)
    nbc = H.nbc
    qx, qz, qt_x, qt_z, qxt = H.qx, H.qz, H.qt_x, H.qt_z, H.qxt
    j_start = is_free_surface ? nbc + 1 : 2

    # --- 1. Pure Edges (1D Absorption) ---

    # Left Edge
    @tturbo for j in j_start:(nz-nbc-1), i in 2:nbc+1
        sum_x = -qx * f[i+1, j] - qt_x * f_old[i, j] - qxt * f_old[i+1, j]
        f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * sum_x
    end

    # Right Edge
    @tturbo for j in j_start:(nz-nbc-1), i in (nx-nbc):nx-1
        sum_x = -qx * f[i-1, j] - qt_x * f_old[i, j] - qxt * f_old[i-1, j]
        f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * sum_x
    end

    # Bottom Edge
    @tturbo for j in (nz-nbc):nz-1, i in (nbc+2):(nx-nbc-1)
        sum_z = -qz * f[i, j-1] - qt_z * f_old[i, j] - qxt * f_old[i, j-1]
        f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * sum_z
    end

    # Top Edge (Skip if Free Surface)
    if !is_free_surface
        @tturbo for j in 2:nbc+1, i in (nbc+2):(nx-nbc-1)
            sum_z = -qz * f[i, j+1] - qt_z * f_old[i, j] - qxt * f_old[i, j+1]
            f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * sum_z
        end
    end

    # --- 2. Corner Coupling (Average of X and Z Absorption) ---

    # Left-Bottom Corner
    @tturbo for j in (nz-nbc):nz-1, i in 2:nbc+1
        sum_x = -qx * f[i+1, j] - qt_x * f_old[i, j] - qxt * f_old[i+1, j]
        sum_z = -qz * f[i, j-1] - qt_z * f_old[i, j] - qxt * f_old[i, j-1]
        f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * 0.5f0 * (sum_x + sum_z)
    end

    # Right-Bottom Corner
    @tturbo for j in (nz-nbc):nz-1, i in (nx-nbc):nx-1
        sum_x = -qx * f[i-1, j] - qt_x * f_old[i, j] - qxt * f_old[i-1, j]
        sum_z = -qz * f[i, j-1] - qt_z * f_old[i, j] - qxt * f_old[i, j-1]
        f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * 0.5f0 * (sum_x + sum_z)
    end

    if !is_free_surface
        # Left-Top Corner
        @tturbo for j in 2:nbc+1, i in 2:nbc+1
            sum_x = -qx * f[i+1, j] - qt_x * f_old[i, j] - qxt * f_old[i+1, j]
            sum_z = -qz * f[i, j+1] - qt_z * f_old[i, j] - qxt * f_old[i, j+1]
            f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * 0.5f0 * (sum_x + sum_z)
        end

        # Right-Top Corner
        @tturbo for j in 2:nbc+1, i in (nx-nbc):nx-1
            sum_x = -qx * f[i-1, j] - qt_x * f_old[i, j] - qxt * f_old[i-1, j]
            sum_z = -qz * f[i, j+1] - qt_z * f_old[i, j] - qxt * f_old[i, j+1]
            f[i, j] = weights[j, i] * f[i, j] + (1.0f0 - weights[j, i]) * 0.5f0 * (sum_x + sum_z)
        end
    end
end

"""
    apply_habc_velocity!(backend, W, H, M)

Apply HABC to velocity fields.
"""
function apply_habc_velocity!(::CPUBackend, W::Wavefield, H::HABCConfig, M::Medium)
    nx, nz = M.nx, M.nz
    apply_habc!(W.vx, W.vx_old, H, H.w_vx, nx, nz, M.is_free_surface)
    apply_habc!(W.vz, W.vz_old, H, H.w_vz, nx, nz, M.is_free_surface)
    return nothing
end

"""
    apply_habc_stress!(backend, W, H, M)

Apply HABC to stress fields.
"""
function apply_habc_stress!(::CPUBackend, W::Wavefield, H::HABCConfig, M::Medium)
    nx, nz = M.nx, M.nz
    apply_habc!(W.txx, W.txx_old, H, H.w_tau, nx, nz, M.is_free_surface)
    apply_habc!(W.tzz, W.tzz_old, H, H.w_tau, nx, nz, M.is_free_surface)
    apply_habc!(W.txz, W.txz_old, H, H.w_tau, nx, nz, M.is_free_surface)
    return nothing
end

# ==============================================================================
# Free Surface Condition (Exact copy from original Solver.jl)
# ==============================================================================

"""
    apply_free_surface!(backend, W, M)

Apply stress-free boundary condition at top (z=0).
"""
function apply_free_surface!(::CPUBackend, W::Wavefield, M::Medium)
    if !M.is_free_surface
        return nothing
    end
    
    nx = M.nx
    j_fs = M.pad + 1
    
    @tturbo for i in 1:nx, j in j_fs-5:j_fs
        W.tzz[i, j] = 0.0f0
        W.txz[i, j] = 0.0f0
    end
    return nothing
end

# ==============================================================================
# GPU Implementations
# ==============================================================================

function backup_boundary!(::CUDABackend, W::Wavefield, H::HABCConfig, M::Medium)
    nx, nz = M.nx, M.nz
    nbc = H.nbc
    is_fs = M.is_free_surface
    threads = (16, 16)
    blocks = (cld(nx, 16), cld(nz, 16))
    
    @cuda threads=threads blocks=blocks _backup_boundary_kernel!(
        W.vx, W.vx_old, W.vz, W.vz_old, W.txx, W.txx_old, 
        W.tzz, W.tzz_old, W.txz, W.txz_old,
        nx, nz, nbc, is_fs
    )
    return nothing
end

function _backup_boundary_kernel!(vx, vx_old, vz, vz_old, txx, txx_old, 
                                   tzz, tzz_old, txz, txz_old,
                                   nx, nz, nbc, is_free_surface)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    
    if i <= nx && j <= nz
        j_top = is_free_surface ? nbc + 1 : 1
        if j >= j_top
            is_strip = (i <= nbc + 2) || (i >= nx - nbc - 1) || 
                       (j <= nbc + 2) || (j >= nz - nbc - 1)
            if is_strip
                @inbounds begin
                    vx_old[i, j] = vx[i, j]
                    vz_old[i, j] = vz[i, j]
                    txx_old[i, j] = txx[i, j]
                    tzz_old[i, j] = tzz[i, j]
                    txz_old[i, j] = txz[i, j]
                end
            end
        end
    end
    return nothing
end

function apply_habc_velocity!(::CUDABackend, W::Wavefield, H::HABCConfig, M::Medium)
    _apply_habc_gpu!(W.vx, W.vx_old, H, H.w_vx, M)
    _apply_habc_gpu!(W.vz, W.vz_old, H, H.w_vz, M)
    return nothing
end

function apply_habc_stress!(::CUDABackend, W::Wavefield, H::HABCConfig, M::Medium)
    _apply_habc_gpu!(W.txx, W.txx_old, H, H.w_tau, M)
    _apply_habc_gpu!(W.tzz, W.tzz_old, H, H.w_tau, M)
    _apply_habc_gpu!(W.txz, W.txz_old, H, H.w_tau, M)
    return nothing
end

function _apply_habc_gpu!(f, f_old, H, weights, M)
    nx, nz = M.nx, M.nz
    threads = (16, 16)
    blocks = (cld(nx, 16), cld(nz, 16))
    
    @cuda threads=threads blocks=blocks _habc_kernel!(
        f, f_old, weights, H.qx, H.qz, H.qt_x, H.qt_z, H.qxt, 
        H.nbc, nx, nz, M.is_free_surface
    )
end

function _habc_kernel!(f, f_old, w, qx, qz, qt_x, qt_z, qxt, nbc, nx, nz, is_free_surface)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    
    j_start = is_free_surface ? nbc + 1 : 2
    
    if i > 1 && i < nx && j >= j_start && j < nz
        in_left = i <= nbc + 1
        in_right = i >= nx - nbc
        in_bottom = j >= nz - nbc
        in_top = !is_free_surface && j <= nbc + 1
        
        # Note: weights are (nz, nx), so index as w[j, i]
        wt = w[j, i]
        
        if in_left && !in_bottom && !in_top
            sum_x = -qx * f[i+1, j] - qt_x * f_old[i, j] - qxt * f_old[i+1, j]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * sum_x
        elseif in_right && !in_bottom && !in_top
            sum_x = -qx * f[i-1, j] - qt_x * f_old[i, j] - qxt * f_old[i-1, j]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * sum_x
        elseif in_bottom && !in_left && !in_right
            sum_z = -qz * f[i, j-1] - qt_z * f_old[i, j] - qxt * f_old[i, j-1]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * sum_z
        elseif in_top && !in_left && !in_right
            sum_z = -qz * f[i, j+1] - qt_z * f_old[i, j] - qxt * f_old[i, j+1]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * sum_z
        elseif in_left && in_bottom
            sum_x = -qx * f[i+1, j] - qt_x * f_old[i, j] - qxt * f_old[i+1, j]
            sum_z = -qz * f[i, j-1] - qt_z * f_old[i, j] - qxt * f_old[i, j-1]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * 0.5f0 * (sum_x + sum_z)
        elseif in_right && in_bottom
            sum_x = -qx * f[i-1, j] - qt_x * f_old[i, j] - qxt * f_old[i-1, j]
            sum_z = -qz * f[i, j-1] - qt_z * f_old[i, j] - qxt * f_old[i, j-1]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * 0.5f0 * (sum_x + sum_z)
        elseif in_left && in_top
            sum_x = -qx * f[i+1, j] - qt_x * f_old[i, j] - qxt * f_old[i+1, j]
            sum_z = -qz * f[i, j+1] - qt_z * f_old[i, j] - qxt * f_old[i, j+1]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * 0.5f0 * (sum_x + sum_z)
        elseif in_right && in_top
            sum_x = -qx * f[i-1, j] - qt_x * f_old[i, j] - qxt * f_old[i-1, j]
            sum_z = -qz * f[i, j+1] - qt_z * f_old[i, j] - qxt * f_old[i, j+1]
            @inbounds f[i, j] = wt * f[i, j] + (1.0f0 - wt) * 0.5f0 * (sum_x + sum_z)
        end
    end
    return nothing
end

function apply_free_surface!(::CUDABackend, W::Wavefield, M::Medium)
    if !M.is_free_surface
        return nothing
    end
    
    nx = M.nx
    j_fs = M.pad + 1
    threads = (256, 1)
    blocks = (cld(nx, 256), 6)
    
    @cuda threads=threads blocks=blocks _free_surface_kernel!(W.tzz, W.txz, nx, j_fs)
    return nothing
end

function _free_surface_kernel!(tzz, txz, nx, j_fs)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j_offset = blockIdx().y - 1
    j = j_fs - 5 + j_offset
    
    if i <= nx && j >= 1 && j <= j_fs
        @inbounds tzz[i, j] = 0.0f0
        @inbounds txz[i, j] = 0.0f0
    end
    return nothing
end
