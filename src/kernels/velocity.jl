# ==============================================================================
# kernels/velocity.jl
#
# Velocity update kernel - dispatched by backend
# ==============================================================================

using LoopVectorization

"""
    update_velocity!(backend, W, M, a, params)

Update velocity fields (vx, vz) based on stress gradients.
Dispatches to CPU or GPU implementation based on backend.
"""
function update_velocity! end

# ==============================================================================
# CPU Implementation
# ==============================================================================

function update_velocity!(::CPUBackend, W::Wavefield, M::Medium, a::Vector{Float32}, p::SimParams)
    nx, nz = M.nx, M.nz
    dtx, dtz = p.dtx, p.dtz
    M_order = p.M
    
    vx, vz = W.vx, W.vz
    txx, tzz, txz = W.txx, W.tzz, W.txz
    rho_vx, rho_vz = M.rho_vx, M.rho_vz
    
    @tturbo for j in (M_order+1):(nz-M_order)
        for i in (M_order+1):(nx-M_order)
            dtxxdx, dtxzdz, dtxzdx, dtzzdz = 0.0f0, 0.0f0, 0.0f0, 0.0f0

            for l in 1:M_order
                dtxxdx += a[l] * (txx[i+l-1, j] - txx[i-l, j])
                dtxzdz += a[l] * (txz[i, j+l-1] - txz[i, j-l])
                dtxzdx += a[l] * (txz[i+l, j] - txz[i-l+1, j])
                dtzzdz += a[l] * (tzz[i, j+l] - tzz[i, j-l+1])
            end

            vx[i, j] += (dtx / rho_vx[i, j]) * dtxxdx + (dtz / rho_vx[i, j]) * dtxzdz
            vz[i, j] += (dtx / rho_vz[i, j]) * dtxzdx + (dtz / rho_vz[i, j]) * dtzzdz
        end
    end
    return nothing
end

# ==============================================================================
# CUDA Implementation
# ==============================================================================

function _update_velocity_kernel!(vx, vz, txx, tzz, txz, rho_vx, rho_vz, a, 
                                   nx, nz, dtx, dtz, M_order)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    
    if i > M_order && i <= nx - M_order && j > M_order && j <= nz - M_order
        dtxxdx, dtxzdz, dtxzdx, dtzzdz = 0.0f0, 0.0f0, 0.0f0, 0.0f0
        
        for l in 1:M_order
            dtxxdx += a[l] * (txx[i+l-1, j] - txx[i-l, j])
            dtxzdz += a[l] * (txz[i, j+l-1] - txz[i, j-l])
            dtxzdx += a[l] * (txz[i+l, j] - txz[i-l+1, j])
            dtzzdz += a[l] * (tzz[i, j+l] - tzz[i, j-l+1])
        end
        
        @inbounds vx[i, j] += (dtx / rho_vx[i, j]) * dtxxdx + (dtz / rho_vx[i, j]) * dtxzdz
        @inbounds vz[i, j] += (dtx / rho_vz[i, j]) * dtxzdx + (dtz / rho_vz[i, j]) * dtzzdz
    end
    return nothing
end

function update_velocity!(::CUDABackend, W::Wavefield, M::Medium, a::CuVector{Float32}, p::SimParams)
    nx, nz = M.nx, M.nz
    threads = (16, 16)
    blocks = (cld(nx, 16), cld(nz, 16))
    
    @cuda threads=threads blocks=blocks _update_velocity_kernel!(
        W.vx, W.vz, W.txx, W.tzz, W.txz, M.rho_vx, M.rho_vz, a,
        nx, nz, p.dtx, p.dtz, p.M
    )
    return nothing
end
