# ==============================================================================
# kernels/stress.jl
#
# Stress update kernel - dispatched by backend
# ==============================================================================

"""
    update_stress!(backend, W, M, a, params)

Update stress fields (txx, tzz, txz) based on velocity gradients.
"""
function update_stress! end

# ==============================================================================
# CPU Implementation
# ==============================================================================

function update_stress!(::CPUBackend, W::Wavefield, M::Medium, a::Vector{Float32}, p::SimParams)
    nx, nz = M.nx, M.nz
    dtx, dtz = p.dtx, p.dtz
    M_order = p.M
    
    vx, vz = W.vx, W.vz
    txx, tzz, txz = W.txx, W.tzz, W.txz
    lam, mu_txx, mu_txz = M.lam, M.mu_txx, M.mu_txz
    
    @tturbo for j in (M_order+1):(nz-M_order)
        for i in (M_order+1):(nx-M_order)
            dvxdx, dvzdz, dvxdz, dvzdx = 0.0f0, 0.0f0, 0.0f0, 0.0f0

            for l in 1:M_order
                dvxdx += a[l] * (vx[i+l, j] - vx[i-l+1, j])
                dvzdz += a[l] * (vz[i, j+l-1] - vz[i, j-l])
                dvxdz += a[l] * (vx[i, j+l] - vx[i, j-l+1])
                dvzdx += a[l] * (vz[i+l-1, j] - vz[i-l, j])
            end

            l_val = lam[i, j]
            m_val = mu_txx[i, j]

            txx[i, j] += (l_val + 2.0f0 * m_val) * (dvxdx * dtx) + l_val * (dvzdz * dtz)
            tzz[i, j] += l_val * (dvxdx * dtx) + (l_val + 2.0f0 * m_val) * (dvzdz * dtz)
            txz[i, j] += mu_txz[i, j] * (dvxdz * dtz + dvzdx * dtx)
        end
    end
    return nothing
end

# ==============================================================================
# CUDA Implementation
# ==============================================================================

function _update_stress_kernel!(txx, tzz, txz, vx, vz, lam, mu_txx, mu_txz, a,
                                 nx, nz, dtx, dtz, M_order)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    
    if i > M_order && i <= nx - M_order && j > M_order && j <= nz - M_order
        dvxdx, dvzdz, dvxdz, dvzdx = 0.0f0, 0.0f0, 0.0f0, 0.0f0
        
        for l in 1:M_order
            dvxdx += a[l] * (vx[i+l, j] - vx[i-l+1, j])
            dvzdz += a[l] * (vz[i, j+l-1] - vz[i, j-l])
            dvxdz += a[l] * (vx[i, j+l] - vx[i, j-l+1])
            dvzdx += a[l] * (vz[i+l-1, j] - vz[i-l, j])
        end
        
        l_val, m_val = lam[i, j], mu_txx[i, j]
        
        @inbounds txx[i, j] += (l_val + 2.0f0 * m_val) * (dvxdx * dtx) + l_val * (dvzdz * dtz)
        @inbounds tzz[i, j] += l_val * (dvxdx * dtx) + (l_val + 2.0f0 * m_val) * (dvzdz * dtz)
        @inbounds txz[i, j] += mu_txz[i, j] * (dvxdz * dtz + dvzdx * dtx)
    end
    return nothing
end

function update_stress!(::CUDABackend, W::Wavefield, M::Medium, a::CuVector{Float32}, p::SimParams)
    nx, nz = M.nx, M.nz
    threads = (16, 16)
    blocks = (cld(nx, 16), cld(nz, 16))
    
    @cuda threads=threads blocks=blocks _update_stress_kernel!(
        W.txx, W.tzz, W.txz, W.vx, W.vz, M.lam, M.mu_txx, M.mu_txz, a,
        nx, nz, p.dtx, p.dtz, p.M
    )
    return nothing
end
