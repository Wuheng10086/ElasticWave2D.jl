# ==============================================================================
# kernels/source_receiver.jl
#
# Source injection and receiver recording kernels
# ==============================================================================

# ==============================================================================
# Source Injection
# ==============================================================================

"""
    inject_source!(backend, W, source, k, dt)

Inject source wavelet at time step k.
"""
function inject_source! end

function inject_source!(::CPUBackend, W::Wavefield, src::Source, k::Int, dt::Float32)
    if k > length(src.wavelet)
        return nothing
    end
    
    wav = src.wavelet[k]
    i, j = src.i, src.j
    
    # Pressure source (explosion)
    W.txx[i, j] += wav
    W.tzz[i, j] += wav
    return nothing
end

function _inject_source_kernel!(txx, tzz, wavelet, i0, j0, k)
    if threadIdx().x == 1 && blockIdx().x == 1
        wav = wavelet[k]
        @inbounds begin
            txx[i0, j0] += wav
            tzz[i0, j0] += wav
        end
    end
    return nothing
end

function inject_source!(::CUDABackend, W::Wavefield, src::Source, k::Int, dt::Float32)
    if k > length(src.wavelet)
        return nothing
    end
    
    @cuda threads=1 blocks=1 _inject_source_kernel!(
        W.txx, W.tzz, src.wavelet, src.i, src.j, k
    )
    return nothing
end

# ==============================================================================
# Receiver Recording
# ==============================================================================

"""
    record_receivers!(backend, W, rec, k)

Record wavefield values at receiver locations for time step k.
"""
function record_receivers! end

function record_receivers!(::CPUBackend, W::Wavefield, rec::Receivers, k::Int)
    for r in 1:length(rec.i)
        ri, rj = rec.i[r], rec.j[r]
        if rec.type == :vz
            rec.data[k, r] = W.vz[ri, rj]
        elseif rec.type == :vx
            rec.data[k, r] = W.vx[ri, rj]
        elseif rec.type == :p
            rec.data[k, r] = (W.txx[ri, rj] + W.tzz[ri, rj]) * 0.5f0
        end
    end
    return nothing
end

function _record_kernel!(data, field, rec_i, rec_j, k, n_rec)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx <= n_rec
        ii, jj = rec_i[idx], rec_j[idx]
        @inbounds data[k, idx] = field[ii, jj]
    end
    return nothing
end

function record_receivers!(::CUDABackend, W::Wavefield, rec::Receivers, k::Int)
    n_rec = length(rec.i)
    
    # Select field based on type
    field = if rec.type == :vz
        W.vz
    elseif rec.type == :vx
        W.vx
    else  # :p - would need special handling, defaulting to vz
        W.vz
    end
    
    @cuda threads=256 blocks=cld(n_rec, 256) _record_kernel!(
        rec.data, field, rec.i, rec.j, k, n_rec
    )
    return nothing
end

# ==============================================================================
# Wavefield Reset
# ==============================================================================

"""
    reset!(backend, W)

Reset all wavefield components to zero.
"""
function reset!(::CPUBackend, W::Wavefield)
    for field in (W.vx, W.vz, W.txx, W.tzz, W.txz,
                  W.vx_old, W.vz_old, W.txx_old, W.tzz_old, W.txz_old)
        fill!(field, 0.0f0)
    end
    return nothing
end

function reset!(::CUDABackend, W::Wavefield)
    for field in (W.vx, W.vz, W.txx, W.tzz, W.txz,
                  W.vx_old, W.vz_old, W.txx_old, W.tzz_old, W.txz_old)
        fill!(field, 0.0f0)
    end
    return nothing
end

# Default dispatch (infer backend from wavefield type)
reset!(W::Wavefield{<:Array}) = reset!(CPU_BACKEND, W)
reset!(W::Wavefield{<:CuArray}) = reset!(CUDA_BACKEND, W)
