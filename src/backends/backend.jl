# ==============================================================================
# backends/backend.jl
#
# Backend abstraction for CPU/CUDA dispatch
# ==============================================================================

abstract type AbstractBackend end

struct CPUBackend <: AbstractBackend end
struct CUDABackend <: AbstractBackend end

# Singleton instances
const CPU_BACKEND = CPUBackend()
const CUDA_BACKEND = CUDABackend()

# ==============================================================================
# Backend Selection
# ==============================================================================

"""
    backend(type::Symbol) -> AbstractBackend

Select the computation backend for the simulation.

# Arguments
- `type::Symbol`: The backend type.
  - `:cpu` - Use CPU (threaded via LoopVectorization).
  - `:cuda` - Use NVIDIA GPU (via CUDA.jl).

# Returns
- `AbstractBackend`: A singleton instance of `CPUBackend` or `CUDABackend`.

# Throws
- `ErrorException`: If `:cuda` is requested but no GPU is available/functional.

# Example
```julia
# Auto-select
be = is_cuda_available() ? backend(:cuda) : backend(:cpu)

# Force CPU
be = backend(:cpu)
```
"""
function backend(type::Symbol)
    if type == :cpu
        return CPUBackend()
    elseif type == :cuda
        if !CUDA_AVAILABLE[]
            error("CUDA backend requested but no GPU available. Use backend(:cpu)")
        end
        return CUDABackend()
    else
        error("Unknown backend: $type. Use :cpu or :cuda")
    end
end

# ==============================================================================
# Device Transfer
# ==============================================================================

"""
    to_device(x, backend) -> device_array

Transfer array to the appropriate device.
"""
# CPU
to_device(x::AbstractArray, ::CPUBackend) = Array(x)
to_device(x::Number, ::CPUBackend) = x

# CUDA
to_device(x::AbstractArray, ::CUDABackend) = CuArray(x)
to_device(x::Number, ::CUDABackend) = x

# ==============================================================================
# Synchronization
# ==============================================================================

"""
    synchronize(backend)

Synchronize computation on the given backend.
"""
synchronize(::CPUBackend) = nothing
synchronize(::CUDABackend) = CUDA.synchronize()
