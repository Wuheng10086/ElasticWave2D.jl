# ==============================================================================
# backends/backend.jl
#
# Backend abstraction for CPU/GPU dispatch
# ==============================================================================

"""
Abstract backend type. All compute operations are dispatched based on backend.
"""
abstract type AbstractBackend end

"""
CPU backend - uses LoopVectorization for SIMD.
"""
struct CPUBackend <: AbstractBackend end

"""
CUDA backend - uses CUDA.jl kernels.
"""
struct CUDABackend <: AbstractBackend end

# Singleton instances
const CPU_BACKEND = CPUBackend()
const CUDA_BACKEND = CUDABackend()

"""
    backend(name::Symbol) -> AbstractBackend
    backend(name::String) -> AbstractBackend

Get backend by name.

# Examples
```julia
b = backend(:cpu)
b = backend(:cuda)
b = backend("gpu")
```
"""
function backend(name::Symbol)
    if name in (:cpu, :CPU)
        return CPU_BACKEND
    elseif name in (:cuda, :CUDA, :gpu, :GPU)
        if !CUDA_AVAILABLE[]
            @warn "CUDA not available, falling back to CPU"
            return CPU_BACKEND
        end
        return CUDA_BACKEND
    else
        error("Unknown backend: $name. Use :cpu or :cuda")
    end
end

backend(name::String) = backend(Symbol(lowercase(name)))

"""
    ArrayType(::AbstractBackend)

Get the array type for a backend.
"""
ArrayType(::CPUBackend) = Array
ArrayType(::CUDABackend) = CuArray

"""
    to_device(data, backend::AbstractBackend)

Move data to the specified backend's device.
"""
to_device(data::Array, ::CPUBackend) = data
to_device(data::Array, ::CUDABackend) = CuArray(data)
to_device(data::CuArray, ::CPUBackend) = Array(data)
to_device(data::CuArray, ::CUDABackend) = data

"""
    synchronize(::AbstractBackend)

Synchronize the backend (no-op for CPU, CUDA.synchronize() for GPU).
"""
synchronize(::CPUBackend) = nothing
synchronize(::CUDABackend) = CUDA.synchronize()
