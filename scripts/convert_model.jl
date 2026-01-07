#!/usr/bin/env julia
# ==============================================================================
# scripts/convert_model.jl
#
# Convert velocity models to JLD2 format
# Supports separate Vp, Vs, Rho files (common workflow)
#
# Usage:
#   # Three separate files (most common)
#   julia convert_model.jl --vp=vp.segy --vs=vs.segy --rho=rho.segy -o model.jld2 --dx=10
#   julia convert_model.jl --vp=vp.bin --vs=vs.bin --rho=rho.bin -o model.jld2 --dx=10 --nx=500 --nz=200
#
#   # Only Vp (Vs and Rho will be estimated)
#   julia convert_model.jl --vp=vp.segy -o model.jld2 --dx=10
#
#   # Single JLD2 file containing all (just copy/validate)
#   julia convert_model.jl --input=model_old.jld2 -o model.jld2
# ==============================================================================

import Pkg
Pkg.activate(dirname(@__DIR__))

using Fomo
using Fomo
using Printf

# Try to load optional packages at startup
const HAS_SEGYIO = try
    using SegyIO
    true
catch
    false
end

const HAS_MAT = try
    using MAT
    true
catch
    false
end

const HAS_NPZ = try
    using NPZ
    true
catch
    false
end

# ==============================================================================
# File Reading Functions
# ==============================================================================

"""
Read a single property file (auto-detect format from extension)
"""
function read_property_file(path::String; nx=nothing, nz=nothing, dtype=Float32, order=:column_major)
    ext = lowercase(splitext(path)[2])
    
    if ext == ".jld2"
        data = JLD2.load(path)
        # Try common key names
        for key in ["vp", "Vp", "VP", "vs", "Vs", "VS", "rho", "Rho", "density", "data"]
            if haskey(data, key)
                return Float32.(data[key])
            end
        end
        error("Cannot find data in JLD2 file: $path")
        
    elseif ext in [".segy", ".sgy"]
        if !HAS_SEGYIO
            error("SegyIO.jl not installed. Run: using Pkg; Pkg.add(\"SegyIO\")")
        end
        block = SegyIO.segy_read(path)
        return Float32.(block.data)
        
    elseif ext == ".bin"
        if nx === nothing || nz === nothing
            error("Binary files require --nx and --nz parameters")
        end
        data = zeros(dtype, nx * nz)
        open(path, "r") do io
            read!(io, data)
        end
        if order == :row_major
            return Float32.(permutedims(reshape(data, nz, nx)))
        else
            return Float32.(reshape(data, nx, nz))
        end
        
    elseif ext == ".npy"
        if !HAS_NPZ
            error("NPZ.jl not installed. Run: using Pkg; Pkg.add(\"NPZ\")")
        end
        return Float32.(NPZ.npzread(path))
        
    elseif ext == ".mat"
        if !HAS_MAT
            error("MAT.jl not installed. Run: using Pkg; Pkg.add(\"MAT\")")
        end
        data = MAT.matread(path)
        # Return first array found
        for (k, v) in data
            if v isa AbstractMatrix
                return Float32.(v)
            end
        end
        error("No matrix data found in MAT file: $path")
        
    else
        error("Unsupported format: $ext")
    end
end

# ==============================================================================
# Argument Parsing
# ==============================================================================

function parse_args(args)
    opts = Dict{Symbol, Any}(
        :vp => nothing,
        :vs => nothing,
        :rho => nothing,
        :output => nothing,
        :dx => nothing,
        :dz => nothing,
        :nx => nothing,
        :nz => nothing,
        :order => :column_major,
        :name => "model",
        :transpose => false,
        :no_transpose => false
    )
    
    i = 1
    while i <= length(args)
        arg = args[i]
        
        if startswith(arg, "--vp=")
            opts[:vp] = arg[6:end]
        elseif startswith(arg, "--vs=")
            opts[:vs] = arg[6:end]
        elseif startswith(arg, "--rho=")
            opts[:rho] = arg[7:end]
        elseif startswith(arg, "--dx=")
            opts[:dx] = parse(Float32, arg[6:end])
        elseif startswith(arg, "--dz=")
            opts[:dz] = parse(Float32, arg[6:end])
        elseif startswith(arg, "--nx=")
            opts[:nx] = parse(Int, arg[6:end])
        elseif startswith(arg, "--nz=")
            opts[:nz] = parse(Int, arg[6:end])
        elseif startswith(arg, "--name=")
            opts[:name] = arg[8:end]
        elseif startswith(arg, "--order=")
            opts[:order] = Symbol(arg[9:end])
        elseif arg == "--transpose"
            opts[:transpose] = true
        elseif arg == "--no-transpose"
            opts[:no_transpose] = true
        elseif arg == "-o" && i < length(args)
            i += 1
            opts[:output] = args[i]
        elseif startswith(arg, "-o")
            opts[:output] = arg[3:end]
        elseif !startswith(arg, "-")
            # Positional argument - assume output if not set
            if opts[:output] === nothing
                opts[:output] = arg
            end
        end
        
        i += 1
    end
    
    # Default dz = dx
    if opts[:dz] === nothing && opts[:dx] !== nothing
        opts[:dz] = opts[:dx]
    end
    
    return opts
end

# ==============================================================================
# Main
# ==============================================================================

function print_usage()
    println("""
Velocity Model Converter - Convert Vp/Vs/Rho files to JLD2 format

Usage:
  julia convert_model.jl --vp=<file> [--vs=<file>] [--rho=<file>] -o <output.jld2> [options]

Required:
  --vp=<file>       Path to Vp (P-wave velocity) file
  -o <output>       Output JLD2 file path

Optional:
  --vs=<file>       Path to Vs file (if not provided: Vs = Vp / 1.73)
  --rho=<file>      Path to density file (if not provided: Gardner relation)
  --dx=<float>      Grid spacing in X (required for SEG-Y, optional for others)
  --dz=<float>      Grid spacing in Z (default: same as dx)
  --nx=<int>        Grid points in X (required for binary files)
  --nz=<int>        Grid points in Z (required for binary files)
  --order=<order>   Memory order for binary: column_major (default) or row_major
  --name=<string>   Model name (default: "model")

Supported formats:
  .segy, .sgy       SEG-Y format (requires SegyIO.jl)
  .bin              Raw binary (requires --nx, --nz)
  .mat              MATLAB format (requires MAT.jl)
  .npy              NumPy format (requires NPZ.jl)
  .jld2             Julia JLD2 format

Examples:
  # Three SEG-Y files
  julia convert_model.jl --vp=vp.segy --vs=vs.segy --rho=rho.segy -o model.jld2 --dx=12.5

  # Three binary files
  julia convert_model.jl --vp=vp.bin --vs=vs.bin --rho=rho.bin -o model.jld2 \\
        --dx=10 --nx=500 --nz=200

  # Only Vp (Vs and Rho will be estimated)
  julia convert_model.jl --vp=velocity.segy -o model.jld2 --dx=10

  # Mixed formats
  julia convert_model.jl --vp=vp.segy --vs=vs.bin --rho=rho.mat -o model.jld2 \\
        --dx=10 --nx=500 --nz=200
""")
end

function main()
    if length(ARGS) < 1 || "--help" in ARGS || "-h" in ARGS
        print_usage()
        return
    end
    
    opts = parse_args(ARGS)
    
    # Validate required arguments
    if opts[:vp] === nothing
        println("Error: --vp is required")
        println("Run with --help for usage information")
        return
    end
    
    if opts[:output] === nothing
        println("Error: Output file (-o) is required")
        println("Run with --help for usage information")
        return
    end
    
    println("=" ^ 60)
    println("  Velocity Model Converter")
    println("=" ^ 60)
    println()
    
    # Read Vp
    println("Reading Vp: $(opts[:vp])")
    vp = read_property_file(opts[:vp]; 
                            nx=opts[:nx], nz=opts[:nz], order=opts[:order])
    nx, nz = size(vp)
    @printf("  Size: %d × %d\n", nx, nz)
    @printf("  Range: %.1f - %.1f m/s\n", minimum(vp), maximum(vp))
    
    # Read or estimate Vs
    if opts[:vs] !== nothing
        println("\nReading Vs: $(opts[:vs])")
        vs = read_property_file(opts[:vs]; 
                                nx=opts[:nx], nz=opts[:nz], order=opts[:order])
        @printf("  Range: %.1f - %.1f m/s\n", minimum(vs), maximum(vs))
    else
        println("\nEstimating Vs = Vp / 1.73")
        vs = vp ./ 1.73f0
        @printf("  Range: %.1f - %.1f m/s\n", minimum(vs), maximum(vs))
    end
    
    # Read or estimate Rho
    if opts[:rho] !== nothing
        println("\nReading Rho: $(opts[:rho])")
        rho = read_property_file(opts[:rho]; 
                                 nx=opts[:nx], nz=opts[:nz], order=opts[:order])
        @printf("  Range: %.1f - %.1f kg/m³\n", minimum(rho), maximum(rho))
    else
        println("\nEstimating Rho using Gardner relation: ρ = 310 × Vp^0.25")
        rho = 310.0f0 .* (vp .^ 0.25f0)
        @printf("  Range: %.1f - %.1f kg/m³\n", minimum(rho), maximum(rho))
    end
    
    # Validate sizes
    if size(vs) != (nx, nz)
        error("Vs size $(size(vs)) doesn't match Vp size ($nx, $nz)")
    end
    if size(rho) != (nx, nz)
        error("Rho size $(size(rho)) doesn't match Vp size ($nx, $nz)")
    end
    
    # Grid spacing
    dx = opts[:dx] !== nothing ? opts[:dx] : 1.0f0
    dz = opts[:dz] !== nothing ? opts[:dz] : dx
    
    if opts[:dx] === nothing
        @warn "dx not specified, using dx=1.0"
    end
    
    # Check dimensions - typically X (horizontal) > Z (depth)
    x_extent = (nx - 1) * dx
    z_extent = (nz - 1) * dz
    
    println("\n" * "-" ^ 60)
    println("  Dimension Check:")
    @printf("    X extent: %.1f m (%.2f km)\n", x_extent, x_extent/1000)
    @printf("    Z extent: %.1f m (%.2f km)\n", z_extent, z_extent/1000)
    
    need_transpose = false
    
    if opts[:transpose]
        # Force transpose
        need_transpose = true
        println("  → Transpose requested via --transpose")
    elseif opts[:no_transpose]
        # Force no transpose
        println("  → Keeping original orientation (--no-transpose)")
    elseif nz > nx * 2
        println()
        println("  ⚠️  WARNING: 维度可能是反的！")
        println("     通常地震模型: X (水平) > Z (深度)")
        println("     当前: X = $(round(x_extent/1000, digits=2))km, Z = $(round(z_extent/1000, digits=2))km")
        println()
        print("  是否转置? [y/N]: ")
        response = lowercase(strip(readline()))
        
        if response == "y" || response == "yes"
            need_transpose = true
        else
            println("  Keeping original orientation.")
        end
    else
        println("  ✓ Dimensions look correct (X > Z)")
    end
    
    if need_transpose
        println("  Transposing...")
        vp = permutedims(vp)
        vs = permutedims(vs)
        rho = permutedims(rho)
        nx, nz = size(vp)
        @printf("  New size: %d × %d\n", nx, nz)
        @printf("  New X extent: %.1f m (%.2f km)\n", (nx-1)*dx, (nx-1)*dx/1000)
        @printf("  New Z extent: %.1f m (%.2f km)\n", (nz-1)*dz, (nz-1)*dz/1000)
    end
    println("-" ^ 60)
    
    # Save
    println("\nSaving to: $(opts[:output])")
    model = VelocityModel(vp, vs, rho, dx, dz; name=opts[:name])
    save_model(opts[:output], model)
    
    println()
    println("=" ^ 60)
    println("  Conversion Complete!")
    println("=" ^ 60)
    model_info(model)
end

main()
