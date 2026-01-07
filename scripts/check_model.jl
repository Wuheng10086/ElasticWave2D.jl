#!/usr/bin/env julia
# ==============================================================================
# scripts/check_model.jl
#
# Check and fix model dimensions
# ==============================================================================

import Pkg
Pkg.activate(dirname(@__DIR__))

using Fomo
using Fomo

function main()
    if length(ARGS) < 1
        println("Usage: julia check_model.jl <model.jld2> [--fix]")
        return
    end
    
    path = ARGS[1]
    do_fix = "--fix" in ARGS
    
    println("=" ^ 60)
    println("  Model Dimension Check")
    println("=" ^ 60)
    
    model = load_model(path)
    
    println("\nCurrent dimensions:")
    println("  nx (X方向) = $(model.nx)  →  $(model.nx * model.dx) m")
    println("  nz (Z方向) = $(model.nz)  →  $(model.nz * model.dz) m")
    
    # Check if dimensions look swapped
    # Typically X (horizontal) > Z (depth) for seismic models
    if model.nz > model.nx * 2
        println("\n⚠️  WARNING: 维度可能是反的！")
        println("  通常地震模型: X (水平) > Z (深度)")
        println("  当前模型:     X = $(model.nx * model.dx)m, Z = $(model.nz * model.dz)m")
        println("\n  Marmousi2 标准尺寸: X ≈ 17km, Z ≈ 3.5km")
        println("  你的模型尺寸:       X ≈ $(round(model.nx * model.dx / 1000, digits=1))km, Z ≈ $(round(model.nz * model.dz / 1000, digits=1))km")
        
        if do_fix
            println("\n🔧 Fixing: 转置模型...")
            
            # Transpose all fields
            vp_t = permutedims(model.vp)
            vs_t = permutedims(model.vs)
            rho_t = permutedims(model.rho)
            
            # Create new model with swapped dimensions
            fixed_model = VelocityModel(vp_t, vs_t, rho_t, model.dx, model.dz; 
                                        name=model.name * "_fixed")
            
            # Save
            out_path = replace(path, ".jld2" => "_fixed.jld2")
            save_model(out_path, fixed_model)
            
            println("\n✅ Fixed model saved to: $out_path")
            println("\nNew dimensions:")
            println("  nx (X方向) = $(fixed_model.nx)  →  $(fixed_model.nx * fixed_model.dx) m")
            println("  nz (Z方向) = $(fixed_model.nz)  →  $(fixed_model.nz * fixed_model.dz) m")
        else
            println("\n💡 To fix, run:")
            println("   julia scripts/check_model.jl $path --fix")
        end
    else
        println("\n✅ 维度看起来正确！")
        println("   X (水平) = $(round(model.nx * model.dx / 1000, digits=1))km")
        println("   Z (深度) = $(round(model.nz * model.dz / 1000, digits=1))km")
    end
    
    println("\n" * "=" ^ 60)
end

main()
