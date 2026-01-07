#!/usr/bin/env julia
# ==============================================================================
# scripts/load_results.jl
#
# Load and inspect simulation results (geometry + gathers)
# Useful for migration and other post-processing workflows
#
# Usage:
#   julia load_results.jl survey_geometry.jld2
#   julia load_results.jl survey_geometry.jld2 --plot
#   julia load_results.jl survey_geometry.jld2 shot_1.bin shot_2.bin --plot
# ==============================================================================

import Pkg
Pkg.activate(dirname(@__DIR__))

using Fomo
using Fomo
using Printf
using CairoMakie

# ==============================================================================
# Load Gather Data
# ==============================================================================

"""
    load_gather_data(path, nt, n_rec) -> Matrix{Float32}

Load gather from binary file.
Returns matrix of shape [nt × n_rec].
"""
function load_gather_data(path::String, nt::Int, n_rec::Int)
    data = zeros(Float32, nt * n_rec)
    open(path, "r") do io
        read!(io, data)
    end
    return reshape(data, nt, n_rec)
end

"""
    load_gather_data(path, geom::SurveyGeometry) -> Matrix{Float32}

Load gather using geometry for dimensions.
"""
function load_gather_data(path::String, geom::SurveyGeometry)
    return load_gather_data(path, geom.nt, geom.n_rec)
end

# ==============================================================================
# Print Functions
# ==============================================================================

function print_geometry_summary(geom::SurveyGeometry)
    println("=" ^ 60)
    println("  Single Shot Geometry - Shot #$(geom.shot_id)")
    println("=" ^ 60)
    println()
    println("  Source (actual discretized position):")
    @printf("    Position:  (%.2f, %.2f) m\n", geom.src_x, geom.src_z)
    @printf("    Grid idx:  (%d, %d)\n", geom.src_i, geom.src_j)
    println()
    println("  Receivers: $(geom.n_rec)")
    @printf("    X range:   %.2f - %.2f m\n", minimum(geom.rec_x), maximum(geom.rec_x))
    @printf("    Z range:   %.2f - %.2f m\n", minimum(geom.rec_z), maximum(geom.rec_z))
    @printf("    Spacing:   %.2f m (approx)\n", 
            geom.n_rec > 1 ? (maximum(geom.rec_x) - minimum(geom.rec_x)) / (geom.n_rec - 1) : 0.0)
    println()
    println("  Time:")
    @printf("    dt:        %.6f s (%.3f ms)\n", geom.dt, geom.dt * 1000)
    @printf("    nt:        %d samples\n", geom.nt)
    @printf("    T_max:     %.4f s\n", geom.t_max)
    println()
    println("  Grid (physical domain):")
    @printf("    Size:      %d × %d\n", geom.nx, geom.nz)
    @printf("    Spacing:   dx=%.2f m, dz=%.2f m\n", geom.dx, geom.dz)
    @printf("    Extent:    %.2f × %.2f m\n", (geom.nx-1)*geom.dx, (geom.nz-1)*geom.dz)
    println("=" ^ 60)
end

function print_geometry_summary(mg::MultiShotGeometry)
    println("=" ^ 60)
    println("  Multi-Shot Survey Geometry")
    println("=" ^ 60)
    println()
    println("  Survey Overview:")
    @printf("    Number of shots: %d\n", mg.n_shots)
    
    src_x = [s.src_x for s in mg.shots]
    src_z = [s.src_z for s in mg.shots]
    @printf("    Source X range:  %.2f - %.2f m\n", minimum(src_x), maximum(src_x))
    @printf("    Source Z range:  %.2f - %.2f m\n", minimum(src_z), maximum(src_z))
    if mg.n_shots > 1
        @printf("    Shot spacing:    %.2f m (approx)\n", 
                (maximum(src_x) - minimum(src_x)) / (mg.n_shots - 1))
    end
    println()
    
    s1 = mg.shots[1]
    println("  Receivers (per shot): $(s1.n_rec)")
    @printf("    X range:   %.2f - %.2f m\n", minimum(s1.rec_x), maximum(s1.rec_x))
    @printf("    Z range:   %.2f - %.2f m\n", minimum(s1.rec_z), maximum(s1.rec_z))
    println()
    
    println("  Time:")
    @printf("    dt:        %.6f s (%.3f ms)\n", mg.dt, mg.dt * 1000)
    @printf("    nt:        %d samples\n", mg.nt)
    @printf("    T_max:     %.4f s\n", mg.dt * mg.nt)
    println()
    
    println("  Grid (physical domain):")
    @printf("    Size:      %d × %d\n", mg.nx, mg.nz)
    @printf("    Spacing:   dx=%.2f m, dz=%.2f m\n", mg.dx, mg.dz)
    println()
    
    println("  Shot List:")
    println("  " * "-" ^ 56)
    @printf("  %6s  %10s  %10s  %6s  %6s  %6s\n", 
            "ID", "src_x(m)", "src_z(m)", "src_i", "src_j", "n_rec")
    println("  " * "-" ^ 56)
    for s in mg.shots
        @printf("  %6d  %10.2f  %10.2f  %6d  %6d  %6d\n",
                s.shot_id, s.src_x, s.src_z, s.src_i, s.src_j, s.n_rec)
    end
    println("  " * "-" ^ 56)
    println("=" ^ 60)
end

function print_gather_summary(gather::Matrix{Float32}, name::String="Gather")
    nt, n_rec = size(gather)
    println()
    println("  $name:")
    @printf("    Shape:     [%d × %d] (nt × n_rec)\n", nt, n_rec)
    @printf("    Min:       %.6e\n", minimum(gather))
    @printf("    Max:       %.6e\n", maximum(gather))
    @printf("    RMS:       %.6e\n", sqrt(sum(gather.^2) / length(gather)))
end

# ==============================================================================
# Plot Functions
# ==============================================================================

function plot_gather(gather::Matrix{Float32}, geom::SurveyGeometry; 
                     output::String="gather_shot_$(geom.shot_id).png",
                     clip_percentile::Float64=99.0)
    nt, n_rec = size(gather)
    t_axis = range(0, geom.t_max, length=nt)
    
    # Clip for visualization
    vmax = percentile(abs.(gather[:]), clip_percentile)
    
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1,1],
        xlabel = "Receiver Index",
        ylabel = "Time (s)",
        title = "Shot Gather #$(geom.shot_id) - Source at ($(geom.src_x), $(geom.src_z)) m"
    )
    
    hm = heatmap!(ax, 1:n_rec, t_axis, gather',
                  colormap=:seismic, colorrange=(-vmax, vmax))
    ax.yreversed = true
    
    Colorbar(fig[1,2], hm, label="Amplitude")
    
    save(output, fig)
    println("  Plot saved: $output")
    return fig
end

function plot_geometry(geom::MultiShotGeometry; output::String="survey_geometry.png")
    fig = Figure(size=(900, 400))
    ax = Axis(fig[1,1],
        xlabel = "X (m)",
        ylabel = "Z (m)",
        title = "Survey Geometry ($(geom.n_shots) shots, $(geom.shots[1].n_rec) receivers/shot)",
        aspect = DataAspect()
    )
    
    # Plot receivers for first shot
    s1 = geom.shots[1]
    scatter!(ax, s1.rec_x, s1.rec_z, 
             marker=:dtriangle, markersize=8, color=:blue, label="Receivers")
    
    # Plot all sources
    src_x = [s.src_x for s in geom.shots]
    src_z = [s.src_z for s in geom.shots]
    scatter!(ax, src_x, src_z,
             marker=:star5, markersize=12, color=:red, label="Sources")
    
    ax.yreversed = true
    axislegend(ax, position=:rt)
    
    save(output, fig)
    println("  Plot saved: $output")
    return fig
end

# Helper
function percentile(arr, p)
    sorted = sort(arr[:])
    idx = clamp(ceil(Int, p/100 * length(sorted)), 1, length(sorted))
    return sorted[idx]
end

# ==============================================================================
# Export for External Use
# ==============================================================================

"""
    load_simulation_results(geometry_path; gather_paths=String[]) -> (geom, gathers)

Load geometry and optionally gather data.

# Returns
- `geom`: SurveyGeometry or MultiShotGeometry
- `gathers`: Vector of gather matrices (empty if no gather_paths provided)

# Example
```julia
include("scripts/load_results.jl")

# Load geometry only
geom, _ = load_simulation_results("survey.jld2")

# Load geometry + gathers
geom, gathers = load_simulation_results("survey.jld2"; 
    gather_paths=["shot_1.bin", "shot_2.bin"])

# Access data
geom.shots[1].src_x      # Source X position
geom.shots[1].rec_x      # Receiver X positions
gathers[1]               # First shot gather [nt × n_rec]
```
"""
function load_simulation_results(geometry_path::String; gather_paths::Vector{String}=String[])
    # Load geometry
    geom = load_geometry(geometry_path)
    
    # Load gathers if provided
    gathers = Matrix{Float32}[]
    if !isempty(gather_paths)
        for (i, gpath) in enumerate(gather_paths)
            if geom isa MultiShotGeometry
                g = load_gather_data(gpath, geom.shots[i])
            else
                g = load_gather_data(gpath, geom)
            end
            push!(gathers, g)
        end
    end
    
    return geom, gathers
end

# ==============================================================================
# Main
# ==============================================================================

function main()
    if length(ARGS) < 1
        println("""
Usage: julia load_results.jl <geometry_file> [gather_files...] [options]

Examples:
  # View geometry only
  julia load_results.jl survey_geometry.jld2

  # View geometry + plot
  julia load_results.jl survey_geometry.jld2 --plot

  # Load geometry + gather files
  julia load_results.jl survey_geometry.jld2 shot_1.bin shot_2.bin

  # Load all and plot
  julia load_results.jl survey_geometry.jld2 shot_1.bin --plot

Options:
  --plot    Generate plots for geometry and gathers

Output files (with --plot):
  - survey_geometry.png    Survey layout
  - gather_shot_N.png      Shot gather visualization

For programmatic use:
  include("scripts/load_results.jl")
  geom, gathers = load_simulation_results("survey.jld2"; 
      gather_paths=["shot_1.bin"])
""")
        return
    end
    
    # Parse arguments
    do_plot = "--plot" in ARGS
    args = filter(a -> !startswith(a, "--"), ARGS)
    
    geometry_path = args[1]
    gather_paths = length(args) > 1 ? args[2:end] : String[]
    
    # Load geometry
    println("\nLoading geometry: $geometry_path")
    geom = load_geometry(geometry_path)
    print_geometry_summary(geom)
    
    # Load gathers
    gathers = Matrix{Float32}[]
    if !isempty(gather_paths)
        println("\nLoading gathers:")
        for (i, gpath) in enumerate(gather_paths)
            println("  [$i] $gpath")
            shot_geom = geom isa MultiShotGeometry ? geom.shots[i] : geom
            g = load_gather_data(gpath, shot_geom)
            push!(gathers, g)
            print_gather_summary(g, "Shot #$(shot_geom.shot_id)")
        end
    end
    
    # Plot if requested
    if do_plot
        println("\nGenerating plots...")
        
        if geom isa MultiShotGeometry
            plot_geometry(geom)
        end
        
        for (i, g) in enumerate(gathers)
            shot_geom = geom isa MultiShotGeometry ? geom.shots[i] : geom
            plot_gather(g, shot_geom)
        end
    end
    
    println("\nDone!")
    
    # Return for interactive use
    return geom, gathers
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
