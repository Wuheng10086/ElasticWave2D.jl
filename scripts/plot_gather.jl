#!/usr/bin/env julia
# ==============================================================================
# scripts/plot_gather.jl
#
# Plot shot gather from binary file.
#
# Usage:
#   julia scripts/plot_gather.jl shot_1.bin 600 100 0.00167
#   julia scripts/plot_gather.jl shot_1.bin 600 100 0.00167 output.png
#
# Arguments:
#   1. filename: Path to binary gather file
#   2. nt: Number of time samples
#   3. n_rec: Number of receivers
#   4. dt: Time step (seconds)
#   5. output (optional): Output image path (default: display)
# ==============================================================================

using CairoMakie
using Printf

function load_gather(filename::String, nt::Int, n_rec::Int)
    data = zeros(Float32, nt, n_rec)
    open(filename, "r") do io
        read!(io, data)
    end
    return data
end

function plot_gather(data::Matrix{Float32}, dt::Real; 
                     title::String="Shot Gather",
                     clip_percentile::Float64=98.0)
    nt, n_rec = size(data)
    t = range(0, step=dt, length=nt)
    
    # Auto clip
    sorted = sort(abs.(vec(data)))
    idx = round(Int, clip_percentile / 100 * length(sorted))
    vmax = sorted[min(idx, length(sorted))]
    vmax = vmax > 0 ? vmax : 1.0f0
    
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], 
              title=title,
              xlabel="Receiver Index",
              ylabel="Time (s)",
              yreversed=true)
    
    hm = heatmap!(ax, 1:n_rec, collect(t), data',
                  colormap=:seismic,
                  colorrange=(-vmax, vmax))
    
    Colorbar(fig[1, 2], hm, label="Amplitude")
    
    return fig
end

function plot_wiggle(data::Matrix{Float32}, dt::Real;
                     title::String="Shot Gather (Wiggle)",
                     trace_skip::Int=1,
                     scale::Float64=1.0)
    nt, n_rec = size(data)
    t = range(0, step=dt, length=nt)
    
    # Normalize traces
    max_amp = maximum(abs.(data))
    data_norm = data ./ (max_amp + 1e-10)
    
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1],
              title=title,
              xlabel="Receiver Index", 
              ylabel="Time (s)",
              yreversed=true)
    
    for r in 1:trace_skip:n_rec
        trace = data_norm[:, r] .* scale .+ r
        lines!(ax, trace, collect(t), color=:black, linewidth=0.5)
    end
    
    return fig
end

# Main
function main()
    if length(ARGS) < 4
        println("Usage: julia plot_gather.jl <filename> <nt> <n_rec> <dt> [output.png]")
        println("Example: julia plot_gather.jl shot_1.bin 600 100 0.00167")
        return
    end
    
    filename = ARGS[1]
    nt = parse(Int, ARGS[2])
    n_rec = parse(Int, ARGS[3])
    dt = parse(Float64, ARGS[4])
    output = length(ARGS) >= 5 ? ARGS[5] : nothing
    
    @info "Loading gather" filename=filename nt=nt n_rec=n_rec dt=dt
    
    data = load_gather(filename, nt, n_rec)
    
    @info "Data statistics" min=minimum(data) max=maximum(data) rms=sqrt(sum(data.^2)/length(data))
    
    fig = plot_gather(data, dt; title="Shot Gather: $(basename(filename))")
    
    if output !== nothing
        save(output, fig)
        @info "Saved to $output"
    else
        display(fig)
        println("Press Enter to exit...")
        readline()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
