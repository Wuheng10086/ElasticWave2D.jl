using ElasticWave2D
using Printf
using Dates

function batch_simulation_demo()
    println("ElasticWave2D.jl - High Performance Batch Simulation Demo")
    println("================================================")

    # 1. Define Model
    # 500 x 200 grid
    nx, nz = 500, 200
    dx, dz = 10.0f0, 10.0f0

    vp = fill(3000.0f0, nz, nx)
    vs = fill(1800.0f0, nz, nx)
    rho = fill(2200.0f0, nz, nx)

    # Add a reflector
    vp[100:end, :] .= 4000.0f0
    vs[100:end, :] .= 2400.0f0
    rho[100:end, :] .= 2500.0f0

    model = VelocityModel(vp, vs, rho, dx, dz; name="batch_model")

    # 2. Define Survey Geometry
    # 20 shots spaced 200m apart
    n_shots = 20
    shot_spacing = 200.0f0
    start_x = 200.0f0
    src_z = 20.0f0

    shots_x = Float32[start_x + (i - 1) * shot_spacing for i in 1:n_shots]
    shots_z = fill(src_z, n_shots)

    # Fixed receiver spread for all shots (full offset)
    # Receivers every 20m across the whole model
    rec_x = Float32.(collect(50:20:4950))
    rec_z = fill(10.0f0, length(rec_x))

    println("Survey Parameters:")
    println("  Model Size: $nx x $nz")
    println("  Shots: $n_shots")
    println("  Receivers per shot: $(length(rec_x))")

    # Ensure output directory exists before saving setup plot
    mkpath("outputs/batch_demo")

    # Plot Setup (Only once)
    plot_setup(model, shots_x, shots_z, rec_x, rec_z;
        title="Batch Simulation Setup",
        output="outputs/batch_demo/setup.png"
    )

    # 3. Initialize Batch Simulator
    println("\nInitializing Simulator...")

    # Simulation config (minimal IO)
    config = SimulationConfig(
        nt=1500,
        f0=25.0f0,
        free_surface=true,
        output_dir="outputs/batch_demo",
        save_gather=false,   # Don't save individual files automatically
        plot_gather=false,   # Don't plot individual files
        show_progress=false  # Disable progress bars for speed
    )

    # Pre-allocate everything on GPU/CPU
    # This is where the heavy lifting happens (once)
    sim = BatchSimulator(model, rec_x, rec_z; config=config)

    println("Simulator ready. Starting batch run...")
    println("----------------------------------------")

    # 4. Run Shots Loop
    t_start = now()

    # Pre-compile / Warm-up
    print("Warm-up shot... ")
    simulate_shot!(sim, shots_x[1], shots_z[1])
    println("Done.")

    total_samples = 0

    # Main Loop
    for i in 1:n_shots
        # Run single shot
        # No info printing inside the loop for max speed
        gather = simulate_shot!(sim, shots_x[i], shots_z[i])

        # Optional: Save gather manually if needed
        # save_gather(gather, joinpath(config.output_dir, "shot_$i.bin"))

        total_samples += length(gather)

        # Minimal progress indicator
        if i % 5 == 0
            print(".")
        end
    end

    t_end = now()
    duration = (t_end - t_start).value / 1000.0 # seconds

    println("\n----------------------------------------")
    println("Batch Run Complete!")
    println(@sprintf("Total Time:     %.2f seconds", duration))
    println(@sprintf("Throughput:     %.2f shots/sec", n_shots / duration))
    println(@sprintf("Avg Time/Shot:  %.2f seconds", duration / n_shots))

    # Save one example gather just to prove it worked
    last_gather = simulate_shot!(sim, shots_x[end], shots_z[end])
    output_dir = config.output_dir
    mkpath(output_dir)

    # Calculate dt manually since config.dt might be nothing (auto-calculated)
    # Or retrieve it from the simulator/model if possible.
    # In this demo, we can just recalculate or assume a safe value for plotting.
    # But better: The simulator object should have the actual dt used.
    # Let's assume standard CFL calculation or just pass a float if we know it.
    # For now, let's recalculate the CFL-based dt for plotting purposes.
    # v_max = 4000.0 (from model definition)
    # dx = 10.0
    # dt <= 0.5 * 10 / 4000 = 0.00125 s = 1.25 ms
    # Let's just use 0.001f0 (1ms) as a reasonable estimate for the plot axis
    # or better, fix the config to have explicit dt.

    dt_plot = 0.001f0 # 1ms

    # Corrected arguments for plot_gather
    plot_gather(last_gather, rec_x, dt_plot;
        output_path=joinpath(output_dir, "example_gather.png"),
        title="Example Shot Gather (Shot $n_shots)"
    )
    println("Saved example gather to: $output_dir/example_gather.png")
end

# Run the demo
batch_simulation_demo()
