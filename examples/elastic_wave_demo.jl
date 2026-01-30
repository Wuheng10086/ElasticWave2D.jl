using ElasticWave2D

function elastic_wave_demo()
    println("ElasticWave2D.jl - Elastic Wave Phenomena Demo (Two-Layer Model)")
    println("======================================================")

    # 1. Define Model (Balanced Resolution)
    # 0.5m spacing => 3200 x 2000 grid points
    # Keep physical size same: 1600m x 1000m
    dx, dz = 0.5f0, 0.5f0
    nx = round(Int, 1600.0 / dx)
    nz = round(Int, 1000.0 / dz)

    println("1. Creating two-layer velocity model...")
    println("   Grid: $nx x $nz (dx=dz=$dx m)")

    # Initialize with top layer properties
    vp = fill(3000.0f0, nz, nx)
    vs = fill(1800.0f0, nz, nx)
    rho = fill(2200.0f0, nz, nx)

    # Define interface depth at 500m
    interface_z_idx = round(Int, 500.0 / dz)

    # Set bottom layer properties (High contrast)
    # Z range: interface to bottom
    vp[interface_z_idx:end, :] .= 4000.0f0
    vs[interface_z_idx:end, :] .= 2400.0f0
    rho[interface_z_idx:end, :] .= 2600.0f0

    model = VelocityModel(vp, vs, rho, dx, dz; name="two_layer_model")

    # 2. Define Source and Receivers
    println("2. Setting up source and receivers...")
    # Source near the top center
    src_x = 800.0f0
    src_z = 20.0f0   # Shallow source

    # Receivers on the surface
    rec_x = Float32.(collect(50:10:1550))
    rec_z = fill(2.0f0, length(rec_x)) # Near surface

    # Ensure output directory exists
    mkpath("outputs/elastic_wave_demo")

    # Plot Setup
    plot_setup(model, [src_x], [src_z], rec_x, rec_z;
        title="Two-Layer Model Setup",
        output="outputs/elastic_wave_demo/setup.png"
    )

    # 3. Run Simulation
    println("3. Running simulation...")
    # Higher frequency source (40Hz)

    config = SimulationConfig(
        nt=12000,           # Sufficient steps for 0.5m grid
        f0=40.0f0,          # Higher frequency
        cfl=0.3f0,          # Safe CFL
        free_surface=true,  # Explicit free surface
        output_dir="outputs/elastic_wave_demo",
        save_gather=true,
        plot_gather=true,
        show_progress=true
    )

    # Adjust video skip to avoid generating too many frames
    # 12000 steps / 20fps / 20 sec video = ~30 steps per frame
    video_config = VideoConfig(
        fields=[:vz],  # Only record Vz
        skip=50,       # Save every 50th step (240 frames total)
        fps=20,        # Slower playback
        colormap=:seismic # Default seismic (Red-White-Blue)
    )

    result = simulate!(model, src_x, src_z, rec_x, rec_z;
        config=config,
        video_config=video_config
    )

    println("\nSimulation complete!")
    println("Results saved to: $(config.output_dir)")
    println(" - Check 'wavefield_vz.mp4' to see reflection/refraction")
end

# Run the demo
elastic_wave_demo()
