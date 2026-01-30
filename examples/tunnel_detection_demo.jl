using ElasticWave2D

"""
    tunnel_detection_demo()

Engineering Application: Detecting underground tunnels using seismic waves.

This demo simulates a common geotechnical scenario: detecting a shallow tunnel
(e.g., subway, utility corridor, or abandoned mine) using surface seismic methods.

Uses the `seismic_survey` API with vacuum formulation for the free surface.
The tunnel cavity is also modeled using vacuum (ρ=0).

Physical setup:
- Surface seismic survey (source and receivers on ground)
- Shallow tunnel at ~30m depth
- Three-layer geological model (soil → weathered rock → bedrock)
"""
function tunnel_detection_demo()
    println("ElasticWave2D.jl - Tunnel Detection Demo")
    println("========================================")
    println("Engineering scenario: Detecting shallow underground tunnel")
    println()

    # ===========================
    # 1. Model Parameters
    # ===========================
    nx, nz = 400, 200          # Grid size
    dx, dz = 1.0f0, 1.0f0      # 1m grid spacing

    println("1. Creating geological model...")
    println("   Grid: $(nx) × $(nz) points")
    println("   Physical size: $(nx*dx)m × $(nz*dz)m")
    println("   Grid spacing: $(dx)m")

    # ===========================
    # 2. Build Layered Geological Model
    # ===========================
    # Layer 1: Soil/Fill (0-15m)
    # Layer 2: Weathered rock (15-50m)
    # Layer 3: Bedrock (50m+)

    vp = zeros(Float32, nz, nx)
    vs = zeros(Float32, nz, nx)
    rho = zeros(Float32, nz, nx)

    for ix in 1:nx
        for iz in 1:nz
            depth = (iz - 1) * dz

            if depth < 15
                # Layer 1: Soil/Fill
                vp[iz, ix] = 800.0f0
                vs[iz, ix] = 400.0f0
                rho[iz, ix] = 1800.0f0
            elseif depth < 50
                # Layer 2: Weathered rock
                vp[iz, ix] = 2000.0f0
                vs[iz, ix] = 1000.0f0
                rho[iz, ix] = 2200.0f0
            else
                # Layer 3: Bedrock
                vp[iz, ix] = 3500.0f0
                vs[iz, ix] = 2000.0f0
                rho[iz, ix] = 2500.0f0
            end
        end
    end

    # ===========================
    # 3. Add Tunnel (Vacuum Region)
    # ===========================
    tunnel_center_x = 200.0f0   # Center at x = 200m
    tunnel_center_z = 30.0f0    # Depth = 30m
    tunnel_width = 10.0f0       # 10m wide
    tunnel_height = 8.0f0       # 8m tall

    println("2. Adding tunnel (vacuum region)...")
    println("   Position: x = $(tunnel_center_x)m, z = $(tunnel_center_z)m")
    println("   Size: $(tunnel_width)m × $(tunnel_height)m")

    # Convert to grid indices
    ix_start = round(Int, (tunnel_center_x - tunnel_width / 2) / dx) + 1
    ix_end = round(Int, (tunnel_center_x + tunnel_width / 2) / dx) + 1
    iz_start = round(Int, (tunnel_center_z - tunnel_height / 2) / dz) + 1
    iz_end = round(Int, (tunnel_center_z + tunnel_height / 2) / dz) + 1

    ix_start = clamp(ix_start, 1, nx)
    ix_end = clamp(ix_end, 1, nx)
    iz_start = clamp(iz_start, 1, nz)
    iz_end = clamp(iz_end, 1, nz)

    # Set tunnel to vacuum (ρ = 0)
    vp[iz_start:iz_end, ix_start:ix_end] .= 0.0f0
    vs[iz_start:iz_end, ix_start:ix_end] .= 0.0f0
    rho[iz_start:iz_end, ix_start:ix_end] .= 0.0f0

    model = VelocityModel(vp, vs, rho, dx, dz; name="tunnel_model")

    # ===========================
    # 4. Survey Geometry
    # ===========================
    # Source and receivers near surface
    # seismic_survey will add vacuum layers and adjust coordinates
    src_x = 100.0f0
    src_z = 2.0f0     # 2m below surface

    rec_x = Float32.(collect(5:2:395))
    rec_z = fill(2.0f0, length(rec_x))

    println("3. Survey geometry...")
    println("   Source: x = $(src_x)m")
    println("   Receivers: $(length(rec_x)) geophones, spacing = 2m")

    # ===========================
    # 5. Output Directory
    # ===========================
    output_dir = "outputs/tunnel_detection"
    mkpath(output_dir)

    plot_setup(model, [src_x], [src_z], rec_x, rec_z;
        title="Tunnel Detection Survey Setup",
        output=joinpath(output_dir, "setup.png")
    )
    println("   Setup plot saved.")

    # ===========================
    # 6. Run Simulation
    # ===========================
    println("4. Running simulation...")
    println("   Using seismic_survey API with vacuum formulation")

    config = SimulationConfig(
        nt=4000,
        f0=50.0f0,            # 50 Hz (typical for engineering)
        nbc=50,
        fd_order=8,
        output_dir=output_dir,
        save_gather=true,
        plot_gather=true,
        show_progress=true
    )

    # Use seismic_survey with vacuum formulation
    result = seismic_survey(
        model,
        (src_x, src_z),
        (rec_x, rec_z);
        surface_method=:vacuum,
        vacuum_layers=10,     # 10m vacuum at top
        config=config
    )

    # ===========================
    # 7. Summary
    # ===========================
    println()
    println("="^50)
    println("Simulation Complete!")
    println("="^50)
    println()
    println("Results saved to: $(output_dir)/")
    println()
    println("Vacuum formulation used for:")
    println("  • Free surface (via seismic_survey)")
    println("  • Tunnel cavity (ρ=0 in model)")
    println()
    println("Key outputs:")
    println("  • setup.png  - Survey geometry with tunnel")
    println("  • gather.png - Shot gather")
    println()
    println("What to look for:")
    println("  1. Diffracted waves from tunnel edges")
    println("  2. Shadow zone behind tunnel")
    println("  3. Surface wave (Rayleigh) propagation")
end

# Run the demo
tunnel_detection_demo()