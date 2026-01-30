using ElasticWave2D

"""
    exploration_seismic_demo()

Exploration Seismology: Imaging an anticlinal structure.

This demo simulates a classic petroleum exploration scenario: imaging a 
subsurface anticline (fold structure) that could serve as a hydrocarbon trap.

Uses the `seismic_survey` API with vacuum formulation for the free surface.

Physical setup:
- Deep geological model (3km depth)
- Anticline structure at ~1.5km depth
- Multi-layer sedimentary sequence
- Surface seismic reflection survey
"""
function exploration_seismic_demo()
    println("ElasticWave2D.jl - Exploration Seismic Demo")
    println("===========================================")
    println("Scenario: Imaging subsurface anticline structure")
    println()

    # ===========================
    # 1. Model Parameters
    # ===========================
    nx, nz = 600, 300          # Grid size
    dx, dz = 10.0f0, 10.0f0    # 10m grid spacing

    println("1. Creating geological model...")
    println("   Grid: $(nx) × $(nz) points")
    println("   Physical size: $(nx*dx/1000)km × $(nz*dz/1000)km")
    println("   Grid spacing: $(dx)m")

    # ===========================
    # 2. Build Sedimentary Basin Model
    # ===========================
    # Layer structure (from top to bottom):
    # - Layer 1: Weathering layer (0-200m)
    # - Layer 2: Shale overburden (200m - interface1)
    # - Layer 3: Sandstone reservoir (interface1 - interface2) ← ANTICLINE
    # - Layer 4: Shale seal (interface2 - interface3)
    # - Layer 5: Basement (interface3+)

    vp = zeros(Float32, nz, nx)
    vs = zeros(Float32, nz, nx)
    rho = zeros(Float32, nz, nx)

    # Define anticline geometry
    function anticline_depth(x, base_depth, amplitude, wavelength)
        center_x = nx * dx / 2
        return base_depth - amplitude * exp(-((x - center_x)^2) / (2 * (wavelength / 2)^2))
    end

    # Interface depths (at model center)
    base_interface1 = 800.0f0    # Top of reservoir
    base_interface2 = 1200.0f0   # Bottom of reservoir
    base_interface3 = 2000.0f0   # Top of basement

    anticline_amplitude = 300.0f0   # 300m structural relief
    anticline_wavelength = 2000.0f0 # 2km wide structure

    for ix in 1:nx
        x = (ix - 1) * dx

        interface1 = anticline_depth(x, base_interface1, anticline_amplitude, anticline_wavelength)
        interface2 = anticline_depth(x, base_interface2, anticline_amplitude, anticline_wavelength)
        interface3 = anticline_depth(x, base_interface3, anticline_amplitude * 0.5f0, anticline_wavelength * 1.5f0)

        for iz in 1:nz
            depth = (iz - 1) * dz

            if depth < 200
                # Layer 1: Weathering layer
                vp[iz, ix] = 2000.0f0
                vs[iz, ix] = 800.0f0
                rho[iz, ix] = 1900.0f0
            elseif depth < interface1
                # Layer 2: Shale overburden
                vp[iz, ix] = 3000.0f0
                vs[iz, ix] = 1500.0f0
                rho[iz, ix] = 2300.0f0
            elseif depth < interface2
                # Layer 3: Sandstone reservoir (anticline)
                vp[iz, ix] = 3800.0f0
                vs[iz, ix] = 2200.0f0
                rho[iz, ix] = 2400.0f0
            elseif depth < interface3
                # Layer 4: Shale seal
                vp[iz, ix] = 3200.0f0
                vs[iz, ix] = 1700.0f0
                rho[iz, ix] = 2350.0f0
            else
                # Layer 5: Basement
                vp[iz, ix] = 5000.0f0
                vs[iz, ix] = 2900.0f0
                rho[iz, ix] = 2700.0f0
            end
        end
    end

    println("2. Geological structure:")
    println("   Layer 1: Weathering layer (Vp=2000 m/s)")
    println("   Layer 2: Shale overburden (Vp=3000 m/s)")
    println("   Layer 3: Sandstone reservoir - ANTICLINE (Vp=3800 m/s)")
    println("   Layer 4: Shale seal (Vp=3200 m/s)")
    println("   Layer 5: Basement (Vp=5000 m/s)")
    println("   Anticline amplitude: $(anticline_amplitude)m")

    model = VelocityModel(vp, vs, rho, dx, dz; name="anticline_model")

    # ===========================
    # 3. Survey Geometry
    # ===========================
    # Source and receivers at surface (z = 0, relative to model top)
    # seismic_survey will add vacuum layers and adjust coordinates automatically
    src_x = 500.0f0
    src_z = 20.0f0     # 20m below surface (will be adjusted by seismic_survey)

    rec_x = Float32.(collect(100:25:5900))
    rec_z = fill(20.0f0, length(rec_x))

    println("3. Survey geometry:")
    println("   Source: x = $(src_x)m")
    println("   Receivers: $(length(rec_x)) geophones, spacing = 25m")

    # ===========================
    # 4. Output Directory
    # ===========================
    output_dir = "outputs/exploration_seismic"
    mkpath(output_dir)

    # Plot setup (original model without vacuum)
    plot_setup(model, [src_x], [src_z], rec_x, rec_z;
        title="Exploration Seismic: Anticline Structure",
        output=joinpath(output_dir, "setup.png")
    )
    println("   Setup plot saved.")

    # ===========================
    # 5. Run Simulation
    # ===========================
    println("4. Running simulation...")
    println("   Using seismic_survey API with vacuum formulation")

    config = SimulationConfig(
        nt=8000,              # Long record for deep reflections
        f0=25.0f0,            # 25 Hz (typical exploration source)
        nbc=60,               # Thick absorbing boundary
        fd_order=8,
        output_dir=output_dir,
        save_gather=true,
        plot_gather=true,
        show_progress=true
    )

    # Use seismic_survey with vacuum formulation
    result = seismic_survey(
        model,
        (src_x, src_z),         # Source position
        (rec_x, rec_z);         # Receiver positions
        surface_method=:vacuum,
        vacuum_layers=10,     # 100m vacuum (10 layers × 10m)
        config=config
    )

    # ===========================
    # 6. Summary
    # ===========================
    println()
    println("="^55)
    println("Simulation Complete!")
    println("="^55)
    println()
    println("Results saved to: $(output_dir)/")
    println()
    println("Key outputs:")
    println("  • setup.png        - Velocity model with anticline")
    println("  • gather.png       - Shot gather with reflections")
    println()
    println("What to look for in the shot gather:")
    println("  1. Direct wave (first arrival)")
    println("  2. Reflection from top of reservoir (anticline shape)")
    println("  3. Reflection from bottom of reservoir")
    println("  4. Reflection from basement")
    println("  5. \"Pull-up\" in reflection times at anticline crest")
    println()
    println("This anticline structure is a classic petroleum trap where")
    println("hydrocarbons migrate upward and accumulate at the crest.")
end

# Run the demo
exploration_seismic_demo()