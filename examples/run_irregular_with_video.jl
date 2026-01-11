# ==============================================================================
# run_irregular_with_video.jl
#
# Irregular free surface simulation - using high-level API
#
# Key point: YOU define the surface shape (z_surface)!
# ==============================================================================

using Fomo

# ==============================================================================
# Example 1: Using surface helper functions
# ==============================================================================

println("="^60)
println("  Example 1: Surface helper functions")
println("="^60)

nx, nz = 400, 200
dx = 10.0f0

vp = fill(3000.0f0, nz, nx)
vs = fill(1800.0f0, nz, nx)
rho = fill(2200.0f0, nz, nx)
vp[100:end, :] .= 4000.0f0
vs[100:end, :] .= 2400.0f0
model = VelocityModel(vp, vs, rho, dx, dx; name="Two-layer model")

# Choose a surface shape:
z_surface = sinusoidal_surface(nx, dx; base_depth=50, amplitude=25, wavelength=1500)

# Other options:
# z_surface = gaussian_valley(nx, dx; base_depth=40, valley_depth=35, width=300)
# z_surface = gaussian_hill(nx, dx; base_depth=80, hill_height=40, width=400)
# z_surface = tilted_surface(nx, dx; depth_left=30, depth_right=70)
# z_surface = random_surface(nx, dx; base_depth=50, amplitude=15, smoothness=8)
# z_surface = combine_surfaces(
#     sinusoidal_surface(nx, dx; base_depth=50, amplitude=20),
#     gaussian_valley(nx, dx; base_depth=0, valley_depth=25, width=250)
# )

println("Surface: $(minimum(z_surface)) - $(maximum(z_surface)) m")

# Video config - separate parameter
video_cfg = VideoConfig(fields=[:vz], skip=10, fps=30)

result = simulate_irregular!(
    model,
    z_surface,
    2000.0f0,
    Float32.(100:20:3900);
    config=IrregularSurfaceConfig(
        nt=3000,
        f0=15.0f0,
        ibm_method=:mirror,
        src_depth=30.0f0,
        output_dir="outputs/example1_sinusoidal"
    ),
    video_config=video_cfg
)


# ==============================================================================
# Example 2: Fully custom surface shape
# ==============================================================================

println("\n" * "="^60)
println("  Example 2: Custom surface shape")
println("="^60)

# Define ANY shape using plain Julia code
x = Float32.((0:nx-1) .* dx)
z_custom = Float32.(
    60.0 .+
    20.0 .* sin.(2π .* x ./ 1500.0) .+
    10.0 .* sin.(2π .* x ./ 300.0) .+
    5.0 .* sin.(2π .* x ./ 80.0)
)

# Add a sharp canyon
for i in 1:nx
    xi = (i - 1) * dx
    if 1300 < xi < 1700
        z_custom[i] += 30.0f0 * (1 - abs(xi - 1500) / 200)
    end
end

println("Custom surface: $(minimum(z_custom)) - $(maximum(z_custom)) m")

result2 = simulate_irregular!(
    model, z_custom, 2500.0f0, Float32.(100:20:3900);
    config=IrregularSurfaceConfig(
        nt=3000,
        ibm_method=:mirror,
        src_depth=40.0f0,
        output_dir="outputs/example2_custom"
    ),
    video_config=VideoConfig(fields=[:vz], skip=10)
)


println("\nAll examples complete!")