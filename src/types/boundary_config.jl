# ==============================================================================
# types/boundary_config.jl
#
# Boundary configuration for ElasticWave2D.jl
# ==============================================================================

"""
    BoundaryConfig(; kwargs...)

Configuration for boundary conditions on all sides of the simulation domain.

# Keyword Arguments
- `top_boundary::Symbol = :image`: Top boundary condition.
  - `:image`: Image method (stress-free surface).
  - `:absorbing`: Absorbing boundary (no surface waves).
  - `:vacuum`: Vacuum formulation (for topography).
- `bottom_boundary::Symbol = :absorbing`: Bottom boundary condition (`:absorbing`, `:image`).
- `left_boundary::Symbol = :absorbing`: Left boundary condition (`:absorbing`, `:periodic`).
- `right_boundary::Symbol = :absorbing`: Right boundary condition (`:absorbing`, `:periodic`).
- `nbc::Int = 50`: Number of absorbing boundary layers (PML/HABC thickness).
- `top_padding::Int = 0`: Additional padding layers at the top (useful if source is very close to boundary).

# Returns
- `BoundaryConfig`: A configuration object.

# Example
```julia
# Standard survey with free surface
bc = BoundaryConfig(top_boundary=:image, nbc=50)

# Full absorbing boundaries (infinite medium)
bc_abs = BoundaryConfig(top_boundary=:absorbing)
```
"""
struct BoundaryConfig
    top_boundary::Symbol    # :absorbing, :image, :vacuum
    bottom_boundary::Symbol # :absorbing, :image
    left_boundary::Symbol   # :absorbing, :periodic
    right_boundary::Symbol  # :absorbing, :periodic
    nbc::Int               # 吸收边界层数
    top_padding::Int       # 顶部额外填充层（用于震源靠近边界的情况）
end

# Default constructor
function BoundaryConfig(;
    top_boundary::Symbol=:image,
    bottom_boundary::Symbol=:absorbing,
    left_boundary::Symbol=:absorbing,
    right_boundary::Symbol=:absorbing,
    nbc::Int=50,
    top_padding::Int=0
)
    BoundaryConfig(top_boundary, bottom_boundary, left_boundary, right_boundary, nbc, top_padding)
end