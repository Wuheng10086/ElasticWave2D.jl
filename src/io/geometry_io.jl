# ==============================================================================
# io/geometry_io.jl
#
# Geometry I/O - Save/Load survey geometry for migration
# 
# IMPORTANT: All coordinates are the ACTUAL discretized positions on the grid,
#            not the originally requested positions!
# ==============================================================================

# JLD2 and Printf are already imported at module level

# ==============================================================================
# Geometry Structure
# ==============================================================================

"""
    SurveyGeometry

Complete survey geometry for a shot gather.
Contains all information needed for migration.

**Note**: All coordinates (`src_x`, `src_z`, `rec_x`, `rec_z`) are the 
ACTUAL positions after grid discretization, not the originally requested values.
The grid indices (`src_i`, `src_j`, `rec_i_idx`, `rec_j_idx`) are also provided.
"""
struct SurveyGeometry
    # Shot info
    shot_id::Int

    # Source position - ACTUAL discretized position
    src_x::Float32          # Source X position (m) - discretized
    src_z::Float32          # Source Z position (m) - discretized
    src_i::Int              # Source grid index in X
    src_j::Int              # Source grid index in Z

    # Receiver info - ACTUAL discretized positions
    n_rec::Int              # Number of receivers
    rec_x::Vector{Float32}  # Receiver X positions (m) - discretized
    rec_z::Vector{Float32}  # Receiver Z positions (m) - discretized
    rec_i_idx::Vector{Int}  # Receiver grid indices in X
    rec_j_idx::Vector{Int}  # Receiver grid indices in Z

    # Time info
    dt::Float32             # Time sampling (s)
    nt::Int                 # Number of time samples
    t_max::Float32          # Total recording time (s)

    # Grid info (physical domain, excluding padding)
    dx::Float32
    dz::Float32
    nx::Int
    nz::Int
end

"""
    MultiShotGeometry

Geometry for multiple shots (entire survey).
"""
struct MultiShotGeometry
    n_shots::Int
    shots::Vector{SurveyGeometry}

    # Common parameters
    dt::Float32
    nt::Int
    dx::Float32
    dz::Float32
    nx::Int
    nz::Int
end

# ==============================================================================
# Create Geometry from Simulation Results
# ==============================================================================

"""
    create_geometry(result::ShotResult, medium::Medium, params::SimParams) -> SurveyGeometry

Create geometry from a single shot result.

The output coordinates are the ACTUAL discretized grid positions,
calculated from the grid indices used during simulation.
"""
function create_geometry(result::ShotResult, medium::Medium, params::SimParams)
    pad = medium.pad

    # Convert grid indices to ACTUAL physical coordinates
    # These are the discretized positions, not the original input!
    src_x = Float32((result.src_i - pad - 1) * medium.dx)
    src_z = Float32((result.src_j - pad - 1) * medium.dz)

    n_rec = length(result.rec_i)
    rec_x = Float32[(result.rec_i[r] - pad - 1) * medium.dx for r in 1:n_rec]
    rec_z = Float32[(result.rec_j[r] - pad - 1) * medium.dz for r in 1:n_rec]

    # Grid indices relative to physical domain (0-based, for external use)
    src_i_phys = result.src_i - pad
    src_j_phys = result.src_j - pad
    rec_i_phys = [result.rec_i[r] - pad for r in 1:n_rec]
    rec_j_phys = [result.rec_j[r] - pad for r in 1:n_rec]

    return SurveyGeometry(
        result.shot_id,
        src_x, src_z, src_i_phys, src_j_phys,
        n_rec,
        rec_x, rec_z, rec_i_phys, rec_j_phys,
        params.dt, params.nt, params.dt * params.nt,
        medium.dx, medium.dz,
        medium.nx - 2 * pad, medium.nz - 2 * pad  # Physical grid size
    )
end

"""
    create_geometry(results::Vector{ShotResult}, medium::Medium, params::SimParams) -> MultiShotGeometry

Create geometry from multiple shot results.
"""
function create_geometry(results::Vector{ShotResult}, medium::Medium, params::SimParams)
    shots = [create_geometry(r, medium, params) for r in results]
    pad = medium.pad

    return MultiShotGeometry(
        length(shots),
        shots,
        params.dt, params.nt,
        medium.dx, medium.dz,
        medium.nx - 2 * pad, medium.nz - 2 * pad
    )
end

# ==============================================================================
# Save Geometry
# ==============================================================================

"""
    save_geometry(path, geom::SurveyGeometry)
    save_geometry(path, geom::MultiShotGeometry)

Save geometry to file. Supports .jld2, .json, .txt formats.

**Note**: All coordinates in the output are ACTUAL discretized grid positions.

# Example
```julia
geom = create_geometry(result, medium, params)
save_geometry("shot_001_geom.jld2", geom)
save_geometry("shot_001_geom.json", geom)
save_geometry("shot_001_geom.txt", geom)
```
"""
function save_geometry(path::String, geom::SurveyGeometry)
    ext = lowercase(splitext(path)[2])

    if ext == ".jld2"
        _save_geometry_jld2(path, geom)
    elseif ext == ".json"
        _save_geometry_json(path, geom)
    elseif ext in [".txt", ".dat", ".geo"]
        _save_geometry_txt(path, geom)
    else
        error("Unsupported format: $ext. Use .jld2, .json, or .txt")
    end
end

function save_geometry(path::String, geom::MultiShotGeometry)
    ext = lowercase(splitext(path)[2])

    if ext == ".jld2"
        _save_multigeom_jld2(path, geom)
    elseif ext == ".json"
        _save_multigeom_json(path, geom)
    elseif ext in [".txt", ".dat", ".geo"]
        _save_multigeom_txt(path, geom)
    else
        error("Unsupported format: $ext. Use .jld2, .json, or .txt")
    end
end

# ==============================================================================
# JLD2 Format (Recommended)
# ==============================================================================

function _save_geometry_jld2(path::String, g::SurveyGeometry)
    jldsave(path;
        # Meta
        shot_id=g.shot_id,
        # Source - actual discretized position
        src_x=g.src_x, src_z=g.src_z,
        src_i=g.src_i, src_j=g.src_j,
        # Receivers - actual discretized positions
        n_rec=g.n_rec,
        rec_x=g.rec_x, rec_z=g.rec_z,
        rec_i=g.rec_i_idx, rec_j=g.rec_j_idx,
        # Time
        dt=g.dt, nt=g.nt, t_max=g.t_max,
        # Grid
        dx=g.dx, dz=g.dz, nx=g.nx, nz=g.nz,
        # Note
        _note="All coordinates are ACTUAL discretized grid positions"
    )
    @info "Geometry saved (actual discretized positions)" path = path
end

function _save_multigeom_jld2(path::String, mg::MultiShotGeometry)
    n = mg.n_shots
    shot_ids = [s.shot_id for s in mg.shots]
    src_x = [s.src_x for s in mg.shots]
    src_z = [s.src_z for s in mg.shots]
    src_i = [s.src_i for s in mg.shots]
    src_j = [s.src_j for s in mg.shots]
    n_rec = [s.n_rec for s in mg.shots]
    rec_x = [s.rec_x for s in mg.shots]
    rec_z = [s.rec_z for s in mg.shots]
    rec_i = [s.rec_i_idx for s in mg.shots]
    rec_j = [s.rec_j_idx for s in mg.shots]

    jldsave(path;
        n_shots=n,
        shot_ids=shot_ids,
        src_x=src_x, src_z=src_z,
        src_i=src_i, src_j=src_j,
        n_rec=n_rec,
        rec_x=rec_x, rec_z=rec_z,
        rec_i=rec_i, rec_j=rec_j,
        dt=mg.dt, nt=mg.nt,
        dx=mg.dx, dz=mg.dz, nx=mg.nx, nz=mg.nz,
        _note="All coordinates are ACTUAL discretized grid positions"
    )
    @info "Multi-shot geometry saved (actual discretized positions)" path = path n_shots = n
end

# ==============================================================================
# JSON Format (Cross-language compatible)
# ==============================================================================

function _save_geometry_json(path::String, g::SurveyGeometry)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"_note\": \"All coordinates are ACTUAL discretized grid positions\",")
        println(io, "  \"shot_id\": $(g.shot_id),")
        println(io, "  \"source\": {")
        println(io, "    \"x\": $(g.src_x), \"z\": $(g.src_z),")
        println(io, "    \"i\": $(g.src_i), \"j\": $(g.src_j)")
        println(io, "  },")
        println(io, "  \"receivers\": {")
        println(io, "    \"n\": $(g.n_rec),")
        println(io, "    \"x\": [$(join(g.rec_x, ", "))],")
        println(io, "    \"z\": [$(join(g.rec_z, ", "))],")
        println(io, "    \"i\": [$(join(g.rec_i_idx, ", "))],")
        println(io, "    \"j\": [$(join(g.rec_j_idx, ", "))]")
        println(io, "  },")
        println(io, "  \"time\": {\"dt\": $(g.dt), \"nt\": $(g.nt), \"t_max\": $(g.t_max)},")
        println(io, "  \"grid\": {\"dx\": $(g.dx), \"dz\": $(g.dz), \"nx\": $(g.nx), \"nz\": $(g.nz)}")
        println(io, "}")
    end
    @info "Geometry saved (actual discretized positions)" path = path format = "JSON"
end

function _save_multigeom_json(path::String, mg::MultiShotGeometry)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"_note\": \"All coordinates are ACTUAL discretized grid positions\",")
        println(io, "  \"n_shots\": $(mg.n_shots),")
        println(io, "  \"common\": {")
        println(io, "    \"dt\": $(mg.dt), \"nt\": $(mg.nt),")
        println(io, "    \"dx\": $(mg.dx), \"dz\": $(mg.dz),")
        println(io, "    \"nx\": $(mg.nx), \"nz\": $(mg.nz)")
        println(io, "  },")
        println(io, "  \"shots\": [")
        for (idx, s) in enumerate(mg.shots)
            comma = idx < mg.n_shots ? "," : ""
            println(io, "    {")
            println(io, "      \"id\": $(s.shot_id),")
            println(io, "      \"src_x\": $(s.src_x), \"src_z\": $(s.src_z),")
            println(io, "      \"src_i\": $(s.src_i), \"src_j\": $(s.src_j),")
            println(io, "      \"n_rec\": $(s.n_rec),")
            println(io, "      \"rec_x\": [$(join(s.rec_x, ", "))],")
            println(io, "      \"rec_z\": [$(join(s.rec_z, ", "))],")
            println(io, "      \"rec_i\": [$(join(s.rec_i_idx, ", "))],")
            println(io, "      \"rec_j\": [$(join(s.rec_j_idx, ", "))]")
            println(io, "    }$comma")
        end
        println(io, "  ]")
        println(io, "}")
    end
    @info "Multi-shot geometry saved (actual discretized positions)" path = path n_shots = mg.n_shots format = "JSON"
end

# ==============================================================================
# Text Format (Human readable, easy to parse)
# ==============================================================================

function _save_geometry_txt(path::String, g::SurveyGeometry)
    open(path, "w") do io
        println(io, "# Survey Geometry - Single Shot")
        println(io, "# Generated by ElasticWave2D")
        println(io, "# NOTE: All coordinates are ACTUAL DISCRETIZED grid positions!")
        println(io, "#" * "="^60)
        println(io, "")
        println(io, "# Shot Info")
        println(io, "shot_id      $(g.shot_id)")
        println(io, "")
        println(io, "# Source (actual discretized position)")
        @printf(io, "src_x        %.4f    # meters\n", g.src_x)
        @printf(io, "src_z        %.4f    # meters\n", g.src_z)
        println(io, "src_i        $(g.src_i)        # grid index (0-based)")
        println(io, "src_j        $(g.src_j)        # grid index (0-based)")
        println(io, "")
        println(io, "# Time Info")
        @printf(io, "dt           %.6f  # seconds\n", g.dt)
        println(io, "nt           $(g.nt)")
        @printf(io, "t_max        %.4f    # seconds\n", g.t_max)
        println(io, "")
        println(io, "# Grid Info (physical domain)")
        @printf(io, "dx           %.4f    # meters\n", g.dx)
        @printf(io, "dz           %.4f    # meters\n", g.dz)
        println(io, "nx           $(g.nx)")
        println(io, "nz           $(g.nz)")
        println(io, "")
        println(io, "# Receivers: $(g.n_rec) total (actual discretized positions)")
        println(io, "# Format: rec_id  x(m)  z(m)  i_idx  j_idx")
        println(io, "n_rec        $(g.n_rec)")
        println(io, "receivers")
        for r in 1:g.n_rec
            @printf(io, "%6d  %12.4f  %12.4f  %6d  %6d\n",
                r, g.rec_x[r], g.rec_z[r], g.rec_i_idx[r], g.rec_j_idx[r])
        end
    end
    @info "Geometry saved (actual discretized positions)" path = path format = "TXT"
end

function _save_multigeom_txt(path::String, mg::MultiShotGeometry)
    open(path, "w") do io
        println(io, "# Survey Geometry - Multi Shot")
        println(io, "# Generated by ElasticWave2D")
        println(io, "# NOTE: All coordinates are ACTUAL DISCRETIZED grid positions!")
        println(io, "#" * "="^60)
        println(io, "")
        println(io, "# Common Parameters")
        println(io, "n_shots      $(mg.n_shots)")
        @printf(io, "dt           %.6f  # seconds\n", mg.dt)
        println(io, "nt           $(mg.nt)")
        @printf(io, "dx           %.4f    # meters\n", mg.dx)
        @printf(io, "dz           %.4f    # meters\n", mg.dz)
        println(io, "nx           $(mg.nx)")
        println(io, "nz           $(mg.nz)")
        println(io, "")
        println(io, "# Shot List (actual discretized positions)")
        println(io, "# Format: shot_id  src_x(m)  src_z(m)  src_i  src_j  n_rec")
        println(io, "shots")
        for s in mg.shots
            @printf(io, "%6d  %12.4f  %12.4f  %6d  %6d  %6d\n",
                s.shot_id, s.src_x, s.src_z, s.src_i, s.src_j, s.n_rec)
        end
        println(io, "")
        println(io, "# Receiver Coordinates per Shot (actual discretized positions)")
        println(io, "# Format: rec_id  x(m)  z(m)  i_idx  j_idx")
        for s in mg.shots
            println(io, "")
            println(io, "# Shot $(s.shot_id) receivers")
            println(io, "shot_$(s.shot_id)_receivers  $(s.n_rec)")
            for r in 1:s.n_rec
                @printf(io, "%6d  %12.4f  %12.4f  %6d  %6d\n",
                    r, s.rec_x[r], s.rec_z[r], s.rec_i_idx[r], s.rec_j_idx[r])
            end
        end
    end
    @info "Multi-shot geometry saved (actual discretized positions)" path = path n_shots = mg.n_shots format = "TXT"
end

# ==============================================================================
# Load Geometry
# ==============================================================================

"""
    load_geometry(path) -> SurveyGeometry or MultiShotGeometry

Load geometry from file.
"""
function load_geometry(path::String)
    ext = lowercase(splitext(path)[2])

    if ext == ".jld2"
        return _load_geometry_jld2(path)
    else
        error("Only JLD2 loading is currently supported. Convert with save_geometry first.")
    end
end

function _load_geometry_jld2(path::String)
    data = JLD2.load(path)

    if haskey(data, "n_shots")
        # Multi-shot geometry
        shots = SurveyGeometry[]
        for i in 1:data["n_shots"]
            push!(shots, SurveyGeometry(
                data["shot_ids"][i],
                data["src_x"][i], data["src_z"][i],
                data["src_i"][i], data["src_j"][i],
                data["n_rec"][i],
                data["rec_x"][i], data["rec_z"][i],
                data["rec_i"][i], data["rec_j"][i],
                data["dt"], data["nt"], data["dt"] * data["nt"],
                data["dx"], data["dz"], data["nx"], data["nz"]
            ))
        end
        return MultiShotGeometry(
            data["n_shots"], shots,
            data["dt"], data["nt"],
            data["dx"], data["dz"], data["nx"], data["nz"]
        )
    else
        # Single shot geometry
        return SurveyGeometry(
            data["shot_id"],
            data["src_x"], data["src_z"],
            data["src_i"], data["src_j"],
            data["n_rec"],
            data["rec_x"], data["rec_z"],
            data["rec_i"], data["rec_j"],
            data["dt"], data["nt"], data["t_max"],
            data["dx"], data["dz"], data["nx"], data["nz"]
        )
    end
end

# ==============================================================================
# Print Geometry Info
# ==============================================================================

function Base.show(io::IO, g::SurveyGeometry)
    println(io, "SurveyGeometry (Shot #$(g.shot_id))")
    println(io, "  Source:    ($(g.src_x), $(g.src_z)) m  [grid: $(g.src_i), $(g.src_j)]")
    println(io, "  Receivers: $(g.n_rec)")
    println(io, "    X range: $(minimum(g.rec_x)) - $(maximum(g.rec_x)) m")
    println(io, "    Z range: $(minimum(g.rec_z)) - $(maximum(g.rec_z)) m")
    println(io, "  Time:      $(g.nt) samples, dt=$(g.dt*1000) ms, T=$(g.t_max) s")
    println(io, "  Grid:      $(g.nx) × $(g.nz), dx=$(g.dx) m")
    println(io, "  Note:      Coordinates are ACTUAL discretized positions")
end

function Base.show(io::IO, mg::MultiShotGeometry)
    println(io, "MultiShotGeometry ($(mg.n_shots) shots)")
    src_x = [s.src_x for s in mg.shots]
    println(io, "  Sources X: $(minimum(src_x)) - $(maximum(src_x)) m")
    println(io, "  Receivers: $(mg.shots[1].n_rec) per shot")
    println(io, "  Time:      $(mg.nt) samples, dt=$(mg.dt*1000) ms")
    println(io, "  Grid:      $(mg.nx) × $(mg.nz), dx=$(mg.dx) m")
    println(io, "  Note:      Coordinates are ACTUAL discretized positions")
end

# ==============================================================================
# Load Gather with Geometry
# ==============================================================================

"""
    load_gather(path::String, geom::SurveyGeometry) -> Matrix{Float32}

Load gather from binary file using geometry for dimensions.

# Example
```julia
geom = load_geometry("survey.jld2")
gather = load_gather("shot_1.bin", geom)  # or geom.shots[1] for multi-shot
```
"""
function load_gather(path::String, geom::SurveyGeometry)
    data = zeros(Float32, geom.nt, geom.n_rec)
    open(path, "r") do io
        read!(io, data)
    end
    return data
end
