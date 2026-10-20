phase_info(::AqueousPhase) = (keyword = "SWAT", mask = 1, slot = 1)
phase_info(::LiquidPhase) = (keyword = "SOIL", mask = 2, slot = 2)
phase_info(::VaporPhase) = (keyword = "SGAS", mask = 4, slot = 3)
phase_info(phase) = throw(ArgumentError("Unsupported phase $(typeof(phase))"))
model_phases(case) = map(phase_info, get_phases(reservoir_model(case).system))
phase_mask(phases) = foldl(|, (p.mask for p in phases); init = 0)

function case_start_date(case::Jutul.JutulCase)
    if !isnothing(case.start_date)
        return DateTime(case.start_date)
    end
    data = case.input_data
    if data isa AbstractDict && haskey(data, "RUNSPEC")
        start = get(data["RUNSPEC"], "START", nothing)
        if !isnothing(start) && !ismissing(start)
            return DateTime(start)
        end
    end
    return DateTime(1970)
end

function case_grid(case::Jutul.JutulCase)
    domain = reservoir_model(case).data_domain
    mesh = Jutul.physical_representation(domain)
    data = case.input_data
    if data isa AbstractDict
        grid = get(data, "GRID", data)
        if haskey(grid, "COORD") || haskey(grid, "DX") || haskey(grid, "DXV")
            grid = Dict{String, Any}(corner_point_data(grid))
            dims = grid_dims(grid)
            grid["ACTNUM"] = get(grid, "ACTNUM", trues(dims))
            if mesh isa Jutul.UnstructuredMesh && mesh.structure isa CartesianIndex
                Tuple(mesh.structure.I) == dims || throw(ArgumentError("Input grid dimensions differ from case mesh"))
                mask = falses(prod(dims))
                mask[mesh.cell_map] .= true
                grid["ACTNUM"] = mask
            end
            count(!iszero, grid["ACTNUM"]) == Jutul.number_of_cells(mesh) ||
                throw(ArgumentError("Input grid active-cell count differs from case mesh"))
            return grid
        end
    end
    return fake_grid(Jutul.number_of_cells(mesh))
end

"""
    write_egrid_jutuldarcy(case, path; kw...)

Write EGRID using `case.input_data` geometry or a synthetic one-dimensional grid.
"""
function write_egrid_jutuldarcy(case::Jutul.JutulCase, path::AbstractString; kw...)
    return GeoEnergyIO.write_egrid(case_grid(case), path; kw...)
end

function reservoir_part(x)
    if haskey(x, :Reservoir)
        return x[:Reservoir]
    end
    return x
end

"""
    write_init_jutuldarcy(case, path; units=:metric, endian=:big)

Write INIT properties from the reservoir domain and parameters.
"""
function write_init_jutuldarcy(case::Jutul.JutulCase, path::AbstractString; units = :metric, endian = :big,
        phases = phase_mask(model_phases(case)))
    units = Symbol(lowercase(string(units)))
    model = reservoir_model(case)
    domain = model.data_domain
    grid = case_grid(case)
    dims = grid_dims(grid)
    active = findall(!iszero, vec(grid["ACTNUM"]))
    nc = length(active)
    kws = standard_headers(dims, nc; units, phases, start_date = case_start_date(case))
    function add(name, vals; full = false)
        if vals isa Number
            vals = fill(vals, nc)
        else
            vals = vec(vals)
        end
        length(vals) == nc || throw(ArgumentError("$name must have $nc cell values"))
        converted = convert_array(name, vals, units; to_si = false)
        if full
            a = zeros(prod(dims))
            a[active] = converted
            converted = a
        end
        push!(kws, BinaryKeyword(name, Float32.(converted)))
    end
    poro = domain[:porosity]
    vols = domain[:volumes]
    if haskey(domain, :net_to_gross)
        ntg = domain[:net_to_gross]
    else
        ntg = ones(nc)
    end
    add("PORO", poro)
    add("NTG", ntg)
    params = reservoir_part(case.parameters)
    pv = get(params, :StaticFluidVolume, get(params, :FluidVolume, vols.*poro.*ntg))
    add("PORV", pv; full = true)
    centroids = domain[:cell_centroids]
    add("DEPTH", centroids[end, :])
    perm = domain[:permeability]
    if perm isa AbstractMatrix
        diagonal = permeability_diagonal_rows(perm, Jutul.dim(Jutul.physical_representation(domain)))
        rows = ntuple(i -> diagonal[min(i, length(diagonal))], 3)
        for (name, row) in zip(("PERMX", "PERMY", "PERMZ"), rows)
            add(name, perm[row, :])
        end
    else
        for name in ("PERMX", "PERMY", "PERMZ")
            add(name, perm)
        end
    end
    # Region arrays are stored as integer cell arrays on the domain.
    for (name, key) in (("SATNUM", :satnum), ("PVTNUM", :pvtnum), ("EQLNUM", :eqlnum), ("FIPNUM", :fipnum))
        haskey(domain, key) && push!(kws, BinaryKeyword(name, Int32.(domain[key])))
    end
    return write_binary(output_path(path, "INIT"), kws; endian)
end

"""
    write_summary_jutuldarcy(case, result, path; units=missing, kw...)

Write SMSPEC and summary data using `summary_result`. Return the written paths
as `(spec, data)`. `units=missing` selects the case's output units.
"""
function write_summary_jutuldarcy(case::Jutul.JutulCase,
        result::ReservoirSimResult, path::AbstractString; units = missing, kw...)
    summary = case_summary(case, result, units)
    return GeoEnergyIO.write_summary(summary, path; kw...)
end

function case_summary(case, result, units)
    summary = summary_result(case, result, units; start_date = case_start_date(case))
    summary["DIMENS"] = grid_dims(case_grid(case))
    return summary
end

"""
    write_restart_jutuldarcy(case, result, path; kw...)
    write_restart_jutuldarcy(case, states, path; times=cumsum(case.dt), kw...)

Write reservoir restart arrays in the case's active-cell and phase ordering.
Times are in seconds. Set `unified=false` for separate Xnnnn files.
"""
function write_restart_jutuldarcy(case::Jutul.JutulCase,
        result::ReservoirSimResult, path::AbstractString; kw...)
    return write_restart_jutuldarcy(case, result.states, path; times = result.time, kw...)
end

function write_restart_jutuldarcy(case::Jutul.JutulCase,
        states::AbstractVector, path::AbstractString; times = cumsum(case.dt), units = :metric,
        unified = true, endian = :big, report_steps = eachindex(states))
    length(states) == length(times) == length(report_steps) || throw(ArgumentError("States, times and report steps must have equal lengths"))
    all(isfinite, times) && all(>=(0), times) && issorted(times) || throw(ArgumentError("Restart times must be finite, nonnegative and sorted"))
    grid = case_grid(case)
    dims = grid_dims(grid)
    nc = count(!iszero, grid["ACTNUM"])
    phases = model_phases(case)
    mask = phase_mask(phases)
    start = case_start_date(case)
    blocks = KeywordFile[]
    for (state, elapsed, report_step) in zip(states, times, report_steps)
        header = standard_headers(dims, nc; units, phases = mask, elapsed, start_date = start, report_step)
        arrays = state_keywords(state, case, nc; units)
        push!(blocks, restart_block(header, arrays; report_step))
    end
    return GeoEnergyIO.write_restart(blocks, path; unified, endian)
end

"""
    write_results_jutuldarcy(case, result, path; units=:metric, unified=true)

Export EGRID, INIT, SMSPEC/UNSMRY and reservoir UNRST solution arrays from a
JutulDarcy case and simulation result. Values are converted from SI using
`summary_result` and deck unit conventions.
"""
function write_results_jutuldarcy(case::Jutul.JutulCase,
        result::ReservoirSimResult, path::AbstractString; units = :metric,
        unified = true, endian = :big)
    grid = write_egrid_jutuldarcy(case, path; units, endian)
    init = write_init_jutuldarcy(case, path; units, endian)
    smry = write_summary_jutuldarcy(case, result, path; units, unified, endian)
    restart = write_restart_jutuldarcy(case, result, path; units, unified, endian)
    return (grid = grid, init = init, summary = smry, restart = restart)
end

"""
    read_restart_jutuldarcy(path, case; steps=missing, extra_out=false)

Read restart arrays as SI state dictionaries using the case's cell, phase and
component ordering. Phase mole fractions are converted to mass fractions.
Missing oil saturation is inferred from water and gas saturations. The
dictionaries can be passed to `setup_reservoir_state`.

`steps` selects report positions. With `extra_out=true`, return `(states,
restart_steps)`, where the second output contains file metadata and native arrays.
"""
function read_restart_jutuldarcy(path::AbstractString, case::Jutul.JutulCase;
        extra_out = false, kw...)
    grid = case_grid(case)
    steps = GeoEnergyIO.read_restart(path; actnum = grid["ACTNUM"], units = :native, kw...)
    nc = count(!iszero, grid["ACTNUM"])
    states = Dict{Symbol, Any}[]
    for step in steps
        h = step["INTEHEAD"]
        Tuple(h[9:11]) == grid_dims(grid) && h[12] == nc || throw(ArgumentError("Restart grid does not match case"))
        push!(states, state_from_restart(step, case, nc))
    end
    return ifelse(extra_out, (states, steps), states)
end
