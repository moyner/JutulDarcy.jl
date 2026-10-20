const SCALAR_FIELDS = (
    (:Pressure, "PRESSURE", :pressure),
    (:Temperature, "TEMP", :temperature),
    (:Rs, "RS", :u_rs),
    (:Rv, "RV", :u_rv),
    (:FluidVolume, "RPORV", :volume),
    (:BulkVolume, "JU_BVOL", :volume),
    (:TotalThermalEnergy, "JU_THERM", :energy),
    (:TotalEnergy, "JU_ENERG", :energy),
    (:PotentialEnergy, "JU_POTEN", :energy),
    (:RockInternalEnergy, "JU_ROCKU", :specific_energy),
    (:RockDensity, "JU_RHO_R", :density),
    (:RockHeatCapacity, "JU_CP_R", :mass_heat_capacity)
)

# Phase keyword tuples use water, oil, gas ordering.
const PHASE_FIELDS = (
    (:Saturations, ("SWAT", "SOIL", "SGAS"), :identity),
    (:PhaseMassDensities, ("DENW", "DENO", "DENG"), :density),
    (:PhaseViscosities, ("VWAT", "VOIL", "VGAS"), :viscosity),
    (:RelativePermeabilities, ("KRW", "KRO", "KRG"), :identity),
    (:ShrinkageFactors, ("1OVERBW", "1OVERBO", "1OVERBG"), :inverse_fvf),
    (:PhasePressures, ("WAT_PRES", "OIL_PRES", "GAS_PRES"), :pressure),
    (:MaxSaturations, ("SWMAX", "SOMAX", "SGMAX"), :identity),
    (:FluidEnthalpy, ("JU_H_W", "JU_H_O", "JU_H_G"), :specific_energy),
    (:FluidInternalEnergy, ("JU_U_W", "JU_U_O", "JU_U_G"), :specific_energy),
    (:PhaseMobilities, ("JU_MOBW", "JU_MOBO", "JU_MOBG"), :mobility),
    (:PhaseMassMobilities, ("JU_MMBW", "JU_MMBO", "JU_MMBG"), :mass_mobility)
)

const INDEXED_FIELDS = (
    (:OverallMoleFractions, "ZMF", :identity),
    (:TotalMasses, "JU_M", :mass),
    (:ComponentHeatCapacity, "JU_CP", :mass_heat_capacity),
    (:CapillaryPressure, "JU_PC", :pressure),
    (:TracerConcentrations, "JU_TR", :identity),
    (:TracerMasses, "JU_TM", :mass)
)

function restart_factor(quantity, units, phase)
    system = unit_system(units)
    specific_energy = system.mass_heat_capacity*system.absolute_temperature_numeric
    if quantity == :specific_energy
        return specific_energy
    elseif quantity == :energy
        return specific_energy*system.mass
    elseif quantity == :mobility
        return inv(system.viscosity)
    elseif quantity == :mass_mobility
        return system.density/system.viscosity
    elseif quantity == :inverse_fvf
        if phase.mask == 4
            return inv(system.gas_formation_volume_factor)
        end
        return inv(system.liquid_formation_volume_factor)
    end
    return getproperty(system, quantity)
end

function restart_values(values, quantity, units; to_si, phase = nothing)
    if quantity == :identity
        return copy(values)
    elseif quantity == :temperature
        return convert_array("TEMP", values, units; to_si)
    end
    factor = restart_factor(quantity, units, phase)
    return values .* ifelse(to_si, factor, inv(factor))
end

function cell_values(values, nc, name)
    if values isa Number
        return fill(values, nc)
    end
    length(values) == nc || throw(ArgumentError("$name must have $nc cell values"))
    return vec(values)
end

function cell_matrix(values, nc, name)
    if values isa AbstractVector
        values = reshape(values, 1, :)
    end
    ndims(values) == 2 && size(values, 2) == nc ||
        throw(ArgumentError("$name must have one column per active cell ($nc)"))
    return values
end

function component_molar_masses(case)
    system = reservoir_model(case).system
    system isa CompositionalSystem || throw(ArgumentError("Phase compositions require a compositional model"))
    return MultiComponentFlash.molar_masses(system.equation_of_state)
end

function convert_fractions(fractions, masses; to_mass)
    size(fractions, 1) == length(masses) || throw(ArgumentError("Composition rows differ from model components"))
    weights = ifelse(to_mass, masses, inv.(masses))
    converted = fractions .* weights
    for col in eachcol(converted)
        total = sum(col)
        if total > 0
            col ./= total
        end
    end
    return converted
end

function state_keywords(state, case, nc; units = :metric)
    units = Symbol(lowercase(string(units)))
    state = reservoir_part(state)
    phases = model_phases(case)
    keywords = BinaryKeyword[]
    function add(name, values, quantity; phase = nothing)
        values = cell_values(values, nc, name)
        values = restart_values(values, quantity, units; to_si = false, phase)
        push!(keywords, BinaryKeyword(name, Float32.(values)))
    end
    for (key, name, quantity) in SCALAR_FIELDS
        if haskey(state, key)
            add(name, state[key], quantity)
        end
    end
    for (key, names, quantity) in PHASE_FIELDS
        if !haskey(state, key)
            continue
        end
        values = cell_matrix(state[key], nc, key)
        size(values, 1) == length(phases) || throw(ArgumentError("$key rows differ from model phases"))
        for (i, phase) in enumerate(phases)
            add(names[phase.slot], values[i, :], quantity; phase)
        end
    end
    if !haskey(state, :Saturations) && haskey(state, :ImmiscibleSaturation)
        system = reservoir_model(case).system
        if system isa CompositionalSystem && has_other_phase(system)
            other = first(phase_indices(system))
            add(phases[other].keyword, state[:ImmiscibleSaturation], :identity)
        end
    end
    for (key, prefix, quantity) in INDEXED_FIELDS
        if haskey(state, key)
            values = cell_matrix(state[key], nc, key)
            for i in axes(values, 1)
                add("$prefix$i", values[i, :], quantity)
            end
        end
    end
    for (key, prefix) in ((:LiquidMassFractions, "XMF"), (:VaporMassFractions, "YMF"))
        if haskey(state, key)
            fractions = cell_matrix(state[key], nc, key)
            fractions = convert_fractions(fractions, component_molar_masses(case); to_mass = false)
            for i in axes(fractions, 1)
                add("$prefix$i", fractions[i, :], :identity)
            end
        end
    end
    return keywords
end

function indexed_array(step, prefix, quantity, units, nc)
    indices = Int[]
    for name in keys(step)
        if startswith(name, prefix)
            index = tryparse(Int, name[length(prefix)+1:end])
            if !isnothing(index)
                push!(indices, index)
            end
        end
    end
    if isempty(indices)
        return nothing
    end
    sort!(indices) == collect(1:length(indices)) || throw(ArgumentError("$prefix indices must be contiguous from 1"))
    values = zeros(length(indices), nc)
    for i in indices
        name = "$prefix$i"
        values[i, :] = restart_values(cell_values(step[name], nc, name), quantity, units; to_si = true)
    end
    return values
end

function state_from_restart(step, case, nc)
    units = file_units(step)
    phases = model_phases(case)
    state = Dict{Symbol, Any}()
    for (key, name, quantity) in SCALAR_FIELDS
        if haskey(step, name)
            values = cell_values(step[name], nc, name)
            state[key] = restart_values(values, quantity, units; to_si = true)
        end
    end
    for (key, names, quantity) in PHASE_FIELDS
        present = [haskey(step, names[p.slot]) for p in phases]
        if !any(present) || (key != :Saturations && !all(present))
            continue
        end
        values = zeros(length(phases), nc)
        for (i, phase) in enumerate(phases)
            name = names[phase.slot]
            if present[i]
                values[i, :] = restart_values(cell_values(step[name], nc, name), quantity, units; to_si = true, phase)
            elseif name == "SOIL"
                values[i, :] .= 1
                for other in ("SWAT", "SGAS")
                    if haskey(step, other)
                        values[i, :] .-= cell_values(step[other], nc, other)
                    end
                end
            end
        end
        state[key] = values
    end
    for (key, prefix, quantity) in INDEXED_FIELDS
        values = indexed_array(step, prefix, quantity, units, nc)
        if !isnothing(values)
            state[key] = values
        end
    end
    for (key, prefix) in ((:LiquidMassFractions, "XMF"), (:VaporMassFractions, "YMF"))
        fractions = indexed_array(step, prefix, :identity, units, nc)
        if !isnothing(fractions)
            state[key] = convert_fractions(fractions, component_molar_masses(case); to_mass = true)
        end
    end
    system = reservoir_model(case).system
    if system isa CompositionalSystem && has_other_phase(system) && haskey(state, :Saturations)
        other = first(phase_indices(system))
        state[:ImmiscibleSaturation] = copy(state[:Saturations][other, :])
    end
    return state
end
