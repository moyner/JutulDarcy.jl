using GeoEnergyIO, JutulDarcy, Jutul, Test, Dates
using MultiComponentFlash

@testset "StandardIO JutulDarcy states and results" begin
    g = UnstructuredMesh(CartesianMesh((3, 2, 1), (30.0, 20.0, 5.0)))
    permeability = repeat(1e-12.*[1.0, 0.0, 0.0, 4.0, 0.0, 9.0], 1, 6)
    domain = reservoir_domain(g; porosity = 0.25, permeability)
    # Reversed phase ordering detects accidental hard-coded saturation rows.
    sys = ImmiscibleSystem((LiquidPhase(), AqueousPhase()); reference_densities = [800.0, 1000.0])
    model, parameters = setup_reservoir_model(domain, sys; extra_out = true)
    state0 = setup_reservoir_state(model; Pressure = 1e7, Saturations = [0.7, 0.3])
    dt = [86400.0, 86400.0]
    case = JutulCase(model, dt; state0, parameters, start_date = DateTime(2020, 1, 1))
    states = [deepcopy(state0[:Reservoir]), deepcopy(state0[:Reservoir])]
    states[2][:Pressure] .*= 1.1
    for state in states
        state[:ShrinkageFactors] = repeat([0.8, 1.1], 1, 6)
        state[:MaxSaturations] = repeat([0.9, 0.5], 1, 6)
        state[:PhasePressures] = repeat([1e7, 0.99e7], 1, 6)
        state[:CapillaryPressure] = fill(-1e5, 1, 6)
        state[:BulkVolume] = fill(500.0, 6)
        state[:FluidVolume] = fill(125.0, 6)
        state[:RockDensity] = fill(2500.0, 6)
        state[:RockHeatCapacity] = fill(1000.0, 6)
        state[:ComponentHeatCapacity] = repeat([2000.0, 4184.0], 1, 6)
        state[:TracerConcentrations] = fill(0.01, 1, 6)
        state[:TracerMasses] = fill(10.0, 1, 6)
    end
    mktempdir() do dir
        base = joinpath(dir, "RESULT")
        @test read_egrid(write_egrid_jutuldarcy(case, base))["cartDims"] == [6, 1, 1]
        griddata = GeoEnergyIO.StandardIO.output.rectangular_grid([10.0, 10.0, 10.0], [10.0, 10.0], [5.0])
        structured_case = JutulCase(model, dt; state0, parameters,
            input_data = Dict("GRID" => griddata), start_date = case.start_date)
        @test read_egrid(write_egrid_jutuldarcy(structured_case, base))["cartDims"] == [3, 2, 1]
        for units in (:metric, :field, :lab)
            write_init_jutuldarcy(case, base; units)
            init = read_init(base; units = :si)
            @test init["PORO"] ≈ fill(0.25, 6)
            @test init["INTEHEAD"][15] == 3
            @test init["PERMX"] ≈ fill(1e-12, 6) rtol=1e-6
            @test init["PERMY"] ≈ fill(4e-12, 6) rtol=1e-6
            @test init["PERMZ"] ≈ fill(9e-12, 6) rtol=1e-6
            @test init["PORV"] ≈ fill(125.0, 6) rtol=1e-6
            for unified in (true, false)
                path = base*"_$(units)_$unified"
                write_restart_jutuldarcy(case, states, path; units, unified)
                back = read_restart_jutuldarcy(path, case)
                @test read_restart(path)[1]["days"] == 1.0
                @test length(back) == 2
                for (a, b) in zip(states, back)
                    @test a[:Pressure] ≈ b[:Pressure] rtol=1e-6
                    @test a[:Saturations] ≈ b[:Saturations] rtol=1e-6
                    for key in (:ShrinkageFactors, :MaxSaturations, :PhasePressures,
                            :CapillaryPressure, :BulkVolume, :FluidVolume, :RockDensity,
                            :RockHeatCapacity, :ComponentHeatCapacity, :TracerConcentrations,
                            :TracerMasses)
                        @test a[key] ≈ b[key] rtol=1e-6
                    end
                end
            end
        end
        # A real simulation exercises ReservoirSimResult and summary_result.
        result = simulate_reservoir(case; info_level = -1)
        paths = write_results_jutuldarcy(case, result, joinpath(dir, "SIM"))
        @test isfile(paths.init)
        @test isfile(paths.grid)
        @test isfile(paths.restart)
        smry = read_summary(joinpath(dir, "SIM"); format = :structured)
        @test smry["TIME"].seconds ≈ result.time
        @test smry["TIME"].start_date == case.start_date
        @test smry["DIMENS"] == Tuple(read_egrid(paths.grid)["cartDims"])
        @test length(read_restart_jutuldarcy(paths.restart, case)) == length(result.states)
        field_paths = write_summary_jutuldarcy(case, result, joinpath(dir, "FIELD"); units = :field)
        field = read_summary(field_paths.spec; format = :structured)
        @test field["UNIT_SYSTEM"] == "FIELD"
        @test field["VALUES"]["FIELD"]["FPR"] ≈ JutulDarcy.summary_result(case, result, :field)["VALUES"]["FIELD"]["FPR"] rtol=1e-6
        write_summary(result.summary, joinpath(dir, "DIRECT_RESULT"))
        direct = read_summary(joinpath(dir, "DIRECT_RESULT"); format = :structured)
        expected = JutulDarcy.summary_result(case, result, :metric)
        @test direct["VALUES"]["FIELD"]["FPR"] ≈ expected["VALUES"]["FIELD"]["FPR"] rtol=1e-6
    end
end

@testset "Compositional and thermal restart output" begin
    mesh = UnstructuredMesh(CartesianMesh((2, 1, 1), (20.0, 10.0, 5.0)))
    domain = reservoir_domain(mesh; porosity = 0.25, permeability = 1e-12)
    co2 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
    c10 = MolecularProperty(0.1420, 2.10e6, 617.7, 6.098e-4, 0.488)
    mixture = MultiComponentMixture([co2, c10]; names = ["CO2", "C10"])
    system = MultiPhaseCompositionalSystemLV(GenericCubicEOS(mixture),
        (VaporPhase(), AqueousPhase(), LiquidPhase()))
    model, parameters = setup_reservoir_model(domain, system; thermal = true,
        energy_formulation = :total, extra_outputs = true, extra_out = true)
    state0 = setup_reservoir_state(model; Pressure = 1e7, Temperature = 350.0,
        OverallMoleFractions = [0.2, 0.8], ImmiscibleSaturation = 0.15)
    case = JutulCase(model, [10.0, 10.0]; state0, parameters)
    result = simulate_reservoir(case; info_level = -1)
    fields = (:Pressure, :Temperature, :OverallMoleFractions, :LiquidMassFractions,
        :VaporMassFractions, :ImmiscibleSaturation, :Saturations, :PhaseMassDensities,
        :PhaseViscosities, :RelativePermeabilities, :PhaseMobilities,
        :PhaseMassMobilities, :FluidEnthalpy, :FluidInternalEnergy,
        :RockInternalEnergy, :TotalThermalEnergy, :TotalEnergy, :PotentialEnergy,
        :TotalMasses)
    mktempdir() do dir
        for units in (:metric, :field, :lab), unified in (true, false)
            path = joinpath(dir, "$(units)_$unified")
            write_restart_jutuldarcy(case, result, path; units, unified)
            restored, raw = read_restart_jutuldarcy(path, case; extra_out = true)
            for (original, back) in zip(result.states, restored), key in fields
                @test back[key] ≈ original[key] rtol=2e-6
            end
            expected_x = [f.liquid.mole_fractions[1] for f in result.states[1][:FlashResults]]
            expected_y = [f.vapor.mole_fractions[1] for f in result.states[1][:FlashResults]]
            @test raw[1]["XMF1"] ≈ expected_x rtol=2e-6
            @test raw[1]["YMF1"] ≈ expected_y rtol=2e-6
            if units == :metric
                @test raw[1]["JU_THERM"] ≈ result.states[1][:TotalThermalEnergy]./1000 rtol=2e-6
                @test raw[1]["JU_H_G"] ≈ result.states[1][:FluidEnthalpy][1, :]./1000 rtol=2e-6
            elseif units == :lab
                @test raw[1]["JU_THERM"] ≈ result.states[1][:TotalThermalEnergy] rtol=2e-6
            end
        end
    end
end
