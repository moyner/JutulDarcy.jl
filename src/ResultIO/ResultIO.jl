module ResultIO
    using Dates: DateTime
    import Jutul
    import GeoEnergyIO
    import MultiComponentFlash
    using GeoEnergyIO: BinaryKeyword, KeywordFile, write_binary
    using GeoEnergyIO.StandardIO: convert_array, output_path, unit_system, file_units
    using GeoEnergyIO.StandardIO.output: corner_point_data, fake_grid, grid_dims,
        standard_headers, restart_block
    using ..JutulDarcy: AqueousPhase, LiquidPhase, VaporPhase, ReservoirSimResult,
        get_phases, reservoir_model, summary_result, permeability_diagonal_rows,
        CompositionalSystem, has_other_phase, phase_indices

    include("restart_fields.jl")
    include("standard_io.jl")

    export write_egrid_jutuldarcy, write_init_jutuldarcy, write_restart_jutuldarcy
    export write_summary_jutuldarcy, write_results_jutuldarcy, read_restart_jutuldarcy
end
