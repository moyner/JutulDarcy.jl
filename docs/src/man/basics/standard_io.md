```@meta
CurrentModule = JutulDarcy
```

# Industry standard format input and output

`JutulDarcy.ResultIO` reads and writes industry standard format files using
reservoir domains, state dictionaries and `summary_result`. Its public
functions are also exported by JutulDarcy.

## Exporting simulation results

```julia
using JutulDarcy, GeoEnergyIO

files = write_results_jutuldarcy(case, result, "RESULT"; units = :metric)
```

The combined export writes EGRID, INIT, SMSPEC/UNSMRY and reservoir UNRST
solution arrays. Use `units=:field` or `units=:lab` for other output systems,
and `unified=false` for separate Xnnnn/Snnnn files.

Individual files can be written with:

```julia
write_egrid_jutuldarcy(case, "RESULT")
write_init_jutuldarcy(case, "RESULT")
write_restart_jutuldarcy(case, result, "RESULT")
write_summary_jutuldarcy(case, result, "RESULT")
```

EGRID uses the grid stored in `case.input_data`, with ACTNUM rebuilt from the
processed mesh. Parsed COORD/ZCORN and separable DX/DY/DZ or DXV/DYV/DZV
geometry are accepted. Without input geometry, a synthetic one-metre cell row
with dimensions `(number_of_cells, 1, 1)` preserves active-cell ordering. The
grid dimensions and start date are kept consistent across all exported files.

INIT maps porosity, permeability, pore volumes, depths and region arrays from
the reservoir domain and parameters.

For state vectors, supply the report times explicitly when they differ from
the case timesteps:

```julia
write_restart_jutuldarcy(case, states, "RESULT"; times = report_times)
```

Summary dictionaries can also be written directly through GeoEnergyIO:

```julia
summary = JutulDarcy.summary_result(case, result, :field)
GeoEnergyIO.write_summary(summary, "RESULT")
GeoEnergyIO.write_summary(result.summary, "RESULT_SI")
```

The summary writer converts SI dictionaries such as `result.summary` to metric
by default. The case/result adapter uses the chosen output system through
`summary_result` and does not modify the simulation results.

```@docs
write_results_jutuldarcy
write_egrid_jutuldarcy
write_init_jutuldarcy
write_restart_jutuldarcy
write_summary_jutuldarcy
```

## Reading reservoir states

```julia
states = read_restart_jutuldarcy("RESULT", case)
state0 = setup_reservoir_state(case.model; states[end]...)
```

The reader checks the grid dimensions and active-cell count, converts physical
arrays to SI, and maps the file arrays back to JutulDarcy state names and phase
ordering. Use `steps=[1, 3]` to select report positions, or `extra_out=true` to
also retrieve the restart dictionaries and their metadata.

Restart export contains reservoir solution arrays for visualization and state
initialization. It does not contain the complete well, control and solver state
needed to resume another simulator.

```@docs
read_restart_jutuldarcy
```

## Restart fields

Restart writers include supported fields present in the state dictionaries.
Use `extra_outputs=true` or a list of variables in `setup_reservoir_model` to
save additional secondary variables. Phase rows follow the model's ordering.

| State field | File keywords |
|:--|:--|
| Pressure, temperature, black-oil ratios | `PRESSURE`, `TEMP`, `RS`, `RV` |
| Saturations | `SWAT`, `SOIL`, `SGAS` |
| Phase densities | `DENW`, `DENO`, `DENG` |
| Phase viscosities | `VWAT`, `VOIL`, `VGAS` |
| Relative permeabilities | `KRW`, `KRO`, `KRG` |
| Shrinkage factors | `1OVERBW`, `1OVERBO`, `1OVERBG` |
| Phase pressures | `WAT_PRES`, `OIL_PRES`, `GAS_PRES` |
| Maximum saturations | `SWMAX`, `SOMAX`, `SGMAX` |
| Overall mole fractions | `ZMF1`, `ZMF2`, ... |
| Liquid/vapor mass fractions | `XMF1`, ... / `YMF1`, ... (mole fractions) |
| Fluid volume | `RPORV` |
| Fluid enthalpy / internal energy | `JU_H_W/O/G` / `JU_U_W/O/G` |
| Thermal / total / potential energy | `JU_THERM` / `JU_ENERG` / `JU_POTEN` |
| Rock internal energy, density, heat capacity | `JU_ROCKU`, `JU_RHO_R`, `JU_CP_R` |
| Bulk volume | `JU_BVOL` |
| Phase mobility / mass mobility | `JU_MOBW/O/G` / `JU_MMBW/O/G` |
| Component masses / heat capacities | `JU_M1`, ... / `JU_CP1`, ... |
| Capillary-pressure rows | `JU_PC1`, ... |
| Tracer concentrations / masses | `JU_TR1`, ... / `JU_TM1`, ... |

The reader converts physical quantities to SI and restores phase and component
ordering. `XMF` and `YMF` are converted to mass fractions using the case's
molecular weights. For compositional models, `ImmiscibleSaturation` is restored
from the corresponding phase saturation. Missing oil saturation is inferred
from water and gas saturations.

See [OPM's component output](https://github.com/OPM/opm-simulators/blob/master/opm/simulators/flow/CompositionalContainer.cpp)
for the indexed mole-fraction arrays.

`JU_` keywords contain auxiliary arrays. Energies use kJ, BTU and J for metric,
field and laboratory files; specific energies use kJ/kg, BTU/lb and J/g.
Tracer concentrations retain their simulation units. `extra_out=true` returns
native file arrays alongside the SI state dictionaries.
