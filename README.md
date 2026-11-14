# Carbon transfer for methane pyrolysis in OpenFOAM v2312

This is a modified OpenFOAM.com v2312 solver for methane pyrolysis. It allows carbon formed by the gas-phase reaction to be tracked as Lagrangian particles.

Carbon is first produced by the gas chemistry as a species named `C`. At a chosen time interval, the injection model converts the carbon in each selected cell into a particle parcel. The parcel starts at a random position inside the cell and uses the local gas velocity and temperature. It represents the carbon mass transferred from that cell.

The gas loses this carbon through source terms in the mass, species, momentum and energy equations. The solver then updates the gas fields in the usual way. Once injected, parcels move through the reactor and exchange momentum and heat with the gas. Particle absorption and emission can also be included in the radiation calculation.

This code and the accompanying tutorial case are associated with our published paper:
Zhang, X., Sun, Z., Chin, R. C., & Nathan, G. J. (2026). An Eulerian–Lagrangian CFD framework for methane pyrolysis environments incorporating reaction-generated carbon particles with radiative heating. Applications in Energy and Combustion Science, 100573. 

**DOI:** [10.1016/j.jaecs.2026.100573](https://doi.org/10.1016/j.jaecs.2026.100573)

## Included code

- `periodicCellCarbonInjection` — creates parcels at regular intervals, without a positions file. Cells with little carbon can be skipped using optional thresholds.
- `reactingParcelFoamCarbonTransfer` — a modified reacting parcel solver with gas-to-particle carbon transfer.

## Compile

This code is for **OpenFOAM.com v2312** and needs a working `wmake` environment and the v2312 source files.

After loading OpenFOAM v2312, run from this folder:

```bash
./Allwmake
```

The build writes the modified `liblagrangianIntermediate.so` to `$FOAM_USER_LIBBIN` and the two solvers to `$FOAM_USER_APPBIN`. It does not change the system OpenFOAM installation, but it may replace files with the same names in your user directories.

To remove the local build files, run `./Allwclean`.

## Set up a case

A complete methane-pyrolysis example is provided separately. It includes the dictionaries, mesh, initial fields and chemistry files.

The main settings are in:

| File | What to set |
| --- | --- |
| `constant/reactingCloud1Properties` | `periodicCellCarbonInjection`, particle properties, tracking and heat-transfer models; set `radiation on` when using particle radiation. |
| `constant/carbonToParticleProperties` | `carbonSpecie C`, first transfer time, transfer interval, and particle energy reference settings. |
| `constant/radiationProperties` | Gas radiation and, when needed, `cloudAbsorptionEmission` for `reactingCloud1`. |
| `system/controlDict` | Solver name, time step, run time and output interval. |

Use the same `firstInjectionTime` and `firstTransferTime`, and the same `injectionInterval` and `transferInterval`. Injection and gas-side carbon removal then occur in the same time step. If using `staggerGroups`, use the same value in both files. Staggered injection has been tested with a fixed time step (`deltaT <= injectionInterval/staggerGroups`); variable time steps have not yet been validated.

With `parcelsPerCell 1`, each selected cell creates one computational parcel per event. Its number of physical particles is set by the mass of carbon in that cell. For the unfiltered method, set `minCarbonMassFraction`, `minCarbonMass` and `minParticlesPerParcel` to zero. The particle diameter and material density are set separately in the cloud dictionary.

For an example that starts with a developed methane flow, map the gas fields into `0/` before running. The supplied example does not contain a precomputed developed field.

## Notes

- The gas-phase `C` field is a temporary representation of the carbon product. It is not intended to represent gaseous atomic carbon.
- The transfer uses finite-volume source terms for mass, carbon species, momentum and energy. It does not directly overwrite the gas fields at the injection event.
- Carbon transfer has been run stably through repeated injection events in an example with mapped initial fields. A full check of the solved gas/particle mass and energy balances, including the thermodynamic reference energies, is still needed.
- This implementation does not add particle growth, agglomeration or wall deposits. Other particle-force and collision options depend on the case settings.
