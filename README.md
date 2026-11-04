# Carbon transfer for methane pyrolysis in OpenFOAM v2312

This is a modified OpenFOAM.com v2312 solver for methane pyrolysis. It allows carbon formed by the gas-phase reaction to be tracked as Lagrangian particles.

Carbon is first produced by the gas chemistry as a species named `C`. At a chosen time interval, the injection model converts the carbon in each selected cell into a particle parcel. The parcel starts at a random position inside the cell and uses the local gas velocity and temperature. It represents the carbon mass transferred from that cell.

The gas loses this carbon through source terms in the mass, species, momentum and energy equations. The solver then updates the gas fields in the usual way. Once injected, parcels move through the reactor and exchange momentum and heat with the gas. Particle absorption and emission can also be included in the radiation calculation.

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

## Run a case

A complete methane-pyrolysis example is provided as a separate download. The main settings are in:

- `constant/reactingCloud1Properties` — injection interval, particle size and properties, tracking, heat transfer and particle radiation.
- `constant/carbonToParticleProperties` — carbon species and gas-side transfer interval.
- `constant/radiationProperties` — gas radiation and particle radiation coupling.
- `system/controlDict` — solver, time step and run time.

Set the injection and gas-side transfer to the **same time and interval**. For one parcel per cell, use `parcelsPerCell 1`. Set the carbon thresholds to zero if all carbon is to be transferred. No cell-centre positions file is needed.

For particle radiation, enable `radiation on` in the cloud settings and include `cloudAbsorptionEmission` in the radiation model.

If your case uses a previously developed gas field, map it into the `0/` folder before running. The example does not include a mapped field.

## Notes

The `C` species is used to hold the carbon product before transfer; it is not gaseous atomic carbon. This version (v1.0) does not include particle growth, agglomeration or carbon deposits on walls.

Repeated carbon injection has been tested in a methane-pyrolysis case. A separate check of the full gas–particle mass and energy balances is still needed.

Released under GPL-3.0-or-later. See `LICENSE`.
