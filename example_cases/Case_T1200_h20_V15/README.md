# Methane pyrolysis — carbon transfer tutorial

This example uses `reactingParcelFoamCarbonTransfer` with the `periodicCellCarbonInjection` model in OpenFOAM.com v2312. Carbon formed by the gas-phase reaction is transferred to Lagrangian parcels every 0.1 s, with the corresponding carbon removed from the gas through source terms.

The case corresponds to the **Figure 7** condition in the paper by Zhang et al. (2026): furnace set-point temperature **1200 °C**, methane inlet flow rate **15 SLPM**, and external heat-transfer coefficient **h_ext = 20 W/(m² K)**. The `tubeWall` temperature profile in `0.org/T` was taken from the corresponding conjugate heat-transfer (CHT) simulation. It contains one temperature (K) for each of the 14,400 wall faces; `h_ext` is used in that earlier CHT calculation, not as a boundary condition in this case.

## Settings specific to this model

- `constant/reactingCloud1Properties` — injection times, staggered groups, carbon thresholds, parcel size, and gas–particle coupling.
- `constant/carbonToParticleProperties` — matching gas-side transfer times, particle energy settings, and an optional mass check (off by default).
- `constant/radiationProperties` — gas radiation and carbon-particle absorption/emission.
- `0.org/T` — nonuniform `tubeWall` temperatures from the CHT simulation.

`firstInjectionTime` and `injectionInterval` set the first injection and the interval for each group. `staggerGroups 20` spreads injections across the 0.1 s cycle. The carbon thresholds are zero, so cells with positive carbon mass can inject. `SOI`, `massTotal`, and `nParticle` are required by the OpenFOAM injection-model interface; the custom model sets its own injection times and parcel weights.

The simulation uses a fixed time step of `0.0001 s`. Injections start at `0.1 s`, followed by the other groups at 0.005 s intervals. The case includes chemistry files, starting fields, and a compressed mesh (`constant/polyMesh.tar.gz`), which `Allrun` unpacks if needed.
