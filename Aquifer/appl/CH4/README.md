# CH4 3D Cylindrical Central-Well Case

TwoPNC simulation of underground CH4 storage with a **true 3D cylindrical domain** (r, θ, z) and a central well modeled as a volumetric source/sink.

## Domain

- **Grid type:** `Dune::UGGrid<3>` created by `CakeGridManager` (DuMux built-in cylindrical mesh generator)
- **Geometry:** 3D cylindrical mesh, axes: `r` (radial), `θ` (angular), `z` (axial/vertical)
- **Radial grading:** Fine cells near the well (r = 0.2 m), progressively coarser outward to r = 500 m. Grading controlled by `Grid.Grading0`.
- **Angular resolution:** 8 sectors (full 360°). Increase `Cells1` for finer angular resolution.
- **Vertical spacing:** Uniform, 4 layers over 22 m height.
- **Cell types:** Hexahedral cells away from center, prism (wedge) cells at the center axis.

Grid parameters in `params.input`:
```
[Grid]
Radial0 = 0.2 5 50 500    # radial breakpoints [m]
Angular1 = 0 360           # full circle [deg]
Axial2 = 0 22              # vertical extent [m]
Cells0 = 4 6 6             # cells per radial zone
Cells1 = 8                 # angular sectors
Cells2 = 4                 # vertical layers
Grading0 = 1.0 1.3 1.3     # radial grading factor (>1 = cells grow outward)
```

## Central Well

- **Implementation:** Volumetric source/sink term in `source()` (not a boundary condition).
- **Well identification:** Cells with `r = sqrt(x² + y²) ≤ Well.Radius` and `z ∈ [Well.PerfBottom, Well.PerfTop]`.
- **Rate distribution:** Total well rate (`InjectionRateOp` / `ProductionRate`, mol/s) is uniformly distributed over `wellSupportVolume_` so each well cell gets rate per bulk volume `[mol/(m³·s)]`.
- **Injection:** Positive source with prescribed composition (`HydrogenInjectionConcentration` for H₂, remainder CO₂).
- **Production:** Negative source proportional to local composition (mobility-weighted).
- **Schedule:** Same injection/idle/production cycle semantics as the original 2D case.

## Boundary Conditions

- **Outer radial boundary (r = R_outer):** Dirichlet — fixed initial pressure.
- **Top / Bottom / Inner well surface:** Neumann — no flow (well handled via source term).

## Key Differences from H2 (2D Axisymmetric)

| Aspect | H2 (2D axisym) | CH4 (3D cylindrical) |
|--------|----------------|----------------------|
| Grid | 2D YaspGrid + RotationalExtrusion | 3D UGGrid via CakeGridManager |
| Mesh in ParaView | Rectangle | Cylinder |
| Well | Neumann BC on left boundary | Volumetric source/sink at r ≈ 0 |
| Permeability | Kx, Ky | Kx = Ky (horizontal), Kz = Kx/10 |
| Linear solver | ILU-BiCGSTAB | AMG-BiCGSTAB (required for UGGrid) |

## Dependencies

- **dune-uggrid** must be installed and found by CMake (provides `Dune::UGGrid`).
- The `CakeGridManager` is a built-in DuMux utility (`dumux/io/grid/cakegridmanager.hh`).

## Building

```bash
make appl_2pnc_box_CH4     # Box discretization (recommended)
make appl_2pnc_tpfa_CH4    # CCTpfa discretization
```

## Smoke Test

```bash
./appl_2pnc_box_CH4 params.input -TimeLoop.TEnd 1
```

## Assumptions & Limitations

- Prism cells at the grid center axis may cause slower linear solver convergence; AMG preconditioner is used.
- Well rate is distributed uniformly over support volume; no Peaceman-type well index correction.
- Dispersion length scale is estimated as an average cell size (approximate for graded meshes).
