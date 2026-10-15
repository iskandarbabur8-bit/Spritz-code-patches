# Spritz patches: z-reflection symmetry + cold_tabulated (QHC21A) EOS support

These are diffs against the Spritz GRMHD thorn (Bruno Giacomazzo et al.,
part of the SpritzCode Cactus arrangement), NOT a copy of Spritz itself.
Spritz's own source carries no explicit license, so it is not
redistributed here -- only the changes are. To apply, obtain your own
copy of the Spritz thorn and, from the arrangement's `Spritz/` directory:

```
patch -p1 < patches/0002-noble-con2prim-cold_tabulated-eps-fix.patch
patch -p1 < patches/0003-prim2con-eps-roundtrip-fix.patch
patch -p1 < patches/0004-minima-eps-roundtrip-fix.patch
patch -p1 < patches/0005-schedule.ccl.patch
patch -p1 < patches/0006-make.code.defn.patch
cp patches/0001-Spritz_set_Avec_gz_symmetries.c src/Spritz_set_Avec_gz_symmetries.c
patch -p1 < patches/0007-interface.ccl.patch
patch -p1 < patches/0008-param.ccl.patch
patch -p1 < patches/0009-schedule.ccl-coldtab-storage.patch
patch -p1 < patches/0010-noble-con2prim-matched-polytrope.patch
patch -p1 < patches/0011-con2prim-mhd-noble-call-sites.patch
patch -p1 < patches/0012-minima-coldtab-calibration.patch
patch -p1 < patches/0013-avec-gz-symmetries-guard.patch
```

Patches 0007-0012 apply on top of 0002-0006 (sequential, not independent -- in
particular 0010 modifies the same function 0002 already patched, and 0012
modifies the same routine 0004 already patched). 0013 applies on top of the
`0001` file copy, not the 0002-0012 sequence.

## What's here

1. **`0001-Spritz_set_Avec_gz_symmetries.c`** (new file, original work, full
   content) -- fixes z-reflection symmetry corrupting Bx/By near
   refinement boundaries. `Avec_x`/`Avec_y`/`Avec_z` are each their own
   one-component Cactus group (required so each can carry its own
   `Prolongation="STAGGERxyz"` tag), so ReflectionSymmetry's generic
   ghostzone fill can't tell them apart and silently treats the two
   z-staggered ones (`Avec_x`, `Avec_y`) as cell-centered when mirroring
   across z=0. This corrupts their z-ghostzones and, through
   `Bx = dAz/dy - dAy/dz` / `By = dAx/dz - dAz/dx`, Bx/By as well. This
   routine does the mirroring correctly and is scheduled after every SYNC
   of Avec (Con2Prim, Boundaries, and once per Runge-Kutta substep in the
   RHS group -- the last one is the fix's critical point, since Carpet
   refreshes ghostzones between substeps too).

2. **`0002-noble-con2prim-cold_tabulated-eps-fix.patch`** --
   `Spritz_Con2PrimM_Noble_pt_EOSOmni.c`'s Noble con2prim routine calls
   `EOS_Omni_EpsFromPress` in its atmosphere-reset branch (on
   Newton-Raphson failure); that function has no `case(5)`
   (`cold_tabulated`) and hard-aborts the run. `EOS_Omni_press`'s own
   `case(5)` already computes both press and eps in one call when given
   `keytemp=1`, so the redundant second call is replaced with a local
   `keytemp=1` override on the existing first call.

3. **`0003-prim2con-eps-roundtrip-fix.patch`** -- `Spritz_Prim2Con.F90`
   unconditionally computed `press = Press(rho, eps)` then
   `eps = EpsFromPress(rho, press)`, regardless of con2prim method. For
   any self-consistent EOS this round trip is a no-op (the second step
   must return exactly the eps the first step's press was computed
   from), so the redundant, cold_tabulated-incompatible second call is
   dropped and eps is left unchanged.

4. **`0004-minima-eps-roundtrip-fix.patch`** -- `Spritz_Minima.cc`'s
   atmosphere `Cs2_min` setup (runs once at startup, also regardless of
   con2prim method) already calls `EOS_Omni_press` with `keytemp=1`
   first, which already writes the cold-curve eps directly into the
   output buffer -- the subsequent `EpsFromPress` call was pure
   redundancy on top of that (two occurrences, same pattern).

5. **`0005-schedule.ccl.patch`** -- schedules
   `Spritz_set_Avec_gz_symmetries` at the three points needed: after
   `Con2Prim`, after `Spritz_Boundaries`, and -- critically -- at the
   start of the `SpritzRHS` group (which runs once per Runge-Kutta
   substep, not once per full timestep).

6. **`0006-make.code.defn.patch`** -- registers the new source file in
   the thorn's build list.

7. **`0007-interface.ccl.patch`** / **`0008-param.ccl.patch`** /
   **`0009-schedule.ccl-coldtab-storage.patch`** -- adds a new
   `Spritz_coldtab_calib` grid-scalar group (`rhomin`, `press_lo`,
   `eps_lo`, `gamma_lo`) and a `shares: EOS_Omni` / `USES STRING
   coldeos_table_name` param, for the matched-polytrope fix below.

8. **`0010-noble-con2prim-matched-polytrope.patch`** -- Spritz's EOS
   engine (EOS_Omni) extrapolates the cold_tabulated EOS below the
   table's own `rho_min` with a **hardcoded Gamma=2** power law,
   regardless of the table's actual local slope (confirmed 1.3853 for
   the QHC21A table in use here), unlike IllinoisGRMHD/WhiskyTHC, which
   both use the table's real local Gamma. The Noble con2prim routine's
   atmosphere-reset branch (patched by 0002 to stop crashing) now uses a
   pure power-law "matched polytrope" calibrated against the table
   itself (see 0012) instead: `eps = eps_lo*ratio^(Gamma_lo-1)`,
   `P = P_lo*ratio^Gamma_lo`, continuous with the table at the matching
   point, computed with **zero EOS_Omni calls** during the con2prim
   failure itself (calibration happens once, at startup -- see 0012).
   An earlier additive form (integrating `deps=-P d(1/rho)` from the
   matching point) was tried first and rejected after it was found to go
   *negative* ~200x below the table's `rho_min`; the multiplicative form
   used here is a pure power law and stays positive by construction --
   matches GRHayL's own `NRPyEOS_hybrid_coldtab_compute_P_cold_and_eps_cold.c`,
   which documents the same finding.

9. **`0011-con2prim-mhd-noble-call-sites.patch`** -- threads the new
   calibration constants through both Fortran call sites of the Noble
   routine (`Spritz_Con2Prim_Noble`'s C signature gained 4 new trailing
   arguments).

10. **`0012-minima-coldtab-calibration.patch`** -- `Spritz_Minima.cc`'s
    one-time atmosphere setup now also calibrates the matched-polytrope
    constants for cold_tabulated: reads the table's own `RhoMin` from its
    header, then calls the existing `EOS_Omni_press` API twice (at
    `rho_min` and `1.1*rho_min`, a WhiskyTHC-style two-point secant) to
    recover the table's real local Gamma there -- the only EOS_Omni
    calls involved, made once at startup rather than during a con2prim
    failure. Also replaces this file's own `Cs2_min` computation for
    cold_tabulated with the same analytic formula, dropping its
    `EOS_Omni_DPressByDEps`/`DPressByDRho` calls too (`dP/deps=0` below
    `rho_min`, confirmed against `EOS_Omni_EOS_DPressByDEps`'s own
    `case(5)`; `dP/drho` is the power law's analytic derivative).

11. **`0013-avec-gz-symmetries-guard.patch`** -- `0001`
    (`Spritz_set_Avec_gz_symmetries`) had no guard at all: it ran
    unconditionally on every call (every RK substep, per 0005's
    scheduling), regardless of whether `ReflectionSymmetry` was even
    active. In a full-domain (`Symmetry="none"`-style) run this would
    still mirror `Avec_x`/`Avec_y` across `z=0` every substep, corrupting
    real (non-symmetry-plane) physics there via
    `Bx = dAz/dy - dAy/dz` / `By = dAx/dz - dAz/dx`. Adds a
    `Spritz_avec_gz_symmetries_is_active()` check (cached after first
    call, since this runs every substep) that reads
    `ReflectionSymmetry::reflection_z` via `CCTK_ParameterValString` --
    that parameter is `PRIVATE` to `ReflectionSymmetry` (no
    `RESTRICTED:`/`GLOBAL:` header in its `param.ccl`), so it can't be
    pulled in with `shares:` and has to be read via this cross-thorn
    introspection call instead, the same technique
    `IllinoisGRMHD_symmetry_gzs_is_vertex_centered()` uses for the
    analogous check. Returns immediately (no mirroring at all) if
    `ReflectionSymmetry` isn't active or `reflection_z="no"`.

## Verification

- z-reflection symmetry: validated against a non-symmetric reference run
  (same physical setup) at matching iteration/time -- `ADMBase::lapse`,
  `HydroBase::rho`, `ML_ADMConstraints::ML_Ham` all agree to 5-9
  significant figures.
- Noble + cold_tabulated (QHC21A): verified past the atmosphere-minima
  startup crash that all three `EpsFromPress` fixes address. Exercising
  the con2prim-failure path specifically (patch 2) still pending -- zero
  con2prim failures were observed in an equivalent GRHayL/IllinoisGRMHD
  production run with the same physical setup, so that branch may not be
  hit in a short smoke test.
- Matched-polytrope calibration (patches 7-12): startup log confirmed
  `rho_min=2.0e-9, P_lo=8.897e-13, eps_lo=8.805e-05, Gamma_lo=1.3175`,
  matching a hand-computed cross-check exactly and landing within 0.1% of
  GRHayL's independently-precomputed `Gamma_lo=1.3192` for its own
  (differently-generated) copy of the same QHC21A table. A full smoke
  test ran cleanly through 26 atmosphere resets with zero crashes and
  zero `Cs2_min is negative` aborts.
- `reflection_z` guard (patch 13): parameter name, type (`BOOLEAN`), and
  `PRIVATE` visibility confirmed directly against
  `CactusNumerical/ReflectionSymmetry/param.ccl`; `CCTK_ParameterValString`/
  `CCTK_IsThornActive` signatures confirmed against the flesh headers
  (`cctk_Parameter.h`, `cctk_ActiveThorns.h`). Not yet exercised with an
  actual `Symmetry="none"` run -- the equatorial-symmetry runs used for the
  z-reflection verification above would not have caught this gap, since
  `reflection_z="yes"` in all of them.
