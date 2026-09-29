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
```

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
