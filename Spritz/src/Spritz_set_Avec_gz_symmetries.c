#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <stdlib.h>
#include <string.h>

// Is z-reflection symmetry actually active? Cached once per run: this
// function is called every substep (see schedule.ccl), so it must not
// re-query/re-parse a parameter on every call. ReflectionSymmetry::reflection_z
// is private to that thorn (no visibility keyword in its param.ccl), so it
// can't be imported via "shares:" -- read it via CCTK_ParameterValString
// instead, the same cross-thorn-read technique used for
// IllinoisGRMHD_symmetry_gzs_is_vertex_centered() in
// GRHayLET/IllinoisGRMHD/src/symmetry_gzs_staggered.c.
static int Spritz_avec_gz_symmetries_is_active(void) {
  static int cached = -1; // -1 == not yet queried
  if (cached < 0) {
    if (CCTK_IsThornActive("ReflectionSymmetry")) {
      char *reflection_z = CCTK_ParameterValString("reflection_z", "ReflectionSymmetry");
      cached = (reflection_z != NULL) && (strcmp(reflection_z, "yes") == 0);
      free(reflection_z);
    } else {
      cached = 0; // no symmetry thorn active -> nothing to mirror
    }
  }
  return cached;
}

/*
 * Direct port of IllinoisGRMHD_set_symmetry_gzs_staggered
 * (GRHayLET/IllinoisGRMHD/src/symmetry_gzs_staggered.c), adapted for Spritz's
 * Avec_x/Avec_y, which are staggered along z:
 *   Avec_x(i, j+1/2, k+1/2)
 *   Avec_y(i+1/2, j, k+1/2)
 *   Avec_z(i+1/2, j+1/2, k)   -- NOT z-staggered
 * (per interface.ccl's own comments and STAGGER011/101/110 Prolongation tags).
 *
 * ReflectionSymmetry's generic ghostzone fill cannot handle this: for a
 * "VECTOR"+"edge"-staggered group it derives which axis a component belongs
 * to via tcomponent = (vi-firstvar)/vectorlength, which only works if all
 * three components live in one 3-member Cactus group. Avec_x/Avec_y/Avec_z
 * are each their own one-component group instead (required so each can
 * carry its own distinct Prolongation="STAGGERxyz" tag -- Carpet's
 * Prolongation tag is per-group). So every one of them looks like
 * "component 0" to ReflectionSymmetry, which then silently treats them as
 * cell-centered when mirroring across the z=0 plane -- corrupting Avec_x's
 * and Avec_y's z-ghostzones, which feeds directly into Bx = dAz/dy - dAy/dz
 * and By = dAx/dz - dAz/dx.
 */
static void Spritz_mirror_gz_staggered_z(
      const cGH *cctkGH,
      const CCTK_REAL *restrict z,
      CCTK_REAL *restrict gridfunc,
      const CCTK_REAL parity_z,
      const int stagger_z) {

  const int lsh[3] = {cctkGH->cctk_lsh[0], cctkGH->cctk_lsh[1], cctkGH->cctk_lsh[2]};

  const CCTK_REAL dz = z[CCTK_GFINDEX3D(cctkGH,0,0,1)] - z[CCTK_GFINDEX3D(cctkGH,0,0,0)];
  const CCTK_REAL z_offset = dz*0.5*stagger_z;

  // Only components whose LOCAL patch actually touches the z=0 plane need
  // any mirroring; for every other (non-boundary) component this loop
  // condition is false at k=0, so num_gzs stays 0 and the fix below is a
  // no-op, same as IllinoisGRMHD_set_symmetry_gzs_staggered.
  int num_gzs=0;
  while( (z[CCTK_GFINDEX3D(cctkGH,0,0,num_gzs)]+z_offset) < -dz*0.1 && num_gzs<lsh[2] ) num_gzs++;
  if(num_gzs*2>=lsh[2]) CCTK_VERROR("ERROR in Spritz_mirror_gz_staggered_z.c");

#pragma omp parallel for
  for(int k=0; k<num_gzs; k++) {
    for(int j=0; j<lsh[1]; j++) {
      for(int i=0; i<lsh[0]; i++) {
        const int index_inside_sym_gz  = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int index_outside_sym_gz = CCTK_GFINDEX3D(cctkGH,i,j,(num_gzs*2-stagger_z)-k);
        gridfunc[index_inside_sym_gz] = parity_z*gridfunc[index_outside_sym_gz];
      }
    }
  }
}

/*
 * Parities cross-checked against Spritz_InitSymBound (Spritz_Boundaries.F90):
 *   Avec_x: sym=(-1,1,1)  -> z-parity = +1 (even, same as vx/scon_x)
 *   Avec_y: sym=(1,-1,1)  -> z-parity = +1 (even, same as vy/scon_y)
 *   Avec_z: sym=(1,1,-1)  -> z-parity = -1 (odd,  same as vz/scon_z)
 */
void Spritz_set_Avec_gz_symmetries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // Nothing to mirror if z-reflection symmetry isn't actually turned on --
  // this was previously missing, so this function ran unconditionally on
  // every substep even in full-domain, Symmetry="none"-style runs, where it
  // could silently corrupt Avec_x/Avec_y near z=0 (or abort outright, if some
  // component's whole local extent happened to look "below" the nonexistent
  // plane). Mirrors IllinoisGRMHD_set_symmetry_gzs_staggered's own
  // `if(CCTK_EQUALS(Symmetry, "none")) return;` guard.
  if (!Spritz_avec_gz_symmetries_is_active()) return;

  Spritz_mirror_gz_staggered_z(cctkGH, z, Avec_x, +1.0, 1);
  Spritz_mirror_gz_staggered_z(cctkGH, z, Avec_y, +1.0, 1);
  // Avec_z isn't z-staggered, so the generic ReflectionSymmetry fill is
  // already correct for it -- included anyway for defense in depth, at
  // negligible cost, using the unstaggered (stagger_z=0) formula.
  Spritz_mirror_gz_staggered_z(cctkGH, z, Avec_z, -1.0, 0);
}
