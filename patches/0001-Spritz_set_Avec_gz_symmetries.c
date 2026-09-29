#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

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

  Spritz_mirror_gz_staggered_z(cctkGH, z, Avec_x, +1.0, 1);
  Spritz_mirror_gz_staggered_z(cctkGH, z, Avec_y, +1.0, 1);
  // Avec_z isn't z-staggered, so the generic ReflectionSymmetry fill is
  // already correct for it -- included anyway for defense in depth, at
  // negligible cost, using the unstaggered (stagger_z=0) formula.
  Spritz_mirror_gz_staggered_z(cctkGH, z, Avec_z, -1.0, 0);
}
