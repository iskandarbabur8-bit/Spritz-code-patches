#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define VERR_DEF_PARAMS __LINE__, __FILE__, CCTK_THORNSTRING

#define IDX(i,j,k) CCTK_GFINDEX3D(cctkGH,(i),(j),(k))

#define XMAX_OB_LINEAR_EXTRAP(FUNC,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imax,j,k)] = 2.0 * FUNC[IDX(imax-1,j,k)] - FUNC[IDX(imax-2,j,k)];
#define YMAX_OB_LINEAR_EXTRAP(FUNC,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmax,k)] = 2.0 * FUNC[IDX(i,jmax-1,k)] - FUNC[IDX(i,jmax-2,k)];
#define ZMAX_OB_LINEAR_EXTRAP(FUNC,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmax)] = 2.0 * FUNC[IDX(i,j,kmax-1)] - FUNC[IDX(i,j,kmax-2)];

#define XMIN_OB_LINEAR_EXTRAP(FUNC,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imin,j,k)] = 2.0 * FUNC[IDX(imin+1,j,k)] - FUNC[IDX(imin+2,j,k)];
#define YMIN_OB_LINEAR_EXTRAP(FUNC,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmin,k)] = 2.0 * FUNC[IDX(i,jmin+1,k)] - FUNC[IDX(i,jmin+2,k)];
#define ZMIN_OB_LINEAR_EXTRAP(FUNC,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmin)] = 2.0 * FUNC[IDX(i,j,kmin+1)] - FUNC[IDX(i,j,kmin+2)];

/*********************************************
 * Apply outer boundary conditions on A_{\mu}
 * Original routine from IllinoisGRMHD
 * by Etienne et al 2015
 ********************************************/
extern "C" void Spritz_outer_boundaries_on_A_mu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  //bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  bool Symmetry_none=true; //No Reflection Symmetry

  int levelnumber = GetRefinementLevel(cctkGH);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(CCTK_EQUALS(Spritz_Extr_OB_Avec_Psi,"no")) {
     //CCTK_INFO("Your choice was NOT to perform linear extrapolation of Avec and Psi at outer boundaries"); 
     return;
  }

  //CCTK_INFO("Your choice was to perform linear extrapolation of Avec and Psi at outer boundaries");

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2]) 
    CCTK_VError(VERR_DEF_PARAMS,"ERROR: Spritz outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {
    int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt; // for cctk_nghostzones==3, this goes {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
    int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;

    int imin=cctk_nghostzones[0]-which_bdry_pt-1; // for cctk_nghostzones==3, this goes {2,1,0}
    int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
    int kmin=cctk_nghostzones[2]-which_bdry_pt-1;

    if(Spritz_mhd_fluxCD_Lorenz != 0)
    {
      if(cctk_bbox[1]) { XMAX_OB_LINEAR_EXTRAP(Avec_x,imax); XMAX_OB_LINEAR_EXTRAP(Avec_y,imax); XMAX_OB_LINEAR_EXTRAP(Avec_z,imax); XMAX_OB_LINEAR_EXTRAP(psi_mhd,imax);}//  } 
      if(cctk_bbox[3]) { YMAX_OB_LINEAR_EXTRAP(Avec_x,jmax); YMAX_OB_LINEAR_EXTRAP(Avec_y,jmax); YMAX_OB_LINEAR_EXTRAP(Avec_z,jmax); YMAX_OB_LINEAR_EXTRAP(psi_mhd,jmax);}//  }
      if(cctk_bbox[5]) { ZMAX_OB_LINEAR_EXTRAP(Avec_x,kmax); ZMAX_OB_LINEAR_EXTRAP(Avec_y,kmax); ZMAX_OB_LINEAR_EXTRAP(Avec_z,kmax); ZMAX_OB_LINEAR_EXTRAP(psi_mhd,kmax);}//  }
      
      if(cctk_bbox[0]) {                    XMIN_OB_LINEAR_EXTRAP(Avec_x,imin); XMIN_OB_LINEAR_EXTRAP(Avec_y,imin); XMIN_OB_LINEAR_EXTRAP(Avec_z,imin); XMIN_OB_LINEAR_EXTRAP(psi_mhd,imin);}//  }
      if(cctk_bbox[2]) {                    YMIN_OB_LINEAR_EXTRAP(Avec_x,jmin); YMIN_OB_LINEAR_EXTRAP(Avec_y,jmin); YMIN_OB_LINEAR_EXTRAP(Avec_z,jmin); YMIN_OB_LINEAR_EXTRAP(psi_mhd,jmin);}//  }
      if((cctk_bbox[4]) && Symmetry_none) { ZMIN_OB_LINEAR_EXTRAP(Avec_x,kmin); ZMIN_OB_LINEAR_EXTRAP(Avec_y,kmin); ZMIN_OB_LINEAR_EXTRAP(Avec_z,kmin); ZMIN_OB_LINEAR_EXTRAP(psi_mhd,kmin);}//  }
    }
    else
    {
      if(cctk_bbox[1]) { XMAX_OB_LINEAR_EXTRAP(Avec_x,imax); XMAX_OB_LINEAR_EXTRAP(Avec_y,imax); XMAX_OB_LINEAR_EXTRAP(Avec_z,imax);}//  }
      if(cctk_bbox[3]) { YMAX_OB_LINEAR_EXTRAP(Avec_x,jmax); YMAX_OB_LINEAR_EXTRAP(Avec_y,jmax); YMAX_OB_LINEAR_EXTRAP(Avec_z,jmax);}//  }
      if(cctk_bbox[5]) { ZMAX_OB_LINEAR_EXTRAP(Avec_x,kmax); ZMAX_OB_LINEAR_EXTRAP(Avec_y,kmax); ZMAX_OB_LINEAR_EXTRAP(Avec_z,kmax);}//  }
      
      if(cctk_bbox[0]) {                    XMIN_OB_LINEAR_EXTRAP(Avec_x,imin); XMIN_OB_LINEAR_EXTRAP(Avec_y,imin); XMIN_OB_LINEAR_EXTRAP(Avec_z,imin);}//  }
      if(cctk_bbox[2]) {                    YMIN_OB_LINEAR_EXTRAP(Avec_x,jmin); YMIN_OB_LINEAR_EXTRAP(Avec_y,jmin); YMIN_OB_LINEAR_EXTRAP(Avec_z,jmin);}//  }
      if((cctk_bbox[4]) && Symmetry_none) { ZMIN_OB_LINEAR_EXTRAP(Avec_x,kmin); ZMIN_OB_LINEAR_EXTRAP(Avec_y,kmin); ZMIN_OB_LINEAR_EXTRAP(Avec_z,kmin);}//  }
    }
  }
}
