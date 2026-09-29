 /*@@
   @file      Spritz_Weights.F90
   @date      Sat Jan 26 01:40:14 2002
   @author    Ian Hawke
   @desc 
   This routine calculates the "weights" of the cells and cell boundaries.
   This is only required if using FishEye.

   There is also a routine that calculates the "physical" velocity.

   @enddesc 
 @@*/
   
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "FishEye.h"

 /*@@
   @routine    Spritz_Weights
   @date       Sat Jan 26 01:41:02 2002
   @author     Ian Hawke
   @desc 
   Calculates the weights from FishEye.
   @enddesc 
   @calls     
   @calledby   
   @history 

   @endhistory 

@@*/

void Spritz_Weights(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL  dx, dy, dz, xbound, ybound, zbound;
  CCTK_REAL  jdet[27], jacobian[9], cell_centre_slice, cell_surface_l;
  CCTK_INT  i, j, k, nx, ny, nz, l, m, n, index, totalsize;

  dx = cctk_delta_space[0] / cctk_levfac[0];
  dy = cctk_delta_space[1] / cctk_levfac[1];
  dz = cctk_delta_space[2] / cctk_levfac[2];
  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];
  totalsize = nx*ny*nz;

  if (use_weighted_fluxes)
  {
    
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {

          /*
            Calculate the jacobian of the transformation at the 
            27 points; 3 points in every direction
          */

          index = CCTK_GFINDEX3D(cctkGH, i, j, k);

          for (n = 0; n < 3; n++)
          {
            for (m = 0; m < 3; m++)
            {
              for (l = 0; l < 3; l++)
              {
                
                xbound = x[index] + 0.5 * (l-2) * dx;
                ybound = y[index] + 0.5 * (m-2) * dy;
                zbound = z[index] + 0.5 * (n-2) * dz;

#ifdef FISHEYE_ACTIVE
                activejacobian(xbound,ybound,zbound,
                               &jacobian[0],&jacobian[1],&jacobian[2],
                               &jacobian[3],&jacobian[4],&jacobian[5],
                               &jacobian[6],&jacobian[7],&jacobian[8]);
#else
                CCTK_WARN(0,"You must compile with FishEye for this to work!");
#endif
                jdet[l+3*m+9*n] = 
                  -jacobian[2]*jacobian[4]*jacobian[6] + 
                  jacobian[1]*jacobian[5]*jacobian[6] + 
                  jacobian[2]*jacobian[3]*jacobian[7] - 
                  jacobian[0]*jacobian[5]*jacobian[7] - 
                  jacobian[1]*jacobian[3]*jacobian[8] + 
                  jacobian[0]*jacobian[4]*jacobian[8];

              }
            }
          }

          cell_jdet[index] = jdet[0];
          
          /* 
             Integrate by Simpsons rule over both directions
          */

          /*
            x
          */

          cell_surface[index] = 
            1.0 / 36.0 * ( (jdet[2] + 4.0 * jdet[5] + jdet[8])+ 
                          4.0 *(jdet[11] + 4.0 * jdet[14] + jdet[17])+
                          (jdet[20] + 4.0 * jdet[23] + jdet[26]) );

          /*
            y
          */

          cell_surface[index + totalsize] = 
            1.0 / 36.0 * ( (jdet[6] + 4.0 * jdet[7] + jdet[8])+
                          4.0 *(jdet[15] + 4.0 * jdet[16] + jdet[17])+
                          (jdet[24] + 4.0 * jdet[25] + jdet[26]) );

          /*
            z
          */

          cell_surface[index + 2 * totalsize] = 
            1.0 / 36.0 * ( (jdet[18] + 4.0 * jdet[19] + jdet[20])+
                          4.0 *(jdet[21] + 4.0 * jdet[22] + jdet[23])+
                          (jdet[24] + 4.0 * jdet[25] + jdet[26]) );
          
          /* 
             scalar for the x direction at the left, needed for
             the volume calculations
          */

          cell_surface_l = 
            1.0 / 36.0 * ( (jdet[0] + 4.0 * jdet[3] + jdet[6])+ 
                          4.0 *(jdet[9] + 4.0 * jdet[12] + jdet[15])+
                          (jdet[18] + 4.0 * jdet[21] + jdet[24]));

          /*
            The slice through the centre holding x fixed.
          */

          cell_centre_slice = 
            1.0 / 36.0 * ( (jdet[1] + 4.0 * jdet[4] + jdet[7])+ 
                          4.0 *(jdet[10] + 4.0 * jdet[13] + jdet[16])+
                          (jdet[19] + 4.0 * jdet[22] + jdet[25]) );

          /*
            The cell volume.
            Again using Simpsons rule, only in the x direction.
          */

          cell_volume[index] = 
            1.0 / 6.0 * ( cell_surface_l + 4.0 * cell_centre_slice + 
                          cell_surface[index] );
          
          
        }
      }
    }
    
  }
  else
  {
    
    for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {

          index = CCTK_GFINDEX3D(cctkGH, i, j, k);
    
          cell_volume[index] = 1.0;
          cell_surface[index] = 1.0;
          cell_surface[index + totalsize] = 1.0;
          cell_surface[index + 2 * totalsize] = 1.0;
          
        }
      }
    }
    
  }
  
  return;
}



 /*@@
   @routine    Spritz_Fisheye_Analysis
   @date       
   @author     Christian Ott, Ian Hawke
   @desc 
   Prepare output of Spritz variables in physical coordinates
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int Spritz_Fisheye_Analysis(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  double jacobian[9];

  int i,j,k, index;  

  /* let's first do the scalars: simple copying */

  for(i=0;
	i<cctkGH->cctk_lsh[0]*cctkGH->cctk_lsh[1]*cctkGH->cctk_lsh[2];
       i++) 
  {

    frho[i] = rho[i];
    fpress[i] = press[i];
    feps[i] = eps[i];

  }    

  /* now the vectors: more complicated... first we need to (numerically) 
     calculate the jacobian */


  for(i=0;i<cctkGH->cctk_lsh[0];i++)
  {
    for(j=0;j<cctkGH->cctk_lsh[1];j++)
    {
      for(k=0;k<cctkGH->cctk_lsh[2];k++)
      {
       
        index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        
#ifdef FISHEYE_ACTIVE
        activejacobian(x[index], y[index], z[index],
                       &jacobian[0],&jacobian[1],&jacobian[2],
                       &jacobian[3],&jacobian[4],&jacobian[5],
                       &jacobian[6],&jacobian[7],&jacobian[8]);
#else
        CCTK_WARN(1,"The output of \"FishEye\" velocities is only meaningful if FishEye is active!");
        jacobian[0]=1.0;
        jacobian[1]=0.0;
        jacobian[2]=0.0;
        jacobian[3]=0.0;
        jacobian[4]=1.0;
        jacobian[5]=0.0;
        jacobian[6]=0.0;
        jacobian[7]=0.0;
        jacobian[8]=1.0;
#endif

        fvelx[index] =
          jacobian[0] * velx[index] +
          jacobian[1] * vely[index] +
          jacobian[2] * velz[index];
        
        fvely[index] =
          jacobian[3] * velx[index] +
          jacobian[4] * vely[index] +
          jacobian[5] * velz[index];
        
        fvelz[index] =
          jacobian[6] * velx[index] +
          jacobian[7] * vely[index] +
          jacobian[8] * velz[index];
        
      }
    }
  }
  
  return;

}
