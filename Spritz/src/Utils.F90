 /*@@
   @file      Utils.F
   @date      Sat Jan 26 02:28:46 2002
   @author    
   @desc 
   Utility functions for other thorns. Calculation of the determinant
   of the spatial metric and the upper metric.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    Spritz_RefinementLevel
   @date       July 2005
   @author     
   @desc 
   Calculates the current refinement level from flesh data
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

subroutine Spritz_RefinementLevel(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS

  Spritz_reflevel = aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))

end subroutine Spritz_RefinementLevel


 /*@@
   @routine    SpatialDeterminant
   @date       Sat Jan 26 02:30:23 2002
   @author     
   @desc 
   Calculates the determinant of the spatial metric.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

subroutine SpatialDeterminant(gxx,gxy,gxz,gyy,gyz,gzz,det)

  implicit none
  
  CCTK_REAL :: det
  CCTK_REAL :: gxx,gxy,gxz,gyy,gyz,gzz

  det = -gxz**2*gyy + 2*gxy*gxz*gyz - gxx*gyz**2 - gxy**2*gzz + gxx*gyy*gzz

  
!!$  Why is this weird order necessary? Search me. It just seemed 
!!$  to make a really odd NaN go away.

!!$  det =  -gxz**2*gyy + 2*gxy*gxz*gyz 
!!$  det = det - gxx*gyz**2 - gxy**2*gzz
!!$  det = det + gxx*gyy*gzz
  
  return
  
end subroutine SpatialDeterminant



/*@@
@routine    UpperMetric
@date       Sat Jan 26 02:32:26 2002
@author     
@desc 
Calculates the upper metric. The determinant is given, not
calculated.
@enddesc 
@calls     
@calledby   
@history 

@endhistory 
@@*/

subroutine UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, &
     det, gxx, gxy, gxz, gyy, gyz, gzz)
  
  implicit none
  
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det, &
       gxx, gxy, gxz, gyy, gyz, gzz
  
  uxx=(-gyz**2 + gyy*gzz)/det
  uxy=(gxz*gyz - gxy*gzz)/det
  uyy=(-gxz**2 + gxx*gzz)/det
  uxz=(-gxz*gyy + gxy*gyz)/det
  uyz=(gxy*gxz - gxx*gyz)/det
  uzz=(-gxy**2 + gxx*gyy)/det
  
  return
  
end subroutine UpperMetric

