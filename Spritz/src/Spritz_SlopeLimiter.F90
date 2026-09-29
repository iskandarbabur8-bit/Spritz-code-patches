 /*@@
   @file      Spritz_SlopeLimiter.F90
   @date      Sat Jan 26 02:00:32 2002
   @author    
   @desc 
   The routine for the more complex slope limiters. 
   See Toros book for most of these.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

 /*@@
   @routine    slopelimiter
   @date       Sat Jan 26 02:01:15 2002
   @author     Luca Baiotti, Ian Hawke
   @desc 
   Given the slope delta and the ratio of the two local slopes,
   limit the slope so the reconstruction is TVD.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Original minmod from GR3D, author Mark Miller. Other limiters from 
   codes of Hawke, Nikiforakis and Toro.
   @endhistory 

@@*/

subroutine slopelimiter(r, delta)
  
  USE Spritz_Scalars
  
  implicit none
  
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_REAL :: phi, r, denor, phir, delta
  
  denor = 1.d0 + r
  phir = 2.d0/denor
  
  if (MC1) then
     
     if (r > 0.d0) then 
        phi = min( 2.d0*r/(1.d0+r), phir )
     else
        phi = 0.d0
     end if
    
  else if (MINMOD2) then
     
     if (r > 1.d0) then 
        phi = min(1.d0, phir )
     else if (r > 0.d0) then
        phi = r
     else
        phi = 0.d0
     end if
    
  else if (MINMOD3) then
     
     if (r > 0.d0) then 
        phi = min(1.d0, 4.d0/(1.d0+r) )
     else
        phi = 0.d0
     end if
    
  else if (SUPERBEE) then
     
     if (r > 1.d0) then
        phi = min(2.d0, r, 2.d0 / (1.d0 + r))
     else if (r > 0.5d0) then
        phi = 1.d0
     else if (r > 0.d0) then
        phi = 2.d0 * r
     else
        phi = 0.d0
     end if
     
  else 
    
     call CCTK_WARN(0, "Type of limiter not recognized")
     
  end if

  delta = delta * phi
  
end subroutine slopelimiter
  
