 /*@@
   @file      Spritz_Minima.F90
   @date      Mon Feb 25 11:43:36 2002
   @author    
   @desc 
   Sets up the scalars used for the atmosphere, before initial data.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    Spritz_Minima_Setup
   @date       Mon Feb 25 11:25:27 2002
   @author     Ian Hawke
   @desc 
   Before initial data, set up the scalar Spritz_rho_min used for the atmosphere.
   This is computed only from parameters.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Modified 30 Aug 2006 by Luca Baiotti
   @endhistory 

@@*/

subroutine Spritz_Rho_Minima_Setup(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (initial_rho_abs_min > 0.0) then
    Spritz_rho_min = initial_rho_abs_min
  else if (initial_rho_rel_min > 0.0) then 
    Spritz_rho_min = Spritz_rho_central * initial_rho_rel_min
  else if (rho_abs_min > 0.d0) then
    Spritz_rho_min = rho_abs_min
  else
    Spritz_rho_min = Spritz_rho_central * rho_rel_min
  end if

  if (initial_atmosphere_factor > 0.0) Spritz_rho_min = Spritz_rho_min * initial_atmosphere_factor

  return

end subroutine Spritz_Rho_Minima_Setup


 /*@@
   @routine    Spritz_Check_Rho_Minimum
   @date       Mon Jul 7 16:35:45 2008
   @author     Luca Baiotti 
   @desc 
   Check whether at some point rho < Spritz_rho_min and print a warning in case.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @endhistory 
@@*/

subroutine Spritz_Check_Rho_Minimum(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT i,j,k
  character(len=100) warnline

  do i=1,cctk_lsh(1)
    do j=1,cctk_lsh(2)
      do k=1,cctk_lsh(3)

        if (rho(i,j,k) < Spritz_rho_min) then
          call CCTK_WARN(2,"rho<Spritz_rho_min!!!")
          write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
          call CCTK_WARN(2,warnline)
          write(warnline,'(a25,g15.6)') 'Spritz_rho_min: ',Spritz_rho_min
          call CCTK_WARN(2,warnline)
          write(warnline,'(a25,g15.6)') 'rho: ',rho(i,j,k)
          call CCTK_WARN(2,warnline)
          write(warnline,'(a25,4g15.6)') 'coordinates: x,y,z,r:',x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
          call CCTK_WARN(2,warnline)
        end if

      end do
    end do
  end do

  return

end subroutine Spritz_Check_Rho_Minimum


 /*@@
   @routine    Spritz__Change_Rho_Minimum_At_Recovery
   @date       Thu Aug 14 17:11:32 2008
   @author     Luca Baiotti 
   @desc 
   Change, via a parameter, the value of Spritz_rho_min at recovery.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @endhistory 
@@*/

subroutine Spritz_Change_Rho_Minimum_At_Recovery(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  Spritz_rho_min = rho_abs_min_after_recovery

  return

end subroutine Spritz_Change_Rho_Minimum_At_Recovery
