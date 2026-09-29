 /*@@
   @file      Spritz_Maxima.F90
   @date      Mon May 17 07:05:55 2004
   @author    Ian Hawke, Bruno Giacomazzo
   @desc 
   Find the maximum of the density over the whole grid.
   @enddesc 
 @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Spritz_FindMaxima
   @date       Mon May 17 07:06:49 2004
   @author     Ian Hawke, Bruno Giacomazzo
   @desc 
   Find the maximum of the density over the whole grid.
   Some of this code is borrowed, with alterations, from
   Peter Dieners EHFinder.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_InitFindMaxima(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  Spritz_maxima_iteration = -1

  maximum_density_plus  = 0.0d0 !maximum of rho in x>0 (this is actually the max of sign(x)*rho)
  maximum_density_minus = 0.0d0 !maximum of rho in x<0 (this is actually the min of sign(x)*rho)

  return

end subroutine Spritz_InitFindMaxima


subroutine Spritz_FindMaxima(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CALL Spritz_FindMaxima_Guts(                                      &
       & cctkGH, cctk_lsh, cctk_iteration, Spritz_maxima_iteration, &
       & maximum_density_plus, maximum_density_minus,               &
       & maxima_x_plus,  maxima_y_plus,  maxima_z_plus,             &
       & maxima_i_plus,  maxima_j_plus,  maxima_k_plus,             &
       & maxima_x_minus, maxima_y_minus, maxima_z_minus,            &
       & maxima_i_minus, maxima_j_minus, maxima_k_minus,            &
       & rho, x, y, z )

  !note maxrho_global is a protected variable which is needed by
  !Spritz_MoveCarpetGrids for the parameter
  !stop_moving_grids_when_rho_max_has_grown_a_factor
  !
  !note that maxrho_global_minus is negative (since it is the minimum of sign(x)*rho)
  maxrho_global_plus = maximum_density_plus
  maxrho_global_minus = maximum_density_minus
  maxrho_global = max(maxrho_global_plus,abs(maxrho_global_minus))

  return

end subroutine Spritz_FindMaxima


!except for not being called from the scheduler and for the declaration of variables, this routine is identical to Spritz_FindMaxima in Spritz_Dev.
subroutine Spritz_FindMaxima_Guts(                                &
     & cctkGH, cctk_lsh, cctk_iteration, Spritz_maxima_iteration, &
     & maximum_density_plus, maximum_density_minus,               &
     & maxima_x_plus,  maxima_y_plus,  maxima_z_plus,             &
     & maxima_i_plus,  maxima_j_plus,  maxima_k_plus,             &
     & maxima_x_minus, maxima_y_minus, maxima_z_minus,            &
     & maxima_i_minus, maxima_j_minus, maxima_k_minus,            &
     & rho, x, y, z )

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS


  CCTK_POINTER, intent(IN) :: cctkGH
  CCTK_INT,  intent(IN), dimension(3) :: cctk_lsh
  CCTK_INT,  intent(IN)    :: cctk_iteration
  CCTK_INT,  intent(INOUT) :: Spritz_maxima_iteration
  CCTK_REAL, intent(INOUT) :: maximum_density_plus, maximum_density_minus
  CCTK_INT,  intent(OUT)   :: maxima_i_plus,  maxima_j_plus,  maxima_k_plus
  CCTK_REAL, intent(OUT)   :: maxima_x_plus,  maxima_y_plus,  maxima_z_plus
  CCTK_INT,  intent(OUT)   :: maxima_i_minus, maxima_j_minus, maxima_k_minus
  CCTK_REAL, intent(OUT)   :: maxima_x_minus, maxima_y_minus, maxima_z_minus
  CCTK_REAL, intent(IN), dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: rho, x, y, z

! local variables
  CCTK_INT :: imin, imax, jmin, jmax, kmin, kmax 
  CCTK_INT :: ierr, max_handle, sum_handle, min_handle
  CCTK_INT, dimension(3) :: maxrho_local_location_plus, maxrho_location_plus
  CCTK_REAL, dimension(3) :: maxrho_local_position_plus, maxrho_position_plus
  CCTK_REAL :: maxrho_local_plus, maxrho_global_plus
  CCTK_INT, dimension(3) :: maxrho_local_location_minus, maxrho_location_minus
  CCTK_REAL, dimension(3) :: maxrho_local_position_minus, maxrho_position_minus
  CCTK_REAL :: maxrho_local_minus, maxrho_global_minus

  maxrho_local_location_plus = 0 
  maxrho_location_plus       = 0 
  maxrho_local_position_plus = 0.d0
  maxrho_position_plus       = 0.d0

  maxrho_local_location_minus = 0 
  maxrho_location_minus       = 0 
  maxrho_local_position_minus = 0.d0
  maxrho_position_minus       = 0.d0

  if (mod(cctk_iteration, DCandMaxima_every) .gt. 0) then
     return
  end if

  if (cctk_iteration > Spritz_maxima_iteration) then
     Spritz_maxima_iteration = cctk_iteration
     maximum_density_plus = 0.d0
     maximum_density_minus = 0.d0
  end if
  ! remember that in fortran array indeces start from 1
  imin = Spritz_stencil + 1
  jmin = Spritz_stencil + 1
  kmin = Spritz_stencil + 1
  imax = cctk_lsh(1) - Spritz_stencil
  jmax = cctk_lsh(2) - Spritz_stencil
  kmax = cctk_lsh(3) - Spritz_stencil

  call CCTK_ReductionArrayHandle (max_handle, 'maximum')
  call CCTK_ReductionArrayHandle (min_handle, 'minimum')
  call CCTK_ReductionArrayHandle (sum_handle, 'sum')

! The below is the new API for reductions but Erik says it is still buggy when used with carpet.
!  call CCTK_LocalArrayReductionHandle (sum_handle, 'sum')


  !find the maximum of rho on x>0 and x<0
  !note that maxrho_local_minus is a negative number
  maxrho_local_location_plus  = maxloc(sign(rho(imin:imax,jmin:jmax,kmin:kmax),x(imin:imax,jmin:jmax,kmin:kmax)))


!!$  write(*,*) "maxrho_local_location_plus", maxrho_local_location_plus
!!$  write(*,*) "maxrho_local_location_minus", maxrho_local_location_minus
!!$  
!!$  write(*,*) "xplus",  x (maxrho_local_location_plus(1) + Spritz_stencil, maxrho_local_location_plus(2) + Spritz_stencil, maxrho_local_location_plus(3) + Spritz_stencil)
!!$  write(*,*) "xminus",  x (maxrho_local_location_minus(1) + Spritz_stencil, maxrho_local_location_minus(2) + Spritz_stencil, maxrho_local_location_minus(3) + Spritz_stencil)

  maxrho_local_plus = sign(rho (&
       maxrho_local_location_plus(1) + Spritz_stencil, &
       maxrho_local_location_plus(2) + Spritz_stencil, &
       maxrho_local_location_plus(3) + Spritz_stencil  ), x (&
       maxrho_local_location_plus(1) + Spritz_stencil, &
       maxrho_local_location_plus(2) + Spritz_stencil, &
       maxrho_local_location_plus(3) + Spritz_stencil  ))


!!$  Only check the x>0 part of the grid. This ensures that if there are 
!!$  two stars either side of x=0 only one will be picked every time.
!!$
!!$  if ( Spritz_useFullGridForMaxima .eq. 0 ) then
!!$    if ( x (&
!!$         maxrho_local_location(1) + Spritz_stencil, &
!!$         maxrho_local_location(2) + Spritz_stencil, &
!!$         maxrho_local_location(3) + Spritz_stencil  ) < 0.d0 ) then
!!$      maxrho_local = 0.d0 ! this is how maxrho_local can be seen to be zero in some output
!!$    end if
!!$  end if
  
  ierr = 0 !initialize ierr so that Valgrind doesn't complain


  !----------------------------------------------------------!
  !--- FIND the location of the maximum in the region x>0 ---!
  !----------------------------------------------------------!


  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
       maxrho_local_plus, maxrho_global_plus, CCTK_VARIABLE_REAL)
  if ( ierr .ne. 0 ) then
     call CCTK_WARN(0, 'Reduction of maxrho failed')
  end if
  
  if (maxrho_local_plus .eq. maxrho_global_plus) then
!!$     write(*,*) "New Test Spritz_Maxima FOUND MAXIMA LOCATION PLUS"
!!$     write(*,*) "maxrho_local_location_plus", maxrho_local_location_plus
!!$     write(*,*) "maxrho_local_location_minus", maxrho_local_location_minus 
     maxrho_local_position_plus(1) = x (&
          maxrho_local_location_plus(1) + Spritz_stencil, &
          maxrho_local_location_plus(2) + Spritz_stencil, &
          maxrho_local_location_plus(3) + Spritz_stencil )
     maxrho_local_position_plus(2) = y (&
          maxrho_local_location_plus(1) + Spritz_stencil, &
          maxrho_local_location_plus(2) + Spritz_stencil, &
          maxrho_local_location_plus(3) + Spritz_stencil )
     maxrho_local_position_plus(3) = z (&
          maxrho_local_location_plus(1) + Spritz_stencil, &
          maxrho_local_location_plus(2) + Spritz_stencil, &
          maxrho_local_location_plus(3) + Spritz_stencil )
  else
     
     !if the processor does not have the global maximum, set maxrho_local_position to zero 
     !(for use in the below "sum" reduction)
     maxrho_local_position_plus = 0.d0 !it is a vector
     
     
     !if the processor does not have the global maximum, set maxrho_local_location to zero 
     !(for use in the below "sum" reduction)
     maxrho_local_location_plus = 0 !it is a vector
     
  end if

!!$  write(*,*)
!!$  write(*,*) "2) On reflevel",aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))
!!$  write(*,*) "Max rho on level",maxrho_local
!!$  write(*,*) "Max rho this time step up to now",maximum_density
!!$  write(*,*) "Mrho_loc_loc",maxrho_local_location
!!$  write(*,*) "Mrho_loc_pos",maxrho_local_position

  ierr = 0 !initialize ierr so that Valgrind doesn't complain
  
  !find maxrho_local_position on the processor that has the maximum on x>0
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
       maxrho_local_position_plus, maxrho_position_plus, 3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
     call CCTK_WARN(0, 'Reduction of maxrho_position failed')
  end if

  !find maxrho_local_location on the processor that has the maximum on x>0
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
       maxrho_local_location_plus, maxrho_location_plus, 3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
     call CCTK_WARN(0, 'Reduction of maxrho_location failed')
  end if

!!$  write(*,*) "Mrho_glob_loc",maxrho_location
!!$  write(*,*) "Mrho_glob_pos",maxrho_position


!!$  If the maximum on this iteration was larger on another level, 
!!$  do not set the position

  if (maxrho_global_plus > maximum_density_plus) then
!!$  Only set values if well above atmosphere level
     if (maxrho_global_plus .ge. rho_max_min ) then
        
        maximum_density_plus = maxrho_global_plus

        maxima_x_plus = maxrho_position_plus(1)
        maxima_y_plus = maxrho_position_plus(2)
        maxima_z_plus = maxrho_position_plus(3)
        
        maxima_i_plus = maxrho_location_plus(1)
        maxima_j_plus = maxrho_location_plus(2)
        maxima_k_plus = maxrho_location_plus(3)
        
     end if
  end if

!!$  write(*,*)
!!$  write(*,*) "CHECK MAXIMA LOCATION MINUS (1)"
!!$  write(*,*) "maxrho_local_location_plus", maxrho_local_location_plus
!!$  write(*,*) "maxrho_local_location_minus", maxrho_local_location_minus 


  !----------------------------------------------------------!
  !--- FIND the location of the maximum in the region x<0 ---!
  !----------------------------------------------------------!

  if (Spritz_useFullGridForMaxima /= 0) then

  maxrho_local_location_minus = minloc(sign(rho(imin:imax,jmin:jmax,kmin:kmax),x(imin:imax,jmin:jmax,kmin:kmax)))


  maxrho_local_minus = sign(rho (&
       maxrho_local_location_minus(1) + Spritz_stencil, &
       maxrho_local_location_minus(2) + Spritz_stencil, &
       maxrho_local_location_minus(3) + Spritz_stencil  ), x (&
       maxrho_local_location_minus(1) + Spritz_stencil, &
       maxrho_local_location_minus(2) + Spritz_stencil, &
       maxrho_local_location_minus(3) + Spritz_stencil  ))


!!$  write(*,*)
!!$  write(*,*) "CHECK MAXIMA LOCATION MINUS (2)"
!!$  write(*,*) "maxrho_local_location_plus", maxrho_local_location_plus
!!$  write(*,*) "maxrho_local_location_minus", maxrho_local_location_minus 
!!$  write(*,*)


  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, min_handle, &
       maxrho_local_minus, maxrho_global_minus, CCTK_VARIABLE_REAL)
  if ( ierr .ne. 0 ) then
     call CCTK_WARN(0, 'Reduction of maxrho failed')
  end if
  
!!$  write(*,*) "New Test Spritz_Maxima"
!!$  write(*,*) "maxrho_local_minus", maxrho_local_minus
!!$  write(*,*) "maxrho_global_minus", maxrho_global_minus
!!$  write(*,*) "maxrho_local_location_plus", maxrho_local_location_plus
!!$  write(*,*) "maxrho_local_location_minus", maxrho_local_location_minus 
!!$  write(*,*)

  if (abs(maxrho_local_minus-maxrho_global_minus).le.1.0d-14) then
!!$     write(*,*) "New Test Spritz_Maxima FOUND MAXIMA LOCATION MINUS"
!!$     write(*,*) "maxrho_local_location_plus", maxrho_local_location_plus
!!$     write(*,*) "maxrho_local_location_minus", maxrho_local_location_minus 
     maxrho_local_position_minus(1) = x (&
          maxrho_local_location_minus(1) + Spritz_stencil, &
          maxrho_local_location_minus(2) + Spritz_stencil, &
          maxrho_local_location_minus(3) + Spritz_stencil )
     maxrho_local_position_minus(2) = y (&
          maxrho_local_location_minus(1) + Spritz_stencil, &
          maxrho_local_location_minus(2) + Spritz_stencil, &
          maxrho_local_location_minus(3) + Spritz_stencil )
     maxrho_local_position_minus(3) = z (&
          maxrho_local_location_minus(1) + Spritz_stencil, &
          maxrho_local_location_minus(2) + Spritz_stencil, &
          maxrho_local_location_minus(3) + Spritz_stencil )

!!$     write(*,*) "maxrho_local_position_minus(1)", maxrho_local_position_minus(1)
!!$     write(*,*) "maxrho_local_position_minus(2)", maxrho_local_position_minus(2)
!!$     write(*,*) "maxrho_local_position_minus(3)", maxrho_local_position_minus(3)
!!$     write(*,*)

  else
     
     !if the processor does not have the global maximum, set maxrho_local_position to zero 
     !(for use in the below "sum" reduction)
     maxrho_local_position_minus = 0.d0 !it is a vector

     
     !if the processor does not have the global maximum, set maxrho_local_location to zero 
     !(for use in the below "sum" reduction)
     maxrho_local_location_minus= 0 !it is a vector
     
  end if

  ierr = 0 !initialize ierr so that Valgrind doesn't complain

  !find maxrho_local_position_minus on the processor that has the maximum in x<0
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
       maxrho_local_position_minus, maxrho_position_minus, 3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
     call CCTK_WARN(0, 'Reduction of maxrho_position failed')
  end if

!!$  write(*,*) "CHECK waht happens after reduction..."
!!$  write(*,*) "maxrho_local_position_minus", maxrho_local_position_minus
!!$  write(*,*) "maxrho_position_minus", maxrho_position_minus



  !find maxrho_local_location_minus on the processor that has the maximum in x<0
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
       maxrho_local_location_minus, maxrho_location_minus, 3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
     call CCTK_WARN(0, 'Reduction of maxrho_location failed')
  end if

!!$  If the maximum on this iteration was larger on another level, 
!!$  do not set the position

  if (maxrho_global_minus < maximum_density_minus) then
!!$  Only set values if well above atmosphere level
     if (abs(maxrho_global_minus) .ge. rho_max_min ) then
        
        maximum_density_minus = maxrho_global_minus
        
        maxima_x_minus = maxrho_position_minus(1)
        maxima_y_minus = maxrho_position_minus(2)
        maxima_z_minus = maxrho_position_minus(3)
        
        maxima_i_minus = maxrho_location_minus(1)
        maxima_j_minus = maxrho_location_minus(2)
        maxima_k_minus = maxrho_location_minus(3)
        
     end if
  end if


!!$  write(*,*) "maxima location",maxima_i,maxima_j,maxima_k

  ! write(*,*) "maxima position plus",maxima_x_plus, maxima_y_plus, maxima_z_plus
  ! write(*,*) "maxima position minus",maxima_x_minus, maxima_y_minus, maxima_z_minus
  ! write(*,*) "maxrho_global_plus, maximum_density_plus", maxrho_global_plus, maximum_density_plus
  ! write(*,*) "maxrho_global_minus, maximum_density_minus", maxrho_global_minus, maximum_density_minus


else

        maximum_density_minus = 0.0
        
        maxima_x_minus = 0.0
        maxima_y_minus = 0.0
        maxima_z_minus = 0.0
        
        maxima_i_minus = 0
        maxima_j_minus = 0
        maxima_k_minus = 0

end if

  return

end subroutine Spritz_FindMaxima_Guts


! /*@@
!   @routine    Spritz_FindWeightedMaxima
!   @date       Mon May 24 23:22:18 2004
!   @author     Ian Hawke
!   @desc 
!   Find the volume weighted centre of mass, 
!            /
!           |       i
!           |  rho x    dV
!           |
!  - i      /
!  x  =     ----------------
!           /             
!           |             
!           |  rho      dV 
!           |              
!           /              
!
!   @enddesc 
!   @calls     
!   @calledby   
!   @history 
! 
!   @endhistory 
!
!@@*/
!
!subroutine Spritz_FindWeightedMaxima(CCTK_ARGUMENTS)
!
!  implicit none
!
!  DECLARE_CCTK_ARGUMENTS
!  DECLARE_CCTK_PARAMETERS
!  DECLARE_CCTK_FUNCTIONS
!
!  CCTK_INT :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
!  CCTK_INT :: ierr, sum_handle
!  CCTK_REAL, dimension(3) :: weightedrho_local, weightedrho_global
!  CCTK_REAL :: rho_local, rho_global
!
!  imin = Spritz_stencil
!  jmin = Spritz_stencil
!  kmin = Spritz_stencil
!  imax = cctk_lsh(1) - Spritz_stencil
!  jmax = cctk_lsh(2) - Spritz_stencil
!  kmax = cctk_lsh(3) - Spritz_stencil
!
!  call CCTK_ReductionArrayHandle (sum_handle, 'sum')
!
!  weightedrho_local = 0.d0
!  rho_local = 0.d0
!
!  do k = kmin, kmax
!    do j = jmin, jmax
!      do i = imin, imax
!
!        weightedrho_local(1) = weightedrho_local(1) + &
!             rho(i, j, k) * x(i, j, k)
!        weightedrho_local(2) = weightedrho_local(2) + &
!             rho(i, j, k) * y(i, j, k)
!        weightedrho_local(3) = weightedrho_local(3) + &
!!!$             rho(i, j, k) * z(i, j, k)
!!!$  To get around problems with symmetry boundaries and because
!!!$  the BNS do not have spin (yet) I just set \bar{z} to zero.             
!             0.d0
!        rho_local = rho_local + rho(i,j,k)
!
!      end do
!    end do
!  end do
!  
!  ierr = 0 !initialize ierr so that Valgrind doesn't complain
!
!  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
!       weightedrho_local, weightedrho_global, 3, CCTK_VARIABLE_REAL )
!  if ( ierr .ne. 0 ) then
!    call CCTK_WARN(0, 'Reduction of weightedrho failed')
!  end if
!  
!  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, sum_handle, &
!       rho_local, rho_global, CCTK_VARIABLE_REAL )
!  if ( ierr .ne. 0 ) then
!    call CCTK_WARN(0, 'Reduction of rho failed')
!  end if
!
!  maxima_x_plus = weightedrho_global(1) / rho_global
!  maxima_y_plus = weightedrho_global(2) / rho_global
!  maxima_z_plus = weightedrho_global(3) / rho_global
!
!  return  
!
!end subroutine Spritz_FindWeightedMaxima
!
!
!!$ /*@@
!!$   @routine    Spritz_SetDCCentroid
!!$   @date       Mon May 17 07:06:49 2004
!!$   @author     Ian Hawke
!!$   @desc 
!!$
!!$   @enddesc 
!!$   @calls     
!!$   @calledby   
!!$   @history 
!!$ 
!!$   @endhistory 
!!$
!!$@@*/
!!$
!!$subroutine Spritz_SetDCCentroid(CCTK_ARGUMENTS)
!!$
!!$  implicit none
!!$
!!$  DECLARE_CCTK_ARGUMENTS
!!$  DECLARE_CCTK_PARAMETERS
!!$  DECLARE_CCTK_FUNCTIONS
!!$
!!$  CCTK_INT :: ierr
!!$
!!$  ierr = 0 !initialize ierr so that Valgrind doesn't complain
!!$
!!$  if (set_dc_centroid .ne. 0) then
!!$
!!$    call CCTK_IsFunctionAliased(ierr, "SetDriftCorrectPosition")
!!$    if ( ierr .ne. 0 ) then
!!$      call SetDriftCorrectPosition ( cctkGH, &
!!$           maxima_x_plus, maxima_y_plus, maxima_z_plus )
!!$    end if
!!$  
!!$  end if
!!$  
!!$  return
!!$
!!$end subroutine Spritz_SetDCCentroid


 /*@@
   @routine    Spritz_FindSeparation
   @date       Thu May 20 12:35:20 2004
   @author     Ian Hawke
   @desc 
   Finds the separation (in coordinate and proper distance)
   between the NS. Equal mass symmetry is assumed so it is
   actually the distance between the origin and the location
   of maximum density. This is along a straight line, not a
   geodesic, so still not perfect.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_FindSeparation(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: separation_dx, separation_dy, separation_dz
  CCTK_REAL, dimension(:), allocatable :: separation_x, separation_y, &
       separation_z
  CCTK_REAL, dimension(:), allocatable :: s_gxx, s_gxy, s_gxz, &
       s_gyy, s_gyz, s_gzz
  CCTK_INT :: ierr, i
  CCTK_INT :: param_table_handle, interp_handle, coord_system_handle
  CCTK_INT :: vindex

  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_INT, dimension(6) :: in_array_indices
  CCTK_POINTER, dimension(6) :: out_arrays
  CCTK_INT, dimension(6) :: out_array_type_codes

  ierr = 0 !initialize ierr so that Valgrind doesn't complain

!!$  Coordinate separation is easy.

  IF ( Spritz_useFullGridForMaxima == 0 ) THEN
     Spritz_separation =                                    &
          & SQRT( maxima_x_plus**2 + maxima_y_plus**2 + maxima_z_plus**2 )
     maxima_x_minus = -maxima_x_plus
     maxima_y_minus = -maxima_y_plus
     maxima_z_minus =  maxima_z_plus
  ELSE
     Spritz_separation =                                    &
          & SQRT(   ( maxima_x_plus - maxima_x_minus )**2   &
          &       + ( maxima_y_plus - maxima_y_minus )**2   &
          &       + ( maxima_z_plus - maxima_z_minus )**2 )
  END IF

!!$  Proper separation requires interpolation

  allocate(separation_x(separation_npoints), &
       separation_y(separation_npoints), &
       separation_z(separation_npoints), STAT=ierr)

  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Failed to allocate separation coordinate arrays")
  end if
  
  separation_dx = maxima_x_plus / dble(separation_npoints)
  separation_dy = maxima_y_plus / dble(separation_npoints)
  separation_dz = maxima_z_plus / dble(separation_npoints)

  do i = 1, separation_npoints
    
    separation_x(i) = i * separation_dx
    separation_y(i) = i * separation_dy
    separation_z(i) = i * separation_dz

  end do

  ierr = 0 !initialize ierr so that Valgrind doesn't complain
  
  allocate(s_gxx(separation_npoints), &
       s_gxy(separation_npoints), &
       s_gxz(separation_npoints), &
       s_gyy(separation_npoints), &
       s_gyz(separation_npoints), &
       s_gzz(separation_npoints), &
       STAT=ierr)

  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Failed to allocate separation metric arrays")
  end if

  param_table_handle = -1
  interp_handle = -1
  coord_system_handle = -1

  call Util_TableCreateFromString (param_table_handle, "order = 2")
  if (param_table_handle .lt. 0) then
    call CCTK_WARN(0,"Cannot create parameter table for interpolator")
  endif

!  call CCTK_InterpHandle (interp_handle, "uniform cartesian")
  call CCTK_InterpHandle (interp_handle, "Lagrange polynomial interpolation")
  if (interp_handle.lt.0) then
    call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
  endif

  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
    call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
  endif

!     fill in the input/output arrays for the interpolator
  interp_coords(1) = CCTK_PointerTo(separation_x)
  interp_coords(2) = CCTK_PointerTo(separation_y)
  interp_coords(3) = CCTK_PointerTo(separation_z)
  
  call CCTK_VarIndex (vindex, "admbase::gxx")
  in_array_indices(1) = vindex
  call CCTK_VarIndex (vindex, "admbase::gyy")
  in_array_indices(2) = vindex
  call CCTK_VarIndex (vindex, "admbase::gzz")
  in_array_indices(3) = vindex
  call CCTK_VarIndex (vindex, "admbase::gxy")
  in_array_indices(4) = vindex
  call CCTK_VarIndex (vindex, "admbase::gxz")
  in_array_indices(5) = vindex
  call CCTK_VarIndex (vindex, "admbase::gyz")
  in_array_indices(6) = vindex
  
  out_arrays(1) = CCTK_PointerTo(s_gxx)
  out_arrays(2) = CCTK_PointerTo(s_gyy)
  out_arrays(3) = CCTK_PointerTo(s_gzz)
  out_arrays(4) = CCTK_PointerTo(s_gxy)
  out_arrays(5) = CCTK_PointerTo(s_gxz)
  out_arrays(6) = CCTK_PointerTo(s_gyz)
  
  out_array_type_codes = CCTK_VARIABLE_REAL
  
!     Interpolation.

  ierr = 0 !initialize ierr so that Valgrind doesn't complain

  call CCTK_InterpGridArrays (ierr, cctkGH, 3, interp_handle,&
       param_table_handle, coord_system_handle,&
       separation_npoints, CCTK_VARIABLE_REAL, interp_coords,&
       6, in_array_indices,&
       6, out_array_type_codes, out_arrays)
  if (ierr < 0) then
    call CCTK_WARN (1, "interpolator call returned an error code");
  endif
  
!     release parameter table
  call Util_TableDestroy (ierr, param_table_handle)
  
!        Integrate using trapezoidal rule.
  
  Spritz_proper_separation = 0.d0

  do i=2, separation_npoints
    
    Spritz_proper_separation = Spritz_proper_separation + 0.5d0 &
         *(sqrt(s_gxx(i-1)*separation_dx**2  + &
                s_gyy(i-1)*separation_dy**2  + &
                s_gzz(i-1)*separation_dz**2 &
        + 2.d0*(s_gxy(i-1)*separation_dx*separation_dy + &
                s_gxz(i-1)*separation_dx*separation_dz + &
                s_gyz(i-1)*separation_dy*separation_dz)) &
         + sqrt(s_gxx(i  )*separation_dx**2  + &
                s_gyy(i  )*separation_dy**2  + &
                s_gzz(i  )*separation_dz**2 &
        + 2.d0*(s_gxy(i  )*separation_dx*separation_dy + &
                s_gxz(i  )*separation_dx*separation_dz + &
                s_gyz(i  )*separation_dy*separation_dz)))
    
  end do

  return

end subroutine Spritz_FindSeparation
