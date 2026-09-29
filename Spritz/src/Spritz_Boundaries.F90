 /*@@
   @file      Spritz_Boundaries.F90
   @date      Sat Jan 26 01:01:14 2002
   @author
   @desc
   The two routines for dealing with boundary conditions.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "util_Table.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)


 /*@@
   @routine    Spritz_InitSymBound
   @date       Sat Jan 26 01:03:04 2002
   @author     Ian Hawke
   @desc
   Sets up the symmetries at the boundaries of the hydrodynamical variables.
   @enddesc
   @calls
   @calledby
   @history
   Direct translation of routines from GR3D, GRAstro_Hydro,
   written by Mark Miller, or WaveToy routines, or...
   @endhistory

@@*/

subroutine Spritz_InitSymBound(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer, dimension(3) :: sym
  integer :: ierr

  sym(1) = 1
  sym(2) = 1
  sym(3) = 1


  call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::rho")
  call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::press")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::dens")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::tau")
  call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::w_lorentz")
  call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::eps")
  !call SetCartSymVN(ierr, cctkGH, sym, "Spritz::Spritz_C2P_failed")

  if(evolve_y_e.ne.0) then
     call SetCartSymGN(ierr, cctkGH, sym, "HydroBase::Y_e")
     call SetCartSymGN(ierr, cctkGH, sym, "Spritz::Y_e_con")
  endif

  if(GetTemp_FromThermo.ne.0) then
      call SetCartSymGN(ierr, cctkGH, sym, "HydroBase::temperature")
      call SetCartSymGN(ierr, cctkGH, sym, "HydroBase::entropy")
  endif


  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::psi_mhd")

  sym(1) = -1
  sym(2) = 1
  sym(3) = 1

  call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::vel[0]")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::scon[0]")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::Avec_x") !vector potential has the same symmetry of v^j


  sym(1) = -1
  sym(2) = 1
  sym(3) = -1
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::Bx")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::dBx")

  sym(1) = 1
  sym(2) = -1
  sym(3) = 1

  call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::vel[1]")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::scon[1]")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::Avec_y") !vector potential has the same symmetry of v^j

  sym(1) = 1
  sym(2) = -1
  sym(3) = -1
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::By")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::dBy")

  sym(1) = 1
  sym(2) = 1
  sym(3) = -1

  call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::vel[2]")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::scon[2]")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::Avec_z") !vector potential has the same symmetry of v^j

  sym(1) = 1
  sym(2) = 1
  sym(3) = 1
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::Bz")
  call SetCartSymVN(ierr, cctkGH, sym, "Spritz::dBz")

end subroutine Spritz_InitSymBound

 /*@@
   @routine    Spritz_Boundaries
   @date       Sat Jan 26 01:04:04 2002
   @author
   @desc
   Calls the appropriate boundary routines
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine Spritz_Boundaries(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: sw
  integer :: ierr = 0
!!$  CCTK_REAL :: pi = 3.141569d0

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1

  sw = Spritz_stencil

!!$Symmetry boundaries first

!!$  call CartSymGN(ierr,cctkGH,"Spritz::Spritz_cons_var")
!!$  call CartSymGN(ierr,cctkGH,"Spritz::Spritz_prim_var")


!!$Flat boundaries if required

  if (CCTK_EQUALS(bound,"flat")) then
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "Spritz::dens", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "Spritz::tau", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "Spritz::scon", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "HydroBase::w_lorentz", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "HydroBase::rho", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "HydroBase::press", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "HydroBase::eps", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "HydroBase::vel", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "Spritz::Spritz_mhd_vector_var", "Flat")
       ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
            "Spritz::Spritz_mhd_evolved_var", "Flat")
       if (CCTK_EQUALS(Spritz_Extr_OB_Avec_Psi,"no")) then
          !call CCTK_INFO("Setting Flat Outer Boundaries for Avec and Psi...")
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
               "Spritz::psi_mhd", "Flat")
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
               "Spritz::Avec_x", "Flat")
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
               "Spritz::Avec_y", "Flat")
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
               "Spritz::Avec_z", "Flat")
       end if

      if(evolve_y_e.ne.0) then
        if (sync_conserved_only .eq. 0) then
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
          "HydroBase::Y_e", "Flat")
        endif
        ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
        "Spritz::Y_e_con", "Flat")
      endif
      if(GetTemp_FromThermo.ne.0) then
        if (sync_conserved_only .eq. 0) then
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
          "HydroBase::temperature", "Flat")
        endif
        ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
        "HydroBase::entropy", "Flat")
      endif

  end if

  if (CCTK_EQUALS(bound,"none")) then
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "Spritz::dens", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "Spritz::tau", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "Spritz::scon", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "HydroBase::w_lorentz", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "HydroBase::rho", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "HydroBase::press", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "HydroBase::eps", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "HydroBase::vel", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "Spritz::Spritz_mhd_vector_var", "None")
      ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
           "Spritz::Spritz_mhd_evolved_var", "None")
      if (CCTK_EQUALS(Spritz_Extr_OB_Avec_Psi,"no")) then
         !call CCTK_INFO("Setting None Outer Boundaries for Avec and Psi...")
         ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
              "Spritz::psi_mhd", "None")
         ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
              "Spritz::Avec_x", "None")
         ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
              "Spritz::Avec_y", "None")
         ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
              "Spritz::Avec_z", "None")
      end if

      if(evolve_y_e.ne.0) then
        if (sync_conserved_only .eq. 0) then
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
          "HydroBase::Y_e", "None")
        endif
        ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
        "Spritz::Y_e_con", "None")
      endif
      if(GetTemp_FromThermo.ne.0) then
        if (sync_conserved_only .eq. 0) then
          ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
          "HydroBase::temperature", "None")
        endif
        ierr = Boundary_SelectGroupForBC(cctkGH, faces, Spritz_stencil, -ione, &
        "HydroBase::entropy", "None")
      endif

  end if

  if (CCTK_EQUALS(bound,"scalar")) then
     call CCTK_WARN(0, "Until somebody uses this I see no reason to support it")
  end if

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine Spritz_Boundaries
