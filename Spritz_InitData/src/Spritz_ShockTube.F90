 /*@@
   @file      Spritz_ShockTube.F90
   @date      
   @author    Bruno Giacomazzo
   @desc 
   Test with Riemann problems
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)

 /*@@
   @routine    Spritz_shocktube
   @date       
   @author     Bruno Giacomazzo
   @desc 
   Initial data for Riemann problems. Either diagonal or parallel to
   a coordinate axis.
   @enddesc 
@@*/

subroutine Spritz_shocktube(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: direction, det
  CCTK_REAL :: rhol, rhor, velxl, velxr, velyl, velyr, &
       velzl, velzr, epsl, epsr, Bxl, Bxr, Byl, Byr, Bzl, Bzr, &
       pressl, pressr
  CCTK_REAL :: dx, dy, dz

  CCTK_REAL,dimension(15)::value
  character(len=200)::dump,input_file_for
  character(len=500)::warnline
  CCTK_INT :: string_length
  CCTK_REAL :: small_val

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)
  CCTK_REAL :: xeos_gamma(1)
  
  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!!

  ! Hard-coded small values for checking the point position
  small_val = 1.0d-6

  if (CCTK_EQUALS(initial_shift,"none")) then
     call CCTK_WARN(0,"You have to set ADMBase::initial_shift=zero")
  end if

  psi_mhd = 0.0d0
  psi_mhd_p = 0.0d0
  if (timelevels == 3) then
     psi_mhd_p_p = 0.0d0
  end if

  if (CCTK_EQUALS(shock_case,"Simple")) then
     rhol = 10.d0
     rhor = 1.d0
     velxl = 0.d0
     velxr = 0.d0
     velyl = 0.d0
     velyr = 0.d0
     velzl = 0.d0
     velzr = 0.d0
     epsl = 2.d0
     epsr = 1.d-6
  else if (CCTK_EQUALS(shock_case,"Sod")) then
     rhol = 1.d0
     rhor = 0.125d0
     velxl = 0.d0
     velxr = 0.d0
     velyl = 0.d0
     velyr = 0.d0
     velzl = 0.d0
     velzr = 0.d0
     epsl = 1.5d0
     epsr = 1.2d0
!!$This line only for polytrope, k=1
!!$          epsr = 0.375d0
  else if (CCTK_EQUALS(shock_case,"Blast")) then
     rhol = 1.d0
     rhor = 1.d0
     velxl = 0.d0
     velxr = 0.d0
     velyl = 0.d0
     velyr = 0.d0
     velzl = 0.d0
     velzr = 0.d0
     epsl = 1500.d0
     epsr = 1.5d-2
  elseif (CCTK_EQUALS(shock_case,"KO_ST2")) then
     rhol  = 1.0d0
     rhor  = 0.1d0
     velxl = 0.0d0
     velxr = 0.0d0
     velyl = 0.0d0
     velyr = 0.0d0
     velzl = 0.0d0
     velzr = 0.0d0
     pressl= 30.0d0 
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl)
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 2.0d0 !1.0d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 0.0d0
     Bxr   = 0.0d0
     Byl   = 20.0d0
     Byr   = 0.0d0
     Bzl   = 0.0d0
     Bzr   = 0.0d0
  else if (CCTK_EQUALS(shock_case,"generic_ST")) then
     rhol  = 1.0d0
     rhor  = 0.01d0
     velxl = 0.1d0
     velxr = 0.5d0
     velyl = 0.3d0
     velyr = 0.4d0
     velzl = 0.4d0
     velzr = 0.3d0
     pressl= 0.01d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl)
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 5000.0d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 0.0d0
     Bxr   = 0.0d0
     Byl   = 6.0d0
     Byr   = 5.0d0
     Bzl   = 2.0d0
     Bzr   = 20.0d0
  else if (CCTK_EQUALS(shock_case,"KO_ST1")) then
     rhol  = 1.0d0
     rhor  = 0.1d0
     velxl = 0.0d0
     velxr = 0.0d0
     velyl = 0.0d0
     velyr = 0.0d0
     velzl = 0.0d0
     velzr = 0.0d0
     pressl= 1000.0d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl)
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 1.0d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 1.0d0
     Bxr   = 1.0d0
     Byl   = 0.0d0
     Byr   = 0.0d0
     Bzl   = 0.0d0
     Bzr   = 0.0d0
  else if (CCTK_EQUALS(shock_case,"KO_Collision")) then
     rhol  = 1.0d0
     rhor  = 1.0d0
     velxl = 5.0d0/sqrt(26.0d0)
     velxr = -5.0d0/sqrt(26.0d0)
     velyl = 0.0d0
     velyr = 0.0d0
     velzl = 0.0d0
     velzr = 0.0d0
     pressl= 1.0d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl)
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 1.0d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 10.0d0
     Bxr   = 10.0d0
     Byl   = 10.0d0
     Byr   = -10.0d0
     Bzl   = 0.0d0
     Bzr   = 0.0d0
  else if (CCTK_EQUALS(shock_case,"Balsara1")) then
     rhol  = 1.0d0
     rhor  = 0.125d0
     velxl = 0.0d0
     velxr = 0.0d0
     velyl = 0.0d0
     velyr = 0.0d0
     velzl = 0.0d0
     velzr = 0.0d0
     pressl= 1.0d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl)
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 0.1d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 0.5d0
     Bxr   = 0.5d0
     Byl   = 1.0d0
     Byr   = -1.0d0
     Bzl   = 0.0d0
     Bzr   = 0.0d0
  else if (CCTK_EQUALS(shock_case,"Balsara2")) then
     rhol  = 1.0d0
     rhor  = 1.0d0
     velxl = 0.0d0
     velxr = 0.0d0
     velyl = 0.0d0
     velyr = 0.0d0
     velzl = 0.0d0
     velzr = 0.0d0
     pressl= 30.0d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl)
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 1.0d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 5.0d0
     Bxr   = 5.0d0
     Byl   = 6.0d0
     Byr   = 0.7d0
     Bzl   = 6.0d0
     Bzr   = 0.7d0
  else if (CCTK_EQUALS(shock_case,"Balsara3")) then
     rhol  = 1.0d0
     rhor  = 1.0d0
     velxl = 0.0d0
     velxr = 0.0d0
     velyl = 0.0d0
     velyr = 0.0d0
     velzl = 0.0d0
     velzr = 0.0d0
     pressl= 1000.0d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl) 
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 0.1d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 10.0d0
     Bxr   = 10.0d0
     Byl   = 7.0d0
     Byr   = 0.7d0
     Bzl   = 7.0d0
     Bzr   = 0.7d0
  else if (CCTK_EQUALS(shock_case,"Balsara4")) then
     rhol  = 1.0d0
     rhor  = 1.0d0
     velxl = 0.999d0
     velxr = -0.999d0
     velyl = 0.0d0
     velyr = 0.0d0
     velzl = 0.0d0
     velzr = 0.0d0
     pressl= 0.1d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl) 
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 0.1d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 10.0d0
     Bxr   = 10.0d0
     Byl   = 7.0d0
     Byr   = -7.0d0
     Bzl   = 7.0d0
     Bzr   = -7.0d0
  else if (CCTK_EQUALS(shock_case,"Balsara5")) then
     rhol  = 1.08d0
     rhor  = 1.0d0
     velxl = 0.4d0
     velxr = -0.45d0
     velyl = 0.3d0
     velyr = -0.2d0
     velzl = 0.2d0
     velzr = 0.2d0
     pressl= 0.95d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl) 
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 1.0d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 2.0d0
     Bxr   = 2.0d0
     Byl   = 0.3d0
     Byr   = -0.7d0
     Bzl   = 0.3d0
     Bzr   = 0.5d0
  else if (CCTK_EQUALS(shock_case,"generic_alfven")) then
     rhol  = 1.0d0
     rhor  = 0.9d0
     velxl = 0.0d0
     velxr = 0.0d0
     velyl = 0.3d0
     velyr = 0.0d0
     velzl = 0.4d0
     velzr = 0.0d0
     pressl= 5.0d0
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl) 
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr= 5.3d0
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr)
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     Bxl   = 1.0d0
     Bxr   = 1.0d0
     Byl   = 6.0d0
     Byr   = 5.0d0
     Bzl   = 2.0d0
     Bzr   = 2.0d0
  else if (CCTK_EQUALS(shock_case,"user")) then
     !an example file can be found in 
     ! Spritz_InitData/src/include
     call CCTK_FortranString(string_length,input_file,input_file_for)     
     write(*,*) "Reding Initial Data from ", input_file_for
     open(UNIT=200,FILE=input_file_for,STATUS='OLD',ACTION='READ')
     do i = 1,15
        read (200,280) dump,value(i)
        write(*,*) dump,value(i)
     enddo
     write(*,*)

280  format(a20,f20.10)
     close (200)

     rhol  = value(1)  
     rhor  = value(2)  
     pressl= value(3)  
     !epsl  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhol, pressl) 
     !!!!!!
     xrho = rhol
     xpress = pressl
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsl = xeps(1)
     !!!!!!
     pressr=value(4)  
     !epsr  = EOS_SpecificIntEnergy(Spritz_eos_handle, rhor, pressr) 
     !!!!!!
     xrho = rhor
     xpress = pressr
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsr = xeps(1)
     !!!!!!
     velxl = value(5)
     velxr = value(6)
     velyl = value(7)
     velyr = value(8)
     velzl = value(9)
     velzr = value(10)
     Bxl   = value(11)
     Bxr   = value(11)
     Byl   = value(12)
     Byr   = value(13)
     Bzl   = value(14)
     Bzr   = value(15)
     close(200)
  else
     call CCTK_WARN(0,"Shock case not recognized")
  end if

  write(warnline,750) "Input values: rhol = ",rhol,", rhor = ",rhor,&
                      ", Pl = ",pressl,", Pr = ",pressr,&
                      ", epsl = ",epsl,", epsr = ",epsr,&
                      ", vxl = ",velxl,", vxr = ",velxr,&
                      ", vyl = ",velyl,", vyr = ",velyr,&
                      ", vzl = ",velzl,", vzr = ",velzr,&
                      ", Bxl = ",Bxl,", Bxr = ",Bxr,&
                      ", Byl = ",Byl,", Byr = ",Byr,&
                      ", Bzl = ",Bzl,", Bzr = ",Bzr 
  750 format(a21,g16.7,a10,g16.7,a7,g16.7,a7,g16.7,a9,g16.7,a9,g16.7,a8,g16.7,a8,g16.7,&
             a8,g16.7,a8,g16.7,a8,g16.7,a8,g16.7,a8,g16.7,a8,g16.7,a8,g16.7,a8,g16.7,&
             a8,g16.7,a8,g16.7)
  call CCTK_INFO(warnline)

  dx = CCTK_DELTA_SPACE(1)	
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  do i=1,nx
     do j=1,ny
        do k=1,nz
  
           gxx(i,j,k) = 1.0d0
           gxy(i,j,k) = 0.0d0
           gxz(i,j,k) = 0.0d0
           gyy(i,j,k) = 1.0d0
           gyz(i,j,k) = 0.0d0
           gzz(i,j,k) = 1.0d0
           alp(i,j,k) = alp_value

           if (CCTK_EQUALS(initial_shift,"zero")) then
              betax(i,j,k) = 0.0d0
              betay(i,j,k) = 0.0d0
              betaz(i,j,k) = 0.0d0
           end if


           if (CCTK_EQUALS(shocktube_type,"diagshock")) then
              direction = x(i,j,k) - shock_xpos + &
                   y(i,j,k) - shock_ypos + z(i,j,k) - shock_zpos
           elseif (CCTK_EQUALS(shocktube_type,"2Ddiagshock")) then
              direction = x(i,j,k) - shock_xpos + &
                   y(i,j,k) - shock_ypos
           else if (CCTK_EQUALS(shocktube_type,"xshock")) then
              direction = x(i,j,k) - shock_xpos
              if (CCTK_EQUALS(initial_shift,"zero")) betax(i,j,k) = beta_value
           else if (CCTK_EQUALS(shocktube_type,"yshock")) then
              direction = y(i,j,k) - shock_ypos
              if (CCTK_EQUALS(initial_shift,"zero")) betay(i,j,k) = beta_value
           else if (CCTK_EQUALS(shocktube_type,"zshock")) then
              direction = z(i,j,k) - shock_zpos
              if (CCTK_EQUALS(initial_shift,"zero")) betaz(i,j,k) = beta_value
           else if (CCTK_EQUALS(shocktube_type,"sphere")) then
              direction = sqrt((x(i,j,k)-shock_xpos)**2+&
                   (y(i,j,k)-shock_ypos)**2+&
                   (z(i,j,k)-shock_zpos)**2)-shock_radius
           else
              call CCTK_WARN(0,"shocktube_type not recognized")
           end if

       if ( ((change_shock_direction==0).and.(direction .le. small_val)).or.&   !.lt.   0.0d0
            ((change_shock_direction==1).and.(direction .ge. (-1.0d0*small_val))) ) then !.gt.   0.0d0
          rho(i,j,k) = rhol
          eps(i,j,k) = epsl
          if (CCTK_EQUALS(shocktube_type,"xshock")) then
             velx(i,j,k) = velxl
             vely(i,j,k) = velyl
             velz(i,j,k) = velzl
             Bx(i,j,k) = Bxl
             By(i,j,k) = Byl
             Bz(i,j,k) = Bzl   
             Avec_x(i,j,k) = (By(i,j,k) * (z(i,j,k)+dz/2.0d0)) - (Bz(i,j,k) * (y(i,j,k)+dy/2.0d0))!!!!!!+0.5d0*CCTK_DELTA_SPACE(3))
             Avec_y(i,j,k) = 0.0d0 !!!!!!(Bzl * x(i,j,k)) - (Bxl * z(i,j,k))+CCTK_DELTA_SPACE(1))!need to have the same value of B at x=0
             Avec_z(i,j,k) = (Bx(i,j,k) * (y(i,j,k)+dy/2.0d0)) !!!!!! ( - (Byl * x(i,j,k)) !+0.5d0*CCTK_DELTA_SPACE(2))


          else if (CCTK_EQUALS(shocktube_type,"yshock")) then
             velx(i,j,k) = velzl !velyl
             vely(i,j,k) = velxl
             velz(i,j,k) = velyl !velzl
             Bx(i,j,k) = Bzl !Byl
             By(i,j,k) = Bxl
             Bz(i,j,k) = Byl !Bzl   
             Avec_x(i,j,k) = (By(i,j,k) * (z(i,j,k)+dz/2.0d0)) - (Bz(i,j,k) * (y(i,j,k)+dy/2.0d0))!!!!!!+0.5d0*CCTK_DELTA_SPACE(3))
             Avec_y(i,j,k) = 0.0d0 !!!!!!(Bzl * x(i,j,k)) - (Bxl * z(i,j,k))+CCTK_DELTA_SPACE(1))!need to have the same value of B at x=0
             Avec_z(i,j,k) = (Bx(i,j,k) * (y(i,j,k)+dy/2.0d0)) !!!!!! ( - (Byl * x(i,j,k)) !+0.5d0*CCTK_DELTA_SPACE(2))

          else if (CCTK_EQUALS(shocktube_type,"zshock")) then
             velx(i,j,k) = velyl
             vely(i,j,k) = velzl
             velz(i,j,k) = velxl
             Bx(i,j,k) = Byl
             By(i,j,k) = Bzl
             Bz(i,j,k) = Bxl   
             Avec_x(i,j,k) = (By(i,j,k) * (z(i,j,k)+dz/2.0d0)) - (Bz(i,j,k) * (y(i,j,k)+dy/2.0d0))!!!!!!+0.5d0*CCTK_DELTA_SPACE(3))
             Avec_y(i,j,k) = 0.0d0 !!!!!!(Bzl * x(i,j,k)) - (Bxl * z(i,j,k))+CCTK_DELTA_SPACE(1))!need to have the same value of B at x=0
             Avec_z(i,j,k) = (Bx(i,j,k) * (y(i,j,k)+dy/2.0d0)) !!!!!! ( - (Byl * x(i,j,k)) !+0.5d0*CCTK_DELTA_SPACE(2))


          else if CCTK_EQUALS(shocktube_type,"diagshock") then
             !NOTE: we want that Bxr and vxr are set along the normal to
             ! the shock front, i.e. along the main diagonal.
             ! So we rotate of 45 degrees around the z axes and then of
             ! 45 degrees around the y axes.
             velx(i,j,k) = 0.5d0*velxl - sqrt(0.5d0)*velyl - 0.5d0*velzl
             vely(i,j,k) = 0.5d0*velxl + sqrt(0.5d0)*velyl - 0.5d0*velzl
             velz(i,j,k) = sqrt(0.5d0)*(velxl+velzl)
             Bx(i,j,k) = 0.5d0*Bxl - sqrt(0.5d0)*Byl - 0.5d0*Bzl !Bxl/sqrt(3.0d0)
             By(i,j,k) = 0.5d0*Bxl + sqrt(0.5d0)*Byl - 0.5d0*Bzl !Bxl/sqrt(3.0d0)
             Bz(i,j,k) = sqrt(0.5d0)*(Bxl+Bzl) !Bxl/sqrt(3.0d0)
             Avec_x(i,j,k) = 0.5 * ((By(i,j,k) * z(i,j,k)) - (Bz(i,j,k) * y(i,j,k))) !!!!!!By(i,j,k) * z(i,j,k)
             Avec_y(i,j,k) = 0.5 * ((Bz(i,j,k) * x(i,j,k)) - (Bx(i,j,k) * z(i,j,k))) !!!!!!Bz(i,j,k) * x(i,j,k)
             Avec_z(i,j,k) = 0.5 * ((Bx(i,j,k) * y(i,j,k)) - (By(i,j,k) * x(i,j,k))) !!!!!!Bx(i,j,k) * y(i,j,k)

          else if CCTK_EQUALS(shocktube_type,"2Ddiagshock") then
             !NOTE: we want that Bxr and vxr are set along the normal to
             ! the shock front, i.e. along the diagonal on the xy plabe.
             ! So we rotate of 45 degrees around the z axes.
             velx(i,j,k) = sqrt(0.5d0)*velxl - sqrt(0.5d0)*velyl
             vely(i,j,k) = sqrt(0.5d0)*velxl + sqrt(0.5d0)*velyl
             velz(i,j,k) = velzl
             Bx(i,j,k) = sqrt(0.5d0)*Bxl - sqrt(0.5d0)*Byl
             By(i,j,k) = sqrt(0.5d0)*Bxl + sqrt(0.5d0)*Byl
             Bz(i,j,k) = Bzl
             Avec_x(i,j,k) = 0.5 * ((By(i,j,k) * z(i,j,k)) - (Bz(i,j,k) * y(i,j,k))) !!!!!!By(i,j,k) * z(i,j,k)
             Avec_y(i,j,k) = 0.5 * ((Bz(i,j,k) * x(i,j,k)) - (Bx(i,j,k) * z(i,j,k))) !!!!!!Bz(i,j,k) * x(i,j,k)
             Avec_z(i,j,k) = 0.5 * ((Bx(i,j,k) * y(i,j,k)) - (By(i,j,k) * x(i,j,k))) !!!!!!Bx(i,j,k) * y(i,j,k)
          end if
       else
          rho(i,j,k) = rhor
          eps(i,j,k) = epsr
          if (CCTK_EQUALS(shocktube_type,"xshock")) then
             velx(i,j,k) = velxr
             vely(i,j,k) = velyr
             velz(i,j,k) = velzr
             Bx(i,j,k) = Bxr
             By(i,j,k) = Byr
             Bz(i,j,k) = Bzr   
             Avec_x(i,j,k) = (By(i,j,k) * (z(i,j,k)+dz/2.0d0)) - (Bz(i,j,k) * (y(i,j,k)+dy/2.0d0))!!!!!!+0.5d0*CCTK_DELTA_SPACE(3))
             Avec_y(i,j,k) = 0.0d0 !!!!!!(Bzl * x(i,j,k)) - (Bxl * z(i,j,k))+CCTK_DELTA_SPACE(1))!need to have the same value of B at x=0
             Avec_z(i,j,k) = (Bx(i,j,k) * (y(i,j,k)+dy/2.0d0)) !!!!!! ( - (Byl * x(i,j,k)) !+0.5d0*CCTK_DELTA_SPACE(2))

          else if (CCTK_EQUALS(shocktube_type,"yshock")) then
             velx(i,j,k) = velzr !velyr
             vely(i,j,k) = velxr
             velz(i,j,k) = velyr !velzr
             Bx(i,j,k) = Bzr !Byr
             By(i,j,k) = Bxr
             Bz(i,j,k) = Byr !Bzr  
             Avec_x(i,j,k) = (By(i,j,k) * (z(i,j,k)+dz/2.0d0)) - (Bz(i,j,k) * (y(i,j,k)+dy/2.0d0))!!!!!!+0.5d0*CCTK_DELTA_SPACE(3))
             Avec_y(i,j,k) = 0.0d0 !!!!!!(Bzl * x(i,j,k)) - (Bxl * z(i,j,k))+CCTK_DELTA_SPACE(1))!need to have the same value of B at x=0
             Avec_z(i,j,k) = (Bx(i,j,k) * (y(i,j,k)+dy/2.0d0)) !!!!!! ( - (Byl * x(i,j,k)) !+0.5d0*CCTK_DELTA_SPACE(2))

          else if (CCTK_EQUALS(shocktube_type,"zshock")) then
             velx(i,j,k) = velyr
             vely(i,j,k) = velzr
             velz(i,j,k) = velxr
             Bx(i,j,k) = Byr
             By(i,j,k) = Bzr
             Bz(i,j,k) = Bxr   
             Avec_x(i,j,k) = (By(i,j,k) * (z(i,j,k)+dz/2.0d0)) - (Bz(i,j,k) * (y(i,j,k)+dy/2.0d0))!!!!!!+0.5d0*CCTK_DELTA_SPACE(3))
             Avec_y(i,j,k) = 0.0d0 !!!!!!(Bzl * x(i,j,k)) - (Bxl * z(i,j,k))+CCTK_DELTA_SPACE(1))!need to have the same value of B at x=0
             Avec_z(i,j,k) = (Bx(i,j,k) * (y(i,j,k)+dy/2.0d0)) !!!!!! ( - (Byl * x(i,j,k)) !+0.5d0*CCTK_DELTA_SPACE(2))

          else if CCTK_EQUALS(shocktube_type,"diagshock") then
             !NOTE: we want that Bxr and vxr are set along the normal to
             ! the shock front, i.e. along the main diagonal.
             ! So we rotate of 45 degrees around the z axes and then of
             ! 45 degrees around the y axes.
             velx(i,j,k) = 0.5d0*velxr - sqrt(0.5d0)*velyr - 0.5d0*velzr
             vely(i,j,k) = 0.5d0*velxr + sqrt(0.5d0)*velyr - 0.5d0*velzr
             velz(i,j,k) = sqrt(0.5d0)*(velxr+velzr)
             Bx(i,j,k) = 0.5d0*Bxr - sqrt(0.5d0)*Byr - 0.5d0*Bzr !Bxr/sqrt(3.0d0)
             By(i,j,k) = 0.5d0*Bxr + sqrt(0.5d0)*Byr - 0.5d0*Bzr !Bxr/sqrt(3.0d0)
             Bz(i,j,k) = sqrt(0.5d0)*(Bxr+Bzr) !Bxr/sqrt(3.0d0)
             Avec_x(i,j,k) = 0.5 * ((By(i,j,k) * z(i,j,k)) - (Bz(i,j,k) * y(i,j,k))) !!!!!!By(i,j,k) * z(i,j,k)
             Avec_y(i,j,k) = 0.5 * ((Bz(i,j,k) * x(i,j,k)) - (Bx(i,j,k) * z(i,j,k))) !!!!!!Bz(i,j,k) * x(i,j,k)
             Avec_z(i,j,k) = 0.5 * ((Bx(i,j,k) * y(i,j,k)) - (By(i,j,k) * x(i,j,k))) !!!!!!Bx(i,j,k) * y(i,j,k)

          else if CCTK_EQUALS(shocktube_type,"2Ddiagshock") then
             !NOTE: we want that Bxr and vxr are set along the normal to
             ! the shock front, i.e. along the diagonal on the xy plabe.
             ! So we rotate of 45 degrees around the z axes.
             velx(i,j,k) = sqrt(0.5d0)*velxr - sqrt(0.5d0)*velyr
             vely(i,j,k) = sqrt(0.5d0)*velxr + sqrt(0.5d0)*velyr
             velz(i,j,k) = velzr
             Bx(i,j,k) = sqrt(0.5d0)*Bxr - sqrt(0.5d0)*Byr
             By(i,j,k) = sqrt(0.5d0)*Bxr + sqrt(0.5d0)*Byr
             Bz(i,j,k) = Bzr
             Avec_x(i,j,k) = 0.5 * ((By(i,j,k) * z(i,j,k)) - (Bz(i,j,k) * y(i,j,k))) !!!!!!By(i,j,k) * z(i,j,k)
             Avec_y(i,j,k) = 0.5 * ((Bz(i,j,k) * x(i,j,k)) - (Bx(i,j,k) * z(i,j,k))) !!!!!!Bz(i,j,k) * x(i,j,k)
             Avec_z(i,j,k) = 0.5 * ((Bx(i,j,k) * y(i,j,k)) - (By(i,j,k) * x(i,j,k))) !!!!!!Bx(i,j,k) * y(i,j,k)
          end if
       end if

        call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

        if (CCTK_EQUALS(Spritz_eos_type,"Polytype")) then
           call prim2conpolytype_mhd(Spritz_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                tau(i,j,k),rho(i,j,k),&
                velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                Bx(i,j,k), By(i,j,k), Bz(i,j,k), &
                dBx(i,j,k), dBy(i,j,k), dBz(i,j,k), alp(i,j,k),&
                betax(i,j,k),betay(i,j,k),betaz(i,j,k))
        else
           call prim2con_mhd(Spritz_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                tau(i,j,k),rho(i,j,k),&
                velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                Bx(i,j,k), By(i,j,k), Bz(i,j,k), &
                dBx(i,j,k), dBy(i,j,k), dBz(i,j,k),alp(i,j,k),&
                betax(i,j,k),betay(i,j,k),betaz(i,j,k), &
                Spritz_coldtab_rhomin,Spritz_coldtab_press_lo, &
                Spritz_coldtab_eps_lo,Spritz_coldtab_gamma_lo, &
                Spritz_coldtab_gammath,Spritz_coldtab_thfac)
        end if

        rho_p(i,j,k) = rho(i,j,k)
        dens_p(i,j,k) = dens(i,j,k)

        scon_p(i,j,k,1)   = sx(i,j,k)
        scon_p(i,j,k,2)   = sy(i,j,k)
        scon_p(i,j,k,3)   = sz(i,j,k)
        vel_p(i,j,k,1) = vel(i,j,k,1)
        vel_p(i,j,k,2) = vel(i,j,k,2)
        vel_p(i,j,k,3) = vel(i,j,k,3)

        tau_p(i,j,k)    = tau(i,j,k)
        press_p(i,j,k)  = press(i,j,k)
        eps_p(i,j,k)    = eps(i,j,k)

        Bx_p(i,j,k) = Bx(i,j,k)
        By_p(i,j,k) = By(i,j,k)
        Bz_p(i,j,k) = Bz(i,j,k)
        dBx_p(i,j,k) = dBx(i,j,k)
        dBy_p(i,j,k) = dBy(i,j,k)
        dBz_p(i,j,k) = dBz(i,j,k)
        Avec_x_p(i,j,k) = Avec_x(i,j,k)
        Avec_y_p(i,j,k) = Avec_y(i,j,k)
        Avec_z_p(i,j,k) = Avec_z(i,j,k)

        if (timelevels==3) then
           rho_p_p(i,j,k) = rho(i,j,k)
           dens_p_p(i,j,k) = dens(i,j,k)

           scon_p_p(i,j,k,1)   = sx(i,j,k)
           scon_p_p(i,j,k,2)   = sy(i,j,k)
           scon_p_p(i,j,k,3)   = sz(i,j,k)
           vel_p_p(i,j,k,1) = vel(i,j,k,1)
           vel_p_p(i,j,k,2) = vel(i,j,k,2)
           vel_p_p(i,j,k,3) = vel(i,j,k,3)

           tau_p_p(i,j,k)    = tau(i,j,k)
           press_p_p(i,j,k)  = press(i,j,k)
           eps_p_p(i,j,k)    = eps(i,j,k)

           Bx_p_p(i,j,k) = Bx(i,j,k)
           By_p_p(i,j,k) = By(i,j,k)
           Bz_p_p(i,j,k) = Bz(i,j,k)
           dBx_p_p(i,j,k) = dBx(i,j,k)
           dBy_p_p(i,j,k) = dBy(i,j,k)
           dBz_p_p(i,j,k) = dBz(i,j,k)
           Avec_x_p_p(i,j,k) = Avec_x(i,j,k)
           Avec_y_p_p(i,j,k) = Avec_y(i,j,k)
           Avec_z_p_p(i,j,k) = Avec_z(i,j,k)
        end if

     enddo
  enddo
enddo

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0
  Bxrhs = 0.0d0
  Byrhs = 0.0d0
  Bzrhs = 0.0d0
  Avec_x_rhs = 0.0d0
  Avec_y_rhs = 0.0d0
  Avec_z_rhs = 0.0d0

  Ex=0.0d0
  Ey=0.0d0
  Ez=0.0d0

  return
  
end subroutine Spritz_shocktube



subroutine change_alp_beta(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT:: nx,ny,nz,i,j,k  

  character(len=100) warnline

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  if (CCTK_EQUALS(initial_shift,"none")) then
     call CCTK_WARN(0,"You have to set ADMBase::initial_shift=zero")
  end if

  do i=1,nx
     do j=1,ny
        do k=1,nz


           if (use_alp .eq. 1) then
              alp(i,j,k) = alp_value
              alp_p(i,j,k) = alp_value
              if (CCTK_EQUALS(shocktube_type,"xshock")) then
                 betax(i,j,k) = beta_value
                 betax_p(i,j,k) = beta_value
              else if (CCTK_EQUALS(shocktube_type,"yshock")) then
                 betay(i,j,k) = beta_value
                 betay_p(i,j,k) = beta_value
              else if (CCTK_EQUALS(shocktube_type,"zshock")) then
                 betaz(i,j,k) = beta_value
                 betaz_p(i,j,k) = beta_value
              else if (CCTK_EQUALS(shocktube_type,"2Ddiagshock")) then
                 betax(i,j,k) = beta_value*sqrt(0.5d0)
                 betax_p(i,j,k) = beta_value*sqrt(0.5d0)
                 betay(i,j,k) = beta_value*sqrt(0.5d0)
                 betay_p(i,j,k) = beta_value*sqrt(0.5d0)
              else if (CCTK_EQUALS(shocktube_type,"diagshock")) then
                 betax(i,j,k) = beta_value/sqrt(3.0d0)
                 betax_p(i,j,k) = beta_value/sqrt(3.0d0)
                 betay(i,j,k) = beta_value/sqrt(3.0d0)
                 betay_p(i,j,k) = beta_value/sqrt(3.0d0)
                 betaz(i,j,k) = beta_value/sqrt(3.0d0)
                 betaz_p(i,j,k) = beta_value/sqrt(3.0d0)
              else
                 call CCTK_WARN(0,"wrong shocktube_type parameter")
              end if
              if (timelevels==3) then
                 alp_p_p(i,j,k) = alp_value
                 betax_p_p(i,j,k) = beta_value
                 betay_p_p(i,j,k) = beta_value
                 betaz_p_p(i,j,k) = beta_value
              end if
           end if

        end do
     end do
  end do

  call CCTK_INFO("Values of the lapse and shift have been changed")
  write(warnline,'(a40,g16.7)') 'The new value for the lapse is: ',alp_value
  call CCTK_INFO(warnline)
  write(warnline,'(a40,g16.7)') 'The new value for the shift is: ',beta_value
  call CCTK_INFO(warnline)

end subroutine change_alp_beta
