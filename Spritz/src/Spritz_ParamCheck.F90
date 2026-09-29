 /*@@
   @file      Spritz_ParamCheck.F90
   @date      Sat Feb  9 23:48:01 2002
   @author
   @desc
   Parameter checking routine.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Spritz_ParamCheck
   @date       Sat Feb  9 23:48:43 2002
   @author     Ian Hawke
   @desc
   Checks the parameters.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine Spritz_ParamCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  character*512 :: message
  CCTK_INT      :: ierr

  if (Spritz_stencil > minval(cctk_nghostzones)) then
    call CCTK_PARAMWARN("The stencil is larger than the number of ghost zones. Answer will be dependent on processor number...")
  end if

  if (Spritz_stencil < 3) then
    call CCTK_PARAMWARN("Spritz requires at least 3 ghostzones. Please fix the value of Spritz_stencil")
  end if
  
  if (Spritz_Avec_Staggered.eq.0) then
     call CCTK_PARAMWARN("Spritz only with a staggered vector potential. Please set Spritz_Avec_Staggered to yes")
  end if


  if (CCTK_EQUALS(recon_method,"ppm")) then
    if (CCTK_EQUALS(ppm_flatten,"stencil_3").and.(Spritz_stencil < 3)) then
      call CCTK_PARAMWARN("The stencil size must be at least 3 to use PPM reconstruction with the stencil-3 variant of the flattening procedure")
    else if (CCTK_EQUALS(ppm_flatten,"stencil_4").and.(Spritz_stencil < 4)) then
      call CCTK_PARAMWARN("The stencil size must be at least 4 to use PPM reconstruction with the stencil-4 (original) flattening procedure")
    end if
  end if



  if (CCTK_EQUALS(Spritz_eos_table,"2D_Polytrope").and.&
       (.not.CCTK_EQUALS(Spritz_eos_type,"Polytype"))) then
    call CCTK_PARAMWARN("When using the 2D_Polytrope EOS you need to set eos_type to Polytype")
  end if

  if (CCTK_EQUALS(Spritz_eos_table,"Ideal_Fluid").and.&
       (.not.CCTK_EQUALS(Spritz_eos_type,"General"))) then
    call CCTK_PARAMWARN("When using the Ideal_Fluid EOS you need to set eos_type to General")
  end if

  if (CCTK_EQUALS(Y_e_evolution_method,"Spritz")) then
    evolve_Y_e = 1
  else
    evolve_Y_e = 0
  end if

  if (CCTK_EQUALS(temperature_evolution_method,"Spritz")) then
    GetTemp_FromThermo = 1
  else
    GetTemp_FromThermo = 0
  end if

  if (CCTK_EQUALS(Spritz_eos_table,"nuc_eos") .and.&
      (.not.CCTK_EQUALS(con2prim_mhd_method,"Palenzuela"))) then
    call CCTK_PARAMWARN("You have to set con2prim_mhd_method to Palenzuela when using nuc_eos.")
  end if
  
  if (CCTK_EQUALS(Spritz_eos_table,"nuc_eos") .and.&
      (Spritz_Y_e_max.eq.1.0)) then
    call CCTK_PARAMWARN("Spritz_Y_e_max has been set to 1, are you sure this is the maximum Ye available in your EOS?")
  end if


  ! Need for checking tabulated eos parameters:
    ! When temperature_evolution_method is Spritz:
      ! nuc_eos must be used
        ! If Spritz_EvolveT_at_time, just evolve after that time
        ! Otherwise, always evolve temperature.
      ! Y_e must be advected
      ! If you choose to constrain temperature (in the whole sim or just the beginning) check uniform temperature.
    ! When temperature_evolution_method is NOT Spritz:
      ! Error if trying to use nuc_eos
    ! If constraining temperature:
      ! you have to use nuc_eos
      ! if Spritz_EvolveT_at_time > 0, there's an error

  ! Thus, requiring GetTemp_FromThermo.eq.1 or Spritz_eos_table.eq.4 is exactly the same

  if (GetTemp_FromThermo.eq.0) then
    call CCTK_INFO("Temperature won't be changed by Spritz.")
    if(CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) then
      call CCTK_PARAMWARN("But you can't use tabulated EoS without dealing with temperature.")
    end if
  else
    if(.not.CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) then
      call CCTK_PARAMWARN("But you have to use nuc_eos as Spritz_eos_table, otherwise it makes no sense.")
    end if
    if(evolve_Y_e.eq.0) then
      call CCTK_PARAMWARN("Are you sure that you don't want Y_e to be advected???")
    end if
    if(Spritz_EvolveT_at_time.gt.0.0d0) then
      write(message,*) "Temperature will be kept fixed till time ", Spritz_EvolveT_at_time
      call CCTK_INFO (message)
    else
      call CCTK_INFO("Temperature will be dealt with using EOS_Omni and nuc_eos.")
    end if
    if(Spritz_EvolveT_at_time.eq.-2.0d0.or.Spritz_EvolveT_at_time.gt.0.0d0) then
      call CCTK_INFO("Temperature will be constrained to the initial value. I'll make sure that the initial value is constant. Otherwise, I'd have to advect it, but that feature hasn't been implemented yet.")
    end if
  end if

  if(Spritz_EvolveT_at_time.eq.-2.0d0.and..not.(CCTK_EQUALS(Spritz_eos_table,"nuc_eos"))) then
    call CCTK_PARAMWARN("You have to use nuc_eos tables in order to be able to constrain temperature.")
  end if

  if(CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) then
    if((Spritz_EvolveT_at_time.ge.-1.0d-6 .and. cctk_time.lt.Spritz_EvolveT_at_time) .or.&
        Spritz_EvolveT_at_time.eq.-2.0d0) then
      Constrain_Temp_now = 1
      call CCTK_INFO("At beginning, I'm CONSTRAINING temperature.")
      call CCTK_ParameterSet(ierr,"Spritz_Constrain_Temperature",CCTK_THORNSTRING,"yes")
      if(ierr.ne.0) then
        call CCTK_WARN(0,"Error in setting Spritz_Constrain_Temperature.")
      end if
    else
      Constrain_Temp_now = 0
      call CCTK_INFO("At beginning, I'm NOT constraining temperature.")
      call CCTK_ParameterSet(ierr,"Spritz_Constrain_Temperature",CCTK_THORNSTRING,"no")
      if(ierr.ne.0) then
        call CCTK_WARN(0,"Error in setting Spritz_Constrain_Temperature.")
      end if
    end if
  end if


  if (.not.(CCTK_EQUALS(metric_type, "Physical").or.&
       CCTK_EQUALS(metric_type, "Static Conformal"))) then
    call CCTK_PARAMWARN("Spritz only knows how to deal with physical metric type at the moment. Complain to the maintainers...")
  end if


  if (.not.CCTK_EQUALS(metric_type, "Static Conformal")) then
    call CCTK_INFO("Spritz will use the physical metric and set conformal_state=0.")
    conformal_state = 0
  end if


  if (CCTK_EQUALS(initial_shift, "none")) then
     call CCTK_PARAMWARN("ADMBase::initial_shift cannot be set to none when using mhd. I need it to compute the energy-momentum tensor.")
  end if

  if ((Spritz_Avec_Staggered.eq.0).and.(Spritz_mhd_fluxCD_Lorenz.ne.0)) then
     call CCTK_PARAMWARN("Spritz only allows (Generalized) Lorenz Gauge when using staggered Avec.")
  end if
  
  if ((approximation_order.eq.4).and.(Spritz_stencil < 4)) then
  	call CCTK_PARAMWARN("If the order of approximation of the Riemann solver is 4, then the Spritz_stencil must be at least 4")
  end if
  
  if ((approximation_order.eq.6).and.(Spritz_stencil < 5)) then
  	call CCTK_PARAMWARN("If the order of approximation of the Riemann solver is 6, then the Spritz_stencil must be at least 5")
  end if


  return

end subroutine Spritz_ParamCheck

/*@@
  @routine    Spritz_TemperatureCheck
  @date       22-07-2019 14:41
  @author     Lorenzo Sala
  @desc
  Checks temperature consistency.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

@@*/

subroutine Spritz_TemperatureCheck(CCTK_ARGUMENTS)

 implicit none

 DECLARE_CCTK_ARGUMENTS
 DECLARE_CCTK_PARAMETERS
 DECLARE_CCTK_FUNCTIONS

 integer :: i,j,k
 CCTK_REAL :: firsttempval
 character(len=300) warnline

  if(Constrain_Temp_now.eq.1) then

    call CCTK_INFO("You want to constrain temperature to the initial value. Thus, I make sure that you have set up a uniform initial value.")

    firsttempval = temperature(1,1,1)
    do k = 1, cctk_lsh(3)  !!!1, nz
      do j = 1, cctk_lsh(2) !!!1, ny
        do i = 1, cctk_lsh(1) !!!1, nx
          if(temperature(i,j,k).ne.firsttempval) then
            call CCTK_WARN(0, "Wait, wait, you don't have uniform temperature but you want to constrain it to a uniform value")
          end if
        end do
      end do
    end do
    write(warnline,'(a20,g16.7,a20)') 'Ok, temp is', firsttempval, "MeV."
    call CCTK_INFO(warnline)
  end if

 return

end subroutine Spritz_TemperatureCheck
