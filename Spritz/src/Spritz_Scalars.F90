 /*@@
   @file      Spritz_Scalars.F90
   @date      Mon Feb 25 11:11:23 2002
   @author    
   @desc 
   Module containing various scalars to avoid having to use
   CCTK_EQUALS on keywords at every cell (slow).
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 module Spritz_Scalars

   implicit none

   LOGICAL :: MINMOD, MINMOD2, MINMOD3, MCB, MC1, MC2, SUPERBEE, PPM3, PPM4, VANLEERBG
   LOGICAL :: ANALYTICAL
   LOGICAL :: FAST

   CCTK_INT :: EOS_RiemannCallPlus, EOS_RiemannCallMinus, EOS_RoeAverageCall
   CCTK_INT :: EOS_Con2PrimCall
   CCTK_INT :: EOS_Prim2ConBndCallPlus, EOS_Prim2ConBndCallMinus
   CCTK_INT :: EOS_Prim2ConCellsCall

!   CCTK_REAL :: eosgeneral_pmin
   
 end module Spritz_Scalars
