 /*@@
   @file      Spritz_WENOZ.F90
   @date      
   @author    Beatrice Giudici
   @desc 
   WENO-Z reconstruction routine. Based on the one implemented in GRHydro code
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


/*@@
	@routine	wenoz_reconstruct
	@date	
	@author	Beatrice Giudici
	@desc
	Performs the reconstruction of a given variable using 5th order WENO-Z method
	(Borges et al., ''An improved weighted essentially non-oscillatory scheme for 
	hyperbolic conservation laws'', 2008)
	@enddesc
	@calls
	@calledby
@@*/

subroutine wenoz_reconstruct_1D(N, v, vminus, vplus, trivial_rp)

	use Spritz_Scalars
	
	implicit none

	DECLARE_CCTK_PARAMETERS
	DECLARE_CCTK_FUNCTIONS
	
	CCTK_INT :: N, i
	CCTK_REAL :: tau
	CCTK_REAL, dimension(N) :: v, vminus, vplus
	CCTK_REAL, dimension(2,3) :: alphaZ, omegaZ !! 1->plus, 2->minus
	CCTK_REAL, dimension(3) :: beta
	logical, dimension(N) :: trivial_rp
	CCTK_INT :: HO_stc
	

!! If high order approximation is enabled, we need a larger stencil
	if (high_order .ne. 0) then
		if (approximation_order == 4) HO_stc = 1
		if (approximation_order == 6) HO_stc = 2
	else
		HO_stc = 0
	end if
	
	trivial_rp = .false.
	
	do i = Spritz_stencil - HO_stc, N - Spritz_stencil + 1 + HO_stc
!! Calculating the smoothness indicators (Borges et al. 2008)
!! Note: these quantities are symmetric in the reflection of the stencil
		beta(1) = (13.0d0/12.0d0) * (v(i-2) - 2.0d0*v(i-1) + v(i))**2 &
			   + (1.0d0/4.0d0) * (v(i-2) - 4.0d0*v(i-1) + 3.0d0*v(i))**2
		
		beta(2) = (13.0d0/12.0d0) * (v(i-1) - 2.0d0*v(i) + v(i+1))**2 &
			   + (1.0d0/4.0d0) * (v(i-1) - v(i+1))**2
		
		beta(3) = (13.0d0/12.0d0) * (v(i) - 2.0d0*v(i+1) + v(i+2))**2 &
			   + (1.0d0/4.0d0) * (3.0d0*v(i) - 4.0d0*v(i+1) + v(i+2))**2

!! Computing the weights according to WENO-Z algorithm
		tau = abs(beta(1) - beta(3))

	!! Unnormalized weights (the optimal weights are chosen according to Del Zanna et al. 2007)
		alphaZ(1,1) = (1.0d0/16.0d0) * (1.0d0 + tau / (beta(1) + weno_eps))
		alphaZ(1,2) = (10.0d0/16.0d0) * (1.0d0 + tau / (beta(2) + weno_eps))
		alphaZ(1,3) = (5.0d0/16.0d0) * (1.0d0 + tau / (beta(3) + weno_eps))
		
	!! Weights for the reconstruction
		omegaZ(1,1) = alphaZ(1,1) / (alphaZ(1,1) + alphaZ(1,2) + alphaZ(1,3))
		omegaZ(1,2) = alphaZ(1,2) / (alphaZ(1,1) + alphaZ(1,2) + alphaZ(1,3))
		omegaZ(1,3) = alphaZ(1,3) / (alphaZ(1,1) + alphaZ(1,2) + alphaZ(1,3))
		
!! Reconstruct the function at the right interface (Del Zanna et al. 2007)
		vplus(i) = omegaZ(1,1)/8.0d0 * (3.0d0*v(i-2) - 10.0d0*v(i-1) + 15.0d0*v(i)) &
			    + omegaZ(1,2)/8.0d0 * (-1.0d0*v(i-1) + 6.0d0*v(i) + 3.0d0*v(i+1)) &
			    + omegaZ(1,3)/8.0d0 * (3.0d0*v(i) + 6.0d0*v(i+1) - 1.0d0*v(i+2))
			    
	!! Unnormalized weights
		alphaZ(2,1) = (5.0d0/16.0d0) * (1.0d0 + tau / (beta(1) + weno_eps))
		alphaZ(2,2) = (10.0d0/16.0d0) * (1.0d0 + tau / (beta(2) + weno_eps))
		alphaZ(2,3) = (1.0d0/16.0d0) * (1.0d0 + tau / (beta(3) + weno_eps))
		
	!! Weights for the reconstruction (Del Zanna et al. 2007)
		omegaZ(2,1) = alphaZ(2,1) / (alphaZ(2,1) + alphaZ(2,2) + alphaZ(2,3))
		omegaZ(2,2) = alphaZ(2,2) / (alphaZ(2,1) + alphaZ(2,2) + alphaZ(2,3))
		omegaZ(2,3) = alphaZ(2,3) / (alphaZ(2,1) + alphaZ(2,2) + alphaZ(2,3))
			    
!! Reconstruct the function at the left interface, simply by swapping the points of the sub-stencils
!! and the corresponding weights (Del Zanna et al. 2007)
		vminus(i) = omegaZ(2,3)/8.0d0 * (3.0d0*v(i+2) - 10.0d0*v(i+1) + 15.0d0*v(i)) &
			     + omegaZ(2,2)/8.0d0 * (-1.0d0*v(i+1) + 6.0d0*v(i) + 3.0d0*v(i-1)) &
			     + omegaZ(2,1)/8.0d0 * (3.0d0*v(i) + 6.0d0*v(i-1) - 1.0d0*v(i-2))
	end do
	
	return
	
end subroutine wenoz_reconstruct_1D
