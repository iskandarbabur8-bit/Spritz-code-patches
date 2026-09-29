#include "cctk.h"

SUBROUTINE zrhqr(a,rtr,rti) 
  USE nrutil, ONLY : assert_eq,nrerror 
  USE nrmodule, ONLY : balanc,hqr
  IMPLICIT NONE 
  CCTK_REAL, DIMENSION(:), INTENT(IN) :: a 
  CCTK_REAL, DIMENSION(:), INTENT(OUT) :: rtr,rti
  !Find all the roots of a polynomial with real coefficients, 
  !\sum_{i=1}^{M+1}a(i)x^{i-1}, given the array of M +1 coefficients a.
  !The method is to construct an upper Hessenberg matrix whose eigenvalues are 
  !the desired roots,and then use the routines balanc and hqr. 
  !The real and imaginary parts of the M roots are returned in 
  !rtr and rti respectively.

  !(C) 1986-1996 Numerical Recipes Software

  CCTK_INT :: k,m 
!  CCTK_INT, DIMENSION(size(rtr)) :: indx 
  CCTK_REAL, DIMENSION(size(a)-1,size(a)-1) :: hess 
  m=assert_eq(size(rtr),size(rti),size(a)-1,'zrhqr')
  if (a(m+1) == 0.0) call & 
       nrerror('zrhqr: Last value of array a must not be 0') 
  hess(1,:)=-a(m:1:-1)/a(m+1) !Construct the matrix. 
  hess(2:m,:)=0.0 
  do k=1,m-1 
     hess(k+1,k)=1.0 
  end do
  call balanc(hess) !Find its eigenvalues. 
  call hqr(hess,rtr,rti) 
  !call indexx(rtr,indx) !Sort roots by their real parts. 
  !rtr=rtr(indx) 
  !rti=rti(indx) 
END SUBROUTINE zrhqr
