#include "cctk.h"

MODULE nrutil
  !(C) 1986-1996 Numerical Recipes Software
  IMPLICIT NONE
  CCTK_INT, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_rv, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE

  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE

  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE

  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r
  END INTERFACE

  INTERFACE diagadd
     MODULE PROCEDURE diagadd_rv,diagadd_r
  END INTERFACE  

  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_i
     !MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE

  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_i
     !MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
  END INTERFACE


contains
  
  SUBROUTINE swap_i(a,b)
    CCTK_INT, INTENT(INOUT) :: a,b
    CCTK_INT :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  SUBROUTINE swap_r(a,b)
    CCTK_REAL, INTENT(INOUT) :: a,b
    CCTK_REAL :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r

  SUBROUTINE swap_rv(a,b)
    CCTK_REAL, DIMENSION(:), INTENT(INOUT) :: a,b
    CCTK_REAL, DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv

  SUBROUTINE masked_swap_rs(a,b,mask)
    CCTK_REAL, INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    CCTK_REAL :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs

  SUBROUTINE masked_swap_rv(a,b,mask)
    CCTK_REAL, DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    CCTK_REAL, DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv

  SUBROUTINE masked_swap_rm(a,b,mask)
    CCTK_REAL, DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    CCTK_REAL, DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm

  FUNCTION imaxloc_r(arr)
    CCTK_REAL, DIMENSION(:), INTENT(IN) :: arr
    CCTK_INT :: imaxloc_r
    CCTK_INT, DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r

  FUNCTION imaxloc_i(iarr)
    CCTK_INT, DIMENSION(:), INTENT(IN) :: iarr
    CCTK_INT, DIMENSION(1) :: imax
    CCTK_INT :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i

  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    CCTK_INT, INTENT(IN) :: n1,n2
    CCTK_INT :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    CCTK_INT, INTENT(IN) :: n1,n2,n3
    CCTK_INT :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3

  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    CCTK_INT, INTENT(IN) :: n1,n2,n3,n4
    CCTK_INT :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    CCTK_INT, DIMENSION(:), INTENT(IN) :: nn
    CCTK_INT :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

  FUNCTION outerprod_r(a,b)
    CCTK_REAL, DIMENSION(:), INTENT(IN) :: a,b
    CCTK_REAL, DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r

  SUBROUTINE diagadd_rv(mat,diag)
    CCTK_REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
    CCTK_REAL, DIMENSION(:), INTENT(IN) :: diag
    CCTK_INT :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  END SUBROUTINE diagadd_rv
  !BL
  SUBROUTINE diagadd_r(mat,diag)
    CCTK_REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
    CCTK_REAL, INTENT(IN) :: diag
    CCTK_INT :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  END SUBROUTINE diagadd_r

  FUNCTION upper_triangle(j,k,extra)
    CCTK_INT, INTENT(IN) :: j,k
    CCTK_INT, OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: upper_triangle
    CCTK_INT :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle

  FUNCTION arth_r(first,increment,n)
    CCTK_REAL, INTENT(IN) :: first,increment
    CCTK_INT, INTENT(IN) :: n
    CCTK_REAL, DIMENSION(n) :: arth_r
    CCTK_INT :: k,k2
    CCTK_REAL :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r
  !BL
  FUNCTION arth_i(first,increment,n)
    CCTK_INT, INTENT(IN) :: first,increment,n
    CCTK_INT, DIMENSION(n) :: arth_i
    CCTK_INT :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  FUNCTION outerdiff_r(a,b)
    CCTK_REAL, DIMENSION(:), INTENT(IN) :: a,b
    CCTK_REAL, DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
  !BL
  FUNCTION outerdiff_i(a,b)
    CCTK_INT, DIMENSION(:), INTENT(IN) :: a,b
    CCTK_INT, DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i



END MODULE nrutil
