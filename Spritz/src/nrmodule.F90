#include "cctk.h"

module nrmodule
  implicit none
contains
  SUBROUTINE lubksb(a,indx,b)
    !Solves the set of N linear equations A·X=. Here the N×N matrix a is input,
    !not as the original matrix A, but rather as its LU decomposition,
    !determined by the routine ludcmp. indx is input as the permutation vector 
    !of length N returned by ludcmp. b is input as the right-hand-side vector B,
    !also of length N, and returns with the solution vector X.
    !a and indx are not modified by this routine and can be left 
    !in place for successive calls with different right-hand sides b. 
    !This routine takes into account the possibility that b will begin 
    !with many zero elements, so it is efficient for use in matrix inversion.

    !(C) 1986-1996 Numerical Recipes Software

    USE nrutil, ONLY : assert_eq 
    IMPLICIT NONE 
    CCTK_REAL, DIMENSION(:,:), INTENT(IN) :: a 
    CCTK_INT, DIMENSION(:), INTENT(IN) :: indx 
    CCTK_REAL, DIMENSION(:), INTENT(INOUT) :: b
    CCTK_INT :: i,n,ii,ll 
    CCTK_REAL :: summ 
    n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    ii=0 
    !When ii is set to a positive value, it will become the index 
    !of the first nonvanishing element of b. 
    !We now do the forward substitution, equation (2.3.6). 
    !The only new wrinkle is to unscramble the permutation as we go. 
    do i=1,n 
       ll=indx(i) 
       summ=b(ll)
       b(ll)=b(i) 
       if (ii /= 0) then 
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1)) 
       else if (summ /= 0.0) then 
          ii=i 
          !A nonzero element was encountered, so from now on 
          !we will have to do the dot product above.
       end if
       b(i)=summ 
    end do
    do i=n,1,-1 !Now we do the backsubstitution,equation (2.3.7). 
       b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i) 
    end do
  END SUBROUTINE lubksb

  SUBROUTINE ludcmp(a,indx,d)
    !Given an N×N input matrix a, 
    !this routine replaces it by the LU decomposition
    !of a rowwise permutation of itself.
    !On output, a is arranged as in equation (2.3.14);
    !indx is an output vector of length N that records the row permutation
    !effected by the partial pivoting;
    !d is output as ±1 depending on whether the number of row 
    !interchanges was even or odd, respectively.
    !This routine is used in combination with lubksb to solve linear equations 
    !or invert a matrix.

    !(C) 1986-1996 Numerical Recipes Software

    USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap 
    IMPLICIT NONE 
    CCTK_REAL, DIMENSION(:,:), INTENT(INOUT) :: a 
    CCTK_INT, DIMENSION(:), INTENT(OUT) :: indx 
    CCTK_REAL, INTENT(OUT) :: d
    CCTK_REAL, DIMENSION(size(a,1)) :: vv 
    !vv stores the implicit scaling of each row. 
    CCTK_REAL, PARAMETER :: TINY=1.0e-20  !A small number. 
    CCTK_INT :: j,n,imax 
    n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
    d=1.0 !No row interchanges yet. 
    vv=maxval(abs(a),dim=2) 
    !Loop over rows to get the implicit scaling information
    if (any(vv == 0.0)) then
       write(*,*) 'matrix=',a
       call nrerror('singular matrix in ludcmp') !There is a row of zeros. 
    end if
    vv=1.0/vv !Save the scaling.
    do j=1,n 
       imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row. 
       if (j /= imax) then !Do we need to interchange rows? 
          call swap(a(imax,:),a(j,:)) !Yes,do so... 
          d=-d !...and change the parity of d 
          vv(imax)=vv(j) !Also interchange the scale factor. 
       end if
       indx(j)=imax 
       if (a(j,j) == 0.0) a(j,j)=TINY 
       !If the pivot element is zero the matrix is singular 
       !(at least to the precision of the algorithm).
       !For some applications on singular matrices, it is desirable 
       !to substitute TINY for zero. 
       a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide y the pivot element. 
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n)) 
       !Reduce remaining submatrix. 
    end do
  END SUBROUTINE ludcmp

SUBROUTINE balanc(a) 
  USE nrutil, ONLY : assert_eq 
  IMPLICIT NONE 
  CCTK_REAL, DIMENSION(:,:), INTENT(INOUT) :: a 
  CCTK_REAL, PARAMETER :: RADX=radix(a),SQRADX=RADX**2
  !Given a NxN matrix a this routine replaces it by a balanced matrix 
  !with identical eigenvalues. A symmetric matrix is already balanced and is 
  !unaffected by this procedure.
  !The parameter RADX is the machine's floating-point radix.

  !(C) 1986-1996 Numerical Recipes Software

  CCTK_INT :: i,last,ndum 
  CCTK_REAL :: c,f,g,r,s 
  ndum=assert_eq(size(a,1),size(a,2),'balanc')
  do 
     last=1 
     do i=1,size(a,1) !Calculate row and colum norms. 
        c=sum(abs(a(:,i)))-a(i,i) 
        r=sum(abs(a(i,:)))-a(i,i) 
        if (c /= 0.0 .and. r /= 0.0) then !If both are nonzero, 
           g=r/RADX 
           f=1.0 
           s=c+r 
           do  
              !find the integer power of the machine radix that comes closest 
              !to balancing the matrix. 
              if (c >= g) exit 
              f=f*RADX 
              c=c*SQRADX 
           end do
           g=r*RADX 
           do 
              if (c <= g) exit 
              f=f/RADX 
              c=c/SQRADX 
           end do
           if ((c+r)/f < 0.95d0*s) then 
              last=0 
              g=1.0d0/f 
              a(i,:)=a(i,:)*g !Apply similarity transformation. 
              a(:,i)=a(:,i)*f 
           end if
        end if
     end do
     if (last /= 0) exit 
  end do
END SUBROUTINE balanc

SUBROUTINE hqr(a,wr,wi) 
  USE nrutil, ONLY : assert_eq,diagadd,nrerror,upper_triangle 
  IMPLICIT NONE 
  CCTK_REAL, DIMENSION(:), INTENT(OUT) :: wr,wi 
  CCTK_REAL, DIMENSION(:,:), INTENT(INOUT) :: a
  !Finds all eigenvalues of an NxN upper Hessenberg matrix a. 
  !On input a can be exactly as output from elmhes chap.11.5;
  !on output it is destroyed.
  !The real and imaginary parts of the N eigenvalues are returned in 
  !wr and wi, respectively.
  CCTK_INT :: i,its,k,l,m,n,nn,mnnk 
  CCTK_REAL:: anorm,p,q,r,s,t,u,v,w,x,y,z 
  CCTK_REAL, DIMENSION(size(a,1)) :: pp 
  n=assert_eq(size(a,1),size(a,2),size(wr),size(wi),'hqr')
  anorm=sum(abs(a),mask=upper_triangle(n,n,extra=2)) 
  !Compute matrix norm for possible use in locating single small 
  !subdiagonal element. 
  nn=n
  t=0.0 !Gets changed only by an exceptional shift. 
  do !Begin search for next eigenvalue: "Do while nn >= 1".
     if (nn < 1) exit 
     its=0
     iterate: do !Begin iteration. 
        small: do l=nn,2,-1 !Look for single small subdiagonal element. 
           s=abs(a(l-1,l-1))+abs(a(l,l)) 
           if (s == 0.0) s=anorm 
           if (abs(a(l,l-1))+s == s) then 
              a(l,l-1)=0.0 
              exit small 
           end if
        end do small
        x=a(nn,nn) 
        if (l == nn) then !One root found. 
           wr(nn)=x+t 
           wi(nn)=0.0 
           nn=nn-1 
           exit iterate !Go back for next eigenvalue. 
        end if
        y=a(nn-1,nn-1) 
        w=a(nn,nn-1)*a(nn-1,nn) 
        if (l == nn-1) then !Two roots found ... 
           p=0.5d0*(y-x) 
           q=p**2+w 
           z=sqrt(abs(q)) 
           x=x+t 
           if (q >= 0.0) then !...a real pair ... 
              z=p+sign(z,p) 
              wr(nn)=x+z 
              wr(nn-1)=wr(nn) 
              if (z /= 0.0) wr(nn)=x-w/z 
              wi(nn)=0.0 
              wi(nn-1)=0.0 
           else !...a complex pair. 
              wr(nn)=x+p 
              wr(nn-1)=wr(nn) 
              wi(nn)=z 
              wi(nn-1)=-z 
           end if
           nn=nn-2 
           exit iterate !Go back for next eigenvalue. 
        end if
        !No roots found. Continue iteration. 
        if (its == 30) call nrerror('too many iterations in hqr') 
        if (its == 10 .or. its == 20) then !Form exceptional shift. 
           t=t+x 
           call diagadd(a(1:nn,1:nn),-x) 
           s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
           x=0.75d0*s 
           y=x 
           w=-0.4375d0*s**2 
        end if
        its=its+1 
        do m=nn-2,l,-1 
           !Form shift and then look for 2 consecutive small subdiagonal elements.
           z=a(m,m) 
           r=x-z 
           s=y-z
           p=(r*s-w)/a(m+1,m)+a(m,m+1) !Equation (11.6.23) of NR. 
           q=a(m+1,m+1)-z-r-s 
           r=a(m+2,m+1) 
           s=abs(p)+abs(q)+abs(r) !Scale to prevent overflow or underflow. 
           p=p/s 
           q=q/s 
           r=r/s 
           if (m == l) exit 
           u=abs(a(m,m-1))*(abs(q)+abs(r)) 
           v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1))) 
           if (u+v == v) exit !Equation (11.6.26) of NR. 
        end do
     do i=m+2,nn 
        a(i,i-2)=0.0 
        if (i /= m+2) a(i,i-3)=0.0 
     end do
     do k=m,nn-1 !Double QR step on rows l to nn and columns m to nn 
        if (k /= m) then 
           p=a(k,k-1) !Begin setup of Householder vector. 
           q=a(k+1,k-1) 
           r=0.0 
           if (k /= nn-1) r=a(k+2,k-1) 
           x=abs(p)+abs(q)+abs(r) 
           if (x /= 0.0) then 
              p=p/x !Scale to prevent overflow or underflow. 
              q=q/x 
              r=r/x 
           end if
        end if
        s=sign(sqrt(p**2+q**2+r**2),p) 
        if (s /= 0.0) then 
           if (k == m) then 
              if (l /= m) a(k,k-1)=-a(k,k-1) 
           else 
              a(k,k-1)=-s*x 
           end if
           p=p+s !Equations (11.6.24) of NR. 
           x=p/s 
           y=q/s 
           z=r/s 
           q=q/p 
           r=r/p !Ready for row modification. 
           pp(k:nn)=a(k,k:nn)+q*a(k+1,k:nn) 
           if (k /= nn-1) then 
              pp(k:nn)=pp(k:nn)+r*a(k+2,k:nn) 
              a(k+2,k:nn)=a(k+2,k:nn)-pp(k:nn)*z 
           end if
           a(k+1,k:nn)=a(k+1,k:nn)-pp(k:nn)*y 
           a(k,k:nn)=a(k,k:nn)-pp(k:nn)*x 
           mnnk=min(nn,k+3) !Column modification. 
           pp(l:mnnk)=x*a(l:mnnk,k)+y*a(l:mnnk,k+1) 
           if (k /= nn-1) then 
              pp(l:mnnk)=pp(l:mnnk)+z*a(l:mnnk,k+2) 
              a(l:mnnk,k+2)=a(l:mnnk,k+2)-pp(l:mnnk)*r 
           end if
           a(l:mnnk,k+1)=a(l:mnnk,k+1)-pp(l:mnnk)*q 
           a(l:mnnk,k)=a(l:mnnk,k)-pp(l:mnnk) 
        end if
     end do
  end do iterate  !Go back for next iteration on current eigenvalue.
end do
END SUBROUTINE hqr



end module nrmodule
