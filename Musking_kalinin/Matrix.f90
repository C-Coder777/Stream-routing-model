!-------------------------------------------------------------------------------------------------
Subroutine Matrix_Det(c, n, det)
implicit none
integer i, j, k, is, js, n
real*8 a(n,n), det, f, q, s, c(n, n)
a=c
f=1.0
det=1.0
do k=1, n-1
  q=0.0
  do i=k, n
	  do j=k, n
	    if (abs(a(i, j)) > q) then
	      q=abs(a(i, j))
	      is=i
	      js=j
	  end if
	end do
  end do
!--
  if (q+1.0==1.0) then
	det=0.0
	return
  end if
!--
  if (is/=k) then
	f=-f
	do j=k, n
	    s=a(k, j)
	    a(k, j)=a(is, j)
	    a(is, j)=s
    end do
  end if
!--
  if (js/=k) then
	f=-f
	do i=k, n
	  s=a(i, js)
	  a(i, js)=a(i, k)
	  a(i, k)=s
    end do
  end if
!--
  det=det*a(k, k)
  do i=k+1, n
    s=a(i, k)/a(k, k)
    do j=k+1,n
      a(i, j)=a(i, j)-s*a(k, j)
    end do
  end do
end do
det=f*det*a(n, n)
End subroutine
!============================================================
Subroutine inverse(aa, c, n)
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n - dimension
! output ...
! c(n,n) - inverse matrix of A
!===========================================================
implicit none 
integer n
double precision a(n, n), c(n, n), aa(n, n)
double precision L(n, n), U(n, n), b(n), d(n), x(n)
double precision coeff
integer i, j, k
a=aa
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
End subroutine
!----------------------------------------------------