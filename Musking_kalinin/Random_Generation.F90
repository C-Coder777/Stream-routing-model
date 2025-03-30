!===========================================================================
! These programs are used to generate some random numbers
!-----------------------------------------------------------------------------
!========Generate multivariate uniform random numbers e(d) from U(-b, b)
subroutine multi_uni(d, b, e)
implicit none
integer i, d
real*8 b, e(d), u
do i=1, d
  call random_number(u)
  e(i)=-b+2.0*b*u
end do
end subroutine
!-----------------------------------------------------------------------------
!========generate n different random integer numbers p(n) from a to b,and c is not include in p(n)
!--r is the seed of this program, r must be odd number, and p(n) is different for different r
subroutine ran_inter(r, a, b, c, n, p)
implicit none
integer n, a, b, c, k, l, i, kk
integer p(n), m, j
real*8 s, x, r
!--
i = int(r)
if (Mod(i, 2) == 0) then
  r = (i + 1)*1.0
end if
!--
kk=0
30 s=b-a+1.0
k=int(log(s-0.5)/log(2.0)+1)
l=1
do i=1,k
  l=2*l
end do
k=1
s=4.0*l
i=1
20 if ((i<=l).and.(k<=n)) then
 r=r+r+r+r+r
 m=int(r/s)
 r=r-m*s
 j=int(a+r/4.0)
 if (j<=b) then
   p(k)=j
   k=k+1
 end if
 i=i+1 
goto 20 
end if
x=1.0
do i=1,n
  x=x*(p(i)-c)
end do
if (x==0.0) then
  r=r/5.0+2.0
  kk=kk+1
  if(kk>=100) then
    write(*,*) 'Fatal error in regenerating p(n) in ran_inter subroutine'
    return
  end if
  goto 30
end if
end subroutine
!---------------------------------------------------------------------------------------------------------------
!========generate n different random integer numbers p(n) from a to b
!--r is the seed of this program, r must be odd number, and p(n) is different for different r
subroutine ran_inter2(r, a, b, n, p)
implicit none
integer n, a, b, k, l, i, kk
integer p(n),m,j
real*8 s,r
!--
i = int(r)
if (Mod(i, 2) == 0) then
  r = (i + 1)*1.0
end if
!--
kk=0
s=b-a+1.0
k=int(log(s-0.5)/log(2.0)+1)
l=1
do i=1,k
  l=2*l
end do
k=1
s=4.0*l
i=1
20 if ((i<=l).and.(k<=n)) then
 r=r+r+r+r+r
 m=int(r/s)
 r=r-m*s
 j=int(a+r/4.0)
 if (j<=b) then
   p(k)=j
   k=k+1
 end if
 i=i+1 
goto 20 
end if
end subroutine
!---------------------------------------------------------------------------------------------------------------