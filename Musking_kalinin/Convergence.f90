!===================================================
!------------------convergence test; SR is the scale reduction score by Gelman and Rubin
!--chains(num, d, n) -> chains2(num/2, d, n2), l3=n.
subroutine GR_test(chains, number, num, d, n, SRd)
implicit none
integer i, j, k1, k2, num, d, n, number
real*8 chains(number, d, n), SRd(d), array(num/2)
real*8 a1, a2, mean(n, d), var(n, d), a3(n), b(d), w(d), g, a4
!--
k1=num/2
k2=num-num/2+1
!--
do i=1, d
  do j=1, n
    array(1:k1)=chains(k2:num, i, j)
    call mean_var(array, k1, a1, a2)
    mean(j, i)=a1
	var(j, i)=a2
  end do
end do
!--------------
do i=1, d
  a3(1:n)=mean(1:n, i)
  call mean_var(a3, n, a1, a2)
  b(i)=a2
end do
!----
do i=1, d
  a3(1:n)=var(1:n, i)
  call mean_var(a3, n, a1, a2)
  w(i)=a1
end do
!--------------
g=k1
do i=1, d
  a4=(g-1.0)/(g*1.0)+((n+1.0)/(n*1.0*g*1.0))*(b(i)/w(i))
  SRd(i)=sqrt(a4)
end do
!------
end subroutine
!--------------------------------------------------Get mean and variance
subroutine mean_var(array, number, ex, ss)
implicit none
integer number,i
real*8 array(number),ex,ss,sum1,sum2
ex=0.0
ss=0.0
sum1=0.0
sum2=0.0
do i=1,number
  sum1=sum1+array(i)
end do
ex=sum1/(number*1.0)
do i=1,number
  sum2=sum2+(array(i)-ex)*(array(i)-ex)
end do
ss=sum2/(number*1.0-1.0)
end subroutine