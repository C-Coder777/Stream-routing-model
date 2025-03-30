!====
! Compute the Likelihood value
!--
Subroutine Likelihood(obs, sim, n, log_lik)
implicit none
integer i, n
real*8 obs(n), sim(n), det, log_lik
real*8 c1(n, 1), c2(1, n), c3(1, n), cov(n, n), cov_1(n, n), c(1, 1)
!--
cov=0.0
cov_1=0.0
!--
do i=1, n
 cov(i, i) = obs(i)*0.1*obs(i)*0.1
!  cov(i, i) = 0.2*0.2
!  cov(i, i) = 0.1*0.1
  !cov(i, i) = obs(i)*0.01
  cov_1(i, i) = 1.0/cov(i, i)
end do
!--
do i=1, n
  c1(i, 1)=obs(i) - sim(i)
end do
!--
c2=transpose(c1)
c3=matmul(c2, cov_1)
c=matmul(c3, c1)
log_lik = - 0.5*c(1,1)
!--
End subroutine