!-----------------------------------------------------------------------------------------------
!----Code for a Markov Chain Monte Carlo simulation method----DREAMzs
!----A self-adaptive Differential Evolution learning strategy within a population-based evolutionary framework.
!----DiffeRential Evolution Adaptive Metropolis (DREAM) combined with snooker updating and sampling from past states.
!
!-----------------------------------------------------------------------------------------------
!  Parameter discription
!--d dimension of parameter set
!--n number of initial Markov Chains (parallel chain)
!--limit is the maximum length of burn-in period
!--nrc user defined parameter for deriving CR probability distribution
!--num number of iternations after burn-in period
!--sig number of chain pairs used for sample generation
!--r1(del),r2(del) are chain number, and they are different
!--b1,b2 are the b and b* in equation 4, Ud(-b,b), Nd(0,b*)
!----the use of b2 is different from Vrugt at here
!--e(d), ibu(d) sampled from Ud(-b,b), Nd(0,b*)
!--ilen is the initial length of burn in chain
!--st the start time to test outlier chains and convergence
!--ipoints(ilen*n, d), z(d) initial point and alternative point
!--p(ncr) probability distribtion of CR values, CR = (m/ncr, m=1,2,...,ncr)
!--CR(ncr) CR values
!--del(ncr) jumping distance for CR(ncr)
!--sd(d) standard devivation for each dimension of current x(n,d)
!--bnum length of burn-in period
!--ar1, ar2 acceptance rate durning and after burn-in
!--CRUpdate is the start time to update CR distribution
!--Status = 0 for ending of burn-in, 1 for after burn-in
!----------------------------------------------------------------------------------
Subroutine DREAM_1(d, ncr, n, sig, limit, st, ilen, CRUpdate, b1, b2, SR, M0, TK, Psk)
!--
Use Model_data1
!--
implicit none
integer d, ncr, n, sig, limit, st, ilen, CRUpdate, TK, M0, im
real*8 b1, b2, SR, SRd(d), Psk, SR_Con(limit, d), SR_Chain(limit, d, n), ZM(M0+limit*n/TK, d)
integer i, j, Lm(ncr), m, d2, sig2, c1, c2, r(sig*2), r1(sig), r2(sig), i1, j1, i2, j2, bnum, k1, k2
real*8 z(d), e(d), ibu(d), jumpsize, CR(ncr), crv, p(ncr), del(ncr), delta, pnom(ncr-1)
real*8 sd(d), x(d), limburn(limit, d+1, n), Q_out1(nn1),Q(Qn)
real*8 dx(d), u, lik, acc_pro, ss, s, ex, ar1, ar2, rr, log_lik
real*8 para(d), ipoints(ilen*n, d), ilik(ilen*n), cov(d, d), mu(d), xnor(d, 1)
integer seed, numnor, numnom, ix(ncr), sn(3), label(d), r_label, Status
real*8 line(2, d), points(2, d), p_points(2, d), acc_coeff
!----
call random_seed()
!--
numnor=1
!----numnom for multinomial
numnom=1
!--Burn-in period--
limburn=0.0
SR_Con=0.0
SR_Chain=0.0
!--
Lm=0
do i=1, ncr
  CR(i)=i*1.0/(ncr*1.0)
  p(i)=1.0/(ncr*1.0)
end do
del=0.0
sig2=2*sig
c1=1
c2=ilen-1
i2=0
j2=0
im=M0
!----cov, numnor, mu for generate multivariate normal number ibu
cov=0.0
do i=1, d
  cov(i, i)=(ranges(i, 2) - ranges(i, 1))**2/120.0
  mu(i)=0.5*(ranges(i, 1)+ranges(i, 2))
end do
!--
!---------------------------------------------
!-(1). obtain initial points ipoints(ilen*n, d) for archive ZM
i=ilen*n
call Pre_samples(ipoints, i, mu, cov, d)
!-(2). obtaining the likelihood values (probability densities) of these vectors
ilik=0.0
!--
do j=1, i
  para(1:d)=ipoints(j, 1:d)
  call Run_Model(para, d, log_lik, Q_out1, Q)
!  ilik(j)=exp(log_lik)
   ilik(j)=log_lik
end do
!-(3). initial points for ZM and burn in chain
do j=1, M0
  ZM(j, :)=ipoints(j, :)
end do
!--
do i=1, ilen
  do j=1, n
    i1=j+(i-1)*n
    limburn(i, 1:d, j)=ipoints(i1, 1:d)
    SR_Chain(i, 1: d, j)=limburn(i, 1: d, j)
    limburn(i, d+1, j)=ilik(i1)
  end do
end do
!----
!========
write(*,*) 'Start of burn in period'
!--the start of burn in period
!--
do i=ilen, limit-1
!--
  do j=1, n
    d2=d
    label=0
    x(1:d)=limburn(i, 1:d, j)
!------obtaining ibu, r1, r2, jumpsize, dx
    call random_number (u)
    rr = int(u*1E8)*1.0
    call random_number (u)
    seed = int(u*1E8)
    call multi_uni(d, b1, e)
	call multinormal_sample (d, numnor, cov, mu, seed, xnor)
	ibu(1:d)=xnor(1:d, 1)*b2
    call ran_inter2(rr, c1, im, sig2, r)
	do i1=1, sig
	  j1=sig+i1
	  r1(i1)=r(i1)
	  r2(i1)=r(j1)
	end do
!----------------------------------------compute d2
!----obtaining m by a multinomial distribution
    do i1=1, ncr-1
      pnom(i1)=p(i1)
    end do
  !--
    call genmul (numnom, pnom, ncr, ix)
    do i1=1, ncr
      if (ix(i1)==1) then
        m=i1
        exit
      end if
    end do
    crv=CR(m)
    Lm(m)=Lm(m)+1
	do i1=1, d
	  call random_number(u)
	  if(u<=(1.0-crv)) then
        label(i1)=1
		d2=d2-1
	  end if
	end do
	if (d2==0) then
	  d2=1
	end if
!----
    i1=j+(i-ilen)*n
	j1=mod(i1, 5)
	if(j1==0) then
	  jumpsize=1.0
	else
      jumpsize=2.38/sqrt(2.0*sig*d2*1.0)
	end if
!--
	dx=0.0
	do i1=1, d
	  do j1=1, sig
	    k1=r1(j1)
		k2=r2(j1)
	    dx(i1)=dx(i1)+ZM(k1, i1)-ZM(k2, i1)
	  end do
	end do
!----obtaining a alternative point   
    do i1=1, d
	  z(i1)=x(i1)+(1.0+e(i1))*jumpsize*dx(i1)+ibu(i1)
    end do
!----replace z by x with probability crv
    do i1=1, d
      if (label(i1) == 1) then
        z(i1)=x(i1)
      end if
    end do
!===========================
!  Snooker update
!----
  i1=Psk*100
  k1=3
  if (mod(i, i1)==0) then
!----
    if (j2==0) then
      call random_number (u)
      rr = int(u*1E8)*1.0
      call ran_inter2(rr, c1, im, k1, sn)
    else if (j2==1) then 
      call random_number (u)
      rr = int(u*1E8)*1.0
      k2=im-n+j
      call ran_inter(rr, c1, im, k2, k1, sn)
    end if
!----
    line(1, 1:d)=x(1:d)
    line(2, 1:d)=ZM(sn(1), 1:d)
    points(1, 1:d)=ZM(sn(2), 1:d)
    points(2, 1:d)=ZM(sn(3), 1:d)
    k1=2
    k2=2
    call Ortho_Pro(line, points, d, k1, p_points, k2)
!----
    do i1=1, d
      z(i1)=x(i1)+1.6829*(p_points(1, i1)-p_points(2, i1))
    end do
!----
    s=0.0
    ss=0.0
    do i1=1, d
      s=s+(z(i1)-line(2, i1))*(z(i1)-line(2, i1))
      ss=ss+(x(i1)-line(2, i1))*(x(i1)-line(2, i1))
    end do
    acc_coeff=(sqrt(s/ss))**(d-1)
  end if
!--
!  End of snooker update
!===========================
!----calculating the likelihood value of zi and the acceptance probability
    call range_jud(ranges, d, col, z, r_label)
    if (r_label == 1) then
	  call Run_Model(z, d, log_lik, Q_out1, Q)
!	  lik=exp(log_lik)
!	  acc_pro=lik/limburn(i, d+1, j)
	  lik=log_lik
	  acc_pro=exp(lik-limburn(i, d+1, j)*1.9)
! 1.1, 1.2, 1.25, 1.3, 1.4, 1.5
	else
!	  lik=0.0
	  lik=-1.0D+15
	  acc_pro=0.0
	end if
!===============
! Design for snooker update
  i1=Psk*100
  if (mod(i, i1)==0) then
    acc_pro=acc_pro*acc_coeff
  end if
!===============
    call random_number(u)
	if(acc_pro>=u) then
	  limburn(i+1, 1:d, j)=z(1:d)
      limburn(i+1, d+1, j)=lik
	  i2=i2+1
    else
	  limburn(i+1, 1:d+1, j)=limburn(i, 1:d+1, j)
    end if
!----calculating squared normalized jumping distance
    delta=0.0
	do i1=1, d
	  s=sum(limburn(1:i+1, i1, j))
	  ex=s/((i+1)*1.0)
	  ss=0.0
	  do j1=1, i+1
	    ss=ss+(limburn(j1, i1, j)-ex)*(limburn(j1, i1, j)-ex)
	  end do
	  ss=ss/((i+1)*1.0)
	  s=(limburn(i+1, i1, j)-limburn(i, i1, j))*(limburn(i+1, i1, j)-limburn(i, i1, j))/ss
	  delta=delta+s
    end do
    del(m)=del(m)+delta
!====(do j=1, n -->> end do)
  end do
!================
!--Update archive ZM
  j2=0
  if (mod(i, TK)==0) then
    do i1=1, n
      ZM(im+i1, :)=limburn(i, 1:d, i1)
    end do
    im=im+n
    j2=1
  end if
!================
!----updating of CR probability distribution p(ncr)
  if (minval(del(1:ncr)) >1.0D-06 .and. minval(Lm(1:ncr)) >1.0D-06) then
    s=sum(del(1:ncr)/Lm(1:ncr))
    if (i>=CRUpdate) then
      do j=1, ncr-1
        p(j)=(del(j)/Lm(j))/s
      end do
      p(ncr)=1.0D+00-sum(p(1:ncr-1))
    end if
  end if
!----
  if (mod(i, 100) == 0) then
!     write(*, fmt='(a40, i6)') 'During Burn-in Period: Iteration time =', i
   write(*, *) 'During Burn-in Period: Iteration time =', i
  end if
  c2=c2+1
!--Convergence test
  do i1=1, n
  	SR_Chain(c2, 1:d, i1) = limburn(i+1, 1:d, i1)
  end do
  i1=limit
  if (i2>=n*5) then
   call GR_test(SR_Chain, i1, c2, d, n, SRd)
   SR_Con(c2, 1:d)=SRd(1:d)
  end if
!----
  if(maxval(SRd(1:d))<=SR .and. i>=st) then
	bnum=i
	write(*,*) 'convergence is obtained!'
	write(*,*) 'Length of burn-in period is', i
	exit
  end if
!----
  if(i==limit-1) then
	bnum=limit
	write(*,*) 'convergence is not obtained !!!'
	write(*,*) 'Length of burn-in period is limit time'
	pause
  end if
!--
end do
ar1=(i2*1.0)/((c2-ilen)*n*1.0)
write(*,*) '>>>> End of burn-in period <<<<'
!-----------------------------------------------------------------------
!--
write(*,*) 'Writing data into file'
!--
open (unit=101, file='DREAMzs_1.txt', recl=500, blank='null', position='rewind', pad='yes')
write(unit=101, fmt='(a20, i6)') 'Length of burn-in:', bnum
write(unit=101, fmt='(a25, f8.4)') 'Accepting rate in Burn_in', ar1
write(unit=101, fmt='(a8, <ncr>f10.6)') 'p(ncr)=', p(1:ncr)
close (101)
!--
open (unit=102, file='Burn_in.txt', recl=500, blank='null', position='rewind', pad='yes') 
write(unit=102, fmt='(a33)') 'parameters and likelihood values'
do i=1, n
  do j=1, bnum
    write(unit=102, fmt='(<d>f14.6, E16.6E3)') limburn(j, 1:d+1, i)
  end do
end do
close(102)
!--
open (unit=103,file='SR_Records.txt', recl=500, blank='null', position='rewind', pad='yes') 
write(unit=103, fmt='(a30)') 'SR records of each interation'
do i=1, c2
    write(unit=103, fmt='(<d>f20.4)') SR_Con(i, 1: d)
end do
close(103)
!--
open (unit=104, file='Arichive_matrix.txt', recl=500, blank='null', position='rewind', pad='yes')
write(unit=104, fmt='(i8)') im
do i = 1, im
  write(unit=104, fmt='(<d>f14.6)') ZM(i, 1: d)
end do
close (104)
!--
open (unit=105, file='Status.txt', recl=50, blank='null', position='rewind', pad='yes')
write(unit=105, fmt='(i1)') 0
write(unit=105, fmt='(i8)') c2
close (105)
!--
open (unit=106, file='Burn_ends.txt', recl=500, blank='null', position='rewind', pad='yes')
do i=1, n
  write(unit=106, fmt='(<d>f14.6, E16.6E3)') limburn(bnum, 1: d+1, i)
end do
close(106)
!--
End subroutine
!============================================================================