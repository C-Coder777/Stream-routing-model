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
!--Status = 0 for ending of burn-in, 1 for after burn-in, -1 for before burn-in
!----------------------------------------------------------------------------------
Subroutine DREAM_2(d, ncr, n, sig, limit, num, b1, b2, M0, TK, Psk, Out_step)
!--
Use Model_data1
!--
implicit none
integer d, ncr, n, num, limit, sig, TK, M0, im, Out_step
real*8 b1, b2, SR, SRd(d), Psk, SR_Con(limit+num, d), SR_Chain(limit+num, d, n), ZM(M0+(limit+num)*n/TK, d)
integer i, j, m, d2, sig2, c1, c2, r(sig*2), r1(sig), r2(sig), i1, j1, i2, j2, bnum, k1, k2
real*8 z(d), e(d), ibu(d), jumpsize, CR(ncr), crv, p(ncr), pnom(ncr-1)
real*8 sd(d), x(d), chains(num, d+1, n), Q_out(num, Qn, n), Q_out1(nn1), Q_pre(num, Qm, n), Q(Qn)
real*8 dx(d), u, lik, acc_pro, ss, s, ex, ar1, ar2, rr, log_lik
real*8 para(d), cov(d, d), mu(d), xnor(d, 1)
integer seed, numnor, numnom, ix(ncr), sn(3), label(d), r_label, Status, i0, num_ZM
real*8 line(2, d), points(2, d), p_points(2, d), acc_coeff, paras_burn(d+1, n), w
Character (len = 20) Head
!--
Status = 0
num_ZM = M0+(limit+num)*n/TK
!--
numnor=1
!----numnom for multinomial
numnom=1
!--
do i=1, ncr
  CR(i)=i*1.0/(ncr*1.0)
end do
sig2=2*sig
c1=1
i2=0
j2=0
!--
SR_Con=0.0
SR_Chain=0.0
chains=0.0
Q_out = 0.0
Q_pre = 0.0
!--
open (unit=101, file='Status.txt', recl=50, blank='null', position='rewind', pad='yes')
read (unit=101, fmt='(i1)') Status
read (unit=101, fmt='(i8)') c2
close (101)
!--
if (Status == 0) then
  call Read_files_1(d, n, ncr, limit, num, bnum, im, p, SR_Chain, SR_Con, ZM, num_ZM, paras_burn)
  i0 = 1
!--
  do i=1, n
    chains(1, 1: d+1, i) = paras_burn(1: d+1, i)
  end do
!--
else if (Status == 1) then
  call Read_files_2 (d, n, ncr, limit, num, c2, i2, bnum, im, p, SR_Chain, SR_Con, ZM, num_ZM, chains)
  i0 = c2 - bnum
!--
end if
!----cov, numnor, mu for generate multivariate normal number ibu
cov=0.0
do i=1, d
  cov(i, i)=(ranges(i, 2) - ranges(i, 1))**2/120.0
  mu(i)=0.5*(ranges(i, 1)+ranges(i, 2))
end do
!-----------------------------------------------------------------------
!--Generating the samples of posterior distribution after burn-in period
!-----------------------------------------------------------------------
!====
write(*,*) 'Start of the period after burn in'
do i=i0, num-1
!--
  do j=1, n
    d2=d
    label=0
    x(1:d)=chains(i, 1:d, j)
!------obtaining ibu, r1,r2,jumpsize,dx
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
!-------------------------------------compute d2
    crv=CR(m)
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
    i1=j+(i-1)*n
	j1=mod(i1,5)
	if(j1==0) then
	  jumpsize=1.0
	else
      jumpsize=2.38/sqrt(2.0*sig*d2*1.0)
	end if
!----
	dx=0.0
	do i1=1, d
	  do j1=1, sig
	    k1=r1(j1)
		k2=r2(j1)
	    dx(i1)=dx(i1)+ZM(k1, i1)-ZM(k2, i1)
	  end do
	end do
!----obtaining a alternative point   
    do i1=1,d
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
      k2=im-n+j
      call random_number (u)
      rr = int(u*1E8)*1.0      
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
!	  acc_pro=lik/chains(i, d+1, j)
	  lik=log_lik
	  acc_pro=exp(lik - chains(i, d+1, j)*1.5)
! 1.0 - 1.5
	else
!	  lik=0.0	  
	  lik=-1.0D+15
	  acc_pro=0.0
	end if
!===============
! Design for snooker update
    i1=Psk*100
    if (mod(i, i1)==0) then
      acc_pro=acc_pro* acc_coeff
    end if
!===============
    call random_number(u)
	if (acc_pro>=u) then
	  chains(i+1, 1:d, j)=z(1:d)
      chains(i+1, d+1, j)=lik
      Q_out(i+1, 1: Qn, j) = Q(1: Qn)
	  i2=i2+1
    else
	  chains(i+1, 1:d+1, j)=chains(i, 1:d+1, j)
      Q_out(i+1, 1: Qn, j) = Q_out(i, 1: Qn, j)
    end if

!----
  end do
  if (mod(i, 100) == 0) then
    write(*, fmt='(a39, i6)') 'After Burn-in Period: Iteration time =', i
  end if
!--Update archive ZM
  if (mod(i, TK)==0) then
    do i1=1, n
      ZM(im+i1, :)=chains(i, 1:d, i1)
    end do
    im=im+n
  end if
!--
!--GR record
  c2=c2+1
  do i1=1, n
  	 SR_Chain(c2, 1:d, i1) = chains(i+1, 1:d, i1)
  end do
  i1=limit+num
  call GR_test(SR_Chain, i1, c2, d, n, SRd)
  SR_Con(c2, 1:d)=SRd(1:d)
!---------------------------------------------------------------------------------------------------------------------------------
  if (i  >= Out_step) then
    if (Mod(i, Out_step) == 0 .or. i == num - 1) then
      write(*,*) 'Writing data into file'
!----
      open (unit=103, file='SR_Records.txt', recl=500, blank='null', position='rewind', pad='yes') 
      write (unit=103, fmt='(a30)') 'SR records of each interation'
      do i1=1, c2
          write (unit=103, fmt='(<d>f20.4)') SR_Con(i1, 1: d)
      end do
      close (103)
!--
      open (unit=104, file='Arichive_matrix.txt', recl=500, blank='null', position='rewind', pad='yes')
      write (unit=104, fmt='(2i8)') im, i2
      do i1 = 1, im
        write (unit=104, fmt='(<d>f14.6)') ZM(i1, 1: d)
      end do
      close (104)
!----
      open (unit=105, file='Status.txt', recl=50, blank='null', position='rewind', pad='yes')
      write (unit=105, fmt='(i1)') 1
      write (unit=105, fmt='(i8)') c2
      close (105)
!-----
      open (unit=107, file='chains.txt', recl=500, blank='null', position='rewind', pad='yes') 
      write (unit=107, fmt='(a33)') 'parameters and likelihood values'
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=107, fmt='(<d>f14.6, E16.6E3)') chains(j1, 1:d+1, i1)
        end do
      end do
      close (107)
!-----
      
      open (unit=108, file='Q_out1.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=108, fmt='(600f8.2)') Q_out(j1, 1: 600, i1)
        end do
      end do
      close (108)  
     
      open (unit=109, file='Q_out2.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=109, fmt='(600f8.2)') Q_out(j1, 601: 1200, i1)
        end do
      end do
      close (109)  
      
      open (unit=112, file='Q_out3.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=112, fmt='(600f8.2)') Q_out(j1, 1201: 1800, i1)
        end do
      end do
      close (112)
     
      open (unit=113, file='Q_out4.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=113, fmt='(600f8.2)') Q_out(j1, 1801: 2400, i1)
        end do
      end do
      close (113)
      
      open (unit=118, file='Q_out5.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=118, fmt='(600f8.2)') Q_out(j1, 2401: 3000, i1)
        end do
      end do
      close (118)
      
      open (unit=119, file='Q_out6.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=119, fmt='(600f8.2)') Q_out(j1, 3001: 3600, i1)
        end do
      end do
      close (119)
      
      open (unit=120, file='Q_out7.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=120, fmt='(600f8.2)') Q_out(j1, 3601: 4200, i1)
        end do
      end do
      close (120)
      
      open (unit=121, file='Q_out8.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=121, fmt='(600f8.2)') Q_out(j1, 4201: 4800, i1)
        end do
      end do
      close (121)
      
    open (unit=122, file='Q_out9.txt', recl=5000, blank='null', position='rewind', pad='yes') 
      do i1=1, n
        do j1=1, c2-bnum
          write (unit=122, fmt='(313f8.2)') Q_out(j1, 4801: 5113, i1)
        end do
      end do
      close (122)
!--
      
       
      
    end if
  end if
!---------------------------------------------------------------------------------------------------------------------
end do
ar2=(i2*1.0)/(num*n*1.0)
!==============
open (unit=110, file='DREAMzs_2.txt', recl=500, blank='null', position='rewind', pad='yes')
write(unit=110, fmt='(a28, f7.4)') 'Accepting rate after Burn_in', ar2
close (110)
!--
Call Prediction(chains, Q_pre, num, d, n, Qm)
!--
open (unit=111, file='Q_pre1.txt', recl=5000, blank='null', position='rewind', pad='yes') 
do i=1, n
  do j=1, num
    write (unit=111, fmt='(600f8.2)') Q_pre(j, 1: 600, i)
  end do
end do
close (111) 
!--
open (unit=114, file='Q_pre2.txt', recl=5000, blank='null', position='rewind', pad='yes') 
do i=1, n
  do j=1, num
    write (unit=114, fmt='(600f8.2)') Q_pre(j, 601: 1200, i)
  end do
end do
close (114) 

open (unit=115, file='Q_pre3.txt', recl=5000, blank='null', position='rewind', pad='yes') 
do i=1, n
  do j=1, num
    write (unit=115, fmt='(600f8.2)') Q_pre(j, 1201: 1800, i)
  end do
end do
close (115) 

open (unit=116, file='Q_pre4.txt', recl=5000, blank='null', position='rewind', pad='yes') 
do i=1, n
  do j=1, num
    write (unit=116, fmt='(600f8.2)') Q_pre(j, 1801: 2192, i)
  end do
end do
close (116) 



write(*,*) '>>>> End of DREAM-MCMC Program <<<<'
!--
!--
End subroutine
!============================================================================