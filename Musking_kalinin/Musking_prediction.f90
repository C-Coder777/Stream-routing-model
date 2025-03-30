
subroutine msk_pre(k,x,n,dt,nt,qi,qo,qc)
implicit none

real*8 k
real*8 x
real dt

integer n
integer nt


real*8 qi(nt), qo(nt), qc(nt)
real*8 qo1, qo2, qi1, qi2

real*8 x1, c0, c1, c2

integer i, j

x1=k-k*x+0.5*dt
c0=(0.5*dt-k*x)/x1
c1=(k*x+0.5*dt)/x1
c2=(k-k*x-0.5*dt)/x1

do  i=1,n
  if (qc(i).le.0.0) qc(i)=qi(1)
end do

if (qo(1).le.0.0) qo(1)=qc(n)

do  i=2,nt
  do j=1,n
	qo1=qc(j)
	if (j.eq.1) then
	  qi1=qi(i-1)
	  qi2=qi(i)
	end if
	qo2=c0*qi2+c1*qi1+c2*qo1
	qi1=qo1
	qi2=qo2
	qc(j)=qo2
  end do
  qo(i)=qo2
end do
return
end subroutine
    
subroutine nash_bias_pre(nstep, qobs,qcal,nash,bias)
implicit none

integer nstep
real*8 qobs(nstep)
real*8 qcal(nstep)
real*8 nash
real*8 bias

integer i
real*8  s_qobs
real*8  s_qcal
real*8  qobs_mean
real*8  qcal_mean
real*8  s1
real*8  s2

s_qobs=0.
s_qcal=0.
s1=0.
s2=0.

do i=1,nstep
  s_qobs=s_qobs+qobs(i)
  s_qcal=s_qcal+qcal(i)
end do

qobs_mean=s_qobs/real(nstep)


do i=1,nstep
  s1=s1+(qobs(i)-qcal(i))**2
  s2=s2+(qobs(i)-qobs_mean)**2
end do      

bias=(s_qcal-s_qobs)/abs(s_qobs)
nash=1-s1/s2

return

end subroutine

!-------------------------------------------------------------------------------------------------Snow Runoff Model
! n is number of days
! m is number of areas
! nn is number of months
! unit transformation: 10000/86400 = 0.11574074
!--
Subroutine Run_model_pre(point, d, Q)
Use Model_data5
Use Model_data6
!--
implicit none
! this program routes the daily surface runoff and base flow that are calculated
! by vic from each grid cells to the watershed outlet through the lineal 
! reservior and muskingum methods.
integer, parameter:: ngrid=25 ! Number of grid cells within the watershed
real, parameter:: dt=24 ! Time step of model simulation
integer, parameter:: startyear=2016 ! Year starting the routing
integer, parameter:: endyear=2021 ! Year ending the routing
character(72), parameter:: gridinfo_dir='./PARAMETER/zones.txt' ! Directory of grid information file
character(72), parameter:: output_daily_dir='./results_CJ/routing_results_day(40)_pre.txt' ! Directory of output daily streamflow file
character(72), parameter:: output_monthly_dir='./results_CJ/routing_results_month(40)_pre.txt' ! Directory of output daily streamflow file
character(72), parameter:: flux_dir='./flux2012-2017/fluxes_' ! Directory of flux files
character(72), parameter:: qobs_dir='./PARAMETER/Qo_pre.txt' ! Directory of observed streamflow data file                
!--
integer  nt, nmus  ! Total number of time steps for model simulation
real*8  s, qrs0, qrg0, qri0
real*8  qgrid(n), qi(n), qj(n), qo(n), qc(n), pmean(n), emean(n), r1mean(n), r2mean(n), qobs(n)
integer  i, j, k, d, year(n), month(n), day(n)
real*8  temp, gridid, area, lat, lon, riverlen
real*8  p, e, rs, rg, rss, rgg, rii, qrs, qrg, qri, u
real*8  pmonth(n), emonth(n), rmonth(n), qobsmonth(n), qcalmonth(n) ! Monthly calculated streamflow
real*8  nash_day, nash_month, bias ! Relative runoff error  
real*8  point(d), Qave(nn), Q(n), Qb(n), Q_out(nn)
real*8  kmus, xmus, cs,sm, ci, cg, kg, ki
character:: chrtemp
character*72 filename
integer nmyfile


!--
!real,parameter:: sm=1.5
!real,parameter:: ki=0.59
!real,parameter:: kg=0.01
!real,parameter:: cs=0.59
!real,parameter:: ci=0.9
!real,parameter:: cg=0.995
real,parameter:: v=6
!real,parameter:: kmus=18
!real,parameter:: xmus=0.3
real, parameter:: qrs00=0. ! Initial surface discharge at each grid cell
real, parameter:: qri00=0. ! Initial interflow discharge at each grid cell
real, parameter:: qrg00=0. ! Initial groundwater discharge at each grid cell
real, parameter:: s0=0.0 ! Initial free water storage at each grid cell

kmus = point(1)
xmus = point(2)
sm=point(3)
cs=point(4)
cg = point(5)
ci = point(6)
kg = point(7)
!v = point(8)
ki = point(8)

Q = 0.0
Qb = 0.0
Qave=0.0
Q_out=0.0
!write (*,'(8F10.4)')kmus,xmus 

! Calculate the total number of time steps for routing

nt=0
do i=startyear, endyear
  if (mod(i,4).eq.0.and.mod(i,100).ne.0 &
	  .or. mod(i,100).eq.0.and.mod(i,400).eq.0) then
	nt=nt+366
  else
	nt=nt+365
  end if
end do

open (unit=6,file=gridinfo_dir, status='old')
open (unit=7,file=output_daily_dir,action='write') 
open (unit=9,file=qobs_dir,action='read')
open (unit=10,file=output_monthly_dir,action='write') 

do i=1,nt
  read(9,*) qobs(i)
end do 
close (9)
	
do i=1,nt
  qj(i)=0.
  pmean(i)=0.
  emean(i)=0.
  r1mean(i)=0.
  r2mean(i)=0.
end do

read(6,*) chrtemp
!write(*,*) ngrid
!pause
do i=1,ngrid
  
  do j=1,nt
	qgrid(j)=0
	qi(j)=0.
	qo(j)=0.
	qc(j)=0.
  end do

  qrs0=qrs00
  qrg0=qrg00
  qri0=qri00
  s=s0

  read (6,*) gridid,lat,lon,area,riverlen
  u=area/(dt*3.6) 
filename=flux_dir 
nmyfile=len_trim(filename)
write (filename(nmyfile+1:nmyfile+6),'(f6.3)') lat
  write (filename(nmyfile+7:nmyfile+7),'(a1)') '_'
write (filename(nmyfile+8:nmyfile+15),'(f6.3)') lon
!write (*,*) filename
!   pause
open (unit=8, file=filename, status='old')
  nmus=nint(riverlen/3600.0/dt/v)
  if (nmus.le.1) nmus=1

  read (8,*) year(1), month(1), day(1), p, e, &
		  rs, rg, temp, temp, temp, temp, temp, &
		  temp, temp, temp, temp, temp, temp, &
		  temp, temp, temp, temp, temp, temp, &
		  temp, temp, temp
  do while (year(1).lt.startyear)
	read (8,*) year(1), month(1), day(1), p, e, &
		  rs, rg, temp, temp, temp, temp, temp, &
		  temp, temp, temp, temp, temp, temp, &
		  temp, temp, temp, temp, temp, temp, &
		  temp, temp, temp
  end do
  s=s+rs+rg
	
  if (s.gt.sm) then
	  rss=s-sm
	  rgg=sm*kg
	  rii=sm*ki
	  s=s-rss-rgg-rii
  else
	  rss=0
	  rgg=s*kg
	  rii=s*ki
	  s=s-rss-rgg-rii
  end if            

  qrs=qrs0*cs+rss*(1-cs)*u
  qrg=qrg0*cg+rgg*(1-cg)*u
  qri=qri0*ci+rii*(1-ci)*u
  qgrid(1)=qrs+qrg+qri
  qrs0=qrs
  qrg0=qrg
  qri0=qri

  r1mean(1)=r1mean(1)+rs/real(ngrid)
  r2mean(1)=r2mean(1)+rg/real(ngrid)
  pmean(1)=pmean(1)+p/real(ngrid)
  emean(1)=emean(1)+e/real(ngrid)

  do j=2,nt
! read flux file
  read (8,*) year(j), month(j), day(j), p, e, &
		  rs, rg, temp, temp, temp, temp, temp, &
		  temp, temp, temp, temp, temp, temp, &
		  temp, temp, temp

! routing within grid cells
	s=s+rs+rg
	
	if (s.gt.sm) then
	  rss=s-sm
	  rgg=sm*kg
	  rii=sm*ki
	  s=s-rss-rgg-rii
	else
	  rss=0
	  rgg=s*kg
	  rii=s*ki
	  s=s-rss-rgg-rii
	end if            

	qrs=qrs0*cs+rss*(1-cs)*u
	qrg=qrg0*cg+rgg*(1-cg)*u
	qri=qri0*ci+rii*(1-ci)*u
	qgrid(j)=qrs+qrg+qri
	qrs0=qrs
	qrg0=qrg
	qri0=qri

	r1mean(j)=r1mean(j)+rs/real(ngrid)
	r2mean(j)=r2mean(j)+rg/real(ngrid)
	pmean(j)=pmean(j)+p/real(ngrid)
	emean(j)=emean(j)+e/real(ngrid)
	
  end do
  close (8)

! river flow routing

  do j=1,nt
	qi(j)=qgrid(j)
  end do
  
  do j=1,nmus
	qc(j)=qi(1)
  end do
  !write(*,*) kmus,xmus,nmus,dt,nt,qi,qo,qc
  !pause
  call msk_pre(kmus,xmus,nmus,dt,nt,qi,qo,qc)
  !write(*,*) xmus
  !pause

  
  do j=1,nt
	qj(j)=qj(j)+qo(j)+Qb(j)
  end do
end do 
  open (unit=26, file='baseflow_pre.txt', recl=100, blank='null', position='rewind', pad='yes')
  do j = 1, n
    read(unit=26, fmt='(f7.2)') Qb(j)
  end do
  close(26)
    do j=1,nt
	qj(j)=qj(j)+Qb(j)
  end do
  
do i=1,nt
  pmonth(i)=0.
  emonth(i)=0.
  rmonth(i)=0.
  qobsmonth(i)=0.
  qcalmonth(i)=0.
end do

pmonth(1)=pmean(1)
emonth(1)=emean(1)
rmonth(1)=r1mean(1)+r2mean(1)
qobsmonth(1)=qobs(1)
qcalmonth(1)=qj(1)
write (10,'(7a10)') 'Year','Month','Prec','Evap','Runoff','Q_obs','Q_cal'
j=1
do i=2,nt
  if (month(i).eq.month(i-1)) then
	pmonth(j)=pmonth(j)+pmean(i)
	emonth(j)=emonth(j)+emean(i)
	rmonth(j)=rmonth(j)+r1mean(i)+r2mean(i)
	qobsmonth(j)=qobsmonth(j)+qobs(i)
	qcalmonth(j)=qcalmonth(j)+qj(i)
  else
	qobsmonth(j)=qobsmonth(j)/real(day(i-1))
	qcalmonth(j)=qcalmonth(j)/real(day(i-1))
	write (10,'(2i10,5f10.3)')year(i-1),month(i-1),&
        pmonth(j),emonth(j),rmonth(j),qobsmonth(j),qcalmonth(j) 
	j=j+1
	pmonth(j)=pmonth(j)+pmean(i)
	emonth(j)=emonth(j)+emean(i)
	rmonth(j)=rmonth(j)+r1mean(i)+r2mean(i)
	qobsmonth(j)=qobsmonth(j)+qobs(i)
	qcalmonth(j)=qcalmonth(j)+qj(i)            
  end if
  if (i.eq.nt) then
	qobsmonth(j)=qobsmonth(j)/real(day(i))
	qcalmonth(j)=qcalmonth(j)/real(day(i))
  end if
  
end do
!write (10,'(2i10,5f10.3)')year(i-1),month(i-1),&
!	  pmonth(j),emonth(j),rmonth(j),qobsmonth(j),qcalmonth(j) 

write (7,'(8a10)') 'Year','Month','Day','Prec','Evap','Runoff','Q_obs','Q_cal'
do i=1,nt
  write (7,'(3i10,5f10.3)')year(i),month(i),&
   day(i), pmean(i),emean(i),r1mean(i)+r2mean(i),qobs(i), &
   qj(i)
end do 




 !Calculate Nash and Bias for daily streamflow simulation
open(unit=101,file='NS_day_pre.txt')
open(unit=102,file='NS_month_pre.txt')
open(unit=103,file='bias_pre.txt')
call nash_bias_pre(nt,qobs,qj,nash_day,bias)

 !Calculate Nash and Bias for monthly streamflow simulation
call nash_bias_pre(j,qobsmonth,qcalmonth,nash_month,bias)
  write (101,'(a110,F10.3)') 'Nash-Sutcliffe model efficiency coefficient for daily streamflow simulation:  ', Nash_day
  write (102,'(a110,F10.3)') 'Nash-Sutcliffe model efficiency coefficient for monthly streamflow simulation:  ', Nash_month
  write (103,'(a110,F10.3)') 'Relative error between Qcal and Qobs:  ', bias        
Q(1:n)=qj(1:nt)
Q_out(1: nn)=qcalmonth(1: nn)
Qave(1: nn) = qcalmonth(1: nn)

close (6)
close (7)
close (10)

return
End