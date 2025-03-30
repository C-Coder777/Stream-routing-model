!====
! Module used for piror information
!--
Module Model_data1
implicit none
integer, parameter :: nn1=168
integer, parameter :: Qn=5113
integer, parameter :: Qm=2192
integer, parameter :: d1=8
integer, parameter :: col=2
real*8 ranges (d1, col)
!--
! x, y, cr, cs, a
!--
data ranges(1, 1:2) /16, 30/
data ranges(2, 1:2) /0.0, 1.0/
data ranges(3, 1:2) /60.00, 100/
data ranges(4, 1:2) /0.1, 1.0/
data ranges(5, 1:2) /0.5, 1.0/
data ranges(6, 1:2) /0.8, 1.0/
data ranges(7, 1:2) /0.05, 0.65/
!data ranges(8, 1:2) /2.0, 20.0/
data ranges(8, 1:2) /0.05, 0.65/
!--
End Module
!====
! Module for Head observations
!--
Module Model_data2
implicit none
!--
integer, parameter :: n3 = 168
real*8 Qo_ave(n3)
integer nd(n3)
!--
data Qo_ave(1: n3) /5.36,5.57,4.31,4.12,14.3,46.7,92.6,108,40.6	,16,8.54,6.43,&
                    5.75,6.11,6.1,7.99,19.8,67.4,131,93.8,61.1,15.4,9.34,7.09,&
                    6.16,6.22, 6.96,9.29,30.8,57.6,92.8,112,47.1,11.4,6.85,5.61,&
                    5.94,5.5,4.89,11.6,33.6,130,126,105,56.6,18.4,9.63,7.59,&
                    6.78,8.33,6.07,13.5,41.8,69.5,133,103,36.9,13,9.42,7.37,&
                    7.18	,6.12,5.1,6.57,24.8,68.1,76.3,76.5,47.2,10.3,7.21,6.03,&
                    5.36, 2.79,1.5,0.584,19.5,84,82.7,56.1,13.4, 1.74,1.8,2.06,&
                    5.81,4.37,3.29,3.1,12.4,54.2,74.5,88.4,24.8,8.99,5.58,5.09,&
                    5.74,5.2,3.37,6.17,36.5,79.6,142,123,45.7,23.8,8.96,5.68,&
                    5.31,5.14	,4.87,10.7,21.7,49.8,84.5,99.4,69.9,12.3,7.33,5.45,&
                    4.7, 4.19, 5.01, 12.7, 34.8, 97.5, 104, 127, 59, 14.2, 7.96, 6.25,&
                    4.69, 4.69, 2.41, 1.94, 35.6, 115, 135, 136, 48.6, 6.84, 2.73, 2.63,&
                   5.23, 4.94, 3.65, 3.3, 8.55, 46.6, 83.3, 83, 35.2, 4.78, 7.9, 5.77,&
                    4.93, 4.92, 3.65, 5.9, 14.1, 51.4, 95.5, 131, 41.7, 13.9, 11, 7.16	/

data nd(1: n3) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,&
                31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!--
End Module
!====
! n1 n2 for 2002, n3, n4 for 2003
!--
Module Model_data3
implicit none
Integer, parameter :: m = 1
Real*8, parameter :: C_Tem = 2.0
Integer, parameter :: n = 5113
Integer, parameter :: nn = 168
!--
End Module

  Module Model_data5
implicit none
Integer, parameter :: m = 1
Real*8, parameter :: C_Tem = 2.0
Integer, parameter :: n = 2192
Integer, parameter :: nn = 72
!--
End Module
!--
!====
! Module for observations
!--
Module Model_data6
implicit none
!--
integer, parameter :: n3 = 72
real*8 Qo_ave(n3)
integer nd(n3)
!--
data Qo_ave(1: n3) /5.86, 5.27, 5.63, 4.89, 22.6, 74.3, 66.3, 91, 65.5, 12.6, 9.4, 5.72,&
                                  5.12,5.36,4.8,9.95,33.6,60,46,120,68	,19,9.53,3.39,&
                                  2.93,2.94,4.97,5.58,21.4,61.1,117,82.2	,49.6,16.6,12,5.87,&
                                  3.09,2.35,4.05,8.01,15.5,44.7,99.7	,103,65,15.4,9.35,7.15,&
                                  6.72,6.38,4.09,6.95,6.9,58.6,98.4,120,36.8,10.2,7.29,6.07,&
                                  5.54,4.41,3.55,7.39,30.9,66.3,139,68.2	,38.6,15.6,10.5,7.43/

data nd(1: n3) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,&
                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!--
End Module  
    !--
!--Parameter description for the DREAMzs module
!-----------------------------------------------------------------------------------------------
!--d, dimension of parameter set
!--n, number of initial Markov Chains (parallel chain)
!--limit, the maximum length of burn-in period
!--nrc, user defined parameter for deriving CR probability distribution
!--num, number of iternations after burn-in period
!--sig, number of chain pairs used for sample generation
!--b1,b2, are the b and b* in equation 4, Ud(-b,b), Nd(0,b*)
!    the use of b2 is different from Vrugt at here
!--st, the start time to test convergence
!--ilen is the initial length of burn_in chain, for obtaining initial points and archive Z
!--CRUpdate is the start time to update CR distribution
!--SR, the convergence criterion.
!--M0, the length of arichive matrix Z[M0, d]
!--TK, thinning rate
!--Psk, the percentage of snooker update, which must has less than or equal 2 decimals.
!--Out_step the step size of output
!--Status = 0 for ending of burn-in, 1 for after burn-in, -1 for before burn-in
!----------------------------------------------------------------------------------
Module Model_data4
implicit none
!----
integer, parameter :: d4=8
integer, parameter :: n4=3
integer, parameter :: limit=5000
integer, parameter :: num=5000
Integer, parameter :: Out_step = num/20
!----
integer, parameter :: ncr=3
integer, parameter :: sig=1
integer, parameter :: st=limit/2
integer, parameter :: M0=10*d4
integer, parameter :: TK=10
integer, parameter :: CRUpdate=limit/10
integer, parameter :: ilen=M0/n4+1
!--
real*8, parameter :: b1=0.05
real*8, parameter :: b2=0.0001
real*8, parameter :: SR=1.2
real*8, parameter :: Psk=0.1D+00
!--
End Module
!====
!=============================================================
Program main
Use Model_data4

!--
implicit none
Integer Status
logical alive
!--
inquire (file = 'Status.txt', exist = alive)
if (alive) then
!--
  open (unit=101, file='Status.txt', recl=50, blank='null', position='rewind', pad='yes')
  read (unit=101, fmt='(i1)') Status
  close (101)
!--
else
  Status = -1
end if
!--
if (Status == -1) then
  Call DREAM_1(d4, ncr, n4, sig, limit, st, ilen, CRUpdate, b1, b2, SR, M0, TK, Psk)
  Call DREAM_2(d4, ncr, n4, sig, limit, num, b1, b2, M0, TK, Psk, Out_step)
else if (Status == 1 .or. Status == 0) then
  Call DREAM_2(d4, ncr, n4, sig, limit, num, b1, b2, M0, TK, Psk, Out_step)
end if
!--
Stop
End Program