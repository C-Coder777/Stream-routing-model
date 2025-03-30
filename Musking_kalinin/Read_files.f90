!----
!--Read files
Subroutine Read_files_1(d, n, ncr, limit, num, bnum, im, p, SR_Chain, SR_Con, ZM, num_ZM, paras_burn)
Implicit none
!--
Integer i, j, ncr, limit, num, d, n, bnum, im, num_ZM
Real*8 p(ncr), SR_Con(limit+num, d), SR_Chain(limit+num, d, n), ZM(num_ZM, d), paras_burn(d+1, n)
Character (len = 25) head
!--
open (unit=101, file='DREAMzs_1.txt', recl=50, blank='null', position='rewind', pad='yes')
read (unit=101, fmt='(a20, i6)') head, bnum
read (unit=101, fmt='(a20)') head
read (unit=101, fmt='(a8, <ncr>f10.6)') head, p(1: ncr)
close (101)
!--
open (unit=102, file='Burn_in.txt', recl=500, blank='null', position='rewind', pad='yes') 
read (unit=102, fmt='(a20)') head
do i = 1, n
  do j=1, bnum
    read (unit=102, fmt='(<d>f14.6)') SR_Chain(j, 1: d, i)
  end do
end do
close(102)
!--
open (unit=103, file='SR_Records.txt', recl=500, blank='null', position='rewind', pad='yes') 
read (unit=103, fmt='(a20)') head
do i=1, bnum
    read (unit=103, fmt='(<d>f20.4)') SR_Con(i, 1: d)
end do
close(103)
!--
open (unit=104, file='Arichive_matrix.txt', recl=500, blank='null', position='rewind', pad='yes')
read (unit=104, fmt='(i8)') im
do i = 1, im
  read (unit=104, fmt='(<d>f14.6)') ZM(i, 1: d)
end do
close (104)
!--
open (unit=106, file='Burn_ends.txt', recl=500, blank='null', position='rewind', pad='yes')
do i=1, n
  read (unit=106, fmt='(<d>f14.6, E16.6E3)') paras_burn(1: d+1, i)
end do
close(106)
!--
End subroutine
!--------------------------------------------------------------------------------------------------------------------------
Subroutine Read_files_2 (d, n, ncr, limit, num, c2, i2, bnum, im, p, SR_Chain, SR_Con, ZM, num_ZM, chains)
Implicit none
!--
Integer i, j, bnum, im, ncr, limit, num, d, n, c2, num_ZM, i2
Real*8 p(ncr), SR_Con(limit+num, d), SR_Chain(limit+num, d, n), ZM(num_ZM, d)
Real*8 chains(num, d+1, n)
Character (len = 25) head
!--
open (unit=101, file='DREAMzs_1.txt', recl=50, blank='null', position='rewind', pad='yes')
read (unit=101, fmt='(a20, i6)') head, bnum
read (unit=101, fmt='(a20)') head
read (unit=101, fmt='(a8, <ncr>f10.6)') head, p(1: ncr)
close (101)
!--
open (unit=102, file='Burn_in.txt', recl=500, blank='null', position='rewind', pad='yes') 
read (unit=102, fmt='(a33)') head
do i = 1, n
  do j=1, bnum
    read (unit=102, fmt='(<d>f14.6)') SR_Chain(j, 1:d, i)
  end do
end do
close(102)
!--
open (unit=103, file='SR_Records.txt', recl=500, blank='null', position='rewind', pad='yes') 
read (unit=103, fmt='(a20)') head
do i=1, c2
    read (unit=103, fmt='(<d>f20.4)') SR_Con(i, 1: d)
end do
close(103)
!--
open (unit=104, file='Arichive_matrix.txt', recl=500, blank='null', position='rewind', pad='yes')
read (unit=104, fmt='(2i8)') im, i2
do i = 1, im
  read (unit=104, fmt='(<d>f14.6)') ZM(i, 1: d)
end do
close (104)
!--
open (unit=107, file='chains.txt', recl=500, blank='null', position='rewind', pad='yes') 
read (unit=107, fmt='(a20)') head
do i=1, n
  do j=1, c2-bnum
    read (unit=107, fmt='(<d>f14.6, E16.6E3)') chains(j, 1: d+1, i)
    SR_Chain(bnum+j, 1: d, i) = chains(j, 1: d, i)
  end do
end do
close (107)
!--
!====
End subroutine