!====
! Used to generated initial samples
!--
Subroutine Pre_samples(samples, n, mean, cov, d)
implicit none
integer i, n, d, seed
real*8 samples(n, d), mean(d), cov(d, d), samples2(d, n)
!--
seed=1024768
!--
call multinormal_sample (d, n, cov, mean, seed, samples2)
!--
do i=1, n
  samples(i, :)=samples2(:, i)
end do
!--
End subroutine
!==========================================
!---------------------------judge whether the alternative point beyonds the range or not
subroutine range_jud(ranges, d, col, point, label)
implicit none
integer i, j, d, col, label, lab(d)
real*8 ranges(d, col), point(d)
label=1
!----
do i=1, d
  if(ranges(i, 1) <= point(i) .and. point(i) <= ranges(i, 2)) then
    lab(i)=1
  else
    lab(i)=0
  end if
end do
!--
do i=1,d
  label=label*lab(i)
end do
end subroutine
!--------------------------------------------------------------------------------


