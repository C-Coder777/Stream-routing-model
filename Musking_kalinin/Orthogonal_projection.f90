!-------------------------------------------------------------------
! Orthogonal projection of points by a line
!----
!--line(2, d) is the start and end points of a line
!--points(n, d) are the points need to be pojected.
!--p_points(n, d) are hte projected points
!--
Subroutine Ortho_Pro(line, points, d, n1, p_points, n2)
implicit none
integer i, j, d, n1, n2
real*8 point(d), line(2, d), points(n2, d), p_points(n2, d), t, s1, s2
!--
do i=1, n2
  s1=0.0
  s2=0.0
  do j=1, d
    s1=s1+(line(1, j)-points(i, j))*(line(2, j)-line(1, j))
    s2=s2+(line(2, j)-line(1, j))*(line(2, j)-line(1, j))
  end do
  t=-s1/s2
  do j=1, d
    p_points(i, j)=line(1, j)+(line(2, j)-line(1, j))*t
  end do
end do
!--
End subroutine