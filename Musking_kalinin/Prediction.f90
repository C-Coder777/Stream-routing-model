!----Program for prediction
!--
Subroutine Prediction(chains, Q_pre, num, d, n, Qm)
Implicit none
!--
Integer i, j, num, d, n, Qm
Real*8 chains(num, d+1, n), Q_pre(num, Qm, n), Q_pre1(Qm), point(d)
!--
Do i = 1, num
  if (i > 100 .and. Mod(i, 100) == 0) then
    write(*, fmt = '(a31, i8)') 'Iteration number of predicion:', i
  end if
  Do j = 1, n
    point(1: d) = chains(i, 1: d, j)
    call Run_model_pre(point, d, Q_pre1)

    Q_pre(i, 1: Qm, j) = Q_pre1(1: Qm)
  End do
End do
!--
End subroutine
