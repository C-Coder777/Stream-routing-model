!================================================================
subroutine multinormal_sample (m, n, a, mu, seed, x)

!*****************************************************************************
!
!    MULTINORMAL_SAMPLE samples a multivariate normal distribution.
!
!    * M, the spatial dimension,
!    * N, the number of points to generate,
!    * SEED, the seed, a positive integer.
!    * MU, the mean vector.
!    * A, the MxM variance-covariance matrix.
!
!  Discussion:
!    The multivariate normal distribution for the M dimensional vector X has the form:
!      pdf(X) = (2*pi*det(A))**(-M/2) * exp(-0.5*(X-MU)'*inverse(A)*(X-MU))
!    where MU is the mean vector, and A is a positive definite symmetric matrix called the variance-covariance matrix.
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!  Modified:
!    07 December 2009
!  Author:
!    John Burkardt
!
!  Parameters:
!    Input, integer ( kind = 4 ) M, the dimension of the space.
!    Input, integer ( kind = 4 ) N, the number of points.
!    Input, real ( kind = 8 ) A(M,M), the variance-covariance 
!    matrix.  A must be positive definite symmetric.
!    Input, real ( kind = 8 ) MU(M), the mean vector.
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!    Output, real ( kind = 8 ) X(M, N), the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) a(m,m)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu(m)
  real ( kind = 8 ) r(m,m)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(m,n)
!
!  Compute the upper triangular Cholesky factor R of the variance-covariance matrix.
!
  r(1:m,1:m) = a(1:m,1:m)

  call r8po_fa ( m, r, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MULTINORMAL_SAMPLE - Fatal error!'
    write ( *, '(a)' ) &
      '  The variance-covariance matrix is not positive definite symmetric.'
    stop
  end if
!
!  Get an MxN matrix of samples of the 1D normal distribution with mean 0 and variance 1.  
!
  call r8vec_normal_01 ( m*n, seed, x(1:m,1:n) )
!
!  Compute R' * X.
!  We actually carry out this computation in the equivalent form X' * R.
!
  do j = 1, n
    x(1:m,j) = mu(1:m) + matmul ( x(1:m,j), r(1:m,1:m) )
  end do

  return
end
!------------------------------------------------------------------------------------------------------------------
subroutine r8po_fa (n, a, info )
!*****************************************************************************
! R8PO_FA factors an R8PO matrix.
!  Discussion:
!    The R8PO storage format is used for a symmetric positive definite matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.) Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero. R8PO storage is used by LINPACK and LAPACK.
!
!    The positive definite symmetric matrix A has a Cholesky factorization of the form:
!    A = R' * R
!    where R is an upper triangular matrix with positive elements on its diagonal.  This routine overwrites the matrix A with its factor R.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!    Modified:  22 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!  Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart, LINPACK User's Guide, SIAM, 1979, ISBN13: 978-0-898711-72-1, LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the matrix in R8PO storage.
!    On output, the Cholesky factor R in R8GE storage.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal return.
!    K, error condition.  The principal minor of order K is not
!    positive definite, and the factorization was not completed.
!--------------------------------------------------------------------------
  implicit none
  integer n, i, info, j, k
  real*8 a(n,n),s
  do j = 1, n
    do k = 1, j-1
      a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
    end do
    s = a(j,j) - sum( a(1:j-1,j)**2 )
    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if
    a(j,j) = sqrt ( s )
  end do
  info = 0
!
!  Since the Cholesky factor is stored in R8GE format, be sure to zero out the lower triangle.
!
  do i = 1, n
    do j = 1, i-1
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
!--------------------------------------------------------------------------------------------------------------------
subroutine r8vec_normal_01 ( n, seed, x )
!*****************************************************************************
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over. In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    17 July 2006
!  Author:
!    John Burkardt
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer MADE, records the number of values that have been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1), and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if SAVED is 1.
!
  implicit none
  integer n, m, seed
  integer, save :: made = 0
  real*8, parameter :: pi = 3.141592653589793D+00
  real*8 r(n+1), r8_uniform_01, x(n)
  integer, save :: saved = 0
  integer x_hi_index, x_lo_index
  real*8, save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =  sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2 * m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) * cos ( 2.0D+00 * pi * r(2:2*m:2) )
    x(x_lo_index+1:x_hi_index:2) = sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) )  * sin ( 2.0D+00 * pi * r(2:2*m:2) )
    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1
    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2 * m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )
    x(x_lo_index+1:x_hi_index:2) = sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) * cos ( 2.0D+00 * pi * r(2*m) )
    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1
    made = made + x_hi_index - x_lo_index + 2

  end if
  return
end
!--------------------------------------------------------------------------------------------------------
function r8_uniform_01 ( seed )
!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!    An R8 is a real ( kind = 8 ) value.
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits, including a sign bit.
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!  Modified:
!    05 July 2006
!  Author:
!    John Burkardt
!  Reference:
!    Paul Bratley, Bennett Fox, Linus Schrage, A Guide to Simulation, Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer, Random Number Generation, in Handbook of Simulation, edited by Jerry Banks, Wiley Interscience, page 95, 1998.
!    Bennett Fox, Algorithm 647: Implementation and Relative Efficiency of Quasirandom Sequence Generators,
!    ACM Transactions on Mathematical Software, Volume 12, Number 4, pages 362-376, 1986.
!    Peter Lewis, Allen Goodman, James Miller A Pseudo-Random Number Generator for the System/360, IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773
  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer, it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
!--------------------------------------------------------------------------------------------------
subroutine r8vec_uniform_01 ( n, seed, r )
!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!    An R8VEC is a vector of real ( kind = 8 ) values.
!    For now, the input quantity SEED is an integer variable.
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!  Modified:
!    05 July 2006
!  Author:
!    John Burkardt
!  Reference:
!    Paul Bratley, Bennett Fox, Linus Schrage, A Guide to Simulation, Springer Verlag, pages 201-202, 1983.
!    Bennett Fox, Algorithm 647: Implementation and Relative Efficiency of Quasirandom Sequence Generators,
!    ACM Transactions on Mathematical Software, Volume 12, Number 4, pages 362-376, 1986.
!    Peter Lewis, Allen Goodman, James Miller A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal, Volume 8, pages 136-143, 1969.
!
!  Parameters:
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should NOT be 0.  On output, SEED has been updated.
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773
    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
!==============================================================
subroutine genmul ( n, p, ncat, ix )
!*****************************************************************************
!! GENMUL generates a multinomial random deviate.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    31 March 2013
!  Author:
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!  Reference:
!    Luc Devroye, Non-Uniform Random Variate Generation, Springer, 1986, ISBN: 0387963057, LC: QA274.D48.
!
!  Parameters:
!    Input, integer ( kind = 4 ) N, the number of events, which will be classified into one of the NCAT categories.
!    Input, real ( kind = 4 ) P(NCAT-1).  P(I) is the probability that an event will be classified into category I.  Thus, each P(I) must be between 
!    0.0 and 1.0.  Only the first NCAT-1 values of P must be defined since  P(NCAT) would be 1.0 minus the sum of the first NCAT-1 P's.
!
!    Input, integer ( kind = 4 ) NCAT, the number of categories.
!    Output, integer ( kind = 4 ) IX(NCAT), a random observation from  the multinomial distribution.  All IX(i) will be nonnegative and their 
!    sum will be N.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icat
  integer ( kind = 4 ) ignbin
  integer ( kind = 4 ) ix(ncat)
  integer ( kind = 4 ) ntot
  real ( kind = 8 ) p(ncat-1)
  real ( kind = 8 ) prob
  real ( kind = 8 ) ptot

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENMUL - Fatal error!'
    write ( *, '(a)' ) '  N < 0'
    stop
  end if

  if ( ncat <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENMUL - Fatal error!'
    write ( *, '(a)' ) '  NCAT <= 1'
    stop
  end if

  do i = 1, ncat - 1

    if ( p(i) < 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENMUL - Fatal error!'
      write ( *, '(a)' ) '  Some P(i) < 0.'
      stop
    end if

    if ( 1.0E+00 < p(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENMUL - Fatal error!'
      write ( *, '(a)' ) '  Some 1 < P(i).'
      stop
    end if

  end do

  ptot = 0.0E+00
  do i = 1, ncat - 1
    ptot = ptot + p(i)
  end do

  if ( 0.99999E+00 < ptot ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENMUL - Fatal error!'
    write ( *, '(a)' ) '  1 < Sum of P().'
    stop
  end if
!
!  Initialize variables.
!
  ntot = n
  ptot = 1.0E+00
  do i = 1, ncat
    ix(i) = 0
  end do
!
!  Generate the observation.
!
  do icat = 1, ncat - 1
    prob = p(icat) / ptot
    ix(icat) = ignbin ( ntot, prob )
    ntot = ntot - ix(icat)
    if ( ntot <= 0 ) then
      return
    end if
    ptot = ptot - p(icat)
  end do

  ix(ncat) = ntot

  return
end
!--------------------------------------------------------------------------------------------------------------------------------
function r4_uniform_01 ()
implicit none
real*8 r4_uniform_01
    call random_number(r4_uniform_01)
return
end
!----------------------------------------------------------------------------------------------------------------------------------
function ignbin ( n, pp )
!*****************************************************************************80
!
!! IGNBIN generates a binomial random deviate.
!
!  Discussion:
!
!    This procedure generates a single random deviate from a binomial distribution whose number of trials is N and whose
!    probability of an event in each trial is P.
!
!    The previous version of this program relied on the assumption that local memory would be preserved between calls.  It set up data
!    one time to be preserved for use over multiple calls.  In the  interests of portability, this assumption has been removed, and
!    the "setup" data is recomputed on every call.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    31 March 2013
!  Author:
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!    Voratas Kachitvichyanukul, Bruce Schmeiser, Binomial Random Variate Generation, Communications of the ACM,
!    Volume 31, Number 2, February 1988, pages 216-222.
!
!  Parameters:
!    Input, integer ( kind = 4 ) N, the number of binomial trials, from which a random deviate will be generated. 0 < N.
!    Input, real ( kind = 4 ) PP, the probability of an event in each trial of the binomial distribution from which a random deviate is to be generated.
!    0.0 < PP < 1.0.
!    Output, integer ( kind = 4 ) IGNBIN, a random deviate from the distribution.
!
  implicit none

  real ( kind = 8 ) al
  real ( kind = 8 ) alv
  real ( kind = 8 ) amaxp
  real ( kind = 8 ) c
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) ffm
  real ( kind = 8 ) fm
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ignbin
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ix1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp
  real ( kind = 8 ) pp
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) p4
  real ( kind = 8 ) q
  real ( kind = 8 ) qn
  real ( kind = 8 ) r
  real ( kind = 8 ) r4_uniform_01
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) w2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xl
  real ( kind = 8 ) xll
  real ( kind = 8 ) xlr
  real ( kind = 8 ) xm
  real ( kind = 8 ) xnp
  real ( kind = 8 ) xnpq
  real ( kind = 8 ) xr
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z
  real ( kind = 8 ) z2

  if ( pp <= 0.0D+00 .or. 1.0D+00 <= pp ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IGNBIN - Fatal error!'
    write ( *, '(a)' ) '  PP is out of range.'
    stop
  end if

  p = min ( pp, 1.0D+00 - pp )
  q = 1.0D+00 - p
  xnp = real ( n, kind = 8 ) * p

  if ( xnp < 30.0D+00 ) then
    qn = q ** n
    r = p / q
    g = r * real ( n + 1, kind = 8 )

    do

      ix = 0
      f = qn
      u = r4_uniform_01 ( )

      do

        if ( u < f ) then
          if ( 0.5D+00 < pp ) then
            ix = n - ix
          end if
          ignbin = ix
          return
        end if

        if ( 110 < ix ) then
          exit
        end if

        u = u - f
        ix = ix + 1
        f = f * ( g / real ( ix, kind = 8 ) - r )

      end do
    end do
  end if

  ffm = xnp + p
  m = ffm
  fm = m
  xnpq = xnp * q
  p1 = int ( 2.195D+00 * sqrt ( xnpq ) - 4.6D+00 * q ) + 0.5D+00
  xm = fm + 0.5D+00
  xl = xm - p1
  xr = xm + p1
  c = 0.134D+00 + 20.5D+00 / ( 15.3D+00 + fm )
  al = ( ffm - xl ) / ( ffm - xl * p )
  xll = al * ( 1.0D+00 + 0.5D+00 * al )
  al = ( xr - ffm ) / ( xr * q )
  xlr = al * ( 1.0D+00 + 0.5D+00 * al )
  p2 = p1 * ( 1.0D+00 + c + c )
  p3 = p2 + c / xll
  p4 = p3 + c / xlr
!
!  Generate a variate.
!
  do

    u = r4_uniform_01 ( ) * p4
    v = r4_uniform_01 ( )
!
!  Triangle
!
    if ( u < p1 ) then
      ix = xm - p1 * v + u
      if ( 0.5D+00 < pp ) then
        ix = n - ix
      end if
      ignbin = ix
      return
    end if
!
!  Parallelogram
!
    if ( u <= p2 ) then

      x = xl + ( u - p1 ) / c
      v = v * c + 1.0D+00 - abs ( xm - x ) / p1

      if ( v <= 0.0D+00 .or. 1.0D+00 < v ) then
        cycle
      end if

      ix = x

    else if ( u <= p3 ) then

      ix = xl + log ( v ) / xll
      if ( ix < 0 ) then
        cycle
      end if
      v = v * ( u - p2 ) * xll

    else

      ix = xr - log ( v ) / xlr
      if ( n < ix ) then
        cycle
      end if
      v = v * ( u - p3 ) * xlr

    end if

    k = abs ( ix - m )

    if ( k <= 20 .or. xnpq / 2.0 - 1.0 <= k ) then

      f = 1.0D+00
      r = p / q
      g = ( n + 1 ) * r

      if ( m < ix ) then
        mp = m + 1
        do i = m + 1, ix
          f = f * ( g / i - r )
        end do
      else if ( ix < m ) then
        ix1 = ix + 1
        do i = ix + 1, m
          f = f / ( g / real ( i, kind = 8 ) - r )
        end do
      end if

      if ( v <= f ) then
        if ( 0.5D+00 < pp ) then
          ix = n - ix
        end if
        ignbin = ix
        return
      end if

    else

      amaxp = ( k / xnpq ) * ( ( k * ( k / 3.0D+00 + 0.625D+00 ) + 0.1666666666666D+00 ) / xnpq + 0.5D+00 )
      ynorm = - real ( k * k, kind = 8 ) / ( 2.0D+00 * xnpq )
      alv = log ( v )

      if ( alv < ynorm - amaxp ) then
        if ( 0.5D+00 < pp ) then
          ix = n - ix
        end if
        ignbin = ix
        return
      end if

      if ( ynorm + amaxp < alv ) then
        cycle
      end if

      x1 = real ( ix + 1, kind = 8 )
      f1 = fm + 1.0D+00
      z = real ( n + 1, kind = 8 ) - fm
      w = real ( n - ix + 1, kind = 8 )
      z2 = z * z
      x2 = x1 * x1
      f2 = f1 * f1
      w2 = w * w

      t = xm * log ( f1 / x1 ) + ( n - m + 0.5D+00 ) * log ( z / w ) + real ( ix - m, kind = 8 ) * log ( w * p / ( x1 * q ) ) &
        + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 - ( 99.0D+00 - 140.0D+00 / f2 ) / f2 ) / f2 ) / f2 ) / f1 / 166320.0D+00 &
        + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 - ( 99.0D+00 - 140.0D+00 / z2 ) / z2 ) / z2 ) / z2 ) / z / 166320.0D+00 &
        + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 - ( 99.0D+00 - 140.0D+00 / x2 ) / x2 ) / x2 ) / x2 ) / x1 / 166320.0D+00 &
        + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 - ( 99.0D+00 - 140.0D+00 / w2 ) / w2 ) / w2 ) / w2 ) / w / 166320.0D+00

      if ( alv <= t ) then
        if ( 0.5D+00 < pp ) then
          ix = n - ix
        end if
        ignbin = ix
        return
      end if
    end if
  end do

  return
end
!--------------------------------------------------------------------------------------------------------------------------------