!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> External routines used for gaussian integration 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module gauss

  implicit none

  private

  public :: legendre_ek_compute, legendre_gw_compute

  contains

subroutine imtql2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors
!    of a symmetric tridiagonal matrix by the implicit QL method.
!    The eigenvectors of a full symmetric matrix can also
!    be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer  N, the order of the matrix.
!
!    Input/output, real  D(N).  On input, the diagonal elements of
!    the input matrix.  On output, the eigenvalues in ascending order.  If an
!    error exit is made, the eigenvalues are correct but
!    unordered for indices 1,2,...,IERR-1.
!
!    Input/output, real  E(N).  On input, the subdiagonal elements
!    of the input matrix in E(2:N).  E(1) is arbitrary.  On output, E is
!    overwritten.
!
!    Input/output, real  Z(N,N).  On input, the transformation
!    matrix produced in the reduction by TRED2, if performed.  If the
!    eigenvectors of the tridiagonal matrix are desired, Z must contain the
!    identity matrix.  On output, Z contains orthonormal eigenvectors of the
!    symmetric tridiagonal (or full) matrix.  If an error exit is made, Z
!    contains the eigenvectors associated with the stored eigenvalues.
!
!    Output, integer  IERR, error flag.
!    0, for normal return,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  use global, only : r_tiny
  integer  :: n

  real     ::  b
  real     ::  c
  real     ::  d(n)
  real     ::  e(n)
  real     ::  f
  real     ::  g
  integer  i
  integer  ::  ierr
  integer  ::  ii
  integer  ::  j
  integer  ::  k
  integer  ::  l
  integer  ::  m
  integer  ::  mml
  real     ::  p
  real     ::  r
  real     ::  s
  real     ::  t(n)
  real     ::  tst1
  real     ::  tst2
  real     ::  z(n,n)

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do
  e(n) = 0.0E+00

  do l = 1, n

    j = 0
!
!  Look for a small sub-diagonal element.
!
105 continue

    do m = l, n

      if ( m == n ) then
        exit
      end if

      tst1 = abs ( d(m) ) + abs ( d(m+1) )
      tst2 = tst1 + abs ( e(m) )

      if (abs(tst2 - tst1) < r_tiny) then
        exit
      end if

    end do

    p = d(l)

    if ( m == l ) then
      cycle
    end if

    if ( 30 <= j ) then
      ierr = l
      return
    end if

    j = j + 1
!
!  Form shift.
!
    g = ( d(l+1) - p ) / ( 2.0E+00 * e(l) )
    r = pythag ( g, 1.0E+00 )
    g = d(m) - p + e(l) / ( g + sign ( r, g ) )
    s = 1.0E+00
    c = 1.0E+00
    p = 0.0E+00
    mml = m - l

    do ii = 1, mml

      i = m - ii
      f = s * e(i)
      b = c * e(i)
      r = pythag ( f, g )
      e(i+1) = r
!
!  Recover from underflow.
!
      if (abs(r) < r_tiny) then
        d(i+1) = d(i+1) - p
        e(m) = 0.0E+00
        go to 105
      end if

      s = f / r
      c = g / r
      g = d(i+1) - p
      r = ( d(i) - g ) * s + 2.0E+00 * c * b
      p = s * r
      d(i+1) = g + p
      g = c * r - b
!
!  Form vector.
!
      do k = 1, n
        f = z(k,i+1)
        z(k,i+1) = s * z(k,i) + c * f
        z(k,i) = c * z(k,i) - s * f
      end do

    end do

    d(l) = d(l) - p
    e(l) = g
    e(m) = 0.0E+00
    go to 105

  end do
!
!  Order eigenvalues and eigenvectors.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      t(1:n)   = z(1:n,i)
      z(1:n,i) = z(1:n,k)
      z(1:n,k) = t(1:n)

    end if

  end do

end subroutine imtql2


subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer  N, the order of the matrix.
!
!    Input/output, real  D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real  E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real  Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  integer  ::  n

  real     ::  b
  real     ::  c
  real     ::  d(n)
  real     ::  e(n)
  real     ::  f
  real     ::  g
  integer  ::  i
  integer  ::  ii
  integer, parameter :: itn = 30
  integer  ::  j
  integer  ::  k
  integer  ::  l
  integer  ::  m
  integer  ::  mml
  real     ::  p
  real     ::  prec
  real     ::  r
  real     ::  s
  real     ::  z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0E+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop 1
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0E+00 * e(l) )
      r =  sqrt ( g * g + 1.0E+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0E+00
      c = 1.0E+00
      p = 0.0E+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0E+00 )
          e(i+1) = f * r
          s = 1.0E+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0E+00 )
          e(i+1) = g * r
          c = 1.0E+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0E+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0E+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

end subroutine imtqlx

subroutine legendre_ek_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_EK_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer  N, the order.
!
!    Output, real  X(N), the abscissas.
!
!    Output, real  W(N), the weights.
!
  integer  ::  n

  real     ::  bj(n)
  integer  ::  i
  real     ::  w(n)
  real     ::  x(n)
  real     ::  zemu
!
!  Define the zero-th moment.
!
  zemu = 2.0E+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * i, kind = kind(0.E0) )   /                            &
          & real ( 4 * i * i - 1, kind = kind(0.E0) )
  end do

  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0E+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0E+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

end subroutine legendre_ek_compute

subroutine legendre_gw_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_GW_COMPUTE: Legendre quadrature rule by the Golub-Welsch method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Gene Golub, John Welsch,
!    Calculation of Gaussian Quadrature Rules,
!    Mathematics of Computation,
!    Volume 23, Number 106, April 1969, pages 221-230.
!
!  Parameters:
!
!    Input, integer  N, the order.
!
!    Output, real  X(N), the abscissas.
!
!    Output, real  W(N), the weights.
!
  integer  ::  n

  real     ::  beta(n)
  real     ::  d(n)
  integer  ::  i
  integer  ::  ierr
  real     ::  w(n)
  real     ::  x(n)
  real     ::  z(n,n)

  do i = 1, n - 1
    beta(i) = 1.0E+00 / real ( 2 * i, kind = kind(0.E0) )
  end do
  beta(1:n-1) = 0.5E+00 &
   / sqrt ( ( 1.0E+00 - beta(1:n-1) ) * ( 1.0E+00 + beta(1:n-1) ) )
!
!  Compute eigenvalues and eigenvectors.
!
  d(1:n) = 0.0E+00
  beta(2:n) = beta(1:n-1)
  beta(1) = 0.0E+00

  z(1:n,1:n) = 0.0E+00
  do i = 1, n
    z(i,i) = 1.0E+00
  end do

  call imtql2 ( n, d, beta, z, ierr )
!
!  X values are eigenvalues.
!
  x(1:n) = d(1:n)
!
!  W is related to first eigenvector.
!
  w(1:n) = 2.0E+00 * z(1,1:n)**2

  w(1:n) = 2.0E+00 * w(1:n) / sum ( w(1:n) )

end subroutine legendre_gw_compute

function pythag ( a, b )

!*****************************************************************************80
!
!! PYTHAG computes SQRT ( A * A + B * B ) carefully.
!
!  Discussion:
!
!    The formula
!
!      pythag = sqrt ( a * a + b * b )
!
!    is reasonably accurate, but can fail if, for example, A*A is larger
!    than the machine overflow.  The formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, real  A, B, the two legs of a right triangle.
!
!    Output, real  PYTHAG, the length of the hypotenuse.
!
  use global, only : r_tiny
  real     ::  a
  real     ::  b
  real     ::  p
  real     ::  pythag
  real     ::  r
  real     ::  s
  real     ::  t
  real     ::  u

  p = max ( abs ( a ), abs ( b ) )

  if (.not. abs(p) < r_tiny) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = 4.0E+00 + r

      if ( abs(t - 4.0E+00) < r_tiny ) then
        exit
      end if

      s = r / t
      u = 1.0E+00 + 2.0E+00 * s
      p = u * p
      r = ( s / u ) * ( s / u ) * r

    end do

  end if

  pythag = p

end function pythag

end module gauss
