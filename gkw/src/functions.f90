!-----------------------------------------------------------------------------
!> Contains mathematical functions required by GKW.
!> Bessel and Gamma functions
!-----------------------------------------------------------------------------
module functions

  private

  public :: gamma_gkw, gamma1_gkw, besselj0_gkw, mod_besselj1_gkw, mod_sinh_gkw
  public :: legendre, leguerre, cheby

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This function calculates the Gamma function needed in Poissons and
!> Ampere's equation
!----------------------------------------------------------------------------

function gamma_gkw(imod,ix,i,is)

  use specfun,    only : expbessi0
  use components, only : mas, vthrat, signz  
  use geom,       only : bn
  use mode,       only : krloc

  integer, intent(in) :: imod, ix, i, is

  real :: gamma_gkw, dum1

  dum1 = 0.5*(mas(is)*vthrat(is)*krloc(imod,ix,i)/(signz(is)*bn(ix,i)))**2 
  gamma_gkw = expbessi0(dum1)

end function gamma_gkw

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This function calculates the first order Gamma function needed in Poissons
!> and Ampere's equation
!----------------------------------------------------------------------------

function gamma1_gkw(imod,ix,i,is)

  use specfun,    only : expbessi1
  use components, only : mas, vthrat, signz  
  use geom,       only : bn
  use mode,       only : krloc

  integer, intent(in) :: imod, ix, i, is

  real :: gamma1_gkw, dum1

  dum1 = 0.5*(mas(is)*vthrat(is)*krloc(imod,ix,i)/(signz(is)*bn(ix,i)))**2
  gamma1_gkw = expbessi1(dum1)

end function gamma1_gkw

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This function calculates the bessel function J0 used in the Gyro average
!----------------------------------------------------------------------------

function besselj0_gkw(imod,ix,i,j,is)

  use specfun,      only : bessj0
  use components,   only : mas, vthrat, signz  
  use velocitygrid, only : mugr
  use geom,         only : bn
  use mode,         only : krloc

  integer, intent(in) :: imod, ix, i, j, is 

  real :: besselj0_gkw

  besselj0_gkw = bessj0(mas(is)*vthrat(is)*krloc(imod,ix,i)*                & 
               & sqrt(2.E0*mugr(j)/bn(ix,i)) / signz(is) )

end function besselj0_gkw


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This function calculates the bessel function 2*J1(x)/x used in the Gyro
!> average
!----------------------------------------------------------------------------

function mod_besselj1_gkw(imod,ix,i,j,is)

  use specfun,      only : bessj1
  use components,   only : mas, vthrat, signz  
  use velocitygrid, only : mugr
  use geom,         only : bn
  use mode,         only : krloc

  integer, intent(in) :: imod, ix, i, j, is 

  real :: mod_besselj1_gkw

  !In case of zero mode the following limit is taken:
  !  lim( 2*J_1(x)/x ) = 1/2 as x->0

  if (abs(krloc(imod,ix,i)) < 1.0E-5) then
    mod_besselj1_gkw = 0.5
  else
    mod_besselj1_gkw = 2.0*bessj1(mas(is)*vthrat(is)*krloc(imod,ix,i) *         & 
               & sqrt(2.E0*mugr(j)/bn(ix,i))/signz(is))/(mas(is)*krloc(imod,ix,i)* & 
               & vthrat(is)*sqrt(2.E0*mugr(j)/bn(ix,i))/signz(is))
  end if

end function mod_besselj1_gkw

function mod_sinh_gkw(input_sinh)
  
  real, intent(in)  :: input_sinh
  real, parameter   :: input_sinh_tol = 1.E-2
  real              :: mod_sinh_gkw

  !This function determines the value of sinh(x)/x. For the limit x->0,
  !a Taylor expansion is used

  if (abs(input_sinh).le.input_sinh_tol) then
     mod_sinh_gkw = 1.E0 + 1.E0/6.E0*input_sinh**2 + 1.E0/120.E0*input_sinh**4
  else
     mod_sinh_gkw = sinh(input_sinh)/input_sinh
  end if
  
end function mod_sinh_gkw

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>Function to evaluate the Legendre polynomial for integer order, n and 
!> argument x
!> Using Bonnets recursion formula for computation
!> (n + 1) P_{n+1} (x) = (2n +1) x P_{n} (x) - n P_{n-1} (x)
!----------------------------------------------------------------------------
function legendre(n, x)

  use general, only : gkw_abort

  integer :: n,j
  real :: x
  real :: legendre
  real :: p0,p1,p2

  if(x.lt.-1.0 .or. x.gt.1.0)then
    call gkw_abort('Argument of Legendre polynomial out of range')
  endif

  if(n.eq.0)then
    legendre = 1.E0
    return
  else if (n.eq.1)then
    legendre = x
    return
  else if(n.lt.0)then
    return
  else
    p0 = 1
    p1 = x
    ! note: j = n+1 compared with the form given above.
    do j=2,n
      p2 = (2*j-1)*x*p1 - (j-1)*p0
      p2 = p2/j  
      p0 = p1
      p1 = p2
    end do
    legendre = p2
    return
  endif
  
end function legendre

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Function to evaluate the Leguerre polynomial of condition alpha to order x
!> with argument x
!> Using the recursive formula
!> \f[ L_{j}^{(\alpha)} (x) = ( (2 j - 1 + \alpha - x) L_{j-1}^{(\alpha)} (x) - (j - 1 + \alpha) L_{j-2}^{(\alpha)} (x) ) / j \f]
!----------------------------------------------------------------------------
function leguerre(n, x, alpha)


  integer :: n,j
  real :: x, alpha
  real :: leguerre
  real :: p0,p1,p2

  if(n.eq.0)then
    leguerre = 1.E0
    return
  else if (n.eq.1)then
    leguerre = 1.0E0 + alpha - x
    return
  else if(n.lt.0)then
    return
  else
    p0 = 1.0
    p1 = 1.0 + alpha - x
    do j=2,n
      p2 = (2*j  - 1 + alpha - x)*p1 - (j-1+alpha)*p0
      p2 = p2/j
      p0 = p1
      p1 = p2
    end do
    leguerre = p2
    return
  endif

end function leguerre

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>Function to evaluate the Chebyschev polynomial of the first
!>kind using its generating function Tn+1=2xTn - Tn-1, for order n
!>and argument x
!----------------------------------------------------------------------------
function cheby(n, x)

  use general, only : gkw_abort

  integer :: n,j
  real :: x
  real :: cheby
  real :: p0,p1,p2

  if(x.lt.-1.0 .or. x.gt.1.0)then
    call gkw_abort('Argument of Chebyshev polynomial out of range')
  endif

  if(n.eq.0)then
    cheby = 1.E0
    return
  else if (n.eq.1)then
    cheby = 1.0E0 - x
    return
  else if(n.lt.0)then
    return
  else
    p0 = 1.0
    p1 = x
    do j=2,n
      p2 = 2*x*p1 - p0
      p0 = p1
      p1 = p2
    end do
    cheby = p2
    return
  endif

end function cheby

end module functions
