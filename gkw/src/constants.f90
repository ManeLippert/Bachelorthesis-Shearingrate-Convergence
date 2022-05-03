!-----------------------------------------------------------------------------
!> Contains some physical and mathematical constants
!-----------------------------------------------------------------------------
module constants 

  implicit none

!
! ========== mathematical constants ==========
!

  !> \f$ \pi \f$
  real, parameter :: pi =  3.1415926535897932384626433832795028841971693993751

  !> The imaginary number i=sqrt(-1) 
  complex, parameter :: ci1 = (0.,1.)
  !> the complex identity
  complex, parameter :: c1 = (1.,0.)
  !> The complex zero
  complex, parameter :: czero = (0.0, 0.0)
  
!
! =========== physical constants ==========
!

  !> The mass of a proton 
  real, parameter :: proton_mass = 1.67262158E-27

  !> \f$ \epsilon_0 \f$, permittivity of free space
  real, parameter :: epsilon_0 = 8.85418782E-12

  !> The unit of charge 
  real, parameter :: unit_charge_si = 1.602176487E-19

end module constants 
