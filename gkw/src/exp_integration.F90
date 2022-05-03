!-----------------------------------------------------------------------------------
!> This module handles the main time integration, it is the core of the solver
!> The code spends most of its time here.
!> It contains a choice of explicit time integrators (RK2, RK4, 3rd order)
!> which all call a generalised calculate_RHS and calculate_feilds
!> It contains the matrix vector multiply of the linear terms
!> This module is optimised, contains OpenMP, and (optional) perflib timings
!----------------------------------------------------------------------------------
module exp_integration

implicit none

#ifdef blas
#ifdef real_precision_default
  external :: blas_copy => ccopy
  external :: blas_axpy => caxpy
#else
  external :: blas_copy => zcopy
  external :: blas_axpy => zaxpy
#endif
#endif

  private

  public :: advance_large_step_explicit, init_explicit, exp_integration_deallocate
  public :: calculate_rhs, set_persistent_mode

  complex, allocatable, public, save :: persistent_mode(:)

  !> parameter used to set the implicitness in some of the schemes
  real, save :: delta

  !> parameters used for the 3rd order (or is it second ??) scheme
  !> that is stable for waves
  real, dimension(3), save :: alf, bet
  real, save :: gam

  !> integers corresponding to the number of temporary solution copies needed
  !> in explicit time integration
  integer, save, public :: isizef, isizer

  !> logical to flip once nonlinear timestep estimator is first used
  logical, save :: clean_dtim=.true.

contains


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Initialization routine for the explicit integration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine init_explicit
use control,   only : meth
use dist,      only : nsolc, fdisk, rhsk
use mpighosts, only : persistent_comm_init
use general,   only : gkw_abort

! for error
integer :: ierr

! This subroutine need also be called by 'EIV'.
!if (method/='EXP')  return

! Set the parameters for the 3rd order scheme
! second order (for completeness is written out here)
alf(1) = 2.E0
alf(2) = -0.5E0
alf(3) = 0.E0
bet(1) = 2.E0
bet(2) = -1.E0
bet(3) = 0.E0
gam = 1.5E0

! Third order (is actually used)
alf(1) = 3.E0
alf(2) = -1.5E0
alf(3) = 1.E0/3.E0
bet(1) = 3.E0
bet(2) = -3.E0
bet(3) = 1.E0
gam = 11.E0/6.E0

! set the implicitness parameter
delta = 0.5E0

! allocate the arrays of the distribution function and the rhs
! note the size depends on the scheme used and therefore on meth
select case(meth)
case(1)
  isizef = 1
  isizer = 1
case(2)
  isizef = 2
  isizer = 1
case(3)
  isizef = 4
  isizer = 3
case(-2,-3,-4,-5,-6,-7,-8)
  isizef = 3
  isizer = 2
case default
  call gkw_abort('exp_integration : Unknown explicit integration scheme')
end select

  ierr = 0
  allocate(fdisk(nsolc,isizef),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate fdisk in exp_integration')

  ! allocate the right hand side
  allocate(rhsk(nsolc,isizer),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate rhsk in exp_integration')

  ! Set up the persistent communication for exp integration
  ! and and nonspectral fields calculation
  ! needs only mpicomms, grid, and dist to be setup (could go higher up)
  ! The actual communication is started in calculate_rhs
  ! and field_solve_nonspec_wrap
  call persistent_comm_init()

end subroutine init_explicit

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine does the explicit time integration
!> Several methods can be used. Their choice is
!> controlled through the 'meth' parameter
!> meth = 1 Midpoint method
!> meth = 2 Fourth order Runga Kutta
!> meth = 3 Third order scheme that is stable for
!>          waves
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine advance_large_step_explicit(itime, call_diagnostics)
use dist,       only : fdisi
use control,    only : meth, dtim, time, naverage,       &
                     & dtim_input, dtim_est, dtim_est_save, non_linear, &
                     & ntotstep, nl_dtim_est, stop_me, dt_min
use control,    only : spectral_radius, last_smallstep_time
use rotation,   only : shear_shift_ky, shear_ky_shift, shear_remap, wavevector_remap
use general,    only : gkw_warn
use diagnostic, only : diagnostic_pre_naverage
use fields,     only : calculate_cons, calculate_fields, field_solve_nonspec_wrap
use mpiinterface, only : mpiallreduce_min, root_processor
use perform,      only : perfloop, perfdo, perf_measure
use control,      only : method
integer, intent(in) :: itime
logical, intent(in) :: call_diagnostics
integer :: iloop
real :: dtim_dum

  ! the loop over the small timesteps
  time_stepping : do iloop = 1, naverage

    if(perf_measure) then
      call perfloop(2)
      if (iloop <= 2 .and. itime == 1.and. method=='EXP') then
        perfdo=.false.  ! Skip the first two loop iterations for code section timings
      else
        perfdo=.true.
      end if
    end if

    !By doing the shearing here we do not need to shift the potential
    !Since the potential is recalculated in calculate_rhs
    if (shear_ky_shift) call shear_shift_ky(fdisi)
    if (shear_remap) call wavevector_remap(fdisi,time)

    ! select the integration method
    method_of_integration : select case (meth)
      case(1); call rk2
      case(2); call rk4
      case(3); call rk3(itime,iloop)
      case(-2,-3,-4,-5,-6,-7); call rkc_2(abs(meth))
      case default
        ! Meth does not match any of the implemented methods
        stop 'Not a proper selection of the numerical method [METH]'
    end select method_of_integration
    
    ! Set persistent mode before calling diagnostics. This ensures that
    ! the last update of the persistent mode in the integration routine
    ! chosen is not considered in the diagnostics output.
    call set_persistent_mode(fdisi)

    ! reset the time step if necessary - will not work for meth=3?
    if (meth /= 3) then
      if (nl_dtim_est) then
        if (dtim_est < dtim_input) then
          !if (root_processor) write(*,*) 'Timestep estimator is off'
          if (root_processor.and.clean_dtim) then
               write (*,*)
               write (*,*) 'WARNING: NL timestep estimator activated'
               write (*,*) 'Timestep was set from dtim=', dtim, &
                  & ' to dtim_est=', dtim_est
               write (*,*)
          end if
          dtim = dtim_est
          clean_dtim=.false.
        else
          dtim = dtim_input
        end if
      end if !nl_dtim_est
    end if !meth

    ! check the time step is not too small
    if (dtim < dt_min) then
      call gkw_warn('dt < dt_min; the run will terminate shortly')
      stop_me = .true.
    end if

    ! advance the time
    time = time + dtim

    if(call_diagnostics) call diagnostic_pre_naverage(iloop, fdisi)

    ! one more timestep done. Update some counters:
    ntotstep = ntotstep + 1
    last_smallstep_time = time
  end do time_stepping

  !If the nonlinear timestep estimator reduction is off, check across processors

  if (non_linear .and. (.not. nl_dtim_est)) then
    call mpiallreduce_min(dtim_est_save,dtim_dum,1)
    if (dtim_dum .lt. dtim) then
       call gkw_warn ('Timestep too big (nl_est)!, aborting')
       !The ideal solution is to make the code go back one large timestep
       !Until this is implemented will abort.
       stop_me=.true.
    end if
  end if


  ! at the end of the timestep, make the electromagnetic field
  ! consistent with the distribution function
  if (spectral_radius) then
    call calculate_fields(fdisi)
  else
    call field_solve_nonspec_wrap(fdisi,0.,.false.)
  endif
  ! Also the collisions conservation for neoclassical fluxes
  call calculate_cons(fdisi)

end subroutine advance_large_step_explicit

!****************************************************************************
!> A Runge Kutta second order timestep
!----------------------------------------------------------------------------
subroutine rk2
  use constants, only : c1
  use dist,      only : fdisi, fdisk, rhsk, nsolc

  integer :: i
  complex :: cdum

  ! first 'half' time step
  ! calculate the rhs
  do i = 1, nsolc
    fdisk(i,1) = fdisi(i)
  end do
  call calculate_rhs(fdisk(:,1),rhsk(:,1))

  ! advance a timestep delta*dtime, calculate f
  ! and store in fdisk(:,2)
  cdum = c1*delta
  do i = 1, nsolc
    fdisk(i,1) = fdisi(i) + cdum*rhsk(i,1)
  end do

  !  second part full time step calculated from fdisk(:,2)
  ! calculate the rhs
  call calculate_rhs(fdisk(:,1),rhsk(:,1))

  ! add the rhs to fdisi
  do i = 1, nsolc
    fdisi(i) = fdisi(i) + rhsk(i,1)
  end do

end subroutine rk2

!****************************************************************************
!> A Runge Kutta fourth order timestep
!----------------------------------------------------------------------------
subroutine rk4

  use global,    only : logical_false
  use control,   only : dtim
  use constants, only : c1
  use dist,      only : fdisi, fdisk, rhsk, nsolc
  use perform,   only : perfon, perfoff, perf_measure

  complex :: cdum
  integer :: i

  if (perf_measure) call perfon('rk4',2)

!BLAS version, double precision
#if defined(blas)

  ! copy fdisi into both buffers
  call blas_copy(nsolc, fdisi, 1, fdisk(:,1),1)
  call blas_copy(nsolc, fdisi, 1, fdisk(:,2),1)
  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1))

  ! first step into solution
  ! fdisi = fdisi + rhsk(:,1)/6.E0
  cdum=c1/6.0
  call blas_axpy(nsolc, cdum, rhsk(:,1), 1,fdisi,1)
  ! second step initialization
  ! fdisk(:,2) = fdisk(:,1) + rhsk(:,1)/2.E0
  cdum=c1/2.0
  call blas_copy(nsolc, fdisk(:,1), 1, fdisk(:,2), 1)
  call blas_axpy(nsolc, cdum, rhsk(:,1), 1,fdisk(:,2),1)
  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1),0.5*dtim)

  ! second step into solution
  ! fdisi = fdisi + rhsk(:,1)/3.E0
  cdum=c1/3.0
  call blas_axpy(nsolc, cdum, rhsk(:,1), 1,fdisi,1)
  ! third step initialization
  ! fdisk(:,2) = fdisk(:,1) + rhsk(:,1)/2.E0
  cdum=c1/2.0
  call blas_copy(nsolc, fdisk(:,1), 1, fdisk(:,2),1)
  call blas_axpy(nsolc, cdum, rhsk(:,1), 1,fdisk(:,2),1)

  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1),0.5*dtim)

  ! third step into solution
  ! fdisi = fdisi + rhsk(:,1)/3.E0
  cdum=c1/3.0
  call blas_axpy(nsolc, cdum, rhsk(:,1), 1,fdisi,1)
  ! fourth step initialization
  ! fdisk(:,2) = fdisk(:,1) + rhsk(:,1)
  cdum=c1
  call blas_copy(nsolc, fdisk(:,1), 1, fdisk(:,2),1)
  call blas_axpy(nsolc, cdum, rhsk(:,1), 1,fdisk(:,2),1)

  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1),dtim)

  ! fourth step into solution
  ! fdisi = fdisi + rhsk(:,1)/6.E0
  cdum=c1/6.0 !promote to complex
  call blas_axpy(nsolc, cdum, rhsk(:,1), 1,fdisi,1)
!NON BLAS version
!OpenMP on these loops seem to have little effect
!These operations may be memory bandwidth limited.
#else

  ! initialize to fdisk(:,1)
  !$omp parallel do schedule(static)
  do i = 1, nsolc
    fdisk(i,1) = fdisi(i)
    fdisk(i,2) = fdisi(i)
  end do
  !$omp end parallel do

  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1))

  ! first step into solution
  cdum = c1 / 6.
  !$omp parallel
  !$omp do schedule(static)
  do i = 1, nsolc
    fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
  end do
  !$omp end do

  ! second step initialization
  cdum = c1 / 2.
  !$omp do schedule(static)
  do i = 1, nsolc
    fdisk(i,2) = fdisk(i,1) + cdum*rhsk(i,1)
  end do
  !$omp end do
  !$omp end parallel

  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1),0.5*dtim)

  ! second step into solution
  cdum = c1 / 3.
  !$omp parallel
  !$omp do schedule(static)
  do i = 1, nsolc
    fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
  end do
  !$omp end do

  ! third step initialization
  cdum = c1 / 2.
  !$omp do schedule(static)
  do i = 1, nsolc
    fdisk(i,2) = fdisk(i,1) + cdum*rhsk(i,1)
  end do
  !$omp end do
  !$omp end parallel

  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1),0.5*dtim)

  ! third step into solution
  cdum = c1 / 3.
  !$omp parallel
  !$omp do schedule(static)
  do i = 1, nsolc
    fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
  end do
  !$omp end do

  ! fourth step initialization
  !$omp do schedule(static)
  do i = 1, nsolc
    fdisk(i,2) = fdisk(i,1) + rhsk(i,1)
  end do
  !$omp end do
  !$omp end parallel

  ! advance a timestep delta*dtime, calculate delta f
  call calculate_rhs(fdisk(:,2),rhsk(:,1),dtim)

  ! fourth step into solution
  cdum = c1 / 6.
  !$omp parallel do schedule(static)
  do i = 1, nsolc
    fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
  end do
  !$omp end parallel do

#endif


  ! keep compiler quiet
  if (logical_false) write (*,*) cdum

  if (perf_measure) call perfoff(2)

end subroutine rk4

!****************************************************************************
!> Runge Kutta third order (midpoint method; stable for waves) timestep. The
!> timestep estimator is not programmed for this method
!----------------------------------------------------------------------------
subroutine rk3(itime,iloop)
  use constants, only : c1
  use dist,      only : fdisi, fdisk, rhsk, nsolc

  integer, intent(in) :: itime,iloop
  integer :: i,j
  complex :: cdum,cdum2

  if (itime == 1) then

    ! To set up the scheme use fourth order Runge Kutta
    if (iloop <= 3) then

      ! initialize to fdisi
      do i = 1, nsolc
        fdisk(i,1) = fdisi(i)
        fdisk(i,4) = fdisi(i)
      end do

      ! advance a timestep delta*dtime, calculate delta f
      call calculate_rhs(fdisk(:,4),rhsk(:,1))

      ! Put in storage
      do i = 1, nsolc
        fdisk(i,4-iloop) = fdisi(i)
        rhsk(i,4-iloop)  = rhsk(i,1)
      end do

      ! no point to calculate the last time step
      if (iloop /= 3) then

        ! first step into solution
        cdum = c1 / 6.
        do i = 1, nsolc
          fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
        end do

        ! second step initialization
        cdum = c1 / 2.
        do i = 1, nsolc
          fdisk(i,4) = fdisk(i,1) + cdum*rhsk(i,1)
        end do

        ! advance a timestep delta*dtime, calculate delta f
        call calculate_rhs(fdisk(:,4),rhsk(:,1))

        ! second step into solution
        cdum = c1 / 3.
        do i = 1, nsolc
          fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
        end do

        ! third step initialization
        cdum = c1 / 2.
        do i = 1, nsolc
          fdisk(i,4) = fdisk(i,1) + cdum*rhsk(i,1)
        end do

        ! advance a timestep delta*dtime, calculate delta f
        call calculate_rhs(fdisk(:,4),rhsk(:,1))

        ! third step into solution
        cdum = c1 / 3.
        do i = 1, nsolc
          fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
        end do

        ! fourth step initialization
        do i = 1, nsolc
          fdisk(i,4) = fdisk(i,1) + rhsk(i,1)
        end do

        ! advance a timestep delta*dtime, calculate delta f
        call calculate_rhs(fdisk(:,4),rhsk(:,1))

        ! fourth step into solution
        cdum = c1 / 6.
        do i = 1, nsolc
          fdisi(i) = fdisi(i) + cdum*rhsk(i,1)
        end do

      end if

    end if

  end if

  if (itime == 1 .and. iloop < 3) return

  ! initialize
  do i = 1, nsolc
    fdisk(i,4) = (0.,0.)
  end do

  ! update the solution
  do j = 1, 3
    cdum = c1*alf(j) ; cdum2 = c1*bet(j)
    do i = 1, nsolc
      fdisk(i,4) = fdisk(i,4) + cdum*fdisk(i,j) &
                & + cdum2*rhsk(i,j)
    end do
  end do
  cdum = c1 / gam
  do i = 1, nsolc
    fdisk(i,4) = cdum*fdisk(i,4)
  end do

  ! resuffle (could be made more efficient)
  do i = 1, nsolc
    fdisk(i,3) = fdisk(i,2)
    fdisk(i,2) = fdisk(i,1)
    fdisk(i,1) = fdisk(i,4)
    fdisi(i) = fdisk(i,4)
    rhsk(i,3) = rhsk(i,2)
    rhsk(i,2) = rhsk(i,1)
  end do

  ! calculate a new right hand side
  call calculate_rhs(fdisk(:,4),rhsk(:,1))

end subroutine rk3

!*****************************************************************************
!> An experimental Runge-Kutta-Chebyshev second-order arbitrary stage timestep.
!> Requires input s, the number of stages. Largely based on the subroutine
!> 'step' of http://www.netlib.org/ode/rkc.f (see other reference below).
!> As is, the scheme has been shown to allow larger timesteps where diffusive
!> terms (e.g. collisions) dominate through improved stability. A primary
!> reference is
!>
!> B. P. Sommeijer, L. F. Shampine and J. G. Verwer (1997),
!> RKC: An explicit solver for parabolic PDEs, J. Comp. Appl. Math. 88, 315--326
!> J. G. Verwer (1980), Explicit Runge-Kutta methods for parabolic partial
!> differential equations, Appl. Numer. Math. 22, 359--379
!>
!>
!> If a diffusive term is dominating the timestep estimate by more than a
!> factor of two, then this scheme will be beneficial over standard RK4.
!> One can also adapt the scheme in use in real time, but it is not clear how
!> useful this is for delta-f gyrokinetics -- new variants of this scheme and
!> similar schemes may exist that are more fruitful.
!-----------------------------------------------------------------------------
subroutine rkc_2(s)
  use control,   only : dtim
  use global,    only : dp
  use constants, only : c1
  use dist,      only : fdisi, fdisk, rhsk, nsolc
  use perform,   only : perfon, perfoff, perf_measure

  !> number of stages
  integer, intent(in) :: s

  !> Parameter which might in principle be varied according to the size of
  !> advection/diffusion terms or the timestep control method.
  real (dp), parameter :: deps = 2._dp/13._dp

  ! Use type double in calculations, then convert to complex when used in
  ! conjunction with fdisi, fdisk or rhsh.
  real (dp) :: a1,a2,a3,w0,w1,bj,b1,b2,tf,tf1,tf2
  real (dp) :: tj,t1,t2,tjp,t1p,t2p,tjpp,t1pp,t2pp
  real (dp) :: mu,nu,mus,gt,mm
  real :: dtdum
  complex   :: h1,h2,h3,h4,h5,cdum
  integer   :: i,j,i1,i2,ii

  if (perf_measure) call perfon('rkc_arb2',2)

  ! Could here precalculate all Chebyshev polynomials and derivatives exactly,
  ! then convert to desired datatype, but it is relatively inexpensive to
  ! calculate and re-calculate them recursively as needed.
  w0 = 1._dp  + deps / real(s**2,dp)
  a1 = w0**2 - 1._dp
  a2 = sqrt(a1)
  a3 = real(s,dp)*log(w0 + a2)
  ! w1 = T_s'(w0) / T_s''(w0)
  w1 = sinh(a3)*a1 / (cosh(a3)*real(s,dp)*a2 - w0*sinh(a3))

  ! intial values
  b1 = 0.25_dp / w0**2
  b2 = b1
  bj = 0._dp    ! b_j       [calculate later]
  tj = 0._dp    ! T_{j}(w0) [calculate later]
  t1 = w0       ! T_{j-1}(w0)
  t2 = 1._dp    ! T_{j-2}(w0)
  tf1 = w1*b1   ! sub-timestep
  tf2 = 0._dp   ! previous sub-timestep
  tjp = 1._dp   ! T'_j(w0)
  t1p = 1._dp   ! T'_{j-1}(w0)
  t2p = 0._dp   ! T'_{j-2}(w0)
  tjpp = 0._dp  ! T''_j(w0)
  t1pp = 0._dp  ! T''_{j-1}(w0)
  t2pp = 0._dp  ! T''_{j-2}(w0)

  ! initialize to fdisi
  do i=1, 3
    fdisk(:,i) = fdisi
  end do

  ! calculate delta f at t
  call calculate_rhs(fdisk(:,2),rhsk(:,1))

  ! store [f_1 = f_n(t,f_n) + dtim * tf1 * df_n(t,f_n) ]
  cdum = c1*tf1
  do i=1, nsolc
    fdisk(i,1) = fdisk(i,3) + cdum*rhsk(i,1)
  end do

  ! first step into solution fdisi; f_n = f_1
  fdisi = fdisk(:,1)

  ! array referer initial values
  i2 = 2 ; i1 = 1 ; ii = 0

  do j=2, s

    ! swap the array references
    ii = i1 ; i1 = i2 ; i2 = ii

    ! Evaluate Chebyshev polynomials and derivatives; caculate coefficients,
    ! convert to complex.
    tj   = 2._dp*w0*t1   -       t2
    tjp  = 2._dp*w0*t1p  + 2._dp*t1   - t2p
    tjpp = 2._dp*w0*t1pp + 4._dp*t1p  - t2pp
    bj   = tjpp / tjp**2
    mu   = 2._dp*bj*w0 / b1    ; h2 = c1*mu
    nu   =      -bj    / b2    ; h1 = c1*nu
    mus  = mu*w1 / w0          ; h4 = c1*mus
    gt   = (b1*t1 - 1._dp)*mus ; h5 = c1*gt
    mm   = 1._dp - mu - nu     ; h3 = c1*mm

    ! calculate delta f_{j-1} at t + tf_{j-1}*dtim
    ! prevent type mismatch with calculate_rhs
    dtdum = tf1*dtim
    call calculate_rhs(fdisi,rhsk(:,2),dtdum)

    do i=1, nsolc
      fdisk(i,i1) = h1 * fdisk(i,i1) &
                & + h2 * fdisk(i,i2) &
                & + h3 * fdisk(i,3)  &
                & + h4 * rhsk(i,2)   &
                & + h5 * rhsk(i,1)
    end do

    ! update fdisi
    fdisi = fdisk(:,i1)

    ! sub-timestep (used for next delta f_{j-1})
    tf  = mu*tf1 + nu*tf2 + mus + gt
    ! re-order for next evaluations
    tf2 = tf1
    tf1 = tf
    b2 = b1
    b1 = bj
    t2pp = t1pp
    t1pp = tjpp
    t2p = t1p
    t1p = tjp
    t2 = t1
    t1 = tj

  end do

  if (perf_measure) call perfoff(2)

end subroutine rkc_2


!--------------------------------------------------------------------
!>  Subroutine that calculates the right hand side of the gyrokinetic equation.
!>  Optimised routine, contains OpenMP
!--------------------------------------------------------------------
subroutine calculate_rhs(fdis,rhs,DPART_IN) ! Optimised
use control,    only : nlapar,non_linear,dtim
use grid,       only : lsendrecv_x
use matdat,     only : mat, mat_maxwll_background, matg2f, add_source
use rotation,   only : shear_real
use dist,       only : nsolc, nf, fdis_tmp, fdis_tmp2
use non_linear_terms, only : add_non_linear_terms
use control,    only : spectral_radius
use mpighosts,  only : gc1, gc1x2, gc2x2_f, gc2x2_pga, gckrook, mpistart, mpiwait
use constants,  only : c1
use dist,       only : iphi_ga, ifdis, iapar_ga
use rotation,   only : shear_periodic_bcs
use fields,     only : calculate_cons, calculate_fields, field_solve_nonspec_wrap
use krook,      only : nlkrook, add_krook
use source_time,only : add_source_time
use perform, only : perfon, perfoff, perffields, perf_measure, perfswitch

use general, only : gkw_abort
use matrix_format, only : usmv
  !> The distribution function and all fields. This is input for this
  !> routine and the values are destroyed afterwards.
  complex, intent(inout) :: fdis(nsolc)
  !>the right hand side. An explict step can be calculated through
  !>fdis = fdis + rhs, the time interval value of the time step is
  !>already included in rhs.
  complex, intent(out) :: rhs(nsolc)
  !> This is the sub time step dt, which is used for anthing in the
  !> RHS which depends on time.
  real, intent(in), optional :: DPART_IN

  complex :: deltatime
  integer :: i, ierr
  real :: dpart

  if (present(DPART_IN)) then
    dpart = DPART_IN
  else
    dpart = 0.
  end if

  ierr = 0

  ! copy mode that should be temporally persistent into fdis
  call set_persistent_mode(fdis)

  if (perf_measure) then
   call perfon('Calc RHS',2)
  end if

  if (nlkrook) then
    ! the calculate_rhs routine is invokled for several different
    ! buffers, namely fdisi or the fdisk columns. The ghost cell
    ! communication must already know which array it will be working
    ! on when it is set up. For this reason, we copy:
    fdis_tmp(1:nf) = fdis(1:nf)
    ! start velocity mirror communications for krook operator from
    ! fdis_tmp into fdis_nvpar
    call mpistart(gckrook)
  end if

  if(perf_measure) then
    call perfon('Calc Fields',2)
    perffields=.true.
  end if

  ! First calculate the electro-static adnd electro-magnetic fields,
  ! from fdis into fdis itself
  if (spectral_radius) then
    call calculate_fields(fdis,DPART=dpart)
  else
    call field_solve_nonspec_wrap(fdis,dpart,non_linear)
  endif

  ! if using radial communications, fdis_tmp2 = fdis already.  now
  ! start communications of gyroaveraged phi ghost points (only) into
  ! fdis_tmp2.  fdisi part was already started before the field solve
  ! by field_solve_nonspec_wrap.
  if (non_linear) call mpistart(gc2x2_pga)

  if (perf_measure) then
    call perfswitch('Copy fdis vector to tmp',2)
    perffields=.false.
  end if

  ! Copy fdis into fdis_tmp. For electro-magnetic runs, undo the A||
  ! correction of g at the same time.
  if (nlapar) then
    !$omp parallel
    !$omp do schedule(static)
    do i = 1, matg2f%nmat
      fdis_tmp(i) = fdis(i) + matg2f%mat(i)*fdis(matg2f%jj(i))
    end do
    !$omp end do
    !$omp do schedule(static)
    do i= matg2f%nmat+1, nsolc
      fdis_tmp(i) = fdis(i)
    end do
    !$omp end do
    !$omp end parallel
  else
#if defined(blas)
    call blas_copy(nsolc, fdis, 1, fdis_tmp,1)
#else
    !$omp parallel do schedule(static)
    do i= 1, nsolc
      fdis_tmp(i) = fdis(i)
    end do
    !$omp end parallel do
#endif
  end if

  ! Send/Recv the distribution function to/from neighbours into fdis_tmp
  ! (vpar, s and mu ghost points only)
  call mpistart(gc1)
  ! Send/Recv the distribution function (x ghost points only) into fdis_tmp
  call mpistart(gc1x2)

  ! Calculate the collisions conservation 'fields' after g2f transform
  ! APS suggests this can be overlapped with the communication
  ! Since no part of this field goes into ghost cells.
  if (perf_measure) call perfswitch('Collisions conservation',2)
  call calculate_cons(fdis_tmp)
  if (perf_measure) call perfswitch('Init RHS to zero',2)

  ! initialise the RHS to zero
  !$omp parallel do schedule(static)
  do i = 1, nsolc
    rhs(i) = (0.,0.)
  end do
  !$omp end parallel do

  if (perf_measure) call perfswitch('MPI radial ghost comms',2)

  ! Wait for the radial derivatives communication (only) to finish into fdis_tmp2.
  ! Both the fields and the distribution are required before nonlinear terms
  if (non_linear) call mpiwait(gc2x2_f)
  if (non_linear) call mpiwait(gc2x2_pga)

  ! apply shear periodic boundary conditions when needed
  call shear_periodic_bcs(fdis_tmp2,ifdis,dpart)
  call shear_periodic_bcs(fdis_tmp2,iphi_ga,dpart)
  call shear_periodic_bcs(fdis_tmp2,iapar_ga,dpart)

  ! The krook operator
  if (nlkrook) then
    if (perf_measure) call perfswitch('Krook operator:   MPI mirror',2)
    call mpiwait(gckrook)
    if (perf_measure) call perfswitch('Krook op: Local + MPI reduce',2)
    
    if (lsendrecv_x) then
      call add_krook(fdis_tmp2,rhs)
    else
      call add_krook(fdis,rhs)
    endif
  endif

  if (perf_measure) call perfswitch('Non linear terms: FFT, No MPI',2)

  ! Call the nonlinear terms routine if necessary.
  ! The nonlinear terms are called with fdis (or fdis_tmp2), which contains the
  ! distribution g = f + Z v\\ A\\ etc., rather than f, which is presently
  ! in fdis_tmp.
  if (non_linear .or. shear_real) then
    if (lsendrecv_x) then
      call add_non_linear_terms(fdis_tmp2,rhs)
    else
      call add_non_linear_terms(fdis,rhs)
    end if
  end if


  ! Call the routine which adds the time dependent source term
  call add_source(rhs, dtim)
  call add_source_time(rhs, dpart, dtim)

  if (perf_measure) call perfswitch('MPI other ghost comms',2)

 ! If necessary, wait for the s vpar mu + x ghosts communication to finish.
 ! fdis and fields all required here
 call mpiwait(gc1)
 call mpiwait(gc1x2)

 ! apply shear periodic boundary conditions when needed
 call shear_periodic_bcs(fdis_tmp,ifdis,dpart)
 call shear_periodic_bcs(fdis_tmp,iphi_ga,dpart)
 call shear_periodic_bcs(fdis_tmp,iapar_ga,dpart)

!APS: In order to have the conserving parts of collisions calculated
!APS: correctly when using parallel velocity grids, the calculation of those
!APS: `fields' must be performed after the communication has ended, rather
!APS: than in the fields calculation subroutine.
! FJC: Does this comment still apply?

  if (perf_measure) call perfswitch('Linear terms matmul, No MPI',2)

  deltatime=dtim*c1
  call usmv(deltatime,mat,fdis_tmp,rhs,ierr)
  call usmv(deltatime,mat_maxwll_background,fdis_tmp,rhs,ierr)
  if(ierr /= 0) call gkw_abort('error in usmv')

  if (perf_measure) call perfoff(2) ! matmul
  if (perf_measure) call perfoff(2) ! calc RH

end subroutine calculate_rhs


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>Deallocate arrays for the module
!-----------------------------------------------------------------------
subroutine exp_integration_deallocate
  use dist, only : fdisk, rhsk

  if (allocated(fdisk))  deallocate(fdisk)
  if (allocated(rhsk))   deallocate(rhsk)

end subroutine exp_integration_deallocate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> force a mode to a constant value.
!-----------------------------------------------------------------------
subroutine set_persistent_mode(fdis)

  use components, only : mode_persist, ind_persist
  use grid,    only : nx, ns, nmu, nvpar, nsp

  complex, intent(inout) :: fdis(:)

  integer :: i

  if (mode_persist) then
    do i = 1, nx*ns*nmu*nvpar*nsp
      fdis(ind_persist(i)) = persistent_mode(i)
    end do
  end if

end subroutine set_persistent_mode


end module exp_integration
