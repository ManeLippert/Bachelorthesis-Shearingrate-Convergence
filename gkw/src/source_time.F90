!-------------------------------------------------------------------------------
!> Provides the modulated source for numerical experiments on turbulence
!> as used in publications by Miglano et al.
!-------------------------------------------------------------------------------
module source_time

  use global, only : lenswitch
#ifdef blas
#ifdef real_precision_default
  external :: blas_copy => ccopy
  external :: blas_axpy => caxpy
#else
  external :: blas_copy => zcopy
  external :: blas_axpy => zaxpy
#endif
#endif
  use global, only : lenswitch
  
  implicit none 

  private 
  
  public :: source_time_allocate
  public :: source_modulation, source_time_read_nml, source_time_write_nml
  public :: source_time_check_params, source_time_bcast_nml
  public :: add_source_time, source_time_ampl

  !------------------------------------------------------------------------------
  !parameter which determine cosine heating profile
  !------------------------------------------------------------------------------

  !> deprecated, do not use
  logical, save :: modulation

  !> source modulation wave length
  integer, save :: source_wave_number

  !> amplitude of the time dependent source
  complex, save :: source_time_ampl

  !> modulation frequency
  real, save :: mod_freq
  
  !------------------------------------------------------------------------------
  !parameter which determine gaussian source sink heating profile
  !------------------------------------------------------------------------------
  
  !>broadness of gaussian
  real, save :: gauss_source_stdev

  !>position of maximum of gaussian (distance from inner radial boundary)
  integer, save :: gauss_source_median
  
  !>begin of source free region (distance from inner radial boundary)
  integer, save :: dsfr

  
  !> Complex array that contains the source in the equation for the evolution
  !> of the distribution function. This source is per definition all the
  !> terms that are independent of either the distribution function or any
  !> of the fields.
  complex, save, allocatable :: source_t(:)
  
  !> source profile (cosine, gauss_source_sink, )
  character(len=lenswitch), save :: source_profile

  
  interface source_time_write_nml
     module procedure source_time_read_nml
  end interface

contains

  !-----------------------------------------------------------------------------
  !>
  !-----------------------------------------------------------------------------
  subroutine source_time_allocate()
    use dist, only : nsolc
    use general, only : gkw_abort
    integer :: ierr, i

    allocate(source_t(nsolc),stat=ierr)
    if (ierr /= 0) call gkw_abort('matdat_allocate: cannot allocate source_t')

    do i = 1, nsolc
      source_t(i) = (0.,0.)
    end do
  end subroutine source_time_allocate

!--------------------------------------------------------------------------------
!> This routine produces the matrix element for a time dependent source.
!> Possible profiles are
!> 
!> cosine:
!>
!>      S=amp*(v^2/v_th^2-3/2)*F_M*exp(i*(x/L-omega*t)) L=lx/source_wave_number
!>
!> gauss: 
!>
!>      S=amp*(v^2/v_th^2-3/2)*F_M* &
!>        & *(exp(-((gx(ix)-gauss_source_median)*dxgr)**2 / (2.0 * gauss_source_stdev**2.0)) &
!>        & -exp(-((gx(ix)-1-(n_x_grid-gauss_source_median))*dxgr)**2 / (2.0 * gauss_source_stdev**2.0))
!>
!> gauss_sfr:
!>
!>
!>                                   | * exp(-((gx(ix)-gauss_source_median)*dxgr)**2 / 
!>                                   | / (2.0 * gauss_source_stdev**2.0)) 
!>                                   |   if gx(ix) <= dsfr
!>                                   |
!>      S=amp*(v^2/v_th^2-3/2)*F_M * | * 0
!>                                   |   if dsfr < gx(ix) <= n_x_grid-dsfr
!>                                   |
!>                                   | * -exp(-((gx(ix)-1-(n_x_grid-gauss_source_median))*dxgr)**2 /
!>                                   | /(2.0 * gauss_source_stdev**2.0))
!>                                   |   if n_x_grid - dsfr < gx(ix)
!>
!>
!> The source term is defined such that when added to the rhs of the gyrokinetic 
!> equation it does not produce particle and parallel momentum sources. 
!> For this purpose intfm (normalized velocity space integral of the Maxwellian: 
!> int(v^2/v_th^2*F_M)*dv ~ 3/2 ) is calculated and used in the matrix element 
!> instead of the 3/2 that appears in the formula above. The source term is 
!> then stored.
!> All parameters (length, amplitude, ecc..) have to be given in the namelist  
!> source_time. The routine is called in linear_terms.f90.  
!--------------------------------------------------------------------------------
subroutine source_modulation
  
  use grid,           only : nx, ns, nmu, nvpar, nsp, lx, n_x_grid
  use grid,           only : parallel_vpar, parallel_mu, gx 
  use dist,           only : fmaxwl, ifdis
  use index_function, only : indx
  use geom,           only : bn, dxgr
  use velocitygrid,   only : vpgr, mugr, intvp, intmu
  use components,     only : tmp, tgrid
  use general,        only : gkw_abort
  use mode,           only : iyzero
  use mpiinterface,   only : mpiallreduce_sum_inplace
  use mpicomms,       only : COMM_VPAR_NE_MU_NE
  use constants,      only : ci1, pi
  use global, only : r_tiny

  ! real for the normalized velocity space integral of the Maxwellian
  ! the integral is calculated at each grid points (species, radial, field)
  real :: intfm, intfm1, intfm2

  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is

  ! matrix element and amplitude of the source
  complex :: mat_elem 

  ! dummy variables
  integer :: iih

  ! normalization factor for gauss_source
  real :: gauss_source_norm
  
  gauss_source_norm = 1.0 / sqrt( 2.0 * pi * gauss_source_stdev**2.0 )
  
  if (abs(source_time_ampl) < r_tiny) return

  ! check if there is a zero mode 
  if (iyzero == 0) call gkw_abort('No zero mode for modulation')
    
  do ix=1,nx; do i=1,ns; do is=1,nsp
    
    ! initialize the integral to zero
    intfm1 = 0.0
    intfm2 = 0.0 

    ! perform the integral in velocity space
    do j=1,nmu ; do k=1,nvpar

      intfm1 = intfm1 &
               & +(vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j))/(tmp(ix,is)/tgrid(is))  &
               & *fmaxwl(ix,i,j,k,is) &
               & *intvp(i,j,k,is)*intmu(j)*bn(ix,i) 

      intfm2 = intfm2 + fmaxwl(ix,i,j,k,is)*intvp(i,j,k,is)*intmu(j)*bn(ix,i) 

    end do ; end do
      
    ! MPI allreduce over COMM_VPAR_NE_MU_NE 
    if (parallel_vpar.or.parallel_mu) then 
      call mpiallreduce_sum_inplace(intfm1,1,COMM_VPAR_NE_MU_NE)
      call mpiallreduce_sum_inplace(intfm2,1,COMM_VPAR_NE_MU_NE)
    endif 

    ! normalize the sum of the energy moment to the sum of the Maxwell such that 
    ! the particle source is exactly zero 
    intfm = intfm1 / intfm2 

    ! calculate and store the matrix element
    do j=1,nmu ; do k=1,nvpar

      ! select the zonal mode
      imod = iyzero
 
      ! index of the perturbed distribution
      iih = indx(ifdis,imod,ix,i,j,k,is) 
 
      ! matrix element 
      mat_elem = ((vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j))/(tmp(ix,is)/tgrid(is))-intfm)  &
               &  *source_time_ampl*fmaxwl(ix,i,j,k,is)
               
      ! select source profile
      ! 'cosine' -> cosine shaped
      ! 'gauss_source_sink' -> gauss shaped source and sink
      ! 'default' -> cosine shaped
      select case(source_profile)
      
            case('cosine')
              mat_elem = mat_elem*exp(2.E0*pi*ci1*source_wave_number*(gx(ix)-1)*dxgr/lx)
              
            case('gauss')
              mat_elem = mat_elem*(exp(-1.0 * ((gx(ix)-gauss_source_median)*dxgr)**2 &
              & / (2.0 * gauss_source_stdev**2.0)) &
              & -exp(-1.0*((gx(ix)-1-(n_x_grid-gauss_source_median))*dxgr)**2) &
              & / (2.0 * gauss_source_stdev**2.0)) * gauss_source_norm
              
            case('gauss_sfr')
              ! source region (inner edge): S = positive gaussian
              if(gx(ix) <= dsfr) then
                mat_elem = mat_elem*(exp(-1.0 * ((gx(ix)-gauss_source_median)*dxgr)**2 &
                & / (2.0 * gauss_source_stdev**2.0)) * gauss_source_norm )
              ! source free region: S = 0
              else 
                if (gx(ix) > dsfr .and. gx(ix) <= n_x_grid-dsfr) then
                  mat_elem = mat_elem * 0.0
                ! sink region (outer edge): S = negative gaussian
                else 
                  mat_elem = mat_elem*(-exp(-1.0 * ((gx(ix)-1-(n_x_grid-gauss_source_median))*dxgr)**2 &
                  & / (2.0 * gauss_source_stdev**2.0) ) * gauss_source_norm)
                endif
              end if
              
            case default
              mat_elem = mat_elem*exp(2.E0*pi*ci1*source_wave_number*(gx(ix)-1)*dxgr/lx)
              
      end select
      ! put the element 
      source_t(iih) = source_t(iih) + mat_elem

    end do ; end do 

  end do  ; end do ; end do

end subroutine source_modulation

!------------------------------------------------------------------------------
!> This subroutine adds the time dependent source term
!> of the gyrokinetic equation
!------------------------------------------------------------------------------
subroutine add_source_time(rhs, DPART_IN, deltatime)
  use constants, only : ci1
  use control,   only : time
  use dist,      only : nsolc
  complex, intent(inout) :: rhs(nsolc)
  real, intent(in) :: DPART_IN
  real, intent(in) :: deltatime

  integer :: i
  complex :: cmodul
  if (abs(source_time_ampl) > 0.0) then

    cmodul = exp(-ci1*mod_freq*(time+DPART_IN))
    !$omp parallel private(i)
    !$omp do
    do i = 1, nsolc
      !rhs(i) = rhs(i) + 0.5*deltatime*(cmodul*source_t(i) + conjg(cmodul*source_t(i)))
      rhs(i) = rhs(i) + deltatime*real(cmodul*source_t(i))
    end do
    !$omp end do
    !$omp end parallel
  end if
end subroutine add_source_time

!------------------------------------------------------------------------------
!> This subroutine reads (or writes) the time dependent source namelist
!------------------------------------------------------------------------------
subroutine source_time_read_nml(lun,io_stat,lwrite)
  use io, only : write_run_parameter
  use constants, only : c1
  use mpiinterface, only : root_processor
  use grid, only : n_x_grid

  integer, intent(in)  :: lun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  namelist /source_time/ modulation, source_wave_number, &
         & source_time_ampl, mod_freq, source_profile, gauss_source_stdev, &
         & gauss_source_median, dsfr

  io_stat = 0
  ! read the input
  if (present(lwrite)) then
    if (.not. lwrite) then
      ! set the default values 
      source_wave_number = 1
      source_time_ampl = 0.0
      mod_freq = 1.0
      source_profile = 'cosine'
      gauss_source_stdev = 1.8
      dsfr = int(n_x_grid/4)
      gauss_source_median = int(n_x_grid/8)

      !deprecated, do not use:
      modulation = .true.

      ! read namelist
      read(lun,NML=source_time,IOSTAT=io_stat) 
    end if
  else
    ! write the namelist
    if(root_processor) write(lun,NML=source_time)
    call write_run_parameter('source_time', 'source_wave_number', source_wave_number)
    call write_run_parameter('source_time', 'source_time_ampl', source_time_ampl)
    call write_run_parameter('source_time', 'mod_freq', mod_freq)
    call write_run_parameter('source_time', 'source_profile', source_profile)
    call write_run_parameter('source_time', 'gauss_source_stdev', gauss_source_stdev)
    call write_run_parameter('source_time', 'gauss_source_median', gauss_source_median)
    call write_run_parameter('source_time', 'dsfr', dsfr)
  end if

end subroutine source_time_read_nml


!------------------------------------------------------------------------------
!> bcast the time dependent source namelist params
!------------------------------------------------------------------------------
subroutine source_time_bcast_nml
  use mpiinterface, only : mpibcast
  call mpibcast(source_wave_number, 1)
  call mpibcast(source_time_ampl,   1)
  call mpibcast(mod_freq,           1)
  call mpibcast(source_profile,     lenswitch)
  call mpibcast(gauss_source_stdev, 1)
  call mpibcast(gauss_source_median,1)
  call mpibcast(dsfr,               1)
  !deprecated:
  call mpibcast(modulation,         1)
end subroutine source_time_bcast_nml

!------------------------------------------------------------------------------
!> put any checks that can be done before memory allocation in here
!------------------------------------------------------------------------------
subroutine source_time_check_params
  use grid, only : n_x_grid
  use control,    only : spectral_radius
  use general, only : gkw_abort, gkw_warn

  if(.not.modulation) then
    call gkw_warn('The modulation paramater in the SOURCE TIME namelist has been&
       & deprecated. Please do not use it.')
  end if
  if(abs(source_time_ampl) > 0.0) then
    if (spectral_radius) then
      call gkw_abort('Must be running non spectral to run modulation')
    end if
    if (gauss_source_median > n_x_grid) then
      call gkw_abort('gauss_source_median must not be larger than n_x_grid')
    end if
    if (2 * dsfr > n_x_grid) then
      call gkw_abort(' The total source free region (2 * dsfr) must not&
         & be larger than n_x_grid')
    end if
  end if


end subroutine source_time_check_params

end module source_time
