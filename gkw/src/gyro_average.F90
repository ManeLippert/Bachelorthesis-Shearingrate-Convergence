!------------------------------------------------------------------------------
!> This module contains the data, functions and subroutines to perform the 
!> gyro-average of the distribution function as well as the fields. 
!> at present also the field equations of the nonspectral case are part of 
!> this module 
!------------------------------------------------------------------------------
module gyro_average
 
  use global,          only : iumf,rumf
  use matrix_format, only : sparse_matrix

  implicit none 

  ! public subroutines routines 
  public :: gyro_average_read_nml, gyro_average_bcast_nml, gyro_average_write_nml
  public :: gyro_average_init, gyro_average_allocate, polarization_init
  public :: write_fdis

  ! public variables used by fields module (could be moved and calculated there)
  public :: matred, matav
  public :: matint, iiscatz, iilocz, nscatz
  public :: control_umf, info, ihavemodes, n_phi, nmod_l
  public :: phi_ga_start, phi_ga_end, apar_ga_start, apar_ga_end
  ! variables for the gyro average in the s direction
  public :: matred_s, matpol, n_s_ga
  public :: first_imod, last_imod, blending_order

  private 

  !> the number of points in the gyro-angle integral
  integer, save :: n_points_ga

  !> Additional number of points for gyroaverage array
  integer, parameter :: n_points_ga_additional = 4

  !> do the polarisation as in orb
  logical, save :: orb_polarize 

  !> number of indices in index_gyro_av 
  integer, allocatable, save :: nind_gyro_av(:,:,:,:,:) 

  !> the ix index array for the gyro-average 
  integer, allocatable, save :: index_gyro_av(:,:,:,:,:,:) 

  !> the coefficients for the gyro-average
  complex, allocatable, save :: coeff_gyro_av(:,:,:,:,:,:) 

  !> number of indices in index_gyro_av in the s direction 
  integer, allocatable, save :: nind_gyro_av_s(:,:,:,:,:) 

  !> the ix index array for the gyro-average  in the s direction
  integer, allocatable, save :: index_ix_gyro_av_s(:,:,:,:,:,:) 

  !> the i index array for the gyro-average in the s direction 
  integer, allocatable, save :: index_i_gyro_av_s(:,:,:,:,:,:) 

  !> the coefficients for the gyro-average in the s direction
  complex, allocatable, save :: coeff_gyro_av_s(:,:,:,:,:,:) 

  !> array for the gyro-average of the field solve 
  type(sparse_matrix) :: matav

  !> s component of polarisation (global in nmod and n_x_grid)
  type(sparse_matrix) :: matpol

  !> matrix for the integral of the distribution 
  type(sparse_matrix) :: matint

  type(sparse_matrix) :: matred
 
  !> parameters needed by UMFPACK
  real(rumf), save :: control_umf(20)
  real(rumf), save :: info (90)

  !> pointers to matrices stored in UMFPACK  
  integer (iumf), save, public :: numeric, symbolic

  !> help integers to set the gyro-averaged fields to zero 
  integer, save :: phi_ga_start, phi_ga_end, apar_ga_start, apar_ga_end

  !> Maxwell global in x 
  real,    allocatable, save :: fmaxwl_Gx(:,:,:,:,:)

  !> Energetic particles dist global in x
  real,    allocatable, save :: fEP_Gx(:,:,:,:,:)

  !> magnetic field global in x 
  real,    allocatable, save :: bn_Gx(:,:)

  !> coefficients of the gyro-average global in x 
  complex, allocatable, save :: coeff_gyro_av_Gx(:,:,:,:,:,:)

  !> index array for gyro average, global in x 
  integer, allocatable, save :: index_gyro_av_Gx(:,:,:,:,:,:)

  !> number of indices in the index array for gyro-average, global in x 
  integer, allocatable, save :: nind_gyro_av_Gx(:,:,:,:,:) 

  !> number of indices in index_gyro_av in the s direction 
  integer, allocatable, save :: nind_gyro_av_s_gx(:,:,:,:,:) 

  !> the ix index array for the gyro-average  in the s direction
  integer, allocatable, save :: index_ix_gyro_av_s_gx(:,:,:,:,:,:) 

  !> the i index array for the gyro-average in the s direction 
  integer, allocatable, save :: index_i_gyro_av_s_gx(:,:,:,:,:,:) 

  !> the coefficients for the gyro-average in the s direction
  complex, allocatable, save :: coeff_gyro_av_s_gx(:,:,:,:,:,:) 

  !> number of indices in index_gyro_cj in the s direction 
  integer, allocatable, save :: nind_gyro_cj_s_gx(:,:,:,:,:) 

  !> the ix index array for the gyro-average  in the s direction
  integer, allocatable, save :: index_ix_gyro_cj_s_gx(:,:,:,:,:,:) 

  !> the i index array for the gyro-average in the s direction 
  integer, allocatable, save :: index_i_gyro_cj_s_gx(:,:,:,:,:,:) 

  !> the coefficients for the conjugated gyro-average in the s direction
  complex, allocatable, save :: coeff_gyro_cj_s_gx(:,:,:,:,:,:) 

  type(sparse_matrix), save :: matred_s

  !> Switch for the use of the conjugated operator 
  logical, save :: use_conj  

  !> Switch for the consistent long wavelength treatment of the gyro-average
  logical, save :: consistent_long_wave
  
  !> Switch to enable / disable parallelism of field solve over toroidal modes
  !> Disable if nonspectral + parallel_mu/sp/vpar + nmod > 1 gives problems
  logical, save :: parallel_mod

  !> ccoefficients of the gyro-average conjugate global in x 
  complex, allocatable, save :: coeff_gyro_cj_Gx(:,:,:,:,:,:)

  !> index array for gyro average cojugate, global in x 
  integer, allocatable, save :: index_gyro_cj_Gx(:,:,:,:,:,:)

  !> number of indices in the index array for gyro-average conjugate, global 
  ! in x 
  integer, allocatable, save :: nind_gyro_cj_Gx(:,:,:,:,:) 

  !> gather_array switch
  logical, save :: ALL_PROCS=.true. 

  !> First and last modes to do on this processor, and number locally
  integer, save :: first_imod, last_imod, nmod_l

  !> Indices for rearranging toroidal modes
  integer, save :: nscatz, n_phi

  !> arrays for rearranging the fields global in toroidal mode
  integer, allocatable, save :: iiscatz(:), iilocz(:)

  !> If this processor has any work to do in parallelisation
  !> of implicit polarisation solve over modes
  logical, save :: ihavemodes = .true.

  !> Logical that determines whether the ions are gyro-averaged 
  logical, save :: gyro_average_ions

  !> Logical that determines whether the electrons are gyro-averaged 
  logical, save :: gyro_average_electrons

  !> Logical that if true enforces the gyro-average operator to be Hermitian 
  logical, save :: mk_gyro_av_hermitian 

  !> number of grid points in s-direction used for gyro average
  !> For s_average = .false. this is not used
  integer, save :: n_s_ga

  !> Integer that determines the order of the blending scheme 
  integer, save :: blending_order

  !> Integer that measures the actual size of n_gav_bound (the number of
  !> radial grid points beyond the boundary in gyroaverage). In
  !> shaped geometries the estimate of this value sometimes fails.
  integer, save :: n_gav_bound_measure = 0 
  
  interface gyro_average_write_nml
    module procedure gyro_average_read_nml
  end interface


contains 

!------------------------------------------------------------------------------
!> This subroutine reads (or writes) the gyro-average namelist
!------------------------------------------------------------------------------

subroutine gyro_average_read_nml(lun,io_stat,lwrite)
  use io, only : write_run_parameter
  use dist, only : n_gav_bound_ex
  use mpiinterface, only : root_processor

  integer, intent(in)  :: lun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  namelist /gyroaverage / n_points_ga, gyro_average_ions, gyro_average_electrons, &
                        & orb_polarize, mk_gyro_av_hermitian, use_conj, &
                        & consistent_long_wave, parallel_mod, blending_order, n_gav_bound_ex
 
  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then 
      ! set the default values 
      n_points_ga            = 32 
      gyro_average_ions      = .true.
      gyro_average_electrons = .true.
      orb_polarize           = .false. 
      mk_gyro_av_hermitian   = .false. 
      use_conj               = .false. 
      consistent_long_wave   = .false.
      parallel_mod           = .true. 
      blending_order         =  2 
      n_gav_bound_ex         = 0
           
      read(lun,NML=gyroaverage,IOSTAT=io_stat) 
    else
      ! do nothing
    end if
  else
    if(root_processor) write(lun,NML=gyroaverage)

    call write_run_parameter('gyroaverage', 'n_points_ga', n_points_ga)
    call write_run_parameter('gyroaverage', 'gyro_average_ions', gyro_average_ions)
    call write_run_parameter('gyroaverage', 'gyro_average_electrons', gyro_average_electrons)
    call write_run_parameter('gyroaverage', 'orb_polarize', orb_polarize)
    call write_run_parameter('gyroaverage', 'mk_gyro_av_hermitian', mk_gyro_av_hermitian)
    call write_run_parameter('gyroaverage', 'use_conj', use_conj)
    call write_run_parameter('gyroaverage', 'consistent_long_wave', consistent_long_wave)
    call write_run_parameter('gyroaverage', 'parallel_mod', parallel_mod)
    call write_run_parameter('gyroaverage', 'blending_order', blending_order)
    call write_run_parameter('gyroaverage', 'n_gav_bound_ex', n_gav_bound_ex)

  end if 

end subroutine gyro_average_read_nml

!------------------------------------------------------------------------------
!> bcast the linear terms namelist params, and do checks
!------------------------------------------------------------------------------

subroutine gyro_average_bcast_nml

  use mpiinterface, only : mpibcast 
  use dist,         only : n_gav_bound_ex
  
  call mpibcast(n_points_ga,            1)
  call mpibcast(gyro_average_ions,      1)
  call mpibcast(gyro_average_electrons, 1)
  call mpibcast(orb_polarize,           1)
  call mpibcast(mk_gyro_av_hermitian,   1) 
  call mpibcast(use_conj,               1) 
  call mpibcast(consistent_long_wave,   1)
  call mpibcast(blending_order,         1)  
  call mpibcast(parallel_mod,           1)
  call mpibcast(n_gav_bound_ex,         1)
   
end subroutine gyro_average_bcast_nml


!------------------------------------------------------------------------------
!> This routine allocates the arrays needed for the gyro-average
!------------------------------------------------------------------------------

subroutine gyro_average_allocate  

  use control,        only : spectral_radius, nlapar
  use grid,           only : nmod, nx, ns, nmu, nvpar, nsp, n_x_grid
  use grid,           only : mpisetup_mod 
  use dist,           only : n_gav_bound
  use components,     only : energetic_particles
  use mpidatatypes,   only : type_mod_send, type_mod_recv 
  use rho_par_switch, only : s_average 
  use general,        only : gkw_abort
  use matrix_format, only : create_matrix, matrix_format_gkwcrs
  
  integer :: nelem, ierr
  !> maximal number of elements
  integer :: max_nmatred, max_nmatred_s, max_nmatav

  ! get first, last and number of toroidal modes to do on this processor
  call mpisetup_mod(first_imod,last_imod,nmod_l,TYPE_MOD_SEND,TYPE_MOD_RECV,parallel_mod)

  allocate(nind_gyro_av(nmod,nx,ns,nmu,nsp),stat=ierr)
  if (ierr /= 0) call gkw_abort('can not allocate nind_gyro_av')

  allocate(index_gyro_av(nmod,nx,ns,nmu,nsp,n_points_ga+n_points_ga_additional),stat=ierr)
  if (ierr /= 0) call gkw_abort('can not allocate index_gyro_av')

  allocate(coeff_gyro_av(nmod,nx,ns,nmu,nsp,n_points_ga+n_points_ga_additional),stat=ierr)
  if (ierr /= 0) call gkw_abort('can not allocate coeff_gyro_av')
  
  max_nmatav = nmod*nx*ns*nmu*nsp*n_points_ga 
  if (s_average)  max_nmatav = n_s_ga * max_nmatav
  if (nlapar) max_nmatav = 2*max_nmatav
  
  matav = create_matrix("matrix for gyroavg of the field", max_nmatav)

  nelem = 6*nmod*nx*ns*nmu*nvpar*nsp
  if (nlapar) nelem = 2*nelem

  matint = create_matrix("matrix for integral of the distribution", nelem, &
     & matrix_format_gkwcrs)

  max_nmatred = nmod*nx*ns*nmu*nsp*n_points_ga
  if (nlapar) max_nmatred = 2*max_nmatred

  matred = create_matrix('matrix for the gyroaveraged charge density calculation',&
     & max_nmatred, matrix_format_gkwcrs)
  
  if (s_average) max_nmatred_s=max_nmatred*n_s_ga

  ! arrays for rearranging toroidal modes after fields calc
  nelem = nmod_l*nx*ns
  if (nlapar) nelem = 2*nelem
  
  allocate(iiscatz(nelem),stat=ierr) 
  if (ierr /= 0) call gkw_abort('can not allocate iiiscatz')
 
  allocate(iilocz(nelem),stat=ierr) 
  if (ierr /= 0) call gkw_abort('can not allocate iiilocz') 

  allocate(bn_Gx(1:n_x_grid,ns), stat = ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate bn_Gx in gyro-average')

  allocate(fmaxwl_Gx(1:n_x_grid,ns,nmu,nvpar,nsp), stat = ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate fmaxwl_Gx in gyro-average')

  if (energetic_particles) then
    allocate(fEP_Gx(1:n_x_grid,ns,nmu,nvpar,nsp), stat = ierr)
    if (ierr /= 0) call gkw_abort('could not allocate fEP_Gx in gyro-average')
  end if

  allocate(coeff_gyro_av_Gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,ns,nmu, &
          &                 nsp,n_points_ga+n_points_ga_additional), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate coeff_gyro_av_Gx in gyro-average')

  allocate(index_gyro_av_Gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,ns,nmu, &
          &                 nsp,n_points_ga+n_points_ga_additional), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate index_gyro_av_Gx in gyro-average')

  allocate(nind_gyro_av_Gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,ns,nmu, &
          &                nsp), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate nind_gyro_av_Gx in gyro-average')

  allocate(coeff_gyro_cj_Gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,ns,nmu, &
          &                 nsp,n_points_ga+n_points_ga_additional), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate coeff_gyro_cj_Gx in gyro-average')

  allocate(index_gyro_cj_Gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,ns,nmu, &
          &                 nsp,n_points_ga+n_points_ga_additional), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate index_gyro_cj_Gx in gyro-average')

  allocate(nind_gyro_cj_Gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,ns,nmu, &
          &                nsp), stat = ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate nind_gyro_cj_Gx in gyro-average')  

  if(s_average) then
    nelem = nmod*n_x_grid*ns*n_x_grid*(n_s_ga-1)
    if (nlapar) nelem = nelem + nmod_l*n_x_grid*ns*n_x_grid

    matpol = create_matrix("matrix for s component of polarisation" , nelem,&
       & matrix_format_gkwcrs)

    allocate(nind_gyro_av_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate nind_gyro_av_s_gx')

    allocate(index_ix_gyro_av_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp,n_points_ga*n_s_ga),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate index_ix_gyro_av_s_gx')

    allocate(index_i_gyro_av_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp,n_points_ga*n_s_ga),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate index_i_gyro_av_s_gx')

    allocate(coeff_gyro_av_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp,n_points_ga*n_s_ga),stat=ierr) 
    if (ierr /= 0) call gkw_abort('can not allocate coeff_gyro_av_s_gx')

    allocate(nind_gyro_cj_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate nind_gyro_av_s_gx')

    allocate(index_ix_gyro_cj_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp,n_points_ga*n_s_ga),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate index_ix_gyro_av_s_gx')

    allocate(index_i_gyro_cj_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp,n_points_ga*n_s_ga),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate index_i_gyro_av_s_gx')

    allocate(coeff_gyro_cj_s_gx(nmod,-n_gav_bound:n_x_grid+n_gav_bound+1,0:ns+1,nmu,nsp,n_points_ga*n_s_ga),stat=ierr) 
    if (ierr /= 0) call gkw_abort('can not allocate coeff_gyro_av_s_gx')


    matred_s = create_matrix('coefficient for the gyro average of the &
       & distribution in the Poisson equation', max_nmatred_s, &
       & matrix_format_gkwcrs)

    allocate(nind_gyro_av_s(nmod,nx,ns,nmu,nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate nind_gyro_av_s')

    allocate(index_ix_gyro_av_s(nmod,nx,ns,nmu,nsp,n_points_ga*n_s_ga),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate index_ix_gyro_av_s')

    allocate(index_i_gyro_av_s(nmod,nx,ns,nmu,nsp,n_points_ga*n_s_ga),stat=ierr)
    if (ierr /= 0) call gkw_abort('can not allocate index_i_gyro_av_s')

    allocate(coeff_gyro_av_s(nmod,nx,ns,nmu,nsp,n_points_ga*n_s_ga),stat=ierr) 
    if (ierr /= 0) call gkw_abort('can not allocate coeff_gyro_av_s')

  endif
end subroutine gyro_average_allocate


!------------------------------------------------------------------------------
!> The initialization routine of the gyro-average 
!------------------------------------------------------------------------------

subroutine gyro_average_init 
  use control,        only : method, nlapar
  use grid,           only : lsendrecv_x, n_procs_s
  use rho_par_switch, only : s_average
  use general,        only : gkw_abort
 
  if (method == 'IMP') then
    call gkw_abort('Need method=EXP for non-spectral')
  end if
  if ( s_average.and. (n_procs_s .ne. 1) ) then
    call gkw_abort ('Parallelisation over s-direction is not allowed for &
       &s_average=T')
  endif
  if ( s_average.and.  lsendrecv_x) then 
    !problem with boundary condition and matpol%ii
    call gkw_abort ('Parrallelisation over x-direction is not allowed for &
       &s_average=T')
  endif
  if ( nlapar .and.s_average) then
    call gkw_abort('apar for s_average is not implemented')
  end if

  ! without s-direction the gyro average is only in ix- and zeta.
  if (.not.s_average) then
    n_s_ga = 1
  else
    n_s_ga = 3
  end if
     
  ! calculate the coefficients 
  call gyro_average_indx 

  call index_diagnostics()
  
  ! setup the arays for rearranging the toroidal modes after the fields solve
  call setup_indx_nmod

  ! set up the arrays for the gyro-average of the field
  ! LOCAL IN NMOD, LOCAL IN X
  call gyro_average_fields_init 

  ! set up the arrays for the integral of the poisson and ampere' equation 
  !GLOBAL in NMOD
  call field_integral_eqs 

  ! set up the polarization
  ! LOCAL IN NMOD, GLOBAL IN X
  call polarization_init
  
end subroutine gyro_average_init

!------------------------------------------------------------------------------
!> Some diagnostics (if desired). This diagnostic is useful
!> when parallelizing over the ix direction, but the number of
!> ghost points has to be predetermined in dist.F90
!------------------------------------------------------------------------------
subroutine index_diagnostics()
  use grid, only : nmod, nx, ns, nmu, nsp, lsendrecv_x
  use mpiinterface, only : mpiallreduce_max, root_processor
  use global, only : root_and_verbose
  use rho_par_switch, only : s_average
  use dist, only : xgp => ghost_points_x
  integer :: imod, ix, i, j, is, m, ind_max_dif, im, ip, iout

  ind_max_dif = 0
  do imod=1,nmod; do ix=1,nx; do i=1,ns; do j=1,nmu; do is =1,nsp
    do m = 1, nind_gyro_av(imod,ix,i,j,is)
      im = abs(index_gyro_av(imod,ix,i,j,is,m)-ix)
      if (.not. lsendrecv_x) then
        ip = abs(im-nx)
        im = min(im,ip)
      end if
      ind_max_dif = max(ind_max_dif,im)
    end do
    if (s_average) then
      do m = 1, nind_gyro_av_s(imod,ix,i,j,is)
        im = abs(index_ix_gyro_av_s(imod,ix,i,j,is,m)-ix)
        if (.not. lsendrecv_x) then
          ip = abs(im-nx)
          im = min(im,ip)
        end if
        ind_max_dif = max(ind_max_dif,im)
      end do
    endif
  end do; end do; end do; end do; end do;
  
  call mpiallreduce_max(ind_max_dif,iout,1)
  ! In some cases iout can exceed n_gav_bound_measure by one - why ?
  
  if (root_and_verbose .or. (lsendrecv_x .and. root_processor)) then
    write(*,*)
    write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*)'Calculated the coefficients for the gyro average'
    write(*,*)'Maximum distance in the x-grid:  ', iout
    write(*,*)'Ghost points in the x-grid:      ', xgp
    write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*)
  endif
  
  ! check we have enough ghost cells for the gyro-average
  ! (would not need repeating if n_gav_bound_measure was correct)
  if (lsendrecv_x .and. iout /= xgp) call xgrid_check(iout,xgp,'x ghost cells')

end subroutine index_diagnostics

!------------------------------------------------------------------------------
!> This routine builds the arrays for the gyro average. Here, simply the 
!> integral over the circular orbit is implemented. 
!------------------------------------------------------------------------------
subroutine gyro_average_indx 
  use global,        only : int2char
  use grid,           only : nmod, nx, ns, nmu, nsp, n_x_grid, gx, gs
  use velocitygrid,   only : mugr 
  use components,     only : mas, tgrid, signz, rhostar
  use geom,           only : bn, metric, efun, dxgr, jacobian_G, metric_G
  use constants,      only : pi 
  use dist,           only : n_gav_bound 
  use mpiinterface,   only : mpiallreduce_max, root_processor
  use global,         only : r_tiny 
  use rho_par_switch, only : s_average 
  use general,        only : gkw_abort

  !> help array,
  !> collects gyro-av element sum for each jj ix point 
  !> for ii point currently being calculated  
  complex, allocatable :: gyro_interp_coeffs(:)
  complex, allocatable :: gyro_interp_coeffs_s(:,:)
  real    :: gyro_angle, rho, delta_zeta, delta_psi, max_val, deviation, dev
  real    :: delta_s, jacix, bnix, jacxp, bnxp 
  real    :: xderiv, yderiv
  logical :: found 
  integer :: imod, ix, i, j, is, m, ierr, nind, nind_s, l, ig
  integer :: idum, ixp, mix, mi
  complex :: cdum

  ! initialize the index and coefficient array 
  nind_gyro_av  = 0
  index_gyro_av = 0
  coeff_gyro_av = 0.E0
  if(s_average) then
    nind_gyro_av_s  = 0 
    index_i_gyro_av_s = 0 
    index_ix_gyro_av_s = 0 
    coeff_gyro_av_s = 0.E0  
  endif

  allocate(gyro_interp_coeffs(-n_gav_bound:nx+n_gav_bound+1),stat=ierr) 
  if(s_average) allocate(gyro_interp_coeffs_s(-n_gav_bound:nx+n_gav_bound+1,0:ns+1),stat=ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate gyro_interp_coeffs') 

  ! the gyroaverage is a function of position space coordinates, species and mu
  do imod = 1,nmod; do ix = 1,nx; do i = 1,ns 

    ! global s point 
    ig = gs(i) 

    ! The second order correction to the centre of the gyro-ring in 
    ! the field aligned hamada coordinates. Using this correction the 
    ! proper long wave length limit of the polarization is obtained. 
    if (consistent_long_wave) then 
      if (gx(ix) == 1) then 
        xderiv = (jacobian_G(gx(ix+1))*metric_G(gx(ix+1),ig,1,1)-   &
               & jacobian_G(gx(ix))*metric_G(gx(ix),ig,1,1))        &
               & / (1.E0*dxgr*jacobian_G(gx(ix))) 
        yderiv = (jacobian_G(gx(ix+1))*metric_G(gx(ix+1),ig,1,2) -  &
               & jacobian_G(gx(ix))*metric_G(gx(ix),ig,1,2))        &
               & / (1.E0*dxgr*jacobian_G(gx(ix))) 
      else 
        if (gx(ix) == n_x_grid) then 
          xderiv = (jacobian_G(gx(ix))*metric_G(gx(ix),ig,1,1) -    &
                 & jacobian_G(gx(ix-1))*metric_G(gx(ix-1),ig,1,1))  &
                 & / (1.E0*dxgr*jacobian_G(gx(ix))) 
          yderiv = (jacobian_G(gx(ix))*metric_G(gx(ix),ig,1,2) -    &
                 & jacobian_G(gx(ix-1))*metric_G(gx(ix-1),ig,1,2))  &
                 & / (1.E0*dxgr*jacobian_G(gx(ix))) 
        else 
          xderiv = (jacobian_G(gx(ix+1))*metric_G(gx(ix+1),ig,1,1)- &
                 & jacobian_G(gx(ix-1))*metric_G(gx(ix-1),ig,1,1))  &
                 & / (2.E0*dxgr*jacobian_G(gx(ix))) 
          yderiv = (jacobian_G(gx(ix+1))*metric_G(gx(ix+1),ig,1,2)- &
                 & jacobian_G(gx(ix-1))*metric_G(gx(ix-1),ig,1,2))  &
                 & / (2.E0*dxgr*jacobian_G(gx(ix))) 
        endif
      endif 
    else 
      xderiv = 0. 
      yderiv = 0. 
    endif 

    do j = 1,nmu; do is = 1,nsp 

      if ((signz(is)<0 .and. gyro_average_electrons) .or. &
          & (signz(is)>0 .and. gyro_average_ions)) then 

        ! to compute the gyroaverage at a given position ix, imod, i, j, is,
        ! an average along a circular path is calculated.
        
        ! the Larmor radius is obtained from background quantities
        rho = sqrt(2.E0*mas(is)*tgrid(is)*mugr(j)/bn(ix,i)) / abs(signz(is))

        ! initialize array for accumulating the coefficients of one
        ! sampling point
        gyro_interp_coeffs = 0.E0
        if(s_average) gyro_interp_coeffs_s = 0.E0 

        ! loop over all the sampling points of the orbit
        do m = 1, n_points_ga 
 
          ! the gyro angle corresponding to the sampling point is 
          gyro_angle = 2.E0*pi*real(m-0.5E0)/real(n_points_ga) 

          ! calculate the distance of the sampling point to the
          ! gyrocenter in both directions
          delta_zeta = cos(gyro_angle) * metric(ix,i,1,2) * rho / sqrt(metric(ix,i,1,1)) + &
           & sin(gyro_angle) * 2.E0 * bn(ix,i) * efun(ix,i,1,2) * rho / sqrt(metric(ix,i,1,1)) &
           & + 0.25 * rho**2 * yderiv 
          delta_psi = cos(gyro_angle) * rho * sqrt(metric(ix,i,1,1)) + 0.25*rho**2*xderiv
          ! if desired calculated the distance in s direction
          if ( s_average ) then
            delta_s  = rhostar* rho/sqrt(metric(ix,i,1,1)) * ( cos(gyro_angle) * metric(ix,i,1,3) + &
                   &         2.E0 * sin(gyro_angle) * bn(ix,i) * efun(ix,i,1,3) )
          else
            delta_s = 0
          endif 

          ! calculate the interpolation coefficients, which are then
          ! used to obtain the field value at the sample point.
          if (s_average) then
            call calculate_gyro_coef(imod,ix,i,delta_zeta,delta_psi,delta_s, &
               & gyro_interp_coeffs,gyro_interp_coeffs_s)
          else
            call calculate_gyro_coef(imod,ix,i,delta_zeta,delta_psi,delta_s, &
               & gyro_interp_coeffs) 
          end if

        end do
        ! the array gyro_interp_coeffs now contains the accumulated
        ! interpolation coefficients for each local radial gridpoint,
        ! needed to compute the average along the orbit around the
        ! current gyrocenter position imod, ix, i, is, j.

!       cdum = 0.
!       do m = - n_gav_bound, nx + n_gav_bound + 1 
!         cdum = cdum + gyro_interp_coeffs(m) 
!       end do 
!       cdum = cdum / n_points_ga 
!       if (abs(cdum) > 1.E0) then 
!         write(*,*)ix,i,j,is,cdum 
!       endif 

        nind = 0
        ! look at all radial gridpoints, incl. necessary ghost cells
        do ixp = -n_gav_bound, nx+n_gav_bound+1
          if (abs(gyro_interp_coeffs(ixp)) > r_tiny) then
            ! if this grid point contributes, then increment the counter...
            nind = nind + 1
            if (nind > n_points_ga+n_points_ga_additional) then
              call gkw_abort('Too few indices for gyro-av: Workaround: &
                 &Increase n_points_gav') 
            end if
            ! ... and store this radial index together with its coefficient
            index_gyro_av(imod,ix,i,j,is,nind) = ixp
            ! this division by number of sampling points is
            ! where the average over all sampling values is done.
            coeff_gyro_av(imod,ix,i,j,is,nind) = gyro_interp_coeffs(ixp) &
               & / n_points_ga
          endif 
        end do
        ! the number of grid points which contribute to the
        ! gyroaverage at that position
        nind_gyro_av(imod,ix,i,j,is) = nind

        ! the coeff for the average in the s direction
        if (s_average) then
          nind_s = 0
          do mi = 0, ns+1
            do mix = -n_gav_bound, nx+n_gav_bound+1
              if (abs(gyro_interp_coeffs_s(mix,mi)) > r_tiny)then
                ! save the coeff in the lowest order first
                nind_s = nind_s + 1
                if (nind_s > n_points_ga*n_s_ga) then
                  call gkw_abort('Too few indices for gyro-av: Workaround: &
                     & Increase n_points_gav')
                end if
                index_ix_gyro_av_s(imod,ix,i,j,is,nind_s) = mix
                index_i_gyro_av_s(imod,ix,i,j,is,nind_s) = mi
                nind_gyro_av_s(imod,ix,i,j,is) = nind_s
                coeff_gyro_av_s(imod,ix,i,j,is,nind_s) = gyro_interp_coeffs_s(mix,mi) / n_points_ga
!              write(1,*) imod,ix,i,j,is, coeff_gyro_av_s(imod,ix,i,j,is,nind_s)
              endif 
            end do 
          end do 
        end if

      else
        ! this species shall not be gyroaveraged.
        ! then the gyroaveraging will just be an identity matrix.
        nind_gyro_av(imod,ix,i,j,is) = 1
        index_gyro_av(imod,ix,i,j,is,1) = ix
        coeff_gyro_av(imod,ix,i,j,is,1) = 1.0
        if(s_average) nind_gyro_av_s(imod,ix,i,j,is) = 0

      endif 

    end do; end do 

  end do; end do; end do 

  ! Find the global maximum n_gav_bound_measure
  call mpiallreduce_max(n_gav_bound_measure,idum,1)
  n_gav_bound_measure = idum

  ! report and check.
  ! WARNING n_gav_bound_measure sometimes underestimates by one
  ! compared to ind_max_dif
  if (root_processor) then
    write(*,*)
    write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*)'Calculated the coefficients for the gyro average'
    write(*,*)'Maximum points over the x-boundary: ', n_gav_bound_measure
    write(*,*)'No. of boundary points allocated:   ', n_gav_bound
    write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) 
  end if
  if (n_points_ga < n_gav_bound_measure) then
    call gkw_abort('n_points_ga must be equal to or greater than the maximum number of&
                  & points over the x-boundary, here '//trim(int2char(n_gav_bound_measure)))
  end if

  if (n_gav_bound_measure > n_gav_bound) then 
    call xgrid_check(n_gav_bound_measure,n_gav_bound,'radial boundary points')      
  end if

  ! deallocate the help array 
  deallocate(gyro_interp_coeffs)
  if (s_average) deallocate(gyro_interp_coeffs_s)

  ! do the global gathering (can be done after the gyro average matrix is calculated
  call gather_global_x_arrays

  ! now manipulate the global gyro-average matrix. First make sure that if points 
  ! across the boundary exist there is an entry that refers back to this point.
  ! This is important only for non-periodic boundary conditions, and allows 
  ! one to do the double gyro-average. 

  ! loop over local position space, mu and species points
  do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do is = 1, nsp
    ! but global x space.
    do ix = 1, n_x_grid
      ! loop over all gridpoints contributing to the gyroavg for that point.
      do m = 1, nind_gyro_av_Gx(imod,ix,i,j,is)
        ! is the point across the lower boundary of the global grid?
        if (index_gyro_av_Gx(imod,ix,i,j,is,m) < 1) then
          ixp = index_gyro_av_Gx(imod,ix,i,j,is,m)
          ! so lets add one more contributing grid point for
          ! that point on the other side...
          nind = nind_gyro_av_Gx(imod,ixp,i,j,is) + 1
          ! ... and that makes nind in total (up to now)
          nind_gyro_av_Gx(imod,ixp,i,j,is) = nind
          ! ... and it is the point we are looking at right now which
          ! contributes.
          index_gyro_av_Gx(imod,ixp,i,j,is,nind) = ix 
          coeff_gyro_av_Gx(imod,ixp,i,j,is,nind) = conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m))
        else if (index_gyro_av_Gx(imod,ix,i,j,is,m) > n_x_grid) then  
          ixp = index_gyro_av_Gx(imod,ix,i,j,is,m) 
          nind = nind_gyro_av_Gx(imod,ixp,i,j,is) + 1 
          nind_gyro_av_Gx(imod,ixp,i,j,is) = nind 
          index_gyro_av_Gx(imod,ixp,i,j,is,nind) = ix 
          coeff_gyro_av_Gx(imod,ixp,i,j,is,nind) = conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m))        
        endif 
      end do 
    end do 
  end do; end do; end do; end do
  ! note that the rest of the orbit sampling points for the ixp added
  ! above give zero contributions, which reflects the Dirichlet
  ! boundary condition.

  if (s_average) then
    do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do is = 1, nsp 
      do ix = 1, n_x_grid 
        do m = 1, nind_gyro_av_s_Gx(imod,ix,i,j,is) 
          ! is the point across the lower boundary 
          if (index_ix_gyro_av_s_Gx(imod,ix,i,j,is,m) < 1) then  
            ixp = index_ix_gyro_av_s_Gx(imod,ix,i,j,is,m) 
            nind_s = nind_gyro_av_s_gx(imod,ixp,i,j,is) + 2
            nind_gyro_av_s_Gx(imod,ixp,i,j,is) = nind_s 
            index_ix_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s) = ix 
            index_ix_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s-1) = ix 
            index_i_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s) = index_i_gyro_av_s_Gx(imod,ixp,i,j,is,m) 
            index_i_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s-1) = index_i_gyro_av_s_Gx(imod,ixp,i,j,is,m) 
            coeff_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s) = conjg(coeff_gyro_av_s_Gx(imod,ix,i,j,is,m))
            coeff_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s-1) = conjg(coeff_gyro_av_s_Gx(imod,ix,i,j,is,m))
          endif 
          if (index_ix_gyro_av_s_Gx(imod,ix,i,j,is,m) > n_x_grid) then  
            ixp = index_ix_gyro_av_s_Gx(imod,ix,i,j,is,m) 
            nind_s = nind_gyro_av_s_Gx(imod,ixp,i,j,is) + 2 
            nind_gyro_av_s_Gx(imod,ixp,i,j,is) = nind_s 
            index_ix_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s) = ix 
            index_ix_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s-1) = ix 
            index_i_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s) = index_i_gyro_av_s_Gx(imod,ixp,i,j,is,m) 
            index_i_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s-1) = index_i_gyro_av_s_Gx(imod,ixp,i,j,is,m) 
            coeff_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s) = conjg(coeff_gyro_av_s_Gx(imod,ix,i,j,is,m))        
            coeff_gyro_av_s_Gx(imod,ixp,i,j,is,nind_s-1) = conjg(coeff_gyro_av_s_Gx(imod,ix,i,j,is,m))        
          endif 
        end do 
      end do 
    end do; end do; end do; end do 
  end if

  ! If requested make the gyro-average operator Hermitian.  That is,
  ! the matrix is equal its transposed&complexconjugated matrix.
  if (mk_gyro_av_hermitian) then 
    if (s_average) call gkw_abort('hermitian gyro average operator and s_average is not supported jet')
    ! Then make the gyro_average matrix Hermitian 
    max_val = 0.
    deviation = 0.
    ! loop over all points of the gyroaveraged field
    do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do is = 1, nsp
      ! including the points across the boundary
      do ix = -n_gav_bound, n_x_grid+n_gav_bound+1
        ! for each point, loop over all the contributing sample points
        do m = 1, nind_gyro_av_Gx(imod,ix,i,j,is)
          ! get the x index of the sampling point
          ixp = index_gyro_av_Gx(imod,ix,i,j,is,m)
          if (ixp == ix) then
            ! just collect some statistics about our matrix
            max_val = max(max_val,abs(coeff_gyro_av_Gx(imod,ix,i,j,is,m)))
            deviation = max(deviation,abs(aimag(coeff_gyro_av_Gx(imod,ix,i,j,is,m))))
            ! if the sampling point is the point itself then further
            ! below this will become a matrix element on the main diagonal.
            ! For a Hermitian matrix this must be a purely real number.
            coeff_gyro_av_Gx(imod,ix,i,j,is,m) = real(coeff_gyro_av_Gx(imod,ix,i,j,is,m))
          endif
          if (ixp > ix) then
            ! if the sampling point is in the upper triangular part then
            found = .false.
            ! now in turn search all contributing points of this sample point
            do l = 1, nind_gyro_av_Gx(imod,ixp,i,j,is)
              if (index_gyro_av_Gx(imod,ixp,i,j,is,l) == ix) then 
                ! and when the sample point is met again then we now
                ! know have two matrix elements which should be the
                ! complex conjugates of each other for a Hermitian
                ! matrix.
                found = .true.
                ! keep collecting some statistics about our matrix
                max_val = max(max_val,abs(coeff_gyro_av_Gx(imod,ixp,i,j,is,l))) 
                dev =   abs(coeff_gyro_av_Gx(imod,ixp,i,j,is,l) - &
                    & conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m)))
                if (dev > deviation) deviation = dev
                ! Now, the modified value has the averaged real and imag part.
                cdum  = 0.5*coeff_gyro_av_Gx(imod,ix,i,j,is,m) + &
                      & 0.5*conjg(coeff_gyro_av_Gx(imod,ixp,i,j,is,l)) 
                coeff_gyro_av_Gx(imod,ixp,i,j,is,l) = 0.5*coeff_gyro_av_Gx(imod,ixp,i,j,is,l)+&
                      & 0.5*conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m))
                coeff_gyro_av_Gx(imod,ix,i,j,is,m) = cdum 
              endif
            end do 
            if (.not.found) then
              ! add a completely new element, whose value is the
              ! complex conjugate.
              nind = nind_gyro_av_Gx(imod,ixp,i,j,is) + 1 
              nind_gyro_av_Gx(imod,ixp,i,j,is) = nind 
              index_gyro_av_Gx(imod,ixp,i,j,is,nind) = ix 
              coeff_gyro_av_Gx(imod,ixp,i,j,is,nind) = conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m))
              deviation = max(deviation,abs(coeff_gyro_av_Gx(imod,ix,i,j,is,m)))
            endif
          else 
            found = .false.
            ! in turn loop all elements contributing to this sample point
            do l = 1, nind_gyro_av_Gx(imod,ixp,i,j,is)
              ! and when the sample point is met again
              if (index_gyro_av_Gx(imod,ixp,i,j,is,l) == ix) then 
                found = .true.
                dev =   abs(coeff_gyro_av_Gx(imod,ixp,i,j,is,l) - &
                    & conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m)))
                if (dev > r_tiny) call gkw_abort('Error in making gyro-average Hermitian') 
              endif 
            end do 
            if (.not.found) then
              ! if it is never met, i.e. this is a matrix element in
              ! the lower triangular part which does not have a
              ! corresponding element in the upper triangle; then add
              ! a new complex conjugated element there
              nind = nind_gyro_av_Gx(imod,ixp,i,j,is) + 1 
              nind_gyro_av_Gx(imod,ixp,i,j,is) = nind 
              index_gyro_av_Gx(imod,ixp,i,j,is,nind) = ix 
              coeff_gyro_av_Gx(imod,ixp,i,j,is,nind) = conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m))
              deviation = max(deviation,abs(coeff_gyro_av_Gx(imod,ix,i,j,is,m)))
            endif           
          endif
        end do 
      end do 
    end do; end do; end do; end do; 

    if (root_processor) then    
      write(*,*) 'Made the Gyro-average matrix Hermitian: correction of order ',& 
               & deviation / max_val 
    endif 

    ! store the local gyro_average array 
    do imod = 1, nmod; do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do is = 1, nsp 
      ixp = gx(ix) 
      nind_gyro_av(imod,ix,i,j,is) = nind_gyro_av_Gx(imod,ixp,i,j,is) 
      do m = 1, nind_gyro_av(imod,ix,i,j,is)
        index_gyro_av(imod,ix,i,j,is,m) = index_gyro_av_Gx(imod,ixp,i,j,is,m) 
        coeff_gyro_av(imod,ix,i,j,is,m) = coeff_gyro_av_Gx(imod,ixp,i,j,is,m) 
      end do 
    end do; end do; end do; end do; end do 
  endif 

  call check_count_index_gyro_av_gx

  ! build the conjugated operator. 
  if (use_conj) then 
    if(s_average) call gkw_abort ('use_conj .and. s_average  is not done' ) 

    ! initialize 
    nind_gyro_cj_Gx  = 0
    index_gyro_cj_Gx = 0
    coeff_gyro_cj_Gx = 0 

    do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do is = 1, nsp 
      do ix = -n_gav_bound, n_x_grid+n_gav_bound+1
        ! Determine jacobian and magnetic field 
        if (ix < 1) then 
          jacix = jacobian_G(1)
          bnix  = bn_Gx(1,i) 
        else if (ix > n_x_grid) then 
          jacix = jacobian_G(n_x_grid) 
          bnix  = bn_Gx(n_x_grid,i) 
        else 
          jacix = jacobian_G(ix) 
          bnix  = bn_Gx(ix,i)
        endif 
        do m = 1, nind_gyro_av_Gx(imod,ix,i,j,is) 
          ixp = index_gyro_av_Gx(imod,ix,i,j,is,m)
          if (ixp < 1) then  
            jacxp = jacobian_G(1)
            bnxp  = bn_Gx(1,i)
          else if (ixp > n_x_grid) then 
            jacxp = jacobian_G(n_x_grid) 
            bnxp  = bn_Gx(n_x_grid,i)
          else 
            jacxp = jacobian_G(ixp) 
            bnxp  = bn_Gx(ixp,i)
          endif
          nind = nind_gyro_cj_Gx(imod,ixp,i,j,is)+1 
          nind_gyro_cj_Gx(imod,ixp,i,j,is) = nind 
          index_gyro_cj_Gx(imod,ixp,i,j,is,nind) = ix 
          coeff_gyro_cj_Gx(imod,ixp,i,j,is,nind) = jacix*bnix/(jacxp*bnxp)* &
               &                    conjg(coeff_gyro_av_Gx(imod,ix,i,j,is,m))
        end do 
      end do 
    end do; end do; end do; end do 

  else

    ! the conjugated operator is identical to the gyro-average operator

    nind_gyro_cj_Gx  = 0
    index_gyro_cj_Gx = 0
    coeff_gyro_cj_Gx = 0 

    do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do is = 1, nsp 
      do ix = -n_gav_bound, n_x_grid+n_gav_bound+1
        nind_gyro_cj_Gx(imod,ix,i,j,is) = nind_gyro_av_Gx(imod,ix,i,j,is) 
        do m = 1, nind_gyro_av_Gx(imod,ix,i,j,is) 
          index_gyro_cj_Gx(imod,ix,i,j,is,m) = index_gyro_av_Gx(imod,ix,i,j,is,m) 
          coeff_gyro_cj_Gx(imod,ix,i,j,is,m) = coeff_gyro_av_Gx(imod,ix,i,j,is,m)
        end do 
      end do 
    end do; end do; end do; end do 

    if (s_average) then

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do is = 1, nsp 
        do ix = -n_gav_bound, n_x_grid+n_gav_bound+1
          nind_gyro_cj_s_Gx(imod,ix,i,j,is) = nind_gyro_av_s_Gx(imod,ix,i,j,is) 
          do m = 1, nind_gyro_av_s_Gx(imod,ix,i,j,is) 
            index_ix_gyro_cj_s_Gx(imod,ix,i,j,is,m) = index_ix_gyro_av_s_Gx(imod,ix,i,j,is,m) 
            index_i_gyro_cj_s_Gx(imod,ix,i,j,is,m) = index_i_gyro_av_s_Gx(imod,ix,i,j,is,m) 
            coeff_gyro_cj_s_Gx(imod,ix,i,j,is,m) = coeff_gyro_av_s_Gx(imod,ix,i,j,is,m)
          end do 
        end do 
      end do; end do; end do; end do  
    endif
  endif 

end subroutine gyro_average_indx 

!------------------------------------------------------------------------------
!> For use_conj=.true., count how often a certain x-index appears in
!> index_gyro_av_gx. If this exceeds n_points_ga + n_points_ga_additional
!> then a segfault would occour, thus abort.
!------------------------------------------------------------------------------
subroutine check_count_index_gyro_av_gx
  use dist, only : n_gav_bound
  use general, only : gkw_abort
  use global, only : int2char
  use grid, only : nmod, ns, nmu, nsp, n_x_grid, gx, gs
  use mpicomms, only : COMM_X_EQ
  use mpiinterface, only : mpiallreduce_max

  integer, allocatable :: counter(:), maxcounter(:)
  integer :: ierr
  integer :: imod, i, j, is, ix, gyroangle

  if (use_conj) then
    allocate(counter(-n_gav_bound:n_x_grid+n_gav_bound+1), stat = ierr)
    if (ierr /= 0) call gkw_abort('could not allocate counter in check_count_index_gyro_av_gx')
    allocate(maxcounter(-n_gav_bound:n_x_grid+n_gav_bound+1), stat = ierr)
    if (ierr /= 0) call gkw_abort('could not allocate maxcounter in check_count_index_gyro_av_gx')

    do imod = 1,nmod
      do i = 1,ns
        do j = 1,nmu
          do is = 1,nsp
            counter = 0
            maxcounter = 0
            do ix = -n_gav_bound, n_x_grid+n_gav_bound+1
              do gyroangle = 1, nind_gyro_av_gx(imod,ix,i,j,is)
                counter(index_gyro_av_gx(imod,ix,i,j,is,gyroangle)) = &
                     & counter(index_gyro_av_gx(imod,ix,i,j,is,gyroangle)) + 1
              end do
            end do

            ! This may seem wrong, as with paralellization values for different
            ! s/mu/sp might be compared, but as only the global maximum is of
            ! importance, it should be working correct.
            call mpiallreduce_max(counter, maxcounter, n_x_grid+2*(n_gav_bound+1), COMM_X_EQ)

            if (maxval(maxcounter) > (n_points_ga + n_points_ga_additional)) then
              call gkw_abort('Number of references ('//int2char(maxval(maxcounter)) &
                      & //') for x-location '&
                      & //int2char(maxloc(maxcounter,1)-n_gav_bound-1)&
                      & //' in index_gyro_av_gx exceeds the array size of &
                      &index_gyro_cj_gx ('&
                      & //int2char(n_points_ga + n_points_ga_additional) &
                      & //'). This would lead to a segfault.')
            end if
          end do
        end do
      end do
    end do

    if (allocated(counter)) deallocate(counter)
    if (allocated(maxcounter)) deallocate(maxcounter)
  end if
end subroutine check_count_index_gyro_av_gx

!------------------------------------------------------------------------------
!> Routine that calculates the contribution to the gyro-averaged field
!> at a single sampling point (imod, ix, i) using blending functions of order 
!> 2 (linear interpolation), 3 and 4. The finite rho* terms is always 
!> calculated using linear interpolation. The routine repeatedly adds 
!> to the output array gyro_interp_coeffs
!------------------------------------------------------------------------------
subroutine calculate_gyro_coef(imod,ix,i,delta_zeta,delta_psi,delta_s,gyro_interp_coeffs,  & 
                    &     gyro_interp_coeffs_s)   !Don't pass an unallocated variable

  use grid,           only : nx, ns
  use mode,           only : krho
  use dist,           only : n_gav_bound
  use constants,      only : ci1 
  use geom,           only : dxgr, sgr_dist
  use general,        only : gkw_abort
  use global, only : r_tiny

  integer, intent(in)    :: imod, ix, i 
  real,    intent(in)    :: delta_zeta, delta_psi, delta_s
  complex, intent(inout) :: gyro_interp_coeffs(-n_gav_bound:nx+n_gav_bound+1) 
  complex, intent(inout), optional :: gyro_interp_coeffs_s(-n_gav_bound:nx+n_gav_bound+1,0:ns+1) 

  integer :: ix_sample, m, stencil_extent
  real    :: ix_sample_xgr_dist, xv, vv 

  integer,parameter         :: ncoeff=6
  complex,dimension(ncoeff) :: coeff

  ! local x-grid location of the sampling point.
  ! if delta_psi is a very small number, this yields ix again...
  ix_sample = floor(ix + delta_psi/dxgr)
  if(ix == ix_sample .and. delta_psi < 0) then
    ! ... then do the floor 'by hand'
    ix_sample = ix - 1
  end if
  

  ! the stencil for higher order interpolation schemes involves
  ! neighbour points after next.
  if (blending_order > 2) then
    stencil_extent = 1
  else
    stencil_extent = 0
  end if
  ! measure, how many grid points beyond the grid boundary are necessary.
  if (ix_sample-stencil_extent < 1) then
    n_gav_bound_measure = max(1-ix_sample + stencil_extent,n_gav_bound_measure)
  endif
  if (ix_sample+stencil_extent > nx) then
    n_gav_bound_measure = max(ix_sample + stencil_extent - nx,n_gav_bound_measure)
  endif
  
  ! this error check is checked in parent routine
  if (n_gav_bound_measure > n_gav_bound) return

  ! ix_sample_xgr_dist is the distance of the sampling point to the corresponding
  ! local x grid point, in terms of grid distance.
  ix_sample_xgr_dist = (delta_psi - dxgr * real(ix_sample - ix))/dxgr

  ! simple check (remove later) 
  if (ix_sample_xgr_dist < 0) then
    write(*,*)imod,ix,i,delta_zeta, delta_psi, ix_sample, ix_sample_xgr_dist, dxgr
    call gkw_abort('ix_sample_xgr_dist < 0. ') 
  end if
  if (abs(delta_s) > sgr_dist) then
    call gkw_abort('rho grad s > sgr_dist : gyro average in s direction fails')
  end if
  
  select case(blending_order) 

  ! second order blending for x direction (linear interpolation) 
  case(2)
    ! the value at that sampling point is interpolated using values of
    ! the 2 neighbouring grid points ix_sample and ix_sample+1.

    ! The field values at these gridpoints for this imod are obtained
    ! in position space representation by multiplying the complex
    ! field with the exp(i k zeta) factor below. Effectively this is
    ! already a (maximum order) interpolation with respect to zeta.
    gyro_interp_coeffs(ix_sample)   = gyro_interp_coeffs(ix_sample) + (1.-ix_sample_xgr_dist) &
       & * exp(ci1*krho(imod)*delta_zeta)
    gyro_interp_coeffs(ix_sample+1) = gyro_interp_coeffs(ix_sample+1) + &
       & ix_sample_xgr_dist * exp(ci1*krho(imod)*delta_zeta)

  ! third order blending for x direction
  case(3) 

    xv = ix_sample_xgr_dist + 1.E0 
    if (xv < 1.5) then 
      vv = (4.E0*xv**2 - 12.E0*xv + 9.E0) / 8.E0  
      gyro_interp_coeffs(ix_sample-1) = gyro_interp_coeffs(ix_sample-1) + vv * exp(ci1*krho(imod)*delta_zeta) 
    endif 
    xv = ix_sample_xgr_dist 
    if (xv <= 0.5) then
      vv = - xv**2 + 3.E0 / 4.E0  
    else 
      vv = (4.E0*xv**2 - 12.E0*xv + 9.E0) / 8.E0  
    endif 
    gyro_interp_coeffs(ix_sample) = gyro_interp_coeffs(ix_sample) + vv * exp(ci1*krho(imod)*delta_zeta) 
    xv = abs(ix_sample_xgr_dist-1.E0)  
    if (xv <= 0.5) then
      vv = - xv**2 + 3.E0 / 4.E0  
    else 
      vv = (4.E0*xv**2 - 12.E0*xv + 9.E0) / 8.E0  
    endif 
    gyro_interp_coeffs(ix_sample+1) = gyro_interp_coeffs(ix_sample+1) + vv * exp(ci1*krho(imod)*delta_zeta) 
    xv = abs(ix_sample_xgr_dist-2.E0)  
    if (xv < 1.5) then 
      vv = (4.E0*xv**2 - 12.E0*xv + 9.E0) / 8.E0  
      gyro_interp_coeffs(ix_sample+2) = gyro_interp_coeffs(ix_sample+2) + vv * exp(ci1*krho(imod)*delta_zeta) 
    endif 

  ! fourth order blending for x direction
  case(4) 

    xv = ix_sample_xgr_dist + 1.E0 
    vv = (2.E0 - xv)**3 / 6.E0   
    gyro_interp_coeffs(ix_sample-1) = gyro_interp_coeffs(ix_sample-1) + vv * exp(ci1*krho(imod)*delta_zeta) 
    xv = ix_sample_xgr_dist 
    vv = (3.E0*xv**3 - 6.E0*xv**2 + 4.E0)/6.E0  
    gyro_interp_coeffs(ix_sample) = gyro_interp_coeffs(ix_sample) + vv * exp(ci1*krho(imod)*delta_zeta) 
    xv = abs(ix_sample_xgr_dist-1.E0)  
    vv = (3.E0*xv**3 - 6.E0*xv**2 + 4.E0)/6.E0  
    gyro_interp_coeffs(ix_sample+1) = gyro_interp_coeffs(ix_sample+1) + vv * exp(ci1*krho(imod)*delta_zeta) 
    xv = abs(ix_sample_xgr_dist-2.E0)  
    vv = (2.E0 - xv)**3 / 6.E0   
    gyro_interp_coeffs(ix_sample+2) = gyro_interp_coeffs(ix_sample+2) + vv * exp(ci1*krho(imod)*delta_zeta) 

  case default 
  
    call gkw_abort('Not implemented blending_order: Options are 2,3,4') 

  end select 

  ! the finite rho* terms are done through linear interpolation in the 
  ! x-direction (i.e. always with second order blending). 
  if (present(gyro_interp_coeffs_s)) then
    coeff(3) = - delta_s / sgr_dist /2.  !ix, i-1
    coeff(4) =   delta_s / sgr_dist /2.   !ix,i+1
    coeff(5) = delta_s/sgr_dist * ix_sample_xgr_dist /2. ! ix +1, i+1
    coeff(4) = coeff(4) -delta_s/sgr_dist * ix_sample_xgr_dist/2. ! ix , i+1
    coeff(6) = -delta_s/sgr_dist * ix_sample_xgr_dist/2. ! ix +1 ,i-1
    coeff(3) = coeff(3)+ delta_s/sgr_dist * ix_sample_xgr_dist /2. ! ix , i -1
    do m=3,ncoeff
      coeff(m)=coeff(m)*exp(ci1*krho(imod)*delta_zeta)
    end do
    gyro_interp_coeffs_s(ix_sample,i-1)   = gyro_interp_coeffs_s(ix_sample,i-1)   + coeff(3)
    gyro_interp_coeffs_s(ix_sample,i+1)   = gyro_interp_coeffs_s(ix_sample,i+1)   + coeff(4)
    gyro_interp_coeffs_s(ix_sample+1,i+1) = gyro_interp_coeffs_s(ix_sample+1,i+1) + coeff(5)
    gyro_interp_coeffs_s(ix_sample+1,i-1) = gyro_interp_coeffs_s(ix_sample+1,i-1) + coeff(6)
  endif

end subroutine calculate_gyro_coef


!------------------------------------------------------------------------------
!> Routine that calculates the contribution to the gyro-averaged field
!> at a single point (imod, ix, i), from the points delta_[xi/ps] away
!> repeatedly adds to the output array gyro_interp_coeffs. 
!> This routine uses a third order polynomal fit. It is at present 
!> not incorporated. The code must be changed to call it. 
!------------------------------------------------------------------------------

subroutine calculate_gyro_coef_higherorder(imod,ix,i,delta_zeta,delta_psi,gyro_interp_coeffs) 

  use grid,       only : nx 
  use mode,       only : krho
  use dist,       only : n_gav_bound
  use constants,  only : ci1 
  use geom,       only : dxgr 
  use general,    only : gkw_abort

  integer, intent(in)    :: imod, ix, i 
  real,    intent(in)    :: delta_zeta, delta_psi
  complex, intent(inout) :: gyro_interp_coeffs(-n_gav_bound:nx+n_gav_bound+1) 

  integer :: ix_sample 
  real    :: ix_sample_xgr_dist

  ! local x-grid location of point of interest
  ix_sample   = floor(ix + delta_psi/dxgr) 

  ! check on ix_sample 
  if ((ix_sample < -n_gav_bound).or.(ix_sample > nx+n_gav_bound+1)) & 
    & call gkw_abort('Index out of range in gyro-average')

  ! distance from local x location
  ix_sample_xgr_dist   = (delta_psi - dxgr * real(ix_sample - ix))/dxgr + 1.0E0  

  gyro_interp_coeffs(ix_sample-1) = gyro_interp_coeffs(ix_sample-1) - &
     & (ix_sample_xgr_dist-1.E0)*(ix_sample_xgr_dist -2.E0)*(ix_sample_xgr_dist-3.E0)* &
     & exp(ci1*krho(imod)*delta_zeta)/ 6.E0
  gyro_interp_coeffs(ix_sample)   = gyro_interp_coeffs(ix_sample)   + &
     & ix_sample_xgr_dist * (ix_sample_xgr_dist-2.E0)*(ix_sample_xgr_dist - 3.E0)*     &
     & exp(ci1*krho(imod)*delta_zeta)/ 2.E0
  gyro_interp_coeffs(ix_sample+1) = gyro_interp_coeffs(ix_sample+1) - &
     & ix_sample_xgr_dist * (ix_sample_xgr_dist-1.E0)*(ix_sample_xgr_dist - 3.E0)*     &
     & exp(ci1*krho(imod)*delta_zeta) / 2.E0
  gyro_interp_coeffs(ix_sample+2) = gyro_interp_coeffs(ix_sample+2) + &
     & ix_sample_xgr_dist * (ix_sample_xgr_dist-1.E0)*(ix_sample_xgr_dist - 2.E0)*     &
     & exp(ci1*krho(imod)*delta_zeta) / 6.E0

  ! The dummy argument 'i' is kept, as it will be needed for finite rho_* terms.
  if (.false.) write(*,*) i

end subroutine calculate_gyro_coef_higherorder


!------------------------------------------------------------------------------
!> This routine sets up the matrix for doing the gyro
!> average of the fields.
!------------------------------------------------------------------------------

subroutine gyro_average_fields_init 

  use control,        only : nlphi, nlapar, nlbpar 
  use grid,           only : nmod, nx, ns, nmu, nsp
  use dist,           only : iphi, iapar, iphi_ga, iapar_ga
  use index_function, only : indx
  use matdat,         only : compress_piece
  use structures,     only : matrix_element
  use matdat,         only : set_indx, connect_rad
  use rho_par_switch, only : s_average  
  use general,        only : gkw_abort
  use matrix_format, only : put_element, finish_matrix
  integer              :: imod, ix, i, j, is, m, nind, nind_ix, nind_i
  integer              :: sind 
  type(matrix_element) :: E 
  logical              :: ingrid 

  if (nlphi) then
    do imod = 1,nmod; do ix = 1,nx; do i=1,ns; do j = 1,nmu; do is = 1,nsp 
    
      ! number of grid points involved in the gyro-average 
      do m = 1, nind_gyro_av(imod,ix,i,j,is)

        ! index of the radial point       
        nind = index_gyro_av(imod,ix,i,j,is,m) 
        
        ! The radial boundary condition  
        call set_indx(E,imod,ix,i,j,1,is) 
        E%ixloc  = nind 
        E%val    = 1. 
        call connect_rad(E,ingrid) 
        nind     = E%ixloc  
        sind     = E%iloc 
        
        ! Only put the element if it is in the grid 
        if (ingrid) then 
          call put_element(matav, &
             & indx(iphi_ga,imod,ix,i,j,is), &
             & indx(iphi,imod,nind,sind), &
             & E%val*coeff_gyro_av(imod,ix,i,j,is,m))
        endif
      end do

      ! rhostar correction to the gyro average.
      ! the gyro average of the potential is not necessary for the polarisation.
      ! I don't have to distinguish between the gyroav in the lowest order and
      ! the rhostar correction and therefore I don't have to reset matav%nmat
      if ( s_average ) then
        do m = 1, nind_gyro_av_s(imod,ix,i,j,is)

          ! index of the radial point       
          nind_ix = index_ix_gyro_av_s(imod,ix,i,j,is,m) 
          ! index of the parallel point       
          nind_i = index_i_gyro_av_s(imod,ix,i,j,is,m) 
          
          ! The boundary condition  
          call set_indx(E,imod,ix,i,j,1,is) 
          E%ixloc  = nind_ix 
          E%iloc  = nind_i
          E%val    = 1. 
          call connect_wrapper(E,ingrid) 
          if (E%iloc < 0 .or. E%iloc > ns) call gkw_abort('error in E%iloc')
            

          ! Only put the element if it is in the grid 
          if (ingrid) then 
            call put_element(matav, indx(iphi_ga,imod,ix,i,j,is), &
               & indx(iphi,imod,E%ixloc,E%iloc), &
               & E%val*coeff_gyro_av_s(imod,ix,i,j,is,m))
          endif 
        enddo
      endif
    end do; end do; end do; end do; end do 
  endif 

  if (nlapar) then  
    do imod = 1,nmod; do ix = 1,nx; do i=1,ns; do j = 1,nmu; do is = 1,nsp 

      ! number of radial grid points involved in the gyro-average  
      do m = 1, nind_gyro_av(imod,ix,i,j,is)
      
        ! index of the radial point 
        nind = index_gyro_av(imod,ix,i,j,is,m) 

        ! apply the boundary condition 
        call set_indx(E,imod,ix,i,j,1,is) 
        E%ixloc  = nind 
        E%val    = 1. 
        call connect_rad(E,ingrid) 
        nind     = E%ixloc  
        sind     = E%iloc 
        
        ! only put the element if it is in the grid 
        if (ingrid) then 
          call put_element(matav, indx(iapar_ga,imod,ix,i,j,is), &
             & indx(iapar,imod,nind,sind), &
             & E%val*coeff_gyro_av(imod,ix,i,j,is,m))
        endif 
 
      end do 

      ! rhostar correction to the gyro average
      if ( s_average) then
        do m = 1, nind_gyro_av_s(imod,ix,i,j,is)

          ! index for the radial point       
          nind_ix = index_ix_gyro_av_s(imod,ix,i,j,is,m) 
          nind_i = index_i_gyro_av_s(imod,ix,i,j,is,m) 
          
          ! The radial boundary condition  
          call set_indx(E,imod,ix,i,j,1,is) 
          E%ixloc  = nind_ix 
          E%iloc   = nind_i
          E%val    = 1. 
          call connect_wrapper(E,ingrid) 
            

          ! Only put the element if it is in the grid 
          if (ingrid) then 
            call put_element(matav, indx(iapar_ga,imod,ix,i,j,is), &
               & indx(iapar,imod,E%ixloc,E%iloc), &
               & E%val*coeff_gyro_av_s(imod,ix,i,j,is,m))
          endif
        enddo
      endif
    end do; end do; end do; end do; end do 
  endif
      
  if (nlbpar) then 
    call gkw_abort('Have not found out how to do nlbpar nonspectral yet')
  endif 

  ! set the integers to initialize the gyro-averaged field to zero 
  phi_ga_start = indx(iphi_ga,1,1,1,1,1) 
  phi_ga_end   = indx(iphi_ga,nmod,nx,ns,nmu,nsp)
  if (nlapar) then
    apar_ga_start = indx(iapar_ga,1,1,1,1,1)
    apar_ga_end = indx(iapar_ga,nmod,nx,ns,nmu,nsp)
  endif
  
  ! compress the matav matrix
  call finish_matrix(matav)

end subroutine gyro_average_fields_init 

!------------------------------------------------------------------------------
!> The field integral equations (similar to poisson_int in linear_terms) 
!------------------------------------------------------------------------------
subroutine field_integral_eqs 

  use control,        only : nlphi, nlapar
  use linear_terms,   only : lampere
  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod, gx 
  use geom,           only : bn 
  use components,     only : signz, vthrat, veta, dgrid
  use dist,           only : ifdis, iphi, iphi_ga, iapar, iapar_ga
  use velocitygrid,   only : intvp, intmu, vpgr 
  use index_function, only : indx
  use matdat,         only : compress_piece 
  use structures,     only : matrix_element
  use matdat,         only : set_indx, connect_rad 
  use rho_par_switch, only : s_average 
  use general,        only : gkw_abort
  use matrix_format, only : put_element, finish_matrix
  
  integer              :: imod, ix, i, j, k, is, m, nind, nind_i, nind_ix
  integer              :: sind 
  type(matrix_element) :: E 
  logical              :: ingrid 

  ! first build the array that calculates the integral over the parallel 
  ! velocity. The integral over vpar is temporarily stored as the gyro-
  ! averaged potential. The ghost cells are required in the next stage.

  ! for the potential 
  if (nlphi) then 
    do imod=1,nmod; do ix=1,nx; do i=1,ns; 
       do j=1,nmu; do k=1,nvpar; do is=1,nsp

        call put_element(matint,indx(iphi_ga,imod,ix,i,j,is),&
           & indx(ifdis,imod,ix,i,j,k,is), cmplx(intvp(i,j,k,is)))

      end do; end do; end do
    end do; end do; end do
  endif 

  ! for the vector potential  
  if (nlapar.and.lampere) then 
    do imod=1,nmod; do ix=1,nx; do i=1,ns; 
       do j=1,nmu; do k=1,nvpar; do is=1,nsp

         call put_element(matint, indx(iapar_ga,imod,ix,i,j,is), &
           & indx(ifdis,imod,ix,i,j,k,is), cmplx(vpgr(i,j,k,is)*intvp(i,j,k,is)))

      end do; end do; end do
    end do; end do; end do 
  endif 

  ! Then build the array that does the gyro-average as well as 
  ! integration over mu and summation of the charge. The solution
  ! is to be stored in a help array of reals (1:2*nmod*nx*ns) 
  ! for the potential
  if (nlphi) then   
    do imod = 1,nmod; do ix = 1, nx; do i = 1, ns 

      do j = 1, nmu; do is = 1, nsp 

        ! number of radial grid points in the gyro-average 
        do m = 1, nind_gyro_cj_Gx(imod,gx(ix),i,j,is)

          ! index of the radial point 
          nind = index_gyro_av(imod,ix,i,j,is,m)
 
          ! apply the boundary condition 
          call set_indx(E,imod,ix,i,j,1,is) 
          E%ixloc  = nind 
          E%val    = 1. 
          call connect_rad(E,ingrid) 
          nind     = E%ixloc  
          sind     = E%iloc 

          ! only put if in grid 
          if (ingrid) then 
            call put_element(matred, 2*(indx(iphi,imod,ix,i)-indx(iphi,1,1,1))+1,&
               & indx(iphi_ga,imod,nind,sind,j,is), &
               & E%val*signz(is)*intmu(j)*bn(ix,i)*dgrid(is)*  &
                             & coeff_gyro_cj_Gx(imod,gx(ix),i,j,is,m))
          end if 
        end do

        if (s_average) then
          do m = 1, nind_gyro_av_s(imod,ix,i,j,is)

            ! index of the radial point       
            nind_ix = index_ix_gyro_av_s(imod,ix,i,j,is,m) 
            ! index of the parallel point       
            nind_i = index_i_gyro_av_s(imod,ix,i,j,is,m) 
            
            ! The boundary condition  
            call set_indx(E,imod,ix,i,j,1,is) 
            E%ixloc  = nind_ix 
            E%iloc  = nind_i
            E%val    = 1. 
            call connect_wrapper(E,ingrid) 
            if (E%iloc < 0 .or. E%iloc > ns) call gkw_abort('error in E%iloc')
              

            ! Only put the element if it is in the grid 
            if (ingrid) then 
              call put_element(matred_s, 2*(indx(iphi,imod,ix,i) - indx(iphi,1,1,1))+1, &
                 & indx(iphi_ga,imod,E%ixloc,E%iloc,j,is), &
                 & E%val*signz(is)*intmu(j)*bn(ix,i)*dgrid(is)*&
                 & coeff_gyro_av_s(imod,ix,i,j,is,m))
            endif
          enddo
        endif


      end do; end do; 

    end do; end do; end do
  endif 
    
  ! for the vector potential 
  if (nlapar.and.lampere) then 
    do imod = 1, nmod; do ix = 1, nx; do i = 1, ns 
  
      do j = 1, nmu; do is = 1, nsp 
  
        ! number of points in the gyro-average 
        do m = 1, nind_gyro_cj_Gx(imod,gx(ix),i,j,is) 
        
          ! radial index of the grid point 
          nind = index_gyro_av(imod,ix,i,j,is,m)

          ! apply the boundary condition 
          call set_indx(E,imod,ix,i,j,1,is) 
          E%ixloc  = nind 
          E%val    = 1. 
          call connect_rad(E,ingrid) 
          nind     = E%ixloc  
          sind     = E%iloc

          ! only put if in grid 
          if (ingrid) then 
            call put_element(matred, 2*(indx(iapar,imod,ix,i)-indx(iphi,1,1,1))+1,&
               & indx(iapar_ga,imod,nind,sind,j,is), &
               & E%val*veta(ix)*signz(is)*intmu(j)*bn(ix,i)*vthrat(is)*  &
                            & dgrid(is)*coeff_gyro_cj_Gx(imod,gx(ix),i,j,is,m))
          end if 
        end do

        if (s_average) then
          do m = 1, nind_gyro_av_s(imod,ix,i,j,is)

            ! index of the radial point       
            nind_ix = index_ix_gyro_av_s(imod,ix,i,j,is,m) 
            ! index of the parallel point       
            nind_i = index_i_gyro_av_s(imod,ix,i,j,is,m) 
            
            ! The boundary condition  
            call set_indx(E,imod,ix,i,j,1,is) 
            E%ixloc  = nind_ix 
            E%iloc  = nind_i
            E%val    = 1. 
            call connect_wrapper(E,ingrid) 
            if (E%iloc < 0 .or. E%iloc > ns) call gkw_abort('error in E%iloc')
              

            ! Only put the element if it is in the grid 
            if (ingrid) then 
              call put_element(matred_s, indx(iapar_ga,imod,ix,i,j,is), &
                 & indx(iapar,imod,E%ixloc,E%iloc), &
                 & E%val*signz(is)*intmu(j)*bn(ix,i)*dgrid(is)*  &
                 & coeff_gyro_av_s(imod,ix,i,j,is,m))
            endif 
          enddo
        endif
      end do; end do; 

    end do; end do; end do
  endif 

  call finish_matrix(matint)
  call finish_matrix(matred)
  call finish_matrix(matred_s)

end subroutine field_integral_eqs 


!------------------------------------------------------------------------------
!> Build the polarization term of the Poisson equation by integrating over 
!> the gyro-average operator and the Maxwellian 
!> Also builds the Ampere's law operator into the same matrix
!> Note the polarisation matrix is constructed global in x, but local in y
!> (that is, after construction, only a subset of modes are put into the 
!> matrix on each proc)
!> In the case of shear periodic boundary conditions (only), this routine
!> may be called repeatedly as a function of time, (indicated by the optional 
!> argument dpart) (this is NOT optimised !) 
!> hence each allocation and deallocation statement must be conditional.
!------------------------------------------------------------------------------

subroutine polarization_init(dpart) 

  use general,         only : gkw_abort
#if defined(umfpack)  
  use control,         only : nlphi, nlapar, time, zonal_adiabatic, shear_periodic
  use linear_terms,    only : lampere
  use control,         only : radial_boundary_conditions
  use grid,            only : ns, nmu, nvpar, nsp, n_x_grid, nmod, lx, lsendrecv_x 
  use grid,            only : parallel_mu, parallel_sp, parallel_vpar, parallel_s
  use grid,            only : n_procs_mu, n_procs_vpar, n_procs_sp, gs
  use velocitygrid,    only : intvp, intmu, vpgr
  use geom,            only : ints, metric_G, jacobian_G, dxgr
  use geom,            only : sgr_dist
  use mode,            only : iyzero, krho 
  use components,      only : iadia, signz, adiabatic_electrons, mas, dgrid, rhostar
  use components,      only : veta_Gx, de_Gx, tmp_Gx, tgrid, types 
  use dist,            only : iphi, iapar, n_phi_start
  use matdat,          only : compress_piece, connect_parallel 
  use global,          only : root_and_verbose, r_tiny, gkw_a_equal_b_accuracy
  use rotation,        only : shear_rate 
  use mpicomms,        only : COMM_S_EQ_X_EQ
  use mpiinterface,    only : mpiallreduce_sum, root_processor
  use constants,       only : ci1, czero
  use structures,      only : matrix_element
  use matdat,          only : set_indx, connect_rad
  use mode,            only : parallel_phase_shift
  use index_function,  only : indx
  use rho_par_switch,  only : s_average
  use matconv,         only : crstoreal_umf
  use global,          only : iumf 
  use general,         only : gkw_warn
  use perform,         only : perfon, perfoff, perffields
  use matrix_format, only : create_matrix, compress_matrix, finalize_matrix
  use matrix_format, only : put_element, matrix_format_gkwcrs

  real, optional, intent(in) :: dpart

  ! This logical determines how the condition of the field equation matrix is 
  ! constrained. If true the integral form is used. If false the potential 
  ! on the first grid point is set to an arbitrary number (determined by the 
  ! charge density on this point) 
  logical :: integration_constraint = .true. 
  logical :: bc_choice = .false.

  integer   :: imod, ix, i, j, k, is, m, l, p, lto, nind, lint, ierr=0
  integer   :: lind_in_x, lind_in_i, nind_out_i, nind_out_x, sind, sint 
  integer   :: nhr, nelem, ig, iv, nmatpol_max, ib, igs 
  integer, save :: icall = 0
  integer(kind=iumf) :: nr, nt  
  real      :: bc_shift, np
  complex   :: phase1, phase2, phase3
  logical   :: ingrid, ingrid1, ingrid2
  type(matrix_element) :: E 
  save :: nelem, nmatpol_max

  ! help arrays for the building of the matrix

  ! help(ix,i,ixref) contains the matrix element connected with the point 
  ! (ix,i) refering to the point (ixref,i). To save memory it is used 
  ! multiple times for each of the toroidal modes. 
  complex, save,      allocatable :: help(:,:,:) 
  ! hilp(ix,i,ixref) contains the matrix element connected with the point 
  ! (ix,i) refering to the point (ixref,isin(i)). To save memory it is used 
  ! multiple times for each of the toroidal modes. The difference with then
  ! help array is that also the i location of the referenced element is 
  ! changed. This can occur in the axis boundary conditions 
  complex, save,      allocatable :: hilp(:,:,:) 
  ! aphe(ix,i,ixref) is the same as help array above, but is used for 
  ! Ampere's law 
  complex, save,      allocatable :: aphe(:,:,:) 
  ! aphi(ix,i,ixref) is the same as the hilp array above, but is used for 
  ! Ampere's law 
  complex, save,      allocatable :: aphi(:,:,:) 
  ! help_za(ix,i,iref) is a help array for the zonal adiabatic correction. 
  ! the latter is the flux surfaced averaged potential correction in the 
  ! adiabatic response. 
  real,    save,      allocatable :: help_za(:,:,:)
  
  complex, save,      allocatable :: help_buf(:,:,:)
  complex, save,      allocatable :: help_s(:,:,:,:), aphe_s(:,:,:,:)
  complex, save,      allocatable :: aprho(:,:,:,:) 
  type(sparse_matrix), save :: mathel
  integer, save,      allocatable :: iachel(:)
  integer(kind=iumf), save,    allocatable :: ap(:), ai(:)
  integer,            save, allocatable :: ixprh(:,:), ixmrh(:,:), iiprh(:,:), iimrh(:,:)
  real(kind=rumf)   , save, allocatable :: ax(:)
  real,               save, allocatable :: F(:,:,:,:,:) ! factor in the coefficients of the laplacian
  real,               save, allocatable :: maxwellint(:,:,:,:), maxwvp2int(:,:,:,:) 

  real    :: jacmetl1, jacmetl2, jacmetp1, jacmetp2
  logical :: apar_rho = .false. 
      
  first: if (.not. present(dpart)) then
      
    ! allocate the help array 
    if (nlphi) then 
      allocate(maxwellint(n_x_grid,ns,nmu,nsp),stat = ierr) 
      if (ierr /= 0) call gkw_abort('could not allocate maxwellint in gyro-average')
      allocate(help(1:n_x_grid,ns,1:n_x_grid), stat = ierr) 
      if (ierr /= 0) call gkw_abort('could not allocate help in gyro-average')
      allocate(hilp(1:n_x_grid,ns,1:n_x_grid), stat = ierr) 
      if (ierr /= 0) call gkw_abort('could not allocate help in gyro-average')
    endif 
    if (nlapar) then 
      allocate(maxwvp2int(n_x_grid,ns,nmu,nsp),stat = ierr) 
      if (ierr /= 0) call gkw_abort('could not allocate maxwvp2int in gyro-average')
      allocate(aphe(1:n_x_grid,ns,1:n_x_grid), stat = ierr) 
      if (ierr /= 0) call gkw_abort('could not allocate aphe in gyro-average')
      allocate(aphi(1:n_x_grid,ns,1:n_x_grid), stat = ierr) 
      if (ierr /= 0) call gkw_abort('could not allocate aphe in gyro-average')
    endif 
    if (apar_rho) then 
      if (parallel_s) call gkw_abort('Finite rho* for Ampere does not work with parallel_s')
      allocate(aprho(n_x_grid,ns,3,3), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate aprho in gyro-average')
      allocate(ixprh(n_x_grid,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate ixprh in gyro-average')
      allocate(ixmrh(n_x_grid,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate ixmrh in gyro-average')
      allocate(iiprh(n_x_grid,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate iiprh in gyro-average')
      allocate(iimrh(n_x_grid,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate iimrh in gyro-average')
    endif 

    ! allocate the matrix
    if (integration_constraint) then 
      nelem = nmod_l*n_x_grid*ns*(4*n_gav_bound_measure+1) + n_x_grid*ns*ns  &
            & + n_x_grid**2*ns
      if (nlapar) nelem = nelem + nmod_l*n_x_grid*ns*(4*n_gav_bound_measure+ &
            & 1) + n_x_grid**2*ns
    else   
      nelem = nmod_l*n_x_grid*ns*(4*n_gav_bound_measure+1) + n_x_grid*ns*ns 
      if (nlapar) nelem = nelem + nmod_l*n_x_grid*ns*(4*n_gav_bound_measure+1)
    endif 

    if (.not. allocated(mathel%mat)) then
      mathel = create_matrix("polarization helper matrix", nelem, &
         & matrix_format_gkwcrs)
    end if

    if ( s_average ) then
      ierr = 0
      if (nlphi) allocate (help_s(n_x_grid,ns,n_x_grid,ns), stat = ierr)
      if (ierr /= 0) call gkw_abort('could not allocate help_s in gyro-average')
      ! is not completly implemented  jet
      if (.false.) allocate (aphe_s(n_x_grid,ns,n_x_grid,ns), stat = ierr)
      if (ierr /= 0) call gkw_abort('could not allocate aphe_s in gyro-average')
      nmatpol_max=nmod*n_x_grid*ns*nsp*n_x_grid*n_s_ga
      if (nlapar) nelem=nelem*2
    endif
  
    ! allocate the help array 
    if (zonal_adiabatic) allocate(help_za(1:n_x_grid,ns,ns), stat = ierr) 
    if (ierr /= 0) call gkw_abort('could not allocate help_za in gyro-average')
        
    bc_shift = 0.
  
  else
    
    bc_shift = -(time + dpart) * shear_rate * lx
    
  end if first

  if(radial_boundary_conditions=='periodic'.or.radial_boundary_conditions &
             &  =='Neuslab') bc_choice = .true.

  if (zonal_adiabatic .and. parallel_s .and. (.not. (iyzero == 0))) then
    call gkw_abort('non-spectral zonal_adiabatic not yet with parallel s')
  end if

  if (perffields) call perfon('Build mathel',2)

  ! set the counter for the mathel matrix to zero
  mathel%nmat = 0

  ! the Maxwell integrated over the parallel velocity. This is used
  ! in both the polarization as well as the Maxwell correction in 
  ! Ampere's law. 
  if (nlphi) then 
    maxwellint(:,:,:,:) = 0. 
    do ix = 1, n_x_grid; do i = 1, ns; do j = 1, nmu; do is = 1, nsp 
      do k = 1, nvpar 
        if (types(is) == 'EP') then
          maxwellint(ix,i,j,is) = maxwellint(ix,i,j,is) +  &
            &        intvp(i,j,k,is)*fEP_Gx(ix,i,j,k,is)
        else
          maxwellint(ix,i,j,is) = maxwellint(ix,i,j,is) +  & 
           &        intvp(i,j,k,is)*fmaxwl_Gx(ix,i,j,k,is)
        endif
      end do 
    end do; end do; end do; end do 
  endif 
  if (nlapar) then 
    maxwvp2int(:,:,:,:) = 0. 
    do ix = 1, n_x_grid; do i = 1, ns; do j = 1, nmu; do is = 1, nsp 
      do k = 1, nvpar 
        if (types(is) == 'EP') then
          maxwvp2int(ix,i,j,is) = maxwvp2int(ix,i,j,is) +  &
           & 2.E0*vpgr(i,j,k,is)**2*intvp(i,j,k,is)*fEP_Gx(ix,i,j,k,is)
        else 
          maxwvp2int(ix,i,j,is) = maxwvp2int(ix,i,j,is) +  & 
           & 2.E0*vpgr(i,j,k,is)**2*intvp(i,j,k,is)*fmaxwl_Gx(ix,i,j,k,is)
        endif
      end do 
    end do; end do; end do; end do 
  endif 
   
  ! one big outer loop for the toroidal modes (this is to save memory)
  do imod = 1, nmod ! imod loop 

    ! Initialize the help arrays 
    if (nlphi) then 
      help(:,:,:)     = (0.E0,0.E0)
      hilp(:,:,:)     = (0.E0,0.E0) 
    endif 
    if (nlapar) then 
      aphe(:,:,:)     = (0.E0,0.E0)
      aphi(:,:,:)     = (0.E0,0.E0) 
    endif 
    if (zonal_adiabatic) help_za(:,:,:)  = 0.E0
    if (apar_rho)        aprho(:,:,:,:)  = (0.E0,0.E0) 
    if (s_average) then
      if (nlphi)         help_s(:,:,:,:) = (0.E0,0.E0)   
      if (nlapar)        aphe_s(:,:,:,:) = (0.E0,0.E0)   
    endif
 
    !-------------------------------------------------------------------------
    ! The double gyro-average that occurs in the Poisson and Ampere's
    ! law equations. In the Poisson equation this represent the 
    ! part of the polarization, in Ampere's law the Maxwell correction
    !-------------------------------------------------------------------------
 
    ! All grid points of the Poisson / Ampere's law equation.
    ! The integral over the parallel velocity is performed before 
    ! the matrix solve 
    do ix = 1, n_x_grid; do i = 1,ns; do j = 1, nmu; do is = 1, nsp 

      ! the number of points used for the gyro-average at each 
      ! grid point. The loop presents the outer gyro-average. 

      ! outer gyro-average  
      do m = 1, nind_gyro_cj_Gx(imod,ix,i,j,is) 

        ! Set the matrix element (to diagonal) 
        call set_indx(E,imod,ix,i,j,1,is)
        ! find the ix-point to which the gyro-average refers. 
        E%ixloc  = index_gyro_cj_Gx(imod,ix,i,j,is,m)
        ! Set the value of the element to 1.0. After caling the 
        ! boundary condition routine the value could contain a 
        ! phase shift or be zero (Dirichlet)
        E%val    = 1.0 
        ! apply the boundary conditions. This can lead to a remap of 
        ! the ix point, and a change in the element E%val. When the 
        ! proper boundary conditions on the axis are used also a 
        ! remap of the s-point is possible 
        call connect_rad(E,ingrid,bc_shift) 
        nind   = E%ixloc ! ix index of the referenced point 
        sind   = E%iloc  ! i  index of the referenced point 
        phase1 = E%val   ! multiplication factor of the boundary cond. 

        ! Then do the inner gyro-average. lto is the number of points
        ! of the inner gyro-average 
        lto = nind_gyro_av_Gx(imod,nind,sind,j,is) 
        do l = 1, lto ! Inner gyro-average
 
          ! find the ix point and apply the boundary conditions 
          ! Set the element (to the gyro-centre referenced by the outer
          ! gyro-average) 
          call set_indx(E,imod,nind,sind,j,1,is)
          ! Find the ix point referenced by the inner gyro-average 
          E%ixloc = index_gyro_av_Gx(imod,nind,sind,j,is,l)
          E%val   = 1.0
          ! apply the boundary conditions 
          call connect_rad(E,ingrid,bc_shift) 
          lint    = E%ixloc
          sint    = E%iloc 
          phase2  = E%val 
          
          ! the phase defined below is the product of the phase shifts of the 
          ! two boundary conditions. Note that the boundary can be crossed 
          ! twice. The exception is the Dirichlet boundary condition, for which 
          ! the product of two boundary crossing would be zero, but the referenced
          ! point after two crossing can lie in the computational domain. In this 
          ! case it is considered in the double gyro-average. 
          if (.not. gkw_a_equal_b_accuracy(phase1, czero)) phase2 = phase1 * phase2

          if (nlphi.and.ingrid) then
            ! ib is set to ix at the moment. In the future it might be possible 
            ! to shift it to the point referenced by the outer gyro-average. 
            ! This has perhaps the correct long wave length limit (with the 
            ! density inside the divergence) but this needs further study
            ib = min(max(1,nind),n_x_grid) ; ib = ix 
            ! In case of the boundary condition on the axis also the i-point 
            ! can be shifted. This is then stored in a separate array. 
            if (sint /= i) then
              ! check (removed later) 
              if (sint /= isin(i)) call gkw_abort('inconsistent isin') 
              hilp(ix,i,lint) = hilp(ix,i,lint) - phase2*bn_Gx(ix,i)*             &
           &  intmu(j)*dgrid(is)*signz(is)**2*maxwellint(ib,i,j,is)/tmp_Gx(ib,is) &
           &  * coeff_gyro_av_Gx(imod,nind,sind,j,is,l) & ! inner gyro-average G
           &  * coeff_gyro_cj_Gx(imod,ix,i,j,is,m)        ! outer gyro-average G^dag
            else
              help(ix,i,lint) = help(ix,i,lint) - phase2*bn_Gx(ix,i)*             &
           &  intmu(j)*dgrid(is)*signz(is)**2*maxwellint(ib,i,j,is)/tmp_Gx(ib,is) &
           &  * coeff_gyro_av_Gx(imod,nind,sind,j,is,l) & ! inner gyro-average G
           &  * coeff_gyro_cj_Gx(imod,ix,i,j,is,m)        ! outer gyro-average G^dag
            endif 
          endif 
          if (nlapar.and.ingrid.and.lampere) then 
            ! ib is set to ix at the moment. In the future it might be possible 
            ! to shift it to the point referenced by the outer gyro-average. 
            ! This has perhaps the correct long wave length limit (with the 
            ! density inside the divergence) but this needs further study
            ib = min(max(1,nind),n_x_grid) ; ib = ix 
            if (sint /= i) then 
              ! check (removed later) 
              if (sint /= isin(i)) call gkw_abort('inconsistent isin') 
              aphi(ix,i,lint) = aphi(ix,i,lint) + phase2*veta_Gx(ib)*intmu(j)*    &
           &  signz(is)**2*dgrid(is)*bn_Gx(ix,i)*maxwvp2int(ib,i,j,is)/(mas(is)   &
           &  * tmp_Gx(ib,is)/tgrid(is))                                          &
           &  * coeff_gyro_av_Gx(imod,nind,sind,j,is,l) & ! inner gyro-average G
           &  * coeff_gyro_cj_Gx(imod,ix,i,j,is,m)        ! outer gyro-average G^dag
            else 
              aphe(ix,i,lint) = aphe(ix,i,lint) + phase2*veta_Gx(ib)*intmu(j)*    &
           &  signz(is)**2*dgrid(is)*bn_Gx(ix,i)*maxwvp2int(ib,i,j,is)/(mas(is)   &
           &  * tmp_Gx(ib,is)/tgrid(is))                                          &
           &  * coeff_gyro_av_Gx(imod,nind,sind,j,is,l) & ! inner gyro-average G
           &  * coeff_gyro_cj_Gx(imod,ix,i,j,is,m)        ! outer gyro-average G^dag
            endif 
          endif 

        end do ! inner gyro-average 
      end do ! outer gyro-average 

    end do; end do; end do; end do ! loop over all grid points. 

    !-------------------------------------------------------------------------
    ! The part of the polarization that is directly proportional (i.e. 
    ! without gyro-average) to the potential. 
    !-------------------------------------------------------------------------
    if (nlphi) then
      do ix = 1, n_x_grid; do i = 1,ns; do j = 1, nmu; do is = 1, nsp
        help(ix,i,ix) = help(ix,i,ix) + maxwellint(ix,i,j,is)*               &
                   & bn_Gx(ix,i)*intmu(j)*signz(is)**2*dgrid(is)/tmp_Gx(ix,is)

        ! If integration_contraint reduce the condition number for zonal mode
        if (imod == iyzero .and. bc_choice .and. integration_constraint) then
          do l = 1, n_x_grid
            help(ix,i,l) = help(ix,i,l) - ints(i)*dxgr
          end do
        end if
      end do; end do; end do; end do ! loop over all grid points. 
    endif

    !-------------------------------------------------------------------------
    ! The zonal adiabatic correction. This makes the problem nonlocal 
    ! in the s-direction. 
    !-------------------------------------------------------------------------
    if (zonal_adiabatic .and. imod == iyzero) then
      do ix = 1, n_x_grid; do i = 1,ns; do p = 1, ns    
        help_za(ix,i,p) = help_za(ix,i,p) - ints(i)*de_Gx(ix,nsp+iadia)     &
                        & /tmp_Gx(ix,nsp+iadia)
      end do; end do; end do 
    end if 
    
    !-------------------------------------------------------------------------
    ! Finite rho* corrections due to the parallel derivatives. 
    ! This part is usually neglected as it has very little influence 
    ! on heat and particle fluxes. The finite rho* parallel derivatives
    ! can generate a momentum flux though
    !-------------------------------------------------------------------------
    if(s_average) then 
      do ix = 1, n_x_grid; do i = 1,ns; do j = 1, nmu; do is = 1, nsp 

        ! first the build the rho* correction for then outer gyro av and then
        ! the rho* correction for the inner gyro av. the rhostar correction
        ! for the inner and outer gyro av \sim rho**2 and neglected.
        
        ! Outer rhostar correction to gyroaverage
        do m = 1, nind_gyro_cj_s_Gx(imod,ix,i,j,is)

          call set_indx(E,imod,ix,i,j,1,is)
          E%ixloc  = index_ix_gyro_av_s_Gx(imod,ix,i,j,is,m)
          E%iloc   =  index_i_gyro_av_s_Gx(imod,ix,i,j,is,m)
          E%val    = 1.0 
          call connect_wrapper(E,ingrid,bc_shift) 
          ! outer index 
          nind_out_x  = E%ixloc 
          nind_out_i   = E%iloc 
          phase1   = E%val 

          lto = nind_gyro_av_Gx(imod,nind_out_x,nind_out_i,j,is) 
          do l = 1, lto ! Inner gyro-average in the lowest order
   
            call set_indx(E,imod,ix,i,j,1,is) 
            E%ixloc = index_gyro_av_Gx(imod,nind_out_x,nind_out_i,j,is,l)
            if (.not. gkw_a_equal_b_accuracy(phase1, czero)) then
              E%val = phase1
            else 
              E%val   = 1.0
            endif 
            call connect_wrapper(E,ingrid,bc_shift) 
            ! inner index; the gyro average in the lowest order refers to the
            ! same i point
            lind_in_x   = E%ixloc 
            lind_in_i   = nind_out_i
            phase2     = E%val 
            !check inner and outer coeff
            if (nlphi.and.ingrid) then
              help_s(ix,i,lind_in_x,lind_in_i) = help_s(ix,i,lind_in_x,lind_in_i) - &
              phase2*bn_Gx(ix,i)*intmu(j)*maxwellint(ix,i,j,is)*&
             &  coeff_gyro_av_Gx(imod,nind_out_x,nind_out_i,j,is,l)*coeff_gyro_cj_s_Gx(imod,ix,i,j,is,m)  &
             &  * dgrid(is)*signz(is)**2/tmp_Gx(ix,is) 
            endif 
!              if (nlapar.and.ingrid.and.lampere) then 
!                aphe_s(ix,i,lint) = aphe_s(ix,i,lint) + phase2*veta_Gx(ix)*signz(is)**2*       & 
!                  & dgrid(is)*bn_Gx(ix,i)*intmu(j)*maxwellint/mas(is)                      &
!                  & *coeff_gyro_av_s_Gx(imod,nind,i,j,is,l)*coeff_gyro_cj_s_Gx(imod,ix,i,j,is,m)
!              endif 

          end do 
        end do !outer

        ! now the inner gyro av is the rho* correction
        
        ! Outer gyroaverage in the lowest order
        do m = 1, nind_gyro_cj_Gx(imod,ix,i,j,is)

          call set_indx(E,imod,ix,i,j,1,is)
          E%ixloc  = index_gyro_av_Gx(imod,ix,i,j,is,m)
          E%val    = 1.0 
          call connect_rad(E,ingrid,bc_shift) 
          nind_out_x  = E%ixloc 
          nind_out_i  = i
          phase1   = E%val

          lto = nind_gyro_av_s_Gx(imod,nind,i,j,is) 
          do l = 1, lto ! Inner rhostar correction gyro-average 
   
            call set_indx(E,imod,ix,i,j,1,is) 
            E%ixloc = index_ix_gyro_av_s_Gx(imod,nind_out_x,nind_out_i,j,is,l)
            E%iloc  = index_i_gyro_av_s_Gx(imod,nind_out_x,nind_out_i,j,is,l)
            if (.not. gkw_a_equal_b_accuracy(phase1, czero)) then
              E%val = phase1
            else 
              E%val   = 1.0
            endif 
            call connect_wrapper(E,ingrid,bc_shift) 
            lind_in_x    = E%ixloc 
            lind_in_i    = E%iloc 
            phase2     = E%val 
            
            if (nlphi.and.ingrid) then
              help_s(ix,i,lind_in_x,lind_in_i) = help_s(ix,i,lind_in_x,lind_in_i) &
             &      - phase2*bn_Gx(ix,i)*intmu(j)*maxwellint(ix,i,j,is)* &
             &  coeff_gyro_av_Gx(imod,nind_out_x,nind_out_i,j,is,l)*coeff_gyro_cj_s_Gx(imod,ix,i,j,is,m)  &
             &  * dgrid(is)*signz(is)**2/tmp_Gx(ix,is) 
            endif 

          end do 
        end do !outer

      end do; end do; end do; end do 
    endif !s average
        
    !-------------------------------------------------------------------------
    ! The nabla^2 term in Ampere's law 
    !-------------------------------------------------------------------------
    if (nlapar.and.lampere) then  

      do ix = 1, n_x_grid; do i = 1,ns 
    
        ! Workaround to undo the reduction (could be removed with some reordering)
        np = n_procs_mu * n_procs_vpar * n_procs_sp

        ig = gs(i) 

        call set_indx(E,imod,ix,i,1,1,1) 
        E%val   = 1.0 
        E%ixloc = ix - 1
        call connect_rad(E,ingrid1,bc_shift) 
        l       = E%ixloc 
        phase1  = E%val
        ! only the ix-1 point can be remapped to refer to another i position when 
        ! axis boundaries are used. 
        igs     = E%iloc 
        
        ! the lower middle grid point Lagrangian (no igs correction at present)
        if (igs /= i) then 
          jacmetl1 = - jacobian_G(m)*metric_G(m,ig,1,1) 
          jacmetl2 = - jacobian_G(m)*metric_G(m,ig,1,2) 
        else 
          if (l >= 1) then 
            jacmetl1 = jacobian_G(l)*metric_G(l,ig,1,1) 
            jacmetl2 = jacobian_G(l)*metric_G(l,ig,1,2) 
          else 
            jacmetl1 = 2*jacobian_G(1)*metric_G(1,ig,1,1) -  &
                     &   jacobian_G(2)*metric_G(2,ig,1,1)
            jacmetl2 = 2*jacobian_G(1)*metric_G(1,ig,1,2) -  &
                     &   jacobian_G(2)*metric_G(2,ig,1,2)
          endif  
        endif 
          
        call set_indx(E,imod,ix,i,1,1,1) 
        E%val   = 1.0 
        E%ixloc = ix 
        call connect_rad(E,ingrid,bc_shift) 
        m       = E%ixloc 
        phase2  = E%val
        

        call set_indx(E,imod,ix,i,1,1,1) 
        E%val   = 1.0 
        E%ixloc = ix + 1
        call connect_rad(E,ingrid2,bc_shift) 
        p       = E%ixloc 
        phase3  = E%val  
        if (p <= n_x_grid) then 
          jacmetp1 = jacobian_G(p)*metric_G(p,ig,1,1) 
          jacmetp2 = jacobian_G(p)*metric_G(p,ig,1,2) 
        else 
          jacmetp1 = 2*jacobian_G(n_x_grid)*metric_G(n_x_grid,ig,1,1) -  &
                   &   jacobian_G(n_x_grid-1)*metric_G(n_x_grid-1,ig,1,1)
          jacmetp2 = 2*jacobian_G(n_x_grid)*metric_G(n_x_grid,ig,1,2) -  &
                   &   jacobian_G(n_x_grid-1)*metric_G(n_x_grid-1,ig,1,2)
        endif  
 
 !      jacmetl1 = jacobian_G(m)*metric_G(m,ig,1,1)
 !      jacmetl2 = jacobian_G(m)*metric_G(m,ig,1,2)
!       jacmetp1 = jacobian_G(m)*metric_G(m,ig,1,1)
!       jacmetp2 = jacobian_G(m)*metric_G(m,ig,1,2)

        ! 1/J (d/dpsi)[ J g^11 d/dpsi] 
        if (ingrid1) then 
          if (igs == i) then 
            aphe(ix,i,l) = aphe(ix,i,l) - 0.5/(jacobian_G(m)*dxgr**2) * phase1 * (          &
              & jacmetl1 + jacobian_G(m)*metric_G(m,ig,1,1)) /np    
          else 
            ! check (removed later) 
            if (igs /= isin(i)) call gkw_abort('inconsistent isin') 
            aphi(ix,i,l) = aphi(ix,i,l) - 0.5/(jacobian_G(m)*dxgr**2) * phase1 * (          &
              & jacmetl1 + jacobian_G(m)*metric_G(m,ig,1,1)) /np    
          endif           
        endif 

        aphe(ix,i,m) = aphe(ix,i,m) +0.5/(jacobian_G(m)*dxgr**2) * phase2 * (               &
          & jacmetl1 + 2*jacobian_G(m)*metric_G(m,ig,1,1) + jacmetp1) / np                                          
  
        if (ingrid2) then 
          aphe(ix,i,p) = aphe(ix,i,p) - 0.5/(jacobian_G(m)*dxgr**2) * phase3 * (              &
            & jacobian_G(m)*metric_G(m,ig,1,1) + jacmetp1) / np    
        endif 


        ! 1/J (d/dpsi) [ J g^12 i k_zeta ]  + i k_zeta g^21 d/dpsi
        if (ingrid1) then 
          if (igs == i) then 
            aphe(ix,i,l) = aphe(ix,i,l) + 0.5/(jacobian_G(m)*dxgr)*ci1*krho(imod)*phase1*     &
                       & (jacobian_G(m)*metric_G(m,ig,2,1) + jacmetl2)/np     
          else 
            aphi(ix,i,l) = aphi(ix,i,l) + 0.5/(jacobian_G(m)*dxgr)*ci1*krho(imod)*phase1*     &
                       & (jacobian_G(m)*metric_G(m,ig,2,1) + jacmetl2)/np     
          endif 
        endif 
        if (ingrid2) then 
          aphe(ix,i,p) = aphe(ix,i,p) - 0.5/(jacobian_G(m)*dxgr)*ci1*krho(imod)*phase3*      &
                       & (jacobian_G(m)*metric_G(m,ig,2,1)+jacmetp2)/np 
        endif 

        ! - k_zeta^2 g^22
        aphe(ix,i,m) = aphe(ix,i,m) + phase2*metric_G(m,ig,2,2)*krho(imod)**2/np  

        ! finite rho* effects in apar (necessary for tearing modes) 
        if (apar_rho) then 


          ixprh(ix,i) = p 
          ixmrh(ix,i) = l 

          call set_indx(E,imod,ix,i,1,1,1)
          E%iloc = i + 1 
          E%val  = parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
          call connect_parallel(E,ingrid)
          iiprh(ix,i) = E%iloc 
          if (ingrid) then 
            iv     = E%iloc 
            aprho(ix,i,3,3) = aprho(ix,i,3,3) -E%val*phase3*rhostar*(jacobian_G(p)*          &
                   & metric_G(p,i,1,3)+jacobian_G(m)*metric_G(m,iv,3,1)) /                   &
                   & (4*np*jacobian_G(m)*sgr_dist*dxgr)
            aprho(ix,i,1,3) =  aprho(ix,i,1,3) + E%val*phase1*rhostar*(jacobian_G(l)*        &
                   & metric_G(l,i,1,3)+jacobian_G(m)*metric_G(m,iv,3,1)) /                   &
                   & (4*np*jacobian_G(m)*sgr_dist*dxgr)
            aprho(ix,i,2,3) = aprho(ix,i,2,3) - E%val*phase2*rhostar*(ci1*krho(imod)*        &
                   & (metric_G(m,i,2,3) +metric_G(m,iv,3,2))/(2*sgr_dist) + rhostar*         &
                   & (metric_G(m,i,3,3) + metric_G(m,iv,3,3))/(2*sgr_dist**2)) /np 
            aprho(ix,i,2,2)=aprho(ix,i,2,2)+rhostar**2*(metric_G(ix,iv,3,3)+                 &
                   & metric_G(ix,i,3,3))/(2.E0*np*sgr_dist**2)
          endif 
          call set_indx(E,imod,ix,i,1,1,1)
          E%iloc = i - 1 
          E%val  = parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
          call connect_parallel(E,ingrid)
          iimrh(ix,i) = E%iloc 
          if (ingrid) then 
            iv     = E%iloc 
            aprho(ix,i,3,1) =  aprho(ix,i,3,1) + E%val*phase3*rhostar*(jacobian_G(p)*        &
                   & metric_G(p,i,1,3)+jacobian_G(m)*metric_G(m,iv,3,1))                     &
                   & / (4*np*jacobian_G(m)*sgr_dist*dxgr) 
            aprho(ix,i,1,1) = aprho(ix,i,1,1)-E%val*phase1*rhostar*(jacobian_G(l)*           &
                   & metric_G(l,i,1,3)+jacobian_G(m)*metric_G(m,iv,3,1)) /                   &
                   & (4*np*jacobian_G(m)*sgr_dist*dxgr)
            aprho(ix,i,2,1) =  aprho(ix,i,2,1) + E%val*phase2*rhostar*(ci1*krho(imod)*       &
                   & (metric_G(m,i,2,3)+metric_G(m,iv,3,2))/(2*sgr_dist) - rhostar*          &
                   & (metric_G(m,i,3,3) +metric_G(m,iv,3,3))/(2*sgr_dist**2))/np 
            aprho(ix,i,2,2)=aprho(ix,i,2,2)+rhostar**2*(metric_G(ix,iv,3,3)+                 &
                   & metric_G(ix,i,3,3))/(2.E0*np*sgr_dist**2)
          endif 


        endif 

      end do; end do     
    endif 
  
    if (perffields) call perfoff(2)
    if (perffields) call perfon('mathel reduce, convert, put',2)
  
    !-------------------------------------------------------------------------
    ! The mpiallreduce_sum for the case of parallelization in the velocity 
    ! space or species. 
    !-------------------------------------------------------------------------
    if (parallel_mu .or. parallel_sp .or. parallel_vpar) then       

      ! allocate the help buffer 
      if (.not. allocated(help_buf)) then
        allocate(help_buf(n_x_grid,ns,1:n_x_grid), stat = ierr) 
        if (ierr /= 0) call gkw_abort('Can not allocate help_buf in gyro_average')
      end if

      if (nlphi) then 

        ! add the elements 
        call mpiallreduce_sum(help,help_buf,n_x_grid,ns,n_x_grid,COMM_S_EQ_X_EQ)

        ! copy back (could be written as an implicit loop!)
        do ix = 1, n_x_grid; do i = 1, ns; do m = 1, n_x_grid
          help(ix,i,m) = help_buf(ix,i,m) 
        end do; end do; end do

        ! add the elements 
        call mpiallreduce_sum(hilp,help_buf,n_x_grid,ns,n_x_grid,COMM_S_EQ_X_EQ)

        ! copy back (could be written as an implicit loop!)
        do ix = 1, n_x_grid; do i = 1, ns; do m = 1, n_x_grid
          hilp(ix,i,m) = help_buf(ix,i,m) 
        end do; end do; end do

        ! when a remapping of the i-position occurs through the axis boundary 
        ! it must always have the same i value independent of imod,ix,j,k,is
        ! However, not all j points must have a Larmor radius large enough for 
        ! to encircle the axis and so remapping must not occur for all grid 
        ! points. 
        !call mpiallreduce_max(isin,isin_buf,ns,COMM_S_EQ_X_EQ) 
        
        ! copy back 
        !do i = 1, ns 
        !  isin(i) = isin_buf(i) 
        !end do 
        
      endif 
 
      if (nlapar.and.lampere) then 

        ! add the elements 
        call mpiallreduce_sum(aphe,help_buf,n_x_grid,ns,n_x_grid,COMM_S_EQ_X_EQ)

        ! copy back (could be written as an implicit loop!)
        do ix = 1, n_x_grid; do i = 1, ns; do m = 1, n_x_grid
          aphe(ix,i,m) = help_buf(ix,i,m) 
        end do; end do; end do

        ! add the elements 
        call mpiallreduce_sum(aphi,help_buf,n_x_grid,ns,n_x_grid,COMM_S_EQ_X_EQ)

        ! copy back (could be written as an implicit loop!)
        do ix = 1, n_x_grid; do i = 1, ns; do m = 1, n_x_grid
          aphi(ix,i,m) = help_buf(ix,i,m) 
        end do; end do; end do

        ! find isin 
        !call mpiallreduce_max(isin,isin_buf,ns,COMM_S_EQ_X_EQ) 
        
        ! copy back 
        !do i = 1, ns 
        !  isin(i) = isin_buf(i) 
        !end do 

      endif 

      ! deallocate the buffer arrays  
      if (.not. shear_periodic) then 
        deallocate(help_buf)
      endif 

    endif 

    !-------------------------------------------------------------------------
    ! The part below is for an approximation of the polarization used, for 
    ! instance, in ORB. This is not the best physics model, but useful for 
    ! benchmarks
    !-------------------------------------------------------------------------
    if (orb_polarize) then 

      if(root_processor)then
        write(*,*)'WARNING: calculating the polarization using'
        write(*,*)'1/jac * div ( (sum_s m_s n_s / 2 B^2 * jac * metric) nabla phi)' 
      endif
      
      ! allocate the F array 
      allocate(F(0:n_x_grid+1,1:ns,1:nsp,3,3), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate the F in gyro-average')
   
      ! Does not work when parallelizing over species 
      if (n_procs_sp .ne. 1) stop 'No parallelization over species' 

      ! first undo the polarization calculation 
      help(:,:,:) = (0.E0,0.E0) 

      ! calculation of 1/jac * m_s n_s / 2 B^2 * jac * metric
      do i = 1, ns; do is = 1, nsp; do l = 1, 3; do m = 1, 3
        ig = gs(i) 

        do ix = 1, n_x_grid 
          F(ix,i,is,l,m) = 0.5*(mas(is)*de_Gx(ix,is)/bn_Gx(ix,i)**2)*  &
                         & jacobian_G(ix)*metric_G(ix,ig,l,m)
        end do 
        F(0,i,is,l,m)        = 2*F(1,i,is,l,m) - F(2,i,is,l,m) 
        F(n_x_grid+1,i,is,l,m) = 2*F(n_x_grid,i,is,l,m) - F(n_x_grid-1,i,is,l,m) 
 
      end do; end do; end do; end do
      

      do ix = 1, n_x_grid; do i = 1, ns 

        do is = 1, nsp 

          ! new laplacian
          ig = gs(i)

          call set_indx(E,imod,ix,i,1,1,1) 
          E%val   = 1. 
          E%ixloc = ix -1
          call connect_rad(E,ingrid1,bc_shift) 
          l       = E%ixloc 
          igs     = E%iloc 
          phase1  = E%val  

          call set_indx(E,imod,ix,i,1,1,1) 
          E%val   = 1. 
          E%ixloc = ix 
          call connect_rad(E,ingrid,bc_shift) 
          m       = E%ixloc 
          phase2  = E%val  

          call set_indx(E,imod,ix,i,1,1,1) 
          E%val   = 1. 
          E%ixloc = ix + 1
          call connect_rad(E,ingrid2,bc_shift) 
          p       = E%ixloc 
          phase3  = E%val  

          if (ingrid1) then
            if (igs /= i) then 
              hilp(ix,i,l) = hilp(ix,i,l) - phase1*0.5/(dxgr*jacobian_G(ix))*(           &
                 & ( F(l,i,is,1,1)+F(m,i,is,1,1) )/dxgr -                                &
                 & 0.5*ci1*krho(imod)*( F(l,i,is,1,2)+F(m,i,is,1,2) ) -                  &
                 & ci1*krho(imod)*F(m,i,is,2,1) )
            else 
              help(ix,i,l) = help(ix,i,l) - phase1*0.5/(dxgr*jacobian_G(ix))*(           &
                 & ( F(l,i,is,1,1)+F(m,i,is,1,1) )/dxgr -                                &
                 & 0.5*ci1*krho(imod)*( F(l,i,is,1,2)+F(m,i,is,1,2) ) -                  &
                 & ci1*krho(imod)*F(m,i,is,2,1) )
            endif 
          endif 

          help(ix,i,m) = help(ix,i,m) - phase2*(1.E0/jacobian_G(ix))*(               &
             & -0.5*( F(p,i,is,1,1) + 2.E0*F(m,i,is,1,1) +F(l,i,is,1,1) )/dxgr**2 +  &
             & 0.25*ci1*krho(imod)*( F(p,i,is,1,2)-F(l,i,is,1,2) ) / dxgr -          &
             & krho(imod)**2*F(m,i,is,2,2) )

          if (ingrid2) then 
            help(ix,i,p) = help(ix,i,p) - phase3*0.5/(dxgr*jacobian_G(ix))*(           &
               & ( F(p,i,is,1,1)+F(m,i,is,1,1) )/dxgr +                                &
               & 0.5*ci1*krho(imod)*( F(p,i,is,1,2)+F(m,i,is,1,2) ) +                  &
               & ci1*krho(imod)*F(m,i,is,2,1) )
          endif      

        end do 
      end do; end do  
    
      if (.not. shear_periodic) deallocate(F)

    endif 

    ! add the adiabatic electron correction (if present) 
    if (nlphi .and. adiabatic_electrons) then 
      do ix = 1,n_x_grid; do i = 1,ns 
        help(ix,i,ix) = help(ix,i,ix) - signz(nsp+iadia)*de_Gx(ix,nsp+iadia)/tmp_Gx(ix,nsp+iadia)
      end do; end do 
    endif   
  
    if (nmod_l < 1) then  ! no modes on this processor, nothing more to do
      ihavemodes = .false.
      continue ! must stay in loop to avoid MPI deadlock
    else
      ihavemodes = .true.
    end if

    !open(43,file = 'effe') 
    !ix = 1 
    !do nind =1, n_x_grid
    !  write(43,500)(help(ix,i,nind),i = 1,ns), (hilp(ix,i,nind), i = 1, ns) 
    !end do 
    !500 format(128(1pe13.5,1x))
    !close(43) 
      
    ! Now build the matrix (perhaps the two steps can be combined to save memory)
    ! Here only a subset of modes are kept
    if (nlphi) then 

      if ((imod >= first_imod).and.(imod <= last_imod)) then 

        do ix = 1, n_x_grid; do i = 1, ns 

          if (imod == iyzero .and. ix == 1 .and. (.not. integration_constraint)) then 

            call put_element(mathel,&
               & indx_gx(iphi,imod,ix,i),indx_gx(iphi,imod,ix,i),cmplx(1.0,0.0))

          else 
       
            do m = 1, n_x_grid 
              if (abs(help(ix,i,m)) > r_tiny) then
                call put_element(mathel,&
                   & indx_gx(iphi,imod,ix,i), indx_gx(iphi,imod,m,i),&
                   & help(ix,i,m))
              endif 
            end do 
   
            do m = 1, n_x_grid 
              if (abs(hilp(ix,i,m)) > r_tiny) then
                call put_element(mathel,&
                   & indx_gx(iphi,imod,ix,i),indx_gx(iphi,imod,m,isin(i)), &
                   & hilp(ix,i,m))
                if ((isin(i) < 1) .or. (isin(i) > ns)) call gkw_abort('error isin')
              endif 
            end do 

            if (imod == iyzero .and. zonal_adiabatic) then
              do p = 1, ns 
                if (abs(help_za(ix,i,p)) > r_tiny) then
                  call put_element(mathel, &
                     & indx_gx(iphi,imod,ix,i), indx_gx(iphi,imod,ix,p), &
                     & cmplx(help_za(ix,i,p)))
                endif 
              end do 
            end if 
            

            !check s_average here (if not integration constraint) or alway
            if(s_average) then
              do m=1,n_x_grid; do p=1,ns
                if(abs(help_s(ix,i,m,p)) > r_tiny) then
                  ! matpol%ii is requiered to build z in field solve non spec. but z
                  ! is real, but not complex.
                  ! how do I get the local imod
                  if ( lsendrecv_x )then
                    ! matpol%ii has to be global in imod and x, not local.
                    call gkw_abort('matpol%ii contains the wrong index')
                  endif
                  ! matpol%jj accesses either phi_dum or fdis(iphi(...)).
                  if (lsendrecv_x) then
                    call put_element(matpol, &
                       & 2*(indx(iphi,imod,ix,i)-indx(iphi,1,1,1))+1, &
                       & indx_gx(iphi,imod,m,p), &
                       & help_s(ix,i,m,p))
                  else
                    call put_element(matpol, &
                       & 2*(indx(iphi,imod,ix,i)-indx(iphi,1,1,1))+1, &
                       & indx_gx(iphi,imod,m,p)+n_phi_start-1, &
                       & help_s(ix,i,m,p))
                  endif
                endif

              enddo;enddo
            endif!s_average

          endif !integration_contsraint
        end do; end do 
      
      endif 
    endif 

    if (nlapar.and.lampere) then 

      if ((imod >= first_imod).and.(imod <= last_imod)) then 

        do ix = 1, n_x_grid; do i = 1, ns 
   
          do m = 1, n_x_grid 
            if (abs(aphe(ix,i,m)) > r_tiny) then
              call put_element(mathel, &
                 & indx_gx(iapar,imod,ix,i), indx_gx(iapar,imod,m,i), &
                 & aphe(ix,i,m))
            endif 
          end do 

          do m = 1, n_x_grid 
            if (abs(aphi(ix,i,m)) > r_tiny) then
              call put_element(mathel, &
                 & indx_gx(iapar,imod,ix,i), indx_gx(iapar,imod,m,isin(i)), &
                 & aphi(ix,i,m))
            endif 
          end do 

        end do; end do 

       ! add the finite rho* terms 
       if (apar_rho) then
     
         do ix = 1, n_x_grid; do i = 1, ns
           call put_element(mathel, &
              & indx_gx(iapar,imod,ix,i), &
              & indx_gx(iapar,imod,ixprh(ix,i),iiprh(ix,i)), &
              & aprho(ix,i,3,3))
           call put_element(mathel, &
              & indx_gx(iapar,imod,ix,i), &
              & indx_gx(iapar,imod,ixmrh(ix,i),iiprh(ix,i)), &
              & aprho(ix,i,1,3))
           call put_element(mathel, &
              & indx_gx(iapar,imod,ix,i), &
              & indx_gx(iapar,imod,ixmrh(ix,i),iimrh(ix,i)), &
              & aprho(ix,i,1,1))
           call put_element(mathel, &
              & indx_gx(iapar,imod,ix,i), &
              & indx_gx(iapar,imod,ixprh(ix,i),iimrh(ix,i)), &
              & aprho(ix,i,3,1))


           call put_element(mathel, &
              & indx_gx(iapar,imod,ix,i), &
              & indx_gx(iapar,imod,ix,iiprh(ix,i)), &
              & aprho(ix,i,2,3))
           call put_element(mathel, &
              & indx_gx(iapar,imod,ix,i), &
              & indx_gx(iapar,imod,ix,i), &
              & aprho(ix,i,2,2))
           call put_element(mathel, &
              & indx_gx(iapar,imod,ix,i), &
              & indx_gx(iapar,imod,ix,iimrh(ix,i)), &
              & aprho(ix,i,2,1))

         end do; end do 

       endif 

      endif 
    endif 

  end do ! loop over toroidal modes 
  
  ! no longer needed   
  if(.not. shear_periodic) then
    if (allocated(bn_Gx))            deallocate(bn_Gx)
    if (allocated(fmaxwl_Gx))        deallocate(fmaxwl_Gx)
    if (allocated(fEP_Gx))           deallocate(fEP_Gx)
    if (allocated(coeff_gyro_av_Gx)) deallocate(coeff_gyro_av_Gx)
    if (allocated(index_gyro_av_Gx)) deallocate(index_gyro_av_Gx)
    if (allocated(help))             deallocate(help)
    if (allocated(help_za))          deallocate(help_za)
    if (allocated(aphe))             deallocate(aphe) 
    if (allocated(aprho))            deallocate(aprho)
    if (allocated(ixprh))            deallocate(ixprh)
    if (allocated(ixmrh))            deallocate(ixmrh)
    if (allocated(iiprh))            deallocate(iiprh)
    if (allocated(iimrh))            deallocate(iimrh)       
    if (allocated(maxwellint))       deallocate(maxwellint)
  end if

  ! sort the mathel array (on columns first, for umfpack)
  call compress_matrix(mathel, on_columns=.true.)

  ! The number of "rows" 
  nhr=maxval(mathel%jj(1:mathel%nmat)) - minval(mathel%jj(1:mathel%nmat)) + 1
  !write(*,*)  maxval(mathel%jj(1:mathel%nmat)), minval(mathel%jj(1:mathel%nmat)) 
      
  ! set up an ia array
  if (.not. allocated(iachel)) then
     allocate(iachel(nhr + 1), stat = ierr)
     if (ierr /= 0) call gkw_abort('could not allocate iachel in gyro-average')
  end if

  !integer array that points to the beginning of the next "row"
  iachel(1) = 1
  i = 1
  do j = 1, nhr-1
    do while (mathel%jj(i) == j)
      i = i + 1
      if (i>mathel%nmat) call gkw_abort('internal problem in gyro_average') 
    end do
    iachel(j+1) = i
  end do
  iachel(nhr+1) = mathel%nmat + 1

  ! Set the array values for the real matrix 
  nr  = 2*nhr
  nt  = 4*mathel%nmat

  ! allocate the real arrays
  if (.not. allocated(ap)) then
    allocate(ap(nr+1),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate ap in gyro_average')
    allocate(ai(nt),stat = ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate ai in gyro_average')
    allocate(ax(nt),stat = ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate ax in gyro_average')
  end if
  
  ap = 0; ai = 0; ax = 0

  ! convert this array to a real array in UMF precision
  call crstoreal_umf(nhr,iachel,mathel,nt,ap,ai,ax)

  ! Complex matrix no longer needed
  if (.not. shear_periodic) then
    call finalize_matrix(mathel)
    deallocate(iachel)
  end if

  ! set up the indices for umfpack 
  ! convert to 0 base (because of the use of C inside UMFPACK) 
  do i = 1, int(nr+1)
    ap(i) = ap(i) - 1
  end do 
  do i = 1, int(ap(nr+1))  ! no -1 here. just subtracted
    ai(i) = ai(i) - 1
  end do 

  if (perffields) call perfoff(2)
  if (perffields) call perfon('mathel umfpack',2)
    
  ! set default parameters
  call umf4def (control_umf)

  control_umf(1) = 2
  if (root_and_verbose) call umf4pcon (control_umf)

  ! pre-order and symbolic analysis
  if (ihavemodes) then
  
    call umf4sym (nr, nr, Ap, Ai, Ax, symbolic, control_umf, info)

    if (root_and_verbose) then 
      ! print statistics computed so far
      ! call umf4pinf (control, info) could also be done.
      write(*,80) info (1), info (16), (info (21) * info (4)) / 2**20, &
        &      (info (22) * info (4)) / 2**20,info (23), info (24), info (25)
 80      format ('symbolic analysis:',/,                            &
      &      '   status:  ', f5.0, /,                               &
      &      '   time:    ', e10.2, ' (sec)'/,                      &
      &      '   estimates (upper bound) for numeric LU:', /,       &
      &      '   size of LU:    ', f10.2, ' (MB)', /,               &
      &      '   memory needed: ', f10.2, ' (MB)', /,               &
      &      '   flop count:    ', e10.2, /                         &
      &      '   nnz (L):       ', f10.0, /                         &
      &      '   nnz (U):       ', f10.0)

    endif 
      
    ! check umf4sym error condition
    if (info (1) < 0) then
      print *, 'Error occurred in umf4sym: ', info (1)
      call gkw_abort('Problem in gyro_average')
    endif

    ! numeric factorization
    call umf4num (Ap, Ai, Ax, symbolic, numeric, control_umf, info)

    if (root_and_verbose) then 
      ! print statistics for the numeric factorization
      ! call umf4pinf (control, info) could also be done.
      write(*,90) info (1), info (66), (info (41) * info (4)) / 2**20,   &
      &      (info (42) * info (4)) / 2**20, info (43), info (44), info (45)
  90      format ('numeric factorization:',/,                       &
      &      '   status:  ', f5.0, /,                              &
      &      '   time:    ', e10.2, /,                             &
      &      '   actual numeric LU statistics:', /,                &
      &      '   size of LU:    ', f10.2, ' (MB)', /,              &
      &      '   memory needed: ', f10.2, ' (MB)', /,              &
      &      '   flop count:    ', e10.2, /                        &
      &      '   nnz (L):       ', f10.0, /                        &
      &      '   nnz (U):       ', f10.0)
    endif 
    
    if (1./info(68) > 1e6) call gkw_warn('Ill conditioned polarisation matrix')
    
    if (root_processor .and. icall == 0) then
      write(*,*)
      write(*,'(" Condition number of polarisation matrix: ", es10.3)') 1./info(68)
      write(*,*)    
    end if
    icall = icall + 1
      
    ! check umf4num error condition
    if (info (1) < 0) then
      print *, 'Error occurred in umf4num: ', info (1)
      call gkw_abort('Problem with umfpack in gyro_average')
    endif

  end if !ihavemodes

  if (perffields) call perfoff(2)

  ! no longer need the real matrix
  if (.not. shear_periodic) then 
    deallocate(ap)
    deallocate(ai)
    deallocate(ax)
  end if
  
#else
  real, optional, intent(in) :: dpart
  call gkw_abort('Non spectral presently requires compilation with umfpack')
#endif

  contains 

  !--------------------------------------------------------------------
  ! Integer function connected with the boundary conditions. Of course,
  ! this should somehow be moved into the boundary conditions, but until
  ! proven that the mapping works it is here
  !---------------------------------------------------------------------
  integer function isin(i)

    use grid, only : ns
    integer, intent(in) :: i

    isin = i - ns/2
    if (isin < 1) isin = isin + ns

    return
  end function isin

end subroutine polarization_init 

!----------------------------------------------------------------------------
!> This routine sets-up indexing arrays for putting the toroidal modes
!> back into the right global location after the fields solve is complete  
!> Likely the same could be accomplished with MPI derived datatypes,
!> and might be faster. As with indx_gx, this also relies on matching index order
!----------------------------------------------------------------------------
subroutine setup_indx_nmod

  use index_function, only : indx
  use dist,           only : iphi, iapar
  use grid,           only : nmod, ns, nx
  use control,        only : nlapar
  use matdat,         only : compress_piece
  use general,        only : gkw_abort
  real, allocatable :: dummy(:)
  integer :: imod, ix, i, ic2, n_apar

  ic2 =0   
  n_phi     = indx(iphi,1,1,1) - 1
  n_apar    = indx(iapar,1,1,1) - 1
  
  do ix = 1, nx;
    do i = 1, ns;
      do imod = first_imod, last_imod
        ic2 = ic2 + 1

        ! These arrays translate z,
        ! from (local in x, global in y) into (local in x, local in y)
        iiscatz(ic2) = 2*(indx(iphi, imod, ix, i) - n_phi) - 1
        iilocz(ic2)  = 2*((imod-first_imod+1 -1) + (i-1)*nmod_l + (ix-1)*ns*nmod_l + 1 ) - 1

        ! safety check
        if (nmod_l == nmod .and. iiscatz(ic2) .ne. iilocz(ic2)) then
          call gkw_abort('Error in setup_indx_nmod')
        end if

      end do;
    end do;
  end do;
 
  if (nlapar) then
    do ix = 1, nx; do i = 1, ns; do imod = first_imod, last_imod 
      ic2 = ic2 + 1
      
      ! These arrays translate z, 
      ! from (local in x, global in y) into (local in x, local in y)
      iiscatz(ic2) = 2*(indx(iapar, imod, ix, i) - n_phi) - 1
      iilocz(ic2)  = 2*((imod-first_imod+1 -1) + (i-1)*nmod_l + (ix-1)*ns*nmod_l + 1 ) - 1 
      iilocz(ic2)  = iilocz(ic2) + 2*nmod_l*ns*nx
               
      ! safety check 
      if (nmod_l == nmod .and. iiscatz(ic2) .ne. iilocz(ic2)) then
        call gkw_abort('Error in setup_indx_nmod')
      end if  
      
    end do; end do; end do; 
  end if
      
  nscatz = ic2

  allocate(dummy(1:size(iilocz)))
  ! sort for cache efficiency
  call compress_piece(1,nscatz,iilocz,iiscatz,dummy)
  deallocate(dummy)
  
end subroutine setup_indx_nmod 

!--------------------------------------------------------------
!> for debugging parallel (remove later)
!--------------------------------------------------------------
subroutine write_fdis(fdis,filelun,z,insw)

  use index_function,  only : indx
  use grid,            only : nmod, nx, ns, nmu, nvpar, nsp, n_x_grid
  use grid,            only : gx, gs
  use dist,            only : iphi, ifdis, n_phi_start, iapar_ga
  use mpiinterface,    only : processor_number
  
  real, optional,  intent(in) :: z(:)
  integer, optional,  intent(in) :: insw
  complex, intent(inout) :: fdis(:)
  integer, intent(in) :: filelun
  
  integer :: ix, j, k, is, imod, i
  
  if (present(z)) then
  
    if (present(insw)) then

      do ix = 1,n_x_grid; do i = 1,ns; do imod=first_imod,last_imod 
        write(filelun+processor_number,*)  imod, gs(i), ix, z(2*(indx_gx(iphi,imod,ix,i))-1)
      end do; end do; end do;
    
    else
      
      do ix = 1,nx; do i = 1,ns; do imod=1,nmod
        write(filelun+processor_number,*)  imod, gs(i), gx(ix), z(2*(indx(iphi,imod,ix,i)-n_phi_start + 1)-1)
      end do; end do; end do;
      
    end if

  else
      
    if (present(insw)) then
      
      do ix = 1,n_x_grid; do i = 1,ns; do imod=first_imod,last_imod 
        write(filelun+processor_number,*)  imod, gs(i), ix, abs(fdis(indx_gx(iphi,imod,ix,i)))
      end do; end do; end do;
    
    else
          
      do ix = 1,nx; do i = 1,ns; do imod=1,nmod 
        write(filelun+processor_number,*)  imod, gs(i), gx(ix), abs(fdis(indx(iphi,imod,ix,i)))
      end do; end do; end do;
    
    end if
        
    do j=nmu,nmu; do k=nvpar,nvpar; do is= 1,nsp
    do imod=1,nmod; do i = 1,ns; do ix = 1,nx
      write(filelun+50+processor_number,*)  imod, j, k, is, gs(i), gx(ix), abs(fdis(indx(ifdis,imod,ix,i,j,k,is)))
    end do; end do; end do
    end do; end do; end do
    
    do j=1,nmu; do is= 1,nsp
    do imod=1,nmod; do i = 1,ns; do ix = 1,nx
      write(filelun+20+processor_number,*)  imod, j, is, gs(i), gx(ix), abs(fdis(indx(iapar_ga,imod,ix,i,j,is)))
    end do; end do; end do
    end do; end do
    
  end if
  
end subroutine write_fdis

!----------------------------------------------------------------------------
!> Mini index function for phi/apar and z (global in the x grid)
!> but local in the nmod grid !
!> Requires consistency with the main index order 
!> The version here is chosen to be easiest to gather without defining
!> a derived datatype, this should be fixed.
!> Should also be adapted to generalise with index_order 
!> and moved into the index function
!> When not parallel in x or nmod, is equivalent to
!> indx(iphi,imod,ix,1) - indx(iphi,1,1,1) + 1 for the potential and 
!> indx(iapar,imod,ix,1) - indx(iapar,1,1,1) + 1 + nmod*n_x_grid*ns
!----------------------------------------------------------------------------
function indx_gx(ifield,imod,ix,i)

  use grid,           only : ns, n_x_grid, nx, nmod
  use dist,           only : n_phi_start, iphi, iapar 
  use index_function, only : indx
  use general,        only : gkw_abort

  integer, intent(in) :: ifield, imod, ix, i
  integer :: indx_gx
  
  !easiest for gathering over x only without defining new derived datatypes
  indx_gx = ((imod-first_imod+1)-1) + (i-1)*nmod_l + (ix-1)*ns*nmod_l +1 
  if (ifield == iapar) then   
    indx_gx = indx_gx + nmod_l*ns*n_x_grid 
  endif 

  !agrees with normal index default order
  !!indx_gx = (imod-1) + (ix-1)*nmod + (i-1)*n_x_grid*nmod +1 

  ! safety check only for the potential 
  if ((n_x_grid == nx).and.(ifield == iphi).and.nmod == nmod_l) then
    if (indx(iphi,imod,ix,i)-n_phi_start+1 .ne. indx_gx) then
      call gkw_abort('Error in indx_phi_gx')
    end if
  end if
  
  if (ix < 1 .or. ix > n_x_grid) call gkw_abort('indx_gx does not handle ghost cells')

end function


!----------------------------------------------------------------------------
!> This routine creates global arrays in x needed for the polarisation solve
!> The values returned by index_gyro_av_Gx are also made global
!----------------------------------------------------------------------------
subroutine gather_global_x_arrays

  use grid,           only : nmod, ns, nx, n_x_grid, nsp, nvpar, nmu, gx
  use mpicomms,       only : COMM_X_NE, COMM_DUMMY
  use mpiinterface,   only : gather_array
  use components,     only : energetic_particles
  use geom,           only : bn 
  use dist,           only : fmaxwl, f_EP 
  use rho_par_switch, only : s_average
  use general,        only : gkw_abort

  integer :: ierr, i, j, k, imod, l, ix, is, nind_max
  integer, allocatable :: indx_dum(:,:)
  
  ! allocate the local helper
  ierr=0
  allocate(indx_dum(nx,nsp), stat = ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate indx_dum in gyro-average')  

  ! initialize 
  nind_gyro_av_Gx  = 0 
  index_gyro_av_Gx = 0
  coeff_gyro_av_Gx = 0. 

  ! gather all arrays that will be needed global in x
                      
  call gather_array(bn_Gx(1:n_x_grid,1:ns), n_x_grid,  ns,               &
                  & bn   (1:nx,      1:ns), nx,        ns,               &
                  & COMM_X_NE, COMM_DUMMY, ALLGATHER = ALL_PROCS)

  do imod = 1,nmod; do i = 1, ns; do j = 1, nmu
    call gather_array(nind_gyro_av_Gx(imod,1:n_x_grid,i,j,:),n_x_grid,nsp,&
                    & nind_gyro_av (imod,1:nx,   i,j,:)   ,nx,      nsp,  &
                    & COMM_X_NE, COMM_DUMMY,  ALLGATHER = ALL_PROCS)

    do k=1,nvpar
      call gather_array(fmaxwl_Gx(1:n_x_grid,i,j,k,1:nsp), n_x_grid, nsp,   &
                      & fmaxwl   (1:nx,      i,j,k,1:nsp), nx,       nsp,   &
                      & COMM_X_NE, COMM_DUMMY, ALLGATHER = ALL_PROCS)
      if (energetic_particles) then
        call gather_array(fEP_Gx(1:n_x_grid,i,j,k,1:nsp), n_x_grid, nsp,   &
                      & f_EP   (1:nx,      i,j,k,1:nsp), nx,       nsp,   &
                      & COMM_X_NE, COMM_DUMMY, ALLGATHER = ALL_PROCS)
      end if
    end do; !nvpar

    do l = 1, n_points_ga+4

      ! translate the x-indices to global x-indices
      do is=1,nsp; do ix=1,nx
        if (l <= nind_gyro_av(imod,ix,   i,j,is)) then
          indx_dum(ix,is)=gx(index_gyro_av(imod,ix,i,j,is,l))
        else
          ! These values should not matter, but better be sure than sorry.
          indx_dum(ix,is)= -100 !index_gyro_av(imod,ix,i,j,is,l)
        end if
      end do; end do

      call gather_array(index_gyro_av_Gx(imod,1:n_x_grid,i,j,:,l),n_x_grid,nsp, &
                      &         indx_dum(1:nx,1:nsp)             ,nx,      nsp, &
                      & COMM_X_NE, COMM_DUMMY,  ALLGATHER = ALL_PROCS)

      call gather_array(coeff_gyro_av_Gx(imod,1:n_x_grid,i,j,:,l),n_x_grid,nsp, &
                      & coeff_gyro_av   (imod,1:nx,      i,j,:,l),nx,      nsp, &
                      & COMM_X_NE, COMM_DUMMY,  ALLGATHER = ALL_PROCS)

    end do ! l

    if (s_average) then 
      nind_max=0
      do ix=1,nx; do is=1,nsp
        nind_max=max(nind_gyro_av_s(imod,ix,i,j,is),nind_max)
      enddo;enddo
      do l=1,nind_max
        ! gather the length of the of the coeff array
        call gather_array(nind_gyro_av_s_Gx(imod,1:n_x_grid,i,j,:),n_x_grid,nsp,&
                        & nind_gyro_av_s (imod,1:nx,   i,j,:)   ,nx,      nsp,  &
                        & COMM_X_NE, COMM_DUMMY,  ALLGATHER = ALL_PROCS)

        ! gather the i indices for the coeff. I stay local in i
        do is=1,nsp; do ix=1,nx
          indx_dum(ix,is)=index_i_gyro_av_s(imod,ix,i,j,is,l)
        end do; end do
        call gather_array(index_i_gyro_av_s_Gx(imod,1:n_x_grid,i,j,:,l),n_x_grid,nsp, &
                        &         indx_dum(1:nx,1:nsp)             ,nx,      nsp, &
                        & COMM_X_NE, COMM_DUMMY,  ALLGATHER = ALL_PROCS)
                        
        ! gather the ix indices for the coeff. I translate to global ix
        do ix=1,nx; do is=1,nsp
          indx_dum(ix,is)=gx(index_ix_gyro_av_s(imod,ix,i,j,is,l))
        end do; end do
        call gather_array(index_ix_gyro_av_s_Gx(imod,1:n_x_grid,i,j,:,l),n_x_grid,nsp, &
                        &         indx_dum(1:nx,1:nsp)             ,nx,      nsp, &
                        & COMM_X_NE, COMM_DUMMY,  ALLGATHER = ALL_PROCS)
        ! and the coeff
        call gather_array(coeff_gyro_av_s_Gx(imod,1:n_x_grid,i,j,:,l),n_x_grid,nsp, &
                        & coeff_gyro_av_s   (imod,1:nx,      i,j,:,l),nx,      nsp, &
                        & COMM_X_NE, COMM_DUMMY,  ALLGATHER = ALL_PROCS)
      enddo
    endif

  end do; end do; end do !nmod, ns, nmu 

  ! deallocate the helping arrays. 
  if (allocated(indx_dum))   deallocate(indx_dum)

end subroutine gather_global_x_arrays

!------------------------------------------------------------------------------
!> This subroutine handles the parallel boundary condition for the gyro average
!> in the s direction. it successively calls connect_rad and connect_parallel
!------------------------------------------------------------------------------
subroutine connect_wrapper(E,ingrid,bc_shift)

  use grid,       only : nx, ns
  use structures, only : matrix_element
  use matdat,     only : connect_rad, connect_parallel 
  use general,    only : gkw_abort

  type (matrix_element), intent(inout)  :: E
  logical ,intent(inout)                :: ingrid
  real,optional, intent(in)             :: bc_shift 


  
  ! if not determined by the connect_(parallel|rad) the point lies in the grid.
  ingrid =.true.
  if(E%iloc <= 0 .and. (E%ixloc <= 0 .or. E%ixloc>nx) ) then
    ! the i direction should be only off the grid by 1
    E%iloc=E%iloc + 1
    if(E%iloc<=0) call gkw_abort( "error 1 in connect_wrapper")
    if(present (bc_shift)) then 
      call connect_rad(E,ingrid,bc_shift)
    else
      call connect_rad(E,ingrid)
    endif
    E%iloc=E%iloc - 1
    call connect_parallel(E,ingrid)
  else if(E%iloc > ns .and. (E%ixloc<=0.or.E%ixloc > nx) ) then
    E%iloc=E%iloc -1 
    if(E%iloc>ns) call gkw_abort("error 3 in connect_wrapper")
    if(present (bc_shift)) then 
      call connect_rad(E,ingrid,bc_shift)
    else
      call connect_rad(E,ingrid)
    endif
    E%iloc=E%iloc + 1
    call connect_parallel(E,ingrid)
  else
    if(present (bc_shift)) then 
      call connect_rad(E,ingrid,bc_shift)
    else
      call connect_rad(E,ingrid)
    endif
    call connect_parallel(E,ingrid)
  endif
end subroutine connect_wrapper

!------------------------------------------------------------------------------
!> Exit / warn with a useful error message if x ghost cells sizes are incorrect 
!> All processors must call with the same values
!------------------------------------------------------------------------------
subroutine xgrid_check(xgc_need, xgc_have,reason)
 
  use grid,    only : n_procs_x, n_x_grid  
  use general, only : max_factor, gkw_warn, gkw_exit
  use global, only : int2char

  character(len=*)  :: reason
  integer, intent(in) :: xgc_need, xgc_have
  integer :: id

  id = xgc_need - xgc_have
  
  ! Too many, could limit the parallelism
  if (id < 0) then
    if (max_factor(n_x_grid,n_x_grid / xgc_need) >  & 
       &  max_factor(n_x_grid,n_x_grid / xgc_have)) then
      call gkw_warn('Decrease n_gav_bound_ex in GYROAVERAGE by '//           &
                & trim(int2char(abs(id),0))//char(10)//                      &
                & ' to increase maximum n_procs_x to '//                     &
                & trim(int2char(max_factor(n_x_grid,n_x_grid / xgc_need),0)))
 
    end if
  ! Too few, can't run
  else if (id > 0) then
    if (max_factor(n_x_grid,n_x_grid / xgc_need) < n_procs_x) then
      call gkw_exit('Insufficient ' //trim(reason)// ' for gyro-average:'//  &
                & char(10)//                                                 & 
                & ' Increase n_gav_bound_ex in GYROAVERAGE namelist by '//   &
                & trim(int2char(id,0))//char(10)//' and reduce n_procs_x to '// &
                & trim(int2char(max_factor(n_x_grid,n_x_grid / xgc_need),0)))    
    else
      call gkw_exit('Insufficient ' //trim(reason)// ' for gyro-average:'//  &
                & char(10)//                                                 & 
                & ' Increase n_gav_bound_ex in GYROAVERAGE namelist by '//   &
                & trim(int2char(id,0)))
    end if              
  end if

end subroutine xgrid_check

end module gyro_average
