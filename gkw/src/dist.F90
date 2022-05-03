!-----------------------------------------------------------------------------
!> This module contains the code solution and sets up how the different parts
!> are arranged in the main array. The index function module is responsible
!> for the mapping between points in the computational domain and the array
!> index and so needs to be initialised accordingly. Communicators used for
!> sending parts of the local array to other parts of the same array on other
!> processors are also initialised from here.
!-----------------------------------------------------------------------------

module dist
  use index_function, only : IS_3D_FIELD, IS_COLL_CONS_FIELD, IS_GYROAVG_FIELD
  use index_function, only : IS_6D_DISTR
  implicit none

  private

  !
  ! public subroutines and functions
  !

  public :: dist_init
  public :: get_apar, get_phi, get_bpar

  !
  ! public parameters
  !

  ! Identifiers to access certain parts of fdisi using the index function.
  ! These identifiers incorporate information on the dimensionality of the
  ! respective data (testable using bitwise-and).
  ! Check with the index function to see how the are used.
  ! Changes will need to be made in index_function if anything here is changed.
  integer, public, parameter :: iphi     = 1 + IS_3D_FIELD      !< potential
  integer, public, parameter :: iphi_ga  = 1 + IS_GYROAVG_FIELD !< gyro-averaged potential 
  integer, public, parameter :: iapar    = 2 + IS_3D_FIELD      !< apar
  integer, public, parameter :: iapar_ga = 2 + IS_GYROAVG_FIELD !< gyro-averaged vector potential 
  integer, public, parameter :: ibpar    = 3 + IS_3D_FIELD      !< bpar
  integer, public, parameter :: ibpar_ga = 3 + IS_GYROAVG_FIELD !< gyro-avared parallel magn. field 
  integer, public, parameter :: i_mom    = 4 + IS_COLL_CONS_FIELD !< momentum change
  integer, public, parameter :: i_ene    = 5 + IS_COLL_CONS_FIELD !< energy change
  integer, public, parameter :: ifdis    = 6 + IS_6D_DISTR      !< distribution function


  !
  ! public variables
  !

  !> Number of fields
  integer, public, save :: number_of_fields
  
  !> A counter, to predict the number of elements in the matrix of the linear terms
  !> (or at least get a very good upper estimate)
  integer, public, save :: ntot
  
  !> The distribution function (NOT including the potential) has a
  !> total of nf elements per processor.
  integer, public, save :: nf
  
  !> The distribution function (including the potential) has a total
  !> of nsolc elements to solve for per processor.
  ! nsolc is the total number of elements of the distribution and the
  ! fields, INCLUDING gyroavg fields and INCLUDING the collisionop-related fields
  integer, public, save :: nsolc
  
  !> The distribution (including the potential) has a total of
  !> msolc ( = nsolc + the ghostcells ) elements in memory.
  integer, public, save :: msolc
  
  !> The matrix for the neoclassics diagnostics has nelem_nc elements
  integer, public, save :: nelem_nc
  
  !> The matrix for the collisions conservations has nelem_cc elements
  integer, public, save :: nelem_cc

  !> The matrix for the g to f conversion has nelem_g2f (= nf) elements
  integer, public, save :: nelem_g2f
  
  !> location of the last point of `regular' fields in fdisi
  integer, public, save :: nregular_fields_end
  
  !> the total number of elements in the collisions conservation fields
  integer, public, save :: nelem_conserve

  !> the offset for the collisions conservation feilds (start of both)
  integer, public, save :: n_conserve

  !> number of radial grid points beyond the boundary in gyroaverage 
  integer, public, save :: n_gav_bound   
  
  !> extra number of radial grid points beyond the boundary in gyroaverage, 
  !> which can be used if estimate fails
  !> this is read and broadcast in gyroaverage namelist
  integer, public, save :: n_gav_bound_ex = 0

  !> the maxwell fmaxwl(nx,ns,nmu,nvp,nsp)
  real, public, allocatable, save :: fmaxwl(:,:,:,:,:)

  !> the EP f_EP(nx,ns,nmu,nvp,nsp) and derivatives
  real, public, allocatable, save :: f_EP(:,:,:,:,:)
  real, public, allocatable, save :: df_EPdv(:,:,:,:,:)
  real, public, allocatable, save :: df_EPdv_sinhc(:,:,:,:,:)

  !> the slowing down alpha particle distribution
  real, public, allocatable, save :: falpha(:,:,:)
  
  !> the distribution g as labeled in docs (NOT f) is stored in fdisi(nsolc)
  complex, public, allocatable, save :: fdisi(:)
  
  !> the time averaged distribution g used if laverage_dist_over_time
  complex, public, allocatable, save :: fdisi_tavg(:)
  
  !> and temporary copy of fdisi with MPI ghost cells 
  !> which can at various points contain either f or g (see exp_integration)
  !> fdis_tmp(msolc)
  complex, public, allocatable, save :: fdis_tmp(:), fdis_tmp2(:)

  !> the distribution function and rhs of the intermediate Runge-Kutta steps
  complex, allocatable, save, dimension(:,:), public :: fdisk
  complex, allocatable, save, dimension(:,:), public :: rhsk

  !> potential phi(nmod,nx,ns)
  complex, public, allocatable, target, save :: phi(:,:,:)
  
  !> The parallel component of the vector potential  apar(nmod,nx,ns)
  complex, public, allocatable, target, save :: apar(:,:,:)
  
  !> The parallel component of the magnetic field perturbation bpar(nmod,nx,ns)
  complex, public, allocatable, target, save :: bpar(:,:,:)
  
  !> The position where phi starts in the solution
  integer, public, save :: n_phi_start

  !> the size of a central difference stencil for the respective coordinate
  !> direction. The total size of the stencil is usually
  !> stencil_side(id_something)*2+1 .
  integer, public, save :: stencil_side(6)
  !> the size of a central difference stencil only for the zonal flow
  !> for all coordinate directions
  integer, public, save :: stencil_side_zf(6)

  real, public, save :: nmat_factor_for_deriv(6)
  
  !> The total number of points to be communicated in some direction.
  integer, public, save :: ghost_size_vpar, ghost_size_s, ghost_size_mu
  integer, public, save :: ghost_size_x, ghost_size_vpar_mu, ghost_size_vpar_s
  integer, public, save :: ghost_size_x_phi, ghost_size_x_pga, ghost_size_x_f
  integer, public, save :: ghost_size_x2_pga, ighost_sbp, ighost_sbn
  integer, public, save :: ighost_vparbp, ighost_vparbn
  integer, public, save :: ighost_xbp, ighost_xbn
  integer, public, save :: ighost_mubp, ighost_mubn
  integer, public, save :: ighost_vparbp_mubp, ighost_vparbn_mubp
  integer, public, save :: ighost_vparbp_mubn, ighost_vparbn_mubn
  integer, public, save :: ighost_vparbp_sbp, ighost_vparbn_sbp
  integer, public, save :: ighost_vparbp_sbn, ighost_vparbn_sbn

  !> The number of ghost points in various directions (at each end).
  integer, public, save  :: ghost_points_vpar, ghost_points_s, ghost_points_mu
  integer, public, save  :: ghost_points_vpar_mu, ghost_points_vpar_s
  integer, public, save  :: ghost_points_x
  
  !> the distribution only needs ghost_points_xf, less than ghost_points_x
  !> used for the fields  
  integer, public, save :: ghost_points_xf 

  !
  ! private variables
  !

  !> size of fields counter
  integer, save :: nfields
  !> to check in initialization routines have been called
  logical :: first_call = .true.
  !> various offsets for each field
  integer, save :: ioffset_phi, ioffset_apar, ioffset_bpar
  integer, save :: ioffset_phi_ga, ioffset_apar_ga, ioffset_bpar_ga  
  integer, save :: ioffset_phi2, ioffset_apar2, ioffset_bpar2
  integer, save :: ioffset_phi_ga2, ioffset_apar_ga2, ioffset_bpar_ga2
  integer, save :: n_phi, n_apar, n_bpar, n_phi_ga, n_apar_ga, n_bpar_ga 
  integer, save :: nelem_phi, nelem_apar, nelem_bpar
  integer, save :: nelem_phi_ga, nelem_apar_ga, nelem_bpar_ga 
  integer, save :: n_mom_conserve, n_ene_conserve  
  !> needed by mpighosts
  public :: ioffset_phi2, ioffset_phi_ga2
  
contains

!****************************************************************************
!> subroutine that intializes everything required to use dist
!> The distribution function is stored in in a single vector as
!> [FFFFFFFFFFFFFFF PPPPPP ffffff ppp ffffff ppp ffff pp ffff pp]
!> where F is g-dist, P is the fields, and the lower cases represent 
!> ghost point blocks, left and right, for a particular direction
!----------------------------------------------------------------------------
subroutine dist_init(mom_conservation,ene_conservation,blending_order)

  use general,        only : gkw_abort, gkw_warn, gkw_exit
  use global,         only : int2char, DISTRIBUTION
  use global,         only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
  use global,         only : id_vpar,id_mu,id_s,id_mod,id_x,id_sp
  use index_function, only : register_offset, index_init,                    &
                           & index_set_ghostpoints, indx
  use grid,           only : nsp, nx, ns, nmu, nvpar, nmod, n_x_grid, lx, mumax,   &
                           & parallel_s, lsendrecv_mu, lsendrecv_x, parallel_vpar, &
                           & n_mu_grid, psih
  use control, only : nlphi, nlapar, nlbpar, vp_trap
  use control, only : order_of_the_scheme, order_of_the_radial_scheme
  use control, only : order_of_the_zf_scheme, lcollisions, ltrapping_arakawa
  use control, only : spectral_radius, method, flux_tube, uniform_mu_grid
  use components,     only : max_rho, rhostar
  use geom,           only : eps, geom_type
  use general,        only : max_factor
  use mpiinterface,   only : root_processor

  logical, intent(in) :: mom_conservation
  logical, intent(in) :: ene_conservation
  integer, intent(in) :: blending_order

  real :: rdum, psim
  integer :: idum
  integer :: iquantity, direction
  logical :: derivs_in_lin_terms(BPAR_GA_FIELD,6,6)

  ! we have got here into dist, so at least an attempt was made to call these
  first_call = .false.

  ! no stencil is needed for spectral directions, or directions
  ! without finite differences (like species)
  stencil_side = 0
  stencil_side_zf = 0

  ! The stencil of terms with finite differences will involve more
  ! points, depending on the scheme:
  select case(order_of_the_scheme)
  case('second_order')
    stencil_side(id_s) = 1
  case('fourth_order')
    stencil_side(id_s) = 2
  case default
    call gkw_abort('order_of_the_scheme: use second_order or fourth_order')
  end select
  
  stencil_side(id_vpar) = stencil_side(id_s)
  ! We only have second order derivatives in mu so far.
  stencil_side(id_mu) = 1
  
  if (spectral_radius) then
    ! computation is spectral in x
    stencil_side(id_x) = 0
  else
    ! finite differences in x, just as in s or different, if stated in
    ! the input file.
    stencil_side(id_x) = stencil_side(id_s)
    select case(order_of_the_radial_scheme)
    case('second_order')
      stencil_side(id_x) = max(1,stencil_side(id_x))
    case('fourth_order')
      stencil_side(id_x) = max(2,stencil_side(id_x))
    case default
      call gkw_abort('order_of_the_radial_scheme: use at most ''fourth_order''')
    end select
  end if

  ! Allow to use higher order finite differences for the zonal
  ! flow.
  ! FIXME Currently, this leads to increased communication and memory
  ! consumption for *all* toroidal modes.
  select case(order_of_the_zf_scheme)
  case('second_order')
    stencil_side_zf(id_s) = 1
  case('fourth_order')
    stencil_side_zf(id_s) = 2
  case('sixth_order')
    stencil_side_zf(id_s) = 3
  case default
    call gkw_abort('order_of_the_zf_scheme: use at most ''sixth_order'' not '''&
       & //order_of_the_zf_scheme//'''')
  end select


  ! The size of the solution, without any fields
  nf = nsp*nx*ns*nmu*nvpar*nmod

  nmat_factor_for_deriv(id_s) = real((nmod-1)*(2*stencil_side(id_s)) &
     & + (2*stencil_side_zf(id_s)))/real(nmod)
  nmat_factor_for_deriv(id_x) = 2*stencil_side(id_x)
  nmat_factor_for_deriv(id_mu) = 2*stencil_side(id_mu)
  nmat_factor_for_deriv(id_vpar) = 2*stencil_side(id_vpar)
  nmat_factor_for_deriv(id_mod) = 2*stencil_side(id_mod)
  nmat_factor_for_deriv(id_sp) = 2*stencil_side(id_sp)


  ! Set the array size of the main matrix 
  ! These numbers are informed based on the terms and checked by tests
  ! But there could be cases where they are insufficient or oversufficient,
  ! in which case you can simply increase / decrease them.

  ! The various linear terms contribute to the elements of the matrix,
  ! but some overlap.

  ! the diagonal elements of the gyrokinetic equation
  ntot = nf

  ! the diagonal elements of the Poisson and Ampere equations
  if(nlphi) ntot = ntot + nf
  if(nlapar) ntot = ntot + nf
  if(nlbpar) ntot = ntot + nf

  ! the central elements of all field derivatives
  if(nlphi) ntot = ntot + nf
  if(nlapar) ntot = ntot + nf
  if(nlbpar) ntot = ntot + nf

  derivs_in_lin_terms = .false.
  
  ! term I : s-deriv of distribution
  derivs_in_lin_terms(DISTRIBUTION,id_s,id_s) = .true.

  ! term II: x and y (and s) deriv of distribution
  derivs_in_lin_terms(DISTRIBUTION,id_x,id_x) = .true.
  derivs_in_lin_terms(DISTRIBUTION,id_mod,id_mod) = .true.
  if(rhostar > 0.0) then
    derivs_in_lin_terms(DISTRIBUTION,id_s,id_s) = .true.
  end if

  ! term III is nonlinear...

  ! term IV: vpar-derivative of distribution
  derivs_in_lin_terms(DISTRIBUTION,id_vpar,id_vpar) = .true.

  ! term V: x and y (and s) deriv of phi_ga and apar_ga and bpar_ga
  derivs_in_lin_terms(PHI_GA_FIELD,id_x,id_x) = .true.
  derivs_in_lin_terms(PHI_GA_FIELD,id_mod,id_mod) = .true.
  if(rhostar > 0.0) then
    derivs_in_lin_terms(PHI_GA_FIELD,id_s,id_s) = .true.
  end if
  if (nlapar) then
    ! x and y and s deriv of apar_ga
    derivs_in_lin_terms(APAR_GA_FIELD,id_x,id_x) = .true.
    derivs_in_lin_terms(APAR_GA_FIELD,id_mod,id_mod) = .true.
    if(rhostar > 0.0) then
      derivs_in_lin_terms(APAR_GA_FIELD,id_s,id_s) = .true.
    end if

  end if
  if (nlbpar) then
    ! x and y and s deriv of apar_ga
    derivs_in_lin_terms(BPAR_GA_FIELD,id_x,id_x) = .true.
    derivs_in_lin_terms(BPAR_GA_FIELD,id_mod,id_mod) = .true.
    if(rhostar > 0.0) then
      derivs_in_lin_terms(BPAR_GA_FIELD,id_s,id_s) = .true.
    end if
  end if

  ! term VI is a pure background term... this is like a source!

  ! term VII: s-deriv of phi_ga
  derivs_in_lin_terms(PHI_GA_FIELD,id_s,id_s) = .true.

  ! term VIII: x and y and s deriv of phi_ga
  derivs_in_lin_terms(PHI_GA_FIELD,id_x,id_x) = .true.
  derivs_in_lin_terms(PHI_GA_FIELD,id_mod,id_mod) = .true.
  if(rhostar > 0.0) then
    derivs_in_lin_terms(PHI_GA_FIELD,id_s,id_s) = .true.
  end if

  ! there is no term IX in the final set of equations...

  ! term X: s deriv of bpar_ga
  derivs_in_lin_terms(BPAR_GA_FIELD,id_s,id_s) = .true.

  ! term XI: x and y and s deriv of bpar_ga
  derivs_in_lin_terms(BPAR_GA_FIELD,id_x,id_x) = .true.
  derivs_in_lin_terms(BPAR_GA_FIELD,id_mod,id_mod) = .true.
  derivs_in_lin_terms(BPAR_GA_FIELD,id_s,id_s) = .true.

  if (ltrapping_arakawa) ntot= ntot + 6*nf !guess
  
  if (lcollisions) then
    derivs_in_lin_terms(DISTRIBUTION,id_vpar,id_vpar) = .true.
    derivs_in_lin_terms(DISTRIBUTION,id_mu,id_mu) = .true.
    ! note that there are mixed derivatives in the collision operator:
    derivs_in_lin_terms(DISTRIBUTION,id_vpar,id_mu) = .true.
  end if

  ! count elements produced by mixed derivatives (this is supposedly
  ! only in the collisionop, but for the sake of generality...)
  do iquantity = 1, size(derivs_in_lin_terms,1)
    do direction = 1, size(derivs_in_lin_terms,2)
      do idum = direction, size(derivs_in_lin_terms,2)
        if(derivs_in_lin_terms(iquantity, direction, idum) .or. &
           & derivs_in_lin_terms(iquantity, idum, direction)) then
          if(direction == idum) then
            ! count the elements produced by unmixed first or higher derivatives
            ntot = ntot + ceiling(nmat_factor_for_deriv(direction)*nf)
          else
            ! the central element an outer derivative overlaps with an
            ! 'unmixed' derivative
            ntot = ntot + ceiling( &
               & nmat_factor_for_deriv(direction)*nmat_factor_for_deriv(idum)*nf)
            if(.not.derivs_in_lin_terms(iquantity, &
               & min(direction,idum), min(direction,idum))) then
              call gkw_abort('This should not be the case.')
            end if
          end if
        end if
      end do
    end do
  end do

  nelem_cc = 0
  if (mom_conservation) then
      !Minimum that will always compress in unlimited time:
      !Exact stencil + compression workspace:
      !nelem_cc = nsp*nx*ns*nmu*nvpar*nmod + 25*(nmu+2)*(nvpar+2)
      !But can be a little more generous to speed compression
      nelem_cc= nelem_cc +int(1.3*nf) +50*(nmu+2)*(nvpar+2)
  end if
  if (ene_conservation) then  !same as mom_conservation
      nelem_cc= nelem_cc +int(1.3*nf) +50*(nmu+2)*(nvpar+2)
  end if

  ! move this calculation to grid so it can be used for checking parallel layout
  if (.not. spectral_radius) then
      ! precalculate maximum rho       
      if (flux_tube) then
        psim=eps
      else
        psim=psih  
      endif
      rdum = 1.0*max_rho*sqrt(mumax*(1.0+psim))

      select case(geom_type)
        ! approximate g11 factor for non circular geometries
        ! might need to be adjusted if gyro_average aborts
        case('chease','miller')
          rdum = rdum*(1.1+psim)
        case('slab','slab_periodic');  
          rdum = 1.0*max_rho*sqrt(mumax)
        case default
      end select

      if (uniform_mu_grid) then
        rdum =rdum*sqrt((real(n_mu_grid)-0.5)/real(n_mu_grid))
      else
        rdum =rdum*(real(n_mu_grid)-0.5)/real(n_mu_grid)
      end if

      ! need an extra point for higher blending orders
      idum = 0
      if (blending_order > 2) idum = 1
      
      n_gav_bound = ceiling(rdum*real(n_x_grid) / lx) + n_gav_bound_ex + idum
  else
      n_gav_bound = 0
  end if

  !The number of elements for the neoclassics diagnostic - only the (0,0) mode
  !is needed
  nelem_nc = 72*nsp*ns*nmu*(nvpar+1)

  nfields = 0

  !
  ! (1) Work out how many grid points to communicate between adjacent
  !     processes. This depends on the order of the scheme.
  !
  
  ghost_points_vpar = 0
  ghost_points_s = 0
  ghost_points_mu = 0
  ghost_points_x = 0
  ghost_points_xf = 0
  ghost_points_vpar_mu = 0
  ghost_points_vpar_s  = 0
  ! using higher order finite differences for the zonal flow
  ! requires more ghost points; currently no distinction between the
  ! modes is made.
  if (parallel_s)    ghost_points_s = max(stencil_side(id_s),stencil_side_zf(id_s))
  if (lsendrecv_x)   ghost_points_xf = stencil_side(id_x)
  if (lsendrecv_mu)  ghost_points_mu = stencil_side(id_mu)
  if (parallel_vpar.and.vp_trap==0) then
    ghost_points_vpar = stencil_side(id_vpar)
    ! only 1 point needed below
    if (lsendrecv_mu) ghost_points_vpar_mu = 1
    ! only 1 point needed below
    if (parallel_s .and. ltrapping_arakawa) then
      ghost_points_vpar_s = 1
    end if
  end if

  ! more ghost points are needed for gyroaverage
  ! move this calculation to grid so it can be used for checking parallel layout
  if (lsendrecv_x) then
    ghost_points_x = max(n_gav_bound,ghost_points_xf)
    
    if (root_processor) then
      write(*,*)
      write(*,*) '*** ghost_points_x:    ', ghost_points_x 
      write(*,*) '*** maximum n_procs_x: ',  &
         &  max_factor(n_x_grid,n_x_grid / ghost_points_x)
      write(*,*)
    end if

    if (ghost_points_x > nx) then
      call gkw_warn('The maximum gyroradius cannot cross more than two x processors')    
      call gkw_warn('Reduce n_procs_x to ' // &
         & trim(int2char(max_factor(n_x_grid,n_x_grid / ghost_points_x),0)))
      call gkw_exit('More x ghosts points than x points per proc')
    end if
  end if
  
  !
  ! (2) Work out the total number of elements of complex datatype to be sent
  !     to the adjacent processors based on the required derivatives.
  !     Additional contributions will come from the fields.
  !
  
  ghost_size_vpar    = ghost_points_vpar*nsp*nx*ns*nmu*nmod
  ghost_size_mu      = ghost_points_mu*nsp*nx*ns*nvpar*nmod
  ! FJC: the index function cannot yet handle different numbers
  ! of ghost points for fields and f, so use larger for now
  ! the distribution only needs ghost_points_xf 
  !ghost_size_x       = ghost_points_xf*nmu*nsp*ns*nmod*nvpar
  ghost_size_x       = ghost_points_x*nmu*nsp*ns*nmod*nvpar
  ghost_size_s       = ghost_points_s*nsp*nx*nvpar*nmu*nmod
  ghost_size_vpar_mu = ghost_points_vpar_mu*nsp*nx*ns*nmod
  ghost_size_vpar_s  = ghost_points_vpar_s*nmu*nsp*nx*nmod
  ghost_size_x_phi   = 0 !< fields only
  ghost_size_x_pga   = 0 !< ga fields only
  ghost_size_x2_pga  = 0 !< ga fields (2 gp) only
  
  !> presently only used for the datatype, not the index function
  ghost_size_x_f     = ghost_points_xf*nmu*nsp*ns*nmod*nvpar
    
  !
  ! (3) Add the fields and their contributions to the ghost elements.
  !

  number_of_fields = 0
  n_phi = 0; n_phi_ga = 0; n_apar = 0; n_apar_ga = 0; n_bpar = 0; n_bpar_ga = 0 
  n_mom_conserve=0; n_ene_conserve=0 

  ! if the electro-static potential is kept increase the size
  if (nlphi) then
    number_of_fields = number_of_fields + 1
    n_phi       = nf + nfields         ! the phi offset (=nf)
    nelem_phi   = nx*ns*nmod           ! number of elements  
    nfields     = nfields + nelem_phi  ! total size of fields
    ! the offsets of phi within the ghost block
    ioffset_phi  = ghost_size_s    
    ioffset_phi2 = ghost_size_x
    ! increase the ghost size; derivatives in phi are of same order as fdisi
    ghost_size_s     = ghost_size_s     + nx*nmod*ghost_points_s
    ghost_size_x     = ghost_size_x     + nmod*ns*ghost_points_x
    ghost_size_x_phi = ghost_size_x_phi + nmod*ns*ghost_points_x
  endif 

  ! if the parallel vector potential is kept increase the size
  if (nlapar) then
    number_of_fields = number_of_fields + 1
    n_apar      = nf + nfields
    nelem_apar  = nx*ns*nmod
    nfields     = nfields + nelem_apar
    ioffset_apar  = ghost_size_s 
    ioffset_apar2 = ghost_size_x
    ghost_size_s     = ghost_size_s     + nx*nmod*ghost_points_s
    ghost_size_x     = ghost_size_x     + nmod*ns*ghost_points_x
    ghost_size_x_phi = ghost_size_x_phi + nmod*ns*ghost_points_x
  endif 

  ! if the parallel magnetic field is kept increase the size
  if (nlbpar) then
    number_of_fields = number_of_fields + 1
    n_bpar      = nf + nfields
    nelem_bpar  = nx*ns*nmod
    nfields     = nfields + nelem_bpar
    ioffset_bpar  = ghost_size_s 
    ioffset_bpar2 = ghost_size_x
    ghost_size_s     = ghost_size_s     + nx*nmod*ghost_points_s
    ghost_size_x     = ghost_size_x     + nmod*ns*ghost_points_x
    ghost_size_x_phi = ghost_size_x_phi + nmod*ns*ghost_points_x
  endif 

  ! nregular_fields_end is the total number of elements of the distribution and
  ! the fields, EXCLUDING gyroavg fields and EXCLUDING the collisionop-related fields
  nregular_fields_end = nf + nfields
  ! N.B. no more actual fields are allowed after here.

  if ((nlphi).and.(.not. spectral_radius)) then 
    n_phi_ga         = nf + nfields 
    nelem_phi_ga     = nmod*nx*ns*nmu*nsp 
    nfields          = nfields + nelem_phi_ga 
    ! the offset of phi_ga within the ghost block 
    ioffset_phi_ga   = ghost_size_s 
    ioffset_phi_ga2  = ghost_size_x
    ! increase the ghost size: derivatives in phi_ga 
    ghost_size_s     = ghost_size_s     + nmod*nx*nmu*nsp*ghost_points_s
    ghost_size_x     = ghost_size_x     + nmod*ns*nmu*nsp*ghost_points_x
    ghost_size_x_pga = ghost_size_x_pga + nmod*ns*nmu*nsp*ghost_points_x
    ghost_size_x2_pga = ghost_size_x2_pga + nmod*ns*nmu*nsp*ghost_points_xf
  endif 

  if ((nlapar).and.(.not. spectral_radius)) then 
    n_apar_ga         = nf + nfields 
    nelem_apar_ga     = nmod*nx*ns*nmu*nsp 
    nfields           = nfields + nelem_apar_ga 
    ! the offset of apar_ga within the ghost block 
    ioffset_apar_ga   = ghost_size_s 
    ioffset_apar_ga2  = ghost_size_x
    ! increase the ghost size: derivatives in phi_ga 
    ghost_size_s     = ghost_size_s     + nmod*nx*nmu*nsp*ghost_points_s
    ghost_size_x     = ghost_size_x     + nmod*ns*nmu*nsp*ghost_points_x
    ghost_size_x_pga = ghost_size_x_pga + nmod*ns*nmu*nsp*ghost_points_x
    ghost_size_x2_pga = ghost_size_x2_pga + nmod*ns*nmu*nsp*ghost_points_xf
  endif 

  if ((nlbpar).and.(.not. spectral_radius)) then 
    n_bpar_ga        = nf + nfields 
    nelem_bpar_ga    = nmod*nx*ns*nmu*nsp 
    nfields          = nfields + nelem_bpar_ga 
    ! the offset of bpar_ga within the ghost block 
    ioffset_bpar_ga   = ghost_size_s 
    ioffset_bpar_ga2  = ghost_size_x
    ! increase the ghost size: derivatives in phi_ga 
    ghost_size_s     = ghost_size_s     + nmod*nx*nmu*nsp*ghost_points_s
    ghost_size_x     = ghost_size_x     + nmod*ns*nmu*nsp*ghost_points_x
    ghost_size_x_pga = ghost_size_x_pga + nmod*ns*nmu*nsp*ghost_points_x
    ghost_size_x2_pga = ghost_size_x2_pga + nmod*ns*nmu*nsp*ghost_points_xf
  endif 
  
  nelem_conserve = 0
  n_conserve = nf + nfields
  
  ! The momentum conserving `field': no derivatives.
  if (mom_conservation) then  
    n_mom_conserve = nf + nfields
    nelem_conserve = nelem_conserve + nx*nmod*ns*nsp
    nfields = nfields + nx*nmod*ns*nsp     
  end if
  
  ! The energy conserving `field': no derivatives.
  if (ene_conservation) then  
    n_ene_conserve = nf + nfields
    nelem_conserve = nelem_conserve + nx*nmod*ns*nsp
    nfields = nfields + nx*nmod*ns*nsp
  end if   
 
  ! nsolc is the total number of elements of the distribution and the
  ! fields, INCLUDING gyroavg fields and INCLUDING the collisionop-related fields
  nsolc = nf + nfields

  ! matrix for f-to-g conversion
  if (nlapar) then
    nelem_g2f = nsolc+nf
    if (method=='IMP') ntot = ntot + nelem_g2f
  else
    nelem_g2f = 0
  end if


  !
  ! (4) Decide where to store the points received from adjacent processors.
  !
  
  ! Start at nsolc. These variables hold the index of the start of the
  ! respective ghost cell block. (Actually the first ghost cell is
  ! then at i_ghost_*[pn] + 1)
  ighost_vparbp = nsolc                            ! previous proc in vpar
  ighost_vparbn = ighost_vparbp + ghost_size_vpar  ! next proc in vpar
  ighost_mubp   = ighost_vparbn + ghost_size_vpar  ! previous proc in mu
  ighost_mubn   = ighost_mubp   + ghost_size_mu    ! next pro in mu
  ighost_vparbp_mubp = ighost_mubn + ghost_size_mu ! prev mu, prev vpar
  ighost_vparbn_mubp = ighost_vparbp_mubp + ghost_size_vpar_mu
  ighost_vparbp_mubn = ighost_vparbn_mubp + ghost_size_vpar_mu 
  ighost_vparbn_mubn = ighost_vparbp_mubn + ghost_size_vpar_mu
  ighost_vparbp_sbp  = ighost_vparbn_mubn + ghost_size_vpar_mu
  ighost_vparbn_sbp  = ighost_vparbp_sbp  + ghost_size_vpar_s
  ighost_vparbp_sbn  = ighost_vparbn_sbp  + ghost_size_vpar_s
  ighost_vparbn_sbn  = ighost_vparbp_sbn  + ghost_size_vpar_s
  ighost_sbp         = ighost_vparbn_sbn  + ghost_size_vpar_s
  ighost_sbn         = ighost_sbp         + ghost_size_s
  ighost_xbp         = ighost_sbn         + ghost_size_s
  ighost_xbn         = ighost_xbp         + ghost_size_x
    
  msolc = ighost_xbn + ghost_size_x
  ! Initialise the index function. msolc is the maximum value
  call index_init(msolc)

  !
  ! (5) Register all the offsets with the index function. This allows the
  !     index function to relate call patterns with offsets.
  !
  
  ! The main part of fdisi
  call register_offset(imod=0,ix=0,is=0,imu=0,ivpar=0,isp=0,ioffset=0,ifield=ifdis)
  
  ! The fields, together with the tags that are used to reference them.
  ! The momentum conserving `field' also has a number of local species
  ! associated with it.
  call register_offset(ioffset=n_phi,    ifield=iphi    )
  call register_offset(ioffset=n_apar,   ifield=iapar   ) 
  call register_offset(ioffset=n_bpar,   ifield=ibpar   ) 
  call register_offset(ioffset=n_phi_ga, ifield=iphi_ga )
  call register_offset(ioffset=n_apar_ga,ifield=iapar_ga) 
  call register_offset(ioffset=n_bpar_ga,ifield=ibpar_ga) 
  call register_offset(ioffset=n_mom_conserve,ifield=i_mom,nsp=nsp)
  call register_offset(ioffset=n_ene_conserve,ifield=i_ene,nsp=nsp)

  ! The ghost cells for fdisi
  call register_offset(ivpar=-1,ioffset=ighost_vparbp,ifield=ifdis) 
  call register_offset(ivpar=+1,ioffset=ighost_vparbn,ifield=ifdis)
  call register_offset(imu=-1,ioffset=ighost_mubp,ifield=ifdis)
  call register_offset(imu=+1,ioffset=ighost_mubn,ifield=ifdis)
  call register_offset(imu=-1,ivpar=-1,ioffset=ighost_vparbp_mubp,ifield=ifdis)
  call register_offset(imu=-1,ivpar=+1,ioffset=ighost_vparbn_mubp,ifield=ifdis)
  call register_offset(imu=+1,ivpar=-1,ioffset=ighost_vparbp_mubn,ifield=ifdis)
  call register_offset(imu=+1,ivpar=+1,ioffset=ighost_vparbn_mubn,ifield=ifdis)
  call register_offset(is=-1,ivpar=-1,ioffset=ighost_vparbp_sbp,ifield=ifdis)
  call register_offset(is=-1,ivpar=+1,ioffset=ighost_vparbn_sbp,ifield=ifdis)
  call register_offset(is=+1,ivpar=-1,ioffset=ighost_vparbp_sbn,ifield=ifdis)
  call register_offset(is=+1,ivpar=+1,ioffset=ighost_vparbn_sbn,ifield=ifdis)
  call register_offset(is=-1,ioffset=ighost_sbp,ifield=ifdis)
  call register_offset(is=+1,ioffset=ighost_sbn,ifield=ifdis)
  call register_offset(ix=-1,ioffset=ighost_xbp,ifield=ifdis)
  call register_offset(ix=+1,ioffset=ighost_xbn,ifield=ifdis)
  
  ! the ghost cells for the fields
  call register_offset(is=-1,ioffset=ighost_sbp+ioffset_phi,ifield=iphi)
  call register_offset(is=-1,ioffset=ighost_sbp+ioffset_apar,ifield=iapar)
  call register_offset(is=-1,ioffset=ighost_sbp+ioffset_bpar,ifield=ibpar)

  call register_offset(is=-1,ioffset=ighost_sbp+ioffset_phi_ga,ifield=iphi_ga)
  call register_offset(is=-1,ioffset=ighost_sbp+ioffset_apar_ga,ifield=iapar_ga)
  call register_offset(is=-1,ioffset=ighost_sbp+ioffset_bpar_ga,ifield=ibpar_ga)

  call register_offset(is=+1,ioffset=ighost_sbn+ioffset_phi,ifield=iphi)
  call register_offset(is=+1,ioffset=ighost_sbn+ioffset_apar,ifield=iapar)
  call register_offset(is=+1,ioffset=ighost_sbn+ioffset_bpar,ifield=ibpar)

  call register_offset(is=+1,ioffset=ighost_sbn+ioffset_phi_ga,ifield=iphi_ga)
  call register_offset(is=+1,ioffset=ighost_sbn+ioffset_apar_ga,ifield=iapar_ga)
  call register_offset(is=+1,ioffset=ighost_sbn+ioffset_bpar_ga,ifield=ibpar_ga)
  
  call register_offset(ix=-1,ioffset=ighost_xbp+ioffset_phi2,ifield=iphi)
  call register_offset(ix=-1,ioffset=ighost_xbp+ioffset_apar2,ifield=iapar)
  call register_offset(ix=-1,ioffset=ighost_xbp+ioffset_bpar2,ifield=ibpar)

  call register_offset(ix=-1,ioffset=ighost_xbp+ioffset_phi_ga2,ifield=iphi_ga)
  call register_offset(ix=-1,ioffset=ighost_xbp+ioffset_apar_ga2,ifield=iapar_ga)
  call register_offset(ix=-1,ioffset=ighost_xbp+ioffset_bpar_ga2,ifield=ibpar_ga)

  call register_offset(ix=+1,ioffset=ighost_xbn+ioffset_phi2,ifield=iphi)
  call register_offset(ix=+1,ioffset=ighost_xbn+ioffset_apar2,ifield=iapar)
  call register_offset(ix=+1,ioffset=ighost_xbn+ioffset_bpar2,ifield=ibpar)

  call register_offset(ix=+1,ioffset=ighost_xbn+ioffset_phi_ga2,ifield=iphi_ga)
  call register_offset(ix=+1,ioffset=ighost_xbn+ioffset_apar_ga2,ifield=iapar_ga)
  call register_offset(ix=+1,ioffset=ighost_xbn+ioffset_bpar_ga2,ifield=ibpar_ga)
  
  !
  ! (6) set up ghost points in the index function
  !
  ! FJC: the index function cannot yet handle different numbers
  ! of ghost points for feilds and f, so use larger for now
  call index_set_ghostpoints(gp_x=ghost_points_x,gp_s=ghost_points_s,gp_mu=ghost_points_mu,&
      & gp_vpar=ghost_points_vpar,gp_vpar_mu=ghost_points_vpar_mu,     &
      & gp_vpar_s=ghost_points_vpar_s)
  
  !
  ! (7) For each set of points to be sent to adjacent processor ghost zones,
  !     create a datatype to aid communication ( i.e. in sending part of fdisi
  !     to the offset point on the receiving processor).
  !
  call setup_types_for_ghost_send
 
! Starting point of the fields
  n_phi_start = n_phi+1

  ! we now allocate everything else required in this module
  call dist_allocate()

end subroutine dist_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Create arrays for copying parts of the solution into communication buffers
!>
!> For example, the first datatype is for communicating the field and fdisi to
!> next processor in s. The call to get_block_bounds returns an array 
!> of starts and ends. The values in the start array
!> (not necessarily in this order) will be (1,1,1,1,1,ns-2+1).
!> The ends in this case would then be (nvpar,nx,nmod,nmu,nsp,ns)
!> If you called the index function indx(), looping over those
!> starts(i),ends(i), you would get the last 2 s-points in fdisi.
!> 
!> Next we call fill_type_index, where we loop over those values.
!> Because the index function can change its order, we use match_indices_with_dims() 
!> to map n,m,l,k,j,i to the correct indices. That allows us to call
!> the index function and recieve the locations of the 2 last s points
!> of fdisi in an order which will correspond to contiguous locations in the
!> ghost points on the receiving processor.
!> The returned location is stored in index_array in sequence, which is only
!> reset when fill_type_index is called for the main part (before doing the
!> next communicator, index_array will be deallocated, reallocated and zeroed)
!> 
!> For the s direction, we then also use get_block_bounds to get the locations
!> of the local fields that we want to send. These points are put into
!> the index array as well.  Note that all these indices are all relative 
!> to the _start_ of fdisi, and (obviously) may not be contiguous.
!> 
!> Then we create the datatype, TYPE_NEXT_S. That corresponds to an array of 
!> pointers to the complex values that we want to send, determined by
!> the index_array. Think of it as an implicit `packing'. This datatype has
!> size size_of_complex*ghost_points_s, i.e. ghost_points_s complex values.
!> 
!> When we use the datatype to mpi send, we are just sending 1 of type
!> TYPE_NEXT_S. That is picking up the right points, reordering and sending 
!> them in the right sized buffer. Because the index_array is created with 
!> this correct reordering certain order to start with, 
!> when the data is received at the next processor, no `unpack'.is needed 
!> because the ghost_point_s complex values are already in the right order. 
! 
!-----------------------------------------------------------------------------
subroutine setup_types_for_ghost_send
!-----------------------------------------------------------------------------

  use general,        only : gkw_abort
  use control,        only : nlapar, nlbpar, spectral_radius 
  use index_function, only : get_block_bounds, get_block_bounds_x_hack
  use index_function, only : IS_6D_DISTR, IS_3D_FIELD, IS_GYROAVG_FIELD
  use mpidatatypes,   only : TYPE_NEXT_MU, TYPE_NEXT_S, TYPE_NEXT_S_NEXT_VPAR
  use mpidatatypes,   only : TYPE_NEXT_S_PREV_VPAR, TYPE_NEXT_VPAR
  use mpidatatypes,   only : TYPE_NEXT_VPAR_NEXT_MU, TYPE_NEXT_VPAR_PREV_MU
  use mpidatatypes,   only : TYPE_PREV_MU, TYPE_PREV_S, TYPE_PREV_S_NEXT_VPAR
  use mpidatatypes,   only : TYPE_PREV_S_PREV_VPAR, TYPE_PREV_VPAR
  use mpidatatypes,   only : TYPE_PREV_VPAR_NEXT_MU, TYPE_PREV_VPAR_PREV_MU
  use mpidatatypes,   only : TYPE_NEXT_X2_F, TYPE_NEXT_X2_F_RECV
  use mpidatatypes,   only : TYPE_NEXT_X2_PGA, TYPE_NEXT_X2_PGA_RECV 
  use mpidatatypes,   only : TYPE_NEXT_X_PGA, TYPE_PREV_X2_F, TYPE_PREV_X_PGA
  use mpidatatypes,   only : TYPE_PREV_X2_F_RECV, TYPE_PREV_X2_PGA
  use mpidatatypes,   only : TYPE_PREV_X2_PGA_RECV, TYPE_NEXT_X_PHI
  use mpidatatypes,   only : TYPE_PREV_X_PHI

  integer, allocatable :: index_array(:)
  integer :: ii, ierr

  integer, dimension(6) :: starts,ends

  !
  ! S, including fields
  !
  
  s : if (ghost_size_s > 0) then
    
    allocate(index_array(ghost_size_s),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for s')
    
    ! main part for next s
    call get_block_bounds(starts,ends,IS_6D_DISTR,gps_next=ghost_points_s)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
  
    ! fields for next s
    call get_block_bounds(starts,ends,IS_3D_FIELD,gps_next=ghost_points_s)
    call fill_type_index(iphi, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar, starts, ends, index_array)
    if (.not.spectral_radius) then 
      call get_block_bounds(starts,ends,IS_GYROAVG_FIELD,gps_next=ghost_points_s)
      call fill_type_index(iphi_ga, starts, ends, index_array) 
      if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array)
      if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array)  
    endif 
    
    ! create the communicator    
    call get_ghost_block_type(index_array,TYPE_NEXT_S)

    ! main part for prev s
    call get_block_bounds(starts,ends,IS_6D_DISTR,gps_prev=ghost_points_s)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)

    ! fields for prev s
    call get_block_bounds(starts,ends,IS_3D_FIELD,gps_prev=ghost_points_s)
    call fill_type_index(iphi, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar, starts, ends, index_array)
    if (.not.spectral_radius) then
      call get_block_bounds(starts,ends,IS_GYROAVG_FIELD,gps_prev=ghost_points_s) 
      call fill_type_index(iphi_ga, starts, ends, index_array)
      if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array) 
      if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array) 
    endif     

    call get_ghost_block_type(index_array,TYPE_PREV_S)
    if (allocated(index_array)) deallocate(index_array)

  end if s
    
  ! X, distribution function (2 gp) only - reduced in size to ghost_points_xf
  ! until the index function can handle different ghost sizes for fields
  ! we create a recv datatype too as a temporary measure
  
  xf : if (ghost_size_x_f > 0) then
    
    allocate(index_array(ghost_size_x_f),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for x')
    
    ! main part for next x
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpx_next=ghost_points_xf)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)          
    call get_ghost_block_type(index_array,TYPE_NEXT_X2_F)

    ! main part for prev x
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpx_prev=ghost_points_xf)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_X2_F)
    
    ! main part for prev x
    call get_block_bounds(starts,ends,IS_6D_DISTR)
    call get_block_bounds_x_hack(starts,ends,ghost_points_xf)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_NEXT_X2_F_RECV)
    
    ! main part for prev x
    call get_block_bounds(starts,ends,IS_6D_DISTR)
    call get_block_bounds_x_hack(starts,ends,-ghost_points_xf)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_X2_F_RECV)

    if (allocated(index_array)) deallocate(index_array)

  end if xf
   
  
  ! X - normal fields only, all ghost points
  
  x2 : if (ghost_size_x_phi > 0) then
    
    allocate(index_array(ghost_size_x_phi),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for x')
    
    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for next x

    call get_block_bounds(starts,ends,IS_3D_FIELD,gpx_next=ghost_points_x)
    call fill_type_index(iphi, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_NEXT_X_PHI)

    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for prev x
    call get_block_bounds(starts,ends,IS_3D_FIELD,gpx_prev=ghost_points_x)
    call fill_type_index(iphi, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_PREV_X_PHI)

    if (allocated(index_array)) deallocate(index_array)

  end if x2
  
  ! X - gyro-averaged fields only, all ghost points
  
  x3 : if (ghost_size_x_pga > 0) then
    
    allocate(index_array(ghost_size_x_pga),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for x')
    
    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for next x
    call get_block_bounds(starts,ends,IS_GYROAVG_FIELD,gpx_next=ghost_points_x)
    call fill_type_index(iphi_ga, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_NEXT_X_PGA)

    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for prev x
    call get_block_bounds(starts,ends,IS_GYROAVG_FIELD,gpx_prev=ghost_points_x)
    call fill_type_index(iphi_ga, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_PREV_X_PGA)

    if (allocated(index_array)) deallocate(index_array)

  end if x3
  
  ! X - gyro-averaged fields (2gp) only.  
  ! Here we create also RECV types for noncontiguous received blocks
  
  x4 : if (ghost_size_x2_pga > 0) then
    
    allocate(index_array(ghost_size_x2_pga),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for x')
    
    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for next x
    call get_block_bounds(starts,ends,IS_GYROAVG_FIELD,gpx_next=ghost_points_xf)
    call fill_type_index(iphi_ga, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_NEXT_X2_PGA)

    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for prev x
    call get_block_bounds(starts,ends,IS_GYROAVG_FIELD,gpx_prev=ghost_points_xf)
    call fill_type_index(iphi_ga, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_PREV_X2_PGA)

    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for next x
    call get_block_bounds(starts,ends,IS_GYROAVG_FIELD)
    call get_block_bounds_x_hack(starts,ends,ghost_points_xf)
    call fill_type_index(iphi_ga, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_NEXT_X2_PGA_RECV)

    ! manual reset index_array counter for this case (since no fdis part)
    ii = 0
    
    ! fields for prev x
    call get_block_bounds(starts,ends,IS_GYROAVG_FIELD)
    call get_block_bounds_x_hack(starts,ends,-ghost_points_xf)
    call fill_type_index(iphi_ga, starts, ends, index_array)
    if (nlapar) call fill_type_index(iapar_ga, starts, ends, index_array) ! same size and call sequence as phi
    if (nlbpar) call fill_type_index(ibpar_ga, starts, ends, index_array)

    call get_ghost_block_type(index_array,TYPE_PREV_X2_PGA_RECV)

    if (allocated(index_array)) deallocate(index_array)

  end if x4
  

  !
  ! VPAR
  !
  
  vpar : if (ghost_size_vpar > 0) then
    
    allocate(index_array(ghost_size_vpar),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for vpar')
  
    ! next vpar
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpvpar_next=ghost_points_vpar)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_NEXT_VPAR)
    
    ! prev vpar
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpvpar_prev=ghost_points_vpar)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_VPAR)
    
    if (allocated(index_array)) deallocate(index_array)

  end if vpar 

  !
  ! MU
  !
  
  mu : if (ghost_size_mu > 0) then
    
    allocate(index_array(ghost_size_mu),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for mu')
  
    ! next mu
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpmu_next=ghost_points_mu)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_NEXT_MU)
    
    ! prev mu
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpmu_prev=ghost_points_mu)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_MU)
    
    if (allocated(index_array)) deallocate(index_array)

  end if mu    
  
  !
  ! VPAR-S N.B. there are 4, FOUR(!) DATATYPES HERE!
  !
  
  vpar_s : if (ghost_size_vpar_s > 0) then
    
    allocate(index_array(ghost_size_vpar_s),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for vpar-s')
    
    ! s next, vpar next
    call get_block_bounds(starts,ends,IS_6D_DISTR,gps_next=ghost_points_vpar_s,gpvpar_next=ghost_points_vpar_s)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_NEXT_S_NEXT_VPAR)
    
    ! s prev, vpar prev
    call get_block_bounds(starts,ends,IS_6D_DISTR,gps_prev=ghost_points_vpar_s,gpvpar_prev=ghost_points_vpar_s)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_S_PREV_VPAR)
    
    ! s next, vpar prev
    call get_block_bounds(starts,ends,IS_6D_DISTR,gps_next=ghost_points_vpar_s,gpvpar_prev=ghost_points_vpar_s)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array) 
   call get_ghost_block_type(index_array,TYPE_NEXT_S_PREV_VPAR)
      
    ! s prev, vpar next
    call get_block_bounds(starts,ends,IS_6D_DISTR,gps_prev=ghost_points_vpar_s,gpvpar_next=ghost_points_vpar_s)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_S_NEXT_VPAR)
    
    if (allocated(index_array)) deallocate(index_array)
    
  end if vpar_s
  
  !
  ! VPAR-MU N.B. there are 4 datatypes here!
  !
  
  vpar_mu : if (ghost_size_vpar_mu > 0) then
    
    allocate(index_array(ghost_size_vpar_mu),stat=ierr)
    if (ierr /= 0) call gkw_abort('dist_allocate: Could not allocate index_array for vpar-mu')
    
    ! mu next, vpar next
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpmu_next=ghost_points_vpar_mu,gpvpar_next=ghost_points_vpar_mu)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_NEXT_VPAR_NEXT_MU)

    ! mu prev, vpar prev
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpmu_prev=ghost_points_vpar_mu,gpvpar_prev=ghost_points_vpar_mu)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_VPAR_PREV_MU)
    
    ! mu next, vpar prev
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpmu_next=ghost_points_vpar_mu,gpvpar_prev=ghost_points_vpar_mu)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_PREV_VPAR_NEXT_MU)
    
    ! mu prev, vpar next
    call get_block_bounds(starts,ends,IS_6D_DISTR,gpmu_prev=ghost_points_vpar_mu,gpvpar_next=ghost_points_vpar_mu)
    ii = 0
    call fill_type_index(ifdis, starts, ends, index_array)
    call get_ghost_block_type(index_array,TYPE_NEXT_VPAR_PREV_MU)
  
    if (allocated(index_array)) deallocate(index_array)
    
  end if vpar_mu

  contains

  !*************************************************************************
    subroutine fill_type_index(field_id, starts, ends, index_array)
      use index_function, only : match_indices_with_dims, indx
    integer, intent(in) :: field_id
    integer, dimension(6), intent(in) :: starts,ends
    integer, dimension(:), intent(inout) :: index_array
    integer :: n, m, l, k, j, i
    integer :: i_s, i_sp, i_mu, i_vpar, i_x, i_mod
  
    do n=starts(6),ends(6)
      do m=starts(5),ends(5)
        do l=starts(4),ends(4)
          do k=starts(3),ends(3)
            do j=starts(2),ends(2)
              do i=starts(1),ends(1)
                call match_indices_with_dims(i,j,k,l,m,n,i_mod,i_x,i_s,i_mu,i_vpar,i_sp)
                ii=ii+1
                select case(field_id) 
                case(iphi,iapar,ibpar)
                  index_array(ii) = indx(field_id,i_mod,i_x,i_s)
                case(iphi_ga,iapar_ga,ibpar_ga) 
                  index_array(ii) = indx(field_id,i_mod,i_x,i_s,i_mu,i_sp)
                case(ifdis) 
                  index_array(ii) = indx(field_id,i_mod,i_x,i_s,i_mu,i_vpar,i_sp)
                case default 
                  call gkw_abort('Unkown field in fill_type_index')
                end select  
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine fill_type_index
  !*************************************************************************
  
end subroutine setup_types_for_ghost_send

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine copies the potential from the distribution into the the phic
!> array. If nlphi = .false. phic is set to zero 
!-----------------------------------------------------------------------------
subroutine get_phi(fdis,phic)
  use general, only : gkw_abort
  use control, only : nlphi, flux_tube
  use grid, only : nmod, nx, ns
  use index_function, only : indx
  complex, dimension(:), intent(in) :: fdis
  !> the potential will be put into it the array passed. This
  !> only works if the array has lbounds 1!
  complex, dimension(nmod,1-ghost_points_x:nx+ghost_points_x,&
     & 1-ghost_points_s:ns+ghost_points_s), intent(out) :: phic
  integer :: ix, i, imod, ghost_s, ghost_x

  if (first_call) then
    call gkw_abort('get_phi: you can not call this before dist_init')
  end if

  if (nlphi) then 

    ! copy phi
    if(size(fdis) == nsolc .or. .not.flux_tube) then
      ! the indx(..) function throws a bug if I just ask for ghost
      ! cells like this with x AND s parallelisation at the same time.
      
      ! without ghostcells
      ghost_s = 0
      ghost_x = 0
    else
      ! with ghostcells
      ghost_s = ghost_points_s
      ghost_x = ghost_points_x
    end if
    do i = 1-ghost_s, ns+ghost_s
      do ix = 1-ghost_x, nx+ghost_x
        do imod = 1, nmod 
          phic(imod,ix,i) = fdis(indx(iphi,imod,ix,i))
        end do
      end do
    end do
  else
    ! phi is not solved for and therefore set to zero
    phic = (0.E0,0.E0)
  end if

end subroutine get_phi

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine copies the potential from the distribution into the the aparc
!> array. If nlapar = .false. aparc is set to zero 
!-----------------------------------------------------------------------------
subroutine get_apar(fdis,aparc)
  use general, only : gkw_abort
  use control, only : nlapar
  use grid, only : nmod, nx, ns
  use index_function, only : indx
  complex, dimension(:), intent(in) :: fdis
  complex, dimension(nmod,1-ghost_points_x:nx+ghost_points_x,&
     & 1-ghost_points_s:ns+ghost_points_s), intent(out) :: aparc
  integer :: ix, i, imod, ghost_s, ghost_x

  if (first_call) then
    call gkw_abort('get_apar: you can not call this before dist_init')
  end if

  if (nlapar) then 

    ! copy apar
    if(size(fdis) == nsolc) then
      ! without ghostcells
      ghost_s = 0
      ghost_x = 0
    else
      ! with ghostcells
      ghost_s = ghost_points_s
      ghost_x = ghost_points_x
    end if
    do i = 1-ghost_s, ns+ghost_s
      do ix = 1-ghost_x, nx+ghost_x
        do imod = 1, nmod 
          aparc(imod,ix,i) = fdis(indx(iapar,imod,ix,i))
        end do
      end do
    end do
  else
    ! apar is not solved for and therefore set to zero
    aparc = (0.E0,0.E0)
  end if

end subroutine get_apar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine copies the potential from the distribution into the the bparc
!> array. If nlbpar = .false. bparc is set to zero 
!-----------------------------------------------------------------------------
subroutine get_bpar(fdis,bparc)
  use general, only : gkw_abort
  use control, only : nlbpar
  use grid, only : nmod, nx, ns
  use index_function, only : indx
  complex, dimension(:), intent(in) :: fdis
  complex, dimension(nmod,1-ghost_points_x:nx+ghost_points_x,&
     & 1-ghost_points_s:ns+ghost_points_s), intent(out) :: bparc
  integer :: ix, i, imod, ghost_s, ghost_x

  if (first_call) then
    call gkw_abort('get_bpar: you can not call this before dist_init')
  end if

  if (nlbpar) then 

    ! copy bpar
    if(size(fdis) == nsolc) then
      ! without ghostcells
      ghost_s = 0
      ghost_x = 0
    else
      ! with ghostcells
      ghost_s = ghost_points_s
      ghost_x = ghost_points_x
    end if
    do i = 1-ghost_s, ns+ghost_s
      do ix = 1-ghost_x, nx+ghost_x
        do imod = 1, nmod 
          bparc(imod,ix,i) = fdis(indx(ibpar,imod,ix,i))
        end do
      end do
    end do
  else
    ! bpar is not solved for and therefore set to zero
    bparc = (0.E0,0.E0)
  end if

end subroutine get_bpar


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine allocates all the arrays connected with dist.
!-----------------------------------------------------------------------------
subroutine dist_allocate()

  use general,    only : gkw_abort
  use grid,       only : ns, nmu, nvpar, nx
  use grid,       only : nmod, nsp, lsendrecv_x
  use components, only : energetic_particles
  use control,    only : laverage_dist_over_time

  !> integer for error status
  integer :: ierr, i
  ! initialize ierr
  ierr= 0

  ! allocate the array that contains the Maxwell with ghost values
  allocate(fmaxwl(nx,-1:ns+2,0:nmu+1,-1:nvpar+2,nsp),stat=ierr)
  if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate fmaxwl')
  fmaxwl(:,:,:,:,:) = 0.0

  if (energetic_particles) then
    ! allocate the array that contains the EP without ghost values
    allocate(f_EP(nx,ns,nmu,nvpar,nsp),stat=ierr)
    if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate f_EP')
    f_EP(:,:,:,:,:) = 0.0

    ! allocate the array that contains the derivative of EP without ghost values
    allocate(df_EPdv(nx,ns,nmu,nvpar,nsp),stat=ierr)
    if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate df_EPdv')
    df_EPdv(:,:,:,:,:) = 0.0

    ! allocate the array that contains the sinhc part
    ! of the derivative of EP without ghost values
    allocate(df_EPdv_sinhc(nx,ns,nmu,nvpar,nsp),stat=ierr)
    if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate df_EPdv_sinhc')
    df_EPdv_sinhc(:,:,:,:,:) = 0.0
  end if

  ! allocate the array that contains the alpha particle slowing down
  ! distribution
  allocate(falpha(ns,nmu,nvpar),stat=ierr)
  if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate falpha')

  ! allocate the array that contains phi
  allocate(phi(nmod,1-ghost_points_x:nx+ghost_points_x,&
     & 1-ghost_points_s:ns+ghost_points_s),stat=ierr)
  if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate phi')

  ! allocate the array that contains apar 
  allocate(apar(nmod,1-ghost_points_x:nx+ghost_points_x,&
     & 1-ghost_points_s:ns+ghost_points_s),stat=ierr)
  if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate apar')

  ! allocate the array that contains bpar 
  allocate(bpar(nmod,1-ghost_points_x:nx+ghost_points_x,&
     & 1-ghost_points_s:ns+ghost_points_s),stat=ierr)
  if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate bpar')

  ! allocate the distribution function 
  allocate(fdisi(nsolc),stat=ierr)
  if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate fdisi')
  !$omp parallel do schedule(static)
  do i = 1, nsolc
    fdisi(i) = 0.0
  end do
  !$omp end parallel do
  
  ! allocate buffer for temporally averaged distribution
  if(laverage_dist_over_time) then
    allocate(fdisi_tavg(nsolc),stat=ierr)
    !$omp parallel do schedule(static)
    do i = 1, nsolc
      fdisi_tavg(i) = 0.0
    end do
    !$omp end parallel do
  end if

  ! allocate tmp space for ghost cell communcations
  allocate(fdis_tmp(msolc),stat=ierr)
  if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate fdis_tmp')
  !$omp parallel do schedule(static)
  do i = 1, msolc
    ! Set fdis_tmp initially to a value that should affect the result if ghost
    ! points are incorrectly referenced (for example, if elem_is_on_vpar_grid() or
    ! connect_parallel were broken).
    fdis_tmp(i) = (4352931.,-9876104.)
  end do
  !$omp end parallel do

  ! allocate tmp2 space for radial ghost cell communications
  if (lsendrecv_x) then
    allocate(fdis_tmp2(msolc),stat=ierr)
    if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate fdis_tmp2')
    !$omp parallel do schedule(static)
    do i = 1, msolc
      fdis_tmp2(i) = (5352931.,-7876104.)
    end do
    !$omp end parallel do
  end if

end subroutine dist_allocate

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine dist_deallocate()
  if(allocated(fmaxwl)) deallocate(fmaxwl)
  if(allocated(f_EP)) deallocate(f_EP)
  if(allocated(df_EPdv)) deallocate(df_EPdv)
  if(allocated(df_EPdv_sinhc)) deallocate(df_EPdv_sinhc)
  if(allocated(falpha)) deallocate(falpha)
  if(allocated(phi)) deallocate(phi)
  if(allocated(apar)) deallocate(apar)
  if(allocated(bpar)) deallocate(bpar)
  if(allocated(fdisi)) deallocate(fdisi)
  if(allocated(fdis_tmp)) deallocate(fdis_tmp)
  if(allocated(fdis_tmp2)) deallocate(fdis_tmp)
  if(allocated(fdisi_tavg)) deallocate(fdisi_tavg)
end subroutine dist_deallocate


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Creates a new datatype to be used for communicating gathered array in 3d
!> (NOT YET) used in the non-spectral polarisation gather
!-----------------------------------------------------------------------------
subroutine create_gather_datatype(type_new)

  use grid,         only : iproc_x, nmod, ns, nx, n_x_grid
  use mpiinterface, only : MPICOMPLEX_X
#if defined(mpi2)
  use mpiinterface, only : MPI_ORDER_FORTRAN
#endif

  integer :: type_old,type_new,ierr
  integer, dimension(3) :: whole_array_size, sub_array_size, starts

  type_old = MPICOMPLEX_X  ! communicated array is a real array twice the size

  ! this depends on the index order
  whole_array_size = (/nmod,n_x_grid,ns/)
  sub_array_size   = (/nmod,nx,ns/)
  starts           = (/1,   iproc_x*nmod*nx, 1/)

#if defined(mpi2)

  call MPI_Type_create_subarray(3, whole_array_size, sub_array_size, starts, &
                        & MPI_ORDER_FORTRAN, type_old, type_new, ierr)

  call MPI_TYPE_COMMIT(type_new,ierr)

#else

  ! not very good, but there are no obvious better choices
  type_new = type_old

#endif

end subroutine create_gather_datatype


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Creates a new datatype for for communicating ghost points. The input is
!> an array of indices `ghost_offset', and its length, `ghost_size'. The old
!> type is assumed to be complex.
!-----------------------------------------------------------------------------
subroutine get_ghost_block_type(ghost_offset,type_new)

  use global,       only : root_and_verbose
  use mpiinterface, only : MPICOMPLEX_X
  
  integer, intent(in), dimension(:) :: ghost_offset
  !> The returned `type_new' will be of the MPI indexed-block type
  !> (can be non-contiguous, though) of complex numbers.
  integer, intent(out) :: type_new

  integer :: blocklen,ierr, ghost_size
  ghost_size = size(ghost_offset)

  blocklen = 1

  if (root_and_verbose) then
    write (*,*) '* creating indexed block datatype of length',ghost_size
  end if

#if defined(mpi2)
  ! reduce offset by 1 for MPI call (offsets start at 0)
  call MPI_TYPE_CREATE_INDEXED_BLOCK(ghost_size,blocklen,ghost_offset-1,       &
      & MPICOMPLEX_X,type_new,ierr)
  call MPI_TYPE_COMMIT(type_new,ierr)

#else

  ! not very good, but there are no obvious better choices
  type_new = MPICOMPLEX_X

#endif

end subroutine get_ghost_block_type

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
function combine_nmat_factors(id1,id2,id3) result(ret)
  integer, intent(in) :: id1,id2
  integer, intent(in),optional :: id3
  real :: ret
  ret = 0
  ret = ret + nmat_factor_for_deriv(id1)-1
  ret = ret + nmat_factor_for_deriv(id2)-1
  if(present(id3)) then
    ret = ret + nmat_factor_for_deriv(id3)-1
  end if
  ! all the finite differences share the central element
  ret = ret + 1
end function combine_nmat_factors

end module dist
