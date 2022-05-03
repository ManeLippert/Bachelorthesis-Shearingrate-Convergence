!----------------------------------------------------------------------------
!> This module contains variables related to toroidal rotation
!> And routines and variables associated with perpendicular shearing.
!> Also contains the routine to calculate the centrifugal potential
!> This module reads the parameters in the ROTATION namelist.
!----------------------------------------------------------------------------
module rotation

  use global, only : lenswitch, root_and_verbose

  implicit none

  private

  public :: rotation_read_nml, rotation_write_nml, rotation_check_params
  public :: rotation_bcast_nml, vcor, ts_uprim, need_fft, grad_pot
  public :: rotation_init, rotation_allocate, parallelize_rotation 
  public :: toroidal_shear, kxshift, kxshift_write_nml
  public :: shear_real, shear_remap, shear_ky_shift, kxshift_bcast_nml
  public :: shear_shift_ky, wavevector_remap, kxshift_read_nml
  public :: coriolis, cf_drift, cf_trap, cfen, cfen_G, rotation_parallelized
  public :: cfen_hfs, cfen_lfs
  public :: dcfen_ds, dcfphi_dpsi, dcfphi_ds, cf_upsrc, cfenh, cfenl
  public :: shear_periodic_bcs, shear_rate

  interface rotation_write_nml
    module procedure rotation_read_nml
  end interface

  interface kxshift_write_nml
    module procedure kxshift_read_nml
  end interface
  
  !> The rotation of the plasma vcor = Vtor / vthref
  !> Will be parallel to toroidal magnetic field if signB=1 
  real, save :: vcor 

  !> Normalised shearing rate for the perp shear added in non linear terms
  !> Only used if perp_shear = true
  !> shear_rate is internally rescaled
  real, save :: shear_rate

  !> uprim when toroidal_shear='use_shear_rate'
  real, save :: ts_uprim=0.E0

  !> time at which ExB shearing begins.
  real, save :: t_shear_begin=0.E0

  !> Selects the options for the perpendicular shear flow
  character (len = lenswitch), save :: shear_profile

  !> Select input overrides for toroidal rotation.
  character (len = lenswitch), save :: toroidal_shear

  logical, save :: perp_shear           !< true if perpendicular shear included
  logical, save :: shear_real           !< true if shear flow added in real space
  logical, save :: shear_remap          !< true for shearing by wavevector remap
  logical, save :: shear_ky_shift       !< true for shearing by shifting theorem in ky
  logical, save :: need_fft             !< true if fft is used
  logical, save :: coriolis             !< true if coriolis drift term is kept
  logical, save :: cf_drift             !< true if centrifugal drift terms are kept
  logical, save :: cf_trap              !< true if centrifugal trapping is kept
  logical, save :: cf_upphi             !< add uprim correction to centrifugal phi
  logical, save :: cf_upsrc             !< add uprim correction to centrifugal source
  logical, save :: cf_qncheck           !< enforce quasineutrality in gradients

  !>Array that tracks number of kx mode shifts: kxshift(nmod)
  integer, allocatable, dimension(:), save :: kxshift

  !>Indexing arrays for shearing wavevector remap:
  !> aindx(nmod*nmu*nvpar*ns*nsp*(nx-1))
  integer, allocatable, dimension(:), save :: aindx, aindx_shift

  !> integers used in remap loops
  integer, save :: ixstart, ixend, ixedge, ixdir

  !> Array for shear flow function potential gradient. grad_pot(mrad,ns)
  !> used for shear_real
  real, allocatable, save :: grad_pot(:,:)

  !> Array for the centrifugal energy, cfen(-1:n_s_grid+2,nsp+iadia)
  real, allocatable, save :: cfen(:,:), cfeno(:,:)
  !> Maximum and minimum of centrifugal energy, which can be used to
  !> check particle trapping
  real, allocatable, save :: cfen_hfs(:), cfen_lfs(:)
  !> Array for the centrifugal energy Global in sp: cfen_G(ns,number_of_species)
  real, allocatable, save :: cfen_G(:,:)
  !> And derivative along s (n_s_grid,nsp)
  real, allocatable, save :: dcfen_ds(:,:)  
  !> Radial derivative of centrifugal potential (n_s_grid)
  real, allocatable, save :: dcfphi_dpsi(:)
  !> parallel derivative of centrifugal potential (n_s_grid)
  real, allocatable, save :: dcfphi_ds(:)
  !> Exact high and low field side centrifugal energy by species (nsp+iadia)
  real, allocatable, save :: cfenh(:), cfenl(:)
  !> ICRH species energy (n_s_grid)
  real, allocatable, save :: cfen_icrh(:)  
  ! > workspace, cf_phi and adjacent flux surfaces cfphi(n_s_grid,-2:2) 
  ! > NOT a public array
  real, allocatable, save :: cfphi(:,:) 

  !Arrays for the ffts, if required.
  complex, ALLOCATABLE, save :: arr(:,:)
  
  !> Index arrays for the ffts - one needed for both forward and back
  integer,    allocatable, save :: jind_nx(:)

  !>Size of padded fft array
  integer, save ::mx

  !>Should be identical to the ones in non_linear_terms
  integer, save :: mrad

  !>Flag to set once rotation is parallelized
  logical, save :: rotation_parallelized=.false.

contains


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reads the plasma toroidal rotation and perpendicular shearing rate, and 
!> shearing method.
!----------------------------------------------------------------------------

subroutine rotation_read_nml(ilun,io_stat,lwrite)
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  
  namelist / rotation / vcor, shear_rate, shear_profile, toroidal_shear, &
               coriolis, cf_trap, cf_drift, t_shear_begin,  & 
               cf_upphi, cf_upsrc, cf_qncheck
  
  integer, intent(in)  :: ilun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then
      
    ! the default is no rotation, no shearing
    vcor = 0
    shear_rate = 0.E0
    ! Set the type of perpendicular shearing
    shear_profile     =  'none' 
    t_shear_begin     =  0.E0
    toroidal_shear    =  'none'
    coriolis          =  .true.
    cf_drift       =  .false.
    cf_trap        =  .false.
    cf_upphi       = .true.
    cf_upsrc       = .true.
    cf_qncheck     = .true.

    ! read namelist and return on error
    read(ilun,NML=rotation,IOSTAT=io_stat)

    end if
  else
    if(root_processor) write(ilun,NML=rotation)

    call write_run_parameter('rotation', 'vcor', vcor)
    call write_run_parameter('rotation', 'shear_rate', shear_rate)
    call write_run_parameter('rotation', 'shear_profile', shear_profile)
    call write_run_parameter('rotation', 'toroidal_shear', toroidal_shear)
    call write_run_parameter('rotation', 'coriolis', coriolis)
    call write_run_parameter('rotation', 'cf_trap', cf_trap)
    call write_run_parameter('rotation', 'cf_drift', cf_drift)
    call write_run_parameter('rotation', 't_shear_begin', t_shear_begin)
    call write_run_parameter('rotation', 'cf_upphi', cf_upphi)
    call write_run_parameter('rotation', 'cf_upsrc', cf_upsrc)
    call write_run_parameter('rotation', 'cf_qncheck', cf_qncheck)
    
  end if

end subroutine rotation_read_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> broadcast the rotation parameters
!----------------------------------------------------------------------------
subroutine rotation_bcast_nml

  use  mpiinterface, only : mpibcast
  
  call mpibcast(vcor, 1)
  call mpibcast(shear_rate, 1) 
  call mpibcast(coriolis, 1)
  call mpibcast(cf_drift, 1)
  call mpibcast(cf_trap, 1)
  call mpibcast(cf_upphi, 1)
  call mpibcast(cf_upsrc, 1)
  call mpibcast(cf_qncheck, 1)
  call mpibcast(t_shear_begin, 1)
  call mpibcast(toroidal_shear, lenswitch)
  call mpibcast(shear_profile, lenswitch)

end subroutine rotation_bcast_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Run some checks to test for appropriate rotation paremters.
!----------------------------------------------------------------------------

subroutine rotation_check_params

  use control,      only : method, meth, non_linear, flux_tube 
  use control,      only : radial_boundary_conditions, lcollisions, vp_trap 
  use control,      only : shear_periodic, spectral_radius
  use components,   only : zeff_sp, nsps, icrh_model
  use geom,         only : R0_loc, vcor_miller
  use mode,         only : mode_box
  use fft,          only : working_fft_library
  use general,      only : gkw_warn, gkw_abort, gkw_exit
  use mpiinterface, only : root_processor
  use global,only : r_tiny
  
  select case(shear_profile)

    case('none')
      perp_shear=.false.
      shear_real=.false.
      shear_remap=.false.
      shear_ky_shift=.false.
      need_fft=.false.
      if (toroidal_shear.ne.'none') then
        toroidal_shear = 'none'
        call gkw_warn('No E x B shear: uprim not overridden')
      end if

      !Warn if shearing rate is not zero
      if(abs(shear_rate) > 1e-10) then
        call gkw_warn('No E x B shear active')
        shear_rate=0.
      end if

    case('wavevector_remap')
      perp_shear=.true.
      shear_real=.false.
      shear_remap=.true.
      shear_ky_shift=.false.
      need_fft=.false.
      if (root_processor) write(*,*) &
        & 'Shear method: wavevector_remap: Check for nx convergence'
      if (.not. spectral_radius) then
        call gkw_abort('wavevector_remap only with spectral_radius')
      end if

    case('ky_shift')
      perp_shear=.true.
      shear_real=.false.
      shear_remap=.false.
      shear_ky_shift=.true.
      need_fft=.true.
      if(root_processor) then
      if (.not. spectral_radius) then
        write(*,*) 'Shear method: ky_shift: Boundary discontinuity in ExB profile'
        write(*,*) 'Shear profile: Linear sawtooth with boundary discontinuity'
      end if
      end if
 
    case('symmetric')
      perp_shear=.true.
      shear_real=.true.
      shear_remap=.false.
      need_fft=.true.
      shear_ky_shift=.false.
      write(*,*) 'Shear method: Triangle saw symmetric '
      call gkw_warn('Check shear rate normalisation for shear_real')
      if (.not. spectral_radius) then
        call gkw_abort('symmetric only with spectral_radius')
      end if
 
    case('linear')
      perp_shear=.true.
      shear_real=.true.
      shear_remap=.false.
      need_fft=.true.
      shear_ky_shift=.false.
      if (root_processor.and. spectral_radius) then
        write(*,*) 'Shear profile: Linear sawtooth with boundary discontinuity'
        call gkw_warn('Check shear rate normalisation for shear_real')
      end if 

    case default 
      call gkw_exit('Unknown shearing option')
 
  end select
  
  if (.not. spectral_radius .and. perp_shear) then
    if (radial_boundary_conditions == 'periodic') then
      shear_periodic = .true.
      if(.not. flux_tube) call gkw_abort('shear_periodic only with flux_tube')
    else  
      call gkw_warn('Use ExB shear with boundary damping or periodic boundaries')
    end if
  end if 
 
  if (non_linear) need_fft=.true.
  if (mode_box) need_fft=.true.

  if (abs(vcor).gt.1E-6 .and. perp_shear) then
     call gkw_warn('shear rate is defined in the rotating frame')
  endif

  !Toroidal shear not consistent for vcor/=0?
  !Or is this ok in rotating frame?
   select case(toroidal_shear)

     case('use_shear_rate')
          call gkw_warn('uprim inputs ignored for toroidal shear')

     case('use_uprim','add_uprim')
        call gkw_warn('shear_rate input overridden by toroidal shear')

     case('none') !Default
        !Do nothing
     
     case default 
       call gkw_exit('Unknown toroidal_shear option')

   end select

  if(cf_trap.and.lcollisions) then
     call gkw_warn('cf_trap: Parallel variation in collision frequency')
     call gkw_warn('cf_trap: Collision freq. taken at R=R_0=R_'//trim(R0_loc))
     call gkw_warn('cf_trap: Parallel variation not kept in Coloumb logarithm')
  end if
  
  ! The non uniform velocity grid is not yet tested by species
  if(vp_trap /= 0 .and. cf_trap .and. (nsps > 2 .or. zeff_sp > 1.)) then
    call gkw_warn('Centrifugal trapping + vp_trap + impurities not yet benchmarked')
  end if
  
  !Run some checks to test for appropriate control paremters
  if (shear_remap) then
    select case (method)
       case('EXP')
          select case(meth)
             case(3) 
                call gkw_exit('Shear not implemented for 3rd order scheme')
          end select
       case default ; call gkw_exit('rotation: This shearing option only '//&
                          &          'works with explicit time integration')
    end select
  end if

  ! The rotating frame is not strictly valid globally
  ! But at the very least the cf_max parameters in componnents need fixing
  ! Also, most of the cf solver uses the species quantities at ix=1 only
  if(.not. flux_tube .and. cf_trap) then
    call gkw_exit('Centrifugal effects only valid in flux tube')
  end if

  !Needs to be tested after reading mode namelist....
  if(perp_shear.and..not.mode_box) then
     call gkw_exit('Shearing only availiable with mode_box')
     !This restriction could be lifted only if mode_box=false 
     !is made to work with multiple nx modes.
  end if

  if (.not. working_fft_library.and.need_fft) then
      call gkw_exit('(shear?) option requires a working FFT library')
  end if

  ! Needed if rotation effects are taken into account in the magnetic equilibrium
  ! through Miller geometry
  vcor_miller = vcor
  
  if(abs(vcor) < r_tiny .and. icrh_model == 0) then
    coriolis = .false.
    cf_drift = .false.
    cf_trap = .false.
    cf_upphi = .false.
    cf_upsrc = .false.
    !cf_qncheck = .false.
  end if
  
end subroutine rotation_check_params 


!--------------------------------------------------------------------
!> Adds ExB shear as using 1D FFt only - boundary discontinuity
!> Uses FFT shifting theorem
!> This routine is significantly faster than adding using shear_real
!> to add the shear as a term in add_non_linear_terms
!> NOT yet optimised
!--------------------------------------------------------------------
subroutine shear_shift_ky(inout)

  use fft,            only : fourcol, FFT_INVERSE, FFT_FORWARD
  use control,        only : dtim, spectral_radius
  use grid,           only : n_x_grid, nx, nmod, nsp, nmu
  use grid,           only : nvpar, ns, lx, gx 
  use mode,           only : krho
  use dist,           only : ifdis, nsolc
  use geom,           only : dxgr
  use index_function, only : indx
  use general,        only : gkw_abort
  use constants,      only : ci1 
  
  complex,dimension(nsolc),intent(inout) :: inout

  integer :: ix,iky, ipar, is, jv, kt, imod, ixplus
  real    :: dx, dum2
  
  if(.not.shear_ky_shift) call gkw_abort('Invalid call to shear_shift_ky')
  
  !arr=0.E0
  
  !Possibly could drop repeated nyquist mode and do nx-1 sized fft. 
  !But check they are the same first.
  if (spectral_radius) then
    dx = lx /(n_x_grid+1)
    ixplus = 1
  else
    dx = dxgr
    ixplus = 0
  end if
    
  dum2=-shear_rate*dx*dtim
  
  do ipar = 1, ns; do is = 1, nsp; do jv = 1, nmu; do kt = 1, nvpar         
    do ix=1,nx
        do imod=1, nmod
          arr(jind_nx(ix),imod)=inout(indx(ifdis,imod,ix,ipar,jv,kt,is))
        end do
    end do

    if (spectral_radius) then 
      call fourcol(arr, FFT_INVERSE)
      arr = arr / ((n_x_grid+1))
    end if

    do iky=1,nmod
        !Add the shear by shifting in ky.
        do ix=1,nx + ixplus
            arr(ix,iky)=exp(ci1*krho(iky)*(gx(ix)-n_x_grid/2)*dum2)*arr(ix,iky)
        end do
    end do

    !Forward FFT of 2D array in x only
    !Intent (INOUT, For dir -1)
    if (spectral_radius) call fourcol(arr, FFT_FORWARD)

    !Do the inverse rearrangment
    do ix=1,nx
      do imod=1, nmod
        inout(indx(ifdis,imod,ix,ipar,jv,kt,is))=arr(jind_nx(ix),imod)
      end do
    end do
    
  end do; end do; end do; end do

end subroutine shear_shift_ky


!--------------------------------------------------------------------
!> This subroutine applies shear periodic boundary condiions to the
!> ghost cells in both f and all fields
!> Not yet optimised
!--------------------------------------------------------------------
subroutine shear_periodic_bcs(fdis,ipart,dpart)

 use control,        only : spectral_radius, shear_periodic
 use control,        only : time, nlapar
 use grid,           only : lsendrecv_x, n_x_grid, nx, nmod, nsp, nmu 
 use grid,           only : nvpar, ns, lx, gx 
 use dist,           only : msolc, ifdis, iphi, iphi_ga, iapar 
 use dist,           only : iapar_ga, xgp => ghost_points_x
 use mode,           only : krho 
 use index_function, only : indx
 use general,        only : gkw_abort
 use constants,      only : ci1
 use global,         only : compiled_with_mpi
 
 complex, intent(inout) :: fdis(:)
 integer, intent(in) :: ipart
 real, intent(in)    :: dpart
 real :: shift
 
 integer :: i, j, k, is, ix, imod, idum

 if (.not. shear_periodic) return
 
 if (size(fdis).ne.msolc) call gkw_abort('shear_periodic_bcs: wrong size fdis')
 if (.not. lsendrecv_x) call gkw_abort('Need x ghost cells (mpi) for shear periodic boundaries')
 if (spectral_radius) call gkw_abort('Need non spectral_radius for shear periodic boundaries')
 ! Require the ghost cells which are populated by MPI (no alternative at present)
 if (.not.compiled_with_mpi) call gkw_abort('Need MPI compilation for shear periodic boundaries')
 if (.not. nlapar .and. (ipart == iapar .or. ipart== iapar_ga)) return
 
 shift = -(time+dpart) * shear_rate * lx
     
  lmod: do imod = 2, nmod 
    ! left boundary
    if (gx(1) == 1) then
      x: do ix = 1-xgp, 0 
    
        if (ipart == ifdis) then
          do i = 1, ns; do is = 1, nsp; do j = 1, nmu; do k = 1, nvpar; 
            idum = indx(ifdis,imod,ix,i,j,k,is)
            fdis(idum)=fdis(idum)*exp(-ci1*krho(imod)*shift)
          end do; end do; end do; end do
        end if 
   
        if (ipart == iphi_ga .or. ipart == iapar_ga) then
          do i = 1, ns; do is = 1, nsp; do j = 1, nmu 
            idum = indx(ipart,imod,ix,i,j,is)
            fdis(idum)=fdis(idum)*exp(-ci1*krho(imod)*shift)
          end do; end do; end do
        end if
        
        if (ipart == iphi .or. ipart == iapar) then
          do i = 1, ns
            idum = indx(ipart,imod,ix,i)
            fdis(idum)=fdis(idum)*exp(-ci1*krho(imod)*shift)
          end do
        end if
        
      end do x
    end if
    
    ! right boundary
    if (gx(nx) == n_x_grid) then
      x2: do ix = nx+1, nx+xgp
        
        if (ipart == ifdis) then
          do i = 1, ns; do is = 1, nsp; do j = 1, nmu; do k = 1, nvpar; 
            idum = indx(ifdis,imod,ix,i,j,k,is)
            fdis(idum)=fdis(idum)*exp(ci1*krho(imod)*shift)
          end do; end do; end do; end do
        end if
    
        if (ipart == iphi_ga .or. ipart == iapar_ga) then
          do i = 1, ns; do is = 1, nsp; do j = 1, nmu 
            idum = indx(ipart,imod,ix,i,j,is)
            fdis(idum)=fdis(idum)*exp(ci1*krho(imod)*shift)
          end do; end do; end do
        end if
        
        if (ipart == iphi .or. ipart == iapar) then
          do i = 1, ns
            idum = indx(ipart,imod,ix,i)
            fdis(idum)=fdis(idum)*exp(ci1*krho(imod)*shift)
          end do
        end if
        
      end do x2
    end if 
  end do lmod
      
end subroutine shear_periodic_bcs


!--------------------------------------------------------------------
!> This routine is used for setting up an index array
!--------------------------------------------------------------------
subroutine shear_remap_init

  use grid,           only : nmod, nx, ns, nmu, nvpar, nsp
  use dist,           only : ifdis, nf
  use index_function, only : indx
  use control,        only : spectral_radius
  use general,        only : gkw_abort
  use mpiinterface,   only : root_processor

  integer ::  imod, ix, jv, kt, is, ipar

  integer, save :: idx = 1

  if (.not. spectral_radius) return

  if (.not.shear_remap) call gkw_abort('Invalid call to shear_remap_init')

    kxshift(:)=0
    aindx=0
    aindx_shift=0

    if (shear_rate.gt.0) then 
        ixdir=1
        ixstart=1
        ixend=nx-1
        ixedge=nx
    else if (shear_rate.lt.0) then
        ixdir=-1
        ixstart=nx
        ixend=2
        ixedge=1
    else          !This might happen if the shear rate is zero
        shear_remap=.false.
        if (root_processor) write(*,*) 'Shear_remap off'
        return    !The index arrays will be blank but should never be used
    end if

  idx=1
  !Setup the indexing arrays for remapping
  !Based on the assumption that index function increases monotonically with ix
  !Direction is reversed for negative shearing rates
  !The order of the loops must exactly match that in wavevector_remap
  do imod = 1, nmod               !Loop over poloidal modes
    do ipar = 1, ns               !Loop along feild line points
      do jv = 1, nmu              !Loop over perpendicular velocity
        do kt = 1, nvpar          !Loop over parallel velocity
          do is = 1, nsp          !Loop over species
            do ix = ixstart, ixend, ixdir       !Loop over radial modes
              aindx(idx) = indx(ifdis,imod,ix,ipar,jv,kt,is)
              aindx_shift(idx)= indx(ifdis,imod,ix+ixdir,ipar,jv,kt,is)
                idx=idx+1         
            end do
          end do
        end do
      end do
    end do
  end do

  !Perform check
  !write(*,*) nf, nmod*ns*nvpar*nmu*nsp, idx
  if (idx.ne.nf-(nmod*ns*nmu*nvpar*nsp)+1) then
    call gkw_abort('Severe error in wavevector remap init: idx')
  end if

end subroutine shear_remap_init


!--------------------------------------------------------------------
!> Perpendicular shearing by remapping of spectral kx wavevector grid
!>
!> This ExB shear implementation depends on the history of kx shifts:
!> On restart, these are read from FDS.dat, in subroutine
!> kxshift_read_nml
!--------------------------------------------------------------------
subroutine wavevector_remap(fdis,t)

  use dist,           only : ifdis, nsolc, nf
  use index_function, only : indx
  use grid,           only : nmod, nx, ns, nmu, nvpar, nsp
  use control,        only : t_init
  use mode,           only : krho, kxspace
  use control,        only : spectral_radius
  use general,        only : gkw_abort, gkw_exit
  use mpiinterface,   only : root_processor

  complex, intent(inout) :: fdis(nsolc)
  real, INTENT(IN) :: t !Time
  real :: shift, tdum
  integer :: dkx, ix, idx=1, imod, jv, kt, is, ipar
  logical, save :: first_call=.true.

  if(.not. shear_remap .or. .not. spectral_radius) then
     call gkw_abort('exp_integration: Invalid call to wavevector_remap')
  end if

  ! Check for restart which adds shear, and set t_shear_begin correctly
  ! Unfortunately this cannot be in shear_remap_init
  ! since t_init is set after reading of restart files 
  if (first_call) then
       !Case of no previous shear
       if (maxval(abs(kxshift)).eq.0) then
          if (t_shear_begin.lt.t_init) t_shear_begin=t_init
       else !previous shear
          !t_shear_begin should have been read correctly from input.dat or FDS.dat
       end if
       first_call=.false.
  end if

  idx=1

  tdum=t-t_shear_begin
  !At negative time there is no shear.
  if (tdum.lt.0.) tdum=0.E0

  do imod = 1, nmod       !Loop over "poloidal" modes

     shift=krho(imod)*shear_rate*(tdum)
     !write(*,*) shift

     !If the wavevectors need to be remapped
     if (nint(shift/kxspace).ne.kxshift(imod)) THEN

        !Calculate the number of kx grid points to shift the wavevectors
        dkx=nint(shift/kxspace)-kxshift(imod) 

        ! The ExB shear implementation depends on the history of the shifts
        ! And for fidelity should not allow shifts > 1 
        if (abs(dkx).gt.1) then
           if (root_processor) then
             write(*,*)
             write(*,*) 'ERROR: ExB shear: shift request ', abs(dkx), ' kx places'
             write(*,*) 'Possible causes, if shift request is >> 1 at restart'
             write(*,*) ' 1) Discontinuity: No FDS.dat restart data (required)'
             write(*,*) ' 2) Discontinuity: Changed shear_rate on restart'
             write(*,*) '   -> remove incompatible FDS.dat, set t_shear_begin >= time'
             write(*,*) ' 3) Discontinuity: Incorrect t_shear_begin for restart'
             write(*,*) '   -> Set t_shear_begin consistently in FDS.dat and input.dat'
             write(*,*) 'Possible causes, if shift request is > 1'
             write(*,*) ' 4) Time resolution insufficient for shearing rate'
             write(*,*) '    -> Decrease timestep, or increase kx spacing'
            end if
            call gkw_exit('Invalid shift in rotation:wavevector_remap')
        end if

        if (ixdir.ne.dkx) then
            write(*,*) "kxspace", kxspace
            write(*,*) "t, imod, ixdir, dkx", t, imod, ixdir, dkx
            write(*,*) 'kxshift', kxshift
            call gkw_abort("Severe error in rotation:wavevector_remap")
        end if

        do ipar = 1, ns                      !Loop along feild line points
          do jv = 1, nmu                     !Loop over perpendicular velocity
             do kt = 1, nvpar                !Loop over parallel velocity
                do is = 1, nsp               !Loop over species
                    !Shift the wavevectors
                    do ix = ixstart, ixend, ixdir    !Loop over radial modes
                        fdis(aindx(idx))=fdis(aindx_shift(idx))
                        idx=idx+1
                    end do !nx
                    !Boundary condition in kx space
                    !This call to index function could also be removed
                    fdis(indx(ifdis,imod,ixedge,ipar,jv,kt,is))=0.E0
                end do !nsp
             end do  !vpar
          end do !nmu
        end do !ns

        kxshift(imod)=kxshift(imod)+dkx

      ELSE
         idx=idx+ns*nmu*nvpar*nsp*(nx-1)
      endIF

   end do !nmod

 !Perform check
 !write(*,*) idx, nf-(nmod*ns*nvpar*nmu*nsp)
 if (idx.ne.nf-(nmod*ns*nvpar*nmu*nsp)+1) then
     call gkw_abort('Severe error in wavevector remap: idx')
 end if

end subroutine wavevector_remap


!--------------------------------------------------------------------
!> This subroutine calculates the poloidal 
!> asymmetries from strong rotation and ICRH poloidal asymmetry.
!> The outputs are the centrifugal energy (cfen), and various 
!> derivatives of the centrifugal potential (the potential itself
!> should not be used anywhere).
!> In the case of iterate_fsa, this may be called multiple times.
!>
!> The bisection method is used to find the centrifugal potential
!> for each point in s, and for two adjacent flux surfaces 
!> to calculate radial derivatives from a finite difference.
!>
!> Must be called BEFORE parallelize_geom and parallelize_rotation:
!> In this routine, s is not yet parallel.
!> The quantities calculated by this routine are saved in arrays 
!> which are species-local, but since these are later gathered by the
!> centrifugal diagnostics, this is now a bit backwards and it might 
!> now be neater to parallelise the species at the end; presently
!> this routine and dependencies cf_* must not use the species-local
!> density and gradients as input, as they are not updated until 
!> after the iterate_fsa loop (dangerous, could be improved).
!> TO DO: Pass data in/out/down with arguments, for clearer intent.
!--------------------------------------------------------------------
subroutine centrifugal_energy

 use components, only : cf_quasineutral, cf_phi_max, cf_mass_weight
 use components, only : mas, signz, tmp, iadia, vp_G
! use components, only : icrh_fmin, signz_G, icrh_rln
 use components, only : icrh_params, icrh_norm, isg_icrh
 use components, only : icrh_model
 use geom, only : jfun, jfunh, jfunl, kfun, geom_parallelized, bn_G 
 use geom, only : dBdpsi, R0_loc, bmin, gradp_type, iterate_fsa, isg_lfs
 use grid, only : n_s_grid, nsp, lsp
 use mpicomms, only : COMM_SP_NE
 use general,  only : gkw_abort, gkw_warn
 use mpiinterface, only : mpiallreduce_sum_inplace, root_processor

 real :: hi, lo, mid, a, b, fa, fmid, icrh_eta, BR0 = 0.0, bmag
 real :: accuracy, mass_weight, dpsi, jfunc, vtor !dummies
 real :: coeff4d1(-2:2)
 integer :: i, is, it, ip, ix 
 integer, parameter :: max_iterations=32000 !should fit normal integer
 logical :: converged

 ! set ix and do all calculations for one surface only (at the moment)
 ix = 1 

 !Should always be initialised as used in field equations and trapping
 cfen(:,:)=0.E0
 dcfen_ds(:,:)=0.E0
 dcfphi_dpsi(:)=0.E0
 dcfphi_ds(:)=0.E0
 cfphi(:,:)=0.E0

 if(geom_parallelized.or.rotation_parallelized) then 
   call gkw_abort('centrifugual_energy: call before parallelize')
 endif

 !write (*,*) icrh_model, isg_icrh, icrh_params

 ! Some checks for the ICRH case
 if (icrh_model > 0) then 
   if (.not. cf_trap) call gkw_abort('Set cf_trap = T to use an anisotropic minority')
   if (gradp_type=='rota_miller') call gkw_abort('Cannot use rota_miller with anisotropic minority')
   !if (R0_loc == 'LFS') call gkw_warn('R0_loc = LFS is imprecise for anisotropic minority gradients')
   if (R0_loc == 'axis') call gkw_warn('R0_loc = axis assumes Baxis=Bref for anisotropic minority gradients')
 end if  
   
 if(.not. cf_trap) return ! put back call to centrifugal diagnostics ?

 accuracy=10.E0**(-precision(1.E0))

 !Converged result should not depend on this value in limit of small dpsi
 !However numerical errors amplify if the value is too small
 !Testing in double precision showed 1e-4 to be optimal
 !Gives fractional error from analytic solution of 1e-12
 dpsi=1.E-4

 if (iterate_fsa) dpsi=1.E-8

 if (precision(1.E0)<10) then 
   call gkw_warn('Accuracy of radial derivative of centrifugal potential &
                  &  reduced in single precision')
   dpsi=1.E-3
 end if 

 !Setup finite differencing coefficients.  4th order 1st derivative
 coeff4d1=(/1.E0,-8.E0,0.E0,8.E0,-1.E0/)/(12.E0*dpsi)
 
 ! Be swedish: no radial gradients of any kind
 ! Need to fix the analytic bits to match this...
 !if (abs(icrh_params(3)) > 999.0) then
 !  coeff4d1=(/0.0,0.0,0.0,0.0,0.0/)
 !end if

 !write(*,*) 'accuracy', accuracy, 'precision', precision(1.E0), 'dpsi', dpsi
 !write(*,*) 'cf_phi_max', cf_phi_max
  
 ! MOVE TO A SINGLE INSTANCE 
 select case(R0_loc)
  ! Radial derivative approximate for even n_s_grid.
  ! however, the differences are small
  case('LFS'); icrh_norm = bmin(ix)
  ! Assumes Bref = Baxis (could run a check on mean(bng))
  case('axis'); icrh_norm = 1.0
  case default; call gkw_abort('Error in R0_loc')
 end select 
 
 if (root_processor.and..not.iterate_fsa) then
   write(*,*)
   write(*,*) 'Solving for centrifugal potential...'
 end if
 
 !Dummy outside s-loop
 mass_weight=cf_mass_weight() 

  n_s_loop: do i=1,n_s_grid

       psi_points: do ip=-2,2
       
          ! MOVE TO A SINGLE INSTANCE
          select case(R0_loc)
            ! Radial derivative approximate for even n_s_grid.
            case('LFS'); BR0 = icrh_norm + ip*dpsi*dBdpsi(1,isg_lfs)
            ! Assumes Bref = Baxis (could run a check on mean(bng))
            case('axis'); BR0 = icrh_norm
            case default; call gkw_abort('Error in R0_loc')
          end select

          bmag  = bn_G(1,i) + dBdpsi(1,i)*ip*dpsi
          jfunc =jfun(ix,i) + kfun(ix,i)*ip*dpsi
          
          
          !FJC need to add some warning about the use of vp_G(1)
          
          if (cf_upphi) then
            vtor=vcor-vp_G(1,1)*ip*dpsi
          else
            vtor=vcor
          end if
                
          !Bisection method
          !Assumption: One of f(a),f(b) is > 0 and the other is < 0
          a=((cf_phi_max(3)+cf_phi_max(4))*vcor*vcor*abs(jfunc))+0.1

          !old estimate, deprecated
          !a=(cf_phi_max(1)+log(abs(vcor*vcor*jfunc)+1.E0))
          !a=a*cf_phi_max(2)+0.1
          b=-a
    
          converged=.false.
          
          !Calculate the function Q(Phi(psi,theta)) start value
          fa=cf_quasineutral(a,vtor*vtor*jfunc,bmag/BR0,ip*dpsi)
    
          if (fa < 0.E0) then
              lo = a
              hi = b
          else
              lo = b
              hi = a
          endif
        
          iterate: do it=1,max_iterations

            mid = lo + (hi-lo)/2.E0
            !Check for convergence
            if(abs(lo-mid) < accuracy.and.abs(hi-mid) < accuracy) then
              !if (abs(fmid) .lt. accuracy) then 
              converged=.true.
              exit iterate
            endif

            !Get midpoint value Q(Phi(psi,theta))
            fmid=cf_quasineutral(mid,vtor*vtor*jfunc,bmag/BR0,ip*dpsi)

            !Do the bisection
            if (fmid < 0.E0) then
              lo = mid
            else
              hi = mid
            endif
          end do iterate

          if (.not.converged.or.abs(mid).gt.abs(a)) then
            if (root_processor) write(*,*) 'ip', ip, 's',i,'it',it,'f',fmid, &
                              & 'phi-l-h-m', lo, hi, mid, 'upper limit', a
            call gkw_abort('centrifugal potential bisection method failed') 
          else
            cfphi(i,ip)=mid

            if(root_and_verbose) write(*,*) 'phi-point: ', ip,  's-point: ',i,' used ',it,' iterations'
            if(root_and_verbose) write(*,*) '  upper bisection limit: ',a
            if(root_and_verbose) write(*,*)   
          endif
   
          !Now repeat solver for adjacent flux surface psi+dpsi
       end do psi_points

   end do n_s_loop  
   
   ! ICRH combination method: input unadjustments complication ?
   ! recombine the density gradients of bulk and minority: REMOVE ME
   if (isg_icrh == 0 .and. icrh_model > 0 .and. .false.) then    
     call gkw_abort('icrh combination method removed')
!       fp_G(:,1) =  ((de_G(:,1)-icrh_fmin*icrh_params(4)/signz_G(1))*fp_G(:,1)  &
!                    &  +  icrh_fmin*icrh_params(4)*icrh_rln/signz_G(1))/de_G(:,1)
! 
!       ! check the gradient QN again                   
!       if (abs(sum(signz_G(:)*de_G(1,:)*fp_G(1,:))) > 1e-10 ) call gkw_abort('ICRH QN error 2a')
! 
!       de(:,is) = de_G(:,1)
!       fp(:,is) = fp_G(:,1)
!       
!       if (lsp(1) >= 1 .and.  lsp(1) <= nsp) then      
!         de(:,lsp(1))        = de_G(:,1)
!         fp(:,lsp(1))        = fp_G(:,1)
!       end if 
   end if

  cfen_icrh(:) = 0.0
   
  ! Calculate cfen and radial derivative of cf_phi
  do i=1,n_s_grid

    ! Density correction for anisotropic minority
    if (icrh_model > 0) then        
      if (icrh_model == 1) then 
        icrh_eta = icrh_params(2) - 1.0
        cfen_icrh(i) = icrh_eta*log(bn_G(1,i)/icrh_norm)
      else if (icrh_model == 2) then
        a = icrh_params(2) + (1.0 - icrh_params(2))/(bn_G(1,i)/icrh_norm)
        cfen_icrh(i) = log(a)
      end if 
    end if

    do is=1,nsp+iadia
      !Note the factor (1/2) is lost in the normalisation
      cfen(i,is)=(signz(is)*cfphi(i,0)-vcor*vcor*jfun(ix,i)*mas(is))/tmp(1,is)
      !cfen(i,is)=signz(is)*cfphi(i,0)-0.5*vcor*vcor*jfun(i)*mas(is)

      ! keep a copy of the cfen unmodified for icrh for use in diags
      cfeno(i,is) =       cfen(i,is)

      ! Add the anisotropic minority correction to the correct species   
      ! Independant method 
      if (isg_icrh /= 0) then   
        if (is == lsp(isg_icrh).and. is <=nsp) then
          cfen(i,is)   = cfen(i,is) + cfen_icrh(i)
          cfen_icrh(i) = cfen(i,is)
        end if      
      ! Combination method:  minority combined with bulk ions (species 1)
      else if (icrh_model > 0) then
        call gkw_abort('icrh_combination method removed')
!         if (is == lsp(1).and. is <=nsp) then
!           a = icrh_params(4)*icrh_fmin/de_G(1)/signz_G(1)
!           b = (de_G(1)-icrh_params(4)*icrh_fmin/signz_G(1))/de_G(1)
!           
!           ! The following allows z_min /= z_i but still assumes masses and temperatures are equal
!           cfen_icrh(i) =  cfen_icrh(i) + & 
!                           (icrh_params(4)*cfphi(i,0)-vcor*vcor*jfun(ix,i)*mas(is))/tmp(1,is) 
!           cfen(i,is)   = -log(a*exp(-cfen_icrh(i)) + b*exp(-cfen(i,is)))
!         end if
      end if

      cfen_hfs(is) = maxval(cfen(:,is))
      cfen_lfs(is) = minval(cfen(:,is))

    end do !species

    !Radial finite difference of potential, fourth order
    do ip=-2,2
       dcfphi_dpsi(i)=dcfphi_dpsi(i)+cfphi(i,ip)*coeff4d1(ip)
    end do
  end do !n_s_grid
  
  ! All processors keep a copy of cfen_icrh
  call mpiallreduce_sum_inplace(cfen_icrh,n_s_grid,COMM_SP_NE)

  
  !Calculate the parallel derivative of cfen 
  do is=1,nsp
    call parallel_derivative(cfen(1:n_s_grid,is),dcfen_ds(:,is))
  end do

  !And parallel deriviatve of cfphi
  call parallel_derivative(cfphi(:,0),dcfphi_ds(:))

  !Report success of iterative solver
  if(root_processor.and..not.iterate_fsa) then
    write(*,*) '                                   ...done'
    write(*,*)
  end if
  
  ! set exact cfenh and cfenl by using linear dependence on jfunc.
  ! FAILS FOR ICRH ?
  ! could also be done by root finding, or interpolation, or analytically
  ! for hydrogenic plasma.
  cfenh(:) = cfen(1,:) * jfunh / jfun(1,1)
  cfenl(:) = cfen(1,:) * jfunl / jfun(1,1)

end subroutine centrifugal_energy
 
!--------------------------------------------------------------------
!> Diagnostics of centrifugal density variation. 
!> Currently gathers over the species
!> Assumes arrays are parallel in species but global in s
!> If called after rotation_parallelize and geom parallelize
!> would also need to gather over s
!--------------------------------------------------------------------
subroutine centrifugal_diagnostics(phi,fsa_dens,RLNFSA,write_out)
use io,         only : get_free_file_unit, xy_fmt, xy_fmt_long
use io,         only : output_array, ascii_fmt, output_enabled
use mpicomms,   only : COMM_SP_NE, COMM_DUMMY
use grid, only : n_s_grid, nsp,  number_of_species, gsp, nperiod
use geom,       only : ints, sgr, R0_loc, pol_angle, jfun, kfun 
use geom,       only : bn_G, dBdpsi, iterate_fsa, isg_lfs
use components, only : fp, tp, isg_icrh, icrh_fmin, de_G, icrh_model 
use components, only : tmp_G, mas_G, vp_G, tp_G, signz_G, fp_G 
use components, only : icrh_params, icrh_norm, adiabatic_electrons
use components, only : cf_mass_weight, cf_dtf1
use general, only : gkw_abort, gkw_warn
use mpiinterface, only : root_processor, gather_array
use control, only : io_legacy
use global, only : dotdat
  
logical, intent(in) :: write_out

real, dimension(number_of_species), intent(out) :: fsa_dens, RLNFSA
real, dimension(1:n_s_grid), intent(in) :: phi

real, dimension(n_s_grid,nsp) :: rlne_l, rlne1_l, n_l !Diagnostic local dummies
real, dimension(nsp) :: fsa_dens_l

real, dimension(n_s_grid,number_of_species) :: rlne, rlne1, n, dcfen_dpsi  !Diagnostic global dummies
real, dimension(number_of_species) :: e_1, e_2, e_3, e_4, F_V_SP, F_V 
real, dimension(n_s_grid) :: bnn, denom, dphi_est_dpsi

integer :: i,is, file_unit
logical :: ALL_PROCS=.true.
real    :: QN_dum, QN_dum2, eta, deta, e_eta = 0.0, tpp, dtpp, dBR0dpsi = 0.0

 if(rotation_parallelized) then 
    call gkw_abort('centrifugal_diagnostics: call before parallelize')
 endif

 dcfen_dpsi(:,:)=0.E0

 ! These diagnostics are of little interest without the centrifugal terms.
 ! if (.not.cf_trap) return
 ! But it makes post processing easier if they are output anyway

  ! density and Flux surface average density
  do is=1,nsp
    !implicit element by element multiplication of conformable arrays
    fsa_dens_l(is)=sum(exp(-cfen(1:n_s_grid,is))*ints(1:n_s_grid))
    n_l(1:n_s_grid,is)=exp(-cfen(1:n_s_grid,is))
  end do

  ! gather over species
  call gather_array(fsa_dens(1:number_of_species), number_of_species,  &
        & fsa_dens_l(1:nsp),nsp,COMM_SP_NE, ALLGATHER = ALL_PROCS)

  call gather_array(n(1:n_s_grid,1:number_of_species), n_s_grid,number_of_species,  &
        & n_l(1:n_s_grid,1:nsp),n_s_grid,nsp,COMM_DUMMY,COMM_SP_NE, ALLGATHER = ALL_PROCS)

  ! calculate the Angioni memo quantities (also calculated by NEO):
  do is = 1, number_of_species

    do i = 1, n_s_grid
      ! Radial derivative of centrifugal energy 
      dcfen_dpsi(i,is)=(signz_G(is)*dcfphi_dpsi(i)-vcor*vcor*kfun(1,i)*mas_G(is))/tmp_G(1,is)
      !Correction of the 2012 Errata to the 2009 papers
      if (cf_upphi) dcfen_dpsi(i,is)=dcfen_dpsi(i,is)+2.*mas_G(is)*jfun(1,i)*vp_G(1,is)*vcor/tmp_G(1,is)
    end do

    !e_0 = fsa_dens

    e_1(is) = sum(n(1:n_s_grid,is)*signz_G(is)*phi(1:n_s_grid)*ints(1:n_s_grid)/tmp_G(1,is))
    if ((abs(icrh_params(3)) > 999.0)) e_1(is) = 0.0

    e_2(is) = sum(n(1:n_s_grid,is)*signz_G(is)*dcfphi_dpsi(1:n_s_grid)*ints(1:n_s_grid)/tmp_G(1,is))

    e_3(is) = sum(n(1:n_s_grid,is)*jfun(1,1:n_s_grid)*ints(1:n_s_grid))

    e_4(is) = sum(n(1:n_s_grid,is)*kfun(1,1:n_s_grid)*ints(1:n_s_grid))

    F_V_SP(is) = -e_2(is) - e_3(is)*mas_G(is)*2.0*vp_G(1,is)*vcor/tmp_G(1,is)  &
               &          + e_4(is)*mas_G(is)*vcor*vcor/tmp_G(1,is)            &
               &          - e_1(is)*tp_G(1,is)                                 &  
               & + e_3(is)*tp_G(1,is)*mas_G(is)*vcor*vcor/tmp_G(1,is)
               
    F_V_SP(is) = F_V_SP(is)/fsa_dens(is)             

    ! This is the factor needed to go from rln_R0 to rln_FSA
    ! as above, but uses the uprim gradient only from species 1 only
    F_V(is)    = -e_2(is) - e_3(is)*mas_G(is)*2.0*vp_G(1,1)*vcor/tmp_G(1,is)   &
               &          + e_4(is)*mas_G(is)*vcor*vcor/tmp_G(1,is)            &
               &          - e_1(is)*tp_G(1,is)                                 &  
               & + e_3(is)*tp_G(1,is)*mas_G(is)*vcor*vcor/tmp_G(1,is)
               
    F_V(is)    = F_V(is)/fsa_dens(is)
    
  end do ! species

  ! Quantity for transformations of the anisotropic minority density gradient
  if (icrh_model > 0 .and. (abs(icrh_params(3)) < 999.0)) then
    ! FJC Fix dBdpsi axis (R0_loc)
    select case(R0_loc)
      case('LFS'); dBR0dpsi = dBdpsi(1,isg_lfs)
      case('axis'); dBR0dpsi = 0.0
      case default; call gkw_abort('R0_loc error') 
    end select

    bnn(:) = (Bn_G(1,1:n_s_grid)/icrh_norm)

    !corrections for anisotropic minorities
    if (icrh_model == 1) then
      eta = icrh_params(2) - 1.0
      deta = icrh_params(3)
  
      if (isg_icrh /=0) then    
        is = isg_icrh
        do i= 1,n_s_grid
          dcfen_dpsi(i,is) = dcfen_dpsi(i,is) +                             &
                    &        (deta*log(bnn(i)) + eta*dBdpsi(1,i)/bn_G(1,i)  &
                    &                         - eta*dBR0dpsi/icrh_norm)
        end do
      else  ! combination method
        ! not yet implemented
        call gkw_abort('ICRH species combination method removed')
      end if

      !implicit element by element multiplication of conformable arrays
      e_eta  =  sum(n(1:n_s_grid,is)*ints(1:n_s_grid)*                     &
                    &  ( -deta*log(bnn(:))                                 &
                    &    -eta*dBdpsi(1,1:n_s_grid)/Bn_g(1,1:n_s_grid)      & 
                    &    +eta*dBR0dpsi/icrh_norm*bnn(:)/bnn(:) )  )

    else if (icrh_model == 2) then
      tpp = icrh_params(2)
      dtpp = icrh_params(3)
      denom(:) = -(tpp*bnn(:)/bnn(:) + (1.0 - tpp)/bnn(:))

      if (isg_icrh /=0) then    
        is = isg_icrh
        do i= 1,n_s_grid
          dcfen_dpsi(i,is) = dcfen_dpsi(i,is) -                           &
                           &  (dtpp +                                     & 
                           &  (- dBdpsi(1,i)/bnn(i)**2.0/icrh_norm        &
                           &   + dBR0dpsi/bnn(i)/icrh_norm)*(1.0- tpp)    &
                           &   - dtpp/bnn(i))/denom(i)                    
                           
        end do
      else ! combination method
        call gkw_abort('ICRH species combination method removed')
      end if

      !implicit element by element multiplication of conformable arrays
      e_eta  =  sum(n(1:n_s_grid,is)*ints(1:n_s_grid)*                    &
                    &  (   dtpp*bnn(:)/bnn(:)                             &
                    &  + (-dBdpsi(1,1:n_s_grid)/bnn(:)**2.0/icrh_norm     &
                    &     +dBR0dpsi/bnn(:)/icrh_norm)*(1.0- tpp)          &
                    &     -dtpp/bnn(:)   ) / denom(:) )
    end if

    if (isg_icrh /= 0) then   ! independent method
      is = isg_icrh
      F_V(is)    = F_V(is)    + e_eta /fsa_dens(is)              
      F_V_SP(is) = F_V_SP(is) + e_eta /fsa_dens(is) 
    else                      ! combination method
      call gkw_abort('ICRH species combination method removed')
    end if   
      
  end if

  ! Density and effective density gradients
  ! Calculate as quantities which are local in species
  do is=1,nsp   
    do i=1,n_s_grid

      !And analytic estimate for singly charged ions (returns zero if none)
      dphi_est_dpsi(i)=cf_dtf1()*vcor*vcor*jfun(1,i)+cf_mass_weight()*vcor*vcor*kfun(1,i)
      !Correction of the 2012 Errata to the 2009 papers
      if (cf_upphi) dphi_est_dpsi(i)=dphi_est_dpsi(i)-2.*cf_mass_weight()*vcor*vp_G(1,1)*jfun(1,i)

      ! Use cfeno because the additional part of the ICRH energy does not have on 1/T
      ! clearer DOCS needed.  Without this, QN tests fail. 
      rlne_l(i,is)=(dcfen_dpsi(i,gsp(is))+tp(1,is)*cfeno(i,is))+fp(1,is)
      rlne1_l(i,is)=rlne_l(i,is)*exp(-cfen(i,is))/fsa_dens_l(is)
    end do
  end do ! species

  ! And gather over species
  call gather_array(rlne(1:n_s_grid,1:number_of_species), n_s_grid,number_of_species,  &
        & rlne_l(1:n_s_grid,1:nsp),n_s_grid,nsp,COMM_DUMMY,COMM_SP_NE, ALLGATHER = ALL_PROCS)

  call gather_array(rlne1(1:n_s_grid,1:number_of_species), n_s_grid,number_of_species,  &
        & rlne1_l(1:n_s_grid,1:nsp),n_s_grid,nsp,COMM_DUMMY,COMM_SP_NE, ALLGATHER = ALL_PROCS)

  ! set outputs to return
  do is=1,number_of_species
    RLNFSA(is) = sum(rlne1(1:n_s_grid,is)*ints(1:n_s_grid))
  end do

  !Write to file and to screen
  if (root_processor .and. write_out) then
    
    call output_array(dotdat('cfdens',io_legacy), 'rotation', reshape( &
       & (/ ( sgr(i), rlne(i,:), n(i,:), pol_angle(1,i) , i=1,n_s_grid) /), &
       & (/ 1+size(rlne,2)+size(n,2)+1, n_s_grid /)), &
       & 'F', xy_fmt_long, ascii_fmt)

    call output_array(dotdat('cfphi',io_legacy), 'rotation', reshape( &
       & (/ ( phi(i), dcfphi_dpsi(i), dphi_est_dpsi(i), dcfphi_ds(i) ,&
       & i=1,n_s_grid) /), &
       & (/ 4, n_s_grid /)), &
       & 'F', xy_fmt, ascii_fmt)
    
    if (.not.cf_trap) return !don't write to screen
  
    if (iterate_fsa) then
      write(*,*) 'Input densities and density gradients defined as FSA values'
      write(*,*) '    reference location used internally is R0=R_'//trim(R0_loc)
    else
      write(*,*) 'Input densities and density gradients defined at R=R0=R_'//trim(R0_loc)
    endif
    if (mod(n_s_grid,2)==0) write(*,*) 'Outboard MP value below is for nearest grid point'
    !Need an odd number of points for output "outboard MP" value to be exactly at LFS
    write(*,*)
    
    write(*,*) 'Densities for kinetic species (relative to n_R0):'
    do is=1,number_of_species !global species number
            write(*,'((A,i2),2(A,f11.6))') ' Species: ', is, ' :  n / n_R0  - &
            &  FSA: ',  fsa_dens(is),                  &
            &  '    LFS: ',  n(isg_lfs,is)
            !Note outboard MP is NOT necessarily the same as Rmax in general geom.
    end do
    write(*,*)

    !if (iterate_fsa) then
      write(*,*) 'Densities for kinetic species (relative to n_ref):'
      do is=1,number_of_species !global species number
              write(*,'((A,i2),2(A,f11.6))') ' Species: ', is, ' :  n / n_ref - &
              &  FSA: ',  fsa_dens(is)*de_G(1,is),                         &
              &  '    LFS: ',  n(isg_lfs,is)*de_G(1,is)
      end do
    !end do
    
    ! Report the exact ICRH minority fraction if applicable
    if (isg_icrh /=0 .and. number_of_species >= 2) then
      write(*,*)
      write(*,'(A,f11.6)') ' ICRH minority fraction (exact FSA): ',  de_G(1,isg_icrh)*fsa_dens(isg_icrh)/de_G(1,2)/fsa_dens(2)
      write(*,*)  'Assumption: Electrons are species 2'
    else if (icrh_fmin > 0.0 .and. number_of_species >= 2) then  
      write(*,*)
      write(*,'(A,f11.6)') ' ICRH minority fraction (exact FSA): ',  &
         & icrh_fmin*sum(exp(-cfen_icrh(1:n_s_grid))*ints(1:n_s_grid))/de_G(1,2)/fsa_dens(2)
      write(*,*) 'Assumptions: Bulk ions + minority are species 1, electrons are species 2'
    end if    
  
    !And effective rln
    !write(*,*)
    !write(*,*) 'Effective local rln for kinetic species:'
    !do is=1,number_of_species !global species number
    !        write(*,'((A,i2),2(A,f11.6))') ' Species: ', is, ' :   rln - &
    !        & Outboard MP: ',  rlne(isg_lfs,is),                         &
    !        & ' :  Flux surface average: ', sum(rlne(1:n_s_grid,is)*ints(1:n_s_grid))
    !        !Note implict element by element multiplication of conformable arrays
    !        !Need an odd number of points for LFS value to be exactly at LFS
    !        !Note outboard MP is NOT necessarily the same as Rmax in general geom.
    !end do
  
    !And effective rln_1
    !write(*,*)
    !write(*,*) 'Effective rln_1 for kinetic species:'
    !do is=1,number_of_species !global species number
    !        write(*,'((A,i2),2(A,f11.6))') ' Species: ', is, ' :  rln1 - &
    !        & Outboard MP: ',   rlne1(isg_lfs,is),                       &
    !        & ' :  Flux surface average: ', sum(rlne1(1:n_s_grid,is)*ints(1:n_s_grid))
    !        !Note implict element by element multiplication of conformable arrays
    !        !Note outboard MP is NOT necessarily the same as Rmax in general geom.
    !end do
    !write(*,*)

    write(*,*)
    write(*,*) 'In-Out Asymmetry: '
    do is=1,number_of_species 
      write(*,'((A,i2),2(A,f11.6))') ' Species: ', is, ' : (n_LFS - n_HFS) / n_FSA', & 
                             &  (n(isg_lfs,is)-n(1,is))/fsa_dens(is)
    end do
    
    !And rln of the flux surface average
    write(*,*)
    write(*,*) 'RLN of the FSA density for kinetic species: (memo transformations)'
    write(*,*) 'Uses uprim of species 1, as used by dPhi/dpsi solve'
    do is=1,number_of_species 
      write(*,'((A,i2),2(A,f11.6))') ' Species: ', is, ' : RLN_FSA : ', & 
                & fp_G(1,is) - F_V(is), ' :      RLN_R0 : ', fp_G(1,is)    
    end do
    write(*,*)

    write(*,*) 'RLN of the FSA density for kinetic species (code diagnostics):'
    write(*,*) 'Uses uprim of each species, as used by RHS'
    do is=1,number_of_species 
      write(*,'((A,i2),2(A,f11.6))') ' Species: ', is, ' : RLN_FSA : ', & 
                & RLNFSA(is), ' :      RLN_R0 : ', fp_G(1,is)    
    end do
    write(*,*)

    if(output_enabled) then
      ! output the Angioni memo quantities (Vcal_spgr is output also by NEO)
      call get_free_file_unit(file_unit)
      open(file_unit,FILE = 'cffsa.dat')

      write(file_unit,*) 'e_0' 
      write(file_unit,fmt = '(256(es13.5,1x))')  fsa_dens(:)    

      write(file_unit,*) 'e_1'   
      write(file_unit,fmt = '(256(es13.5,1x))') e_1(:)

      write(file_unit,*) 'e_2'   
      write(file_unit,fmt = '(256(es13.5,1x))') e_2(:)

      write(file_unit,*) 'e_3'   
      write(file_unit,fmt = '(256(es13.5,1x))') e_3(:)

      write(file_unit,*) 'e_4'   
      write(file_unit,fmt = '(256(es13.5,1x))') e_4(:)

      e_4(:) = 0.0 
      if (isg_icrh /= 0) e_4(isg_icrh) = e_eta
      write(file_unit,*) 'e_eta'   
      write(file_unit,fmt = '(256(es13.5,1x))') e_4(:)

      write(file_unit,*) 'Vcal_spgr'    
      write(file_unit,fmt = '(256(es13.5,1x))') F_V_SP(:)

      write(file_unit,*) 'Vcal_iongr_uprim_only'    
      write(file_unit,fmt = '(256(es13.5,1x))') F_V(:)

      write(file_unit,*) 'RLN_FSA'
      write(file_unit,fmt = '(256(es13.5,1x))') fp_G(1,1:number_of_species) - F_V(:)

      if (iterate_fsa) then

        write(file_unit,*) 'RLN_FSA_code'
        write(file_unit,fmt = '(256(es13.5,1x))') RLNFSA(:)

        write(file_unit,*) 'n_R0'
        write(file_unit,fmt = '(256(es13.5,1x))') de_G(1,:)

        write(file_unit,*) 'RLN_R0'
        write(file_unit,fmt = '(256(es13.5,1x))') fp_G(1,:)

      end if

      close(file_unit)
    end if

    ! Check the quasi-neutrality of the FSA gradients 
    ! This is a very strong consistency check on the QN solver and diagnostics
    if (.not. adiabatic_electrons) then
      QN_dum = sum( signz_G(1:number_of_species)*de_G(1,1:number_of_species)*  &
                &   fsa_dens(:)*(fp_G(1,1:number_of_species) - F_V(:)) )

      QN_dum2 = sum( signz_G(1:number_of_species)*de_G(1,1:number_of_species)*  &
                &   fsa_dens(:)*(RLNFSA(:)) )

      !write(*,*) ' Quasi-neutral FSA gradients ?: ',  QN_dum, QN_dum2, F_V
      if (abs(QN_dum) > 1e-9 .or. abs(QN_dum2) > 1e-9) then
        write(*,*) ' Quasi-neutral FSA gradients (memo, code) ?: ',  QN_dum, QN_dum2 
        call gkw_warn('NO QUASI-NEUTRALITY IN THE FSA GRADIENTS')  
        call gkw_warn('With CF asymmetry, all non-trace species must have same uprim')
        call gkw_warn('To override, set cf_qncheck = .false. in ROTATION')        
        if (cf_qncheck) call gkw_abort('DO NOT USE THIS RUN FOR PARTICLE TRANSPORT')      
        call gkw_warn('DO NOT USE THIS RUN FOR PARTICLE TRANSPORT')
      end if
    end if  
  
  end if ! root processor

end subroutine centrifugal_diagnostics


!--------------------------------------------------------------------
!> This subroutine initialises anything required for shearing
!> or centrifugal force
!> Depending on the shearing method, different init routines are called
!> This routine also recaluates some init quantities from non_linear terms
!> This is done to avoid a circular dependency
!--------------------------------------------------------------------
subroutine rotation_init

  use grid,         only : nx, number_of_species,n_s_grid
  use mode,         only : ixzero
  use geom,         only : efun, signB, iterate_fsa
  use components,   only : vp_G, de_G, fp_G, de, fp
  use grid,         only : nsp, lsp
  use control,      only : spectral_radius
  use general,      only : gkw_abort, gkw_warn
  use mpiinterface, only : root_processor
  use constants,    only : pi 

  ! all global in species quantities
  real :: defsa_nR0(number_of_species), fp_fsa(number_of_species)
  real :: de_in(number_of_species), fp_in(number_of_species)

  real :: accuracy

  integer ix, il, is
  integer, parameter :: maxit = 50

  ! set the radial surface and do all calculations for one surface (at 
  ! the moment 
  ix = 1 

    !Solve for the centrifugal energy
    !Must call to intialise energy even if not keeping centrifugal terms
    !since cfen is used in field equations.
    call centrifugal_energy
    if (cf_trap) call centrifugal_diagnostics(cfphi(1:n_s_grid,0),defsa_nR0,fp_fsa,(.not.iterate_fsa))

    ! In the case the user wants to input FSA vlues for dens and rln,
    ! chosen by the geom input switch R0_loc='FSA',
    ! the asymmetry solver uses R0_loc='LFS' and we iterate the LFS values 
    ! to converge on desired FSA values
    ! MOVE to a subroutine ?
    iterate_fsa_cond: if (iterate_fsa .and. cf_trap) then

      !Store input values (desired FSA values)
      de_in=de_G(1,:)
      fp_in=fp_G(1,:)

      accuracy=10.E0**(-precision(1.E0))

      ! iterate on density
      iterate: do il=1,maxit

        !write(*,*) 'de_fsa     ',  de_G(1,:)*defsa_nR0

        ! intuitive (overshoots and corrects fast)
        de_G(1,:) = de_in(:)/defsa_nR0(:)
        ! ugly hack to keep all global variables consitent
        call distribute_species_arrays

        !write(*,*) 'de_G   ', de_G(1,:)

        call centrifugal_energy ! uses de_G values inside components
        call centrifugal_diagnostics(cfphi(1:n_s_grid,0),defsa_nR0,fp_fsa,.false.)

        ! condition to stop iterating
        if (sum(abs(de_G(1,:)*defsa_nR0(:)-de_in(:))) < accuracy) exit iterate

      end do iterate
      
      if (root_processor) write(*,*) 'n_fsa loop', il, sum(abs(de_G(1,:)*defsa_nR0(:)-de_in(:)))
      if (root_processor) write(*,*)

      !if (il >= maxit) then 
      if (sum(abs(de_G(1,:)*defsa_nR0(:)-de_in(:))) > accuracy) then 
        call gkw_warn('electrons and ions require same uprim for R0_loc=FSA')
        call gkw_abort('R0_loc=FSA density iteration failed to converge')
      end if

      ! iterate on gradients
      iterate2: do il=1,maxit

        !write(*,*) 'fp_G   ', fp_G(1,:)
        !write(*,*) 'fp_in     ',  fp_in
        !write(*,*) 'fp_fsa     ',  fp_fsa
       
        ! intuitive, but prone to oscillatory solutions
        fp_G(1,:) = fp_in(:) - fp_fsa(:) + fp_G(1,:)

        !small quasineutral kicks to prevent oscillating solutions (does not work)
        !fp_G(1,:) = fp_G(1,:) + 0.8*rosc*ran_kiss()*signz_G(:)*de_G(1,:)*(fp_fsa(:)-fp_in(:))
        ! grows slowly towards 1      
        !rosc = -abs(rosc)**0.75

        ! Most cases converge in about 5 iterations.
        ! A few cases get stuck oscillating around 1e-7 convergence level.
        ! Kludge: add non quasineutral grit to jump out of stable oscillating solutions
        ! In all cases tested this allows convergence towards machine precision, 1e-15
        if (il == 10) fp_G(1,1) = fp_G(1,1) - 1e1*(fp_in(1) - fp_fsa(1))
        if (il == 15) fp_G(1,1) = fp_G(1,1) + 1e2*(fp_in(1) - fp_fsa(1))
        if (il == 22) fp_G(1,1) = fp_G(1,1) - 2e2*(fp_in(1) - fp_fsa(1))
        if (il == 30) fp_G(1,1) = fp_G(1,1) + 3e2*(fp_in(1) - fp_fsa(1))
        if (il == 40) fp_G(1,1) = fp_G(1,1) - 1e3*(fp_in(1) - fp_fsa(1))

        ! ugly hack to keep all global variables consistent
        call distribute_species_arrays

        !write(*,*) 'fp_G (new)', fp_G(1,:)

        call centrifugal_energy ! uses de_G values inside components
        call centrifugal_diagnostics(cfphi(1:n_s_grid,0),defsa_nR0,fp_fsa,.false.)

        ! condition to stop iterating (not so accurate due to numerical derivatives)
        ! write(*,*) il, sum(abs(fp_in(:) - fp_fsa(:))), (fp_in(1:2) - fp_fsa(1:2))

        if (sum(abs(fp_in(:) - fp_fsa(:))) < 1e1*accuracy) exit iterate2

      end do iterate2

      if (root_processor) write(*,*) 'RLN_fsa loop', il, sum(abs(fp_in(:) - fp_fsa(:)))
      if (root_processor) write(*,*)

      if (il >= maxit .and. sum(abs(fp_in(:) - fp_fsa(:))) > 1e-7) then 
        if (root_processor) write(*,*) sum(abs(fp_in(:) - fp_fsa(:))), 1e-7
        call gkw_abort('R0_loc=FSA gradient iteration failed to converge')
      end if 

      !write(*,*) de_G(1,:)*defsa_nR0
      call centrifugal_diagnostics(cfphi(1:n_s_grid,0),defsa_nR0,fp_fsa,.true.)

    end if iterate_fsa_cond

    if(.not.perp_shear) return

    !Override input shear rate for toroidal shear
    !Ideally these values would be written to input.out
    !Due to the order of init this is impractical 
    !Since all namelists are read / written before geom_init
    select case(toroidal_shear)

      case ('use_uprim') ! Use uprim of first species, ignore input
        !Valid for general geometry
        shear_rate = signB*vp_G(1,1)/(4.E0*pi*efun(ix,1,1,2))
        if (root_processor) then
        !In future this info should go to input.out
            write(*,*)
            write(*,*) 'shear_rate used : ', shear_rate
            write(*,*)
        end if
  
      case ('add_uprim') ! Input ExB shear is from poloidal flow
        !Valid for general geometry
        shear_rate = shear_rate + signB*vp_G(1,1)/(4.E0*pi*efun(ix,1,1,2))
        if (root_processor) then
        !In future this info should go to input.out
            write(*,*)
            write(*,*) 'shear_rate used : ', shear_rate
            write(*,*)
        end if
  
      case ('use_shear_rate') ! Overrides all input values of uprim
        !Valid for general geometry
        ts_uprim=4.E0*pi*signB*efun(ix,1,1,2)*shear_rate
        !In future this info should go to input.out
        if (root_processor) then
            write(*,*)
            write(*,*) 'uprim used : ', ts_uprim
            write(*,*)
        end if

      case ('none') ! Each is controlled independently
        ! Do nothing

      case default
       call gkw_abort('uknown option for toroidal_shear') 

    end select

    !Now account for the different coordinate normalisations of kx and ky
    !The factor 2 arises because v_s & v_\chi have different normalisations
    !Due to definiton of thermal velocity
    !Note this routine is now called after parallelize geom
    !Rescaling works since efun(i,1,2) is a flux function anyway.
    !because of this rescaling, shear_rate should not be public
    shear_rate=2.*efun(ix,1,1,2)*shear_rate

    !Equivalent to: 
    !*(lx/ly box picture) shear_rate=s_j*shear_rate*q/(2.*pi*eps)
    !*(normalisation picture) shear_rate=s_j*shear_rate*kthnorm
    !This is normalisation FORWARDS from rho ref units to angular units

    if(shear_real) call shear_real_init

    if(shear_remap) call shear_remap_init

    if(shear_ky_shift) then
    !With odd number of modes
    !The nyquist frequency is repeated - incorrectly?
    !Setup index arrays for ky_shift
      if (spectral_radius) then
        do ix = 1, nx
          jind_nx(ix) = ix
        end do
      else
        do ix = ixzero, nx
          jind_nx(ix) = ix - ixzero + 1
        end do 
 
        do ix = ixzero-1, 1, -1 
          jind_nx(ix) = nx + ix - ixzero + 2
        end do
      end if
    endif

contains

  !--------------------------------------------------------------------
  !> broadcast / distribute everything back across ix and local species arrays
  !--------------------------------------------------------------------
  subroutine distribute_species_arrays

      do is= 1,number_of_species
        de_G(:,is) = de_G(1,is)
        fp_G(:,is) = fp_G(1,is)

        ! if species on local processor, copy back local values
        if (lsp(is) >= 1 .and.  lsp(is) <= nsp) then
          de(:,lsp(is)) = de_G(:,is)
          fp(:,lsp(is)) = fp_G(:,is)
        end if 
      end do

  end subroutine distribute_species_arrays

end subroutine rotation_init


!--------------------------------------------------------------------
!> This subroutine initialises the grad_pot array for the ExB shearing
!> When added directly in add_non_linear terms, with boundary discontinuity
!--------------------------------------------------------------------
subroutine shear_real_init

use grid,         only : ns, lx 
use control,      only : spectral_radius
use general,      only : gkw_abort

integer j, ipar, qmrad

if (.not. spectral_radius) return

if(.not.shear_real) call gkw_abort('Invalid call to shear_real_init')

!Beware integer division
qmrad = mrad/4

!!NEED TO CHECK NORMALISATION OF SHEAR RATE IF USING THIS METHOD

  if(4*qmrad.ne.mrad) then
      write(*,*) 'Warning: mrad not a multiple of 4' 
      write(*,*)  '-> small discontinuity in shear'
  endif
  !Could do this with qmrad=real(mrad/4).

  !if(root_and_verbose) write(*,*) 'mrad', mrad, 'qmrad', qmrad, 'lx', lx 

  select case(shear_profile)

  case('linear') !Linear shear discontinuous in periodic radial boundary
            if(root_and_verbose) write(*,*) 'Discontinuous linear shearing' 

            shear_rate=shear_rate*lx/mrad
 
            do ipar = 1, ns  
                do j = 1, mrad
                    !grad_pot(j, ipar)= (-(mrad/2)+j)*bn(ipar)*shear_rate
                    grad_pot(j, ipar) =(-(mrad/2)+j)*shear_rate
                end do
            end do

  case('symmetric') !Symmetric profile continuous at boundary
            !The factor of B_N in geom has been removed.
            !Factor of two to cancel out efun divisor of two???
            shear_rate=shear_rate*lx/mrad

            if(root_and_verbose) write(*,*) 'Symmetric shearing profile' 
            do ipar = 1, ns
                !mrad will always be a multiple of 2
                do j = 1, mrad/2
                  !grad_pot(j, ipar) =(-qmrad+j)*bn(ipar)*shear_rate
                  grad_pot(j, ipar) =(-qmrad+j)*shear_rate
                end do

                do j = mrad/2+1, mrad
                  !grad_pot(j, ipar) =(3*qmrad-j)*bn(ipar)*shear_rate
                  grad_pot(j, ipar) =(3*qmrad-j)*shear_rate 
                end do
            end do
  
  case default
        call gkw_abort ('shear_real_init called incorrectly')
  
  end select

end subroutine shear_real_init

!--------------------------------------------------------------------
!> This subroutine allocates the arrays needed for the rotation module
!--------------------------------------------------------------------
subroutine rotation_allocate

  use grid, only : n_s_grid,nx,nmod,nsp,nmu,nvpar, ns
  use components, only : iadia, nsps
  use general, only : gkw_abort

  integer ierr
  real dum

  !Note this must be identical to the one in non_linear_terms
  !Must also be calculated here to avoid a circular dependecy
  dum  = 1.5*real(nx+1) 
  mrad = int(log(dum)/log(2.E0) + 1.E0)
  mrad = 2**mrad 
  mx=2*nx

  !Always allocate the following since it may be read from restart file.
  !But only used with shear_remap
  ierr=0
  allocate(kxshift(nmod),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate kxshift')
  !initialise
  kxshift=0

  !Allocate the array for the centrifugal energy
  ierr=0
  ! 2 extra points, parallel or otherwise
  allocate(cfen(1-2:n_s_grid+2,nsp+iadia),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfen')
  !initialise to be safe
  cfen(:,:)=0.

  allocate(cfen_hfs(nsp+iadia),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfen_hfs')
  allocate(cfen_lfs(nsp+iadia),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfen_lfs')

  !Allocate the workspace array for cf_phi
  ierr=0
  allocate(cfphi(n_s_grid,-2:2),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate (CF) phi')
  !initialise to be safe
  cfphi(:,:) = 0.0

  ierr=0
  allocate(cfeno(1:n_s_grid,nsp+iadia),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfeno')
  !initialise to be safe
  cfeno(:,:)=0.

  !Allocate the array for the centrifugal energy high field side values
  ierr=0
  allocate(cfenh(nsp+iadia),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfenh')
  !initialise to be safe
  cfenh(:)=0.

  !Allocate the array for the centrifugal energy low field side values
  ierr=0
  allocate(cfenl(nsp+iadia),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfenl')
  !initialise to be safe
  cfenl(:)=0.

  !Allocate the array for the icrh energy
  allocate(cfen_icrh(n_s_grid),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfen_icrh')
  !initialise to be safe
  cfen_icrh(:)=0.

  !Allocate the array for the centrifugal energy, global in species
  ierr=0
  allocate(cfen_G(ns,nsps),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate cfen_G')
  !initialise to be safe
  cfen_G(:,:)=0.

  !Derivative arrays do no need to keep the adiabatic species data
  !s derivative of cfen
  ierr=0
  allocate(dcfen_ds(n_s_grid,nsp),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate dcfen_ds')
  !initialise to be safe
  dcfen_ds(:,:)=0.

  !Radial derivative of cfphi
  ierr=0
  allocate(dcfphi_dpsi(n_s_grid),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate dcfphi_dpsi')
  !initialise to be safe
  dcfphi_dpsi(:)=0.

  !Parallel derivative of cfphi
  ierr=0
  allocate(dcfphi_ds(n_s_grid),stat=ierr)
  if (ierr.ne.0) call gkw_abort('Could not allocate dcfphi_ds')
  !initialise to be safe
  dcfphi_ds(:)=0.

  if(.not.perp_shear) return

  if(shear_real) then
    ierr = 0
    allocate(grad_pot(mrad,ns), stat = ierr)
    if (ierr.ne.0) then
      write(*,*)'Could not allocate grad_pot in module rotation'
      stop 1
    end if
  end if

  if(shear_ky_shift) then
    ierr = 0
    allocate(arr(nx+1,nmod), stat = ierr)
    if (ierr.ne.0) then
      write(*,*)'Could not allocate arr in module rotation'
      stop 1
    end if

    ierr = 0 
    allocate(jind_nx(nx), stat = ierr) 
    if (ierr.ne.0) then 
      write(*,*)'Could not allocate jind_nx in module rotation'
      stop 1
    end if
  end if

  if (shear_remap) then
    ierr=0
    allocate(aindx(nmod*nmu*nvpar*ns*nsp*(nx-1)),stat=ierr)
    if (ierr.ne.0) call gkw_abort('Could not allocate aindx in rotation')

    ierr=0
    allocate(aindx_shift(nmod*nmu*nvpar*ns*nsp*(nx-1)),stat=ierr)
    if (ierr.ne.0) call gkw_abort('Could not allocate aindx_shift in rotation')

  end if

end subroutine rotation_allocate

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> read and broadcast or write the kxshift restart namelist
!> For restart consistency, kxshift must be stored.
!> although there may be a cleaner way to rewrite without this
!----------------------------------------------------------------------------
subroutine kxshift_read_nml(ilun,io_stat,lwrite)
  use grid, only : nmod
  use general, only : gkw_abort, gkw_warn
  use global, only : r_tiny

  integer, intent(in) :: ilun
  integer,intent(out) :: io_stat
  logical, optional, intent(in) :: lwrite
  integer, parameter :: nmmx = 512
  integer, dimension(nmmx) :: kxshft
  integer :: imod
  real :: t_dum

  !kxshift is always the same on all processors
  namelist /kx_shift/ t_shear_begin, kxshft
  if (nmod.gt.nmmx) call GKW_ABORT('Reset kxshft nmmx and recompile')
  kxshft(:)=0

  if (present(lwrite)) then !read
    if ((.not. lwrite)) then
      ! read nml
      t_dum = t_shear_begin

      read(ilun,NML=kx_shift,IOSTAT=io_stat)
      do imod=1,nmod
         kxshift(imod)=kxshft(imod)
      end do

      ! case where t_shear_begin is changed in input.dat at restart
      if (t_dum .gt. t_shear_begin + r_tiny) then
          t_shear_begin = t_dum
          call gkw_warn('t_shear_begin used from input.dat, not restart file')
      end if 
     
      ! case where restart overwrites t_shear_begin (e.g. not set in input.dat)
      if (abs(t_dum - t_shear_begin) .gt. r_tiny) then
          call gkw_warn('t_shear_begin used from restart file, not input.dat')
      end if 
    else
      ! do nothing
    end if
  else !write
      do imod=1,nmod
        kxshft(imod)=kxshift(imod)
      end do
      write(ilun,NML=kx_shift)
  end if

end subroutine kxshift_read_nml


!---------------------------------------------------------------------
!>
!---------------------------------------------------------------------
subroutine kxshift_bcast_nml
  
  use mpiinterface, only : mpibcast
  use grid,         only : nmod

  call mpibcast(kxshift,nmod)
  call mpibcast(t_shear_begin,1)

end subroutine kxshift_bcast_nml

!---------------------------------------------------------------------
!> Takes the fourth order parallel first derivative
!> operates on real global arrays of size n_s_grid
!> uses periodicity at end of n_s_grid
!> must be called before parallelize_geom
!---------------------------------------------------------------------
subroutine parallel_derivative(in,out)

use grid,    only : n_s_grid
use geom,    only : sgr_dist, geom_parallelized
use general, only : gkw_abort

real, intent(in), dimension(n_s_grid) :: in
real, intent(out), dimension(n_s_grid) :: out
real, dimension(-1:n_s_grid+2) :: dum

integer :: is, ic
real :: coeff4d1(-2:2)

 if(geom_parallelized.or.rotation_parallelized) then
     call gkw_abort('call parallel_derivative before parallelize_rotation')
 end if

 !intialise
 out(:)=0.E0

 !copy in to dum
 dum(1:n_s_grid)=in(1:n_s_grid)
 
 !wrap ghost cells
 dum(0)=in(n_s_grid)
 dum(-1)=in(n_s_grid-1)
 dum(n_s_grid+1)=in(1)
 dum(n_s_grid+2)=in(2)

 !Setup finite differencing coeffcients.  4th order 1st derivative
 coeff4d1=(/1.E0,-8.E0,0.E0,8.E0,-1.E0/)/(12.E0*sgr_dist)

 !write(*,*) coeff4d1

 do is=1,n_s_grid
   !Finite difference
   do ic=-2,2
       out(is)=out(is)+dum(is+ic)*coeff4d1(ic)
   end do
 end do

end subroutine parallel_derivative

!---------------------------------------------------------------------
!> In order to solve the local problem in s, copy the local elements to
!> the begining of each array used elsewhere in the code. 
!> See also parallelize geom.
!>
!> Before this routine, cfen is global in n_s_grid
!> and cfen_G is NOT initialisied.
!> After this routine both cfen_G and cfen are local in ns
!> And cfen_G is global in species.
!> i.e:
!> cfen(n_s_grid,nsp+iadia)->cfen(ns,nsp+iadia)
!> Then cfen()->cfen_G(ns,nsp)
!>
!> It is a little unfortunate that there are two cfen arrays global
!> in opposite directions; this is due to the initialisation order.
!---------------------------------------------------------------------
subroutine parallelize_rotation

use grid,         only : ns, gs, n_s_grid, nsp,  number_of_species
use mpicomms,     only : COMM_SP_NE, COMM_DUMMY
use components,   only : iadia, nsps, adiabatic_electrons
use general,      only : gkw_abort
use mpiinterface, only : gather_array

integer :: i1, i2
real   :: G(1:2,nsp+iadia) !dummy
logical :: ALL_PROCS

!i2=global s point
!i1=local s point

if (rotation_parallelized) then
    call gkw_abort('parallelize rotation may only be called once')
end if

!Initialize G with silly values
G(:,:)= -3.32176E+21

!cfen is periodic in n_s_grid
!two ghost points of cfen (only) used by arakawa scheme
G(1:2,:)=cfen(1:2,:)

  !left local s boundary
  do i1=-1,0
    i2 = gs(i1)
    if (i2 .lt. 1) then !outside n_s_grid
      cfen(i1,:)=cfen(i2+n_s_grid,:)
    else                !adjacent s processor
      cfen(i1,:)=cfen(i2,:)
    endif
  enddo

  ! all cases can use the following
  do i1=1, ns
    i2 = gs(i1)
    cfen(i1,:)=cfen(i2,:)
    dcfphi_dpsi(i1)=dcfphi_dpsi(i2)
    dcfphi_ds(i1)=dcfphi_ds(i2)
    dcfen_ds(i1,:)=dcfen_ds(i2,:)
  enddo

  !right local s boundary
  do i1=ns+1,ns+2
    i2 = gs(i1)
    if (i2 .gt. n_s_grid) then !outside n_s_grid, wrap
      !cfen(i1,:)=cfen(i2-n_s_grid,:)
      cfen(i1,:)=G(i2-n_s_grid,:)
    else                       !adjacent s processor
      cfen(i1,:)=cfen(i2,:)
    endif
  enddo

  !Gather local cfen(ns,nsp) across species into cfen_G(ns,number_of_species)
  !species ordering after the gather should be as input and tmp_G, de_G, etc.
  !cfen_G is only used for collisions
  ALL_PROCS=.true.
  call gather_array(cfen_G(1:ns,1:number_of_species), ns,number_of_species,  &
      & cfen(1:ns,1:nsp),ns,nsp,COMM_DUMMY,COMM_SP_NE, ALLGATHER = ALL_PROCS)

  !Put in the adiabatic species part of cfen_G, all procs have same so no MPI.
  if (adiabatic_electrons) then
     cfen_G(1:ns,nsps)=cfen(1:ns,nsp+iadia)
  endif

rotation_parallelized=.true.

end subroutine parallelize_rotation

end module rotation
