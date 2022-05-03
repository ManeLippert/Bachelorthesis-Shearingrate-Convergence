!------------------------------------------------------------------------------
!> Module contains most of the information on the modes kept in the
!> simulation. Sets up 2D mode grid and connections.
!------------------------------------------------------------------------------
module mode

  use global, only : lenswitch
  
  implicit none

  private

  public :: mode_read_nml, mode_write_nml, mode_bcast_nml, mode_init
  public :: mode_check_params, krbal, kgrid, mode_box_recon, mode_allocate
  public :: kgrid_nonspec_rad, parallel_phase_shift 

  !> The array of 'toroidal' (zeta) wave vectors, krho(nmod).
  !> These wavevectors are perpendicular to the field line in the plane of
  !> the flux surface. Confusion can arise because they are called both
  !> 'toroidal' and 'poloidal' at various points! This direction is also often
  !> referred to as the y direction.
  real, allocatable, public, save :: krho(:)

  !> The perpendicular wave vector (as a function along the field line)
  real, allocatable, public, save :: krloc(:,:,:)

  !> The array of radial (psi) wave vectors, kxrh(nx).
  real, allocatable, public, save :: kxrh(:)

  !> Integers that determine to which kx mode the mode is connected to
  !> through the parallel boundary conditions. ixplus(nmod,nx) is for the
  !> boundary connection in the positive s-direction?
  integer, allocatable, public, save :: ixplus(:,:)
  integer, allocatable, public, save :: ixminus(:,:) !< Opposite boundary to ixplus.

  !> Logical that determines if there is a 2D grid of ky,kx. Used in
  !> conjunction with nperiod = 1 (necessary for non linear runs).
  logical, public, save :: mode_box

  !> For mode_box, the integer that determines the spacing between the
  !> different kx modes.
  integer, public, save :: ikxspace
  
  !> LX/LY for zero shear case, spectral only
  real, save :: rkxspace

  !> The location of the kx=0 mode.
  integer, public, save :: ixzero

  !> The location of the ky=0 mode.  Needed by linear terms.
  integer, public, save :: iyzero

  !> The spacing of the toroidal mode numbers. Only used in global simulations 
  !> when set to a nonzero value (otherwise the code uses krhomax
  integer, public, save :: n_spacing 

  !> The size of the perpendicular box in real (non-spectral) space.
  real, public, save :: lxinv, lyinv
  real, save :: ly
  
  !> The real space box sizes in units of rhoref = ly/kthnorm
  real, public, save :: lyn, lxn

  !> maximum kx value of the kgrid
  real, public, save :: kxmax, kymax !< maximum ky value of the kgrid

  !> minimum kx value calculated by mode_box_recon
  real, public, save :: kxmin

  !> Spacing between kxmodes
  real, public, save :: kxspace

  !> The poloidal shift of the ballooning transform (mode_box=F)
  real, save :: chin

  !> The radial wave vector (mode_box=F)
  real, save :: krrho

  !> Switch for the radial wave vector / poloidal shift  input 
  character (len=lenswitch), save :: kr_type

  !> Logical for treating the special case when shat=0.
  logical, public, save :: lshat_zero = .false.

  !> For mode_box, the maximum ky used in the simulations (input).
  real, save :: krhomax

  !> Number of distinct modes in a linear run (i.e. not connected through the
  !> parallel boundary conditions).
  integer, public, save :: nmodes_l = 0
  integer, public, save :: nmodes_G = 0

  !> Map to one of the distinct modes in a linear run from each radial and
  !> toroidal point. Each value in the array is an integer between 1 and
  !> nmodes. This can be used for normalisations, calculating the individual
  !> mode frequencies and growth rates etc.
  integer, allocatable, public, save :: mode_label_l(:,:), mode_label_G(:,:)

  !> Size of kthrho array.
  integer, parameter :: nmmx = 512

  !> Input values for k_theta rho when mode_box=F
  real, save :: kthrho(nmmx)

  !> list of wavenumbers (indices, i.e. the zero-mode corresponds to 1),
  !> from/to which nonlinear transfer shall be suppressed.
  integer, save :: no_transfer_to(nmmx)
  integer, save :: no_transfer_from(nmmx)
  integer, allocatable, public, save :: no_transfer_to_modes(:,:)
  integer, allocatable, public, save :: no_transfer_from_modes(:,:)
  logical, public, save :: erase_any_transfer_to, erase_any_transfer_from

  !> list of wavenumbers (indices, i.e. the zero-mode corresponds to 1),
  !> from/to which linear drive shall be suppressed.
  integer, save :: no_drive_of(nmmx)
  integer, allocatable, public, save :: no_drive_of_modes(:)
  logical, public, save :: erase_any_drive

  interface mode_write_nml
    module procedure mode_read_nml
  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

!-----------------------------------------------------------------------------
!> read (or write) the mode namelist
!-----------------------------------------------------------------------------

subroutine mode_read_nml(ilun,io_stat,lwrite)
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  use grid,    only : nmod
  use general, only : gkw_abort

  integer, intent(in)  :: ilun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  namelist /mode/ kthrho,     & ! the poloidal wave vector times rho
                & chin,       & ! poloidal angle shift
                & krrho,      & ! the radial wave vector times rho
                & kr_type,    & ! switch to use krrho or chin as an input
                & mode_box,   & ! true if a 2D grid of modes is used
                & krhomax,    & ! maximum krho used in 2D box
                & ikxspace,   & ! the spacing of the kx modes if shat/=0
                & rkxspace,   & ! the spacing of the kx modes if shat=0
                & n_spacing,  & ! toroidal mode number spacing (global runs)
                & no_transfer_to,  & !
                & no_transfer_from, & !
                & no_drive_of !
  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then

      ! test of nx is not too large
      if (nmod > nmmx) then
        call gkw_abort('nmod > nmmx in mode.F90. Reset nmmx and recompile')
      end if

      ! Set default values
      kr_type    = 'chin'
      chin       = 0.E0
      krrho      = 0.E0
      kthrho     = 0.E0
      mode_box   = .false.
      krhomax    = 0.E0
      ikxspace   = 0
      rkxspace   = 0.
      n_spacing  = -1
      no_transfer_to = 0
      no_transfer_from = 0
      no_drive_of = 0
      read(ilun,NML=mode,IOSTAT=io_stat)

    else
      ! do nothing
    end if
  else
    if(root_processor) write(ilun,NML=mode)

    ! write metadata
    call write_run_parameter('mode', 'kthrho', kthrho)
    call write_run_parameter('mode', 'chin', chin)
    call write_run_parameter('mode', 'krrho', krrho)
    call write_run_parameter('mode', 'kr_type', kr_type)
    call write_run_parameter('mode', 'mode_box', mode_box)
    call write_run_parameter('mode', 'krhomax', krhomax)
    call write_run_parameter('mode', 'ikxspace', ikxspace)
    call write_run_parameter('mode', 'rkxspace', rkxspace)
    call write_run_parameter('mode', 'n_spacing', n_spacing)
    call write_run_parameter('mode', 'no_transfer_to', no_transfer_to)
    call write_run_parameter('mode', 'no_transfer_from', no_transfer_from)
    call write_run_parameter('mode', 'no_drive_of', no_drive_of)
    
  end if

end subroutine mode_read_nml


!-----------------------------------------------------------------------------
!> broadcast the input parameters for mode to other processors
!-----------------------------------------------------------------------------

subroutine mode_bcast_nml()

  use mpiinterface, only : mpibcast
  use grid,         only : nmod

  call mpibcast(kthrho,       nmod)
  call mpibcast(kr_type, lenswitch)
  call mpibcast(chin,            1)
  call mpibcast(krrho,           1)
  call mpibcast(mode_box,        1)
  call mpibcast(krhomax,         1)
  call mpibcast(ikxspace,        1)
  call mpibcast(rkxspace,        1)
  call mpibcast(n_spacing,       1)
  call mpibcast(no_transfer_to, nmmx)
  call mpibcast(no_transfer_from, nmmx)
  call mpibcast(no_drive_of, nmmx)

end subroutine mode_bcast_nml


!-----------------------------------------------------------------------------
!> Check the parameters for mode. In the case when the chease interface is
!> used, at present this routine may be called twice.
!-----------------------------------------------------------------------------

subroutine mode_check_params(icall)

  use mpiinterface, only : root_processor
  use general,      only : gkw_warn, gkw_exit, gkw_abort
  use fft,          only : working_fft_library
  use grid,         only : nmod, nx, n_x_grid, nperiod
  use control,      only : non_linear, disp_x, disp_y, spectral_radius
  use geom,         only : shat, geom_type
  use constants,    only : pi
  use global,       only : r_tiny
  use components,   only : finit
  
  integer, intent(in) :: icall
  integer :: i

  ! do a few checks on the given input
  if (mode_box) then
    if (.not. working_fft_library) then
      call gkw_exit('mode_box requires a working FFT library')
    end if
    if (ikxspace <= 0 .and. abs(shat) > r_tiny .and. spectral_radius) then
      call gkw_exit('Unreasonable value of ikxspace')
    end if
    if (nperiod /= 1 .and. geom_type/='slab') then
      call gkw_warn('With mode_box = .true. only nperiod=1 is globally consistent!')
    end if
    if (root_processor .and. icall == 1) then
      write(*,*)
      write(*,*) 'With mode_box input value(s) of kthrho are ignored'
      write(*,*)
    end if

    ! Do some checks on value of shat with mode_box.
    ! These are repeated in case chease changed shat.
    if(abs(shat) < 1e-5) then  ! The zero shear case is treated differently
      lshat_zero = .true.
      ! Write actual values to file input.out
      shat = 0.
      
      if (spectral_radius) then
        call gkw_warn('with shat=0, rkxspace = Lx/Ly')

        if (abs(rkxspace) < 1e-5) then           
          call gkw_warn('Using rkxspace = Lx/Ly = 1')
          rkxspace = 1.
        end if
      end if
      
      if (ikxspace /= 0) then
           call gkw_warn('ikxspace not used with shat=0')
           ikxspace = 0
      end if
      
      if(root_processor) write(*,*) 'zero shear case selected'
    else if (abs(shat) < 0.05) then ! Near zero shear not implemented
      call gkw_exit('Magnetic shear must be exactly zero, '//               &
                    &'case close to zero not possible in flux tube.')
    else if (abs(shat) < 0.1) then
      call gkw_warn('Small magnetic shear requires v.large nx for mode_box')
    else
      ! In the "normal" case, do nothing
    end if
    
    if ((.not. lshat_zero) .and. abs(rkxspace) > 1e-5) then
      call gkw_warn('rkxspace not used unless shat=0')
      rkxspace = 0.
    end if 

    ! kthrho not used, say so in input.out
    kthrho = 0.

    if (abs(chin) > 1e-5) then
      call gkw_warn('chin not implemented for mode_box')
      chin = 0.
    end if

    if (abs(krrho) > 1e-5) then
      call gkw_warn('krrho not implemented for mode_box')
      krrho = 0.
    end if
    
    if ( (finit == 'zonal') .and. (spectral_radius) ) then
      if (nx.le.2) then
        call gkw_abort('NX must be at least 3 for init zonal')
      endif
    endif

  else ! not mode_box

    if (non_linear) call gkw_exit('mode_box must be true for nonlinear runs')

    if (.not.spectral_radius) &
      & call gkw_exit('mode_box must be true for spectral_radius = .false.')

    ! If mode_box=false, all kxrh=0, since kgrid returns.
    if (nx > 1) call gkw_exit('nx=1 should be used for mode_box=.false.')

    if (ikxspace /= 0) then
       call gkw_warn('ikxspace not used with mode_box=.false.')
       ikxspace = 0
    end if

    if (kr_type/='chin' .and. kr_type/='kr') then
      call gkw_abort('mode: unkown kr_type option, available are "chin" and "kr"')
    end if

    if (kr_type=='chin' .and. abs(chin)>r_tiny .and. nmod > 1) then
      call gkw_abort('mode: kr_type="chin" option only available for NMOD=1')
    end if

    if (kr_type=='chin' .and. abs(krrho)>r_tiny) then
      call gkw_warn('krrho value not used for kr_type="chin", krrho set to 0')
      krrho = 0.
    end if

    if (kr_type=='chin' .and. abs(chin) > pi) then
      call gkw_warn('chin value should be between -pi and pi, chin value '//  &
                   &'changed')
      chin = mod(chin,pi)
    end if

    if (kr_type=='kr' .and. abs(chin)>r_tiny) then
      call gkw_warn('chin value not used for kr_type="kr", chin set to 0')
      chin = 0.
    end if

    if (root_processor .and. icall == 1) then
      write(*,*)
      write(*,*) 'With mode_box off input value of krhomax is ignored'
      write(*,*) 'No 2D diagnostics will be written'
      write(*,*) 'Ensure you have ', nmod, ' value(s) in kthrho list'
      write(*,*)
    end if
    
    if (disp_x > 0. .or. disp_y > 0.) then
       call gkw_warn('Perpendicular dissipation is only with mode_box')
       disp_x = 0.
       disp_y = 0.
    end if
    
    if ( (finit == 'zonal') .and. (spectral_radius) ) then
      call gkw_abort('init zonal requires mode_box=.true.')
    endif

    ! krhomax not used, say so in input.out
    krhomax = 0.

  end if ! mode_box

  do i = 1, size(no_transfer_to), 2
    if(no_transfer_to(i) > nmod .or. no_transfer_to(i+1) > n_x_grid) then
      call gkw_abort('Invalid values in no_transfer_to list, >nmod &
         &or >n_x_grid')
    end if
    if(.not.spectral_radius .and. no_transfer_to(i+1) > 0) then
      call gkw_abort('Cannot suppress transfer to kx modes for nonspectral.')
    end if
  end do
  do i = 1, size(no_transfer_from), 2
    if(no_transfer_from(i) > nmod .or. no_transfer_from(i+1) > n_x_grid) then
      call gkw_abort('Invalid values in no_transfer_from list, >nmod &
         &or >n_x_grid')
    end if
    if(.not.spectral_radius .and. no_transfer_from(i+1) > 0) then
      call gkw_abort('Cannot suppress transfer from kx modes for nonspectral.')
    end if
  end do



end subroutine mode_check_params


!-----------------------------------------------------------------------------
!> Allocation of the arrays for mode
!-----------------------------------------------------------------------------

subroutine mode_allocate()

  use general, only : gkw_abort
  use grid,    only : nmod, nx, n_x_grid, ns, parallel_s

  integer :: ierr
  integer :: i,n

  ! initialize the error integer
  ierr = 0

  ! allocate the krloc array
  if (parallel_s) then
    allocate(krloc(nmod,nx,-1:ns+2),stat=ierr)
  else
    allocate(krloc(nmod,nx,ns),stat=ierr)
  end if
  if (ierr /= 0) call gkw_abort('Could not allocate krloc in mode')

  ! allocate the krho array
  allocate(krho(nmod),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate krho in mode')

  ! allocate the kxrh array
  allocate(kxrh(nx),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate kxrh in mode')

  ! allocate the arrays with the integers for the parallel boundary conditions
  allocate(ixplus(nmod,nx),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate ixplus in mode')
  allocate(ixminus(nmod,nx),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate ixminus in mode')

  ! allocate the mode label array
  allocate(mode_label_l(nmod,nx),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate mode_label_l in mode')
  allocate(mode_label_G(nmod,n_x_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate mode_label_G in mode')

  ! allocate the no_transfer_to/from_modes arrays:
  ! At first count all elements from the arrays no_transfer_to/from
  ! which are in the range of admissible indices.
  n = 0
  do i = 1, size(no_transfer_to), 2
    if(no_transfer_to(i) <= 0 .and. no_transfer_to(i+1) <= 0) cycle
    n = n + 1
  end do
  ! Now n contains the number of mode indices given by the input file.
  allocate(no_transfer_to_modes(2,n),stat=ierr)
  if (ierr /= 0) then
    call gkw_abort('Could not allocate no_transfer_to_modes in mode')
  end if
  
  n = 0
  do i = 1, size(no_transfer_from), 2
    if(no_transfer_from(i) <= 0 .and. no_transfer_from(i+1) <= 0) cycle
    n = n + 1
  end do
  allocate(no_transfer_from_modes(2,n),stat=ierr)
  if (ierr /= 0) then
    call gkw_abort('Could not allocate no_transfer_from_modes in mode')
  end if

  n = 0
  do i = 1, size(no_drive_of)
    if(no_drive_of(i) >= 1 .and. no_drive_of(i) <= nmod) then
      n = n + 1
    end if
  end do
  allocate(no_drive_of_modes(n),stat=ierr)
  if (ierr /= 0) then
    call gkw_abort('Could not allocate no_drive_of_modes in mode')
  end if

end subroutine mode_allocate


!-----------------------------------------------------------------------------
!> This subroutine initializes the parameters of mode. This subroutine could
!> be merged with another, perhaps kgrid.
!-----------------------------------------------------------------------------

subroutine mode_init()
  use grid, only : nmod, n_y_grid, nymod

  integer :: imod, i, j

  ! initialisations of grid sizes. These are saved in the grid module,
  ! but can only be determined if krhomax is known. Unfortunately it
  ! had been decided historically that krhomax is placed in the mode
  ! namelist.  For now, this means that the following grid sizes can be
  ! initialised only here instead of in the grid module:
  if(nmod == 1 .and. krhomax > 0.0) then
    nymod = 2
  else
    nymod = nmod
  end if
  n_y_grid = (nymod - 1) * 2 + 1
  ! because of nmod = floor(n_y_grid/2) + 1
  ! there are two position space gridsizes
  ! which lead to nmod modes in the kspace buffer. If one chooses
  ! the position space size without +1,
  ! then there is trouble with the smallest scale mode if one fourier-
  ! transforms e.g. a Poten000000n dataset and compares it to Spc3d000000n

  ! further initialisations:

  do imod = 1, nmod
    krho(imod) = kthrho(imod)
  end do

  erase_any_transfer_to = .false.
  erase_any_transfer_from = .false.
  erase_any_drive = .false.
  j = 0
  do i = 1, size(no_transfer_to), 2
    if(no_transfer_to(i) <= 0 .and. no_transfer_to(i+1) <= 0) cycle

    j = j + 1
    no_transfer_to_modes(1,j) = no_transfer_to(i)
    no_transfer_to_modes(2,j) = no_transfer_to(i+1)
    erase_any_transfer_to = .true.
  end do

  j = 0
  do i = 1, size(no_transfer_from), 2
    if(no_transfer_from(i) <= 0 .and. no_transfer_from(i+1) <= 0) cycle

    j = j + 1
    no_transfer_from_modes(1,j) = no_transfer_from(i)
    no_transfer_from_modes(2,j) = no_transfer_from(i+1)
    erase_any_transfer_from = .true.
  end do

  j = 0
  do i = 1, size(no_drive_of)
    if(no_drive_of(i) >= 1 .and. no_drive_of(i) <= nmod) then
      j = j + 1
      erase_any_drive = .true.
      no_drive_of_modes(j) = no_drive_of(i)
    end if
  end do
  
end subroutine mode_init

!-----------------------------------------------------------------------------
!> Calculates the grids for the 2D case and determines the integers necessary
!> for the parallel boundary conditions. Must be called after geom_init_grids
!> but before parallelize_geom. This routine can only be used for the local 
!> spectral flux tube version.
!> This routine must be called before parallelize geom due to way pol_angle
!> array is used
!-----------------------------------------------------------------------------
subroutine kgrid()
  use general,      only : gkw_abort
  use grid,         only : nmod, nx, n_s_grid, lx, n_x_grid, nxmod
  use geom,         only : q, eps, shat, kthnorm, signB, signJ
  use geom,         only : metric_G, pol_angle, geom_type, kxnorm
  use constants,    only : pi
  use mpiinterface, only : root_processor
  use global,       only : r_tiny
  use grid,         only : nperiod

  integer :: imod, ix, i, i_chin, j_chin, ierr, ixnext, swap_temp_int
  real :: kx, ky, kxplus, kxminus, kxhalf
  real, allocatable :: dum_s(:), dum_angle(:)

  ! Check that shat is correctly intialised
  if (shat > 1.e4) call gkw_abort('shat not correctly intialised. '//      &
      & 'geom_init_grids must be called before kgrid for geom type chease')
  ! Check that q is correctly intialised
  if (q > 1.e4) call gkw_abort('q not correctly intialised. '//            &
      & 'geom_init_grids must be called before kgrid for geom type chease')
!   ! Check that kthnorm is correctly intialised
  if (kthnorm > 1.e4) call gkw_abort('kthnorm not correctly intialised. '//&
      & 'geom_init_grids must be called before kgrid')

  ! initialize
  ixplus  = 0
  ixminus = 0
  kxrh    = 0.
  ixzero  = 0
  iyzero  = 0

  ! temporary arrays
  allocate(dum_s(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dum_s in mode')
  allocate(dum_angle(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dum_angle in mode')

  ! For mode_box = .false. the only function of this routine is to
  ! properly normalise the wave vectors given in input.
  if (.not. mode_box) then
    krho = krho / kthnorm
    kymax = maxval(krho)
    if (minval(abs(krho)) < r_tiny) iyzero = minloc(krho,1)

    select case(kr_type)
     case('kr')
       if (abs(krrho) .gt. r_tiny) then
         kxrh(1)= krrho  / kxnorm
       else ! (kr =0.)
         kxrh(1) = 0.
       end if

     case('chin')
       ! compute the value of kx to minimise k_perp @ pol_angle=CHIN
       ! In s-alpha corresponds to
       ! kxrh(1) = - chin* abs(q*shat*krho(imod=1)) / (2*pi*eps).
       if (abs(chin) .gt. r_tiny) then
         do i =  1, n_s_grid
          dum_angle(i) = pol_angle(1,i) - chin
          dum_s(i) = -metric_G(1,i,1,2) / metric_G(1,i,1,1)*krho(1)
         end do
         i_chin = minloc(abs(dum_angle),1)
         if (abs(pol_angle(1,i_chin) - chin) < r_tiny) then
           kxrh(1) = dum_s(i_chin)
         else
           if (pol_angle(1,i_chin) > chin) then
            j_chin = i_chin - 1
           else
            j_chin = i_chin + 1
           end if
           kxrh(1) = (dum_s(i_chin) - dum_s(j_chin)) /                        &
                & (pol_angle(1,i_chin) - pol_angle(1,j_chin)) *                &
                & (chin - pol_angle(1,i_chin)) + dum_s(i_chin)
                ! (linear interpolation)
         end if
       else ! (chin =0.)
         kxrh(1) = 0.
       end if
       krrho = kxrh(1)* kxnorm

     case default
        call gkw_abort('mode: unknown kr_type option')
     end select

    ! Incorrect if chin /= 0 !
    ixzero = 1 ! used (for example) in tearingmodes, diagnostic
    kxmax = 1. ! appears in kxrh/kxmax in hyper dissipation term
               ! This value is meaningless, but prevents div by 0
    iyzero = 0

    ! label the modes
    nmodes_l = 0
    do ix = 1, nx
      do imod = 1, nmod
        nmodes_l = nmodes_l + 1
        mode_label_l(imod,ix) = nmodes_l
      end do
    end do
    ! there is no x-parallelisation, hence these
    ! are the same:
    nmodes_G = nmodes_l
    mode_label_G = mode_label_l

    ! avoid unitialised values:
    lxinv = 0.0
    lyinv = 0.0
    lxn = 0.0
    lyn = 0.0

  else ! (mode_box)

    ! Calculate the 'toroidal' wave numbers. Note only nmod modes are
    ! used in the time intgration.
    if (nmod > 1) then
      do imod = 1, nmod
        krho(imod) = krhomax*(imod-1)/(real(nmod-1)*kthnorm)
      end do
    else
      krho(1) = krhomax/kthnorm
    end if
    if (abs(krho(1)) < r_tiny) then 
      iyzero = 1
    else
      ! then there is no kyzero mode
      iyzero = 0
    endif

    if (nmod > 1) then
      if (lshat_zero) then
        ! Lx / Ly = rkxspace
        kxspace = krho(2)*kthnorm/rkxspace
      else
        ! In order not to reverse the mode ordering need abs()
        ! In theory shat < 0 is the only part that could be negative.
        kxspace = abs(q*shat*krho(2) / (eps*real(ikxspace)))
      end if
    else
      if (lshat_zero) call gkw_abort('shat=0 not implemented for nmod=1')
      
      kxspace = abs(q*shat*krho(1) / (eps*real(ikxspace)))
    end if

    kxhalf = kxspace / 2.

    ! Determine the kx modes again. Not all these modes are used in time
    ! integration.

    ! NX must always be odd
    if (mod(nx,2) == 0) call gkw_abort('mode: NX must be odd')

    kxrh(1) = - real(nx-1)*kxspace / 2.
    ! the radial kspace grid has spacing kxspace; the zeromode is in
    ! the middle of the grid (i.e. not the first element as FFT
    ! libraries expect it)
    ixzero  = (n_x_grid+1)/2
    ! nxmod is the number of radial modes, with +k and -k counted as
    ! one, not two modes. Hence nxmod is analogous to nmod.
    nxmod = (n_x_grid+1)/2
    ! is the same as
    !nxmod = floor(n_x_grid/2) + 1 ! (cf FFT library docs)
    !nxmod = n_x_grid - ixzero - 1)

    do ix = 2, nx
      kxrh(ix) = kxrh(ix-1) + kxspace
    end do

    ! Calculate the box size in real space.
    ! WARNING: still need to look at nmod = 1 case.
    if (nmod > 1) then
      ly = 2.*pi/krho(2)
      lx = 2.*pi/kxrh(ixzero+1)
    else
      lx = 1.
      ly = 1.
    end if

    lxinv = 1./lx
    lyinv = 1./ly
    
    ! box sizes in real space at LFS in units of rho_ref
    lyn=ly/kthnorm
    lxn=lx/kxnorm
   
    if (root_processor) then
      write(*,*)
      write(*,'(A,F8.2,A,F8.2)') ' LFS real box size: lxn/rho_ref ', lxn, ' lyn/rho_ref ', lyn
      write(*,*)
    end if

    kxmax = maxval(kxrh)
    kxmin = minval(kxrh)
    kymax = maxval(krho)

    ! Make the integer connections for the parallel boundary conditions.
    do imod = 1, nmod
      do ix = 1, nx
        ky = krho(imod)
        kx = kxrh(ix)

        if (geom_type == 'slab') then !No shear-periodicity in the slab
          ixminus(imod,ix) = 0
          ixplus(imod,ix) = 0

        ! ky = 0 mode is always treated differently
        else if (imod == iyzero) then
          ! the ky = 0 mode is always periodic, connects to itself
          ixminus(imod,ix) = ix
          ixplus(imod,ix) = ix

        else

          ! kx value after nperiod poloidal turns
          kxplus = kx + (2*nperiod-1)*abs(q*shat*ky/eps)
          if (kxplus > kxmax + kxhalf) then
            ixplus(imod,ix) = 0
          else

            ! inefficient programming, but very general ...
            i = 1
            find_ixplus : do
              if (abs(kxplus-kxrh(i)) <= 0.5*kxhalf) exit find_ixplus
              i = i + 1
              if (i > nx) call gkw_abort('Severe internal error: ixplus')
            end do find_ixplus
            ixplus(imod,ix) = i

          end if

          ! kx value before nperiod poloidal turns
          kxminus = kx - (2*nperiod-1)*abs(q*shat*ky/eps)
          if (kxminus < kxmin - kxhalf) then
            ixminus(imod,ix) = 0
          else

            ! inefficient programming, but very general ...
            i = 1
            find_ixminus : do
              if (abs(kxminus-kxrh(i)) <= 0.5*kxhalf) exit find_ixminus
              i = i + 1
              if (i > nx) call gkw_abort('Severe internal error: ixminus')
            end do find_ixminus
            ixminus(imod,ix) = i

          end if

          ! Swap the connections over for dqdpsi*ky < 0.
          ! Here dqdpsi=signB*signJ*q*shat/eps is a signed quantity
          if (signB*signJ*q*shat*ky/eps < 0.) then
            swap_temp_int = ixminus(imod,ix)
            ixminus(imod,ix) = ixplus(imod,ix)
            ixplus(imod,ix) = swap_temp_int
          end if

        end if

      end do
    end do

    ! label the independent modes locally unambiguously
    nmodes_l = 0
    mode_label_l(:,:) = 0
    do imod = 1, nmod
      do ix = 1, nx

        if (mode_label_l(imod,ix) == 0) then
          nmodes_l = nmodes_l + 1
          mode_label_l(imod,ix) = nmodes_l
          ixnext = ixplus(imod,ix)

          ! mode can be connected to itself
          if (ixnext /= ix) then
            get_next_x : do
              if (ixnext > 0) then
                mode_label_l(imod,ixnext) = nmodes_l
                ixnext = ixplus(imod,ixnext)
              else
                exit get_next_x
              end if
            end do get_next_x
          end if

        end if

      end do
    end do

    ! for spectral runs, there is no x-parallelisation, hence these
    ! are the same:
    nmodes_G = nmodes_l
    mode_label_G = mode_label_l

  end if ! mode_box

  ! deallocate the temporary arrays
  deallocate(dum_s)
  deallocate(dum_angle)

end subroutine kgrid


!-----------------------------------------------------------------------------
!> Calculates the grids for the 2D case and determines the integers necessary
!> for the parallel boundary conditions. Must be called after geom_init_grids
!> but before parallelize_geom (due to way qx and shatx are used)!
!-----------------------------------------------------------------------------
subroutine kgrid_nonspec_rad()

  use control,      only : flux_tube, radial_boundary_conditions
  use general,      only : gkw_abort, gkw_warn
  use grid,         only : nmod, nx, lx, n_x_grid, nxmod
  use geom,         only : qx, eps, shatx, kthnorm
  use geom,         only : geom_parallelized, kxnorm
  use constants,    only : pi
  use components,   only : rhostar
  use mpiinterface, only : root_processor
  use global,       only : r_tiny

  integer :: imod, ix, ikxspace_floor, ikxspace_ceil
  real    :: q, dum, shat, krho_min = 0.0, lx_kx_shift = 0.0

  ! this can be fixed by not using x-global qx and shatx grids below
  if (geom_parallelized) call gkw_abort('kgrid_nonspec_rad too late')

  do ix = 1, n_x_grid 

    ! Check that shat is correctly intialised
    if (shatx(ix) > 1.22e4) call gkw_abort('shat not correctly intialised. '//&
        & 'geom_init_grids must be called before kgrid for geom type chease')

    ! Check that q is correctly intialised
    if (qx(ix) > 1.22e4) call gkw_abort('q not correctly intialised. '//      &
        & 'geom_init_grids must be called before kgrid for geom type chease')

  end do 

  ! Check that kthnorm is correctly intialised
  if (kthnorm > 1.22e4) call gkw_abort('kthnorm not correctly intialised. '// &
      & 'geom_init_grids must be called before kgrid')

  ! initialize
  ixplus  = 0
  ixminus = 0
  kxrh    = 0.E0

  ! set the q and shear value to the value at the centre of the box 
  ! works for even and odd n_x_grid
  q    = (qx((n_x_grid+1)/2)+qx((n_x_grid+2)/2))/2.
  shat = (shatx((n_x_grid+1)/2)+shatx((n_x_grid+2)/2))/2.

  ! If spectral_radius = .false. the code must be run with modebox is true. 
  if (.not. mode_box) call gkw_abort('spectral_radius=T requires mode_box = T')

  ! Calculate the 'toroidal' wave numbers. Note only nmod modes are
  ! used in the time integration.
  if ((n_spacing < 0).or.(flux_tube)) then 
    if (nmod > 1) then
      do imod = 1, nmod
        krho(imod) = krhomax*(imod-1)/real(nmod-1)/kthnorm
      end do
      iyzero = 1
    else
      krho(1) = krhomax/kthnorm
      !write(*,*) krho(1)
      if (abs(krho(1)) < r_tiny) then 
        iyzero = 1 
      else 
        iyzero = 0
      endif
    end if
  else 
    ! use actual mode numbers
    if (nmod == 1) then 
      krho(1) = 2*pi*rhostar*n_spacing 
      if (n_spacing == 0) then 
        iyzero = 1 
      else 
        iyzero = 0 
      endif 
    else 
      do imod = 1, nmod
        krho(imod) = 2*pi*rhostar*n_spacing * (imod - 1)  
      end do
      iyzero = 1
    endif 
  endif 

  ! set the mode box in the y-direction 
  ly = 1.
  if (nmod > 1) then
    ly = 2.*pi/krho(2)
  end if
  lyinv = 1./ly
    
  ! Binormal box size in units of rho_ref
  lyn=ly/kthnorm
  
  ! for the nonspectral case there is no zero radial mode as the
  ! fields are in position space represenation with respect to the
  ! radial direction.
  ixzero = 1
  ! If one would do a fft of a field with respect to x, then:
  ! nxmod the number of radial modes, with +k and -k counted as one, not two
  ! modes. Hence nxmod is analogous to nmod.
  !ixzero = (n_x_grid+1)/2 ! but some code may rely on ixzero = 1 for nonspec
  nxmod = (n_x_grid+1)/2
  
  if (nmod > 1 ) krho_min = krho(2)
  if (nmod == 1) krho_min = krho(1)

  ! For flux tube, lx must still obey ballooning periodicity ?
  ! connection criterion, in order for radial boundary conditions
  ! to work correctly
  if (ikxspace < 0) then
    !There is a problem with this: LX cannot be reset at this
    !stage because it has already been used in geom and dist.
    call gkw_abort('Negative ikxspace feature needs fixing')
    ! use ikxspace when lx <= 0
    if (lx <=0. .and. flux_tube .and. .not. lshat_zero) then
      kxspace = abs(q*shat*krho_min / (eps*real(abs(ikxspace))))
      lx = 2.*pi/kxspace  ! CAN'T DO THIS NOW
    ! otherwise set ikxspace appropriately
    else if (lx >0. .and. flux_tube .and. .not. lshat_zero) then
      dum = abs(q*shat*krho_min / eps)
      ikxspace = nint(lx*dum/(2.*pi))
      kxspace = dum / ikxspace
      lx = 2.*pi/kxspace  ! CAN'T DO THIS NOW
    end if    
  end if  

  ! check if the value of lx has been properly set 
  ! if (lx <= 0.) call gkw_abort('Zero or negative value of lx')
  
  ! Check if radial boxsize is consistent with the toroidal and poloidal
  ! periodicity constraint in field aligned Hamada coordinates with 
  ! sheared background magnetic field and radial periodic boundary 
  ! conditions.
  if((flux_tube) .and. (.not. lshat_zero) .and. &
    & (radial_boundary_conditions == 'periodic')) then
    
    ! Wavelength connected to radial wave vector shift at parallel 
    ! boundary.
    lx_kx_shift = 2*pi*eps/abs(q*shat*krho_min)
    
    ! Nearest (to lx) integer multiples of the wavelength connected to 
    ! the radial wave vector shift. The former must be an integer 
    ! multiple of the latter to be consistent with periodicity 
    ! constraint.
    ikxspace_floor = floor(lx/lx_kx_shift)
    ikxspace_ceil = ceiling(lx/lx_kx_shift)
    
  
    ! Check if radial wave vector shift over parallel boundary is 
    ! consistent with the radial boxsize specified as input parameter.
    ! The accuracy is set to 1e-5.
    if(abs(lx-ikxspace_floor * lx_kx_shift) > 1e-5 .and. &
      & abs(lx-ikxspace_ceil * lx_kx_shift) > 1e-5) then
      
      call gkw_warn('The given radial boxsize lx is not consistent with&
      & the toroidal and poloidal periodicity constraint in field&
      & aligned Hamada coordinates with sheared background magnetic&
      & field and periodic radial boundary conditions. ')
      
      if(root_processor) then
        write(*,*)
        write(*,*) 'To prevent radial discontinuities set lx to an&
        & integer multiple of: ', 2*pi*eps/abs(q*shat*krho_min)
      end if
      
    end if
  
  end if
 
  ! rescale input lx psi coordinate space to real space
  lxn = lx / kxnorm  
  lxinv = 1./lx 

  ! lyn=ly_real, lyn=lx_real.  TO DO: These quantities could be usefully written to geom.dat
  if (root_processor) then
    write(*,*)
    write(*,'(A,F8.2,A,F8.2)') ' LFS real box size / rho_ref: lxn', lxn, ', lyn', lyn
  ! write(*,*) 'ikxspace: ', ikxspace
    write(*,*)
  end if

  kymax = maxval(krho)

  ! Make the integer connections for the parallel boundary conditions.
  ! Every radial mode connects to itself
  do imod = 1, nmod
    do ix = 1, nx
      ixminus(imod,ix) = ix
      ixplus(imod,ix) = ix
    end do 
  end do 
  
  ! label the independent modes, globally
  nmodes_G = 0
  do imod = 1, nmod
    do ix = 1, n_x_grid
      nmodes_G = nmodes_G + 1
      mode_label_G(imod,ix) = nmodes_G
    end do
  end do

  ! label the independent modes, locally
  nmodes_l = 0
  do imod = 1, nmod
    do ix = 1, nx
      nmodes_l = nmodes_l + 1
      ! correct would be mode_label_l(imod,ix) =
      !mode_label_G(imod,gx(ix)) but there is a bug in
      !diagnos_growth_freq and the above would cause a
      !segfault. diagnos_nonlin_transfer expects this, too, as it uses
      !it as an index for a local array:
      mode_label_l(imod,ix) = nmodes_l
    end do
  end do
  ! now mode_label_l still holds nmodes_l different numbers, but they
  ! do not go from 1 to nmodes_l

end subroutine kgrid_nonspec_rad


!-----------------------------------------------------------------------------
!> This subroutine calculates some quantities that are necessary for the
!> ballooning transform. krloc is the local perpendicular wave vector. 
!> This routine must be called after parallelize_geom and kgrid.
!-----------------------------------------------------------------------------

subroutine krbal()

  use general, only : gkw_abort
  use grid,    only : nx, ns, nmod, parallel_s
  use geom,    only : metric, geom_parallelized

  integer :: ix, i, imod, i1, i2

  if (.not. geom_parallelized) call gkw_abort('krbal before parallelize geom')

  if (parallel_s) then
    i1= -1 ; i2 = ns + 2
  else
    i1=  1 ; i2 = ns
  end if

  ! Calculate the peperpendicular wave vector
  do imod = 1, nmod
    do ix = 1, nx
      do i = i1, i2
        krloc(imod,ix,i) = krho(imod)**2*metric(ix,i,2,2) + 2.*krho(imod)*      &
                         & kxrh(ix)*metric(ix,i,1,2) + kxrh(ix)**2*metric(ix,i,1,1)
        krloc(imod,ix,i) = sqrt(krloc(imod,ix,i))
      end do
    end do
  end do

end subroutine krbal


!-----------------------------------------------------------------------------
!> Small routine that allows one to plot quantities along the field line when
!> mode_box = true is used. At present puts on file krloc_G (the perpendicular
!> wave number), the curvature, and the coriolis operator. This must be called
!> before parallelize geom shifts the tensors in s.
!-----------------------------------------------------------------------------

subroutine mode_box_recon()

  use general,      only : gkw_abort
  use grid,         only : nmod, n_s_grid
  use geom,         only : metric_G, dfun, hfun, geom_parallelized
  use mpiinterface, only : root_processor
  use control,      only : spectral_radius
  
  integer :: imod, i, ixmin, ixref
  real :: curv, cori, krloc_G

  if (geom_parallelized) then
    call gkw_abort('call mode_box_recon before parallelize_rotation')
  end if

  if (.not. root_processor) return

  if (.not. spectral_radius) return

  if (mode_box) then

    ! Find the minimum kx
    kxmin = minval(kxrh)
    ixmin = minloc(kxrh,1)

    ! Curvature and kperp and coriolis for modes connected to kx=0 only
    ! Same can be achieved made using the mode_label formalism?
    open(13,file = 'par.dat')
    do imod = 1, nmod
      write(13,11) imod
      11 format('The toroidal mode ',I4)
      write(13,12)
      12 format(' Perpend. k-vec  Curvature        Coriolis')

      ! start at minimum kx
      ixref = ixmin
      do
        do i = 1, n_s_grid
          curv = dfun(1,i,1)*kxrh(ixref) + dfun(1,i,2)*krho(imod)
          cori = hfun(1,i,1)*kxrh(ixref) + hfun(1,i,2)*krho(imod)
          krloc_G = krho(imod)**2*metric_G(ixref,i,2,2) + 2.*krho(imod)*             &
             & kxrh(ixref)*metric_G(ixref,i,1,2) + kxrh(ixref)**2*metric_G(ixref,i,1,1)
          krloc_G = sqrt(krloc_G)
          write(13,fmt = '(30(es16.8,1X))') krloc_G,curv,cori
        end do
        if (ixplus(imod,ixref) /= 0 .and. ixref /= ixplus(imod,ixref)) then
          ixref = ixplus(imod,ixref)
        else
          exit
        end if
      end do
    end do
    close(13)

!   open(13,file = 'lxly.dat')
!     write(13,*) 's_point, kxmin, kymin, lx_perp, ly_perp, lx, ly'
!     do i = 1, n_s_grid
!       write(13,21) i, krloc_G(1,ixzero+1,i), krloc_G(2,ixzero,i), 2.*pi/krloc_G(1,ixzero+1,i), 2.*pi/krloc_G(2,ixzero,i), lx, ly
!       21 format(i,30(es12.4,1X))
!
!     end do
!   close(13)

  end if

end subroutine mode_box_recon

!------------------------------------------------------------------------------
!> This function calculates the transformation from one grid point along the 
!> field to another. A phase shift is needed
!>    * in the case of shifted metric whereas the multiplication factor
!>      is 1.0 for the non shifted case.
!>    * due to periodicity when a point is outside the s-grid.
!> As this routine needs the information if the considered point is outside the 
!> s-grid, it must be called before connect_parallel().
!>
!> Must be called after parallelize_geom
!------------------------------------------------------------------------------
function parallel_phase_shift(imod,ix,ip,il) 

  use control,   only : spectral_radius  
  use geom,      only : alphak, shift_end_grid 
  use constants, only : ci1 
  use grid,      only : n_s_grid, gs
  !> index of the binormal mode 
  integer,intent(in) :: imod
  !> index of the radial surface 
  integer,intent(in) :: ix
  !> parallel index at which one wants to calculate the distribution function
  integer,intent(in) :: ip
  !> parallel index at which the distribution function is given
  integer,intent(in) :: il
  complex :: parallel_phase_shift 

  if (spectral_radius) then 

    ! No phase shift if not using shifted metric 
    parallel_phase_shift = (1.E0,0.E0)

  else 
    
    ! Shift due to shifted metric (if not shift_metric, alphak= 0.):
    parallel_phase_shift = exp(ci1*krho(imod)*(alphak(ix,ip) & 
                                       &     - alphak(ix,il)))

    ! is a shift across the boundary necessary ? 
    if (gs(il).gt.n_s_grid) then 
      parallel_phase_shift = parallel_phase_shift * & 
                           & exp( -ci1*krho(imod)*shift_end_grid(ix) ) 
    endif 

    if (gs(il).lt.1) then 
      parallel_phase_shift = parallel_phase_shift * & 
                           & exp( +ci1*krho(imod)*shift_end_grid(ix) ) 
    endif 

  endif 

end function parallel_phase_shift 

!------------------------------------------------------------------------------

end module mode
