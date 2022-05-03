!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Provide the local grid sizes and associated quantities.
!> Decides how the global grid will be distributed across the processors
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module grid

  implicit none

  private

  !
  ! public procedures
  !

  public :: grid_read_nml, grid_bcast_nml, grid_check_params, grid_write_nml
  public :: setup_grid, gsp, gmu, gvpar, gs, gx, lsp, lmu, lvpar, ls, lrx
  public :: mpisetup_mod, proc_subset
  public :: jind_flexible, jinv_flexible

  !
  ! public variables
  !

  ! the *global* sizes of the grids (only for parallelizable directions)
  integer, save, public :: number_of_species !< total number of species
                                  ! (N.B. do not count the adiabatic species)
  integer, save, public :: n_mu_grid    !< total number of magnetic moment grid points
  integer, save, public :: n_vpar_grid  !< total number of vparallel grid points
  integer, save, public :: n_s_grid     !< total number of s grid points
  integer, save, public :: n_x_grid     !< total number of x grid points

  integer, save, public :: n_total      !< total number of grid points.
  integer, save, public :: n_field_total!< total number of points for the fields.
  
  integer, save, public :: nsg_pt       !< the number of s points per turn
  ! the sizes of the *local* grids (per processor)
  ! These are derived from the global sizes in setup_grid
  integer, save, public :: ns    !< number of grid points along the field line
  integer, save, public :: nmu   !< number of grid points in the magnetic moment
  integer, save, public :: nvpar !< number of grid points for the parallel velocity
  integer, save, public :: nx    !< number of radial wave vectors
  integer, save, public :: nmod  !< total number of toroidal modes
  integer, save, public :: nsp   !< the number of species

  !> total number grid points of the full binormal grid, (corresponds
  !> to n_x_grid, just for y direction)
  integer, save, public :: n_y_grid

  !> number of grid points on binormal k-grid, always including zero-mode
  integer, save, public :: nymod

  !> number of radial modes (corresponds to nymod, just for x direction)
  integer, save, public :: nxmod

  !> value specifies the length of the field line, i.e. the field line makes
  !> 2*nperiod - 1 poloidal turns. For nonlinear runs nperiod should be 1.
  integer, save, public :: nperiod

  !> The number of grid points in the trapping region with positive parallel
  !> velocity (used when vp_trap = 1). The total number in the trapped domain
  !> is 2*n_trapped.
  integer, save, public :: n_trapped

  !>The maximum value of the parallel velocity grid 
  real, save, public :: vpmax=3.E0

  !> The maximum value fo the mu grid 
  real, save, public :: mumax=3**2/2.E0

  !> The lower boundary of the radial domain (global only) 
  real, save, public :: psil 

  !> The upper boundary of the radial domain (global only) 
  real, save, public :: psih 

  !> The size of the radial domain normalized to the Larmor radius (for flux 
  !> tube simulations only, overwritten in global simulations) 
  real, save, public :: lx 

  ! The number of processors to be used in each parallel direction;
  ! these are either set automatically (if possible) in setup_grid
  ! (based on the total number of processors), or are prescribed via the
  ! GRIDSIZE namelist (previously in CONTROL).
  integer, save, public :: n_procs_s      !< # of procs to divide up s-direction
  integer, save, public :: n_procs_x      !< # of procs to divide up x-direction
  integer, save, public :: n_procs_mu     !< # of procs to divide up mu grid
  integer, save, public :: n_procs_sp     !< # of procs to divide up species
  integer, save, public :: n_procs_vpar   !< # of procs to divide up parallel velocity

  ! Logicals to know if we are parallel over some dimension; these are
  ! set in setup_grid
  logical, save, public :: parallel_sp    !< T if the code parallelizes over species
  logical, save, public :: parallel_mu    !< T if the code parallelizes over mu grid
  logical, save, public :: parallel_s     !< T if the code parallelizes over s grid
  logical, save, public :: parallel_vpar  !< T if the code parallelizes over vpar
  logical, save, public :: parallel_x     !< T if the code parallelizes over x

  !> rank of the next processor in the vpar direction
  integer, save, public :: proc_vpar_next
  !> rank of the previous processor in the vpar direction
  integer, save, public :: proc_vpar_prev
  !> rank of the next processor in the mu direction
  integer, save, public :: proc_mu_next
  !> rank of the previous processor in the mu direction
  integer, save, public :: proc_mu_prev
  !> rank of the next processor in the s direction
  integer, save, public :: proc_s_next
  !> rank of the previous processor in the s direction
  integer, save, public :: proc_s_prev
  !> rank of the previous processor in the x direction  
  integer, save, public :: proc_x_next
  !> rank of the previous processor in the x direction
  integer, save, public :: proc_x_prev
  !> rank of the next processor in the vpar and mu
  integer, save, public :: proc_vpar_next_mu_next
  !> rank of the previous processor in the vpar direction and next in mu
  integer, save, public :: proc_vpar_prev_mu_next
  !> rank of the next processor in the vpar and prev in mu
  integer, save, public :: proc_vpar_next_mu_prev
  !> rank of the previous processor in the vpar direction and prev in mu
  integer, save, public :: proc_vpar_prev_mu_prev
  !> rank of the next processor in the vpar and s
  integer, save, public :: proc_vpar_next_s_next
  !> rank of the previous processor in the vpar direction and next in s
  integer, save, public :: proc_vpar_prev_s_next
  !> rank of the next processor in the vpar and prev in s
  integer, save, public :: proc_vpar_next_s_prev
  !> rank of the previous processor in the vpar direction and prev in s
  integer, save, public :: proc_vpar_prev_s_prev
  !> rank of the opposite processor in the vpar direction in COMM_VPAR_NE
  integer, save, public :: proc_vpar_opposite
  

  !> send and recv boundary in the s direction
  logical, save, public :: lsendrecv_s = .false.
  !> send and recv boundary in the x direction (non spectral only)
  logical, save, public :: lsendrecv_x = .false.
  !> send and recv boundary in the vpar direction
  logical, save, public :: lsendrecv_vpar = .false.
  !> send and recv boundary in the mu direction (used with e.g. collisions)
  logical, save, public :: lsendrecv_mu = .false.

  !> position of the processor in the x direction 
  integer, save, public :: iproc_x

  !
  ! locals
  !

  !> number of possible parallel directions
  integer, parameter :: nplan = 5

  !> plan ordering (to optimise MPI cartesian topology)
  !> Tests on CRAY XT6 (Hector phase 2b) with N_vpar_grid=48=2*cores/node and
  !> ivpar=1, isp=2, imu=3, is=4, n_procs=1344
  !> Come out 30% slower than the default configuration below.
  !> Reversing the order below was 80% slower.
  !> The bottleneck is in the mu-sp reduction and the order
  !> below is probably optimal for most situations. If altered, it may be
  !> desirable to also reorder the index function correspondingly.
  !> Most machines place the first dimensions on local nodes where possible
  !> Optimal performance therefore also depends on gridsizes
  !> With the combination below: number_of_species*n_mu_grid = cores per node 
  !> might be the ideal
  
  integer, parameter :: isp = 1, imu = 2, ivpar = 3, is = 5, ix = 4 

  !> Convenient way of dealing with the plan in this module only.
  !> Each plan corresponds to a parallel direction.
  type :: plan
    logical :: parallel = .false.      !< true if parallel
    integer :: procs    = 0            !< number of procs
    integer :: points   = 0            !< number of points per proc
    integer :: min_points   = 0        !< minimum number of points per proc
    integer :: max_points   = 0        !< max. # of points per proc in parallel
    character (len=8) :: name = 'none' !< name
    integer :: ipb                     !< first point on this proc
    integer :: ipe                     !< last point on this proc
    integer :: iproc                   !< processor in the ? direction
    logical :: cart = .false.          !< if the direction needs local comm
    logical :: periodic = .false.      !< if the direction is periodic
    integer :: direction = -1          !< direction for MPI calls in cartesian
  end type plan

  !> structure to hold parallel plan details and associated indicies
  type (plan), save, dimension(nplan) :: pp !< pp is the plan structure

  !> re-order allowed in creating sub-communicators
  logical, parameter :: reorder = .true.

  !> ixpb contains the first x point on the local processor
  integer, save :: ixpb
  !> ixpe contains the last x point on the local processor
  integer, save :: ixpe
  !> isppb contains the first of the species number on the local processor
  integer, save :: isppb
  !> isppe contains the last of the species number on the local processor
  integer, save :: isppe
  !> imupb contains the first mu grid point on the local processor
  integer, save :: imupb
  !> imupe contains the last mu grid point on the local processor
  integer, save :: imupe
  !> ispb contains the first s grid point on the local processor
  integer, save :: ispb
  !> ispe contains the last s grid point on the local processor
  integer, save :: ispe
  !> ivparpb contains the first vpar grid point on the local processor.
  !> For vp_trap = 1, the processor will be responsible for a mirror point
  !> for every point apparently in the grid.
  integer, save :: ivparpb
  !> ivparpe contains the last vpar grid point on the local processor.
  !> For vp_trap = 1, the processor will be responsible for a mirror point
  !> for every point apparently in the grid.
  integer, save :: ivparpe
  
  !> arrays for the mpi_gatherv over nmod in the nonspectral field solve
  integer, save, allocatable, public :: gathv_mod_disps(:), gathv_mod_recvno(:) 

  ! other locals
  integer, save :: n_dims
  integer, save :: i
  integer, save :: ierr = 0
  logical, save, dimension (nplan) :: periodic
  integer, save, dimension (nplan) :: dims,coords

  !
  ! interfaces
  !

  ! use the same routine for reading and writing
  interface grid_write_nml
    module procedure grid_read_nml
  end interface

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> In this subroutine it is determined how the calculation is divided
!> among the different processors. This can be done automatically, by
!> setting all the processor numbers to 1 in the input file (or not
!> including them) or manually, by specifying n_procs_[mu|sp|s|vpar] in
!> the control namelist. Manual setup is implied if *any* of n_procs_*
!> are found to be > 1.
!> Additional communicators are set up for a mpi cartesian topology if necessary.
!> 
!> Nota bene this routine acts by side-effect:
!> It sets the variables nsp,ns,nmu and nvpar. These are the local grid sizes
!> on each processor.
!>  (which describe the size of the subproblem of the local process)
!----------------------------------------------------------------------------
subroutine setup_grid

  use control,      only : order_of_the_scheme, vp_trap, lcollisions
  use control,      only : spectral_radius 
  use control,      only : shear_periodic
  use global,       only : i_huge
  use general,      only : gkw_abort
  use mpiinterface, only : number_of_processors, root_processor
  
  ! make nothing parallel by default
  pp(:)%parallel   = .false.

  ! put some value in by default
  pp(:)%max_points = i_huge

  !
  ! put the input parameters into the plan structure
  !

  ! species:
  pp(isp)%name       = 'species'
  pp(isp)%points     = number_of_species
  pp(isp)%procs      = n_procs_sp
  pp(isp)%min_points = 1

  ! mu grid:
  pp(imu)%name       = 'mu'
  pp(imu)%points     = n_mu_grid
  pp(imu)%procs      = n_procs_mu
  pp(imu)%min_points = 1
  if (lcollisions) pp(imu)%min_points = 2

  ! vpar grid:
  pp(ivpar)%name       ='vpar'
  pp(ivpar)%points     = n_vpar_grid
  pp(ivpar)%procs      = n_procs_vpar
  if (vp_trap == 1) then
    pp(ivpar)%min_points = n_vpar_grid
  else
    select case(order_of_the_scheme)
      case('second_order')
        pp(ivpar)%min_points = 1
        pp(ivpar)%max_points = n_vpar_grid
      case('fourth_order')
        pp(ivpar)%min_points = 2
        pp(ivpar)%max_points = n_vpar_grid
      case default
        pp(ivpar)%min_points = i_huge
    end select
  end if

  ! s grid:
  pp(is)%name       ='s'
  pp(is)%points     = n_s_grid
  pp(is)%procs      = n_procs_s
  pp(is)%min_points = 2
  
  ! x grid:
  pp(ix)%name       ='x'
  pp(ix)%points     = n_x_grid
  pp(ix)%procs      = n_procs_x
  pp(ix)%min_points = n_x_grid
  if (.not. spectral_radius) pp(ix)%min_points = 2

  lsendrecv_x  = shear_periodic
  
  !write(*,*) pp(:)%procs
  ! Use the above to determine the number of processors to be used in each
  ! direction.
  call parallel_plan

  ! Copy the number of processors obtained into the variables use by the
  ! rest of the code and check the layout before doing any more.
  n_procs_s   = pp(is)%procs
  n_procs_vpar= pp(ivpar)%procs
  n_procs_mu  = pp(imu)%procs
  n_procs_sp  = pp(isp)%procs
  n_procs_x   = pp(ix)%procs

  ! set nsp, nmu, nvpar, ns
  nsp   = pp(isp)%points
  nmu   = pp(imu)%points
  nvpar = pp(ivpar)%points
  ns    = pp(is)%points
  nx    = pp(ix)%points

  ! set parallel_* logicals
  parallel_sp      = pp(isp)%parallel
  parallel_mu      = pp(imu)%parallel
  parallel_vpar    = pp(ivpar)%parallel
  parallel_s       = pp(is)%parallel
  parallel_x       = pp(ix)%parallel

  !Also check allowed inputs with parallel options
  !Do not move above logicals above.
  call check_parallel_layout

  ! enable sendrecv dist in the mu-direction with collisions
  lsendrecv_mu = parallel_mu .and. lcollisions
  lsendrecv_s  = parallel_s
  lsendrecv_x  = parallel_x .or. lsendrecv_x
  lsendrecv_vpar = parallel_vpar .and. (vp_trap == 0)

  ! Configure the processors for nearest neighbour communication of the
  ! boundaries in some of the parallel directions and determine the first
  ! and last points the local processor is responsible for.
  call setup_proc_coords
  call setup_communicators
  call get_local_ranks
  call setup_local_grid_ranges
  call set_comm_sizes()
  call openmpi_fix

  n_total      =number_of_species*n_x_grid*nmod*n_s_grid*n_mu_grid*n_vpar_grid
  n_field_total=number_of_species*n_x_grid*nmod*n_s_grid

  ! explain (via stdout) how things will be parallelized
  if (number_of_processors > 1) call parallel_report(pp)

  nsg_pt = n_s_grid/(2*nperiod-1)

  if (root_processor .neqv. proc_subset(1,1,1,1,1)) then
    call gkw_abort('This grid decomposition could lead to I/O problems')
  end if

end subroutine setup_grid

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Read the gridsize namelist
!----------------------------------------------------------------------------
subroutine grid_read_nml(ifile,io_stat,lwrite)
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  integer, intent(in)  :: ifile
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  !> *request* non-blocking boundary sendrecv in vpar direction
  !> Is obsolete now.
  !APS: remove non_blocking_vpar when the namelists can be cleanly updated
  logical :: non_blocking_vpar = .true.

  namelist /gridsize/ nx, nperiod, n_vpar_grid, nmod, number_of_species,     &
                    & n_mu_grid, n_trapped, n_s_grid, n_procs_s, n_procs_sp, &
                    & n_procs_vpar, n_procs_mu, non_blocking_vpar, vpmax,    &
                    & mumax, psil, psih, lx, n_procs_x !,n_x_grid 

  ! read nml, if not writing
  if (present(lwrite)) then

    if (.not. lwrite) then

      ! Set the default sizes of the grid to zero
      number_of_species = 0
      n_mu_grid         = 0
      n_vpar_grid       = 0
      n_s_grid          = 0
      n_x_grid          = 0

      lx   = 0
      psih = 0
      psil = 0

      ! Set other default values
      nperiod           = 1
      n_trapped         = 0
      nx                = 1
      nmod              = 1

      ! import any values set in control for the defaults
      n_procs_s         = 1
      n_procs_x         = 1
      n_procs_vpar      = 1
      n_procs_mu        = 1
      n_procs_sp        = 1

      ! read nml
      read(ifile,NML=gridsize,IOSTAT=io_stat)

    end if

  else

    ! Historically, 'nx' is the radial global grid size in the input
    ! file.  When the input.out file and the metadat are written here,
    ! the variable nx has already become the *local* grid size for every process.
    
    ! Workaround, to have nx in input.dat and input.out with the same
    ! value: temporarily overwrite nx with the global grid size
    nx = n_x_grid
    
    ! write nml
    if(root_processor) write(ifile,NML=gridsize)
    
    nx = n_x_grid / n_procs_x

    ! Do it better for the new metadata: Note that the parameter 'nx'
    ! as output to the metadata is not the value of the input file
    ! (which is the global radial grid size), but the actual local
    ! radial grid size on every process:
    call write_run_parameter('grid', 'nx', nx)
    call write_run_parameter('grid', 'n_x_grid', n_x_grid)
    call write_run_parameter('grid', 'nperiod', nperiod)
    call write_run_parameter('grid', 'n_vpar_grid', n_vpar_grid)
    call write_run_parameter('grid', 'nmod', nmod)
    call write_run_parameter('grid', 'number_of_species', number_of_species)
    call write_run_parameter('grid', 'n_mu_grid', n_mu_grid)
    call write_run_parameter('grid', 'n_trapped', n_trapped)
    call write_run_parameter('grid', 'n_s_grid', n_s_grid)
    call write_run_parameter('grid', 'n_procs_s', n_procs_s)
    call write_run_parameter('grid', 'n_procs_sp', n_procs_sp)
    call write_run_parameter('grid', 'n_procs_vpar', n_procs_vpar)
    call write_run_parameter('grid', 'n_procs_mu', n_procs_mu)
    call write_run_parameter('grid', 'non_blocking_vpar', non_blocking_vpar)
    call write_run_parameter('grid', 'vpmax', vpmax)
    call write_run_parameter('grid', 'mumax', mumax)
    call write_run_parameter('grid', 'psil', psil)
    call write_run_parameter('grid', 'psih', psih)
    call write_run_parameter('grid', 'lx', lx)
    call write_run_parameter('grid', 'n_procs_x', n_procs_x)

  end if

end subroutine grid_read_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> distribute the (global) grid sizes over the different processors
!----------------------------------------------------------------------------
subroutine grid_bcast_nml

  use mpiinterface, only : mpibcast

  call mpibcast(number_of_species,1)
  call mpibcast(n_mu_grid,        1)
  call mpibcast(n_vpar_grid,      1)
  call mpibcast(n_s_grid,         1)
  call mpibcast(n_x_grid,         1)
  call mpibcast(nx,               1)
  call mpibcast(nmod,             1)
  call mpibcast(n_trapped,        1)
  call mpibcast(nperiod,          1)
  call mpibcast(n_procs_s,        1)
  call mpibcast(n_procs_x,        1)
  call mpibcast(n_procs_vpar,     1)
  call mpibcast(n_procs_mu,       1)
  call mpibcast(n_procs_sp,       1)
  call mpibcast(vpmax,            1)
  call mpibcast(mumax,            1)
  call mpibcast(psil,             1)
  call mpibcast(psih,             1)
  call mpibcast(lx,               1)

end subroutine grid_bcast_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> perform checks on the grid params
!----------------------------------------------------------------------------
subroutine grid_check_params

  use general, only : gkw_exit
  use control, only : vp_trap, spectral_radius, flux_tube

  ! NX must always be odd
  if(mod(nx,2).eq.0 .and. spectral_radius) then
     call gkw_exit('grid_size: nx must be odd')
  end if
  
  ! n_x_grid = nx for now (nx may be reset to local value)
  ! do not remove nx from input namelist because many previous
  ! input files that have it. This avoids having both in namelist
  n_x_grid = nx

  ! check we have positive values for the n_procs_*
  if (n_procs_s <= 0 .or. n_procs_mu <= 0 .or. n_procs_sp <= 0 .or.          &
      &                      n_procs_x <= 0 .or. n_procs_vpar <= 0) then
    call gkw_exit('control: '//                                             &
        &          'n_procs_[s|mu|sp|vpar|x] <= 0 not allowed!')
  endif

  ! nperiod:
  if (nperiod <= 0) then
    call gkw_exit('grid_size: unreasonable value of nperiod')
  end if

  ! if vp_trap we require an even number of points in n_vpar_grid
  if (vp_trap == 1 .and. mod(n_vpar_grid,2) /= 0) then
    call gkw_exit('grid_size: '//                                           &
        &          ' for vp_trap we require n_vpar_grid to be even')
  end if

  ! check if the grid sizes are set
  if (number_of_species <= 0) then
    call gkw_exit('grid_size: '//                                           &
        &          'The number of species is zero; set "number_of_species".')
  end if
  if (n_mu_grid <= 0 .or. n_vpar_grid <= 0 .or. n_s_grid <= 0) then
    call gkw_exit('grid_size: '//                                           &
        &          'There are no points in one of n_[s|mu|vpar]_grid; '//    &
        &          'check the input file')
  end if

  if (nmod <= 0 .or. nx <= 0) then
     call gkw_exit('grid_size: Must have at least one mode')
  end if

  if (flux_tube .and. (.not. spectral_radius) .and. lx <= 0) then
     call gkw_exit('lx: radial boxsize must be specified for nonspectral &
        & fluxtube runs')
  end if

  if ((.not.flux_tube) .and. (.not. spectral_radius) .and. (psil <= 0 .and. &
     & psih <= 0)) then
     call gkw_exit('psil, psih: radial boxsize must be specified for nonspectral &
        & global runs')
  end if
end subroutine grid_check_params

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> find out how many processors we shall use in each direction
!----------------------------------------------------------------------------
subroutine parallel_plan

  use mpiinterface, only : number_of_processors
  use general,      only : gkw_warn, gkw_abort
  use global,       only : root_and_verbose

  integer, dimension(nplan) :: itest, procs
  integer :: i
  logical :: done

  ! In order to obtain the maximum number of procs in each direction we
  ! must increase the minimum number of points in that direction until
  ! it divides the total number of points (for efficiency).
  ! In the case where we have less points than the minimum, we set the
  ! minimum to the total number of points, then this will get done on 1
  ! processor in this direction.

  !FJC: this loop could be skipped for single processor jobs
  do i=1, nplan
    if (pp(i)%points < pp(i)%min_points) then
      pp(i)%min_points = pp(i)%points
      if (root_and_verbose) then
        write(*,*) '* reducing the min_points in ',pp(i)%name,' to',&
                  & pp(i)%points
      end if
    else
      addpoints : do
        if ( mod(pp(i)%points,pp(i)%min_points) == 0 ) exit addpoints
        pp(i)%min_points = pp(i)%min_points + 1
        if (root_and_verbose) then
          write(*,*) '* increasing min_points in ',pp(i)%name,' to',&
                    & pp(i)%min_points
        end if
      end do addpoints
    end if
  end do

  ! check various possibilities
  done = .false.

  ! check for single processor -- override any input number of processors
  ! and print a warning message
  if (number_of_processors == 1) then
    done = .true.
    if (any(pp(:)%procs > 1)) then
       call gkw_warn('setup_grid:&
                    & n_procs_<something> has been set > 1 for a single&
                    & processor run')
    end if
    pp(:)%procs = 1
    pp(:)%parallel = .false.
  end if
  
  if (lsendrecv_x) pp(ix)%parallel=.true.

  ! Check if processor numbers have been specified manually (and that they are
  ! consistent with the total number of processors). Avoid the other checks
  ! below if this is the case.
  if (.not. done) then
    do i=1, nplan

      ! check if the number of processors in one direction has been specified
      if ( pp(i)%procs > 1 ) then
        pp(i)%parallel = .true.

        ! check that this number of processors divides the total
        ! number of points
        if (mod(pp(i)%points,pp(i)%procs) == 0) then
          pp(i)%points = pp(i)%points / pp(i)%procs
        else
          call gkw_abort('setup_grid:&
                        & mismatch of points/procs for '//pp(i)%name)
        end if

        ! check if the specified number is too big
        if (pp(i)%points < pp(i)%min_points) then
          call gkw_abort('setup_grid:'//                                     &
                        &' too few points per proc (or too many processors'//&
                        &' specified ) in '//pp(i)%name)
        end if
        done = .true.
      end if

    end do

  end if

  ! For the most efficient parallelisation, the user should explicitly set
  ! and test various combinations. The tests below check through all the
  ! possible combinations in a specific order, which in most cases will
  ! result in a reasonably efficient combination.
  if (.not. done) then
    itest(1)=isp
    call get_parallel_combination(number_of_processors,1,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,1)
  end if
  if (.not. done) then
    itest(1)=imu
    call get_parallel_combination(number_of_processors,1,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,1)
  end if
  if (.not. done) then
    itest(1)=is
    call get_parallel_combination(number_of_processors,1,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,1)
  end if
  if (.not. done) then
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,1,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,1)
  end if

  if (.not. done) then
    itest(2)=isp
    itest(1)=imu
    call get_parallel_combination(number_of_processors,2,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,2)
  end if
  if (.not. done) then
    itest(2)=is
    itest(1)=isp
    call get_parallel_combination(number_of_processors,2,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,2)
  end if
  if (.not. done) then
    itest(2)=imu
    itest(1)=is
    call get_parallel_combination(number_of_processors,2,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,2)
  end if

  if (.not. done) then
    itest(3)=isp
    itest(2)=imu
    itest(1)=is
    call get_parallel_combination(number_of_processors,3,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,3)
  end if

  if (.not. done) then
    itest(2)=isp
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,2,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,2)
  end if
  if (.not. done) then
    itest(2)=imu
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,2,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,2)
  end if

  if (.not. done) then
    itest(2)=is
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,2,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,2)
  end if

  if (.not. done) then
    itest(3)=isp
    itest(2)=imu
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,3,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,3)
  end if

  if (.not. done) then
    itest(4)=is
    itest(3)=isp
    itest(2)=imu
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,4,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,4)
  end if

  if (.not. done) then
    itest(3)=is
    itest(2)=imu
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,3,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,3)
  end if
  if (.not. done) then
    itest(3)=is
    itest(2)=isp
    itest(1)=ivpar
    call get_parallel_combination(number_of_processors,3,itest,procs,done)
    if (done) call set_parallel_plan(itest,procs,3)
  end if

  ! give up if no possible combinations are found
  if (.not. done) call gkw_abort('parallel_plan:'//                          &
                  &' no reasonable processor and gridsize combinations '//   &
                  &' possible for this number of processors')

end subroutine parallel_plan

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine checks if it might be possible to parallelize over the
!> first direction contained in itest. Calling recursively tests all the
!> possibilities for the next directions.
!----------------------------------------------------------------------------
recursive subroutine get_parallel_combination(&
                       & total_procs, n, itest, nprocs, found)

  use general, only : gkw_abort

  !> number of processors available for the decomposition
  integer, intent(in) :: total_procs
  !> number of directions to test
  integer, intent(in) :: n
  !> index of the directions to test
  integer, dimension(nplan), intent(in) :: itest

  !> number of processors in each direction
  integer, dimension(nplan), intent(inout) :: nprocs
  !> denotes combination found successfully
  logical, intent(out) :: found

  !>the maximum number of processors we can use for this combination
  integer :: m_max_procs
  !> maxmimum processors in each direction
  integer, dimension(nplan) :: max_procs

  ! local variables
  integer :: i, j, procs_i, procs_remain
  integer, dimension(nplan) :: itest_next
  integer, dimension(nplan) :: nprocs_next
  logical :: found_next, good
  integer :: s_vpar_count,vpar_mu_count

  ! defaults
  found = .false.
  found_next = .false.

  !if (root_and_verbose) write (*,*) '* itest: n=',n
  !if (root_and_verbose) write (*,*) '* testing: ',pp(itest(1:n))%name

  s_vpar_count = 0 ; vpar_mu_count = 0
  ! checks - probably not necessary (could just return?)
  do i=1, n
    if (itest(i) == is) s_vpar_count = s_vpar_count + 1
    if (itest(i) == ivpar ) then
      s_vpar_count = s_vpar_count + 1
      vpar_mu_count = vpar_mu_count + 1
    end if
    if (itest(i) == imu ) vpar_mu_count = vpar_mu_count + 1
    do j=1, n
      if (i .ne. j) then
        if (itest(i) .eq. itest(j)) then
          call gkw_abort('get_parallel_combination: bad input')
        end if
      end if
    end do
  end do

  !if (lcollisions .and. vpar_mu_count >= 2) then
  !  if (root_and_verbose) write (*,*) '* skipping both ivpar and imu with collisions'
  !  return
  !end if
  !if (ltrapping_arakawa .and. s_vpar_count >= 2) then
  !  if (root_and_verbose) write (*,*) '* skipping both ivpar and is for Arakawa'
  !  return
  !end if

  ! set the max procs in each direction
  do i=1, n
    max_procs(i) = pp(itest(i))%points / pp(itest(i))%min_points
  end do

  ! set the maximum number of procs we can possibly use
  m_max_procs = 1
  do i=1, n
    m_max_procs = m_max_procs * max_procs(i)
  end do

  ! check if we have too many processors
  if (total_procs > m_max_procs) then
    return
  end if

  ! loop over the possible number of processors in the 1st direction
  loop_over_dim_procs : do procs_i = 2, max_procs(1)

    good = .false.

    ! check if this number of processors divides the total
    good = mod(total_procs,procs_i) == 0

    ! check that there would be the same number of points per processor in
    ! that direction
    good = good .and. ( mod(pp(itest(1))%points,procs_i) == 0)

    this_one_is_good : if (good) then

      ! this value works, so either stop if we are at the last level or call
      ! the routine again to check the next dimension
      procs_remain = total_procs / procs_i

      ! check that the number of points per proc would not be greater than the
      ! maximum number of points in a parallel scheme (e.g. 2 for the
      ! fourth-order scheme with vpar).
      if ((pp(itest(1))%points / procs_i) <= pp(itest(1))%max_points) &
        & then

        ! terminate the recursion if we get to the end and pass back found or
        ! found_next together with the number of procs
        if (n == 1 .and. procs_remain == 1) then
           nprocs(1) = procs_i
           found = .true.
           return
        end if

        ! check the next direction
        if (n > 1) then
          do j=1, n-1
            itest_next(j) = itest(j+1)
            nprocs_next(j)= nprocs(j+1)
          end do
          call get_parallel_combination(&
              & procs_remain, n-1, itest_next, nprocs_next, found_next)
        end if

        ! If the next direction works, then all directions must work so unpack
        ! the number of processors for that direction into nprocs. The next
        ! level up will do the same again until it returns to the original
        ! calling routine.
        if (found_next) then
          nprocs(1) = procs_i
          do j=1, n-1
            nprocs(j+1) = nprocs_next(j)
          end do
          found = .true.
          return
        end if

      end if

    end if this_one_is_good

  end do loop_over_dim_procs

  ! We have already looped over all possibilities for the first direction,
  ! so there is no match.
  found = .false.
  return

end subroutine get_parallel_combination

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine fills in the parallel plan based on the number of
!> processors in each direction in iplan
!----------------------------------------------------------------------------
subroutine set_parallel_plan(itest, procs, n_parallel)

  !> array holding the indicies of the things to parallelize over
  integer, dimension(nplan), intent(in) :: itest
  !> the number of procs for the parallel directions
  integer, dimension(nplan), intent(in) :: procs
  !> the number of elements in itest over which we shall parallelize
  integer, intent(in) :: n_parallel

  integer :: i, j
  logical :: fill

  ! for all parallel directions, set the grid sizes and number of processors
  do i=1, n_parallel
    pp(itest(i))%parallel = .true.
    pp(itest(i))%procs    = procs(i)
    pp(itest(i))%points   = pp(itest(i))%points / pp(itest(i))%procs
  end do

  ! for all other directions, set the number of processors to 1
  do i=1, nplan
    fill = .true.
    do j=1, n_parallel
      if (i .eq. itest(j)) fill = .false.
    enddo
    if (fill) then
      pp(i)%parallel    = .false.
      pp(i)%procs       = 1
    end if
  end do

end subroutine set_parallel_plan

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine checks that the parallel layout is consistent for the
!> variables provided to the rest of the code
!----------------------------------------------------------------------------
subroutine check_parallel_layout

  use general,      only : gkw_abort
  use mpiinterface, only : number_of_processors
  use control,      only : vp_trap, spectral_radius

  ! -- not implemented or tested--
  if (parallel_vpar) then
    if (vp_trap == 1) then
      call gkw_abort('check_parallel_layout:'//                              &
          & 'n_procs_vpar > 1 not tested with vp_trap=1?')
    end if
  end if !parallel v_par

  if (parallel_s .and. vp_trap == 1) call gkw_abort('vp_trap not implemented for parallel_s')
  
  if (parallel_x .and. spectral_radius) then
     call gkw_abort('Parallel x is not possible with spectral_radius')
  end if !parallel_x

  ! -- consistency checks --
  ! total number of procs selected
  if(n_procs_mu*n_procs_s*n_procs_sp*n_procs_vpar*n_procs_x /=               &
                                   & number_of_processors) then
    write(*,*) 'mu',n_procs_mu
    write(*,*) 'sp',n_procs_sp
    write(*,*) ' s',n_procs_s
    write(*,*) ' x',n_procs_x
    write(*,*) 'vp',n_procs_vpar
    write(*,*) 'nn',number_of_processors
    call gkw_abort('check_parallel_layout:'//                                &
                  &' n_procs_mu*n_procs_s*n_procs_sp*n_procs_vpar /='//      &
                  &' number_of_procs !')
  end if

  ! too many procs in particular direction
  if (n_procs_mu > n_mu_grid) then
    call gkw_abort('check_parallel_layout: too many procs for mu_grid')
  end if
  if (n_procs_sp > number_of_species) then
    call gkw_abort('check_parallel_layout:'//                                &
        &          ' too many procs for number of species')
  end if
  if (n_procs_s > n_s_grid) then
    call gkw_abort('check_parallel_layout: too many procs in s')
  end if
  if (n_procs_x > n_x_grid) then
    call gkw_abort('check_parallel_layout: too many procs in x')
  end if
  if (n_procs_vpar > n_vpar_grid) then
    call gkw_abort('check_parallel_layout: too many procs in vpar')
  end if

  ! -- efficiency checks --
  ! all procs should have the same amount of work so quit if that is not true
  if (mod(n_mu_grid, n_procs_mu) /= 0) then
    call gkw_abort('check_parallel_layout:'//                                &
                  &' n_mu_grid is not a multiple of n_procs_mu')
  end if
  if (mod(n_s_grid, n_procs_s) /= 0)  then
    call gkw_abort('check_parallel_layout:'//                                &
                  &' ns is not a multiple of n_procs_s')
  end if
  if (mod(n_x_grid, n_procs_x) /= 0)  then
    call gkw_abort('check_parallel_layout:'//                                &
                  &' ns is not a multiple of n_procs_x')
  end if
  if (mod(number_of_species, n_procs_sp) /= 0) then
    call gkw_abort('check_parallel_layout:'//                                &
                  &' number_of_species is not a multiple of n_procs_sp')
  end if
  if (mod(n_vpar_grid, n_procs_vpar) /= 0) then
    call gkw_abort('check_parallel_layout:'//                                &
                  &' n_vpar_grid is not a multiple of n_procs_vpar')
  end if

end subroutine check_parallel_layout

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> print the parallel plan to stdout
!----------------------------------------------------------------------------
subroutine parallel_report(a)

  use mpiinterface, only : root_processor, number_of_processors
  use control,      only : vp_trap
  use global,       only : root_and_verbose

  type (plan), intent(in), dimension(nplan) :: a

  integer :: i

  if (root_processor) then
    write (*,*) ' --- parallel report --- '
    write (*,*)
    write (*,*) 'given:'
    write (*,*) ' number of processors =',number_of_processors
    write (*,*)
    write (*,*) ' number of species    =',number_of_species
    write (*,*) ' n_mu_grid            =',n_mu_grid
    write (*,*) ' n_vpar_grid          =',n_vpar_grid
    write (*,*) ' n_s_grid             =',n_s_grid
    write (*,*) ' n_x_grid             =',n_x_grid
    write (*,*)
    write (*,*) 'The code will:'
    write (*,*)

    do i=1, nplan

     if(a(i)%parallel) then
       write (*,*) 'parallelise over ',a(i)%name
       write (*,*) '   with'
       write (*,*) 'procs           =',a(i)%procs
       write (*,*) 'points per proc =',a(i)%points
       write (*,*)
     end if

    end do

    if (root_and_verbose) then
      write (*,*) ' per processor values:'
      write (*,*) '   nsp  =',nsp
      write (*,*) '   nmu  =',nmu
      write (*,*) '  nvpar =',nvpar
      write (*,*) '    ns  =',ns
      write (*,*) '    nx  =',nx
      write (*,*) '         '
      write (*,*) ' arrays for the first and last points on the root processor'
      write (*,*) ' isppb     =',isppb
      write (*,*) ' isppe     =',isppe
      write (*,*) ' ixppb     =',ixpb
      write (*,*) ' ixppe     =',ixpe
      write (*,*) ' isppb     =',ispb
      write (*,*) ' isppe     =',ispe
      write (*,*) ' imupb     =',imupb
      write (*,*) ' imupe     =',imupe
      write (*,*) ' ivparpb   =',ivparpb
      write (*,*) ' ivparpe   =',ivparpe
    end if

    if (vp_trap .eq. 1) then
      write (*,*) ' (using vp_trap)'
    end if

    write (*,*)
    write (*,*) ' --- end of parallel report ---'

  end if

end subroutine parallel_report


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Set up additional communicators for various hyperslabs in the cartesian
!> topology. When using MPI, any additional communicator should be correctly
!> initialised, irrespective of if it will be used or not.
!----------------------------------------------------------------------------
subroutine setup_communicators

  use mpiinterface, only : MPI_COMM_SELF, number_of_processors
  use mpicomms,     only : COMM_ALL_EQ, COMM_VPAR_NE, COMM_VPAR_EQ 
  use mpicomms,     only : COMM_CART, COMM_S_NE, COMM_S_EQ, COMM_X_NE
  use mpicomms,     only : COMM_X_EQ, COMM_MU_NE, COMM_MU_EQ, COMM_SP_NE
  use mpicomms,     only : COMM_SP_EQ, COMM_SP_EQ_S_EQ, COMM_SP_EQ_X_EQ
  use mpicomms,     only : COMM_S_EQ_X_EQ, COMM_VPAR_NE_MU_NE
  use mpicomms,     only : COMM_VPAR_EQ_MU_EQ, COMM_VPAR_EQ_X_EQ
  use mpicomms,     only : COMM_SP_NE_MU_NE, COMM_S_NE_X_NE, COMM_S_NE_MU_NE
  use mpicomms,     only : COMM_SP_NE_X_NE, COMM_SP_NE_S_NE
  use mpicomms,     only : COMM_ROOT

  logical, dimension (nplan) :: remain_dims
  integer :: comm_cart_group, only_root_group
  
  ! ! debugging:
  !integer :: rank

  ! Default values; these assume there is no parallelisation in a particular
  ! direction e.g. COMM_VPAR_NE should contain only the calling processor by
  ! default, as all the vpar points are local (there are no other vpar
  ! points); COMM_VPAR_EQ should contain all processors as they all solve for
  ! the same part of the vpar grid. MPI_COMM_SELF and COMM_CART are suitable
  ! communicators for these examples; these are provided by either MPI or a
  ! dummy value in mpiinterface.
  COMM_ALL_EQ     = MPI_COMM_SELF
  COMM_VPAR_NE    = MPI_COMM_SELF
  COMM_VPAR_EQ    = COMM_CART
  COMM_S_NE       = MPI_COMM_SELF
  COMM_S_EQ       = COMM_CART
  COMM_X_NE       = MPI_COMM_SELF
  COMM_X_EQ       = COMM_CART
  COMM_MU_NE      = MPI_COMM_SELF
  COMM_MU_EQ      = COMM_CART
  COMM_SP_EQ      = COMM_CART
  COMM_SP_NE      = MPI_COMM_SELF
  COMM_SP_EQ_S_EQ = COMM_CART
  COMM_SP_EQ_X_EQ = COMM_CART
  COMM_S_EQ_X_EQ  = COMM_CART
  COMM_VPAR_NE_MU_NE = MPI_COMM_SELF
  COMM_VPAR_EQ_MU_EQ = COMM_CART
  COMM_VPAR_EQ_X_EQ  = COMM_CART
  COMM_SP_NE_MU_NE   = MPI_COMM_SELF
  COMM_S_NE_X_NE     = MPI_COMM_SELF
  COMM_S_NE_MU_NE     = MPI_COMM_SELF
  COMM_SP_NE_X_NE     = MPI_COMM_SELF
  COMM_SP_NE_S_NE     = MPI_COMM_SELF

#if defined(mpi2)
  call mpi_comm_group(COMM_CART, comm_cart_group, ierr)
  call mpi_group_incl(comm_cart_group, 1, (/0/), only_root_group, ierr)
  call mpi_comm_create(COMM_CART, only_root_group, COMM_ROOT, ierr)
#endif
  
  ! Sub-communicators for processors parallel (NE) and perpendicular (EQ) in
  ! the processor grid.
  remain_dims(:) = .false.
#if defined(mpi2)
  if (number_of_processors > 1 .or. lsendrecv_x) then
    ! vpar
    remain_dims(:) = .false.
    if (pp(ivpar)%parallel) remain_dims(pp(ivpar)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_VPAR_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_VPAR_NE)
    remain_dims(:) = (.not. remain_dims(:))
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_VPAR_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_VPAR_EQ)

    ! s
    remain_dims(:) = .false.
    if (pp(is)%parallel) remain_dims(pp(is)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_S_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_S_NE)
    remain_dims(:) = (.not. remain_dims(:))
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_S_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_S_EQ)
    
    ! x
    remain_dims(:) = .false.
    if (pp(ix)%parallel) remain_dims(pp(ix)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_X_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_X_NE)
    remain_dims(:) = (.not. remain_dims(:))
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_X_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_X_EQ)

    ! mu
    remain_dims(:) = .false.
    if (pp(imu)%parallel) remain_dims(pp(imu)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_MU_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_MU_NE)
    remain_dims(:) = (.not. remain_dims(:))
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_MU_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_MU_EQ)

    ! species
    remain_dims(:) = .false.
    if (pp(isp)%parallel) remain_dims(pp(isp)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_SP_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_SP_NE)
    remain_dims(:) = (.not. remain_dims(:))
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_SP_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_SP_EQ)

    ! not equal mu and vpar
    remain_dims(:) = .false.
    if (pp(ivpar)%parallel) remain_dims(pp(ivpar)%direction) = .true.
    if (pp(imu)%parallel)   remain_dims(pp(imu)%direction)   = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_VPAR_NE_MU_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_VPAR_NE_MU_NE)

    ! equal mu and vpar
    remain_dims(:) = .true.
    if (pp(imu)%parallel)   remain_dims(pp(imu)%direction)    = .false.
    if (pp(ivpar)%parallel) remain_dims(pp(ivpar)%direction)  = .false.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_VPAR_EQ_MU_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_VPAR_EQ_MU_EQ)

    ! equal x and vpar
    remain_dims(:) = .true.
    if (pp(ix)%parallel)   remain_dims(pp(ix)%direction)    = .false.
    if (pp(ivpar)%parallel) remain_dims(pp(ivpar)%direction)  = .false.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_VPAR_EQ_X_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_VPAR_EQ_X_EQ)
    
    ! not equal species and mu
    remain_dims(:) = .false.
    if (pp(isp)%parallel) remain_dims(pp(isp)%direction) = .true.
    if (pp(imu)%parallel)   remain_dims(pp(imu)%direction)   = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_SP_NE_MU_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_SP_NE_MU_NE)

    ! equal s and species
    remain_dims(:) = .true.
    if (pp(isp)%parallel) remain_dims(pp(isp)%direction) = .false.
    if (pp(is)%parallel)  remain_dims(pp(is)%direction)  = .false.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_SP_EQ_S_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_SP_EQ_S_EQ)

    ! equal x and species
    remain_dims(:) = .true.
    if (pp(ix)%parallel) remain_dims(pp(ix)%direction) = .false.
    if (pp(isp)%parallel) remain_dims(pp(isp)%direction) = .false.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_SP_EQ_X_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_SP_EQ_X_EQ)
    
    ! equal s and x
    remain_dims(:) = .true.
    if (pp(is)%parallel) remain_dims(pp(is)%direction) = .false.
    if (pp(ix)%parallel)  remain_dims(pp(ix)%direction)  = .false.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_S_EQ_X_EQ,ierr)
    call mpicartsubfix(remain_dims,COMM_S_EQ_X_EQ)

    ! not equal s and x
    remain_dims(:) = .false.
    if (pp(is)%parallel) remain_dims(pp(is)%direction) = .true.
    if (pp(ix)%parallel) remain_dims(pp(ix)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_S_NE_X_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_S_NE_X_NE)

    ! not equal s and mu
    remain_dims(:) = .false.
    if (pp(is)%parallel) remain_dims(pp(is)%direction) = .true.
    if (pp(imu)%parallel) remain_dims(pp(imu)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_S_NE_MU_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_S_NE_MU_NE)

    ! ! debugging
    ! call mpi_comm_rank(COMM_S_NE_X_NE,rank,ierr)
    ! call MPI_CART_COORDS(COMM_S_NE_X_NE,rank,nplan,coords,ierr)
    ! write (*,*) "COMM_S_NE_X_NE rank=",rank,"coords=",coords

    ! not equal species and x
    remain_dims(:) = .false.
    if (pp(isp)%parallel) remain_dims(pp(isp)%direction) = .true.
    if (pp(ix)%parallel) remain_dims(pp(ix)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_SP_NE_X_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_SP_NE_X_NE)

    ! not equal species and s
    remain_dims(:) = .false.
    if (pp(isp)%parallel) remain_dims(pp(isp)%direction) = .true.
    if (pp(is)%parallel) remain_dims(pp(is)%direction) = .true.
    call MPI_CART_SUB(COMM_CART,remain_dims(1:n_dims),COMM_SP_NE_S_NE,ierr)
    call mpicartsubfix(remain_dims,COMM_SP_NE_S_NE)

  end if
#endif

end subroutine setup_communicators

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Use COMM_S_EQ_X_EQ (or subset) to distribute fields solve work over 
!> the toroidal modes and setup mpi datatypes for gathering them again.
!> This is currently only used for nonspectral fields.
!> TO DO: Fix to be independent of index function order
!> TO DO: More flexible choice of communicator subset using mpi_comm_split
!----------------------------------------------------------------------------
subroutine mpisetup_mod(first_imod,last_imod,nmod_l,type_send,type_recv,enable)

  use mpiinterface, only : MPI_COMM_SELF, mpi_address_kind, MPICOMPLEX_X
  use mpiinterface, only : number_of_processors, root_processor 
  use mpiinterface, only : gather_array
  use mpicomms,     only : COMM_MOD, COMM_S_EQ_X_EQ, COMM_SP_NE_MU_NE
  use mpicomms,     only : COMM_MU_NE
  use control,      only : nlapar
  use general,      only : gkw_abort
  
  !> first, last and number of modes to do on this proc
  integer, intent(out) :: first_imod, last_imod, nmod_l
  integer, intent(out) :: type_send, type_recv
  logical, intent(in), optional :: enable ! set F to disable the parallelism
  integer :: rank_cart, ierr, n_procs_mod, nmod_pp, dummax, i, size
  integer :: type_dum2, type_dum, fields
  integer, allocatable :: blocklens(:), disps(:)

#if defined(mpi2)  
  integer(kind=mpi_address_kind) :: start, extent
#endif

  ! First choose a communicator of a sensible size
  ! For maximum work distribution over the chosen communicator.
  ! The mu and sp communicators are preferred due to their 
  ! ordering in the cartesian topology.
  ! (The choices made here to some extent reflect typical grid sizes,
  ! there are probably more elegant general solutions which mpi_comm_split
  ! the most local possible communicator into the ideal size)
  
  ! This assumes the gain from work sharing outweighs communication overheads. 
  ! Timing tests have indicated the gather is typically quite efficient.
  ! A smaller nmod_pp also minimises the communications in the x_gather.

  n_procs_mod = n_procs_vpar * n_procs_mu * n_procs_sp
  COMM_MOD = COMM_S_EQ_X_EQ  
      
  if (n_procs_mod > nmod) then  
    n_procs_mod =  n_procs_mu * n_procs_sp
    COMM_MOD = COMM_SP_NE_MU_NE
  end if  

  if (n_procs_mod > nmod) then   
    n_procs_mod =  n_procs_mu
    COMM_MOD = COMM_MU_NE
  end if 
     
  ! Other options that may be desireable for atypical grids / decompositions
  ! n_procs_mod = n_procs_vpar * n_procs_mu 
  ! COMM_MOD = COMM_VPAR_NE_MU_NE
   
  ! n_procs_mod =  n_procs_vpar
  ! COMM_MOD = COMM_VPAR_NE
  
  ! n_procs_mod =  n_procs_sp
  ! COMM_MOD = COMM_SP_NE
  
  if (nmod == 1) then
    n_procs_mod = 1
    COMM_MOD = MPI_COMM_SELF
  end if

  first_imod = 0
  last_imod = 0

#if defined(mpi2)
  !Get my rank in COMM_MOD
  call MPI_COMM_RANK(COMM_MOD,rank_cart,ierr)
#else
  rank_cart = 0
#endif
  
  !Some checks
  if (rank_cart < 0 .or. rank_cart > n_procs_mod - 1) then
    call gkw_abort('You have not understood COMM_S_EQ_X_EQ rank')
  end if
  
  !Split the work as evenly as possible
  nmod_pp = ceiling(real(nmod)/real(n_procs_mod))
  first_imod = rank_cart * nmod_pp + 1
  dummax = (rank_cart + 1) * nmod_pp
  if (rank_cart == n_procs_mod -1 .and. nmod > dummax) then
    call gkw_abort('Error in mpisetup_mod')
  end if
  last_imod  = min(dummax,nmod)
  
  ! final checks
  if (first_imod <= 0 .or. last_imod == 0 .or. last_imod > nmod) then
    call gkw_abort('Error in mpisetup_mod')
  end if
  
  nmod_l = last_imod - first_imod + 1
  if (nmod_l < 1) nmod_l = 0
  
  ! create a datatype corresponding to the data received for a single mode
  ! DEPENDS ON INDEX ORDER
  ! The datatype setup is explained very nicely here:
  ! http://stackoverflow.com/questions/5530648/  
  ! http://stackoverflow.com/questions/6511356/  
#if defined(mpi2)
  fields = 1; if (nlapar) fields = 2
  allocate(blocklens(fields*nx*ns))
  allocate(disps(fields*nx*ns)) 
  
  blocklens(:)=1
  disps(1) = 0
  do i = 2, fields*nx*ns
    disps(i) = disps(i-1) + nmod
  end do
  ! Might create more efficent datastructure if used mpi_type_vector ?
  call MPI_TYPE_INDEXED(fields*nx*ns, blocklens, disps, MPICOMPLEX_X, type_dum, ierr)
  call MPI_TYPE_size(MPICOMPLEX_X,size,ierr)        
  ! redefine the upper and lower bounds for receiving multiple modes
  start = 0
  extent = size
  call MPI_TYPE_create_resized(type_dum, start, extent, type_recv, ierr);  
  call MPI_TYPE_COMMIT(type_recv,ierr)

  !create a datatype corresponding to the data sent for a single mode  
  disps(1) = 0
  do i = 2, fields*nx*ns
    disps(i) = disps(i-1) + nmod_l
  end do
  call MPI_TYPE_INDEXED(fields*nx*ns, blocklens, disps, MPICOMPLEX_X, type_dum2, ierr)
  ! redefine the upper and lower bounds for sending multiple modes
  call MPI_TYPE_create_resized(type_dum2, start, extent, type_send, ierr);
  call MPI_TYPE_COMMIT(type_send,ierr)      
#else
  type_recv = MPICOMPLEX_X
  type_send = MPICOMPLEX_X
#endif

  !setup arrays for gatherv
  allocate(gathv_mod_recvno(n_procs_mod),stat=ierr)
  allocate(gathv_mod_disps(n_procs_mod),stat=ierr)
  
  ! use as dummy array of size 1
  gathv_mod_disps(1)=nmod_l
    
  call gather_array(gathv_mod_recvno, n_procs_mod, gathv_mod_disps(1:1), 1, &
                    & COMM_MOD, ALLGATHER = .true.)
  
  gathv_mod_disps(1) = 0
  do i = 2,n_procs_mod
    gathv_mod_disps(i) = gathv_mod_disps(i-1) + nmod_pp
  end do
  
  ! TO DISABLE PARALLEL TOROIDAL and make the gatherv a simple copy
  if (present(enable)) then 
    if (.not.enable) then 
      first_imod = 1
      last_imod = nmod
      nmod_l = nmod
      n_procs_mod = 1
      nmod_pp = nmod
      COMM_MOD = MPI_COMM_SELF
      type_recv = type_send
      gathv_mod_disps(:) = 0
      gathv_mod_recvno(:) = nmod   
    end if
  end if 
  
  !write(*,*) 'rank', rank_cart, 'first_imod', first_imod, 'last_imod', last_imod, 'nmod_l', nmod_l
  if (nmod > 1 .and. root_processor .and. number_of_processors > 1) then
    write(*,*) 
    write(*,*) 'Field solve modes per proc:', nmod_pp
    write(*,*) 'Size of modes communicator:', n_procs_mod 
    write(*,*)
  end if
  
  !if (root_processor) write(*,*) gathv_mod_recvno 
  !if (root_processor) write(*,*) gathv_mod_disps

end subroutine mpisetup_mod


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Possible fix for MPI_CART_SUB problem on some systems where a sub
!> cartesian communicator containing none of the original dimensions is
!> assigned MPI_COMM_NULL rather than MPI_COMM_SELF. In the first observed
!> instance of this behaviour, MPI_COMM_SELF was given to the root processor
!> but MPI_COMM_NULL was given to the others.
!----------------------------------------------------------------------------
subroutine mpicartsubfix(rdims,comm)

  use mpiinterface, only : MPI_COMM_SELF

  logical, dimension(:), intent(in) :: rdims
  integer, intent(inout) :: comm

  if (count(rdims(1:n_dims)) == 0) comm = MPI_COMM_SELF

end subroutine mpicartsubfix


!-----------------------------------------------------------------------------
!> Work out the number of processors associated with some sub-communicators
!> (*_nprocsets) and give the communicators a label (*_procset) in the range
!> 1 <= procset <= nprocsets.
!-----------------------------------------------------------------------------
subroutine set_comm_sizes()

  use mpicomms, only : COMM_SP_EQ_nprocsets, COMM_SP_EQ_procset
  use mpicomms, only : COMM_S_EQ_nprocsets, COMM_S_EQ_procset
  use mpicomms, only : COMM_X_EQ_nprocsets, COMM_X_EQ_procset
  use mpicomms, only : COMM_SP_EQ_S_EQ_nprocsets, COMM_SP_EQ_S_EQ_procset

  COMM_SP_EQ_nprocsets = n_procs_sp
  COMM_SP_EQ_procset   = pp(isp)%iproc + 1

  COMM_S_EQ_nprocsets = n_procs_s
  COMM_S_EQ_procset   = pp(is)%iproc + 1
  
  COMM_X_EQ_nprocsets = n_procs_x
  COMM_X_EQ_procset   = pp(ix)%iproc + 1

  COMM_SP_EQ_S_EQ_nprocsets = n_procs_s*n_procs_sp
  COMM_SP_EQ_S_EQ_procset   = pp(is)%iproc + (pp(isp)%iproc*n_procs_s) + 1

end subroutine set_comm_sizes

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Provide the ranks of neighbouring processors for directions in which
!> send/recv communications are required. This requires the additional
!> communicators of mpicomms to be initialised.
!----------------------------------------------------------------------------
subroutine get_local_ranks

  use mpiinterface, only : number_of_processors, MPI_PROC_NULL
  use mpicomms,     only : COMM_MU_NE, COMM_S_NE, COMM_X_NE, COMM_VPAR_NE
  use general,      only : gkw_abort
  
  integer :: my_rank, my_coord(1), opp_coord(1)

  ! Default values; MPI_PROC_NULL - a dummy processor
  proc_vpar_prev          = MPI_PROC_NULL
  proc_vpar_next          = MPI_PROC_NULL
  proc_s_prev             = MPI_PROC_NULL
  proc_s_next             = MPI_PROC_NULL
  proc_x_prev             = MPI_PROC_NULL
  proc_x_next             = MPI_PROC_NULL
  proc_mu_prev            = MPI_PROC_NULL
  proc_mu_next            = MPI_PROC_NULL
  proc_vpar_prev_mu_prev  = MPI_PROC_NULL
  proc_vpar_prev_mu_next  = MPI_PROC_NULL
  proc_vpar_next_mu_prev  = MPI_PROC_NULL
  proc_vpar_next_mu_next  = MPI_PROC_NULL
  proc_vpar_prev_s_prev   = MPI_PROC_NULL
  proc_vpar_prev_s_next   = MPI_PROC_NULL
  proc_vpar_next_s_prev   = MPI_PROC_NULL
  proc_vpar_next_s_next   = MPI_PROC_NULL
  proc_vpar_opposite      = MPI_PROC_NULL

#ifdef mpi2
  ! next and previous processor for mu, s, x and vpar
  if (pp(imu)%parallel) then
    call MPI_CART_SHIFT(COMM_MU_NE,0,1,proc_mu_prev,proc_mu_next,ierr)
  end if
  ! If this was always made true, even not with parallel s 
  ! Then many of the if statements could be removed in connect_parallel
  if (pp(is)%parallel)  then
    call MPI_CART_SHIFT(COMM_S_NE,0,1,proc_s_prev,proc_s_next,ierr)
  end if
  if (pp(ix)%parallel)  then
    call MPI_CART_SHIFT(COMM_X_NE,0,1,proc_x_prev,proc_x_next,ierr)
  end if  
  if (pp(ivpar)%parallel) then
    call MPI_CART_SHIFT(COMM_VPAR_NE,0,1,proc_vpar_prev,proc_vpar_next,ierr) 
    
    ! Only needed for krook operator option 1
    ! get my rank in COMM_VPAR_NE (0 based)
    call MPI_COMM_RANK(COMM_VPAR_NE,my_rank,ierr)   
    ! Needed in unlikely case that ranks /= coords in 1D communicator
    call MPI_CART_COORDS(COMM_VPAR_NE,my_rank,1,my_coord,ierr) 
    ! get opposite coord in COMM_VPAR_NE (0 based)
    opp_coord = pp(ivpar)%procs - my_coord - 1
    if (opp_coord(1) < 0 .or. opp_coord(1) > pp(ivpar)%procs - 1) then 
      call gkw_abort('vpar unexpected opp_coord')
    end if
    call MPI_CART_RANK(COMM_VPAR_NE,opp_coord,proc_vpar_opposite,ierr)    
    !write(*,*) my_rank, my_coord, opp_coord, proc_vpar_opposite
    
  end if
#endif

  ! The processors along diagonals in the grid can always be found by
  ! systematically going through the processors till the one responsible for
  ! the right part of the grid is found. However, we do something a little
  ! more complicated here...

  if (number_of_processors > 1) then
    proc_vpar_next_mu_next = cart_rank(dim1='vpar',s1=+1,dim2='mu',s2=+1)
    proc_vpar_next_mu_prev = cart_rank(dim1='vpar',s1=+1,dim2='mu',s2=-1)
    proc_vpar_prev_mu_next = cart_rank(dim1='vpar',s1=-1,dim2='mu',s2=+1)
    proc_vpar_prev_mu_prev = cart_rank(dim1='vpar',s1=-1,dim2='mu',s2=-1)
    proc_vpar_next_s_next = cart_rank(dim1='vpar',s1=+1,dim2='s',s2=+1)
    proc_vpar_next_s_prev = cart_rank(dim1='vpar',s1=+1,dim2='s',s2=-1)
    proc_vpar_prev_s_next = cart_rank(dim1='vpar',s1=-1,dim2='s',s2=+1)
    proc_vpar_prev_s_prev = cart_rank(dim1='vpar',s1=-1,dim2='s',s2=-1)
  end if

end subroutine get_local_ranks


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> find the rank in COMM_CART of the proc shifted in different directions
!----------------------------------------------------------------------------
function cart_rank(dim1,s1,dim2,s2,dim3,s3,dim4,s4,dim5,s5)

  use mpiinterface, only : MPI_PROC_NULL
#if defined(mpi2)
  use mpicomms,     only : COMM_CART
  
  integer :: ierr 

#endif 

  character (len=*), intent(in) :: dim1,dim2,dim3,dim4,dim5
  integer, intent(in) :: s1,s2,s3,s4,s5
  optional :: dim2,s2,dim3,s3,dim4,s4,dim5,s5
  integer :: cart_rank

  integer, dimension(nplan) :: coords_in
  logical :: bad_shift, bad_shift_test
  integer :: j

  coords_in(:) = coords(:)
  j=1

  bad_shift_test = .false.
  bad_shift = .false.
  !
  ! The plan dims do not correspond to the parallel dims so each must
  ! be checked in sequence. Periodic shifts are wrapped around.

  check_shift : do i=1,nplan

      if (pp(i)%parallel) then

      if (trim(pp(i)%name)==trim(dim1)) call test_periodic_shift(i,j,s1,bad_shift_test)
      bad_shift = bad_shift .or. bad_shift_test
      if (present(dim2) .and. present(s2)) then
        if (trim(pp(i)%name)==trim(dim2)) call test_periodic_shift(i,j,s2,bad_shift_test)
        bad_shift = bad_shift .or. bad_shift_test
      end if
      if (present(dim3) .and. present(s3)) then
        if (pp(i)%name==dim3) call test_periodic_shift(i,j,s3,bad_shift_test)
        bad_shift = bad_shift .or. bad_shift_test
      end if
      if (present(dim4) .and. present(s4)) then
        if (pp(i)%name==dim4) call test_periodic_shift(i,j,s4,bad_shift_test)
        bad_shift = bad_shift .or. bad_shift_test
      end if
      if (present(dim5) .and. present(s5)) then
        if (pp(i)%name==dim5) call test_periodic_shift(i,j,s5,bad_shift_test)
        bad_shift = bad_shift .or. bad_shift_test
      end if

      if (bad_shift) exit check_shift

      j=j+1

    end if
  enddo check_shift

  ! default
  cart_rank = MPI_PROC_NULL

  ! This avoids the error of attempting to obtain the rank if the coords are
  ! not in the processor grid.
#ifdef mpi2
  if (.not. bad_shift) call MPI_CART_RANK(COMM_CART,coords_in,cart_rank,ierr)
#endif

  contains

  !*************************************************************************
  !> If the direction is periodic, shift the coordinate appropriately.
  !> Otherwise, note a bad shift so that we get a null proc at the end.
  !-------------------------------------------------------------------------
  subroutine test_periodic_shift(i2,j2,shift,bad)

    logical, intent(out) :: bad
    integer, intent(in) :: i2, j2, shift

    coords_in(j2) = coords_in(j2) + shift
    bad = .false.

    if (coords_in(j2) < 0) then
      if (pp(i2)%periodic) then
        coords_in(j2) = coords_in(j2) + pp(i2)%procs
        bad = .false.
      else
        bad = .true.
      end if
    end if

    if (coords_in(j2) >= pp(i2)%procs) then
      if (pp(i2)%periodic) then
        coords_in(j2) = coords_in(j2) - pp(i2)%procs
        bad = .false.
      else
        bad = .true.
      end if
    end if

  end subroutine test_periodic_shift
  !*************************************************************************

end function cart_rank


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Create a cartesian processor topology (communicator COMM_CART) and get the
!> coordinates of the local processor within the processor grid (coords).
!----------------------------------------------------------------------------
subroutine setup_proc_coords

  use general,      only : gkw_abort
  use mpicomms,     only : COMM_CART
#if defined(mpi2)
  use mpiinterface, only : number_of_processors, mpicart_create
  use mpiinterface, only : mpicopy_world

  integer :: i,j,parallel_dim_count, rank_cart
  
#else 
  use global,       only : i_huge

  integer :: i,j,parallel_dim_count
#endif 


  ! state which directions should be periodic or not
  pp(:)%periodic      =   .false.
  pp(is)%periodic     =   .true.

  ! For flux tube only !!
  pp(ix)%periodic     =   .true.

  ! find out how many dimensions the cartesian domain has
  n_dims = count(pp(:)%parallel)

  ! consistency check counter
  parallel_dim_count = 0

  do i=1, nplan
    if (pp(i)%parallel) then
      ! MPI direction indexing starts at 1
      parallel_dim_count = parallel_dim_count + 1

      pp(i)%direction = parallel_dim_count
      dims(parallel_dim_count)     = pp(i)%procs
      periodic(parallel_dim_count) = pp(i)%periodic
    end if
  end do

  if (parallel_dim_count /= n_dims) then
    call gkw_abort('setup_proc_coords: bad parallel_dim_count')
  end if

  coords = 0

#if defined(mpi2)
  ! create the topology; give the communicator a value in all cases
  if (number_of_processors > 1 .or. lsendrecv_x) then
    call mpicart_create(n_dims,dims(1:n_dims),      &
        &                periodic(1:n_dims),reorder,COMM_CART)

    ! get the rank of the local processor in this new grid
    call mpi_comm_rank(COMM_CART,rank_cart,ierr)

    ! get the coordinates of the local processor in this new grid
    call mpi_cart_coords(COMM_CART,rank_cart,nplan,coords,ierr)
  else
    ! equivalent to COMM_CART = MPI_COMM_WORLD! 
    call mpicopy_world(COMM_CART)
  end if
#else
   COMM_CART = -i_huge+432523
#endif

  ! copy the coords into the plan structure
  j=1
  do i = 1, nplan
    if (pp(i)%parallel) then
      pp(i)%iproc = coords(j)
      j = j + 1
    else
      pp(i)%iproc = 0
    end if
  enddo

  ! provide positions in the processor grid
  iproc_x = pp(ix)%iproc

end subroutine setup_proc_coords

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Create integers corresponding to the first and last points in each
!> direction that the local processor is responsible for. Provide logicals
!> indicating if a processor is responsible for a upper/lowel boundary in the
!> computation domain.
!----------------------------------------------------------------------------
subroutine setup_local_grid_ranges

  use control, only : vp_trap

  integer :: i,j

  !
  ! Use the coords array and the number of points in each direction to
  ! determine the start and end points.
  !

  j=1
  do i=1,nplan
    if (pp(i)%parallel) then
      if (i == ivpar .and. vp_trap == 1) then
        pp(i)%ipb=coords(j)*(pp(i)%points/2) + 1
        pp(i)%ipe=(coords(j)+1)*(pp(i)%points/2)
      else
        pp(i)%ipb=coords(j)*pp(i)%points + 1
        pp(i)%ipe=(coords(j)+1)*pp(i)%points
      end if
      j=j+1
    else
      pp(i)%ipb=1
      pp(i)%ipe=pp(i)%points
    end if
  enddo

  !
  ! copy the values into the variables used elsewhere in the code
  !

  isppb   = pp(isp)%ipb
  isppe   = pp(isp)%ipe
  imupb   = pp(imu)%ipb
  imupe   = pp(imu)%ipe
  ivparpb = pp(ivpar)%ipb
  ivparpe = pp(ivpar)%ipe
  ispb    = pp(is)%ipb
  ispe    = pp(is)%ipe
  ixpb    = pp(ix)%ipb
  ixpe    = pp(ix)%ipe

end subroutine setup_local_grid_ranges

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> It appears that sometimes when performing MPI_CART_SHIFT on a sub cartesian
!> topology, openmpi forgets which dimensions should be periodic. This routine
!> does not fix the topology; it just identifies the processors required for
!> communication in periodic directions. This is done by communicating the
!> ranks of the first and last processors in the grid to eachother.
!----------------------------------------------------------------------------
subroutine openmpi_fix

  use mpicomms,     only : COMM_S_NE, COMM_X_NE
  use mpiinterface, only : mpiallreduce_max, mpicomm_rank

  integer :: rank, flag, sumflag

  ! s-direction
  if (pp(is)%periodic .and. pp(is)%parallel) then

    ! get my rank
    call mpicomm_rank(COMM_S_NE,rank)

    ! I am at the start of the grid; the processor at the end obtains my rank.
    flag = -1
    if (pp(is)%ipb == 1) flag = rank
    call mpiallreduce_max(flag,sumflag,1,COMM_S_NE)
    if (pp(is)%ipe == n_s_grid) proc_s_next = sumflag

    ! I am at the end of the grid; the processor at the start obtains  my rank.
    flag = -1
    if (pp(is)%ipe == n_s_grid) flag = rank
    call mpiallreduce_max(flag,sumflag,1,COMM_S_NE)
    if (pp(is)%ipb == 1) proc_s_prev = sumflag

  end if
  
  ! x-direction
  if (pp(ix)%periodic .and. pp(ix)%parallel) then

    ! get my rank
    call mpicomm_rank(COMM_X_NE,rank)

    ! I am at the start of the grid; the processor at the end obtains my rank.
    flag = -1
    if (pp(ix)%ipb == 1) flag = rank
    call mpiallreduce_max(flag,sumflag,1,COMM_X_NE)
    if (pp(ix)%ipe == n_x_grid) proc_x_next = sumflag

    ! I am at the end of the grid; the processor at the start obtains  my rank.
    flag = -1
    if (pp(ix)%ipe == n_x_grid) flag = rank
    call mpiallreduce_max(flag,sumflag,1,COMM_X_NE)
    if (pp(ix)%ipb == 1) proc_x_prev = sumflag

  end if

end subroutine openmpi_fix


!****************************************************************************
!> get the global species number
!----------------------------------------------------------------------------
function gsp(local_species_number)

  integer, intent(in) :: local_species_number
  integer :: gsp

  gsp = isppb + local_species_number - 1

end function gsp

!****************************************************************************
!> get the global v_parallel grid point number
!----------------------------------------------------------------------------
function gvpar(local_point)

  integer, intent(in) :: local_point
  integer :: gvpar

  gvpar = ivparpb + local_point - 1

end function gvpar

!****************************************************************************
!> get the global mu grid point number
!----------------------------------------------------------------------------
function gmu(local_point)

  integer, intent(in) :: local_point
  integer :: gmu

  gmu = imupb + local_point - 1

end function gmu

!****************************************************************************
!> get the global s point number
!----------------------------------------------------------------------------
function gs(local_point)

  integer, intent(in) :: local_point
  integer :: gs

  gs = ispb + local_point - 1

end function gs

!****************************************************************************
!> get the global x point number
!----------------------------------------------------------------------------
function gx(local_point)

  integer, intent(in) :: local_point
  integer :: gx

  gx = ixpb + local_point - 1

end function gx

!****************************************************************************
!> get the local species number
!----------------------------------------------------------------------------
function lsp(global_species_number)

  integer, intent(in) :: global_species_number
  integer :: lsp

  lsp = global_species_number - isppb + 1

end function lsp

!****************************************************************************
!> get the local v_parallel grid point number
!----------------------------------------------------------------------------
function lvpar(global_point)

  integer, intent(in) :: global_point
  integer :: lvpar

  lvpar = global_point - ivparpb + 1

end function lvpar

!****************************************************************************
!> get the local mu grid point number
!----------------------------------------------------------------------------
function lmu(global_point)

  integer, intent(in) :: global_point
  integer :: lmu

  lmu = global_point - imupb + 1

end function lmu

!****************************************************************************
!> get the local s point number
!----------------------------------------------------------------------------
function ls(global_point)

  integer, intent(in) :: global_point
  integer :: ls

  ls = global_point - ispb + 1

end function ls

!****************************************************************************
!> get the local x point number (named to avoid clash with Lx box size)
!----------------------------------------------------------------------------
function lrx(global_point)

  integer, intent(in) :: global_point
  integer :: lrx

  lrx = global_point - ixpb + 1

end function lrx


!----------------------------------------------------------------------
!> this helper function calculates indexes which can be used to
!> shift Fourier space arrays (cf. for example the fftshift function
!> in Octave/Matlab)
!>
!> (/ (jind_flexible(nx, mrad, i), i = 1,nx)/) reproduces jind from
!> non_linear_terms
!----------------------------------------------------------------------
pure function jind_flexible(ingridsize, outgridsize, i) result(jind)
  use control, only : spectral_radius
  integer, intent(in) :: ingridsize, outgridsize, i
  integer :: jind

  integer :: izero

  if(spectral_radius) then
    jind = 0
    izero = (ingridsize+1)/2

    if(i < izero) then
      jind = outgridsize + i - izero + 1;
    else
      jind = i - izero + 1;
    end if
  else
    jind = i
  end if
end function jind_flexible

!----------------------------------------------------------------------
!> 
!> (jinv_flexible(nx, mrad, i), i = 1,mrad)
!> reproduces jinv from non_linear_terms
!----------------------------------------------------------------------
pure function jinv_flexible(ingridsize, outgridsize, i) result(jinv)
  use control, only : spectral_radius
  integer, intent(in) :: ingridsize, outgridsize, i
  integer :: jinv

  if(spectral_radius) then
    ! The jinv array provides the inverse translation
    do jinv = 1, ingridsize
      if (jind_flexible(ingridsize, outgridsize, jinv) == i) then
        return
      end if
    end do
    jinv = 0
  else
    jinv = i
  end if
end function jinv_flexible

!***************************************************************************
!> Return logical true for processors which contain a specific global grid 
!> point or contain some of a slice of global points.  
!> Use to select a single processor for writing, or a subset for 
!> a particlar integration direction calculation (e.g. moments).
!> No MPI is required (fast). Generalised and replaces previous xy_init 
!> functionality, but is intended to be called when needed, not at init.
!>
!> If called with 0 in an argument, all processors in the direction of
!> that dimension will be returned (e.g. for moment calculations).
!> 
!> For already integrated directions it does not matter which point is selected.
!> All the inputs are global grid points.
!> If no arguments are zero, will always return true on a single processor.
!----------------------------------------------------------------------------
function proc_subset(ixg,isg,imug,ivpg,ispg)
  use general, only : gkw_abort 
  
  integer, intent(in) :: ixg,isg,imug,ivpg,ispg
  logical :: proc_subset, proc_s, proc_mu, proc_vp, proc_sp, proc_x
  
  proc_s = .false.
  proc_mu = .false.
  proc_vp = .false.
  proc_sp = .false.
  proc_x = .false.
  
  !Error checks on input params
  if (isg  <0 .or. isg  > n_s_grid) call gkw_abort('proc_subset: isg ')
  if (imug <0 .or. imug > n_mu_grid) call gkw_abort('proc_subset: imug')
  if (ivpg <0 .or. ivpg > n_vpar_grid) call gkw_abort('proc_subset: ivpg')
  if (ispg <0 .or. ispg > number_of_species) call gkw_abort('proc_subset: ispg')
  if (ixg  <0 .or.  ixg > n_x_grid) call gkw_abort('proc_subset: ixg')

  !Decide which points / directions are needed
  if (isg == 0  .or. (ls(isg)     >=1 .and. ls(isg)<=ns ))         proc_s  = .true. 
  if (imug == 0 .or. (lmu(imug)   >=1 .and. lmu(imug) <=nmu ))     proc_mu = .true.
  if (ivpg == 0 .or. (lvpar(ivpg) >=1 .and. lvpar(ivpg) <=nvpar )) proc_vp = .true. 
  if (ispg == 0 .or. (lsp(ispg)   >=1 .and. lsp(ispg) <=nsp ))     proc_sp = .true.
  if (ixg == 0  .or. (lrx(ixg)    >=1 .and. lrx(ixg) <=nx ))       proc_x  = .true.
  
  !Take intersection of dimensions
  proc_subset = proc_s .and. proc_mu .and. proc_vp .and. proc_sp .and. proc_x

end function proc_subset

end module grid
