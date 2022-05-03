!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Control contains all the top level switches of the code.  
!> Control is the top most main code module.
!> (it is below only mpiinterface and ompinterface containing the basic mpi 
!> parameters, and general, which contains general purpose routines).
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module control

  use global, only : lenswitch
  
  implicit none

  private

  !
  ! publicly available procedures
  !

  public :: control_init, control_initt, control_read_nml, control_bcast_nml
  public :: control_check_params, control_write_nml
  public :: get_start_time, get_max_runtime, max_runtime_reached

  !
  ! publicly available variables
  !

  !> True if collisions are to be used
  logical, save, public :: lcollisions
  logical, save :: collisions

  !< True if the neoclassical effects are to be calculated 
  logical, save, public :: neoclassics
  !> True if the zonal flows are used in the equation for the adiabatic response
  logical, save, public :: zonal_adiabatic
  !< true for including the nonlinear terms
  logical, save, public :: non_linear
  !> True if the electrostatic potential is kept in the equations
  logical, save, public :: nlphi
  !> True if A|| is kept in the equations 
  logical, save, public :: nlapar
  !> True if B|| is kept in the equations 
  logical, save, public :: nlbpar
  !> True if a spectral method is used for the radial direction 
  logical, save, public :: spectral_radius
  !> True for the flux tube case 
  logical, save, public :: flux_tube 
  !> True for the parallel velocity nonlinearity (experimental) 
  logical, save, public :: lpar_vel_nl = .false. 
  !> True if one the simulation does not include the background 
  logical, save, public :: fullf_wo_Fm = .false. 
  !> True if the code uses the shifted metric only for the non spectral case
  logical, save, public :: shift_metric = .false. 

  !> True if the energetics diagnostic is switched on.
  logical, save, public :: lcalc_energetics = .false.
  
  !> Putting the two definitions below into control is necessary to 
  !> allow for memory optimization.
  !> True if the zonal evolution diagnostic is switched on.
  logical, save, public :: lcalc_zonal_evo = .false.
  !> True if zonal evolution diagnostic should provide detailed output. 
  logical, save, public :: zevo_detail = .false.
  
  
  !> Some diagnostics use individual linear terms matrices. The switches
  !> below are used for memory optimization.
  !> True if matd is needed by diagnostics
  logical, save, public :: l_matd = .false.
  !> True if matvpd is needed by diagnostics
  logical, save, public :: l_matvpd = .false.
  !> True if matperpd is needed by diagnostics
  logical, save, public :: l_matperpd = .false.
  !> True if matcoll is needed by diagnostics
  logical, save, public :: l_matcoll = .false.
  !> True if matoutflow is needed by diagnostics
  logical, save, public :: l_matoutflow = .false.
  !> True if mat_vpar_grad_df is needed by diagnostics
  logical, save, public :: l_mat_vpar_grad_df = .false.
  !> True if mat_vdgradf is needed by diagnostics
  logical, save, public :: l_mat_vdgradf = .false.
  !> True if mat_trapdf_4d is needed by diagnostics
  logical, save, public :: l_mat_trapdf_4d = .false.
  !> True if mat_ve_grad_fm is needed by diagnostics
  logical, save, public :: l_mat_ve_grad_fm = .false.
  !> True if mat_vpar_grd_phi is needed by diagnostics
  logical, save, public :: l_mat_vpar_grd_phi = .false.
  !> True if mat_vd_grad_phi_fm is needed by diagnostics
  logical, save, public :: l_mat_vd_grad_phi_fm = .false.
  !> True if mat_vpar_grd_bpar is needed by diagnostics
  logical, save, public :: l_mat_vpar_grd_bpar = .false.
  !> True if mat_vd_grad_bpar_fm is needed by diagnostics
  logical, save, public :: l_mat_vd_grad_bpar_fm = .false.
  
  
  !> True if distribution function is averaged over time
  logical, save, public :: laverage_dist_over_time

  !> True if the code decides it needs to stop 
  logical, save, public :: stop_me = .false.
  logical, save, public :: nan_stop
  !> True if the nonlinear timestep estimator is to be used (nonlinear only) 
  logical, save, public :: nl_dtim_est
  !> True if shear-periodic boundaries are used in non spectral case
  logical, save, public :: shear_periodic = .false.  
  !> Factors to multiply the timesteps estimators by
  real, save, public :: fac_dtim_est, fac_dtim_nl
  !> Disipation / upwind parameter for the parallel (to the field) derivatives
  real, save, public :: disp_par 
  !> Dissipation / upwind parameter for the parallel velocity 
  real, save, public :: disp_vp
  !> Perpendicular dissipation in x
  real, save, public :: disp_x
  !> Perpendicular dissipation in y
  real, save, public :: disp_y

  !> Integer that determines how the parallel velocity grid is set up. 
  !> If = 0, the parallel velocity of a grid point is constant along the field
  !> line and the grid is uniform. 
  !> If = 1, the parallel velocity follows the trapping condition. 
  integer, save, public :: vp_trap
  !Makes the spacing between mu grid points equal -> Useful for collisional runs.
  logical, save, public :: uniform_mu_grid
  !'Flappy' boundary conditions in velocity space switch.  
  logical, save, public :: lflapv
  !> Switch for the parallel boundary conditions, can be "Dirichlet" or "open" (recommended)
  character (len = lenswitch), save, public :: parallel_boundary_conditions
  !> Radial boundary conditions 'periodic', 'Dirichlet' or 'Neu-Dir' 
  character (len = lenswitch), save, public :: radial_boundary_conditions
  !> Selects the order of the numerical scheme (accuracy). Allowed are
  !> 'second_order' and 'fourth_order'
  character (len = lenswitch), save, public :: order_of_the_scheme
  character (len = lenswitch), save, public :: order_of_the_radial_scheme
  character (len = lenswitch), save, public :: order_of_the_zf_scheme

  !> one of several sparse matrix formats
  character (len = lenswitch), save, public :: matrix_format

  !> True if a normalization is applied for the distrubtion function. Nonlinear
  !> runs are not normalized.
  logical, save, public :: normalized
  !> normalise per toroidal mode
  logical, save, public :: normalize_per_toroidal_mode
  !> Use Arakawa type differencing for trapping terms
  logical, save, public :: ltrapping_arakawa
  !> another switch to finetune screen output: no printing of the
  !> matrix compression information. 'silent' appears in the namelist,
  !> root_and_not_silent should be used in conditionals.
  logical, save, public :: silent
  logical, save, public :: root_and_not_silent
  !> set to true only when testing code MPI performance
  logical, save, public :: testing = .false.

  ! time-related parameters:
  
  !> total number of large timesteps to run in this run
  integer, save, public :: ntime
  !> number of small timesteps per large timestep
  integer, save, public :: naverage
  !> number of large timesteps between checkpoint dumps
  integer, save, public :: ndump_ts
  !> timestep size in normalized time units
  real,    save, public :: dtim
  !> estimated timestep size for nonlinear terms stability
  real,    save, public :: dtim_est
  !> keeps local minimum timestep for nonlinear terms
  real,    save, public :: dtim_est_save
  !> the original input time step
  real,    save, public :: dtim_input

  ! counters:
  
  !> Total time
  real, save, public :: time
  !> number of current large timestep, since the beginning of the
  !> current run
  integer, save, public :: itime
  !> number of large timesteps, since the last restart file has been
  !> written
  integer, save, public :: itime_rst
  !> number of current small timestep, including previous runs; note
  !> that int(ntotstep/naverage) is not the same as nt_complete+itime
  !> for restarted runs whose naverage has been changed.
  integer, save, public :: ntotstep

  !> Total time at the last large time step
  real, save, public :: last_largestep_time
  !> Total time at the last small time step
  real, save, public :: last_smallstep_time

  !> time completed by earlier runs (i.e. time at which the current
  !> run began); this number is updated at restart and does not change
  !> during the run.
  real, save, public :: t_init
  !> number of timesteps completed by earlier runs; this number is 
  !> updated when restart files are written, including both restart 
  !> files that are written in the end of a run and automatic restart 
  !> dump files written during the run.
  integer, save, public :: nt_complete=0

  !> method for solving
  character (len=lenswitch), save, public :: method
  !> choice of algorithm for method
  integer, save, public :: meth

  !> if true, code attempts to resume a run from FDS file
  logical, save, public :: read_file
  !> code will restart from a checkpoint file DMP if present.
  !> DMP file (if present) overrides FDS.
  logical, save, public :: auto_restart
  !> if restarted code is allowed to use a new gridsize
  logical, save, public :: lrestart_new_grid
  
  !> Maximum number of seconds for a run
  real,    save, public :: max_seconds
  !> Maximum number of second for a run (integer)
  integer, save, public :: max_sec
  
  !> Tolerance in gamma that stops the code.
  real,    save, public :: gamatol
  !> Tolerance in es fluxes that stops code.
  real,    save, public :: fluxtol
  !> Which column of es fluxes is used to stop code.
  integer, save, public :: ifluxtol
  !> Tolerance in neoclassical heat flux that stops code.
  real,    save, public :: ncqtol
  !> Minimum growth rate that stops the code.
  real,    save, public :: min_gr
  !> Maximum growth rate that stops the code.
  real,    save, public :: max_gr
  !> mimimum value of dt, below which code stops.
  real,    save, public :: dt_min

  !> restart file version (0 = no files written)
  integer, save, public :: restart_file_version
  !> current restart file version
  integer, parameter, public :: restart_file_current_version = 2
  !> run number; used for (optional) sequential number of restart files. 
  integer, save, public :: irun
  !> a guess for the time needed to write the final restart file.
  !> this value is updated when dump files are written.
  real, save, public :: max_t_fdis_write = 0.0

  !> a file, whose existence makes the iteration stop and gkw exit cleanly
  character(len=8), parameter, public :: stop_filename = 'gkw.stop'


  !> The file format for data output, e.g. data from diagnostics.
  character (len=lenswitch), save, public :: io_format

  !> If this switch is set to true, some diagnostics output their data
  !> in the old, legacy, traditional ordering, rather than a
  !> (intentionally) more natural ordering.
  logical, save, public :: io_legacy

  !> This switch causes very small numbers in the output data to be
  !> rounded to zero. This is useful to make testcases pass more
  !> robustly on different machines.
  logical, save, public :: io_testdata

  !> Integer that determines for which part of the code the performance is 
  !> measured
  integer, save, public :: iperform_set

  !
  ! interfaces
  !

  ! use one routine for reading and writing
  interface control_write_nml
    module procedure control_read_nml
  end interface

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Anything needed for initialisation of control
!----------------------------------------------------------------------------
subroutine control_init


end subroutine control_init

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reads from the input file the various control switches; writes the
!> namelist to file if the optional switch is *not* present.
!----------------------------------------------------------------------------
subroutine control_read_nml(ilun,io_stat,lwrite)
  use global, only : lverbose, logical_false, compiled_with_hdf5
  use global, only : compiled_with_librsb, compiled_with_mkl
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  use general, only : gkw_warn
  
  integer, intent(in)  :: ilun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite
  
  namelist /control/ disp_par, disp_vp, disp_x, disp_y,       &  !Numerical dissipation
      & dtim, fac_dtim_est, fac_dtim_nl, nl_dtim_est,         &  !Timestep 
      & lverbose, silent,                                     &  !Screen output
      & meth, method, order_of_the_scheme, spectral_radius,   &  !Numerical scheme
      & ntime, naverage, ndump_ts, nlapar, nlbpar, nlphi,     &  !Loops counts, Physics
      & collisions, neoclassics, zonal_adiabatic,             &  !Physics
      & non_linear, flux_tube, shift_metric, lpar_vel_nl,     &  !Physics 
      & normalized, normalize_per_toroidal_mode,              &  !Normalization
      & read_file, restart_file_version, auto_restart, irun,  &  !Restart options
      & lrestart_new_grid, vp_trap, uniform_mu_grid, lflapv,  &  !Velocity grids
      & max_sec, max_seconds, gamatol, fluxtol, ifluxtol, ncqtol,   &  !Stop tolerances
      & max_gr, min_gr, dt_min,                               &  !Stop tolerances
      & ltrapping_arakawa,                                    &  !Experimental numerics
      & radial_boundary_conditions,                           &  !boundary conditions 
      & order_of_the_radial_scheme, order_of_the_zf_scheme,   &
      & io_format, io_legacy, io_testdata,                    &  !input/output
      & iperform_set,                                         &  !performance measure
      & matrix_format,                                        &  !Sparse Matrix Format
      & testing, parallel_boundary_conditions,                &  !Obsolete
      & laverage_dist_over_time
      
  ! keep compiler quiet
  if (logical_false) write(*,*) matrix_format

  io_stat = 0
  ! read the input; this needs the switch
  if (present(lwrite)) then

    if (.not. lwrite) then
      ! Set the default values for the control parameters; the default
      ! corresponds to a linear run without resume.
      max_seconds       = -1.
      lverbose          = .false.
      nl_dtim_est       = .true.
      fac_dtim_est      =  0.95
      fac_dtim_nl       =  1.0
      read_file         = .false.
      silent            = .false. 
      zonal_adiabatic   = .false.
      collisions        = .false. 
      nlphi             = .true. 
      nlapar            = .false.
      nlbpar            = .false.
      spectral_radius   = .true. 
      shift_metric      = .false.
      neoclassics       = .false. 
      non_linear        = .false.
      flux_tube         = .true. 
      lpar_vel_nl       = .false. 
      normalized        = .true.
      normalize_per_toroidal_mode = .false.
      auto_restart      = .false.
      lrestart_new_grid = .false.
      ndump_ts          = 0
      vp_trap           = 0 
      uniform_mu_grid   = .false.
      lflapv            = .false.
      ltrapping_arakawa = .false.
      radial_boundary_conditions = 'periodic'
      restart_file_version = restart_file_current_version
      gamatol           = 0.
      fluxtol           = 0.
      ifluxtol          = 2
      ncqtol            = 0.
      dt_min            = 1e-6
      min_gr            = 0.01
      max_gr            = 100.0
      
      parallel_boundary_conditions = 'open'
      ! The differentiation is by default 4th order 
      order_of_the_scheme = 'fourth_order'
      order_of_the_radial_scheme = 'unstated'
      order_of_the_zf_scheme = 'unstated'
      
      laverage_dist_over_time = .false.

      irun       = 0 ! run number
      ntime      = 0 ! number of large timesteps
      naverage   = 0 ! number of small time steps. Total number is ntime*naverage 
      dtim       = 0.005 ! time step 
      disp_par   = 0.2E0 ! The dissipation coefficient for parallel derivatives 
      disp_vp    = 0.2E0 ! The dissipation coefficient for parallel velocity space 
      disp_x     = 0.0E0 ! 'Radial' perpendicular dissipation coeffcient
      disp_y     = 0.0E0 ! 'Poloidal' perpendicular dissipation coeffcient
      method     = 'EXP' ! method of integration
      meth       = 2     ! switch between different schemes
      max_sec    = -1    ! maximum seconds for a run (<=0) is infinite 
      ! initialize also the total number of steps taken 
      ntotstep = 0

      iperform_set = 2

      if(compiled_with_hdf5) then
        io_format = 'hdf5+ascii'
      else
        io_format = 'mixed'
      end if
      io_legacy = .true.
      io_testdata = .false.

      ! automatically use a certain library if available
      ! if(compiled_with_mkl) then
      !   matrix_format = 'mkl-crs'
      ! else if(compiled_with_librsb) then
      !   matrix_format = 'librsb'
      ! else
      !   matrix_format = 'gkw-crs'
      ! end if
      ! always use inhouse sparse matrix routines by default
      matrix_format = 'gkw-crs'

      ! read namelist
      read(ilun,NML=control,IOSTAT=io_stat)

      if(trim(matrix_format) == 'complex') then
        call gkw_warn('matrix_format = ''complex'' is outdated (matrix is always&
        & complex). Default matrix format is used.')
        if(compiled_with_mkl) then
          matrix_format = 'mkl-crs'
        else if(compiled_with_librsb) then
          matrix_format = 'librsb'
        else
          matrix_format = 'gkw-crs'
        end if
      end if

      if(.not. spectral_radius .and. matrix_format /= 'gkw-crs') then
        call gkw_warn('For nonspectral runs, only matrix_format = ''gkw-crs'' is&
           & allowed. Forcing this.')
        ! This is because in fields, usmv is called at several place
        ! with the identical vector in two arguments. This is not
        ! allowed for the external libraries.
         matrix_format = 'gkw-crs'
      end if

      if(order_of_the_radial_scheme == 'unstated') then
        order_of_the_radial_scheme = order_of_the_scheme
      end if
      if(order_of_the_zf_scheme == 'unstated') then
        order_of_the_zf_scheme = order_of_the_scheme
      end if

      if(testing) then
        ! For mpi testing, disable output to files
        io_format = 'none'
        write (*,*) "* output disabled:", &
           & "io_format = 'none' because testing is enabled"

        ! For mpi testing, need only one large loop
        ntime = 1
      end if

    end if

  else

    ! write to input.out if called without the switch; this is the default

    ! write the namelist
    if(root_processor) write(ilun,NML=control)

    ! further checks

    if(normalize_per_toroidal_mode .and. normalized) then
      ! normalize=T means normalization with a single factor.
      ! normalize_per_toroidal_mode is more specialised and thus
      ! takes precedence if both are used.
      normalized = .false.
    end if

    ! write metadata
    call write_run_parameter('control', 'disp_par', disp_par)
    call write_run_parameter('control', 'disp_vp', disp_vp)
    call write_run_parameter('control', 'disp_x', disp_x)
    call write_run_parameter('control', 'disp_y', disp_y)
    call write_run_parameter('control', 'dtim', dtim)
    call write_run_parameter('control', 'fac_dtim_est', fac_dtim_est)
    call write_run_parameter('control', 'nl_dtim_est', nl_dtim_est)
    call write_run_parameter('control', 'lverbose', lverbose)
    call write_run_parameter('control', 'silent', silent)
    call write_run_parameter('control', 'meth', meth)
    call write_run_parameter('control', 'method', method)
    call write_run_parameter('control', 'order_of_the_scheme', order_of_the_scheme)
    call write_run_parameter('control', 'spectral_radius', spectral_radius)
    call write_run_parameter('control', 'ntime', ntime)
    call write_run_parameter('control', 'naverage', naverage)
    call write_run_parameter('control', 'ndump_ts', ndump_ts)
    call write_run_parameter('control', 'nlapar', nlapar)
    call write_run_parameter('control', 'nlbpar', nlbpar)
    call write_run_parameter('control', 'nlphi', nlphi)
    call write_run_parameter('control', 'collisions', collisions)
    call write_run_parameter('control', 'neoclassics', neoclassics)
    call write_run_parameter('control', 'normalized', normalized)
    call write_run_parameter('control', 'normalize_per_toroidal_mode', &
       & normalize_per_toroidal_mode)
    call write_run_parameter('control', 'zonal_adiabatic', zonal_adiabatic)
    call write_run_parameter('control', 'non_linear', non_linear)
    call write_run_parameter('control', 'flux_tube', flux_tube)
    call write_run_parameter('control', 'lpar_vel_nl', lpar_vel_nl)
    call write_run_parameter('control', 'shift_metric', shift_metric)
    call write_run_parameter('control', 'read_file', read_file)
    call write_run_parameter('control', 'restart_file_version', restart_file_version)
    call write_run_parameter('control', 'auto_restart', auto_restart)
    call write_run_parameter('control', 'irun', irun)
    call write_run_parameter('control', 'lrestart_new_grid', lrestart_new_grid)
    call write_run_parameter('control', 'vp_trap', vp_trap)
    call write_run_parameter('control', 'uniform_mu_grid', uniform_mu_grid)
    call write_run_parameter('control', 'lflapv', lflapv)
    call write_run_parameter('control', 'max_sec', max_sec)
    call write_run_parameter('control', 'max_seconds', max_seconds)
    call write_run_parameter('control', 'gamatol', gamatol)
    call write_run_parameter('control', 'fluxtol', fluxtol)
    call write_run_parameter('control', 'ifluxtol', ifluxtol)
    call write_run_parameter('control', 'ncqtol', ncqtol)
    call write_run_parameter('control', 'max_gr', max_gr)
    call write_run_parameter('control', 'min_gr', min_gr)
    call write_run_parameter('control', 'dt_min', dt_min)
    call write_run_parameter('control', 'ltrapping_arakawa', ltrapping_arakawa)
    call write_run_parameter('control', 'radial_boundary_conditions', radial_boundary_conditions)
    call write_run_parameter('control', 'order_of_the_radial_scheme', order_of_the_radial_scheme)
    call write_run_parameter('control', 'order_of_the_zf_scheme', order_of_the_zf_scheme)
    call write_run_parameter('control', 'io_format', io_format)
    call write_run_parameter('control', 'io_legacy', io_legacy)
    call write_run_parameter('control', 'io_testdata', io_testdata)
    call write_run_parameter('control', 'testing', testing)
    call write_run_parameter('control', 'parallel_boundary_conditions', parallel_boundary_conditions)
    call write_run_parameter('control', 'iperform_set', iperform_set)
    call write_run_parameter('control', 'matrix_format', matrix_format)
    call write_run_parameter('control', 'laverage_dist_over_time', laverage_dist_over_time)
  end if

end subroutine control_read_nml


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Broadcast control details to other processors
!----------------------------------------------------------------------------
subroutine control_bcast_nml

  use mpiinterface, only : mpibcast
  use global,       only : lverbose

  call mpibcast(irun,1)
  call mpibcast(dt_min,1)
  call mpibcast(ntime,1)
  call mpibcast(naverage,1)
  call mpibcast(ndump_ts,1)
  call mpibcast(max_sec,1)
  call mpibcast(vp_trap,1)
  call mpibcast(uniform_mu_grid,1)
  call mpibcast(lflapv,1)
  call mpibcast(dtim,1)
  call mpibcast(max_seconds,1)
  call mpibcast(disp_par,1)
  call mpibcast(disp_vp,1)
  call mpibcast(disp_x,1)
  call mpibcast(disp_y,1)
  call mpibcast(method,lenswitch)
  call mpibcast(meth,1)  
  call mpibcast(nlphi,          1) 
  call mpibcast(nlapar,         1) 
  call mpibcast(spectral_radius,1)
  call mpibcast(shift_metric,   1)
  call mpibcast(flux_tube,      1) 
  call mpibcast(lpar_vel_nl,    1)
  call mpibcast(nlbpar,         1) 
  call mpibcast(normalized,     1) 
  call mpibcast(normalize_per_toroidal_mode,1)
  call mpibcast(read_file,      1) 
  call mpibcast(neoclassics,    1) 
  call mpibcast(collisions,     1) 
  call mpibcast(non_linear,     1)
  call mpibcast(nl_dtim_est,    1)
  call mpibcast(fac_dtim_est,   1)
  call mpibcast(fac_dtim_nl,    1)
  call mpibcast(silent,         1)
  call mpibcast(zonal_adiabatic,1) 
  call mpibcast(ltrapping_arakawa,1) 
  call mpibcast(auto_restart,1)
  call mpibcast(lrestart_new_grid,1)
  call mpibcast(parallel_boundary_conditions, lenswitch) 
  call mpibcast(radial_boundary_conditions,   lenswitch) 
  call mpibcast(order_of_the_scheme,          lenswitch)
  call mpibcast(order_of_the_radial_scheme,   lenswitch)
  call mpibcast(order_of_the_zf_scheme,       lenswitch)
  call mpibcast(lverbose,1)
  call mpibcast(restart_file_version,1)
  call mpibcast(gamatol,1)
  call mpibcast(fluxtol,1)
  call mpibcast(ifluxtol,1)
  call mpibcast(ncqtol,1)
  call mpibcast(min_gr,1)
  call mpibcast(max_gr,1)
  call mpibcast(testing,1)
  
  call mpibcast(laverage_dist_over_time,1)

  call mpibcast(io_format, lenswitch)
  call mpibcast(io_legacy, 1)
  call mpibcast(io_testdata, 1)

  call mpibcast(iperform_set,1)

  call mpibcast(matrix_format,lenswitch)


  lcollisions = collisions

end subroutine control_bcast_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine control_initt

  ! Set switches for allocation of individual linear terms matrices.
  ! Some of them are used by diagnos_energetics.
  if(lcalc_energetics) then
    l_matd = .true.
    l_matvpd = .true.
    l_matperpd = .true.
    l_matcoll = .true.
    l_matoutflow = .true.
  end if
  
  ! Set switches for allocation of individual linear terms matrices.
  ! Some of them are used by diagnos_zonal_evo.
  if(lcalc_zonal_evo .and. zevo_detail) then
  
    l_mat_vpar_grad_df = .true.
    l_mat_vdgradf = .true.
    l_mat_trapdf_4d = .true.
    l_mat_ve_grad_fm = .true.
    l_mat_vpar_grd_phi = .true.
    l_mat_vd_grad_phi_fm = .true.
    
    l_matd = .true.
    l_matvpd = .true.
    l_matperpd = .true.
    
    if(nlbpar) then
      l_mat_vpar_grd_bpar = .true.
      l_mat_vd_grad_bpar_fm = .true.
    end if
    
    if(lcollisions) then
      l_matcoll = .true.
    end if
    
  end if
  

end subroutine control_initt

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Check the input and initialise everything necessary for allocation.
!----------------------------------------------------------------------------
subroutine control_check_params

  use fft,          only : working_fft_library
  use general,      only : gkw_exit, gkw_warn
  use mpiinterface, only : root_processor
  use global,       only : lverbose, root_and_verbose

  root_and_verbose = root_processor .and. lverbose
  root_and_not_silent = root_processor .and. .not.silent

  if (root_and_verbose) write (*,*) '* verbose output enabled'

  !bad values for vp_trap
  if (.not. ( vp_trap == 0 .or. vp_trap == 1) ) then
    call gkw_exit('Control: vp_trap must be 1 or 0')
  endif

  ! Do we need a FFT library?
  if (non_linear .and. (.not. working_fft_library)) then
    call gkw_exit('non_linear option requires a working FFT '//&
        &          'library!')
  end if

  !method checks
  select case(method)
  case('EXP')
    select case(meth)
    case(1)
      if (non_linear) call gkw_warn('Meth=1 Nonlinear timestep estimator off?')
    case(2,99) !RK4 This is recommended
      !No warnings
    case(3) 
      if (non_linear) call gkw_warn('Meth=3 Nonlinear timestep estimator off?')
    case(-4,-6,-7)
      call gkw_warn('Experimental rkc time integration selected')
    case default
      call gkw_exit('Unknown value of meth for method EXP')
    end select !meth

  case('IMP')
    call gkw_warn('Implicit scheme under development &
        & many options not supported: (nonlinear, parallel, ...)')
    !Should add throw outs here after more development before release.

  case('EIV')

  case default 
    if (root_processor) then
      write(*,*)'control: You have given ',method,' as method'
      write(*,*)'allowed are: EIV, EXP, IMP'
    endif
    call gkw_exit('(see message above)')
  end select  !method

  ! parallel boundary conditions
  select case(parallel_boundary_conditions)
  case('zero') ! 
    call gkw_warn('Control: parallel_boundary_conditions reset to "open"')
    call gkw_warn('Control: The "zero" option is no longer recognised')    
    parallel_boundary_conditions = 'open'    
  case('open')      ! The default  
  case('Dirichlet') ! not recommended due to oscillations at boundary  
  case default
    call gkw_warn('Control: parallel_boundary_conditions unknown option:')
    call gkw_exit('Reset to either Dirichlet, or open (recommended)')
  end select !parallel boundary conditions

  ! radial boundary conditions 
  select case(radial_boundary_conditions)
  case('Dirichlet') 
  case('periodic')  ! The default
  case('Neu-Dir')
  case('Neuslab')
  case default
    call gkw_warn('Control: radial_boundary_conditions unknown option:')
    call gkw_exit('Reset to either Dirichlet, periodic or Neu-Dir')
  end select !radial boundary conditions

  ! order of the scheme
  select case(order_of_the_scheme)

  case('second_order')

    if (ltrapping_arakawa) then
      call gkw_exit('second_order not fully implemented for ltrapping_'// &
          &          'arakawa')
    end if

    ! Need to correct normalizations and add the trapping due to Apar F_M G
    ! term
    !APS: what is this above comment about?
  case('fourth_order')   
  case default
    if (root_processor) then 
      write(*,*) 'Control: you specified ',order_of_the_scheme,' as input'
      write(*,*) 'to order_of_the_scheme.'
      write(*,*) 'Only known options are: second_order and fourth_order'
    endif
    call gkw_exit('(see message above)')
  end select !order of the scheme

  ! store the value given in the input file 
  dtim_input = dtim
  !The estimate value must be initialised larger than the input value. 
  dtim_est_save=dtim+1.
  !Do not use the nonlinear timestep estimator for linear runs
  if (.not.non_linear) nl_dtim_est=.false. 

  ! intialize the time to zero 
  time = 0.
  last_largestep_time = 0 ! or better -dtim ?
  last_smallstep_time = 0 ! or better -dtim ?
  t_init=0.

  ! restart file version
  if (restart_file_version < 0 .or.                                          &
      &  restart_file_version > restart_file_current_version) then
    call gkw_exit('bad restart_file_version, must be between 0 and 2')
  end if

  ! check for vp_trap = 0 and collisions, not implemented
  if (vp_trap /= 0 .and. collisions) then
    call gkw_exit('control: vp_trap=1 does not work with collisions')
  endif

  ! this combination makes no sense
  if (vp_trap /= 0 .and. ltrapping_arakawa) then
   call gkw_exit('control: vp_trap=1 incompatible with arakawa scheme')
  endif

  ! vp_trap is currently numerically unstable in nonlinear case
  if (vp_trap /= 0 .and. non_linear .and. method .eq. 'EXP') then
!   call gkw_exit('control: vp_trap unstable for nonlinear explicit runs')
  endif

  if (min_gr > max_gr) then
     call gkw_warn('Control: min_gr cannot be larger than max_gr, override')
     min_gr=-10.0
     max_gr=100.0
  end if
  
  if (dtim < dt_min) then
    call gkw_warn('Control: Reducing dt_min to value of input dtim')
    dt_min = dtim_input
  end if 
  
  if (dt_min < 1e-6 ) then
     call gkw_warn('Control: The code is unlikely to work correctly if ' // &
                  & 'the timestep gets smaller than 1e-6, but try normalized=.false.')
  end if
  
  ! remedy for legacy inputs with fac_dtim_est > 1.0
  if (fac_dtim_est > 1.0) then
     call gkw_warn('Control: fac_dtim_est > 1.0 always unstable, since r3553')
     call gkw_warn('Control: resetting fac_dtim_est = 0.98')
     call gkw_warn('Control: To override, use negative fac_dtim_est')
     fac_dtim_est = 0.98
  else
     fac_dtim_est = abs(fac_dtim_est)
  end if
  
  if (neoclassics) then
     min_gr=-100.
  end if
  
  !Negative dissipation coeffcients are not allowed
  if (disp_x < 0. .or. disp_y < 0.) then
    call gkw_warn('Control: Negative x/y dissipation treated as k^2 dissipation')
  endif
  if (disp_vp < 0. .or. disp_par < 0.) then
    call gkw_exit('Control: Negative vp/par dissipation coeffcients not allowed')
  endif
  
  !Small dissipation coeffcients are not recommended
  if (abs(disp_x) + abs(disp_y) < 0.1 .and. non_linear) then
    call gkw_warn('Control: x/y dissipation advised NL runs: Check high k spectra')
  endif

  ! set the maximum number of seconds
  if (max_seconds > 0. .and. max_sec > 0) then
    call gkw_warn('max_seconds will be used, not max_sec')
  end if
  if (max_seconds > 0.) then
    ! make them appear equal in input.out
    max_sec=int(max_seconds)
  else
    max_seconds=real(max_sec)
  end if

  ! shift_metric not valid for spectral_radius
  if (shift_metric .and. spectral_radius) then
    call gkw_warn('shift metric not possible for spectral_radius')
    shift_metric = .false.
  end if
  
  if ((.not. flux_tube) .and. spectral_radius) then
    call gkw_exit('Global runs require spectral_radius = F')
  end if
  
  if (radial_boundary_conditions == 'Dirichlet' .and. spectral_radius) then 
    call gkw_exit('Dirichlet BC only possible for spectral_radius = F')
  end if
  
  !if ((.not.spectral_radius).and. non_linear .and. nlapar) then
  !  call gkw_warn('Electromagnetic flutter not yet tested in nonspectral NL terms')
  !end if
  
  if (ifluxtol < 1) ifluxtol = 1

end subroutine control_check_params

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function get_start_time()
  use mpiinterface, only : mpiwtime
  double precision, save :: t_begin = 0.0
  double precision :: get_start_time
  if(t_begin <= 0.0) then
    ! note the start time of the run
    t_begin = mpiwtime()
  end if
  get_start_time = t_begin
end function get_start_time


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function get_max_runtime(deltat_ahead_of_max)
  real, intent(in) :: deltat_ahead_of_max
  real :: get_max_runtime
  real :: extra_time

  ! make the run stop at least 5 minutes before walltime, if we have
  ! no idea of how long it takes to write a restart file.  if it turns
  ! out that writing a restart file takes longer than this, you better
  ! have requested more time (than precisely max_seconds) from the
  ! queueing system.
  !
  ! this also fails, if the eigenvalue solver needs too much tume
  ! between calls of the custom stopping test.
  !
  ! Thus, be aware that in particular for large runs, it is necessary
  ! to request more time than just max_seconds from the queueing
  ! system, to safely finish the run.
  
  extra_time = 5 * 60
  if(max_seconds < 2*extra_time) then
    ! but make it work somehow to set short values like max_seconds = 3 minutes
    extra_time = 0
  end if
  get_max_runtime = max_seconds - deltat_ahead_of_max - max_t_fdis_write - extra_time

end function get_max_runtime
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function max_runtime_reached(deltat_ahead_of_max)
  use mpiinterface, only : mpiwtime
  real, intent(in) :: deltat_ahead_of_max
  logical :: max_runtime_reached
  
  double precision :: t_current
  if (max_seconds > 0.) then
    t_current = mpiwtime()
    if (t_current - get_start_time() > get_max_runtime(deltat_ahead_of_max)) then
      max_runtime_reached = .true.
      return
    end if
  end if

  max_runtime_reached = .false.
end function max_runtime_reached

end module control
