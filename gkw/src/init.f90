!-----------------------------------------------------------------------------
!> 
!> This module reads the input file, and checks for errors
!> It then calls the initialisation routines
!> from all other modules in the correct order.
!>
!-----------------------------------------------------------------------------
module init
  
  implicit none

  private
  ! public subroutines
  public :: initialize, deallocate_runtime_arrays, finalize

contains

! ----------------------------------------------------------------------------
!> This subroutine should do the initalization of all 
!> the quantities
!-----------------------------------------------------------------------------

subroutine initialize(check_input)

  use general,         only : gkw_exit
  use control,         only : spectral_radius, flux_tube
  use control,         only : read_file, auto_restart, testing 
  use control,         only : method, control_initt, irun
  use control,         only : io_format, io_legacy, io_testdata
  use grid,            only : setup_grid, n_x_grid 
  use mode,            only : mode_init, mode_box_recon, mode_check_params
  use mode,            only : kgrid, krbal, kgrid_nonspec_rad
  use normalise,       only : normalise_init, normalise_fdisi
  use diagnostic,      only : diagnostic_allocate, diagnostic_read_last_data
  use diagnostic,      only : diagnostic_init, diagnostic_check_params
  use diagnostic,      only : abort_if_old_data_exists
  use restart,         only : restart_init, read_restart_file, restarted
  use dist,            only : dist_init
  use dist,            only : dist_init, fdisi, nsolc
  use components,      only : components_input_species, components_profiles
  use components,      only : components_beta_cal, rhostar
  use geom,            only : geom_init_grids, parallelize_geom
  use geom,            only : geom_output, geom_type
  use geom,            only : shift_end_grid
  use exp_integration, only : init_explicit
  use eiv_integration, only : init_eiv
  use linear_terms,    only : init_linear_terms, calc_linear_terms
  use non_linear_terms,only : nonlinear_init 
  use rotation,        only : rotation_init, need_fft, parallelize_rotation
  use velocitygrid,    only : velgrid_init
  use tearingmodes,    only : tearingmodes_init
  use collisionop,     only : mom_conservation,ene_conservation
  use gyro_average,    only : gyro_average_init, gyro_average_allocate
  use gyro_average,    only : blending_order
  use krook,           only : krook_init, nlkrook 
  use mpiinterface,    only : processor_number
  use io,              only : init_io
  use matrix_format,   only : matrix_format_init
  use components,      only : imod_init
  use exp_integration, only : set_persistent_mode
  use fields, only : calculate_fields, field_solve_nonspec_wrap
  use control, only : spectral_radius 
  use components, only : mode_persist

  
  ! for input checking, call this routine with check_input=.true.
  logical, optional, intent(in) :: check_input
  
  integer :: ix

  ! Read, broadcast and check all the namelist items; initialize anything
  ! necessary before allocation.
  call read_params

  ! initialize the IO interfaces
  call init_io(processor_number, io_format, io_testdata, &
     & read_file .or. auto_restart, &
     & io_legacy)

  ! anything else in control
  call control_initt

  call matrix_format_init()

  call setup_grid

  ! Exit here if the code is compiled only to check the input
  if (present(check_input)) then
    if (check_input) then
      call gkw_exit('Input file read sucessfully, input.out written')
    end if
  end if

  ! This routine also sets up and allocates everything required in
  ! dist 
  call dist_init(mom_conservation,ene_conservation,blending_order)

  call restart_init
  
  ! Allocate all the most fundamental arrays, necessary for the computations
  call allocate_core

  ! need to reread the species into the now allocated arrays
  call components_input_species

  ! Output all the namelist items and output some other metadata
  call write_params

  ! calculate the proper value of beta and beta_prime 
  ! (called here in case the species beta_prime is required by miller geom)
  call components_beta_cal 

  ! Read the mode information
  ! this should be after allocate
  call mode_init

  ! Initialize the grids
  call geom_init_grids

  ! normalize the first for the global case
  ! FJC: suggest rhostar be moved to geom, and these lines to geom
  if (.not.flux_tube) then 
    do ix = 1, n_x_grid
      shift_end_grid(ix) = shift_end_grid(ix) / rhostar
    end do 
  endif 

  !Since chease read can modify q and shat, some checks need repeating
  !APS: ideally this hack should be removed.
  if(geom_type=='chease') then 
    call mode_check_params(2)
  endif

  ! if 2D array of modes is used one must call kgrid 
  ! Must be called before krbal. For a single mode one 
  ! must still call this routine since it normalises 
  ! the wave vectors used in the code 
  if (spectral_radius) then 
    call kgrid
  else 
    call kgrid_nonspec_rad
  endif

  ! Further allocate arrays
  call allocate_more

  ! Write geom quantities to file before parallelize_geom!
  ! Note these cannot be written with other write_run_params.
  if (.not. testing) call geom_output()
  if (.not. testing) call mode_box_recon

  !Must occur after components_init, and after geom_init_grids
  !Uses global tensors, so must be before parallelize_geom
  call rotation_init

  !Warning DO NOT change the order of this and the above routines
  call parallelize_geom
  call parallelize_rotation

  ! for the global case recalculate the profiles 
  call components_profiles 

  ! calculate the proper value of beta and beta_prime 
  ! (called again in case the chease inputs are required)
  call components_beta_cal

  ! calculate the parameters along the field line 
  call krbal

  !Allocates and initialises the velocity grid (uniform or not)
  !For vptrap case must be called after parallelize_geom
  call velgrid_init

  ! Allocate the arrays of diagnostic
  call diagnostic_check_params
  call diagnostic_allocate()

  ! set up the grids for the distribution function
  !call dist_grid_setup <! done by dist

  ! Initialise the array containing the Maxwellian
  call init_fmaxwl

  ! Initialise the gyro_average module
  if (.not.spectral_radius) then
    call gyro_average_allocate
    call gyro_average_init
  end if
  
  ! read from file (restart) if requested
  if (read_file .or. auto_restart) then
    call read_restart_file(irun)
  end if

  if (restarted) then
    call diagnostic_read_last_data
  else
    call abort_if_old_data_exists
  endif

  !The sizes of the real space box and location of the modes are calculated
  !This may be needed even if non_linear = false
  !This must be called after kgrid
  if (need_fft) call nonlinear_init
    
  ! calculate the linear terms, setup the matrix
  call init_linear_terms
  call calc_linear_terms
    
  ! if krook operator is required do the initialization (before init_explicit)
  if (nlkrook) call krook_init 
    
  ! Initialize the parameters for the explicit/implicit time step 
  select case(method) 
  case('EIV')
    call init_eiv
  case('EXP') 
    ! initialise the explict time integration and persistent communication
    ! needs only control, dist, mpicomms, and grid (could go higher up)
    call init_explicit        
  case('IMP') 
    ! currently imp_int called from gkw.f90 contains the initialisation.
    ! If the initialisation is split off, it should be called here
  case default 
    call gkw_exit('init: unknown integration scheme (method)')
  end select
      
  !The following case is for the study of tearing modes.  The magnetic island
  !is initialised as a perturbation in the parallel vector potential.
  call tearingmodes_init

  ! and initialise the distribution (but should not change the fields)
  ! in the cases where the g2f corretion is applied
  ! uses the linear terms / gyro_average matrices since it calls
  ! the fields solver, and in the cases also with parallel_x
  ! also uses the persistent_comm_init setup inside init_explicit
  call init_fdisi

  ! initialize normalise
  call normalise_init()

  ! set up the diagnostics (must be after nonlinear init)
  ! done last to reduce incidence of opening files then aborting
  call diagnostic_init

  ! not absolutely necessary, but
  ! the distribution is normalized
  call normalise_fdisi(fdisi,nsolc)
  
  ! If there is a bi-normal mode that should be persistent in time 
  ! reset it to its constant value.
  if(mode_persist) then
    call set_persistent_mode(fdisi)
    
    ! Calculate the fields with the correct persistent mode
    if (spectral_radius) then
      call calculate_fields(fdisi)
    else
      call field_solve_nonspec_wrap(fdisi,0.,.false.)
    endif
  end if
  
  
end subroutine initialize


!-----------------------------------------------------------------------------
!> Read the parameters from the input file. Call the corresponding routines
!> that do the parameter checking.
!-----------------------------------------------------------------------------
subroutine read_params()
  use version,         only : gkw_info
  use mpiinterface,    only : root_processor
  use control,         only : control_bcast_nml, control_read_nml
  use control,         only : control_check_params
  use grid,            only : grid_read_nml
  use grid,            only : grid_bcast_nml, grid_check_params
  use diagnostic,      only : diagnostic_read_nml
  use diagnostic,      only : diagnostic_bcast_nml
  use rotation,        only : rotation_read_nml
  use rotation,        only : rotation_bcast_nml, rotation_check_params
  use geom,            only : geom_read_nml
  use geom,            only : geom_bcast_nml, geom_check_params
  use mode,            only : mode_read_nml
  use mode,            only : mode_bcast_nml, mode_check_params
  use components,      only : components_read_nml_spcg
  use components,      only : components_bcast_nml_spcg
  use components,      only : components_check_params_spcg
  use components,      only : components_read_nml_spec  
  use components,      only : components_bcast_nml_spec
  use components,      only : components_check_params_spec
  use components,      only : n_spec
  use collisionop,     only : collisionop_read_nml
  use collisionop,     only : collisionop_bcast_nml, collisionop_check_params
  use linear_terms,    only : linear_terms_read_nml, linear_terms_bcast_nml
  use linear_terms,    only : linear_terms_check_params
  use krook,           only : krook_read_nml, krook_bcast_nml
  use gyro_average,    only : gyro_average_read_nml, gyro_average_bcast_nml 
  use source_time,     only : source_time_read_nml, source_time_bcast_nml
  use source_time,     only : source_time_check_params
  use eiv_integration, only : eiv_integration_read_nml
  use eiv_integration, only : eiv_integration_bcast_nml
  use eiv_integration, only : eiv_integration_check_params
  use rho_par_switch,  only : rho_par_switch_read_nml
  use rho_par_switch,  only : rho_par_switch_bcast_nml 
  use rho_par_switch,  only : rho_par_switch_check_params 

  integer, parameter :: file_unit = 91
  
  integer :: io_stat, i
  logical :: input_exists

  ! exit for now when no input file exists
  inquire(FILE='input.dat',EXIST=input_exists)
  if (.not. input_exists) call gkw_info()

  ! Read and check the various namelists
  ! Set io_stat to zero; only root_processor will obtain a different value.
  io_stat = 0
  
  ! The order here matters for the checks, and should follow the
  ! dependency hierarcy

  open(file_unit,file='input.dat',FORM='formatted',STATUS='old', &
         POSITION='rewind', ACTION='read')
  ! control
  if (root_processor) then
    rewind(file_unit)
    call control_read_nml(file_unit, io_stat, .false.)
  end if
  call namelist_error_check('control',io_stat)
  call control_bcast_nml
  call control_check_params

  ! grid
  if (root_processor) then
    rewind(file_unit)
    call grid_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('gridsize',io_stat)
  call grid_bcast_nml
  call grid_check_params

  ! geom
  if (root_processor) then
    rewind(file_unit)
    call geom_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('geom',io_stat)
  call geom_bcast_nml
  call geom_check_params(1) !may be called again later

  ! mode
  if (root_processor) then
    rewind(file_unit)
    call mode_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('mode',io_stat)
  call mode_bcast_nml
  call mode_check_params(1) !may be called again later

  ! components spcgeneral
  if (root_processor) then
    rewind(file_unit)
    call components_read_nml_spcg(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('spcgeneral',io_stat)
  call components_bcast_nml_spcg
  call components_check_params_spcg
  
  ! components species
  ! *** The species can be read for checking, but not initialised till after
  ! *** allocate.
  do i=1,n_spec
    ! Rewind the file when the first species namelist is read.
    if (root_processor) then
      if (i == 1) rewind(file_unit)
      call components_read_nml_spec(file_unit, io_stat)
    end if
    call namelist_error_check('species',io_stat,i)
    call components_bcast_nml_spec
    call components_check_params_spec
  end do

  ! rotation
  if (root_processor) then
    rewind(file_unit)
    call rotation_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('rotation',io_stat)
  call rotation_bcast_nml
  call rotation_check_params

  ! eiv_integration (optional)
  if (root_processor) then
    rewind(file_unit)
    call eiv_integration_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('eiv_integration',io_stat)
  call eiv_integration_bcast_nml
  call eiv_integration_check_params

  ! linear terms (optional)
  if (root_processor) then
    rewind(file_unit)
    call linear_terms_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('linear_terms',io_stat)
  call linear_terms_bcast_nml
  call linear_terms_check_params

  ! collisions (optional)
  if (root_processor) then
    rewind(file_unit)
    call collisionop_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('collisions',io_stat)
  call collisionop_bcast_nml
  call collisionop_check_params

  ! krook (optional)
  if (root_processor) then
    rewind(file_unit)
    call krook_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('krook',io_stat)
  call krook_bcast_nml

  ! gyro_average (optional / nonspectral only)
  if (root_processor) then
    rewind(file_unit)
    call gyro_average_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('gyroaverage',io_stat)
  call gyro_average_bcast_nml

  ! source time (optional)
  if (root_processor) then
    rewind(file_unit)
    call source_time_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('source_time',io_stat)
  call source_time_bcast_nml
  call source_time_check_params

  ! source time (optional)
  if (root_processor) then
    rewind(file_unit)
    call rho_par_switch_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('finite_rho_parallel',io_stat)
  call rho_par_switch_bcast_nml
  call rho_par_switch_check_params

  ! diagnostic (optional)
  if (root_processor) then
    rewind(file_unit)
    call diagnostic_read_nml(file_unit,io_stat, .false.)
  end if
  call namelist_error_check('diagnostic',io_stat)
  call diagnostic_bcast_nml
  ! check of diagnostic namelist parameters is done later, because at
  ! this stage many other things are not initialised yet and
  ! cannot be tested for.
  
  close(file_unit)
  
  ! The parameters which were read from input.dat into namelists are metadata.
  ! We want to attach them also to the output data, be it a simple ascii file
  ! or HDF5 attributes or other.
  ! Therefore the IO interfaces must be initialized before those metadata can
  ! be output. The IO interfaces can be initialized only *after* the input.dat
  ! has been read, because they must know the desired format and whether the
  ! run is a restarted one or not.

end subroutine read_params


!-----------------------------------------------------------------------------
!> Write the parameters in the desired output
!> format, e.g. to the file input.out .
!-----------------------------------------------------------------------------
subroutine write_params()
  use mpiinterface, only : root_processor, mpibarrier
  use version,         only : output_header
  use control,         only : control_write_nml
  use grid,            only : grid_write_nml
  use diagnostic,      only : diagnostic_write_nml
  use rotation,        only : rotation_write_nml
  use geom,            only : geom_write_nml
  use mode,            only : mode_write_nml
  use components,      only : components_write_nml_spcg
  use components,      only : n_spec, components_write_nml_spec
  use collisionop,     only : collisionop_write_nml
  use linear_terms,    only : linear_terms_write_nml
  use krook,           only : krook_write_nml
  use gyro_average,    only : gyro_average_write_nml
  use source_time,     only : source_time_write_nml
  use eiv_integration, only : eiv_integration_write_nml
  use rho_par_switch,  only : rho_par_switch_write_nml
  
  integer, parameter :: out_file_unit = 92
  integer :: io_stat, i
  
  if (root_processor) then
    open(UNIT=out_file_unit, FILE='input.out', FORM='formatted', &
       & STATUS='replace', POSITION='rewind')
  end if
  
  ! write the header for the input.out file
  call output_header(out_file_unit)

  call control_write_nml(out_file_unit,io_stat)
  call grid_write_nml(out_file_unit,io_stat)
  call diagnostic_write_nml(out_file_unit,io_stat)
  call geom_write_nml(out_file_unit,io_stat)
  call mode_write_nml(out_file_unit,io_stat)
  call components_write_nml_spcg(out_file_unit,io_stat)
  do i=1,n_spec
    call components_write_nml_spec(out_file_unit, i, io_stat)
  end do
  call rotation_write_nml(out_file_unit,io_stat)
  call eiv_integration_write_nml(out_file_unit,io_stat)
  call linear_terms_write_nml(out_file_unit,io_stat)
  call collisionop_write_nml(out_file_unit,io_stat)
  call krook_write_nml(out_file_unit,io_stat)
  call gyro_average_write_nml(out_file_unit,io_stat)
  call source_time_write_nml(out_file_unit,io_stat)
  call rho_par_switch_write_nml(out_file_unit,io_stat)
  if (root_processor) then
    close(out_file_unit)
  end if
  ! (formerly: don't want to use input.out till it is written) Now,
  ! input.out is not to be read in anymore, but input.dat is read
  ! twice.
  call mpibarrier()

  
end subroutine write_params

! !-----------------------------------------------------------------------------
! !> processor the namelist given in read_namelist_from_input
! !> FJC: This routine wrapper should be removed, it does not add anything
! !> except confusion.The readwrite namelist routines should be called directly 
! !> with no interface since the read/write argument has not need to be optional.
! !> If necessary, they can be renamed to module_readwrite_namelist
! !-----------------------------------------------------------------------------
! subroutine read_nml(read_namelist_from_input,file_unit,io_stat,i_spec)
  
!   use control,       only : control_read_nml
!   use grid,          only : grid_read_nml
!   use diagnostic,    only : diagnostic_read_nml
!   use rotation,      only : rotation_read_nml
!   use geom,          only : geom_read_nml
!   use mode,          only : mode_read_nml
!   use components,    only : components_read_nml_spcg
!   use components,    only : components_read_nml_spec
!   use collisionop,   only : collisionop_read_nml
  
!   interface
!     subroutine read_namelist_from_input(i_lun,i_stat,l_write)
!       integer, intent(in)  :: i_lun
!       integer, intent(out) :: i_stat
!       logical, optional, intent(in)  :: l_write
!     end subroutine read_namelist_from_input
!   end interface

!   integer, intent(in)           :: file_unit
!   integer, intent(out)          :: io_stat
!   integer, optional, intent(in) :: i_spec

!   logical, parameter :: l_write = .false.

!   if (present(i_spec)) then
!     ! The input file may contain several 'species' namelists.
!     ! In order to read them consecutively, the file has to stay open while they
!     ! are read.
!     if (i_spec == 1) then
!       ! Open the file anew when the first species namelist is read.
! #error      open(file_unit,file='input.dat',FORM='formatted',STATUS='old')
!     end if
!   else
!     ! By rewinding the file pointer before any other namelist is read,
!     ! we allow for an arbitrary order of the namelists.
!     open(file_unit,file='input.dat',FORM='formatted',STATUS='old', &
!        & POSITION='rewind')
!   end if

!   ! call the read_input routine in the appropriate module
!   call read_namelist_from_input(file_unit,io_stat,l_write)
  
!   ! Close the file, unless we have more species to read
!   if (present(i_spec)) then
!     if (i_spec == n_spec) close(file_unit)
!   else
!     close(file_unit)
!   end if
    
! end subroutine read_nml


!-----------------------------------------------------------------------------
!> Check io_stat for a read then abort if any non-zero values are found.
!-----------------------------------------------------------------------------
subroutine namelist_error_check(list_name,io_stat,ispc)

  use mpiinterface, only : root_processor, mpibcast
  use general,      only : gkw_exit

  character (len=*), intent(in) :: list_name
  integer, intent(in)           :: io_stat
  integer, optional, intent(in) :: ispc

  integer :: my_io_stat

  my_io_stat = 0

  ! root_processor first reports any error
  if (root_processor) then
    my_io_stat = io_stat

    if (io_stat > 0) then
      write(*,*) '* Error in namelist '//list_name//'.'
      if (present(ispc)) write(*,*) '* for species number: ',ispc
    end if

    if (io_stat < 0) then
      if (check_if_optional_namelist(list_name)) then
        write(*,*) '* Optional namelist '//list_name//' MISSING.'
        ! For optional namelists, ignore missing status
        my_io_stat = 0
      else
        write(*,*) '* Namelist '//list_name//' MISSING.'
        if (present(ispc)) write(*,*) '* for species number: ',ispc
      end if
    end if

  end if

  ! broadcast io_stat from root_processor
  call mpibcast(my_io_stat,1)

  ! abort on non-zero io_stat
  if (my_io_stat /= 0) then
    call gkw_exit('namelist_error_check: problem with '//list_name)
  end if

end subroutine namelist_error_check

!-----------------------------------------------------------------------------
!> Function that returns true if the name given is that of an optional
!> namelist and false otherwise, i.e. also if the name given is no namelist at
!> all.
!-----------------------------------------------------------------------------
function check_if_optional_namelist(list_name)
  character (len=*), intent(in) :: list_name
  logical :: check_if_optional_namelist

  if (     list_name == 'collisions'  &
    & .or. list_name == 'diagnostic' &
    & .or. list_name == 'eiv_integration' &
    & .or. list_name == 'finite_rho_parallel' &
    & .or. list_name == 'gyroaverage' &
    & .or. list_name == 'krook' &
    & .or. list_name == 'linear_terms' &
    & .or. list_name == 'rotation' &
    & .or. list_name == 'source_time' ) then
    check_if_optional_namelist = .true.
  else
    check_if_optional_namelist = .false.
  end if
end function check_if_optional_namelist


!-----------------------------------------------------------------------------
!> Top level initialization routine that selects the initialization method
!> which is then done by the soubroutine init_fdis.
!> It allows for the initialization of a specific mode set by imod_init. 
!> Furthermore, it sets up the persisten_mode (for mode_persist = .true.)
!> that is set constant during the time integration. 
!-----------------------------------------------------------------------------
subroutine init_fdisi
  use components, only : amp_init, amp_init_imag
  use mpiinterface, only : root_processor
  use dist, only : fdisi, ifdis, fdis_tmp, nf
  use fields, only : get_averaged_apar
  use fields, only : calculate_fields, field_solve_nonspec_wrap
  use components, only : finit, finit_imod
  use components, only : amp_imod, amp_imod_imag, amp_zon, amp_zon_imag
  use components, only : imod_init, ind_persist, mode_persist
  use mode,       only : iyzero
  use exp_integration, only : persistent_mode
  use grid,            only : nmod,ns,nx,nvpar,nmu,nsp
  use index_function,  only : indx
  use components,      only : init_coef
  
  integer :: ix, i, j, k, is, iiii, imod

  if(root_processor) then
    write(*,*) 'For initialization of the distribution funciton use amplitude: ', amp_init, amp_init_imag
  end if
  
  ! in case of restart and imod_init > 0, backup fdisi in fdis_tmp 
  if(imod_init > 0) then
    fdis_tmp(1:nf) = fdisi(1:nf)
  end if
  
  ! set the initial conditions for fdisi 
  call init_fdis(fdisi, amp_init, amp_init_imag, finit)

  ! set the initial conditions for fdis_tmp that is used as buffer to treat
  ! the specific mode set by imod_init and persisten mode later 
  if(imod_init > 0) then
  
    if(root_processor) then
      write(*,*) 'For initialization of the mode specified by imod_init use amplitude: ', amp_imod, amp_imod_imag
    end if
	  
	  ! initialize fdis_tmp
    call init_fdis(fdis_tmp, amp_imod, amp_imod_imag, finit_imod)
	
    ! now combine fdis_tmp into fdisi by copying mode imod_init from fdis_tmp to fdisi:
    do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
      fdisi(indx(ifdis,imod_init,ix,i,j,k,is)) = fdis_tmp(indx(ifdis,imod_init,ix,i,j,k,is))
    end do; end do; end do; end do; end do
		
  end if
  
	
  ! scale the zonal mode according to amp_zon
  if(iyzero > 0) then
  
    if(root_processor) then
      write(*,*) 'For scaling of the zonal mode use factor: ', amp_zon, amp_zon_imag
    end if
    
    do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
      fdisi(indx(ifdis,iyzero,ix,i,j,k,is)) = cmplx(amp_zon, amp_zon_imag) * &
	    & fdisi(indx(ifdis,iyzero,ix,i,j,k,is))
    end do; end do; end do; end do; end do
  end if


  ! now backup the mode imod_init according to mode_persist
  if(mode_persist .and. imod_init > 0) then
  
	  ! Allocate space for a copy of the modes of the distribution function 
	  ! that should be temporally persistent and the corresponding index array 
	  ! for copying it into fdisi later on
    allocate(ind_persist(ns*nx*nvpar*nmu*nsp))
    allocate(persistent_mode(ns*nx*nvpar*nmu*nsp))
	  
	  ! set up the persistent mode
    iiii=0
    do imod=1,nmod; do i=1,ns; do j=1,nmu; do k=1,nvpar; do is=1,nsp; do ix=1,nx
      
      ! save mode specified by imod_init in persistent_mode and 
      ! save corresponding index
      if(imod==imod_init) then
        iiii=iiii+1
        ind_persist(iiii)=indx(ifdis,imod,ix,i,j,k,is)
        persistent_mode(iiii) = fdisi(ind_persist(iiii))
      end if
	    
    end do ; end do ; end do  ; end do ; end do ; end do
    
  end if
  
end subroutine init_fdisi



!-----------------------------------------------------------------------------
!> Subroutine that initialises the distribution function with a 
!> choice of intial conditions, and the maxwellian
!> should not change the fields !
!> Note: fdisi contains g, but since the A|| perturbation is intially
!> zero, initialising f is the same as g (except for tearingmodes)
!> GLOBAL: Needs a 1/dgrid for all the fdisi init options
!> and some tgrids also
!-----------------------------------------------------------------------------
subroutine init_fdis(fdis, amp_init_real, amp_init_imag, finit)
  use marsaglia,      only : ran_kiss
  use control,        only : spectral_radius, fullf_wo_Fm
  use restart,        only : restarted
  use grid,           only : nmod, nx, ns, nsp, nmu, nvpar, n_x_grid, lx
  use grid,           only : n_vpar_grid, n_mu_grid, n_s_grid, nxmod
  use grid,           only : ls, lmu, lvpar, lsp, lrx, gx, number_of_species
  use grid,           only : parallel_vpar, parallel_mu
  use dist,           only : ifdis, nelem_conserve, fmaxwl, n_conserve
  use dist,           only : f_EP
  use fields,         only : get_averaged_apar, calculate_fields 
  use fields,         only : field_solve_nonspec_wrap
  use index_function, only : indx
  use geom,           only : bn, sgr, kthnorm, dxgr
  use mode,           only : mode_box, ixzero, iyzero, ikxspace 
  use mode,           only : kxrh, krho
  use components,     only : quench_modes, rhostar, types 
  use components,     only : signz, de, kwid_ini
  use components,     only : vthrat, isl_ls, isl_mode, tmp, tearingmode, lfinit_radial_dirichlet
  use components,     only : quench_switch, n_quench
  use velocitygrid,   only : vpgr, mugr, intvp, intmu
  use constants,      only : ci1, pi 
  use mpiinterface,   only : root_processor, mpiallreduce_sum_inplace
  use general,        only : gkw_abort
  use global,         only : lenswitch
  use mpicomms,       only : COMM_VPAR_NE_MU_NE
  use components,     only : init_coef, tgrid
  use general,        only : gkw_warn
  
  complex, intent(out) :: fdis(:)
  real, intent(in) :: amp_init_real, amp_init_imag
  character(len=lenswitch), intent(in) :: finit

  ! integers for the loop over all grid points 
  integer :: imod, ix ,i, j, k, is, line, i_s, j_mu, k_vpar, i_sp, i_x, idx
  
 
  ! Dummy variables 
  real    :: random1, random2, vperp, random3, random4, a
  complex :: amp_ini, apar_ga
  integer :: idum, ikxspace_local, imodisland 
  integer :: iplusmode,iminusmode, en = 0, p, modnum
  real    :: balloonpos, dumbuf, dum
  character(len=lenswitch) :: finit_restart
  
   ! Variables for cos-zonal-tmp init case
  real :: intfm, intfm1, intfm2, intv2, v2norm
  real :: n_zonal, dens_phase, amp_tmp, tmp_phase, no_tmp
  complex :: tmp_norm
  
  n_zonal = init_coef(1)
  dens_phase = init_coef(2)
  amp_tmp = init_coef(3)
  tmp_phase = init_coef(4)
  no_tmp = init_coef(5)
  

  ! initialise
  k = -999 ; is = -999

  !Initialise the complex initial amplitude
  amp_ini=cmplx(amp_init_real,amp_init_imag)  
  
  !Initialise collisions conservation fields to zero 
  !(probably not needed)
  !(could this field be real rather than complex?)
  if (nelem_conserve > 0) then 
    fdis(n_conserve+1:n_conserve+nelem_conserve) = (0.E0,0.E0)
  end if
    
  ! There may be ways to initialise the distribution which first read
  ! a restart file and then modify it.  For those, it is checked
  ! whether this run is indeed a restarted one.
  if (( index(finit,'restart+').EQ.1 ) .and. .not. restarted) then
    call gkw_abort("The initial conditions '"//trim(finit)//"' can only be used when restarting a run.")
  end if
  

  if ( restarted ) then
  
    ! if a restart-initial-condition is desired:
    if ( index(finit,'restart+').EQ.1 ) then
      ! These initial conditions differs from the ordinary set below:
      ! They assumes that the run is restarted and fdis is already filled.
      
      finit_restart = trim(finit(scan(finit,'+')+1:))
      
      if(root_processor) then
        write(*,*)
        write(*,*) '------------------------------------------------------'
        write(*,*)"'restart+' chosen as finit, altering distribution function.&
        & Using '" // trim(finit_restart) // "' to alter."
      end if
      
      select case(finit_restart)
      
      case('quench_ky')
        ! With this initial condition the zonal flow eigenmode is set to zero.
        
        ! For electromagnetic runs, A_par has to be added to f in fdis. For this, f is filled into fdis,
        ! then A_par is calculated with it and added to it.
        
        do idx = 1, size(quench_modes)
          if(1 <= quench_modes(idx) .and. quench_modes(idx) <= nmod) then
            do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar;
              do is = 1, nsp
                ! shifted cosine prevents radial boundary artefacts
                fdis(indx(ifdis,quench_modes(idx),ix,i,j,k,is)) = amp_ini * &
                & (cos(2*pi*sgr(i))+1.0) 
              end do
            end do; end do; end do; end do;
          end if
        end do
      
      case('quench_turbulence')
        ! With this initial condition the turbulent modes can be set to a small value,
        ! only the zonal flow is kept.
        
        do imod = 1, nmod
          if (imod.ne.iyzero) then
						! shifted cosine prevents radial boundary artefacts
            do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
              fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini * (cos(2*pi*sgr(i))+1.0) 
            end do; end do; end do; end do; end do
          end if
        end do
        
        
      case('quench_turb_and_some_ZF')
        ! With this initial condition the turbulent modes can be set to a small value and
        ! furthermore only a part of the zonal flow is kept.
        
        if ( n_quench < 0 .or. n_quench > nx/2 ) then
          call gkw_abort('Cannot quench a nonexistent mode, n_quench has to be in [0,floor(nx/2)].')
        end if

        do imod = 1, nmod
        
          if (imod.ne.iyzero) then
            do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp 
              ! shifted cosine prevents radial boundary artefacts
              fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini * (cos(2*pi*sgr(i))+1.0) 
            end do; end do; end do; end do; end do
            
          else
            
            select case (quench_switch)
            
            case('just')
              
              do ix = 1, nx
                if ( ix.eq.nxmod-n_quench .or. ix.eq.nxmod+n_quench) then
                  do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
                    fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini
                  end do; end do; end do; end do
                endif
              end do
          
            case('except')
            
              do ix = 1, nx
                if ( ix.ne.nxmod-n_quench .and. ix.ne.nxmod+n_quench) then
                  do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
                    fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini
                  end do; end do; end do; end do
                endif
              end do
            
            case('above')
            
              do ix = 1, nx
                if ( ix.lt.nxmod-n_quench .or. ix.gt.nxmod+n_quench) then
                  do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
                    fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini 
                  end do; end do; end do; end do
                endif
              end do
            
            case('below')
            
              do ix = 1, nx
                if ( ix.gt.nxmod-n_quench .and. ix.lt.nxmod+n_quench) then
                  do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
                    fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini
                  end do; end do; end do; end do
                endif
              end do
          
            case default
            ! abort if an unknown condition is specified  
              call gkw_abort('Unknown quench_switch: '//trim(quench_switch))
              
            end select
            
          end if
            
        end do
        
        if(root_processor) then
          write(*,*) 'Use amplitude: ', amp_ini, 'to quench ', trim(quench_switch), n_quench
        end if
        
        
      case('quench_all_except_SI')
      
        ! With this initial condition the turbulent modes can be set to a small value, 
        ! only the self-interacting modes in the zonal flow are kept.
        
        do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
          if (imod.ne.iyzero) then
            do ix = 1, nx
              ! shifted cosine prevents radial boundary artefacts
              fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini * (cos(2*pi*sgr(i))+1.0) 
            end do
          else
            do ix = 1, nx  
              if ( modulo(ix-(nx/2)-1,ikxspace).ne.0) then
                fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini
              endif
            end do
            fdis(indx(ifdis,imod,nxmod,i,j,k,is)) = amp_ini
          end if
        end do; end do; end do; end do; end do
        
        
        if(root_processor) then
          write(*,*) 'Use amplitude: ', amp_ini
        end if
        
      case('scale_turbulence')
        ! With this initial condition the turbulent modes can be scaled to a small value,
        ! only the zonal flow is kept.
        
        do imod = 1, nmod; do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp
            if (imod.ne.iyzero) then
              fdis(indx(ifdis,imod,ix,i,j,k,is)) = fdis(indx(ifdis,imod,ix,i,j,k,is)) * amp_ini
            end if
        end do; end do; end do; end do; end do; end do
        
        if(root_processor) then
          write(*,*) 'Use scaling-factor: ', amp_ini
        end if
        
      case('scale_zonal_mode')
        ! With this initial condition individual radial modes of the 
        ! zonal mode (ky=0 mode) are scaled.
        ! The ky/=0 modes are kept.
        
        do idx = 1, size(quench_modes)
          if(1 <= quench_modes(idx) .and. quench_modes(idx) <= (nx-1)/2) then
            do ix = 1, nx; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar;
              do is = 1, nsp
                if(ix == ixzero+quench_modes(idx)) then
                  fdis(indx(ifdis,iyzero,ix,i,j,k,is)) &
                  & = fdis(indx(ifdis,iyzero,ix,i,j,k,is)) * amp_ini
                end if
                if(ix == ixzero-quench_modes(idx)) then
                  fdis(indx(ifdis,iyzero,ix,i,j,k,is)) &
                  & = fdis(indx(ifdis,iyzero,ix,i,j,k,is)) * amp_ini
                end if
              end do
            end do; end do; end do; end do;
          end if
        end do
      
      case default
        ! abort if an unknown condition is specified
        
        call gkw_abort('Unknown initilization switch finit: '//trim(finit))
        
      end select
      
      if(root_processor) then
        write(*,*) '------------------------------------------------------'
        write(*,*)
      end if
      
      
    endif
    
    ! if no alteration is desired, do nothing
    
    
  else
  
  
    ! The following initial conditions are applicable for new runs,
    ! i.e. not restarted ones.
    select case(finit)
    case('noise','gnoise') ! grid scale uniform (or gaussian) noise

      do imod = 1, nmod; do i_s = 1, n_s_grid; do j_mu = 1, n_mu_grid 
        do k_vpar = 1, n_vpar_grid; do i_sp = 1, number_of_species
          do i_x = 1, n_x_grid

            ! always generate 2 numbers on interval [0 1]
            random1 = ran_kiss() ; random2 = ran_kiss()

            ! apply box-muller transformation to make gaussian
            if (finit=='gnoise') then 
              do while (random1 < 1.e-8 .or. random1 > 1.-1.e-8)
                random1 = ran_kiss()
              end do
              random3 = sqrt(-2.*log(random1))*cos(2*pi*random2)
              random4 = sqrt(-2.*log(random1))*sin(2*pi*random2)
              random1 = random3; random2=random4
            else ! shift onto [-1 1]
              random1 = 1.-2.*random1
              random2 = 1.-2.*random2
            end if

            ! get the local index
            i= ls(i_s); j= lmu(j_mu); k=lvpar(k_vpar)
            is= lsp(i_sp); ix=lrx(i_x)

            ! only initialise if we are on the local proc
            if (i >= 1 .and. i <= ns .and. j >= 1 .and. j <= nmu &
               & .and. k >= 1 .and. k <= nvpar .and. is >= 1    &
               & .and. is <= nsp .and. ix >= 1 .and. ix <= nx) then

              ! No initialisation of the zonal flow 
              if (imod.eq.iyzero) then
                fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  
              else 
                fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini *   &
                   &            (random1 + random2*ci1)
              end if

            end if

          end do
        end do; end do
      end do; end do; end do 

    case('hybrid') ! gaussian noise shaped towards something physical

      ! modal peak in ky maxwellian
      a=sqrt(2.)*kwid_ini/kthnorm

      do imod = 1, nmod; do i_s = 1, n_s_grid; do j_mu = 1, n_mu_grid
        do k_vpar = 1, n_vpar_grid; do i_sp = 1, number_of_species
          do i_x = 1, n_x_grid

            ! always generate 2 numbers on interval [0 1]
            random1 = ran_kiss() ; random2 = ran_kiss()

            ! apply box-muller transformation to make gaussian
            do while (random1 < 1.e-8 .or. random1 > 1.-1.e-8)
              random1 = ran_kiss()
            end do
            random3 = sqrt(-2.*log(random1))*cos(2*pi*random2)
            random4 = sqrt(-2.*log(random1))*sin(2*pi*random2)
            random1 = random3; random2=random4

            ! get the local index
            i= ls(i_s); j= lmu(j_mu); k=lvpar(k_vpar)
            is= lsp(i_sp); ix=lrx(i_x)

            ! only initialise if we are on the local proc
            if (i >= 1 .and. i <= ns .and. j >= 1 .and. j <= nmu &
               & .and. k >= 1 .and. k <= nvpar .and. is >= 1    &
               & .and. is <= nsp .and. ix >= 1 .and. ix <= nx) then

              ! No initialisation of the zonal flow 
              if (imod.eq.iyzero) then
                fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  
              else 
                fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini *        &                                              
                                ! random imaginary number, uniform on [-1,1]
                   (random1 + random2*ci1) *                            &
                                !ballooning, gaussian centered on kx=0
                   cos(2*pi*sgr(i))*exp(-kxrh(ix)*kxrh(ix)/(20.*a*a))*  &
                                ! maxwellian in velocity space
                   fmaxwl(ix,i,j,k,is) *                                &
                                ! maxwellian in ky, with a modal peak at sqrt(2)*a
                   exp(-krho(imod)*krho(imod)/(2*a*a)) *                &
                   krho(imod)*krho(imod)/(a*a*a)
              end if

            end if

          end do
        end do; end do
      end do; end do; end do
      
      
      
      
      

    case('cosine6')

      a=sqrt(2.)*kwid_ini/kthnorm

      ! Maxwellian in velocity space
      ! Ballooning and gaussian centered about kx = 0
      ! Like hybrid but without the random element or ky dependence - for MTM
      ! Like cosine3 and 4, but with better kx dependence
      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx
          ! No initialisation of the zonal flow 
          if (imod == iyzero .and. mode_box .and. nmod > 1) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)
          ! or radial streamers
          elseif ((ix == ixzero).and.mode_box) then          
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)
          else
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini                      &
               & * (cos(2*pi*sgr(i)) + 1.0)*exp(-kxrh(ix)*kxrh(ix)/(20.*a*a))  &
               & * exp(-(vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j)))
          endif
        end do; end do
      end do; end do; end do; end do

    case('cosine5')

      a=sqrt(2.)*kwid_ini/kthnorm

      ! Maxwellian in velocity space
      ! Ballooning and gaussian centered about kx = 0, excluding ixzero modes
      ! Like hybrid but without the random element or ky dependence - for MTM
      ! Like cosine3 and 4, but with better kx dependence
      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx
          ! No initialisation of the zonal flow 
          if (imod == iyzero .and. mode_box .and. nmod > 1) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)
          else
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini                      &
               & * (cos(2*pi*sgr(i)) + 1.0)*exp(-kxrh(ix)*kxrh(ix)/(20.*a*a))  &
               & * exp(-(vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j)))
          endif
        end do; end do
      end do; end do; end do; end do

    case('cosine4')
      ! This makes sense particularly with the spectral method: The
      ! initial condition is bell-shaped in velocity space and has in
      ! position space a cosine-shaped hill at the center of the
      ! field line, i.e. at the central radial wavenumber k_psi=0
      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx
          ! No initialisation of the zonal flow 
          if ((imod.eq.1).and.mode_box.and.(nmod.ne.1)) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)
          elseif(ix == ixzero) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini              &
               & * (cos(2*pi*sgr(i)) + 1.0)                             &
               & * exp(-(vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j)))
          else
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  
          endif
        end do; end do
      end do; end do; end do; end do

    case('cosine3')
      ! resembles the option cosine2 but is bell-shaped in velocity space:
      ! the distribution function becomes smaller towards the borders
      ! of the velocity space

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx 

          ! No initialisation of the zonal flow 
          if (imod == iyzero .and. mode_box .and. nmod > 1) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  
          else 
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini              &
               & * (cos(2*pi*sgr(i)) + 1.0)                            &
               & * exp(-(vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j)))
          endif

        end do; end do
      end do; end do; end do; end do

    case('cosine2') !default 

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx

          ! No initialisation of the zonal flow
          if (imod == iyzero .and. mode_box .and. nmod > 1) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)
          else
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini*(cos(2*pi*sgr(i))+1.0) 
          endif

        end do; end do
      end do; end do; end do; end do

    case('cosine') !deprecated - see issue 31.

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx 

          ! No initialisation of the zonal flow 
          if (imod == iyzero) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  
          else 
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini * cos(2*pi*sgr(i)) 
          endif

        end do; end do
      end do; end do; end do; end do

    case('sine')  

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu;  do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx 

          ! No initialisation of the zonal flow 
          if ((imod.eq.1).and.mode_box.and.(nmod.ne.1)) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  
          else 
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini*de(ix,is)    &
               &      *(sin(2*pi*sgr(i))+1.0) 
          endif

        end do; end do
      end do; end do; end do; end do
      
      
    case('cos-zonal-tmp')
    
      ! Initialize sinusoidal (with respect to radial coordinate) zonal 
      ! density and temperature perturbation, both having radial wave number n_zonal.
      ! The amplitude of the density perturbation is controlled by amp_init_real
      ! and the phase of the sinusoidal perturbation is controlled by dens_phase.
      ! The amplitude of the temperature perturbation is controlled by amp_tmp
      ! and the phas is set by tmp_phase.
      ! The turbulent modes are initialized to zero.
      
      if(spectral_radius) then
        call gkw_warn('Note: The phase of the zonal density and temperature &
        & perturbation is ignored in spectral runs.')
      endif
    
      do ix = 1, nx; do i = 1, ns; do is = 1, nsp
      
        ! initialize the integral to zero
        intfm1 = 0.0
        intfm2 = 0.0
        intv2 = 0.0 
        
        ! Perform the integral in velocity space
        ! necessary for normalisation, such that only a temperature
        ! perturbation is generated and no density or momentum.
        do j=1, nmu; do k=1, nvpar
        
          intfm1 = intfm1 &
               & +(vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j))/(tmp(ix,is)/tgrid(is))  &
               & *fmaxwl(ix,i,j,k,is) &
               & *intvp(i,j,k,is)*intmu(j)*bn(ix,i) 
        
          intfm2 = intfm2 + fmaxwl(ix,i,j,k,is)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)
          
          intv2 = intv2 + &
               & (vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j))/(tmp(ix,is)/tgrid(is))  &
               & *intvp(i,j,k,is)*intmu(j)*bn(ix,i) 
        
        end do; end do
      
      
        ! finish velocity space integral
         if (parallel_vpar .or. parallel_mu) then 
           call mpiallreduce_sum_inplace(intfm1,1,COMM_VPAR_NE_MU_NE)
           call mpiallreduce_sum_inplace(intfm2,1,COMM_VPAR_NE_MU_NE)
           call mpiallreduce_sum_inplace(intv2,1,COMM_VPAR_NE_MU_NE)
         endif 
         
        ! normalize the sum of the energy moment to the sum of the Maxwell such that 
        ! the particle source is exactly zero 
        intfm = intfm1 / intfm2
        
        ! normalize such that initial density perturbation does not produce a
        ! temperature perturbation
        v2norm = intfm1 / intv2
      

        ! initialization of the distribution function
        do imod = 1, nmod; do j = 1, nmu; do k = 1, nvpar
        
          ! initialize the zonal
          if (imod.eq.iyzero) then
          
            ! factor that ensures temperature perturbation only
            tmp_norm = ((vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j))/(tmp(ix,is)/tgrid(is))-intfm)
           
            if (spectral_radius) then
        
              if (ixzero.eq.0) call gkw_abort('No kx = 0 mode found in init zonal')
              if ( (ixzero.le.1) .or. (ixzero.ge.nx) ) &
                  & call gkw_abort('No finite kx mode found in init')
              
                fdis(indx(ifdis,imod,ixzero-int(n_zonal),i,j,k,is)) = -ci1 * &
                  & ((fmaxwl(ixzero-int(n_zonal),i,j,k,is) - no_tmp * v2norm) * amp_init_real + &
                  & fmaxwl(ixzero-int(n_zonal),i,j,k,is) * tmp_norm * amp_tmp)/2.
                fdis(indx(ifdis,imod,ixzero+int(n_zonal),i,j,k,is)) = +ci1 * &
                  & ((fmaxwl(ixzero+int(n_zonal),i,j,k,is) - no_tmp * v2norm) * amp_init_real + &
                  & fmaxwl(ixzero+int(n_zonal),i,j,k,is) * tmp_norm * amp_tmp)/2.
              else
                fdis(indx(ifdis,imod,ix,i,j,k,is)) = (fmaxwl(ix,i,j,k,is) - no_tmp * v2norm) * &
                  & amp_init_real * sin(2*pi*n_zonal*gx(ix)/n_x_grid + dens_phase * pi) + &
                  & fmaxwl(ix,i,j,k,is) * tmp_norm * amp_tmp * sin(2*pi*n_zonal*gx(ix)/n_x_grid + tmp_phase * pi)
              endif
            else
              ! initialize turbulent modes to zero
              fdis(indx(ifdis,imod,ix,i,j,k,is)) = 0.0
          endif

        end do; end do; end do
      end do; end do; end do



    case('gam')

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx

          ! Initialisation of the mode (m,n)=(1,0)
          if ((imod.eq.1).and.mode_box.and.(nmod.ne.1)) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini*de(ix,is)    &
               &      *(sin(sgr(i)))
          else
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = 0.E0
          endif

        end do; end do
      end do; end do; end do; end do

    case('gauss')

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx 

          i_x = gx(ix) + int(n_x_grid/2.)
          if (i_x .gt. n_x_grid ) i_x = i_x - n_x_grid

          ! No initialisation of the zonal flow 
          if (imod.eq.iyzero) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  
          else if(spectral_radius) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini *            & 
               & exp(-(real(imod)/real(2*nmod))**2)*         &
               & exp(-(real(i_x)/real(2*n_x_grid))**2)
          else
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini *            & 
               & exp(-(real(imod)/real(2*nmod))**2)*         &
               & exp(-(real(i_x)**2/(2.*dxgr)**2))
          endif

        end do; end do 
      end do; end do; end do; end do

    case('zero')

      do imod = 1, nmod; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        do is = 1, nsp; do ix = 1, nx 

          ! No initialisation of the dist function
          ! For use in conjunction with the tearing mode the
          ! presence of an island is a sufficient perturbation 
          ! Note the g2f correction is applied below !!
          fdis(indx(ifdis,imod,ix,i,j,k,is)) = (0.E0,0.E0)  

        end do; end do
      end do; end do; end do; end do

    case('zonal')
      ! For the Rosenbluth-Hinton-test, all non-zonal components are zero and
      ! the zonal component is used to produce a radial electric field
 
      fdis(:) = (0.E0,0.E0)

      do imod = 1, nmod; do i = 1, ns;  do j = 1, nmu 
        do k = 1, nvpar; do is = 1, nsp 
 
          if (iyzero == imod .and. signz(is) > 0.0) then
            if (spectral_radius) then
              
              ! find central mode
              idum = 0
              do ix = 1, nx
                if (ixzero == ix) idum = ix
              enddo
              
              if (idum.eq.0) call gkw_abort('No kx = 0 mode found in init zonal')
              if ( (idum.le.1) .or. (idum.ge.nx) )                            &
                &            call gkw_abort('No finite kx mode found in init')
              
              fdis(indx(ifdis,imod,ixzero-1,i,j,k,is)) = -ci1 * amp_ini *   &
                & fmaxwl(ixzero-1,i,j,k,is)/2.
              fdis(indx(ifdis,imod,ixzero+1,i,j,k,is)) = ci1 *amp_ini *     &
                & fmaxwl(ixzero+1,i,j,k,is)/2.
              
            else ! (.not. spectral_radius)
              
              do ix = 1, nx
                if (types(is) .eq. 'EP') then
                  fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini *         &
                    & f_EP(ix,i,j,k,is) * (sin(2*pi*gx(ix)/n_x_grid))
                else
                  fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini *         &
                        & fmaxwl(ix,i,j,k,is) * (sin(2*pi*gx(ix)/n_x_grid))
                endif
              enddo
              
            endif
          endif

        end do; end do
      end do; end do; end do

    case('sgauss') !For single mode runs only

      do is = 1, nsp; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        vperp = sqrt(2.E0*bn(1,i)*mugr(j))
        fdis(indx(ifdis,1,1,i,j,k,is)) = amp_ini*exp(-((vperp-2.0)/0.5)**2) &
           &                        * exp(-((vpgr(i,j,k,is)+2.E0)/0.5)**2)
      end do; end do; end do; end do

    case('sgauss2')

      ! positive part
      do is = 1, nsp; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        vperp = sqrt(2.E0*bn(1,i)*mugr(j))
        fdis(indx(ifdis,1,1,i,j,k,is)) = amp_ini*                     &
           & exp(-((vpgr(i,j,k,is)+2.E0)/0.5)**2)*                     &
           & exp(-((vperp-2.0)/0.5)**2)
      end do; end do; end do; end do

      ! negative part
      do is = 1, nsp; do i = 1, ns; do j = 1, nmu; do k = 1, nvpar
        vperp = sqrt(2.E0*bn(1,i)*mugr(j))
        fdis(indx(ifdis,1,1,i,j,k,is)) = fdis(indx(ifdis,1,1,i,j,k,is)) &
           & -amp_ini * exp(-((vpgr(i,j,k,is)-2.E0)/0.5)**2)*            &
           & exp(-((vperp-2.0)/0.5)**2)
      end do; end do; end do; end do

    case('helical')
      !A helical density perturbation resonant with the magnetic island
      !to study GAM dynamnics in the presence of a magnetic island.
      !ikxspace_local = 1
      if (mode_box) then
        ikxspace_local = ikxspace
      else !if not mode_box ikxspace has no meaning in the rest of the code
        ikxspace_local = 1
        !ixzero set in mode.
      end if

      if (nmod.eq.1) then
        imodisland=1
        en = 1
      else if (nmod.gt.1) then       
        imodisland=isl_mode
        en = imodisland-1
      end if
      if(imodisland.gt.nmod) then
        call gkw_abort('Tearing mode: The poloidal mode the island is &
           &  defined as must be less than nmod')
      end if

      if(root_processor) then
        write(*,*)'Helical density perturbation chosen to be resonant &
           & with magnetic island'
      end if

      do is = 1,nsp; do j = 1,nmu; do k = 1,nvpar; do i=1,ns
        balloonpos = ikxspace_local*en*sgr(i)
        if (abs(sgr(i)).lt.1e-10) then
          dumbuf = pi*fmaxwl(ixzero,i,j,k,is)*exp(-(balloonpos/isl_ls)**2)
          fdis(indx(ifdis,imodisland,ixzero,i,j,k,is)) = dumbuf*ci1
        else
          dumbuf = sin(pi*balloonpos)*exp(-(balloonpos/isl_ls)**2)     &
             & *fmaxwl(ixzero,i,j,k,is)/balloonpos
          fdis(indx(ifdis,imodisland,ixzero,i,j,k,is)) = dumbuf*ci1
        end if

        !number of non-zero positive frequency radial modes
        modnum = (nx-1)/2

        if (nx.gt.1) then
          do p = 1,modnum
            iplusmode = ixzero+p
            iminusmode = ixzero-p
            balloonpos = (ikxspace_local*en*sgr(i)+1.E0*p)
            dumbuf = fmaxwl(iminusmode,i,j,k,is)*sin(pi*balloonpos)*   &
               &  exp(-(balloonpos/isl_ls)**2)/balloonpos
            fdis(indx(ifdis,imodisland,iminusmode,i,j,k,is)) = dumbuf*ci1
            balloonpos = (ikxspace_local*en*sgr(i)-1.E0*p)
            dumbuf = fmaxwl(iplusmode,i,j,k,is)*sin(pi*balloonpos)*    &
               &  exp(-(balloonpos/isl_ls)**2)/balloonpos
            fdis(indx(ifdis,imodisland,iplusmode,i,j,k,is)) = dumbuf*ci1
          end do
        end if

      end do; end do; end do; end do

    case('kxzero') !A single mode with kx=0, useful for testing

      if (ixzero.eq.0) call gkw_abort("No kx = 0 mode")
      fdis=(0.E0,0.E0)
      do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp 
        fdis(indx(ifdis,3,ixzero,i,j,k,is))=amp_ini
      end do; end do; end do; end do

    case('kyzero') !A single mode with ky=0, useful for testing

      fdis=(0.E0,0.E0)
      do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp 
        fdis(indx(ifdis,1,ixzero+3,i,j,k,is))=amp_ini
      end do; end do; end do; end do

    case('line') !A line of gaussians (in real space)

      fdis=(0.E0,0.E0)
      do i = 1, ns; do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp 
        do ix = 1, nx; do imod = 1, nmod
          i_x = gx(ix)

          a = 0.0
          if (spectral_radius) a = kxrh(ix)

          do line=-49,50 !100 Gaussians shifted by 0.01
            idx=indx(ifdis,imod,ix,i,j,k,is)
            fdis(idx)=fdis(idx) + amp_ini * (                       &
               & exp(-(real(i_x-ixzero)/real(2*n_x_grid))**2) +       &
               &              exp(-(real(imod)/(2*nmod)**2)) )        &
               * exp(ci1*a*line*0.01*lx) 
          end do
        end do; end do
      end do; end do; end do; end do

    case('ecurrent') !Initialisation of the electron distribution function
      !so that there is a radial current profile so as to induce a tearing mode.
      !Velocity space is a shifted Maxwellian.
      fdis=(0.E0,0.E0)
      a=0.5  
      do imod = 1,nmod; do ix = 1,nx; do i = 1, ns; do j = 1, nmu
        do k = 1, nvpar; do is = 1, nsp 
          !Parallel current in electrons to trigger a tearing mode.  
          !Perturbation is uniform toroidally and uniform along s,
          !with a gaussian shape radially.
          if ((signz(is).lt.0).and.(imod.eq.1)) then
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = -amp_ini             &
               &  *vthrat(is) *vpgr(i,j,k,is) *fmaxwl(ix,i,j,k,is)   &
               &  *sqrt(1.E0*pi/a) *exp(-pi*pi*kxrh(ix)**2)
          else
            fdis(indx(ifdis,imod,ix,i,j,k,is)) = amp_ini              &
               &      *de(ix,is)*(cos(2*pi*sgr(i))+1.0) 
          end if
        end do; end do
      end do; end do; end do; end do

    case default
      call gkw_abort('Unknown initilization switch finit: '//trim(finit))

    end select
  end if

  ! Multiply the initial distribution with some factor, that will make it zero
  ! near the axis and the plasma edge.
  if (lfinit_radial_dirichlet) then
    do imod = 1,nmod; do ix = 1,nx; do i = 1, ns; do j = 1, nmu
      do k = 1, nvpar; do is = 1, nsp
        i_x = gx(ix)
        ! \note The product (\tanh(+x) + 1)(\tanh(-x) + 1) \leq 1, thus there
        !       is no factor 1/2 needed in front.
        fdis(indx(ifdis,imod,ix,i,j,k,is)) = fdis(indx(ifdis,imod,ix,i,j,k,is)) &
             &      *(tanh(0.1/tan(pi*i_x/n_x_grid)) + 1.E0)*(tanh(-0.1/tan(pi*i_x/n_x_grid)) + 1.E0)
      end do; end do
    end do; end do; end do; end do
  end if

  ! Calculate the fields
  if (spectral_radius) then
    call calculate_fields(fdis)
  else
    call field_solve_nonspec_wrap(fdis,0.,.false.)
  endif


  ! In electromagnetic cases where the initial apar field is not zero
  ! (e.g. tearingmodes)
  ! correct the initialised distribution initialised by init_fdis
  ! to correspond to f, not g

  if (tearingmode.and..not.restarted) then
    do imod = 1, nmod; do ix = 1, nx; do i = 1, ns; do is = 1, nsp
      do j = 1, nmu; do k = 1, nvpar

        apar_ga = get_averaged_apar(imod,ix,i,j,is,fdis)
        
        fdis(indx(ifdis,imod,ix,i,j,k,is)) =           &
            & fdis(indx(ifdis,imod,ix,i,j,k,is))+      &
            & 2.E0*signz(is)*vthrat(is)*vpgr(i,j,k,is)  &
            & *fmaxwl(ix,i,j,k,is)/tmp(ix,is)*apar_ga

      end do; end do
    end do; end do; end do; end do

  end if
  
  ! Add the Maxwell in the case of a run without a background 
  ! distribution 
  if (fullf_wo_Fm) then 
  
    imod = iyzero 
    if (imod == 0) call gkw_abort('No zero mode to store FM (fullf_wo_Fm = .true.)')
    
    do ix = 1, nx; do i = 1, ns; 
    
      ! The total charge is subtracted 
      !WARNING needs to be parallelized 
      dum = 0.E0 
      do is = 1, nsp; do j = 1, nmu; do k = 1, nvpar 
        dum = dum + signz(is)*bn(ix,i)*intvp(i,j,k,is)*intmu(j)*fmaxwl(ix,i,j,k,is) 
      end do; end do; end do; 

      do is = 1, nsp; do j = 1, nmu; do k= 1, nvpar         
        fdis(indx(ifdis,imod,ix,i,j,k,is)) = fdis(indx(ifdis,imod,ix,i,j,k,is)) + & 
             (fmaxwl(ix,i,j,k,is) - dum) / rhostar              
      end do; end do; end do 
      
    end do; end do 
    
  endif 

end subroutine init_fdis

!-----------------------------------------------------------------------------
!> Initialise the array containing the Maxwellian
!> (or slowing down distribution)
!-----------------------------------------------------------------------------
subroutine init_fmaxwl

  use grid,         only : nx, ns, nsp, nmu, nvpar
  use dist,         only : fmaxwl, falpha
  use dist,         only : f_EP, df_EPdv, df_EPdv_sinhc
  use dist,         only : ghost_points_mu
  use geom,         only : bn
  use components,   only : de, types, pbg 
  use components,   only : dgrid, tmp, tgrid
  use components,   only : energetic_particles
  use components,   only : vpar_mean
  use rotation,     only : cfen
  use constants,    only : pi
  use global,       only : r_tiny
  use velocitygrid, only : vpgr, mugr!, intmu, intvp
  use functions,    only : mod_sinh_gkw
  use control,      only : order_of_the_scheme  

  ! integers for the loop over all grid points 
  integer :: ix ,i, j, k, is
  integer :: is_extra, ivpar_extra 

  ! Dummy variables 
  logical :: alpha_dist
  logical :: energetic_dist
  real    :: vn, par!, norm

  !A bit of a hack so that the extended (ghost points)
  !Maxwellian is initialised even when not parallelised.
  is_extra=2
  ivpar_extra=2
  if(order_of_the_scheme.eq.'second_order')then
    is_extra=1
    ivpar_extra=1
  endif

  ! calculate the normalized maxwellian w**3 / n_grid Fm
  do i = 1-is_extra, ns+is_extra
    do j = 1-ghost_points_mu, nmu+ghost_points_mu
      do k = 1-ivpar_extra, nvpar+ivpar_extra
        do ix = 1, nx; do is = 1, nsp
          if (types(is)=='EP') then
            fmaxwl(ix,i,j,k,is) = 0.E0
          else
            fmaxwl(ix,i,j,k,is) = exp(-(vpgr(i,j,k,is)**2+             &
                     &  2.E0*bn(ix,i)*mugr(j))/(tmp(ix,is)/tgrid(is))) &
                     &  / (sqrt(tmp(ix,is)*pi/tgrid(is))**3)  
          end if
          if (dgrid(is) > r_tiny) then
            fmaxwl(ix,i,j,k,is) = de(ix,is)*fmaxwl(ix,i,j,k,is) / dgrid(is)
          end if
          !Centrifugal parallel density variation,  cfen=0 if vcor=0
          fmaxwl(ix,i,j,k,is)=fmaxwl(ix,i,j,k,is)*exp(-cfen(i,is))
        end do; end do 
      end do
    end do
  end do

  ! the non Maxwellian distributions
 
  alpha_dist     = .false.
  energetic_dist = .false.
  par            = -1.
  do is = 1, nsp 
    if (types(is).eq.'alpha') then 
      alpha_dist = .true.
      par        = pbg(1,is)
    endif
    if ((types(is).eq.'EP').and.(energetic_particles)) then
      energetic_dist = .true.
    end if
  end do
  !Slowing-down distribution for alpha particles
  if (alpha_dist) then 
    ! WARNING this part works only for the vpar_max = 3. 
    ! specify as temperature E_alpha / (9 T_ref) 
    ! the parameter is pbg = 9 E_c / E_alpha 
    do i = 1, ns; do j = 1, nmu ; do k = 1, nvpar 
      vn = sqrt(vpgr(i,j,k,is)**2+2.E0*bn(ix,i)*mugr(j)) 
      if (vn.lt.3) then 
        falpha(i,j,k) = 3.E0 / (4.*pi*log(1.E0 + 27.E0*par**(-1.5))* & 
                      & (par**1.5 + vn**3))
      endif
    end do; end do; end do
  end if
  !Double shifted Maxwellian for NBI-like distribution
  if (energetic_dist) then
    do i = 1, ns; do j = 1,nmu; do k = 1, nvpar; do ix = 1, nx; do is = 1, nsp
      if (types(is) == 'EP') then
        f_EP(ix,i,j,k,is) = exp(-2.E0*bn(ix,i)*mugr(j)/                      &
            & (tmp(ix,is)/tgrid(is)))*                                       &
            & (                                                              &
            & exp(-((vpgr(i,j,k,is)-vpar_mean)**2)/(tmp(ix,is)/tgrid(is)))   &
            & +                                                              &
            & exp(-((vpgr(i,j,k,is)+vpar_mean)**2)/(tmp(ix,is)/tgrid(is)))   &
            & )*                                                             &
            & 0.5E0/(sqrt(tmp(ix,is)*pi/tgrid(is))**3)
        df_EPdv(ix,i,j,k,is) = - 1.E0/(tmp(ix,is)/tgrid(is)) *                         &
            & exp(-2.E0*bn(ix,i)*mugr(j)/                                    &
            & (tmp(ix,is)/tgrid(is)))*                                       &
            & ( (vpgr(i,j,k,is)-vpar_mean) *                                 &
            & exp(-((vpgr(i,j,k,is)-vpar_mean)**2)/(tmp(ix,is)/tgrid(is)))   &
            & + (vpgr(i,j,k,is)+vpar_mean) *                                 &
            & exp(-((vpgr(i,j,k,is)+vpar_mean)**2)/(tmp(ix,is)/tgrid(is)))   &
            & )*                                                             &
            & 1.E0/(sqrt(tmp(ix,is)*pi/tgrid(is))**3)
        df_EPdv_sinhc(ix,i,j,k,is) = exp(-(vpgr(i,j,k,is)**2 +               &
            &  2.E0*bn(ix,i)*mugr(j))/(tmp(ix,is)/tgrid(is)))                &
            &  / (sqrt(tmp(ix,is)*pi/tgrid(is))**3) *                        &
            & exp(-vpar_mean**2/(tmp(ix,is)/tgrid(is))) *                    &
            & 4.E0*vpar_mean/(tmp(ix,is)/tgrid(is))*                         &
            & mod_sinh_gkw(2.E0*vpar_mean*vpgr(i,j,k,is)/(tmp(ix,is)))*      &
            & 0.5E0/(sqrt(tmp(ix,is)*pi/tgrid(is))**3)
      else
        f_EP(ix,i,j,k,is)          = 0.E0
        df_EPdv(ix,i,j,k,is)       = 0.E0
        df_EPdv_sinhc(ix,i,j,k,is) = 0.E0
      end if
      if (dgrid(is) > r_tiny) then 
        f_EP(ix,i,j,k,is) = de(ix,is)*f_EP(ix,i,j,k,is) / dgrid(is) 
        df_EPdv(ix,i,j,k,is) = de(ix,is)*df_EPdv(ix,i,j,k,is) / dgrid(is) 
        df_EPdv_sinhc(ix,i,j,k,is) = de(ix,is)*df_EPdv_sinhc(ix,i,j,k,is)    &
                                   & / dgrid(is) 
      end if
    end do; end do; end do; end do; end do
  end if
  ! Renormalize the Maxwell to ensure the correct density.
  ! This is switched off at the moment since this is incomplete
  ! when parallizing over the magnetic moment.
  ! do i = 1, ns; do is =1,nsp ; do ix = 1, nx  
  !    norm = 0.E0
  !    do j = 1, nmu ; do k = 1, nvpar 
  !      norm = norm + bn(ix,i)*intmu(j)*intvp(i,j,k,is)*fmaxwl(ix,i,j,k,is)
  !    end do; end do
  ! end do; end do; end do 


end subroutine init_fmaxwl
  

!--------------------------------------------------------------------
!> This routine calls the routines that allocate the most fundamental arrays in
!> several modules.
!--------------------------------------------------------------------
subroutine allocate_core

  use matdat,           only : matdat_allocate 
  use components,       only : components_allocate 
  use geom,             only : geom_allocate 
  use mode,             only : mode_allocate
  use rotation,         only : rotation_allocate
  use fields,           only : fields_allocate

  ! Allocate the arrays of matdat.f90 
  call matdat_allocate
    
  ! Allocate the arrays of fields.F90 
  call fields_allocate

  ! Allocate the arrays of components.f90 
  call components_allocate

  ! Allocate the arrays of geom.f90 
  call geom_allocate

  ! Allocate the arrays of mode.f90  
  call mode_allocate  

  ! Allocate the arrays in rotation.f90
  call rotation_allocate
 
end subroutine allocate_core

!--------------------------------------------------------------------
!> This routine calls the routines that allocate other arrays.
!> (These arrays may have sizes which had to be computed before. And
!> may that computation may need an array allocated in the
!> routine allocate_core(). This is why a second round of allocation can
!> make sense.)
!--------------------------------------------------------------------
subroutine allocate_more
  use non_linear_terms, only : nonlinear_allocate
  use source_time, only : source_time_allocate

  ! Allocate the arrays of non_linear_terms.f90 Note 
  ! that these arrays are also used for the recon-
  ! struction of linear modes if mode_box is true 
  call nonlinear_allocate

  call source_time_allocate

end subroutine allocate_more

!-----------------------------------------------------------------------------
!> This routine calls the various routines that deallocate the arrays. 
!> Many of the arrays are, at present, not deallocated at the end of the 
!> run. Of course, this is not strictly necessary.
!-----------------------------------------------------------------------------
subroutine deallocate_runtime_arrays

  use exp_integration, only : exp_integration_deallocate
  use fields,          only : fields_deallocate
  use mpighosts, only : persistent_comm_end

  call exp_integration_deallocate

  call persistent_comm_end()
  call fields_deallocate

end subroutine deallocate_runtime_arrays

subroutine finalize()
  use components,      only : components_deallocate
  use diagnostic, only : diagnostic_finalize
  use io, only    : finalize_io
  use matrix_format, only : matrix_format_finalize

  call components_deallocate

  ! clean up any diagnostics and close files
  call diagnostic_finalize()
  
  ! close the IO interfaces
  call finalize_io()

  ! close any sparse matrix libraries
  call matrix_format_finalize()

end subroutine finalize

end module init
