!------------------------------------------------------------------------------
!> This module provides switches and some utility routines which are
!> seen and used by several of the specific diagnos_* modules.
!------------------------------------------------------------------------------
module diagnos_generic
  use global, only : BPAR_GA_FIELD

  implicit none

  private
 
  public :: set_default_nml_values
  public :: init, finalize
  public :: bcast, check, allocate_mem
  public :: allocate_spec_cmpx_buffers
  public :: compute_fft_sizes
  public :: four2real_2D
  public :: dfieldds, dfielddx, d2fielddx2
  public :: dphigadx_xy_point_get, dphigadzeta_xy_point_get
  public :: dphigads_xy_point_get, phiga_xy_point_get
  public :: dapargadx_xy_point_get, dapargadzeta_xy_point_get
  public :: dapargads_xy_point_get

  public :: velocity_slice_output
  public :: kyx_output_array, kyxsp_output_array
  public :: xy_output_array, xysp_output_array
  public :: xs_kyzero_output_array

  public :: inv_fft_xy_slice
  public :: get_tr_ps_mask
  
  public attach_metadata_grid

  interface attach_metadata_grid
    module procedure attach_metadata_grid_1d
    module procedure attach_metadata_grid_2d
    module procedure attach_metadata_grid_3d
    module procedure attach_metadata_grid_4d
    module procedure attach_metadata_grid_5d
    module procedure attach_metadata_grid_6d
    module procedure attach_metadata_grid_7d
    module procedure attach_metadata_grid_1d_name
    module procedure attach_metadata_grid_2d_name
    module procedure attach_metadata_grid_3d_name
    module procedure attach_metadata_grid_4d_name
  end interface attach_metadata_grid


  !> output switches which are seen by and used in several of the specific
  !> diagnos_* modules:
  
  !> Meta switch for all per timestep output (fluxes, time traces, 2D)
  logical, save, public :: lwrite_output1
  !> logical to output XY slices *every* time
  !> step, otherwise just output at the end of the run:
  logical, save, public :: xy_estep
  
  !> switch that sets whether radial fluxes and background profiles
  !> are written out
  logical, save, public :: lradial_profile
  
  !> Logical if you want to follow the magnetic island rotation around
  logical, save, public :: lisl_follow

  !> 2d diagnostics will output slices at the global s location xy_slice_ipar_
  integer, save, public :: xy_slice_ipar

  !> Logicals to select various XY diagnostic outputs
  logical, save, public :: xy_fluxes
  logical, save, public :: xy_fluxes_fsa
  logical, save, public :: xy_fluxes_k
  logical, save, public :: xy_fluxes_bi
  logical, save, public :: xy_fluxes_p
  logical, save, public :: xy_fluxes_v
  logical, save, public :: xy_fluxes_em
  logical, save, public :: xy_fluxes_bpar
  logical, save, public :: xy_dens
  logical, save, public :: xy_temp
  logical, save, public :: xy_current
  logical, save, public :: xy_current2
  !> If you want the xy data also outputted per mode rather than in 
  !> real space.  
  logical, save, public :: xy_spec

  !> logicals for 3D output
  logical, save, public :: apa3d
  logical, save, public :: spc3d
  logical, save, public :: phi3d
  logical, save, public :: apc3d
  logical, save, public :: den3d
  logical, save, public :: ene3d
  logical, save, public :: bpa3d
  logical, save, public :: bpc3d
  !> any_3d_output is the logical OR of the *3d switches above.
  logical, save, public :: any_3d_output
  
  !> Interval of large time steps for output of 3d data.
  integer, save, public :: out3d_interval
  

  !----------------------------------------------------------------
  ! Deprecated, obsolete diagnostic namelist entries:
  !----------------------------------------------------------------
  logical, save, public :: lphi_diagnostics
  ! This switch does not switch
  ! "phi diagnostics" as the name suggests, but acts as a meta
  ! switch for 2d diagnostics.
  !----------------------------------------------------------------

  !> Suppression factor for zonal flows in the 3d output
  real, save, public :: zonal_scale_3d

  !> complex slice in xy for mpi reduction
  complex, save, allocatable, dimension(:,:), public :: c_xy_local
  complex, save, allocatable, public :: out3dbuf_local_spec_cmpx(:,:,:)
  complex, save, allocatable, public :: out3dbuf_global_spec_cmpx(:,:,:)
  !> real slice for outputing xy
  real, save, allocatable, dimension(:,:), public :: r_xy_local
  real, save, allocatable, dimension(:,:,:), public :: r_xysp_local


  !> a mpi datatype for a perpendicular slice
  integer, save, public :: mpi_dtype_yx_slice
  integer, save, public :: mpi_dtype_yxsp_slice
  integer, save, public :: mpi_dtype_kyx_slice
  integer, save, public :: mpi_dtype_kyxsp_slice

  !> mpi datatypes for 3d real-valued subarrays, for realspace and kspace
  integer, save, public :: mpi_dtype_real_yxs
  integer, save, public :: mpi_dtype_real_spec_yxs

  !> mpi datatypes for 3d complex-valued subarrays, for realspace and kspace
  integer, save, public :: mpi_dtype_cpx_yxs
  integer, save, public :: mpi_dtype_cpx_spec_yxs

  !> mpi datatype for 2d real-valued velocity space
  integer, save, public :: mpi_dtype_real_vparmu

  !> buffer of the size of the global velocity space grid
  real, allocatable, dimension(:,:), save :: global_vspace

  !> a buffer for xy_output_array
  real, save, allocatable, public :: r_xy_global(:,:)
  real, save, allocatable, public :: r_xysp_global(:,:,:)
  
  
  !> real slice for outputting sx
  real, save, allocatable, dimension(:,:), public :: r_xs_local
  
  !> a buffer for xs_output_array
  real, save, allocatable, public :: r_xs_global(:,:)
  complex, save, allocatable, dimension(:,:), public :: c_xs_local

  
  !> mpi datatype for 2d xs-slice
  integer, save, public :: mpi_dtype_xs_slice
  
  
  !----------------------------------------------------------------
  ! other useful stuff for diagnostics
  !----------------------------------------------------------------
  !> The Parseval correction factor, documented in manual Sec. 10.3  
  !> = 2 for non-zero modes, = 1 for the zero mode
  !> (Hermitian symmetry of binormal Fourier representation)
  integer, save, allocatable, dimension(:), public :: parseval_correction
  integer, save, public :: mphit, mphiw3t, mrad_G, mrad_l

  !> a string variable, which holds either the name of the spectral
  !> radial grid (in case of spectral runs) or the name of the
  !> position space grid (in case of nonspectral runs). This is useful
  !> when it comes to setting grid metadata for spectral fields.
  character (len=10), save, public :: rad_gridname

  integer, parameter, public :: PARTICLE_FLUX = 1, ENERGY_FLUX = 2
  integer, parameter, public :: MOMENTUM_FLUX = 3

  ! diagnostics that want to iterate over both fields and moments may need this:
  integer, parameter, public :: START_MOMENTS = BPAR_GA_FIELD
  
  integer, parameter, public :: DENSITY_MOMENT = 1
  integer, parameter, public :: TEMPERATURE_MOMENT = 2
  integer, parameter, public :: CURRENT_MOMENT = 3
  integer, parameter, public :: CURRENT_SQ_MOMENT = 4
  integer, parameter, public :: VPERP_SQ_MOMENT = 5
  integer, parameter, public :: PAR_TEMPERATURE_MOMENT = 6
  integer, parameter, public :: PERP_TEMPERATURE_MOMENT = 7
  integer, parameter, public :: QPAR_MOMENT = 8
  integer, parameter, public :: M12_MOMENT = 9
  integer, parameter, public :: M24_MOMENT = 10
  integer, parameter, public :: PASSING_DENSITY_MOMENT = 11
  integer, parameter, public :: TRAPPED_DENSITY_MOMENT = 12
  integer, parameter, public :: DENSITY_GA_MOMENT = 13
  integer, parameter, public :: DENSITY_POLAR_MOMENT = 14
  integer, parameter, public :: CURRENT_GA_MOMENT = 15
  integer, parameter, public :: PAR_TEMPERATURE_GA_MOMENT = 16
  integer, parameter, public :: PERP_TEMPERATURE_GA_MOMENT = 17
  integer, parameter, public :: QPAR_GA_MOMENT = 18
  integer, parameter, public :: M12_GA_MOMENT = 19
  integer, parameter, public :: M24_GA_MOMENT = 20
  integer, parameter, public :: PERP_TEMPERATURE_J1_MOMENT = 21
  ! gyro-averaged phi moment of background and perturbed distribution
  integer, parameter, public :: PHI_GA_FM_MOMENT = 22
  integer, parameter, public :: PHI_GA_DELTAF_MOMENT = 23
 
  integer, parameter, public :: N_MOMENTS = 23


  integer, parameter, public :: max_token_length = 10
  character(len=max_token_length), parameter, public :: PHI_FIELD_TOKEN = 'phi'
  character(len=max_token_length), parameter, public :: APAR_FIELD_TOKEN = 'Apar'
  character(len=max_token_length), parameter, public :: BPAR_FIELD_TOKEN = 'Bpar'
  character(len=max_token_length), parameter, public :: DISTRIBUTION_TOKEN = 'fdis'
  character(len=max_token_length), parameter, public :: PHI_GA_FIELD_TOKEN = 'phi_ga'
  character(len=max_token_length), parameter, public :: APAR_GA_FIELD_TOKEN = 'Apar_ga'
  character(len=max_token_length), parameter, public :: BPAR_GA_FIELD_TOKEN = 'Bpar_ga'
  character(len=max_token_length), parameter, public :: DENSITY_MOMENT_TOKEN = 'dens'
  character(len=max_token_length), parameter, public :: TEMPERATURE_MOMENT_TOKEN = 'T'
  character(len=max_token_length), parameter, public :: CURRENT_MOMENT_TOKEN = 'vpar'
  character(len=max_token_length), parameter, public :: CURRENT_SQ_MOMENT_TOKEN = 'vparsq'
  character(len=max_token_length), parameter, public :: VPERP_SQ_MOMENT_TOKEN = 'vperpsq'
  character(len=max_token_length), parameter, public :: PAR_TEMPERATURE_MOMENT_TOKEN = 'Tpar'
  character(len=max_token_length), parameter, public :: PERP_TEMPERATURE_MOMENT_TOKEN = 'Tperp'
  character(len=max_token_length), parameter, public :: QPAR_MOMENT_TOKEN = 'Qpar'
  character(len=max_token_length), parameter, public :: M12_MOMENT_TOKEN = 'M12'
  character(len=max_token_length), parameter, public :: M24_MOMENT_TOKEN = 'M24'
  character(len=max_token_length), parameter, public :: PASSING_DENSITY_MOMENT_TOKEN = 'dens_ps'
  character(len=max_token_length), parameter, public :: TRAPPED_DENSITY_MOMENT_TOKEN = 'dens_tr'
  character(len=max_token_length), parameter, public :: DENSITY_GA_MOMENT_TOKEN = 'dens_ga'
  character(len=max_token_length), parameter, public :: DENSITY_POLAR_MOMENT_TOKEN = 'dens_polar'
  character(len=max_token_length), parameter, public :: CURRENT_GA_MOMENT_TOKEN = 'vpar_ga'
  character(len=max_token_length), parameter, public :: PAR_TEMP_GA_MOMENT_TOKEN = 'Tpar_ga'
  character(len=max_token_length), parameter, public :: PERP_TEMP_GA_MOMENT_TOKEN = 'Tperp_ga'
  character(len=max_token_length), parameter, public :: QPAR_GA_MOMENT_TOKEN = 'Qpar_ga'
  character(len=max_token_length), parameter, public :: M12_GA_MOMENT_TOKEN = 'M12_ga'
  character(len=max_token_length), parameter, public :: M24_GA_MOMENT_TOKEN = 'M24_ga'
  character(len=max_token_length), parameter, public :: PERP_TEMP_J1_MOMENT_TOKEN = 'Tperp_J1'

  !> a list, to obtain the token from the iquantity number
  character(len=max_token_length), parameter, public :: &
     & tokens(START_MOMENTS+N_MOMENTS-2) = (/ &
     & PHI_FIELD_TOKEN, APAR_FIELD_TOKEN, BPAR_FIELD_TOKEN, &
     & DISTRIBUTION_TOKEN, &
     & PHI_GA_FIELD_TOKEN, APAR_GA_FIELD_TOKEN, BPAR_GA_FIELD_TOKEN, &
     & DENSITY_MOMENT_TOKEN, TEMPERATURE_MOMENT_TOKEN, &
     & CURRENT_MOMENT_TOKEN, CURRENT_SQ_MOMENT_TOKEN, &
     & VPERP_SQ_MOMENT_TOKEN, &
     & PAR_TEMPERATURE_MOMENT_TOKEN, PERP_TEMPERATURE_MOMENT_TOKEN, &
     & QPAR_MOMENT_TOKEN, M12_MOMENT_TOKEN, M24_MOMENT_TOKEN, &
     & PASSING_DENSITY_MOMENT_TOKEN, TRAPPED_DENSITY_MOMENT_TOKEN, &
     & DENSITY_GA_MOMENT_TOKEN, DENSITY_POLAR_MOMENT_TOKEN, &
     & CURRENT_GA_MOMENT_TOKEN, PAR_TEMP_GA_MOMENT_TOKEN, &
     & PERP_TEMP_GA_MOMENT_TOKEN, &
     & QPAR_GA_MOMENT_TOKEN, M12_GA_MOMENT_TOKEN, M24_GA_MOMENT_TOKEN, &
     & PERP_TEMP_J1_MOMENT_TOKEN &
     & /)

  !> constants needed to indicate what a particular diagnostic needs
  !> to be provided with.
  integer, parameter, public :: LOCAL_DATA = 1, S_GHOSTCELLS = 2
  integer, parameter, public :: X_GHOSTCELLS = 3, VPAR_GHOSTCELLS = 4
  integer, parameter, public :: MU_GHOSTCELLS = 5

contains

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init()
    use control, only : spectral_radius, io_legacy
    use grid, only : proc_subset, n_y_grid, n_s_grid
    use mpidatatypes, only : create_subarray_datatype
    use mpiinterface, only : MPIREAL_X, MPICOMPLEX_X
    use global, only : id_mod, id_x, id_sp, id_s, id_vpar, id_mu

    if(spectral_radius) then
      rad_gridname = 'kxrh'
    else
      rad_gridname = 'xgr'
    end if

    !if (proc_subset(0,0,0,1,1)) then ! would be better, but then
    !other procs must not try to resize the datatype in the gather
    !routine
    if(io_legacy) then
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_yx_slice, &
         & id_mod, id_x, &
         & global_and_local_ysize=mphit)
    else
      ! create subarray dataype for 2d real space field.
      ! the correct size is (as nx==n_x_grid for spectral runs)
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_yx_slice, &
         & id_mod, id_x, &
         & global_and_local_ysize=n_y_grid)
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_yxsp_slice, &
         & id_mod, id_x, id_sp, &
         & global_and_local_ysize=n_y_grid)
      !if (proc_subset(0,0,0,1,1)) then ! would be better, but then
      !certain procs must not try to resize this dtype inside the
      !gather routine.  FIXME maybe create dtype only on relevant s
      !coord procs
        ! create subarray dataype for 2d real space field:
        ! the correct size is (as nx==n_x_grid for spectral runs)
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_kyxsp_slice, &
         & id_mod, id_x, id_sp)
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_kyx_slice, &
         & id_mod, id_x)
    end if

    if(.not. io_legacy) then
      ! create subarray dataype for 3d real space field
      call create_subarray_datatype(MPICOMPLEX_X, mpi_dtype_cpx_yxs, &
         & id_mod, id_x, id_s, &
         & global_and_local_ysize=n_y_grid)
  
    end if
    ! create subarray MPI datatype for 3d spectral field
    call create_subarray_datatype(MPICOMPLEX_X, mpi_dtype_cpx_spec_yxs, &
       & id_mod, id_x, id_s)


    if (proc_subset(0,0,1,1,0)) then
      ! create subarray MPI datatype for 3d spectral field
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_real_spec_yxs, &
         & id_mod, id_x, id_s)

     if(io_legacy) then
        ! create subarray dataype for 3d real space field
        if (spectral_radius) then
          ! this is actually twice as large as it needs to be and
          ! wastes space
          call create_subarray_datatype(MPIREAL_X, mpi_dtype_real_yxs, &
             & id_mod, id_x, id_s,&
             & global_and_local_xsize=mrad_l, global_and_local_ysize=mphit)
        else
          call create_subarray_datatype(MPIREAL_X, mpi_dtype_real_yxs, &
             & id_mod, id_x, id_s,&
             & global_and_local_ysize=mphit)
        end if
      else
        ! instead, the correct size is (as nx==n_x_grid for spectral runs)
        call create_subarray_datatype(MPIREAL_X, mpi_dtype_real_yxs, &
           & id_mod, id_x, id_s,&
           & global_and_local_ysize=n_y_grid)
      end if
    endif

    call create_subarray_datatype(MPIREAL_X, mpi_dtype_real_vparmu, &
       & id_vpar, id_mu)
       
    call create_subarray_datatype(MPIREAL_X, mpi_dtype_xs_slice, &
       & id_x, id_s)

  end subroutine init



  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    use control, only : flux_tube
    use geom, only : isg_lfs
    lwrite_output1 = .true.

    lradial_profile=(.not. flux_tube)
    lphi_diagnostics = .true.

    xy_slice_ipar = isg_lfs
    xy_estep = .true.

    lisl_follow = .true.

    zonal_scale_3d = 1.0

    xy_fluxes = .false.
    xy_fluxes_fsa = .false.
    xy_fluxes_k = .false.
    xy_fluxes_bi = .false.
    xy_fluxes_p  = .false.
    xy_fluxes_v  = .false.
    xy_fluxes_em = .false.
    xy_fluxes_bpar = .false.
    xy_dens = .false.
    xy_temp = .false.
    xy_current = .false.
    xy_current2 = .false.

    xy_spec = .false.

    apa3d = .false.
    spc3d = .false.
    phi3d = .false.
    apc3d = .false.
    den3d = .false.
    ene3d = .false.
    bpa3d = .false.
    bpc3d = .false.
    
    out3d_interval = 1

  end subroutine set_default_nml_values


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lwrite_output1,1)

    call mpibcast(xy_slice_ipar,1)
    call mpibcast(lradial_profile, 1)

    call mpibcast(XY_fluxes,1)
    call mpibcast(XY_fluxes_fsa,1)
    call mpibcast(XY_fluxes_k,1)
    call mpibcast(XY_fluxes_bi,1)
    call mpibcast(XY_fluxes_V,1)
    call mpibcast(XY_fluxes_P,1)
    call mpibcast(XY_fluxes_em,1)
    call mpibcast(XY_fluxes_bpar,1)
    call mpibcast(XY_dens,1)
    call mpibcast(XY_temp,1)
    call mpibcast(XY_spec,1)
    call mpibcast(XY_current,1)
    call mpibcast(XY_current2,1)
    call mpibcast(apa3d,1)
    call mpibcast(bpa3d,1)
    call mpibcast(den3d,1)
    call mpibcast(ene3d,1)
    call mpibcast(spc3d,1)
    call mpibcast(apc3d,1)
    call mpibcast(bpc3d,1)
    call mpibcast(phi3d,1)
    call mpibcast(xy_estep,1)
    call mpibcast(out3d_interval,1)

    call mpibcast(lphi_diagnostics,1)

    call mpibcast(zonal_scale_3d,1)
    
    call mpibcast(lisl_follow,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine check()
    use mode, only : mode_box
    use grid, only : n_s_grid
    use control, only : spectral_radius, method
    use general, only : gkw_warn
    use geom, only : isg_lfs

    if (method == 'EIV') then
      call gkw_warn('lwrite_output1, the per-timestep diagnostics metaswitch,&
         & is used in a particular way for eigenmode diagnostics.&
         & Forced it to .true.')
      lwrite_output1 = .true.
    end if
    
    ! Select the global s-grid location at which we should make the xy slices.
    if (.not.(xy_slice_ipar >= 1 .and. xy_slice_ipar <= n_s_grid)) then
      ! default is the middle of the grid
      xy_slice_ipar = isg_lfs
      call gkw_warn('xy_slice_ipar out of range, will use central point')
    end if

    ! decide if there will be any 3D output, then print a warning
    any_3d_output = (apa3d .or. spc3d .or. phi3d .or. apc3d .or. den3d .or. &
       & ene3d .or. bpa3d .or. bpc3d)
    if (any_3d_output) then
      call gkw_warn('3D output data can be very large')

      if (.not. mode_box) then
        call gkw_warn('3D diagnostics: To get 3D output data, '// &
           & 'mode_box must also be set to true')
        apa3d = .false.
        spc3d = .false.
        phi3d = .false.
        apc3d = .false.
        den3d = .false.
        ene3d = .false.
        bpa3d = .false.
        bpc3d = .false.
        any_3d_output = .false.
      end if
    end if

    if (spectral_radius .and. lradial_profile) then
      call gkw_warn('The radial profiles have no physical meaning running spectral.')
    endif
    
    
    if(out3d_interval < 1) then
      call gkw_warn('out3d_interval less then one makes no sense. &
                   & Setting it to one!')
      out3d_interval = 1
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use grid, only : nmod, nx, nsp, n_x_grid, n_y_grid, number_of_species
    use grid, only : n_s_grid, ns
    use mode, only : iyzero
    use general, only : gkw_abort
    use control, only : io_legacy
    integer :: ierr

    call compute_fft_sizes(nmod, n_x_grid, nx, &
     & mphit, mphiw3t, mrad_G, mrad_l)
    
    allocate(parseval_correction(nmod), stat = ierr)
    if (ierr.ne.0) then
      call gkw_abort('Could not allocate parseval_correction diagnos_generic')
    endif

    ! might be allocated in vain - but if it does exist, it is a
    ! helpful buffer at several places in diagnostics.
    if(io_legacy) then
      allocate(r_xy_local(mphit,mrad_l),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xy_local')
    else
      allocate(r_xy_local(n_y_grid,nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xy_local')
      allocate(r_xysp_local(n_y_grid,nx,nsp),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xysp_local')
    end if

    parseval_correction = 2
    if(iyzero > 0) then
      ! then one mode is the zero mode
      parseval_correction(iyzero) = 1
    end if

    ! FIXME allocate this only if it is really needed
    if(io_legacy) then
      allocate(r_xy_global(mphit,mrad_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xy_global')
    else
      allocate(r_xy_global(n_y_grid,n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xy_global')
      allocate(r_xysp_global(n_y_grid, n_x_grid, number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xysp_global')
    end if
    
    allocate(r_xs_local(mrad_l,ns),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xs_local')
    
    allocate(r_xs_global(mrad_G,n_s_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_generic :: r_xs_global')
    
  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_spec_cmpx_buffers()
    use grid, only : nmod, nx, ns, n_x_grid, n_s_grid, proc_subset
    use mpiinterface, only : root_processor
    use general, only : gkw_abort
    integer :: ierr
    if(.not.allocated(out3dbuf_local_spec_cmpx)) then
      if(proc_subset(0,0,0,0,0))then
        allocate(out3dbuf_local_spec_cmpx(nmod,nx,ns),stat=ierr)
        if (ierr /= 0) call gkw_abort('diagnostic :: out3dbuf_local_spec_cmpx')
      else
        allocate(out3dbuf_local_spec_cmpx(1,1,1),stat=ierr)
      end if
    end if
    if (.not.allocated(out3dbuf_global_spec_cmpx)) then
      if(root_processor) then
        allocate(out3dbuf_global_spec_cmpx(nmod,n_x_grid,n_s_grid),&
           & stat=ierr)
        if (ierr /= 0) call gkw_abort('diagnostic :: out3dbuf_global_spec_cmpx')
      else
        allocate(out3dbuf_global_spec_cmpx(1,1,1),stat=ierr)
      end if
    end if
  end subroutine allocate_spec_cmpx_buffers
  

  !****************************************************************************
  !> Generic function to output velocity space data, handles MPI.
  !> Give it a local (real) 2D velocity space array and will gather and write the
  !> global one to a file.
  !> Call this routine with *all* processes.  isg and ispg are GLOBAL species 
  !> and s points for which you want to write. (if array is not dependent on
  !> these then choose 1,1 for example)
  !----------------------------------------------------------------------------
  subroutine velocity_slice_output(groupname,local_array,luname,ixg,isg,ispg, &
     & preference_if_mixed, tag)
    use io,           only : output_array
    use io,           only : xy_fmt
    use grid,         only : nmu,nvpar,n_vpar_grid,n_mu_grid
    use grid,         only : proc_subset, nx, ns, nsp
    use mpiinterface, only : gather_array, root_processor
    use mpicomms,     only : COMM_VPAR_NE_MU_NE

    character (len=*), intent(in) :: groupname
    real, dimension(nvpar,nmu), intent(in) :: local_array
    integer, intent(in) :: isg, ispg, ixg
    character (len=*), intent(in) :: luname
    character (len=*), intent(in) :: preference_if_mixed
    integer, intent(in) :: tag


    ! allocate array to contain the full slice
    if (.not. allocated(global_vspace)) then
      if(root_processor) then
        allocate(global_vspace(n_vpar_grid,n_mu_grid))
      else
        allocate(global_vspace(1,1))
      end if
    end if

    ! gather the data
    if (proc_subset(ixg,isg,0,0,ispg) .or. root_processor) then
      call gather_array(global_vspace,local_array,mpi_dtype_real_vparmu,&
         & COMM_VPAR_NE_MU_NE, .true., &
         & ixg <= nx .and. isg <= ns .and. ispg <= nsp, tag)
    end if

    ! write with 1 processor
    if (root_processor) then
      call output_array(trim(luname), groupname, &
         & global_vspace, 'F', xy_fmt, preference_if_mixed)
    end if

  end subroutine velocity_slice_output


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine compute_fft_sizes(nmod, n_x_grid, nx, &
     & mphit, mphiw3t, mrad_G, mrad_l)
    use control, only : spectral_radius
    use non_linear_terms, only : mphi, mphiw3, mrad
    use non_linear_terms, only : get_extended_firstdim_fft_size
    integer, intent(in) :: nmod, n_x_grid, nx
    integer, intent(out) :: mphit, mphiw3t, mrad_G, mrad_l
    
    !If the number of nmod is small the XY diagnostic lacks points in
    !the poloidal direction...this pads the ffts with
    !zeros so there are more points in that direction
    if(nmod < 4)then
      ! Calculate the size of the grid for the FFT
      call get_extended_firstdim_fft_size(4, mphit, mphiw3t)
    else
      mphit = mphi
      mphiw3t = mphiw3
    end if

    if (.not. spectral_radius) then
      mrad_G = n_x_grid
      mrad_l = nx ! /= mrad in nonlinear terms which contains ghost cells
    else
      mrad_l = mrad
      mrad_G = mrad
    end if

  end subroutine compute_fft_sizes

  !----------------------------------------------------------------------------
  !> This function uses c_xy_local as a helper buffer!
  !----------------------------------------------------------------------------
  subroutine inv_fft_xy_slice(c_kxxky_local, r_xy_local, zonal_scale)
    use components, only : tearingmode, isl_mode
    use tearingmodes, only : omega_rot
    use grid, only : nmod, nx, nymod
    use control, only : time, io_legacy
    use non_linear_terms, only : jind
    use general, only : gkw_abort
    use grid, only : jind_flexible
    use fft, only : FFT_INVERSE
    
    !> the input argument is a perpendicular slice in fourier representation,
    !> so either in kxky or in xky .
    complex, dimension(nmod,nx), intent(inout) :: c_kxxky_local
    !> output buffer size depends on io_legacy, because of backwards compatibility
    !> (mphit,mrad_l) or (n_y_grid,nx)
    real, dimension(:,:), intent(out) :: r_xy_local
    logical, optional, intent(in) :: zonal_scale

    integer :: ix, imod, ierr

    if(.not. allocated(c_xy_local)) then
      if(io_legacy) then
        allocate(c_xy_local(mphiw3t,mrad_l),stat=ierr)
      else
        allocate(c_xy_local(nymod,nx),stat=ierr)
      end if
      if (ierr /= 0) call gkw_abort('diagnos_generic :: c_xy_local')
    end if

    ! if run is radially spectral, reorder with respect to the radial
    ! dimension.
    c_xy_local = (0.,0.)
    do ix = 1, nx
      do imod = 1, nmod
        if(io_legacy) then
          c_xy_local(imod,jind(ix)) = c_kxxky_local(imod,ix)
        else
          if(nmod == 1) then
            c_xy_local(imod+1,jind_flexible(nx,nx,ix)) = c_kxxky_local(imod,ix)
          else
            c_xy_local(imod,jind_flexible(nx,nx,ix)) = c_kxxky_local(imod,ix)
          end if
        end if
      end do
    end do
    if (present(zonal_scale)) then
      if (zonal_scale) then
        ! Cheat on the zonal flow appearance (useful for illustrative pictures)
        c_xy_local(1,:) = zonal_scale_3d * c_xy_local(1,:)
      end if
    end if

    ! are we following an island?
    ! FIXME there is an 'imod' lying
    ! around here, whose value is equal to nmod. Is this looks as if
    ! it was not intended?
    if (tearingmode .and. lisl_follow) then
      c_xy_local = c_xy_local*exp(cmplx(0.,(imod-1)*omega_rot*time)/(isl_mode-1))
    end if

    ! call the inverse FFT routine
    call four2real_2D(r_xy_local,c_xy_local,FFT_INVERSE)
    
  end subroutine inv_fft_xy_slice
  
  !****************************************************************************
  !> This routine outputs a real field living on a global xy slice of position
  !> space.
  !> 
  !> It expects a complex array in local x-ky or kx-ky space.
  !> The Fourier Transform and global gather in x are performed if needed
  !> (depending on parallelisation scheme and spectral/nonspectral setting),
  !> and in y in any case.
  !> 
  !> This routine has to be called with all procs,
  !> and the proc_is_in_slice parameter must be set to true for the processes of
  !> the xy-slice, and false for all others.
  !----------------------------------------------------------------------------
  subroutine xy_output_array(array, zonal_scale, lun, &
     & proc_is_in_slice, root_is_in_slice)
    use io,               only : xy_fmt, binary_fmt
    use grid,             only : nx, nmod
    use mpiinterface,     only : gather_array, root_processor, mpicomm_rank
    use mpiinterface,     only : send_to_root, mpibarrier
    use mpicomms,         only : COMM_DUMMY, COMM_X_NE
!    use mpicomms,         only : COMM_X_EQ

    ! for output in a single file:
    use io, only : append_chunk
    use control, only : io_legacy
    
    logical, intent(in) :: proc_is_in_slice, root_is_in_slice
    logical, intent(in) :: zonal_scale
    !> the local xy-slice, value does only matter on processes with
    !> proc_is_in_slice == .true. and who are with root in COMM_VPAR_EQ_MU_EQ
    complex, dimension(nmod,nx), intent(inout) :: array
    integer, intent(in) :: lun

    integer :: rank
    
    if (proc_is_in_slice) then
      ! fourier transform the local slice
      call inv_fft_xy_slice(array, r_xy_local, zonal_scale)
    end if

    ! gather arrays with respect to x direction
    if(io_legacy) then
      if(proc_is_in_slice) then
        call gather_array(r_xy_global(1:mphit,1:mrad_G),mphit,mrad_G, &
         & r_xy_local(1:mphit,1:mrad_l),mphit,mrad_l, &
         & COMM_DUMMY,COMM_X_NE,ALLGATHER=.true.)
        ! now every process on the chosen perp slice has the global
        ! position space array.
      else
        ! the processes on other slices hold zeros.
        r_xy_global = 0
      end if
      ! if(proc_is_with_root_in(COMM_X_EQ)) then
      !   ! sum over the communicator which is complementary to COMM_X_NE,
      !   ! this effectively sends the array
      !   ! to the root process.
      !   call mpireduce_sum_inplace(r_xy_global, &
      !    & shape(r_xy_global), &
      !    & COMM_X_EQ)
      ! end if
      call mpicomm_rank(COMM_X_NE, rank)
      if((proc_is_in_slice .and. rank == 0 .and. .not. root_processor) &
         & .or. (root_processor .and. .not. (proc_is_in_slice .and. rank == 0))) then
        ! if the root process of COMM_X_NE is not the root of
        ! COMM_WORLD, then these two distinct processes will do this:
        call send_to_root(r_xy_global, shape(r_xy_global))
      end if
    else
      if(proc_is_in_slice .or. root_processor) then
        call gather_array(r_xy_global, &
           & r_xy_local, mpi_dtype_yx_slice, &
           & COMM_X_NE, .true., root_is_in_slice)
      end if
    end if

    ! as long as in general only the root process can do IO: send
    ! output data to root process and
    ! let root do the IO.

    if (root_processor) then

      ! NOTE with the rewrite of the io, both the binary
      ! and the ascii output produced by the following line are in 'F'
      ! order. Before, the ascii output used to be in 'C' order.
      call append_chunk(lun, r_xy_global, xy_fmt, binary_fmt)
    end if
    call mpibarrier()

  end subroutine xy_output_array
  
  
  
  !****************************************************************************
  !> This routine outputs the ky=0 component of a field living on a global 
  !> xs slice of position space.
  !----------------------------------------------------------------------------
  subroutine xs_kyzero_output_array(array, lun, &
     & proc_is_in_slice, root_is_in_slice)
    use io,               only : xy_fmt, binary_fmt, ascii_fmt
    use grid,             only : nx, ns, n_s_grid
    use mpiinterface,     only : gather_array, root_processor, mpicomm_rank
    use mpiinterface,     only : send_to_root, mpibarrier
    use mpicomms,         only : COMM_DUMMY, COMM_S_NE_X_NE, COMM_S_NE
    use mpicomms,         only : COMM_X_NE
    use control,          only : spectral_radius
    use fft,              only : FFT_INVERSE, fourcol
    use non_linear_terms, only : jind
    use general,          only : gkw_abort
    

    ! for output in a single file:
    use io,      only : append_chunk
    use control, only : io_legacy
    
    logical, intent(in) :: proc_is_in_slice, root_is_in_slice
    !> the local xs-slice, value does only matter on processes with
    !> proc_is_in_slice == .true. and who are with root 
    complex, dimension(nx,ns), intent(in) :: array
    integer, intent(in) :: lun
    
    integer :: is, ix, ierr
    integer :: rank1, rank2
    

    if(proc_is_in_slice) then
    
      ! since zonal mode, i.e., ky=0, and non-spectral the field is real
      if(.not. spectral_radius) then
        r_xs_local = real(array)
        
      ! if spectral radius make inverse fft in radial direction
      else
      
        if(.not. allocated(c_xs_local)) then
          if(io_legacy) then
            allocate(c_xs_local(mrad_l,ns),stat=ierr)
          else
            allocate(c_xs_local(nx,ns),stat=ierr)
          end if
          if (ierr /= 0) call gkw_abort('diagnos_generic :: could not allocate c_xs_local')
        end if
        
        
        ! reorder with respect to the radial dimension.
        c_xs_local = (0.,0.)
        do ix = 1, nx
          do is = 1, ns
            c_xs_local(jind(ix),is) = array(ix,is)
          end do
        end do
        
        ! make inverse fft
        call fourcol(c_xs_local, FFT_INVERSE)
        
        ! take the real part only, since iyzero is considered
        r_xs_local = real(c_xs_local)
            
      end if
      
    end if

    ! gather arrays with respect to x and s direction
    if(io_legacy) then
    
      if(proc_is_in_slice) then
        call gather_array(r_xs_global(1:mrad_G,1:n_s_grid),mrad_G,n_s_grid, &
         & r_xs_local(1:mrad_l,1:ns),mrad_l,ns, &
         & COMM_X_NE,COMM_S_NE,ALLGATHER=.true.)
        ! now every process on the chosen xs-slice has the global
        ! position space array.
      else
        ! the processes on other slices hold zeros.
        r_xs_global = 0
      end if
  
      call mpicomm_rank(COMM_S_NE, rank1)
      call mpicomm_rank(COMM_X_NE, rank2)
      
      if((proc_is_in_slice .and. rank1 == 0 .and. rank2 == 0 .and. .not. root_processor) &
         & .or. (root_processor .and. .not. (proc_is_in_slice .and. rank1 == 0 .and. rank2 == 0))) then
        ! if the root process of COMM_X_NE is not the root of
        ! COMM_WORLD, then these two distinct processes will do this:
        call send_to_root(r_xs_global, shape(r_xs_global))
      end if
    
    else
    
      if(proc_is_in_slice .or. root_processor) then
    
        ! gather local (in x and s) data in global array in x and s
        call gather_array(r_xs_global, r_xs_local, &
          & mpi_dtype_xs_slice, &
          & COMM_S_NE_X_NE, .true., root_is_in_slice)
        
      end if
      
    end if

    if (root_processor) then

      ! NOTE with the rewrite of the io, both the binary
      ! and the ascii output produced by the following line are in 'F'
      ! order. Before, the ascii output used to be in 'C' order.
      call append_chunk(lun, r_xs_global, xy_fmt, binary_fmt)

    end if
    
  end subroutine xs_kyzero_output_array
  
  

  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine xysp_output_array(array, lun, &
     & proc_is_in_slice, root_is_in_slice)
    use io,               only : xy_fmt, binary_fmt
    use grid,             only : nx,nsp, nmod
    use mpiinterface,     only : gather_array, root_processor
    use mpicomms,         only : COMM_SP_NE_X_NE
    
    ! for output in a single file:
    use io, only : append_chunk
    
    logical, intent(in) :: proc_is_in_slice, root_is_in_slice
    !> the local xy-slice, value does only matter on processes with
    !> proc_is_in_slice == .true. and who are with root in COMM_VPAR_EQ_MU_EQ
    complex, dimension(nmod,nx,nsp), intent(inout) :: array
    integer, intent(in) :: lun
    integer :: isp
    
    if (proc_is_in_slice) then
      ! fourier transform the local slice
      do isp = 1, nsp
        call inv_fft_xy_slice(array(:,:,isp), r_xysp_local(:,:,isp), .false.)
      end do
    end if
    if(proc_is_in_slice .or. root_processor) then
      call gather_array(r_xysp_global, &
         & r_xysp_local, mpi_dtype_yxsp_slice, &
         & COMM_SP_NE_X_NE, .true., root_is_in_slice)
    end if
    
    ! as long as in general only the root process can do IO: send
    ! output data to root process (this is done by a reduce_sum) and
    ! let root do the IO.

    if (root_processor) then

      ! NOTE with the rewrite of the io, both the binary
      ! and the ascii output produced by the following line are in 'F'
      ! order. Before, the ascii output used to be in 'C' order.
      call append_chunk(lun, r_xysp_global, xy_fmt, binary_fmt)
    end if

  end subroutine xysp_output_array

  !****************************************************************************
  !> This routine outputs real array in fourier space.
  !> It includes the gather in x when required.
  !> It produces different output for spectral (ky(kx)) and nonspectral (ky(xgr)).
  !> 
  !> Call with *all* procs, and pass proc_is_in_slice=T on those processes
  !> which hold their local data in the array argument.
  !----------------------------------------------------------------------------
  subroutine kyx_output_array(array, lun, &
     & proc_is_in_slice, root_is_in_slice)
    use io,           only : xy_fmt, binary_fmt
    use grid,         only : nmod,nx, n_x_grid
    use mpiinterface, only : gather_array, mpireduce_sum_inplace, root_processor
    use mpiinterface, only : proc_is_with_root_in
    use mpicomms,     only : COMM_DUMMY, COMM_X_NE, COMM_X_EQ

    ! for output in a single file:
    use io, only : append_chunk
    use control, only : io_legacy

    integer, intent(in) :: lun
    real, intent(in) :: array(nmod,nx)
    logical, intent(in) :: proc_is_in_slice, root_is_in_slice
    !> a buffer for the global ky-x slice
    real :: kyx_buf_global(nmod,n_x_grid)
    
    if(io_legacy) then
      if(proc_is_in_slice) then
        !gather global array in x (allgather)
        call gather_array(kyx_buf_global(1:nmod,1:n_x_grid),nmod,n_x_grid, &
           & array(1:nmod,1:nx),nmod,nx, &
           & COMM_DUMMY,COMM_X_NE,ALLGATHER=.true.)
      else
        kyx_buf_global = 0
      end if

      ! sum over the communicator which is complementary to COMM_X_NE,
      ! this effectively sends the array
      ! to the root process.
      if(proc_is_with_root_in(COMM_X_EQ)) then
        call mpireduce_sum_inplace(kyx_buf_global, &
           & shape(kyx_buf_global), &
           & COMM_X_EQ)
      end if
    else
      if(proc_is_in_slice .or. root_processor) then
        call gather_array(kyx_buf_global, &
           & array, mpi_dtype_kyx_slice, &
           & COMM_X_NE, .true., root_is_in_slice)
      end if
    end if
    if (root_processor) then

      ! NOTE with the rewrite of the io, both the binary
      ! and the ascii output produced by the following line are in 'F'
      ! order. Before, the ascii output used to be in 'C' order.
      call append_chunk(lun, kyx_buf_global, xy_fmt, binary_fmt)

    end if
  end subroutine kyx_output_array


  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine kyxsp_output_array(array, lun, &
     & proc_is_in_slice, root_is_in_slice)
    use io,           only : xy_fmt, binary_fmt
    use grid,         only : nmod,nx, nsp, n_x_grid, number_of_species
    use mpiinterface, only : gather_array, root_processor
    use mpicomms,     only : COMM_SP_NE_X_NE

    ! for output in a single file:
    use io, only : append_chunk

    integer, intent(in) :: lun
    real, intent(in) :: array(nmod,nx,nsp)
    logical, intent(in) :: proc_is_in_slice, root_is_in_slice
    !> a buffer for the global ky-x slice
    real :: kyxsp_buf_global(nmod,n_x_grid,number_of_species)

    if(proc_is_in_slice .or. root_processor) then
      call gather_array(kyxsp_buf_global, &
         & array, mpi_dtype_kyxsp_slice, &
         & COMM_SP_NE_X_NE, .true., root_is_in_slice)
    end if
    if (root_processor) then

      !kyx_buf_global = kyx_buf_global / n_procs_x

      ! NOTE with the rewrite of the io, both the binary
      ! and the ascii output produced by the following line are in 'F'
      ! order. Before, the ascii output used to be in 'C' order.
      call append_chunk(lun, kyxsp_buf_global, xy_fmt, binary_fmt)

    end if
  end subroutine kyxsp_output_array



  !****************************************************************************
  !> Transforms 2D arrays in x,y from k space to real space
  !> Selects the correct transform based on spectral_radius or not
  !> Note that the input array is destroyed by the FFT.
  !----------------------------------------------------------------------------
  subroutine four2real_2D(r_out,c_in,dir)

    use control, only  : spectral_radius
    use fft,     only  : four2D_real, four1D_real

    !> complex array (input for forward)
    complex, intent(inout)  :: c_in(:,:)
    !> real_array (output for forward)
    real,    intent(inout)  :: r_out(:,:)
    !> direction of transform (FFT_INVERSE or FFT_FORWARD)
    integer, intent(in)     :: dir
    integer :: ix

    if (spectral_radius) then
      call four2D_real(r_out,c_in,dir) 
      !scale the FFT by sqrt(n), not needed, values are independent of gridsize
      !r_out=r_out/sqrt(real(mphit*mrad_G))
    else
      !  FJC: replace with fourrow_real once working in nl_terms
      do ix = 1, size(r_out,2)!mrad_l
        call four1D_real(r_out(:,ix),c_in(:,ix),dir)
      end do
      !scale the FFT by sqrt(n), not needed, values are independant of gridsize
      !r_out=r_out/sqrt(real(mphit))
    end if

  end subroutine four2real_2D

  !--------------------------------------------------------------------------
  !> This routine calculates the s derivative of the gyroaveraged or
  !> not-gyroaveraged fields. it
  !> uses the differential scheme introduced in linear_terms.F90.
  !>
  !> Important: each diagnostic which uses this has to specify
  !> requirements(PHI_GA, LOCAL_DATA) = .true
  !> requirements(PHI_GA, S_GHOSTCELLS) = .true
  !> in its init() routine (and analogous for other fields).
  !-------------------------------------------------------------------------
  subroutine dfieldds(ifield,imod,ix,j,is,dfield)
    use linear_terms, only : differential_scheme
    use matdat, only : pos_par_grid, connect_parallel, set_indx
    use structures, only : matrix_element
    use fields, only : get_averaged_apar,get_averaged_bpar,get_averaged_phi
    use dist, only : iphi, iapar, ibpar, fdis_tmp, stencil_side, stencil_side_zf
    use index_function, only : indx
    use mode, only : parallel_phase_shift
    use geom, only : sgr_dist
    use grid, only : ns
    use general,only : gkw_abort
    use global, only : r_tiny, id_s
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    !>  is one of PHI_FIELD,PHI_GA_FIELD,APAR_FIELD,etc.
    integer, intent(in) :: ifield
    !> coordinates of of the field line. mu index j and species index
    !> are only relevant if a gyroavg is requested and ignored
    !> otherwise.
    integer, intent(in) :: imod,ix,j,is
    ! the s derivative of the field
    complex, intent(out) :: dfield(ns)
    ! Note that this buffer stores the local field and a few
    ! cells at both ends. These cells are not exactly ghost cells
    ! (which would not be needed if there was no parallelisation), but
    ! they are convenient:
    ! In case of spectral_radius=T the open
    ! boundary at the end of the line is implemented using one-sided
    ! stencils next to it. In contrast to this, in the case of
    ! spectral_radius=F there are no open boundaries and only
    ! symmetric stencils are used. Then a buffer like fields(:) needs 
    ! neighbourcells everywhere, regardless of the parallelisation.
    complex, allocatable :: field(:)
    type (matrix_element) :: elem
    integer :: ist, ipw, id, ierr, m, i, k
    logical ::ingrid
    complex :: dum
    real, allocatable :: w(:)
    allocate(w(1 + 4*max(stencil_side(id_s), stencil_side_zf(id_s))))

    if(.not. allocated(field)) then
      allocate(field(1-stencil_side(id_s):ns+stencil_side(id_s)), stat=ierr)
      if(ierr.ne.0) call gkw_abort('could not allocate field')
    end if

    !Fills an array with the field along the field-line
    !Array extended as it needs the ghost points for the
    !four point stencil
    do i=1-stencil_side(id_s), ns+stencil_side(id_s)
      call set_indx(elem,imod,ix,i,j,k,is)
      call connect_parallel(elem,ingrid)

      if(ingrid)then
        select case (ifield)
        case(PHI_FIELD)
          field(i) = fdis_tmp(indx(iphi,elem%imod,elem%ixloc,elem%iloc))
        case(APAR_FIELD)
          field(i) = fdis_tmp(indx(iapar,elem%imod,elem%ixloc,elem%iloc))
        case(BPAR_FIELD)
          field(i) = fdis_tmp(indx(ibpar,elem%imod,elem%ixloc,elem%iloc))
        case(PHI_GA_FIELD)
          field(i) = get_averaged_phi(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(APAR_GA_FIELD)
          field(i) = get_averaged_apar(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(BPAR_GA_FIELD)
          field(i) = get_averaged_bpar(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case default
          call gkw_abort('unknown ifield in diagnos_generic::dfieldds')
        end select
      end if
    end do

    ipw = -1
    id = + 1

    do i=1,ns
      dfield(i) = (0.0, 0.0)
      
      ! select the scheme 
      ist = pos_par_grid(imod,ix,i,k)
      call differential_scheme(ist,ipw,id,w)
      
      do m = 1, size(w)
        if (abs(w(m)) > r_tiny) then 
          call set_indx(elem,imod,ix,i,j,k,is)
          elem%iloc= i + m - ((size(w)+1)/2)
          dum = w(m) * field(elem%iloc) / sgr_dist * &
             & parallel_phase_shift(imod,ix,i,elem%iloc)
          ! note that connect_parallel was already used to fill the
          ! field(:) buffer. Call it again to find out if the current
          ! point is in the grid.
          call connect_parallel(elem,ingrid)
          if(ingrid) dfield(i) = dfield(i) + dum
        endif
      enddo
    enddo

    deallocate(w)
    
  end subroutine dfieldds

  !****************************************************************************
  !> return the real space value of gyro-averaged electrostatic
  !> potential at a specific psi zeta point 
  !> the value of <phi> array for vpar, mu, s, sp, zeta, psi is returned
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !> The points i,j,k,is are local
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function phiga_xy_point_get(psi,zeta,i,j,is)
    use dist,             only : fdisi
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : spectral_radius, flux_tube
    use grid,             only : nx, nmod, lx, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use matdat,           only : get_f_from_g
    use general,          only : gkw_abort
    use fields,           only : get_averaged_phi
    integer, intent(in) :: i,j,is
    real, intent(in) :: psi,zeta

    real                  :: phiga_xy_point_get
    !real :: ix_mid
    integer               :: ix, ixg, imod, helpint
    complex :: dum
    complex               :: dfield

    dfield=(0.0E0,0.0E0)
    do imod = 1, nmod

      if (spectral_radius) then    !Manually do the fourier sum in 2D
        do ix = 1, nx
          dum = get_averaged_phi(imod,ix,i,j,is,fdisi)

          ! This remaps the island to 0,0
          helpint=ix-1-(nx-1)/2

          ! should have used krho here
          dfield = dfield + 2 * dum *&
             & exp(ci1*2*pi/lyn*(imod-1)*zeta + &
             & ci1*2*pi*helpint/lx*psi)
        end do
      else
        ! Manually do the Fourier sum in 1D, choose the nearest point in x .
        ! Note that island is at centre in nonspectral.
        !ix_mid = real((n_x_grid+1)*0.5E0)      
        !ix=nint(ix_mid+psi/lx*real(n_x_grid))
        if (flux_tube) then
          ! We use the global number of radial points n_x_grid to determine the
          ! global ix for the chosen radial position
          ixg = nint(psi/lxn*real(n_x_grid))
        else
          ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
        end if

        if (ixg < 1 .or. ixg > n_x_grid) then      
          call gkw_abort('nonspectral: Error in phiga_xy_point_get')
        end if

        ix  = lrx(ixg)

        if (proc_subset(ixg,0,0,0,0)) then
          ! dum is the distribution without A_par contribution
          ! The derivative is calculated from spectral representation &
          ! to real space via an inverse Fourier transform
          dum = get_averaged_phi(imod,ix,i,j,is,fdisi)
          ! should have used krho here
          dfield = dfield + 2*dum*exp((ci1*2*pi/lyn*(imod-1))*zeta)   
        else ! psi point not on this processor
          dfield = 0.0
        end if
      end if
    end do

    phiga_xy_point_get = real(dfield)

  end function phiga_xy_point_get

  !****************************************************************************
  !> return the real space value of grad of gyro-averaged electrostatic
  !> potential at a specific psi zeta point.
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !>
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dphigads_xy_point_get(psi,zeta,i,j,is)
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : flux_tube
    use grid,             only : ns
    use grid,             only : nmod, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use general,          only : gkw_abort
    use global, only : PHI_GA_FIELD

    !> local parallel, mu, and species index
    integer, intent(in) :: i,j,is
    !> radial and binormal position space coordinates (not indices!)
    real, intent(in) :: psi,zeta
    !> the value of d<phi>/ds at the specified position is returned
    real :: dphigads_xy_point_get

    !real ::ix_mid
    integer               :: ix, ixg, imod
    complex               :: dfield
    complex               :: dfieldds_line(ns)

    dfield= (0.0E0, 0.0E0)

    do imod = 1, nmod
      !ix_mid = real((n_x_grid+1)*0.5E0)
      !ix=nint(ix_mid+psi/lx*real(n_x_grid))
      if (flux_tube) then
        ! We use the global number of radial points n_x_grid to determine the
        ! global ix for the chosen radial position
        ixg = nint(psi/lxn*real(n_x_grid))
      else
        ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
      end if

      ix = lrx(ixg)

      call dfieldds(PHI_GA_FIELD,imod,ix,j,is,dfieldds_line)
      
      if (ixg < 1 .or. ixg > n_x_grid) then
        call gkw_abort('nonspectral: Error in dphigads_xy_point_get')
      end if

      if (proc_subset(ixg,0,0,0,0)) then
        dfield = dfield + 2*dfieldds_line(i) * &
           & exp(ci1*2*pi/lyn*(imod-1)*zeta)
      else ! psi point not on this processor
        dfield = 0.0
      end if
    end do

    dphigads_xy_point_get = real(dfield)

  end function dphigads_xy_point_get

  !****************************************************************************
  !> return the real space value of grad of gyro-averaged electrostatic
  !> potential at a specific psi zeta point.
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !>
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dphigadx_xy_point_get(psi,zeta,i,j,is)
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : spectral_radius, flux_tube
    use grid,             only : nx, nmod, lx, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use general,          only : gkw_abort
    use global, only : PHI_GA_FIELD

    !> local parallel, mu, and species index
    integer, intent(in) :: i,j,is
    !> radial and binormal position space coordinates (not indices!)
    real, intent(in) :: psi,zeta
    !> the value of d<phi>/dx at the specified position is returned
    real :: dphigadx_xy_point_get

    !real ::ix_mid
    integer               :: ix, ixg, imod, helpint
    complex               :: dfield
    complex :: dfielddx_line(nx)

    dfield= (0.0E0, 0.0E0)

    do imod = 1, nmod
      call dfielddx(PHI_GA_FIELD,imod,i,j,is,dfielddx_line)

      if (spectral_radius) then
        !Manually do the fourier sum in 2D
        do ix = 1, nx
          ! This remaps the island to 0,0
          helpint=ix-1-(nx-1)/2

          dfield = dfield + 2*dfielddx_line(ix) *&
             & exp(ci1*2*pi/lyn*(imod-1)*zeta + &
             & ci1*2*pi*helpint/lx*psi)
        end do
      else
        ! Manually do the Fourier sum in 1D, choose the nearest point in x .
        ! Note that island is at centre in nonspectral.
        !ix_mid = real((n_x_grid+1)*0.5E0)
        !ix=nint(ix_mid+psi/lx*real(n_x_grid))
        if (flux_tube) then
          ! We use the global number of radial points n_x_grid to determine the
          ! global ix for the chosen radial position
          ixg = nint(psi/lxn*real(n_x_grid))
        else
          ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
        end if

        ix = lrx(ixg)

        if (ix < 1 .or. ix > n_x_grid) then
          call gkw_abort('nonspectral: Error in dphigadx_xy_point_get')
        end if

        if (proc_subset(ixg,0,0,0,0)) then
          dfield = dfield + 2*dfielddx_line(ix) * &
             & exp(ci1*2*pi/lyn*(imod-1)*zeta)
        else ! psi point not on this processor
          dfield = 0.0
        end if
      end if
    end do

    dphigadx_xy_point_get = real(dfield)

  end function dphigadx_xy_point_get

  !****************************************************************************
  !> return the real space value of grad of gyro-averaged electrostatic
  !> potential at a specific psi zeta point 
  !> the value of d<phi>/dzeta array for vpar, mu, s, sp, zeta, psi is returned
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !> The points i,j,k,is are local
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dphigadzeta_xy_point_get(psi,zeta,i,j,is)
    use dist,             only : fdisi
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : spectral_radius, flux_tube
    use grid,             only : nx, nmod, lx, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use general,          only : gkw_abort
    use fields,           only : get_averaged_phi

    integer, intent(in) :: i,j,is
    real, intent(in) :: psi,zeta

    real                  :: dphigadzeta_xy_point_get
    !real                  :: ix_mid
    integer               :: ix, ixg, imod, helpint
    complex               :: dum
    complex               :: dfield

    dfield= (0.0E0, 0.0E0)
    do imod = 1, nmod

      if (spectral_radius) then    !Manually do the fourier sum in 2D
        do ix = 1, nx
          dum = get_averaged_phi(imod,ix,i,j,is,fdisi)

          ! This remaps the island to 0,0
          helpint=ix-1-(nx-1)/2

          ! should have used krho here
          dfield = dfield + 2 * dum * (ci1*2*pi/lyn*(imod-1)) *&
             & exp(ci1*2*pi/lyn*(imod-1)*zeta + &
             & ci1*2*pi*helpint/lx*psi)
        end do
      else
        ! Manually do the Fourier sum in 1D, choose the nearest point in x .
        ! Note that island is at centre in nonspectral.
        !ix_mid = real((n_x_grid+1)*0.5E0)      
        !ix=nint(ix_mid+psi/lx*real(n_x_grid))
        if (flux_tube) then
          ! We use the global number of radial points n_x_grid to determine the
          ! global ix for the chosen radial position
          ixg = nint(psi/lxn*real(n_x_grid))
        else
          ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
        end if

        if (ixg < 1 .or. ixg > n_x_grid) then      
          call gkw_abort('nonspectral: Error in dphigadzeta_xy_point_get')
        end if

        ix  = lrx(ixg)

        if (proc_subset(ixg,0,0,0,0)) then
          ! dum is the distribution without A_par contribution
          ! The derivative is calculated from spectral representation &
          ! to real space via an inverse Fourier transform
          dum = get_averaged_phi(imod,ix,i,j,is,fdisi)
          ! should have used krho here
          dfield = dfield + (ci1*2*pi/lyn*(imod-1))*2*dum*exp((ci1*2*pi/lyn*(imod-1))*zeta)   
        else ! psi point not on this processor
          dfield = 0.0
        end if
      end if
    end do

    dphigadzeta_xy_point_get = real(dfield)

  end function dphigadzeta_xy_point_get

  !****************************************************************************
  !> return the real space value of parallel derivative of gyro-averaged A||
  !> potential at a specific psi zeta point.
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !>
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dapargads_xy_point_get(psi,zeta,i,j,is)
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : flux_tube
    use grid,             only : ns
    use grid,             only : nmod, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use matdat,           only : get_f_from_g
    use general,          only : gkw_abort
    use fields,           only : get_averaged_phi
    use matdat,           only : connect_rad
    use structures,       only : matrix_element
    use linear_terms,     only : differential_scheme
    use global,           only : APAR_GA_FIELD
    !> local parallel, mu, and species index
    integer, intent(in) :: i,j,is
    !> radial and binormal position space coordinates (not indices!)
    real, intent(in) :: psi,zeta
    !> the value of d<phi>/ds at the specified position is returned
    real :: dapargads_xy_point_get

    !real :: ix_mid
    integer               :: ix, ixg, imod
    complex               :: dfield
    complex               :: dfieldds_line(ns)

    dfield=(0.0E0,0.0E0)

    do imod = 1, nmod
      !ix_mid = real((n_x_grid+1)*0.5E0)
      !ix=nint(ix_mid+psi/lx*real(n_x_grid))
      if (flux_tube) then
        ! We use the global number of radial points n_x_grid to determine the
        ! global ix for the chosen radial position
        ixg = nint(psi/lxn*real(n_x_grid))
      else
        ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
      end if

      ix = lrx(ixg)

      call dfieldds(APAR_GA_FIELD,imod,ix,j,is,dfieldds_line)
      
      if (ixg < 1 .or. ixg > n_x_grid) then
        call gkw_abort('nonspectral: Error in dapargads_xy_point_get')
      end if

      if (proc_subset(ixg,0,0,0,0)) then
        dfield = dfield + 2*dfieldds_line(i) * &
           & exp(ci1*2*pi/lyn*(imod-1)*zeta)
      else ! psi point not on this processor
        dfield = 0.0
      end if
    end do

    dapargads_xy_point_get = dfield

  end function dapargads_xy_point_get

  !****************************************************************************
  !> return the real space value of grad of gyro-averaged A||
  !> potential at a specific psi zeta point.
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !>
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dapargadx_xy_point_get(psi,zeta,i,j,is)
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : spectral_radius, flux_tube
    use grid,             only : nx, nmod, lx, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use matdat,           only : get_f_from_g
    use general,          only : gkw_abort
    use fields,           only : get_averaged_phi
    use matdat,           only : connect_rad
    use structures,       only : matrix_element
    use linear_terms,     only : differential_scheme
    use global,           only : APAR_GA_FIELD
    !> local parallel, mu, and species index
    integer, intent(in) :: i,j,is
    !> radial and binormal position space coordinates (not indices!)
    real, intent(in) :: psi,zeta
    !> the value of d<phi>/dx at the specified position is returned
    real :: dapargadx_xy_point_get

    !real :: ix_mid
    integer               :: ix, ixg, imod, helpint
    complex               :: dfield
    complex :: dfielddx_line(nx)

    dfield=(0.0E0,0.0E0)

    do imod = 1, nmod
      call dfielddx(APAR_GA_FIELD,imod,i,j,is,dfielddx_line)

      if (spectral_radius) then
        !Manually do the fourier sum in 2D
        do ix = 1, nx

          ! This remaps the island to 0,0
          helpint=ix-1-(nx-1)/2

          dfield = dfield + 2*dfielddx_line(ix) *&
             & exp(ci1*2*pi/lyn*(imod-1)*zeta + &
             & ci1*2*pi*helpint/lx*psi)
        end do
      else
        ! Manually do the Fourier sum in 1D, choose the nearest point in x .
        ! Note that island is at centre in nonspectral.
        !ix_mid = real((n_x_grid+1)*0.5E0)
        !ix=nint(ix_mid+psi/lx*real(n_x_grid))
        if (flux_tube) then
          ! We use the global number of radial points n_x_grid to determine the
          ! global ix for the chosen radial position
          ixg = nint(psi/lxn*real(n_x_grid))
        else
          ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
        end if

        ix = lrx(ixg)

        if (ix < 1 .or. ix > n_x_grid) then
          call gkw_abort('nonspectral: Error in dapargadx_xy_point_get')
        end if

        if (proc_subset(ixg,0,0,0,0)) then
          dfield = dfield + 2*dfielddx_line(ix) * &
             & exp(ci1*2*pi/lyn*(imod-1)*zeta)
        else ! psi point not on this processor
          dfield = 0.0
        end if
      end if
    end do

    dapargadx_xy_point_get = dfield

  end function dapargadx_xy_point_get

  !****************************************************************************
  !> return the real space value of grad of gyro-averaged A||
  !> potential at a specific psi zeta point 
  !> the value of d<phi>/dzeta array for vpar, mu, s, sp, zeta, psi is returned
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !> The points i,j,k,is are local
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dapargadzeta_xy_point_get(psi,zeta,i,j,is)
    use dist,             only : fdisi
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : spectral_radius, flux_tube
    use grid,             only : nx, nmod, lx, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use matdat,           only : get_f_from_g
    use general,          only : gkw_abort
    use fields,           only : get_averaged_apar
    integer, intent(in) :: i,j,is
    real, intent(in) :: psi,zeta
    real                  :: dapargadzeta_xy_point_get
    !real :: ix_mid
    integer               :: ix, ixg, imod, helpint
    complex :: dum
    complex               :: dfield
    

    dfield=(0.0E0,0.0E0)
    do imod = 1, nmod

      if (spectral_radius) then    !Manually do the fourier sum in 2D
        do ix = 1, nx
          dum = get_averaged_apar(imod,ix,i,j,is,fdisi)

          ! This remaps the island to 0,0
          helpint=ix-1-(nx-1)/2

          ! should have used krho here
          dfield = dfield + 2 * dum * (ci1*2*pi/lyn*(imod-1)) *&
             & exp(ci1*2*pi/lyn*(imod-1)*zeta + &
             & ci1*2*pi*helpint/lx*psi)
        end do
      else
        ! Manually do the Fourier sum in 1D, choose the nearest point in x .
        ! Note that island is at centre in nonspectral.
        !ix_mid = real((n_x_grid+1)*0.5E0)      
        !ix=nint(ix_mid+psi/lx*real(n_x_grid))
        if (flux_tube) then
          ! We use the global number of radial points n_x_grid to determine the
          ! global ix for the chosen radial position
          ixg = nint(psi/lxn*real(n_x_grid))
        else
          ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
        end if

        if (ixg < 1 .or. ixg > n_x_grid) then      
          call gkw_abort('nonspectral: Error in daparga_xy_point_get')
        end if

        ix  = lrx(ixg)

        if (proc_subset(ixg,0,0,0,0)) then
          ! dum is the distribution without A_par contribution
          ! The derivative is calculated from spectral representation &
          ! to real space via an inverse Fourier transform
          dum = get_averaged_apar(imod,ix,i,j,is,fdisi)
          ! should have used krho here
          dfield = dfield + (ci1*2*pi/lyn*(imod-1))*2*dum*exp((ci1*2*pi/lyn*(imod-1))*zeta)   
        else ! psi point not on this processor
          dfield = 0.0
        end if
      end if
    end do

    dapargadzeta_xy_point_get = real(dfield)

  end function dapargadzeta_xy_point_get

  !--------------------------------------------------------------------------
  !> This routine calculates the 1st x derivative of the gyroaveraged or
  !> not-gyroaveraged fields. it
  !> uses the differential scheme introduced in linear_terms.F90.
  !>
  !> Important: each diagnostic which uses this has to specify
  !> requirements(PHI_GA, LOCAL_DATA) = .true
  !> requirements(PHI_GA, X_GHOSTCELLS) = .true
  !> in its init() routine (and analogous for other fields).
  !-------------------------------------------------------------------------
  subroutine dfielddx(ifield,imod,i,j,is,dfield)
    use linear_terms, only : differential_scheme
    use matdat, only : set_indx, connect_rad
    use structures, only : matrix_element
    use fields, only : get_averaged_apar,get_averaged_bpar,get_averaged_phi
    use dist, only : iphi, iapar, ibpar, fdisi, fdis_tmp
    use dist, only : stencil_side, stencil_side_zf
    use index_function, only : indx
    use mode, only : kxrh
    use geom, only : dxgr
    use grid, only : nx, gx
    use general,only : gkw_abort
    use global, only : r_tiny, id_x
    use constants, only : ci1
    use control, only : spectral_radius, order_of_the_radial_scheme
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    
    !>  is one of PHI_FIELD,PHI_GA_FIELD,APAR_FIELD,etc.
    integer, intent(in) :: ifield
    !> coordinates of of the field line. mu index j and species index
    !> are only relevant if a gyroavg is requested and ignored
    !> otherwise.
    integer, intent(in) :: imod,i,j,is
    ! the 2nd x derivative of the field
    complex, intent(out):: dfield(nx)
    complex, allocatable :: field(:)
    type (matrix_element) :: elem
    integer :: ist, ipw, id, ierr, m, k
    integer :: ix
    logical ::ingrid
    complex :: dum
    real, allocatable :: w(:)

    allocate(w(1 + 4*max(stencil_side(id_x), stencil_side_zf(id_x))))
    
    if(.not. allocated(field)) then
      allocate(field(1-stencil_side(id_x):nx+stencil_side(id_x)), stat=ierr)
      if(ierr.ne.0) call gkw_abort('could not allocate field')
    end if
    !Fills an array with the field along the field-line
    !Array extended as it needs the ghost points for the
    !four point stencil
    if (spectral_radius) then
      do ix = 1, nx
        call set_indx(elem,imod,ix,i,j,k,is)
        select case (ifield)
        case(PHI_FIELD)
          field(ix) = fdisi(indx(iphi,elem%imod,elem%ixloc,elem%iloc))
        case(APAR_FIELD)
          field(ix) = fdisi(indx(iapar,elem%imod,elem%ixloc,elem%iloc))
        case(BPAR_FIELD)
          field(ix) = fdisi(indx(ibpar,elem%imod,elem%ixloc,elem%iloc))
        case(PHI_GA_FIELD)
          field(ix)=get_averaged_phi(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(APAR_GA_FIELD)
          field(ix)=get_averaged_apar(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(BPAR_GA_FIELD)
          call gkw_abort('check if fdis_tmp or fdis_tmp2 is correct&
             & in diagnos_generic::dfieldds')
          field(ix)=get_averaged_bpar(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(9999)
          ! test
          field(ix) = 1.0
        case default
          call gkw_abort('unknown ifield in diagnos_generic::dfielddx')
        end select
      end do
    else
      do ix=1-stencil_side(id_x), nx+stencil_side(id_x)
        call set_indx(elem,imod,ix,i,j,k,is)
        elem%val=(0.,0.)
        call connect_rad(elem,ingrid)

        if(ingrid)then
          select case (ifield)
          case(PHI_FIELD)
            field(ix) = fdis_tmp(indx(iphi,elem%imod,elem%ixloc,elem%iloc))
          case(APAR_FIELD)
            field(ix) = fdis_tmp(indx(iapar,elem%imod,elem%ixloc,elem%iloc))
          case(BPAR_FIELD)
            field(ix) = fdis_tmp(indx(ibpar,elem%imod,elem%ixloc,elem%iloc))
          case(PHI_GA_FIELD)
            field(ix)=get_averaged_phi(elem%imod,elem%ixloc,elem%iloc, &
               & elem%jloc,elem%isloc,fdis_tmp)
          case(APAR_GA_FIELD)
            field(ix)=get_averaged_apar(elem%imod,elem%ixloc,elem%iloc, &
               & elem%jloc,elem%isloc,fdis_tmp)
          case(BPAR_GA_FIELD)
            call gkw_abort('check if fdis_tmp or fdis_tmp2 is correct&
               & in diagnos_generic::dfieldds')
            field(ix)=get_averaged_bpar(elem%imod,elem%ixloc,elem%iloc, &
               & elem%jloc,elem%isloc,fdis_tmp)
          case(9999)
            ! test
            field(ix) = 3*gx(ix)*dxgr
          case default
            call gkw_abort('unknown ifield in diagnos_generic::dfielddx')
          end select
        end if
      end do
    end if

    ipw = -1
    ! pick 1st deriv
    id = 1
    ist = 0
    do ix=1,nx
      if (spectral_radius) then
        call set_indx(elem,imod,ix,i,j,k,is)
        dfield(ix) = (ci1*kxrh(ix)) * field(elem%ixloc)
      else
        dfield(ix) = 0.0
        ! select the scheme
        call differential_scheme(ist,ipw,id,w,order_of_the_radial_scheme)

        do m = 1, size(w)
          if (abs(w(m)) > r_tiny) then
            call set_indx(elem,imod,ix,i,j,k,is)
            elem%ixloc= ix + m - ((size(w)+1)/2)
            ! note that connect_rad was already used to fill the
            ! field(:) buffer. Call it again to find out if the current
            ! point is in the grid.
            call connect_rad(elem,ingrid)
            if(ingrid) then
              dum = w(m) * field(elem%ixloc) / dxgr
              dfield(ix) = dfield(ix) + dum
             endif
          endif
        enddo
      end if
    enddo

    deallocate(w)
  end subroutine dfielddx

  !--------------------------------------------------------------------------
  !> This routine calculates the 2nd x derivative of the gyroaveraged or
  !> not-gyroaveraged fields. it
  !> uses the differential scheme introduced in linear_terms.F90.
  !>
  !> Important: each diagnostic which uses this has to specify
  !> requirements(PHI_GA, LOCAL_DATA) = .true
  !> requirements(PHI_GA, X_GHOSTCELLS) = .true
  !> in its init() routine (and analogous for other fields).
  !-------------------------------------------------------------------------
  subroutine d2fielddx2(ifield,imod,i,j,is,dfield)
    use linear_terms, only : differential_scheme
    use matdat, only : set_indx, connect_rad
    use structures, only : matrix_element
    use fields, only : get_averaged_apar,get_averaged_bpar,get_averaged_phi
    use dist, only : iphi, iapar, ibpar, fdisi, fdis_tmp
    use dist, only : stencil_side, stencil_side_zf
    use index_function, only : indx
    use mode, only : kxrh
    use geom, only : dxgr
    use grid, only : nx, gx
    use general,only : gkw_abort
    use global, only : r_tiny, id_x
    use constants, only : ci1
    use control, only : spectral_radius, order_of_the_radial_scheme
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    !>  is one of PHI_FIELD,PHI_GA_FIELD,APAR_FIELD,etc.
    integer, intent(in) :: ifield
    !> coordinates of of the field line. mu index j and species index
    !> are only relevant if a gyroavg is requested and ignored
    !> otherwise.
    integer, intent(in) :: imod,i,j,is
    ! the 2nd x derivative of the field
    complex, intent(out):: dfield(nx)
    complex, allocatable :: field(:)
    type (matrix_element) :: elem
    integer :: ist, ipw, id, ierr, m, k
    integer :: ix
    logical ::ingrid
    complex :: dum
    real, allocatable :: w(:)
    allocate(w(1 + 4*max(stencil_side(id_x), stencil_side_zf(id_x))))
    
    if(.not. allocated(field)) then
      allocate(field(1-stencil_side(id_x):nx+stencil_side(id_x)), stat=ierr)
      if(ierr.ne.0) call gkw_abort('could not allocate field')
    end if
    !Fills an array with the field along the field-line
    !Array extended as it needs the ghost points for the
    !four point stencil
    if (spectral_radius) then
      do ix = 1, nx
        call set_indx(elem,imod,ix,i,j,k,is)
        select case (ifield)
        case(PHI_FIELD)
          field(ix) = fdisi(indx(iphi,elem%imod,elem%ixloc,elem%iloc))
        case(APAR_FIELD)
          field(ix) = fdisi(indx(iapar,elem%imod,elem%ixloc,elem%iloc))
        case(BPAR_FIELD)
          field(ix) = fdisi(indx(ibpar,elem%imod,elem%ixloc,elem%iloc))
        case(PHI_GA_FIELD)
          field(ix)=get_averaged_phi(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(APAR_GA_FIELD)
          call gkw_abort('check if fdis_tmp or fdis_tmp2 is correct&
             & in diagnos_generic::dfieldds')
          field(ix)=get_averaged_apar(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(BPAR_GA_FIELD)
          call gkw_abort('check if fdis_tmp or fdis_tmp2 is correct&
             & in diagnos_generic::dfieldds')
          field(ix)=get_averaged_bpar(elem%imod,elem%ixloc,elem%iloc, &
             & elem%jloc,elem%isloc,fdis_tmp)
        case(9999)
          ! test
          field(ix) = 1.0
        case default
          call gkw_abort('unknown ifield in diagnos_generic::d2fielddx2')
        end select
      end do
    else
      do ix=1-stencil_side(id_x), nx+stencil_side(id_x)
        call set_indx(elem,imod,ix,i,j,k,is)
        call connect_rad(elem,ingrid)

        if(ingrid)then
          select case (ifield)
          case(PHI_FIELD)
            field(ix) = fdis_tmp(indx(iphi,elem%imod,elem%ixloc,elem%iloc))
          case(APAR_FIELD)
            field(ix) = fdis_tmp(indx(iapar,elem%imod,elem%ixloc,elem%iloc))
          case(BPAR_FIELD)
            field(ix) = fdis_tmp(indx(ibpar,elem%imod,elem%ixloc,elem%iloc))
          case(PHI_GA_FIELD)
            field(ix)=get_averaged_phi(elem%imod,elem%ixloc,elem%iloc, &
               & elem%jloc,elem%isloc,fdis_tmp)
          case(APAR_GA_FIELD)
            call gkw_abort('check if fdis_tmp or fdis_tmp2 is correct&
               & in diagnos_generic::dfieldds')
            field(ix)=get_averaged_apar(elem%imod,elem%ixloc,elem%iloc, &
               & elem%jloc,elem%isloc,fdis_tmp)
          case(BPAR_GA_FIELD)
            call gkw_abort('check if fdis_tmp or fdis_tmp2 is correct&
               & in diagnos_generic::dfieldds')
            field(ix)=get_averaged_bpar(elem%imod,elem%ixloc,elem%iloc, &
               & elem%jloc,elem%isloc,fdis_tmp)
          case(9999)
            ! test
            field(ix) = 3*(gx(ix)*dxgr)**2
          case default
            call gkw_abort('unknown ifield in diagnos_generic::d2fielddx2')
          end select
        end if
      end do
    end if

    ipw = -1
    ! do not pick 4th deriv, but always 2nd deriv, regardless of the order of the scheme.
    id = 3
    ist = 0
    do ix=1,nx
      if (spectral_radius) then
        call set_indx(elem,imod,ix,i,j,k,is)
        dfield(ix) = (ci1*kxrh(ix))**2 * field(elem%ixloc)
      else
        dfield(ix) = 0.0
        ! select the scheme
        call differential_scheme(ist,ipw,id,w,order_of_the_radial_scheme)

        do m = 1, size(w)
          if (abs(w(m)) > r_tiny) then
            call set_indx(elem,imod,ix,i,j,k,is)
            elem%ixloc= ix + m - ((size(w)+1)/2)
            ! note that connect_rad was already used to fill the
            ! field(:) buffer. Call it again to find out if the current
            ! point is in the grid.
            call connect_rad(elem,ingrid)
            if(ingrid) then
              dum = w(m) * field(elem%ixloc) / (dxgr**2)
              dfield(ix) = dfield(ix) + dum
            endif
          endif
        enddo
      end if
    enddo

    deallocate(w)
  end subroutine d2fielddx2

  !--------------------------------------------------------------------
  !> This function returns the integer that is used in various diagnostics
  !> to distinguish trapped and passing particles.
  !>
  !> Put a factor (-mask+1) into 'passing' quantities,
  !> and a factor (mask) into trapped ones.
  !>
  !--------------------------------------------------------------------
  pure function get_tr_ps_mask(ix,i,j,k,is) result(mask)
    use velocitygrid, only : mugr, vpgr
    use geom, only : bn, bmax
    use rotation, only : cfen_hfs, cfen
    !> radial,parallel,mu,vpar and species index
    integer, intent(in) :: ix,i,j,k,is
    integer :: mask
    ! real :: b_max_this_fieldline
    ! b_max_this_fieldline = maxval(bn_G(ix,:))
    
    if(mugr(j) > &
       & (vpgr(i,j,k,is)**2 +cfen(i,is)-cfen_hfs(is))/(2*(bmax(ix)-bn(ix,i)))) then
      ! trapped
      mask= 1
    else
      ! passing
      mask= 0
    end if

  end function get_tr_ps_mask


  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()

    if (allocated(c_xy_local)) deallocate(c_xy_local)
    if (allocated(r_xy_local)) deallocate(r_xy_local)
    if (allocated(r_xysp_local)) deallocate(r_xysp_local)
    if (allocated(global_vspace)) deallocate(global_vspace)

    if (allocated(r_xy_global)) deallocate(r_xy_global)
    if (allocated(r_xysp_global)) deallocate(r_xysp_global)

    if (allocated(parseval_correction)) deallocate(parseval_correction)
    
  end subroutine finalize

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_1d(lun, gridname1, preference_if_mixed)
    use io, only : attach_metadata
    integer, intent(in) :: lun
    character(len=*), intent(in) :: gridname1
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(lun, 'grid1', gridname1, preference_if_mixed)

  end subroutine attach_metadata_grid_1d

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_2d(lun, gridname1, gridname2, &
     & preference_if_mixed)
    use io, only : attach_metadata
    integer, intent(in) :: lun
    character(len=*), intent(in) :: gridname1, gridname2
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(lun, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(lun, 'grid2', gridname2, preference_if_mixed)

  end subroutine attach_metadata_grid_2d

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_3d(lun, gridname1, gridname2, gridname3, &
     & preference_if_mixed)
    use io, only : attach_metadata
    integer, intent(in) :: lun
    character(len=*), intent(in) :: gridname1, gridname2, gridname3
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(lun, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(lun, 'grid2', gridname2, preference_if_mixed)
    call attach_metadata(lun, 'grid3', gridname3, preference_if_mixed)

  end subroutine attach_metadata_grid_3d

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_4d(lun, &
     & gridname1, gridname2, gridname3, gridname4, &
     & preference_if_mixed)
    use io, only : attach_metadata
    integer, intent(in) :: lun
    character(len=*), intent(in) :: gridname1, gridname2, gridname3, gridname4
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(lun, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(lun, 'grid2', gridname2, preference_if_mixed)
    call attach_metadata(lun, 'grid3', gridname3, preference_if_mixed)
    call attach_metadata(lun, 'grid4', gridname4, preference_if_mixed)

  end subroutine attach_metadata_grid_4d

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_5d(lun, &
     & gridname1, gridname2, gridname3, gridname4, &
     & gridname5, &
     & preference_if_mixed)
    use io, only : attach_metadata
    integer, intent(in) :: lun
    character(len=*), intent(in) :: gridname1, gridname2, gridname3, gridname4
    character(len=*), intent(in) :: gridname5
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(lun, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(lun, 'grid2', gridname2, preference_if_mixed)
    call attach_metadata(lun, 'grid3', gridname3, preference_if_mixed)
    call attach_metadata(lun, 'grid4', gridname4, preference_if_mixed)
    call attach_metadata(lun, 'grid4', gridname5, preference_if_mixed)

  end subroutine attach_metadata_grid_5d

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_6d(lun, &
     & gridname1, gridname2, gridname3, gridname4, &
     & gridname5, gridname6, &
     & preference_if_mixed)
    use io, only : attach_metadata
    integer, intent(in) :: lun
    character(len=*), intent(in) :: gridname1, gridname2, gridname3, gridname4
    character(len=*), intent(in) :: gridname5, gridname6
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(lun, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(lun, 'grid2', gridname2, preference_if_mixed)
    call attach_metadata(lun, 'grid3', gridname3, preference_if_mixed)
    call attach_metadata(lun, 'grid4', gridname4, preference_if_mixed)
    call attach_metadata(lun, 'grid4', gridname5, preference_if_mixed)
    call attach_metadata(lun, 'grid6', gridname6, preference_if_mixed)

  end subroutine attach_metadata_grid_6d

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_7d(lun, &
     & gridname1, gridname2, gridname3, gridname4, &
     & gridname5, gridname6, gridname7, &
     & preference_if_mixed)
    use io, only : attach_metadata
    integer, intent(in) :: lun
    character(len=*), intent(in) :: gridname1, gridname2, gridname3, gridname4
    character(len=*), intent(in) :: gridname5, gridname6, gridname7
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(lun, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(lun, 'grid2', gridname2, preference_if_mixed)
    call attach_metadata(lun, 'grid3', gridname3, preference_if_mixed)
    call attach_metadata(lun, 'grid4', gridname4, preference_if_mixed)
    call attach_metadata(lun, 'grid4', gridname5, preference_if_mixed)
    call attach_metadata(lun, 'grid6', gridname6, preference_if_mixed)
    call attach_metadata(lun, 'grid7', gridname7, preference_if_mixed)

  end subroutine attach_metadata_grid_7d

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_1d_name(luname, groupname, gridname1, preference_if_mixed)
    use io, only : attach_metadata
    character(len=*), intent(in) :: luname, groupname
    character(len=*), intent(in) :: gridname1
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(luname, groupname, 'grid1', gridname1, preference_if_mixed)

  end subroutine attach_metadata_grid_1d_name

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_2d_name(luname, groupname, gridname1, gridname2, &
     & preference_if_mixed)
    use io, only : attach_metadata
    character(len=*), intent(in) :: luname, groupname
    character(len=*), intent(in) :: gridname1, gridname2
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(luname, groupname, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(luname, groupname, 'grid2', gridname2, preference_if_mixed)

  end subroutine attach_metadata_grid_2d_name

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_3d_name(luname, groupname, gridname1, gridname2, gridname3, &
     & preference_if_mixed)
    use io, only : attach_metadata
    character(len=*), intent(in) :: luname, groupname
    character(len=*), intent(in) :: gridname1, gridname2, gridname3
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(luname, groupname, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(luname, groupname, 'grid2', gridname2, preference_if_mixed)
    call attach_metadata(luname, groupname, 'grid3', gridname3, preference_if_mixed)

  end subroutine attach_metadata_grid_3d_name

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine attach_metadata_grid_4d_name(luname, groupname, &
     & gridname1, gridname2, gridname3, gridname4, &
     & preference_if_mixed)
    use io, only : attach_metadata
    character(len=*), intent(in) :: luname, groupname
    character(len=*), intent(in) :: gridname1, gridname2, gridname3, gridname4
    character (len=*), intent(in) :: preference_if_mixed

    call attach_metadata(luname, groupname, 'grid1', gridname1, preference_if_mixed)
    call attach_metadata(luname, groupname, 'grid2', gridname2, preference_if_mixed)
    call attach_metadata(luname, groupname, 'grid3', gridname3, preference_if_mixed)
    call attach_metadata(luname, groupname, 'grid4', gridname4, preference_if_mixed)

  end subroutine attach_metadata_grid_4d_name

end module diagnos_generic
