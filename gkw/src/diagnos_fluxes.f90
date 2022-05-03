!------------------------------------------------------------------------------
!> Calculates the fluxes of particles, energy and  toroidal angular momentum.
!> Calculates both the anomalous as well as the neo-classical fluxes. 
!> The anomalous flux is furthermore split into contributions due to the 
!> ExB velocity and the magnetic flutter and magnetic compression.
!------------------------------------------------------------------------------
module diagnos_fluxes

  implicit none

  private

  logical, save, public :: lcalc_fluxes

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output, screen_output

  public :: pflux_es, eflux_es, vflux_es

  !> the range of tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl

  ! logicals to control the fluxes output; default to T
  logical, save, public :: lfluxes_spectra
  logical, save, public :: lfluxes_em_spectra

  logical, save, public :: flux3d

  !> logical unit numbers of various output files (new format/order/grouping)
  integer, save :: i_pflux_es = -1, i_eflux_es = -1, i_vflux_es = -1
  integer, save :: i_pflux_apar = -1, i_eflux_apar = -1, i_vflux_apar = -1
  integer, save :: i_pflux_bpar = -1, i_eflux_bpar = -1, i_vflux_bpar = -1
  integer, save :: i_pflux_es_lab = -1, i_eflux_es_lab = -1, i_vflux_es_lab = -1
  integer, save :: i_pflux_apar_lab = -1, i_eflux_apar_lab = -1, i_vflux_apar_lab = -1
  integer, save :: i_pflux_bpar_lab = -1, i_eflux_bpar_lab = -1, i_vflux_bpar_lab = -1

  !> logical unit numbers for output files (in legacy format/order/grouping)
  integer, save :: i_fluxes_legacy = -1, i_fluxes_em = -1, i_neoclass = -1
  integer, save :: i_fluxes_legacy_lab = -1, i_fluxes_em_lab = -1
  integer, save :: i_neoclass_old = -1
  integer, save :: i_efluxmag = -1, i_pfluxmag = -1
  integer, save :: i_efluxmag_apar = -1, i_pfluxmag_apar = -1
  integer, save :: i_fluxes_bpar = -1, i_fluxes_bpar_lab=-1
  integer, save :: i_efluxspec = -1, i_pfluxspec = -1, i_vfluxspec = -1
  integer, save :: i_efluxxspec = -1, i_pfluxxspec = -1, i_vfluxxspec = -1
  integer, save :: i_pflux_rad_es = -1, i_vflux_rad_es = -1, i_eflux_rad_es = -1 
  integer, save :: i_pflux_rad_apar = -1, i_vflux_rad_apar = -1, i_eflux_rad_apar = -1 
  integer, save :: i_pflux_rad_bpar = -1, i_vflux_rad_bpar = -1, i_eflux_rad_bpar = -1 
  integer, save :: i_efluxspec_apar = -1, i_pfluxspec_apar = -1, i_vfluxspec_apar = -1
  integer, save :: i_efluxxspec_apar = -1, i_pfluxxspec_apar = -1, i_vfluxxspec_apar = -1

  integer, save :: i_deltaprime =-1
  integer, save :: i_torque =-1

  integer, save, allocatable, dimension(:) :: i_eflux_rad_xy, i_eflux_pol_xy

  integer, save, allocatable, dimension(:,:,:,:,:) :: lun_flux_xy
  integer, save, allocatable, dimension(:,:,:) :: lun_flux_k2


  !> number of slots in fluxes write array.
  !> Typically we want all the fluxes to be output.
  integer, parameter :: nfluxes = 3
  
  !> the particle flux per mode : pflux(nmod,nx,number_of_species) due to the 
  !> ExB motion 
  real, save, allocatable, dimension(:,:,:) :: pflux_es

  !> the energy flux per mode : eflux(nmod,nx,number_of_species) due to the 
  !> ExB motion 
  real, save, allocatable, dimension(:,:,:) :: eflux_es
 
  !> Lab frame correction to the energy flux per mode : eflux(nmod,nx,number_of_species) due to the 
  !> ExB motion 
  real, save, allocatable, dimension(:,:,:) :: deflux_es

  !> the toroidal angular momentum flux per mode : vflux(nmod,nx,number_of_species)
  !> due to the ExB motion 
  real, save, allocatable, dimension(:,:,:) :: vflux_es

  !> Lab frame correction to the toroidal angular momentum flux per mode : vflux(nmod,nx,number_of_species)
  !> due to the ExB motion 
  real, save, allocatable, dimension(:,:,:) :: dvflux_es

  !> the particle flux per mode pflux(nmod,nx,number_of_species) due to the
  !> flutter of the field  
  real, save, allocatable, dimension(:,:,:) :: pflux_apar

  !> the energy flux per mode : eflux(nmod,nx,number_of_species) due to the 
  !> flutter of the field 
  real, save, allocatable, dimension(:,:,:) :: eflux_apar

  !> Lab frame correction to the energy flux per mode : eflux(nmod,nx,number_of_species) due to the 
  !> flutter of the field 
  real, save, allocatable, dimension(:,:,:) :: deflux_apar

  !> the toroidal angular momentum flux per mode : vflux(nmod,nx,number_of_species)
  !> due to the flutter of the field  
  real, save, allocatable, dimension(:,:,:) :: vflux_apar

  !> Lab frame correction to the toroidal angular momentum flux per mode : vflux(nmod,nx,number_of_species)
  !> due to the flutter of the field  
  real, save, allocatable, dimension(:,:,:) :: dvflux_apar

  !> arrays for island torque and deltap diagnostics
  real, save, allocatable, dimension(:,:) :: deltap_isl
  real, save, allocatable, dimension(:,:) :: torque_isl
   
  !> the particle flux per mode pflux(nmod,nx,number_of_species) due to the
  !> compression of the field  
  real, save, allocatable, dimension(:,:,:) :: pflux_bpar

  !> the energy flux per mode : eflux(nmod,nx,number_of_species) due to the 
  !> compression of the field 
  real, save, allocatable, dimension(:,:,:) :: eflux_bpar

  !> Lab frame correction to the energy flux per mode : eflux(nmod,nx,number_of_species) due to the 
  !> compression of the field 
  real, save, allocatable, dimension(:,:,:) :: deflux_bpar

  !> the toroidal angular momentum flux per mode : vflux(nmod,nx,number_of_species)
  !> due to the compression of the field  
  real, save, allocatable, dimension(:,:,:) :: vflux_bpar
  
  !> Lab frame correction to the toroidal angular momentum flux per mode : vflux(nmod,nx,number_of_species)
  !> due to the compression of the field  
  real, save, allocatable, dimension(:,:,:) :: dvflux_bpar

  !> total particle flux per species : pflux_tot(number_of_species)
  !> electrostatic, electromagnetic, compression
  real, save, allocatable, dimension(:) :: pflux_tot_es
  real, save, allocatable, dimension(:) :: pflux_tot_apar
  real, save, allocatable, dimension(:) :: pflux_tot_bpar
  
  !> the total energy flux per species : eflux_tot(number_of_species) 
  !> electrostatic, electromagnetic and compression contributions 
  real, save, allocatable, dimension(:) :: eflux_tot_es
  real, save, allocatable, dimension(:) :: eflux_tot_apar
  real, save, allocatable, dimension(:) :: eflux_tot_bpar
  real, save, allocatable, dimension(:) :: eflux_tot_es_lab
  real, save, allocatable, dimension(:) :: eflux_tot_apar_lab
  real, save, allocatable, dimension(:) :: eflux_tot_bpar_lab
  
  !> total toroidal angular momentum flux per species : vflux_tot(number_of_species)
  !> electrostatic, electromagnetic and compression contributions 
  real, save, allocatable, dimension(:) :: vflux_tot_es
  real, save, allocatable, dimension(:) :: vflux_tot_apar
  real, save, allocatable, dimension(:) :: vflux_tot_bpar
  real, save, allocatable, dimension(:) :: vflux_tot_es_lab
  real, save, allocatable, dimension(:) :: vflux_tot_apar_lab
  real, save, allocatable, dimension(:) :: vflux_tot_bpar_lab
  
  !> array for writing lines of the fluxes file
  real, save, allocatable, dimension(:) :: flux_tot_es
  real, save, allocatable, dimension(:) :: flux_tot_apar
  real, save, allocatable, dimension(:) :: flux_tot_bpar
  real, save, allocatable, dimension(:) :: flux_tot_es_lab
  real, save, allocatable, dimension(:) :: flux_tot_apar_lab
  real, save, allocatable, dimension(:) :: flux_tot_bpar_lab
  
  !> The neoclassical particle flux.  pflux_nc(number_of_species)
  real, save, allocatable, dimension(:) :: pflux_nc
  
  !> The neoclassical energy flux. eflux_nc(number_of_species)
  real, save, allocatable, dimension(:) :: eflux_nc

  !> The neoclassical momentum flux. vflux_nc(number_of_species)
  real, save, allocatable, dimension(:) :: vflux_nc
  
  !> The neoclassical particle flux.  pflux_nc_old(number_of_species)
  real, save, allocatable, dimension(:) :: pflux_nc_old
  
  !> The neoclassical energy flux. eflux_nc_old(number_of_species)
  real, save, allocatable, dimension(:) :: eflux_nc_old
  !> The neoclassical potential energy flux. eflux_nc(number_of_species)
  real, save, allocatable, dimension(:) :: eflux2_nc_old

  !> The neoclassical momentum flux. vflux_nc_old(number_of_species)
  real, save, allocatable, dimension(:) :: vflux_nc_old
  
  !> The buffer for communication. fluxbuf(nmod,nx,number_of_species)
  real, save, allocatable, dimension(:,:,:) :: fluxbuf(:,:,:)
  
  !> The buffer for communication. fluxncbuf(number_of_species)
  real, save, allocatable, dimension(:) :: fluxncbuf

  !> arrays for islands diagnostics
  real, save, allocatable, dimension(:) :: deltaprime
  real, save, allocatable, dimension(:) :: torque

  !> spectral particle flux per species : pflux_spec(nmod,number_of_species)
  !> electrostatic and electro-magnetic contributions 
  real, save, allocatable, dimension(:,:) :: pflux_spec
  real, save, allocatable, dimension(:,:) :: pflux_xspec !<(nx,number_of_species)
  real, save, allocatable, dimension(:,:) :: pflux_apar_spec
  real, allocatable, dimension(:,:) :: pflux_apar_xspec !<(nx,number_of_species)

  !> spectral energy flux per species : eflux_spec(nmod,number_of_species)
  !> electrostatic and electromagnetic contributions 
  real, save, allocatable, dimension(:,:) :: eflux_spec
  real, save, allocatable, dimension(:,:) :: eflux_xspec !<(nx,number_of_species)
  real, save, allocatable, dimension(:,:) :: eflux_apar_spec
  real, save, allocatable, dimension(:,:) :: eflux_apar_xspec !<(nx,number_of_species)

  !> the spectral toroidal angular momentum per species : 
  !> vflux_spec(nmod,number_of_species)
  !> electrostatic and electromagnetic contributions 
  real, save, allocatable, dimension(:,:) :: vflux_spec
  real, save, allocatable, dimension(:,:) :: vflux_xspec !<(nx,number_of_species)
  real, save, allocatable, dimension(:,:) :: vflux_apar_spec
  real, save, allocatable, dimension(:,:) :: vflux_apar_xspec !<(nx,number_of_species)

  ! the derivation of the fields with respect to s
  complex, allocatable :: dphi_gads(:), dapar_gads(:), dbpar_gads(:)

  !> Private FFT arrays.
  real, save, allocatable :: rdum(:,:), rdum_yxsp(:,:,:)
  complex, save, allocatable :: a(:,:), b(:,:)
  real, save, allocatable :: ar(:,:), br(:,:), cr(:,:)
  real, save, allocatable :: cr_yxsp(:,:,:)

  !> complex slice in xy for mpi reduction
  complex, save, allocatable, dimension(:,:) :: c_xy
  complex, save, allocatable, dimension(:,:,:) :: c_xysp

  integer, parameter :: flux_history_size = 6
  
  integer, parameter :: RAD_DIRECTION = 1, POL_DIRECTION = 2
  integer, parameter :: FSAVG = 1, NO_FSAVG = 2
  logical, save, allocatable, dimension(:,:,:,:,:) :: flux_xy_enabled

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()

    lfluxes_spectra = .true.
    lfluxes_em_spectra = .false.
    lcalc_fluxes = .true.
    flux3d = .false.

  end subroutine set_default_nml_values


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    
    call mpibcast(lcalc_fluxes,1)
    call mpibcast(lfluxes_spectra,1)
    call mpibcast(lfluxes_em_spectra,1)
    call mpibcast(flux3d,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : nlapar, flux_tube, io_legacy, nlapar, nlbpar
    use grid, only : parallel_x
    use general, only : gkw_warn
    use diagnos_generic, only : xy_fluxes
    use diagnos_generic, only : xy_fluxes_em, xy_fluxes_bpar

    ! output of the EM fluxes spectra only if EM run
    if (.not. nlapar) then
      if (lfluxes_em_spectra) then
        call gkw_warn('Output of the EM fluxes spectra not possible with nlapar=.false.,&
           & force lfluxes_em_spectra=.false.')
        lfluxes_em_spectra = .false.
      end if
      if (xy_fluxes_em) then
        call gkw_warn('Output of the XY EM fluxes not possible with nlapar=.false.,&
           & force xy_fluxes_em=.false.')
        xy_fluxes_em = .false.
      end if
    end if
    if (.not. nlbpar) then
      if (xy_fluxes_bpar) then
        call gkw_warn('Output of the XY EM fluxes not possible with nlbpar=.false.,&
           & force xy_fluxes_em=.false.')
        xy_fluxes_bpar = .false.
      end if
    end if


    if (io_legacy .and. (.not. flux_tube) .and. parallel_x .and. xy_fluxes) then
      call gkw_warn('Output of 2d perp. slices of fluxes in realspace &
         & is incorrect with &
         & io_legacy=T and global runs and therefore disabled! (Please &
         & do not parallelise radially, or use io_legacy=F)')
      ! the reason is that the old code fourier-transforms to higher
      ! resolution, but the B field is not available at that
      ! resolution.
      xy_fluxes = .false.
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use grid, only : number_of_species, n_x_grid, nmod, n_y_grid
    use io, only : open_real_lu, ascii_fmt, binary_fmt
    use io, only : attach_metadata, description_key, comments_key
    use io, only : phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use global, only : int2char_zeros
    use control, only : io_legacy, io_format
    use control, only : nlphi, nlapar, neoclassics, nlbpar
    use mode, only : mode_box
    use mpiinterface, only : root_processor, register_tag_range
    use diagnos_generic, only : mphit, mrad_G
    use diagnos_generic, only : lradial_profile
    use diagnos_generic, only : lphi_diagnostics
    use diagnos_generic, only : xy_fluxes, xy_fluxes_bi, xy_fluxes_bpar
    use diagnos_generic, only : xy_fluxes_em, xy_fluxes_fsa, xy_fluxes_p
    use diagnos_generic, only : xy_fluxes_v, xy_fluxes_k
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD, DISTRIBUTION
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    use diagnos_generic, only : PARTICLE_FLUX, ENERGY_FLUX, MOMENTUM_FLUX
    use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS
    use diagnos_generic, only : VPAR_GHOSTCELLS, MU_GHOSTCELLS
    use components, only : tearingmode, tm_drive
    use global, only : dotdat
    use rho_par_switch, only : lflux_rhostar
    use collisionop, only : conservation
    logical, intent(inout) :: requirements(:,:)

    character(len=25), dimension(number_of_species*nfluxes) :: labels
    integer :: is, iflux
    integer :: fsa
    logical :: eflux_write, pflux_write, vflux_write
    integer :: lun_flux_xy_is
    integer, allocatable :: xyslice_shape(:)

    ! the number 3 appears because there are 3 fluxes (particle, energy, momentum)
    call register_tag_range(3 * number_of_species, &
       & tag_range_start, tag_range_end_inkl)

    requirements(DISTRIBUTION,LOCAL_DATA) = .true.
    if (neoclassics) then
      requirements(DISTRIBUTION,VPAR_GHOSTCELLS) = .true.
      requirements(DISTRIBUTION,MU_GHOSTCELLS) = .true.
    end if

    if(lflux_rhostar) then
      if(nlphi) then
        requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
        requirements(PHI_GA_FIELD,S_GHOSTCELLS) = .true.
      end if
      if(nlapar) then
        requirements(APAR_GA_FIELD,LOCAL_DATA) = .true.
        requirements(APAR_GA_FIELD,S_GHOSTCELLS) = .true.
      end if
      if(nlbpar) then
        requirements(BPAR_GA_FIELD,LOCAL_DATA) = .true.
        requirements(BPAR_GA_FIELD,S_GHOSTCELLS) = .true.
      end if
    end if

    ! default values for flags:
    flux_xy_enabled = .false.

    if(root_processor) then
      ! open the file to write the fluxes
      if(io_legacy) then
        do is = 1, number_of_species
          iflux=1
          labels(iflux+(is-1)*nfluxes) = 'pflux_species'// &
             & trim(adjustl(int2char_zeros(is,2)))
          iflux=iflux+1
          labels(iflux+(is-1)*nfluxes) = 'eflux_species'// &
             & trim(adjustl(int2char_zeros(is,2)))
          iflux=iflux+1         
          labels(iflux+(is-1)*nfluxes) = 'vflux_species'// &
             & trim(adjustl(int2char_zeros(is,2)))
        end do

        if (nlphi) then
          call open_real_lu(dotdat('fluxes',io_legacy), &
             & 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species*nfluxes /), &
             & ascii_fmt, i_fluxes_legacy, labels)
          call open_real_lu(dotdat('fluxes_lab',io_legacy), &
             & 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species*nfluxes /), &
             & ascii_fmt, i_fluxes_legacy_lab, labels)
        end if

        if (nlapar) then
          call open_real_lu(dotdat('fluxes_em', io_legacy), &
             & 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species*nfluxes /), &
             & ascii_fmt, i_fluxes_em, labels)
          call open_real_lu(dotdat('fluxes_em_lab', io_legacy), &
             & 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species*nfluxes /), &
             & ascii_fmt, i_fluxes_em_lab, labels)
        end if

        if (nlbpar) then
          call open_real_lu(dotdat('fluxes_bpar', io_legacy), &
             & 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species*nfluxes /), &
             & ascii_fmt, i_fluxes_bpar, labels)
          call open_real_lu(dotdat('fluxes_bpar_lab', io_legacy), &
             & 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species*nfluxes /), &
             & ascii_fmt, i_fluxes_bpar_lab, labels)
        end if

      else

        if (nlphi) then
          call open_real_lu('pflux_es', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_pflux_es)
          call attach_metadata_grid(i_pflux_es, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_pflux_es, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_es, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_es, comments_key, not_avail, ascii_fmt)

          call open_real_lu('eflux_es', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_eflux_es)
          call attach_metadata_grid(i_eflux_es, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_eflux_es, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_es, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_es, comments_key, not_avail, ascii_fmt)

          call open_real_lu('vflux_es', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_vflux_es)
          call attach_metadata_grid(i_vflux_es, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_vflux_es, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_es, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_es, comments_key, not_avail, ascii_fmt)

          call open_real_lu('pflux_es_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_pflux_es_lab)
          call attach_metadata_grid(i_pflux_es_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_pflux_es_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_es_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_es_lab, comments_key, not_avail, ascii_fmt)

          call open_real_lu('eflux_es_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_eflux_es_lab)
          call attach_metadata_grid(i_eflux_es_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_eflux_es_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_es_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_es_lab, comments_key, not_avail, ascii_fmt)

          call open_real_lu('vflux_es_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_vflux_es_lab)
          call attach_metadata_grid(i_vflux_es_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_vflux_es_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_es_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_es_lab, comments_key, not_avail, ascii_fmt)
        end if

        if(nlapar) then
          call open_real_lu('pflux_apar', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_pflux_apar)
          call attach_metadata_grid(i_pflux_apar, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_pflux_apar, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_apar, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_apar, comments_key, not_avail, ascii_fmt)

          call open_real_lu('eflux_apar', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_eflux_apar)
          call attach_metadata_grid(i_eflux_apar, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_eflux_apar, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_apar, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_apar, comments_key, not_avail, ascii_fmt)

          call open_real_lu('vflux_apar', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_vflux_apar)
          call attach_metadata_grid(i_vflux_apar, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_vflux_apar, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_apar, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_apar, comments_key, not_avail, ascii_fmt) 

          call open_real_lu('pflux_apar_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_pflux_apar_lab)
          call attach_metadata_grid(i_pflux_apar_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_pflux_apar_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_apar_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_apar_lab, comments_key, not_avail, ascii_fmt)

          call open_real_lu('eflux_apar_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_eflux_apar_lab)
          call attach_metadata_grid(i_eflux_apar_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_eflux_apar_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_apar_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_apar_lab, comments_key, not_avail, ascii_fmt)

          call open_real_lu('vflux_apar_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_vflux_apar_lab)
          call attach_metadata_grid(i_vflux_apar_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_vflux_apar_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_apar_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_apar_lab, comments_key, not_avail, ascii_fmt) 
        end if

        if(nlbpar) then
          call open_real_lu('pflux_bpar', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_pflux_bpar)
          call attach_metadata_grid(i_pflux_bpar, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_pflux_bpar, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_bpar, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_bpar, comments_key, not_avail, ascii_fmt)

          call open_real_lu('eflux_bpar', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_eflux_bpar)
          call attach_metadata_grid(i_eflux_bpar, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_eflux_bpar, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_bpar, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_bpar, comments_key, not_avail, ascii_fmt)

          call open_real_lu('vflux_bpar', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_vflux_bpar)
          call attach_metadata_grid(i_vflux_bpar, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_vflux_bpar, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_bpar, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_bpar, comments_key, not_avail, ascii_fmt) 

          call open_real_lu('pflux_bpar_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_pflux_bpar_lab)
          call attach_metadata_grid(i_pflux_bpar_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_pflux_bpar_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_bpar_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pflux_bpar_lab, comments_key, not_avail, ascii_fmt)

          call open_real_lu('eflux_bpar_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_eflux_bpar_lab)
          call attach_metadata_grid(i_eflux_bpar_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_eflux_bpar_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_bpar_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_eflux_bpar_lab, comments_key, not_avail, ascii_fmt)

          call open_real_lu('vflux_bpar_lab', 'diagnostic/diagnos_fluxes', &
             & (/ number_of_species /), ascii_fmt, i_vflux_bpar_lab)
          call attach_metadata_grid(i_vflux_bpar_lab, 'time', 'number_of_species', ascii_fmt)
          call attach_metadata(i_vflux_bpar_lab, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_bpar_lab, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vflux_bpar_lab, comments_key, not_avail, ascii_fmt) 
        end if

      end if

      if (neoclassics) then
        call open_real_lu(dotdat('fluxes_nc_old', io_legacy), &
           & 'diagnostic/diagnos_fluxes', (/ 4*number_of_species /), &
           & ascii_fmt, i_neoclass_old)
        call attach_metadata_grid(i_neoclass_old, 'time', '3*number_of_species', ascii_fmt)
        call attach_metadata(i_neoclass_old, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_neoclass_old, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_neoclass_old, comments_key, not_avail, ascii_fmt)

        if(conservation) then
          call open_real_lu(dotdat('fluxes_nc', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ 3*number_of_species /), &
             & ascii_fmt, i_neoclass)
          call attach_metadata_grid(i_neoclass, 'time', '3*number_of_species', ascii_fmt)
          call attach_metadata(i_neoclass, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_neoclass, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_neoclass, comments_key, not_avail, ascii_fmt)
        end if
      end if


      ! radial profiles of fluxes
      if (lradial_profile) then       
        ! electrostatic flux due to ExB-drift
        if(nlphi) then
          if(io_format == 'hdf5') then
            call open_real_lu('pflux_rad_es', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_pflux_rad_es)
            call attach_metadata_grid(i_pflux_rad_es, 'time', 'n_x_grid' , 'number_of_species', ascii_fmt)
            call attach_metadata(i_pflux_rad_es, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_es, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_es, comments_key, not_avail, ascii_fmt)

            call open_real_lu('vflux_rad_es', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_vflux_rad_es)
            call attach_metadata_grid(i_vflux_rad_es, 'time', 'n_x_grid' , 'number_of_species', ascii_fmt)
            call attach_metadata(i_vflux_rad_es, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_es, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_es, comments_key, not_avail, ascii_fmt)

            call open_real_lu('eflux_rad_es', 'diagnostic/diagnos_fluxes', &
               &(/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_eflux_rad_es)
            call attach_metadata_grid(i_eflux_rad_es, 'time', 'n_x_grid' , 'number_of_species', ascii_fmt)
            call attach_metadata(i_eflux_rad_es, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_es, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_es, comments_key, not_avail, ascii_fmt)
          else
            call open_real_lu('pflux_rad_es', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_pflux_rad_es)
            call attach_metadata_grid(i_pflux_rad_es, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_pflux_rad_es, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_es, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_es, comments_key, not_avail, ascii_fmt)

            call open_real_lu('vflux_rad_es', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_vflux_rad_es)
            call attach_metadata_grid(i_vflux_rad_es, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_vflux_rad_es, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_es, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_es, comments_key, not_avail, ascii_fmt)

            call open_real_lu('eflux_rad_es', 'diagnostic/diagnos_fluxes', &
               &(/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_eflux_rad_es)
            call attach_metadata_grid(i_eflux_rad_es, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_eflux_rad_es, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_es, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_es, comments_key, not_avail, ascii_fmt)
          end if
        end if
        
        
        ! electromagnetic flux due to magnetic flutter 
        if(nlapar) then
          if(io_format == 'hdf5') then
            call open_real_lu('pflux_rad_apar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_pflux_rad_apar)
            call attach_metadata_grid(i_pflux_rad_apar, 'time', 'n_x_grid', 'number_of_species', ascii_fmt)
            call attach_metadata(i_pflux_rad_apar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_apar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_apar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('vflux_rad_apar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_vflux_rad_apar)
            call attach_metadata_grid(i_vflux_rad_apar, 'time', 'n_x_grid', 'number_of_species', ascii_fmt)
            call attach_metadata(i_vflux_rad_apar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_apar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_apar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('eflux_rad_apar', 'diagnostic/diagnos_fluxes', &
               &(/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_eflux_rad_apar)
            call attach_metadata_grid(i_eflux_rad_apar, 'time', 'n_x_grid', 'number_of_species', ascii_fmt)
            call attach_metadata(i_eflux_rad_apar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_apar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_apar, comments_key, not_avail, ascii_fmt)
          else
            call open_real_lu('pflux_rad_apar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_pflux_rad_apar)
            call attach_metadata_grid(i_pflux_rad_apar, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_pflux_rad_apar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_apar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_apar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('vflux_rad_apar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_vflux_rad_apar)
            call attach_metadata_grid(i_vflux_rad_apar, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_vflux_rad_apar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_apar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_apar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('eflux_rad_apar', 'diagnostic/diagnos_fluxes', &
               &(/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_eflux_rad_apar)
            call attach_metadata_grid(i_eflux_rad_apar, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_eflux_rad_apar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_apar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_apar, comments_key, not_avail, ascii_fmt)
          end if
        end if
        
        
        ! electromagnetic flux due to magnetic compression
        if(nlbpar) then
          if(io_format == 'hdf5') then
            call open_real_lu('pflux_rad_bpar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_pflux_rad_bpar)
            call attach_metadata_grid(i_pflux_rad_bpar, 'time', 'n_x_grid', 'number_of_species', ascii_fmt)
            call attach_metadata(i_pflux_rad_bpar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_bpar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_bpar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('vflux_rad_bpar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_vflux_rad_bpar)
            call attach_metadata_grid(i_vflux_rad_bpar, 'time', 'n_x_grid', 'number_of_species', ascii_fmt)
            call attach_metadata(i_vflux_rad_bpar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_bpar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_bpar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('eflux_rad_bpar', 'diagnostic/diagnos_fluxes', &
               &(/ n_x_grid , number_of_species /), &
               & ascii_fmt, i_eflux_rad_bpar)
            call attach_metadata_grid(i_eflux_rad_bpar, 'time', 'n_x_grid', 'number_of_species', ascii_fmt)
            call attach_metadata(i_eflux_rad_bpar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_bpar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_bpar, comments_key, not_avail, ascii_fmt)
          else
            call open_real_lu('pflux_rad_bpar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_pflux_rad_bpar)
            call attach_metadata_grid(i_pflux_rad_bpar, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_pflux_rad_bpar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_bpar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_pflux_rad_bpar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('vflux_rad_bpar', 'diagnostic/diagnos_fluxes', &
               & (/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_vflux_rad_bpar)
            call attach_metadata_grid(i_vflux_rad_bpar, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_vflux_rad_bpar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_bpar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_vflux_rad_bpar, comments_key, not_avail, ascii_fmt)

            call open_real_lu('eflux_rad_bpar', 'diagnostic/diagnos_fluxes', &
               &(/ n_x_grid*number_of_species /), &
               & ascii_fmt, i_eflux_rad_bpar)
            call attach_metadata_grid(i_eflux_rad_bpar, 'time', 'n_x_grid*number_of_species', ascii_fmt)
            call attach_metadata(i_eflux_rad_bpar, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_bpar, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_bpar, comments_key, not_avail, ascii_fmt)
          end if
        end if
      endif

      if (lfluxes_spectra) then
        call open_real_lu(dotdat('eflux_sup', io_legacy), &
           & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species /), &
           & ascii_fmt, i_efluxmag)
        call attach_metadata_grid(i_efluxmag, 'time', 'nmod*number_of_species', ascii_fmt)
        call attach_metadata(i_efluxmag, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_efluxmag, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_efluxmag, comments_key, not_avail, ascii_fmt)

        call open_real_lu(dotdat('pflux_sup', io_legacy), &
           & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species /), &
           & ascii_fmt, i_pfluxmag)
        call attach_metadata_grid(i_pfluxmag, 'time', 'nmod*number_of_species', ascii_fmt)
        call attach_metadata(i_pfluxmag, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_pfluxmag, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_pfluxmag, comments_key, not_avail, ascii_fmt)

        if (nmod > 1) then
          call open_real_lu(dotdat('eflux_spectra', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species /), &
             & ascii_fmt, i_efluxspec)
          call attach_metadata_grid(i_efluxspec, 'time', 'nmod*number_of_species', ascii_fmt)
          call attach_metadata(i_efluxspec, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_efluxspec, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_efluxspec, comments_key, not_avail, ascii_fmt)

          call open_real_lu(dotdat('pflux_spectra', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species /), &
             & ascii_fmt, i_pfluxspec)
          call attach_metadata_grid(i_pfluxspec, 'time', 'nmod*number_of_species', ascii_fmt)
          call attach_metadata(i_pfluxspec, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pfluxspec, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pfluxspec, comments_key, not_avail, ascii_fmt)

          call open_real_lu(dotdat('vflux_spectra', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species /), &
             & ascii_fmt, i_vfluxspec)
          call attach_metadata_grid(i_vfluxspec, 'time', 'nmod*number_of_species', ascii_fmt)
          call attach_metadata(i_vfluxspec, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vfluxspec, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vfluxspec, comments_key, not_avail, ascii_fmt)

        end if
        if (n_x_grid > 1) then
          call open_real_lu(dotdat('eflux_xspec', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ n_x_grid*number_of_species/), &
             & ascii_fmt, i_efluxxspec)
          call attach_metadata_grid(i_efluxxspec, 'time', 'n_x_grid*number_of_species', ascii_fmt)
          call attach_metadata(i_efluxxspec, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_efluxxspec, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_efluxxspec, comments_key, not_avail, ascii_fmt)

          call open_real_lu(dotdat('pflux_xspec', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ n_x_grid*number_of_species/), &
             & ascii_fmt, i_pfluxxspec)
          call attach_metadata_grid(i_pfluxxspec, 'time', 'n_x_grid*number_of_species', ascii_fmt)
          call attach_metadata(i_pfluxxspec, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_pfluxxspec, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_pfluxxspec, comments_key, not_avail, ascii_fmt)

          call open_real_lu(dotdat('vflux_xspec', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ n_x_grid*number_of_species/), &
             & ascii_fmt, i_vfluxxspec)
          call attach_metadata_grid(i_vfluxxspec, 'time', 'n_x_grid*number_of_species', ascii_fmt)
          call attach_metadata(i_vfluxxspec, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_vfluxxspec, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_vfluxxspec, comments_key, not_avail, ascii_fmt)
        end if
      end if

      if (lfluxes_em_spectra) then
        call open_real_lu(dotdat('eflux_em_sup', io_legacy), &
           & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species /), &
           & ascii_fmt, i_efluxmag_apar)

        call open_real_lu(dotdat('pflux_em_sup', io_legacy), &
           & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species /), &
           & ascii_fmt, i_pfluxmag_apar)
        if (nmod > 1) then
          call open_real_lu(dotdat('eflux_em_spectra', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species/), &
             & ascii_fmt, i_efluxspec_apar)
          call open_real_lu(dotdat('pflux_em_spectra', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species/), &
             & ascii_fmt, i_pfluxspec_apar)
          call open_real_lu(dotdat('vflux_em_spectra', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ nmod*number_of_species/), &
             & ascii_fmt, i_vfluxspec_apar)
        end if
        if (n_x_grid > 1) then
          call open_real_lu(dotdat('eflux_em_xspec', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ n_x_grid*number_of_species /), &
             & ascii_fmt, i_efluxxspec_apar)
          call open_real_lu(dotdat('pflux_em_xspec', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ n_x_grid*number_of_species /), &
             & ascii_fmt, i_pfluxxspec_apar)
          call open_real_lu(dotdat('vflux_em_xspec', io_legacy), &
             & 'diagnostic/diagnos_fluxes', (/ n_x_grid*number_of_species /), &
             & ascii_fmt, i_vfluxxspec_apar)
        end if
      end if


      if((tearingmode.or.tm_drive).and. nlapar) then
        call open_real_lu(dotdat('deltaprime', io_legacy), &
           & 'diagnostic/diagnos_fluxes', (/ number_of_species /), &
           & ascii_fmt, i_deltaprime)
        call attach_metadata_grid(i_deltaprime, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_deltaprime, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_deltaprime, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_deltaprime, comments_key, not_avail, ascii_fmt)

        call open_real_lu(dotdat('torque', io_legacy), &
           & 'diagnostic/diagnos_fluxes', (/ number_of_species /), &
           & ascii_fmt, i_torque)
        call attach_metadata_grid(i_torque, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_torque, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_torque, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_torque, comments_key, not_avail, ascii_fmt)
      end if
      
      if (mode_box .and. lphi_diagnostics) then
        if(xy_fluxes) then
          ! flux magnitude radial
          if(io_legacy) then
            do is = 1, number_of_species
              call open_real_lu(trim('flmgr')//trim(int2char_zeros(is,2)), &
                 & 'diagnostic/diagnos_fluxes', (/ mphit, mrad_G /), &
                 & binary_fmt, i_eflux_rad_xy(is))
            end do
          else
            call open_real_lu(trim('flmgr'), &
               & 'diagnostic/diagnos_fluxes', &
               & (/ n_y_grid, n_x_grid, number_of_species /), &
               & binary_fmt, i_eflux_rad_xy(1))
            call attach_metadata_grid(i_eflux_rad_xy(1), 'time', 'n_y_grid', &
               & 'n_x_grid', 'number_of_species', ascii_fmt)
            call attach_metadata(i_eflux_rad_xy(1), phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_xy(1), description_key, not_avail, ascii_fmt)
            call attach_metadata(i_eflux_rad_xy(1), comments_key, not_avail, ascii_fmt)
          end if

          if(xy_fluxes_bi) then
            ! flux magnitude poloidal
            if(io_legacy) then
              do is = 1, number_of_species
                call open_real_lu(trim('flmgp')//trim(int2char_zeros(is,2)), &
                   & 'diagnostic/diagnos_fluxes', (/ mphit, mrad_G /), &
                   & binary_fmt, i_eflux_pol_xy(is))
              end do
            else
              call open_real_lu(trim('flmgp'), &
                 & 'diagnostic/diagnos_fluxes', &
                 & (/ n_y_grid, n_x_grid, number_of_species /), &
                 & binary_fmt, i_eflux_pol_xy(1))
              call attach_metadata_grid(i_eflux_pol_xy(1), 'time', 'n_y_grid', &
                 & 'n_x_grid', 'number_of_species', ascii_fmt)
              call attach_metadata(i_eflux_pol_xy(1), phys_unit_key, not_avail, ascii_fmt)
              call attach_metadata(i_eflux_pol_xy(1), description_key, not_avail, ascii_fmt)
              call attach_metadata(i_eflux_pol_xy(1), comments_key, not_avail, ascii_fmt)
            end if
          end if
        end if
      end if

    end if

    eflux_write = xy_fluxes
    pflux_write = xy_fluxes .and. xy_fluxes_p
    vflux_write = xy_fluxes .and. xy_fluxes_v

    if (mode_box .and. lphi_diagnostics) then
      do is = 1, number_of_species

        if((.not. io_legacy) .and. is > 1) exit
        ! because then all the data goes into one lu, instead of one
        ! lu per species.
        
        ! loop over two types of output
        fsavg_or_nofsavg: do fsa = FSAVG, NO_FSAVG
          select case(fsa)
          case(FSAVG)
            ! slice will be flux surface averaged
            if(.not.xy_fluxes_fsa) then
              cycle fsavg_or_nofsavg
            end if
          end select
          if(io_legacy) then
            if(.not. allocated(xyslice_shape)) allocate(xyslice_shape(2))
            xyslice_shape(1) = mphit
            xyslice_shape(2) = mrad_G
            ! for each species, a separate lu is opened, species
            ! index goes into the lu name
            lun_flux_xy_is = is
          else
            if(.not. allocated(xyslice_shape)) allocate(xyslice_shape(3))
            xyslice_shape(1) = n_y_grid
            xyslice_shape(2) = n_x_grid
            xyslice_shape(3) = number_of_species
            ! data of all species goes into one lu, species index
            ! becomes a dimension
            lun_flux_xy_is = 1
          end if

          if (eflux_write) then
            flux_xy_enabled(RAD_DIRECTION,PHI_FIELD,ENERGY_FLUX,fsa,is) = .true.
            if(root_processor) then
              call open_real_lu( &
                 & get_xy_flux_luname(RAD_DIRECTION,PHI_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), &
                 & 'diagnostic/diagnos_fluxes', &
                 & xyslice_shape, &
                 & binary_fmt, &
                 & lun_flux_xy(RAD_DIRECTION,PHI_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is))
              if(.not.io_legacy) then
                call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                   & PHI_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), 'time', &
                   & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & ENERGY_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & ENERGY_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & ENERGY_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
              end if
            end if

            if (XY_fluxes_em) then
              flux_xy_enabled(RAD_DIRECTION,APAR_FIELD,ENERGY_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(RAD_DIRECTION,APAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(RAD_DIRECTION,APAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                     & APAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if
            end if
            if (XY_fluxes_bpar) then
              flux_xy_enabled(RAD_DIRECTION,BPAR_FIELD,ENERGY_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(RAD_DIRECTION,BPAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                     & BPAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if
            end if

            if (XY_fluxes_bi) then
              flux_xy_enabled(POL_DIRECTION,PHI_FIELD,ENERGY_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(POL_DIRECTION,PHI_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(POL_DIRECTION,PHI_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                     & PHI_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & ENERGY_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if

              if (XY_fluxes_em) then
                flux_xy_enabled(POL_DIRECTION,APAR_FIELD,ENERGY_FLUX,fsa,is) = .true.
                if(root_processor) then
                  call open_real_lu( &
                     & get_xy_flux_luname(POL_DIRECTION,APAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), &
                     & 'diagnostic/diagnos_fluxes', &
                     & xyslice_shape, &
                     & binary_fmt, &
                     & lun_flux_xy(POL_DIRECTION,APAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is))
                  if(.not.io_legacy) then
                    call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                       & APAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), 'time', &
                       & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & ENERGY_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & ENERGY_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & ENERGY_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                  end if
                end if
              end if
              if (XY_fluxes_bpar) then    
                flux_xy_enabled(POL_DIRECTION,BPAR_FIELD,ENERGY_FLUX,fsa,is) = .true.
                if(root_processor) then
                  call open_real_lu( &
                     & get_xy_flux_luname(POL_DIRECTION,BPAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), &
                     & 'diagnostic/diagnos_fluxes', &
                     & xyslice_shape, &
                     & binary_fmt, &
                     & lun_flux_xy(POL_DIRECTION,BPAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is))
                  if(.not.io_legacy) then
                    call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                       & BPAR_FIELD,ENERGY_FLUX,fsa,lun_flux_xy_is), 'time', &
                       & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & ENERGY_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & ENERGY_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & ENERGY_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                  end if
                end if
              end if
            end if
          end if

          if (pflux_write) then
            flux_xy_enabled(RAD_DIRECTION,PHI_FIELD,PARTICLE_FLUX,fsa,is) = .true.
            if(root_processor) then
              call open_real_lu( &
                 & get_xy_flux_luname(RAD_DIRECTION,PHI_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), &
                 & 'diagnostic/diagnos_fluxes', &
                 & xyslice_shape, &
                 & binary_fmt, &
                 & lun_flux_xy(RAD_DIRECTION,PHI_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is))
              if(.not.io_legacy) then
                call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                   & PHI_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), 'time', &
                   & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & PARTICLE_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & PARTICLE_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & PARTICLE_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
              end if
            end if

            if (XY_fluxes_em) then
              flux_xy_enabled(RAD_DIRECTION,APAR_FIELD,PARTICLE_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(RAD_DIRECTION,APAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(RAD_DIRECTION,APAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                     & APAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if
            end if
            if (XY_fluxes_bpar) then
              flux_xy_enabled(RAD_DIRECTION,BPAR_FIELD,PARTICLE_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(RAD_DIRECTION,BPAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                     & BPAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if
            end if
            if (XY_fluxes_bi) then
              flux_xy_enabled(POL_DIRECTION,PHI_FIELD,PARTICLE_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(POL_DIRECTION,PHI_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(POL_DIRECTION,PHI_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                     & PHI_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & PARTICLE_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if

              if (XY_fluxes_em) then
                flux_xy_enabled(POL_DIRECTION,APAR_FIELD,PARTICLE_FLUX,fsa,is) = .true.
                if(root_processor) then
                  call open_real_lu( &
                     & get_xy_flux_luname(POL_DIRECTION,APAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), &
                     & 'diagnostic/diagnos_fluxes', &
                     & xyslice_shape, &
                     & binary_fmt, &
                     & lun_flux_xy(POL_DIRECTION,APAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is))
                  if(.not.io_legacy) then
                    call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                       & APAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), 'time', &
                       & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & PARTICLE_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & PARTICLE_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & PARTICLE_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                  end if
                end if
              end if
              if (XY_fluxes_bpar) then
                flux_xy_enabled(POL_DIRECTION,BPAR_FIELD,PARTICLE_FLUX,fsa,is) = .true.
                if(root_processor) then
                  call open_real_lu( &
                     & get_xy_flux_luname(POL_DIRECTION,BPAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), &
                     & 'diagnostic/diagnos_fluxes', &
                     & xyslice_shape, &
                     & binary_fmt, &
                     & lun_flux_xy(POL_DIRECTION,BPAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is))
                  if(.not.io_legacy) then
                    call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                       & BPAR_FIELD,PARTICLE_FLUX,fsa,lun_flux_xy_is), 'time', &
                       & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & PARTICLE_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & PARTICLE_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & PARTICLE_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                  end if
                end if
              end if
            end if
          end if

          if (vflux_write) then
            flux_xy_enabled(RAD_DIRECTION,PHI_FIELD,MOMENTUM_FLUX,fsa,is) = .true.
            if(root_processor) then
              call open_real_lu( &
                 & get_xy_flux_luname(RAD_DIRECTION,PHI_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), &
                 & 'diagnostic/diagnos_fluxes', &
                 & xyslice_shape, &
                 & binary_fmt, &
                 & lun_flux_xy(RAD_DIRECTION,PHI_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is))
              if(.not.io_legacy) then
                call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                   & PHI_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), 'time', &
                   & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & MOMENTUM_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & MOMENTUM_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                call attach_metadata(lun_flux_xy(RAD_DIRECTION,PHI_FIELD,&
                   & MOMENTUM_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
              end if
            end if

            if (XY_fluxes_em) then
              flux_xy_enabled(RAD_DIRECTION,APAR_FIELD,MOMENTUM_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(RAD_DIRECTION,APAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(RAD_DIRECTION,APAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                     & APAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,APAR_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if
            end if
            if (XY_fluxes_bpar) then
              flux_xy_enabled(RAD_DIRECTION,BPAR_FIELD,MOMENTUM_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(RAD_DIRECTION,BPAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(RAD_DIRECTION,&
                     & BPAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(RAD_DIRECTION,BPAR_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if
            end if

            if (XY_fluxes_bi) then
              flux_xy_enabled(POL_DIRECTION,PHI_FIELD,MOMENTUM_FLUX,fsa,is) = .true.
              if(root_processor) then
                call open_real_lu( &
                   & get_xy_flux_luname(POL_DIRECTION,PHI_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), &
                   & 'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, &
                   & binary_fmt, &
                   & lun_flux_xy(POL_DIRECTION,PHI_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                     & PHI_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), 'time', &
                     & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                  call attach_metadata(lun_flux_xy(POL_DIRECTION,PHI_FIELD,&
                     & MOMENTUM_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                end if
              end if

              if (XY_fluxes_em) then
                flux_xy_enabled(POL_DIRECTION,APAR_FIELD,MOMENTUM_FLUX,fsa,is) = .true.
                if(root_processor) then
                  call open_real_lu( &
                     & get_xy_flux_luname(POL_DIRECTION,APAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), &
                     & 'diagnostic/diagnos_fluxes', &
                     & xyslice_shape, &
                     & binary_fmt, &
                     & lun_flux_xy(POL_DIRECTION,APAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is))
                  if(.not.io_legacy) then
                    call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                       & APAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), 'time', &
                       & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & MOMENTUM_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & MOMENTUM_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,APAR_FIELD,&
                       & MOMENTUM_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                  end if
                end if
              end if
              if (XY_fluxes_bpar) then
                flux_xy_enabled(POL_DIRECTION,BPAR_FIELD,MOMENTUM_FLUX,fsa,is) = .true.
                if(root_processor) then
                  call open_real_lu( &
                     & get_xy_flux_luname(POL_DIRECTION,BPAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), &
                     & 'diagnostic/diagnos_fluxes', &
                     & xyslice_shape, &
                     & binary_fmt, &
                     & lun_flux_xy(POL_DIRECTION,BPAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is))
                  if(.not.io_legacy) then
                    call attach_metadata_grid(lun_flux_xy(POL_DIRECTION,&
                       & BPAR_FIELD,MOMENTUM_FLUX,fsa,lun_flux_xy_is), 'time', &
                       & 'n_y_grid', 'n_x_grid', 'number_of_species', ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & MOMENTUM_FLUX,fsa,lun_flux_xy_is), phys_unit_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & MOMENTUM_FLUX,fsa,lun_flux_xy_is), description_key, not_avail, ascii_fmt)
                    call attach_metadata(lun_flux_xy(POL_DIRECTION,BPAR_FIELD,&
                       & MOMENTUM_FLUX,fsa,lun_flux_xy_is), comments_key, not_avail, ascii_fmt)
                  end if
                end if
              end if
            end if
          end if
        end do fsavg_or_nofsavg
      end do

      if(root_processor) then
        ! Open lus for the 2d spectral fields, too:
        
        if(io_legacy) then
          if(.not. allocated(xyslice_shape)) allocate(xyslice_shape(2))
          xyslice_shape(1) = nmod
          xyslice_shape(2) = n_x_grid
          ! for each species, a separate lu is opened, species
          ! index goes into the lu name
        else
          if(.not. allocated(xyslice_shape)) allocate(xyslice_shape(3))
          xyslice_shape(1) = nmod
          xyslice_shape(2) = n_x_grid
          xyslice_shape(3) = number_of_species
          ! data of all species goes into one lu, species index
          ! becomes a dimension
        end if
        
        do is = 1, number_of_species

          if(.not. io_legacy .and. is > 1) exit
          ! because then all the data goes into one lu, instead of one
          ! lu per species.
          
          if (xy_fluxes_k) then
            call open_real_lu(get_k_flux_luname(PHI_FIELD, ENERGY_FLUX, is), &
               & 'diagnostic/diagnos_fluxes', &
               & xyslice_shape, binary_fmt, &
               & lun_flux_k2(PHI_FIELD, ENERGY_FLUX, is))
            if(.not.io_legacy) then
              call attach_metadata_grid(lun_flux_k2(PHI_FIELD,&
                 & ENERGY_FLUX,is), 'time', &
                 & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
              call attach_metadata(lun_flux_k2(PHI_FIELD,&
                 & ENERGY_FLUX,is), phys_unit_key, not_avail, binary_fmt)
              call attach_metadata(lun_flux_k2(PHI_FIELD,&
                 & ENERGY_FLUX,is), description_key, not_avail, binary_fmt)
              call attach_metadata(lun_flux_k2(PHI_FIELD,&
                 & ENERGY_FLUX,is), comments_key, not_avail, binary_fmt)
            end if
            
            if (xy_fluxes_P) then
              call open_real_lu(get_k_flux_luname(PHI_FIELD, PARTICLE_FLUX, is), &
                 & 'diagnostic/diagnos_fluxes', &
                 & xyslice_shape, binary_fmt, &
                 & lun_flux_k2(PHI_FIELD, PARTICLE_FLUX, is))
              if(.not.io_legacy) then
                call attach_metadata_grid(lun_flux_k2(PHI_FIELD,&
                   & PARTICLE_FLUX,is), 'time', &
                   & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
                call attach_metadata(lun_flux_k2(PHI_FIELD,&
                   & PARTICLE_FLUX,is), phys_unit_key, not_avail, binary_fmt)
                call attach_metadata(lun_flux_k2(PHI_FIELD,&
                   & PARTICLE_FLUX,is), description_key, not_avail, binary_fmt)
                call attach_metadata(lun_flux_k2(PHI_FIELD,&
                   & PARTICLE_FLUX,is), comments_key, not_avail, binary_fmt)
              end if
            end if

          if (xy_fluxes_v) then
            call open_real_lu(get_k_flux_luname(PHI_FIELD, MOMENTUM_FLUX, is), &
               & 'diagnostic/diagnos_fluxes', &
               & xyslice_shape, binary_fmt, &
               & lun_flux_k2(PHI_FIELD, MOMENTUM_FLUX, is))
            if(.not.io_legacy) then
              call attach_metadata_grid(lun_flux_k2(PHI_FIELD,&
                 & MOMENTUM_FLUX,is), 'time', &
                 & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
              call attach_metadata(lun_flux_k2(PHI_FIELD,&
                 & MOMENTUM_FLUX,is), phys_unit_key, not_avail, binary_fmt)
              call attach_metadata(lun_flux_k2(PHI_FIELD,&
                 & MOMENTUM_FLUX,is), description_key, not_avail, binary_fmt)
              call attach_metadata(lun_flux_k2(PHI_FIELD,&
                 & MOMENTUM_FLUX,is), comments_key, not_avail, binary_fmt)
            end if
          end if

          if (xy_fluxes_em) then
            call open_real_lu(get_k_flux_luname(APAR_FIELD, ENERGY_FLUX, is), &
               & 'diagnostic/diagnos_fluxes', &
               & xyslice_shape, binary_fmt, &
               & lun_flux_k2(APAR_FIELD, ENERGY_FLUX, is))
            if(.not.io_legacy) then
              call attach_metadata_grid(lun_flux_k2(APAR_FIELD,&
                 & ENERGY_FLUX,is), 'time', &
                 & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
              call attach_metadata(lun_flux_k2(APAR_FIELD,&
                 & ENERGY_FLUX,is), phys_unit_key, not_avail, binary_fmt)
              call attach_metadata(lun_flux_k2(APAR_FIELD,&
                 & ENERGY_FLUX,is), description_key, not_avail, binary_fmt)
              call attach_metadata(lun_flux_k2(APAR_FIELD,&
                 & ENERGY_FLUX,is), comments_key, not_avail, binary_fmt)
            end if
            
              if (xy_fluxes_P) then
                call open_real_lu(get_k_flux_luname(APAR_FIELD, PARTICLE_FLUX, &
                   & is),'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, binary_fmt, &
                   & lun_flux_k2(APAR_FIELD, PARTICLE_FLUX, is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_k2(APAR_FIELD,&
                     & PARTICLE_FLUX,is), 'time', &
                     & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
                  call attach_metadata(lun_flux_k2(APAR_FIELD,&
                     & PARTICLE_FLUX,is), phys_unit_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(APAR_FIELD,&
                     & PARTICLE_FLUX,is), description_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(APAR_FIELD,&
                     & PARTICLE_FLUX,is), comments_key, not_avail, binary_fmt)
                end if
              end if
              if (xy_fluxes_v) then
                call open_real_lu(get_k_flux_luname(APAR_FIELD, MOMENTUM_FLUX, &
                   & is),'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, binary_fmt, &
                   & lun_flux_k2(APAR_FIELD, MOMENTUM_FLUX, is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_k2(APAR_FIELD,&
                     & MOMENTUM_FLUX,is), 'time', &
                     & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
                  call attach_metadata(lun_flux_k2(APAR_FIELD,&
                     & MOMENTUM_FLUX,is), phys_unit_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(APAR_FIELD,&
                     & MOMENTUM_FLUX,is), description_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(APAR_FIELD,&
                     & MOMENTUM_FLUX,is), comments_key, not_avail, binary_fmt)
                end if
              end if
            end if

            if (xy_fluxes_bpar) then
              call open_real_lu(get_k_flux_luname(BPAR_FIELD, ENERGY_FLUX, &
                 & is),'diagnostic/diagnos_fluxes', &
                 & xyslice_shape, binary_fmt, &
                 & lun_flux_k2(BPAR_FIELD, ENERGY_FLUX, is))
              if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_k2(BPAR_FIELD,&
                     & ENERGY_FLUX,is), 'time', &
                     & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & ENERGY_FLUX,is), phys_unit_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & ENERGY_FLUX,is), description_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & ENERGY_FLUX,is), comments_key, not_avail, binary_fmt)
                end if
                
              if (xy_fluxes_P) then
                call open_real_lu(get_k_flux_luname(BPAR_FIELD, PARTICLE_FLUX, &
                   & is),'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, binary_fmt, &
                   & lun_flux_k2(BPAR_FIELD, PARTICLE_FLUX, is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_k2(BPAR_FIELD,&
                     & PARTICLE_FLUX,is), 'time', &
                     & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & PARTICLE_FLUX,is), phys_unit_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & PARTICLE_FLUX,is), description_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & PARTICLE_FLUX,is), comments_key, not_avail, binary_fmt)
                end if
              end if
              if (xy_fluxes_v) then
                call open_real_lu(get_k_flux_luname(BPAR_FIELD, MOMENTUM_FLUX, &
                   & is),'diagnostic/diagnos_fluxes', &
                   & xyslice_shape, binary_fmt, &
                   & lun_flux_k2(BPAR_FIELD, MOMENTUM_FLUX, is))
                if(.not.io_legacy) then
                  call attach_metadata_grid(lun_flux_k2(BPAR_FIELD,&
                     & MOMENTUM_FLUX,is), 'time', &
                     & 'nmod', 'n_x_grid', 'number_of_species', binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & MOMENTUM_FLUX,is), phys_unit_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & MOMENTUM_FLUX,is), description_key, not_avail, binary_fmt)
                  call attach_metadata(lun_flux_k2(BPAR_FIELD,&
                     & MOMENTUM_FLUX,is), comments_key, not_avail, binary_fmt)
                end if
              end if
            end if
          end if
        end do
      end if

      if(allocated(xyslice_shape)) deallocate(xyslice_shape)

    end if

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use grid, only : nmod,nymod,nx,nsp,number_of_species,n_x_grid, n_y_grid
    use grid, only : ns
    use general, only : gkw_abort
    use components, only : tearingmode, tm_drive
    use diagnos_generic, only : mphit, mphiw3t, mrad_l, mrad_G
    use diagnos_generic, only : xy_fluxes, xy_fluxes_bi
    use diagnos_generic, only : allocate_spec_cmpx_buffers
    use control, only : io_legacy, nlphi, nlapar, nlbpar
    use rho_par_switch, only : lflux_rhostar
    use diagnos_fluxes_vspace, only : allocate_flux_det
    integer :: ierr

    ierr=0

    call allocate_flux_det

    if(io_legacy) then
      allocate(rdum(mphit,mrad_G), stat = ierr)
    else
      allocate(rdum_yxsp(n_y_grid,n_x_grid, number_of_species), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('Could not allocate rdum in diagnostic')
    endif
    
    if(io_legacy) then
      allocate(a(mphiw3t,mrad_l), stat = ierr)
    else
      allocate(a(nymod,nx), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('diagnos_fluxes: Could not allocate a in diagnostic')
    endif

    if(io_legacy) then
      allocate(b(mphiw3t,mrad_l), stat = ierr) 
    else
      allocate(b(nymod,nx), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('Could not allocate b in diagnostic')
    endif

    if(io_legacy) then
      allocate(ar(mphit,mrad_l), stat = ierr) 
    else
      allocate(ar(n_y_grid,nx), stat = ierr) 
    end if
    if (ierr.ne.0) then 
      call gkw_abort('Could not allocate ar in diagnostic')
    endif

    if(io_legacy) then
      allocate(br(mphit,mrad_l), stat = ierr) 
    else
      allocate(br(n_y_grid,nx), stat = ierr) 
    end if
    if (ierr.ne.0) then 
      call gkw_abort('Could not allocate br in diagnostic')
    endif

    if(io_legacy) then
      allocate(cr(mphit,mrad_l), stat = ierr)
    else
      allocate(cr_yxsp(n_y_grid,nx,nsp), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('Could not allocate cr in diagnostic')
    endif

    if(xy_fluxes) then
      if(io_legacy) then
        allocate(c_xy(nmod,nx),stat=ierr)
        if (ierr /= 0) call gkw_abort('diagnostic :: c_xy')
      else
        allocate(c_xysp(nymod,nx,nsp),stat=ierr)
        if (ierr /= 0) call gkw_abort('diagnostic :: c_xysp')
      end if
    end if

    if(flux3d) then
      call allocate_spec_cmpx_buffers()
    end if

    allocate(pflux_es(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_es')
    allocate(eflux_es(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_es')
    allocate(vflux_es(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_es')
    allocate(deflux_es(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: deflux_es')
    allocate(dvflux_es(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: dvflux_es')
    allocate(pflux_apar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_apar')
    allocate(eflux_apar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_apar')
    allocate(vflux_apar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_apar')
    allocate(deflux_apar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: deflux_apar')
    allocate(dvflux_apar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: dvflux_apar')
    allocate(pflux_bpar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_bpar')
    allocate(eflux_bpar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_bpar') 
    allocate(vflux_bpar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_bpar') 
    allocate(deflux_bpar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: deflux_bpar') 
    allocate(dvflux_bpar(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: dvflux_bpar') 
    allocate(fluxbuf(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: fluxbuf')  
    allocate(pflux_tot_es(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_tot_es')
    allocate(pflux_tot_apar(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_tot_apar')
    allocate(pflux_tot_bpar(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_tot_bpar')
    allocate(eflux_tot_es(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_tot_es')
    allocate(eflux_tot_apar(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_tot_apar')
    allocate(eflux_tot_bpar(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_tot_bpar')
    allocate(eflux_tot_es_lab(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_tot_es_lab')
    allocate(eflux_tot_apar_lab(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_tot_apar_lab')
    allocate(eflux_tot_bpar_lab(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_tot_bpar_lab')
    allocate(vflux_tot_es(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_tot_es')
    allocate(vflux_tot_apar(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_tot_apar')
    allocate(vflux_tot_bpar(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_tot_bpar')
    allocate(vflux_tot_es_lab(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_tot_es_lab')
    allocate(vflux_tot_apar_lab(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_tot_apar_lab')
    allocate(vflux_tot_bpar_lab(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_tot_bpar_lab')
    allocate(pflux_nc(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_nc')
    allocate(eflux_nc(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_nc')
    allocate(vflux_nc(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_nc')
    allocate(pflux_nc_old(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pflux_nc_old')
    allocate(eflux_nc_old(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux_nc_old')
    allocate(eflux2_nc_old(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: eflux2_nc_old')
    allocate(vflux_nc_old(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vflux_nc_old')
    allocate(fluxncbuf(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: fluxncbuf')
    allocate(flux_tot_es(nfluxes*number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: flux_tot_es')
    allocate(flux_tot_apar(nfluxes*number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: flux_tot_apar')
    allocate(flux_tot_bpar(nfluxes*number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: flux_tot_bpar')
    allocate(flux_tot_es_lab(nfluxes*number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: flux_tot_es_lab')
    allocate(flux_tot_apar_lab(nfluxes*number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: flux_tot_apar_lab')
    allocate(flux_tot_bpar_lab(nfluxes*number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: flux_tot_bpar_lab')

    ! Arrays for spectral fluxes
    if (lfluxes_spectra) then
      allocate(pflux_spec(nmod,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: pflux_spec')
      allocate(eflux_spec(nmod,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: eflux_spec')
      allocate(vflux_spec(nmod,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: vflux_spec')
      allocate(pflux_xspec(n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: pflux_xspec')
      allocate(eflux_xspec(n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: eflux_xspec')
      allocate(vflux_xspec(n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: vflux_xspec')
    end if

    if (lfluxes_em_spectra) then  
      allocate(pflux_apar_spec(nmod,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: pflux_apar_spec')
      allocate(eflux_apar_spec(nmod,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: eflux_apar_spec')
      allocate(vflux_apar_spec(nmod,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: vflux_apar_spec')
      allocate(pflux_apar_xspec(nx,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: pflux_apar_xspec')
      allocate(eflux_apar_xspec(nx,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: eflux_apar_xspec')
      allocate(vflux_apar_xspec(nx,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: vflux_apar_xspec')
    end if

    if ((tearingmode.or.tm_drive).and.nlapar) then
      allocate(deltaprime(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: deltaprime')
      allocate(torque(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: torque')
      allocate(deltap_isl(n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: deltap_isl')
      allocate(torque_isl(n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: torque_isl')
    end if

    if(xy_fluxes) then
      allocate(i_eflux_rad_xy(number_of_species), stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: i_eflux_rad_xy')
      
      if(xy_fluxes_bi) then
        allocate(i_eflux_pol_xy(number_of_species), stat=ierr)
        if (ierr /= 0) call gkw_abort('diagnostic :: i_eflux_pol_xy')
      end if
      
    end if

    allocate(lun_flux_xy(2,3,3,2,number_of_species), stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: lun_flux_xy')
    allocate(lun_flux_k2(3,3,number_of_species), stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: lun_flux_k2')
    allocate(flux_xy_enabled(2,3,3,2,number_of_species), stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: lun_flux_enabled')

    ! rho* correction of fluxes requires s-derivative of the fields
    if (lflux_rhostar) then
      if (nlphi) allocate(dphi_gads(ns))
      if (nlapar) allocate(dapar_gads(ns))
      if (nlbpar) allocate(dbpar_gads(ns))
    endif

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    use rho_par_switch, only : lflux_rhostar
    if(allocated(rdum)) deallocate(rdum)
    if(allocated(rdum_yxsp)) deallocate(rdum_yxsp)

    if(allocated(a)) deallocate(a)
    if(allocated(ar)) deallocate(ar)
    if(allocated(b)) deallocate(b)
    if(allocated(br)) deallocate(br)
    if(allocated(cr)) deallocate(cr)
    if(allocated(cr_yxsp)) deallocate(cr_yxsp)
    
    if(allocated(c_xy)) deallocate(c_xy)
    if(allocated(c_xysp)) deallocate(c_xysp)

    if (lflux_rhostar)  then
      if(allocated(dphi_gads)) deallocate(dphi_gads)
      if(allocated(dapar_gads)) deallocate(dapar_gads)
      if(allocated(dbpar_gads)) deallocate(dbpar_gads)
    endif

    ! ... and many others ...
    

  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()

  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine final_output(number)
    use diagnos_generic, only : xy_estep
    integer, intent(in) :: number

    if (.not.xy_estep) then
      call output_xy(number)
    end if
    
    if (.not.xy_estep) then
      call output_3d(number)
    end if

  end subroutine final_output

  !--------------------------------------------------------------------------
  !>
  !> This routine calculates the fluxes of particles, energy and  
  !> toroidal angular momentum
  !>
  !>
  !> Normalisation is such that the flux in the gradient of psi 
  !> (R_ref \nabla \psi, because psi is normalised) is 
  !>
  !> Gamma (R_r nabla psi)    = n_s rho^*2 v_thref pflux 
  !> Q_s (R_r nabla psi)      = n_s rho^*2 v_thref T_s eflux 
  !> Pi_phi(R_r nabla psi) = n_s rho^*2 v_thref m_s R_ref v_ths vflux 
  !>
  !> So to get the real flux one still needs to multiply with the 
  !> density and or temperature of the particular species 
  !>
  !> The normalization, however, has the advantage that 
  !>
  !> D       = rho^*2 Ln vthref pflux 
  !> chi     = rho^*2 Lt vthref eflux 
  !> chi_phi = rho^*2 Lu vthref v_R vflux 
  !>
  !> where the L refers to the gradient length. Note that the momentum 
  !> flux must be multiplied with the relative velocity 
  !>
  !> The routine calculates both the anomalous as well as the neo-
  !> classical fluxes. The anomalous flux is furthermore split in 
  !> the contributions due to the ExB velocity and the magnetic 
  !> flutter. 
  !--------------------------------------------------------------------------
  subroutine calc_fluxes(magnitude)
    use control,        only : nlapar
    use control,        only : nlphi, nlbpar
    use grid,           only : nx, ns, nmu, nvpar, nsp, n_x_grid 
    use grid,           only : nmod, number_of_species, gsp, gx
    use dist,           only : fdisi
    use dist,           only : iapar
    use geom,           only : ints, bn, efun, signB, bt_frac, Rfun
    use mode,           only : krho, kxrh
    use components,     only : tmp, vthrat, signz, mas, tgrid
    use components,     only : tearingmode, rhostar, tm_drive
    use rotation,       only : vcor
    use velocitygrid,   only : intmu, intvp, mugr, vpgr
    use constants,      only : ci1
    use matdat,         only : get_f_from_g 
    use index_function, only : indx 
    use fields,         only : get_averaged_phi, get_averaged_apar
    use fields,         only : get_averaged_bpar 
    use rho_par_switch, only : lflux_rhostar
    use mpiinterface,   only : mpiallreduce_sum_inplace
    use mpiinterface,   only : number_of_processors
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    use diagnos_generic, only : dfieldds, parseval_correction
    use global,         only : gkw_a_equal_b_accuracy

    !> The argument is a switch for calculating flux suprema assuming 90
    !> degree phase angle when this result is used to normalise the
    !> actual fluxes.
    !> It extracts the phase angle (e.g. for quasilinear comparisions) .
    !> In future, may also want to compute the total flux as complex number
    !> so that average phase can be computed AFTER all the averages.
    logical, intent(in) :: magnitude

    ! integers for the loop over all grid points 
    integer :: imod, ix, i, j, k, is 

    ! the gyro-averaged fields 
    complex :: phi_ga, apar_ga, bpar_ga 

    ! toggles the rhostar correction in the momentum flux
    ! TODO make input parameter

    ! Dummy variables 
    complex :: dum, dumes1, dumes2, dumem1, dumem2, dumbpar1, dumbpar2
    complex :: dumes_rs, dumem_rs, dumbpar_rs
    complex :: dumem1_cos
    complex :: fdis, aparallel
    real    :: omega_times_time
    real :: d3v

    ! The global species (ix) index is in isglb (ixg)
    integer :: isglb, ixg

    ! Initialize the fluxes to zero 
    deflux_es = 0.    ; dvflux_es = 0.
    deflux_apar = 0.  ; dvflux_apar   = 0.
    deflux_bpar = 0.  ; dvflux_bpar   = 0.

    if ((tm_drive.or.tearingmode).and. nlapar) then
      deltap_isl=0.
      torque_isl=0.
    end if

    !  Calculate the particle flux
    nmod1: do imod = 1, nmod 
      nx1:  do ix = 1, nx 
        nsp1: do is = 1, nsp

          ! the actual (global) species index
          isglb = gsp(is)
          ! the actual (global) x index
          ixg = gx(ix)

          ! Integral over the velocity space 
          nmu3: do j = 1, nmu
            ! rho* correction requires derivative of the fields 
            if (lflux_rhostar ) then
              if (nlphi)  call dfieldds(PHI_GA_FIELD,imod,ix,j,is,dphi_gads)
              if (nlapar) call dfieldds(APAR_GA_FIELD,imod,ix,j,is,dapar_gads)
              if (nlbpar) call dfieldds(BPAR_GA_FIELD,imod,ix,j,is,dbpar_gads)
            endif
            nvpar3: do k = 1, nvpar

              ! Do the average over the flux surface
              ns3: do i = 1, ns

                ! the gyro-averaged fields 
                phi_ga  = get_averaged_phi(imod,ix,i,j,is,fdisi) 
                apar_ga = get_averaged_apar(imod,ix,i,j,is,fdisi)
                bpar_ga = get_averaged_bpar(imod,ix,i,j,is,fdisi) 

                ! fdis is the distribution without A_par contribution  
                fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)

                ! If assuming 90 degree phase angle in fluxes,
                ! take the magnitudes of the complex numbers  
                if (magnitude) then
                  phi_ga  = abs(phi_ga)
                  apar_ga = abs(apar_ga) 
                  bpar_ga = abs(bpar_ga)
                  fdis    = -ci1*abs(fdis)
                end if


                ! in the implicit scheme fdis can be NaN for intvp = 0 
                if (gkw_a_equal_b_accuracy(intvp(i,j,k,is), 0.0)) fdis = 0.
                ! the 1 in the last index of efun indicates that this
                ! is the radial component of the respective fluxes
                d3v = bn(ix,i)*intvp(i,j,k,is)*intmu(j)
                dum  = d3v*parseval_correction(imod)*ints(i)*&
                   & (efun(ix,i,1,1)*kxrh(ix) +&    ! YC: E^psi^psi=0 by definition. 
                                                    ! so this line could be removed 
                   & efun(ix,i,2,1)*krho(imod))*fdis

                dumes1 = dum*conjg(phi_ga)
                dumes2 = dum*conjg(phi_ga)*bn(ix,i)
                dumem1 = -2.E0*vthrat(is)*vpgr(i,j,k,is)*dum*              &
                   & conjg(apar_ga)
                dumem2 = -2.E0*vthrat(is)*vpgr(i,j,k,is)*dum*              &
                   & conjg(apar_ga)*bn(ix,i)
                dumbpar1 = 2.E0*mugr(j)*tmp(ix,is)*dum*                 &
                   & conjg(bpar_ga)/signz(is)
                dumbpar2 = 2.E0*mugr(j)*tmp(ix,is)*dum*                 &
                   & conjg(bpar_ga)*bn(ix,i)/signz(is)

                if (lflux_rhostar) then
                  if (nlphi) then
                    dumes_rs = d3v*parseval_correction(imod)*ints(i)* &
                       & efun(ix,i,3,1) * fdis &
                       & * conjg(dphi_gads(i))
                    !the fluxes requires the real part of the rho*
                    !correction. below the imaginary part of dumes[12]
                    !is used. therefore I multiply with ci1
                    dumes1 = dumes1 + dumes_rs * rhostar * ci1
                    dumes2 = dumes2 + rhostar * dumes_rs*bn(ix,i) * ci1
                  endif
                  if( nlapar) then
                    dumem_rs = -d3v*parseval_correction(imod)*ints(i)* &
                       & efun(ix,i,1,3) * fdis *&
                       & vthrat(is)*vpgr(i,j,k,is) * conjg(dapar_gads(i))
                    dumem1 = dumem1 + rhostar * dumem_rs * ci1
                    dumem2 = dumem2 + rhostar * dumem_rs * bn(ix,i) * ci1
                  endif
                  if (nlbpar) then
                    ! there was no bn in dumbpar_rs, in contrast to
                    ! the other two _rs variables, therefore divide by
                    ! bn.
                    dumbpar_rs = d3v*parseval_correction(imod)*ints(i)* &
                       & efun(ix,i,1,3) * fdis * &
                       & mugr(j)* tmp(ix,is)* conjg(dbpar_gads(i))&
                       & /(signz(is)*bn(ix,i))
                    dumbpar1 = dumbpar1 + rhostar*dumbpar_rs * ci1
                    dumbpar2 = dumbpar2 + rhostar*dumbpar_rs*bn(ix,i) * ci1
                  endif
                endif

                deflux_es(imod,ixg,isglb) = deflux_es(imod,ixg,isglb) +       &
                   & aimag(dumes1)*vpgr(i,j,k,is)*Rfun(ix,i)*bt_frac(ix,i)*2*vcor/vthrat(is) + &
                   & aimag(dumes1)/2*mas(is)*(Rfun(ix,i)*vcor)**2/tgrid(is)
                dvflux_es(imod,ixg,isglb)=dvflux_es(imod,ixg,isglb) +        &
                   & aimag(dumes1)*(Rfun(ix,i)*bt_frac(ix,i))**2*signB*vcor/vthrat(is)
                deflux_apar(imod,ixg,isglb) = deflux_apar(imod,ixg,isglb) +       &
                   & aimag(dumem1)*vpgr(i,j,k,is)*Rfun(ix,i)*bt_frac(ix,i)*2*vcor/vthrat(is) + &
                   & aimag(dumem1)/2*mas(is)*(Rfun(ix,i)*vcor)**2/tgrid(is)
                dvflux_apar(imod,ixg,isglb)=dvflux_apar(imod,ixg,isglb) +        &
                   & aimag(dumem1)*(Rfun(ix,i)*bt_frac(ix,i))**2*signB*vcor/vthrat(is)

                deflux_bpar(imod,ixg,isglb) = deflux_bpar(imod,ixg,isglb) +       &
                   & aimag(dumbpar1)*vpgr(i,j,k,is)*Rfun(ix,i)*bt_frac(ix,i)*2*vcor/vthrat(is) + &
                   & aimag(dumbpar1)/2*mas(is)*(Rfun(ix,i)*vcor)**2/tgrid(is)
                dvflux_bpar(imod,ixg,isglb)=dvflux_bpar(imod,ixg,isglb) +        &
                   & aimag(dumbpar1)*(Rfun(ix,i)*bt_frac(ix,i))**2*signB*vcor/vthrat(is)

                if (nlapar .and. (tearingmode.or.tm_drive) &
                   .and. (imod.eq.2)) then

                  !IMPORTANT!! This is the leading term, having confused the
                  !periodicidy angle of the island xi with the 
                  !binormal coordinate of GKW, may be incorrect for large islands 
                  aparallel             = fdisi(indx(iapar,imod,ix,i))
                  omega_times_time      = atan2(-aimag(aparallel),real(aparallel))
                  dumem1_cos            = signz(is)*vthrat(is)* &
                     & vpgr(i,j,k,is)*fdis* &
                     & EXP(CMPLX(0,omega_times_time))*ints(i)
                  deltap_isl(ixg,isglb) = deltap_isl(ixg,isglb)+ &
                     & real(dumem1_cos)*intmu(j)* &
                     & intvp(i,j,k,is)*bn(ix,i)
                  torque_isl(ixg,isglb) = torque_isl(ixg,isglb)- &
                     & aimag(dumem1_cos)*intmu(j)* &
                     & intvp(i,j,k,is)*bn(ix,i)
                end if

              end do ns3
            end do nvpar3
          end do nmu3
        end do nsp1
      end do nx1
    end do nmod1

    ! only when run on more than one processor (saves two copies)
    if (number_of_processors.gt.1) then 

      ! The 3D arrays ([k]y, k[x], species) are MPI reduced over s, x, sp and
      ! velocity. As the arrays have global size in x and species, the
      ! reduction has the effect of "gathering" over x and species.
      call mpiallreduce_sum_inplace(deflux_es, nmod, n_x_grid, number_of_species)
      call mpiallreduce_sum_inplace(deflux_apar, nmod, n_x_grid, number_of_species)
      call mpiallreduce_sum_inplace(deflux_bpar, nmod, n_x_grid, number_of_species)
      call mpiallreduce_sum_inplace(dvflux_es, nmod, n_x_grid, number_of_species)
      call mpiallreduce_sum_inplace(dvflux_apar, nmod, n_x_grid, number_of_species)
      call mpiallreduce_sum_inplace(dvflux_bpar, nmod, n_x_grid, number_of_species)
      if ((tearingmode.or.tm_drive).and. nlapar) then
        call mpiallreduce_sum_inplace(deltap_isl, n_x_grid, number_of_species)
        call mpiallreduce_sum_inplace(torque_isl, n_x_grid, number_of_species)
      end if

    end if

  end subroutine calc_fluxes

  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine calc_fluxes_yxsp(pflux,eflux,vflux)
    use mpiinterface, only : mpiallreduce_sum_inplace
    use grid, only : nsp, nmod, nx, ns, nmu, nvpar, n_x_grid, number_of_species
    use grid, only : gx, gsp
    use diagnos_fluxes_vspace, only : pflux_det, eflux_det, vflux_det
    real, intent(out), dimension(:,:,:) :: pflux,eflux,vflux
    integer :: imod, ix, i, j, k, is
    pflux = 0.
    eflux = 0.
    vflux = 0.

    !  Integrate
    do is = 1, nsp
      do imod = 1, nmod
        do ix = 1, nx
          do i = 1, ns

            ! velocity space
            do j = 1, nmu
              do k = 1, nvpar
                pflux(imod,gx(ix),gsp(is)) = pflux(imod,gx(ix),gsp(is))&
                   & + pflux_det(k,j,i,ix,imod,is)

                eflux(imod,gx(ix),gsp(is)) = eflux(imod,gx(ix),gsp(is))&
                   & + eflux_det(k,j,i,ix,imod,is)

                vflux(imod,gx(ix),gsp(is)) = vflux(imod,gx(ix),gsp(is))&
                   & + vflux_det(k,j,i,ix,imod,is)

              end do
            end do
          end do
        end do
      end do
    end do

    ! The 3D arrays ([k]y, k[x], species) are MPI reduced over s, x, sp and
    ! velocity. As the arrays have global size in x and species, the
    ! reduction has the effect of "gathering" over x and species.
    call mpiallreduce_sum_inplace(pflux, nmod, n_x_grid, number_of_species)
    call mpiallreduce_sum_inplace(eflux, nmod, n_x_grid, number_of_species)
    call mpiallreduce_sum_inplace(vflux, nmod, n_x_grid, number_of_species)
    
    ! if ((tearingmode.or.tm_drive).and. nlapar) then
    !   call mpiallreduce_sum_inplace(deltap_isl, n_x_grid, number_of_species)
    !   call mpiallreduce_sum_inplace(torque_isl, n_x_grid, number_of_species)
    ! end if

  end subroutine calc_fluxes_yxsp


  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine sum_total_fluxes
    use grid, only : number_of_species, nmod, n_x_grid
    use control, only : spectral_radius, flux_tube, nlapar
    use components, only : tearingmode, tm_drive
    integer :: is, imod, ix
    integer :: iflux


    species: do is = 1, number_of_species 
      pflux_tot_es(is) = 0. 
      eflux_tot_es(is) = 0. 
      vflux_tot_es(is) = 0. 
      eflux_tot_es_lab(is) = 0. 
      vflux_tot_es_lab(is) = 0. 
      pflux_tot_apar(is) = 0. 
      eflux_tot_apar(is) = 0. 
      vflux_tot_apar(is) = 0. 
      eflux_tot_apar_lab(is) = 0. 
      vflux_tot_apar_lab(is) = 0. 
      pflux_tot_bpar(is) = 0. 
      eflux_tot_bpar(is) = 0. 
      vflux_tot_bpar(is) = 0. 
      eflux_tot_bpar_lab(is) = 0. 
      vflux_tot_bpar_lab(is) = 0. 

      if ((tearingmode.or.tm_drive) .and. nlapar) then
        deltaprime(is)=0.
        torque(is)=0.

        do ix = 1, n_x_grid
          deltaprime(is) = deltaprime(is)+ deltap_isl(ix,is) 
          torque(is) = torque(is)+ torque_isl(ix,is)
        end do
      end if

      do imod = 1, nmod

        if(lfluxes_spectra) then
          pflux_spec(imod,is) = 0. 
          eflux_spec(imod,is) = 0. 
          vflux_spec(imod,is) = 0. 
        end if
        if(lfluxes_em_spectra) then
          pflux_apar_spec(imod,is) = 0. 
          eflux_apar_spec(imod,is) = 0. 
          vflux_apar_spec(imod,is) = 0. 
        end if

        do ix = 1, n_x_grid

          if(lfluxes_spectra) then
            pflux_spec(imod,is) = pflux_spec(imod,is) + pflux_es(imod,ix,is)
            eflux_spec(imod,is) = eflux_spec(imod,is) + eflux_es(imod,ix,is)
            vflux_spec(imod,is) = vflux_spec(imod,is) + vflux_es(imod,ix,is)
          end if
          if(lfluxes_em_spectra) then
            pflux_apar_spec(imod,is) = pflux_apar_spec(imod,is) + pflux_apar(imod,ix,is)
            eflux_apar_spec(imod,is) = eflux_apar_spec(imod,is) + eflux_apar(imod,ix,is)
            vflux_apar_spec(imod,is) = vflux_apar_spec(imod,is) + vflux_apar(imod,ix,is)
          end if


          ! co-moving frame fluxes
          pflux_tot_es(is) = pflux_tot_es(is) + pflux_es(imod,ix,is)
          eflux_tot_es(is) = eflux_tot_es(is) + eflux_es(imod,ix,is)
          vflux_tot_es(is) = vflux_tot_es(is) + vflux_es(imod,ix,is)

          pflux_tot_apar(is) = pflux_tot_apar(is) + pflux_apar(imod,ix,is)
          eflux_tot_apar(is) = eflux_tot_apar(is) + eflux_apar(imod,ix,is) 
          vflux_tot_apar(is) = vflux_tot_apar(is) + vflux_apar(imod,ix,is) 

          pflux_tot_bpar(is) = pflux_tot_bpar(is) + pflux_bpar(imod,ix,is)
          eflux_tot_bpar(is) = eflux_tot_bpar(is) + eflux_bpar(imod,ix,is) 
          vflux_tot_bpar(is) = vflux_tot_bpar(is) + vflux_bpar(imod,ix,is) 

          ! Lab frame fluxes
          eflux_tot_es_lab(is) = eflux_tot_es_lab(is) + eflux_es(imod,ix,is) &
             & + deflux_es(imod,ix,is)
          vflux_tot_es_lab(is) = vflux_tot_es_lab(is) &
             & + vflux_es(imod,ix,is) &
             & + dvflux_es(imod,ix,is)

          eflux_tot_apar_lab(is) = eflux_tot_apar_lab(is) + &
             & eflux_apar(imod,ix,is) + deflux_apar(imod,ix,is)
          vflux_tot_apar_lab(is) = vflux_tot_apar_lab(is) + &
             & vflux_apar(imod,ix,is) + dvflux_apar(imod,ix,is)
 
          eflux_tot_bpar_lab(is) = eflux_tot_bpar_lab(is) + &
             & eflux_bpar(imod,ix,is) + deflux_bpar(imod,ix,is)
          vflux_tot_bpar_lab(is) = vflux_tot_bpar_lab(is) + &
             & vflux_bpar(imod,ix,is) + dvflux_bpar(imod,ix,is)

        end do !n_x_grid

        ! for the non-spectral case (but still flux tube) pflux_tot (as well) 
        ! as the spectrum must be divided by n_x_grid
        if (.not.spectral_radius .and. flux_tube) then 

          if(lfluxes_spectra) then
            pflux_spec(imod,is) = pflux_spec(imod,is) / n_x_grid 
            eflux_spec(imod,is) = eflux_spec(imod,is) / n_x_grid
            vflux_spec(imod,is) = vflux_spec(imod,is) / n_x_grid
          end if
          if(lfluxes_em_spectra) then
            pflux_apar_spec(imod,is) = pflux_apar_spec(imod,is) / n_x_grid 
            eflux_apar_spec(imod,is) = eflux_apar_spec(imod,is) / n_x_grid
            vflux_apar_spec(imod,is) = vflux_apar_spec(imod,is) / n_x_grid
          end if

        endif

      end do !nmod

      ! modify the total fluxes for the non-spectral case 
      if (.not.spectral_radius) then 

        pflux_tot_es(is) = pflux_tot_es(is) / n_x_grid
        eflux_tot_es(is) = eflux_tot_es(is) / n_x_grid
        eflux_tot_es_lab(is) = eflux_tot_es_lab(is) / n_x_grid
        vflux_tot_es(is) = vflux_tot_es(is) / n_x_grid
        vflux_tot_es_lab(is) = vflux_tot_es_lab(is) / n_x_grid

        pflux_tot_apar(is) = pflux_tot_apar(is) / n_x_grid
        eflux_tot_apar(is) = eflux_tot_apar(is) / n_x_grid
        eflux_tot_apar_lab(is) = eflux_tot_apar_lab(is) / n_x_grid
        vflux_tot_apar(is) = vflux_tot_apar(is) / n_x_grid
        vflux_tot_apar_lab(is) = vflux_tot_apar_lab(is) / n_x_grid

        pflux_tot_bpar(is) = pflux_tot_bpar(is) / n_x_grid
        eflux_tot_bpar(is) = eflux_tot_bpar(is) / n_x_grid
        eflux_tot_bpar_lab(is) = eflux_tot_bpar_lab(is) / n_x_grid
        vflux_tot_bpar(is) = vflux_tot_bpar(is) / n_x_grid
        vflux_tot_bpar_lab(is) = vflux_tot_bpar_lab(is) / n_x_grid

        if (tearingmode.or.tm_drive) then
          deltaprime(is) = deltaprime(is) / n_x_grid  
          torque(is) = torque(is) / n_x_grid
        end if

      endif

      !kx_spectra
      ! FJC_NON_SPECTRAL need FFT
      if (lfluxes_spectra) then
        do ix = 1, n_x_grid
          pflux_xspec(ix,is) = 0. 
          eflux_xspec(ix,is) = 0. 
          vflux_xspec(ix,is) = 0.
          do imod = 1, nmod 
            pflux_xspec(ix,is) = pflux_xspec(ix,is) + pflux_es(imod,ix,is)
            eflux_xspec(ix,is) = eflux_xspec(ix,is) + eflux_es(imod,ix,is)
            vflux_xspec(ix,is) = vflux_xspec(ix,is) + vflux_es(imod,ix,is)
          end do !nmod
        end do !nx
      end if !lfluxes_spectra

      if (lfluxes_em_spectra) then
        do ix = 1, n_x_grid
          pflux_apar_xspec(ix,is) = 0. 
          eflux_apar_xspec(ix,is) = 0. 
          vflux_apar_xspec(ix,is) = 0.
          do imod = 1, nmod 
            pflux_apar_xspec(ix,is) = pflux_apar_xspec(ix,is) + pflux_apar(imod,ix,is)
            eflux_apar_xspec(ix,is) = eflux_apar_xspec(ix,is) + eflux_apar(imod,ix,is)
            vflux_apar_xspec(ix,is) = vflux_apar_xspec(ix,is) + vflux_apar(imod,ix,is)
          end do !nmod
        end do !nx
      end if !lfluxes_em_spectra

      !if(io_legacy) then
      ! order the required fluxes for legacy output
      iflux = 1
      flux_tot_es(iflux + (is-1)*nfluxes) = pflux_tot_es(is)
      flux_tot_apar(iflux + (is-1)*nfluxes) = pflux_tot_apar(is)
      flux_tot_bpar(iflux + (is-1)*nfluxes) = pflux_tot_bpar(is)
      flux_tot_es_lab(iflux + (is-1)*nfluxes) = pflux_tot_es(is)
      flux_tot_apar_lab(iflux + (is-1)*nfluxes) = pflux_tot_apar(is)
      flux_tot_bpar_lab(iflux + (is-1)*nfluxes) = pflux_tot_bpar(is)
      iflux = iflux + 1
      flux_tot_es(iflux + (is-1)*nfluxes) = eflux_tot_es(is)
      flux_tot_apar(iflux + (is-1)*nfluxes) = eflux_tot_apar(is)
      flux_tot_bpar(iflux + (is-1)*nfluxes) = eflux_tot_bpar(is)
      flux_tot_es_lab(iflux + (is-1)*nfluxes) = eflux_tot_es_lab(is)
      flux_tot_apar_lab(iflux + (is-1)*nfluxes) = eflux_tot_apar_lab(is)
      flux_tot_bpar_lab(iflux + (is-1)*nfluxes) = eflux_tot_bpar_lab(is)
      iflux = iflux + 1
      flux_tot_es(iflux + (is-1)*nfluxes) = vflux_tot_es(is)
      flux_tot_apar(iflux + (is-1)*nfluxes) = vflux_tot_apar(is)
      flux_tot_bpar(iflux + (is-1)*nfluxes) = vflux_tot_bpar(is)
      flux_tot_es_lab(iflux + (is-1)*nfluxes) = vflux_tot_es_lab(is)
      flux_tot_apar_lab(iflux + (is-1)*nfluxes) = vflux_tot_apar_lab(is)
      flux_tot_bpar_lab(iflux + (is-1)*nfluxes) = vflux_tot_bpar_lab(is)
      !end if

    end do species

  end subroutine sum_total_fluxes

  !--------------------------------------------------------------------------
  !> check if any convergence conditions are fulfilled
  !--------------------------------------------------------------------------
  subroutine check_exit_conditions_nc
    use control, only : neoclassics, stop_me, ncqtol
    use mpiinterface, only : root_processor
    real :: std_ncq, mean_ncq
    real, dimension(flux_history_size), save :: nc_flux_history=(/6,5,4,3,2,1/)
    integer, save :: incq = 0

    integer :: i

    if (ncqtol > 0. .and. neoclassics) then
      ! fill in the array storing the nc flux history
      incq = mod(incq+1,flux_history_size)
      nc_flux_history(incq+1) = eflux_nc(1)

      ! standard deviation of the nc flux history
      std_ncq = 0.
      mean_ncq = sum(nc_flux_history)/flux_history_size
      do i=1,flux_history_size
        std_ncq = std_ncq + (nc_flux_history(i)-mean_ncq)**2
      end do
      std_ncq = sqrt(std_ncq/flux_history_size)

      stop_me = (std_ncq < ncqtol) .or. stop_me
      if (stop_me .and. root_processor) then
        write(*,*) 'Convergence of neoclassical ion heat flux reached: stop'
      end if
    end if
  end subroutine check_exit_conditions_nc

  !--------------------------------------------------------------------------
  !> check if any convergence conditions are fulfilled
  !> or if there is a NaN.
  !--------------------------------------------------------------------------
  subroutine check_exit_conditions
    use general, only : gkw_is_nan
    use control, only : nan_stop, stop_me
    use control, only : ifluxtol
    use control, only : fluxtol
    use mpiinterface, only : root_processor
    real :: std_flux, mean_flux
    !> current index for flux history
    integer, save :: influx = 0
    real, dimension(flux_history_size), save :: flux_history=(/6,5,4,3,2,1/)

    integer :: i


    if (fluxtol > 0.) then
      ! increment the index
      influx = mod(influx+1,flux_history_size)
      ! fill in the array storing the flux history
      flux_history(influx+1) = flux_tot_es(min(ifluxtol,size(flux_tot_es)))

      ! standard deviation of the flux history
      std_flux = 0.
      mean_flux = sum(flux_history)/flux_history_size
      do i=1,flux_history_size
        std_flux = std_flux + (flux_history(i)-mean_flux)**2
      end do
      std_flux = sqrt(std_flux/flux_history_size)

      ! do a check for flux convergence
      if (.not. stop_me) then
        stop_me = (std_flux/mean_flux < fluxtol) .or. stop_me
        if (stop_me .and. root_processor) write(*,*) 'Fluxes convergence reached: stop'
      end if
    end if

    ! attempt to stop the code if things are going badly wrong
    if (gkw_is_nan(eflux_tot_es(1))) then
      if (root_processor) then
        write(*,*) ' ******* Hit a NaN in fluxes: stop *******'
        write(*,*)
      end if
      nan_stop = .true.
      stop_me = .true.
    end if
  end subroutine check_exit_conditions

  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine calc_fluxes_nc
    use grid, only : nmod, nx, nsp, ns, nmu, nvpar, number_of_species
    use grid, only : gsp, gx
    use control, only : neoclassics
    use collisionop, only : mom_conservation, ene_conservation, conservation
    use index_function, only : indx
    use mode, only : ixzero, iyzero
    use velocitygrid, only : intmu, intvp, mugr, vpgr
    use dist, only : fdisi, fdis_tmp, fmaxwl, i_mom, i_ene
    use geom, only : ints, bn, efun, signB, bt_frac, Rfun 
    use geom, only : jfun
    use rotation, only : cf_trap, cf_drift, cfen, vcor
    use constants, only : pi
    use components, only : tmp, signz, mas
    use linear_terms, only : lpoisson, lneorotsource
    use matdat, only : matn, matn_e, matn_v
    use mpiinterface, only : mpiallreduce_sum_inplace
    ! integers for the loop over all grid points 
    integer :: imod, ix, i, j, k, is
    complex :: dum
    real :: dumnc
    integer :: l, imomidx
    ! The global species (ix) index is in isglb (ixg)
    integer :: isglb, ixg
    
    real :: velshift
    complex :: phi_ga
    if(.not.(neoclassics .and. conservation)) return
    pflux_nc = 0.
    eflux_nc = 0.
    vflux_nc = 0.
    phi_ga = 0.0

    !  Calculate the particle flux
    nmod1: do imod = 1, nmod 
      nx1:  do ix = 1, nx 
        nsp1: do is = 1, nsp

          ! the actual (global) species index
          isglb = gsp(is)
          ! the actual (global) x index
          ixg = gx(ix)

          ! New neoclassical diagnostics (momentum conserving term)
          ! FJC_NON_SPECTRAL: ix /= ixzero has no meaning
          neo: if (ix == ixzero .and. imod == iyzero) then
            dumnc = 0.
            ns2: do i = 1, ns 
              nmu2: do j = 1, nmu 
                nvpar2: do k = 1, nvpar
                  dumnc = 0.E0
                  dum = vpgr(i,j,k,is)**2 +  2.E0*mugr(j)*bn(ix,i) 

                  if(mom_conservation) then
                    imomidx = indx(i_mom,imod,ix,i,is) 
                    dumnc = real(fdisi(imomidx)) * fmaxwl(ix,i,j,k,is) * vpgr(i,j,k,is)
                  end if

                  if(ene_conservation) then
                    imomidx = indx(i_ene,imod,ix,i,is)
                    dumnc = dumnc + real(fdisi(imomidx) * fmaxwl(ix,i,j,k,is) * (dum-3.0E0/2.0E0))
                  end if

                  velshift = signB*bt_frac(ix,i)*rfun(ix,i)*vpgr(i,j,k,is)
                  if(.not.lneorotsource) then
                    velshift = velshift + vcor*jfun(ix,i)
                  endif
                  dumnc = - dumnc*velshift*4*pi*efun(ix,i,1,2)*bn(ix,i)/signz(is)*&
                     & sqrt(mas(is)*tmp(ix,is))*intmu(j)* intvp(i,j,k,is)*ints(i)

                  !Correction for when the centrifugal drift is included.
                  if(cf_drift.or.cf_trap) then
                    dum = dum + cfen(i,is) 
                  endif
                  if(lpoisson) then
                    dum = dum + signz(is)*phi_ga/tmp(ix,is)
                  end if

                  pflux_nc(isglb) = pflux_nc(isglb) + dumnc  
                  eflux_nc(isglb) = eflux_nc(isglb) + dumnc * real(dum-5.0E0/2.0E0)
                  vflux_nc(isglb) = vflux_nc(isglb) + dumnc * 0.5E0*velshift
                end do nvpar2
              end do nmu2
            end do ns2

            !write(*,*) 'Neoclassical particle flux (momcon)', pflux_nc(isglb)
            !write(*,*) 'Neoclassical heat flux (momcon)', eflux_nc(isglb)
            !write(*,*) 'Neoclassical momentum flux (momcon)', vflux_nc(isglb)
          end if neo


          ! need to communicate the derivatives before doing this
          ! this is done in exp integration
          ! FJC_NON_SPECTRAL: ix /= ixzero has no meaning
          neo2: if (neoclassics .and. (ix .eq. ixzero).and.(imod .eq. iyzero)) then
            dumnc = 0.
            nmatm: do l = 1, matn%nmat
              if (matn%jj(l) .eq. 0) then 
                write(*,*) l, matn%nmat
              end if
              !if (jjn(l).gt.nf) then
              !  write(*,*) 'ghost cell'
              !end if
              ! the minus is part of the factors
              ! fdis_tmp here contains f, not g (see exp_int)
              dumnc = - real(matn%mat(l) * fdis_tmp(matn%jj(l)))
              pflux_nc(isglb) = pflux_nc(isglb) + dumnc 
              dumnc = - real(matn_e%mat(l) * fdis_tmp(matn_e%jj(l)))
              eflux_nc(isglb) = eflux_nc(isglb) + dumnc
              dumnc = - real(matn_v%mat(l) * fdis_tmp(matn_v%jj(l)))
              vflux_nc(isglb) = vflux_nc(isglb) + dumnc  
            end do nmatm !loops over species!!, nmu, nvpar, and ns

            !write(*,*) 'Neoclassical particle flux (differential part)', pflux_nc(isglb)
            !write(*,*) 'Neoclassical heat flux (differential part)', eflux_nc(isglb)
            !write(*,*) 'Neoclassical momentum flux (differential part)', vflux_nc(isglb)
          end if neo2
        end do nsp1
      end do nx1
    end do nmod1

    call mpiallreduce_sum_inplace(pflux_nc,number_of_species)
    call mpiallreduce_sum_inplace(eflux_nc,number_of_species)
    call mpiallreduce_sum_inplace(vflux_nc,number_of_species)

  end subroutine calc_fluxes_nc



  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine calc_fluxes_nc_old
    use grid, only : nmod, nx, nsp, ns, nmu, nvpar, number_of_species
    use grid, only : gsp, gx
    use control, only : neoclassics
    use mode, only : ixzero, iyzero
    use velocitygrid, only : intmu, intvp, mugr, vpgr
    use dist, only : fdisi, fmaxwl
    use geom, only : ints, bn, signB, bt_frac, Rfun
    use geom, only : jfun
    use matdat, only : get_f_from_g 
    use rotation, only : cf_trap, cf_drift, cfen, vcor, coriolis
    use components, only : tmp, signz
    use linear_terms, only : drift, lpoisson, lvpgrphi, lneorotsource
    use mpiinterface, only : mpiallreduce_sum_inplace
    use fields, only : get_averaged_phi
    use global, only : gkw_a_equal_b_accuracy
    
    ! integers for the loop over all grid points 
    integer :: imod, ix, i, j, k, is
    complex :: dum
    real :: dumnc
    complex :: fdis
    real :: drift_x, drift_y, drift_z, dumbs

    ! The global species (ix) index is in isglb (ixg)
    integer :: isglb, ixg
    ! the gyro-averaged fields 
    complex :: phi_ga
    real :: velshift

    if(.not. neoclassics) return

    pflux_nc_old = 0.
    eflux_nc_old = 0.
    vflux_nc_old = 0.
    eflux2_nc_old = 0.
    
    !  Calculate the particle flux
    nmod1: do imod = 1, nmod 
      nx1:  do ix = 1, nx 
        nsp1: do is = 1, nsp

          ! the actual (global) species index
          isglb = gsp(is)
          ! the actual (global) x index
          ixg = gx(ix)

          !Neoclassical fluxes old version
          ! check if this is the 0,0 mode for which the neoclassical
          ! fluxes are calculated 
          ! FJC_NON_SPECTRAL: ix /= ixzero has no meaning
          if (ix == ixzero .and. imod == iyzero) then
            ns1: do i = 1, ns 
              nmu1: do j = 1, nmu 
                nvpar1: do k = 1, nvpar 

                  if(lneorotsource) then
                    call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,.false.,cf_drift,.true.)       
                  else
                    call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
                  endif
                  ! Distribution function without A|| correction 
                  fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)
                  if (gkw_a_equal_b_accuracy(intvp(i,j,k,is), 0.0)) fdis = 0.

                  ! The potential 
                  phi_ga  = get_averaged_phi(imod,ix,i,j,is,fdisi) 

                  !Following is a correction for the Landau term
                  if(lpoisson.and.lvpgrphi) then
                    fdis = fdis + signz(is)*fmaxwl(ix,i,j,k,is)*phi_ga/tmp(ix,is)
                  end if
                  velshift = signB*bt_frac(ix,i)*rfun(ix,i)*vpgr(i,j,k,is)
                  if(.not.lneorotsource) then
                    velshift = velshift + vcor*jfun(ix,i)
                  endif

                  ! common factors 
                  dumnc = drift_x * (bn(ix,i)*intmu(j)*intvp(i,j,k,is)/signz(is))*  &
                     & real(fdis)*ints(i)
                  dumbs = (bn(ix,i)**2*intmu(j)*intvp(i,j,k,is)*ints(i))*     &
                     & real(fdis)*vpgr(i,j,k,is)*signz(is)
                  dum = (vpgr(i,j,k,is)**2 +  2.E0*mugr(j)*bn(ix,i))            

                  !Correction for when the centrifugal drift is included.
                  if(cf_drift.or.cf_trap) then
                    dum = dum + cfen(i,is)
                  endif
                  !Further correction if a self-consistent potential is included.
                  if(lpoisson) then
                    dum = dum + signz(is)*phi_ga/tmp(ix,is)
                  end if

                  ! the fluxes
                  pflux_nc_old(isglb) = pflux_nc_old(isglb) + dumnc  
                  eflux_nc_old(isglb) = eflux_nc_old(isglb) + dumnc*real(dum)
                  vflux_nc_old(isglb) = vflux_nc_old(isglb) + dumnc * velshift
                  eflux2_nc_old(isglb) = eflux2_nc_old(isglb) + dumbs

                end do nvpar1
              end do nmu1
            end do ns1
            !write(*,*) 'Neoclassical particle flux (old version)', pflux_nc(isglb)
            !write(*,*) 'Neoclassical heat flux (old version)', eflux_nc(isglb)
            !write(*,*) 'Neoclassical momentum flux (old version)', vflux_nc(isglb)

          end if
        end do nsp1
      end do nx1
    end do nmod1

    call mpiallreduce_sum_inplace(pflux_nc_old,number_of_species)
    call mpiallreduce_sum_inplace(eflux_nc_old,number_of_species)
    call mpiallreduce_sum_inplace(eflux2_nc_old,number_of_species)
    call mpiallreduce_sum_inplace(vflux_nc_old,number_of_species)

  end subroutine calc_fluxes_nc_old

  !-------------------------------------------------
  ! Calculate the 2D fluxes vd*f s-averaged
  ! Please document me better !
  !-------------------------------------------------
  subroutine output_neoclass_efluxes_perp2d(direction)
    use dist,             only : fdisi
    use velocitygrid,     only : intvp,intmu,vpgr,mugr
    use grid,             only : nmod,nx,nmu,nvpar,ns, nsp
    use grid,             only : proc_subset
    use geom,             only : bn,ints,dfun
    use components,       only : tmp,signz  
    use global,           only : gkw_a_equal_b_accuracy
    use matdat,           only : get_f_from_g
    use mpiinterface,     only : mpiallreduce_sum_inplace
    use mpicomms,         only : COMM_SP_EQ_X_EQ
    use general,          only : gkw_abort
    use diagnos_generic,  only : xysp_output_array

    integer, intent(in) :: direction
    complex :: fdis, scalar
    integer :: j,k,imod,ix,is,i
    real :: vdrift
    integer :: lun
    
    select case(direction)
    case(RAD_DIRECTION)
      ! flux magnitude radial
      !prefix='flmgr'
      lun = i_eflux_rad_xy(1)

    case(POL_DIRECTION)
      ! flux magnitude poloidal
      !prefix='flmgp'
      lun = i_eflux_pol_xy(1)

    case default
      call gkw_abort('output_neoclass_efluxes_perp2d unknown direction')
    end select

    if(nmod == 1) c_xysp(1,:,:) = 0.0
    
    species_loop: do is = 1, nsp

      ! loop over the slice points for this species
      local_slice: do imod=1,nmod
        do ix=1,nx

          ! integrate over the velocity space for this point
          scalar=(0.E0,0.E0)

          velocity_space: do j=1,nmu
            do k=1,nvpar
              ! N.B. this is probably slow!
              flux_surf_avg: do i=1,ns

                fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)
                if (gkw_a_equal_b_accuracy(intvp(i,j,k,is), 0.0)) fdis = 0.

                ! WARNING: Should call the drift function of linear terms
                ! here instead call drift - this is incomplete with rotation

                vdrift = (vpgr(i,j,k,is)**2+(mugr(j)*bn(ix,i)))* &
                   & dfun(ix,i,direction)*tmp(ix,is)/signz(is)
                scalar = scalar + vdrift*bn(ix,i)*intvp(i,j,k,is)* &
                   & intmu(j)* fdis * ints(i)

              end do flux_surf_avg
            end do
          end do velocity_space
          
          ! store the scalar in a temporary slice for reduction later
          if(nmod == 1) then
            c_xysp(imod+1,ix, is) = scalar
          else
            c_xysp(imod,ix, is) = scalar
          end if

        end do
      end do local_slice
    end do species_loop

    !Sum the slice over other processors if necessary
    call mpiallreduce_sum_inplace(c_xysp,shape(c_xysp), &
       & COMM_SP_EQ_X_EQ)

    ! All processes now have the x-sp-local array.
    ! Only on one x slice (with global s index 1): FFT, gather and
    ! output to a logical unit.
    ! (and generate a luname based on the global species)

    call xysp_output_array( &
       & c_xysp, &
       & lun, &
       & proc_subset(0,1,1,1,0), .true.)

  end subroutine output_neoclass_efluxes_perp2d

  !-------------------------------------------------
  ! 
  !-------------------------------------------------
  subroutine outp_nc_efluxes_perp2d_legacy(direction)
    use dist,             only : fdisi
    use velocitygrid,     only : intvp,intmu,vpgr,mugr
    use grid,             only : nmod,nx,nmu,nvpar,ns, lsp, nsp
    use grid,             only : proc_subset, number_of_species
    use geom,             only : bn,ints,dfun
    use components,       only : tmp,signz  
    use global,           only : gkw_a_equal_b_accuracy
    use matdat,           only : get_f_from_g
    use mpiinterface,     only : mpiallreduce_sum_inplace
    use mpicomms,         only : COMM_SP_EQ_X_EQ
    use general,          only : gkw_abort
    use diagnos_generic,  only : xy_output_array

    integer, intent(in) :: direction
    complex :: fdis, scalar
    integer :: j,k,imod,ix,is,i
    real :: vdrift
    integer :: lun

    ! loop over the global species
    species_loop: do is = 1, number_of_species

      select case(direction)
      case(RAD_DIRECTION)
        ! flux magnitude radial
        !prefix='flmgr'
        lun = i_eflux_rad_xy(is)

      case(POL_DIRECTION)
        ! flux magnitude poloidal
        !prefix='flmgp'
        lun = i_eflux_pol_xy(is)

      case default
        call gkw_abort('output_neoclass_efluxes_perp2d unknown direction')
      end select

      if(proc_subset(0,0,0,0,is)) then

        ! loop over the slice points for this species
        local_slice: do imod=1,nmod
          do ix=1,nx

            ! integrate over the velocity space for this point
            scalar=(0.E0,0.E0)

            velocity_space: do j=1,nmu
              do k=1,nvpar
                ! N.B. this is probably slow!
                flux_surf_avg: do i=1,ns

                  fdis = get_f_from_g(imod,ix,i,j,k,lsp(is),fdisi)
                  if (gkw_a_equal_b_accuracy(intvp(i,j,k,lsp(is)), 0.0)) fdis = 0.

                  ! WARNING: Should call the drift function of linear terms
                  ! here instead call drift - this is incomplete with rotation

                  vdrift = (vpgr(i,j,k,lsp(is))**2+(mugr(j)*bn(ix,i)))* &
                     & dfun(ix,i,direction)*tmp(ix,lsp(is))/signz(lsp(is))
                  scalar = scalar + vdrift*bn(ix,i)*intvp(i,j,k,lsp(is))* &
                     & intmu(j)* fdis * ints(i)

                end do flux_surf_avg
              end do
            end do velocity_space
            ! store the scalar in a temporary slice for reduction later
            c_xy(imod,ix) = scalar

          end do
        end do local_slice

        !Sum the slice over other processors if necessary
        call mpiallreduce_sum_inplace(c_xy,shape(c_xy), &
           & COMM_SP_EQ_X_EQ)

      end if

      ! All processes now have the x-local array.
      ! Only on one x slice (with global s index 1): FFT, gather and
      ! output to a logical unit.
      ! (and generate a luname based on the global species)
      call xy_output_array( &
         & c_xy, .false., &
         & lun, &
         & proc_subset(0,1,1,1,is), is <= nsp)

    end do species_loop
  end subroutine outp_nc_efluxes_perp2d_legacy
  
!-----------------------------------------------------------------------------
!> Calculation of fluxes in real space for 2D output
!> This in independent from the spectral flux calculation as it operates 
!> in real space.  This is necessary because the spectral fluxes calculation 
!> loses the phase (and therefore spatial) information when it multiplies the
!> complex conjugates. For the flux surface average, the two methods agree
!> i.e. the average of the XY_flux exactly matches the total flux calculation

!> This is a generalised routine for any direction and any field and moment
!> and for a single slice and for the flux surface average
!> The input arguments are used to select which flux should be calculated:
!>   direc  selects between radial and poloidal directions 
!>   fluxf  selects between phi, apa or bpar.
!>   moment selects between particle, energy, velocity flux
!>   fs_av  selects between flux surface average or single slice 
!>
!> The way this routine is called repeatedly is a bit inefficient, since
!> the FFTs are repeated unnecessarily.  Note however, the diagnostics timings 
!> will warn the user is if they become a significant fraction of runtime
!-----------------------------------------------------------------------------
subroutine xy_flux_calculation(fdis,direc,field,moment,fsa,gsp)

  use dist,             only : nsolc
  use velocitygrid,     only : intvp,intmu,vpgr,mugr
  use grid,             only : nmod,nx,nmu,nvpar,ls, nsp,ns,ls, lsp
  use grid,             only : proc_subset, n_y_grid
  use control,          only : time
  use geom,             only : bn, efun, ints
  use mode,             only : krho, kxrh
  use components,       only : tmp,vthrat,signz, tearingmode, isl_mode
  use io,               only : xy_fmt, binary_fmt
  use tearingmodes,     only : omega_rot
  use non_linear_terms, only : jind 
  use constants,        only : ci1
  use matdat,           only : get_f_from_g 
  use fields,           only : get_averaged_phi, get_averaged_apar 
  use fields,           only : get_averaged_bpar
  use mpiinterface,     only : gather_array, mpiallreduce_sum_inplace
  use mpiinterface,     only : root_processor
  use mpiinterface,     only : send_to_root, mpibarrier
  use mpicomms,         only : COMM_VPAR_NE_MU_NE, COMM_SP_EQ_X_EQ 
  use mpicomms,         only : COMM_X_NE, COMM_DUMMY, COMM_SP_NE_X_NE
!  use mpicomms,         only : COMM_X_EQ
  use general,          only : gkw_abort
  use diagnos_generic,  only : xy_slice_ipar, lisl_follow
  use diagnos_generic,  only : four2real_2D
  use fft, only : FFT_INVERSE
  use grid,             only : jind_flexible
  use diagnos_generic,  only : mphit, mrad_G, mrad_l, mpi_dtype_yxsp_slice
  use control, only : io_legacy
  use io, only : append_chunk, xy_fmt, binary_fmt
  use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
  use diagnos_generic, only : PARTICLE_FLUX, ENERGY_FLUX, MOMENTUM_FLUX
  use global,          only : gkw_a_equal_b_accuracy

  complex, dimension(nsolc), intent(in) :: fdis
  integer, intent(in) :: direc     !< poloidal or radial
  integer, intent(in) :: field     !< which field
  integer, intent(in) :: moment    !< which moment
  integer, intent(in) :: fsa     !< if flux surface average
  !> global species index, remove this when io_legacy is removed.
  integer, intent(in) :: gsp
  integer :: i,j,k,imod,ix, iy, isp
  !> a helper index. Since for nmod == 1 the complex data must not be
  !> stored into index imod=1 of the buffer, but into imod=2, because
  !> imod=1 is the zero-mode (which is not simulated if nmod == 1)
  integer :: imod_in_cbuf
  integer :: i_begin, i_end, isp_begin, isp_end
  real    :: vsqr = 0.0, tens1 = 0.0, tens2 = 0.0
  complex :: fdisi
  complex :: phi_ga, apar_ga, bpar_ga
  integer :: reordered_ix
  logical :: root_is_in_comm
  integer :: comm
  logical :: proc_is_involved_in_integr
  logical :: proc_receives_gather, proc_in_slice

  if (io_legacy) then
    isp_begin = lsp(gsp)
    isp_end = lsp(gsp)
    cr = 0
  else
    if(gsp /= 0) call gkw_abort('gsp should be zero! Coding mistake. &
       & species loop is inside this routine, not outside.')
    isp_begin = 1
    isp_end = nsp
    cr_yxsp = 0

  end if
  
  ! Select the points on the field line needed
  if (fsa == NO_FSAVG) then
    ! if no flux surface avg is demanded, pick a single s point
    ! local s index
    i_begin = ls(xy_slice_ipar) 
    i_end = ls(xy_slice_ipar)
    proc_is_involved_in_integr = proc_subset(0,xy_slice_ipar,0,0,gsp)
    proc_receives_gather = proc_subset(1,xy_slice_ipar,1,1,gsp)
    proc_in_slice = proc_subset(0,xy_slice_ipar,1,1,gsp)

    root_is_in_comm = (xy_slice_ipar <= ns)
  else
    ! Setting for flux surface average
    i_begin = 1
    i_end = ns
    proc_is_involved_in_integr = proc_subset(0,0,0,0,gsp)
    proc_receives_gather = proc_subset(1,1,1,1,gsp)
    proc_in_slice = proc_subset(0,1,1,1,gsp)

    root_is_in_comm = .true.
  end if

  procs_involved: if (proc_is_involved_in_integr) then
    species_loop: do isp = isp_begin, isp_end

      spoints: do i = i_begin, i_end
        velspace: do k = 1, nvpar
          do j = 1, nmu
            a = (0.,0.)
            b = (0.,0.)
            ar= 0.0
            br= 0.0

            ! loop over the slice points for this species
            perp_slice: do imod=1,nmod
              if(nmod == 1) then
                imod_in_cbuf = imod + 1
              else
                imod_in_cbuf = imod
              end if
              do ix=1,nx

                select case(direc)
                case(RAD_DIRECTION)
                  tens1 = efun(ix,i,1,1)
                  tens2 = efun(ix,i,2,1)
                case(POL_DIRECTION)
                  tens1 = efun(ix,i,1,2)
                  tens2 = 0.
                end select

                ! fdis is the distribution without A_par contribution  
                fdisi = get_f_from_g(imod,ix,i,j,k,isp,fdis) 
                if (gkw_a_equal_b_accuracy(intvp(i,j,k,isp), 0.0)) fdisi = 0.
                !Just the distribution function
                if(io_legacy) then
                  reordered_ix = jind(ix)
                else
                  reordered_ix = jind_flexible(nx,nx,ix)
                end if
                b(imod_in_cbuf,reordered_ix) = fdisi

                select case(field)
                case(PHI_FIELD)
                  ! here is a factor of two less than the spectral calculation
                  ! this is because only the spectral calculation needs
                  ! the complex conjugate conversion
                  ! |f|^2 = (1/2) |f||f*|
                  phi_ga = get_averaged_phi(imod,ix,i,j,isp,fdis)
                  a(imod_in_cbuf,reordered_ix)=ci1*(tens1*kxrh(ix) +           &
                     & tens2*krho(imod))*bn(ix,i)*phi_ga
                case(APAR_FIELD)
                  apar_ga = get_averaged_apar(imod,ix,i,j,isp,fdis) 
                  a(imod_in_cbuf,reordered_ix)=-2.E0*ci1*(tens1*kxrh(ix) +          &
                     & tens2*krho(imod))*  &
                     & bn(ix,i)*vthrat(isp)*vpgr(i,j,k,isp)*apar_ga  
                case(BPAR_FIELD)
                  bpar_ga = get_averaged_bpar(imod,ix,i,j,isp,fdis) 
                  a(imod_in_cbuf,reordered_ix)=2.E0*ci1*(tens1*kxrh(ix) +           &
                     & tens2*krho(imod))*  &
                     & bn(ix,i)*tmp(ix,isp)*mugr(j)*bpar_ga/signz(isp)
                case default
                  call gkw_abort('Wrong field in xy fluxes calculation')  
                end select

                if(tearingmode.and.lisl_follow)then
                  b(imod_in_cbuf,reordered_ix) = b(imod_in_cbuf,reordered_ix) * &
                     & exp(cmplx(0.0, (imod-1)*omega_rot*time/(isl_mode-1)))
                  a(imod_in_cbuf,reordered_ix) = a(imod_in_cbuf,reordered_ix) * &
                     & exp(cmplx(0.0, (imod-1)*omega_rot*time/(isl_mode-1)))
                end if

              end do
            end do perp_slice

            ! transform to position space
            call four2real_2D(ar,a,FFT_INVERSE)
            call four2real_2D(br,b,FFT_INVERSE)


            if(io_legacy) then
              ! Choose the moment to calculate
              select case(moment)
              case(PARTICLE_FLUX)
                vsqr = 1.0

              case(ENERGY_FLUX)
                vsqr = vpgr(i,j,k,isp)**2+(2.0*mugr(j)*bn(1,i))

              case(MOMENTUM_FLUX)
                !! FIXME: Still needs RBt/B factors to agree quantatively with total
                vsqr = vpgr(i,j,k,isp)
              case default
                call gkw_abort('Wrong moment in xy fluxes calculation')
              end select

              do ix = 1, mrad_l
                do iy = 1, mphit
                  cr(iy, ix) = cr(iy, ix)+ar(iy,ix)*br(iy,ix)*  &
                     & vsqr*intvp(i,j,k,isp)*intmu(j)
                end do
              end do
            else
              do ix = 1, nx
                ! Choose the moment to calculate
                select case(moment)
                case(PARTICLE_FLUX)
                  vsqr = 1.0

                case(ENERGY_FLUX)
                  vsqr = vpgr(i,j,k,isp)**2+(2.0*mugr(j)*bn(ix,i))

                case(MOMENTUM_FLUX)
                  !! FIXME: Still needs RBt/B factors to agree quantatively with total
                  vsqr = vpgr(i,j,k,isp)
                case default
                  call gkw_abort('Wrong moment in xy fluxes calculation')
                end select

                do iy = 1, n_y_grid
                  cr_yxsp(iy, ix,isp) = cr_yxsp(iy, ix, isp)+ar(iy,ix)*br(iy,ix)*  &
                     & vsqr*intvp(i,j,k,isp)*intmu(j)
                end do
              end do
            end if

          end do
        end do velspace
      end do spoints
    end do species_loop

    ! finish the velspace integral and evtl. the flux surface avg
    if (fsa == FSAVG) then
      comm = COMM_SP_EQ_X_EQ
      if(io_legacy) then
        cr = cr*ints(1)
      else
        cr_yxsp = cr_yxsp*ints(1)
      end if
    else
      comm = COMM_VPAR_NE_MU_NE
    end if

    if(io_legacy) then
      call mpiallreduce_sum_inplace(cr,size(cr,1), size(cr,2),comm)
    else
      call mpiallreduce_sum_inplace(cr_yxsp, &
         & size(cr_yxsp,1), size(cr_yxsp,2), size(cr_yxsp,3),comm)
    end if
  end if procs_involved

  if(io_legacy) then
    ! gather global array in x
    !if(proc_is_involved_in_integr) then ! FIXME WRONG?
    if(proc_in_slice) then
      call gather_array(rdum,mphit,mrad_G,  &
         & cr,mphit,mrad_l, &
         & COMM_DUMMY,COMM_X_NE,ALLGATHER=.true.)
    else
      rdum = 0
    end if
    if((root_processor .and. .not. proc_receives_gather) &
       & .or. (proc_receives_gather .and. .not. root_processor)) then
      call send_to_root(rdum, shape(rdum))
    end if
    if (root_processor) then
      call append_chunk(lun_flux_xy(direc, field, moment, fsa, gsp), &
         & rdum, xy_fmt, binary_fmt)
    end if
  else
    !gather in x and  species
    if(proc_in_slice .or. root_processor) then
      call gather_array(rdum_yxsp, cr_yxsp, mpi_dtype_yxsp_slice, COMM_SP_NE_X_NE, &
         & .true., root_is_in_comm)
    end if
    if(root_processor) then
      call append_chunk(lun_flux_xy(direc, field, moment, fsa, 1), &
         & rdum_yxsp, xy_fmt, binary_fmt)
    end if
  end if
  call mpibarrier()
end subroutine xy_flux_calculation

!--------------------------------------------------------------------
!>
!--------------------------------------------------------------------
subroutine output(magnitude,file_count)
  use diagnos_generic, only : xy_estep
  
  logical, intent(in) :: magnitude
  integer, intent(in) :: file_count
  
  call output_regular(magnitude)

  if (.not. magnitude .and. xy_estep) then
    call output_xy(file_count)
  end if

  
  if(flux3d .and. xy_estep) then
    call output_3d(file_count)
  end if

end subroutine output



!--------------------------------------------------------------------
!>
!--------------------------------------------------------------------
subroutine output_regular(magnitude)
  use collisionop, only : conservation
  use io, only : append_chunk, xy_fmt, ascii_fmt
  use control, only : io_legacy, method
  use control, only : nlphi, nlapar, neoclassics, nlbpar
  use diagnos_generic, only : lradial_profile
  use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
  use grid, only : n_x_grid, nmod, number_of_species
  use components, only : tearingmode, tm_drive
  use mpiinterface, only : root_processor
  use diagnos_fluxes_vspace, only : calc_fluxes_full_detail
  use control, only : io_format
  ! to use the other routine:
  
  logical, intent(in) :: magnitude
  integer :: is
  if(.not. lcalc_fluxes) return
  if(magnitude .and..not. lfluxes_spectra)  return
  
  ! either:
  call calc_fluxes(magnitude)
  
  ! or use the other routine:
  call calc_fluxes_full_detail(PHI_FIELD, magnitude)
  call calc_fluxes_yxsp(pflux_es, eflux_es, vflux_es)
  call calc_fluxes_full_detail(APAR_FIELD, magnitude)
  call calc_fluxes_yxsp(pflux_apar, eflux_apar, vflux_apar)
  call calc_fluxes_full_detail(BPAR_FIELD, magnitude)
  call calc_fluxes_yxsp(pflux_bpar, eflux_bpar, vflux_bpar)

  call sum_total_fluxes
  
  magnitude_cond: if (magnitude) then
    ! Recompute fluxes suprema.

    ! WARNING: do use any results computed
    ! in subroutine fluxes except fluxes suprema after the call of
    ! calc_fluxes(.true.)
    ! above
    
    ! output only fluxes suprema
    if (root_processor) then
      if(lfluxes_spectra) then
        call append_chunk(i_efluxmag, &
           & reshape(eflux_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_pfluxmag, &
           & reshape(pflux_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
      end if

      if (lfluxes_em_spectra) then
        call append_chunk(i_efluxmag_apar, &
           & reshape(eflux_apar_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_pfluxmag_apar, &
           & reshape(pflux_apar_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
      end if
    end if
  else
    call calc_fluxes_nc_old
    call calc_fluxes_nc

    if(method /= 'EIV') then
      call check_exit_conditions
      call check_exit_conditions_nc
    end if
    
    if (neoclassics) then
      call append_chunk(i_neoclass_old, &
         & (/ (pflux_nc_old(is), eflux_nc_old(is), vflux_nc_old(is), &
         & eflux2_nc_old(is), is = 1, number_of_species) /), xy_fmt, ascii_fmt)
      if(conservation) then
        call append_chunk(i_neoclass, &
           & (/ (pflux_nc(is), eflux_nc(is), vflux_nc(is), &
           & is = 1, number_of_species) /), xy_fmt, ascii_fmt)
      end if
    end if

    ! for non-fluxtube cases write the complete radial profile to file   
    if (lradial_profile) then 
    
      ! electrostatic contribution to flux due to ExB-drift
      if(nlphi) then
        if(io_format == 'hdf5') then
          call append_chunk(i_pflux_rad_es, &
             & reshape(sum(pflux_es,1),(/n_x_grid , number_of_species/)), &
             & xy_fmt, ascii_fmt)
          call append_chunk(i_vflux_rad_es, &
             & reshape(sum(vflux_es,1),(/n_x_grid , number_of_species/)), &
             & xy_fmt, ascii_fmt)
          call append_chunk(i_eflux_rad_es, &
             & reshape(sum(eflux_es,1),(/n_x_grid , number_of_species/)), &
             & xy_fmt, ascii_fmt)
        else
          call append_chunk(i_pflux_rad_es, &
             & reshape(sum(pflux_es,1),(/n_x_grid * number_of_species/)), &
             & xy_fmt, ascii_fmt)
          call append_chunk(i_vflux_rad_es, &
             & reshape(sum(vflux_es,1),(/n_x_grid * number_of_species/)), &
             & xy_fmt, ascii_fmt)
          call append_chunk(i_eflux_rad_es, &
             & reshape(sum(eflux_es,1),(/n_x_grid * number_of_species/)), &
             & xy_fmt, ascii_fmt)
        end if
      end if
         
      ! electromagnetic contribution to flux due to magnetic flutter
      if(nlapar) then
        if(io_format == 'hdf5') then
          call append_chunk(i_pflux_rad_apar, &
            & reshape(sum(pflux_apar,1),(/n_x_grid , number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_vflux_rad_apar, &
            & reshape(sum(vflux_apar,1),(/n_x_grid , number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_eflux_rad_apar, &
            & reshape(sum(eflux_apar,1),(/n_x_grid , number_of_species/)), &
            & xy_fmt, ascii_fmt)
        else
          call append_chunk(i_pflux_rad_apar, &
            & reshape(sum(pflux_apar,1),(/n_x_grid * number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_vflux_rad_apar, &
            & reshape(sum(vflux_apar,1),(/n_x_grid * number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_eflux_rad_apar, &
            & reshape(sum(eflux_apar,1),(/n_x_grid * number_of_species/)), &
            & xy_fmt, ascii_fmt)
        end if
      end if
      
      ! electromagnetic contribution to flux due to magnetic compression
      if(nlbpar) then
        if(io_format == 'hdf5') then
          call append_chunk(i_pflux_rad_bpar, &
            & reshape(sum(pflux_bpar,1),(/n_x_grid , number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_vflux_rad_bpar, &
            & reshape(sum(vflux_bpar,1),(/n_x_grid , number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_eflux_rad_bpar, &
            & reshape(sum(eflux_bpar,1),(/n_x_grid , number_of_species/)), &
            & xy_fmt, ascii_fmt)
        else
          call append_chunk(i_pflux_rad_bpar, &
            & reshape(sum(pflux_bpar,1),(/n_x_grid * number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_vflux_rad_bpar, &
            & reshape(sum(vflux_bpar,1),(/n_x_grid * number_of_species/)), &
            & xy_fmt, ascii_fmt)
          call append_chunk(i_eflux_rad_bpar, &
            & reshape(sum(eflux_bpar,1),(/n_x_grid * number_of_species/)), &
            & xy_fmt, ascii_fmt)
        end if
      end if
      
    endif

    if(io_legacy) then
      if (nlphi) then
        call append_chunk(i_fluxes_legacy, flux_tot_es, xy_fmt, ascii_fmt)
        call append_chunk(i_fluxes_legacy_lab, flux_tot_es_lab, xy_fmt, ascii_fmt)
      end if
      if (nlapar) then
        call append_chunk(i_fluxes_em, flux_tot_apar, xy_fmt, ascii_fmt)
        call append_chunk(i_fluxes_em_lab, flux_tot_apar_lab, xy_fmt, ascii_fmt)
      end if
      if (nlbpar) then
        call append_chunk(i_fluxes_bpar, flux_tot_bpar, xy_fmt, ascii_fmt)
        call append_chunk(i_fluxes_bpar_lab, flux_tot_bpar_lab, xy_fmt, ascii_fmt)
      end if
    else
      if (nlphi) then
        call append_chunk(i_pflux_es, pflux_tot_es, xy_fmt, ascii_fmt)
        call append_chunk(i_eflux_es, eflux_tot_es, xy_fmt, ascii_fmt)
        call append_chunk(i_vflux_es, vflux_tot_es, xy_fmt, ascii_fmt)
        call append_chunk(i_pflux_es_lab, pflux_tot_es, xy_fmt, ascii_fmt)
        call append_chunk(i_eflux_es_lab, eflux_tot_es_lab, xy_fmt, ascii_fmt)
        call append_chunk(i_vflux_es_lab, vflux_tot_es_lab, xy_fmt, ascii_fmt)
      end if
      if (nlapar) then
        call append_chunk(i_pflux_apar, pflux_tot_apar, xy_fmt, ascii_fmt)
        call append_chunk(i_eflux_apar, eflux_tot_apar, xy_fmt, ascii_fmt)
        call append_chunk(i_vflux_apar, vflux_tot_apar, xy_fmt, ascii_fmt)
        call append_chunk(i_pflux_apar_lab, pflux_tot_apar, xy_fmt, ascii_fmt)
        call append_chunk(i_eflux_apar_lab, eflux_tot_apar_lab, xy_fmt, ascii_fmt)
        call append_chunk(i_vflux_apar_lab, vflux_tot_apar_lab, xy_fmt, ascii_fmt)
      end if
      if (nlbpar) then
        call append_chunk(i_pflux_bpar, pflux_tot_bpar, xy_fmt, ascii_fmt)
        call append_chunk(i_eflux_bpar, eflux_tot_bpar, xy_fmt, ascii_fmt)
        call append_chunk(i_vflux_bpar, vflux_tot_bpar, xy_fmt, ascii_fmt)
        call append_chunk(i_pflux_bpar_lab, pflux_tot_bpar, xy_fmt, ascii_fmt)
        call append_chunk(i_eflux_bpar_lab, eflux_tot_bpar_lab, xy_fmt, ascii_fmt)
        call append_chunk(i_vflux_bpar_lab, vflux_tot_bpar_lab, xy_fmt, ascii_fmt)
      end if
    end if

    if ((tearingmode.or.tm_drive).and.nlapar) then
      call append_chunk(i_deltaprime, deltaprime, xy_fmt, ascii_fmt)
      call append_chunk(i_torque, torque, xy_fmt, ascii_fmt)
    end if

    ! fluxes spectra
    fluxes_spectra : if(lfluxes_spectra) then
      !Must have 1024 > nmod * number_of_species
      !Possibility of a fortran column limit?

      if (nmod > 1) then
        call append_chunk(i_efluxspec, &
           & reshape(eflux_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_pfluxspec, &
           & reshape(pflux_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_vfluxspec, &
           & reshape(vflux_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
      end if
      if (n_x_grid > 1) then 
        call append_chunk(i_efluxxspec, &
           & reshape(eflux_xspec,(/n_x_grid * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_pfluxxspec, &
           & reshape(pflux_xspec,(/n_x_grid * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_vfluxxspec, &
           & reshape(vflux_xspec,(/n_x_grid * number_of_species/)), &
           & xy_fmt, ascii_fmt)
      end if

    end if fluxes_spectra

    ! EM fluxes spectra
    fluxes_em_spectra : if(lfluxes_em_spectra) then

      !Must have 1024 > nmod * number_of_species
      !Possibility of a fortran column limit?
      if (nmod > 1) then
        call append_chunk(i_efluxspec_apar, &
           & reshape(eflux_apar_spec,(/nmod * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_pfluxspec_apar, &
           & reshape(pflux_apar_spec,(/nmod * number_of_species/)), &
           xy_fmt, ascii_fmt)
        call append_chunk(i_vfluxspec_apar, &
           & reshape(vflux_apar_spec,(/nmod * number_of_species/)), &
           xy_fmt, ascii_fmt)
      end if
      if (n_x_grid > 1) then 
        call append_chunk(i_efluxxspec_apar, &
           & reshape(eflux_apar_xspec,(/n_x_grid * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_pfluxxspec_apar, &
           & reshape(pflux_apar_xspec,(/n_x_grid * number_of_species/)), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_vfluxxspec_apar, &
           & reshape(vflux_apar_xspec,(/n_x_grid * number_of_species/)), &
           & xy_fmt, ascii_fmt)
      end if

    end if fluxes_em_spectra
  end if magnitude_cond

end subroutine output_regular

!--------------------------------------------------------------------
!>
!--------------------------------------------------------------------
subroutine output_3d(file_count)
  use grid, only : number_of_species
  integer, intent(in) :: file_count
  integer :: imoment
  
  if(flux3d) then
      do imoment = 1,3
          call output_fluxes_kykxs_parallel(imoment, &
             & file_count, &
             & tag_range_start + (imoment-1)*number_of_species)
        end do
        !TO DO attach metadata
  end if

end subroutine output_3d

!--------------------------------------------------------------------
!>
!--------------------------------------------------------------------
subroutine output_fluxes_kykxs_parallel(imoment, file_count, tag)
    use io,               only : mpi_output_array, xy_fmt, binary_fmt
    use grid,             only : nmod, nx, ns, proc_subset, nmu, nvpar
    use grid,             only : lsp, number_of_species, nsp
    use mpicomms,         only : COMM_S_NE_X_NE, COMM_VPAR_NE_MU_NE
    use mpiinterface,     only : mpiallreduce_sum_inplace
    use global,           only : int2char_zeros
    use diagnos_generic,  only : mpi_dtype_cpx_spec_yxs, out3dbuf_local_spec_cmpx
    use diagnos_generic,  only : out3dbuf_global_spec_cmpx
    use diagnos_fluxes_vspace, only : pflux_det, eflux_det, vflux_det
    use diagnos_generic, only : PARTICLE_FLUX, ENERGY_FLUX, MOMENTUM_FLUX
    use diagnos_fluxes_vspace, only : calc_fluxes_full_detail
    use global,                only : EVERY_FIELD

    integer, intent(in) :: imoment, file_count, tag

    integer :: imod,ix,isg, isp
    
    character (len=30) :: luname
    character (len=15) :: prefix
    logical :: root_is_in_comm
    integer :: i, j, k
    real, pointer :: flux(:,:,:,:,:,:)
    
    ! calculate p/e/vflux_det, otherwise the pointer flux hold zeros only
    call calc_fluxes_full_detail(EVERY_FIELD, .false.)

    prefix='?Fluxk'
    select case(imoment)
    case(PARTICLE_FLUX)
      prefix(1:1) = 'P'
      flux => pflux_det
    case(ENERGY_FLUX)
      prefix(1:1) = 'E'
      flux => eflux_det
    case(MOMENTUM_FLUX)
      prefix(1:1) = 'V'
      flux => vflux_det
    end select

    species_loop : do isg = 1, number_of_species
      if(proc_subset(0,0,0,0,isg)) then
        out3dbuf_local_spec_cmpx = 0.0
        isp = lsp(isg)
        !  Integrate
        do imod = 1, nmod
          do ix = 1, nx
            do i = 1, ns
              ! velocity space
              do j = 1, nmu
                do k = 1, nvpar
                  out3dbuf_local_spec_cmpx(imod,ix,i) = &
                     & out3dbuf_local_spec_cmpx(imod,ix,i) + &
                     & flux(k,j,i,ix,imod,isp)
                end do
              end do
            end do
          end do
        end do
        ! reduce over processes which work on the same spatial
        ! point+species;
        call mpiallreduce_sum_inplace(out3dbuf_local_spec_cmpx, &
           & shape(out3dbuf_local_spec_cmpx), COMM_VPAR_NE_MU_NE)
      end if

      ! create a luname
      luname = trim(prefix)//trim(int2char_zeros(isg,2))//'_'// &
         & trim(int2char_zeros(file_count,6))

      !on an arbitrary process, find out if root (of comm_cart)
      !works on this global species.
      root_is_in_comm = (isg <= nsp)
      ! *all* processes go into IO...
      call mpi_output_array(luname, 'diagnostic/diagnos_fluxes', &
         & out3dbuf_local_spec_cmpx, mpi_dtype_cpx_spec_yxs, &
         & out3dbuf_global_spec_cmpx, &
         & COMM_S_NE_X_NE, xy_fmt, binary_fmt, proc_subset(0,0,1,1,isg), &
         & root_is_in_comm, tag+isg-1)

    end do species_loop
end subroutine output_fluxes_kykxs_parallel


!--------------------------------------------------------------------
!>
!--------------------------------------------------------------------
subroutine output_xy(file_count)
  use mode, only : mode_box
  use dist, only : fdisi
  use grid, only : number_of_species
  use control, only : io_legacy
  use diagnos_generic, only : lphi_diagnostics
  use diagnos_generic, only : xy_fluxes, xy_fluxes_bi
  use diagnos_generic, only : xy_fluxes_v, xy_fluxes_k
  use diagnos_generic, only : xy_fluxes_p, xy_fluxes_bpar, xy_fluxes_em
  use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
  use diagnos_generic, only : PARTICLE_FLUX, ENERGY_FLUX, MOMENTUM_FLUX
  use io, only : append_chunk, binary_fmt, xy_fmt
  integer, intent(in) :: file_count

  integer :: direction, moment, field, is, fsa

  ! To keep the compiler quiet.
  if (file_count > 0) continue

  if (mode_box .and. lphi_diagnostics) then

    if (xy_fluxes) then
      if(io_legacy) then
        call outp_nc_efluxes_perp2d_legacy(RAD_DIRECTION)
      else
        call output_neoclass_efluxes_perp2d(RAD_DIRECTION)
      end if
    end if

    if (xy_fluxes_bi .and. xy_fluxes) then
      if(io_legacy) then
        call outp_nc_efluxes_perp2d_legacy(POL_DIRECTION)
      else
        call output_neoclass_efluxes_perp2d(POL_DIRECTION)
      end if
    end if

    do direction = 1, 2
      do field = 1, 3
        do moment = 1, 3
          do fsa = FSAVG, NO_FSAVG
            if(io_legacy) then
              do is = 1, number_of_species
                if(flux_xy_enabled(direction, field, moment, fsa, is)) then
                  call xy_flux_calculation(fdisi, direction, field, moment, &
                     & fsa, is)
                end if
              end do
            else
              if(flux_xy_enabled(direction, field, moment, fsa, 1)) then
                ! pass a zero as last arg, do the species loop inside
                ! the routine:
                call xy_flux_calculation(fdisi, direction, field, moment, &
                   & fsa, 0)
              end if
            end if
          end do
        end do
      end do
    end do

    if(io_legacy) then
      do is = 1, number_of_species

        ! 2D spectral outputs (these are all flux surface averages)
        if (xy_fluxes_k) then
          call append_chunk(&
             & lun_flux_k2(PHI_FIELD, ENERGY_FLUX, is), &
             & eflux_es(:,:,is), xy_fmt, binary_fmt)
          if (xy_fluxes_p) call append_chunk(&
             & lun_flux_k2(PHI_FIELD, PARTICLE_FLUX, is), &
             & pflux_es(:,:,is), xy_fmt, binary_fmt)
          if (xy_fluxes_v) call append_chunk(&
             & lun_flux_k2(PHI_FIELD, MOMENTUM_FLUX, is), &
             & vflux_es(:,:,is), xy_fmt, binary_fmt)

          if(xy_fluxes_em) then
            call append_chunk(&
               & lun_flux_k2(APAR_FIELD, ENERGY_FLUX, is), &
               & eflux_es(:,:,is), xy_fmt, binary_fmt)
            if (xy_fluxes_p) call append_chunk(&
               & lun_flux_k2(APAR_FIELD, PARTICLE_FLUX, is), &
               & pflux_es(:,:,is), xy_fmt, binary_fmt)
            if (xy_fluxes_v) call append_chunk(&
               & lun_flux_k2(APAR_FIELD, MOMENTUM_FLUX, is), &
               & vflux_es(:,:,is), xy_fmt, binary_fmt)
          end if

          if(xy_fluxes_bpar) then
            call append_chunk(&
               & lun_flux_k2(BPAR_FIELD, ENERGY_FLUX, is), &
               & eflux_es(:,:,is), xy_fmt, binary_fmt)
            if (xy_fluxes_p) call append_chunk(&
               & lun_flux_k2(BPAR_FIELD, PARTICLE_FLUX, is), &
               & pflux_es(:,:,is), xy_fmt, binary_fmt)
            if (xy_fluxes_v) call append_chunk(&
               & lun_flux_k2(BPAR_FIELD, MOMENTUM_FLUX, is), &
               & vflux_es(:,:,is), xy_fmt, binary_fmt)
          end if
        end if
      end do
    else
      ! write the data into one lu per quantity, species becomes a dimension

      ! 2D spectral outputs (these are all flux surface averages)
      if (xy_fluxes_k) then
        call append_chunk(&
           & lun_flux_k2(PHI_FIELD, ENERGY_FLUX, 1), &
           & eflux_es(:,:,:), xy_fmt, binary_fmt)
        if (xy_fluxes_P) call append_chunk(&
           & lun_flux_k2(PHI_FIELD, PARTICLE_FLUX, 1), &
           & pflux_es(:,:,:), xy_fmt, binary_fmt)
        if (xy_fluxes_v) call append_chunk(&
           & lun_flux_k2(PHI_FIELD, MOMENTUM_FLUX, 1), &
           & vflux_es(:,:,:), xy_fmt, binary_fmt)

        if(xy_fluxes_em) then
          call append_chunk(&
             & lun_flux_k2(APAR_FIELD, ENERGY_FLUX, 1), &
             & eflux_es(:,:,:), xy_fmt, binary_fmt)
          if (xy_fluxes_P) call append_chunk(&
             & lun_flux_k2(APAR_FIELD, PARTICLE_FLUX, 1), &
             & pflux_es(:,:,:), xy_fmt, binary_fmt)
          if (xy_fluxes_v) call append_chunk(&
             & lun_flux_k2(APAR_FIELD, MOMENTUM_FLUX, 1), &
             & vflux_es(:,:,:), xy_fmt, binary_fmt)
        end if

        if(xy_fluxes_bpar) then
          call append_chunk(&
             & lun_flux_k2(BPAR_FIELD, ENERGY_FLUX, 1), &
             & eflux_es(:,:,:), xy_fmt, binary_fmt)
          if (xy_fluxes_P) call append_chunk(&
             & lun_flux_k2(BPAR_FIELD, PARTICLE_FLUX, 1), &
             & pflux_es(:,:,:), xy_fmt, binary_fmt)
          if (xy_fluxes_v) call append_chunk(&
             & lun_flux_k2(BPAR_FIELD, MOMENTUM_FLUX, 1), &
             & vflux_es(:,:,:), xy_fmt, binary_fmt)
        end if
      end if
    end if
  end if
end subroutine output_xy

!--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine screen_output()
    use mpiinterface, only : root_processor
    use collisionop, only : conservation
    use control,   only : non_linear, neoclassics
    use grid,      only : number_of_species, nmod, nx
    use mode,      only : mode_box

    integer :: is, imod, ix

    ! Only root process shall write to terminal
    if(.not.root_processor) return

    if(non_linear) then
      write(*,*)
      do is = 1, number_of_species
        write(*,*)'Species ',is
        write(*,*)'Particle flux : ',pflux_tot_es(is)
        write(*,*)'Energy flux   : ',eflux_tot_es(is)
        write(*,*)'Momentum flux : ',vflux_tot_es(is)
        write(*,*)
      end do

    else
      if(mode_box) then
        do is = 1, number_of_species
          write(*,10)is 
          write(*,20)pflux_tot_es(is)
          write(*,30)eflux_tot_es(is)
          write(*,40)vflux_tot_es(is)
        end do
10      format('Species  : ',I4)
20      format('Total particle flux :',es13.5)
30      format('Total energy flux   :',es13.5)
40      format('Total momentum flux :',es13.5)

        if(neoclassics .and. conservation) then
          do is = 1, number_of_species
            write(*,10)is
            write(*,34)pflux_nc(is)
            write(*,35)eflux_nc(is)
            write(*,36)vflux_nc(is)
          end do
34        format('Total neoclassical particle flux :',es13.5)
35        format('Total neoclassical energy flux   :',es13.5)
36        format('Total neoclassical momentum flux :',es13.5)
        end if
      else
        mod : do imod = 1, nmod 
          x : do ix = 1, nx 
            species : do is = 1, number_of_species  
              !Note this output block generates zeros for nlapar = false
              !Since the electromagnetic fluxes are zero
              write(*,1)imod,ix,is
1             format('Toroidal mode ',I3,' Radial mode ',I3,' Species ',I2) 
              write(*,2)pflux_es(imod,ix,is), pflux_apar(imod,ix,is), pflux_bpar(imod,ix,is) 
2             format('The particle flux  (ES/APAR/BPAR) : ',3(es13.5,1X)) 
              write(*,3)eflux_es(imod,ix,is), eflux_apar(imod,ix,is), eflux_bpar(imod,ix,is) 
3             format('The energy flux    (ES/APAR/BPAR) : ',3(es13.5,1X)) 
              write(*,4)vflux_es(imod,ix,is), vflux_apar(imod,ix,is),vflux_bpar(imod,ix,is)
4             format('The momentum flux  (ES/APAR/BPAR) : ',3(es13.5,1X))
              write(*,*)
            end do species
          end do x
        end do mod
      end if
    end if

  end subroutine screen_output


  !--------------------------------------------------------------------
  !> generate a luname based on the global species
  !--------------------------------------------------------------------
  function get_xy_flux_luname(direction, field, moment, fsa, gsp)
    use general, only : gkw_abort
    use global, only : int2char_zeros
    use control, only : io_legacy
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use diagnos_generic, only : PARTICLE_FLUX, ENERGY_FLUX, MOMENTUM_FLUX
    integer, intent(in) :: field, direction, moment, fsa, gsp
    character (len=15) :: get_xy_flux_luname

    select case(direction)

    case(RAD_DIRECTION)

      !Choses the tensor components
      !tens1 = efun(ipar,1,1)
      !tens2 = efun(ipar,2,1)
      !Changes the name according to the field
      select case(field)
      case(PHI_FIELD)
        get_xy_flux_luname = 'FFlesr'
      case(APAR_FIELD)
        get_xy_flux_luname = 'FFlemr'
      case(BPAR_FIELD)
        get_xy_flux_luname = 'FFlbpr'
      case default

        call gkw_abort('Wrong field in xy fluxes calculation')

      end select

    case(POL_DIRECTION)

      !Choses the tensor components
      !tens1 = efun(ipar,1,2)
      !tens2 = 0.E0
      select case(field)
      case(PHI_FIELD)
        get_xy_flux_luname = 'FFlesp'
      case(APAR_FIELD)
        get_xy_flux_luname = 'FFlemp'
      case(BPAR_FIELD)
        get_xy_flux_luname = 'FFlbpp'
      case default

        call gkw_abort('Wrong field in xy fluxes calculation')

      end select

    case default

      call gkw_abort('Wrong direction choice in xy fluxes calculation')

    end select

    ! Overwrite the first letter of the flux based on moment
    select case(moment)
    case(PARTICLE_FLUX)  ! Particle (PFlesr, PFlesp, PFlemr, PFlemp)
      get_xy_flux_luname(1:1) = 'P'
    case(ENERGY_FLUX)  ! Energy   (EFlesr, EFlesp, EFlemr, EFlemp)
      get_xy_flux_luname(1:1) = 'E'
    case(MOMENTUM_FLUX)  ! Momentum (VFlesr, VFlesp, VFlemr, VFlemp)
      get_xy_flux_luname(1:1) = 'V'
    case default
      call gkw_abort('Wrong moment in xy fluxes calculation')

    end select

    ! Overwrite the third letter of the flux for fsa, gives lunames:
    ! PFAesr, PFAesp, PFAemr, PFAemp
    ! EFAesr, EFAesp, EFAemr, EFAemp 
    ! VFAesr, VFAesp, VFAemr, VFAemp
    if (fsa == FSAVG) get_xy_flux_luname(3:3) = 'A'

    if(io_legacy) then
      get_xy_flux_luname = trim(get_xy_flux_luname)//int2char_zeros(gsp, 4)
    end if

  end function get_xy_flux_luname

  !--------------------------------------------------------------------
  !> generate a luname based on the global species
  !--------------------------------------------------------------------
  function get_k_flux_luname(field, moment, gsp)
    use general, only : gkw_abort
    use global, only : int2char_zeros
    use control, only : io_legacy
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use diagnos_generic, only : PARTICLE_FLUX, ENERGY_FLUX, MOMENTUM_FLUX
    integer, intent(in) :: field, moment, gsp
    character (len=15) :: get_k_flux_luname

    get_k_flux_luname = '?Fl??k'

    ! Overwrite the first letter of the flux based on moment
    select case(moment)
    case(PARTICLE_FLUX)
      get_k_flux_luname(1:1) = 'P'
    case(ENERGY_FLUX)
      get_k_flux_luname(1:1) = 'E'
    case(MOMENTUM_FLUX)
      get_k_flux_luname(1:1) = 'V'
    case default
      call gkw_abort('Wrong moment in k fluxes calculation')
    end select
    
    select case(field)
    case(PHI_FIELD)
      get_k_flux_luname(4:5) = 'es'
    case(APAR_FIELD)
      get_k_flux_luname(4:5) = 'em'
    case(BPAR_FIELD)
      get_k_flux_luname(4:5) = 'bp'
    case default
      call gkw_abort('Wrong field in k fluxes calculation')
    end select

    if(io_legacy) then
      get_k_flux_luname = trim(get_k_flux_luname)//int2char_zeros(gsp, 4)
    end if
    
  end function get_k_flux_luname



end module diagnos_fluxes
