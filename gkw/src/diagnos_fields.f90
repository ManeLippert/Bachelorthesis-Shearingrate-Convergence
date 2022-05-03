!------------------------------------------------------------------------------
!>
!> This diagnostic calculates and outputs data which is closely
!> related to the electromagnetic fields, i.e. the fluctuations of
!> the electrostatic potential,the vector potential and the
!> magnetic field.
!>
!> (A) kx- and ky-spectra of phi^2 and other perturbed fields
!> (B) A spectrum suitable especially for the Rosenbluth-Hinton test
!> (C) kx- and ky-spectra of the vorticity
!> (D) Parallel structure of the fields
!> (D) The 3D fields for the whole flux tube
!>
!------------------------------------------------------------------------------
module diagnos_fields

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output

  public :: gather_field_3d_serial

  !> the range of tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl

  logical, save, public :: xy_phi
  logical, save, public :: xy_apar
  logical, save, public :: xy_bpar
  logical, save, public :: xy_vort

  !> a switch to produce data typically used for the Rosenbluth-Hinton test
  logical, save, public :: field_fsa_kyzero
 
  !> logicals to select various x s diagnostic outputs
  !> (to be completed)
  logical, save, public :: xs_phi
  integer, parameter :: nmodemax = 256
  !> A list of toroidal mode numbers. The requested field is going to
  !> be output for each entry.
  integer, save, public :: imod_field(nmodemax)
  integer, save, dimension(:), allocatable :: imod_list

  !> Logicals to select ky,kx,s fields outputs
  !> Works for mode_box=T and F
  logical, save, public :: kykxs_phi
  logical, save, public :: kykxs_apar
  logical, save, public :: kykxs_bpar

  !> Obsolete namelist item: Number of toroidal points to consider (if
  !> -1 outputs all)
  integer, save, public :: nmodepoints


  !> logical unit numbers for output files
  integer, save :: i_kxvort = -1, i_kyvort = -1

  integer, save :: i_kyspec = -1, i_kxspec = -1
  integer, save :: i_kyspec_em = -1, i_kxspec_em = -1, i_rhtest = -1
  integer, save :: lun_phi_xs
  integer, save :: i_apar_fsa_kyzero= -1
  
  !> Private buffers, for FFT of perpendicular slices
  complex, save, allocatable :: a(:,:)
  real, save, allocatable :: ar(:,:)

  !> The fields spectra. ky_spec(nmod), kx_spec(nx)
  real, save, allocatable, dimension(:) :: ky_spec, kx_spec

  !> 3d local real space buffer
  real, save, allocatable :: out3dbuf_local_real(:,:,:)
  real, save, allocatable :: out3dbuf_global_real(:,:,:)
  real, save, allocatable :: out3dbuf_global_spec(:,:,:)
  complex, save, allocatable :: out3dbuf_global_cpx_spec(:,:,:)


  !> Buffer for MPI reductions, buffer(max(nx,nmod))
  real, save, allocatable, dimension(:) :: buffer

  !> arrays for the vorticity(nmod,nx)
  real, save, allocatable, dimension(:,:) :: vort
  
  !> arrays for parallel_phi, apar and bpar: field_par(n_s_grid)
  complex, save, allocatable, dimension(:) :: field_par
  complex, save, allocatable, dimension(:) :: field_par_global

  integer :: mphit, mrad_l, mrad_G
  integer :: mphiw3t

  !> logical unit numbers for 2d diagnostics
  integer, save :: lun_field_xy(3), lun_field_spec(4)

  integer, parameter, public :: VORTICITY_FIELD = 4

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    use control, only : spectral_radius, nlapar, nlbpar, nlphi
    use components, only : finit
    xy_phi = .false.
    xy_apar = .false.
    xy_bpar = .false.
    xy_vort = .false.
   
    kykxs_phi = .false.
    kykxs_apar = .false.
    kykxs_bpar = .false.
    
    xs_phi      = .false.
    imod_field = 0

    if (finit == 'zonal' .or. finit == 'gam'  .or. finit == 'alfven') then
       field_fsa_kyzero = .true.
    else
       field_fsa_kyzero = .false.
    end if
 
  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(xy_phi,1)
    call mpibcast(xy_apar,1)
    call mpibcast(xy_bpar,1)
    call mpibcast(xy_vort,1)

    call mpibcast(kykxs_phi,1)
    call mpibcast(kykxs_apar,1)
    call mpibcast(kykxs_bpar,1)

    call mpibcast(field_fsa_kyzero,1)

    call mpibcast(xs_phi,1)
    call mpibcast(imod_field,nmodemax)

    call mpibcast(field_fsa_kyzero,1)

    call mpibcast(nmodepoints,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : nlapar, spectral_radius, nlbpar, nlphi, io_format
    use general, only : gkw_warn
    use io, only : ascii_fmt
    use mode, only : iyzero
    use diagnos_generic, only : phi3d, spc3d, apa3d, apc3d, bpa3d, bpc3d
    
    if ((.not. spectral_radius) .and. (kykxs_phi .or. kykxs_apar .or. kykxs_bpar)) then
      call gkw_warn('kykxs_xxx diagnostics only implemented for spectral_radius=.true.')       
      kykxs_phi=.false.
      kykxs_apar=.false.
      kykxs_bpar=.false.
    end if

    if ((io_format == ascii_fmt).and. (kykxs_phi .or. kykxs_apar .or. kykxs_bpar)) then
      call gkw_warn('kykxs_xxx diagnostics not implemented for ascii io_format')       
      kykxs_phi=.false.
      kykxs_apar=.false.
      kykxs_bpar=.false.
    end if

    if (xs_phi .and. spectral_radius) then
      call gkw_warn('xs_phi: &
         & not implemented yet for the spectral version')
      xs_phi = .false.
    end if

    if (.not. nlphi) then
      if(xy_phi) then
        call gkw_warn('xy_phi: &
          & phi is not calculated and will not be output')
        xy_phi = .false.
      end if
      if(phi3d .or. spc3d) then
        call gkw_warn('phi3d, spc3d : phi is not calculated and will not be output')
        phi3d = .false.
        spc3d = .false.
      end if
      if (kykxs_phi) then
        call gkw_warn('kykxs_phi: phi is not calculated and will not be output')
        kykxs_phi=.false.
      end if
      if(xy_vort) then
      call gkw_warn('xy_vort: &
         & phi is not calculated and vorticity will not be output')
      xy_vort = .false.
      end if
    end if
    
    if (.not. nlapar) then
      if(xy_apar) then
        call gkw_warn('xy_apar: &
          & Apar is not calculated and will not be output')
        xy_apar = .false.
      end if
      if(apa3d .or. apc3d) then
        call gkw_warn('apa3d, apc3d : Apar is not calculated and will not be output')
        apa3d = .false.
        apc3d = .false.
      end if
      if (kykxs_apar) then
        call gkw_warn('kykxs_apar: Apar is not calculated and will not be output')
        kykxs_apar=.false.
      end if
    end if
    
    if (.not. nlbpar) then
      if(xy_bpar) then
        call gkw_warn('xy_bpar: &
           & Bpar is not calculated and will not be output')
        xy_bpar = .false.
      end if
      if(bpa3d .or. bpc3d) then
        call gkw_warn('bpa3d, bpc3d : Bpar is not calculated and will not be output')
        bpa3d = .false.
        bpc3d = .false.
      end if
      if (kykxs_bpar) then
        call gkw_warn('kykxs_bpar: Bpar is not calculated and will not be output')
        kykxs_bpar=.false.
      end if
  end if

    if (iyzero <= 0) then
      call gkw_warn('Cannot use the field_fsa_kyzero (Rosenbluth-Hinton)&
         & diagnostic: There is no k_y=0 mode.')
      field_fsa_kyzero=.false.
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor, register_tag_range
    use io, only : open_real_lu, open_complex_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use io, only : binary_fmt
    use diagnos_generic, only : attach_metadata_grid, rad_gridname
    use diagnos_generic, only : lphi_diagnostics, LOCAL_DATA
    use control, only : normalized, nlphi, nlapar
    use control, only : flux_tube, spectral_radius, io_legacy
    use mode, only : mode_box
    use grid, only : nmod, n_x_grid, n_s_grid, n_y_grid
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    logical, intent(inout) :: requirements(:,:)

    integer :: xyslice_shape(2)

    requirements(PHI_FIELD,LOCAL_DATA) = .true.
    requirements(APAR_FIELD,LOCAL_DATA) = .true.
    requirements(BPAR_FIELD,LOCAL_DATA) = .true.

    call register_tag_range(9, &
       & tag_range_start, tag_range_end_inkl)

    if (.not.root_processor) return

    if(io_legacy) then
      xyslice_shape(1) = mphit
      xyslice_shape(2) = mrad_G
    else
      xyslice_shape(1) = n_y_grid
      xyslice_shape(2) = n_x_grid
    end if

    if (.not. (mode_box .and. lphi_diagnostics)) return
    
    !if (non_linear .or. mode_box) then! legacy switch
    if (mode_box .and. lphi_diagnostics) then
      call open_real_lu('kyvort', 'diagnostic/diagnos_fields', (/ nmod /), &
         & ascii_fmt, i_kyvort)
      call open_real_lu('kxvort', 'diagnostic/diagnos_fields', (/ n_x_grid /), &
         & ascii_fmt, i_kxvort)
    end if

    if (field_fsa_kyzero .and. nlphi) then
      call open_complex_lu('rhtest', 'diagnostic/diagnos_fields', (/ n_x_grid /), &
         & ascii_fmt, i_rhtest)
      if(spectral_radius) then
         call attach_metadata_grid(i_rhtest, 'time', 'kxrh', ascii_fmt)
      else
         call attach_metadata_grid(i_rhtest, 'time', 'n_x_grid', ascii_fmt)
      end if
      call attach_metadata(i_rhtest, phys_unit_key, 'T_{ref}/e', ascii_fmt)
      call attach_metadata(i_rhtest, description_key, &
         & 'The ky=0 mode of the es. potential, flux-surface-averaged.', ascii_fmt)
      call attach_metadata(i_rhtest, comments_key, &
         & 'See also apar_fsa_kyzero', ascii_fmt)
    end if
    if (field_fsa_kyzero .and. nlapar) then
      call open_complex_lu('apar_fsa_kyzero', 'diagnostic/diagnos_fields', &
         & (/ n_x_grid /), ascii_fmt, i_apar_fsa_kyzero)
      if(spectral_radius) then
         call attach_metadata_grid(i_apar_fsa_kyzero, 'time', 'kxrh', ascii_fmt)
      else
         call attach_metadata_grid(i_apar_fsa_kyzero, 'time', 'n_x_grid', &
            & ascii_fmt)
      end if
      call attach_metadata(i_apar_fsa_kyzero, phys_unit_key, 'B_{ref}R_{ref}', &
         & ascii_fmt)
      call attach_metadata(i_apar_fsa_kyzero, description_key, &
         & 'The ky=0 mode of the parall. component of the vector potential, &
         &flux-surface-averaged.', ascii_fmt)
      call attach_metadata(i_apar_fsa_kyzero, comments_key, &
         & 'see also apar_fsa_kyzero', ascii_fmt)
    end if

    if (.not. normalized) then
      call open_real_lu('kyspec', 'diagnostic/diagnos_fields', (/ nmod /), &
         & ascii_fmt, i_kyspec)
      if(mode_box .and. flux_tube) then
        call attach_metadata_grid(i_kyspec, 'time', 'krho', ascii_fmt)
      else
        call attach_metadata_grid(i_kyspec, 'time', 'kzeta', ascii_fmt)
      end if
      call attach_metadata(i_kyspec, phys_unit_key, '(T_{ref}/e)^2', ascii_fmt)
      call attach_metadata(i_kyspec, description_key, &
         & 'binormal power spectrum of the electrostatic potential', ascii_fmt)
      call attach_metadata(i_kyspec, comments_key, &
         & 'Traditionally, the modulus square of the es. potential is used as &
         & a measure of the mode intensity.', ascii_fmt)
      
      call open_real_lu('kxspec', 'diagnostic/diagnos_fields', (/ n_x_grid /), &
         & ascii_fmt, i_kxspec)
      if(mode_box .and. spectral_radius) then
        call attach_metadata_grid(i_kxspec, 'time', 'kxrh', ascii_fmt)
      else
        call attach_metadata_grid(i_kxspec, 'time', not_avail, ascii_fmt)
      end if
      call attach_metadata(i_kxspec, phys_unit_key, '(T_{ref}/e)^2', ascii_fmt)
      call attach_metadata(i_kxspec, description_key, &
         & 'radial power spectrum of the electrostatic potential', ascii_fmt)
      call attach_metadata(i_kxspec, comments_key, &
         & 'Traditionally, the modulus square of the es. potential is used as &
         & a measure of the mode intensity.', ascii_fmt)

      if(nlapar) then
        call open_real_lu('kyspec_em', 'diagnostic/diagnos_fields', (/ nmod /), &
           & ascii_fmt, i_kyspec_em)
        call open_real_lu('kxspec_em', 'diagnostic/diagnos_fields', (/ n_x_grid /), &
           & ascii_fmt, i_kxspec_em)
      end if
    end if

    if (xy_phi) then
      !better name:luname_spec = 'phi_xky'//file_count_suffix
      call open_real_lu('spc', &
         & 'diagnostic/diagnos_fields', (/ nmod, n_x_grid /), &
         & binary_fmt, lun_field_spec(PHI_FIELD))
      call attach_metadata(lun_field_spec(PHI_FIELD), &
           & description_key, &
           & 'electrostatic potential on the complete 3d computational spatial &
           & domain, in spectral representation with respect to the 1st &
           & (binormal) dimension.', binary_fmt)
        call attach_metadata(lun_field_spec(PHI_FIELD), &
           & comments_key, not_avail, binary_fmt)
        call attach_metadata(lun_field_spec(PHI_FIELD), &
         & phys_unit_key, not_avail, &
         & binary_fmt)
        call attach_metadata_grid(lun_field_spec(PHI_FIELD), &
           & 'time', 'krho', 'xphi', binary_fmt)
        
      !better name:luname = 'phi_xy'//file_count_suffix
      call open_real_lu('phi', &
         & 'diagnostic/diagnos_fields', xyslice_shape, &
         & binary_fmt, lun_field_xy(PHI_FIELD))
      call attach_metadata(lun_field_xy(PHI_FIELD), &
           & description_key, &
           & 'electrostatic potential on the complete 3d computational spatial&
           & domain, in position space representation', binary_fmt)
        call attach_metadata(lun_field_xy(PHI_FIELD), &
           & comments_key, not_avail, binary_fmt)
        call attach_metadata(lun_field_xy(PHI_FIELD), &
           & phys_unit_key, not_avail, &
           & binary_fmt)
        call attach_metadata_grid(lun_field_xy(PHI_FIELD), &
           & 'time', 'yphi', 'xphi', binary_fmt)
    end if

    if(xy_bpar) then
      !better name:luname_spec = 'Bpar_xky'//file_count_suffix
      call open_real_lu('bpc', &
         & 'diagnostic/diagnos_fields', (/ nmod, n_x_grid /), &
         & binary_fmt, lun_field_spec(BPAR_FIELD))
      call attach_metadata(lun_field_spec(BPAR_FIELD), &
         & description_key, &
         & 'parallel component of the magn. field on the complete 3d&
         & computational spatial domain, in spectral representation with &
         & respect to the 1st (binormal) dimension.', binary_fmt)
      call attach_metadata(lun_field_spec(BPAR_FIELD), &
         & comments_key, not_avail, binary_fmt)
      call attach_metadata(lun_field_spec(BPAR_FIELD), &
         & phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata_grid(lun_field_spec(BPAR_FIELD), &
         & 'time', 'krho', 'xphi', binary_fmt)

      ! better name:luname = 'Bpar_xy'//file_count_suffix
      call open_real_lu('bpa', &
         & 'diagnostic/diagnos_fields', xyslice_shape, &
         & binary_fmt, lun_field_xy(BPAR_FIELD))
      call attach_metadata(lun_field_xy(BPAR_FIELD), &
         & description_key, &
         & 'parallel component of the magn. field on the complete 3d&
         & computational spatial domain, in position space representation', &
         & binary_fmt)
      call attach_metadata(lun_field_xy(BPAR_FIELD), &
         & comments_key, not_avail, binary_fmt)
      call attach_metadata(lun_field_xy(BPAR_FIELD), &
         & phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata_grid(lun_field_xy(BPAR_FIELD), &
         & 'time', 'yphi', 'xphi', binary_fmt)
    end if

    if(xy_apar) then
      !better name:luname = 'Apar_xy'//file_count_suffix
      !better name:luname_spec = 'Apar_xky'//file_count_suffix
      call open_real_lu('sac', &
         & 'diagnostic/diagnos_fields', (/ nmod, n_x_grid /), &
         & binary_fmt, lun_field_spec(APAR_FIELD))
      call attach_metadata(lun_field_spec(APAR_FIELD), &
         & description_key, &
         & 'parallel component of the vector potential on the complete 3d&
         & computational spatial domain, in spectral representation with &
         & respect to the 1st (binormal) dimension.', binary_fmt)
      call attach_metadata(lun_field_spec(APAR_FIELD), &
         & comments_key, not_avail, binary_fmt)
      call attach_metadata(lun_field_spec(APAR_FIELD), &
         & phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata_grid(lun_field_spec(APAR_FIELD), &
         & 'time', 'krho', 'xphi', binary_fmt)

      call open_real_lu('apa', &
         & 'diagnostic/diagnos_fields', xyslice_shape, &
         & binary_fmt, lun_field_xy(APAR_FIELD))
      call attach_metadata(lun_field_xy(APAR_FIELD), &
         & description_key, &
         & 'parallel component of the vector potential on the complete 3d&
         & computational spatial domain, in position space representation', &
         & binary_fmt)
      call attach_metadata(lun_field_xy(APAR_FIELD), &
         & comments_key, not_avail, binary_fmt)
      call attach_metadata(lun_field_xy(APAR_FIELD), &
         & phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata_grid(lun_field_xy(APAR_FIELD), &
         & 'time', 'yphi', 'xphi', binary_fmt)
    end if

    if(xy_vort) then
      !better name:luname_spec = 'vort_xky'//file_count_suffix
      call open_real_lu('vok', &
         & 'diagnostic/diagnos_fields', (/ nmod, n_x_grid /), &
         & binary_fmt, lun_field_spec(VORTICITY_FIELD))
      call attach_metadata(lun_field_spec(VORTICITY_FIELD), &
         & description_key, &
         & 'vorticity on the complete 3d&
         & computational spatial domain, in spectral representation with &
         & respect to the 1st (binormal) dimension.', binary_fmt)
      call attach_metadata(lun_field_spec(VORTICITY_FIELD), &
         & comments_key, not_avail, binary_fmt)
      call attach_metadata(lun_field_spec(VORTICITY_FIELD), &
         & phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata_grid(lun_field_spec(VORTICITY_FIELD), &
         & 'time', 'krho', 'xphi', binary_fmt)
    end if

    if(xs_phi) then
      call open_complex_lu('phi_xs', &
         & 'diagnostic/diagnos_fields', &
         & (/ n_x_grid, n_s_grid, size(imod_list) /), &
         & binary_fmt, lun_phi_xs)
      call attach_metadata(lun_phi_xs, &
         & description_key, &
         & 'electrostatic potential slice for the binormal modes specified by&
         & the parameter imod_list', binary_fmt)
      call attach_metadata(lun_phi_xs, &
         & comments_key, not_avail, binary_fmt)
      call attach_metadata(lun_phi_xs, &
         & phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata_grid(lun_phi_xs, &
         & 'time', rad_gridname, 'sgrid', &
         & 'imod_list', binary_fmt)
    end if
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem(allocate_3d_buffers_)
    use grid, only : nymod, nmod, n_x_grid, nx, n_y_grid, n_s_grid, ns, proc_subset
    use diagnos_generic, only : compute_fft_sizes
    use diagnos_generic, only : phi3d, apa3d, bpa3d, spc3d, apc3d, bpc3d
    use general, only : gkw_abort
    use mpiinterface, only : root_processor
    use control, only : io_legacy
    !> by calling with this optional argument .true., the
    !> diagnos_mode_struct diagnostic can make sure that the buffers necessary
    !> for the gather operation are allocated.
    logical, optional, intent(in) :: allocate_3d_buffers_
    integer :: ierr
    integer :: n, i
    logical :: allocate_3d_buffers
    if(present(allocate_3d_buffers_)) then
      allocate_3d_buffers = allocate_3d_buffers_
    else
      allocate_3d_buffers = .false.
    end if
    
    call compute_fft_sizes(nmod, n_x_grid, nx, &
       & mphit, mphiw3t, mrad_G, mrad_l)

    if (phi3d .or. apa3d .or. bpa3d .or. allocate_3d_buffers) then
      if(.not.allocated(out3dbuf_local_real)) then
        if(proc_subset(0,0,1,1,1)) then
          if(io_legacy) then
            allocate(out3dbuf_local_real(mphit,mrad_l,ns),stat=ierr)
          else
            allocate(out3dbuf_local_real(n_y_grid,nx,ns),stat=ierr)
          end if
          if (ierr /= 0) call gkw_abort('diagnostic :: out3dbuf')
        else
          allocate(out3dbuf_local_real(1,1,1),stat=ierr)
        end if
      end if
      if(.not.allocated(out3dbuf_global_real)) then
        if(root_processor) then
          if(io_legacy) then
            allocate(out3dbuf_global_real(mphit,mrad_G,n_s_grid),&
               & stat=ierr)
          else
            allocate(out3dbuf_global_real(n_y_grid,n_x_grid,n_s_grid),&
               & stat=ierr)
          end if
          if (ierr /= 0) call gkw_abort('diagnostic :: out3dbuf')
        else
          allocate(out3dbuf_global_real(1,1,1),&
             & stat=ierr)
        end if
      end if
    end if
    if(spc3d .or. apc3d .or. bpc3d .or. allocate_3d_buffers) then
      if(.not.allocated(out3dbuf_global_spec)) then
        if(root_processor) then
          allocate(out3dbuf_global_spec(nmod,n_x_grid,n_s_grid),&
             & stat=ierr)
          if (ierr /= 0) call gkw_abort('diagnostic :: out3dbuf')
        else
          allocate(out3dbuf_global_spec(1,1,1),&
             & stat=ierr)
        end if
      end if
    end if
    if(kykxs_phi .or. kykxs_apar .or. kykxs_bpar .or. allocate_3d_buffers) then
      if(.not.allocated(out3dbuf_global_cpx_spec)) then
        if(root_processor) then
          allocate(out3dbuf_global_cpx_spec(nmod,n_x_grid,n_s_grid),&
             & stat=ierr)
          if (ierr /= 0) call gkw_abort('diagnostic :: out3dbuf')
        else
          allocate(out3dbuf_global_cpx_spec(1,1,1),&
             & stat=ierr)
        end if
      end if
    end if

    ! if this is called by diagnos_mode_struct, then there is nothing
    ! else to do - the rest is already allocated.
    if(allocate_3d_buffers) return
    
    ! allocate the imod_field arrays:
    ! At first count all elements from the arrays no_transfer_to/from
    ! which are in the range of admissible indices.
    n = 0
    listloop: do i = 1, size(imod_field)
      if(imod_field(i) >= 1 .and. imod_field(i) <= nmod) then
        n = n + 1
      else if(imod_field(i) > nmod) then
        n = nmod
        imod_field(1) = imod_field(i)
        exit listloop
      end if
    end do listloop
    ! Now n contains the number of mode indices given by the input file.
    ! Allocate the list of requested mode numbers:
    allocate(imod_list(n),stat=ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate imod_list in diagnos_fields')
    end if
    ! fill this list: either will all mode numbers...
    if(imod_field(1) > nmod) then
      imod_list = (/ (i, i=1,nmod) /)
    else
      ! ... or just with a few, as given in the input namelist.
      n = 0
      do i = 1, size(imod_field)
        if(imod_field(i) >= 1 .and. imod_field(i) <= nmod) then
          n = n + 1
          imod_list(n) = imod_field(i)
        end if
      end do
    end if

    ! Arrays for FFT of a perpendicular slice
    if(io_legacy) then
      allocate(a(mphiw3t,mrad_l), stat = ierr)
    else
      allocate(a(nymod,nx), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('diagnos_fields: Could not allocate a in diagnostic')
    endif
    if(io_legacy) then
      allocate(ar(mphit,mrad_l), stat = ierr)
    else
      allocate(ar(n_y_grid,nx), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('diagnos_fields: Could not allocate ar in diagnostic')
    endif

    allocate(vort(nmod,nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: vort')

    ! array for MPI reduction
    ierr = 0 
    allocate(buffer(max(nx,nmod)), stat = ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: buffer')



  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    use mpiinterface, only : root_processor
    use io, only : close_lu, ascii_fmt, binary_fmt
    use control, only : normalized
    use control, only : nlphi, nlapar
    use mode, only : mode_box
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    
    if(allocated(a)) deallocate(a)
    if(allocated(ar)) deallocate(ar)
    if(allocated(vort)) deallocate(vort)
    if(allocated(field_par)) deallocate(field_par)
    if(allocated(field_par_global)) deallocate(field_par_global)
    if(allocated(buffer)) deallocate(buffer)
    if(allocated(out3dbuf_local_real)) then
      deallocate(out3dbuf_local_real)
    end if
    if(allocated(out3dbuf_global_real)) then
      deallocate(out3dbuf_global_real)
    end if
    if(allocated(out3dbuf_global_spec)) then
      deallocate(out3dbuf_global_spec)
    end if
    if(allocated(out3dbuf_global_cpx_spec)) then
      deallocate(out3dbuf_global_spec)
    end if

    if(.not. root_processor) return

    ! be nice and close all logical units
    if (mode_box) then 
      call close_lu(i_kyvort, ascii_fmt)
      call close_lu(i_kxvort, ascii_fmt)
    end if

    if (field_fsa_kyzero) then
       if(nlphi) call close_lu(i_rhtest, ascii_fmt)
       if(nlapar) call close_lu(i_apar_fsa_kyzero, ascii_fmt)
    end if

    if (.not. normalized) then
      call close_lu(i_kyspec, ascii_fmt)
      call close_lu(i_kxspec, ascii_fmt)

      if(nlapar) then
        call close_lu(i_kyspec_em, ascii_fmt)
        call close_lu(i_kxspec_em, ascii_fmt)
      end if
    end if

    if(xy_vort) call close_lu(lun_field_spec(VORTICITY_FIELD), binary_fmt)
    if(xy_phi) then
      call close_lu(lun_field_spec(PHI_FIELD), binary_fmt)
      call close_lu(lun_field_xy(PHI_FIELD), binary_fmt)
    end if
    if(xy_apar) then
      call close_lu(lun_field_spec(APAR_FIELD), binary_fmt)
      call close_lu(lun_field_xy(APAR_FIELD), binary_fmt)
    end if
    if(xy_bpar) then
      call close_lu(lun_field_spec(BPAR_FIELD), binary_fmt)
      call close_lu(lun_field_xy(BPAR_FIELD), binary_fmt)
    end if

    if(xs_phi) then
      call close_lu(lun_phi_xs, binary_fmt)
    end if

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
    use diagnos_generic, only : xy_slice_ipar, xy_estep
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use global, only : int2char_zeros
    integer, intent(in) :: number
    character (len=6), save :: number_suffix
    
    number_suffix=trim(int2char_zeros(number,6))
    
    if (.not.xy_estep) then
      call calc_largestep()
      
      call output_field_kykxs(PHI_FIELD,number,tag_range_start+6)
      call output_field_kykxs(APAR_FIELD,number,tag_range_start+7)
      call output_field_kykxs(BPAR_FIELD,number,tag_range_start+8)

      if (xy_apar) call output_field_perp2d(APAR_FIELD,xy_slice_ipar,number)
      if (xy_bpar) call output_field_perp2d(BPAR_FIELD,xy_slice_ipar,number)
      if (xy_phi) call output_field_perp2d(PHI_FIELD,xy_slice_ipar,number)
      if (xy_vort) call output_field_perp2d(VORTICITY_FIELD,1,number)

      if (xs_phi) then
        call output_field_xs(PHI_FIELD, number)
        !call output_field_xs(APAR_FIELD, number) ! uncomment, if you want it
        !call output_field_xs(BPAR_FIELD, number)
      end if
    end if
  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_largestep()
    call calc_vorticity
  end subroutine calc_largestep
  
  !--------------------------------------------------------------------
  !> This routine calculates the vorticity
  !> \f[
  !> \omega = \sum_{sp}\int d^3v \int ds n_{R,sp} Z_{sp} |f_{sp}|^2
  !> \f]
  !--------------------------------------------------------------------
  subroutine calc_vorticity()
    use grid,           only : nx, ns, nmu, nvpar, nsp
    use grid,           only : nmod, gx
    use dist,           only : fdisi
    use geom,           only : ints, bn
    use components,     only : de, signz
    use velocitygrid,   only : intmu, intvp
    use matdat,         only : get_f_from_g 
    use mpiinterface,   only : mpiallreduce_sum_inplace, number_of_processors
    use mpicomms,   only : COMM_X_EQ
    use global,         only : gkw_a_equal_b_accuracy

    ! integers for the loop over all grid points 
    integer :: imod, ix, i, j, k, is 

    ! ! Dummy variables 
    ! complex :: dum, dumes1, dumes2, dumem1, dumem2, dumbpar1, dumbpar2
    ! complex :: dumes_rs, dumem_rs, dumbpar_rs
    ! complex :: dumem1_cos
    complex :: fdis
    ! real    :: omega_times_time
    ! ! real phi2, apa2, bpa2
    ! real drift_x, drift_y, drift_z, dumnc, dumbs
    ! integer l, imomidx

    ! ! index for flux output array
    ! integer :: iflux

    ! Initialize to zero 
    vort = 0. 

    nmod1: do imod = 1, nmod 
      nx1: do ix = 1, nx 
        nsp1: do is = 1, nsp

          ! Integral over the velocity space 
          nmu3: do j = 1, nmu
            nvpar3: do k = 1, nvpar

              ! Do the average over the flux surface
              ns3: do i = 1, ns

                ! fdis is the distribution without A_par contribution  
                fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)

                ! in the implicit scheme fdis can be NaN for intvp = 0 
                if (gkw_a_equal_b_accuracy(intvp(i,j,k,is), 0.0)) fdis = 0.

                !Vorticity
                vort(imod,ix) = vort(imod,ix)+ signz(is)*de(ix,is)* &
                   & abs(fdis)**2 &
                   & *2.E0*bn(ix,i)*ints(i)*intvp(i,j,k,is)*intmu(j)
              end do ns3
            end do nvpar3
          end do nmu3
        end do nsp1
      end do nx1
    end do nmod1

    ! only when run on more than one processor (saves two copies)
    if (number_of_processors > 1) then

      ! FIXME This routine better should not sum a global array with a
      ! lot of zeros but gather local arrays. This should save
      ! communication and memory.
      call mpiallreduce_sum_inplace(vort,nmod,nx, COMM_X_EQ)

    end if

  end subroutine calc_vorticity

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output(file_count)
    use mode, only : mode_box
    use diagnos_generic, only : lphi_diagnostics, xy_estep
    use diagnos_generic, only : xy_slice_ipar, out3d_interval
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use control, only : nlphi, nlapar, normalized, itime_rst
    use control, only : nt_complete
    use global, only : int2char_zeros
    use mpiinterface, only : root_processor, gather_array, mpireduce_sum_inplace
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use grid, only : n_x_grid, nx, nmod
    use mpicomms, only : COMM_X_NE
    
    integer, intent(in) :: file_count

    character (len=6), save :: file_count_suffix
    real :: buf_global_x(n_x_grid), buf_y(nmod)

    ! Complex 3D fields on a (kx,ky,s) grid
    ! WARNING, data volume can be huge.
    if (xy_estep) then
      call output_field_kykxs(PHI_FIELD,file_count,tag_range_start)
      call output_field_kykxs(APAR_FIELD,file_count,tag_range_start+1)
      call output_field_kykxs(BPAR_FIELD,file_count,tag_range_start+2)
    end if

    if (.not. (mode_box .and. lphi_diagnostics)) return
    
    file_count_suffix=trim(int2char_zeros(file_count,6))

    call calc_largestep

    if (xy_estep) then
      if (xy_apar) call output_field_perp2d(APAR_FIELD,xy_slice_ipar, file_count)
      if (xy_bpar) call output_field_perp2d(BPAR_FIELD,xy_slice_ipar, file_count) 
      if (xy_phi) call output_field_perp2d(PHI_FIELD,xy_slice_ipar, file_count)
      if (xy_vort) call output_field_perp2d(VORTICITY_FIELD,1, file_count)

      if (xs_phi) then
        call output_field_xs(PHI_FIELD, file_count)
        !call output_field_xs(APAR_FIELD, file_count) ! uncomment, if you want it
        !call output_field_xs(BPAR_FIELD, file_count)
      end if
    end if

    ! fields spectra
    if (.not. normalized) then
       if(nlphi) then
          call output_fieldspec_kykx(PHI_FIELD,i_kyspec,i_kxspec)
          if(field_fsa_kyzero) call output_fsa_potential_zeromode(PHI_FIELD, &
             & i_rhtest)
       end if
       if (nlapar) then
          call output_fieldspec_kykx(APAR_FIELD,i_kyspec_em,i_kxspec_em)
          if(field_fsa_kyzero) call output_fsa_potential_zeromode(APAR_FIELD, &
             & i_apar_fsa_kyzero)
       end if
    end if

    ! If you want to add the 'Bpar' field output, you need to create files
    ! and luns then call
    !if (mode_box.and.nlbpar) then
    !  call output_fieldspec_kykx(BPAR_FIELD,i_kyspec_em2,i_kxspec_em2)
    !end if

    ! if(fluxes_spectra) then !legacy switch
      ! vorticity
      if (mode_box) then
        buf_y = sum(vort, 2)
        call mpireduce_sum_inplace(buf_y,shape(buf_y),COMM_X_NE)
        if(root_processor) then
          call append_chunk(i_kyvort, buf_y, xy_fmt, ascii_fmt)
        end if
        
        call gather_array(buf_global_x, n_x_grid, sum(vort, 1), nx, COMM_X_NE)
        if(root_processor) then
          call append_chunk(i_kxvort, buf_global_x, xy_fmt, ascii_fmt)
        end if
      end if
    ! end if

    ! Outputs 3D files for the whole flux tube every out3d_interval'th
    ! large timestep.
    ! WARNING, data volume can be huge.
    if(modulo(nt_complete+itime_rst,out3d_interval) == 0) then
      call output_field_3d(PHI_FIELD,file_count,tag_range_start+3)
      call output_field_3d(APAR_FIELD,file_count,tag_range_start+4)
      call output_field_3d(BPAR_FIELD,file_count,tag_range_start+5)
    end if
  end subroutine output


  !-----------------------------------------------------------------------------
  !> Generic 2D xy slice output of any field variable.
  !----------------------------------------------------------------------------
  subroutine output_field_perp2d(fieldname,isg,file_count)
    use grid,    only : nmod, nx, ls, proc_subset, ns
    use dist,    only : phi, apar, bpar
    use general, only : gkw_abort
    use diagnos_generic, only : xy_output_array, kyx_output_array
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD

    !> String that defines the type of the field.
    integer, intent(in) :: fieldname
    !> Global s index of the slice to output.
    integer, intent(in) :: isg
    integer, intent(in) :: file_count
    logical :: proc_is_in_slice

    ! an integer which contains
    !   * for procs containing the desired global s point: the
    !     corresponding local s index
    !   * for all other procs: the value 1, in order not to commit a memory
    !     access fault. (for those procs, the values passed to the
    !     *_output_array routines do not matter)
    integer :: is

    !> violates the letter of the law against pointers
    !> if you object you can define a new dummy field array
    !> or find a cleverer solution
    complex, pointer :: field(:,:,:)

    ! keep compiler quiet
    if (file_count > 0) continue

    if(fieldname == VORTICITY_FIELD) then
      ! output vok, including x gather
      proc_is_in_slice = proc_subset(0,1,1,1,1)
      call kyx_output_array( &
         & vort(1:nmod,1:nx),lun_field_spec(fieldname), &
         & proc_is_in_slice, .true.)

    else
      proc_is_in_slice = proc_subset(0,isg,1,1,1)
      if(proc_is_in_slice) then
        is = ls(isg)
      else
        is = 1
      end if

      ! fill the buffer
      select case (fieldname)
      case (PHI_FIELD)
        field => phi
      case(APAR_FIELD)
        field => apar
      case(BPAR_FIELD)
        field => bpar
      case default
        call gkw_abort('output_field_perp2d, wrong field labels') 
      end select

      ! output potential spectrum (modulus), including x gather
      call kyx_output_array( &
         & abs(field(1:nmod,1:nx,is)),lun_field_spec(fieldname), &
         & proc_is_in_slice, isg <= ns)

      ! output 2d slice of phi (in real space, including x gather)
      ! xy_output_array must be called from all x processes.
      call xy_output_array( &
         & field(1:nmod,1:nx,is), .true., lun_field_xy(fieldname), &
         & proc_is_in_slice, isg <= ns)

    end if

  end subroutine output_field_perp2d

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine output_field_xs(fieldname, file_count)
    use grid,         only : nx, n_s_grid, n_x_grid, ns, proc_subset
    use dist,         only : phi, apar, bpar
    use mpiinterface, only : gather_array
    use mpicomms,     only : COMM_S_NE, COMM_X_NE
    use mpiinterface, only : gather_array
    use io,           only : xy_fmt, binary_fmt
    use general,      only : gkw_abort
    use io, only : append_chunk
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD

    
    integer, intent(in) :: fieldname
    integer, intent(in) :: file_count
    integer ix, is, imode
    complex, dimension(nx,ns)              :: field_xs_local
    complex, dimension(n_x_grid, n_s_grid, &
       & size(imod_list)) :: field_xs_global

    !> violates the letter of the law against pointers
    !> if you object you can define a new dummy field array
    !> or find a cleverer solution
    complex, pointer :: field(:,:,:)

    integer :: lun

    ! keep compiler quiet
    if (file_count > 0) continue

    ! avoid some unnecessary gather-communication
    if (.not.proc_subset(0,0,1,1,1)) return

    ! fill the buffer
    select case (fieldname)
    case (PHI_FIELD)
      field => phi
      lun = lun_phi_xs
    case(APAR_FIELD)
      field => apar
      !lun = lun_apar_xs
      call gkw_abort('output_field_xs: apar is not implemented.')
    case(BPAR_FIELD)
      field => bpar
      !lun = lun_bpar_xs
      call gkw_abort('output_field_xs: bpar is not implemented.')
    case default
      call gkw_abort('output_field_xs: another field is requested.')
    end select
    
    do imode = 1, size(imod_list)
      ! fill the buffer
      do ix = 1, nx
        do is = 1, ns
          field_xs_local(ix,is) = field(imod_list(imode),ix,is)
        end do
      end do

      ! parallel_s and parallel_x gather
      call gather_array(field_xs_global(:,:,imode),n_x_grid,n_s_grid, &
         & field_xs_local,nx,ns, &
         & COMM_X_NE,COMM_S_NE,ALLGATHER=.true.)

    end do

    if (proc_subset(1,1,1,1,1)) then
      ! let the IO module deal with the real/imag parts:
      call append_chunk(lun, field_xs_global, xy_fmt, binary_fmt)
    end if
  end subroutine output_field_xs


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine output_fieldspec_kykx(fieldname,lun_ky,lun_kx) 

    use grid,         only : nmod, nx, n_s_grid, ls, n_x_grid, ns
    use grid,         only : proc_subset
    use dist,         only : phi, apar, bpar
    use mpicomms,     only : COMM_S_NE_X_NE
    use control,      only : spectral_radius
    use fft,          only : fourcol
    use geom,         only : ints
    use mpiinterface, only : gather_array, root_processor
    use mpiinterface, only : mpireduce_sum_inplace
    use mpicomms,     only : COMM_S_NE, COMM_X_NE, COMM_DUMMY
    use io,           only : append_chunk, xy_fmt, ascii_fmt
    use general,      only : gkw_abort
    use diagnos_generic, only : parseval_correction
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    integer, intent(in) :: fieldname
    integer, intent(in) :: lun_kx, lun_ky

    integer ::  imod, ix, is
    logical, save :: initialised = .false.
    complex, pointer :: fil(:,:,:)
    complex, allocatable, save :: ax(:,:), agx(:,:)
    integer :: ierr

    if (.not. proc_subset(0,0,1,1,1)) return
    
    if (.not.initialised) then 
      
      ierr = 0 
      allocate(ky_spec(nmod), stat = ierr)
      if (ierr /= 0) call gkw_abort('could not allocate ky_spec')

      ierr = 0 
      allocate(kx_spec(n_x_grid), stat = ierr)
      if (ierr /= 0) call gkw_abort('could not allocate kx_spec')
 
      if (spectral_radius) then
        ierr = 0 
        allocate(agx(n_x_grid,nmod), stat = ierr)
        if (ierr /= 0) call gkw_abort('could not allocate agx')
      end if
 
      if (.not. spectral_radius) then
        ierr = 0 
        allocate(ax(nx,nmod), stat = ierr)
        if (ierr /= 0) call gkw_abort('could not allocate ax')

        ierr = 0 
        allocate(agx(n_x_grid,nmod), stat = ierr)
        if (ierr /= 0) call gkw_abort('could not allocate agx')
      end if

      initialised = .true. 

    end if

    ! alternatively could just pass it to the function
    select case (fieldname)
    case (PHI_FIELD)
      fil => phi
    case (APAR_FIELD) 
      fil => apar
    case (BPAR_FIELD)
      fil => bpar
    case default
      call gkw_abort('unknown fieldname')
    end select

    kx_spec(:)=0.
    ky_spec(:)=0.

    ! kyspec
    do imod = 1, nmod 
      do is = 1, n_s_grid
        if(proc_subset(0,is,0,0,0)) then !For parallel_s
          do ix = 1, nx 
            ky_spec(imod) = ky_spec(imod) + abs(fil(imod,ix,ls(is)))**2 *ints(ls(is))
          end do
        end if
      end do
    end do

    ! PARALLEL_S or X
    call mpireduce_sum_inplace(ky_spec, shape(ky_spec), COMM_S_NE_X_NE)

    if (spectral_radius) then
      ! kxspec, buffer is reused
      do ix = 1, nx 
        do is = 1, n_s_grid 
          if(proc_subset(0,is,0,0,0)) then !For parallel_s
            do imod = 1, nmod
              kx_spec(ix) = kx_spec(ix) + parseval_correction(imod)*abs(fil(imod,ix,ls(is)))**2 *ints(ls(is))
            end do
          end if
        end do
      end do

      ! PARALLEL_S
      call mpireduce_sum_inplace(kx_spec, shape(kx_spec), COMM_S_NE)

    else
      do is = 1,ns
        ! get a slice
        do imod = 1, nmod
          
          do ix = 1, nx
            ax(ix, imod) = fil(imod,ix,is)
          end do
        end do

        ! gather in x direction
        call gather_array(agx, n_x_grid, nmod, &
           & ax, nx, nmod, &
           & COMM_X_NE, COMM_DUMMY)
        
        if(proc_subset(1,0,1,1,1)) then

          ! fouriertrafo in x direction (= the first dimension of the array)
          call fourcol(agx,-1)

          ! integrate the |phi|^2 locally
          do imod = 1, nmod
            do ix = 1, n_x_grid
              kx_spec(ix) = kx_spec(ix) + abs(agx(ix, imod))**2 * &
                 & ints(is) * parseval_correction(imod)
            end do
          end do
        end if
      end do

      if(proc_subset(1,0,1,1,1)) then
        ! finish off FSA by mpireduce_sum in s direction
        call mpireduce_sum_inplace(kx_spec, shape(kx_spec), COMM_S_NE)
      end if

    end if

    if (root_processor) then
      call append_chunk(lun_ky, ky_spec, xy_fmt, ascii_fmt)
      call append_chunk(lun_kx, kx_spec, xy_fmt, ascii_fmt)
    end if
      
  end subroutine output_fieldspec_kykx

  
  !----------------------------------------------------------------------
  !> Produce the radial field spectrum used in the Rosenbluth-Hinton Test
  !----------------------------------------------------------------------
  subroutine output_fsa_potential_zeromode(fieldname, lun)
    use mpiinterface, only : mpireduce_sum_inplace, gather_array, root_processor
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use mpicomms, only : COMM_S_NE, COMM_X_NE
    use grid, only : proc_subset, nx, n_x_grid, ns
    use geom, only : ints
    use mode, only : iyzero
    use general, only : gkw_abort
    use dist, only : phi, apar, bpar
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD

    integer, intent(in) :: fieldname
    integer, intent(in) :: lun
    !> a line buffers to hold a profile or a spectrum; are only used on a few processes
    complex :: x_array_local(1:nx), x_array_global(1:n_x_grid)
    complex, pointer :: fil(:,:,:)

    integer :: is, ix
    
    if(proc_subset(0,0,1,1,1)) then

      select case (fieldname)
      case (PHI_FIELD)
        fil => phi
      case (APAR_FIELD) 
        fil => apar
      case (BPAR_FIELD)
        fil => bpar
      case default
        call gkw_abort('unknown fieldname')
      end select

      x_array_local = 0.0e0
      do is = 1, ns
        do ix = 1, nx
          x_array_local(ix) = x_array_local(ix) + fil(iyzero,ix,is)*ints(is)
        end do
      end do
      call mpireduce_sum_inplace(x_array_local,shape(x_array_local),COMM_S_NE)
      if(proc_subset(0,1,1,1,1)) call gather_array(x_array_global, n_x_grid, &
         & x_array_local, nx, COMM_X_NE)
      ! do the fourier trafo in x when analysing the data, if you want it
      if(root_processor) call append_chunk(lun, x_array_global, xy_fmt, ascii_fmt)
    end if

  end subroutine output_fsa_potential_zeromode


  !----------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------
  subroutine output_field_3d(fieldname, file_count, tag)
    use control, only : spectral_radius, time, io_legacy
    use io, only : binary_fmt, not_avail, phys_unit_key, description_key
    use io, only : attach_metadata, comments_key, lu_exists
    use global, only : int2char_zeros
    use general, only : gkw_abort
    use diagnos_generic, only : phi3d, spc3d, apa3d, apc3d, bpa3d, bpc3d
    use diagnos_generic, only : attach_metadata_grid, rad_gridname
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use mpiinterface, only : root_processor

    integer, intent(in) :: fieldname
    integer, intent(in) :: file_count
    integer, intent(in) :: tag
    character (len=13) :: luname_spec, luname_real
    character (len=30) :: desc_snippet
    character (len=50) :: name
    logical :: output_spec_flag, output_real_flag

    if(spectral_radius) then
      desc_snippet = 'and the 3rd (radial)'
    else
      desc_snippet = ''
    end if

    select case(fieldname)
    case(PHI_FIELD)
      luname_spec = 'Spc3d'//trim(int2char_zeros(file_count,8))
      luname_real = 'Poten'//trim(int2char_zeros(file_count,8))
      name = 'electrostatic potential'
      output_spec_flag = spc3d
      output_real_flag = phi3d
    case(APAR_FIELD)
      luname_spec = 'Apc3d'//trim(int2char_zeros(file_count,8))
      luname_real = 'Apara'//trim(int2char_zeros(file_count,8))
      name = 'parallel component of the vector potential'
      output_spec_flag = apc3d
      output_real_flag = apa3d
    case(BPAR_FIELD)
      luname_spec = 'Bpc3d'//trim(int2char_zeros(file_count,8))
      luname_real = 'Bpara'//trim(int2char_zeros(file_count,8))
      name = 'parallel component of the magn. field'
      output_spec_flag = bpc3d
      output_real_flag = bpa3d
    case default
      call gkw_abort('Wrong field selected in binary output')
    end select


    call output_field_3d_parallel(fieldname, luname_spec, luname_real, &
       & output_spec_flag, output_real_flag, tag)

    ! Workaround for issue #256
    return

    !FIXME correct x grid metadata for spectral runs!
    
    !FIXME metadata is not nice in the gkwdata.meta file, would have
    !to split the next block into 2 routines and call just after output or so.
    if(output_spec_flag  .and. lu_exists(luname_spec, &
       & 'diagnostic/diagnos_fields', binary_fmt) &
       & .and. root_processor) then
      call attach_metadata_grid(luname_spec, &
         & 'diagnostic/diagnos_fields', 'krho', rad_gridname, 'sgrid', &
         & binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', description_key, &
         & trim(name)//' on the complete 3d computational spatial&
         & domain, in spectral representation with respect to the 1st &
         & (binormal) '//trim(desc_snippet)//' dimension.', binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', comments_key, not_avail, binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', 'time', time, binary_fmt)
    end if

    if(output_real_flag .and. &
       & lu_exists(luname_real, 'diagnostic/diagnos_fields', binary_fmt) &
       & .and. root_processor) then
      ! may not exist, as serial 3d real space output is not yet implemented
      if(io_legacy) then
        call attach_metadata_grid(luname_real, &
           & 'diagnostic/diagnos_fields', 'yphi', 'xphi', 'sgrid', binary_fmt)
      else
        ! FIXME what is the grid which is equivalent to yphi and xphi,
        ! but is just n_y_grid and n_x_grid elements long?
        call attach_metadata_grid(luname_real, &
           & 'diagnostic/diagnos_fields', 'n_y_grid', 'n_x_grid', 'sgrid', binary_fmt)
      end if
      call attach_metadata(luname_real, &
         & 'diagnostic/diagnos_fields', phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata(luname_real, &
         & 'diagnostic/diagnos_fields', description_key, &
         & trim(name)//' on the complete 3d computational spatial&
         & domain, in position space representation',&
         & binary_fmt)
      call attach_metadata(luname_real, &
         & 'diagnostic/diagnos_fields', comments_key, not_avail, binary_fmt)
      call attach_metadata(luname_real, &
         & 'diagnostic/diagnos_fields', 'time', time, binary_fmt)
    end if
  end subroutine output_field_3d


  !----------------------------------------------------------------------
  !> Outputs the full three dimensional field by calling a routine from the
  !> IO module. which either (i) uses MPI-IO and the subarray datatype to
  !> output raw binary data, or (ii) uses the subarray datatype to gather
  !> the global field to one process and output the data serially.
  !----------------------------------------------------------------------
  subroutine output_field_3d_parallel(fieldname, luname_spec, luname_real, &
     & output_spec_flag, output_real_flag, tag)

    use io,               only : mpi_output_array, xy_fmt, binary_fmt
    use grid,             only : nmod,nx,ns, proc_subset, n_y_grid
    use mpicomms,         only : COMM_S_NE_X_NE
    use dist,             only : phi,apar,bpar
    use control,          only : io_legacy
    use non_linear_terms, only : jind
    use general,          only : gkw_abort
    use diagnos_generic,  only : zonal_scale_3d
    use diagnos_generic,  only : four2real_2D
    use fft, only : FFT_INVERSE
    use grid,             only : jind_flexible
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use diagnos_generic, only : mpi_dtype_real_spec_yxs, mpi_dtype_real_yxs
    
    integer, intent(in) :: fieldname
    character (len=*), intent(in) :: luname_spec, luname_real
    logical, intent(in) :: output_spec_flag, output_real_flag
    integer, intent(in) :: tag

    integer :: i,j,imod,ix,ipar

    !> violates the rule against pointers, but convenient
    complex, pointer :: buf3d_local_cmplx(:,:,:)
    logical, parameter :: root_is_in_comm = .true.
    logical :: proc_is_in_slab = .false.

    if(.not. output_spec_flag .and. .not. output_real_flag) return
    
    if(proc_subset(0,0,1,1,1)) then
      proc_is_in_slab = .true.

      !Obtain the local field
      select case(fieldname)
      case(PHI_FIELD)
        buf3d_local_cmplx => phi
      case(APAR_FIELD)
        buf3d_local_cmplx => apar
      case(BPAR_FIELD)
        buf3d_local_cmplx => bpar
      case default
        call gkw_abort('Wrong field selected in MPI-IO binary output')
      end select
    else
      ! just anything, so that it can be passed to the subroutine below
      buf3d_local_cmplx => phi
    end if


    if(output_spec_flag) then
      !Write the spectral field
      call mpi_output_array(luname_spec, 'diagnostic/diagnos_fields', &
         & abs(buf3d_local_cmplx(1:nmod, &
         & 1:nx, &
         & 1:ns)), mpi_dtype_real_spec_yxs, &
         & out3dbuf_global_spec, &
         & COMM_S_NE_X_NE, xy_fmt, binary_fmt, proc_is_in_slab, root_is_in_comm,&
         & tag)
    end if

    if(.not. output_real_flag) return

    if(proc_subset(0,0,1,1,1)) then
      ! fourier transform every perpendicular slice of the local array,
      ! to obtain the real space representation of the local array
      local_s_grid : do ipar = 1, ns
        if(io_legacy) then
        ! fill the array for the potential 
          a = (0.,0.)
          do ix = 1, nx;
            do imod = 1, nmod
              a(imod,jind(ix)) = buf3d_local_cmplx(imod,ix,ipar)
            end do
          end do
        else
          if(nmod == 1) a(1,:) = (0.,0.)
          do ix = 1, nx;
            do imod = 1, nmod
              if(nmod == 1) then
                a(imod+1,jind_flexible(nx, nx, ix)) = buf3d_local_cmplx(imod,ix,ipar)
              else
                a(imod,jind_flexible(nx, nx, ix)) = buf3d_local_cmplx(imod,ix,ipar)
              end if
            end do
          end do
        end if

        !Cheat on the zonal flow appearance
        a(1,:)=zonal_scale_3d*a(1,:)

        !Do the inverse FFT
        call four2real_2D(ar,a,FFT_INVERSE)
        
        ! copy slice into write buffer
        if(io_legacy) then
          do j=1,mrad_l
            do i=1,mphit
              out3dbuf_local_real(i,j,ipar) = ar(i,j)
            end do
          end do
        else
          do j=1,nx
            do i=1,n_y_grid
              out3dbuf_local_real(i,j,ipar) = ar(i,j)
            end do
          end do
        end if
      end do local_s_grid
    end if
    ! Write the real space field
    call mpi_output_array(luname_real, 'diagnostic/diagnos_fields', &
       & out3dbuf_local_real, mpi_dtype_real_yxs, &
       & out3dbuf_global_real, &
       & COMM_S_NE_X_NE, xy_fmt, binary_fmt, proc_is_in_slab, root_is_in_comm,&
       & tag)

  end subroutine output_field_3d_parallel

  !----------------------------------------------------------------------
  !> Output the 3D fields for spectral runs on a ky,kx,s grid
  !----------------------------------------------------------------------
  subroutine output_field_kykxs(fieldname, file_count, tag)
    use control, only : time
    use io, only : binary_fmt, not_avail, phys_unit_key, description_key
    use io, only : attach_metadata, comments_key, lu_exists
    use global, only : int2char_zeros
    use general, only : gkw_abort
    use diagnos_generic, only : attach_metadata_grid
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    
    integer, intent(in) :: fieldname
    integer, intent(in) :: file_count
    integer, intent(in) :: tag
    character (len=17) :: luname_spec
    character (len=50) :: name
    logical :: output_spec_flag


    select case(fieldname)
    case(PHI_FIELD)
      luname_spec = 'Phi_kykxs'//trim(int2char_zeros(file_count,8))
      name = 'electrostatic potential'
      output_spec_flag = kykxs_phi
    case(APAR_FIELD)
      luname_spec = 'Apa_kykxs'//trim(int2char_zeros(file_count,8))
      name = 'parallel component of the vector potential'
      output_spec_flag = kykxs_apar
    case(BPAR_FIELD)
      luname_spec = 'Bpa_kykxs'//trim(int2char_zeros(file_count,8))
      name = 'parallel component of the magn. field'
      output_spec_flag = kykxs_bpar
    case default
      call gkw_abort('Wrong field selected in binary output')
    end select

    call output_field_kykxs_parallel(fieldname, luname_spec, output_spec_flag, &
       & tag)

    if(output_spec_flag  .and. lu_exists(luname_spec, &
       & 'diagnostic/diagnos_fields', binary_fmt)) then
      call attach_metadata_grid(luname_spec, &
         & 'diagnostic/diagnos_fields', 'krho', 'kxrh', 'sgrid', &
         & binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', phys_unit_key, not_avail, &
         & binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', description_key, &
         & trim(name)//' on the complete 3d computational spatial&
         & domain, in spectral representation with respect to the &
         & binormal, radial and parallel dimensions', binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', comments_key, not_avail, binary_fmt)
      call attach_metadata(luname_spec, &
         & 'diagnostic/diagnos_fields', 'time', time, binary_fmt)
    end if

  end subroutine output_field_kykxs


  !----------------------------------------------------------------------
  !> Outputs the full three dimensional field by calling a routine from the
  !> IO module. which either (i) uses MPI-IO and the subarray datatype to
  !> output raw binary data, or (ii) uses the subarray datatype to gather
  !> the global field to one process and output the data serially.
  !----------------------------------------------------------------------
  subroutine output_field_kykxs_parallel(fieldname, luname_spec, &
     & output_spec_flag, tag)
    use io,               only : mpi_output_array, xy_fmt, binary_fmt
    use grid,             only : nmod,nx,ns, proc_subset
    use mpicomms,         only : COMM_S_NE_X_NE
    use dist,             only : phi,apar,bpar
    use general,          only : gkw_abort
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use diagnos_generic, only : mpi_dtype_real_spec_yxs
    
    integer, intent(in) :: fieldname
    character (len=*), intent(in) :: luname_spec
    logical, intent(in) :: output_spec_flag
    integer, intent(in) :: tag

    !> violates the rule against pointers, but convenient
    complex, pointer :: buf3d_local_cmplx(:,:,:)
    logical, parameter :: root_is_in_comm = .true.
    logical :: proc_is_in_slab = .false.

    if(.not. output_spec_flag) return
    
    if(proc_subset(0,0,1,1,1)) then
      proc_is_in_slab = .true.

      !Obtain the local field
      select case(fieldname)
      case(PHI_FIELD)
        buf3d_local_cmplx => phi
      case(APAR_FIELD)
        buf3d_local_cmplx => apar
      case(BPAR_FIELD)
        buf3d_local_cmplx => bpar
      case default
        call gkw_abort('Wrong field selected in MPI-IO binary output')
      end select
    else
      ! just anything, so that it can be passed to the subroutine below
      buf3d_local_cmplx => phi
    end if


    if(output_spec_flag) then
      !Write the spectral field
      !Warning: the mpi_dtype needs to be of REAL type as ultimately
      !two files of real are written 
      call mpi_output_array(luname_spec, 'diagnostic/diagnos_fields', &
         & buf3d_local_cmplx(1:nmod, &
         & 1:nx, &
         & 1:ns), mpi_dtype_real_spec_yxs, &
         & out3dbuf_global_cpx_spec, &
         & COMM_S_NE_X_NE, xy_fmt, binary_fmt, proc_is_in_slab, root_is_in_comm,&
         & tag)

    end if


  end subroutine output_field_kykxs_parallel

  !----------------------------------------------------------------------
  !> gather the 3D global field to root of comm_cart
  !----------------------------------------------------------------------
  subroutine gather_field_3d_serial(fieldname, field_buf_global, &
     & get_real_space_field, rotate)
    use grid,           only : nmod, ns, nx, proc_subset
    use dist,           only : phi, apar, bpar
    use mpiinterface,   only : gather_array
    use mpicomms, only : COMM_S_NE_X_NE
    use general,          only : gkw_abort, gkw_warn
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use diagnos_generic, only : mpi_dtype_cpx_yxs, mpi_dtype_cpx_spec_yxs
    
    integer, intent(in) :: fieldname
    complex, dimension(:,:,:), intent(out) :: field_buf_global
    logical, intent(in) :: get_real_space_field
    complex, optional, intent(in) :: rotate(nmod)

    complex, dimension(nmod, nx, ns) :: field_buf
    integer :: mpi_dtype
    integer :: imod

    logical, parameter :: to_root_of_commcart = .true.
    logical, parameter :: root_is_in_comm = .true.

    proc_in_slab: if(proc_subset(0,0,1,1,1)) then
      select case(fieldname)
      case(PHI_FIELD)
        field_buf = phi(:,1:nx,1:ns)
      case(APAR_FIELD)
        field_buf = apar(:,1:nx,1:ns)
      case(BPAR_FIELD)
        field_buf = bpar(:,1:nx,1:ns)
      case default
        call gkw_abort('Wrong moment called in gather_field_3d_serial')
      end select

      rotate_is_present: if(present(rotate)) then
        do imod = 1, nmod
          ! for all other fields and moments, the rotation factor is now
          ! applied, as it is given as an argument to this routine:
          ! Rotate all values in the complex plane relative to maximum potential
          ! defined to be (1.,0.).  Effectively, normalise out frequency.
          field_buf(imod, :, :) = field_buf(imod, :, :) * rotate(imod)
        end do
      end if rotate_is_present

      if(get_real_space_field) then
        ! first, perform an FFT on the local array to obtain the real
        !space repr...  TODO

        mpi_dtype = mpi_dtype_cpx_yxs
        ! THIS IS NOT YET IMPLEMENTED - and not used at the moment
        call gkw_warn('not implemented: a real space field is requested from&
           & gather_field_3d_serial.')

      else

        mpi_dtype = mpi_dtype_cpx_spec_yxs

      end if

      ! (In order to create a subarray datatype, every process must
      ! know the shape of the global field buffer, so hardcode that
      ! shape.)  gather the data to the root process(es) of the
      ! communicator(s) COMM_S_NE_X_NE:
      call gather_array(field_buf_global, &
         & field_buf, &
         & mpi_dtype, COMM_S_NE_X_NE, to_root_of_commcart, root_is_in_comm)

    end if proc_in_slab

  end subroutine gather_field_3d_serial

end module diagnos_fields

