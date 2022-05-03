!------------------------------------------------------------------------------
!> Output moments (density, flow, temperature perturbations) in 3D and in
!> 2D perpendicular slices (spectral and real).
!------------------------------------------------------------------------------
module diagnos_moments
  use diagnos_generic, only : DENSITY_MOMENT, TEMPERATURE_MOMENT
  use diagnos_generic, only : PAR_TEMPERATURE_MOMENT, PERP_TEMPERATURE_MOMENT
  use diagnos_generic, only : PAR_TEMPERATURE_GA_MOMENT, PERP_TEMPERATURE_GA_MOMENT
  use diagnos_generic, only : PERP_TEMPERATURE_J1_MOMENT
  use diagnos_generic, only : CURRENT_MOMENT, CURRENT_SQ_MOMENT, VPERP_SQ_MOMENT
  use diagnos_generic, only : CURRENT_GA_MOMENT
  use diagnos_generic, only : PASSING_DENSITY_MOMENT, TRAPPED_DENSITY_MOMENT
  use diagnos_generic, only : DENSITY_GA_MOMENT, DENSITY_POLAR_MOMENT
  use diagnos_generic, only : PHI_GA_FM_MOMENT, PHI_GA_DELTAF_MOMENT
  use diagnos_generic, only : QPAR_MOMENT, M12_MOMENT, M24_MOMENT
  use diagnos_generic, only : QPAR_GA_MOMENT, M12_GA_MOMENT, M24_GA_MOMENT

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output

  public :: gather_moment_4d_serial
  public :: get_4d_moment

  !> the range of tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl

  !> For output of moments spectra
  real, save, allocatable, dimension(:,:) :: moment_spec

  !> logical unit numbers
  integer, save :: i_denspec = -1, i_enespec = -1

  !> logical unit numbers for 2d output of various moments
  !> Note: please remove the second dimension if io_legacy is ever removed
  integer, save, allocatable, dimension(:,:) :: lun_moment_spec
  integer, save, allocatable, dimension(:,:) :: lun_moment_xy
  integer, save, allocatable, dimension(:,:) :: lun_moment_xs

  !> complex slice in xy for mpi reduction
  complex, save, allocatable, dimension(:,:) :: c_xy
  complex, save, allocatable, dimension(:,:,:) :: c_xysp
  
  !> complex slice in xs for mpi reduction
  complex, save, allocatable, dimension(:,:) :: mom_xs_local
  
  !> A buffer for communication. buf(nmod,nx,number_of_species)
  real, save, allocatable, dimension(:,:,:) :: buf(:,:,:)
  
  !> (testing) buffer for 3D outputs
  real, save, allocatable :: out3dbuf_local_real(:,:,:)
  real, save, allocatable :: out3dbuf_global_real(:,:,:)

  !> Private FFT arrays.
  complex, save, allocatable :: a(:,:)
  real, save, allocatable :: ar(:,:)
  
  !> Logicals to select ky,kx,s moments output
  !> Works for mode_box=T and F
  logical, save, public :: kykxs_moments
  logical, save, public :: kykxs_j0_moments
  logical, save, public :: kykxs_j1_moments

  !> swtiches for xs_output
  logical, save, public :: xs_kyzero_dens
  logical, save, public :: xs_kyzero_ene
  logical, save, public :: xs_kyzero_ene_par
  logical, save, public :: xs_kyzero_ene_perp
  logical, save, public :: xs_kyzero_current
  logical, save, public :: xs_kyzero_current2
  logical, save, public :: xs_kyzero_phi_ga_fm
  logical, save, public :: xs_kyzero_phi_ga_deltaf

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
 
     kykxs_moments = .false.
     kykxs_j0_moments = .false.
     kykxs_j1_moments = .false.
     
    xs_kyzero_dens = .false.
    xs_kyzero_ene = .false.
    xs_kyzero_ene_par = .false.
    xs_kyzero_ene_perp = .false.
    xs_kyzero_current = .false.
    xs_kyzero_current2 = .false.
    xs_kyzero_phi_ga_fm = .false.
    xs_kyzero_phi_ga_deltaf = .false.
  
  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(kykxs_moments,1)
    call mpibcast(kykxs_j0_moments,1)
    call mpibcast(kykxs_j1_moments,1)
    
    call mpibcast(xs_kyzero_dens,1)
    call mpibcast(xs_kyzero_ene,1)
    call mpibcast(xs_kyzero_ene_par,1)
    call mpibcast(xs_kyzero_ene_perp,1)
    call mpibcast(xs_kyzero_current,1)
    call mpibcast(xs_kyzero_current2,1)
    call mpibcast(xs_kyzero_phi_ga_fm,1)
    call mpibcast(xs_kyzero_phi_ga_deltaf,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : spectral_radius, io_format
    use general, only : gkw_warn
    use io, only : ascii_fmt

    if ((.not. spectral_radius) .and. (kykxs_moments .or. kykxs_j0_moments .or. kykxs_j1_moments)) then
      call gkw_warn('kykxs_xxx diagnostics only implemented for spectral_radius=.true.')       
      kykxs_moments=.false.
      kykxs_j0_moments=.false.
      kykxs_j1_moments=.false.
    end if

    if ((io_format == ascii_fmt).and. (kykxs_moments .or. kykxs_j0_moments .or. kykxs_j1_moments)) then
      call gkw_warn('kykxs_xxx diagnostics not implemented for ascii io_format')       
      kykxs_moments=.false.
      kykxs_j0_moments=.false.
      kykxs_j1_moments=.false.
    end if
    

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use io, only : open_real_lu, ascii_fmt, binary_fmt, open_complex_lu
    use grid, only : nmod, number_of_species, n_x_grid
    use grid, only : n_y_grid, n_s_grid
    use global, only : int2char_zeros
    use mpiinterface, only : root_processor, register_tag_range
    use control, only : io_legacy
    use diagnos_generic, only : mphit, mrad_G
    use diagnos_generic, only : xy_current, xy_current2, xy_dens, xy_temp
    use diagnos_generic, only : xy_spec, lphi_diagnostics, xy_estep
    use mode, only : mode_box
    logical, intent(inout) :: requirements(:,:)

    integer :: is
    integer :: xyslice_shape(2)
    integer :: xsslice_shape(2)
    
    ! To keep the compiler quiet, as the array is not used so far (in this routine).
    if (.false.) write(*,*) requirements

    ! tags for den3d/ene3d (2*number_of_species and for kykxs diagnostics (7+7+1)*number_of_species
    call register_tag_range(2 * number_of_species &
       & + 15 * number_of_species, &
       & tag_range_start, tag_range_end_inkl)

    if(io_legacy) then
      xyslice_shape(1) = mphit
      xyslice_shape(2) = mrad_G
      
      xsslice_shape(1) = mrad_G
      xsslice_shape(2) = n_s_grid
      
    else
      xyslice_shape(1) = n_y_grid
      xyslice_shape(2) = n_x_grid
      
      xsslice_shape(1) = n_x_grid
      xsslice_shape(2) = n_s_grid
      
    end if

    if(root_processor) then
      if(io_legacy) then
        call open_real_lu('den_spectra.dat', 'diagnostic/diagnos_moments', &
           & (/ nmod*number_of_species /), &
           & ascii_fmt, i_denspec)
        call open_real_lu('ene_spectra.dat', 'diagnostic/diagnos_moments', &
           & (/ nmod*number_of_species /), &
           & ascii_fmt, i_enespec)
      else
        call open_real_lu('den_spectra', 'diagnostic/diagnos_moments', &
           & (/ nmod, number_of_species /), &
           & ascii_fmt, i_denspec)
        call open_real_lu('ene_spectra', 'diagnostic/diagnos_moments', &
           & (/ nmod, number_of_species /), &
           & ascii_fmt, i_enespec)
      end if
    end if
 
   
    if (.not.xy_estep .or. .not. (mode_box .and. lphi_diagnostics)) return
    
    if(io_legacy) then
      do is = 1, number_of_species
        if(xy_dens) then
          if(xy_spec) then
            call open_real_lu('des'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
               & (/ nmod, n_x_grid /), &
               & binary_fmt, lun_moment_spec(DENSITY_MOMENT,is))
          end if

          call open_real_lu('den'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xyslice_shape, &
             & binary_fmt, lun_moment_xy(DENSITY_MOMENT,is))
        end if
        if(xy_temp) then
          if(xy_spec) then
            call open_real_lu('ens'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
               & (/ nmod, n_x_grid /), &
               & binary_fmt, lun_moment_spec(TEMPERATURE_MOMENT,is))
          end if
          call open_real_lu('ene'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xyslice_shape, &
             & binary_fmt, lun_moment_xy(TEMPERATURE_MOMENT,is))

        end if
        if(xy_current) then
          if(xy_spec) then
            call open_real_lu('pas'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
               & (/ nmod, n_x_grid /), &
               & binary_fmt, lun_moment_spec(CURRENT_MOMENT,is))
          end if
          call open_real_lu('pac'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xyslice_shape, &
             & binary_fmt, lun_moment_xy(CURRENT_MOMENT,is))

        end if
        if(xy_current2) then
          if(xy_spec) then
            call open_real_lu('p2s'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
               & (/ nmod, n_x_grid /), &
               & binary_fmt, lun_moment_spec(CURRENT_SQ_MOMENT,is))
          end if
          call open_real_lu('p2c'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xyslice_shape, &
             & binary_fmt, lun_moment_xy(CURRENT_SQ_MOMENT,is))
        end if
      end do
      
      
    else
      if(xy_dens) then
        if(xy_spec) then
          call open_real_lu('dens_xky', 'diagnostic/diagnos_moments', &
             & (/ nmod, n_x_grid, number_of_species /), &
             & binary_fmt, lun_moment_spec(DENSITY_MOMENT,1))
        end if

        call open_real_lu('dens_xy', 'diagnostic/diagnos_moments', &
           & (/ n_y_grid, n_x_grid, number_of_species /), &
           & binary_fmt, lun_moment_xy(DENSITY_MOMENT,1))

      end if
      if(xy_temp) then
        if(xy_spec) then
          call open_real_lu('temp_xky', 'diagnostic/diagnos_moments', &
             & (/ nmod, n_x_grid, number_of_species /), &
             & binary_fmt, lun_moment_spec(TEMPERATURE_MOMENT,1))
        end if
        call open_real_lu('temp_xy', 'diagnostic/diagnos_moments', &
           & (/ n_y_grid, n_x_grid, number_of_species /), &
           & binary_fmt, lun_moment_xy(TEMPERATURE_MOMENT,1))
      end if
      if(xy_current) then
        if(xy_spec) then
          call open_real_lu('current_xky', 'diagnostic/diagnos_moments', &
             & (/ nmod, n_x_grid, number_of_species /), &
             & binary_fmt, lun_moment_spec(CURRENT_MOMENT,1))
        end if
        call open_real_lu('current_xy', 'diagnostic/diagnos_moments', &
           & (/ n_y_grid, n_x_grid, number_of_species /), &
           & binary_fmt, lun_moment_xy(CURRENT_MOMENT,1))
      end if
      if(xy_current2) then
        if(xy_spec) then
          call open_real_lu('currentsq_xky', 'diagnostic/diagnos_moments', &
             & (/ nmod, n_x_grid, number_of_species /), &
             & binary_fmt, lun_moment_spec(CURRENT_SQ_MOMENT,1))
        end if
        call open_real_lu('currentsq_xy', 'diagnostic/diagnos_moments', &
           & (/ n_y_grid, n_x_grid, number_of_species /), &
           & binary_fmt, lun_moment_xy(CURRENT_SQ_MOMENT,1))
      end if
    end if
    
    
    ! xs-output
      do is = 1, number_of_species
        if(xs_kyzero_dens) then
          call open_real_lu('dens_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(DENSITY_MOMENT,is))
        end if
        if(xs_kyzero_ene) then
          call open_real_lu('ene_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(TEMPERATURE_MOMENT,is))
        end if
        if(xs_kyzero_ene_par) then
          call open_real_lu('ene_par_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(PAR_TEMPERATURE_MOMENT,is))
        end if
        if(xs_kyzero_ene_perp) then
          call open_real_lu('ene_perp_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(PERP_TEMPERATURE_MOMENT,is))
        end if
        if(xs_kyzero_current) then
          call open_real_lu('current_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(CURRENT_MOMENT,is))
        end if
        if(xs_kyzero_current2) then
          call open_real_lu('currentsq_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(CURRENT_SQ_MOMENT,is))
        end if
        if(xs_kyzero_phi_ga_fm) then
          call open_real_lu('phi_ga_fm_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(PHI_GA_FM_MOMENT,is))
        end if
        if(xs_kyzero_phi_ga_deltaf) then
          call open_real_lu('phi_ga_deltaf_kyzero_xs'//trim(int2char_zeros(is,2)), 'diagnostic/diagnos_moments', &
             & xsslice_shape, &
             & binary_fmt, lun_moment_xs(PHI_GA_DELTAF_MOMENT,is))
        end if
      end do
    
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem(allocate_3d_buffers_)
    use general, only : gkw_abort
    use grid, only : nmod,nymod,nx,number_of_species,n_x_grid,ns, n_s_grid, n_y_grid
    use grid, only : proc_subset
    use diagnos_generic, only : den3d, ene3d, xy_dens, xy_temp
    use diagnos_generic, only : xy_current, xy_current2
    use diagnos_generic, only : mphiw3t, mphit, mrad_l, mrad_G
    use diagnos_generic, only : allocate_spec_cmpx_buffers
    use diagnos_generic, only : N_MOMENTS
    use mpiinterface, only : root_processor
    use control, only : io_legacy
    !> by calling with this optional argument .true., the
    !> diagnos_mode_struct diagnostic can make sure that the buffers necessary
    !> for the gather operation are allocated.
    logical, optional, intent(in) :: allocate_3d_buffers_
    integer :: ierr
    logical :: allocate_3d_buffers
    if(present(allocate_3d_buffers_)) then
      allocate_3d_buffers = allocate_3d_buffers_
    else
      allocate_3d_buffers = .false.
    end if

    if (den3d .or. ene3d .or. allocate_3d_buffers) then
      if(.not.allocated(out3dbuf_local_real)) then
        if(proc_subset(0,0,1,1,0))then
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
      if (.not.allocated(out3dbuf_global_real)) then
        if(root_processor) then
          if(io_legacy) then
            allocate(out3dbuf_global_real(mphit,mrad_G,n_s_grid), &
               & stat=ierr)
          else
            allocate(out3dbuf_global_real(n_y_grid,n_x_grid,n_s_grid),&
               & stat=ierr)
          end if
          if (ierr /= 0) call gkw_abort('diagnostic :: out3dbuf')
        else
          allocate(out3dbuf_global_real(1,1,1),stat=ierr)
        end if
      end if
    end if
    ! if this is called by diagnos_mode_struct, then there is nothing
    ! else to do - the rest is already allocated.
    if(allocate_3d_buffers) return
    if (kykxs_moments .or. kykxs_j0_moments .or. kykxs_j1_moments) then
      call allocate_spec_cmpx_buffers()
    end if
   

    allocate(moment_spec(nmod,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: moment_spec')
    allocate(buf(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: fluxbuf')

    ! Arrays for FFT of a perpendicular slice
    if(io_legacy) then
      allocate(a(mphiw3t,mrad_l), stat = ierr)
    else
      allocate(a(nymod,nx), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('diagnos_moments: Could not allocate a in diagnostic')
    endif
    if(io_legacy) then
      allocate(ar(mphit,mrad_l), stat = ierr)
    else
      allocate(ar(n_y_grid,nx), stat = ierr)
    end if
    if (ierr.ne.0) then 
      call gkw_abort('diagnos_moments: Could not allocate ar in diagnostic')
    endif
    
    if(den3d.or.ene3d.or.xy_dens.or.xy_temp.or. &
       & xy_current.or.xy_current2) then
      if(io_legacy) then
        allocate(c_xy(mphiw3t,mrad_l),stat=ierr)
        if (ierr /= 0) call gkw_abort('diagnostic :: c_xy')
      else
        allocate(c_xysp(nmod,nx,ns),stat=ierr)
        if (ierr /= 0) call gkw_abort('diagnostic :: c_xysp')
      end if
    end if

    allocate(lun_moment_spec(N_MOMENTS,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: lun_moment_spec')
    allocate(lun_moment_xy(N_MOMENTS,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: lun_moment_xy')
    allocate(lun_moment_xs(N_MOMENTS,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: lun_moment_xs')

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()

    if(allocated(moment_spec)) deallocate(moment_spec)
    if(allocated(out3dbuf_local_real)) then
      deallocate(out3dbuf_local_real)
    end if
    if(allocated(out3dbuf_global_real)) then
      deallocate(out3dbuf_global_real)
    end if
    

  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()
    
    
    ! One could add generic switch linitial_output
    !if (xs_kyzero_dens) call output_moment_xs_kyzero(DENSITY_MOMENT)
    !if (xs_kyzero_ene) call output_moment_xs_kyzero(TEMPERATURE_MOMENT)
    !if (xs_kyzero_ene_par) call output_moment_xs_kyzero(PAR_TEMPERATURE_MOMENT)
    !if (xs_kyzero_ene_perp) call output_moment_xs_kyzero(PERP_TEMPERATURE_MOMENT)
    !if (xs_kyzero_current) call output_moment_xs_kyzero(CURRENT_MOMENT)
    !if (xs_kyzero_current2) call output_moment_xs_kyzero(CURRENT_SQ_MOMENT)
    !if (xs_kyzero_phi_ga_fm) call output_moment_xs_kyzero(PHI_GA_FM_MOMENT)
    !if (xs_kyzero_phi_ga_deltaf) call output_moment_xs_kyzero(PHI_GA_DELTAF_MOMENT)
    
    if(kykxs_moments) then
      call output_moments_kykxs(0)
    end if

  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine final_output(number)
    use diagnos_generic, only : xy_estep
    use dist, only : fdisi, nsolc
    integer, intent(in), optional :: number

    ! To keep the compiler quiet.
    if (.false.) write(*,*) number

    if(.not.xy_estep) then
      ! moments spectra
      !if (lfluxes_spectra) then ! legacy switch
        call output_momentspec_ky(fdisi(1:nsolc),DENSITY_MOMENT,i_denspec)
        call output_momentspec_ky(fdisi(1:nsolc),TEMPERATURE_MOMENT,i_enespec)
      !end if
      if(kykxs_moments .or. kykxs_j0_moments .or. kykxs_j1_moments) then
         call output_moments_kykxs(number)
      end if
    end if

  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output(file_count)
    use dist, only : fdisi, nsolc
    use mode, only : mode_box
    use diagnos_generic, only : xy_estep, lphi_diagnostics
    use diagnos_generic, only : den3d, ene3d, xy_dens, xy_temp
    use diagnos_generic, only : xy_current, xy_current2
    use diagnos_generic, only : xy_slice_ipar, out3d_interval
    use control, only : io_legacy, spectral_radius, nt_complete
    use control, only : itime_rst
    use grid, only : number_of_species
    integer, intent(in) :: file_count

    ! Shouldn't we get the local f once for all here and pass it to the 
    ! various moment routines instead of using get_f_from_g everywhere?

    if(xy_estep) then
      ! Complex 3D moments on a (kx,ky,s) grid
      ! WARNING, data volume can be huge.
      if(kykxs_moments .or. kykxs_j0_moments .or. kykxs_j1_moments) then
         call output_moments_kykxs(file_count)
      end if

      ! moments spectra
      !if (lfluxes_spectra) then
        call output_momentspec_ky(fdisi(1:nsolc),DENSITY_MOMENT,i_denspec)
        call output_momentspec_ky(fdisi(1:nsolc),TEMPERATURE_MOMENT,i_enespec)
      !end if
  
      if (.not. (mode_box .and. lphi_diagnostics)) return

      ! Outputs 3D files for the whole flux tube  every 
      ! out3d_interval'th large timestep.
      ! WARNING, data volume can be huge.
      if(modulo(nt_complete+itime_rst,out3d_interval) == 0) then
        if(den3d) call output_moment_3d(DENSITY_MOMENT,file_count,tag_range_start)
        if(ene3d) call output_moment_3d(TEMPERATURE_MOMENT,file_count,tag_range_start+number_of_species)
      end if
      
      if(io_legacy) then
        if (xy_dens) call output_moment_perp2d_legacy(DENSITY_MOMENT,xy_slice_ipar)
        if (xy_temp) call output_moment_perp2d_legacy(TEMPERATURE_MOMENT,xy_slice_ipar)
        if (xy_current) call output_moment_perp2d_legacy(CURRENT_MOMENT,xy_slice_ipar)
        if (xy_current2) call output_moment_perp2d_legacy(CURRENT_SQ_MOMENT,xy_slice_ipar)
      else
        if (xy_dens) call output_moment_perp2d(DENSITY_MOMENT,xy_slice_ipar)
        if (xy_temp) call output_moment_perp2d(TEMPERATURE_MOMENT,xy_slice_ipar)
        if (xy_current) call output_moment_perp2d(CURRENT_MOMENT,xy_slice_ipar)
        if (xy_current2) call output_moment_perp2d(CURRENT_SQ_MOMENT,xy_slice_ipar)
      end if

    end if
    
    ! xs-output of iyzero mode
    if (xs_kyzero_dens) call output_moment_xs_kyzero(DENSITY_MOMENT)
    if (xs_kyzero_ene) call output_moment_xs_kyzero(TEMPERATURE_MOMENT)
    if (xs_kyzero_ene_par) call output_moment_xs_kyzero(PAR_TEMPERATURE_MOMENT)
    if (xs_kyzero_ene_perp) call output_moment_xs_kyzero(PERP_TEMPERATURE_MOMENT)
    if (xs_kyzero_current) call output_moment_xs_kyzero(CURRENT_MOMENT)
    if (xs_kyzero_current2) call output_moment_xs_kyzero(CURRENT_SQ_MOMENT)
    if (xs_kyzero_phi_ga_fm) call output_moment_xs_kyzero(PHI_GA_FM_MOMENT)
    if (xs_kyzero_phi_ga_deltaf) call output_moment_xs_kyzero(PHI_GA_DELTAF_MOMENT)

  end subroutine output

  !---------------------------------------------------------------------------
  !> Generic routine to write 2D perpendicular slices of a distribution
  !> function moment over all velocity space for each species
  !> individually at global parallel grid point ig.
  !>
  !> This routine outputs position space data, as well as spectral data.
  !>
  !---------------------------------------------------------------------------
  subroutine output_moment_perp2d(momentname, ig)
    use dist,            only : fdisi
    use velocitygrid,    only : intvp,intmu,vpgr,mugr
    use grid,            only : nmod,nx,nmu,nvpar, ls, ns, nsp
    use grid,            only : proc_subset
    use geom,            only : bn
    use general,         only : gkw_abort
    use matdat,          only : get_f_from_g 
    use mpiinterface,    only : mpiallreduce_sum_inplace
    use mpicomms,        only : COMM_VPAR_NE_MU_NE
    use diagnos_generic, only : xy_spec
    use diagnos_generic, only : xysp_output_array, kyxsp_output_array


    !> String defines the type of the moment.
    integer, intent(in) :: momentname
    !> global s point to output
    integer, intent(in) :: ig
    complex :: mom_scalar
    integer :: j,k,imod,ix,isp
    real    :: integrand = 0.0
    
    if(proc_subset(0,ig,0,0,0)) then
      species_loop : do isp = 1, nsp

        ! loop over the slice points for this species
        local_slice : do imod = 1, nmod
          do ix = 1, nx

            ! integrate over the velocity space for this point
            mom_scalar=(0.E0,0.E0)
            velocity_space : do j= 1, nmu
              do k = 1, nvpar

                select case(momentname)
                case(DENSITY_MOMENT)
                  integrand = 1.E0
                case(TEMPERATURE_MOMENT)
                  integrand = vpgr(ls(ig),j,k,isp)**2 + 2.E0*mugr(j)*bn(ix,ls(ig))
                case(CURRENT_MOMENT)
                  integrand = vpgr(ls(ig),j,k,isp)
                case(CURRENT_SQ_MOMENT)
                  integrand = vpgr(ls(ig),j,k,isp)**2
                case(VPERP_SQ_MOMENT)
                  integrand = sqrt(2.E0*mugr(j)*bn(ix,ls(ig)))
                case default
                  call gkw_abort('output_moment_perp2d, wrong moment labels')
                end select
                ! N.B. this is probably slow!
                mom_scalar = mom_scalar + integrand*bn(ix,ls(ig))* &
                   & intvp(ls(ig),j,k,isp)*intmu(j)* &
                   & get_f_from_g(imod,ix,ls(ig),j,k,isp,fdisi)
              end do
            end do velocity_space

            ! store the scalar into the slice
            c_xysp(imod,ix,isp) = mom_scalar

          end do
        end do local_slice

      end do species_loop
      ! finish the velocity space integral
      call mpiallreduce_sum_inplace(c_xysp, shape(c_xysp), &
         & COMM_VPAR_NE_MU_NE)
    else
      ! processes not working on parallel coordinate ig
      ! will just hold zeros.
      c_xysp = 0
    end if

    ! choose the logical unit according to the global species
    call xysp_output_array( &
       & c_xysp, &
       & lun_moment_xy(momentname,1), &
       & proc_subset(0,ig,1,1,0), ig <= ns)

    if(xy_spec) then
      call kyxsp_output_array( &
         & abs(c_xysp), &
         & lun_moment_spec(momentname,1), &
         & proc_subset(0,ig,1,1,0), ig <= ns)
    end if


  end subroutine output_moment_perp2d
  
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> 
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine output_moment_xs_kyzero(momentname)
  
    use mode,            only : iyzero
    use dist,            only : fdisi, fmaxwl
    use velocitygrid,    only : intvp,intmu,vpgr,mugr
    use grid,            only : nx,nmu,nvpar, lsp, ns, nsp
    use grid,            only : number_of_species, proc_subset
    use geom,            only : bn
    use general,         only : gkw_abort
    use matdat,          only : get_f_from_g 
    use mpiinterface,    only : mpiallreduce_sum_inplace
    use mpicomms,        only : COMM_VPAR_NE_MU_NE
    use diagnos_generic, only : xs_kyzero_output_array
    use fields,          only : get_averaged_phi


    complex :: mom_scalar
    integer :: j,k,imod,ix,isp,is
    complex    :: integrand = (0.E0,0.E0)
    integer :: lun, ierr
    !> String defines the type of the moment.
    integer, intent(in) :: momentname
    
    
    if(.not. allocated(mom_xs_local)) then
      allocate(mom_xs_local(nx,ns),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_moments :: Cannot allocate mom_xs_local!')
    end if

    ! only for ky zero mode
    imod = iyzero

    species_loop : do isp = 1, number_of_species
      
      velspace_integration: if(proc_subset(0,0,0,0,isp)) then
      
        mom_xs_local = (0.E0, 0.E0)
        
        ! loop over the slice points for this species
        local_slice : do ix = 1, nx; do is = 1, ns

          ! integrate over the velocity space for this point
          mom_scalar = (0.E0, 0.E0)

          velocity_space : do j= 1, nmu; do k = 1, nvpar

            select case(momentname)
            case(DENSITY_MOMENT)
              integrand = 1.E0 * get_f_from_g(imod,ix,is,j,k,lsp(isp),fdisi)
            case(TEMPERATURE_MOMENT)
              integrand = (vpgr(is,j,k,lsp(isp))**2 + 2.E0*mugr(j)*bn(ix,is)) &
                        & * get_f_from_g(imod,ix,is,j,k,lsp(isp),fdisi)
            case(CURRENT_MOMENT)
              integrand = vpgr(is,j,k,lsp(isp)) * get_f_from_g(imod,ix,is,j,k,lsp(isp),fdisi)
            case(CURRENT_SQ_MOMENT)
              integrand = vpgr(is,j,k,lsp(isp))**2 * get_f_from_g(imod,ix,is,j,k,lsp(isp),fdisi)
            case(PAR_TEMPERATURE_MOMENT)
              integrand = vpgr(is,j,k,lsp(isp))**2 * get_f_from_g(imod,ix,is,j,k,lsp(isp),fdisi)
            case(PERP_TEMPERATURE_MOMENT)
              integrand = 2.E0*mugr(j)*bn(ix,is) * get_f_from_g(imod,ix,is,j,k,lsp(isp),fdisi)
            case(PHI_GA_FM_MOMENT)
              integrand = get_averaged_phi(imod,ix,is,j,lsp(isp),fdisi) &
                        & * fmaxwl(ix,is,j,k,lsp(isp))
            case(PHI_GA_DELTAF_MOMENT)
              integrand = get_averaged_phi(imod,ix,is,j,lsp(isp),fdisi) &
                        & * get_f_from_g(imod,ix,is,j,k,lsp(isp),fdisi)
            case default
              call gkw_abort('output_moment_xs_kyzero, wrong moment label')
            end select
            
            ! multiply with velocity space jacobian and volumen element
            mom_scalar = mom_scalar + integrand * bn(ix,is) * &
              & intvp(is,j,k,lsp(isp)) * intmu(j)
              
          end do; end do velocity_space

          ! store the scalar into the slice
          mom_xs_local(ix,is) = mom_xs_local(ix,is) + mom_scalar

        end do; end do local_slice
        
        ! finish the velocity space integral
        call mpiallreduce_sum_inplace(mom_xs_local(1:nx,1:ns), (/ nx, ns /), &
           & COMM_VPAR_NE_MU_NE)

      end if velspace_integration

    
      ! get output unit
      lun = lun_moment_xs(momentname,isp)
      
      
      ! output xs-moment
      call xs_kyzero_output_array( &
        & mom_xs_local(1:nx,1:ns), &
        & lun_moment_xs(momentname, isp), &
        & proc_subset(0,0,1,1,isp), isp <= nsp)

    end do species_loop
    
  end subroutine output_moment_xs_kyzero
  
  
  

  !---------------------------------------------------------------------------
  !> Generic routine to write 2D perpendicular slices of a distribution
  !> function moment over all velocity space for each species
  !> individually at global parallel grid point ig.
  !>
  !> This routine outputs position space data, as well as spectral data.
  !>
  !---------------------------------------------------------------------------
  subroutine output_moment_perp2d_legacy(momentname, ig)
    use dist,            only : fdisi
    use velocitygrid,    only : intvp,intmu,vpgr,mugr
    use grid,            only : nmod,nx,nmu,nvpar, ls, lsp, ns, nsp
    use grid,            only : number_of_species, proc_subset
    use geom,            only : bn
    use general,         only : gkw_abort
    use matdat,          only : get_f_from_g 
    use mpiinterface,    only : mpiallreduce_sum_inplace
    use mpicomms,        only : COMM_VPAR_NE_MU_NE
    use diagnos_generic, only : xy_spec
    use diagnos_generic, only : xy_output_array, kyx_output_array

    !> String defines the type of the moment.
    integer, intent(in) :: momentname
    !> global s point to output
    integer, intent(in) :: ig
    !FJC future: integer, intent(in) :: isg !< global species number to output
    complex :: mom_scalar
    integer :: j,k,imod,ix,isp
    real    :: integrand = 0.0

    ! Return if the point ig is not on the local processor.
    !if (.not. proc_subset(0,ig,0,0,0)) return
    !FJC future: if .not. proc_subset(0,ig,0,0,is) return

    species_loop : do isp = 1, number_of_species
      c_xy = 0
      
      if(proc_subset(0,ig,0,0,isp)) then

        ! loop over the slice points for this species
        local_slice : do imod = 1, nmod
          do ix = 1, nx

            ! integrate over the velocity space for this point
            mom_scalar=(0.E0,0.E0)
            velocity_space : do j= 1, nmu
              do k = 1, nvpar

                select case(momentname)
                case(DENSITY_MOMENT)
                  integrand = 1.E0
                case(TEMPERATURE_MOMENT)
                  integrand = vpgr(ls(ig),j,k,lsp(isp))**2 + 2.E0*mugr(j)*bn(ix,ls(ig))
                case(CURRENT_MOMENT)
                  integrand = vpgr(ls(ig),j,k,lsp(isp))
                case(CURRENT_SQ_MOMENT)
                  integrand = vpgr(ls(ig),j,k,lsp(isp))**2
                case(VPERP_SQ_MOMENT)
                  integrand = sqrt(2.E0*mugr(j)*bn(ix,ls(ig)))
                case default
                  call gkw_abort('output_moment_perp2d, wrong moment labels')
                end select
                ! N.B. this is probably slow!
                mom_scalar = mom_scalar + integrand*bn(ix,ls(ig))* &
                   & intvp(ls(ig),j,k,lsp(isp))*intmu(j)* &
                   & get_f_from_g(imod,ix,ls(ig),j,k,lsp(isp),fdisi)
              end do
            end do velocity_space

            ! store the scalar into the slice
            c_xy(imod,ix) = mom_scalar

          end do
        end do local_slice

        ! finish the velocity space integral
        call mpiallreduce_sum_inplace(c_xy(1:nmod,1:nx), (/ nmod, nx /), &
           & COMM_VPAR_NE_MU_NE)
      end if
      ! processes not working on parallel coordinate ig and species isp
      ! will just hold zeros.

      ! choose the logical unit according to the global species
      call xy_output_array( &
         & c_xy(1:nmod,1:nx), .false., &
         & lun_moment_xy(momentname, isp), &
         & proc_subset(0,ig,1,1,isp), ig <= ns .and. isp <= nsp)

      ! Select one processor x slice to write
      if(xy_spec) then
        call kyx_output_array( &
           & abs(c_xy(1:nmod,1:nx)), &
           & lun_moment_spec(momentname, isp), &
           & proc_subset(0,ig,1,1,isp), ig <= ns .and. isp <= nsp)
      end if

    end do species_loop

  end subroutine output_moment_perp2d_legacy

  !--------------------------------------------------------------------
  !> Outputs ky spectra of flux surface averaged moments.
  !--------------------------------------------------------------------
  subroutine output_momentspec_ky(fdis,momentname,lun)
    use mpiinterface,     only : mpiallreduce_sum, root_processor
    !use mpicomms,         only : COMM_SP_EQ
    use grid,             only : nmod, nx, ns, nmu, nvpar, proc_subset
    use grid,             only : lsp, number_of_species
    use dist,             only : nsolc
    use geom,             only : bn, ints
    use velocitygrid,     only : intvp,intmu,vpgr,mugr
    use matdat,           only : get_f_from_g
    use io,               only : append_chunk, xy_fmt, ascii_fmt
    use general,          only : gkw_abort
    use control,          only : io_legacy

    complex, dimension(nsolc), intent(in) :: fdis
    integer, intent(in) :: momentname
    integer, intent(in)  :: lun

    real :: integrand = 0.0, scalar, fdisi
    integer :: ipar,j,k,imod,ix,is

    moment_spec(:,:) = 0.0
    buf(:,:,:) = 0.0

    species : do is=1,number_of_species
      procset : if (proc_subset(0,0,0,0,is)) then
        mod: do imod=1,nmod
          scalar= 0.0
          ! loop through the local field line points and radial wavevectors
          sx : do ipar = 1, ns
            do ix = 1, nx
              ! integrate over the velocity space for this point
              velocity_space : do j=1,nmu
                do k=1,nvpar

                  select case(momentname)
                  case(DENSITY_MOMENT)
                    integrand = 1.E0
                  case(TEMPERATURE_MOMENT)
                    integrand = vpgr(ipar,j,k,lsp(is))**2 + &
                       & 2.E0*mugr(j)*bn(ix,ipar)
                  case(CURRENT_MOMENT)
                    integrand = vpgr(ipar,j,k,lsp(is))
                  case default
                    call gkw_abort('Wrong moment called in output_momentspec_ky')
                  end select

                  ! N.B. this is probably slow!
                  fdisi = abs(get_f_from_g(imod,ix,ipar,j,k,lsp(is),fdis))
                  scalar = scalar + bn(ix,ipar)*intvp(ipar,j,k,lsp(is))* &
                     & intmu(j)*fdisi*integrand*ints(ipar)

                end do;
              end do velocity_space
            end do
          end do sx

          ! store the moment in a temporary slice for reduction later
          buf(imod,1,is) = scalar
        end do mod
      end if procset
    end do species

    ! complete the integrals over processors, "gather" over species:
    ! sum over all processors
    call mpiallreduce_sum(buf(:,1,:),moment_spec(:,:), &
       & nmod,number_of_species)
    !SRG thinks the following should be correct, rather than the preceeding
    !line, but why do testcases pass then?! (Remove this comment
    !when you have checked that)
    ! call mpiallreduce_sum(buf(:,1,:),moment_spec(:,:), &
    !    & nmod,number_of_species, COMM_SP_EQ)

    if (root_processor) then
      if(io_legacy) then
        !call append_chunk(lun, &
        !   & (/ (moment_spec(1:nmod,is),is = 1, number_of_species) /), &
        !   & xy_fmt, ascii_fmt)

        ! issue 214#40 bug with PGI compiler here
        ! the shape( (/ (moment_spec(1:nmod,is),is = 1, number_of_species) /) ) 
        ! is not correct inside the passed routine!

        call append_chunk(lun, &
           & reshape(moment_spec, (/ nmod*number_of_species /)) , &
           & xy_fmt, ascii_fmt)

      else
        ! note, for ascii, this causes interleaving of species data in a single column
        ! therefore, io_legacy is somehow always required for ascii
        call append_chunk(lun, &
           & moment_spec(1:nmod,:), &
           & xy_fmt, ascii_fmt)
      end if
    end if

  end subroutine output_momentspec_ky


  !----------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------
  subroutine output_moment_3d(momentname, file_count, tag)
    use general, only : gkw_abort
    
    integer, intent(in) :: momentname
    integer, intent(in) :: file_count
    integer, intent(in) :: tag
    character (len=13) :: prefix

    select case(momentname)
    case(DENSITY_MOMENT)
      prefix = 'D3d'
    case(TEMPERATURE_MOMENT)
      prefix = 'E3d'
    case(CURRENT_MOMENT)
      prefix = 'P3d'
    case default
      call gkw_abort('Wrong moment selected in output_moment_3d')
    end select

    call output_moment_3d_parallel(momentname, prefix, file_count, tag)

  end subroutine output_moment_3d


  !--------------------------------------------------------------------
  !> Calculates moments on the local domain and
  !> passes it to mpi_output_array (which uses MPI-IO and subarray
  !> datatype to write in parallel) in order to write the complete 3d
  !> moment field.
  !> 
  !> Output in hdf5 and ascii format will be done serially, by another
  !> routine.
  !--------------------------------------------------------------------
  subroutine output_moment_3d_parallel(momentname, prefix, file_count, tag)

    use mpiinterface,     only : mpireduce_sum_inplace
    use io,               only : mpi_output_array, xy_fmt, binary_fmt
    use grid,             only : nmod, nx, ns, nmu, nvpar, proc_subset
    use grid,             only : lsp, number_of_species, nsp, n_y_grid
    use mpicomms,         only : COMM_S_NE_X_NE, COMM_VPAR_NE_MU_NE
    use dist,             only : fdisi
    use geom,             only : bn
    use non_linear_terms, only : jind
    use global,           only : int2char_zeros
    use velocitygrid,     only : intvp,intmu,vpgr,mugr
    use matdat,           only : get_f_from_g
    use general,          only : gkw_abort
    use diagnos_generic,  only : zonal_scale_3d, four2real_2D, mphit, mrad_l
    use diagnos_generic,  only : mpi_dtype_real_yxs
    use fft, only : FFT_INVERSE
    use grid,             only : jind_flexible
    use control, only : io_legacy

    character (len=*), intent(in) :: prefix
    integer, intent(in) :: momentname, file_count, tag

    real :: integrand = 0.0
    complex :: fdis,dens
    integer :: i,j,k,imod,ix,ipar,is
    
    character (len=30) :: luname
    logical :: root_is_in_comm

    species_loop : do is = 1, number_of_species

      ! loop through the local field line points 
      local_s_grid : do ipar = 1, ns
        
        procset: if (proc_subset(0,0,0,0,is)) then

          ! loop over the slice points for this species
          local_slice : do imod = 1, nmod
            do ix = 1, nx

              ! integrate over the velocity space for this point
              dens=(0.E0,0.E0)

              velocity_space : do j = 1, nmu
                do k = 1, nvpar

                  select case(momentname)
                  case(DENSITY_MOMENT)
                    integrand = 1.E0
                  case(TEMPERATURE_MOMENT)
                    integrand = vpgr(ipar,j,k,lsp(is))**2 + &
                       & 2.E0*mugr(j)*bn(ix,ipar)
                  case(CURRENT_MOMENT)
                    integrand = vpgr(ipar,j,k,lsp(is))
                  case default
                    call gkw_abort('Wrong moment called in output_moment_3d')
                  end select

                  ! N.B. this is probably slow!
                  fdis = get_f_from_g(imod,ix,ipar,j,k,lsp(is),fdisi)
                  dens = dens + bn(ix,ipar)*intvp(ipar,j,k,lsp(is))* &
                     & intmu(j)*fdis*integrand
                  
                end do
              end do velocity_space

              if(io_legacy) then
                ! store the density in a temporary slice for reduction later
                c_xy(imod,ix) = dens
              else
                c_xysp(imod, ix, 1) = dens
              end if
            end do
          end do local_slice

          ! finish the velspace integration by summing over other processors
          if(io_legacy) then
            call mpireduce_sum_inplace(c_xy(1:nmod,1:nx), &
               & (/nmod, nx/), COMM_VPAR_NE_MU_NE)
          else
            call mpireduce_sum_inplace(c_xysp(1:nmod,1:nx,1), &
               & (/nmod, nx/), COMM_VPAR_NE_MU_NE)
          end if
          ! Note that here it is assumed that the topology is such
          ! that the root processes of the COMM_VPAR_NE_MU_NE
          ! communicators all work on the (vpar, mu) = (1,1)
          ! gridpoints.
        end if procset
        if(proc_subset(0,0,1,1,is)) then
          if(io_legacy) then
            ! reorder the radial direction if this is a spectral run,
            ! because the fft libs expects the zeromode at index 1,
            ! not in the middle of the grid
            a = (0.,0.)
            do ix = 1, nx
              do imod = 1, nmod
                a(imod,jind(ix)) = c_xy(imod,ix)
              end do
            end do
          else
            a = (0.,0.)
            do ix = 1, nx
              do imod = 1, nmod
                if(nmod == 1) then
                  a(imod+1,jind_flexible(nx,nx,ix)) = c_xysp(imod,ix,1)
                else
                  a(imod,jind_flexible(nx,nx,ix)) = c_xysp(imod,ix,1)
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
        end if

      end do local_s_grid

      ! create a luname
      luname = trim(prefix)//trim(int2char_zeros(is,2))//'_'// &
         & trim(int2char_zeros(file_count,6))

      !on an arbitrary process, find out if root (of comm_cart)
      !works on this global species.
      root_is_in_comm = (is <= nsp)
      ! *all* processes go into IO...
      call mpi_output_array(luname, 'diagnostic/diagnos_moments', &
         & out3dbuf_local_real, mpi_dtype_real_yxs, &
         & out3dbuf_global_real, &
         & COMM_S_NE_X_NE, xy_fmt, binary_fmt, proc_subset(0,0,1,1,is), &
         & root_is_in_comm, tag+is-1)
    end do species_loop

  end subroutine output_moment_3d_parallel


  !----------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------
  subroutine output_moments_kykxs(file_count)
    use general, only : gkw_abort
    use grid, only : number_of_species
    use diagnos_generic, only : tokens
    use diagnos_generic, only : START_MOMENTS

    integer, intent(in) :: file_count
    integer :: momentname(7), i_moment
    character (len=14) :: prefix

    if (kykxs_moments) then
      
      momentname(1)=DENSITY_MOMENT
      momentname(2)=CURRENT_MOMENT
      momentname(3)=PAR_TEMPERATURE_MOMENT
      momentname(4)=PERP_TEMPERATURE_MOMENT
      momentname(5)=QPAR_MOMENT
      momentname(6)=M12_MOMENT
      momentname(7)=M24_MOMENT
     
      do i_moment=1,7
        prefix=trim(tokens(START_MOMENTS+momentname(i_moment)))//'_kykxs'
        call output_moments_kykxs_parallel(momentname(i_moment),prefix, &
        & file_count, tag_range_start+2*number_of_species+ &
        & (i_moment-1)*number_of_species)
        !TO DO, return the luname and attach metadata
      end do

    end if
      
    if (kykxs_j0_moments) then
      
      momentname(1)=DENSITY_GA_MOMENT
      momentname(2)=CURRENT_GA_MOMENT
      momentname(3)=PAR_TEMPERATURE_GA_MOMENT
      momentname(4)=PERP_TEMPERATURE_GA_MOMENT
      momentname(5)=QPAR_GA_MOMENT
      momentname(6)=M12_GA_MOMENT
      momentname(7)=M24_GA_MOMENT

      do i_moment=1,7
        prefix=trim(tokens(START_MOMENTS+momentname(i_moment)))//'_kykxs'
        call output_moments_kykxs_parallel(momentname(i_moment),prefix, &
           & file_count, tag_range_start+2*number_of_species+ &
           & 7*number_of_species+(i_moment-1)*number_of_species)
        !TO DO, return the luname and attach metadata
      end do
    end if

    if (kykxs_j1_moments) then
      
      momentname(1)=PERP_TEMPERATURE_J1_MOMENT
     
      i_moment=1
        prefix=trim(tokens(START_MOMENTS+momentname(i_moment)))//'_kykxs'
        call output_moments_kykxs_parallel(momentname(i_moment),prefix, &
           & file_count, tag_range_start+2*number_of_species+ &
           & 14*number_of_species+(i_moment-1)*number_of_species)
        !TO DO, return the luname and attach metadata
     end if

 end subroutine output_moments_kykxs


  !--------------------------------------------------------------------
  !> Calculates moments on the local domain and
  !> passes it to mpi_output_array (which uses MPI-IO and subarray
  !> datatype to write in parallel) in order to write the complete 3d
  !> moment field.
  !> 
  !> Output in hdf5 and ascii format will be done serially, by another
  !> routine.
  !--------------------------------------------------------------------
  subroutine output_moments_kykxs_parallel(momentname, prefix, file_count, tag)

    use io,               only : mpi_output_array, xy_fmt, binary_fmt
    use grid,             only : nmod, nx, ns, proc_subset
    use grid,             only : lsp, number_of_species, nsp
    use mpicomms,         only : COMM_S_NE_X_NE, COMM_VPAR_NE_MU_NE
    use mpiinterface,     only : mpiallreduce_sum_inplace
    use global,           only : int2char_zeros
    use diagnos_generic,  only : mpi_dtype_real_spec_yxs, out3dbuf_local_spec_cmpx
    use diagnos_generic,  only : out3dbuf_global_spec_cmpx

    character (len=*), intent(in) :: prefix
    integer, intent(in) :: momentname, file_count, tag

    integer :: imod,ix,ipar,is

    character (len=30) :: luname
    logical :: root_is_in_comm

    complex, dimension(nsp, nmod, ns, nx) :: moment_buf

    moment_buf = get_4d_moment(momentname)

    ! reduce over processes which work on the same spatial
    ! point+species;
    call mpiallreduce_sum_inplace(moment_buf, &
       & shape(moment_buf), COMM_VPAR_NE_MU_NE)

    species_loop : do is = 1, number_of_species

      if(proc_subset(0,0,1,1,is)) then

        ! copy slice into write buffer
          do ipar=1,ns
           do ix=1,nx
            do imod=1,nmod
              out3dbuf_local_spec_cmpx(imod,ix,ipar) = &
                 & moment_buf(lsp(is),imod,ipar,ix)
            end do
          end do
        end do

      end if

      ! create a luname
      luname = trim(prefix)//trim(int2char_zeros(is,2))//'_'// &
         & trim(int2char_zeros(file_count,6))

      !on an arbitrary process, find out if root (of comm_cart)
      !works on this global species.
      root_is_in_comm = (is <= nsp)
      ! *all* processes go into IO...
      call mpi_output_array(luname, 'diagnostic/diagnos_moments', &
         & out3dbuf_local_spec_cmpx, mpi_dtype_real_spec_yxs, &
         & out3dbuf_global_spec_cmpx, &
         & COMM_S_NE_X_NE, xy_fmt, binary_fmt, proc_subset(0,0,1,1,is), &
         & root_is_in_comm, tag+is-1)

    end do species_loop

  end subroutine output_moments_kykxs_parallel


  !----------------------------------------------------------------------
  !> gather the 4D global moment field to root of comm_cart
  !----------------------------------------------------------------------
  subroutine gather_moment_4d_serial(iquantity, moment_buf_global, &
    & get_real_space_field, rotate)
    use grid,           only : nsp, nmod, nx, ns
    use grid,           only : number_of_species, n_x_grid, n_s_grid, n_y_grid
    use grid,           only : proc_subset
    use mpicomms,       only : COMM_VPAR_NE_MU_NE, COMM_VPAR_EQ_MU_EQ
    use control, only : non_linear
    use mpiinterface, only : mpiallreduce_sum_inplace
    use mpiinterface, only : gather_array, MPICOMPLEX_X
    use mpidatatypes, only : create_subarray_datatype
    use general, only : gkw_warn
    use global, only : id_sp, id_mod, id_x, id_s
        
    integer, intent(in) :: iquantity
    complex, dimension(number_of_species, nmod, n_s_grid, n_x_grid), intent(out) :: moment_buf_global
    logical, intent(in) :: get_real_space_field
    complex, optional, intent(in) :: rotate(nmod)

    integer :: imod
    complex, dimension(nsp, nmod, ns, nx) :: moment_buf

    integer :: mpi_dtype

    logical, parameter :: to_root_of_commcart = .true.
    logical, parameter :: root_is_in_comm = .true.


    moment_buf = get_4d_moment(iquantity)
    
    loop_nmod: do imod = 1, nmod
      ! Rotate all values in the complex plane relative to maximum potential
      ! defined to be (1.,0.).  Effectively, normalise out frequency.
      if (.not.non_linear) then
        if (present(rotate)) then
          moment_buf(:,imod,:,:) = moment_buf(:,imod,:,:) * rotate(imod)
        endif
      endif
    end do loop_nmod
    
    ! reduce over processes which work on the same spatial
    ! point+species;
    call mpiallreduce_sum_inplace(moment_buf, &
       & shape(moment_buf), COMM_VPAR_NE_MU_NE)

    if(get_real_space_field) then
      ! first, perform an FFT on the local array...
      
      !TODO (needs change of moment_buf_global shape, too)

      ! create subarray dataype for 4d real space field
      call create_subarray_datatype(MPICOMPLEX_X, mpi_dtype, &
         & id_sp, id_mod, id_s, id_x, &
         & global_and_local_ysize=n_y_grid)
      
      ! THIS IS NOT YET IMPLEMENTED - and not used at the moment
      call gkw_warn('not implemented: a real space field is requested from&
         & gather_moment_4d_serial.')

    else
      ! create subarray MPI datatype for 4d spectral field
      call create_subarray_datatype(MPICOMPLEX_X, mpi_dtype, &
         & id_sp, id_mod, id_s, id_x)

    end if

    ! gather and return the 4D global moment
    !if(proc_is_with_root_in(COMM_VPAR_EQ_MU_EQ)) then
    if(proc_subset(0,0,1,1,0)) then
      ! we want to gather data from the processes which are in
      ! COMM_VPAR_EQ_MU_EQ, together with the root process.
      ! (there as many different COMM_VPAR_EQ_MU_EQ communicators
      ! as there are processes in COMM_VPAR_NE_MU_NE)
      call gather_array(moment_buf_global, &
         & moment_buf, &
         & mpi_dtype, COMM_VPAR_EQ_MU_EQ, to_root_of_commcart, root_is_in_comm)
      
    end if
    
  end subroutine gather_moment_4d_serial

  !----------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------
  function get_4d_moment(iquantity) result(moment_buf)
    use global, only : r_tiny
    use control, only : spectral_radius
    use dist, only : fdisi, phi, fmaxwl
    use matdat, only : get_f_from_g
    use grid, only : nsp, nmod, ns, nx, nvpar, nmu, gsp
    use geom, only : bn, bmax, ints
    use velocitygrid, only : intmu, intvp, mugr, vpgr
    use general, only : gkw_abort
    use functions,      only : besselj0_gkw,  mod_besselj1_gkw
    use components,     only : adiabatic_electrons, tmp, signz
    integer, intent(in) :: iquantity
    
    complex :: moment_buf(nsp,nmod,ns,nx)
    complex :: fdis, phi_fsa
    integer :: imod, ix, i, is, k, j
    integer :: mask
    moment_buf = 0.0

    ! If we want to calculate a moment of the distribution,
    ! then we have to compute a velocity space integral.

    loop_nmod: do imod = 1, nmod
      loop_nx: do ix = 1, nx
        ! compute the fsa of a potential line
        phi_fsa = 0.0
        do i = 1, ns
          phi_fsa = phi_fsa + phi(imod,ix,i)*ints(i)
        end do
        
        loop_ns: do i = 1, ns
          loop_nsp: do is = 1, nsp
            loop_nvpar: do k = 1, nvpar
              loop_nmu: do j = 1, nmu
                ! The distribution function
                fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)

                ! There might be gridcells where intvp = 0 ?
                ! Still, I do not see why it is necessary to set fdis = 0
                ! as it used to be done in the parallel_output() routine.
                !if (intvp(i,j,k,is).eq.0.) fdis = 0.
                
                ! mask for trapped passing particles
                if (abs(bmax(ix)-bn(ix,i)) > r_tiny) then
                  if(sqrt(2*mugr(j)*bn(ix,i)) > abs(vpgr(i,j,k,is)*sqrt(bn(ix,i) &
                     & /(bmax(ix)-bn(ix,i))))) then
                    mask= 1
                  else
                    mask = 0
                  end if
                else
                  mask= 0
                end if

                select case(iquantity)

                 ! The following four moments are as they were in 
                 ! the 'old' parallel.dat output
                case(DENSITY_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & fdis
                case(CURRENT_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & vpgr(i,j,k,is)* &
                     & fdis
                case(PAR_TEMPERATURE_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & vpgr(i,j,k,is)**2* &
                     & fdis
                case(PERP_TEMPERATURE_MOMENT)
                  ! with T=T_//+2*T_perp
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & mugr(j)*bn(ix,i)* &
                     & fdis
                case(QPAR_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & vpgr(i,j,k,is)* &
                     & (vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i)) *&
                     & fdis
                case(M12_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & vpgr(i,j,k,is)* &
                     & 2.E0*mugr(j)*bn(ix,i)*&
                     & fdis
                case(M24_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & 2.E0*mugr(j)*bn(ix,i)*&
                     & (vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i)) *&
                     & fdis

                  ! These ones include the Bessel function
                  ! They can be used to compute fluxes or non-linear terms
                case(DENSITY_GA_MOMENT)
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       ! the adjoint gyroavg operator (in the local
                       ! limit identical to the gyroavg operator),
                       ! acting on f
                       & fdis * besselj0_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged density moment is not &
                       & implemented in get_4d_moment')
                  end if
                case(CURRENT_GA_MOMENT)
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & vpgr(i,j,k,is)* &
                       & fdis * besselj0_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged parallel velocity moment is not &
                       & implemented in get_4d_moment')
                  end if
                case(PAR_TEMPERATURE_GA_MOMENT)
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & vpgr(i,j,k,is)**2* &
                       & fdis * besselj0_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged parallel temperature moment is not &
                       & implemented in get_4d_moment')
                  end if
                case(PERP_TEMPERATURE_GA_MOMENT)
                  ! with T=T_//+2*T_perp
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & mugr(j)*bn(ix,i)* &
                       & fdis * besselj0_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged perpendicular temperature moment is not &
                       & implemented in get_4d_moment')
                  end if
                 case(QPAR_GA_MOMENT)
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & vpgr(i,j,k,is)* &
                       & (vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i)) *&
                       & fdis * besselj0_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged parallel heat flux moment is not &
                       & implemented in get_4d_moment')
                  end if
                 case(M12_GA_MOMENT)
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & vpgr(i,j,k,is)* &
                       & 2.E0*mugr(j)*bn(ix,i)*&
                       & fdis * besselj0_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged vpar*vperp^2 moment is not &
                       & implemented in get_4d_moment')
                  end if
                 case(M24_GA_MOMENT)
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & 2.E0*mugr(j)*bn(ix,i)*&
                       & (vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i)) *&
                       & fdis * besselj0_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged vperp^2*v^2 moment is not &
                       & implemented in get_4d_moment')
                  end if


               case(PERP_TEMPERATURE_J1_MOMENT)
                  ! with T=T_//+2*T_perp
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & mugr(j)*bn(ix,i)* &
                       & fdis * mod_besselj1_gkw(imod,ix,i,j,is)
                  else
                    call gkw_abort('Gyroaveraged perpendicular temperature moment is not &
                       & implemented in get_4d_moment')
                  end if



                  ! Other useful moments
               case(PASSING_DENSITY_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & fdis * (-mask+1)
                case(TRAPPED_DENSITY_MOMENT)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & fdis * mask
                case(DENSITY_POLAR_MOMENT)
                  if(spectral_radius) then
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & (signz(is)/tmp(ix,is))*fmaxwl(ix,i,j,k,is)*&
                       & (besselj0_gkw(imod,ix,i,j,is)**2*phi(imod,ix,i) &
                       & - phi(imod,ix,i))
                  else
                    call gkw_abort('Polarisation density moment is not &
                       & implemented in get_4d_moment')
                  end if
                  if(gsp(is) == 1 .and. adiabatic_electrons) then
                    ! the adiabatic electron correction to the
                    ! polarisation density.
                    moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                       & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                       & (signz(is)/tmp(ix,is))* &
                       & (phi(imod,ix,i) - phi_fsa)
                  end if
                case(TEMPERATURE_MOMENT)
                  !(same integrand as in output_moment_3d_parallel)
                  moment_buf(is, imod, i, ix) = moment_buf(is, imod, i, ix) + &
                     & bn(ix,i)*intvp(i,j,k,is)*intmu(j)* &
                     & (vpgr(i,j,k,is)**2 + &
                     & 2.E0*mugr(j)*bn(ix,i)) *&
                     & fdis
  

                case default
                  call gkw_abort('Unknown iquantity passed to &
                     & get_4d_moment')
                end select
              end do loop_nmu
            end do loop_nvpar
          end do loop_nsp
        end do loop_ns
      end do loop_nx

    end do loop_nmod
  end function get_4d_moment

end module diagnos_moments
