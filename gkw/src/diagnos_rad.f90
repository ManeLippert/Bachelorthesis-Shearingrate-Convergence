!------------------------------------------------------------------------------
!> Outputs radial profiles for non-spectral runs.
!> A mixture of fields, moments, energetics, and background profiles.
!> Could be logically re-organised into other relevant modules.
!------------------------------------------------------------------------------
module diagnos_rad

  implicit none

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, output

  private

  !> the luns for the radial profile files
  integer, save, allocatable :: i_delta_n(:), i_delta_v(:), i_delta_t(:), i_delta_e(:)
  integer, save :: i_fields
  !> luns for turbulence and zonal mode intensity
  integer, save :: i_f_int, i_zf_int
  !> luns for turbulence and zonal mode entropy
  integer, save :: i_f_entropy, i_zf_entropy
  !> luns for  turbulence and zonal mode entropy fields
  integer, save :: i_f_entropy_fields, i_zf_entropy_fields
  !> luns for kparallel diagnostic and kparallel^2 and theta_0
  integer, save :: i_kpar, i_kpar2, i_theta0, i_kr, i_kr2

  !> zonal mode entropy flux ExB particles and fields real part
  integer, save :: i_zf_entropy_flux_re, i_zf_entropy_flux_fields_re
  !> zonal mode entropy flux ExB particles and fields imag part
  integer, save :: i_zf_entropy_flux_im, i_zf_entropy_flux_fields_im
  !> perturbations entropy flux ExB 
  integer, save :: i_f_entropy_flux
  !> zonal mode entropy flux drift particles and fields 
  integer, save :: i_zf_entropy_flux_drift, i_zf_entropy_flux_fields_drift
  !> perturbations entropy flux drift particles and fields 
  integer, save :: i_f_entropy_flux_drift, i_f_entropy_flux_fields_drift

  !> Help arrays
  real, save, allocatable  :: dum_G(:)            !< array global in x
  real, save, allocatable  :: dum(:),dum2(:),dum_buf(:)  !< array local in x

  complex, save, allocatable :: cdum_G(:)         !< array global in x
  complex, save, allocatable :: cdum(:)           !< array local in x
  complex, save, allocatable :: cdum2(:)           !< array local in x
  complex, save, allocatable :: cdum3(:)           !< array local in x

  !> Diagnostic switches
  !> turbulence intensity
  logical, save, public :: lrad_tint
  logical, save, public :: lrad_moment
  logical, save, public :: lrad_field
  logical, save, public :: lrad_entropy
  logical, save, public :: lrad_kpar
  !> radial profiles
  logical, save, public :: lradial_entropy

  !> logical unit numbers for output files (in legacy format/order/grouping)
  integer, save :: i_radial_entropy = -1, i_entropy_flux  = -1

  !> Entropy results
  real, allocatable :: rad_entr(:,:)

  !> mpi datatypes 
  integer, save :: mpi_dtype_radial_line
  
contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides.
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    use control, only : flux_tube
    lrad_field =(.not. flux_tube)
    lrad_moment = .false.
    lrad_tint = .false.
    lrad_entropy = .false.
    lradial_entropy= .false.
    lrad_kpar = .false.
  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lradial_entropy, 1)
    call mpibcast(lrad_field, 1)
    call mpibcast(lrad_moment, 1)
    call mpibcast(lrad_tint, 1)
    call mpibcast(lrad_entropy, 1)
    call mpibcast(lrad_kpar, 1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : nlapar, nlbpar, spectral_radius
    use general, only : gkw_warn

    if (lradial_entropy) then
      if (spectral_radius) then
        call gkw_warn('No radial entropy flux in the spectral version')
        lradial_entropy = .false.
      endif
      if (nlapar) then
        call gkw_warn('No radial entropy flux in EM runs')
        lradial_entropy = .false.
      endif
      if (nlbpar) then
        call gkw_warn('No radial entropy flux in EM runs')
        lradial_entropy = .false.
      endif
    endif

    if (spectral_radius .and. (lrad_field.or.lrad_moment.or. &
      & lrad_tint.or.lrad_entropy))then
      call gkw_warn('The radial profiles have no physical meaning running spectral.')
    endif

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use grid, only : n_x_grid, nx, number_of_species
    use general, only : gkw_abort
    integer :: ierr

    if (lradial_entropy) then
      allocate(rad_entr(n_x_grid,2),stat = ierr)
      if (ierr /= 0) call gkw_abort('Cannot allocate rad_entr in diagnostics')
    endif

    if (lrad_field.or.lrad_moment.or.lrad_tint.or.lrad_entropy.or.lrad_kpar)then
      ! allocate help array
      ierr = 0
      allocate(dum_G(n_x_grid), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate dum_G in diagnos_rad')
      allocate(dum(nx), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate dum in diagnos_rad')
      allocate(dum2(nx), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate dum2 in diagnos_rad')
      allocate(dum_buf(nx), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate dum_buf in diagnos_rad')
      allocate(cdum_G(n_x_grid), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate cdum_G in diagnos_rad')
      allocate(cdum(nx), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate cdum in diagnos_rad')
      allocate(cdum2(nx), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate cdum2 in diagnos_rad')
      allocate(cdum3(nx), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate cdum2 in diagnos_rad')
    end if

    if(lrad_moment) then
      allocate(i_delta_n(number_of_species), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate i_delta_n in diagnos_rad')
      allocate(i_delta_v(number_of_species), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate i_delta_v in diagnos_rad')
      allocate(i_delta_t(number_of_species), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate i_delta_t in diagnos_rad')
      allocate(i_delta_e(number_of_species), stat = ierr)
      if (ierr /= 0) call gkw_abort('unable to allocate i_delta_e in diagnos_rad')
    end if

  end subroutine allocate_mem


  !---------------------------------------------------------------------------
  !> This subroutine opens the necessary files
  !---------------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor, MPIREAL_X
    use mpidatatypes, only : create_subarray_datatype
    use io,           only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use grid,         only : n_x_grid, number_of_species, proc_subset
    use global,       only : int2char_zeros, dotdat, id_x, PHI_FIELD
    use control,      only : io_legacy, spectral_radius, flux_tube
    use mode, only : iyzero
    use general, only : gkw_warn
    use diagnos_generic, only : attach_metadata_grid, LOCAL_DATA, X_GHOSTCELLS
    use diagnos_generic, only : S_GHOSTCELLS
    logical, intent(inout) :: requirements(:,:)
    integer :: is
    character (len=18) ::  filename

    if (iyzero == 0) then
      ! warn and disable here, because iyzero is not initialised when
      ! check() is called
      call gkw_warn('No zero mode exists. Radial diagnostics are disabled.')
      lrad_entropy = .false.
      lrad_field = .false.
      lrad_moment = .false.
      lrad_tint = .false.
      !lrad_kpar=T is useful linearly
    end if

    if(lrad_kpar) then
      requirements(PHI_FIELD,LOCAL_DATA) = .true.
      if(.not.spectral_radius) then
        requirements(PHI_FIELD,X_GHOSTCELLS) = .true.
      end if
      !requirements(PHI_FIELD,S_GHOSTCELLS) = .true.
    else if(lrad_field) then
      requirements(PHI_FIELD,LOCAL_DATA) = .true.
    end if
    
    if (root_processor) then
      if (lradial_entropy) then 
        call open_real_lu(dotdat('entropy_rad',io_legacy), &
           & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_radial_entropy)

        call open_real_lu(dotdat('entropy_flr',io_legacy), &
           & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_entropy_flux)
      endif
    end if

    if (lrad_moment) then
      ! open one file per species, on relevant processors
      ! FIXME This probably does not work with the current serial HDF5 IO!!
      do is = 1, number_of_species
        if (root_processor) then
        
        
          filename="delta_n"//trim(int2char_zeros(is,2))
          call open_real_lu(dotdat(trim(filename),io_legacy), &
             & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
             & ascii_fmt, i_delta_n(is))
          call attach_metadata_grid(i_delta_n(is), 'time', 'xgr', ascii_fmt)
          if(flux_tube) then
            ! in the local_limit the distribution function is normalized by \rho_\ast n_{R0}/v_{th}^3
            ! and, hence, the density is given in units of \rho_\ast n_{R0},
            ! where n_{R0} is the background density of the species at the center of the flux tube
            call attach_metadata(i_delta_n(is), phys_unit_key, '\rho_\ast n_{R0}', ascii_fmt)
          else
            ! in the global case the distribution function is normalized by \rho_\ast n_{G}/w^3
            ! and, hence, the density is given in units of \rho_\ast n_{G},
            ! where n_{G} is the 'grid density'
            call attach_metadata(i_delta_n(is), phys_unit_key, 'n_{G}', ascii_fmt)
          endif
          call attach_metadata(i_delta_n(is), description_key,'This file has &
          & n_x_grid columns and contains the radial profile of the perturbed &
          & density, i. e., the zeroth order velocity space momenent over the &
          & pertrubed distribution function.', ascii_fmt)
          call attach_metadata(i_delta_n(is), comments_key, not_avail, ascii_fmt)


          filename="delta_v"//trim(int2char_zeros(is,2))
          call open_real_lu(dotdat(trim(filename),io_legacy), &
             & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
             & ascii_fmt, i_delta_v(is))
          call attach_metadata_grid(i_delta_v(is), 'time', 'xgr', ascii_fmt)
          if(flux_tube) then
            ! in the local_limit the distribution function is normalized by \rho_\ast n_{R0}/v_{th}^3
            call attach_metadata(i_delta_v(is), phys_unit_key, '\rho_\ast n_{R0} v_{th}', ascii_fmt)
          else
            ! in the global case the distribution function is normalized by \rho_\ast n_{G}/w^3
            call attach_metadata(i_delta_v(is), phys_unit_key, 'n_{G} w', ascii_fmt)
          endif
          call attach_metadata(i_delta_v(is), description_key,'This file has &
          & n_x_grid columns and contains the radial profile of the perturbed &
          & parallel rotation velocity times the background density, i. e., &
          & the first order velocity space momenent of parallel velocity over &
          & the pertrubed distribution function.', ascii_fmt)
          call attach_metadata(i_delta_v(is), comments_key, not_avail, ascii_fmt)
          
          
          filename="delta_e"//trim(int2char_zeros(is,2))
          call open_real_lu(dotdat(trim(filename),io_legacy), &
             & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
             & ascii_fmt, i_delta_e(is))
          call attach_metadata_grid(i_delta_e(is), 'time', 'xgr', ascii_fmt)
          if(flux_tube) then
            ! in the local_limit the distribution function is normalized by \rho_\ast n_{R0}/v_{th}^3
            ! T = m v_{th}^2/2 is the background temperature of the species.
            call attach_metadata(i_delta_e(is), phys_unit_key, '\rho_\ast n_{R0} T', ascii_fmt)
          else
            ! in the global case the distribution function is normalized by \rho_\ast n_{G}/w^3
            ! T_{G} is the grid temperature.
            call attach_metadata(i_delta_e(is), phys_unit_key, 'n_{G} T_{G}', ascii_fmt)
          endif
          call attach_metadata(i_delta_e(is), description_key,'This file has &
          & n_x_grid columns and contains the radial profile of the perturbed &
          & energy density, i. e., the second order moment of the squared modulus of the velocity &
          & over the pertrubed distribution function.', ascii_fmt)
          call attach_metadata(i_delta_e(is), comments_key, not_avail, ascii_fmt)


          filename="delta_t"//trim(int2char_zeros(is,2))
          call open_real_lu(dotdat(trim(filename),io_legacy), &
             & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
             & ascii_fmt, i_delta_t(is))
          call attach_metadata_grid(i_delta_t(is), 'time', 'xgr', ascii_fmt)
          if(flux_tube) then
            ! in the local_limit the distribution function is normalized by \rho_\ast n_{R0}/v_{th}^3
            ! T = m v_{th}^2/2 is the background temperature of the species.
            call attach_metadata(i_delta_t(is), phys_unit_key, '\rho_\ast  T', ascii_fmt)
          else
            ! in the global case the distribution function is normalized by \rho_\ast n_{G}/w^3
            ! T_{G} is the grid temperature.
            call attach_metadata(i_delta_t(is), phys_unit_key, '\rho_\ast T_{G}', ascii_fmt)
          endif
          call attach_metadata(i_delta_t(is), description_key,'This file has &
          & n_x_grid columns and contains the radial profile of the perturbed &
          & temperature. The perturbed temperature is related to the perturbed density and &
          & the perturbed energy through the equation of state of an ideal gas.', ascii_fmt)
          call attach_metadata(i_delta_t(is), comments_key, not_avail, ascii_fmt)
          
          
          
          
        end if
      end do
    end if

    if (root_processor) then
      if (lrad_tint)then
        call open_real_lu(dotdat('f_int',io_legacy), &
           & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_f_int)
        call open_real_lu(dotdat('zf_int',io_legacy), &
           & 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_zf_int)
      endif
      if (lrad_field)then
        call open_real_lu('fields_rad', 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_fields)
      endif
      if (lrad_kpar)then
      call open_real_lu('kpar_rad', 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_kpar)
      call open_real_lu('kpar2_rad', 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_kpar2)
      call open_real_lu('theta0_rad', 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_theta0)
      call open_real_lu('kr_r_rad', 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_kr)
      call open_real_lu('kr_i_rad', 'diagnostic/diagnos_rad', (/ n_x_grid /), &
           & ascii_fmt, i_kr2)
      endif
      if (lrad_entropy)then
        call open_real_lu(dotdat('f_entropy_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_f_entropy)
        call open_real_lu(dotdat('zf_entropy_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy)
        call open_real_lu(dotdat('f_entropy_fields_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_f_entropy_fields)
        call open_real_lu(dotdat('zf_entropy_fields_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy_fields)
        call open_real_lu(dotdat('zf_entropy_flux_re_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy_flux_re)
        call open_real_lu(dotdat('zf_entropy_flux_im_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy_flux_im)
        call open_real_lu(dotdat('zf_entropy_flux_fields_re_radial_profile', &
           & io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy_flux_fields_re)
        call open_real_lu(dotdat('zf_entropy_flux_fields_im_radial_profile', &
           & io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy_flux_fields_im)
        call open_real_lu(dotdat('f_entropy_flux_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_f_entropy_flux)
        call open_real_lu(dotdat('zf_entropy_flux_drift_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy_flux_drift)
        call open_real_lu(dotdat('zf_entropy_flux_fields_drift_radial_profile', &
           & io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_zf_entropy_flux_fields_drift)
        call open_real_lu(dotdat('f_entropy_flux_drift_radial_profile',io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_f_entropy_flux_drift)
        call open_real_lu(dotdat('f_entropy_flux_fields_drift_radial_profile', &
           & io_legacy), &
           & 'diagnostic/diagnos_rad', &
           & (/ n_x_grid /), &
           & ascii_fmt, i_f_entropy_flux_fields_drift)
      endif
    endif

    if (proc_subset(0,1,1,1,0)) then
      ! create subarray MPI datatype, used for gathering the radial
      ! profiles to root
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_radial_line, &
         & id_x)
    end if
  end subroutine init

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()

    call output_background_profiles

  end subroutine initial_output

  !------------------------------------------------------------------------------
  !> write the background profiles on file 
  !------------------------------------------------------------------------------
  subroutine output_background_profiles 
    use io,           only : output_array, xy_fmt, ascii_fmt
    use geom,         only : xgr, qx, shatx
    use grid,         only : nsp, nx, number_of_species, n_x_grid
    use mpiinterface, only : root_processor, gather_array 
    use components,   only : de, tmp, fp, tp, vp, iadia, tm_drive, ecurx
    use components,   only : n_spec
    use mpicomms,     only : COMM_X_NE, COMM_SP_NE
    use general,      only : gkw_abort
    use control,      only : io_legacy
    use global,       only : dotdat
    use io, only : attach_metadata
    use io, only : phys_unit_key, description_key, comments_key, not_avail
    use diagnos_generic, only : attach_metadata_grid

    integer :: ix, ierr, dum
    real, allocatable  :: den_G(:,:), temp_G(:,:), rlt_G(:,:), rln_G(:,:)
    real, allocatable  :: ecurx_G(:), uprim_G(:,:)

    ! allocate help array
    ierr = 0
    allocate(den_G(n_x_grid,n_spec), stat = ierr)
    if (ierr /= 0) call gkw_abort('unable to allocate den_G in diagnos_rad')
    allocate(temp_G(n_x_grid,n_spec), stat = ierr)
    if (ierr /= 0) call gkw_abort('unable to allocate temp_G in diagnos_rad')
    allocate(rlt_G(n_x_grid,n_spec), stat = ierr)
    if (ierr /= 0) call gkw_abort('unable to allocate rlt_G in diagnos_rad')
    allocate(rln_G(n_x_grid,n_spec), stat = ierr)
    if (ierr /= 0) call gkw_abort('unable to allocate rln_G in diagnos_rad')
    allocate(uprim_G(n_x_grid,n_spec), stat = ierr)
    if (ierr /= 0) call gkw_abort('unable to allocate uprim_G in diagnos_rad')
    allocate(ecurx_G(n_x_grid), stat = ierr)
    if (ierr /= 0) call gkw_abort('unable to allocate ecurx_G in diagnos_rad')

    ! gather the density/.. information onto the root processor
    call gather_array(den_G(  1:n_x_grid, 1:number_of_species), n_x_grid,      &
       & number_of_species, de( 1:nx, 1:nsp), nx, nsp, COMM_X_NE, COMM_SP_NE)
    call gather_array(temp_G( 1:n_x_grid, 1:number_of_species), n_x_grid,      &
       & number_of_species, tmp(1:nx, 1:nsp), nx, nsp, COMM_X_NE, COMM_SP_NE)
    call gather_array(rln_G(  1:n_x_grid, 1:number_of_species), n_x_grid,      &
       & number_of_species, fp( 1:nx, 1:nsp), nx, nsp, COMM_X_NE, COMM_SP_NE)
    call gather_array(rlt_G(  1:n_x_grid, 1:number_of_species), n_x_grid,      &
       & number_of_species, tp( 1:nx, 1:nsp), nx, nsp, COMM_X_NE, COMM_SP_NE)
    call gather_array(uprim_G(1:n_x_grid, 1:number_of_species), n_x_grid,      &
       & number_of_species, vp( 1:nx, 1:nsp), nx, nsp, COMM_X_NE, COMM_SP_NE)

    ! If there is an adiabatic species, these information also has to be collected.
    if (0 /= iadia) then
      call gather_array(den_G(  1:n_x_grid, n_spec), n_x_grid,&
         & de( 1:nx, nsp+iadia), nx, COMM_X_NE)
      call gather_array(temp_G( 1:n_x_grid, n_spec), n_x_grid,&
         & tmp(1:nx, nsp+iadia), nx, COMM_X_NE)
      call gather_array(rln_G(  1:n_x_grid, n_spec), n_x_grid,&
         & fp( 1:nx, nsp+iadia), nx, COMM_X_NE)
      call gather_array(rlt_G(  1:n_x_grid, n_spec), n_x_grid,&
         & tp( 1:nx, nsp+iadia), nx, COMM_X_NE)
      call gather_array(uprim_G(1:n_x_grid, n_spec), n_x_grid,&
         & vp( 1:nx, nsp+iadia), nx, COMM_X_NE)
    end if

    if(tm_drive) then
      call gather_array(ecurx_G(1:n_x_grid), n_x_grid,  &
         & ecurx(1:nx),nx,COMM_X_NE, ALLGATHER = .true.)
      dum = number_of_species
      do ix=1, n_x_grid
        !Currently a dirty hack as this wont work with multiple ion species
        !Assumes electrons are the last species.  
        uprim_G(ix,dum) = uprim_G(ix,dum)+ecurx_G(ix)
      enddo
    endif

    if(root_processor) then
      !> Format for the file 'prof_back'.  This file has 3 (x, q,
      !> shat) + 5 * n_spec columns, where n_spec is the total number
      !> of species, adiabatic and kinetic.  The columns for the
      !> species are respectively (n, T_s, R/L_n, R/L_T, u'), where
      !> first all the densities are given, then the temperature, etc.

      call output_array(dotdat('prof_back', io_legacy) , &
         & 'diagnostic/diagnos_rad', reshape(&
         & (/ (xgr(ix), qx(ix), shatx(ix), den_G(ix,:), temp_G(ix,:), &
         &  rln_G(ix,:), rlt_G(ix,:), uprim_G(ix,:), ix = 1, n_x_grid) /), &
         & (/ 3+5*n_spec, n_x_grid /)), &
         & 'F', xy_fmt, ascii_fmt)
      call attach_metadata_grid('prof_back', &
         & 'diagnostic/diagnos_rad', 'xgr', ascii_fmt)
      call attach_metadata('prof_back', &
         & 'diagnostic/diagnos_rad', phys_unit_key, &
         & not_avail, &
         & ascii_fmt)
      call attach_metadata('prof_back', &
         & 'diagnostic/diagnos_rad', description_key, &
         & 'This file has 3 (x, q, shat) + 5 * n_spec columns, where n_spec is &
         & the total number of species, adiabatic and kinetic.  The columns for &
         & the species are respectively (n, T_s, R/L_n, R/L_T, u-prime), where &
         & first all the densities are given, then the temperature, etc.' , ascii_fmt)
      call attach_metadata('prof_back', &
         & 'diagnostic/diagnos_rad', comments_key, &
         & not_avail, ascii_fmt)
      
      if(.not.io_legacy) then

        ! the q profile appears also in geom, but the traditional
        ! format of geom cannot be read by scripts so easily...
        call output_array('qx', &
           & 'diagnostic/diagnos_rad/background_profiles', qx, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata_grid('qx', &
           & 'diagnostic/diagnos_rad/background_profiles', 'xgr', ascii_fmt)
        call attach_metadata('qx', &
           & 'diagnostic/diagnos_rad/background_profiles', phys_unit_key, &
           & not_avail, &
           & ascii_fmt)
        call attach_metadata('qx', &
           & 'diagnostic/diagnos_rad/background_profiles', description_key, &
           & 'radial profile of the safety factor q' , ascii_fmt)
        call attach_metadata('qx', &
           & 'diagnostic/diagnos_rad/background_profiles', comments_key, &
           & not_avail, ascii_fmt)

        ! the magnetic shear profile appears also in geom
        call output_array('shatx', &
           & 'diagnostic/diagnos_rad/background_profiles', shatx, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata_grid('shatx', &
           & 'diagnostic/diagnos_rad/background_profiles', 'xgr', ascii_fmt)
        call attach_metadata('shatx', &
           & 'diagnostic/diagnos_rad/background_profiles', phys_unit_key, &
           & not_avail, &
           & ascii_fmt)
        call attach_metadata('shatx', &
           & 'diagnostic/diagnos_rad/background_profiles', description_key, &
           & 'magnetic shear profile', ascii_fmt)
        call attach_metadata('shatx', &
           & 'diagnostic/diagnos_rad/background_profiles', comments_key, &
           & 'This quantity appears as shat in geom', ascii_fmt)

        call output_array('background_dens', &
           & 'diagnostic/diagnos_rad/background_profiles', den_G, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata_grid('background_dens', &
           & 'diagnostic/diagnos_rad/background_profiles', 'xgr', 'species', &
           & ascii_fmt)
        call attach_metadata('background_dens', &
           & 'diagnostic/diagnos_rad/background_profiles', phys_unit_key, &
           & not_avail, &
           & ascii_fmt)
        call attach_metadata('background_dens', &
           & 'diagnostic/diagnos_rad/background_profiles', description_key, &
           & 'background radial profile of the density', ascii_fmt)
        call attach_metadata('background_dens', &
           & 'diagnostic/diagnos_rad/background_profiles', comments_key, &
           & not_avail, ascii_fmt)

        call output_array('background_temp', &
           & 'diagnostic/diagnos_rad/background_profiles', temp_G, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata_grid('background_temp', &
           & 'diagnostic/diagnos_rad/background_profiles', 'xgr', 'species', &
           & ascii_fmt)
        call attach_metadata('background_temp', &
           & 'diagnostic/diagnos_rad/background_profiles', phys_unit_key, &
           & not_avail, &
           & ascii_fmt)
        call attach_metadata('background_temp', &
           & 'diagnostic/diagnos_rad/background_profiles', description_key, &
           & 'background radial profile of the temperature', ascii_fmt)
        call attach_metadata('background_temp', &
           & 'diagnostic/diagnos_rad/background_profiles', comments_key, &
           & not_avail, ascii_fmt)

        call output_array('background_rln', &
           & 'diagnostic/diagnos_rad/background_profiles', rln_G, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata_grid('background_rln', &
           & 'diagnostic/diagnos_rad/background_profiles', 'xgr', 'species', &
           & ascii_fmt)
        call attach_metadata('background_rln', &
           & 'diagnostic/diagnos_rad/background_profiles', phys_unit_key, &
           & not_avail, &
           & ascii_fmt)
        call attach_metadata('background_rln', &
           & 'diagnostic/diagnos_rad/background_profiles', description_key, &
           & 'background radial profile of the radial density gradient', &
           & ascii_fmt)
        call attach_metadata('background_rln', &
           & 'diagnostic/diagnos_rad/background_profiles', comments_key, &
           & not_avail, ascii_fmt)

        call output_array('background_rlt', &
           & 'diagnostic/diagnos_rad/background_profiles', rlt_G, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata_grid('background_rlt', &
           & 'diagnostic/diagnos_rad/background_profiles', 'xgr', 'species',  &
           & ascii_fmt)
        call attach_metadata('background_rlt', &
           & 'diagnostic/diagnos_rad/background_profiles', phys_unit_key, &
           & not_avail, &
           & ascii_fmt)
        call attach_metadata('background_rlt', &
           & 'diagnostic/diagnos_rad/background_profiles', description_key, &
           & 'background radial profile of the radial temperature gradient', &
           & ascii_fmt)
        call attach_metadata('background_rlt', &
           & 'diagnostic/diagnos_rad/background_profiles', comments_key, &
           & not_avail, ascii_fmt)

        call output_array('background_uprim', &
           & 'diagnostic/diagnos_rad/background_profiles', uprim_G, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata_grid('background_uprim', &
           & 'diagnostic/diagnos_rad/background_profiles', 'xgr', 'species', &
           & ascii_fmt)
        call attach_metadata('background_uprim', &
           & 'diagnostic/diagnos_rad/background_profiles', phys_unit_key, &
           & not_avail, &
           & ascii_fmt)
        call attach_metadata('background_uprim', &
           & 'diagnostic/diagnos_rad/background_profiles', description_key, &
           & 'background radial profile of the negative toroidal rotation &
           & gradient', &
           & ascii_fmt)
        ! FIXME fix this or warn properly
        call attach_metadata('background_uprim', &
           & 'diagnostic/diagnos_rad/background_profiles', comments_key, &
           & 'Currently a dirty hack as this wont work with multiple ion &
           & species!', ascii_fmt)

      end if

    endif

    ! deallocate the help arrays
    deallocate(den_G)
    deallocate(temp_G)
    deallocate(rln_G)
    deallocate(rlt_G)
    deallocate(uprim_G)
    deallocate(ecurx_G)

  end subroutine output_background_profiles

  !------------------------------------------------------------------------------
  ! Subroutine that writes the radial profiles
  ! Not optimized 
  !------------------------------------------------------------------------------
  subroutine output_radial_profile(lrad_tint,lrad_moment,lrad_field,lrad_entropy,lrad_kpar) 

    use control,        only : nlapar
    use geom,           only : ints, ffun, sgr
    use grid,           only : nx, n_x_grid, ns, nmod, nvpar, nmu
    use grid,           only : nsp, lsp, gs
    use grid,           only : proc_subset, number_of_species
    use dist,           only : iphi, iapar, fdisi, ifdis, fmaxwl
    use index_function, only : indx
    use mpiinterface,   only : root_processor, gather_array
    use mpiinterface,   only : mpireduce_sum_inplace
    use mpicomms,       only : COMM_X_NE, COMM_S_NE
    use mpicomms,       only : COMM_X_EQ, COMM_SP_EQ_X_EQ
    use components,     only : rhostar, de, dgrid
    use geom,           only : bn, efun
    use velocitygrid,   only : intvp, intmu, vpgr, mugr 
    use mode,           only : iyzero, krho
    use fields,         only : get_averaged_phi
    use matdat,         only : get_f_from_g
    use components,     only : signz, tmp, tgrid
    use constants,      only : ci1
    use global,         only : r_tiny, r_huge, PHI_FIELD
    use io,             only : append_chunk, xy_fmt, ascii_fmt
    use linear_terms,   only : drift
    use rotation,       only : coriolis, cf_drift
    use diagnos_generic, only : dfieldds,  dfielddx, parseval_correction
    ! loop integers
    integer :: ix, i, imod, j, k, is, imod_dum
    logical :: lrad_tint,lrad_moment,lrad_field,lrad_entropy,lrad_kpar

    ! real for non-global rhostar
    real :: dum_rhostar
    
    ! dummy variables
    real :: drift_x, drift_y, drift_z
    complex :: par(ns), rad(nx)

    ! real for the normalized velocity space integral of the Maxwellian
    ! the integral is calculated at each grid points (species, radial, field)
    !real :: intfm, intfm1, intfm2, intfm_sum

    ! matrix element and amplitude of the source
    !complex :: mat_elem 

    ! check on rhostar
    if (abs(rhostar) < r_tiny) then
      dum_rhostar = 1.
    else
      dum_rhostar = rhostar
    endif

    ! THIS IS A MESS !
    if(lrad_field)then
      call gather_reduce_write_field(iphi)
      if (nlapar) call gather_reduce_write_field(iapar)
    endif

    if(lrad_moment) then
      ! The density, toroidal velocity, and temperature perturbations
      imod = iyzero 
      loop_global_species: do is = 1, number_of_species

        if (proc_subset(0,0,0,0,is)) then
          ! calculate the density profile
          dum(:) = 0.E0
          do ix = 1, nx; do i = 1, ns 
            do j = 1, nmu; do k = 1, nvpar 
              dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,lsp(is))*intmu(j)*bn(ix,i) &
                 & * dum_rhostar * real(get_f_from_g(imod,ix,i,j,k,lsp(is),fdisi))
            end do; end do;
          end do; end do;

          call mpireduce_sum_inplace(dum,shape(dum),COMM_SP_EQ_X_EQ)

        end if

        ! then gather the array on the root of COMM_CART (not
        ! necessarily working on the species is)
        if(proc_subset(0,1,1,1,is) .or. root_processor) then
          call gather_array(dum_G, dum, mpi_dtype_radial_line, COMM_X_NE, .true., is <= nsp)
        end if

        ! write on file 
        if (root_processor) then 
          call append_chunk(i_delta_n(is), real(dum_G), xy_fmt, ascii_fmt)
        endif

        if (proc_subset(0,0,0,0,is)) then
        ! calculate the parallel velocity profile
          dum(:) = 0.0E0
          do ix = 1, nx; do i = 1, ns 
            do j = 1, nmu; do k = 1, nvpar 
              dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,lsp(is))*intmu(j)*bn(ix,i) &
                 & * dum_rhostar * real(get_f_from_g(imod,ix,i,j,k,lsp(is),fdisi)) &
                 & * vpgr(i,j,k,lsp(is))
            end do; end do;
          end do; end do;

          call mpireduce_sum_inplace(dum,shape(dum),COMM_SP_EQ_X_EQ)

        end if

        ! then gather the array on the root of COMM_CART (not
        ! necessarily working on the species is)
        if(proc_subset(0,1,1,1,is) .or. root_processor) then
          call gather_array(dum_G, dum, mpi_dtype_radial_line, COMM_X_NE, .true., is <= nsp)
        end if

        ! write on file 
        if (root_processor) then 
          call append_chunk(i_delta_v(is), real(dum_G), xy_fmt, ascii_fmt)
        endif

        if (proc_subset(0,0,0,0,is)) then
          ! calculate the energy profile
          dum(:) = 0.E0
          do ix = 1, nx ; do i = 1, ns 
            do j = 1, nmu; do k = 1, nvpar 
              dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,lsp(is))*intmu(j)*bn(ix,i) &
                 & * dum_rhostar * real(get_f_from_g(imod,ix,i,j,k,lsp(is),fdisi)) &
                 & * (vpgr(i,j,k,lsp(is))**2 + 2.E0*mugr(j)*bn(ix,i)) 
            end do; end do; 
          end do; end do;

          call mpireduce_sum_inplace(dum,shape(dum),COMM_SP_EQ_X_EQ)

        end if

        ! then gather the array on the root of COMM_CART (not
        ! necessarily working on the species is)
        if(proc_subset(0,1,1,1,is) .or. root_processor) then
          call gather_array(dum_G, dum, mpi_dtype_radial_line, COMM_X_NE, .true., is <= nsp)
        end if

        ! write on file from single processor containing this species
        if (root_processor) then 
          call append_chunk(i_delta_e(is), real(dum_G), xy_fmt, ascii_fmt)
        endif
        
        
        if (proc_subset(0,0,0,0,is)) then
          ! calculate the temperature profile, i.e., the energy moment * (2./3.)
          ! with the density correction - 1.*delta_n
          ! The prefactors containting de, dgrid, tmp, tgrid consider the 
          ! normalization with the grid density in the global case
          dum(:) = 0.E0
          do ix = 1, nx ; do i = 1, ns 
            do j = 1, nmu; do k = 1, nvpar 
              dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,lsp(is))*intmu(j)*bn(ix,i)* &
                 & dum_rhostar * real(get_f_from_g(imod,ix,i,j,k,lsp(is),fdisi)) &
                 & * ((vpgr(i,j,k,lsp(is))**2 + 2.E0*mugr(j)*bn(ix,i)) * 2.E0/3.E0 &
                 & - 1.E0*tmp(ix,lsp(is))/(tgrid(lsp(is))))*dgrid(lsp(is))/de(ix,lsp(is))
            end do; end do; 
          end do; end do;

          call mpireduce_sum_inplace(dum,shape(dum),COMM_SP_EQ_X_EQ)

        end if

        ! then gather the array on the root of COMM_CART (not
        ! necessarily working on the species is)
        if(proc_subset(0,1,1,1,is) .or. root_processor) then
          call gather_array(dum_G, dum, mpi_dtype_radial_line, COMM_X_NE, .true., is <= nsp)
        end if

        ! write on file from single processor containing this species
        if (root_processor) then 
          call append_chunk(i_delta_t(is), real(dum_G), xy_fmt, ascii_fmt)
        endif

      end do loop_global_species
    end if

    ! turbulence intensity
    if(lrad_tint)then
      if(proc_subset(0,0,1,1,1)) then

        dum(:) = 0.
        do imod=1, nmod
          if (imod .ne. iyzero) then
            do ix = 1, nx
              do i = 1, ns
                dum(ix) = dum(ix) + &
                   & ints(i)*abs(fdisi(indx(iphi,imod,ix,i)))**2 &
                   & * parseval_correction(imod)
              end do
            end do
          end if
        end do

        ! finish the flux surface average
        call mpireduce_sum_inplace(dum,shape(dum),COMM_S_NE)

        if(proc_subset(0,1,1,1,1)) then
          ! then gather the array (no intercomm gather necessary)
          call gather_array(dum_G, n_x_grid, dum, nx, COMM_X_NE)
        end if

        ! write on file 
        if (root_processor) then
          call append_chunk(i_f_int, dum_G, xy_fmt, ascii_fmt)
        endif

        ! the zonal mode
        dum(:) = 0
        imod = iyzero
        do ix = 1, nx
          do i = 1, ns 
            dum(ix) = dum(ix) + & 
               & ints(i)*abs(fdisi(indx(iphi,imod,ix,i)))**2
          end do
        end do

        call mpireduce_sum_inplace(dum, shape(dum),COMM_S_NE)

        if(proc_subset(0,1,1,1,1)) then
          ! then gather the array (no intercomm gather necessary)
          call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
        end if

        ! write on file 
        if (root_processor) then 
          call append_chunk(i_zf_int, dum_G, xy_fmt, ascii_fmt)
        endif
      end if
    endif

    ! radial profile of various spectrally averaged quantities useful
    ! in determining the asymmetry in the mode structure.
    if(lrad_kpar)then
      if(proc_subset(0,0,1,1,1)) then

        dum(:) = 0. 
        dum2(:) = 0.
        cdum(:)=0.
        cdum2(:)=0.
        cdum3(:)=0.
        
        do ix = 1, nx
          do imod=1, nmod
            if(imod .ne. iyzero) then
              !Get parallel derivative of phi, inputs 4 and 5 are simply
              !1 because they are ignored for non gyroaveraged fields.
              call dfieldds(PHI_FIELD,imod,ix,1,1,par)

              do i = 1, ns
                call dfielddx(PHI_FIELD,imod,i,1,1,rad)

                dum(ix) = dum(ix) + &
                   & ints(i)*abs(fdisi(indx(iphi,imod,ix,i)))**2 &
                   & * parseval_correction(imod)

                dum2(ix) = dum2(ix) + &
                   & ints(i)*sgr(i)*abs(fdisi(indx(iphi,imod,ix,i)))**2 &
                   & * parseval_correction(imod)

                cdum(ix) = cdum(ix) + &
                   & ffun(ix,i)*par(i)*conjg(fdisi(indx(iphi,imod,ix,i)))* &
                   & ints(i) * parseval_correction(imod)

                cdum2(ix) = cdum2(ix) + &
                   & (ffun(ix,i)*par(i)*conjg(fdisi(indx(iphi,imod,ix,i))))**2 &
                   & *ints(i)* parseval_correction(imod)

                cdum3(ix) = cdum3(ix) + &
                   & rad(ix)*conjg(fdisi(indx(iphi,imod,ix,i)))* &
                   & ints(i)* parseval_correction(imod)

              end do
            end if
          end do
        end do

        ! finish the flux surface average
        call mpireduce_sum_inplace(dum,shape(dum),COMM_S_NE)
        call mpireduce_sum_inplace(dum2,shape(dum2),COMM_S_NE)
        call mpireduce_sum_inplace(cdum,shape(cdum),COMM_S_NE)
        call mpireduce_sum_inplace(cdum2,shape(cdum2),COMM_S_NE)
        call mpireduce_sum_inplace(cdum3,shape(cdum3),COMM_S_NE)

        do ix = 1, nx
          if(abs(dum(ix)) > r_tiny) then
            cdum(ix) = cdum(ix)/dum(ix)
            cdum2(ix) = cdum2(ix)/dum(ix)
            cdum3(ix) = cdum3(ix)/dum(ix)
            dum2(ix) = dum2(ix)/dum(ix)
          else
            cdum(ix) = r_huge
            cdum2(ix) = r_huge
            cdum3(ix) = r_huge
            dum2(ix) = r_huge
          end if
        end do

        if(proc_subset(0,1,1,1,1)) then
          ! then gather the array (no intercomm gather necessary)
          call gather_array(cdum_G, n_x_grid, cdum, nx, COMM_X_NE)
        end if

        ! write on file 
        if (root_processor) then
          call append_chunk(i_kpar, aimag(cdum_G), xy_fmt, ascii_fmt)
        endif
        
        if(proc_subset(0,1,1,1,1)) then
          ! then gather the array (no intercomm gather necessary)
          call gather_array(cdum_G, n_x_grid, cdum2, nx, COMM_X_NE)
        end if
        
        ! write on file 
        if (root_processor) then
          call append_chunk(i_kpar2, aimag(cdum_G), xy_fmt, ascii_fmt)
        endif

        if(proc_subset(0,1,1,1,1)) then
          ! then gather the array (no intercomm gather necessary)
          call gather_array(dum_G, n_x_grid, dum2, nx, COMM_X_NE)
        end if
        
        ! write on file 
        if (root_processor) then
          call append_chunk(i_theta0, dum_G, xy_fmt, ascii_fmt)
        endif
        
        if(proc_subset(0,1,1,1,1)) then
          ! then gather the array (no intercomm gather necessary)
          call gather_array(cdum_G, n_x_grid, cdum3, nx, COMM_X_NE)
        end if
        
        ! write on file 
        if (root_processor) then
          call append_chunk(i_kr, aimag(cdum_G), xy_fmt, ascii_fmt)
        endif

        ! write on file 
        if (root_processor) then
          call append_chunk(i_kr2, real(cdum_G), xy_fmt, ascii_fmt)
        endif

      end if
    endif

    ! entropy profile
    if (lrad_entropy) then

      ! entropy particles in the perturbations
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                   & 2.E0*abs(get_f_from_g(imod,ix,i,j,k,is,fdisi))**2/(2.E0*fmaxwl(ix,i,j,k,is)) )
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum), COMM_X_EQ)
      if(proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_f_entropy, dum_G, xy_fmt, ascii_fmt)
      endif

      ! entropy particles in the zonal mode
      imod = iyzero
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
               & abs(get_f_from_g(imod,ix,i,j,k,is,fdisi))**2/(2.E0*fmaxwl(ix,i,j,k,is)) )
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if(proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy, dum_G, xy_fmt, ascii_fmt)
      endif

      ! entropy fileds in the perturbations
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                   & 2.E0*(signz(is)**2*fmaxwl(ix,i,j,k,is))/(2.E0*tmp(ix,is)**2)* &
                   & (abs(fdisi(indx(iphi,imod,ix,i)))**2- &
                   & abs(get_averaged_phi(imod,ix,i,j,is,fdisi))**2))
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx, COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_f_entropy_fields, dum_G, xy_fmt, ascii_fmt)
      endif

      ! entropy fields in the zonal mode
      imod = iyzero
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
               & (signz(is)**2*fmaxwl(ix,i,j,k,is))/(2.E0*tmp(ix,is)**2)* &
               & (abs(fdisi(indx(iphi,imod,ix,i)))**2- &
               & abs(get_averaged_phi(imod,ix,i,j,is,fdisi))**2))
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum,nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy_fields, dum_G, xy_fmt, ascii_fmt)
      endif

      ! entropy flux ExB particles zonal mode real part
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*2.E0*real( &
                      & (1.0/fmaxwl(ix,i,j,k,is))* &
                      & get_f_from_g(iyzero,ix,i,j,k,is,fdisi)*get_f_from_g(imod,ix,i,j,k,is,fdisi)* &
                      & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy_flux_re, dum_G, xy_fmt, ascii_fmt)
      endif


      ! entropy flux ExB particles zonal mode imag part
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*2.E0*aimag( &
                      & (1.0/fmaxwl(ix,i,j,k,is))* &
                      & get_f_from_g(iyzero,ix,i,j,k,is,fdisi)*get_f_from_g(imod,ix,i,j,k,is,fdisi)* &
                      & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy_flux_im, dum_G, xy_fmt, ascii_fmt)
      endif


      ! entropy flux ExB fields zonal mode real part
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*2.E0*real( &
                      & (signz(is)/tmp(ix,is))* &
                      & get_averaged_phi(iyzero,ix,i,j,is,fdisi)*get_f_from_g(imod,ix,i,j,k,is,fdisi)* &
                      & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy_flux_fields_re, dum_G, xy_fmt, ascii_fmt)
      endif


      ! entropy flux ExB fields zonal mode imag part
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*2.E0*aimag( &
                      & (signz(is)/tmp(ix,is))* &
                      & get_averaged_phi(iyzero,ix,i,j,is,fdisi)*get_f_from_g(imod,ix,i,j,k,is,fdisi)* &
                      & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy_flux_fields_im, dum_G, xy_fmt, ascii_fmt)
      endif


      ! entropy flux ExB perturbations
      cdum(:) = (0.0, 0.0)
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=-nmod, nmod
              do imod_dum=-nmod, nmod
                if (imod .ne. 0 .AND. imod_dum .ne. 0) then
                  if (imod .ne. iyzero .AND. imod_dum .ne. iyzero .AND. &
                     & imod .ne. -iyzero .AND. imod_dum .ne. -iyzero) then
                    if (abs(mode_difference(imod,imod_dum)) <= nmod .AND. imod .ne. imod_dum) then
                      if (imod > 0 .AND. imod_dum > 0 .AND. imod-imod_dum > 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (get_f_from_g(imod_dum,ix,i,j,k,is,fdisi) / &
                           & (2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)* &
                           & get_averaged_phi(imod_dum,ix,i,j,is,fdisi) )* &
                           & get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi)* &
                           & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))
                        
                      else if (imod > 0 .AND. imod_dum < 0 .AND. imod-imod_dum > 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (conjg(get_f_from_g(-imod_dum,ix,i,j,k,is,fdisi))/ &
                           & (2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)* &
                           & conjg(get_averaged_phi(-imod_dum,ix,i,j,is,fdisi)) )* &
                           & get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi)* &
                           & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))

                      else if (imod > 0 .AND. imod_dum > 0 .AND. imod-imod_dum < 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (get_f_from_g(imod_dum,ix,i,j,k,is,fdisi) / &
                           & (2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)*get_averaged_phi(imod_dum,ix,i,j,is,fdisi) )* &
                           & conjg(get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi))* &
                           & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))

                      else if (imod > 0 .AND. imod_dum < 0 .AND. imod-imod_dum < 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (conjg(get_f_from_g(-imod_dum,ix,i,j,k,is,fdisi))/(2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)*conjg(get_averaged_phi(-imod_dum,ix,i,j,is,fdisi)) )* &
                           & conjg(get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi))* &
                           & conjg(ci1*krho(imod)*efun(ix,i,2,1)*get_averaged_phi(imod,ix,i,j,is,fdisi)))

                      else if (imod < 0 .AND. imod_dum > 0 .AND. imod-imod_dum > 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (get_f_from_g(imod_dum,ix,i,j,k,is,fdisi)/(2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)*get_averaged_phi(imod_dum,ix,i,j,is,fdisi) )* &
                           & get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi)* &
                           & conjg(-ci1*krho(-imod)*efun(ix,i,2,1)*conjg(get_averaged_phi(-imod,ix,i,j,is,fdisi))))

                      else if (imod < 0 .AND. imod_dum < 0 .AND. imod-imod_dum > 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (conjg(get_f_from_g(-imod_dum,ix,i,j,k,is,fdisi))/(2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)*conjg(get_averaged_phi(-imod_dum,ix,i,j,is,fdisi)) )* &
                           & get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi)* &
                           & conjg(-ci1*krho(-imod)*efun(ix,i,2,1)*conjg(get_averaged_phi(-imod,ix,i,j,is,fdisi))))

                      else if (imod < 0 .AND. imod_dum > 0 .AND. imod-imod_dum < 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (get_f_from_g(imod_dum,ix,i,j,k,is,fdisi)/(2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)*get_averaged_phi(imod_dum,ix,i,j,is,fdisi) )* &
                           & conjg(get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi))* &
                           & conjg(-ci1*krho(-imod)*efun(ix,i,2,1)*conjg(get_averaged_phi(-imod,ix,i,j,is,fdisi))))

                      else if (imod < 0 .AND. imod_dum < 0 .AND. imod-imod_dum < 0) then
                        cdum(ix) = cdum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                           & (conjg(get_f_from_g(-imod_dum,ix,i,j,k,is,fdisi))/(2.E0*fmaxwl(ix,i,j,k,is)) + &
                           & signz(is)/tmp(ix,is)*conjg(get_averaged_phi(-imod_dum,ix,i,j,is,fdisi)) )* &
                           & conjg(get_f_from_g(abs(mode_difference(imod,imod_dum)),ix,i,j,k,is,fdisi))* &
                           & conjg(-ci1*krho(-imod)*efun(ix,i,2,1)*conjg(get_averaged_phi(-imod,ix,i,j,is,fdisi))))
                      endif
                    endif
                  endif
                endif
              end do;
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(cdum, shape(cdum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        call gather_array(cdum_G, n_x_grid, cdum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_f_entropy_flux, real(cdum_G), xy_fmt, ascii_fmt)
      endif


      ! entropy flux drift particles zonal mode
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
            dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                    & abs(get_f_from_g(iyzero,ix,i,j,k,is,fdisi))**2/(2.E0*fmaxwl(ix,i,j,k,is))* &
                    & drift_x)
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy_flux_drift, dum_G, xy_fmt, ascii_fmt)
      endif


      ! entropy flux drift fields zonal mode
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
            dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                    & (signz(is)/tmp(ix,is))**2*fmaxwl(ix,i,j,k,is)* & 
                    & abs(get_averaged_phi(iyzero,ix,i,j,is,fdisi))**2/2.E0* &
                    & drift_x)
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zf_entropy_flux_fields_drift, dum_G, xy_fmt, ascii_fmt)
      endif


      ! entropy flux drift particles perturbations
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                        & 2.E0*abs(get_f_from_g(imod,ix,i,j,k,is,fdisi))**2/(2.E0*fmaxwl(ix,i,j,k,is))* &
                        & drift_x)
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_f_entropy_flux_drift, dum_G, xy_fmt, ascii_fmt)
      endif


      ! entropy flux drift fields perturbations
      dum(:) = 0.
      do is = 1, nsp
        do ix = 1, nx ; do i = 1, ns 
          do j = 1, nmu; do k = 1, nvpar
            do imod=1, nmod
              if (imod .ne. iyzero) then
                call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
                dum(ix) = dum(ix) + ints(i)*intvp(i,j,k,is)*intmu(j)*bn(ix,i)*( &
                        & 2.E0*(signz(is)/tmp(ix,is))**2*fmaxwl(ix,i,j,k,is)* & 
                        & abs(get_averaged_phi(imod,ix,i,j,is,fdisi))**2/2.E0* &
                        & drift_x)
              endif
            end do;
          end do; end do; 
        end do; end do;
      end do;

      call mpireduce_sum_inplace(dum, shape(dum),COMM_X_EQ)
      if (proc_subset(0,1,1,1,1)) then
        ! then gather the array
        call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)
      end if

      ! write on file 
      if (root_processor) then 
        call append_chunk(i_f_entropy_flux_fields_drift, dum_G, xy_fmt, ascii_fmt)
      endif

    endif ! lrad_entropy

  end subroutine output_radial_profile

  !--------------------------------------------------------------------------
  !> given two binormal mode indices in the interval [-nmod ... -1 1 ... nmod ]
  !> (nota bene the zero is not allowed)
  !> this function returns the difference (imod1 - imod2), corrects the
  !> returned index so that it fits to the fact that the zero mode has index 1
  !--------------------------------------------------------------------------
  pure function mode_difference(imod1, imod2)
    integer, intent(in) :: imod1, imod2
    integer :: mode_difference

    integer :: i1, i2
    
    i1 = imod1 - imod1/abs(imod1)
    i2 = imod2 - imod2/abs(imod2)
    mode_difference = i1 - i2
    if(mode_difference == 0) then
      mode_difference = 1
    else
      mode_difference = mode_difference + mode_difference/abs(mode_difference)
    end if

  end function mode_difference

  subroutine gather_reduce_write_field(id)
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use mpicomms, only : COMM_X_NE, COMM_S_NE
    use dist, only : fdisi
    use index_function, only : indx
    use grid, only : nx, n_x_grid, ns, nmod, parallel_s
    use geom, only : ints
    use mpiinterface, only : root_processor, gather_array, mpiallreduce_sum

    integer :: id
    integer :: imod, ix, i

    ! loop over binormal modes 
    do imod = 1, nmod 

      dum(:) = 0. 
      do ix = 1, nx 
        do i = 1, ns 
          dum(ix) = dum(ix) + ints(i)*abs(fdisi(indx(id,imod,ix,i))) 
        end do
      end do

      ! sum for global s 
      if (parallel_s) then 
        call mpiallreduce_sum(dum,dum_buf,nx,COMM_S_NE)
        dum(:) = dum_buf(:)
      endif

      ! then gather the array
      call gather_array(dum_G, n_x_grid, dum, nx,COMM_X_NE)

      if (root_processor) then 
        call append_chunk(i_fields, dum_G, xy_fmt, ascii_fmt)
      endif

    end do ! binormal

  end subroutine gather_reduce_write_field

  !--------------------------------------------------------------------------
  !> Close the files and deallocates the arrays
  !--------------------------------------------------------------------------
  subroutine finalize
    use io, only : close_lu, ascii_fmt
    use grid, only : number_of_species
    integer :: is
    ! as this diagnostic opens many logical units, depending on
    ! various switches, it is easier to rely on the io module to close
    ! them.

    if(lrad_moment) then
      do is = 1, number_of_species
        call close_lu(i_delta_n(is), ascii_fmt)
        call close_lu(i_delta_v(is), ascii_fmt)
        call close_lu(i_delta_t(is), ascii_fmt)
        call close_lu(i_delta_e(is), ascii_fmt)
      end do
      if(allocated(i_delta_n)) deallocate(i_delta_n)
      if(allocated(i_delta_v)) deallocate(i_delta_v)
      if(allocated(i_delta_t)) deallocate(i_delta_t)
      if(allocated(i_delta_e)) deallocate(i_delta_e)
    end if
    
    if(allocated(rad_entr)) deallocate(rad_entr)

    ! deallocate the help arrays 
    if(allocated(dum_G)) deallocate(dum_G)
    if(allocated(dum)) deallocate(dum) 
    if(allocated(dum2)) deallocate(dum2) 
    if(allocated(dum_buf)) deallocate(dum_buf)
    if(allocated(cdum)) deallocate(cdum)
    if(allocated(cdum2)) deallocate(cdum2)
    if(allocated(cdum3)) deallocate(cdum3)
    if(allocated(cdum_G)) deallocate(cdum_G)
    
  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output()
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use mpiinterface, only : root_processor
    use non_linear_terms, only : entropy_radial

    ! calculate the entropy flux 
    if (lradial_entropy) call entropy_radial(rad_entr)

    ! The radial profiles of phi / delta_n / delta_v / delta_t + background 
    call output_radial_profile(lrad_tint,lrad_moment, &
       & lrad_field,lrad_entropy,lrad_kpar)

    if(root_processor) then
      if (lradial_entropy) then 
        call append_chunk(i_radial_entropy, &
           & rad_entr(:,1), xy_fmt, ascii_fmt)

        call append_chunk(i_entropy_flux, &
           & rad_entr(:,2), xy_fmt, ascii_fmt)
      endif
    end if

  end subroutine output


end module diagnos_rad 
