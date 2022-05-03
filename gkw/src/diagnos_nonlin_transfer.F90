!------------------------------------------------------------------------------
!>
!> This diagnostic calculates the energy transfer by nonlinear triad
!> interactions.
!------------------------------------------------------------------------------
module diagnos_nonlin_transfer

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: output

  public :: lnonlin_transfer, lnonlin_transfer_fsa, nonlin_transfer_interval

  ! namelist switch
  !> True if the nonlin_transfer diagnostic is switched on which
  !> outputs detailed data on nonlinear interaction
  logical, save :: lnonlin_transfer
  logical, save :: lnonlin_transfer_fsa

  integer, save :: lun_time
  integer, save :: lun_distfluct
  integer, save :: lun_tempfluct

  integer, save :: nonlin_transfer_interval

  !> array for the data of the nonlin_transfer diagnostic
  real, save, allocatable :: nonlin_transfer(:,:,:,:)
  real, save, allocatable :: nonlin_transfer_neg(:,:,:,:)

  !> indices of the modes, used to iterate over possible triad interactions.
  integer, parameter :: NONLIN_TRANSFER_PLUSPLUS = 1
  integer, parameter :: NONLIN_TRANSFER_PLUSMINUS = 2
  integer, parameter :: NONLIN_TRANSFER_MINUSMINUS = 3
  integer, parameter :: NONLIN_TRANSFER_MINUSPLUS = 4

  !> the range of MPI tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    lnonlin_transfer =  .false.
    lnonlin_transfer_fsa =  .false.
    nonlin_transfer_interval = 1
  end subroutine set_default_nml_values


  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    
    call mpibcast(lnonlin_transfer,1)
    call mpibcast(lnonlin_transfer_fsa,1)
    call mpibcast(nonlin_transfer_interval,1)
    
  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use general, only : gkw_warn
    use control, only : non_linear, nlapar, nlbpar, flux_tube
    
    if(.not.lnonlin_transfer) return

    if (.not.non_linear) then
      call gkw_warn('nonlinear transfer diagnostic can only be used for &
         & nonlinear runs')
      lnonlin_transfer = .false.
    endif

    if (.not.flux_tube) then
      call gkw_warn('nonlinear transfer diagnostic shows diffs between&
      & serial/parallel data for global runs and entropy balance is &
         & derived for fluxtube system anyway, thus this diagnostic is&
         & disabled')
      lnonlin_transfer = .false.
    endif
    
    if (nonlin_transfer_interval < 1) then
      call gkw_warn('nonlin_transfer_interval: only values >= 1 make sense. &
         & Value is set to 1.')
      nonlin_transfer_interval = 1
    endif

    if(nlapar .or. nlbpar) then
      call gkw_warn('nonlinear transfer diagnostic implements only &
         & electrostatic contributions at the moment!')
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use control, only : spectral_radius
    use mpiinterface, only : root_processor, register_tag_range
    use io, only : open_real_lu, open_real_lu, binary_fmt, xy_fmt
    use io, only : attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use io, only : output_array, ascii_fmt
    use diagnos_generic, only : attach_metadata_grid, TEMPERATURE_MOMENT_TOKEN
    use diagnos_generic, only : DISTRIBUTION_TOKEN, LOCAL_DATA, X_GHOSTCELLS
    use diagnos_generic, only : N_MOMENTS
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, DISTRIBUTION
    use grid, only : nmod, nxmod, n_x_grid, n_y_grid
    use general, only : gkw_abort
    logical, intent(inout) :: requirements(:,:)
    integer, allocatable :: tmp(:,:)
    integer :: ierr, prod_case, mode1, mode2, mode_flat, ix
    integer :: iyzero, iy
    character(len=40) :: dsetname

    if (.not.lnonlin_transfer) return

    call register_tag_range(2 * (N_MOMENTS + 1), &
       & tag_range_start, tag_range_end_inkl)

    requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
    requirements(PHI_GA_FIELD,X_GHOSTCELLS) = .true.
    requirements(APAR_GA_FIELD,LOCAL_DATA) = .true.
    requirements(DISTRIBUTION,X_GHOSTCELLS) = .true.

    if(root_processor) then

      ! open logical units
      if(lnonlin_transfer_fsa) then
        dsetname = 'nonlin_transfer_fsa_'
      else
        dsetname = 'nonlin_transfer_'
      end if
      if(spectral_radius) then
        call open_real_lu(trim(dsetname)//trim(TEMPERATURE_MOMENT_TOKEN), &
           & 'diagnostic/diagnos_nonlin_transfer', &
           & shape(nonlin_transfer), &
           & binary_fmt, lun_tempfluct)
        call attach_metadata_grid(lun_tempfluct, 'nonlin_transfer_time', 'y_pair', &
           & 'x_pair', 'x_triads','y_triads', binary_fmt)
        call attach_metadata(lun_tempfluct, description_key, &
           & 'The nonlin_transfer, computed from momenta of the distribution &
           & function.', &
           & binary_fmt)

        call open_real_lu(trim(dsetname)//trim(DISTRIBUTION_TOKEN), &
           & 'diagnostic/diagnos_nonlin_transfer', &
           & shape(nonlin_transfer), &
           & binary_fmt, lun_distfluct)
        call attach_metadata_grid(lun_distfluct, 'nonlin_transfer_time', 'y_pair', &
           & 'x_pair', 'x_triads','y_triads', binary_fmt)
        call attach_metadata(lun_distfluct, description_key, &
           & 'The nonlin_transfer, computed from the distribution &
           & function.', &
           & binary_fmt)
      else
        call open_real_lu(trim(dsetname)//trim(TEMPERATURE_MOMENT_TOKEN), &
           & 'diagnostic/diagnos_nonlin_transfer', &
           & shape(nonlin_transfer(:,:,1,:)), &
           & binary_fmt, lun_tempfluct)
        call attach_metadata_grid(lun_tempfluct, 'nonlin_transfer_time', 'y_pair', &
           & 'xgr', 'y_triads', binary_fmt)
        call attach_metadata(lun_tempfluct, description_key, &
           & 'The nonlin_transfer, computed from momenta of the distribution &
           & function.', &
           & binary_fmt)

        call open_real_lu(trim(dsetname)//trim(DISTRIBUTION_TOKEN), &
           & 'diagnostic/diagnos_nonlin_transfer', &
           & shape(nonlin_transfer(:,:,1,:)), &
           & binary_fmt, lun_distfluct)
        call attach_metadata_grid(lun_distfluct, 'nonlin_transfer_time', 'y_pair', &
           & 'xgr', 'y_triads', binary_fmt)
        call attach_metadata(lun_distfluct, description_key, &
           & 'The nonlin_transfer, computed from the distribution &
           & function.', &
           & binary_fmt)
      end if

      
      call open_real_lu('nonlin_transfer_time', &
         & 'diagnostic/diagnos_nonlin_transfer', &
         & (/ 1 /), &
         & ascii_fmt, lun_time)
      call attach_metadata_grid(lun_time, 'nonlin_transfer_time', ascii_fmt)
      call attach_metadata(lun_time, description_key, &
         & 'The timestamps of the nonlin_transfer data.', &
         & ascii_fmt)
      call attach_metadata(lun_time, comments_key, &
         & 'In principle these values can be computed from /grid/time and &
         & nonlin_transfer_interval, but this can be a pain.', &
         & ascii_fmt)

      ! output helper arrays:
      allocate(tmp(nmod**2,4),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate temporary array')
      tmp = 0
      do mode2 = 1,nmod
        do mode1 = 1,nmod
          mode_flat = mode1 + (mode2-1)*nmod

          iyzero = nmod

          prod_case = NONLIN_TRANSFER_PLUSPLUS
          iy = iyzero + (mode1-1) + (mode2-1)
          if(iy < 1 .or. iy > n_y_grid) iy = 0
          tmp(mode_flat,prod_case) = iy

          prod_case = NONLIN_TRANSFER_PLUSMINUS
          iy = iyzero + (mode1-1) - (mode2-1)
          if(iy < 1 .or. iy > n_y_grid) iy = 0
          tmp(mode_flat,prod_case) = iy

          prod_case = NONLIN_TRANSFER_MINUSMINUS
          iy = iyzero - (mode1-1) - (mode2-1)
          if(iy < 1 .or. iy > n_y_grid) iy = 0
          tmp(mode_flat,prod_case) = iy

          prod_case = NONLIN_TRANSFER_MINUSPLUS
          iy = iyzero - (mode1-1) + (mode2-1)
          if(iy < 1 .or. iy > n_y_grid) iy = 0
          tmp(mode_flat,prod_case) = iy
        end do
      end do
      call output_array('y_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & real(tmp), 'F', xy_fmt, ascii_fmt)
      call attach_metadata('y_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata('y_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & description_key, &
         & 'This translates the penultimate dimension of nonlin_transfer into&
         & binormal modes (in the range 1:nmod), for each mode pair.', &
         & ascii_fmt)
      call attach_metadata('y_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & comments_key, not_avail, ascii_fmt)
      deallocate(tmp)


      allocate(tmp(nxmod**2,4),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate temporary array')
      tmp = 0
      do mode2 = 1,nxmod
        do mode1 = 1,nxmod
          ! mode1 is always >= mode2
          do ix = 1,n_x_grid
            prod_case = nonlin_transfer_get_ixsparse(ix,mode1,mode2)
            if(prod_case /= 0) then
              mode_flat = mode1 + (mode2-1)*nxmod
              tmp(mode_flat,prod_case) = ix
            end if
          end do
        end do
      end do
      call output_array('x_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & real(tmp), 'F', xy_fmt, ascii_fmt)
      call attach_metadata('x_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata('x_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & description_key, &
         & 'This translates the last dimension of nonlin_transfer into&
         & radial modes (in the range 1:n_x_grid), for each mode pair index.', &
         & ascii_fmt)
      call attach_metadata('x_triads', '/diagnostic/diagnos_nonlin_transfer', &
         & comments_key, not_avail, ascii_fmt)
      deallocate(tmp)


      allocate(tmp(2, nmod**2),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate temporary array')
      tmp = 0
      do mode2 = 1,nmod
        do mode1 = 1,nmod
          mode_flat = mode1 + (mode2-1)*nmod
          ! mode1 is always >= mode2
          tmp(1,mode_flat) = mode1
          tmp(2,mode_flat) = mode2
        end do
      end do
      call output_array('y_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & real(tmp), 'F', xy_fmt, ascii_fmt)
      call attach_metadata('y_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata('y_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & description_key, &
         & 'This translates the flat mode pair index into two mode indices in &
         & the range 1:nmod', &
         & ascii_fmt)
      call attach_metadata('y_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & comments_key, not_avail, ascii_fmt)
      deallocate(tmp)

      allocate(tmp(2, nxmod**2),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate temporary array')
      tmp = 0
      do mode2 = 1,nxmod
        do mode1 = 1,nxmod
          ! mode1 is always >= mode2
          mode_flat = mode1 + (mode2-1)*nxmod
          tmp(1,mode_flat) = mode1
          tmp(2,mode_flat) = mode2
        end do
      end do
      call output_array('x_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & real(tmp), 'F', xy_fmt, ascii_fmt)
      call attach_metadata('x_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata('x_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & description_key, &
         & 'This translates the flat mode pair index into two mode indices in &
         & the range 1:nxmod', &
         & ascii_fmt)
      call attach_metadata('x_pairs', '/diagnostic/diagnos_nonlin_transfer', &
         & comments_key, not_avail, ascii_fmt)
      deallocate(tmp)

    end if

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use control, only : spectral_radius
    use grid, only : nmod, nxmod, n_x_grid
    use general, only : gkw_abort

    integer :: ierr

    if(.not.lnonlin_transfer) return

    ! for binormal modes, there are 4 combinations
    ! (k1+k2), (-k1-k2),
    ! (-k1+k2), (k1-k2)
    if(spectral_radius) then
      ! and for radial modes the same.
      allocate(nonlin_transfer(nmod**2,nxmod**2, &
         & 4,4), &
         & stat = ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate nonlin_transfer')
      allocate(nonlin_transfer_neg(nmod**2,nxmod**2, &
         & 4,4), &
         & stat = ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate nonlin_transfer')
    else
      ! and for radial modes one has n_x_grid position space points
      allocate(nonlin_transfer(nmod**2,n_x_grid, &
         & 1,4), &
         & stat = ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate nonlin_transfer')
      allocate(nonlin_transfer_neg(nmod**2,n_x_grid, &
         & 1,4), &
         & stat = ierr)
      if (ierr /= 0) call gkw_abort('diagnos_nonlin_transfer :: could not &
         & allocate nonlin_transfer')
    end if

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : close_lu, binary_fmt, ascii_fmt
    use mpiinterface, only : root_processor

    if (.not.lnonlin_transfer) return

    if(allocated(nonlin_transfer)) deallocate(nonlin_transfer)
    if(allocated(nonlin_transfer_neg)) deallocate(nonlin_transfer_neg)

    ! be nice and close all logical units
    if(root_processor) then
      call close_lu(lun_tempfluct, binary_fmt)
      call close_lu(lun_distfluct, binary_fmt)

      call close_lu(lun_time, ascii_fmt)
    end if
  end subroutine finalize

  !--------------------------------------------------------------------
  !> The routine calc_largestep() is to be called repeatedly, after every
  !> large timestep, and here the diagnostic should calculate its
  !> quantities.
  !> 
  !> Splitting calculation and output is to be encouraged, particularly when
  !>  - the output is higher dimensional, and requires MPI or MPI-IO
  !>  - same data is output in different dimensional slices
  !--------------------------------------------------------------------
  subroutine calc_largestep()
    use mpiinterface, only : root_processor
    
    if (.not.lnonlin_transfer) return

    if(root_processor) then
      write (*,*) "diagnos_nonlin_transfer: Measuring nonlinear interaction &
       &(this can take a long time)..."
    end if
    
  end subroutine calc_largestep

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output()
    use control, only : time
    use mpiinterface, only : root_processor
    use control, only : nt_complete, itime
    use diagnos_generic, only : START_MOMENTS, TEMPERATURE_MOMENT
    use global, only : DISTRIBUTION
    use io, only : append_chunk, xy_fmt, ascii_fmt
    integer :: tag
    
    if (.not.lnonlin_transfer) return

    ! as this is very expensive, only do something every
    ! nonlin_transfer_interval'th large timestep
    if (modulo(nt_complete+itime,nonlin_transfer_interval) /= 0) return

    if(root_processor) then
      write(*,*) 'Nonlinear transfer diagnostic is computing, this may &
         &take a while...'
    end if

    tag = TEMPERATURE_MOMENT
    call calc_output_allsp(START_MOMENTS + TEMPERATURE_MOMENT, &
       & lun_tempfluct,tag)
    tag = 0
    call calc_output_allsp(DISTRIBUTION, lun_distfluct,tag)

    if(root_processor) then
      call append_chunk(lun_time, (/ time /), xy_fmt, ascii_fmt)
      write(*,*) 'Nonlinear transfer diagnostic is done.'
    end if
    
  end subroutine output

  !-----------------------------------------------------------------------------
  !> 
  !-----------------------------------------------------------------------------
  subroutine calc_output_allsp(iquantity, lun,tag)
    use mpiinterface, only : mpireduce_sum_inplace, root_processor, send_to_root
    use io, only : append_chunk, xy_fmt, ascii_fmt, binary_fmt
    use grid, only : number_of_species, proc_subset
    use mpicomms, only : COMM_X_EQ, COMM_CART, COMM_S_EQ_X_EQ, COMM_SP_NE
    use mpicomms, only : COMM_S_EQ, COMM_SP_NE_S_NE
    use global, only : DISTRIBUTION
    use diagnos_generic, only : xy_slice_ipar
    use control, only : spectral_radius
    integer, intent(in) :: iquantity, lun, tag
    integer :: gipar, isg, comm
    logical :: proc_is_involved
    if(lnonlin_transfer_fsa) then
      gipar = 0
    else
      gipar = xy_slice_ipar
    end if
    nonlin_transfer = 0.0
    nonlin_transfer_neg = 0.0
    
    do isg = 1, number_of_species
      if(spectral_radius) then
        if(iquantity == DISTRIBUTION) then
          call calc_nl_transfer_entr(iquantity,isg, gipar)
        else
          call calc_nl_transfer_moments(iquantity,isg, gipar)
        end if
      else
        if(iquantity == DISTRIBUTION) then
          call calc_nl_transfer_entr_nonspec(iquantity,isg, gipar)
        else
          ! not implemented yet
          return
        end if
      end if
    end do

    if(lnonlin_transfer_fsa) then
      if(iquantity == DISTRIBUTION) then
        if(spectral_radius) then
          ! complete the flux surface avg and the velspace integral,
          ! and the species sum
          ! (this is actually equal to COMM_CART, yes)
          comm = COMM_X_EQ
          proc_is_involved = proc_subset(0,0,0,0,0)
        else
          ! complete the flux surface avg and the velspace integral,
          ! and the species sum, and the x dependency
          comm = COMM_CART
          proc_is_involved = proc_subset(0,0,0,0,0)
        end if
      else
        if(spectral_radius) then
          ! complete the flux surface avg and the species sum
          comm = COMM_SP_NE_S_NE
          proc_is_involved = proc_subset(0,0,1,1,0)
        else
          !not implemented
          return
        end if
      end if
    else
      if(iquantity == DISTRIBUTION) then
        if(spectral_radius) then
          ! complete the velspace integral and the species sum
          comm = COMM_S_EQ_X_EQ
          proc_is_involved = proc_subset(0,gipar,0,0,0)
        else
          ! complete the velspace integral and the species sum and the
          ! x dependency
          comm = COMM_S_EQ
          proc_is_involved = proc_subset(0,gipar,0,0,0)
        end if
      else
        if(spectral_radius) then
          ! complete the species sum
          comm = COMM_SP_NE
          proc_is_involved = proc_subset(0,gipar,1,1,0)
        else
          !not implemented
          return
        end if
      end if
    end if
    if(proc_is_involved) then
      call mpireduce_sum_inplace(nonlin_transfer, &
         & shape(nonlin_transfer), comm)
      call mpireduce_sum_inplace(nonlin_transfer_neg, &
         & shape(nonlin_transfer_neg), comm)
    end if
    if(.not.lnonlin_transfer_fsa) then
      ! that means gipar /= 0
      
      ! if the transfer is not fsaed, but computed on a slice that is
      ! not on the root process, then send the result to the root.
      if((root_processor .and. .not. proc_subset(1,gipar,1,1,1)) &
         & .or. (proc_subset(1,gipar,1,1,1) .and. .not. root_processor)) then
        call send_to_root(nonlin_transfer, shape(nonlin_transfer), &
           & COMM_CART, tag*2)
        call send_to_root(nonlin_transfer_neg, shape(nonlin_transfer_neg), &
           & COMM_CART, tag*2+1)
      end if
    end if
    if(root_processor) then
      if(spectral_radius) then
        call append_chunk(lun, &
           & nonlin_transfer + nonlin_transfer_neg, &
           & xy_fmt, binary_fmt)
      else
        call append_chunk(lun, &
           & nonlin_transfer(:,:,1,:) + nonlin_transfer_neg(:,:,1,:), &
           & xy_fmt, binary_fmt)
      end if
    end if
  end subroutine calc_output_allsp


  !-----------------------------------------------------------------------------
  !> 
  !-----------------------------------------------------------------------------
  subroutine calc_nl_transfer_entr(iquantity, ispg, gipar)
    use dist,           only : fdisi, phi, apar, fmaxwl
    use control, only : nlapar
    use matdat,         only : get_f_from_g
    use components,     only : signz, tmp, vthrat
    use mode,           only : krho, kxrh, ixzero
    use geom,           only : efun, bn, ints
    use grid,           only : nmod, nxmod, n_x_grid, ns, nmu, nvpar, n_y_grid
    use grid,           only : gs, proc_subset, lsp
    use velocitygrid,   only : intmu, intvp, vpgr
    use constants,      only : ci1
    use functions,      only : besselj0_gkw
    use dist, only : ifdis
    use index_function, only : indx

    use general, only : gkw_abort

    use diagnos_cross_phase, only : get_field_or_moment
    integer, intent(in) :: iquantity
    !> global species index. This is used to access moment_buf and
    !> gyroavg fields.
    integer, intent(in) :: ispg
    !> global parallel index. If this is 0, then a flux surf avg will
    !> be computed.
    integer, intent(in) :: gipar
    complex :: cdum
    integer :: imod, ix, iy, jv, kt, is, ipar
    real    :: dum
    integer :: ix_sparse, iy_sparse, ymode_flat, xmode_flat
    !> loop boundaries, will cause loops or no loops
    integer :: nonlin_transfer_ymode1, nonlin_transfer_ymode2
    integer :: nonlin_transfer_xmode1, nonlin_transfer_xmode2
    integer :: xmode, ymode, iyzero
    integer :: ix_mode1, ix_mode2, iy_mode1, iy_mode2
    integer :: iy_sign, ix_sign, xmode_sign1, xmode_sign2, ymode_sign1, ymode_sign2
    integer :: ix_theothersign
    complex, dimension(n_y_grid,n_x_grid) :: v_exb_x, v_exb_y, grad_x_buf, grad_y_buf

    ! this is NOT the iyzero from module mode! This is the index of
    ! the ky=0 mode in a grid of n_y_grid modes, where the zeromode is
    ! in the centre (analogous to the kx grid we are working with).
    iyzero = nmod

    procs_on_gipar_ispg: if(proc_subset(0,gipar,0,0,ispg)) then

      is = lsp(ispg)
      loop_ns: do ipar = 1, ns

        ! setting nonlin_transfer_ipar=0 gives a flux-surface-avg of the transfer,
        ! while >= 1 gives transfer in the chosen slice.
        if(gipar /= 0 .and. gs(ipar) /= gipar) cycle

        loop_mu: do jv = 1, nmu
          loop_vpar: do kt = 1, nvpar
            ! fill buffers with the "complete" k-space fields. All
            ! values are stored explicitely, despite of hermitian
            ! symmetry, to make life easier later.
            do xmode = 1, nxmod
              do ymode = 1, nmod
                iy_sign = +1
                iy = iyzero + iy_sign*(ymode-1)

                ix_sign = +1
                ix = ixzero + ix_sign*(xmode-1)
                ! efun is actually
                ! symmetric and the diagnonal elements are
                ! zero, but for readability I write it out
                cdum = phi(ymode,ix,ipar)
                if(nlapar) then
                  cdum = cdum - 2.0 * vthrat(is) * vpgr(ipar,jv,kt,is)&
                   & * apar(ymode,ix,ipar)
                end if
                v_exb_x(iy,ix) = &
                   & (efun(ix,ipar,1,1)*ci1*kxrh(ix) + &
                   & efun(ix,ipar,2,1)*ci1*krho(ymode)) * &
                   & besselj0_gkw(ymode,ix,ipar,jv,is) * cdum
                v_exb_y(iy,ix) = &
                   & (efun(ix,ipar,1,2)*ci1*kxrh(ix) + &
                   & efun(ix,ipar,2,2)*ci1*krho(ymode)) * &
                   & besselj0_gkw(ymode,ix,ipar,jv,is) * cdum
                grad_x_buf(iy,ix) = ci1 * kxrh(ix) * &
                   & fdisi(indx(ifdis,ymode, &
                   & ix,ipar,jv,kt,is))
                grad_y_buf(iy,ix) = ci1 * krho(ymode) * &
                   & fdisi(indx(ifdis,ymode, &
                   & ix,ipar,jv,kt,is))

                if(xmode /= 1) then
                  ix_sign = -1
                  ix = ixzero + ix_sign*(xmode-1)
                  cdum = phi(ymode,ix,ipar)
                  if(nlapar) then
                    cdum = cdum - 2.0 * vthrat(is) * vpgr(ipar,jv,kt,is)&
                       & * apar(ymode,ix,ipar)
                  end if
                  v_exb_x(iy,ix) = &
                     & (efun(ix,ipar,1,1)*ci1*kxrh(ix) + &
                     & efun(ix,ipar,2,1)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix,ipar,jv,is) * cdum
                  v_exb_y(iy,ix) = &
                     & (efun(ix,ipar,1,2)*ci1*kxrh(ix) + &
                     & efun(ix,ipar,2,2)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix,ipar,jv,is) * cdum
                  grad_x_buf(iy,ix) = ci1 * kxrh(ix) * &
                     & fdisi(indx(ifdis,ymode, &
                     & ix,ipar,jv,kt,is))
                  grad_y_buf(iy,ix) = ci1 * krho(ymode) * &
                     & fdisi(indx(ifdis,ymode, &
                     & ix,ipar,jv,kt,is))
                end if

                if(ymode /= 1) then
                  !use hermitian symmetry:
                  ! f(-kx, -ky) = f*(kx,ky)
                  !and
                  ! f(kx, -ky)  = f*(-kx,ky)
                  iy_sign = -1
                  iy = iyzero + iy_sign*(ymode-1)

                  ix_sign = +1
                  ix = ixzero + ix_sign*(xmode-1)
                  ix_theothersign = ixzero - ix_sign*(xmode-1)
                  cdum = phi(ymode,ix_theothersign,ipar)
                  if(nlapar) then
                    cdum = cdum - 2.0 * vthrat(is) * vpgr(ipar,jv,kt,is)&
                       & * apar(ymode,ix_theothersign,ipar)
                  end if
                  v_exb_x(iy,ix) = &
                     & conjg((efun(ix_theothersign,ipar,1,1)*ci1*kxrh(ix_theothersign) + &
                     & efun(ix_theothersign,ipar,2,1)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * cdum)
                  v_exb_y(iy,ix) = &
                     & conjg((efun(ix_theothersign,ipar,1,2)*ci1*kxrh(ix_theothersign) + &
                     & efun(ix_theothersign,ipar,2,2)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * cdum)
                  grad_x_buf(iy,ix) = conjg(ci1 * kxrh(ix_theothersign) * &
                     & fdisi(indx(ifdis,ymode, &
                     & ix_theothersign,ipar,jv,kt,is)))
                  grad_y_buf(iy,ix) = conjg(ci1 * krho(ymode) * &
                     & fdisi(indx(ifdis,ymode, &
                     & ix_theothersign,ipar,jv,kt,is)))

                  if(xmode /= 1) then
                    ix_sign = -1
                    ix = ixzero + ix_sign*(xmode-1)
                    ix_theothersign = ixzero - ix_sign*(xmode-1)
                    cdum = phi(ymode,ix_theothersign,ipar)
                    if(nlapar) then
                      cdum = cdum - 2.0 * vthrat(is) * vpgr(ipar,jv,kt,is)&
                         & * apar(ymode,ix_theothersign,ipar)
                    end if
                    v_exb_x(iy,ix) = &
                       & conjg((efun(ix_theothersign,ipar,1,1)*ci1*kxrh(ix_theothersign) + &
                       & efun(ix_theothersign,ipar,2,1)*ci1*krho(ymode)) * &
                       & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * cdum)
                    v_exb_y(iy,ix) = &
                       & conjg((efun(ix_theothersign,ipar,1,2)*ci1*kxrh(ix_theothersign) + &
                       & efun(ix_theothersign,ipar,2,2)*ci1*krho(ymode)) * &
                       & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * cdum)
                    grad_x_buf(iy,ix) = conjg(ci1 * kxrh(ix_theothersign) * &
                       & fdisi(indx(ifdis,ymode, &
                       & ix_theothersign,ipar,jv,kt,is)))
                    grad_y_buf(iy,ix) = conjg(ci1 * krho(ymode) * &
                       & fdisi(indx(ifdis,ymode, &
                       & ix_theothersign,ipar,jv,kt,is)))
                  end if
                end if
              end do
            end do

            ! Measure the interaction of different modes

            ! iterate over all mode pairs: mode2 is the "from-mode", mode1
            ! is the "mediating mode" and the other one is the "to-mode"
            loop_ymode2: do nonlin_transfer_ymode2 = 1, nmod
              do ymode_sign2 = +1, -1, -2
                if(nonlin_transfer_ymode2 == 1 .and. ymode_sign2 == -1) then
                  cycle
                end if
                loop_ymode1: do nonlin_transfer_ymode1 = 1, nmod
                  do ymode_sign1 = +1, -1, -2
                    if(nonlin_transfer_ymode1 == 1 .and. ymode_sign1 == -1) then
                      cycle
                    end if
                    iy_mode1 = iyzero + ymode_sign1*(nonlin_transfer_ymode1-1)
                    iy_mode2 = iyzero + ymode_sign2*(nonlin_transfer_ymode2-1)
                    ymode_flat = nonlin_transfer_ymode1 + (nonlin_transfer_ymode2-1)*nmod

                    iy = iyzero + ymode_sign1*(nonlin_transfer_ymode1-1) &
                       & + ymode_sign2*(nonlin_transfer_ymode2-1)
                    if(iy < 1 .or. iy > n_y_grid) cycle

                    if(ymode_sign1 == +1 .and. ymode_sign2 == +1) then
                      iy_sparse = NONLIN_TRANSFER_PLUSPLUS
                    elseif(ymode_sign1 == +1 .and. ymode_sign2 == -1) then
                      iy_sparse = NONLIN_TRANSFER_PLUSMINUS
                    elseif(ymode_sign1 == -1 .and. ymode_sign2 == +1) then
                      iy_sparse = NONLIN_TRANSFER_MINUSPLUS
                    elseif(ymode_sign1 == -1 .and. ymode_sign2 == -1) then
                      iy_sparse = NONLIN_TRANSFER_MINUSMINUS
                    end if

                    loop_xmode2: do nonlin_transfer_xmode2 = 1, nxmod
                      do xmode_sign2 = +1, -1, -2
                        if(nonlin_transfer_xmode2 == 1 .and. xmode_sign2 == -1) then
                          cycle
                        end if
                        loop_xmode1: do nonlin_transfer_xmode1 = 1,nxmod
                          do xmode_sign1 = +1, -1, -2
                            if(nonlin_transfer_xmode1 == 1 .and. xmode_sign1 == -1) then
                              cycle
                            end if
                            ix_mode1 = ixzero + xmode_sign1*(nonlin_transfer_xmode1-1)
                            ix_mode2 = ixzero + xmode_sign2*(nonlin_transfer_xmode2-1)
                            xmode_flat = nonlin_transfer_xmode1 + (nonlin_transfer_xmode2-1)*nxmod

                            ix = ixzero + xmode_sign1*(nonlin_transfer_xmode1-1) &
                               & + xmode_sign2*(nonlin_transfer_xmode2-1)
                            ix_theothersign = ixzero - (ix-ixzero)
                            if(ix < 1 .or. ix > n_x_grid) cycle

                            if(xmode_sign1 == +1 .and. xmode_sign2 == +1) then
                              ix_sparse = NONLIN_TRANSFER_PLUSPLUS
                            elseif(xmode_sign1 == +1 .and. xmode_sign2 == -1) then
                              ix_sparse = NONLIN_TRANSFER_PLUSMINUS
                            elseif(xmode_sign1 == -1 .and. xmode_sign2 == +1) then
                              ix_sparse = NONLIN_TRANSFER_MINUSPLUS
                            elseif(xmode_sign1 == -1 .and. xmode_sign2 == -1) then
                              ix_sparse = NONLIN_TRANSFER_MINUSMINUS
                            end if

                            imod = abs(iy - iyzero) + 1

                            if(iy >= iyzero) then
                              cdum = (get_f_from_g(imod,ix,ipar,jv,kt,is,fdisi) &
                                 & / fmaxwl(ix, ipar, jv, kt, is) + signz(is)/tmp(ix,is) &
                                 & * besselj0_gkw(imod,ix,ipar,jv,is) * phi(imod,ix,ipar))
                            else
                              cdum = conjg(get_f_from_g(imod,ix_theothersign,ipar,jv,kt,is,fdisi) &
                                 & / fmaxwl(ix_theothersign, ipar, jv, kt, is) + signz(is)/tmp(ix_theothersign,is) &
                                 & * besselj0_gkw(imod,ix_theothersign,ipar,jv,is) * phi(imod,ix_theothersign,ipar))
                            end if

                            dum = - ints(ipar) * bn(ix,ipar) * intmu(jv) * intvp(ipar, jv, kt, is) &
                               & * real( &
                               & conjg(cdum) &
                               & * ( &
                               & v_exb_x(iy_mode1, ix_mode1)*&
                               & grad_x_buf(iy_mode2, ix_mode2) &
                               & + &
                               & v_exb_y(iy_mode1, ix_mode1)*&
                               & grad_y_buf(iy_mode2, ix_mode2) &
                               &) &
                               & )

                            if(dum > 0) then
                              nonlin_transfer(ymode_flat, &
                                 & xmode_flat, ix_sparse, iy_sparse) = &
                                 & nonlin_transfer(ymode_flat, &
                                 & xmode_flat, ix_sparse, iy_sparse) &
                                 & + dum
                            else
                              nonlin_transfer_neg(ymode_flat, &
                                 & xmode_flat, ix_sparse, iy_sparse) = &
                                 & nonlin_transfer_neg(ymode_flat, &
                                 & xmode_flat, ix_sparse, iy_sparse) &
                                 & + dum
                            end if
                          end do
                        end do loop_xmode1
                      end do
                    end do loop_xmode2
                  end do
                end do loop_ymode1
              end do
            end do loop_ymode2
          end do loop_vpar
        end do loop_mu

      end do loop_ns
    end if procs_on_gipar_ispg
  end subroutine calc_nl_transfer_entr

  !-----------------------------------------------------------------------------
  !> 
  !-----------------------------------------------------------------------------
  subroutine calc_nl_transfer_moments(iquantity, ispg, gipar)
    use dist,           only : phi
    use matdat,         only : get_f_from_g
    use mode,           only : krho, kxrh, ixzero
    use geom,           only : efun, bn, ints
    use grid,           only : nmod, nxmod, n_x_grid, ns, nmu, nvpar, nsp, n_y_grid
    use grid,           only : gs, proc_subset, lsp
    use velocitygrid,   only : intmu, intvp
    use constants,      only : ci1
    use functions,      only : besselj0_gkw
    use dist, only : ifdis
    use mpicomms, only : COMM_VPAR_NE_MU_NE
    use mpiinterface, only : mpireduce_sum_inplace

    use general, only : gkw_abort

    ! use diagnos_cross_phase, only : get_field_or_moment
    use diagnos_moments, only : get_4d_moment
    integer, intent(in) :: iquantity
    !> global species index. This is used to access moment_buf and
    !> gyroavg fields.
    integer, intent(in) :: ispg
    !> global parallel index. If this is 0, then a flux surf avg will
    !> be computed.
    integer, intent(in) :: gipar
    complex :: cdum
    integer :: imod, ix, iy, jv, kt, is, ipar
    real    :: dum, d3v
    integer :: ix_sparse, iy_sparse, ymode_flat, xmode_flat
    !> loop boundaries, will cause loops or no loops
    integer :: nonlin_transfer_ymode1, nonlin_transfer_ymode2
    integer :: nonlin_transfer_xmode1, nonlin_transfer_xmode2
    integer :: xmode, ymode, iyzero
    integer :: ix_mode1, ix_mode2, iy_mode1, iy_mode2
    integer :: iy_sign, ix_sign, xmode_sign1, xmode_sign2, ymode_sign1, ymode_sign2
    integer :: ix_theothersign
    !> this is local in species and local in parallel direction
    complex, dimension(nsp,nmod,ns,n_x_grid) :: buf
    complex, dimension(n_y_grid,n_x_grid) :: v_exb_x, v_exb_y, grad_x_buf, grad_y_buf

    ! this is NOT the iyzero from module mode! This is the index of
    ! the ky=0 mode in a grid of n_y_grid modes, where the zeromode is
    ! in the centre (analogous to the kx grid we are working with).
    iyzero = nmod

    procs_on_gipar_ispg: if(proc_subset(0,gipar,0,0,ispg)) then

      buf = get_4d_moment(iquantity)
      call mpireduce_sum_inplace(buf, shape(buf), COMM_VPAR_NE_MU_NE)


      is = lsp(ispg)
      loop_ns: do ipar = 1, ns

        ! setting nonlin_transfer_ipar=0 gives a flux-surface-avg of the transfer,
        ! while >= 1 gives transfer in the chosen slice.
        if(gipar /= 0 .and. gs(ipar) /= gipar) cycle

        v_exb_x = 0.0
        v_exb_y = 0.0

        ! fill buffers with the "complete" k-space fields. All
        ! values are stored explicitely, despite of hermitian
        ! symmetry, to make life easier later.
        do xmode = 1, nxmod
          do ymode = 1, nmod
            iy_sign = +1
            iy = iyzero + iy_sign*(ymode-1)

            ix_sign = +1
            ix = ixzero + ix_sign*(xmode-1)
            ! efun is actually
            ! symmetric and the diagnonal elements are
            ! zero, but for readability I write it out
            do jv = 1, nmu
              do kt = 1, nvpar
                d3v = bn(ix,ipar) * intmu(jv) * intvp(ipar, jv, kt, is)
                v_exb_x(iy,ix) = v_exb_x(iy,ix) + d3v * &
                   & (efun(ix,ipar,1,1)*ci1*kxrh(ix) + &
                   & efun(ix,ipar,2,1)*ci1*krho(ymode)) * &
                   & besselj0_gkw(ymode,ix,ipar,jv,is) * &
                   & phi(ymode,ix,ipar)
                v_exb_y(iy,ix) = v_exb_y(iy,ix) + d3v * &
                   & (efun(ix,ipar,1,2)*ci1*kxrh(ix) + &
                   & efun(ix,ipar,2,2)*ci1*krho(ymode)) * &
                   & besselj0_gkw(ymode,ix,ipar,jv,is) * &
                   & phi(ymode,ix,ipar)
              end do
            end do
            if(proc_subset(0,gipar,1,1,ispg)) then
              grad_x_buf(iy,ix) = ci1 * kxrh(ix) * &
                 & buf(is,ymode, ipar, ix)
              grad_y_buf(iy,ix) = ci1 * krho(ymode) * &
                 & buf(is,ymode, ipar, ix)
            end if

            if(xmode /= 1) then
              ix_sign = -1
              ix = ixzero + ix_sign*(xmode-1)
              do jv = 1, nmu
                do kt = 1, nvpar
                  d3v = bn(ix,ipar) * intmu(jv) * intvp(ipar, jv, kt, is)
                  v_exb_x(iy,ix) = v_exb_x(iy,ix) + d3v *&
                     & (efun(ix,ipar,1,1)*ci1*kxrh(ix) + &
                     & efun(ix,ipar,2,1)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix,ipar,jv,is) * &
                     & phi(ymode,ix,ipar)
                  v_exb_y(iy,ix) = v_exb_y(iy,ix) + d3v *&
                     & (efun(ix,ipar,1,2)*ci1*kxrh(ix) + &
                     & efun(ix,ipar,2,2)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix,ipar,jv,is) * &
                     & phi(ymode,ix,ipar)
                end do
              end do
              if(proc_subset(0,gipar,1,1,ispg)) then
                grad_x_buf(iy,ix) = ci1 * kxrh(ix) * &
                   & buf(is,ymode, ipar, ix)
                grad_y_buf(iy,ix) = ci1 * krho(ymode) * &
                   & buf(is,ymode, ipar, ix)
              end if

            end if

            if(ymode /= 1) then
              !use hermitian symmetry:
              ! f(-kx, -ky) = f*(kx,ky)
              !and
              ! f(kx, -ky)  = f*(-kx,ky)
              iy_sign = -1
              iy = iyzero + iy_sign*(ymode-1)

              ix_sign = +1
              ix = ixzero + ix_sign*(xmode-1)
              ix_theothersign = ixzero - ix_sign*(xmode-1)
              do jv = 1, nmu
                do kt = 1, nvpar
                  d3v = bn(ix,ipar) * intmu(jv) * intvp(ipar, jv, kt, is)
                  v_exb_x(iy,ix) = v_exb_x(iy,ix) + d3v * &
                     & conjg((efun(ix_theothersign,ipar,1,1)*ci1*kxrh(ix_theothersign) + &
                     & efun(ix_theothersign,ipar,2,1)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * &
                     & phi(ymode,ix_theothersign,ipar))
                  v_exb_y(iy,ix) = v_exb_y(iy,ix) + d3v *&
                     & conjg((efun(ix_theothersign,ipar,1,2)*ci1*kxrh(ix_theothersign) + &
                     & efun(ix_theothersign,ipar,2,2)*ci1*krho(ymode)) * &
                     & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * &
                     & phi(ymode,ix_theothersign,ipar))
                end do
              end do
              if(proc_subset(0,gipar,1,1,ispg)) then
                grad_x_buf(iy,ix) = conjg(ci1 * kxrh(ix_theothersign) * &
                   & buf(is,ymode, ipar, ix_theothersign))
                grad_y_buf(iy,ix) = conjg(ci1 * krho(ymode) * &
                   & buf(is,ymode, ipar, ix_theothersign))
              end if

              if(xmode /= 1) then
                ix_sign = -1
                ix = ixzero + ix_sign*(xmode-1)
                ix_theothersign = ixzero - ix_sign*(xmode-1)
                do jv = 1, nmu
                  do kt = 1, nvpar
                    d3v = bn(ix,ipar) * intmu(jv) * intvp(ipar, jv, kt, is)
                    v_exb_x(iy,ix) = v_exb_x(iy,ix) + d3v * &
                       & conjg((efun(ix_theothersign,ipar,1,1)*ci1*kxrh(ix_theothersign) + &
                       & efun(ix_theothersign,ipar,2,1)*ci1*krho(ymode)) * &
                       & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * &
                       & phi(ymode,ix_theothersign,ipar))
                    v_exb_y(iy,ix) = v_exb_y(iy,ix) + d3v * &
                       & conjg((efun(ix_theothersign,ipar,1,2)*ci1*kxrh(ix_theothersign) + &
                       & efun(ix_theothersign,ipar,2,2)*ci1*krho(ymode)) * &
                       & besselj0_gkw(ymode,ix_theothersign,ipar,jv,is) * &
                       & phi(ymode,ix_theothersign,ipar))
                  end do
                end do
                if(proc_subset(0,gipar,1,1,ispg)) then
                  grad_x_buf(iy,ix) = conjg(ci1 * kxrh(ix_theothersign) * &
                     & buf(is,ymode, ipar, ix_theothersign))
                  grad_y_buf(iy,ix) = conjg(ci1 * krho(ymode) * &
                     & buf(is,ymode, ipar, ix_theothersign))
                end if

              end if
            end if
          end do
        end do
        ! end debug
        call mpireduce_sum_inplace(v_exb_x, shape(v_exb_x), COMM_VPAR_NE_MU_NE)
        call mpireduce_sum_inplace(v_exb_y, shape(v_exb_y), COMM_VPAR_NE_MU_NE)


        if(.not. proc_subset(0,gipar,1,1,ispg)) cycle
        
        ! Measure the interaction of different modes

        ! iterate over all mode pairs: mode2 is the "from-mode", mode1
        ! is the "mediating mode" and the other one is the "to-mode"
        loop_ymode2: do nonlin_transfer_ymode2 = 1, nmod
          do ymode_sign2 = +1, -1, -2
            if(nonlin_transfer_ymode2 == 1 .and. ymode_sign2 == -1) then
              cycle
            end if
            loop_ymode1: do nonlin_transfer_ymode1 = 1, nmod
              do ymode_sign1 = +1, -1, -2
                if(nonlin_transfer_ymode1 == 1 .and. ymode_sign1 == -1) then
                  cycle
                end if
                iy_mode1 = iyzero + ymode_sign1*(nonlin_transfer_ymode1-1)
                iy_mode2 = iyzero + ymode_sign2*(nonlin_transfer_ymode2-1)
                ymode_flat = nonlin_transfer_ymode1 + (nonlin_transfer_ymode2-1)*nmod

                iy = iyzero + ymode_sign1*(nonlin_transfer_ymode1-1) &
                   & + ymode_sign2*(nonlin_transfer_ymode2-1)
                if(iy < 1 .or. iy > n_y_grid) cycle

                if(ymode_sign1 == +1 .and. ymode_sign2 == +1) then
                  iy_sparse = NONLIN_TRANSFER_PLUSPLUS
                elseif(ymode_sign1 == +1 .and. ymode_sign2 == -1) then
                  iy_sparse = NONLIN_TRANSFER_PLUSMINUS
                elseif(ymode_sign1 == -1 .and. ymode_sign2 == +1) then
                  iy_sparse = NONLIN_TRANSFER_MINUSPLUS
                elseif(ymode_sign1 == -1 .and. ymode_sign2 == -1) then
                  iy_sparse = NONLIN_TRANSFER_MINUSMINUS
                end if

                loop_xmode2: do nonlin_transfer_xmode2 = 1, nxmod
                  do xmode_sign2 = +1, -1, -2
                    if(nonlin_transfer_xmode2 == 1 .and. xmode_sign2 == -1) then
                      cycle
                    end if
                    loop_xmode1: do nonlin_transfer_xmode1 = 1,nxmod
                      do xmode_sign1 = +1, -1, -2
                        if(nonlin_transfer_xmode1 == 1 .and. xmode_sign1 == -1) then
                          cycle
                        end if
                        ix_mode1 = ixzero + xmode_sign1*(nonlin_transfer_xmode1-1)
                        ix_mode2 = ixzero + xmode_sign2*(nonlin_transfer_xmode2-1)
                        xmode_flat = nonlin_transfer_xmode1 + (nonlin_transfer_xmode2-1)*nxmod

                        ix = ixzero + xmode_sign1*(nonlin_transfer_xmode1-1) &
                           & + xmode_sign2*(nonlin_transfer_xmode2-1)
                        ix_theothersign = ixzero - (ix-ixzero)
                        if(ix < 1 .or. ix > n_x_grid) cycle

                        if(xmode_sign1 == +1 .and. xmode_sign2 == +1) then
                          ix_sparse = NONLIN_TRANSFER_PLUSPLUS
                        elseif(xmode_sign1 == +1 .and. xmode_sign2 == -1) then
                          ix_sparse = NONLIN_TRANSFER_PLUSMINUS
                        elseif(xmode_sign1 == -1 .and. xmode_sign2 == +1) then
                          ix_sparse = NONLIN_TRANSFER_MINUSPLUS
                        elseif(xmode_sign1 == -1 .and. xmode_sign2 == -1) then
                          ix_sparse = NONLIN_TRANSFER_MINUSMINUS
                        end if

                        imod = abs(iy - iyzero) + 1

                        if(iy >= iyzero) then
                          cdum = buf(is,imod, ipar, ix)
                        else
                          cdum = conjg(buf(is,imod, ipar, ix_theothersign))
                        end if

                        ! only a s-integral here
                        dum = - ints(ipar) &
                           & * real( &
                           ! CHECK factor is conjugated or not?
                           & conjg(cdum) &
                           & * ( &
                           & v_exb_x(iy_mode1, ix_mode1)*&
                           & grad_x_buf(iy_mode2, ix_mode2) &
                           & + &
                           & v_exb_y(iy_mode1, ix_mode1)*&
                           & grad_y_buf(iy_mode2, ix_mode2) &
                           &) &
                           & )

                        if(dum > 0) then
                          nonlin_transfer(ymode_flat, &
                             & xmode_flat, ix_sparse, iy_sparse) = &
                             & nonlin_transfer(ymode_flat, &
                             & xmode_flat, ix_sparse, iy_sparse) &
                             & + dum
                        else
                          nonlin_transfer_neg(ymode_flat, &
                             & xmode_flat, ix_sparse, iy_sparse) = &
                             & nonlin_transfer_neg(ymode_flat, &
                             & xmode_flat, ix_sparse, iy_sparse) &
                             & + dum
                        end if

                      end do
                    end do loop_xmode1
                  end do
                end do loop_xmode2
              end do
            end do loop_ymode1
          end do
        end do loop_ymode2
      end do loop_ns
    end if procs_on_gipar_ispg
  end subroutine calc_nl_transfer_moments

  !-----------------------------------------------------------------------------
  !> 
  !-----------------------------------------------------------------------------
  subroutine calc_nl_transfer_entr_nonspec(iquantity, ispg, gipar)
    use dist,           only : fdisi, fdis_tmp, fmaxwl
    use dist, only : ghost_points_xf, stencil_side
    use matdat,         only : get_f_from_g
    use components,     only : signz, tmp
    use mode,           only : krho
    use geom,           only : efun, bn, ints
    use grid,           only : nmod, ns, nmu, nvpar, n_y_grid, nx
    use grid,           only : gs, proc_subset, lsp, gx
    use velocitygrid,   only : intmu, intvp
    use constants,      only : ci1
    use global, only : PHI_GA_FIELD, id_x
    use dist, only : ifdis
    use index_function, only : indx

    use general, only : gkw_abort
    use mpiinterface, only : mpireduce_sum_inplace
    use diagnos_cross_phase, only : get_field_or_moment
    use diagnos_generic, only : dfielddx
    use fields, only : get_averaged_phi
    use non_linear_terms, only : mbd, ibd

    use grid,       only : lsendrecv_x
    integer, intent(in) :: iquantity
    !> global species index. This is used to access moment_buf and
    !> gyroavg fields.
    integer, intent(in) :: ispg
    !> global parallel index. If this is 0, then a flux surf avg will
    !> be computed.
    integer, intent(in) :: gipar
    complex :: cdum
    integer :: imod, ix, iy, jv, kt, is, ipar
    real    :: dum
    integer :: iy_sparse, ymode_flat
    !> loop boundaries, will cause loops or no loops
    integer :: nonlin_transfer_ymode1, nonlin_transfer_ymode2
    integer :: ymode, iyzero
    integer :: iy_sign, iy_mode1, iy_mode2
    integer :: ymode_sign1, ymode_sign2
    complex, dimension(n_y_grid) :: v_exb_x, v_exb_y, grad_x_buf, grad_y_buf
    complex, dimension(nx,nmod) :: dphigadx
    complex, dimension(nmod,1-stencil_side(id_x):nx+stencil_side(id_x)) :: buf
    integer :: iy_theothersign

    ! this is NOT the iyzero from module mode! This is the index of
    ! the ky=0 mode in a grid of n_y_grid modes, where the zeromode is
    ! in the centre (analogous to the kx grid we are working with).
    iyzero = nmod

    if(proc_subset(0,gipar,0,0,ispg)) then
      
      is = lsp(ispg)
      loop_ns: do ipar = 1, ns
        ! setting nonlin_transfer_ipar=0 gives a flux-surface-avg of the transfer,
        ! while >= 1 gives transfer in the chosen slice.
        if(gipar /= 0 .and. gs(ipar) /= gipar) cycle


        loop_mu: do jv = 1, nmu
          loop_vpar: do kt = 1, nvpar
            do imod = 1, nmod
              call dfielddx(PHI_GA_FIELD, imod, ipar, jv, is, dphigadx(:,imod))
            end do
            
            ! copy out the distribution function 
            do ix = 1-ghost_points_xf, nx+ghost_points_xf
              do imod = 1, nmod
                if (lsendrecv_x) then
                  buf(imod,ix) = fdis_tmp(indx(ifdis,imod,ix,ipar,jv,kt,is))
                  !buf(imod,ix,th) = fdis_tmp(lindx3(idx,i,j,is))
                else
                  buf(imod,ix) = fdisi(indx(ifdis,imod,ix,ipar,jv,kt,is))
                  !buf(imod,ix,th) = fdisi(lindx3(idx,i,j,is))
                end if
                !idx = idx + 1
              end do
            end do
            ! the boundary conditions 
            do imod = 1, nmod
              do ix = 1, stencil_side(id_x)
                buf(imod,1-ix) = mbd(1-ix,imod,ipar)*buf(imod,ibd(1-ix,imod,ipar))
                buf(imod,nx+ix) = mbd(ix,imod,ipar)*buf(imod,ibd(ix,imod,ipar))
              end do
            end do

            loop_nx: do ix = 1, nx

              v_exb_x = 0.0
              v_exb_y = 0.0
              grad_x_buf = 0.0
              grad_y_buf = 0.0

              ! fill buffers with the "complete" k-space fields. All
              ! values are stored explicitely, despite of hermitian
              ! symmetry, to make life easier later.
              do ymode = 1, nmod
                iy_sign = +1
                iy = iyzero + iy_sign*(ymode-1)

                ! efun is actually
                ! symmetric and the diagnonal elements are
                ! zero, but for readability I write it out
                v_exb_x(iy) = v_exb_x(iy) + &
                   & efun(ix,ipar,1,1)*dphigadx(ix,ymode) + &
                   & efun(ix,ipar,2,1)*ci1*krho(ymode) * &
                   & get_averaged_phi(ymode,ix,ipar, jv,is,fdisi)
                v_exb_y(iy) = v_exb_y(iy) + &
                   & efun(ix,ipar,1,2)*dphigadx(ix,ymode) + &
                   & efun(ix,ipar,2,2)*ci1*krho(ymode) * &
                   & get_averaged_phi(ymode,ix,ipar, jv,is,fdisi)
                ! FIXME for now, just hardcode the radial finite difference
                ! of the distribution function
                grad_x_buf(iy) = buf(ymode,ix-2)-8.0E0*buf(ymode,ix-1) &
                   & + 8.0E0*buf(ymode,ix+1)-buf(ymode,ix+2)
                grad_y_buf(iy) = grad_y_buf(iy) + ci1 * krho(ymode) * &
                   & fdisi(indx(ifdis,ymode, &
                   & ix,ipar,jv,kt,is))

                if(ymode /= 1) then
                  !use hermitian symmetry:
                  ! f(-ky) = f*(ky)
                  iy_sign = -1
                  iy = iyzero + iy_sign*(ymode-1)
                  iy_theothersign = iyzero - iy_sign*(ymode-1)

                  v_exb_x(iy) = conjg(v_exb_x(iy_theothersign))
                  v_exb_y(iy) = conjg(v_exb_y(iy_theothersign))
                  grad_x_buf(iy) = conjg(grad_x_buf(iy_theothersign))
                  grad_y_buf(iy) = conjg(grad_y_buf(iy_theothersign))
                end if
              end do

              ! Measure the interaction of different modes

              ! iterate over all mode pairs: mode2 is the "from-mode", mode1
              ! is the "mediating mode" and the other one is the "to-mode"
              loop_ymode2: do nonlin_transfer_ymode2 = 1, nmod
                do ymode_sign2 = +1, -1, -2
                  if(nonlin_transfer_ymode2 == 1 .and. ymode_sign2 == -1) then
                    cycle
                  end if
                  loop_ymode1: do nonlin_transfer_ymode1 = 1, nmod
                    do ymode_sign1 = +1, -1, -2
                      if(nonlin_transfer_ymode1 == 1 .and. ymode_sign1 == -1) then
                        cycle
                      end if
                      iy_mode1 = iyzero + ymode_sign1*(nonlin_transfer_ymode1-1)
                      iy_mode2 = iyzero + ymode_sign2*(nonlin_transfer_ymode2-1)
                      ymode_flat = nonlin_transfer_ymode1 + (nonlin_transfer_ymode2-1)*nmod

                      iy = iyzero + ymode_sign1*(nonlin_transfer_ymode1-1) &
                         & + ymode_sign2*(nonlin_transfer_ymode2-1)
                      if(iy < 1 .or. iy > n_y_grid) cycle

                      if(ymode_sign1 == +1 .and. ymode_sign2 == +1) then
                        iy_sparse = NONLIN_TRANSFER_PLUSPLUS
                      elseif(ymode_sign1 == +1 .and. ymode_sign2 == -1) then
                        iy_sparse = NONLIN_TRANSFER_PLUSMINUS
                      elseif(ymode_sign1 == -1 .and. ymode_sign2 == +1) then
                        iy_sparse = NONLIN_TRANSFER_MINUSPLUS
                      elseif(ymode_sign1 == -1 .and. ymode_sign2 == -1) then
                        iy_sparse = NONLIN_TRANSFER_MINUSMINUS
                      end if

                      imod = abs(iy - iyzero) + 1
                      cdum = (get_f_from_g(imod,ix,ipar,jv,kt,is,fdisi) &
                         & / fmaxwl(ix, ipar, jv, kt, is) + signz(is)/tmp(ix,is) &
                         & * get_averaged_phi(imod,ix,ipar,jv,is,fdisi))
                      if(iy < iyzero) then
                        cdum = conjg(cdum)
                      end if

                      dum = - ints(ipar) * bn(ix,ipar) * intmu(jv) * intvp(ipar, jv, kt, is) &
                         & * real( &
                         & conjg(cdum) &
                         & * ( &
                         & v_exb_x(iy_mode1)*&
                         & grad_x_buf(iy_mode2) &
                         & + &
                         & v_exb_y(iy_mode1)*&
                         & grad_y_buf(iy_mode2) &
                         &) &
                         & )

                      if(dum > 0) then
                        nonlin_transfer(ymode_flat, &
                           & gx(ix), 1, iy_sparse) = &
                           & nonlin_transfer(ymode_flat, &
                           & gx(ix), 1, iy_sparse) &
                           & + dum
                      else
                        nonlin_transfer_neg(ymode_flat, &
                           & gx(ix), 1, iy_sparse) = &
                           & nonlin_transfer_neg(ymode_flat, &
                           & gx(ix), 1, iy_sparse) &
                           & + dum
                      end if
                    end do
                  end do loop_ymode1
                end do
              end do loop_ymode2
            end do loop_nx
          end do loop_vpar
        end do loop_mu
      end do loop_ns
    end if
  end subroutine calc_nl_transfer_entr_nonspec
  
  !------------------------------------------------------------------------------
  !>
  !------------------------------------------------------------------------------
  pure function nonlin_transfer_get_imodsparse(imod,ymode1,ymode2) &
     & result(iysparse)
    integer, intent(in) :: imod, ymode1, ymode2
    integer :: iysparse

    ! because it is simpler here, this function calculates examines a
    ! zero-based index (ymode1=ymode2=1 yields m = 0, not 1)

    ! second summand is surely >= 0
    if(imod-1 == (ymode1-1) + (ymode2-1)) then
      ! this is (k1+k2):
      iysparse = NONLIN_TRANSFER_PLUSPLUS
    else if(imod-1 == (ymode1-1) -(ymode2-1)) then
      ! this is (k1-k2):
      iysparse = NONLIN_TRANSFER_PLUSMINUS
      ! else if(imod == iyzero + abs(-(ymode1-1) +(ymode2-1))) then
      !   ! this is (-k1+k2)
      !   iysparse = 3
    else
      iysparse = 0
    end if

  end function nonlin_transfer_get_imodsparse

  !------------------------------------------------------------------------------
  !>
  !------------------------------------------------------------------------------
  pure function nonlin_transfer_get_ixsparse(ix,xmode1,xmode2) result(ixsparse)
    use mode, only : ixzero
    !  use grid, only : n_x_grid
    integer, intent(in) :: ix, xmode1, xmode2
    integer :: ixsparse

    ! second summand is surely >= 0
    if(ix == ixzero + ((xmode1-1) + (xmode2-1))) then
      ! this is (k1+k2):
      ixsparse = NONLIN_TRANSFER_PLUSPLUS

      ! second summand is surely < 0 (it would be ==0 for xmode12==1, but this has
      ! already been catched by the first test)
      !else if(ix == (n_x_grid + 1) + (-(xmode1-1) -(xmode2-1))) then
    else if(ix == ixzero + (-(xmode1-1) -(xmode2-1))) then
      ! this is (-k1-k2):
      ixsparse = NONLIN_TRANSFER_MINUSMINUS

      ! second summand is surely >= 0
    else if(ix == ixzero + ((xmode1-1) -(xmode2-1))) then
      ! this is (k1-k2):
      ixsparse = NONLIN_TRANSFER_PLUSMINUS

      ! second summand is surely < 0 (it would be ==0 for xmode1==xmode2, but this has
      ! already been catched by the preceeding test)
    else if(ix == ixzero + (-(xmode1-1) + (xmode2-1))) then
      ! this is (-k1+k2)
      ixsparse = NONLIN_TRANSFER_MINUSPLUS

    else
      ixsparse = 0
    end if

  end function nonlin_transfer_get_ixsparse


end module diagnos_nonlin_transfer
