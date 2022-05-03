!------------------------------------------------------------------------------
!> Add module description here
!------------------------------------------------------------------------------
module diagnos_grid

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output
  public :: read_last_data
  public :: calc_smallstep
  public :: output

  integer, save :: lun_time, lun_timefine, lun_filecount

  real, save, allocatable :: timefine_buf(:)
contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor
    use control, only : io_legacy
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    logical, intent(inout) :: requirements(:,:)

    ! To keep the compiler quiet, as the array is not used so far (in this routine).
    if (.false.) write(*,*) requirements

    if(root_processor) then
      if(.not. io_legacy) then
        call open_real_lu('time', 'grid', &
           & (/ 1 /), ascii_fmt, lun_time)
        call attach_metadata(lun_time, phys_unit_key, &
           & 'R_{ref}/v_{th,ref}', ascii_fmt)
        call attach_metadata(lun_time, description_key, &
           & 'The time grid, including only large timesteps.', &
           & ascii_fmt)
        call attach_metadata(lun_time, comments_key, &
           & 'GKW performs naverage "small timesteps" of length dtim per &
           & "large timestep". Most timedependent diagnostic output data &
           & refers to this time grid.', ascii_fmt)
      end if
      
      call open_real_lu('time_fine', 'grid', &
         & (/ 1 /), ascii_fmt, lun_timefine)
      call attach_metadata(lun_timefine, phys_unit_key, &
         & 'R_{ref}/v_{th,ref}', ascii_fmt)
      call attach_metadata(lun_timefine, description_key, &
         & 'The time grid, including all small timesteps.', &
         & ascii_fmt)
      call attach_metadata(lun_timefine, comments_key, &
         & 'GKW performs naverage "small timesteps" of length dtim per &
         & "large timestep". (For every higher order explicit integration &
         & a scheme with substeps is used, but those do not appear anywhere in &
         & the output data.)', ascii_fmt)

      call open_real_lu('file_count', 'grid', &
         & (/ 1 /), ascii_fmt, lun_filecount)
      call attach_metadata(lun_timefine, description_key, &
         & 'File number used in per time step diagnostic outputs.', &
         & ascii_fmt)

    end if
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use control, only : naverage
    use mpiinterface, only : root_processor
    use general, only : gkw_abort
    integer :: ierr
    
    if(.not. root_processor) return
    allocate(timefine_buf(naverage),stat=ierr)
    if (ierr /= 0) &
         & call gkw_abort('diagnos_grid :: could not allocate &
         & timefine_buf')

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : ascii_fmt, close_lu
    use mpiinterface, only : root_processor

    if(root_processor) then
      call close_lu(lun_time, ascii_fmt)
      call close_lu(lun_timefine, ascii_fmt)
      call close_lu(lun_filecount, ascii_fmt)
      deallocate(timefine_buf)
    end if
  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()
    use mode, only : mode_label_G
    use mpiinterface, only : root_processor
    use io, only : output_array, ascii_fmt, attach_metadata
    use io, only : phys_unit_key, description_key, comments_key, not_avail

    ! write the model_label array to a file. It gives information about
    ! which modes are coupled over the boundary conditions.
    if(root_processor) then
      call output_array('mode_label', 'diagnostic/diagnos_grid', &
         & mode_label_G*1.0, 'F', &
         & '(4096(F5.0,1x))', ascii_fmt)
      call attach_metadata('mode_label', 'diagnostic/diagnos_grid', &
         & phys_unit_key, not_avail, &
         & ascii_fmt)
      call attach_metadata('mode_label', 'diagnostic/diagnos_grid', &
         & description_key, not_avail, &
         & ascii_fmt)
      call attach_metadata('mode_label', 'diagnostic/diagnos_grid', &
         & comments_key, not_avail, &
         & ascii_fmt)
    end if

    call write_xykxky_grids
    call write_kx_connections
    
    call write_box_parameters
    
    call write_s_grid

    call write_velspace_grids
    call write_nonuni_velspace_grids
    
  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine write_kx_connections
    use mpiinterface, only : root_processor
    use grid, only : nmod, nx
    use mode, only : ixplus, ixminus, kxrh, ikxspace, mode_box
    use geom, only : q, eps, shat
    use io, only : get_free_file_unit, clean0r, output_enabled
    integer :: lun
    integer :: imod, ix
    ! Write the kx connections to file
    if (.not.root_processor) return
    
    if (mode_box .and. output_enabled) then
      call get_free_file_unit(lun)
      open(lun,FILE="kx_connect.dat", FORM='formatted', POSITION='asis')
22    format('shat: ',f6.2,2x,'q: ',f6.2,2x,'eps: ',f6.2,2x,'ikxspace: ', i4)
23    format(i6,f10.6,i6,i6,i6)
      write(lun,22) shat, q, eps, ikxspace
      write(lun,*) 'imod, kxrh(ix), ix, ixplus(imod,ix), ixminus(imod,ix)'
      do imod = 1, nmod
        do ix = 1, nx
          write(lun,23) imod, clean0r(kxrh(ix)), ix, &
             & ixplus(imod,ix), ixminus(imod,ix)
        end do
      end do
      close(lun)
    end if
  end subroutine write_kx_connections

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine read_last_data()
    use control,      only : nt_complete, io_legacy
    use control,      only : last_largestep_time
    use control,      only : time
    use io,           only : lu_exists, ascii_fmt, close_lu
    use io,           only : read_last_chunk, open_real_lu
    use mpiinterface, only : root_processor, mpibcast
    use general,      only : gkw_warn
    integer :: lun, nlast
    real, dimension(1) :: last_chunk

    last_largestep_time = 0.0

    if(.not. io_legacy) then
      if(root_processor) then
        if(.not. lu_exists('time', 'grid', ascii_fmt)) then
          write(*,*) 'file or dataset grid/time not found. Do not read last value.'
          ! Otherwise the time value set in the restart
          ! module is kept.
        else
          call open_real_lu('time', 'grid', &
           & (/ 1 /), ascii_fmt, lun)
          call read_last_chunk(lun, '(3(es13.5))', last_chunk(1:1), &
             & nlast, ascii_fmt)
          call close_lu(lun, ascii_fmt)

          time = last_chunk(1)

          if(nt_complete /= nlast) then
            call gkw_warn('The logical unit "time" '// &
               & 'is shorter or longer than expected.');
            write (*,*) "nt_complete = ", nt_complete
            write (*,*) "number of chunks in logical unit = ", nlast
          end if
        end if
      end if
      
      call mpibcast(time,1)
    end if
    
  end subroutine read_last_data


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> 
  !---------------------------------------------------------------------------
  subroutine write_xykxky_grids
    use io, only : output_array, xy_fmt, ascii_fmt, attach_metadata
    use io, only : phys_unit_key, description_key, comments_key, not_avail
    use mode,             only : lyn, lxn, kxrh, krho, mode_box, krloc
    use grid,             only : ns, n_s_grid, nx, n_x_grid, nmod, proc_subset
    use grid,             only : n_y_grid
    use geom,             only : kthnorm, xgr
    use control,          only : spectral_radius, flux_tube, non_linear
    use control, only : io_legacy
    use non_linear_terms, only : mphiw3
    use mpiinterface,     only : root_processor, gather_array
    use mpicomms, only : COMM_S_NE_X_NE
    use diagnos_generic, only :  mphit, mrad_G, mpi_dtype_real_spec_yxs
    use general, only : gkw_abort

    integer :: i,j,ix,imod
    integer :: ierr
    real, allocatable :: krloc_G(:,:,:)

    if(proc_subset(0,0,1,1,1)) then
      if(root_processor) then
        allocate(krloc_G(nmod, n_x_grid, n_s_grid),stat=ierr)
      else
        allocate(krloc_G(1,1,1),stat=ierr)
      end if
      if (ierr /= 0) call gkw_abort('Could not allocate krloc_G in diagnos_grid')
      call gather_array(krloc_G, &
         & krloc(:,1:nx,1:ns), &
         & mpi_dtype_real_spec_yxs, COMM_S_NE_X_NE, .true., .true.)
      
      if(root_processor) then
        call output_array('krloc', 'grid', krloc_G, &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata('krloc', 'grid', phys_unit_key, not_avail, &
           & ascii_fmt)
        call attach_metadata('krloc', 'grid', description_key, &
           & 'local perpendicular wavenumber', ascii_fmt)
        call attach_metadata('krloc', 'grid', comments_key, &
           & 'kzeta', &
           & ascii_fmt)

      end if
      deallocate(krloc_G)
    end if

    if(.not. root_processor) return

    ! phi diagnostic related outputs
    modebox : if (mode_box) then
      
      if(flux_tube) then
        call output_array('xgr', 'grid', &
           & reshape( (/ (real(j-0.5)*lxn/real(n_x_grid)-lxn/2.0E0, &
           & j = 1, n_x_grid) /),(/ n_x_grid /)), &
           & 'F', xy_fmt, ascii_fmt)
        call attach_metadata('xgr', &
          & 'grid', phys_unit_key, &
          & '\rho_{ref}', &
          & ascii_fmt)
        call attach_metadata('xgr', &
          & 'grid', description_key, &
          & 'radial coordinate grid', ascii_fmt)
        call attach_metadata('xgr', &
          & 'grid', comments_key, &
          & 'In the local limit the radial extend of the box is  &
          & lxn (in units of \rho_{ref}). The radial grid is then &
          & set by lxn devided by the number of radial grid &
          & points (n_x_grid). Here, the zero point of the radial grid &
          & is chosen to be the center of the flux tube.', ascii_fmt)
      else
        call output_array('xgr', &
          & 'grid', xgr, &
          & 'F', xy_fmt, ascii_fmt)
        call attach_metadata('xgr', &
          & 'grid', phys_unit_key, &
          & 'dimensionless', &
          & ascii_fmt)
        call attach_metadata('xgr', &
          & 'grid', description_key, &
          & 'radial coordinate grid', ascii_fmt)
        call attach_metadata('xgr', &
          & 'grid', comments_key, &
          & 'In GKW the radial coordinate is epsilon, the inverse aspect ratio. &
          & This quantity appears &
          & also as eps in geom. Some people call the radial coord. the &
          & flux tube label.', ascii_fmt)
      endif

      ! the xphi grid goes from 0 to lxn. This is only useful for
      ! fluxtube runs. A global radial grid must go from psil to psih.
      if(io_legacy) then
        call output_array('xphi', 'grid', reshape( &
           & (/ ((real(j-1)*lxn/real(mrad_G), j = 1, mrad_G), i = 1, mphit) /), &
           & (/ mrad_G, mphit /)), &
           & 'F', xy_fmt, ascii_fmt)
      else
        call output_array('xphi', 'grid', reshape( &
           & (/ ((real(j-1)*lxn/real(n_x_grid), j = 1, n_x_grid), i = 1, n_y_grid) /), &
           & (/ n_x_grid, n_y_grid /)), &
           & 'F', xy_fmt, ascii_fmt)
      end if
      call attach_metadata('xphi', 'grid', phys_unit_key, '\rho_{ref}', &
         & ascii_fmt)
      call attach_metadata('xphi', 'grid', description_key, &
         & 'radial coordinate grid', ascii_fmt)
      call attach_metadata('xphi', 'grid', comments_key, &
         & 'The motivation for the xphi,yphi outputs was that when you do &
         & contourf(xphi,yphi,some_2d_data) &
         & You get an output in which the units of both x and y directions &
         & are the same, i.e. in units of rho_ref. ', &
         & ascii_fmt)

      if(io_legacy) then
        call output_array('yphi', 'grid', reshape( &
           & (/ ((real(i-1)*lyn/real(mphit), j = 1, mrad_G), i = 1, mphit) /), &
           & (/ mrad_G, mphit /)), &
           & 'F', xy_fmt, ascii_fmt)
      else
        call output_array('yphi', 'grid', reshape( &
           & (/ ((real(i-1)*lyn/real(n_y_grid), j = 1, n_x_grid), i = 1, n_y_grid) /), &
           & (/ n_x_grid, n_y_grid /)), &
           & 'F', xy_fmt, ascii_fmt)
      end if
      call attach_metadata('yphi', 'grid', phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata('yphi', 'grid', description_key, not_avail, ascii_fmt)
      call attach_metadata('yphi', 'grid', comments_key, &
         ' Actually the grid yphi is not the real space coordinate zeta but &
         & its projection onto theta, i.e. 2 pi/krho_theta .', &
         & ascii_fmt)
    end if modebox

    if (spectral_radius) then
      call output_array('kxrh', 'grid', reshape( &
         & (/ (kxrh, imod = 1, nmod) /), &
         & (/ n_x_grid, nmod /)), &
         & 'F', xy_fmt, ascii_fmt)
      call attach_metadata('kxrh', 'grid', phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata('kxrh', 'grid', description_key, not_avail, ascii_fmt)
      call attach_metadata('kxrh', 'grid', comments_key, not_avail, ascii_fmt)
    end if

    if(flux_tube) then
      call output_array('krho', 'grid', reshape (&
         & (/ ((krho(imod)*kthnorm, ix = 1,n_x_grid), imod = 1,nmod) /), &
         & (/ n_x_grid, nmod /)), &
         & 'F', xy_fmt, ascii_fmt)

      if(non_linear .and. nmod > 1) then
        ! in order to compute the nonlinear term, a finer spatial
        ! grid/a larger k-space grid is used. Output the wavevectors
        ! of that grid:
        ! FIXME this is a bit sloppy
        call output_array('krho_extended', 'grid', reshape (&
           & (/ (krho(2)*kthnorm* imod, imod = 0,mphiw3) /), &
           & (/ mphiw3 /)), &
           & 'F', xy_fmt, ascii_fmt)
      end if
    end if

    call output_array('kzeta', 'grid', reshape (&
       & (/ ((krho(imod), ix = 1,n_x_grid), imod = 1,nmod) /), &
       & (/ n_x_grid, nmod /)), &
       & 'F', xy_fmt, ascii_fmt)


  end subroutine write_xykxky_grids

  !-------------------------------------------------------------------------
  !> Output velocity grids at s point xy_slice_ipar
  !>
  !-------------------------------------------------------------------------
  subroutine write_velspace_grids
    use grid,           only : nmu,nvpar, ls, gx, lsp
    use grid,           only : proc_subset
    use velocitygrid,   only : intmu,intvp, vpgr, mugr
    use mode,           only : ixzero
    use diagnos_generic, only : velocity_slice_output, xy_slice_ipar
    use io, only : ascii_fmt
    use control, only : io_legacy
    use global, only : dotdat
    use mpiinterface, only : register_tag_range
    integer :: i, j, isl, ixg, ispg
    real, allocatable, dimension(:,:) :: local_vpar_mu
    integer :: tag_range_start, tag_range_end_inkl

    call register_tag_range(4, &
       & tag_range_start, tag_range_end_inkl)

    ! set ix here to ixzero 
    ixg = ixzero
    ispg = 1

    ! allocate arrays to contain the fine slice and local slice
    allocate(local_vpar_mu(nvpar,nmu))

    !USE xy_slice_ipar as for xy-slices
    !local s
    isl=ls(xy_slice_ipar)

    if(proc_subset(ixg,xy_slice_ipar,0,0,ispg)) then
      ! copy the local array into the buffer
      do j=1,nmu
        do i=1,nvpar
          local_vpar_mu(i,j) = intvp(isl,j,i,lsp(ispg))
        end do
      end do
    end if

    call velocity_slice_output('diagnostic/diagnos_grid', &
       & local_vpar_mu,dotdat('intvp',io_legacy),ixg,xy_slice_ipar,ispg, &
       & ascii_fmt, tag_range_start)

    if(proc_subset(ixg,xy_slice_ipar,0,0,ispg)) then
      ! copy the local array into the buffer
      do j=1,nmu
        ! bn should not feature here: the velocity grid is normalised to Bref
        ! The velocity of a particle is directly the value of the grid
        !  local_vpar_mu(i,j)=sqrt(2.*bn(ix,isl)*mugr(j))
        local_vpar_mu(:,j)=sqrt(2.*mugr(j))
      end do
    end if
    ! former filename: distr2.dat
    call velocity_slice_output('grid', &
       & local_vpar_mu,dotdat('vperp', io_legacy),ixg,xy_slice_ipar,ispg, &
       & ascii_fmt, tag_range_start+1)


    if(proc_subset(ixg,xy_slice_ipar,0,0,ispg)) then
      ! copy the local array into the buffer
      do j=1,nmu
        local_vpar_mu(:,j)=intmu(j)
      end do
    end if

    call velocity_slice_output('diagnostic/diagnos_grid', &
       & local_vpar_mu,dotdat('intmu',io_legacy),ixg,xy_slice_ipar,ispg, &
       & ascii_fmt, tag_range_start+2)

    if(proc_subset(ixg,xy_slice_ipar,0,0,ispg)) then
      ! copy the local array into the buffer
      do j=1,nmu
        do i=1,nvpar
          local_vpar_mu(i,j) = vpgr(isl,j,i,lsp(ispg))
        end do
      end do
    end if

    !Reorder array to do this
    !local_vpar_mu(:,:)=vpgr(isl,:,:)

    !former filename: distr1.dat
    call velocity_slice_output('grid', &
       & local_vpar_mu, dotdat('vpgr',io_legacy),ixg,xy_slice_ipar,ispg, &
       & ascii_fmt, tag_range_start+3)

  end subroutine write_velspace_grids

  !-------------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------------
  subroutine write_nonuni_velspace_grids
    use mpiinterface, only : root_processor
    use grid, only : nmu, nvpar, ns
    use velocitygrid, only : mugr, vpgr, intvp
    use control, only : vp_trap
    integer :: j, k

    if(.not. root_processor) return
    
    ! check whether the parallel velocity grid is set up with
    ! the parallel velocity following the trapping condition
    if (vp_trap.ne.1) return
    
    ! write the grid
    open(18,file = 'velgrid') !The grid at the high field point
    do j = 1, nmu
      do k =  1, nvpar
        write(18,*) sqrt(mugr(j)),vpgr((ns+1)/2,j,k,1)
      end do
    end do
    
    open(18,file = 'exintvp') !The grid at the high field point
    do j = 1, nmu
      do k =  1, nvpar
        write(18,*) intvp((ns+1)/2,j,k,1)
      end do
    end do

  end subroutine write_nonuni_velspace_grids


  !-------------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------------
  subroutine write_box_parameters
    use mode,             only : lyn, lxn, mode_box
    use non_linear_terms, only : mphi, mphiw3
    use mpiinterface,     only : root_processor
    use diagnos_generic,  only : mrad_G, mrad_l
    use io, only : output_array, xy_fmt, ascii_fmt, attach_metadata
    use io, only : phys_unit_key, description_key, comments_key, not_avail

    if(.not. root_processor) return

    if(mode_box) then
    call output_array('lxn','diagnostic/diagnos_grid', (/ lxn /), 'F', &
       & xy_fmt, ascii_fmt)
    call attach_metadata('lxn', 'diagnostic/diagnos_grid', &
       & phys_unit_key, not_avail, &
       & ascii_fmt)
    call attach_metadata('lxn', 'diagnostic/diagnos_grid', &
       & description_key, 'radial length of the computational box', &
       & ascii_fmt)
    call attach_metadata('lxn', 'diagnostic/diagnos_grid', &
       & comments_key, not_avail, &
       & ascii_fmt)

    call output_array('lyn','diagnostic/diagnos_grid', (/ lyn /), 'F', &
       & xy_fmt, ascii_fmt)
    call attach_metadata('lyn', 'diagnostic/diagnos_grid', &
       & phys_unit_key, not_avail, &
       & ascii_fmt)
    call attach_metadata('lyn', 'diagnostic/diagnos_grid', &
       & description_key, 'binormal length of the computational box', &
       & ascii_fmt)
    call attach_metadata('lyn', 'diagnostic/diagnos_grid', &
       & comments_key, not_avail, &
       & ascii_fmt)
    end if

    call output_array('mrad_G','diagnostic/diagnos_grid', (/ mrad_G*1.0 /), &
       & 'F', xy_fmt, ascii_fmt)
    call attach_metadata('mrad_G', 'diagnostic/diagnos_grid', &
       & description_key, 'global radial size of the higher resolved grid which&
       & is used in the computation of nonlinear terms', &
       & ascii_fmt)

    call output_array('mrad_l','diagnostic/diagnos_grid', (/ mrad_l*1.0 /), &
       & 'F', xy_fmt, ascii_fmt)
    call attach_metadata('mrad_l', 'diagnostic/diagnos_grid', &
       & description_key, 'local radial size of the higher resolved grid which&
       & is used in the computation of nonlinear terms', &
       & ascii_fmt)

    call output_array('mphi','diagnostic/diagnos_grid', (/ mphi*1.0 /), &
       & 'F', xy_fmt, ascii_fmt)
    call attach_metadata('mphi', 'diagnostic/diagnos_grid', &
       & description_key, 'binormal size of the higher resolved grid, which &
       & is used in the computation of the nonlinear terms', &
       & ascii_fmt)

    call output_array('mphiw3','diagnostic/diagnos_grid', (/ mphiw3*1.0 /), &
       & 'F', xy_fmt, ascii_fmt)
    call attach_metadata('mphiw3', 'diagnostic/diagnos_grid', &
       & description_key, 'binormal size of the higher resolved grid, which &
       & is used in the computation of the nonlinear terms, in k-space.', &
       & ascii_fmt)

  end subroutine write_box_parameters


  !-------------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------------
  subroutine write_s_grid
    use grid,             only : n_s_grid
    use geom,             only : sgr
    use io, only : output_array, xy_fmt, ascii_fmt, attach_metadata
    use io, only : phys_unit_key, description_key, comments_key
    use mpiinterface, only : root_processor
    
    if(.not.root_processor) return

    call output_array('sgrid','diagnostic/diagnos_grid', sgr(1:n_s_grid), 'F', &
       & xy_fmt, ascii_fmt)
    call attach_metadata('sgrid', 'diagnostic/diagnos_grid', &
       & phys_unit_key, 'dimensionless', &
       & ascii_fmt)
    call attach_metadata('sgrid', 'diagnostic/diagnos_grid', &
       & description_key, 'parallel coordinate s', &
       & ascii_fmt)
    call attach_metadata('sgrid', 'diagnostic/diagnos_grid', &
       & comments_key, 'The s coordinate is zero at the LFS midplane and increasing upward from this position', &
       & ascii_fmt)
    
 
  end subroutine write_s_grid

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_smallstep(i_smallstep)
    use control, only : time
    use mpiinterface, only : root_processor
    integer, intent(in) :: i_smallstep
    if(.not. root_processor) return

    ! at every small timestep:
    timefine_buf(i_smallstep) = time

  end subroutine calc_smallstep

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output(file_count)
    use control, only : time, naverage, io_legacy
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use mpiinterface, only : root_processor
    integer :: i
    integer, intent(in) :: file_count
    if(.not. root_processor) return

    if(.not. io_legacy) then
      call append_chunk(lun_time, (/ time /), '(1(es13.5))', ascii_fmt)
    end if

    do i = 1, naverage
      call append_chunk(lun_timefine, (/ timefine_buf(i) /), xy_fmt, ascii_fmt)
    end do

    call append_chunk(lun_filecount,(/ real(file_count) /), '(f10.0)', ascii_fmt)

  end subroutine output


end module diagnos_grid
