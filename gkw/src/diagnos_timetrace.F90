!------------------------------------------------------------------------------
!> 
!------------------------------------------------------------------------------
module diagnos_timetrace

  implicit none

  private

  !>==== DESCRIPTION ==================================================
  !>
  !> This diagnostic outputs a time trace, which can be subject to a
  !> signal analysis technique at postprocessing such as
  !> time-frequency analysis, higher-order spectral analysis,
  !> wavelet analysis or whatever makes sense.
  !>
  !>---- FURTHER NOTES ------------------------------------------------
  !>
  !> see also figure 5a in
  !> http://iopscience.iop.org/0029-5515/50/5/054005/
  !> or
  !> D. Told et al., PoP 20, 122312 (2013)
  !>   http://dx.doi.org/10.1063/1.4858899
  !>
  !>---- LIMITATIONS --------------------------------------------------
  !>
  !> This diagnostic does not care about at-runtime refinement of the
  !> timestep by any of the timestep estimators, and the analysis
  !> scripts probably have a hard time handling such cases (see
  !> comments there).
  !>
  !>===================================================================

  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: finalize
  public :: initial_output, read_last_data
  public :: output

  public :: calc_smallstep

  !> The general on/off switch for this diagnostic
  logical, save, public :: lcomplex_timetrace

  !> array for quick look up of the fields
  integer, save, dimension(:), allocatable :: lindx

  !> buffers
  complex, save, dimension(:,:), allocatable :: cbuf

  integer :: lun_csignal

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()

    lcomplex_timetrace = .false.

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcomplex_timetrace,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use general, only : gkw_warn

    if(.not.lcomplex_timetrace) return

  end subroutine check

  !---------------------------------------------------------------------------
  !> construct the array for quick look up of the fields
  !---------------------------------------------------------------------------
  function get_field_index_lookup_array() result(lindx)
    use grid, only : nx, nmod, ns
    use dist, only : iphi, iapar, ibpar, number_of_fields
    use index_function, only : indx
    use control, only : nlphi, nlapar, nlbpar
    integer, dimension(number_of_fields*ns*nx*nmod) :: lindx
    integer, save, allocatable, dimension(:) :: field_id
    integer :: nfields, k, ix, imod, is, idx

    ! construct the array for quick look up of the fields
    allocate(field_id(number_of_fields))
    nfields = 0
    if (nlphi) then
      nfields = nfields + 1
      field_id(nfields) = iphi
    end if
    if (nlapar) then
      nfields = nfields + 1
      field_id(nfields) = iapar
    end if
    if (nlbpar) then
      nfields = nfields + 1
      field_id(nfields) = ibpar
    end if

    idx = 0
    lindx(:) = 0
    do k = 1, number_of_fields
      do ix = 1, nx
        do imod = 1, nmod
          do is = 1, ns
            idx=idx+1
            lindx(idx) = indx(field_id(k),imod,ix,is)
          end do
        end do
      end do
    end do
    deallocate(field_id)

  end function get_field_index_lookup_array
  
  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor
    use io, only : open_complex_lu, ascii_fmt, binary_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use io, only : output_array, open_real_lu
    use io, only : xy_fmt
    use diagnos_generic, only : attach_metadata_grid
    use grid, only : nmod

    use constants, only : pi
    use general, only : gkw_abort
    logical, intent(inout) :: requirements(:,:)
    
    if(.not.lcomplex_timetrace) return

    if(root_processor) then
      call open_complex_lu('complex_signal', 'diagnostic/diagnos_timetrace', &
         & (/ nmod /), &
         & binary_fmt, lun_csignal)
      call attach_metadata_grid(lun_csignal, &
         & 'time_fine', 'krho', binary_fmt)
      call attach_metadata(lun_csignal, phys_unit_key, &
         & not_avail, binary_fmt)
      call attach_metadata(lun_csignal, description_key, &
         & 'Time trace for every binormal mode. &
         & ', &
         & binary_fmt)
      call attach_metadata(lun_csignal, comments_key, &
         & not_avail, binary_fmt)
    end if

    lindx = get_field_index_lookup_array()
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use general, only : gkw_abort
    use control, only : naverage
    use grid, only : ns, nx, nmod
    use mpiinterface, only : root_processor
    use dist, only : number_of_fields
    integer :: ierr

    if(.not.lcomplex_timetrace) return

    allocate(lindx(number_of_fields*ns*nx*nmod),stat=ierr)
    if (ierr /= 0) &
         & call gkw_abort('diagnos_timetrace :: could not allocate &
         & lindx')

    if(root_processor) then

      allocate(cbuf(nmod, naverage),stat=ierr)
      if (ierr /= 0) &
         & call gkw_abort('diagnos_timetrace :: could not allocate &
         & cbuf')
    end if
  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine read_last_data()

    if(.not.lcomplex_timetrace) return

  end subroutine read_last_data

  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : close_lu, ascii_fmt, binary_fmt

    if(.not.lcomplex_timetrace) return

    ! deallocate all arrays of this diagnostic
    if(allocated(cbuf)) deallocate(cbuf)
    ! be nice and close all logical units
    call close_lu(lun_csignal, binary_fmt)

  end subroutine finalize

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()

    if(.not.lcomplex_timetrace) return

  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_smallstep(fdis)
    use control, only : naverage
    use grid, only : nmod, ns, nx, ls, lrx, proc_subset
    use dist, only : fdisi, number_of_fields
    use mpiinterface, only : mpireduce_sum_inplace, root_processor
    use mpicomms, only : COMM_S_NE_X_NE
    complex, intent(in) :: fdis(:)
    
    complex :: values(nmod)
    integer :: idx, j, imod, ipar
    integer, save :: buf_index = 0

    if(.not.lcomplex_timetrace) return

    ! get a timeseries, per ky mode, similar to what is used in the
    ! computations in diagnos_growth_freq. Here, we just output a
    ! complex signal which can be analysed outside of GKW.
    values = 0.0
    idx = 0
    do j = 1, number_of_fields*nx
      do imod = 1, nmod
        do ipar = 1, ns
          idx=idx+1
          values(imod) = values(imod) + fdisi(lindx(idx))
        end do
      end do
    end do

    call mpireduce_sum_inplace(values, shape(values), COMM_S_NE_X_NE)

    if(root_processor) then
      buf_index = modulo(buf_index, naverage)+1
      cbuf(:, buf_index) = values
    end if

  end subroutine calc_smallstep

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()
    use io, only : append_chunk, xy_fmt, binary_fmt
    use control, only : naverage
    use mpiinterface, only : root_processor
    integer :: n

    if(.not.lcomplex_timetrace) return

    if(root_processor) then
      do n = 1, naverage
        call append_chunk(lun_csignal, &
             & cbuf(:,n), &
             & xy_fmt, binary_fmt)
      end do
    end if

  end subroutine output

end module diagnos_timetrace
