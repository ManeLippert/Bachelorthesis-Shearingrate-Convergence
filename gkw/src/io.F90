!===========================================================================
!> This module provides the abstracted interface for data output.
!> It is an interface wrapper around swappable modules io_* which 
!> perform the actual file handling for different output options.
!>
!> The recognized output options (selected by input variable io_format) are
!> 'none'    : do not output data to files at all
!> 'ascii'   : write *all* data in ascii format
!> 'binary'  : write *all* data in binary format
!> 'hdf5'    : write *all* data in binary format
!> 'mixed'   : write some data in ascii and some in binary format.
!>             Diagnostics can specify which format they want to write
!>             in the 'mixed' case with the preference_if_mixed argument.
!> 'hdf5+ascii' : Write the same file as with 'hdf5', but additionally
!>                write low dimensional data to ascii files.
!===========================================================================
module io
  use global, only : lenswitch, gkw_throw_abort, gkw_throw_warn, gkw_check_err
  use global, only : gkw_warn_any_proc
  ! use format strings provided by io_ascii, so that diagnostics can
  ! use them from here.
  use io_ascii, only : xy_fmt, xy_fmt_long

  ! ToDo:
  !  * Make sure the hdf5 file is properly finalized if GKW halts
  !    preliminarily.
  !  * maybe name the routines in io_hdf5 not specially but like the
  !    others; maybe say 'record' instead of 'chunk'
  !  * Think about moving the restart IO to the io_* modules -
  !    probably one needs to get rid of the IO dependency of the MPI
  !    wrappers, to make this possible

  implicit none

  private

  ! the following line makes it possible to write 'use io, only : xy_fmt'
  public :: xy_fmt, xy_fmt_long

  !> Defines the output format (options described in module description)
  character (len=lenswitch), save :: io_format

  !> IMPORTANT: Note that at the moment the diagnostic itself has to make
  !>            sure that always the same preference_if_mixed argument is
  !>            used for each logical unit.
  !> The following io_format names can be used by diagnostics to announce
  !> how they want to output their data in case of 'mixed'.
  character(len=lenswitch), parameter, public :: ascii_fmt = "ascii"
  character(len=lenswitch), parameter, public :: binary_fmt = "binary"
  character(len=lenswitch), parameter, public :: hdf5_fmt = "hdf5"
  !> The io_format='mixed' and 'hdf5+ascii' are special: one can set
  !> them in the input file, but no diagnostic should
  !> make any use of these tags.
  !> These are only public because they need to be tested exceptionally when
  !>  MPI-parallel and serial writing is distinguished, as in diagnos_fields.
  character(len=lenswitch), parameter, public :: mixed_fmt = "mixed"
  character(len=lenswitch), parameter, public :: hdf5_mixed_fmt = "hdf5+ascii"
  character(len=lenswitch), parameter, public :: none_fmt = "none"

  !> a switch, to disable output at some remaining places, were
  !> historic/debugging output is done using Fortran formatted output
  !> (i.e. one file per quantity)
  logical, save, public :: output_enabled


  !> Keys for the recommended metadata which every logical unit should
  !> be equipped with. If a value is not yet known, please give it the
  !> not_avail value.
  character(len=11), parameter, public :: description_key = "description"
  character(len=8),  parameter, public :: comments_key = "comments"
  character(len=13), parameter, public :: phys_unit_key = "physical unit"
  character(len=4),  parameter, public :: not_avail = "n.a."

  !> A switch which enables rounding of small real numbers to 0.0 .
  !> This is useful for testcases.
  logical, save :: io_testdata

  !> if MPI IO does not fully support derived datatypes, this should be T
  logical, save, public :: lmpi_broken_io = .false.

  ! utility routine, only made public to fix problems with +-0.0 and tiny
  ! number in particular files like kx_connect.dat and parfun.dat.
  public :: clean0r, clean0c

  public :: init_io, finalize_io
  public :: finalize_and_abort, finalize_and_exit
  public :: flush_all_open_files
  public :: lu_exists
  ! only temporarily a public routine, as long as there is no better solution:
  public :: get_free_file_unit
  
  !> * The basic item of this module is the 'Logical Unit' (lu)
  !> * A logical unit may refer to an ascii file, a simple binary file or
  !>   some structure in a more complicated format, for example
  !>   an HDF5 dataset.
  !> * A logical unit may even refer to 2 things (see the
  !>   duplicated_lus list) Then an output call produces e.g. both
  !>   a file and a HDF5 dataset.
  !> * A logical unit is identified by its 'Logical Unit Number' (lun).
  !> * A logical unit number is simply an integer number.

  ! subroutines to open, write (as often as desired) and close
  ! extendible datasets:
  public :: open_real_lu, open_complex_lu
  public :: close_lu, close_all_lus
  public :: append_chunk, seek_to_chunk
  public :: read_last_chunk
  public :: attach_metadata

  ! subroutines which open, write (once and for all) and close fixed
  ! size datasets:
  public :: output_array
  public :: mpi_output_array

  ! another subroutine which deals with metadata:
  public :: write_run_parameter

  ! subroutines which deal with logging:
  ! (not implemented)

  logical, save :: is_initialized = .false.
  integer, save :: proc_number
  logical, save :: is_restarted

  integer, parameter :: max_duplicated_lus = 1024
  integer, parameter :: UNUSED_LUN = -1
  integer, dimension(max_duplicated_lus, 2) :: duplicated_lus = UNUSED_LUN
  
  interface output_array
    module procedure output_array_1d_real
    module procedure output_array_2d_real
    module procedure output_array_3d_real
    module procedure output_array_4d_real
    module procedure output_array_5d_real
    module procedure output_array_6d_real
    module procedure output_array_1d_complex
    module procedure output_array_2d_complex
    module procedure output_array_3d_complex
    module procedure output_array_6d_complex
  end interface output_array

  interface append_chunk
    module procedure append_chunk_1d_real
    module procedure append_chunk_2d_real
    module procedure append_chunk_3d_real
    module procedure append_chunk_4d_real
    module procedure append_chunk_1d_complex
    module procedure append_chunk_2d_complex
    module procedure append_chunk_3d_complex
    module procedure append_chunk_4d_complex
    module procedure append_chunk_6d_complex
  end interface append_chunk

  interface read_last_chunk
    module procedure read_last_chunk_1d_real
  end interface read_last_chunk

  interface attach_metadata
    module procedure attach_metadata_string
    module procedure attach_metadata_real
    module procedure attach_metadata_integer
    module procedure attach_metadata_logical
    module procedure attach_metadata_r_array
    module procedure attach_metadata_i_array
    module procedure attach_metadata_string_name
    module procedure attach_metadata_real_name
    module procedure attach_metadata_integer_name
    module procedure attach_metadata_logical_name
    module procedure attach_metadata_r_array_name
    module procedure attach_metadata_i_array_name
  end interface attach_metadata
  
  interface write_run_parameter
    module procedure write_run_param_string
    module procedure write_run_param_real
    module procedure write_run_param_complex
    module procedure write_run_param_integer
    module procedure write_run_param_logical
    module procedure write_run_param_r_array
    module procedure write_run_param_i_array
  end interface write_run_parameter

  interface mpi_output_array
    module procedure mpi_output_array_1d_real
    module procedure mpi_output_array_3d_real
    module procedure mpi_output_array_3d_complex
    module procedure mpi_output_array_6d_real
  end interface mpi_output_array

  interface mpi_serial_output_array
    module procedure mpi_serial_output_3d_real
    module procedure mpi_serial_output_3d_cmplx
    module procedure mpi_serial_output_6d_real
  end interface mpi_serial_output_array

  ! These have been added to avoid warnings of the type "'xyz' defined but not used".
  public roundr

contains

  !---------------------------------------------------------------------------
  !> round numbers close to zero to a clean positive zero for the test cases.
  !> This function works also on arrays.
  !---------------------------------------------------------------------------
  elemental function clean0r(input) result(output)
    real, intent(in) :: input
    real :: output
    ! cast all values to zero whose abs value is smaller than some tiny number.
    output = merge(input,+0.0,abs(input) > 1.e-13 .or. .not.io_testdata)

    ! alternative: round to a certain precision
    !output = merge(input,anint(input*1.0e10)*1.0e-10,.not.io_testdata)
  end function clean0r
  elemental function clean0c(input) result(output)
    complex, intent(in) :: input
    complex :: output
    output = cmplx(clean0r(real(input)), &
       & clean0r(aimag(input)))
  end function clean0c
  
  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  elemental function roundr(input) result(output)
    real, intent(in) :: input
    real :: output
    ! alternative: round to a certain precision
    output = merge(input,anint(input*1.0e10)*1.0e-10,.not.io_testdata)
  end function roundr


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine init_io(proc_number_, io_format_, io_testdata_, is_restarted_, &
     & io_legacy)
    use io_ascii, only : init_ascii => init
    use io_binary, only : init_binary => init
#if HAVE_HDF5
    use io_hdf5, only : init_hdf5 => init
#endif
    integer, intent(in) :: proc_number_
    character (len=*), intent(in) :: io_format_
    logical, intent(in) :: io_testdata_
    logical, intent(in) :: is_restarted_
    logical, intent(in) :: io_legacy

    proc_number = proc_number_
    io_format = io_format_
    io_testdata = io_testdata_

    ! the IO interfaces need to know whether they shall
    !   a) overwrite existing logical units [is_restarted == .false.]
    !   b) append to existing logical units [is_restarted == .true. ]
    is_restarted = is_restarted_

    select case(io_format)
    case(none_fmt)
      call init_ascii(proc_number, is_restarted, .false., io_legacy)
    case(ascii_fmt)
      call init_ascii(proc_number, is_restarted, .true., io_legacy)
    case(binary_fmt, mixed_fmt)
      call init_ascii(proc_number, is_restarted, .true., io_legacy)
      call init_binary(proc_number, is_restarted, io_legacy)
#if HAVE_HDF5
    case(hdf5_mixed_fmt)
      call init_hdf5(proc_number, is_restarted)
      call init_ascii(proc_number, is_restarted, .true., io_legacy)
    case(hdf5_fmt)
      call init_hdf5(proc_number, is_restarted)
      call init_ascii(proc_number, is_restarted, .false., io_legacy)
#else
    case(hdf5_fmt)
      call gkw_throw_abort('GKW was not compiled with HDF5, io_format="hdf5" cannot be chosen.')
    case (hdf5_mixed_fmt)
      call gkw_throw_abort('GKW was not compiled with HDF5, io_format="hdf5+ascii" cannot be chosen.')
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('init_io')

    output_enabled = (io_format /= none_fmt)

    is_initialized = .true.
  end subroutine init_io


  !---------------------------------------------------------------------------
  !> This routine closes the concrete IO interface(s). It must be called
  !> after close_all_lus().
  !---------------------------------------------------------------------------
  subroutine finalize_io()
    use io_ascii, only : finalize_ascii => finalize
    use io_binary, only : finalize_binary => finalize
#if HAVE_HDF5
    use io_hdf5, only : finalize_hdf5 => finalize
#endif

    if(.not. is_initialized) return
    
    call close_all_lus()

    select case(io_format)
    case(none_fmt)
      call finalize_ascii()
    case(ascii_fmt)
      call finalize_ascii()
    case(binary_fmt)
      call finalize_binary()
    case(mixed_fmt)
      call finalize_ascii()
      call finalize_binary()
#if HAVE_HDF5
    case(hdf5_fmt)
      ! but at the moment we still need some normal file units, too,
      ! therefore:
      call finalize_ascii()
      call finalize_hdf5()
    case (hdf5_mixed_fmt)
      call finalize_ascii()
      call finalize_hdf5()
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    is_initialized = .false.

    ! I cannot call abort here, because abort calls this and then I land in
    ! an infinite loop...
    if(gkw_check_err()) then
      !at least warn
      call gkw_warn_any_proc('A critical error occured in finalize_io().')
    end if
  end subroutine finalize_io

  !---------------------------------------------------------------------------
  !> This routine is used if an abort is desired but one cannot make sure
  !> that all processors call it.
  !---------------------------------------------------------------------------
  subroutine finalize_and_abort(abort_message)
    use mpiinterface, only : mpiabort

    character (len=*),intent(in), optional :: abort_message

    ! SRG At the moment I see no way to properly close files if an
    ! arbitrary processor calls this. This only really closes files if the
    ! root processor does it:
    call finalize_io()
    
    ! Quit safely without MPI deadlock
    call mpiabort(abort_message)

  end subroutine finalize_and_abort

  !---------------------------------------------------------------------------
  !> This routine is to be called by *all processors* and used if a run is
  !> to be ended (preliminarily).
  !---------------------------------------------------------------------------
  subroutine finalize_and_exit()
    use mpiinterface, only : mpifinalize

    call finalize_io()

    call mpifinalize()
    stop 1
  end subroutine finalize_and_exit

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine flush_all_open_files()
    use io_ascii, only : flush_all_open_files_ascii => flush_all_open_files
    use io_binary, only : flush_all_open_files_binary => flush_all_open_files
#if HAVE_HDF5
    use io_hdf5, only : flush_all_open_files_hdf5 => flush_all_open_files
#endif
    use mpiinterface, only : root_processor

    if(.not. root_processor) return

    select case(io_format)
    case(none_fmt)
      call flush_all_open_files_ascii()
    case(ascii_fmt)
      call flush_all_open_files_ascii()
    case(binary_fmt)
      call flush_all_open_files_binary()
    case(mixed_fmt)
      call flush_all_open_files_ascii()
      call flush_all_open_files_binary()
#if HAVE_HDF5
    case(hdf5_fmt)
      call flush_all_open_files_hdf5()
    case(hdf5_mixed_fmt)
      call flush_all_open_files_hdf5()
      call flush_all_open_files_ascii()
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('flush_all_open_files')
  end subroutine flush_all_open_files

  !---------------------------------------------------------------------------
  !> This routine provides a file_unit number. This is the kind of number
  !> which is used in Fortran to access files. Note that this different from a
  !> logical unit number (lun), in terms of this module. A lun is one
  !> abstraction layer higher. A lun refers to a file unit in this IO module,
  !> but in another IO module, it may refer e.g. to a dataset ID.
  !---------------------------------------------------------------------------
  subroutine get_free_file_unit(file_unit)
    use io_ascii, only : get_free_file_unit_ascii => get_free_file_unit
    use io_binary, only : get_free_file_unit_binary => get_free_file_unit

    integer, intent(out) :: file_unit

    select case(io_format)
    case(none_fmt)
      call get_free_file_unit_ascii(file_unit)
    case(ascii_fmt)
      call get_free_file_unit_ascii(file_unit)
    case(binary_fmt)
      call get_free_file_unit_binary(file_unit)
    case(mixed_fmt)
      call get_free_file_unit_ascii(file_unit)
#if HAVE_HDF5
    case(hdf5_fmt)
      call get_free_file_unit_ascii(file_unit)
    case(hdf5_mixed_fmt)
      call get_free_file_unit_ascii(file_unit)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('get_free_file_unit')
  end subroutine get_free_file_unit

  !---------------------------------------------------------------------------
  !> return corresponding other ("duplicated") logical unit number, by
  !> looking it up in the duplicated_lus list.
  !---------------------------------------------------------------------------
  function get_duplicate_lun(lun)
    integer, intent(in) :: lun
    integer :: get_duplicate_lun

    integer :: i

    get_duplicate_lun = -1
    
    do i = 1, size(duplicated_lus, 1)
      if(duplicated_lus(i,1) == lun) then
        get_duplicate_lun = duplicated_lus(i,2)
        return
      end if
    end do
  end function get_duplicate_lun

#if HAVE_HDF5
  !---------------------------------------------------------------------------
  !> put an entry into the duplicated_lus list, with a pair of lus which
  !> are used for the same quantity but go into two different formats.
  !---------------------------------------------------------------------------
  subroutine set_duplicate_lun(lun1, lun2)
    integer, intent(in) :: lun1,lun2
    integer :: i

    ! write (*,*) "HDF5 lun ", lun1, " is duplicated to ascii lun ", lun2
    
    do i = 1, size(duplicated_lus, 1)
      if(duplicated_lus(i,1) == UNUSED_LUN) then
        duplicated_lus(i,1) = lun1
        duplicated_lus(i,2) = lun2
        return
      end if
    end do
    ! if program flow gets here this means that the array must be set bigger
    call gkw_throw_abort('not enough logical unit numbers: &
       & increase max_duplicated_lus in module io: ')
  end subroutine set_duplicate_lun

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine remove_duplicate_lun(lun)
    integer, intent(in) :: lun
    integer :: i

    do i = 1, size(duplicated_lus, 1)
      if(duplicated_lus(i,1) == lun) then
        duplicated_lus(i,1) = UNUSED_LUN
        return
      end if
    end do
  end subroutine remove_duplicate_lun
#endif

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine open_real_lu(luname, groupname, hyperslab_dims, &
     & preference_if_mixed, lun, hyperslab_labels)
    use io_ascii, only : open_real_lu_ascii => open_real_lu
    use io_binary, only : open_real_lu_binary => open_real_lu
#if HAVE_HDF5
    use io_hdf5, only : open_real_lu_hdf5 => open_real_dataset_hdf5
    integer :: duplicate_lun
#endif
    character (len=*), intent(in) :: luname, groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    character (len=*), intent(in) :: preference_if_mixed
    integer, intent(out) :: lun
    character (len=*), dimension(:), intent(in), optional :: hyperslab_labels
#define NONE_LUN 99999

    select case(io_format)
    case(none_fmt)
      lun = NONE_LUN
    case(ascii_fmt)
      call open_real_lu_ascii(luname, groupname, hyperslab_dims, lun)
    case(binary_fmt)
      call open_real_lu_binary(luname, groupname, hyperslab_dims, lun)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call open_real_lu_ascii(luname, groupname, hyperslab_dims, lun)
      case(binary_fmt)
        call open_real_lu_binary(luname, groupname, hyperslab_dims, lun)
      case default
        call gkw_throw_abort("open_real_lu: unknown value for preference format: '"&
           & //preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      if(present(hyperslab_labels)) then
        call open_real_lu_hdf5(luname, groupname, hyperslab_dims, lun, &
           & hyperslab_labels)
      else
        call open_real_lu_hdf5(luname, groupname, hyperslab_dims, lun)
      end if
    case (hdf5_mixed_fmt)
      if(present(hyperslab_labels)) then
        call open_real_lu_hdf5(luname, groupname, hyperslab_dims, lun, &
           & hyperslab_labels)
      else
        call open_real_lu_hdf5(luname, groupname, hyperslab_dims, lun)
      end if

      ! if this is a lowdimensional dataset, then duplicate it to an ascii file
      if(size(hyperslab_dims) == 1) then
        call open_real_lu_ascii(luname, groupname, &
           & hyperslab_dims, duplicate_lun)
        ! Put an entry to the list of duplicated lus.
        call set_duplicate_lun(lun, duplicate_lun)
      end if
#else
    if (.false.) write(*,*) hyperslab_labels
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('open_real_lu')
  end subroutine open_real_lu

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine open_complex_lu(luname, groupname, hyperslab_dims, &
     & preference_if_mixed, lun)
    use io_ascii, only : open_complex_lu_ascii => open_complex_lu
    use io_binary, only : open_complex_lu_binary => open_complex_lu
#if HAVE_HDF5
    use io_hdf5, only : open_complex_lu_hdf5 => open_complex_dataset_hdf5
    integer :: duplicate_lun
#endif

    character (len=*), intent(in) :: luname, groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    character (len=*), intent(in) :: preference_if_mixed
    integer, intent(out) :: lun

    select case(io_format)
    case(none_fmt)
      lun = NONE_LUN
    case(ascii_fmt)
      call open_complex_lu_ascii(luname, groupname, hyperslab_dims, lun)
    case(binary_fmt)
      call open_complex_lu_binary(luname, groupname, hyperslab_dims, lun)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call open_complex_lu_ascii(luname, groupname, hyperslab_dims, lun)
      case(binary_fmt)
        call open_complex_lu_binary(luname, groupname, hyperslab_dims, lun)
      case default
        call gkw_throw_abort("open_complex_lu: unknown value for preference format: '"//&
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call open_complex_lu_hdf5(luname, groupname, hyperslab_dims, lun)
    case(hdf5_mixed_fmt)
      call open_complex_lu_hdf5(luname, groupname, hyperslab_dims, lun)
      ! if this is a lowdimensional dataset, then duplicate it to an ascii file
      if(size(hyperslab_dims) == 1) then
        call open_complex_lu_ascii(luname, groupname, hyperslab_dims, duplicate_lun)
        ! Put an entry to the list of duplicated lus.
        call set_duplicate_lun(lun, duplicate_lun)
      end if
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('open_complex_lu')
  end subroutine open_complex_lu

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine close_lu(lun, preference_if_mixed)
    use io_ascii, only : close_lu_ascii => close_lu
    use io_binary, only : close_lu_binary => close_lu
#if HAVE_HDF5
    use io_hdf5, only :  close_lu_hdf5 => close_dataset_hdf5
    integer :: duplicate_lun
#endif

    integer, intent(in) :: lun
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call close_lu_ascii(lun)
    case(binary_fmt)
      call close_lu_binary(lun)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call close_lu_ascii(lun)
      case(binary_fmt)
        call close_lu_binary(lun)
      case default
        call gkw_throw_abort("close_lu: unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call close_lu_hdf5(lun)
    case(hdf5_mixed_fmt)
      call close_lu_hdf5(lun)
      ! Find out if the given lu is duplicaded to another format.
      duplicate_lun = get_duplicate_lun(lun)
      if(duplicate_lun > 0) then
        ! Close also that other lu.
        call close_lu_ascii(duplicate_lun)
        ! Remove the entry from the list of duplicated lus.
        call remove_duplicate_lun(lun)
      end if
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select
    
    if(gkw_check_err()) call finalize_and_abort('close_lu')
  end subroutine close_lu

  !---------------------------------------------------------------------------
  !> This routine closes all open logical units, so that a subsequent call
  !> of finalize_io() leaves everything in a clean state.
  !---------------------------------------------------------------------------
  subroutine close_all_lus()
    use io_ascii, only : close_all_lus_ascii => close_all_lus
    use io_binary, only : close_all_lus_binary => close_all_lus
#if HAVE_HDF5
    use io_hdf5, only :  close_all_lus_hdf5 => close_all_lus
#endif
    use mpiinterface, only : root_processor

    if (.not.root_processor) return

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call close_all_lus_ascii()
    case(binary_fmt)
      call close_all_lus_binary()
    case(mixed_fmt)
      call close_all_lus_ascii()
      call close_all_lus_binary()
#if HAVE_HDF5
    case(hdf5_fmt)
      call close_all_lus_hdf5()
    case(hdf5_mixed_fmt)
      call close_all_lus_hdf5()
      call close_all_lus_ascii()
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    ! I cannot call abort here, because abort calls this and then I land in
    ! an infinite loop...
    if(gkw_check_err()) then
      !at least warn
      call gkw_warn_any_proc('A critical error occured in finalize_io().')
    end if
  end subroutine close_all_lus

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  function lu_exists(luname, groupname, preference_if_mixed)
    use io_ascii, only : lu_exists_ascii => lu_exists
    use io_binary, only : lu_exists_binary => lu_exists
#if HAVE_HDF5
    use io_hdf5, only :  lu_exists_hdf5 => dataset_exists_hdf5
#endif
    use mpiinterface, only : mpibcast

    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: preference_if_mixed
    logical :: lu_exists

    select case(io_format)
    case(none_fmt)
      lu_exists = .false.
    case(ascii_fmt)
      lu_exists = lu_exists_ascii(luname, groupname)
    case(binary_fmt)
      lu_exists = lu_exists_binary(luname, groupname)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        lu_exists = lu_exists_ascii(luname, groupname)
      case(binary_fmt)
        lu_exists = lu_exists_binary(luname, groupname)
      case default
        call gkw_throw_abort("lu_exists: unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      lu_exists = lu_exists_hdf5(luname, groupname)
    case(hdf5_mixed_fmt)
      lu_exists = lu_exists_hdf5(luname, groupname) .or. &
         & lu_exists_ascii(luname, groupname)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select
    
    ! all processes get the result of the root process, to make sure
    ! this function returns the same value for all processes.
    call mpibcast(lu_exists,1)

    if(gkw_check_err()) call finalize_and_abort('lu_exists')
  end function lu_exists

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_1d_real(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use mpiinterface, only : root_processor
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
    integer :: duplicate_lun
#endif
    real, dimension(:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    ! for PGI compiler reshaping bug issue 214#50
    ! write(*,*) 'data shape', shape(data)

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0r(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0r(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0r(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0r(data))
      case default
        call gkw_throw_abort("append_chunk_1d_real: unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0r(data))
    case(hdf5_mixed_fmt)
      ! While 2D chunks are output only in HDF5, 1D chunks are
      ! duplicated to textfiles.
      call append_chunk_hdf5(lun, clean0r(data))
      duplicate_lun = get_duplicate_lun(lun)
      if(root_processor .and. duplicate_lun < 0) then
        call gkw_throw_abort('WHAT? Never ever!')
      else
        call append_chunk_ascii(duplicate_lun, clean0r(data), fmt)
      end if
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_1d_real')
  end subroutine append_chunk_1d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_2d_real(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
#endif

    real, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0r(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0r(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0r(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0r(data))
      case default
        call gkw_throw_abort("append_chunk_2d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0r(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0r(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_2d_real')
  end subroutine append_chunk_2d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_3d_real(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
#endif

    real, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0r(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0r(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0r(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0r(data))
      case default
        call gkw_throw_abort("append_chunk_3d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0r(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0r(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_3d_real')
  end subroutine append_chunk_3d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_4d_real(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
#endif

    real, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0r(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0r(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0r(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0r(data))
      case default
        call gkw_throw_abort("append_chunk_4d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0r(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0r(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_4d_real')
  end subroutine append_chunk_4d_real
  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_1d_complex(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use mpiinterface, only : root_processor
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
    integer :: duplicate_lun
#endif
    complex, dimension(:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0c(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0c(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0c(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0c(data))
      case default
        call gkw_throw_abort("append_chunk_1d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
      duplicate_lun = get_duplicate_lun(lun)
      if(root_processor .and. duplicate_lun < 0) then
        call gkw_throw_abort('WHAT? Never ever!')
      else
        call append_chunk_ascii(duplicate_lun, clean0c(data), fmt)
      end if
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_1d_complex')
  end subroutine append_chunk_1d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_2d_complex(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
#endif

    complex, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0c(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0c(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0c(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0c(data))
      case default
        call gkw_throw_abort("append_chunk_2d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_2d_complex')
  end subroutine append_chunk_2d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_3d_complex(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
#endif

    complex, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0c(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0c(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0c(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0c(data))
      case default
        call gkw_throw_abort("append_chunk_3d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_3d_complex')
  end subroutine append_chunk_3d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_4d_complex(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
#endif

    complex, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0c(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0c(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0c(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0c(data))
      case default
        call gkw_throw_abort("append_chunk_4d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_4d_complex')
  end subroutine append_chunk_4d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_6d_complex(lun, data, fmt, preference_if_mixed)
    use io_ascii, only : append_chunk_ascii => append_chunk
    use io_binary, only : append_chunk_binary => append_chunk
#if HAVE_HDF5
    use io_hdf5, only :  append_chunk_hdf5 => append_chunk
#endif

    complex, dimension(:,:,:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call append_chunk_ascii(lun, clean0c(data), fmt)
    case(binary_fmt)
      call append_chunk_binary(lun, clean0c(data))
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call append_chunk_ascii(lun, clean0c(data), fmt)
      case(binary_fmt)
        call append_chunk_binary(lun, clean0c(data))
      case default
        call gkw_throw_abort("append_chunk_6d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
    case(hdf5_mixed_fmt)
      call append_chunk_hdf5(lun, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('append_chunk_6d_complex')
  end subroutine append_chunk_6d_complex

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine read_last_chunk_1d_real(lun, fmt, data, nlast, preference_if_mixed)
    use io_ascii, only : read_last_chunk_ascii => read_last_chunk
    use io_binary, only : read_last_chunk_binary => read_last_chunk
#if HAVE_HDF5
    use io_hdf5, only :  read_last_chunk_hdf5 => read_last_chunk
#endif

    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt, preference_if_mixed
    real, dimension(:), intent(out) :: data
    integer, intent(out) :: nlast

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call read_last_chunk_ascii(lun, fmt, data, nlast)
    case(binary_fmt)
      call read_last_chunk_binary(lun, data, nlast)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call read_last_chunk_ascii(lun, fmt, data, nlast)
      case(binary_fmt)
        call read_last_chunk_binary(lun, data, nlast)
      case default
        call gkw_throw_abort("read_last_chunk_1d_real: unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call read_last_chunk_hdf5(lun, data, nlast)
    case(hdf5_mixed_fmt)
      call read_last_chunk_ascii(lun, fmt, data, nlast)
      if(gkw_check_err()) then
        ! then the file did probably not exist. try to read from the hdf5 data.
        call read_last_chunk_hdf5(lun, data, nlast)
      end if
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('read_last_chunk_1d_real')
  end subroutine read_last_chunk_1d_real

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine seek_to_chunk(lun, index, preference_if_mixed)
    use io_ascii, only : seek_to_chunk_ascii => seek_to_chunk
    use io_binary, only : seek_to_chunk_binary => seek_to_chunk
#if HAVE_HDF5
    use io_hdf5, only : seek_to_chunk_hdf5 => seek_to_chunk
#endif

    integer, intent(in) :: lun
    integer, intent(in) :: index
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call seek_to_chunk_ascii(lun, index)
    case(binary_fmt)
      call seek_to_chunk_binary(lun, index)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call seek_to_chunk_ascii(lun, index)
      case(binary_fmt)
        call seek_to_chunk_binary(lun, index)
      case default
        call gkw_throw_abort("seek_to_chunk: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call seek_to_chunk_hdf5(lun, index)
    case(hdf5_mixed_fmt)
      call seek_to_chunk_hdf5(lun, index)
      call seek_to_chunk_ascii(get_duplicate_lun(lun), index)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('seek_to_chunk')
  end subroutine seek_to_chunk


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_1d_real(luname, groupname, data, order, fmt, &
     & preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    real, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0r(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0r(data), order)
      case default
        call gkw_throw_abort("output_array_1d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_1d_real')
  end subroutine output_array_1d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_2d_real(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    real, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0r(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0r(data), order)
      case default
        call gkw_throw_abort("output_array_2d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_2d_real')
  end subroutine output_array_2d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_3d_real(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    real, dimension(:,:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0r(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0r(data), order)
      case default
        call gkw_throw_abort("output_array_3d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_3d_real')
  end subroutine output_array_3d_real

    !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_4d_real(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    real, dimension(:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0r(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0r(data), order)
      case default
        call gkw_throw_abort("output_array_4d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_4d_real')
  end subroutine output_array_4d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_5d_real(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    real, dimension(:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0r(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0r(data), order)
      case default
        call gkw_throw_abort("output_array_5d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_5d_real')
  end subroutine output_array_5d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_6d_real(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    real, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0r(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0r(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0r(data), order)
      case default
        call gkw_throw_abort("output_array_5d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0r(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_5d_real')
  end subroutine output_array_6d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_1d_complex(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    complex, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0c(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0c(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0c(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0c(data), order)
      case default
        call gkw_throw_abort("output_array_1d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_1d_complex')
  end subroutine output_array_1d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_2d_complex(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    complex, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0c(data), order,fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0c(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0c(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0c(data), order)
      case default
        call gkw_throw_abort("output_array_2d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format)
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_2d_complex')
  end subroutine output_array_2d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_3d_complex(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    complex, dimension(:,:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0c(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0c(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0c(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0c(data), order)
      case default
        call gkw_throw_abort("output_array_3d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_3d_complex')
  end subroutine output_array_3d_complex


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_array_6d_complex(luname, groupname, data, order, fmt, preference_if_mixed)
    use io_ascii, only : output_array_ascii => output_array
    use io_binary, only : output_array_binary => output_array
#if HAVE_HDF5
    use io_hdf5, only :  output_array_hdf5 => output_array
#endif

    complex, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order
    character (len=*), intent(in) :: preference_if_mixed

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      call output_array_ascii(luname, groupname, clean0c(data), order, fmt)
    case(binary_fmt)
      call output_array_binary(luname, groupname, clean0c(data), order)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        call output_array_ascii(luname, groupname, clean0c(data), order, fmt)
      case(binary_fmt)
        call output_array_binary(luname, groupname, clean0c(data), order)
      case default
        call gkw_throw_abort("output_array_3d_complex: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
    case(hdf5_mixed_fmt)
      call output_array_hdf5(luname, groupname, clean0c(data))
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('output_array_3d_complex')
  end subroutine output_array_6d_complex


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_string(lun, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt, binary_fmt, mixed_fmt)
        call attach_metadata_ascii(lun, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(lun, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(lun, key, val)
        call attach_metadata_ascii(get_duplicate_lun(lun), key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_string')
  end subroutine attach_metadata_string

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_real(lun, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    real, intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt, binary_fmt, mixed_fmt)
        call attach_metadata_ascii(lun, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(lun, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(lun, key, val)
        call attach_metadata_ascii(get_duplicate_lun(lun), key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_real')
  end subroutine attach_metadata_real

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_integer(lun, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    integer, intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt, binary_fmt, mixed_fmt)
        call attach_metadata_ascii(lun, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(lun, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(lun, key, val)
        call attach_metadata_ascii(get_duplicate_lun(lun), key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_integer')
  end subroutine attach_metadata_integer


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_logical(lun, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    logical, intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt, binary_fmt, mixed_fmt)
        call attach_metadata_ascii(lun, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(lun, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(lun, key, val)
        call attach_metadata_ascii(get_duplicate_lun(lun), key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_logical')
  end subroutine attach_metadata_logical


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_r_array(lun, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt, binary_fmt, mixed_fmt)
        call attach_metadata_ascii(lun, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(lun, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(lun, key, val)
        call attach_metadata_ascii(get_duplicate_lun(lun), key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_r_array')
  end subroutine attach_metadata_r_array


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_i_array(lun, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt, binary_fmt, mixed_fmt)
        call attach_metadata_ascii(lun, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(lun, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(lun, key, val)
        call attach_metadata_ascii(get_duplicate_lun(lun), key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_i_array')
  end subroutine attach_metadata_i_array


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_string_name(luname, groupname, key, val, &
     & preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt,binary_fmt,mixed_fmt)
        call attach_metadata_ascii(luname, groupname, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
        call attach_metadata_ascii(luname, groupname, key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_string_name')
  end subroutine attach_metadata_string_name

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_real_name(luname, groupname, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    real, intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt,binary_fmt,mixed_fmt)
        call attach_metadata_ascii(luname, groupname, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
        call attach_metadata_ascii(luname, groupname, key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_real_name')
  end subroutine attach_metadata_real_name

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_integer_name(luname, groupname, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    integer, intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt,binary_fmt,mixed_fmt)
        call attach_metadata_ascii(luname, groupname, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
        call attach_metadata_ascii(luname, groupname, key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_integer_name')
  end subroutine attach_metadata_integer_name

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_logical_name(luname, groupname, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    logical, intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt,binary_fmt,mixed_fmt)
        call attach_metadata_ascii(luname, groupname, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
        call attach_metadata_ascii(luname, groupname, key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_logical_name')
  end subroutine attach_metadata_logical_name

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_r_array_name(luname, groupname, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt,binary_fmt,mixed_fmt)
        call attach_metadata_ascii(luname, groupname, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
        call attach_metadata_ascii(luname, groupname, key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_r_array_name')
  end subroutine attach_metadata_r_array_name


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_i_array_name(luname, groupname, key, val, preference_if_mixed)
    use io_ascii, only : attach_metadata_ascii => attach_metadata
#if HAVE_HDF5
    use io_hdf5, only : attach_metadata_hdf5 => attach_metadata
#endif
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val
    character (len=*), intent(in) :: preference_if_mixed

    ! Silence compiler until it is decided what to do with this argument
    if (len(preference_if_mixed) > 0) continue

    select case(io_format)
      case(none_fmt)
        continue
      case(ascii_fmt,binary_fmt,mixed_fmt)
        call attach_metadata_ascii(luname, groupname, key, val)
#if HAVE_HDF5
      case(hdf5_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
      case(hdf5_mixed_fmt)
        call attach_metadata_hdf5(luname, groupname, key, val)
        call attach_metadata_ascii(luname, groupname, key, val)
#endif
    end select
    if(gkw_check_err()) call finalize_and_abort('attach_metadata_i_array_name')
  end subroutine attach_metadata_i_array_name

  
  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_string(group, key, val)
    use mpiinterface, only : check_equal_on_all_processors
#if HAVE_HDF5
    use io_hdf5, only : write_run_parameter_hdf5 => write_run_parameter
#endif
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val

    select case(io_format)
#if HAVE_HDF5
    case(hdf5_fmt)
      call write_run_parameter_hdf5(group, key, val)
    case(hdf5_mixed_fmt)
      call write_run_parameter_hdf5(group, key, val)
#endif
    end select

    if(gkw_check_err()) call finalize_and_abort('write_run_param_string')
    if(.not.check_equal_on_all_processors(val)) then
      write (*,*) group, key, val
      call finalize_and_abort(key//' in '//group//' is not equal on all &
         & processes, has not been broadcast.')
    end if
  end subroutine write_run_param_string


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_real(group, key, val)
    use mpiinterface, only : check_equal_on_all_processors
#if HAVE_HDF5
    use io_hdf5, only : write_run_parameter_hdf5 => write_run_parameter
#endif
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    real, intent(in) :: val

    select case(io_format)
#if HAVE_HDF5
    case(hdf5_fmt)
      call write_run_parameter_hdf5(group, key, val)
    case(hdf5_mixed_fmt)
      call write_run_parameter_hdf5(group, key, val)
#endif
    end select

    if(gkw_check_err()) call finalize_and_abort('write_run_param_real')
    if(.not.check_equal_on_all_processors(val)) then
      write (*,*) group, key, val
      call finalize_and_abort(key//' in '//group//' is not equal on all &
         & processes, has not been broadcast.')
    end if
  end subroutine write_run_param_real


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_complex(group, key, val)
    use mpiinterface, only : check_equal_on_all_processors
#if HAVE_HDF5
    use io_hdf5, only : write_run_parameter_hdf5 => write_run_parameter
#endif
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    complex, intent(in) :: val

    select case(io_format)
#if HAVE_HDF5
    case(hdf5_fmt)
      call write_run_parameter_hdf5(group, key, val)
    case(hdf5_mixed_fmt)
      call write_run_parameter_hdf5(group, key, val)
#endif
    end select

    if(gkw_check_err()) call finalize_and_abort('write_run_param_complex')
    if(.not.check_equal_on_all_processors(val)) then
      write (*,*) group, key, val
      call finalize_and_abort(key//' in '//group//' is not equal on all &
         & processes, has not been broadcast.')
    end if
  end subroutine write_run_param_complex


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_integer(group, key, val)
    use mpiinterface, only : check_equal_on_all_processors
#if HAVE_HDF5
    use io_hdf5, only : write_run_parameter_hdf5 => write_run_parameter
#endif
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    integer, intent(in) :: val

    select case(io_format)
#if HAVE_HDF5
    case(hdf5_fmt)
      call write_run_parameter_hdf5(group, key, val)
    case(hdf5_mixed_fmt)
      call write_run_parameter_hdf5(group, key, val)
#endif
    end select

    if(gkw_check_err()) call finalize_and_abort('write_run_param_integer')
    if(.not.check_equal_on_all_processors(val)) then
      write (*,*) group, key, val
      call finalize_and_abort(key//' in '//group//' is not equal on all &
         & processes, has not been broadcast.')
    end if
  end subroutine write_run_param_integer


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_logical(group, key, val)
    use mpiinterface, only : check_equal_on_all_processors
#if HAVE_HDF5
    use io_hdf5, only : write_run_parameter_hdf5 => write_run_parameter
#endif
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    logical, intent(in) :: val

    select case(io_format)
#if HAVE_HDF5
    case(hdf5_fmt)
      call write_run_parameter_hdf5(group, key, val)
    case(hdf5_mixed_fmt)
      call write_run_parameter_hdf5(group, key, val)
#endif
    end select

    if(gkw_check_err()) call finalize_and_abort('write_run_param_logical')
    if(.not.check_equal_on_all_processors(val)) then
      write (*,*) group, key, val
      call finalize_and_abort(key//' in '//group//' is not equal on all &
         & processes, has not been broadcast.')
    end if
  end subroutine write_run_param_logical


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_r_array(group, key, val)
    use mpiinterface, only : check_equal_on_all_processors
#if HAVE_HDF5
    use io_hdf5, only : write_run_parameter_hdf5 => write_run_parameter
#endif
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val

    select case(io_format)
#if HAVE_HDF5
    case(hdf5_fmt)
      call write_run_parameter_hdf5(group, key, val)
    case(hdf5_mixed_fmt)
      call write_run_parameter_hdf5(group, key, val)
#endif
    end select

    if(gkw_check_err()) call finalize_and_abort('write_run_param_r_array')
    if(.not.check_equal_on_all_processors(val)) then
      write (*,*) group, key, val
      call finalize_and_abort(key//' in '//group//' is not equal on all &
         & processes, has not been broadcast.')
    end if
  end subroutine write_run_param_r_array


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_i_array(group, key, val)
    use mpiinterface, only : check_equal_on_all_processors
#if HAVE_HDF5
    use io_hdf5, only : write_run_parameter_hdf5 => write_run_parameter
#endif
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    select case(io_format)
#if HAVE_HDF5
    case(hdf5_fmt)
      call write_run_parameter_hdf5(group, key, val)
    case(hdf5_mixed_fmt)
      call write_run_parameter_hdf5(group, key, val)
#endif
    end select

    if(gkw_check_err()) call finalize_and_abort('write_run_param_i_array')
    if(.not.check_equal_on_all_processors(val)) then
      write (*,*) group, key, val
      call finalize_and_abort(key//' in '//group//' is not equal on all &
         & processes, has not been broadcast.')
    end if
  end subroutine write_run_param_i_array


  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine mpi_output_array_1d_real(luname, groupname, array, mpi_dtype, &
       & comm, fmt, preference_if_mixed)
    use io_binary, only : binary_mpi_output_array => mpi_output_array

    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: mpi_dtype
    real, dimension(:), intent(in) :: array
    integer, intent(in), optional :: comm

    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=*), intent(in) :: preference_if_mixed

    if (.false.) write(*,*) fmt, ' ', groupname

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      !ascii mpi-io is not possible, do it serially

      !or do it binary:
      call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
         & comm)
    case(binary_fmt)
      call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
         & comm)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        !ascii mpi-io is not possible, do it serially

        !or do it binary:
        call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
           & comm)
      case(binary_fmt)
        call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
           & comm)
      case default
        call gkw_throw_abort( &
     & "mpi_output_array_1d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      !parallel hdf5 io is not implemented yet, do it serially

      !or do it binary:
      call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
         & comm)
    case(hdf5_mixed_fmt)
      !parallel hdf5 io is not implemented yet, do it serially

      !or do it binary:
      call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
         & comm)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('mpi_output_array_1d_real')
  end subroutine mpi_output_array_1d_real

  !--------------------------------------------------------------------------
  !> should be called by *all* processes of COMM_CART.  (although at
  !> the moment just the procs working on the slab (plus or including
  !> root) are sufficient)
  !--------------------------------------------------------------------------
  subroutine mpi_output_array_3d_real(luname, groupname, array, mpi_dtype, &
     & buf_global, &
     & comm, fmt, preference_if_mixed, proc_is_in_slab, root_is_in_comm, tag)
    use io_binary, only : binary_mpi_output_array => mpi_output_array
    
    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: mpi_dtype
    real, dimension(:,:,:), intent(in) :: array
    real, dimension(:,:,:), intent(out) :: buf_global
    integer, intent(in) :: comm, tag

    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=*), intent(in) :: preference_if_mixed
    logical, intent(in) :: proc_is_in_slab, root_is_in_comm

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      !ascii mpi-io is not possible, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
           & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
           & root_is_in_comm, tag)
    case(binary_fmt)
      ! do truely parallel binary MPI-IO
      if(proc_is_in_slab) call binary_mpi_output_array(luname, clean0r(array),&
          & mpi_dtype, comm)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        !ascii mpi-io is not possible, do it serially
        call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
           & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
           & root_is_in_comm, tag)
      case(binary_fmt)
        ! do truely parallel binary MPI-IO
        if(proc_is_in_slab) call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
           & comm)
      case default
        call gkw_throw_abort( &
     & "mpi_output_array_1d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      !parallel hdf5 io is not implemented yet, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
         & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
         & root_is_in_comm, tag)
    case(hdf5_mixed_fmt)
      !parallel hdf5 io is not implemented yet, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
         & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
         & root_is_in_comm, tag)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('mpi_output_array_3d_real')

  end subroutine mpi_output_array_3d_real

  !--------------------------------------------------------------------------
  !> should be called by *all* processes of COMM_CART.  (although at
  !> the moment just the procs working on the slab (plus or including
  !> root) are sufficient)
  !--------------------------------------------------------------------------
  subroutine mpi_output_array_6d_real(luname, groupname, array, mpi_dtype, &
     & buf_global, &
     & comm, fmt, preference_if_mixed, proc_is_in_slab, root_is_in_comm, tag)
    use io_binary, only : binary_mpi_output_array => mpi_output_array

    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: mpi_dtype
    real, dimension(:,:,:,:,:,:), intent(in) :: array
    real, dimension(:,:,:,:,:,:), intent(out) :: buf_global
    integer, intent(in) :: comm, tag

    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=*), intent(in) :: preference_if_mixed
    logical, intent(in) :: proc_is_in_slab, root_is_in_comm

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      !ascii mpi-io is not possible, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
         & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
         & root_is_in_comm, tag)
    case(binary_fmt)
      ! do truely parallel binary MPI-IO
      if(proc_is_in_slab) call binary_mpi_output_array(luname, clean0r(array),&
         & mpi_dtype, comm)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        !ascii mpi-io is not possible, do it serially
        call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
           & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
           & root_is_in_comm, tag)
      case(binary_fmt)
        ! do truely parallel binary MPI-IO
        if(proc_is_in_slab) call binary_mpi_output_array(luname, clean0r(array), mpi_dtype, &
           & comm)
      case default
        call gkw_throw_abort( &
           & "mpi_output_array_6d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      !parallel hdf5 io is not implemented yet, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
         & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
         & root_is_in_comm, tag)
    case(hdf5_mixed_fmt)
      !parallel hdf5 io is not implemented yet, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
         & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
         & root_is_in_comm, tag)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('mpi_output_array_6d_real')

  end subroutine mpi_output_array_6d_real


  !--------------------------------------------------------------------------
  !> should be called by *all* processes of COMM_CART.  (although at
  !> the moment just the procs working on the slab (plus or including
  !> root) are sufficient)
  !--------------------------------------------------------------------------
  subroutine mpi_output_array_3d_complex(luname, groupname, array, mpi_dtype, &
     & buf_global, &
     & comm, fmt, preference_if_mixed, proc_is_in_slab, root_is_in_comm, tag)
    use io_binary, only : binary_mpi_output_array => mpi_output_array
    
    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: mpi_dtype
    complex, dimension(:,:,:), intent(in) :: array
    complex, dimension(:,:,:), intent(out) :: buf_global
    integer, intent(in) :: comm, tag

    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=*), intent(in) :: preference_if_mixed
    logical, intent(in) :: proc_is_in_slab, root_is_in_comm

    select case(io_format)
    case(none_fmt)
      continue
    case(ascii_fmt)
      !ascii mpi-io is not possible, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
           & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
           & root_is_in_comm, tag)
    case(binary_fmt)
      ! do truely parallel binary MPI-IO
      if(proc_is_in_slab) call binary_mpi_output_array(luname, clean0c(array),&
          & mpi_dtype, comm)
    case(mixed_fmt)
      select case(preference_if_mixed)
      case(ascii_fmt)
        !ascii mpi-io is not possible, do it serially
        call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
           & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
           & root_is_in_comm, tag)
      case(binary_fmt)
        ! do truely parallel binary MPI-IO
        if(proc_is_in_slab) call binary_mpi_output_array(luname, clean0c(array), mpi_dtype, &
           & comm)
      case default
        call gkw_throw_abort( &
     & "mpi_output_array_1d_real: Unknown value for preference format: '"// &
           & preference_if_mixed//"'")
      end select
#if HAVE_HDF5
    case(hdf5_fmt)
      !parallel hdf5 io is not implemented yet, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
         & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
         & root_is_in_comm, tag)
    case(hdf5_mixed_fmt)
      !parallel hdf5 io is not implemented yet, do it serially
      call mpi_serial_output_array(luname, groupname, array, mpi_dtype, &
         & buf_global, comm, fmt, preference_if_mixed, proc_is_in_slab, &
         & root_is_in_comm, tag)
#endif
    case default
      call gkw_throw_abort('Unknown value for io_format: "'//io_format//'"')
    end select

    if(gkw_check_err()) call finalize_and_abort('mpi_output_array_3d_complex')

  end subroutine mpi_output_array_3d_complex

  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine mpi_serial_output_3d_cmplx(luname, groupname, array, &
     & mpi_dtype, buf_global, &
     & comm, fmt, preference_if_mixed, proc_is_in_slab, root_is_in_comm, tag)
    use mpiinterface, only : gather_array, root_processor

    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: mpi_dtype
    complex, dimension(:,:,:), intent(out) :: buf_global
    complex, dimension(:,:,:), intent(in) :: array
    integer, intent(in) :: comm, tag

    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=*), intent(in) :: preference_if_mixed
    logical, intent(in) :: proc_is_in_slab, root_is_in_comm
    
    logical, parameter :: to_root_of_commcart = .true.

    if (.false.) write(*,*) fmt

    if(root_processor .and. to_root_of_commcart .and. &
       & .not.proc_is_in_slab .and. root_is_in_comm) then
      ! if this is root, and we want to gather to root, and root is
      ! not in the slab and the flag to signal this is not set either...
      
      ! this is not good. Have coded wrong.
      write(*,*) root_processor, to_root_of_commcart, proc_is_in_slab, root_is_in_comm
      call gkw_throw_abort('mpi_serial_output_3d_cmplx: root is not in slab &
         & but root_is_in_comm is .true.')
      return
    end if

    if(proc_is_in_slab .or. (to_root_of_commcart .and. root_processor)) then
      ! gather the global array to root
      call gather_array(buf_global, &
         & array, &
         & mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
    end if

    ! serially output the global array from root
    if(root_processor) then
      call output_array(luname, groupname, &
         & buf_global, &
         & 'F', xy_fmt, preference_if_mixed)
    end if
  end subroutine mpi_serial_output_3d_cmplx

  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine mpi_serial_output_3d_real(luname, groupname, array, &
     & mpi_dtype, buf_global, &
     & comm, fmt, preference_if_mixed, proc_is_in_slab, root_is_in_comm, tag)
    use mpiinterface, only : gather_array, root_processor

    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: mpi_dtype
    real, dimension(:,:,:), intent(out) :: buf_global
    real, dimension(:,:,:), intent(in) :: array
    integer, intent(in) :: comm, tag
    logical, intent(in) :: proc_is_in_slab, root_is_in_comm

    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=*), intent(in) :: preference_if_mixed
    logical, parameter :: to_root_of_commcart = .true.

    if (.false.) write(*,*) fmt

    ! either this routine is called by a group of processes that
    ! happens to include root, or root calls this *additionally*
    
    if(root_processor .and. to_root_of_commcart .and. &
       & .not.proc_is_in_slab .and. root_is_in_comm) then
      ! if this is root, and we want to gather to root, and root is
      ! not in the slab and the flag to signal this is not set either...

      ! this is not good. Have coded wrong.
      write(*,*) root_processor, to_root_of_commcart, proc_is_in_slab, root_is_in_comm
      call gkw_throw_abort('mpi_serial_output_3d_real: root is not in slab')
      return
    end if

    ! gather the global array to root (with respect to COMM_CART)
    if(proc_is_in_slab .or. (to_root_of_commcart .and. root_processor)) then
      call gather_array(buf_global, &
         & array, &
         & mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
    end if

    ! serially output the global array from root
    if(root_processor) then
      call output_array(luname, groupname, &
         & buf_global, &
         & 'F', xy_fmt, preference_if_mixed)
    end if
  end subroutine mpi_serial_output_3d_real

  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine mpi_serial_output_6d_real(luname, groupname, array, &
     & mpi_dtype, buf_global, &
     & comm, fmt, preference_if_mixed, proc_is_in_slab, root_is_in_comm, tag)
    use mpiinterface, only : gather_array, root_processor

    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: mpi_dtype
    real, dimension(:,:,:,:,:,:), intent(out) :: buf_global
    real, dimension(:,:,:,:,:,:), intent(in) :: array
    integer, intent(in) :: comm, tag
    logical, intent(in) :: proc_is_in_slab, root_is_in_comm

    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=*), intent(in) :: preference_if_mixed
    logical, parameter :: to_root_of_commcart = .true.

    if (.false.) write(*,*) fmt

    ! either this routine is called by a group of processes that
    ! happens to include root, or root calls this *additionally*

    if(root_processor .and. to_root_of_commcart .and. &
       & .not.proc_is_in_slab .and. root_is_in_comm) then
      ! if this is root, and we want to gather to root, and root is
      ! not in the slab and the flag to signal this is not set either...

      ! this is not good. Have coded wrong.
      write(*,*) root_processor, to_root_of_commcart, proc_is_in_slab, root_is_in_comm
      call gkw_throw_abort('mpi_serial_output_6d_real: root is not in slab')
      return
    end if

    ! gather the global array to root (with respect to COMM_CART)
    if(proc_is_in_slab .or. (to_root_of_commcart .and. root_processor)) then
      call gather_array(buf_global, &
         & array, &
         & mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
    end if

    ! serially output the global array from root
    if(root_processor) then
      call output_array(luname, groupname, &
         & buf_global, &
         & 'F', xy_fmt, preference_if_mixed)
    end if
  end subroutine mpi_serial_output_6d_real

end module io
