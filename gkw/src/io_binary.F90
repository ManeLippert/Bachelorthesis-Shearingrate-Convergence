!------------------------------------------------------------------------------
!> Fortran binary output routines for the swappable interfaces in the IO module.
!------------------------------------------------------------------------------
module io_binary
  use global, only : gkw_throw_warn, root_and_verbose
  use io_ascii, only : write_lu_metadata
  
  implicit none

  private

  ! In terms of this module, a logical unit refers to one binary file,
  ! for both real and complex data.

  ! general
  public :: init, finalize
  public :: flush_all_open_files, flush_file
  public :: lu_exists
  ! FIXME this is only temporarily public, as long as we do not deal
  ! e.g. with 3DOutputParam.dat and the restart mechanism in a cleaner way
  ! (hence nobody should use it for diagnostic code):
  public :: get_free_file_unit

  ! reading and writing restart data
  ! FIXME this should be transferred here, thinks SRG

  ! subroutines to open, write (as often as desired) and close
  ! extendible datasets:
  public :: open_real_lu, open_complex_lu
  public :: close_lu, close_all_lus
  public :: append_chunk, seek_to_chunk
  public :: read_last_chunk

  ! subroutines which open, write (once and for all) and close fixed
  ! size datasets:
  public :: output_array
  public :: mpi_output_array

  !
  public :: mpifile_open
  public :: mpifile_close
  
  character(len=4), parameter :: suffix = ".bin"
  
  character(len=5), parameter :: real_suffix = "_real"
  character(len=5), parameter :: imag_suffix = "_imag"

  integer, parameter :: max_file_units = 256
  integer, save, dimension(max_file_units) :: file_unit_list
  logical, save, dimension(max_file_units) :: unit_is_complex
  
  !> number of files kept open (some of them may already be closed again)
  integer, save      :: n_open_file_units = 0

  !> gets a copy of mpiinterface value
  integer, save :: proc_number

  logical, save :: is_restarted
  logical, save :: io_legacy

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
    module procedure append_chunk_6d_real
    module procedure append_chunk_1d_complex
    module procedure append_chunk_2d_complex
    module procedure append_chunk_3d_complex
    module procedure append_chunk_4d_complex
    module procedure append_chunk_6d_complex
  end interface append_chunk

  interface read_last_chunk
    module procedure read_last_chunk_1d_real
  end interface read_last_chunk

  interface mpi_output_array
    module procedure mpi_output_array_1d_real
    module procedure mpi_output_array_3d_real
    module procedure mpi_output_array_6d_real
    module procedure mpi_output_array_3d_complex
  end interface mpi_output_array
  
contains

  !---------------------------------------------------------------------------
  !> A tiny helper function which is useful to make sure that only the
  !> root process does IO.
  !---------------------------------------------------------------------------
  function is_not_root_proc()
    logical :: is_not_root_proc
    is_not_root_proc = (proc_number/=0)
    if(is_not_root_proc) then
      ! call gkw_throw_warn("Parallel binary IO is not implemented yet.")
    end if
  end function is_not_root_proc
  

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine init(proc_number_, is_restarted_, io_legacy_)
    integer, intent(in) :: proc_number_
    logical, intent(in) :: is_restarted_, io_legacy_
    
    proc_number = proc_number_
    is_restarted = is_restarted_
    io_legacy = io_legacy_
  end subroutine init

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine finalize()
    ! There is nothing to finalize.
    continue
  end subroutine finalize

  !---------------------------------------------------------------------------
  !> flush all output files that were opened by this module
  !---------------------------------------------------------------------------
  subroutine flush_all_open_files()
    integer :: i
    do i = 1, n_open_file_units
      call flush_file(i)
    end do
  end subroutine flush_all_open_files


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine open_lu(luname, groupname, hyperslab_dims, lun)
    use global, only : int2char
    character (len=*), intent(in) :: luname, groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    integer, intent(out) :: lun
    
    integer :: iostatus

    call get_free_lun(lun)

    if (root_and_verbose) then
      write (*,*) '* opening file '//luname//', LOGICAL UNIT='//int2char(lun)//','
      write (*,*) '* FILE UNIT='//int2char(file_unit_list(lun))//','
      write (*,*) '  FORM=unformatted, POSITION=append'
    end if

    if(lu_exists(luname, groupname) .and. is_restarted) then
      ! Open an existing file and put the pointer at the end.
      if (root_and_verbose) write(*,*) luname, ' is found. Will append to it.'
      open(file_unit_list(lun), FILE=luname, FORM='unformatted', &
         & STATUS='old', POSITION='append', IOSTAT=iostatus)
    else
      ! Overwrite an existing file or create a new one.
      if (root_and_verbose) write(*,*) luname, &
         &' not found or this is not a restarted run. Will create new.'
      
      ! FIXME Should one specify recl? Do we want ACCSS=direct?
      !open(lun,FILE=luname,FORM='unformatted',ACCESS='direct',recl=??,...)
      open(file_unit_list(lun), FILE=luname, FORM='unformatted', &
         & STATUS='replace', POSITION='rewind', IOSTAT=iostatus)
    end if

    if (.false.) write(*,*) hyperslab_dims

  end subroutine open_lu
    

  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine open_real_lu(luname, groupname, hyperslab_dims, lun)
    use global, only : ss
    character (len=*), intent(in) :: luname, groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    integer, intent(out) :: lun

    if(is_not_root_proc()) return

    if(io_legacy) then
      call open_lu(luname, groupname, hyperslab_dims, lun)
      call write_lu_metadata(luname, groupname, 'real', hyperslab_dims, .true.)
    else
      call open_lu(trim(ss(luname))//suffix, groupname, hyperslab_dims, lun)
      call write_lu_metadata(trim(ss(luname))//suffix, groupname, 'real', &
         & hyperslab_dims, .true.)
    end if

  end subroutine open_real_lu

  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine open_complex_lu(luname, groupname, hyperslab_dims, real_lun)
    use global, only : ss
    character (len=*), intent(in) :: luname, groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    ! the returned lun will be the one corresponding to the real part
    integer, intent(out) :: real_lun
    integer :: imag_lun
    if(is_not_root_proc()) return

    if(io_legacy) then
      call open_lu(luname, groupname, hyperslab_dims, real_lun)
      call write_lu_metadata(luname, groupname, 'complex', hyperslab_dims, &
         & .true.)
    else
      ! real and imag part into one file:
      ! call open_lu(trim(ss(luname))//suffix, groupname, hyperslab_dims, real_lun)
      ! call write_lu_metadata(trim(ss(luname))//suffix, groupname, 'complex', &
      !    & hyperslab_dims, .true.)
      ! ... or into two separate files:
      call open_lu(trim(ss(luname))//real_suffix//suffix, groupname, &
        & hyperslab_dims, real_lun)
      call open_lu(trim(ss(luname))//imag_suffix//suffix, groupname, &
        & hyperslab_dims, imag_lun)
      call write_lu_metadata(trim(ss(luname))//real_suffix//suffix// &
         & ' and '//trim(ss(luname))//imag_suffix//suffix, &
         & groupname, 'complex', hyperslab_dims, .true.)

      unit_is_complex(imag_lun) = .false.
    end if
    ! set an entry which can be tested when closing the dataset.
    unit_is_complex(real_lun) = .true.


  end subroutine open_complex_lu

  !----------------------------------------------------------------------------
  !> Close the file associated to the given logical unit number.
  !----------------------------------------------------------------------------
  subroutine close_lu(lun)

    integer, intent(in) :: lun

    if(is_not_root_proc()) return

    if(unit_is_complex(lun)) then
      call close_complex_lu(lun)
    else
      call close_real_lu(lun)
    end if

  end subroutine close_lu

  !---------------------------------------------------------------------------
  !> This subroutine closes a complex logical unit.
  !---------------------------------------------------------------------------
  subroutine close_real_lu(lun)
    use global, only : int2char
    integer, intent(in) :: lun
    logical :: is_open

    ! unit numbers may not be negative
    ! If negative, means the file was never opened
    ! The inquire statement may fail with negative lun
    ! For the IBM compiler, unit zero is always open but can never be closed
    
    ! check for invalid argument:
    if (lun > n_open_file_units .or. lun <= 0) then
      call gkw_throw_warn("Invalid argument to ascii close_real_lu() ")
      return
    end if

    inquire(unit=file_unit_list(lun), opened=is_open)
    if(is_open) then
      close(file_unit_list(lun))

      if (root_and_verbose) write(*,*) '* closed file unit '//int2char(lun)
    else
      if (root_and_verbose) call gkw_throw_warn("Attempted to close file unit "// &
         & int2char(file_unit_list(lun))// &
         & " associated to logical unit "//int2char(lun)//" more than once.")
    end if

  end subroutine close_real_lu

  !---------------------------------------------------------------------------
  !> This subroutine closes a complex logical unit.
  !---------------------------------------------------------------------------
  subroutine close_complex_lu(lun)
    integer, intent(in) :: lun

    if(is_not_root_proc()) return

    ! End access to the dataset with the real part.
    call close_real_lu(lun)
    if(io_legacy) then
      ! Close also the associated dataset with the imaginary part.
      call close_real_lu(lun+1)
    end if
  end subroutine close_complex_lu

  !----------------------------------------------------------------------------
  !> Close all open files.
  !----------------------------------------------------------------------------
  subroutine close_all_lus
    use global, only : int2char

    integer :: i
    logical :: is_open

    if(is_not_root_proc()) return

    do i = 1, n_open_file_units
      inquire(unit=file_unit_list(i),opened=is_open)
      if(is_open) then
        close(file_unit_list(i))
        if(root_and_verbose) write(*,*) '* closed file unit '// &
           & int2char(file_unit_list(i))
      end if
    end do
  end subroutine close_all_lus

  !---------------------------------------------------------------------------
  !> Check wether a file with the given name exists.
  !---------------------------------------------------------------------------
  function lu_exists(luname, groupname)
    use global, only : ss
    character (len=*), intent(in) :: luname, groupname
    logical :: lu_exists

    if(is_not_root_proc()) then
      lu_exists = .false.
    else
      if(io_legacy) then
        inquire(file=trim(luname), exist = lu_exists)
      else
        inquire(file=trim(ss(luname))//suffix, exist = lu_exists)
      end if
    end if

    if (.false.) write(*,*) groupname

  end function lu_exists


  !----------------------------------------------------------------------------
  !> Obtain a logical unit number which is associated to a  free file unit
  !> for input/output.
  !> For safety, do not reuse the same file units on different processors.
  !> (Note that here a logical unit number is not the same as a file unit.)
  !>
  !> A local variable proc_number is used to avoid a circular dependency.
  !>
  !----------------------------------------------------------------------------
  subroutine get_free_lun(lun)
    use global, only : gkw_throw_abort
    integer, intent(out) :: lun
    integer :: file_unit

    if(is_not_root_proc()) return

    if(n_open_file_units < max_file_units - 1) then

      call get_free_file_unit(file_unit)
      
      n_open_file_units = n_open_file_units + 1
      
      ! Register this file_unit in the list of open file units.
      lun = n_open_file_units
      file_unit_list(lun) = file_unit
      return
    else
      call gkw_throw_abort('binary IO: Increase max_file_units, or use less diagnostics')
      return
    end if
    
  end subroutine get_free_lun

  subroutine get_free_file_unit(file_unit)
    use global, only : gkw_throw_abort
    integer, intent(out) :: file_unit
    logical :: is_open

    ! In restart.F90 this routine is called by 'the last process'.
    !! if(is_not_root_proc()) return
    
    
    ! every processor has a number of max_file_units logical units at
    ! disposition.
    do file_unit = (110 + max_file_units*proc_number), &
       & (110 + max_file_units*(proc_number+1))
      inquire(unit=file_unit, opened=is_open)
      if (.not. is_open) then
        ! A lun was found which is not yet associated to an open file.
        return
      end if
    end do

    file_unit = -1
    
    ! Program flow reaches this point if
    !  a) no free file unit could be found.
    !  b) max_file_units was too small.

    write(*,*) 'binary IO: Increase max_file_units, or use less diagnostics'
    call gkw_throw_abort('binary IO: Increase max_file_units, or use less diagnostics')

  end subroutine get_free_file_unit


  !----------------------------------------------------------------------------
  !> This function attempts to imitate the flush intrinsic in Fortran 2003. 
  !> It if causes difficulties for a given compiler, the contents can simply be
  !> commented and the function made empty; flushing files to disk is only for 
  !> convenience and the code will run fine without it.
  !----------------------------------------------------------------------------
  subroutine flush_file(lun)
    use global, only : i8 

    integer, intent(in) :: lun
    integer, parameter :: max_filename_len = 512
    character (len=max_filename_len) :: fname
    ! The integer below is kind=i8 for the IBM compiler (see r3446, Issue 203)
    ! However, the Intel compiler complains that this violates the Fortran95
    ! standard. Unless an alterative completely portable solution is found 
    ! it is better to break the Intel interpretation of the standard than 
    ! to have something that won't run with a known, in-use compiler.
    ! You can remove the kind specification in your local copy if required.
    integer (kind=i8) :: recll

    if(is_not_root_proc()) return

    if (file_is_open(lun)) then
#ifdef STD2003_FLUSH
      flush(file_unit_list(lun))
#else
      inquire(UNIT=file_unit_list(lun),NAME=fname)
      inquire(UNIT=file_unit_list(lun),RECL=recll)
      close(file_unit_list(lun))

      ! workaround for gfortran bug:
      ! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53796
      if (recll < 1) then
        open(file_unit_list(lun), FILE=fname, &
           & FORM='unformatted', POSITION='append')
      else
        open(file_unit_list(lun), FILE=fname, &
           & FORM='unformatted', POSITION='append', RECL=recll)
      end if
#endif
    end if

  end subroutine flush_file

  !----------------------------------------------------------------------------
  !> 
  !----------------------------------------------------------------------------
  function file_is_open(lun)
    integer, intent(in) :: lun

    logical :: file_is_open

    if(is_not_root_proc()) then
      file_is_open = .false.
      return
    end if

    if (lun <= 0 .or. lun > n_open_file_units) then
      file_is_open = .false.
    else
      inquire(unit=file_unit_list(lun),opened=file_is_open)
    end if
  end function file_is_open


  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine open_file_to_append_or_replace(filename, groupname, file_unit)
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    integer, intent(in) :: file_unit
    integer :: ios

    character(len=128) :: full_filename

    if(io_legacy) then
      full_filename = filename
    else
      full_filename = trim(filename)//suffix
    end if
    ! Overwrite the file called full_filename
    
    ! FIXME The only difference between output_slice and output_array
    ! seems to be the presence/absence of the RECL argument. Should
    ! one specify it?
    
    if (root_and_verbose) write(*,*) full_filename, &
       &' will be overwritten, if it exists.'
    
    open(file_unit, FILE=full_filename, FORM='unformatted', &
       & STATUS='replace', POSITION='rewind', IOSTAT=ios)
    if(ios /= 0) then
      write(*,*) 'Error while opening ',full_filename,', iostat=',ios
      write(*,*) 'file_unit=', file_unit
    end if

    if (.false.) write(*,*) groupname

  end subroutine open_file_to_append_or_replace

  !----------------------------------------------------------------------------
  !> Output a real 1D array data(:) to file, as unformatted raw binary.
  !> The order 'C' is not implemented for binary output (no demand yet).
  !----------------------------------------------------------------------------
  subroutine output_array_1d(filename, groupname, data,order)
    real, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    integer :: file_unit
    !integer :: record_length

    if(order == 'C' .and. root_and_verbose) then
      call gkw_throw_warn("The order 'C' is not implemented for binary output.")
    end if

    call get_free_file_unit(file_unit)

    ! FIXME The only difference between output_slice and output_array
    ! seems to be the presence/absence of the RECL argument
    !inquire(iolength=record_length) data
    call open_file_to_append_or_replace(filename, groupname, file_unit)
    write(file_unit) data
    close(file_unit)

  end subroutine output_array_1d

  !----------------------------------------------------------------------------
  !> Output a real 1D array data(:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_1d_real(filename, groupname, data,order)
    real, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_1d(filename, groupname, data,order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)

  end subroutine output_array_1d_real

  !----------------------------------------------------------------------------
  !> Output a real 2D array data(:,:) to file, as unformatted raw binary.
  !> The order 'C' is not implemented for binary output (no demand yet).
  !----------------------------------------------------------------------------
  subroutine output_array_2d(filename, groupname, data,order)
    real, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    integer :: file_unit
    !integer :: record_length

    if(order == 'C' .and. root_and_verbose) then
      call gkw_throw_warn("The order 'C' is not implemented for binary output.")
    end if
    
    call get_free_file_unit(file_unit)

    ! FIXME The only difference between output_slice and output_array
    ! seems to be the presence/absence of the RECL argument
    !inquire(iolength=record_length) data
    call open_file_to_append_or_replace(filename, groupname, file_unit)
    write(file_unit) data
    close(file_unit)

  end subroutine output_array_2d

  !----------------------------------------------------------------------------
  !> Output a real 2D array data(:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_2d_real(filename, groupname, data,order)
    real, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_2d(filename, groupname, data,order)

    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)

  end subroutine output_array_2d_real

  !----------------------------------------------------------------------------
  !> Output a real 3D array data(:,:,:) to file, as unformatted raw binary.
  !> The order 'C' is not implemented for binary output (no demand yet).
  !----------------------------------------------------------------------------
  subroutine output_array_3d(filename, groupname, data,order)
    real, dimension(:, :,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    integer :: file_unit
    !integer :: record_length

    if(order == 'C' .and. root_and_verbose) then
      call gkw_throw_warn("The order 'C' is not implemented for binary output.")
    end if

    call get_free_file_unit(file_unit)

    ! FIXME The only difference between output_slice and output_array
    ! seems to be the presence/absence of the RECL argument
    !inquire(iolength=record_length) data
    call open_file_to_append_or_replace(filename, groupname, file_unit)
    write(file_unit) data
    close(file_unit)

  end subroutine output_array_3d

  !----------------------------------------------------------------------------
  !> Output a real 3D array data(:,:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_3d_real(filename, groupname, data, order)
    real, dimension(:, :,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_3d(filename, groupname, data, order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)

  end subroutine output_array_3d_real

  !----------------------------------------------------------------------------
  !> Output a real 4D array data to file, as unformatted raw binary.
  !> The order 'C' is not implemented for binary output (no demand yet).
  !----------------------------------------------------------------------------
  subroutine output_array_4d(filename, groupname, data,order)
    real, dimension(:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    integer :: file_unit
    !integer :: record_length

    if(order == 'C'.and. root_and_verbose) then
      call gkw_throw_warn("The order 'C' is not implemented for binary output.")
    end if

    call get_free_file_unit(file_unit)

    ! FIXME The only difference between output_slice and output_array
    ! seems to be the presence/absence of the RECL argument
    !inquire(iolength=record_length) data
    call open_file_to_append_or_replace(filename, groupname, file_unit)
    write(file_unit) data
    close(file_unit)

  end subroutine output_array_4d

  !----------------------------------------------------------------------------
  !> Output a real 4D array data(:,:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_4d_real(filename, groupname, data, order)
    real, dimension(:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_4d(filename, groupname, data, order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)

  end subroutine output_array_4d_real

  !----------------------------------------------------------------------------
  !> Output a real 5D array data to file, as unformatted raw binary.
  !> The order 'C' is not implemented for binary output (no demand yet).
  !----------------------------------------------------------------------------
  subroutine output_array_5d(filename, groupname, data,order)
    real, dimension(:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    integer :: file_unit
    !integer :: record_length

    if(order == 'C'.and. root_and_verbose) then
      call gkw_throw_warn("The order 'C' is not implemented for binary output.")
    end if

    call get_free_file_unit(file_unit)

    ! FIXME The only difference between output_slice and output_array
    ! seems to be the presence/absence of the RECL argument
    !inquire(iolength=record_length) data
    call open_file_to_append_or_replace(filename, groupname, file_unit)
    write(file_unit) data
    close(file_unit)

  end subroutine output_array_5d

  !----------------------------------------------------------------------------
  !> Output a real 5D array data(:,:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_5d_real(filename, groupname, data, order)
    real, dimension(:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_5d(filename, groupname, data, order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)

  end subroutine output_array_5d_real

  !----------------------------------------------------------------------------
  !> Output a real 6D array data to file, as unformatted raw binary.
  !> The order 'C' is not implemented for binary output (no demand yet).
  !----------------------------------------------------------------------------
  subroutine output_array_6d(filename, groupname, data,order)
    real, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    integer :: file_unit
    !integer :: record_length

    if(order == 'C'.and. root_and_verbose) then
      call gkw_throw_warn("The order 'C' is not implemented for binary output.")
    end if

    call get_free_file_unit(file_unit)

    ! FIXME The only difference between output_slice and output_array
    ! seems to be the presence/absence of the RECL argument
    !inquire(iolength=record_length) data
    call open_file_to_append_or_replace(filename, groupname, file_unit)
    write(file_unit) data
    close(file_unit)

  end subroutine output_array_6d
  
  !----------------------------------------------------------------------------
  !> Output a real 6D array data(:,:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_6d_real(filename, groupname, data, order)
    real, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_6d(filename, groupname, data, order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)

  end subroutine output_array_6d_real

  !----------------------------------------------------------------------------
  !> Output a complex 1D array data(:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_1d_complex(filename, groupname, data,order)
    complex, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_1d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order)
    call output_array_1d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(trim(filename)//real_suffix//suffix//' and ' &
       & //trim(filename)//imag_suffix//suffix, &
       & groupname, 'complex', shape(data), .false.)

  end subroutine output_array_1d_complex

  !----------------------------------------------------------------------------
  !> Output a complex 2D array data(:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_2d_complex(filename, groupname, data,order)
    complex, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_2d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order)
    call output_array_2d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(trim(filename)//real_suffix//suffix//' and ' &
       & //trim(filename)//imag_suffix//suffix, &
       & groupname, 'complex', shape(data), .false.)

  end subroutine output_array_2d_complex

  !----------------------------------------------------------------------------
  !> Output a complex 3D array data(:,:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_3d_complex(filename, groupname, data,order)
    complex, dimension(:, :,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_3d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order)
    call output_array_3d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(trim(filename)//real_suffix//suffix//' and ' &
       & //trim(filename)//imag_suffix//suffix, &
       & groupname, 'complex', shape(data), .false.)

  end subroutine output_array_3d_complex

  !----------------------------------------------------------------------------
  !> Output a complex 6D array data(:,:,:,:,:,:) to file, as unformatted raw binary.
  !----------------------------------------------------------------------------
  subroutine output_array_6d_complex(filename, groupname, data,order)
    complex, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename, groupname
    character (len=1), intent(in) :: order

    call output_array_6d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order)
    call output_array_6d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(trim(filename)//real_suffix//suffix//' and ' &
       & //trim(filename)//imag_suffix//suffix, &
       & groupname, 'complex', shape(data), .false.)

  end subroutine output_array_6d_complex

  !----------------------------------------------------------------------------
  !> This subroutine appends a 1D vector data(:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_1d_real(lun, data)
    real, dimension(:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! FIXME which one is better?
    !write(file_unit_list(lun), RECL=1) data
    write(file_unit_list(lun)) data

  end subroutine append_chunk_1d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 2D matrix data(:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_2d_real(lun, data)
    real, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! FIXME which one is better?
    !write(file_unit_list(lun), RECL=1) data
    write(file_unit_list(lun)) data

  end subroutine append_chunk_2d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 3D matrix data(:,:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_3d_real(lun, data)
    real, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! FIXME which one is better?
    !write(file_unit_list(lun), RECL=1) data
    write(file_unit_list(lun)) data

  end subroutine append_chunk_3d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 4D matrix data(:,:,:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_4d_real(lun, data)
    real, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! FIXME which one is better?
    !write(file_unit_list(lun), RECL=1) data
    write(file_unit_list(lun)) data

  end subroutine append_chunk_4d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 6D matrix to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_6d_real(lun, data)
    real, dimension(:,:,:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! FIXME which one is better?
    !write(file_unit_list(lun), RECL=1) data
    write(file_unit_list(lun)) data

  end subroutine append_chunk_6d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 1D vector data(:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_1d_complex(lun, data)
    complex, dimension(:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    if(unit_is_complex(lun)) then
      if(io_legacy) then
        write(file_unit_list(lun)) data
      else
        call append_chunk_1d_real(lun, real(data))
        call append_chunk_1d_real(lun+1, aimag(data))
      end if
    else
      call gkw_throw_warn("Complex data should not be appended to a real &
         & valued logical unit.")
    end if

  end subroutine append_chunk_1d_complex

  !----------------------------------------------------------------------------
  !> This subroutine appends a 2D matrix data(:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_2d_complex(lun, data)
    complex, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      if(io_legacy) then
        write(file_unit_list(lun)) data
      else
        call append_chunk_2d_real(lun, real(data))
        call append_chunk_2d_real(lun+1, aimag(data))
      end if
    else
      call gkw_throw_warn("Complex data should not be appended to a real &
         & valued logical unit.")
    end if

  end subroutine append_chunk_2d_complex

  !----------------------------------------------------------------------------
  !> This subroutine appends a 3D matrix data(:,:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_3d_complex(lun, data)
    complex, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      if(io_legacy) then
        write(file_unit_list(lun)) data
      else
        call append_chunk_3d_real(lun, real(data))
        call append_chunk_3d_real(lun+1, aimag(data))
      end if
    else
      call gkw_throw_warn("Complex data should not be appended to a real &
         & valued logical unit.")
    end if

  end subroutine append_chunk_3d_complex

  !----------------------------------------------------------------------------
  !> This subroutine appends a 4D matrix data(:,:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_4d_complex(lun, data)
    complex, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      if(io_legacy) then
        write(file_unit_list(lun)) data
      else
        call append_chunk_4d_real(lun, real(data))
        call append_chunk_4d_real(lun+1, aimag(data))
      end if
    else
      call gkw_throw_warn("Complex data should not be appended to a real &
         & valued logical unit.")
    end if

  end subroutine append_chunk_4d_complex

  !----------------------------------------------------------------------------
  !> This subroutine appends a 6D matrix to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_6d_complex(lun, data)
    complex, dimension(:,:,:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    if(is_not_root_proc()) return
    
    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      if(io_legacy) then
        write(file_unit_list(lun)) data
      else
        call append_chunk_6d_real(lun, real(data))
        call append_chunk_6d_real(lun+1, aimag(data))
      end if
    else
      call gkw_throw_warn("Complex data should not be appended to a real &
         & valued logical unit.")
    end if

  end subroutine append_chunk_6d_complex
  
  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine seek_to_chunk(lun, index)
    use global, only : gkw_throw_abort
    integer, intent(in) :: lun
    integer, intent(in) ::  index
    call gkw_throw_abort("The routine seek_to_chunk was not implemented for binary.")
    return
    
    ! no demand yet
    if (.false.) write(*,*) lun, index
    
  end subroutine seek_to_chunk

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine read_last_chunk_1d_real(lun, data, nlast)
    use global, only : gkw_throw_abort
    integer, intent(in) :: lun
    real, dimension(:), intent(out) :: data
    integer, intent(out) :: nlast
    
    if(is_not_root_proc()) return
    
    nlast = 0
    data = 0

    call gkw_throw_abort("The routine read_last_chunk was not implemented for binary.")
    return
    
    ! no demand yet
    if (.false.) write(*,*) lun
    
  end subroutine read_last_chunk_1d_real

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine mpi_output_array_1d_real(luname, array, view_type, &
       & comm)

    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: view_type
    real, dimension(:), intent(in) :: array
    integer, intent(in), optional :: comm

#if defined(mpi2)
    if(present(comm)) then
      call mpi_output_array_generic(luname, array, size(array), view_type, &
         & comm)
    else
      call mpi_output_array_generic(luname, array, size(array), view_type &
         & )
    end if
#else
    call output_array_1d_real(luname, 'mpi_tmp', array, 'F')
#endif
  end subroutine mpi_output_array_1d_real


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine mpi_output_array_3d_real(luname, array, view_type, &
       & comm)
    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: view_type
    real, dimension(:,:,:), intent(in) :: array
    integer, intent(in), optional :: comm

#if defined(mpi2)
    if(present(comm)) then
      call mpi_output_array_generic(luname, array, size(array), view_type, &
         & comm)
    else
      call mpi_output_array_generic(luname, array, size(array), view_type &
         & )
    end if
#else
    call output_array_3d_real(luname, 'mpi_tmp', array, 'F')
#endif
  end subroutine mpi_output_array_3d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine mpi_output_array_6d_real(luname, array, view_type, &
       & comm)
    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: view_type
    real, dimension(:,:,:,:,:,:), intent(in) :: array
    integer, intent(in), optional :: comm

#if defined(mpi2)
    if(present(comm)) then
      call mpi_output_array_generic(luname, array, size(array), view_type, &
         & comm)
    else
      call mpi_output_array_generic(luname, array, size(array), view_type &
         & )
    end if
#else
    call output_array_6d_real(luname, 'mpi_tmp', array, 'F')
#endif
  end subroutine mpi_output_array_6d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine mpi_output_array_3d_complex(luname, array, view_type, &
       & comm)
    character(len=*), intent(in)  :: luname
    integer, intent(in)           :: view_type
    complex, dimension(:,:,:), intent(in) :: array
    integer, intent(in), optional :: comm

#if defined(mpi2)
    if(present(comm)) then
      call mpi_output_array_generic(trim(luname)//real_suffix, real(array), &
         & size(array), view_type, comm)
      call mpi_output_array_generic(trim(luname)//imag_suffix, aimag(array), &
         & size(array), view_type, comm)
    else
      call mpi_output_array_generic(trim(luname)//real_suffix, real(array), &
         & size(array), view_type)
      call mpi_output_array_generic(trim(luname)//imag_suffix, aimag(array), &
         & size(array), view_type)
    end if
#else
    call output_array_3d_complex(luname, 'mpi_tmp', array, 'F')
#endif
  end subroutine mpi_output_array_3d_complex

  !---------------------------------------------------------------------------
  !> This routine opens a file with MPI-IO, in write mode.
  !---------------------------------------------------------------------------
  subroutine mpifile_open(filename,file_unit,comm,append_to_file)
    use mpiinterface, only : MPI_INFO_NULL, MPI_MODE_WRONLY, MPI_MODE_CREATE
    use mpiinterface, only : MPI_MODE_APPEND
    use mpiinterface, only : processor_number, mpiabort
    use global, only : gkw_throw_abort
    character(len=*), intent(in)  :: filename
    integer, intent(out)          :: file_unit
    integer, intent(in) :: comm
    logical, intent(in) :: append_to_file

#if defined(mpi2)
    integer :: ierr,amode

    amode = MPI_MODE_WRONLY + MPI_MODE_CREATE
    if(append_to_file) then
      amode = amode + MPI_MODE_APPEND
    end if
    call mpi_file_open(comm,filename,amode, &
       & MPI_INFO_NULL,file_unit,ierr)
    if (ierr /= 0) then
      write(*,*) 'Error calling MPI_FILE_OPEN', processor_number
      !FIXME choose one - how should one deal with an error here?
      call gkw_throw_abort('Error calling mpi_file_open')
      call mpiabort
    end if
#else
    write(*,*) 'Cannot do mpi IO without mpi!'
    stop 1
#endif

  end subroutine mpifile_open

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine mpifile_close(file_unit)
    use mpiinterface, only : processor_number, mpiabort
    use global, only : gkw_throw_abort
    integer, intent(in) :: file_unit
    integer :: ierr
    
#if defined(mpi2)
    call mpi_file_close(file_unit,ierr)
    if (ierr /= 0) then
      write(*,*) 'Error calling MPI_FILE_CLOSE', processor_number
      !FIXME choose one - how should one deal with an error here?
      call gkw_throw_abort('Error calling mpi_file_open')
      call mpiabort
    end if
#else
    write(*,*) 'Cannot do mpi IO without mpi!'
    stop 1
#endif

  end subroutine mpifile_close

  !---------------------------------------------------------------------------
  !> this routine must be called by all processes of the given
  !> communicator comm. The given local arrays with 'nelem' elements are
  !> written at the position specified by view_type and view_disp.
  !> 
  !---------------------------------------------------------------------------
  subroutine mpi_output_array_generic(filename, array, nelem, view_type, &
     & comm)
    use mpiinterface, only : MPI_OFFSET_KIND
    use mpiinterface, only : MPI_INFO_NULL
    use mpiinterface, only : MPIREAL_X, statusmpi
    use mpicomms, only : COMM_CART
    use global, only : gkw_throw_warn
    character(len=*), intent(in)  :: filename
    integer, intent(in)           :: nelem, view_type
    real, intent(in)              :: array(nelem)
    integer, intent(in), optional :: comm
#if defined(mpi2)
    integer :: file_unit
    integer :: communicator, ierr 
    integer(kind=MPI_OFFSET_KIND), parameter :: view_disp = 0

    communicator = COMM_CART
    if (present(comm)) communicator = comm

    call mpifile_open(filename, file_unit, communicator, .false.)
    call mpi_file_set_view(file_unit, view_disp, MPIREAL_X, view_type, &
       & "native", MPI_INFO_NULL, ierr)
    call mpi_file_write_all(file_unit, array, nelem, MPIREAL_X, &
       & statusmpi, ierr)
    call mpi_file_close(file_unit,ierr)

#else
    ! This should never ever be the case.
    call gkw_throw_warn("No MPI, but mpi_output_array_generic() called!")
#endif

  end subroutine mpi_output_array_generic

end module io_binary
