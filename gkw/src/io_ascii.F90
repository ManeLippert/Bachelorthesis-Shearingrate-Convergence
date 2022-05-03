!------------------------------------------------------------------------------
!> ASCII output routines for the swappable interfaces of the IO module.
!------------------------------------------------------------------------------
module io_ascii
  use global, only : gkw_throw_warn, root_and_verbose
  
  implicit none

  private

  ! In terms of this module, a logical unit refers to one ascii file,
  ! for both real and complex data.

  ! general
  public :: init, finalize
  public :: flush_all_open_files, flush_file
  public :: lu_exists
  ! this should not be used for diagnostic code (!) :
  public :: get_free_file_unit
  
  ! subroutines to open, write (as often as desired) and close
  ! extendible datasets:
  public :: open_real_lu, open_complex_lu
  public :: close_lu, close_all_lus
  public :: append_chunk, seek_to_chunk
  public :: read_last_chunk

  public :: write_lu_metadata
  public :: attach_metadata

  ! subroutines which open, write (once and for all) and close fixed
  ! size datasets:
  public :: output_array
  
  ! other
  public :: xy_fmt, xy_fmt_long

  ! Historically, there were 2 different kinds of data files:
  !  1.) Files which contained a single n-dimensional dataset, like time.dat,
  !      parallel.dat and entr.kyspec
  !  2.) Files which contained a multitude of datasets, of different
  !      dimensionality respectively, like geom.dat
  !
  ! It is difficult for this module to reproduce both kinds as before,
  ! without a too complicated interface.
  ! SRG thinks sooner or later we should deprecate the 2nd kind of files
  ! and get rid of them by writing several files of the 1st kind.
  !
  ! In the meanwhile, this module will
  !    * produce files of the 2nd kind only if the given
  !      groupname appears in the list file_2nd_kind_groupnames
  !    * produce files of the 1st kind otherwise.
  ! Other parts of the code will not have any control over which kind of
  ! file is produced.
  !
  ! One simplification:
  ! Historically, data was not appended piecewise to files of the 2nd kind,
  ! but "at a single blow". Therefore the special handling of the groupname
  ! is not implemented into the open_*_lu, close_lu and append_chunk_*
  ! routines, but only into the output_array_* routines.

  integer, parameter :: luname_maxsize = 256
  character(len=4), parameter :: suffix = ".dat"
  character(len=12), parameter :: meta_filename = "gkwdata.meta"
  integer, parameter :: meta_file_unit = 2
  
  integer, parameter :: groupname_maxsize = 40
  character(len=groupname_maxsize), dimension(1) :: file_2nd_kind_groupnames
  logical, dimension(size(file_2nd_kind_groupnames)) :: &
     & file_2nd_kind_is_touched = .false.
  
  integer, parameter :: max_file_units = 512
  integer, save, dimension(max_file_units) :: file_unit_list
  logical, save, dimension(max_file_units) :: unit_is_complex

  !> number of files kept open (some of them may already be closed again)
  integer, save      :: n_open_file_units = 0

  character(len=5), parameter :: real_suffix = "_real"
  character(len=5), parameter :: imag_suffix = "_imag"

  !> max. number of columns in ascii output.
  integer, parameter :: max_cols = 4096
  !> size of each output number (according to xy_fmt) in characters
  !> (i.e. bytes) in the file (sign, significand,'E',sign,exponent,blank)
  integer, parameter :: fmt_wid = 1+18+1+2+1
  character (len=18), parameter :: xy_fmt = '(65536(es13.5,1x))'
  ! an alternative format to xy_fmt, to print some more decimal places
  ! for debugging purposes:
  character (len=18), parameter :: xy_fmt_long = '(65536(es18.8,1x))'

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
  
contains

  !---------------------------------------------------------------------------
  !> A tiny helper function which is useful to make sure that only the
  !> root process does IO.
  !---------------------------------------------------------------------------
  function is_not_root_proc()
    logical :: is_not_root_proc
    is_not_root_proc = (proc_number/=0)
    if(is_not_root_proc) then
      ! call gkw_throw_warn("Parallel ascii IO is not implemented yet.")
    end if
  end function is_not_root_proc
  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine init(proc_number_, is_restarted_, create_metadata_file, io_legacy_)
    integer, intent(in) :: proc_number_
    logical, intent(in) :: is_restarted_
    logical, intent(in) :: create_metadata_file, io_legacy_
    integer :: iostatus
    proc_number = proc_number_
    is_restarted = is_restarted_
    io_legacy = io_legacy_

    file_2nd_kind_groupnames(1) = 'geom'

    if(is_not_root_proc()) return

    if(create_metadata_file) then
      ! open the metadata file, just to discard previous content.
      open(meta_file_unit, FILE=meta_filename, FORM='formatted', &
         & STATUS='replace', POSITION='rewind', IOSTAT=iostatus, &
         & ACTION='write')
      close(meta_file_unit)
    end if

  end subroutine init


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine finalize()
    
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
  subroutine open_metadata_file()
    integer :: iostatus
    ! instead of requesting a free file_unit each time with
    !call get_free_file_unit(meta_file_unit)
    ! we take a hardcoded file_unit because the metadata file needs to be
    ! frequently opened and closed by different processes.
    open(meta_file_unit, FILE=meta_filename, FORM='formatted', &
       & STATUS='old', POSITION='append', IOSTAT=iostatus, &
       & ACTION='write')
    
  end subroutine open_metadata_file
  
  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine open_real_lu(luname, groupname, hyperslab_dims, lun)
    use global, only : ss

    character (len=*), intent(in) :: luname
    ! the ascii IO ignores the groupname here
    character (len=*), intent(in) :: groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    integer, intent(out) :: lun

    if(is_not_root_proc()) return

    if(io_legacy) then
      call open_lu(trim(luname), groupname, hyperslab_dims, lun)
      call write_lu_metadata(luname, groupname, 'real', hyperslab_dims, .true.)
    else
      call open_lu(trim(ss(luname))//suffix, groupname, hyperslab_dims, lun)
      call write_lu_metadata(trim(ss(luname))//suffix, groupname, 'real', hyperslab_dims, .true.)
    end if
    
  end subroutine open_real_lu

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine open_lu(luname, groupname, hyperslab_dims, lun)
    use global, only : int2char, gkw_its_critical

    character (len=*), intent(in) :: luname
    ! the ascii IO ignores the groupname here
    character (len=*), intent(in) :: groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    integer, intent(out) :: lun
    integer :: iostatus

    integer :: record_length

    record_length = hyperslab_dims(1) * fmt_wid
    
    call get_free_lun(lun)
    if(gkw_its_critical()) return

    if (root_and_verbose) then
      write (*,*) '* opening file '//luname//', LOGICAL UNIT='//int2char(lun)//','
      write (*,*) '* FILE UNIT='//int2char(file_unit_list(lun))//','
      write (*,*) '  FORM=formatted, POSITION=append'
    end if

    if(lu_exists(luname, groupname) .and. is_restarted) then
      ! Open an existing file and put the pointer at the end.
      if (root_and_verbose) then
        write(*,*) luname, ' is found. Will append to it.'
      end if
      open(file_unit_list(lun), FILE=luname, FORM='formatted', &
         & STATUS='old', POSITION='append', IOSTAT=iostatus, &
         & ACTION='readwrite', RECL=record_length)

    else
      if (root_and_verbose) then
        ! Overwrite an existing file or create a new one.
        write(*,*) luname, &
           &' not found or this is not a restarted run. Will create new.'
      end if

      open(file_unit_list(lun), FILE=luname, FORM='formatted', &
         & STATUS='replace', POSITION='rewind', IOSTAT=iostatus, &
         & ACTION='write',RECL=record_length)!, ACCESS='sequential')
    end if


  end subroutine open_lu

  
  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine write_lu_metadata(luname, groupname, type, dims, ishyperslab)
    character (len=*), intent(in) :: luname, groupname, type
    integer, dimension(:), intent(in) :: dims
    logical, intent(in) :: ishyperslab
    
    call open_metadata_file
    write(meta_file_unit,*)
    write(meta_file_unit,*)
    write(meta_file_unit,'(A,A)') 'logical unit name : ', luname
    write(meta_file_unit,'(A,A)') 'group : ', groupname
    write(meta_file_unit,'(A,A)') 'type : ', type
    if(ishyperslab) then
      write(meta_file_unit,'(A,1x,A,1x,(512(I5,1x)))') 'dimensions : ', &
         & 'n', dims
    else
      write(meta_file_unit,'(A,1x,(512(I5,1x)))') 'dimensions : ', &
         & dims
    end if
    close(meta_file_unit)
    
  end subroutine write_lu_metadata
    


  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine open_complex_lu(luname, groupname, hyperslab_dims, lun)
    use global, only : ss
    character (len=*), intent(in) :: luname
    ! the ascii IO ignores the groupname here
    character (len=*), intent(in) :: groupname
    integer, dimension(:), intent(in) :: hyperslab_dims
    integer, intent(out) :: lun
    integer :: real_lun, imag_lun
    if(is_not_root_proc()) return

    if(io_legacy) then
      call open_lu(trim(luname)//real_suffix, groupname, hyperslab_dims, &
         & real_lun)
      call open_lu(trim(luname)//imag_suffix, groupname, hyperslab_dims, &
         & imag_lun)
      call write_lu_metadata(trim(luname)//real_suffix//' and ' &
         & //trim(luname)//imag_suffix, &
         & groupname, 'complex', hyperslab_dims, .true.)
    else
      call open_lu(trim(ss(luname))//real_suffix//suffix, groupname, &
         & hyperslab_dims, real_lun)
      call open_lu(trim(ss(luname))//imag_suffix//suffix, groupname, &
         & hyperslab_dims, imag_lun)
      call write_lu_metadata(trim(ss(luname))//real_suffix//suffix// &
         & ' and '//trim(ss(luname))//imag_suffix//suffix, &
         & groupname, 'complex', hyperslab_dims, .true.)
    end if

    ! the returned lun will be the one corresponding to the real part
    lun = real_lun

    ! set an entry which can be tested when closing the dataset.
    unit_is_complex(lun) = .true.
    unit_is_complex(imag_lun) = .false.

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
    ! Close also the associated dataset with the imaginary part.
    call close_real_lu(lun+1)
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
      inquire(unit=file_unit_list(i), opened=is_open)
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
    character (len=*), intent(in) :: luname
    ! the ascii IO ignores the groupname here
    character (len=*), intent(in) :: groupname
    logical :: lu_exists

    if (len(groupname) /= 0) continue

    if(is_not_root_proc()) then
      lu_exists = .false.
    else
      if(io_legacy) then
        inquire(file=trim(luname), exist = lu_exists)
      else
        inquire(file=trim(ss(luname))//suffix, exist = lu_exists)
      end if
    end if

  end function lu_exists

  !----------------------------------------------------------------------------
  !> Obtain a logical unit number which is associated to a free file unit
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
      call gkw_throw_abort('ascii IO: Increase max_file_units, or use less diagnostics')
      return
    end if
    
  end subroutine get_free_lun

  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine get_free_file_unit(file_unit)
    use global, only : gkw_throw_abort
    integer, intent(out) :: file_unit
    logical :: is_open

    ! In restart.F90 this routine is called by 'the last process'.
    ! Some diagnostics (distr{3,4}.dat, fluxes_vspace..) write from
    ! non-root-processes.
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

    call gkw_throw_abort('ascii IO: Increase max_file_units, or use less diagnostics')

  end subroutine get_free_file_unit


  !----------------------------------------------------------------------------
  !> This function attempts to imitate the flush intrinsic in Fortran 2003. 
  !> It if causes difficulties for a given compiler, the contents can simply be
  !> commented and the function made empty; flushing files to disk is only for 
  !> convenience and the code will run fine without it.
  !----------------------------------------------------------------------------
  subroutine flush_file(lun)
    use global, only : i8 

    integer, intent(in)              :: lun
    integer, parameter :: max_filename_len = 512
    character (len=max_filename_len) :: fname
    ! The integer below is kind=i8 for the IBM compiler (see r3446, Issue 203)
    ! However, the Intel compiler complains that this violates the Fortran95
    ! standard. Unless an alterative completely portable solution is found 
    ! it is better to break the Intel interpretation of the standard than 
    ! to have something that won't run with a known, in-use compiler.
    ! You can remove the kind specification in your local copy if required.
    integer (kind=i8)                :: recll

    if(is_not_root_proc()) return

    if (file_is_open(lun)) then
#ifdef STD2003_FLUSH
      flush(file_unit_list(lun))
#else
      inquire(UNIT=file_unit_list(lun), NAME=fname)
      inquire(UNIT=file_unit_list(lun), RECL=recll)
      close(file_unit_list(lun))

      ! workaround for gfortran bug:
      ! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53796
      !write (*,*) "ascii file", fname, "is flushed, recll=", recll
      if (recll < 1) then
        open(file_unit_list(lun), FILE=fname, &
           & FORM='formatted', POSITION='append')
      else
        open(file_unit_list(lun), FILE=fname, &
           & FORM='formatted', POSITION='append', RECL=recll)
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
      inquire(unit=file_unit_list(lun), opened=file_is_open)
    end if
  end function file_is_open

  
  !----------------------------------------------------------------------------
  !> This module has an array file_2nd_kind_groupnames where all files are
  !> which are to be realised as "files of 2nd kind" (see comments on top).
  !----------------------------------------------------------------------------
  function get_file_2nd_kind_index(groupname)
    character (len=*), intent(in) :: groupname
    integer :: get_file_2nd_kind_index
    integer :: i

    get_file_2nd_kind_index = -1
    do i = 1, size(file_2nd_kind_groupnames)
      if(file_2nd_kind_groupnames(i) == groupname) then
        get_file_2nd_kind_index = i
        exit
      end if
    end do

  end function get_file_2nd_kind_index


  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine open_file_to_append_or_replace(filename, groupname, file_unit)
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    integer, intent(in) :: file_unit
    
    logical :: file_exists
    integer :: i, ios

    ! FIXME the fluxes_vspace diagnostic writes from non-root processes
    !if(is_not_root_proc()) return
    
    i = get_file_2nd_kind_index(groupname)
    if(i > 0) then
      ! This is a file of the 2nd kind.

      inquire(file=trim(groupname)//'.dat', exist = file_exists)
      
      if(file_exists .and. .not.file_2nd_kind_is_touched(i)) then
        ! The file exists already and was not yet written to in this
        ! run. We overwrite it.
        if (root_and_verbose) then
          write(*,*) filename, ' is in ', groupname, &
             & ' which is found. Will overwrite it.'
        end if
        open(file_unit, FILE=trim(groupname)//suffix, FORM='formatted', &
           & STATUS='replace', POSITION='rewind')
      elseif(file_exists .and. file_2nd_kind_is_touched(i)) then
        ! The file exists already but was already written to in this
        ! run. We append to it.
        if (root_and_verbose) then
          write(*,*) filename, ' is in ', groupname, &
             & ' which is found. Will append to it.'
        end if
        open(file_unit, FILE=trim(groupname)//suffix, FORM='formatted', &
           & STATUS='old', POSITION='append')
      else
        ! The file does not exist yet. We create a new one.
        if (root_and_verbose) then
          write(*,*) filename, ' is in ', groupname, &
             & ' which is not found. Will create new.'
        end if
        open(file_unit, FILE=trim(groupname)//suffix, FORM='formatted', &
           & STATUS='new', POSITION='rewind')
      end if

      ! anyway, now it is clear that this file shall not be
      ! overwritten later during this run:
      file_2nd_kind_is_touched(i) = .true.
      
    else
      ! This is a file of the 1st kind:
      ! Overwrite the file called filename

      ! FIXME Does the RECL argument here actually make sense?
      ! open(file_unit, FILE=filename, FORM='formatted', &
         ! & STATUS='replace', POSITION='rewind', RECL=max_cols*fmt_wid)
      ! FIXME Before SRG rewrote this, the only difference between
      ! output_slice and output_array was the presence/absence of RECL
      if(io_legacy) then
        if (root_and_verbose) then
          write(*,*) trim(filename),' will be overwritten, if it exists.'
        end if
        open(file_unit, FILE=filename, FORM='formatted', &
           & STATUS='replace', POSITION='rewind',IOSTAT=ios)
        if(ios /= 0) then
          write(*,*) 'Error while opening ',filename,', iostat=',ios
          write(*,*) 'file_unit=', file_unit
        end if
      else
        if (root_and_verbose) then
          write(*,*) trim(filename)//suffix,' will be overwritten, if it exists.'
        end if
        open(file_unit, FILE=trim(filename)//suffix, FORM='formatted', &
           & STATUS='replace', POSITION='rewind',IOSTAT=ios)
        if(ios /= 0) then
          write(*,*) 'Error while opening ',trim(filename)//suffix,', iostat=',ios
          write(*,*) 'file_unit=', file_unit
        end if
      end if

    end if

  end subroutine open_file_to_append_or_replace

  
  !----------------------------------------------------------------------------
  !> Output a real 1D array data(:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_1d(filename, groupname, data, order, fmt)
    real, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    integer :: file_unit
    integer :: i
    integer :: dim1

    if (len(order) /= 0) continue

    call get_free_file_unit(file_unit)

    dim1=size(data,1)

    call open_file_to_append_or_replace(filename, groupname, file_unit)
    i = get_file_2nd_kind_index(groupname)
    if(i > 0) then
      ! write (*,*) 'file of 2nd kind:', filename, groupname, file_unit
      write(file_unit, FMT='(A)') filename
    end if
    
    write(file_unit, FMT=fmt) (data(i), i=1, dim1)
    close(file_unit)

  end subroutine output_array_1d

  !----------------------------------------------------------------------------
  !> Output a real 1D array data(:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_1d_real(filename, groupname, data, order, fmt)
    real, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_1d(filename, groupname, data, order, fmt)
    
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), &
       & .false.)
  end subroutine output_array_1d_real

  !----------------------------------------------------------------------------
  !> Output a real 2D array data(:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_2d(filename, groupname, data, order, fmt)
    real, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    integer :: file_unit
    integer :: i, j
    integer :: dim1, dim2

    call get_free_file_unit(file_unit)

    dim1=size(data,1)
    dim2=size(data,2)

    call open_file_to_append_or_replace(filename, groupname, file_unit)
    i = get_file_2nd_kind_index(groupname)
    if(i > 0) then
      ! write (*,*) 'file of 2nd kind:', filename, groupname, file_unit
      write(file_unit, FMT='(A)') filename
    end if
    
    
    if (order == 'C') then
      do i=1, dim1
        write(file_unit, FMT=fmt) (data(i, j), j=1, dim2)
      end do
    else
      do j=1, dim2
        write(file_unit, FMT=fmt) (data(i, j), i=1, dim1)
      end do
    end if

    close(file_unit)

  end subroutine output_array_2d

  !----------------------------------------------------------------------------
  !> Output a real 2D array data(:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_2d_real(filename, groupname, data, order, fmt)
    real, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_2d(filename, groupname, data, order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)
  end subroutine output_array_2d_real

  !----------------------------------------------------------------------------
  !> Output a real 3D array data(:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_3d(filename, groupname, data, order, fmt)
    real, dimension(:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    integer :: file_unit
    integer :: dim1, dim2, dim3
    integer :: i, j, k

    call get_free_file_unit(file_unit)

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)

    call open_file_to_append_or_replace(filename, groupname, file_unit)
    i = get_file_2nd_kind_index(groupname)
    if(i > 0) then
      write(file_unit, FMT='(A)') filename
    end if
    
    if (order == 'C') then
      do i=1, dim1
        do j=1, dim2
          write(file_unit, FMT=fmt) (data(i, j, k), k=1, dim3)
        end do
      end do
    else
      do k=1, dim3
        do j=1, dim2
          write(file_unit, FMT=fmt) (data(i, j, k), i=1, dim1)
        end do
      end do
    end if

    close(file_unit)


  end subroutine output_array_3d

  !----------------------------------------------------------------------------
  !> Output a real 3D array data(:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_3d_real(filename, groupname, data, order, fmt)
    real, dimension(:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_3d(filename, groupname, data, order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)
  end subroutine output_array_3d_real

  !----------------------------------------------------------------------------
  !> Output a real 4D array data(:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_4d(filename, groupname, data, order, fmt)
    real, dimension(:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    integer :: file_unit
    integer :: dim1, dim2, dim3, dim4
    integer :: i, j, k, l

    call get_free_file_unit(file_unit)

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)
    dim4=size(data,4)

    call open_file_to_append_or_replace(filename, groupname, file_unit)
    i = get_file_2nd_kind_index(groupname)
    if(i > 0) then
      write(file_unit, FMT='(A)') filename
    end if
    
    if (order == 'C') then
      do i=1, dim1
        do j=1, dim2
          do k=1, dim3
            write(file_unit, FMT=fmt) (data(i, j, k, l), l=1, dim4)
          end do
        end do
      end do
    else
      do l=1, dim4
        do k=1, dim3
          do j=1, dim2
            write(file_unit, FMT=fmt) (data(i, j, k, l), i=1, dim1)
          end do
        end do
      end do
    end if

    close(file_unit)


  end subroutine output_array_4d

  !----------------------------------------------------------------------------
  !> Output a real 4D array data(:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_4d_real(filename, groupname, data, order, fmt)
    real, dimension(:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_4d(filename, groupname, data, order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)
  end subroutine output_array_4d_real

  !----------------------------------------------------------------------------
  !> Output a real 5D array data(:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_5d(filename, groupname, data, order, fmt)
    real, dimension(:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    integer :: file_unit
    integer :: dim1, dim2, dim3, dim4, dim5
    integer :: i, j, k, l, m

    call get_free_file_unit(file_unit)

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)
    dim4=size(data,4)
    dim5=size(data,5)

    call open_file_to_append_or_replace(filename, groupname, file_unit)
    i = get_file_2nd_kind_index(groupname)
    if(i > 0) then
      write(file_unit, FMT='(A)') filename
    end if
    
    if (order == 'C') then
      do i=1, dim1
        do j=1, dim2
          do k=1, dim3
            do l=1, dim4
              write(file_unit, FMT=fmt) (data(i, j, k, l, m), m=1, dim5)
            end do
          end do
        end do
      end do
    else
      do m=1, dim5
        do l=1, dim4
          do k=1, dim3
            do j=1, dim2
              write(file_unit, FMT=fmt) (data(i, j, k, l, m), i=1, dim1)
            end do
          end do
        end do
      end do
    end if

    close(file_unit)


  end subroutine output_array_5d

  !----------------------------------------------------------------------------
  !> Output a real 5D array data(:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_5d_real(filename, groupname, data, order, fmt)
    real, dimension(:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_5d(filename, groupname, data, order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)
  end subroutine output_array_5d_real

  !----------------------------------------------------------------------------
  !> Output a real 6D array to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_6d(filename, groupname, data, order, fmt)
    real, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    integer :: file_unit
    integer :: dim1, dim2, dim3, dim4, dim5, dim6
    integer :: i, j, k, l, m, n

    call get_free_file_unit(file_unit)

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)
    dim4=size(data,4)
    dim5=size(data,5)
    dim6=size(data,6)

    call open_file_to_append_or_replace(filename, groupname, file_unit)
    i = get_file_2nd_kind_index(groupname)
    if(i > 0) then
      write(file_unit, FMT='(A)') filename
    end if
    
    if (order == 'C') then
      do i=1, dim1
        do j=1, dim2
          do k=1, dim3
            do l=1, dim4
              do m=1, dim5
                write(file_unit, FMT=fmt) (data(i, j, k, l, m, n), n=1, dim6)
              end do
            end do
          end do
        end do
      end do
    else
      do n=1, dim6
        do m=1, dim5
          do l=1, dim4
            do k=1, dim3
              do j=1, dim2
                write(file_unit, FMT=fmt) (data(i, j, k, l, m, n), i=1, dim1)
              end do
            end do
          end do
        end do
      end do
    end if

    close(file_unit)

  end subroutine output_array_6d

  !----------------------------------------------------------------------------
  !> Output a real 6D array file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_6d_real(filename, groupname, data, order, fmt)
    real, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_6d(filename, groupname, data, order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(filename, groupname, 'real', shape(data), .false.)
  end subroutine output_array_6d_real

  
  !----------------------------------------------------------------------------
  !> Output a complex 1D array data(:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_1d_complex(filename, groupname, data, order, fmt)
    complex, dimension(:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_1d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order, fmt)
    call output_array_1d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(trim(filename)//real_suffix//suffix//' and ' &
       & //trim(filename)//imag_suffix//suffix, &
       & groupname, 'complex', shape(data), .false.)

  end subroutine output_array_1d_complex

  !----------------------------------------------------------------------------
  !> Output a complex 2D array data(:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_2d_complex(filename, groupname, data, order, fmt)
    complex, dimension(:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_2d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order, fmt)
    call output_array_2d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(trim(filename)//real_suffix//suffix//' and ' &
       & //trim(filename)//imag_suffix//suffix, &
       & groupname, 'complex', shape(data), .false.)

  end subroutine output_array_2d_complex

  !----------------------------------------------------------------------------
  !> Output a complex 3D array data(:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_3d_complex(filename, groupname, data, order, fmt)
    complex, dimension(:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_3d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order, fmt)
    call output_array_3d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order, fmt)
    ! Spit metadata into a simple text file, too.
    call write_lu_metadata(trim(filename)//real_suffix//suffix//' and ' &
       & //trim(filename)//imag_suffix//suffix, &
       & groupname, 'complex', shape(data), .false.)

  end subroutine output_array_3d_complex


  !----------------------------------------------------------------------------
  !> Output a complex 6D array data(:,:,:,:,:,:) to file, as formatted ascii,
  !> with order 'C' or 'F'.
  !> If the file already exists, it is overwritten.
  !----------------------------------------------------------------------------
  subroutine output_array_6d_complex(filename, groupname, data, order, fmt)
    complex, dimension(:,:,:,:,:,:), intent(in) :: data
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: groupname
    character (len=*), intent(in) :: fmt
    character (len=1), intent(in) :: order

    call output_array_6d(trim(filename)//real_suffix//suffix, groupname, &
       & real(data), order, fmt)
    call output_array_6d(trim(filename)//imag_suffix//suffix, groupname, &
       & aimag(data), order, fmt)
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
  subroutine append_chunk_1d_real(lun, data, fmt)
    real, dimension(:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt
    integer :: i
    integer :: dim1

    if(is_not_root_proc()) return
    
    dim1 = size(data,1)
    
    write(file_unit_list(lun), FMT=fmt) (data(i), i=1, dim1)

  end subroutine append_chunk_1d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 2D matrix data(:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_2d_real(lun, data, fmt)
    real, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt
    integer :: i, j
    integer :: dim1, dim2

    if(is_not_root_proc()) return
    
    dim1 = size(data,1)
    dim2 = size(data,2)

    do j = 1, dim2
      write(file_unit_list(lun), FMT=fmt) (data(i, j), i=1, dim1)
    end do

  end subroutine append_chunk_2d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 3D matrix data(:,:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_3d_real(lun, data, fmt)
    real, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt
    integer :: i, j, k
    integer :: dim1, dim2, dim3

    if(is_not_root_proc()) return
    
    dim1 = size(data,1)
    dim2 = size(data,2)
    dim3 = size(data,3)

    do k = 1, dim3
      do j = 1, dim2
        write(file_unit_list(lun), FMT=fmt) (data(i, j, k), i=1, dim1)
      end do
    end do

  end subroutine append_chunk_3d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 4D matrix data(:,:,:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_4d_real(lun, data, fmt)
    real, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt
    integer :: i, j, k, l
    integer :: dim1, dim2, dim3, dim4

    if(is_not_root_proc()) return
    
    dim1 = size(data,1)
    dim2 = size(data,2)
    dim3 = size(data,3)
    dim4 = size(data,4)

    do l = 1, dim4
      do k = 1, dim3
        do j = 1, dim2
          write(file_unit_list(lun), FMT=fmt) (data(i, j, k, l), i=1, dim1)
        end do
      end do
    end do

  end subroutine append_chunk_4d_real

  !----------------------------------------------------------------------------
  !> This subroutine appends a 6D matrix to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_6d_real(lun, data, fmt)
    real, dimension(:,:,:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt
    integer :: i, j, k, l, m, n
    integer :: dim1, dim2, dim3, dim4, dim5, dim6

    if(is_not_root_proc()) return
    
    dim1 = size(data,1)
    dim2 = size(data,2)
    dim3 = size(data,3)
    dim4 = size(data,4)
    dim5 = size(data,5)
    dim6 = size(data,6)

    do n = 1, dim6
      do m = 1, dim5
        do l = 1, dim4
          do k = 1, dim3
            do j = 1, dim2
              write(file_unit_list(lun), FMT=fmt) &
                 & (data(i, j, k, l, m, n), i=1, dim1)
            end do
          end do
        end do
      end do
    end do

  end subroutine append_chunk_6d_real
  
  !----------------------------------------------------------------------------
  !> This subroutine appends a 1D vector data(:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_1d_complex(lun, data, fmt)
    use global, only : int2char
    complex, dimension(:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt

    if(is_not_root_proc()) return

    if(unit_is_complex(lun)) then
      call append_chunk_1d_real(lun, real(data), fmt)
      call append_chunk_1d_real(lun+1, aimag(data), fmt)
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if
    
  end subroutine append_chunk_1d_complex

  !----------------------------------------------------------------------------
  !> This subroutine appends a 2D matrix data(:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_2d_complex(lun, data, fmt)
    use global, only : int2char
    complex, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt

    if(is_not_root_proc()) return

    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_2d_real(lun, real(data), fmt)
      call append_chunk_2d_real(lun+1, aimag(data), fmt)
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if

  end subroutine append_chunk_2d_complex

!----------------------------------------------------------------------------
  !> This subroutine appends a 2D matrix data(:,:) to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_3d_complex(lun, data, fmt)
    use global, only : int2char
    complex, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt

    if(is_not_root_proc()) return

    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_3d_real(lun, real(data), fmt)
      call append_chunk_3d_real(lun+1, aimag(data), fmt)
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if

  end subroutine append_chunk_3d_complex

!----------------------------------------------------------------------------
  !> This subroutine appends a matrix to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_4d_complex(lun, data, fmt)
    use global, only : int2char
    complex, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt

    if(is_not_root_proc()) return
    
    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_4d_real(lun, real(data), fmt)
      call append_chunk_4d_real(lun+1, aimag(data), fmt)
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if

  end subroutine append_chunk_4d_complex

  !----------------------------------------------------------------------------
  !> This subroutine appends a matrix to a file.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_complex_lu().
  !----------------------------------------------------------------------------
  subroutine append_chunk_6d_complex(lun, data, fmt)
    use global, only : int2char
    complex, dimension(:,:,:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt

    if(is_not_root_proc()) return
    
    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_6d_real(lun, real(data), fmt)
      call append_chunk_6d_real(lun+1, aimag(data), fmt)
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if

  end subroutine append_chunk_6d_complex
  
  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine read_last_chunk_1d_real(lun, fmt, data, nlast)
    use global, only : int2char, gkw_throw_error
    
    integer, intent(in) :: lun
    character (len=*), intent(in) :: fmt
    real, dimension(:), intent(out) :: data
    integer, intent(out) :: nlast
    integer:: ios

    if(is_not_root_proc()) return

    nlast = 0

    ! rewind the file pointer to the beginning of the file
    rewind file_unit_list(lun)
    ! read records and count how many there are
    do
      read(UNIT=file_unit_list(lun), FMT=fmt, IOSTAT=ios) data
      if(ios < 0) then
        ! negative ios means end-of-file
        if(nlast == 0) then
          call gkw_throw_error("Could not read last line from logical unit " &
             &//int2char(lun)//" because the file is empty.")
          return
        else
          exit
        end if
      end if
      if(ios/=0) then
        call gkw_throw_error(&
           &"An error occurred when reading the last line from logical unit " &
           &//int2char(lun)//" (line "//int2char(nlast)//")")
        return
      end if

      nlast = nlast + 1
    end do
    ! backspace to before end-of-file record
    backspace file_unit_list(lun)
    ! Read the last record
    read(UNIT=file_unit_list(lun), FMT=fmt, IOSTAT=ios) data
  end subroutine read_last_chunk_1d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine seek_to_chunk(lun, index)
    integer, intent(in) :: lun
    integer, intent(in) :: index

    if(.false.) write(*,*) lun, ' ', index

    continue

  end subroutine seek_to_chunk

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_string(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val

    call open_metadata_file
    write(meta_file_unit,'(A,A)') key//' : ', val
    close(meta_file_unit)

    if (.false.) write(*,*) lun

  end subroutine attach_metadata_string

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_real(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    real, intent(in) :: val
    
    call open_metadata_file
    write(meta_file_unit,'(A,es13.5)') key//' : ', val
    close(meta_file_unit)

    if (.false.) write(*,*) lun

  end subroutine attach_metadata_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_integer(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    integer, intent(in) :: val

    call open_metadata_file
    write(meta_file_unit,'(A,I8)') key//' : ', val
    close(meta_file_unit)

    if (.false.) write(*,*) lun

  end subroutine attach_metadata_integer

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_logical(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    logical, intent(in) :: val

    call open_metadata_file
    if(val) then
      write(meta_file_unit,'(A,A)') key//' : ', 'T'
    else
      write(meta_file_unit,'(A,A)') key//' : ', 'F'
    end if
    close(meta_file_unit)

    if (.false.) write(*,*) lun

  end subroutine attach_metadata_logical

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_r_array(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val

    call open_metadata_file
    write(meta_file_unit,'(A,(512(es13.5,1x)))') key//' : ', val
    close(meta_file_unit)

    if (.false.) write(*,*) lun

  end subroutine attach_metadata_r_array

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_i_array(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    call open_metadata_file
    write(meta_file_unit,'(A,(512(I8,1x)))') key//' : ', val
    close(meta_file_unit)

    if (.false.) write(*,*) lun

  end subroutine attach_metadata_i_array

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_string_name(luname, groupname, key, val)
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val

    call attach_metadata_string(-1, key, val)

    if (.false.) write(*,*) luname, ' ', groupname

  end subroutine attach_metadata_string_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_real_name(luname, groupname, key, val)
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    real, intent(in) :: val

    call attach_metadata_real(-1, key, val)

    if (.false.) write(*,*) luname, ' ', groupname

  end subroutine attach_metadata_real_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_integer_name(luname, groupname, key, val)
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    integer, intent(in) :: val

    call attach_metadata_integer(-1, key, val)

    if (.false.) write(*,*) luname, ' ', groupname

  end subroutine attach_metadata_integer_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_logical_name(luname, groupname, key, val)
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    logical, intent(in) :: val

    call attach_metadata_logical(-1, key, val)

    if (.false.) write(*,*) luname, ' ', groupname

  end subroutine attach_metadata_logical_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_r_array_name(luname, groupname, key, val)
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val

    call attach_metadata_r_array(-1, key, val)

    if (.false.) write(*,*) luname, ' ', groupname

  end subroutine attach_metadata_r_array_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_i_array_name(luname, groupname, key, val)
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    call attach_metadata_i_array(-1, key, val)

    if (.false.) write(*,*) luname, ' ', groupname

  end subroutine attach_metadata_i_array_name

end module io_ascii

