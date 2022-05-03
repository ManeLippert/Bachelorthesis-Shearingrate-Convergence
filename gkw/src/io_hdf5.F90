!------------------------------------------------------------------------------
!> HDF5 output routines for the swappable interfaces in the IO module.
!>
!> In terms of this module, a logical unit for real data refers to
!> * one HDF5 dataset
!> or
!>  * an arbitrary number of HDF5 datasets ("sub-logical units")
!> while a logical unit for complex data refers to
!>  * exactly two HDF5 datasets.
!>
!> If one creates a real logical unit using open_real_dataset_hdf5
!> and passes a list of labels as an argument, then the logical unit
!> will be split up along the first dimension: For every hyperslab e.g.
!> data(i, :, :) a separate dataset is created.
!>
!> Notes About HDF5:
!>  * 'Attaching metadata' to a dataset or group means creating HDF5
!>    attributes.  An HDF5 attribute is a metadata object describing
!>    the nature and/or intended usage of a primary data object (a
!>    dataset, group, ...). For various reasons (e.g. there is no
!>    compression), attributes should be *very small*.
!>  * Octave can read in the whole data using
!>        data=load('-hdf5', 'gkwdata.h5')
!>    One can load only certain parts with
!>        data=load('-hdf5', 'gkwdata.h5', 'diagnostic', 'geom')
!>  * If octave complains that the HDF5 libs you compiled yourself
!>    are more recent than the system libs it was compiled against, try
!>    starting it with something like
!>        LD_LIBRARY_PATH=/usr/lib64/ yourscript.m
!>  * Octave does not seem to be able to read HDF5 attributes at the
!>    moment.
!------------------------------------------------------------------------------
module io_hdf5
  use global, only : ss
#ifdef HAVE_HDF5

! define Macro for comparing versions, not introduced until 1.8.7
#define H5_VERSION_GE(Maj,Min,Rel)  \
        (((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR==Min) && (H5_VERS_RELEASE>=Rel)) || \
         ((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR>Min)) || \
         (H5_VERS_MAJOR>Maj)) 

  ! UPPERCASE USE to deliberately exclude from mkdeps script
  USE hdf5

  implicit none

  private

  ! general
  public :: get_libversion, init, finalize
  public :: flush_all_open_files, flush_file
  public :: dataset_exists_hdf5

  ! subroutines to open, write (as often as desired) and close
  ! extendible datasets:
  public :: open_real_dataset_hdf5, open_complex_dataset_hdf5
  public :: close_dataset_hdf5
  public :: close_all_lus
  public :: append_chunk, seek_to_chunk
  public :: read_last_chunk
  public :: attach_metadata


  ! subroutines which open, write (once and for all) and close fixed
  ! size datasets:
  public :: output_array

  ! subroutines which deal with metadata:
  public :: write_run_parameter


  logical, save :: is_initialized = .false.
  
  logical, save :: with_szip_compression = .false.
  integer, parameter :: MIN_PIXELS_PER_BLOCK = 8

  integer(HID_T), save :: APPROPRIATE_HDF5_REAL_TYPE
  
  integer, parameter :: max_dsets = 256
  integer, save :: n_open_dsets = 0
  integer, parameter :: I_DSET = 1
  integer, parameter :: I_GROUP = 2
  integer(HID_T), save, dimension(max_dsets, 2) :: dset_id_list
  logical, save, dimension(max_dsets) :: unit_is_complex

  !> A logical unit may refer to several other real datasets. This
  !> datastructure is used for this purpose: (splitted_luns_list is
  !> like a ragged matrix of integers, i.e. the rows have different
  !> length)
  type ragged
    ! unfortunately this is Fortran 2003:
    !integer, dimension(:), allocatable :: lun_sublist
    ! this is Fortran 95:
    integer, dimension(:), pointer :: lun_sublist
  end type ragged
  type(ragged), dimension(max_dsets) :: splitted_luns_list

  ! File identifier and name
  integer(HID_T), save :: main_output_file_id
  character(len=10), parameter :: main_output_file_name = "gkwdata.h5"

  ! Data group identifier and name
  character(len=11), parameter :: standard_group_name = "diagnostic"

  character(len=5), parameter :: input_group_name = "input"

  character(len=5), parameter :: real_suffix = "_real"
  character(len=5), parameter :: imag_suffix = "_imag"
  
  character(len=128) :: dset_name_buf
  integer(size_t) :: dset_name_buf_size = 128, name_size

  !> gets a copy of mpiinterface value
  integer, save :: proc_number
  logical, save :: is_restarted

  interface open_real_dataset_hdf5
    module procedure open_real_single_dataset_hdf5
    module procedure open_real_splitted_dataset_hdf5
  end interface open_real_dataset_hdf5
  
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

contains

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  function get_libversion()
    use global, only : int2char
    character (len=10) :: get_libversion
    integer :: majnum, minnum, relnum, err
    call h5get_libversion_f(majnum, minnum, relnum, err)
    get_libversion = trim(int2char(majnum,2))//'.'//trim(int2char(minnum,2))//'.'&
       & //trim(int2char(relnum,2))
  end function get_libversion

  !---------------------------------------------------------------------------
  !> A tiny helper function which is useful to make sure that only the
  !> root process does IO.
  !---------------------------------------------------------------------------
  function is_not_root_proc()
    logical :: is_not_root_proc
    is_not_root_proc = (proc_number/=0)
    if(is_not_root_proc) then
      !call gkw_throw_warn("Parallel HDF5 IO is not implemented yet.")
      continue
      ! actually, this warning should pop up rather seldom, if the
      ! rest of the code lets only the root proc do IO.
    end if
  end function is_not_root_proc
  
  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine handle_hdf5_error(err, filename, line)
    use global, only : int2char, gkw_throw_warn
    integer, intent(in) :: err, line
    character(len=*) :: filename
    ! an integer in case of yet another internal error
    integer :: ierr
    if(err<0) then
      write (*,*) trim(filename)//":"//int2char(line,4) &
         & //" HDF5 err = "//int2char(err,2)
#define HDF5_ERROR_LOGFILE "gkw_hdf5_errors.txt"
      call h5eprint_f(ierr, HDF5_ERROR_LOGFILE)
      write (*,*) 'This error has been logged to', HDF5_ERROR_LOGFILE
      ! call gkw_throw_warn(trim(filename)//":"//int2char(line,4) &
      !    & //" HDF5 err = "//int2char(err,2))
    end if
  end subroutine handle_hdf5_error
!#define check_err(err)  if(err<0) write (*,*) " HDF5 err at ",__FILE__,__LINE__
#define check_err(err)  if(err<0) call handle_hdf5_error(err,__FILE__,__LINE__)

  !---------------------------------------------------------------------------
  !> This subroutine initializes the HDF5 interface by
  !>  a) opening the output file if it already exists or
  !>  b) creating a new file otherwise. Moreover, some HDF5 groups
  !>     are created.
  !> 
  !---------------------------------------------------------------------------
  subroutine init(proc_number_, is_restarted_)
    use global, only : gkw_throw_abort
    use global, only : have_double_precision, root_and_verbose
    integer, intent(in) :: proc_number_
    logical, intent(in) :: is_restarted_
    ! error flag
    integer :: err
    logical :: status

    !> the size in bytes of each file member. This size will be saved
    !> in file when the property list fapl_id is used to create a new
    !> file. If fapl_id is used to open an existing file, memb_size has
    !> to be equal to the original size saved in file.
    ! integer(HSIZE_T), parameter :: family_file_member_size = &
    !   & H5F_FAMILY_DEFAULT_F
    ! integer(HSIZE_T), parameter :: family_file_member_size = &
    !   & 1024*1024*1024
    integer(HID_T) :: group_id
    integer(HID_T) :: access_plist
    integer :: filter_info_bitfield, filter_info_both
    
    
    proc_number = proc_number_
    is_restarted = is_restarted_

    if(is_not_root_proc()) return

    if(.not. is_initialized) then
      ! Initialize the Fortran interface.
      call h5open_f(err)

      ! suppress hdf5 error output, to avoid some unimportant noise in the log
      call h5eset_auto_f(0,err);

      !  Check if SZIP compression is available and can be used for both
      !  compression and decompression. This filter is an
      !  optional part of the hdf5 library.
      call h5zfilter_avail_f(H5Z_FILTER_SZIP_F, status, err)

      if (.not.status) then
        if(root_and_verbose) write(*,*) "SZIP compression is not&
           & available in the HDF5 library."
      else
        call h5zget_filter_info_f(H5Z_FILTER_SZIP_F, filter_info_bitfield, err)

        filter_info_both = ior(H5Z_FILTER_ENCODE_ENABLED_F, &
           & H5Z_FILTER_DECODE_ENABLED_F)
        if (filter_info_bitfield .ne. filter_info_both) then
          if (root_and_verbose) write (*,*) "SZIP compression is not&
             & available for encoding and/or decoding."
        else
          ! Comment out this line in order to switch off SZIP
          ! compression even if it would be available.
          !with_szip_compression = .true.
          
          if(with_szip_compression) then
            write (*,*) "SZIP compression is available and is used."
          else
            write (*,*) "SZIP compression is available but is not used."
          end if
        end if
      end if
      ! restore old setting
      call h5eset_auto_f(1,err);

      ! We have to choose the type used for floating point
      ! numbers, according to the precision GKW was compiled with.
      if(have_double_precision) then
        if (root_and_verbose) then
          write (*,*) "HDF5 double precision data is treated as H5T_NATIVE_DOUBLE"
        end if
        call h5tcopy_f(H5T_NATIVE_DOUBLE, APPROPRIATE_HDF5_REAL_TYPE, err)
        check_err(err)
      else
        if (root_and_verbose) then
          write (*,*) "HDF5 single precision data is treated as H5T_NATIVE_REAL"
        end if
        call h5tcopy_f(H5T_NATIVE_REAL, APPROPRIATE_HDF5_REAL_TYPE, err)
        check_err(err)
      end if

      call h5pcreate_f(H5P_FILE_ACCESS_F, access_plist, err)
      ! Choose a file driver underlying the HDF5 library by
      ! uncommenting *one* of the following:
      !   1) Unbuffered output to a single file.
      call h5pset_fapl_sec2_f(access_plist, err)

      !   2) Buffered output to a single file.
      ! call h5pset_fapl_stdio_f(access_plist, err)

      !   3) (Un/buffered??) output to a single HDF5, consisting of
      !      several physical files. *Use this to avoid files
      !      larger than 2GB*
      ! call h5pset_fapl_family_f(access_plist, family_file_member_size, &
      !    & H5P_DEFAULT_F, err)

      ! Check if there is already a file, then open it or create a new one.
      inquire(file=main_output_file_name,exist=status)
      if (status .and. is_restarted) then
        ! Test if the existing file is in HDF5 format
        call h5fis_hdf5_f(main_output_file_name, status, err)
        if(.not. status) then
          call gkw_throw_abort(main_output_file_name//" is not a valid HDF5 file")
          return
        else
          ! Open the file
          call h5fopen_f(main_output_file_name, H5F_ACC_RDWR_F, &
             & main_output_file_id, err, access_plist)
          call h5pclose_f(access_plist, err)
        end if
      else
        ! Create a new file.
        ! An existing file with the same name is overwritten.
        call h5fcreate_f(main_output_file_name, H5F_ACC_TRUNC_F, &
           & main_output_file_id, err, H5P_DEFAULT_F, access_plist)

        call h5pclose_f(access_plist, err)
        
        ! Force the creation of certain groups. This is the occasion to add
        ! comments or attributes to them.
        
        ! Create a group in the file which will contain the diagnostic
        ! datasets.
        call h5gcreate_f(main_output_file_id, standard_group_name, &
           & group_id, err)
#if H5_VERSION_GE(1,8,11)
         call h5oset_comment_f(group_id, &
            & "This general group contains datasets produced by various diagnostics.", &
            & err)
#endif
        call h5gclose_f(group_id, err)

        call h5gcreate_f(main_output_file_id, "geom", &
           & group_id, err)
#if H5_VERSION_GE(1,8,11)
         call h5oset_comment_f(group_id, &
            & "This group contains geometry tensors.", &
            & err)
#endif
        call h5gclose_f(group_id, err)

        call h5gcreate_f(main_output_file_id, input_group_name, &
           & group_id, err)
#if H5_VERSION_GE(1,8,11)
         call h5oset_comment_f(group_id, &
            & "This group contains the input parameters, additionally as datasets.", &
            & err)
#endif
        call h5gclose_f(group_id, err)

      end if
      is_initialized = .true.
      ! write (*,*) "HDF 5 initialized"
    end if
  end subroutine init


  !---------------------------------------------------------------------------
  !> This subroutine finalizes the HDF5 interface and tries to close properly
  !> all remaining open objects.
  !---------------------------------------------------------------------------
  subroutine finalize()
    use global, only : gkw_throw_warn
    ! error flag
    integer :: err
    integer(SIZE_T) :: obj_count

    if(is_not_root_proc()) return

    if(is_initialized) then
      ! It is a bit unclear how we should finalize HDF5 properly: In
      ! principle h5close_f() causes a general shutdown of the
      ! library: "all data is written to disk, all identifiers are
      ! closed, and all memory used by the library is cleaned up" (say
      ! the docs).  However, in general it is often a good idea to
      ! close all objects individually, 'manually'. Therefore,
      ! each diagnostic should call its close_lu().

      ! We check here if some piece of code has forgotten to close its
      ! logical unit.  This amounts here to checking a list of open datasets
      ! and a list of open groups.

      call h5fget_obj_count_f(main_output_file_id, H5F_OBJ_DATASET_F, &
         & obj_count, err)
      if(obj_count/=0) then
        call gkw_throw_warn("HDF5 is to be finalized but there are still open datasets.")

        ! Close all datasets which are still open, because someone has
        ! forgotten to close them.
        call close_all_datasets_hdf5()
      end if

      ! Close any other open group objects.
      call h5fget_obj_count_f(main_output_file_id, H5F_OBJ_GROUP_F, &
         & obj_count, err)
      if(obj_count/=0) then
        call gkw_throw_warn("HDF5 is to be finalized but there are still open groups.")

        ! Close all groups which are still open.
        call close_all_groups_hdf5()
      end if

      ! Close the data type
      call h5tclose_f(APPROPRIATE_HDF5_REAL_TYPE, err)

      ! Close the file.
      call h5fclose_f(main_output_file_id, err)
      ! Close the Fortran interface.
      call h5close_f(err)
      is_initialized = .false.

      ! write (*,*) "HDF 5 finalized"
    end if
  end subroutine finalize


  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine close_all_lus()

    if(is_not_root_proc()) return

    call close_all_datasets_hdf5()
    call close_all_groups_hdf5()

  end subroutine close_all_lus

  !---------------------------------------------------------------------------
  !> Close all datasets.
  !---------------------------------------------------------------------------
  subroutine close_all_datasets_hdf5()
    use global, only : gkw_throw_abort, root_and_verbose
    integer :: err
    integer(HID_T), dimension(:), allocatable :: obj_ids
    integer(SIZE_T) :: obj_count
    ! FIXME Can this be implemented without a hardcoded buffersize?
    integer(SIZE_T), parameter :: dsetname_bufsize = 128
    integer(SIZE_T) :: dsetname_size
    character (len=dsetname_bufsize) :: dsetname
    integer :: i

    if(is_not_root_proc()) return

    call h5fget_obj_count_f(main_output_file_id, H5F_OBJ_DATASET_F, &
       & obj_count, err)
    if(obj_count/=0) then
      ! There are open datasets.
      allocate(obj_ids(obj_count), stat=err)
      if (err /= 0) then
        call gkw_throw_abort('could not allocate obj_ids')
        return
      end if
      ! Get a list of all open datasets.
      call h5fget_obj_ids_f(main_output_file_id, H5F_OBJ_DATASET_F, &
         & obj_count, obj_ids, err)
      check_err(err)
      do i = 1,size(obj_ids)
        ! Close each one and log it.
        call h5iget_name_f(obj_ids(i), dsetname, dsetname_bufsize, &
           & dsetname_size, err)
        call h5dclose_f(obj_ids(i), err)
        if (root_and_verbose) then
          write (*,*) "The dataset '"//dsetname(1:dsetname_size)//"' was closed."
        end if
      end do
      deallocate(obj_ids)
    end if

  end subroutine close_all_datasets_hdf5

  !---------------------------------------------------------------------------
  !> Close all groups.
  !---------------------------------------------------------------------------
  subroutine close_all_groups_hdf5()
    use global, only : gkw_throw_abort, root_and_verbose
    integer :: err
    integer(HID_T), dimension(:), allocatable :: obj_ids
    integer(SIZE_T) :: obj_count
    ! FIXME Can this be implemented without a hardcoded buffersize?
    integer(SIZE_T), parameter :: groupname_bufsize = 128
    integer(SIZE_T) :: groupname_size
    character (len=groupname_bufsize) :: groupname
    integer :: i

    if(is_not_root_proc()) return

    call h5fget_obj_count_f(main_output_file_id, H5F_OBJ_GROUP_F, &
       & obj_count, err)
    if(obj_count/=0) then
      ! There are open groups.
      allocate(obj_ids(obj_count), stat=err)
      if (err /= 0) then
        call gkw_throw_abort('could not allocate obj_ids')
        return
      end if
      ! Get a list of all open groups.
      call h5fget_obj_ids_f(main_output_file_id, H5F_OBJ_GROUP_F, &
         & obj_count, obj_ids, err)
      check_err(err)
      do i = 1,size(obj_ids)
        ! Close each one and log it.
        call h5iget_name_f(obj_ids(i), groupname, groupname_bufsize, &
           & groupname_size, err)
        call h5gclose_f(obj_ids(i), err)
        if (root_and_verbose) then
          write (*,*) "The group '"//groupname(1:groupname_size)//"' was closed."
        end if
      end do
      deallocate(obj_ids)
    end if

  end subroutine close_all_groups_hdf5

  !---------------------------------------------------------------------------
  !> Check wether a dataset with the given name exists.
  !---------------------------------------------------------------------------
  function dataset_exists_hdf5(dsetname, groupname)
    use global, only : ss
    character (len=*), intent(in) :: dsetname
    character (len=*), intent(in) :: groupname
    logical :: dataset_exists_hdf5
    integer :: err
    integer(HID_T) :: group_id

    if(is_not_root_proc()) then
      dataset_exists_hdf5 = .false.
    elseif(.not.group_exists_hdf5(ss(groupname))) then
      dataset_exists_hdf5 = .false.
    else
      ! open the specified group
      call h5gopen_f(main_output_file_id, ss(groupname), &
         & group_id, err)
      check_err(err)

      ! check if a link with the name given by dsetname exists in the
      ! specified group.  (We assume that the link points to a
      ! dataset - this is not checked.)
      call h5lexists_f(group_id, ss(dsetname), &
         & dataset_exists_hdf5, err)
      check_err(err)

      ! FIXME Should one check here if its datatype/precision
      ! matches APPROPRIATE_HDF5_REAL_TYPE ?

      call h5gclose_f(group_id, err)
    end if

  end function dataset_exists_hdf5

  !---------------------------------------------------------------------------
  !> Check wether a group with the given name exists.
  !---------------------------------------------------------------------------
  function group_exists_hdf5(groupname)
    use global, only : ss
    character (len=*), intent(in) :: groupname
    logical :: group_exists_hdf5
    integer :: err

    if(is_not_root_proc()) then
      group_exists_hdf5 = .false.
    else
      call h5lexists_f(main_output_file_id, ss(groupname), &
         & group_exists_hdf5, err)
      check_err(err)
    end if

  end function group_exists_hdf5


  !---------------------------------------------------------------------------
  !> This subroutine
  !>  a) creates a new dataset 'dsetname' in the diagnostics group of the 
  !>     output HDF5 file, if that dataset does not yet exist
  !>  b) opens the dataset 'dsetname' preserving the old data, if it does
  !>     already exist.
  !> To create the dataset, the desired datatype must already be known, hence
  !> the _real_ and the _complex_ version.
  !> This module remembers the chosen datatype and therefore provides a single
  !> routine close_dataset_hdf5().
  !>
  !> This simplifying output wrapper expects the array hyperslab_dims.
  !> It creates an extendible dataset with 1+size(hyperslab_dims) dimensions.
  !> Using the subroutine append_chunk(), one can append chunks of the
  !> shape described by hyperslab_dims.
  !---------------------------------------------------------------------------
  subroutine open_real_single_dataset_hdf5(dsetname, groupname, hyperslab_dims, &
     & lun)
    use global, only : gkw_throw_abort, gkw_its_critical, root_and_verbose
    use global, only : int2char
    character (len=*), intent(in) :: dsetname
    integer, dimension(:), intent(in) :: hyperslab_dims
    character (len=*), intent(in) :: groupname
    integer, intent(out) :: lun

    integer :: rank
    ! error flag
    integer :: err
    logical :: link_exists
    integer(HSIZE_T), dimension(:), allocatable :: maxdims
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(HSIZE_T), dimension(:), allocatable :: dims_chunk

    integer(HID_T) :: crp_list  ! Dataset creation property identifier
    integer(HID_T) :: dspace_id ! Dataspace identifier

    integer :: pixels_per_block

    if(is_not_root_proc()) return

    ! FIXME Check if it causes trouble when accessing the data from octave,
    ! if the datasetname contains a dot like 'time.dat'

    rank = size(hyperslab_dims)+1
    if (root_and_verbose) then
      write (*,*) "Create extendible HDF5 dataset ",groupname//'/'//dsetname, &
         & ", rank ", rank
    end if

    allocate(maxdims(rank),stat=err)
    if (err /= 0) then
      call gkw_throw_abort('could not allocate integer buffer maxdims')
      return
    end if
    allocate(dims(rank),stat=err)
    if (err /= 0) then
      call gkw_throw_abort('could not allocate integer buffer dims')
      return
    end if
    allocate(dims_chunk(rank),stat=err)
    if (err /= 0) then
      call gkw_throw_abort('could not allocate integer buffer dims_chunk')
      return
    end if

    ! Create the data space with unlimited max. first dimension.  dims is
    ! a one-dimensional array of size rank specifying the current size of each
    ! dimension of the dataset. maxdims is an array of the same
    ! size specifying the upper limit on the size of each dimension.
    maxdims(1) = H5S_UNLIMITED_F
    maxdims(2:) = hyperslab_dims

    !#define METHODA 1
#define METHODB 1
#if defined(METHODA)
    ! Either (METHOD A) we set the initial extent of the dataset to a
    ! precomputed value (but in which state is the data left if GKW
    ! then exits before completing the run as expected?)

    ! One would have to compute here how much data this run will add to the
    ! dataset.
#define EXPECTED_NHYPERSLABS 1000
    dims(1) = EXPECTED_NHYPERSLABS

    ! We would also need a counter for every dataset which tells how
    ! much data was already written.
#error "METHOD A to create extendible HDF5 datasets is not implemented."

#elif defined(METHODB)
    ! or we (METHOD B) extend the dataset size by one hyperslab each
    ! time we append something to the dataset.

    ! return an integer which can be used to write data into that
    ! dataset.
    call get_free_lun(lun)
    if(gkw_its_critical()) return
    unit_is_complex(lun) = .false.

    ! check if a group with the given name exists already in the HDF5
    ! root group
    call h5lexists_f(main_output_file_id, ss(groupname), &
       & link_exists, err)
    if(link_exists) then
      ! From docs: HDF5 objects may be opened more than once at the same
      ! time. For example, an application might open a group twice,
      ! receiving two identifiers.

      ! We handle one dataset lun=12 in a group through one identifier
      ! dset_id_list(12, I_GROUP), the information from another
      ! dataset lun=34 in the group through a different
      ! identifier dest_id_list(34, I_GROUP).

      call h5gopen_f(main_output_file_id, ss(groupname), &
         & dset_id_list(lun, I_GROUP), err)
    else
      ! create a new group
      call h5gcreate_f(main_output_file_id, ss(groupname), &
         & dset_id_list(lun, I_GROUP), err)
    end if

    ! check if a dataset with the given name exists already in the
    ! group
    call h5lexists_f(dset_id_list(lun, I_GROUP), ss(dsetname), &
       & link_exists, err)
    check_err(err)
    if(link_exists) then
      call h5dopen_f(dset_id_list(lun, I_GROUP), ss(dsetname), &
         & dset_id_list(lun, I_DSET), err)
      ! FIXME Should one check here if its datatype/precision
      ! matches APPROPRIATE_HDF5_REAL_TYPE ?
      ! FIXME maybe check here if the existing dataset has indeed an infinite
      ! max. first dimension
      ! FIXME maybe check here if the existing
      ! dataset is indeed chunked
    else

      dims = maxdims
      dims(1) = 0
      call h5screate_simple_f(rank, dims, dspace_id, err, maxdims)
      check_err(err)
      ! Any dataset with an unlimited dimension must be chunked;
      ! (Generally, a dataset must be chunked if current dims does not
      ! equal maxdims)

      ! Get a dataset creation property list.
      call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, err)
      
      check_err(err)
      ! Modify dataset creation properties, i.e. enable chunking: Set
      ! the size of the chunks used to store a chunked layout dataset.
      ! The choice of the chunk size affects performance. But the
      ! library also buffers internally a bit.
      dims_chunk = maxdims
      dims_chunk(1) = 1
      call h5pset_layout_f(crp_list, H5D_CHUNKED_F, err)
      check_err(err)
      call h5pset_chunk_f(crp_list, rank, dims_chunk, err)
      check_err(err)

      if(with_szip_compression) then
        pixels_per_block = get_pixels_per_block(int(hyperslab_dims, SIZE_T))
        if(pixels_per_block > 1) then
          call h5pset_szip_f(crp_list, H5_SZIP_EC_OM_F, pixels_per_block, err)
          if(root_and_verbose) write (*,*) 'Use '//int2char(pixels_per_block) &
             & //' pixels per block to compress '//dsetname
        else
          if(root_and_verbose) write (*,*) 'Could not determine pixels&
             & per block for '//dsetname, hyperslab_dims
        end if
      end if

      ! With this property list crp_list in hand, create a dataset.
      call h5dcreate_f(dset_id_list(lun, I_GROUP), ss(dsetname), &
         & APPROPRIATE_HDF5_REAL_TYPE, dspace_id, &
         & dset_id_list(lun, I_DSET), err, crp_list)
      check_err(err)

      ! The dataspace can be closed, because later it won't be valid anymore.
      call h5sclose_f(dspace_id, err)
      check_err(err)
      call h5pclose_f(crp_list, err)
      check_err(err)
    end if

#endif
    deallocate(maxdims)
    deallocate(dims)
    deallocate(dims_chunk)

  end subroutine open_real_single_dataset_hdf5

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  function get_creation_prop_list(dims) result(crp_list)
    integer(HSIZE_T), intent(in) :: dims(:)
    integer(HID_T) :: crp_list      ! Dataset creation property list
    integer :: err
    integer :: pixels_per_block, rank
    
    ! Create the new dataset creation property list, add the szip
    ! compression filter and set a chunk size for the dataset.
    call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, err)

    ! If compression is possible, overwrite crp_list with other
    ! properties.
    if(with_szip_compression) then

      ! First examine the shape of the data array and try to find a
      ! reasonable value pixels_per_block.
      ! Rudimentarily:
      pixels_per_block = get_pixels_per_block(int(dims, SIZE_T))
      if(pixels_per_block > 1) then
        rank = size(dims)
        call h5pset_layout_f(crp_list, H5D_CHUNKED_F, err)
        call h5pset_chunk_f(crp_list, rank, dims, err)
        call h5pset_szip_f(crp_list, H5_SZIP_EC_OM_F, pixels_per_block, err)
      end if
    end if

  end function get_creation_prop_list

  !---------------------------------------------------------------------------
  !> The aim of this subroutine is to provide a means to split a
  !> dataset like 'time.dat' into several distinct datasets, each containing
  !> a column of the legacy 'time.dat'.
  !>
  !> 
  !---------------------------------------------------------------------------
  subroutine open_real_splitted_dataset_hdf5(dsetname, groupname, hyperslab_dims, &
     & lun, hyperslab_labels)
    use global, only : gkw_throw_abort, gkw_its_critical, ss

    character (len=*), intent(in) :: dsetname
    integer, dimension(:), intent(in) :: hyperslab_dims
    character (len=*), intent(in) :: groupname
    integer, intent(out) :: lun
    character (len=*), dimension(:), intent(in) :: hyperslab_labels

    integer :: i

    ! get a lun for the whole thing
    call get_free_lun(lun)
    if(gkw_its_critical()) return
    unit_is_complex(lun) = .false.
    
    ! Now we register this lun as the mother of several luns, referring to
    ! a sub-dataset, respectively.
    ! We have to associate each element of the hyperslab with a lun.

    ! The fact that splitted_luns_list(lun)%lun_sublist is allocated
    ! serves as a flag that this logical unit has several sub-units.
    if(.not.associated(splitted_luns_list(lun)%lun_sublist)) then
      allocate(splitted_luns_list(lun)%lun_sublist(hyperslab_dims(1)))
    else
      call gkw_throw_abort("Error: sub-unit list was already allocated")
      ! This should never be the case. If so, then there is a
      ! programming mistake.
    end if

    do i = 1, hyperslab_dims(1)
      if(size(hyperslab_dims) > 1) then
        call open_real_single_dataset_hdf5(hyperslab_labels(i), ss(groupname), &
           & hyperslab_dims(2:), &
           & splitted_luns_list(lun)%lun_sublist(i))
      else
        call open_real_single_dataset_hdf5(hyperslab_labels(i), ss(groupname), &
         & (/ 1 /), &
         & splitted_luns_list(lun)%lun_sublist(i))
      end if
    end do

  end subroutine open_real_splitted_dataset_hdf5

  !---------------------------------------------------------------------------
  !> This subroutine opens a complex dataset. As HDF5 does not support
  !> Fortrans intrinsic complex datatype, this subroutine at the
  !> moment actually creates 2 datasets, one for the real part and
  !> another for the imaginary part.
  !>
  !> The name of the real part dataset is dsetname//real_suffix .
  !> The name of the imaginary part dataset is dsetname//imag_suffix .
  !---------------------------------------------------------------------------
  subroutine open_complex_dataset_hdf5(dsetname, groupname, hyperslab_dims, &
     & lun)
    use global, only : ss
    character (len=*), intent(in) :: dsetname
    integer, dimension(:), intent(in) :: hyperslab_dims
    character (len=*), intent(in) :: groupname
    integer, intent(out) :: lun
    integer :: real_lun, imag_lun
    if(is_not_root_proc()) return

    call open_real_dataset_hdf5(trim(ss(dsetname))//real_suffix, &
       & ss(groupname), hyperslab_dims, real_lun)
    call open_real_dataset_hdf5(trim(ss(dsetname))//imag_suffix, &
       & ss(groupname), hyperslab_dims, imag_lun)
    lun = real_lun

    ! set an entry which can be tested when closing the dataset
    unit_is_complex(lun) = .true.

    ! there should never be the case that this value is queried by
    ! anyone, but set it nevertheless:
    unit_is_complex(lun+1) = .false.

  end subroutine open_complex_dataset_hdf5

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine get_free_lun(lun)
    use global, only : gkw_throw_abort
    integer, intent(out) :: lun
    if(n_open_dsets < max_dsets - 2) then
      n_open_dsets = n_open_dsets + 1
      ! return an integer which can be used to write data into that
      ! dataset.
      lun = n_open_dsets
    else
      ! In this case one has to increase max_dsets
      call gkw_throw_abort("Not enough HDF5 logical units available.")
      return
    end if
  end subroutine get_free_lun

  !---------------------------------------------------------------------------
  !> This subroutine closes a dataset.
  !---------------------------------------------------------------------------
  subroutine close_dataset_hdf5(lun)
    integer, intent(in) :: lun
    integer :: i
    if(is_not_root_proc()) return

    if(unit_is_complex(lun)) then
      call close_complex_dataset_hdf5(lun)
    else
      ! check if this lu has sub-units
      if(associated(splitted_luns_list(lun)%lun_sublist)) then
        ! This lu has sub-units. Close them all.
        do i = 1, size(splitted_luns_list(lun)%lun_sublist)
          call close_real_dataset_hdf5(splitted_luns_list(lun)%lun_sublist(i))
        end do
        ! Then deallocate the list.
        deallocate(splitted_luns_list(lun)%lun_sublist)
      else
        ! This lu does not have sub-units.
        call close_real_dataset_hdf5(lun)
      end if
    end if
  end subroutine close_dataset_hdf5

  !---------------------------------------------------------------------------
  !> This subroutine closes a real dataset.
  !---------------------------------------------------------------------------
  subroutine close_real_dataset_hdf5(lun)
    use global, only : gkw_throw_warn, int2char
    integer, intent(in) :: lun
    ! error flag
    integer :: err
    logical :: is_valid
    call h5iis_valid_f(dset_id_list(lun, I_DSET), is_valid, err)
    if(is_valid) then
      ! End access to the dataset and release resources used by it.
      call h5dclose_f(dset_id_list(lun, I_DSET), err)
    else
      call gkw_throw_warn("An HDF5 dataset is closed twice, lun="//int2char(lun))
    end if

    call h5iis_valid_f(dset_id_list(lun, I_GROUP), is_valid, err)
    if(is_valid) then
      ! End access to the group which contains this dataset.  Note that
      ! other identifiers associated to this group may continue to be
      ! open.
      call h5gclose_f(dset_id_list(lun, I_GROUP), err)
    end if
  end subroutine close_real_dataset_hdf5

  !---------------------------------------------------------------------------
  !> This subroutine closes a complex dataset.
  !---------------------------------------------------------------------------
  subroutine close_complex_dataset_hdf5(lun)
    integer, intent(in) :: lun

    if(is_not_root_proc()) return

    ! End access to the dataset with the real part.
    call close_real_dataset_hdf5(lun)
    ! Close also the associated dataset with the imaginary part.
    call close_real_dataset_hdf5(lun+1)
  end subroutine close_complex_dataset_hdf5

  !---------------------------------------------------------------------------
  !> This subroutine causes buffered data to be flushed to the
  !> physical file which contains the given logical unit (i.e. dataset in the
  !> sense of this module). Note that HDF5 does not have complete control
  !> about the flushing at low level.
  !---------------------------------------------------------------------------
  subroutine flush_file(lun)
    integer, intent(in) :: lun
    ! error flag
    integer :: err
    integer(HID_T), save :: file_id

    if(is_not_root_proc()) return

    call h5fflush_f(main_output_file_id, H5F_SCOPE_GLOBAL_F, err)
    check_err(err)

    ! just in case not everything is stored into a single hdf5 file:
    call h5iget_file_id_f(dset_id_list(lun,I_DSET), file_id, err)
    check_err(err)
    call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, err)
    check_err(err)

  end subroutine flush_file

  !---------------------------------------------------------------------------
  !> This subroutine causes buffered data to be flushed to the
  !> physical file. Note that HDF5 does not have complete control
  !> about the flushing at low level.
  !---------------------------------------------------------------------------
  subroutine flush_all_open_files()
    ! error flag
    integer :: err
    if(is_not_root_proc()) return

    call h5fflush_f(main_output_file_id, H5F_SCOPE_GLOBAL_F, err)
    check_err(err)

  end subroutine flush_all_open_files


  !---------------------------------------------------------------------------
  !> This is a place to put a more or less clever algorithm to find a
  !> block size for SZIP compression.
  !---------------------------------------------------------------------------
  function get_pixels_per_block(dims)
    integer(SIZE_T), dimension(:), intent(in) :: dims
    integer :: get_pixels_per_block

    ! block size is passed in the parameter pixels_per_block and must
    ! be even and not greater than 32, with typical values being 8,
    ! 10, 16, or 32.
    integer, dimension(7), parameter :: block_sizes = &
       & (/ 8, 10, 12, 16, 20, 24, 32 /)
    integer :: i

    ! To achieve optimal performance for SZIP compression, it is
    ! recommended that a chunk's fastest-changing dimension be
    ! equal to N times pixels_per_block where N is the maximum
    ! number of blocks per scan line allowed by the SZIP
    ! library. In the current version of SZIP, N is set to 128.

    get_pixels_per_block = -1

    ! At the moment this is really rudimentary (does not even run
    ! over all elements)
    do i = size(block_sizes), 1, -1
      if(block_sizes(i) <= int(dims(1)) .and. &
         & mod(int(dims(1)), block_sizes(i)) == 0) then
        get_pixels_per_block = block_sizes(i)
        return
      end if
    end do
    do i = size(block_sizes), 1, -1
      if(block_sizes(i) <= int(dims(1))) then
        get_pixels_per_block = block_sizes(i)
        return
      end if
    end do

    
  end function get_pixels_per_block


  !---------------------------------------------------------------------------
  !> output a real 1-D array data(:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_1d_real(dsetname, groupname, data, order)
    real, dimension(:), intent(in)        :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order

    integer, parameter :: rank = 1
    integer :: err
    logical :: link_exists
    integer(HID_T) :: group_id      ! Group identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier
    integer(HID_T) :: crp_list      ! Dataset creation property list
    integer(HSIZE_T), dimension(rank) :: dims ! Dataset dimensions
    if(is_not_root_proc()) return

    ! To keep the compiler quiet.
    if (present(order)) continue

    ! check if the desired group is there. 
    call h5lexists_f(main_output_file_id, ss(groupname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    else
      ! create the group
      call h5gcreate_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    end if

    ! check if the dataset is already there
    call h5lexists_f(group_id, ss(dsetname), link_exists, err)
    check_err(err)
    if(link_exists) then
      call h5dopen_f(group_id, ss(dsetname), dset_id, err)
      check_err(err)
    else
      ! Create the dataspace.
      dims = shape(data)
      call h5screate_simple_f(rank, dims, dspace_id, err)

      crp_list = get_creation_prop_list(dims)

      ! Create the dataset.
      call h5dcreate_f(group_id, ss(dsetname), &
         & APPROPRIATE_HDF5_REAL_TYPE, dspace_id, &
         & dset_id, err, crp_list)
      ! The data space object and the propertly list are not needed any more.
      call h5sclose_f(dspace_id, err)
      call h5pclose_f(crp_list, err)
      
    end if

    
    ! Write the data
    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, data, dims, err)

    ! End access
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)

  end subroutine output_array_1d_real

  !---------------------------------------------------------------------------
  !> output a real 2-D array data(:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_2d_real(dsetname, groupname, data, order)
    real, dimension(:,:), intent(in)        :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order

    integer, parameter :: rank = 2
    integer :: err
    logical :: link_exists
    integer(HID_T) :: group_id      ! Group identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier
    integer(HID_T) :: crp_list      ! Dataset creation property list
    integer(HSIZE_T), dimension(rank) :: dims ! Dataset dimensions
    if(is_not_root_proc()) return

    ! To keep the compiler quiet.
    if (present(order)) continue

    ! check if the desired group is there. 
    call h5lexists_f(main_output_file_id, ss(groupname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    else
      ! create the group
      call h5gcreate_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    end if

    call h5lexists_f(group_id, ss(dsetname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! Open the dataset
      call h5dopen_f(group_id, ss(dsetname), dset_id, err)
      check_err(err)

    else
      ! Create the dataspace.
      dims = shape(data)
      call h5screate_simple_f(rank, dims, dspace_id, err)

      crp_list = get_creation_prop_list(dims)

      ! Create the dataset.
      call h5dcreate_f(group_id, ss(dsetname), &
         & APPROPRIATE_HDF5_REAL_TYPE, dspace_id, &
         & dset_id, err, crp_list)

      ! The data space object and the propertly list are not needed any more.
      call h5sclose_f(dspace_id, err)
      call h5pclose_f(crp_list, err)

    end if

    ! Write the data
    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, data, dims, err)

    ! End access
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)

  end subroutine output_array_2d_real

  !---------------------------------------------------------------------------
  !> output a real 3-D array data(:,:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_3d_real(dsetname, groupname, data, order)
    real, dimension(:,:,:), intent(in)        :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order

    integer, parameter :: rank = 3
    integer :: err
    logical :: link_exists
    integer(HID_T) :: group_id      ! Group identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier
    integer(HID_T) :: crp_list      ! Dataset creation property list
    integer(HSIZE_T), dimension(rank) :: dims ! Dataset dimensions
    if(is_not_root_proc()) return

    ! To keep the compiler quiet.
    if (present(order)) continue

    ! check if the desired group is there. 
    call h5lexists_f(main_output_file_id, ss(groupname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    else
      ! create the group
      call h5gcreate_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    end if

    ! check if the dataset is already there
    call h5lexists_f(group_id, ss(dsetname), link_exists, err)
    check_err(err)
    if(link_exists) then
      call h5dopen_f(group_id, ss(dsetname), dset_id, err)
      check_err(err)
    else
      ! Create the dataspace.
      dims = shape(data)
      call h5screate_simple_f(rank, dims, dspace_id, err)

      crp_list = get_creation_prop_list(dims)

      ! Create the dataset.
      call h5dcreate_f(group_id, ss(dsetname), &
         & APPROPRIATE_HDF5_REAL_TYPE, dspace_id, &
         & dset_id, err, crp_list)

      ! The data space object and the propertly list are not needed any more.
      call h5sclose_f(dspace_id, err)
      call h5pclose_f(crp_list, err)
    end if

    ! Write the data
    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, data, dims, err)

    ! End access
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)

  end subroutine output_array_3d_real


  !---------------------------------------------------------------------------
  !> output a real 4D array data(:,:,:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_4d_real(dsetname, groupname, data, order)
    real, dimension(:,:,:,:), intent(in)        :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order

    integer, parameter :: rank = 4
    integer :: err
    logical :: link_exists
    integer :: pixels_per_block
    integer(HID_T) :: group_id      ! Group identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier
    integer(HID_T) :: crp_list      ! Dataset creation property list
    integer(HSIZE_T), dimension(rank) :: dims ! Dataset dimensions
    if(is_not_root_proc()) return

    ! To keep the compiler quiet.
    if (present(order)) continue

    ! check if the desired group is there. 
    call h5lexists_f(main_output_file_id, ss(groupname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    else
      ! create the group
      call h5gcreate_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    end if

    ! check if the dataset is already there
    call h5lexists_f(group_id, ss(dsetname), link_exists, err)
    check_err(err)
    if(link_exists) then
      call h5dopen_f(group_id, ss(dsetname), dset_id, err)
      check_err(err)
    else
      ! Create the dataspace.
      dims = shape(data)
      call h5screate_simple_f(rank, dims, dspace_id, err)

      ! A priori, use the default properties.
      call h5pcopy_f(H5P_DEFAULT_F, crp_list, err)
      ! If compression is possible, overwrite crp_list with other
      ! properties.
      if(with_szip_compression) then
        ! Rudimentarily:
        if(dims(1) > MIN_PIXELS_PER_BLOCK) then
          pixels_per_block = get_pixels_per_block(int(dims, SIZE_T))
          ! Close the default property list.
          call h5pclose_f(crp_list, err)
          ! Create the new dataset creation property list, add the szip
          ! compression filter and set a chunk size for the dataset.
          call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, err)
          call h5pset_layout_f(crp_list, H5D_CHUNKED_F, err)
          call h5pset_chunk_f(crp_list, rank, dims, err)
          call h5pset_szip_f(crp_list, H5_SZIP_EC_OM_F, pixels_per_block, err)
        end if
      end if

      ! Create the dataset.
      call h5dcreate_f(group_id, ss(dsetname), &
         & APPROPRIATE_HDF5_REAL_TYPE, dspace_id, &
         & dset_id, err, crp_list)

      ! The data space object and the propertly list are not needed any more.
      call h5sclose_f(dspace_id, err)
      call h5pclose_f(crp_list, err)
    end if

    ! Write the data
    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, data, dims, err)

    ! End access
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)

  end subroutine output_array_4d_real


  !---------------------------------------------------------------------------
  !> output a real 5D array data(:,:,:,:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_5d_real(dsetname, groupname, data, order)
    real, dimension(:,:,:,:,:), intent(in)        :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order

    integer, parameter :: rank = 5
    integer :: err
    logical :: link_exists
    integer :: pixels_per_block
    integer(HID_T) :: group_id      ! Group identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier
    integer(HID_T) :: crp_list      ! Dataset creation property list
    integer(HSIZE_T), dimension(rank) :: dims ! Dataset dimensions
    if(is_not_root_proc()) return

    ! To keep the compiler quiet.
    if (present(order)) continue

    ! check if the desired group is there. 
    call h5lexists_f(main_output_file_id, ss(groupname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    else
      ! create the group
      call h5gcreate_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    end if

    ! check if the dataset is already there
    call h5lexists_f(group_id, ss(dsetname), link_exists, err)
    check_err(err)
    if(link_exists) then
      call h5dopen_f(group_id, ss(dsetname), dset_id, err)
      check_err(err)
    else
      ! Create the dataspace.
      dims = shape(data)
      call h5screate_simple_f(rank, dims, dspace_id, err)

      ! A priori, use the default properties.
      call h5pcopy_f(H5P_DEFAULT_F, crp_list, err)
      ! If compression is possible, overwrite crp_list with other
      ! properties.
      if(with_szip_compression) then
        ! Rudimentarily:
        if(dims(1) > MIN_PIXELS_PER_BLOCK) then
          pixels_per_block = get_pixels_per_block(int(dims, SIZE_T))
          ! Close the default property list.
          call h5pclose_f(crp_list, err)
          ! Create the new dataset creation property list, add the szip
          ! compression filter and set a chunk size for the dataset.
          call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, err)
          call h5pset_layout_f(crp_list, H5D_CHUNKED_F, err)
          call h5pset_chunk_f(crp_list, rank, dims, err)
          call h5pset_szip_f(crp_list, H5_SZIP_EC_OM_F, pixels_per_block, err)
        end if
      end if

      ! Create the dataset.
      call h5dcreate_f(group_id, ss(dsetname), &
         & APPROPRIATE_HDF5_REAL_TYPE, dspace_id, &
         & dset_id, err, crp_list)

      ! The data space object and the propertly list are not needed any more.
      call h5sclose_f(dspace_id, err)
      call h5pclose_f(crp_list, err)
    end if

    ! Write the data
    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, data, dims, err)

    ! End access
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)

  end subroutine output_array_5d_real

  !---------------------------------------------------------------------------
  !> output a real 6D array data(:,:,:,:,:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_6d_real(dsetname, groupname, data, order)
    real, dimension(:,:,:,:,:,:), intent(in)        :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order

    integer, parameter :: rank = 6
    integer :: err
    logical :: link_exists
    integer :: pixels_per_block
    integer(HID_T) :: group_id      ! Group identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier
    integer(HID_T) :: crp_list      ! Dataset creation property list
    integer(HSIZE_T), dimension(rank) :: dims ! Dataset dimensions
    if(is_not_root_proc()) return

    ! To keep the compiler quiet.
    if (present(order)) continue

    ! check if the desired group is there. 
    call h5lexists_f(main_output_file_id, ss(groupname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    else
      ! create the group
      call h5gcreate_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    end if

    ! check if the dataset is already there
    call h5lexists_f(group_id, ss(dsetname), link_exists, err)
    check_err(err)
    if(link_exists) then
      call h5dopen_f(group_id, ss(dsetname), dset_id, err)
      check_err(err)
    else
      ! Create the dataspace.
      dims = shape(data)
      call h5screate_simple_f(rank, dims, dspace_id, err)

      ! A priori, use the default properties.
      call h5pcopy_f(H5P_DEFAULT_F, crp_list, err)
      ! If compression is possible, overwrite crp_list with other
      ! properties.
      if(with_szip_compression) then
        ! Rudimentarily:
        if(dims(1) > MIN_PIXELS_PER_BLOCK) then
          pixels_per_block = get_pixels_per_block(int(dims, SIZE_T))
          ! Close the default property list.
          call h5pclose_f(crp_list, err)
          ! Create the new dataset creation property list, add the szip
          ! compression filter and set a chunk size for the dataset.
          call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, err)
          call h5pset_layout_f(crp_list, H5D_CHUNKED_F, err)
          call h5pset_chunk_f(crp_list, rank, dims, err)
          call h5pset_szip_f(crp_list, H5_SZIP_EC_OM_F, pixels_per_block, err)
        end if
      end if

      ! Create the dataset.
      call h5dcreate_f(group_id, ss(dsetname), &
         & APPROPRIATE_HDF5_REAL_TYPE, dspace_id, &
         & dset_id, err, crp_list)

      ! The data space object and the propertly list are not needed any more.
      call h5sclose_f(dspace_id, err)
      call h5pclose_f(crp_list, err)
    end if

    ! Write the data
    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, data, dims, err)

    ! End access
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)

  end subroutine output_array_6d_real

  !---------------------------------------------------------------------------
  !> output a complex 1-D array data(:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_1d_complex(dsetname, groupname, data, order)
    complex, dimension(:), intent(in)     :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order
    if(is_not_root_proc()) return

    ! From HDF5 docs: "The HDF5 Fortran Library before version 1.8.8
    ! could handle only INTEGER, REAL, and CHARACTER types, and an
    ! obsolete DOUBLE PRECISION type. It could not support COMPLEX and
    ! LOGICAL types because there is no support for the corresponding
    ! C types in the HDF5 C Library."  New Fortran2003 based features
    ! in HDF5 allow more datatypes to be written and read.

    ! For now, I seems that we have to store real and imaginary part
    ! separately, into one real dataset respectively.
    call output_array_1d_real(trim(ss(dsetname))//real_suffix, ss(groupname), &
       & real(data), order)
    call output_array_1d_real(trim(ss(dsetname))//imag_suffix, ss(groupname), &
       & aimag(data), order)

  end subroutine output_array_1d_complex

  !---------------------------------------------------------------------------
  !> output a complex 2-D array data(:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_2d_complex(dsetname, groupname, data, order)
    complex, dimension(:,:), intent(in)     :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order
    if(is_not_root_proc()) return

    call output_array_2d_real(trim(ss(dsetname))//real_suffix, ss(groupname), &
       & real(data), order)
    call output_array_2d_real(trim(ss(dsetname))//imag_suffix, ss(groupname), &
       & aimag(data), order)

  end subroutine output_array_2d_complex

  !---------------------------------------------------------------------------
  !> output a complex 3-D array data(:,:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_3d_complex(dsetname, groupname, data, order)
    complex, dimension(:,:,:), intent(in)     :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order
    if(is_not_root_proc()) return

    call output_array_3d_real(trim(ss(dsetname))//real_suffix, ss(groupname), &
       & real(data), order)
    call output_array_3d_real(trim(ss(dsetname))//imag_suffix, ss(groupname), &
       & aimag(data), order)

  end subroutine output_array_3d_complex

  !---------------------------------------------------------------------------
  !> output a complex 6-D array data(:,:,:,:,:,:) to file
  !---------------------------------------------------------------------------
  subroutine output_array_6d_complex(dsetname, groupname, data, order)
    complex, dimension(:,:,:,:,:,:), intent(in)     :: data
    character (len=*), intent(in)           :: dsetname, groupname
    character (len=1), optional, intent(in) :: order
    if(is_not_root_proc()) return

    call output_array_6d_real(trim(ss(dsetname))//real_suffix, ss(groupname), &
       & real(data), order)
    call output_array_6d_real(trim(ss(dsetname))//imag_suffix, ss(groupname), &
       & aimag(data), order)

  end subroutine output_array_6d_complex

  !---------------------------------------------------------------------------
  !> This subroutine appends a 1D line of data data(:) to a dataset.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_dataset_hdf5().
  !>
  !---------------------------------------------------------------------------
  recursive subroutine append_chunk_1d_real(lun, data)
    use global, only : gkw_throw_abort
    use global, only : int2char

    real, dimension(:), intent(in) :: data
    integer, intent(in) :: lun

    !flag to check operation success
    integer :: err

    integer, parameter :: rank = 2
    integer :: r
    !integer(HSIZE_T) :: npoints ! just for curiosity
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier 
    integer(HID_T) :: mem_dspace_id ! Memory dataspace identifier 
    integer(HID_T) :: crp_list      ! Dataset creation property identifier 

    integer(HSIZE_T), dimension(rank) :: dims, maxdims
    integer(HSIZE_T), dimension(rank) :: dims_chunk
    integer(HSIZE_T), dimension(rank) :: newdims
    integer(HSIZE_T), dimension(rank) :: offset
    integer(HSIZE_T), dimension(rank), parameter :: ones = 1
    if(is_not_root_proc()) return

    ! check if the lun refers to a splitted logical unit
    if(associated(splitted_luns_list(lun)%lun_sublist)) then
      do r = 1, size(data, 1)
        call append_chunk(splitted_luns_list(lun)%lun_sublist(r), (/ data(r) /))
      end do
      return
    end if

    dset_id = dset_id_list(lun, I_DSET)
    call h5iget_name_f(dset_id, dset_name_buf, dset_name_buf_size, &
       & name_size, err)

#if defined(METHODA)
    ! (METHOD A) If necessary, extend the dataset.
    ! This could be necessary, e.g. in case of a restarted run.

    !if (.false.) then
    !  newdims(1) = newdims + 500
    !  call h5dset_extent_f(dset_id, newdims, err)
    !end if
#error "METHOD A for extendible HDF5 datasets is not implemented."

#elif defined(METHODB)
    ! (METHOD B) extend the dataset each time we append a chunk of data.

    ! Get the dataspace.
    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)

    ! Get some info, do some checks.
    call h5sget_simple_extent_ndims_f(dspace_id, r, err)
    check_err(err)
    if(r /= rank) then
      call gkw_throw_abort("The queried HDF5 dataset "// &
         & dset_name_buf(1:name_size)//" has rank "//int2char(r,2) &
         & //" but "//int2char(rank,2)//" is needed.")
      return
    end if

    ! call h5sget_simple_extent_npoints_f(dspace_id, npoints, err)
    ! FIXME FIXME When h5sget_simple_extent_npoints_f is called first, this
    ! routine returns the correct npoints==0 but err==-1 !
    ! Is this a flaw in the library or do I misunderstand something?
    ! check_err(err)
    ! write (*,*) "npoints = ", npoints

    ! Find the current extent of the dataset, i.e. out how much data is
    ! already there.
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    check_err(err)
    if(any(dims(2:) /= shape(data))) then
      write (*,*) "hyperslab shape of the queried HDF5 dataset:", dims(2:)
      write (*,*) "data shape:", shape(data)
      ! FIXME Is it necessary and practicable to check if dims(1,:) has
      ! the same shape as data(:)?
      call gkw_throw_abort("Wrong data shape for "//dset_name_buf(1:name_size))
      return
    end if

    ! Find the chunk size associated to this dataset. Of course, this
    ! should match the shape of the data
    call h5dget_create_plist_f(dset_id, crp_list, err)
    check_err(err)
    call h5pget_chunk_f(crp_list, rank, dims_chunk, err)
    check_err(err)
    ! Note that for the chunks size(dims_chunk,1)==1
    if(any(dims_chunk(2:) /= shape(data))) then
      write (*,*) "dims_chunk = ", dims_chunk
      write (*,*) "shape(data) = ", shape(data)
      call gkw_throw_abort("HDF5 chunk size does not match data shape.")
      return
    end if
    call h5pclose_f(crp_list, err)
    check_err(err)

    ! We have set the initial size of the dataset to 0 x
    ! shape(hyperslab) (e.g. 0 x ncolumns), thus we need to extend it
    ! by one hyperslab before we write.  Note that we extend the
    ! dataset itself, not its dataspace.
    newdims = dims
    newdims(1) = dims(1) + 1
    call h5dset_extent_f(dset_id, newdims, err)
    check_err(err)

    ! Create a dataspace mem_dspace_id which indicates the size of our
    ! buffer in memory. (Actually, data does have 1 dimension less, but
    ! there is a reshape(data) below)
    call h5screate_simple_f (rank, dims_chunk, mem_dspace_id, err)
    check_err(err)

    ! We now select the hyperslab of the dataspace (e.g a line-shaped
    ! part in a 2D dataspace) describing the part of the dataset we
    ! want to write into. For this, we first get the dataspace.
    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)
    ! Replace the existing dataspace selection with a new
    ! selection.
    offset = 0
    offset(1) = newdims(1)-1

    ! Select a row in that dataspace.
    ! call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
    !    & offset, dims_chunk, err, (/1,1/), (/1,1/))
    ! If SRG understands h5sselect_hyperslab_f correctly, this
    ! call selects a total of dims_chunk(1)*dims_chunk(2) blocks of 1 element
    ! each, without stride. This would be equivalent to writing
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
       & offset, ones, err, ones, dims_chunk)
    check_err(err)
    ! which selects 1 block with dims_chunk(1)*dims_chunk(2) elements.

    ! Just for curiosity: Now mem_space_id and dspace_id should have
    ! the same number of elements selected.
    ! FIXME FIXME Those information calls would throw errors! Is this a flaw
    ! in the library or do I misunderstand something?
    ! call H5Sget_select_elem_npoints_f(dspace_id, npoints, err)
    ! check_err(err)
    ! write (*,*) "dspace npoints = ", npoints
    ! call H5Sget_select_elem_npoints_f(mem_dspace_id, npoints, err)
    ! check_err(err)
    ! write (*,*) "memspace npoints = ", npoints

    ! Write a hyperslab of data. 
    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, &
       & reshape(data, dims_chunk), &
       & dims_chunk, err, &
       & mem_dspace_id, dspace_id)
    check_err(err)

    call h5sclose_f(dspace_id, err)
    ! FIXME We could also reuse the memory dataspace instead of
    ! creating a new one each time.
    call h5sclose_f(mem_dspace_id, err)
#endif

  end subroutine append_chunk_1d_real

  !---------------------------------------------------------------------------
  !> This subroutine appends a 2D matrix of data data(:,:) to a dataset.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_dataset_hdf5().
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_2d_real(lun, data)
    use global, only : gkw_throw_abort
    use global, only : int2char

    real, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun

    !flag to check operation success
    integer :: err

    integer, parameter :: rank = 3
    integer :: r
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier 
    integer(HID_T) :: mem_dspace_id ! Memory dataspace identifier 
    integer(HID_T) :: crp_list      ! Dataset creation property identifier 

    integer(HSIZE_T), dimension(rank) :: dims, maxdims
    integer(HSIZE_T), dimension(rank) :: dims_chunk
    integer(HSIZE_T), dimension(rank) :: newdims
    integer(HSIZE_T), dimension(rank) :: offset
    integer(HSIZE_T), dimension(rank), parameter :: ones = 1
    if(is_not_root_proc()) return

    ! check if the lun refers to a splitted logical unit
    if(associated(splitted_luns_list(lun)%lun_sublist)) then
      do r = 1, size(data, 1)
        call append_chunk(splitted_luns_list(lun)%lun_sublist(r), data(r, :))
      end do
      return
    end if

    dset_id = dset_id_list(lun, I_DSET)
    call h5iget_name_f(dset_id, dset_name_buf, dset_name_buf_size, name_size, err)

#if defined(METHODA)
#error "METHOD A for extendible HDF5 datasets is not implemented."
#elif defined(METHODB)
    ! (METHOD B) extend the dataset each time we append a chunk of data.

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)

    call h5sget_simple_extent_ndims_f(dspace_id, r, err)
    check_err(err)
    if(r /= rank) then
      call gkw_throw_abort("The queried HDF5 dataset "// &
         & dset_name_buf(1:name_size)//" has rank "//int2char(r,2) &
         & //" but "//int2char(rank,2)//" is needed.")
      return
    end if

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    check_err(err)
    if(any(dims(2:) /= shape(data))) then
      write (*,*) "hyperslab shape of the queried HDF5 dataset:", dims(2:)
      write (*,*) "data shape:", shape(data)
      ! FIXME Is it necessary and practicable to check if dims(1,:) has
      ! the same shape as data(:)?
      call gkw_throw_abort("Wrong data shape for "//dset_name_buf(1:name_size))
      return
    end if

    call h5dget_create_plist_f(dset_id, crp_list, err)
    check_err(err)
    call h5pget_chunk_f(crp_list, rank, dims_chunk, err)
    check_err(err)
    if(any(dims_chunk(2:) /= shape(data))) then
      write (*,*) "dims_chunk = ", dims_chunk
      write (*,*) "shape(data) = ", shape(data)
      call gkw_throw_abort("HDF5 chunk size does not match data shape.")
      return
    end if
    call h5pclose_f(crp_list, err)
    check_err(err)

    newdims = dims
    newdims(1) = dims(1) + 1
    call h5dset_extent_f(dset_id, newdims, err)
    check_err(err)

    call h5screate_simple_f (rank, dims_chunk, mem_dspace_id, err)
    check_err(err)

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)
    offset = 0
    offset(1) = newdims(1)-1

    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
       & offset, ones, err, ones, dims_chunk)
    check_err(err)

    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, &
       & reshape(data, dims_chunk), &
       & dims_chunk, err, &
       & mem_dspace_id, dspace_id)
    check_err(err)

    call h5sclose_f(dspace_id, err)
    call h5sclose_f(mem_dspace_id, err)
#endif

  end subroutine append_chunk_2d_real

  !---------------------------------------------------------------------------
  !> This subroutine appends a 3D matrix of data data(:,:,:) to a dataset.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_dataset_hdf5().
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_3d_real(lun, data)
    use global, only : gkw_throw_abort
    use global, only : int2char

    real, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    !flag to check operation success
    integer :: err

    integer, parameter :: rank = 4
    integer :: r
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier 
    integer(HID_T) :: mem_dspace_id ! Memory dataspace identifier 
    integer(HID_T) :: crp_list      ! Dataset creation property identifier 

    integer(HSIZE_T), dimension(rank) :: dims, maxdims
    integer(HSIZE_T), dimension(rank) :: dims_chunk
    integer(HSIZE_T), dimension(rank) :: newdims
    integer(HSIZE_T), dimension(rank) :: offset
    integer(HSIZE_T), dimension(rank), parameter :: ones = 1
    if(is_not_root_proc()) return

    ! check if the lun refers to a splitted logical unit
    if(associated(splitted_luns_list(lun)%lun_sublist)) then
      do r = 1, size(data,1)
        call append_chunk(splitted_luns_list(lun)%lun_sublist(r), data(r, :, :))
      end do
      return
    end if

    dset_id = dset_id_list(lun, I_DSET)
    call h5iget_name_f(dset_id, dset_name_buf, dset_name_buf_size, name_size, err)

#if defined(METHODA)
#error "METHOD A for extendible HDF5 datasets is not implemented."
#elif defined(METHODB)
    ! (METHOD B) extend the dataset each time we append a chunk of data.

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)

    call h5sget_simple_extent_ndims_f(dspace_id, r, err)
    check_err(err)
    if(r /= rank) then
      call gkw_throw_abort("The queried HDF5 dataset "// &
         & dset_name_buf(1:name_size)//" has rank "//int2char(r,2) &
         & //" but "//int2char(rank,2)//" is needed.")
      return
    end if

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    check_err(err)
    if(any(dims(2:) /= shape(data))) then
      write (*,*) "hyperslab shape of the queried HDF5 dataset:", dims(2:)
      write (*,*) "data shape:", shape(data)
      ! FIXME Is it necessary and practicable to check if dims(1,:) has
      ! the same shape as data(:)?
      call gkw_throw_abort("Wrong data shape for "//dset_name_buf(1:name_size))
      return
    end if

    call h5dget_create_plist_f(dset_id, crp_list, err)
    check_err(err)
    call h5pget_chunk_f(crp_list, rank, dims_chunk, err)
    check_err(err)
    if(any(dims_chunk(2:) /= shape(data))) then
      write (*,*) "dims_chunk = ", dims_chunk
      write (*,*) "shape(data) = ", shape(data)
      call gkw_throw_abort("HDF5 chunk size does not match data shape.")
      return
    end if
    call h5pclose_f(crp_list, err)
    check_err(err)

    newdims = dims
    newdims(1) = dims(1) + 1
    call h5dset_extent_f(dset_id, newdims, err)
    check_err(err)

    call h5screate_simple_f (rank, dims_chunk, mem_dspace_id, err)
    check_err(err)

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)
    offset = 0
    offset(1) = newdims(1)-1

    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
       & offset, ones, err, ones, dims_chunk)
    check_err(err)

    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, &
       & reshape(data, dims_chunk), &
       & dims_chunk, err, &
       & mem_dspace_id, dspace_id)
    check_err(err)

    call h5sclose_f(dspace_id, err)
    call h5sclose_f(mem_dspace_id, err)
#endif

  end subroutine append_chunk_3d_real

  !---------------------------------------------------------------------------
  !> This subroutine appends a 4D matrix of data data(:,:,:,:) to a dataset.
  !> The first parameter is a lun. The caller obtains a
  !> valid lun by using open_real_dataset_hdf5().
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_4d_real(lun, data)
    use global, only : gkw_throw_abort
    use global, only : int2char

    real, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    !flag to check operation success
    integer :: err

    integer, parameter :: rank = 5
    integer :: r
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier 
    integer(HID_T) :: mem_dspace_id ! Memory dataspace identifier 
    integer(HID_T) :: crp_list      ! Dataset creation property identifier 

    integer(HSIZE_T), dimension(rank) :: dims, maxdims
    integer(HSIZE_T), dimension(rank) :: dims_chunk
    integer(HSIZE_T), dimension(rank) :: newdims
    integer(HSIZE_T), dimension(rank) :: offset
    integer(HSIZE_T), dimension(rank), parameter :: ones = 1
    if(is_not_root_proc()) return

    ! check if the lun refers to a splitted logical unit
    if(associated(splitted_luns_list(lun)%lun_sublist)) then
      do r = 1, size(data,1)
        call append_chunk(splitted_luns_list(lun)%lun_sublist(r), &
           & data(r, :, :, :))
      end do
      return
    end if

    dset_id = dset_id_list(lun, I_DSET)
    call h5iget_name_f(dset_id, dset_name_buf, dset_name_buf_size, name_size, err)

#if defined(METHODA)
#error "METHOD A for extendible HDF5 datasets is not implemented."
#elif defined(METHODB)
    ! (METHOD B) extend the dataset each time we append a chunk of data.

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)

    call h5sget_simple_extent_ndims_f(dspace_id, r, err)
    check_err(err)
    if(r /= rank) then
      call gkw_throw_abort("The queried HDF5 dataset "// &
         & dset_name_buf(1:name_size)//" has rank "//int2char(r,2) &
         & //" but "//int2char(rank,2)//" is needed.")
      return
    end if

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    check_err(err)
    if(any(dims(2:) /= shape(data))) then
      write (*,*) "hyperslab shape of the queried HDF5 dataset:", dims(2:)
      write (*,*) "data shape:", shape(data)
      ! FIXME Is it necessary and practicable to check if dims(1,:) has
      ! the same shape as data(:)?
      call gkw_throw_abort("Wrong data shape for "//dset_name_buf(1:name_size))
      return
    end if

    call h5dget_create_plist_f(dset_id, crp_list, err)
    check_err(err)
    call h5pget_chunk_f(crp_list, rank, dims_chunk, err)
    check_err(err)
    if(any(dims_chunk(2:) /= shape(data))) then
      write (*,*) "dims_chunk = ", dims_chunk
      write (*,*) "shape(data) = ", shape(data)
      call gkw_throw_abort("HDF5 chunk size does not match data shape.")
      return
    end if
    call h5pclose_f(crp_list, err)
    check_err(err)

    newdims = dims
    newdims(1) = dims(1) + 1
    call h5dset_extent_f(dset_id, newdims, err)
    check_err(err)

    call h5screate_simple_f (rank, dims_chunk, mem_dspace_id, err)
    check_err(err)

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)
    offset = 0
    offset(1) = newdims(1)-1

    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
       & offset, ones, err, ones, dims_chunk)
    check_err(err)

    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, &
       & reshape(data, dims_chunk), &
       & dims_chunk, err, &
       & mem_dspace_id, dspace_id)
    check_err(err)

    call h5sclose_f(dspace_id, err)
    call h5sclose_f(mem_dspace_id, err)
#endif

  end subroutine append_chunk_4d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_6d_real(lun, data)
    use global, only : gkw_throw_abort
    use global, only : int2char

    real, dimension(:,:,:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun

    !flag to check operation success
    integer :: err

    integer, parameter :: rank = 7
    integer :: r
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier 
    integer(HID_T) :: mem_dspace_id ! Memory dataspace identifier 
    integer(HID_T) :: crp_list      ! Dataset creation property identifier 

    integer(HSIZE_T), dimension(rank) :: dims, maxdims
    integer(HSIZE_T), dimension(rank) :: dims_chunk
    integer(HSIZE_T), dimension(rank) :: newdims
    integer(HSIZE_T), dimension(rank) :: offset
    integer(HSIZE_T), dimension(rank), parameter :: ones = 1
    if(is_not_root_proc()) return

    ! ! check if the lun refers to a splitted logical unit
    ! if(associated(splitted_luns_list(lun)%lun_sublist)) then
    !   do r = 1, size(data,1)
    !     call append_chunk(splitted_luns_list(lun)%lun_sublist(r), &
    !        & data(r, :, :, :, :, :))
    !   end do
    !   return
    ! end if

    dset_id = dset_id_list(lun, I_DSET)
    call h5iget_name_f(dset_id, dset_name_buf, dset_name_buf_size, name_size, err)

#if defined(METHODA)
#error "METHOD A for extendible HDF5 datasets is not implemented."
#elif defined(METHODB)
    ! (METHOD B) extend the dataset each time we append a chunk of data.

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)

    call h5sget_simple_extent_ndims_f(dspace_id, r, err)
    check_err(err)
    if(r /= rank) then
      call gkw_throw_abort("The queried HDF5 dataset "// &
         & dset_name_buf(1:name_size)//" has rank "//int2char(r,2) &
         & //" but "//int2char(rank,2)//" is needed.")
      return
    end if

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    check_err(err)
    if(any(dims(2:) /= shape(data))) then
      write (*,*) "hyperslab shape of the queried HDF5 dataset:", dims(2:)
      write (*,*) "data shape:", shape(data)
      ! FIXME Is it necessary and practicable to check if dims(1,:) has
      ! the same shape as data(:)?
      call gkw_throw_abort("Wrong data shape for "//dset_name_buf(1:name_size))
      return
    end if

    call h5dget_create_plist_f(dset_id, crp_list, err)
    check_err(err)
    call h5pget_chunk_f(crp_list, rank, dims_chunk, err)
    check_err(err)
    if(any(dims_chunk(2:) /= shape(data))) then
      write (*,*) "dims_chunk = ", dims_chunk
      write (*,*) "shape(data) = ", shape(data)
      call gkw_throw_abort("HDF5 chunk size does not match data shape.")
      return
    end if
    call h5pclose_f(crp_list, err)
    check_err(err)

    newdims = dims
    newdims(1) = dims(1) + 1
    call h5dset_extent_f(dset_id, newdims, err)
    check_err(err)

    call h5screate_simple_f (rank, dims_chunk, mem_dspace_id, err)
    check_err(err)

    call h5dget_space_f(dset_id, dspace_id, err)
    check_err(err)
    offset = 0
    offset(1) = newdims(1)-1

    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
       & offset, ones, err, ones, dims_chunk)
    check_err(err)

    call h5dwrite_f(dset_id, APPROPRIATE_HDF5_REAL_TYPE, &
       & reshape(data, dims_chunk), &
       & dims_chunk, err, &
       & mem_dspace_id, dspace_id)
    check_err(err)

    call h5sclose_f(dspace_id, err)
    call h5sclose_f(mem_dspace_id, err)
#endif

  end subroutine append_chunk_6d_real


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_1d_complex(lun, data)
    use global, only : gkw_throw_warn,int2char
    complex, dimension(:), intent(in) :: data
    integer, intent(in) :: lun
    if(is_not_root_proc()) return

    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_1d_real(lun, real(data))
      call append_chunk_1d_real(lun+1, aimag(data))
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if

  end subroutine append_chunk_1d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_2d_complex(lun, data)
    use global, only : gkw_throw_warn,int2char
    complex, dimension(:,:), intent(in) :: data
    integer, intent(in) :: lun
    if(is_not_root_proc()) return

    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_2d_real(lun, real(data))
      call append_chunk_2d_real(lun+1, aimag(data))
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if
  end subroutine append_chunk_2d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_3d_complex(lun, data)
    use global, only : gkw_throw_warn,int2char
    complex, dimension(:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    if(is_not_root_proc()) return

    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_3d_real(lun, real(data))
      call append_chunk_3d_real(lun+1, aimag(data))
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if
  end subroutine append_chunk_3d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_4d_complex(lun, data)
    use global, only : gkw_throw_warn,int2char
    complex, dimension(:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    if(is_not_root_proc()) return

    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_4d_real(lun, real(data))
      call append_chunk_4d_real(lun+1, aimag(data))
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if
  end subroutine append_chunk_4d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine append_chunk_6d_complex(lun, data)
    use global, only : gkw_throw_warn,int2char
    complex, dimension(:,:,:,:,:,:), intent(in) :: data
    integer, intent(in) :: lun
    if(is_not_root_proc()) return

    ! check if the lun is associated to a complex valued logical unit
    if(unit_is_complex(lun)) then
      call append_chunk_6d_real(lun, real(data))
      call append_chunk_6d_real(lun+1, aimag(data))
    else
      call gkw_throw_warn("Complex data should not be appended to a real valued &
         & logical unit. lun="//int2char(lun))
    end if
  end subroutine append_chunk_6d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine seek_to_chunk(lun, index)
    integer, intent(in) :: lun
    integer, intent(in) :: index

    if(is_not_root_proc()) return

    ! To keep the compiler quiet.
    if (lun > 0) continue
    if (index > 0) continue

    ! This is probably not necessary, as append_chunk() acts so by
    ! itself.

  end subroutine seek_to_chunk

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  recursive subroutine read_last_chunk_1d_real(lun, data, nlast)
    use global, only : gkw_throw_abort
    integer, intent(in) :: lun
    real, dimension(:), intent(out) :: data
    integer, intent(out) :: nlast

    integer, parameter :: rank = 2
    integer :: err
    integer(HSIZE_T), dimension(rank-1)  :: mem_dims
    integer(HID_T) :: dspace_id, mem_dspace_id, crp_list
    integer(HSIZE_T), dimension(rank) :: dims, maxdims
    integer(HSIZE_T), dimension(rank) :: offset
    integer(HSIZE_T), dimension(rank) :: dims_chunk
    integer(HSIZE_T), dimension(rank), parameter :: ones = 1
    integer :: r

    ! Evtl. check if this is indeed not a complex dataset
    ! [not implemented]
    if(is_not_root_proc()) return

    ! Check if this is a splitted dataset
    if(associated(splitted_luns_list(lun)%lun_sublist)) then
      do r = 1, size(splitted_luns_list(lun)%lun_sublist)
        call read_last_chunk_1d_real(splitted_luns_list(lun)%lun_sublist(r), &
           & data(r:r), nlast)
      end do
      return
    end if

    mem_dims = shape(data)

    call h5dget_space_f(dset_id_list(lun,I_DSET), dspace_id, err);
    call h5sget_simple_extent_ndims_f(dspace_id, r, err)
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, err)
    if(any(dims(2:) /= mem_dims)) then
      write (*,*) "hyperslab shape of the queried HDF5 dataset:", dims(2:)
      write (*,*) "expected data shape:", mem_dims
    end if
    ! write (*,*) "rank",rank,"dimensions",(dims_out[0]), &
    !    & " x ", (dims_out[1]);

    ! Find the chunk size associated to this dataset.
    call h5dget_create_plist_f(dset_id_list(lun,I_DSET), crp_list, err)
    check_err(err)
    call h5pget_chunk_f(crp_list, rank, dims_chunk, err)

    ! The chunk size and the shape of the
    ! data array in memory should be equal.
    if(any(dims_chunk(2:) /= mem_dims)) then
      write (*,*) "dims_chunk = ", dims_chunk
      write (*,*) "mem_dims = ", mem_dims
      call gkw_throw_abort("HDF5 chunk size does not match data shape in memory.")
      return
    end if

    ! Define hyperslab in the dataset. This routine wants to read
    ! the very last chunk, so we select that.
    offset = 0
    offset(1) = dims(1)-1
    nlast = int(dims(1))
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
         & offset, ones, err, ones, dims_chunk)
    check_err(err)

    ! Define the memory dataspace.
    call h5screate_simple_f(rank-1, dims_chunk(2:), mem_dspace_id, err)
    check_err(err)

    ! Define memory hyperslab.
    offset = 0
    call h5sselect_hyperslab_f(mem_dspace_id, H5S_SELECT_SET_F, &
         & offset(2:), ones(2:), err, ones(2:), dims_chunk(2:))
    check_err(err)
    
    ! Read data from hyperslab in the file into the hyperslab in
    ! memory.
    call h5dread_f(dset_id_list(lun, I_DSET), APPROPRIATE_HDF5_REAL_TYPE, &
       & data, mem_dims, err, &
       & mem_dspace_id, dspace_id, H5P_DEFAULT_F)
    check_err(err)

  end subroutine read_last_chunk_1d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_string(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val
    if(is_not_root_proc()) return

    call attach_attribute_string(dset_id_list(lun, I_DSET), key, val)
    if(unit_is_complex(lun)) then
      call attach_attribute_string(dset_id_list(lun+1, I_DSET), key, val)
    end if
    
  end subroutine attach_metadata_string


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_real(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    real, intent(in) :: val
    if(is_not_root_proc()) return

    call attach_attribute_real(dset_id_list(lun, I_DSET), key, val)
    if(unit_is_complex(lun)) then
      call attach_attribute_real(dset_id_list(lun+1, I_DSET), key, val)
    end if
    
  end subroutine attach_metadata_real


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_integer(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    integer, intent(in) :: val

    if(is_not_root_proc()) return

    call attach_attribute_integer(dset_id_list(lun, I_DSET), key, val)
    if(unit_is_complex(lun)) then
      call attach_attribute_integer(dset_id_list(lun+1, I_DSET), key, val)
    end if
    
  end subroutine attach_metadata_integer

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_logical(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    logical, intent(in) :: val

    character(len=1) :: val_char
    if(is_not_root_proc()) return

    ! There is no support for the fortran logical type in HDF5.
    
    if(val) then
      val_char = 'T'
    else
      val_char = 'F'
    end if
    call attach_metadata_string(lun, key, val_char)
    if(unit_is_complex(lun)) then
      call attach_attribute_string(dset_id_list(lun+1, I_DSET), key, val_char)
    end if
    
  end subroutine attach_metadata_logical

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_r_array(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val

    if(is_not_root_proc()) return

    call attach_attribute_r_array(dset_id_list(lun, I_DSET), key, val)
    if(unit_is_complex(lun)) then
      call attach_attribute_r_array(dset_id_list(lun+1, I_DSET), key, val)
    end if
    
  end subroutine attach_metadata_r_array

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_metadata_i_array(lun, key, val)
    integer, intent(in) :: lun
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    if(is_not_root_proc()) return

    call attach_attribute_i_array(dset_id_list(lun, I_DSET), key, val)
    if(unit_is_complex(lun)) then
      call attach_attribute_i_array(dset_id_list(lun+1, I_DSET), key, val)
    end if
    
  end subroutine attach_metadata_i_array

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine open_dataset_by_name(luname, groupname, dset_id)
    use global, only : gkw_throw_abort
    character (len=*), intent(in) :: luname, groupname
    integer(HID_T), intent(out) :: dset_id

    integer(HID_T) :: group_id
    logical :: link_exists
    integer :: err

    ! check if the desired group is there. 
    call h5lexists_f(main_output_file_id, ss(groupname), link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, ss(groupname), group_id, err)
      check_err(err)
    else
      call gkw_throw_abort(luname//" is not found")
      return
    end if
      
    ! check if a dataset with the given name exists already in the HDF5
    ! root group
    call h5lexists_f(group_id, ss(luname), &
       & link_exists, err)
    check_err(err)
    if(link_exists) then
      !get the dataset id corresponding to luname
      call h5dopen_f(group_id, ss(luname), &
         & dset_id, err)
      check_err(err)
    else
      call gkw_throw_abort(luname//" is not found.")
      return
    end if

    call h5gclose_f(group_id, err)
    check_err(err)
    
  end subroutine open_dataset_by_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  recursive subroutine attach_metadata_string_name(luname, groupname, key, val)
    use global, only : gkw_its_critical
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val

    integer(HID_T) :: dset_id
    integer :: err
    if(is_not_root_proc()) return

    if(dataset_exists_hdf5(luname//real_suffix, groupname) &
       & .and. dataset_exists_hdf5(luname//imag_suffix, groupname)) then
      call attach_metadata_string_name(luname//real_suffix, groupname, key, val)
      if(gkw_its_critical()) return
      call attach_metadata_string_name(luname//imag_suffix, groupname, key, val)
    else
      call open_dataset_by_name(luname, groupname, dset_id)
      if(gkw_its_critical()) return

      call attach_attribute_string(dset_id, key, val)
      call h5dclose_f(dset_id, err)
      check_err(err)
    end if

  end subroutine attach_metadata_string_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  recursive subroutine attach_metadata_integer_name(luname, groupname, key, val)
    use global, only : gkw_its_critical
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    integer, intent(in) :: val

    integer(HID_T) :: dset_id
    integer :: err
    if(is_not_root_proc()) return

    if(dataset_exists_hdf5(luname//real_suffix, groupname) &
       & .and. dataset_exists_hdf5(luname//imag_suffix, groupname)) then
      call attach_metadata_integer_name(luname//real_suffix, groupname, key, val)
      if(gkw_its_critical()) return
      call attach_metadata_integer_name(luname//imag_suffix, groupname, key, val)
    else
      call open_dataset_by_name(luname, groupname, dset_id)
      if(gkw_its_critical()) return

      call attach_attribute_integer(dset_id, key, val)
      call h5dclose_f(dset_id, err)
      check_err(err)
    end if

  end subroutine attach_metadata_integer_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  recursive subroutine attach_metadata_real_name(luname, groupname, key, val)
    use global, only : gkw_its_critical
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    real, intent(in) :: val

    integer(HID_T) :: dset_id
    integer :: err
    if(is_not_root_proc()) return

    if(dataset_exists_hdf5(luname//real_suffix, groupname) &
       & .and. dataset_exists_hdf5(luname//imag_suffix, groupname)) then
      call attach_metadata_real_name(luname//real_suffix, groupname, key, val)
      if(gkw_its_critical()) return
      call attach_metadata_real_name(luname//imag_suffix, groupname, key, val)
    else
      call open_dataset_by_name(luname, groupname, dset_id)
      if(gkw_its_critical()) return

      call attach_attribute_real(dset_id, key, val)
      call h5dclose_f(dset_id, err)
      check_err(err)
    end if

  end subroutine attach_metadata_real_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  recursive subroutine attach_metadata_logical_name(luname, groupname, key, val)
    use global, only : gkw_its_critical
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    logical, intent(in) :: val

    character(len=1) :: val_char

    if(is_not_root_proc()) return

    if(dataset_exists_hdf5(luname//real_suffix, groupname) &
       & .and. dataset_exists_hdf5(luname//imag_suffix, groupname)) then
      call attach_metadata_logical_name(luname//real_suffix, groupname, key, val)
      if(gkw_its_critical()) return
      call attach_metadata_logical_name(luname//imag_suffix, groupname, key, val)
    else
      ! There is no support for the fortran logical type in HDF5.
      if(val) then
        val_char = 'T'
      else
        val_char = 'F'
      end if
      call attach_metadata_string_name(luname, groupname, key, val_char)
    end if

  end subroutine attach_metadata_logical_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  recursive subroutine attach_metadata_r_array_name(luname, groupname, key, val)
    use global, only : gkw_its_critical
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val

    integer(HID_T) :: dset_id
    integer :: err
    if(is_not_root_proc()) return

    if(dataset_exists_hdf5(luname//real_suffix, groupname) &
       & .and. dataset_exists_hdf5(luname//imag_suffix, groupname)) then
      call attach_metadata_r_array_name(luname//real_suffix, groupname, key, val)
      if(gkw_its_critical()) return
      call attach_metadata_r_array_name(luname//imag_suffix, groupname, key, val)
    else
      call open_dataset_by_name(luname, groupname, dset_id)
      if(gkw_its_critical()) return

      call attach_attribute_r_array(dset_id, key, val)
      call h5dclose_f(dset_id, err)
      check_err(err)
    end if

  end subroutine attach_metadata_r_array_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  recursive subroutine attach_metadata_i_array_name(luname, groupname, key, val)
    use global, only : gkw_its_critical
    character (len=*), intent(in) :: luname, groupname
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    integer(HID_T) :: dset_id
    integer :: err
    if(is_not_root_proc()) return

    if(dataset_exists_hdf5(luname//real_suffix, groupname) &
       & .and. dataset_exists_hdf5(luname//imag_suffix, groupname)) then
      call attach_metadata_i_array_name(luname//real_suffix, groupname, key, val)
      if(gkw_its_critical()) return
      call attach_metadata_i_array_name(luname//imag_suffix, groupname, key, val)
    else
      call open_dataset_by_name(luname, groupname, dset_id)
      if(gkw_its_critical()) return

      call attach_attribute_i_array(dset_id, key, val)
      call h5dclose_f(dset_id, err)
      check_err(err)
    end if

  end subroutine attach_metadata_i_array_name

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_string(group, key, val)
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val
    
    if(is_not_root_proc()) return

    call attach_attribute_string(main_output_file_id, group//'.'//key, val)

    ! In addition, put a dataset into a special group.
    ! This is useful, because
    !   1. it may be available in the file even if the run crashes
    !   2. some analysis tools may not be able to read HDF5 attributes

    call output_string_dataset(key, input_group_name//'/'//group, (/val/))
    
  end subroutine write_run_param_string

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_string_dataset(dsetname, groupname, data)
    character (len=*), intent(in) :: dsetname, groupname
    character (len=*), dimension(1), intent(in) :: data
    logical :: link_exists
    integer, parameter :: rank = 1
    integer(HID_T) :: dtype, dspace_id, dset_id, group_id
    integer :: err
    integer(HSIZE_T), dimension(rank) :: dims ! Dataset dimensions
    integer(SIZE_T), dimension(1) :: str_len
    
    if(is_not_root_proc()) return

    ! Create file and memory datatypes.  We will save
    ! the strings as C variable length strings, H5T_STRING is defined
    ! as a variable length string.
    call h5tcopy_f(H5T_STRING, dtype, err)
    check_err(err)
    call h5tset_strpad_f(dtype, H5T_STR_NULLPAD_F, err)
    check_err(err)
    ! Create dataspace.
    dims = shape(data)
    call h5screate_simple_f(rank, dims, dspace_id, err)
    check_err(err)
    
    ! check if the desired group is there.
    call h5lexists_f(main_output_file_id, groupname, link_exists, err)
    check_err(err)
    if(link_exists) then
      ! open the group
      call h5gopen_f(main_output_file_id, groupname, group_id, err)
      check_err(err)
    else
      ! create the group
      call h5gcreate_f(main_output_file_id, groupname, group_id, err)
      check_err(err)
    end if


    ! check if the dataset is already there
    call h5lexists_f(group_id, dsetname, link_exists, err)
    check_err(err)
    if(link_exists) then
      call h5dopen_f(group_id, dsetname, dset_id, err)
      check_err(err)
    else
      ! Create the dataset and write the variable-length string data to
      ! it.
      call h5dcreate_f(group_id, dsetname, dtype, dspace_id, dset_id, err)
      check_err(err)
    end if
    
    str_len = len_trim(data)
    call h5dwrite_vl_f(dset_id, dtype, data, &
       & int((/ len(data(1)), 1/), HSIZE_T), str_len, err, dspace_id)
    check_err(err)

    ! Close and release resources.
    call h5dclose_f(dset_id, err)
    check_err(err)
    call h5sclose_f(dspace_id, err)
    check_err(err)
    call h5tclose_f(dtype, err)
    check_err(err)
    call h5gclose_f(group_id, err)
    check_err(err)
  end subroutine output_string_dataset

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_real(group, key, val)
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    real, intent(in) :: val

    if(is_not_root_proc()) return

    call attach_attribute_real(main_output_file_id, group//'.'//key, val)

    ! In addition, put a dataset into a special group.
    call output_array(key, input_group_name//'/'//group, &
       & (/val/), 'F')
    
  end subroutine write_run_param_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_complex(group, key, val)
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    complex, intent(in) :: val
    if(is_not_root_proc()) return
    

    call attach_attribute_complex(main_output_file_id, group//'.'//key, val)

    ! In addition, put a dataset into a special group.
    call output_array(key, input_group_name//'/'//group, &
       & (/ real(val), aimag(val) /), 'F')

  end subroutine write_run_param_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_integer(group, key, val)
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    integer, intent(in) :: val

    if(is_not_root_proc()) return

    call attach_attribute_integer(main_output_file_id, group//'.'//key, val)

    ! In addition, put a dataset into a special group.
    call output_array(key, input_group_name//'/'//group, &
       & (/real(val)/), 'F')
    
  end subroutine write_run_param_integer

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_logical(group, key, val)
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    logical, intent(in) :: val

    character(len=1) :: val_char
    if(is_not_root_proc()) return

    ! There is no support for the fortran logical type in HDF5.
    
    if(val) then
      val_char = 'T'
    else
      val_char = 'F'
    end if
    call write_run_param_string(group, key, val_char)

    ! In addition, put a dataset into a special group.
    call output_string_dataset(key, input_group_name//'/'//group, (/val_char/))
    
  end subroutine write_run_param_logical

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_r_array(group, key, val)
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val
    if(is_not_root_proc()) return
    
    call attach_attribute_r_array(main_output_file_id, group//'.'//key, val)

    ! In addition, put a dataset into a special group.
    call output_array(key, input_group_name//'/'//group, &
       & val, 'F')
    
  end subroutine write_run_param_r_array

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine write_run_param_i_array(group, key, val)
    character (len=*), intent(in) :: group
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val
    if(is_not_root_proc()) return
    
    call attach_attribute_i_array(main_output_file_id, group//'.'//key, val)
    
    ! In addition, put a dataset into a special group.
    call output_array(key, input_group_name//'/'//group, &
       & real(val), 'F')
    
  end subroutine write_run_param_i_array

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_attribute_integer(loc_id, key, val)
    integer(HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: key
    integer, intent(in) :: val
    
    integer :: err
    integer(HID_T) :: attr_id, aspace_id
    integer, parameter :: val_rank = 1
    integer(HSIZE_T), dimension(1), parameter :: val_dims = (/ 1 /)

    logical :: attr_exists
    
    if(is_not_root_proc()) return

    ! Check if an attribute with the given name exists already at that object.
    call h5aexists_f(loc_id, key, attr_exists, err)
    if(attr_exists) then
      ! Delete the existing attribute
      call h5adelete_f(loc_id, key, err)
    end if

    ! Create data space for the attribute.
    call h5screate_simple_f(val_rank, val_dims, aspace_id, err)

    ! Create the attribute.
    call h5acreate_f(loc_id, key, &
       & H5T_NATIVE_INTEGER, aspace_id, attr_id, err)
    ! Write the attribute data.
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, val, val_dims, err)
    
    ! Terminate access
    call h5aclose_f(attr_id, err)
    call h5sclose_f(aspace_id, err)
    
  end subroutine attach_attribute_integer


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_attribute_real(loc_id, key, val)
    integer(HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: key
    real, intent(in) :: val

    integer :: err
    integer(HID_T) :: attr_id, aspace_id
    integer, parameter :: val_rank = 1
    integer(HSIZE_T), dimension(1), parameter :: val_dims = (/ 1 /)

    logical :: attr_exists

    if(is_not_root_proc()) return

    ! Check if an attribute with the given name exists already at that object.
    call h5aexists_f(loc_id, key, attr_exists, err)
    if(attr_exists) then
      ! Delete the existing attribute
      call h5adelete_f(loc_id, key, err)
    end if

    ! Create data space for the attribute.
    call h5screate_simple_f(val_rank, val_dims, aspace_id, err)

    ! Create the attribute.
    call h5acreate_f(loc_id, key, &
       & APPROPRIATE_HDF5_REAL_TYPE, aspace_id, attr_id, err)
    ! Write the attribute data.
    call h5awrite_f(attr_id, APPROPRIATE_HDF5_REAL_TYPE, (/val/), val_dims, err)
    
    ! Terminate access
    call h5aclose_f(attr_id, err)
    call h5sclose_f(aspace_id, err)
    
  end subroutine attach_attribute_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_attribute_string(loc_id, key, val)
    integer(HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: key
    character (len=*), intent(in) :: val

    integer :: err
    integer(HID_T) :: attr_id, aspace_id, mem_type_id
    integer, parameter :: val_rank = 1
    integer(HSIZE_T), dimension(1), parameter :: val_dims = (/ 1 /)
    integer(SIZE_T) :: string_len

    logical :: attr_exists
    
    if(is_not_root_proc()) return

    ! Check if an attribute with the given name exists already at that object.
    call h5aexists_f(loc_id, key, attr_exists, err)
    if(attr_exists) then
      ! Delete the existing attribute
      call h5adelete_f(loc_id, key, err)
    end if

    ! Create a new attribute.

    ! Create data space for the attribute.
    call h5screate_simple_f(val_rank, val_dims, aspace_id, err)

    ! Create datatype for the attribute.
    call h5tcopy_f(H5T_NATIVE_CHARACTER, mem_type_id, err)
    string_len = len(val)
    call h5tset_size_f(mem_type_id, string_len, err)

    ! Create the attribute.
    call h5acreate_f(loc_id, key, &
       & mem_type_id, aspace_id, attr_id, err)

    ! Write the attribute data.
    call h5awrite_f(attr_id, mem_type_id, val, val_dims, err)
    
    ! Terminate access
    call h5aclose_f(attr_id, err)
    call h5sclose_f(aspace_id, err)
    
    
  end subroutine attach_attribute_string

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_attribute_complex(loc_id, key, val)
    integer(HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: key
    complex, intent(in) :: val

    integer :: err
    integer(HID_T) :: attr_id, aspace_id
    integer, parameter :: val_rank = 1
    integer(HSIZE_T), dimension(1), parameter :: val_dims = (/ 2 /)

    logical :: attr_exists
    
    if(is_not_root_proc()) return

    ! Check if an attribute with the given name exists already at that object.
    call h5aexists_f(loc_id, key, attr_exists, err)
    if(attr_exists) then
      ! Delete the existing attribute
      call h5adelete_f(loc_id, key, err)
    end if

    ! Create data space for the attribute.
    call h5screate_simple_f(val_rank, val_dims, aspace_id, err)

    ! Create the attribute.
    call h5acreate_f(loc_id, key, &
       & APPROPRIATE_HDF5_REAL_TYPE, aspace_id, attr_id, err)
    ! Write the attribute data.
    call h5awrite_f(attr_id, APPROPRIATE_HDF5_REAL_TYPE, &
       & (/real(val), aimag(val)/), val_dims, err)
    
    ! Terminate access
    call h5aclose_f(attr_id, err)
    call h5sclose_f(aspace_id, err)

  end subroutine attach_attribute_complex


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_attribute_r_array(loc_id, key, val)
    integer(HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: key
    real, dimension(:), intent(in) :: val

    integer :: err
    integer(HID_T) :: attr_id, aspace_id
    integer, parameter :: val_rank = 1
    integer(HSIZE_T) :: val_dim1

    logical :: attr_exists
    
    if(is_not_root_proc()) return

    ! Check if an attribute with the given name exists already at that object.
    call h5aexists_f(loc_id, key, attr_exists, err)
    if(attr_exists) then
      ! Delete the existing attribute
      call h5adelete_f(loc_id, key, err)
    end if

    val_dim1 = size(val)

    ! Create data space for the attribute.
    call h5screate_simple_f(val_rank, (/val_dim1/), aspace_id, err)

    ! Create the attribute.
    call h5acreate_f(loc_id, key, &
       & APPROPRIATE_HDF5_REAL_TYPE, aspace_id, attr_id, err)
    ! Write the attribute data.
    call h5awrite_f(attr_id, APPROPRIATE_HDF5_REAL_TYPE, val, (/val_dim1/), err)

    ! Terminate access
    call h5aclose_f(attr_id, err)
    call h5sclose_f(aspace_id, err)
  end subroutine attach_attribute_r_array
  

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine attach_attribute_i_array(loc_id, key, val)
    integer(HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    integer :: err
    integer(HID_T) :: attr_id, aspace_id
    integer, parameter :: val_rank = 1
    integer(HSIZE_T) :: val_dim1
    logical :: attr_exists
    
    if(is_not_root_proc()) return

    ! Check if an attribute with the given name exists already at that object.
    call h5aexists_f(loc_id, key, attr_exists, err)
    if(attr_exists) then
      ! Delete the existing attribute
      call h5adelete_f(loc_id, key, err)
    end if

    ! Create data space for the attribute.
    val_dim1 = size(val)
    call h5screate_simple_f(val_rank, (/val_dim1/), aspace_id, err)

    ! Create the attribute.
    call h5acreate_f(loc_id, key, &
       & H5T_NATIVE_INTEGER, aspace_id, attr_id, err)
    ! Write the attribute data.
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, val, (/val_dim1/), err)
    
    ! Terminate access
    call h5aclose_f(attr_id, err)
    call h5sclose_f(aspace_id, err)
    
  end subroutine attach_attribute_i_array


#endif
end module io_hdf5

