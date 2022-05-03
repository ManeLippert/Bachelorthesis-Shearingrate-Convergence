!===========================================================================
!> This module provides the abstracted interface for the matrix format.
!> It is an interface wrapper around swappable modules matrix_format_* which
!> perform the actual sparse matrix computations for different formats.
!>
!> The recognized output options (selected by input variable io_format) are
!> gkw-crs     : self implemented compressed row storage (called CRS,
!>               CSR or Yale format)
!> librsb      : the publicly available library written by M.Martone
!>
!===========================================================================
module matrix_format

#define FIRST_TOUCH

#define CRS_FASTDIAGONAL_FORMAT

#ifdef HAVE_MKL
! UPPERCASE USE to deliberately exclude from mkdeps script
  ! USE f95_precision, only: wp => dp

  !lowercase 'use', to include mkl_sblas module as an additional
  ! dependency.  SRG: it seems to me that installing the MKL does not
  ! compile a .o and .mod files for sparse BLAS. One has to compile
  ! them on it's own.
  USE mkl_spblas

  ! USE ISO_C_BINDING
  ! INCLUDE 'mkl_spblas.fi'
  ! INCLUDE 'mkl_sparse_handle.fi'
#endif

  implicit none

  public

  !> coordinate format. This is the simplest and most flexible format.
  integer, parameter :: FORMAT_GKWCOO = 0

  character(len=10) :: matrix_format_gkwcrs = 'gkw-crs'
  integer, parameter :: FORMAT_GKWCRS = 1
  character(len=10) :: matrix_format_librsb = 'librsb'
  integer, parameter :: FORMAT_LIBRSB = 2
  character(len=10) :: matrix_format_mklcrs = 'mkl-crs'
  integer, parameter :: FORMAT_MKLCRS = 3
  
  type :: sparse_matrix
    ! Initially, every sparse matrix is constructed using the arrays below
    ! and a simple coordinate format. As soon as the matrix is filled,
    ! it can be transformed to a (potentially) more efficient sparse
    ! matrix format.

    ! FIXME find out if it does have any significance to put the
    ! elements of this struct in a certain order (e.g. scalars first)
    ! (maybe to help for good alignment or so).

    !> the total number of nonzero elements
    integer :: nmat = 0
    !> counter for the number of compressions
    integer :: compressions = 0
    !> the matrix name is use for error messages, for example
    character(len=80) :: name
    !> the format that this matrix is going to be converted to, as
    !> soon as all values are filled in. This can be influenced by the
    !> inputfile, at least for some matrices.
    character(len=10) :: format
    !> the format this matrix has currently. If finish_matrix has
    !> never been called, this matrix is still in coordinate (GKWCOO)
    !> format.
    integer :: current_format = FORMAT_GKWCOO
    !> this is true as soon as the finish_matrix() routine has been
    !> called one time.
    logical :: is_finished = .false.
    !> Do not store values in memory but only count the number of
    !> elements.  (Afterwards, the nmat can then be used for
    !> allocation) This extra flag is used because it is said that checks
    !> like if(allocated(mat)) may not always work as expected.
    logical :: count_only = .true.

    !> this array contains the values of the nonzero matrix elements,
    !> and will be initialised for good data locality
    complex, allocatable :: mat(:)
    !> Integer array that determines the row of the corresponding
    !> element in mat
    integer, allocatable :: ii(:)
    !> Integer array that determines the column of the corresponding
    !> element in mat, and will be initialised for good data locality
    integer, allocatable :: jj(:)

    !! further datastructures needed by gkw-crs format:
    !> row of first nonzero element
    integer :: start_row
    !> row of last nonzero element
    integer :: end_row
    !> number of rows between first and last nonzero element
    integer :: nrows

    !> row of first nonzero element
    integer :: start_col
    !> row of last nonzero element
    integer :: end_col
    !> number of rows between first and last nonzero element
    integer :: ncols
    
    !> array for the starting element of each matrix row
    integer, allocatable :: irs(:)
#ifdef CRS_FASTDIAGONAL_FORMAT
    !> array which holds alternately an row start index for a purely
    !> diagonal block, followed by a row start index for a crs format
    !> block.  the indices are row indices.
    integer, allocatable :: diag_block_row_start(:)
    !> the number of valid blocks specified by
    !> diag_block_row_start. (Due to the algorithm, the array is
    !> allocated with a larger size in general.)
    integer :: ndiag_blocks = 0

    !> this array describes if the diagonal blocks are unity matrix
    !> blocks, respectively. For simplicity, this has the same
    !> structure as diag_block_row_start, i.e. only the odd indices
    !> belong to diagonal blocks and are relevant here.
    logical, allocatable :: diag_block_is_unity(:)
#endif

    !> true if the matrix is purely diagonal, without a single
    !> off-diagonal element. There may still be zeros on the main diagonal,
    !> of course.
    logical :: is_diagonal = .false.

    !> for debugging: number of invalid elements that have been rejected.
    integer :: nrejected_elements

#ifdef HAVE_LIBRSB
    !> a librsb matrix handle
    integer :: librsb_handle
#endif

#ifdef HAVE_MKL
    ! Matrix descriptor
    type(matrix_descr) :: mkl_desc
    ! CSR matrix representation
    type(sparse_matrix_t) :: mkl_handle
#endif
  end type sparse_matrix

  interface create_matrix
    module procedure create_matrix_dry
    module procedure create_matrix_from_mat
    module procedure create_matrix_from_name
  end interface create_matrix

  interface put_element
    module procedure put_element_simple
    module procedure put_element_struct
    module procedure put_element_matrix
  end interface put_element

  interface usmv
    module procedure usmv_xy
    module procedure usmv_yupdate
  end interface usmv
  
  contains

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine matrix_format_init()
      use control, only : matrix_format
      use global, only : compiled_with_librsb
      use general, only : gkw_abort
#ifdef HAVE_LIBRSB
      ! UPPERCASE USE to deliberately exclude from mkdeps script
      USE rsb
      USE blas_sparse
      type(C_PTR),parameter :: init_null = C_NULL_PTR
#endif
      integer :: istat

      if(.not.compiled_with_librsb .and. matrix_format==matrix_format_librsb) then
        matrix_format = matrix_format_gkwcrs
      end if

      if(matrix_format == matrix_format_librsb) then
#ifdef HAVE_LIBRSB
        ! librsb initialization
        ! call rsb_lib_init_np(istat)
        !istat = rsb_lib_init(RSB_NULL_INIT_OPTIONS)
        istat = rsb_lib_init(init_null)
#else
        istat = 1
#endif
        if(istat /= 0) call gkw_abort('librsb could not be initialised')
      end if

    end subroutine matrix_format_init

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine matrix_format_finalize()
      use control, only : matrix_format
      use general, only : gkw_warn
#ifdef HAVE_LIBRSB
      ! UPPERCASE USE to deliberately exclude from mkdeps script
      USE rsb
      USE blas_sparse
      type(C_PTR),parameter :: exit_null = C_NULL_PTR
#endif
      integer :: istat

      if(matrix_format == matrix_format_librsb) then
#ifdef HAVE_LIBRSB
        !istat = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS)
        istat = rsb_lib_init(exit_null)
#else
        istat = 1
#endif
        if(istat /= 0) call gkw_warn('librsb could not be finalized cleanly')
      end if


    end subroutine matrix_format_finalize

    !--------------------------------------------------------------------------
    !> Create a matrix struct, but do not allocate arrays. This can be
    !> used to only count elements.
    !--------------------------------------------------------------------------
    function create_matrix_dry(name, force_format) result(mat)
      use control, only : matrix_format
      use general, only : gkw_abort
      character(len=*), intent(in) :: name
      character(len=10), intent(in), optional :: force_format
      type(sparse_matrix) :: mat

      if(present(force_format)) then
        if(force_format /= matrix_format_gkwcrs) then
          ! programming mistake
          call gkw_abort('A format other than gkw-crs has been enforced.  &
             & This is not allowed.')
        end if
      end if

      mat = create_matrix(name, 0, force_format, .true.)
    end function create_matrix_dry

    !--------------------------------------------------------------------------
    !> Given an sparse_matrix type struct, this function returns a struct with
    !> the arrays allocated.
    !> This can be used for first dry counting of the number of elements, then
    !> calling this routine, then actually put the elements.
    !--------------------------------------------------------------------------
    function create_matrix_from_mat(origmat) result(mat)
      use general, only : gkw_warn
      type(sparse_matrix), intent(in) :: origmat
      type(sparse_matrix) :: mat

      if(origmat%nmat == 0) then
        call gkw_warn('Matrix will be empty: '//trim(origmat%name))
      end if
      mat = create_matrix(origmat%name, origmat%nmat, origmat%format, .false.)
    end function create_matrix_from_mat
    

    !--------------------------------------------------------------------------
    !> This function returns a sparse_matrix struct with the arrays allocated,
    !> ready to put elements.
    !--------------------------------------------------------------------------
    function create_matrix_from_name(name, max_nelements, force_format, &
       & count_only) result(mat)
      use control, only : root_and_not_silent
      use control, only : matrix_format
      character(len=*), intent(in) :: name
      !> pass max_nelements=0 to create a matrix which only counts
      !> elements, without allocated data.
      !>
      !> pass max_nelements>0 to get a matrix ready to store matrix elements.
      integer, intent(in) :: max_nelements
      !> some matrices are used in ways which are not (yet) clearly
      !> brought into patters which can be computed using the usmv
      !> routines (or will never be).  For such cases the builtin CRS
      !> format can be enforced.
      character(len=10), intent(in), optional :: force_format
      logical, intent(in), optional :: count_only
      type(sparse_matrix) :: mat
      if(present(count_only)) then
        mat%count_only = count_only
      else
        mat%count_only = .false.
      end if

      if(.not. mat%count_only) then
        if(allocated(mat%ii)) deallocate(mat%ii)
        if(allocated(mat%jj)) deallocate(mat%jj)
        if(allocated(mat%mat)) deallocate(mat%mat)

        allocate(mat%ii(max(0,max_nelements)))
        allocate(mat%jj(max(0,max_nelements)))
        allocate(mat%mat(max(0,max_nelements)))
      end if
      mat%name = name
      mat%nmat = 0
      mat%compressions = 0
      mat%nrejected_elements = 0
      mat%start_row = huge(1)
      mat%end_row = 0
      mat%start_col = huge(1)
      mat%end_col = 0
      mat%is_finished = .false.
      ! this is a workaround: because creating the librsb matrix will
      ! destroy simple coordinate data (which we need for several
      ! other matrices), the more sophisticated matrix formats are used
      ! only for certain matrices.
      if(present(force_format)) then
        mat%format = force_format
      else
        mat%format = matrix_format
      end if

      if (root_and_not_silent) then
        write(*,*) '********************************************'
        write(*,*) 'Matrix created'
        write(*,*) '(',trim(mat%name),')'
        write(*,*) max_nelements,'max. elements'
        write(*,*) trim(mat%format),' format'
        write(*,*) '********************************************'
      end if

    end function create_matrix_from_name

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine detect_diagonal(mat)
      type(sparse_matrix), intent(inout) :: mat
      integer :: ir

      mat%is_diagonal = .true.
      do ir = mat%start_row, mat%end_row
        if(.not. &
           ! either exactly one element in this row, and this element
           ! is on the main diagonal
           & ((mat%irs(ir+1)-mat%irs(ir) == 1  .and. mat%jj(mat%irs(ir)) == ir) &
           ! or zero elements in this row
           & .or. mat%irs(ir+1)-mat%irs(ir) == 0)) then
          mat%is_diagonal = .false.
          exit
        end if
      end do
    end subroutine detect_diagonal

#ifdef CRS_FASTDIAGONAL_FORMAT
    !--------------------------------------------------------------------------
    !> This routine detects diagonal and non-diagonal blocks in the matrix.
    !> This information is stored in mat%diag_block_row_start and can be
    !> used to make matrix-vector calculations in CRS format faster.
    !--------------------------------------------------------------------------
    subroutine find_diagonal_blocks(mat)
      use global, only : r_tiny
      type(sparse_matrix), intent(inout) :: mat
      integer :: nelems
      integer :: ir, i, j
      integer, parameter :: MIN_DIAG_BLOCK_SIZE = 100
      nelems = 1

      outer_loop: do while (.true.)
        if(allocated(mat%diag_block_row_start)) then
          deallocate(mat%diag_block_row_start)
          deallocate(mat%diag_block_is_unity)
        end if
        nelems = nelems * 2
        allocate(mat%diag_block_row_start(nelems+1))
        allocate(mat%diag_block_is_unity(nelems+1))
        mat%diag_block_is_unity = .false.

        mat%ndiag_blocks = 0
        ir = mat%start_row
        do while(ir <= mat%end_row)
          ! first look for a diagonal block
          mat%ndiag_blocks = mat%ndiag_blocks + 1
          mat%diag_block_is_unity(mat%ndiag_blocks) = .true.
          if(mat%ndiag_blocks > nelems) cycle outer_loop

          mat%diag_block_row_start(mat%ndiag_blocks) = ir

          ! while there is only one element in the row,
          ! and this element is on the diagonal
          do while(mat%irs(ir)+1 == mat%irs(ir+1) .and. &
             & mat%jj(mat%irs(ir)) == ir)
            if(abs(mat%mat(ir) - 1.0) < r_tiny) then
              mat%diag_block_is_unity(mat%ndiag_blocks) = .false.
            end if

            ir = ir + 1
            ! fortran logic operators may or may not short circuit.
            ! therefore, I cannot write
            !    while(ir <= mat%end_row .and. ...)
            ! because this would still let it access element irs(ir+1)
            ! (out of bounds).
            ! Instead we must catch the case of the last line separately:
            if(ir > mat%end_row) then
              exit
            end if
          end do

          ! now find the size of the unstructured sparse block
          mat%ndiag_blocks = mat%ndiag_blocks + 1
          if(mat%ndiag_blocks > nelems) cycle outer_loop

          mat%diag_block_row_start(mat%ndiag_blocks) = ir
          if(ir <= mat%end_row) then
            ! while there is not exactly one element in the row or
            ! there is at least one off-diagonal element
            do while(.not. (mat%irs(ir)+1 == mat%irs(ir+1).and. &
               & mat%jj(mat%irs(ir)) == ir))
              ir = ir + 1

              ! fortran logic operators may or may not short circuit.
              if(ir > mat%end_row) then
                exit
              end if
            end do
          end if
        end do
        exit outer_loop
      end do outer_loop

      ! write (*,*) trim(mat%name), mat%ndiag_blocks, &
      !    & 'diag_block_row_start:', mat%diag_block_row_start(1:mat%ndiag_blocks)

      ! now there are likely some diagonal blocks which consist only
      ! of very few elements.
      do i = 3, mat%ndiag_blocks, 2
        ! a block with uneven number is a diagonal block.

        if(mat%diag_block_row_start(i) - mat%diag_block_row_start(i+1) &
           & < MIN_DIAG_BLOCK_SIZE) then
          ! this block is too small.

          ! merge this block to the preceeding one by simply deleting
          ! it. it is deleted by shifting all consecutive blocks by
          ! one position to the left (if there are any).
          do j = i, mat%ndiag_blocks-1-2
            mat%diag_block_row_start(j) = mat%diag_block_row_start(j+2)
          end do
          mat%ndiag_blocks = mat%ndiag_blocks - 2
        end if
      end do

      mat%diag_block_row_start(mat%ndiag_blocks+1) = mat%end_row + 1

      ! write (*,*) trim(mat%name), mat%ndiag_blocks, &
      !    & 'optdiag_block_row_start:', mat%diag_block_row_start(1:mat%ndiag_blocks)


    end subroutine find_diagonal_blocks

#endif

    !--------------------------------------------------------------------------
    !> As soon as all elements have been added to the matrix, it can be
    !> 'finished'. Finishing means that the data is reorganised to a
    !> (hopefully) more efficient datastructure, i.e. a specific sparse matrix
    !> format.
    !>
    !> If a matrix is never finished, it is in the most simple coordinate format
    !> where mat(i) is the value of the element at indices ii(i) and jj(i).
    !--------------------------------------------------------------------------
    subroutine finish_matrix(mat)
      use general, only : gkw_abort
      use global, only : i_huge
#ifdef HAVE_LIBRSB
      ! UPPERCASE USE to deliberately exclude from mkdeps script
      USE rsb
      USE blas_sparse
#endif
#ifdef HAVE_MKL
      !lowercase 'use', to include mkl_sblas module as an additional
      ! dependency.
      USE mkl_spblas
#endif
      type(sparse_matrix), intent(inout) :: mat
#ifdef HAVE_LIBRSB
      integer :: nmat
#endif
#ifdef HAVE_MKL
      type(sparse_matrix_t) :: coo_handle
#endif
#if (defined(HAVE_LIBRSB) || defined(HAVE_MKL))
      integer :: istat
#endif
#ifdef _OPENMP
#ifdef FIRST_TOUCH
      !> this array contains the values of the nonzero matrix elements
      !> during the construction of the matrix
      complex, allocatable :: mat_(:)
      !> Integer array that determines the column of the corresponding
      !> element in mat, during construction of the matrix
      integer, allocatable :: jj_(:)
#endif
#endif

      if(mat%is_finished) call gkw_abort('finish_matrix: is called multiple times for ' //trim(mat%name))
      ! do not finish, while we are dry counting matrix elements.
      if(mat%count_only) return

      call compress_matrix(mat)

      if(mat%nmat > 0) then
        mat%start_row = minval(mat%ii(1:mat%nmat))
        mat%end_row = maxval(mat%ii(1:mat%nmat))

        mat%start_col = minval(mat%jj(1:mat%nmat))
        mat%end_col = maxval(mat%jj(1:mat%nmat))
      else
        mat%end_row = 0
        mat%start_row = 1
        mat%end_col = 0
        mat%start_col = 1
      end if
      mat%nrows = mat%end_row - mat%start_row + 1
      mat%ncols = mat%end_col - mat%start_col + 1

      if(mat%format == matrix_format_gkwcrs) then

        ! + 1 for convenience, otherwise one would need to treat the
        ! last index exceptionally
        allocate(mat%irs(mat%start_row:mat%end_row + 1))
        call crs_find_row_index(mat)
        call detect_diagonal(mat)
#ifdef CRS_FASTDIAGONAL_FORMAT
        ! detect the number of purely diagonal and not-purely-diagonal
        ! blocks in this matrix.
        call find_diagonal_blocks(mat)
#endif

#ifdef _OPENMP
#ifdef FIRST_TOUCH
        ! allocate temporary buffers
        allocate(mat_(1:size(mat%mat)))
        allocate(jj_(1:size(mat%mat)))
        ! put the data into the buffers
        mat_ = mat%mat
        jj_ = mat%jj
        ! allocate the memory anew
        deallocate(mat%mat)
        deallocate(mat%jj)
        allocate(mat%mat(1:size(mat_)))
        allocate(mat%jj(1:size(jj_)))
        ! and this time do the first touch of every memory location in
        ! the same pattern as the usmv routine will access it later
        call crs_first_touch(mat, mat_, jj_)
        ! free the buffers again
        deallocate(mat_)
        deallocate(jj_)
#endif
#endif

        ! if(root_and_not_silent) then
        !   write(*,*)'finish ',mat%name
        !   write(*,*)'irs',mat%irs
        !   write(*,*)'ii',mat%ii(1:mat%nmat)
        !   write(*,*)'jj',mat%jj(1:mat%nmat)
        !   write(*,*)'mat', mat%mat(1:mat%nmat)
        !   write(*,*)'end finish ',mat%name
        ! end if
        mat%current_format = FORMAT_GKWCRS
      elseif(mat%format == matrix_format_mklcrs) then
#ifdef HAVE_MKL
          ! Create COO matrix
          istat = mkl_sparse_z_create_coo(coo_handle,sparse_index_base_one, &
             & mat%end_row,mat%end_col,mat%nmat,mat%ii,mat%jj,mat%mat)
          ! convert to CSR format
          istat = mkl_sparse_convert_csr(coo_handle, &
             & sparse_operation_non_transpose, mat%mkl_handle)
        ! create matrix descriptor
        mat%mkl_desc%type = sparse_matrix_type_general
        !mat%mkl_desc%mode and mat%mkl_desc%diag are not of interest
        !for a general matrix

        !estimate of number and type of upcoming matrix-vector operations.
        ! istat = mkl_sparse_set_mv_hint (mat%mkl_handle, &
        !    & SPARSE_OPERATION_NON_TRANSPOSE, mat%desc, i_huge)
        ! are these other hints necessary?
        ! istat = mkl_sparse_set_sv_hint (mat%mkl_handle, &
        !    & SPARSE_OPERATION_NON_TRANSPOSE, mat%desc, 0)
        ! I am not sure what the layout parameter means
        !mkl_sparse_set_mm_hint (A, operation, descr, SPARSE_LAYOUT_COLUMN_MAJOR,&
           ! & mat%nmat, 0)

        ! Allow the routine to allocate memory up to the size of matrix A
        ! for converting into the appropriate sparse format.
        ! istat = mkl_sparse_set_memory_hint(mat%mkl_handle, &
        !    & SPARSE_MEMORY_AGGRESSIVE)
        
        ! chose proper kernels and workload balancing strategy based
        ! on hints given above
        istat = mkl_sparse_optimize(mat%mkl_handle)
        mat%current_format = FORMAT_MKLCRS
#endif
      else if(mat%format == matrix_format_librsb) then
#ifdef HAVE_LIBRSB
        ! begin matrix creation
        call zuscr_begin(mat%end_row,mat%end_col,mat%librsb_handle,istat)
        if(istat /= 0) then
          call gkw_abort('librsb matrix could not be created for '//mat%name)
        end if

        ! set symmetry property
        !call ussp(m12_rsb,blas_lower_symmetric,istat)

        nmat = mat%nmat
        ! insert some nonzeroes
        call uscr_insert_entries(mat%librsb_handle,mat%nmat,mat%mat(1:nmat),&
           & mat%ii(1:nmat),mat%jj(1:nmat),istat)
        if(istat /= 0) then
          call gkw_abort('librsb matrix could not be filled for ' //mat%name)
        end if

        ! end matrix creation
        call uscr_end(mat%librsb_handle,istat)
        if(istat /= 0) then
          call gkw_abort('librsb matrix could not be finished for ' //mat%name)
        end if
        mat%current_format = FORMAT_LIBRSB
#else
        call gkw_abort('finish_matrix: unsupported matrix format librsb ' //trim(mat%name))
#endif
      else
        call gkw_abort('code error: unsupported matrix format ''' //trim(mat%name)//''', format '//trim(mat%format))
      end if
      mat%is_finished = .true.

      call print_matrix_report(mat)
    end subroutine finish_matrix

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine print_matrix_report(mat)
      use control, only : root_and_not_silent
      use mpiinterface, only : root_processor
      type(sparse_matrix), intent(in) :: mat
      if (root_processor) then
        write(*,*) 'matrix "' // trim(mat%name) // '" has format '// trim(mat%format)
      end if
      if (root_and_not_silent) then
        write(*,*) '------------------------------------------------------'
        write(*,*) 'Matrix report:', trim(mat%name)
        write(*,*) 'current format:', mat%current_format
        write(*,*) 'final format:', trim(mat%format)
300     format(I8,' nonzero elements. Max ', I8,' elements')
224     format(' first row ',I8, ',last row ', I8,', in total ',I8,' rows. Full compression #: ',I8)
225     format(' first col ',I8, ',last col ', I8,', in total ',I8,' columns.')
        write(*,300) mat%nmat,size(mat%mat)
        write(*,224) mat%start_row, mat%end_row, mat%nrows, mat%compressions
        write(*,225) mat%start_col, mat%end_col, mat%ncols
        write(*,*) 'rejected elements:', mat%nrejected_elements
        write(*,*) 'is purely diagonal:', mat%is_diagonal
        write(*,*) 'is finished:', mat%is_finished
        write(*,*) 'is only dry counting:', mat%count_only
#ifdef CRS_FASTDIAGONAL_FORMAT
        write(*,*) 'blocks:', mat%ndiag_blocks
        if(mat%ndiag_blocks > 0) then
          write(*,*) 'diagonal blocks being unity:', &
             & mat%diag_block_is_unity(1:mat%ndiag_blocks:2)
          write(*,*) 'blocks row starts:', &
             & mat%diag_block_row_start(1:mat%ndiag_blocks)
        end if
#endif
        write(*,*) '------------------------------------------------------'
        write(*,*)
      end if
    end subroutine print_matrix_report

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine finalize_matrix(mat)
#ifdef HAVE_MKL
      !lowercase 'use', to include mkl_sblas module as an additional
      ! dependency.
      USE mkl_spblas
#endif
      type(sparse_matrix), intent(inout) :: mat
#if HAVE_MKL
      integer :: istat
#endif

      deallocate(mat%ii)
      deallocate(mat%jj)
      deallocate(mat%mat)

      if(mat%current_format == FORMAT_GKWCRS) then
        if(allocated(mat%irs)) deallocate(mat%irs)
      else if(mat%current_format == FORMAT_MKLCRS) then
#if HAVE_MKL
        ! Release internal representation of CSR matrix
        istat = mkl_sparse_destroy(mat%mkl_handle)
#endif
      end if
    end subroutine finalize_matrix

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine put_element_struct(mat,E)
      use structures, only : matrix_element
      use general, only : gkw_abort
      type(sparse_matrix), intent(inout) :: mat
      type (matrix_element), intent(in) :: E

      if (E%term == '') call gkw_abort('elem%term not set')
      if (E%ideriv == - 100) call gkw_abort('elem%ideriv not set ' // E%term)
      if (E%itloc == -200) call gkw_abort('elem%itloc not set ' // E%term)
      if (E%itype == -300) call gkw_abort('elem%itype not set ' // E%term)

      if (abs(E%val) > 1.E10) then
        write(*,*) 'Very large matrix element from term ' // trim(E%term)
      end if

      call put_element_simple(mat, get_elem_ii(E), get_elem_jj(E), E%val)

    end subroutine put_element_struct

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine put_element_simple(mat,i,j,val)
      use general, only : gkw_abort, gkw_warn
      use global, only : r_tiny
      type(sparse_matrix), intent(inout) :: mat
      integer, intent(in) :: i,j
      complex, intent(in) :: val
      if(abs(val) < r_tiny) then
        if(i == j) then
          if(mat%name == "field matrix, diag part") then
            ! this zero valued element is supposed to be on the main
            ! diagonal.  If this is intended for one of the purely
            ! diagonal matrices abort: since the non-csr algorithm
            !do i = start, end
            !  fdis(ii(i)) = mat(i)*fdis(ii(i))
            !end do
            ! that was used for diagonal matrices
            ! implicitely assumes unity for empty matrix rows, a 0.0
            ! element not being added to this matrix implicitely
            ! becomes 1.0 value. This bug may never play a role if we never
            ! attempt to add zeros and later use this algorithm:

            ! abort, to make sure this case does not go unnoticed.

            ! if this aborts for you, then your old runs of this input
            ! file may evtl. suffer from a small bug. You may replace
            ! the abort line with the warning and rerun. This will
            ! then explicitely add the 0.0 to the matrix. Please check
            ! what difference it makes to the old data and report
            ! back.
            call gkw_abort('A 0.0 value matrix element on the main diagonal was &
             &added to '//trim(mat%name))
            ! call gkw_warn('A 0.0 value matrix element on the main diagonal was &
            !  &added to '//trim(mat%name))

            ! note that this is a workaround which should be removed as
            ! soon as the fast non-csr usmv_yupdate() algorithm is
            ! corrected and the respective diagonal matrices are
            ! equipped with explicit unity blocks.
          else
            return
          end if
        else
          ! if not on the diagonal then do not add this zero valued
          ! element, because the CRS machinery does not have the
          ! problem of implicit assumption of unity elements described above.
          return
        end if
      end if

      ! Check if main matrix is getting too big
      if (.not. mat%count_only .and. mat%nmat+1 > size(mat%mat)) then
        call compress_matrix(mat)
        if (mat%nmat+1 > size(mat%mat)) then
          call gkw_abort('matrix '//trim(mat%name)//' is too small')
        end if
      end if
      ! do only put elements with reasonable indices
      if(i > 0 .and. j > 0) then
        mat%nmat = mat%nmat + 1
        if(.not. mat%count_only) then
          mat%ii(mat%nmat) = i
          mat%jj(mat%nmat) = j
          mat%mat(mat%nmat) = val
        end if
        ! update these on the fly, too
        mat%start_row = min(mat%start_row, i)
        mat%end_row = max(mat%end_row, i)
        mat%start_col = min(mat%start_col, j)
        mat%end_col = max(mat%end_col, j)
        mat%nrows = mat%end_row - mat%start_row + 1
        mat%ncols = mat%end_col - mat%start_col + 1
      else
        mat%nrejected_elements = mat%nrejected_elements + 1
      end if
    end subroutine put_element_simple

    !--------------------------------------------------------------------------
    !>
    !--------------------------------------------------------------------------
    subroutine put_element_matrix(mat,matval)
      use general, only : gkw_abort
      type(sparse_matrix), intent(inout) :: mat
      type(sparse_matrix), intent(in) :: matval
      integer :: ir, i
      if(matval%current_format == FORMAT_LIBRSB) then
#ifdef HAVE_LIBRSB
        call gkw_abort('it is not possible to access single matrix elements')
#endif
      else if(matval%current_format == FORMAT_MKLCRS) then
#if HAVE_MKL
        call gkw_abort('it is not implemented to access single matrix elements')
#endif
      else if(matval%current_format == FORMAT_GKWCRS) then
        do ir=matval%start_row, matval%end_row
          !loop over elements in row:
          do i = matval%irs(ir),matval%irs(ir+1)-1
            call put_element_simple(mat, ir, matval%jj(i), matval%mat(i))
          end do
        end do

      else if(matval%current_format == FORMAT_GKWCOO) then
        do i = 1,matval%nmat
          call put_element_simple(mat, matval%ii(i), matval%jj(i), matval%mat(i))
        end do
      end if

    end subroutine put_element_matrix



    !--------------------------------------------------------------------------
    !> compute y = y + a * mat * x
    !>
    !> usmv means:
    !> us means unstructured sparse
    !> mv means matrix vector multiplication
    !>
    !> usaxpy means:
    !> us means unstructured sparse
    !> axpy means scalar a times vector x plus vector y
    !>
    !> NOTE: it is dangerous to pass the same vector to x and y.
    !> In General the result will be wrong then. It is only allowed,
    !> If the result updates indices of y which are never used in x.
    !> NOTE: In general, external libraries will not allow to pass the same
    !> vector to x and y.
    !--------------------------------------------------------------------------
    subroutine usmv_xy(a,mat,x,y, ierr)
      !use control, only : root_and_not_silent
#ifdef HAVE_LIBRSB
      !include 'blas_sparse.fi' It is a module, can't put it here
      ! UPPERCASE USE to deliberately exclude from mkdeps script
      USE blas_sparse, only : usmv_librsb => usmv, blas_no_trans
#endif
#ifdef HAVE_MKL
      !lowercase 'use', to include mkl_sblas module as an additional
      ! dependency.
      USE mkl_spblas
#endif
      complex, intent(in) :: a
      type(sparse_matrix), intent(in) :: mat
      complex, intent(in) :: x(:)
      complex, intent(inout) :: y(:)
      integer, intent(out) :: ierr
#ifdef HAVE_LIBRSB
      integer :: incX, incB
#endif
      integer :: i, ir
#ifdef CRS_FASTDIAGONAL_FORMAT
      integer :: idiag_block, idiag_elem_start
      integer :: diag_block_row_end
      integer :: offset
#endif
      complex :: acc
      ierr = 0

      ! if(root_and_not_silent) then
      !   !for debugging
      !   write (*,*) 'multiplication with ', mat%name, mat%nrows
      ! end if

      if(mat%current_format == FORMAT_LIBRSB) then
#ifdef HAVE_LIBRSB
        incX = 1
        incB = 1
        call usmv_librsb(blas_no_trans,a,mat%librsb_handle,x,incX,y,incB,ierr)
#else
        ierr = 1
#endif
      else if(mat%current_format == FORMAT_MKLCRS) then
#if HAVE_MKL
        ! Compute y = alpha * A * x + beta * y
        ierr = mkl_sparse_z_mv(sparse_operation_non_transpose,a,mat%mkl_handle,&
           & mat%mkl_desc,x,cmplx(1.0,0.0),y)
#endif
      else if(mat%current_format == FORMAT_GKWCRS) then
        ! CSR format version (faster: cache access and multi-threading
        ! compared to most simplistic coordinate format) the multiplication
        ! is off diagonal and jj < ii for all i=n2+1, n3 this must be
        ! the case to avoid a datarace

        ! this computation is wrong, if x and y refer to the same memory

        !$omp parallel 
#ifdef CRS_FASTDIAGONAL_FORMAT
        do idiag_block = 1, mat%ndiag_blocks, 2
          ! the first block is a purely diagonal block
          idiag_elem_start = mat%irs(mat%diag_block_row_start(idiag_block))
          diag_block_row_end = mat%diag_block_row_start(idiag_block+1)-1
          offset = idiag_elem_start - 1 - mat%diag_block_row_start(idiag_block) + 1
          !$omp do private(ir,i,acc) schedule(static)
          do ir = mat%diag_block_row_start(idiag_block), diag_block_row_end
            y(ir) = y(ir) + a*mat%mat(offset + ir)*x(ir)
          end do
          !$omp end do
          diag_block_row_end = mat%diag_block_row_start(idiag_block+2)-1
#endif
        !$omp do private(ir,i,acc) schedule(static)
#ifdef CRS_FASTDIAGONAL_FORMAT
          do ir=mat%diag_block_row_start(idiag_block+1), diag_block_row_end
#else
        do ir=mat%start_row, mat%end_row
#endif
          acc = 0.0
          !loop over elements in row:
          do i = mat%irs(ir),mat%irs(ir+1)-1
            acc = acc + a*mat%mat(i)*x(mat%jj(i))
          end do
          y(ir) = y(ir) + acc
        end do
        !$omp end do
#ifdef CRS_FASTDIAGONAL_FORMAT
      end do
#endif
        !$omp end parallel
      else if(mat%current_format == FORMAT_GKWCOO) then
        ! this computation is wrong, if x and y refer to the same memory
        do i = 1,mat%nmat
          y(mat%ii(i)) = y(mat%ii(i)) + a*mat%mat(i)*x(mat%jj(i))
        end do
      else
        ierr = 1
      end if

    end subroutine usmv_xy

    !--------------------------------------------------------------------------
    !> This routine initialises the arrays of matrix mat with the
    !> values given in the other two arguments. The initialisation is
    !> done in parallel and the access pattern should be the same
    !> as in usmv_xy() .
    !--------------------------------------------------------------------------
    subroutine crs_first_touch(mat, mat_, jj_)
      type(sparse_matrix), intent(inout) :: mat
      complex, intent(in) :: mat_(:)
      integer, intent(in) :: jj_(:)
      integer :: ir, i
#ifdef CRS_FASTDIAGONAL_FORMAT
      integer :: idiag_elem_start, diag_block_row_end
      integer :: idiag_block
#endif
      !$omp parallel
#ifdef CRS_FASTDIAGONAL_FORMAT
      do idiag_block = 1, mat%ndiag_blocks, 2
        ! the first block is a purely diagonal block
        idiag_elem_start = mat%irs(mat%diag_block_row_start(idiag_block))
        diag_block_row_end = mat%diag_block_row_start(idiag_block+1)-1
        !$omp do private(ir) schedule(static)
        do ir = mat%diag_block_row_start(idiag_block), diag_block_row_end
          mat%mat(idiag_elem_start-1 + ir) = mat_(idiag_elem_start-1 + ir)
        end do
        !$omp end do
        diag_block_row_end = mat%diag_block_row_start(idiag_block+2)-1
#endif
      !$omp do private(ir,i) schedule(static)
#ifdef CRS_FASTDIAGONAL_FORMAT
        do ir=mat%diag_block_row_start(idiag_block+1), diag_block_row_end
#else
          do ir=mat%start_row, mat%end_row
#endif
        !loop over elements in row:
        do i = mat%irs(ir),mat%irs(ir+1)-1
          mat%mat(i) = mat_(i)
          mat%jj(i) = jj_(i)
        end do
      end do
      !$omp end do
#ifdef CRS_FASTDIAGONAL_FORMAT
    end do
#endif
      !$omp end parallel
     end subroutine crs_first_touch

    !--------------------------------------------------------------------------
    !> compute y = a * mat * y
    !> with an implicit unity matrix in
    !> the matrix block 1:mat%start_row-1 . This means, if this function is
    !> invoked with the fdis array and there are only matrix elements
    !> acting on the fields part, then the distribution part will be
    !> unmodified (instead of nullified)!
    !>
    !> NOTE: it is dangerous to pass the same vector to x and y.
    !> In General the result will be wrong then. It is only allowed,
    !> If the result updates indices of y which are never used in x.
    !--------------------------------------------------------------------------
    subroutine usmv_yupdate(a,mat,y,ierr)
      use general, only : gkw_abort
#ifdef HAVE_LIBRSB
      !include 'blas_sparse.fi' It is a module, can't put it here
      ! UPPERCASE USE to deliberately exclude from mkdeps script
      USE blas_sparse, only : usmv_librsb => usmv, blas_no_trans
#endif
      complex, intent(in) :: a
      type(sparse_matrix), intent(in) :: mat
      complex, intent(inout) :: y(:)
      integer, intent(out) :: ierr
#ifdef HAVE_LIBRSB
      integer :: incX, incB
#endif
      integer :: i, ir
      complex :: acc

      if(mat%current_format == FORMAT_LIBRSB) then
#ifdef HAVE_LIBRSB
        incX = 1
        incB = 1
        call gkw_abort('pick the routine computing y=a*mat*y from sparse lib &
           &and remove this.')
        !There is no such routine. There is only the usmv_xy
        !variant. One has to make the matrix such that it does the
        !same thing with usmv_xy.
#else
        ierr = 1
#endif
      else if(mat%current_format == FORMAT_GKWCRS) then

        ! For single-threaded execution, the algorithm below yields
        ! the correct result, if the matrix is a lower triangular
        ! matrix: the condition jj <= ii for all matrix
        ! elements must be fulfilled.

        ! For thread-parallel execution this is generally a wrong
        ! algorithm and causes a datarace for any matrix.  This is
        ! because elements of y are updated, and subsequently used
        ! to update other elements of y...

        ! The datarace can be avoided by any of these methods:
        !   * making sure that the matrix is purely diagonal
        !   * introducing a local buffer vector
        !   * rewriting the calling code to use the usmv_xy() routine.

        if(mat%is_diagonal) then
          !$omp parallel private(ir,i,acc)
          !$omp do
          do ir=mat%start_row, mat%end_row
            acc = 0.0
            !loop over elements in row:
            do i = mat%irs(ir),mat%irs(ir+1)-1
              acc = acc + a*mat%mat(i)*y(mat%jj(i))
            end do
            y(ir) = acc
          end do
          !$omp end do
          !$omp end parallel
        else
          !$omp master
          do ir=mat%start_row, mat%end_row
            acc = 0.0
            !loop over elements in row:
            do i = mat%irs(ir),mat%irs(ir+1)-1
              acc = acc + a*mat%mat(i)*y(mat%jj(i))
            end do
            y(ir) = acc
          end do
          !$omp end master
        end if
        ierr = 0
      else
        ierr = 1
        call gkw_abort('usmv_yupdate is not implemented for the chosen matrix &
           &format.')
      end if

    end subroutine usmv_yupdate

    


    !------------------------------------------------------------------------------
    !> This function returns the matrix row (ii) value of a given element 
    !------------------------------------------------------------------------------
    function get_elem_ii(E)

      use dist,           only : iapar, iphi, ibpar, i_mom, i_ene, ifdis 
      use dist,           only : iphi_ga, iapar_ga, ibpar_ga 
      use index_function, only : indx
      use general,        only : gkw_abort
      use structures,     only : matrix_element

      type (matrix_element), intent(in) :: E

      integer :: get_elem_ii

      ! Calculate iih 
      select case(E%itype)
      case(iphi,iapar,ibpar)
        get_elem_ii = indx(E%itype,E%imod,E%ix,E%i)
      case(iphi_ga,iapar_ga,ibpar_ga) 
        get_elem_ii = indx(E%itype,E%imod,E%ix,E%i,E%j,E%is)
      case(i_mom,i_ene) 
        get_elem_ii = indx(E%itype,E%imod,E%ix,E%i,E%is)
      case(ifdis)
        get_elem_ii = indx(E%itype,E%imod,E%ix,E%i,E%j,E%k,E%is)
      case default
        get_elem_ii = 0
        call gkw_abort('get_elem_ii: undefined itype in term '//E%term)
      end select

    end function get_elem_ii


    !---------------------------------------------------------------------------
    !> This function returns the matrix column (jj) value of a given element
    !---------------------------------------------------------------------------
    function get_elem_jj(E)
      use dist,           only : iapar, iphi, ibpar, i_mom, i_ene, ifdis
      use dist,           only : iphi_ga, iapar_ga, ibpar_ga
      use index_function, only : indx
      use general,        only : gkw_abort
      use structures,     only : matrix_element
      type (matrix_element), intent(in) :: E
      integer :: get_elem_jj

      select case(E%itloc)
      case(iphi,iapar,ibpar)
        get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc)
      case(iphi_ga,iapar_ga,ibpar_ga)
        get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc,E%jloc,E%isloc)
      case(i_mom,i_ene)
        get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc,E%isloc)
      case(ifdis)
        get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc,E%jloc,E%kloc,E%isloc)
      case default
        get_elem_jj = 0
        call gkw_abort('get_elem_jj: undefined itloc in term '//E%term)
      end select

    end function get_elem_jj

    !-----------------------------------------------------------------------------
    !> Compress complex matrix piece between istart and iend.
    !> Same as compress_piece_real, matched via interface to compress_piece.
    !>
    !> iend returns the updated compressed location of the segment end.
    !> Ranges must always be compressed sequentially to avoid gaps and junk data.
    !> All of the matrix below istart must already have been compressed
    !> It is acceptable however to recompress a range with a new piece at the end.
    !>
    !> If sections are compressed together after the per-section compression 
    !> this results in a re-sorting but not a size reduction because there is no 
    !> overlap between the four matrix sections n1,n2,n3,n4
    !-----------------------------------------------------------------------------
    subroutine compress_matrix(mat, on_columns)
      use control, only : root_and_not_silent
      type(sparse_matrix), intent(inout):: mat
      logical, optional, intent(in) :: on_columns
      integer :: i, ireduced

      if(mat%nmat <= 1) then
        return
      end if
      if(present(on_columns)) then
        if(on_columns) then
          call swap(mat%ii,mat%jj)
        end if
      end if

      ! if (root_and_not_silent) then
      !   write(*,*) "compressing ", trim(mat%name), "..."
      ! end if

      ! first sort the matrix
      call sort_matrix(mat)

      ireduced = 1
      ! compress the matrix (not robust if the sort was incorrect)
      do i = 2, mat%nmat
        if (mat%ii(i) == mat%ii(ireduced) .and. &
           & mat%jj(i) == mat%jj(ireduced)) then
          mat%mat(ireduced) = mat%mat(ireduced) + mat%mat(i)
        else
          ireduced = ireduced + 1
          mat%ii(ireduced) = mat%ii(i)
          mat%jj(ireduced) = mat%jj(i)
          mat%mat(ireduced) = mat%mat(i)
        end if
      end do

      if (root_and_not_silent) then
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*) 'Matrix compression successfully completed'
        write(*,*) '(',trim(mat%name),')'
        write(*,223) mat%nmat,ireduced,size(mat%mat)
223     format(' Before ',I8,' elements. New ',I8,' elements. Max ', I8,' elements')
224     format(' first row ',I8, ',last row ', I8,', in total ',I8,' rows. Full compression #: ',I8)
        ! Count the full section compressions
        mat%compressions = mat%compressions + 1
        write(*,224) mat%start_row, mat%end_row, mat%nrows, mat%compressions
        if(mat%nrejected_elements > 0) then
          write(*,*)'rejected elements:',mat%nrejected_elements
        end if
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)
      end if

      mat%nmat=ireduced

      if(present(on_columns)) then
        if(on_columns) then
          call swap(mat%ii,mat%jj)
        end if
      end if

    contains
      pure subroutine swap(a, b)
        integer, dimension(:),intent(inout) :: a, b
        integer, dimension(size(a)) :: work
        work = a
        a = b
        b = work

        ! ! XOR swap of column and row indices
        ! a = ieor(a,b)
        ! b = ieor(b,a)
        ! a = ieor(a,b)

      end subroutine swap
    end subroutine compress_matrix

    

    !-----------------------------------------------------------------------------
    !> Sort the matrix on ii and then jj. Use in-place recursive quicksort, then
    !> switch to insertion sort when the partitions are of size min_size or less.
    !> (min_size = 9 was found to be optimal in most cases tested)
    !-----------------------------------------------------------------------------
    subroutine sort_matrix(mat)
      type(sparse_matrix), intent(inout) :: mat

      call qsort(mat,1,mat%nmat)

    contains

      !---------------------------------------------------------------------
      !> Quicksort
      !---------------------------------------------------------------------
      recursive subroutine qsort(mat,b,e)
        type(sparse_matrix), intent(inout) :: mat
        integer, intent(in) :: b,e
        integer :: ind
        !> The smallest sized partition on which to perform quicksort
        integer, parameter :: min_size = 9

        if (e - b > min_size) then
          ind = median_of_3_ind(mat,b,e)
          call partition(mat,b,e,ind)
          call qsort(mat,b,ind-1)
          call qsort(mat,ind+1,e)
        else
          call isort(mat,b,e)
        end if

      end subroutine qsort
      !---------------------------------------------------------------------
      subroutine partition(mat,b,e,ind)
        type(sparse_matrix), intent(inout) :: mat
        integer, intent(in) :: b, e
        integer, intent(inout) :: ind
        integer :: ival, jval, i

        ! the pivot values
        ival = mat%ii(ind)
        jval = mat%jj(ind)

        ! swap the pivot point with the end point
        call swap_elements(mat,ind,e)

        ! put elements to the right or left of the pivot value
        ind = b
        do i = b, e-1
          if (lesseq(mat%ii(i),mat%jj(i),ival,jval)) then
            call swap_elements(mat,i,ind)
            ind = ind + 1
          end if
        end do

        ! put the pivot value back at new ind
        call swap_elements(mat,ind,e)

      end subroutine partition
      !---------------------------------------------------------------------------
      subroutine swap_elements(mat,i,j)
        type(sparse_matrix), intent(inout) :: mat
        integer, intent(in) :: i, j
        integer :: itmp
        complex :: ctmp

        itmp = mat%ii(i)
        mat%ii(i) = mat%ii(j)
        mat%ii(j) = itmp

        itmp = mat%jj(i)
        mat%jj(i) = mat%jj(j)
        mat%jj(j) = itmp

        ctmp = mat%mat(i)
        mat%mat(i) = mat%mat(j)
        mat%mat(j) = ctmp
      end subroutine swap_elements
      !---------------------------------------------------------------------------
      function median_of_3_ind(mat,a,c)
        type(sparse_matrix), intent(in) :: mat
        integer, intent(in) :: a, c
        integer :: b, median_of_3_ind

        b = (a + c) / 2
        if (larger(mat%ii(a),mat%jj(a),mat%ii(b),mat%jj(b))) then
          if (larger(mat%ii(c),mat%jj(c),mat%ii(a),mat%jj(a))) then
            median_of_3_ind = a
          else if (larger(mat%ii(c),mat%jj(c),mat%ii(b),mat%jj(b))) then
            median_of_3_ind = c
          else
            median_of_3_ind = b
          end if
        else if (larger(mat%ii(b),mat%jj(b),mat%ii(c),mat%jj(c))) then
          if (larger(mat%ii(c),mat%jj(c),mat%ii(a),mat%jj(a))) then
            median_of_3_ind = c
          else
            median_of_3_ind = a
          end if
        else
          median_of_3_ind = b
        end if

      end function median_of_3_ind
      !-------------------------------------------------------------------------
      !> Insertion sort
      !-------------------------------------------------------------------------
      subroutine isort(mat,b,e)
        type(sparse_matrix), intent(inout) :: mat
        integer, intent(in) :: b, e
        integer :: ival, jval, i, j
        complex :: val = (0.0, 0.0)

        do i = b+1, e
          ival = mat%ii(i)
          jval = mat%jj(i)
          val  = mat%mat(i)
          
          ajgtval : do j = i-1, 1, -1
            if (less(ival,jval,mat%ii(j),mat%jj(j))) then
              mat%mat(j+1) = mat%mat(j)
              mat%ii(j+1)  = mat%ii(j)
              mat%jj(j+1)  = mat%jj(j)
            else
              exit ajgtval
            end if
          end do ajgtval
          mat%mat(j+1) = val
          mat%ii(j+1)  = ival
          mat%jj(j+1)  = jval
        end do

      end subroutine isort
      !---------------------------------------------------------------------------
      pure function larger(ii1,jj1,ii2,jj2)
        integer, intent(in) :: ii1, ii2, jj1, jj2
        logical :: larger
        larger = (ii1 > ii2) .or. (ii1 == ii2 .and. jj1 > jj2)
      end function larger
      !---------------------------------------------------------------------------
      pure function lesseq(ii1,jj1,ii2,jj2)
        integer, intent(in) :: ii1, ii2, jj1, jj2
        logical :: lesseq
        lesseq =  (ii1 < ii2) .or. (ii1 == ii2 .and. jj1 <= jj2)
      end function lesseq
      !---------------------------------------------------------------------------
      pure function less(ii1,jj1,ii2,jj2)
        integer, intent(in) :: ii1, ii2, jj1, jj2
        logical :: less
        less = (ii1 < ii2) .or. (ii1 == ii2 .and. jj1 < jj2)
      end function less

    end subroutine sort_matrix

    !-----------------------------------------------------------------------------
    !> This routine finds the starting element of each row in the matrix, and
    !> stores the index it in irs(:).
    !> The matrix must be fully sorted and compressed.
    !> The indices found are used for faster matrix-vector multiply, with
    !> the explicit scheme.
    !>
    !> The nth element of irs is the index of first element of the n'th
    !> matrix row with respect to the total number of elements of the
    !> matrix. The elements of the n'th matrix row are then
    !> mat(irs(n):(irs(n+1)-1)) .
    !-----------------------------------------------------------------------------
    subroutine crs_find_row_index(mat)
      type(sparse_matrix), intent(inout) :: mat
      integer :: i, current_row

#ifdef _OPENMP
#ifdef FIRST_TOUCH
      ! parallel first touch of irs to achieve good data locality.
      !$omp parallel do private(i) schedule(static)
      do i=mat%start_row, mat%end_row
        mat%irs(i) = 0
      end do
      !$omp end parallel do
#endif
#endif

      current_row = 0
      ! iterate over all matrix elements
      do i = 1, mat%nmat

        if (current_row == mat%ii(i)) then
          ! if that element belongs to the row indicated by the
          ! current_row counter, then just proceed
        else if (mat%ii(i) > current_row ) then
          ! if that element belongs to the next row, then
          do
            ! increment the counter..
            current_row = current_row + 1
            if (mat%ii(i) == current_row) then
              ! until it matches the row of the element
              mat%irs(current_row) = i
              exit
            else if(current_row >= mat%start_row) then
              ! produces an empty row
              mat%irs(current_row) = i
              ! this is needed in some cases (e.g. n_s_grid = 1, k=0 mode)
            end if
          end do
        end if
      end do

      ! Extra line marking row beyond end of the matrix for convenience
      mat%irs(current_row+1) = mat%nmat+1

    end subroutine crs_find_row_index

    !--------------------------------------------------------------------------
    !> usmm means:
    !> unstructured sparse matrix matrix multiplication.
    !>
    !> C = a A B + C
    !>
    !--------------------------------------------------------------------------
    subroutine usmm(a,matA,matB,matC,ierr)
      use general, only : gkw_abort
      type(sparse_matrix), intent(in) :: matA
      type(sparse_matrix), intent(in) :: matB
      type(sparse_matrix), intent(inout) :: matC
      complex, intent(in) :: a
      integer, intent(out) :: ierr
      integer :: irA, irB, irC, iA, iB, iC

      if(matC%is_finished) then
        ! programming error
        call gkw_abort('usmm cannot be used, because '//trim(matC%name)//&
           & ' has already been finished.')
      end if

      if(matA%count_only .and. matB%count_only) then
        ! this is crude!
        matC%nmat = matA%nrows*matB%ncols
      else

        if(matA%current_format == FORMAT_GKWCRS .and. &
           & matB%current_format == FORMAT_GKWCRS) then
          irC = 1
          iC = 1
          do irA = matA%start_row, matA%end_row;
            ! loop over elements in row of A:
            do iA = matA%irs(irA), matA%irs(irA+1)-1

              do irB = matB%start_row, matB%end_row;
                ! loop over elements in row of B:
                do iB = matB%irs(irB), matB%irs(irB+1)-1
                  if(matA%jj(iA) == irB) then
                    ! call put_element_simple(matC, irA, matB%jj(iB), &
                    !    & a * matA%mat(iA) * matB%mat(iB))
                    call put_element_simple(matC, irA, matB%jj(iB), &
                       & a * matA%mat(iA) * matB%mat(iB))
                  end if
                end do
              end do

            end do
          end do
        end if
      end if

    end subroutine usmm

end module matrix_format
