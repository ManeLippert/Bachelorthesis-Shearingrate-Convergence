!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Provide an interface to the mpi parameters and functions.
!> The "use mpi" fortran 90 interface is less portable than the mpif.h header
!> but sometimes gives better error reporting (see issue 140) 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module mpiinterface

#if defined(mpif90_interface)
  ! UPPERCASE USE to deliberately exclude from mkdeps script
  USE mpi
#endif
#if !defined(mpi2)
  use global, only : I_RUN_WITHOUT_MPI
#endif 

  implicit none

  private :: get_basic_mpi_params, ierr
  
  interface mpibarrier
    module procedure mpibarrier_comm
  end interface mpibarrier

  interface mpibcast
    module procedure mpibcast_logical_scalar
    module procedure mpibcast_real_scalar
    module procedure mpibcast_cmplx_scalar
    module procedure mpibcast_real_array1
    module procedure mpibcast_real_array2
    module procedure mpibcast_real_array3
    module procedure mpibcast_real_array4
    module procedure mpibcast_real_array5
    module procedure mpibcast_integer_scalar
    module procedure mpibcast_integer_array1
    module procedure mpibcast_character_scalar
  end interface mpibcast

  interface mpireduce_sum_inplace
    module procedure mpireduce_sum_r_inpl
    module procedure mpireduce_sum_r_array_1_inpl
    module procedure mpireduce_sum_r_array_2_inpl
    module procedure mpireduce_sum_r_array_3_inpl
    module procedure mpireduce_sum_r_array_4_inpl
    module procedure mpireduce_sum_c_array_1_inpl
    module procedure mpireduce_sum_c_array_2_inpl
    module procedure mpireduce_sum_c_array_3_inpl
    module procedure mpireduce_sum_c_array_4_inpl
    module procedure mpireduce_sum_c_array_6_inpl
  end interface mpireduce_sum_inplace
  
  interface mpireduce_sum
    module procedure mpireduce_sum_r_array_2
    module procedure mpireduce_sum_r_array_3
    module procedure mpireduce_sum_c_array_2
    module procedure mpireduce_sum_c_array_4
  end interface mpireduce_sum

  interface mpiallreduce_sum
    module procedure mpiallreduce_sum_r_array_1
! used by imp integration
#if defined(real_precision_default)
    module procedure mpiallreduce_sum_d_array_1
#endif
    module procedure mpiallreduce_sum_r_array_2
    module procedure mpiallreduce_sum_r_array_3
    module procedure mpiallreduce_sum_r_array_4
    module procedure mpiallreduce_sum_c_scalar
    module procedure mpiallreduce_sum_r_scalar
    module procedure mpiallreduce_sum_int_scalar
    module procedure mpiallreduce_sum_c_array_1
    module procedure mpiallreduce_sum_c_array_2
    module procedure mpiallreduce_sum_c_array_3
    module procedure mpiallreduce_sum_c_array_4
    module procedure mpiallreduce_sum_c_array_5
  end interface mpiallreduce_sum
  
  interface mpiallreduce_sum_inplace
    module procedure mpiallreduce_sum_i_inpl
    module procedure mpiallreduce_sum_r_inpl
    module procedure mpiallreduce_sum_c_inpl
    module procedure mpiallreduce_sum_r_array_1_inpl
    module procedure mpiallreduce_sum_r_array_2_inpl
    module procedure mpiallreduce_sum_r_array_3_inpl
    module procedure mpiallreduce_sum_r_array_4_inpl
    module procedure mpiallreduce_sum_c_array_1_inpl
    module procedure mpiallreduce_sum_c_array_2_inpl
    module procedure mpiallreduce_sum_c_array_3_inpl
    module procedure mpiallreduce_sum_c_array_4_inpl
  end interface mpiallreduce_sum_inplace

  interface mpiallreduce_max
    module procedure mpiallreduce_max_int_scalar
    module procedure mpiallreduce_max_int_array_1
    module procedure mpiallreduce_max_real_scalar
  end interface mpiallreduce_max

  interface mpiallreduce_min
    module procedure mpiallreduce_min_int_scalar
    module procedure mpiallreduce_min_r_scalar
  end interface mpiallreduce_min

  interface mpiallreduce_or
    module procedure mpiallreduce_or_scalar
  end interface mpiallreduce_or

  interface mpiallreduce_or_inplace
    module procedure mpiallreduce_or_scalar_inpl
  end interface mpiallreduce_or_inplace

  interface mpiallreduce_and
    module procedure mpiallreduce_and_logical_scalar
  end interface mpiallreduce_and

  interface mpicart_coords_self
    module procedure mpicart_coords_self_1d
  end interface mpicart_coords_self

!APS: Could we have a simpler interface which does not care how many dims your
!APS: arrays are, but instead provide a list of the distributed dimensions,
!APS: together with their communicators (or alternative just a list of id's, so
!APS: that you don't need the communicator even...)????
!APS: So the call would be gather_array(array,  [real, dim(:,...,:)]
!APS:                                   dims,   [integer]
!APS:                     communicated_dimms,   [logical, dim(dims)] { could use zeros to
!APS:                                  comms,   [integer, dim(dims)]   merge these two }
!APS:                            local_array,   [real, dim(:,...,:)]
!APS:                        ALLGATHER=[T|F],   [logical, optional ] )
!APS: (and also one for different datatypes) Would it be easy to do????
!FJC: It might not be easy, but it might make gyro-average / diagnostic a bit neater, 
!FJC: since at present we are looping over lower dimensional gathers
!FJC: Surely someone has already written this kind of thing elsewhere?

  interface gather_array
    module procedure gather_1d_int
    module procedure gather_1d_real
    module procedure gather_1d_real_intercomm
    module procedure gather_1d_cmplx
    module procedure gather_1d_cmplx_intercomm
    module procedure gather_2d_real
    module procedure gather_2d_real_intercomm
    module procedure gather_2d_cmplx
    module procedure gather_2d_int
    module procedure gather_3d_real
    module procedure gather_3d_cmplx
    module procedure gather_4d_real
    module procedure gather_4d_cmplx
    module procedure gather_6d_real
  end interface gather_array

  interface check_equal_on_all_processors
    module procedure chk_eq_string
    module procedure chk_eq_real
    module procedure chk_eq_complex
    module procedure chk_eq_integer
    module procedure chk_eq_logical
    module procedure chk_eq_r_array
    module procedure chk_eq_i_array
  end interface check_equal_on_all_processors

  interface send_to_root
    module procedure send_scalar_int
    module procedure send_scalar_matrix_element
    module procedure send_1d_int
    module procedure send_1d_complex
    module procedure send_2d_real
    module procedure send_4d_complex
    module procedure send_4d_real
  end interface send_to_root


! Swap variables and routines depending on if we use mpi or not.

#if defined(mpi2)
 
#if defined(mpif90_interface)
  ! included at top of module
#else
  ! UPPERCASE INCLUDE to deliberately exclude from mkdeps script
  INCLUDE 'mpif.h'
#endif     
  
  ! lines below stops us accidently using these in the code;
  ! eventually make everything in mpi private by default.
  ! If you need a new usage of an MPI routine not already defined,
  ! please add a new wrapper routine to this module (see issues 56 and
  ! 124).
  private :: MPI_DOUBLE_PRECISION, MPI_REAL, MPI_COMPLEX
  ! Please DO NOT make MPI_COMM_WORLD public, it is private to prevent code
  ! appearing that does not compile without MPI. Usually COMM_CART can be 
  ! used instead of MPI_COMM_WORLD.  
  private :: MPI_COMM_WORLD

  character (len=4), parameter, public :: mpistandard='MPI2'

#else

  ! THESE SHOULD BE made into proper dummies
  !> dummy MPI_STATUS_SIZE
  integer, parameter :: MPI_STATUS_SIZE = 1
  !> dummy MPI_PROC_NULL
  integer, parameter :: MPI_PROC_NULL = -1
  !> default offset kind
  integer, parameter :: MPI_OFFSET_KIND = KIND(1.0D0)

  character (len=4), parameter, public :: mpistandard='NONE'
#endif

  !> total number of processors
  integer, save :: number_of_processors
  !> local processor rank
  integer, save :: processor_number 
  !> for mpi error 
  integer, save :: ierr
  !> root processor
  logical, save :: root_processor
  !> the processor with the largest rank
  logical, save:: last_processor
  !> Is MPI ready to use?
  logical, save :: lmpi_is_ready = .false.
  !> data representation
  character (len=6), parameter :: data_representation = 'native'
  ! status
  integer, dimension(MPI_STATUS_SIZE), save :: statusmpi

  !>
  integer, save, private :: tag_range_upper_limit = 0
  integer, parameter, private :: tag_range_size = 1000


! set MPIREAL_X etc. to some MPI precision determined at compile time

#if defined(mpi2)


#if defined(real_precision_8)
  !> the MPI type corresponding to the default real type 
  integer, parameter :: MPIREAL_X = MPI_REAL8
  integer, parameter :: MPI2REAL_X = MPI_2DOUBLE_PRECISION
  !> the MPI type corresponding to the default complex type used
  integer, parameter :: MPICOMPLEX_X = MPI_COMPLEX16
#endif     

#if defined(real_precision_4)
  integer, parameter :: MPIREAL_X = MPI_REAL4
  integer, parameter :: MPI2REAL_X = MPI_2REAL
  integer, parameter :: MPICOMPLEX_X = MPI_COMPLEX8
#endif

#if defined(real_precision_double)
  integer, parameter :: MPIREAL_X = MPI_DOUBLE_PRECISION
  integer, parameter :: MPI2REAL_X = MPI_2DOUBLE_PRECISION
  integer, parameter :: MPICOMPLEX_X = MPI_DOUBLE_COMPLEX
#endif

#if defined(real_precision_default)
  integer, parameter :: MPIREAL_X = MPI_REAL
  integer, parameter :: MPI2REAL_X = MPI_2REAL
  integer, parameter :: MPICOMPLEX_X = MPI_COMPLEX
#endif

#else
  integer, parameter :: MPI_MODE_WRONLY   = I_RUN_WITHOUT_MPI
  integer, parameter :: MPI_MODE_RDONLY   = I_RUN_WITHOUT_MPI
  integer, parameter :: MPI_MODE_CREATE   = I_RUN_WITHOUT_MPI
  integer, parameter :: MPI_MODE_APPEND   = I_RUN_WITHOUT_MPI
  integer, parameter :: MPI_INFO_NULL     = I_RUN_WITHOUT_MPI
  integer, parameter :: MPI_INTEGER       = I_RUN_WITHOUT_MPI
  integer, parameter :: mpi_address_kind  = I_RUN_WITHOUT_MPI
  integer, parameter :: MPICOMPLEX_X      = I_RUN_WITHOUT_MPI
  integer, parameter :: MPIREAL_X         = I_RUN_WITHOUT_MPI 
  integer, parameter :: MPI2REAL_X        = I_RUN_WITHOUT_MPI
  integer, parameter :: MPI_COMM_SELF     = I_RUN_WITHOUT_MPI 
#endif


  !> a MPI datatype which is used to communicate complete matrix elements
  integer, private, save :: matrix_element_type


contains


!-----------------------------------------------------------------------------
!> The basic mpi parameters are obtained here. This provides the number of
!> processors and the processor number to the rest of the code.
!-----------------------------------------------------------------------------
subroutine get_basic_mpi_params()

#if defined(mpi2)

  call mpi_comm_size(MPI_COMM_WORLD, number_of_processors, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, processor_number, ierr)

#else

  number_of_processors = 1
  processor_number     = 0

#endif 

end subroutine get_basic_mpi_params

!-----------------------------------------------------------------------------
!> find the coordinate in one direction
!-----------------------------------------------------------------------------
subroutine mpicart_coords_self_1d(ICOMM,coord)

  integer, intent(in)  :: ICOMM
  integer, intent(out) :: coord(1)

#if defined(mpi2)

  integer            :: rank,comm
  integer, parameter :: maxdims = 1

  comm=ICOMM
  if (number_of_processors > 1) then
    call mpicomm_rank(comm,rank)
    call mpi_cart_coords(comm,rank,maxdims,coord,ierr)
  else
    coord = 0
  end if

#else

  coord = 0

#endif

end subroutine mpicart_coords_self_1d

!-----------------------------------------------------------------------------
!> mpi_wtime() substitute if no mpi implementation is available
!-----------------------------------------------------------------------------
function mpiwtime()

  use global, only : dp

  real (kind=dp) :: mpiwtime
#if defined(mpi2)
  mpiwtime = MPI_WTIME()
#else
  integer :: ccount, ccount_rate

  call system_clock(ccount, ccount_rate)

  if (ccount_rate == 0) then
    mpiwtime = 0.
  else
    mpiwtime = (1.*ccount) / ccount_rate
  end if
#endif

end function mpiwtime


!-----------------------------------------------------------------------------
!> mpi_wtick() substitute if no mpi implementation is available
!-----------------------------------------------------------------------------
function mpiwtick()

  use global, only : dp

  real (kind=dp) :: mpiwtick
#if defined(mpi2)
  mpiwtick = MPI_WTICK()
#else
  integer :: ccount_rate

  call system_clock(COUNT_RATE=ccount_rate)

  if (ccount_rate == 0) then
    mpiwtick = 0.
  else
    mpiwtick = 1./ccount_rate
  end if
#endif

end function mpiwtick

!-----------------------------------------------------------------------------
!> match LOGICAL SCALAR for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_logical_scalar(lscalar,nscalar,PROC)
 
  use global, only : logical_false

  logical, intent(inout)         :: lscalar
  integer, intent(in)            :: nscalar
  integer, optional,  intent(in) :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(lscalar,nscalar,MPI_LOGICAL,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) lscalar,nscalar,PROC

#endif

end subroutine mpibcast_logical_scalar

!-----------------------------------------------------------------------------
!> match INTEGER SCALAR for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_integer_scalar(iscalar,nscalar,PROC)

  use global, only : logical_false

  integer, intent(inout)         :: iscalar
  integer, intent(in)            :: nscalar
  integer, optional,  intent(in) :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(iscalar,nscalar,MPI_INTEGER,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) iscalar,nscalar,PROC

#endif

end subroutine mpibcast_integer_scalar

!-----------------------------------------------------------------------------
!> match INTEGER ARRAY for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_integer_array1(rscalar,nscalar,PROC)

  use global, only : logical_false

  integer, dimension(:), intent(inout) :: rscalar
  integer, intent(in)                  :: nscalar
  integer, optional,  intent(in)       :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(rscalar,nscalar,MPI_INTEGER,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC

#endif

end subroutine mpibcast_integer_array1

!-----------------------------------------------------------------------------
!> match REAL SCALAR for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_real_scalar(rscalar,nscalar,PROC)

  use global, only : logical_false

  real, intent(inout)            :: rscalar
  integer, intent(in)            :: nscalar
  integer, optional,  intent(in) :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(rscalar,nscalar,MPIREAL_X,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC

#endif

end subroutine mpibcast_real_scalar

!-----------------------------------------------------------------------------
!> match COMPLEX SCALAR for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_cmplx_scalar(rscalar,nscalar,PROC, comm)

  use global, only : logical_false

  complex, intent(inout)         :: rscalar
  integer, intent(in)            :: nscalar
  integer, optional,  intent(in) :: PROC
  integer, optional,  intent(in) :: comm

#if defined(mpi2)

  integer :: iproc_bcast
  integer :: icom

  ! keep compiler quiet
  if (logical_false) continue

  icom = MPI_COMM_WORLD
  if (present(comm)) icom = comm
  
  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  
  call MPI_BCAST(rscalar,nscalar,MPICOMPLEX_X,iproc_bcast,icom,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC

#endif

end subroutine mpibcast_cmplx_scalar

!-----------------------------------------------------------------------------
!> match REAL ARRAY for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_real_array1(rscalar,nscalar,PROC)

  use global, only : logical_false

  real, dimension(:), intent(inout) :: rscalar
  integer, intent(in)               :: nscalar
  integer, optional,  intent(in)    :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(rscalar,nscalar,MPIREAL_X,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC
 
#endif

end subroutine mpibcast_real_array1

!-----------------------------------------------------------------------------
!> match REAL ARRAY for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_real_array2(rscalar,nscalar,PROC)

  use global, only : logical_false

  real, dimension(:,:), intent(inout) :: rscalar
  integer, intent(in)                 :: nscalar
  integer, optional,  intent(in)      :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(rscalar,nscalar,MPIREAL_X,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC
 
#endif

end subroutine mpibcast_real_array2

!-----------------------------------------------------------------------------
!> match REAL ARRAY for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_real_array3(rscalar,nscalar,PROC)

  use global, only : logical_false

  real, dimension(:,:,:), intent(inout) :: rscalar
  integer, intent(in)                   :: nscalar
  integer, optional,  intent(in)        :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(rscalar,nscalar,MPIREAL_X,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC

#endif

end subroutine mpibcast_real_array3

!-----------------------------------------------------------------------------
!> match REAL ARRAY for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_real_array4(rscalar,nscalar,PROC)

  use global, only : logical_false

  real, dimension(:,:,:,:), intent(inout) :: rscalar
  integer, intent(in)                     :: nscalar
  integer, optional,  intent(in)          :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(rscalar,nscalar,MPIREAL_X,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC

#endif

end subroutine mpibcast_real_array4

!-----------------------------------------------------------------------------
!> match REAL ARRAY for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_real_array5(rscalar,nscalar,PROC)

  use global, only : logical_false

  real, dimension(:,:,:,:,:), intent(inout) :: rscalar
  integer, intent(in)                     :: nscalar
  integer, optional,  intent(in)          :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(rscalar,nscalar,MPIREAL_X,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) rscalar,nscalar,PROC

#endif

end subroutine mpibcast_real_array5

!-----------------------------------------------------------------------------
!> match CHARACTER SCALAR for MPI_BCAST
!----------------------------------------------------------------------------
subroutine mpibcast_character_scalar(cscalar,nscalar,PROC)

  use global, only : logical_false

  character (len=*), intent(inout) :: cscalar
  integer, intent(in)              :: nscalar
  integer, optional,  intent(in)   :: PROC

#if defined(mpi2)

  integer :: iproc_bcast

  ! keep compiler quiet
  if (logical_false) continue

  iproc_bcast = 0
  if (present(PROC)) iproc_bcast = PROC
  call MPI_BCAST(cscalar,nscalar,MPI_CHARACTER,iproc_bcast,MPI_COMM_WORLD,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write(*,*) cscalar,nscalar,PROC

#endif

end subroutine mpibcast_character_scalar

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpifinalize()

  call mpi_type_free(matrix_element_type, ierr)

#if defined(mpi2)
  call mpibarrier()
  call mpi_finalize(ierr)
#endif

  lmpi_is_ready = .false.

end subroutine mpifinalize

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiinit()

#if defined(mpi2)

#ifdef _OPENMP
  integer :: support
  !$omp parallel
  !$omp master
  !If the processor is multithreaded, only the thread that called
  !mpi_init_thread() will make MPI calls.
  call mpi_init_thread(MPI_THREAD_FUNNELED, support, ierr)
  !$omp end master
  !$omp end parallel
#else
  call mpi_init(ierr)
#endif

#endif

  ! set processor_number and number_of_processors
  call get_basic_mpi_params

  ! set the root_processor and last_processor logicals
  root_processor = ( processor_number == 0 )
  last_processor = ( processor_number == number_of_processors - 1 )
  if (root_processor) then
    write(*,*)
    write(*,*)'Running with ',number_of_processors,' MPI processes'
#ifdef _OPENMP
    write(*,*) 'MPI has been initialised for multithreading.', MPI_THREAD_FUNNELED, &
       & support
#endif
    write(*,*)
  endif

  lmpi_is_ready = .true.


  call create_matrix_element_type()
  
end subroutine mpiinit

!-----------------------------------------------------------------------------
!> call mpi_barrier() with communicator icomm, or just MPI_COMM_WORLD
!-----------------------------------------------------------------------------
subroutine mpibarrier_comm(COMM)
  
  use global, only : logical_false

  integer, optional, intent(in) :: COMM

#if defined(mpi2)
  integer :: icom
  icom = MPI_COMM_WORLD
  if (present(COMM)) icom = COMM
  call MPI_BARRIER(icom,ierr)
#endif

  ! keep compiler quiet
  if (logical_false) write(*,*) COMM

end subroutine mpibarrier_comm

!-----------------------------------------------------------------------------
!> call MPI_COMM_RANK
!-----------------------------------------------------------------------------
subroutine mpicomm_rank(icomm,rank)

  use global, only : logical_false

  integer, intent(in)  :: icomm
  integer, intent(out) :: rank

#if defined(mpi2)
  ! keep compiler quiet
  if (logical_false) continue

  call MPI_COMM_RANK(icomm,rank,ierr)
#else
  rank = 0
  ! keep compiler quiet
  if (logical_false) write (*,*) icomm
#endif

end subroutine mpicomm_rank


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_r_inpl(a, COMM)
  real, intent(inout)  :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  call mpi_comm_rank(i_comm, rank, ierr)
  
  ! sum the elements of all processes and store the resulting sum only in
  ! the root process.
  if (rank == root_proc_rank) then
    call mpi_reduce(MPI_IN_PLACE, a, 1, MPIREAL_X, MPI_SUM, &
       & root_proc_rank, i_comm, ierr)
  else
    call mpi_reduce(a, MPI_IN_PLACE, 1, MPIREAL_X, MPI_SUM, &
       & root_proc_rank, i_comm, ierr)
  end if
#else

  continue

#endif

end subroutine mpireduce_sum_r_inpl


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_r_array_1_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  real, dimension(dims(1)), intent(inout)  :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
  
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if (rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else
  
  continue

#endif

end subroutine mpireduce_sum_r_array_1_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_c_array_1_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(dims(1)), intent(inout)  :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else
  
  continue

#endif

end subroutine mpireduce_sum_c_array_1_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_c_array_2_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(dims(1), dims(2)), intent(inout) :: a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else
  
  continue

#endif

end subroutine mpireduce_sum_c_array_2_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_c_array_3_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(dims(1), dims(2), dims(3)), intent(inout) :: a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else
  
  continue

#endif

end subroutine mpireduce_sum_c_array_3_inpl


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_c_array_4_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(dims(1), dims(2), dims(3), dims(4)), intent(inout) :: a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else
  
  continue

#endif

end subroutine mpireduce_sum_c_array_4_inpl


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_c_array_6_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)), &
     & intent(inout) :: a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPICOMPLEX_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else
  
  continue

#endif

end subroutine mpireduce_sum_c_array_6_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_r_array_2(a, sum_a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  real, dimension(dims(1),dims(2)), intent(in)  :: a
  real, dimension(dims(1),dims(2)), intent(out) :: sum_a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  
  nelem = product(dims)
  if (nelem > 0) then
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    call mpi_reduce(a, sum_a, nelem, MPIREAL_X, MPI_SUM, &
       & root_proc_rank, i_comm, ierr)
  end if
#else

  sum_a = a

#endif

end subroutine mpireduce_sum_r_array_2

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_c_array_2(a, sum_a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(dims(1),dims(2)), intent(in)  :: a
  complex, dimension(dims(1),dims(2)), intent(out) :: sum_a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  
  nelem = product(dims)
  if (nelem > 0) then
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    call mpi_reduce(a, sum_a, nelem, MPICOMPLEX_X, MPI_SUM, &
       & root_proc_rank, i_comm, ierr)
  end if
#else

  sum_a = a

#endif

end subroutine mpireduce_sum_c_array_2

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_r_array_2_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  real, dimension(dims(1),dims(2)), intent(inout)  :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  
  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else

  continue

#endif

end subroutine mpireduce_sum_r_array_2_inpl



!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_r_array_3(a, sum_a, n_elements, m_elements, &
   & k_elements, COMM)
  use global, only : logical_false

  integer, intent(in)                 :: n_elements, m_elements, k_elements
  real, dimension(:,:,:), intent(in)  :: a
  real, dimension(n_elements,m_elements,k_elements), &
          & intent(out)               :: sum_a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, i_elements 
  integer, parameter :: root_proc_rank = 0

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  
  i_elements = n_elements*m_elements*k_elements 
  if (i_elements > 0) then
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    call mpi_reduce(a, sum_a, i_elements, MPIREAL_X, MPI_SUM, &
       & root_proc_rank, i_comm, ierr)
  end if
  
#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,k_elements, &
      & COMM

#endif

end subroutine mpireduce_sum_r_array_3


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_r_array_3_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  real, dimension(dims(1),dims(2),dims(3)), intent(inout)  :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  
  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else

  continue

#endif

end subroutine mpireduce_sum_r_array_3_inpl


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_r_array_4_inpl(a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  real, dimension(dims(1),dims(2),dims(3),dims(4)), intent(inout)  :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, rank
  integer :: nelem
  integer, parameter :: root_proc_rank = 0

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  
  nelem = product(dims)
  if (nelem > 0) then

    call mpi_comm_rank(i_comm, rank, ierr)
    
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    if(rank == root_proc_rank) then
      call mpi_reduce(MPI_IN_PLACE, a, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    else
      call mpi_reduce(a, MPI_IN_PLACE, nelem, MPIREAL_X, MPI_SUM, &
         & root_proc_rank, i_comm, ierr)
    end if
  end if
#else

  continue

#endif

end subroutine mpireduce_sum_r_array_4_inpl



!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpireduce_sum_c_array_4(a, sum_a, dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(:,:,:,:), intent(in)  :: a
  complex, dimension(dims(1),dims(2),dims(3),dims(4)), &
          & intent(out)               :: sum_a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)
  integer, parameter :: root_proc_rank = 0
  integer :: nelem
  nelem = product(dims)
  if (nelem > 0) then
    ! sum the elements of all processes and store the resulting sum only in
    ! the root process.
    call mpi_reduce(a, sum_a, nelem, MPICOMPLEX_X, MPI_SUM, &
       & root_proc_rank, COMM, ierr)
  end if
  
#else
  sum_a = a
#endif
end subroutine mpireduce_sum_c_array_4


!-----------------------------------------------------------------------------
!> wrapper for mpi_allreduce MAXLOC
!> fortran mapping of reals to ranks is unfortunate, but that is how 
!> it is implemented in MPI.  An cleaner wrapper could be made however
!-----------------------------------------------------------------------------
subroutine mpiallreduce_maxloc(a,a_max,n,COMM)

  use global, only : logical_false
  ! the first element of a is the local scalar value, the second
  ! element of a is the local process rank.
  real, dimension(2), intent(in)  :: a
  ! the first element of a_max is the global maximum value,
  ! the second element is the rank of the process which holds it.
  real, dimension(2), intent(out) :: a_max
  integer, intent(in)             :: n
  integer, optional, intent(in)   :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call mpi_allreduce(a,a_max,1,MPI2REAL_X,MPI_MAXLOC,icomm,ierr)

#else

  a_max = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n,COMM

#endif

end subroutine mpiallreduce_maxloc


!-----------------------------------------------------------------------------
!> wrapper for logical scalar .or. MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_or_scalar(a,or_a,n_a,COMM)

  use global, only : logical_false

  logical, intent(in)           :: a
  integer, intent(in)           :: n_a
  logical, intent(out)          :: or_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n_a > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call MPI_ALLREDUCE(a,or_a,1,MPI_LOGICAL,MPI_LOR,icomm,ierr)

#else

  or_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_a,COMM

#endif

end subroutine mpiallreduce_or_scalar

!-----------------------------------------------------------------------------
!> wrapper for logical scalar .or. MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_or_scalar_inpl(a,COMM)
  use global, only : logical_false
  logical, intent(in)           :: a
  integer, optional, intent(in) :: COMM
#if defined(mpi2)
  integer :: icomm

  if (logical_false) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call mpi_allreduce(MPI_IN_PLACE,a,1,MPI_LOGICAL,MPI_LOR,icomm,ierr)

#else
  ! keep compiler quiet
  if (logical_false) write (*,*) a,COMM
#endif

end subroutine mpiallreduce_or_scalar_inpl


!-----------------------------------------------------------------------------
!> wrapper for logical scalar .and. MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_and_logical_scalar(a,and_a,n_a,COMM)

  use global, only : logical_false

  logical, intent(in)           :: a
  integer, intent(in)           :: n_a
  logical, intent(out)          :: and_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n_a > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call MPI_ALLREDUCE(a,and_a,1,MPI_LOGICAL,MPI_LAND,icomm,ierr)

#else

  and_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_a,COMM

#endif

end subroutine mpiallreduce_and_logical_scalar


!-----------------------------------------------------------------------------
!> wrapper for integer scalar max MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_max_int_scalar(a,max_a,n_a,COMM)

  use global, only : logical_false

  integer, intent(in)           :: a,n_a
  integer, intent(out)          :: max_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n_a > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call MPI_ALLREDUCE(a,max_a,1,MPI_INTEGER,MPI_MAX,icomm,ierr)

#else

  max_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_a,COMM

#endif

end subroutine mpiallreduce_max_int_scalar

!-----------------------------------------------------------------------------
!> wrapper for integer array max MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_max_int_array_1(a,max_a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)           :: n_elements
  integer, intent(in)           :: a(n_elements)
  integer, intent(out)          :: max_a(n_elements)
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n_elements > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call MPI_ALLREDUCE(a,max_a,n_elements,MPI_INTEGER,MPI_MAX,icomm,ierr)

#else

  max_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) a,max_a,n_elements,COMM

#endif

end subroutine mpiallreduce_max_int_array_1

!-----------------------------------------------------------------------------
!> wrapper for real scalar MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_max_real_scalar(a,max_a,n_a,COMM)

  use global, only : logical_false

  real, intent(in)              :: a
  integer, intent(in)           :: n_a
  real, intent(out)             :: max_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n_a > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call MPI_ALLREDUCE(a,max_a,1,MPIREAL_X,MPI_MAX,icomm,ierr)

#else

  max_a=a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_a,COMM

#endif

end subroutine mpiallreduce_max_real_scalar


!-----------------------------------------------------------------------------
!> wrapper for integer scalar MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_min_int_scalar(a,min_a,n_a,COMM)

  use global, only : logical_false

  integer, intent(in)           :: a, n_a
  integer, intent(out)          :: min_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n_a > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call MPI_ALLREDUCE(a,min_a,1,MPI_INTEGER,MPI_MIN,icomm,ierr)

#else

  min_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_a,COMM

#endif

end subroutine mpiallreduce_min_int_scalar


!-----------------------------------------------------------------------------
!> wrapper for real scalar MPI_ALLREDUCE
!-----------------------------------------------------------------------------
subroutine mpiallreduce_min_r_scalar(a,min_a,n_a,COMM)

  use global, only : logical_false

  real, intent(in)              :: a
  integer, intent(in)           :: n_a
  real, intent(out)             :: min_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: icomm

  if (logical_false) continue
  if (n_a > 0) continue

  icomm = MPI_COMM_WORLD
  if (present(COMM)) icomm = COMM
  call MPI_ALLREDUCE(a,min_a,1,MPIREAL_X,MPI_MIN,icomm,ierr)

#else

  min_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_a,COMM

#endif

end subroutine mpiallreduce_min_r_scalar

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_scalar(a,sum_a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)           :: n_elements
  complex, intent(in)           :: a
  complex, intent(out)          :: sum_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue
  if (n_elements > 0) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(a,sum_a,1,MPICOMPLEX_X,MPI_SUM,i_comm,ierr)

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM

#endif

end subroutine mpiallreduce_sum_c_scalar

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_int_scalar(a,sum_a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)           :: n_elements
  integer, intent(in)           :: a
  integer, intent(out)          :: sum_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue
  if (n_elements > 0) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(a,sum_a,1,MPI_INTEGER,MPI_SUM,i_comm,ierr)

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM

#endif

end subroutine mpiallreduce_sum_int_scalar

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_scalar(a,sum_a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)           :: n_elements
  real, intent(in)              :: a
  real, intent(out)             :: sum_a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue
  if (n_elements > 0) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(a,sum_a,1,MPIREAL_X,MPI_SUM,i_comm,ierr)

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM

#endif

end subroutine mpiallreduce_sum_r_scalar

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_i_inpl(a,n_elements,COMM)


  integer, intent(inout)  :: a
  !> the argument n_elements is ignored
  integer, intent(in) :: n_elements
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(MPI_IN_PLACE,a,1,MPI_INTEGER,MPI_SUM,i_comm,ierr)

#else

  continue

#endif

  if (n_elements > 0) continue

end subroutine mpiallreduce_sum_i_inpl


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_inpl(a,n_elements,COMM)

  use global, only : logical_false

  real, intent(inout)  :: a
  !> the argument n_elements is ignored
  integer, intent(in) :: n_elements
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue
  if (n_elements > 0) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(MPI_IN_PLACE,a,1,MPIREAL_X,MPI_SUM,i_comm,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM,a

#endif

end subroutine mpiallreduce_sum_r_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_inpl(a,n_elements,COMM)

  use global, only : logical_false

  complex, intent(inout)  :: a
  !> the argument n_elements is ignored
  integer, intent(in) :: n_elements
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue
  if (n_elements > 0) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(MPI_IN_PLACE,a,1,MPICOMPLEX_X,MPI_SUM,i_comm,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM,a

#endif

end subroutine mpiallreduce_sum_c_inpl

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_1(a,sum_a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                         :: n_elements
  complex, dimension(:), intent(in)           :: a
  complex, dimension(n_elements), intent(out) :: sum_a
  integer, optional, intent(in)               :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(a,sum_a,n_elements,MPICOMPLEX_X,MPI_SUM,i_comm,ierr)

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM

#endif

end subroutine mpiallreduce_sum_c_array_1

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_1_inpl(a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                   :: n_elements
  complex, dimension(:), intent(inout)  :: a
  integer, optional, intent(in)         :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  call MPI_ALLREDUCE(MPI_IN_PLACE,a,n_elements,MPICOMPLEX_X,MPI_SUM,i_comm,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM,a

#endif

end subroutine mpiallreduce_sum_c_array_1_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_2_inpl(a,dims,COMM)
  integer, dimension(2), intent(in) :: dims
  complex, dimension(dims(1), dims(2)), intent(inout) :: a
  integer, optional, intent(in) :: COMM
#if defined(mpi2)
  integer :: i_comm
  integer :: nelem
  
  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then
    call mpi_allreduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
       & i_comm, ierr)
  end if
#else

  continue
  
#endif

end subroutine mpiallreduce_sum_c_array_2_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_3_inpl(a, dims, COMM)
  integer, dimension(3), intent(in) :: dims
  complex, dimension(dims(1),dims(2),dims(3)), &
          & intent(inout)               :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)
  integer :: nelem
  integer :: i_comm

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then
    ! sum the elements of all processes
    call mpi_allreduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
       & i_comm, ierr)
  end if
#else

  continue

#endif
end subroutine mpiallreduce_sum_c_array_3_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_4_inpl(a, dims, COMM)
  integer, dimension(4), intent(in) :: dims
  complex, dimension(dims(1),dims(2),dims(3),dims(4)), &
          & intent(inout)               :: a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)
  integer :: nelem
  integer :: i_comm

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM

  nelem = product(dims)
  if (nelem > 0) then
    ! sum the elements of all processes
    call mpi_allreduce(MPI_IN_PLACE, a, nelem, MPICOMPLEX_X, MPI_SUM, &
       & i_comm, ierr)
  end if
#else

  continue

#endif
end subroutine mpiallreduce_sum_c_array_4_inpl


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_2(a,sum_a,n_elements,m_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                                    :: n_elements
  integer, intent(in)                                    :: m_elements
  complex, dimension(:,:), intent(in)                    :: a
  complex, dimension(n_elements,m_elements), intent(out) :: sum_a
  integer, optional, intent(in)                          :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  if (n_elements*m_elements > 0) then
    call MPI_ALLREDUCE(a,sum_a,n_elements*m_elements,MPICOMPLEX_X,MPI_SUM,   &
        &              i_comm,ierr)
  end if

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,COMM

#endif

end subroutine mpiallreduce_sum_c_array_2

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_3(a,sum_a,n_elements,m_elements, & 
                                    &  k_elements, COMM)

  use global, only : logical_false

  integer, intent(in)                    :: n_elements,m_elements,k_elements
  complex, dimension(:,:,:), intent(in)  :: a
  complex, dimension(n_elements,m_elements,k_elements), &
          & intent(out)                  :: sum_a
  integer, optional, intent(in)          :: COMM

#if defined(mpi2)

  integer :: i_comm, i_elements

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  i_elements = n_elements*m_elements*k_elements 
  if (i_elements > 0) then
    call MPI_ALLREDUCE(a,sum_a,i_elements,MPICOMPLEX_X,MPI_SUM,   &
        &              i_comm,ierr)
  end if

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,k_elements, &
      & COMM

#endif

end subroutine mpiallreduce_sum_c_array_3

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_4(a,sum_a,dims, COMM)
  integer, dimension(:), intent(in) :: dims
  complex, dimension(dims(1),dims(2),dims(3),dims(4)), intent(out) :: sum_a
  complex, dimension(dims(1),dims(2),dims(3),dims(4)), intent(in) :: a
  integer, optional, intent(in)            :: COMM

#if defined(mpi2)

  integer :: i_comm, i_elements 

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  i_elements = product(dims)
  if (i_elements > 0) then
    call mpi_allreduce(a,sum_a,i_elements,MPICOMPLEX_X,MPI_SUM, &
       & i_comm,ierr)
  end if

#else

  sum_a = a
  
#endif

end subroutine mpiallreduce_sum_c_array_4


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_c_array_5(a,sum_a,n_elements,m_elements, & 
                                  &  k_elements, p_elements, l_elements, COMM)

  use global, only : logical_false

  integer, intent(in)                        :: n_elements,m_elements,k_elements
  integer, intent(in)                        :: p_elements,l_elements
  complex, dimension(:,:,:,:,:), intent(in)  :: a
  complex, dimension(n_elements,m_elements,k_elements, p_elements,l_elements), &
          & intent(out)                      :: sum_a
  integer, optional, intent(in)              :: COMM

#if defined(mpi2)

  integer :: i_comm, i_elements

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  i_elements = n_elements*m_elements*k_elements*p_elements*l_elements 
  if (i_elements > 0) then
    call MPI_ALLREDUCE(a,sum_a,i_elements,MPICOMPLEX_X,MPI_SUM,   &
        &              i_comm,ierr)
  end if

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,k_elements, &
      & p_elements, l_elements, COMM

#endif

end subroutine mpiallreduce_sum_c_array_5

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_4(a,sum_a,n_elements,m_elements, & 
                                    &  k_elements, p_elements, COMM)

  use global, only : logical_false

  integer, intent(in)                   :: n_elements,m_elements
  integer, intent(in)                   :: k_elements,p_elements
  real, dimension(:,:,:,:), intent(in)  :: a
  real, dimension(n_elements,m_elements,k_elements, p_elements), &
          & intent(out)                 :: sum_a
  integer, optional, intent(in)         :: COMM

#if defined(mpi2)

  integer :: i_comm, i_elements

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  i_elements = n_elements*m_elements*k_elements*p_elements 
  if (i_elements > 0) then
    call MPI_ALLREDUCE(a,sum_a,i_elements,MPIREAL_X,MPI_SUM,   &
        &              i_comm,ierr)
  end if

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,k_elements, &
      & p_elements, COMM

#endif

end subroutine mpiallreduce_sum_r_array_4

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_3(a,sum_a,n_elements,m_elements, & 
                                    &  k_elements, COMM)

  use global, only : logical_false

  integer, intent(in)                 :: n_elements,m_elements,k_elements
  real, dimension(:,:,:), intent(in)  :: a
  real, dimension(n_elements,m_elements,k_elements), &
          & intent(out)               :: sum_a
  integer, optional, intent(in)       :: COMM

#if defined(mpi2)

  integer :: i_comm, i_elements 

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  i_elements = n_elements*m_elements*k_elements 
  if (i_elements > 0) then
    call MPI_ALLREDUCE(a,sum_a,i_elements,MPIREAL_X,MPI_SUM,   &
        &              i_comm,ierr)
  end if

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,k_elements, &
      & COMM

#endif

end subroutine mpiallreduce_sum_r_array_3

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_3_inpl(a,n_elements,m_elements, & 
                                    &  k_elements, COMM)

  use global, only : logical_false
    
  integer, intent(in)           :: n_elements,m_elements,k_elements
  real, dimension(n_elements,m_elements,k_elements), intent(inout)  :: a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)

  integer :: i_comm, i_elements

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  i_elements = n_elements*m_elements*k_elements 
  if (i_elements > 0) then
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,i_elements,MPIREAL_X,MPI_SUM,   &
        &              i_comm,ierr)
  end if

#else

  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,k_elements, &
      & COMM, a

#endif

end subroutine mpiallreduce_sum_r_array_3_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_4_inpl(a, n_elements, m_elements, &
   & k_elements, l_elements, COMM)
  use global, only : logical_false
  integer, intent(in) :: n_elements, m_elements, k_elements, l_elements
  real, dimension(n_elements, m_elements, k_elements, l_elements), &
     & intent(inout) :: a
  integer, optional, intent(in) :: COMM

#if defined(mpi2)
  integer :: i_comm, i_elements

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  i_elements = n_elements*m_elements*k_elements*l_elements
  if (i_elements > 0) then
    call mpi_allreduce(MPI_IN_PLACE, a, i_elements, MPIREAL_X, MPI_SUM, &
       & i_comm,ierr)
  end if

#else

  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,k_elements,l_elements, &
      & COMM, a

#endif

end subroutine mpiallreduce_sum_r_array_4_inpl


!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_2(a,sum_a,n_elements,m_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                                 :: n_elements,m_elements
  real, dimension(:,:), intent(in)                    :: a
  real, dimension(n_elements,m_elements), intent(out) :: sum_a
  integer, optional, intent(in)                       :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  if (n_elements*m_elements > 0) then
    call MPI_ALLREDUCE(a,sum_a,n_elements*m_elements,MPIREAL_X,MPI_SUM,      &
        &              i_comm,ierr)
  end if

#else

  sum_a = a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,COMM

#endif

end subroutine mpiallreduce_sum_r_array_2

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_2_inpl(a,n_elements,m_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                                  :: n_elements,m_elements
  real, dimension(n_elements,m_elements), intent(inout):: a
  integer, optional, intent(in)                        :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
  if (n_elements*m_elements > 0) then
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,n_elements*m_elements,MPIREAL_X,MPI_SUM,      &
        &              i_comm,ierr)
  end if

#else

  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,m_elements,COMM,a

#endif

end subroutine mpiallreduce_sum_r_array_2_inpl

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_1(a,sum_a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                      :: n_elements
  real, dimension(n_elements), intent(in)  :: a
  real, dimension(n_elements), intent(out) :: sum_a
  integer, optional, intent(in)            :: COMM

#if defined(mpi2)

  integer :: i_comm

  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm=COMM
  call MPI_ALLREDUCE(a,sum_a,n_elements,MPIREAL_X,MPI_SUM,i_comm,ierr)

#else

  sum_a=a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM

#endif

end subroutine mpiallreduce_sum_r_array_1

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------

#if defined(real_precision_default)
subroutine mpiallreduce_sum_d_array_1(a,sum_a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                                  :: n_elements
  double precision, dimension(n_elements), intent(in)  :: a
  double precision, dimension(n_elements), intent(out) :: sum_a
  integer, optional, intent(in)                        :: COMM

#if defined(mpi2)

  integer :: i_comm

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm=COMM
  call MPI_ALLREDUCE(a,sum_a,n_elements,MPI_DOUBLE_PRECISION,MPI_SUM,i_comm,ierr)

#else

  sum_a=a
  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM

#endif

end subroutine mpiallreduce_sum_d_array_1
#endif

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine mpiallreduce_sum_r_array_1_inpl(a,n_elements,COMM)

  use global, only : logical_false

  integer, intent(in)                         :: n_elements
  real, dimension(n_elements), intent(inout)  :: a
  integer, optional, intent(in)               :: COMM

#if defined(mpi2)

  integer :: i_comm

  ! keep compiler quiet about unused parameter.
  if (logical_false) continue

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm=COMM
  call MPI_ALLREDUCE(MPI_IN_PLACE,a,n_elements,MPIREAL_X,MPI_SUM,i_comm,ierr)

#else

  ! keep compiler quiet
  if (logical_false) write (*,*) n_elements,COMM,a

#endif

end subroutine mpiallreduce_sum_r_array_1_inpl

!-----------------------------------------------------------------------------
!> Get the displacement _idisp_ from the beginning of the file associated
!> with _file_unit_ of the file pointer on processor _iproc_, for the communicator
!> _COMM_. This is useful for setting the next file view if only 1 processor
!> has written data to a file open for parallel write.
!-----------------------------------------------------------------------------
subroutine mpigetfiledisp(file_unit,iproc,idisp,COMM)

  use global, only : logical_false

  integer, intent(in) :: iproc, file_unit
    
#if defined(mpi2)

  integer (KIND=MPI_OFFSET_KIND), intent(out) :: idisp
  integer (KIND=MPI_OFFSET_KIND)              :: offset
  integer, optional, intent(in)               :: COMM
  integer :: i_comm

  i_comm = MPI_COMM_WORLD
  if (present(COMM)) i_comm = COMM
    
  if (processor_number == iproc) then
    call MPI_FILE_GET_POSITION(file_unit,offset,ierr)
    call MPI_FILE_GET_BYTE_OFFSET(file_unit,offset,idisp,ierr)
  end if
  
  ! 8 integers should be long enough?
  call MPI_BCAST(idisp,8,MPI_INTEGER,iproc,i_comm,ierr)
    
#else

  integer, intent(in) :: idisp
  integer, optional, intent(in) :: COMM
  ! do nothing

#endif

  ! keep compiler quiet
  if (logical_false) write (*,*) file_unit,iproc,idisp,COMM

end subroutine mpigetfiledisp
  
!-----------------------------------------------------------------------------
!> Gather 1d integer array into a global array. This routine should only be
!> called by processors responsible for the different parts of the array.
!> The local array, local and global array sizes, communicators are input.
!> A global array is returned (on the processor with rank = 0 in the comm).
!----------------------------------------------------------------------------
subroutine gather_1d_int(global_x,n_x_grid,local_x,nx,comm_x,ALLGATHER)

  use global, only : logical_false

  integer, intent(in)                         :: n_x_grid,nx,comm_x
  integer, dimension(n_x_grid), intent(inout) :: global_x
  integer, dimension(nx), intent(in)          :: local_x
  integer, dimension(nx)                      :: tmp_x
  logical, optional, intent(in)               :: ALLGATHER
  
  ! use tmp_x for the MPI calls
  tmp_x(:) = local_x(:)
  
  ! gather in the X-direction
  if (n_x_grid > nx .and. mod(n_x_grid,nx) == 0) then
#if defined(mpi2)
  if (present(ALLGATHER)) then
    if (ALLGATHER) then
      call MPI_ALLGATHER(tmp_x,nx,MPI_INTEGER,global_x,nx,MPI_INTEGER,comm_x,ierr)
    else
      call MPI_GATHER(tmp_x,nx,MPI_INTEGER,global_x,nx,MPI_INTEGER,0,comm_x,ierr)
    end if
  else
    call MPI_GATHER(tmp_x,nx,MPI_INTEGER,global_x,nx,MPI_INTEGER,0,comm_x,ierr)
  end if
#endif
  else if (n_x_grid == nx) then
    global_x(:) = tmp_x(:)
  else
    call mpiabort('gather_1d_real: global and local array size mismatch')
  end if

  ! keep compiler quiet
  if (logical_false) write (*,*) comm_x,ALLGATHER

end subroutine gather_1d_int

!-----------------------------------------------------------------------------
!> Gather 1d real array into a global array. This routine should only be
!> called by processors responsible for the different parts of the array.
!> The local array, local and global array sizes, communicators are input.
!> A global array is returned (on the processor with rank = 0 in the comm).
!----------------------------------------------------------------------------
subroutine gather_1d_real(global_x,n_x_grid,local_x,nx,comm_x,ALLGATHER)

  use global, only : logical_false

  integer, intent(in)                      :: n_x_grid,nx,comm_x
  real, dimension(n_x_grid), intent(inout) :: global_x
  real, dimension(nx), intent(in)          :: local_x
  real, dimension(nx)                      :: tmp_x
  logical, optional, intent(in)            :: ALLGATHER
  
  ! use tmp_x for the MPI calls
  tmp_x(:) = local_x(:)
  
  ! gather in the X-direction
  if (n_x_grid > nx .and. mod(n_x_grid,nx) == 0) then
#if defined(mpi2)
  if (present(ALLGATHER)) then
    if (ALLGATHER) then
      call MPI_ALLGATHER(tmp_x,nx,MPIREAL_X,global_x,nx,MPIREAL_X,comm_x,ierr)
    else
      call MPI_GATHER(tmp_x,nx,MPIREAL_X,global_x,nx,MPIREAL_X,0,comm_x,ierr)
    end if
  else
    call MPI_GATHER(tmp_x,nx,MPIREAL_X,global_x,nx,MPIREAL_X,0,comm_x,ierr)
  end if
#endif
  else if (n_x_grid == nx) then
    global_x(:) = tmp_x(:)
  else
    call mpiabort('gather_1d_real: global and local array size mismatch')
  end if

  ! keep compiler quiet
  if (logical_false) write (*,*) comm_x,ALLGATHER

end subroutine gather_1d_real

!-----------------------------------------------------------------------------
!> This routine must be called by all processes of the communicator comm.
!>
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> IMPORTANT: If to_root_of_commcart is .true., then this routine must 
!> additionally be called by the root process (with respect to COMM_CART). The
!> comm and local_array arguments are ignored by root if root is not in the
!> communicator of the other processes.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator).
!----------------------------------------------------------------------------
subroutine gather_1d_real_intercomm(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
  use mpicomms, only : COMM_CART, COMM_ROOT
  use global, only : gkw_warn_any_proc

  integer, parameter :: data_ndims = 1
  real, dimension(:), &
     & intent(out) :: global_array
  real, dimension(:), &
     & intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart, root_is_in_comm
  integer, intent(in),optional :: tag

#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  real :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start

#if defined(mpi3)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif

  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  integer :: remote_leader_rank
  integer :: intercomm, tag_
  integer :: status_arr (MPI_STATUS_SIZE)
  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
    global_array = local_array
  else
    lshape = shape(local_array)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 0
    end if


    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)

    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 2d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if

    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(to_root_of_commcart .and. root_processor) then
      allocate(displacements(num_procs))
    else if(.not. to_root_of_commcart .and. cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if
    
    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPIREAL_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)

    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    !write (*,*) "bytes of a real number:", elem_bytes
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)

    if(.not. to_root_of_commcart .or. root_is_in_comm) then
      ! Collectively gather using an intracommunicator
      
      ! every process sends its displacement to root
      call gather_array(displacements, num_procs, displacements, 1, comm)

      counts = 1
      call mpi_gatherv(local_array, 1, local_subarray_type, &
         & global_array, counts, &
         & displacements, resized_subarray_type, 0, comm, ierr)
    else
      ! That means that the gather is to root of commcart and that
      ! root is not in comm.
      
      ! Collectively gather using an intercommunicator
      
      if(root_processor) then
        
        ! root needs to know the rank (with respect to COMM_CART) of
        ! the remote group's leader. Send this via point-to-point
        ! communication.
        call mpi_recv(remote_leader_rank, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, COMM_CART, status_arr, ierr);
        call mpi_intercomm_create(COMM_ROOT, 0, COMM_CART, &
           & remote_leader_rank, tag_, intercomm, ierr)
      else
        if (cart_rank == 0) then
          ! if this is root (i.e. here used as the leader), then
          ! send its rank (with respect to COMM_CART) to root (with
          ! respect to COMM_CART)
          call mpi_send(processor_number, 1, MPI_INTEGER, 0, tag_, COMM_CART, ierr);
        end if
        call mpi_intercomm_create(comm, 0, COMM_CART, &
           & 0, tag_, intercomm, ierr)
      end if

      ! now use the obtained intercommunicator for a gather
      if(root_processor) then
        ! if there were other processes in this group, they would
        ! have to specify MPI_PROC_NULL instead of MPI_ROOT

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, MPI_ROOT, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, MPI_ROOT, intercomm, ierr)
      else
        ! the receiver rank parameter here is the root in the remote group

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, 0, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, 0, intercomm, ierr)
      end if

      call mpi_comm_free(intercomm, ierr);
    end if

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)
    
  end if
#else
  global_array = local_array
#endif
  
end subroutine gather_1d_real_intercomm

!-----------------------------------------------------------------------------
!> This routine must be called by all processes of the communicator comm.
!>
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> IMPORTANT: If to_root_of_commcart is .true., then this routine must
!> additionally be called by the root process (with respect to COMM_CART). The
!> comm and local_array arguments are ignored by root if root is not in the
!> communicator of the other processes.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator).
!----------------------------------------------------------------------------
subroutine gather_1d_cmplx_intercomm(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
  use mpicomms, only : COMM_CART, COMM_ROOT
  use global, only : gkw_warn_any_proc

  integer, parameter :: data_ndims = 1
  complex, dimension(:), intent(out) :: global_array
  complex, dimension(:), intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart, root_is_in_comm
  integer, intent(in),optional :: tag

#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  complex :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start

#if defined(mpi3)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif

  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  integer :: remote_leader_rank
  integer :: intercomm, tag_
  integer :: status_arr (MPI_STATUS_SIZE)
  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
    global_array = local_array
  else
    lshape = shape(local_array)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 0
    end if


    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)

    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 2d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if

    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(to_root_of_commcart .and. root_processor) then
      allocate(displacements(num_procs))
    else if(.not. to_root_of_commcart .and. cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if

    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPIREAL_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)

    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    !write (*,*) "bytes of a real number:", elem_bytes
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)

    if(.not. to_root_of_commcart .or. root_is_in_comm) then
      ! Collectively gather using an intracommunicator

      ! every process sends its displacement to root
      call gather_array(displacements, num_procs, displacements, 1, comm)

      counts = 1
      call mpi_gatherv(local_array, 1, local_subarray_type, &
         & global_array, counts, &
         & displacements, resized_subarray_type, 0, comm, ierr)
    else
      ! That means that the gather is to root of commcart and that
      ! root is not in comm.

      ! Collectively gather using an intercommunicator

      if(root_processor) then

        ! root needs to know the rank (with respect to COMM_CART) of
        ! the remote group's leader. Send this via point-to-point
        ! communication.
        call mpi_recv(remote_leader_rank, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, COMM_CART, status_arr, ierr);
        call mpi_intercomm_create(COMM_ROOT, 0, COMM_CART, &
           & remote_leader_rank, tag_, intercomm, ierr)
      else
        if (cart_rank == 0) then
          ! if this is root (i.e. here used as the leader), then
          ! send its rank (with respect to COMM_CART) to root (with
          ! respect to COMM_CART)
          call mpi_send(processor_number, 1, MPI_INTEGER, 0, tag_, COMM_CART, ierr);
        end if
        call mpi_intercomm_create(comm, 0, COMM_CART, &
           & 0, tag_, intercomm, ierr)
      end if

      ! now use the obtained intercommunicator for a gather
      if(root_processor) then
        ! if there were other processes in this group, they would
        ! have to specify MPI_PROC_NULL instead of MPI_ROOT

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, MPI_ROOT, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, MPI_ROOT, intercomm, ierr)
      else
        ! the receiver rank parameter here is the root in the remote group

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, 0, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, 0, intercomm, ierr)
      end if

      call mpi_comm_free(intercomm, ierr);
    end if

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)

  end if
#else
  global_array = local_array
#endif

end subroutine gather_1d_cmplx_intercomm


!-----------------------------------------------------------------------------
!> Gather 1d complex array into a global array. This routine should only be
!> called by processors responsible for the different parts of the array.
!> The local array, local and global array sizes, communicators are input.
!> A global array is returned (on the processor with rank = 0 in the comm).
!----------------------------------------------------------------------------
subroutine gather_1d_cmplx(global_x,n_x_grid,local_x,nx,comm_x,ALLGATHER)

  use global, only : logical_false

  integer, intent(in)                         :: n_x_grid,nx,comm_x
  complex, dimension(n_x_grid), intent(inout) :: global_x
  complex, dimension(nx), intent(in)          :: local_x
  complex, dimension(nx)                      :: tmp_x
  logical, optional, intent(in)               :: ALLGATHER
  
  ! use tmp_x for the MPI calls
  tmp_x(:) = local_x(:)
  
  ! gather in the X-direction
  if (n_x_grid > nx .and. mod(n_x_grid,nx) == 0) then
#if defined(mpi2)
  if (present(ALLGATHER)) then
    if (ALLGATHER) then
      call MPI_ALLGATHER(tmp_x,nx,MPICOMPLEX_X,global_x,nx,MPICOMPLEX_X,comm_x,ierr)
    else
      call MPI_GATHER(tmp_x,nx,MPICOMPLEX_X,global_x,nx,MPICOMPLEX_X,0,comm_x,ierr)
    end if
  else
    call MPI_GATHER(tmp_x,nx,MPICOMPLEX_X,global_x,nx,MPICOMPLEX_X,0,comm_x,ierr)
  end if
#endif
  else if (n_x_grid == nx) then
    global_x(:) = tmp_x(:)
  else
    call mpiabort('gather_1d_cmplx: global and local array size mismatch')
  end if

  ! keep compiler quiet
  if (logical_false) write (*,*) comm_x,ALLGATHER

end subroutine gather_1d_cmplx

!-----------------------------------------------------------------------------
!> Gather 2d real slices in a global slice. This routine should (usually)
!> only be called by processors responsible for the slice. The local array,
!> local and global array sizes, communicators are input. A
!> global array is returned (on the processor with rank = 0 in both comms).
!----------------------------------------------------------------------------
subroutine gather_2d_real(global_x_y,n_x_grid,n_y_grid,local_x_y,nx,ny,      &
    &                           comm_x,comm_y,ALLGATHER)

  use global, only : logical_false

  integer, intent(in)                             :: n_x_grid,n_y_grid,nx,ny
  integer, intent(in)                             :: comm_x,comm_y
  real, dimension(n_x_grid,n_y_grid), intent(out) :: global_x_y
  real, dimension(nx,ny), intent(in)              :: local_x_y
  real, dimension(nx,ny)                          :: tmp_x_y
  logical, optional, intent(in)                   :: ALLGATHER
  real, dimension(nx,n_y_grid)                    :: part_x_y  
#if defined(mpi2)
  real, dimension(n_x_grid*n_y_grid)              :: global_x_y_tmp
  integer                                         :: ind,i,j,k
#endif

  ! use tmp_x_y as the array for MPI calls
  tmp_x_y(:,:) = local_x_y(:,:)
  
  ! first gather in the Y-direction
  if (n_y_grid > ny .and. mod(n_y_grid,ny) == 0) then
#if defined(mpi2)
    gather_y : if (present(ALLGATHER)) then
      if (ALLGATHER) then
        call MPI_ALLGATHER(tmp_x_y,nx*ny,MPIREAL_X,part_x_y,nx*ny,MPIREAL_X,   &
            &           comm_y,ierr)
      else
        call MPI_GATHER(tmp_x_y,nx*ny,MPIREAL_X,part_x_y,nx*ny,MPIREAL_X,0,    &
            &           comm_y,ierr)
      end if
    else
      call MPI_GATHER(tmp_x_y,nx*ny,MPIREAL_X,part_x_y,nx*ny,MPIREAL_X,0,      &
          &           comm_y,ierr)
    end if gather_y
#endif
  else if (n_y_grid == ny) then
    part_x_y(:,:) = tmp_x_y(:,:)
  else
    call mpiabort('gather_2d_real: dimension 2 grid size mismatch')
  end if

  ! gather in the X-direction
  if (n_x_grid > nx .and. mod(n_x_grid,nx) == 0) then
#if defined(mpi2)
    gather_x : if (present(ALLGATHER)) then
      if (ALLGATHER) then
        call MPI_ALLGATHER(part_x_y,nx*n_y_grid,MPIREAL_X,global_x_y_tmp,      &
            & nx*n_y_grid,MPIREAL_X,comm_x,ierr)
      else  
        call MPI_GATHER(part_x_y,nx*n_y_grid,MPIREAL_X,global_x_y_tmp,         &
            & nx*n_y_grid,MPIREAL_X,0,comm_x,ierr)
      end if
    else
      call MPI_GATHER(part_x_y,nx*n_y_grid,MPIREAL_X,global_x_y_tmp,           &
          & nx*n_y_grid,MPIREAL_X,0,comm_x,ierr)
    end if gather_x
    ! Need to re-order the gathered data from the second gather; the
    ! first gather is in the correct order.
    ind=0
    re_order : do k=1,(n_x_grid/nx)
      do j=1,n_y_grid
        do i=1+(k-1)*nx,k*nx
          ind=ind+1
          global_x_y(i,j) = global_x_y_tmp(ind)
        end do
      end do
    end do re_order
#endif
  else if (n_x_grid == nx) then
    global_x_y(:,:) = part_x_y(:,:)
  else
    call mpiabort('gather_2d_real: dimension 1 gridsize mismatch')
  end if

  ! keep compiler quiet
  if (logical_false) write (*,*) comm_x,comm_y,ALLGATHER

end subroutine gather_2d_real

!-----------------------------------------------------------------------------
!> Gather 2d complex slices in a global slice. This routine should (usually)
!> only be called by processors responsible for the slice. The local array,
!> local and global array sizes, communicators are input. A
!> global array is returned (on the processor with rank = 0 in both comms).
!----------------------------------------------------------------------------
subroutine gather_2d_cmplx(global_x_y,n_x_grid,n_y_grid,local_x_y,nx,ny,      &
    &                           comm_x,comm_y,ALLGATHER)

  use global, only : logical_false
        
  integer, intent(in)                                :: n_x_grid,n_y_grid,nx,ny
  integer, intent(in)                                :: comm_x,comm_y
  complex, dimension(n_x_grid,n_y_grid), intent(out) :: global_x_y
  complex, dimension(nx,ny), intent(in)              :: local_x_y
  complex, dimension(nx,ny)                          :: tmp_x_y
  logical, optional, intent(in)                      :: ALLGATHER
  complex, dimension(nx,n_y_grid)                    :: part_x_y  
#if defined(mpi2)
  complex, dimension(n_x_grid*n_y_grid)              :: global_x_y_tmp
  integer                                            :: ind,i,j,k
#endif 

  ! use tmp_x_y as the array for MPI calls
  tmp_x_y(:,:) = local_x_y(:,:)
  
  ! first gather in the Y-direction
  if (n_y_grid > ny .and. mod(n_y_grid,ny) == 0) then
#if defined(mpi2)
    gather_y : if (present(ALLGATHER)) then
      if (ALLGATHER) then
        call MPI_ALLGATHER(tmp_x_y,nx*ny,MPICOMPLEX_X,part_x_y,nx*ny,MPICOMPLEX_X,   &
            &           comm_y,ierr)
      else
        call MPI_GATHER(tmp_x_y,nx*ny,MPICOMPLEX_X,part_x_y,nx*ny,MPICOMPLEX_X,0,    &
            &           comm_y,ierr)
      end if
    else
      call MPI_GATHER(tmp_x_y,nx*ny,MPICOMPLEX_X,part_x_y,nx*ny,MPICOMPLEX_X,0,      &
          &           comm_y,ierr)
    end if gather_y
#endif
  else if (n_y_grid == ny) then
    part_x_y(:,:) = tmp_x_y(:,:)
  else
    call mpiabort('gather_2d_cmplx: dimension 2 grid size mismatch')
  end if

  ! gather in the X-direction
  if (n_x_grid > nx .and. mod(n_x_grid,nx) == 0) then
#if defined(mpi2)
    gather_x : if (present(ALLGATHER)) then
      if (ALLGATHER) then
        call MPI_ALLGATHER(part_x_y,nx*n_y_grid,MPICOMPLEX_X,global_x_y_tmp,      &
            & nx*n_y_grid,MPICOMPLEX_X,comm_x,ierr)
      else
        call MPI_GATHER(part_x_y,nx*n_y_grid,MPICOMPLEX_X,global_x_y_tmp,         &
            & nx*n_y_grid,MPICOMPLEX_X,0,comm_x,ierr)
      end if
    else
      call MPI_GATHER(part_x_y,nx*n_y_grid,MPICOMPLEX_X,global_x_y_tmp,           &
          & nx*n_y_grid,MPICOMPLEX_X,0,comm_x,ierr)
    end if gather_x
    ! Need to re-order the gathered data from the second gather; the
    ! first gather is in the correct order.
    ind=0
    re_order : do k=1,(n_x_grid/nx)
      do j=1,n_y_grid
        do i=1+(k-1)*nx,k*nx
          ind=ind+1
          global_x_y(i,j) = global_x_y_tmp(ind)
        end do
      end do
    end do re_order
#endif
  else if (n_x_grid == nx) then
    global_x_y(:,:) = part_x_y(:,:)
  else
    call mpiabort('gather_2d_cmplx: dimension 1 gridsize mismatch')
  end if

  ! keep compiler quiet
  if (logical_false) write (*,*) comm_x,comm_y,ALLGATHER

end subroutine gather_2d_cmplx

!-----------------------------------------------------------------------------
!> Gather 2d integer slices in a global slice. This routine should (usually)
!> only be called by processors responsible for the slice. The local array,
!> local and global array sizes, communicators are input. A
!> global array is returned (on the processor with rank = 0 in both comms).
!----------------------------------------------------------------------------
subroutine gather_2d_int(global_x_y,n_x_grid,n_y_grid,local_x_y,nx,ny,      &
    &                           comm_x,comm_y,ALLGATHER)

  integer, intent(in)                                :: n_x_grid,n_y_grid,nx,ny
  integer, intent(in)                                :: comm_x,comm_y
  integer, dimension(n_x_grid,n_y_grid), intent(out) :: global_x_y
  integer, dimension(nx,ny), intent(in)              :: local_x_y
  logical, optional, intent(in)                      :: ALLGATHER
  
  integer, dimension(nx,ny)                          :: tmp_x_y
  integer, dimension(nx,n_y_grid)                    :: part_x_y  
#if defined(mpi2)
  integer, dimension(n_x_grid*n_y_grid)              :: global_x_y_tmp
  integer                                            :: ind,i,j,k
#endif 

  ! use tmp_x_y as the array for MPI calls
  tmp_x_y(:,:) = local_x_y(:,:)
  
  ! first gather in the Y-direction
  if (n_y_grid > ny .and. mod(n_y_grid,ny) == 0) then
#if defined(mpi2)
    gather_y : if (present(ALLGATHER)) then
      if (ALLGATHER) then
        call MPI_ALLGATHER(tmp_x_y,nx*ny,MPI_INTEGER,part_x_y,nx*ny,MPI_INTEGER,   &
            &           comm_y,ierr)
      else
        call MPI_GATHER(tmp_x_y,nx*ny,MPI_INTEGER,part_x_y,nx*ny,MPI_INTEGER,0,    &
            &           comm_y,ierr)
      end if
    else
      call MPI_GATHER(tmp_x_y,nx*ny,MPI_INTEGER,part_x_y,nx*ny,MPI_INTEGER,0,      &
          &           comm_y,ierr)
    end if gather_y
#endif
  else if (n_y_grid == ny) then
    part_x_y(:,:) = tmp_x_y(:,:)
  else
    call mpiabort('gather_2d_int: dimension 2 grid size mismatch')
  end if

  ! gather in the X-direction
  if (n_x_grid > nx .and. mod(n_x_grid,nx) == 0) then
#if defined(mpi2)
    gather_x : if (present(ALLGATHER)) then
      if (ALLGATHER) then
        call MPI_ALLGATHER(part_x_y,nx*n_y_grid,MPI_INTEGER,global_x_y_tmp,      &
            & nx*n_y_grid,MPI_INTEGER,comm_x,ierr)
      else
        call MPI_GATHER(part_x_y,nx*n_y_grid,MPI_INTEGER,global_x_y_tmp,         &
            & nx*n_y_grid,MPI_INTEGER,0,comm_x,ierr)
      end if
    else
      call MPI_GATHER(part_x_y,nx*n_y_grid,MPI_INTEGER,global_x_y_tmp,           &
          & nx*n_y_grid,MPI_INTEGER,0,comm_x,ierr)
    end if gather_x
    ! Need to re-order the gathered data from the second gather; the
    ! first gather is in the correct order.
    ind=0
    re_order : do k=1,(n_x_grid/nx)
      do j=1,n_y_grid
        do i=1+(k-1)*nx,k*nx
          ind=ind+1
          global_x_y(i,j) = global_x_y_tmp(ind)
        end do
      end do
    end do re_order
#endif
  else if (n_x_grid == nx) then
    global_x_y(:,:) = part_x_y(:,:)
  else
    call mpiabort('gather_2d_int: dimension 1 gridsize mismatch')
  end if

end subroutine gather_2d_int


function proc_is_with_root_in(comm)
  integer, intent(in) :: comm
  logical :: proc_is_with_root_in
  integer :: i

  if(root_processor) then
    i = 1
  else
    i = 0
  end if
  call mpiallreduce_sum_inplace(i, 1,comm)
  proc_is_with_root_in = (i == 1)

end function proc_is_with_root_in


!-----------------------------------------------------------------------------
!> This routine must be called by all processes of the communicator comm.
!>
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> IMPORTANT: If to_root_of_commcart is .true., then this routine must 
!> additionally be called by the root process (with respect to COMM_CART). The
!> comm and local_array arguments are ignored by root if root is not in the
!> communicator of the other processes.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator).
!----------------------------------------------------------------------------
subroutine gather_2d_real_intercomm(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
  use mpicomms, only : COMM_CART, COMM_ROOT
  use global, only : gkw_warn_any_proc

  integer, parameter :: data_ndims = 2
  real, dimension(:,:), &
     & intent(out) :: global_array
  real, dimension(:,:), &
     & intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart, root_is_in_comm
  integer, intent(in),optional :: tag

#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  real :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start

#if defined(mpi3)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif

  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  integer :: remote_leader_rank
  integer :: intercomm, tag_
  integer :: status_arr (MPI_STATUS_SIZE)
  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
    global_array = local_array
  else
    lshape = shape(local_array)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 0
    end if

    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)

    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 2d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if

    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(to_root_of_commcart .and. root_processor) then
      allocate(displacements(num_procs))
    else if(.not. to_root_of_commcart .and. cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if
    
    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPIREAL_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)

    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    !write (*,*) "bytes of a real number:", elem_bytes
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)

    if(.not. to_root_of_commcart .or. root_is_in_comm) then
      ! Collectively gather using an intracommunicator
      
      ! every process sends its displacement to root
      call gather_array(displacements, num_procs, displacements, 1, comm)

      counts = 1
      call mpi_gatherv(local_array, 1, local_subarray_type, &
         & global_array, counts, &
         & displacements, resized_subarray_type, 0, comm, ierr)
    else
      ! That means that the gather is to root of commcart and that
      ! root is not in comm.
      
      ! Collectively gather using an intercommunicator
      
      if(root_processor) then
        
        ! root needs to know the rank (with respect to COMM_CART) of
        ! the remote group's leader. Send this via point-to-point
        ! communication.
        call mpi_recv(remote_leader_rank, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, COMM_CART, status_arr, ierr);
        call mpi_intercomm_create(COMM_ROOT, 0, COMM_CART, &
           & remote_leader_rank, tag_, intercomm, ierr)
      else
        if (cart_rank == 0) then
          ! if this is root (i.e. here used as the leader), then
          ! send its rank (with respect to COMM_CART) to root (with
          ! respect to COMM_CART)
          call mpi_send(processor_number, 1, MPI_INTEGER, 0, tag_, COMM_CART, ierr);
        end if
        call mpi_intercomm_create(comm, 0, COMM_CART, &
           & 0, tag_, intercomm, ierr)
      end if

      ! now use the obtained intercommunicator for a gather
      if(root_processor) then
        ! if there were other processes in this group, they would
        ! have to specify MPI_PROC_NULL instead of MPI_ROOT

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, MPI_ROOT, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, MPI_ROOT, intercomm, ierr)
      else
        ! the receiver rank parameter here is the root in the remote group

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, 0, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, 0, intercomm, ierr)
      end if

      call mpi_comm_free(intercomm, ierr);
    end if

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)
    
  end if
#else
  global_array = local_array
#endif
  
end subroutine gather_2d_real_intercomm

!-----------------------------------------------------------------------------
!> This routine must be called by all processes of the communicator comm.
!>
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> IMPORTANT: If to_root_of_commcart is .true., then this routine must 
!> additionally be called by the root process (with respect to COMM_CART). The
!> comm and local_array arguments are ignored by root if root is not in the
!> communicator of the other processes.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator if to_root_of_commcart=F, or root of COMM_CART
!> if it is T).
!----------------------------------------------------------------------------
subroutine gather_3d_real(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
  use mpicomms, only : COMM_CART, COMM_ROOT
  use global, only : gkw_warn_any_proc

  integer, parameter :: data_ndims = 3
  real, dimension(:,:,:), &
     & intent(out) :: global_array
  real, dimension(:,:,:), &
     & intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart, root_is_in_comm
  integer, intent(in),optional :: tag
  
#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  real :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start

#if defined(mpi3)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif

  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  integer :: remote_leader_rank
  integer :: intercomm, tag_
  integer :: status_arr (MPI_STATUS_SIZE)
  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
    global_array = local_array
  else
    lshape = shape(local_array)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 0
    end if


    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)

    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 3d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if

    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(to_root_of_commcart .and. root_processor) then
      allocate(displacements(num_procs))
    else if(.not. to_root_of_commcart .and. cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if
    
    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPIREAL_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)


    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)

    if(.not. to_root_of_commcart .or. root_is_in_comm) then
      ! Collectively gather using an intracommunicator
      
      ! every process sends its displacement to root
      call gather_array(displacements, num_procs, displacements, 1, comm)

      counts = 1
      ! write (*,*) tag_, cart_rank, "lshape:", lshape
      ! write (*,*) tag_, cart_rank, "starts:", starts
      ! write (*,*) tag_, cart_rank, "bytes of a real number:", elem_bytes
      ! write (*,*) tag_, cart_rank, "displacements:", displacements
      ! write (*,*) tag_, cart_rank, "num_procs:", num_procs
      ! write (*,*) tag_, cart_rank, "counts:", counts
      ! write (*,*) tag_, cart_rank, "global_array shape:", shape(global_array)
      call mpi_gatherv(local_array, 1, local_subarray_type, &
         & global_array, counts, &
         & displacements, resized_subarray_type, 0, comm, ierr)
    else
      ! That means that the gather is to root of commcart and that
      ! root is not in comm.
      
      ! Collectively gather using an intercommunicator
      
      ! Why this: The problem is that with serial output, only the root
      ! process can write files.
      ! However, the root process does not necessarily work on the
      ! hyperslab which shall be output, for example, one wants to
      ! output the 3d density for species 2, but there is
      ! parallelisation on the species and root works on species
      ! 1. Hence, to output serially, one has to gather over the
      ! communicator COMM_S_NE_X_NE (of processes working on the slab), but
      ! not to the root of that communicator, but to the root of
      ! COMM_CART!
      
      ! To achieve this, MPI also supports inter-communication: communication
      ! between two non-overlapping groups of processes.
      
      ! The group containing a process that initiates an
      ! inter-communication operation is called the "local group"
      ! that is, the sender in a send and the receiver in a receive.

      ! As in intra-communication, the target process is specified
      ! using a (communicator, rank) pair. Unlike
      ! intra-communication, the rank is relative to a second,
      ! remote group (called peer group in the docs).

      if(root_processor) then
        
        ! root needs to know the rank (with respect to COMM_CART) of
        ! the remote group's leader. Send this via point-to-point
        ! communication.
        call mpi_recv(remote_leader_rank, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, COMM_CART, status_arr, ierr);
        call mpi_intercomm_create(COMM_ROOT, 0, COMM_CART, &
           & remote_leader_rank, tag_, intercomm, ierr)
      else
        if (cart_rank == 0) then
          ! if this is root (i.e. here used as the leader), then
          ! send its rank (with respect to COMM_CART) to root (with
          ! respect to COMM_CART)
          call mpi_send(processor_number, 1, MPI_INTEGER, 0, tag_, COMM_CART, ierr);
        end if
        call mpi_intercomm_create(comm, 0, COMM_CART, &
           & 0, tag_, intercomm, ierr)
      end if

      ! now use the obtained intercommunicator for a gather
      if(root_processor) then
        ! if there were other processes in this group, they would
        ! have to specify MPI_PROC_NULL instead of MPI_ROOT

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, MPI_ROOT, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, MPI_ROOT, intercomm, ierr)
      else
        ! the receiver rank parameter here is the root in the remote group

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, 0, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, 0, intercomm, ierr)
      end if

      call mpi_comm_free(intercomm, ierr);
    end if

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)
    
  end if
#else
  global_array = local_array
#endif
  
end subroutine gather_3d_real

!-----------------------------------------------------------------------------
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator).
!----------------------------------------------------------------------------
subroutine gather_3d_cmplx(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
  use mpicomms, only : COMM_CART, COMM_ROOT
  use global, only : gkw_warn_any_proc

  integer, parameter :: data_ndims = 3
  complex, dimension(:,:,:), &
     & intent(out) :: global_array
  complex, dimension(:,:,:), &
     & intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart
  logical, intent(in) :: root_is_in_comm
  integer, intent(in),optional :: tag

#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  complex :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start
#if defined(mpi3)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif
  
  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  integer :: remote_leader_rank
  integer :: intercomm, tag_
  integer :: status_arr (MPI_STATUS_SIZE)

  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
    global_array = local_array
  else
    lshape = shape(local_array)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 0
    end if

    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)
    
    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 3d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if
        
    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(to_root_of_commcart .and. root_processor) then
      allocate(displacements(num_procs))
    else if(.not. to_root_of_commcart .and. cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if

    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPICOMPLEX_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)

    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    !call mpi_type_get_true_extent(resized_subarray_type, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)

    if(.not. to_root_of_commcart .or. root_is_in_comm) then
      ! Collectively gather using an intracommunicator
      
      ! gather all displacements on root
      call gather_array(displacements, num_procs, displacements, 1, comm)

      ! gather the data, using the displacements for the data coming
      ! from the respective process
      counts = 1
      ! write (*,*) tag_, cart_rank, "lshape=", lshape
      ! write (*,*) tag_, cart_rank, "starts=", starts
      ! write (*,*) tag_, cart_rank, "bytes of a complex number=", elem_bytes
      ! write (*,*) tag_, cart_rank, "displacements=", displacements
      ! write (*,*) tag_, cart_rank, "num_procs=", num_procs
      ! write (*,*) tag_, cart_rank, "counts=", counts
      ! write (*,*) tag_, cart_rank, "global_array shape=", shape(global_array)
      call mpi_gatherv(local_array, 1, local_subarray_type, &
         & global_array, counts, &
         & displacements, resized_subarray_type, 0, comm, ierr)
    else
      ! That means that the gather is to root of commcart and that
      ! root is not in comm.
      
      ! Collectively gather using an intercommunicator

      if(root_processor) then
        
        ! root needs to know the rank (with respect to COMM_CART) of
        ! the remote group's leader. Send this via point-to-point
        ! communication.
        call mpi_recv(remote_leader_rank, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, COMM_CART, status_arr, ierr);
        call mpi_intercomm_create(COMM_ROOT, 0, COMM_CART, &
           & remote_leader_rank, tag_, intercomm, ierr)
      else
        if (cart_rank == 0) then
          ! if this is root (i.e. here used as the leader), then
          ! send its rank (with respect to COMM_CART) to root (with
          ! respect to COMM_CART)
          call mpi_send(processor_number, 1, MPI_INTEGER, 0, tag_, COMM_CART, ierr);
        end if
        call mpi_intercomm_create(comm, 0, COMM_CART, &
           & 0, tag_, intercomm, ierr)
      end if

      ! now use the obtained intercommunicator for a gather
      if(root_processor) then
        ! if there were other processes in this group, they would
        ! have to specify MPI_PROC_NULL instead of MPI_ROOT

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, MPI_ROOT, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, MPI_ROOT, intercomm, ierr)
      else
        ! the receiver rank parameter here is the root in the remote group

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, 0, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, 0, intercomm, ierr)
      end if

      call mpi_comm_free(intercomm, ierr);
    end if

    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
  end if
#else
  global_array = local_array
#endif
  
end subroutine gather_3d_cmplx

!-----------------------------------------------------------------------------
!> This routine must be called by all processes of the communicator comm.
!>
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> IMPORTANT: If to_root_of_commcart is .true., then this routine must 
!> additionally be called by the root process (with respect to COMM_CART). The
!> comm and local_array arguments are ignored by root if root is not in the
!> communicator of the other processes.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator).
!----------------------------------------------------------------------------
subroutine gather_4d_real(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
  use mpicomms, only : COMM_CART, COMM_ROOT
  use global, only : gkw_warn_any_proc

  integer, parameter :: data_ndims = 4
  real, dimension(:,:,:,:), &
     & intent(out) :: global_array
  real, dimension(:,:,:,:), &
     & intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart, root_is_in_comm
  integer, intent(in),optional :: tag

#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  real :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start

#if defined(mpi3)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif

  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  integer :: remote_leader_rank
  integer :: intercomm, tag_
  integer :: status_arr (MPI_STATUS_SIZE)
  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
    global_array = local_array
  else
    lshape = shape(local_array)

    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)

    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 3d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if

    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(to_root_of_commcart .and. root_processor) then
      allocate(displacements(num_procs))
    else if(.not. to_root_of_commcart .and. cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if
    
    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPIREAL_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)

    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    !write (*,*) "bytes of a real number:", elem_bytes
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)

    if(.not. to_root_of_commcart .or. root_is_in_comm) then
      ! Collectively gather using an intracommunicator
      
      ! every process sends its displacement to root
      call gather_array(displacements, num_procs, displacements, 1, comm)

      counts = 1
      call mpi_gatherv(local_array, 1, local_subarray_type, &
         & global_array, counts, &
         & displacements, resized_subarray_type, 0, comm, ierr)
    else
      ! That means that the gather is to root of commcart and that
      ! root is not in comm.
      
      ! Collectively gather using an intercommunicator
      
      ! Why this: The problem is that with serial output, only the root
      ! process can write files.
      ! However, the root process does not necessarily work on the
      ! hyperslab which shall be output, for example, one wants to
      ! output the 3d density for species 2, but there is
      ! parallelisation on the species and root works on species
      ! 1. Hence, to output serially, one has to gather over the
      ! communicator COMM_S_NE_X_NE (of processes working on the slab), but
      ! not to the root of that communicator, but to the root of
      ! COMM_CART!
      
      ! To achieve this, MPI also supports inter-communication: communication
      ! between two non-overlapping groups of processes.
      
      ! The group containing a process that initiates an
      ! inter-communication operation is called the "local group"
      ! that is, the sender in a send and the receiver in a receive.

      ! As in intra-communication, the target process is specified
      ! using a (communicator, rank) pair. Unlike
      ! intra-communication, the rank is relative to a second,
      ! remote group (called peer group in the docs).

      if(root_processor) then
        
        ! root needs to know the rank (with respect to COMM_CART) of
        ! the remote group's leader. Send this via point-to-point
        ! communication.
        call mpi_recv(remote_leader_rank, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, COMM_CART, status_arr, ierr);
        call mpi_intercomm_create(COMM_ROOT, 0, COMM_CART, &
           & remote_leader_rank, tag_, intercomm, ierr)
      else
        if (cart_rank == 0) then
          ! if this is root (i.e. here used as the leader), then
          ! send its rank (with respect to COMM_CART) to root (with
          ! respect to COMM_CART)
          call mpi_send(processor_number, 1, MPI_INTEGER, 0, tag_, COMM_CART, ierr);
        end if
        call mpi_intercomm_create(comm, 0, COMM_CART, &
           & 0, tag, intercomm, ierr)
      end if

      ! now use the obtained intercommunicator for a gather
      if(root_processor) then
        ! if there were other processes in this group, they would
        ! have to specify MPI_PROC_NULL instead of MPI_ROOT

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, MPI_ROOT, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, MPI_ROOT, intercomm, ierr)
      else
        ! the receiver rank parameter here is the root in the remote group

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, 0, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, 0, intercomm, ierr)
      end if

      call mpi_comm_free(intercomm, ierr);
    end if

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)
    
  end if
#else
  global_array = local_array
#endif
  
end subroutine gather_4d_real


!-----------------------------------------------------------------------------
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator).
!----------------------------------------------------------------------------
subroutine gather_4d_cmplx(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm)
  use global, only : gkw_warn_any_proc
  
  integer, parameter :: data_ndims = 4
  complex, dimension(:,:,:,:), intent(out) :: global_array
  complex, dimension(:,:,:,:), intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart
  logical, intent(in) :: root_is_in_comm

#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  complex :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start
#if defined(mpi3)
  ! Just a note: in case GKW is ever going to use MPI3, note that the
  ! types of some arguments change and if this is not respected then
  ! spurious memory corruption happens (e.g. handles may suddenly
  ! become invalid)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif
  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  ! To silence compiler about unsused dummy argument.
  if (to_root_of_commcart) continue

  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
      global_array = local_array
  else
    lshape = shape(local_array)

    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)

    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 4d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if
    
    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if
    
    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPICOMPLEX_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)

    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    !write (*,*) "bytes of a complex number:", elem_bytes
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)
    ! every process sends its displacement to root
    call gather_array(displacements, num_procs, displacements, 1, comm)

    counts = 1
    call mpi_gatherv(local_array, 1, local_subarray_type, &
       & global_array, counts, &
       & displacements, resized_subarray_type, 0, comm, ierr)

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)
    
  end if
#else
  global_array = local_array
#endif
  
end subroutine gather_4d_cmplx

!-----------------------------------------------------------------------------
!> This routine must be called by all processes of the communicator comm.
!>
!> On every process the local array must be given. On rank=0, additional the
!> global_array buffer must be allocated.
!>
!> IMPORTANT: If to_root_of_commcart is .true., then this routine must 
!> additionally be called by the root process (with respect to COMM_CART). The
!> comm and local_array arguments are ignored by root if root is not in the
!> communicator of the other processes.
!>
!> A global array is returned (on the processor with rank = 0
!> in the given communicator).
!----------------------------------------------------------------------------
subroutine gather_6d_real(global_array, &
   & local_array, mpi_dtype, comm, to_root_of_commcart, root_is_in_comm, tag)
  use mpicomms, only : COMM_CART, COMM_ROOT
  use global, only : gkw_warn_any_proc

  integer, parameter :: data_ndims = 6
  real, dimension(:,:,:,:,:,:), &
     & intent(out) :: global_array
  real, dimension(:,:,:,:,:,:), &
     & intent(in) :: local_array
  integer, intent(in) :: mpi_dtype, comm
  logical, intent(in) :: to_root_of_commcart, root_is_in_comm
  integer, intent(in),optional :: tag

#if defined(mpi2)
  integer, dimension(data_ndims) :: lshape
  integer :: status
  real :: elem
  integer :: num_procs, cart_rank
  integer :: local_subarray_type, resized_subarray_type
  integer(kind=mpi_address_kind) :: start

#if defined(mpi3)
  integer(kind=mpi_count_kind) :: lowbound, true_extent, elem_bytes_extent
#else
  integer(kind=mpi_address_kind) :: lowbound, true_extent, elem_bytes_extent
#endif

  integer :: elem_bytes, ierr
  integer, dimension(:), allocatable :: counts, displacements
  integer, dimension(:), allocatable :: starts

  integer :: remote_leader_rank
  integer :: intercomm, tag_
  integer :: status_arr (MPI_STATUS_SIZE)
  call mpi_topo_test(comm, status, ierr)
  if(status == MPI_UNDEFINED .and. root_is_in_comm) then
    ! the trivial case: no parallelisation in any of the s or x directions.
    global_array = local_array
  else
    lshape = shape(local_array)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 0
    end if


    call mpi_comm_size(comm, num_procs, ierr)
    call mpi_comm_rank(comm, cart_rank, ierr)

    if(root_processor .and. num_procs > 256) then
      call gkw_warn_any_proc('gather 3d: there are very many processes &
         & involved, this is probably much less performant than binary MPI-IO')
    end if

    allocate(starts(data_ndims))
    allocate(counts(num_procs))
    if(to_root_of_commcart .and. root_processor) then
      allocate(displacements(num_procs))
    else if(.not. to_root_of_commcart .and. cart_rank == 0) then
      allocate(displacements(num_procs))
    else
      allocate(displacements(1))
    end if
    
    starts = 0
    call mpi_type_create_subarray(data_ndims, lshape, lshape, &
       & starts, MPI_ORDER_FORTRAN, &
       & MPIREAL_X, local_subarray_type, ierr)
    call mpi_type_commit(local_subarray_type, ierr)

    ! get the size of a complex number in bytes
    inquire(iolength=elem_bytes) elem
    !write (*,*) "bytes of a real number:", elem_bytes
    ! implicit type conversion to avoid compiler complaints
    elem_bytes_extent = elem_bytes

    ! resize the size of this type to make calculation of the
    ! displacement easier
    start = 0
    call mpi_type_create_resized(mpi_dtype, start, elem_bytes_extent, &
       & resized_subarray_type, ierr);
    call mpi_type_commit(resized_subarray_type, ierr);

    ! the lower bound returned by this function is the index of the
    ! first non-gap element in this datatype
    call mpi_type_get_true_extent(mpi_dtype, lowbound, true_extent, ierr)
    displacements = int(lowbound/elem_bytes)

    if(.not. to_root_of_commcart .or. root_is_in_comm) then
      ! Collectively gather using an intracommunicator
      
      ! every process sends its displacement to root
      call gather_array(displacements, num_procs, displacements, 1, comm)

      counts = 1
      call mpi_gatherv(local_array, 1, local_subarray_type, &
         & global_array, counts, &
         & displacements, resized_subarray_type, 0, comm, ierr)
    else
      ! That means that the gather is to root of commcart and that
      ! root is not in comm.
      
      ! Collectively gather using an intercommunicator
      
      ! Why this: The problem is that with serial output, only the root
      ! process can write files.
      ! However, the root process does not necessarily work on the
      ! hyperslab which shall be output, for example, one wants to
      ! output the 3d density for species 2, but there is
      ! parallelisation on the species and root works on species
      ! 1. Hence, to output serially, one has to gather over the
      ! communicator COMM_S_NE_X_NE (of processes working on the slab), but
      ! not to the root of that communicator, but to the root of
      ! COMM_CART!
      
      ! To achieve this, MPI also supports inter-communication: communication
      ! between two non-overlapping groups of processes.
      
      ! The group containing a process that initiates an
      ! inter-communication operation is called the "local group"
      ! that is, the sender in a send and the receiver in a receive.

      ! As in intra-communication, the target process is specified
      ! using a (communicator, rank) pair. Unlike
      ! intra-communication, the rank is relative to a second,
      ! remote group (called peer group in the docs).

      if(root_processor) then
        
        ! root needs to know the rank (with respect to COMM_CART) of
        ! the remote group's leader. Send this via point-to-point
        ! communication.
        call mpi_recv(remote_leader_rank, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, COMM_CART, status_arr, ierr);
        call mpi_intercomm_create(COMM_ROOT, 0, COMM_CART, &
           & remote_leader_rank, tag_, intercomm, ierr)
      else
        if (cart_rank == 0) then
          ! if this is root (i.e. here used as the leader), then
          ! send its rank (with respect to COMM_CART) to root (with
          ! respect to COMM_CART)
          call mpi_send(processor_number, 1, MPI_INTEGER, 0, tag_, COMM_CART, ierr);
        end if
        call mpi_intercomm_create(comm, 0, COMM_CART, &
           & 0, tag_, intercomm, ierr)
      end if

      ! now use the obtained intercommunicator for a gather
      if(root_processor) then
        ! if there were other processes in this group, they would
        ! have to specify MPI_PROC_NULL instead of MPI_ROOT

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, MPI_ROOT, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, MPI_ROOT, intercomm, ierr)
      else
        ! the receiver rank parameter here is the root in the remote group

        ! every process sends its displacement to root
        call mpi_gather(displacements, 1, MPI_INTEGER, &
           & displacements, 1, &
           & MPI_INTEGER, 0, intercomm, ierr)

        counts = 1
        call mpi_gatherv(local_array, 1, local_subarray_type, &
           & global_array, counts, &
           & displacements, resized_subarray_type, 0, intercomm, ierr)
      end if

      call mpi_comm_free(intercomm, ierr);
    end if

    deallocate(starts)
    deallocate(counts)
    deallocate(displacements)
    call mpi_type_free(resized_subarray_type, ierr)
    call mpi_type_free(local_subarray_type, ierr)
    
  end if
#else
  global_array = local_array
#endif
  
end subroutine gather_6d_real



!-----------------------------------------------------------------------------
!> Abort routine for this module (would not be necessary if some
!> restructuring is done).
!-----------------------------------------------------------------------------
subroutine mpiabort(abort_message)
  character (len=*),intent(in), optional :: abort_message
#if defined(mpi2)
  integer :: ierrorcode = 0
#endif
  
  if (root_processor) then
    if (present(abort_message)) then
      write (*,'(A,A,A)') ' ABORT -- ', abort_message
    else
      write (*,'(A,A,A)') ' ABORT -- ', '(no message given!)'
    end if
  end if

#if defined(mpi2)
  call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
#endif
  stop 1

end subroutine mpiabort



!-----------------------------------------------------------------------------
!> basic wrapper for mpi_cart_create to allow MPI_COMM_WORLD to be private
!-----------------------------------------------------------------------------
subroutine mpicart_create(n_dims,dims,periodic,reorder,COMM_CART)

  integer, intent(in)  :: n_dims, dims(n_dims)
  logical, intent(in)  :: periodic(n_dims), reorder
  integer, intent(out) :: COMM_CART

#if defined(mpi2)
  call mpi_cart_create(MPI_COMM_WORLD,n_dims,dims(1:n_dims),&
     & periodic(1:n_dims),reorder,COMM_CART,ierr)
#else
  COMM_CART = I_RUN_WITHOUT_MPI
#endif

end subroutine mpicart_create

!-----------------------------------------------------------------------------
!> workaround to copy MPI_COMM_WORLD for grid/setup_proc_coords
!-----------------------------------------------------------------------------
subroutine mpicopy_world(COMM)

  integer, intent(out) :: COMM

#if defined(mpi2)
  COMM = MPI_COMM_WORLD 
#else
  COMM = I_RUN_WITHOUT_MPI
#endif

end subroutine mpicopy_world

!------------------------------------------------------------------------------
!> Pick a communicator to read/write the full fdisi.
!>
!> For large problems with more than about 64 processors, there can be problems
!> when reading/writing restart files in parallel, depending on the filesystem 
!> and the MPI implementation. This was observed on some Cray XT4 machines
!> where it was necessary to modify environment variables to successfully write 
!> restart files. The communciator setup here provides a workaround which must
!> be set at compile time via adjustments to two variables in switches.F90: 
!> Rather than writing to file with all processors simultaneously, this
!> communicator allows sequential writes with subsets of processors, with the 
!> file closed and re-opened between subsets. Increasing the number of subsets 
!> reduces overheads and can allow files to be written successfully. This might
!> also improve the I/O performance on machines where such operations already 
!> work. However, if the number of subsets is too large, the I/O performance 
!> will suffer. The variable min_fdisi_io subsets selects the minimum number 
!> of processor subsets to use. When problems are encountered, it is suggested 
!> to first try setting this value to the value of n_procs_sp, then to 
!> n_procs_sp*n_procs_s if problems persist.  Using larger values will likely
!> give very poor performance, but may be necessary. The other variable, 
!> min_procs_io split, sets the minimum number of processors for which the 
!> read/write split should be used.
!------------------------------------------------------------------------------
subroutine set_comm_fdisi_io()

  use mpicomms, only : COMM_FDISI_IO, COMM_FDISI_IO_nprocsets
  use mpicomms, only : COMM_FDISI_IO_procset, COMM_SP_EQ, COMM_S_EQ
  use mpicomms, only : COMM_SP_EQ_nprocsets, COMM_SP_EQ_procset
  use mpicomms, only : COMM_CART, COMM_SP_EQ_S_EQ, COMM_SP_EQ_S_EQ_nprocsets
  use mpicomms, only : COMM_SP_EQ_S_EQ_procset, COMM_ALL_EQ
  use mpicomms, only : COMM_S_EQ_nprocsets, COMM_S_EQ_procset
  use switches, only : min_fdisi_io_subsets, min_procs_io_split

#if defined(mpi2)
  integer :: my_key,my_color,ierr
#endif

  if (min_fdisi_io_subsets <= 1 .or.                                         &
      & number_of_processors < min_procs_io_split) then

    ! all processors involved in 1 write
    COMM_FDISI_IO           = COMM_CART
    COMM_FDISI_IO_nprocsets = 1
    COMM_FDISI_IO_procset   = 1

  else if (COMM_SP_EQ_nprocsets >= min_fdisi_io_subsets) then

    ! split writes over different parts of the species
    COMM_FDISI_IO           = COMM_SP_EQ
    COMM_FDISI_IO_nprocsets = COMM_SP_EQ_nprocsets
    COMM_FDISI_IO_procset   = COMM_SP_EQ_procset

  else if (COMM_SP_EQ_nprocsets == 1 .and.                                   &
      & COMM_S_EQ_nprocsets >= min_fdisi_io_subsets) then

    ! split writes over different parts of the s-grid
    COMM_FDISI_IO           = COMM_S_EQ
    COMM_FDISI_IO_nprocsets = COMM_S_EQ_nprocsets
    COMM_FDISI_IO_procset   = COMM_S_EQ_procset

  else if (COMM_SP_EQ_S_EQ_nprocsets >= min_fdisi_io_subsets) then

    ! split writes over different parts of the species and s-grid
    COMM_FDISI_IO           = COMM_SP_EQ_S_EQ
    COMM_FDISI_IO_nprocsets = COMM_SP_EQ_S_EQ_nprocsets
    COMM_FDISI_IO_procset   = COMM_SP_EQ_S_EQ_procset
#if defined(mpi2)
  else if (min_fdisi_io_subsets <= number_of_processors / 2) then

    !
    ! In this case, create some new communicators for now...
    !

    ! pick a comm which I would like to belong to
    my_color = mod(processor_number,min_fdisi_io_subsets)

    ! For this purpose I don't care what rank I get in the new comm.
    my_key = 0

    ! split comm_world
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,my_color,my_key,COMM_FDISI_IO,ierr)

    ! use the prescribed number of subsets
    COMM_FDISI_IO_nprocsets = min_fdisi_io_subsets

    ! pick my procset using the coloring
    COMM_FDISI_IO_procset = my_color + 1

#endif
  else

    ! only the local processor at a time - will be very slow!
    COMM_FDISI_IO           = COMM_ALL_EQ
    COMM_FDISI_IO_nprocsets = number_of_processors
    COMM_FDISI_IO_procset   = processor_number + 1

  end if

end subroutine set_comm_fdisi_io

  !---------------------------------------------------------------------------
  !> This function checks if a scalar integer is the same on all processors.
  !---------------------------------------------------------------------------
  function chk_eq_integer(param) result(is_eq)
    integer, intent(in) :: param
    logical :: is_eq
    integer :: param_ref

    ! obtain the value of an arbitrary other process as a reference
    param_ref = param
    call mpibcast(param_ref, 1)
    ! compare equality across processes
    call mpiallreduce_and(param == param_ref, &
       & is_eq,1)
  end function chk_eq_integer

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  function chk_eq_real(param) result(is_eq)
    use global, only : r_tiny
    real, intent(in) :: param
    logical :: is_eq
    real :: param_ref

    ! obtain the value of an arbitrary other process as a reference
    param_ref = param
    call mpibcast(param_ref, 1)
    ! compare equality across processes
    call mpiallreduce_and(abs(param - param_ref) < r_tiny, &
       & is_eq,1)

  end function chk_eq_real

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  function chk_eq_complex(param) result(is_eq)
    use global, only : r_tiny
    complex, intent(in) :: param
    logical :: is_eq
    complex :: param_ref

    ! obtain the value of an arbitrary other process as a reference
    param_ref = param
    call mpibcast(param_ref, 1)
    ! compare equality across processes
    call mpiallreduce_and(abs(real(param - param_ref)) < r_tiny .and. &
       & abs(aimag(param-param_ref)) < r_tiny, &
       & is_eq,1)

  end function chk_eq_complex

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  function chk_eq_string(param) result(is_eq)
    character(len=*), intent(in) :: param
    logical :: is_eq
    character(len=len(param)) :: param_ref

    ! obtain the value of an arbitrary other process as a reference
    param_ref = param
    call mpibcast(param_ref, len(param_ref))
    ! compare equality across processes
    call mpiallreduce_and(param .eq. param_ref, &
       & is_eq,1)

  end function chk_eq_string

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  function chk_eq_i_array(param) result(is_eq)
    integer, dimension(:), intent(in) :: param
    logical :: is_eq
    integer, dimension(size(param,1)) :: param_ref

    ! obtain the value of an arbitrary other process as a reference
    param_ref = param
    call mpibcast(param_ref, size(param_ref,1))
    ! compare equality across processes
    call mpiallreduce_and(all(param == param_ref,1), &
       & is_eq,1)

  end function chk_eq_i_array

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  function chk_eq_r_array(param) result(is_eq)
    use global, only : r_tiny
    real, dimension(:), intent(in) :: param
    logical :: is_eq
    real, dimension(size(param,1)) :: param_ref

    ! obtain the value of an arbitrary other process as a reference
    param_ref = param
    call mpibcast(param_ref, size(param_ref,1))
    ! compare equality across processes
    call mpiallreduce_and(all(abs(param - param_ref) < r_tiny,1), &
       & is_eq,1)

  end function chk_eq_r_array

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  function chk_eq_logical(param) result(is_eq)
    logical, intent(in) :: param
    logical :: is_eq
    logical :: param_ref

    ! obtain the value of an arbitrary other process as a reference
    param_ref = param
    call mpibcast(param_ref, 1)
    ! compare equality across processes
    call mpiallreduce_and(param .eqv. param_ref, &
       & is_eq,1)

  end function chk_eq_logical


  !---------------------------------------------------------------------------
  !> This routine can help to prevent MPI tag collisions, without
  !> hardcoding anything. Every part of the code which needs
  !> to use MPI tags (e.g. to avoid dataraces when data
  !> is sent to root by different processes using send_to_root)
  !> should call this initially, to make sure that no other part of the code
  !> will use those same tags.
  !---------------------------------------------------------------------------
  subroutine register_tag_range(range_size, range_start, range_end_inkl)
    integer, intent(in) :: range_size
    integer, intent(out) :: range_start, range_end_inkl

    range_start = tag_range_upper_limit + 1
    range_end_inkl = range_start + range_size
    tag_range_upper_limit = range_end_inkl

  end subroutine register_tag_range
  
  !---------------------------------------------------------------------------
  !> This routine must be called by both the sending process and root
  !> of the given communicator. It is illegal to call this routine with more than
  !> two processes.
  !---------------------------------------------------------------------------
  subroutine send_scalar_int(data, comm, tag)
    integer, intent(inout) :: data
    integer, intent(in), optional :: comm, tag
    integer :: nelem
#if defined(mpi2)
    integer :: status_arr (MPI_STATUS_SIZE)
    integer :: c, tag_
    if(present(comm)) then
      c = comm
    else
      c = MPI_COMM_WORLD
    end if

    nelem = 1

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 1
    end if
    if (nelem > 0) then
      if(root_processor) then
        ! Receive the array via point-to-point communication. As the
        ! rank of the sending process is unknown, use MPI_ANY_SOURCE.
        call mpi_recv(data, nelem, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, c, status_arr, ierr);
      else
        call mpi_send(data, nelem, MPI_INTEGER, 0, tag_, c, ierr);
      end if
    end if
#endif
  end subroutine send_scalar_int

  !---------------------------------------------------------------------------
  !> initialise the type used to communicate matrix elements.
  !---------------------------------------------------------------------------
  subroutine create_matrix_element_type()
    integer, parameter :: nblocks = 4
    integer, parameter, dimension(nblocks) :: array_of_blocklengths = &
       & (/ 15, 1, 1, 64/)
    integer(kind=MPI_ADDRESS_KIND), dimension(nblocks) :: &
       & array_of_displacements
    integer, dimension(nblocks) :: array_of_types = &
       & (/ MPI_INTEGER, MPI_LOGICAL, MPI_COMPLEX, MPI_CHARACTER /)

    integer(kind=MPI_ADDRESS_KIND) :: bytes_of_int, bytes_of_complex
    integer(kind=MPI_ADDRESS_KIND) :: lbound, bytes_of_logical
    integer(kind=MPI_ADDRESS_KIND) :: bytes_of_real
    call mpi_type_get_extent(MPI_INTEGER, lbound, bytes_of_int, ierr)
    call mpi_type_get_extent(MPI_COMPLEX, lbound, bytes_of_complex, ierr)
    call mpi_type_get_extent(MPI_REAL, lbound, bytes_of_real, ierr)
    call mpi_type_get_extent(MPI_LOGICAL, lbound, bytes_of_logical, ierr)

    array_of_displacements(1) = 0
    array_of_displacements(2) = bytes_of_int*array_of_blocklengths(1)&
       & +array_of_displacements(1)
    array_of_displacements(3) = bytes_of_logical*array_of_blocklengths(2)&
       & +array_of_displacements(2)
    array_of_displacements(4) = bytes_of_complex*array_of_blocklengths(3)&
       & +array_of_displacements(3)

    call mpi_type_create_struct(nblocks, array_of_blocklengths, &
       & array_of_displacements, array_of_types, matrix_element_type, ierr)
    call mpi_type_commit(matrix_element_type, ierr)

  end subroutine create_matrix_element_type

  !---------------------------------------------------------------------------
  !> This routine must be called by both the sending process and root
  !> of the given communicator. It is illegal to call this routine with more than
  !> two processes.
  !---------------------------------------------------------------------------
  subroutine send_scalar_matrix_element(data, comm, tag)
    use structures, only : matrix_element
    type(matrix_element), intent(inout) :: data
    integer, intent(in), optional :: comm, tag
#if defined(mpi2)
    integer :: status_arr (MPI_STATUS_SIZE)
    integer :: c, tag_

    if(present(comm)) then
      c = comm
    else
      c = MPI_COMM_WORLD
    end if

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 1
    end if

    if(root_processor) then
      ! Receive the array via point-to-point communication. As the
      ! rank of the sending process is unknown, use MPI_ANY_SOURCE.
      call mpi_recv(data, 1, matrix_element_type, MPI_ANY_SOURCE, &
         & tag_, c, status_arr, ierr);
    else
      call mpi_send(data, 1, matrix_element_type, 0, tag_, c, ierr);
    end if
#endif
  end subroutine send_scalar_matrix_element

  !---------------------------------------------------------------------------
  !> This routine must be called by both the sending process and root
  !> of the given communicator. It is illegal to call this routine with more than
  !> two processes.
  !---------------------------------------------------------------------------
  subroutine send_1d_complex(data, dims, comm, tag)
    integer, dimension(:), intent(in) :: dims
    complex, dimension(dims(1)), intent(inout)  :: data
    integer, intent(in), optional :: comm, tag
    integer :: nelem
#if defined(mpi2)
    integer :: status_arr (MPI_STATUS_SIZE)
    integer :: c, tag_
    if(present(comm)) then
      c = comm
    else
      c = MPI_COMM_WORLD
    end if
    
    nelem = product(dims)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 1
    end if
    if (nelem > 0) then
      if(root_processor) then
        ! Receive the array via point-to-point communication. As the
        ! rank of the sending process is unknown, use MPI_ANY_SOURCE.
        call mpi_recv(data, nelem, MPICOMPLEX_X, MPI_ANY_SOURCE, &
           & tag_, c, status_arr, ierr);
      else
        call mpi_send(data, nelem, MPICOMPLEX_X, 0, tag_, c, ierr);
      end if
    end if
#endif
  end subroutine send_1d_complex


  !---------------------------------------------------------------------------
  !> This routine must be called by both the sending process and root
  !> of the given communicator. It is illegal to call this routine with more than
  !> two processes.
  !---------------------------------------------------------------------------
  subroutine send_1d_int(data, dims, comm, tag)
    integer, dimension(:), intent(in) :: dims
    integer, dimension(dims(1)), intent(inout)  :: data
    integer, intent(in), optional :: comm, tag
    integer :: nelem
#if defined(mpi2)
    integer :: status_arr (MPI_STATUS_SIZE)
    integer :: c, tag_
    if(present(comm)) then
      c = comm
    else
      c = MPI_COMM_WORLD
    end if
    
    nelem = product(dims)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 1
    end if
    if (nelem > 0) then
      if(root_processor) then
        ! Receive the array via point-to-point communication. As the
        ! rank of the sending process is unknown, use MPI_ANY_SOURCE.
        call mpi_recv(data, nelem, MPI_INTEGER, MPI_ANY_SOURCE, &
           & tag_, c, status_arr, ierr);
      else
        call mpi_send(data, nelem, MPI_INTEGER, 0, tag_, c, ierr);
      end if
    end if
#endif
  end subroutine send_1d_int


  !---------------------------------------------------------------------------
  !> This routine must be called by both the sending process and root
  !> of the given communicator. It is illegal to call this routine with more than
  !> two processes.
  !---------------------------------------------------------------------------
  subroutine send_2d_real(data, dims, comm, tag)
    integer, dimension(:), intent(in) :: dims
    real, dimension(dims(1), dims(2)), intent(inout)  :: data
    integer, optional, intent(in) :: comm, tag
    integer :: nelem
#if defined(mpi2)
    integer :: status_arr (MPI_STATUS_SIZE)
    integer :: c, tag_
    if(present(comm)) then
      c = comm
    else
      c = MPI_COMM_WORLD
    end if
    
    nelem = product(dims)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 0
    end if
    if (nelem > 0) then
      if(root_processor) then
        ! Receive the array via point-to-point communication. As the
        ! rank of the sending process is unknown, use MPI_ANY_SOURCE.
        call mpi_recv(data, nelem, MPIREAL_X, MPI_ANY_SOURCE, &
           & tag_, c, status_arr, ierr);
      else
        call mpi_send(data, nelem, MPIREAL_X, 0, tag_, c, ierr);
      end if
    end if
#endif
  end subroutine send_2d_real

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine send_4d_complex(data, dims, comm, tag)
    integer, dimension(:), intent(in) :: dims
    complex, dimension(dims(1), dims(2), dims(3), dims(4)), intent(inout)  :: data
    integer, optional, intent(in) :: comm, tag
    integer :: nelem
#if defined(mpi2)
    integer :: status_arr (MPI_STATUS_SIZE)
    integer :: c, tag_
    if(present(comm)) then
      c = comm
    else
      c = MPI_COMM_WORLD
    end if
    
    nelem = product(dims)

    if(present(tag)) then
      tag_ = tag
    else
      tag_ = 1
    end if
    if (nelem > 0) then
      if(root_processor) then
        ! Receive the array via point-to-point communication. As the
        ! rank of the sending process is unknown, use MPI_ANY_SOURCE.
        call mpi_recv(data, nelem, MPICOMPLEX_X, MPI_ANY_SOURCE, &
           & tag_, c, status_arr, ierr);
      else
        call mpi_send(data, nelem, MPICOMPLEX_X, 0, tag_, c, ierr);
      end if
    end if
#endif
  end subroutine send_4d_complex

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine send_4d_real(data, dims, comm, tag)
    integer, dimension(:), intent(in) :: dims
    real, dimension(dims(1), dims(2), dims(3), dims(4)), intent(inout)  :: data
    integer, optional, intent(in) :: comm, tag
    integer :: nelem
#if defined(mpi2)
    integer :: status_arr (MPI_STATUS_SIZE)
    integer :: c, t
    if(present(comm)) then
      c = comm
    else
      c = MPI_COMM_WORLD
    end if

    nelem = product(dims)

    if(present(tag)) then
      t = tag
    else
      t = 1
    end if
    if (nelem > 0) then
      if(root_processor) then
        ! Receive the array via point-to-point communication. As the
        ! rank of the sending process is unknown, use MPI_ANY_SOURCE.
        call mpi_recv(data, nelem, MPIREAL_X, MPI_ANY_SOURCE, &
           & t, c, status_arr, ierr);
      else
        call mpi_send(data, nelem, MPIREAL_X, 0, t, c, ierr);
      end if
    end if
#endif
  end subroutine send_4d_real

end module mpiinterface

