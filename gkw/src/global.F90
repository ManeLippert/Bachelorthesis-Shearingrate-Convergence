!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Various unchanging parameters needed in most parts of the code.
!> Furthermore some very small functions and subroutines which do not depend
!> on any other module and are also needed in many parts of the code.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module global

  implicit none

  private

  !
  ! Parameters for capabilities.
  !
#ifdef HAVE_SLEPC
  logical, parameter, public :: compiled_with_slepc = .true.
#else
  logical, parameter, public :: compiled_with_slepc = .false.
#endif

#ifdef HAVE_HDF5
  logical, parameter, public :: compiled_with_hdf5 = .true.
#else
  logical, parameter, public :: compiled_with_hdf5 = .false.
#endif

#ifdef HAVE_LIBRSB
  logical, parameter, public :: compiled_with_librsb = .true.
#else
  logical, parameter, public :: compiled_with_librsb = .false.
#endif

#ifdef HAVE_MKL
  logical, parameter, public :: compiled_with_mkl = .true.
#else
  logical, parameter, public :: compiled_with_mkl = .false.
#endif

#if !defined(mpi2)
  !> parameter for setting some dummy integers when the code runs without MPI
  integer, parameter, public :: I_RUN_WITHOUT_MPI = -64523
  logical, parameter, public :: compiled_with_mpi = .false.
#else
  logical, parameter, public :: compiled_with_mpi = .true.
#endif

#ifdef _OPENMP
  !> reporting label
  logical, parameter, public :: compiled_with_openmp = .true.
  integer, parameter, public :: openmp_date = _OPENMP
#else
  !> reporting label
  logical, parameter, public :: compiled_with_openmp = .false.
  integer, parameter, public :: openmp_date = 0
#endif

 !> reporting label
#ifdef umfpack
#ifdef kernel32bit
  integer, parameter, public :: compiled_with_umfpack = 32
#else
  integer, parameter, public :: compiled_with_umfpack = 64
#endif
#else
   integer, parameter, public :: compiled_with_umfpack = 0
#endif


  !
  ! Precision
  !
  !> i8 is for long integer
  integer, parameter, public :: i8 = selected_int_kind(R=18)
  !> dp is for doubles
  integer, parameter, public :: dp = kind(0.d0)
  !> rp is the default real precision which should be used almost everywhere

  ! Umfpack pointers
#ifdef kernel32bit
  integer, parameter, public :: iumf = kind(int(1))
#else
  !This could use the int64  ISO_FORTRAN_ENV but not until fortran 2008
  integer, parameter, public :: iumf = i8
#endif

  ! Umfpack always works in double precision 
  ! if the rest of the code is in single precision, conversions are used
  ! all type conversions should occur implicitly in the conversion
  ! routine in matconv.f90, and where the outputs of 
  ! umf4sol are copied back to GKW arrays (in imp_integration and fields) 
  integer, parameter, public :: rumf = kind(0.d0)

#ifdef real_precision_default
  integer, parameter, public :: rp = kind(0.e0)
  integer, parameter, public :: cp = kind((0.e0,0.e0))
  ! Kind numbers for single precision integer container
  integer, parameter :: Single = selected_int_kind(precision(1.e0))
  ! Single precision IEEE NaN bit value
  integer(Single), save, public :: iNaN  !  = Z"7FC00000"
  data iNaN/Z"7FC00000"/
#else
  integer, parameter, public :: rp = kind(0.d0)
  integer, parameter, public :: cp = kind((0.d0,0.d0))
  ! Kind numbers for double precision integer container
  integer, parameter :: Double = selected_int_kind(precision(1.d0))
  ! double precision IEEE NaN bit value
  integer(Double), save, public :: iNaN  !  = Z"7FF8000000000000"
  data iNaN/Z"7FF8000000000000"/
#endif

  ! Determine real_precision from rp
  logical, parameter, public :: have_double_precision = (rp == dp)

  !
  ! Error and message handling:
  !
  integer, parameter, public :: LEVEL_NORMAL_OPERATION = 0
  integer, parameter, public :: LEVEL_INFO = 10
  integer, parameter, public :: LEVEL_WARN = 20
  integer, parameter, public :: LEVEL_ERROR = 30
  integer, parameter, public :: LEVEL_ABORT = 40
  type :: error
    logical :: err = .false.
    integer :: level = LEVEL_NORMAL_OPERATION
    character (len=256) :: msg = ''
    ! type(error) , pointer :: next => null()
  end type error
  type(error), save :: latest_error

  integer, parameter, public :: PHI_FIELD = 1, APAR_FIELD = 2, BPAR_FIELD = 3
  integer, parameter, public :: DISTRIBUTION = 4, PHI_GA_FIELD = 5
  integer, parameter, public :: APAR_GA_FIELD = 6, BPAR_GA_FIELD = 7
  integer, parameter, public :: EVERY_FIELD = 0


  !
  ! Identifiers for coordinate directions which are used at various
  ! places, e.g. in the index function module:
  !
  !> An id to indicate vpar (parallel velocity direction)
  integer, parameter, public :: id_vpar  = 1
  !> An id to indicate mu (magn. moment direction)
  integer, parameter, public :: id_mu    = 2
  !> An id to indicate s (parallel direction)
  integer, parameter, public :: id_s     = 3
  !> An id to indicate y modes (toroidal modes)
  integer, parameter, public :: id_mod   = 4
  !> An id to indicate x (radial direction)
  integer, parameter, public :: id_x     = 5
  !> An id to indicate species direction
  integer, parameter, public :: id_sp    = 6
  !> dummy id
  integer, parameter, public :: id_dummy = -7

  !
  ! Various other settings:
  !

  !> parameter which is false
  logical, parameter, public :: logical_false = .false.
  !> length for characters
  integer, parameter, public :: lenswitch = 32
  !> large integer
  integer, parameter, public :: i_huge = huge(1)/5
  !> some large integer
  integer, parameter, public :: i_huge_tag = 62768-huge(1)/5
  !> large real
  real (kind=rp), parameter, public :: r_huge = 0.2_rp*huge(1._rp)
  !> small real
  real (kind=rp), parameter, public :: r_tiny = 5._rp*tiny(1._rp)

  logical, save, public :: lverbose = .false.
  logical, save, public :: root_and_verbose
  
  interface Kahan_sum_step
    module procedure Kahan_sum_step_r
    module procedure Kahan_sum_step_c
  end interface Kahan_sum_step

  interface Kahan_sum
    module procedure Kahan_sum_r
    module procedure Kahan_sum_c
  end interface Kahan_sum

  !
  ! Global small helper functions and routines which do *not depend* on any
  ! other module.
  !

  public :: int2char, int2char_zeros, real2char, ss, dotdat
  ! Error handling:
  public :: gkw_check_err, gkw_its_critical, gkw_throw_abort, gkw_throw_error
  public :: gkw_throw_warn, gkw_throw_info
  public :: gkw_warn_any_proc

  public :: gkw_a_equal_b_accuracy

  interface gkw_a_equal_b_accuracy
    module procedure gkw_a_equal_b_accuracy_real
    module procedure gkw_a_equal_b_accuracy_complex
  end interface gkw_a_equal_b_accuracy

  public ::Kahan_sum, Kahan_sum_step

contains
  
  !-----------------------------------------------------------------------------
  !> Take an integer, i, and write it into a character of length j. The total
  !> length is always 8, so it should usually be trimmed.
  !------------------------------------------------------------------------------
  function int2char(number,ndigits)

    integer, intent(in)           :: number
    integer, optional, intent(in) :: ndigits

    character (len=8) :: int2char
    character (len=6) :: ifmt

    if (present(ndigits)) then
      select case(ndigits)
      case(0) ; ifmt='(I4)'
      case(1) ; ifmt='(I1)'
      case(2) ; ifmt='(I2)'
      case(3) ; ifmt='(I3)'
      case(4) ; ifmt='(I4)'
      case(5) ; ifmt='(I5)'
      case(6) ; ifmt='(I6)'
      case(7) ; ifmt='(I7)'
      case default ; ifmt='(I8)'
      end select
      !write(int2char,'(A8)') '        '
    else
      ifmt='(I8)'
    end if

    write(int2char,ifmt) number

  end function int2char

  function int2char_zeros(number, ndigits)
    integer, intent(in)           :: number
    integer, optional, intent(in) :: ndigits

    character (len=8) :: int2char_zeros
    character (len=6) :: ifmt

    if (present(ndigits)) then
      select case(ndigits)
      case(0) ; ifmt='(I4)'
      case(1) ; ifmt='(I1.1)'
      case(2) ; ifmt='(I2.2)'
      case(3) ; ifmt='(I3.3)'
      case(4) ; ifmt='(I4.4)'
      case(5) ; ifmt='(I5.5)'
      case(6) ; ifmt='(I6.6)'
      case(7) ; ifmt='(I7.7)'
      case default ; ifmt='(I8.8)'
      end select
      write(int2char_zeros,'(A8)') '        '
    else
      ifmt='(I8.8)'
    end if

    write(int2char_zeros,ifmt) number

  end function int2char_zeros

  !-----------------------------------------------------------------------------
  !> Take a real number, and write it into a string. This is meant for
  !> use in gkw_warn.
  !------------------------------------------------------------------------------
  function real2char(number)

    real, intent(in)           :: number

    character (len=13) :: real2char

    write(real2char,'(es13.5)') number

  end function real2char

  !---------------------------------------------------------------------------
  !> A tiny helper function: Naturally one would call this
  !> 'strip_suffix', but in order not to lengthen much the lines where
  !> this function is called, its name is very short.
  !>
  !> this function strips '.dat' from 'foobar.dat'
  !> and '.kyspec' from 'foofoo.kyspec'
  !> and so on..
  !---------------------------------------------------------------------------
  function ss(str)
    !> the input string
    character(len=*), intent(in) :: str
    !> The return value is the input string without the trailing
    !> suffix, if there was any.
    character(len=len(str)) :: ss
    integer :: i

    ! find the position of the dot, from the end of the string
    i = index(str, '.', .true.)
    !i = index(str, '.dat', .true.)

    ! strip everything after the last dot, including the dot.
    if(i > 0) then
      ss = str(1:i-1)
    else
      ss = str
    end if

    ! replace all remaining dots with underscores
    do i = 1, len(ss)
      if(ss(i:i) == '.') then
        ss(i:i) = '_'
      end if
    end do

  end function ss

  function dotdat(str, is_legacy)
    character(len=*), intent(in) :: str
    logical, intent(in) :: is_legacy
    character(len=len(str)+4) :: dotdat
    if(is_legacy) then
      dotdat = trim(str)//'.dat'
    else
      dotdat = trim(str)
    end if
  end function dotdat


  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine gkw_throw(level, message)
    integer, intent(in) :: level
    character(len=*), intent(in) :: message

    if(level > latest_error%level) then
      latest_error%err = .true.
      latest_error%level = level
      latest_error%msg = message
    end if
  end subroutine gkw_throw

  !---------------------------------------------------------------------------
  !> This function checks, if the state of the error struct signals
  !> that an error has recently occurred. It prints the corresponding message.
  !>
  !> The return value is .true. if GKW shall be aborted,
  !> and .false. otherwise.
  !> The intention behind this is that the functions to call to properly abort
  !> a run are part of the IO, MPI and OpenMP modules. Errors can occur in
  !> any of them, but they cannot depend all from each other.
  !> 
  !---------------------------------------------------------------------------
  function gkw_check_err()
    logical :: gkw_check_err
    
    ! ToDo: Print all messages on a fifo queue
    ! Now: Print the (latest) error message
    ! and cause an abort, if requested.

    gkw_check_err = .false.
    
    if(latest_error%err) then
      select case(latest_error%level)
      case(LEVEL_INFO)
        write (*, '(A)') latest_error%msg
      case(LEVEL_WARN)
        call gkw_warn_any_proc(latest_error%msg)
      case(LEVEL_ERROR)
        write (*, '(A,A)') '*** ERROR (not critical) *** ', latest_error%msg
      case(LEVEL_ABORT)
        write (*, '(A,A)') '*** ABORT *** ', latest_error%msg
        gkw_check_err = .true.
      end select
    end if

    latest_error%err = .false.
    latest_error%level = LEVEL_NORMAL_OPERATION
    latest_error%msg = ''
    
  end function gkw_check_err

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  function gkw_its_critical()
    logical :: gkw_its_critical
    
    gkw_its_critical = .false.
    
    if(latest_error%err) then
      select case(latest_error%level)
      case(LEVEL_ABORT)
        gkw_its_critical = .true.
      end select
    end if
    
  end function gkw_its_critical

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine gkw_throw_abort(message)
    character(len=*), intent(in) :: message
    call gkw_throw(LEVEL_ABORT, message)
  end subroutine gkw_throw_abort

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine gkw_throw_error(message)
    character(len=*), intent(in) :: message
    call gkw_throw(LEVEL_ERROR, message)
  end subroutine gkw_throw_error

  !---------------------------------------------------------------------------
  !> This routine prints a warning message (if called by any processor, as
  !> opposed to gkw_warn in the module general).
  !---------------------------------------------------------------------------
  subroutine gkw_throw_warn(message)
    character(len=*), intent(in) :: message
    call gkw_throw(LEVEL_WARN, message)
  end subroutine gkw_throw_warn

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine gkw_throw_info(message)
    character(len=*), intent(in) :: message
    call gkw_throw(LEVEL_INFO, message)
  end subroutine gkw_throw_info

  !-----------------------------------------------------------------------------
  !> This routine prints a warning message. Any process who calls this spits
  !> out a message.
  !-----------------------------------------------------------------------------
  subroutine gkw_warn_any_proc(message)  
    character(len=*), intent(in) :: message 

    write (*,'(A,A)') '*** WARNING *** ', message

  end subroutine gkw_warn_any_proc

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  elemental function gkw_a_equal_b_accuracy_real(a, b, accuracy)
    real, intent(in) :: a, b
    real, optional, intent(in) :: accuracy

    logical :: gkw_a_equal_b_accuracy_real

    real :: local_accuracy

    if (present(accuracy)) then
      local_accuracy = accuracy
    else
      local_accuracy = r_tiny
    end if

    gkw_a_equal_b_accuracy_real = (abs(a - b) < local_accuracy)
  end function gkw_a_equal_b_accuracy_real

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  elemental function gkw_a_equal_b_accuracy_complex(a, b, accuracy)
    complex, intent(in) :: a, b
    real, optional, intent(in) :: accuracy

    logical :: gkw_a_equal_b_accuracy_complex

    real :: local_accuracy

    if (present(accuracy)) then
      local_accuracy = accuracy
    else
      local_accuracy = r_tiny
    end if

    gkw_a_equal_b_accuracy_complex = (abs(a - b) < local_accuracy)
  end function gkw_a_equal_b_accuracy_complex

  !> Implements a single step of the Kahan algorithm
  !> https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  !> See also http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
  !> Written with capital K, because the algorithm is named after a person.
  !>
  !> Use this if you do not have an input array, e.g. because the elements are
  !> computed.
  !> Example:
  !> real :: c = 0.0, kahan = 0.0
  !> do i = 1, 100
  !>   do j = 1, 100
  !>     call Kahan_sum_step(a(i,j)*b(j)*c*x(i,j), kahan, c)
  !>   end do
  !> end do
  !---------------------------------------------------------------------------
  subroutine Kahan_sum_step_r(input, Kahan_sum, compensator)
    real, intent(in) :: input
    real, intent(inout) :: Kahan_sum, compensator

    ! Help/temporary variables.
    real :: intermediate_sum
    real :: intermediate_element
    intermediate_sum = 0.0
    intermediate_element = 0.0

    intermediate_element = input - compensator
    intermediate_sum = Kahan_sum + intermediate_element
    compensator = (intermediate_sum - Kahan_sum) - intermediate_element
    Kahan_sum = intermediate_sum
  end subroutine Kahan_sum_step_r

  subroutine Kahan_sum_step_c(input, Kahan_sum, compensator)
    complex, intent(in) :: input
    complex, intent(inout) :: Kahan_sum, compensator

    ! Help/temporary variables.
    complex :: intermediate_sum
    complex :: intermediate_element
    intermediate_sum = 0.0
    intermediate_element = 0.0

    intermediate_element = input - compensator
    intermediate_sum = Kahan_sum + intermediate_element
    compensator = (intermediate_sum - Kahan_sum) - intermediate_element
    Kahan_sum = intermediate_sum
  end subroutine Kahan_sum_step_c

  !---------------------------------------------------------------------------
  !> Implementing https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  !> see also http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
  !> Written with capital K, because the algorithm is named after a person.
  !>
  !> Use this if you want to sum over an array with improved accuracy.
  !> Example:
  !> real, dimension(100) :: a
  !> a = ...
  !> write(*,*) Kahan_sum(a)
  !---------------------------------------------------------------------------
  function Kahan_sum_r(input)
    real, dimension(:), intent(in) :: input
    real :: Kahan_sum_r
    ! Keeps track of fraction that are lost in summing.
    real :: compensator
    integer :: i
    compensator = 0.0
    Kahan_sum_r = 0.0

    do i = 1, size(input)
       call Kahan_sum_step_r(input(i), Kahan_sum_r, compensator)
    end do

  end function Kahan_sum_r

  function Kahan_sum_c(input)
    complex, dimension(:), intent(in) :: input
    complex :: Kahan_sum_c

    ! Keeps track of fraction that are lost in summing.
    complex :: compensator
    integer :: i
    compensator = 0.0
    Kahan_sum_c = 0.0

    do i = 1, size(input)
       call Kahan_sum_step_c(input(i), Kahan_sum_c, compensator)
    end do

  end function Kahan_sum_c

end module global
