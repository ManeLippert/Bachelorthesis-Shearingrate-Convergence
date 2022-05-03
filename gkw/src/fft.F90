!----------------------------------------------------------------------
!> Contains generalised FFT interfaces for fftw3 and fftw3
!> Also contains a lot of stuff that is not currently used by GKW
!----------------------------------------------------------------------
module fft
#if defined(FFT_MKL)
  ! UPPERCASE USE to deliberately exclude from mkdeps script
  USE mkl_dfti
#endif
  
  use global, only : i8
  implicit none

  private

  public :: fourcol, four1D_real, four2D_real

  integer, parameter, public :: FFT_FORWARD = -1
  integer, parameter, public :: FFT_INVERSE = 1

  ! define the maximum number of plans/work arrays/descriptor handles.
  integer, parameter :: MXPLAN = 64
  integer, parameter :: N_FFT_ROUTINES = 9

  interface fourcol
    module procedure fourcol_ca
    module procedure fourcol_caa
    module procedure fourcol_ra_ca
    module procedure fourcol_raa_caa
  end interface fourcol

#if defined(FFT_FFTW3)
!=============================================================================
! FFTW v3
!=============================================================================

  !> if this interface compiles, FFTW should work
  logical, parameter, public :: WORKING_FFT_LIBRARY=.true.
  character (len=5), parameter, public :: fftlib = 'FFTW3'

  ! UPPERCASE INCLUDE to deliberately exclude from mkdeps script
  INCLUDE 'fftw3.f'

  !> plans for FFT transforms
  integer (kind=i8), dimension(MXPLAN,N_FFT_ROUTINES), save :: plans
  integer, dimension(MXPLAN,N_FFT_ROUTINES,4), save :: plan_properties = 1
  !> number of plans saved
  integer, dimension(N_FFT_ROUTINES), save :: nplans=0

contains

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  pure function get_plan_id(arr_shape, isign, num) result(id)
    integer, dimension(:), intent(in) :: arr_shape
    integer, intent(in) :: isign, num
    integer :: id
    integer :: i, j
    id = -1

    search_plan: do i = 1, nplans(num)
      if (plan_properties(i,num,1) == isign) then
        do j = 1,size(arr_shape)
          if(plan_properties(i,num,1+j) /= arr_shape(j)) cycle search_plan
        end do
        do j = size(arr_shape)+1,3
          if(plan_properties(i,num,1+j) /= 1) cycle search_plan
        end do
        ! a plan is considered to fit if it was created for the same
        ! trafo direction and for a real-space array of the same size
        id = i
        return
      end if
    end do search_plan

  end function get_plan_id

  !--------------------------------------------------------------------
  !> \attention Use of complex-to-real FFT (inverse/backeard) overwrites the
  !>   input array (vec_ca).
  !>   Using FFTW_PRESERVE_INPUT, will keep the input intact, but reduces
  !>   performance.
  !--------------------------------------------------------------------
  subroutine four1D_real(vec_ra, vec_ca, isign)
    real,    dimension(:), intent(inout) :: vec_ra
    complex, dimension(:), intent(inout) :: vec_ca
    integer, intent(in) :: isign

    integer, parameter :: NUM=1

    integer :: k
    integer :: dim1_ra, dim1_ca, id, istat
    real,    dimension(:), allocatable :: vec_ra_tmp
    complex, dimension(:), allocatable :: vec_ca_tmp

    dim1_ra = size(vec_ra,1)

    k = isign*dim1_ra
    id = get_plan_id(shape(vec_ra), isign, NUM)

    select case (id)
    case (-1)

      ! test if the maximal number of plans is alredy reached.

      if (nplans(NUM) == MXPLAN) then
        write(*,*) 'four1D_real: MXPLAN too small! Increase it and recompile'
        stop
      end if

      allocate(vec_ra_tmp(dim1_ra), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: Allocation of  vec_ra_tmp  failed!'
        stop
      end if

      dim1_ca = size(vec_ca)
      allocate(vec_ca_tmp(dim1_ca), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: Allocation of  vec_ca_tmp  failed!'
        stop
      end if

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      id = nplans(NUM)

      select case (isign)
      case (FFT_FORWARD)
#if defined(real_precision_default)
        call sfftw_plan_dft_r2c_1d(plans(id,NUM), dim1_ra, vec_ra_tmp(1), &
           vec_ca_tmp(1), FFTW_MEASURE)
#else
        call dfftw_plan_dft_r2c_1d(plans(id,NUM), dim1_ra, vec_ra_tmp(1), &
           vec_ca_tmp(1), FFTW_MEASURE)
#endif
      case (FFT_INVERSE)
#if defined(real_precision_default)
        call sfftw_plan_dft_c2r_1d(plans(id,NUM), dim1_ra, vec_ca_tmp(1), &
           vec_ra_tmp(1), FFTW_MEASURE)
#else
        call dfftw_plan_dft_c2r_1d(plans(id,NUM), dim1_ra, vec_ca_tmp(1), &
           vec_ra_tmp(1), FFTW_MEASURE)
#endif
      end select

      deallocate(vec_ra_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: Dellocation of  vec_ra_tmp  failed!'
        stop
      end if

      deallocate(vec_ca_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: Dellocation of  vec_ca_tmp  failed!'
        stop
      end if

    end select

    ! Using the Guru execution of plans.

    select case (isign)
    case (FFT_FORWARD)
#if defined(real_precision_default)
      call sfftw_execute_dft_r2c(plans(id,NUM), vec_ra(1), vec_ca(1))
#else
      call dfftw_execute_dft_r2c(plans(id,NUM), vec_ra(1), vec_ca(1))
#endif
    case (FFT_INVERSE)
#if defined(real_precision_default)
      call sfftw_execute_dft_c2r(plans(id,NUM), vec_ca(1), vec_ra(1))
#else
      call dfftw_execute_dft_c2r(plans(id,NUM), vec_ca(1), vec_ra(1))
#endif
    end select

  end subroutine four1D_real

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ra_ca(arr_ra, arr_ca, isign)
    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign

    integer, parameter :: NUM=2, RANK=1 

    integer :: k
    integer :: dim1_ra, dim2_ra, dim1_ca, dim2_ca
    integer :: idist, odist, howmany, id, istat
    integer, dimension(RANK)             :: n_arr, inembed, onembed
    real,    dimension(:,:), allocatable :: arr_ra_tmp
    complex, dimension(:,:), allocatable :: arr_ca_tmp

    dim1_ra = size(arr_ra,1)
    dim2_ra = size(arr_ra,2)
    howmany = dim2_ra

    ! test if a plan that fits is already created.
    k = isign*dim1_ra
    id = get_plan_id(shape(arr_ra), isign, NUM)

    select case (id)
    case (-1)

      ! test if the maximal number of plans is alredy reached.

      if (nplans(NUM) == MXPLAN) then
        write(*,*) 'FOURCOL_RA_CA: MXPLAN too small! Increase it and recompile'
        stop
      end if

      allocate(arr_ra_tmp(dim1_ra, dim2_ra), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: Allocation of  arr_ra_tmp  failed!'
        stop
      end if

      dim1_ca = size(arr_ca,1)
      dim2_ca = size(arr_ca,2)
      allocate(arr_ca_tmp(dim1_ca, dim2_ca), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: Allocation of  arr_ca_tmp  failed!'
        stop
      end if

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      plan_properties(nplans(NUM),NUM,3) = dim2_ra
      id = nplans(NUM)

      select case (isign)
      case (FFT_FORWARD)
        n_arr(1) = dim1_ra
        howmany  = dim2_ra
        inembed(1) = size(arr_ra)
        onembed(1) = size(arr_ca)
        idist = dim1_ra
        odist = dim1_ca
#if defined(real_precision_default)
        call sfftw_plan_many_dft_r2c(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ra_tmp(1,1), inembed, 1, idist, &
           arr_ca_tmp(1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#else
        call dfftw_plan_many_dft_r2c(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ra_tmp(1,1), inembed, 1, idist, &
           arr_ca_tmp(1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#endif
      case (FFT_INVERSE)
        n_arr(1) = dim1_ra
        howmany  = dim2_ca
        inembed(1) = size(arr_ca)
        onembed(1) = size(arr_ra)
        idist = dim1_ca
        odist = dim1_ra
#if defined(real_precision_default)
        call sfftw_plan_many_dft_c2r(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ca_tmp(1,1), inembed, 1, idist, &
           arr_ra_tmp(1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#else
        call dfftw_plan_many_dft_c2r(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ca_tmp(1,1), inembed, 1, idist, &
           arr_ra_tmp(1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#endif
      end select

      deallocate(arr_ra_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: Dellocation of  arr_ra_tmp  failed!'
        stop
      end if

      deallocate(arr_ca_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: Dellocation of  arr_ca_tmp  failed!'
        stop
      end if

    end select

    ! Using the Guru execution of plans.

    select case (isign)
    case (FFT_FORWARD)
#if defined(real_precision_default)
      call sfftw_execute_dft_r2c(plans(id,NUM), arr_ra(1,1), arr_ca(1,1))
#else
      call dfftw_execute_dft_r2c(plans(id,NUM), arr_ra(1,1), arr_ca(1,1))
#endif
    case (FFT_INVERSE)
#if defined(real_precision_default)
      call sfftw_execute_dft_c2r(plans(id,NUM), arr_ca(1,1), arr_ra(1,1))
#else
      call dfftw_execute_dft_c2r(plans(id,NUM), arr_ca(1,1), arr_ra(1,1))
#endif
    end select

  end subroutine fourcol_ra_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_raa_caa(arr_raa, arr_caa, isign)
    real,    dimension(:,:,:), intent(inout) :: arr_raa
    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign

    integer, parameter :: NUM=3, RANK=1

    integer :: k
    integer :: dim1_raa, dim2_raa, dim3_raa, dim1_caa, dim2_caa, dim3_caa
    integer :: idist, odist, howmany, id, istat
    integer, dimension(RANK)               :: n_arr, inembed, onembed
    real,    dimension(:,:,:), allocatable :: arr_raa_tmp
    complex, dimension(:,:,:), allocatable :: arr_caa_tmp

    dim1_raa = size(arr_raa,1)
    dim2_raa = size(arr_raa,2)
    dim3_raa = size(arr_raa,3)
    howmany = dim2_raa*dim3_raa

    ! test if a plan that fits is already created.
    k = isign*dim1_raa
    id = get_plan_id(shape(arr_raa), isign, NUM)

    select case (id)
    case (-1)

      ! test if the maximal number of plans is alredy reached.

      if (nplans(NUM) == MXPLAN) then
        write(*,*) 'FOURCOL_RAA_CAA: MXPLAN too small! Increase it and recompile'
        stop
      end if

      allocate(arr_raa_tmp(dim1_raa, dim2_raa, dim3_raa), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: Allocation of  arr_raa_tmp  failed!'
        stop
      end if

      dim1_caa = size(arr_caa,1)
      dim2_caa = size(arr_caa,2)
      dim3_caa = size(arr_caa,3)
      allocate(arr_caa_tmp(dim1_caa, dim2_caa, dim3_caa), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: Allocation of  arr_caa_tmp  failed!'
        stop
      end if

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_raa
      plan_properties(nplans(NUM),NUM,3) = dim2_raa
      plan_properties(nplans(NUM),NUM,4) = dim3_raa
      id = nplans(NUM)

      select case (isign)
      case (FFT_FORWARD)
        n_arr(1) = dim1_raa
        howmany  = dim2_raa*dim3_raa
        inembed(1) = size(arr_raa)
        onembed(1) = size(arr_caa)
        idist = dim1_raa
        odist = dim1_caa
#if defined(real_precision_default)
        call sfftw_plan_many_dft_r2c(plans(id,NUM), RANK, n_arr, howmany, &
           arr_raa_tmp(1,1,1), inembed, 1, idist, &
           arr_caa_tmp(1,1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#else
        call dfftw_plan_many_dft_r2c(plans(id,NUM), RANK, n_arr, howmany, &
           arr_raa_tmp(1,1,1), inembed, 1, idist, &
           arr_caa_tmp(1,1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#endif
      case (FFT_INVERSE)
        n_arr(1) = dim1_raa
        howmany  = dim2_caa*dim3_caa
        inembed(1) = size(arr_caa)
        onembed(1) = size(arr_raa)
        idist = dim1_caa
        odist = dim1_raa
#if defined(real_precision_default)
        call sfftw_plan_many_dft_c2r(plans(id,NUM), RANK, n_arr, howmany, &
           arr_caa_tmp(1,1,1), inembed, 1, idist, &
           arr_raa_tmp(1,1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#else
        call dfftw_plan_many_dft_c2r(plans(id,NUM), RANK, n_arr, howmany, &
           arr_caa_tmp(1,1,1), inembed, 1, idist, &
           arr_raa_tmp(1,1,1), onembed, 1, odist, &
           FFTW_MEASURE)
#endif
      end select

      deallocate(arr_raa_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: Dellocation of  arr_raa_tmp  failed!'
        stop
      end if

      deallocate(arr_caa_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: Dellocation of  arr_caa_tmp  failed!'
        stop
      end if

    end select

    ! Using the Guru execution of plans.

    select case (isign)
    case (FFT_FORWARD)
#if defined(real_precision_default)
      call sfftw_execute_dft_r2c(plans(id,NUM), arr_raa(1,1,1), arr_caa(1,1,1))
#else
      call dfftw_execute_dft_r2c(plans(id,NUM), arr_raa(1,1,1), arr_caa(1,1,1))
#endif
    case (FFT_INVERSE)
#if defined(real_precision_default)
      call sfftw_execute_dft_c2r(plans(id,NUM), arr_caa(1,1,1), arr_raa(1,1,1))
#else
      call dfftw_execute_dft_c2r(plans(id,NUM), arr_caa(1,1,1), arr_raa(1,1,1))
#endif
    end select

  end subroutine fourcol_raa_caa


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four2D_real(arr_ra, arr_ca, isign)
    ! the position space array: a field of real numbers
    real,    dimension(:,:), intent(inout) :: arr_ra
    ! the fourier transformed array: a field of complex numbers
    complex, dimension(:,:), intent(inout) :: arr_ca
    ! the direction of the transformation:
    ! isign = 1 inverse fourier transformation: input arr_ca returns arr_ra
    ! isign = -1        fourier transformation: input arr_ra returns arr_ca
    integer,                 intent(in)    :: isign

    integer, parameter :: NUM=9

    integer :: k
    integer :: dim1_ra, dim2_ra, dim1_ca, dim2_ca, id, istat
    real,    dimension(:,:), allocatable :: arr_ra_tmp
    complex, dimension(:,:), allocatable :: arr_ca_tmp

    dim1_ra = size(arr_ra,1)
    dim2_ra = size(arr_ra,2)

    ! test if a plan that fits is already created.
    k = isign*dim1_ra
    id = get_plan_id(shape(arr_ra), isign, NUM)

    select case (id)
    case (-1)
      ! test if the maximal number of plans is alredy reached.
      if (nplans(NUM) == MXPLAN) then
        write(*,*) 'four2D_real: MXPLAN too small! Increase it and recompile'
        stop
      end if

      allocate(arr_ra_tmp(dim1_ra, dim2_ra), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: Allocation of  arr_ra_tmp  failed!'
        stop
      end if

      ! allocate a temporary array just as big as the k space array arr_ca
      ! given as an argument
      dim1_ca = size(arr_ca,1)
      dim2_ca = size(arr_ca,2)
      allocate(arr_ca_tmp(dim1_ca, dim2_ca), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: Allocation of  arr_ca_tmp  failed!'
        stop
      end if

      ! save the plan that will now be created for later reuse:
      ! remember the size of the real space array and
      ! the trafo direction
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      plan_properties(nplans(NUM),NUM,3) = dim2_ra
      id = nplans(NUM)

      select case (isign)
      case (FFT_FORWARD)
        ! create a plan to do a (forward) fft

        ! allocate a temporary array just as big as the real space
        ! array.  further down, this array is passed to the routine
        ! which creates the fft plan. Doing its tests the plan
        ! creation routine will overwrite values in both arrays
        ! passed as arguments.

#if defined(real_precision_default)
        call sfftw_plan_dft_r2c_2d(plans(id,NUM), dim1_ra, dim2_ra, &
           arr_ra_tmp(1,1), arr_ca_tmp(1,1), FFTW_MEASURE)
#else
        call dfftw_plan_dft_r2c_2d(plans(id,NUM), dim1_ra, dim2_ra, &
           arr_ra_tmp(1,1), arr_ca_tmp(1,1), FFTW_MEASURE)
#endif
      case (FFT_INVERSE)
        ! create a plan to do an inverse fft
#if defined(real_precision_default)
        call sfftw_plan_dft_c2r_2d(plans(id,NUM), dim1_ra, dim2_ra, &
           arr_ca_tmp(1,1), arr_ra_tmp(1,1), FFTW_MEASURE)
#else
        call dfftw_plan_dft_c2r_2d(plans(id,NUM), dim1_ra, dim2_ra, &
           arr_ca_tmp(1,1), arr_ra_tmp(1,1), FFTW_MEASURE)
#endif
      end select

      deallocate(arr_ra_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: Dellocation of  arr_ra_tmp  failed!'
        stop
      end if

      deallocate(arr_ca_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: Dellocation of  arr_ca_tmp  failed!'
        stop
      end if
    end select

    ! Using the basic execution of plans.
    select case (isign)
    case (FFT_FORWARD)
#if defined(real_precision_default)
      call sfftw_execute_dft_r2c(plans(id,NUM), arr_ra(1,1), arr_ca(1,1))
#else
      call dfftw_execute_dft_r2c(plans(id,NUM), arr_ra(1,1), arr_ca(1,1))
#endif
    case (FFT_INVERSE)
#if defined(real_precision_default)
      call sfftw_execute_dft_c2r(plans(id,NUM), arr_ca(1,1), arr_ra(1,1))
#else
      call dfftw_execute_dft_c2r(plans(id,NUM), arr_ca(1,1), arr_ra(1,1))
#endif
    end select

  end subroutine four2D_real


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ca(arr_ca, isign)
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign

    integer, parameter :: NUM=6, RANK=1 

    integer :: k
    integer :: dim1_ca, dim2_ca, howmany, id, istat, n
    integer, dimension(RANK)             :: n_arr, nembed
    complex, dimension(:,:), allocatable :: arr_ca_tmp

    dim1_ca = size(arr_ca,1)
    dim2_ca = size(arr_ca,2)
    howmany = size(arr_ca,2)

    ! test if a plan that fits is already created.
    k = isign*dim1_ca
    id = get_plan_id(shape(arr_ca), isign, NUM)

    select case (id)
    case (-1)

      ! test if the maximal number of plans is alredy reached.

      if (nplans(NUM) == MXPLAN) then
        write(*,*) 'FOURCOL_CA: MXPLAN too small! Increase it and recompile'
        stop
      end if

      allocate(arr_ca_tmp(dim1_ca, dim2_ca), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CA: Allocation of  arr_ca_tmp  failed!'
        stop
      end if

      nembed(1) = size(arr_ca)
      n_arr(1)  = dim1_ca
      n         = n_arr(1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ca
      plan_properties(nplans(NUM),NUM,3) = dim2_ca
      id = nplans(NUM)

      select case (isign)
      case (FFT_FORWARD)
#if defined(real_precision_default)
        call sfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           FFTW_FORWARD, FFTW_MEASURE)
#else
        call dfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           FFTW_FORWARD, FFTW_MEASURE)
#endif
      case (FFT_INVERSE)
#if defined(real_precision_default)
        call sfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           FFTW_BACKWARD, FFTW_MEASURE)
#else
        call dfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           arr_ca_tmp(1,1), nembed, 1, n, &
           FFTW_BACKWARD, FFTW_MEASURE)
#endif
      end select

      deallocate(arr_ca_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CA: Dellocation of  arr_ca_tmp  failed!'
        stop
      end if

    end select

#if defined(real_precision_default)
    call sfftw_execute_dft(plans(id,NUM), arr_ca(1,1), arr_ca(1,1))
#else
    call dfftw_execute_dft(plans(id,NUM), arr_ca(1,1), arr_ca(1,1))
#endif

  end subroutine fourcol_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_caa(arr_caa, isign)
    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign

    integer, parameter :: NUM=7, RANK=1

    integer :: k
    integer :: dim1_caa, dim2_caa, dim3_caa, howmany, id, istat, n
    integer, dimension(RANK)               :: n_arr, nembed
    complex, dimension(:,:,:), allocatable :: arr_caa_tmp

    dim1_caa = size(arr_caa,1)
    dim2_caa = size(arr_caa,2)
    dim3_caa = size(arr_caa,3)
    howmany = dim2_caa*dim3_caa

    ! test if a plan that fits is already created.
    k = isign*dim1_caa
    id = get_plan_id(shape(arr_caa), isign, NUM)

    select case (id)
    case (-1)

      ! test if the maximal number of plans is alredy reached.

      if (nplans(NUM) == MXPLAN) then
        write(*,*) 'FOURCOL_CAA: MXPLAN too small! Increase it and recompile'
        stop
      end if

      allocate(arr_caa_tmp(dim1_caa, dim2_caa, dim3_caa), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CAA: Allocation of  arr_caa_tmp  failed!'
        stop
      end if

      nembed(1) = size(arr_caa)
      n_arr(1)  = dim1_caa
      n         = n_arr(1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_caa
      plan_properties(nplans(NUM),NUM,3) = dim2_caa
      plan_properties(nplans(NUM),NUM,4) = dim3_caa
      id = nplans(NUM)

      select case (isign)
      case (FFT_FORWARD)
#if defined(real_precision_default)
        call sfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           FFTW_FORWARD, FFTW_MEASURE)
#else
        call dfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           FFTW_FORWARD, FFTW_MEASURE)
#endif
      case (FFT_INVERSE)
#if defined(real_precision_default)
        call sfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           FFTW_BACKWARD, FFTW_MEASURE)
#else
        call dfftw_plan_many_dft(plans(id,NUM), RANK, n_arr, howmany, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           arr_caa_tmp(1,1,1), nembed, 1, n, &
           FFTW_BACKWARD, FFTW_MEASURE)
#endif
      end select

      deallocate(arr_caa_tmp, stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CAA: Dellocation of  arr_caa_tmp  failed!'
        stop
      end if

    end select


#if defined(real_precision_default)
    call sfftw_execute_dft(plans(id,NUM), arr_caa(1,1,1), arr_caa(1,1,1))
#else
    call dfftw_execute_dft(plans(id,NUM), arr_caa(1,1,1), arr_caa(1,1,1))
#endif

  end subroutine fourcol_caa


#elif defined(FFT_ESSL)
!=============================================================================
! ESSL
!=============================================================================

  logical, parameter, public :: WORKING_FFT_LIBRARY=.true.

  ! string to copy error list entry

  CHARACTER (len=8), save :: S2015

  EXTERNAL :: ENOTRM

  ! auxilary arrays

  real, dimension(15) :: aux1
  real, dimension(1)  :: aux2


  ! work arrays for the ESSL routine

  type pointer_ra
    real, dimension(:), POINTER :: poi_ra
  end type pointer_ra

  type(pointer_ra), dimension(:,:), allocatable, save :: aux1d_poi1, aux1d_poi2
  type(pointer_ra), dimension(:,:), allocatable, save :: aux2d_poi1, aux2d_poi2

  integer, dimension(MXPLAN,N_FFT_ROUTINES,4), save :: plan_properties = 1

  ! number of plans saved

  integer, dimension(N_FFT_ROUTINES), save :: nplans=0

  interface fourrow_real
    module procedure fourrow_ra_ca
  end interface fourrow_real

  interface four1D
    module procedure four1D_ca
  end interface four1D

  interface fourrow
    module procedure fourrow_ca
  end interface fourrow

contains

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four2D_real(arr_ra, arr_ca, isign)

    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=9
    integer :: dim1_ra, dim2_ra, dim1_ca, dim2_ca, n1_ra, i, id, istat, &
       k, naux1, naux2

    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux2d_poi1(MXPLAN,1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: 1. Allocation of  aux2d_poi1  failed!'
        stop
      end if

      allocate(aux2d_poi2(MXPLAN,1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: 1. Allocation of  aux2d_poi2  failed!'
        stop
      end if

    end if

    dim1_ra = size(arr_ra,1)
    dim2_ra = size(arr_ra,2)
    dim1_ca = size(arr_ca,1)
    dim2_ca = size(arr_ca,2)

    ! the number of rows of data - that is, the length of the columns in
    ! array arr_ra involved in the computation.

    n1_ra = (dim1_ca-1)*2

    if (dim1_ra < n1_ra+2) then
      write(*,*) 'four2D_real: The stride of the array arr_ra is too small!'
      stop
    end if


    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_ra), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      plan_properties(nplans(NUM),NUM,3) = dim2_ra
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux1 = 1
      naux2 = size(aux2)

      select case (isign)
      case(FFT_FORWARD)
        call drcft2(1, arr_ra(1,1), dim1_ra,  &
           arr_ca(1,1), dim1_ca, &
           n1_ra, dim2_ra, -isign, 1.0, aux1, naux1, aux2, naux2)
      case(FFT_INVERSE)
        call dcrft2(1, arr_ca(1,1), dim1_ca,  &
           arr_ra(1,1), dim1_ra, &
           n1_ra, dim2_ra, -isign, 1.0, aux1, naux1, aux2, naux2)
      end select

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arrays.

      allocate(aux2d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: 2. Allocation of  aux2d_poi1  failed!'
        stop
      endif

      allocate(aux2d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four2D_real: 2. Allocation of  aux2d_poi2  failed!'
        stop
      endif

      select case (isign)
      case(FFT_FORWARD)
        call drcft2(1, arr_ra(1,1), dim1_ra, &
           arr_ca(1,1), dim1_ca, &
           n1_ra, dim2_ra, -isign, 1.0, &
           aux2d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux2d_poi2(id,NUM)%poi_ra(1), naux2)
      case(FFTW_INVERSE)
        call dcrft2(1, arr_ca(1,1), dim1_ca, &
           arr_ra(1,1), dim1_ra, &
           n1_ra, dim2_ra, -isign, 1.0, &
           aux2d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux2d_poi2(id,NUM)%poi_ra(1), naux2)
      end select

    end select

    select case (isign)
    case(FFT_FORWARD)
      call drcft2(0, arr_ra(1,1), dim1_ra, &
         arr_ca(1,1), dim1_ca, &
         n1_ra, dim2_ra, -isign, 1.0,  &
         aux2d_poi1(id,NUM)%poi_ra(1), size(aux2d_poi1(id,NUM)%poi_ra), &
         aux2d_poi2(id,NUM)%poi_ra(1), size(aux2d_poi2(id,NUM)%poi_ra))
    case(FFT_INVERSE)
      call dcrft2(0, arr_ca(1,1), dim1_ca, &
         arr_ra(1,1), dim1_ra, &
         n1_ra, dim2_ra, -isign, 1.0,  &
         aux2d_poi1(id,NUM)%poi_ra(1), size(aux2d_poi1(id,NUM)%poi_ra), &
         aux2d_poi2(id,NUM)%poi_ra(1), size(aux2d_poi2(id,NUM)%poi_ra))
    end select

  end subroutine four2D_real

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four1D_real(vec_ra, vec_ca, isign)

    real,    dimension(:), intent(inout) :: vec_ra
    complex, dimension(:), intent(inout) :: vec_ca
    integer,               intent(in)    :: isign
    integer, parameter :: NUM=1
    integer :: dim1_ra, dim1_ca, i, id, istat, k, naux1, naux2

    logical, save :: initflag=.true.


    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux1d_poi1(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: 1. Allocation of  aux1d_poi1  failed!'
        stop
      end if

      allocate(aux1d_poi2(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: 1. Allocation of  aux1d_poi2  failed!'
        stop
      end if

    end if

    dim1_ra = size(vec_ra)
    dim1_ca = size(vec_ca)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(vec_ra), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux2 = size(aux2)

      select case (isign)
      case(FFT_FORWARD)
        naux1 = 15
        call drcft(1, vec_ra(1), dim1_ra,  &
           vec_ca(1), dim1_ca, &
           dim1_ra, 1, -isign, 1.0, aux1, naux1, aux2, naux2)
      case(FFT_INVERSE)
        naux1 = 14
        call dcrft(1, vec_ca(1), dim1_ca,  &
           vec_ra(1), dim1_ra, &
           dim1_ra, 1, -isign, 1.0, aux1, naux1, aux2, naux2)
      end select

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arrays.

      allocate(aux1d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: 2. Allocation of  aux1d_poi1  failed!'
        stop
      endif

      allocate(aux1d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'four1D_real: 2. Allocation of  aux1d_poi2  failed!'
        stop
      endif

      select case (isign)
      case(FFT_FORWARD)
        call drcft(1, vec_ra(1), dim1_ra, &
           vec_ca(1), dim1_ca, &
           dim1_ra, 1, -isign, 1.0, &
           aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux1d_poi2(id,NUM)%poi_ra(1), naux2)
      case(FFT_INVERSE)
        call dcrft(1, vec_ca(1), dim1_ca, &
           vec_ra(1), dim1_ra, &
           dim1_ra, 1, -isign, 1.0, &
           aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux1d_poi2(id,NUM)%poi_ra(1), naux2)
      end select

    end select

    select case (isign)
    case(FFT_FORWARD)
      call drcft(0, vec_ra(1), dim1_ra, &
         vec_ca(1), dim1_ca, &
         dim1_ra, 1, -isign, 1.0,  &
         aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
         aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))
    case(FFT_INVERSE)
      call dcrft(0, vec_ca(1), dim1_ca, &
         vec_ra(1), dim1_ra, &
         dim1_ra, 1, -isign, 1.0,  &
         aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
         aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))
    end select

  end subroutine four1D_real

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ra_ca(arr_ra, arr_ca, isign)

    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=2
    integer :: dim1_ra, dim1_ca, howmany, i, id, istat, k, naux1, naux2

    logical, save :: initflag=.true.

    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux1d_poi1(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: 1. Allocation of  aux1d_poi1  failed!'
        stop
      end if

      allocate(aux1d_poi2(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: 1. Allocation of  aux1d_poi2  failed!'
        stop
      end if

    end if

    dim1_ra = size(arr_ra,1)
    howmany = size(arr_ra,2)
    dim1_ca = size(arr_ca,1)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(vec_ra), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      plan_properties(nplans(NUM),NUM,3) = size(arr_ra,2)
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux2 = size(aux2)

      select case (isign)
      case(FFT_FORWARD)
        naux1 = 15
        call drcft(1, arr_ra(1,1), dim1_ra,  &
           arr_ca(1,1), dim1_ca, &
           dim1_ra, howmany, -isign, 1.0, aux1, naux1, aux2, naux2)
      case(FFT_INVERSE)
        naux1 = 14
        call dcrft(1, arr_ca(1,1), dim1_ca,  &
           arr_ra(1,1), dim1_ra, &
           dim1_ra, howmany, -isign, 1.0, aux1, naux1, aux2, naux2)
      end select

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arrays.

      allocate(aux1d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: 2. Allocation of  aux1d_poi1  failed!'
        stop
      endif

      allocate(aux1d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RA_CA: 2. Allocation of  aux1d_poi2  failed!'
        stop
      endif

      select case (isign)
      case(FFT_FORWARD)
        call drcft(1, arr_ra(1,1), dim1_ra, &
           arr_ca(1,1), dim1_ca, &
           dim1_ra, howmany, -isign, 1.0, &
           aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux1d_poi2(id,NUM)%poi_ra(1), naux2)
      case(FFT_INVERSE)
        call dcrft(1, arr_ca(1,1), dim1_ca, &
           arr_ra(1,1), dim1_ra, &
           dim1_ra, howmany, -isign, 1.0, &
           aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux1d_poi2(id,NUM)%poi_ra(1), naux2)
      end select

    end select

    select case (isign)
    case(FFT_FORWARD)
      call drcft(0, arr_ra(1,1), dim1_ra, &
         arr_ca(1,1), dim1_ca, &
         dim1_ra, howmany, -isign, 1.0,  &
         aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
         aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))
    case(FFT_INVERSE)
      call dcrft(0, arr_ca(1,1), dim1_ca, &
         arr_ra(1,1), dim1_ra, &
         dim1_ra, howmany, -isign, 1.0,  &
         aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
         aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))
    end select

  end subroutine fourcol_ra_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_raa_caa(arr_raa, arr_caa, isign)

    real,    dimension(:,:,:), intent(inout) :: arr_raa
    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign
    integer, parameter :: NUM=3
    integer :: dim1_raa, dim1_caa, howmany, i, id, istat, k, naux1, naux2

    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux1d_poi1(MXPLAN,NUM), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: 1. Allocation of  aux1d_poi1  failed!'
        stop
      end if

      allocate(aux1d_poi2(MXPLAN,NUM), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: 1. Allocation of  aux1d_poi2  failed!'
        stop
      end if

    end if

    dim1_raa = size(arr_raa,1)
    howmany  = size(arr_raa,2)*size(arr_raa,3)
    dim1_caa = size(arr_caa,1)

    ! test if a plan that fits is already created.

    id = get_plan_id(shape(arr_raa), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_raa
      plan_properties(nplans(NUM),NUM,3) = dim2_raa
      plan_properties(nplans(NUM),NUM,4) = dim3_raa
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux2 = size(aux2)

      select case (isign)
      case(FFT_FORWARD)
        naux1 = 15
        call drcft(1, arr_raa(1,1,1), dim1_raa,  &
           arr_caa(1,1,1), dim1_caa, &
           dim1_raa, howmany, -isign, 1.0, aux1, naux1, aux2, naux2)
      case(FFT_INVERSE)
        naux1 = 14
        call dcrft(1, arr_caa(1,1,1), dim1_caa,  &
           arr_raa(1,1,1), dim1_raa, &
           dim1_raa, howmany, -isign, 1.0, aux1, naux1, aux2, naux2)
      end select

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arrays.

      allocate(aux1d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: 2. Allocation of  aux1d_poi1  failed!'
        stop
      endif

      allocate(aux1d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_RAA_CAA: 2. Allocation of  aux1d_poi2  failed!'
        stop
      endif

      select case (isign)
      case(FFT_FORWARD)
        call drcft(1, arr_raa(1,1,1), dim1_raa, &
           arr_caa(1,1,1), dim1_caa, &
           dim1_raa, howmany, -isign, 1.0, &
           aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux1d_poi2(id,NUM)%poi_ra(1), naux2)
      case(FFT_INVERSE)
        call dcrft(1, arr_caa(1,1,1), dim1_caa, &
           arr_raa(1,1,1), dim1_raa, &
           dim1_raa, howmany, -isign, 1.0, &
           aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
           aux1d_poi2(id,NUM)%poi_ra(1), naux2)
      end select

    end select

    select case (isign)
    case(FFT_FORWARD)
      call drcft(0, arr_raa(1,1,1), dim1_raa, &
         arr_caa(1,1,1), dim1_caa, &
         dim1_raa, howmany, -isign, 1.0,  &
         aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
         aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))
    case(FFT_INVERSE)
      call dcrft(0, arr_caa(1,1,1), dim1_caa, &
         arr_raa(1,1,1), dim1_raa, &
         dim1_raa, howmany, -isign, 1.0,  &
         aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
         aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))
    end select

  end subroutine fourcol_raa_caa

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourrow_ra_ca(arr_ra, arr_ca, isign)

    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign

    select case (isign)
    case(FFT_FORWARD)
      write(*,*) 'FOURROW_RA_CA: Real-to-Complex FFT of rows'// &
         ' not implementable in ESSL v4.2'
    case(FFT_INVERSE)
      write(*,*) 'FOURROW_RA_CA: Complex-to-Real FFT of rows'// &
         ' not implementable in ESSL v4.2'
    end select

    stop

  end subroutine fourrow_ra_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four1D_ca(vec_ca, isign)

    complex, dimension(:), intent(inout) :: vec_ca
    integer,               intent(in)    :: isign
    integer, parameter :: NUM=5
    integer :: dim1_ca, i, id, istat, k, n, naux1, naux2

    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux1d_poi1(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOUR1D_CA: 1. Allocation of  aux1d_poi1  failed!'
        stop
      end if

      allocate(aux1d_poi2(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOUR1D_CA: 1. Allocation of  aux1d_poi2  failed!'
        stop
      end if

    end if

    dim1_ca = size(vec_ca)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(vec_ca), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ca
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux1 = 8
      naux2 = size(aux2)

      call dcft(1, vec_ca(1), 1, dim1_ca,  vec_ca(1), 1, dim1_ca, &
         dim1_ca, 1, -isign, 1.0, aux1, naux1, aux2, naux2)

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arrays.

      allocate(aux1d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOUR1D_CA: 2. Allocation of  aux1d_poi1  failed!'
        stop
      endif

      allocate(aux1d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOUR1D_CA: 2. Allocation of  aux1d_poi2  failed!'
        stop
      endif

      call dcft(1, vec_ca(1), 1, dim1_ca, vec_ca(1), 1, dim1_ca, &
         dim1_ca, 1, -isign, 1.0, &
         aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
         aux1d_poi2(id,NUM)%poi_ra(1), naux2)

    end select

    call dcft(0, vec_ca(1), 1, dim1_ca, vec_ca(1), 1, dim1_ca, &
       dim1_ca, 1, -isign, 1.0,  &
       aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
       aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))

  end subroutine four1D_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ca(arr_ca, isign)

    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=6
    integer :: dim1_ca, howmany, i, id, istat, k, naux1, naux2

    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux1d_poi1(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CA: 1. Allocation of  aux1d_poi1  failed!'
        stop
      end if

      allocate(aux1d_poi2(MXPLAN,8), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CA: 1. Allocation of  aux1d_poi2  failed!'
        stop
      end if

    end if

    dim1_ca = size(arr_ca,1)
    howmany = size(arr_ca,2)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_ca), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ca
      plan_properties(nplans(NUM),NUM,3) = size(arr_ca,2)
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux1 = 8
      naux2 = size(aux2)

      call dcft(1, arr_ca(1,1), 1, dim1_ca,  arr_ca(1,1), 1, dim1_ca, &
         dim1_ca, howmany, -isign, 1.0, aux1, naux1, aux2, naux2)

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arrays.

      allocate(aux1d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CA: 2. Allocation of  aux1d_poi1  failed!'
        stop
      endif

      allocate(aux1d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CA: 2. Allocation of  aux1d_poi2  failed!'
        stop
      endif

      call dcft(1, arr_ca(1,1), 1, dim1_ca, arr_ca(1,1), 1, dim1_ca, &
         dim1_ca, howmany, -isign, 1.0, &
         aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
         aux1d_poi2(id,NUM)%poi_ra(1), naux2)

    end select

    call dcft(0, arr_ca(1,1), 1, dim1_ca, arr_ca(1,1), 1, dim1_ca, &
       dim1_ca, howmany, -isign, 1.0, &
       aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
       aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))

  end subroutine fourcol_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_caa(arr_caa, isign)

    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign
    integer, parameter :: NUM=7
    integer :: dim1_caa, howmany, i, id, istat, k, naux1, naux2

    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux1d_poi1(MXPLAN,NUM), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CAA: 1. Allocation of  aux1d_poi1  failed!'
        stop
      end if

      allocate(aux1d_poi2(MXPLAN,NUM), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CAA: 1. Allocation of  aux1d_poi2  failed!'
        stop
      end if

    end if

    dim1_caa = size(arr_caa,1)
    howmany  = size(arr_caa,2)*size(arr_caa,3)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_caa), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_caa
      plan_properties(nplans(NUM),NUM,3) = dim2_caa
      plan_properties(nplans(NUM),NUM,4) = dim3_caa
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux1 = 8
      naux2 = size(aux2)

      call dcft(1, arr_caa(1,1,1), 1, dim1_caa, arr_caa(1,1,1), 1, dim1_caa, &
         dim1_caa, howmany, -isign, 1.0, aux1, naux1, aux2, naux2)

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arrays.

      allocate(aux1d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CAA: 2. Allocation of  aux1d_poi1  failed!'
        stop
      endif

      allocate(aux1d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURCOL_CAA: 2. Allocation of  aux1d_poi2  failed!'
        stop
      endif

      call dcft(1, arr_caa(1,1,1), 1, dim1_caa, arr_caa(1,1,1), 1, dim1_caa, &
         dim1_caa, howmany, -isign, 1.0, &
         aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
         aux1d_poi2(id,NUM)%poi_ra(1), naux2)

    end select

    call dcft(0, arr_caa(1,1,1), 1, dim1_caa, arr_caa(1,1,1), 1, dim1_caa, &
       dim1_caa, howmany, -isign, 1.0, &
       aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
       aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))

  end subroutine fourcol_caa

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourrow_ca(arr_ca, isign)

    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=8
    integer :: dim1_ca, dim2_ca, i, id, istat, k, naux1, naux2

    if (initflag) then
      initflag = .false.
      call EINFO(0)
      call ERRSAV(2015,S2015)

      allocate(aux1d_poi1(MXPLAN,NUM), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURROW_CA: 1. Allocation of  aux1d_poi1  failed!'
        stop
      end if

      allocate(aux1d_poi2(MXPLAN,NUM), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURROW_CA: 1. Allocation of  aux1d_poi2  failed!'
        stop
      end if

    end if

    dim1_ca = size(arr_ca,1)
    dim2_ca = size(arr_ca,2)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_ca), isign, NUM)

    select case (id)
    case (-1)

      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ca
      plan_properties(nplans(NUM),NUM,3) = dim2_ca
      id = nplans(NUM)

      call ERRSET(2015,0,-1,1,ENOTRM,0)

      naux1 = 8
      naux2 = size(aux2)

      call dcft(1, arr_ca(1,1), dim2_ca, 1,  arr_ca(1,1), dim2_ca, 1, &
         dim1_ca, dim2_ca, -isign, 1.0, aux1, naux1, aux2, naux2)

      call ERRSTR(2015,S2015)

      ! dynamic allocation of the work arr_caays.

      allocate(aux1d_poi1(id,NUM)%poi_ra(naux1), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURROW_CA: 2. Allocation of  aux1d_poi1  failed!'
        stop
      endif

      allocate(aux1d_poi2(id,NUM)%poi_ra(naux2), stat=istat)
      if (istat /= 0) then
        write(*,*) 'FOURROW_CA: 2. Allocation of  aux1d_poi2  failed!'
        stop
      endif

      call dcft(1, arr_ca(1,1), dim1_ca, 1, arr_ca(1,1), dim1_ca, 1, &
         dim1_ca, dim2_ca, -isign, 1.0, &
         aux1d_poi1(id,NUM)%poi_ra(1), naux1, &
         aux1d_poi2(id,NUM)%poi_ra(1), naux2)

    end select

    call dcft(0, arr_ca(1,1), dim1_ca, 1, arr_ca(1,1), dim1_ca, 1, &
       dim1_ca, dim2_ca, -isign, 1.0,  &
       aux1d_poi1(id,NUM)%poi_ra(1), size(aux1d_poi1(id,NUM)%poi_ra), &
       aux1d_poi2(id,NUM)%poi_ra(1), size(aux1d_poi2(id,NUM)%poi_ra))

  end subroutine fourrow_ca


#elif defined(FFT_MKL)
!=============================================================================
! MKL
!=============================================================================

  public :: fourrow_real
  public :: four1D, fourrow

  logical, parameter, public :: WORKING_FFT_LIBRARY=.true.

  type pointer_r
    type(DFTI_DESCRIPTOR), pointer :: desc_handle
  end type pointer_r

  !> descriptor handles for FFT transforms
  type(pointer_r), dimension(MXPLAN,N_FFT_ROUTINES), save :: handles
  integer, dimension(MXPLAN,N_FFT_ROUTINES,4), save :: plan_properties = 1
  !> number of descriptor handles saved
  integer, dimension(N_FFT_ROUTINES), save :: nplans=0

  interface fourrow_real
    module procedure fourrow_ra_ca
  end interface fourrow_real

  interface four1D
    module procedure four1D_ca
  end interface four1D

  interface fourrow
    module procedure fourrow_ca
  end interface fourrow

contains

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four2D_real(arr_ra, arr_ca, isign)

    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=9
    logical :: init_flag
    integer :: dim1_ra, dim2_ra, dim1_ca, id, i, k, status

    dim1_ra = size(arr_ra,1)
    dim2_ra = size(arr_ra,2)
    dim1_ca = size(arr_ca,1)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(vec_ra), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      plan_properties(nplans(NUM),NUM,3) = dim2_ra
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call four2D_mkl_ra_ca(arr_ra(1,1), arr_ca(1,1), dim1_ra, dim2_ra, &
       dim1_ca, isign, init_flag, id, NUM)

  end subroutine four2D_real

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four1D_real(vec_ra, vec_ca, isign)

    integer,               intent(in)    :: isign
    real,    dimension(:), intent(inout) :: vec_ra
    complex, dimension(:), intent(inout) :: vec_ca
    integer, parameter :: NUM=1
    logical :: init_flag
    integer :: dim1_ra, dim1_ca, id, i, k, status

    dim1_ra = size(vec_ra,1)
    dim1_ca = size(vec_ca,1)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(vec_ra), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourcol_mkl_ra_ca(vec_ra(1), vec_ca(1), dim1_ra, dim1_ca, &
       1, isign, init_flag, id, NUM)

  end subroutine four1D_real

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ra_ca(arr_ra, arr_ca, isign)

    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=2
    logical :: init_flag
    integer :: dim1_ra, dim1_ca, howmany, id, i, k, status

    dim1_ra = size(arr_ra,1)
    howmany = size(arr_ra,2)
    dim1_ca = size(arr_ca,1)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_ra), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ra
      plan_properties(nplans(NUM),NUM,3) = size(arr_ra,2)
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourcol_mkl_ra_ca(arr_ra(1,1), arr_ca(1,1), dim1_ra, dim1_ca, &
       howmany, isign, init_flag, id, NUM)

  end subroutine fourcol_ra_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_raa_caa(arr_raa, arr_caa, isign)

    real,    dimension(:,:,:), intent(inout) :: arr_raa
    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign
    integer, parameter :: NUM=3
    logical :: init_flag
    integer :: dim1_raa, dim1_caa, howmany, id, i, k, status

    dim1_raa = size(arr_raa,1)
    howmany  = size(arr_raa,2)*size(arr_raa,3)
    dim1_caa = size(arr_caa,1)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_raa), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_raa
      plan_properties(nplans(NUM),NUM,3) = dim2_raa
      plan_properties(nplans(NUM),NUM,4) = dim3_raa
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourcol_mkl_ra_ca(arr_raa(1,1,1), arr_caa(1,1,1), dim1_raa, dim1_caa, &
       howmany, isign, init_flag, id, NUM)

  end subroutine fourcol_raa_caa

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourrow_ra_ca(arr_ra, arr_ca, isign)

    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=4
    logical :: init_flag
    integer :: dim2_ra, dim2_ca, howmany, id, i, k, status

    howmany = size(arr_ra,1)
    dim2_ra = size(arr_ra,2)
    dim2_ca = size(arr_ca,2)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_ra), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = size(arr_ra,1)
      plan_properties(nplans(NUM),NUM,3) = size(arr_ra,2)
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourrow_mkl_ra_ca(arr_ra(1,1), arr_ca(1,1), howmany, &
       dim2_ra, dim2_ca, isign, init_flag, id, NUM)

  end subroutine fourrow_ra_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four1D_ca(vec_ca, isign)

    integer,               intent(in)    :: isign
    complex, dimension(:), intent(inout) :: vec_ca
    integer, parameter :: NUM=5
    logical :: init_flag
    integer :: dim1_ca, id, i, k, status

    dim1_ca = size(vec_ca,1)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(vec_ca), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim1_ca
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourcol_mkl_ca(vec_ca(1), dim1_ca, 1, isign, init_flag, id, NUM)

  end subroutine four1D_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ca(arr_ca, isign)

    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=6
    logical :: init_flag
    integer :: dim_ca, howmany, id, i, k, status

    dim_ca  = size(arr_ca,1)
    howmany = size(arr_ca,2)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_ca), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dim_ca
      plan_properties(nplans(NUM),NUM,3) = size(arr_ca,2)
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourcol_mkl_ca(arr_ca(1,1), dim_ca, howmany, isign, init_flag, id, NUM)

  end subroutine fourcol_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_caa(arr_caa, isign)

    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign
    integer, parameter :: NUM=7
    logical :: init_flag
    integer :: dim_caa, howmany, id, i, k, status

    dim_caa = size(arr_caa,1)
    howmany = size(arr_caa,2)*size(arr_caa,3)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_caa), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = dimm_caa
      plan_properties(nplans(NUM),NUM,3) = size(arr_caa,2)
      plan_properties(nplans(NUM),NUM,4) = size(arr_caa,3)
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourcol_mkl_ca(arr_caa(1,1,1), dim_caa, howmany, isign, init_flag, id, NUM)

  end subroutine fourcol_caa

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourrow_ca(arr_ca, isign)

    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign
    integer, parameter :: NUM=8
    logical :: init_flag
    integer :: dim2_ca, howmany, id, i, k, status

    howmany = size(arr_ca,1)
    dim2_ca = size(arr_ca,2)

    ! test if a plan that fits is already created.
    id = get_plan_id(shape(arr_ca), isign, NUM)

    select case (id)
    case (-1)
      init_flag = .true.
      nplans(NUM) = nplans(NUM)+1
      plan_properties(nplans(NUM),NUM,1) = isign
      plan_properties(nplans(NUM),NUM,2) = size(arr_ca,1)
      plan_properties(nplans(NUM),NUM,3) = size(arr_ca,2)
      id = nplans(NUM)

    case default
      init_flag = .false.
    end select

    call fourrow_mkl_ca(arr_ca(1,1), howmany, dim2_ca, isign, init_flag, id, NUM)

  end subroutine fourrow_ca
  
  !--------------------------------------------------------------------
  !> COMMENT: This subroutine is necessary to prevent the Lahey/Fujitsu
  !>          compiler from making a copy of array arr_ra when passing
  !>          arguments.
  !--------------------------------------------------------------------
  subroutine fourcol_mkl_ra_ca(arr_ra, arr_ca, dim1_ra, dim1_ca, howmany, &
     isign, init_flag, id, num)

    real,    dimension(*), intent(inout) :: arr_ra
    complex, dimension(*), intent(inout) :: arr_ca
    integer,               intent(in)    :: dim1_ra
    integer,               intent(in)    :: dim1_ca
    integer,               intent(in)    :: howmany
    integer,               intent(in)    :: isign
    logical,               intent(in)    :: init_flag
    integer,               intent(in)    :: id
    integer,               intent(in)    :: num
    integer :: i, status

    if (init_flag) then

      status = DftiCreateDescriptor(handles(id,num)%desc_handle, &
         DFTI_double, DFTI_real, 1, dim1_ra)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 0:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 1:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_NUMBER_OF_TRANSFORMS, howmany)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 2:', DftiErrorMessage(status)

      select case (isign)
      case (FFT_FORWARD)
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_INPUT_DISTANCE, dim1_ra)
        if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 3:', DftiErrorMessage(status)
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_OUTPUT_DISTANCE, dim1_ca)
        if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 4:', DftiErrorMessage(status)
      case(FFT_INVERSE)
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_INPUT_DISTANCE, dim1_ca)
        if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 5:', DftiErrorMessage(status)
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_OUTPUT_DISTANCE, dim1_ra)
        if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 6:', DftiErrorMessage(status)
      end select

      status = DftiCommitDescriptor(handles(id,num)%desc_handle)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_RA_CA 7:', DftiErrorMessage(status)

    end if

    select case (isign)
    case (FFT_FORWARD)
      status = DftiComputeForward(handles(id,num)%desc_handle, arr_ra, arr_ca)
    case (FFT_INVERSE)
      status = DftiComputeBackward(handles(id,num)%desc_handle, arr_ca, arr_ra)
    end select

  end subroutine fourcol_mkl_ra_ca

  !--------------------------------------------------------------------
  !> COMMENT: This subroutine is necessary to prevent the Lahey/Fujitsu
  !>          compiler from making a copy of array arr_ra when passing
  !>          arguments.
  !--------------------------------------------------------------------
  subroutine fourrow_mkl_ra_ca(arr_ra, arr_ca, howmany, dim2_ra, dim2_ca, &
     isign, init_flag, id, num)

    real,    dimension(*), intent(inout) :: arr_ra
    complex, dimension(*), intent(inout) :: arr_ca
    integer,               intent(in)    :: howmany
    integer,               intent(in)    :: dim2_ra
    integer,               intent(in)    :: dim2_ca
    integer,               intent(in)    :: isign
    logical,               intent(in)    :: init_flag
    integer,               intent(in)    :: id
    integer,               intent(in)    :: num
    integer :: i, status
    integer, dimension(2) :: stride

    if (init_flag) then

      stride(1) = 0
      stride(2) = howmany

      status = DftiCreateDescriptor(handles(id,num)%desc_handle, &
         DFTI_double, DFTI_real, 1, dim2_ra)
      if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 0:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 1:', DftiErrorMessage(status)

      select case (isign)
      case (FFT_FORWARD)
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_NUMBER_OF_TRANSFORMS, dim2_ra)
        if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 2:', DftiErrorMessage(status)
      case (FFT_INVERSE)
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_NUMBER_OF_TRANSFORMS, dim2_ca)
        if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 3:', DftiErrorMessage(status)
      end select

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_INPUT_DISTANCE, 1)
      if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 4:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_INPUT_STRIDES, stride)
      if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 5:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_OUTPUT_DISTANCE, 1)
      if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 6:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_OUTPUT_STRIDES, stride)
      if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 7:', DftiErrorMessage(status)

      status = DftiCommitDescriptor(handles(id,num)%desc_handle)
      if (status /= 0) write(*,*) 'FOURROW_MKL_RA_CA 8:', DftiErrorMessage(status)

    end if

    select case (isign)
    case (FFT_FORWARD)
      status = DftiComputeForward(handles(id,num)%desc_handle, arr_ra, arr_ca)
    case (FFT_INVERSE)
      status = DftiComputeBackward(handles(id,num)%desc_handle, arr_ca, arr_ra)
    end select

  end subroutine fourrow_mkl_ra_ca

  !--------------------------------------------------------------------
  !> COMMENT: This subroutine is necessary to prevent the Lahey/Fujitsu
  !>          compiler from making a copy of array arr_ra when passing
  !>          arguments.
  !--------------------------------------------------------------------
  subroutine four2D_mkl_ra_ca(arr_ra, arr_ca, dim1_ra, dim2_ra, dim1_ca, &
     isign, init_flag, id, num)

    real,    dimension(*), intent(inout) :: arr_ra
    complex, dimension(*), intent(inout) :: arr_ca
    integer,               intent(in)    :: dim1_ra
    integer,               intent(in)    :: dim2_ra
    integer,               intent(in)    :: dim1_ca
    integer,               intent(in)    :: isign
    logical,               intent(in)    :: init_flag
    integer,               intent(in)    :: id
    integer,               intent(in)    :: num
    integer :: i, status
    integer, dimension(2) :: length
    integer, dimension(3) :: strides_in, strides_out

    if (init_flag) then

      length(1) = dim1_ra
      length(2) = dim2_ra
      status = DftiCreateDescriptor(handles(id,num)%desc_handle, &
         DFTI_double, DFTI_real, 2, length)
      if (status /= 0) write(*,*) 'FOUR2D_MKL_RA_CA 0:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_CONJUGATE_EVEN_STORAGE, DFTI_complex_complex)
      if (status /= 0) write(*,*) 'FOUR2D_MKL_RA_CA 1:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      if (status /= 0) write(*,*) 'FOUR2D_MKL_RA_CA 2:', DftiErrorMessage(status)

      select case (isign)
      case (FFT_FORWARD)
        strides_out(1) = 0
        strides_out(2) = 1
        strides_out(3) = dim1_ca
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_OUTPUT_STRIDES, strides_out)
        if (status /= 0) write(*,*) 'FOUR2D_MKL_RA_CA 3:', DftiErrorMessage(status)

      case (FFT_INVERSE)
        strides_in(1) = 0
        strides_in(2) = 1
        strides_in(3) = dim1_ca
        status = DftiSetValue(handles(id,num)%desc_handle, &
           DFTI_INPUT_STRIDES, strides_in)
        if (status /= 0) write(*,*) 'FOUR2D_MKL_RA_CA 4:', DftiErrorMessage(status)
      end select

      status = DftiCommitDescriptor(handles(id,num)%desc_handle)
      if (status /= 0) write(*,*) 'FOUR2D_MKL_RA_CA 5:', DftiErrorMessage(status)

    end if

    select case (isign)
    case (FFT_FORWARD)
      status = DftiComputeForward(handles(id,num)%desc_handle, arr_ra, arr_ca)
    case (FFT_INVERSE)
      status = DftiComputeBackward(handles(id,num)%desc_handle, arr_ca, arr_ra)
    end select

  end subroutine four2D_mkl_ra_ca

  !--------------------------------------------------------------------
  !> COMMENT: This subroutine is necessary to prevent the Lahey/Fujitsu
  !>          compiler from making a copy of array arr_ca when passing
  !>          arguments.
  !--------------------------------------------------------------------
  subroutine fourcol_mkl_ca(arr_ca, dim1_ca, howmany, isign, init_flag, id, num)

    complex, dimension(*), intent(inout) :: arr_ca
    integer,               intent(in)    :: dim1_ca
    integer,               intent(in)    :: howmany
    integer,               intent(in)    :: isign
    logical,               intent(in)    :: init_flag
    integer,               intent(in)    :: id
    integer,               intent(in)    :: num
    integer :: i, status

    if (init_flag) then

      status = DftiCreateDescriptor(handles(id,num)%desc_handle, &
         DFTI_double, DFTI_complex, 1, dim1_ca)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_CA 0:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_NUMBER_OF_TRANSFORMS, howmany)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_CA 1:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_INPUT_DISTANCE, dim1_ca)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_CA 2:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_OUTPUT_DISTANCE, dim1_ca)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_CA 3:', DftiErrorMessage(status)

      status = DftiCommitDescriptor(handles(id,num)%desc_handle)
      if (status /= 0) write(*,*) 'FOURCOL_MKL_CA 4:', DftiErrorMessage(status)

    end if

    select case (isign)
    case (FFT_FORWARD)
      status = DftiComputeForward(handles(id,num)%desc_handle, arr_ca)
    case (FFT_INVERSE)
      status = DftiComputeBackward(handles(id,num)%desc_handle, arr_ca)
    end select

  end subroutine fourcol_mkl_ca

  !--------------------------------------------------------------------
  !>This subroutine is necessary to prevent the Lahey/Fujitsu
  !> compiler from making a copy of array arr_ca when passing
  !> arguments.
  !--------------------------------------------------------------------
  subroutine fourrow_mkl_ca(arr_ca, howmany, dim2_ca, isign, init_flag, id, num)
    complex, dimension(*), intent(inout) :: arr_ca
    integer,               intent(in)    :: howmany
    integer,               intent(in)    :: dim2_ca
    integer,               intent(in)    :: isign
    logical,               intent(in)    :: init_flag
    integer,               intent(in)    :: id
    integer,               intent(in)    :: num

    integer :: i, status
    integer, dimension(2) :: stride

    if (init_flag) then

      stride(1) = 0
      stride(2) = howmany

      status = DftiCreateDescriptor(handles(id,num)%desc_handle, &
         DFTI_double, DFTI_complex, 1, dim2_ca)
      if (status /= 0) write(*,*) 'FOURROW_MKL_CA 0:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_NUMBER_OF_TRANSFORMS, dim2_ca)
      if (status /= 0) write(*,*) 'FOURROW_MKL_CA 1:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_INPUT_DISTANCE, 1)
      if (status /= 0) write(*,*) 'FOURROW_MKL_CA 2:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_INPUT_STRIDES, stride)
      if (status /= 0) write(*,*) 'FOURROW_MKL_CA 3:', DftiErrorMessage(status)

      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_OUTPUT_DISTANCE, 1)
      if (status /= 0) write(*,*) 'FOURROW_MKL_CA 4:', DftiErrorMessage(status)
      status = DftiSetValue(handles(id,num)%desc_handle, &
         DFTI_OUTPUT_STRIDES, stride)
      if (status /= 0) write(*,*) 'FOURROW_MKL_CA 5:', DftiErrorMessage(status)

      status = DftiCommitDescriptor(handles(id,num)%desc_handle)
      if (status /= 0) write(*,*) 'FOURROW_MKL_CA 6:', DftiErrorMessage(status)

    end if

    select case (isign)
    case (FFT_FORWARD)
      status = DftiComputeForward(handles(id,num)%desc_handle, arr_ca)
    case (FFT_INVERSE)
      status = DftiComputeBackward(handles(id,num)%desc_handle, arr_ca)
    end select

  end subroutine fourrow_mkl_ca


#else
!=============================================================================
! DUMMY FFT
!=============================================================================

  logical, parameter, public :: WORKING_FFT_LIBRARY=.false.
  
  character (len=4), parameter, public :: fftlib = 'none'  


  interface four1D
    module procedure four1D_ca
  end interface four1D

contains

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four1D_real(vec_ra, vec_ca, isign)
    real,    dimension(:), intent(inout) :: vec_ra
    complex, dimension(:), intent(inout) :: vec_ca
    integer,               intent(in)    :: isign

    ! keep the compiler quiet regarding unused dummy arguments
    if (WORKING_FFT_LIBRARY) then
      write (*,*) vec_ra, vec_ca, isign
    end if

  end subroutine four1D_real

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ra_ca(arr_ra, arr_ca, isign)
    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign

    ! keep the compiler quiet regarding unused dummy arguments
    if (WORKING_FFT_LIBRARY) then
      write (*,*) arr_ra, arr_ca, isign
    end if

  end subroutine fourcol_ra_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_raa_caa(arr_raa, arr_caa, isign)
    real,    dimension(:,:,:), intent(inout) :: arr_raa
    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign

    ! keep the compiler quiet regarding unused dummy arguments
    if (WORKING_FFT_LIBRARY) then
      write (*,*) arr_raa, arr_caa, isign
    end if


  end subroutine fourcol_raa_caa


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four2D_real(arr_ra, arr_ca, isign)
    real,    dimension(:,:), intent(inout) :: arr_ra
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign

    ! keep the compiler quiet regarding unused dummy arguments
    if (WORKING_FFT_LIBRARY) then
      write (*,*) arr_ra, arr_ca, isign
    end if


  end subroutine four2D_real

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine four1D_ca(vec_ca, isign)
    complex, dimension(:), intent(inout) :: vec_ca
    integer,               intent(in)    :: isign

    ! keep the compiler quiet regarding unused dummy arguments
    if (WORKING_FFT_LIBRARY) then
      write (*,*) vec_ca, isign
    end if

  end subroutine four1D_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_ca(arr_ca, isign)
    complex, dimension(:,:), intent(inout) :: arr_ca
    integer,                 intent(in)    :: isign

    ! keep the compiler quiet regarding unused dummy arguments
    if (WORKING_FFT_LIBRARY) then
      write (*,*) arr_ca, isign
    end if

  end subroutine fourcol_ca

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine fourcol_caa(arr_caa, isign)
    complex, dimension(:,:,:), intent(inout) :: arr_caa
    integer,                   intent(in)    :: isign

    ! keep the compiler quiet regarding unused dummy arguments
    if (WORKING_FFT_LIBRARY) then
      write (*,*) arr_caa, isign
    end if

  end subroutine fourcol_caa


#endif
end module fft
