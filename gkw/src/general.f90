!*****************************************************************************
!> General purpose routines, having no dependencies other than things
!> in mpiinterface, should be in here.
!*****************************************************************************
module general

  use global, only : int2char

  implicit none

  private

  public :: general_init
  public :: gkw_abort, gkw_warn, gkw_exit
  public :: matout
  public :: gkw_is_nan
  public :: max_factor

  contains

!-----------------------------------------------------------------------------
!> initialize anything for general
!-----------------------------------------------------------------------------
subroutine general_init()

end subroutine general_init

!-----------------------------------------------------------------------------
!> Abort the code safely (mpiabort) with an optional error message.
!> This routine can be used to abort from a module higher in the hierarchy
!> than IO.
!-----------------------------------------------------------------------------
subroutine gkw_abort(abort_message)
  use mpiinterface, only : processor_number
  use io, only : finalize_and_abort
  
  character (len=*),intent(in), optional :: abort_message

  if (present(abort_message)) then
    write (*,'(A,A,A,I4)') ' GKW_ABORT -- ', abort_message,                  &
                         & '  ', processor_number
  else
    write (*,'(A,A,A,I4)') ' GKW_ABORT -- ', '(no message given!)',          &
                         & '  ', processor_number
  end if

  call finalize_and_abort

end subroutine gkw_abort

!-----------------------------------------------------------------------------
!> Clean abort that must only be called from ALL procs. 
!> This is only to be used when we are sure that all processors
!> will call the routine (check params routines) 
!> Otherwise, gkw_abort should be used as a safe option to avoid MPI deadlock
!-----------------------------------------------------------------------------
subroutine gkw_exit(abort_message)
  use mpiinterface, only : root_processor
  use io, only : finalize_and_exit
  
  character(len=*), intent(in), optional :: abort_message

  if (root_processor) then
    if (present(abort_message)) then
      write (*,'(A,A,A)') ' GKW_EXIT -- ', abort_message, '  '
    else
      write (*,'(A,A,A)') ' GKW_EXIT -- ', '(no message given!)', '  '
    end if
  end if

  call finalize_and_exit

end subroutine gkw_exit

!-----------------------------------------------------------------------------
!> Function that should return true if the passed real is a NaN 
!> Attempts to be robust to all compilers, arcitectures and optimisations
!> but probably not as good as ieee_is_nan in 2003 standard
!-----------------------------------------------------------------------------
function gkw_is_nan(real_to_test)

  use mpiinterface, only : root_processor

  logical            :: gkw_is_nan
  logical            :: isnan1, isnan2, isnan3, isnan4, isnan5, isnan6, isnan7
  real, intent(in)   :: real_to_test
  character(len=128) :: mychar

  mychar='000'
  write(mychar,*) real_to_test
  mychar = adjustl(mychar)

  isnan1 = (trim(mychar)==trim(adjustl('NaN'))) 
  isnan2 = (trim(mychar)==trim(adjustl('NAN'))) 
  isnan3 = (trim(mychar)==trim(adjustl('nan'))) 
  isnan4 = (real_to_test /= real_to_test)
  isnan5 = ((real_to_test > 0.0) .EQV. (real_to_test <= 0.0))
  isnan6 = bittest_is_nan(real_to_test)

  !http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2004-05/0882.html
  if (exponent(real_to_test) > maxexponent(real_to_test)) then 
    if(fraction(real_to_test) == 0.5) then
      isnan7 = .false.  ! 'infinity' -> might want another function for this 
    else
      isnan7 = .true. 
    endif 
  else 
    isnan7 = .false.
  endif 

  gkw_is_nan = isnan1 .or. isnan2 .or. isnan3 .or. isnan4 .or. isnan5
  gkw_is_nan = gkw_is_nan .or. isnan6 .or. isnan7

  if (gkw_is_nan .and. root_processor) then
      write(*,*) 'isnan1', isnan1
      write(*,*) 'isnan2', isnan2
      write(*,*) 'isnan3', isnan3
      write(*,*) 'isnan4', isnan4
      write(*,*) 'isnan5', isnan5
      write(*,*) 'isnan6', isnan6
      write(*,*) 'isnan7', isnan7
   end if

end function gkw_is_nan

!-----------------------------------------------------------------------------
!> bit test for NaN on selected real precision
!>
!> The function operates by transferring bit pattern from a real variable to 
!> an integer container. The value is exclusive ORed with the value being 
!> tested for  (the ieee NAN bit representation hardcoded in global.F90)
!> The integer result of the IEOR function is converted to a logical result 
!> by comparing it to zero.  Note, other NaNs may also exist which do NOT
!> match the bit pattern, but these will not be ieee conforming.
!>
!> If the argument is array valued, the function returns a conformable logical
!> array, suitable for use with the ANY function, or as a logical mask.
!>
!> This function (only) is adapted for GKW from infnan.f90 
!> which is available from http://www.lahey.com/code.htm:
!!
!! Copyright(c) 2003, Lahey Computer Systems, Inc.
!! Copies of this source code, or standalone compiled files 
!! derived from this source may not be sold without permission
!! from Lahey Computers Systems. All or part of this [function] may be 
!! freely incorporated into executable programs which are offered
!! for sale. Otherwise, distribution of all or part of this file is
!! permitted, provided this copyright notice and header are included.
!----------------------------------------------------------------------------
elemental function bittest_is_nan(real_to_test) result(res)
    
  use global, only : iNaN
  
  ! Position of the sign bit (Intel) bit numbering starts at zero
  integer, parameter :: PSB = bit_size(iNaN) - 1
  real, intent(in)   :: real_to_test
  logical            :: res
    
  res = ieor(ibclr(transfer(real_to_test,iNaN),PSB), iNaN) == 0

end function

!-----------------------------------------------------------------------------
!> Check if io_stat is nonzero on all processes.
!>
!> FIXME This function is probably obsolete? At least no place uses it?
!----------------------------------------------------------------------------
function check_iostat(io_stat,abort_message)
  use mpiinterface, only : mpiallreduce_max, mpiallreduce_min
  character (len=*),intent(in), optional :: abort_message
  integer, intent(in)                    :: io_stat
  integer :: io_stat_max,io_stat_min
  logical :: check_iostat

  if (.false.) write(*,*) abort_message

  call mpiallreduce_max(io_stat,io_stat_max,1)
  call mpiallreduce_min(io_stat,io_stat_min,1)

  if (io_stat_min == 0 .and. io_stat_max == 0) then
    check_iostat = .false.
  else
    write (*,*) 'iostat=',io_stat
    write (*,*) 'iostat_min=',io_stat_min, 'iostat_max=',io_stat_max
    check_iostat = .true.
  end if

end function check_iostat

!-----------------------------------------------------------------------------
!> This routine prints a warning message, if called by the MPI root processor
!-----------------------------------------------------------------------------
subroutine gkw_warn(warning_message)  
  use mpiinterface, only : root_processor
  use global, only : gkw_warn_any_proc
  character(len=*), intent(in) :: warning_message 

  if (root_processor) then
    call gkw_warn_any_proc(warning_message)
  end if

end subroutine gkw_warn

!-----------------------------------------------------------------------------
!> Simple output routine that is only used for testing the matix elements 
!------------------------------------------------------------------------------
subroutine matout(iih,jjh,mat_elem) 

  integer, intent(in) :: iih, jjh 
  complex, intent(in) :: mat_elem 
  
  integer, save :: file_unit = 43423
  
  if (file_unit == 0) then
    open(file_unit, file = 'effe.dat') 
  endif 
  
  write(file_unit,fmt='(I6,1X,I6,1X,2(1e13.5,1X))') iih, jjh, mat_elem
  
end subroutine matout 

!-----------------------------------------------------------------------------
!> Very simple way to find the maximum factor of an integer up to a given no.,
!> used for grid decomposition checks.  Not intended to hunt for large primes!
!-----------------------------------------------------------------------------
pure function max_factor(num,max_fac)

  integer, intent(in) :: num, max_fac
  integer :: max_factor, i

  max_factor = 0

  do i = 1, max_fac
    if (mod(num, i) == 0) max_factor = i
  end do 

end function max_factor

end module general
