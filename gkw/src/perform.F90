!-----------------------------------------------------------------------------
!
!>  Module that measures the performance of individual pieces of code.
!!
!!  Usage  : Call perfon just before the piece of code with as
!!           argument a string that contains the name used as
!!           indentifier for the code (maximum 30 charachters)
!!           After the loop call perfoff
!!           One can call the routine several times with the
!!           same name. The time spent will then be aded. Also
!!           loops can be nested
!!           The final output is obtained (on the screen and
!!           in file perform.dat) by calling perfout.
!!
!!           Example:
!!           use perform, only : perfon, perfoff, perfout
!!           call perfon('fft',4) ! as part of timing set 4
!!           call fft(a,ar,1)
!!           call perfoff(4)
!!           call perfout
!!
!! The output is a table which has the following format:
!! name_of_the_piece   number_of_calls  total_time_spent  percentage
!! where percentage is always measured with respect to name of the first
!! call in the code to perfon
!!
!! A.G. Peeters April 2007
!<-----------------------------------------------------------------------------
module perform

  implicit none

  private

  ! public subroutines
  public :: perfon, perfoff, perfswitch, perfout, perfloop
  
  ! public variables / parameters 
  public :: perffields, perfdo, perf_measure

  !> This parameter determines the maximum pieces of code that can be tracked.
  integer, parameter :: nregion = 500

  !> The maximum length of the character string
  integer, parameter :: maxchar = 30

  !> Number of loop iterations to record individually
  integer, parameter :: number_loop = 5000

  !> array in which the names of the indentifiers is stored
  character (len = maxchar), dimension(nregion), save :: perf_name

  !> This array keeps track of the names of the nested indentifiers
  integer, dimension(nregion), save :: perf_stack

  !> This integer keeps track of how many calls of perf are nested
  integer, save :: perf_stack_count

  !> Array in which the begin time for every call to perfon is stored.
  double precision, dimension(nregion), save :: perf_time

  !> Array in which the iteration time for every call to perfloop is stored.
  double precision, dimension(0:number_loop), save :: perf_loop

  !> Real array in which the sum of time spent in a certain piece is stored.
  double precision, dimension(nregion), save :: perf_total_time

  !> logical that is true if the perfon routine is called for a certain
  !> identifier. Calling perfon twice with the same name without closing in
  !> between with a call to perfoff results in an error message.
  logical, dimension(nregion), save :: perf_on

  !> Array that contains the counter for the number of calls with
  !> a certain identifier.
  integer, dimension(nregion), save :: perf_icount

  !> Warning message are stored in this string 
  character (len = 80)  :: perf_warn = '' 
  
  !> The current total of different indentifiers
  integer, save :: actual_region

  !> Loop counter
  integer, save :: iloop = 0

  !> Logical that determines if the routine has been called before
  logical, save :: initialised = .false.
  
  !> Logical to control fields timings
  logical, save :: perffields = .false.

  !> Logical to switch on / off all timings (use to avoid first loop)
  logical, save :: perfdo=.true.

  !> Switch, which enables/disables performance timing
  !> preferred over preprocessor statemenents for neatness
  !> This must remain a parameter, so the compiler can remove the perf calls 
  !> (they have a measureable performance impact in some linear runs)
#if defined(perf)
  logical, parameter :: perf_measure = .true.
#else
  logical, parameter :: perf_measure = .false.
#endif
  
contains

!-----------------------------------------------------------------------------
!> This subroutine start the measurement of the time for indentifier name.
!> If name does not exist, it will be added to the list. The name will also be
!> added to the stack. Two calls with the same name without a closing call in
!> between will result in an error message.
!-----------------------------------------------------------------------------
subroutine perfon(name,iperform)

  use control,      only : iperform_set
  use mpiinterface, only : mpiwtime 
  use general,      only : gkw_abort 
  
  ! The name of the piece of code
  integer,   intent(in) :: iperform
  character, intent(in) :: name*(*)

  ! local variables
  integer :: i

  ! is this part of the code followed ?
  if (.not.(iperform_set < 0 .or. iperform == iperform_set)) return 
  
  if (.not. perfdo) return

  ! if not called before initialise the counters
  if (.not. initialised) then
    actual_region = 0
    perf_stack_count = 0
    initialised = .true.
    perf_total_time(:) = 0.
  end if

  ! First check if the name already exists
  do i = 1, actual_region

    if (name == perf_name(i)) then

      ! check if the loop is open already
      if (perf_on(i)) then
        write(*,*)'perfon called twice for ',name
        call gkw_abort('PERFORM: internal error') 
      else
        perf_on(i) = .true.
      end if
      ! add 1 to the counter for the number of calls
      perf_icount(i) = perf_icount(i) + 1
      ! add the name to the stack
      perf_stack_count = perf_stack_count + 1
      if (perf_stack_count > nregion) then
        perf_warn = 'PERFORM ERROR : stack overflow for perf'
        perf_stack_count = perf_stack_count - 1
      end if
      perf_stack(perf_stack_count) = i

      ! call the timing routine. NOTE in perfon it is the last statement
      perf_time(i) = mpiwtime()

      ! jump out of the routine
      return

    end if

  end do

  ! if this part of the code is reached the name did not yet exist.
  ! First add the name to the list
  actual_region = actual_region + 1
  if (actual_region > nregion) then
    call gkw_abort('Too many pieces of code are being followed or &
                   & string is too long')
  end if

  ! add the name to the stack
  perf_stack_count = perf_stack_count + 1
  if (perf_stack_count > nregion) then
    call gkw_abort('stack overflow for perf')
  end if
  perf_stack(perf_stack_count) = i

  ! Store the name, identify that it is open and set the counter
  perf_name(actual_region) = name
  perf_on(i) = .true.
  perf_icount(actual_region) = 1

  ! call the time routine (Note this is the last statement
  perf_time(i) = mpiwtime()

end subroutine perfon

!-----------------------------------------------------------------------------
!> This routine ends the measurement. 
!> It always closes the last of the open loops and stores the time spent.
!-----------------------------------------------------------------------------
subroutine perfoff(iperform)

  use control,      only : iperform_set
  use mpiinterface, only : mpiwtime
  use general,      only : gkw_abort
  
  integer, intent(in) :: iperform
  
  integer          :: n
  double precision :: dum

  ! Follow this piece of code ? 
  if (.not.(iperform_set < 0 .or. iperform == iperform_set)) return 
  
  if (.not. perfdo) return

  ! call the timer (Note in perfoff this is the first statment)
  dum = mpiwtime()

  ! Store the current number
  n = perf_stack(perf_stack_count)

  ! reduce the stack
  perf_stack_count = perf_stack_count - 1

  ! check that perfoff is not called too often
  if (perf_stack_count < 0) then
    call gkw_abort('Error: perfoff is called more often than perfon')
  end if

  ! Store the time spent
  perf_total_time(n) = perf_total_time(n) + (dum - perf_time(n))

  ! Set that the loop is no longer active
  perf_on(n) = .false.

end subroutine perfoff

!-----------------------------------------------------------------------------
!> End the last measurement and immediately start the next one.
!-----------------------------------------------------------------------------
subroutine perfswitch(name,iperform)
  ! The name of the piece of code
  integer,   intent(in) :: iperform
  character, intent(in) :: name*(*)
  
  call perfoff(iperform)
  call perfon(name, iperform)
  
end subroutine perfswitch

!----------------------------------------------------------------------------- 
!> This routine stores the time for each loop iteration, up to number_loop
!-----------------------------------------------------------------------------
subroutine perfloop(iperform)

  use control,      only : iperform_set 
  use mpiinterface, only : mpiwtime
  
  integer, intent(in) :: iperform
  
  double precision, save :: last = 0.0
  double precision       :: dum
  
  ! follow this piece of code 
  if (.not. (iperform_set < 0 .or. iperform == iperform_set)) return 
  
  if (iloop <= number_loop) then
    dum = mpiwtime()
    perf_loop(iloop) = dum-last
    last = dum 
    iloop = iloop + 1
  else 
    ! Since these numbers are used to get an idea of machine 
    ! state intermittency, this routine is expected to discard later loops
    ! perf_warn = 'PERFORM WARNING: later perfloop timings strorage discarded'
  end if

end subroutine perfloop


!-----------------------------------------------------------------------------
!> This routine writes to the screen the output of the timing routine.
!> A table is made which has the following format:
!>
!>   name_of_the_piece   number_of_calls  total_time_spent  percentage
!>
!> where percentage is always measured with respect to name of the first
!> call in the code to perfon.
!-----------------------------------------------------------------------------
subroutine perfout(itime)

  use io,           only : get_free_file_unit
  use control,      only : iperform_set
  use mpiinterface, only : mpiwtick
  use global, only : r_tiny
  
  integer :: i, file_unit
  double precision :: percentage, perf_time
  integer, optional, intent(in) :: itime
  character(len=32) :: file_sec, file_loop

  perf_time = 0.E0

  if (iperform_set == 0) return

  ! it may happen that no timing was done at the first timestep, if
  ! naverage is equal to 2. Then we should avoid a division by zero
  ! in the output below.
  if (abs(perf_total_time(1)) < r_tiny) return

  ! check if all loops are closed
  if (perf_stack_count /= 0) then
    write(*,*)'WARNING not all loops were closed with perfoff'
    write(*,*)'Results may be corrupted'
  end if
  
  file_sec='perform.dat'
  file_loop='perfloop.dat'
  
  if (present(itime)) then 
    file_sec = 'perform_first.dat'
    file_loop = 'perfloop_first.dat'
  end if
  
  ! File for easy read with matlab: importdata('perform.dat')
  !oldio: call open_file(file_unit,trim(file_sec))
  call get_free_file_unit(file_unit)
  open(UNIT=file_unit,FILE=file_sec, FORM='formatted', &
     & STATUS='replace', POSITION='rewind')
  
  if (perf_warn /= '')  write(*,300) perf_warn
  write(*,300)'=====================================================',&
         &    '========='
  write(*,100)'Number of code pieces tracked ',actual_region
  write(*,300)'=====================================================',&
         &    '========='
  write(*,300)'Name                     Number of calls   time spent',&
         &    ' percent'
  write(*,300)'=====================================================',&
         &    '========='
  do i = 1, actual_region
    percentage = 1.E2*perf_total_time(i) / perf_total_time(1)
    if (i > 1) perf_time = perf_time + perf_total_time(i)
    write(*,200) perf_name(i), perf_icount(i), perf_total_time(i),    &
         & percentage
    write(file_unit,200) perf_name(i), perf_icount(i), perf_total_time(i),    &
         & percentage
  end do
  perf_time=perf_total_time(1)-perf_time

  write(*,250) 'perflib / other               ','-', perf_time,         &
     &  1.E2*perf_time / perf_total_time(1)
  write(*,*)
  write(*,400)'   (timing error: ', mpiwtick(), 's)'
  write(*,300)'=====================================================',&
         &    '========='
  100 format(A31,I5)
  200 format(A30,I10,E13.5,F9.2)
  250 format(A30,A10,E13.5,F9.2)
  300 format(A53,A9)
  400 format(A18,F9.6,A2)
  
  close(file_unit)

  call get_free_file_unit(file_unit)
  !oldio: call open_file(file_unit,trim(file_loop))
  open(UNIT=file_unit,FILE=file_loop, FORM='formatted', &
     & STATUS='replace', POSITION='rewind')
  do i = 1, iloop - 1
    write(file_unit,*) i, perf_loop(i)
  end do
  close(file_unit)
  
end subroutine perfout

end module perform
