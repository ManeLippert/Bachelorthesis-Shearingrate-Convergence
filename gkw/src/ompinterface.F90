!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Provide an interface to the openmp routines
!> This module is probably unncessary, if !$ were used on the library routines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module ompinterface
#ifdef _OPENMP
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE omp_lib
#endif

implicit none

public :: ompreport_and_init, ompfinal_report, update_thread_statistics
public :: ompget_thread_num, ompget_max_threads, ompget_num_threads

private

integer, save :: nl_threads_min
integer, save :: nl_threads_max

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> 
  !---------------------------------------------------------------------------
  subroutine ompreport_and_init()
    use global, only : i_huge
    use mpiinterface, only: root_processor
#ifdef _OPENMP
    integer (kind=omp_sched_kind) :: schedule
    !> chunk size for 'guided' schedule
    integer :: modifier
#endif

    nl_threads_max = 0
    nl_threads_min = i_huge
#ifdef _OPENMP
    call omp_set_dynamic(.false.)
    ! only static scheduling makes sense for GKW, since it is memory bound.
    ! chunk size -1 means default.
    call omp_set_schedule(omp_sched_static, -1)
#endif

    if(root_processor) then
      ! set the env variable OMP_NUM_THREADS to control this
      write (*,*) "Number of threads requested when encountering parallel &
         & regions (on each process):", &
         & ompget_max_threads()
#ifdef _OPENMP
      write (*,*) "Dynamic adjustment of the number of threads: ", &
         & omp_get_dynamic()
      call omp_get_schedule(schedule, modifier)
      write (*,*) "Thread scheduling strategy: ", &
         & schedule, modifier 
      ! static 1
      ! dynamic 2
      ! guided 3
      ! auto 4
      write (*,*) "Nested thread parallelism: ", &
         & omp_get_nested()
      write (*,*) "Max. active levels of parallelism: ", &
         & omp_get_max_active_levels()
      ! set OMP_THREAD_LIMIT to control this:
      write (*,*) "Max. number of threads per process at the same time: ", &
         & omp_get_thread_limit()
#endif
      write (*,*)
    end if
  end subroutine ompreport_and_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> get info how many threads actually were generated for the current team
  !---------------------------------------------------------------------------
  subroutine update_thread_statistics()
    nl_threads_min = min(nl_threads_min, ompget_num_threads())
    nl_threads_max = max(nl_threads_max, ompget_num_threads())
  end subroutine update_thread_statistics
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> 
  !-----------------------------------------------------------------------------
  subroutine ompfinal_report()
    use global, only : int2char
    use general, only : gkw_warn
    use mpiinterface, only: root_processor

    if(root_processor .and. nl_threads_max /= 0) then
      if(nl_threads_min == nl_threads_max) then
        write (*,*) 'Throughout the run, '//int2char(nl_threads_max)// &
           & ' threads were for computation of nonlinear terms.'
      else
        ! will most likely never be the case
        call gkw_warn('Between '//int2char(nl_threads_min)//' and '// &
           & int2char(nl_threads_max)//' threads&
           & were used for computation of nonlinear terms.')
      end if
      write(*,*)
    end if

  end subroutine ompfinal_report

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> wrapper for omp_get_thread_num
  !-----------------------------------------------------------------------------
  function ompget_thread_num()
    integer :: ompget_thread_num

#ifdef _OPENMP
    ompget_thread_num = omp_get_thread_num()
#else
    ompget_thread_num = 0
#endif
  end function ompget_thread_num


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> wrapper for omp_get_max_threads
  !-----------------------------------------------------------------------------
  function ompget_max_threads()
    integer :: ompget_max_threads

#ifdef _OPENMP
    ompget_max_threads = omp_get_max_threads()
#else
    ompget_max_threads = 1
#endif
  end function ompget_max_threads

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>
  !---------------------------------------------------------------------------fg
  function ompget_num_threads()
    integer :: ompget_num_threads

#ifdef _OPENMP
    ompget_num_threads = omp_get_num_threads()
#else
    ompget_num_threads = 1
#endif
  end function ompget_num_threads

end module ompinterface
