!-----------------------------------------------------------------------------
!> Various switches to be set at compile time, which enable the code to be
!> run in unusual ways. These may be for debugging purposes, but can be for
!> essentially anything that will be changed infrequently.
!-----------------------------------------------------------------------------
module switches

  implicit none

  private

  !
  ! The following switches are for use in the normalise-general branch
  !

  !> When using imod_init, option to use only the zero radial mode.
  logical, parameter, public :: lmod_init_zerorad = .false.

  !> Split the read/write of fdisi between at least min_fdisi_io_subsets
  !> distinct open+read/write+close operations. Increase this value if
  !> problems are encountered with write/read of large restart files when
  !> running on many processors.  See subroutine set_comm_fdisi_io()
  integer, parameter, public :: min_fdisi_io_subsets = 1

  !> Do not split the read/write of fdisi via min_fdisi_io_subsets unless the
  !> number of processors is greater than min_procs_io_split.
  integer, parameter, public :: min_procs_io_split = 64



end module switches
