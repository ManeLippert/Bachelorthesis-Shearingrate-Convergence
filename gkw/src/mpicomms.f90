!-----------------------------------------------------------------------------
!> This module keeps the additional MPI communicators
!-----------------------------------------------------------------------------
module mpicomms

  implicit none

  public

  !> a communicator incorporating all processors
  integer, save :: COMM_CART

  !> a communicator incorporating only the root process
  integer, save :: COMM_ROOT

  !> procs for the same points in all directions other than s
  integer, save :: COMM_S_NE
  !> procs for the same points in all directions other than x
  integer, save :: COMM_X_NE
  !> procs for the same points in all directions other than x and s
  integer, save :: COMM_S_NE_X_NE
  !> procs for the same points in all directions other than mu and s
  integer, save :: COMM_S_NE_MU_NE
  !> procs for the same points in all directions other than x and sp
  integer, save :: COMM_SP_NE_X_NE
  !> procs for the same points in all directions other than s and sp
  integer, save :: COMM_SP_NE_S_NE

  !> procs for the same points in all directions other than vpar
  integer, save :: COMM_VPAR_NE
  !> procs for the same points in all directions other than vpar and mu
  integer, save :: COMM_VPAR_NE_MU_NE
  !> procs for the same point in velocity space (aka COMM_SP_NE_S_NE_X_NE)
  integer, save :: COMM_VPAR_EQ_MU_EQ
  !> procs for the same points in all directions other than species and mu
  integer, save :: COMM_SP_NE_MU_NE
  !> all procs responsible for the same points in vpar
  integer, save :: COMM_VPAR_EQ
  !> all procs responsible for the same points in vpar and x
  integer, save :: COMM_VPAR_EQ_X_EQ

  !> procs for the same points in all directions other than mu
  integer, save :: COMM_MU_NE
  !> all procs responsible for the same points in mu
  integer, save :: COMM_MU_EQ

  !> procs for the same points in all directions, but different species
  integer, save :: COMM_SP_NE

  !> all procs responsible for the same species
  integer, save :: COMM_SP_EQ
  integer, save :: COMM_SP_EQ_nprocsets  !< number of `colors' in the comm split
  integer, save :: COMM_SP_EQ_procset    !< local `color' from comm split

  !> all procs responsible for the same part of the s-grid
  integer, save :: COMM_S_EQ
  integer, save :: COMM_S_EQ_nprocsets  !< number of `colors' in the comm split
  integer, save :: COMM_S_EQ_procset    !< local `color' from comm split

  !> all procs responsible for the same part of the x-grid
  integer, save :: COMM_X_EQ
  integer, save :: COMM_X_EQ_nprocsets  !< number of `colors' in the comm split
  integer, save :: COMM_X_EQ_procset    !< local `color' from comm split
  
  !> intersection of COMM_SP_EQ and COMM_S_EQ
  integer, save :: COMM_SP_EQ_S_EQ
  integer, save :: COMM_SP_EQ_S_EQ_nprocsets !< number of `colors' in the comm split
  integer, save :: COMM_SP_EQ_S_EQ_procset   !< local `color' from comm split

  !> intersection of COMM_S_EQ and COMM_X_EQ
  integer, save :: COMM_S_EQ_X_EQ
  
  !> intersection of COMM_SP_EQ and COMM_X_EQ
  integer, save :: COMM_SP_EQ_X_EQ
  
  !> used for toroidal mode work division in nonspectral field solve ONLY
  integer, save :: COMM_MOD
  
  !> a communicator for self
  integer, save :: COMM_ALL_EQ

  !> commumicator primarily for reading/writing restart files
  integer, save :: COMM_FDISI_IO
  integer, save :: COMM_FDISI_IO_nprocsets  !< number of splits from COMM_CART
  integer, save :: COMM_FDISI_IO_procset    !< a label for the local set

  !> a dummy communicator
  integer, parameter :: COMM_DUMMY = -93456

end module mpicomms
