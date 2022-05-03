!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> GKW Input Checking Program.  This calls the initialization routines, with
!> the optional CHECK_INPUT set to .true.  This will catch any problems with
!> namelist formatting or variable names, and many problems of incompatible 
!> settings that are checked in the [namelist]_check_params routines.
!>
!> The checker is compiled with "make checker".  The checker is run in the
!> same way as GKW, in a directory with input.dat.   Multiple input files can 
!> be checked using the wrapper script "gkw_check_input_new file[1,2] input*".  
!>
!> However, if an input file passes the input checker, there is no guarantee
!> that GKW will run correctly, as many checks take place later than the 
!> initialisation stage.  The only way to check if GKW will run is to run GKW!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program gkw_input_check

  use mpiinterface, only : mpiinit
  use init,         only : initialize

  implicit none

  ! initialize mpi
  call mpiinit()

  ! call the initialization only for checking
  call initialize(check_input=.true.)

end program gkw_input_check

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
