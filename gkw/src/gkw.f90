!----------------------------------------------------------------------------
!! Doxygen frontpage
!> \mainpage
!! 
!! GKW is a gyro-kinetic simulation code for the study of turbulence 
!! in magnetised plasmas developed from the original linear code LINART.
!!
!! This is the Doxygen generated <b>reference dictionary</b> for GKW.
!! The quickest way to lookup a variable or routine is to use the <b>search
!! box above.</b>
!! 
!! The project source homepage, is at <b>http://bitbucket.org/gkw/gkw</b>
!! where you can:  
!!  - Find the full latex documentation and the source
!!  - Contribute to the issue tracker
!!  - Review and comment on changes to the source
!!  - Obtain new versions of the code
!!
!!<hr> 
!!   Copyright (C) 2003, 2004, 2005
!!     A.G. Peeters, D. Strintzi
!! 
!!   Copyright (C) 2007, 2008, 2009, 2010
!!     A.G. Peeters, Y. Camenen, F.J. Casson, W.A. Hornsby, A.P. Snodin,
!!     D. Strintzi, G. Szepesi
!!
!!  Copyright (C) 2012, 2013, 2014, 2015
!!    A.G. Peeters, Y. Camenen, F.J. Casson, W.A. Hornsby, A.P. Snodin,
!!    D. Strintzi, G. Szepesi, R. Buchholz, S. Grosshauser, P. Manas, 
!!    P. Migliano, M. Siccinio, T. Sung,  D. Zarzoso 
!!<hr> 
!! 
!! A detailed description of the code, how to build it, and how run it
!! can be found in the associated paper [1].
!! 
!! If you use GKW (or some results obtained from it) in any publication, we
!! politely request that you cite paper [1] below, in which the code is
!! comprehensively described. If you wish, you could also cite the original
!! LINART paper [2] and/or the first paper in which GKW was used [3].
!! 
!!  - [1] A.G. Peeters Y. Camenen, F.J. Casson, W.A. Hornsby, A.P. Snodin,
!!      D. Strintzi, and G. Szepesi,
!!      Computer Physics Communications, <b>180</b>, 2650 (2009)
!!      http://dx.doi.org/10.1016/j.cpc.2009.07.001
!! 
!!  - [2] A.G. Peeters, D. Strintzi, Phys. Plasmas, <b>11</b>, 3748 (2004)
!!      http://dx.doi.org/10.1063/1.1762876
!! 
!!  - [3] A.G. Peeters, C. Angioni, D. Strintzi,
!!      Phys. Rev. Lett. <b>98</b>, 265003 (2007) 
!!      http://dx.doi.org/10.1103/PhysRevLett.98.265003
!! 
!! GKW is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!! 
!! GKW is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License (http://www.gnu.org/licenses/gpl-3.0.txt)
!! for more details.
!< End Doxygen frontpage
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Main Program.
!> Top level for GKW
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
program gkw
!----------------------------------------------------------------------------
! Normalization: all length scales are normalized
! to the major radius. All velocities are normalized
! the the thermal velocity of the ions. 
! NOTE vth = sqrt(2 T / m )
!--------------------------------------------------------------------

  use mpiinterface,    only : mpiinit, mpiwtime, root_processor
  use mpiinterface,    only : mpibarrier, mpibcast, mpifinalize
  use ompinterface,    only : ompreport_and_init, ompfinal_report
  use control,         only : ntime, max_seconds,method, ntotstep
  use control,         only : stop_me, ndump_ts, itime, irun
  use control,         only : last_smallstep_time, last_largestep_time, time
  use control,         only : stop_filename
  use diagnostic,      only : diagnostic_naverage, diagnostic_initial_output, &
                            & diagnostic_final_output
  use restart,         only : write_restart_file
  use exp_integration, only : advance_large_step_explicit, init_explicit
  use init,            only : initialize, finalize, deallocate_runtime_arrays
  use imp_integration, only : imp_int
  use eiv_integration, only : eiv_solve
  use general,         only : gkw_exit
  use perform,         only : perfout, perf_measure
  use normalise, only : normalise_after_timestep
  use control, only : get_start_time, max_runtime_reached, get_max_runtime
  use exp_integration, only : set_persistent_mode
  use fields, only : calculate_fields, field_solve_nonspec_wrap
  use control, only : spectral_radius 
  use components, only : mode_persist
  use dist, only : fdisi
  use control, only : laverage_dist_over_time, itime_rst
  use diagnos_f, only : tavg_end

  implicit none

  integer :: i, idump=0
  !> time at the end of the run
  double precision :: t_begin, t_end
  double precision :: t_1, t_predict
  double precision :: t_begin_main, t_end_main
  double precision :: t_diag_b, t_diag, t_step_b, t_step
  real    :: t_1r
  logical :: exstop, exdump

  ! Exit GKW already before anywork is done, if a certain file is
  ! found.
  exstop = .false.
  inquire(file=stop_filename, exist = exstop)
  if(exstop) then
    write(*,*) 'External stop, '//stop_filename//' is present.'
    open (9, FILE = stop_filename)
    close (9, STATUS='delete')
    stop 1
  end if

  ! initialize mpi
  call mpiinit()
  call ompreport_and_init()
  
  ! set the start time of the run
  t_begin = get_start_time()

  ! call the initialization
  call initialize



  ! get the time at the start of the main loops
  t_begin_main = mpiwtime()
  
  ! check if there is anything to output at the beginning of the run
  call diagnostic_initial_output()


  select case(method) 
  case('EIV')
    call eiv_solve
    ! get the main loop end time
    t_end_main = mpiwtime()

  case default

    last_largestep_time = time
    ! the setting here does not matter with
    ! current diagnostics:
    last_smallstep_time = time

    if(root_processor) then
      write (*,*)
      write (*,*) "Enter main loop (time integration)..."
      write (*,*)
    end if
    
    
    ! Set timesteps since last restart file has been written to zero.
    itime_rst = 0

    ! loop over large time steps
    large_time_steps : do i = 1, ntime

      ! store loop value in control
      itime = i
      
      ! Update timesteps since last restart file has been written.
      ! Do this before diagnostics are called, since numbering of output
      ! files depends on itime_rst.
      itime_rst = itime_rst + 1

      t_step_b = mpiwtime() 

      select case(method) 
      case('EXP')
        ! do the time integration naverage times 
        call advance_large_step_explicit(i, .true.)
      case('IMP') 
        call imp_int
      case default 
        call gkw_exit('Unknown method of integration')
      end select

      call normalise_after_timestep()
      
      ! If there is a bi-normal mode that should be persistent in time 
      ! reset it to its constant value after normalisation.
      ! F.Rath: It would be cleaner to put this into the normalise module.
      ! However, this is not possible due to the modules dependencies.
      if(mode_persist) then
        call set_persistent_mode(fdisi)
    
        ! Calculate the fields with the correct persistent mode
        if (spectral_radius) then
          call calculate_fields(fdisi)
        else
          call field_solve_nonspec_wrap(fdisi,0.,.false.)
        endif
      end if
  
       
      t_step = mpiwtime() - t_step_b

      ! calculate and output diagnostics
      t_diag_b = mpiwtime()
      call diagnostic_naverage()
      t_diag = mpiwtime() - t_diag_b

      ! Warn if the diagnostics are slow.
      ! This can happen if there are filesystem problems
      if (t_diag/t_step > 0.1) then
        if (root_processor) then
          write(*,*) 'WARNING: Diagnostics are slow:'
          write(*,'(A,es13.5)') ' Diagnostics time (last call):     ', t_diag  
          write(*,'(A,es13.5)') ' Non-diagnostics time (last step): ', t_step
          write(*,*)
        end if
      end if

      if (root_processor) then
200     format('Time step : ',I7,' Normalised time : ',es13.5)
        write(*,200) ntotstep, time
        write(*,*)
        write(*,*)
      end if
      
      
      ! If temporal averag of distribution finished , write it to restart file
      if(laverage_dist_over_time .and. itime == tavg_end) then
        ! write restart file (named FTA)
        call write_restart_file(.true.,0,.true.)
      end if

      ! check for external stop criterium
      exstop = .false.
      if (root_processor) inquire(file=stop_filename,exist = exstop)
      call mpibcast(exstop,1) ! broadcast from root
      if (exstop) then
        if (root_processor) then
          write(*,*) 'External stop, '//stop_filename//' is present.'
          open (9, FILE = stop_filename)
          close (9, STATUS='delete')
        end if
        exit large_time_steps
      end if

      !Predict total runtime after first iteration
      if (i == 1) then
        if(root_processor) then

          !Time for one iteration in seconds
          t_1 = mpiwtime()-t_begin_main
          t_1r = t_1
          !Total predicted runtime in seconds
          t_predict=t_1*ntime

          if (t_predict > 60.) then
            write(*,*)
            write(*,*) 'Iteration 1 completed successfully.'
            write(*,*) 'Predicted runtime: ', nint(t_predict/60.), ' minutes'
            write(*,*)
          end if

          if (t_predict > max_seconds .and. max_seconds > 0.) then
            write(*,*) 'WARNING: Run likely to terminate from max_seconds input'
            write(*,*) 'Stop will occur after ', nint(get_max_runtime(t_1r)/60.), 'minutes' 
            write(*,*) 'Iterations expected: ', nint(get_max_runtime(t_1r)/t_1), 'of ',  &
               &       ntime,' requested'
            write(*,*)
          end if

          if(perf_measure) then
            call perfout(i)
          end if

        end if

        call mpibcast(t_1r,1)
      end if !time predictor

      if(max_runtime_reached(t_1r)) then
        if (root_processor) then
          write(*,*) 'max_sec: stop' 
        end if
        exit large_time_steps
      end if

      ! check if the code needs to stop (convergence, timestep, ...)
      if (stop_me) then 
        if (root_processor) then 
          write(*,*) 'Internal: stop' 
        endif
        exit large_time_steps 
      endif
      

      ! check for dump criterion
      ! If the code hit a stop condition above, the last dump file 
      ! is not overwritten, because FDS will be written on exit anyway
      idump = idump + 1
      if (root_processor) inquire(file='gkw.dump',exist = exdump)
      ! every process waits for root to check if that file exists
      call mpibcast(exdump,1)
      if ((exdump.or.(idump == ndump_ts)).and. i /= ntime) then
        idump=0
        if (root_processor.and.exdump) then
          write(*,*)'External dump '
          open (9, FILE = 'gkw.dump')
          close (9, STATUS='delete')              
        end if

        call write_restart_file(.false., irun, .false.)
        
        ! Set timesteps since last restart file has been written to 
        ! zero.
        itime_rst = 0
        
      end if

      last_largestep_time = time

    end do large_time_steps

    ! get the main loop end time
    t_end_main = mpiwtime()

    ! write restart file first just in case other diagnostics fail
    call write_restart_file(.true., irun, .false.)

    ! check if there is anything to output at the end of the run
    call diagnostic_final_output()

  end select

  ! deallocate arrays used only in the main time loop
  call deallocate_runtime_arrays

  ! determine the CPU time used 
  t_end = mpiwtime()

  if (root_processor) then
    ! Please do not remove the success message from the following line
    ! unless the test script is modified accordingly.
    write(*,*)'Run successfully completed, Run time :',t_end-get_start_time()
    write(*,*)'                    (main loop time):',t_end_main-t_begin_main
    write(*,*)
    call ompfinal_report()
  end if

  if (root_processor .and. perf_measure) call perfout()

  ! finalize in general
  call finalize()

  ! finalize mpi 
  call mpifinalize()

  ! exit status zero indicates successful run
  stop 0

end program gkw
