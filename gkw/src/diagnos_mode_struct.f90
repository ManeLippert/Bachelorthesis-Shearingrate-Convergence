!------------------------------------------------------------------------------
!> Outputs field perturbations and moments for each mode.
!>
!------------------------------------------------------------------------------
module diagnos_mode_struct

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output

  
  !> switch for the parallel output 
  logical, save, public :: lparallel_output = .true.

  !> a list of timestamps can be provided in the namelist, to get more than
  !> one exemplaire of the parallel mode structure without rerunning the code
  real, save, public :: parallel_output_timestamps(512)
  
  !> switch for the rotation of the parallel output such that Re[phi(s=0)] = 1
  logical, save, public :: lrotate_parallel

  integer, dimension(7), save :: mode_struct_lus = -1
  integer, dimension(7), save :: mode_struct_ids
  character(len=10), dimension(7), save :: mode_struct_names

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
      lparallel_output = .true.
      parallel_output_timestamps = -1
      lrotate_parallel = .true.
  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lparallel_output,1)
    call mpibcast(parallel_output_timestamps,512)
    call mpibcast(lrotate_parallel,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use diagnos_generic, only : DENSITY_MOMENT, PAR_TEMPERATURE_MOMENT
    use diagnos_generic, only : PERP_TEMPERATURE_MOMENT, CURRENT_MOMENT
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD

    use diagnos_generic, only : DENSITY_MOMENT_TOKEN
    use diagnos_generic, only : PAR_TEMPERATURE_MOMENT_TOKEN, PERP_TEMPERATURE_MOMENT_TOKEN
    use diagnos_generic, only : CURRENT_MOMENT_TOKEN
    use diagnos_generic, only : PHI_FIELD_TOKEN, APAR_FIELD_TOKEN, BPAR_FIELD_TOKEN

    use diagnos_generic, only : LOCAL_DATA
    logical, intent(inout) :: requirements(:,:)

    if(.not. lparallel_output) return
    
    requirements(PHI_FIELD,LOCAL_DATA) = .true.
    requirements(APAR_FIELD,LOCAL_DATA) = .true.
    requirements(BPAR_FIELD,LOCAL_DATA) = .true.

    ! initialize the list which holds the names for
    ! the three fields...
    mode_struct_names(1) = PHI_FIELD_TOKEN
    mode_struct_ids(1) = PHI_FIELD
    mode_struct_names(2) = APAR_FIELD_TOKEN
    mode_struct_ids(2) = APAR_FIELD
    mode_struct_names(3) = BPAR_FIELD_TOKEN
    mode_struct_ids(3) = BPAR_FIELD
    ! ...and any moments which shall be output.
    mode_struct_names(4) = DENSITY_MOMENT_TOKEN  ! density fluctuations
    mode_struct_ids(4) = DENSITY_MOMENT
    mode_struct_names(5) = PAR_TEMPERATURE_MOMENT_TOKEN  ! parallel temperature fluctuations
    mode_struct_ids(5) = PAR_TEMPERATURE_MOMENT
    mode_struct_names(6) = PERP_TEMPERATURE_MOMENT_TOKEN ! parallel temperature fluct.
    mode_struct_ids(6) = PERP_TEMPERATURE_MOMENT
    mode_struct_names(7) = CURRENT_MOMENT_TOKEN ! parallel velocity field fluct.
    mode_struct_ids(7) = CURRENT_MOMENT

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine open_lus()
    use mpiinterface, only : root_processor
    use control, only : io_legacy, nlphi, nlapar, nlbpar
    use grid, only : number_of_species, nmod, n_s_grid, n_x_grid
    use io, only : open_complex_lu, binary_fmt
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD


    integer :: iquantity
    logical, save :: lus_are_opened = .false.
    
    if(root_processor .and. .not. lus_are_opened) then

      if(lparallel_output) then
        if(.not.io_legacy) then
          do iquantity = 1, size(mode_struct_names)
            if(iquantity == PHI_FIELD .and. .not.nlphi) cycle
            if(iquantity == BPAR_FIELD .and. .not.nlbpar) cycle
            if(iquantity == APAR_FIELD .and. .not.nlapar) cycle

            if(iquantity <= 3) then
              call open_complex_lu(mode_struct_names(iquantity), &
                 & 'diagnostic/diagnos_mode_struct', &
                 & (/ nmod, n_x_grid, n_s_grid /), &
                 & binary_fmt, mode_struct_lus(iquantity))
            else
              call open_complex_lu(mode_struct_names(iquantity), &
                 & 'diagnostic/diagnos_mode_struct', &
                 & (/ number_of_species, nmod, n_s_grid, n_x_grid /), &
                 & binary_fmt, mode_struct_lus(iquantity))
            end if
          end do
        end if
      end if

    end if
    lus_are_opened = .true.

  end subroutine open_lus

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use diagnos_moments, only : allocate_mem_moments => allocate_mem
    use diagnos_fields, only : allocate_mem_fields => allocate_mem

    logical, parameter :: allocate_3d_buffers = .true.

    if(lparallel_output .or. any(parallel_output_timestamps >= 0)) then
      call allocate_mem_moments(allocate_3d_buffers)
      call allocate_mem_fields(allocate_3d_buffers)
    end if

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()

  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()

    call mode_struct_output_at_times()

  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine final_output(number)
    use control, only : io_legacy
    !> an argument to number the output, is used to output eigenmode structures
    integer, intent(in), optional :: number
    
    if (lparallel_output) then
      call open_lus()
      
      if (present(number)) then
        if(io_legacy) then
          call parallel_output(number)
        else
          !call parallel_output(number)
          call output_mode_structure(number)
        end if
      else
        if(io_legacy) then
          call parallel_output()
        else
          !call parallel_output()
          call output_mode_structure()
        end if
      end if
    end if
  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output()
    
    call mode_struct_output_at_times()

  end subroutine output

  !----------------------------------------------------------------------------
  !> This routine checks if the time of the simulation has overstepped the next
  !> timestamp given in the input.dat file. Then a file is produced by
  !> parallel_output() .
  !----------------------------------------------------------------------------
  subroutine mode_struct_output_at_times()
    use control, only : time, last_largestep_time, io_legacy
    use global, only : r_tiny
    integer :: i
    if (lparallel_output) then
      if(any(parallel_output_timestamps >= 0)) then
        call open_lus()
      end if

      do i = 1, 512
        if(parallel_output_timestamps(i) >= 0) then
          !check if one of the timestamps in the list was passed
          if (abs(time) <= r_tiny .and. &
             & abs(time - parallel_output_timestamps(i)) < r_tiny) then
            if(io_legacy) then
              call parallel_output(i)
            else
              call output_mode_structure(i)
              !call parallel_output(i)
            end if
          elseif(time >= parallel_output_timestamps(i) .and. &
             & last_largestep_time < parallel_output_timestamps(i)) then
            if(io_legacy) then
              call parallel_output(i)
            else
              call output_mode_structure(i)
              !call parallel_output(i)
            end if
          end if
        end if
      end do
    end if
  end subroutine mode_struct_output_at_times

  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine output_mode_structure(inumber)
    use grid,           only : nmod
    use grid,           only : number_of_species, n_x_grid, n_s_grid
    use control,        only : non_linear, nlphi, nlapar, nlbpar
    use io,             only : append_chunk, xy_fmt
    use io,             only : binary_fmt
    use mpiinterface,   only : root_processor
    use diagnos_fields, only : gather_field_3d_serial
    use diagnos_moments, only : gather_moment_4d_serial
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use normalise, only : get_cmplx_rotation_factor

    integer, intent(in), optional :: inumber
    
    integer :: iquantity
    complex, dimension(number_of_species, &
       & nmod, n_s_grid, n_x_grid) :: moment_buf_global
    complex, dimension(nmod, n_x_grid, n_s_grid) :: field_buf_global
    character(len=20) :: name
    complex :: rotate(nmod)
    logical, parameter :: get_real_space_field = .false.

    ! TODO allocate the global arrays only on the root process

    ! The mode structure (i.e. fields and moments) may be output more than once:
    !  * for the eigenvalue solver, there may be multiple modes found, that
    !    need to be stored
    !  * it is possible to specify a list of timestamps where the mode should
    !    be stored

    ! To keep the compiler quiet.
    if (present(inumber)) continue
    
    rotate = 1.0
    
    do iquantity = 1,size(mode_struct_names)

      ! we would like to avoid numbers in the dataset name
      !name = trim(mode_struct_names(iquantity))//trim(adjustl(int2char(inumber,3)))
      ! therefore:
      name = trim(mode_struct_names(iquantity))

      if(iquantity == PHI_FIELD .and. .not.nlphi) cycle
      if(iquantity == BPAR_FIELD .and. .not.nlbpar) cycle
      if(iquantity == APAR_FIELD .and. .not.nlapar) cycle

      if(iquantity == PHI_FIELD) then
        rotate = get_cmplx_rotation_factor()
      end if
      
      ! At the moment, the IO is serial. The root process gathers the
      ! data from the other processes and does the IO.
      if(iquantity <= 3) then
        if(.not. non_linear .and. lrotate_parallel) then
          call gather_field_3d_serial(mode_struct_ids(iquantity), &
             & field_buf_global, get_real_space_field, rotate)
        else
          call gather_field_3d_serial(mode_struct_ids(iquantity), &
             & field_buf_global, get_real_space_field)
          rotate = 1.0
        end if
        if (root_processor) then
          write(*,*) "output mode struct ", name
          call append_chunk(mode_struct_lus(iquantity), &
             & field_buf_global, xy_fmt, binary_fmt)
        end if
      else
        if(.not. non_linear .and. lrotate_parallel) then
          call gather_moment_4d_serial(mode_struct_ids(iquantity), &
             & moment_buf_global, get_real_space_field, rotate)
        else
          call gather_moment_4d_serial(mode_struct_ids(iquantity), &
             & moment_buf_global, get_real_space_field)
        end if
        if (root_processor) then
          write(*,*) "output mode struct ", name
          call append_chunk(mode_struct_lus(iquantity), &
             & moment_buf_global, xy_fmt, binary_fmt)
        end if
      end if
    end do ! iquantity

  end subroutine output_mode_structure


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>  This subroutine writes some quantities connected with the parallel 
  !>  mode structure. For every mode the following quantities are written
  !>  to file
  !>
  !>  sgr         the length along the field line 
  !>  phi         the potential 
  !>  apar        the parallel component of the vector potential 
  !>  dens        the perturbed density 
  !>  epar        the perturbed parallel energy 
  !>  eperp       the perturbed perpendicular energy 
  !>  parallel_velocity       the perturbed parallel flow velocity
  !>  bpar        the perturbed parallel (compressional) magnetic field
  !>
  !>  Note only sgr is real. The other quantities are complex and fill 2
  !>  columns in the output dataset.
  !>
  !----------------------------------------------------------------------------
  subroutine parallel_output(inumber)
    use grid,           only : nmod, nx, ns, nmu, nvpar, lsp, proc_subset
    use grid,           only : number_of_species, n_x_grid, n_s_grid
    use control,        only : non_linear, flux_tube, io_legacy
    use control,        only : nlphi, nlapar, nlbpar
    use dist,           only : fdisi, phi, apar, bpar
    use geom,           only : sgr, bn
    use mode,           only : ixzero
    use velocitygrid,   only : intmu, intvp, mugr, vpgr
    use io,             only : open_real_lu, close_lu, append_chunk, ascii_fmt
    use global,         only : r_tiny
    use matdat,         only : get_f_from_g
    use mpiinterface,   only : mpiallreduce_maxloc, mpibcast, mpiallreduce_sum
    use mpiinterface,   only : mpibarrier, gather_array, root_processor
    use mpiinterface,   only : number_of_processors, processor_number
    use mpicomms,       only : COMM_S_EQ_X_EQ, COMM_X_NE, COMM_S_NE
    use global,        only : int2char
    use global,         only : dotdat, gkw_a_equal_b_accuracy

    integer, optional, intent(in) :: inumber

    ! integers for the loop over all grid points 
    integer :: imod, ix, i, j, k, is, ihelp, il, i_s
    integer, dimension(2) :: ihelp2

    real, dimension(2) :: phimax_local,phimax_global
    ! for reduction:
    complex, dimension(4) :: arr, rarr
    complex :: phid, apad, bpad, epar, eperp, dens, parallel_velocity
    complex :: fdis, rotate
    complex, dimension(1:nx,1:ns,8) :: local_write
    complex, dimension(1:n_x_grid,1:n_s_grid,8) :: global_write

    integer :: lun

    ! A number can added to the logical unit name as
    !  * for the eigenvalue solver, there may be multiple modes found, that
    !    need to be stored
    !  * it is possible to specify a list of timestamps where the mode should
    !    be stored
    if (root_processor) then
      if (present(inumber)) then
        call open_real_lu(dotdat( &
           &'parallel'//trim(adjustl(int2char(inumber,5))), io_legacy), &
           & 'diagnostic/diagnos_mode_struct', (/ 15 /), &
           & ascii_fmt, lun)
      else
        call open_real_lu(dotdat('parallel',io_legacy), &
           & 'diagnostic/diagnos_mode_struct', (/ 15 /), &
           & ascii_fmt, lun)
      end if
    end if

    do i_s = 1, number_of_species
      do imod = 1, nmod

        if (.not.non_linear .and. lrotate_parallel) then
          ! Find the maximum of the potential and with (mpi) maxloc 
          ! Find a (single) processor on which that point exists
          if (flux_tube) then
            ihelp=maxloc(abs(phi(imod,ixzero,1:ns)),1)
            ! beware using maxloc with array indices not starting at 1
            phimax_local(1)=abs(phi(imod,ixzero,ihelp))
            rotate=phi(imod,ixzero,ihelp)
          else
            ihelp2 = maxloc(abs(phi(imod,1:nx,1:ns)))
            !beware using maxloc with array indices not starting at 1
            rotate = phi(imod,ihelp2(1),ihelp2(2))
            phimax_local(1) = abs(rotate)
          end if

          ! use real(processor_number) for maxloc interface
          ! (could make a more Fortran friendly wrapper here)
          phimax_local(2)=real(1.*processor_number)
          call mpiallreduce_maxloc(phimax_local,phimax_global,1)
          ! convert back to integer for comparison
          ihelp=int(phimax_global(2))
          !Broadcast the point with maximum potential to all.
          call mpibcast(rotate,2,PROC=ihelp)
          !Avoid dividing by zero
          if (abs(rotate) < r_tiny) rotate = (1.,0.)
        else
          rotate=(1.,0.)
        endif

        do ix = 1, nx
          do i = 1, ns

            epar  = 0.
            eperp = 0.
            dens  = 0.
            parallel_velocity = 0.

            is = lsp(i_s)
            ! When using parallel species, only perform the velocity space sums
            ! when the outer loop reaches the local species.
            if (proc_subset(0,0,0,0,i_s)) then
              do k = 1, nvpar
                do j = 1, nmu

                  ! The distribution function 
                  fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)
                  if (gkw_a_equal_b_accuracy(intvp(i,j,k,is), 0.0)) fdis = 0.
                  
                  ! perturbed density
                  dens = dens + bn(ix,i)*intvp(i,j,k,is)*intmu(j)*fdis
                  
                  ! perturbed parallel velocity
                  parallel_velocity= parallel_velocity + bn(ix,i)*intvp(i,j,k,is)*intmu(j)*vpgr(i,j,k,is)*   &
                     & fdis
                     
                  ! perturbed parallel energy
                  epar = epar + bn(ix,i)*intvp(i,j,k,is)*intmu(j)*vpgr(i,j,k,is)**2* &
                     & fdis
                     
                  ! perturbed perpendicular energy
                  eperp= eperp + bn(ix,i)**2*intvp(i,j,k,is)*intmu(j)*2.E0*mugr(j)*    &
                     & fdis

                end do
              end do
            end if
            ! reduce over points of equal s; any procs not responsible for the
            ! current species will use zero values.
            if (number_of_processors > 1) then
              arr(1) = dens ; arr(2) = epar ; arr(3) = eperp ; arr(4) = parallel_velocity
              call mpiallreduce_sum(arr,rarr,4,COMM_S_EQ_X_EQ)
              dens = rarr(1) ; epar = rarr(2); eperp = rarr(3) ; parallel_velocity = rarr(4)
            end if

            !Rotate all values in the complex plane relative to maximum potential
            !defined to be (1.,0.).  Effectively, normalise out frequency.
            if(nlphi) then
              phid=phi(imod,ix,i)
            else
              phid = 0.0
            end if
            if(nlapar) then
              apad=apar(imod,ix,i)
            else
              apad = 0.0
            end if
            if(nlbpar) then
              bpad=bpar(imod,ix,i)
            else
              bpad = 0.0
            end if
            if (.not.non_linear .and. lrotate_parallel) then
              apad=apad/rotate
              phid=phid/rotate
              bpad=bpad/rotate
              dens=dens/rotate
              epar=epar/rotate
              eperp=eperp/rotate
              parallel_velocity=parallel_velocity/rotate
            endif

            ! Store for the columns of the file         
            local_write(ix,i,1) = sgr(i)
            local_write(ix,i,2) = phid
            local_write(ix,i,3) = apad
            local_write(ix,i,4) = dens
            local_write(ix,i,5) = epar
            local_write(ix,i,6) = eperp
            local_write(ix,i,7) = parallel_velocity
            local_write(ix,i,8) = bpad

            ! debug output
            !local_write(ix,i,:) = processor_number

            ! this barrier appears to be necessary ??
            call mpibarrier()      

          end do ! ns
        end do ! nx

        ! gather (only for the root processor) for all 8 complex columns of the file
        do il = 1,8     
          call gather_array(global_write(:,:,il),n_x_grid, & 
             & n_s_grid,local_write(:,:,il),nx,ns,COMM_X_NE,COMM_S_NE)
        end do

        if (root_processor) then
          ! output
          ! real: sgr(i), complex: phid, apad, dens, epar, eperp, parallel_velocity, bpad

          do ix = 1, n_x_grid
            do i = 1, n_s_grid
              !       write(ifile,fmt = '(15(1x,e13.5))') real(global_write(ix,i,1)), &
              !                          & (global_write(ix,i,il), il = 2,8)
              call append_chunk(lun, &
                 & (/ real(global_write(ix,i,1)), &
                 & (real(global_write(ix,i,il)), aimag(global_write(ix,i,il)), il=2,8) /), &
                 & '(15(1x,e13.5))', ascii_fmt)
            end do
          end do
        end if

      end do ! nmod
    end do  ! species

    call close_lu(lun, ascii_fmt)

    if ((.not.non_linear).and.root_processor.and.lrotate_parallel) then
      write(*,*) 
      write(*,*) 'parallel.dat: normalized (rotated in complex plane) by: '
      write(*,*) rotate 
      write(*,*)
    end if

  end subroutine parallel_output



end module diagnos_mode_struct
