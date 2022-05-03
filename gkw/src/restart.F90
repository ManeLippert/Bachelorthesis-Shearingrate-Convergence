!------------------------------------------------------------------------------
!> This module contains routines for reading and writing the restart file
!------------------------------------------------------------------------------
module restart

  implicit none

  private

  public :: restart_init, read_restart_file, write_restart_file

  !> Set to true if code has successfully restarted
  logical, public, save :: restarted=.false.

  !> restart file gridsizes if different from input file
  !> (lrestart_new_grid only)
  integer, save :: nx_rst, nmod_rst, number_of_species_rst

  !> buffer for writing fdisi in the right order for restart
  complex, allocatable, save, dimension(:) :: fdisiiobuf

  integer, public, save :: TYPE_RW_FDISI          !< fdisi to/from file
  
  !use one routine for reading and writing
  interface restart_write_nml
    module procedure restart_read_nml
  end interface


contains

  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine restart_init()
    use mpiinterface, only : set_comm_fdisi_io
    ! setup communicator to read / write the restart file (after setup_grid) 
    call set_comm_fdisi_io()
    ! Generate a new datatype TYPE_FDISI_IO for a mapping between the (variable)
    ! memory layout and the (presently fixed) restart file layout.
    call create_fdisi_io_type()

  end subroutine restart_init


  !-----------------------------------------------------------------------------
  !> This routine returns the prefix for restart filenames, based on an
  !> enumeration index
  !-----------------------------------------------------------------------------
  subroutine restartfilename(prefix,ival,ltavg)
    integer, intent(in)            :: ival
    logical, intent(in)            :: ltavg
    character (len=3), intent(out) :: prefix

    if(ltavg) then
      prefix='FTA'
    else
      if (ival == -345) then
        prefix='INI'
      else if (ival.lt.1) then
        prefix='FDS'
      else if (ival.lt.10) then
        write(prefix,'("FD",I1.1)') ival
      else if (ival.lt.100) then
        write(prefix,'("F",I2.2)') ival
      else
        prefix='FDS' 
      end if
    end if

  end subroutine restartfilename


  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine write_restart_file(lfinal, irun, ltavg)
    logical, intent(in) :: lfinal, ltavg
    integer, intent(in) :: irun

    call write_restart_file_binary(lfinal, irun, ltavg)

  end subroutine write_restart_file

  !--------------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------------
  subroutine read_restart_file(irun)
    use control, only : read_file, auto_restart
    use general, only : gkw_abort
    integer, intent(in) :: irun
    
    if (.not.(read_file .or. auto_restart)) then
      call gkw_abort('Invalid call to read_restart_file')
    end if

    call read_restart_file_binary(restarted, irun)

  end subroutine read_restart_file


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Write file(s) for restarting.
!>
!> In the latest version (2) this is to a single
!> file "FDS" via MPI IO. This allows restarts on a different number of
!> processors.
!> In the original version (1), one file per processor is output as
!> a fortran binary file, so this must be run again on the same number of
!> processors.
!>
!> Some MPI implementations do not support writing of the derived
!> datatypes used in the code, so in that case the original version
!> (restart_file_version = 1) must be selected in the control input namelist.
!----------------------------------------------------------------------------
subroutine write_restart_file_binary(lfinal, irun, ltavg)

  use global,         only : id_vpar,id_mu,id_s,id_x,id_mod,id_sp
  use mpiinterface,   only : mpibarrier, number_of_processors, statusmpi
  use mpiinterface,   only : last_processor, data_representation
  use mpiinterface,   only : processor_number, root_processor
  use mpiinterface,   only : MPICOMPLEX_X, MPI_OFFSET_KIND
  use mpiinterface,   only : MPI_MODE_WRONLY, MPI_INFO_NULL
  use mpiinterface,   only : MPI_MODE_CREATE, mpiwtime
  use mpidatatypes,   only : create_subarray_datatype
  use mpidatatypes,   only : free_datatype
  use mpicomms,       only : COMM_FDISI_IO, COMM_FDISI_IO_nprocsets
  use mpicomms,       only : COMM_FDISI_IO_procset
  use grid,           only : nsp,ns,nvpar,nmu,nmod,nx
  use dist,           only : ifdis, fdisi, nf, fdisi_tavg
  use control,        only : restart_file_version, max_t_fdis_write
  use io,             only : lmpi_broken_io, get_free_file_unit
  use index_function, only : indx
  use control,        only : nan_stop, nt_complete
  use mpiinterface,   only : mpibcast, send_to_root
  
    
  logical, intent(in) :: lfinal, ltavg
  integer, intent(in) :: irun

  integer (kind=MPI_OFFSET_KIND) :: idisp
  character (len=3) :: root
  character (len=9) :: filename, filename2
  integer :: np,file_unit,ierr,new_type
  integer :: i,j,k,l,m,n,ind,ierr_t,iprocset
  integer, save :: idump=0
  real :: t_beg, t_write

  ierr=0
  ierr_t=0

  t_beg = mpiwtime()

  ! Determine the output file name
  if (lfinal) then   
    ! generate the file name root, 'FDS' is default if irun=0
    ! FD1 if irun=1, etc
    call restartfilename(root,irun,ltavg)
  else
    !alternate dump files, 1 and 2 - autorestart uses 1 for now
    idump = idump + 1
    root = 'DM1'
    if (mod(idump,2)==0) root='DM2'
  end if

  if(nan_stop)restart_file_version=0

  !do something based on the restart file version
  restartfileversion : select case(restart_file_version)
  
    case(0) 
      
      return   ! write no restart files
    
    case(1)   ! restart_file_version 1: multiple fortran raw binary files
      
      call genfilename(root,filename,processor_number)

      ! In this method a barrier is used to prevent all processors outputing
      ! the data in one go. This might not always be necessary, and there
      ! should be more elegant solutions.
      call mpibarrier()

      do np = 1, number_of_processors

        if (np-1 == processor_number) then
          call get_free_file_unit(file_unit)
          open(file_unit,FILE = filename, STATUS= 'unknown', FORM = 'unformatted')
          write(file_unit)(fdisi(i), i = 1, nf)
          close(file_unit)
        end if
        call mpibarrier()

      end do 

    case(2)   ! restart_file_version 2: single binary file "FDS"

#if defined(mpi2)

  !
  ! Create a sub-array dataype so that the local part can be written to the
  ! right part of the global array. This should be consistent with the call
  ! in the routine reading the restart file.
  ! The order used in the call here is consistent with the order used in
  ! create_fdisi_io_type, since the data is always ordered in the same way in
  ! the file (although this restriction is not entirely necessary).
  !

  call create_subarray_datatype(MPICOMPLEX_X,new_type,                       &
      &                         id_vpar,id_mu,id_s,id_x,id_mod,id_sp)

  !
  ! Open and write to a file with all processors, or alternatively perform an
  ! open+write+close with some subset of the processors at a time (which can
  ! be useful when there are otherwise problems with MPI IO). COMM_FDISI_IO is
  ! a communcator over which the file can be opened with a subset of processors
  ! at a time. The number of independent subsets is COMM_FDISI_IO_nprocsets
  ! and the local processor belongs to subset number COMM_FDISI_IO_procset.
  ! This communicator will actually be COMM_S_EQ or COMM_SP_EQ_S_EQ or some
  ! other communicator that corresponds to a subdomain of the problem
  ! decomposition.  The communicator size is determined in set_comm_fdisi_io()
  !

  do iprocset = 1, COMM_FDISI_IO_nprocsets

    ! global barrier - need to wait till file is closed before re-open
    call mpibarrier()

    ! wait for my turn to open the file and write
    my_turn : if (COMM_FDISI_IO_procset == iprocset) then

      ! First processor set creates a new file; other sets just write.
      if (iprocset == 1) then
        call MPI_FILE_OPEN(COMM_FDISI_IO,root,MPI_MODE_WRONLY+MPI_MODE_CREATE&
            & ,MPI_INFO_NULL,file_unit,ierr)
      else
        call MPI_FILE_OPEN(COMM_FDISI_IO,root,MPI_MODE_WRONLY,               &
            &  MPI_INFO_NULL,file_unit,ierr)
      end if

      ! Do something on file open error?
      if (ierr /= 0) write(*,*) processor_number,                            &
          &                     ': WARNING: Error on MPI file OPEN'
      ierr_t = ierr + ierr_t


#if defined(FORCE_IO_ATOMICITY)
      ! force sequential I/O to prevent deadlock
      call MPI_File_set_atomicity(file_unit,.true.)
#endif


      ! Set the file displacement to be used with the file view to be zero
      idisp = 0

      ! I only get to see the bits of the full global array that I need to
      ! write to.
      call MPI_FILE_SET_VIEW(file_unit,idisp,MPICOMPLEX_X,new_type,                &
          &                  data_representation,MPI_INFO_NULL,ierr)

      ! Write to file: some MPI implementations do not support writing derived
      ! datatypes; if so we re-order the data to write in a temporary array
      ! before writing, so that it can be written as a standard type.
      if (lmpi_broken_io) then

        allocate (fdisiiobuf(nf),stat=ierr)

        ! Re-order the data; the order of these loops should not be changed!
        ind = 0
        do n=1,nsp; do m=1,nmod; do l=1,nx; do k=1,ns; do j=1,nmu
          do i=1,nvpar
            ind = ind + 1
            if(ltavg) then
              fdisiiobuf(ind) = fdisi_tavg(indx(ifdis,m,l,k,j,i,n))
            else
              fdisiiobuf(ind) = fdisi(indx(ifdis,m,l,k,j,i,n))
            end if
          end do
        end do ;    end do ;     end do ;   end do ;   end do

        ! write the temporary array as nf complex values
        call MPI_FILE_WRITE_ALL(file_unit,fdisiiobuf(1),nf,MPICOMPLEX_X,statusmpi,ierr)

        deallocate (fdisiiobuf)

      else

        ! Usually I can just write and re-order fdisi without difficulty.
        call MPI_FILE_WRITE_ALL(file_unit,fdisi,1,TYPE_RW_FDISI,statusmpi,ierr)

      end if

      ! Do something on write close error?
      if (ierr /= 0) write(*,*) processor_number,                            &
          &                     ': WARNING: Error on MPI file write'
      ierr_t = ierr + ierr_t

      ! barrier for the writing processors
      call mpibarrier(COMM_FDISI_IO)

      ! close the file
      call MPI_FILE_CLOSE(file_unit,ierr)

      ! Do something on file closer error?
      if (ierr /= 0) write(*,*) processor_number,                            &
          &                     ': WARNING: Error on MPI file close'
      ierr_t=ierr+ierr_t

    end if my_turn

  end do

  ! Really there is no need to create and free these everytime. It should be
  ! done once on the first call or before time integration.
  call free_datatype(new_type)

!    Check for errors accross all processors
!    call mpiallreduce_max(ierr_t,ierr,1)
!    if (ierr /=0) then retry write once or twice

#endif

  end select restartfileversion

  ! Now write associated data about FDS dump file
  if (last_processor) then

    filename2=root//'.dat'
    filename2=trim(filename2)
    call write_restart_info(filename2, root, ierr)
        
  end if
  
  
  ! The lines below make sure that nt_complete is transferred to all 
  ! processes. First, it is send from last_processor (restart file info
  ! has been written by last_processor) to root_processor. Second, 
  ! nt_complete is broadcast from root_processor to all processors.
  ! This is necessary for consistent numbering of diagnostic output 
  ! files, which relies on nt_complete.
  if((last_processor .and. .not. root_processor) .or. &
    &(.not. last_processor .and. root_processor)) then
    call send_to_root(nt_complete)
  end if
  call mpibcast(nt_complete,1)
  

  t_write = mpiwtime() - t_beg
  max_t_fdis_write = max(t_write,max_t_fdis_write)

  if (root_processor) then
    write(*,*) '* wrote restart file ',root,', count=', nf
    write(*,*) '* restart_file_version=', restart_file_version
    write(*,'(A,es13.5,A)') ' * in ', t_write, ' seconds'
    write(*,*) 
  end if

end subroutine write_restart_file_binary

!-------------------------------------------------------------------------
!> Write a file with additional information associated to the
!> dump/checkpoint/restart file.
!-------------------------------------------------------------------------
subroutine write_restart_info(filename, dumpfilename, ierr)
  use io, only : get_free_file_unit
  use control, only : nan_stop
  use rotation, only : kxshift_write_nml
  use version, only : GKW_REV, GKW_EXE
  use mpiinterface, only : number_of_processors
  use control, only : restart_file_version

  character (len=*), intent(in) :: filename, dumpfilename
  integer, intent(in) :: ierr
  integer :: dat_file_unit, io_stat
  
  call get_free_file_unit(dat_file_unit)
  open(dat_file_unit,FILE=filename,FORM='formatted',STATUS='unknown')

  !Write file header.
  write(dat_file_unit,*) '!Dump filename: ', dumpfilename
  write(dat_file_unit,*) '!Written with restart file version: ', &
     & restart_file_version
  if (restart_file_version==2) then
    write(dat_file_unit,*) '!Write status error code: ', ierr
  end if
  write(dat_file_unit,*) '!Restart file written with GKW version: ', GKW_REV
  write(dat_file_unit,*) '!Executable name: ',GKW_EXE
  write(dat_file_unit,*) '!Run on ', number_of_processors, 'processor(s)'
  write(dat_file_unit,*) '!Values in this file will override input.dat'
  if(nan_stop)then
    write(dat_file_unit,*) '!No FDS file (over)written as code is unstable (NaNs)'
  endif
  !Write restart namelist 
  call restart_write_nml(dat_file_unit,io_stat,.true.)
  !And store kxshift from rotation
  call kxshift_write_nml(dat_file_unit,io_stat)
  close(dat_file_unit)
end subroutine write_restart_info



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Read the restart files(s).
!----------------------------------------------------------------------------
subroutine read_restart_file_binary(success,irun)

  use general,        only : gkw_abort, gkw_warn
  use global,         only : id_mod, id_mu, id_s, id_sp, id_sp, id_vpar, id_x
  use io,             only : get_free_file_unit, lmpi_broken_io
  use dist,           only : ifdis, fdisi,nf
  use index_function, only : indx
  use control,        only : auto_restart, lrestart_new_grid
  use grid,           only : nx,ns,nmod,nvpar,nmu,nsp, number_of_species
  use grid,           only : parallel_sp, parallel_x, n_x_grid
  use rotation,       only : kxshift_read_nml, kxshift_bcast_nml
  use mpiinterface,   only : root_processor, processor_number, statusmpi
  use mpiinterface,   only : number_of_processors, mpibarrier
  use mpiinterface,   only : data_representation
  use mpiinterface,   only : MPICOMPLEX_X, MPI_MODE_RDONLY
  use mpiinterface,   only : MPI_OFFSET_KIND, MPI_INFO_NULL
  use mpidatatypes,   only : create_subarray_datatype
  use mpicomms,       only : COMM_FDISI_IO_nprocsets, COMM_FDISI_IO_procset
  use mpicomms,       only : COMM_FDISI_IO

  logical, intent(out) :: success
  integer, intent(in) :: irun

  integer (kind=MPI_OFFSET_KIND) :: idisp
  integer :: NEW_TYPE, new_type_rst
  logical lread_input
  character (len=3) :: root 
  character (len=9) :: filename2, name 
  integer :: np,ierr,ierr_t,io_stat,file_unit,version,dum,iprocset
  integer :: nsp_dum, nx_dum, nmod_dum, nf_dum, ix0_df

  integer :: i,j,k,l,m,n,ind

  success = .false.

!Still a potential problem if incorrect input files are provided
!This is hard to make foolproof but should be solved by the rewrite 
!of the restarts allowing changed parallel configurations

  
if (lrestart_new_grid .and. .not. lmpi_broken_io) then
  call gkw_abort('Restart on new grid requires lmpi_broken_io') 
end if 

!Decide which file and type to read. restart_files_exists returns 0,1, or 2
!Version 1 takes priority over 2 if it exists
!DMP files take priority over FDS files if they exist.
version=0
if (auto_restart) then !DMP file takes preference in case of auto_restart
   version=restart_files_exist('DM1') 
end if

if (version==1 .or. version==2) then 
  root='DM1'
else if (version==0) then ! There is no dump file
  ! Generate restart filename from irun. FDS is default.
  call restartfilename(root,-345,.false.)
  ! look for FDS file
  version=restart_files_exist(root)
  if (version==0) then
    ! Generate restart filename from irun. FDS is default.
    call restartfilename(root,irun-1,.false.)
    ! look for FDS file
    version=restart_files_exist(root)
  end if
  if (version==0) then  ! No file has been found
    if(root_processor) then
      write(*,*) 'No restart file ',root,' found, intialising new run'
    end if
    success = .false.
    return
  end if
end if

if (lrestart_new_grid .and. version/=2) then 
  call gkw_abort('Restart on new grid requires restart file of version 2')
end if

if (root_processor) then
   write(*,*) 'Reading restart file: ', root
   write(*,*) 'Data will be appended to existing files.'
end if

! READ restart namelist if file exists, otherwise ignore.
  !write(filename2,*) root,'.dat'
  filename2=root//'.dat'
  filename2=trim(filename2)
  inquire(FILE=filename2,EXIST=lread_input)
  if (lread_input) then
    if (root_processor) then
      call get_free_file_unit(file_unit)
      open(file_unit,file=filename2,FORM='formatted',STATUS='unknown')
      call restart_read_nml(file_unit,io_stat,.false.)
      call kxshift_read_nml(file_unit,io_stat,.false.)
      close(file_unit)
    end if
    !Line below causes circular dependeny, omit for now.
    !call namelist_error_check('restart',io_stat)
    call restart_bcast_nml
    call kxshift_bcast_nml
    !if (last_processor) call restart_write_nml(ilun_out,io_stat,.true.)
  else if (root_processor) then
    write(*,*) 'WARNING: No companion file ', filename2, ' found, using defaults' 
  end if

call mpibarrier()

select case(version)

  case(1) ! Individual restart files
    do np = 1 , number_of_processors 

     dum = processor_number
     ! Existence of all files has already been checked for
     ! in restart_files_exist
     call genfilename(root,name,dum)
  
      if (np-1.eq.processor_number) then 
        open(11,file = name, status = 'unknown', form = 'unformatted')
          read(11)(fdisi(i), i = 1, nf)
          write(*,*)nf
        close(11)
      end if 

      ! At present a barrier is used to prevent all processors of reading 
      ! in all the data in one go. This might not always be necessary, 
      ! and there should be more elegant solutions       
      call mpibarrier()
  
    end do 
  
  case(2) ! Single restart file

#if defined(mpi2)

  !
  ! Create a sub-array dataype so that the local part can be read from the
  ! right part of the global array. This should be consistent with the call
  ! in the routine writing the restart file.
  ! The order used in the call here is consistent with the order used in
  ! create_fdisi_io_type, since the data is always ordered in the same way in
  ! the file (although this restriction is not entirely necessary).
  !

  call create_subarray_datatype(MPICOMPLEX_X,new_type,                       &
      &                         id_vpar,id_mu,id_s,id_x,id_mod,id_sp)
      
      
  ! if run is restarted on a new grid, create a sub-array datatype which is
  ! compatible with the global (unchanged) grid layout of the old run read 
  ! from restart namelist.
  if(lrestart_new_grid) then
    call create_subarray_datatype(MPICOMPLEX_X,new_type_rst,        &
        &                         id_vpar,id_mu,id_s,id_x,id_mod,id_sp, &
        & global_spsize=number_of_species_rst, global_xsize=nx_rst,     &
        & global_ysize=nmod_rst)
  end if
  

  !
  ! Open and read from a file with all processors, or alternatively perform an
  ! open+read+close with some subset of the processors at a time (which can
  ! be useful when there are otherwise problems with MPI IO). COMM_FDISI_IO is
  ! a communcator over which the file can be opened with a subset of processors
  ! at a time. The number of independent subsets is COMM_FDISI_IO_nprocsets
  ! and the local processor belongs to subset number COMM_FDISI_IO_procset.
  ! This communicator will actually be COMM_S_EQ or COMM_SP_EQ_S_EQ or some
  ! other communicator that corresponds to a subdomain of the problem
  ! decomposition.
  !

  ierr_t=0

  do iprocset = 1, COMM_FDISI_IO_nprocsets

    ! global barrier - need to wait till file is closed before re-open
    call mpibarrier()

    ! wait for my turn to open the file and write
    my_turn : if (COMM_FDISI_IO_procset == iprocset) then
      ierr=0
      call MPI_FILE_OPEN(COMM_FDISI_IO,root,MPI_MODE_RDONLY,MPI_INFO_NULL,   &
          &              file_unit,ierr)

      ! Do something on file open error?
      if (ierr /= 0) write(*,*) processor_number,                            &
          &                     ': WARNINNG: Error on MPI file OPEN'
      ierr_t = ierr + ierr_t

      ! Set the file displacement to be used with the file view to be zero
      idisp = 0

      ! I only get to see the bits of the full global array that I need to
      ! read from. If run is restarted on new grid use datatype layout 
      ! corresponding to old global grid to set file view.
      if(lrestart_new_grid) then
        call MPI_FILE_SET_VIEW(file_unit,idisp,MPICOMPLEX_X,new_type_rst,      &
            &                  data_representation,MPI_INFO_NULL,ierr)
      else
        call MPI_FILE_SET_VIEW(file_unit,idisp,MPICOMPLEX_X,new_type,          &
            &                  data_representation,MPI_INFO_NULL,ierr)
      end if
      
      ! Read from file: some MPI implementations do not support reading
      ! derived datatypes; if so we read the data into a temporary
      ! array before re-ordering, so that it can be read as a standard type.
      if (lmpi_broken_io) then

        if (lrestart_new_grid) then ! grid is different size in restart file
          if (parallel_sp.and. number_of_species_rst/=number_of_species) then
             call gkw_abort('Cannot restart on this new grid with parallel species')
          end if
          
          if (parallel_x.and. nx_rst/=n_x_grid) then
             call gkw_abort('Cannot restart on this new grid with parallel x')
          end if
          
!          nsp_dum = number_of_species_rst
          if(parallel_sp) then
            ! the above gkw_abort makes sure that with parallel_sp the current
            ! global species grid has to agree with the original one, then:
            nsp_dum = nsp
          else
            nsp_dum = number_of_species_rst
          end if
          if(parallel_x) then
            ! the above gkw_abort makes sure that with parallel_x the current
            ! global xgrid has to agree with the original one, then:
            nx_dum = nx
          else
            nx_dum  = nx_rst
          end if
          nmod_dum = nmod_rst
          ! The size of the restart solution, without any fields
          nf_dum = nsp_dum*nx_dum*ns*nmu*nvpar*nmod_dum
          !shift in the zero x mode
          ix0_df  = (nx_dum+1)/2 - (nx+1)/2                
       
       
          ! Since grid dimensions specified in restart namelist are global
          ! grid dimensions compare to current global grid dimensions to 
          ! identify, if grid has been modified.
          if (number_of_species_rst > number_of_species) &
          & call gkw_warn('Restarting with excess species dropped')
          if (nx_rst > n_x_grid) call gkw_warn('Restarting with smaller x grid')
          if (nmod_rst > nmod) call gkw_warn('Restarting with smaller y grid')
        
          if (number_of_species_rst < number_of_species) &
          & call gkw_warn('Restarting with new species added')
          if (nx_rst < n_x_grid) call gkw_warn('Restarting with enlarged x grid')
          if (nmod_rst < nmod) call gkw_warn('Restarting with enlarged y grid')
          
          if (nmod /= nmod_dum .and. (nmod==1 .or. nmod_dum == 1)) then
            call gkw_abort('Cannot restart when nmod changed to or from 1')
          end if
           
        else ! normal restart
          nsp_dum = nsp; nx_dum = nx; nmod_dum = nmod; nf_dum = nf; ix0_df = 0
        end if
        
        ! allocate the buffer used to read the part of the restart file that
        ! correspond to current processor subset
        allocate (fdisiiobuf(nf_dum),stat=ierr)

        ! read the temporary array as nf_dum complex values
        call MPI_FILE_READ_ALL(file_unit,fdisiiobuf(1),nf_dum,MPICOMPLEX_X,statusmpi,ierr)
          
        ! Re-order the data; the order of these loops should not be changed!
        ind = 0
        do n=1,nsp_dum; do m=1,nmod_dum; do l=1,nx_dum; do k=1,ns; do j=1,nmu
          do i=1,nvpar
            ind = ind + 1
            !Line below allows dropping of grid point in reduced grid restart
            if (n > nsp .or. l-ix0_df > nx .or. m > nmod .or. l-ix0_df < 1) cycle
            fdisi(indx(ifdis,m,l-ix0_df,k,j,i,n)) = fdisiiobuf(ind)
          end do
        end do ;    end do ;     end do ;   end do ;   end do

        deallocate (fdisiiobuf)

      else

        ! Usually I can just write and re-order fdisi without difficulty.
        call MPI_FILE_READ_ALL(file_unit,fdisi,1,TYPE_RW_FDISI,statusmpi,ierr)

      end if

      ! Do something on read close error?
      if (ierr /= 0) write(*,*) processor_number,                            &
          &                     ': WARNING: Error on MPI file read'
      ierr_t = ierr + ierr_t

      ! barrier for the writing processors
      call mpibarrier(COMM_FDISI_IO)

      ! close the file
      call MPI_FILE_CLOSE(file_unit,ierr)
      !call free_datatype(new_type)
      ! Do something on file close error?
      if (ierr /= 0) write(*,*) processor_number,                            &
          &                     ': WARNING: Error on MPI file close'
      ierr_t=ierr+ierr_t

    end if my_turn

  end do

#else
    call gkw_abort('Single restart file requires MPI')
#endif
  case default
    call gkw_abort('function: restart_files_exist: error') 

end select !File version.

  !restarted = .true.
  success = .true.

end subroutine read_restart_file_binary


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Checks for correct number and existance of restart files 
!> with given prefix name. 
!> returns 0 if no restart files exist
!> returns 1 if type 1 exists, returns 2 if only type 2 exists
!> Aborts if incorrect number of single files are found.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function restart_files_exist(root)

  use mpiinterface, only : processor_number, mpiallreduce_and
  use general,      only : gkw_abort

  integer :: restart_files_exist

  character (len=3), intent(IN) ::  root  !The filename (or prefix)

  character (len=9) :: firstfile, name
  logical :: docont, individual_restart_files
  logical :: all_files_present

  restart_files_exist=0

  ! 1) Check for first of individual restart files
  call genfilename(root,firstfile,0)
  inquire(file=firstfile,exist = individual_restart_files)
  if (individual_restart_files) restart_files_exist=1

  ! 2) Check for the single restart file
  if (.not. individual_restart_files) then
    inquire(file=root,exist = docont)
    if (docont) then
      restart_files_exist=2
      individual_restart_files = .false.
    end if
  end if

  ! for version 1 restarts, check all the files exist
  if (individual_restart_files) then
    call genfilename(root,name,processor_number)

    ! find out if the file exists for this processor
    docont=.false.
    inquire(file=name,exist = docont)

    ! check that all processors have a restart file
    call mpiallreduce_and(docont,all_files_present,1)

    if (.not. all_files_present) then
      call gkw_abort('Not enough input files for restart')
    end if

  else
    ! do nothing -- only 1 file; nothing to check
  end if

end function restart_files_exist


!----------------------------------------------------------------------------
!> read or write the restart namelist
!----------------------------------------------------------------------------
subroutine restart_read_nml(file_unit,io_stat,lwrite)
  use control, only : ntime, time, itime, nt_complete, dtim, method
  use control, only : itime_rst
  use global, only : r_tiny, r_huge, int2char, real2char
  use general, only : gkw_warn
  use grid, only : n_x_grid, nmod, number_of_species
  use io, only : open_real_lu, ascii_fmt, xy_fmt, attach_metadata, append_chunk
  use io, only : close_lu, description_key


  
  integer, intent(in) :: file_unit
  integer,intent(out) :: io_stat
  logical, intent(in) :: lwrite
  real :: dtim_current
  integer :: nt_remain, file_count

  integer :: lun
  
  namelist /restart/ dtim, &
     & time, &
     & nt_complete, nt_remain, & !,irun
     & nx_rst, nmod_rst, number_of_species_rst, &
     ! obsolete:
     & file_count

  !At present irun has to be manually incremented in input.dat for next run.
  !This preserves original behaviour of overwriting the restart file
  !If irun is not manually set by the user

  !For automatic restarts to be seamless, code would need to
  !seek back to correct line of time trace files
  
  if (.not. lwrite) then
    time = -r_huge
    
    
    ! save the current runs dtim
    dtim_current = dtim
    ! read the namelist
    read(file_unit,NML=restart,IOSTAT=io_stat)

    if (nt_remain > 0) then
      ntime=nt_remain
      call gkw_warn('This is a restart: ntime='//int2char(nt_remain)//' set to &
         & complete remaining time of previous run')
    end if

    if (abs(dtim-dtim_current) > 100.*r_tiny) then
      call gkw_warn('dtim was different: '//real2char(dtim)//' in previous run')
    end if
    !don't overwrite input dtim with the previous one
    dtim = dtim_current

    if(abs(time - (-r_huge)) < r_tiny) then
      call gkw_warn("The restart parameters lack the 'time' value. It is set to 0.")
      time = 0.0
    end if
    write(*,*) 'Resuming at normalised time: ', time

    call open_real_lu('nt_complete', 'restart', (/ 1 /), &
       & ascii_fmt, lun)
    call append_chunk(lun, (/ real(nt_complete) /), xy_fmt, ascii_fmt)
    call attach_metadata(lun, &
       & description_key, &
       & 'Each entry in this dataset is the value nt_complete as it appears in &
       & the restart input file, for the past restarts of this run. &
       & Example: Values 4, 10, 17 mean that the first run ended after 4 &
       & timesteps, then it was restarted and ran for another 6 timesteps, &
       & restarted again and ran for another 7 steps (and then ran again, but &
       & values are appended only at restart). Note that the last entry &
       & appears also as metadata attribute restart_parameters.nt_complete .' &
       & , ascii_fmt)
    call close_lu(lun, ascii_fmt)

    call open_real_lu('nt_remain', 'restart', (/ 1 /), &
       & ascii_fmt, lun)
    call append_chunk(lun, (/ real(nt_remain) /), xy_fmt, ascii_fmt)
    call attach_metadata(lun, &
       & description_key, &
       & 'Each entry in this dataset is the value nt_remain as it appears in &
       & the restart input file, for the past restarts of this run. &
       & Example: Consider ntime to be fixed to 12. Values 8, 2, 7 mean that &
       & the first run ended after 4 &
       & timesteps (8 steps before completing the run), then it was restarted &
       & and ran for another 6 timesteps (2 steps left), &
       & restarted again and ran for another 7 steps (and 7 steps left) &
       & (and then ran again, but &
       & values are appended only at restart). Note that the last entry &
       & appears also as metadata attribute restart_parameters.nt_remain .'&
       & , ascii_fmt)
    call close_lu(lun, ascii_fmt)

    ! a different name than 'time' is needed, because otherwise the
    ! time lu of diagnos_grid would be overwritten for ascii or binary
    call open_real_lu('time_complete', 'restart', (/ 1 /), &
       & ascii_fmt, lun)
    call append_chunk(lun, (/ time /), xy_fmt, ascii_fmt)
    call attach_metadata(lun, &
       & description_key, &
       & 'Each entry in this dataset is the value time as it appears in &
       & the restart input file, for the past restarts of this run. The real &
       & values of time_complete correspond to the integer timestep &
       & values nt_complete.', &
       & ascii_fmt)
    call close_lu(lun, ascii_fmt)

    call open_real_lu('dtim', 'restart', (/ 1 /), &
       & ascii_fmt, lun)
    call append_chunk(lun, (/ dtim /), xy_fmt, ascii_fmt)
    call attach_metadata(lun, &
       & description_key, &
       & 'Each entry in this dataset is the value dtim as it appears in &
       & the restart input file, for the past restarts of this run.', &
       & ascii_fmt)
    call close_lu(lun, ascii_fmt)

  else !write the namelist
    nx_rst=n_x_grid
    nmod_rst=nmod
    number_of_species_rst=number_of_species
    !nperiod_rst=nperiod
    !n_s_grid_rst=n_s_grid


    ! If restart file is written by Eigenvaluesolver, use old 
    ! implementation.
    if(method == 'EIV') then
      nt_complete = itime + nt_complete
      
    ! Otherwise: itime_rst counts the number of large time steps since 
    ! the last restart file has been written. This method allows for a
    ! consistent numbering, when dump files are automatically written 
    ! during runtime.
    else
      nt_complete = itime_rst + nt_complete
    end if
        
    
    nt_remain=ntime-itime
    !increment the run number for the next run
    !irun=irun+1
    

    !obsolete:
    file_count = nt_complete
    
    write(file_unit,NML=restart)
    !but keep irun unchanged for this run.
    !irun=irun-1
  end if

end subroutine restart_read_nml


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> broadcast the restart namelist
!----------------------------------------------------------------------------
subroutine restart_bcast_nml

  use control,      only : irun, ntime, dtim, time, t_init, nt_complete 
  use mpiinterface, only : mpibcast 

  call mpibcast(irun,1)
  call mpibcast(ntime,1)
  call mpibcast(dtim,1)
  call mpibcast(time,1)
  call mpibcast(nt_complete,1)
  call mpibcast(nx_rst,1)
  call mpibcast(nmod_rst,1)
  call mpibcast(number_of_species_rst,1)
  !Note file_count read from the file is not used
  t_init=time

end subroutine restart_bcast_nml

!-----------------------------------------------------------------------------
!> this routine returns a filename `filename', starting with the first 3
!> characters of `start', followed by a 6 digit representation of the input
!> integer `ival'. ival is inc. by 1 on return
!-----------------------------------------------------------------------------
subroutine genfilename(start,filename,ival)

  character (len=3), intent(in)  :: start
  integer, intent(in) :: ival
  character (len=9), intent(out) :: filename

  write(filename,'(A,I6.6)') start, ival

end subroutine genfilename

!****************************************************************************
!> Create a datatype for reading/writing the distribution function from/to
!! file in a particular order. The distribution function can be stored in
!! any order in memory, but for the purpose of restarting a run with an
!! alternative memory layout it is useful to always have the file data in the
!! same order. The order of the data in the file could also be recorded in the
!! file so that it can always be correctly read, if for some reason is useful
!! to have differently ordered outputs.
!!
!! Presently the input `type_old' for this routine is MPICOMPLEX_X and the
!! returned type `type_new' is the then the type to write a the full local
!! This is particularly useful for writing the distribution function
!! APS: subroutine get_reordered_type(ndims,dims,type_in,type_out)
!! returns old_type for non-MPI run
!<----------------------------------------------------------------------------
subroutine create_fdisi_io_type()
  use dist, only : nf, ifdis
  use global,         only : root_and_verbose 
  use general,        only : gkw_abort
  use index_function, only : indx
  use mpiinterface,   only : MPICOMPLEX_X
  use grid,           only : nsp,nvpar,nmu,ns,nx,nmod
  integer, allocatable, dimension(:) :: offset
  integer :: blocklen
  integer, dimension(6) :: iend
  integer :: ierr, new_index
  integer :: ispecies, ivpar,imu,is,ix,imod
  integer :: i,j,k,l,m,n
  integer :: u_count

  u_count = 0
  blocklen = 1
  allocate(offset(nf),stat=ierr) 
  if (ierr /= 0) call gkw_abort('get_reordered_type: cannot allocate offset')

  iend(:) = 1
  !APS  iend(1:n_dims) = dims(1:n_dims)
  iend(1) = nsp
  iend(2) = nmod
  iend(3) = nx
  iend(4) = ns
  iend(5) = nmu
  iend(6) = nvpar

  new_index = 0
  do n= 1, iend(1)
    do m= 1, iend(2)
      do l= 1, iend(3)
        do k= 1, iend(4)
          do j = 1, iend(5)
            do i = 1, iend(6)
              !
              ispecies = n
              imod     = m
              ix       = l
              is       = k
              imu      = j
              ivpar    = i
              !ivpar,imu,is,ix,imod,isp
              new_index = new_index + 1
              ! starts from zero in mpi call

              offset(new_index) = indx(ifdis, imod, ix, is, imu, ivpar, ispecies) - 1
              if (offset(new_index) + 1 == new_index) u_count = u_count + 1
              !blocklen(new_index) = 1
            end do
          end do
        end do
      end do
    end do
  end do

  if (.not. u_count == nf) then
    if (root_and_verbose) write (*,*) '* re-ordered memory layout'
#if defined(mpi2)
    call MPI_TYPE_CREATE_INDEXED_BLOCK(nf,blocklen,offset,MPICOMPLEX_X,          &
       &    TYPE_RW_FDISI,ierr)
#endif
  else
    if (root_and_verbose) write (*,*) '* default memory layout'
#if defined(mpi2)
    call MPI_TYPE_CONTIGUOUS(nf,MPICOMPLEX_X,TYPE_RW_FDISI,ierr)
#endif
  end if

#if defined(mpi2)
  call MPI_TYPE_COMMIT(TYPE_RW_FDISI,ierr)
#else
  ! not very good, but there are no obvious better choices
  TYPE_RW_FDISI = TYPE_OLD
#endif

  deallocate(offset)
  ! deallocate(blocklen)

end subroutine create_fdisi_io_type


end module restart
