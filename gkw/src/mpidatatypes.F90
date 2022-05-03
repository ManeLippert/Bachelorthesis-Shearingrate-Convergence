!-----------------------------------------------------------------------------
!> provide derived datatypes for communication of processor boundaries and
!> parallel io etc.
!> 
!> Derived datatypes for ghost cell communciation are initialised in dist.
!> The communications using these datatypes are managed by module mpighosts.
!-----------------------------------------------------------------------------

module mpidatatypes
  
  implicit none

  private

  ! --- Public Routines ---
  public :: create_subarray_datatype,free_datatype
  
  ! --- Public --- !

  ! In the comments below, "fields" refers only to NON gyro-averaged fields.
  ! All ghost cells have two points, used for finite difference derivatives.

  ! Non-spectral scheme only:
  ! "ga-fields" refers to gyro-averaged fields, which exist only for the non-spectral scheme.
  ! The x ghost cells may have > 2 points, used only for taking gyro-averages.

  integer, public, save :: TYPE_NEXT_MU     !< dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_PREV_MU     !< dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_NEXT_VPAR   !< dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_PREV_VPAR   !< dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_NEXT_S      !< dist., fields, and ga-fields.  Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_PREV_S      !< dist., fields, and ga-fields.  Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_NEXT_X_PHI  !< fields only (all gp). Communicated in fdis_tmp2 using gc2x_phi
  integer, public, save :: TYPE_PREV_X_PHI  !< fields only (all gp). Communicated in fdis_tmp2 using gc2x_phi
  integer, public, save :: TYPE_NEXT_X_PGA  !< ga-fields only (all gp). Communicated in fdis_tmp2 using gc2x_pga
  integer, public, save :: TYPE_PREV_X_PGA  !< ga-fields only (all gp). Communicated in fdis_tmp2 using gc2x_pga
  integer, public, save :: TYPE_NEXT_X2_F    !< dist. (2 gp) only. Communicated in fdis_tmp using gc1x2 and in fdis_tmp2 using gc2x2_f
  integer, public, save :: TYPE_PREV_X2_F    !< dist. (2 gp) only. Communicated in fdis_tmp using gc1x2 and in fdis_tmp2 using gc2x2_f
  integer, public, save :: TYPE_NEXT_X2_PGA  !< ga-fields (2 gp) only. Communicated in fdis_tmp2 using gc2x2_pga
  integer, public, save :: TYPE_PREV_X2_PGA  !< ga-fields (2 gp) only. Communicated in fdis_tmp2 using gc2x2_pga
  integer, public, save :: TYPE_NEXT_X2_PGA_RECV  !< recv version of the above. Communicated in fdis_tmp2 using gc2x2_pga
  integer, public, save :: TYPE_PREV_X2_PGA_RECV  !< recv version of the above. Communicated in fdis_tmp2 using gc2x2_pga
  integer, public, save :: TYPE_NEXT_S_NEXT_VPAR  !< Arakawa corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_NEXT_S_PREV_VPAR  !< Arakawa corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_PREV_S_NEXT_VPAR  !< Arakawa corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_PREV_S_PREV_VPAR  !< Arakawa corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_NEXT_VPAR_NEXT_MU  !< Collision operator corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_PREV_VPAR_NEXT_MU  !< Collision operator corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_NEXT_VPAR_PREV_MU  !< Collision operator corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_PREV_VPAR_PREV_MU  !< Collision operator corner points, dist. only. Communicated in fdis_tmp using gc1
  integer, public, save :: TYPE_MOD_SEND   !< gather over toroidal modes in nonspectral field solve
  integer, public, save :: TYPE_MOD_RECV   !< gather over toroidal modes in nonspectral field solve

  ! the following recv datatypes are used as temporary workaround until
  ! the index function can handle different numbers of ghost cells for
  ! the fields, at that point they could be removed.
  integer, public, save :: TYPE_NEXT_X2_F_RECV  !< dist. only. Communicated in fdis_tmp using gc1x2 and in fdis_tmp2 using gc2x2_f
  integer, public, save :: TYPE_PREV_X2_F_RECV  !< dist. only. Communicated in fdis_tmp using gc1x2 and in fdis_tmp2 using gc2x2_f
  
contains


!----------------------------------------------------------------------------
!> Create a new sub-array datatype given an existing atom datatype and up to 6
!> identifiers for the various dimensions (found in the global module, all
!> prefixed `id_').
!> Optionally, x and y gridsizes may be specified. FIXME
!>
!> Optionally, global species-, y- and x-size deviating from current global
!> grid sizes can be specified. 
!>
!> Typically, the new datatype is used to write a local array
!> through a `slot' into a global array on file, or to read local arrays
!> from a global array on file.
!----------------------------------------------------------------------------  
subroutine create_subarray_datatype(atom_type, new_type, &
   & i1, i2, i3, i4, i5, i6, global_and_local_xsize, global_and_local_ysize, &
  & global_and_local_spsize, global_and_local_ssize, global_spsize, &
  & global_ysize, global_xsize)

#if defined(mpi2)
  use mpiinterface, only : MPI_ORDER_FORTRAN
#else 
  use global,       only : I_RUN_WITHOUT_MPI
#endif
  use general,      only : gkw_abort
  use global,       only : id_s,id_sp,id_vpar,id_mu,id_x,id_mod,id_dummy
  use grid,         only : n_x_grid,nmod,n_s_grid,ns,n_mu_grid,nmu,n_vpar_grid
  use grid,         only : nx,nvpar,number_of_species,nsp,ls,lvpar,lmu,lsp,lrx
  use grid,         only : n_procs_x, n_procs_sp
  use mpicomms,     only : COMM_SP_EQ_procset, COMM_X_EQ_procset

  integer, intent(in) :: atom_type
  integer, intent(out) :: new_type
  integer, intent(in) :: i1
  integer, optional, intent(in) :: i2,i3,i4,i5,i6
  !> real space array sizes
  integer, optional, intent(in) :: global_and_local_xsize
  integer, optional, intent(in) :: global_and_local_ysize
  integer, optional, intent(in) :: global_and_local_spsize
  integer, optional, intent(in) :: global_and_local_ssize
  !> real space global array size
  integer, optional, intent(in) :: global_spsize 
  integer, optional, intent(in) :: global_ysize
  integer, optional, intent(in) :: global_xsize
 
  integer :: ndims,i,ierr
  integer :: ng_x, nl_x, start_x
  integer :: ng_y
  integer :: ng_sp, nl_sp, start_sp
  integer :: ng_s, nl_s, start_s
  integer, parameter :: maxdims=6
  !> an array to hold the IDs of the chosen dimensions
  integer, dimension(maxdims) :: id_in
  !> global array shape and shape of the local subarray
  integer, dimension(maxdims) :: gshape, lshape
  !> start point in the respective dimension
  integer, dimension(maxdims) :: starts
  logical :: ldum
  
  nl_x = nx
  ng_x = n_x_grid
  start_x = -lrx(0)
  if(present(global_and_local_xsize)) then
    ! this is intended for spectral_radius=T runs where there is no
    ! parallelisation in x-direction.
    !if (.not. spectral_radius) call gkw_abort('No radial fft in subarray_datatype')
    ng_x = global_and_local_xsize
    nl_x = global_and_local_xsize
    start_x = 0
  end if 
  if(present(global_xsize)) then
    ! this is intended for restarted runs, where the gridsize in x has been
    ! changed compared to the old grid size.
    ng_x = global_xsize
    ! the local grid size derived from to the global_xsize passed as function 
    ! argument together with the current processor setup
    nl_x = global_xsize / n_procs_x
    ! the start index of the current processor subset with respect to the 
    ! global x-grid is determined by explicitely using the local 'color'
    ! of the communication split
    start_x = nl_x * (COMM_X_EQ_procset - 1)
  end if 

  ng_y = nmod
  if(present(global_and_local_ysize)) ng_y = global_and_local_ysize
  if(present(global_ysize)) then
    ng_y = global_ysize
  end if

  nl_sp = nsp
  ng_sp = number_of_species
  start_sp = -lsp(0)
  if(present(global_and_local_spsize)) then
    nl_sp = global_and_local_spsize
    ng_sp = global_and_local_spsize
  end if
  if(present(global_spsize)) then
    ! this is intended for restarted runs, where the number of species has been
    ! changed compared to the old number.
    ng_sp = global_spsize
    ! the local grid size derived from the global_spsize passed as function 
    ! argument together with the current processor setup
    nl_sp = global_spsize / n_procs_sp
    ! the start index of the current processor subset with respect to the 
    ! global species grid is determined by explicitely using the local 'color'
    ! of the communication split
    start_sp = nl_sp * (COMM_SP_EQ_procset - 1)
  end if 

  nl_s = ns
  ng_s = n_s_grid
  start_s = -ls(0)
  if(present(global_and_local_ssize)) then
    nl_s = global_and_local_ssize
    ng_s = global_and_local_ssize
    start_s = 0
  end if

  ! put the given IDs into an array
  id_in(:) = id_dummy
  id_in(1) = i1
  if (present(i2)) id_in(2) = i2
  if (present(i3)) id_in(3) = i3
  if (present(i4)) id_in(4) = i4
  if (present(i5)) id_in(5) = i5
  if (present(i6)) id_in(6) = i6
  ! now id_in holds valid indices up to a certain element, and
  ! id_dummy behind that.

  ndims = 0
  get_dims : do i = 1, maxdims
  
    ldum = .false.
    select case (id_in(i))
      !  The starts(i) expression is a slight hack to return start
      !  point = ispb - 1
      case(id_s)
        gshape(i)=ng_s;              lshape(i)=nl_s;  starts(i)=start_s
      case(id_vpar)
        gshape(i)=n_vpar_grid;       lshape(i)=nvpar; starts(i)=-lvpar(0)
      case(id_mu)
        gshape(i)=n_mu_grid;         lshape(i)=nmu;   starts(i)=-lmu(0)
      case(id_x)
        gshape(i)=ng_x;              lshape(i)=nl_x;  starts(i)=start_x
      case(id_mod)
        gshape(i)=ng_y;              lshape(i)=ng_y;  starts(i)=0
      case(id_sp)
        gshape(i)=ng_sp;             lshape(i)=nl_sp; starts(i)=start_sp
      case(id_dummy)
        ldum = .true.
      case default
        call gkw_abort('bad id in mpidatatypes')
    end select
    
    if (ldum) exit get_dims
    ndims=ndims+1

  end do get_dims

  ! create the subarray datatype
#if defined(mpi2)
  ierr=0
  call mpi_type_create_subarray(ndims,gshape,lshape,starts, &
     & MPI_ORDER_FORTRAN,atom_type, &
     & new_type,ierr)
  call mpi_type_commit(new_type,ierr)
#else
  
  ! safe defaults
  new_type = I_RUN_WITHOUT_MPI
  ierr = 0
  
#endif

end subroutine create_subarray_datatype


!****************************************************************************
!> Free MPI datatypes no longer needed, perhaps to reduce MPI overheads.
!----------------------------------------------------------------------------
subroutine free_datatype(old_type)

  integer, intent(in) :: old_type

  integer :: ierr

#if defined(mpi2)
  call MPI_TYPE_FREE(old_type,ierr)
#endif

  continue

end subroutine free_datatype


end module mpidatatypes
