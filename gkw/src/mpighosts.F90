!------------------------------------------------------------------------------
!> This module sets-up communication of MPI ghost cells
!> (formerly in exp_integration)
!>
!> The solution vector (distribution and fields) excluding ghost cells is 
!> stored at various stages of explicit integration in the arrays 
!> fdis(1:nsolc), fdis_tmp(1:nsolc) and fdis_tmp2(1:nsolc), 
!> which belong to the dist module.  Ghost cells are only ever stored in 
!> fdis_tmp(nsolc+1:msolc) and fdis_tmp2(nsolc+1:msolc); the index function
!> handles the access to these cells in non-optimised routines.
!> 
!> The ghost cells are only ever communicated WITHIN a single array 
!> (NOT between arrays) so the range (1:nsolc) must be updated before
!> the communication starts.  The MPI communications are initialised and 
!> grouped into sets in this module.  The derived datatypes used for these 
!> communications belong to the mpidatatypes module (and are commented 
!> there at declaration); these are initialised in the dist module.  
!------------------------------------------------------------------------------
module mpighosts

  use structures, only : ghostcomm

  implicit none

  private

  public :: persistent_comm_init, persistent_comm_end, mpistart, mpiwait

  !> mpi status and request info for all normal non-blocking communications

  !> In the comments below, "fields" refers only to NON gyro-averaged fields.
  !> "ga-fields" refers to gyro-averaged fields, which exist only for the non-spectral scheme.

  !> The first communication set is used in ALL schemes and includes
  !> all required non-x ghost cells, communicated in fdis_tmp:
  !>  * The distribution for all (non-x) ghost points and corner points
  !>  * The fields and ga-fields for s ghost points only
  type (ghostcomm), public, save :: gc1      

  !> Mirror of the distribution in velocity space, used only for the krook
  !> operator with parallel velocity parallelisation, communicated into fdis_nvpar
  type (ghostcomm), public, save ::  gckrook

  !> The remaining communications are used only for the non-spectral scheme with
  !> radial parallelisation.  The sub-divisions are needed to minimise the 
  !> communiciation at various stages, for an efficient non-spectral field solve.
  !> The x ghost cells may have > 2 points, used only for taking gyro-averages.  
  type (ghostcomm), public, save ::  gc1x2     !< f + ga-fields only (2 x gp) in fdis_tmp
  type (ghostcomm), public, save ::  gc1x_phi  !< fields only (all x gp) in fdis_tmp (for diags only)
  type (ghostcomm), public, save ::  gc2x_phi  !< fields only (all x gp) in fdis_tmp2 
  type (ghostcomm), public, save ::  gc2x_pga  !< ga-fields only (all x gp) in fdis_tmp2
  type (ghostcomm), public, save ::  gc2x2_pga !< ga-fields (2 x gp) in fdis_tmp2
  type (ghostcomm), public, save ::  gc2x2_f   !< f only (2 x gp)     in fdis_tmp2
    

contains

!-----------------------------------------------------------------------------
!> start a particular communication
!-----------------------------------------------------------------------------
subroutine mpistart(gcin)

 use general, only : gkw_abort
 
type (ghostcomm), intent(inout) :: gcin

integer :: ierr

#ifdef mpi2
    if (gcin%nreq > 0) then
      if (.not. gcin%running) then
        call MPI_STARTALL(gcin%nreq,gcin%preq,ierr)      
        gcin%running = .true.
      else
        call gkw_abort('Cannot start same communication twice')
      end if 
    end if
#else
   ! do nothing 
#endif

end subroutine mpistart

!-----------------------------------------------------------------------------
!> wait for a particular communication to finish
!-----------------------------------------------------------------------------
subroutine mpiwait(gcin)

#if defined(mpi2)
  use mpiinterface, only : MPI_STATUSES_IGNORE
#endif
  use general,      only : gkw_abort
  
type (ghostcomm), intent(inout) :: gcin
integer :: ierr

#ifdef mpi2
  if (gcin%nreq > 0) then
    if (gcin%running) then  
      call MPI_WAITALL(gcin%nreq,gcin%preq,MPI_STATUSES_IGNORE,ierr)
      gcin%running = .false.     
    else
      call gkw_abort('Cannot wait for a request not running')
    end if
  endif
#else
   ! do nothing 
#endif

end subroutine mpiwait


!-----------------------------------------------------------------------------
!> allocate for the persistent communication
!-----------------------------------------------------------------------------
subroutine persistent_comm_allocate

use grid,      only : parallel_s,parallel_vpar,lsendrecv_mu, lsendrecv_x,    &
                  & proc_mu_prev, proc_mu_next, proc_vpar_prev_mu_prev,      &
                  & proc_s_next,proc_s_prev,proc_vpar_next,proc_vpar_prev,   &
                  & proc_mu_next,proc_mu_prev,proc_vpar_next_s_next,         &
                  & proc_vpar_prev_s_next,proc_vpar_next_s_prev,             &
                  & proc_vpar_prev_s_prev,proc_vpar_next_mu_next,            &
                  & proc_vpar_prev_mu_next,proc_vpar_next_mu_prev,           &
                  & proc_x_next,proc_x_prev

use dist,      only : ghost_size_vpar_mu,ghost_size_vpar_s
use krook,     only : nlkrook, krook_option
use mpiinterface, only : MPI_PROC_NULL

  integer :: ierr

 ! setup the mpi communication information arrays
  gc1%nreq = 0
  gc1x2%nreq = 0
  gc1x_phi%nreq = 0
  gc2x_phi%nreq = 0
  gc2x_pga%nreq = 0
  gc2x2_pga%nreq = 0
  gc2x2_f%nreq = 0
  gckrook%nreq = 0

#ifdef mpi2
  ! fdisi_tmp ghost cells in vpar, s, mu
  if (parallel_s) then
    if (proc_s_next            /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_s_prev            /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
  end if  
  if (parallel_vpar) then
    if (proc_vpar_next         /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_vpar_prev         /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (nlkrook .and. ((krook_option ==1).or.(krook_option==4))) gckrook%nreq = 2
  end if
  if (lsendrecv_mu) then
    if (proc_mu_next           /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_mu_prev           /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
  end if
  if (ghost_size_vpar_s > 0) then
    if (proc_vpar_next_s_next  /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_vpar_prev_s_next  /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_vpar_next_s_prev  /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_vpar_prev_s_prev  /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
  end if
  if (ghost_size_vpar_mu > 0) then
    if (proc_vpar_next_mu_next /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_vpar_prev_mu_next /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_vpar_next_mu_prev /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
    if (proc_vpar_prev_mu_prev /= MPI_PROC_NULL) gc1%nreq = gc1%nreq + 2
  end if


  if (lsendrecv_x) then
    ! fdisi_tmp2 f ghost cells in x
    if (proc_x_next           /= MPI_PROC_NULL) gc2x2_f%nreq = gc2x2_f%nreq + 2
    if (proc_x_prev           /= MPI_PROC_NULL) gc2x2_f%nreq = gc2x2_f%nreq + 2
    
    gc1x2%nreq    = 2*gc2x2_f%nreq   ! fdisi_tmp f and phi_ga ghost cells in x
    gc1x_phi%nreq  = gc2x2_f%nreq   ! fdisi_tmp1 phi ghost cells in x
    gc2x_phi%nreq  = gc2x2_f%nreq   ! fdisi_tmp2 phi ghost cells in x
    gc2x_pga%nreq  = gc2x2_f%nreq   ! fdisi_tmp2 phi_ga ghost cells in x    
    gc2x2_pga%nreq = gc2x2_f%nreq   ! fdisi_tmp2 phi_ga ghost cells in x

     
  end if

#endif

  if (gc1%nreq > 0) allocate(gc1%preq(gc1%nreq),stat=ierr)
  if (gc1x2%nreq > 0) allocate(gc1x2%preq(gc1x2%nreq),stat=ierr)
  if (gc2x2_pga%nreq > 0) allocate(gc2x2_pga%preq(gc2x2_pga%nreq),stat=ierr)
  if (gc1x_phi%nreq > 0) allocate(gc1x_phi%preq(gc1x_phi%nreq),stat=ierr)
  if (gc2x_phi%nreq > 0) allocate(gc2x_phi%preq(gc2x_phi%nreq),stat=ierr)
  if (gc2x_pga%nreq > 0) allocate(gc2x_pga%preq(gc2x_pga%nreq),stat=ierr)
  if (gc2x2_f%nreq > 0) allocate(gc2x2_f%preq(gc2x2_f%nreq),stat=ierr)
  if (gckrook%nreq > 0) allocate(gckrook%preq(gckrook%nreq),stat=ierr) 

end subroutine persistent_comm_allocate

!-----------------------------------------------------------------------------
!> setup the persistent communication
!-----------------------------------------------------------------------------
subroutine persistent_comm_init()
 
  use grid, only : proc_vpar_next, proc_vpar_prev, proc_mu_prev, proc_mu_next,&
                     & proc_s_next, proc_s_prev,proc_vpar_prev_mu_prev,       &
                     & proc_vpar_prev_mu_next, proc_vpar_next_mu_next,        &
                     & proc_vpar_next_mu_prev,proc_vpar_next_s_next,          &
                     & proc_vpar_prev_s_next, proc_vpar_prev_s_prev,          &
                     & proc_vpar_next_s_prev, proc_x_next, proc_x_prev,       &
                     & proc_vpar_opposite
  use mpicomms, only : COMM_MU_NE, COMM_VPAR_NE, COMM_CART, COMM_X_NE
  use mpicomms, only : COMM_VPAR_NE, COMM_S_NE, COMM_MU_NE
  use dist,     only : ghost_size_vpar, ghost_size_s, ighost_sbp, &
                     & ighost_sbn,ighost_vparbp,ighost_vparbn,ighost_mubp,   &
                     & ighost_mubn, ighost_vparbp_mubp, ighost_vparbn_mubp,  &
                     & ighost_vparbp_mubn, ighost_vparbn_mubn,               &
                     & ghost_size_vpar_mu, ighost_vparbp_sbn, ghost_size_mu, &
                     & ighost_vparbp_sbp, ighost_vparbn_sbn,                 &
                     & ighost_vparbn_sbp, ghost_size_vpar_s, fdis_tmp,       &
                     & ighost_xbn, ighost_xbp, fdis_tmp2, ghost_size_x_phi,  &
                     & ghost_size_x_pga, ioffset_phi2, ioffset_phi_ga2, nf 
   use krook,        only : fdis_nvpar
   use general,      only : gkw_abort
   use mpiinterface, only : MPI_PROC_NULL, MPICOMPLEX_X, register_tag_range
   use mpicomms,     only : COMM_VPAR_NE
   use mpidatatypes, only : TYPE_PREV_S, TYPE_NEXT_S, TYPE_PREV_VPAR 
   use mpidatatypes, only : TYPE_NEXT_VPAR, TYPE_PREV_MU, TYPE_NEXT_MU
   use mpidatatypes, only : TYPE_PREV_VPAR_PREV_MU, TYPE_PREV_VPAR_NEXT_MU
   use mpidatatypes, only : TYPE_NEXT_VPAR_PREV_MU, TYPE_NEXT_VPAR_NEXT_MU
   use mpidatatypes, only : TYPE_PREV_S_PREV_VPAR, TYPE_PREV_S_NEXT_VPAR
   use mpidatatypes, only : TYPE_NEXT_S_PREV_VPAR, TYPE_NEXT_S_NEXT_VPAR
   use mpidatatypes, only : TYPE_NEXT_X2_F, TYPE_NEXT_X2_F_RECV
   use mpidatatypes, only : TYPE_NEXT_X2_PGA, TYPE_NEXT_X2_PGA_RECV
   use mpidatatypes, only : TYPE_NEXT_X_PGA, TYPE_PREV_X2_F
   use mpidatatypes, only : TYPE_PREV_X2_F_RECV, TYPE_PREV_X2_PGA
   use mpidatatypes, only : TYPE_PREV_X2_PGA_RECV, TYPE_PREV_X_PGA
   use mpidatatypes, only : TYPE_NEXT_X_PHI, TYPE_PREV_X_PHI
   
  integer :: itag, ierr

  call persistent_comm_allocate

  !
  ! communication in the s direction (gc1)
  !
  call mpiinit_pair(gc1,fdis_tmp,COMM_S_NE,                                &
                  & ghost_size_s, ighost_sbp + 1, ighost_sbn + 1,          &
                  & 1, 1, TYPE_PREV_S, TYPE_NEXT_S, PROC_S_PREV, PROC_S_NEXT)
  
  !
  ! communication in the parallel velocity direction (gc1)
  !
  call mpiinit_pair(gc1,fdis_tmp,COMM_VPAR_NE,                             &
                  & ghost_size_vpar, ighost_vparbp + 1, ighost_vparbn + 1, &
                  & 1, 1, TYPE_PREV_VPAR, TYPE_NEXT_VPAR,                  &
                  &       PROC_VPAR_PREV, PROC_VPAR_NEXT                   )
  
  !
  ! Communicate with procs along the mu direction; used for collisions (gc1)
  !
  call mpiinit_pair(gc1,fdis_tmp,COMM_MU_NE,                               &
                  & ghost_size_mu, ighost_mubp + 1, ighost_mubn + 1,       &
                  & 1, 1, TYPE_PREV_MU, TYPE_NEXT_MU,                      &
                  &       PROC_MU_PREV, PROC_MU_NEXT                       )

  !
  ! These are only used when running with collisions, where a point is
  ! required from a processor shifted in both the vpar and mu directions (when
  ! parallelising over both vpar and mu). (gc1)
  !
   call mpiinit_pair(gc1,fdis_tmp,COMM_CART,                               &
                  & ghost_size_vpar_mu,                                    &
                  & ighost_vparbp_mubp + 1, ighost_vparbn_mubn + 1,        &
                  & 1, 1, TYPE_PREV_VPAR_PREV_MU, TYPE_NEXT_VPAR_NEXT_MU,  &
                  &       proc_vpar_prev_mu_prev, proc_vpar_next_mu_next   )
                  
   call mpiinit_pair(gc1,fdis_tmp,COMM_CART,                               &
                  & ghost_size_vpar_mu,                                    &
                  & ighost_vparbp_mubn + 1, ighost_vparbn_mubp + 1,        &
                  & 1, 1, TYPE_PREV_VPAR_NEXT_MU, TYPE_NEXT_VPAR_PREV_MU,  &
                  &       proc_vpar_prev_mu_next,  proc_vpar_next_mu_prev  )
  
  !
  ! These will only be enabled if Arakawa type finite-differencing is used,
  ! where a point is required on a processor shifted in both vpar and s
  ! directions (when parallelising over both vpar and s). (gc1)
  !
   call mpiinit_pair(gc1,fdis_tmp,COMM_CART,                             &
                  & ghost_size_vpar_s,                                   &
                  & ighost_vparbp_sbp + 1, ighost_vparbn_sbn + 1,        &
                  & 1, 1, TYPE_PREV_S_PREV_VPAR, TYPE_NEXT_S_NEXT_VPAR,  &
                  &       proc_vpar_prev_s_prev, proc_vpar_next_s_next   )
                  
   call mpiinit_pair(gc1,fdis_tmp,COMM_CART,                             &
                  & ghost_size_vpar_s,                                   &
                  & ighost_vparbp_sbn + 1, ighost_vparbn_sbp + 1,        &
                  & 1, 1, TYPE_NEXT_S_PREV_VPAR, TYPE_PREV_S_NEXT_VPAR,  &
                  &       proc_vpar_prev_s_next,  proc_vpar_next_s_prev  )
                  
  
  ! All following communicator sets are only used for the nonspectral 
  ! version of the code with parallel_x
  
  ! 
  ! communication of f (2 gp) in the x direction into fdis_tmp (gc1x2) 
  ! THIS DIFFERS FROM THE OTHERS BY USING THE RECV DATATYPE OPTIONAL ARGUMENTS
  call mpiinit_pair(gc1x2,fdis_tmp,COMM_X_NE,                              &
                  & 1, 1, 1,                                               &
                  & 1, 1, TYPE_PREV_X2_F, TYPE_NEXT_X2_F,                  &
                  &       PROC_X_PREV,   PROC_X_NEXT,                      &
                  &       TYPE_PREV_X2_F_RECV, TYPE_NEXT_X2_F_RECV         ) 
  
  ! communication of phi_ga (2 gp) in the x direction into fdis_tmp (gc1x2) 
  ! THIS DIFFERS FROM THE OTHERS BY USING THE RECV DATATYPE OPTIONAL ARGUMENTS
  call mpiinit_pair(gc1x2,fdis_tmp,COMM_X_NE,                              &
                  & 1, 1, 1,                                               &
                  & 1, 1, TYPE_PREV_X2_PGA, TYPE_NEXT_X2_PGA,              &
                  &       PROC_X_PREV,   PROC_X_NEXT,                      &
                  &       TYPE_PREV_X2_PGA_RECV, TYPE_NEXT_X2_PGA_RECV     ) 
    
  if (allocated(fdis_tmp2)) then 
     
  ! 
  ! communication of f (2 gp) only in the x direction into fdis_tmp2 (gc2x2_f) 
  ! THIS DIFFERS FROM THE OTHERS BY USING THE RECV DATATYPE OPTIONAL ARGUMENTS
    call mpiinit_pair(gc2x2_f,fdis_tmp2,COMM_X_NE,                         &
                  & 1, 1, 1,                                               &
                  & 1, 1, TYPE_PREV_X2_F, TYPE_NEXT_X2_F,                  &
                  &       PROC_X_PREV,   PROC_X_NEXT,                      &
                  &       TYPE_PREV_X2_F_RECV, TYPE_NEXT_X2_F_RECV         )  
                  
  ! 
  ! communication of fields only in the x direction in fdis_tmp2 (gc2x_phi) 
  !
    call mpiinit_pair(gc2x_phi,fdis_tmp2,COMM_X_NE,                         &
                  & ghost_size_x_phi, ioffset_phi2 + ighost_xbp + 1,        &
                                      ioffset_phi2 + ighost_xbn + 1,        &
                  & 1, 1, TYPE_PREV_X_PHI, TYPE_NEXT_X_PHI,                 &
                  &       PROC_X_PREV,     PROC_X_NEXT                      )
                  
  ! 
  ! communication of fields only in the x direction in *fdis_tmp* (gc1x_phi) 
  ! used by diagnostics only:  fdis_tmp is always allocated if fdis_tmp2 is allocated
  !
    call mpiinit_pair(gc1x_phi,fdis_tmp,COMM_X_NE,                          &
                  & ghost_size_x_phi, ioffset_phi2 + ighost_xbp + 1,        &
                                      ioffset_phi2 + ighost_xbn + 1,        &
                  & 1, 1, TYPE_PREV_X_PHI, TYPE_NEXT_X_PHI,                 &
                  &       PROC_X_PREV,     PROC_X_NEXT                      )
  
  ! 
  ! communication of ga fields only in x direction in fdis_tmp2 (gc2x_pga) 
  !                  
    call mpiinit_pair(gc2x_pga,fdis_tmp2,COMM_X_NE,                         &
                  & ghost_size_x_pga, ioffset_phi_ga2 + ighost_xbp + 1,     &
                  &                   ioffset_phi_ga2 + ighost_xbn + 1,     &
                  & 1, 1, TYPE_PREV_X_PGA, TYPE_NEXT_X_PGA,                 &
                  &       PROC_X_PREV,     PROC_X_NEXT                      )

  ! 
  ! communication of ga fields (2 gp) only in x in fdis_tmp2 (gc2x2_pga) 
  ! THIS DIFFERS FROM THE OTHERS BY USING THE RECV DATATYPE OPTIONAL ARGUMENTS
    call mpiinit_pair(gc2x2_pga,fdis_tmp2,COMM_X_NE,                       &
                  & 1, 1, 1,                                               &
                  & 1, 1, TYPE_PREV_X2_PGA, TYPE_NEXT_X2_PGA,              &
                  &       PROC_X_PREV,   PROC_X_NEXT,                      &
                  &       TYPE_PREV_X2_PGA_RECV, TYPE_NEXT_X2_PGA_RECV     )


   end if !allocated fdis_tmp2
#ifdef mpi2   
   if (allocated(fdis_nvpar).and. allocated(fdis_tmp)) then
     !
     ! mirror communication in the parallel velocity direction (gckrook)
     !                  
     itag = get_new_tag()

     if (proc_vpar_opposite == MPI_PROC_NULL) then 
       call gkw_abort('Severe error in mpighosts: gckrook')
     end if  
       
     gckrook%irequest = gckrook%irequest + 1
     call MPI_SEND_INIT(                                               &
                     &          fdis_tmp(1),             nf,           &
                     &          MPICOMPLEX_X,     proc_vpar_opposite,  &
                     &          itag,             COMM_VPAR_NE,        &
                     &          gckrook%preq(gckrook%irequest),        ierr)

     gckrook%irequest = gckrook%irequest + 1
     call MPI_RECV_INIT(                                               &
                    &         fdis_nvpar(1),           nf,             &
                    &         MPICOMPLEX_X,     proc_vpar_opposite,    &
                    &         itag,             COMM_VPAR_NE,          &
                    &         gckrook%preq(gckrook%irequest),         ierr)
  end if
#endif 

  ! check that the requests are allocated the right size
  call check_ghostcomm(gc1)
  call check_ghostcomm(gc1x2)
  call check_ghostcomm(gc1x_phi)
  call check_ghostcomm(gc2x_phi)
  call check_ghostcomm(gc2x_pga)
  call check_ghostcomm(gc2x2_f)
  call check_ghostcomm(gc2x2_pga)
  call check_ghostcomm(gckrook)  
                    
  ! TOP level ROUTINE END

  contains 
  
  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine mpiinit_pair(gcin, fdis_dum, comm,                             &
                            recv_size, recv_loc_p, recv_loc_n,              &
                            send_loc, send_size,                            &
                            type_prev, type_next,  proc_prev, proc_next,    &
                            type_prev_recv_op, type_next_recv_op            )

  use dist, only  : msolc

  type (ghostcomm), intent(inout) :: gcin  
  integer, intent(in) :: recv_loc_p, recv_loc_n, recv_size 
  integer, intent(in) :: send_loc, send_size
  integer, intent(in) :: type_prev, type_next, proc_prev, proc_next
  integer, intent(in) :: comm
  integer, intent(in), optional :: type_prev_recv_op, type_next_recv_op
  complex, intent(inout) :: fdis_dum(msolc)
  
  integer :: tag_next, tag_prev, ierr, type_prev_recv, type_next_recv 
  
  tag_next = get_new_tag()
  tag_prev = get_new_tag()
  
  if (recv_size < 1) return

  type_next_recv = MPICOMPLEX_X
  type_prev_recv = MPICOMPLEX_X

  if (present(type_prev_recv_op)) type_prev_recv = type_prev_recv_op   
  if (present(type_next_recv_op)) type_next_recv = type_next_recv_op
  
#ifdef mpi2    
  if (proc_prev /= MPI_PROC_NULL) then
    gcin%irequest = gcin%irequest + 1
    call MPI_RECV_INIT(                                                      &
                   &         fdis_dum(recv_loc_p),           recv_size,      &
                   &                   type_prev_recv,       proc_prev,      &
                   &                   tag_next,             COMM,           & 
                   &                   gcin%preq(gcin%irequest),         ierr)
  end if
  if (proc_next /= MPI_PROC_NULL) then
    gcin%irequest = gcin%irequest + 1
    call MPI_RECV_INIT(                                                      &
                   &         fdis_dum(recv_loc_n),           recv_size,      &
                   &                   type_next_recv,       proc_next,      &
                   &                   tag_prev,             COMM,           &
                   &                   gcin%preq(gcin%irequest),         ierr)
  end if
  if (proc_prev /= MPI_PROC_NULL) then
    gcin%irequest = gcin%irequest + 1
    call MPI_SEND_INIT(                                                      &
                   &          fdis_dum(send_loc),             send_size,     &
                   &                    type_prev,            proc_prev,     &
                   &                    tag_prev,             COMM,          &
                   &                    gcin%preq(gcin%irequest),        ierr)
  end if
  if (proc_next /= MPI_PROC_NULL) then
    gcin%irequest = gcin%irequest + 1
    call MPI_SEND_INIT(                                                      &
                   &          fdis_dum(send_loc),             send_size,     &
                   &                    type_next,            proc_next,     &
                   &                    tag_next,             COMM,          &
                   &                    gcin%preq(gcin%irequest),        ierr)
  end if
  
#endif
 
  end subroutine mpiinit_pair
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !> return a tag; this must be called collectively by all processes, and in
  !> the same order, for the purpose of the parent module/routine this is
  !> fullfilled.
  !----------------------------------------------------------------------------
  function get_new_tag()
    use mpiinterface, only : register_tag_range
    integer, save :: tag = 0
    integer, save :: mpi_tag_end_inkl = 0
    integer :: get_new_tag

    if(tag == mpi_tag_end_inkl) then
      ! register tags in chunks of 10
      call register_tag_range(10,tag, &
         & mpi_tag_end_inkl)
    end if
    
    tag = tag + 1  
    get_new_tag = tag 
  
  end function get_new_tag
  !-------------------------------------
  
  !-------------------------------------
  subroutine check_ghostcomm(gcin)
  
    type (ghostcomm), intent(inout) :: gcin  
  
    if (gcin%irequest < gcin%nreq) then
      call gkw_abort('a ghostcomm structure is allocated too large') 
    else if (gcin%irequest > gcin%nreq) then
      call gkw_abort('a ghostcomm structure is allocated too small') 
    end if  
    gcin%initialised = .true.

  end subroutine check_ghostcomm
  !-------------------------------------
  
end subroutine persistent_comm_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!> end persistent communications
!-----------------------------------------------------------------------
subroutine persistent_comm_end()

  call deallocate_free(gc1)
  call deallocate_free(gc1x2)
  call deallocate_free(gc1x_phi)
  call deallocate_free(gc2x_phi)
  call deallocate_free(gc2x_pga)
  call deallocate_free(gc2x2_f)
  call deallocate_free(gc2x2_pga)
  call deallocate_free(gckrook)
  
  contains
  
  !--------------------------------------------------------
  subroutine deallocate_free(gcin)
  
    type (ghostcomm), intent(inout) :: gcin  
    
      integer :: ierr, i
      
      gcin%initialised=.false.
      
#ifdef mpi2
     if (gcin%nreq > 0) then
      do i = 1, gcin%nreq
        call MPI_REQUEST_FREE(gcin%preq(i),ierr)
      end do
     end if
#endif

     ! should also free (most) mpi datatypes ?

     if (associated(gcin%preq))   deallocate(gcin%preq)
     !if (allocated(gcin%preq))   deallocate(gcin%preq)
     !if (allocated(gcin%pstat))  deallocate(gcin%pstat)
  
  end subroutine deallocate_free 
  !---------------------------------------------------------

end subroutine persistent_comm_end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

end module mpighosts
