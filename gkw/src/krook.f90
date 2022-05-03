!-----------------------------------------------------------------------------
!> Module containing the source connected with a global run 
!-----------------------------------------------------------------------------
module krook 

  implicit none 

  private 

  public :: krook_read_nml, krook_bcast_nml 
  public :: krook_write_nml, add_krook, krook_init 
  public :: nlkrook, nlbound, bwidth, gammab, krook_option, fdis_nvpar 

  !> Logical that determines if the Krook operator is used 
  logical, save :: nlkrook 

  !> Damping rate used in the Krook operator 
  real, save :: gammak 
  
  !> For krook option 4, prefactor of n=0 mode damping rate.  I.e
  !> can selectively chose to damp the zonal mode stronger than the 
  !> turbulence
  real, save :: gamkpre

  !> Logical that determines in damping is applied in the boundary layers 
  logical, save :: nlbound 

  !> width of the boundary layers (in grid points) 
  real, save :: bwidth 

  !> damping rate of the boundary layers 
  real, save :: gammab 

  !> option for the various forms of the krook operator 
  integer, save :: krook_option 

  interface krook_write_nml
    module procedure krook_read_nml
  end interface

  ! arrays for the krook operator 
  integer, save              :: nmats 
  integer, save, allocatable :: iss(:), jss(:), id1(:), id2(:)
  real,    save, allocatable :: mats1(:), mats2(:), mats3(:), mats4(:)
  real,    save, allocatable :: pref(:), prefgam(:)
  complex, save, allocatable :: sss(:), sst(:), fdis_nvpar(:) !, sss_tmp(:)

contains 

!------------------------------------------------------------------------------
!> Subroutine that reads the namelist 
!------------------------------------------------------------------------------
subroutine krook_read_nml(lun,io_stat,lwrite)
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor

  integer, intent(in)  :: lun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  namelist / krook / nlkrook, gammak, nlbound, bwidth, gammab, krook_option, gamkpre

  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then

      ! default values 
      nlkrook = .false. 
      gammak = 0.0 ! please check if this default value makes sense
                   ! and remove comment
      nlbound = .false.
      bwidth = 0.0 ! please check if this default value makes sense
                   ! and remove comment
      gammab = 0.0 ! please check if this default value makes sense
                   ! and remove comment
      krook_option = 1
      gamkpre = 1.0E0

      read(lun,NML=krook,IOSTAT=io_stat)
    else
      ! do nothing
    end if
  else
    if(root_processor) write(lun,NML=krook)

    call write_run_parameter('krook', 'nlkrook', nlkrook)
    call write_run_parameter('krook', 'gammak', gammak)
    call write_run_parameter('krook', 'nlbound', nlbound)
    call write_run_parameter('krook', 'bwidth', bwidth)
    call write_run_parameter('krook', 'gamkpre', gamkpre)
    call write_run_parameter('krook', 'gammab', gammab)
    call write_run_parameter('krook', 'krook_option', krook_option)

  end if

end subroutine krook_read_nml

!------------------------------------------------------------------------------
!> bcast the source terms namelist params (and check)
!------------------------------------------------------------------------------
subroutine krook_bcast_nml

  use mpiinterface, only : mpibcast
  use general,      only : gkw_exit
  use control,      only : spectral_radius
  
  call mpibcast(nlkrook,      1)
  call mpibcast(gammak,       1)
  call mpibcast(nlbound,      1)
  call mpibcast(bwidth,       1)
  call mpibcast(gamkpre,      1)
  call mpibcast(gammab,       1)
  call mpibcast(krook_option, 1)
  
  if ((nlkrook .or. nlbound).and.spectral_radius) then
    call gkw_exit('Use krook / boundary operator with spectral_radius = F')
  end if   

end subroutine krook_bcast_nml


!------------------------------------------------------------------------------
!> Set up the arrays needed for the calculation in krook. Note that this is 
!> really a field and should likely be moved into the matrix as such. First 
!> we will however estabilish its usefulness
!> NOTE TOO that this field is very closely related to the gyro-centre density 
!> that is anyway calculated to obtain the potential. It is likely a good idea
!> to combine the two and too save computational time.  
!------------------------------------------------------------------------------
subroutine krook_init

  use grid,           only : nmod, nx, ns, nvpar, nmu, nsp, parallel_vpar 
  use dist,           only : ifdis, nsolc, fmaxwl, nf 
  use index_function, only : indx 
  use velocitygrid,   only : intvp, intmu, vpgr
  use geom,           only : bn 
  use components,     only : de, dgrid 
  use general,        only : gkw_abort 
  use matdat,         only : compress_piece
  use mode,           only : iyzero

  integer itel, itel2, imod, ix, i, j, k, is, ierr 
  integer, allocatable :: iid(:), jjd(:)  ! dummies to allow repeated sort

  real :: int_max
  
  if (.not. nlkrook) return

  ! allocate the arrays 
  ierr = 0 
  allocate(iss(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate iss in source')

  allocate(jss(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate jss in source')
  
  allocate(iid(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate iid in source')

  allocate(jjd(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate jjd in source')  

  allocate(mats1(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate mats1 in source')

  allocate(mats2(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate mats2 in source')

  allocate(mats3(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate mats3 in source')

  allocate(id1(nsolc), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate id1 in source')

  !Pref is the prefactor that weights the symmetrisation and prefgam is a
  !prefactor to increase or deacrease the damping rate for the n=0 mode.
  !i.e this allows one to damp away the equilibrium parallel flow
  if (krook_option == 4)then
    allocate(pref(nmod*nx*ns*nsp), stat = ierr) 
    if (ierr.ne.0) call gkw_abort('Could not allocate pref in source')
    allocate(prefgam(nmod*nx*ns*nsp), stat = ierr) 
    if (ierr.ne.0) call gkw_abort('Could not allocate pref in source')
  endif
 
  if ((krook_option == 1).or.(krook_option == 3).or.(krook_option == 4)) then 
    allocate(id2(nsolc), stat = ierr) 
    if (ierr.ne.0) call gkw_abort('Could not allocate is in source')
    
    if ((krook_option == 1).or.(krook_option == 4)) then
      if (parallel_vpar) allocate(fdis_nvpar(nf), stat = ierr)
      if (ierr.ne.0) call gkw_abort('Could not allocate fdis_nvpar in source')
    endif

  endif 
      
  allocate(sss(nmod*nx*ns*nsp), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Could not allocate sss in source')

  !allocate(sss_tmp(nmod*nx*ns*nsp), stat = ierr) 
  !if (ierr.ne.0) call gkw_abort('Could not allocate sss_tmp in source')

  if (krook_option == 2) then 

    allocate(mats4(nsolc), stat = ierr) 
    if (ierr.ne.0) call gkw_abort('Could not allocate mats4 in source')

    allocate(sst(nmod*nx*ns*nsp), stat = ierr) 
    if (ierr.ne.0) call gkw_abort('Could not allocate sst in source')

  endif 

  if ((krook_option==1).or.(krook_option==3).or.(krook_option==4)) then    
   
    itel = 0 
    itel2 = 0
    do imod = 1, nmod; do ix = 1, nx; do i = 1, ns; do is = 1, nsp 
    
       !This is here as these factors are only spatially varying
       itel2 = itel2 + 1
       if(krook_option==4)then
          if(imod==iyzero)then
            pref(indx_sss(imod,ix,i,is)) = 1.E0
            prefgam(indx_sss(imod,ix,i,is)) = gamkpre
          else
            pref(indx_sss(imod,ix,i,is)) = 0.5E0
            prefgam(indx_sss(imod,ix,i,is)) = 1.0E0 
          endif        
      endif

      do j = 1, nmu; do k = 1, nvpar
        itel = itel + 1 
        iss(itel) = indx_sss(imod,ix,i,is)
        jss(itel) = indx(ifdis,imod,ix,i,j,k,is) 
        mats1(itel) = intvp(i,j,k,is)*intmu(j)*bn(ix,i) 
        mats2(itel) = fmaxwl(ix,i,j,k,is) * dgrid(is) / de(ix,is)  
        id1(itel) = indx(ifdis,imod,ix,i,j,k,is) 
        id2(itel) = indx(ifdis,imod,ix,i,j,nvpar-k+1,is)  

      end do; end do 

    end do; end do; end do ; end do

    ! store the number of elements 
    nmats = itel 
    
    ! sort the matrices for basic cache efficiency, reset row and columns except last time
    ! since no sort is implemented for integer matrices, use mats3 as a dummy
    ! saves 25% of time in krook routine (krook routine is ~5% of total)
    iid=iss; jjd=jss; nmats = itel 
    call compress_piece(1,nmats,iid,jjd,mats1)
    iid=iss; jjd=jss; nmats = itel 
    call compress_piece(1,nmats,iid,jjd,mats2)
    iid=iss; jjd=jss; nmats = itel; mats3=id1 
    call compress_piece(1,nmats,iid,jjd,mats3)    
    id1(1:nmats) = nint(mats3(1:nmats));
    mats3(1:nmats) = real(id2(1:nmats));
    call compress_piece(1,nmats,iss,jss,mats3)     
    id2(1:nmats) = nint(mats3(1:nmats))
    deallocate(mats3)

  endif  
 
  if (krook_option ==2) then 

    itel = 0 
    do imod = 1, nmod; do ix = 1, nx; do i = 1, ns 
    
      do j = 1, nmu; do is = 1, nsp  

        int_max = 0. 
        do k = 1, nvpar
          int_max = int_max + intvp(i,j,k,is)*fmaxwl(ix,i,j,k,is)* & 
                  & vpgr(i,j,k,is)**2
        end do  

        do k = 1, nvpar
          itel = itel + 1 
          iss(itel) = indx_sss(imod,ix,i,is)
          jss(itel) = indx(ifdis,imod,ix,i,j,k,is) 
          mats1(itel) = intvp(i,j,k,is)*intmu(j)*bn(ix,i) 
          mats3(itel) = intvp(i,j,k,is)*vpgr(i,j,k,is) 
          mats2(itel) = fmaxwl(ix,i,j,k,is) * dgrid(is) / de(ix,is)  
          mats4(itel) = vpgr(i,j,k,is)*fmaxwl(ix,i,j,k,is)/int_max
          id1(itel) = indx(ifdis,imod,ix,i,j,k,is) 
        end do
      end do; end do 

    end do; end do; end do 

    ! store the number of elements 
    nmats = itel 
    
    ! sort the matrices for basic cache efficiency (not enabled because no test case yet)
    !iid=iss; jjd=jss; nmats = itel 
    !call compress_piece(1,nmats,iid,jjd,mats1)
    !iid=iss; jjd=jss; nmats = itel 
    !call compress_piece(1,nmats,iid,jjd,mats2)
    !iid=iss; jjd=jss; nmats = itel 
    !call compress_piece(1,nmats,iid,jjd,mats3)
    !iid=iss; jjd=jss; nmats = itel 
    !call compress_piece(1,nmats,iid,jjd,mats4) 
    !also need a dummy real array for id1   
    !call compress_piece(1,nmats,iss,jss,id1)  


  endif 
  
  deallocate(iid)
  deallocate(jjd)

end subroutine krook_init 

!------------------------------------------------------------------------------
!> This subroutine contains the Krook operator. Note that this operator uses 
!> an integral over delta f and can therefore be introduced in the code as a 
!> field. This is not done at present since it first has to prove its 
!> usefulness
!> NOTE TOO that this field is very closely related to the gyro-centre density 
!> that is anyway calculated to obtain the potential. It is likely a good idea
!> to combine the two and too save computational time.  
!------------------------------------------------------------------------------
subroutine add_krook(fdis,rhs)

  use control,      only : dtim 
  use grid,         only : parallel_mu, parallel_vpar
  use grid,         only : nmod, nx, ns, nsp 
  use dist,         only : nsolc 
  use constants,    only : c1 
  use mpiinterface, only : mpiallreduce_sum_inplace
  use mpicomms,     only : COMM_VPAR_NE_MU_NE 

  complex, intent(in)    :: fdis(nsolc)
  complex, intent(inout) :: rhs(nsolc)

  ! complex timestep 
  complex :: dtim_cmplx 

  ! loop variables 
  integer :: i

  ! use complex dtim
  dtim_cmplx = c1*dtim

  if (krook_option == 1) then 

    ! initialize sss to zero 
    sss(:) = 0. 

    ! calculate the integral 
    do i = 1, nmats 
      sss(iss(i)) = sss(iss(i)) + mats1(i)*fdis(jss(i)) 
    end do 

    ! mpi reduce over velocity space if needed
    if (parallel_mu .or. parallel_vpar) then

      ! do the sum over velocity space    
      call mpiallreduce_sum_inplace(sss,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)

      ! call mpiallreduce_sum(sss,sss_tmp,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)
      ! sss(:) = sss_tmp(:) 

    endif 

    ! add the damping to the solution 
    if (parallel_vpar) then
      do i = 1, nmats 
        ! the fdis_nvpar array holds the distribution on the opposite processor
        ! it is filled in exp_integration using mpistart(gckrook) and mpiwait(gckrook)
        rhs(jss(i)) = - gammak * dtim_cmplx * ( 0.5*( fdis(id1(i))+fdis_nvpar(id2(i)) ) - &
                    &   mats2(i)*sss(iss(i)) )
      end do 
    else
      do i = 1, nmats 
        rhs(jss(i)) = - gammak * dtim_cmplx * ( 0.5*( fdis(id1(i))+fdis(id2(i)) ) - &
                    &   mats2(i)*sss(iss(i)) )
      end do     
    end if

  endif

  if (krook_option == 3) then 

    ! initialize sss to zero 
    sss(:) = 0. 

    ! calculate the integral 
    do i = 1, nmats 
      sss(iss(i)) = sss(iss(i)) + mats1(i)*fdis(jss(i)) 
    end do 

    ! mpi reduce over velocity space if needed
    if (parallel_mu .or. parallel_vpar) then

      ! do the sum over velocity space    
      call mpiallreduce_sum_inplace(sss,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)

      ! call mpiallreduce_sum(sss,sss_tmp,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)
      ! sss(:) = sss_tmp(:) 

    endif 

    ! add the damping to the solution 
    do i = 1, nmats 
        ! the fdis_nvpar array holds the distribution on the opposite processor
        ! it is filled in exp_integration using mpistart(gckrook) and mpiwait(gckrook)
        rhs(jss(i)) = - gammak * dtim_cmplx * ( fdis(id1(i)) - &
                    &   mats2(i)*sss(iss(i)) )
    end do 

  endif

  !This option is the same of krook=1 but, with krook=3 for the n=0 toroidal mode
  if (krook_option == 4) then 

    ! initialize sss to zero 
    sss(:) = 0. 

    ! calculate the integral 
    do i = 1, nmats 
      sss(iss(i)) = sss(iss(i)) + mats1(i)*fdis(jss(i)) 
    end do 

    ! mpi reduce over velocity space if needed
    if (parallel_mu .or. parallel_vpar) then

      ! do the sum over velocity space    
      call mpiallreduce_sum_inplace(sss,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)

      ! call mpiallreduce_sum(sss,sss_tmp,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)
      ! sss(:) = sss_tmp(:) 

    endif 

    ! add the damping to the solution 
    if (parallel_vpar) then
      do i = 1, nmats 
        ! the fdis_nvpar array holds the distribution on the opposite processor
        ! it is filled in exp_integration using mpistart(gckrook) and mpiwait(gckrook)
        rhs(jss(i)) = - gammak*prefgam(iss(i))* dtim_cmplx * (pref(iss(i))*fdis(id1(i))+  & 
           & (1.E0-pref(iss(i)))*fdis_nvpar(id2(i)) - mats2(i)*sss(iss(i)) )
      end do 
    else
      do i = 1, nmats 
        rhs(jss(i)) = - gammak * prefgam(iss(i)) * dtim_cmplx * (pref(iss(i))*fdis(id1(i))+ &
           & (1.E0-pref(iss(i)))*fdis(id2(i)) - mats2(i)*sss(iss(i)) )

      end do     
    end if

  endif

  if (krook_option == 2) then 

    ! initialize sss to zero 
    sss(:) = 0.; sst(:) = 0. 

    ! calculate the integral 
    do i = 1, nmats 
      sss(iss(i)) = sss(iss(i)) + mats1(i)*fdis(jss(i)) 
      sst(iss(i)) = sst(iss(i)) + mats3(i)*fdis(jss(i))
    end do 

    ! mpi reduce over mu, and vpar if needed
    if (parallel_mu .or. parallel_vpar) then

      ! do the sum over velocity space
      call mpiallreduce_sum_inplace(sss,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)
      call mpiallreduce_sum_inplace(sst,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)      
      
      !call mpiallreduce_sum(sss,sss_tmp,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)
      !sss(:) = sss_tmp(:) 
      !call mpiallreduce_sum(sst,sss_tmp,nmod*nx*ns*nsp,COMM_VPAR_NE_MU_NE)
      !sst(:) = sss_tmp(:) 

    endif 

    ! add the damping to the solution 
    do i = 1, nmats 
      rhs(jss(i)) = - gammak * dtim_cmplx * ( fdis(id1(i)) -  &
                  &   mats2(i)*sss(iss(i)) - mats4(i)*sst(iss(i)) )
    end do 

  endif

end subroutine add_krook 

!------------------------------------------------------------------------------
! Simple index function for the sss array. (no ghost cells need to be 
! considered.
!------------------------------------------------------------------------------
function indx_sss(imod,ix,i,is) 

  use grid, only : nmod, nx, ns 

  integer :: indx_sss 
  integer, intent(in) :: imod, ix, i, is 

  indx_sss = imod + (ix-1)*nmod + (i-1)*nmod*nx + (is-1)*nmod*nx*ns 

end function indx_sss 


end module krook 
