! ----------------------------------------------------------------------------
!> Module that contains the implicit integration using UMFPACK 
!> Developmental: Many things are not yet implemented or tested.
!-----------------------------------------------------------------------------
module imp_integration

  use global, only : iumf, rumf
  use matrix_format, only : sparse_matrix

  implicit none 

  private 
  
  public :: imp_int

#if defined(umfpack)
  !> store separately the size of the matrix nrcol
  integer(kind=iumf), save :: nrcol  ! size of the real matrix (nrcol * nrcol) 

  !> The total number of elements in the real matrix (in integer(kind=iumf)) 
  integer(kind=iumf), save :: ntelem

  !> The number of diagonal elelments connected with f (not the fields) 
  integer(kind=iumf), save :: ndia 

  !> The index arrays (note that they have to be in integer(kind=iumf) 
  integer(kind=iumf), allocatable, save :: ap(:), ai(:)
  
  complex, allocatable, save :: rnl(:), rnlold(:)

  !> The real matrix and vectors 
  real(rumf), allocatable, save :: ax(:), x(:), b(:), c(:), y(:), z(:)

  !> help arrays that contain the distribution 
  real(rumf), allocatable, save :: fdis1(:), fdis2(:), sourcr(:)

  !> The implicit time step 
  real(rumf), save :: dtime_imp, dtime 

  !> parameters needed by UMFPACK
  real(rumf), save :: control_umf(20), control_umf2(20)
  real(rumf), save :: info (90), info2(90) 

  !> pointers in UMFPACK
  integer (iumf), save :: numeric, symbolic, numeric2, symbolic2

  !> help arrays for the second method 
  ! complex, allocatable, save :: mth(:)
  ! integer, allocatable   , save :: ih(:)
  integer, allocatable   , save :: iah(:)
  type(sparse_matrix), save :: mth

  !> Matrix for the solution of the potential 
  complex, allocatable, save :: gmint(:), gmat_buf(:)
  ! complex, save :: gdum 
  integer, allocatable   , save :: ig(:)
  ! complex, allocatable, save :: gmat(:)
  ! integer, allocatable   , save :: iig(:)
  type(sparse_matrix), save :: gmat
  ! integer, save :: idum 

  !> real matrix for the potential inversion 
  integer(kind=iumf)                     , save :: gsize
  real(kind=rumf),           allocatable , save :: gmatr(:)
  integer(kind=iumf),        allocatable , save :: igr(:), iigr(:)

  !> loading of phi 
  !> The number of columns in the field matrix (different on each proc.)
  integer, save :: n_phi_elem_pc
  ! The column numbers solved on a particular processor 
  integer, allocatable, save :: isel_phi_proc(:)
  ! The reduction array for the field solve 
  integer, allocatable, save :: ireduce_field(:)

  real(kind=rumf), parameter    :: half = 0.5_rumf
  complex(kind=rumf), parameter :: ci1  = (0.0_rumf,1.0_rumf)

#endif

contains 


!------------------------------------------------------------------------------
!> Routine that does the initialization necessary for the inversion of the 
!> matrix using method1. This method treats the whole matrix and can be 
!> potentially very slow, and furthermore can run out of memory very quickly 
!------------------------------------------------------------------------------
subroutine imp_init_meth1

  ! all of this only if umfpack is linked, and only for double precision
  use general, only : gkw_abort
#if defined(umfpack)
  use matdat,  only : iac, mat, source, mat_field_diag
  use dist,    only : nsolc, n_phi_start
  use control, only : dtim
  use global,  only : root_and_verbose
  use matconv, only : crstoreal_umf

  integer(kind=iumf) :: i, j
  integer :: ierr
  logical :: diagonal_elem

  ! FIXME check what n4 is meant to be
  integer :: n4
  n4 = mat_field_diag%nmat

  ! Set the array values for the real matrix 
  nrcol  = 2*nsolc 
  ntelem = 4*n4

  ! allocate the column array
  allocate(ap(nrcol+1),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate ap in imp_integration')

  ! allocate the row array 
  allocate(ai(ntelem),stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate ai in imp_integration')

  ! allocate the real matrix
  allocate(ax(ntelem),stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate ax in imp_integration')

  ! allocate the real matrix
  allocate(sourcr(nrcol),stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate sourcr in imp_integration')

  ! Then convert the complex matrix to a real matrix  
  call crstoreal_umf(nsolc,iac,mat,ntelem,ap,ai,ax)

  ! prepare it for an implicit step 
  ndia      = 2*(n_phi_start-1)    ! the number of diagonal elements of f 
  dtime_imp = 0.5_rumf*dtim          ! the implicit time step 
  call crstoimp(nrcol,ndia,dtime_imp,ap,ai,ax)

  ! simple check on the matrix
  do j = 1, nrcol 
    diagonal_elem = .false. 
    if (ap(j).eq.0) then 
      write(*,*)j, 'zero ap elem '
      call gkw_abort('Problem with the real matrix')
    endif 
    if (ap(j).eq.ap(j+1)) then 
      write(*,*)j, 'Value of j never referenced'
      call gkw_abort('Problem with the real matrix')
    endif
    do i = ap(j),ap(j+1)-1
      if (ai(i).eq.j) diagonal_elem = .true.
      if (ai(i).eq.0) then 
        write(*,*)i,j, ' error zero i' 
        call gkw_abort('Problem with the real matrix')
      endif
      if (ai(i).gt.nrcol) then 
        write(*,*)i,j, ' error i > nrcol ', nrcol 
        call gkw_abort('Problem with the real matrix')
      endif
    end do 
    if (.not.diagonal_elem) then 
      write(*,*) j, 'no diagnonal element'
      call gkw_abort('Problem with the real matrix')
    endif 
  end do 

  ! convert to 0 base (because of the use of C inside UMFPACK) 
  do i = 1, nrcol+1
    Ap(i) = ap(i) - 1
  end do 
  do i = 1, ap(nrcol+1)  ! no -1 here. just subtracted 
    Ai(i) = ai(i) - 1
  end do 

  ! set default parameters
  call umf4def (control_umf)

  control_umf(1) = 2
  call umf4pcon (control_umf)

  ! pre-order and symbolic analysis
  call umf4sym (nrcol, nrcol, Ap, Ai, Ax, symbolic, control_umf, info)

  if (root_and_verbose) then 

    ! print statistics computed so far
    ! call umf4pinf (control, info) could also be done.
    write(*,80) info (1), info (16), (info (21) * info (4)) / 2**20, &
         &      (info (22) * info (4)) / 2**20,info (23), info (24), info (25)
80      format ('symbolic analysis:',/,                            &
     &      '   status:  ', f5.0, /,                               &
     &      '   time:    ', e10.2, ' (sec)'/,                      &
     &      '   estimates (upper bound) for numeric LU:', /,       &
     &      '   size of LU:    ', f10.2, ' (MB)', /,               &
     &      '   memory needed: ', f10.2, ' (MB)', /,               &
     &      '   flop count:    ', e10.2, /                         &
     &      '   nnz (L):       ', f10.0, /                         &
     &      '   nnz (U):       ', f10.0)

  endif 

  ! check umf4sym error condition
  if (info (1) .lt. 0) then
    print *, 'Error occurred in umf4sym: ', info (1)
    call gkw_abort('Problem in imp_integration')
  endif

  ! numeric factorization
  call umf4num (Ap, Ai, Ax, symbolic, numeric, control_umf, info)

  ! print statistics for the numeric factorization
  ! call umf4pinf (control, info) could also be done.

  if (root_and_verbose) then 
    write(*,90) info (1), info (66), (info (41) * info (4)) / 2**20,   &
         &      (info (42) * info (4)) / 2**20, info (43), info (44), info (45)
90      format ('numeric factorization:',/,                       &
     &      '   status:  ', f5.0, /,                              &
     &      '   time:    ', e10.2, /,                             &
     &      '   actual numeric LU statistics:', /,                &
     &      '   size of LU:    ', f10.2, ' (MB)', /,              &
     &      '   memory needed: ', f10.2, ' (MB)', /,              &
     &      '   flop count:    ', e10.2, /                        &
     &      '   nnz (L):       ', f10.0, /                        &
     &      '   nnz (U):       ', f10.0)

  endif 

  ! check umf4num error condition
  if (info (1) .lt. 0) then
    print *, 'Error occurred in umf4num: ', info (1)
    call gkw_abort('Probem encountered in numerc factorization with UMFPACK')
  endif


  ! The matrix Ax is no longer needed now that the LU decomposition is 
  ! in memory. It is reconstructed for the Explicit time step 
  call crstoreal_umf(nsolc,iac,mat,ntelem,ap,ai,ax)


  ! allocate the help arrays for the numerical integration 
  allocate(fdis1(nrcol),stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate fdis1 in imp_integration')

  allocate(fdis2(nrcol),stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate fdis2 in imp_integration')

  ! Finally deal with the source 
  ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
  do i = 1, nsolc 
    sourcr(2*i-1) = real(source(i))
    sourcr(2*i)   = aimag(source(i))
  end do 

! umfpack linked, and double precision
#else
  call gkw_abort('You selected method = IMP but umfpack is not linked')
#endif

end subroutine imp_init_meth1


!-----------------------------------------------------------------------------
!> Routine that does the initialization necessary for the inversion of the 
!> matrix base on method 2. This method treats the fields differently 
!> 
!> The system of equations can be written as 
!> 
!>  M F^n+1 + N G^n+1 = F^n 
!>  P F^n+1 + Q G^n+1 = 0 
!> 
!> where F is the distribution and G is the combination of all field 
!> variables. The system allows for the solution of G 'separtely' 
!> through 
!>
!> [ P M^-1 N - Q ] G^n+1 = P M^-1 F^n
!>
!> After which F^n+1 follows from 
!>
!> F^n+1 = M^-1 ( F^n - N G^n+1 ) 
!> 
!> The gain lies in the fact that the Q matrix is much smaller (ns*ns) and 
!> The right hand side of the first equation can be parallelized over 
!> mu (and nvpar in the case of vp_trap = 1). 
!-----------------------------------------------------------------------------
subroutine imp_init_meth2
  use global,          only : r_tiny
  use general,         only : gkw_abort
#if defined(umfpack)
  use matdat,          only : iac, mat, mat_field_diag
  use dist,            only : nsolc, n_phi_start
  use control,         only : dtim, non_linear
  use rotation,        only : shear_real
  use matconv,         only : crstoreal_umf
  use mpiinterface,    only : mpiallreduce_sum
  use matrix_format, only : create_matrix, finalize_matrix
  integer(kind=iumf)   :: i
  integer   :: j, k, jg, itel, ierr, npr, nm, gmatnm, mm, npelement, j_el
  integer(kind=iumf)   :: sys
  real      :: diagonal = 0.0, phitol
  logical   :: new_j_element
  ! FIXME check what n4 is meant to be
  integer :: n4
  n4 = mat_field_diag%nmat

  ! set the error variable to zero 
  ierr = 0 

  ! set the tollerance of the potential 
  phitol = 0.

  !----------------------------------------------------------------------------
  ! Handle the M matrix : 
  !----------------------------------------------------------------------------

  ! first allocate 
  allocate(iah(n_phi_start),stat=ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate iah in imp_integration')
  endif
  ! allocate(ih(n4),stat=ierr)
  ! if (ierr.ne.0) then 
  !   call gkw_abort('Could not allocate ih in imp_integration')
  ! endif
  ! allocate(mth(n4),stat=ierr)
  ! if (ierr.ne.0) then 
  !   call gkw_abort('Could not allocate mth in imp_integration')
  ! endif
  mth = create_matrix("FIXME I have no idea", n4)

  ! save in math block 1 of the matrix (all elements with both 
  ! ih and jh < n_phi_start 
  itel = 0 
  iah(1) = 1 
  do j = 1, n_phi_start-1 
    do i = iac(j), iac(j+1)-1
      if (mat%ii(i).lt.n_phi_start) then 
        itel = itel + 1 
        mth%ii(itel) = mat%ii(i)
        mth%mat(itel) = mat%mat(i)
      endif 
    end do 
    iah(j+1) = itel + 1
  end do 

  ! Now convert this matrix to a real matrix 
  nrcol  = 2*(n_phi_start-1)
  ntelem = 4*itel 

  ! first do the allocation 
  allocate(ap(nrcol+1),stat=ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate ap in imp_integration')
  endif
  ! allocate the row array  
  allocate(ai(ntelem),stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate ai in imp_integration')
  endif
  ! allocate the real matrix
  allocate(ax(ntelem),stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate ax in imp_integration')
  endif


  ! The conversion (npr is the number of columns, nm the total 
  ! number of elements of the complex matrix) 
  npr = n_phi_start - 1 
  nm  = itel 
  call crstoreal_umf(npr,iah,mth,ntelem,ap,ai,ax)

  ! prepare it for an implicit step 
  ndia      = 2*(n_phi_start-1)    ! the number of diagonal elements of f 
  dtime_imp = 0.6_rumf*dtim          ! the implicit time step 
  call crstoimp(nrcol,ndia,dtime_imp,ap,ai,ax)


  ! convert to 0 base (because of the use of C inside UMFPACK) 
  do i = 1, nrcol+1
    Ap(i) = ap(i) - 1
  end do 
  do i = 1, ntelem
    Ai(i) = ai(i) - 1
  end do 

  ! set default parameters
  call umf4def (control_umf)

  control_umf(1) = 2
  call umf4pcon (control_umf)

  ! pre-order and symbolic analysis
  call umf4sym (nrcol, nrcol, Ap, Ai, Ax, symbolic, control_umf, info)

  ! print statistics computed so far
  ! call umf4pinf (control, info) could also be done.
  write(*,80) info (1), info (16), (info (21) * info (4)) / 2**20, &
       &      (info (22) * info (4)) / 2**20,info (23), info (24), info (25)
80      format ('symbolic analysis:',/,                            &
     &      '   status:  ', f5.0, /,                               &
     &      '   time:    ', e10.2, ' (sec)'/,                      &
     &      '   estimates (upper bound) for numeric LU:', /,       &
     &      '   size of LU:    ', f10.2, ' (MB)', /,               &
     &      '   memory needed: ', f10.2, ' (MB)', /,               &
     &      '   flop count:    ', e10.2, /                         &
     &      '   nnz (L):       ', f10.0, /                         &
     &      '   nnz (U):       ', f10.0)

  ! check umf4sym error condition
  if (info (1) .lt. 0) then
    print *, 'Error occurred in umf4sym: ', info (1)
    write(*,*)nrcol, ap(200), ai(200)
    !stop 'imp 1.0'
    stop 1
  endif

  ! numeric factorization
  call umf4num (Ap, Ai, Ax, symbolic, numeric, control_umf, info)

  ! print statistics for the numeric factorization
  !  call umf4pinf (control, info) could also be done.
  write(*,90) info (1), info (66), (info (41) * info (4)) / 2**20,   &
       &      (info (42) * info (4)) / 2**20, info (43), info (44), info (45)
90      format ('numeric factorization:',/,                       &
     &      '   status:  ', f5.0, /,                              &
     &      '   time:    ', e10.2, /,                             &
     &      '   actual numeric LU statistics:', /,                &
     &      '   size of LU:    ', f10.2, ' (MB)', /,              &
     &      '   memory needed: ', f10.2, ' (MB)', /,              &
     &      '   flop count:    ', e10.2, /                        &
     &      '   nnz (L):       ', f10.0, /                        &
     &      '   nnz (U):       ', f10.0)

  ! check umf4num error condition
  if (info (1) .lt. 0) then
    print *, 'Error occurred in umf4num: ', info (1)
    !stop 'imp 1.2'
    stop 1
  endif

  ! The arrays of ai, ap and ax are no longer needed (but are used later
  ! again for the explicit step) Here they are deallocated. 
  deallocate(ax,stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not deallocate ax in imp_integration')
  endif
  deallocate(ai,stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not deallocate ai in imp_integration')
  endif
  deallocate(ap,stat=ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not deallocate ap in imp_integration')
  endif
  deallocate(iah,stat=ierr) 
  if (ierr.ne.0) then 
    call gkw_abort('Could not deallocate iah in imp_integration')
  endif
  ! deallocate(ih,stat=ierr)
  ! if (ierr.ne.0) then 
  !   call gkw_abort('Could not deallocate ih in imp_integration')
  ! endif
  ! deallocate(mth,stat=ierr)
  ! if (ierr.ne.0) then 
  !   call gkw_abort('Could not deallocate mth in imp_integration')
  ! endif
  call finalize_matrix(mth)

  !----------------------------------------------------------------------------
  ! start building the gmat matrix = P M^(-1) N - Q for the calculation of the 
  ! fields 
  !----------------------------------------------------------------------------

  ! Build the plan for parallelization  
  call parallel_phi_solve 

  ! allocate 
  if (n_phi_elem_pc.ne.0) then 
    ! allocate(gmat(n_phi_elem_pc**2),stat = ierr)
    ! if (ierr.ne.0) then 
    !   call gkw_abort('Could not allocate gmat in imp_integration')
    ! endif
    ! FIXME what is the right size estimate?
    gmat = create_matrix("FIXME I have no idea", max(n_phi_elem_pc**2,jg**2))
    
  endif 
  jg = nsolc - n_phi_start + 1 
  ! index array for gmat 
  allocate(ig(jg+1),stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate ig in imp_integration')
  endif
  ! ! index array ii for gmat 
  ! allocate(iig(jg*jg),stat = ierr)
  ! if (ierr.ne.0) then 
  !   call gkw_abort('Could not allocate jjg in imp_integration')
  ! endif
  ! allocate x array (the solution of Ax = b)
  allocate(x(nrcol),stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate x in imp_integration')
  endif
  ! allocate b array (the solution of Ax = b)
  allocate(b(nrcol),stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate b in imp_integration')
  endif
  ! help array to store the intermediate results of gmat 
  allocate(gmint(jg), stat = ierr) 
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate gmint in imp_integration')
  endif
  allocate(gmat_buf(jg),stat = ierr)
  if (ierr.ne.0) then 
    call gkw_abort('Could not allocate gmat_buf in imp_integration')
  endif


  ! intialize 
  ! if (allocated(gmat)) gmat = 0. 
  gmatnm = 0 
  ig(1) = 1
  j_el = 1 

  k_loop : do k = n_phi_start, nsolc 

    ! first find the diagonal element of the Q matrix (for 
    ! comparison only 
    do i = iac(k), iac(k+1)-1 
      if (mat%ii(i).eq.k) then 
        diagonal = abs(mat%mat(i)) 
      endif 
    end do 
 
    ! index in the gmat matrix 
    jg = k -n_phi_start + 1 

    ! Take one column of the matrix N 
    b = 0. 
    do i = iac(k), iac(k+1)-1 
      if (mat%ii(i).lt.n_phi_start) then 
        b(2*mat%ii(i)-1) = -dtime_imp*real(mat%mat(i)) 
        b(2*mat%ii(i))   = -dtime_imp*aimag(mat%mat(i))
      endif 
    end do 

    ! then do the inverse 
    sys = 0
    x = 0 
    call umf4sol (sys, x, b, numeric, control_umf, info)
    if (info (1) .lt. 0) then
      print *, 'Error occurred in umf4sol: ', info (1)
      call gkw_abort('Prolbem in imp_integration')
    endif
 
    ! then contruct P M^-1 N 

  ! initialize 
  gmint = 0. 

  do j = 1, n_phi_start - 1 
    do i = iac(j), iac(j+1) -1 

      ! test if element is to be stored 
      if (mat%ii(i).ge.n_phi_start) then

        ! test if the ii element maybe exists already 
        !npelement = 0 
        !do mm = ig(jg), gmatnm 
        ! if (iig(mm).eq.mat%ii(i)-n_phi_start+1) then 
        !   npelement = mm
        ! endif 
        !end do 

        ! the does not exist 
        !if (npelement.eq.0) then 

          ! is it big enough to be stored 
        !  if ((abs(mat%mat(i)*(x(2*j-1)+(0.,1.)*x(2*j))).gt.phitol*diagonal).or. &
         !    &  (mat%ii(i)-n_phi_start+1.eq.jg)) then 
         !   gmatnm = gmatnm + 1 
         !   if (gmatnm.gt.(nsolc-n_phi_start+1)**2) then 
         !     call gkw_abort('too many elements in the matrix gmat') 
         !   endif         
         !   gmat(gmatnm) = mat%mat(i)*(x(2*j-1) + (0.,1.)*x(2*j))
         !   iig(gmatnm)  = ii(i) - n_phi_start + 1 
         ! endif 
        !else ! the element exists 
         ! gmat(npelement) = gmat(npelement) + mat%mat(i)*(x(2*j-1) + (0.,1.)*x(2*j))
        !endif 

        gmint(mat%ii(i)-n_phi_start+1) = gmint(mat%ii(i)-n_phi_start+1) + mat%mat(i)*(x(2*j-1) + (0.,1.)*x(2*j))
      endif 

    end do 
  end do 
  
! sum over all processors 
  call mpiallreduce_sum(gmint,gmat_buf,nsolc-n_phi_start+1)

  do i = 1, nsolc-n_phi_start+1 
    gmint(i) = gmat_buf(i) 
  end do 

  ! Does this processor store anything? 

  new_j_element  =  .false. 

  if (n_phi_elem_pc.ne.0) then 
  
    do i = 1, nsolc - n_phi_start + 1 

      ! is the element big enough to be kept ?? 
      if ((abs(gmint(i)).gt.phitol*diagonal).and.abs(gmint(i)) > r_tiny) then 
        ! stored on this processor ?? 
        if (ireduce_field(i+n_phi_start-1).ne.0) then 
          new_j_element = .true. 
          gmatnm = gmatnm + 1 
          gmat%mat(gmatnm) = gmint(i) 
          gmat%ii(gmatnm) = ireduce_field(i+n_phi_start-1)
        endif 
      endif 
    
    end do 
  
  endif 

  if (new_j_element) then 
    j_el = j_el + 1 
    ig(j_el) = gmatnm+1 
  endif 

end do  k_loop 

write(*,*)gmatnm, (nsolc-n_phi_start+1)**2, real(gmatnm) / (nsolc-n_phi_start+1)**2

write(*,*)j_el-1, n_phi_elem_pc

! sort the array (bubble sort MUST be changed) 
!do j = 1, nsolc-n_phi_start+1 
!  do i = ig(j), ig(j+1)-1 
!    npelement = i  
!    do k = i+1, ig(j+1)-1 
!      if (iig(k).lt.iig(npelement)) then 
!        npelement = k 
!      endif 
!    end do 
!    if (npelement.ne.i) then 
!      idum = iig(i) 
!      iig(i) = iig(npelement)
!      iig(npelement) = idum 
!      gdum = gmat(i) 
!      gmat(i) = gmat(npelement)
!      gmat(npelement) = gdum 
!    endif
!  end do 
!end do 


!#if defined(mpi2) 

! first allocate the buffer array 
!allocate(gmat_buf(gmatnm),stat = ierr)
!if (ierr.ne.0) then 
!  call gkw_abort('Could not allocate gmat_buf in imp_integration')
!endif
!
!call MPI_ALLREDUCE(gmat,gmat_buf,gmatnm,MPICOMPLEX_X, MPI_SUM, &
!                  & MPI_COMM_WORLD, ierr) 
!
!do i = 1, gmatnm  
!  gmat(i) = gmat_buf(i) 
!end do 

! buffer array no longer needed 
!deallocate(gmat_buf,stat = ierr)
!if (ierr.ne.0) then 
!  call gkw_abort('Could not deallocate gmat_buf in imp_integration')
!endif

!#endif

if (n_phi_elem_pc.gt.0) then 

! finally add the - Q matrix 
do k = 1, n_phi_elem_pc  
  ! the j-index of q 
  j = isel_phi_proc(k) 
  write(*,*)'k,j ',k,j
  do i = iac(j), iac(j+1)-1 
    if (mat%ii(i).ge.n_phi_start) then 
       npelement = 0 
       do mm = ig(k), ig(k+1)-1 
         if (isel_phi_proc(gmat%ii(mm)).eq.mat%ii(i)) then 
           npelement = mm 
         endif 
       end do 
       if (npelement.ne.0) then 
         gmat%mat(npelement) = gmat%mat(npelement) - mat%mat(i) 
       endif 
!      gmat(mat%ii(i)-n_phi_start+1,j-n_phi_start+1) = gmat(mat%ii(i)-n_phi_start+1,j-n_phi_start+1) - mat%mat(i)
    endif 
  end do 
end do 

! the gmat matrix is finished. 
!open(19, file = 'gmat') 
!write(*,*)'wrote the file' 
!do i = 1, nsolc-n_phi_start+1 
!  write(19,100)(abs(gmat(j,i)), j = 1, nsolc-n_phi_start+1)
!  100 format(500(1pe13.5,1X))
!end do

endif 

! fill the index arrays. 
jg = nsolc - n_phi_start + 1 


write(*,*)'did the gmatrix' 

! Note both x and b are used in the main routine b is deallocate
deallocate(b,stat=ierr) 
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate b in imp_integration')
endif


! allocate the real arrays 
allocate(igr(2*n_phi_elem_pc+1),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate ig in imp_integration')
endif
allocate(iigr(4*n_phi_elem_pc**2),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate jjg in imp_integration')
endif
allocate(gmatr(4*n_phi_elem_pc**2),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate jjg in imp_integration')
endif
allocate(b(2*jg),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate b in imp_integration')
endif
allocate(c(2*jg),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate c in imp_integration')
endif

write(*,*)'second round' 
if (n_phi_elem_pc.ne.0) then 

! The conversion 
gsize = 2*n_phi_elem_pc
call crstoreal_umf(n_phi_elem_pc,ig,gmat,gsize*gsize,igr,iigr,gmatr)

do i = 1, gsize + 1 
  igr(i) = igr(i)-1 
end do 
do i = 1, gsize*gsize
  iigr(i) = iigr(i) - 1 
end do 

! set default parameters
call umf4def (control_umf2)

 control_umf2(1) = 2
call umf4pcon (control_umf2)

! pre-order and symbolic analysis
call umf4sym (gsize, gsize, igr, iigr, gmatr, symbolic2, control_umf2, info2)

! print statistics computed so far
! call umf4pinf (control, info) could also be done.
write(*,80) info2(1), info2(16), (info2(21) * info2(4)) / 2**20, &
     &      (info2(22) * info2(4)) / 2**20,info2(23), info2(24), info2(25)

! check umf4sym error condition
if (info2(1) .lt. 0) then
  print *, 'Error occurred in umf4sym: ', info (1)
  stop 1
endif

! numeric factorization
call umf4num (igr, iigr, gmatr, symbolic2, numeric2, control_umf2, info2)

! print statistics for the numeric factorization
! call umf4pinf (control, info) could also be done.
write(*,90) info2(1), info2(66), (info2(41) * info2(4)) / 2**20,   &
     &      (info2(42) * info2(4)) / 2**20, info2(43), info2(44), info2(45)

! check umf4num error condition
if (info2(1) .lt. 0) then
  print *, 'Error occurred in umf4num: ', info2(1)
  stop 1
endif

endif 

! deallocate the arrays no longer used and reuse some for the time 
! integration 
! start building the gmat matrix for the solution of phi 
if (n_phi_elem_pc.ne.0) then 
! deallocate(gmat,stat = ierr)
! if (ierr.ne.0) then 
!   call gkw_abort('Could not deallocate gmat in imp_integration')
! endif
call finalize_matrix(gmat)
endif 
deallocate(ig,stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate ig in imp_integration')
endif
! deallocate(iig,stat = ierr)
! if (ierr.ne.0) then 
!   call gkw_abort('Could not deallocate iig in imp_integration')
! endif
deallocate(x,stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate x in imp_integration')
endif
deallocate(b,stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate b in imp_integration')
endif
deallocate(igr,stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate igr in imp_integration')
endif
deallocate(iigr,stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate jjg in imp_integration')
endif
deallocate(gmatr,stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate gmatr in imp_integration')
endif
deallocate(c,stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate c in imp_integration')
endif


!----------
! Prepare the matrices and vectors for the time integration 
!----------

! allocate x array (the solution of Ax = b)
jg = nsolc - n_phi_start + 1 
ndia = 2*(n_phi_start -1) 
nrcol = 2*nsolc 
allocate(x(nrcol),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate x in imp_integration')
endif
allocate(y(nrcol),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate y in imp_integration')
endif
allocate(z(nrcol),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate z in imp_integration')
endif
jg = nsolc - n_phi_start + 1 
allocate(b(2*jg),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate b in imp_integration')
endif
allocate(c(2*jg),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate c in imp_integration')
endif

! ap, ai, ax will contain the M+N array for an explicit step 
itel = 0 
do j = 1, nsolc 
  do i = iac(j), iac(j+1)-1 
    if (mat%ii(i).lt.n_phi_start) then 
      itel = itel + 1 
    endif 
  end do 
end do 

! first do the allocation 
ntelem = 4*itel 
allocate(ap(2*nsolc+1),stat=ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate ap in imp_integration')
endif
! allocate the row array 
allocate(ai(ntelem),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate ai in imp_integration')
endif
! allocate the real matrix
allocate(ax(ntelem),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate ax in imp_integration')
endif

! gmatr, ig, iig will contain the N+P arrays for the implicit 
! time step 
itel = 0
do j = 1, nsolc 
  do i = iac(j), iac(j+1)-1 
    if (j.lt.n_phi_start) then 
      if (mat%ii(i).ge.n_phi_start) then 
        itel = itel + 1 
      endif 
    else 
      if (mat%ii(i).lt.n_phi_start) then 
        itel = itel + 1 
      endif 
    endif 
  end do 
end do 

! allocate the real arrays 
allocate(igr(2*nsolc+1),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate ig in imp_integration')
endif
allocate(iigr(4*itel),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate jjg in imp_integration')
endif
allocate(gmatr(4*itel),stat = ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate jjg in imp_integration')
endif

! then allocate the complex help matrix (and index arrays) 
allocate(iah(nsolc+1),stat=ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate iah in imp_integration')
endif
mth = create_matrix("FIXME I have no idea", n4)
! allocate(ih(n4),stat=ierr)
! if (ierr.ne.0) then 
!   call gkw_abort('Could not allocate ih in imp_integration')
! endif
! allocate(mth(n4),stat=ierr)
! if (ierr.ne.0) then 
!   call gkw_abort('Could not allocate mth in imp_integration')
! endif

! fill the ap, ai, ax arrays 
itel = 0 
iah(1) = 1 
do j = 1, nsolc 

  do i = iac(j), iac(j+1)-1 
    if (mat%ii(i).lt.n_phi_start) then 
      itel = itel + 1 
      mth%ii(itel) = mat%ii(i) 
      mth%mat(itel) = mat%mat(i) 
    endif 
  end do 
  iah(j+1) = itel + 1 
end do 

! convert this matrix 
call crstoreal_umf(nsolc,iah,mth,ntelem,ap,ai,ax)


! fill the gmatr ig and iig arrays 
itel = 0 
iah(1) = 1 
do j = 1, n_phi_start-1 
  do i = iac(j), iac(j+1)-1 
    if (mat%ii(i).ge.n_phi_start) then 
      itel = itel + 1 
      mth%ii(itel) = mat%ii(i)
      mth%mat(itel) = mat%mat(i) 
    endif 
  end do 
  iah(j+1) = itel + 1 
end do 


do j = n_phi_start, nsolc
 
  do i = iac(j), iac(j+1)-1 
    if (mat%ii(i).lt.n_phi_start) then 
      itel = itel+1 
      mth%ii(itel) = mat%ii(i) 
      mth%mat(itel) = mat%mat(i) 
    endif 
  end do 
  iah(j+1) = itel+1 
end do  

! convert this array 
ntelem = 4*itel 
call crstoreal_umf(nsolc,iah,mth,ntelem,igr,iigr,gmatr)

! convert to implicit array 
call crstoimp(nrcol,ndia,dtime_imp,igr,iigr,gmatr)

! finally subtract 2*(n_phi_start-1) from iigr for storage in a shorter 
! vector 
do i = 1, igr(nrcol+1)-1 
  if (iigr(i).gt.ndia) then 
    iigr(i) = iigr(i) - ndia 
  endif 
end do 

! the complex matrix is no longer needed. 
deallocate(iah,stat=ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate iah in imp_integration')
endif
! deallocate(ih,stat=ierr)
! if (ierr.ne.0) then 
!   call gkw_abort('Could not deallocate ih in imp_integration')
! endif
! deallocate(mth,stat=ierr)
! if (ierr.ne.0) then 
!   call gkw_abort('Could not deallocate mth in imp_integration')
! endif
call finalize_matrix(mth)

! finally allocate the arrays necessary for the nonlinear terms 
if (non_linear.or.shear_real) then 
allocate(rnl(nsolc),stat=ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate rnl in imp_integration')
endif
allocate(rnlold(nsolc),stat=ierr)
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate rnlold in imp_integration')
endif
endif 

#else 
  call gkw_abort('You selected method = IMP but umfpack is not linked')
#endif

end subroutine imp_init_meth2


!------------------------------------------------------------------------------
!> Implicit time integration. Call the correct initialization routine 
!------------------------------------------------------------------------------
subroutine imp_int 

  use control, only : meth

  select case(meth) 

  case(1) 
  ! the full matrix is taken into an implicit step this does not parallelize 
  ! over anything and is potentially very slow. 
  call imp_int_meth1 

  case(2) 
  ! The potential is solved for independently 
  call imp_int_meth2

  end select 
 
end subroutine imp_int 

!------------------------------------------------------------------------------
!> Routine that does the time stepping for method 1. 
!------------------------------------------------------------------------------
subroutine imp_int_meth1 

  use general,         only : gkw_abort
! all of this only if umfpack is linked 
#if defined(umfpack)
  use control,         only : naverage, time, dtim, ntotstep
  use dist,            only : fdisi, nsolc
  use diagnostic,      only : diagnostic_pre_naverage

  logical, save :: initialized = .false.  
  integer(kind=iumf) :: sys, i, j
  integer            :: iloop 
    
  dtime = real(dtim,rumf)

  ! When called the first time the matrix needs to be manipulated 
  if (.not.initialized) then 
    call imp_init_meth1 
    initialized = .true.
  endif 


  ! inside this routine use the real representation 
  ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
  do i = 1, nsolc 
    fdis1(2*i-1) = real(fdisi(i))
    fdis1(2*i)   = aimag(fdisi(i))
  end do 


  time_stepping : do iloop = 1, naverage

    ! Store the present distribution function 
    do i = 1, nrcol 
      fdis2(i) = fdis1(i)
    end do 

    ! do an explicit time step 
    do j = 1, nrcol;  do i = ap(j), ap(j+1)-1 
      fdis2(ai(i)) = fdis2(ai(i)) + half*dtime*ax(i)*fdis1(j)
    end do;  end do 

    ! Add the source 
    do j = 1, nrcol 
      fdis2(j) = fdis2(j) + dtime*sourcr(j)
    end do 

    ! before going into the implicit solver set the fields to zero 
    do i = ndia+1, nrcol  
      fdis2(i) = 0.0_rumf 
    end do 

    ! solve Ax=b, without iterative refinement
    sys = 0
    call umf4sol (sys, fdis1, fdis2, numeric, control_umf, info)
    if (info (1) .lt. 0) then
      print *, 'Error occurred in umf4sol: ', info (1)
      call gkw_abort('Prolbem in imp_integration')
    endif

    ! advance the time 
    time=time+dtim*1.

  ! Some diagnostics (such as mode frequencies) need values from the
  ! previous time step, so need to be updated before naverage steps have
  ! been completed.
  !if (iloop > naverage - 2) then
    !copy back into fdisi
    ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
    do i = 1, nsolc
      fdisi(i) = fdis1(2*i-1) + ci1*fdis1(2*i)
    end do
    call diagnostic_pre_naverage(iloop, fdisi)
  !end if
   
    ! one more timestep done 
    ntotstep = ntotstep + 1 

  end do time_stepping

  ! copy the result back into the complex array fdisi 
  ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
  do i = 1, nsolc
    fdisi(i) = fdis1(2*i-1) + ci1*fdis1(2*i)
  end do 

  ! umfpack not linked 
#else
  call gkw_abort('you have selected implicit but umfpack is not linked')
#endif

end subroutine imp_int_meth1


!------------------------------------------------------------------------------
!> Routine that does the time stepping
!------------------------------------------------------------------------------
subroutine imp_int_meth2

use general,          only : gkw_abort, gkw_warn
! all of this only if umfpack is linked
#if defined(umfpack) 
use control,          only : naverage, time, dtim, ntotstep,  &
                           & non_linear, dtim_est,  stop_me, dt_min
use dist,             only : fdisi, nsolc, n_phi_start
use rotation,         only : shear_real
use non_linear_terms, only : add_non_linear_terms
use diagnostic,       only : diagnostic_pre_naverage
use mpiinterface,     only : mpiallreduce_sum, mpibarrier


logical, save :: initialized = .false.  
integer(kind=iumf) :: sys, i, j  
integer   :: iloop, ihelp 

dtime = real(dtim,rumf)

! When called the first time the matrix needs to be manipulated 
if (.not.initialized) then 
  call imp_init_meth2 
  initialized = .true.
endif 

do i = 1, nsolc
  ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
  y(2*i-1) = real(fdisi(i)) 
  y(2*i)   = aimag(fdisi(i))
end do 

time_stepping : do iloop = 1, naverage


  ! nonlinear terms if desired 
  if (non_linear.or.shear_real) then 
    do i = 1, nsolc 
      rnlold(i) = rnl(i)
      rnl(i) = 0.
      ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
      fdisi(i) = y(2*i-1) + ci1*y(2*i) 
    end do 
    call add_non_linear_terms(fdisi,rnl) 
    write(*,*)'1234567890 ',dtim_est
  endif 


! if (nlapar) then 
!    call gkw_abort('not in this version')
!    do i = 1, nmata
!      y(2*i-1) = y(2*i-1) + dreal(mata(i)*(y(2*jja(i)-1) + (0.,1.)*y(2*jja(i))))
!      y(2*i)   = y(2*i)   + dimag(mata(i)*(y(2*jja(i)-1) + (0.,1.)*y(2*jja(i))))
!    end do
!  endif 

  ! do half an explicit step
  do i = 1, nrcol 
    z(i) = y(i) 
  end do 
  do j = 1, nrcol  
    do i = ap(j), ap(j+1)-1 
      z(ai(i)) = z(ai(i)) + (dtime-dtime_imp)*ax(i)*y(j)  !use to be 0.5*dtim*ax*y
    end do 
  end do 

   
  ! first calculate the right hand side of the equation for the 
  ! fields 
  sys = 0
  ! \attention Here might be a zero passed/input that leads to a zero; by this
  !            value is then divided.
  call umf4sol (sys, x, z, numeric, control_umf, info)
  if (info (1) .lt. 0) then
    print *, 'Error occurred in umf4sol: ', info (1)
    call gkw_abort('Problem in imp_integration : inversion of fdisi')
  endif

  !write(*,*)'at the barrier', processor_number, sum(x)
  call mpibarrier()
  !write(*,*)'before filling b' 

  b = 0. 
  do j = 1, ndia
    do i = igr(j), igr(j+1)-1 
      !if (processor_number.eq.0) then 
      !  write(*,*)i,j,igr(j),iigr(i),gmatr(i)
      !endif
      b(iigr(i)) = b(iigr(i)) + gmatr(i)*x(j)
    end do 
  end do 

  !write(*,*)processor_number, 'before global sum xxx', sum(b)
  ! Do a global sum 
  c = 0.
  !write(*,*)'at the barrier 2' 
  !call mpi_barrier(mpi_comm_world,ierr)

  call mpiallreduce_sum(b,c,2*(nsolc-n_phi_start+1))

  !write(*,*)processor_number,'at the barrier 3 xxx' 
  !call mpi_barrier(mpi_comm_world,ierr)

  !write(*,*)sum(c)
  ! prepare for the solve (if necessary) 
  b = 0.0_rumf
  if (n_phi_elem_pc.ne.0) then 

    ! do the reduction here 
    do i = 1, nsolc-n_phi_start+1 
      ihelp = ireduce_field(i+n_phi_start-1)
      if (ihelp.gt.0) then 
        b(2*ihelp-1) = c(2*i-1)
        b(2*ihelp) = c(2*i)
      endif 
    end do 

    
   ! write(*,*)processor_number,'at the 4 xxx',sum(b), sum(c)

    sys = 0
    call umf4sol (sys, c, b, numeric2, control_umf2, info2)
    if (info2(1) .lt. 0) then
      print *, 'Error occurred in umf4sol: ', info2(1)
      call gkw_abort('Prolbem in imp_integration : field solve')
    endif

    !write(*,*)processor_number, 'solved', sum(c)

    ! undo the compression 
    b = 0.0_rumf 
    do i = 1, n_phi_elem_pc 
      ihelp = isel_phi_proc(i)-n_phi_start+1
      b(2*ihelp-1) = c(2*i-1)
      b(2*ihelp) = c(2*i) 
    end do 

    !write(*,*)sum(b)
  endif 

  ! After the total reduction, c contains the implicitly calculted 
  ! fields: field^n+1 
  call mpiallreduce_sum(b,c,2*(nsolc-n_phi_start+1))
 
  ! calculate the right hand side of the f equation 
  ! Load in x the result of the half explicit step 
  do i = 1, nrcol 
    x(i) = z(i)
  end do 
  do j = ndia+1, nrcol 
    ! x <-- g_^n - N field^n+1   (both g^n and field^n+1 are calculated 
    ! from the half explicit step 
    do i = igr(j), igr(j+1) -1 
      x(iigr(i)) = x(iigr(i)) - gmatr(i)*c(j-ndia) 
    end do 
  end do 
  ! Add the nonlinear terms 
  if (non_linear.or.shear_real) then 
    do i = 1, n_phi_start - 1 
      x(2*i-1) = x(2*i-1) + real(1.5*rnl(i) - 0.5*rnlold(i))
      x(2*i)   = x(2*i)   + aimag(1.5*rnl(i) - 0.5*rnlold(i))
    end do 
  endif 

  ! first calculate the right hand side of the equation for the 
  ! fields 
  sys = 0
  call umf4sol (sys, y, x, numeric, control_umf, info)
  if (info (1) .lt. 0) then
    print *, 'Error occurred in umf4sol: ', info (1)
    call gkw_abort('Prolbem in imp_integration')
  endif

  ! copy back in the potential solution 
  do i = ndia+1, nrcol 
    y(i) = c(i-ndia)
  end do 

  ! advance the time 
  time=time+dtim*1.

  ! Some diagnostics (such as mode frequencies) need values from the
  ! previous time step, so need to be updated before naverage steps have
  ! been completed.
  !if (iloop > naverage - 2) then
    ! copy back into fdisi 
    ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
    do i = 1, nsolc 
      fdisi(i) = y(2*i-1) + (0.,1.)*y(2*i) 
    end do 
    call diagnostic_pre_naverage(iloop, fdisi)
  !end if
   
  ! one more timestep done 
  ntotstep = ntotstep + 1 

end do time_stepping

! check if the timestep is too small; if so, end the run
  if (dtim < dt_min) then
    call gkw_warn('dt < dt_min; the run will terminate shortly')
    stop_me = .true.
  end if
  
! copy back into fdisi 
! WARNING: IMPLICIT UMFPACK TYPE CONVERSION
do i = 1, nsolc 
  fdisi(i) = y(2*i-1) + (0.,1.)*y(2*i) 
end do 

#else 
  call gkw_abort('IMP is selected but UMFPACK is not linked')
#endif

end subroutine imp_int_meth2

#if defined(umfpack)
!------------------------------------------------------------------------------
!> Calculate the residual
!------------------------------------------------------------------------------
subroutine resid (n, nz, Ap, Ai, Ax, x, b, r)


 
integer(kind=iumf) :: n, nz, Ap (n+1), Ai (nz), j, i, p
real(kind=rumf)    :: Ax (nz), x (n), b (n), r (n), rmax, aij

do i = 1, n
  r (i) = -b (i)
end do 

do j = 1,n
  do p = Ap (j) + 1, Ap (j+1)
    i = Ai (p) + 1
    aij = Ax (p)
    r (i) = r (i) + aij * x (j)
  end do 
end do 

rmax = 0
do i = 1, n
  rmax = max (rmax, r (i))
end do 

print *, 'norm (A*x-b): ', rmax
end subroutine resid

!-----------------------------------------------------------------
!>     this routine converts the real matrix to the implicit format 
!>
!>     input 
!>     np     number of colums in the matrix 
!>     ndia   number of elements in the matrix that need the 
!>            diagonal element 1. added 
!>     apk    index array (np + 1) 
!>     aik    index array (ia(np+1)-1) 
!>     axk    Real matrix. (should be ordered) 
!>
!>     output 
!>    axk    Real matrix prepared for implicit run (compressed)
!> 
!>     NOTE call this routine with ap(1) = 1 (not with the 
!>          c-convention) 
!------------------------------------------------------------------
subroutine crstoimp(np,nd,dtime,apk,aik,axk)

integer(kind=iumf)        :: np, nd
integer(kind=iumf)        :: apk(np+1) 
integer(kind=iumf)        :: aik(ap(np+1)-1)  
real(kind=rumf)    :: axk(ap(np+1)-1), dtime 

integer(kind=iumf)        :: i, j 

do j = 1, np 
  do i = apk(j), apk(j+1)-1
    if (aik(i).le.nd) then 
      axk(i) = -dtime*axk(i) 
      if (aik(i).eq.j) then 
        axk(i) = axk(i) + 1. 
      endif 
    endif 
  end do 
end do  

end subroutine crstoimp

!------------------------------------------------------------------------------
!> Modifies the real matrix in order to solve for the potential 
!------------------------------------------------------------------------------
subroutine strip_real(np,nd,ia,ja,a)
integer(kind=iumf)        :: np, nd 
integer(kind=iumf)        :: ia(np+1) 
integer(kind=iumf)        :: ja(ia(np+1)-1) 
real(kind=rumf)           :: a(ia(np+1)-1)

integer(kind=rumf) :: i, j, idiag
real(kind=rumf)    :: diael

! modify the potential part into an explict from 
do i = nd+1, np 

  ! first find the diagonal element 
  idiag = 0
  do j = ia(i), ia(i+1)-1 
    if (ja(j).eq.i) then 
      idiag = j 
    endif 
  end do 

  if (idiag.eq.0) then 
    !stop 'error'
    write(*,*) 'error'
    stop 1
  endif 

  ! recalculate the matrix elements 
  diael = a(idiag) 
  do j = ia(i), ia(i+1)-1 
    a(j) = -a(j) / diael 
    if (ja(j).eq.idiag) then 
      a(j) = 0. 
    endif 
  end do 

end do   

end subroutine strip_real 

!------------------------------------------------------------------------------
!> Calulate the fields for an implicit step 
!------------------------------------------------------------------------------
subroutine calculate_fields_imp(fdis,np,nd,ia,ja,a)
integer(kind=iumf)        :: np, nd
integer(kind=iumf)        :: ia(np+1)
integer(kind=iumf)        :: ja(ia(np+1)-1) 
real(kind=rumf)           :: fdis(np), a(ia(np)-1) 

integer(kind=rumf) :: i, j

do i = nd + 1, np 
  fdis(i) = 0 
  do j = ia(i), ia(i+1)-1 
    fdis(i) = fdis(i) + a(j)*fdis(ja(j))
  end do 
end do 

end subroutine calculate_fields_imp 

!-------------------------------------------------------------------------------
!> This subroutine makes a plan for the solution of the fields in parallel 
!> The minimum amount of work in the field solve corresponds to the solution 
!> for the mode with the longest field line. Subdividing in smaller parts is 
!> not useful since it will not reduce the cost. The maximum number of processors
!> for the field solve is then nmod * ikxspace 
!-------------------------------------------------------------------------------
subroutine parallel_phi_solve 

use mpiinterface,   only : number_of_processors, processor_number, mpibarrier
use grid,           only : nmod, nx, ns
use dist,           only : iphi, nsolc, iapar
use index_function, only : indx
use mode,           only : ikxspace
use general,        only : gkw_abort
use control,        only : nlapar 

integer, allocatable :: ikpc(:) 
integer :: ikp, i, ierr, iktel, imod, ikpm, it, ix, xs, n_phi_c

! Set the error parameter
ierr = 0 

! allocate the storage array 
allocate(ikpc(number_of_processors+1), stat = ierr) 
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate ikpc in imp_integration')
endif 

! number of submodes on each processor (note ikp is rounded) 
ikp = nmod*ikxspace / number_of_processors 

! fill the array that determines the start value on each processor
ikpc(1) = 0
do i = 1, number_of_processors 
  ikpc(i+1) = ikp 
end do 
! left over 
ikpm = mod(nmod*ikxspace,number_of_processors)
do i = 1, ikpm 
  ikpc(number_of_processors-i+2) = ikpc(number_of_processors-i+2) + 1 
end do 
it = 1 
do while (ikpc(it+1).eq.0) 
  it = it + 1 
end do 
ikpc(it) = 1
do i = it+1, number_of_processors+1
  ikpc(i) = ikpc(i-1) + ikpc(i) 
end do 


! initialize
n_phi_elem_pc = 0 

! simple but effective 
iktel = 0 
do imod = 1, nmod 
  do ix = 1, ikxspace
    iktel = iktel + 1 

    ! check if the value is in the range for this processor 
    if ((iktel.ge.ikpc(processor_number+1)).and. &
      & (iktel.lt.ikpc(processor_number+2))) then 

      ! count the number of elements 
      n_phi_elem_pc = n_phi_elem_pc +  ((nx-ix)/ikxspace + 1)*ns

      if (nlapar) then 
        n_phi_elem_pc = n_phi_elem_pc + ((nx-ix)/ikxspace+1)*ns 
      endif 

    endif 
  end do 
end do 

write(*,*)'before allocation'
if (n_phi_elem_pc.eq.0) goto 10 

! allocate the index arrays 
allocate(isel_phi_proc(n_phi_elem_pc), stat = ierr) 
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate isel_phi_proc in imp_integration')
endif 
allocate(ireduce_field(nsolc), stat = ierr) 
if (ierr.ne.0) then 
  call gkw_abort('Could not allocate ireduce_field in imp_integration')
endif 



! set the element number to zero 
n_phi_c = 0 
iktel = 0 
do imod = 1, nmod 
  do ix = 1, ikxspace
    iktel = iktel + 1 

    ! check if the value is in the range for this processor 
    if ((iktel.ge.ikpc(processor_number+1)).and. &
      & (iktel.lt.ikpc(processor_number+2))) then 

        xs = ix 
        do while(xs.le.nx)  
          do i = 1, ns 
            n_phi_c = n_phi_c + 1 
            isel_phi_proc(n_phi_c) = indx(iphi,imod,xs,i)
          end do 
          xs = xs + ikxspace 
        end do 
     
        if (nlapar) then 
          xs = ix 
          do while(xs.le.nx)  
            do i = 1, ns 
              n_phi_c = n_phi_c + 1 
              isel_phi_proc(n_phi_c) = indx(iapar,imod,xs,i)
            end do 
            xs = xs + ikxspace 
          end do 
        endif 

    endif 
  end do 
end do 

if (n_phi_c.ne.n_phi_elem_pc) then 
  call gkw_abort('Internal error in imp_integration')
endif 

! deallocate the unneeded arrays 
deallocate(ikpc, stat = ierr) 
if (ierr.ne.0) then 
  call gkw_abort('Could not deallocate ikpc in imp_integration')
endif 

! sort the array 
call sort(n_phi_elem_pc,isel_phi_proc) 

write(*,*)'do ireduce' 

! make the reduction array 
ireduce_field = 0 
do i = 1, n_phi_elem_pc 
  ireduce_field(isel_phi_proc(i)) = i 
end do 


! a barrier just to be on the safe side (can be removed) 
10 continue 
call mpibarrier()

end subroutine parallel_phi_solve 

#endif

!------------------------------------------------------------------------------
!> Simple sorting algorithm
!------------------------------------------------------------------------------
subroutine sort(n,ra)
integer :: n, ra(n), l, ir, j, rra, i

l = n/2+1 
ir = n 

10 continue 
  if (l.gt.1) then 
    l = l-1
    rra= ra(l)
  else 
    rra = ra(ir)
    ra(ir) = ra(1) 
    ir = ir -1 
    if (ir.eq.1) then 
      ra(1) = rra
      return 
    endif 
  endif
  i = l 
  j = l+l 
  20 continue 
    if (j.le.ir) then 
      if (j.lt.ir) then 
        if (ra(j).lt.ra(j+1)) j = j + 1 
      endif 
      if (rra.lt.ra(j)) then 
        ra(i) = ra(j)
        i = j 
        j = j + j 
      else 
        j = ir + 1 
      endif 
      goto 20 
    endif 
  ra(i) = rra 
goto 10 

end subroutine sort 

end module imp_integration 
