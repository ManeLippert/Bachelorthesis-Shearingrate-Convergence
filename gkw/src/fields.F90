!------------------------------------------------------------------------------
!> Module containing routines for returning gyroaveraged fields
!> Will contain more soon, if not could be merged with gyro_average
!------------------------------------------------------------------------------
module fields 
  
  use global, only : rumf 
  
  implicit none 

  public :: get_averaged_phi, get_averaged_apar, get_averaged_bpar 
  public :: fields_allocate, fields_deallocate, calculate_fields
  public :: field_solve_nonspec_wrap, calculate_cons, gyro_average_fields

  private 

  !> help array for the mpi
  complex, allocatable, save :: bufphi(:)

  !> help array for the nonspectral polarisation of the potential 
  real(rumf), allocatable, save :: z_buf(:),  z(:)
  complex(rumf), allocatable :: phi_dum(:)

  !>  help array to store the additional potential which is of the order rhostar
  real(rumf), allocatable :: phi1(:)

contains 

!------------------------------------------------------------------------------
!> This function gives the value of the gyro-averaged potential for a given 
!> grid point. If nlphi = .false. a zero is returned 
!------------------------------------------------------------------------------
function get_averaged_phi(imod,ix,i,j,is,fdis)
  use control,        only : nlphi, spectral_radius 
  use index_function, only : indx 
  use functions,      only : besselj0_gkw
  use dist,           only : iphi, iphi_ga 

  integer, intent(in)  :: imod, ix, i, j, is
  complex, intent(in)  :: fdis(:)
  complex :: get_averaged_phi 

  if (nlphi) then 
    if (spectral_radius) then 
      get_averaged_phi = besselj0_gkw(imod,ix,i,j,is) * &
                       & fdis(indx(iphi,imod,ix,i)) 
    else 
      get_averaged_phi = fdis(indx(iphi_ga,imod,ix,i,j,is)) 
    endif 
  else 
    get_averaged_phi = (0.E0,0.E0) 
  endif 
 
end function get_averaged_phi 

!------------------------------------------------------------------------------
!> This function gives the value of the gyro-averaged vector potential for a  
!> given grid point. If nlapar = .false. a zero is returned 
!------------------------------------------------------------------------------
function get_averaged_apar(imod,ix,i,j,is,fdis) 

  use control,        only : nlapar, spectral_radius 
  use index_function, only : indx 
  use functions,      only : besselj0_gkw
  use dist,           only : iapar, iapar_ga

  integer, intent(in)  :: imod, ix, i, j, is
  complex, intent(in)  :: fdis(:)  
  complex :: get_averaged_apar 

  if (nlapar) then 
    if (spectral_radius) then 
      get_averaged_apar = besselj0_gkw(imod,ix,i,j,is) * &
                        & fdis(indx(iapar,imod,ix,i)) 
    else 
      get_averaged_apar = fdis(indx(iapar_ga,imod,ix,i,j,is)) 
    endif 
  else 
    get_averaged_apar = (0.E0,0.E0) 
  endif 

end function get_averaged_apar 

!------------------------------------------------------------------------------
!> This function gives the value of the averaged parallel magnetic field  
!> perturbation for a grid point. If nlbpar = .false. a zero is returned 
!------------------------------------------------------------------------------
function get_averaged_bpar(imod,ix,i,j,is,fdis) 

  use control,        only : nlbpar, spectral_radius 
  use index_function, only : indx 
  use functions,      only : mod_besselj1_gkw 
  use dist,           only : ibpar, ibpar_ga

  integer, intent(in)  :: imod, ix, i, j, is
  complex, intent(in)  :: fdis(:)  
  complex :: get_averaged_bpar 

  if (nlbpar) then 
    if (spectral_radius) then 
      get_averaged_bpar = mod_besselj1_gkw(imod,ix,i,j,is) * &
                        & fdis(indx(ibpar,imod,ix,i)) 
    else 
      get_averaged_bpar = fdis(indx(ibpar_ga,imod,ix,i,j,is)) 
    endif 
  else 
    get_averaged_bpar = (0.E0,0.E0) 
  endif 


end function get_averaged_bpar 

!---------------------------------------------------------------------
!> This routine caclulates the integrated quantities of the perturbed
!> distribution function: the potential and the vector potential 
!> 
!> It is an optimized routine, meaning that explicit assumptions are 
!> made on the the memory layout of the code. 
!>
!>. calls MPI_ALL_REDUCE to sum and distribute the integrals to all 
!>  processors. 
!---------------------------------------------------------------------
subroutine calculate_fields(fdis,DPART)  ! OPTIMISED routine 

  use dist,           only : nsolc, n_phi_start, nregular_fields_end, iapar
  use matdat,         only : maty,matz
  use matdat, only : mat_poisson, mat_field_diag
  use control,        only : nlphi, nlapar, nlbpar, zonal_adiabatic, time
  use grid,           only : nx, ns, nmod, parallel_s
  use mpiinterface,   only : number_of_processors, mpiallreduce_sum
  use mpiinterface,   only : mpiallreduce_sum_inplace
  use mpicomms,       only : COMM_S_EQ, COMM_S_NE
  use components,     only : adiabatic_electrons, tearingmode, tear_zero_epar
  use components,     only : isl_shear, isl_ls, wstar, tm_start, tm_sat
  use tearingmodes,   only : islandstruct,islandindx, omega_rot
  use tearingmodes,   only : isl_phi_indx, isl_phi, imodisland
  use constants,      only : ci1, pi
  use mode,           only : krho, kxrh, lshat_zero, ixzero
  use index_function, only : indx
  use perform,        only : perfon, perfoff, perffields, perf_measure
  use perform,        only : perfswitch
  use matrix_format, only : usmv

  complex, intent(inout) :: fdis(nsolc)

  integer :: ix, i, nelem, k
  integer :: ierr

  real :: phase, tm_fac = 0.0, dumt
  real, optional, intent(in) :: DPART

  ! warning #1 here an assumption of the postion of phi in 
  ! the solution is made. I.e. the index function is not 
  ! used. This is of course faster.
!APS: This could be replaced with another appropriate integer denoting the
!APS: index of the first field element.

  !Initialise to zero (required if fields are not solved: tearingmodes, tracers)
  !$omp parallel do schedule(static) 
  zero : do i = n_phi_start, nregular_fields_end
    fdis(i) = (0.E0,0.E0)
  end do zero
  !$omp end parallel do

  ! no point in calculating the potential if it is to 
  ! be zero. (important also for the testcases that do not use phi)
  if ((.not. nlphi) .and. (.not.nlapar) .and. (.not.nlbpar)) then
    return
  end if

  if (perffields .and. perf_measure) then
    call perfon('fields - poisson mat multiply',2)
  end if
  
  ! first calculate the contribution of f 
  ! NOTE: Calling this with x=fdis and y=fdis is only possible because
  ! in y the gyoraveraged fields are updated whereas in x the
  ! not-gyroaveraged fields are read.
  call usmv(cmplx(1.0,0.0), mat_poisson,fdis,ierr)

  if (perffields .and. perf_measure) call perfswitch('fields - MPI allreduce',2)

  if (number_of_processors > 1) then
    ! finish the species sum and the velspace integral
    call mpiallreduce_sum_inplace(fdis(n_phi_start:nregular_fields_end),&
       nregular_fields_end-n_phi_start+1,COMM=COMM_S_EQ)
  end if

  ! if the zonal flow correction is present 
  zonal_adiabat_correction : if (adiabatic_electrons .and. zonal_adiabatic) then
    if (perffields .and. perf_measure) call perfswitch('fields - zonal adiabatic',2)

    !Calculate {(1/A)poisson_int} since poisson_int is now in fdis
    bufphi = 0.0
    call usmv(cmplx(1.0,0.0), matz,fdis,bufphi,ierr)

    ! Finish off the flux surface averaging
    if (parallel_s) then
      call mpiallreduce_sum_inplace(bufphi(1:nx),nx,COMM_S_NE)
    end if

    call usmv(cmplx(1.0,0.0), maty, bufphi, fdis, ierr)
  end if zonal_adiabat_correction
  
  if (perffields .and. perf_measure) call perfswitch('fields - n4 renorm diag',2)

  ! Then divide by A: division by fields diagonal part, all elements
  ! are independent. Since we know that mat_field_diag has a unity
  ! block in the distribution part of the matrix, we can skip that and
  ! begin at n_phi.
  !!$omp parallel do schedule(static)
  normalise : do i = n_phi_start, mat_field_diag%nmat
  !normalise : do i = 1, mat_field_diag%nmat
    fdis(mat_field_diag%ii(i)) = fdis(mat_field_diag%ii(i))*mat_field_diag%mat(i)
    !fdis(mat_field_diag%ii(i)) = fdis(mat_field_diag%jj(i))*mat_field_diag%mat(i)
  end do normalise
  !!$omp end parallel do

  tearingmod : if(tearingmode) then
    if (perffields .and. perf_measure) call perfswitch('fields - tearingmod block',2)

    if (present(DPART)) then
      phase = (time + DPART)*omega_rot
      dumt = time + DPART
    else
      phase = time*omega_rot
      dumt = time
    end if

    ! Re-apply the perturbation that is the magnetic island.
    if (dumt > tm_start) then
      ! Grow the mode gradually:
      if (dumt < tm_sat) then
        ! the prefactor tm_fac follows a linear slope for dumt < tm_sat
        tm_fac = (dumt - tm_start)/(tm_sat-tm_start)
        tm_fac = min(tm_fac,1.0)
        tm_fac = max(tm_fac,0.0)
      else
        tm_fac = 1.0
      end if

!!$omp parallel do private(i,ix) schedule(static)
      do i=1,ns
        do ix=1,nx
          nelem = islandindx(ix,i)
          fdis(nelem)=fdis(nelem)+islandstruct(ix,i)*exp(-ci1*phase)*tm_fac
        end do
      end do
!!$omp end parallel do
    end if

    ! Tearing modes phi
    ! Re-apply the the island phi perturbation
    tear2 : if (tear_zero_epar .and. time > tm_start) then
      if (perffields .and. perf_measure) call perfswitch('fields - tear2 block',2)
      
!!$omp parallel do private(i,ix) schedule(static)
      do i=1,ns
        do k=1,nmod; do ix=1,nx
          nelem = isl_phi_indx(k,ix,i)
          fdis(nelem)=fdis(nelem)+isl_phi(k,ix,i)*exp(-ci1*phase*krho(k)*tm_fac/krho(imodisland))
        end do; end do
      end do
!!$omp end parallel do
    end if tear2

    ! Tearing modes shearless case only
    ! Re-apply the perturbation that is the island in ky=0 mode
    tear3 : if (lshat_zero .and. time > tm_start) then
      if (perffields .and. perf_measure) call perfswitch('fields - tear3 block',2)
      
!!$omp parallel do private(i,ix) schedule(static)
      do i=1,ns
        do ix=1,nx
          if (ix == ixzero) cycle

          if (wstar > 0) then
            !fdis(indx(iapar,1,ix,i)) = (-1)**(ix-ixzero)*4*isl_shear*tm_fac*&
            ! & exp(-((ix-ixzero)/isl_ls)**2)/kxrh(ix)**2
            fdis(indx(iapar,1,ix,i)) = (-1)**(ix-ixzero+1) * isl_shear * tm_fac &
               & * exp(-((ix-ixzero)/isl_ls)**2) / (kxrh(ix)**2)
          else ! 2 periodic islands
            if (ix - ixzero == 1) then       
              fdis(indx(iapar,1,ix,i)) =  pi * isl_shear * tm_fac / (kxrh(ix)**2)
            else if (ix - ixzero == -1) then
              fdis(indx(iapar,1,ix,i)) =  pi * isl_shear * tm_fac / (kxrh(ix)**2)
            end if
          end if

        end do
      end do
!!$omp end parallel do
    end if tear3
  end if tearingmod

  if (perffields .and. perf_measure)  call perfoff(2) 

end subroutine calculate_fields


!----------------------------------------------------------------
!> wrapper routine to do MPI for nonspec field solve (if needed)
!> two gyro-averages in the calculation (one fdis, once iphi)
!> imply that two sets of communications of radial ghost
!> points will always be required
!----------------------------------------------------------------
subroutine field_solve_nonspec_wrap(fdis,dpart,start_gc2x)

  use dist,         only : fdis_tmp2, nsolc 
  use grid,         only : lsendrecv_x
  use mpighosts,    only : mpistart, gc2x2_f
  
  complex, intent(inout) :: fdis(nsolc)
  real,    intent(in)    :: dpart
  !> start communication of x ghost cells before solving
  logical, intent(in) :: start_gc2x
  
  if (lsendrecv_x) then
    ! Copy (maybe this can be removed later)
    ! These excess memory copies can take up to 2% of total non MPI time   
    fdis_tmp2(1:nsolc)=fdis(1:nsolc)
    
    ! start communication of x ghost cells (f only) into fdis_tmp2
    ! this will be waited for in calculate_rhs
    if (start_gc2x) call mpistart(gc2x2_f)

    call field_solve_nonspec(fdis_tmp2,dpart)
   
    ! copy back again - can optimise later
    ! These excess memory copies can take up to 2% of total non MPI time   
    fdis(1:nsolc)=fdis_tmp2(1:nsolc) 
 
  else
    call field_solve_nonspec(fdis,dpart)
    
  end if

end subroutine

!------------------------------------------------------------------------------
!> Routine that calculates the nonspectral fields 
!> using the matrices calculated and stored in gyro_average
!------------------------------------------------------------------------------
subroutine field_solve_nonspec(fdis,dpart)
  use general,        only : gkw_abort
#if defined(umfpack)
  use control,        only : nlapar, shear_periodic, time
  use dist,           only : msolc, nsolc, n_phi_start
  use grid,           only : nmod, nx, ns, iproc_x, n_x_grid, lx
  use grid,           only : parallel_sp, parallel_vpar, parallel_mu, gx
  use grid,           only : gathv_mod_disps, gathv_mod_recvno, lsendrecv_x
  use mpiinterface,   only : mpiallreduce_sum, MPIREAL_X
  use mpicomms,       only : COMM_S_EQ_X_EQ, COMM_X_NE, COMM_MOD
  use mpighosts,      only : mpistart, mpiwait, gc2x_phi, gc2x_pga
  use mpidatatypes,   only : TYPE_MOD_SEND, TYPE_MOD_RECV
  use dist,           only : iphi_ga, iphi, iapar_ga, iapar
  use rotation,       only : shear_periodic_bcs
  use constants,      only : ci1
  use gyro_average,   only : matred, matred_s
  use gyro_average,   only : matint, matpol
  use gyro_average,   only : iilocz, iiscatz, nmod_l, n_phi, nscatz
  use gyro_average,   only : numeric, control_umf, ihavemodes, info
  use gyro_average,   only : phi_ga_start, phi_ga_end, apar_ga_start, apar_ga_end
  use gyro_average,   only : polarization_init, first_imod
  use rho_par_switch, only : s_average
  use components,     only : isl_shear, wstar, tm_start, tm_sat, tearingmode
  use tearingmodes,   only : islandstruct,islandindx, omega_rot
  use constants,      only : ci1
  use geom,           only : bmin
  use mode,           only : lshat_zero
  use index_function, only : indx
  use global,         only : iumf 
  use perform,        only : perfon, perfoff, perffields, perf_measure, perfswitch
  use linear_terms,   only : lpoisson
  use matrix_format, only : usmv

  complex, intent(inout) :: fdis(:)
  real, intent(in)       :: dpart  
  real                   :: phase, tm_fac, dumt, ix_mid, psi_tmp
  integer   :: i, nelem, ixoffset, ierr, nelem_l
  integer(kind=iumf) :: sys
  integer :: ix, ixstart
    
  ! checks on the size of passed fdis (could be removed)
  if (lsendrecv_x) then
    if (size(fdis) /= msolc) call gkw_abort('field_solve wrong size fdis') 
  else
    if (size(fdis) /= nsolc) call gkw_abort('field_solve wrong fdis size') 
  end if

  if (perffields .and. perf_measure) call perfon('Fields Int',2)

  ! zero all parts of the ga_array about to be used
  fdis(phi_ga_start:phi_ga_end) = (0.E0,0.E0)
  if (nlapar) fdis(apar_ga_start:apar_ga_end) = (0.E0,0.E0)

  ! do the local integral of f over the parallel velocity
  ! temporarily stored in the gyro-average field
  call usmv(cmplx(1.0), matint, fdis, fdis, ierr)

  if (perffields .and. perf_measure) call perfswitch('Fields Int radial comms',2)

  ! Communicate the x ghost cells into fdis_tmp2 (only ga-phi is needed)
  ! assuming this routine was called by reference with fdis_tmp2
  call mpistart(gc2x_pga)
  ! Wait for the communication to finish.
  call mpiwait(gc2x_pga)
  
  ! shift the ghost cells
  call shear_periodic_bcs(fdis,iphi_ga,dpart)
  call shear_periodic_bcs(fdis,iapar_ga,dpart)
  
  if (perffields .and. perf_measure) call perfswitch('Fields Int GA',2)

  ! gyro-average and do the integral over mu and the summation over is 
  ! jj(iphi_ga) accesses the ghosts, ii(iphi) does not
  z(:) = 0. 
  do i = 1, matred%nmat
    z(matred%ii(i)) = z(matred%ii(i)) + real(matred%mat(i)*fdis(matred%jj(i)))
    z(matred%ii(i)+1) = z(matred%ii(i)+1) + aimag(matred%mat(i)*fdis(matred%jj(i)))
  end do 
  
  if (perffields .and. perf_measure) call perfswitch('Fields Int MPIallreduce',2)

  ! Global sum 
  if (parallel_mu .or. parallel_sp .or. parallel_vpar) then 

    ! add the elements 
    nelem = 2*nmod*nx*ns
    if (nlapar) nelem = 2*nelem
    call mpiallreduce_sum(z(1:nelem),z_buf(1:nelem),nelem,COMM_S_EQ_X_EQ)

  if (perffields .and. perf_measure) call perfswitch('Fields nmod copy / scatter',2)

  ! copy back (local values of nmod only) 
    do i = 1, nscatz 
      z(iilocz(i))   = z_buf(iiscatz(i))
      z(iilocz(i)+1) = z_buf(iiscatz(i)+1)
    ! z is now  local in mod global in ix, local in i
    end do
    
  endif

  if (perffields .and. perf_measure) call perfswitch('Fields x gather',2)
  
  sys = 0
  ! This does the polarization, using a matrix stored inside umfpack
  ! which is referenced by the pointer numeric
  ! the result is stored in the non gyro-averaged field iphi
  nelem = 2*nmod_l*n_x_grid*ns
  nelem_l = 2*nmod_l*nx*ns

  phi_dum(:) = 0.0

  if (lsendrecv_x) then
#if defined(mpi2)

    ! This gather over x points relies on index_order indexing last in ix
    ! direct MPI faster than gather_array which allocates temporary arrays
    ! argument order REVERSED from gather_array, size arguments differ 
    call MPI_ALLGATHER(z,     nelem_l,  MPIREAL_X, &
                  &    z_buf, nelem_l,  MPIREAL_X, COMM_X_NE, ierr)
                  
    if (nlapar) then
      call MPI_ALLGATHER(z(nelem_l+1), nelem_l,  MPIREAL_X, &
                  &    z_buf(nelem+1), nelem_l,  MPIREAL_X, COMM_X_NE, ierr)
    end if
    
    !call gather_array(z_buf(1:nelem),         nelem,          &
    !              &       z(1:nelem_l ),      nelem_l,        &
    !                      COMM_X_NE, ALLGATHER = ALL_PROCS)
#else
    if (nlapar) nelem = 2*nelem
    z_buf(1:nelem) =  z(1:nelem)
#endif

  if (perffields .and. perf_measure) call perfswitch('Fields Polarize umfpack',2)
    
    ! rebuild the polarisation matrix when boundaries are time dependant
    ! VERY SLOW - not optimised !
    if (shear_periodic) call polarization_init(dpart)

    ! global polarisation solve in x, local in y, using matrix stored in umfpack 
    if (ihavemodes) call umf4sol (sys, phi_dum(1), z_buf(1), numeric, control_umf, info)

  else  !no lsendrecvx 

    ! uses phi_dum for possible double -> single precision conversion from umfpack
    if (ihavemodes) call umf4sol (sys, phi_dum(1), z(1), numeric, control_umf, info)
    ! Solves Mat(numeric)*fdis(n_phi_start) = z  
    ! breaks with runtime seqgfault in single precision: type mismatch
    !if (ihavemodes) call umf4sol (sys, fdis(n_phi_start), z(1), numeric, control_umf, info)

  end if !lsendrecvx

  ! This copy is fast, removing it does not gain much
  ! "scatter" back to local x, also relies on index order
  ixoffset = iproc_x*nx*ns*nmod_l - n_phi_start + 1
  do i = n_phi_start,n_phi_start + nmod_l*nx*ns -1
    ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION: DO NOT REMOVE THIS COPY
    if (lpoisson) then
      fdis(i) = phi_dum(i+ixoffset) 
    else
      fdis(i) = (0.E0,0.E0)
    end if
  end do
  
  ! possibly not the cleanest way to do this...
  if (nlapar) then
    ixoffset = iproc_x*nx*ns*nmod_l + n_x_grid*ns*nmod_l - nmod_l*nx*ns - n_phi_start + 1
    do i = n_phi_start + nmod_l*nx*ns, n_phi_start + 2*nmod_l*nx*ns -1
      ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION: DO NOT REMOVE THIS COPY
      fdis(i) = phi_dum(i+ixoffset) 
    end do
  end if

  if (info (1) .lt. 0) then
    print *, 'Error occurred in umf4sol: ', info (1)
    call gkw_abort('Problem in gyro_average : inversion of phi')
  endif
  
  if (perffields .and. perf_measure) call perfswitch('Fields nmod gather',2)

    ! solve polarisation order rhostar
    ! the gyro av and integration over the distribution is done local in ix  and global in mod
    ! the polarisation is done global in ix and local in mod
    sav:if (s_average) then
      phi1=0
      z(:)=0
      ! integrate over the velocity and calculate the gyroaverage along the s-direction of the distribution function
      ! do z=\bar f

      do i = 1, matred_s%nmat
        z(matred_s%ii(i)) = z(matred_s%ii(i)) + real(matred_s%mat(i)*fdis(matred_s%jj(i)))
        z(matred_s%ii(i)+1) = z(matred_s%ii(i)+1) + aimag(matred_s%mat(i)*fdis(matred_s%jj(i)))
      end do 

      nelem = 2*nmod_l*n_x_grid*ns
      nelem_l = 2*nmod_l*nx*ns
      if (lsendrecv_x) then
        call gkw_abort('parallel x and s_average not done jet')
#if defined(mpi2)
        ! This gather over x points relies on index_order indexing last in ix
        ! direct MPI faster than gather_array which allocates temporary arrays
        ! argument order REVERSED from gather_array, size arguments differ 
        call MPI_ALLGATHER(z,     nelem_l,  MPIREAL_X, &
                      &    z_buf, nelem_l,  MPIREAL_X, COMM_X_NE, ierr)
                      
        if (nlapar) then
          call MPI_ALLGATHER(z(nelem_l+1), nelem_l,  MPIREAL_X, &
                      &    z_buf(nelem+1), nelem_l,  MPIREAL_X, COMM_X_NE, ierr)
        end if
#else
        if (nlapar) nelem = 2*nelem
        z_buf(1:nelem) =  z(1:nelem)
#endif

        ! z is now local in mod and global in ix
        ! but we haven't finished yet: the contribution of potential O(0) to the
        ! equation is still missing.
#ifdef DEBUG
        if (any(matpol%jj(1:matpol%nmat) > nelem)) then
          call gkw_abort("fields::fields_solve_nonspec:  error in matpol%jj")
        endif
#endif
        do i=1,matpol%nmat
          z_buf(matpol%ii(i))      = z_buf(matpol%ii(i)  ) - real (matpol%mat(i)*phi_dum(matpol%jj(i)))
          z_buf(matpol%ii(i)+1)    = z_buf(matpol%ii(i)+1) - aimag(matpol%mat(i)*phi_dum(matpol%jj(i)))
        enddo

        nelem = 2*nmod*nx*ns
        if (nlapar) nelem = 2*nelem
        call mpiallreduce_sum(z_buf(1:nelem),z(1:nelem),nelem,COMM_S_EQ_X_EQ)      
        do i = 1, nscatz 
          z_buf(iilocz(i))   = z(iiscatz(i))
          z_buf(iilocz(i)+1) = z(iiscatz(i)+1)
        end do


        ! solve the equation for the potential O(rhostar)
        if (ihavemodes) call umf4sol (sys, phi1(1), z_buf(1), numeric, control_umf, info)

        if (info (1) .lt. 0) then
          print *, 'error occurred in umf4sol: ', info (1)
          call gkw_abort('problem in gyro_average')
        endif


        ! "scatter" back to local x, also relies on index order
        ixoffset = iproc_x*nx*ns*nmod_l - n_phi_start + 1
        do i = n_phi_start,n_phi_start + nmod_l*nx*ns -1
          ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION: DO NOT REMOVE THIS COPY
          if (lpoisson) then
            fdis(i) = fdis(i)+phi1(2*(i+ixoffset)-1)+ci1*phi1(2*(i+ixoffset))
          else
            fdis(i) = (0.E0,0.E0)
          end if
        end do
        
        ! possibly not the cleanest way to do this...
        if (nlapar) then
          ixoffset = iproc_x*nx*ns + n_x_grid*ns*nmod_l - nmod_l*nx*ns - n_phi_start + 1
          do i = n_phi_start + nx*ns, n_phi_start + nmod_l*nx*ns -1
            ! WARNING: IMPLICIT UMFPACK TYPE CONVERSION: DO NOT REMOVE THIS COPY
            fdis(i) = fdis(i)+phi1(2*(i+ixoffset)-1)+ci1*phi1(2*(i+ixoffset))
          end do
        end if


      else !lsendrecv_x
        ! since the umfpack solver can only deal with real numbers, we
        ! flatten our complex field into buffer z.
        do i=1,matpol%nmat
          ! matpol%ii refers to the 'global' first_imod position.
          z(matpol%ii(i))      = z(matpol%ii(i))- real(matpol%mat(i)*fdis(matpol%jj(i)))
          z(matpol%ii(i)+1)    = z(matpol%ii(i)+1) - aimag(matpol%mat(i)*fdis(matpol%jj(i)))
        enddo
        if (parallel_mu .or. parallel_sp .or. parallel_vpar ) then !parallelize over the mods
          ! add the elements 
          nelem = 2*nmod_l*nx*ns
          if (nlapar) nelem = 2*nelem
          call mpiallreduce_sum(z(1:nelem),z_buf(1:nelem),nelem,COMM_S_EQ_X_EQ)      
          do i = 1, nscatz 
            z(iilocz(i))   = z_buf(iiscatz(i))
            z(iilocz(i)+1) = z_buf(iiscatz(i)+1)
          end do
        endif
        if (ihavemodes) call umf4sol (sys, phi1, z(1), numeric, control_umf, info)
        ixoffset = -n_phi_start+1 
        ixstart=n_phi_start+(first_imod-1)*nx*ns-1
        do i=ixstart+1,ixstart+ nmod_l*nx*ns
          if (lpoisson) then
            fdis(i) = fdis(i)+phi1(2*(i+ixoffset)-1)+ci1*phi1(2*(i+ixoffset))
          else
            fdis(i) = (0.E0,0.E0)
          end if
        enddo
!        if ( nlapar ) then
!          call gkw_abort( ' not implemented')
!          ixoffset = -indx(iapar,1,1,1) + 1 
!          do i=ixstart,ixstart+ nmod_l*nx*ns
!            fdis(i) = fdis(i)+phi1(2*(i+ixoffset)-1)+ci1*phi1(2*(i+ixoffset))
!          enddo
!        endif
      endif
    endif sav  ! s_average

    ! distribute the local nmod back to global locations
    ! then reduce to get all modes on all procs (mpi scatter would be better)
    ! Maybe it can be done faster with derived datatypes and mpi scatter
#if defined(mpi2)
    
    ! copy: could also be combined with the "x scatter"
    nelem = nx*ns*nmod_l
    if (nlapar) nelem = 2*nelem  
    do i = 1, nelem
      bufphi(i) = fdis(i+n_phi)
    end do
 
    ! Gather the distributed mode solutions back to all procs
    ! If this line causes a segfault, look at issue 137. 
    call MPI_ALLGATHERV(bufphi(1),1*nmod_l,TYPE_MOD_SEND,fdis(n_phi+1), &
       &      gathv_mod_recvno(1),gathv_mod_disps(1),TYPE_MOD_RECV,COMM_MOD,ierr)
    
#endif            
   
    if (tearingmode) then

       phase = time*omega_rot
       dumt = time
  
       ix_mid = real((n_x_grid+1)*0.5E0)
       ! Grow the mode gradually
       if (dumt < tm_sat) then
         tm_fac = (dumt - tm_start)/(tm_sat-tm_start)
         tm_fac = min(tm_fac,1.0)
         tm_fac = max(tm_fac,0.0)
       else
          tm_fac = 1.0
       end if
       do i=1,ns
         do ix=1,nx
           nelem = islandindx(ix,i)
           fdis(nelem)=islandstruct(ix,i)*exp(-ci1*phase)*tm_fac
           if (lshat_zero) then
             if (wstar.gt.0) then
               psi_tmp = lx*(real(gx(ix))-ix_mid)/n_x_grid
               fdis(indx(iapar,1,ix,i)) = &
                      & -tm_fac*0.5E0*isl_shear*bmin(ix)* &
                      & (psi_tmp**2-0.25E0*lx**2)* &
                      !& 0.25E0*(1+tanh((psi_tmp+psi_0)/delta_psi_0))* &
                      !& (1-tanh((psi_tmp-psi_0)/delta_psi_0))* &
                      & (1.E0,0.E0)
             else ! 2 periodic islands
               call gkw_abort('Tearing mode: 2 periodic &
&                        islands not implemented in the &
&                        non spectral case.')
             end if
           end if
         end do
       end do      
    end if
    
  if (perffields .and. perf_measure) call perfswitch('Fields solved radial comms',2)
    
    ! Communicate the x ghost cells into fdis_tmp2 (only iphi is needed).
    call mpistart(gc2x_phi)
    ! Wait for the communication to finish.
    call mpiwait(gc2x_phi)
    
    call shear_periodic_bcs(fdis,iphi,dpart)
    call shear_periodic_bcs(fdis,iapar,dpart)

    if (perffields .and. perf_measure) call perfoff(2)
    if (perffields .and. perf_measure) call perfon('Fields gyroaverage',2)
         
    ! calculate the gyro_average field from the non-gyro-averaged field
    ! jj acceses the ghosts in iphi to give ii in iphi_ga (no ghosts)
    call gyro_average_fields(fdis)

    if (perffields .and. perf_measure) call perfoff(2)
     
#else
  complex, intent(inout) :: fdis(:)
  real, intent(in)       :: dpart
  call gkw_abort('Non spectral presently requires compilation with umfpack')
#endif 

end subroutine field_solve_nonspec


!------------------------------------------------------------------------------
!> This routine performs the gyro-average of the potential, parallel vector
!> potential and perturbed parallel magnetic field.
!> input  fdis : contains the fields (phi, apar, bpar)
!> output fdis : the fields that contain the gyro-averages of the field are
!>               recomputed (phi_ga, apar_ga, bpar_ga)
!------------------------------------------------------------------------------
subroutine gyro_average_fields(fdis) 
  use control,      only : nlapar
  use control,      only : nlapar
  use gyro_average, only : matav
  use gyro_average, only : phi_ga_start, phi_ga_end, apar_ga_start, apar_ga_end
  use matrix_format, only : usmv
  !> size can be either msolc or nsolc
  complex, intent (inout) :: fdis(:)
  integer :: ierr

  fdis(phi_ga_start:phi_ga_end) = (0.E0,0.E0)
  if (nlapar) fdis(apar_ga_start:apar_ga_end) = (0.E0,0.E0)

  ! This accesses ghostpoints of the x-argument array.
  ! NOTE: Calling this with x=fdis and y=fdis is only possible because
  ! in y the gyoraveraged fields are updated whereas in x the
  ! not-gyroaveraged fields are read.
  call usmv(cmplx(1.0), matav, fdis, fdis, ierr)
end subroutine gyro_average_fields 


!-----------------------------------------------------------------------------
!> Calculate the collisions conservation fields (should be optimised routine)
!> Lacks OpenMP.  
!-----------------------------------------------------------------------------
subroutine calculate_cons(fdis)

    use matdat,       only : matm
    use dist,         only : nelem_conserve, n_conserve, nsolc
    use mpiinterface, only : number_of_processors, mpiallreduce_sum_inplace
    use mpicomms,     only : COMM_VPAR_NE_MU_NE
    use perform,      only : perfon, perfoff
    use matrix_format, only : usmv
    
    complex, intent(inout) :: fdis(nsolc)
    integer :: i, ierr
   
    if (matm%nmat < 1) return
    
    call perfon('Coll cons: local',2)
 
    zero : do i = n_conserve+1,n_conserve+nelem_conserve
      fdis(i) = (0.E0,0.E0)
    end do zero    
    
    ! first calculate the contribution of f.
    ! NOTE: Calling this with x=fdis and y=fdis is only possible because
    ! in y the gyoraveraged fields are updated whereas in x the
    ! not-gyroaveraged fields are read.
    call usmv(cmplx(1.0), matm, fdis, fdis, ierr)
    
    call perfoff(2)
    call perfon('Coll cons: MPI allreduce',2)

    ! reduces fields over same species and s points
    if (number_of_processors > 1) then
      call mpiallreduce_sum_inplace(                              &
          & fdis(n_conserve+1:n_conserve+nelem_conserve),         &
          & nelem_conserve,COMM=COMM_VPAR_NE_MU_NE)
      
      !call mpiallreduce_sum(                                       &
      !    & fdis(n_conserve+1:n_conserve+nelem_conserve),          &
      !    & bufphi(1:nelem_conserve),                              &
      !    & nelem_conserve,COMM=COMM_VPAR_NE_MU_NE)

    
      !!$OMP PARALLEL DO
      !do i = n_conserve+1,n_conserve+nelem_conserve
      !   fdis(i) = bufphi(i-n_conserve)
      !end do
      !!$OMP END PARALLEL DO
    end if
    
    call perfoff(2)
   
end subroutine calculate_cons


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Allocate arrays for the module
!-----------------------------------------------------------------------
subroutine fields_allocate

  use control,        only : nlapar, spectral_radius
  use grid,           only : nx,ns, n_x_grid, nmod, ns
  use dist,           only : n_phi_start, nsolc
  use gyro_average,   only : n_s_ga
  use rho_par_switch, only : s_average 
  use general,        only : gkw_abort
  
  integer :: ierr, nelem 

  ! allocate the help array for the use of mpi
  if (nsolc-n_phi_start+1 < nx*2) then
    allocate(bufphi(2*nx), stat = ierr)
  else
    allocate(bufphi(nsolc-n_phi_start+1), stat = ierr)
  end if
  if (ierr.ne.0) call gkw_abort('Could not allocate bufphi in fields')

  if (.not.spectral_radius) then

    nelem = 2*nmod*n_x_grid*ns 
    if (nlapar) nelem = 2*nelem
    if(s_average) nelem=nelem
    
    allocate(z(nelem), stat = ierr) 
    if (ierr /= 0) call gkw_abort('can not allocate z')

    allocate(z_buf(nelem), stat = ierr) 
    if (ierr /= 0) call gkw_abort('can not allocate z_buf')

    nelem = nmod*n_x_grid*ns
    if (nlapar) nelem = 2*nelem

    allocate(phi_dum(nelem), stat = ierr) 
    if (ierr /= 0) call gkw_abort('can not allocate phi_dum')

  if(s_average) then
    nelem = 2*nmod*n_x_grid*ns*n_s_ga
    if (nlapar) nelem = nelem + 2*nmod*n_x_grid*ns
    allocate(phi1(nelem), stat = ierr)
    if (ierr /= 0) call gkw_abort('could not allocate phi1 in fields')
  end if
 end if
  
end subroutine fields_allocate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>Deallocate arrays for the module
!-----------------------------------------------------------------------
subroutine fields_deallocate
 
  if (allocated(bufphi)) deallocate(bufphi)
  if (allocated(z)) deallocate(z)
  if (allocated(z_buf)) deallocate(z_buf)
  if (allocated(phi1)) deallocate(phi1)

end subroutine fields_deallocate

end module fields
