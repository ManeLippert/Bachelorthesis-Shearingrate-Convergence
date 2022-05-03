!-----------------------------------------------------------------------------
!> This module calculates all the linear terms and puts them in the matrix
!> All terms I-XI except III and IX.
!> only used during the initialisation of the code.
!> The matrix is stored in matdat.  The matrix vector multiply is
!> performed during the time integration. The vector is stored in dist.
!-----------------------------------------------------------------------------
module linear_terms

  implicit none

  private 

  public :: init_linear_terms, calc_linear_terms, linear_terms_read_nml 
  public :: linear_terms_bcast_nml, linear_terms_check_params 
  public :: linear_terms_write_nml, differential_scheme, drift
  
  public :: lpoisson, lvpgrphi, lampere, lneorotsource
  public :: lvdgradf, lvpar_grad_df, lneo_equil_switch_default

  public :: dmaxwel
  
  !> switches for hard disabling of calls
  logical, parameter :: linear_term_switch_default = .true.
  logical, parameter :: lneo_equil_switch_default = .false.
  integer, dimension(128), save :: apply_on_imod = 0
  logical, save :: lvpar_grad_df          = linear_term_switch_default
  logical, save :: lvdgradf               = linear_term_switch_default
  logical, save :: ltrapdf                = linear_term_switch_default
  logical, save :: lve_grad_fm            = linear_term_switch_default
  logical, save :: lvd_grad_phi_fm        = linear_term_switch_default
  logical, save :: lvpgrphi               = linear_term_switch_default
  logical, save :: lpoisson               = linear_term_switch_default
  logical, save :: lg2f_correction        = linear_term_switch_default
  logical, save :: lpoisson_zf            = linear_term_switch_default
  logical, save :: lampere                = linear_term_switch_default
  logical, save :: lneoclassical          = linear_term_switch_default
  logical, save :: lbpar                  = linear_term_switch_default
  logical, save :: lneorotsource          = lneo_equil_switch_default
  logical, save, public :: lneo_equil     = lneo_equil_switch_default
  logical, save :: lneo_trap              = linear_term_switch_default
  logical, save :: lneo_rad               = linear_term_switch_default
  character (len = 3), save ::  neo_fsource = ''

  !> Dissipation switch idisp for parallel and v|| derivatives 
  !> 1: use the absolute value of the velocity
  !> 2: use the RMS velocity (no velocity space dependence)
  !> Positive: Use 4th derivative at 2nd order (fourth_order scheme only)
  !> Negative: Use 2nd derivative at 4th order (fourth_order scheme only)
  !>
  !> Early on it was found that option 2 was more stable for EM runs
  !> However, recently option 1 has proved to be more stable for EM runs at low ky
  !> when vpmax is large or beta is large (see issue 201)
  !> The reason is not understood. To test: is there an impact of geometry? 
  integer, save :: idisp = 2
 
  interface linear_terms_write_nml
    module procedure linear_terms_read_nml
  end interface

contains

!------------------------------------------------------------------------------
!> This subroutine reads (or writes) the linear terms namelist
!------------------------------------------------------------------------------
subroutine linear_terms_read_nml(lun,io_stat,lwrite)
  
  use io, only : write_run_parameter
  use neoequil, only : neo_equil_parse_sp_seq
  use mpiinterface, only : root_processor


  integer, intent(in)  :: lun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  namelist /linear_term_switches/ lvpar_grad_df, lvdgradf, ltrapdf, &
       & lve_grad_fm, lvd_grad_phi_fm, lvpgrphi, lpoisson, lampere,  & 
       & lg2f_correction,  lpoisson_zf, lneoclassical, lbpar,  &
       & lneorotsource, lneo_equil, lneo_trap, neo_fsource, lneo_rad, &
       & neo_equil_parse_sp_seq, &
       & idisp, apply_on_imod

  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then       
      read(lun,NML=linear_term_switches,IOSTAT=io_stat) 
    end if
  else
    if(root_processor) write(lun,NML=linear_term_switches)

    call write_run_parameter('linear_term_switches', 'lvpar_grad_df', lvpar_grad_df)
    call write_run_parameter('linear_term_switches', 'lvdgradf', lvdgradf)
    call write_run_parameter('linear_term_switches', 'ltrapdf', ltrapdf)
    call write_run_parameter('linear_term_switches', 'lve_grad_fm', lve_grad_fm)
    call write_run_parameter('linear_term_switches', 'lvd_grad_phi_fm', lvd_grad_phi_fm)
    call write_run_parameter('linear_term_switches', 'lvpgrphi', lvpgrphi)
    call write_run_parameter('linear_term_switches', 'lpoisson', lpoisson)
    call write_run_parameter('linear_term_switches', 'lampere', lampere)
    call write_run_parameter('linear_term_switches', 'lg2f_correction', lg2f_correction)
    call write_run_parameter('linear_term_switches', 'lpoisson_zf', lpoisson_zf)
    call write_run_parameter('linear_term_switches', 'lneoclassical', lneoclassical)
    call write_run_parameter('linear_term_switches', 'lbpar', lbpar)
    call write_run_parameter('linear_term_switches', 'lneorotsource', lneorotsource)
    call write_run_parameter('linear_term_switches', 'lneo_equil', lneo_equil)
    call write_run_parameter('linear_term_switches', 'neo_equil_parse_sp_seq', &
       & neo_equil_parse_sp_seq)
    call write_run_parameter('linear_term_switches', 'lneo_trap', lneo_trap)
    call write_run_parameter('linear_term_switches', 'lneo_rad', lneo_rad)  
    call write_run_parameter('linear_term_switches', 'idisp', idisp)
    call write_run_parameter('linear_term_switches', 'apply_on_imod', apply_on_imod)

  end if 

end subroutine linear_terms_read_nml

!------------------------------------------------------------------------------
!> bcast the linear terms namelist params
!------------------------------------------------------------------------------
subroutine linear_terms_bcast_nml
  use mpiinterface, only : mpibcast
  use neoequil, only : neo_equil_parse_sp_seq
  
  call mpibcast(lvpar_grad_df,   1)
  call mpibcast(lvdgradf,        1)
  call mpibcast(ltrapdf,         1)
  call mpibcast(lve_grad_fm,     1)
  call mpibcast(lvd_grad_phi_fm, 1)
  call mpibcast(lvpgrphi,        1)
  call mpibcast(lpoisson,        1)
  call mpibcast(lg2f_correction, 1)
  call mpibcast(lampere,         1)
  call mpibcast(lpoisson_zf,     1)
  call mpibcast(lneoclassical,   1)
  call mpibcast(lbpar,           1)
  call mpibcast(lneorotsource,   1)
  call mpibcast(lneo_equil,      1)
  call mpibcast(neo_equil_parse_sp_seq, size(neo_equil_parse_sp_seq))
  call mpibcast(lneo_trap,       1)
  call mpibcast(lneo_rad,        1)
  call mpibcast(neo_fsource,     3)
  call mpibcast(idisp,           1)
  call mpibcast(apply_on_imod,   size(apply_on_imod))

end subroutine linear_terms_bcast_nml

!------------------------------------------------------------------------------
!> put any checks that can be done before memory allocation in here
!------------------------------------------------------------------------------
subroutine linear_terms_check_params

  use components, only : tearingmode 
  use control,    only : zonal_adiabatic, fullf_wo_Fm, ltrapping_arakawa
  use general,    only : gkw_warn, gkw_exit
  use grid,       only : n_s_grid
  use neoequil,   only : neoequil_check_params 
 
  if (zonal_adiabatic.and.tearingmode) then
    zonal_adiabatic = .false.
    call gkw_warn('Zonal adiabatic switched off for tearingmode')
  end if

  if (lampere.and.tearingmode) then
    lampere = .false.
    call gkw_warn('Amperes law disabled with static imposed island')
  end if 

  if (n_s_grid < 2) then
    call gkw_warn('Parallel derivative terms switched off for single s point')
    lvpar_grad_df=.false.
    lvpgrphi=.false.
  end if
  
  if(lneorotsource .and. .not.lneoclassical) then
    call gkw_exit('Must be running neoclassical to run Hinton Wong source term')
  end if

  if(fullf_wo_Fm) then 
    if (lve_grad_fm.or.lvd_grad_phi_fm.or.lvpgrphi) then 
      call gkw_exit('Switch off the linear background terms when running with &
                    & fullf_wo_Fm')
    endif 
  endif

  if (ltrapping_arakawa .and. ((.not. lvpar_grad_df) .or. (.not. ltrapdf))) then
    call gkw_exit('When using arakawa scheme for terms I and IV, switching&
       & these of (lvpar_grad_df and ltrapdf) is not possible.')
  end if

  if(lneo_equil) call neoequil_check_params

end subroutine linear_terms_check_params

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine init_linear_terms
  use neoequil, only : neoequil_init

  ! set up the neoclassical distribution function correction
  ! by reading data from NEO or GKW or an analytical form.
  if(lneo_equil .or. lneo_equil_switch_default) call neoequil_init

end subroutine init_linear_terms

!------------------------------------------------------------------------------
!> This routine calls in sequence all the subroutines
!> that put the linear terms in the equation:
!> First part the linear terms of the perturbed distribution
!> Second part the Maxwell background
!> Third part the field equations
!------------------------------------------------------------------------------
subroutine calc_linear_terms

  use control,      only : dtim, nlapar, nlbpar
  use control,      only : lcollisions, vp_trap, neoclassics, disp_par, disp_vp 
  use control,      only : disp_x, disp_y, ltrapping_arakawa, uniform_mu_grid 
  use control,      only : fac_dtim_est, method, dtim_input
  use matdat,       only : finish_matrix_section, abort_if_bad_element
  use matdat,       only : mat, mat_maxwll_background, get_estimated_timestep
  use matrix_format, only : finish_matrix
  use collisionop,  only : collision_operator_setup, conservation
  use collisionop,  only : coll_mom_change_int, coll_mom_change_int_numu
  use general,      only : gkw_warn 
  use mpiinterface, only : root_processor
  use mode,         only : mode_box
  use structures,   only : matrix_element
  use source_time,  only : source_modulation
  use global,       only : gkw_a_equal_b_accuracy
 
  real    :: dtim_est, dtim_dum
  logical :: lreduced_dt = .false.

  ! To large dissipation at low velocity drives the implicit scheme unstable 
  if (method == 'IMP') idisp = 1 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! First part the linear terms of the perturbed distribution
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !Add terms I and IV as Poisson bracket
  arakawa: if (ltrapping_arakawa) then 
    call igh(disp_par,disp_vp)

  !Or add terms I and IV separately
  else 

    ! add the convection parallel to the field (Term I in the manual)
    call vpar_grd_df
    ! with v_par dissipation.  Both 2nd order and 4th order.
    call parallel_dissipation(disp_par) 

    ! The trapping term 
    if ((vp_trap.ne.1)) then 
      call dfdvp_trap  
      call dfdvp_dissipation(disp_vp)
    endif 
  
  end if arakawa

  ! add the part due to the drift in the gradient of the eikonal
  ! (Term II in the manual)
  call vdgradf

  ! Perpendicular hyperdissipation
  if (mode_box.and.( (.not. gkw_a_equal_b_accuracy(disp_x, 0.0)) &
              & .or. (.not. gkw_a_equal_b_accuracy(disp_y, 0.0)))) then
    call hyper_disp_perp(disp_x,disp_y)
  end if

  ! The collision opeator
  if (lcollisions) call collision_operator_setup

  ! the boundary region damping 
  call krook_bound

  call finish_matrix(mat)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Second part the Maxwell background
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! the ExB in the Maxwell background (Term V in the manual)
  call ve_grad_fm
  
  ! The drift in the gradient of phi times the velocity derivative of the 
  ! Maxwellian. (Term VIII in the manual)
  call vd_grad_phi_fm
  
  ! Landau damping (both 4th and 2nd order) (Term VII in the manual)
  call vpar_grd_phi
  
  call source_modulation
  
  ! Retrieve minimum timestep estimate for ALL terms (excluding fields), 
  dtim_est = get_estimated_timestep()
  
  dtim_dum = dtim_input

  ! reduce the timestep to the estimated value times user factor
  if (fac_dtim_est*dtim_est < dtim   &
    & .and. (method=='EXP' .or. method == 'EIV')) then

    dtim=min(fac_dtim_est*dtim_est,dtim_input)
    if (dtim < dtim_input) then
      call gkw_warn('Linear timestep automatically reduced')
      lreduced_dt=.true.
      dtim_input=dtim !Avoid NL estimator resetting to input value
    end if
  end if

  if (root_processor) then
    write(*,*)
    write(*,'(A,es13.5)') ' Maximum linear timestep estimate:', dtim_est
    write(*,'(A,es13.5)') ' Input timestep:                  ', dtim_dum
    write(*,'(A,es13.5)') ' Timestep in use:                 ', dtim
    write(*,*)
  end if
  
  ! store the value of nmat
  call finish_matrix(mat_maxwll_background)
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Third part: Integral (over the perturbed distribution) part of the 
! field equations
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! place a unity matrix block in the fdis->fdis part of the matrix
  !call put_unity_block('poisson_int',cmplx(1.0))
  
  ! The poisson equation
  call poisson_int
  
  ! Electro-magnetic corrections
  if (nlapar) then
    ! The correction to be added to fdisi to generate the
    ! distribution without A|| correction
    call g2f_correct

    ! The integral part of ampere's law
    call ampere_int
  endif

  if (nlbpar) call ampere_bpar_int

  ! Initialisation of the integrals required for momentum conservation
  ! Now go into a separate matrix.
  if (lcollisions.and.conservation) then
    if (uniform_mu_grid)then
      call coll_mom_change_int
    else
      call coll_mom_change_int_numu
    endif
  endif
  
  call finish_matrix_section(3)
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Fourth part: Diagonal part of the field equations 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  call put_unity_block('poisson_dia',cmplx(1.0))
  
  ! diagonal part of the poisson equation
  call poisson_dia

  ! calculate the zonal flow matrices when necessary
  call poisson_zf

  ! diagonal part of Ampere's law
  if (nlapar) call ampere_dia

  ! number of elements (n4)
  call finish_matrix_section(4)
 
  ! Finally, if neoclassical effects are to be kept, call the source routine.
  ! (Term VI in the manual)
  if (neoclassics) call neoclassical

  ! uncomment the following line to output the matrix
  !call output_matrix()
 
  ! Cleanly abort if any bad elements have been put into the matrix.
  call abort_if_bad_element
 
end subroutine calc_linear_terms

!------------------------------------------------------------------------------
!> This routine puts the motion along the field line (Term I in the manual)
!!  \f$
!! -v_R v_{parallel N} {\cal F} (d f / d s) \f$
!!
!! in the matrix. The parallel boundary conditions are implemented throug a 
!! call to connect_parallel, the differential scheme is set by the routine 
!! differential_scheme and an estimate of the critial time step is made for 
!! every species
!<-----------------------------------------------------------------------------
subroutine vpar_grd_df

  use structures,     only : matrix_element
  use matdat,         only : set_indx, add_element
  use matdat,         only : pos_par_grid, connect_parallel
  use grid,           only : nx,ns,nmu,nvpar,nsp,nmod, n_s_grid, gs
  use velocitygrid,   only : vpgr
  use geom,           only : ffun, sgr_dist
  use components,     only : vthrat, rhostar_linear
  use mode,           only : parallel_phase_shift, ixplus, ixminus, iyzero
  use dist,           only : ifdis, stencil_side, stencil_side_zf
  use global,         only : id_s
  use rotation,       only : coriolis, cf_drift
  use rho_par_switch, only : lvdgradf_rhostar
  use control,        only : lcalc_energetics
  use global,         only : gkw_a_equal_b_accuracy
  use control,        only : order_of_the_zf_scheme

  !> The integers for the loop over all grid points
  integer :: imod, ix, i, j, k ,is

  !> the variables used in the parallel boundary conditions
  integer :: ist
  logical :: ingrid

  !> the element to attempt
  type (matrix_element) :: elem

  !> dummy variables
  real    :: dum, drift_x, drift_y, drift_z
  real, allocatable :: w(:)
  integer :: ierr, ipw, id, m


  ! if only one parallel grid point (2D case) return
  if (n_s_grid .lt. 2) return

  ! Identifier of the term 
  elem%term = 'I: vpar_grad_df'
  elem%itype = ifdis 
  elem%itloc = ifdis
  elem%ideriv = 1

  ! provide a buffer large enough to hold the finite difference coefficients of a
  ! completely one-sided stencil, for the chosen order of the scheme
  allocate(w(1 + 4*max(stencil_side(id_s), stencil_side_zf(id_s))))

  do is = 1, nsp

    do imod=1,nmod
    
      if (all(apply_on_imod == 0)) then
        if(.not. lvpar_grad_df) cycle
      else
        if (any(apply_on_imod == imod)) then
          if(.not. lvpar_grad_df) cycle
        else
          if(.not.linear_term_switch_default) cycle
        endif
      endif
      
      do ix=1,nx ; do j = 1,nmu ; do k=1,nvpar ; do i=1,ns

      ! parallel velocity
      dum = -ffun(ix,i)*vthrat(is)*vpgr(i,j,k,is)
      
      if ( lvdgradf_rhostar ) then
        ! parallel derivative part of term II 
        call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
        dum = dum - rhostar_linear * drift_z 
      endif
      ! the direction of the parallel motion 
      if (dum > 0) then 
        ipw = 1 
      else 
        ipw = -1 
      end if 
      id = + 1 

      ! select the scheme 
      ist = pos_par_grid(imod,ix,i,k)
      if (imod == iyzero) then
        call differential_scheme(ist,ipw,id,w,order_of_the_zf_scheme)
      else
        call differential_scheme(ist,ipw,id,w) 
      endif
      
      
      do m = 1, size(w)
        
        if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
          call set_indx(elem,imod,ix,i,j,k,is) 
          elem%iloc = i + m - ((size(w)+1)/2) ! (max(stencil_side(id_s), stencil_side_zf(id_s)) + 1)
          elem%val  = w(m) * dum / sgr_dist  
          elem%val = elem%val * parallel_phase_shift(elem%imod, & 
            elem%ix,elem%i,elem%iloc)
          call connect_parallel(elem,ingrid)
          
          if (lcalc_energetics) then
            !------------------------------------------------------------------
            ! This part of the routine calculates the outflow through the
            ! open boundary condition at the end of a field line.
            ! Conceptionally, it should calculate the same advection term
            ! as vpar_grd_df, but on the first cell *outside* the grid. As
            ! there are no ghost cells at the global s-boundary this could
            ! maybe be done via extrapolation,
            ! resulting in an extreme one-sided upwind differentiation?
            !
            !   use values of those cells     +- to calculate the derivative 
            !                                 |                     here
            !          *     *     *    *     v
            !   ... | n-3 | n-2 | n-1 | n | (n+1)
            !
            ! Implemented is however:
            ! Assuming that the solution is rather smooth, the parallel 
            ! advection term on the last cell of the grid next to an outflow
            ! boundary is used as an approximation of the outflow.
            !------------------------------------------------------------------
            
            ! FIXME What about outflow in the nonspectral case? This is written
            ! for the spectral method and won't make sense in nonspectral?

            if ((ixplus(imod,ix) == 0 .and. ipw == 1 .and. gs(i) == n_s_grid) .or. &
              & (ixminus(imod,ix) == 0 .and. ipw == -1 .and. gs(i) == 1)) then
              
              ! Label outflow matrix elements
              ! add_element can then store it in a separate matrix 
              ! the outflow can be calculated with.
              ! See also the energetics diagnostic.
              elem%outflow = .true.
              
            end if
          end if

          if (ingrid) call add_element(elem,ierr)
          elem%outflow = .false.
        end if

      end do

    end do ; end do ; end do ; end do ;
  end do
  end do

  deallocate(w)
  
end subroutine vpar_grd_df


!------------------------------------------------------------------------------
!> This routine calculates the parallel disipation. The differential scheme 
!> is selected through the use of the routine differential_scheme. 
!>
!> Input  disp   : real, the dissipation coefficient 
!> CLEANUP: This routine replicates subroutine diffus and should be merged
!------------------------------------------------------------------------------
subroutine parallel_dissipation(disp)

  use structures,   only : matrix_element
  use matdat,       only : set_indx, add_element
  use matdat,       only : connect_parallel, pos_par_grid  
  use geom,         only : ffun, sgr_dist
  use grid,         only : nx,ns,nmu,nvpar,nsp,nmod, n_s_grid
  use velocitygrid, only : vpgr, vpgr_rms
  use components,   only : vthrat
  use general,      only : gkw_abort
  use dist,         only : ifdis, stencil_side, stencil_side_zf
  use global,       only : gkw_a_equal_b_accuracy, id_s
  use mode,         only : parallel_phase_shift, iyzero
  use control,      only : order_of_the_zf_scheme

  !> The dissipation coefficient
  real, intent(in) :: disp

  !> The integers for the loop over all grid points
  integer :: imod, ix, i, j, k ,is

  !> Variables to deal with the parallel boundary conditions 
  integer :: ist
  logical :: ingrid

  !> the element to attempt
  type (matrix_element) :: elem

  !> dummy variables
  real    :: dum, dum2
  real, allocatable :: w(:)
  integer :: id, ipw, m, ierr


  ! if only one parallel grid point (2D case) return
  if (n_s_grid .lt. 2) return

  allocate(w(1 + 4*max(stencil_side(id_s), stencil_side_zf(id_s))))

  ! Set the string to identify the term 
  elem%term  = 'parallel dissipation' 

  ! type of the elements 
  elem%itype = ifdis
  elem%itloc = ifdis
  elem%ideriv = 4

  id = + 2                ! use 4th derivative at 2nd order
  if (idisp < 0) id = -2  ! use 2nd derivative at 4th order

  do is = 1, nsp

    do imod=1,nmod

      if (all(apply_on_imod == 0)) then
        if(.not. lvpar_grad_df) cycle
      else
        if (any(apply_on_imod == imod)) then
          if(.not. lvpar_grad_df) cycle
        else
          if(.not.linear_term_switch_default) cycle
        endif
      endif

      do ix=1,nx ; do j = 1,nmu ; do k=1,nvpar ; do i=1,ns

      ! the dissipation 
      dum = -ffun(ix,i)*vthrat(is)*vpgr(i,j,k,is)
      select case(idisp)
      case(1,-1) 
        dum2 = dum 
      case(2,-2) 
        dum2  = ffun(ix,i)*vthrat(is)*vpgr_rms      
      case default 
        dum2 = 0.
        call gkw_abort('parallel_dissipation: unknown idisp')
      end select 

      ! direction of the parallel motion 
      if (dum > 0) then 
        ipw = 1 
      else 
        ipw = -1 
      end if 

      ! select the scheme 
      ist = pos_par_grid(imod,ix,i,k)
      
      if (imod == iyzero) then 
        ! Calculate the damping on the zonal mode with a higher order scheme.
        call differential_scheme(ist,ipw,id,w,order_of_the_zf_scheme)
      else
        ! returns (for ist = 0) 
        ! w = dfoc00   /-1.E0,   4.E0,  -6.E0,   4.E0,  -1.E0 /  12
        ! second order fourth derivative
        call differential_scheme(ist,ipw,id,w) 
      endif
      

      do m = 1, size(w)

        if ((.not. gkw_a_equal_b_accuracy(w(m), 0.0)) &
           & .and. ((.not. gkw_a_equal_b_accuracy(dum, 0.0)).or.(ist.eq.0))) then
          call set_indx(elem,imod,ix,i,j,k,is)
          elem%iloc = i + m - ((size(w)+1)/2)
          elem%val = w(m) * abs(dum2) * disp / sgr_dist  
          elem%val = elem%val * parallel_phase_shift(elem%imod, & 
          & elem%ix,elem%i,elem%iloc)
          call connect_parallel(elem,ingrid)
          ! if (.not. ingrid .and. zonal_flow_sixth_order_FD) then
          !   call gkw_abort('Something is wrong')
          ! endif
          if (ingrid) call add_element(elem,ierr)
        endif 

      end do 


    end do ; end do ; end do ; end do ; end do

  end do

  deallocate(w)
end subroutine parallel_dissipation

!-----------------------------------------------------------------------------
!> This routine puts the trapping (Term IV in the manual)
!!
!! \f$  + v_R \mu_N B_N {\cal G} (d f / d v_{\parallel} N)        \f$
!!
!! in the matrix. Boundary conditions are that f is zero outside the parallel
!! velocity grid - could be made to be "flappy" as in term I
!<-----------------------------------------------------------------------------
subroutine dfdvp_trap 

  use structures,     only : matrix_element
  use matdat,         only : set_indx, add_element
  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod
  use velocitygrid,   only : mugr, dvp, vpgr
  use velocitygrid, only : get_vpar_stencil, elem_is_on_vpar_grid
  use geom,           only : bn, gfun, ffun, dfun, hfun, efun
  use components,     only : vthrat, tgrid, tmp, rhostar_linear, signz, mas
  use components,     only : veta_prime
  use rotation,       only : dcfen_ds, vcor
  use dist,           only : ifdis, stencil_side, stencil_side_zf
  use global,         only : id_vpar
  use rho_par_switch, only : ltrapdf_rhostar
  use neoequil,       only : neo_nsp
  use global,         only : gkw_a_equal_b_accuracy

  character(len=64) :: term='IV: trapdf_4d'

  ! The integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is, ist

  ! element to add
  type (matrix_element) :: elem

  ! Dummy variables
  real    :: dum,  drift_z, dBdeps, dBdzeta, dBds
  real, allocatable :: w(:)
  integer :: ierr, ipw, id, m 
  logical :: ingrid

  allocate(w(1 + 4*max(stencil_side(id_vpar), stencil_side_zf(id_vpar))))

  ! The same for all elements 
  elem%term  = term 

  dBdzeta=0
  ! type of the elements 
  elem%itype = ifdis
  elem%itloc = ifdis
  elem%ideriv = 1

  do imod=1,nmod
    
    if (all(apply_on_imod == 0)) then
      if(.not. ltrapdf) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. ltrapdf) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif

    do is = 1,nsp ; do ix=1,nx ; do i=1,ns ; do j=1,nmu ; do k=1,nvpar

      call set_indx(elem,imod,ix,i,j,k,is)

      dum = vthrat(is)*mugr(j)*bn(ix,i)*gfun(ix,i)

      !Add the centrifugal correction to the trapping
      !cfen=0 if vcor=0 or cf_trap=.false.
      dum=dum+0.5*vthrat(is)*tmp(ix,is)*ffun(ix,i)*dcfen_ds(i,is)/tgrid(is)

      if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_trap ) .or. & 
      & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. linear_term_switch_default ) ) then
        if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_equil ) .or. & 
        & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
          if(is <= neo_nsp) then
            dum=dum+0.5*signz(is)*vthrat(is)*tmp(ix,is)*ffun(ix,i)* &
              & dphineods(imod,ix,i)/tgrid(is)
          endif
        endif
      endif

      ! add the rhostar term to the trapping
      if (rhostar_linear.gt.0.and. ltrapdf_rhostar) then
        dBds = gfun(ix,i)/ffun(ix,i)*bn(ix,i)
        !use dfun(ix,i,3) to calculate dBdeps
        dBdeps = dfun(ix,i,3) *bn(ix,i)/efun(ix,i,3,1)
        ! driftz_over_vpar would be more precise
        drift_z = (2E0* tmp(ix,is)*vpgr(i,j,k,is)/bn(ix,i)*  &
                    (efun(ix,i,1,1)*dBdeps*dBdeps + efun(ix,i,2,2)*dBdzeta*dBdzeta &
                      +efun(ix,i,3,3)*dBds*dBds)&
                   +(2E0-2E0) * tmp(ix,is)*vpgr(i,j,k,is)/bn(ix,i)* efun(ix,i,1,2)*dBdeps * dBdzeta &
                   +(2E0-2E0) * tmp(ix,is)*vpgr(i,j,k,is)/bn(ix,i)* efun(ix,i,1,3)*dBdeps * dBds    &
                   +(2E0-2E0) * tmp(ix,is)*vpgr(i,j,k,is)/bn(ix,i)* efun(ix,i,2,3)*dBdzeta* dBds    &
                  ) * mugr(j)
        drift_z = drift_z + tmp(ix,is)/bn(ix,i) * vpgr(i,j,k,is) * veta_prime(ix)/bn(ix,i)**2 &
                             * efun(ix,i,1,3) * dBds * mugr(j)
        drift_z = drift_z / signz(is) 
        drift_z = drift_z + 2.E0*mas(is)*vcor*( hfun(ix,i,1) *dBdeps+hfun(ix,i,3)*dBds  * mugr(j) )
        dum=dum+rhostar_linear * drift_z
      end if
      !Use function to find the value of ist which determines which 
      !stencil to use.
      ist = get_vpar_stencil(elem)
      ingrid = elem_is_on_vpar_grid(elem)
      
      if (dum > 0) then 
        ipw = 1 
      else 
        ipw = -1 
      end if 
      id = + 1 

      ! select the scheme 
      call differential_scheme(ist,ipw,id,w) 

      do m = 1, size(w)
 
        if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
          call set_indx(elem,imod,ix,i,j,k,is)
          elem%kloc = k + m - ((size(w)+1)/2)
          elem%val  = w(m) * dum / dvp
          call add_element(elem,ierr)

        endif 

      end do 

    end do ; end do ; end do ; end do ; end do

  end do

  deallocate(w)
end subroutine dfdvp_trap

!-----------------------------------------------------------------------------
!> This routine puts the dissipation on the trapping term (Term IV in the 
!! manual) in the matrix. Boundary conditions are that f is zero outside the 
!! parallel velocity grid 
!> CLEANUP: This routine replicates subroutine diffus and should be merged.
!<-----------------------------------------------------------------------------
subroutine dfdvp_dissipation(disp) 

  use structures,   only : matrix_element
  use matdat,       only : set_indx, add_element
  use grid,         only : nx, ns, nmu, nvpar, nsp, nmod
  use velocitygrid, only : mugr, dvp, mugr_rms
  use velocitygrid, only : get_vpar_stencil, elem_is_on_vpar_grid
  use geom,         only : bn, gfun, ffun
  use components,   only : vthrat, tgrid, tmp
!  use components,   only : signz
  use rotation,     only : dcfen_ds
  use dist,         only : ifdis, stencil_side, stencil_side_zf
  use global,       only : id_vpar, gkw_a_equal_b_accuracy

  ! The dissipation coefficient 
  real, intent(in) :: disp 

  ! The integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is, ist

  ! element to add
  type (matrix_element) :: elem

  ! Dummy variables
  real    :: dum, dum2
  real, allocatable :: w(:)
  integer :: ierr, ipw, id, m 
  logical :: ingrid

  allocate(w(1 + 4*max(stencil_side(id_vpar), stencil_side_zf(id_vpar))))

  ! The same for all elements 
  elem%term  = 'parallel velocity dissipation' 

  ! type of the elements 
  elem%itype = ifdis 
  elem%itloc = ifdis
  elem%ideriv = 4
  if (disp<0) elem%ideriv = 2

  id = + 2                  ! use 4th derivative at 2nd order
  if (idisp < 0) id = -2    ! use 2nd derivative at 4th order

  do imod=1,nmod
    
    if (all(apply_on_imod == 0)) then
      if(.not. ltrapdf) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. ltrapdf) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
  
    do is = 1,nsp ; do ix=1,nx ; do i=1,ns ; do j=1,nmu ; do k=1,nvpar

      call set_indx(elem,imod,ix,i,j,k,is)

      dum = vthrat(is)*mugr(j)*bn(ix,i)*gfun(ix,i)

      !Add the centrifugal correction to the trapping
      !cfen=0 if vcor=0 or cf_trap=.false.
      dum=dum+0.5*vthrat(is)*ffun(ix,i)*dcfen_ds(i,is)

      !if(lneo_equil.and.lneo_trap.and.is <= neo_nsp) then
      !  dum=dum+0.5*signz(is)*vthrat(is)*ffun(ix,i)*dphineods(imod,ix,i)
      !endif

      select case(idisp)
        case(1,-1)
          dum2 = dum
        case(2,-2) ! use the mugr rms value only
          dum2 = vthrat(is)*bn(ix,i)*gfun(ix,i)*mugr_rms
          dum2 = dum2+0.5*vthrat(is)*tmp(ix,is)*ffun(ix,i)*dcfen_ds(i,is) & 
               & / tgrid(is)
        case default 
          dum2 = 0. 
      end select

      !Use function to find the value of ist which determines which 
      !stencil to use.
      ist = get_vpar_stencil(elem)
      ingrid = elem_is_on_vpar_grid(elem)

      if (dum > 0) then 
        ipw = 1 
      else 
        ipw = -1 
      end if 

      ! select the scheme 
      call differential_scheme(ist,ipw,id,w) 

      do m = 1, size(w)
 
        if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then

          call set_indx(elem,imod,ix,i,j,k,is)
          elem%kloc = k + m - ((size(w)+1)/2)
          elem%val  = w(m) * disp * abs(dum2) / dvp
          call add_element(elem,ierr)

        endif 

      end do 

    end do ; end do ; end do ; end do ; end do

  end do

  deallocate(w)
  
end subroutine dfdvp_dissipation 



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Add terms I and IV as one term : \f$ v_R {\cal F} \{H,g_N\} \f$
!! where
!!         \f$ H  = (1/2) v_{\parallel N}^2 + \mu_N B_N + (1/2) {\cal E}_R = H(s,\mu,v_{\parallel}) \f$,
!! and
!!         \f$ \{H,g_N\} = (d H/d s)(d g_N/d v) - (d H/d v)(d g_N/d s). \f$
!!
!! At the boundaries in the s-direction, we apply something similar to the
!! wills=1 scheme to do upwinding. Note that the routines that this calls to
!! deal with points near the boundaries do not presently work in second order.
!!
!! See jhg_interior for the main bit of these terms.
!!
!! ( Note that above we dropped the second index with respect to `j'(for \mu)
!!  in H; the second index is always j. Below we use `HH' for `H')
!<----------------------------------------------------------------------------
subroutine igh(disp_par,disp_vp)

  use structures,   only : matrix_element
  use components,   only : vthrat
  use grid,         only : nx,ns,nmu,nvpar,nsp,nmod
  use general,      only : gkw_abort 
  use geom,         only : ffun, sgr_dist, gfun, bn
  use velocitygrid, only : dvp,vpgr,mugr,vpgr_rms,mugr_rms
  use rotation,     only : dcfen_ds
  use matdat,       only : set_indx, pos_par_grid
  use dist,         only : ifdis

  real, intent(in) :: disp_vp, disp_par
  real    :: dum, dum2, disp_v_dum, disp_s_dum
  integer :: is, imod, ix, i ,j ,k, ist

  type (matrix_element) :: elem
  
  disp_v_dum =  0.
  disp_s_dum  =  0.

  ! the term 
  elem%term  = 'I + IV: Arawaka'

  ! type of the elements 
  elem%itype = ifdis
  elem%itloc = ifdis
  elem%ideriv = 1
  
  do is = 1, nsp

    do imod=1,nmod ; do ix=1,nx ; do j = 1,nmu ; do k=1,nvpar ; do i=1,ns
            
      ! common factor of the mat elements
      dum = vthrat(is)*ffun(ix,i)/(sgr_dist*dvp)

      ! dum2 is for the advection sign and (standard) diffusion in s
      dum2 = -ffun(ix,i)*vthrat(is)*vpgr(i,j,k,is)

      ! disp_v_dum is for diffusion in vpar
      ! See subroutine trapdf_4d for analouge
      select case(idisp)
        case(1) ; disp_v_dum = vthrat(is)*(mugr(j)*bn(ix,i)*gfun(ix,i)+0.5*ffun(ix,i)*dcfen_ds(i,is))
        case(2) ; disp_v_dum = vthrat(is)*(mugr_rms*bn(ix,i)*gfun(ix,i)+0.5*ffun(ix,i)*dcfen_ds(i,is))
        case default; call gkw_abort('Arakawa does not yet implement this idisp')
      end select

      ! disp_s_dum is for diffusion in s
      ! See subroutine vpar_grad_df for analouge
      select case(idisp)
        case(1) ; disp_s_dum = dum2
        case(2) ; disp_s_dum = ffun(ix,i)*vthrat(is)*vpgr_rms
        case default; call gkw_abort('Arakawa does not yet implement this idisp')
      end select
 
      ! Check if the point is near the boundary in s; 0 denotes somewhere in
      ! the middle and +/- 1, +/- 2 correspond to the points nearest the
      ! boundary.
      !
      !      ist:  -2 -1     0    1  2                  
      ! location:   |  |  -  -  - |  | 
      !
      ist = pos_par_grid(imod,ix,i,k) 


      ! set the element 
      call set_indx(elem,imod,ix,i,j,k,is) 

      ! Add different term depending on where the point is on the grid *and*
      ! the sign of the advection velocity.
      select case(ist)
        case(-2,2)
          
          if (dum2*ist < 0.) then
            call igh_zero_two(elem,dum,ist,dum2)
          else if (dum2*ist > 0.) then
            call jhg_interior(elem,dum)
            ! apply any hyperdiffusion at second order only
            if (disp_par > 0.) call diffus(elem,disp_par,disp_s_dum,'s',2)
          else
            ! do nothing; vanishing term
          end if
          
        case(-1,1)
          
          if (dum2*ist < 0.) then
            call igh_two(elem,dum,ist,dum2)
          else if (dum2*ist > 0.) then
            call jhg_interior(elem,dum)
            ! apply any hyperdiffusion at fourth order
            if (disp_par > 0.) call diffus(elem,disp_par,disp_s_dum,'s',4)
          else
            ! do nothing; vanishing term
          end if
          
        case(0)
          
          ! original Arakawa differencing of the right order.
          call jhg_interior(elem,dum)
          
          ! apply any hyperdiffusion to the interior part of the grid, only
          ! for the s-direction (vpar-direction is applied for all cases).
          if (disp_par > 0.) call diffus(elem,disp_par,disp_s_dum,'s',4)
          
        case default
          
          call gkw_abort('igh: bad case of ist')
          
      end select
 
      ! (any) parallel velocity dissipation is applied everywhere
      if (disp_vp  > 0.) call diffus(elem,disp_vp, disp_v_dum,'vpar',4)
      
    end do ; end do ; end do ; end do ; end do

  end do

end subroutine igh


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine is called from igh in order to deal with points at the end of
!! the s grid with upwinding.
!! We difference 
!!  \f$        \{H,g_N\} = (d H/d s)(d g_N/d v) - (d H/d v)(d g_N/d s)  \f$
!! as
!!
!!  \f$      d Z(i)/d s = 1/(2 \delta s) [3 Z(i) - 4 Z(i-1) + Z(i-2)]
!!           d Z(j)/d v = 1/(2 \delta v) [Z(i+1) - Z(i-1)]              \f$ 
!!
!! For dum = vpgr*vthrat/(d v d s) and case(ist), we have the following;
!!
!! case(-2) .and. dum > 0:
!!
!!     elems = 0.25*dum*[
!!          + {  -3.*H(i,k) + 4.*H(i+1,k) - H(i+2,k) }    g(i,k+1)
!!          - {  -3.*H(i,k) + 4.*H(i+1,k) - H(i+2,k) }    g(i,k-1)
!!          - {      H(i,k+1) - H(i,k-1)             } -3*g(i,k) 
!!          - {      H(i,k+1) - H(i,k-1)             }  4*g(i+1,k)
!!          - {      H(i,k+1) - H(i,k-1)             }   -g(i+2,k)
!!                       ]
!!
!! case(-2) .and. dum < 0 .or. case(2) .and. dum > 0:
!!
!!     (apply J_1 or J_2 with regular bcs; this routine is not called)
!!
!! case(2) .and. dum < 0:
!!
!!     elems = 0.25*dum*[
!!          + {   3.*H(i,k) - 4.*H(i-1,k) + H(i-2,k) }    g(i,k+1)
!!          - {   3.*H(i,k) - 4.*H(i-1,k) + H(i-2,k) }    g(i,k-1)
!!          - {      H(i,k+1) - H(i,k-1)             }  3*g(i,k) 
!!          - {      H(i,k+1) - H(i,k-1)             } -4*g(i-1,k)
!!          - {      H(i,k+1) - H(i,k-1)             }    g(i-2,k)
!!                       ]
!<----------------------------------------------------------------------------
subroutine igh_zero_two(elem,dum,ist,dum2)

  use structures,   only : matrix_element
  use matdat,       only : add_element, connect_parallel 
  use general,      only : gkw_abort
  use mode,         only : parallel_phase_shift 

  type (matrix_element), intent(in) :: elem
  integer, intent(in) :: ist
  real, intent(in) :: dum,dum2

  type(matrix_element) :: E 
  integer :: ix, i, j, k, is, ierr
  real    :: val
  logical :: ingrid 
  real, parameter :: fac=0.25

  ix = elem%ix
  i = elem%i
  j = elem%j
  k = elem%k
  is = elem%is

  if (ist == 2 .and. dum2 < 0.) then
    
    val = fac*dum*(3.*HH(ix,i,j,k,is) - 4.*HH(ix,i-1,j,k,is) + 1.*HH(ix,i-2,j,k,is))
    
    E = elem 
    E%iloc = i
    E%kloc = k-1
    E%val  = -1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
    
    E = elem 
    E%iloc = i
    E%kloc = k+1
    E%val  =  1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    val = fac*dum*( 1.*HH(ix,i,j,k+1,is) - 1.*HH(ix,i,j,k-1,is))
    
    E = elem 
    E%iloc = i
    E%kloc = k
    E%val  = -3.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
  
    E = elem 
    E%iloc = i-1
    E%kloc = k
    E%val  =  4.* val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc) 
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i-2
    E%kloc = k
    E%val  = -1.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

  else if (ist == -2 .and. dum2 > 0) then
    
    val = fac*dum*(-3.*HH(ix,i,j,k,is) + 4.*HH(ix,i+1,j,k,is) - 1.*HH(ix,i+2,j,k,is))
    
    E = elem 
    E%iloc = i
    E%kloc = k-1
    E%val  = -1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
    
    E = elem 
    E%iloc = i
    E%kloc = k+1
    E%val  =  1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    val = fac*dum*(1.*HH(ix,i,j,k+1,is) - 1.*HH(ix,i,j,k-1,is))
    
    E = elem 
    E%iloc = i
    E%kloc = k
    E%val  = 3.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
  
    E = elem 
    E%iloc = i+1
    E%kloc = k
    E%val  = -4.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i+2
    E%kloc = k
    E%val  = 1.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

  else
    
    call gkw_abort('igh_zero_two: bad call to this routine')
    
  end if

end subroutine igh_zero_two


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This is a similar routine to igh_zero_two; differences
!>   
!>  \f$     \{H,g_N\} = (d H/d s)(d g_N/d v) - (d H/d v)(d g_N/d s)      \f$ 
!>
!> at the *second* from end point of the s direction with upwinds
!> dum = vpgr*vthrat/(d v d s).
!>
!> d Z(i)/d s = 1/(6 \delta s)  [2 Z(i+1) + 3 Z(i)   - 6 Z(i-1) + Z(i-2)]
!>      *or*    1/(6 \delta s) [- Z(i+2) + 6 Z(i+1) - 3 Z(i)   - 2 z(i-1)]
!> d Z(j)/d v = 1/(12 \delta v) [Z(i-2) - 8 Z(i-1) + 8 Z(i+1) - Z(i+2) ]
!>
!> case(-1) .and. dum > 0:
!>
!>   elems = (1/72)*dum*[
!>     + { - 2.*H(i-1,k) - 3.*H(i,k)   + 6.*H(i+1,k) - H(i+2,k) }    g(i,k-2)
!>     - { - 2.*H(i-1,k) - 3.*H(i,k)   + 6.*H(i+1,k) - H(i+2,k) }  8*g(i,k-1)
!>     + { - 2.*H(i-1,k) - 3.*H(i,k)   + 6.*H(i+1,k) - H(i+2,k) }  8*g(i,k+1)
!>     - { - 2.*H(i-1,k) - 3.*H(i,k)   + 6.*H(i+1,k) - H(i+2,k) }    g(i,k+2)
!>   -   {      H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1) - H(i,k+2) } -2*g(i-1,k)
!>   -   {      H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1) - H(i,k+2) } -3*g(i,k)
!>   -   {      H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1) - H(i,k+2) } +6*g(i+1,k)
!>   -   {      H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1) - H(i,k+2) } -1*g(i+2,k)
!>                       ]
!>
!> case(-1) .and. dum < 0 .or. case(1) .and. dum > 0:
!>
!>   (apply J_1 or J_2 with regular bcs; this routine is not called)
!>
!> case(1) .and. dum < 0:
!>
!>   elems = (1/72)*dum*[
!>     + {  H(i-2,k) - 6.*H(i-1,k) + 3.*H(i,k) + 2.*H(i+1,k) }    g(i,k-2)
!>     - {  H(i-2,k) - 6.*H(i-1,k) + 3.*H(i,k) + 2.*H(i+1,k) }  8*g(i,k-1)
!>     + {  H(i-2,k) - 6.*H(i-1,k) + 3.*H(i,k) + 2.*H(i+1,k) }  8*g(i,k+1)
!>     - {  H(i-2,k) - 6.*H(i-1,k) + 3.*H(i,k) + 2.*H(i+1,k) }    g(i,k+2)
!>   -   {  H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1)  - H(i,k+2) } +1*g(i-2,k)
!>   -   {  H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1)  - H(i,k+2) } -6*g(i-1,k)
!>   -   {  H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1)  - H(i,k+2) } +3*g(i,k)
!>   -   {  H(i,k-2) - 8.*H(i,k-1) + 8.*H(i,k+1)  - H(i,k+2) } +2*g(i+1,k)
!>                      ]
!----------------------------------------------------------------------------
subroutine igh_two(elem,dum,ist,dum2)

  use structures,   only : matrix_element
  use matdat,       only : add_element, connect_parallel 
  use general,      only : gkw_abort 
  use mode,         only : parallel_phase_shift

  type (matrix_element), intent(in) :: elem
  integer, intent(in) :: ist
  real, intent(in)    :: dum,dum2

  type (matrix_element) :: E 
  real    :: val
  integer :: ix, i, j, k, is, ierr
  logical :: ingrid 
  real, parameter :: fac=1./72.

  ix = elem%ix
  i = elem%i
  j = elem%j
  k = elem%k
  is = elem%is

  if (ist == -1 .and. dum2 > 0) then
    
    val = fac*dum*                                                            &
      & (-2.*HH(ix,i-1,j,k,is)-3.*HH(ix,i,j,k,is)+6.*HH(ix,i+1,j,k,is)-HH(ix,i+2,j,k,is))
    
    E = elem 
    E%iloc = i
    E%kloc = k-2
    E%val  = +1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
    
    E = elem 
    E%iloc = i
    E%kloc = k-1
    E%val  = -8.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i
    E%kloc = k+1
    E%val  = +8.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
  
    E = elem 
    E%iloc = i
    E%kloc = k+2
    E%val  = -1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    val =  -1.*fac*dum*                                                       &
     & (1.*HH(ix,i,j,k-2,is)-8.*HH(ix,i,j,k-1,is)+8.*HH(ix,i,j,k+1,is)-1.*HH(ix,i,j,k+2,is))
    
    E = elem 
    E%iloc = i-1
    E%kloc = k
    E%val  = -2.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i
    E%kloc = k
    E%val  = -3.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i+1
    E%kloc = k
    E%val  = +6.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i+2
    E%kloc = k
    E%val  = -1.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

  else if (ist == 1 .and. dum2 < 0) then

    val = fac*dum*                                                            &
     &  (HH(ix,i-2,j,k,is)-6.*HH(ix,i-1,j,k,is)+3.*HH(ix,i,j,k,is)+ 2.*HH(ix,i+1,j,k,is))
    
    E = elem 
    E%iloc = i
    E%kloc = k-2
    E%val  = +1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
    
    E = elem 
    E%iloc = i
    E%kloc = k-1
    E%val  = -8.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i
    E%kloc = k+1
    E%val  = +8.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
  
    E = elem 
    E%iloc = i
    E%kloc = k+2
    E%val  = -1.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    val =  -1.*fac*dum*                                                       &
     & (1.*HH(ix,i,j,k-2,is)-8.*HH(ix,i,j,k-1,is)+8.*HH(ix,i,j,k+1,is)-1.*HH(ix,i,j,k+2,is))

    E = elem 
    E%iloc = i-2
    E%kloc = k
    E%val  = +1.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i-1
    E%kloc = k
    E%val  = -6.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i
    E%kloc = k
    E%val  = +3.*val
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)

    E = elem 
    E%iloc = i+1
    E%kloc = k
    E%val  = +2.*val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
  
  else
    
    call gkw_abort('bad something wrong linear terms')
    
  end if

end subroutine igh_two


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Difference {H,g} in the interior, which includes the boundary in the vpar
!! direction.
!! For 2nd order we difference {H,g_N} = J as J = J_1 *OR* J = J_2, with
!! J_1 and J_2 as defined below.  For a 4th order scheme, we combine J_1 and
!! J_2; J = 2*J_1 - J_2 (Arakawa, JcompPhys,1,1,_119_ 1966).
!!
!! J_1 =
!!   1/3 (1 / (4 \delta s \delta v_{\par}) )
!!   *[
!!        {              H(i-1,k)   - H(i,k-1)                } g(i-1,k-1)
!!      + { H(i-1,k+1) - H(i-1,k-1) - H(i,k-1)   + H(i,k+1)   } g(i-1,k)
!!      + {              H(i,k+1)   - H(i-1,k)                } g(i-1,k+1)
!!      + { H(i-1,k-1) + H(i-1,k)   - H(i+1,k-1) - H(i+1,k)   } g(i,k-1)
!!      + { H(i+1,k)   + H(i+1,k+1) - H(i-1,k)   - H(i-1,k+1) } g(i,k+1)
!!      + {              H(i,k-1)   - H(i+1,k)                } g(i+1,k-1)
!!      + { H(i,k-1)   + H(i+1,k-1) - H(i,k+1)   - H(i+1,k+1) } g(i+1,k)
!!      + {              H(i+1,k)   - H(i,k+1)                } g(i+1,k+1)
!!    ]
!! J_2 =
!!   1/3 (1 / (8 \delta s \delta v_{\par}) )
!!   *[
!!        {            - H(i-1,k-1) + H(i-1,k+1)              } g(i-2,k)
!!      + { H(i-2,k)   + H(i-1,k+1) - H(i,k-2)   - H(i+1,k-1) } g(i-1,k-1)
!!      + { H(i+1,k+1) - H(i-1,k-1) - H(i-2,k)   + H(i,k+2)   } g(i-1,k+1)
!!      + {              H(i-1,k-1) - H(i+1,k-1)              } g(i,k-2)
!!      + {              H(i+1,k+1) - H(i-1,k+1)              } g(i,k+2)
!!      + { H(i-1,k-1) - H(i+1,k+1) - H(i+2,k)   + H(i,k-2)   } g(i+1,k-1)
!!      + { H(i+1,k-1) - H(i-1,k+1) + H(i+2,k)   - H(i,k+2)   } g(i+1,k+1)
!!      + {              H(i+1,k-1) - H(i+1,k+1)              } g(i+2,k)
!!    ]
!!
!! At the boundaries in the s-direction, we apply something similar to the
!! wills=1 scheme to do upwinding.
!!
!! (Above we drop the second index with respect to `j' (for \mu) in H; the
!!  second index is always j. Below we use `HH' for `H')
!<----------------------------------------------------------------------------

subroutine jhg_interior(elem,dum,scheme)

  use structures, only : matrix_element
  use control,    only: order_of_the_scheme
  use matdat,     only : add_element, connect_parallel 
  use general,    only : gkw_abort 
  use mode,       only : parallel_phase_shift

  type (matrix_element), intent(in) :: elem
  real, intent(in)                  :: dum

  real                  :: d1, d2
  integer               :: ix,i,j,k,is,ierr
  logical               :: ingrid 
  type (matrix_element) :: E 

  character (len=1), optional, intent(in) :: scheme
  character (len=1) :: second_order_scheme 

  if (present(scheme)) then
    second_order_scheme = scheme
  else
    second_order_scheme = "+" ! or "x"
  endif

  ix=elem%ix
  i=elem%i
  j=elem%j
  k=elem%k
  is=elem%is
  
  d1 = 0. ; d2 = 0.

  select case(order_of_the_scheme)

    case('second_order') 
      select case(second_order_scheme)
        case("+") ; d1 =  1. ; d2 =  0.
        case("x") ; d1 =  0. ; d2 =  1.
        case default ; call gkw_abort('unknown case of second_order_scheme')
      endselect

    case('fourth_order') 
      if (present(scheme)) then
        select case(second_order_scheme)
          case("+") ; d1 =  1. ; d2 =  0.
          case("x") ; d1 =  0. ; d2 =  1.
          case default ; call gkw_abort('unknown case of second_order_scheme')
        endselect
      else
        d1 =  2. ; d2 = -1.
      endif
      
    case default ; call gkw_abort('linear_terms: unknown case for J')

  end select

  d1 = d1*dum/12.
  d2 = d2*dum/24.

  ! g(i-2,k); J_2
  E = elem 
  E%iloc = i-2
  E%kloc = k
  E%val  =   d2*(       HH(ix,i-1,j,k+1,is) - HH(ix,i-1,j,k-1,is)      )
  E%val = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i-1,k-1); J_1, J_2
  E = elem 
  E%iloc = i-1
  E%kloc = k-1
  E%val  =   d1*(       HH(ix,i-1,j,k,is)   - HH(ix,i,j,k-1,is)         )    &
   &        +d2*(       HH(ix,i-2,j,k,is) + HH(ix,i-1,j,k+1,is)              &
   &                  - HH(ix,i,j,k-2,is) - HH(ix,i+1,j,k-1,is)         )
  E%val = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i-1,k); J_1
  E = elem 
  E%iloc = i-1
  E%kloc = k
  E%val  = d1*(         HH(ix,i-1,j,k+1,is) - HH(ix,i-1,j,k-1,is)            &
   &                  - HH(ix,i,j,k-1,is) + HH(ix,i,j,k+1,is))
  E%val  = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i-1,k+1); J_1, J_2
  E = elem 
  E%iloc = i-1
  E%kloc = k+1
  E%val  =   d1*(       HH(ix,i,j,k+1,is)   - HH(ix,i-1,j,k,is)          )   &
   &        +d2*(       HH(ix,i+1,j,k+1,is) - HH(ix,i-1,j,k-1,is)            &
   &                  - HH(ix,i-2,j,k,is) + HH(ix,i,j,k+2,is)            )
  E%val  = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i,k-2); J_2
  E = elem 
  E%iloc = i
  E%kloc = k-2
  E%val  =   d2*(       HH(ix,i-1,j,k-1,is) - HH(ix,i+1,j,k-1,is)   )
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i,k-1); J_1
  E = elem 
  E%iloc = i
  E%kloc = k-1
  E%val  =   d1*(      HH(ix,i-1,j,k-1,is) + HH(ix,i-1,j,k,is)          &
   &                 - HH(ix,i+1,j,k-1,is) - HH(ix,i+1,j,k,is)  )
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i,k+1); J_1
  E = elem 
  E%iloc = i
  E%kloc = k+1
  E%val  =   d1*( HH(ix,i+1,j,k,is) + HH(ix,i+1,j,k+1,is)  &
   &            - HH(ix,i-1,j,k,is) - HH(ix,i-1,j,k+1,is))
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i,k+2); J_2
  E = elem 
  E%iloc = i
  E%kloc = k+2
  E%val  =   d2*(       HH(ix,i+1,j,k+1,is) - HH(ix,i-1,j,k+1,is)   )
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i+1,k-1); J_1, J_2
  E = elem 
  E%iloc = i+1
  E%kloc = k-1
  E%val  =   d1*(       HH(ix,i,j,k-1,is)   - HH(ix,i+1,j,k,is)         )   &
   & +       d2*(       HH(ix,i-1,j,k-1,is) - HH(ix,i+1,j,k+1,is)           &
   &                  - HH(ix,i+2,j,k,is) + HH(ix,i,j,k-2,is)           )
  E%val  = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i+1,k); J_1
  E = elem 
  E%iloc = i+1
  E%kloc = k
  E%val  = d1*(          HH(ix,i,j,k-1,is)   + HH(ix,i+1,j,k-1,is)          &
   &                   - HH(ix,i,j,k+1,is) - HH(ix,i+1,j,k+1,is)     )
  E%val  = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i+1,k+1); J_1, J_2
  E = elem 
  E%iloc = i+1
  E%kloc = k+1
  E%val  =   d1*(       HH(ix,i+1,j,k,is)   - HH(ix,i,j,k+1,is)         )   &
   &       + d2*(       HH(ix,i+1,j,k-1,is) - HH(ix,i-1,j,k+1,is)           &
   &                  + HH(ix,i+2,j,k,is)   - HH(ix,i,j,k+2,is)         )
  E%val  = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

  ! g(i+2,k); J_2
  E = elem 
  E%iloc = i+2
  E%kloc = k
  E%val  =   d2*(       HH(ix,i+1,j,k-1,is) - HH(ix,i+1,j,k+1,is)   )
  E%val  = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
  call connect_parallel(E,ingrid)
  if (ingrid) call add_element(E,ierr)

end subroutine jhg_interior
  

!-----------------------------------------------------------------------------
!> adds vdrif grad delta f into the matrix
!! (Term II in the manual)
!!
!! \f$ -(1/Z) i [ T_R E_D {\cal D}^\alpha + T_R v_{\parallel N}^2 \beta^\prime x
!!    {\cal E}^{\psi \alpha} + 2 m_R v_R v_{\parallel N} \Omega {\cal H}^\alpha
!!     + m_R {\cal I}^\alpha \Omega^2 ] k_\alpha  \f$
!!
!!  NB: vcor = \f$ \Omega  \f$
!!
!! The parallel term proportional to the parallel derivative is neglected. 
!< 
!-----------------------------------------------------------------------------

subroutine vdgradf

  use structures,     only : matrix_element
  use control,        only : spectral_radius, order_of_the_radial_scheme
  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod, gx, n_x_grid
  use mode,           only : krho, kxrh
  use matdat,         only : add_element, set_indx
  use matdat,         only : connect_rad
  use constants,      only : ci1
  use dist,           only : ifdis, stencil_side, stencil_side_zf
  use global,         only : id_x
  use geom,           only : dxgr, efun 
  use rotation,       only : shear_real, shear_rate, coriolis, cf_drift
  use neoequil,       only : gradneophi, neo_nsp
  use global,         only : gkw_a_equal_b_accuracy

  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is

  ! The matrix element 
  type (matrix_element) :: elem

  ! error variable 
  integer :: ierr 

  ! dummy 
  integer :: id, ist, ipw, m 
  real    :: drift_x, drift_y, drift_z,ekapka
  real, allocatable :: w(:)
  logical :: ingrid

  allocate(w(1 + 4*max(stencil_side(id_x), stencil_side_zf(id_x))))

  ! Setting which term 
  elem%term  = 'II: vdgradf'
  
  ! type of the elements 
  elem%itype = ifdis
  elem%itloc = ifdis

  if (spectral_radius) then
    elem%ideriv = 0
  else
    elem%ideriv = 1  
  end if
  
  do imod = 1, nmod
    
    if (all(apply_on_imod == 0)) then
      if(.not. lvdgradf) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lvdgradf) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
    
    do is = 1, nsp  ; do ix = 1, nx ; do i = 1, ns ; do j = 1, nmu ; do k = 1, nvpar

        ! store the indices 
        call set_indx(elem,imod,ix,i,j,k,is) 

        ! calculate the drift 
        call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)

        ! if spectral 
        if (spectral_radius) then 

          ! Matrix element : Multiply with -ci1 and wave vectore 
          elem%val = -ci1*(drift_x*kxrh(ix) + drift_y*krho(imod)) 

          !Here is a neoclassical correction term to the Nonlinear term (3)
          !however it is linear (phi is fixed) and therefore the best place 
          !to put it seems to be here.
          if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_rad ) .or. & 
          & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. linear_term_switch_default ) ) then
            if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_equil ) .or. & 
            & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
              if(is <= neo_nsp) then
                ekapka = efun(ix,i,1,1)*kxrh(ix) + efun(ix,i,1,2)*krho(imod)
                elem%val = elem%val - ci1*ekapka*gradneophi(i)

                ekapka = efun(ix,i,3,1)*kxrh(ix) + efun(ix,i,3,2)*krho(imod)
                elem%val = elem%val - ci1*ekapka*dphineods(imod,ix,i)
              endif
            endif
          endif

          
          ! store the element
          call add_element(elem,ierr) 

        else 

          ! background ExB drift (in spectral case is in nonlinear terms)
          ! FIX ME FOR EVEN NX ?
          if (shear_real) then
            drift_y = drift_y + (gx(ix) - real(n_x_grid)/2.)*shear_rate*dxgr
          end if

          ! the differential scheme 
          id  = +1
          ipw = -1 
          if (drift_x.lt.0) ipw = +1  
          ist = 0
          call differential_scheme(ist,ipw,id,w,order_of_the_radial_scheme) 

          elem%val = - ci1 * drift_y * krho(imod) 
          call add_element(elem,ierr) 

          do m = 1, size(w)
            if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
              elem%ixloc = ix + m - ((size(w)+1)/2)
              elem%val   =  - w(m) * drift_x / dxgr
              call connect_rad(elem,ingrid)
              if (ingrid) call add_element(elem,ierr)
            endif
          end do 
 
        endif

      end do ; end do ; end do ; end do ; end do

  end do
  
  deallocate(w)

end subroutine vdgradf

!-----------------------------------------------------------------------------
!> This routine adds perpendicular (hyper)dissipation, coefficients in control
!> Dissipation only on f, not on fields.
!> This dissipation can prevent spectral pile up in nonlinear runs.
!> Some argue dissipation should be isotropic for a good saturation.
!> To compare with a finite difference code may need to put dissipation.
!> In term II, proportional to V_D
!-----------------------------------------------------------------------------
subroutine hyper_disp_perp(disp_x,disp_y)

  use structures,     only : matrix_element
  use control,        only : spectral_radius, order_of_the_radial_scheme
  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod
  use mode,           only : krho, kxrh, kxmax, kymax
  use matdat,         only : add_element, set_indx
  use matdat,         only : connect_rad 
  use general,        only : gkw_abort
  use global,         only : r_tiny, id_x
  use dist,           only : ifdis, stencil_side, stencil_side_zf
  use geom,           only : metric, dxgr
  use constants,      only : ci1
  use mpiinterface,   only : root_processor
  use global,         only : gkw_a_equal_b_accuracy

  real, intent(in) :: disp_x, disp_y

  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is, kpowx, kpowy

  ! dummies 
  integer :: ipw, id, ist, m
  real, allocatable :: w(:)
  logical :: ingrid

  ! the matrix element 
  type (matrix_element) :: elem 

  ! Dummy variables
  real    :: dspx, dspy
  integer :: ierr

  if(.not. spectral_radius) then
    allocate(w(1 + 4*max(stencil_side(id_x), stencil_side_zf(id_x))))
  end if

  dspx=abs(disp_x)
  dspy=abs(disp_y)

  ! type of dissipation depends on the sign of the input coefficient
  if (disp_x < 0.E0) then
    kpowx = 2   ! dissipation
  else
    kpowx = 4   ! hyper-dissipation
  end if
  
  if (disp_y < 0.E0) then
    kpowy = 2   ! dissipation
  else
    kpowy = 4   ! hyper-dissipation
  end if    

  if (spectral_radius) then 
    if (kxmax < r_tiny .or. kymax < r_tiny) then
      call gkw_abort('Invalid call to hyper_disp_perp')
    end if
  else 
    if ((kymax < r_tiny).and.(dspy > r_tiny)) then 
      call gkw_abort('Invalid call to hyper_disp_perp')
    endif 
  endif 

  ! Identifier for the term 
  elem%term = 'hyper_disp_perp'

  ! type of the elements 
  elem%itype = ifdis 
  elem%itloc = ifdis
  
  if (spectral_radius) then
    elem%ideriv = 0
  else
    if (disp_x < 0) then
     elem%ideriv = 2
    else
     elem%ideriv = 4
    end if
  end if

  ! select the scheme. 
  ipw = 1; id = +2;
  ! for the radial direction there are boundary cells, thus one can
  ! use central differences everywhere. No need to use one-sided
  ! stencils near the boundaries.
  ist = 0
  if(.not. spectral_radius) then
    call differential_scheme(ist,ipw,id,w,order_of_the_radial_scheme)
  end if

  if ((root_processor).and.(.not.spectral_radius).and.(disp_x<0.)) &
    & write(*,*)'WARNING: disp_y is set equal to disp_x'

  do is = 1, nsp; do imod = 1, nmod ; do ix = 1, nx ; do i = 1, ns

    do j = 1, nmu ; do k = 1, nvpar

      ! set the indices of the matrix element 
      call set_indx(elem,imod,ix,i,j,k,is)

      if (spectral_radius) then ! hyper-dissipation
      
        elem%val = -(dspy*(krho(imod)/kymax)**kpowy + dspx*(kxrh(ix)/kxmax)**kpowx)
 
        ! The line below is for normal (not hyper) dissipiation 
        ! mat_elem = -(disp_y*krho(imod)**2 + disp_x*kxrh(ix)**2)
        call add_element(elem,ierr) 
   
      else ! dissipation

        ! dissipation in the binormal direction 
        if (dspy.gt.r_tiny) then 
          elem%val = -dspy*(krho(imod)/kymax)**kpowy 
          call add_element(elem,ierr) 
        endif

        ! dissipation in the radial direction 
        if(disp_x.gt.0) then
          ! Fourth order dissipation
          do m = 1, size(w)
            if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
              call set_indx(elem,imod,ix,i,j,k,is)
              elem%ixloc = ix + m - ((size(w)+1)/2)
              elem%val  =  w(m) * dspx / dxgr
              call connect_rad(elem,ingrid)
              if (ingrid) call add_element(elem,ierr)
            endif 
          end do 

        else
          ! Second order dissipation 

          ! Note: (1) the radial variation of the metric elements is neglected
          !       (2) the jacobian is set equal to one

          call set_indx(elem,imod,ix,i,j,k,is)
          elem%ixloc = ix-1
          elem%val = dspx*dxgr*( metric(ix,i,1,1)/(dxgr**2)         &
                   & - metric(ix,i,1,2)*ci1*krho(imod)/dxgr )    
          call connect_rad(elem,ingrid)
          if (ingrid) call add_element(elem,ierr)

          call set_indx(elem,imod,ix,i,j,k,is)
          elem%ixloc = ix
          elem%val = dspx*dxgr*( -2.E0*metric(ix,i,1,1)/(dxgr**2)   &
                   & - metric(ix,i,2,2)*krho(imod)**2 ) 
          call connect_rad(elem,ingrid)
          if (ingrid) call add_element(elem,ierr)
 
          call set_indx(elem,imod,ix,i,j,k,is)
          elem%ixloc = ix+1
          elem%val = dspx*dxgr*( metric(ix,i,1,1)/(dxgr**2)         &
                   & + metric(ix,i,1,2)*ci1*krho(imod)/dxgr )     
          call connect_rad(elem,ingrid)
          if (ingrid) call add_element(elem,ierr)

        endif  

      endif 

    end do ; end do ; end do ; end do ; end do ; end do

    if(.not. spectral_radius) deallocate(w)

end subroutine hyper_disp_perp

!------------------------------------------------------------------------------
!> This routine adds the damping of the solution in the boundary regions into 
!> the matrix. 
!------------------------------------------------------------------------------
subroutine krook_bound

  use structures, only : matrix_element
  use matdat,     only : add_element, set_indx
  use dist,       only : ifdis 
  use krook,      only : nlbound, gammab, bwidth 
  use grid,       only : gx, n_x_grid, nmod, nx, ns, nmu, nvpar, nsp 
  use mode,       only : iyzero

  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is

  ! the matrix element 
  type (matrix_element) :: elem 

  ! dummy 
  integer :: ierr 
 
  if (.not.nlbound) return 

  ! Identifier for the term (not in timestep estimator since no derivative)
  elem%term = 'Krook operator for boundary'

  ! type of the elements 
  elem%itype = ifdis 
  elem%itloc = ifdis
  elem%ideriv = 0

  do ix = 1, nx 

    do imod = 1, nmod; do i = 1, ns; 
      do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp 

        ! set the indices of the matrix element (diagonal, no derivative)
        call set_indx(elem,imod,ix,i,j,k,is)

        ! default zero damping 
        elem%val = 0. 

        if (bwidth > 0) then

          ! inner boundary 
          if (gx(ix).le.bwidth) then  
            elem%val = -gammab*(exp(-(gx(ix)-1)/bwidth) - exp(-1.)) / & 
                     & (1 - exp(-1.)) 
          endif 

          ! outer boundary 
          if (gx(ix).gt.n_x_grid-bwidth) then 
            elem%val = -gammab*(exp(-(n_x_grid-gx(ix))/bwidth) - exp(-1.))/ &
                     &  (1. - exp(-1.)) 
          endif
        
        else  ! heaviside damping layer on 0 mode only               
           
          if ((gx(ix).le.abs(bwidth) .or. gx(ix).gt.n_x_grid-abs(bwidth)).and. imod == iyzero) then  
            elem%val = -gammab 
          endif   
          
        end if

        ! store the element 
        call add_element(elem,ierr) 

      end do; end do; end do  
    end do; end do

  end do 

end subroutine krook_bound

!------------------------------------------------------------------------------
!> This routine adds the \f$  {\bf v}_\chi \nabla F_M \f$ term
!! Term V in the manual
!! This routine would be better named vchi_grad_fm as
!! it now also includes v_del_B_perp
!!
!! \f$  +  {\cal I} E^(\alpha \psi) k_\alpha chi 
!! [ 1/L_N + E_T 1/L_T + 2 v_{\parallel N} R_N B_t/B u^\prime / v_R   \f$
!!
!! in the matrix.
!!
!! Note that because the diagonal components of efun are zero, the components
!! containing kxrh are actually zero, hence for the nonspectral method no 
!! radial derivatives need to be calculated.
!< 
!------------------------------------------------------------------------------
subroutine ve_grad_fm

  use structures,     only : matrix_element
  use control,        only : nlapar, nlbpar, spectral_radius
  use control,        only : order_of_the_radial_scheme
  use grid,           only : nx,ns,nmu,nvpar,nsp,nmod
  use velocitygrid,   only : vpgr, mugr
  use dist,           only : fmaxwl,falpha, iphi, iphi_ga, iapar, iapar_ga 
  use dist,           only : ibpar, ibpar_ga, ifdis, f_EP, df_EPdv
  use dist,           only : stencil_side, stencil_side_zf
  use components,     only : fp, tp, vthrat, types, pbg, tmp, signz, types
  use components,     only : tgrid, vpar_mean
  use mode,           only : krho, kxrh
  use matdat,         only : add_element, set_indx
  use geom,           only : efun, bn, dfun
  use constants,      only : ci1, pi
  use neoequil,       only : neof, gradneof, neo_nsp
  use index_function, only : indx
  use matdat,         only : connect_rad
  use mode,           only : erase_any_drive, no_drive_of_modes
  use global,         only : gkw_a_equal_b_accuracy, id_x

  ! indices integers 
  integer :: ix, i, j, k, is, imod
  integer :: id, ipw, m, ist
  ! The matrix element 
  type (matrix_element) :: elem, elem2, elem3  

  complex :: mat_elem_electrostatic
  real    :: dum, b, b1_mod, ekapka, vn, ET, dneofds, dum2
  real, allocatable :: w(:)
  integer :: ierr  
 
  ! Switch for correction due to bpar in term V
  logical :: l_term5_bpar = .true.
  logical :: ingrid
  logical :: erase_drive_this_imod 

  ! vE_grad_fm
  elem%ideriv = 0
  elem%term  = 'V: ve_grad_fm' 
  elem%itype = ifdis 
  if (spectral_radius) then 
    elem%itloc = iphi 
  else 
    elem%itloc = iphi_ga 
  endif 

  ! Electromagnetic correction, v_del_B_perp
  elem2%ideriv = 0
  elem2%term  = 'V: v_del_B_perp_grad_fm' 
  elem2%itype = ifdis 
  if (spectral_radius) then 
    elem2%itloc = iapar 
  else
    elem2%itloc = iapar_ga
  endif 

  ! Magnetic field compression correction del_v_gradB 
  elem3%ideriv = 0
  elem3%term  = 'V: V_del_gradB_grad_fm'   
  elem3%itype = ifdis 
  if (spectral_radius) then 
    elem3%itloc = ibpar 
  else 
    elem3%itloc = ibpar_ga 
  endif 

  ! provide a buffer large enough to hold the finite difference coefficients of a
  ! completely one-sided stencil, for the chosen order of the scheme
  allocate(w(1 + 4*max(stencil_side(id_x), stencil_side_zf(id_x))))

  ! calculate the terms due to the Maxwell background
  binormal_loop: do imod=1,nmod;
    
    if (all(apply_on_imod == 0)) then
      if(.not. lve_grad_fm) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lve_grad_fm) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
  
    erase_drive_this_imod = .false.
    if(erase_any_drive) then
      ! find out if this imod is in the list of modes which shall not be
      ! driven
      do i = 1, size(no_drive_of_modes)
        ! if it is found, then do not add matrix elements for this
        ! term
        if(imod == no_drive_of_modes(i)) erase_drive_this_imod = .true.
      end do
    end if
    
    do is=1,nsp; do ix=1,nx; do i=1,ns; do j=1,nmu; do k=1,nvpar

    ! E^(alpha psi) k_alpha
    ! (first component is always zero because diagonal of efun is zero)
    ekapka = efun(ix,i,1,1)*kxrh(ix) + efun(ix,i,2,1)*krho(imod)
    
    if (types(is) .eq. 'EP') then
      dum = dmaxwel(ix,i,j,k,is,imod)*f_EP(ix,i,j,k,is)*ekapka +                  &
         & 0.5E0* (                                                          &
         &  (vpgr(i,j,k,is)-vpar_mean)**2/(tmp(ix,is)/tgrid(is))*            &
         &     exp(-(vpgr(i,j,k,is)-vpar_mean)**2/(tmp(ix,is)/tgrid(is)))+   &
         &  (vpgr(i,j,k,is)+vpar_mean)**2/(tmp(ix,is)/tgrid(is))*            &
         &     exp(-(vpgr(i,j,k,is)+vpar_mean)**2/(tmp(ix,is)/tgrid(is)))    &
         &  )                                                                &
         & * exp(-2.E0*bn(ix,i)*mugr(j)/(tmp(ix,is)/tgrid(is)))              &
         & *tp(ix,is)*0.5E0/(sqrt(tmp(ix,is)*pi/tgrid(is))**3)*ekapka
    else
      dum = dmaxwel(ix,i,j,k,is,imod)*fmaxwl(ix,i,j,k,is)*ekapka
    endif 

    if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_equil ) .or. & 
    & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
      if(is <= neo_nsp)then
        !Adds the radial gradient of the neoclassical correction 
        !to the maxwellian background distribution function
        !The radial derivative of Fm x h is done in two parts.
        !This adds hdFm/dr
        if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_rad ) .or. & 
        & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. linear_term_switch_default ) ) then
          dum = dum*(1.E0 + neof(3,is,i,j,k))
          !this is Fmdh/dr
          dum = dum -  efun(ix,i,2,1)*krho(imod)*fmaxwl(ix,i,j,k,is)*gradneof(is,i,j,k)
          !This term is due to the radial derivative of the magnetic field that appears when you take a
          !derivative of a Maxwellian
          ekapka = dfun(ix,i,1)*kxrh(ix) + dfun(ix,i,2)*krho(imod)
          dum = dum - ekapka*fmaxwl(ix,i,j,k,is)*neof(3,is,i,j,k)*mugr(j)*bn(ix,i)/tmp(ix,is)
        endif

        !The i is multiplied below
        
        !The neoclassical correction related to the parallel derivative of the
        !background distribution function
        ekapka = efun(ix,i,1,3)*kxrh(ix) + efun(ix,i,2,3)*krho(imod)
        !Calculate the parallel derivative
        dneofds = neo_dFds(ix,i,j,k,is)
        dum = dum - ekapka*dneofds
      endif
    endif

    ! for an alpha particle distribution change to
    if (types(is) .eq. 'alpha') then
      vn = sqrt(vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i))
      ET =   3.E0 / 2.E0 * (1.E0 / (log(1.E0 + 27.E0*pbg(ix,is)**(-1.5))*        &
          & (1.E0+pbg(ix,is)**(1.5)/27E0)) - pbg(ix,is)**1.5/ (pbg(ix,is)**1.5 + &
          & vn**3))
      ekapka = efun(ix,i,1,1)*kxrh(ix) + efun(ix,i,2,1)*krho(imod) 
      dum = (fp(ix,is) + ET*tp(ix,is))*falpha(i,j,k)*ekapka
    endif

   
    ! For energetic particle distribution, a term is added as a sum of vchi_gradFeq
    ! and vchi_gradB_dFeq_dvpar terms
    ! This part adds the electrostatic contribution to the common factor of matrix elements
    if (types(is) .eq. 'EP') then 
      dum2 = - mugr(j) * bn(ix,i) / (tmp(ix,is)/tgrid(is)) * (f_EP(ix,i,j,k,is) &
        & + 0.5E0*(tmp(ix,is)/tgrid(is)/vpgr(i,j,k,is))*df_EPdv(ix,i,j,k,is))
      dum = dum + dfun(ix,i,2)*krho(imod)*dum2

      if (spectral_radius) then
        dum = dum + dfun(ix,i,1)*kxrh(ix)*dum2 ! &
        !  &  mugr(j) * bn(ix,i) / (tmp(ix,is)/tgrid(is)) * (f_EP(ix,i,j,k,is) &
        !  & + 0.5E0*(tmp(ix,is)/tgrid(is)/vpgr(i,j,k,is))*df_EPdv(ix,i,j,k,is))
      else
        !the differencial scheme
        id  = +1
        ipw = -1
        if (dum2 .lt. 0) ipw = +1
        ist=0
        call differential_scheme(ist,ipw,id,w,order_of_the_radial_scheme)

        call set_indx(elem,imod,ix,i,j,k,is)
          elem%val = dfun(ix,i,1) * dum2 ! mugr(j) * bn(ix,i) / &
            !& (tmp(ix,is)/tgrid(is)) * (f_EP(ix,i,j,k,is) &
            !& + 0.5E0*(tmp(ix,is)/tgrid(is)/vpgr(i,j,k,is))*df_EPdv(ix,i,j,k,is))
        call add_element(elem,ierr)
       
        do m = 1, size(w)
          if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
            call set_indx(elem,imod,ix,i,j,k,is)
            elem%ixloc = ix + m - ((size(w)+1)/2)
            elem%val   =  w(m) * elem%val 
            call connect_rad(elem,ingrid)
            if (ingrid) call add_element(elem,ierr)
          endif
        end do
  
        if (nlapar) then  
          call set_indx(elem2,imod,ix,i,j,k,is)
          elem2%val = - 2.0E0*vthrat(is)*vpgr(i,j,k,is) * dfun(ix,i,1) * dum2!dfun(ix,i,1) &
            !& * mugr(j) * bn(ix,i) / (tmp(ix,is)/tgrid(is)) * ( f_EP(ix,i,j,k,is) &
            !& + 0.5E0*(tmp(ix,is)/tgrid(is)/vpgr(i,j,k,is))*df_EPdv(ix,i,j,k,is) )
          call add_element(elem2,ierr)  

          do m = 1, size(w)
            if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
              call set_indx(elem2,imod,ix,i,j,k,is)
              elem%ixloc = ix + m - ((size(w)+1)/2)
              elem%val   =  w(m) * elem2%val
              call connect_rad(elem2,ingrid)
              if (ingrid) call add_element(elem2,ierr)
            endif
          end do
        end if ! nlapar
               
      endif ! spectral_radius 
    
    endif ! types(is)='EP'

    ! common factor of the matrix elements
    mat_elem_electrostatic = ci1*dum

    ! the bessel functions 
    b = wrap_beslj0_gkw(imod,ix,i,j,is)
    b1_mod = wrap_mod_besj1_gkw(imod,ix,i,j,is)

    ! vE_grad_fm
    ! WARNING THIS IF STATEMENT ERASES ALSO nlapar but not nlbpar
    if(.not. erase_drive_this_imod) then
      call set_indx(elem,imod,ix,i,j,k,is) 
      elem%val   = mat_elem_electrostatic*b
      call add_element(elem,ierr)
    end if

    ! Electromagnetic correction, v_del_B_perp
    if (nlapar) then     
      call set_indx(elem2,imod,ix,i,j,k,is)
      elem2%val = - 2.*elem%val*vthrat(is)*vpgr(i,j,k,is)
      call add_element(elem2,ierr)
    endif

    ! Magnetic field compression correction del_v_gradB 
    if (nlbpar .and. l_term5_bpar) then
      call set_indx(elem3,imod,ix,i,j,k,is)
      elem3%val = mat_elem_electrostatic*b1_mod*2.0*mugr(j)*tmp(ix,is)/signz(is)
      call add_element(elem3,ierr)
    end if  

   
      end do;end do; end do ; end do ; end do ;
    end do binormal_loop

 end subroutine ve_grad_fm

!------------------------------------------------------------------------------
!> This routine puts the drift in the gradient of phi times the maxwel
!! Term VIII in the manual
!! \f$
!! - [E_D {\cal D}^\alpha + \beta^\prime v_{\parallel N}^2 {\cal E}^{\psi \alpha}
!! + 2 (m_R v_R / T_R) v_{\parallel N} {\cal H}^\alpha \Omega + m_R/T_R 
!!  {\cal I}^\alpha \Omega^2 ] k_\alpha \chi_N F_{MN}  \f$
!<-----------------------------------------------------------------------------
subroutine vd_grad_phi_fm

  use structures, only : matrix_element
  use control,        only : nlbpar, spectral_radius, order_of_the_radial_scheme 
  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod
  use velocitygrid,   only : mugr, vpgr 
  use geom,           only : bn 
  use dist,           only : fmaxwl, falpha
  use dist,           only : df_EPdv
  use dist,           only : iphi, iphi_ga 
  use dist,           only : ibpar, ibpar_ga, ifdis, stencil_side, stencil_side_zf
  use global,         only : id_x
  use mode,           only : krho, kxrh
  use components,     only : types, signz, pbg, tmp
  use matdat,         only : add_element, set_indx
  use matdat,         only : connect_rad 
  use constants,      only : ci1
  use rotation,       only : coriolis, cf_drift
  use geom,           only : dxgr 
  use neoequil,       only : neo_nsp
  use global,         only : gkw_a_equal_b_accuracy

  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is

  ! the matrix element 
  type (matrix_element) :: elem, elem2 

  ! reference values and matrix element
  complex :: mat_elem_electrostatic = (0.0, 0.0)

  ! Dummy variables
  integer :: ipw, id, ist, m 
  real    :: b, b1_mod, dum, vn, drift_x, drift_y, drift_z
  real, allocatable :: w(:) 
  real    :: dneofdv
  integer :: ierr 

  !Switch for the correction due to bpar in term VIII
  logical :: l_term8_bpar = .true.
  logical :: ingrid 

  allocate(w(1 + 4*max(stencil_side(id_x), stencil_side_zf(id_x))))
  
  ! identifier of the terms 
  elem%term  = 'VIII: vd_grad_phi_fm'
  elem%itype = ifdis   
  if (spectral_radius) then 
    elem%ideriv = 0
    elem%itloc = iphi 
  else 
    elem%ideriv = 1
    elem%itloc = iphi_ga
  endif

  !Correction due to bpar
  elem2%term  = 'XI: vd_grad_bpar_fm'
  elem2%itype = ifdis 
  if (spectral_radius) then 
    elem2%ideriv = 0
    elem2%itloc = ibpar
  else 
    elem2%ideriv = 1
    elem2%itloc = ibpar_ga
  endif   

  do imod = 1, nmod
  
    if (all(apply_on_imod == 0)) then
      if(.not. lvd_grad_phi_fm) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lvd_grad_phi_fm) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
  
    do is = 1, nsp ; do ix = 1, nx ; do i = 1, ns; do j = 1, nmu ; do k = 1, nvpar

      ! The bessel function for gyro-averaging
      b = wrap_beslj0_gkw(imod,ix,i,j,is)
      b1_mod = wrap_mod_besj1_gkw(imod,ix,i,j,is)

      ! calculate the drift 
      call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.) 

      ! the matrix elements 
      call set_indx(elem,imod,ix,i,j,k,is) 
      call set_indx(elem2,imod,ix,i,j,k,is) 

      if (spectral_radius) then 

        dum = signz(is)*(drift_x*kxrh(ix) + drift_y*krho(imod))

        mat_elem_electrostatic = - ci1*dum*fmaxwl(ix,i,j,k,is)/tmp(ix,is)

        if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_equil ) .or. & 
        & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
          if(is <= neo_nsp)then
            !The gradient drift is removed here as it doesnt cancel with another term elsewhere.
            call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.false.) 
            dum = signz(is)*(drift_x*kxrh(ix) + drift_y*krho(imod))

            dneofdv = neo_dFdvpar(ix,i,j,k,is)
            mat_elem_electrostatic = mat_elem_electrostatic  &
              & + 0.5E0*ci1*dum*dneofdv/(tmp(ix,is)*vpgr(i,j,k,is))
          endif
        endif
        !if (energetic_particles) then
        !  mat_elem_electrostatic = mat_elem_electrostatic &
        !                       & - ci1*dum*f_EP(ix,i,j,k,is)/T_EP &
        !                       & + ci1*dum*df_EPdv_sinhc(ix,i,j,k,is)/T_EP
        !end if
        elem%val = b*mat_elem_electrostatic
        if (types(is) .eq. 'EP') then
          mat_elem_electrostatic = 0.5E0*ci1*dum* &
            & df_EPdv(ix,i,j,k,is)/vpgr(i,j,k,is)/tmp(ix,is)
          elem%val = b*mat_elem_electrostatic
        endif  
       ! if (types(is) .eq. 'EP') then
       !  elem%val = b*(-ci1*dum*f_EP(ix,i,j,k,is) &
       !                      & + ci1*dum*df_EPdv_sinhc(ix,i,j,k,is))
       !
       ! endif
        if (types(is) .eq. 'alpha') then
          vn = sqrt(vpgr(i,j,k,is)**2+2.E0*mugr(j)*bn(ix,i))
          elem%val = -ci1*dum*b*falpha(i,j,k)*(3./2.)*vn / (pbg(ix,is)**1.5 + vn**3)
        endif

        call add_element(elem,ierr)

      else

        ! the differential scheme 
        id  = +1
        ipw = -1 
        if (drift_x.lt.0) ipw = +1  
        ist = 0 
        call differential_scheme(ist,ipw,id,w,order_of_the_radial_scheme) 

        call set_indx(elem,imod,ix,i,j,k,is)
        elem%val = - ci1 * signz(is) *  drift_y * krho(imod) * fmaxwl(ix,i,j,k,is) & 
                 &   / tmp(ix,is)! / tgrid(is)

      !  if (types(is) .eq. 'EP') then
      !     elem%val = &
      !           & - ci1 * signz(is) *  drift_y * krho(imod) * f_EP(ix,i,j,k,is) &
      !           &   / tmp(ix,is) &
      !           & + ci1 * signz(is) *  drift_y * krho(imod) * df_EPdv_sinhc(ix,i,j,k,is) &
      !           &   / tmp(ix,is)
      !  end if

       if (types(is) .eq. 'EP') then 
         elem%val = 0.5E0 * ci1 * signz(is) * drift_y * krho(imod) * df_EPdv(ix,i,j,k,is) &
               & / vpgr(i,j,k,is)
       end if

        call add_element(elem,ierr)

        do m = 1, size(w)
          if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
            call set_indx(elem,imod,ix,i,j,k,is)      
            elem%ixloc = ix + m - ((size(w)+1)/2)
            elem%val   =  - w(m) * drift_x / dxgr 
            elem%val   =  elem%val * signz(is) * &
                         &  fmaxwl(ix,i,j,k,is) / tmp(ix,is)!/ tgrid(is)
            if (types(is) .eq. 'EP') then
              elem%val = &
                       &  w(m) * drift_x / dxgr               &
                       &  * signz(is) * df_EPdv(ix,i,j,k,is)    &
                       &  * 0.5E0 / vpgr(i,j,k,is)
            end if
            call connect_rad(elem,ingrid)
            if (ingrid) call add_element(elem,ierr)
          endif
        end do 

      endif 

      if (spectral_radius) then 

        !Correction due to bpar
        if (nlbpar .and. l_term8_bpar) then
          elem2%val = mat_elem_electrostatic*b1_mod*2.0*tmp(ix,is)*mugr(j)/signz(is)
          call add_element(elem2,ierr)
        endif

      else 

        if (nlbpar.and.l_term8_bpar) then
          write(*,*) 'This hasnt been programmed correctly yet'
          stop 1
    
          ! the differential scheme 
          id  = +1
          ipw = -1 
          if (drift_x.lt.0) ipw = +1  
          ist = 0 
          call differential_scheme(ist,ipw,id,w,order_of_the_radial_scheme) 

          do m = 1, size(w)
            if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
              call set_indx(elem2,imod,ix,i,j,k,is)
              elem2%ixloc = ix + m - ((size(w)+1)/2)
              elem2%val   =  - w(m) * drift_x / dxgr 
              ! if (m.eq.3) elem2%val = elem2%val - ci1 * drift_y*krho(imod) 
              elem2%val   =  2.0 * elem2%val * fmaxwl(ix,i,j,k,is) * mugr(j) 
              call connect_rad(elem2,ingrid)
              if (ingrid) call add_element(elem2,ierr)
            endif
          end do 

        endif 

      endif 

    end do ; end do ; end do ; end do ; end do

  end do

  deallocate(w)
end subroutine vd_grad_phi_fm

!------------------------------------------------------------------------------
!> This routine puts the Landau damping term (Term VII in the manual)
!! \f$ 
!!
!! - (Z/T_R) v_R v_{\parallel N} {\cal F} (d <\phi> / d s) F_{MN}      \f$ 
!!
!! in the matrix. As well as the trapping term due to the perturbed parallel 
!! magnetic field.  
!!
!! The differential scheme is set through the use of the differential_scheme
!! routine.  
!! The paralllel boundary conditions are implemented through calls to
!! connect_parallel (which are made from add_element).
!<------------------------------------------------------------------------------
subroutine vpar_grd_phi

  use structures,     only : matrix_element
  use control,        only : nlbpar, nlapar, spectral_radius 
  use dist,           only : fmaxwl, f_EP, df_EPdv, df_EPdv_sinhc
  use dist,           only : falpha, iphi, iphi_ga, f_EP, df_EPdv
  use dist,           only : iapar, iapar_ga, ibpar, ibpar_ga, ifdis
  use geom,           only : ffun, bn, sgr_dist,efun, dfun
  use matdat,         only : set_indx, add_element
  use matdat,         only : connect_parallel, pos_par_grid
  use grid,           only : nx,ns,nmu,nvpar,nsp,nmod, gx
  use velocitygrid,   only : vpgr, mugr
  use components,     only : tmp, signz, vthrat, types, pbg, tgrid
  use components,     only : rhostar_linear, mas, rhostar
  use dist,           only : ifdis, stencil_side, stencil_side_zf
  use global,         only : id_s
  use mode,           only : parallel_phase_shift
  use rotation,       only : coriolis, cf_drift
  use rho_par_switch, only : lvdgrad_phi_fm_rhostar, lve_grad_fm_rhostar
  use neoequil,       only : neo_nsp
  use global,         only : gkw_a_equal_b_accuracy

  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is

  ! variables for the parallel boundary conditions
  integer :: ist, ierr
  logical :: ingrid

  ! element to test/add to the matrix
  type (matrix_element) :: elem, elem2, elem3

  ! Dummy variables
  real    :: b, b1_mod, dum, vn, drift_x, drift_y, drift_z
  real    :: term7, term8, term5, term9, type_alpha = 0.0, efem, dneofdv
  real, allocatable :: w(:)
  integer :: ipw, id, m 

  if (ns .lt. 2) return

  allocate(w(1 + 4*max(stencil_side(id_s), stencil_side_zf(id_s))))

  elem%term = 'VII: vpar_grd_phi (Landau damping)'
  if (spectral_radius) then 
    elem%itloc = iphi
  else 
    elem%itloc = iphi_ga
  endif
  elem%itype = ifdis
  elem%ideriv = 1

  elem2%term = 'X: Trapping due to the perturbed magnetic field'
  if (spectral_radius) then 
    elem2%itloc = ibpar
  else 
    elem2%itloc = ibpar_ga
  endif 
  elem2%itype = ifdis   
  elem2%ideriv = 1

  elem3%term = 'V: v_del_B_perp_grad_fm (rhostar)'
  if (spectral_radius) then 
    elem3%itloc = iapar
  else 
    elem3%itloc = iapar_ga
  endif 
  elem3%itype = ifdis 
  elem3%ideriv = 1

  do imod = 1,nmod
    
    if (all(apply_on_imod == 0)) then
      if(.not. lvpgrphi) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lvpgrphi) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif

    do is = 1,nsp ; do ix = 1,nx ; do j = 1,nmu ; do k = 1,nvpar ; do i = 1,ns

      ! parallel derivative part of term VII vpar grad phi FM
      efem = fmaxwl(ix,i,j,k,is)

      term7 = -signz(is)*ffun(ix,i)*vthrat(is)*vpgr(i,j,k,is)*efem/tmp(ix,is)

      if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_equil ) .or. & 
      & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
        if(is <= neo_nsp)then
          dneofdv = neo_dFdvpar(ix,i,j,k,is)
          term7 = term7 + 0.5E0*signz(is)*ffun(ix,i)*vthrat(is)*dneofdv/tmp(ix,is)
        endif
      endif

      if (types(is) .eq. 'EP') then
        term7 = 0.5E0*signz(is)*ffun(ix,i)*vthrat(is)*df_EPdv(ix,i,j,k,is)/tmp(ix,is)
      end if

      ! parallel derivative part of term VIII vd grad phi FM
      call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
      if ( lvdgrad_phi_fm_rhostar ) then
        term8 = - drift_z * rhostar_linear * signz(is) * fmaxwl(ix,i,j,k,is) / tmp(ix,is)
        if (types(is) .eq. 'EP') then
          term8 = &
              & - drift_z * rhostar_linear * signz(is) * f_EP(ix,i,j,k,is) / (tmp(ix,is)) &
              & + drift_z * rhostar_linear * signz(is) * df_EPdv_sinhc(ix,i,j,k,is) / (tmp(ix,is))
        end if
      else
        term8=0
      endif
      ! common factor of parallel derivative part of term V vchi grad FM
      if (lve_grad_fm_rhostar ) then
        term5 = rhostar_linear * dmaxwel(ix,i,j,k,is,imod) * fmaxwl(ix,i,j,k,is) * efun(ix,i,3,1)
      else
        term5=0
      endif
     
      !Parallel derivative part of the term added when the distribution function is not a maxwellian
      if (types(is) .eq. 'EP') then
        term9 = - rhostar*mugr(j)/(tmp(ix,is)/tgrid(is))*bn(ix,i)*( f_EP(ix,i,j,k,is) + &
          & 0.5E0*(tmp(ix,is)/tgrid(is)/vpgr(i,j,k,is))*df_EPdv(ix,i,j,k,is) )*dfun(ix,i,3)
      !write(*,*) term9
      else
        term9=0
      end if

      dum = term5 + term7 + term8 + term9

      ! for an alpha particle distribution replace
      if (types(is) .eq. 'alpha') then
        vn = sqrt(vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i))
        dum = -signz(is)*ffun(ix,i)*vthrat(is)*vpgr(i,j,k,is)*falpha(i,j,k) &
              & *3.E0*vn/(2.E0*tmp(ix,is)*(pbg(ix,is)**1.5+vn**3))
        type_alpha = dum
      endif

      ! dum * dphi/ds 
      ! direction of the parallel motion 
      if (dum < 0) then 
        ipw = 1 
      else 
        ipw = -1 
      end if 
      id = + 1 

      ! select the scheme 
      ist = pos_par_grid(imod,ix,i,k)
      call differential_scheme(ist,ipw,id,w) 

      do m = 1, size(w)
 
        if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
 
          ! grd_phi 
          call set_indx(elem,imod,ix,i,j,k,is)       
 
          elem%iloc  = i + m - ((size(w)+1)/2)
          elem%val   = w(m) * dum / sgr_dist  
          elem%val   = elem%val * parallel_phase_shift(elem%imod, & 
             elem%ix,elem%i,elem%iloc)
          call connect_parallel(elem,ingrid)
          if (ingrid) then 
            b = wrap_beslj0_gkw(elem%imloc,elem%ixloc,elem%iloc,elem%jloc,elem%isloc)
            elem%val   = elem%val * b 
            call add_element(elem,ierr)
          endif
        endif
      end do 

      if ( nlbpar ) then

        dum = 2 * mugr(j) * tmp(ix,is) * term5 / signz(is)
        dum =  dum + term7 + term8

        ! for an alpha particle distribution replace
        if (types(is) .eq. 'alpha') then
          dum = type_alpha
        endif


        ! dum * dB1par/ds 
        ! direction of the parallel motion 
        if (dum < 0) then 
          ipw = 1 
        else 
          ipw = -1 
        end if 
        id = + 1 

        ! select the scheme 
        ist = pos_par_grid(imod,ix,i,k)
        call differential_scheme(ist,ipw,id,w) 

        ! parallel rho* effects - electromagnetic compressional
        ! somehow it would be nice to make it clear that these 
        ! terms are not generally in use - perhaps hide them 
        ! their own wrapper function
        do m = 1, size(w)
          if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then

            ! trapping due to magnetic field compression 
            call set_indx(elem2,imod,ix,i,j,k,is)
                  
            elem2%iloc  = i + m - ((size(w)+1)/2)
            elem2%val   = w(m) * dum / sgr_dist 
            elem2%val   = elem2%val * parallel_phase_shift(elem2%imod, & 
               elem2%ix,elem2%i,elem2%iloc)
            call connect_parallel(elem2,ingrid)

            if (ingrid) then 
              b1_mod = wrap_mod_besj1_gkw(elem2%imloc,elem2%ixloc,elem2%iloc,elem2%jloc,elem2%isloc)
              b          = 2.0*b1_mod*mugr(j)*tmp(ix,is)/signz(is)
              elem2%val   = elem2%val * b
              call add_element(elem2,ierr)
            endif 

          endif 
        end do
      end if 

      ! parallel rho* effects - electromagnetic
      ! somehow it would be nice to make it clear that these 
      ! terms are not generally in use - perhaps hide them 
      ! their own wrapper function
      if ( nlapar .and. rhostar_linear > 0.) then
        !parallel derivation of term V (electrodynamik potential)
        dum= -2 * tmp(ix,is) / vthrat(is) / mas(is) * vpgr(i,j,k,is) * (term5+term9) / signz(is)

        ! for an alpha particle distribution replace
        if (types(is) .eq. 'alpha') then
          dum = type_alpha
        endif

        ! dum * dapar/ds
        ! direction of the parallel motion 
        if (dum < 0) then 
          ipw = 1 
        else 
          ipw = -1 
        end if 
        id = + 1 

        ! select the scheme 
        ist = pos_par_grid(imod,ix,i,k)
        call differential_scheme(ist,ipw,id,w) 

        if ( nlapar ) then   

        ! select the scheme 
        ist = pos_par_grid(imod,ix,i,k)
        call differential_scheme(ist,ipw,id,w)

        do m = 1, size(w)
          if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then

            ! rho_star effect term V 
            call set_indx(elem3,imod,ix,i,j,k,is)

            elem3%iloc  = i + m - ((size(w)+1)/2)
            elem3%val   = w(m) * dum / sgr_dist 
            elem3%val   = elem3%val * parallel_phase_shift(elem3%imod, & 
               elem3%ix,elem3%i,elem3%iloc)
            call connect_parallel(elem3,ingrid)
            if (ingrid) then 
              b = wrap_beslj0_gkw(elem3%imloc,elem3%ixloc,elem3%iloc,elem3%jloc,elem3%isloc)
              elem3%val   = elem3%val * b 
              call add_element(elem3,ierr)
            endif 

          endif 
        end do
      endif
    end if 

    end do ; end do ; end do ; end do ; end do  
  end do 

  deallocate(w)
end subroutine vpar_grd_phi


subroutine put_unity_block(mat_name, val)
  use structures,     only : matrix_element
  use grid,           only : nx,ns,nmu,nvpar,nsp,nmod
  use matdat,         only : add_element, set_indx
  use dist,           only : ifdis
  use control, only : spectral_radius
  character(len=*),intent(in) :: mat_name
  complex, intent(in) :: val
  integer :: imod, ix, i, j, k, is
  type (matrix_element) :: elem 
  integer :: ierr
  ! non spectral case is done elsewhere (in gyro_average)
  if (.not. spectral_radius) return

  elem%term = mat_name
  elem%itype = ifdis
  elem%itloc = ifdis
  elem%ideriv = 0

    do imod = 1, nmod
      do ix = 1, nx
        do i = 1, ns
          do j = 1, nmu
            do k = 1, nvpar
              do is = 1, nsp
                call set_indx(elem,imod,ix,i,j,k,is)

                elem%val = val
                call add_element(elem,ierr)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine put_unity_block

!------------------------------------------------------------------------------
!> This routine adds the integral part of the Poisson equation in the matrix
!> for the spectral case
!------------------------------------------------------------------------------
subroutine poisson_int

  use structures,     only : matrix_element
  use grid,           only : nx,ns,nmu,nvpar,nsp,nmod, number_of_species 
  use mpicomms,       only : COMM_SP_NE
  use components,     only : de, signz, tmp, veta
  use geom,           only : bn
  use matdat,         only : add_element, set_indx
  use velocitygrid,   only : intvp, intmu, mugr
  use control,        only : nlbpar, spectral_radius
  use mode,           only : krloc
  use rotation,       only : cfen
  use dist,           only : ifdis, iphi 
  use functions,      only : gamma1_gkw, gamma_gkw, besselj0_gkw, mod_besselj1_gkw 
  use mpiinterface,   only : mpiallreduce_sum_inplace
  use global,         only : gkw_a_equal_b_accuracy

  real :: gamma_diff, bes_j0, mod_bes_j1

  !coefficients used if bpar is on:
  real :: I_sp1, I_sp2, B_sp1, B_sp2

  !dummy for the mpi allreduce
  real :: dum

  ! indices for the distribution 
  integer :: imod, ix, i, j, k, is

  ! the matrix element 
  type (matrix_element) :: elem 

  ! error variable 
  integer :: ierr

  ! non spectral case is done elsewhere (in gyro_average)
  if (.not. spectral_radius) return

  ! The term 
  elem%term = 'poisson_int' 

  ! type of the elements 
  elem%itype = iphi 
  elem%itloc = ifdis 
  elem%ideriv = 0

  do imod = 1, nmod
    
    if (all(apply_on_imod == 0)) then
      if(.not. lpoisson) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lpoisson) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
  
    if (nlbpar) then 

      do ix = 1, nx ; do i = 1, ns ; do j = 1, nmu ; do k = 1, nvpar

        I_sp1 = 0.E0
        I_sp2 = 0.E0
        B_sp1 = 0.E0
        B_sp2 = 0.E0

        !calculating coefficients B_sp1, B_sp2, F_sp1, F_sp2 locally:
        !in case of zero mode the coefficients are simplified since gamma(0)=1
        !and gamma_diff(0)=1.

        do is = 1, nsp

          if (abs(krloc(imod,ix,i)) < 1.0E-5) then !zero mode
            !Gamma functions together with correction due to high rotation
            gamma_diff = 1.0*exp(-cfen(i,is))
          else !non-zero modes
            !Gamma functions together with correction due to high rotation
            gamma_diff = (gamma_gkw(imod,ix,i,is) - gamma1_gkw(imod,ix,i,is))*exp(-cfen(i,is))
          end if

          B_sp1 = B_sp1 + signz(is)*de(ix,is)*gamma_diff/(bn(ix,i))
          B_sp2 = B_sp2 + tmp(ix,is)*de(ix,is)*veta(ix)*gamma_diff/(bn(ix,i)**2)

        end do

        ! calculating coefficients B_sp1, B_sp2, F_sp1, F_sp2 globally
        if (number_of_species > nsp) then

          call mpiallreduce_sum_inplace(B_sp1,1,COMM_SP_NE)
          call mpiallreduce_sum_inplace(B_sp2,1,COMM_SP_NE)

        end if

        do is = 1, nsp

          !calculating coefficients I_sp1, I_sp2. Their sum over velocity space and species
          !is done in exp_integration.

          call set_indx(elem,imod,ix,i,j,k,is)

          bes_j0 = besselj0_gkw(imod,ix,i,j,is)
          !modified J_1 Bessel fundtion: mod_bes_j1 = 2*J_1(k_perp rho)/(k_perp rho):
          !the limit in case of zero mode is dealt with in functions
          mod_bes_j1 = mod_besselj1_gkw(imod,ix,i,j,is)

          I_sp1 = signz(is)*de(ix,is)*bn(ix,i)*bes_j0*intvp(i,j,k,is)*intmu(j)
          I_sp2 = veta(ix)*bn(ix,i)*tmp(ix,is)*de(ix,is)*intvp(i,j,k,is)*intmu(j)*   &
  &               mugr(j)*mod_bes_j1

          elem%val = I_sp1*(1.+B_sp2) - I_sp2*B_sp1
          call add_element(elem,ierr)

        end do

      end do ; end do ; end do; end do

    else !bpar is off

      do is = 1, nsp ; do ix = 1, nx ; do i = 1, ns

        do j = 1, nmu ; do k = 1, nvpar

        dum = besselj0_gkw(imod,ix,i,j,is)
        call set_indx(elem,imod,ix,i,j,k,is) 
        elem%val = signz(is)*de(ix,is)*intmu(j)*intvp(i,j,k,is)*dum*bn(ix,i)
        if (.not. gkw_a_equal_b_accuracy(intvp(i,j,k,is), 0.0)) then
          call add_element(elem,ierr)
        endif 

        end do ; end do

      end do ; end do ; end do

    end if
    
  end do

end subroutine poisson_int

!------------------------------------------------------------------------------
!> Adds the part of the Ampere's equation that is related with the
!> integral over the distribution function
!------------------------------------------------------------------------------
subroutine ampere_int

  use structures,     only : matrix_element
  use control,        only : nlapar, spectral_radius
  use grid,           only : nx,ns,nmu,nvpar,nsp,nmod
  use components,     only : de, signz, vthrat, veta
  use matdat,         only : add_element, set_indx
  use geom,           only : bn
  use velocitygrid,   only : vpgr, intvp, intmu
  use dist,           only : ifdis, iapar 
  use functions,      only : besselj0_gkw 

  real                  :: bes
  integer               :: ix, i, j, k, is, imod, ierr 
  type (matrix_element) :: elem 

  ! non spectral case is done elsewhere (in gyro_average)
  if (.not. spectral_radius) return

  ! Indentifier for the term 
  elem%term = 'ampere_int'

  ! Type of iih and jjh 
  elem%itype = iapar 
  elem%itloc = ifdis 
  elem%ideriv = 0
 
  if (.not. nlapar) return

  do imod = 1, nmod

    if (all(apply_on_imod == 0)) then
      if(.not. lampere) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lampere) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
  
    do ix = 1, nx ; do i = 1, ns ; do j = 1, nmu ; do k = 1, nvpar ; do is = 1, nsp
      
      call set_indx(elem,imod,ix,i,j,k,is) 

      ! Bessel function for gyro-average
      bes = besselj0_gkw(imod,ix,i,j,is)

      elem%val = signz(is)*de(ix,is)*veta(ix)*intvp(i,j,k,is)*intmu(j)      &
          &     *vthrat(is)*bn(ix,i)*vpgr(i,j,k,is)*bes

      call add_element(elem,ierr)

    end do ; end do ; end do ; end do ; end do
    
  end do
end subroutine ampere_int

!------------------------------------------------------------------------------
!> This routine adds the integral part of the perpendicular Ampere's law.
!------------------------------------------------------------------------------
subroutine ampere_bpar_int

  use structures,     only : matrix_element
  use grid,           only : nx,ns,nmu,nvpar,nsp,nmod, number_of_species
  use mpicomms,       only : COMM_SP_NE
  use components,     only : de, signz, tmp, adiabatic_electrons, iadia, veta
  use geom,           only : bn
  use matdat,         only : add_element, set_indx
  use velocitygrid,   only : intvp, intmu, mugr
  use control,        only : nlbpar, spectral_radius
  use mode,           only : krloc
  use rotation,       only : cfen
  use dist,           only : ifdis, ibpar 
  use functions,      only : gamma1_gkw, gamma_gkw, mod_besselj1_gkw 
  use functions,      only : besselj0_gkw  
  use mpiinterface,   only : mpiallreduce_sum_inplace

  real :: gamma, gamma_diff, bes_j0, mod_bes_j1
  real :: I_sp1, I_sp2, F_sp1, F_sp2

  integer :: ix, i, j, k, is, imod, ierr
  type (matrix_element) :: elem 

  ! non spectral case is done elsewhere (in gyro_average)
  if (.not. spectral_radius) return

  ! Identifier of the term 
  elem%term = 'ampere_bpar_int' 

  ! type of iih and jjh 
  elem%itype = ibpar 
  elem%itloc = ifdis
  elem%ideriv = 0

  if (nlbpar) then

    do imod = 1, nmod

      if (all(apply_on_imod == 0)) then
        if(.not. lbpar) cycle
      else
        if (any(apply_on_imod == imod)) then
          if(.not. lbpar) cycle
        else
          if(.not.linear_term_switch_default) cycle
        endif
      endif
      
      do ix = 1, nx ; do i = 1, ns ; do j = 1, nmu ; do k = 1, nvpar

        I_sp1 = 0.E0
        I_sp2 = 0.E0
        F_sp1 = 0.E0
        F_sp2 = 0.E0

        !calculating coefficients B_sp1, B_sp2, F_sp1, F_sp2 locally:
        !in case of zero mode the coefficients are simplified since gamma(0)=1
        !and gamma_diff(0)=1.

        do is = 1, nsp

          if (abs(krloc(imod,ix,i)) < 1.0E-5) then !zero mode
            !Gamma functions together with correction due to high rotation
            gamma_diff = 1.0*exp(-cfen(i,is))
            gamma = 1.0*exp(-cfen(i,is))
          else !non-zero modes
            !Gamma functions together with correction due to high rotation
            gamma_diff = (gamma_gkw(imod,ix,i,is) - gamma1_gkw(imod,ix,i,is))*exp(-cfen(i,is))
            gamma = gamma_gkw(imod,ix,i,is)*exp(-cfen(i,is))
          end if

          F_sp1 = F_sp1 + signz(is)**2*de(ix,is)*(gamma-1)/tmp(ix,is)
          F_sp2 = F_sp2 + signz(is)*veta(ix)*de(ix,is)*gamma_diff/(2.0*bn(ix,i))

        end do

        ! calculating coefficients F_sp1, F_sp2 globally
        if (number_of_species > nsp) then

          call mpiallreduce_sum_inplace(F_sp1,1,COMM_SP_NE)
          call mpiallreduce_sum_inplace(F_sp2,1,COMM_SP_NE)

        end if

        !if one of the species is adiabatic there is an extra term in F_sp1
        if (adiabatic_electrons) then 
          F_sp1 = F_sp1 - signz(nsp+iadia)*de(ix,nsp+iadia)/tmp(ix,nsp+iadia)
        end if

        do is = 1, nsp

          !calculating coefficients I_sp1, I_sp2. Their sum over velocity space and species
          !is done in exp_integration.

          bes_j0 = besselj0_gkw(imod,ix,i,j,is)
          !modified j_1 bessel function: mod_bes_j1 = 2*j_1(k_perp rho)/(k_perp rho):
          !the limit in case of zero mode is dealt with in functions
          mod_bes_j1 = mod_besselj1_gkw(imod,ix,i,j,is)

          i_sp1 = signz(is)*de(ix,is)*bn(ix,i)*bes_j0*intvp(i,j,k,is)*intmu(j)
          i_sp2 = veta(ix)*bn(ix,i)*tmp(ix,is)*de(ix,is)*intvp(i,j,k,is)*intmu(j)*   &
&               mugr(j)*mod_bes_j1

          elem%val = i_sp2*f_sp1 - i_sp1*f_sp2

          call set_indx(elem,imod,ix,i,j,k,is) 
          call add_element(elem,ierr)

        end do

      end do ; end do ; end do ; end do
    
    end do

  else
     return
  end if

end subroutine ampere_bpar_int

!-----------------------------------------------------------------------------
!> add the diagonal part of the poisson equation (and the perpendicular
!> ampere's law if bpar is on) into the matrix
!-----------------------------------------------------------------------------
subroutine poisson_dia

  use structures,     only : matrix_element
  use grid,           only : nx,ns,nsp,nmod,number_of_species
  use mpicomms,       only : comm_sp_ne
  use control,        only : nlbpar, spectral_radius
  use mode,           only : krloc, ixzero, iyzero
  use geom,           only : bn
  use components,     only : de, tmp, signz, adiabatic_electrons, iadia, veta
  use matdat,         only : add_element, set_indx
  use rotation,       only : cfen
  use dist,           only : ibpar, iphi 
  use functions,      only : gamma1_gkw, gamma_gkw
  use mpiinterface,   only : mpiallreduce_sum_inplace

  integer :: imod, ix, i, is, idum, ierr
  real    :: gamma, gamma_diff
  type (matrix_element) :: elem, elem2 

  ! coefficients used if bpar is on:
  real :: f_sp1, f_sp2, b_sp1, b_sp2

  ! non spectral case is done elsewhere (in gyro_average - polarisation term)
  if (.not. spectral_radius) return

  ! identifier of the term 
  elem%term = 'poisson_dia'
  elem%ideriv = 0
  ! type of iih and jjh and indices 
  elem%itype = iphi 
  elem%itloc = iphi   

  ! NOTE: Before r4094, the B|| itloc was iphi
  ! which meant the diagonal terms for B|| had the wrong
  ! jj value in the matrix.  This did not cause problems 
  ! since in explicit integration the diagonal terms 
  ! use only the ii value under the assumption they are diagonal
  ! However, it would have been a bug if used by the implict scheme
  elem2%term = 'poisson_bpar'
  elem2%ideriv = 0
  ! type of iih and jjh and indices 
  elem2%itype = ibpar
  elem2%itloc = ibpar

  ! set the dummy 
  idum = 0 

  do imod = 1, nmod

    if (all(apply_on_imod == 0)) then
      if(.not. lpoisson) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lpoisson) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
  
    if (nlbpar) then

      do i = 1, ns; do ix = 1, nx

        f_sp1 = 0.e0
        f_sp2 = 0.e0
        b_sp1 = 0.e0
        b_sp2 = 0.e0

        ! calculating coefficients b_sp1, b_sp2, f_sp1, f_sp2 locally
        ! in case of zero mode the coefficients are simplified since gamma(0)=1
        ! and gamma_diff(0)=1.

        do is = 1, nsp

          if (abs(krloc(imod,ix,i)) < 1e-5) then !zero mode
            ! gamma functions together with correction due to high rotation
            gamma_diff = 1.0*exp(-cfen(i,is))
            gamma = 1.0*exp(-cfen(i,is))
          else ! non-zero modes
            ! gamma functions together with correction due to high rotation
            gamma_diff = (gamma_gkw(imod,ix,i,is) - gamma1_gkw(imod,ix,i,is))*exp(-cfen(i,is))
            gamma = gamma_gkw(imod,ix,i,is)*exp(-cfen(i,is))
          end if

            f_sp1 = f_sp1 + signz(is)**2*de(ix,is)*(gamma-1)/tmp(ix,is)
            b_sp1 = b_sp1 + signz(is)*de(ix,is)*gamma_diff/(bn(ix,i))
            f_sp2 = f_sp2 + signz(is)*veta(ix)*de(ix,is)*gamma_diff/(2.0*bn(ix,i))
            b_sp2 = b_sp2 + tmp(ix,is)*de(ix,is)*veta(ix)*gamma_diff/(bn(ix,i)**2)

        end do

        ! calculating coefficients b_sp1, b_sp2, f_sp1, f_sp2 globally
        ! (The same reductions are repeated and could be moved out of this loop -
        ! but it is has not been demonstrated to date that initialisation ever
        ! takes too long on many cores)
        if (number_of_species > nsp) then

          call mpiallreduce_sum_inplace(b_sp1,1,comm_sp_ne)
          call mpiallreduce_sum_inplace(b_sp2,1,comm_sp_ne)
          call mpiallreduce_sum_inplace(f_sp1,1,comm_sp_ne)
          call mpiallreduce_sum_inplace(f_sp2,1,comm_sp_ne)

        end if

        ! if one of the species is adiabatic there is an extra term in f_sp1
        if (adiabatic_electrons) then 
          f_sp1 = f_sp1 - signz(nsp+iadia)*de(ix,nsp+iadia)                &
                          * exp(-cfen(i,nsp+iadia))/tmp(ix,nsp+iadia)
        end if

        call set_indx(elem,imod,ix,i,idum,idum,idum)

        elem%val = f_sp1*(1+b_sp2) - f_sp2*b_sp1

        if (imod==iyzero.and.ix==ixzero .and. .not. adiabatic_electrons) then
          elem%val = 1.
        end if

        ! put the element to iphi position.  Take the inverse to save
        ! some cycles later because multiplication is faster than
        ! division.
        elem%val = -1.0/elem%val
        call add_element(elem,ierr)

        ! the diagonal part of the perpendicular ampere's law is the same,
        ! if bpar is calculated, put the same element to ibpar position, too.        
        if (all(apply_on_imod == 0)) then
          if(lbpar) then 
            call set_indx(elem2,imod,ix,i,idum,idum,idum)
            elem2%val = elem%val
            call add_element(elem2,ierr)
          endif
        else
          if (any(apply_on_imod == imod)) then
            if(lbpar) then
              call set_indx(elem2,imod,ix,i,idum,idum,idum)
              elem2%val = elem%val
              call add_element(elem2,ierr)
            endif
          else
            if(linear_term_switch_default) then
              call set_indx(elem2,imod,ix,i,idum,idum,idum)
              elem2%val = elem%val
              call add_element(elem2,ierr)
            endif
          endif
        endif

      end do; end do

    else !bpar is off

      do ix = 1, nx ; do i = 1, ns

        ! reference the element of the potential
        call set_indx(elem,imod,ix,i,idum,idum,idum)

        ! initialize the mat_element
        elem%val = (0.,0.)
    
        ! detect the (0,0) mode
        if (imod==iyzero.and.ix==ixzero) then
          ! the (0,0) mode does not contain any turbulent physics, but 
          ! for adiabatic electrons can still be considered for the 
          ! neoclassical physics.
          if (adiabatic_electrons) then 
            elem%val = elem%val + signz(nsp+iadia)* exp(-cfen(i,nsp+iadia))  &
                                            *de(ix,nsp+iadia)/tmp(ix,nsp+iadia)
          else
            elem%val = elem%val + 1.
          endif 
        else
      
            ! sum local species contributions, then non-local
            do is = 1, nsp
              elem%val = elem%val + exp(-cfen(i,is))*                           &
                & (signz(is)**2)*de(ix,is)*(gamma_gkw(imod,ix,i,is)-1)/tmp(ix,is)
          end do

          if (number_of_species > nsp) then
            call mpiallreduce_sum_inplace(elem%val,1,COMM_SP_NE)
          end if
          
          ! add adiabatic electrons contribution
          if (adiabatic_electrons) then
            elem%val = elem%val + signz(nsp+iadia)*exp(-cfen(i,nsp+iadia))  &
                                        *de(ix,nsp+iadia) / tmp(ix,nsp+iadia)
          end if
          
        end if
      
        ! put the inverse element into the matrix. This saves some cycles
        ! later because multiplication is faster than division.
        elem%val = -1.0/elem%val
        call add_element(elem,ierr)

      end do ; end do

    end if
    
  end do

end subroutine poisson_dia

!------------------------------------------------------------------------------
!> add the zonal flow corrections
!------------------------------------------------------------------------------
subroutine poisson_zf

  use control,        only : zonal_adiabatic, spectral_radius
  use grid,           only : nx,ns,nsp,number_of_species, parallel_s
  use mpicomms,       only : COMM_SP_NE, COMM_S_NE
  use mode,           only : iyzero, ixzero
  use components,     only : de, tmp, signz, adiabatic_electrons, iadia
  use matdat,         only : maty, matz
  use matrix_format, only : put_element
  use index_function, only : indx
  use geom,           only : ints
  use rotation,       only : cfen
  use dist,           only : iphi 
  use general,        only : gkw_abort 
  use functions,      only : gamma_gkw
  use mpiinterface,   only : mpiallreduce_sum_inplace
  use structures,     only : matrix_element

  real     :: dum2
  integer  :: ix, i, is, imod
  complex  :: dum_elem
  complex :: val
  ! non spectral case is done elsewhere (in gyro_average - polarisation term)
  if (.not. spectral_radius) return

  ! If zonal_adiabatic = F then no zonal flow correction is used
  if (.not. zonal_adiabatic) return

  ! The correction is of importance only for adiabatic electrons
  if (.not.adiabatic_electrons) return

  if (iadia.eq.0) call gkw_abort('severe error in possion_zf')
    
  ! No krho = 0 mode, i.e. no zonal flow correction
  if (iyzero == 0) return

  ! Only put the correction in the zonal mode
  imod = iyzero

  if (all(apply_on_imod == 0)) then
    if(.not. lpoisson_zf) return
  else
    if (any(apply_on_imod == imod)) then
      if(.not. lpoisson_zf) return
    else
      if(.not.linear_term_switch_default) return
    endif
  endif

  x_grid : do ix = 1, nx

    ! initialize the dummy element
    dum_elem = (0.E0,0.E0)

    s_grid : do i = 1, ns
      ! initialize the matrix element to zero
      val = (0.E0,0.E0)

      do is = 1, nsp
        dum2 = signz(is)*(gamma_gkw(imod,ix,i,is)-1.)*exp(-cfen(i,is))/tmp(ix,is) & 
                         -exp(-cfen(i,nsp+iadia))/tmp(ix,nsp+iadia)
        val = val + signz(is)*de(ix,is)*dum2
      end do

      ! sum all the diagonal contributions
      if (nsp < number_of_species) then
        call mpiallreduce_sum_inplace(val,1,COMM_SP_NE)
      endif

      val = -ints(i) / val
      
      ! put first element in matz
      ! -{ds /A}, summation in exp int
      call put_element(matz, ix, indx(iphi,imod,ix,i), val)
      
      !Also keep {exp(-Ee)/A}
      dum_elem = dum_elem - val*exp(-cfen(i,nsp+iadia))

    end do s_grid

    ! sum all the dum elements over the s-direction?
    if (parallel_s) then
      call mpiallreduce_sum_inplace(dum_elem,1,COMM_S_NE)
    endif

    s_grid2 : do i = 1, ns

      val = tmp(ix,nsp+iadia) / (de(ix,nsp+iadia)*exp(-cfen(i,nsp+iadia)))  &
                 + dum_elem / exp(-cfen(i,nsp+iadia))

      ! Ignore zero zero mode
      if (ix==ixzero) then
        val = (1.E0,0.E0)
      endif

      ! take the inverse to save cycles later: floating point
      ! multiplication is faster than division.
      val = 1.0/val

      ! put element into the matrix maty.
      call put_element(maty,&
         & indx(iphi,imod,ix,i),&
         & ix, &
         & val)
      
    end do s_grid2

  end do x_grid

end subroutine poisson_zf

!-----------------------------------------------------------------------------
!> Add the diagonal part of the Ampere's equation
!-----------------------------------------------------------------------------
subroutine ampere_dia

  use structures,     only : matrix_element
  use control,        only : nlapar, spectral_radius
  use grid,           only : nx,ns,nsp,nmod,nmu,nvpar
  use mpicomms,       only : COMM_S_EQ
  use dist,           only : fmaxwl, iapar 
  use components,     only : de, signz, mas, veta
  use matdat,         only : add_element, set_indx
  use geom,           only : bn
  use mode,           only : krloc
  use rotation,       only : cfen
  use velocitygrid,   only : intvp, intmu
  use functions,      only : gamma_gkw, besselj0_gkw 
  use mpiinterface,   only : mpiallreduce_sum_inplace

  real    :: gamma, gamma_num, b, dum
  integer :: imod, ix, i, j, k, is, idum, ierr
  complex :: mat_elem
  type (matrix_element) :: elem 

  ! non spectral case is done elsewhere (in gyro_average - integral term)
  if (.not. spectral_radius) return

  if (.not. nlapar) return

  ! indentifier of the term 
  elem%term = 'ampere_dia' 
  elem%itype = iapar 
  elem%itloc = iapar
  elem%ideriv = 0

  ! set the dummy 
  idum = 0 

  do imod = 1, nmod 

    if (all(apply_on_imod == 0)) then
      if(.not. lampere) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lampere) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
    
    do ix = 1, nx ; do i = 1, ns

      ! reference Apar
      call set_indx(elem,imod,ix,i,idum,idum,idum)

      ! The nabla^2 term
      mat_elem =  - krloc(imod,ix,i)**2

      ! calculate the Maxwell correction
      dum = 0.
      do is = 1, nsp

        ! The gamma function of the species
        gamma = gamma_gkw(imod,ix,i,is)*exp(-cfen(i,is))

        ! numerical calculation of the gamma function 
        !(includes the strong rotation correction implictly in the maxwellian)
        gamma_num = 0.
        do j = 1, nmu ; do k = 1, nvpar
          b = besselj0_gkw(imod,ix,i,j,is)
          gamma_num = gamma_num + bn(ix,i)*intmu(j)*intvp(i,j,k,is)*b**2*fmaxwl(ix,i,j,k,is)
  !       The implementation below is perhaps more consistent (stability issues?) 
  !       gamma_num = gamma_num + 2.E0*bn(ix,i)*intmu(j)*intvp(i,j,k,is)*b**2* &
  !                 & vpgr(i,j,k,is)**2*fmaxwl(ix,i,j,k,is)
        end do ; end do

        ! The 'Maxwell correction'
        dum = dum - veta(ix)*signz(is)**2*de(ix,is)*gamma_num / mas(is)

      end do !nsp

      ! MPI sum of the Maxwell correction, over species, mu and vpar.
      call mpiallreduce_sum_inplace(dum,1,COMM_S_EQ)

      ! put the inverse element into the matrix. This saves some cycles
      ! later because multiplication is faster than division.
      elem%val = -1.0/(mat_elem + dum)
      call add_element(elem,ierr)

    end do ; end do
    
  end do

end subroutine ampere_dia


!------------------------------------------------------------------------------
!> Routine that implements a boundary damping 
!------------------------------------------------------------------------------
subroutine boundary_damping 

end subroutine boundary_damping


!------------------------------------------------------------------------------
!> This routine adds the neoclassical source terms to the source vector.
!! Term VI in the manual
!!
!! \f$ (1/Z)(T_R E_D {\cal D}^\psi + 2 m_R v_R v_\parallel {\cal H}^\psi
!! \Omega + m_R \Omega^2 {\cal I}^\psi)*
!! (1/L_n + E_T 1/L_T + 2 v_par u_{\parallel}^\prime / v_R)    \f$
!!
!<-----------------------------------------------------------------------------
subroutine neoclassical

  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod
  use geom,           only : ints, bn
  use dist,           only : fmaxwl, ifdis, f_EP, df_EPdv
  use index_function, only : indx
  use mode,           only : ixzero, iyzero
  use matdat,         only : put_source
  use velocitygrid,   only : intvp, intmu, vpgr, mugr
  use geom,           only : ints, bn, bt_frac, rfun
  use geom,           only : ffun, gfun
  use rotation,       only : coriolis, cf_drift
  use mpiinterface,   only : mpiallreduce_sum_inplace, root_processor
  use mpicomms,       only : COMM_SP_EQ
  use control,        only : spectral_radius
  use neoequil,       only : gradneof, neo_nsp
  use components,     only : types, tgrid, tmp, vthrat, tp, vpar_mean
  use constants,      only : pi

  !character(len=64) :: term='VI: neoclassical'
  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is
  complex :: mat_elem
  ! Dummy arguments for integration. While the result should be real, the
  ! variables are complex, to be able to verify this.
  complex :: dum1, dum2, dum3
  real    :: dum
  integer :: iih
  real    :: drift_x, drift_y, drift_z
  
  dum1 = (0.E0, 0.E0)
  dum2 = (0.E0, 0.E0)
  dum3 = (0.E0, 0.E0)

  ! calculate the terms due to the Maxwell background
  do imod = 1, nmod

    if (all(apply_on_imod == 0)) then
      if(.not. lneoclassical) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lneoclassical) cycle
      else
        if(.not.linear_term_switch_default) cycle
      endif
    endif
    
    do is=1,nsp ; do ix=1,nx ; do i=1,ns ; do j=1,nmu ; do k=1,nvpar

      ! Neoclassical terms only for the (0,0) mode in the spectral case, and only for 
      ! the imod = iyzero mode for the nonspectral case 
      if ( (spectral_radius.and.(ixzero==ix).and.(iyzero==imod)) .or. &
        &  ((.not.spectral_radius).and.(iyzero==imod)) ) then

        ! calculate the drift
        if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneorotsource ) .or. & 
        & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
          !This version simply removes the coriolis drift from the source drift as
          !per Hinton Wong derivation
          call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,.false.,cf_drift,.true.) 
        else
          call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
        end if

        ! the matrix element 
        iih = indx(ifdis,imod,ix,i,j,k,is)
        if (types(is) .eq. 'EP') then
          mat_elem = dmaxwel(ix,i,j,k,is,imod)*f_EP(ix,i,j,k,is)*drift_x +          &
          & 0.5E0*(                                                           &
          &  (vpgr(i,j,k,is)-vpar_mean)**2/(tmp(ix,is)/tgrid(is))*            &
          &     exp(-(vpgr(i,j,k,is)-vpar_mean)**2/(tmp(ix,is)/tgrid(is)))+   &
          &  (vpgr(i,j,k,is)+vpar_mean)**2/(tmp(ix,is)/tgrid(is))*            &
          &     exp(-(vpgr(i,j,k,is)+vpar_mean)**2/(tmp(ix,is)/tgrid(is)))    &
          &  )                                                                &
          & * exp(-2.E0*bn(ix,i)*mugr(j)/(tmp(ix,is)/tgrid(is)))              &
          & *tp(ix,is)*0.5E0/(sqrt(tmp(ix,is)*pi/tgrid(is))**3)*drift_x
                
        else      
          mat_elem = dmaxwel(ix,i,j,k,is,imod)*fmaxwl(ix,i,j,k,is)*drift_x
        endif

        if (types(is) .eq. 'EP') then
          !In the case of a non maxwellian one of the terms that get added is
          ! 2*vthrat*mu/T * b \cdot \gradB *(Feq + T/(2*vpar)* dFeqdvpar)
          !This term comes from the sum of vpar grad_Feq and vpar mu gradB dFeqdvpar
          dum = 2.0 * vthrat(is) * mugr(j) / (tmp(ix,is)/tgrid(is)) *&
          & bn(ix,i) * ffun(ix,i) * gfun(ix,i) * ( f_EP(ix,i,j,k,is) + &
          & (tmp(ix,is)/tgrid(is)) / (2* vpgr(i,j,k,is)) * df_EPdv(ix,i,j,k,is))
          mat_elem = mat_elem + dum
        write(*,*) dum
        endif

        if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneo_equil ) .or. & 
        & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
          if(is <= neo_nsp)then
            !Adds the radial gradient of the neoclassical correction 
            !to the maxwellian background distribution function
            !Not 100% sure it makes ordering sense.
            mat_elem = mat_elem + drift_x*gradneof(is,i,j,k)
            !TO DO:  The s derivative also?
          endif
        endif
    
        !To avoid confusion.  These are just integral checks and are not
        !placed into the matrix.  The reduction is done below and output.
        dum1 = dum1 + mat_elem*intvp(i,j,k,is)*intmu(j)*ints(i)*bn(ix,i)
        dum2 = dum2 + mat_elem*rfun(ix,i)*bt_frac(ix,i)*vpgr(i,j,k,is)* &
          & intvp(i,j,k,is)*intmu(j)*ints(i)*bn(ix,i)
        dum3 = dum3 + mat_elem*(vpgr(i,j,k,is)**2 + 2.E0*bn(ix,i)*mugr(j))* &
          &  intvp(i,j,k,is)*intmu(j)*ints(i)*bn(ix,i)

        ! put the element 
        call put_source(iih,mat_elem)

      endif
    

    end do ; end do ; end do ; end do  ; end do

  end do

  call mpiallreduce_sum_inplace(dum1,1,COMM_SP_EQ)
  call mpiallreduce_sum_inplace(dum2,1,COMM_SP_EQ)
  call mpiallreduce_sum_inplace(dum3,1,COMM_SP_EQ)
  if(root_processor)then
    write(*,*)'Neoclassical source integrals', dum1, dum2, dum3
  endif

  call calc_correction_fluxes

end subroutine neoclassical

!-----------------------------------------------------------------------------
!> This routine adds the neoclassical source terms as outlined by Hinton and Wong
!! WORK IN PROGRESS
!<-----------------------------------------------------------------------------
subroutine neoclassical_rot

  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod
  use dist,           only : fmaxwl, ifdis
  use index_function, only : indx
  use mode,           only : ixzero, iyzero
  use matdat,         only : put_source
  use velocitygrid,   only : vpgr, mugr
  use geom,           only : bn, gfun, bt_frac, rfun
  use components,     only : signz, tmp, tgrid, tp, vp
  use rotation,       only : cfen, dcfen_ds, toroidal_shear, ts_uprim
  
!  character(len=64) :: term='VI: neoclassical_rot'
  
  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is
  complex :: mat_elem
  real :: dum, dum2, dum3, vel, ET
 
  ! dummy variables
  integer :: iih
  real    :: up_prim, prefac1

  !elem%itloc = ifdis
  !elem%term  = term
  dum = 0.E0
  dum2 = 0.E0
  dum3 = 0.E0

  ! calculate the terms due to the Maxwell background
  do is=1,nsp ; do imod = 1, nmod ; do ix=1,nx ; do i=1,ns ; do j=1,nmu ; do k=1,nvpar

    ! Neoclassical terms only for the (0,0) mode
    if (ixzero==ix.and.iyzero==imod) then

      ! the matrix element 
      iih = indx(ifdis,imod,ix,i,j,k,is)
      
      ! term in front of the temperature gradient
      ET = (vpgr(i,j,k,is)**2 + mugr(j)*bn(ix,i)) / (tmp(ix,is) / tgrid(is) ) &
      &  - 1.5E0 + cfen(i,is)
      
      if (toroidal_shear.eq.'use_shear_rate') then
        up_prim=bt_frac(ix,i)*ts_uprim * rfun(ix,i)
      else !vp is uprim
        up_prim=bt_frac(ix,i)*vp(ix,is) * rfun(ix,i)
      end if

      vel = vpgr(i,j,k,is)**2 + bn(ix,i)*mugr(j)
      prefac1 = bt_frac(ix,i) * rfun(ix,i)/(signz(is)*tmp(ix,is))
      
      ! return the final value
      !The prefactor before the temperature gradient
      mat_elem = prefac1*ET*tp(ix,is)*(vel*gfun(ix,i) + dcfen_ds(i,is))
      !Add the prefactor before the rotation gradient term1
      mat_elem = mat_elem + prefac1*up_prim*vpgr(i,j,k,is)*(vel*gfun(ix,i) + dcfen_ds(i,is))
      !Note there is no density (or pressure) gradient contribution as
      !that is removed by the transformation.  This gives only two drives
      !terms.     
      mat_elem = -mat_elem*fmaxwl(ix,i,j,k,is)

      ! put the element 
      call put_source(iih,mat_elem)

    endif

  end do ; end do ; end do ; end do  ; end do ; end do

end subroutine neoclassical_rot

!------------------------------------------------------------------------------
!>If Hinton Wong transformation is used, this section calculates the correction 
!!in the fluxes due to the transformation
!<------------------------------------------------------------------------------
subroutine calc_correction_fluxes

  use grid,           only : nx, ns, nmu, nvpar, nsp, nmod
  use grid,           only : gsp, gx
  use geom,           only : ints, bn, bt_frac, jfun, rfun, jfun
  use dist,           only : fmaxwl
  use velocitygrid,   only : intvp, intmu, vpgr, mugr
  use geom,           only : ints, bn, signB
  use rotation,       only : coriolis, cf_drift, vcor, cf_trap, cfen
  use components,     only : vthrat, signz
  use mpiinterface,   only : mpiallreduce_sum_inplace, root_processor
  use mpicomms,       only : COMM_SP_EQ

  ! integers for the loop over all grid points
  integer :: imod, ix, i, j, k, is
  real :: dum, dum3
  real :: dum_cor, dum2_cor, dum3_cor
  real :: sum1, sum2, sum3
  integer :: isglb, ixg
  real    :: drift_x, drift_y, drift_z, dumnc
  logical :: is_in_effect

  sum1 = 0.E0
  sum2 = 0.E0
  sum3 = 0.E0
  is_in_effect = .false.

  nmod1: do imod = 1, nmod
    
    if (all(apply_on_imod == 0)) then
      if(.not. lneorotsource) cycle
    else
      if (any(apply_on_imod == imod)) then
        if(.not. lneorotsource) cycle
      else
        if(.not.lneo_equil_switch_default) cycle
      endif
    endif
    
    is_in_effect = .true.
    
    nx1:   do ix = 1, nx 
      nsp1:  do is = 1, nsp
 
       ! the actual (global) species index 
       isglb = gsp(is)
       ! the acutal (global) x index
       ixg = gx(ix)!
        
       !Neoclassical fluxes old version
       ! check if this is the 0,0 mode for which the neoclassical
       ! fluxes are calculated 
       ! FJC_NON_SPECTRAL: ix /= ixzero has no meaning
       ns1: do i = 1, ns 
           nmu1:do j = 1, nmu 
             nvpar1: do k = 1, nvpar 
              
               call drift(ix,i,j,k,is,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)                
               ! common factors 
               dumnc = drift_x * (bn(ix,i)*intmu(j)*intvp(i,j,k,is)/signz(is))*     &
                     & fmaxwl(ix,i,j,k,is)*ints(i)

               dum_cor = signB*vpgr(i,j,k,is)*rfun(ix,i)*bt_frac(ix,i) + vcor*jfun(ix,i)
               dum2_cor = vcor*jfun(ix,i)*(vpgr(i,j,k,is)**2 +  2.E0*mugr(j)*bn(ix,i) - 5.E0/2.E0)
               dum3_cor = (vcor*jfun(ix,i)/vthrat(is)**2)*(vpgr(i,j,k,is)* &
                     & rfun(ix,i)*bt_frac(ix,i) + vcor*jfun(ix,i))
               
               dum = (vpgr(i,j,k,is)**2 +  2.E0*mugr(j)*bn(ix,i) - 5.E0/2.E0)
               !Correction for when the centrifugal drift is included.
               if(cf_drift.or.cf_trap)then
                 dum = dum + cfen(i,is)
               endif
               dum3 = vpgr(i,j,k,is)*Rfun(ix,i)*bt_frac(ix,i)*signB

               sum1 = sum1 + dum_cor*dumnc
               sum1 = sum1 + dum2_cor*dumnc
               sum1 = sum1 + dum3_cor*dumnc

               sum2 = sum2 + dum_cor*dumnc*dum
               sum2 = sum2 + dum2_cor*dumnc*dum
               sum2 = sum2 + dum3_cor*dumnc*dum

               sum3 = sum3 + dum_cor*dumnc*dum3
               sum3 = sum3 + dum2_cor*dumnc*dum3
               sum3 = sum3 + dum3_cor*dumnc*dum3
    
             end do nvpar1 
           end do nmu1
         end do ns1          

       end do nsp1 
     end do nx1
   end do nmod1

  if(is_in_effect) then
    call mpiallreduce_sum_inplace(sum1,1,COMM_SP_EQ)
    call mpiallreduce_sum_inplace(sum2,1,COMM_SP_EQ)
    call mpiallreduce_sum_inplace(sum3,1,COMM_SP_EQ)
    if(root_processor)then
      write(*,*) 'Neo correction - Particle',sum1
      write(*,*) 'Neo correction - Energy',sum2
      write(*,*) 'Neo correction - Parallel Momentum',sum3
    endif
  endif

end subroutine calc_correction_fluxes

!------------------------------------------------------------------------------
! The subroutine below calculates the drift due to the magnetic field 
! inhomogeneity as well as the plasma rotation 
!------------------------------------------------------------------------------
subroutine drift(ix,i,j,k,is,drift_x,drift_y,drift_z,lcoriolis,lcf_drift,lgradb) 

  use components,     only : signz, vthrat, mas, tgrid, veta_prime
  use geom,           only : dfun, efun, hfun, bn, ifun, dpfdpsi, signB,    &
                      &  gradp_type, dpds_rot, dpdpsi_rot, beta_miller,     &
                      &  curv_effect, de_miller, mas_miller_i
  use velocitygrid,   only : vpgr, mugr
  use rotation,       only : vcor!, coriolis, cf_drift
  use rotation,       only : dcfphi_ds, dcfphi_dpsi

  logical, intent(in)  :: lcoriolis, lcf_drift, lgradb
  ! local index for radial-, s-, mu-, vpar-direction and species, respectively.
  integer, intent(in)  :: ix,i,j,k,is
  real,    intent(out) :: drift_x, drift_y, drift_z

  real :: ED, grdp1, grdp2, grdp3

  ED = vpgr(i,j,k,is)**2
  if(lgradb)then
    ED = ED + bn(ix,i)*mugr(j)
  endif

  ! the B\times \nabla B component of the drift
  drift_x = tgrid(is)*ED*dfun(ix,i,1)
  drift_y = tgrid(is)*ED*dfun(ix,i,2) 
  drift_z = tgrid(is)*ED*dfun(ix,i,3) 

  
  if (gradp_type == 'rota_miller' .and. curv_effect) then
    ! change for vcor =/= 0 (used with miller)
   
    ! Derivative of the pressure toward psi            
    grdp1 = 2.*dpdpsi_rot(ix,i) * dpfdpsi(ix)   

    ! Derivative of the pressure toward s
    grdp2 = dpds_rot(ix,i)*signB

    ! Additionnal term due to non zero toroidal velocities
    grdp3 = mas_miller_i*de_miller(ix,1)*vcor**2.*signB*beta_miller(ix)

    !Terms due to toroidal rotation in pressure gradient
    drift_x = drift_x + tgrid(is)*vpgr(i,j,k,is)**2*grdp1*efun(ix,i,1,1)/ bn(ix,i)**2
    drift_y = drift_y + tgrid(is)*vpgr(i,j,k,is)**2*grdp1*efun(ix,i,1,2)/ bn(ix,i)**2
    drift_z = drift_z + tgrid(is)*vpgr(i,j,k,is)**2*grdp1*efun(ix,i,1,3)/ bn(ix,i)**2

    ! Terms due to modification of force balance equation for a 
    ! two fluid model
    drift_x = drift_x + tgrid(is)*vpgr(i,j,k,is)**2*(grdp3-grdp2)*ifun(ix,i,1)/ bn(ix,i)**2
    drift_y = drift_y + tgrid(is)*vpgr(i,j,k,is)**2*(grdp3-grdp2)*ifun(ix,i,2)/ bn(ix,i)**2
    drift_z = drift_z + tgrid(is)*vpgr(i,j,k,is)**2*(grdp3-grdp2)*ifun(ix,i,3)/ bn(ix,i)**2

   else

    ! The finite beta correction of the curvature
    drift_x = drift_x + tgrid(is)*vpgr(i,j,k,is)**2*veta_prime(ix)*efun(ix,i,1,1)/ bn(ix,i)**2
    drift_y = drift_y + tgrid(is)*vpgr(i,j,k,is)**2*veta_prime(ix)*efun(ix,i,1,2)/ bn(ix,i)**2
    drift_z = drift_z + tgrid(is)*vpgr(i,j,k,is)**2*veta_prime(ix)*efun(ix,i,1,3)/ bn(ix,i)**2  

   end if   

  ! The coriolis drift correction (uses m_R v_R = T_R / v_R)
  if (lcoriolis) then 
    drift_x = drift_x + 2.E0*mas(is)*vthrat(is)*vpgr(i,j,k,is)*vcor*hfun(ix,i,1)
    drift_y = drift_y + 2.E0*mas(is)*vthrat(is)*vpgr(i,j,k,is)*vcor*hfun(ix,i,2)
    drift_z = drift_z + 2.E0*mas(is)*vthrat(is)*vpgr(i,j,k,is)*vcor*hfun(ix,i,3)
  endif 

  ! The centrifugal drift
  if (lcf_drift) then 
    drift_x = drift_x + vcor*vcor*mas(is)*ifun(ix,i,1)
    drift_y = drift_y + vcor*vcor*mas(is)*ifun(ix,i,2)
    drift_z = drift_z + vcor*vcor*mas(is)*ifun(ix,i,3)
  endif   

  ! common factor (1/Z)
  drift_x = drift_x / signz(is)
  drift_y = drift_y / signz(is)
  drift_z = drift_z / signz(is)

  ! drift from centrifugal potential
  drift_x = drift_x + efun(ix,i,1,1)*dcfphi_dpsi(i)+efun(ix,i,3,1)*dcfphi_ds(i)
  drift_y = drift_y + efun(ix,i,1,2)*dcfphi_dpsi(i)+efun(ix,i,3,2)*dcfphi_ds(i)
  drift_z = drift_z + efun(ix,i,1,3)*dcfphi_dpsi(i)+efun(ix,i,3,3)*dcfphi_ds(i)

end subroutine drift 


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> For the arakawa scheme
!> \f$   (1/2) v_{\parallel}^2 + \mu B_N + (1/2) {\cal E}_R  \f$ 
!> Note H is periodic in n_s_grid
!----------------------------------------------------------------------------
function HH(ix, i,j,k,is)

  use geom, only : bn
  use velocitygrid, only : vpgr, mugr
  use rotation, only : cfen

  real :: HH
  integer, intent(in) :: ix,i,j,k,is

  !cfen and bn are periodic in n_s_grid and have ghost points as needed.
  HH = 0.5*vpgr(i,j,k,is)**2 + mugr(j)*bn(ix,i) + 0.5*cfen(i,is)

end function HH


!------------------------------------------------------------------------------
!> This routine determines the differential scheme used. The selection is 
!> based on the order of the scheme, the position of the point on the  
!> field line (i.e. end point or not), and up or down wind differentiation 
!------------------------------------------------------------------------------
subroutine differential_scheme(ist,ipw,id,w,optscheme) 

  use control,  only : order_of_the_scheme
  use global,   only : lenswitch
  use general,  only : gkw_abort
  !> an integer indicating if the position is near the grid end
  integer, intent(in) :: ist
  !> +1 / -1 up / down wind
  integer, intent(in) :: ipw
  !> id=1 means 1st deriv
  !> id=2 means dissipation term (resulting in 2nd or 4th deriv)
  integer, intent(in) :: id
  ! a buffer large enough to hold the finite difference coefficients of a
  ! completely one-sided stencil, for the chosen order of the scheme
  real,     intent(out) :: w(:) 
  character (len = lenswitch), optional :: optscheme
  character (len = lenswitch) :: scheme_order

  real, dimension(3) :: ord1_d1_r1
  real, dimension(3) :: ord1_d1_l1
  real, dimension(5) :: ord2_null
  real, dimension(5) :: ord2_d1_c
  real, dimension(5) :: ord2_d2_c
  real, dimension(5) :: ord2_d1_c_ghost2
  real, dimension(5) :: ord2_d1_c_ghost4
  real, dimension(5) :: ord2_d1_r2
  real, dimension(5) :: ord2_d2_c_ghost2
  real, dimension(5) :: ord2_d4_c
  real, dimension(5) :: ord2_d2_c_ghost4
  real, dimension(5) :: ord2_d1_l2
  real, dimension(9) :: ord3_null
  real, dimension(9) :: ord3_d1_l1
  real, dimension(9) :: ord3_d1_r1
  real, dimension(9) :: ord4_d1_c_ghost3
  real, dimension(9) :: ord4_d2_c_ghost3
  real, dimension(9) :: ord4_d1_c
  real, dimension(9) :: ord4_d2_c
  real, dimension(9) :: ord4_d1_c_ghost7
  real, dimension(9) :: ord4_d2_c_ghost7
  real, dimension(13) :: ord6_d1_c
  real, dimension(13) :: ord6_d2_c
  
  
  !---- First order derivatives ------
  ! 1st order, 1st derivative, right-leaning by 1 cell
  ! ('upwind' and 'downwind' refer to flow direction, not grid direction)
  data ord1_d1_r1 /  0.E0, -12.E0,  12.E0 /
  ! 1st order, 1st derivative, left-leaning by 1 cell
  data ord1_d1_l1 /-12.E0,  12.E0,   0.E0 /
  
  !---- Second order derivatives ------
  data ord2_null /  0.E0,   0.E0,   0.E0,   0.E0,   0.E0 / 
  ! 2nd order, 1st deriv, central difference
  data ord2_d1_c   / 0.E0,   -6.E0,   0.E0,   6.E0,   0.E0 / 
  ! 2nd order, 2nd deriv, central difference
  data ord2_d2_c   / 0.E0,   12.E0, -24.E0,  12.E0,   0.E0 /
  ! 2nd order, 1st deriv, central difference, with one cell zeroed
  data ord2_d1_c_ghost2 / 0.E0,   0.E0,   0.E0,   6.E0,   0.E0 / 
  data ord2_d1_c_ghost4 / 0.E0,   -6.E0,   0.E0,   0.E0,   0.E0 /
  ! 2nd order, 1st derivative, right-leaning by 2 cells
  data ord2_d1_r2 / 0.E0,   0.E0, -18.E0,  24.E0,  -6.E0 /
  data ord2_d2_c_ghost2 / 0.E0,   0.E0, -24.E0,  12.E0,   0.E0 /
  data ord2_d4_c   /-1.E0,   4.E0,  -6.E0,   4.E0,  -1.E0 /
  data ord2_d2_c_ghost4 / 0.E0,  12.E0, -24.E0,   0.E0,   0.E0 / 
  data ord2_d1_l2 / 6.E0, -24.E0,  18.E0,   0.E0,   0.E0 / 

  !---- Third order scheme -------
  data ord3_null / 0.E0,   0.E0,   0.E0,   0.E0,   0.E0,   0.E0,   0.E0,   0.E0,   0.E0 /
  data ord3_d1_l1 / 0.E0,   0.E0,   2.E0, -12.E0,   6.E0,   4.E0,   0.E0,   0.E0,   0.E0 /
  data ord3_d1_r1 / 0.E0,   0.E0,   0.E0,  -4.E0,  -6.E0,  12.E0,  -2.E0,   0.E0,   0.E0 / 
  
  !---- Fourth order scheme -------
  data ord4_d1_c_ghost3 / 0.E0,   0.E0,   0.E0,  -8.E0,   0.E0,   8.E0,  -1.E0,   0.E0,   0.E0 /
  data ord4_d2_c_ghost3 / 0.E0,   0.E0,   0.E0,   4.E0,  -6.E0,   4.E0,  -1.E0,   0.E0,   0.E0 /
  data ord4_d1_c  / 0.E0,   0.E0,   1.E0,  -8.E0,   0.E0,   8.E0,  -1.E0,   0.E0,   0.E0/
  data ord4_d2_c  / 0.E0,   0.E0,  -1.E0,  16.E0, -30.E0,  16.E0,  -1.E0,   0.E0,   0.E0/
  data ord4_d1_c_ghost7 / 0.E0,   0.E0,   1.E0,  -8.E0,   0.E0,   8.E0,   0.E0,   0.E0,   0.E0 / 
  data ord4_d2_c_ghost7 / 0.E0,   0.E0,  -1.E0,   4.E0,  -6.E0,   4.E0,   0.E0,   0.E0,   0.E0 / 
  
  !---- Sixth order scheme -------
  ord6_d1_c  =  (/ 0.E0, 0.E0, 0.E0, -1.E0,   9.E0, -45.E0,   0.E0,  45.E0,  -9.E0,   1.E0, 0.E0, 0.E0, 0.E0 /) / 5.0E0
  ord6_d2_c  =  (/ 0.E0, 0.E0, 0.E0, 1.E0,   -6.E0,  15.E0, -20.E0,  15.E0,  -6.E0,   1.E0, 0.E0, 0.E0, 0.E0 /) / 5.0E0
  
  scheme_order = order_of_the_scheme
  if (present(optscheme)) scheme_order = optscheme

  w = 0

  if (id.eq.1) then
    ! First derivatives

    select case(scheme_order) 
    
    case('sixth_order')

      ! decrease the order near the border, to reduce reflections.

      ! put the right stencil into w, around the centre.
      select case(ist)
      case(-2)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      case(-1)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      case(0)
        if (ipw.eq.1) w = put_centered(w,ord6_d1_c)
        if (ipw.eq.-1) w = put_centered(w,ord6_d1_c)
      case(1)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      case(2)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      end select

    case('fourth_order')

      ! decrease the order near the border, to reduce reflections
      select case(ist) 
      case(-2) 
        if (ipw.eq.1) w = put_centered(w,ord2_d1_r2)
        if (ipw.eq.-1) w = put_centered(w,ord2_d1_c_ghost2)
      case(-1)
        if (ipw.eq.1) w = put_centered(w,ord3_d1_r1)
        if (ipw.eq.-1) w = put_centered(w,ord4_d1_c_ghost3) 
      case(0)
        if (ipw.eq.1) w = put_centered(w,ord4_d1_c)
        if (ipw.eq.-1) w = put_centered(w,ord4_d1_c) 
      case(1)
        if (ipw.eq.1) w = put_centered(w,ord4_d1_c_ghost7)
        if (ipw.eq.-1) w = put_centered(w,ord3_d1_l1) 
      case(2) 
        if (ipw.eq.1) w = put_centered(w,ord2_d1_c_ghost4)
        if (ipw.eq.-1) w = put_centered(w,ord2_d1_l2) 
      end select 

    case('second_order') 
  
      select case(ist) 
      case(-2) 
        if (ipw.eq.1) w = put_centered(w,ord1_d1_r1)
        if (ipw.eq.-1) w = put_centered(w,ord2_d1_c_ghost2) 
      case(-1,0,1) 
        if (ipw.eq.1) w = put_centered(w,ord2_d1_c)
        if (ipw.eq.-1) w = put_centered(w,ord2_d1_c) 
      case(2) 
        if (ipw.eq.1) w = put_centered(w,ord2_d1_c_ghost4)
        if (ipw.eq.-1) w = put_centered(w,ord1_d1_l1) 
      end select 
      
    case default
    
       call gkw_abort('linear_terms: unknown differential scheme_order')  

    end select
    
  else if (id == 2) then     ! standard numerical (hyper)dissipation 

    select case(scheme_order) 
    
    case('sixth_order')
    
      select case(ist)
      case(-2)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      case(-1)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      case(0)
        if (ipw.eq.1) w = put_centered(w,ord6_d2_c)
        if (ipw.eq.-1) w = put_centered(w,ord6_d2_c)
      case(1)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      case(2)
        call gkw_abort('linear_terms: sixth order not yet implemented here')
      end select

    case('fourth_order')     ! 4th derivatives at 2nd order

      select case(ist) 
      case(-2) 
        if (ipw.eq.1) w = put_centered(w,ord2_null)
        if (ipw.eq.-1) w = put_centered(w,ord2_d2_c_ghost2) 
      case(-1)
        if (ipw.eq.1) w = put_centered(w,ord3_null)
        if (ipw.eq.-1) w = put_centered(w,ord4_d2_c_ghost3) 
      case(0)
        if (ipw.eq.1) w = put_centered(w,ord2_d4_c)
        if (ipw.eq.-1) w = put_centered(w,ord2_d4_c) 
      case(1)
        if (ipw.eq.1) w = put_centered(w,ord4_d2_c_ghost7)
        if (ipw.eq.-1) w = put_centered(w,ord3_null) 
      case(2) 
        if (ipw.eq.1) w = put_centered(w,ord2_d2_c_ghost4)
        if (ipw.eq.-1) w = put_centered(w,ord2_null) 
      end select 

    case('second_order')     ! 2nd derivatives at 4th order
 
      select case(ist) 
      case(-2) 
        if (ipw.eq.1) w = put_centered(w,ord2_null)
        if (ipw.eq.-1) w = put_centered(w,ord2_null) 
      case(-1,0,1) 
        if (ipw.eq.1) w = put_centered(w,ord2_d2_c)
        if (ipw.eq.-1) w = put_centered(w,ord2_d2_c) 
      case(2) 
        if (ipw.eq.1) w = put_centered(w,ord2_null)
        if (ipw.eq.-1) w = put_centered(w,ord2_null) 
      end select
      
    case default
    
      call gkw_abort('linear_terms: unknown differential scheme_order')  

    end select
    
  else if (id == 3) then     ! second derivative, regardless of the order
    ! this can be used in diagnostics
    select case(scheme_order)
    case('fourth_order')
      select case(ist) 
      case(0) 
        if (ipw.eq.1) w = put_centered(w,ord4_d2_c) ! 4th order 2nd deriv
        if (ipw.eq.-1) w = put_centered(w,ord4_d2_c)
      case default
        call gkw_abort('linear_terms: differential scheme for second derivative &
           & and ist/=0 is not chosen.')
      end select
    case('second_order')
      select case(ist) 
      case(0) 
        if (ipw.eq.1) w = put_centered(w,ord2_d2_c) ! 2nd order 2nd deriv
        if (ipw.eq.-1) w = put_centered(w,ord2_d2_c)
      case default
        call gkw_abort('linear_terms: differential scheme for second derivative &
           & and ist/=0 is not chosen.')
      end select
    end select
    
  else if (id == -2) then  
  
    select case(scheme_order)  ! alternative numerical dissipation

    case('fourth_order')  ! 2nd derivatives at 4th order
  
      select case(ist) 

        case(-2)  
          if (ipw.eq.1) w = put_centered(w,ord2_null)
          if (ipw.eq.-1) w = put_centered(w,ord2_d2_c_ghost2)
        case(-1)      
          if (ipw.eq.1) w = put_centered(w,ord3_null)
          if (ipw.eq.-1) w = put_centered(w,ord4_d2_c_ghost3)
        case(0)  
          if (ipw.eq.1)  w = put_centered(w,ord4_d2_c)
          if (ipw.eq.-1) w = put_centered(w,ord4_d2_c)
        case(1)       
          if (ipw.eq.1) w = put_centered(w,ord4_d2_c_ghost7)
          if (ipw.eq.-1) w = put_centered(w,ord3_null)
        case(2)            
          if (ipw.eq.1) w = put_centered(w,ord2_d2_c_ghost4)
          if (ipw.eq.-1) w = put_centered(w,ord2_null)

      end select    
        
    case('second_order')
    
       call gkw_abort('linear_terms: Use positive idisp for second_order scheme')   
    
    case default
    
       call gkw_abort('linear_terms: unknown differential scheme_order')  

    end select 
    
  else  
  
    call gkw_abort('linear_terms: unknown differential derivative')  

  end if 

  ! Do the proper normalization 
  w = w / 12.E0  

contains

  !-----------------------------------------------------------------------------
  !> copy the stencil array into w, so that it is centered.
  !-----------------------------------------------------------------------------
  function put_centered(w, stencil) result(wout)
    use general, only : gkw_abort
    real, dimension(:), intent(in) :: w, stencil
    real, dimension(size(w)) :: wout
    integer :: delta
    delta = (size(w)-size(stencil))/2
    if(delta < 0) call gkw_abort("put_centered: delta < 0, probably w(:) is too small.")
    wout = w
    wout(1+delta:size(w)-delta) = stencil
  end function put_centered

end subroutine differential_scheme



!-----------------------------------------------------------------------------
!> This routine calculates the correction necessary to go from the
!> distribution g (which includes the correction of the parallel
!> vector potential) to the distribution f
!-----------------------------------------------------------------------------
subroutine g2f_correct

  use control,        only : spectral_radius
  use index_function, only : indx,get_block_bounds,match_indices_with_dims
  use index_function, only : IS_6D_DISTR
  use dist,           only : fmaxwl, iapar, iapar_ga
  use dist,           only : df_EPdv
  use components,     only : signz, vthrat, tmp, types
  use matdat,         only : matg2f
  use matrix_format, only : put_element
  use velocitygrid,   only : vpgr
  integer :: jjh
  real :: b0
  complex :: mat_elem

  integer, dimension(6) :: starts,ends
  integer :: i,j,k,l,m,n
  integer :: i_mod, i_x, i_s, i_mu, i_vpar, i_sp

  integer :: ii
  ii = 0
  
  call get_block_bounds(starts,ends,is_distr_or_field=IS_6D_DISTR)
  do n=starts(6),ends(6)
    do m=starts(5),ends(5)
      do l=starts(4),ends(4)
        do k=starts(3),ends(3)
          do j=starts(2),ends(2)
            do i=starts(1),ends(1)
              call match_indices_with_dims(i,j,k,l,m,n,i_mod,i_x,i_s,i_mu,i_vpar,i_sp)
              b0 = wrap_beslj0_gkw(i_mod,i_x,i_s,i_mu,i_sp)
              mat_elem = -2.E0*signz(i_sp)*vthrat(i_sp)*vpgr(i_s,i_mu,i_vpar,i_sp)*&
                       &  b0*fmaxwl(i_x,i_s,i_mu,i_vpar,i_sp)/tmp(i_x,i_sp)
              if (types(i_sp) .eq. 'EP') then
                mat_elem = signz(i_sp)*vthrat(i_sp)/tmp(i_x,i_sp)*b0* &
                         & df_EPdv(i_x,i_s,i_mu,i_vpar,i_sp)
              endif
              
              if (all(apply_on_imod == 0)) then
                if(.not. lg2f_correction) mat_elem = 0.0
              else
                if (any(apply_on_imod == i_mod)) then
                  if(.not. lg2f_correction) mat_elem = 0.0
                else
                  if(.not.linear_term_switch_default) mat_elem = 0.0
                endif
              endif
              
              if (spectral_radius) then 
                jjh = indx(iapar,i_mod,i_x,i_s)
              else 
                jjh = indx(iapar_ga,i_mod,i_x,i_s,i_mu,i_sp) 
              endif
              ii = ii+1
              call put_element(matg2f,ii,jjh,mat_elem)
            end do
          end do
        end do
      end do
    end do
  end do
              
end subroutine g2f_correct

! !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! THE ROUTINE BELOW COULD BE REMOVED AFTER SOME CODING 
!
! DIFFUS : used for dissipation in Arakawa scheme
!          could be replaced by the implementation of the parallel and velocity
!          space dissipation which should be the same 
! 
! !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!------------------------------------------------------------------------------
!> add 4th or 2nd order diffusion in s or vpar (dx = dvp or sgr_dist)
!> CLEANUP: This routine replicates parallel_dissipation and should be merged
!------------------------------------------------------------------------------
subroutine diffus(E_in,kdiff,dum,direction,iorder)

  use geom,         only : sgr_dist
  use velocitygrid, only : dvp
  use structures,   only : matrix_element
  use matdat,       only : add_element, connect_parallel
  use general,      only : gkw_abort 
  use mode,         only : parallel_phase_shift

  type(matrix_element), intent(in) :: E_in
  real, intent(in) :: kdiff,dum
  character (len=*), intent(in) :: direction
  integer, intent(in) :: iorder

  type(matrix_element) :: E
  real :: d,dx
  real, dimension(5) :: st
  integer :: ierr, jj
  logical :: ingrid 

  dx = 1. ; d = 0.

  E%ideriv = 4
  E%itloc = E_in%itloc
  E%itype = E_in%itype

  select case(direction)
    case('vpar')
      dx = dvp
      E%term='parallel velocity dissipation (ARA)'
    case('s')
      dx = sgr_dist
      E%term='parallel dissipation (ARA)'
    case default
      call gkw_abort('diffus: bad case of direction')
  end select
  
  select case (iorder)
    case (4)
      st = (/ 1., -4., 6., -4., 1./)
      d=-kdiff*abs(dum)/(12.*dx)
    case (2)
      st = (/ 0., -1., 2., -1., 0./)
      d=-kdiff*abs(dum)/(dx)
    case default
      call gkw_abort('diffus: bad case of scheme order')
  end select
       
  do jj=1,5
 
    E     = E_in  
    E%val = d*st(jj)

    select case(direction)
      case('vpar')
        E%iloc = E%i
        E%kloc = E%k + (jj - 3)
      case('s')
        E%iloc = E%i + (jj - 3)
        E%kloc = E%k
        E%val  = E%val * parallel_phase_shift(E%imod,E%ix,E%i,E%iloc)
    end select
    
    call connect_parallel(E,ingrid)
    if (ingrid) call add_element(E,ierr)
    
  end do

end subroutine diffus

!-----------------------------------------------------------------------------
!> returns the terms appearing in front of the the radial derivative 
!> of the maxwellian for a given local species input 
!-----------------------------------------------------------------------------
function dmaxwel(ix,i,j,k,is,imod)

 use velocitygrid, only : vpgr, mugr
 use components, only : tgrid, tmp, fp, vp, mas, tp, tm_drive, signz, ecurx
 use components, only : types
 use geom,       only : jfun, lfun, bt_frac, bn, rfun
 use rotation,   only : vcor, cfen, ts_uprim, cf_trap, cf_upsrc, toroidal_shear
  
 real :: dmaxwel
 integer, intent(in) :: ix, i,j,k,is,imod !< local species number

 real :: ET, cf_dum, up_prim !< dummies
 
  ! term in front of the temperature gradient - To avoid future confusion,
  ! the factor of 2 in front of the mugr term comes from the normalisations of
  ! the magnetic moment
  if (types(is) .eq. 'EP') then
    ET = (2.E0*mugr(j)*bn(ix,i) + cfen(i,is) )                &
        &/ (tmp(ix,is) / tgrid(is) ) - 1.5E0
  else 
    ET = (vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i) + cfen(i,is) ) &
        & / (tmp(ix,is) / tgrid(is) ) - 1.5E0 
  endif 
 
 !calculate u_parallel^prime
  if (toroidal_shear.eq.'use_shear_rate') then
    up_prim=bt_frac(ix,i)*ts_uprim * rfun(ix,i)
  else !vp is uprim
    up_prim=bt_frac(ix,i)*vp(ix,is) * rfun(ix,i)
  end if

  if(tm_drive)then
    !Exact form of drive goes here -> Only present in electrons
    if(signz(is).lt.0)then
      !Negative sign here 
      up_prim = up_prim + ecurx(ix)
      !write(*,*) is,ix,ecurx(ix)
    endif
  endif

  ! correction for cf_trap (radial derivative in density gradient)
  ! EXPERIMENTAL
  if (cf_trap) then
    cf_dum=mas(is)*vcor*vcor*lfun
    if (cf_upsrc)then
      if ( ( ( all(apply_on_imod == 0) .or. any(apply_on_imod == imod) ) .and. lneorotsource ) .or. & 
      & ( ( .not.all(apply_on_imod == 0) .and. .not.any(apply_on_imod == imod) ) .and. lneo_equil_switch_default ) ) then
        cf_dum=cf_dum
      else
        cf_dum=cf_dum+2.*vp(ix,is)*mas(is)*vcor*jfun(ix,i)
      endif
    end if
    cf_dum=cf_dum/tmp(ix,is)
  else 
    cf_dum=0.0
  end if
 
  ! return the final value 
  if (types(is) .eq. 'EP') then
    dmaxwel=fp(ix,is)+ET*tp(ix,is)
  else 
    dmaxwel=fp(ix,is)+ET*tp(ix,is)+2.E0*sqrt(mas(is)*tgrid(is))*up_prim*vpgr(i,j,k,is) &
         & /tmp(ix,is) + cf_dum
  endif

end function

!--------------------------------------------------------------------------
!> Wrapper function for the bessel functions. In the case of the spectral
!> method this function returns the Bessel function and the product with 
!> the potential gives the gyro-averaged potential. In the case of 
!> spectral_radius = .false. the function returns 1 since all elements 
!> refer directly to the gyro-averaged potential, not the potential 
!--------------------------------------------------------------------------
function wrap_beslj0_gkw(imod,ix,i,j,is) 

  use control,   only : spectral_radius
  use functions, only : besselj0_gkw 

  integer, intent(in) :: imod, ix, i, j, is 
  real                :: wrap_beslj0_gkw 

  if (spectral_radius) then 
    wrap_beslj0_gkw = besselj0_gkw(imod,ix,i,j,is) 
  else 
    wrap_beslj0_gkw = 1.0E0 
  endif 

end function wrap_beslj0_gkw 

!--------------------------------------------------------------------------
!> Similar to the wrapper function wrap_beslj0_gkw. For spectral_radius = 
!> .false. this function returns 1., since the elements refer directly to 
!> the gyro-average of Bpar
!--------------------------------------------------------------------------
function wrap_mod_besj1_gkw(imod,ix,i,j,is) 

  use control,   only : spectral_radius 
  use functions, only : mod_besselj1_gkw 

  integer, intent(in) :: imod, ix, i, j, is
  real                :: wrap_mod_besj1_gkw 

  if (spectral_radius) then 
    wrap_mod_besj1_gkw = mod_besselj1_gkw(imod,ix,i,j,is) 
  else 
    wrap_mod_besj1_gkw = 1.0E0 
  endif 
 
end function wrap_mod_besj1_gkw 

!-------------------------------------------------------------------------
!> Calculates the parallel (s) derivative of the neo background
!> distribution function
!> At the moment assumes the local model. -> Possible expanded in future
!-------------------------------------------------------------------------
function neo_dFds(ix,i,j,k,is)

  use neoequil, only : neof, mid_r_point
  use geom, only : ffun, sgr_dist
  use dist,   only : fmaxwl, stencil_side, stencil_side_zf
  use global, only : gkw_a_equal_b_accuracy, id_s

  integer :: ist, m, ii
  integer, intent(in) :: ix, i, j, k, is
  real    :: dum, neo_dFds
  real, allocatable :: w(:)

  allocate(w(1 + 4*max(stencil_side(id_s), stencil_side_zf(id_s))))

  dum = 0.E0
  ! select the scheme 
  ist = 0!Always in the grid as it is periodic !pos_par_grid(1,ix,i,k)
  call differential_scheme(ist,1,1,w)
  do m = 1, size(w)
    if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
      !Only ix=2 in neof as it is local 
      ii = i + m - ((size(w)+1)/2)
      dum = dum + w(m)*ffun(ix,i)* &
        & fmaxwl(ix,i,j,k,is)*neof(mid_r_point,is,ii,j,k)/sgr_dist
    endif
  enddo
  neo_dFds = dum

  deallocate(w)
  
end function neo_dFds

!-------------------------------------------------------------------------
!> Calculates the parallel (s) derivative of the neo background
!> distribution function
!> At the moment assumes the local model. -> Possible expanded in future
!-------------------------------------------------------------------------
function dphineods(imod,ix,i)

  use neoequil, only : neophi, mid_r_point
  use geom, only : sgr_dist
  use dist, only : stencil_side, stencil_side_zf
  use global, only : gkw_a_equal_b_accuracy, id_s

  integer :: ist, m, ii
  integer, intent(in) :: imod, i, ix
  real    :: dum, dphineods
  real, allocatable :: w(:)
  allocate(w(1 + 4*max(stencil_side(id_s), stencil_side_zf(id_s))))

  dphineods = 0.E0
  
  if (all(apply_on_imod == 0)) then
    if(.not. lneo_equil) return
  else
    if (any(apply_on_imod == imod)) then
      if(.not. lneo_equil) return
    else
      if(.not.lneo_equil_switch_default) return
    endif
  endif

  dum = 0.E0
  ! select the scheme 
  ist = 0!Always in the grid as it is periodic !pos_par_grid(1,ix,i,k)
  call differential_scheme(ist,1,1,w)
  do m = 1, size(w)
    if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
      !Only ix=2 in neof as it is local 
      ii = i + m - ((size(w)+1)/2)
      dum = dum + w(m)* &
        & neophi(mid_r_point,ii)/sgr_dist
    endif
  enddo
  dphineods = dum

  deallocate(w)
  
end function dphineods

!--------------------------------------------------------------------
!> Calculates the parallel velocity (vpar) derivative of the
!> neo background distribution function 
!--------------------------------------------------------------------
function neo_dFdvpar(ix,i,j,k,is)

  use neoequil, only : neof, mid_r_point
  use velocitygrid, only : dvp
  use dist,   only : fmaxwl, stencil_side, stencil_side_zf
  use global, only : gkw_a_equal_b_accuracy, id_s

  integer :: ist, m, ii
  integer, intent(in) :: ix, i, j, k, is
  real    :: dum, neo_dFdvpar
  real, allocatable :: w(:)

  allocate(w(1 + 4*max(stencil_side(id_s), stencil_side_zf(id_s))))

  dum = 0.E0

  ist = 0
  call differential_scheme(ist,1,1,w)
  do m = 1, size(w)
    if (.not. gkw_a_equal_b_accuracy(w(m), 0.0)) then
      ii = k + m - ((size(w)+1)/2)
      dum = dum + w(m)*fmaxwl(ix,i,j,k,is)*neof(mid_r_point,is,i,j,ii)/dvp
    endif
  enddo
  neo_dFdvpar = dum

  deallocate(w)

end function neo_dFdvpar

end module linear_terms
