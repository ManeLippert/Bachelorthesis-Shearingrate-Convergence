!----------------------------------------------------------------------
!> Contains calculation of all the geometry quantities (tensors)
!! for six models:
!! 1) s-alpha
!! 2) circular
!! 3) General geometry read from CHEASE output file.
!! 4) slab
!! 5) Miller
!! 6) Fourier
!! This module is used in the initialisation phase.
!<---------------------------------------------------------------------

module geom

  use global,       only : lenswitch

  implicit none

  private

  public :: geom_read_nml, geom_write_nml, geom_bcast_nml
  public :: geom_check_params, parallelize_geom 
  public :: geom_output, q_profile, qprof_coef

  public :: bn, bn_G, dfun, efun, eps, ffun, geom_allocate, bt_frac, q, shatx
  public :: geom_init_grids, gfun, hfun, ints, kthnorm, kxnorm, metric, qx, sgr
  public :: sgr_dist, shat, bmin, bmax, metric_G, jacobian_G, ifun, jfun, kfun 
  public :: Rfun, beta_eq, betaprime_eq, Rref, Bref, pol_angle, lfun, jfunh 
  public :: jfunl, dBdpsi, alphak, alphak_xbnd, shift_end_grid, dxgr, xgr
  public :: betaprime_miller, tmp_miller, de_miller, vp_miller_i, fp_miller
  public :: tp_miller, mas_miller_i, vcor_miller, beta_miller, dpfdpsi, R0 
  public :: gradp_type, dpdpsi_rot, dpds_rot, curv_effect, iterate_fsa 
  public :: interpquad
  public :: isg_lfs
  
  !Ideally the ones below should not need to be public.
  public :: geom_type, signB, signJ, geom_parallelized, R0_loc

  interface geom_write_nml
    module procedure geom_read_nml
  end interface

  !> local radial coordinate (read from file) 
  !> Note that 'psi' rather than 'eps' is used to mean the radial coordinate in 
  !> variable names (as in the GKW manual).
  !> In short, psi=eps, except briefly at initialisation when eps_type=2.
  !> Note also that to avoid confusion, 'pf' is used to mean the poloidal flux 
  real, save :: eps

  !> radial coordinate grid, this is a constant value for a flux tube
  !> run, and equidistant for a global run.
  real, allocatable, save :: xgr(:)
  
  !> Local safety factor. This value is read from input. It is used in all 
  !> routines that require the flux tube solution.  
  real, save :: q

  !> q profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: qx(:)

  !> Local value of the magnetic shear (read from input)
  real, save :: shat

  !> Magnetic shear profile. For global runs 
  real, allocatable, save :: shatx(:) 

  !> Switch for the geometry (private)
  character (len = lenswitch), save :: geom_type  

  !> Length of chararacter string of filename 
  integer, parameter :: lenfile = 180 

  !> File to read the equilibrium related quantities (private)
  character (len = lenfile), save :: eqfile

  !> Sign of B.grad_phi, the toroidal component of the magnetic field
  integer, save :: signB

  !> Sign of j.grad_phi, the toroidal component of the plasma current
  integer, save :: signJ

  !> Radial coordinate used to specify the chosen FS
  !> 1=eps, 2=rho_pf (sqrt of poloidal flux)
  integer, save :: eps_type

  !> The plasma beta corresponding to the MHD equilibrium
  real, save :: beta_eq 

  !> and betaprime
  real, save :: betaprime_eq 

  !> The reference major radius 
  real, save :: Rref 

  !> The reference magnetic field
  real, save :: Bref 

  !> Grid distance along the field line 
  real, save :: sgr_dist 

  !> Grid distance in the radial direction 
  real, save :: dxgr 

  !> the normfactor for calculating the k_zeta values in the code 
  !> initial value given to catch possible errors
  real, save :: kthnorm=1.23e4
  
  !> normfactor for calculating LX values in real space at LFS
  !> (and calculating k_psi in the code when kr_switch='kr') 
  real, save :: kxnorm=1.0

  !> The minimum and maximum magnetic field on the flux surface 
  real, allocatable, save :: bmin(:), bmax(:)

  !> the normalized magnetic field strength : bn_G (GLOBAL)
  !> bn is normalised to Bref at the magnetic axis.
  real, allocatable, save :: bn_G(:,:)

  !> the normalized magnetic field strength : bn (LOCAL)
  real, allocatable, save :: bn(:,:)

  !> The covariant component of the magnetic field bups(n_x_grid) 
  real, allocatable, save :: bups(:)

  !> The coordinate along the field line : sgr(0:n_s_grid+1)
  real, allocatable, save :: sgr(:)

  !> and associated poloidal angle: pol_angle(n_x_grid,n_s_grid)
  real, allocatable, save :: pol_angle(:,:)

  !> this array is used for the flux surface average.
  real, allocatable, save :: ints(:)

  !> the curvature operator  dfun(n_x_grid,n_s_grid,3) 
  real, allocatable, save :: dfun(:,:,:)

  !> the ExB operator efun(n_x_grid,n_s_grid,3,3) 
  real, allocatable, save :: efun(:,:,:,:)

  !> the factor in front of the parallel derivative ffun(n_x_grid,n_s_grid)
  real, allocatable, save :: ffun(:,:) 

  !> high and low field side values of jfun (exact, not nearest point)
  real, save :: jfunh, jfunl

  !> The trapping operator gfun(n_x_grid,n_s_grid)
  real, allocatable, save :: gfun(:,:)

  !> The tensor that determines the Coriolis drift: hfun(n_x_grid,n_s_grid,3) 
  real, allocatable, save :: hfun(:,:,:)

  !> The array connected to the centrifugal drift: ifun(n_x_grid,n_s_grid,3) 
  real, allocatable, save :: ifun(:,:,:)

  !> centrifugal trapping
  !> The quantity R^2-(R_0)^2.  Peeters PoP 09
  real, allocatable, save :: jfun(:,:)

  !> centrifugal trapping
  !> 2 R*dR/dpsi at constant s
  real, allocatable, save :: kfun(:,:)

  !> Local major radius R for angular momentum flux calculation
  !> Note: R=1 in s-alpha Rfun(n_x_grid,n_s_grid)
  real, allocatable, save :: Rfun(:,:)

  !> |Bt| / B bt_frac(n_x_grid,n_s_grid)
  real, allocatable, save :: bt_frac(:,:)

  !> Normalised radius at which input density and gradient are defined
  !> Only used with centrifugal terms.
  !> Value set according to input R0_loc
  real, save :: R0 

  !> d(R0^2)/dpsi (at constant theta)
  real, save :: lfun

  !> Switch to specify location of R0 definition
  character (len = lenswitch), save :: R0_loc 

  !> the metric (GLOBAL in s and x)
  real, allocatable, save :: metric_G(:,:,:,:)

  !> the jacobian (GLOBAL in x)
  real, allocatable, save :: jacobian_G(:)

  !> the metric (LOCAL in s and x)
  real, allocatable, save :: metric(:,:,:,:)

  !> Flag to set after geom has been parallelized
  logical, save :: geom_parallelized=.false.

  !> if rotation needs to iterate on poloidal asymmetries
  logical, save :: iterate_fsa=.false.

  !> The radial derivative of the magnetic field (at constant s)
  real, allocatable , save :: dBdpsi(:,:)
 
  !> The derivative of the magnetic field towards the parallel coordinate
  real, allocatable, save :: dBds(:,:) 

  !> The radial derivative of the Z coordinate of the flux surface (at const s) 
  real, allocatable , save :: dZdpsi(:,:)

  !> The parallel derivative of the Z coordinate of the flux surface 
  real, allocatable, save :: dZds(:,:) 

  !> The z-coordinate 
  real, allocatable, save :: Z_FS(:,:)

  !> the radial derivative of the major radius (at constant s)
  real, allocatable , save :: dRdpsi(:,:)

  !> the parallel derivative of the major radius   
  real, allocatable, save :: dRds(:,:) 

  !> first derivative of poloidal flux towards the radial coordinate (at const s)  
  real, allocatable, save :: dpfdpsi(:)

  !> Shift of the binormal coordinate (GLOBAL in s and x)
  real, allocatable, save :: alphakp(:,:) 

  !> Shift of the binormal coordinate (LOCAL in s and x)
  real, allocatable, save :: alphak(:,:) 

  !> Shift of the binormal coordinate at x boundary
  real, allocatable, save :: alphak_xbnd(:) 

  !> Shift in phase of toroidal mode when reconnection along s, for nonspectral
  real, allocatable, save :: shift_end_grid(:)

  !> Switch for various q-profile functions 
  character (len = lenswitch), save :: prof_type  
  
  !> Coefficients for the analytic q-profile as a function of radius 
  real, save :: qprof_coef(5) 

  !> Below parameters used in miller geometry
  !> R = Rmil + r*cos(theta + arcsin(delta)*sin(theta))
  !> Z = Zmil + kappa*r*sin(theta + zeta*sin(2*theta)) 
  !> Elongation
  real, save :: kappa

  !> Elongation profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: kappax(:)
  
  !> Triangularity
  real, save :: delta

  !> Triangularity profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: deltax(:)
  
  !> Squareness
  real, save :: square

  !> Squareness profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: squarex(:)
  
  !> Radial derivative of kappa
  real, save :: skappa

  !> radial elongation derivative profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: skappax(:)
  
  !> Radial derivative of delta (definition used in Miller et al. PoP 5 973 (1998)
  real, save :: sdelta

  !> radial Triangularity derivative profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: sdeltax(:)
  
  !> Radial derivative of zeta
  real, save :: ssquare
 
  !> radial Squareness derivative profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: ssquarex(:)

  !> Elevation
  real, save :: Zmil

  !> radial elevation profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: zmilx(:)

  !> Radial derivative of the effective major radius
  real, save :: dRmil
  
  !> radial major raidus derivative profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: dRmilx(:)

  !> Radial derivative of elevation
  real, save :: dZmil

  !> radial elevation derivative profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: dZmilx(:)
  
  !> Pressure gradient (with respect to poloidal flux, SI units)
  real, save :: gradp

  !> radial pressure gradient profile. This is used in all routines that are called in the global
  !> version of the code. Note that this involves also routines that are 
  !> called in the local flux tube version 
  real, allocatable, save :: gradpx(:)
  
  !> Switch to determine which definition is taken for the input gradp
  character (len = lenswitch), save :: gradp_type
  
  !> Betaprime if used by miller (value provided by the components module)
  real, allocatable, save :: betaprime_miller(:)
  
  !> Normalised temperature for e- and main ions if used by miller 
  !> (value provided by the components module)
  real, allocatable, save :: tmp_miller(:,:)
  
  !> Normalised density for e- and main ions if used by miller 
  !> (value provided by the components module)
  real, allocatable, save :: de_miller(:,:)
  
  !> Toroidal rotation gradient for main ion species if used by miller 
  !> (value provided by the components module)
  real, allocatable, save :: vp_miller_i(:)
  
  !> Density gradient for e- and main ions if used by miller 
  !> (value provided by the components module)
  real, allocatable, save :: fp_miller(:,:)
  
  !> Temperature gradient for e- and main ions if used by miller 
  !> (value provided by the components module)
  real, allocatable, save :: tp_miller(:,:)  
  
  !> Mass for the main ions if used by miller 
  !> (value provided by the components module)
  real, save :: mas_miller_i
  
  !> Beta if used by miller (value provided by the components module)
  real, allocatable, save :: beta_miller(:)
  
  !> Toroidal rotation of the plasma if used by miller
  !> (value provided by the components module)
  real, save :: vcor_miller = 1.23e4
  
  !> Radial derivative of the pressure when toroidal rotation is taken 
  !> into account in geom miller
  real, allocatable, save :: dpdpsi_rot(:,:)
  
  !> Parallel derivative of the pressure when toroidal rotation is taken 
  !> into account in geom miller
  real, allocatable, save :: dpds_rot(:,:)

  !> Effect of toroidal rotation on the curvature drift when geom Miller is used
  logical, save :: curv_effect = .true.

  !> Switch to determine which beta is used for gradp_type = 'rota_miller'
  character (len = lenswitch), save :: beta_rota_miller_type

  !> Used in rota_miller if beta_rota_miller_type = 'geom'
  real, save :: beta_rota_miller

  !> Below parameters used in fourier geometry
  !> R = Rf + a*cos(theta) and Z = Zf + a*sin(theta)
  !> a = sum_n c*cos(n*theta) + s*sin(n*theta) 
  !> Number of Fourier moments
  integer, save :: N_shape
  !> Max number of Fourier moments
  integer, parameter :: N_shape_max=50
  !> cosine components of the minor radius
  real, save :: c(N_shape_max)
  !> sine components of the minor radius
  real, save :: s(N_shape_max)
  !> radial derivative of cosine components 
  real, save :: c_prime(N_shape_max)
  !> radial derivative of the sine components
  real, save :: s_prime(N_shape_max)

  ! IF you add a new geometry tensor, it must appear in:
  ! 1) calc_geom_tensors
  ! 2) chease calculation (and bcast)
  ! 3) paralleize_geom (_radial)
  ! 4) geom_output

  !> the global index (with respect to the parallel coordinate) of the point
  !> at the low field side (LFS). This is right in the middle of the
  !> grid at the moment.
  integer, save :: isg_lfs

contains

!------------------------------------------------------------------------------
!> read (or write) geom namelist
!------------------------------------------------------------------------------

subroutine geom_read_nml(ifile,io_stat,lwrite)
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  integer, intent(in)  :: ifile
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  namelist /geom/ shat, q, eps, eps_type, geom_type, eqfile, signB,       &
  &  signJ, R0_loc, prof_type, qprof_coef, kappa, skappa, delta, sdelta,  &
  &  square, ssquare, Zmil, dRmil, dZmil, gradp_type, gradp, curv_effect, &
  &  beta_rota_miller_type, beta_rota_miller, &
  &  N_shape, c, s, c_prime, s_prime

  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then

      ! Set the standard values  
      q    = 1
      shat = 1.23e4
      eps  = 1
      eps_type = 1
      geom_type = 's-alpha'
      eqfile = "not_needed_for_s-alpha_option"
      signB = 1
      signJ = 1
      R0_loc = 'LFS'
      prof_type = 'parabolic' 
      qprof_coef = 0.E0 
      qprof_coef(1) = 0.8854
      qprof_coef(2) = 19.656
      kappa = 1.23e4
      skappa = 1.23e4
      delta = 1.23e4
      sdelta = 1.23e4
      square = 1.23e4
      ssquare = 1.23e4
      Zmil = 1.23e4
      dRmil = 1.23e4
      dZmil = 1.23e4
      gradp_type = 'alpha_mhd'
      gradp = 1.23e4
      curv_effect = .true.
      beta_rota_miller_type = 'spc'
      beta_rota_miller = 1.23e4
      N_shape = 1
      c = 1.23e4
      s= 1.23e4
      c_prime = 1.23e4
      s_prime = 1.23e4

      ! read nml
      read(ifile,NML=geom,IOSTAT=io_stat)

    else
      ! do nothing
    end if
  else
    ! write nml
    if(root_processor) write(ifile,NML=geom)

    ! write metadata
    call write_run_parameter('geom', 'shat', shat)
    call write_run_parameter('geom', 'q', q)
    call write_run_parameter('geom', 'eps', eps)
    call write_run_parameter('geom', 'eps_type', eps_type)
    call write_run_parameter('geom', 'geom_type', geom_type)
    call write_run_parameter('geom', 'eqfile', eqfile)
    call write_run_parameter('geom', 'signB', signB)
    call write_run_parameter('geom', 'signJ', signJ)
    call write_run_parameter('geom', 'R0_loc', R0_loc)
    call write_run_parameter('geom', 'prof_type', prof_type)
    call write_run_parameter('geom', 'qprof_coef', qprof_coef)
    call write_run_parameter('geom', 'kappa', kappa)
    call write_run_parameter('geom', 'skappa', skappa)
    call write_run_parameter('geom', 'delta', delta)
    call write_run_parameter('geom', 'sdelta', sdelta)
    call write_run_parameter('geom', 'square', square)
    call write_run_parameter('geom', 'ssquare', ssquare)
    call write_run_parameter('geom', 'Zmil', Zmil)
    call write_run_parameter('geom', 'dRmil', dRmil)
    call write_run_parameter('geom', 'dZmil', dZmil)
    call write_run_parameter('geom', 'gradp_type', gradp_type)
    call write_run_parameter('geom', 'gradp', gradp)
    call write_run_parameter('geom', 'curv_effect', curv_effect)
    call write_run_parameter('geom', 'beta_rota_miller_type', beta_rota_miller_type)
    call write_run_parameter('geom', 'beta_rota_miller', beta_rota_miller)
    call write_run_parameter('geom', 'N_shape', N_shape)
    call write_run_parameter('geom', 'c', c)
    call write_run_parameter('geom', 's', s)
    call write_run_parameter('geom', 'c_prime', c_prime)
    call write_run_parameter('geom', 's_prime', s_prime)
   
  end if
  
end subroutine geom_read_nml
               

!------------------------------------------------------------------------------
!> broadcast geom input parameters
!------------------------------------------------------------------------------

subroutine geom_bcast_nml

  use mpiinterface, only : mpibcast

  call mpibcast(q,                             1)
  call mpibcast(shat,                          1)
  call mpibcast(eps,                           1)
  call mpibcast(eps_type,                      1)
  call mpibcast(geom_type,             lenswitch)
  call mpibcast(eqfile,                  lenfile)
  call mpibcast(signB,                         1)
  call mpibcast(signJ,                         1)
  call mpibcast(R0_loc,                lenswitch)
  call mpibcast(prof_type,             lenswitch)
  call mpibcast(qprof_coef,                    5) 
  call mpibcast(kappa,                         1)
  call mpibcast(skappa,                        1)
  call mpibcast(delta,                         1)
  call mpibcast(sdelta,                        1)
  call mpibcast(square,                        1)
  call mpibcast(ssquare,                       1)
  call mpibcast(Zmil,                          1)
  call mpibcast(dRmil,                         1)
  call mpibcast(dZmil,                         1)
  call mpibcast(gradp,                         1)
  call mpibcast(gradp_type,            lenswitch)
  call mpibcast(curv_effect,                   1)
  call mpibcast(beta_rota_miller_type, lenswitch)
  call mpibcast(beta_rota_miller,              1)
  call mpibcast(N_shape,                       1)
  call mpibcast(c,                   N_shape_max)
  call mpibcast(s,                   N_shape_max)
  call mpibcast(c_prime,             N_shape_max)
  call mpibcast(s_prime,             N_shape_max)

end subroutine geom_bcast_nml


!------------------------------------------------------------------------------
!> check geom input parameters
!> in the case of chease this routine may be called twice
!------------------------------------------------------------------------------
subroutine geom_check_params(icall)

  use grid,    only : n_s_grid, nperiod
  use general, only : gkw_abort, gkw_warn, gkw_exit
  use control, only : vp_trap, flux_tube
  use global,  only : r_tiny


  integer, intent(in) :: icall
  integer             :: i
  logical             :: eqfile_exists

  if (abs(eps-1.23e4) < 1e-4) then
    call gkw_exit('geom_check: You have not specified eps in the input')
  endif

  if (eps < 0.) then
    call gkw_exit('geom_check: I do not understand negative aspect ratio')
  end if

  if(q < 0.) then
    call gkw_exit('geom_check: To run negative q, please reverse signJ or signB')
  end if

  ! note to self, the dbdpsi_LFS in also inexact for even n_s_grid
  if (icall == 1) then
    if (R0_loc=='FSA') then
      R0_loc='LFS'
      call gkw_warn('input densities and gradients interpreted as FSA')
      call gkw_warn('requires iterative solution of poloidal assymmetries')
      call gkw_warn('R0_loc = LFS used internally')
      iterate_fsa = .true.
    else 
      iterate_fsa = .false.
    end if
  end if

  if (mod(n_s_grid,2*nperiod-1).ne.0) then
    call gkw_exit('N_s_grid/(2*NPERIOD-1) should be an integer')
    stop 1
  end if

  if (geom_type == 'miller') then

      if (curv_effect .eqv. .false.) then
        call gkw_warn('Effect of toroidal rotation not taken into account in  &
        &  the curvature drift: should be used only for testing')
      end if

      ! only test the parameters if they are not read from file 
      if (prof_type /= 'file') then 
      
        if (abs(kappa-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified kappa in the input')
        endif

        if (kappa < 0.) then
          call gkw_exit('geom_check: I do not understand negative elongation')
        end if

        if (abs(delta-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified delta in the input')
        endif
      
        if (abs(delta) > 1.) then
          call gkw_exit('geom_check: delta must be between -1 and 1')
        end if

        if (abs(square-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified square in the input')
        endif

        if (abs(skappa-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified skappa in the input')
        endif
      
        if (skappa < -1) then
          call gkw_exit('geom_check: No miller eqm. for skappa < -1 ')
        endif

        if (abs(sdelta-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified sdelta in the input')
        endif

        if (abs(ssquare-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified ssquare in the input')
        endif

        if (abs(Zmil-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified Zmil in the input')
        endif

        if (abs(dRmil-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified dRmil in the input')
        endif
      
        if (abs(dRmil - 1.0) < r_tiny) then
          ! that means dRmil equal to 1.
          call gkw_warn('Bn diverges at HFS for dRmil = 1 in Miller geom')
        end if

        if (abs(dZmil-1.23e4) < 1e-4) then
          call gkw_exit('geom_check: You have not specified dZmil in the input')
        endif

        if (abs(gradp-1.23e4) < 1e-4 .and. .not. gradp_type=='beta_prime'      &
        &  .and. .not. gradp_type=='rota_miller') then
          call gkw_exit('geom_check: You have not specified gradp in the input')
        endif
  
      endif 
      
      !bmax, bmin, jfunh, jfunl are inexact
      if (vp_trap == 1 .and. R0_loc == 'axis') then
        call gkw_warn('Miller geom has inexact trapping grid for R0_loc=axis')
      endif 

      if (abs(beta_rota_miller-1.23e4) < 1e-4 .and. beta_rota_miller_type == 'geom') then
        call gkw_exit('geom_check: You have not specified beta_rota_miller in the input')
      endif

      if (.not. flux_tube) then
        if (gradp_type == 'rota_miller') then
          call gkw_exit('geom_check: gradp=rota_miller not available for global runs')
        endif
      endif
  end if
  
  if (geom_type == 'fourier') then

      if (abs(gradp-1.23e4) < 1e-4 .and. .not. gradp_type=='beta_prime') then
        call gkw_exit('geom_check: You have not specified gradp in the input')
      endif

      if (gradp_type == 'rota_miller') then
        call gkw_exit('geom_check: gradp=rota_miller not available for fourier parametrisation')
      endif

      if (prof_type=='file') then
        call gkw_exit('geom_check: prof_type=file not yet implemented for Fourier parametrisation')
      endif

      if (N_shape > N_shape_max) then
        call gkw_exit('geom_check: max number of shape moments is 50')
      endif

    do i=1,N_shape
      if (abs(c(i)-1.23e4) < 1e-4) then
        call gkw_exit('geom_check: You have not specified all c values in the input')
      endif

      if (abs(s(i)-1.23e4) < 1e-4) then
        call gkw_exit('geom_check: You have not specified all s values in the input')
      endif

      if (abs(c_prime(i)-1.23e4) < 1e-4) then
        call gkw_exit('geom_check: You have not specified all c_prime values in the input')
      endif

      if (abs(s_prime(i)-1.23e4) < 1e-4) then
        call gkw_exit('geom_check: You have not specified all s_prime values in the input')
      endif
    enddo

    if (N_shape < N_shape_max) then    ! set the unused Fourier components to zero
      do i=N_shape+1,N_shape_max
        c(i)=0.
        s(i)=0.
        c_prime(i)=0.
        s_prime(i)=0.
      enddo
    endif

  end if

  if (geom_type /= 'chease') then
    if (abs(q-1.23e4) < 1e-4) then
      call gkw_exit('geom_check: You have not specified q in the input')
    endif
    if (abs(shat-1.23e4) < 1e-4) then
      call gkw_exit('geom_check: You have not specified shat in the input')
    endif
    eqfile = 'not_needed_for_s-alpha_option'
  else if (geom_type == 'chease') then
    if (eqfile == 'not_needed_for_s-alpha_option') then
      call gkw_exit('geom_check: You have not specified eqfile in the input')
    endif
    inquire(file=eqfile,EXIST=eqfile_exists)
    if (.not. eqfile_exists) call gkw_abort('geom: not found: '//eqfile)
    !Namelist input values of q and shat are not used for chease case
    !The values are checked again in kgrid and krbal to make sure they have
    !been read from the chease input and set correctly in geom_init_grids
    if(icall==1) then
      if(q<1.22e4) then
        call gkw_warn('Namelist input q not used for geom type chease')
        q = 1.23e4
      endif
      if(shat<1.22e4) then
        call gkw_warn('Namelist input shat not used for geom type chease')
        shat = 1.23e4
      endif
    endif !call
    
    if (mod(n_s_grid,2*nperiod-1).ne.0) then
      call gkw_exit( &
         & 'N_s_grid/(2*NPERIOD-1) has to be an integer for geom_type=chease')
    end if
        
  endif !geom type

  ! make sure that abs(signB)=1 and abs(signJ)=1
  if (signB >= 0) then
    signB = 1
  else
    signB = -1
  endif
  if (signJ >= 0) then
    signJ = 1
  else
    signJ = -1
  endif

end subroutine geom_check_params


!------------------------------------------------------------------------------
!> allocate the arrays of geom
!------------------------------------------------------------------------------

subroutine geom_allocate

  use grid,    only : n_s_grid, n_x_grid, ns, parallel_s, nx
  use general, only : gkw_abort
  
  ! local parameters 
  integer ierr

  ! intialize the error parameter
  ierr=0

  ! allocate the radial grid array 
  allocate(xgr(n_x_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate xgr in geom')
  
  ! allocate the safety factor profile 
  allocate(qx(n_x_grid), stat=ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate qx in geom')

  ! allocate the magnetic shear profile 
  allocate(shatx(n_x_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate shatx in geom')

  ! allocate the geom/component communication array betaprime_miller
  allocate(betaprime_miller(n_x_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate betaprime_miller in geom')
  betaprime_miller(:) = 1E8

  ! allocate the geom/component communication array beta_miller
  allocate(beta_miller(n_x_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate beta_miller in geom') 
  beta_miller(:)=  1.23E4

  ! For the miller geometry allocate the radial profiles 
  if (geom_type == 'miller') then 
  
    ! allocate the elongation profile 
    allocate(kappax(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate kappax in geom')

    ! allocate the Triangularity profile 
    allocate(deltax(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate deltax in geom')

    ! allocate the Squareness profile 
    allocate(squarex(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate squarex in geom')

    ! allocate the elongation derivative profile 
    allocate(skappax(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate skappax in geom')
  
    ! allocate the Triangularity derivative profile 
    allocate(sdeltax(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate sdeltax in geom')

    ! allocate the Squareness derivative profile 
    allocate(ssquarex(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate ssquarex in geom')
  
    ! allocate the elevation profile 
    allocate(Zmilx(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate Zmilx in geom')

    ! allocate the major radius derivative profile 
    allocate(dRmilx(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate dRmilx in geom')

    ! allocate the elevation derivative profile 
    allocate(dZmilx(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate dZmilx in geom')

    ! allocate the pressure gradient profile 
    allocate(gradpx(n_x_grid), stat=ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate gradpx in geom')

  endif 
  
  ! allocate the magnetic field array
  allocate(bn_G(n_x_grid,1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate bn_G in geom')

  ! allocate the magnetic field array
  ! 2 extra points, parallel or otherwise 
  allocate(bn(n_x_grid,1-3:ns+3),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate bn in geom')

  ! allocate the array for the minimum of the magn field on the
  ! respective flux surface
  allocate(bmin(n_x_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate bmin in geom')

  ! allocate the array for the maximum of the magn field on the
  ! respective flux surface
  allocate(bmax(n_x_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate bmax in geom')

  ! allocate the field line length array
  allocate(sgr(0:n_s_grid+1),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate sgr in geom')

  ! allocate the poloidal angle array
  allocate(pol_angle(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate pol_angle in geom')

  ! allocate the array for integration along the 
  ! field line 
  allocate(ints(n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate ints in geom')

  ! allocate the array that contains the function 
  ! connected witht the curvature operator 
  allocate(dfun(n_x_grid,n_s_grid,3), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate dfun in geom')

  ! allocate the array that contains the function 
  ! connected with the ExB velocity 
  allocate(efun(n_x_grid,n_s_grid,3,3), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate efun in geom')

  ! allocate the array that contains the function in 
  ! front of the parallel derivative 
  allocate(ffun(n_x_grid,n_s_grid), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate ffun in geom')

  ! allocate the array that contains the function 
  ! connected with the trapping terms
  allocate(gfun(n_x_grid,n_s_grid), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate gfun in geom')

  ! allocate the array that contains the function 
  ! connected with the coriolis terms 
  allocate(hfun(n_x_grid,n_s_grid,3), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate hfun in geom')

  ! allocate the array that contains the function 
  ! connected with the centrifugal terms 
  allocate(ifun(n_x_grid,n_s_grid,3), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate ifun in geom')

  ! allocate the array that contains the function 
  ! R^2-R_0^2
  allocate(jfun(n_x_grid,n_s_grid), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate jfun in geom')

  ! allocate the array that contains the local major radius 
  allocate(Rfun(n_x_grid,n_s_grid), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate Rfun in geom') 

  ! allocate the array that contains the toroidal field fraction
  allocate(bt_frac(n_x_grid,n_s_grid), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate bt_frac in geom')

  ! allocate the array that contains the radial derivative of R
  allocate(kfun(n_x_grid,n_s_grid), stat = ierr)
  if (ierr /= 0) call gkw_abort('could not allocate kfun in geom') 

  ! allocate the metric_G (GLOBAL)
  allocate(metric_G(n_x_grid,n_s_grid,3,3),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate the metric_G in geom')

  ! allocate the jacobian_G (GLOBAL)
  allocate(jacobian_G(n_x_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate the jacobian_G in geom')

  ! allocate the metric (LOCAL)
  if (parallel_s) then
    allocate(metric(n_x_grid,1-3:ns+3,3,3),stat=ierr)
  else
    allocate(metric(n_x_grid,ns,3,3),stat=ierr)
  endif
  if (ierr /= 0) call gkw_abort('Could not allocate the metric in geom')

  allocate(dBdpsi(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dBdpsi in geom')

  allocate(dBds(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dBds in geom')

  allocate(Z_FS(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate Z_FS in geom')

  allocate(dRdpsi(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dRdpsi in geom')

  allocate(dRds(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dRds in geom')

  allocate(dZdpsi(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dZdpsi in geom')

  allocate(dZds(n_x_grid,n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dZds in geom')

  allocate(bups(n_x_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate bups in geom')

  allocate(dpfdpsi(n_x_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dpfdpsi in geom')

  allocate(alphakp(n_x_grid,n_s_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could ot allocate alphakp in geom') 

  allocate(alphak(n_x_grid,1-3:ns+3),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate alphak in geom')

  allocate(alphak_xbnd(n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate alphak_xbnd in geom')

  allocate(shift_end_grid(n_x_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate shift_end_grid in geom')

  allocate(tmp_miller(nx,2), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate tmp_miller_i in geom')

  allocate(de_miller(nx,2), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate de_miller_i in geom')
  
  allocate(vp_miller_i(nx), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate vp_miller_i in geom')
  
  allocate(fp_miller(nx,2), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate fp_miller in geom')
  
  allocate(tp_miller(nx,2), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate tp_miller in geom')
  
  allocate(dpdpsi_rot(n_x_grid,n_s_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate dpdpsi_rot in geom')
  
  allocate(dpds_rot(n_x_grid,n_s_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate dpds_rot in geom')

end subroutine geom_allocate


!------------------------------------------------------------------------------
!> Initializes the grids of geom
!------------------------------------------------------------------------------

subroutine geom_init_grids

  use control,      only : flux_tube, neoclassics
  use control,      only : radial_boundary_conditions, shift_metric
  use grid,         only : n_s_grid, n_x_grid, nperiod, psil, psih, lx 
  use constants,    only : pi
  use general,      only : gkw_warn
  use mpiinterface, only : root_processor

  real    :: sgrmax, R02
  integer :: i,ix
  real    :: sumdumD,sumdumH,sumdumI,sumdumE,sumdumF,sumdumG


  ! gridsize in 's' direction (along the field line),
  ! currently the same for all radial modes.
  sgrmax = real(nperiod) - 0.5E0

  ! The s-grid is the same for all geom types 
  do i  = 0, n_s_grid+1
    sgr(i) = -sgrmax + 2E0*sgrmax*(i-0.5E0)/n_s_grid
  end do 

  ! ints=ds/(\int_-sgrmax^sgrmax ds)=ds/(2*nperiod-1) is used for the flux surface average 
  ! WARNING: some parts of the code (e.g. normalise) use ints(1) on the 
  ! assumption that the grid spacings are all equal.  
  ! If this is changed, those sections will need to be updated
  do i = 1, n_s_grid
    ints(i) = 0.5*(sgr(i+1)-sgr(i-1))
  end do
  ints(:) = ints(:)/(2*nperiod - 1)

  ! Set the grid distance (equi-distant and the same for all radial
  ! surfaces.  
  sgr_dist = sgr(2)-sgr(1)

  ! initialise a variable that can be used when the global index at
  ! the LFS is needed
  isg_lfs = n_s_grid/2 + 1

  ! set up the radial grid
  if (flux_tube) then
    do ix = 1, n_x_grid
      xgr(ix) = eps
    end do
  else
    if ((psil <= 0) .and. (radial_boundary_conditions == 'Dirichlet')) then
      psil = psih / (2.E0*n_x_grid + 1.E0)
    endif
    do ix = 1, n_x_grid
      ! In the special case above, the following equation reduces to
      ! \f[
      !   xgr(ix) = psih \frac{ix}{n_xgrid + \frac{1}{2}}
      ! \f]
      xgr(ix) = psil + (psih-psil)*(ix-0.5E0)/real(n_x_grid)
    end do
  endif

  ! always set the radial grid size
  ! this version should be used everywhere (e.g. gyro_average, mode, linear_terms)
  dxgr = lx / n_x_grid

  ! calculate the q and shear profile
  call geom_profiles 
    
  ! Initialises Bref, Rref, beta and betaprime
  Rref          = -100.E0
  Bref          = -100.E0
  beta_eq       = -100.E0
  betaprime_eq  = -100.E0
  bt_frac(:,:)  = -100.E0
  r02 = -999.

  ! chease /slab calculations are done for one radial grid point only
  ix = 1

  ! metric and associated elements (magnetic field, beta,...)
  select case(geom_type)

  case('s-alpha')

    call geom_s_alpha

    ! apply shifted metric when non spectral
    call geom_shift_metric

    !                     (alleps,  gfun_num)
    call calc_geom_tensors(.false., .true.) 
    
    ! jump out (could deallocate helper arrays first) 
    !return    

  ! Experimental: s-alpha with "some" finite epsilon effects
  ! Warning: this geometry is not self-consistent and not recommended
  ! Flux surface averages are incorrect, 
  ! e_eps_zeta is not a flux function
  case('s-alpha-eps')

    call geom_s_alpha

    do i = 1,n_s_grid
      ! The derivatives of the magnetic field, put back higher order dependence
      dBdpsi(ix,i) = dBdpsi(ix,i) * (1 + eps*cos(2 * pi * sgr(i)))**(-2.0)
      dBds(ix,i)   = dBds(ix,i) * (1 + eps*cos(2 * pi * sgr(i)))**(-2.0)

      ! put back higher order terms in eps in metric

      ! the s s element 
      metric_G(ix,i,3,3) = 1.0/(2.E0*pi*xgr(ix))**2
    
      ! the zeta s element 
      metric_G(ix,i,2,3) = qx(ix)/(2*pi*xgr(ix))**2

      ! for the other elements symmetry applies 
      metric_G(ix,i,3,2) = metric_G(ix,i,2,3) 

    end do
    
    ! use the rotation values from circ
    select case(R0_loc)
      case('axis') !axis=place on flux surface at which R=Raxis
        R0 = 1.E0
        lfun = 0.
        jfunl = 2.*eps + eps*eps
        jfunh = -2.*eps + eps*eps 
                          
      case('LFS','Rmax') !Low field side in the plane of magnetic axis
        !Rmax deprecated
        R0 = 1.E0+eps
        lfun = 2.E0*(1.E0+eps)
        jfunl = 0. 
        jfunh = -4.*eps
                          
    end select

    ! apply shifted metric when non spectral
    call geom_shift_metric

    !                     (alleps,  gfun_num, test_efun)
    call calc_geom_tensors(.true., .false., .false.) 
    
    ! jump out (could deallocate helper arrays first) 
    !return 

  case('circ') 

    call geom_circ 

    ! apply shifted metric when non spectral 
    call geom_shift_metric

    ! the tensors         (all_eps, gfun_num)
    call calc_geom_tensors(.true., .false.) 

    ! jump out (could deallocate helper arrays first) 
    !return
  
  case('miller','fourier')
    
    call geom_miller
    
    call geom_shift_metric
    
    call calc_geom_tensors(.true., .false.)
    
    !return
    
  case('slab','slab_periodic') 
  
    ! Both are identical, except for the boundary conditions in the mode_box case
    
    call geom_slab
    
    call geom_shift_metric

    call calc_geom_tensors(.true., .false.)        
      
    if (shift_metric) call gkw_warn('slab geom not yet tested with shift_metric') 

  case('chease')   
    call geom_chease
    

  case default   
    ! No known option specified 
    write(*,*) 'You specified ',geom_type, & 
    & 'for geom_type in the namelist GEOM'
    write(*,*) 'Only known options: s-alpha, circ,' 
    write(*,*) 'chease, slab, slab_periodic, miller, fourier' 
    stop 1

  end select
  
  if(neoclassics)then
    sumdumD = 0.E0
    sumdumH = 0.E0
    sumdumI = 0.E0
    sumdumE = 0.E0
    sumdumF = 0.E0
    sumdumG = 0.E0
    do i = 1, n_s_grid
      sumdumD = sumdumD + ints(i)*dfun(1,i,1)
      sumdumH = sumdumH + ints(i)*hfun(1,i,1)
      sumdumI = sumdumI + ints(i)*ifun(1,i,1)
      sumdumE = sumdumE + ints(i)*efun(1,i,1,3)
      sumdumF = sumdumF + ints(i)*ffun(1,i)
      sumdumG = sumdumG + ints(i)*gfun(1,i)
    enddo
    if(root_processor)then
      write(*,*)'D integrated',sumdumD
      write(*,*)'H integrated',sumdumH
      write(*,*)'I integrated',sumdumI
      write(*,*)'E integrated',sumdumE
      write(*,*)'F integrated',sumdumF
      write(*,*)'G integrated',sumdumG
   endif 
   do i = 1, n_s_grid
      dfun(1,i,1) = dfun(1,i,1) - sumdumD
      hfun(1,i,1) = hfun(1,i,1) - sumdumH
      ifun(1,i,1) = ifun(1,i,1) - sumdumI
      gfun(1,i) = gfun(1,i) - sumdumG
    enddo
    sumdumD = 0.E0
    sumdumH = 0.E0
    sumdumI = 0.E0
    sumdumE = 0.E0
    sumdumF = 0.E0
    sumdumG = 0.E0
    do i = 1, n_s_grid
      sumdumD = sumdumD + ints(i)*dfun(1,i,1)
      sumdumH = sumdumH + ints(i)*hfun(1,i,1)
      sumdumI = sumdumI + ints(i)*ifun(1,i,1)
      sumdumE = sumdumE + ints(i)*efun(1,i,1,3)
      sumdumF = sumdumF + ints(i)*ffun(1,i)
      sumdumG = sumdumG + ints(i)*gfun(1,i)
    enddo
    if(root_processor)then
      write(*,*)'D integrated',sumdumD
      write(*,*)'H integrated',sumdumH
      write(*,*)'I integrated',sumdumI
      write(*,*)'E integrated',sumdumE
      write(*,*)'F integrated',sumdumF
      write(*,*)'G integrated',sumdumG
    endif
  endif
  ! copy the values into the arrays that are a function of radius 
  if (flux_tube) call distribute_geom_flux_tube(2)

end subroutine geom_init_grids


!------------------------------------------------------------------------------
!> This routine determines all the quantities necessary for the s-alpha 
!> geometry 
!> ----------------------------------------------------------------------------

subroutine geom_s_alpha 

  use grid,      only : n_x_grid, n_s_grid  
  use constants, only : pi 
  use general,   only : gkw_exit

  integer :: ix, i
  real    :: dum


  ! Calculate the normfactor for k_zeta  (for global runs this is evaluated
  ! at the centre of the box. equivalent to sqrt(g_zeta_zeta at s=0)
  ix = (n_x_grid + 1) / 2 
  kthnorm = qx(ix)  / ( 2 * pi * xgr(ix))
  ! Calculate the normfactor for k_psi
  kxnorm = sqrt(1.0)

  ! Loop over all the flux surfaces and grid points along the field 
  do ix = 1, n_x_grid 
  
    ! set the minimum and maximum magnetic field strength 
    ! these should be made a function of the radius 
    bmin(ix) = 1.E0 / (1.E0 + xgr(ix)) 
    bmax(ix) = 1.E0 / (1.E0 - xgr(ix)) 

    ! the contra variant component of the field 
    bups(ix) = signJ / ( 2.E0* pi * qx(ix) ) 

    ! radial derivative of poloidal flux 
    dpfdpsi(ix) = xgr(ix) / qx(ix) 

    do i = 1, n_s_grid 

      ! poloidal angle    
      pol_angle(ix,i) =  2.*pi*sgr(i)

      ! Magnetic field strength 
      dum = 1. + xgr(ix) * cos(2.E0*pi*sgr(i))
      bn_G(ix,i) = 1.E0 / dum 

      ! Set the toroidal fraction of the field (bt_frac is always positive)
      bt_frac(ix,i) = 1.
    
      ! The major radius (note this will be reset later to 1.0)  
      rfun(ix,i) = 1. + xgr(ix)*cos(2.0*pi*sgr(i))

      ! the Z-coordinate (not needed) 
      z_fs(ix,i) = xgr(ix)*sin(2.0*pi*sgr(i))

      ! the derivatives of the coordinates
      dRdpsi(ix,i) = cos(2.0*pi*sgr(i))
      dRds(ix,i)   = - xgr(ix)*2.0*pi*sin(2.0*pi*sgr(i))
      dZdpsi(ix,i) = sin(2.0*pi*sgr(i))
      dZds(ix,i)   = xgr(ix)*2.0*pi*cos(2.0*pi*sgr(i))


      ! The derivatives of the magnetic field  
      dBdpsi(ix,i) = - cos(2 * pi * sgr(i))
      dBds(ix,i)   = 2.E0 * pi * xgr(ix) * sin(2.0 * pi * sgr(i))

    
      ! Normalised: metric(i,2,2)=g_zeta_zeta*Rref**2
    
      ! the psi psi element 
      metric_G(ix,i,1,1) = 1.E0 

      ! the psi zeta element 
      metric_G(ix,i,1,2) = qx(ix) * shatx(ix) * sgr(i) / xgr(ix)*signB*signJ

      ! the psi s element 
      metric_G(ix,i,1,3) = 0. 

      ! the zeta zeta element 
      metric_G(ix,i,2,2) = ( qx(ix) / (2 * pi * xgr(ix) ))**2 * (1 +     &
                         & (2.E0*pi*sgr(i)*shatx(ix))**2)
      
      ! the zeta s element 
      metric_G(ix,i,2,3) = qx(ix) / (2 * pi * xgr(ix))**2 *signB*signJ

      ! the s s element 
      metric_G(ix,i,3,3) = 1.E0 / (2 * pi * xgr(ix))**2


      ! for the other elements symmetry applies 
      metric_G(ix,i,2,1) = metric_G(ix,i,1,2) 
      metric_G(ix,i,3,1) = metric_G(ix,i,1,3) 
      metric_G(ix,i,3,2) = metric_G(ix,i,2,3) 
   
    end do 
  end do 

  ! Set the low field side position and the drivative 
  select case(R0_loc)
  case('axis') !axis=place on flux surface at which R=Raxis
    R0   = 1.E0
    lfun = 0.E0
    jfunl =  2.E0*eps  ! = 0 in LFS case, 2*eps in axis case 
    jfunh = -2.E0*eps  ! = -4*eps in LFS, -2*eps in axis case

  case('LFS','Rmax') !Low field side in the plane of magnetic axis
    ! Rmax deprecated...
    ! R0 for output purposes only.
    R0=1.E0+eps
    lfun=2.E0
    jfunl =  0.E0      ! = 0 in LFS case, 2*eps in axis case 
    jfunh = -4.E0*eps  ! = -4*eps in LFS, -2*eps in axis case

  case default
    call gkw_exit('Geom: unknown R0_loc option. &
                  & Allowed are: axis, LFS')
  end select

end subroutine geom_s_alpha 


!------------------------------------------------------------------------------
!> This routine implements the analytical circular equilibrium as described in 
!> Lapillonne [PoP 16 032308 (2009). Note that the starting assumption for the 
!> model is that the poloidal flux psi(x) = psi(r) for circular flux surfaces
!> This is valid neglecting terms of order eps^2 in the case of small eps and 
!> beta. So in the metric tensors, only terms up to order eps should be kept 
!> for consistency. This is however a bit tricky as the results can involve 
!> multiplication of terms of order 1/eps with terms of order eps^2 (or even 
!> 1/eps^2 with eps^3). For simplicity the metric terms are not expanded in 
!> the parameter eps (i.e. the full expression is kept). 
!> Just remember that the results are only valid up to order eps !!!
!------------------------------------------------------------------------------

subroutine geom_circ 

  use grid,      only : n_x_grid, n_s_grid 
  use constants, only : pi 
  use general,   only : gkw_abort, gkw_exit

  real, allocatable :: dzetadeps(:) 

  integer :: ix, i, j, ierr
  real    :: dum, shift 

  ! allocate the help array 
  allocate(dzetadeps(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dzetadeps in geom')

  ! Calculate the normfactor for k_zeta (for global runs evaluated at the 
  ! centre of the box)  equivalent to sqrt(g_zeta_zeta at s=0)
  ix = (n_x_grid+1) / 2 
  kthnorm = (1.E0 /(2.E0*pi)/(1.E0+xgr(ix)))*(1.E0 +              &
          & (qx(ix)/xgr(ix))**2*(1-xgr(ix)**2))**0.5
  ! Calculate the normfactor for k_psi
  kxnorm = sqrt(1.0)

  do ix = 1, n_x_grid  

    ! set the minimum and maximum magnetic field strength 
    ! At magnetic axis, bn_G=1 (these must be made a function of radius 
    ! at some point) 
    dum = (1.E0+xgr(ix)**2/qx(ix)**2/(1.-xgr(ix)**2))**0.5
    bmin(ix) = dum / (1.E0 + xgr(ix)) 
    bmax(ix) = dum / (1.E0 - xgr(ix)) 

    ! the contra variant component of the field 
    bups(ix) = signJ / ( 2.E0* pi * qx(ix) * sqrt(1.0E0-xgr(ix)**2) ) 

    ! radial derivative of poloidal flux 
    dpfdpsi(ix) = xgr(ix) / ( qx(ix) * sqrt(1.-xgr(ix)**2) )  
 
    ! the poloidal angle and derivative of zeta towards eps are calculated first 
    do i = 1, n_s_grid 
      pol_angle(ix,i) = 2.0E0*pi*sgr(i) 
      do j = 1,5 !solves numerically theta+eps*sin(theta)= 2pi*s
        pol_angle(ix,i)=2E0*pi*sgr(i) - xgr(ix) * sin(pol_angle(ix,i))
      end do
    end do 
    if (n_s_grid.lt.2) then 
      dzetadeps(1) = 0.
    else 
      dum = ((1.E0-xgr(ix))/(1.E0+xgr(ix)))**0.5
      dzetadeps(1) = atan(dum*tan(pol_angle(ix,1)/2.E0))
      do i = 2, n_s_grid 
        ! This can cause a floating point underflow but for normal compilation
        ! it gets ignored and seems not to cause a slowdown
        dzetadeps(i) = atan(dum*tan(pol_angle(ix,i)/2.E0))
        do while (dzetadeps(i)<dzetadeps(i-1))
          dzetadeps(i) = dzetadeps(i) + pi
        end do
      end do 
      shift = pi*floor((dzetadeps(1)-pol_angle(ix,1)/2.E0)/pi)
      do i = 1, n_s_grid 
        dzetadeps(i) = dzetadeps(i) - shift
      end do 
    endif 
    do i = 1, n_s_grid
      dum = tan(pol_angle(ix,i)/2.E0)
      dzetadeps(i) = signB*signJ/pi * qx(ix)/xgr(ix) * (shatx(ix)*dzetadeps(i) - &
      & xgr(ix)/(1-xgr(ix)**2)**0.5 * dum / (1.E0+ dum**2 + xgr(ix)*(1.E0-dum**2)))
    end do

    do i = 1, n_s_grid 

      ! Magnetic field strength 
      dum = (1.E0+xgr(ix)**2/qx(ix)**2/(1.-xgr(ix)**2))**0.5
      bn_G(ix,i) = dum / (1.E0 + xgr(ix) * cos(pol_angle(ix,i)))

      !Set the toroidal fractions of the field
      !bt_frac should always be positive
      dum = (1.E0+xgr(ix)**2/qx(ix)**2/(1.-xgr(ix)**2))**0.5
      bt_frac(ix,i) = 1.E0 / dum

      ! the major radius 
      Rfun(ix,i)=1.E0 + xgr(ix)*cos(pol_angle(ix,i))

      ! derivative of bn with respect to eps and theta(=pol_angle)
      dBdpsi(ix,i) = bn_G(ix,i)*(dum*xgr(ix)*(1.E0-(1.E0-xgr(ix)**2)*shatx(ix)/ & 
                   & qx(ix))/qx(ix)**2 - cos(pol_angle(ix,i)) / (1.E0 + xgr(ix) *  &
                   & cos(pol_angle(ix,i))))
      dBds(ix,i)   = bn_G(ix,i)*xgr(ix)*sin(pol_angle(ix,i)) & 
                   & / (1.E0 + xgr(ix) * cos(pol_angle(ix,i)))
      ! Transform to (psi,s) 
      dBdpsi(ix,i) = dBdpsi(ix,i) - sin(pol_angle(ix,i)) * dBds(ix,i)              &
                   & / (1 + xgr(ix)*cos(pol_angle(ix,i))) 
      dBds(ix,i)   = 2.E0 * pi * dBds(ix,i) / (1. + xgr(ix)*cos(pol_angle(ix,i))) 


      ! the derivative of the major radius towards (psi,theta) 
      dRdpsi(ix,i) = cos(pol_angle(ix,i))
      dRds(ix,i)   = - xgr(ix)*sin(pol_angle(ix,i))
      ! Transform to (psi,s) 
      dRdpsi(ix,i) = dRdpsi(ix,i) - sin(pol_angle(ix,i)) * dRds(ix,i)              &
                   & / (1 + xgr(ix)*cos(pol_angle(ix,i))) 
      dRds(ix,i)   = 2.E0 * pi * dRds(ix,i) / (1. + xgr(ix)*cos(pol_angle(ix,i))) 

      ! The derivatives of the Z-coordinate [first (psi,theta) space]
      dZdpsi(ix,i) = sin(pol_angle(ix,i))
      dZds(ix,i)   = xgr(ix)*cos(pol_angle(ix,i))
      ! Transform to (psi,s) space
      dZdpsi(ix,i) = dZdpsi(ix,i) - sin(pol_angle(ix,i)) * dZds(ix,i)        &
                   & / (1 + xgr(ix)*cos(pol_angle(ix,i))) 
      dZds(ix,i)   = 2.E0 * pi * dZds(ix,i) / (1. + xgr(ix)*cos(pol_angle(ix,i))) 


      ! Normalised: metric(i,2,2)=g_zeta_zeta*Rref**2 
    
      ! the psi psi element 
      metric_G(ix,i,1,1) = 1.E0 
      
      ! the psi zeta element 
      metric_G(ix,i,1,2) = dzetadeps(i)

      ! the psi s element 
      metric_G(ix,i,1,3) = sin(pol_angle(ix,i))/(2.E0*pi)
      
      ! the zeta zeta element 
      metric_G(ix,i,2,2) = (1.E0 /(2.E0*pi)/(1.E0+xgr(ix)*cos(pol_angle(ix,i))))**2 * &
                       & (1.E0 + (1.E0-xgr(ix)**2)*(qx(ix)/xgr(ix))**2)+dzetadeps(i)**2
      
      ! the zeta s element 
      metric_G(ix,i,2,3) = qx(ix)/(2*pi*xgr(ix))**2 *signB*signJ*(1.E0-xgr(ix)**2)**0.5 &
                      & + dzetadeps(i)*sin(pol_angle(ix,i))/(2.E0*pi) 

      ! the s s element 
      metric_G(ix,i,3,3) = 1.E0 / (2.E0 * pi)**2 * &
                        & ((1.E0/xgr(ix) + cos(pol_angle(ix,i)))**2 + sin(pol_angle(ix,i))**2)
      
      ! for the other elements symmetry applies 
      metric_G(ix,i,2,1) = metric_G(ix,i,1,2) 
      metric_G(ix,i,3,1) = metric_G(ix,i,1,3) 
      metric_G(ix,i,3,2) = metric_G(ix,i,2,3) 
    
      z_fs(ix,i) = xgr(ix)*sin(pol_angle(ix,i))
    end do 
  end do    

  ! Centrifugal trapping
  ! jfun = R_N^2-R0_N^2
  ! kfun= djfun/dpsi AT CONSTANT s !!
  ! R0 (in Peeters PoP09 on centrifugal force) 
  ! specifies location of density gradient definitions
  select case(R0_loc)
  case('axis') !axis=place on flux surface at which R=Raxis
    R0 = 1.E0
    lfun = 0.
    jfunl = 2.*eps + eps*eps
    jfunh = -2.*eps + eps*eps 
                      
  case('LFS','Rmax') !Low field side in the plane of magnetic axis
    !Rmax deprecated
    R0 = 1.E0+eps
    lfun = 2.E0*(1.E0+eps)
    jfunl = 0. 
    jfunh = -4.*eps
                      
  case default
    call gkw_exit('Geom: unknown R0_loc option. Allowed are: axis, LFS')
  end select

end subroutine geom_circ 

!------------------------------------------------------------------------------
!> This routine implements the Miller and the Fourier  parametrisation following 
!> the methods used by J.Candy [Plasma Phys. Control. Fusion 51 (2009) 105009]
!------------------------------------------------------------------------------
subroutine geom_miller

  use grid,         only : n_x_grid, n_s_grid, nperiod
  use constants,    only : pi
  use general,      only : gkw_abort, gkw_exit, gkw_warn
  use control,      only : flux_tube
  use mpiinterface, only : root_processor

  integer :: ix, i, ierr, isa, N, ith, insh, dum1, dum2 = 0
  integer :: is, it, in, itest, n_x_loop
  integer :: ithetazero, nperiod_extended

  real    :: asindelta, x1, x2, x3, x4, dBdl, dBdrho, integ3, integ4
  real    :: integ5, integ6, integ7, Fprime, F, vol
  real    :: dvoldpsi, integ10, integ11, grdp = 0.0, dum3, dum4, dR0dpsi = 0.0
  real    :: qkappa = 0.0, qdelta = 0.0, qsquare = 0.0, qskappa = 0.0
  real    :: qsdelta = 0.0, qssquare = 0.0, qZmil = 0.0, qdRmil = 0.0
  real    :: qdZmil = 0.0, qgradp = 0.0
  real    :: an, dandpsi, dandth, d2andth2, d2andpsidth

  real, allocatable :: dRdpsi2(:), dZdpsi2(:), dRdth(:), dZdth(:), d2Rdth(:)
  real, allocatable :: d2Zdth(:), rc(:), d2Rdpsidth(:), d2Zdpsidth(:)
  real, allocatable :: gdwnpsipsi(:), gdwnpsith(:), gdwnthth(:), guppsipsi(:)
  real, allocatable :: guppsith(:), gupthth(:), abs_gradpsi(:), drhodpsi(:) 
  real, allocatable :: dldpsi(:), dldth(:), Bp(:), pf1(:), pf2(:) 
  real, allocatable :: dabs_gradpsi_dth(:), dpf1dl(:), zeta1(:), integ1(:) 
  real, allocatable :: integ2(:), integ8(:), Jpsi(:), dsdth(:), dzetadpsi(:)
  real, allocatable :: dzetadth(:), Bt(:), dsdpsi(:), cosu(:), sinu(:)
  real, allocatable :: theta_grid(:), s_grid(:), dBdpsi2(:), dRds2(:), dBds2(:)
  real, allocatable :: bn_G2(:), bt_frac2(:), rfun2(:), z_fs2(:) 
  real, allocatable :: metric_G2(:,:,:), dZds2(:), integ9(:), dsdpf(:), f1tmp(:)
  real, allocatable :: f2tmp(:),f3tmp(:), gradp_rot1(:), gradp_rot2(:)

  ierr = 0
  
  if (root_processor) write(*,*)
  if (root_processor) write(*,*) 'Performing numerical integrals for Miller/Fourier...'

  ! N: number of points in the theta_grid (must be odd number)
  ! 51 was suffcient after testing the convergence of the integrals used in 
  ! this routine (simpson's scheme), but 501 was needed to allow
  ! a reproducible gkw test case
  
  nperiod_extended=nperiod+1
  N = 501*(2*nperiod_extended - 1)
  dum3 = 0.
  dum4 = 0.

  allocate(f1tmp(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate f1tmp in geom')
  
  allocate(f2tmp(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate f2tmp in geom')
  
  allocate(f3tmp(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate f3tmp in geom')
  
  allocate(dsdpf(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dsdpf in geom')

  allocate(integ9(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate integ9 in geom')

  allocate(metric_G2(N,3,3),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate the metric_G2 in geom')

  allocate(Jpsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate Jpsi in geom')

  allocate(s_grid(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate s_grid in geom')

  allocate(theta_grid(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate theta_grid in geom')

  allocate(cosu(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate cosu in geom')

  allocate(sinu(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate sinu in geom')

  allocate(dldth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dldth in geom')

  allocate(dsdpsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dsdpsi in geom')

  allocate(Bt(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate Bt in geom')

  allocate(rc(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate rc in geom')

  allocate(dRdpsi2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dRdpsi in geom')
    
  allocate(dZdpsi2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dZdpsi2 in geom')

  allocate(dRdth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dRdth in geom')

  allocate(dZdth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dZdth in geom')

  allocate(d2Rdth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate d2Rdth in geom')

  allocate(d2Zdth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate d2Zdth in geom')

  allocate(d2Rdpsidth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate d2Rdpsidth in geom')

  allocate(d2Zdpsidth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate d2Zdpsidth in geom')

  allocate(gdwnpsipsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate gdwnpsipsi in geom')

  allocate(gdwnpsith(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate gdwnpsith in geom')

  allocate(gdwnthth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate gdwnthth in geom')

  allocate(guppsipsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate guppsipsi in geom')

  allocate(guppsith(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate guppsith in geom')

  allocate(gupthth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate gupthth in geom')

  allocate(abs_gradpsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate abs_gradpsi in geom')

  allocate(drhodpsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate drhodpsi in geom')

  allocate(dldpsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dldpsi in geom')

  allocate(Bp(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate Bp in geom')

  allocate(pf1(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate pf1 in geom')

  allocate(pf2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate pf2 in geom')

  allocate(dpf1dl(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dpf1dl in geom')

  allocate(zeta1(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate zeta1 in geom')

  allocate(dabs_gradpsi_dth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dabs_gradpsi_dth in geom')

  allocate(integ1(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate integ1 in geom')

  allocate(integ2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate integ2 in geom')

  allocate(integ8(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate integ8 in geom')

  allocate(dsdth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dsdth in geom')

  allocate(dzetadpsi(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dzetadpsi in geom')

  allocate(dzetadth(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dzetadth in geom')

  allocate(rfun2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate rfun2 in geom')

  allocate(z_fs2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate z_fs2 in geom')

  allocate(bt_frac2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate bt_frac2 in geom')    

  allocate(bn_G2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate bn_G2 in geom')    

  allocate(dBds2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dBds2 in geom')    

  allocate(dRds2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dRds2 in geom')

  allocate(dZds2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dZds2 in geom')

  allocate(dBdpsi2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dBdpsi2 in geom')
  
  allocate(gradp_rot1(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate gradp_rot1 in geom')

  allocate(gradp_rot2(1:N),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate gradp_rot2 in geom')

  if (flux_tube) then  ! all x_points are the same
    n_x_loop = 1       ! avoid repeating the Miller integrals
  else
    n_x_loop = n_x_grid    ! store the values of kappa, delta etc. This may not be necessary but 
    ! is done to prevent problems connected with this parameters being 
    ! redefined 
    if (geom_type == 'miller') then
      qkappa = kappa 
      qdelta = delta 
      qsquare = square 
      qskappa = skappa 
      qsdelta = sdelta 
      qssquare = ssquare 
      qZmil    = Zmil 
      qdRmil   = dRmil 
      qdZmil   = dZmil 
      qgradp   = gradp 
    endif
  end if

  do ix = 1, n_x_loop

    ! Integrals used to switch between the two gradp_type
    ! Allows calculation of the volume described by the flux surface of
    ! consideration and its derivative with respect to the poloidal flux
    integ3 = 0.    
    integ4 = 0.
    integ5 = 0.
    integ6 = 0.
    integ7 = 0.
    integ10 = 0.
    integ11 = 0.

    if (.not. flux_tube) then 
      kappa   = kappax(ix) 
      delta   = deltax(ix) 
      square  = squarex(ix) 
      skappa  = skappax(ix) 
      sdelta  = sdeltax(ix) 
      ssquare = ssquarex(ix)
      Zmil    = Zmilx(ix)
      dRmil   = dRmilx(ix)
      dZmil   = dZmilx(ix)
      gradp   = gradpx(ix) 
    endif 

    ! Initialisation of theta_grid 
    ! Here, use nperiod_extended=nperiod+1 to define a domain always larger than s_grid 
    ! and avoid extrapolation that can otherwise happen for up-down asymmetric flux surfaces
    do ith = 1, N

      theta_grid(ith) = -real(2*nperiod_extended-1)*(pi) + 2.*real(2*nperiod_extended-1)    &
      &  * (pi) * (ith-1) / (N-1) 

      f1tmp(ith) = exp(theta_grid(ith))

      gradp_rot1(ith) = 0.

      gradp_rot2(ith) = 0.
      
    end do

    ! Set beta_miller to beta_rota_miller or beta_ref
    if (geom_type=='miller' .and. beta_rota_miller_type == 'geom') then
   
      beta_miller(:) = beta_rota_miller

    end if

    ! integrate a simple exp function. this should yield the same exp
    ! again and allows us to check the accuracy
    call simpson_numerical_integ (theta_grid,f1tmp,N,integ1)

    ! Test if numerical integrals are accurate enough. 
    if (abs(integ1(N)-integ1(1) - (exp(theta_grid(N)) - exp(theta_grid(1))))   &
    &  > (2E-04 *(exp(theta_grid(N)) - exp(theta_grid(1))))) then
      ! The diff between integrated and original diff of first and
      ! last value, relative to the original diff, is larger than a
      ! small number.
      call gkw_warn ('Numerical integrals in geom miller/fourier are not accurate')
    end if

    asindelta = asin(delta)

    do ith = 1, N

      integ1(ith) = 0.
     if (geom_type=='miller') then
      ! precompute terms.
      ! the argument of the cos in the R expression:
      x1 = theta_grid(ith) + asindelta * sin(theta_grid(ith))
      x2 = 1. + asindelta * cos(theta_grid(ith))
      x3 = 1. + 2. * square * cos(2.*theta_grid(ith))
      ! the argument of the sin in the Z expression:
      x4 = theta_grid(ith) + square * sin(2.*theta_grid(ith))

      ! Parametrisation of the flux surface with R = rfun2 and Z = z_fs2
      rfun2(ith) = 1. + xgr(ix) * cos(x1)

      z_fs2(ith) = Zmil + kappa * xgr(ix) * sin(x4)

      ! First derivatives of rfun2 and z_fs2
      dRdpsi2(ith) =  dRmil + cos(x1) - sdelta * sin(theta_grid(ith)) * sin(x1)
    
      dZdpsi2(ith) = dZmil + kappa * sin(x4) * (1. + skappa) +  kappa * &
        & ssquare * sin(2.*theta_grid(ith)) * cos(x4)
      
      dRdth(ith) = -xgr(ix) * x2 * sin(x1)

      dZdth(ith) = kappa * xgr(ix) * x3 * cos(x4)

      ! Second derivatives of rfun2 and z_fs2
      d2Rdth(ith) = xgr(ix) * asindelta * sin(theta_grid(ith)) * sin(x1) -           &
      &  xgr(ix) * (x2)**2. * cos(x1)

      d2Zdth(ith) = -4. * kappa * xgr(ix) * square * sin(2.*theta_grid(ith))    &
      &  * cos(x4) - kappa * xgr(ix) * (x3)**2. * sin(x4)

      d2Zdpsidth(ith) = kappa * (1. + skappa) * x3 * cos(x4) + 2. *          &
      &  kappa * ssquare * cos(2. * theta_grid(ith)) * cos(x4) -             &
      &  kappa * ssquare * sin(2.*theta_grid(ith)) * x3 * sin(x4)

      d2Rdpsidth(ith) = -x2 * sin(x1) - sdelta * cos(theta_grid(ith)) *      &
      &  sin(x1) - sdelta * x2 * sin(theta_grid(ith)) * cos(x1)
    endif

    if (geom_type=='fourier') then
      ! computes the minor radius and its derivatives      
      an=0
      dandpsi=0
      dandth=0
      d2andpsidth=0
      d2andth2=0
      do insh = 1,N_shape
        an = an + c(insh) * cos((insh-1)*theta_grid(ith)) + s(insh)*sin((insh-1)*theta_grid(ith))
        dandpsi = dandpsi + c_prime(insh)*cos((insh-1)*theta_grid(ith))             &
        &                 + s_prime(insh)*sin((insh-1)*theta_grid(ith))
        dandth  = dandth  - (insh-1)*c(insh)* sin((insh-1)*theta_grid(ith))    &
        &                 + (insh-1)*s(insh)* cos((insh-1)*theta_grid(ith)) 
        d2andpsidth = d2andpsidth - (insh-1)*c_prime(insh)* sin((insh-1)*theta_grid(ith))    &
        &                         + (insh-1)*s_prime(insh)* cos((insh-1)*theta_grid(ith)) 
        d2andth2    = d2andth2    - (insh-1)*(insh-1)*c(insh)* cos((insh-1)*theta_grid(ith))    &
        &                         - (insh-1)*(insh-1)*s(insh)* sin((insh-1)*theta_grid(ith))
      enddo

      ! Parametrisation of the flux surface with R = rfun2 and Z = z_fs2
      rfun2(ith) = 1 + an*cos(theta_grid(ith))
      z_fs2(ith) = an*sin(theta_grid(ith))

      ! Zmil used to compute the poloidal angle
      Zmil = 0

      ! First derivatives of rfun2 and z_fs2
      dRdpsi2(ith) = dandpsi*cos(theta_grid(ith))
    
      dZdpsi2(ith) = dandpsi*sin(theta_grid(ith))
      
      dRdth(ith) = dandth*cos(theta_grid(ith)) - an*sin(theta_grid(ith))

      dZdth(ith) = dandth*sin(theta_grid(ith)) + an*cos(theta_grid(ith))

      ! Second derivatives of rfun2 and z_fs2
      d2Rdth(ith) = d2andth2*cos(theta_grid(ith)) - 2*dandth*sin(theta_grid(ith)) - an*cos(theta_grid(ith))

      d2Zdth(ith) = d2andth2*sin(theta_grid(ith)) + 2*dandth*cos(theta_grid(ith)) - an*sin(theta_grid(ith))

      d2Rdpsidth(ith) = d2andpsidth*cos(theta_grid(ith)) - dandpsi*sin(theta_grid(ith))

      d2Zdpsidth(ith) = d2andpsidth*sin(theta_grid(ith)) + dandpsi*cos(theta_grid(ith))

    endif

      ! Definition of the Jacobian (psi, theta, phi)      
      Jpsi(ith) = rfun2(ith) * (dRdpsi2(ith) * dZdth(ith) - dRdth(ith) *     &
      &  dZdpsi2(ith))
      
      ! Definition of the metric tensor's elements
      gdwnpsipsi(ith) = dRdpsi2(ith)**2. + dZdpsi2(ith)**2.
      gdwnthth(ith) = dRdth(ith)**2. + dZdth(ith)**2.
      gdwnpsith(ith) = dRdpsi2(ith) * dRdth(ith) + dZdpsi2(ith) * dZdth(ith)

      ! Contravariant elements of the metric tensor
      guppsipsi(ith) = (gdwnthth(ith) * rfun2(ith)**2.) / (Jpsi(ith))**2.
      gupthth(ith) = (gdwnpsipsi(ith) * rfun2(ith)**2.) / (Jpsi(ith))**2.
      guppsith(ith) = - (gdwnpsith(ith) * rfun2(ith)**2.) / (Jpsi(ith))**2.

      ! Some relations between the Mercier-Luc coordinate system and 
      ! (psi, theta, phi)
      dldth(ith) = sqrt(gdwnthth(ith))

      cosu(ith) = dZdth(ith) / dldth(ith)
      sinu(ith) = - dRdth(ith) / dldth(ith)

      drhodpsi(ith) = cosu(ith) * dRdpsi2(ith) + sinu(ith) * dZdpsi2(ith)
      dldpsi(ith) = cosu(ith) * dZdpsi2(ith) - sinu(ith) * dRdpsi2(ith)
      
      ! Radius of curvature
      rc(ith) = gdwnthth(ith)**(1.5) / (dRdth(ith) * d2Zdth(ith) -           &
      &  dZdth(ith) * d2Rdth(ith))
  
      abs_gradpsi(ith) = sqrt(gdwnthth(ith)) / (dRdpsi2(ith) * dZdth(ith) -  &
      &  dRdth(ith) * dZdpsi2(ith))
      
    end do

    ! Definition of the integrals needed for the relation between alpha and p'  
    do ith=1, N

      f1tmp(ith) = z_fs2(ith) * dRdth(ith)

      f2tmp(ith) = rfun2(ith) * sqrt(dRdth(ith)**2. + dZdth(ith)**2.)

      f3tmp(ith) = sqrt(dRdth(ith)**2. + dZdth(ith)**2.)

    end do
    
    call simpson_numerical_integ (theta_grid,f1tmp,N,integ1)
    call simpson_numerical_integ (theta_grid,f2tmp,N,integ2)
    call simpson_numerical_integ (theta_grid,f3tmp,N,integ8)
    
    integ3 = abs((-integ1(1)+integ1(N))/real(2*nperiod_extended-1))
    integ4 = (-integ2(1)+integ2(N))/real(2*nperiod_extended-1)
    integ10 = (-integ8(1)+integ8(N))/real(2*nperiod_extended-1)
    integ4 = integ4 / integ10

    do ith = 1, N
    
      f1tmp(ith) = dRdpsi2(ith) * sqrt(dRdth(ith)**2. + dZdth(ith)**2.)
      f2tmp(ith) = rfun2(ith) * (d2Rdpsidth(ith) * dRdth(ith) +             &
      &  d2Zdpsidth(ith) * dZdth(ith)) / sqrt(dRdth(ith)**2. + dZdth(ith)**2.)
      f3tmp(ith) =  dZdpsi2(ith) * dRdth(ith) + z_fs2(ith) * d2Rdpsidth(ith)

    end do

    call simpson_numerical_integ (theta_grid,f1tmp,N,integ1)
    call simpson_numerical_integ (theta_grid,f2tmp,N,integ2)
    call simpson_numerical_integ (theta_grid,f3tmp,N,integ8)

    integ5 = (-integ1(1)+integ1(N))/real(2*nperiod_extended-1)
    integ6 = (-integ2(1)+integ2(N))/real(2*nperiod_extended-1)
    integ7 = (-integ8(1)+integ8(N))/real(2*nperiod_extended-1)

    do ith = 1, N    
      
      f1tmp(ith) = Jpsi(ith) / rfun2(ith)**2.
      f2tmp(ith) = (d2Rdpsidth(ith) * dRdth(ith) + d2Zdpsidth(ith) *         &
      &  dZdth(ith)) / sqrt(dRdth(ith)**2. + dZdth(ith)**2.)      
 
    end do
  
    call simpson_numerical_integ (theta_grid,f2tmp,N,integ8)
    
    ! Integration of Jpsi and Jpsi/R^2 used for to determine s_grid and dpfdpsi
    call simpson_numerical_integ (theta_grid,Jpsi,N,integ1)
    call simpson_numerical_integ (theta_grid,f1tmp,N,integ2)
     
    integ11 = (-integ8(1)+integ8(N))/real(2*nperiod_extended-1)   

    ! Definition of s_grid which corresponds to theta_grid
    do ith = 1, N

      s_grid(ith) = integ1(ith) * real(2*nperiod_extended-1) / (-integ1(1) +          &
      &  integ1(N))

    end do

    ! Normalised F  
    F =  1.

    ! Volume defined by the flux surface  
    vol = 2.*pi*integ3 * integ4

    ! Radial derivative of the volume  
    dvoldpsi = 2.*pi* ((integ5 + integ6)/integ10-integ11*integ4/integ10) *   &
    &  integ3 + 2.*pi* integ4 * integ7

    dpfdpsi(ix) =  F * (-integ2(1)+integ2(N)) / (real(2*nperiod_extended-1)*          &
    &  2.*pi * qx(ix))   
    call interpquad (s_grid, rfun2, N, n_s_grid, sgr(1:n_s_grid), rfun(ix,:))
    call interpquad (s_grid, z_fs2, N, n_s_grid, sgr(1:n_s_grid), z_fs(ix,:))

    ! Definition of the poloidal angle (R0 is needed for gradp_type = rota_miller)
    do ith = 1, n_s_grid

      pol_angle(ix,ith) = atan((z_fs(ix,ith)-Zmil)/(rfun(ix,ith)-1.))

    end do

    ! Adjustment of the poloidal angle [-pi/2, pi/2] --> [-(2nperiod - 1)pi, (2nperiod -1)pi]
    itest = 1
    dum1 = 1
    in = -int(2*nperiod-1)

    do is = 1, 2*int(2*nperiod-1)

      do while (pol_angle(ix,itest)<= 0.)
        itest = itest + 1
        if (itest > n_s_grid) call gkw_abort('Need larger n_s_grid ?')
      end do

      do while (pol_angle(ix,itest) > 0.)

        dum2 = minloc(abs(pol_angle(ix,dum1:itest) - pi/2.),1)
        
        if (itest == n_s_grid) then

          exit

        end if
        itest = itest + 1

      end do

      dum2 = dum2 + dum1-1

      do it = dum1, dum2

        pol_angle(ix,it) = pol_angle(ix,it) + real(in * pi)

      end do
      in = in + 1
      dum1 = dum2 + 1

    end do

    do it = dum1, n_s_grid

      pol_angle(ix,it) = pol_angle(ix,it) + real(in * pi)

    end do 
    
    select case(R0_loc)
    case('axis') 
      isa = minloc(abs(pol_angle(ix,1:n_s_grid)-pi/2.),1)
      R0 = rfun(ix,isa)
      dR0dpsi = 0.
      call gkw_warn('Axis location inexact with miller_geom')

    case('LFS')      
      isa = minloc(abs(theta_grid(1:N)),1)      
      R0 = rfun2(isa)
      dR0dpsi = dRdpsi2(isa)

    case default
      call gkw_exit('Geom: unknown R0_loc option. Allowed are: axis, LFS')
    end select

    ! Switch between the different possible inputs for gradp (alpha, p',beta_prime,
    ! alpha_mhd, rota_miller) 
    select case(gradp_type)
      case('alpha')     ! From Miller paper
        grdp = (gradp * (4.*pi**2.)*dpfdpsi(ix)*sqrt(2.*pi**2.)) /       &
        &  (dvoldpsi * sqrt(vol)) 

      case('alpha_mhd') ! alpha = -q^2 * R * dbeta^e_dr - see e.g. Candy 2009
         grdp = - gradp *dpfdpsi(ix) / eps**2.   
     
      case('pprime')   ! Pressure derivative wrt poloidal flux
         if (gradp > 0.) then
           call gkw_warn('gradp > 0 corresponds to a hollow pressure profile ! ')
         end if
         grdp = gradp
                 
      case('rota_miller')     
           do ith = 1, N 
             
             dum3 = tmp_miller(ix,1)+tmp_miller(ix,2)
             dum4 = tmp_miller(ix,1)*tp_miller(ix,1)+tmp_miller(ix,2)*tp_miller(ix,2)

             gradp_rot1(ith)=(de_miller(ix,1)*dum3*(mas_miller_i/                     &            
             &  dum3*((rfun2(ith)**2.-R0**2.)*(-2.*vp_miller_i(ix)*vcor_miller+        &
             &  vcor_miller**2.*dum4/dum3)-2.*R0*vcor_miller**2.*dR0dpsi))*beta_miller(ix) &
             &  + betaprime_miller(ix))/(2.*dpfdpsi(ix))*                                  &
             &  exp(mas_miller_i*vcor_miller**2.*(rfun2(ith)**2.-R0**2.)/dum3)
         end do

      case('beta_prime')  ! use betaprime consistent with components
         if (betaprime_miller(ix) > 1e6) call gkw_abort('error in betaprime_miller')
         grdp = betaprime_miller(ix) / (2.*dpfdpsi(ix))
 
      case('beta_prime_input') !  Interpret gradp as beta_prime
         grdp = gradp / (2.*dpfdpsi(ix))
         
      case default
        call gkw_abort('unkown gradp_type: use "alpha", "pprime" or "beta_prime"')                      
    end select

    bups(ix) = signJ * real(2*nperiod_extended-1) * dpfdpsi(ix) / (-integ1(1) +       &
    &  integ1(N)) 

    do ith = 1, N

      ! Toroidal magnetic field
      Bt(ith) =  F / (rfun2(ith))

      ! Poloidal magnetic field
      Bp(ith) = (dpfdpsi(ix) * abs_gradpsi(ith)) / (rfun2(ith))

      ! Normalised magnetic field    
      bn_G2(ith) = sqrt(Bp(ith)**2. + Bt(ith)**2.)

      pf1(ith) = rfun2(ith) * Bp(ith)
    
      dabs_gradpsi_dth(ith) = (d2Rdth(ith) * dRdth(ith) + d2Zdth(ith) *     &
      &  dZdth(ith)) * rfun2(ith) / (Jpsi(ith) * dldth(ith)) -              &
      &  dldth(ith) * rfun2(ith)**2. * (d2Rdpsidth(ith) * dZdth(ith) +      &
      &  dRdpsi2(ith) * d2Zdth(ith) - d2Rdth(ith) * dZdpsi2(ith) -          &
      &  dRdth(ith) * d2Zdpsidth(ith)) / Jpsi(ith)**2.


      dpf1dl(ith) = (dpfdpsi(ix) * dabs_gradpsi_dth(ith) / dldth(ith))
      
      ! Definition of the integrals needed to determine F' and then grad(zeta)
      f1tmp(ith) = dldth(ith)* signB*signJ * (1./(rc(ith) * rfun2(ith) *    &
          &  Bp(ith)*pi) - cosu(ith) / (rfun2(ith)**2.* Bp(ith)*pi)) *F /   &
          &  (rfun2(ith)**2. * Bp(ith))
          
      f2tmp(ith) = dldth(ith) * signB*signJ  * (F/(rfun2(ith)**2.*Bp(ith)**2.) &
          &  + 1./F)/(2.*pi*rfun2(ith)**2.*Bp(ith))
      
      select case (gradp_type)
        case ('rota_miller')
      
          f3tmp(ith) = dldth(ith)* signB*signJ * F/ (rfun2(ith)**2.*         &
          &  Bp(ith)**3.*2.*pi) * gradp_rot1(ith)

        case default
      
          f3tmp(ith) = dldth(ith)* signB*signJ * F / (rfun2(ith)**2.*         &
          &  Bp(ith)**3.*2.*pi)
      end select
      
    end do

    call simpson_numerical_integ (theta_grid,f1tmp,N,integ8)
    call simpson_numerical_integ (theta_grid,f2tmp,N,integ9)
    call simpson_numerical_integ (theta_grid,f3tmp,N,integ1)
    
    select case (gradp_type)

      case ('rota_miller')
        ! Relation between F', the magnetic shear and the pressure gradient
        Fprime = (signB*signJ*real(2*nperiod_extended-1)*qx(ix) * shatx(ix) /          &
        &  (xgr(ix) * dpfdpsi(ix)) - (-integ8(1)+integ8(N)) - (-integ1(1) +   &
        &  integ1(N))) / ((-integ9(1)+integ9(N)) *F)
    
      case default    
        
        Fprime = (signB*signJ*real(2*nperiod_extended-1)*qx(ix) * shatx(ix) /          &
        &  (xgr(ix) * dpfdpsi(ix)) - (-integ8(1)+integ8(N)) - (-integ1(1) +   &
        &  integ1(N)) * grdp) / ((-integ9(1)+integ9(N)) *F)

    end select

    do ith = 1, N

      dzetadth(ith) = signB*signJ*F*dldth(ith)/(rfun2(ith)*pf1(ith)*2.*pi)

      select case (gradp_type)

        case ('rota_miller')
          zeta1(ith) = pf1(ith) * (integ8(ith) + integ9(ith) * Fprime * F +  &
          &  integ1(ith))

          pf2(ith) = 0.5 * (Bp(ith) * (cosu(ith) - rfun2(ith) / rc(ith)) -   &
          &  (rfun2(ith)**2.) * gradp_rot1(ith) - F*Fprime)
      
        case default
          zeta1(ith) = pf1(ith) * (integ8(ith) + integ9(ith) * Fprime * F +   &
          &  integ1(ith) * grdp)

          pf2(ith) = 0.5 * (Bp(ith) * (cosu(ith) - rfun2(ith) / rc(ith)) -    &
          &  (rfun2(ith)**2.) * grdp - F*Fprime)
      
      end select

      ! Definition of the integrals needed for grad(s)
      f1tmp(ith) = signJ*dldth(ith)*rfun2(ith)/pf1(ith)
      
      f2tmp(ith) = signJ*dldth(ith)*(cosu(ith)+rfun2(ith)/rc(ith) -       &
          &  2.*rfun2(ith)*pf2(ith)/pf1(ith))/pf1(ith)**2.
      
    end do

    call simpson_numerical_integ (theta_grid,f1tmp,N,integ9)
    call simpson_numerical_integ (theta_grid,f2tmp,N,integ8)    

    do ith = 1, N

      dsdth(ith) = signJ*rfun2(ith)*dldth(ith)*real(2*nperiod_extended-1) /        &
      &  (pf1(ith) * (-integ9(1)+integ9(N)))

      if (gradp_type == 'rota_miller') then
        
          gradp_rot2(ith) = mas_miller_i*vcor_miller**2.* &
          &  de_miller(ix,1)*exp(mas_miller_i*(rfun2(ith)**2.-R0**2.)*  &
          &  vcor_miller**2./dum3)*beta_miller(ix)  

      end if

      dsdpf(ith) = (integ8(ith)*real(2*nperiod_extended-1)/(-integ9(1) +           &
      &  integ9(N))-(-integ8(1)+integ8(N))*integ9(ith)*real(2*nperiod_extended-1)/ &
      &  (-integ9(1)+integ9(N))**2.)

      dzetadpsi(ith) = dzetadth(ith) * dldpsi(ith)/dldth(ith) + zeta1(ith) * drhodpsi(ith)

      dsdpsi(ith) = dsdpf(ith)*dpfdpsi(ix) + dsdth(ith)*dldpsi(ith)/dldth(ith)
  
      ! Definition of the contravariant elements of the metric tensor (psi, zeta, s)  
      ! the psi zeta element 
      metric_G2(ith,1,2) = dzetadpsi(ith) * guppsipsi(ith) + dzetadth(ith) *  &
      &  guppsith(ith)

      ! the psi s element 
      metric_G2(ith,1,3) = dsdth(ith) * guppsith(ith) + dsdpsi(ith) *         &
      &  guppsipsi(ith)

      ! the zeta zeta element
      metric_G2(ith,2,2) = (dzetadpsi(ith)**2.) * guppsipsi(ith) +            & 
      &  (dzetadth(ith)**2.) * gupthth(ith) + 1./(rfun2(ith)**2.*4.*pi**2.)   &
      &  +  2.*dzetadpsi(ith) * dzetadth(ith)*guppsith(ith) 

      ! the zeta s element 
      metric_G2(ith,2,3) = dzetadpsi(ith) * dsdpsi(ith) * guppsipsi(ith) +    &
      &  dzetadth(ith) * dsdth(ith) * gupthth(ith) + (dsdpsi(ith) *           &
      &  dzetadth(ith) + dsdth(ith) * dzetadpsi(ith)) * guppsith(ith)

      ! the s s element 
      metric_G2(ith,3,3) = dsdpsi(ith)**2. * guppsipsi(ith) + dsdth(ith)**2. *&
      &  gupthth(ith) + 2.*dsdpsi(ith)*dsdth(ith)*guppsith(ith) 

      ! the psi psi element 
      metric_G2(ith,1,1) = guppsipsi(ith)

      ! for the other elements symmetry applies 
      metric_G2(ith,2,1) = metric_G2(ith,1,2)
      metric_G2(ith,3,1) = metric_G2(ith,1,3)
      metric_G2(ith,3,2) = metric_G2(ith,2,3)

      dBdl = 0.5 * (2.*pf1(ith) * dpf1dl(ith) + 2.*sinu(ith) *          &
      &  pf1(ith)**2. / rfun2(ith) -2.*F**2. * dRdth(ith) /             & 
      &  (dldth(ith) * rfun2(ith))) / (bn_G2(ith) * rfun2(ith)**2.) 

      dBds2(ith) = dBdl * dldth(ith) / dsdth(ith)   

      bt_frac2(ith) = Bt(ith) / sqrt(Bp(ith)**2. + Bt(ith)**2.)

      dBdrho = 0.5 *(-2.*F**2. * cosu(ith) / rfun2(ith)**3. + 2.* F *   &
      &  Fprime * Bp(ith) / rfun2(ith) - 2.*Bp(ith)**2. * cosu(ith) /   &
      &  rfun2(ith) + 4.*Bp(ith) * pf2(ith) / rfun2(ith)) / bn_G2(ith)

      dBdpsi2(ith) = drhodpsi(ith) * dBdrho + dldpsi(ith) * dBdl
      
      dBdpsi2(ith) = dBdpsi2(ith) - dBds2(ith) * dsdpsi(ith)

      dRds2(ith) =  dRdth(ith) / dsdth(ith)

      dZds2(ith) = dZdth(ith) / dsdth(ith)
      
      dRdpsi2(ith) = dRdpsi2(ith) - dRds2(ith) * dsdpsi(ith)
      
      dZdpsi2(ith) = dZdpsi2(ith) - dZds2(ith) * dsdpsi(ith)             
      
    end do

    ! kthrho is the projection at the LFS, so take these at theta=0
    ithetazero = minloc(abs(theta_grid),1)
    kthnorm = sqrt(metric_G2(ithetazero,2,2))
    kxnorm = sqrt(metric_G2(ithetazero,1,1))
    
    ! Quadratic interpolation 

    ! Interpolation of the outputs to match s_grid

    call interpquad (s_grid, bt_frac2, N,n_s_grid, sgr(1:n_s_grid), bt_frac(ix,:))
    call interpquad (s_grid, bn_G2, N,n_s_grid, sgr(1:n_s_grid), bn_G(ix,:))
    call interpquad (s_grid,dBds2, N,n_s_grid, sgr(1:n_s_grid), dBds(ix,:))
    call interpquad (s_grid, dRds2, N,n_s_grid, sgr(1:n_s_grid), dRds(ix,:))
    call interpquad (s_grid, dZds2, N,n_s_grid, sgr(1:n_s_grid), dZds(ix,:))
    call interpquad (s_grid, dBdpsi2, N,n_s_grid, sgr(1:n_s_grid), dBdpsi(ix,:))
    call interpquad (s_grid, metric_G2(:,1,2), N,n_s_grid, sgr(1:n_s_grid), metric_G(ix,:,1,2))
    call interpquad (s_grid, metric_G2(:,1,3), N,n_s_grid, sgr(1:n_s_grid), metric_G(ix,:,1,3))
    call interpquad (s_grid, metric_G2(:,2,2), N,n_s_grid, sgr(1:n_s_grid), metric_G(ix,:,2,2))
    call interpquad (s_grid, metric_G2(:,2,3), N,n_s_grid, sgr(1:n_s_grid), metric_G(ix,:,2,3))
    call interpquad (s_grid, metric_G2(:,3,3), N,n_s_grid, sgr(1:n_s_grid), metric_G(ix,:,3,3))
    call interpquad (s_grid, metric_G2(:,1,1), N,n_s_grid, sgr(1:n_s_grid), metric_G(ix,:,1,1))
    call interpquad (s_grid, dZdpsi2, N,n_s_grid, sgr(1:n_s_grid), dZdpsi(ix,:))
    call interpquad (s_grid, dRdpsi2, N,n_s_grid, sgr(1:n_s_grid), dRdpsi(ix,:))
    call interpquad (s_grid, gradp_rot1, N,n_s_grid, sgr(1:n_s_grid), dpdpsi_rot(ix,:))
    call interpquad (s_grid, gradp_rot2, N,n_s_grid, sgr(1:n_s_grid), dpds_rot(ix,:))

    do i = 1, n_s_grid

      metric_G(ix,i,2,1) = metric_G(ix,i,1,2) 
      metric_G(ix,i,3,1) = metric_G(ix,i,1,3) 
      metric_G(ix,i,3,2) = metric_G(ix,i,2,3) 
      
    end do

    lfun = 2.*R0*dR0dpsi 
    select case(R0_loc)
    case('axis') 
      jfunl = maxval(rfun2(1:N)*rfun2(1:N)-R0*R0)
      jfunh = minval(rfun2(1:N)*rfun2(1:N)-R0*R0)
      call gkw_warn('Axis location inexact with miller_geom')
    case('LFS')   
      jfunl = 0.
      jfunh = minval(rfun2(1:N)*rfun2(1:N)-R0*R0)  
    case default
      call gkw_exit('Geom: unknown R0_loc option. Allowed are: axis, LFS')
    end select

    if (geom_type == 'miller') then !this test assumes up-down symmetry, not applicable for fourier
      !bmin(ix) = sqrt((1./(1.+eps))**2.+(dpfdpsi(ix)/((1.+eps)*(dRmil+1)))**2.)
      bmin(ix) = sqrt((1./(1.+xgr(ix)))**2.+(dpfdpsi(ix)/((1.+xgr(ix))*maxval(dRdpsi2(1:N))))**2.)
      if(bmin(ix) > minval(bn_G(ix,:))) then
        write(*,*) 'bmin:', ix, bmin(ix), minval(bn_G(ix,:))
        bmin(ix) = minval(bn_G(ix,:))
        call gkw_warn('Geom: wrong bmin for miller')
      end if

      !bmax(ix) = sqrt((1./(1.-eps))**2.+(dpfdpsi(ix)/((1.-eps)*(dRmil-1)))**2.)
      bmax(ix) = sqrt((1./(1.-xgr(ix)))**2.+(dpfdpsi(ix)/((1.-xgr(ix))*minval(dRdpsi2(1:N))))**2.)
      if(bmax(ix) < maxval(bn_G(ix,:))) then
        write(*,*) 'bmax:', ix, bmax(ix), maxval(bn_G(ix,:))
        bmax(ix) = maxval(bn_G(ix,:))
        call gkw_warn('Geom: wrong bmax for miller')
      end if
    else if (geom_type == 'fourier') then 
      bmin(ix) = minval(bn_G2(:))
      bmax(ix) = maxval(bn_G2(:))
    end if
  end do

  deallocate(dsdpf);     deallocate(integ9);     deallocate(metric_G2)
  deallocate(dBdpsi2);   deallocate(dZds2);      deallocate(dRds2)
  deallocate(dBds2);     deallocate(bn_G2);      deallocate(bt_frac2)
  deallocate(rfun2);     deallocate(z_fs2);      deallocate(s_grid)
  deallocate(cosu);      deallocate(sinu);       deallocate(integ8)
  deallocate(dldth);     deallocate(dsdpsi);     deallocate(Bt)
  deallocate(rc);        deallocate(dRdpsi2);    deallocate(dZdpsi2)
  deallocate(dRdth);     deallocate(dZdth);      deallocate(d2Rdth)
  deallocate(d2Zdth);    deallocate(d2Rdpsidth); deallocate(d2Zdpsidth)
  deallocate(gdwnpsipsi);deallocate(gdwnpsith);  deallocate(gdwnthth)
  deallocate(guppsipsi); deallocate(guppsith);   deallocate(gupthth)
  deallocate(drhodpsi);  deallocate(dldpsi);     deallocate(abs_gradpsi)
  deallocate(Bp);        deallocate(pf1);        deallocate(pf2)
  deallocate(dpf1dl);    deallocate(zeta1);      deallocate(dabs_gradpsi_dth)
  deallocate(integ1);    deallocate(integ2);     deallocate(Jpsi)
  deallocate(dsdth);     deallocate(dzetadpsi);  deallocate(dzetadth)
  deallocate(gradp_rot1);deallocate(gradp_rot2)

  ! copy the values into the arrays that are a function of radius 
  if (flux_tube) then 
    call distribute_geom_flux_tube(1)
    bmin = bmin(1)
    bmax = bmax(1)
  else 
   if ( geom_type == 'miller') then
     kappa   = qkappa 
     delta   = qdelta 
     square  = qsquare 
     skappa  = qskappa 
     sdelta  = qsdelta 
     ssquare = qssquare
     Zmil    = qZmil
     dRmil   = qdRmil 
     dZmil   = qdZmil 
     gradp   = qgradp 
   endif
  endif 
  
  if (root_processor) then
    write(*,*) '                                           ...done'
  end if
        
end subroutine geom_miller



!------------------------------------------------------------------------------
!> This routine implements a simple slab geometry
!> Here we use the sheared slab coordinates
!> x' = x
!> y' = y - shat*x*z
!> s = z
!> As used in Newton et. al PPCF (2010) and Hammett et. al. PPCF (1993)
!> The eps and q inputs no longer have any meaning
!> while s is the shat input number (which no longer has a q in
!> its meaning)
!> however eps is still used in mode.f90 to determine kx spacing (spectral case)
!------------------------------------------------------------------------------
subroutine geom_slab

  use grid,      only : n_s_grid, n_x_grid
  use constants, only : pi
  use general,   only : gkw_warn, gkw_exit
  use global,    only : r_tiny

  integer :: i, ix
  real    :: R02

  call gkw_warn('With slab geometry the parameters eps and q are redundant')

  if (abs(q - 1.E0) > r_tiny) then
    call gkw_exit('q must be set to unity when using slab geometry')
  end if

  if (abs(eps - 1.E0) > r_tiny) then
    call gkw_exit('eps must be set to unity when using slab geometry')
  end if

  ! Magnetic field strength (same for all radial modes)
  ! and poloidal angle and s      
  do ix = 1,n_x_grid
     do i = 1, n_s_grid
        bn_G(ix,i) = 1.E0
        pol_angle(ix,i)=2E0*pi*sgr(i)
        
        !Set the toroidal fraction of the field
        !bt_frac should always be positive
        !Whats this do? It projects the rotation from toroidal to parallel
        !Used in source uprim term and in calc of momentum flux
        bt_frac(ix,i)=1.
     end do
  end do

  ! set the minimum and maximum magnetic field strength 
  bmin = 1.E0
  bmax = 1.E0

  !There is no choice to be made here any more 
  R0   = 1.E0
  R02  = R0*R0
  lfun = 0.E0
  !These are set to zero just in case.
  do ix = 1,n_x_grid
     do i = 1, n_s_grid
        Rfun(ix,i)=1.E0
     end do
  end do

  ! Calculate the normfactor for k_zeta 
  kthnorm = 1.E0
  kxnorm = 1.E0

  ! metric elements
  ! Normalised: metric(i,2,2)=g_zeta_zeta*Rref**2 -> What is Rref?
  do ix = 1,n_x_grid
     do i = 1, n_s_grid 

        !In slab 
        
        ! the x' x' element 
        metric_G(ix,i,1,1) = 1.E0 
        
        ! the x' y' element 
        metric_G(ix,i,1,2) = shat*sgr(i)
        
        ! the x' s element 
        metric_G(ix,i,1,3) = 0. 
        
        ! the y' y' element 
        metric_G(ix,i,2,2) = 1 + (sgr(i)*shat)**2

        ! the y' s element 
        metric_G(ix,i,2,3) = 0.
        
        ! the s s element 
        metric_G(ix,i,3,3) = 1.E0
        
        ! for the other elements symmetry applies 
        metric_G(ix,i,2,1) = metric_G(ix,i,1,2) 
        metric_G(ix,i,3,1) = metric_G(ix,i,1,3) 
        metric_G(ix,i,3,2) = metric_G(ix,i,2,3) 
        
     end do
  end do
  
  dpfdpsi(:)   = 1./(2.*pi)
  bups(:)      = real(signJ) 
  dBdpsi(:,:)  = 0.0 
  dBds(:,:)    = 0.0
  dRdpsi(:,:)  = 0.0
  dRds(:,:)    = 0.0
  dZdpsi(:,:)  = 0.0
  dZds(:,:)    = 0.0

end subroutine geom_slab


!------------------------------------------------------------------------------
!> This routine performs the transformation necessary for the shifted 
!> metric. An alpha' is calculated that makes the g^zeta-psi element 
!> zero. Mathematically alpha is a function only of the radius. 
!> This routine must always be called for all nonspectral cases
!> since it also calculates the parallel boundary connection
!> This routine currently has two functions which should be split as follows
!>   geom_nonspectral_setup: the lines before the second return
!>   geom_shift_metric:      the lines after the second return
!------------------------------------------------------------------------------

subroutine geom_shift_metric 
 
  use control, only : spectral_radius, flux_tube, shift_metric
  use grid,    only : n_x_grid, n_s_grid

  integer :: ix, i
  real    :: dum1, dum2, ix_mid

  ! At present the metric is only shifted if spectral_radius is .false. 
  if (spectral_radius) return  
  
  ! Initialization values 
  do i = 1, n_s_grid    !initialization of alphak_xbnd and alphakp
    alphak_xbnd(i) = 0.E0
    do ix = 1, n_x_grid
      alphakp(ix,i) = 0.E0
    end do
  end do 
  
  ix_mid = real(n_x_grid + 1) / 2.
      
  ! Deterimine the shift end grid (the same for shift_metric = T and F) 
  do ix = 1, n_x_grid 
  
    if (flux_tube) then 
      shift_end_grid(ix) = shatx(ix)*qx(ix)/xgr(ix) * dxgr * (real(ix)-ix_mid)
    else   
      ! NOTE the shift must still be normalized to rhostar for the global case 
      ! this is done in init !!! 
      shift_end_grid(ix) = qx(ix)  
    endif 
  
  end do

  ! return if the metric is not shifted
  if (.not. shift_metric) return
  
  ! shift the metric
  do ix = 1, n_x_grid; do i = 1, n_s_grid 
  
    alphakp(ix,i) = metric_G(ix,i,2,1) / metric_G(ix,i,1,1) 
  
    metric_G(ix,i,2,2) = metric_G(ix,i,2,2) - 2*alphakp(ix,i)*          &
                       & metric_G(ix,i,2,1) + alphakp(ix,i)**2 *        &
                       & metric_G(ix,i,1,1) 
    metric_G(ix,i,2,1) = 0.E0 
    metric_G(ix,i,2,3) = metric_G(ix,i,2,3) - alphakp(ix,i)*            &
                       & metric_G(ix,i,1,3) 
  
    metric_G(ix,i,1,2) = metric_G(ix,i,2,1) 
    metric_G(ix,i,3,2) = metric_G(ix,i,2,3) 
  
  end do; end do 
  
  ! Integrate the alpha^prime towards radius (trapezium rule) 
  do i = 1, n_s_grid 
  
    dum1 = alphakp(1,i) 
  
    alphakp(1,i) = 0.5*dxgr*alphakp(1,i) 
    do ix = 2, n_x_grid 
      dum2 = alphakp(ix,i)
  
      alphakp(ix,i) = alphakp(ix-1,i) + dxgr * 0.5 * (dum2+dum1)
      dum1 = dum2
    end do 
  
    alphak_xbnd(i) = alphakp(n_x_grid,i) + 0.5 * dxgr * dum1 
  
  end do 
  
  ! alphakp array now contains alphak !
  ! choose the middle of the box to have zero alpha    
  do i = 1, n_s_grid
    ! this works for both odd and even n_x_grid
    dum1 = (alphakp((n_x_grid+1)/2,i)+alphakp((n_x_grid+2)/2,i))/2.

    do ix = 1, n_x_grid 
      alphakp(ix,i) = alphakp(ix,i) - dum1 
    end do 
  end do 
  
end subroutine geom_shift_metric


!----------------------------------------------------------------------------
!>
!----------------------------------------------------------------------------
subroutine geom_chease()
  use grid, only : n_s_grid, nperiod,n_x_grid
  use io, only : get_free_file_unit
  use constants, only : pi
  use mpiinterface, only : root_processor, mpibcast
  use general, only : gkw_abort, gkw_warn, gkw_exit
  use control, only : spectral_radius, flux_tube
  use control, only : shift_metric
  integer, allocatable :: s_indx(:)  ! indexes for the s_grid
  real, allocatable    :: R_FS(:), dzetadpf(:), dzetadchi(:)
  real, allocatable    :: dRdpf(:), dZdpf(:)
  real, allocatable    :: pf_dum(:)  ! 1D (pf) arrays
  real, allocatable    :: pf_s_dum(:,:)  ! 2D (pf,s) arrays

  ! variables specific to geom_type='chease'
  ! radial coordinate in chease is pf, the poloidal flux
  integer :: igeom, ierr, ns_c, npf_c, npol, s_coeff
  real    :: dpsidpf, p, dpdpsi, jac, dqdpf, dqdpf_dum, F, Raxis
  real    :: g_zz_LFS, g_pp_LFS, R_LFS, dRdpf_LFS
  integer :: pf_indx = 0 ! index for the pf-grid
  character (len=20) tdum
  integer :: ios=0

  integer :: i,ix,j,isa
  real :: dum
  
  ! WARNING - before changing the chease geom to the new system
  ! add test cases with signB/J = -1 for an up-down asymmetric case, and also with mode_box.

  ! chease /slab calculations are done for one radial grid point only
  ix = 1

  !number of points per poloidal turn for the GKW s-grid
  npol = n_s_grid/(2*nperiod-1)

  ! the CHEASE file is read by root processor only

  if (root_processor) then 
    ! 1) OPEN
    call get_free_file_unit(igeom)
    open(igeom, status='old', action='read', file=eqfile)
    rewind igeom

    ! 2) READ dimensions 
    ! number of points for pf-grid and s-grid
    read(igeom,*,IOSTAT=ios) tdum, npf_c, tdum, ns_c
    call check_read(ios)
    ! reference R and B (used for normalisation in GKW)
    ! taken to be R0EXP and B0EXP used for normalisation in CHEASE 
    ! R0EXP and B0EXP close to the magnetic axis values (but not exactly) 
    read(igeom,*,IOSTAT=ios) tdum, Rref, tdum, Bref, tdum, Raxis
    call check_read(ios)
  endif ! root processor

  call mpibcast(npf_c,1)
  call mpibcast(ns_c,1)
  call mpibcast(Rref,1)
  call mpibcast(Bref,1)
  call mpibcast(Raxis,1)

  ! 3) more checks
  ! check that CHEASE and GKW s-grid are compatible and that
  ! CHEASE grid is at least two times denser than GKW grid
  if (mod(ns_c,2*npol) /= 0) then
    call gkw_exit( &
       &'NCHI in chease has to be a multiple of 2*N_s_grid/(2*NPERIOD-1)')
  end if
  s_coeff=ns_c/npol ! how much the s-grid is denser in chease

  ! 4) ALLOCATE the arrays
  allocate(s_indx(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate s_indx in geom')
  allocate(R_FS(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate R_FS in geom')
  allocate(dRdpf(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dRdpf in geom')
  allocate(dZdpf(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dZdpf in geom')
  allocate(dzetadpf(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dzetadpf in geom')
  allocate(dzetadchi(1:n_s_grid),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate dzetadchi in geom')
  allocate(pf_dum(1:npf_c),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate pf_dum in geom')
  allocate(pf_s_dum(1:npf_c,1:ns_c),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate pf_s_dum in geom')

  if (root_processor) then 

    ! 5) READ 1D arrays

    ! find the index for the radial grid
    select case (eps_type)
    case (1) ! use eps to select FS

      ! pf-grid, s-grid, Rgeom (not used)
      read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,npf_c) !poloidal flux grid
      call check_read(ios)
      read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,ns_c) !s 
      call check_read(ios)
      read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,npf_c) !Rgeom
      call check_read(ios)

      ! amin=(Rmax-Rmin)/2 
      read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) !amin
      call check_read(ios)

      ! find the index for the pf grid
      pf_indx=1
      do while ( abs(pf_dum(pf_indx+1)/Rref-eps)<abs(pf_dum(pf_indx)/Rref-eps) &
         & .and. pf_indx.LT.npf_c-2 ) 
        pf_indx = pf_indx+1
      end do
      if (abs(pf_dum(pf_indx)/Rref-eps) .GE. 0.002) then
        call gkw_exit('Radial grid too coarse in CHEASE: increase NPSI')
      endif

    case (2)! use rho_pf to select FS (rho_pf is the sqrt of the poloidal flux)

      ! pf-grid
      read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) !poloidal flux grid
      call check_read(ios)

      ! find the index for the pf grid
      pf_indx=1
      do while (abs(sqrt(pf_dum(pf_indx+1)/pf_dum(npf_c))-eps) < &
         & abs(sqrt(pf_dum(pf_indx)/pf_dum(npf_c))-eps) &
         & .and. pf_indx < npf_c-2)
        pf_indx = pf_indx+1
      end do
      if (abs(sqrt(pf_dum(pf_indx)/pf_dum(npf_c))-eps) .GE. 0.008) then
        call gkw_exit('Radial grid too coarse in CHEASE: increase NPSI')
      endif

      ! s-grid, Rgeom (not used)
      read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,ns_c) !s 
      call check_read(ios)
      read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,npf_c) !Rgeom
      call check_read(ios)

      ! amin=(Rmax-Rmin)/2 
      read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) !amin
      call check_read(ios)

    case default
      call gkw_exit('geom: value not allowed for eps_type')
    end select

    ! value of eps=amin/Rref (used in mode.F90 for kxspace)
    eps = pf_dum(pf_indx)/Rref

    ! dpsidpf 
    ! -> to go from chease (pf,zeta,s) to GKW (psi,zeta,s)
    ! Reminder: pf = poloidal flux,  psi = eps
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! damindpf
    call check_read(ios)
    dpsidpf = pf_dum(pf_indx)/Rref

    ! d2amindpf2 (not yet used) 
    ! 1 / second radial derivative of poloidal flux
    read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,npf_c) ! d2amindpf2
    call check_read(ios)

    ! Bmax, Bmin, q, shat
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! Bmax
    call check_read(ios)
    bmax = pf_dum(pf_indx)/Bref
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! Bmin
    call check_read(ios)
    bmin = pf_dum(pf_indx)/Bref
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! q
    call check_read(ios)
    q = pf_dum(pf_indx)
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! dqdpf
    call check_read(ios)
    dqdpf = pf_dum(pf_indx)
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! dqdpf_chk
    call check_read(ios)
    dqdpf_dum = pf_dum(pf_indx)
    shat = eps * dqdpf_dum / q / dpsidpf

    ! p and dpdpsi (not normalised) 
    ! -> to be used for beta and betaprime
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! p
    call check_read(ios)
    p = pf_dum(pf_indx)
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! dpdpf
    call check_read(ios)
    dpdpsi = pf_dum(pf_indx)/dpsidpf

    ! jacobian J_pf_zeta_s
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c) ! jacobian
    call check_read(ios)
    jac=pf_dum(pf_indx)

    ! not used (don't change the order!)
    read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,npf_c) ! djacdpf
    call check_read(ios)

    !chease F=RB_t>0 (not the same as GKW F tensor which is ffun)
    read(igeom,*,IOSTAT=ios) tdum, (pf_dum(i),i=1,npf_c)
    call check_read(ios)
    F=pf_dum(pf_indx)

    !not used
    read(igeom,*,IOSTAT=ios) tdum, (dum,i=1,npf_c) ! dFdpf
    call check_read(ios)

    ! 6) READ 2D arrays
    ! build the index correspondance for the s-grid
    ! CHEASE s-array is sgr_c(j)= (j-1) / ns_c 
    ! with ns_c = NCHI and sgr_c=0 at the LFS midplane (sgr_c increases counterclockwise)
    do i = 1, n_s_grid
      s_indx(i) = nint(modulo(sgr(i),1.)*ns_c + 1)
    end do

    ! metric_G elements, (pf,zeta,s) coordinates
    ! g_pf_pf
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) !g11
    call check_read(ios)
    do i = 1, n_s_grid
      metric_G(ix,i,1,1)=pf_s_dum(pf_indx,s_indx(i))
    end do
    g_pp_LFS=pf_s_dum(pf_indx,1)
    ! g_pf_zeta
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) !g12
    call check_read(ios)
    do i = 1, n_s_grid
      metric_G(ix,i,1,2)=pf_s_dum(pf_indx,s_indx(i))
      metric_G(ix,i,2,1)=pf_s_dum(pf_indx,s_indx(i))
    end do
    ! g_pf_s
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) !g13
    call check_read(ios)
    do i = 1, n_s_grid
      metric_G(ix,i,1,3)=pf_s_dum(pf_indx,s_indx(i))
      metric_G(ix,i,3,1)=pf_s_dum(pf_indx,s_indx(i))
    end do
    ! g_zeta_zeta
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) !g22
    call check_read(ios)
    do i = 1, n_s_grid
      metric_G(ix,i,2,2)=pf_s_dum(pf_indx,s_indx(i))
    end do
    g_zz_LFS=pf_s_dum(pf_indx,1)

    ! g_zeta_s
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) !g23
    call check_read(ios)
    do i = 1, n_s_grid
      metric_G(ix,i,2,3)=pf_s_dum(pf_indx,s_indx(i))
      metric_G(ix,i,3,2)=pf_s_dum(pf_indx,s_indx(i))
    end do
    ! g_s_s
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) !g33
    call check_read(ios)
    do i = 1, n_s_grid
      metric_G(ix,i,3,3)=pf_s_dum(pf_indx,s_indx(i))
    end do

    ! magnetic field
    ! norm(B)
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) !B
    call check_read(ios)
    do i = 1, n_s_grid
      bn_G(ix,i)=pf_s_dum(pf_indx,s_indx(i))/Bref
    end do
    ! dBdpsi (not normalized)
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dBdpf
    call check_read(ios)
    do i = 1, n_s_grid
      dBdpsi(ix,i)=pf_s_dum(pf_indx,s_indx(i))/dpsidpf
      !Need dBdpsi_LFS for ICRH to be exactl correct
    end do
    ! dBds (not normalized)
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dBds
    call check_read(ios)
    do i = 1, n_s_grid
      dBds(ix,i)=pf_s_dum(pf_indx,s_indx(i))
    end do

    ! R and Z
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! R
    call check_read(ios)
    do i = 1, n_s_grid
      R_FS(i)=pf_s_dum(pf_indx,s_indx(i)) ! in [m], divided later by Rref (rfun=R_FS/Rref written in geom.dat)
    end do
    R_LFS=pf_s_dum(pf_indx,1) !value of R at the LFS and Z=Zaxis
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! Z
    call check_read(ios)
    ! dimensionless, not used in the code so far (only written in geom.dat)
    do i = 1, n_s_grid
      Z_FS(ix,i)=pf_s_dum(pf_indx,s_indx(i))/Rref
    end do

    ! R and Z gradients
    ! dRdpf
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dRdpf
    call check_read(ios)
    do i = 1, n_s_grid
      dRdpf(i)=pf_s_dum(pf_indx,s_indx(i))
    end do
    dRdpf_LFS=pf_s_dum(pf_indx,1)
    ! dRds
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dRds
    call check_read(ios)
    do i = 1, n_s_grid
      dRds(ix,i)=pf_s_dum(pf_indx,s_indx(i))
    end do
    ! dZdpf
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dZdpf
    call check_read(ios)
    do i = 1, n_s_grid
      dZdpf(i)=pf_s_dum(pf_indx,s_indx(i))
    end do
    ! dZds
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dZds
    call check_read(ios)
    do i = 1, n_s_grid
      dZds(ix,i)=pf_s_dum(pf_indx,s_indx(i))
    end do

    ! poloidal angle
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! theta
    call check_read(ios)
    do i = 1, n_s_grid
      dum=pf_s_dum(pf_indx,s_indx(i))
      if (dum.gt.pi) then
        dum = dum - 2E0*pi
      endif
      pol_angle(ix,i) = dum + 2E0*pi*(floor(real(i-1)/real(npol))-real(nperiod)+1)
    end do

    ! not used
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((dum,i=1, npf_c),j=1, ns_c) ! Ah
    call check_read(ios)
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((dum,i=1, npf_c),j=1, ns_c) ! dAhdpf
    call check_read(ios)

    ! dzetadpf
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dzetadpf
    call check_read(ios)
    do i = 1, n_s_grid
      dzetadpf(i)=pf_s_dum(pf_indx,s_indx(i))
    end do

    ! dzetadchi
    read(igeom,'(A)',IOSTAT=ios) tdum
    call check_read(ios)
    read(igeom,'(5ES20.10)',IOSTAT=ios) ((pf_s_dum(i,j),i=1, npf_c),j=1, ns_c) ! dzetadchi
    call check_read(ios)
    do i = 1, n_s_grid
      dzetadchi(i)=pf_s_dum(pf_indx,s_indx(i))
    end do

    close(igeom)

  endif ! root processor

  call mpibcast(eps,        1)  
  call mpibcast(dpsidpf,    1)  
  call mpibcast(bmax,       n_x_grid)  
  call mpibcast(bmin,       n_x_grid)  
  call mpibcast(q,          1)  
  call mpibcast(dqdpf,     1)  
  call mpibcast(dqdpf_dum, 1)  
  call mpibcast(shat,       1)  
  call mpibcast(p,          1)  
  call mpibcast(dpdpsi,     1)  
  call mpibcast(jac,        1)  
  call mpibcast(F,          1)
  ! just send a part, as most of the elements are undefined
  call mpibcast(metric_G(1,:,:,:),1*n_s_grid*3*3)
  call mpibcast(bn_G(1,:),       1*n_s_grid)
  call mpibcast(dBdpsi(1,:),     1*n_s_grid)
  call mpibcast(dBds(1,:),       1*n_s_grid)
  call mpibcast(Z_FS(1,:),       1*n_s_grid)
  call mpibcast(dRds(1,:),       1*n_s_grid)
  call mpibcast(dZds(1,:),       1*n_s_grid)
  call mpibcast(pol_angle(1,:),  1*n_s_grid)
  call mpibcast(g_zz_LFS,   1)
  call mpibcast(g_pp_LFS,   1)
  call mpibcast(R_FS,       n_s_grid)
  call mpibcast(R_LFS,      1)
  call mpibcast(dRdpf,      n_s_grid)
  call mpibcast(dRdpf_LFS,  1)
  call mpibcast(dZdpf,      n_s_grid)
  call mpibcast(dzetadpf,   n_s_grid)
  call mpibcast(dzetadchi,  n_s_grid)

  ! 7) Compute the missing elements

  ! Correction to the metric_G elements involving dzetadpf
  ! because dzetadpf is not periodic in s:
  ! dzetadpf(1) = dqdpf_dum 
  ! dzetadpf(s) = dzetadpf(s_0) + dum * dqdpf
  !    with  s=s_0+ dum  and  0<=s_0<1 
  do i = 1, n_s_grid 
    dum = int(i/npol) - nperiod
    if (mod(i,npol).GT.real(npol/2)) then
      dum = dum + 1 
    endif
    metric_G(ix,i,1,2) = metric_G(ix,i,1,2) + &
       & dum * dqdpf_dum * metric_G(ix,i,1,1)
    metric_G(ix,i,2,1) = metric_G(ix,i,2,1) + &
       & dum * dqdpf_dum * metric_G(ix,i,1,1)
    metric_G(ix,i,2,3) = metric_G(ix,i,2,3) + &
       & dum * dqdpf_dum * metric_G(ix,i,1,3)
    metric_G(ix,i,3,2) = metric_G(ix,i,3,2) + &
       & dum * dqdpf_dum * metric_G(ix,i,1,3)
    metric_G(ix,i,2,2) = metric_G(ix,i,2,2) + (dum * dqdpf_dum)**2 * metric_G(ix,i,1,1) + &
       & 2 * dum * dqdpf_dum * dzetadpf(i) * metric_G(ix,i,1,1) + & 
       & 4 * pi * dum * dqdpf_dum * dzetadchi(i) * metric_G(ix,i,1,3)
  end do

  ! beta and beta prime (does not have the bn_G dependence)
  beta_eq = 2 * 1.25663706144E-6 * p / Bref**2
  betaprime_eq = 2 * 1.25663706144E-6 * dpdpsi / Bref**2

  ! calculate the array for the parallel derivative (ffun)
  do i = 1, n_s_grid
    ffun(ix,i) = 2*pi*signJ*Rref/bn_G(ix,i)/Bref/jac
  end do

  ! the function connected with the trapping terms (gfun)
  do i = 1, n_s_grid 
    gfun(ix,i) = ffun(ix,i)*dBds(ix,i)/bn_G(ix,i)/Bref
  end do

  ! calculate the array connected with the ExB velocity (efun)
  do i = 1, n_s_grid 

    ! the diagonal components are zero 
    efun(ix,i,1,1) = 0.
    efun(ix,i,2,2) = 0. 
    efun(ix,i,3,3) = 0. 

    ! the psi zeta component 
    efun(ix,i,1,2) =  signJ*pi*Rref**2/bn_G(ix,i)**2/Bref*  &
       & (metric_G(ix,i,1,1)*metric_G(ix,i,2,2) -  &
       &  metric_G(ix,i,1,2)**2)*dpsidpf

    ! the psi s component 
    efun(ix,i,1,3) =  signJ*pi*Rref**2/bn_G(ix,i)**2/Bref*  &
       & (metric_G(ix,i,1,1)*metric_G(ix,i,2,3) -  &
       &  metric_G(ix,i,1,2)*metric_G(ix,i,1,3))*dpsidpf*signB*signJ

    ! the zeta s component 
    efun(ix,i,2,3) =  signJ* pi*Rref**2/bn_G(ix,i)**2/Bref* &
       & (metric_G(ix,i,1,2)*metric_G(ix,i,2,3) -  &
       &  metric_G(ix,i,2,2)*metric_G(ix,i,1,3))

    ! the other components are anti-symmetric_G 
    efun(ix,i,2,1) = - efun(ix,i,1,2) 
    efun(ix,i,3,1) = - efun(ix,i,1,3) 
    efun(ix,i,3,2) = - efun(ix,i,2,3) 

  end do

  ! calculate the curvature function  (dfun)
  do i = 1, n_s_grid 

    ! the psi component
    dfun(ix,i,1) = -2/bn_G(ix,i)/Bref*efun(ix,i,1,3)*dBds(ix,i)

    ! the zeta component 
    dfun(ix,i,2) = -2/bn_G(ix,i)/Bref*(efun(ix,i,2,1)*dBdpsi(ix,i) + &
       &  efun(ix,i,2,3)*dBds(ix,i))

    ! the s component 
    dfun(ix,i,3) = -2/bn_G(ix,i)/Bref*efun(ix,i,3,1)*dBdpsi(ix,i)

  end do

  ! calculate the array connected with the coriolis drift (hfun)
  ! Omega/vthref=Cte assumed (rigid body)
  ! As VCOR=vtor(R=Rref)/vthref, Omega/vthref=VCOR/Rref
  ! term VCOR=vtor(R=Rref)/vthref not included here (added in linear_terms)
  do i = 1, n_s_grid 

    ! the psi component
    hfun(ix,i,1) = - signB*Rref**2/bn_G(ix,i)*(dZdpf(i)*metric_G(ix,i,1,1) + &
       & dZds(ix,i)*metric_G(ix,i,3,1))*dpsidpf / Rref

    ! the zeta component 
    hfun(ix,i,2) = - signB*Rref**2/bn_G(ix,i)* &
       & (dZdpf(i)*metric_G(ix,i,1,2)*signB*signJ + &
       & dZds(ix,i)*metric_G(ix,i,3,2)*signB*signJ) / Rref

    ! the s component 
    hfun(ix,i,3) = - signB*Rref**2/bn_G(ix,i)*(dZdpf(i)*metric_G(ix,i,1,3) + &
       & dZds(ix,i)*metric_G(ix,i,3,3) - dZds(ix,i)*ffun(ix,i)**2/Rref**2) &
       & / Rref

  end do

  ! Centrifugal trapping
  ! R0 (in Peeters PoP09 on centrifugal force) 
  ! specifies location of density gradient definitions
  select case(R0_loc)
  case('axis') ! axis=place on flux surface at which R=Raxis
    !R0 = Raxis/Rref

    ! Some s_point (isa) nearest to the axis must be found. 
    ! Taken as 'top' intersection of the axis with the flux surface
    isa=minloc(abs(pol_angle(ix,1:n_s_grid)-pi/2.),1)
    R0=R_FS(isa)/Rref

    lfun = 2.E0*dRdpf(isa)*R_FS(isa)/(dpsidpf*Rref**2)
    call gkw_warn('Geom: R0_loc=axis option not exact for CHEASE')
    call gkw_warn('For heavy impurities densities may be inaccurate')
    if (root_processor) write(*,'(A,i3,A,f9.4)') ' s point ', isa,  &
       & ' used with poloidal angle=', pol_angle(1,isa)

  case('LFS')  !Low field side in the plane of magnetic axis
    R0 = R_LFS/Rref
    lfun = 2.E0*dRdpf_LFS*R_LFS/(dpsidpf*Rref**2)     
  case default
    call gkw_exit('Geom: unknown R0_loc option. &
       & Allowed are: axis, LFS')
  end select

  ! calculate the arrays connected with centrifugal drift (ifun)
  ! does not include the (Rref*Omega/vthref)**2 term
  do i = 1, n_s_grid 

    ! the psi component
    ifun(ix,i,1) = 2*R_FS(i)*(efun(ix,i,1,3)*dRds(ix,i))/ Rref**2

    ! the zeta component 
    ifun(ix,i,2) = 2*R_FS(i)*(efun(ix,i,2,1)*dRdpf(i)/dpsidpf+ &
       & efun(ix,i,2,3)*dRds(ix,i))/ Rref**2
    ! the s component 
    ifun(ix,i,3) = 2*R_FS(i)*(efun(ix,i,3,1)*dRdpf(i)/dpsidpf)/ Rref**2

    !Centrifugal trapping 
    ! R_N^2-R0_N^2
    jfun(ix,i) = R_FS(i)*R_FS(i)/Rref**2 - R0*R0

    ! djfun/dpsi (FJC implemented)
    kfun(ix,i)=2.E0*dRdpf(i)*R_FS(i)/(dpsidpf*Rref**2)
    kfun(ix,i)=kfun(ix,i)-lfun

    ! local major radius
    Rfun(ix,i)=R_FS(i)/Rref

  end do

  ! not exact if no point at LFS
  ! YC: to update
  jfunl=maxval(jfun(ix,1:n_s_grid))
  jfunh=minval(jfun(ix,1:n_s_grid))

  ! Calculate the normfactor for k_zeta 
  kthnorm = Rref * sqrt(g_zz_LFS) 
  ! Calculate the normfactor for k_pf 
  kxnorm = Rref * sqrt(g_pp_LFS) * dpsidpf


  ! metric_G elements for the (psi,zeta,s) coordinates
  ! all multiplied by Rref**2 for normalisation
  do i=  1, n_s_grid
    metric_G(ix,i,1,1) = metric_G(ix,i,1,1) * dpsidpf**2 * Rref**2
    metric_G(ix,i,1,2) = signB*signJ*metric_G(ix,i,1,2) * dpsidpf * Rref**2
    metric_G(ix,i,2,1) = signB*signJ*metric_G(ix,i,2,1) * dpsidpf * Rref**2
    metric_G(ix,i,1,3) = metric_G(ix,i,1,3) * dpsidpf * Rref**2
    metric_G(ix,i,3,1) = metric_G(ix,i,3,1) * dpsidpf * Rref**2
    metric_G(ix,i,2,2) = metric_G(ix,i,2,2) * Rref**2
    metric_G(ix,i,2,3) = signB*signJ*metric_G(ix,i,2,3) * Rref**2
    metric_G(ix,i,3,2) = signB*signJ*metric_G(ix,i,3,2) * Rref**2
    metric_G(ix,i,3,3) = metric_G(ix,i,3,3) * Rref**2
  end do

  !Compute the toroidal fraction of the field
  !and derivatives of R
  do i=  1, n_s_grid
    !bt_frac should always be positive
    bt_frac(ix,i)=F/(R_FS(i)*Bref*bn_G(ix,i)) 
    !F=RB_t 
    !Rref normalises R_FS everywhere else, but here R_FS cancels with R in F.

    !Some checks
    if(bt_frac(ix,i)<0..or.bt_frac(ix,i)>1.0) then
      call gkw_abort('error in bt_frac')
    end if

  end do

  ! Temporary hacks to allow chease geom to work with (some) nonspectral cases
  if (.not. shift_metric .and. flux_tube) then 
    call distribute_geom_flux_tube(2)   

    ! cheat for the jacobian for now
    ! select an arbitrary i point 
    i = 1       
    bups(:)=ffun(:,i)*bn_G(:,i) 
    do ix = 1, n_x_grid
      jacobian_G(ix) = 1.E0 / (2.E0*abs(bups(ix)*efun(1,i,1,2)))               
    end do

    ! alternatively ??      
    ! jacobian_G(:) = jac/dpsidpf

    ! need to setup alphak_xbnd
    call geom_shift_metric
    if (.not. spectral_radius) then
      call gkw_warn('Chease geom is not verified for nonspectral!')
    end if

  else
    call gkw_abort('chease geom cannot yet be used for this nonspectral run')
  end if

  call geom_check_params(2)

end subroutine geom_chease


!------------------------------------------------------------------------------
!> This subroutine calculates the tensors needed by GKW from the various 
!> geometry quantities. Before calling this routine the following quantities 
!> must be specified 
!> 
!> 1. The magnetic field strength as a function of radius and coordinate 
!>    along the field line [bn_G(ix,i)]
!> 2. B^s the contra-variant component of the magnetic field. [bups(ix)]
!> 3. The metric tensor [metric_G(ix,i,3,3)] 
!> 4. The derivative of the magnetic field towards the radial coordinate 
!>    [dBdpsi(ix,i)] and towards the coordinate along the magnetic field 
!>    [dBds(ix,i)]  
!> 5. The value of R [Rfun(ix,i)] and the derivative of R [dRdpsi,dRds] and 
!>    Z [dZdpsi,dZds]. NOTE: Rfun will be reset if no finite epsilon 
!>    effects are taken into acount to the value 1.0   
!> 6. The first derivative of the poloidal flux towards the radial 
!>    coordinate [dpfdpsi(ix)] 
!> 7. The sign of the magnetic field and plasma current [sign_b, sign_j]
!> 8. The location of the centrifugal reference radius (R0) and the 
!>    derivative of its square towards the radial coordinate (lfun) 
!>
!> Other quantities that are not used by this routine (but elsewhere):
!> and must also be calculated for each geometry:
!> 9. bt_frac, pol_angle, jfunh, jfunl, bmax, bmin, kthnorm, kxnorm
!> 
!> Input to this routine: 
!> finite_epsilon : logical that determines if finite epsilon effects are 
!>                  retained. This is to be set .false. in the case of the 
!>                  s-alpha geometry 
!> gfun_num:      : Logical that determines whether the gfun tensor is 
!>                  calculated through numerical differentiation. If it 
!>                  is, it is a little inconsistent with some other tensor 
!>                  quantities, but it is to perserve the relation with 
!>                  previous implementations.
!> opt_efun:        optional to ignore check on efun flux_function
!------------------------------------------------------------------------------

subroutine calc_geom_tensors(finite_epsilon, gfun_num, opt_efun) 

  use grid,         only : n_x_grid, n_s_grid
  use constants,    only : pi
  use general,      only : gkw_abort
  use mpiinterface, only : root_processor

  logical, intent(in)           :: finite_epsilon, gfun_num
  logical, intent(in), optional :: opt_efun

  integer :: ix, i, j
  real    :: rtest  
  logical :: test_efun = .true.  

  rtest=0.

  if (present(opt_efun)) test_efun = opt_efun
  
  ! the F tensor (derivative along the field line) 
  do ix = 1, n_x_grid; do i = 1, n_s_grid 
    ffun(ix,i) = bups(ix) 
    if (finite_epsilon) ffun(ix,i) = ffun(ix,i) / bn_G(ix,i) 
  end do; end do 

  ! The G tensor (related to trapping) 
  if (gfun_num) then   
    call logbderiv    
  else 
    do ix = 1, n_x_grid; do i = 1, n_s_grid 
      gfun(ix,i) = ffun(ix,i) * dBds(ix,i) / bn_G(ix,i)  
    end do; end do 
  endif 

  ! the E tensor (connected with the ExB)
  do ix = 1, n_x_grid; 
  
    do i = 1, n_s_grid 

      ! anti-symmetric 
      efun(ix,i,1,1) = 0. 
      efun(ix,i,2,2) = 0. 
      efun(ix,i,3,3) = 0.

      efun(ix,i,1,2) = metric_G(ix,i,1,1)*metric_G(ix,i,2,2) - & 
                       metric_G(ix,i,2,1)*metric_G(ix,i,1,2) 
      efun(ix,i,1,3) = metric_G(ix,i,1,1)*metric_G(ix,i,2,3) - & 
                       metric_G(ix,i,2,1)*metric_G(ix,i,1,3) 
      efun(ix,i,2,3) = metric_G(ix,i,1,2)*metric_G(ix,i,2,3) - & 
                       metric_G(ix,i,2,2)*metric_G(ix,i,1,3) 

      ! anti-symmetric 
      efun(ix,i,2,1) = - efun(ix,i,1,2) 
      efun(ix,i,3,1) = - efun(ix,i,1,3)
      efun(ix,i,3,2) = - efun(ix,i,2,3) 

      ! common factors 
      efun(ix,i,:,:) = signJ*pi*dpfdpsi(ix)*efun(ix,i,:,:)
      if (finite_epsilon) efun(ix,i,:,:) = efun(ix,i,:,:) / bn_G(ix,i)**2
 
    end do 

    ! test on numerical equilibria (efun_eps_zeta should be a flux function)
    rtest = max(rtest,maxval(efun(ix,:,1,2)) - minval(efun(ix,:,1,2)))

    ! Take the mean value of efun_eps_zeta for all points
    ! Does it help numerical stability ?
    ! efun(ix,:,1,2)=sum(efun(ix,:,1,2))/n_s_grid
    ! efun(ix,:,2,1)=-efun(ix,:,1,2)
    
  end do

  if (root_processor) write(*,*)
  if (root_processor) write(*,'(" Geometry numerical accuracy: ",es10.3)') rtest
  if (root_processor) write(*,*)
  if (abs(rtest)>5e-4 .and. test_efun) then      
    if (root_processor) write(*,*) 'check geometry inputs ' // & 
                       & 'or increase points in Miller integrals' // &
                       & '(integer N in geom_miller in geom.f90)'
    call gkw_abort('Geom: Numerical equilibrium not converged')
  end if
  
  ! The D tensor (connected with the drift)
  do ix = 1, n_x_grid; do i = 1, n_s_grid; do j = 1, 3 
    dfun(ix,i,j)=-2.0*efun(ix,i,j,1)*dBdpsi(ix,i)-2.0*efun(ix,i,j,3)*dBds(ix,i) 
  end do; end do; end do 
  if (finite_epsilon) then 
    do ix = 1, n_x_grid; do i = 1, n_s_grid; do j = 1, 3 
      dfun(ix,i,j) = dfun(ix,i,j) / bn_G(ix,i) 
    end do; end do; end do 
  endif 

  ! the Hfun tensor (connected with the Coriolis drift) 
  do ix = 1, n_x_grid; do i = 1, n_s_grid
    do j = 1, 3 
      hfun(ix,i,j) = - signB * ( metric_G(ix,i,j,1) * dZdpsi(ix,i) +  &
                   &             metric_G(ix,i,j,3) * dZds(ix,i)     ) 
    end do  
    if (finite_epsilon) then 
      hfun(ix,i,3) = hfun(ix,i,3) + signB * bups(ix)**2 * dZds(ix,i)  &
                   & / (bn_G(ix,i)**2)  
    endif 
  end do; end do 

  if (finite_epsilon) then 
    do ix = 1, n_x_grid; do i = 1, n_s_grid; do j = 1, 3 
      hfun(ix,i,j) = hfun(ix,i,j) / bn_G(ix,i) 
    end do; end do; end do 
  endif 

  ! do the J-tensor first such that one can redfine Rfun for the case 
  ! where finite epsilon effects are neglected 
  ! Related to centrifugal trapping: jfun = R_N^2-R0_N^2, kfun= djfun/dpsi  
  do ix = 1, n_x_grid; do i = 1, n_s_grid
    if (finite_epsilon) then 
      jfun(ix,i) = rfun(ix,i)**2 - R0**2 
    else 
      jfun(ix,i) = 2.E0*(rfun(ix,i)-R0) 
    endif 
  end do; end do 

  ! Redefine rfun for the case in which finite epsilon effects are neglected
  if (.not. finite_epsilon) then 
    do ix = 1, n_x_grid; do i = 1, n_s_grid 
      rfun(ix,i) = 1.0 
    end do; end do 
  endif 

  ! The ifun tensor (connected with the centrifugal force) 
  do ix = 1, n_x_grid; do i = 1, n_s_grid
    do j = 1, 3 
      ifun(ix,i,j) = +  ( efun(ix,i,j,1) * dRdpsi(ix,i) +  &
                   &      efun(ix,i,j,3) * dRds(ix,i)        )
      ifun(ix,i,j) = 2.0E0 * rfun(ix,i) * ifun(ix,i,j) 
    end do 
  end do; end do 

  do ix = 1, n_x_grid; do i = 1, n_s_grid 
    kfun(ix,i) = 2.0*rfun(ix,i)*dRdpsi(ix,i) - lfun 
  end do; end do 

  ! select an arbitrary i point 
  i = 1 
  do ix = 1, n_x_grid
    jacobian_G(ix) = 1.E0 / (2.E0*abs(bups(ix)*efun(ix,i,1,2)))     
    !jacobian_G(ix) = 1.0E0           
  end do
  
end subroutine calc_geom_tensors 


!------------------------------------------------------------------------------
!> This subroutine simply copies local values into the arrays that have an 
!> dependence on the radial coordinate.  It is a temporary measure for the 
!> geometries that have not yet been transferred to the new format.
!------------------------------------------------------------------------------
subroutine distribute_geom_flux_tube(switch)

  use grid, only : n_x_grid  
  
  ! switch = 1 only for tensors needed before before calc_geom_tensors
  ! switch = 2 for all tensors
  integer, intent(in) :: switch

  integer :: ix

  do ix = 1, n_x_grid 
    qx(ix)             = q 
    shatx(ix)          = shat
    xgr(ix)            = eps  
    bups(ix)           = bups(1)
    bn_G(ix,:)         = bn_G(1,:) 
    metric_G(ix,:,:,:) = metric_G(1,:,:,:)
    pol_angle(ix,:)    = pol_angle(1,:)
    bt_frac(ix,:)      = bt_frac(1,:)
    rfun(ix,:)         = rfun(1,:)

    if (switch == 1) then
      dBdpsi(ix,:)     = dBdpsi(1,:)
      dBds(ix,:)       = dBds(1,:)
      dRdpsi(ix,:)     = dRdpsi(1,:)
      dRds(ix,:)       = dRds(1,:)
      dZdpsi(ix,:)     = dZdpsi(1,:)
      dZds(ix,:)       = dZds(1,:)
      dpfdpsi(ix)      = dpfdpsi(1)
    end if 

    if (switch < 2) cycle

    ! following two are set by parallelize_geom
    !metric(ix,:,:,:)   = metric(1,:,:,:)
    !bn(ix,:)           = bn(1,:)
    ffun(ix,:)         = ffun(1,:)
    gfun(ix,:)         = gfun(1,:)
    jfun(ix,:)         = jfun(1,:)    
    kfun(ix,:)         = kfun(1,:)
    dfun(ix,:,:)       = dfun(1,:,:)
    hfun(ix,:,:)       = hfun(1,:,:)
    ifun(ix,:,:)       = ifun(1,:,:)
    efun(ix,:,:,:)     = efun(1,:,:,:) 
  end do 

end subroutine distribute_geom_flux_tube


!------------------------------------------------------------------------------
!> This subroutine writes the equilibrium dependent quantities 
!> Must be called before parallelize_geom
!------------------------------------------------------------------------------
subroutine geom_output()

  use control,      only : flux_tube
  use grid,         only : n_s_grid, n_x_grid 
  use general,      only : gkw_abort
  use io,           only : output_array, ascii_fmt
  use mpiinterface, only : root_processor
 
  integer :: i, ix

  if (geom_parallelized) call gkw_abort('geom_output after parallelize') 

  if (.not.root_processor) return 
    
  ! for the flux tube none of the quantities depends on radius 
  if (flux_tube) then 

    ! scalar quantities
    call output_array('NS', 'geom', (/ n_s_grid*1.0 /), 'C', '(F8.0)', ascii_fmt)
    call output_array('Rref', 'geom', (/ Rref /), 'C', '(ES13.5)', ascii_fmt)
    call output_array('R0', 'geom', (/ R0/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('lfun', 'geom', (/ lfun/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('Bref', 'geom', (/ Bref/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('eps', 'geom', (/ eps/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('q', 'geom', (/ q/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('shat', 'geom', (/ shat/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('bmin', 'geom', (/ bmin(1)/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('bmax', 'geom', (/ bmax(1)/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('kthnorm', 'geom', (/ kthnorm/), 'C', '(ES13.5)', &
       & ascii_fmt) ! = Rref*sqrt(metric(isg_lfs,2,2))
    call output_array('krnorm', 'geom', (/ kxnorm/), 'C', '(ES13.5)', &
       & ascii_fmt) ! = Rref*sqrt(metric(isg_lfs,1,1))
    call output_array('beta_eq', 'geom', (/ beta_eq/), 'C', '(ES13.5)', &
       & ascii_fmt) ! = 2*mu0*p/Bref**2
    call output_array('betaprime_eq', 'geom', (/ betaprime_eq/), 'C', &
       & '(ES13.5)', ascii_fmt) ! = 2*mu0*dpdeps/Bref**2
    call output_array('jfunh', 'geom', (/ jfunh/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('jfunl', 'geom', (/ jfunl/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('Jacobian', 'geom', (/ clean0(Jacobian_G(1))/), 'C', &
       & '(5ES13.5)', ascii_fmt)

    ! s-dependent quantities
    call output_array('s_grid', 'geom', &
       & sgr(1:n_s_grid), 'C', '(5ES13.5)', ascii_fmt)
    !bn = B/Bref
    call output_array('bn', 'geom', &
       & bn_G(1,:), 'C', '(5ES13.5)', ascii_fmt) ! = B/Bref
    call output_array('poloidal_angle', 'geom', &
       & clean0(pol_angle(1,:)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_eps_eps', 'geom', &
       & clean0(metric_G(1,:,1,1)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_eps_zeta', 'geom', &
       & clean0(metric_G(1,:,1,2)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_eps_s', 'geom', &
       & clean0(metric_G(1,:,1,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_zeta_zeta', 'geom', &
       & clean0(metric_G(1,:,2,2)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_zeta_s', 'geom', &
       & clean0(metric_G(1,:,2,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_s_s', 'geom', &
       & clean0(metric_G(1,:,3,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('D_eps', 'geom', &
       & clean0(dfun(1,:,1)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('D_zeta', 'geom', &
       & clean0(dfun(1,:,2)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('D_s', 'geom', &
       & clean0(dfun(1,:,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('E_eps_zeta', 'geom', &
       & clean0(efun(1,:,1,2)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('E_eps_s', 'geom', &
       & clean0(efun(1,:,1,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('E_zeta_s', 'geom', &
       & clean0(efun(1,:,2,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('F', 'geom', &
       & clean0(ffun(1,:)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('G', 'geom', &
       & clean0(gfun(1,:)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('H_eps', 'geom', &
       & clean0(hfun(1,:,1)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('H_zeta', 'geom', &
       & clean0(hfun(1,:,2)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('H_s', 'geom', &
       & clean0(hfun(1,:,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('I_eps', 'geom', &
       & clean0(ifun(1,:,1)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('I_zeta', 'geom', &
       & clean0(ifun(1,:,2)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('I_s', 'geom', &
       & clean0(ifun(1,:,3)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('J', 'geom', &
       & clean0(jfun(1,:)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('K', 'geom', &
       & clean0(kfun(1,:)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('Bt_frac', 'geom', &
       & clean0(Bt_frac(1,:)), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('R', 'geom', &
       & clean0(Rfun(1,:)), 'C', '(5ES13.5)', ascii_fmt)
    !temporarily removed until test cases are updated 
    call output_array('Z', 'geom', &
       & clean0(Z_fs(1,:)), 'C', '(5ES13.5)', ascii_fmt)  

    
  else ! quantities depending on the radius 

    call output_array('NS', 'geom', &
       & (/n_s_grid*1.0/), 'C', '(F8.0)', ascii_fmt)
    call output_array('NX', 'geom', &
       & (/n_x_grid*1.0/) , 'C', '(F8.0)', ascii_fmt)
    call output_array('Rref', 'geom', &
       & (/Rref/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('R0', 'geom', &
       & (/R0/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('lfun', 'geom', &
       & (/lfun/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('Bref', 'geom', &
       & (/Bref/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('eps', 'geom', &
       & xgr, 'C', '(5ES13.5)', ascii_fmt)
    call output_array('q', 'geom', &
       & qx, 'C', '(5ES13.5)', ascii_fmt)
    call output_array('shat', 'geom', &
       & shatx , 'C', '(5ES13.5)', ascii_fmt)
    call output_array('Jacobian', 'geom', &
       & clean0(Jacobian_G), 'C', '(5ES13.5)', ascii_fmt)
    call output_array('bmin', 'geom', &
       & bmin, 'C', '(ES13.5)', ascii_fmt)
    call output_array('bmax', 'geom', &
       & bmax, 'C', '(ES13.5)', ascii_fmt)
    ! kthnorm = Rref*sqrt(metric(isg_lfs,2,2))
    call output_array('kthnorm', 'geom', &
       & (/kthnorm/), 'C', '(ES13.5)', ascii_fmt)
    !beta_eq = 2*mu0*p/Bref**2
    call output_array('beta_eq', 'geom', &
       & (/beta_eq/), 'C', '(ES13.5)', ascii_fmt)
    !betaprime_eq = 2*mu0*dpdeps/Bref**2
    call output_array('betaprime_eq', 'geom', &
       & (/betaprime_eq/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('jfunh', 'geom', &
       & (/jfunh/), 'C', '(ES13.5)', ascii_fmt)
    call output_array('jfunl', 'geom', &
       & (/jfunl/), 'C', '(ES13.5)', ascii_fmt)

    ! s-dependent quantities
    call output_array('s_grid', 'geom', &
       & clean0(sgr(1:n_s_grid)), 'C', '(5ES13.5)', ascii_fmt)
    !bn = B/Bref
    call output_array('bn', 'geom', &
       & (/ ((clean0(bn_G(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /), &
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('poloidal_angle', 'geom', &
       & (/ ((clean0(pol_angle(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_eps_eps', 'geom', &
       & (/ ((clean0(metric_G(ix,i,1,1)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_eps_zeta', 'geom', &
       & (/ ((clean0(metric_G(ix,i,1,2)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_eps_s', 'geom', &
       & (/ ((clean0(metric_G(ix,i,1,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_zeta_zeta', 'geom', &
       & (/ ((clean0(metric_G(ix,i,2,2)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_zeta_s', 'geom', &
       & (/ ((clean0(metric_G(ix,i,2,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('g_s_s', 'geom', &
       & (/ ((clean0(metric_G(ix,i,3,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('D_eps', 'geom', &
       & (/ ((clean0(dfun(ix,i,1)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('D_zeta', 'geom', &
       & (/ ((clean0(dfun(ix,i,2)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('D_s', 'geom', &
       & (/ ((clean0(dfun(ix,i,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('E_eps_zeta', 'geom', &
       & (/ ((clean0(efun(ix,i,1,2)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('E_eps_s', 'geom', &
       & (/ ((clean0(efun(ix,i,1,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('E_zeta_s', 'geom', &
       & (/ ((clean0(efun(ix,i,2,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('F', 'geom', &
       & (/ ((clean0(ffun(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('G', 'geom', &
       & (/ ((clean0(gfun(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('H_eps', 'geom', &
       & (/ ((clean0(hfun(ix,i,1)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('H_zeta', 'geom', &
       & (/ ((clean0(hfun(ix,i,2)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('H_s', 'geom', &
       & (/ ((clean0(hfun(ix,i,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('I_eps', 'geom', &
       & (/ ((clean0(ifun(ix,i,1)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('I_zeta', 'geom', &
       & (/ ((clean0(ifun(ix,i,2)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('I_s', 'geom', &
       & (/ ((clean0(ifun(ix,i,3)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('J', 'geom', &
       & (/ ((clean0(jfun(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('K', 'geom', &
       & (/ ((clean0(kfun(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('Bt_frac', 'geom', &
       & (/ ((clean0(Bt_frac(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('R', 'geom', &
       & (/ ((clean0(Rfun(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
    call output_array('Z', 'geom', &
       & (/ ((clean0(Z_fs(ix,i)),i=1, n_s_grid), ix = 1, n_x_grid) /),&
       & 'C', '(5ES13.5)', ascii_fmt)
  endif

  contains
    !> round numbers close to zero to a clean positive zero for the test cases
    !> works also on arrays
    elemental function clean0(input) result(output)
      real, intent(in) :: input
      real :: output
      output = merge(input,+0.0,abs(input) > 1.e-15)
    end function clean0

end subroutine geom_output


!------------------------------------------------------------------------------
!> This subroutine calculates the parallel derivative of log(B).
!> it is assumed that the field line grid is periodic in the 
!> poloidal direction
!------------------------------------------------------------------------------
subroutine logbderiv
  use grid,    only : n_x_grid,n_s_grid 
  integer :: i, ix  
  real    :: lbm2, lbm1, lbp1, lbp2 

  do ix = 1, n_x_grid 
    if (n_s_grid == 1) then
      gfun(ix,1) = 0.
    else 
      do i = 1, n_s_grid
        ! periodicity with respect to i index is done with a modulo.
        lbm2 = bn_G(ix,mod(n_s_grid+(i-1)-2, n_s_grid)+1)
        lbm1 = bn_G(ix,mod(n_s_grid+(i-1)-1, n_s_grid)+1)
        lbp1 = bn_G(ix,mod(n_s_grid+(i-1)+1, n_s_grid)+1)
        lbp2 = bn_G(ix,mod(n_s_grid+(i-1)+2, n_s_grid)+1)
        gfun(ix,i) = ffun(ix,i)*(lbm2 - 8E0*lbm1 + 8E0*lbp1 - lbp2)/ &
           & (12E0*sgr_dist) / bn_G(ix,i)
      end do 
    endif 

  end do

end subroutine logbderiv 


!------------------------------------------------------------------------------
!> This subroutine calculates atan(sqrt((1-eps)/1+eps)*tan(theta/2))
!> needed in the 'circ' equilibrium model
!------------------------------------------------------------------------------
subroutine bangle(pol_mod)

  use grid,      only : n_s_grid 
  use constants, only : pi
 
  integer :: i
  real    :: pol_mod(n_s_grid), dum, shift

  if (n_s_grid.lt.2) then 

    pol_mod(1) = 0.

  else 

    dum = ((1.E0-eps)/(1.E0+eps))**0.5
    pol_mod(1) = atan(dum*tan(pol_angle(1,1)/2.E0))
    do i = 2, n_s_grid 
      ! This can cause a floating point underflow but for normal compilation
      ! it gets ignored and seems not to cause a slowdown
      pol_mod(i) = atan(dum*tan(pol_angle(1,i)/2.E0))
      do while (pol_mod(i)<pol_mod(i-1))
        pol_mod(i) = pol_mod(i) + pi
      end do
    end do 
    shift = pi*floor((pol_mod(1)-pol_angle(1,1)/2.E0)/pi)
    do i = 1, n_s_grid 
      pol_mod(i) = pol_mod(i) - shift
    end do 
  endif 

end subroutine bangle

!------------------------------------------------------------------------------
!> routine to check chease input file reads.  Could also do broadcast
!------------------------------------------------------------------------------

subroutine check_read(ios)

  use general, only : gkw_abort
  use global, only : int2char

  integer, intent(in) :: ios
  integer             :: ios_sum
  integer, save       :: called = 0

  called = called + 1

  ios_sum = abs(ios)
  if (ios_sum > 0) call gkw_abort('Chease input parse / read error '//       &
                        & 'in block ' // trim(int2char(called,3)) //         &
                        & '. Use CHEASE_90_9_3 or later')

end subroutine check_read


!------------------------------------------------------------------------------
!> In order to solve the local problem in s, copy the local elements to the
!> begining of each array used elsewhere in the code. The alternatives are to
!> either change this module to work with ns or to change all the other
!> modules using this to use nspb/e and friends. geom_output and krbal must be
!> called before this. There is a number "2" in this routine which corresponds
!> to the maximum number of ghost points; this should be fixed somewhere else
!> in the code. A similar routine to this appears in rotation.
!------------------------------------------------------------------------------

subroutine parallelize_geom

  use grid,    only : ns, n_s_grid, gs, parallel_s
  use general, only : gkw_abort

  integer :: i1,i2,ib,ie

  if (geom_parallelized) call gkw_abort('parallelize geom called twice!')

  ! i2=global s point; i1=local s point
  ! bn_G is periodic in n_s_grid; the ghost points are used by some schemes.
  ! The metric is NOT "periodic" in n_s_grid and krloc is NOT "periodic" in
  ! n_s_grid. With mode_box, metric and krloc values are copied as if
  ! "periodic" and accessed by connect_parallel in construction of the
  ! "ballooning peroidicity".  This occurs only for calculation of
  ! besselj0_gkw in term VII via wrapper routine bessel_j0.


  ! for parallel_s we need to deal with extra points first
  ib = 1-2 ; ie = ns+2

  ! left local s boundary
  do i1=ib, 1-1
    i2 = gs(i1)
    if (i2 < 1) then ! outside n_s_grid
      if (parallel_s) then
        metric(:,i1,:,:) = metric_G(:,i2+n_s_grid,:,:) 
      end if
      ! if n_s_grid < maximum number of ghost points (presently 2), then map 1
      ! to 1
      if (n_s_grid < 2) then
        bn(:,i1)     = bn_G(:,1)
        alphak(:,i1) = alphakp(:,1) 
      else
        bn(:,i1)     = bn_G(:,i2+n_s_grid)
        alphak(:,i1) = alphakp(:,i2+n_s_grid) 
      end if
    else                ! adjacent s processor
      if (parallel_s) then
        metric(:,i1,:,:)= metric_G(:,i2,:,:)
        alphak(:,i1)    = alphakp(:,i2) 
      end if
      bn(:,i1)     = bn_G(:,i2)
      alphak(:,i1) = alphakp(:,i2)
    end if
  end do

  ! right local s boundary
  do i1=ns+1, ie
    i2 = gs(i1)
    if (i2 > n_s_grid) then ! outside n_s_grid
      if (parallel_s) then
        metric(:,i1,:,:)= metric_G(:,i2-n_s_grid,:,:)
      end if
      ! if n_s_grid < maximum number of ghost points (presently 2), then map 1
      ! to 1
      if (n_s_grid < 2) then
        bn(:,i1)     = bn_G(:,1)  
        alphak(:,i1) = alphakp(:,1)
      else
        bn(:,i1)     = bn_G(:,i2-n_s_grid)
        alphak(:,i1) = alphakp(:,i2-n_s_grid) 
      end if
    else                      ! adjacent s processor
      if (parallel_s) then
        metric(:,i1,:,:)= metric_G(:,i2,:,:)
        alphak(:,i1)    = alphakp(:,i2) 
      end if
      bn(:,i1)     = bn_G(:,i2)
      alphak(:,i1) = alphakp(:,i2)
    end if
  end do

  ! all cases can use the following
  do i1=1, ns
    i2 = gs(i1)
    !This isnt use else where so its commented out
    !here to keep the poloidal angle 'global'. WAH
    !pol_angle(:,i1)    = pol_angle(:,i2)
    ints(i1)         = ints(i2)
    dfun(:,i1,:)     = dfun(:,i2,:)
    efun(:,i1,:,:)   = efun(:,i2,:,:)
    ffun(:,i1)       = ffun(:,i2)
    gfun(:,i1)       = gfun(:,i2)
    hfun(:,i1,:)     = hfun(:,i2,:)
    ifun(:,i1,:)     = ifun(:,i2,:)
    jfun(:,i1)       = jfun(:,i2)
    kfun(:,i1)       = kfun(:,i2)
    Rfun(:,i1)       = Rfun(:,i2)
    metric(:,i1,:,:) = metric_G(:,i2,:,:)
    alphak(:,i1)     = alphakp(:,i2) 
    alphak_xbnd(i1)  = alphak_xbnd(i2)
    bn(:,i1)         = bn_G(:,i2)
    bt_frac(:,i1)    = bt_frac(:,i2)
    dpdpsi_rot(:,i1) = dpdpsi_rot(:,i2)
    dpds_rot(:,i1)   = dpds_rot(:,i2)

  end do

  do i1=0, ns+1
    i2 = gs(i1)
    sgr(i1)=sgr(i2)
  end do

  call parallelize_geom_radial

  geom_parallelized = .true.


end subroutine parallelize_geom

!------------------------------------------------------------------------------
!> In order to solve the local problem in x, copy the local elements to the
!> begining of each array used elsewhere in the code, as for s above
!> A similar routine will be required in components if profiles are read
!------------------------------------------------------------------------------
subroutine parallelize_geom_radial

  use grid,    only : nx, gx

  integer :: i2, i1

  ! all cases can use the following
  do i1=1, nx
    i2 = gx(i1)
    alphak(i1,:)     = alphak(i2,:)
    shift_end_grid(i1) = shift_end_grid(i2)
    
    ! For global profiles
    ! actually never used after this routine, but parallelise just to be safe.
    ! It would be better to make them private to geom and provide only the
    ! mid point values to mode.
    !qx(i1)    =  qx(i2)
    !shatx(i1) =  shatx(i2)
    !xgr(i1)  = xgr(i2)
    !WAH at the moment the above doesnt need to be performed
    
    ! those below here should not be needed unless geom is radially global
    pol_angle(i1,:)  = pol_angle(i2,:)
    dfun(i1,:,:)     = dfun(i2,:,:)
    efun(i1,:,:,:)   = efun(i2,:,:,:)
    ffun(i1,:)       = ffun(i2,:)
    gfun(i1,:)       = gfun(i2,:)
    hfun(i1,:,:)     = hfun(i2,:,:)
    ifun(i1,:,:)     = ifun(i2,:,:)
    jfun(i1,:)       = jfun(i2,:)
    kfun(i1,:)       = kfun(i2,:)
    Rfun(i1,:)       = Rfun(i2,:)
    metric(i1,:,:,:) = metric(i2,:,:,:)
    bn(i1,:)         = bn(i2,:)
    bt_frac(i1,:)    = bt_frac(i2,:)

    bmin(i1) = bmin(i2)
    bmax(i1) = bmax(i2)
  end do

end subroutine parallelize_geom_radial


!------------------------------------------------------------------------------
!> The q-profile function
!> For the flux_tube the local values are copied. For the non-flux tube the
!> profile is calculated according to the switch prof_type and the coefficients
!> qprof_type. The various options are (qc() = qprof_coef())
!>
!> prof_type = 'parabolic'  
!>   q(psi)    = ( qc(1) + qc(2) * psi^2 ) 
!>   shat(psi) = 2 * qc(2) * psi^2 / (qc(1) + qc(2)*psi^2) 
!>   qc(1) = q(0)   qc(2) = (q(a) - q(0))*(R/a)^2
!>   For reasons why this has not been merged with parabolic2, see there.
!>
!> prof_type = 'parabolic2'
!>   qq = qprof_coef(1) + qprof_coef(2) * psi + qprof_coef(3) * psi ** 2
!>   ss = (qprof_coef(2) * psi + 2 * qprof_coef(3) * psi ** 2) / qq
!>   This is similar to parabolic, just with an addtional linear term.
!>   \note This is not merged with parabolic, because adding there a linear term
!>   would result either in a unlogic order of the parameters or one would
!>   have to break compability with previous versions (old input files). This
!>   error would also be hard to notice.
!>
!> prof_type = 'orb'  (for circular geometry) 
!>   q(psi)    = (qc(1) + qc(2) * psi^2) / sqrt(1 - psi^2) 
!>   shat(psi) = 2.*qc(2)*psi^2 / (qc(1) + qc(2)*psi^2) - psi^2 / (1 - psi^2)
!>   In ORB5, NEMORB one specifies qbar = q * sqrt(1-(r/R)^2) the meaning of
!>   the coefficents qc is in terms of qbar 
!>   qc(1) = q(0)  qc(2) = (qbar(a) - q(0))*(R/a)^2
!>
!> prof_type = 'wesson'
!>   nu = qc(2)
!>   q0 = qc(3)
!>   qq = q0*(psi^2)*qc(1)/(1 - exp((nu+1)*log(1- qc(1)*psi*^2)))
!>   ss = 2*(1-(qq/q0)*(nu+1)*exp(nu*log(1-qc(1)*psi^2)))
!>   \warning For this profile qprof_coef(1)*psi*psi < 1.E0 has to hold,
!>   otherwise the term in the logarithm gets negative. If this happens, the
!>   simulation will abort.
!>
!> prof_type = 'rexp'
!>   qq = qc(2)*psi*exp(qc(1) * psi)
!>   ss = 1 + psi * qc(1)
!>   qc(3) us a cutoff, values of qq lower than this will be set to qc(3) and
!>   ss at these position(s) to zero.
!>
!> prof_type = 'mishchenko'
!>   qq = qc(1) + (1 - qc(1)) (psi / qc(2)) ** qc(3)
!>   ss = (qc(3) * (1 - qc(1)) (psi / qc(2)) ** qc(3)) / qq
!>   This has for example been used in Mishchenko and Zocco, Phys. Plasmas 19,
!>   122104 (2012).
!>   The corespondence to the names used in this paper is
!>   qc(1) <-> q_0
!>   qc(2) <-> r_c (in our case a \f$ \psi_0 \f$)
!>   qc(3) <-> p_q (set to 1 in the paper).
!>   \warning Value combinations of qc(1) > 1 and/or qc(3) < 1 might result
!>   in a decreasing q and negative shear. Be careful what you want, when using
!>   these parameterranges.
!------------------------------------------------------------------------------
subroutine q_profile(psi,qq,ss)
  use global, only : r_tiny
  use control, only : flux_tube
  use general, only : gkw_abort
  use mpiinterface, only : root_processor

  real, intent(in)  :: psi
  real, intent(out) :: qq, ss

  real :: lq, qv, nu, q0, prefac

  if (flux_tube) then
    ! copy in the local values
    qq = q
    ss = shat
    return
  else

    select case(prof_type)
    case('parabolic')
      qq = qprof_coef(1) + qprof_coef(2) * psi ** 2
      ss = 2.* qprof_coef(2)* psi**2 / (qprof_coef(1) + qprof_coef(2)*psi**2)
    case('parabolic2')
      qq = qprof_coef(1) + psi*qprof_coef(2) + qprof_coef(3)*psi**2
      ss = (2.* qprof_coef(3)* psi**2 + qprof_coef(2)*psi)/ (qprof_coef(1) + psi*qprof_coef(2) + qprof_coef(3)*psi**2)
    case('orb')
      qq = (qprof_coef(1) + qprof_coef(2) * psi **2)/sqrt(1-psi**2)
      ss = 2.* qprof_coef(2)* psi**2 / (qprof_coef(1) + qprof_coef(2)*psi**2) &
         & - psi**2 / (1 - psi**2)
    case('wesson')
      nu = qprof_coef(2)
      q0 = qprof_coef(3)
      prefac = q0/(nu+1)
      qq = q0*(psi*psi)*qprof_coef(1)/(1 - (1- &
         & qprof_coef(1)*psi*psi)**(nu+1))
      ss = 2*(1-(qq/q0)*(nu+1)*(1-qprof_coef(1)*psi*psi)**nu)
      if(qprof_coef(1)*psi*psi > 1.E0 .and. mod(nu+1.E0,2.E0) > r_tiny) then
        if(root_processor)then
          call gkw_abort('The parameter 1-a^2psi^2 becomes negative in Wesson q-profile. &
            & Set the (R/a)^2 parameter to be smaller.')
        endif
      endif
    case('rexp')
      lq = qprof_coef(1)
      qv = qprof_coef(2)
      q0 = qprof_coef(3)
      qq = qv*psi*exp(psi*lq)
      if (qq.lt.q0) then
        qq = q0
        ss = 0.E0
      else
        ss = 1 + psi * lq
      endif

    case('mishchenko')
      qq = qprof_coef(1) + (1 - qprof_coef(1))*(psi / qprof_coef(2)) ** qprof_coef(3)
      ss = (qprof_coef(3)                                                    &
         & * (1 - qprof_coef(1))*(psi / qprof_coef(2)) ** qprof_coef(3)) / qq
    case default 
      call gkw_abort('unknown q-profile type')
    end select

    return

  endif

end subroutine q_profile 


!------------------------------------------------------------------------------
!> Calculate the geometry profiles. This includes the q shat profiles (for which 
!> a separate routine is called as well as the parameters of the Miller geom. 
!------------------------------------------------------------------------------
subroutine geom_profiles

  use grid,         only : n_x_grid 
  use general,      only : gkw_abort
  use io,           only : get_free_file_unit 
  use mpiinterface, only : root_processor, mpibcast
  
  integer :: ix, file_unit 
  real    :: xdum
  character(300) :: line 

  !----------------------------------------------------------------------------
  ! Geometry parameters are read from file 
  !----------------------------------------------------------------------------
  if (prof_type == 'file') then 
    
    if(root_processor)then
      call get_free_file_unit(file_unit)
      open(file_unit,file='input.prof',status="old")
      line = ''
      do while (line(1:9) /= '#Geometry')  
        read(file_unit,'(A)',ERR = 100,END=100) line 
      end do 
      do ix = 1, n_x_grid 
        if (geom_type == 'miller') then 
          read(file_unit,*,ERR=100)xdum,qx(ix),shatx(ix),kappax(ix),deltax(ix), & 
            & squarex(ix), skappax(ix), sdeltax(ix),ssquarex(ix),Zmilx(ix), &
            & dRmilx(ix), dZmilx(ix), gradpx(ix)
          if (abs(xdum-xgr(ix)).gt.1e-5) call gkw_abort('Inconsistent xgr in input.prof')
        else 
          read(file_unit,*)xdum,qx(ix),shatx(ix)
          if (abs(xdum-xgr(ix)).gt.1e-5) call gkw_abort('Inconsistent xgr in input.prof')
        endif 
      end do 
      close(file_unit) 
    endif 
 
    call mpibcast(qx,    n_x_grid) 
    call mpibcast(shatx, n_x_grid)
    
    if (geom_type == 'miller') then 
      call mpibcast(kappax,   n_x_grid) 
      call mpibcast(deltax,   n_x_grid) 
      call mpibcast(squarex,  n_x_grid) 
      call mpibcast(skappax,  n_x_grid)
      call mpibcast(sdeltax,  n_x_grid) 
      call mpibcast(ssquarex, n_x_grid) 
      call mpibcast(Zmilx,    n_x_grid) 
      call mpibcast(dRmilx,   n_x_grid)
      call mpibcast(dZmilx,   n_x_grid) 
      call mpibcast(gradpx,   n_x_grid) 
    endif 
                
    return 
    100 call gkw_abort('Error while reading file input.prof') 
 
  endif
 
  !----------------------------------------------------------------------------
  ! Profiles are analytic or constant 
  !----------------------------------------------------------------------------
  do ix = 1, n_x_grid 
    call q_profile(xgr(ix),qx(ix),shatx(ix)) 
  end do   

  ! The miller parameters 
  if (geom_type == 'miller') then 
    
    do ix = 1, n_x_grid 
      kappax(ix)   = kappa 
      deltax(ix)   = delta 
      squarex(ix)  = square
      skappax(ix)  = skappa 
      sdeltax(ix)  = sdelta 
      ssquarex(ix) = ssquare 
      Zmilx(ix)    = Zmil 
      dRmilx(ix)   = dRmil 
      dZmilx(ix)   = dZmil 
      gradpx(ix)   = gradp 
    end do 

  endif 
  
end subroutine geom_profiles


!-----------------------------------------------------------------------------
!> Quadratic interpolation
!> x = x-axis, y = y-axis, n1 = number of points in x and y
!> x_interp = interpolation of x-axis, y_interp = interpolation of y-axis  
!> n2 = number of points in x_interp and y_interp
!> nearly all the CPU time of Miller is spent in this routine
!-----------------------------------------------------------------------------
subroutine interpquad (x, y, n1,n2, x_interp,y_interp)

  integer :: n1, tmp, i, ith,n2

  real :: x(n1), y(n1), diff(n1)
  real :: x_interp(n2)
  real :: y_interp(n2)
  
  do i = 1, n2

    do ith = 1, n1

      diff(ith) = abs(x_interp(i) - x(ith))
      !Could find the minloc in here
    end do

    tmp = minloc(diff,1)

    if (tmp == n1) then
      tmp = n1-1
    end if

    if (tmp == 1) then
      tmp = 2
    end if

    y_interp(i) = y(tmp-1)*(x_interp(i) - x(tmp)) * (x_interp(i) -        &
    &  x(tmp+1)) / ((x(tmp-1) - x(tmp))*(x(tmp-1)-x(tmp+1))) +          &
    &  y(tmp)*(x_interp(i) - x(tmp+1)) * (x_interp(i) -               &
    &  x(tmp-1)) / ((x(tmp) - x(tmp+1))*(x(tmp)-x(tmp-1))) +            &
    &  y(tmp+1)*(x_interp(i) - x(tmp-1)) * (x_interp(i) -             &
    &  x(tmp)) / ((x(tmp+1) - x(tmp-1))*(x(tmp+1)-x(tmp)))

  end do
  
end subroutine interpquad

!-----------------------------------------------------------------------------
!> Subroutine that calculates numerical integrals.
!> This is used for the construction of Miller geometry.
!-----------------------------------------------------------------------------
subroutine simpson_numerical_integ(x,y,N,integ)
  
  use constants, only : pi
  use global, only : r_tiny
  !> number of points in x and y
  integer, intent(in) :: N
  !> x-axis
  real, intent(in) :: x(N)
  !> y-axis
  real, intent(in) :: y(N)
  !> calculated value of the integral \int_0^x.
  !> Note that the value x=0 is in the centre of the grid.
  real, intent(out) :: integ(N)
  real :: x_interp(2*N-1),y_interp(2*N-1)
  integer :: N2,ith,ii
  
  !>  the index of the zero in the x grid. This is used to separate
  !>  the domain of integration between negative and positive range
  integer :: ixzero

  ixzero = minloc(abs(x),1)
  
  N2 = 2*N-1
  do ith = 1, N2
    x_interp(ith) = x(1) + (x(N)-x(1)) * (ith-1)/(N2-1)
  end do
  call interpquad (x,y,N,N2,x_interp,y_interp)
  
  integ = 0.0
  
  do ith =1, N
    if (ith == ixzero) then
      integ(ith) = 0.
    else if (ith < ixzero) then
      ! integrate from x=zero to a finite x<0.
      do ii = ith, ixzero-1
        integ(ith) = integ(ith) - (x(ii+1)-x(ii))*(y(ii)+4.*y_interp(2*ii)        &
        + y(ii+1))/6.
      end do
    else if (ith > ixzero) then
      ! integrate from x=zero to a finide x>0.
      do ii = ixzero, ith-1
        integ(ith) = integ(ith) + (x(ii+1)-x(ii))*(y(ii)+4.*y_interp(2*ii)        &
        + y(ii+1))/6.
      end do
    end if
  end do

end subroutine simpson_numerical_integ

end module geom

