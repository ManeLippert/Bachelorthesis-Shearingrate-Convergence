!------------------------------------------------------------------------------
!> Module containing the species data (and profiles in global case).
!> Also provides functions to rotation for calculating centrifugal potential.
!------------------------------------------------------------------------------
module components

  use global, only : lenswitch

  implicit none

  private

  ! public subroutines 
  public :: components_read_nml_spcg, components_bcast_nml_spcg
  public :: components_check_params_spcg, components_write_nml_spcg
  public :: components_read_nml_spec, components_bcast_nml_spec
  public :: components_check_params_spec, components_write_nml_spec
  public :: components_input_species, components_profiles
  public :: components_beta_cal
  public :: components_deallocate

  ! public variables 
  public :: adiabatic_electrons, amp_init, amp_init_imag, kwid_ini
  public :: components_allocate, de, finit, lfinit_radial_dirichlet, cf_phi_max
  public :: cf_mass_weight, fp, mas, pbg, signz, cf_quasineutral, cf_dtf1
  public :: zeff_sp, tmp, tp, types, vp, vthrat, rhostar, rhostar_linear
  public :: tgrid, tgrid_G, de_Gx, tmp_Gx
  public :: tmp_G, signz_G, vthrat_G, mas_G, de_G, vp_G, fp_G, tp_G
  public :: nsps, iadia, tearingmode, isl_mode, tear_zero_epar
  public :: dgrid, dgrid_G, isl_shear, max_rho
  public :: veta, veta_prime, veta_Gx, betaprime_type
  public :: wstar, isl_Ls, isl_rot_freq  
  public :: tm_drive, ecurx, tm_start, tm_sat
  public :: psi_0, delta_psi_0 
  public :: energetic_particles
  public :: icrh_params, isg_icrh, icrh_norm, icrh_fmin, icrh_rln 
  public :: icrh_model
  public :: rlt_gauss
  public :: vpar_mean
  public :: quench_switch, n_quench
  public :: finit_imod, imod_init, mode_persist, ind_persist
  public :: amp_imod, amp_imod_imag, amp_zon, amp_zon_imag
  public :: init_coef

  public :: n_spec

  interface components_write_nml_spcg
    module procedure components_read_nml_spcg
  end interface

  !> Ratio of the grid temperature (nsp) to the reference temperature 
  !> In the local model this is equivalent to the temperature of the species.
  real, save, allocatable :: tgrid(:)  
  ! and for ALL species tgrid_G(:)
  real, save, allocatable :: tgrid_G(:)

  !> Ratio of the grid density (nsp) to the reference density 
  !> In the local model this is equivalent to the density of the species.
  real, save, allocatable :: dgrid(:)  
  ! and for ALL species tgrid_G(:)
  real, save, allocatable :: dgrid_G(:)

  !> the ratio of the mass to the reference value: local proc: mas(nsp+iadia)
  real, save, allocatable :: mas(:)
  !> and for ALL species mas_G(nsps)
  real, save, allocatable :: mas_G(:)

  !> The charge number of the particles (Z): local proc: signz(nsp+iadia)
  real, save, allocatable :: signz(:)
  !> and for ALL species signz_G(nsps)
  real, save, allocatable :: signz_G(:)

  !> The ratio of the temperature to the reference temp: local: tmp(nsp+iadia)
  real, save, allocatable :: tmp(:,:)
  !> and for ALL species tmp_G(nsps)
  real, save, allocatable :: tmp_G(:,:)
  !> Temperature global in x tmp_Gx(1:n_x_grid,nsp+iadia)
  real, save, allocatable :: tmp_Gx(:,:)

  !> The density of the species: local: de(nx,nsp+iadia)
  real, save, allocatable :: de(:,:)
  !> and for ALL species de_G(nx,nsps)
  real, save, allocatable :: de_G(:,:)
  !> Density global in x de_Gx(1:n_x_grid,nsp+iadia)
  real, save, allocatable :: de_Gx(:,:)

  
  !> The temperature gradient \f$ 1/LT = \nabla T / T \f$ (normalized to the local 
  !> temperature 
  real, save, allocatable :: tp(:,:)
  !> and for ALL species de_G(nsps)
  real, save, allocatable :: tp_G(:,:)

  !> The density gradient of the species 1/Ln (normalized to the local density)
  real, save, allocatable :: fp(:,:)
  !> and for ALL species fp_G(nsps)
  real, save, allocatable :: fp_G(:,:)

  !> The gradient of the toroidal velocity R grad U / v_thref
  real, save, allocatable :: vp(:,:)
  !> and for ALL species vp_G(nsps)
  real, save, allocatable :: vp_G(:,:)
  !> The gradient of the electron current R grad U / v_thref
  real, save, allocatable :: ecurx(:)

  !>The density profile type (global) 
  character (len=lenswitch), save, allocatable :: de_prof_type(:)
  
  !>The temperature profile type (global) 
  character (len=lenswitch), save, allocatable :: tmp_prof_type(:)

  !>The density profile coefficients 
  real, save, allocatable :: de_prof_coef(:,:)

  !>The temperature profile coefficients 
  real, save, allocatable :: tmp_prof_coef(:,:) 

  !> Amplitude and width of a gaussianwise background temperature gradient
  real, save, allocatable :: rlt_gauss_params(:,:) 

  !> the type of the species (only necessary for non-maxwellian 
  character (len=7), save, allocatable :: types(:)
  !> and for ALL species types_G(nsps)
  character (len=7), save, allocatable :: types_G(:)

  !> A parameter that can be used to define the background (only 
  !> necessary for non-maxwelian backgrounds 
  real, save, allocatable :: pbg(:,:)
  !> and for ALL species pbg_G(nsps)
  real, save, allocatable :: pbg_G(:,:)

  !> vthrat is the v_R from the manual, i.e. the normalised thermal velocity
  !> of the species.
  !> The ratio of the 'grid' thermal velocity [sqrt(2 T_G / m)] to the 
  !> reference thermal velocity [2 sqrt(2 Tref / mref]. In the local 
  !> model the grid temperature is equal to the temperature of 
  !> the species, and vthrat is the ratio of the termal velocity of 
  !> the species to the reference velocity. 
  real, save, allocatable :: vthrat(:)
  !> and for ALL species vthrat_G(nsps)
  real, save, allocatable :: vthrat_G(:)

  !> The reference plasma beta 
  real, save :: beta_ref 

  !> maximum species larmor radius
  real, save :: max_rho = 0

  ! Note that beta_ref is used in the Ampere's law
  ! and is defined as beta=2*mu0*n_ref*T_ref/B_ref^2
  ! It is different from beta_eq=2*mu0*p/B_ref^2 that is
  ! calculated in geom.f90

  !> Switch for beta_ref input
  character (len=lenswitch), save :: beta_type 

  !> The reference beta'
  real, save :: betaprime_ref

  !> Switch for beta' input
  character (len=lenswitch), save :: betaprime_type 

  !> Reference beta value used in the equations. Value will depend 
  !> on beta_type [either from input (beta_ref or beta) or from the 
  !> equilibrium]  
  real, save, allocatable :: veta(:) 

  !> only used in input to be backwards compatible 
  real, save :: beta 

  !> value of beta^prime used in the equations. Value will depend 
  !> on the beta_prime_type. 
  real, save, allocatable :: veta_prime(:)  

  !> Logical that is true if the electrons are adiabatic 
  logical, save :: adiabatic_electrons

  !> Total number of species, incl. adiabatic electrons.
  !> (This specifies in particular how many species namelists are read from the
  !> input file)
  integer, save :: n_spec

  !> The initial amplitude of the distribution function (real and imaginary)
  real, save :: amp_init, amp_init_imag

  !> a mode which is initialised in a different way
  integer, save :: imod_init
  !> real amplitude to initialise the mode imod_init with
  real, save :: amp_imod, amp_imod_imag
  !> real amplitude factor, to initialise the zonal mode
  real, save :: amp_zon, amp_zon_imag
  !> to fix a mode with the initial condition
  logical, save :: mode_persist
  !> indices to access the persistent mode
  integer, save, allocatable :: ind_persist(:)
  
  !> list of coefficients used for specific initialization (cos-zonal-tmp)
  real, save :: init_coef(5)

  !> parameter for adjusting the width of initial gaussians in kspace
  !> (used only with cosine5, cosine6 and hybrid initial conditions)
  real, save :: kwid_ini

  !> If tearingmodes are to be used.
  logical, save :: tearingmode

  !> If tearing modes are driven self-consistently from current profile
  logical, save :: tm_drive

  !> If tearingmodes addition to phi is to be used, to give zero E_parallel
  logical, save :: tear_zero_epar

  !> For experimental global effects (ONLY), for flux tube, rhostar = 0. 
  real, save :: rhostar

  !> For tearing mode studies the normalised island width
  real, save :: wstar

  !> For tearing mode studies the perturbation damping length;
  !> is the damping length scale to minimize the jump at the boundary
  real, save :: isl_Ls
  
  !> The tearing_mode shear (ONLY when geometric shear is zero, and 2D) 
  real, save :: isl_shear

  !> For tearing mode studies the magnetic island rotation frequency
  real, save :: isl_rot_freq

  !> For tearing mode studies the magnetic island rotation frequency
  integer, save :: isl_mode
  
  !> Times between which to grow the imposed tearing mode gradually
  real, save :: tm_start, tm_sat

  !> Radial position and width of the tanh to force radial boundary
  !> conditions of A|| to zero in the non spectral version
  real, save :: psi_0, delta_psi_0

  !> Parameters for the definition of the EP distribution function
  logical, save :: energetic_particles
  
  !> Z effective of the species inputs (distinct from collisions inputs)
  real, save :: zeff_sp

  !> Estimate of maximum possible centrifugal potential
  real, save :: cf_phi_max(4)
  
  !> variables describing anisotropic ICRH heated species (experimental)
  !> see M L Reinke et al, Plasma Phys. Control. Fusion 54 045004 (2012)
  !> global species number of ICRH species (only if independent species)
  integer, save :: isg_icrh = 0
  ! fminority, tperp_tpar, dtperp_tpar, z_min
  real, save :: icrh_params(4) 
  ! select between Reinke (1), Bilato (2) anisotropy models
  integer, save :: icrh_model = 0
  
  real, save :: icrh_norm        ! B_R0
  real, save :: icrh_fmin = 0.0  ! nmR0/nref
  real, save :: icrh_rln  = 0.0  ! RLN of minority at R0
  real, save :: tperp_tpar, dtperp_tpar !< dummies
  
  !> list of wavenumbers (indices, i.e. the zero-mode corresponds to 1),
  !> from/to which nonlinear transfer shall be suppressed.
  integer, save, public :: quench_modes(64)
  
  !> switch to specify how to quench modes
  character (len = lenswitch), save :: quench_switch
  
  !> integer to specify initial mode quenching
  integer, save :: n_quench

  !> This charater string determines how the distribution is initialized
  !> allowed are 'noise' and 'cosine2', 'hybrid' 
  character (len = lenswitch), save :: finit
  
  !> This charater string determines how the distribution is initialized
  !> when a mode is specified using imod_init
  character (len = lenswitch), save :: finit_imod

  !> If true, then the distribution is multiplied by a factor that makes it
  !> zero at the radial boundaries.
  logical, save :: lfinit_radial_dirichlet

  !> total global number of species, including adiabatic ones
  integer, save :: nsps

  !> number of adiabatic species - this is either 0 or 1
  integer, save :: iadia = 1

  !>Dummies for reading species namelists
  real, save    :: mass, z, temp, dens, rlt, rln, uprim, param(2)
  character (len=lenswitch) :: dens_prof_type, temp_prof_type
  real, save    :: dens_prof_coef(5), temp_prof_coef(5) 
  character (len=7), save :: background
  real, parameter :: spbad = -9001.0

  !> dummies for reading gausswise gradients for every species
  real, save    :: rlt_gauss(2)
  real, save    :: vpar_mean

  namelist /species/ mass, z, temp, dens, rlt, rln, uprim, background, &
     & param, dens_prof_type, dens_prof_coef, &
     & temp_prof_type, temp_prof_coef, rlt_gauss
  
  !> rhostar used in linear_terms.F90 to switch off effects from parallel
  !> derivative of terms 2,5,7,8 in the manual.
  real, save :: rhostar_linear
  
  !> Plasma beta global in x 
  real,    allocatable, save :: veta_Gx(:)

contains


!------------------------------------------------------------------------------
!> read (or write) the species general namelist
!------------------------------------------------------------------------------

subroutine components_read_nml_spcg(file_unit,io_stat,lwrite)
  use grid, only : lx
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  integer, intent(in)  :: file_unit
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  ! Local dummy for read only, will be copied into isl_Ls, in order
  ! not to break compatibility with old tearing mode input files; what
  ! should have been 'isl_Ls' is called 'Ls' in the input file.
  real :: Ls

  namelist /spcgeneral/ beta, beta_ref, beta_type, betaprime_ref,   &
       & betaprime_type, adiabatic_electrons,                       &
       & amp_init, amp_init_imag, finit, lfinit_radial_dirichlet,   &
       & kwid_ini, quench_modes,                                    &
       & rhostar, wstar, Ls, isl_rot_freq, tearingmode, isl_mode,   &
       & tear_zero_epar, isl_shear, tm_drive, tm_start, tm_sat,     &
       & psi_0, delta_psi_0, icrh_params,                           &
       & energetic_particles, vpar_mean, quench_switch, n_quench,   &
       & imod_init, amp_imod, amp_imod_imag, amp_zon, amp_zon_imag, &
       & mode_persist, finit_imod, init_coef

  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then
    
      ! Set the default values
      adiabatic_electrons = .false.
  
      ! zero beta is the default. Note that beta can affect electro-
      ! static runs because of the finite beta contribution to the 
      ! drift 
      beta_ref = 0.E0

      ! deprecated, kept for backward compatibility
      beta = 0.E0

      ! uses input given in beta_ref
      beta_type = 'ref'

      ! zero betaprime by default
      betaprime_ref = 0.E0

      ! default: use the reference value betaprime_refs
      betaprime_type = 'ref'

      ! The initial amplitude (real and imaginary)
      amp_init = 1e-3 
      amp_init_imag = 0.0
      kwid_ini = 0.3

      ! The choice of initialization 
      finit = 'cosine2'
      lfinit_radial_dirichlet = .false.

      quench_modes = 0
      
      quench_switch = 'except'
      n_quench = 0
      
      ! rhostar is zero for a flux tube
      rhostar = 0

      ! tearingmode
      tearingmode=.false.
      tear_zero_epar=.false.
      isl_mode=2
      tm_drive=.false.
      tm_start=-10.0
      tm_sat=0.0

      ! The magnetic island width
      wstar = 0.E0
      psi_0 = lx/3.E0
      delta_psi_0 = lx/20.E0
      
      ! No ICRH anisotropy by default
      icrh_params(1) = 0.0
      icrh_params(2) = 1.0
      icrh_params(3) = 0.0
      icrh_params(4) = 1.0
      
      !Damping length in the tearing mode that damps away the 
      !island mode to prevent discontinuous boundaries
      Ls = 2.E0
      
      !> The tearing_mode shear (ONLY when geometric shear is zero, and 2D) 
      isl_shear = 0.E0

      !The rotation frequency of the magnetic island
      isl_rot_freq = 0.E0

      !The EP distribution function parameters
      energetic_particles = .false.
  
      ! default vpar_mean is zero
      vpar_mean = 0.E0

			! specification of mode that is used for initializaton if fdisi
      imod_init = 0
      amp_imod = amp_init
      amp_imod_imag = amp_init_imag
      
      ! scaling factor of zonal mode used for intialization of fdisi
      amp_zon = 1.0
      amp_zon_imag = 0.0
      
      ! switch for setting mode specified by imod_init persistent
      mode_persist = .false.
      
      ! The choice of initialization (when imod_init = .true.)
      finit_imod = 'cosine2'
      
      ! Init coefficient list contains: 
      ! 1. Radial wave number of the zonal sinusoidal perturbation
      ! 2. Phase in the argument of the sinusoidal perturbation
      ! 3. Amplitude of the zonal sinusoidal tempearture perturbation
      ! 4. Phase in the argument of the sinusoidal temperature perturbation
      ! 5. Normlization of the inital zonal (density) perturbation such that
      !    no temperature perturbation is introdcued
      init_coef(1) = 1.0
      init_coef(2) = 0.0
      init_coef(3) = 0.0
      init_coef(4) = 0.0
      init_coef(5) = 0.0

      ! read the general parameters 
      read(file_unit,NML=spcgeneral,IOSTAT=io_stat)

      isl_Ls=Ls

    end if
  else
    Ls=isl_Ls  
    if(root_processor) write(file_unit,NML=spcgeneral)

    call write_run_parameter('spcgeneral', 'beta', beta)
    call write_run_parameter('spcgeneral', 'beta_ref', beta_ref)
    call write_run_parameter('spcgeneral', 'beta_type', beta_type)
    call write_run_parameter('spcgeneral', 'betaprime_ref', betaprime_ref)
    call write_run_parameter('spcgeneral', 'betaprime_type', betaprime_type)
    call write_run_parameter('spcgeneral', &
       & 'adiabatic_electrons', adiabatic_electrons)
    call write_run_parameter('spcgeneral', 'amp_init', amp_init)
    call write_run_parameter('spcgeneral', 'finit', finit)
    call write_run_parameter('spcgeneral', 'lfinit_radial_dirichlet', lfinit_radial_dirichlet)
    call write_run_parameter('spcgeneral', 'quench_modes', quench_modes)
    call write_run_parameter('spcgeneral', 'rhostar', rhostar)
    call write_run_parameter('spcgeneral', 'wstar', wstar)
    call write_run_parameter('spcgeneral', 'Ls', Ls)
    call write_run_parameter('spcgeneral', 'isl_rot_freq', isl_rot_freq)
    call write_run_parameter('spcgeneral', 'tearingmode', tearingmode)
    call write_run_parameter('spcgeneral', 'isl_mode', isl_mode)
    call write_run_parameter('spcgeneral', 'tear_zero_epar', tear_zero_epar)
    call write_run_parameter('spcgeneral', 'isl_shear', tm_drive)
    call write_run_parameter('spcgeneral', 'tm_start', tm_sat)
    call write_run_parameter('spcgeneral', 'psi_0', delta_psi_0)
    call write_run_parameter('spcgeneral', 'icrh_params', icrh_params)
    call write_run_parameter('spcgeneral', &
       & 'energetic_particles', energetic_particles)
    call write_run_parameter('spcgeneral', 'vpar_mean', vpar_mean)
    call write_run_parameter('spcgeneral', 'quench_switch', quench_switch)
    call write_run_parameter('spcgeneral', 'n_quench', n_quench)
    call write_run_parameter('spcgeneral', 'imod_init', imod_init)
    call write_run_parameter('spcgeneral', 'amp_imod', amp_imod)
    call write_run_parameter('spcgeneral', 'amp_imod_imag', amp_imod_imag)
    call write_run_parameter('spcgeneral', 'amp_zon', amp_zon)
    call write_run_parameter('spcgeneral', 'amp_zon_imag', amp_zon_imag)
    call write_run_parameter('spcgeneral', 'mode_persist', mode_persist)
    call write_run_parameter('spcgeneral', 'finit_imod', finit_imod)
    call write_run_parameter('spcgeneral', 'init_coef', init_coef)

  end if

  ! used in linear_terms.F90 for effects from parallel derivation
  rhostar_linear = rhostar

end subroutine components_read_nml_spcg


!------------------------------------------------------------------------------
!> broadcast species general namelist of components
!------------------------------------------------------------------------------

subroutine components_bcast_nml_spcg

  use mpiinterface, only : mpibcast 

  call mpibcast(beta,                     1)
  call mpibcast(beta_ref,                 1)
  call mpibcast(betaprime_ref,            1)
  call mpibcast(beta_type,        lenswitch)
  call mpibcast(betaprime_type,   lenswitch)
  call mpibcast(amp_init,                 1)
  call mpibcast(amp_init_imag,            1)
  call mpibcast(kwid_ini,                 1)
  call mpibcast(finit,            lenswitch)
  call mpibcast(lfinit_radial_dirichlet,  1)
  call mpibcast(quench_modes, size(quench_modes))
  call mpibcast(adiabatic_electrons,      1)
  call mpibcast(tearingmode,              1)
  call mpibcast(tm_drive,                 1)
  call mpibcast(tm_start,                 1)
  call mpibcast(tm_sat,                   1)
  call mpibcast(psi_0,                    1)
  call mpibcast(delta_psi_0,              1)
  call mpibcast(tear_zero_epar,           1)
  call mpibcast(wstar,                    1)
  call mpibcast(rhostar,                  1)
  call mpibcast(rhostar_linear,           1)
  call mpibcast(isl_Ls,                   1)
  call mpibcast(isl_rot_freq,             1)
  call mpibcast(isl_mode,                 1)
  call mpibcast(isl_shear,                1)
  call mpibcast(icrh_params,              4)
  call mpibcast(energetic_particles,      1)
  call mpibcast(quench_switch,    lenswitch)
  call mpibcast(n_quench,                 1)
  call mpibcast(amp_imod,                 1)
  call mpibcast(amp_imod_imag,            1)
  call mpibcast(amp_zon,                  1)
  call mpibcast(amp_zon_imag,             1)
  call mpibcast(imod_init,                1)
  call mpibcast(mode_persist,             1)
  call mpibcast(finit_imod,       lenswitch)
  call mpibcast(init_coef,  size(init_coef))


end subroutine components_bcast_nml_spcg


!------------------------------------------------------------------------------
!> check species general params
!> must be called after geom_check_params
!------------------------------------------------------------------------------

subroutine components_check_params_spcg

  use control,      only : zonal_adiabatic, nlapar, nlbpar, flux_tube
  use geom,         only : geom_type, shat, gradp_type
  use general,      only : gkw_warn, gkw_exit
  use global,       only : r_tiny
  use grid,         only : number_of_species, lx, psil, psih

  ! Legacy use
  if (abs(beta) > r_tiny) then
    beta_ref = beta
    call gkw_warn('spcgeneral: beta_ref value has been overriden by beta. &
                        & "beta" input deprecated, use "beta_ref" instead')
    beta =0.
  end if
  ! beta should not appear beyond this point (only beta_ref)

  if (beta_ref < 0.) then
    call gkw_exit('spcgeneral: I do not understand negative beta_ref')
  end if

  if (betaprime_ref > 0.) then
    call gkw_warn('spcgeneral: betaprime_ref>0 corresponds to a hollow&
                 & pressure profile. Be sure it is what you wish...')
  end if

  if (beta_type /= 'ref' .and. beta_type /= 'eq' .and.                  &
     & beta_type /= 'file') then
     call gkw_exit('spcgeneral: unknown beta_type option, available are&
                  & "ref", "file" and "eq".')
  end if

  if (betaprime_type /= 'ref' .and. betaprime_type /= 'sp' .and.        &
     & betaprime_type /= 'eq' .and. betaprime_type /= 'file') then
     call gkw_exit('spcgeneral: unknown betaprime_type option, available&
                  & are "ref", "file", "sp" and "eq".')
  end if

  if (geom_type /= 'chease' .and. beta_type == 'eq') then
      beta_ref = 0.E0
      beta_type = 'ref'
      call gkw_warn('Option beta_type="eq" requires geom_type="chease".&
                      & Will use beta_ref=0 instead.')
  end if

  if(beta_type == 'eq' .and. beta_ref > r_tiny) then
     call gkw_warn('beta_ref is not used with beta_type == "eq"')
     beta_ref=0.E0
  end if

  if (geom_type /= 'chease' .and. betaprime_type == 'eq') then
      betaprime_ref = 0.E0
      betaprime_type = 'ref'
      call gkw_warn('Option betaprime_type="eq" requires geom_type="chease".&
                    & Will use betaprime_ref=0 instead.')
  end if

  if (.not. flux_tube .and. geom_type == 'miller') then
      if (gradp_type == 'beta_prime' .and. betaprime_type == 'sp') then
          call gkw_exit('geom_check: use betaprime_type=file or gradp_type=betaprime_input &
                       &to specify profiles for beta_prime in global runs with miller geometry')
      endif
  endif

  if (adiabatic_electrons .and. (beta_type == 'eq' .or.                    &
     & betaprime_type == 'sp')) then
      call gkw_warn('Be aware that the parameters of the adiabatic species&
                   & are used to calculate beta or beta_prime.')
  end if

  if (tearingmode .and. tm_drive)then
      call gkw_warn('tearingmode is for static islands  &
            &  and is incompatible with tm_drive to drive an island')
  end if

  if (.not.nlapar .and. tm_drive)then
      call gkw_exit('Need electromagnetic and tm_drive to drive an island')
  end if

  if(tm_drive.and.(abs(beta_ref).lt.r_tiny))then
      call gkw_exit('Beta needs to be finite for tm_drive to work')
  end if

  if (adiabatic_electrons .and. nlapar) then
      call gkw_warn('Electromagnetic effects with adiabatic electrons make &
                   & no sense')
  end if
  
  if (adiabatic_electrons .and. nlbpar) then
      call gkw_warn('Cannot use adiabatic electrons with compressional&
                   & magnetic field')
  end if
  
  !Zonal_adiabatic is meaningless for kinetic eletrons
  if (.not.adiabatic_electrons .and. zonal_adiabatic) then
    call gkw_warn('spcgeneral: zonal_adiabatic always false for kinetic&
                 & electrons')
    zonal_adiabatic=.false.
  endif

  !Legacy use
  if (finit=='island') then
    finit='cosine2'
    tearingmode=.true.
    call gkw_warn('finit="island" option is deprecated: &
               & Instead set logical tearingmode=T in spcgeneral namelist')
    call gkw_warn('For better convergence use the initialisation "zero"')
  end if
  
  if (tear_zero_epar.and..not.tearingmode) then
    call gkw_warn('tear_zero_epar requires tearingmode')
    tear_zero_epar=.false.
  end if

  !Ignore tiny beta except for tearing modes.
  if((beta_type=='ref'.and.beta_ref < r_tiny).and.(nlapar.or.nlbpar).and. &
    & (.not.tearingmode)) then
    call gkw_warn('spcgeneral: beta_ref=0 input, this run is electrostatic')
    nlapar=.false.
    nlbpar = .false.
  end if

  !Tearing modes must be electromagnetic
  if (tearingmode.and..not.nlapar) then
    call gkw_warn('tearing modes must be electromagnetic')
    nlapar=.true.
  endif

  ! Warn if beta is not zero for electrostatic runs
  if(beta_ref > r_tiny .and..not.nlapar .and..not.nlbpar) then
    if (beta_type == 'ref') then
      call gkw_warn('spcgeneral_check: input value of beta_ref ignored for&
                 & this electrostatic run')
      !beta_ref=0.
    end if
  end if
  
  if (abs(isl_shear) > r_tiny .and. abs(shat) > r_tiny) then
    call gkw_warn('isl_shear not used when shat /= 0')
    isl_shear = 0.
  end if  
  
  if (tm_sat < tm_start + 0.001) then
    call gkw_warn('tm_sat must be strictly greater than tm_start')
    tm_sat = tm_start + 0.001
  end if
  
  ! Ideally here should report on beta and betaprime and set for input.out
  ! but for now they only written to screen by linear terms because 
  ! they might be overwritten by geom at a later point.

  n_spec = number_of_species
  if (adiabatic_electrons) then
     iadia=1
     n_spec = n_spec + 1
  else
     iadia=0
  end if
  
  if (icrh_params(1) < 0) then
    call gkw_warn('Aniso input density < 0 transform switch deprecated')
    icrh_params(1) = abs(icrh_params(1))
  end if

  ! before calling the geometry, or dist_init routines, set lx for the
  ! global runs
  if (.not. flux_tube) then
    if (rhostar < r_tiny) call gkw_exit('Run is global, but rhostar is set to 0.')
    lx =  (psih - psil) / rhostar
  endif
  
  ! initialization with kinetic electrons
  if (finit=='cosine2' .and. .not. adiabatic_electrons) then
    call gkw_warn('spcgeneral: When adiabatic_electrons=.false. the&
    & initialization option finit="cosine2" (default) may lead to&
    & initial bursts. The option finit="cosine5" may be used&
    & instead.')
  end if
  

  ! more intuitive aliases for icrh inputs (could be removed)
  icrh_fmin   = icrh_params(1)  ! minority fraction  ~ <nm>/<nref>
  tperp_tpar  = abs(icrh_params(2))  ! Tp/T||
  dtperp_tpar = icrh_params(3)  ! -d (Tp/T||) / d psi
  
end subroutine components_check_params_spcg


!------------------------------------------------------------------------------
!> Read *one* species namelist.  Called by one processor
!> only for checking purposes.  The data is read a second time 
!> in components_input_species
!------------------------------------------------------------------------------
subroutine components_read_nml_spec(file_unit, io_stat)

  integer, intent(in)  :: file_unit
  integer, intent(out) :: io_stat

  ! namelist /species/ mass, z, temp, dens, rlt, rln, uprim, background, &
  !                    & param, dens_prof_type, dens_prof_coef,          &
  !                    & temp_prof_type, temp_prof_coef, rlt_gauss 

  ! Set unreasonable defaults for key species parameters
  mass=spbad; dens=spbad; z=spbad; rln=spbad; rlt=spbad
  
  ! default uprim is zero  
  uprim=0.0

  ! default temp  
  temp=1.0
  
  ! initialize the type to be a Maxwellian 
  background = 'maxwell'
  
  ! intialize the profile types
  dens_prof_type = 'none'      
  temp_prof_type = 'none'
  
  ! a parameter that can be used to define the background 
  param(:) = 0.
  
  ! initialisation of the gausswise temperature gradient.
  ! If the list rlt_gauss consists only of zeros or is omitted,
  ! fall back to the standard: a const gradient provided by rlt
  rlt_gauss(:) = 0.

  ! read the namelist and fill the dummy values
  read (file_unit,NML=species,IOSTAT=io_stat)

end subroutine components_read_nml_spec

!------------------------------------------------------------------------------
!> Write *one* species namelist.
!------------------------------------------------------------------------------
subroutine components_write_nml_spec(file_unit, i_spec, io_stat)
  use global, only : int2char_zeros
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor

  integer, intent(in)  :: file_unit, i_spec
  integer, intent(out) :: io_stat

  ! namelist /species/ mass, z, temp, dens, rlt, rln, uprim, background, &
  !                    & param, dens_prof_type, dens_prof_coef,          &
  !                    & temp_prof_type, temp_prof_coef, rlt_gauss

  ! fill the dummy values in the namelist with the actual values of
  ! species i_spec
  mass = mas_G(i_spec)
  z = signz_G(i_spec)
  dens = dgrid_G(i_spec)
  temp = tmp_G(1,i_spec)
  rlt = tp_G(1, i_spec)
  rln = fp_G(1, i_spec)
  uprim = vp_G(1, i_spec)
  background = types_G(i_spec)
  param(1) = pbg_G(1, i_spec)
  !FIXME param(2) is read from the namelist, but param(1) is used.. really??
  param(2) = 0
  dens_prof_type = de_prof_type(i_spec)
  dens_prof_coef = de_prof_coef(:, i_spec)
  temp_prof_type = tmp_prof_type(i_spec)
  temp_prof_coef = tmp_prof_coef(:, i_spec)
  !FIXME have to create an array to save rlt_gauss, if not it is lost
  !after reading
  rlt_gauss = rlt_gauss_params(:, i_spec)

  ! write the namelist to a file
  if(root_processor) write (file_unit, NML=species, IOSTAT=io_stat)

  ! write metadata
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'mass', mass)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'z', z)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'dens', dens)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'temp', temp)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'rlt', rlt)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'rln', rln)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'uprim', uprim)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), &
     & 'background', background)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), 'param', param)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), &
     & 'dens_prof_type', dens_prof_type)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), &
     & 'dens_prof_coef', dens_prof_coef)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), &
     & 'temp_prof_type', temp_prof_type)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), &
     & 'temp_prof_coef', temp_prof_coef)
  call write_run_parameter('species'//trim(int2char_zeros(i_spec,2)), &
     & 'rlt_gauss', rlt_gauss)

end subroutine components_write_nml_spec


!------------------------------------------------------------------------------
!> broadcast 1 species namelist
!------------------------------------------------------------------------------

subroutine components_bcast_nml_spec

  use mpiinterface, only : mpibcast 

  call mpibcast(mass,           1        )
  call mpibcast(temp,           1        )
  call mpibcast(dens,           1        )
  call mpibcast(z,              1        )
  call mpibcast(rln,            1        )
  call mpibcast(rlt,            1        )
  call mpibcast(uprim,          1        )
  call mpibcast(background,     7        )
  call mpibcast(param,          2        )
  call mpibcast(dens_prof_type, lenswitch)
  call mpibcast(temp_prof_type, lenswitch) 
  call mpibcast(dens_prof_coef, 5        )
  call mpibcast(temp_prof_coef, 5        )
  call mpibcast(rlt_gauss,      2        )
  call mpibcast(vpar_mean,      1        )

end subroutine components_bcast_nml_spec


!------------------------------------------------------------------------------
!> check species namelist; called until all species are read
!> called by all processors after all species data has been broadcast
!------------------------------------------------------------------------------

subroutine components_check_params_spec

  use general,      only : gkw_exit
  use grid,         only : number_of_species
  use global,       only : r_tiny
  use mpiinterface, only : root_processor
  use control,      only : flux_tube

  real, save    :: ntot = 0., nqas = 0., ngrad = 0., ne=0.
  integer, save :: call_count = 0
  
  if (call_count == 0) then

    nsps = number_of_species
    if (adiabatic_electrons) nsps = nsps + 1

    ! set the parameters for the check on quasi-neutrality 
    ntot  = 0.
    nqas  = 0.
    ngrad = 0.
    zeff_sp= 0.
    
  end if

  max_rho = max(max_rho,sqrt(2.*temp*mass)/z)

  call_count = call_count + 1
  if (call_count > nsps) then
    call gkw_exit('components_check_species called too many times')
  end if

  ! Check species inputs have been set (correctly)
  if (mass < 0.)          call gkw_exit('Species: Undefined mass forbidden')
  if (z <= spbad+r_tiny)  call gkw_exit('Species: Unreasonable charge Z')
  
  if (flux_tube) then
    if (temp < r_tiny)    call gkw_exit('Species: Undefined temp forbidden')
    if (dens < 0.)        call gkw_exit('Species: Undefined dens forbidden')
    if (rln <= spbad+r_tiny) call gkw_exit('Species: Unreasonable RLN')
    if (rlt <= spbad+r_tiny) call gkw_exit('Species: Unreasonable RLT')
    if (uprim <= spbad+r_tiny) call gkw_exit('Species: Unreasonable uprim')

    ! To check quasi-neutrality 
    ntot  = ntot  +   dens
    nqas  = nqas  + z*dens
    ngrad = ngrad + z*dens*rln

    ! calculate z_eff
    if (z > 0) then 
      zeff_sp = zeff_sp + dens*z*z
    else !electrons
      ne = ne + dens
    endif

    ! M L Reinke et al, Plasma Phys. Control. Fusion 54 045004 (2012) 
    ! Store parameters for the ICRH heated minority distribution function
    ! in the case in which it is treated as a separate species      
    if (trim(background) == 'ICRH') then 
      if (isg_icrh /=0) call gkw_exit('only one ICRH species currently allowed')
      isg_icrh = call_count
      ! For now the species is modelled as a Maxwellian, except in
      ! cf_quasineutral.
      ! For GK stability effects should be treated as
      ! a bi-Maxwellian background='Bi-maxwell'
      background='maxwell'
      icrh_params(1) = dens
      icrh_params(4) = z
    end if

    last_species : if (call_count == nsps) then
      ! Stop if quasi-neutrality is not satisfied
      if (abs(nqas) > 1e-6*abs(ntot)) then
        write(*,*)'sum z_s * n_s = ',nqas, ' number of species = ',nsps
        call gkw_exit('You do not satisfy quasi-neutrality')
      end if
      if (abs(ngrad) > 1e-6*abs(ntot)) then
        write(*,*)'sum z_s * n_s * dn_s/dr = ',ngrad, ' number of species = ',nsps
        call gkw_exit('You do not satisfy quasi-neutrality for the gradients')
      end if

      ! output Z_eff to screen
      zeff_sp = zeff_sp / ne
      if (nsps>2 .and. root_processor) then 
        write(*,*)
        write(*,*) 'Z_eff (of species inputs) = ', zeff_sp
        write(*,*)
      endif

      ! ICRH params (see manual)
      ! M L Reinke et al, Plasma Phys. Control. Fusion 54 045004 (2012) 

      icrh_model = 0
      if (isg_icrh > 0 .or. icrh_params(1) > 0.0) then
        if (isg_icrh == 0 .and. root_processor) then
           write(*,*)
           write(*,*) 'ANISO COMBINATION METHOD:' 
           write(*,*) 'Anisotropic minority added to species 1'        
        end if
   
        if (icrh_params(2) < 0.0) then
          icrh_model = 2
          if (root_processor) then
            write(*,*)
            write(*,*) 'Negative T_||/T_Perp: Bilato anisotropy model selected'
            write(*,*)
          end if
          icrh_params(2) = abs(icrh_params(2))
        else if (icrh_params(2) > 0.0) then
          icrh_model = 1
          if (root_processor) then
            write(*,*)
            write(*,*) 'Reinke anisotropy model selected'
            write(*,*)
          end if
        end if
      end if

    end if last_species
  
  end if ! flux_tube

end subroutine components_check_params_spec


!------------------------------------------------------------------------------
!> This function returns centrifugal quasineutrality test function 
!> Q(Phi(psi,theta)) = 
!>             Sum_s [(Z_s exp(-Z_S*phi/T_s)*exp(m_s*vcor^2*jfun/(T_s))]
!> on the flux surface with radial coordinate psi+dpsi, and also includes
!> the density variation from an (ICRH) species with anisotropic pressure.
!>
!> Note the factor (1/2) is lost in the normalisation
!> Called from rotation to intialise cfphi via bisection root finding method
!>
!> Since the centrifugal force can only be used in the local 
!> model, this function returns the value calculated at ix = 1
!>
!> Also implements the ICRH anisotropy driven poloidal asymmetry,
!> If this is used, the bulk ions must be the first species, and the electrons
!> should be the second species. Assumption: T_|| = T_m input
!>
!> The checks and setup required for
!> this routine are performed in rotation: centrifugal_energy
!------------------------------------------------------------------------------

function cf_quasineutral(phi,vcor2jfun,bn_bR0,dpsi)

  use general, only : gkw_abort 

  real, intent(in) :: phi, dpsi, bn_bR0   !Last one: B / B_R0
  real, intent(in) :: vcor2jfun
  real             :: cf_quasineutral

  real    :: dum, func, temp, dens, eta, tp_tplfs
  integer :: is

  !calculate the function Q(Phi(psi,theta))
  func = 0.E0
 
  do is=1,nsps

    eta = 0.0     
    temp = tmp_G(1,is)*(1.E0-tp_G(1,is)*dpsi)
    dens = de_G(1,is)*(1.E0-fp_G(1,is)*dpsi)

    ! Anisotropic minority parameters
    ! Independent method (own species):
    if (icrh_model > 0 .and. isg_icrh == is) then
      if (icrh_model == 1) then      ! Reinke
          eta = tperp_tpar+dtperp_tpar*dpsi- 1.0  
      else if (icrh_model == 2) then ! Bilato
          dum = tperp_tpar+dtperp_tpar*dpsi
          tp_tplfs= 1.0/(dum + (1.0 - dum)/bn_bR0)
          dens = dens*tp_tplfs
      end if
    ! Combination method (species 1 contains the minority): 
    ! Correct main ion density 
    else if (is == 1 .and. icrh_model > 0 .and. isg_icrh == 0) then
      call gkw_abort('ICRH species combination method removed')
      !dens = (dens-icrh_fmin*icrh_params(4)/signz_G(1))*(1.E0-fp_G(1,is)*dpsi)   
    end if 
    
    dum  = (vcor2jfun*mas_G(is)-signz_G(is)*phi)/temp

    ! Protect against floating point overflow: If triggered, the starting limits
    ! for the bisection algorithm probably need adjusting.
    ! If you are using single precision, try double precision.
    if (abs(dum) > maxexponent(dum)*0.8 ) call gkw_abort('cf_quasineutral & 
            & overflow due to extreme species data')

    func=func+signz_G(is)*dens*(bn_bR0)**(-eta)*exp(dum)

  end do

  ! Alternative: Add the anisotropic minority as an additonal extra species 
  ! Centrifugal effects on minority species assumed equal to bulk ion (is = 1)
  ! and T|| assumed to be equal to T_i
  if (isg_icrh == 0 .and. icrh_model > 0) then
      call gkw_abort('ICRH species combination method removed')
!   
!     is = 1
!     eta = 0.0
!     dens = icrh_fmin*(1.E0-icrh_rln*dpsi)
!     temp = tmp_G(1,is)*(1.E0-tp_G(1,is)*dpsi)
! 
!     if (icrh_model == 1) then      ! Reinke
!         eta = tperp_tpar+dtperp_tpar*dpsi- 1.0  
!     else if (icrh_model == 2) then ! Bilato
!         dum = tperp_tpar+dtperp_tpar*dpsi
!         tp_tplfs= 1.0/(dum + (1.0 - dum)/bn_bR0)
!         dens = dens*tp_tplfs
!     end if
! 
!     dum  = (vcor2jfun*mas_G(is)-icrh_params(4)*phi)/temp
!     
!     if (abs(dum) > maxexponent(dum)*0.8 ) call gkw_abort('cf_quasineutral & 
!             & overflow due to extreme species data')
! 
!     func=func+icrh_params(4)*dens*(bn_bR0)**(-eta)*exp(dum)
! 
   end if
  
  cf_quasineutral=func
  
  ! ICRH: WORK IN PROGRESS
  ! DONE: Fix the simple method with non-trace imp
  ! DONE: Fix the alternative method with non-trace imp
  ! DONE: Check how it works correctly with strong rotation / R0loc, miller geom
  ! DONE: check sensitivity to fmin (linear), zmin (linear), tperp_tpar (linear in eta), tpar (none)
  ! DONE: remove T|| and dT|| inputs: neglible impact even at large fmin and tpar
  ! DONE: Fix and test for Z_m and Z_i /= 1
  ! DONE: Document
  ! DONE: Add input checks: check that the species are in the right places.
  ! DONE: Some sanity checks on the radial derivatives and sensitivities
  ! DONE: Fix and check radial derivatives diagnostics
  ! DONE: Add CLA memo quantities, and use to check FSA gradients quasineutrality
  ! DONE: Update docs and memo on radial derivatives
  ! DONE: Test in parallel, add a test case (or two), 
  ! DONE: Calculate QN RLN for Bilato model and check against analtic (update docs)
  ! TO DO: test with NX>1
  ! TO DO: Fix for R0_loc=axis ? 
  ! TO DO: retest with temperatures /= 1.0, m /= 1
  ! STARTED: Comment, cleanup, make subroutine in rotation

  ! WONTFIX Test and fix the combined method for both models (but not input transformations)
  ! Remove dual usage of icrh_params and tperp_tpar etc 
  ! SIMPLIFY THE input parameters if combined method is to be permanently removed
  
  ! STARTED: adjust rln of inputs to keep constant minority fraction
  ! WONTFIX: the above for the combined species method


end function cf_quasineutral 


!------------------------------------------------------------------------------
!> This function returns [TeTi/(Te+Ti)](mi/Ti-me/Te)
!>                        =(TeMi-TiMe)/(Te+Ti)
!> Note the factor (1/2) is lost in the normalisation
!> If no singly charged ions, returns 0.
!> Uses dominant density single charged species as the bulk ions.
!> Only used to calculate the analytic comparison for pure hydrogenic plasmas
!------------------------------------------------------------------------------

function cf_mass_weight()
  use global, only : r_tiny
  use general, only : gkw_warn 

  real :: cf_mass_weight
  
  integer :: is, iel, iion, isw = 0, isw2 = 0
  real :: de_max, dum
  
  !Which species are electrons and ions?
  iel  = 0
  iion = 0

  !Max ion density, to find bulk ions.
  de_max = 0.E0
  
  !Find the species we want to use
  do is=1,nsps

    if (abs(signz_G(is) - (-1)) < r_tiny) then
      if (iel /= 0 .and. isw == 0) call gkw_warn('cf_mass_weight not setup for >1 electron species')
      iel = is
      isw = 1
    else 
      if (abs(signz_G(is) - 1) < r_tiny) then
        if (de_G(1,is).gt.de_max) then 
          de_max=de_G(1,is)
          iion=is
        end if
      else !Do nothing
      end if
    endif 
  end do

  if (iel*iion ==0) then
    if (isw2 == 0) call gkw_warn('cf_mass_weight: No singly charged ions')
    cf_mass_weight=0.0
    isw2 = 1
    return
  else  
    dum=tmp_G(1,iel)*tmp_G(1,iion)/(tmp_G(1,iel)+tmp_G(1,iion))
    cf_mass_weight=dum*(mas_G(iion)/tmp_G(1,iion)-mas_G(iel)/tmp_G(1,iel))
    return
  end if 

end function cf_mass_weight


!------------------------------------------------------------------------------
!> This function returns complicated radial derivative thing
!> for the singly charged ions case
!> If no singly charged ions, returns 0.
!> Uses dominant density single charged species as the bulk ions.
!> Only used to calculate the analytic comparison for pure hydrogenic plasmas
!------------------------------------------------------------------------------

function cf_dtf1()
  use global, only : r_tiny
  use general, only : gkw_warn 

  real :: cf_dtf1
  
  integer :: is, iel, iion, isw = 0, isw2 = 0
  real    :: de_max, dum, bot, dti, dte
  
  cf_dtf1 = 0.E0

  !Which species are electrons and ions?
  iel  = 0
  iion = 0

  !Max ion density, to find bulk ions.
  de_max = 0.E0
  
  !Find the species we want to use
  do is=1,nsps

    if (abs(signz_G(is) - (-1)) < r_tiny) then
      if (iel /= 0 .and. isw == 0) call gkw_warn('cf_dtf1 not setup for >1 electron species') 
      iel = is
      isw = 1
    else 
      if (abs(signz_G(is) - 1) < r_tiny) then
        if (de_G(1,is).gt.de_max) then 
          de_max=de_G(1,is)
          iion=is
         end if
      else !Do nothing
      endif 
    end if

  end do

  if (iel*iion ==0) then
    if (isw2 == 0) call gkw_warn('cf_dtf1: No singly charged ions')
    cf_dtf1=92.0
    isw2 = 1
    return
  else  
    dti=-tmp_G(1,iion)*tp_G(1,iion)
    dte=-tmp_G(1,iel)*tp_G(1,iel)
    bot=tmp_G(1,iion)+tmp_G(1,iel)
    cf_dtf1=(mas_G(iion)*dte-mas_G(iel)*dti)/bot
    dum=(mas_G(iel)*tmp_G(1,iion)-mas_G(iion)*tmp_G(1,iel))/(bot*bot)
    cf_dtf1=cf_dtf1+dum*(dti+dte)
    return
  end if 

end function cf_dtf1


!------------------------------------------------------------------------------
!> This subroutine allocates the arrays of the module components 
!>
!> mas(nsp+iadia)        the normalised mass 
!> signz(nsp+iadia)      the charge number
!> types(nsp+iaia)       the type of background  
!> tmp(nx,nx,nsp+iadia)  the normalised temperature 
!> de(nx,nx,nsp+iadia)   the normalised density 
!> tp(nx,nx,nsp+iadai)   the temperature gradient 
!> fp(nx,nsp+iadia)      the density gradient 
!> vp(nx,nsp+iadai)      the parallel velocity gradient
!> vthrat(nx,nsp)        the ratio of the thermal / reference velocity
! 
!------------------------------------------------------------------------------

subroutine components_allocate

  use grid,    only : nsp, nx, n_x_grid
  use general, only : gkw_abort

  integer ierr

  if (adiabatic_electrons) then
    iadia=1
  else
    iadia=0
  end if

  ! initilize the error parameter
  ierr = 0

  ! allocate the mass array
  allocate(mas(nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate mas in components')

  ! allocate the charge array
  allocate(signz(nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate signz in components')

  ! allocate the type of the species 
  allocate(types(nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate vp in components')

  ! allocate the grid temperature array 
  allocate(tgrid(nsp+iadia),stat=ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate tgrid in components')

  ! allocate the grid density array 
  allocate(dgrid(nsp+iadia),stat=ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate ngrid in components')

  ! allocate the temperature array
  allocate(tmp(nx,nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate tmp in components')

  ! allocate the density array
  allocate(de(nx,nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate de in components')

  ! allocate the temperature gradient array
  allocate(tp(nx,nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate tp in components')

  ! allocate the density gradient array
  allocate(fp(nx,nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate fp in components')

  ! allocate the velocity gradient array
  allocate(vp(nx,nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate vp in components')

  ! allocate the parameter for the background array 
  allocate(pbg(nx,nsp+iadia),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate vp in components')

  ! the ratio of thermal and reference velocity 
  allocate(vthrat(nsp),stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate vthrat in components')
    
  ! the beta value 
  allocate(veta(nx), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate veta in components')

  ! the beta-prime value 
  allocate(veta_prime(nx), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate veta_prime in components')

  !GLOBAL SPECIES ARRAYS

  ! allocate the mass array
  allocate(mas_G(nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate mas_G in components')

  ! allocate the charge array
  allocate(signz_G(nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate signz_G in components')

  ! allocate
  allocate(types_G(nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate types_G in components')

  ! allocate
  allocate(pbg_G(nx,nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate pbg_G in components')

  ! allocate the grid temperature array 
  allocate(tgrid_G(nsps),stat=ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate tgrid_G in components')

  ! allocate the grid temperature array 
  allocate(dgrid_G(nsps),stat=ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate ngrid_G in components')

  ! allocate the temperature array
  allocate(tmp_G(nx,nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate tmp_G in components')

  ! allocate the density array
  allocate(de_G(nx,nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate de_G in components')

  ! allocate the velocity array
  allocate(vthrat_G(nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate vthrat_G in components')

  ! allocate the temperature gradient array
  allocate(tp_G(nx,nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate tp_G in components')

  ! allocate the density gradient array
  allocate(fp_G(nx,nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate fp_G in components')

  ! allocate
  allocate(vp_G(nx,nsps),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate up_G in components')

  ! allocate
  allocate(ecurx(nx),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate ecurx in components')

  ! The density profile type  
  allocate(de_prof_type(nsps), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Could not allocate de_prof_type in components')
  
  ! The temperature profile type  
  allocate(tmp_prof_type(nsps), stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate tmp_prof_type in components')

  ! The density profile coefficients 
  allocate(de_prof_coef(5,nsps), stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate de_prof_coef in components')

  ! The temperature profile coefficients 
  allocate(tmp_prof_coef(5,nsps), stat = ierr)  
  if (ierr /= 0) call gkw_abort('Could not allocate tmp_prof_coef in components')

  ! Parameters for a gaussian background temperature gradient
  allocate(rlt_gauss_params(2,nsps), stat = ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate rlt_gauss_params in components')

  ! Global arrays in the radial direction 

  allocate(veta_Gx(1:n_x_grid), stat = ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate veta_Gx in gyro-average')

  ! initialize veta_Gx to zero (ES case of the timestep estimator) 
  veta_Gx(:) = 0.E0 

  ! allocate the arrays for gathering in x (different from those in components)
  allocate(tmp_Gx(1:n_x_grid,nsp+iadia), stat = ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate tmp_Gx in gyro-average')

  allocate(de_Gx(1:n_x_grid,nsp+iadia), stat = ierr) 
  if (ierr /= 0) call gkw_abort('could not allocate de_Gx in gyro-average')

return
end subroutine components_allocate


!------------------------------------------------------------------------------
!> This subroutine reads the species data again (after the checks)
!> and puts the data into the arrays.
!> This routine must be called by all processors.
!>
!> Formerly, this routine read from input.out, instead of input.dat.
!>
!> SRG does not know why exactly and hopes that input.dat is also okay:
!>
!> Now, this routine reads from input.dat, because input.out is
!> written *after* this routine initializes the species arrays.
!> This is because metadata in other formats shall also  be written when
!> input.out is produced, and for this, the IO interfaces need to be
!> initialized. And those require some information from input.dat .
!------------------------------------------------------------------------------
subroutine components_input_species
  use global,       only : r_tiny
  use io,           only : get_free_file_unit
  use grid,         only : nsp, lsp, nx, gx, n_x_grid
  use general,      only : gkw_abort
  use mpiinterface, only : root_processor
  use geom,         only : mas_miller_i,tmp_miller,de_miller,fp_miller,tp_miller, &
                         & vp_miller_i 
  use mpicomms,     only : COMM_X_NE, COMM_DUMMY
  use mpiinterface, only : gather_array

  !> gather_array switch
  logical, save :: ALL_PROCS=.true. 
  
  integer :: i, file_unit, io_stat, it1, ix, itempo
  integer, save :: icheck = 0 
  real    :: dum
  
  io_stat        = 0
  it1            = 0
  itempo         = 1
  cf_phi_max(:)  = 0.
  
  if (root_processor) then
    !oldio: call get_free_lun(ilun)
    call get_free_file_unit(file_unit)
    open(UNIT=file_unit, FILE='input.dat', FORM='formatted', STATUS='old', &
       & POSITION='rewind')
  end if
  
  do i=1,nsps
    if (root_processor) then
      ! read the namelist and fill the dummy values
      call components_read_nml_spec(file_unit, io_stat)
    end if
    call components_bcast_nml_spec
    !Adiabatic species always come last.
    if ( z < 0 .and. adiabatic_electrons) then  

      if (iadia.eq.0) call gkw_abort('Severe error in components')
      if (icheck /= 0) call gkw_abort('Only one adiabatic species is allowed')

      mas(nsp+iadia)   = mass
      signz(nsp+iadia) = z
      types(nsp+iadia) = background
      tgrid(nsp+iadia) = temp 
      dgrid(nsp+iadia) = dens
      do ix = 1, nx 
        tmp(ix,nsp+iadia)   = temp
        de(ix,nsp+iadia)    = dens
        fp(ix,nsp+iadia)    = rln
        ! setup the background temperature gradient gaussianwise
        if(.not. (abs(rlt_gauss(1)) < r_tiny .and. &
           & abs(rlt_gauss(2)) < r_tiny)) then
            tp(ix,nsp+iadia) = rlt_gauss(1) * &
               & exp(-(gx(ix)-n_x_grid/2.0)**2/(2*rlt_gauss(2)**2))
        else
          ! temperature gradient is const. on the box
          tp(ix,nsp+iadia)    = rlt
        end if
        vp(ix,nsp+iadia)    = uprim
        pbg(ix,nsp+iadia)   = param(1)
      end do 

      !global arrays of data on ALL species
      mas_G(nsps)   = mass      
      signz_G(nsps) = z
      types_G(nsps) = background
      tgrid_G(nsps) = temp 
      dgrid_G(nsps) = dens 
      vthrat_G(nsps) = sqrt(tgrid_G(nsps)/mas_G(nsps))
      do ix = 1, nx 
        tmp_G(ix,nsps)        = temp
        de_G(ix,nsps)         = dens 
        vp_G(ix,nsps)         = uprim
        fp_G(ix,nsps)         = rln
        ! setup the background temperature gradient gaussianwise
        if(.not. (abs(rlt_gauss(1)) < r_tiny .and. &
           & abs(rlt_gauss(2)) < r_tiny )) then
            tp_G(ix,nsp) = rlt_gauss(1) * &
               & exp(-(gx(ix)-n_x_grid/2.0)**2/(2*rlt_gauss(2)**2))
        else
          ! temperature gradient is const. on the box
          tp_G(ix,nsps)         = rlt 
        endif 
        tmp_miller(ix,2)      = temp
        de_miller(ix,2)       = dens
        fp_miller(ix,2)       = rln
        tp_miller(ix,2)       = rlt
        pbg_G(ix,nsps) = param(1)
      end do 
      de_prof_type(nsps) = dens_prof_type
      tmp_prof_type(nsps) = temp_prof_type
      de_prof_coef(:,nsps) = dens_prof_coef
      tmp_prof_coef(:,nsps) = temp_prof_coef
      rlt_gauss_params(:, nsps) = rlt_gauss
      icheck = 1  

    else !kinetic species
      it1 = it1 + 1
      
      !Local to this processor      
      if ( lsp(it1) >= 1 .and.  lsp(it1) <= nsp ) then 
        mas(lsp(it1))   = mass
        signz(lsp(it1)) = z
        types(lsp(it1)) = background
        tgrid(lsp(it1)) = temp 
        dgrid(lsp(it1)) = dens
        rlt_gauss_params(:, lsp(it1)) = rlt_gauss
        do ix = 1, nx 
          tmp(ix,lsp(it1))   = temp
          de(ix,lsp(it1))    = dens
          fp(ix,lsp(it1))    = rln
          ! setup the background temperature gradient gaussianwise
          if(.not. (abs(rlt_gauss(1)) < r_tiny .and. &
             & abs(rlt_gauss(2)) < r_tiny)) then
            tp(ix,lsp(it1)) = rlt_gauss(1) * &
               & exp(-(gx(ix)-n_x_grid/2.0)**2/(2*rlt_gauss(2)**2))
          else
            ! temperature gradient is const. on the box
          tp(ix,lsp(it1))    = rlt
          endif 
          vp(ix,lsp(it1))    = uprim
          pbg(ix,lsp(it1))   = param(1)
        end do 

      end if !Species on this processor

      !Global arrays of data on ALL species (used for rotation and collisions)
      mas_G(it1)    = mass
      signz_G(it1)  = z
      types_G(it1) = background
      tgrid_G(it1)  = temp
      dgrid_G(it1)  = dens
      vthrat_G(it1) = sqrt(tgrid_G(it1)/mas_G(it1))
      do ix = 1, nx
        tmp_G(ix,it1)    = temp
        de_G(ix,it1)     = dens
        vp_G(ix,it1)     = uprim
        fp_G(ix,it1)     = rln
        ! setup the background temperature gradient gaussianwise
        if(.not. (abs(rlt_gauss(1)) < r_tiny .and. &
           & abs(rlt_gauss(2)) < r_tiny)) then
          tp_G(ix,it1) = rlt_gauss(1)*exp(-(gx(ix)-n_x_grid/2.0)**2/(2*rlt_gauss(2)**2))
        else
          ! temperature gradient is const. on the box
          tp_G(ix,it1)     = rlt
        endif
        !Only the main ion and electrons parameters are needed for rota_miller
        !in miller geometry (Two species model to account for the effect of toroidal
        !rotation on the magnetic equilibrium)
        !Find the main ion species
        if (z > 0 .and. it1 == 1) then
          vp_miller_i(ix)    = uprim
          tmp_miller(ix,1)   = temp
          de_miller(ix,1)    = dens
          tp_miller(ix,1)    = rlt
          fp_miller(ix,1)    = rln
        else if (z < 0 .and. it1 == 1) then
          de_miller(ix,1) = 0.
        end if

        if (it1 > 1 .and. z > 0 .and. dens > de_miller(ix,1)) then
          vp_miller_i(ix)    = uprim
          tmp_miller(ix,1)   = temp
          de_miller(ix,1)    = dens
          tp_miller(ix,1)    = rlt
          fp_miller(ix,1)    = rln
          itempo = it1
        end if

        if (z < 0) then
          tmp_miller(ix,2)   = temp
          de_miller(ix,2)    = dens
          tp_miller(ix,2)    = rlt
          fp_miller(ix,2)    = rln
        end if
        pbg_G(ix,it1)   = param(1)
      end do

      mas_miller_i = mas_G(itempo)

      de_prof_type(it1) = dens_prof_type
      tmp_prof_type(it1) = temp_prof_type
      de_prof_coef(:,it1) = dens_prof_coef
      tmp_prof_coef(:,it1) = temp_prof_coef
      rlt_gauss_params(:, it1) = rlt_gauss

    end if

    ! Estimate limits for centrifugal potential bisection
    if(abs(z*dens)>-1.0) cf_phi_max(3)=max(log(1.0+abs(z*dens)),cf_phi_max(3))
    if(dens > -1.0) cf_phi_max(4)=max(mass*log(1.0+dens)/temp,cf_phi_max(4))

  end do

  if (root_processor) close(file_unit)

  ! Calculate the normalizing coefficients 
  do i = 1, nsp; do ix = 1, nx 
    vthrat(i) = sqrt(tgrid(i)/mas(i))
  end do; end do   

  !Estimate limits for centrifugal potential bisection
  dum=abs(maxval(signz_G)*maxval(mas_G)/minval(tmp_G))
  cf_phi_max(1)=log(dum)
  cf_phi_max(2)=maxval(tmp_G)

  ! gather the arrays for global in the radial direction 
  
  call gather_array(de_Gx(1:n_x_grid,:), n_x_grid,  nsp+iadia,           &
                  & de   (1:nx,      :), nx,        nsp+iadia,           &
                  & COMM_X_NE, COMM_DUMMY, ALLGATHER = ALL_PROCS)
                  
  call gather_array(tmp_Gx(1:n_x_grid,:), n_x_grid, nsp+iadia,           &
                  & tmp   (1:nx,      :), nx,       nsp+iadia,           &
                  & COMM_X_NE, COMM_DUMMY, ALLGATHER = ALL_PROCS) 

  
end subroutine components_input_species


!------------------------------------------------------------------------------
!> Select beta and beta_prime according to the namelist switches
!------------------------------------------------------------------------------
subroutine components_beta_cal 

  use control,      only : flux_tube, nlapar
  use mpicomms,     only : COMM_X_NE
  use mpiinterface, only : gather_array
  use grid,         only : nx, n_x_grid, gx 
  use geom,         only : beta_eq, betaprime_eq, betaprime_miller, beta_miller
  use general,      only : gkw_abort, gkw_warn 
  use mpiinterface, only : root_processor, mpibcast
  use io,           only : get_free_file_unit 
  

  logical, save :: ALL_PROCS=.true.

  real, allocatable :: vetadum(:),vetaprimedum(:)
  real    :: dum, xdum
  character(200) :: line 
  integer :: ix, is, file_unit, ierr

  ! initialize beta and beta_prime
  veta       = 0.
  veta_prime = 0.

  ! read from file if needed
  if (beta_type == 'file' .or. betaprime_type == 'file') then 
    
    ierr = 0 
    allocate(vetadum(n_x_grid), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate vetadum') 
    allocate(vetaprimedum(n_x_grid), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate vetaprimedum') 

    if(root_processor)then
      call get_free_file_unit(file_unit)
      open(file_unit,file='input.prof',status="old")
      line = ''
      do while (line(1:11) /= '#Spcgeneral')  
        read(file_unit,'(A)',ERR = 100,END=100) line 
      end do 
      do ix = 1, n_x_grid 
        read(file_unit,*)xdum,vetadum(ix),vetaprimedum(ix)
        !if (abs(xdum-xgr(ix)).gt.1e-5) call gkw_abort('Inconsistent xgr in input.prof')
        !can't test the radial grid since xgr is filled in yet
      end do 
      close(file_unit) 
    endif 
 
    call mpibcast(vetadum,    n_x_grid) 
    call mpibcast(vetaprimedum, n_x_grid)
    
  endif

  ! Fill in beta
  select case(beta_type)

    case('ref') ! user specified
      do ix = 1, nx 
        veta(ix) = beta_ref
      end do 

    case('file') ! read from file
      do ix = 1, nx 
        veta(ix) = vetadum(gx(ix))
      end do 

    case('eq') ! from actual MHD equilibrium
      call gkw_warn('WARNING: Beta_eq has not yet been made global')
      do ix = 1, nx 
        dum =0.
        do is = 1, nsps
          dum = dum + de_G(ix,is)*tmp_G(ix,is)
        end do
        veta(ix) = beta_eq / dum
      end do 

    case default
      call gkw_abort('Linear_terms: unknown beta_type option')

  end select


  ! Fill in beta_prime
  select case(betaprime_type)

    case('ref') ! user specified
      do ix = 1, nx 
        veta_prime(ix) = betaprime_ref
      end do 

    case('file') ! read from file
      do ix = 1, nx 
        veta_prime(ix) = vetaprimedum(gx(ix))
      end do 

    case('sp') ! consistent with species parameters and beta
      do ix = 1, nx 
        dum =0.
        do is = 1, nsps
          dum = dum - de_G(ix,is)*tmp_G(ix,is)*(tp_G(ix,is)+fp_G(ix,is))
        end do
        veta_prime(ix) = veta(ix) * dum
      end do
      if (veta(1) < 5e-4) call gkw_warn('beta_prime using small beta value')

    case('eq') ! from actual MHD equilibrium
      call gkw_warn('betaprime_eq has not been made global')
      do ix = 1, nx 
        veta_prime(ix) = betaprime_eq
      end do 

    case default
      call gkw_abort('Linear_terms: unknown betaprime_type option')

  end select
 
  ! Fill the global (in x) arrays of beta_miller and betaprime_miller
  call gather_array(beta_miller(1:n_x_grid), n_x_grid,                &
                  & veta(1:nx)             , nx,                      &
                  & COMM_X_NE, ALLGATHER = ALL_PROCS)   
  call gather_array(betaprime_miller(1:n_x_grid), n_x_grid,           &
                  & veta_prime(1:nx)            , nx,                 &
                  & COMM_X_NE, ALLGATHER = ALL_PROCS)   

  ! Report values used if different from refs
  ! ideally these values should be written to input.out, but this is
  ! not possible until geom reads are moved.
  if (flux_tube) then 
    if (beta_type /= 'ref' .and. root_processor) then
      write(*,*) 
      write(*,*) '***** Using beta = ', veta(1)
      write(*,*)
    end if 
  
    if (betaprime_type /= 'ref' .and. root_processor) then
      write(*,*)
      write(*,*) '***** Using beta_prime = ', veta_prime(1)
      write(*,*)
    end if  
  else 
  endif 

  ! Fill the global array (in x) of veta
  if (nlapar) call gather_array(veta_Gx(1:n_x_grid), n_x_grid,        &
                  &             veta   (1:nx      ), nx,              &
                  & COMM_X_NE, ALLGATHER = ALL_PROCS)

  return 
  100 call gkw_abort('Error while reading file input.prof in components_beta_cal') 
   
end subroutine components_beta_cal


!------------------------------------------------------------------------------
!> This subroutine initialises the radial background profiles of density and
!> temperature.
!>
!> This must be called after parallelize_geom.
!>
!> WARNING: a single local uprim is still used in the global version and
!> there is no rotation profile implemented.
!------------------------------------------------------------------------------
subroutine components_profiles
  use global,       only : r_tiny
  use control,      only : flux_tube 
  use grid,         only : nsp, nx, number_of_species, gsp, nsp, n_x_grid, gx 
  use geom,         only : geom_parallelized
  use general,      only : gkw_abort
  use mpicomms,     only : COMM_X_NE, COMM_DUMMY
  use mpiinterface, only : gather_array, root_processor, mpibcast
  use io,           only : get_free_file_unit
  logical :: readfile
  logical, allocatable :: read_prof(:,:) 
  real, allocatable :: ndum(:),tdum(:),fpdum(:),tpdum(:)
  
  integer        :: is, ix, file_unit, ierr, iss, igsp
  real           :: val_norm, ntot, nqas, ngrad, dndum, dtdum, xdum 
  real           :: mdum, zdum 
  character(200) :: line 
  
  ! background profiles are only calculated for the global version of the code. 
  if (flux_tube) return 

  if (.not. geom_parallelized) call gkw_abort('components_profiles too early')

  ! intialize the profiles to zero 
  de_G  = 0 
  tmp_G = 0 
  fp_G  = 0 
  tp_G  = 0 
  
  ! Check if any of the profiles needs to be read from file
  readfile = .false.
  do is = 1, number_of_species+iadia 
    if (de_prof_type(is)  == 'file') readfile = .true.
    if (tmp_prof_type(is) == 'file') readfile = .true.
  end do 
  
  ! Read the profiles from file if necessary 
  if (readfile) then 
    if (root_processor) then 
      call get_free_file_unit(file_unit) 
      open(file_unit,file='input.prof',ERR = 100,status="old")
      line = '' 
    endif 
 
    ierr = 0 
    allocate(read_prof(number_of_species+iadia,2), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Cound not allocate read_prof')
    allocate(ndum(n_x_grid), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate ndum') 
    allocate(tdum(n_x_grid), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate tdum') 
    allocate(fpdum(n_x_grid), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could no allocate fpdum') 
    allocate(tpdum(n_x_grid), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate tpdum') 
    
    do is = 1, number_of_species+iadia 
      read_prof(is,1) = de_prof_type(is)  == 'file'
      read_prof(is,2) = tmp_prof_type(is) == 'file'
    end do
    ndum  = 0. 
    tdum  = 0. 
    fpdum = 0. 
    tpdum = 0.

    do while (readfile) 
    
      ! read the data 
      if (root_processor) then 
        do while (line(1:8) /= '#Species')
          read(file_unit,'(A)',ERR=100,END=100)line 
        end do 
        read(file_unit,*)mdum,zdum,dndum,dtdum 
        read(file_unit,'(A)')line 
 
        ! find the species with the correct mass and charge 
        iss = 0 
        do is = 1, number_of_species+iadia 
          if (abs(mas_G(is)-mdum) < r_tiny .and. &
             & abs(signz_G(is)-zdum) < r_tiny) then
            if (read_prof(is,1) .or. read_prof(is,2)) then 
              iss = is 
            endif 
          endif 
        end do 
        do ix = 1, n_x_grid 
          read(file_unit,*,ERR = 100)xdum,ndum(ix),tdum(ix),fpdum(ix),tpdum(ix) 
        end do 
      endif

      ! broadcast to all processors 
      call mpibcast(iss,         1) 
      call mpibcast(dndum,       1) 
      call mpibcast(dtdum,       1) 
      call mpibcast(ndum, n_x_grid) 
      call mpibcast(tdum, n_x_grid) 
      call mpibcast(fpdum,n_x_grid) 
      call mpibcast(tpdum,n_x_grid) 
  
      
      ! if iss == 0 then the corresponding data in input.prof is not found 
      ! in the input file. This is allowed: one can have a input.prof file 
      ! with many species and run GKW only with a few of these species 
      if (iss /= 0) then 
      
        ! distribute when parallelizing in the radial direction 
        if (read_prof(iss,1)) then 
          read_prof(iss,1) = .false. 
          do ix = 1, nx 
            de_G(ix,iss) = ndum(gx(ix)) 
            fp_G(ix,iss) = fpdum(gx(ix)) 
          end do 
          dgrid_G(iss) = dndum 
        endif 
        if (read_prof(iss,2)) then 
          read_prof(iss,2) = .false. 
          do ix = 1, nx 
            tmp_G(ix,iss) = tdum(gx(ix)) 
            tp_G(ix,iss) = tpdum(gx(ix)) 
          end do 
          tgrid_G(iss) = dtdum 
          vthrat_G(iss) = sqrt(tgrid_G(iss)/mas_G(iss)) 
        endif 
      
      endif 
      
      ! more to read ? 
      readfile = .false. 
      do is = 1, number_of_species+iadia 
        readfile = (readfile .or. read_prof(is,1)) 
        readfile = (readfile .or. read_prof(is,2)) 
      end do 
      
    end do

    if (root_processor) then
      close(file_unit)
    end if
  end if 
  
  ! Calculate the analytic global profiles 
  do is = 1, number_of_species + iadia 
  
    if (de_prof_type(is) /= 'file')  then 
    
      ! density profile 
      call call_prof(de_G(1,is),fp_G(1,is),de_prof_type(is),de_prof_coef(1,is),val_norm)

      ! Grid density value 
      dgrid_G(is) = val_norm 

    endif 
    
    if (tmp_prof_type(is) /= 'file') then 
    
      ! temperature profile 
      call call_prof(tmp_G(1,is),tp_G(1,is),tmp_prof_type(is),tmp_prof_coef(1,is),val_norm) 

      ! Grid temperature value 
      tgrid_G(is) = val_norm

      ! Set the normalizing velocity 
      vthrat_G(is) = sqrt(tgrid_G(is)/mas_G(is))  

    endif 
    
  end do

  !If tearing mode is driven from background current profile
  !This must be calculated from the q profile  
  if(tm_drive)then
    !if(root_processor)write(*,*)'Current profile calc'
    call calc_ecur_prof
  endif

  ! To check quasi-neutrality
  do ix = 1, nx 
    ntot = 0; nqas = 0.; ngrad = 0; 
    do is = 1, number_of_species + iadia 
      ntot  = ntot  +  de_G(ix,is) 
      nqas  = nqas  +  signz_G(is)*de_G(ix,is) 
      ngrad = ngrad +  signz_G(is)*de_G(ix,is)*fp_G(ix,is)
    end do 

    ! Stop if quasi-neutrality is not satisfied
    ! gkw_abort is safer here in case of parallel_x
    if (abs(nqas) > 1e-6*abs(ntot)) then
      write(*,*)'sum z_s * n_s = ',nqas, ' (number of species+adiabatic: ', &
             & number_of_species + iadia, ')'
      call gkw_abort('You do not satisfy quasi-neutrality')
    end if
    if (abs(ngrad) > 1e-6*abs(ntot)) then
      write(*,*)'sum z_s * n_s * dn_s/dr= ',ngrad
      call gkw_abort('You do not satisfy quasi-neutrality for the gradients')
    end if
  end do

  ! copy values from global to local array 
  do is = 1, nsp
    igsp = gsp(is)

    dgrid(is)  = dgrid_G(igsp)
    tgrid(is)  = tgrid_G(igsp)
    vthrat(is) = vthrat_G(igsp)
    do ix = 1, nx
      de(ix,is)  = de_G(ix,igsp)
      tmp(ix,is) = tmp_G(ix,igsp)
      tp(ix,is)  = tp_G(ix,igsp)
      fp(ix,is)  = fp_G(ix,igsp)
      vp(ix,is)  = vp_G(ix,igsp)
    end do
  end do

  ! Adiabatic species
  ! Special treatment as this is present on all processors.
  do is = number_of_species+1, number_of_species + iadia 
    dgrid(nsp+iadia) = dgrid_G(is)
    tgrid(nsp+iadia) = tgrid_G(is) 
    do ix = 1, nx
      de(ix,nsp+iadia)  = de_G(ix,is) 
      tmp(ix,nsp+iadia) = tmp_G(ix,is) 
      tp(ix,nsp+iadia)  = tp_G(ix,is) 
      fp(ix,nsp+iadia)  = fp_G(ix,is)
      vp(ix,nsp+iadia)  = vp_G(ix,is)
    end do
  end do

  ! If reading profiles from file will need something akin to
  ! parallelize_geom for the the global x grids 

  ! gather the arrays for global in the radial direction 
  
  call gather_array(de_Gx(1:n_x_grid,:), n_x_grid,  nsp+iadia,           &
                  & de   (1:nx,      :), nx,        nsp+iadia,           &
                  & COMM_X_NE, COMM_DUMMY, ALLGATHER = .true.)
                  
  call gather_array(tmp_Gx(1:n_x_grid,:), n_x_grid, nsp+iadia,           &
                  & tmp   (1:nx,      :), nx,       nsp+iadia,           &
                  & COMM_X_NE, COMM_DUMMY, ALLGATHER = .true.)
  
  return 
  100 call gkw_abort('Error in reading the input.prof file in components') 
  
end subroutine components_profiles

!------------------------------------------------------------------------------
! This subroutine calculates the radial profiles according to analytic 
! functions. The functional form is chosen through the switch prof_type
! and influenced by the coefficients stored in pf(5) 
! 
! the options are 
! prof_type = 'tanh' 
!   profile form from tanh 
!   pf(1)  = density at reference point 
!   pf(2)  = density gradient at reference point 
!   pf(3)  = psi value of the reference point 
!   pf(4)  = width of the region of finite gradient 
!   pf(5)  = width of the region in which the gradient increases 
!------------------------------------------------------------------------------
subroutine call_prof(prof,dprof,prof_type,pf,val_norm) 

  use grid,    only : nx, gx, n_x_grid  
  use geom,    only : xgr, q_profile 
  use general, only : gkw_abort 

  real, intent(out) :: prof(nx), dprof(nx), val_norm
  character (len=lenswitch), intent(in) :: prof_type
  real, intent(in) :: pf(5)

  real, allocatable :: phelp(:) 
  real    :: ddd, rld, xgr_mp, width, delta, psie, qedge, qnul, a, b
  real    :: sval, dum1, dum2, widtG, xxx  
  integer :: ix, ierr 

  select case(prof_type) 
  case('file') 
    ! nothing is done here. The data is read from file 
    return 
    
  case('cosh2') 

    ddd    = pf(1)   ! the density / temperature at x = x0 
    rld    = pf(2)   ! R/L_G at x = x_0  
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R  
    width  = pf(4)   ! the normalized width  Delta x / R_0   
    widtG  = pf(5)   ! Width that determines the rise of the rld prof Delta_G / R_0 

    do ix = 1, nx 
      xxx = xgr(gx(ix)) 
      if (xxx.lt.xgr_mp - width/2) xxx = xgr_mp-width/2. 
      if (xxx.gt.xgr_mp + width/2) xxx = xgr_mp+width/2. 

      dprof(ix) = rld*(1.0E0 - 1.0E0 / cosh((xxx-(xgr_mp-width/2.E0))/widtG)**2 &
                &            - 1.0E0 / cosh((xxx-(xgr_mp+width/2.E0))/widtG)**2)
      prof(ix) = ddd*exp(-rld*(xxx-xgr_mp) +    &
               &              rld*widtG*tanh((xxx-(xgr_mp-width/2.E0))/widtG) + &
               &              rld*widtG*tanh((xxx-(xgr_mp+width/2.E0))/widtG))
    end do 
   
    val_norm = ddd 
  
  case('const')
    ddd = pf(1)
  
    do ix = 1, nx 
      prof(ix)   = ddd
      dprof(ix)  = 0.E0 
    end do 

    val_norm = ddd  

  case('tanh') 

    ddd    = pf(1)
    rld    = pf(2) 
    xgr_mp = pf(3) 
    width  = pf(4)  
    delta =  pf(5)  

    ! Profile and derivative 
    do ix = 1, nx 
      prof(ix)   = ddd*(1.E0 - rld*int_tanh_prof(xgr(gx(ix)),xgr_mp,width,delta)) 
      dprof(ix)  = rld*tanh_prof(xgr(gx(ix)),xgr_mp,width,delta) 
    end do 

    ! Value of the profile at xgr_mp 
    val_norm = ddd*(1.E0 - rld*int_tanh_prof(xgr_mp,xgr_mp,width,delta))

  case('exp_tanh') 

    ddd    = pf(1)   ! the density / temperature at x = x0 
    rld    = pf(2)   ! R/L_G at x = x_0  
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R  
    width  = pf(4)   ! the normalized width of the function = Delta x / R_0   
    psie   = pf(5)   ! unused 

    do ix = 1, nx 
      prof(ix) = ddd*exp(-rld*width*tanh((xgr(gx(ix))-xgr_mp)/width))
      dprof(ix) = rld/(cosh((xgr(gx(ix))-xgr_mp)/width)**2)
    end do

    val_norm = ddd 

  case('1_m_exp_tanh')

    ddd    = pf(1)   ! the density / temperature at x = x0
    rld    = pf(2)   ! R/L_G at x = x_0
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R
    width  = pf(4)   ! the normalized width of the function = Delta x / R_0
    psie   = pf(5)   ! unused

    do ix = 1, nx
      prof(ix) = 1.E0 - ddd*exp(-rld*width*tanh((xgr(gx(ix))-xgr_mp)/width))
      dprof(ix) = -(rld/(cosh((xgr(gx(ix))-xgr_mp)/width)**2))*(1-prof(ix))/prof(ix)
    end do

    val_norm = 1.E0 - ddd

  case('nth_p_exp_tanh')

    ddd    = pf(1)   ! the density / temperature at x = x0
    rld    = pf(2)   ! R/L_G at x = x_0
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R
    width  = pf(4)   ! the normalized width of the function = Delta x / R_0
    psie   = pf(5)   ! fraction of the thermal species

    do ix = 1, nx
      prof(ix) = psie + ddd*exp(-rld*width*tanh((xgr(gx(ix))-xgr_mp)/width))
      dprof(ix) = (rld/(cosh((xgr(gx(ix))-xgr_mp)/width)**2))*(prof(ix)-psie)/prof(ix)
    end do

    val_norm = psie + ddd

  case ( 'lin' )
  
    ddd    = pf(1)   ! the density / temperature at x = x0 
    rld    = pf(2)   ! R/L_G at x = x_0  
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R 
    width  = pf(4)   ! width of the function (Lin = 0.28*r_0)
    psie   = pf(5)   ! Value of psih (0.36)

    ierr = 0 
    allocate(phelp(n_x_grid),stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate  phelp in call_prof')          
             
    ! the gradient can be readily calculated also when parallelizing of the 
    ! radial domain. 
    do ix = 1, nx   
      xxx = xgr(gx(ix))
      dprof(ix) = rld*exp(-((xxx-xgr_mp)/(width*psie))**6)
    end do 

    ! always integrate the whole domain such that the normaization is known 
    ! on all processors (a more elegant solution is possible but requires a 
    ! bit of rewritting of the code) 
    phelp(1) = 0. 
    do ix = 2, n_x_grid  
      phelp(ix) = phelp(ix-1) - (xgr(ix)-xgr(ix-1))*rld*                  &
           & exp(-(((xgr(ix-1)+xgr(ix))/2.-xgr_mp)/(width*psie))**6)
      if ( (xgr(ix)-0.500001*(xgr(ix)-xgr(ix-1))) < xgr_mp                &
           &  .AND. xgr_mp <= (xgr(ix)+0.500001*(xgr(ix+1)-xgr(ix))) ) then 
           val_norm = exp(phelp(ix)) 
      end if
    end do
    
    ! calculate the profile and normalize 
    val_norm = ddd / val_norm 
    do ix = 1, nx 
      prof(ix) = val_norm * exp(phelp(gx(ix))) 
    end do 
       
    val_norm = ddd       
    
    ! deallocate the help array 
    deallocate(phelp) 
    
  case('orb') 

    ddd    = pf(1)   ! the temperature at s = s0 
    rld    = pf(2)   ! the gradient at s = s0 
    xgr_mp = pf(3)   ! the value of s0 
    width  = pf(4)   ! the width of the function 
    psie   = pf(5)   ! a / R0  

    ! the q - profile values (necessary to calculate s) 
    call q_profile(psie,qedge,dum1)
    ! transform q to qbar (which is the q actually used by ORB) 
    qedge = qedge * sqrt(1.0E0 - psie**2)
    call q_profile(0.E0,qnul,dum1)
    
    ! Rescale the gradient length to have exactly the input value at s0 
    dum1 = cosh(0.) 
    dum2 = cosh(xgr_mp/width) 
    dum1 = (1/dum1**2 - 1 / dum2**2) 
    rld  = rld / dum1 
 
    ! calculate the profiles 
    do ix = 1, nx  

      ! calculate s 
      sval = sqrt(log(1.E0 + (qedge - qnul)*xgr(gx(ix))**2/(qnul*psie**2)) / &
           &      log(1.E0 + (qedge - qnul)/ qnul ) ) 

      dum1 = cosh((sval-xgr_mp)/width) 
      dum2 = cosh(xgr_mp/width) 

      dprof(ix) = + rld*(1/dum1**2 - 1 / dum2**2) 
      prof(ix) =  (dum1/dum2)**(2.*rld*width**2) * exp ( rld * sval**2 / dum2**2 - &
               &  2.*rld*width*sval*tanh((sval-xgr_mp)/width)) 

      ! This write can be used to compare the profiles 
      !write(*,*) xgr(ix),sval

    end do 

    ! calculate the profile at peak width
    dum2 = cosh(xgr_mp/width)
    dum1 = (1.0/dum2)**(2.*rld*width**2) * exp ( rld * xgr_mp**2 / dum2**2 )
    ! normalize
    do ix = 1, nx
      prof(ix) = ddd*prof(ix)/dum1
    end do
 
    ! finally transform the gradient towards the normalized poloidal flux
    ! used in ORB into a gradient towards the radial coordinate
    do ix = 1, nx
      dum1 = 2.*xgr(gx(ix))*(qedge-qnul)/(log(qedge/qnul)*(qnul*psie**2 + &
           & (qedge - qnul)*xgr(gx(ix))**2))
      dprof(ix) = dprof(ix)*dum1
    end do 

    val_norm = ddd

  case('orb3') 

    ddd    = pf(1)   ! the density / temperature at x = x0 
    rld    = pf(2)   ! R/L_G at x = x_0  
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R  
    width  = pf(4)   ! the normalized width  Delta x / R_0   
    widtG  = pf(5)   ! Width that determines the rise of the rld prof Delta_G / R_0

    do ix = 1, nx
      xxx = xgr(gx(ix))
      ! if (xxx.lt.xgr_mp - width/2) xxx = xgr_mp-width/2. 
      ! if (xxx.gt.xgr_mp + width/2) xxx = xgr_mp+width/2. 

      dprof(ix) = 0.5E0*rld*(tanh((xxx-xgr_mp+width)/widtG) -                  &
                &            tanh((xxx-xgr_mp-width)/widtG))
      prof(ix)  = ddd * exp( - 0.5E0 * rld * widtG * log(                      &
                & cosh((xxx-xgr_mp+width)/widtG) /                             &
                & cosh((xxx-xgr_mp-width)/widtG) ) ) 

    end do 
   
    ! value at the centre 
    val_norm = ddd 
  case('exp_poly3')
    ! Results in a quadratic profile of the gradient length.
    ! Limiting cases of linear (for a = 0) and constant (for a=0 and b=0) are
    ! also possible.
    ! Some care has to be taken to get physical profiles.
    ! - To get a gradient length that is positive over the whole range, the
    !   inequality b^2 < 4a has to be fullfilled.
    !  prof(x) = ddd \cdot exp \left[ - rld \left(  \frac{a}{3} (x - xgr_{mp})^3
    !                                             + \frac{b}{2} (x - xgr_{mp})^2
    !                                             +             (x - xgr_{mp}) \right) \right]
    ! dprof(x) = rld \left[ a (x - xgr_{mp})^2 + b (x - xgr_{mp}) + 1 \right]

    ddd    = pf(1)   ! Profile value at x = xgr_{mp}.
    rld    = pf(2)   ! Gradient length value at x = xgr_{mp}.
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R.
    a      = pf(4)   ! Constant for quadratic term.
    b      = pf(5)   ! Constant for linear term.

    do ix = 1,nx
      xxx = xgr(gx(ix))
      prof(ix)  = ddd*exp(-rld*(  a/3.0*(xxx-xgr_mp)**3    &
                &               + b/2.0*(xxx-xgr_mp)**2    &
                &               +       (xxx-xgr_mp)))
      dprof(ix) = rld*(a*(xxx-xgr_mp)**2 + b*(xxx-xgr_mp) + 1)
    end do

    val_norm = ddd

  case('exp_poly6')
    ! Results in a fifth and third order term as well as a constant for the
    ! gradient length.
    !  prof(x) = ddd \cdot exp \left[ - rld \left(  \frac{a}{6} (\frac{x - xgr_{mp}}{width})^6
    !                                             + \frac{1}{4} (\frac{x - xgr_{mp}}{width})^4
    !                                             +             (x - xgr_{mp}) \right) \right]
    ! dprof(x) = rld \left[  \frac{a}{width} (\frac{x - xgr_{mp}}{width})^5
    !                      + \frac{1}{width} (\frac{x - xgr_{mp}}{width})^3
    !                      + 1 \right]
    ddd    = pf(1)   ! Profile value at x = xgr_{mp}.
    rld    = pf(2)   ! Gradient length value at x = xgr_{mp}.
    xgr_mp = pf(3)   ! the value of x0 = r_0 / R.
    a      = pf(4)   ! the normalized width  Delta x / R_0
    b      = pf(5)   ! Constant for the term to the fifth power.

    do ix = 1,nx
      xxx = xgr(gx(ix))
      prof(ix)  = ddd*exp(-rld*(  a/6.0*(xxx-xgr_mp)**6   &
                &               + b/4.0*(xxx-xgr_mp)**4   &
                &               +       (xxx-xgr_mp)))
      dprof(ix) = rld*(a*(xxx-xgr_mp)**5 + b*(xxx-xgr_mp)**3 + 1.0)
    end do

    val_norm = ddd

  case default
 
    call gkw_abort('Unknown profile option')

  end select  

end subroutine call_prof

!------------------------------------------------------------------------------
! This subroutine calculates the radial current profile according to  
! the q-profile calculated in geom.  This is applied only to the electrons.
!------------------------------------------------------------------------------
subroutine calc_ecur_prof

  use grid,    only : nx, gx, n_x_grid
  use geom,    only : dxgr, xgr, qx
  use geom,    only : signB, signJ
  
  integer :: ix, p, n, is
  real :: dum, gppR2np, gppR2n, gppR2nm, invRp, invRm 
  real :: psiip, psiim, bpolp, bpol, bpolm

  gppR2np = 0.E0
  gppR2n  = 0.E0
  gppR2nm = 0.E0
  invRp = 0.E0
  invRm = 0.E0
  is = 1
  
  do is = 1,nsps
    if(signz(is).lt.0)then
  
      do ix= 1,nx 
      
        if(gx(ix).eq.1)then
          n = 1
        else
          n = gx(ix - 1)
        endif   
     
        if(gx(ix).eq.n_x_grid)then
          p = n_x_grid
        else
          p = gx(ix + 1)
        endif
    
    !do i=1,ns
    !  !Flux surface averages of the the quantities g_psi_psi/R^{2} and
    !  !1/R (Which are interpolated to half grid points linearly)
    !  gppR2n = gppR2n + metric_G(ix,i,1,1)*ints(i)/(Rfun(ix,i)*Rfun(ix,i))
    !  gppR2np = gppR2np + metric_G(p,i,1,1)*ints(i)/(Rfun(p,i)*Rfun(p,i))
    !  gppR2nm = gppR2nm + metric_G(n,i,1,1)*ints(i)/(Rfun(n,i)*Rfun(n,i))
    !  invRp = invRp + 0.5E0*ints(i)/(Rfun(p,i)+Rfun(ix,i))
    !  invRm = invRm + 0.5E0*ints(i)/(Rfun(ix,i)+Rfun(n,i))
    !enddo

    !bs = ffun*bn
     
        !At the moment this doesnt work in general geometry.  
        !Only for circular geometry. 
        !Radial coordinate linearly interpolated to mid point plus and minus
        psiip = 0.5E0*(xgr(p)+xgr(gx(ix)))
        psiim = 0.5E0*(xgr(gx(ix))+xgr(n))

        bpolp = xgr(p)*xgr(p)/qx(p)
        bpol = xgr(gx(ix))*xgr(gx(ix))/qx(gx(ix))
        bpolm = xgr(n)*xgr(n)/qx(n)

        dum = ((bpolp-bpol)/psiip) - ((bpol-bpolm)/psiim)

        !ecurx is the radial gradient of the flow -R_0/vthref)du/dr
        dum = dum/(rhostar*rhostar*dxgr*dxgr)
        if(gx(ix).eq.n_x_grid)then
          ecurx(ix)=0.E0
        else if(gx(ix).eq.1)then
          ecurx(ix)=0.E0
        else
          ecurx(ix) = -signB*signJ*dum*rhostar/(signz_G(is)*de_G(ix,is)*beta_ref)
        endif
      enddo
    endif
  enddo
  
  !The inner boundary condition-> Uprime is linearly
  !interpolated from psi=0 (axis) where uprime is zero.
  !This is put here because it requires the 2nd point
  if(gx(1).eq.1)then
    ecurx(1)=(ecurx(gx(2))/xgr(gx(2)))*xgr((gx(1)))
  endif
  if(gx(nx).eq.n_x_grid)then
    n=nx
    ecurx(n) = (2.E0*ecurx(n-1)-ecurx(n-2))
  endif
end subroutine calc_ecur_prof

!------------------------------------------------------------------------------
! Generic function that introduced the tanh profiles
! tanh_prof(psi) = 0.5*(tanh(psi+/delta) - tanh(psi-/delta)) 
! where psi+ = psi - psi0 + width/2 
!       psi- = psi - psi0 - width/2
! The profile is centered around psi0, and has a width (width). The quantity 
! delta determines how fast the profile rises at the edges. 
!------------------------------------------------------------------------------
function tanh_prof(psi,psi0,width,delta) 

  real :: tanh_prof
  real, intent(in) :: psi, psi0, width, delta 

  tanh_prof = 0.5E0*(tanh((psi-psi0+width/2.E0)/delta) -  & 
            &        tanh((psi-psi0-width/2.E0)/delta) ) 

end function tanh_prof 

!------------------------------------------------------------------------------
! Genereric function that determines the integral of the tanh_prof functions
! for parameters see function tanh_prof 
!------------------------------------------------------------------------------
function int_tanh_prof(psi,psi0,width,delta) 

  real :: int_tanh_prof
  real, intent(in) :: psi, psi0, width, delta 

  int_tanh_prof = 0.5*delta*(log(cosh((psi-psi0+width/2.0E0)/delta)) - & 
                &            log(cosh((psi-psi0-width/2.0E0)/delta))) 
  
end function int_tanh_prof



!------------------------------------------------------------------------------
!> This subroutine deallocates the arrays of components 
!------------------------------------------------------------------------------

subroutine components_deallocate 

  if (allocated(mas))      deallocate(mas)
  if (allocated(signz))    deallocate(signz)
  if (allocated(tmp))      deallocate(tmp)
  if (allocated(de))       deallocate(de)
  if (allocated(tp))       deallocate(tp)
  if (allocated(fp))       deallocate(fp)
  if (allocated(vp))       deallocate(vp)
  if (allocated(vthrat))   deallocate(vthrat)
  if (allocated(mas_G))    deallocate(mas_G)
  if (allocated(signz_G))  deallocate(signz_G)
  if (allocated(tmp_G))    deallocate(tmp_G)
  if (allocated(de_G))     deallocate(de_G)
  if (allocated(vthrat_G)) deallocate(vthrat_G)
  if (allocated(vp_G))     deallocate(vp_G)
  if (allocated(tp_G))     deallocate(tp_G)
  if (allocated(fp_G))     deallocate(fp_G)

end subroutine components_deallocate 


end module components
