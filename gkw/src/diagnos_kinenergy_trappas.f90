!------------------------------------------------------------------------------
!> This diagnostic computes the various contributions to the time
!> derivative of the kinetic energy (i.e the work of the particles in
!> the electric field), similar to ORB5.
!>
!> It outputs the potential and magnetic energy, and the time variation of 
!> kinetic energy for all the velocity contributions: gradB drifts,
!> parallel velocity...
!> The growth rate can be computed from these outputs. Thus,
!> this diagnostic can reveal details about the nature of an
!> instability.
!>
!> For a linear and electrostatic instability the fluctuations grow
!> exponentially.
!>   phi(k,t) = phi_0(k) exp[- i omega t + gamma t]
!>   f(k,t) = f_0(k) exp[- i omega t + gamma t]
!> where phi_0 and f_0 are complex amplitudes, so that phi and f may have a
!> phase difference.
!> In what follows, an asterisk means complex conjugation.
!>
!> The potential energy, being quadratic in fluctuations,
!> growths accordingly:
!>   Epot = \sum_k \int dv Z e f_0 exp[- i omega t + gamma t] phi_0* exp[+ i omega t + gamma t]
!>        = \sum_k \int dv Z e f_0 phi_0* exp[2 gamma t]
!> A typical linear run with nmod=1 does not have a zeromode, thus
!>        = 2 real (\int dv Z e f_0 phi_0* exp[2 gamma t] )
!>
!> Note that Epot is real (see also the section 'Hermitian symmetry of binormal
!> spectra' in the GKW manual).
!> Its time derivative is then
!> 
!>   dEpot/dt = 2 gamma Epot
!>
!> We assume then that total energy is conserved (again, Ekin is a real valued
!> quantity)
!>   dEpot/dt = -dEkin/dt
!>
!> Total energy conservation yields that the growth rate is given by:
!>   gamma = (dEpot/dt) / (2 Epot)  = (-dEkin/dt) / (2 Epot)
!>
!> The change in kinetic energy can be computed: it is the dotproduct
!> of current and E-field.
!>
!> dEkin/dt = \sum_{species} \sum_k \int dv
!>   Z_s e f (v_||\mathbf{b} + \mathbf{v}_D ) \cdot \mathbf{E}
!>
!> The E-field is determined by the electrostatic potential and the parallel
!> component of the vector potential
!>   \mathbf{E} = -\ga{\nabla\phi} - \mathbf{b} d\ga{Apar}/dt
!>              = -nabla\ga{\phi}  - \mathbf{b} d\ga{Apar}/dt
!>
!> and the field energy is actually the electrostatic potential energy plus the magnetic field
!> energy.
!>   Epot + Emag = \int dX dv (Ze f \ga{\phi} + 1/(2 mu0) |grad Apar|^2)
!>
!> Normalising the expression for Emag makes the plasma beta appear,
!> as a factor 1/beta_ref.
!>
!> We assume that the mode has already converged and
!> evolves according to
!>       \ga{Apar}(k) = \ga{Apar}_0(k) exp[-i omega t + gamma t]
!>
!> For a given wavenumber k, we have then
!>
!> dEkin(k)/dt = 2 real( \sum_{species} \int dv
!>   Z_s e f(k)   (v_||\mathbf{b} + \mathbf{v}_D + \mathbf{v}_\Chi) \cdot (-nabla\ga{\phi})* +
!> + Z_s e f(k)   (v_||\mathbf{b}) \cdot (-\mathbf{b} d\ga{Apar}/dt)* )
!>
!> dEkin(k)/dt = 2 real( \sum_{species} \int dv
!>   Z_s e f(k)   (v_||\mathbf{b} + \mathbf{v}_D + \mathbf{v}_\Chi) \cdot (-nabla\ga{\phi})* +
!> + Z_s e f(k)   (v_||\mathbf{b}) \cdot (-\mathbf{b} d\ga{Apar}/dt)* )
!>
!> A couple of abbreviations come in handy here:
!>
!> dEkin(k)/dt = ES(k) + I_m(k)                        [eq.1]
!>
!> where ES are all the electrostatic terms
!> ES(k) = 2 real( \sum_{species} \int dv
!>   Z_s e f(k)   (v_||\mathbf{b} + \mathbf{v}_D + \mathbf{v}_\Chi) \cdot (-nabla\ga{\phi})* )
!>
!> and I_m denotes the integral with the magnetic vector potential.
!>
!>   I_m = 2 real( \sum_{species} \int dv Z_s e f(k) v_|| (-(-i omega + gamma)\ga{Apar(k)})* )
!>       = 2 real( \sum_{species} \int dv Z_s e      v_||  -(+i omega + gamma)\ga{Apar(k)}*f(k) )
!>
!> In this form, to be able to solve the above equation for gamma, the real frequency omega must
!> be known. This is, because due to the imaginary part of (\ga{Apar(k)}* f(k)) it
!> contributes to the real valued I_m.

!> Although the overall phase of the converged linear simulation is arbitrary
!> and it is legitimate to rotate f and all fields in the complex plane by the same amount, such a rotation
!> cannot change the phase of a product like (\ga{Apar(k)}* f(k))
!> because exp(-i alpha)exp(+i alpha) = 1

!> As it was already assumed that the mode has very well converged, we
!> may use growthrate gamma and real frequency omega as computed by
!> the diagnos_growth_freq diagnostic.
!>
!> Finally then, we can compute
!> the electromagnetic contribution I_m to dEkin/dt.
!>
!>   dEkin/dt = ES + I_m
!>
!> and
!>
!>   gamma = -dEkin/dt / (2 (Epot + Emag))
!>         = gamma_contrib_vpar + gamma_contrib_gradb + ... + gamma_contrib_vpar_apar
!>
!>
!> Known issues:
!> * very large number of NPERIOD required to get equivalence between
!> the new and old growth rate calculation -> likely due to the free
!> streaming boundary condition at the end of the field line (this
!> should be included in the calculation as for the entropy diagnostic)
!> * good parallel resolution required
!> * low dissipation is also beneficial to get equivalence between the new
!> and old growth rate calculation -> same as above, need to be
!> included in the calculation
!> * consistent magnetic geometry required
!>
!> To do:
!> * include boundary loss and dissipative loss
!>
!> For more informations:
!> A. Bottino, PhD thesis, EPFL 2938 (2004)
!------------------------------------------------------------------------------
module diagnos_kinenergy_trappas

  implicit none

  private

  logical, save, public :: lcalc_kinenergy_trappas

  logical, save, public :: lcalc_kinenergy_trappas_sx
  logical, save, public :: lcalc_kinenergy_trappas_xy

  logical, parameter :: mag_en_with_fdis = .true.

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: output

  !> the variation of kinetic energy per mode : (nmod,nx,number_of_species) for the 
  !> different contributions
  complex, save, allocatable, dimension(:,:,:) :: kin_vpar_ps,kin_gradb_ps
  complex, save, allocatable, dimension(:,:,:) :: kin_betaprime_ps,kin_coriolis_ps
  complex, save, allocatable, dimension(:,:,:) :: kin_centrifugal_ps,kin_centpot_ps
  complex, save, allocatable, dimension(:,:,:) :: kin_vpar_apar_ps

  complex, save, allocatable, dimension(:,:,:) :: kin_vpar_tr,kin_gradb_tr
  complex, save, allocatable, dimension(:,:,:) :: kin_betaprime_tr,kin_coriolis_tr
  complex, save, allocatable, dimension(:,:,:) :: kin_centrifugal_tr,kin_centpot_tr
  complex, save, allocatable, dimension(:,:,:) :: kin_vpar_apar_tr

  complex, save, allocatable, dimension(:,:,:) :: kin_vpar_ps_sx,kin_gradb_ps_sx
  complex, save, allocatable, dimension(:,:,:) :: kin_vpar_tr_sx,kin_gradb_tr_sx

  !> the potential energy (E-field energy) per mode
  complex, save, allocatable, dimension(:,:,:) :: pot_en
  complex, save, allocatable, dimension(:,:,:) :: pot_en_sx
  !> the magnetic field energy per mode
  complex, save, allocatable, dimension(:,:) :: mag_en, mag_en2

  complex, save, allocatable, dimension(:) :: kin_tot_vpar_ps
  complex, save, allocatable, dimension(:) :: kin_tot_gradb_ps,kin_tot_betaprime_ps
  complex, save, allocatable, dimension(:) :: kin_tot_coriolis_ps,kin_tot_centrifugal_ps
  complex, save, allocatable, dimension(:) :: kin_tot_centpot_ps
  complex, save, allocatable, dimension(:) :: kin_tot_vpar_apar_ps

  complex, save, allocatable, dimension(:) :: kin_tot_vpar_tr
  complex, save, allocatable, dimension(:) :: kin_tot_gradb_tr,kin_tot_betaprime_tr
  complex, save, allocatable, dimension(:) :: kin_tot_coriolis_tr,kin_tot_centrifugal_tr
  complex, save, allocatable, dimension(:) :: kin_tot_centpot_tr
  complex, save, allocatable, dimension(:) :: kin_tot_vpar_apar_tr

  complex, save, allocatable, dimension(:) :: pot_tot_en
  complex, save :: mag_tot_en, mag_tot_en2

  integer, save :: i_pot_en
  integer, save :: i_vpar_ps, i_vpar_tr, i_gradb_ps, i_gradb_tr
  integer, save :: i_pot_en_sx
  integer, save :: i_pot_en_xy
  integer, save :: i_vpar_ps_sx, i_vpar_tr_sx, i_gradb_ps_sx, i_gradb_tr_sx
  integer, save :: i_vpar_ps_xy, i_vpar_tr_xy, i_gradb_ps_xy, i_gradb_tr_xy

  integer, save :: i_coriolis_ps, i_coriolis_tr, i_centrifugal_ps
  integer, save :: i_centrifugal_tr, i_centpot_ps, i_centpot_tr
  
  integer, save :: i_betaprime_ps, i_betaprime_tr
  integer, save :: i_vpar_apar_ps, i_vpar_apar_tr
  integer, save :: i_mag_en
  
contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()

   lcalc_kinenergy_trappas = .false.
   lcalc_kinenergy_trappas_sx = .false.
   lcalc_kinenergy_trappas_xy = .false.

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcalc_kinenergy_trappas,1)
    call mpibcast(lcalc_kinenergy_trappas_sx,1)
    call mpibcast(lcalc_kinenergy_trappas_xy,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : spectral_radius, non_linear, flux_tube, nlapar
    use general, only : gkw_warn
    use geom, only : geom_type
    use diagnos_generic, only : xy_estep
    use components, only : veta, betaprime_type, veta_prime
    use global, only : r_tiny
    use grid, only : nmod

    if(.not. lcalc_kinenergy_trappas) return

    lcalc_kinenergy_trappas_sx = lcalc_kinenergy_trappas
    lcalc_kinenergy_trappas_xy = lcalc_kinenergy_trappas .and. xy_estep

    if(.not.spectral_radius) then
      call gkw_warn("lcalc_kinenergy_trappas is only implemented for spectral &
         & fluxtube runs. Disabled it.")
      ! because of the derivatives
      lcalc_kinenergy_trappas = .false.
    end if

    if(.not. flux_tube) then
      call gkw_warn("lcalc_kinenergy_trappas is only implemented for spectral &
         & fluxtube runs. Disabled it.")
      ! because of normalisation factors
      lcalc_kinenergy_trappas = .false.
    end if

    if(geom_type == 's-alpha') then
      call gkw_warn("kinenergy_trappas diagnostic does hardly work with s-alpha &
         & geometry, energy conservation is too bad then. Disabled it.")
      lcalc_kinenergy_trappas = .false.
    end if

    if(nlapar .and. betaprime_type == 'ref') then
      if(all(abs(veta_prime) < r_tiny) .and. all(abs(veta) < r_tiny)) then
        call gkw_warn("kinenergy_trappas diagnostic: betaprime_ref = 0 is surely&
           & not consistent with beta_ref /= 0. This may lead to (small) &
           & inconsistency, compared to the ordinary growth rate diagnostic.")
        ! lcalc_kinenergy_trappas = .false.
      else
        call gkw_warn("kinenergy_trappas diagnostic might not be consistent &
           & if case betaprime is not consistently chosen (you might want to use&
           & betaprime_type='sp')")
      end if
    end if

    if(non_linear) then
      call gkw_warn("lcalc_kinenergy_trappas does not consider the nonlinear ExB&
         & drift.")
    end if

    if(nmod > 1) then
      ! FIXME remove this when fixed
      call gkw_warn("diagnos_kinenergy_trappas is not really made for many &
         & linear modes in one run at the moment.")
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use grid, only : number_of_species, n_x_grid, nmod, n_s_grid
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use global, only : PHI_GA_FIELD, APAR_FIELD, APAR_GA_FIELD
    use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS, X_GHOSTCELLS
    use global, only : int2char_zeros
    use control, only : nlapar
    use mpiinterface, only : root_processor
    use global, only : dotdat
    use rotation, only : coriolis, cf_drift, cf_trap
    logical, intent(inout) :: requirements(:,:)

    if(.not. lcalc_kinenergy_trappas) return

    requirements(PHI_GA_FIELD, LOCAL_DATA) = .true.
    requirements(PHI_GA_FIELD, S_GHOSTCELLS) = .true.
    requirements(APAR_GA_FIELD, LOCAL_DATA) = .true.
    requirements(APAR_FIELD, LOCAL_DATA) = .true.
    requirements(APAR_FIELD, X_GHOSTCELLS) = .true.

    if(root_processor) then

      call open_real_lu('gamma_contrib_pot_en', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ number_of_species /), &
         & ascii_fmt, i_pot_en)
      call attach_metadata_grid(i_pot_en, 'time', 'number_of_species', ascii_fmt)
      call attach_metadata(i_pot_en, phys_unit_key, not_avail, ascii_fmt)
      if(nlapar) then
      call attach_metadata(i_pot_en, description_key, &
         & 'The field energy 2*E_pot (also potential energy) plus the magnetic&
         & field energy.', ascii_fmt)
      else
        call attach_metadata(i_pot_en, description_key, &
         & 'The field energy 2*E_pot (also potential energy).', ascii_fmt)
      end if
      call attach_metadata(i_pot_en, comments_key, 'This includes a factor 2&
         &, i.e. This is the same as ene_e2 (plus ene_m, if electromagnetic) &
         & from diagnos_energetics, apart from the factor 2.&
         & This can be used to compute the&
         & growthrate from the contributions to time deriv. of the kinetic Energy:&
         & gamma = -1/(2*Epot) dEkin/dt.&
         & ', ascii_fmt)

      call open_real_lu('gamma_contrib_pot_en_sx', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ n_s_grid, n_x_grid, number_of_species /), &
         & ascii_fmt, i_pot_en_sx)

      call open_real_lu('gamma_contrib_pot_en_xy', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ nmod, n_x_grid, number_of_species /), &
         & ascii_fmt, i_pot_en_xy)
      

      call open_real_lu('gamma_contrib_vpar_ps', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ number_of_species /), &
         & ascii_fmt, i_vpar_ps)
      call attach_metadata_grid(i_vpar_ps, 'time', 'number_of_species', ascii_fmt)
      call attach_metadata(i_vpar_ps, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_vpar_ps, description_key, 'vpar contribution for &
         & passing particles', ascii_fmt)
      call attach_metadata(i_vpar_ps, comments_key, 'This is one contribution to&
         & the time derivative of the kinetic energy.',&
         & ascii_fmt)

      call open_real_lu('gamma_contrib_vpar_ps_sx', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ n_s_grid, n_x_grid, number_of_species /), &
         & ascii_fmt, i_vpar_ps_sx)

      call open_real_lu('gamma_contrib_vpar_ps_xy', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ nmod, n_x_grid, number_of_species /), &
         & ascii_fmt, i_vpar_ps_xy)
      

      call open_real_lu('gamma_contrib_vpar_tr', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ number_of_species /), &
         & ascii_fmt, i_vpar_tr)
      call attach_metadata_grid(i_vpar_tr, 'time', 'number_of_species', ascii_fmt)
      call attach_metadata(i_vpar_tr, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_vpar_tr, description_key, 'vpar contribution for&
         & trapped particles', ascii_fmt)
      call attach_metadata(i_vpar_tr, comments_key, 'This is one contribution to&
         & the time derivative of the kinetic energy.', ascii_fmt)

      call open_real_lu('gamma_contrib_vpar_tr_sx', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ n_s_grid, n_x_grid, number_of_species /), &
         & ascii_fmt, i_vpar_tr_sx)

      call open_real_lu('gamma_contrib_vpar_tr_xy', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ nmod, n_x_grid, number_of_species /), &
         & ascii_fmt, i_vpar_tr_xy)
      

      call open_real_lu('gamma_contrib_gradb_ps', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ number_of_species /), &
         & ascii_fmt, i_gradb_ps)
      call attach_metadata_grid(i_gradb_ps, 'time', 'number_of_species', ascii_fmt)
      call attach_metadata(i_gradb_ps, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_gradb_ps, description_key, 'gradB drift contribution for&
         & passing particles', ascii_fmt)
      call attach_metadata(i_gradb_ps, comments_key, 'This is one contribution to&
         & the time derivative of the kinetic energy.', ascii_fmt)

      call open_real_lu('gamma_contrib_gradb_ps_sx', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ n_s_grid, n_x_grid, number_of_species /), &
         & ascii_fmt, i_gradb_ps_sx)

      call open_real_lu('gamma_contrib_gradb_ps_xy', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ nmod, n_x_grid, number_of_species /), &
         & ascii_fmt, i_gradb_ps_xy)
      

      call open_real_lu('gamma_contrib_gradb_tr', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ number_of_species /), &
         & ascii_fmt, i_gradb_tr)
      call attach_metadata_grid(i_gradb_tr, 'time', 'number_of_species', ascii_fmt)
      call attach_metadata(i_gradb_tr, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_gradb_tr, description_key, 'gradB drift contribution for&
         & trapped particles', ascii_fmt)
      call attach_metadata(i_gradb_tr, comments_key, 'This is one contribution to&
         & the time derivative of the kinetic energy.', ascii_fmt)

      call open_real_lu('gamma_contrib_gradb_tr_sx', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ n_s_grid, n_x_grid, number_of_species /), &
         & ascii_fmt, i_gradb_tr_sx)

      call open_real_lu('gamma_contrib_gradb_tr_xy', &
         & 'diagnostic/diagnos_kinenergy_trappas', &
         & (/ nmod, n_x_grid, number_of_species /), &
         & ascii_fmt, i_gradb_tr_xy)

      
      if(coriolis) then
        call open_real_lu('gamma_contrib_coriolis_ps', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_coriolis_ps)
        call attach_metadata_grid(i_coriolis_ps, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_coriolis_ps, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_coriolis_ps, description_key, 'coriolis &
           & drift contribution for passing particles', ascii_fmt)
        call attach_metadata(i_coriolis_ps, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)

        call open_real_lu('gamma_contrib_coriolis_tr', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_coriolis_tr)
        call attach_metadata_grid(i_coriolis_tr, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_coriolis_tr, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_coriolis_tr, description_key, 'coriolis &
           & drift contribution for trapped particles', ascii_fmt)
        call attach_metadata(i_coriolis_tr, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)
      end if

      if(cf_drift) then
        call open_real_lu('gamma_contrib_centrifugal_ps', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_centrifugal_ps)
        call attach_metadata_grid(i_centrifugal_ps, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_centrifugal_ps, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_centrifugal_ps, description_key, 'centrifugal &
           & drift contribution for passing particles', ascii_fmt)
        call attach_metadata(i_centrifugal_ps, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)

        call open_real_lu('gamma_contrib_centrifugal_tr', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_centrifugal_tr)
        call attach_metadata_grid(i_centrifugal_tr, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_centrifugal_tr, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_centrifugal_tr, description_key, 'centrifugal &
           & drift contribution for trapped particles', ascii_fmt)
        call attach_metadata(i_centrifugal_tr, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)
      end if

      ! if(cf_trap) then
        call open_real_lu('gamma_contrib_centpot_ps', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_centpot_ps)
        call attach_metadata_grid(i_centpot_ps, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_centpot_ps, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_centpot_ps, description_key, 'centrifugal &
           & potential drift contribution for passing particles', ascii_fmt)
        call attach_metadata(i_centpot_ps, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)

        call open_real_lu('gamma_contrib_centpot_tr', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_centpot_tr)
        call attach_metadata_grid(i_centpot_tr, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_centpot_tr, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_centpot_tr, description_key, 'centrifugal &
           & potential drift contribution for trapped particles', ascii_fmt)
        call attach_metadata(i_centpot_tr, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)
      ! end if

      if(nlapar) then
        call open_real_lu('gamma_contrib_betaprime_ps', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_betaprime_ps)
        call attach_metadata_grid(i_betaprime_ps, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_betaprime_ps, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_betaprime_ps, description_key, 'betaprime &
           & drift contribution for passing particles', ascii_fmt)
        call attach_metadata(i_betaprime_ps, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)

        call open_real_lu('gamma_contrib_betaprime_tr', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_betaprime_tr)
        call attach_metadata_grid(i_betaprime_tr, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_betaprime_tr, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_betaprime_tr, description_key, 'betaprime &
           & drift contribution for trapped particles', ascii_fmt)
        call attach_metadata(i_betaprime_tr, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)

        call open_real_lu('gamma_contrib_vpar_apar_ps', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_vpar_apar_ps)
        call attach_metadata_grid(i_vpar_apar_ps, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_vpar_apar_ps, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_vpar_apar_ps, description_key, 'vpar &
           &  magnetic contribution for passing particles', ascii_fmt)
        call attach_metadata(i_vpar_apar_ps, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)

        call open_real_lu('gamma_contrib_vpar_apar_tr', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ number_of_species /), &
           & ascii_fmt, i_vpar_apar_tr)
        call attach_metadata_grid(i_vpar_apar_tr, 'time', 'number_of_species', ascii_fmt)
        call attach_metadata(i_vpar_apar_tr, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_vpar_apar_tr, description_key, 'vpar &
           & magnetic contribution for trapped particles', ascii_fmt)
        call attach_metadata(i_vpar_apar_tr, comments_key, 'This is one contribution to&
           & the time derivative of the kinetic energy.', ascii_fmt)

        call open_real_lu('gamma_contrib_mag_en', &
           & 'diagnostic/diagnos_kinenergy_trappas', &
           & (/ 2 /), &
           & ascii_fmt, i_mag_en)
        call attach_metadata_grid(i_mag_en, 'time', ascii_fmt)
        call attach_metadata(i_mag_en, phys_unit_key, not_avail, ascii_fmt)
        if(nlapar) then
          call attach_metadata(i_mag_en, description_key, &
             & 'The magnetic&
             & field energy.', ascii_fmt)
        else
          call attach_metadata(i_mag_en, description_key, &
             & 'The magnetic field energy).', ascii_fmt)
        end if
        call attach_metadata(i_mag_en, comments_key, 'This includes a factor 2&
           & ', ascii_fmt)
      end if
    end if

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use control, only : nlapar
    use rotation, only : coriolis, cf_drift, cf_trap
    use grid, only : nmod,number_of_species,n_s_grid,n_x_grid
    use general, only : gkw_abort
    integer :: ierr

    if(.not. lcalc_kinenergy_trappas) return
    
    allocate(kin_vpar_ps(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_vpar_ps')
    allocate(kin_gradb_ps(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_gradb_ps')
    
    allocate(kin_vpar_tr(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_vpar_tr')
    allocate(kin_gradb_tr(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_gradb_tr')
    
    allocate(kin_tot_vpar_ps(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_vpar_ps')
    allocate(kin_tot_gradb_ps(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_gradb_ps')
    
    allocate(kin_tot_vpar_tr(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_vpar_tr')
    allocate(kin_tot_gradb_tr(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_gradb_tr')

    if(coriolis) then
      allocate(kin_coriolis_ps(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_coriolis_ps')
      allocate(kin_coriolis_tr(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_coriolis_tr')
      allocate(kin_tot_coriolis_ps(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_coriolis_ps')
      allocate(kin_tot_coriolis_tr(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_coriolis_tr')
    end if
    
    if(cf_drift) then
      allocate(kin_centrifugal_ps(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_centrifugal_ps')
      allocate(kin_centrifugal_tr(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_centrifugal_tr')
      allocate(kin_tot_centrifugal_ps(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_centrifugal_ps')
      allocate(kin_tot_centrifugal_tr(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_centrifugal_tr')
    end if

    ! if(cf_trap) then
      allocate(kin_centpot_ps(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_centpot_ps')
      allocate(kin_centpot_tr(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_centpot_tr')
      allocate(kin_tot_centpot_ps(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_centpot_ps')
      allocate(kin_tot_centpot_tr(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_centpot_tr')
    ! end if
    
    if(nlapar) then
      allocate(kin_betaprime_ps(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_betaprime_ps')
      allocate(kin_betaprime_tr(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_betaprime_tr')
      
      allocate(kin_tot_betaprime_ps(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_betaprime_ps')
      allocate(kin_tot_betaprime_tr(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_betaprime_tr')
      
      allocate(kin_vpar_apar_ps(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_vpar_apar_ps')
      allocate(kin_vpar_apar_tr(nmod,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_vpar_apar_tr')
      allocate(kin_tot_vpar_apar_ps(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_vpar_apar_ps')
      allocate(kin_tot_vpar_apar_tr(number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_tot_vpar_apar_tr')

      allocate(mag_en(nmod,n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: mag_en')
      allocate(mag_en2(nmod,n_x_grid),stat=ierr)
    end if

    allocate(pot_en(nmod,n_x_grid,number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pot_en')
    allocate(pot_tot_en(number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: pot_tot_en')
    
    if(lcalc_kinenergy_trappas_sx) then
      allocate(kin_vpar_ps_sx(n_s_grid,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_vpar_ps')
      allocate(kin_gradb_ps_sx(n_s_grid,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_gradb_ps')

      allocate(kin_vpar_tr_sx(n_s_grid,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_vpar_tr')
      allocate(kin_gradb_tr_sx(n_s_grid,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: kin_gradb_tr')

      allocate(pot_en_sx(n_s_grid,n_x_grid,number_of_species),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: pot_en')
    end if

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    use control, only : nlapar
    use rotation, only : coriolis, cf_drift, cf_trap

    if(.not. lcalc_kinenergy_trappas) return

    if(allocated(kin_vpar_ps)) deallocate(kin_vpar_ps)
    if(allocated(kin_gradb_ps)) deallocate(kin_gradb_ps)

    if(allocated(kin_vpar_tr)) deallocate(kin_vpar_tr)
    if(allocated(kin_gradb_tr)) deallocate(kin_gradb_tr)

    if(allocated(kin_tot_vpar_ps)) deallocate(kin_tot_vpar_ps)
    if(allocated(kin_tot_gradb_ps)) deallocate(kin_tot_gradb_ps)

    if(allocated(kin_tot_vpar_tr)) deallocate(kin_tot_vpar_tr)
    if(allocated(kin_tot_gradb_tr)) deallocate(kin_tot_gradb_tr)

    if(coriolis) then
      if(allocated(kin_coriolis_ps)) deallocate(kin_coriolis_ps)
      if(allocated(kin_coriolis_tr)) deallocate(kin_coriolis_tr)
      if(allocated(kin_tot_coriolis_ps)) deallocate(kin_tot_coriolis_ps)
      if(allocated(kin_tot_coriolis_tr)) deallocate(kin_tot_coriolis_tr)
    end if

    if(cf_drift) then
      if(allocated(kin_centrifugal_ps)) deallocate(kin_centrifugal_ps)
      if(allocated(kin_centrifugal_tr)) deallocate(kin_centrifugal_tr)
      if(allocated(kin_tot_centrifugal_ps)) deallocate(kin_tot_centrifugal_ps)
      if(allocated(kin_tot_centrifugal_tr)) deallocate(kin_tot_centrifugal_tr)
    end if

    ! if(cf_trap) then
      if(allocated(kin_centpot_ps)) deallocate(kin_centpot_ps)
      if(allocated(kin_centpot_tr)) deallocate(kin_centpot_tr)
      if(allocated(kin_tot_centpot_ps)) deallocate(kin_tot_centpot_ps)
      if(allocated(kin_tot_centpot_tr)) deallocate(kin_tot_centpot_tr)
    ! end if

    if(nlapar) then
      if(allocated(kin_betaprime_ps)) deallocate(kin_betaprime_ps)
      if(allocated(kin_betaprime_tr)) deallocate(kin_betaprime_tr)
      if(allocated(kin_tot_betaprime_ps)) deallocate(kin_tot_betaprime_ps)
      if(allocated(kin_tot_betaprime_tr)) deallocate(kin_tot_betaprime_tr)
      if(allocated(kin_vpar_apar_ps)) deallocate(kin_vpar_apar_ps)
      if(allocated(kin_vpar_apar_tr)) deallocate(kin_vpar_apar_tr)
      if(allocated(kin_tot_vpar_apar_ps)) deallocate(kin_tot_vpar_apar_ps)
      if(allocated(kin_tot_vpar_apar_tr)) deallocate(kin_tot_vpar_apar_tr)
      if(allocated(mag_en)) deallocate(mag_en)
    end if

    if(allocated(pot_en)) deallocate(pot_en)
    if(allocated(pot_tot_en)) deallocate(pot_tot_en)

    if(lcalc_kinenergy_trappas_sx) then
      if(allocated(kin_vpar_ps_sx)) deallocate(kin_vpar_ps_sx)
      if(allocated(kin_gradb_ps_sx)) deallocate(kin_gradb_ps_sx)

      if(allocated(kin_vpar_tr_sx)) deallocate(kin_vpar_tr_sx)
      if(allocated(kin_gradb_tr_sx)) deallocate(kin_gradb_tr_sx)
      if(allocated(pot_en_sx)) deallocate(pot_en_sx)
    end if

  end subroutine finalize

  !--------------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------------
  subroutine calc_kin_en()

    use control, only : nlapar
    use rotation, only : coriolis, cf_drift, cf_trap
    use grid,           only : nx, ns, nmu, nvpar, nsp, n_x_grid
    use grid,           only : nmod, number_of_species, gsp, gx, gs, proc_subset
    use dist,           only : fdisi, iapar
    use geom,           only : ints, bn, efun, ffun
    use geom,           only : hfun, dfun, ifun, metric
    use mode,           only : krho, kxrh, mode_label_G
    use components,     only : vthrat, signz, mas, veta_prime
    use components,     only : tgrid, veta, de
    use velocitygrid,   only : intmu, intvp, mugr, vpgr
    use global,         only : r_tiny
    use constants,      only : ci1
    use matdat,         only : get_f_from_g 
    use fields,         only : get_averaged_phi, get_averaged_apar
    use linear_terms,   only : drift
    use rotation,       only : vcor
    use rotation,       only : dcfphi_dpsi, dcfphi_ds
    use mpiinterface,   only : mpireduce_sum_inplace, root_processor
    use mpicomms, only : COMM_S_NE_X_NE
    use io,             only : append_chunk, xy_fmt, ascii_fmt
    use diagnos_generic, only : dfieldds, dfielddx
    use diagnos_generic, only : get_tr_ps_mask
    use global, only : PHI_GA_FIELD, APAR_FIELD
    use index_function, only : indx
    use diagnos_growth_freq, only : growth_rates, real_frequency
    use linear_terms, only : lvdgradf, lvpar_grad_df

    ! integers for the loop over all grid points 
    integer :: imod, ix, i, j, k, is 

    ! the gyro-averaged field
    complex :: phi_ga, apar_ga, apar
    ! and their deviation to s
    complex :: dphi_gads(ns)
    complex :: dapardx(nx)

    complex :: fdis
    real    :: ED
    integer :: mask
    real :: drift_x, drift_y
    complex :: dum, dum_energy

    ! The global species (ix) index is in isglb (ixg)
    integer :: isglb, ixg, ig


    if(.not. lcalc_kinenergy_trappas) return

    kin_vpar_ps = 0.
    kin_gradb_ps = 0.
    kin_vpar_tr = 0.
    kin_gradb_tr = 0.

    if(coriolis) then
      kin_coriolis_ps = 0.
      kin_coriolis_tr = 0.
    end if
    if(cf_drift) then
      kin_centrifugal_ps = 0.
      kin_centrifugal_tr = 0.
    end if
    ! if(cf_trap) then
      kin_centpot_ps = 0.
      kin_centpot_tr = 0.
    ! end if
    if(nlapar) then
      kin_betaprime_ps = 0.
      kin_betaprime_tr = 0.
      kin_vpar_apar_ps = 0.
      kin_vpar_apar_tr = 0.
      mag_en = 0.
      mag_en2 = 0.
    end if
    pot_en = 0.
    if(lcalc_kinenergy_trappas_sx) then
      kin_vpar_ps_sx = 0.
      kin_vpar_tr_sx = 0.
      kin_gradb_ps_sx = 0.
      kin_gradb_tr_sx = 0.
      pot_en_sx = 0.
    end if

    !  Calculates the potential energy and all the contributions to 
    ! time variation of kinetic energy
    nmod1: do imod = 1, nmod 
      ns3: do i = 1, ns
        ! the global s index
        ig = gs(i)

        ! FIXME the mu index j=1 and the species=1 is passed but
        ! unused, because Apar (without gyroavg) does not depend on mu
        ! or species
        call dfielddx(APAR_FIELD,imod,i,1,1,dapardx)

        nx1: do ix = 1, nx
          ! the actual (global) x index - blanks are left elsewhere
          ixg = gx(ix)
          nsp1: do is = 1, nsp
            ! the actual (global) species index - blanks are left elsewhere
            isglb = gsp(is)
            ! Integral over the velocity space 
            nmu3: do j = 1, nmu
              call dfieldds(PHI_GA_FIELD,imod,ix,j,is,dphi_gads)

              nvpar3: do k = 1, nvpar

                ! the gyro-averaged fields
                phi_ga  = get_averaged_phi(imod,ix,i,j,is,fdisi)
                apar_ga = get_averaged_apar(imod,ix,i,j,is,fdisi)

                ! fdis is the distribution without A_par contribution  
                fdis = de(ix,is) * get_f_from_g(imod,ix,i,j,k,is,fdisi)

                ! mask for trapped or passing particles
                mask = get_tr_ps_mask(ix,i,j,k,is)

                ! in the implicit scheme fdis can be NaN for intvp = 0 
                if (abs(intvp(i,j,k,is)) < r_tiny) fdis = 0.


                ! Time variation of kinetic energy due to the parallel
                ! velocity of trapped particles
                dum = fdis*conjg(-dphi_gads(i)) &
                   & *ffun(ix,i)*vthrat(is)*vpgr(i,j,k,is)*bn(ix,i)* &
                   & intvp(i,j,k,is)*intmu(j)*ints(i)*signz(is)*de(ix,is)
                if(lvpar_grad_df) then
                  kin_vpar_tr(imod,ixg,isglb)=kin_vpar_tr(imod,ixg,isglb)+dum*mask
                  kin_vpar_tr_sx(ig,ixg,isglb)=kin_vpar_tr_sx(ig,ixg,isglb)+dum*mask

                  ! Time variation of kinetic energy due to the parallel
                  ! velocity of passing particles
                  kin_vpar_ps(imod,ixg,isglb)=kin_vpar_ps(imod,ixg,isglb) + dum*(-mask+1)
                  kin_vpar_ps_sx(ig,ixg,isglb)=kin_vpar_ps_sx(ig,ixg,isglb) + dum*(-mask+1)
                end if

                if(lvdgradf) then

                  ED = vpgr(i,j,k,is)**2 + bn(ix,i)*mugr(j)
                  drift_x = tgrid(is)*ED*dfun(ix,i,1)
                  drift_y = tgrid(is)*ED*dfun(ix,i,2)

                ! Time variation of kinetic energy due to the gradB
                ! drift velocity of trapped particles
                dum = (fdis ) *conjg(-phi_ga*ci1)* &
                   & (drift_x*kxrh(ix)+drift_y*krho(imod))*bn(ix,i)*intvp(i,j,k,is) &
                   & *intmu(j)*ints(i)*de(ix,is)
                kin_gradb_tr(imod,ixg,isglb)=kin_gradb_tr(imod,ixg,isglb)+dum*mask
                kin_gradb_tr_sx(ig,ixg,isglb)=kin_gradb_tr_sx(ig,ixg,isglb)+dum*mask

                ! Time variation of kinetic energy due to the gradB
                ! drift velocity of passing particles
                kin_gradb_ps(imod,ixg,isglb)=kin_gradb_ps(imod,ixg,isglb)+dum*(-mask+1)
                kin_gradb_ps_sx(ig,ixg,isglb)=kin_gradb_ps_sx(ig,ixg,isglb)+dum*(-mask+1)


                if(nlapar) then
                  drift_x = tgrid(is)*vpgr(i,j,k,is)**2*veta_prime(ix)*efun(ix,i,1,1)/ bn(ix,i)**2
                  drift_y = tgrid(is)*vpgr(i,j,k,is)**2*veta_prime(ix)*efun(ix,i,1,2)/ bn(ix,i)**2

                  dum = fdis*conjg(-phi_ga*ci1)* &
                     & (drift_x*kxrh(ix)+drift_y*krho(imod))*bn(ix,i)*intvp(i,j,k,is) &
                     & *intmu(j)*ints(i) *de(ix,is)
                  ! Time variation of kinetic energy due to the
                  ! Betaprime drift velocity of trapped particles
                  kin_betaprime_tr(imod,ixg,isglb)=kin_betaprime_tr(imod,ixg,isglb) + &
                     & dum *mask

                  ! Time variation of kinetic energy due to the
                  ! Betaprime drift velocity of passing particles
                  kin_betaprime_ps(imod,ixg,isglb)=kin_betaprime_ps(imod,ixg,isglb)+ &
                     & dum * (-mask+1)

                end if

                if(coriolis) then
                  drift_x = 2.E0*mas(is)*vthrat(is)*vpgr(i,j,k,is)*vcor*hfun(ix,i,1)
                  drift_y = 2.E0*mas(is)*vthrat(is)*vpgr(i,j,k,is)*vcor*hfun(ix,i,2)

                  dum = (fdis*conjg(-phi_ga*ci1))* &
                     & (drift_x*kxrh(ix)+drift_y*krho(imod))*bn(ix,i)*intvp(i,j,k,is) &
                     & *intmu(j)*ints(i)*de(ix,is)
                  ! Time variation of kinetic energy due to the
                  ! coriolis drift velocity of trapped particles
                  kin_coriolis_tr(imod,ixg,isglb)=kin_coriolis_tr(imod,ixg,isglb)+dum * mask 

                  ! Time variation of kinetic energy due to the
                  ! coriolis drift velocity of passing particles
                  kin_coriolis_ps(imod,ixg,isglb)=kin_coriolis_ps(imod,ixg,isglb)+dum * (-mask+1)

                end if

                if(cf_drift) then
                  drift_x = vcor*vcor*mas(is)*ifun(ix,i,1)
                  drift_y = vcor*vcor*mas(is)*ifun(ix,i,2)
                  dum = (fdis*conjg(-phi_ga*ci1))* &
                     & (drift_x*kxrh(ix)+drift_y*krho(imod))*bn(ix,i)*intvp(i,j,k,is) &
                     & *intmu(j)*ints(i) *de(ix,is)
                  
                  ! Time variation of kinetic energy due to the
                  ! centrifugal drift velocity of trapped particles
                  kin_centrifugal_tr(imod,ixg,isglb)=kin_centrifugal_tr(imod,ixg,isglb)+&
                     & dum*mask 

                  ! Time variation of kinetic energy due to the
                  ! centrifugal drift velocity of passing particles
                  kin_centrifugal_ps(imod,ixg,isglb)=kin_centrifugal_ps(imod,ixg,isglb)+&
                     & dum*(-mask+1)

                end if

                ! if(cf_trap) then
                  drift_x = efun(ix,i,1,1)*dcfphi_dpsi(i)+efun(ix,i,3,1)*dcfphi_ds(i)
                  drift_y = efun(ix,i,1,2)*dcfphi_dpsi(i)+efun(ix,i,3,2)*dcfphi_ds(i)

                  dum = (fdis*conjg(-phi_ga*ci1))* &
                     & (drift_x*kxrh(ix)+drift_y*krho(imod))*bn(ix,i)*intvp(i,j,k,is) &
                     & *intmu(j)*ints(i)*signz(is) *de(ix,is)
                  ! Time variation of kinetic energy due to the
                  ! centrifugal potential drift of trapped particles
                  kin_centpot_tr(imod,ixg,isglb)=kin_centpot_tr(imod,ixg,isglb)+ dum * mask

                  ! ... of passing particles
                  kin_centpot_ps(imod,ixg,isglb)=kin_centpot_ps(imod,ixg,isglb)+&
                     & dum * (-mask+1)

                ! end if
                end if
                
                if(lvpar_grad_df) then
                  if(nlapar) then
                    dum = fdis*conjg(-apar_ga) &
                       & *vthrat(is)*vpgr(i,j,k,is)*bn(ix,i)* &
                       & intvp(i,j,k,is)*intmu(j)*ints(i)*signz(is)*de(ix,is) &
                       & *2
                    ! this factor 2 appears due to normalisation:
                    ! m_ref v_thref^2 = 2 T_ref

                    kin_vpar_apar_tr(imod,ixg,isglb)=kin_vpar_apar_tr(imod,ixg,isglb)+&
                       dum*mask
                    kin_vpar_apar_ps(imod,ixg,isglb)=kin_vpar_apar_ps(imod,ixg,isglb)+&
                       & dum*(-mask+1)
                  end if
                end if

                ! This is the same as ene_e2 from diagnos_energetics.
                dum = intmu(j) * 0.5&
                   & * bn(ix,i)*intvp(i,j,k,is) &
                   & *ints(i)*fdis*conjg(phi_ga)*signz(is)*de(ix,is)
                pot_en(imod,ixg,isglb) = pot_en(imod,ixg,isglb) + dum
                pot_en_sx(ig,ixg,isglb) = pot_en_sx(ig,ixg,isglb) + dum

                if(nlapar) then
                  ! This is the same as ene_m from diagnos_energetics.
                  dum =  &
                     & ints(i) * bn(ix,i)*intvp(i,j,k,is) * intmu(j)&
                     &*signz(is) &
                     &*fdis*conjg(apar_ga)*vthrat(is)*vpgr(i,j,k,is)*de(ix,is)
                  mag_en2(imod,ixg) = mag_en2(imod,ixg) + dum
                end if


              end do nvpar3
            end do nmu3
          end do nsp1
          nlapar1: if(nlapar .and. proc_subset(0,0,1,1,1)) then
            apar = fdisi(indx(iapar,imod,ix,i))

            ! integrate 1/(2 mu_0) |grad
            ! apar_ga|^2 .  Note that 1/(2 mu_0) is contained
            ! in beta_ref, i.e. veta here.
            dum = &
               & ints(i)* (metric(ix,i,1,1)*abs(dapardx(ix))**2 + &
               & metric(ix,i,1,2)*(ci1*krho(imod)*apar)*conjg(dapardx(ix)) + &
               & metric(ix,i,2,1)*dapardx(ix)*conjg(ci1*krho(imod)*apar) + &
               & metric(ix,i,2,2)*abs(ci1*krho(imod)*apar)**2)/veta(ix)
            mag_en(imod,ixg) = mag_en(imod,ixg) + dum

          end if nlapar1
        end do nx1
      end do ns3
    end do nmod1

    ! The 3D arrays ([k]y, k[x], species)
    ! are MPI reduced over s and velocity and "gathered" over x and species
    call mpireduce_sum_inplace(kin_vpar_ps, shape(kin_vpar_ps))
    call mpireduce_sum_inplace(kin_vpar_tr, shape(kin_vpar_tr))
    call mpireduce_sum_inplace(kin_gradb_ps, shape(kin_gradb_ps))
    call mpireduce_sum_inplace(kin_gradb_tr, shape(kin_gradb_tr))

    call mpireduce_sum_inplace(kin_vpar_ps_sx, shape(kin_vpar_ps_sx))
    call mpireduce_sum_inplace(kin_vpar_tr_sx, shape(kin_vpar_tr_sx))
    call mpireduce_sum_inplace(kin_gradb_ps_sx, shape(kin_gradb_ps_sx))
    call mpireduce_sum_inplace(kin_gradb_tr_sx, shape(kin_gradb_tr_sx))
    if(coriolis) then
      call mpireduce_sum_inplace(kin_coriolis_ps, shape(kin_coriolis_ps))
      call mpireduce_sum_inplace(kin_coriolis_tr, shape(kin_coriolis_tr))
    end if
    if(cf_drift) then
      call mpireduce_sum_inplace(kin_centrifugal_ps, shape(kin_centrifugal_ps))
      call mpireduce_sum_inplace(kin_centrifugal_tr, shape(kin_centrifugal_tr))
    end if
    ! if(cf_trap) then
      call mpireduce_sum_inplace(kin_centpot_ps, shape(kin_centpot_ps))
      call mpireduce_sum_inplace(kin_centpot_tr, shape(kin_centpot_tr))
    ! end if
    if(nlapar) then
      call mpireduce_sum_inplace(kin_betaprime_ps, shape(kin_betaprime_ps))
      call mpireduce_sum_inplace(kin_betaprime_tr, shape(kin_betaprime_tr))

      call mpireduce_sum_inplace(kin_vpar_apar_ps, shape(kin_vpar_apar_ps))
      call mpireduce_sum_inplace(kin_vpar_apar_tr, shape(kin_vpar_apar_tr))
      call mpireduce_sum_inplace(mag_en2, shape(mag_en2))
      if(proc_subset(0,0,1,1,1)) then
        call mpireduce_sum_inplace(mag_en, shape(mag_en), COMM_S_NE_X_NE)
      end if
    end if
    call mpireduce_sum_inplace(pot_en, shape(pot_en))
    call mpireduce_sum_inplace(pot_en_sx, shape(pot_en_sx))
    
    if(root_processor) then
      kin_tot_vpar_ps = 0.
      kin_tot_vpar_tr=0.
      kin_tot_gradb_ps=0.
      kin_tot_gradb_tr=0.
      if(coriolis) then
        kin_tot_coriolis_ps=0.
        kin_tot_coriolis_tr=0.
      end if
      if(cf_drift) then
        kin_tot_centrifugal_ps=0.
        kin_tot_centrifugal_tr=0.
      end if
      ! if(cf_trap) then
        kin_tot_centpot_ps=0.
        kin_tot_centpot_tr=0.
      ! end if
      if(nlapar) then
        kin_tot_betaprime_ps=0.
        kin_tot_betaprime_tr=0.

        kin_tot_vpar_apar_ps=0.
        kin_tot_vpar_apar_tr=0.
      end if
      pot_tot_en = 0.
      mag_tot_en = 0.
      mag_tot_en2 = 0.

      do imod = 1, nmod ! FIXME experimental
        do ix = 1, n_x_grid
          if(nlapar) then
            mag_tot_en = mag_tot_en + mag_en(imod,ix)
            mag_tot_en2 = mag_tot_en2 + mag_en2(imod,ix)
          end if
        end do
      end do

      species: do is = 1, number_of_species 
        do imod = 1, nmod
          do ix = 1, n_x_grid
            kin_tot_vpar_ps(is) = kin_tot_vpar_ps(is) + &
               & kin_vpar_ps(imod,ix,is)
            kin_tot_vpar_tr(is) = kin_tot_vpar_tr(is) + &
               & kin_vpar_tr(imod,ix,is)
            kin_tot_gradb_ps(is) = kin_tot_gradb_ps(is) + &
               & kin_gradb_ps(imod,ix,is)
            kin_tot_gradb_tr(is) = kin_tot_gradb_tr(is) + &
               & kin_gradb_tr(imod,ix,is)
            if(coriolis) then
              kin_tot_coriolis_ps(is) = kin_tot_coriolis_ps(is) + &
                 & kin_coriolis_ps(imod,ix,is)
              kin_tot_coriolis_tr(is) = kin_tot_coriolis_tr(is) + &
                 & kin_coriolis_tr(imod,ix,is)
            end if
            if(cf_drift) then
              kin_tot_centrifugal_ps(is) = kin_tot_centrifugal_ps(is) + &
                 & kin_centrifugal_ps(imod,ix,is)
              kin_tot_centrifugal_tr(is) = kin_tot_centrifugal_tr(is) + &
                 & kin_centrifugal_tr(imod,ix,is)
            end if
            ! if(cf_trap) then
              kin_tot_centpot_ps(is) = kin_tot_centpot_ps(is) + &
                 & kin_centpot_ps(imod,ix,is)
              kin_tot_centpot_tr(is) = kin_tot_centpot_tr(is) + &
                 & kin_centpot_tr(imod,ix,is)
            ! end if
            pot_tot_en(is) = pot_tot_en(is) + &
               & pot_en(imod,ix,is)

            kin_tot_betaprime_ps(is) = kin_tot_betaprime_ps(is) + &
               & kin_betaprime_ps(imod,ix,is)
            kin_tot_betaprime_tr(is) = kin_tot_betaprime_tr(is) + &
               & kin_betaprime_tr(imod,ix,is)
            if(nlapar) then
              kin_tot_vpar_apar_ps(is) = kin_tot_vpar_apar_ps(is) + &
                 & kin_vpar_apar_ps(imod,ix,is)*conjg(-ci1*real_frequency(mode_label_G(imod,ix)) + &
                 & growth_rates(mode_label_G(imod,ix)))
              kin_tot_vpar_apar_tr(is) = kin_tot_vpar_apar_tr(is) + &
                 & kin_vpar_apar_tr(imod,ix,is)*conjg(-ci1*real_frequency(mode_label_G(imod,ix)) + &
                 & growth_rates(mode_label_G(imod,ix)))
            end if
          end do !n_x_grid
        end do !nmod
      end do species

      ! The result of this diagnostic are the contributions to the
      ! growthrate gamma:
      !
      ! gamma = -dEkin/dt / (2 (Epot + Emag))
      !       = gamma_contrib_vpar + gamma_contrib_gradb+...+gamma_contrib_vpar_apar
      !
      ! to make it this simple (see [eq.3] in the header comment),
      ! divide with the pot.&mag. energy and take the real part:

      if(mag_en_with_fdis) then
        dum_energy = 2*(sum(pot_tot_en) + mag_tot_en2)
      else
        dum_energy = 2*(sum(pot_tot_en) + mag_tot_en)
      end if

      kin_tot_vpar_ps = -real(kin_tot_vpar_ps)/real(dum_energy)
      kin_tot_vpar_tr = -real(kin_tot_vpar_tr)/real(dum_energy)
      kin_tot_gradb_ps = -real(kin_tot_gradb_ps)/real(dum_energy)
      kin_tot_gradb_tr = -real(kin_tot_gradb_tr)/real(dum_energy)
      if(nlapar) then
        kin_tot_betaprime_ps = -real(kin_tot_betaprime_ps)/real(dum_energy)
        kin_tot_betaprime_tr = -real(kin_tot_betaprime_tr)/real(dum_energy)
      end if
      if(coriolis) then
        kin_tot_coriolis_ps = -real(kin_tot_coriolis_ps)/real(dum_energy)
        kin_tot_coriolis_tr = -real(kin_tot_coriolis_tr)/real(dum_energy)
      end if
      if(cf_drift) then
        kin_tot_centrifugal_ps = -real(kin_tot_centrifugal_ps)/real(dum_energy)
        kin_tot_centrifugal_tr = -real(kin_tot_centrifugal_tr)/real(dum_energy)
      end if
      ! if(cf_trap) then
        kin_tot_centpot_ps = -real(kin_tot_centpot_ps)/real(dum_energy)
        kin_tot_centpot_tr = -real(kin_tot_centpot_tr)/real(dum_energy)
      ! end if
      if(nlapar) then
        kin_tot_vpar_apar_ps = -real(kin_tot_vpar_apar_ps)/real(dum_energy)
        kin_tot_vpar_apar_tr = -real(kin_tot_vpar_apar_tr)/real(dum_energy)
      end if
    end if

  end subroutine calc_kin_en


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_regular()
    use control, only : nlapar
    use rotation, only : coriolis, cf_drift!, cf_trap
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use mpiinterface, only : root_processor

    if(.not. lcalc_kinenergy_trappas) return

    call calc_kin_en()

    if(root_processor) then
      call append_chunk(i_pot_en, real(pot_tot_en), xy_fmt, ascii_fmt)
      call append_chunk(i_vpar_ps, real(kin_tot_vpar_ps), xy_fmt, ascii_fmt)
      call append_chunk(i_vpar_tr, real(kin_tot_vpar_tr), xy_fmt, ascii_fmt)
      call append_chunk(i_gradb_ps, real(kin_tot_gradb_ps), xy_fmt, ascii_fmt)
      call append_chunk(i_gradb_tr, real(kin_tot_gradb_tr), xy_fmt, ascii_fmt)

      if(coriolis) then
        call append_chunk(i_coriolis_ps, real(kin_tot_coriolis_ps), xy_fmt, ascii_fmt)
        call append_chunk(i_coriolis_tr, real(kin_tot_coriolis_tr), xy_fmt, ascii_fmt)
      end if
      if(cf_drift) then
        call append_chunk(i_centrifugal_ps, real(kin_tot_centrifugal_ps), xy_fmt, ascii_fmt)
        call append_chunk(i_centrifugal_tr, real(kin_tot_centrifugal_tr), xy_fmt, ascii_fmt)
      end if
      ! if(cf_trap) then
        call append_chunk(i_centpot_ps, real(kin_tot_centpot_ps), xy_fmt, ascii_fmt)
        call append_chunk(i_centpot_tr, real(kin_tot_centpot_tr), xy_fmt, ascii_fmt)
      ! end if
      if(nlapar) then
        call append_chunk(i_betaprime_ps, real(kin_tot_betaprime_ps), xy_fmt, ascii_fmt)
        call append_chunk(i_betaprime_tr, real(kin_tot_betaprime_tr), xy_fmt, ascii_fmt)
        call append_chunk(i_vpar_apar_ps, real(kin_tot_vpar_apar_ps), xy_fmt, ascii_fmt)
        call append_chunk(i_vpar_apar_tr, real(kin_tot_vpar_apar_tr), xy_fmt, ascii_fmt)
        call append_chunk(i_mag_en, (/ real(mag_tot_en), real(mag_tot_en2) /), xy_fmt, ascii_fmt)
      end if

      if(lcalc_kinenergy_trappas_sx) then
        call append_chunk(i_pot_en_sx, real(pot_en_sx), xy_fmt, ascii_fmt)
        call append_chunk(i_vpar_ps_sx, real(kin_vpar_ps_sx), xy_fmt, ascii_fmt)
        call append_chunk(i_vpar_tr_sx, real(kin_vpar_tr_sx), xy_fmt, ascii_fmt)
        call append_chunk(i_gradb_ps_sx, real(kin_gradb_ps_sx), xy_fmt, ascii_fmt)
        call append_chunk(i_gradb_tr_sx, real(kin_gradb_tr_sx), xy_fmt, ascii_fmt)
      end if

      if(lcalc_kinenergy_trappas_xy) then
        call append_chunk(i_pot_en_xy, real(pot_en), xy_fmt, ascii_fmt)
        call append_chunk(i_vpar_ps_xy, real(kin_vpar_ps), xy_fmt, ascii_fmt)
        call append_chunk(i_vpar_tr_xy, real(kin_vpar_tr), xy_fmt, ascii_fmt)
        call append_chunk(i_gradb_ps_xy, real(kin_gradb_ps), xy_fmt, ascii_fmt)
        call append_chunk(i_gradb_tr_xy, real(kin_gradb_tr), xy_fmt, ascii_fmt)
      end if
    end if

  end subroutine output_regular

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output(file_count)
    integer, intent(in) :: file_count

    ! To keep the compiler quiet.
    if (file_count > 0) continue

    if(.not. lcalc_kinenergy_trappas) return
    
    call output_regular()
  end subroutine output

end module diagnos_kinenergy_trappas
