!------------------------------------------------------------------------------
!> Energetics and entropy diagnostics:
!>
!>  1) The free energies spectra
!>     As in Scott PoP 17, 102306 (2010)
!>  2) The various terms in the entropy balance
!>     As in Candy PoP 13, 032310 (2006)
!>  3) The total energy (conserved only in global simulations 
!>     with parallel velocity nonlinearity 
!>
!> to 1)
!> Calculate the free energies spectra 
!> Column 1 is kinetic free energy         (ene_f)   Eq 40
!> Column 2 is ExB energy                  (ene_e)   Eq 41
!> The polarisation is done numerically and not with Gamma function
!> Column 3 is magnetic fluctuation energy (ene_m)   Eq 42
!> 
!> to 2)
!> Derivation of the terms is given in more detail in the GKW manual
!> and in the Diploma thesis of Stefan Grosshauser, Bayreuth.
!>
!>---- FURTHER NOTES ------------------------------------------------
!>
!> This diagnostic makes use of some results of the diagnos_fluxes
!> diagnostic and hence must be called after that one.
!------------------------------------------------------------------------------
module diagnos_energetics
  use control, only : lcalc_energetics

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem, calc_smallstep
  public :: output

  !> The switch for this diagnostic must be defined in the module control,
  !> for dependency reasons.
  !> The following line makes it possible
  !> to write 'use diagnos_energetics, only : lcalc_energetics' .
  public :: lcalc_energetics

  !> total energy switch 
  public :: lcalc_tot_energy 

  integer, save :: i_ene_e = -1, i_ene_m = -1, i_ene_f, i_ene_e2 = -1
  integer, save :: i_dt_entr = -1, i_dt_entr_field = -1
  integer, save :: i_dt_entr_field_adiacorr = -1
  integer, save :: i_dt_entr_tr = -1, i_dt_entr_ps = -1
  integer, save :: i_entr = -1, i_entr_field = -1, i_entr_field_adiacorr = -1
  integer, save :: i_entr_tr = -1, i_entr_ps = -1
  integer, save :: i_entr_coll = -1, i_entr_outflow = -1
  integer, save :: i_entr_num_dis = -1, i_entr_num_vp = -1, i_entr_num_perp =-1
  integer, save :: i_entr_temp_src = -1, i_entr_landau = -1
  integer, save :: i_entr_src01 = -1,i_entr_src02 = -1,i_entr_src03 = -1
  integer, save :: i_entr_src04 = -1,i_entr_src05 = -1,i_entr_src06 = -1
  integer, save :: i_entr_src_kyx = -1, i_entr_src_kyx_fsa = -1
  integer, save :: i_tot_energy = -1
  integer, save :: i_entr_kyx = -1, i_entr_kyx_fsa = -1
  
  !> time between call of calc_smallstep() and output()
  real, save :: delta_time_energetics = 1.
  
  !> last entr and last entr_field, for real calculation of the
  !> time derivative of the entropy.
  real, save, allocatable, dimension(:) :: last_entr
  real, save, allocatable, dimension(:,:) :: last_entr_tr, last_entr_ps
  complex, save, allocatable, dimension(:) :: last_entr_field
  complex, save, allocatable, dimension(:) :: last_entr_field_adiacorr

  !> Switch for the calculation of the total energy (not entropy) 
  logical, save :: lcalc_tot_energy 

  real, allocatable :: ene_f(:), ene_e(:), ene_m(:), ene_e2(:)
  real, allocatable :: dt_entr(:)
  complex, allocatable :: dt_entr_field(:)
  complex, allocatable :: dt_entr_field_adiacorr(:)
  real, allocatable :: dt_entr_tr(:,:)
  real, allocatable :: dt_entr_ps(:,:)
  real, allocatable :: entr(:), entr_kyx_local(:,:), entr_kyx_fsa_local(:,:)
  real, allocatable :: entr_ps(:,:), entr_tr(:,:), kysp_global_buf(:,:)
  complex, allocatable :: entr_field(:), entr_field_adiacorr(:)
  real, allocatable :: entr_num_dis(:), entr_num_vp(:), entr_num_perp(:)
  real, allocatable :: entr_coll(:)
  real, allocatable :: entr_src01(:), entr_src02(:), entr_src03(:)
  real, allocatable :: entr_src04(:), entr_src05(:), entr_src06(:)
  real, allocatable :: entr_src_kyx_local(:,:), entr_src_kyx_fsa_local(:,:)
  real, allocatable :: entr_temp_src(:), entr_outflow(:), entr_landau(:)
  complex, allocatable :: num_disp01(:), num_disp02(:), num_disp03(:)
  complex, allocatable :: collisionop(:), S(:), outflow(:)

  real :: tot_energy(2)
contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    lcalc_energetics = .false.
    lcalc_tot_energy = .false. 

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    call mpibcast(lcalc_energetics,1)
    call mpibcast(lcalc_tot_energy,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : nlapar, naverage, flux_tube
    use general, only : gkw_warn

    if (lcalc_tot_energy) then 
 
      ! if (.not. lpar_vel_nl) then 
      !   call gkw_warn('Total energy diagnostic is called without &
      !       & the parallel velocity non-linearity. Energy should &
      !       & not be conserved.')

    endif 

    if (.not.lcalc_energetics) return

    if (.not. flux_tube) then 
      call gkw_warn('The entropy balance diagnostic does not make much sense &
         & without the local limit.')
    end if
    if (nlapar) then
      call gkw_warn('The entropy balance is incomplete for the &
         & electromagnetic model.')
    end if
    if (naverage == 1) then
      call gkw_warn ('energetics diagnostic does not work with naverage == 1 &
         & because it needs to compute a d/dt!')
      lcalc_energetics = .false.
    end if
 
  end subroutine check
  
  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use control, only : nlapar, io_legacy, spectral_radius, disp_x, disp_vp
    use control, only : lcollisions
    use mode, only : mode_box
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use io, only : binary_fmt
    use diagnos_generic, only : attach_metadata_grid
    use global, only : PHI_FIELD, DISTRIBUTION, PHI_GA_FIELD
    use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS, VPAR_GHOSTCELLS
    use diagnos_generic, only : X_GHOSTCELLS, MU_GHOSTCELLS
    use grid, only : nmod, n_x_grid
    use mpiinterface, only : root_processor
    use global, only : r_tiny

    logical, intent(inout) :: requirements(:,:)
    character(len=7) :: repr

    if (lcalc_energetics .or. lcalc_tot_energy) then
      requirements(PHI_FIELD,LOCAL_DATA) = .true.
      requirements(DISTRIBUTION,LOCAL_DATA) = .true.
      requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
    end if
    
    if(lcalc_energetics) then
      requirements(DISTRIBUTION,S_GHOSTCELLS) = .true.
      requirements(PHI_GA_FIELD,S_GHOSTCELLS) = .true.
      if(lcollisions) then
        requirements(DISTRIBUTION,VPAR_GHOSTCELLS) = .true.
        requirements(DISTRIBUTION,MU_GHOSTCELLS) = .true.
      end if
      if(disp_vp > 0) requirements(DISTRIBUTION,VPAR_GHOSTCELLS) = .true.
      if (mode_box .and. abs(disp_x) > r_tiny .and..not. spectral_radius) then
        ! the perp dissipation is implemented with finite differences
        requirements(DISTRIBUTION,X_GHOSTCELLS) = .true.
      end if
    end if

    if(io_legacy) then
      repr = '.kyspec'
    else
      repr = '_ky'
    end if

    if(root_processor) then
      if (lcalc_tot_energy) then 
        call open_real_lu('total_energy', 'diagnostic/diagnos_energetics', (/ 2 /), &
           & ascii_fmt, i_tot_energy)
        call attach_metadata_grid(i_tot_energy, 'time', ascii_fmt)
        call attach_metadata(i_tot_energy, phys_unit_key, &
           & 'Normalized to background energy', ascii_fmt)
        call attach_metadata(i_tot_energy, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_tot_energy, comments_key, &
           & 'the energy of the system (integrated over the entire radial &
           & domain). The energy is split in kinetic and field part (they &
           & should balance). Both energies are normalized to the total &
           & kinetic energy in the background Maxwell. (energy levels much &
           & smaller than 1. then roughly indicate that f is small compared &
           & to F_M.', ascii_fmt)
      end if

      if (.not.lcalc_energetics) return

      call open_real_lu('ene_e'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_ene_e)
      call attach_metadata_grid(i_ene_e, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_ene_e, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_ene_e, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_ene_e, comments_key, not_avail, ascii_fmt)

      if (nlapar) then
        call open_real_lu('ene_m'//repr, 'diagnostic/diagnos_energetics', &
           & (/ nmod /), ascii_fmt, i_ene_m)
        call attach_metadata_grid(i_ene_m, 'time', 'krho', ascii_fmt)
        call attach_metadata(i_ene_m, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_ene_m, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_ene_m, comments_key, not_avail, ascii_fmt)
      end if

      call open_real_lu('ene_e2'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_ene_e2)
      call attach_metadata_grid(i_ene_e2, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_ene_e2, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_ene_e2, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_ene_e2, comments_key, not_avail, ascii_fmt)

      call open_real_lu('ene_f'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_ene_f)
      call attach_metadata_grid(i_ene_f, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_ene_f, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_ene_f, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_ene_f, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_entr)
      call attach_metadata_grid(i_entr, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr, description_key, &
         & 'entropy', ascii_fmt)
      call attach_metadata(i_entr, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_tr'//repr, 'diagnostic/diagnos_energetics', &
         & shape(kysp_global_buf), ascii_fmt, i_entr_tr)
      call attach_metadata_grid(i_entr_tr, 'time', 'krho', 'number_of_species', &
         & ascii_fmt)
      call attach_metadata(i_entr_tr, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_tr, description_key, &
         & 'entropy of the trapped population', ascii_fmt)
      call attach_metadata(i_entr_tr, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_ps'//repr, 'diagnostic/diagnos_energetics', &
         & shape(kysp_global_buf), ascii_fmt, i_entr_ps)
      call attach_metadata_grid(i_entr_ps, 'time', 'krho', 'number_of_species', &
         & ascii_fmt)
      call attach_metadata(i_entr_ps, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_ps, description_key, &
         & 'entropy of the passing population', ascii_fmt)
      call attach_metadata(i_entr_ps, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_kyx', 'diagnostic/diagnos_energetics', &
         & (/ nmod, n_x_grid /), binary_fmt, i_entr_kyx)
      if(spectral_radius) then
        call attach_metadata_grid(i_entr_kyx, 'time', &
           & 'krho', 'kxrh', binary_fmt)
      else
        call attach_metadata_grid(i_entr_kyx, 'time', &
           & 'krho', 'xphi', binary_fmt)
      end if
      call attach_metadata(i_entr_kyx, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', binary_fmt)
      call attach_metadata(i_entr_kyx, description_key, &
         & 'perpendicular entropy slice, at global s position xy_slice_ipar', &
         & binary_fmt)
      call attach_metadata(i_entr_kyx, comments_key, &
         & not_avail, binary_fmt)

      call open_real_lu('entr_kyx_fsa', 'diagnostic/diagnos_energetics', &
         & (/ nmod, n_x_grid /), binary_fmt, i_entr_kyx_fsa)
      if(spectral_radius) then
        call attach_metadata_grid(i_entr_kyx_fsa, 'time', &
           & 'krho', 'kxrh', binary_fmt)
      else
        call attach_metadata_grid(i_entr_kyx_fsa, 'time', &
           & 'krho', 'xphi', binary_fmt)
      end if
      call attach_metadata(i_entr_kyx_fsa, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', binary_fmt)
      call attach_metadata(i_entr_kyx_fsa, description_key, &
         & 'perpendicular entropy slice, flux surface averaged', &
         & binary_fmt)
      call attach_metadata(i_entr_kyx_fsa, comments_key, &
         & not_avail, binary_fmt)

      call open_real_lu('entr_src_kyx', 'diagnostic/diagnos_energetics', &
         & (/ nmod, n_x_grid /), binary_fmt, i_entr_src_kyx)
      if(spectral_radius) then
        call attach_metadata_grid(i_entr_src_kyx, 'time', &
           & 'krho', 'kxrh', binary_fmt)
      else
        call attach_metadata_grid(i_entr_src_kyx, 'time', &
           & 'krho', 'xphi', binary_fmt)
      end if
      call attach_metadata(i_entr_src_kyx, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', binary_fmt)
      call attach_metadata(i_entr_src_kyx, description_key, &
         & 'perpendicular slice of the entropy rate due to linear drive, at&
         & global s position xy_slice_ipar', binary_fmt)
      call attach_metadata(i_entr_src_kyx, comments_key, &
         & not_avail, binary_fmt)

      call open_real_lu('entr_src_kyx_fsa', 'diagnostic/diagnos_energetics', &
         & (/ nmod, n_x_grid /), binary_fmt, i_entr_src_kyx_fsa)
      if(spectral_radius) then
        call attach_metadata_grid(i_entr_src_kyx_fsa, 'time', &
           & 'krho', 'kxrh', binary_fmt)
      else
        call attach_metadata_grid(i_entr_src_kyx_fsa, 'time', &
           & 'krho', 'xphi', binary_fmt)
      end if
      call attach_metadata(i_entr_src_kyx_fsa, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', binary_fmt)
      call attach_metadata(i_entr_src_kyx_fsa, description_key, &
         & 'perpendicular slice of the entropy rate due to linear drive,&
         & flux surface averaged', binary_fmt)
      call attach_metadata(i_entr_src_kyx_fsa, comments_key, &
         & not_avail, binary_fmt)

      call open_real_lu('entr_field'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_entr_field)
      call attach_metadata_grid(i_entr_field, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_field, description_key, &
         & 'entropy in the field', ascii_fmt)
      call attach_metadata(i_entr_field, comments_key, not_avail, ascii_fmt)
      
      call open_real_lu('entr_field_adiacorr'//repr, &
         & 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_entr_field_adiacorr)
      call attach_metadata_grid(i_entr_field_adiacorr, 'time', 'krho', &
         & ascii_fmt)
      call attach_metadata(i_entr_field_adiacorr, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_field_adiacorr, description_key, &
         & 'a correction to entr_field, for adiabatic electrons', ascii_fmt)
      call attach_metadata(i_entr_field_adiacorr, comments_key, not_avail, &
         & ascii_fmt)
      
      call open_real_lu('dt_entr'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_dt_entr)
      call attach_metadata_grid(i_dt_entr, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_dt_entr, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_dt_entr, description_key, &
         & 'time derivative of the entropy', ascii_fmt)
      call attach_metadata(i_dt_entr, comments_key, not_avail, ascii_fmt)

      call open_real_lu('dt_entr_field'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_dt_entr_field)
      call attach_metadata_grid(i_dt_entr_field, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_dt_entr_field, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_dt_entr_field, description_key, &
         & 'time derivative of the entropy in the field', ascii_fmt)
      call attach_metadata(i_dt_entr_field, comments_key, not_avail, ascii_fmt)

      call open_real_lu('dt_entr_field_adiacorr'//repr, &
         & 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_dt_entr_field_adiacorr)
      call attach_metadata_grid(i_dt_entr_field_adiacorr, 'time', 'krho', &
         & ascii_fmt)
      call attach_metadata(i_dt_entr_field_adiacorr, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_dt_entr_field_adiacorr, description_key, &
         & 'time derivative of the correction term with adiabatic electrons', &
         & ascii_fmt)
      call attach_metadata(i_dt_entr_field_adiacorr, comments_key, not_avail, &
         & ascii_fmt)

      call open_real_lu('dt_entr_tr'//repr, 'diagnostic/diagnos_energetics', &
         & shape(kysp_global_buf), ascii_fmt, i_dt_entr_tr)
      call attach_metadata_grid(i_dt_entr_tr, 'time', 'krho', &
         & 'number_of_species', ascii_fmt)
      call attach_metadata(i_dt_entr_tr, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_dt_entr_tr, description_key, &
         & 'time derivative of the entropy, of the trapped population', ascii_fmt)
      call attach_metadata(i_dt_entr_tr, comments_key, not_avail, ascii_fmt)

      call open_real_lu('dt_entr_ps'//repr, 'diagnostic/diagnos_energetics', &
         & shape(kysp_global_buf), ascii_fmt, i_dt_entr_ps)
      call attach_metadata_grid(i_dt_entr_ps, 'time', 'krho', &
         & 'number_of_species', ascii_fmt)
      call attach_metadata(i_dt_entr_ps, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_dt_entr_ps, description_key, &
         & 'time derivative of the entropy, of the passing population', ascii_fmt)
      call attach_metadata(i_dt_entr_ps, comments_key, not_avail, ascii_fmt)


      call open_real_lu('entr_coll'//repr, 'diagnostic/diagnos_energetics', &
         & (/ nmod /), ascii_fmt, i_entr_coll)
      call attach_metadata_grid(i_entr_coll, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_coll, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_coll, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_coll, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_num_dis'//repr, 'diagnostic/diagnos_energetics',&
         & (/ nmod /), ascii_fmt, i_entr_num_dis)
      call attach_metadata_grid(i_entr_num_dis, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_num_dis, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_num_dis, description_key, &
         & 'Numerical dissipation from parallel derivatives', ascii_fmt)
      call attach_metadata(i_entr_num_dis, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_num_vp'//repr, 'diagnostic/diagnos_energetics',&
         & (/ nmod /), ascii_fmt, i_entr_num_vp)
      call attach_metadata_grid(i_entr_num_vp, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_num_vp, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_num_vp, description_key, &
         & 'Numerical dissipation from v_\parallel derivatives', ascii_fmt)
      call attach_metadata(i_entr_num_vp, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_num_perp'//repr, &
         & 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_num_perp)
      call attach_metadata_grid(i_entr_num_perp, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_num_perp, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_num_perp, description_key, &
         & 'Numerical dissipation from perpendicular derivatives', ascii_fmt)
      call attach_metadata(i_entr_num_perp, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_temp_src'//repr, &
         & 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_temp_src)
      call attach_metadata_grid(i_entr_temp_src, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_temp_src, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_temp_src, description_key, &
         & 'Entropy source/sink associated with the timedependent temperature source', &
         & ascii_fmt)
      call attach_metadata(i_entr_temp_src, comments_key, not_avail, ascii_fmt)
      
      call open_real_lu('entr_outflow'//repr, &
         & 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_outflow)
      call attach_metadata_grid(i_entr_outflow, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_outflow, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_outflow, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_outflow, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_source01'//repr, &
         & 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_src01)
      call attach_metadata_grid(i_entr_src01, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_src01, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_src01, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_src01, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_source02'//repr, 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_src02)
      call attach_metadata_grid(i_entr_src02, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_src02, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_src02, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_src02, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_source03'//repr, 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_src03)
      call attach_metadata_grid(i_entr_src03, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_src03, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_src03, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_src03, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_source04'//repr, 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_src04)
      call attach_metadata_grid(i_entr_src04, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_src04, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_src04, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_src04, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_source05'//repr, 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_src05)
      call attach_metadata_grid(i_entr_src05, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_src05, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_src05, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_src05, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_source06'//repr, 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_src06)
      call attach_metadata_grid(i_entr_src06, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_src06, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_src06, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_src06, comments_key, not_avail, ascii_fmt)

      call open_real_lu('entr_landau'//repr, 'diagnostic/diagnos_energetics', (/ nmod /), &
         & ascii_fmt, i_entr_landau)
      call attach_metadata_grid(i_entr_landau, 'time', 'krho', ascii_fmt)
      call attach_metadata(i_entr_landau, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}/R_{ref}', ascii_fmt)
      call attach_metadata(i_entr_landau, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_entr_landau, comments_key, not_avail, ascii_fmt)
    end if
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use grid, only : nmod, nx, nsp, number_of_species
    use general, only : gkw_abort
    use dist, only : nsolc
    use diagnos_fluxes_vspace, only : allocate_flux_det
    integer :: ierr 
    if (.not.lcalc_energetics) return

    allocate(last_entr(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: last_entr')
    allocate(last_entr_tr(nmod, nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: last_entr_tr')
    allocate(last_entr_ps(nmod, nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: last_entr_ps')
    allocate(last_entr_field(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: last_entr_field')
    allocate(last_entr_field_adiacorr(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: last_entr_field_adiacorr')

    allocate(entr_kyx_local(nmod,nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: entr_kyx_local')
    allocate(entr_src_kyx_local(nmod,nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: entr_src_kyx_local')

    allocate(entr_kyx_fsa_local(nmod,nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: entr_kyx_fsa_local')
    allocate(entr_src_kyx_fsa_local(nmod,nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: entr_src_kyx_fsa_local')

    call allocate_flux_det

    allocate(ene_f(nmod), ene_e(nmod), ene_m(nmod), ene_e2(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(dt_entr(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(dt_entr_field(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(dt_entr_field_adiacorr(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(dt_entr_tr(nmod, nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(dt_entr_ps(nmod, nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_tr(nmod, nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_ps(nmod, nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_field(nmod), entr_field_adiacorr(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_num_dis(nmod), entr_num_vp(nmod), entr_num_perp(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_coll(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_src01(nmod), entr_src02(nmod), entr_src03(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_src04(nmod), entr_src05(nmod), entr_src06(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_landau(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(entr_temp_src(nmod), entr_outflow(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(num_disp01(nsolc), num_disp02(nsolc), num_disp03(nsolc),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(collisionop(nsolc), S(nsolc), outflow(nsolc),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')

    allocate(kysp_global_buf(nmod, number_of_species),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    if(allocated(last_entr)) deallocate(last_entr)
    if(allocated(last_entr_tr)) deallocate(last_entr_tr)
    if(allocated(last_entr_ps)) deallocate(last_entr_ps)
    if(allocated(last_entr_field)) deallocate(last_entr_field)
    if(allocated(last_entr_field_adiacorr)) deallocate(last_entr_field_adiacorr)

    if(allocated(entr_kyx_local)) deallocate(entr_kyx_local)
    if(allocated(entr_src_kyx_local)) deallocate(entr_src_kyx_local)
    if(allocated(entr_kyx_fsa_local)) deallocate(entr_kyx_fsa_local)
    if(allocated(entr_src_kyx_fsa_local)) deallocate(entr_src_kyx_fsa_local)
    
    if(allocated(ene_f)) deallocate(ene_f)
    if(allocated(ene_e)) deallocate(ene_e)
    if(allocated(ene_m)) deallocate(ene_m)
    if(allocated(ene_e2)) deallocate(ene_e2)
    if(allocated(dt_entr)) deallocate(dt_entr)
    if(allocated(dt_entr_field)) deallocate(dt_entr_field)
    if(allocated(dt_entr_field_adiacorr)) deallocate(dt_entr_field_adiacorr)
    if(allocated(dt_entr_tr)) deallocate(dt_entr_tr)
    if(allocated(dt_entr_ps)) deallocate(dt_entr_ps)
    if(allocated(entr)) deallocate(entr)
    if(allocated(entr_field)) deallocate(entr_field)
    if(allocated(entr_field_adiacorr)) deallocate(entr_field_adiacorr)
    if(allocated(entr_tr)) deallocate(entr_tr)
    if(allocated(entr_ps)) deallocate(entr_ps)
    if(allocated(entr_num_dis)) deallocate(entr_num_dis)
    if(allocated(entr_num_vp)) deallocate(entr_num_vp)
    if(allocated(entr_num_perp)) deallocate(entr_num_perp)
    if(allocated(entr_coll)) deallocate(entr_coll)
    if(allocated(entr_src01)) deallocate(entr_src01)
    if(allocated(entr_src02)) deallocate(entr_src02)
    if(allocated(entr_src03)) deallocate(entr_src03)
    if(allocated(entr_src04)) deallocate(entr_src04)
    if(allocated(entr_src05)) deallocate(entr_src05)
    if(allocated(entr_src06)) deallocate(entr_src06)
    if(allocated(entr_landau)) deallocate(entr_landau)
    if(allocated(entr_temp_src)) deallocate(entr_temp_src)
    if(allocated(entr_outflow)) deallocate(entr_outflow)
    if(allocated(num_disp01)) deallocate(num_disp01)
    if(allocated(num_disp02)) deallocate(num_disp02)
    if(allocated(num_disp03)) deallocate(num_disp03)
    if(allocated(collisionop)) deallocate(collisionop)
    if(allocated(S)) deallocate(S)
    if(allocated(outflow)) deallocate(outflow)

    if(allocated(kysp_global_buf)) deallocate(kysp_global_buf)

  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_smallstep(i_smallstep)
    use control,          only : naverage, time
    integer, intent(in) :: i_smallstep
    real, save :: last_time_pre_energetics = 0.

    if (.not.lcalc_energetics) return

    ! if (i_smallstep == n_smallstep - 1) then
    ! else if (i_smallstep == n_smallstep) then
    ! end if

    ! FIXME This calculation of a time derivative presumably does
    ! not work for naverage == 1

    if (i_smallstep == naverage - 1) then
      call calc_entr(last_entr, last_entr_field, &
         & last_entr_field_adiacorr,last_entr_tr,last_entr_ps)
      ! store the time this subroutine was last called in a local variable
      last_time_pre_energetics = time

    else if (i_smallstep == naverage) then
      ! save time interval between this step and the
      ! last call of calc_entr().
      ! this intervall is needed by several subroutines of this module
      delta_time_energetics = time - last_time_pre_energetics
    end if

  end subroutine calc_smallstep
  
  !--------------------------------------------------------------------
  !> calc_energy calculates the energy of the system (integrated over
  !> the entire radial domain). The energy is split in kinetic and 
  !> field part (they should balance). Both energies are normalized 
  !> to the total kinetic energy in the background Maxwell. (energy 
  !> levels much smaller than 1. then roughly indicate that f is small
  !> compared to F_M.  
  !--------------------------------------------------------------------
  subroutine calc_tot_energy 

    use components,   only : rhostar, signz, tmp, tgrid, dgrid
    use components,   only : adiabatic_electrons, iadia, de
    use geom,         only : jacobian_G, bn, ints
    use dist,         only : fmaxwl, fdisi, phi
    use fields,       only : get_averaged_phi 
    use grid,         only : nmod, nx, ns, nmu, nvpar, nsp, gx 
    use matdat,       only : get_f_from_g
    use mode,         only : iyzero
    use velocitygrid, only : vpgr, intvp, mugr, intmu 
    use mpiinterface, only : mpiallreduce_sum 
    use mpicomms,     only : COMM_CART
    use mpiinterface,    only : root_processor

    integer       :: imod, ix, i, j, k, is 
    logical, save :: initialized = .false. 
    real, save    :: energy_back

    real    :: dfac, dum, energy_adi, phisq
    complex :: phiav, cdum
    !> Total kinetic energy in the perturbed distribution normalized to the
    !> total kinetic energy in the background
    complex :: energy_kin
    !> Total energy in the field normalized to the kinetic energy in the
    !> background
    real, save :: energy_fld

    if (.not. initialized) then 
      energy_back = 0.E0 
      ! integrate the kinetic energy  in the background distribution. 
      ! note that the grid distance in the radial direction is missing
      ! This is o.k. as long as it does not appear in any of the other 
      ! integrations and the grid is uniform in psi. Also the integration
      ! over zeta has been surpressed. 
      do ix = 1, nx; do i = 1, ns; do j = 1, nmu; 
        do k = 1, nvpar; do is = 1, nsp 
          energy_back = energy_back + ints(i)*jacobian_G(gx(ix))*intmu(j)* &
                      & bn(ix,i)*intvp(i,j,k,is)* tgrid(is) * dgrid(is) *  &
                      & (vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i))*       &
                      & fmaxwl(ix,i,j,k,is) 
        end do; end do 
      end do; end do; end do 
      initialized = .true.

      ! Global sum over 
      call mpiallreduce_sum(energy_back, dum, 1, COMM_CART) 
      energy_back = dum
      
      ! The output file (should be moved)
      if (root_processor) then 
      endif 
      
    endif 

    ! initialize the energy to zero 
    energy_kin = 0.E0  
    energy_fld = 0.E0 

    ! only the zero mode caries net kinetic energy
    imod = iyzero 

    ! if there is a zero mode in the system 
    if (imod.ne.0) then 
      do ix = 1, nx; do i = 1, ns; do j = 1, nmu; 
        do k = 1, nvpar; do is = 1, nsp  
          energy_kin = energy_kin + ints(i)*jacobian_G(gx(ix))*intmu(j)*    &
                     & bn(ix,i)*intvp(i,j,k,is)* tgrid(is) * dgrid(is) *    & 
                     & (vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(i,j)) *         &
                     & get_f_from_g(imod,ix,i,j,k,is,fdisi)
        end do; end do
      end do; end do; end do 
      energy_kin = rhostar*energy_kin / energy_back 

      ! global sum 
      call mpiallreduce_sum(energy_kin, cdum, 1, COMM_CART)
      energy_kin = cdum
      
    endif 

    do imod = 1, nmod; do ix = 1, nx; do i = 1, ns; 
      do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp 
         dfac = 2.E0; if (imod.eq.iyzero) dfac = 1.E0 
         energy_fld = energy_fld + dfac * signz(is)**2 * tgrid(is) *         &
                    & dgrid(is)* intvp(i,j,k,is) * intmu(j) * bn(ix,i)*      &
                    & ints(i) * jacobian_G(gx(ix)) *                         &
                    & (abs(phi(imod,ix,i))**2 -                              &
                    &  abs(get_averaged_phi(imod,ix,i,j,is,fdisi))**2)       &
                    & * fmaxwl(ix,i,j,k,is) / tmp(ix,is)
      end do; end do; end do; 
    end do; end do; end do 

    ! global sum 
    call mpiallreduce_sum(energy_fld, dum, 1, COMM_CART) 
    energy_fld = dum 
    
    ! For the adiabatic electron case there is only a contribution to 
    ! the field energy.  
    if (adiabatic_electrons) then
      
      energy_adi = 0.E0 
      
      do imod = 1, nmod; do ix = 1, nx 
        if (imod.eq.iyzero) then 
         phiav = 0. 
         phisq = 0. 
         do i = 1, ns 
           phiav = phiav + ints(i)*phi(imod,ix,i)
           phisq = phisq + ints(i)*abs(phi(imod,ix,i))**2 
         end do 
         energy_adi = energy_adi + de(ix,nsp+iadia) / tmp(ix,nsp+iadia) *    &
                    & jacobian_G(gx(ix)) * (phisq - abs(phiav)**2) /         & 
                    & (4.E0*tgrid(nsp+iadia))   
        else 
          do i = 1, ns 
            energy_adi = energy_adi + de(ix,nsp+iadia) / tmp(ix,nsp+iadia) * &
                       & ints(i)*jacobian_G(gx(ix))*                         &
                       & abs(phi(imod,ix,i))**2 / (2.E0 * tgrid(nsp+iadia)) 
          end do 
        endif 
      end do; end do 
      
      ! global sum 
      energy_fld = energy_fld + energy_adi 
      
    endif 
 
    energy_fld = rhostar**2*energy_fld / energy_back 

    tot_energy(1) = real(energy_kin)
    tot_energy(2) = energy_fld
    
  end subroutine calc_tot_energy

  !--------------------------------------------------------------------
  !> The routine calc_largestep() is to be called repeatedly, after every
  !> large timestep, and here the diagnostic should calculate its
  !> quantities.
  !> 
  !> Splitting calculation and output is to be encouraged, particularly when
  !>  - the output is higher dimensional, and requires MPI or MPI-IO
  !>  - same data is output in different dimensional slices
  !--------------------------------------------------------------------
  subroutine calc_largestep()
    use grid,             only : nmod, nx, ns, nvpar, nmu, nsp
    use grid,             only : gsp, gx, gs, proc_subset
    use fields,           only : get_averaged_phi, get_averaged_apar
    use matdat,           only : get_f_from_g
    use matdat,           only : matd, matvpd, matperpd, matcoll, matoutflow
    use matrix_format,    only : usmv
    use dist,             only : fdisi, fmaxwl, ifdis, fdis_tmp, iphi
    use geom,             only : ints, bn, lfun, ffun
    use velocitygrid,     only : intmu, intvp, vpgr
    use components,       only : signz, de, tmp, mas, vthrat, vp, tp, fp
    use index_function,   only : indx
    use constants,        only : c1
    use global,           only : r_tiny
    use normalise,        only : fnorm1d
    use rotation,         only : cfen, vcor
    use source_time,      only : add_source_time
    use matdat, only : add_source
    use control,          only : zonal_adiabatic, normalize_per_toroidal_mode
    use components,       only : adiabatic_electrons
    use mpiinterface,     only : mpireduce_sum_inplace, root_processor
    use mpicomms,         only : COMM_SP_NE_X_NE, COMM_X_EQ, COMM_S_EQ_X_EQ
    use mpicomms,         only : COMM_SP_EQ
    use diagnos_fluxes, only : pflux_es, eflux_es, vflux_es
    use diagnos_fluxes_vspace, only : calc_fluxes_full_detail
    use diagnos_fluxes_vspace, only : pflux_det, eflux_det, vflux_det
    use global,  only : PHI_GA_FIELD, PHI_FIELD
    use diagnos_generic,  only : parseval_correction, xy_slice_ipar, dfieldds

    integer :: i, j, k, is, imod, ix
    integer :: iih, ierr
    integer :: ixg, isglb
    complex :: fdis, phi, phi_ga, apar_ga
    real :: fdis_sq, phi_sq, phi_ga_sq
    real :: d2X, d3X, dumint, d3v, dum
    real :: acc_normfac(nmod), inverse_fnorm
    real :: acc_normfac_imod
    complex :: dphids(ns)

    if (.not.lcalc_energetics) return

    num_disp01 = 0.0
    num_disp02 = 0.0
    num_disp03 = 0.0
    collisionop = 0.0
    outflow = 0.0
    S = 0.0

    acc_normfac = 1.0

    ! If you want the result to be un-normalised (independently of the
    ! 'normalized' switch), then uncomment this line:
    ! acc_normfac = accumulated_normfactor

    ! dissipation elements are stored in the matrices matd, matvpd, matperpd.
    call usmv(c1,matd,fdis_tmp,num_disp01,ierr)
    call usmv(c1,matvpd,fdis_tmp,num_disp02,ierr)
    call usmv(c1,matperpd,fdis_tmp,num_disp03,ierr)
    call usmv(c1,matcoll,fdis_tmp,collisionop,ierr)
    call usmv(c1,matoutflow,fdis_tmp,outflow,ierr)
    
    call calc_entr(entr, entr_field, entr_field_adiacorr, &
       & entr_tr, entr_ps, &
       & entr_kyx_local, entr_kyx_fsa_local)
    entr = acc_normfac**2 * entr
    entr_field = acc_normfac**2 * entr_field
    entr_field_adiacorr = acc_normfac**2 * entr_field_adiacorr
    if(normalize_per_toroidal_mode) then
      do imod = 1,nmod
        entr_kyx_local(imod,:) = acc_normfac(imod)**2 * entr_kyx_local(imod,:)
        entr_kyx_fsa_local(imod,:) = acc_normfac(imod)**2 * entr_kyx_fsa_local(imod,:)
      end do
    else
      entr_kyx_local = acc_normfac(1)**2 * entr_kyx_local
      entr_kyx_fsa_local = acc_normfac(1)**2 * entr_kyx_fsa_local
    end if

    call add_source(S, 1.0)
    call add_source_time(S, 0.0, 1.0)
    call calc_fluxes_full_detail(PHI_FIELD, .false.)

    ene_f = 0.; ene_e = 0.; ene_e2 = 0.; ene_m = 0.
    entr_num_dis = 0.; entr_num_vp = 0.; entr_num_perp = 0.
    entr_coll = 0.; entr_temp_src = 0.; entr_outflow = 0.
    entr_src01 = 0.;entr_src02 = 0.;entr_src03 = 0.
    entr_src04 = 0.;entr_src05 = 0.;entr_src06 = 0.
    entr_src_kyx_local = 0.
    entr_src_kyx_fsa_local = 0.
    entr_landau = 0.

    do imod=1,nmod
      if(normalize_per_toroidal_mode) then
        acc_normfac_imod = acc_normfac(imod)
      else
        acc_normfac_imod = acc_normfac(1)
      end if

      d2X = 1.0

      do ix=1,nx
        ! the actual (global) x index - blanks are left elsewhere
        ixg = gx(ix)

        do is=1,nsp
          if(de(ix,is) < r_tiny) cycle
          
          ! the actual (global) species index - blanks are left elsewhere
          isglb = gsp(is)

          ! entropy sources (or sinks) from currents\cdot gradients:
          ! NOTE: The general diagnostic module *must* make the fluxes
          ! be calculated before this energetics diagnostic.
          ! The flux calculated there can therefore be used here.
          ! The fluxes contain products of the distribution and the
          ! es. potential (in v_\chi). The normalisation correction
          ! is therefore accumulated_normfactor**2 .
          ! The calculation of the fluxes involves already what is here called
          ! the parseval_correction. This is why this factor is not put into 
          ! the following expressions.
          entr_src01(imod) = entr_src01(imod) + d2X * acc_normfac_imod**2 &
             & * de(ix,is) &
             & * pflux_es(imod,ixg,isglb) * fp(ix,is)
          entr_src02(imod) = entr_src02(imod) + d2X * acc_normfac_imod**2 &
             & * de(ix,is) &
             & * mas(is)*vthrat(is)**2 * eflux_es(imod,ixg,isglb) * tp(ix,is)
          entr_src03(imod) = entr_src03(imod) - d2X * acc_normfac_imod**2 &
             & * de(ix,is) &
             & * 1.5 * pflux_es(imod,ixg,isglb) * tp(ix,is)
          ! Formally, cfen depends on parallel coord., but this
          ! balance is only valid in the local limit anyway. This is why
          ! cfen(1,is) can be used instead.
          ! entr_src04(imod) = entr_src04(imod) - cfen(i,is)/tmp(ix,is) ...
          entr_src04(imod) = entr_src04(imod) + d2X * acc_normfac_imod**2 &
             & * de(ix,is) &
             & * cfen(1,is)/tmp(ix,is) * pflux_es(imod,ixg,isglb) * vp(ix,is)
          !FIXME the factor 2 is questionable in term 05
          entr_src05(imod) = entr_src05(imod) + d2X * acc_normfac_imod**2 &
             & * de(ix,is) &
             & * mas(is) * vthrat(is) * vflux_es(imod,ixg,isglb) * vp(ix,is) &
             & * 2.0 / tmp(ix,is)
          entr_src06(imod) = entr_src06(imod) + d2X * acc_normfac_imod**2 &
             & * de(ix,is) &
             & * mas(is)/tmp(ix,is) * vcor*vcor * pflux_es(imod,ixg,isglb) * lfun

          do i=1,ns
            phi  = acc_normfac_imod * fdisi(indx(iphi,imod,ix,i))
            phi_sq = abs(phi)**2
            

            ! ints has n_s_grid elements but in parallelize_geom some elements
            ! are copied so that every processor can index his by ints(1:ns)
            d3X = ints(i) * d2X 

            do j=1,nmu
              phi_ga  = acc_normfac_imod * get_averaged_phi(imod,ix,i,j,is,fdisi)
              phi_ga_sq = abs(phi_ga)**2
              apar_ga = acc_normfac_imod * get_averaged_apar(imod,ix,i,j,is,fdisi)
              call dfieldds(PHI_GA_FIELD,imod,ix,j,is,dphids)

              ! bn is the magnetic field modulus. The 2*pi which is
              ! furthermore contained in the velocity-space Jacobian is
              ! defined into intmu.
              dumint = intmu(j)*bn(ix,i)

              do k=1,nvpar
                ! fdisi stores g, but this function returns f:
                fdis = acc_normfac_imod * get_f_from_g(imod,ix,i,j,k,is,fdisi)
                fdis_sq = abs(fdis)**2

                d3v = dumint * intvp(i,j,k,is)

                ene_f(imod) = ene_f(imod) + parseval_correction(imod) * &
                   & de(ix,is) * & 
                   & 0.5 * fdis_sq * tmp(ix,is) / fmaxwl(ix,i,j,k,is) * d3v * d3X
                ene_e(imod) = ene_e(imod) + parseval_correction(imod) * &
                   & 0.5 * signz(is)**2 * (phi_sq - phi_ga_sq) * &
                   & (de(ix,is) * fmaxwl(ix,i,j,k,is)/ tmp(ix,is)**2) * d3v * d3X
                ene_m(imod) = ene_m(imod) + parseval_correction(imod) * &
                   & de(ix,is) * &
                   & 0.5 * real(fdis*conjg(apar_ga))*vpgr(i,j,k,is)*signz(is)*d3v * d3X
                ! should be equivalent to ene_e - cross check
                ene_e2(imod) = ene_e2(imod) + parseval_correction(imod) * &
                   & de(ix,is) * &
                   & 0.5 * real(fdis*conjg(phi_ga))*signz(is)*d3v * d3X

                iih = indx(ifdis,imod,ix,i,j,k,is) 
                ! collisional dissipation
                entr_coll(imod) = entr_coll(imod) + acc_normfac_imod**2 &
                   & * parseval_correction(imod) * &
                   & de(ix,is) * &
                   & real(conjg(fdis)/fmaxwl(ix,i,j,k,is) &
                   & * collisionop(iih)) * d3v * d3X

                ! numerical parallel dissipation:
                entr_num_dis(imod)  = entr_num_dis(imod) + acc_normfac_imod**2 &
                   & * de(ix,is) &
                   & * d3X * d3v * real(conjg(fdis) / fmaxwl(ix,i,j,k,is) &
                   & * num_disp01(iih)) * parseval_correction(imod)
                ! numerical parallel velocity dissipation
                entr_num_vp(imod)   = entr_num_vp(imod)  + acc_normfac_imod**2 &
                   & * de(ix,is) &
                   & * d3X * d3v * real(conjg(fdis) / fmaxwl(ix,i,j,k,is) &
                   & * num_disp02(iih)) * parseval_correction(imod)
                ! numerical dissipation from perpendicular (hyper)diffusion
                entr_num_perp(imod) = entr_num_perp(imod) + acc_normfac_imod**2 &
                   & * de(ix,is) &
                   & * d3X * d3v * real(conjg(fdis) / fmaxwl(ix,i,j,k,is) &
                   & * num_disp03(iih)) * parseval_correction(imod)

                ! contribution of the time dependent temperature source
                entr_temp_src(imod) = entr_temp_src(imod) + acc_normfac_imod**2 &
                   & * de(ix,is) &
                   & * d3X * d3v * real(conjg(fdis) / fmaxwl(ix,i,j,k,is) &
                   & * S(iih)) * parseval_correction(imod)
                
                !The entropy in the Landau damping term (VII)
                entr_landau(imod) = entr_landau(imod) + d3X * d3v * &
                   & real(conjg(fdis)*dphids(i)) * signz(is) * vthrat(is)* &
                   & vpgr(i,j,k,is) * ffun(ix,i) / tmp(ix,is) &
                   & * parseval_correction(imod)
               
                ! dissipation from outflow at downstream s boundary
                entr_outflow(imod) = entr_outflow(imod) + acc_normfac_imod**2 &
                   & * de(ix,is) &
                   & * d3X * d3v * real(conjg(fdis) / fmaxwl(ix,i,j,k,is) &
                   & * parseval_correction(imod) &
                   & * outflow(iih))

                dum = &
                   & acc_normfac_imod**2 &
                   & * de(ix,is) * ( pflux_det(k,j,i,ix,imod,is) * fp(ix,is) &
                   & + mas(is)*vthrat(is)**2 * eflux_det(k,j,i,ix,imod,is) * tp(ix,is) &
                   & - 1.5 * pflux_det(k,j,i,ix,imod,is) * tp(ix,is) &
                   & + cfen(1,is)/tmp(ix,is) * pflux_det(k,j,i,ix,imod,is) * vp(ix,is) &
                   & + mas(is) * vthrat(is) * vflux_det(k,j,i,ix,imod,is) * vp(ix,is) &
                   & * 2.0 / tmp(ix,is) &
                   & + mas(is)/tmp(ix,is) * vcor*vcor * pflux_det(k,j,i,ix,imod,is) &
                   & * lfun)
                entr_src_kyx_fsa_local(imod, ix) = &
                   & entr_src_kyx_fsa_local(imod, ix) + dum
                if(gs(i) == xy_slice_ipar) then
                  entr_src_kyx_local(imod, ix) = entr_src_kyx_local(imod, ix) + &
                     dum
                end if

              end do !vpar
            end do !mu
          end do !ns
        end do !nsp
      end do !nx
    end do !nmod

    ! mpi reduce 1D arrays:
    call mpireduce_sum_inplace(ene_e,shape(ene_e))
    call mpireduce_sum_inplace(ene_f,shape(ene_f))
    call mpireduce_sum_inplace(ene_m,shape(ene_m))
    call mpireduce_sum_inplace(ene_e2,shape(ene_e2))
    call mpireduce_sum_inplace(entr,shape(entr))
    call mpireduce_sum_inplace(entr_tr,shape(entr_tr), COMM_SP_EQ)
    call mpireduce_sum_inplace(entr_ps,shape(entr_ps), COMM_SP_EQ)
    call mpireduce_sum_inplace(entr_field,shape(entr_field))
    call mpireduce_sum_inplace(last_entr,shape(last_entr))
    call mpireduce_sum_inplace(last_entr_tr,shape(last_entr_tr), COMM_SP_EQ)
    call mpireduce_sum_inplace(last_entr_ps,shape(last_entr_ps), COMM_SP_EQ)
    call mpireduce_sum_inplace(last_entr_field,shape(last_entr_field))
    if (adiabatic_electrons .and. zonal_adiabatic) then
      ! This mpireduction completes the s-integral in in entr_field and
      ! also the other integrals.  It is only after this reduction, that
      ! entr_field_adiacorr is practically a real number (the imag. part
      ! being small).
      call mpireduce_sum_inplace(entr_field_adiacorr, &
         & shape(entr_field_adiacorr))
      call mpireduce_sum_inplace(last_entr_field_adiacorr, &
         & shape(last_entr_field_adiacorr))
    end if
    call mpireduce_sum_inplace(entr_coll,shape(entr_coll))
    call mpireduce_sum_inplace(entr_num_dis,shape(entr_num_dis))
    call mpireduce_sum_inplace(entr_num_vp,shape(entr_num_vp))
    call mpireduce_sum_inplace(entr_num_perp,shape(entr_num_perp))
    call mpireduce_sum_inplace(entr_temp_src,shape(entr_temp_src))
    call mpireduce_sum_inplace(entr_outflow,shape(entr_outflow))
    call mpireduce_sum_inplace(entr_landau,shape(entr_landau))

    if(proc_subset(0,1,1,1,0)) then
      call mpireduce_sum_inplace(entr_src01,shape(entr_src01),COMM_SP_NE_X_NE)
      call mpireduce_sum_inplace(entr_src02,shape(entr_src02),COMM_SP_NE_X_NE)
      call mpireduce_sum_inplace(entr_src03,shape(entr_src03),COMM_SP_NE_X_NE)
      call mpireduce_sum_inplace(entr_src04,shape(entr_src04),COMM_SP_NE_X_NE)
      call mpireduce_sum_inplace(entr_src05,shape(entr_src05),COMM_SP_NE_X_NE)
      call mpireduce_sum_inplace(entr_src06,shape(entr_src06),COMM_SP_NE_X_NE)
    end if
    if(proc_subset(0,xy_slice_ipar,0,0,0)) then
      call mpireduce_sum_inplace(entr_kyx_local,shape(entr_kyx_local), &
         & COMM_S_EQ_X_EQ)
      call mpireduce_sum_inplace(entr_src_kyx_local,shape(entr_src_kyx_local),&
         & COMM_S_EQ_X_EQ)
    end if
    call mpireduce_sum_inplace(entr_kyx_fsa_local,shape(entr_kyx_fsa_local), &
       & COMM_X_EQ)
    call mpireduce_sum_inplace(entr_src_kyx_fsa_local,shape(entr_src_kyx_fsa_local),&
       & COMM_X_EQ)
    
    if(root_processor) then
      do imod=1, nmod
        if (adiabatic_electrons .and. zonal_adiabatic) then
          ! add the adiabatic-electron correction to entr_field

          !entr_field(imod) = entr_field(imod) + entr_field_adiacorr(imod)
          !last_entr_field(imod) = last_entr_field(imod) + &
          !   & last_entr_field_adiacorr(imod)
          continue
        end if

        inverse_fnorm = 1.0/fnorm1d(imod)
        if (abs(entr(imod) - inverse_fnorm**2 * last_entr(imod)) < r_tiny) then
          dt_entr(imod) = 0.0
        else
          dt_entr(imod) = (entr(imod) - &
             & inverse_fnorm**2 * last_entr(imod))/delta_time_energetics
        end if

        if (abs(entr_field(imod) - &
           & inverse_fnorm**2 * last_entr_field(imod)) < r_tiny) then
          dt_entr_field(imod) = 0.0
        else
          dt_entr_field(imod) = (entr_field(imod) - &
             & inverse_fnorm**2 * last_entr_field(imod))/delta_time_energetics
        end if

        if (abs(entr_field_adiacorr(imod) - &
           & inverse_fnorm**2 * last_entr_field_adiacorr(imod)) < r_tiny) then
          dt_entr_field_adiacorr = 0.0
        else
          dt_entr_field_adiacorr(imod) = (entr_field_adiacorr(imod) - &
             & inverse_fnorm**2 * last_entr_field_adiacorr(imod))/delta_time_energetics
        end if

      end do

      where(entr_outflow < r_tiny)
        entr_outflow = 0
      end where
    end if

    if(proc_subset(1,1,1,1,0)) then
      do imod = 1, nmod
        inverse_fnorm = 1.0/fnorm1d(imod)
        do is = 1, nsp
          if (abs(entr_tr(imod,is) - inverse_fnorm**2 * last_entr_tr(imod,is)) &
             & < r_tiny) then
            dt_entr_tr(imod,is) = 0.0
          else
            dt_entr_tr(imod,is) = (entr_tr(imod,is) - &
               & inverse_fnorm**2 * last_entr_tr(imod,is))/delta_time_energetics
          end if

          if (abs(entr_ps(imod,is) - inverse_fnorm**2 * last_entr_ps(imod,is)) &
             & < r_tiny) then
            dt_entr_ps(imod,is) = 0.0
          else
            dt_entr_ps(imod,is) = (entr_ps(imod,is) - &
               & inverse_fnorm**2 * last_entr_ps(imod,is))/delta_time_energetics
          end if
        end do
      end do
    end if

  end subroutine calc_largestep

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> The energetics subroutine is interested in calculating the time derivative
  !> of the quantities entr and entr_field. Do do this, these have to be
  !> calculated and stored in the timestep that preceeds the one the energetics
  !> subroutine is called in, as well as in the subroutine energetics itself.
  !> This routine here is called from both places.
  !----------------------------------------------------------------------------
  subroutine calc_entr(entr, entr_field, entr_field_adiacorr, &
     & entr_tr, entr_ps, entr_kyx_local, entr_kyx_fsa_local)
    use grid,           only : nmod, nx, ns, nvpar, nmu, nsp, parallel_s
    use grid,           only : proc_subset, gs
    use geom,           only : ints, bn
    use mode,           only : ixzero
    use velocitygrid,   only : intmu, intvp
    use dist,           only : fdisi, fmaxwl, iphi
    use index_function, only : indx
    use matdat,         only : get_f_from_g
    use components,     only : signz, tmp, de, iadia, adiabatic_electrons
    use fields,         only : get_averaged_phi
    use control,        only : zonal_adiabatic, spectral_radius
    use mpicomms,       only : COMM_S_NE
    use mpiinterface,   only : mpiallreduce_sum_inplace
    use rotation,       only : cfen
    use diagnos_generic, only : parseval_correction, xy_slice_ipar
    use diagnos_generic, only : get_tr_ps_mask
    use global,         only : r_tiny

    real, intent(out) :: entr(nmod)
    complex, intent(out) :: entr_field(nmod)
    complex, intent(out) :: entr_field_adiacorr(nmod)
    real, intent(out), optional :: entr_kyx_local(nmod,nx)
    real, intent(out), optional :: entr_kyx_fsa_local(nmod,nx)
    real, intent(out), optional :: entr_tr(nmod,nsp)
    real, intent(out), optional :: entr_ps(nmod,nsp)
    integer :: i, j, k, is, imod, ix
    integer :: mask
    real :: d2X, d3X, dumint, d3v
    complex :: fdis, phi, phi_ga, phi_fsa
    real :: fdis_sq, phi_sq, phi_ga_sq

    entr = 0.
    entr_field = 0.
    entr_field_adiacorr = 0.
    if(present(entr_kyx_local)) entr_kyx_local = 0.
    if(present(entr_kyx_fsa_local)) entr_kyx_fsa_local = 0.
    if(present(entr_tr)) entr_tr = 0.
    if(present(entr_ps)) entr_ps = 0.

    do imod=1,nmod
      d2X = 1.0
      do ix=1,nx

        phi_fsa = 0.0
        do i=1,ns
          phi  = fdisi(indx(iphi,imod,ix,i))
          ! write(*,*) 'fdisi', phi
          phi_fsa = phi_fsa + phi * ints(i)
        end do
        if(parallel_s) then
          ! An mpi-reduction is necessary to complete the integral over
          ! the s direction of the potential.
          call mpiallreduce_sum_inplace(phi_fsa, 1, COMM=COMM_S_NE)
        end if

        ! complex conjugate because of Parsefal's theorem
        phi_fsa=conjg(phi_fsa)

        ! sum over species must be inside s-integral
        do i=1,ns
          phi  = fdisi(indx(iphi,imod,ix,i))
          phi_sq = abs(phi)**2
          d3X = ints(i) * d2X
          do is=1,nsp
            if(de(ix,is) < r_tiny) cycle
            do j=1,nmu
              phi_ga  = get_averaged_phi(imod,ix,i,j,is,fdisi)
              phi_ga_sq = abs(phi_ga)**2
              dumint = intmu(j)*bn(ix,i)
              do k=1,nvpar
                fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)
                fdis_sq = abs(fdis)**2
                d3v = dumint * intvp(i,j,k,is)

                ! Important: entr is the entropy S.
                ! The negative of it is often called the fluctuation
                ! intensity H = -S

                ! entr_field is the "field entropy" (-W).
                entr(imod) = entr(imod) - parseval_correction(imod) * 0.5 &
                   & * de(ix,is) &
                   & * fdis_sq / fmaxwl(ix,i,j,k,is) * d3v * d3X

                if(present(entr_tr)) then
                  mask = get_tr_ps_mask(ix,i,j,k,is)
                  entr_tr(imod,is) = entr_tr(imod,is) - parseval_correction(imod) * 0.5 &
                     & * de(ix,is) &
                     & * fdis_sq / fmaxwl(ix,i,j,k,is) * d3v * d3X * mask
                  if(present(entr_ps)) then
                    entr_ps(imod,is) = entr_ps(imod,is) - parseval_correction(imod) * 0.5 &
                       & * de(ix,is) &
                       & * fdis_sq / fmaxwl(ix,i,j,k,is) * d3v * d3X * (-mask+1)
                  end if
                end if

                if(present(entr_kyx_local) .and. gs(i)==xy_slice_ipar) then
                  entr_kyx_local(imod,ix) = entr_kyx_local(imod,ix) &
                     & - parseval_correction(imod) * 0.5 &
                     & * de(ix,is) &
                     & * fdis_sq / fmaxwl(ix,i,j,k,is) * d3v * d3X
                end if
                if(present(entr_kyx_fsa_local)) then
                  entr_kyx_fsa_local(imod,ix) = entr_kyx_fsa_local(imod,ix) &
                     & - parseval_correction(imod) * 0.5 &
                     & * de(ix,is) &
                     & * fdis_sq / fmaxwl(ix,i,j,k,is) * d3v * d3X
                end if

                entr_field(imod) = entr_field(imod) &
                   & - parseval_correction(imod) * 0.5 * signz(is)**2 &
                   & * de(ix,is) &
                   & * fmaxwl(ix,i,j,k,is)*(phi_sq - phi_ga_sq)/tmp(ix,is)**2 &
                   & * d3v * d3X

              end do !vpar
            end do !mu
          end do !nsp
          
          ! use proc_subset to compute the following only once, not
          ! n_procs_sp times.
          if (proc_subset(0,0,1,1,1) .and. adiabatic_electrons .and. zonal_adiabatic) then
            ! The zonal flow correction term affects the field-entropy
            ! because there is a special term for adiabatic electrons in
            ! the poisson equation.

            ! Because of the presence of phi, this expression is
            ! complex-valued at this point. The calculation is done with a
            ! dummy variable.

            ! iadia is 1 for adiabatic electrons, and 0 otherwise, but
            ! in the nonadiabatic case this block is not computed anyway.
            entr_field_adiacorr(imod) = entr_field_adiacorr(imod) &
               & - parseval_correction(imod) * 0.5 * de(ix, nsp+iadia) &
               & * (phi_sq - phi*phi_fsa)/(tmp(ix,nsp+iadia)**2) &
               & * d3X

            ! In a rotating frame the electron contribution needs a further
            ! correction, due to the presence of the centrifugal potential
            ! in the electron term in the poisson equation.
            if(imod == 1) then
              if((spectral_radius .and. ix == ixzero) &
                 & .or. .not. spectral_radius) then
                entr_field_adiacorr(imod) = entr_field_adiacorr(imod) &
                   & - de(ix, nsp+iadia) &
                   & * phi * cfen(i,nsp+iadia) &
                   & / (tmp(ix,nsp+iadia)**2) &
                   & * d3X
              end if
            end if

          end if
        end do !ns
      end do !nx

      ! After the final mpireduction (this is done in the subroutine energetics)
      ! of entr_field, phi*phi_fsa becomes phi_fsa**2, and only then
      ! the complex variable entr_field should be pretty much real, the imag part being
      ! very tiny.
      ! (Maybe this is true as long as the rotational adiabatic-electrons correction
      ! is commented out).
      ! If parallel_s is true, the imag part should be considerable here:
      ! write (*,*) entr_field(imod)
    end do !nmod
  end subroutine calc_entr

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()
    use io, only : append_chunk, xy_fmt, xy_fmt_long, ascii_fmt
    use mpiinterface, only : root_processor
    use control, only : nlapar
    use grid, only : proc_subset, ns
    use diagnos_generic, only : kyx_output_array, xy_slice_ipar

    ! calculation of the total energy 
    if (lcalc_tot_energy) then 
      call calc_tot_energy
      ! Write data
      if (root_processor) then
        call append_chunk(i_tot_energy, tot_energy, xy_fmt, ascii_fmt)
      end if
    endif

    if (.not.lcalc_energetics) return

    ! calculation of the entropy balance terms
    if (lcalc_energetics) then
      call calc_largestep
    end if

    ! gather and output species-dependent data
    if(proc_subset(1,1,1,1,0)) then
      call kysp_output_array(dt_entr_tr, i_dt_entr_tr)
      call kysp_output_array(dt_entr_ps, i_dt_entr_ps)
      call kysp_output_array(entr_tr, i_entr_tr)
      call kysp_output_array(entr_ps, i_entr_ps)
    end if

    ! Write data
    if (root_processor) then

      call append_chunk(i_ene_e, ene_e, xy_fmt, ascii_fmt)
      call append_chunk(i_ene_f, ene_f, xy_fmt, ascii_fmt)
      if (nlapar) then
        call append_chunk(i_ene_m, ene_m, xy_fmt, ascii_fmt)
      end if
      call append_chunk(i_ene_e2, ene_e2, xy_fmt, ascii_fmt)
      call append_chunk(i_dt_entr, dt_entr, xy_fmt, ascii_fmt)
      call append_chunk(i_dt_entr_field, real(dt_entr_field), xy_fmt, ascii_fmt)
      call append_chunk(i_dt_entr_field_adiacorr, real(dt_entr_field_adiacorr), &
         & xy_fmt, ascii_fmt)
      call append_chunk(i_entr, entr, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_field, real(entr_field), xy_fmt, ascii_fmt)
      call append_chunk(i_entr_field_adiacorr, real(entr_field_adiacorr), &
         & xy_fmt, ascii_fmt)
      call append_chunk(i_entr_coll, entr_coll, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_num_dis, entr_num_dis, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_num_vp, entr_num_vp, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_num_perp, entr_num_perp, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_temp_src, entr_temp_src, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_outflow, entr_outflow, xy_fmt_long, ascii_fmt)

      call append_chunk(i_entr_src01, entr_src01, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_src02, entr_src02, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_src03, entr_src03, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_src04, entr_src04, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_src05, entr_src05, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_src06, entr_src06, xy_fmt, ascii_fmt)
      call append_chunk(i_entr_landau, entr_landau, xy_fmt, ascii_fmt)

    end if
    ! output entropy spectrum, including x gather
    call kyx_output_array(entr_kyx_local, i_entr_kyx, &
       & proc_subset(0,xy_slice_ipar,1,1,1), xy_slice_ipar <= ns)
    call kyx_output_array(entr_src_kyx_local, i_entr_src_kyx, &
       & proc_subset(0,xy_slice_ipar,1,1,1), xy_slice_ipar <= ns)
    
    call kyx_output_array(entr_kyx_fsa_local, i_entr_kyx_fsa, &
       & proc_subset(0,1,1,1,1), .true.)
    call kyx_output_array(entr_src_kyx_fsa_local, i_entr_src_kyx_fsa, &
       & proc_subset(0,1,1,1,1), .true.)

  end subroutine output

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine kysp_output_array(array, lun)
    use grid, only : nmod, nsp, number_of_species
    use mpiinterface, only : root_processor, gather_array
    use io, only : xy_fmt, ascii_fmt, append_chunk
    use mpicomms, only : COMM_DUMMY, COMM_SP_NE
    real, dimension(nmod,nsp), intent(in) :: array
    integer, intent(in) :: lun

    call gather_array(kysp_global_buf,nmod,number_of_species, &
       & array,nmod,nsp, &
       & COMM_DUMMY,COMM_SP_NE,ALLGATHER=.false.)

    if(root_processor) then
      call append_chunk(lun, kysp_global_buf, xy_fmt, ascii_fmt)
    end if


  end subroutine kysp_output_array


end module diagnos_energetics
