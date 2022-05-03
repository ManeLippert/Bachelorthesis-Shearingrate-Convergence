!------------------------------------------------------------------------------
!> This module coordinates all the diagnostics, the sub-modules diagnos_*
!> contain individual diagnostics, grouped mainly by related physical quantities.
!> All calculations that do not affect the solution, and data gathering required  
!> for output should appear in the diagnos_* modules.
!> Data output is handled by the io_* modules called from diagnos_* modules
!-------------------------------------------------------------------------------
module diagnostic
 
  implicit none

  private

  public :: diagnostic_read_nml, diagnostic_write_nml, diagnostic_read_last_data
  public :: diagnostic_bcast_nml, diagnostic_check_params
  public :: diagnostic_allocate, diagnostic_initial_output
  public :: diagnostic_final_output, diagnostic_naverage
  public :: diagnostic_init,diagnostic_finalize, diagnostic_pre_naverage
  public :: abort_if_old_data_exists

  interface diagnostic_write_nml
    module procedure diagnostic_read_nml
  end interface

  !> flush the fluxes and time etc. every nflush_ts timesteps (default is 10).
  integer, save :: nflush_ts = 0

  !> output switches which are only used in this general diagnostic module
  logical, save, public :: screen_output,lfinal_output
  
  !> file count and suffix
  integer, save :: file_count = 0

  !> a matrix of flags which is modified by init() in the particular
  !> diagnostics
  logical, save, allocatable :: requirements(:,:)

  !--------- DEPRECATED SWITCHES: -----------
  logical, save, public :: lcalc_freq
  logical, save, public :: lpflux,leflux,lvflux
  logical, save, public :: xy_bin
  logical, save, public :: lparallel_phi, lparallel_apar, lparallel_bpar

contains 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> read or write the diagnostic nml if ldiagnostic_namelist is set in control
!----------------------------------------------------------------------------
subroutine diagnostic_read_nml(ilun,io_stat,lwrite)
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  use diagnos_generic, only : set_defaults_generic => set_default_nml_values, &
     & lwrite_output1, lphi_diagnostics, xy_estep, &
     & lisl_follow, zonal_scale_3d, &
     & spc3d, phi3d, den3d, apa3d, apc3d, ene3d, bpa3d, bpc3d, &
     & xy_fluxes, xy_fluxes_bi, xy_fluxes_p, xy_fluxes_v, xy_fluxes_fsa, &
     & xy_fluxes_em, xy_fluxes_k, &
     & xy_fluxes_bpar, xy_slice_ipar, &
     & xy_dens, xy_temp, xy_current, xy_current2, &
     & xy_spec, out3d_interval, &
     & lradial_profile 
  use diagnos_fluxes, only : set_defaults_fluxes => set_default_nml_values, &
     & lcalc_fluxes, flux3d, &
     & lfluxes_spectra,lfluxes_em_spectra
  use diagnos_kinenergy_trappas, only : set_defaults_kinenergy_trappas => &
     & set_default_nml_values, lcalc_kinenergy_trappas
  use diagnos_fluxes_vspace, only : &
     & set_defaults_fluxes_vspace => set_default_nml_values, &
     & lfluxes_detail, lfluxes_vspace, &
     & lfluxes_vspace_phi,lfluxes_vspace_em,lfluxes_vspace_bpar
  use diagnos_energetics, only : &
     & set_defaults_energetics => set_default_nml_values, &
     & lcalc_energetics, lcalc_tot_energy
  use diagnos_corr, only : &
     & set_defaults_corr => set_default_nml_values, &
     & lcalc_corr, imod_corr
  use diagnos_fields, only : set_defaults_fields => set_default_nml_values, &
     & xy_phi,xy_apar, xy_bpar, xy_vort, &
     & xs_phi, nmodepoints, imod_field, field_fsa_kyzero, &
     & kykxs_phi, kykxs_apar, kykxs_bpar
  use diagnos_stresses, only : set_defaults_stresses => set_default_nml_values, &
        & lcalc_stresses
  use diagnos_eng, only : set_defaults_eng => set_default_nml_values, &
     & lmode_energy
  use diagnos_growth_freq, only : &
     & set_defaults_growth_freq => set_default_nml_values, &
     & lamplitudes,lgrowth_rates,lfrequencies
  use diagnos_mode_struct, only : &
     & set_defaults_mode_struct => set_default_nml_values, &
     & lparallel_output, parallel_output_timestamps, &
     & lrotate_parallel
  use diagnos_velspace, only : set_defaults_velspace => set_default_nml_values, &
     & lvelspace_output, psi_velspace,zeta_velspace, npointsvel,lfinvel
  use diagnos_jdote, only : set_defaults_jdote => set_default_nml_values, &
     & lcalc_jdote, lcalc_jdote_fs
  use diagnos_rad, only : set_defaults_rad => set_default_nml_values, &
     & lrad_moment, lrad_entropy, lradial_entropy, lrad_field,  lrad_tint, lrad_kpar
  use diagnos_f, only : set_defaults_f => set_default_nml_values, &
     & tavg_start, tavg_end 
  use diagnos_grid, only : set_defaults_grid => set_default_nml_values
  use diagnos_moments, only : set_defaults_moments => set_default_nml_values, &
     & kykxs_moments, kykxs_j0_moments,  kykxs_j1_moments , &
     & xs_kyzero_dens, xs_kyzero_current, xs_kyzero_current2, &
     & xs_kyzero_ene, xs_kyzero_ene_perp, xs_kyzero_ene_par, &
     & xs_kyzero_phi_ga_fm, xs_kyzero_phi_ga_deltaf
  use diagnos_zfshear, only : set_defaults_zfshear => set_default_nml_values, &
     & lcalc_zfshear
  use diagnos_nonlin_transfer, only : &
     & set_defaults_nonlin_transfer => set_default_nml_values, &
     & lnonlin_transfer, lnonlin_transfer_fsa, nonlin_transfer_interval
  use diagnos_zonal_evo, only : &
     & set_defaults_zonal_evo => set_default_nml_values, &
     & lcalc_zonal_evo, zevo_detail, zevo_xy_spec
  use diagnos_timetrace, only : set_defaults_timetrace => &
     & set_default_nml_values, lcomplex_timetrace
  use diagnos_cross_phase, only : set_defaults_cross_phase => &
     & set_default_nml_values, cross_phases, cross_phase_timetrace
  use diagnos_neoequil, only : set_defaults_neoequil => &
     & set_default_nml_values
  use diagnos_matrix, only : set_defaults_matrix => &
     & set_default_nml_values, loutput_matrix
  use io, only : lmpi_broken_io
  ! deprecated:
  use diagnos_generic, only : lphi_diagnostics



  integer, intent(in) :: ilun
  integer, intent(out) :: io_stat
  logical, optional, intent(in) :: lwrite

  namelist /diagnostic/ &
     ! from module diagnostic:
     & nflush_ts, screen_output, lfinal_output, &
     ! from module diagnos_generic:
     & lwrite_output1, &
     & lphi_diagnostics, &
     & xy_estep, out3d_interval, &
     & lisl_follow, zonal_scale_3d, &
     & spc3d, phi3d, den3d, apa3d, apc3d, ene3d, bpa3d, bpc3d, &
     & xy_fluxes, xy_fluxes_bi, xy_fluxes_p, xy_fluxes_v, xy_fluxes_fsa, &
     & xy_fluxes_em, xy_fluxes_k, &
     & xy_fluxes_bpar, xy_slice_ipar, &
     & xy_dens, xy_temp, xy_current, xy_current2, &
     & xy_spec, &
     & xs_kyzero_dens, xs_kyzero_current, xs_kyzero_current2, &
     & xs_kyzero_ene, xs_kyzero_ene_perp, xs_kyzero_ene_par, &
     & xs_kyzero_phi_ga_fm, xs_kyzero_phi_ga_deltaf, &
     & lradial_profile, &
     ! from module diagnos_fluxes:
     & lcalc_fluxes, flux3d, &
     & lfluxes_spectra,lfluxes_em_spectra, &
     ! from module diagnos_kinenergy_trappas:
     & lcalc_kinenergy_trappas, &
     ! from module diagnos_fluxes_vspace:
     & lfluxes_detail, lfluxes_vspace, &
     & lfluxes_vspace_phi,lfluxes_vspace_em,lfluxes_vspace_bpar, &
     ! from module diagnos_energetics:
     & lcalc_energetics, lcalc_tot_energy, &
     ! from module diagnos_corr:
     & lcalc_corr, imod_corr , & 
     ! from module diagnos_fields:
     & lparallel_phi, lparallel_apar, lparallel_bpar,&
     & xy_phi,xy_apar, xy_bpar, xy_vort, &
     & xs_phi, nmodepoints, imod_field, field_fsa_kyzero, &
     & kykxs_phi, kykxs_apar, kykxs_bpar, &
     ! from module diagnos_moments:
     & kykxs_moments, kykxs_j0_moments, kykxs_j1_moments, &
     ! from module diagnos_f:
     & tavg_start, tavg_end, &
     ! from module diagnos_grid:
     ! from module diagnos_stresses:
     & lcalc_stresses, &
     ! from module diagnos_eng:
     & lmode_energy, & 
     ! from module diagnos_growth_freq:
     & lamplitudes,lgrowth_rates,lfrequencies, &
     ! from module diagnos_mode_struct:
     & lparallel_output, parallel_output_timestamps, &
     & lrotate_parallel, &
     ! from module diagnos_velspace:
     & lvelspace_output, psi_velspace,zeta_velspace, npointsvel,lfinvel, &
     ! from module diagnos_jdote:
     & lcalc_jdote, lcalc_jdote_fs, &
     ! from module diagnos_rad:
     & lrad_moment, lrad_entropy, lradial_entropy, lrad_field,  lrad_tint, & 
     & lrad_kpar, &
     ! from module diagnos_zfshear:
     & lcalc_zfshear, &
     ! from module diagnos_nonlin_transfer:
     & lnonlin_transfer, lnonlin_transfer_fsa, nonlin_transfer_interval, &
     ! from module diagnos_zonal_evo:
     & lcalc_zonal_evo, zevo_detail, zevo_xy_spec, &
     ! from module diagnos_cross_phase:
     & cross_phases, cross_phase_timetrace, &
     ! from module diagnos_matrix:
     & loutput_matrix, &
     ! from module diagnos_timetrace:
     & lcomplex_timetrace, &
     ! from module io:
     & lmpi_broken_io, &
     ! obsolete or deprecated:
     & xy_bin, &
     & lpflux,leflux,lvflux, &
     & lcalc_freq

  if (present(lwrite)) then
    if (.not. lwrite) then
      ! set any defaults to on/true/non-zero etc (all false at declaration)

      lfinal_output = .true.
      screen_output = .true.
      lmpi_broken_io = .true.
      nflush_ts = 10

      call set_defaults_generic()
      call set_defaults_energetics()
      call set_defaults_corr()
      call set_defaults_growth_freq()
      call set_defaults_eng()
      call set_defaults_f()
      call set_defaults_fluxes_vspace()
      call set_defaults_grid()
      call set_defaults_mode_struct()
      call set_defaults_rad()
      call set_defaults_stresses()
      call set_defaults_fields()
      call set_defaults_moments()
      call set_defaults_fluxes()
      call set_defaults_kinenergy_trappas
      call set_defaults_velspace()
      call set_defaults_jdote()
      call set_defaults_zfshear()
      call set_defaults_nonlin_transfer()
      call set_defaults_zonal_evo()
      call set_defaults_timetrace()
      call set_defaults_cross_phase()

      call set_defaults_neoequil()
      call set_defaults_matrix()

      ! set default values of deprecated variables and switches here
      !...
  
      ! read nml
      io_stat = 0
      read(ilun,NML=diagnostic,IOSTAT=io_stat)
    else
      ! do nothing
    end if
  else
    if(root_processor) write(ilun,NML=diagnostic)

    ! from module diagnostic:
    call write_run_parameter('diagnostic', 'nflush_ts', nflush_ts)
    call write_run_parameter('diagnostic', 'screen_output', screen_output)
    call write_run_parameter('diagnostic', 'lfinal_output', lfinal_output) 
    ! from module diagnos_generic:
    call write_run_parameter('diagnostic', 'lwrite_output1', lwrite_output1) 
    call write_run_parameter('diagnostic', 'lphi_diagnostics', lphi_diagnostics) 
    call write_run_parameter('diagnostic', 'xy_estep', xy_estep) 
    call write_run_parameter('diagnostic', 'out3d_interval', out3d_interval)
    call write_run_parameter('diagnostic', 'lisl_follow', lisl_follow)
    call write_run_parameter('diagnostic', 'zonal_scale_3d', zonal_scale_3d) 
    call write_run_parameter('diagnostic', 'spc3d', spc3d)
    call write_run_parameter('diagnostic', 'phi3d', phi3d)
    call write_run_parameter('diagnostic', 'den3d', den3d)
    call write_run_parameter('diagnostic', 'apa3d', apa3d)
    call write_run_parameter('diagnostic', 'apc3d', apc3d)
    call write_run_parameter('diagnostic', 'ene3d', ene3d)
    call write_run_parameter('diagnostic', 'bpa3d', bpa3d)
    call write_run_parameter('diagnostic', 'bpc3d', bpc3d) 
    call write_run_parameter('diagnostic', 'xy_fluxes', xy_fluxes)
    call write_run_parameter('diagnostic', 'xy_fluxes_bi', xy_fluxes_bi)
    call write_run_parameter('diagnostic', 'xy_fluxes_p', xy_fluxes_p)
    call write_run_parameter('diagnostic', 'xy_fluxes_v', xy_fluxes_v)
    call write_run_parameter('diagnostic', 'xy_fluxes_fsa', xy_fluxes_fsa) 
    call write_run_parameter('diagnostic', 'xy_fluxes_em', xy_fluxes_em)
    call write_run_parameter('diagnostic', 'xy_fluxes_k', xy_fluxes_k) 
    call write_run_parameter('diagnostic', 'xy_fluxes_bpar', xy_fluxes_bpar)
    call write_run_parameter('diagnostic', 'xy_slice_ipar', xy_slice_ipar) 
    call write_run_parameter('diagnostic', 'xy_dens', xy_dens)
    call write_run_parameter('diagnostic', 'xy_temp', xy_temp)
    call write_run_parameter('diagnostic', 'xy_current', xy_current)
    call write_run_parameter('diagnostic', 'xy_current2', xy_current2) 
    call write_run_parameter('diagnostic', 'xs_kyzero_dens', xs_kyzero_dens)
    call write_run_parameter('diagnostic', 'xs_kyzero_current', xs_kyzero_current)
    call write_run_parameter('diagnostic', 'xs_kyzero_current2', xs_kyzero_current2)
    call write_run_parameter('diagnostic', 'xs_kyzero_ene', xs_kyzero_ene)
    call write_run_parameter('diagnostic', 'xs_kyzero_ene_perp', xs_kyzero_ene_perp)
    call write_run_parameter('diagnostic', 'xs_kyzero_ene_par', xs_kyzero_ene_par)
    call write_run_parameter('diagnostic', 'xs_kyzero_phi_ga_fm', xs_kyzero_phi_ga_fm)
    call write_run_parameter('diagnostic', 'xs_kyzero_phi_ga_deltaf', xs_kyzero_phi_ga_deltaf)
    call write_run_parameter('diagnostic', 'xy_spec', xy_spec) 
    call write_run_parameter('diagnostic', 'lradial_profile', lradial_profile)
    ! from module diagnos_fluxes:
    call write_run_parameter('diagnostic', 'lcalc_fluxes', lcalc_fluxes)
    call write_run_parameter('diagnostic', 'flux3d', flux3d) 
    call write_run_parameter('diagnostic', 'lfluxes_spectra', lfluxes_spectra)
    call write_run_parameter('diagnostic', 'lfluxes_em_spectra', &
       & lfluxes_em_spectra) 
   ! from module diagnos_kinenergy_trappas:
    call write_run_parameter('diagnostic', 'lcalc_kinenergy_trappas', lcalc_kinenergy_trappas) 
    ! from module diagnos_fluxes_vspace:
    call write_run_parameter('diagnostic', 'lfluxes_detail', lfluxes_detail)
    call write_run_parameter('diagnostic', 'lfluxes_vspace', lfluxes_vspace) 
    call write_run_parameter('diagnostic', 'lfluxes_vspace_phi', &
       & lfluxes_vspace_phi)
    call write_run_parameter('diagnostic', 'lfluxes_vspace_em', &
       & lfluxes_vspace_em)
    call write_run_parameter('diagnostic', 'lfluxes_vspace_bpar', &
       & lfluxes_vspace_bpar) 
    ! from module diagnos_energetics:
    call write_run_parameter('diagnostic', 'lcalc_energetics', &
       & lcalc_energetics)
    call write_run_parameter('diagnostic', 'lcalc_tot_energy', &
       & lcalc_tot_energy) 
    ! from module diagnos_corr:
    call write_run_parameter('diagnostic', 'lcalc_corr', &
       & lcalc_corr) 
    call write_run_parameter('diagnostic', 'imod_corr', &
       & imod_corr)   
    ! from module diagnos_fields:
    call write_run_parameter('diagnostic', 'xy_phi', xy_phi)
    call write_run_parameter('diagnostic', 'xy_apar', xy_apar)
    call write_run_parameter('diagnostic', 'xy_bpar', xy_bpar)
    call write_run_parameter('diagnostic', 'xy_vort', xy_vort) 
    call write_run_parameter('diagnostic', 'xs_phi', xs_phi)
    call write_run_parameter('diagnostic', 'nmodepoints', nmodepoints)
    call write_run_parameter('diagnostic', 'imod_field', imod_field)
    call write_run_parameter('diagnostic', 'field_fsa_kyzero', field_fsa_kyzero)
    call write_run_parameter('diagnostic', 'kykxs_phi', kykxs_phi)
    call write_run_parameter('diagnostic', 'kykxs_apar', kykxs_apar)
    call write_run_parameter('diagnostic', 'kykxs_bpar', kykxs_bpar)
    ! from module diagnos_moments:
    call write_run_parameter('diagnostic', 'kykxs_moments', kykxs_moments)
    call write_run_parameter('diagnostic', 'kykxs_j0_moments', kykxs_j0_moments)
    call write_run_parameter('diagnostic', 'kykxs_j1_moments', kykxs_j1_moments)
    ! from module diagnos_f:
    ! from module diagnos_grid:
    ! from module diagnos_stresses:
    call write_run_parameter('diagnostic', 'lcalc_stresses', lcalc_stresses) 
    ! from module diagnos_eng:
    call write_run_parameter('diagnostic', 'lmode_energy', lmode_energy) 
    ! from module diagnos_growth_freq:
    call write_run_parameter('diagnostic', 'lamplitudes', lamplitudes)
    call write_run_parameter('diagnostic', 'lgrowth_rates', lgrowth_rates)
    call write_run_parameter('diagnostic', 'lfrequencies', lfrequencies) 
    ! from module diagnos_mode_struct:
    call write_run_parameter('diagnostic', 'lparallel_output', &
       & lparallel_output)
    call write_run_parameter('diagnostic', 'parallel_output_timestamps', &
       & parallel_output_timestamps) 
    call write_run_parameter('diagnostic', 'lrotate_parallel', &
       & lrotate_parallel) 
    ! from module diagnos_velspace:
    call write_run_parameter('diagnostic', 'lvelspace_output', &
       & lvelspace_output)
    call write_run_parameter('diagnostic', 'psi_velspace', psi_velspace)
    call write_run_parameter('diagnostic', 'zeta_velspace', zeta_velspace)
    call write_run_parameter('diagnostic', 'npointsvel', npointsvel)
    call write_run_parameter('diagnostic', 'lfinvel', lfinvel) 
    ! from module diagnos_jdote:
    call write_run_parameter('diagnostic', 'lcalc_jdote', lcalc_jdote)
    call write_run_parameter('diagnostic', 'lcalc_jdote_fs', lcalc_jdote_fs)
    ! from module diagnos_rad:
    call write_run_parameter('diagnostic', 'lrad_moment', lrad_moment)
    call write_run_parameter('diagnostic', 'lrad_entropy', lrad_entropy)
    call write_run_parameter('diagnostic', 'lradial_entropy', lradial_entropy)
    call write_run_parameter('diagnostic', 'lrad_field', lrad_field)
    call write_run_parameter('diagnostic', 'lrad_tint', lrad_tint) 
    call write_run_parameter('diagnostic', 'lrad_kpar', lrad_kpar)
    ! from module diagnos_zfshear:
    call write_run_parameter('diagnostic', 'lcalc_zfshear', lcalc_zfshear) 
    ! from module diagnos_nonlin_transfer:
    call write_run_parameter('diagnostic', 'lnonlin_transfer', &
       & lnonlin_transfer)
    call write_run_parameter('diagnostic', 'lnonlin_transfer_fsa', &
       & lnonlin_transfer_fsa)
    call write_run_parameter('diagnostic', 'nonlin_transfer_interval', &
       & nonlin_transfer_interval)
    ! from module diagnos_zonal_evo
    call write_run_parameter('diagnostic', 'lcalc_zonal_evo', lcalc_zonal_evo)
    call write_run_parameter('diagnostic', 'zevo_detail', zevo_detail)
    call write_run_parameter('diagnostic', 'zevo_xy_spec', zevo_xy_spec)
    ! from module diagnos_cross_phase:
    ! FIXME string array is not implemented:
    ! call write_run_parameter('diagnostic', 'cross_phases', &
    !    & cross_phases)
    call write_run_parameter('diagnostic', 'cross_phase_timetrace', &
       & cross_phase_timetrace)
    ! from module io:
    call write_run_parameter('diagnostic', 'lmpi_broken_io', lmpi_broken_io) 
    ! obsolete or deprecated:
    ! ..(do not output, cause only confusion)

  end if

end subroutine diagnostic_read_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Broadcast the diagnostic namelist: every item should be communicated to
!> all processes.
!----------------------------------------------------------------------------
subroutine diagnostic_bcast_nml

  use mpiinterface, only : mpibcast
  use diagnos_generic, only : bcast_generic => bcast
  use diagnos_grid, only : bcast_grid => bcast
  use diagnos_growth_freq, only : bcast_growth_freq => bcast
  use diagnos_f, only : bcast_f => bcast
  use diagnos_moments, only : bcast_moments => bcast
  use diagnos_mode_struct, only : bcast_mode_struct => bcast
  use diagnos_fluxes, only : bcast_fluxes => bcast
  use diagnos_kinenergy_trappas, only : bcast_kinenergy_trappas => bcast
  use diagnos_fluxes_vspace, only : bcast_fluxes_vspace => bcast
  use diagnos_fields, only : bcast_fields => bcast
  use diagnos_velspace, only : bcast_velspace => bcast
  use diagnos_jdote, only : bcast_jdote => bcast
  use diagnos_energetics, only : bcast_energetics => bcast
  use diagnos_corr, only : bcast_corr => bcast
  use diagnos_nonlin_transfer, only : bcast_nonlin_transfer => bcast
  use diagnos_zonal_evo, only : bcast_zonal_evo => bcast
  use diagnos_eng, only : bcast_eng => bcast
  use diagnos_rad, only : bcast_rad => bcast
  use diagnos_stresses, only : bcast_stresses => bcast
  use diagnos_zfshear, only : bcast_zfshear => bcast
  use diagnos_timetrace, only : bcast_timetrace => bcast
  use diagnos_cross_phase, only : bcast_cross_phase => bcast
  
  use diagnos_neoequil, only : bcast_neoequil => bcast
  use diagnos_matrix, only : bcast_matrix => bcast

  use io, only : lmpi_broken_io

  ! broadcast items of the namelist which are defined in this module here:
  call mpibcast(nflush_ts,1)
  call mpibcast(lfinal_output,1)
  call mpibcast(screen_output,1)
  
  ! (historic) items of the diagnostic namelist which are defined in
  ! other modules but not broadcasted by those:
  call mpibcast(lmpi_broken_io,1)

  call bcast_generic
  call bcast_energetics
  call bcast_corr
  call bcast_nonlin_transfer
  call bcast_eng
  call bcast_f
  call bcast_fields
  call bcast_fluxes
  call bcast_kinenergy_trappas
  call bcast_fluxes_vspace
  call bcast_grid
  call bcast_growth_freq
  call bcast_mode_struct
  call bcast_moments
  call bcast_rad
  call bcast_stresses
  call bcast_velspace
  call bcast_jdote
  call bcast_zfshear
  call bcast_zonal_evo
  call bcast_timetrace
  call bcast_cross_phase
  call bcast_neoequil
  call bcast_matrix

end subroutine diagnostic_bcast_nml


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Check the diagnostic parameters and calculate things needed before
!> allocation of memory will be done.
!>
!----------------------------------------------------------------------------
subroutine diagnostic_check_params

  use diagnos_energetics, only : check_energetics => check
  use diagnos_corr, only : check_corr => check
  use diagnos_eng, only : check_eng => check
  use diagnos_f, only : check_f => check
  use diagnos_fields, only : check_fields => check
  use diagnos_fluxes, only : check_fluxes => check
  use diagnos_kinenergy_trappas, only : check_kinenergy_trappas => check
  use diagnos_fluxes_vspace, only : check_fluxes_vspace => check
  use diagnos_generic, only : check_generic => check
  use diagnos_grid, only : check_grid => check
  use diagnos_growth_freq, only : check_growth_freq => check
  use diagnos_mode_struct, only : check_mode_struct => check
  use diagnos_moments, only : check_moments => check
  use diagnos_nonlin_transfer, only : check_nonlin_transfer => check
  use diagnos_rad, only : check_rad => check
  use diagnos_stresses, only : check_stresses => check
  use diagnos_velspace, only : check_velspace => check
  use diagnos_jdote, only : check_jdote => check
  use diagnos_zfshear, only : check_zfshear => check
  use diagnos_zonal_evo, only : check_zonal_evo => check
  use diagnos_timetrace, only : check_timetrace => check
  use diagnos_cross_phase, only : check_cross_phase => check
  use diagnos_neoequil, only : check_neoequil => check
  use diagnos_matrix, only : check_matrix => check

  call check_generic()

  call check_energetics()
  call check_corr()
  call check_eng()
  call check_f()
  call check_fields()
  call check_fluxes()
  call check_kinenergy_trappas()
  call check_fluxes_vspace()
  call check_grid()
  call check_growth_freq()
  call check_mode_struct()
  call check_moments()
  call check_nonlin_transfer()
  call check_rad()
  call check_stresses()
  call check_velspace()
  call check_jdote()
  call check_zfshear()
  call check_zonal_evo()
  call check_timetrace()
  call check_cross_phase()
  call check_neoequil()
  call check_matrix()
  
end subroutine diagnostic_check_params


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Allocates the arrays of the diagnostic module that are used at runtime;
!> other diagnostics at the end of the run can manage their own memory after
!> most of the runtime memory has been deallocated.
!> Must be called after nonlinear allocate, and kgrid, and is therefore
!> called later than allocate_everything routine.
!----------------------------------------------------------------------------
subroutine diagnostic_allocate
  use diagnos_energetics, only : allocate_energetics => allocate_mem
  use diagnos_corr, only : allocate_corr => allocate_mem
  use diagnos_eng, only : allocate_eng => allocate_mem
  use diagnos_f, only : allocate_f => allocate_mem
  use diagnos_fields, only : allocate_fields => allocate_mem
  use diagnos_fluxes, only : allocate_fluxes => allocate_mem
  use diagnos_kinenergy_trappas, only : allocate_kinenergy_trappas => allocate_mem
  use diagnos_fluxes_vspace, only : allocate_fluxes_vspace => allocate_mem
  use diagnos_generic, only : allocate_generic => allocate_mem
  use diagnos_grid, only : allocate_grid => allocate_mem
  use diagnos_growth_freq, only : allocate_growth_freq => allocate_mem
  use diagnos_mode_struct, only : allocate_mode_struct => allocate_mem
  use diagnos_moments, only : allocate_moments => allocate_mem
  use diagnos_rad, only : allocate_rad => allocate_mem
  use diagnos_stresses, only : allocate_stresses => allocate_mem
  use diagnos_velspace, only : allocate_velspace => allocate_mem
  use diagnos_jdote, only : allocate_jdote => allocate_mem
  use diagnos_zfshear, only : allocate_zfshear => allocate_mem
  use diagnos_zonal_evo, only : allocate_zonal_evo => allocate_mem
  use diagnos_nonlin_transfer, only : allocate_nonlin_transfer => allocate_mem
  use diagnos_timetrace, only : allocate_timetrace => allocate_mem
  use diagnos_cross_phase, only : allocate_cross_phase => allocate_mem

  use global, only : BPAR_GA_FIELD
  use diagnos_generic, only : MU_GHOSTCELLS
  use general, only : gkw_abort
  integer :: ierr
  
  allocate(requirements(BPAR_GA_FIELD,MU_GHOSTCELLS), stat=ierr)
  if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate requirements')
  requirements = .false.
  
  call allocate_generic
  
  call allocate_energetics
  call allocate_corr
  call allocate_eng
  call allocate_f
  call allocate_fields
  call allocate_fluxes
  call allocate_kinenergy_trappas
  call allocate_fluxes_vspace
  call allocate_grid
  call allocate_growth_freq
  call allocate_mode_struct
  call allocate_moments
  call allocate_rad
  call allocate_stresses
  call allocate_velspace
  call allocate_jdote
  call allocate_zfshear
  call allocate_zonal_evo
  call allocate_nonlin_transfer
  call allocate_timetrace
  call allocate_cross_phase

end subroutine diagnostic_allocate

!----------------------------------------------------------------------------
!> output routine called every NAVERAGE timsteps
!----------------------------------------------------------------------------
subroutine diagnostic_naverage()

  use control, only : nt_complete, itime, itime_rst
  use diagnos_grid, only : output_grid => output
  use diagnos_growth_freq, only : output_growth_freq => output, &
     & calc_growth_freq => calc_largestep
  use diagnos_energetics, only : output_energetics => output
  use diagnos_corr, only : output_corr => output
  use diagnos_nonlin_transfer, only : output_nonlin_transfer => output
  use diagnos_fluxes, only : output_fluxes => output
  use diagnos_kinenergy_trappas, only : output_kinenergy_trappas => output
  use diagnos_fluxes_vspace, only : output_fluxes_vspace => output
  use diagnos_eng, only : output_eng => output
  use diagnos_rad, only : output_rad => output
  use diagnos_stresses, only : output_stresses => output
  use diagnos_velspace, only : output_velspace => output
  use diagnos_jdote, only : output_jdote => output
  use diagnos_fields, only : output_fields => output
  use diagnos_moments, only : output_moments => output
  use diagnos_f, only : output_f => output
  use diagnos_mode_struct, only : output_mode_struct => output
  use diagnos_zfshear, only : output_zfshear => output
  use diagnos_zonal_evo, only : output_zonal_evo => output
  use diagnos_timetrace, only : output_timetrace => output
  use diagnos_cross_phase, only : output_cross_phase => output
  use io, only : flush_all_open_files
  use perform, only : perfon, perfoff, perfswitch, perf_measure

  use diagnos_generic, only : lwrite_output1

  !> counter, in order to flush files occasionally
  integer, save :: next_flush_ts = 0

  ! The following is used to achieve consistent numbering after restarts
  ! for diagnostics which produce a series of numbered files.
  ! itime_rst counts the number of large time steps since the last
  ! restart file has been written (and therefore nt_complete has been 
  ! updated). It is required to allow for a consistent numbering with 
  ! automatically written dump-files during runtime.
  file_count = nt_complete + itime_rst

  if(perf_measure) then
    call perfon("diagnostics total",3)
    call perfon("diagnos_growth_freq",3)
  end if

  ! CHECK is this really necessary? and the first call to this routine
  ! from eiv_integration?
  call calc_growth_freq
  if(perf_measure) call perfoff(3)

  if(lwrite_output1) then

    call fill_buffers()
    if(perf_measure) call perfon("diagnos_growth_freq",3)
    call output_growth_freq
    if(perf_measure) call perfswitch("diagnos_grid",3)
    call output_grid(file_count)
    if(perf_measure) call perfswitch("diagnos_fluxes",3)
    call output_fluxes(.false., file_count)
    if(perf_measure) call perfswitch("diagnos_kinenergy_trappas",3)
    call output_kinenergy_trappas(file_count)
    if(perf_measure) call perfswitch("diagnos_energetics",3)
    call output_energetics()
    if(perf_measure) call perfswitch("diagnos_corr",3)
    call output_corr()
    if(perf_measure) call perfswitch("diagnos_nonlin_transfer",3)
    call output_nonlin_transfer()
    if(perf_measure) call perfswitch("diagnos_stresses",3)
    call output_stresses()
    if(perf_measure) call perfswitch("diagnos_fields",3)
    call output_fields(file_count)
    if(perf_measure) call perfswitch("diagnos_moments",3)
    call output_moments(file_count)
    if(perf_measure) call perfswitch("diagnos_velspace",3)
    call output_velspace(file_count)
    if(perf_measure) call perfswitch("diagnos_jdote",3)
    call output_jdote(file_count)
    if(perf_measure) call perfswitch("diagnos_mode_struct",3)
    call output_mode_struct
    if(perf_measure) call perfswitch("diagnos_f",3)
    call output_f(file_count)
    if(perf_measure) call perfswitch("diagnos_eng",3)
    call output_eng
    if(perf_measure) call perfswitch("diagnos_rad",3)
    call output_rad
    if(perf_measure) call perfswitch("diagnos_fluxes_vspace",3)
    call output_fluxes_vspace
    if(perf_measure) call perfswitch("diagnos_zfshear",3)
    call output_zfshear
    if(perf_measure) call perfswitch("diagnos_zonal_evo",3)
    call output_zonal_evo
    if(perf_measure) call perfswitch("diagnos_cross_phase",3)
    call output_cross_phase
    if(perf_measure) call perfswitch("diagnos_timetrace",3)
    call output_timetrace

    if(perf_measure) call perfoff(3)

    ! all processes enter the screen_output routines, but at the lowest
    ! level only one will write to the terminal. This is done in this
    ! way to do it analogously to the logical-unit IO, and thus to
    ! simplify things.
    call write_screen_output

    ! Recompute fluxes suprema
    ! WARNING: do use any results computed in output_fluxes after this call
    if(perf_measure) call perfon("diagnos_fluxes",3)
    call output_fluxes(.true., file_count)
    
    if(perf_measure) call perfswitch("diagnostic flush",3)
    ! occasionally flush files.
    ! flush the main ts files if nflush_ts > 0
    next_flush_ts = next_flush_ts + 1 
    if (next_flush_ts == nflush_ts .and. nflush_ts > 0) then
      next_flush_ts = 0
      call flush_all_open_files
    end if
    if(perf_measure) call perfoff(3)

  end if
  
  if(perf_measure) call perfoff(3)

end subroutine diagnostic_naverage

!----------------------------------------------------------------------------
!>
!----------------------------------------------------------------------------
subroutine fill_buffers()
  use dist, only : fdisi,fdis_tmp
  use dist, only : get_phi, get_bpar, get_apar, apar, bpar, phi
  use control, only : nlphi, nlapar, nlbpar
  use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD, DISTRIBUTION
  use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
  use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS, X_GHOSTCELLS
  use diagnos_generic, only : VPAR_GHOSTCELLS, MU_GHOSTCELLS
  use mpighosts, only : mpistart, mpiwait, gc1, gc1x2, gc1x_phi
  use general, only : gkw_abort

  logical :: fdis_tmp_is_updated
  fdis_tmp_is_updated = .false.

  if (requirements(DISTRIBUTION,LOCAL_DATA)) then
    ! this is trivially available, always.
  end if

  if(any(requirements(PHI_FIELD:BPAR_FIELD, VPAR_GHOSTCELLS:MU_GHOSTCELLS)) .or.&
     & any(requirements(PHI_GA_FIELD:BPAR_GA_FIELD, VPAR_GHOSTCELLS)) &
     & ) then
    call gkw_abort("diagnostics: this requirement does not make sense")
  end if

  ! copy data and communicate

  if(any(requirements(:,S_GHOSTCELLS)) .or. &
     & any(requirements(:,VPAR_GHOSTCELLS)) .or. &
     & any(requirements(:,MU_GHOSTCELLS)) &
     & ) then
    ! Send/Recv the distribution function to/from neighbours into fdis_tmp
    ! (vpar, s and mu ghost points only)
    call fill_fdis_tmp(fdis_tmp,fdis_tmp_is_updated)
    call mpistart(gc1)
    call mpiwait(gc1)
  end if

  if(requirements(DISTRIBUTION, X_GHOSTCELLS) .or. &
     & requirements(PHI_GA_FIELD, X_GHOSTCELLS)) then
    ! Send/Recv the distribution function and phi_ga (two x ghost
    ! points only) into fdis_tmp
    call fill_fdis_tmp(fdis_tmp,fdis_tmp_is_updated)
    call mpistart(gc1x2)
    call mpiwait(gc1x2)
  end if

  if(any(requirements(PHI_FIELD:BPAR_FIELD, X_GHOSTCELLS)) .or. &
     & any(requirements(PHI_GA_FIELD:BPAR_GA_FIELD, X_GHOSTCELLS))) then
    ! Send/Recv phi (two x ghost
    ! points only) into fdis_tmp
    call fill_fdis_tmp(fdis_tmp,fdis_tmp_is_updated)
    call mpistart(gc1x_phi)
    call mpiwait(gc1x_phi)
  end if

  if(any(requirements(APAR_GA_FIELD, VPAR_GHOSTCELLS:MU_GHOSTCELLS)) .or. &
     & any(requirements(BPAR_GA_FIELD, VPAR_GHOSTCELLS:MU_GHOSTCELLS)) &
     & ) then
    ! this is not communicated by gc1
    call gkw_abort("diagnostics: this requirement is not yet implemented")
  end if

  ! fill the field buffers

  if (nlphi) then
    if(requirements(PHI_FIELD,X_GHOSTCELLS) .or. &
       & requirements(PHI_FIELD,S_GHOSTCELLS)) then
      call get_phi(fdis_tmp,phi)
    else if(requirements(PHI_FIELD,LOCAL_DATA)) then
      call get_phi(fdisi,phi)
    end if
  end if
  if (nlapar) then
    if(requirements(APAR_FIELD,X_GHOSTCELLS) .or. &
       & requirements(APAR_FIELD,S_GHOSTCELLS)) then
      call get_apar(fdis_tmp,apar)
    else if(requirements(APAR_FIELD,LOCAL_DATA)) then
      call get_apar(fdisi,apar)
    end if
  end if
  if (nlbpar) then
    if(requirements(BPAR_FIELD,X_GHOSTCELLS) .or. &
       & requirements(BPAR_FIELD,S_GHOSTCELLS)) then
      call get_bpar(fdis_tmp,bpar)
    else if(requirements(BPAR_FIELD,LOCAL_DATA)) then
      call get_bpar(fdisi,bpar)
    end if
  end if

  contains
    subroutine fill_fdis_tmp(fdis_tmp_, is_updated)
      use matdat, only : matg2f
      use dist, only : fdisi, nsolc
      complex, intent(inout) :: fdis_tmp_(:)
      logical, intent(inout) :: is_updated
      integer :: i

      if(is_updated) return

      fdis_tmp_ = (0.,0.)
      if (nlapar) then
        ! For electro-magnetic runs, undo the A|| correction of g at the
        ! same time.
        do i = 1, matg2f%nmat
          fdis_tmp_(i) = fdisi(i) + matg2f%mat(i)*fdisi(matg2f%jj(i))
        end do
        do i= matg2f%nmat+1, nsolc
          fdis_tmp_(i) = fdisi(i)
        end do
      else
        fdis_tmp_(1:nsolc) = fdisi(1:nsolc)
      end if
      
      is_updated = .true.
      
    end subroutine fill_fdis_tmp

end subroutine fill_buffers


!----------------------------------------------------------------------------
!> Some quantities may be interesting to know right after the initialisation.
!> They can be output here.
!----------------------------------------------------------------------------
subroutine diagnostic_initial_output()
  use diagnos_grid, only : initial_output_grid => initial_output
  use diagnos_growth_freq, only : initial_output_growth_freq => initial_output
  use diagnos_fields, only : initial_output_fields => initial_output
  use diagnos_mode_struct, only : &
     & initial_output_mode_struct => initial_output
  use diagnos_cross_phase, only : &
     & initial_output_cross_phase => initial_output
  use diagnos_rad, only : initial_output_rad => initial_output
  use diagnos_neoequil, only : initial_output_neoequil => initial_output
  use diagnos_matrix, only : initial_output_matrix => initial_output

  call fill_buffers()

  call initial_output_grid
  call initial_output_growth_freq
  call initial_output_fields
  call initial_output_cross_phase
  call initial_output_neoequil
  call initial_output_matrix
    
  ! make the rest dependend of the lfinal_output switch for the moment.
  if (.not. lfinal_output) return

  call initial_output_mode_struct
  call initial_output_rad

end subroutine diagnostic_initial_output


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> deal with final output after time integration has been completed 
!----------------------------------------------------------------------------
subroutine diagnostic_final_output(number)
  use diagnos_velspace, only : final_output_velspace => final_output
  use diagnos_fields, only : final_output_fields => final_output
  use diagnos_moments, only : final_output_moments => final_output
  use diagnos_fluxes, only : final_output_fluxes => final_output
  use diagnos_fluxes_vspace, only : final_output_fluxes_vspace => final_output
  use diagnos_mode_struct, only : final_output_mode_struct => final_output
  use diagnos_f, only : final_output_f => final_output
  use diagnos_cross_phase, only : final_output_cross_phase => final_output
  
  !> The optional integer number can be given to produce numbered files.
  integer, optional, intent(in) :: number

  if (.not. lfinal_output) return

  call fill_buffers()

  if (present(number)) then
    call final_output_fluxes(number)
    call final_output_fields(number)
    call final_output_moments(number)
    call final_output_f(number)
    call final_output_mode_struct(number)
    call final_output_cross_phase(number)
  else
    call final_output_fluxes(file_count)
    call final_output_fields(file_count)
    call final_output_moments(file_count)
    call final_output_f()
    call final_output_mode_struct()
    call final_output_cross_phase()
  end if
  call final_output_velspace()
  call final_output_fluxes_vspace
end subroutine diagnostic_final_output


!----------------------------------------------------------------------------
!> perform any initialisation that needs to be done before the main time loop
!----------------------------------------------------------------------------
subroutine diagnostic_init
  use control, only : testing, spectral_radius
  use diagnos_generic, only : init_generic => init
  use diagnos_grid, only : init_grid => init
  use diagnos_growth_freq, only : init_growth_freq => init
  use diagnos_fluxes, only : init_fluxes => init
  use diagnos_kinenergy_trappas, only : init_kinenergy_trappas => init
  use diagnos_fluxes_vspace, only : init_fluxes_vspace => init
  use diagnos_fields, only : init_fields => init
  use diagnos_velspace, only : init_velspace => init
  use diagnos_jdote, only : init_jdote => init
  use diagnos_energetics, only : init_energetics => init
  use diagnos_corr, only : init_corr => init
  use diagnos_rad, only : init_rad => init
  use diagnos_mode_struct, only : init_mode_struct => init
  use diagnos_moments, only : init_moments => init
  use diagnos_nonlin_transfer, only : init_nonlin_transfer => init
  use diagnos_eng, only : init_eng => init
  use diagnos_stresses, only : init_stresses => init
  use diagnos_zfshear, only : init_zfshear => init
  use diagnos_zonal_evo, only : init_zonal_evo => init
  use diagnos_cross_phase, only : init_cross_phase => init
  use diagnos_timetrace, only : init_timetrace => init
  use diagnos_f, only : init_f => init

  use diagnos_generic, only : lwrite_output1
  use global, only : PHI_FIELD, PHI_GA_FIELD, BPAR_FIELD, BPAR_GA_FIELD
  use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS, X_GHOSTCELLS
  
  
  ! no diagnostic files if testing MPI
  if (testing) then
    ! In order not to output any data to files, the io_format = 'none'
    ! is already set.

    ! additional flags spare us the computations of some diagnostics
    lwrite_output1 = .false.
    lfinal_output  = .false.
  end if

  call init_generic()
  call init_growth_freq(requirements)
  call init_grid(requirements)
  call init_fluxes(requirements)
  call init_kinenergy_trappas(requirements)
  call init_fluxes_vspace(requirements)
  call init_fields(requirements)
  call init_mode_struct(requirements)
  call init_moments(requirements)
  call init_nonlin_transfer(requirements)
  call init_velspace(requirements)
  call init_jdote(requirements)
  call init_energetics(requirements)
  call init_corr(requirements)
  call init_eng(requirements)
  call init_rad(requirements)
  call init_stresses(requirements)
  call init_zfshear(requirements)
  call init_zonal_evo(requirements)
  call init_cross_phase(requirements)
  call init_timetrace(requirements)
  call init_f(requirements)
  

  if(spectral_radius) then
    ! There is no phi_ga space in memory but the phi_ga field is
    ! simply computed by bessel*field. Ghostcells of phi_ga are
    ! ghostcells of phi.
    requirements(PHI_FIELD:BPAR_FIELD,LOCAL_DATA) = requirements(PHI_FIELD:BPAR_FIELD,LOCAL_DATA) &
       & .or. requirements(PHI_GA_FIELD:BPAR_GA_FIELD,LOCAL_DATA)
    requirements(PHI_FIELD:BPAR_FIELD,S_GHOSTCELLS) = requirements(PHI_FIELD:BPAR_FIELD,S_GHOSTCELLS) &
       & .or. requirements(PHI_GA_FIELD:BPAR_GA_FIELD,S_GHOSTCELLS)
    requirements(PHI_FIELD:BPAR_FIELD,X_GHOSTCELLS) = requirements(PHI_FIELD:BPAR_FIELD,X_GHOSTCELLS) &
       & .or. requirements(PHI_GA_FIELD:BPAR_GA_FIELD,X_GHOSTCELLS)
    requirements(PHI_GA_FIELD:BPAR_GA_FIELD,:) = .false.
  end if

end subroutine diagnostic_init

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> If earlier data files exist the run is terminated.
!> This prevents confusion of appending / overwriting data of an older run.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine abort_if_old_data_exists
  use general, only : gkw_abort
  use mpiinterface, only : mpibarrier

  if (check_for_old_data()) then
    call gkw_abort('Existing data found. &
       & Please remove previous run data (script gkw_clean_run) &
       & or set read_file=.true. to restart') 
  end if

  call mpibarrier

end subroutine abort_if_old_data_exists


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Write output to the screen.
!----------------------------------------------------------------------------
subroutine write_screen_output
  use diagnos_growth_freq, only : screen_output_growth_freq => screen_output
  use diagnos_fluxes, only : screen_output_fluxes => screen_output
  
  if (screen_output) then
    call screen_output_growth_freq
    call screen_output_fluxes
  end if

end subroutine write_screen_output


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Let diagnostics deallocate memory, close the logical units, etc.
!----------------------------------------------------------------------------
subroutine diagnostic_finalize
  use io, only : close_all_lus
  use diagnos_generic, only : finalize_generic => finalize
  use diagnos_eng, only : finalize_eng => finalize
  use diagnos_rad, only : finalize_rad => finalize
  use diagnos_energetics, only : finalize_energetics => finalize
   use diagnos_corr, only : finalize_corr => finalize
  use diagnos_growth_freq, only : finalize_growth_freq => finalize
  use diagnos_velspace, only : finalize_velspace => finalize
  use diagnos_jdote, only : finalize_jdote => finalize
  use diagnos_zfshear, only : finalize_zfshear => finalize
  use diagnos_zonal_evo, only : finalize_zonal_evo => finalize
  use diagnos_nonlin_transfer, only : finalize_nonlin_transfer => finalize
  use diagnos_timetrace, only : finalize_timetrace => finalize
  use diagnos_cross_phase, only : finalize_cross_phase => finalize

  call finalize_generic()
  call finalize_eng()
  call finalize_energetics()
  call finalize_corr()
  call finalize_growth_freq()
  call finalize_rad()
  call finalize_velspace()
  call finalize_jdote()
  call finalize_zfshear()
  call finalize_zonal_evo()
  call finalize_nonlin_transfer()
  call finalize_timetrace()
  call finalize_cross_phase()

  ! close all logical units, in case the diagnostics are lazy and do
  ! not do it themselves.
  call close_all_lus

end subroutine diagnostic_finalize


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Let diagnostics read in the last chunk of data which they produced
!> in the preceeding run.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine diagnostic_read_last_data
!  use diagnos_growth_freq, only : read_last_data_growth_freq => read_last_data
!  use diagnos_grid, only : read_last_data_grid => read_last_data

  ! call read_last_data_grid
  ! call read_last_data_growth_freq
  
end subroutine diagnostic_read_last_data

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine checks for the existence of a sample of earlier data
!> files.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function check_for_old_data()
  use io, only : lu_exists, ascii_fmt
  use control, only : io_legacy
  logical :: check_for_old_data
  logical, dimension(9) :: checks

  checks(:) = .false.
  if(io_legacy) then
    checks(1) = lu_exists('time', '/diagnostic/diagnos_growth_freq', &
       & ascii_fmt)
    checks(2) = lu_exists('fluxes', '/diagnostic/diagnos_fluxes', ascii_fmt)
  else
    checks(1) = lu_exists('time', '/grid/', ascii_fmt)
    checks(2) = lu_exists('eflux', '/diagnostic/diagnos_fluxes', ascii_fmt)
  end if
  
  checks(3) = lu_exists('eflux_spectra', '/diagnostic/diagnos_fluxes', &
     & ascii_fmt)
  checks(4) = lu_exists('vflux_spectra', '/diagnostic/diagnos_fluxes', &
     & ascii_fmt)
  checks(5) = lu_exists('pflux_spectra', '/diagnostic/diagnos_fluxes', &
     & ascii_fmt)
  checks(6) = lu_exists('fluxes_nc', '/diagnostic/diagnos_fluxes', ascii_fmt)
  checks(7) = lu_exists('kyspec', '/diagnostic/diagnos_fields', ascii_fmt)
  checks(8) = lu_exists('kxspec', '/diagnostic/diagnos_fields', ascii_fmt)
  checks(9) = lu_exists('parallel_phi', '/diagnostic/diagnos_fields', &
     & ascii_fmt)

  check_for_old_data = any(checks)
end function check_for_old_data


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> some quantities need to be calculated at 2 consecutive timesteps, so this
!> routine is called the timestep before naverage, AND on naverage
!> in order to update such quantities
!> Further remarks:
!> At the moment fdis_tmp is updated in advance_large_step_explicit() just before
!> diagnostics_naverage() is called. If one ever needs fdis_tmp in one of
!> the pre_diagnostics, this must be changed.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine diagnostic_pre_naverage(i_smallstep, fdisk)
  use control, only : naverage, spectral_radius
  use diagnos_grid, only : calc_smallstep_grid => calc_smallstep
  use diagnos_growth_freq, only : calc_smallstep_growth_freq => calc_smallstep
  use diagnos_energetics, only : calc_smallstep_energetics => calc_smallstep
  use diagnos_zonal_evo, only : calc_smallstep_zonal_evo => calc_smallstep
  use diagnos_jdote, only : calc_smallstep_jdote => calc_smallstep
  use diagnos_timetrace, only : lcomplex_timetrace
  use diagnos_timetrace, only : calc_smallstep_timetrace => calc_smallstep
  use dist, only : fdisi
  use fields, only : calculate_fields, field_solve_nonspec_wrap

  integer, intent(in) :: i_smallstep
  complex, intent(in) :: fdisk(:)

  
  if (i_smallstep > naverage - 2) then
    ! Some diagnostics (such as computation of frequencies) need
    ! values only from the previous small time step...
    if (spectral_radius) then
      call calculate_fields(fdisi)
    else
      call field_solve_nonspec_wrap(fdisi,0.,.false.)
    endif
    ! Note that for naverage >= 2 this block will be done in two
    ! consecutive small time steps
    call calc_smallstep_growth_freq(i_smallstep)
    call calc_smallstep_energetics(i_smallstep)
    call calc_smallstep_jdote(i_smallstep)
    call calc_smallstep_zonal_evo(i_smallstep)
    
  else if(lcomplex_timetrace) then
    ! for simplicity, make an exception here and use the switch from
    ! the diagnostic to decide if fields need to be updated here
    ! for every small timestep.

    !...  whereas others have to be called at all small time steps.
    if (spectral_radius) then
      call calculate_fields(fdisi)
    else
      call field_solve_nonspec_wrap(fdisi,0.,.false.)
    endif

  end if

  call calc_smallstep_grid(i_smallstep)
  call calc_smallstep_timetrace(fdisk)

end subroutine diagnostic_pre_naverage

end module diagnostic

