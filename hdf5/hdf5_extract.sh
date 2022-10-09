#!/bin/bash

INPUTFILE="data.h5"
OUTPUTFILE="test.h5"

DATASETS=(
#"diagnostic/diagnos_fields/kxspec"
#"diagnostic/diagnos_fields/kxvort"
#"diagnostic/diagnos_fields/kyspec"
#"diagnostic/diagnos_fields/kyvort"
#"diagnostic/diagnos_fields/phi"
#"diagnostic/diagnos_fields/spc"
#"diagnostic/diagnos_fluxes/EFlesr0001"
#"diagnostic/diagnos_fluxes/eflux_species01"
#"diagnostic/diagnos_fluxes/eflux_spectra"
#"diagnostic/diagnos_fluxes/eflux_sup"
#"diagnostic/diagnos_fluxes/eflux_xspec"
#"diagnostic/diagnos_fluxes/flmgr01"
#"diagnostic/diagnos_fluxes/pflux_species01"
#"diagnostic/diagnos_fluxes/pflux_spectra"
#"diagnostic/diagnos_fluxes/pflux_sup"
#"diagnostic/diagnos_fluxes/pflux_xspec"
#"diagnostic/diagnos_fluxes/vflux_species01"
#"diagnostic/diagnos_fluxes/vflux_spectra"
#"diagnostic/diagnos_fluxes/vflux_xspec"
#"diagnostic/diagnos_grid/intmu"
#"diagnostic/diagnos_grid/intvp"
#"diagnostic/diagnos_grid/lxn"
#"diagnostic/diagnos_grid/lyn"
#"diagnostic/diagnos_grid/mode_label"
#"diagnostic/diagnos_grid/mphi"
#"diagnostic/diagnos_grid/mphiw3"
#"diagnostic/diagnos_grid/mrad_G"
#"diagnostic/diagnos_grid/mrad_l"
#"diagnostic/diagnos_grid/sgrid"
#"diagnostic/diagnos_growth_freq/frequencies"
#"diagnostic/diagnos_growth_freq/frequencies_all_modes"
#"diagnostic/diagnos_growth_freq/growth"
#"diagnostic/diagnos_growth_freq/growth_rates_all_modes"
#"diagnostic/diagnos_growth_freq/time"
#"diagnostic/diagnos_mode_struct/parallel"
#"diagnostic/diagnos_moments/den01"
#"diagnostic/diagnos_moments/den_spectra"
#"diagnostic/diagnos_moments/ene01"
#"diagnostic/diagnos_moments/ene_spectra"
#"diagnostic/diagnos_rad/prof_back"
#"geom/Bref"
#"geom/Bt_frac"
#"geom/D_eps"
#"geom/D_s"
#"geom/D_zeta"
#"geom/E_eps_s"
#"geom/E_eps_zeta"
#"geom/E_zeta_s"
#"geom/F"
#"geom/G"
#"geom/H_eps"
#"geom/H_s"
#"geom/H_zeta"
#"geom/I_eps"
#"geom/I_s"
#"geom/I_zeta"
#"geom/J"
#"geom/Jacobian"
#"geom/K"
#"geom/NS"
#"geom/R"
#"geom/R0"
#"geom/Rref"
#"geom/Z"
#"geom/beta_eq"
#"geom/betaprime_eq"
#"geom/bmax"
#"geom/bmin"
#"geom/bn"
#"geom/eps"
#"geom/g_eps_eps"
#"geom/g_eps_s"
#"geom/g_eps_zeta"
#"geom/g_s_s"
#"geom/g_zeta_s"
#"geom/g_zeta_zeta"
#"geom/jfunh"
#"geom/jfunl"
#"geom/krnorm"
#"geom/kthnorm"
#"geom/lfun"
#"geom/poloidal_angle"
#"geom/q"
#"geom/s_grid"
#"geom/shat"
#"grid/file_count"
#"grid/krho"
#"grid/krho_extended"
#"grid/krloc"
#"grid/kxrh"
#"grid/kzeta"
#"grid/time_fine"
#"grid/vperp"
#"grid/vpgr"
#"grid/xgr"
#"grid/xphi"
#"grid/yphi"
#"input/collisions/coll_freq"
#"input/collisions/cons_type"
#"input/collisions/en_scatter"
#"input/collisions/ene_conservation"
#"input/collisions/freq_override"
#"input/collisions/friction_coll"
#"input/collisions/lorentz"
#"input/collisions/mass_conserve"
#"input/collisions/mom_conservation"
#"input/collisions/nref"
#"input/collisions/pitch_angle"
#"input/collisions/rref"
#"input/collisions/selfcollcon"
#"input/collisions/tref"
#"input/collisions/zeff"
#"input/control/auto_restart"
#"input/control/collisions"
#"input/control/disp_par"
#"input/control/disp_vp"
#"input/control/disp_x"
#"input/control/disp_y"
#"input/control/dt_min"
#"input/control/dtim"
#"input/control/fac_dtim_est"
#"input/control/flux_tube"
#"input/control/fluxtol"
#"input/control/gamatol"
#"input/control/ifluxtol"
#"input/control/io_format"
#"input/control/io_legacy"
#"input/control/io_testdata"
#"input/control/iperform_set"
#"input/control/irun"
#"input/control/laverage_dist_over_time"
#"input/control/lflapv"
#"input/control/lpar_vel_nl"
#"input/control/lrestart_new_grid"
#"input/control/ltrapping_arakawa"
#"input/control/lverbose"
#"input/control/matrix_format"
#"input/control/max_gr"
#"input/control/max_sec"
#"input/control/max_seconds"
#"input/control/meth"
#"input/control/method"
#"input/control/min_gr"
#"input/control/naverage"
#"input/control/ncqtol"
#"input/control/ndump_ts"
#"input/control/neoclassics"
#"input/control/nl_dtim_est"
#"input/control/nlapar"
#"input/control/nlbpar"
#"input/control/nlphi"
#"input/control/non_linear"
#"input/control/normalize_per_toroidal_mode"
#"input/control/normalized"
#"input/control/ntime"
#"input/control/order_of_the_radial_scheme"
#"input/control/order_of_the_scheme"
#"input/control/order_of_the_zf_scheme"
#"input/control/parallel_boundary_conditions"
#"input/control/radial_boundary_conditions"
#"input/control/read_file"
#"input/control/restart_file_version"
#"input/control/shift_metric"
#"input/control/silent"
#"input/control/spectral_radius"
#"input/control/testing"
#"input/control/uniform_mu_grid"
#"input/control/vp_trap"
#"input/control/zonal_adiabatic"
#"input/diagnostic/apa3d"
#"input/diagnostic/apc3d"
#"input/diagnostic/bpa3d"
#"input/diagnostic/bpc3d"
#"input/diagnostic/cross_phase_timetrace"
#"input/diagnostic/den3d"
#"input/diagnostic/ene3d"
#"input/diagnostic/field_fsa_kyzero"
#"input/diagnostic/flux3d"
#"input/diagnostic/imod_corr"
#"input/diagnostic/imod_field"
#"input/diagnostic/kykxs_apar"
#"input/diagnostic/kykxs_bpar"
#"input/diagnostic/kykxs_j0_moments"
#"input/diagnostic/kykxs_j1_moments"
#"input/diagnostic/kykxs_moments"
#"input/diagnostic/kykxs_phi"
#"input/diagnostic/lamplitudes"
#"input/diagnostic/lcalc_corr"
#"input/diagnostic/lcalc_energetics"
#"input/diagnostic/lcalc_fluxes"
#"input/diagnostic/lcalc_jdote"
#"input/diagnostic/lcalc_jdote_fs"
#"input/diagnostic/lcalc_kinenergy_trappas"
#"input/diagnostic/lcalc_stresses"
#"input/diagnostic/lcalc_tot_energy"
#"input/diagnostic/lcalc_zfshear"
#"input/diagnostic/lcalc_zonal_evo"
#"input/diagnostic/lfinal_output"
#"input/diagnostic/lfinvel"
#"input/diagnostic/lfluxes_detail"
#"input/diagnostic/lfluxes_em_spectra"
#"input/diagnostic/lfluxes_spectra"
#"input/diagnostic/lfluxes_vspace"
#"input/diagnostic/lfluxes_vspace_bpar"
#"input/diagnostic/lfluxes_vspace_em"
#"input/diagnostic/lfluxes_vspace_phi"
#"input/diagnostic/lfrequencies"
#"input/diagnostic/lgrowth_rates"
#"input/diagnostic/lisl_follow"
#"input/diagnostic/lmode_energy"
#"input/diagnostic/lmpi_broken_io"
#"input/diagnostic/lnonlin_transfer"
#"input/diagnostic/lnonlin_transfer_fsa"
#"input/diagnostic/lparallel_output"
#"input/diagnostic/lphi_diagnostics"
#"input/diagnostic/lrad_entropy"
#"input/diagnostic/lrad_field"
#"input/diagnostic/lrad_kpar"
#"input/diagnostic/lrad_moment"
#"input/diagnostic/lrad_tint"
#"input/diagnostic/lradial_entropy"
#"input/diagnostic/lradial_profile"
#"input/diagnostic/lrotate_parallel"
#"input/diagnostic/lvelspace_output"
#"input/diagnostic/lwrite_output1"
#"input/diagnostic/nflush_ts"
#"input/diagnostic/nmodepoints"
#"input/diagnostic/nonlin_transfer_interval"
#"input/diagnostic/npointsvel"
#"input/diagnostic/out3d_interval"
#"input/diagnostic/parallel_output_timestamps"
#"input/diagnostic/phi3d"
#"input/diagnostic/psi_velspace"
#"input/diagnostic/screen_output"
#"input/diagnostic/spc3d"
#"input/diagnostic/xs_kyzero_current"
#"input/diagnostic/xs_kyzero_current2"
#"input/diagnostic/xs_kyzero_dens"
#"input/diagnostic/xs_kyzero_ene"
#"input/diagnostic/xs_kyzero_ene_par"
#"input/diagnostic/xs_kyzero_ene_perp"
#"input/diagnostic/xs_kyzero_phi_ga_deltaf"
#"input/diagnostic/xs_kyzero_phi_ga_fm"
#"input/diagnostic/xs_phi"
#"input/diagnostic/xy_apar"
#"input/diagnostic/xy_bpar"
#"input/diagnostic/xy_current"
#"input/diagnostic/xy_current2"
#"input/diagnostic/xy_dens"
#"input/diagnostic/xy_estep"
#"input/diagnostic/xy_fluxes"
#"input/diagnostic/xy_fluxes_bi"
#"input/diagnostic/xy_fluxes_bpar"
#"input/diagnostic/xy_fluxes_em"
#"input/diagnostic/xy_fluxes_fsa"
#"input/diagnostic/xy_fluxes_k"
#"input/diagnostic/xy_fluxes_p"
#"input/diagnostic/xy_fluxes_v"
#"input/diagnostic/xy_phi"
#"input/diagnostic/xy_slice_ipar"
#"input/diagnostic/xy_spec"
#"input/diagnostic/xy_temp"
#"input/diagnostic/xy_vort"
#"input/diagnostic/zeta_velspace"
#"input/diagnostic/zevo_detail"
#"input/diagnostic/zevo_xy_spec"
#"input/diagnostic/zonal_scale_3d"
#"input/eiv_integration/freq"
#"input/eiv_integration/growthrate"
#"input/eiv_integration/mat_vec_routine"
#"input/eiv_integration/max_iterations"
#"input/eiv_integration/nr_column_vec"
#"input/eiv_integration/number_eigenvalues"
#"input/eiv_integration/tolerance"
#"input/eiv_integration/type_extraction"
#"input/eiv_integration/type_solver"
#"input/eiv_integration/which_eigenvalues"
#"input/finite_rho_parallel/lflux_rhostar"
#"input/finite_rho_parallel/lnonlinear_rhostar"
#"input/finite_rho_parallel/ltrapdf_rhostar"
#"input/finite_rho_parallel/lvdgrad_phi_fm_rhostar"
#"input/finite_rho_parallel/lvdgradf_rhostar"
#"input/finite_rho_parallel/lve_grad_fm_rhostar"
#"input/finite_rho_parallel/s_average"
#"input/geom/N_shape"
#"input/geom/R0_loc"
#"input/geom/Zmil"
#"input/geom/beta_rota_miller"
#"input/geom/beta_rota_miller_type"
#"input/geom/c"
#"input/geom/c_prime"
#"input/geom/curv_effect"
#"input/geom/dRmil"
#"input/geom/dZmil"
#"input/geom/delta"
#"input/geom/eps"
#"input/geom/eps_type"
#"input/geom/eqfile"
#"input/geom/geom_type"
#"input/geom/gradp"
#"input/geom/gradp_type"
#"input/geom/kappa"
#"input/geom/prof_type"
#"input/geom/q"
#"input/geom/qprof_coef"
#"input/geom/s"
#"input/geom/s_prime"
#"input/geom/sdelta"
#"input/geom/shat"
#"input/geom/signB"
#"input/geom/signJ"
#"input/geom/skappa"
#"input/geom/square"
#"input/geom/ssquare"
#"input/grid/lx"
#"input/grid/mumax"
#"input/grid/n_mu_grid"
#"input/grid/n_procs_mu"
#"input/grid/n_procs_s"
#"input/grid/n_procs_sp"
#"input/grid/n_procs_vpar"
#"input/grid/n_procs_x"
#"input/grid/n_s_grid"
#"input/grid/n_trapped"
#"input/grid/n_vpar_grid"
#"input/grid/n_x_grid"
#"input/grid/nmod"
#"input/grid/non_blocking_vpar"
#"input/grid/nperiod"
#"input/grid/number_of_species"
#"input/grid/nx"
#"input/grid/psih"
#"input/grid/psil"
#"input/grid/vpmax"
#"input/gyroaverage/blending_order"
#"input/gyroaverage/consistent_long_wave"
#"input/gyroaverage/gyro_average_electrons"
#"input/gyroaverage/gyro_average_ions"
#"input/gyroaverage/mk_gyro_av_hermitian"
#"input/gyroaverage/n_gav_bound_ex"
#"input/gyroaverage/n_points_ga"
#"input/gyroaverage/orb_polarize"
#"input/gyroaverage/parallel_mod"
#"input/gyroaverage/use_conj"
#"input/header/compiler"
#"input/header/gkw_executable_name"
#"input/header/gkw_version"
#"input/header/number_of_processors"
#"input/header/timestamp"
#"input/krook/bwidth"
#"input/krook/gamkpre"
#"input/krook/gammab"
#"input/krook/gammak"
#"input/krook/krook_option"
#"input/krook/nlbound"
#"input/krook/nlkrook"
#"input/linear_term_switches/apply_on_imod"
#"input/linear_term_switches/idisp"
#"input/linear_term_switches/lampere"
#"input/linear_term_switches/lbpar"
#"input/linear_term_switches/lg2f_correction"
#"input/linear_term_switches/lneo_equil"
#"input/linear_term_switches/lneo_rad"
#"input/linear_term_switches/lneo_trap"
#"input/linear_term_switches/lneoclassical"
#"input/linear_term_switches/lneorotsource"
#"input/linear_term_switches/lpoisson"
#"input/linear_term_switches/lpoisson_zf"
#"input/linear_term_switches/ltrapdf"
#"input/linear_term_switches/lvd_grad_phi_fm"
#"input/linear_term_switches/lvdgradf"
#"input/linear_term_switches/lve_grad_fm"
#"input/linear_term_switches/lvpar_grad_df"
#"input/linear_term_switches/lvpgrphi"
#"input/linear_term_switches/neo_equil_parse_sp_seq"
#"input/mode/chin"
#"input/mode/ikxspace"
#"input/mode/kr_type"
#"input/mode/krhomax"
#"input/mode/krrho"
#"input/mode/kthrho"
#"input/mode/mode_box"
#"input/mode/n_spacing"
#"input/mode/no_drive_of"
#"input/mode/no_transfer_from"
#"input/mode/no_transfer_to"
#"input/mode/rkxspace"
#"input/rotation/cf_drift"
#"input/rotation/cf_qncheck"
#"input/rotation/cf_trap"
#"input/rotation/cf_upphi"
#"input/rotation/cf_upsrc"
#"input/rotation/coriolis"
#"input/rotation/shear_profile"
#"input/rotation/shear_rate"
#"input/rotation/t_shear_begin"
#"input/rotation/toroidal_shear"
#"input/rotation/vcor"
#"input/source_time/dsfr"
#"input/source_time/gauss_source_median"
#"input/source_time/gauss_source_stdev"
#"input/source_time/mod_freq"
#"input/source_time/source_profile"
#"input/source_time/source_time_ampl"
#"input/source_time/source_wave_number"
#"input/spcgeneral/Ls"
#"input/spcgeneral/adiabatic_electrons"
#"input/spcgeneral/amp_imod"
#"input/spcgeneral/amp_imod_imag"
#"input/spcgeneral/amp_init"
#"input/spcgeneral/amp_zon"
#"input/spcgeneral/amp_zon_imag"
#"input/spcgeneral/beta"
#"input/spcgeneral/beta_ref"
#"input/spcgeneral/beta_type"
#"input/spcgeneral/betaprime_ref"
#"input/spcgeneral/betaprime_type"
#"input/spcgeneral/energetic_particles"
#"input/spcgeneral/finit"
#"input/spcgeneral/finit_imod"
#"input/spcgeneral/icrh_params"
#"input/spcgeneral/imod_init"
#"input/spcgeneral/init_coef"
#"input/spcgeneral/isl_mode"
#"input/spcgeneral/isl_rot_freq"
#"input/spcgeneral/isl_shear"
#"input/spcgeneral/lfinit_radial_dirichlet"
#"input/spcgeneral/mode_persist"
#"input/spcgeneral/n_quench"
#"input/spcgeneral/psi_0"
#"input/spcgeneral/quench_modes"
#"input/spcgeneral/quench_switch"
#"input/spcgeneral/rhostar"
#"input/spcgeneral/tear_zero_epar"
#"input/spcgeneral/tearingmode"
#"input/spcgeneral/tm_start"
#"input/spcgeneral/vpar_mean"
#"input/spcgeneral/wstar"
#"input/species01/background"
#"input/species01/dens"
#"input/species01/dens_prof_coef"
#"input/species01/dens_prof_type"
#"input/species01/mass"
#"input/species01/param"
#"input/species01/rln"
#"input/species01/rlt"
#"input/species01/rlt_gauss"
#"input/species01/temp"
#"input/species01/temp_prof_coef"
#"input/species01/temp_prof_type"
#"input/species01/uprim"
#"input/species01/z"
#"input/species02/background"
#"input/species02/dens"
#"input/species02/dens_prof_coef"
#"input/species02/dens_prof_type"
#"input/species02/mass"
#"input/species02/param"
#"input/species02/rln"
#"input/species02/rlt"
#"input/species02/rlt_gauss"
#"input/species02/temp"
#"input/species02/temp_prof_coef"
#"input/species02/temp_prof_type"
#"input/species02/uprim"
#"input/species02/z"
#"restart/dtim"
#"restart/nt_complete"
#"restart/nt_remain"
#"restart/time_complete"
#"evaluation/second_derivative_phi"
#"evaluation/zonalflow_potential"
#"evaluation/shearing_rate"
#"evaluation/shearing_rate_maximum"
"evaluation/derivative_stepsize"
)

for VAL in ${DATASETS[@]}; do
   h5copy -p -i "$INPUTFILE" -o "$OUTPUTFILE" -s "$VAL" -d "$VAL" 
done