##############################################################################
### GKW object files.
##############################################################################

F90OBJ =  gkw.o init.o perform.o index_function.o constants.o mpicomms.o  \
          specfun.o marsaglia.o components.o fft.o switches.o   \
          exp_integration.o geom.o control.o general.o grid.o io.o  \
          matdat.o matrix_format.o \
          functions.o dist.o mode.o linear_terms.o non_linear_terms.o        \
          imp_integration.o normalise.o mpiinterface.o mpidatatypes.o        \
          collisionop.o global.o velocitygrid.o rotation.o tearingmodes.o    \
          diagnostic.o \
          diagnos_generic.o diagnos_growth_freq.o \
          diagnos_rad.o diagnos_eng.o diagnos_energetics.o \
          diagnos_corr.o diagnos_velspace.o diagnos_jdote.o \
          diagnos_fluxes.o diagnos_fluxes_vspace.o \
          diagnos_kinenergy_trappas.o \
          diagnos_fields.o diagnos_moments.o diagnos_mode_struct.o \
          diagnos_grid.o \
          diagnos_f.o diagnos_stresses.o diagnos_zfshear.o \
          diagnos_nonlin_transfer.o diagnos_timetrace.o \
          diagnos_cross_phase.o diagnos_neoequil.o diagnos_matrix.o \
          diagnos_zonal_evo.o \
          ompinterface.o gauss.o gyro_average.o fields.o \
          structures.o restart.o mpighosts.o \
          matconv.o krook.o source_time.o eiv_integration.o rho_par_switch.o \
          version.o io_ascii.o io_binary.o io_hdf5.o neoequil.o

## objects required for the main program
OBJLIST = $(F90OBJ) $(UMF_INTERFACE) $(MKL_SPBLAS_MODULE)

## objects required for the input checking executable
CHECKOBJLIST = input_check.o $(filter-out gkw.o,$(OBJLIST))
