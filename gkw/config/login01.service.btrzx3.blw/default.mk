## -- Bayreuth Rechenzentrum Cluster btrzx3 using Intel compiler

# need to load the following modules before compilation
# (the following are not enforced, just suggested):
REQUIREDMODULES = intel_parallel_studio_xe_2018_update2 libs/hdf5/5-1.10.5_serial_intel2018_2 libs/fftw/intel libs/atlas+lapack/

OPTFLAGS = on 
DEBUG = off 
PERF=perf

# Uncomment if you need the eigenvalue solver
SLEPC = HAVE_SLEPC
PETSC_DIR=/panasas/data/peeters/shared/petsc-current
SLEPC_DIR=/panasas/data/peeters/shared/slepc-current
FFLAGS_INC_SLEPC= -I$(PETSC_DIR)/include -I$(SLEPC_DIR)/include
LDFLAGS_SLEPC=-Wl,-rpath,${SLEPC_DIR}/lib -L${SLEPC_DIR}/lib -lslepc -Wl,-rpath,${PETSC_DIR}/lib -L${PETSC_DIR}/lib -lpetsc -lflapack -lfblas

# Uncomment if you need the HDF5 format - the paths do not yet work for everybody
IO_LIB = HAVE_HDF5
# The compiler flags for HDF5 are taken from the wrapper provided by HDF5
FFLAGS_INC_H5= $(shell h5fc -show | sed 's/^[^ ]*//')
LDFLAGS_H5=$(shell h5fc -show | sed 's/^[^ ]*//')

## Include the MKL library with (potentially faster)
## sparse linear algebra routines
#FFTLIB = FFT_MKL
#LIBMKL = HAVE_MKL
# MKL_HOME environment variables points to the directory where the MKL
# library is located. It should be set automatically when the module is loaded.
#FFLAGS_INC_MKL= -I$(MKL_HOME)/include/intel64/lp64  -I${MKL_HOME}/include -module ./
#LDFLAGS_MKL = -fpp -L "${MKL_HOME}/lib/intel64" -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -L "${MKL_HOME}/../compiler/lib/intel64" -liomp5 -lm -module ./
# The problem here is that there does not seem to be a mkl_spblas.mod
# file ready to use. We must compile it together with GKW. See also
# the corresponding target in gkw.mk .
#MKL_SPBLAS_MODULE_SRC = ${MKL_HOME}/include/mkl_spblas.f90

## Include the librsb library with (potentially faster)
## sparse linear algebra routines
#LIBRSB = HAVE_LIBRSB
## Those flags were copied from $(librsb-config --ldflags --extra_libs) from a local machine
#LDFLAGS_LIBRSB =  -L/panasas/data/peeters/shared/librsb-current/lib -Wl,-rpath -Wl,/panasas/data/peeters/shared/librsb-current/lib -lrsb -lm #-lgfortran
#FFLAGS_RSB = $(shell librsb-config --I_opts --static --ldflags --extra_libs)
# -L/panasas/data/peeters/shared/librsb-current/lib64 
#FFLAGS_INC_LIBRSB = -I/panasas/data/peeters/shared/librsb-current/include

# use the FFTW3 wrapper provided by MKL, i.e. use the MKL FFT, but
# with the bindings like FFTW3
# FFTLIB = FFT_FFTW3
# # Directory where the fft library is located
# FFLAGS_INC_MKL= -I$(MKL_HOME)/include/intel64/lp64
# LDFLAGS_MKL=-L${MKL_HOME}/lib/intel64/lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl
# FFLAGS_INC_FFT= -I${MKL_HOME}/include/fftw
# LDFLAGS_FFT=-L$(MKL_HMOE)/lib -lfftw3


FFTLIB = FFT_FFTW3
# Directory where the fft library is located, include and link paths
# are taken from the environment and set by the module system
FFLAGS_INC_FFT=
LDFLAGS_FFT=-lfftw3




# For performance analysis one can load
#    module load apps/scalasca
# on btrzx3, set
#    FC = scalasca -instrument $(FC)
# and compile as usual.
# This enables scalasca (http://scalasca.org) instrumentation. Default builds
# (without the scalasca wrapper) remain fully optimized
# and without instrumentation.


CC = mpiicc
FC = mpiifort
# To compile with the HDF5 lib:
LD = $(FC)

MPI = mpi2

IMPLICIT=umfpack

#SMP = OPENMP
#FFLAGS_OMP = -openmp


FFLAGS_DOUBLE = -r8
FFLAGS_OPT    = -O3 -no-prec-div
# further interesting optimisation options for the intel compiler:
# -ip -ipo
FFLAGS_DEBUG  = -g -O0 -traceback -ftrapuv -check pointers,stack,uninit,bounds -fpe0 -fp-stack-check
# further interesting debugging options for the intel compiler:
# -warn all -fpe0 -fpe3 -check bounds -WB -check uninit
# for the 2015 version of the intel compiler:
FFLAGS_OTHER  = -qopt-report=4 -qopt-report-phase=vec -warn all
# for the 2013 version of the intel compiler:
#FFLAGS_OTHER  = -vec-report0

FFLAGS_INC = -I. \
             $(FFLAGS_INC_SLEPC) $(FFLAGS_INC_H5) $(FFLAGS_INC_FFT) $(FFLAGS_INC_MKL) $(FFLAGS_INC_LIBRSB)
LDFLAGS    = $(LDFLAGS_SLEPC) $(LDFLAGS_FFT) $(LDFLAGS_H5) $(LDFLAGS_MKL) $(LDFLAGS_LIBRSB) \
             -L../../libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L../../libs/AMD/Lib/ -lamd

