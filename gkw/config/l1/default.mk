## -- Bayreuth Rechenzentrum Cluster btrzx1 using Intel compiler

OPTFLAGS = on 
DEBUG = off 
PERF = perf


# The first three modules are mandatory, while the latter two are optional (needed if HAVE_HDF5 and/or HAVE_MKL is set)
REQUIREDiMODULES = intel/20.4.304 fftw/3.3.10 mpich/3.3.1 hdf5/1.12.1 inteloneapi/mkl/2022.0.1


# Uncomment if you need the eigenvalue solver
#SLEPC = HAVE_SLEPC
#PETSC_DIR=
#SLEPC_DIR=
#FFLAGS_INC_SLEPC= -I$(PETSC_DIR)/include -I$(SLEPC_DIR)/include
#LDFLAGS_SLEPC=-Wl,-rpath,${SLEPC_DIR}/lib -L${SLEPC_DIR}/lib -lslepc -Wl,-rpath,${PETSC_DIR}/lib -L${PETSC_DIR}/lib -lpetsc -lflapack -lfblas


# Uncomment if you need the HDF5 format - the paths do not yet work for everybody
IO_LIB = HAVE_HDF5
# The compiler flags for HDF5 are taken from the wrapper provided by HDF5
FFLAGS_INC_H5= $(shell h5fc -show | sed 's/^[^ ]*//')
LDFLAGS_H5=$(shell h5fc -show | sed 's/^[^ ]*//')


# Inlcude FFTW3 library
FFTLIB = FFT_FFTW3
# Directory where the fft library is located
FFTW3_DIR = /opt/ohpc/pub/libs/intel/fftw/fftw-3.3.10
# include and library flags
FFLAGS_INC_FFT=-I$(FFTW3_DIR)/include
LDFLAGS_FFT=-L$(FFTW3_DIR)/lib -lfftw3


# # Include the MKL library with (potentially faster)
# # sparse linear algebra routines
# LIBMKL = HAVE_MKL
# # the library path
# MKL_HOME=/opt/ohpc/pub/intel2022_1/mkl/2022.0.1
# # MKL_HOME environment variables points to the directory where the MKL
# # library is located. It should be set automatically when the module is loaded.
# FFLAGS_INC_MKL= -I$(MKL_HOME)/include/intel64/lp64  -I${MKL_HOME}/include -module ./
# LDFLAGS_MKL = -fpp -L "${MKL_HOME}/lib/intel64" -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -L "${MKL_HOME}/../compiler/lib/intel64" -liomp5 -lm -module ./
# # # The problem here is that there does not seem to be a mkl_spblas.mod
# # file ready to use. We must compile it together with GKW. See also
# # the corresponding target in gkw.mk .
# MKL_SPBLAS_MODULE_SRC = ${MKL_HOME}/include/mkl_spblas.f90


# Compiler
CC = mpiicc
FC = mpiifort
# To compile with the HDF5 lib:
LD = $(FC)

MPI = mpi2

IMPLICIT=umfpack


#SMP = OPENMP
#FFLAGS_OMP = -qopenmp


FFLAGS_DOUBLE = -r8
FFLAGS_OPT    = -O3 -no-prec-div
# further interesting optimisation options for the intel compiler:
# -ip -ipo
FFLAGS_DEBUG  = -g -O0 -traceback -ftrapuv -check pointers,stack,uninit,bounds -fpe0
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

