## -- Bayreuth Rechenzentrum Cluster btrzx2 using Intel compiler

OPTFLAGS = on 
DEBUG = off 
PERF=perf

# Uncomment if you need the eigenvalue solver
SLEPC = HAVE_SLEPC
PETSC_DIR=/data/peeters/shared/petsc-current
SLEPC_DIR=/data/peeters/shared/slepc-current
FFLAGS_INC_SLEPC= -I$(PETSC_DIR)/include -I$(SLEPC_DIR)/include
LDFLAGS_SLEPC=-Wl,-rpath,${SLEPC_DIR}/lib -L${SLEPC_DIR}/lib -lslepc -Wl,-rpath,${PETSC_DIR}/lib -L${PETSC_DIR}/lib -lpetsc -lflapack -lfblas

# Uncomment if you need the HDF5 format - the paths do not yet work for everybody
IO_LIB = HAVE_HDF5
# The compiler flags for HDF5 are taken from the wrapper provided by HDF5
FFLAGS_INC_H5= $(shell h5fc -show | sed 's/^[^ ]*//')
LDFLAGS_H5=$(shell h5fc -show | sed 's/^[^ ]*//')

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


# the following are not enforced, just suggested.
REQUIREDMODULES = intel-cluster-studio-2016 mpi/intel/5.1.3.210 lib/hdf5/5-1.8.18_serial_intel lib/scalapack

FFLAGS_DOUBLE = -r8
FFLAGS_OPT    = -O3 -no-prec-div
# further interesting optimisation options for the intel compiler:
# -ip -ipo
FFLAGS_DEBUG  = -g -O0 -traceback -ftrapuv -check pointers,stack,uninit,bounds -fpe0
# further interesting debugging options for the intel compiler:
# -warn all -fpe0 -fpe3 -check bounds -WB -check uninit
# for the 2015 version of the intel compiler:
FFLAGS_OTHER  = -qopt-report=4 -qopt-report-phase=vec
# for the 2013 version of the intel compiler:
#FFLAGS_OTHER  = -vec-report0

FFLAGS_INC = -I. \
             $(FFLAGS_INC_SLEPC) $(FFLAGS_INC_H5) $(FFLAGS_INC_FFT)
LDFLAGS    = $(LDFLAGS_SLEPC) $(LDFLAGS_FFT) $(LDFLAGS_H5) \
             -L../../libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L../../libs/AMD/Lib/ -lamd

