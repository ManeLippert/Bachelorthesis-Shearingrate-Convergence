## ---------- RZG Hydra Cluster IPP Garching --------------
# Eigenvalue solver 

# 
MPI = mpi2
FFTLIB =FFT_FFTW3
REQUIREDMODULES = intel mpi.ibm fftw/3.2.2
CC = cc
FC = mpiifort
FFLAGS_OPT = -O3 -no-prec-div -xHost
FFLAGS_DEBUG = -g -C -traceback
FFLAGS_DOUBLE = -r8
FFLAGS_OMP = -openmp
FFLAGS_INC = -I/usr/include -I$(FFTW_HOME)/include -I./
FFLAGS_OTHER = -vec-report0
GKW_LIBS=/u/fjc/gkw/libs
LDFLAGS = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -L/$(GKW_LIBS)/UMFPACK/Lib/ -lumfpack_gkw64_NB -L/${GKW_LIBS}/AMD/Lib/ -lamd
IMPLICIT=umfpack
SMP=OPENMP
PERF=perf

# Comment the lines below if you don't need the eigenvalue solver
REQUIREDMODULES+= mkl/11.0 petsc-cplx/3.4.2 slepc-cplx/3.4.2
SLEPC = HAVE_SLEPC
LDFLAGS+= -L${PETSC_DIR}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -L/${MKL_HOME}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lslepc -lpetsc -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -mkl 
LDFLAGS+= -Wl,-rpath,${PETSC_DIR}/lib -Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib, -Wl,-rpath,${MKL_HOME}/lib/intel64
FFLAGS_INC+= -I/${PETSC_DIR}/include -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include
