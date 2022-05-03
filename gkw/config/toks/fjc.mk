## ---------- TOK cluster IPP Garching --------------
# 
MPI = usempi2
FFTLIB =FFT_FFTW3
REQUIREDMODULES = impi/4.1.3  intel/14.0 fftw/3.3.3
CC = cc
FC = mpiifort
#FFLAGS_OPT = -O3 -no-prec-div
DEBUG = off # fails with intel 11.1
FFLAGS_OPT = -O2
FFLAGS_DEBUG = -g -C -traceback -check noarg_temp_created #-fpe0
FFLAGS_DOUBLE = -r8
FFLAGS_OMP = -openmp
FFLAGS_INC = -I/usr/include -I$(FFTW_HOME)/include -I./
FFLAGS_OTHER = -vec-report0 #-std95 #-warn stderrors
LDFLAGS = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -L/afs/ipp/home/f/fjc/gkw/libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L/afs/ipp/home/f/fjc/gkw/libs/AMD/Lib/ -lamd
IMPLICIT=umfpack
PERF=perf
#SMP=OPENMP

# Comment the lines below if you don't need the eigenvalue solver
# Needs the exact correct modules to be loaded, and intel 14.0 compiler
REQUIREDMODULES+= mkl/10.3  slepc-cplx/3.6.0  petsc-cplx/3.6.0
SLEPC = HAVE_SLEPC
LDFLAGS+= -L${PETSC_DIR}/lib -L${SLEPC_DIR}/intel14.0-impi413-double-cplx/lib -L/${MKL_HOME}/lib/intel64 -lslepc -lpetsc -mkl #-llapack -lblas -lX11
LDFLAGS+= -Wl,-rpath,${PETSC_DIR}/lib -Wl,-rpath,${SLEPC_DIR}/intel14.0-impi413-double-cplx/lib, -Wl,-rpath,${MKL_HOME}/lib/intel64
FFLAGS_INC+= -I/${PETSC_DIR}/include -I${SLEPC_DIR}/intel14.0-impi413-double-cplx/include -I${SLEPC_DIR}/include 
