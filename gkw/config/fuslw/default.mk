## ---------- CCFE Fusion Unix Network  --------------
# 
MPI = mpi2
FFTLIB =FFT_FFTW3
REQUIREDMODULES = openmpi/1.4.2  ifort/10.0.023
CC = cc
FC = mpif90

#FFLAGS_OPT = -O3 -no-prec-div # These optimisations break the miller geometry test case with intel compiler 12.1
FFLAGS_OPT = -O2
FFLAGS_DEBUG = -g -C -traceback
FFLAGS_DOUBLE = -r8
FFLAGS_OMP = -openmp
FFLAGS_INC = -I/usr/include -I$(FFTW_HOME)/include -I./
FFLAGS_OTHER = -vec-report0

GKW_LIBS=/common/projects/codes/turbulence/GKW/gkw/libs
LDFLAGS = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -L/$(GKW_LIBS)/UMFPACK/Lib/ -lumfpack_gkw64_NB -L/$(GKW_LIBS)/AMD/Lib/ -lamd
IMPLICIT=umfpack
SMP=OPENMP
PERF=perf
