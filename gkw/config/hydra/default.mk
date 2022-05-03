## ---------- RZG Hydra Cluster IPP Garching --------------
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
GKW_LIBS=$(GKW_HOME)/libs
LDFLAGS = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f #-L/$(GKW_LIBS)/UMFPACK/Lib/ -lumfpack_gkw64_NB -L/${GKW_LIBS}/AMD/Lib/ -lamd
#IMPLICIT=umfpack
SMP=OPENMP
PERF=perf
